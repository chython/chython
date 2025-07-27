# cython: undeclared_check_usage=error
# cython: warn.undeclared=True
# cython: warn.unused=True
# cython: warn.unused_arg=True
# cython: warn.maybe_uninitialized=True
# cython: boundscheck=False
# cython: wraparound=False

#
#  Copyright 2025 Ramil Nugmanov <nougmanoff@protonmail.com>
#  This file is part of chython.
#
#  chython is free software; you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with this program; if not, see <https://www.gnu.org/licenses/>.
#
from cpython.tuple cimport PyTuple_New, PyTuple_SET_ITEM
from cpython.mem cimport PyMem_Malloc, PyMem_Free
from libc.limits cimport UINT_MAX
from libc.string cimport memset, memcpy


cdef extern from "Python.h":
    dict _PyDict_NewPresized(Py_ssize_t minused)


cdef unsigned int ** alloc_pid(size_t n_nodes, size_t max_paths):
    cdef unsigned int **pid = <unsigned int **> PyMem_Malloc(n_nodes * n_nodes * max_paths * sizeof(unsigned int *))
    memset(pid, 0, n_nodes * n_nodes * max_paths * sizeof(unsigned int *))
    return pid


cdef unsigned int * alloc_dist(size_t n_nodes):
    cdef unsigned int *dist = <unsigned int *> PyMem_Malloc(n_nodes * n_nodes * sizeof(unsigned int))
    memset(dist, 255, n_nodes * n_nodes * sizeof(unsigned int))
    return dist


cdef void free_pid(unsigned int **pid, size_t n_nodes, size_t max_paths):
    cdef size_t i
    for i in range(n_nodes * n_nodes * max_paths):
        PyMem_Free(pid[i])
    PyMem_Free(pid)


cdef void free_paths(unsigned int **pid, size_t ij, size_t max_paths):
    cdef size_t i

    ij *= max_paths
    for i in range(ij, ij + max_paths):
        if pid[i] == NULL:
            break
        PyMem_Free(pid[i])
        pid[i] = NULL


cdef void free_rings(unsigned int **rings, size_t size):
    cdef size_t i
    for i in range(size):
        PyMem_Free(rings[i])
    PyMem_Free(rings)


cdef void free_hashes(unsigned long long **hashes, size_t size):
    cdef size_t i
    for i in range(size):
        PyMem_Free(hashes[i])
    PyMem_Free(hashes)


cdef void move_paths(unsigned int **pid0, unsigned int **pid1, size_t ij, size_t max_paths):
    cdef size_t i

    ij *= max_paths
    for i in range(ij, ij + max_paths):
        # run max_path times to make sure all pid1 paths overridden
        PyMem_Free(pid1[i])
        pid1[i] = pid0[i]
        pid0[i] = NULL


cdef void append_path(unsigned int **pid, size_t ij, size_t max_paths, unsigned int *path):
    cdef size_t i

    ij *= max_paths
    for i in range(ij, ij + max_paths):
        if pid[i] == NULL:
            pid[i] = path
            return
    raise MemoryError('Reached max allowed paths')


cdef unsigned int * concatenate_paths(unsigned int **pid, unsigned int *dist, size_t i, size_t j, size_t max_paths):
    cdef unsigned int *path
    cdef unsigned int size1, size2

    size1 = dist[i]
    size2 = dist[j] + 1  # sizes are edges counts in a path, not nodes, thus, +1

    path = <unsigned int *> PyMem_Malloc((size1 + size2) * sizeof(unsigned int))

    memcpy(path, pid[i * max_paths], size1 * sizeof(unsigned int))  # dropping last node
    memcpy(path + size1, pid[j * max_paths], size2 * sizeof(unsigned int))
    return path


cdef int has_overlap(unsigned int **pid0, unsigned int **pid1, unsigned int *dist,
                     size_t i, size_t j, size_t k, size_t max_paths, size_t shift):
    cdef size_t n
    cdef unsigned int n1, n2
    cdef unsigned int d
    cdef unsigned int *path

    n1 = pid0[i * max_paths][1]
    n2 = pid0[j * max_paths][dist[j] - 1]

    d = dist[k] - shift
    k *= max_paths
    for n in range(k, k + max_paths):
        path = pid1[n]
        if path == NULL:
            return 0
        if path[1] == n1 or path[d] == n2:
            return 1
    return 0


cdef unsigned int * build_ring(unsigned int *path1, unsigned int *path2, size_t size1, size_t size2):
    cdef size_t i
    cdef unsigned int *path

    # size 1 and 2 are edges counts. paths have 2 overlapped atoms
    path = <unsigned int *> PyMem_Malloc((size1 + size2) * sizeof(unsigned int))

    memcpy(path, path1, size1 * sizeof(unsigned int))  # drop last node

    # inverse second path
    for i in range(size2, 0, -1):  # drop first node
        path[size1] = path2[i]
        size1 += 1
    return path


cdef unsigned long long * build_hash(unsigned int *path1, unsigned int *path2, size_t size1, size_t size2, size_t hash_size):
    cdef size_t i
    cdef unsigned int node
    cdef unsigned long long *hash = <unsigned long long *> PyMem_Malloc(hash_size * sizeof(unsigned long long))

    memset(hash, 0, hash_size * sizeof(unsigned long long))

    for i in range(size1):  # size is edge count; thus, the last node is dropped
        node = path1[i]
        hash[node // 64] |= (1 << (node % 64))
    for i in range(size2, 0, -1):  # 1st node is dropped
        node = path2[i]
        hash[node // 64] |= (1 << (node % 64))
    return hash


cdef size_t get_rank(unsigned int *rings, size_t size, size_t n_rings):
    cdef size_t i
    for i in range(n_rings):
        if size < rings[i]:
            return i
    return UINT_MAX


cdef int compare_rings(unsigned long long *hash1, unsigned long long *ring2, size_t hash_size):
    cdef size_t i
    for i in range(hash_size):
        if hash1[i] ^ ring2[i] != 0:
            return 1  # doesn't match
    return 0


cdef int check_ring_existence(unsigned long long *ring, unsigned long long **rings, unsigned int *ring_sizes,
                              size_t ring_size, size_t hash_size, size_t n_rings):
    cdef size_t i
    cdef unsigned int size
    for i in range(n_rings):
        size = ring_sizes[i]
        if size == ring_size:
            if compare_rings(ring, rings[i], hash_size) == 0:
                return 0
        elif size > ring_size:
            return 1
    return 1


cdef void push_ring(unsigned int **rings, unsigned long long **ring_hashes, unsigned int *ring_sizes,
                    unsigned int *ring, unsigned long long *hash, unsigned int size, size_t rank, size_t n_rings):
    cdef size_t i, i1, n1

    n1 = n_rings - 1
    PyMem_Free(rings[n1])
    PyMem_Free(ring_hashes[n1])

    for i in range(n1, rank, -1):
        i1 = i - 1
        rings[i] = rings[i1]
        ring_hashes[i] = ring_hashes[i1]
        ring_sizes[i] = ring_sizes[i1]

    # Insert the new ring at the rank position
    rings[rank] = ring
    ring_hashes[rank] = hash
    ring_sizes[rank] = size


cdef void add_ring_if_unique(unsigned int *path1, unsigned int *path2,
                             unsigned int **rings, unsigned int *ring_sizes, unsigned long long **ring_hashes,
                             size_t rank, size_t size, size_t size1, size_t size2, size_t hash_size, size_t n_rings):
    cdef unsigned int *ring
    cdef unsigned long long * hash = build_hash(path1, path2, size1, size2, hash_size)

    if check_ring_existence(hash, ring_hashes, ring_sizes, size, hash_size, n_rings) == 1:
        ring = build_ring(path1, path2, size1, size2)
        push_ring(rings, ring_hashes, ring_sizes, ring, hash, size, rank, n_rings)
    else:
        PyMem_Free(hash)


cdef tuple convert_array_to_tuple(unsigned int *array, unsigned int *node_mapping, size_t size):
    cdef size_t i
    cdef tuple output = PyTuple_New(size)

    for i in range(size):
        PyTuple_SET_ITEM(output, i, node_mapping[array[i]])
    return output


cdef void build_pid(unsigned int **pid0, unsigned int **pid1, unsigned int *dist, size_t n_nodes, size_t max_paths):
    cdef size_t i, j, k, sk, si, ij, kj, ik
    cdef unsigned int d, dki, dkj, dij
    cdef unsigned int *path

    for k in range(n_nodes):
        sk = k * n_nodes
        for i in range(n_nodes):
            if i == k: continue
            dki = dist[sk + i]
            if dki == UINT_MAX: continue

            si = i * n_nodes
            ik = si + k
            for j in range(n_nodes):
                if j == k or j == i: continue
                kj = sk + j
                dkj = dist[kj]
                if dkj == UINT_MAX: continue

                ij = si + j
                dij = dist[ij]
                d = dki + dkj
                if d < dij:  # shorter pid0 path found
                    dist[ij] = d

                    if d == dij - 1 and not has_overlap(pid0, pid0, dist, ik, kj, ij, max_paths, 1):
                        move_paths(pid0, pid1, ij, max_paths)
                    else:  # override old paths
                        free_paths(pid0, ij, max_paths)
                        free_paths(pid1, ij, max_paths)

                    path = concatenate_paths(pid0, dist, ik, kj, max_paths)
                    append_path(pid0, ij, max_paths, path)
                elif d == dij:  # new pid0 path
                    if not has_overlap(pid0, pid0, dist, ik, kj, ij, max_paths, 1):
                        path = concatenate_paths(pid0, dist, ik, kj, max_paths)
                        append_path(pid0, ij, max_paths, path)
                elif d == dij + 1:  # new pid1 path
                    if not has_overlap(pid0, pid0, dist, ik, kj, ij, max_paths, 1):
                        if not has_overlap(pid0, pid1, dist, ik, kj, ij, max_paths, 0):
                            path = concatenate_paths(pid0, dist, ik, kj, max_paths)
                            append_path(pid1, ij, max_paths, path)


cdef void find_rings(unsigned int **rings, unsigned int *ring_sizes,
                     unsigned int **pid0, unsigned int **pid1, unsigned int *dist,
                     size_t n_nodes, size_t max_paths, size_t hash_size, size_t n_rings):
    cdef size_t i, j, k, si, ij, size, rank
    cdef unsigned int d, d1
    cdef unsigned int *path1
    cdef unsigned int *path2
    cdef unsigned long long **ring_hashes = <unsigned long long **> PyMem_Malloc(n_rings * sizeof(unsigned long long *))

    memset(ring_hashes, 0, n_rings * sizeof(unsigned long long *))
    memset(ring_sizes, 255, n_rings * sizeof(unsigned int))
    memset(rings, 0, n_rings * sizeof(unsigned int *))

    for i in range(n_nodes):
        si = i * n_nodes
        for j in range(n_nodes):
            if i == j: continue
            ij = si + j
            d = dist[ij]
            if d == UINT_MAX: continue  # different components

            ij *= max_paths
            path2 = pid0[ij + 1]
            if path2 == NULL:  # is odd?
                path2 = pid1[ij]
                if path2 == NULL: continue  # not ring
                size = 2 * d + 1
                rank = get_rank(ring_sizes, size, n_rings)
                if rank == UINT_MAX: continue

                d1 = d + 1
                path1 = pid0[ij]
                add_ring_if_unique(path1, path2, rings, ring_sizes, ring_hashes, rank, size, d, d1, hash_size, n_rings)

                for k in range(ij + 1, ij + max_paths):
                    path2 = pid1[k]
                    if path2 == NULL: break
                    add_ring_if_unique(path1, path2, rings, ring_sizes, ring_hashes, rank, size, d, d1, hash_size, n_rings)
            else:  # is even
                size = 2 * d
                rank = get_rank(ring_sizes, size, n_rings)
                if rank == UINT_MAX: continue  # is not smaller than we have already

                path1 = pid0[ij]
                add_ring_if_unique(path1, path2, rings, ring_sizes, ring_hashes, rank, size, d, d, hash_size, n_rings)

                for k in range(ij + 2, ij + max_paths):
                    path2 = pid0[k]
                    if path2 == NULL: break
                    add_ring_if_unique(path1, path2, rings, ring_sizes, ring_hashes, rank, size, d, d, hash_size, n_rings)
    free_hashes(ring_hashes, n_rings)


def sssr(dict graph, size_t n_rings, size_t max_paths=10):
    cdef size_t n_nodes = len(graph), i, j, si, ij
    cdef size_t hash_size = (n_nodes + 63) // 64
    cdef object n, m, mb
    cdef dict reverse_mapping
    cdef list output = []
    cdef unsigned int *path

    cdef unsigned int *node_mapping = <unsigned int *> PyMem_Malloc(n_nodes * sizeof(unsigned int))
    cdef unsigned int *dist = alloc_dist(n_nodes)
    cdef unsigned int **pid0 = alloc_pid(n_nodes, max_paths)
    cdef unsigned int **pid1 = alloc_pid(n_nodes, max_paths)

    cdef unsigned int *ring_sizes = <unsigned int *> PyMem_Malloc(n_rings * sizeof(unsigned int))
    cdef unsigned int **rings = <unsigned int **> PyMem_Malloc(n_rings * sizeof(unsigned int *))

    reverse_mapping = _PyDict_NewPresized(n_nodes)
    for i, n in enumerate(graph):
        node_mapping[i] = n
        reverse_mapping[n] = i

    for n, mb in graph.items():
        i = reverse_mapping[n]
        si = i * n_nodes
        for m in mb:
            j = reverse_mapping[m]
            ij = si + j
            dist[ij] = 1
            path = <unsigned int *> PyMem_Malloc(2 * sizeof(unsigned int))
            path[0] = i
            path[1] = j
            append_path(pid0, ij, max_paths, path)

    build_pid(pid0, pid1, dist, n_nodes, max_paths)
    find_rings(rings, ring_sizes, pid0, pid1, dist, n_nodes, max_paths, hash_size, n_rings)

    for i in range(n_rings):
        output.append(convert_array_to_tuple(rings[i], node_mapping, ring_sizes[i]))

    PyMem_Free(node_mapping)
    PyMem_Free(dist)
    PyMem_Free(ring_sizes)
    free_rings(rings, n_rings)
    free_pid(pid0, n_nodes, max_paths)
    free_pid(pid1, n_nodes, max_paths)

    return output
