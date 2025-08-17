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
from cpython.mem cimport PyMem_Malloc, PyMem_Free, PyMem_Realloc
from libc.limits cimport UINT_MAX, USHRT_MAX
from libc.string cimport memset, memcpy


cdef struct paths_t:
    # store PID units. paths store only non-terminal nodes
    # in the case of distances equals 1 - don't store any paths
    unsigned short distance  # bonds count on the shortest path between nodes
    unsigned char num_pid0  # shortest paths count
    unsigned char num_pid1  # shortest+1 paths count
    unsigned short **pid0  # array of pointers to the shortest paths arrays
    unsigned short **pid1  # array of pointers to the shortest+1 paths arrays


cdef struct dist_matrix_t:
    unsigned short n_nodes
    unsigned short *mapping  # idx to node label mapping
    paths_t *data


cdef dist_matrix_t *alloc_dist_matrix(unsigned short n_nodes):
    cdef size_t i, total = n_nodes * n_nodes
    cdef dist_matrix_t *matrix = <dist_matrix_t *> PyMem_Malloc(sizeof(dist_matrix_t))

    matrix.n_nodes = n_nodes
    matrix.mapping = <unsigned short *> PyMem_Malloc(n_nodes * sizeof(unsigned short))
    matrix.data = <paths_t *> PyMem_Malloc(total * sizeof(paths_t))
    # set in paths_t all attrs to 0/NULL
    memset(matrix.data, 0, total * sizeof(paths_t))

    # reset distances to max
    for i in range(total):
        matrix.data[i].distance = USHRT_MAX
    return matrix


cdef void clear_pid0(paths_t *paths):
    cdef size_t i
    if paths.num_pid0 == 0: return
    for i in range(paths.num_pid0):
        PyMem_Free(paths.pid0[i])
    PyMem_Free(paths.pid0)
    paths.num_pid0 = 0
    paths.pid0 = NULL


cdef void clear_pid1(paths_t *paths):
    cdef size_t i
    if paths.num_pid1 == 0: return
    for i in range(paths.num_pid1):
        PyMem_Free(paths.pid1[i])
    PyMem_Free(paths.pid1)
    paths.num_pid1 = 0
    paths.pid1 = NULL


cdef void clear_pids(paths_t *paths):
    clear_pid0(paths)
    clear_pid1(paths)


cdef void free_dist_matrix(dist_matrix_t *matrix):
    cdef size_t i, total = matrix.n_nodes * matrix.n_nodes
    for i in range(total):
        clear_pids(&matrix.data[i])
    PyMem_Free(matrix.mapping)
    PyMem_Free(matrix.data)
    PyMem_Free(matrix)


cdef void move_paths(paths_t *paths):
    clear_pid1(paths)

    paths.num_pid1 = paths.num_pid0
    paths.pid1 = paths.pid0
    paths.num_pid0 = 0
    paths.pid0 = NULL


cdef unsigned short *concatenate_paths(paths_t *paths_ik, paths_t *paths_kj, unsigned short k):
    cdef unsigned short *path
    cdef size_t size1, size2

    path = <unsigned short *> PyMem_Malloc((paths_ik.distance + paths_kj.distance - 1) * sizeof(unsigned short))

    size1 = paths_ik.distance - 1
    size2 = paths_kj.distance - 1
    path[size1] = k  # put node k to form continuous i-k-j path excluding i,j terminals
    if size1:
        memcpy(path, paths_ik.pid0[0], size1 * sizeof(unsigned short))  # copy nodes between i-k
    if size2:
        memcpy(path + paths_ik.distance, paths_kj.pid0[0], size2 * sizeof(unsigned short))
    return path


cdef void append_pid0(paths_t *paths_ij, paths_t *paths_ik, paths_t *paths_kj, unsigned short k):
    cdef size_t i = paths_ij.num_pid0
    paths_ij.num_pid0 += 1
    paths_ij.pid0 = <unsigned short **> PyMem_Realloc(paths_ij.pid0, paths_ij.num_pid0 * sizeof(unsigned short *))
    paths_ij.pid0[i] = concatenate_paths(paths_ik, paths_kj, k)


cdef void append_pid1(paths_t *paths_ij, paths_t *paths_ik, paths_t *paths_kj, unsigned short k):
    cdef size_t i = paths_ij.num_pid1
    paths_ij.num_pid1 += 1
    paths_ij.pid1 = <unsigned short **> PyMem_Realloc(paths_ij.pid1, paths_ij.num_pid1 * sizeof(unsigned short *))
    paths_ij.pid1[i] = concatenate_paths(paths_ik, paths_kj, k)


cdef tuple ring_to_tuple(unsigned short *ring, size_t size):
    cdef size_t i
    return tuple(ring[i] for i in range(size))


cdef void free_hashes(unsigned long long **hashes, size_t size):
    cdef size_t i
    for i in range(size):
        PyMem_Free(hashes[i])
    PyMem_Free(hashes)


cdef int has_overlap(paths_t *paths_ik, paths_t *paths_kj, paths_t *paths_ij, int test_pid1):
    cdef size_t i, j
    cdef unsigned short n1, n2
    cdef unsigned short *path

    n1 = paths_ik.pid0[0][1]
    n2 = paths_kj.pid0[0][paths_kj.distance - 1]

    j = paths_ij.distance - 1
    for i in range(paths_ij.num_pid0):
        path = paths_ij.pid0[i]
        if path[1] == n1 or path[j] == n2:
            return 1
    if test_pid1:
        for i in range(paths_ij.num_pid1):
            path = paths_ij.pid1[i]
            if path[1] == n1 or path[paths_ij.distance] == n2:
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


cdef void build_pid(dist_matrix_t *matrix):
    cdef size_t i, j, k, sk, si
    cdef unsigned int d
    cdef paths_t *paths_kj
    cdef paths_t *paths_ij
    cdef paths_t *paths_ik

    for k in range(matrix.n_nodes):
        sk = k * matrix.n_nodes
        for i in range(matrix.n_nodes):
            if i == k: continue
            si = i * matrix.n_nodes
            paths_ik = &matrix.data[si + k]
            if paths_ik.distance == USHRT_MAX: continue

            for j in range(matrix.n_nodes):
                if j == k or j == i: continue
                paths_kj = &matrix.data[sk + j]
                if paths_kj.distance == USHRT_MAX: continue

                paths_ij = &matrix.data[si + j]
                d = paths_ik.distance + paths_kj.distance
                if d < paths_ij.distance:  # shorter pid0 path found
                    if d == paths_ij.distance - 1: # and not has_overlap(paths_ik, paths_kj, paths_ij, 0):
                        move_paths(paths_ij)
                    else:  # drop old paths
                        clear_pids(paths_ij)

                    paths_ij.distance = d
                    append_pid0(paths_ij, paths_ik, paths_kj, matrix.mapping[k])
                elif d == paths_ij.distance:  # new pid0 path
                    #if not has_overlap(paths_ik, paths_kj, paths_ij, 0):
                    append_pid0(paths_ij, paths_ik, paths_kj, matrix.mapping[k])
                elif d == paths_ij.distance + 1:  # new pid1 path
                    #if not has_overlap(paths_ik, paths_kj, paths_ij, 1):
                    append_pid1(paths_ij, paths_ik, paths_kj, matrix.mapping[k])


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


def sssr(dict graph, size_t n_rings):
    cdef size_t si
    cdef unsigned short i, n, m, n_nodes = len(graph)
    cdef unsigned short [USHRT_MAX] reverse_mapping  # 128kb
    cdef object mb
    cdef list output = []

    cdef dist_matrix_t *matrix = alloc_dist_matrix(n_nodes)

    for i, n in enumerate(graph):
        matrix.mapping[i] = n
        reverse_mapping[n] = i

    for n, mb in graph.items():
        si = n_nodes * reverse_mapping[n]
        for m in mb:
            matrix.data[si + reverse_mapping[m]].distance = 1

    build_pid(matrix)
    #
    # DEBUG
    cdef paths_t p
    for i in range(n_nodes):
        for j in range(n_nodes):
            if i == j: continue
            p = matrix.data[i * n_nodes + j]
            if p.distance == 1:
                print('!', matrix.mapping[i], matrix.mapping[j])
            for k in range(p.num_pid0):
                print('?', matrix.mapping[i], matrix.mapping[j], ring_to_tuple(p.pid0[k], p.distance - 1))
            for k in range(p.num_pid1):
                print('$', matrix.mapping[i], matrix.mapping[j], ring_to_tuple(p.pid1[k], p.distance))

    #
    # find_rings(rings, ring_sizes, pid0, pid1, dist, n_nodes, max_paths, hash_size, n_rings)
    #
    # for i in range(n_rings):
    #     output.append(convert_array_to_tuple(rings[i], node_mapping, ring_sizes[i]))
    #
    free_dist_matrix(matrix)
    return output
