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
from libc.limits cimport USHRT_MAX
from libc.string cimport memset, memcpy


cdef enum:
    ULL_BITS = sizeof(unsigned long long) * 8


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


cdef struct ring_t:
    unsigned short n_nodes
    unsigned short *nodes
    unsigned long long *hash


cdef struct rings_t:
    unsigned short n_rings
    unsigned short n_allocated
    unsigned short hash_size
    ring_t *rings


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


cdef rings_t *alloc_rings(unsigned short n_rings, size_t hash_size):
    cdef size_t i
    cdef rings_t *rings = <rings_t *> PyMem_Malloc(sizeof(rings_t))

    rings.n_rings = n_rings
    rings.n_allocated = 0
    rings.hash_size = hash_size
    rings.rings = <ring_t *> PyMem_Malloc(n_rings * sizeof(ring_t))
    memset(rings.rings, 0, n_rings * sizeof(ring_t))

    for i in range(n_rings):
        rings.rings[i].n_nodes = USHRT_MAX
    return rings


cdef void free_pid0(paths_t *paths):
    cdef size_t i
    if paths.num_pid0 == 0: return
    for i in range(paths.num_pid0):
        PyMem_Free(paths.pid0[i])
    PyMem_Free(paths.pid0)
    paths.num_pid0 = 0
    paths.pid0 = NULL


cdef void free_pid1(paths_t *paths):
    cdef size_t i
    if paths.num_pid1 == 0: return
    for i in range(paths.num_pid1):
        PyMem_Free(paths.pid1[i])
    PyMem_Free(paths.pid1)
    paths.num_pid1 = 0
    paths.pid1 = NULL


cdef void free_pids(paths_t *paths):
    free_pid0(paths)
    free_pid1(paths)


cdef void free_dist_matrix(dist_matrix_t *matrix):
    cdef size_t i, total = matrix.n_nodes * matrix.n_nodes
    for i in range(total):
        free_pids(&matrix.data[i])
    PyMem_Free(matrix.mapping)
    PyMem_Free(matrix.data)
    PyMem_Free(matrix)


cdef void free_ring(ring_t *ring):
    PyMem_Free(ring.hash)
    PyMem_Free(ring.nodes)
    ring.hash = NULL
    ring.nodes = NULL


cdef void free_rings(rings_t *rings):
    cdef size_t i
    for i in range(rings.n_rings):
        free_ring(&rings.rings[i])
    PyMem_Free(rings.rings)
    PyMem_Free(rings)


cdef void move_paths(paths_t *paths):
    free_pid1(paths)

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


cdef int has_not_overlap(paths_t *paths_ik, paths_t *paths_kj, paths_t *paths_ij, unsigned short k, int test_pid1):
    cdef size_t i, j
    cdef unsigned short n1, n2
    cdef unsigned short *path

    if paths_ik.distance == 1:
        n1 = k
    else:
        n1 = paths_ik.pid0[0][0]
    if paths_kj.distance == 1:
        n2 = k
    else:
        n2 = paths_kj.pid0[0][paths_kj.distance - 2]

    j = paths_ij.distance - 2
    for i in range(paths_ij.num_pid0):
        path = paths_ij.pid0[i]
        if path[0] == n1 or path[j] == n2:
            return 0
    if test_pid1:
        j = paths_ij.distance - 1
        for i in range(paths_ij.num_pid1):
            path = paths_ij.pid1[i]
            if path[0] == n1 or path[j] == n2:
                return 0
    return 1


cdef unsigned long long *build_hash(unsigned short *path1, unsigned short *path2,
                                    unsigned short i, unsigned short j, size_t size1, size_t size2, size_t hash_size):
    cdef size_t k
    cdef unsigned int node
    cdef unsigned long long *hash = <unsigned long long *> PyMem_Malloc(hash_size * sizeof(unsigned long long))

    memset(hash, 0, hash_size * sizeof(unsigned long long))
    hash[i // ULL_BITS] |= (1 << (i % ULL_BITS))  # add linker node
    hash[j // ULL_BITS] |= (1 << (j % ULL_BITS))

    for k in range(size1):
        node = path1[k]
        hash[node // ULL_BITS] |= (1 << (node % ULL_BITS))
    for k in range(size2):
        node = path2[k]
        hash[node // ULL_BITS] |= (1 << (node % ULL_BITS))
    return hash


cdef unsigned short *build_ring(unsigned short *path1, unsigned short *path2,
                                unsigned short i, unsigned short j, size_t size1, size_t size2):
    cdef size_t k, s = size1 + 1
    cdef unsigned short *nodes

    nodes = <unsigned short *> PyMem_Malloc((size1 + size2 + 2) * sizeof(unsigned short))
    nodes[0] = i
    nodes[s] = j
    memcpy(nodes + 1, path1, size1 * sizeof(unsigned short))

    # inverse second path
    for k in range(size2 - 1, -1, -1):
        s += 1
        nodes[s] = path2[k]
    return nodes


cdef void push_ring(ring_t ring, rings_t *rings):
    cdef size_t i

    if rings.n_allocated == rings.n_rings:
        # drop last ring
        free_ring(&rings.rings[rings.n_rings - 1])
    else:
        rings.n_allocated += 1

    for i in range(rings.n_allocated - 1, rank, -1):
        rings.rings[i] = rings.rings[i - 1]

    # insert the new ring at the rank position
    rings.rings[rank] = ring


cdef int is_unique(ring_t *ring, rings_t *rings):
    cdef size_t i, j
    cdef int hash_match
    cdef ring_t *other

    for i in range(rings.n_rings):
        other = &rings.rings[i]
        if ring.n_nodes == other.n_nodes:
            hash_match = 1
            for j in range(rings.hash_size):
                if ring.hash[j] ^ other.hash[j] != 0:
                    hash_match = 0
                    break
            if hash_match:
                return 0
        elif ring.n_nodes < other.n_nodes:
            return 1
    return 0


cdef void build_pid(dist_matrix_t *matrix):
    cdef size_t i, j, k, sk, si
    cdef unsigned short rk
    cdef unsigned int d
    cdef paths_t *paths_kj
    cdef paths_t *paths_ij
    cdef paths_t *paths_ik

    for k in range(matrix.n_nodes):
        sk = k * matrix.n_nodes
        rk = matrix.mapping[k]
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
                    if d == paths_ij.distance - 1 and has_not_overlap(paths_ik, paths_kj, paths_ij, rk, 0):
                        move_paths(paths_ij)
                    else:  # drop old paths
                        free_pids(paths_ij)

                    paths_ij.distance = d
                    append_pid0(paths_ij, paths_ik, paths_kj, rk)
                elif d == paths_ij.distance:  # new pid0 path
                    if has_not_overlap(paths_ik, paths_kj, paths_ij, rk, 0):
                        append_pid0(paths_ij, paths_ik, paths_kj, rk)
                elif d == paths_ij.distance + 1:  # new pid1 path
                    if has_not_overlap(paths_ik, paths_kj, paths_ij, rk, 1):
                        append_pid1(paths_ij, paths_ik, paths_kj, rk)


cdef void find_rings(dist_matrix_t *matrix, rings_t *rings):
    cdef size_t i, j, k, si
    cdef unsigned short ri, rj, d
    cdef ring_t ring
    cdef paths_t *path

    for i in range(matrix.n_nodes):
        si = i * matrix.n_nodes
        ri = matrix.mapping[i]
        for j in range(matrix.n_nodes):
            if i == j: continue
            path = &matrix.data[si + j]
            d = path.distance
            if d == USHRT_MAX: continue  # different components
            rj = matrix.mapping[j]

            if d == 1:  # triangles
                if not path.num_pid1: continue
                ring.n_nodes = 3
                for k in range(path.num_pid1):
                    ring.hash = build_hash(NULL, path.pid1[k], ri, rj, 0, 1, rings.hash_size)

                    if is_unique(&ring, rings):
                        ring.nodes = build_ring(path.pid1[k], NULL, ri, rj, 1, 0)
                        push_ring(ring, rings)
                    else:
                        PyMem_Free(ring.hash)
            elif path.num_pid0 == 1:  # is odd?
                if not path.num_pid1: continue
                ring.n_nodes = 2 * d + 1
                for k in range(path.num_pid1):
                    ring.hash = build_hash(path.pid0[0], path.pid1[k], ri, rj, d - 1, d, rings.hash_size)

                    if is_unique(&ring, rings):
                        ring.nodes = build_ring(path.pid0[0], path.pid1[k], ri, rj, d - 1, d)
                        push_ring(ring, rings)
                    else:
                        PyMem_Free(ring.hash)
            else:  # is even
                ring.n_nodes = 2 * d
                d -= 1
                for k in range(1, path.num_pid0):
                    ring.hash = build_hash(path.pid0[k - 1], path.pid0[k], ri, rj, d, d, rings.hash_size)

                    if is_unique(&ring, rings):
                        ring.nodes = build_ring(path.pid0[k - 1], path.pid0[k], ri, rj, d, d)
                        push_ring(ring, rings)
                    else:
                        PyMem_Free(ring.hash)


def sssr(dict graph, size_t n_rings):
    cdef size_t si
    cdef unsigned short i, n, m, n_max = 0, n_nodes = len(graph)
    cdef unsigned short [USHRT_MAX] reverse_mapping  # 128kb
    cdef ring_t ring
    cdef rings_t *rings
    cdef object mb
    cdef list output = []

    cdef dist_matrix_t *matrix = alloc_dist_matrix(n_nodes)

    for i, n in enumerate(graph):
        matrix.mapping[i] = n
        reverse_mapping[n] = i
        if n > n_max: n_max = n

    for n, mb in graph.items():
        si = n_nodes * reverse_mapping[n]
        for m in mb:
            matrix.data[si + reverse_mapping[m]].distance = 1

    rings = alloc_rings(n_rings, (n_max + ULL_BITS - 1) // ULL_BITS)

    # run PID matrix calculation
    build_pid(matrix)

    # DEBUG
    # cdef paths_t p
    # for i in range(n_nodes):
    #     for j in range(n_nodes):
    #         if i == j: continue
    #         p = matrix.data[i * n_nodes + j]
    #         if p.distance == 1:
    #             print('!', matrix.mapping[i], matrix.mapping[j])
    #         for k in range(p.num_pid0):
    #             print('?', matrix.mapping[i], matrix.mapping[j], ring_to_tuple(p.pid0[k], p.distance - 1))
    #         for k in range(p.num_pid1):
    #             print('$', matrix.mapping[i], matrix.mapping[j], ring_to_tuple(p.pid1[k], p.distance))

    find_rings(matrix, rings)
    for i in range(n_rings):
        ring = rings.rings[i]
        output.append(ring_to_tuple(ring.nodes, ring.n_nodes))

    free_dist_matrix(matrix)
    free_rings(rings)
    return output
