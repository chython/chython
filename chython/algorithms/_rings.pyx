# -*- coding: utf-8 -*-
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
from cpython.mem cimport PyMem_Malloc, PyMem_Free
from libc.limits cimport UINT_MAX
from libc.string cimport memset, memcpy


cdef extern from "Python.h":
    dict _PyDict_NewPresized(Py_ssize_t minused)


cdef unsigned int ** alloc_pid(size_t n_nodes, size_t max_paths):
    cdef unsigned int **pid = <unsigned int **> PyMem_Malloc(n_nodes * n_nodes * max_paths * sizeof(unsigned int *))
    memset(pid, 0, n_nodes * n_nodes * max_paths * sizeof(unsigned int *))
    return pid


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


def make_pids(dict graph, size_t max_paths=10):
    cdef size_t n_nodes = len(graph), i, j, k, sk, si, ij, kj, ik
    cdef unsigned int d, dki, dkj, dij
    cdef object n, m, mb
    cdef dict reverse_mapping
    cdef unsigned int *path
    cdef unsigned int *path1
    cdef unsigned int *path2

    cdef unsigned int *node_mapping = <unsigned *> PyMem_Malloc(n_nodes * sizeof(unsigned int))
    cdef unsigned int *dist = <unsigned int *> PyMem_Malloc(n_nodes * n_nodes * sizeof(unsigned int))
    cdef unsigned int **pid0 = alloc_pid(n_nodes, max_paths)
    cdef unsigned int **pid1 = alloc_pid(n_nodes, max_paths)

    memset(dist, 255, n_nodes * n_nodes * sizeof(unsigned int))

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
            path[0] = n
            path[1] = m
            append_path(pid0, ij, max_paths, path)

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

    # DEBUG
    # for i in range(n_nodes):
    #     for j in range(n_nodes):
    #         if i == j: continue
    #         for k in range((i * n_nodes + j) * max_paths, (i * n_nodes + j + 1) * max_paths):
    #             path = pid0[k]
    #             if path != NULL:
    #                 print('!', node_mapping[i], node_mapping[j], [path[d] for d in range(dist[i * n_nodes + j] + 1)])
    #             path = pid1[k]
    #             if path != NULL:
    #                 print('?', node_mapping[i], node_mapping[j], [path[d] for d in range(dist[i * n_nodes + j] + 2)])

    PyMem_Free(node_mapping)
    PyMem_Free(dist)
    free_pid(pid0, n_nodes, max_paths)
    free_pid(pid1, n_nodes, max_paths)
