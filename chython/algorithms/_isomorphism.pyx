# -*- coding: utf-8 -*-
#
#  Copyright 2021, 2022 Ramil Nugmanov <nougmanoff@protonmail.com>
#  Copyright 2021 Aleksandr Sizov <murkyrussian@gmail.com>
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
cimport cython
from cpython.mem cimport PyMem_Malloc, PyMem_Free
from libc.string cimport memset

cdef extern from "Python.h":
    dict _PyDict_NewPresized(Py_ssize_t minused)


@cython.boundscheck(False)
@cython.wraparound(False)
def get_mapping(unsigned long[::1] q_numbers not None, unsigned int[::1] q_back not None,
                unsigned long long[::1] q_masks1 not None, unsigned long long[::1] q_masks2 not None,
                unsigned long long[::1] q_masks3 not None, unsigned long long[::1] q_masks4 not None,
                unsigned int[::1] q_closures not None, unsigned int[::1] q_from not None,
                unsigned int[::1] q_to not None, unsigned int[::1] q_indices not None,
                unsigned long long[::1] q_bonds not None, unsigned long[::1] o_numbers not None,
                unsigned long long[::1] o_bits1 not None, unsigned long long[::1] o_bits2 not None,
                unsigned long long[::1] o_bits3 not None, unsigned long long[::1] o_bits4 not None,
                unsigned long long[::1] o_bonds not None, unsigned int[::1] o_from not None,
                unsigned int[::1] o_to not None, unsigned int[::1] o_indices not None,
                unsigned int[::1] scope not None):
    # expected less than 2^16 atoms in structure.
    cdef unsigned int stack = 0, path_size = 0, q_size, q_size_dec, o_size, depth, front, back, closures_num
    cdef unsigned int n, m, o, i, j, closures_counter
    cdef unsigned long long q_mask1, q_mask2, q_mask3, q_mask4, o_bond, c_bond
    cdef dict mapping

    q_size = len(q_numbers)
    q_size_dec = q_size - 1
    o_size = len(o_numbers)
    cdef unsigned int *path = <unsigned int *> PyMem_Malloc(q_size_dec * sizeof(unsigned int))
    cdef unsigned int *stack_index = <unsigned int *> PyMem_Malloc(2 * o_size * sizeof(unsigned int))
    cdef unsigned int *stack_depth = <unsigned int *> PyMem_Malloc(2 * o_size * sizeof(unsigned int))
    cdef bint *matched = <bint *> PyMem_Malloc(o_size * sizeof(bint))
    cdef unsigned long long *o_closures = <unsigned long long *> PyMem_Malloc(o_size * sizeof(unsigned long long))

    if not path or not stack_index or not stack_depth or not matched or not o_closures:
        raise MemoryError()

    memset(matched, 0, o_size * sizeof(bint))
    memset(o_closures, 0, o_size * sizeof(unsigned long long))

    # find entry-points.
    q_mask1 = q_masks1[0]
    q_mask2 = q_masks2[0]
    q_mask3 = q_masks3[0]
    q_mask4 = q_masks4[0]
    for n in range(o_size):
        if (scope[n] and
            q_mask1 & o_bits1[n] and  # o_bits1 doesn't contain bond bits.
            q_mask2 & o_bits2[n] == o_bits2[n] and
            q_mask3 & o_bits3[n] == o_bits3[n] and
            q_mask4 & o_bits4[n]):

            stack_index[stack] = n
            stack_depth[stack] = 0
            stack += 1

    try:
        while stack:
            stack -= 1
            depth = stack_depth[stack]
            n = stack_index[stack]

            if depth == q_size_dec:
                mapping = _PyDict_NewPresized(q_size)
                for i in range(depth):
                    mapping[q_numbers[i]] = o_numbers[path[i]]
                mapping[q_numbers[depth]] = o_numbers[n]
                yield mapping
            else:
                if path_size != depth:  # dead end reached
                    for i in range(depth, path_size):
                        matched[path[i]] = False  # mark unmatched
                    path_size = depth

                matched[n] = True
                path[path_size] = n
                path_size += 1

                front = depth + 1
                back = q_back[front]
                if back != depth:  # branch
                    n = path[back]

                # load next query atom
                q_mask1 = q_masks1[front]
                q_mask2 = q_masks2[front]
                q_mask3 = q_masks3[front]
                q_mask4 = q_masks4[front]
                closures_num = q_closures[front]

                for i in range(o_from[n], o_to[n]):
                    o_bond = o_bonds[i]
                    m = o_indices[i]
                    if (scope[m] and not matched[m] and
                        q_mask1 & o_bond == o_bond and  # bond order, in ring mark and atom bit should match.
                        q_mask2 & o_bits2[m] == o_bits2[m] and
                        q_mask3 & o_bits3[m] == o_bits3[m] and
                        q_mask4 & o_bits4[m]):

                        if closures_num:  # candidate atom should have same closures.
                            closures_counter = 0
                            # make a map of closures for o_n atom
                            # an index is a neighbor atom and a value is a bond between o_n and the neighbor
                            for j in range(o_from[m], o_to[m]):
                                o = o_indices[j]
                                if o != n and matched[o]:
                                    o_closures[o] = o_bonds[j]
                                    closures_counter += 1

                            if closures_counter == closures_num:
                                for j in range(q_from[front], q_to[front]):
                                    c_bond = o_closures[path[q_indices[j]]]
                                    if not c_bond or q_bonds[j] & c_bond != c_bond:  # compare order and ring bits
                                        break
                                else:
                                    stack_index[stack] = m
                                    stack_depth[stack] = front
                                    stack += 1

                            # fill an array with nulls
                            for j in range(o_from[m], o_to[m]):
                                o_closures[o_indices[j]] = 0
                        else:  # candidate atom should not have closures.
                            for j in range(o_from[m], o_to[m]):
                               o = o_indices[j]
                               if o != n and matched[o]:
                                   break  # found closure
                            else:
                                stack_index[stack] = m
                                stack_depth[stack] = front
                                stack += 1
    finally:
        PyMem_Free(path)
        PyMem_Free(matched)
        PyMem_Free(stack_index)
        PyMem_Free(stack_depth)
        PyMem_Free(o_closures)
