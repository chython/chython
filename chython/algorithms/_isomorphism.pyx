# -*- coding: utf-8 -*-
#
#  Copyright 2021 Ramil Nugmanov <nougmanoff@protonmail.com>
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

# long I:
# bond: single, double, triple, aromatic, special = 5 bit
# atom: H-Ce: 58 bit
# transfer bit

# long II:
# atom Pr-Og: 60 bit
# hybridizations: 1-4 = 4 bit

# long III:
# isotope: isotope - common_isotope -8 - +8 = 17 bit
# is_radical: 2 bit
# charge: -4 - +4: 9 bit
# implicit_hydrogens: 0-5 = 6 bit
# neighbors: 0-14 = 15 bit
# heteroatoms: 0-14 = 15 bit

# long IV:
# ring_sizes: not-in-ring bit, 3-atom ring, 4-...., 65-atom ring


@cython.boundscheck(False)
@cython.wraparound(False)
def get_mapping(unsigned long[::1] q_numbers not None, unsigned int[::1] q_back not None,
                unsigned long long[::1] q_masks1 not None, unsigned long long[::1] q_masks2 not None,
                unsigned long long[::1] q_masks3 not None, unsigned long long[::1] q_masks4 not None,
                unsigned int[::1] q_closures not None,
                unsigned int[::1] q_from not None, unsigned int[::1] q_to not None,

                unsigned long[::1] o_numbers not None,
                unsigned long long[::1] o_bits1 not None, unsigned long long[::1] o_bits2 not None,
                unsigned long long[::1] o_bits3 not None, unsigned long long[::1] o_bits4 not None,
                unsigned long long[::1] o_bonds not None,
                unsigned int[::1] o_from not None, unsigned int[::1] o_to not None,
                unsigned int[::1] o_indices not None,
                unsigned int[::1] scope not None):
    # expected less than 2^16 atoms in structure.
    cdef unsigned int stack = 0, path_size = 0, q_size, o_size, n, o_n, o_m, i, j, depth, front, back, has_closures
    cdef unsigned long long q_mask1, q_mask2, q_mask3, q_mask4, o_bond
    cdef dict mapping

    q_size = len(q_numbers)
    o_size = len(o_numbers)
    cdef int *stack_index = <int *> PyMem_Malloc(2 * q_size * sizeof(int))
    cdef int *stack_depth = <int *> PyMem_Malloc(2 * q_size * sizeof(int))
    cdef bint *matched = <bint *> PyMem_Malloc(o_size * sizeof(bint))
    q_size -= 1
    cdef int *path = <int *> PyMem_Malloc(q_size * sizeof(int))
    if not path or not stack_index or not stack_depth or not matched:
        raise MemoryError()
    memset(&matched[0], 0, o_size * sizeof(bint))

    # find entry-points.
    q_mask1 = q_masks1[0]
    q_mask2 = q_masks2[0]
    q_mask3 = q_masks3[0]
    q_mask4 = q_masks4[0]
    for n in range(o_size):
        if (scope[n] and
            q_mask1 & o_bits1[n] and
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

            if depth == q_size:
                mapping = _PyDict_NewPresized(q_size + 1)
                mapping[q_numbers[depth]] = o_numbers[n]
                for i in range(q_size):
                    mapping[q_numbers[i]] = o_numbers[path[i]]
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
                has_closures = q_closures[front]

                for i in range(o_from[n], o_to[n]):
                    o_bond = o_bonds[i]
                    o_n = o_indices[i]
                    if (scope[o_n] and not matched[o_n] and
                        q_mask1 & o_bond == o_bond and
                        q_mask2 & o_bits2[o_n] == o_bits2[o_n] and
                        q_mask3 & o_bits3[o_n] == o_bits3[o_n] and
                        q_mask4 & o_bits4[o_n]):

                        if has_closures:  # candidate atom should have same closures.
                            #for j in range(q_from[o_n], q_to[o_n]):
                            #    ...
                            # o_closures = o_bonds[o_n].keys() & reversed_mapping.keys()
                            # o_closures.discard(n)
                            #if o_closures == {mapping[m] for m, _ in query_closures[s_n]}:
                            #    obon = o_bonds[o_n]
                            #    if all(bond == obon[mapping[m]] for m, bond in query_closures[s_n]):
                            stack_index[stack] = o_n
                            stack_depth[stack] = front
                            stack += 1
                        else:  # candidate atom should not have closures.
                            for j in range(o_from[o_n], o_to[o_n]):
                               o_m = o_indices[j]
                               if o_m != n and matched[o_m]:
                                   break  # found closure
                            else:
                                stack_index[stack] = o_n
                                stack_depth[stack] = front
                                stack += 1
    finally:
        PyMem_Free(path)
        PyMem_Free(matched)
        PyMem_Free(stack_index)
        PyMem_Free(stack_depth)
