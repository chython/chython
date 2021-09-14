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
cimport cpython.array
from cpython.mem cimport PyMem_Malloc, PyMem_Free
from libc.string cimport memset

# bond:
# single, double, triple, aromatic, special = 5 bit
# o_bond & q_bond == o_bond

# atom:
# isotope: isotope - common_isotope -8 - +8 = 17 bit
# is_radical: 2 bit
# charge: -4 - +4: 9 bit
# implicit_hydrogens: 0-5 = 6 bit
# neighbors: 0-14 = 15 bit
# heteroatoms: 0-14 = 15 bit

# hybridizations: 1-4 = 4 bit
# ring_sizes: not-in-ring bit, 3-atom ring, 4-...., 60 atom ring
# exact match
# element: 7 bit


def get_mapping(unsigned long[:] q_numbers not None, unsigned int[:] q_elements not None,
                unsigned int[:] q_isotopes not None, unsigned long long[:] q_masks not None,
                unsigned long long[:] q_rings not None, unsigned int[:] q_bonds not None,
                unsigned int[:] q_back not None,
                #query_closures,
                unsigned long[:] o_numbers not None, unsigned int[:] o_elements not None,
                unsigned int[:] o_isotopes not None, unsigned long long[:] o_masks not None,
                unsigned long long[:] o_rings not None,
                unsigned int[:] o_shifts not None, unsigned int[:] o_neighbors not None,
                unsigned int[:] o_indices not None, unsigned int[:] o_bonds not None,
                bint[:] scope not None, int[:] groups not None):
    # expected less than 2^16 atoms in structure.
    cdef unsigned int stack = 0, path_size = 0, q_size, i, j, depth, front, back, q_element, q_bond, o_bond, q_isotope
    cdef unsigned long q_number, o_number
    cdef unsigned long long q_mask, q_ring
    cdef dict mapping

    q_size = len(q_elements) - 1
    cdef int *path = <int *> PyMem_Malloc(q_size * sizeof(int))
    cdef int *stack_index = <int *> PyMem_Malloc(2 * len(o_elements) * sizeof(int))
    cdef int *stack_depth = <int *> PyMem_Malloc(2 * len(o_elements) * sizeof(int))
    cdef bint *matched = <bint *> PyMem_Malloc(len(o_elements) * sizeof(bint))
    if not path or not stack_index or not stack_depth or not matched:
        raise MemoryError()
    memset(&matched[0], 0, len(o_elements) * sizeof(bint))

    q_element = q_elements[0]
    q_isotope = q_isotopes[0]
    q_mask = q_masks[0]
    q_ring = q_rings[0]
    for i in range(len(o_elements)):
        if scope[i] and q_element == o_elements[i] and q_mask & o_masks[i] == o_masks[i] and \
                q_ring & o_rings[i] == o_rings[i]:
            if q_isotope:
                if q_isotope != o_isotopes[i]:
                    continue
            stack_index[stack] = i
            stack_depth[stack] = 0
            stack += 1

    try:
        while stack:
            stack -= 1
            depth = stack_depth[stack]
            i = stack_index[stack]

            if depth == q_size:
                mapping = {q_numbers[depth]: o_numbers[i]}
                for i in range(q_size):
                    mapping[q_numbers[i]] = o_numbers[path[i]]
                yield mapping
            else:
                if path_size != depth:  # dead end reached
                    for j in range(depth, path_size):
                        matched[path[j]] = False  # mark unmatched
                    path_size = depth

                matched[i] = True
                path[path_size] = i
                path_size += 1

                front = depth + 1
                back = q_back[front]
                if back != depth:  # branch
                    i = path[back]

                # load next query atom
                q_element = q_elements[front]
                q_mask = q_masks[front]
                q_ring = q_rings[front]
                q_bond = q_bonds[front]

                for j in range(o_shifts[i], o_shifts[i] + o_neighbors[i]):
                    o_bond = o_bonds[j]
                    j = o_indices[j]  # now j is atom index
                    if scope[j] and not matched[j] and q_bond & o_bond == o_bond:  # and groups[j] not in uniq:
                        # uniq.add(groups[o_n])

                        if q_element == o_elements[j] and q_mask & o_masks[j] == o_masks[j] and \
                            q_ring & o_rings[j] == o_rings[j]:
                            # check closures equality
                            # o_closures = o_bonds[o_n].keys() & reversed_mapping.keys()
                            # o_closures.discard(n)
                            #if o_closures == {mapping[m] for m, _ in query_closures[s_n]}:
                            #    obon = o_bonds[o_n]
                            #    if all(bond == obon[mapping[m]] for m, bond in query_closures[s_n]):
                            stack_index[stack] = j
                            stack_depth[stack] = front
                            stack += 1
    finally:
        PyMem_Free(path)
        PyMem_Free(matched)
        PyMem_Free(stack_index)
        PyMem_Free(stack_depth)
