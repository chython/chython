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


# bond:
# single, double, triple, aromatic, special = 5 bit
# o_bond & q_bond == o_bond

# atom:
# masked match
# neighbors: 0-14
# implicit_hydrogens: 0-14
# heteroatoms: 0-14
# hybridizations: 1-4
# is_radical: 2 bit
# charge: -4 - +4: 9 bit
# ring_sizes: not-in-ring bit, 3-atom ring, 4-...., 63+ atom ring
# exact match
# element: 7 bit
# isotope: 9 bit (0 - not specified)


def get_mapping(unsigned long[:] q_numbers not None, unsigned int[:] q_elements not None,
                unsigned long long[:] q_masks not None, unsigned long long[:] q_rings not None,
                unsigned int[:] q_bonds not None,
                unsigned long[:] q_back not None,  # ?
                #query_closures,
                unsigned long[:] o_numbers not None, unsigned int[:] o_elements not None,
                unsigned long long[:] o_masks not None, unsigned long long[:] o_rings not None,
                # o_bonds,
                bint[:] scope, int[:] groups):
    # expected less than 2^16 atoms in structure.
    cdef unsigned int stack = 0, path_size = 0, q_size, i, depth, front, back, q_element, q_bond
    cdef unsigned long q_number, o_number
    cdef unsigned long long q_mask, q_ring
    cdef dict mapping

    q_size = len(q_numbers) - 1
    cdef int *path = <int *> PyMem_Malloc(q_size * sizeof(int))
    cdef int *stack_index = <int *> PyMem_Malloc(2 * len(o_numbers) * sizeof(int))
    cdef int *stack_depth = <int *> PyMem_Malloc(2 * len(o_numbers) * sizeof(int))
    if not path or not stack_index or not stack_depth:
        raise MemoryError()

    q_element = q_elements[0]
    q_mask = q_masks[0]
    q_ring = q_rings[0]
    for i in range(len(o_numbers)):
        if scope[i] and q_element == o_elements[i] and q_mask & o_masks[i] == o_masks[i] and \
                q_ring & o_rings[i] == o_rings[i]:
            stack_index[stack] = i
            stack_depth[stack] = 0
            stack += 1

    while stack:
        stack -= 1
        depth = stack_depth[stack]

        if depth == q_size:
            mapping = {q_numbers[depth]: o_numbers[stack_index[stack]]}
            for i in range(q_size):
                mapping[q_numbers[i]] = o_numbers[path[i]]
            yield mapping
        else:
            i = stack_index[stack]
            if path_size != depth:  # dead end reached
                path_size = depth  # reset to

            path[path_size] = i
            path_size += 1

            #s_n, back, s_atom, s_bond = linear_query[depth]

            front = depth + 1
            back = q_back[front]
            if back != depth:  # branch
                i = path[back]

            uniq = set()
            for o_n, o_bond in o_bonds[n].items():
                if o_n in scope and o_n not in reversed_mapping and s_bond == o_bond and groups[o_n] not in uniq:
                    uniq.add(groups[o_n])
                    if s_atom == o_atoms[o_n]:
                        # check closures equality
                        o_closures = o_bonds[o_n].keys() & reversed_mapping.keys()
                        o_closures.discard(n)
                        if o_closures == {mapping[m] for m, _ in query_closures[s_n]}:
                            obon = o_bonds[o_n]
                            if all(bond == obon[mapping[m]] for m, bond in query_closures[s_n]):
                                stack.append((o_n, front))

    PyMem_Free(path)
    PyMem_Free(stack_index)
    PyMem_Free(stack_depth)
