# -*- coding: utf-8 -*-
#
#  Copyright 2021-2025 Ramil Nugmanov <nougmanoff@protonmail.com>
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

cdef packed struct atom_t:
    unsigned long long bits1
    unsigned long long bits2
    unsigned long long bits3
    unsigned long long bits4
    unsigned int from_
    unsigned int to_
    unsigned int mapping

cdef packed struct bond_t:
    unsigned long long bond
    unsigned int index

cdef packed struct molecule_t:
    unsigned int atoms_count
    atom_t *atoms
    bond_t *bonds

cdef packed struct q_atom_t:
    unsigned long long mask1
    unsigned long long mask2
    unsigned long long mask3
    unsigned long long mask4
    unsigned int back
    unsigned int closure
    unsigned int from_
    unsigned int to_
    unsigned int mapping

cdef packed struct query_t:
    unsigned int atoms_count
    q_atom_t *atoms
    bond_t *bonds


@cython.boundscheck(False)
@cython.wraparound(False)
def get_mapping(const unsigned char[::1] q_buffer not None, const unsigned char[::1] m_buffer not None,
                const unsigned int[::1] scope not None):
    # expected less than 2^16 atoms in structure.
    cdef unsigned int stack = 0, path_size = 0, depth, front, q_decrement
    cdef unsigned int n, m, i, j, closures_counter
    cdef unsigned long long c_bond
    cdef dict mapping
    cdef query_t query
    cdef molecule_t molecule
    cdef q_atom_t q_atom
    cdef atom_t n_atom, m_atom
    cdef bond_t i_bond, j_bond

    query.atoms_count = (<unsigned int*> &q_buffer[0])[0]
    query.atoms = <q_atom_t*> (&q_buffer[0] + 4)
    query.bonds = <bond_t*> (&query.atoms[0] + query.atoms_count)
    q_decrement = query.atoms_count - 1

    molecule.atoms_count = (<unsigned int*> &m_buffer[0])[0]
    molecule.atoms = <atom_t*> (&m_buffer[0] + 4)
    molecule.bonds = <bond_t*> (&molecule.atoms[0] + molecule.atoms_count)

    cdef unsigned int *path = <unsigned int *> PyMem_Malloc(q_decrement * sizeof(unsigned int))
    cdef unsigned int *stack_index = <unsigned int *> PyMem_Malloc(2 * molecule.atoms_count * sizeof(unsigned int))
    cdef unsigned int *stack_depth = <unsigned int *> PyMem_Malloc(2 * molecule.atoms_count * sizeof(unsigned int))
    cdef bint *matched = <bint *> PyMem_Malloc(molecule.atoms_count * sizeof(bint))
    cdef unsigned long long *closures = <unsigned long long *> PyMem_Malloc(molecule.atoms_count * sizeof(unsigned long long))

    if not path or not stack_index or not stack_depth or not matched or not closures:
        raise MemoryError()

    memset(matched, 0, molecule.atoms_count * sizeof(bint))
    memset(closures, 0, molecule.atoms_count * sizeof(unsigned long long))

    q_atom = query.atoms[0]
    for n in range(molecule.atoms_count):
        n_atom = molecule.atoms[n]
        if (scope[n] and
            q_atom.mask1 & n_atom.bits1 and  # bits1 doesn't contain bond bits.
            q_atom.mask2 & n_atom.bits2 == n_atom.bits2 and
            q_atom.mask3 & n_atom.bits3 == n_atom.bits3 and
            q_atom.mask4 & n_atom.bits4):

            stack_index[stack] = n
            stack_depth[stack] = 0
            stack += 1

    try:
        while stack:
            stack -= 1
            depth = stack_depth[stack]
            n = stack_index[stack]

            if depth == q_decrement:
                mapping = _PyDict_NewPresized(query.atoms_count)
                for i in range(depth):
                    mapping[query.atoms[i].mapping] = molecule.atoms[path[i]].mapping
                mapping[query.atoms[depth].mapping] = molecule.atoms[n].mapping
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
                # load next query atom
                q_atom = query.atoms[front]
                if q_atom.back != depth:  # branch
                    n = path[q_atom.back]

                n_atom = molecule.atoms[n]
                for i in range(n_atom.from_, n_atom.to_):
                    i_bond = molecule.bonds[i]
                    m = i_bond.index
                    m_atom = molecule.atoms[m]
                    if (scope[m] and not matched[m] and
                        q_atom.mask1 & i_bond.bond == i_bond.bond and  # bond order, in ring mark and atom bit should match.
                        q_atom.mask2 & m_atom.bits2 == m_atom.bits2 and
                        q_atom.mask3 & m_atom.bits3 == m_atom.bits3 and
                        q_atom.mask4 & m_atom.bits4):

                        if q_atom.closure:  # candidate atom should have same closures.
                            closures_counter = 0
                            # make a map of closures for o_n atom
                            # an index is a neighbor atom and a value is a bond between o_n and the neighbor
                            for j in range(m_atom.from_, m_atom.to_):
                                j_bond = molecule.bonds[j]
                                if j_bond.index != n and matched[j_bond.index]:
                                    closures[j_bond.index] = j_bond.bond
                                    closures_counter += 1

                            if closures_counter == q_atom.closure:
                                for j in range(q_atom.from_, q_atom.to_):
                                    j_bond = query.bonds[j]
                                    c_bond = closures[path[j_bond.index]]
                                    if not c_bond or j_bond.bond & c_bond != c_bond:  # compare order and ring bits
                                        break
                                else:
                                    stack_index[stack] = m
                                    stack_depth[stack] = front
                                    stack += 1

                            # fill an array with nulls
                            for j in range(m_atom.from_, m_atom.to_):
                                j_bond = molecule.bonds[j]
                                closures[j_bond.index] = 0
                        else:  # candidate atom should not have closures.
                            for j in range(m_atom.from_, m_atom.to_):
                               j_bond = molecule.bonds[j]
                               if j_bond.index != n and matched[j_bond.index]:
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
        PyMem_Free(closures)
