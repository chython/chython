#
#  Copyright 2026 Ramil Nugmanov <nougmanoff@protonmail.com>
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
# The dict-graph entry point into core ring perception, for a caller that holds an adjacency dict
# and not a `MoleculeContainer` -- which reads `MoleculeContainer.rings` instead.
#
# `_rings.pxi` ends at a minimum cycle basis: `_select_basis` reaches full rank over the relevant
# cycles taking the shortest independent ones first. A candidate-ring-set filter of the CSET kind is
# not an equivalent answer -- it can return a set that is not a cycle basis at all, which
# test_a_cage_whose_smallest_ring_set_is_not_a_cycle_basis pins down.
#
# So this is an adapter and no second algorithm: it takes a skin graph -- atom number to neighbour
# set, bridges already stripped by the caller -- stages it in an arena, runs perception, and hands
# SEG_RELEVANT_RINGS back under the caller's own atom numbers.


def sssr(object graph not None):
    """
    A minimum cycle basis of `graph`, given as {atom number: neighbour atom numbers}.

    :return: list of tuples of atom numbers, each in cyclic order, shortest rings first
    """
    cdef uint32_t n = len(graph)
    cdef list labels = list(graph)
    cdef dict index_of = {}
    cdef list bonds = []
    cdef uint32_t i, j, k, count, end
    cdef object a, b
    cdef Structure structure
    cdef uint32_t *r
    cdef list out = []
    cdef list row

    if not n:
        return out
    for i in range(n):
        index_of[labels[i]] = i
    for a in graph:
        i = <uint32_t> index_of[a]
        for b in graph[a]:
            j = <uint32_t> index_of[b]
            if i < j:
                bonds.append((i, j, 1))
    if not bonds:
        return out

    structure = _build_csr_from_list(n, bonds)
    with nogil:
        mark_bridges(structure)
    perceive_rings(structure)

    # SEG_RELEVANT_RINGS is [count][offset per ring][total][atom indices]
    r = structure_rings(structure)
    count = r[0]
    for i in range(count):
        end = r[2 + i]
        row = []
        for k in range(r[1 + i], end):
            row.append(labels[r[2 + count + k]])
        out.append(tuple(row))
    return out
