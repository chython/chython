# -*- coding: utf-8 -*-
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
from chython.core import MoleculeContainer


# C60-Ih edge list, 90 bonds over 60 vertices, generated from the truncated
# icosahedron: 12 pentagonal and 20 hexagonal faces.
C60_EDGES = [
    (0, 1), (0, 2), (0, 4), (1, 3), (1, 5), (2, 6), (2, 10), (3, 7), (3, 11),
    (4, 8), (4, 12), (5, 9), (5, 13), (6, 7), (6, 14), (7, 15), (8, 9),
    (8, 16), (9, 17), (10, 12), (10, 18), (11, 13), (11, 19), (12, 20),
    (13, 21), (14, 22), (14, 23), (15, 22), (15, 24), (16, 25), (16, 27),
    (17, 26), (17, 27), (18, 23), (18, 28), (19, 24), (19, 29), (20, 25),
    (20, 30), (21, 26), (21, 31), (22, 32), (23, 33), (24, 34), (25, 35),
    (26, 36), (27, 37), (28, 30), (28, 38), (29, 31), (29, 39), (30, 40),
    (31, 41), (32, 42), (32, 43), (33, 38), (33, 42), (34, 39), (34, 43),
    (35, 40), (35, 44), (36, 41), (36, 45), (37, 44), (37, 45), (38, 46),
    (39, 47), (40, 48), (41, 49), (42, 50), (43, 51), (44, 52), (45, 53),
    (46, 48), (46, 54), (47, 49), (47, 55), (48, 56), (49, 57), (50, 51),
    (50, 54), (51, 55), (52, 53), (52, 56), (53, 57), (54, 58), (55, 59),
    (56, 58), (57, 59), (58, 59)]


def test_c60_has_thirty_two_faces_at_rank_thirty_one():
    # 32 faces at circuit rank 31 is exactly why the per-atom descriptors are read off the
    # relevant-cycle prototypes rather than off `rings`: one real hexagonal face is the GF(2)
    # sum of the other 31 faces, so no cycle basis can contain all 32. Rank is not a
    # completeness test, and a descriptor filled from the basis would lose that face.
    m = MoleculeContainer()
    with m.edit():
        ids = [m.add_atom(6) for _ in range(60)]
        for i, j in C60_EDGES:
            m.add_bond(ids[i], ids[j], 1)
    assert m.bond_count == 90
    assert m.atom_count == 60
    assert all(m.atom(s).degree == 3 for s in ids)
    assert m.rings_count == 31           # circuit rank is 90 - 60 + 1 = 31
    # the 12 pentagons are edge-disjoint and therefore independent; one hexagon is redundant
    assert sorted(len(r) for r in m.rings) == [5] * 12 + [6] * 19

    # every vertex of C60 lies on one pentagon and two hexagons, and all 32 faces are relevant
    for s in ids:
        assert m.ring_sizes_of(s) == frozenset({5, 6})
        assert m.ring_count_of(s) == 3
