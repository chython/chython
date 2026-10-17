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
"""The ring-shape post-steps and the plane predicates the controller shares."""
from math import cos, hypot, pi, sin

import pytest

from chython.core import smiles
from chython.depict.layout.lattice import L, project
from chython.depict.layout.quality import _polygon_rmsd
from chython.depict.layout.relax import _arc, _ct_broken, _cycles, _frames, _near, count_defects, fair
from chython.depict.layout.rings import planar_block, ring_blocks


def _rmsd(plane, ring):
    return _polygon_rmsd([plane[x] for x in ring])


def test_a_substituent_on_the_double_bond_line_breaks_its_cis_trans():
    mol = smiles('C/C=C(/F)C')                              # C1 trans to F4, cis to C5
    good = {1: (-1., -1.), 2: (0., 0.), 3: (1., 0.), 4: (2., 1.), 5: (2., -1.)}
    frames = _frames(mol, good)
    assert len(frames) == 1
    assert _ct_broken(mol, good, frames, set(good)) == 0
    assert _ct_broken(mol, {**good, 4: (2., -1.), 5: (2., 1.)}, frames, set(good)) == 1
    assert _ct_broken(mol, {**good, 5: (2., 0.)}, frames, set(good)) == 1          # on the bond's line
    assert _ct_broken(mol, {**good, 5: (2.5, 1.5)}, frames, set(good)) == 1        # beside F


@pytest.mark.parametrize('s', ['C12C3C4C1C5C2C3C45', 'C1C2CC3CC1CC(C2)C3', 'CN1CCC23C4C1CC5=C2C(=C(C=C5)O)OC3C(C=C4)O'])
def test_a_planar_drawing_evens_its_bonds(s):
    mol = smiles(s)
    rank = mol.canonical_order()
    for atoms, rings in ring_blocks(mol):
        pts = {x: project(p) for x, p in planar_block(mol, atoms, rings, rank).items()}
        ls = [hypot(pts[b.n][0] - pts[b.m][0], pts[b.n][1] - pts[b.m][1]) for b in mol.bonds()
              if b.n in atoms and b.m in atoms]
        assert max(ls) / min(ls) < 3


def test_the_local_window_spans_a_moved_atom_s_long_bond():
    mol = smiles('CCCC')                                    # 1-2-3-4; 1 stretched along 2-3, over 3
    plane = {2: (0., 0.), 3: (L, 0.), 4: (1.5 * L, L), 1: (5 * L, 0.)}
    atoms, bonds = list(plane), [(b.n, b.m) for b in mol.bonds()]
    near, local = _near((plane,), atoms, bonds, {1})
    assert count_defects(plane, near, local, {1})[2] == 1


def test_a_spike_already_on_its_arc_is_left_alone():
    # cyclononane on an open circle, atom 0 folded in between its neighbours: the arc is where it already is
    mol = smiles('C1CCCCCCCC1')
    ring = list(mol)
    plane = {x: (1.5 * cos(.2 + (2 * pi - .4) * i / 7), 1.5 * sin(.2 + (2 * pi - .4) * i / 7))
             for i, x in enumerate(ring[1:])}
    p, q = plane[ring[8]], plane[ring[1]]
    plane[ring[0]] = (.5, 0.)
    for _ in range(9):
        side = (hypot(p[0] - plane[ring[0]][0], p[1] - plane[ring[0]][1]) +
                hypot(q[0] - plane[ring[0]][0], q[1] - plane[ring[0]][1])) / 2
        plane[ring[0]], = _arc(p, q, 2, side, False)
    old = dict(plane)
    assert fair(mol, plane, [tuple(ring)], {x: i for i, x in enumerate(ring)}) == 0
    assert plane == old


def test_the_rings_shorten_keeps_clear_are_every_small_cycle_not_one_sssr():
    mol = smiles('C1CC2CCC1C2')                             # norbornane: an SSSR holds two of its three rings
    cycles = _cycles(mol, dict.fromkeys(mol, (0., 0.)), mol.canonical_order())
    assert sorted(len(c) for c in cycles) == [5, 5, 6]
