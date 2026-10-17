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
"""Lattice predicates and ring blocks: the geometry the controller's embed variants come from."""
from math import atan2, degrees, hypot

import pytest

from chython.core import smiles
from chython.depict.layout.lattice import (DIRS, L, ON_BOND, SYMMETRIES, Board, _seg_dist, apply, direction, project,
                                          reflect, rotate)
from chython.depict.layout.quality import _polygon_rmsd
from chython.depict.layout.rings import MOVES, NEAR, embed_block, ring_blocks
from chython.depict.layout.templates import template_block


# --- lattice ---------------------------------------------------------------------------------------------

@pytest.mark.parametrize('k', range(6))
def test_each_direction_is_one_bond_at_its_angle(k):
    x, y = project(DIRS[k])
    assert hypot(x, y) == pytest.approx(L)
    assert degrees(atan2(y, x)) % 360 == pytest.approx(30 + 60 * k)


def test_rotation_and_mirror_permute_the_directions():
    assert [rotate(d, 1) for d in DIRS] == list(DIRS[1:] + DIRS[:1])
    assert [reflect(d) for d in DIRS] == list(DIRS[::-1])
    assert len({tuple(apply(d, g) for d in DIRS) for g in SYMMETRIES}) == 12


def test_direction_reads_straight_steps_only():
    assert direction((3, 0)) == (0, 3)
    assert direction((-2, 2)) == (2, 2)
    assert direction((0, -1)) == (4, 1)
    assert direction((1, 1)) is None


def test_board_refuses_an_occupied_point_a_bond_through_an_atom_and_a_crossing():
    b = Board()
    b.put(0, (0, 0))
    b.put(1, (1, 0))
    b.link(0, 1)
    assert b.bad_atom((0, 0)) == 1
    assert b.bad_bond((-1, 0), (2, 0)) == 2                  # a double step through both atoms
    b.put(2, (0, 1))
    b.put(3, (1, -1))
    b.link(2, 3)                                            # crosses 0-1 at its middle
    assert b.bad_bond((0, 1), (1, -1)) >= 1
    assert b.bad_bond((0, 0), (0, 1)) == 0                  # a unit edge sharing an end point


# --- ring blocks -----------------------------------------------------------------------------------------

def _embed(s):
    mol = smiles(s)
    rank = mol.canonical_order()
    (atoms, rings), = ring_blocks(mol)
    return mol, embed_block(mol, atoms, rings, rank)


@pytest.mark.parametrize('s', ['c1ccccc1', 'c1ccc2ccccc2c1', 'c1ccc2[nH]ccc2c1', 'C1CC1', 'C1CC2CCC1C2',
                               'C1COCCOCCOCCOCCOCCO1', 'c1cc2ccc3cccc4ccc(c1)c2c34', 'C1CCCCCCCCCCC1'])
def test_a_block_embeds_on_the_lattice_without_conflict(s):
    mol, (pts, conflicts) = _embed(s)
    assert conflicts == 0
    assert all(isinstance(v, int) for p in pts.values() for v in p)
    assert len(set(pts.values())) == len(pts)


def test_benzene_is_a_lattice_hexagon():
    mol, (pts, _) = _embed('c1ccccc1')
    for b in mol.bonds():
        d = (pts[b.n][0] - pts[b.m][0], pts[b.n][1] - pts[b.m][1])
        assert d in DIRS or (-d[0], -d[1]) in DIRS


def test_the_near_points_of_a_move_answer_on_bond_with_margin():
    # an ear walk tests ON_BOND by lookup: every lattice point off NEAR is well clear of the threshold
    for (v, _), near in zip(MOVES, NEAR):
        for q in range(-6, 7):
            for r in range(-6, 7):
                if (q, r) in ((0, 0), v):
                    continue
                d = _seg_dist(project((q, r)), (0., 0.), project(v))
                assert (d < ON_BOND - .01 * L) if (q, r) in near else (d > ON_BOND + .01 * L), (v, (q, r), d)


def test_cyclopropane_is_an_equilateral_lattice_triangle():
    mol, (pts, _) = _embed('C1CC1')
    xy = [project(p) for p in pts.values()]
    sides = {round(hypot(a[0] - b[0], a[1] - b[1]), 9) for a in xy for b in xy if a != b}
    assert sides == {round(L, 9)}


def test_a_cage_basis_ignores_atom_order():
    # cubane has six faces and any five are a minimum cycle basis
    mol = smiles('C12C3C4C1C5C2C3C45')
    ref = sorted(sorted(mol.canonical_order()[x] for x in r) for r in ring_blocks(mol)[0][1])
    for _ in range(5):
        other = smiles(format(mol, 'r'))
        rank = other.canonical_order()
        assert sorted(sorted(rank[x] for x in r) for r in ring_blocks(other)[0][1]) == ref


def test_blocks_and_rings_come_in_rank_order():
    mol = smiles('C1CC2(CO1)CCOC2')
    rank = mol.canonical_order()
    blocks = ring_blocks(mol)
    assert [sorted(rank[x] for x in a) for a, _ in blocks] == sorted(sorted(rank[x] for x in a) for a, _ in blocks)
    for _, rings in blocks:
        for r in rings:
            assert rank[r[0]] == min(rank[x] for x in r) and rank[r[1]] < rank[r[-1]]


# --- assembly --------------------------------------------------------------------------------------------


# --- linear atoms ----------------------------------------------------------------------------------------


# --- crosses ---------------------------------------------------------------------------------------------


# --- convex small rings ----------------------------------------------------------------------------------


# --- macrocycles of rigid units --------------------------------------------------------------------------


# PubChem 490721: a rifamycin, its naphthofuranone spanned by a 17-atom ansa chain
RIFAMYCIN = ('CCOC1=CC=C(C=C1)NC=C2C(=C3C4=C(C2=NC(=O)C(=CC=C[C@@H]([C@@H]([C@H]([C@H]([C@H]([C@@H]([C@@H]([C@H](C=CO'
             '[C@@]5(C(=O)C3=C(O5)C(=C4O)C)C)OC)C)OC(=O)C)C)O)C)O)C)C)O)O')


def test_a_unit_of_a_macrocycle_keeps_its_honeycomb():
    mol = smiles(RIFAMYCIN)
    rank = mol.canonical_order()
    atoms, rings = next(b for b in ring_blocks(mol) if max(map(len, b[1])) > 8)
    pts, bad = embed_block(mol, atoms, rings, rank)
    assert not bad
    assert all(_polygon_rmsd([project(pts[x]) for x in r]) < 1e-9 for r in rings if len(r) == 6)


# PubChem 5284616
RAPAMYCIN = ('C[C@@H]1CC[C@H]2C[C@@H](/C(=C/C=C/C=C/[C@H](C[C@H](C(=O)[C@@H]([C@@H](/C(=C/[C@H](C(=O)C[C@H](OC(=O)'
             '[C@@H]3CCCCN3C(=O)C(=O)[C@@]1(O2)O)[C@H](C)C[C@@H]4CC[C@H]([C@@H](C4)OC)O)C)/C)O)OC)C)C)/C)OC')


def test_a_long_ear_closes_on_the_lattice():
    mol = smiles(RAPAMYCIN)
    rank = mol.canonical_order()
    atoms, rings = max(ring_blocks(mol), key=lambda b: len(b[0]))
    pts, conflicts = embed_block(mol, atoms, rings, rank)
    assert conflicts == 0
    assert all(type(v) is int for p in pts.values() for v in p)


# --- spheres ---------------------------------------------------------------------------------------------


# --- templates -------------------------------------------------------------------------------------------


def test_a_block_beyond_a_template_is_not_matched():
    # PubChem 644073: buprenorphine, an epoxymorphinan with a 6,14-ethano bridge
    mol = smiles('CC(C)(C)C(C)(C1CC23CCC1(C4C25CCN(C3CC6=C5C(=C(C=C6)O)O4)CC7CC7)OC)O')
    rank = mol.canonical_order()
    block = max(ring_blocks(mol), key=lambda b: len(b[0]))
    assert template_block(mol, *block, rank) is None


def test_a_block_offers_its_best_distinct_polyhexes():
    mol = smiles('c1ccc2c(c1)ccc1ccccc12')                # phenanthrene: an angular and a mirrored angular polyhex
    rank = mol.canonical_order()
    (atoms, rings), = ring_blocks(mol)
    out = embed_block(mol, atoms, rings, rank, k=4)
    assert isinstance(out, list) and len(out) >= 2
    assert sum(not bad for _, bad in out) >= 2 and [bad for _, bad in out] == sorted(bad for _, bad in out)
    assert len({tuple(sorted(pts.items())) for pts, _ in out}) == len(out)
    assert out[0] == embed_block(mol, atoms, rings, rank)
