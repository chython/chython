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
"""A nucleic-acid strand laid out level: every unit the same drawing, its base up, no link long."""
from math import hypot

import pytest

from chython.core import smiles
from chython.core.nucleotides import nucleotides
from chython.core.test.test_nucleotides import STRANDS
from chython.depict.layout.lattice import L
from chython.depict.layout.molecule import layout2d
from chython.depict.layout.nucleic import PERIOD


@pytest.mark.parametrize('name', sorted(STRANDS))
def test_a_strand_lies_level_with_its_bases_up(name):
    mol = smiles(STRANDS[name])
    strand = nucleotides(mol)
    plane = layout2d(mol, engine='chython', force=True, chain=True)
    ps = [plane[r.path[0]] for r in strand.residues[1:]]
    assert max(y for _, y in ps) - min(y for _, y in ps) < 1e-6
    assert all(abs(q[0] - p[0] - PERIOD * L) < 1e-6 for p, q in zip(ps, ps[1:]))
    lengths = [hypot(plane[b.n][0] - plane[b.m][0], plane[b.n][1] - plane[b.m][1]) for b in mol.bonds()]
    assert max(lengths) < 1.05 * L and min(lengths) > .95 * L
    for r in strand.residues:
        c4 = r.path[-3]
        (c1,) = {y for x in r.atoms if mol.atom(x).atomic_symbol == 'O' and x not in r.path
                 for y in mol.neighbors_of(x) if y in r.atoms and y != c4}
        (n,) = [y for y in mol.neighbors_of(c1) if mol.atom(y).atomic_symbol == 'N']
        assert plane[n][1] > plane[c1][1] > plane[c4][1]                 # the base above its sugar


def test_every_unit_is_the_same_drawing():
    mol = smiles(STRANDS['dT4'])
    strand = nucleotides(mol)
    plane = layout2d(mol, engine='chython', force=True, chain=True)
    shapes = []
    for r in strand.residues[1:]:
        p0 = plane[r.path[0]]
        shapes.append([(round(plane[x][0] - p0[0], 6), round(plane[x][1] - p0[1], 6)) for x in r.path])
    assert all(s == shapes[0] for s in shapes)


@pytest.mark.parametrize('name', sorted(STRANDS))
def test_every_phosphorus_is_a_square_cross(name):
    """O3' and O5' opposite on the link's axis, the other two atoms square to it: four gaps of 90 degrees."""
    from math import atan2, degrees
    mol = smiles(STRANDS[name])
    plane = layout2d(mol, engine='chython', force=True, chain=True)
    for p in (x for x in mol if mol.atom(x).atomic_symbol == 'P'):
        t = sorted(degrees(atan2(plane[y][1] - plane[p][1], plane[y][0] - plane[p][0])) for y in mol.neighbors_of(p))
        assert [round((b - a) % 360) for a, b in zip(t, t[1:] + t[:1])] == [90] * 4
