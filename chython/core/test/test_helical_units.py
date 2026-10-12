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
"""Helicenes: one `SU_HELICAL` unit per ortho-fused chain of five or more angular six-rings.

| field  | value                                                                                    |
| ------ | ---------------------------------------------------------------------------------------- |
| owners | the two fjord atoms, the terminal rings' inner-rim CH atoms; not bonded                   |
| anchor | the lower-slot fjord atom, or the other one where that slot is taken                     |
| refs   | (A's inner-rim neighbour, A's other ring neighbour, B's inner-rim, B's other)             |
| parity | the signed volume of the refs in that order: 1 is the M helix, 2 the P helix              |

The P reference is the conformer below, whose inner-rim torsions are all positive (a right-handed
helix); the Blue Book pair VS010 (P) / VS011 (M) is read in ``chython/test/test_stereo_bluebook.py``.
Compounds: hexahelicene (PubChem CID 98863), pentahelicene (CID 67444), heptahelicene (CID 140151),
2-bromohexahelicene (CID 5256406), benzo[c]phenanthrene (CID 9136), pentacene (CID 8671), picene.
"""

import math
from itertools import product

from pytest import mark

from .. import read_smiles, write_smiles
from ..wedge import assign_parities


SU_ATROPISOMER = 3
SU_HELICAL = 4
HELICENES = {
    '[5]helicene': 'c1ccc2c(c1)ccc1ccc3ccc4ccccc4c3c12',
    '[6]helicene': 'c1ccc2c(c1)ccc1ccc3ccc4ccc5ccccc5c4c3c12',
    '[7]helicene': 'c1ccc2c(c1)ccc1ccc3ccc4ccc5ccc6ccccc6c5c4c3c12',
    '2-bromo[6]helicene': 'Brc1ccc2ccc3ccc4ccc5ccc6ccccc6c5c4c3c2c1',
}
NOT_HELICAL = {
    'benzo[c]phenanthrene': 'c1ccc2c(c1)ccc1ccc3ccccc3c12',
    'pentacene': 'c1ccc2cc3cc4cc5ccccc5cc4cc3cc2c1',
    'picene': 'c1ccc2c(c1)ccc1c2ccc2c3ccccc3ccc21',
}
H6 = HELICENES['[6]helicene']
# ETKDGv3 + MMFF heavy-atom conformer of H6, in its atom order
P_H6 = [(-3.36, -0.28, -1.87), (-2.37, -1.24, -1.72), (-1.33, -1.02, -0.82), (-1.2, 0.17, -0.06),
        (-2.28, 1.08, -0.17), (-3.33, 0.86, -1.09), (-2.34, 2.19, 0.67), (-1.4, 2.35, 1.67),
        (-0.33, 1.47, 1.8), (0.49, 1.6, 2.92), (1.51, 0.69, 3.15), (1.84, -0.26, 2.18),
        (2.84, -1.18, 2.48), (3.24, -2.11, 1.54), (2.79, -2.02, 0.23), (3.31, -2.92, -0.71),
        (3.02, -2.79, -2.06), (2.28, -1.7, -2.51), (2.12, -1.52, -3.9), (1.49, -0.39, -4.41),
        (1.03, 0.58, -3.53), (1.18, 0.4, -2.15), (1.75, -0.77, -1.58), (1.87, -1.01, -0.15),
        (1.18, -0.28, 0.92), (-0.1, 0.44, 0.85)]
RIM = (3, 4, 26, 25, 24, 23, 22)                # fjord atom to fjord atom along the inner rim of H6
BIHELICENYL = 'c1cc(-c2cccc3c2c2c(cc3)ccc3ccc4ccc5ccccc5c4c32)c2c(c1)ccc1ccc3ccc4ccc5ccccc5c4c3c21'


def _helices(mol):
    return [u for u in mol.stereo_units() if u['kind'] == SU_HELICAL]


def _dihedral(p0, p1, p2, p3):
    b0 = [p0[i] - p1[i] for i in range(3)]
    b1 = [p2[i] - p1[i] for i in range(3)]
    b2 = [p3[i] - p2[i] for i in range(3)]
    norm = math.sqrt(sum(x * x for x in b1))
    b1 = [x / norm for x in b1]
    v = [b0[i] - sum(x * y for x, y in zip(b0, b1)) * b1[i] for i in range(3)]
    w = [b2[i] - sum(x * y for x, y in zip(b2, b1)) * b1[i] for i in range(3)]
    c = (b1[1] * v[2] - b1[2] * v[1], b1[2] * v[0] - b1[0] * v[2], b1[0] * v[1] - b1[1] * v[0])
    return math.degrees(math.atan2(sum(x * y for x, y in zip(c, w)), sum(x * y for x, y in zip(v, w))))


def _from_conformer(mirror):
    mol = read_smiles(H6)
    z = {}
    with mol.edit():
        for n, (x, y, h) in zip(mol.atom_numbers, P_H6):
            mol.set_xy(n, x, y)
            z[n] = -h if mirror else h
    log = []
    assign_parities(mol, z=z, log=log)
    assert log == [], log
    return mol


def _configured(smiles, parity):
    mol = read_smiles(smiles)
    (unit,) = _helices(mol)
    mol.set_parity(unit['anchor'], parity)
    return mol


@mark.parametrize('name', sorted(HELICENES))
def test_a_helicene_has_exactly_one_stereogenic_helix(name):
    mol = read_smiles(HELICENES[name])
    (unit,) = _helices(mol)
    assert unit['stereogenic']
    a, b = unit['owners']
    assert unit['anchor'] in (a, b) and b not in mol.neighbors_of(a)


@mark.parametrize('name', sorted(NOT_HELICAL))
def test_four_angular_rings_a_linear_acene_and_a_broken_rim_have_no_helix(name):
    assert _helices(read_smiles(NOT_HELICAL[name])) == []


def test_the_reference_conformer_is_a_right_handed_helix():
    torsions = [_dihedral(*[P_H6[n - 1] for n in RIM[k:k + 4]]) for k in range(len(RIM) - 3)]
    assert all(t > 0 for t in torsions), torsions


def test_the_p_helix_is_parity_2_and_its_mirror_image_is_parity_1():
    p, m = _from_conformer(False), _from_conformer(True)
    (up,), (um,) = _helices(p), _helices(m)
    assert (p.parity_of(up['anchor']), m.parity_of(um['anchor'])) == (2, 1)
    assert p != m


@mark.parametrize('name', sorted(HELICENES))
@mark.parametrize('parity', (1, 2))
def test_a_helix_survives_a_cx_round_trip(name, parity):
    mol = _configured(HELICENES[name], parity)
    log = []
    text = write_smiles(mol, log=log)
    assert ' |w' in text and log == [], (text, log)
    back = read_smiles(text)
    assert list(back.log) == [], list(back.log)
    assert back == mol
    assert back != _configured(HELICENES[name], 3 - parity)


def test_plain_smiles_logs_the_helix_it_cannot_write():
    log = []
    text = write_smiles(_configured(H6, 2), '!x', log=log)
    assert '|' not in text
    assert [r.rule for r in log] == ['smiles:stereo-not-written'], log


@mark.parametrize('parity', (1, 2))
def test_a_helix_survives_pach(parity):
    mol = _configured(H6, parity)
    back = type(mol).unpack(mol.pack())
    assert back == mol
    (unit,) = _helices(back)
    assert back.parity_of(unit['anchor']) == parity


def test_inchi_omits_the_helix_and_logs_it_once():
    mol = _configured(H6, 2)
    first = mol.inchi
    assert first == read_smiles(H6).inchi
    assert mol.inchi == first
    rules = [r.rule for r in mol.log]
    assert rules.count('inchi:stereo-kind-not-represented') == 1, rules


def test_two_helices_and_the_axis_between_them_are_independent():
    """1,1'-bi[6]helicenyl: the axis joins the two fjord atoms, so each helix anchors at its other one.

    The C2 exchange of the halves keeps the axis and swaps the helices: 8 labellings, 6 compounds.
    """
    mol = read_smiles(BIHELICENYL)
    units = [u for u in mol.stereo_units() if u['stereogenic']]
    assert sorted(u['kind'] for u in units) == [SU_ATROPISOMER, SU_HELICAL, SU_HELICAL]
    anchors = [u['anchor'] for u in units]
    assert len(set(anchors)) == 3
    distinct = []
    for parities in product((1, 2), repeat=3):
        m = read_smiles(BIHELICENYL)
        for anchor, parity in zip(anchors, parities):
            m.set_parity(anchor, parity)
        assert [m.parity_of(a) for a in anchors] == list(parities)
        back = read_smiles(write_smiles(m))
        assert back == m, parities
        if not any(m == d for d in distinct):
            distinct.append(m)
    assert len(distinct) == 6
