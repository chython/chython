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
"""Planar chirality: one `SU_PLANAR` unit per substituted ring atom that a second layer screens.

| shape       | anchor                                    | refs                                              |
| ----------- | ----------------------------------------- | ------------------------------------------------- |
| metallocene | a substituted η5 ring atom                | substituent, metal, both ring neighbours          |
| cyclophane  | a substituted decked ring atom            | substituent, bridge atom, bridgehead, other one   |
| twisted     | the later C=C terminal in an 8- or 9-ring | its ring neighbour and H, then its partner's      |

Parity is the signed volume of the refs in that order.  A twisted unit is dead unless its double bond
is trans.  chython has no planar CIP descriptor, so each reference states its sign, not an Rp/Sp label.

| compound                                  | CID             | reference conformer                            |
| ----------------------------------------- | --------------- | ---------------------------------------------- |
| 4-carboxy[2.2]paracyclophane              | 4614666         | PubChem 3-D, heavy atoms                       |
| 4,16- and 4,12-dibromo[2.2]paracyclophane | 5200652         | PubChem 3-D of CID 74210, Br on the replaced H |
| (E)-cyclooctene                           | 5463599         | RDKit ETKDG (E)-cyclononene less one CH2, MMFF |
| PPFA                                      | 129629204       | by hand: η5 rings 1.65 Å from Fe               |
| [2.2]paracyclophane                       | 74210           | none                                           |
| ferrocene, acetylferrocene                | 7611, 228207    | none                                           |
| (Z)-cyclooctene, (E)-cyclodecene          | 638079, 5364362 | none                                           |

RDKit ETKDG embeds no (E)-cyclooctene; PubChem's CID 5463599 conformer has a -1 degree C=C torsion, so
it reads as cis.  PPFA is 1-[2-(diphenylphosphino)ferrocenyl]-N,N-dimethylethylamine.
"""

import math

from pytest import mark

from .. import read_smiles, write_smiles
from ..wedge import assign_parities


SU_CIS_TRANS = 1
SU_PLANAR = 5
WEDGE_UP = 1
WEDGE_DOWN = 2

PCP_COOH = 'OC(=O)c1cc2ccc1CCc1ccc(cc1)CC2'
PCP_COOH_XYZ = [(2.12, -0.89, 2.06), (2.47, -0.76, 0.75), (3.63, -0.53, 0.4), (1.34, -0.94, -0.15),
                (0.3, -1.78, 0.26), (-0.9, -1.81, -0.47), (-0.92, -1.27, -1.76), (0.11, -0.43, -2.16),
                (1.16, -0.13, -1.27), (1.85, 1.22, -1.36), (1.37, 2.27, -0.29), (-0.04, 2.03, 0.24),
                (-0.21, 1.41, 1.49), (-1.33, 0.61, 1.72), (-2.28, 0.46, 0.7), (-2.26, 1.35, -0.39),
                (-1.14, 2.14, -0.62), (-3.07, -0.84, 0.59), (-2.2, -2.11, 0.25)]
DIBROMO = 'Brc1cc2ccc1CCc1ccc(c(Br)c1)CC2'  # each Br ortho to a different bridge
_DECK = [(-1.39, 1.54, 2.06), (-0.79, 1.62, 1.15), (0.6, 1.62, 1.25), (1.39, 1.45, 0.1), (0.79, 1.62, -1.15),
         (-0.6, 1.62, -1.25), (-1.39, 1.45, -0.1), (-2.74, 0.76, -0.22), (-2.74, -0.76, 0.22), (-1.39, -1.45, 0.1)]
PSEUDO_PARA_XYZ = _DECK + [(-0.6, -1.62, 1.25), (0.79, -1.62, 1.15), (1.39, -1.45, -0.1), (0.6, -1.62, -1.25),
                           (1.05, -1.54, -2.24), (-0.79, -1.62, -1.15), (2.74, -0.76, -0.22), (2.74, 0.76, 0.22)]
PSEUDO_ORTHO_XYZ = _DECK + [(-0.79, -1.62, -1.15), (0.6, -1.62, -1.25), (1.39, -1.45, -0.1), (0.79, -1.62, 1.15),
                            (1.39, -1.55, 2.06), (-0.6, -1.62, 1.25), (2.74, -0.76, -0.22), (2.74, 0.76, 0.22)]
TCO = 'C1CCC/C=C/CC1'
TCO_XYZ = [(-0.5635, -0.6661, 0.4763), (-0.6378, 0.4492, -0.5998), (0.2451, 1.7007, -0.3766),
           (1.6704, 1.6378, -0.9693), (2.2978, 0.3079, -0.7011), (2.3312, -0.1538, 0.5592),
           (1.8791, -1.5525, 0.8308), (0.4695, -1.7957, 0.2471)]
# PPFA: dative, ten Fe-C bonds, separate ions
PPFA = ('CN(C)C(C)c%30~%11c(P(c1ccccc1)c1ccccc1)~%12c~%13c~%14c~%15%30.'
        'c%31~%21c~%22c~%23c~%24c~%25%31.[Fe]%11%12%13%14%15%21%22%23%24%25')
PPFA_BONDED = ('CN(C)C(C)C%30%11=C(P(c1ccccc1)c1ccccc1)%12C%13=C%14C%15%30.'
               'C%31%21=C%22C%23=C%24C%25%31.[Fe]%11%12%13%14%15%21%22%23%24%25')
PPFA_IONIC = 'CN(C)C(C)[c-]1c(P(c2ccccc2)c2ccccc2)ccc1.[cH-]1cccc1.[Fe+2]'
ACETYLFERROCENE = ('CC(=O)c%30~%11c~%12c~%13c~%14c~%15%30.'
                   'c%31~%21c~%22c~%23c~%24c~%25%31.[Fe]%11%12%13%14%15%21%22%23%24%25')
FERROCENE = 'c%30~%11c~%12c~%13c~%14c~%15%30.c%31~%21c~%22c~%23c~%24c~%25%31.[Fe]%11%12%13%14%15%21%22%23%24%25'
NOT_PLANAR = {
    '[2.2]paracyclophane (CID 74210)': 'C1Cc2ccc(cc2)CCc2ccc1cc2',
    'ferrocene (CID 7611)': FERROCENE,
    '(E)-cyclodecene (CID 5364362)': 'C1CCCC/C=C/CCC1',
}


def _planar(mol):
    return [u for u in mol.stereo_units() if u['kind'] == SU_PLANAR]


def _sub(a, b):
    return [a[i] - b[i] for i in range(3)]


def _cross(a, b):
    return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]


def _unit(v):
    norm = math.sqrt(sum(x * x for x in v))
    return [x / norm for x in v]


def _ppfa_xyz():
    """PPFA by hand: Fe at the origin, eclipsed-staggered η5 rings of radius 1.2 Å at z = ±1.65 Å."""
    xyz = {29: (0.0, 0.0, 0.0)}
    for k, n in enumerate((6, 7, 21, 22, 23)):
        a = math.radians(72 * k)
        xyz[n] = (1.2 * math.cos(a), 1.2 * math.sin(a), 1.65)
    for k, n in enumerate(range(24, 29)):
        a = math.radians(72 * k + 36)
        xyz[n] = (1.2 * math.cos(a), 1.2 * math.sin(a), -1.65)
    xyz[4] = (2.7, 0.0, 1.8)
    xyz[2] = (3.2, 1.2, 2.3)
    xyz[1] = (4.6, 1.2, 2.4)
    xyz[3] = (2.6, 2.3, 1.6)
    xyz[5] = (3.2, -0.9, 2.8)
    p = math.radians(72)
    xyz[8] = (2.9 * math.cos(p), 2.9 * math.sin(p), 1.9)
    for first, centre in ((9, (0.2, 5.0, 1.2)), (15, (1.8, 3.4, 4.3))):
        for k in range(6):
            a = math.radians(60 * k)
            xyz[first + k] = (centre[0] + 1.4 * math.cos(a), centre[1], centre[2] + 1.4 * math.sin(a))
    return [xyz[n] for n in range(1, 30)]


PPFA_XYZ = _ppfa_xyz()
SOLIDS = {
    '4-carboxy[2.2]paracyclophane': (PCP_COOH, PCP_COOH_XYZ, 1),
    'pseudo-ortho dibromo[2.2]paracyclophane': (DIBROMO, PSEUDO_ORTHO_XYZ, 1),
    '(E)-cyclooctene': (TCO, TCO_XYZ, 2),
    'PPFA': (PPFA, PPFA_XYZ, 2),
}


def _from_conformer(smiles, coordinates, mirror=False):
    mol = read_smiles(smiles)
    z = {}
    with mol.edit():
        for n, (x, y, h) in zip(mol.atom_numbers, coordinates):
            mol.set_xy(n, x, y)
            z[n] = -h if mirror else h
    log = []
    assign_parities(mol, z=z, log=log)
    assert log == [], log
    return mol


def _configured(smiles, parities):
    mol = read_smiles(smiles)
    for unit, parity in zip(_planar(mol), parities):
        mol.set_parity(unit['anchor'], parity)
    return mol


def _frame(mol, unit, xyz):
    """A plane normal for the drawing: the screened ring's, or the twisted C=C's."""
    anchor, refs = unit['anchor'], unit['refs']
    if refs[2] not in mol.neighbors_of(anchor):
        (near,) = [m for m in mol.neighbors_of(anchor) if refs[2] in mol.neighbors_of(m)]
        t, c = xyz[anchor], xyz[near]
        return _cross(_sub(c, t), _sub(_sub(xyz[refs[0]], t), _sub(xyz[refs[2]], c)))
    return _cross(_sub(xyz[refs[2]], xyz[anchor]), _sub(xyz[refs[3]], xyz[anchor]))  # both ring neighbours


@mark.parametrize('name', sorted(SOLIDS))
def test_a_conformer_and_its_mirror_image_read_opposite_and_compare_unequal(name):
    smiles, coordinates, parity = SOLIDS[name]
    mol, mirror = _from_conformer(smiles, coordinates), _from_conformer(smiles, coordinates, True)
    units = _planar(mol)
    assert units and all(u['stereogenic'] for u in units)
    assert mol.parity_of(units[0]['anchor']) == parity
    assert [mirror.parity_of(u['anchor']) for u in units] == [3 - mol.parity_of(u['anchor']) for u in units]
    assert mol != mirror


@mark.parametrize('name', sorted(SOLIDS))
def test_a_wedge_on_the_face_on_drawing_states_what_the_conformer_does(name):
    smiles, coordinates, _ = SOLIDS[name]
    solid = _from_conformer(smiles, coordinates)
    for unit in _planar(solid):
        xyz = dict(zip(solid.atom_numbers, coordinates))
        w = _unit(_frame(solid, unit, xyz))
        u = _unit(_cross(w, [1.0, 0.0, 0.0] if abs(w[0]) < 0.9 else [0.0, 1.0, 0.0]))
        v = _cross(w, u)
        flat = {n: tuple(sum(p[i] * e[i] for i in range(3)) for e in (u, v, w)) for n, p in xyz.items()}
        anchor, refs = unit['anchor'], unit['refs']
        if refs[2] not in solid.neighbors_of(anchor):
            narrow, wide = anchor, refs[0]
        elif refs[1] in solid.neighbors_of(anchor):
            narrow, wide = anchor, refs[1]
        else:
            narrow, wide = refs[2], refs[1]
        mol = read_smiles(smiles)
        with mol.edit():
            for n, (x, y, _) in flat.items():
                mol.set_xy(n, x, y)
        up = flat[wide][2] > flat[narrow][2]
        for code, parity in ((WEDGE_UP, solid.parity_of(anchor)), (WEDGE_DOWN, 3 - solid.parity_of(anchor))):
            drawn = mol.copy()
            drawn.set_wedge(narrow, wide, code if up else 3 - code)
            assign_parities(drawn)
            assert drawn.parity_of(anchor) == parity, (name, anchor, code)


@mark.parametrize('name', sorted(SOLIDS))
@mark.parametrize('order', ('', 'i', 'r'))
def test_a_planar_unit_survives_a_cx_round_trip(name, order):
    smiles, coordinates, _ = SOLIDS[name]
    for mirror in (False, True):
        mol = _from_conformer(smiles, coordinates, mirror)
        log = []
        text = write_smiles(mol, order, log=log)
        assert ' |' in text and log == [], (text, log)
        back = read_smiles(text)
        assert back == mol, text
        assert back != _from_conformer(smiles, coordinates, not mirror), text


@mark.parametrize('name', sorted(SOLIDS))
def test_a_planar_unit_survives_pach(name):
    smiles, coordinates, _ = SOLIDS[name]
    mol = _from_conformer(smiles, coordinates)
    back = type(mol).unpack(mol.pack())
    assert back == mol
    assert [back.parity_of(u['anchor']) for u in _planar(back)] == [mol.parity_of(u['anchor']) for u in _planar(mol)]


def test_a_monosubstituted_paracyclophane_has_exactly_one_planar_unit():
    (unit,) = _planar(read_smiles(PCP_COOH))
    assert unit['stereogenic'] and unit['anchor'] == 4


@mark.parametrize('name', sorted(NOT_PLANAR))
def test_an_unsubstituted_deck_ferrocene_and_a_ten_ring_have_no_planar_unit(name):
    assert _planar(read_smiles(NOT_PLANAR[name])) == []


def test_a_monosubstituted_ferrocene_has_a_mirror_plane():
    (unit,) = _planar(read_smiles(ACETYLFERROCENE))
    assert not unit['stereogenic']


@mark.parametrize('smiles', ('C1CCC/C=C\\CC1', 'C1CCCC=CCC1'))
def test_a_cis_or_unstated_cyclooctene_leaves_its_planar_unit_dead(smiles):
    mol = read_smiles(smiles)
    (unit,) = _planar(mol)
    assert mol.parity_of(unit['anchor']) == 0


def test_e_cyclooctene_is_a_cis_trans_unit_and_a_planar_unit():
    mol = _from_conformer(TCO, TCO_XYZ)
    kinds = sorted(u['kind'] for u in mol.stereo_units() if u['stereogenic'])
    assert kinds == [SU_CIS_TRANS, SU_PLANAR]
    assert write_smiles(mol, '!x') == write_smiles(read_smiles(TCO), '!x')
    assert mol != read_smiles(TCO)


def test_pseudo_para_is_its_own_mirror_image_and_pseudo_ortho_is_not():
    para, ortho = (_from_conformer(DIBROMO, c) for c in (PSEUDO_PARA_XYZ, PSEUDO_ORTHO_XYZ))
    assert para == _from_conformer(DIBROMO, PSEUDO_PARA_XYZ, True)
    assert ortho != _from_conformer(DIBROMO, PSEUDO_ORTHO_XYZ, True)
    assert para != ortho


def test_of_four_dibromo_labellings_three_are_distinct():
    """(1, 2) and (2, 1) are the one achiral pseudo-para; (1, 1) and (2, 2) the pseudo-ortho pair."""
    configured = {p: _configured(DIBROMO, p) for p in ((1, 1), (1, 2), (2, 1), (2, 2))}
    assert configured[(1, 2)] == configured[(2, 1)]
    assert configured[(1, 1)] != configured[(2, 2)]
    assert configured[(1, 1)] != configured[(1, 2)]


def test_the_three_ferrocene_spellings():
    """Dative and ten-bond PPFA perceive the same units; the ionic spelling has no metallocene to perceive."""
    shapes = []
    for smiles in (PPFA, PPFA_BONDED):
        shapes.append([(u['anchor'], u['refs'], u['stereogenic']) for u in _planar(read_smiles(smiles))])
    assert shapes[0] == shapes[1] and len(shapes[0]) == 2
    assert _planar(read_smiles(PPFA_IONIC)) == []


def test_separate_ions_near_a_ring_face_log_that_the_hapto_bonds_are_missing():
    ferrocene = read_smiles('[cH-]1cccc1.[cH-]1cccc1.[Fe+2]')
    xyz = {11: (0.0, 0.0, 0.0)}
    for k in range(5):
        a = math.radians(72 * k)
        xyz[k + 1] = (1.2 * math.cos(a), 1.2 * math.sin(a), 1.65)
        xyz[k + 6] = (1.2 * math.cos(a + 0.6), 1.2 * math.sin(a + 0.6), -1.65)
    with ferrocene.edit():
        for n, (x, y, _) in xyz.items():
            ferrocene.set_xy(n, x, y)
    log = []
    assign_parities(ferrocene, z={n: p[2] for n, p in xyz.items()}, log=log)
    assert 'wedge:hapto-not-bonded' in {r.rule for r in log}, log


def test_inchi_omits_the_planar_unit_and_logs_it_once():
    mol = _from_conformer(PCP_COOH, PCP_COOH_XYZ)
    first = mol.inchi
    assert first == read_smiles(PCP_COOH).inchi
    assert mol.inchi == first
    rules = [r.rule for r in mol.log]
    assert rules.count('inchi:stereo-kind-not-represented') == 1, rules
