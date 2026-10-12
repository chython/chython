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
"""SMILES stereo marks that state no configuration: conflicting directions and the CXSMILES `w:` and
`ctu:` fields.  Structures: bromochlorofluoromethane, but-2-ene, 1-bromo-2-chloro-2-fluoroethene.
"""

from pytest import mark

from .. import STEREO_AND, read_smiles


def _cis_trans(mol):
    (unit,) = [u for u in mol.stereo_units() if u['kind'] == 1]
    return mol.parity_of(unit['anchor'])


def _rules(mol):
    return [r.rule for r in mol.log]


@mark.parametrize('text', ['C(/F)(/Cl)=C/Br', r'F/C(\Cl)=C/Br'])
def test_two_substituents_on_one_side_leave_the_double_bond_unset(text):
    mol = read_smiles(text)
    assert _cis_trans(mol) == 0
    assert 'smiles:stereo-same-side' in _rules(mol), mol.log


def test_opposite_substituents_still_configure():
    mol = read_smiles(r'C(/F)(\Cl)=C/Br')
    assert _cis_trans(mol)
    assert 'smiles:stereo-same-side' not in _rules(mol), mol.log


def test_a_wavy_bond_unsets_a_centre():
    assert read_smiles('F[C@H](Cl)Br').parity_of(2)
    mol = read_smiles('F[C@H](Cl)Br |w:1.0|')
    assert mol.parity_of(2) == 0
    assert 'smiles:cx-wavy-bond' in _rules(mol), mol.log


def test_a_wavy_bond_at_a_terminal_unsets_the_double_bond():
    mol = read_smiles('C/C=C/C |w:1.0|')
    assert _cis_trans(mol) == 0
    assert 'smiles:cx-cis-trans-unknown' in _rules(mol), mol.log


def test_ctu_unsets_the_double_bond_it_names():
    assert _cis_trans(read_smiles('C/C=C/C')) == 1
    mol = read_smiles('C/C=C/C |ctu:1|')
    assert _cis_trans(mol) == 0
    assert 'smiles:cx-cis-trans-unknown' in _rules(mol), mol.log


@mark.parametrize('text, rule', [('F[C@H](Cl)Br |w:9.0|', 'smiles:cx-wavy-bad-index'),
                                 ('C/C=C/C |ctu:7|', 'smiles:cx-ctu-bad-index'),
                                 ('C/C=C/C |ctu:x|', 'smiles:cx-field-not-applied'),
                                 ('F[C@H](Cl)Br |w:1|', 'smiles:cx-field-not-applied')])
def test_a_bad_field_is_logged_and_the_rest_is_read(text, rule):
    mol = read_smiles(text)
    assert rule in _rules(mol), mol.log
    units = [u for u in mol.stereo_units() if u['stereogenic']]
    assert all(mol.parity_of(u['anchor']) for u in units), 'the marks the field did not name stand'


def test_a_wedge_on_a_tetrahedral_centre_is_not_applied():
    mol = read_smiles('F[C@H](Cl)Br |wU:1.0|')
    assert mol.parity_of(2)
    assert 'smiles:cx-wedge-not-applied' in _rules(mol), mol.log


def test_group_zero_is_renumbered():
    mol = read_smiles('F[C@H](Cl)Br |&0:1|')
    assert mol.stereo_groups() == {(STEREO_AND, 1): [2]}
    assert 'smiles:stereo-group-renumbered' in _rules(mol), mol.log
