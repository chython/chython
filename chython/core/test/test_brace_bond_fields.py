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
"""The `{...}` tail's bond fields: `C<i>,<j>=@`/`@@`, an atropisomer axis's configuration, and
`B<i>,<j>=<letter>`, a bond's stated descriptor.

| field         | frame / mapping                                                                  |
| ------------- | -------------------------------------------------------------------------------- |
| `C<i>,<j>=@@` | parity 1, each pivot's ortho atoms in string order, the lower-indexed pivot first |
| `C<i>,<j>=@`  | parity 2 in the same frame                                                       |
| `B` `e` / `z` | stored E / Z                                                                     |
| `B` `r` / `s` | an axis's Ra / Sa, stored as M / P                                               |
| `B` `q`       | no determination: INFO, nothing stored                                           |
| `B` other     | dropped, LOST                                                                    |

Compounds: BINOL (CID 11762), 6,6'-dinitrodiphenic acid (CID 95234), (E)-but-2-ene.
"""

from pytest import mark

from .. import read_smiles, write_smiles


SU_ATROPISOMER = 3
COMPOUNDS = {
    'BINOL': 'Oc1ccc2ccccc2c1-c1c(O)ccc2ccccc12',
    'dinitrodiphenic acid': 'OC(=O)c1cccc(c1-c1c(cccc1[N+](=O)[O-])C(=O)O)[N+](=O)[O-]',
}


def _axis(mol):
    (unit,) = [u for u in mol.stereo_units() if u['kind'] == SU_ATROPISOMER and u['stereogenic']]
    return unit


def _rules(mol):
    return [r.rule for r in mol.log]


@mark.parametrize('name', sorted(COMPOUNDS))
def test_an_axis_reads_its_configuration_from_the_c_field(name):
    text = COMPOUNDS[name]
    unit = _axis(read_smiles(text))
    i, j = unit['owners'][0] - 1, unit['owners'][1] - 1
    read = {sign: read_smiles(f'{text} {{C{i},{j}={sign}}}') for sign in ('@', '@@')}
    assert read['@'] != read['@@']
    for sign, parity in (('@@', 1), ('@', 2)):
        mol = read[sign]
        assert _rules(mol) == [], _rules(mol)
        assert mol.parity_of(unit['anchor']) == parity
        configured = read_smiles(text)
        configured.set_parity(unit['anchor'], parity)
        assert mol == configured
        assert read_smiles(write_smiles(mol)) == mol       # the same configuration through chython's CX wedges
    assert read_smiles(f'{text} {{C{j},{i}=@@}}') == read['@@'], 'the pivot named first does not matter'


def test_the_registry_spelling_of_an_axis_keeps_its_group_and_descriptor():
    text = COMPOUNDS['BINOL']
    unit = _axis(read_smiles(text))
    a, b = unit['owners']
    mol = read_smiles(f'{text} {{B{a - 1},{b - 1}=r;C{a - 1},{b - 1}=@@;o1:{a - 1}}}')
    assert _rules(mol) == []
    assert mol.stereo_groups() == {(2, 1): [(a, b)]}
    assert mol.parity_of(unit['anchor']) == 1 and mol.stereo_is_live((a, b))
    assert mol.bond_cip_of(a, b) == 'M'


def test_a_c_field_on_a_bond_that_is_no_axis_stores_nothing():
    mol = read_smiles(COMPOUNDS['BINOL'] + ' {C0,1=@@}')
    assert _rules(mol) == ['smiles:brace-axis-not-a-unit']
    assert mol == read_smiles(COMPOUNDS['BINOL'])


@mark.parametrize('letter, stored, rule', [
    ('e', 'E', None), ('z', 'Z', None),
    ('q', None, 'smiles:cip-undetermined'),
    ('x', None, 'smiles:cip-field-bad-descriptor'),
])
def test_a_b_field_states_a_bond_descriptor(letter, stored, rule):
    mol = read_smiles(f'C/C=C/C {{B1,2={letter}}}')
    assert mol.bond_cip_of(2, 3) == stored
    assert _rules(mol) == ([rule] if rule else [])


def test_a_b_field_on_atoms_that_share_no_bond_is_dropped():
    mol = read_smiles('C/C=C/C {B0,3=e}')
    assert _rules(mol) == ['smiles:brace-bond-not-a-bond']
