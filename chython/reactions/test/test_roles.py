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
"""Every `roles.tsv` row composes, names a real group, states its cap, and fires on its example."""
from pytest import raises
from chython.core import read_smiles as smiles
from chython.reactions._tables import ROLE_CAP, Role, functional_rules, roles


def test_every_row_names_a_known_group():
    known = functional_rules()
    for name, rows in roles().items():
        for row in rows:
            assert row.group in known, f'{row.id} ({name}) names group {row.group!r}'


def test_every_row_builds_its_cap():
    for name, rows in roles().items():
        for row in rows:
            assert f'#0:{ROLE_CAP}]' in row.product, \
                f'{row.id} ({name}): {row.product!r} states no `[#0:{ROLE_CAP}]`'


def test_every_row_fires_on_its_example():
    known = functional_rules()
    for name, rows in roles().items():
        for row in rows:
            probe = row.example or known[row.group].example
            mol = smiles(probe)
            mol.canonicalize()
            outcomes = list(row.template(mol, report=True))
            assert outcomes, f'{row.id} ({name}) does not fire on {probe!r}'
            for reaction, where in outcomes:
                assert ROLE_CAP in where, f'{row.id} ({name}) built no atom for :{ROLE_CAP}'
                marker = where[ROLE_CAP]
                product = next(p for p in reaction.products if marker in p.atom_numbers)
                assert product.atom(marker).is_r, f'{row.id} ({name}) capped with a real element'
                assert len(product.neighbors_of(marker)) == 1, \
                    f'{row.id} ({name}) built a cap with more than one bond'


def test_no_row_fires_on_a_decoy():
    # An empty `decoys` falls back to the group's own, the way an empty `example` does: a role's reactant
    # side IS its group's SMARTS, so the group's decoys are the ones that must stay unmatched.  Without
    # the fallback this test is vacuous, every row's `decoys` cell being empty.
    known = functional_rules()
    for name, rows in roles().items():
        for row in rows:
            for decoy in row.decoys or known[row.group].decoys:
                mol = smiles(decoy)
                mol.canonicalize()
                assert not list(row.template(mol)), f'{row.id} ({name}) fires on decoy {decoy!r}'


def test_rows_are_grouped_by_role():
    assert isinstance(roles()['aryl_halide'], tuple)
    assert {row.group for row in roles()['aryl_halide']} == \
           {'aryl_chloride', 'aryl_bromide', 'aryl_iodide'}


def test_ids_are_table_qualified():
    for rows in roles().values():
        for row in rows:
            assert row.id.startswith('roles:')


def test_a_product_stating_no_cap_is_refused_at_load():
    from chython.reactions._tables import _cap_in_product

    _cap_in_product(f'[A:1][#0:{ROLE_CAP}]', 'roles:probe')
    with raises(ValueError, match='states no cap'):
        _cap_in_product('[A:1]-[A:2]', 'roles:probe')
    # The cap is the marker under ONE number table-wide, so `where[ROLE_CAP]` names it without a scan:
    # an unmapped marker and a differently numbered one are both refused.
    with raises(ValueError, match='states no cap'):
        _cap_in_product('[A:1][#0]', 'roles:probe')
    with raises(ValueError, match='states no cap'):
        _cap_in_product('[A:1][#0:2]', 'roles:probe')


def test_the_glossary_is_complete():
    assert len(roles()) == 53
    assert sum(len(rows) for rows in roles().values()) == 87


def test_every_role_family_is_present():
    assert set(roles()) == {
        'aryl_halide', 'alkyl_halide', 'alkenyl_halide', 'alkynyl_halide',
        'aryl_fluoride', 'alkyl_fluoride', 'alkenyl_fluoride', 'alkynyl_fluoride',
        'aryl_sulfonate', 'alkyl_sulfonate',
        'aryl_boron', 'alkyl_boron', 'alkenyl_boron', 'alkynyl_boron',
        'aryl_magnesium', 'alkyl_magnesium', 'alkenyl_magnesium',
        'aryl_zinc', 'alkyl_zinc', 'alkenyl_zinc',
        'aryl_stannane', 'alkyl_stannane', 'alkenyl_stannane',
        'aryl_silane', 'alkenyl_silane', 'alkynyl_silane',
        'alkyl_acyl', 'aryl_acyl', 'alkenyl_acyl', 'alkynyl_acyl',
        'acyl_halide', 'carbamoyl_halide',
        'alkyl_amine', 'aryl_amine', 'amide_nitrogen', 'amidine_nitrogen', 'azole_nitrogen',
        'alkyl_thiol', 'aryl_thiol', 'alkyl_hydroxyl', 'aryl_hydroxyl', 'acid_hydroxyl',
        'alkynyl_terminal', 'sulfonyl',
        'alkyl_deoxy', 'aryl_deoxy', 'carbonyl_electrophile',
        'alkyl_decarboxy', 'aryl_decarboxy', 'alkenyl_decarboxy', 'alkynyl_decarboxy',
        'alkyl_deamino', 'aryl_deamino'}


def test_acid_hydroxyl_reaches_the_hydroxyl_oxygen():
    # `3` is written out rather than read from `row.cap`: what this pins is `functional.tsv`'s
    # `carboxylic_acid` numbering its hydroxyl O, which is the only reason this role can cap it.
    row = roles()['acid_hydroxyl'][0]
    mol = smiles('OC(=O)c1ccccc1')
    mol.canonicalize()
    reaction, where = next(row.template(mol, report=True))
    product = next(iter(reaction.products))
    assert product.atom(where[3]).atomic_symbol == 'O'
    assert list(product.neighbors_of(where[ROLE_CAP])) == [where[3]]


def test_the_deoxy_roles_drop_the_oxygen():
    row = next(r for r in roles()['alkyl_deoxy'] if r.group == 'primary_alcohol')
    mol = smiles('CCO')
    mol.canonicalize()
    reaction, where = next(row.template(mol, report=True))
    product = next(iter(reaction.products))
    assert 'O' not in product.brutto
    assert product.atom(where[2]).atomic_symbol == 'C'
