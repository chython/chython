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
    assert len(roles()) == 58
    assert sum(len(rows) for rows in roles().values()) == 94


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
        'acyl_halide', 'carbamoyl_halide', 'acyl_ester', 'acyl_anhydride', 'carbamoyl',
        'michael_acceptor',
        'alkyl_amine', 'aryl_amine', 'amide_nitrogen', 'amidine_nitrogen', 'azole_nitrogen',
        'hydrazine_nitrogen',
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


def _cuts(text, role):
    mol = smiles(text)
    mol.canonicalize()
    out = []
    for fragment in mol.sticky_fragments(role):
        cut = smiles(fragment.canonical_smiles)
        cut.canonicalize()
        out.append(cut)
    return out


def _same(cuts, *expected):
    want = []
    for text in expected:
        mol = smiles(text)
        mol.canonicalize()
        want.append(mol)
    return len(cuts) == len(want) and all(any(c == w for c in cuts) for w in want)


def test_a_sulfonate_of_any_sulfonyl_group_is_cut():
    for text in ('Cc1ccc(cc1)S(=O)(=O)Oc1ccccc1', 'O=S(=O)(Oc1ccccc1)c1ccc(cc1)[N+](=O)[O-]',
                 'O=S(=O)(Oc1ccccc1)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F'):
        assert _same(_cuts(text, 'aryl_sulfonate'), '[R]c1ccccc1'), text
    assert _same(_cuts('CCS(=O)(=O)OCC', 'alkyl_sulfonate'), '[R]CC')


def test_an_ester_gives_the_acid_s_acyl_cut():
    assert _cuts('COC(=O)c1ccccc1', 'acyl_ester') == _cuts('OC(=O)c1ccccc1', 'aryl_acyl')
    assert _same(_cuts('CC(=O)ON1C(=O)CCC1=O', 'acyl_ester'), '[R]C(C)=O')
    assert not _cuts('O=C1CCCO1', 'acyl_ester')                  # a lactone is not an ester row


def test_a_cyclic_anhydride_keeps_its_other_acyl_as_the_acid():
    assert _same(_cuts('O=C1CCC(=O)O1', 'acyl_anhydride'), '[R]C(=O)CCC(=O)O')


def test_an_isocyanate_nitrogen_takes_a_hydrogen():
    assert _same(_cuts('O=C=Nc1ccccc1', 'carbamoyl'), '[R]C(=O)Nc1ccccc1')
    assert _same(_cuts('CNC(=O)Oc1ccccc1', 'carbamoyl'), '[R]C(=O)NC')


def test_a_michael_acceptor_is_cut_at_its_saturated_beta_carbon():
    assert _same(_cuts('C=CC(=O)N1CCCC1', 'michael_acceptor'), '[R]CCC(=O)N1CCCC1')
    assert _same(_cuts('C=CC#N', 'michael_acceptor'), '[R]CCC#N')
    assert _same(_cuts('C=CS(C)(=O)=O', 'michael_acceptor'), '[R]CCS(C)(=O)=O')


def test_an_aryl_hydrazine_is_cut_at_its_terminal_nitrogen():
    assert _same(_cuts('NNc1ccccc1', 'hydrazine_nitrogen'), '[R]NNc1ccccc1')
