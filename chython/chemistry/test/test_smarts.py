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
"""What `compile_smarts` must get right, stated as behaviour rather than as a token stream.

`compile_smarts` calls the core lexer and then numbers the atoms chython 2's way, so most of this
file pins that the TABLES' dialect is still the dialect the tables were written in.  The two
deliberate divergences from chython 2, `z` and `[M]`'s charge, are asserted by name.
"""
import pytest
from chython.chemistry._smarts import SmartsSyntaxError, compile_smarts
from chython.core import read_smiles


def _atom_sets(text, smiles):
    """Every embedding of `text` into `smiles`, as a set of frozensets of molecule atom numbers."""
    query, _, _ = compile_smarts(text)
    molecule = read_smiles(smiles)
    return {frozenset(m.values()) for m in query.get_mapping(molecule)}


# --- the numbering contract -------------------------------------------------------------------- #

def test_unmapped_atoms_are_numbered_by_declaration_order():
    _, numbers, _ = compile_smarts('[O;D1;z1][N;D3;z1][C,N;z1]')
    assert sorted(numbers) == [1, 2, 3]
    assert numbers[1] < numbers[2] < numbers[3], 'stable ids rise with declaration order'


def test_explicit_map_numbers_win_over_position():
    # numbering must be B(1) H(3) B(2) H(4): the patch table's ((1,3,8),(2,4,8)) addresses the right
    # bonds only under it.
    _, numbers, _ = compile_smarts('[B;z1:1]1[H;D2:3][B;z1:2][H;D2:4]1')
    assert sorted(numbers) == [1, 2, 3, 4]
    assert numbers[1] != numbers[3]


def test_a_mixed_pattern_skips_numbers_a_map_already_claimed():
    # atom 1 is mapped :2, so the unmapped atoms take 1 and 3 -- never 2 twice.
    _, numbers, _ = compile_smarts('[C][N:2][O]')
    assert sorted(numbers) == [1, 2, 3]


def test_two_atoms_cannot_share_a_map_number():
    with pytest.raises(SmartsSyntaxError):
        compile_smarts('[C:1][N:1]')


# --- primitives -------------------------------------------------------------------------------- #

def test_an_element_matches_only_that_element():
    assert len(_atom_sets('[N]', 'CCN')) == 1
    assert not _atom_sets('[N]', 'CCO')


def test_an_element_list_matches_any_member():
    assert len(_atom_sets('[N,O]', 'NCCO')) == 2


def test_degree_counts_heavy_neighbours_only():
    assert len(_atom_sets('[C;D1]', 'CC(C)C')) == 3
    assert len(_atom_sets('[C;D3]', 'CC(C)C')) == 1


def test_comma_or_inside_one_primitive_type():
    assert len(_atom_sets('[C;D1,D3]', 'CC(C)C')) == 4


def test_implicit_hydrogen_count():
    assert len(_atom_sets('[C;h3]', 'CC(C)C')) == 3
    assert len(_atom_sets('[C;h1]', 'CC(C)C')) == 1


def test_heteroatom_count_counts_non_carbon_neighbours():
    assert len(_atom_sets('[C;x2]', 'OCN')) == 1
    assert not _atom_sets('[C;x2]', 'CCC')


def test_ring_size_membership():
    assert len(_atom_sets('[C;r6]', 'C1CCCCC1')) == 6
    assert not _atom_sets('[C;r5]', 'C1CCCCC1')


def test_not_in_a_ring_is_acyclic():
    assert len(_atom_sets('[C;!R]', 'C1CCCCC1CC')) == 2


def test_charge():
    assert len(_atom_sets('[N;+]', 'C[N+](C)(C)C')) == 1
    assert not _atom_sets('[N;+]', 'CN(C)C')


def test_a_two_digit_charge_spelling():
    assert len(_atom_sets('[B;+3]', '[B+3]')) == 1
    assert len(_atom_sets('[B;+++]', '[B+3]')) == 1, "'+++' and '+3' are the same charge"


def test_atomic_number_is_an_element():
    assert _atom_sets('[#6]', 'CCN') == _atom_sets('[C]', 'CCN')


def test_any_atom_matches_everything():
    assert len(_atom_sets('[A]', 'CCN')) == 3


def test_metal_matches_a_metal_and_nothing_else():
    assert len(_atom_sets('[M]', '[Fe]CCO')) == 1
    assert not _atom_sets('[M]', 'CCO')


def test_a_charged_metal_is_matched_by_the_spelling_the_tables_now_use():
    """`[M]` is neutral and non-radical, so the metal-organic rules spell the widening out loud.

    An unmentioned charge span is neutral for `[M]` exactly as for `[C]`, so `standardize_metals.tsv`
    writes `[M;*;^,!^]`: `*` withdraws the charge default, `^,!^` is "radical or not".
    """
    assert not _atom_sets('[M]', '[Na+]'), 'a bare [M] means what every other bracket means'
    for smiles in ('[Na+]', '[Ti+4]', '[Fe+2]', '[Fe]', '[Fe-]'):
        assert len(_atom_sets('[M;*;^,!^]', smiles)) == 1, f'the table spelling must match {smiles}'
    assert len(_atom_sets('[M;*;^,!^;D1]', '[Na+]C')) == 1, 'and D still constrains it'
    assert not _atom_sets('[M;*;^,!^;D2]', '[Na+]C')


def test_a_metal_combined_with_a_primitive_chython_2_ignored_now_means_it():
    # `[M]` is an element list like any other, so a primitive combined with it applies rather than
    # being ignored the way chython 2 ignored it.
    assert len(_atom_sets('[M;x0]', '[Fe]CC')) == 1
    assert not _atom_sets('[M;x0]', '[Fe]OC'), 'x0 constrains the metal like it constrains a carbon'


def test_an_ordinary_atom_does_default_to_neutral_and_non_radical():
    assert not _atom_sets('[C]', '[CH3-]')
    assert len(_atom_sets('[C;-]', '[CH3-]')) == 1


def test_aromatic_flag_is_the_aromatic_hybridization():
    assert _atom_sets('[C;a]', 'c1ccccc1C') == _atom_sets('[C;z4]', 'c1ccccc1C')
    assert len(_atom_sets('[C;a]', 'c1ccccc1C')) == 6


def test_isotope_needs_its_element_and_matches_only_it():
    assert len(_atom_sets('[13C]', '[13CH4]')) == 1
    assert not _atom_sets('[13C]', 'C')


def test_a_radical_is_read_from_the_cxsmarts_suffix():
    # the suffix is how the rule tables say 'this atom carries a radical'
    query, numbers, _ = compile_smarts('[O;D1;z1][N;D3;z1][C,N;z1] |^1:0,2|')
    assert len(numbers) == 3
    assert query.atom_count == 3


# --- bonds ------------------------------------------------------------------------------------- #

def test_an_absent_bond_is_single_only_and_does_not_match_aromatic():
    assert len(_atom_sets('[C][C]', 'CCC')) == 2
    assert not _atom_sets('[C][C]', 'c1ccccc1'), 'the single most common SMARTS mistake'


def test_a_colon_bond_matches_aromatic():
    assert len(_atom_sets('[C]:[C]', 'c1ccccc1')) == 6


def test_bond_orders():
    assert len(_atom_sets('[C]=[C]', 'C=CC')) == 1
    assert len(_atom_sets('[C]#[C]', 'C#CC')) == 1


def test_a_bond_or_list():
    assert len(_atom_sets('[C]-,=[C]', 'C=CC')) == 2


def test_ring_closure_builds_the_ring_bond():
    # one atom SET, found six times over the triangle's symmetry; propane matching nothing is the
    # actual assertion, since the ring bond is what makes cyclopropane match at all
    assert _atom_sets('[C]1[C][C]1', 'C1CC1') == {frozenset({1, 2, 3})}
    assert not _atom_sets('[C]1[C][C]1', 'CCC')


def test_a_ring_closure_can_carry_its_own_bond_order():
    # the '=' rides the closure digit, so it is the RING bond that must be double
    assert _atom_sets('[C]=1[C][C]1', 'C1=CC1') == {frozenset({1, 2, 3})}
    assert not _atom_sets('[C]=1[C][C]1', 'C1CC1'), 'no double bond to close on'


def test_branches_attach_to_the_atom_that_opened_them():
    assert len(_atom_sets('[C]([O])[N]', 'OCN')) == 1


# --- refusals: a typo in a knowledge file must fail where the typo is -------------------------- #

@pytest.mark.parametrize('bad', ['', '[C', '[C]1', '[C])', '[Xx]', '[C;Q2]',
                                 '[C;D]', '[:1]', '[C]=', '[C]%12[C]12'])
def test_a_malformed_pattern_raises_at_compile_time(bad):
    with pytest.raises(SmartsSyntaxError):
        compile_smarts(bad)


@pytest.mark.parametrize('text', ['([C]).([N])', '([C].[N])', '[R]', '[C&D1]'])
def test_the_lexer_accepts_what_the_subset_compiler_could_not(text):
    """Four spellings a table author may write.

    `([C]).([N])` demands the two atoms in DIFFERENT molecules and `([C].[N])` in ONE, `[R]` is any
    ring atom whatever its element, and `&` is Daylight's high AND, binding tighter than `,`.
    """
    query, numbers, _ = compile_smarts(text)
    assert len(numbers) == query.atom_count


def test_z3_deliberately_diverges_from_chython_2():
    """The one named behaviour change, asserted rather than left implicit.

    chython 2's `z` saturates at 3, so it calls an allene's central carbon sp.  The core reports 5 --
    two cumulated doubles, no triple -- and reserves 3 for sp alone.  A rule ported without
    translating this primitive silently stops matching, so the translation is done once per rule at
    extraction time and recorded in the knowledge file.
    """
    assert not _atom_sets('[C;z3]', 'C=C=C'), 'an allene carbon is not z3: z3 is sp and nothing else'
    assert len(_atom_sets('[C;z5]', 'C=C=C')) == 1, 'it is z5 there'
