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
"""The molecule-level aggregates over the atom table: `is_radical` and the brutto formula.

Both are folds -- an `any` and a `Counter` -- and both are chython 2 names restored with chython 2's
meaning, which for the formula includes an atom ORDER that is not Hill's.
"""
from chython.core import MoleculeContainer, read_smiles


# --- is_radical ----------------------------------------------------------------------------------

def test_is_radical_is_false_on_a_closed_shell_molecule():
    assert read_smiles('CCO').is_radical is False


def test_is_radical_is_true_when_one_atom_carries_the_flag():
    # the ethyl radical.  THE CXSMILES RADICAL FIELD AND NOT `[CH2]`: an absent hydrogen term inside
    # brackets means zero here, so nothing about a bracket atom implies an unpaired electron, and
    # chython 2 -- which re-derives the count from valence and calls the difference a radical -- reads
    # this string as a radical for a reason this reader deliberately does not have.
    assert read_smiles('C[CH2] |^1:1|').is_radical is True


def test_is_radical_is_true_for_a_biradical():
    # dioxygen written as a biradical rather than as O=O
    mol = read_smiles('[O][O] |^1:0,1|')
    assert [mol.radical_of(n) for n in mol.atom_numbers] == [True, True]
    assert mol.is_radical is True


def test_is_radical_sees_a_radical_in_any_component():
    # TEMPO's aminoxyl next to a spectator ion: the fold is over atoms, not over components
    assert read_smiles('CC1(C)CCCC(C)(C)N1[O].[Na+] |^1:10|').is_radical is True


def test_is_radical_follows_set_radical():
    mol = read_smiles('CCO')
    assert mol.is_radical is False
    with mol.edit():
        mol.set_radical(1, True)
    assert mol.is_radical is True
    with mol.edit():
        mol.set_radical(1, False)
    assert mol.is_radical is False


def test_is_radical_of_an_empty_molecule_is_false():
    assert MoleculeContainer().is_radical is False


# --- brutto and brutto_formula -------------------------------------------------------------------

def test_brutto_folds_implicit_hydrogens_into_h():
    # ethanol has no hydrogen ATOM; all six are implicit counts
    mol = read_smiles('CCO')
    assert mol.element_counts == {6: 2, 8: 1}
    assert mol.brutto == {'C': 2, 'H': 6, 'O': 1}


def test_brutto_is_keyed_by_symbol_and_ordered_c_h_o_n_b_then_by_atomic_number():
    # 4-nitroaniline: C and H and O and N lead, in that order, whatever their atomic numbers
    mol = read_smiles('Nc1ccc(cc1)[N+](=O)[O-]')
    assert list(mol.brutto) == ['C', 'H', 'O', 'N']
    # trifluoromethanesulfonic acid: O before F before S, the tail sorted by atomic number
    assert list(read_smiles('FC(F)(F)S(=O)(=O)O').brutto) == ['C', 'H', 'O', 'F', 'S']
    # and B is fifth in the lead, after N, not sorted among the tail
    assert list(read_smiles('B(O)(O)c1ccccc1').brutto) == ['C', 'H', 'O', 'B']


def test_brutto_omits_an_element_that_is_absent():
    # no seeded key survives with a zero: the dict holds only what the molecule has
    assert read_smiles('[Na+]').brutto == {'Na': 1}


def test_brutto_ignores_isotopes():
    # heavy water is H2O; a formula counts elements, and the isotope is on the atom
    assert read_smiles('[2H]O[2H]').brutto == {'H': 2, 'O': 1}
    assert read_smiles('[13CH4]').brutto == {'C': 1, 'H': 4}


def test_brutto_ignores_charge():
    # the ammonium ion is H4N and the acetate C2H3O2; neither formula carries the charge
    assert read_smiles('[NH4+]').brutto == {'H': 4, 'N': 1}
    assert read_smiles('CC(=O)[O-]').brutto == {'C': 2, 'H': 3, 'O': 2}


def test_brutto_counts_a_hydrogen_atom_and_an_implicit_hydrogen_alike():
    # methane written three ways: no hydrogen atoms, four of them, and a mixture
    assert read_smiles('C').brutto == read_smiles('[H]C([H])([H])[H]').brutto
    assert read_smiles('[H]C').brutto == {'C': 1, 'H': 4}


def test_brutto_covers_every_component():
    assert read_smiles('CC(=O)[O-].[Na+]').brutto == {'C': 2, 'H': 3, 'O': 2, 'Na': 1}


def test_brutto_of_an_empty_molecule_is_an_empty_dict():
    assert MoleculeContainer().brutto == {}


def test_brutto_formula_drops_a_count_of_one():
    assert read_smiles('CC(=O)Oc1ccccc1C(=O)O').brutto_formula == 'C9H8O4'   # aspirin
    assert read_smiles('c1ccc2[nH]ccc2c1').brutto_formula == 'C8H7N'          # indole
    assert read_smiles('[Na+]').brutto_formula == 'Na'


def test_brutto_formula_keeps_bruttos_order_rather_than_hills():
    # Hill would write CHF3O3S; chython 2 puts O before F because O is in the seeded lead
    assert read_smiles('FC(F)(F)S(=O)(=O)O').brutto_formula == 'CHO3F3S'
    # and a salt trails its metal, because sodium sorts after oxygen by atomic number
    assert read_smiles('CC(=O)[O-].[Na+]').brutto_formula == 'C2H3O2Na'


def test_brutto_formula_of_an_empty_molecule_is_an_empty_string():
    assert MoleculeContainer().brutto_formula == ''


def test_brutto_formula_html_subscripts_every_count_above_one():
    assert read_smiles('CC(=O)Oc1ccccc1C(=O)O').brutto_formula_html == \
        'C<sub>9</sub>H<sub>8</sub>O<sub>4</sub>'


def test_brutto_formula_html_leaves_a_count_of_one_bare_and_keeps_bruttos_order():
    assert read_smiles('FC(F)(F)S(=O)(=O)O').brutto_formula_html == 'CHO<sub>3</sub>F<sub>3</sub>S'


def test_brutto_formula_html_of_an_empty_molecule_is_an_empty_string():
    assert MoleculeContainer().brutto_formula_html == ''


def test_brutto_says_nothing_about_a_hydrogen_count_it_does_not_have():
    # an atom whose implicit count is the sentinel contributes no H, the same silence as `float(mol)`
    mol = MoleculeContainer()
    with mol.edit():
        mol.add_atom('C', implicit_h=None)
    assert mol.unknown_h_count == 1
    assert mol.brutto == {'C': 1}


# --- the chython 2 witness -----------------------------------------------------------------------
#
# Not an oracle: every answer above is stated outright.  What this catches is the one thing a stated
# answer cannot -- that the ORDER and the hydrogen folding are chython 2's on a molecule nobody
# thought to write a case for.  Reached through `oracle`, an installed chython 2 in another
# interpreter, so this file imports no chython 2.
#
# THE RADICAL IS SPELT IN CXSMILES ON BOTH SIDES.  Bare `C[CH2]` is a radical to chython 2 and is not
# one here, and that divergence belongs to the READERS -- chython 2 re-derives a bracket atom's
# hydrogen count from valence rules and calls the shortfall a radical, where an absent count inside
# brackets means zero here.  Comparing it would compare two different molecules and say nothing about
# the fold under test; `|^1:1|` states the radical outright and both readers agree on it.

COMPOUNDS = ('CCO', 'c1ccccc1', 'CC(=O)Oc1ccccc1C(=O)O', '[NH4+].[Cl-]', 'CC(=O)[O-].[Na+]',
             '[2H]O[2H]', '[13CH4]', 'C[CH2] |^1:1|', '[O-][N+](=O)c1ccccc1',
             'FC(F)(F)S(=O)(=O)O',
             'B(O)(O)c1ccccc1', '[Fe+2].[Cl-].[Cl-]', 'N', 'O', '[Na+]', 'CC[Si](C)(C)C',
             'c1ccc2[nH]ccc2c1', 'O=[U](=O)([O-])[O-]', '[H][H]',
             'CCCCCCCCCCCCCCCCCC(=O)O')

V2_VALUES = """
from chython import smiles

out = []
for smi in _payload:
    mol = smiles(smi)
    out.append([list(mol.brutto.items()), mol.brutto_formula, mol.is_radical])
_emit(out)
"""


def test_the_aggregates_agree_with_chython_two():
    from .oracle import ask

    answers = ask(V2_VALUES, list(COMPOUNDS))
    assert len(answers) == len(COMPOUNDS)
    for smi, (brutto, formula, radical) in zip(COMPOUNDS, answers):
        mol = read_smiles(smi)
        # `list(...items())` on both sides: the order is part of the answer, and comparing dicts
        # would pass on a molecule whose formula chython 2 spells in a different sequence
        assert list(mol.brutto.items()) == [tuple(x) for x in brutto], smi
        assert mol.brutto_formula == formula, smi
        assert mol.is_radical == radical, smi
