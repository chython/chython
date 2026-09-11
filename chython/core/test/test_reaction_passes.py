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
"""The seven reaction-level passes that need nothing above `core`, and the log rule they all obey.

A pass over a reaction is a loop over molecules plus an answer to "which molecule was that?".  The
second half is the whole risk: a `LogRecord`'s `atoms` are stable ids IN ONE CONTAINER, so a log that
pooled three sides' records without saying which molecule each came from would hand back numbers that
name a different atom depending on which molecule you happened to read them against.  `LogRecord` has
a field for exactly that -- `subject`, whose docstring in `chython/core/_log.py` is "which molecule of
a reaction a record is about" -- and these passes fill it with the molecule's location.
`test_kekule_reports_an_unresolved_system_against_its_molecule` is the ratchet, and it reads the ids
back off the container `subject` names, which is the only reading under which they mean anything.

The sibling ratchet is `test_a_component_keeps_its_own_records_and_the_reaction_gets_a_copy`.  No pass
takes a `log=`: `rxn.log` is the destination and the component's own `.log` holds the same records
unstamped, so the two ends answer the same question.  The provenance must not leak into `rule`, because
one fact with two spellings is the defect the substrate was built to remove.

THE OTHER FOUR PASSES ARE TESTED IN `chython/chemistry/test/test_reaction_passes.py` AND NOT HERE, and
the split is not tidiness.  `standardize`, `canonicalize` and `neutralize` are registered onto
`MoleculeContainer` by injection when `chython.chemistry` is imported, and
`test_no_chython_two_imports.py` forbids a file under `chython/core/` from importing it -- because
`import chython.chemistry` runs `chython/__init__.py`, and a `core` suite that needs the façade is not
provable in isolation.  A reaction-level test of a
chemistry-level pass therefore belongs one layer up, next to the pass it exercises.  What stays here
is everything a bare `core` can answer: `kekule`, `thiele`, `reset_mapping`, `contract_ions`,
`remove_reagents`, `clean_isotopes`, `clean_stereo`, and the half of `explicify_hydrogens` that belongs
to the reaction rather than to any molecule -- giving a hydrogen added on the left and the matching one
added on the right the same map number.  That half is tested directly, through
`reaction_number_new_hydrogens`, which is why it is a function of its own: it is testable without the
molecule pass that feeds it.

The last two report rather than log, and their tests are the only ones here that assert `rxn.log` stays
EMPTY -- `_mirror` copies what a molecule wrote, and these two molecule methods write nothing.
"""
from pytest import raises

from chython.core import LOST, Log, read_reaction_smiles, read_smiles
from chython.core._reaction_passes import reaction_number_new_hydrogens


def _add_explicit_h(molecule, heavy):
    """One explicit hydrogen on `heavy`, unmapped, as `explicify_hydrogens` would leave it."""
    with molecule.edit():
        h = molecule.add_atom('H')
        molecule.add_bond(heavy, h, 1)
    return h


# --------------------------------------------------------------------------------------------------
# kekule and thiele

def test_kekule_and_thiele_are_inverse_over_the_whole_reaction():
    r = read_reaction_smiles('c1ccccc1>>c1ccncc1')
    assert r.kekule() is True
    assert r.reactants[0].smiles == 'C1=CC=CC=C1'
    assert r.kekule() is False, 'a second call changes nothing'
    assert r.thiele() is True
    assert r.reactants[0].smiles == 'c1ccccc1'


def test_kekule_reports_an_unresolved_system_against_its_molecule():
    """`C[n+]1cccc1` has no Kekule form, and the report has to say which molecule that was.

    `subject` is the field that answers it -- not a prefix on `rule`.  A record's `atoms` are stable ids
    in ONE container, so without this the (2, 3, 4, 5, 6) below names nothing in particular.

    The fixture was `Cn1cc[nH]c1` until the kekuliser learned to drop a surplus hydrogen; that is now
    repaired to 1-methylimidazole rather than refused.  A cationic three-coordinate nitrogen is
    must-match, so nothing can relax it and five must-match atoms stay an odd count.
    """
    r = read_reaction_smiles('CC>>C[n+]1cccc1')
    r.kekule()
    log = r.log
    assert log, 'an aromatic system with no Kekule form is an event'
    assert {x.subject for x in log} == {'products[0]'}, 'and the reactant said nothing'
    assert log.by_subject('products[0]') == list(log)
    assert {x.stage for x in log} == {'kekule'}

    unresolved = log.lost()
    assert len(unresolved) == 1
    assert unresolved[0].atoms == (2, 3, 4, 5, 6)
    assert unresolved[0].severity == LOST
    # INVERTED, was `== 'reaction:kekule'`: the kekuliser names its own rule now, and `absorb` fills
    # only blank provenance.  The stage asserted above is what says which pipeline it ran in.
    assert unresolved[0].rule.startswith('kekule:')
    # the atoms are ids in the molecule `subject` names, and reading them there is the whole point
    assert all(r.products[0].element_of(n) in (6, 7) for n in unresolved[0].atoms)
    # One record carries the sentence AND the ids, so there is no prose twin of this line: a summary
    # record synthesised on top of the fold would report one event twice.
    assert len(log) == 1


def test_a_component_keeps_its_own_records_and_the_reaction_gets_a_copy():
    """Both ends answer, and only the reaction-level copy carries a `subject`.

    The component's own records name atoms in the component, which is the container they mean something
    in, so `subject` there would be noise.  On `rxn.log` it is the only thing that makes the ids
    readable at all.
    """
    r = read_reaction_smiles('CC>>C[n+]1cccc1')
    r.kekule()
    product = r.products[0]
    assert product.log, 'the component holds what happened to it'
    assert all(x.subject == '' for x in product.log), 'and needs no subject to say which molecule'
    assert [str(x) for x in r.log] == [str(x) for x in product.log], 'the same records, copied'
    assert all(x.subject == 'products[0]' for x in r.log), 'stamped only on the copy'


def test_thiele_reports_the_rewrite_and_refuses_nothing():
    """One line for the ring it aromatised, subject-stamped, and no refusal."""
    r = read_reaction_smiles('C1=CC=CC=C1>>CC')
    r.thiele()
    assert [(x.rule, x.stage, x.subject) for x in r.log] == [('thiele:aromatized', 'thiele', 'reactants[0]')]
    assert r.log.refused() == []


# --------------------------------------------------------------------------------------------------
# the hydrogen numbering the reaction layer owns, tested without the molecule pass that feeds it

def test_a_new_hydrogen_pairs_across_the_arrow_by_its_heavy_atoms_number():
    """The one piece of `explicify_hydrogens` that is the REACTION's and not a molecule's.

    A hydrogen added to `[CH3:1]` on the left and one added to `[CH3:1]` on the right are the same
    hydrogen, and they have to be given the same map number or the reaction comes out with a mapping
    that says a C-H bond was broken and an identical one formed.
    """
    r = read_reaction_smiles('[CH3:1][OH:2]>>[CH3:1][NH2:3]')
    left = _add_explicit_h(r.reactants[0], 1)
    right = _add_explicit_h(r.products[0], 1)

    assert reaction_number_new_hydrogens(r, [{left}, {right}]) == 2
    log = r.log
    assert {x.subject for x in log} == {'reactants[0]', 'products[0]'}
    assert {x.stage for x in log} == {'number_new_hydrogens'}
    assert any('from the matching new hydrogen' in x for x in log), 'the pairing is an event'
    number = r.reactants[0].map_number_of(left)
    assert number > 3, 'a fresh number, above everything the record already used'
    assert r.products[0].map_number_of(right) == number


def test_an_unpaired_hydrogen_gets_a_number_of_its_own():
    r = read_reaction_smiles('[CH3:1][OH:2]>>[CH3:1][NH2:3]')
    left = _add_explicit_h(r.reactants[0], 2)      # on the oxygen, which the product does not have
    right = _add_explicit_h(r.products[0], 2)      # on the nitrogen, which the reactant does not

    assert reaction_number_new_hydrogens(r, [{left}, {right}]) == 2
    assert r.reactants[0].map_number_of(left) != r.products[0].map_number_of(right)


def test_an_unmapped_molecules_hydrogens_are_left_unmapped():
    """Numbering them would invent a mapping, and a partially mapped molecule cannot even be packed."""
    r = read_reaction_smiles('CC>>CO')
    left = _add_explicit_h(r.reactants[0], 1)
    assert reaction_number_new_hydrogens(r, [{left}, set()]) == 0
    assert r.reactants[0].map_number_of(left) == 0


# --------------------------------------------------------------------------------------------------
# reset_mapping

def test_reset_mapping_makes_every_number_unique_across_the_reaction():
    r = read_reaction_smiles('[CH3:1][CH3:2]>>[CH3:1][OH:2]')
    assert r.reset_mapping() is True
    numbers = [a.map_number for m in r.molecules() for a in m.atoms()]
    assert sorted(numbers) == [1, 2, 3, 4]


def test_reset_mapping_leaves_an_already_unique_numbering_alone():
    r = read_reaction_smiles('[CH3:1][CH3:2]>>[CH3:3][OH:4]')
    assert r.reset_mapping() is False
    assert [a.map_number for m in r.molecules() for a in m.atoms()] == [1, 2, 3, 4]


def test_reset_mapping_numbers_a_reaction_that_had_no_mapping_at_all():
    r = read_reaction_smiles('CC>>CO')
    assert r.reset_mapping() is True
    assert sorted(a.map_number for m in r.molecules() for a in m.atoms()) == [1, 2, 3, 4]


# --------------------------------------------------------------------------------------------------
# contract_ions

def test_contract_ions_merges_a_cation_and_an_anion_into_one_molecule():
    r = read_reaction_smiles('[Na+].[OH-].CC>>CO')
    assert r.contract_ions() is True
    assert len(r.reactants) == 2
    salt = [m for m in r.reactants if m.connected_components_count == 2][0]
    assert int(salt) == 0
    assert salt.smiles in ('[OH-].[Na+]', '[Na+].[OH-]')


def test_contract_ions_refuses_an_ambiguous_side():
    """Two different cations and one anion: nothing says which pairs with which."""
    r = read_reaction_smiles('[Na+].[K+].[OH-]>>CO')
    assert r.contract_ions() is False
    assert len(r.reactants) == 3


def test_contract_ions_leaves_a_side_with_no_ions_alone():
    assert read_reaction_smiles('CC>>CO').contract_ions() is False


# --------------------------------------------------------------------------------------------------
# remove_reagents

MAPPED = '[Na+:1].[OH-:2].[CH3:7][O:5][C:4]([CH3:3])=[O:6]>>[CH3:3][C:4]([OH:8])=[O:6]'


def test_remove_reagents_moves_a_molecule_with_no_reaction_centre_to_the_agents():
    r = read_reaction_smiles(MAPPED)
    assert r.remove_reagents(keep_reagents=True) is True
    assert [m.smiles for m in r.reactants] == ['O=C(C)OC']
    assert sorted(m.smiles for m in r.agents) == ['[Na+]', '[OH-]']
    assert len(r.products) == 1


def test_remove_reagents_drops_them_when_it_is_not_asked_to_keep_them():
    r = read_reaction_smiles(MAPPED)
    assert r.remove_reagents() is True
    assert r.agents == ()
    assert len(r.reactants) == 1


def test_remove_reagents_says_false_when_every_molecule_is_in_the_reaction():
    r = read_reaction_smiles('[CH3:1][CH3:2]>>[CH3:1][OH:2]')
    assert r.remove_reagents() is False


def test_remove_reagents_refuses_an_unmapped_reaction_and_names_the_other_door():
    r = read_reaction_smiles('CCO.CC(=O)O>>CC(=O)OCC.O')
    with raises(ValueError) as e:
        r.remove_reagents()
    assert 'mapping=False' in str(e.value)


def test_the_rule_based_door_moves_a_molecule_that_appears_on_both_sides():
    r = read_reaction_smiles('CCO.CC(=O)O>>CC(=O)OCC.CCO')
    assert r.remove_reagents(mapping=False, keep_reagents=True) is True
    assert [m.smiles for m in r.reactants] == ['C(C)(=O)O']
    assert [m.smiles for m in r.products] == ['C(C)OC(C)=O']
    assert [m.smiles for m in r.agents] == ['C(C)O']


def test_a_molecule_on_both_sides_yields_one_agent_and_not_two():
    """It passed through the flask once, so the agent side says once.

    The obvious implementation demotes per occurrence -- one copy off the left, one off the right, two
    agents -- and reports two equivalents of a solvent where the record shows one.  Pooling reagents
    into a `set` goes the other way and loses a genuinely consumed second equivalent.  The count is
    `min(left, right)`, which is neither.
    """
    r = read_reaction_smiles('CCO.CCO.CC(=O)O>>CC(=O)OCC.CCO')
    assert r.remove_reagents(mapping=False, keep_reagents=True) is True
    assert [m.smiles for m in r.agents] == ['C(C)O'], 'one pass-through copy'
    # and the second equivalent stays a reactant, because the products account for only one
    assert sorted(m.smiles for m in r.reactants) == ['C(C)(=O)O', 'C(C)O']
    assert [m.smiles for m in r.products] == ['C(C)OC(C)=O']


def test_the_rule_based_door_takes_the_common_reagents_as_data():
    """The predefined-solvent list is chemistry knowledge, so it is an argument and not a constant."""
    r = read_reaction_smiles('CCO.CC(=O)O>>CC(=O)OCC')
    assert r.remove_reagents(mapping=False, keep_reagents=True, common=[read_smiles('CCO')]) is True
    assert [m.smiles for m in r.reactants] == ['C(C)(=O)O']
    assert [m.smiles for m in r.agents] == ['C(C)O']


def test_the_rule_based_door_rolls_back_rather_than_empty_a_side():
    """Everything on the left is a common reagent -- so the reaction is left as it was."""
    r = read_reaction_smiles('CCO>>CC=O')
    assert r.remove_reagents(mapping=False, common=[read_smiles('CCO')]) is False
    assert len(r.reactants) == 1


# --------------------------------------------------------------------------------------------------
# the two wipes, which report instead of logging

def test_clean_isotopes_drops_every_label_on_every_side():
    r = read_reaction_smiles('[13CH3]C>>[13CH3]O')
    assert r.clean_isotopes() is True
    assert [m.smiles for m in r.molecules()] == ['CC', 'CO']
    assert r.clean_isotopes() is False, 'a second call has nothing to drop'


def test_clean_isotopes_says_false_and_writes_nothing_when_no_side_carries_one():
    """Both molecule methods REPORT and neither logs, so an empty pass leaves `rxn.log` empty --
    `_mirror` copies what a molecule wrote, and here that is nothing."""
    r = read_reaction_smiles('CC>>CO')
    assert r.clean_isotopes() is False
    assert not r.log


def test_clean_isotopes_takes_a_stranded_parity_with_the_label():
    """`C[C@H](F)[13CH3]` differs at the two methyls only by the label, so the parity goes too.

    The molecule method borrows `validate_stereo` for this; what the reaction level owes is that the
    borrowing happens per molecule and not once over a bag of atoms from three sides.
    """
    r = read_reaction_smiles('CC>>C[C@H](F)[13CH3]')
    assert r.clean_isotopes() is True
    assert r.products[0].smiles == 'C(C)(F)C', 'no isotope and no configuration'


def test_clean_stereo_wipes_every_side_and_keys_the_report_by_molecule():
    """The molecule's report is kept, keyed by location -- reducing five readers to a bool would make
    this method strictly weaker than the loop a caller would write instead."""
    r = read_reaction_smiles('C[C@H](F)Cl.CC>>C[C@@H](F)Br')
    assert r.clean_stereo() == {'reactants[0]': {'parities': [2]}, 'products[0]': {'parities': [2]}}
    assert [m.smiles for m in r.molecules()] == ['C(C)(F)Cl', 'CC', 'C(C)(F)Br']


def test_a_molecule_with_no_stereo_is_absent_from_the_report_rather_than_empty_in_it():
    """The same rule as the molecule's own report, one level up: a key means state was wiped.

    The key is also how a caller addresses the molecule again, so the index is per side and an empty
    molecule earlier on the side does not shift it.
    """
    r = read_reaction_smiles('CC.C[C@H](F)Cl>>CO')
    report = r.clean_stereo()
    assert list(report) == ['reactants[1]']
    # the ids in a value are ids in the molecule the key names, and reading them there is the point
    assert report['reactants[1]'] == {'parities': [2]}
    assert r.reactants[1].parity_of(2) == 0


def test_clean_stereo_reports_nothing_for_a_reaction_that_has_no_stereo_at_all():
    r = read_reaction_smiles('CC>>CO')
    assert r.clean_stereo() == {}
    assert not r.log
