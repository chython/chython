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
"""`retro()`: every one-step precursor set the corpus turns into the target, verified forward."""
import pytest
from .._enumerate import react, retro
from ...core import MoleculeContainer, read_smiles


def _target(smiles):
    molecule = read_smiles(smiles)
    molecule.canonicalize()
    return molecule


def _answers(smiles, reaction):
    return [(out.rule_id, format(out.reaction)) for out in retro(_target(smiles), reaction=reaction)]


def test_an_amide_comes_from_every_acyl_source_the_corpus_has():
    assert _answers('CC(=O)Nc1ccccc1', 'amidation') == [
        ('reactions:2', 'C(C)(=O)O.c1c(N)cccc1>>c1ccc(NC(=O)C)cc1'),
        ('reactions:514', 'C(C)(=O)Cl.c1c(N)cccc1>>c1ccc(NC(=O)C)cc1'),
        ('reactions:520', 'C(C)(=O)Br.c1c(N)cccc1>>c1ccc(NC(=O)C)cc1'),
        ('reactions:524', 'C(C)(=O)F.c1c(N)cccc1>>c1ccc(NC(=O)C)cc1'),
        ('reactions:528', 'C(C)([O-])=O.c1c(N)cccc1>>c1ccc(NC(=O)C)cc1')]


def test_every_leaving_group_alternative_is_its_own_outcome_and_symmetry_is_one():
    """Biphenyl has one biaryl bond seen from two ends; row 8 answers once per halide, not twice."""
    row = [answer for rule_id, answer in _answers('c1ccc(-c2ccccc2)cc1', 'suzuki') if rule_id == 'reactions:8']
    assert row == ['B(O)(c1ccccc1)O.c1c(Cl)cccc1>>c1cc(-c2ccccc2)ccc1',
                   'B(O)(c1ccccc1)O.c1c(Br)cccc1>>c1cc(-c2ccccc2)ccc1',
                   'B(O)(c1ccccc1)O.c1c(I)cccc1>>c1cc(-c2ccccc2)ccc1']


def test_a_lactam_opens_to_one_precursor():
    assert _answers('O=C1CCCN1', 'amidation') == [('reactions:1', 'C(O)(CCCN)=O>>C1CCC(N1)=O'),
                                                  ('reactions:4', 'C(Cl)(CCCN)=O>>C1CCC(N1)=O')]


def test_an_inverting_row_hands_back_the_other_configuration():
    answers = _answers('CC[C@@H](C)Oc1ccccc1', 'williamson')
    assert ('reactions:17', 'C(C)[C@@H](Br)C.c1c(O)cccc1>>c1ccccc1O[C@H](C)CC') in answers


def test_an_unconfigured_centre_accepts_what_the_forward_row_makes():
    assert _answers('CC(O)CC', 'ketone_to_alcohol') == [
        ('reactions:46', 'C(C)C(C)=O>>C(C)[C@H](O)C |&1:7|')]


def test_a_configured_centre_must_come_back_configured():
    """`ketone_to_alcohol` makes a racemate, so it cannot make one drawn enantiomer."""
    assert _answers('CC[C@H](C)O', 'ketone_to_alcohol') == []
    assert _answers('CC[C@@H](O)C |&1:2|', 'ketone_to_alcohol') == [
        ('reactions:46', 'C(C)C(C)=O>>C(C)[C@H](O)C |&1:7|')]


def test_an_untouched_component_rides_with_the_first_precursor():
    assert _answers('CC(=O)NCC.[Na+].[Cl-]', 'amidation')[0] == (
        'reactions:1', 'C(C)(=O)O.[Na+].[Cl-].C(C)N>>C(=O)(NCC)C.[Cl-].[Na+] |f:0.1.2|')


def test_every_answer_is_one_react_gives_back():
    target = _target('CC(=O)Nc1ccccc1')
    for out in retro(target):
        first, *others = out.reaction.reactants
        forward = {format(o.reaction) for o in react(first, others, reaction=out.name)}
        assert format(out.reaction) in forward, out


def test_the_target_is_not_touched():
    target = _target('CC(=O)Nc1ccccc1')
    before = format(target)
    list(retro(target))
    assert format(target) == before


def test_an_unknown_reaction_raises_naming_the_corpus():
    with pytest.raises(ValueError, match='unknown reaction'):
        list(retro(_target('CCO'), reaction='no_such_reaction'))


# --- inputs no row above draws -------------------------------------------------------------------

BASE = _answers('CC(=O)Nc1ccccc1', 'amidation')


def _drawn(smiles, reaction):
    """As drawn: no `canonicalize()`, which would fold the hydrogens and aromatize before `retro()` does."""
    return [(out.rule_id, format(out.reaction)) for out in retro(read_smiles(smiles), reaction=reaction)]


def test_explicit_hydrogens_change_nothing():
    assert _drawn('[H]N(C(=O)C([H])([H])[H])c1ccccc1', 'amidation') == BASE


def test_a_kekule_drawing_changes_nothing():
    assert _drawn('CC(=O)NC1=CC=CC=C1', 'amidation') == BASE


def test_map_numbers_change_nothing():
    mapped = '[CH3:1][C:2](=[O:3])[NH:4][c:5]1[cH:6][cH:7][cH:8][cH:9][cH:10]1'
    assert _drawn(mapped, 'amidation') == BASE


def test_nothing_reachable_yields_nothing():
    assert list(retro(MoleculeContainer())) == []
    assert list(retro(_target('[Ar]'))) == []


def test_an_isotope_label_survives_wherever_it_is():
    for labelled in ('[13CH3]C(=O)Nc1ccccc1', 'C[13C](=O)Nc1ccccc1'):
        answers = list(retro(_target(labelled), reaction='amidation'))
        assert len(answers) == len(BASE)
        assert all(any(atom.isotope for molecule in out.reaction.reactants for atom in molecule.atoms())
                   for out in answers)


def test_the_container_method_is_the_function():
    target = _target('CC(=O)Nc1ccccc1')
    assert [format(o.reaction) for o in target.retro(reaction='amidation')] == [answer for _, answer in BASE]


def test_the_method_takes_no_positional_partner():
    """A partner is `react()`'s question; `retro()` has one molecule and asks what made it."""
    with pytest.raises(TypeError):
        _target('CCO').retro('amidation')


# --- final review --------------------------------------------------------------------------------

def test_a_deuterium_at_the_centre_never_gives_a_hypervalent_precursor():
    """A drawn `[2H]` stays an atom through `implicify_hydrogens()`, so a reverse step raising a bond order
    beside it would make a five-valent carbon the forward row still accepts."""
    from ...core._core import valence_report
    for drawn in ('[2H]C(C)(O)CC', '[2H]OC(C)CC', 'CC([2H])(NC)c1ccccc1'):
        for out in retro(read_smiles(drawn)):
            for molecule in out.reaction.reactants:
                molecule = molecule.copy()
                molecule.kekule()
                assert not [v for v in valence_report(molecule) if v[1] == 'violation'], (drawn, format(out.reaction))


def test_only_a_stereogenic_unit_is_free():
    """A CH2 or CH3 is a stereo candidate and never a centre; counting it free walks every automorphism."""
    from .._enumerate import _Target
    target = _Target(read_smiles('C[C@H](O)CC(CC(C)(C)C)CC(C)(C)C'))
    assert target.configured and target.free == [] and target.position is None
