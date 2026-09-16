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
from chython.core import read_smiles


def test_assign_labels_bromochlorofluoromethane_as_written():
    """(S)-bromochlorofluoromethane.  Br > Cl > F > H by rule 1a alone, and this is the case that
    calibrates even-parity-means-R: a leading bracket atom takes its implicit hydrogen as the FIRST
    neighbour, so the SMILES order is (H, F, Cl, Br), descending priority is (Br, Cl, F, H), the
    permutation between them is even, and `@`'s odd parity therefore stays odd -- odd is S."""
    mol = read_smiles('[C@H](F)(Cl)Br')
    assert mol.assign_cip() is True
    assert mol.atom_cips() == {1: 'S'}


def test_assign_flips_the_letter_with_the_parity():
    a = read_smiles('[C@H](F)(Cl)Br')
    b = read_smiles('[C@@H](F)(Cl)Br')
    a.assign_cip()
    b.assign_cip()
    assert (a.atom_cips()[1], b.atom_cips()[1]) == ('S', 'R')


def test_assign_overrides_a_stated_descriptor_and_logs_the_disagreement():
    mol = read_smiles('[C@H](F)(Cl)Br')
    mol.set_atom_cip(1, 'R')
    assert mol.assign_cip() is True
    assert mol.atom_cips() == {1: 'S'}, 'the computed answer overrides the stated claim'
    assert any(r.rule == 'cip:disagreed' for r in mol.log.by_stage('cip'))


def test_assign_refuses_a_tied_centre_and_says_so():
    """An unresolved rule-1/2 tie is a refusal with a log line, never a guess."""
    mol = read_smiles('C[C@H]1CC[C@@H](C)CC1')
    mol.assign_cip()
    assert mol.atom_cips() == {}
    assert any(r.rule == 'cip:undecided' for r in mol.log.by_stage('cip'))


def test_assign_refuses_an_aromatic_ligand_rather_than_guessing():
    """1-phenylethanol.  A ring bond of order 4 needs the mancude mean this phase does not compute."""
    mol = read_smiles('C[C@H](O)c1ccccc1')
    assert mol.assign_cip() is False
    assert mol.atom_cips() == {}
    assert any(r.rule == 'cip:undecided' for r in mol.log.by_stage('cip'))


def test_assign_says_nothing_about_an_unconfigured_centre():
    """A ranking without a configuration is not a descriptor, and a molecule that states nothing is not
    a molecule with a problem: no descriptor and no line."""
    mol = read_smiles('CC(O)F')
    assert mol.assign_cip() is False
    assert mol.atom_cips() == {}
    assert [r for r in mol.log.by_stage('cip') if r.rule == 'cip:undecided'] == []


def test_assign_is_idempotent_and_the_second_run_disagrees_with_nothing():
    """The pass reads constitution and parity, never a stored descriptor, so the second run cannot
    disagree with the first -- and if it does, something read the nibble it must not read."""
    mol = read_smiles('[C@H](F)(Cl)Br')
    mol.assign_cip()
    first = mol.atom_cips()
    mol.assign_cip()
    assert mol.atom_cips() == first
    assert sum(1 for r in mol.log.by_stage('cip') if r.rule == 'cip:disagreed') == 0


def test_assign_labels_both_centres_of_a_two_centre_molecule():
    """L-threonine's backbone, written without the acid so rule 1a alone settles both sites."""
    mol = read_smiles('C[C@@H](O)[C@@H](N)CO')
    assert mol.assign_cip() is True
    assert set(mol.atom_cips()) == {2, 4}
