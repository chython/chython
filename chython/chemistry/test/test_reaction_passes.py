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
"""The reaction-level passes whose molecule-level pass is injected by this package.

The loops live in `chython/core/_reaction_passes.py`, but the method each iteration calls is injected
onto `MoleculeContainer` on `import chython.chemistry`, so these tests need this package.  They cannot
move under `chython/core/test/`, where importing `chython.chemistry` is forbidden.  The passes needing
nothing above `core` are tested in `chython/core/test/test_reaction_passes.py`.
"""
import chython.chemistry  # noqa: F401  -- injects standardize/canonicalize onto MoleculeContainer
from chython.core import Log, LogRecord, read_reaction_smiles


def test_standardize_runs_on_every_side():
    """The iron carbonyl is the differential suite's own witness that a patch fires."""
    r = read_reaction_smiles('C(=O)[Fe](C=O)C=O>>CC.C(=O)[Fe](C=O)C=O')
    assert r.standardize() is True
    assert r.reactants[0].smiles == '[O+]#[C-]~[Fe](~[C-]#[O+])~[C-]#[O+]'
    assert r.products[1].smiles == '[O+]#[C-]~[Fe](~[C-]#[O+])~[C-]#[O+]'


def test_standardize_says_false_when_nothing_was_mis_drawn():
    assert read_reaction_smiles('CC>>CO').standardize() is False


def test_the_log_says_which_molecule_each_record_came_from():
    """A stable id means nothing without its container, so `subject` names the container."""
    r = read_reaction_smiles('CC>C(=O)[Fe](C=O)C=O>C(=O)[Fe](C=O)C=O')
    assert r.standardize()
    log = r.log
    assert all(isinstance(x, LogRecord) for x in log)
    assert {x.subject for x in log} == {'agents[0]', 'products[0]'}, 'the ethane said nothing'
    assert {x.stage for x in log} == {'standardize'}
    # the molecule pass's own rule id is untouched -- no location was encoded into it
    assert all(x.rule.startswith('metals:') for x in log)
    # and the ids in each record are readable against the molecule its subject names
    for record in log.by_subject('products[0]'):
        assert all(n in r.products[0].atom_numbers for n in record.atoms)


def test_canonicalize_stamps_its_own_stages_under_the_molecule_it_ran_on():
    """`canonicalize` is itself a pipeline, so the two provenance fields must compose.

    The reaction layer stamps `subject` only, leaving each molecule pass's own `stage` in place, or the
    inner stage names are lost and `log.by_stage('standardize')` answers nothing.
    """
    r = read_reaction_smiles('CC>>C(=O)[Fe](C=O)C=O')
    assert r.canonicalize() is True
    log = r.log
    assert {x.subject for x in log} == {'products[0]'}
    assert 'standardize' in {x.stage for x in log}, 'the inner stage survived the outer scope'


def test_the_records_are_there_without_anyone_asking():
    """There is no `log=` to pass and no way to switch recording off; `rxn.log` is the destination."""
    r = read_reaction_smiles('C(=O)[Fe](C=O)C=O>>CC')
    assert r.standardize() is True
    assert r.log, 'a repair with nobody watching is still a repair'
    assert r.reactants[0].log, 'and the component holds its own copy'


def test_canonicalize_runs_over_every_side():
    assert read_reaction_smiles('C1=CC=CC=C1>>CC').canonicalize() is True


# neutralize

def test_neutralize_runs_on_every_side():
    # `|f:0.1|` and not a bare `.`: in a reaction SMILES a dot separates MOLECULES, and the fragment
    # grouping is what says the ammonium and the chloride are one recorded compound.
    r = read_reaction_smiles('C[NH3+].[Cl-]>>[NH3+]CC(=O)[O-] |f:0.1|')
    assert r.neutralize() is True
    assert r.reactants[0].smiles == 'CN.Cl'
    assert r.products[0].smiles == 'C(CN)(=O)O'


def test_neutralize_says_false_when_there_is_nothing_to_pair():
    assert read_reaction_smiles('CN.Cl>>CC').neutralize() is False


def test_neutralize_never_pairs_across_a_molecule_boundary():
    """Two ions written as separate components of one reactant pair; written as two reactants they do
    not, because a proton crossing that boundary would change what each recorded compound is."""
    r = read_reaction_smiles('C[NH3+].[Cl-]>>CC |f:0.1|')
    assert r.neutralize() is True

    r = read_reaction_smiles('C[NH3+].[Cl-]>>CC')
    assert r.neutralize() is False
    assert r.reactants[0].smiles == 'C[NH3+]'


def test_neutralize_forwards_keep_charge():
    r = read_reaction_smiles('C[NH3+]>>CC')
    assert r.neutralize(keep_charge=False) is True
    assert r.reactants[0].smiles == 'CN'


def test_neutralize_records_are_stamped_with_the_molecule_they_came_from():
    r = read_reaction_smiles('C[NH3+].[Cl-]>>CC |f:0.1|')
    r.neutralize()
    records = r.log.by_stage('neutralize')
    assert len(records) == 1
    assert records[0].subject == 'reactants[0]'
    assert not r.reactants[0].log.by_stage('neutralize')[0].subject, \
        'a molecule keeps its own unstamped copy; only the reaction pool needs a subject'


# the hydrogen pair, end to end -- the numbering is tested on its own in the core suite

def test_implicify_hydrogens_counts_across_the_whole_reaction():
    r = read_reaction_smiles('[H]CC>>CC[H]')
    assert r.implicify_hydrogens() == 2


def test_explicify_hydrogens_counts_across_the_whole_reaction():
    r = read_reaction_smiles('C>>C')
    assert r.explicify_hydrogens() == 8


def test_explicify_hydrogens_leaves_a_mapping_that_says_nothing_happened_to_the_hydrogens():
    """Methanol to methylamine: three C-H bonds survive, so the three new hydrogens on each side must
    carry the same three map numbers, or the mapping claims three C-H bonds broke and reformed.
    """
    r = read_reaction_smiles('[CH3:1][OH:2]>>[CH3:1][NH2:3]')
    assert r.explicify_hydrogens() == 4 + 5

    def hydrogens_on(molecule, heavy_map):
        return {molecule.map_number_of(n) for atom in molecule.atoms() if atom.element == 1
                for n in [atom.n]
                if molecule.map_number_of(next(iter(molecule.neighbors_of(n)))) == heavy_map}

    left = hydrogens_on(r.reactants[0], 1)
    right = hydrogens_on(r.products[0], 1)
    assert len(left) == 3 and left == right, 'the methyl hydrogens are the same three hydrogens'
    assert 0 not in left, 'a mapped molecule must not come out partially mapped'
    # the hydroxyl and amine hydrogens are NOT the same hydrogen, and do not share a number
    assert not hydrogens_on(r.reactants[0], 2) & hydrogens_on(r.products[0], 3)
    # and no molecule reuses a number inside itself -- a number repeated across the arrow is the
    # mapping doing its job, a number repeated within one molecule is a collision
    for molecule in r.molecules():
        numbers = [a.map_number for a in molecule.atoms()]
        assert len(set(numbers)) == len(numbers), [m.smiles for m in r.molecules()]
