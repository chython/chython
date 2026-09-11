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
"""`expand_abbreviations`: the contracted group a drawing wrote on one atom, turned into atoms.

Every expansion is checked against the whole compound built from SMILES, canonical form to canonical
form, so a wrong hydrogen count or a lost charge inside the fragment fails here and not later.
"""
import pytest

from .._abbreviations import expand_abbreviations
from .._tables import abbreviation_row, abbreviations_by_label, abbreviations_rows
from .._implicit import check_valence
from ...core import read_smiles


#: label -> the whole compound the graft onto benzene must equal.  Public compounds, one per class of
#: row: a plain alkyl, a charged group, an aromatic fragment, a ring closed inside the fragment, a
#: heteroatom attachment, a silicon and a group whose attachment carries hydrogens.
GRAFTS = {
    'Me': 'Cc1ccccc1',
    'tBu': 'CC(C)(C)c1ccccc1',
    'Cy': 'C1CCCCC1c1ccccc1',
    'Ph': 'c1ccccc1-c1ccccc1',
    'NO2': '[O-][N+](=O)c1ccccc1',
    'NH2': 'Nc1ccccc1',
    'NHMe': 'CNc1ccccc1',
    'OMe': 'COc1ccccc1',
    'CN': 'N#Cc1ccccc1',
    'CHO': 'O=Cc1ccccc1',
    'COOH': 'OC(=O)c1ccccc1',
    'Ts': 'Cc1ccc(cc1)S(=O)(=O)c1ccccc1',
    'Boc': 'CC(C)(C)OC(=O)c1ccccc1',
    'TMS': 'C[Si](C)(C)c1ccccc1',
    'SO2NH2': 'NS(=O)(=O)c1ccccc1',
    'Mor': 'C1COCCN1c1ccccc1',
}


def labelled(label, smiles='c1ccccc1*'):
    """`smiles` with `label` written on its one R atom -- what a reader leaves for this pass."""
    mol = read_smiles(smiles)
    mol.set_aliases({next(a.n for a in mol.atoms() if a.is_r): label})
    mol.log.clear()
    return mol


def canonical(smiles):
    mol = read_smiles(smiles)
    mol.canonicalize()
    return format(mol)


# --- the graft ----------------------------------------------------------------- #

@pytest.mark.parametrize('label', sorted(GRAFTS))
def test_a_labelled_atom_expands_to_the_whole_compound(label):
    mol = labelled(label)
    assert expand_abbreviations(mol)
    mol.canonicalize()
    assert format(mol) == canonical(GRAFTS[label])


@pytest.mark.parametrize('row', abbreviations_rows(), ids=lambda r: r.label)
def test_every_row_grafts_to_a_countable_valid_structure(row):
    mol = labelled(row.label)
    assert expand_abbreviations(mol)
    assert all(a.implicit_h is not None for a in mol.atoms()), 'a grafted atom has no hydrogen count'
    mol.log.clear()
    check_valence(mol)
    assert [r for r in mol.log if r.severity == 'info'] == list(mol.log)


@pytest.mark.parametrize('spelling', ['OMe', 'MeO', 'OCH3', 'CH3O', 'ome', 'OME'])
def test_a_synonym_and_a_case_fold_reach_the_same_row(spelling):
    mol = labelled(spelling)
    assert expand_abbreviations(mol)
    mol.canonicalize()
    assert format(mol) == canonical('COc1ccccc1')


def test_the_alias_is_dropped_and_the_expansion_is_logged_once():
    mol = labelled('OMe')
    assert expand_abbreviations(mol)
    assert not mol.aliases, 'the label is the structure now, so it is no longer a label'
    assert [r.rule for r in mol.log] == ['abbreviations:OMe']
    record = mol.log[0]
    assert record.severity == 'repaired'
    assert record.stage == 'abbreviations'
    assert '*OC' in record


def test_a_label_the_table_does_not_know_is_left_with_its_alias():
    mol = labelled('Q7')
    assert not expand_abbreviations(mol)
    assert mol.aliases == {next(a.n for a in mol.atoms() if a.is_r): b'Q7'}
    assert not mol.log


def test_a_molecule_with_no_alias_is_not_touched():
    mol = read_smiles('c1ccccc1C')
    mol.log.clear()
    assert not expand_abbreviations(mol)
    assert not mol.log


def test_several_labels_expand_in_one_pass():
    mol = read_smiles('*c1ccc(*)cc1')
    mol.set_aliases({a.n: text for a, text in zip((a for a in mol.atoms() if a.is_r), ('OMe', 'NO2'))})
    mol.log.clear()
    assert expand_abbreviations(mol)
    assert sorted(r.rule for r in mol.log) == ['abbreviations:NO2', 'abbreviations:OMe']
    mol.canonicalize()
    assert format(mol) == canonical('COc1ccc(cc1)[N+]([O-])=O')


def test_the_labelled_atom_keeps_its_id_and_a_neighbouring_parity_survives():
    # The reason the pass transmutes instead of deleting: atom 2's parity is stated over a frame of
    # neighbour ids that includes the labelled atom, and a delete-and-add would give the replacement a
    # new id at the end of that frame -- a different configuration, silently.  The group takes the
    # place the label held, so `*[C@@H]` becomes `MeO[C@@H]` and not its mirror image.
    mol = read_smiles('*[C@@H](N)CBr')
    marker = next(a.n for a in mol.atoms() if a.is_r)
    mol.log.clear()
    before = mol.parity_of(2)
    assert expand_abbreviations(mol) is False, 'no alias yet, so there is nothing to expand'
    mol.set_aliases({marker: 'OMe'})
    assert expand_abbreviations(mol)
    assert mol.parity_of(2) == before
    assert mol.atom(marker).element == 8, 'the marker atom became the fragment attachment in place'
    mol.canonicalize()
    assert format(mol) == canonical('CO[C@@H](N)CBr')


def test_the_grafted_atoms_take_the_labelled_atom_coordinates():
    mol = labelled('OMe')
    marker = next(a.n for a in mol.atoms() if a.is_r)
    with mol.edit():
        mol.set_xy(marker, 1.25, -3.5)
    grown = mol.atoms_count
    assert expand_abbreviations(mol)
    assert mol.atoms_count == grown + 1
    assert mol.xy_of(max(a.n for a in mol.atoms())) == (1.25, -3.5)
    assert '2D clean' in mol.log[0]


# --- what is refused ----------------------------------------------------------- #

def test_a_label_on_an_atom_with_two_neighbours_is_refused():
    mol = read_smiles('CC(C)C')
    mol.set_aliases({2: 'OMe'})
    mol.log.clear()
    assert not expand_abbreviations(mol)
    assert mol.aliases == {2: b'OMe'}
    assert [r.rule for r in mol.log] == ['abbreviations:attachment']
    assert mol.log[0].severity == 'refused'


def test_a_label_reached_by_a_double_bond_is_refused():
    mol = read_smiles('CC=C')
    mol.set_aliases({3: 'OMe'})
    mol.log.clear()
    assert not expand_abbreviations(mol)
    assert [r.rule for r in mol.log] == ['abbreviations:attachment']


def test_a_charge_on_the_labelled_atom_is_a_conflict_and_is_refused():
    mol = read_smiles('c1ccccc1[*-]')
    marker = next(a.n for a in mol.atoms() if a.is_r)
    mol.set_aliases({marker: 'OMe'})
    mol.log.clear()
    assert not expand_abbreviations(mol)
    assert mol.aliases == {marker: b'OMe'}
    assert [r.rule for r in mol.log] == ['abbreviations:stated-atom']
    assert mol.log[0].severity == 'refused'


def test_one_refused_site_does_not_stop_another_from_expanding():
    mol = read_smiles('*c1ccc(*)cc1')
    first, second = (a.n for a in mol.atoms() if a.is_r)
    with mol.edit():
        mol.set_charge(first, -1)
    mol.set_aliases({first: 'OMe', second: 'NO2'})
    mol.log.clear()
    assert expand_abbreviations(mol)
    assert mol.aliases == {first: b'OMe'}
    assert sorted(r.rule for r in mol.log) == ['abbreviations:NO2', 'abbreviations:stated-atom']


# --- the table ----------------------------------------------------------------- #

def test_every_row_has_one_attachment_marked_by_one_single_bond():
    for row in abbreviations_rows():
        assert row.id == f'abbreviations:{row.label}'
        assert [a.n for a in row.fragment.atoms() if a.is_r] == [row.marker]
        assert list(row.fragment.neighbors_of(row.marker)) == [row.attachment]
        assert row.fragment.order_of(row.marker, row.attachment) == 1


def test_no_row_states_a_configuration():
    # A parity in a fragment is a statement about an order of references the graft does not reproduce,
    # so the table may not hold one until the pass can carry it across.
    for row in abbreviations_rows():
        assert '@' not in row.smiles and '/' not in row.smiles and '\\' not in row.smiles


def test_no_spelling_is_claimed_twice():
    spellings = abbreviations_by_label()
    assert len(spellings) == sum(1 + len(row.synonyms) for row in abbreviations_rows())
    assert all(abbreviation_row(spelling) is row for spelling, row in spellings.items())


def test_a_label_absent_from_the_table_stays_absent():
    # The record-dependent labels: a marker or a family, not a structure, and R atoms are their home.
    for label in ('R', 'R1', 'X', 'Pol', 'Ar', 'PEG', 'Resin', 'A', 'Q'):
        assert abbreviation_row(label) is None, f'{label} names no one structure'
