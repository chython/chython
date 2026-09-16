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
from chython.core._core import _cip_rank_probe, _cip_rank_word_probe, _cip_ring_systems_probe


def test_ring_systems_number_fused_rings_as_one_and_leave_chains_at_zero():
    """Decalin's two fused rings are ONE system; the methyl on it is on none."""
    mol = read_smiles('CC1CCC2CCCCC2C1')
    systems = _cip_ring_systems_probe(mol)
    ring_ids = {v for v in systems.values() if v}
    assert len(ring_ids) == 1, 'fused rings must share a system'
    assert systems[1] == 0, 'the methyl carbon is on no ring'
    assert sum(1 for v in systems.values() if v) == 10


def test_ring_systems_separate_two_rings_joined_by_a_bond():
    """Bicyclohexyl is two systems: the bond between the rings is on no ring, so it fuses nothing."""
    mol = read_smiles('C1CCCCC1C1CCCCC1')
    systems = _cip_ring_systems_probe(mol)
    assert len({v for v in systems.values() if v}) == 2
    assert all(systems.values()), 'every atom is on a ring'


def test_rank_word_orders_by_atomic_number_first():
    mol = read_smiles('CNOF')
    words = [_cip_rank_word_probe(mol, n) for n in (1, 2, 3, 4)]
    assert words == sorted(words), 'C < N < O < F by atomic number'


def test_rank_word_ranks_a_real_atom_above_its_own_duplicate():
    mol = read_smiles('CC')
    assert _cip_rank_word_probe(mol, 1) > _cip_rank_word_probe(mol, 1, duplicate=True)


def test_rank_word_ranks_a_duplicate_of_a_nearer_atom_above_one_of_a_further_atom():
    """Rule 1b, encoded relative: `back` counts spheres from the duplicate up to the atom it duplicates,
    so at one sphere a bigger `back` is a shorter root distance and ranks higher."""
    mol = read_smiles('CC')
    near = _cip_rank_word_probe(mol, 1, duplicate=True, back=4)
    far = _cip_rank_word_probe(mol, 1, duplicate=True, back=1)
    assert near > far


def test_rank_word_separates_isotopes_of_one_element():
    """Rule 2, and it must not disturb rule 1a: 13C still ranks below any N."""
    mol = read_smiles('[13CH4].[12CH4].N')
    heavy, light, nitrogen = (_cip_rank_word_probe(mol, n) for n in (1, 2, 3))
    assert heavy > light
    assert nitrogen > heavy


def test_ranking_orders_the_four_halogens_of_a_bromochlorofluoro_centre():
    """Br > Cl > F > H, by atomic number alone."""
    mol = read_smiles('FC(Cl)Br')
    assert _cip_rank_probe(mol, 2) == (4, 3, 1)


def test_ranking_puts_the_implicit_hydrogen_last():
    mol = read_smiles('C[CH](O)F')
    assert _cip_rank_probe(mol, 2) == (4, 3, 1), 'F > O > C, and the implicit H is not a named direction'


def test_ranking_reports_a_tie_for_two_constitutionally_equivalent_ligands():
    """cis-1,4-dimethylcyclohexane: the two ring branches from C2 are the same walk in two directions, so
    rules 1a/1b/2 tie and phase 1 declines.  The site anchors a unit -- it is a candidate, which is what
    makes this a refusal rather than a KeyError; a site with no unit at all is the KeyError."""
    mol = read_smiles('C[C@H]1CC[C@@H](C)CC1')
    assert _cip_rank_probe(mol, 2) is None


def test_ranking_reaches_past_the_first_sphere():
    """Both ligands start with carbon; the difference is one sphere out."""
    mol = read_smiles('CC(CO)CC')
    ranked = _cip_rank_probe(mol, 2)
    assert ranked[0] == 3, 'the CH2-OH branch outranks both alkyls'


def test_ranking_ranks_an_aldehyde_branch_above_an_alcohol_branch():
    """Rule 1a at sphere 2: the CH=O carbon's children are (O, O-duplicate, H), the CH2-OH carbon's are
    (O, H, H).  Only expansion clause 1's duplicate makes the first outrank the second."""
    mol = read_smiles('OC(C=O)CO')
    ranked = _cip_rank_probe(mol, 2)
    assert ranked[0] == 1, 'the hydroxyl oxygen is the highest by rule 1a'
    assert ranked[1] == 3, 'CH=O outranks CH2OH'


def test_ranking_refuses_an_aromatic_ligand_in_this_phase():
    """A ring bond of order 4 needs the mancude mean phase 1 does not compute, so 1-phenylethanol is a
    refusal rather than a guess."""
    mol = read_smiles('C[C@H](O)c1ccccc1')
    assert _cip_rank_probe(mol, 2) is None
