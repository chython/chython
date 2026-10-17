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
"""Chains: the segmentation every chain layout reads, a peptide the first of them."""
import pytest

from chython.core import smiles
from chython.core.chains import Chain, Circle, cycles, perceive
from chython.core.monomers import Peptide, monomers
from .peptides import PEPTIDES, molecule


def test_peptide_is_a_chain():
    assert Peptide is Chain


@pytest.mark.parametrize('name', PEPTIDES)
def test_perceive_takes_a_peptide_as_monomers_does(name):
    mol = molecule(name)
    chain = perceive(mol)
    assert chain == monomers(mol) and chain.kind == 'peptide' and chain.beads == ()


ALPHA_CD = ('OC[C@H]1O[C@@H]2O[C@H]3[C@H](O)[C@@H](O)[C@H](O[C@@H]3CO)O[C@H]3[C@H](O)[C@@H](O)[C@H](O[C@@H]3CO)'
            'O[C@H]3[C@H](O)[C@@H](O)[C@H](O[C@@H]3CO)O[C@H]3[C@H](O)[C@@H](O)[C@H](O[C@@H]3CO)O[C@H]3[C@H](O)'
            '[C@@H](O)[C@H](O[C@@H]3CO)O[C@H]1[C@H](O)[C@H]2O')                                     # alpha-cyclodextrin
CALIX4 = 'Oc1c2Cc3cccc(Cc4cccc(Cc5cccc(Cc1ccc2)c5O)c4O)c3O'                                     # calix[4]arene


def beads(chain):
    return sum(r.kind == 'bead' for r in chain.residues)


def test_a_cyclodextrin_is_a_closed_chain_of_glucose_beads():
    chain = perceive(smiles(ALPHA_CD))
    assert chain.kind == 'cycle' and chain.head_to_tail and beads(chain) == 6 and len(chain.beads) == 6


def test_a_calixarene_is_a_chain_of_four_aryl_arcs():
    chain = perceive(smiles(CALIX4))
    assert [r.kind for r in chain.residues].count('arc') == 4 and chain.beads == () and len(chain.backbone) == 16


@pytest.mark.parametrize('text', ['C1COCCOCCOCCOCCOCCO1', 'C1COCCOc2ccccc2OCCOCCOc2ccccc2O1', 'C1CNCCNCCNCCNCCNCCN1'],
                         ids=['18-crown-6', 'dibenzo-18-crown-6', 'hexaaza-18-crown-6'])
def test_a_crown_is_a_circle(text):
    out = cycles(smiles(text))
    assert isinstance(out, Circle) and out.period in (3, 4)


@pytest.mark.parametrize('text', ['C1CCCCCCCCCCCCCCC1', 'CCC1OC(=O)C(C)C(O)C(C)C(O)C(C)(O)CC(C)C(=O)C(C)C(O)C1(C)O'],
                         ids=['cyclohexadecane', 'erythronolide'])
def test_any_other_plain_macrocycle_is_a_chain_of_atoms(text):
    out = cycles(smiles(text))
    assert isinstance(out, Chain) and out.beads == () and out.head_to_tail


@pytest.mark.parametrize('text', ['C1=CC2=CC3=CC=C(N3)C=C4C=CC(=N4)C=C5C=CC(=N5)C=C1N2', 'C1CCCCC1', 'C1CCCCCCCCCCC1',
                                  'C1CC2=NC1=CC3=NC(=CC4=NC(C5=NC(=C2)CC5)CC4)CC3'],
                         ids=['porphine', 'cyclohexane', 'cyclododecane', 'corrin'])
def test_no_cycle(text):
    assert cycles(smiles(text)) is None


def test_a_cyclic_peptide_stays_a_peptide():
    assert perceive(molecule('cyclosporin A')).kind == 'peptide'


def test_a_benzyl_on_a_glucose_is_a_group_not_a_bead():
    mol = smiles(ALPHA_CD.replace('OC[C@H]1', 'c9ccccc9COC[C@H]1', 1))
    chain = perceive(mol)
    assert beads(chain) == 6 and not any(mol.hybridization_of(x) == 4 for x in chain.backbone)


def test_a_cycle_ignores_atom_order():
    mol = smiles(ALPHA_CD)
    ref = perceive(mol)
    rank = mol.canonical_order()
    for _ in range(5):
        other = smiles(format(mol, 'r'))
        got = perceive(other)
        r2 = other.canonical_order()
        assert [r2[x] for x in got.backbone] == [rank[x] for x in ref.backbone]
        assert got.beads == ref.beads


def test_a_mixture_is_asked_per_component():
    assert cycles(smiles(ALPHA_CD + '.[Na+].[Cl-]')) is None
    assert isinstance(cycles(smiles('C1COCCOCCOCCOCCOCCO1')), Circle)
