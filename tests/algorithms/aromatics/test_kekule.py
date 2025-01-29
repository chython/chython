# -*- coding: utf-8 -*-
#
#  Copyright 2025 Ramil Nugmanov <nougmanoff@protonmail.com>
#  Copyright 2025 Tagir Akhmetshin <tagirshin@gmail.com>
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
from chython import smiles
import pytest
from chython.exceptions import InvalidAromaticRing


def test_kekule_basic():
    # Test basic aromatic ring conversion
    mol = smiles('c1ccccc1')  # benzene
    assert mol.kekule()  # should return True for aromatic rings
    
    # Verify alternating single and double bonds
    bonds = mol._bonds
    double_bonds = sum(1 for n, ms in bonds.items() for m, b in ms.items() if b.order == 2 and m > n)
    assert double_bonds == 3  # benzene should have 3 double bonds


def test_kekule_pyridine():
    # Test pyridine and its derivatives
    mol = smiles('n1ccccc1')  # pyridine
    assert mol.kekule()
    
    # Test protonated pyridine
    mol_protonated = smiles('[nH+]1ccccc1')
    assert mol_protonated.kekule()


def test_kekule_pyrrole():
    # Test pyrrole and its derivatives
    mol = smiles('[nH]1cccc1')  # pyrrole
    assert mol.kekule()
    
    # Test N-methylpyrrole
    mol_methyl = smiles('Cn1cccc1')
    assert mol_methyl.kekule()


def test_kekule_furan_thiophene():
    # Test oxygen and sulfur containing aromatics
    mol_furan = smiles('o1cccc1')
    assert mol_furan.kekule()
    
    mol_thiophene = smiles('s1cccc1')
    assert mol_thiophene.kekule()


def test_kekule_complex_systems():
    # Test fused ring systems
    mol_naphthalene = smiles('c1ccc2ccccc2c1')
    assert mol_naphthalene.kekule()
    
    # Test indole
    mol_indole = smiles('c1ccc2[nH]ccc2c1')
    assert mol_indole.kekule()


def test_kekule_enumeration():
    # Test enumeration of Kekulé structures
    mol = smiles('c1ccccc1')  # benzene
    forms = list(mol.enumerate_kekule())
    assert len(forms) == 2  # benzene has 2 Kekulé forms


def test_kekule_invalid_structures():
    # Test invalid aromatic structures
    with pytest.raises(InvalidAromaticRing):
        mol = smiles('c1cccc1')  # 5-membered carbon ring (invalid aromatic)
        mol.kekule()
    
    with pytest.raises(InvalidAromaticRing):
        mol = smiles('c1ccc2c1c3ccccc3cc2')  # acenaphthalene (invalid aromatic form)
        mol.kekule()
    
    with pytest.raises(InvalidAromaticRing):
        mol = smiles('c1cccc1C(=O)c1cccc1')  # cyclopentadiene with carbonyl (invalid aromatic)
        mol.kekule()


def test_kekule_charged_species():
    # Test charged aromatic species
    mol_pyridinium = smiles('[n+]1ccccc1')
    assert mol_pyridinium.kekule()
    
    mol_cyclopentadienyl = smiles('[cH-]1cccc1')
    assert mol_cyclopentadienyl.kekule()


def test_kekule_multiple_rings():
    # Test molecules with multiple aromatic rings
    mol_biphenyl = smiles('c1ccccc1-c2ccccc2')
    assert mol_biphenyl.kekule()
    
    # Test phenylpyridine
    mol_phenylpyridine = smiles('c1ccccc1-c2ccccn2')
    assert mol_phenylpyridine.kekule()


def test_kekule_heteroatoms():
    # Test various heteroatoms in aromatic rings
    mol_pyrazine = smiles('n1ccncc1')  # two nitrogens
    assert mol_pyrazine.kekule()
    
    mol_oxazole = smiles('o1cncc1')  # oxygen and nitrogen
    assert mol_oxazole.kekule()
    
    mol_thiazole = smiles('s1cncc1')  # sulfur and nitrogen
    assert mol_thiazole.kekule()


def test_kekule_buffer_size():
    # Test buffer size parameter for complex heterocycles
    mol1 = smiles('c1ccc2[nH]ccc2c1')  # indole
    assert mol1.kekule(buffer_size=1)  # small buffer
    
    mol2 = smiles('c1ccc2[nH]ccc2c1')  # fresh indole instance
    assert mol2.kekule(buffer_size=10)  # large buffer


def test_kekule_radical_species():
    # Test radical aromatic species
    mol_phenoxy = smiles('[O]c1ccccc1')
    assert mol_phenoxy.kekule()
    
    # Test radical cation
    mol_benzene_radical = smiles('[c]1ccccc1')
    assert mol_benzene_radical.kekule()


def test_kekule_quinones():
    # Test quinone-like structures
    mol_benzoquinone = smiles('O=C1C=CC(=O)C=C1')
    assert not mol_benzoquinone.kekule()  # not aromatic
    
    # Test semiquinone
    mol_semiquinone = smiles('O=C1C=CC(O)C=C1')
    assert not mol_semiquinone.kekule()  # not aromatic 