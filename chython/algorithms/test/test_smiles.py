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


def test_basic_smiles():
    # Test basic SMILES generation
    mol = smiles('CCO')  # ethanol
    assert 'C' in str(mol) and 'O' in str(mol)  # check presence of atoms
    
    mol = smiles('c1ccccc1')  # benzene
    assert 'c1ccccc1' in str(mol)  # aromatic representation


def test_format_options():
    # Test different format options
    mol = smiles('c1ccccc1')
    
    # Test asymmetric closures
    assert mol.__format__('a').startswith('c')
    
    # Test disable stereo
    chiral_mol = smiles('C[C@H](O)CC')
    assert '@' not in chiral_mol.__format__('!s')
    
    # Test aromatic bonds
    kekulized = mol.__format__('A')
    assert 'c' not in kekulized  # should not contain aromatic atoms
    
    # Test atom mapping
    assert ':' in mol.__format__('m')  # atom mapping numbers present
    
    # Test random ordering
    mol_str = str(mol)
    random_smiles = mol.__format__('r')
    assert len(random_smiles) > 0  # valid SMILES generated


def test_smiles_atoms_order():
    # Test atoms order property
    mol = smiles('CCO')
    order = mol.smiles_atoms_order
    assert isinstance(order, tuple)
    assert len(order) == 3  # number of atoms
    assert all(isinstance(x, int) for x in order)


def test_molecule_smiles():
    # Test MoleculeSmiles specific functionality
    mol = smiles('CCO')
    atoms = list(mol._atoms.keys())  # get actual atom indices
    
    # Test sticky smiles generation
    sticky = mol.sticky_smiles(atoms[0])  # fix first atom
    assert sticky and isinstance(sticky, str)
    
    # Test sticky smiles with both ends
    sticky_both = mol.sticky_smiles(atoms[0], atoms[-1])  # fix first and last atoms
    assert sticky_both and isinstance(sticky_both, str)


def test_complex_structures():
    # Test complex molecular structures
    mol = smiles('C1CC(=O)NC(=O)C1')  # cyclic peptide
    assert all(x in str(mol) for x in ('C', 'N', '=O'))  # check for expected fragments
    
    mol = smiles('C[C@H](N)C(=O)O')  # amino acid
    assert '@' in str(mol)  # stereo information preserved


def test_charged_species():
    # Test charged molecules
    mol = smiles('[NH4+]')  # ammonium
    assert '+' in str(mol)
    
    mol = smiles('[OH-]')  # hydroxide
    assert '-' in str(mol)


def test_radical_species():
    # Test radical species
    mol = smiles('[CH3]')
    assert '[' in str(mol) and ']' in str(mol)  # bracketed form
    
    # Test with format options
    assert '[' in mol.__format__('h')  # show hydrogens


def test_cgr_smiles():
    # Test CGR SMILES functionality
    mol = smiles('CC>>CCC')  # dynamic transformation
    assert '>' in str(mol)
    
    # Test dynamic bonds
    mol = smiles('C=C>>CC')
    assert '=' in str(mol)


def test_query_smiles():
    # Test basic query atoms
    mol = smiles('[C]')  # carbon atom
    assert len(mol) == 1
    
    mol = smiles('[N]')  # nitrogen atom
    assert len(mol) == 1
    
    mol = smiles('[O]')  # oxygen atom
    assert len(mol) == 1
    
    mol = smiles('[H]')  # hydrogen atom
    assert len(mol) == 1


def test_smiles_comparison():
    # Test SMILES comparison functionality
    mol1 = smiles('CCO')
    mol2 = smiles('CCO')
    mol3 = smiles('CCC')
    
    assert mol1 == mol2  # same molecules
    assert mol1 != mol3  # different molecules
    assert hash(mol1) == hash(mol2)  # same hash for same molecules


def test_cxsmiles_extensions():
    # Test CXSMILES extensions
    mol = smiles('[CH3]')  # radical
    assert mol.smiles  # valid SMILES generated
    
    # Test without CXSMILES
    assert mol.__format__('!x')  # valid SMILES without extensions


def test_special_cases():
    # Test special cases and edge cases
    mol = smiles('[H][H]')  # hydrogen molecule
    assert '[H]' in str(mol)
    
    mol = smiles('C#N')  # triple bond
    assert '#' in str(mol)
    
    mol = smiles('C~C')  # any bond
    assert '~' in str(mol) 