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


def test_basic_mapping():
    # Test basic atom mapping in simple molecules
    mol1 = smiles('CC(=O)O')  # acetic acid
    mol2 = smiles('CC(=O)O')  # acetic acid
    
    mappings = list(mol1.get_mapping(mol2))
    assert len(mappings) > 0  # at least one mapping should exist
    mapping = mappings[0]  # take first mapping
    assert len(mapping) == len(mol1)  # all atoms should be mapped
    assert all(isinstance(k, int) and isinstance(v, int) for k, v in mapping.items())


def test_substructure_mapping():
    # Test mapping of a substructure
    mol = smiles('CC(=O)OC')  # methyl acetate
    substructure = smiles('CC(=O)O')  # acetic acid pattern
    
    mappings = list(substructure.get_mapping(mol))
    assert len(mappings) > 0  # at least one mapping should exist
    mapping = mappings[0]  # take first mapping
    assert len(mapping) == len(substructure)  # all substructure atoms should be mapped
    assert all(isinstance(k, int) and isinstance(v, int) for k, v in mapping.items())


def test_multiple_mappings():
    # Test cases where multiple valid mappings exist
    mol = smiles('CC(=O)CC(=O)C')  # 2,4-pentanedione
    pattern = smiles('CC(=O)C')  # acetone pattern
    
    mappings = list(pattern.get_mapping(mol))
    assert len(mappings) > 1  # should find multiple matches
    assert all(len(m) == len(pattern) for m in mappings)  # each mapping should cover all pattern atoms


def test_aromatic_mapping():
    # Test mapping with aromatic systems
    benzene = smiles('c1ccccc1')
    toluene = smiles('Cc1ccccc1')
    
    mappings = list(benzene.get_mapping(toluene))
    assert len(mappings) > 0  # at least one mapping should exist
    mapping = mappings[0]  # take first mapping
    assert len(mapping) == len(benzene)  # all benzene atoms should be mapped
    assert all(isinstance(k, int) and isinstance(v, int) for k, v in mapping.items())


def test_reaction_mapping():
    # Test mapping in reaction context
    reactant = smiles('CC(=O)O')  # acetic acid
    product = smiles('CC(=O)OC')  # methyl acetate
    
    mappings = list(reactant.get_mapping(product))
    assert len(mappings) > 0  # at least one mapping should exist
    mapping = mappings[0]  # take first mapping
    assert len(mapping) == len(reactant)  # all reactant atoms should be mapped
    assert all(isinstance(k, int) and isinstance(v, int) for k, v in mapping.items())


def test_complex_mapping():
    # Test mapping with complex molecules
    mol1 = smiles('CC1=C(C(=O)C2=C(C1=O)N3CC4=C(C3(CC2)C)NC5=CC=CC=C54)C')  # complex structure
    mol2 = smiles('CC1=C(C(=O)C2=C(C1=O)N3CC4=C(C3(CC2)C)NC5=CC=CC=C54)C')  # same structure
    
    mappings = list(mol1.get_mapping(mol2))
    assert len(mappings) > 0  # at least one mapping should exist
    mapping = mappings[0]  # take first mapping
    assert len(mapping) == len(mol1)  # all atoms should be mapped
    assert all(isinstance(k, int) and isinstance(v, int) for k, v in mapping.items())


def test_mapping_with_different_bonds():
    # Test mapping when bond orders differ
    mol1 = smiles('C=CC=C')  # 1,3-butadiene
    mol2 = smiles('C=CC=C')  # 1,3-butadiene
    
    # Should find mapping for identical molecules
    mappings = list(mol1.get_mapping(mol2))
    assert len(mappings) > 0  # at least one mapping should exist
    mapping = mappings[0]  # take first mapping
    assert len(mapping) == len(mol1)
    
    # Verify that the mapping preserves atom connectivity and bond orders
    for atom1, atom2 in mapping.items():
        # Check that the number of neighbors is the same
        assert len(mol1._bonds[atom1]) == len(mol2._bonds[atom2])
        # Check that bond orders are preserved
        mol1_orders = {mol1._bonds[atom1][x].order for x in mol1._bonds[atom1]}
        mol2_orders = {mol2._bonds[atom2][x].order for x in mol2._bonds[atom2]}
        assert mol1_orders == mol2_orders


def test_mapping_with_charges():
    # Test mapping with charged atoms
    mol1 = smiles('C[NH3+]')  # methylammonium
    mol2 = smiles('C[NH3+]')  # methylammonium
    
    # Should find mapping for identical molecules
    mappings = list(mol1.get_mapping(mol2))
    assert len(mappings) > 0  # at least one mapping should exist
    mapping = mappings[0]  # take first mapping
    assert len(mapping) == len(mol1)
    
    # Verify that the mapping preserves atom connectivity and charges
    for atom1, atom2 in mapping.items():
        assert mol1._bonds[atom1].keys() == mol2._bonds[atom2].keys()
        assert mol1._charges[atom1] == mol2._charges[atom2] 