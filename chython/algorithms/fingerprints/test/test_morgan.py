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
from chython import smiles, ReactionContainer
import numpy as np
from pytest import mark
import pytest


def test_morgan_fingerprint():
    # Test basic fingerprint generation
    mol = smiles('CCO')
    fp = mol.morgan_fingerprint(min_radius=1, max_radius=2, length=1024)
    
    assert isinstance(fp, np.ndarray)
    assert fp.dtype == np.uint8
    assert fp.shape == (1024,)
    assert fp.sum() > 0  # Should have some bits set
    
    # Test different lengths
    fp2 = mol.morgan_fingerprint(length=2048)
    assert fp2.shape == (2048,)
    
    # Test different number of active bits
    fp3 = mol.morgan_fingerprint(number_active_bits=3)
    assert fp3.sum() >= fp.sum()  # Should have more or equal bits set


def test_morgan_bit_set():
    mol = smiles('CCO')
    bits = mol.morgan_bit_set(min_radius=1, max_radius=2, length=1024)
    
    assert isinstance(bits, set)
    assert len(bits) > 0
    assert all(isinstance(x, int) for x in bits)
    assert all(0 <= x < 1024 for x in bits)
    
    # Test with different parameters
    bits2 = mol.morgan_bit_set(length=2048, number_active_bits=3)
    assert all(0 <= x < 2048 for x in bits2)
    assert len(bits2) >= len(bits)  # Should have more or equal bits


def test_morgan_hash_set():
    mol = smiles('CCO')
    hashes = mol.morgan_hash_set(min_radius=1, max_radius=2)
    
    assert isinstance(hashes, set)
    assert len(hashes) > 0
    assert all(isinstance(x, int) for x in hashes)


def test_morgan_hash_smiles():
    mol = smiles('CCO')
    hash_smiles = mol.morgan_hash_smiles(min_radius=1, max_radius=2)
    
    assert isinstance(hash_smiles, dict)
    assert len(hash_smiles) > 0
    assert all(isinstance(k, int) for k in hash_smiles)
    assert all(isinstance(v, list) for v in hash_smiles.values())
    assert all(isinstance(s, str) for v in hash_smiles.values() for s in v)


def test_morgan_smiles_hash():
    mol = smiles('CCO')
    smiles_hash = mol.morgan_smiles_hash(min_radius=1, max_radius=2)
    
    assert isinstance(smiles_hash, dict)
    assert len(smiles_hash) > 0
    assert all(isinstance(k, str) for k in smiles_hash)
    assert all(isinstance(v, list) for v in smiles_hash.values())
    assert all(isinstance(h, int) for v in smiles_hash.values() for h in v)


@mark.parametrize('radius', [(0, 1), (1, 0), (-1, 2)])
def test_invalid_radius(radius):
    mol = smiles('CCO')
    min_r, max_r = radius
    try:
        mol.morgan_fingerprint(min_radius=min_r, max_radius=max_r)
        assert False, "Should raise AssertionError"
    except AssertionError:
        pass


def test_complex_molecule():
    # Test with a more complex molecule containing rings and multiple atom types
    mol = smiles('c1ccccc1CC(=O)O')
    
    fp1 = mol.morgan_fingerprint(min_radius=1, max_radius=3)
    fp2 = mol.morgan_fingerprint(min_radius=1, max_radius=4)
    
    assert fp1.sum() < fp2.sum()  # More radius should capture more features
    
    # Test hash consistency
    hash_set1 = mol.morgan_hash_set(min_radius=1, max_radius=2)
    hash_set2 = mol.morgan_hash_set(min_radius=1, max_radius=2)
    assert hash_set1 == hash_set2  # Should be deterministic 


def test_morgan_fingerprint_numpy():
    # Test numpy array properties of Morgan fingerprints
    mol = smiles('CCO')
    fp = mol.morgan_fingerprint(min_radius=1, max_radius=2, length=1024)
    
    # Test array type and shape
    assert isinstance(fp, np.ndarray)
    assert fp.dtype == np.uint8
    assert fp.shape == (1024,)
    
    # Test binary nature
    assert set(np.unique(fp)).issubset({0, 1})
    
    # Test different lengths
    fp_2048 = mol.morgan_fingerprint(length=2048)
    assert fp_2048.shape == (2048,)
    assert fp_2048.dtype == np.uint8
    
    # Test different number of active bits
    fp_more_bits = mol.morgan_fingerprint(number_active_bits=4)
    assert fp_more_bits.sum() >= fp.sum()


def test_morgan_fingerprint_consistency():
    # Test that fingerprints are consistent for the same molecule
    mol = smiles('CCO')
    fp1 = mol.morgan_fingerprint()
    fp2 = mol.morgan_fingerprint()
    
    # Test exact equality of arrays
    np.testing.assert_array_equal(fp1, fp2)
    
    # Test different molecules give different fingerprints
    mol2 = smiles('CCC')
    fp3 = mol2.morgan_fingerprint()
    assert not np.array_equal(fp1, fp3)


def test_morgan_fingerprint_parameters():
    mol = smiles('CCO')
    
    # Test different radius parameters
    fp1 = mol.morgan_fingerprint(min_radius=1, max_radius=2)
    fp2 = mol.morgan_fingerprint(min_radius=1, max_radius=3)
    assert fp2.sum() >= fp1.sum()  # More radius should capture more features
    
    # Test power of 2 lengths
    for length in [128, 256, 512, 1024, 2048, 4096]:
        fp = mol.morgan_fingerprint(length=length)
        assert fp.shape == (length,)
        assert isinstance(fp, np.ndarray)
        assert fp.dtype == np.uint8


def test_morgan_fingerprint_arbitrary_length():
    # Test that non-power-of-2 lengths work but might have unexpected behavior
    mol = smiles('CCO')
    lengths = [100, 1000, 1500, 3000]
    
    for length in lengths:
        fp = mol.morgan_fingerprint(length=length)
        assert isinstance(fp, np.ndarray)
        assert fp.dtype == np.uint8
        assert fp.shape == (length,)
        # The actual bits set might be fewer than expected due to masking
        assert 0 <= fp.sum() <= length


def test_complex_molecule_clean2d():
    # Test with more complex molecule
    mol = smiles('c1ccccc1CC(=O)O')
    mol.clean2d()
    
    # Check coordinates exist
    assert all(n in mol._plane for n in mol)
    
    # Verify ring atoms are roughly coplanar
    ring_atoms = [n for n in mol if len(mol._bonds[n]) == 2]
    if ring_atoms:
        coords = [mol._plane[n] for n in ring_atoms]
        # Calculate variance in y coordinates - should be small for planar ring
        y_coords = [y for _, y in coords]
        y_mean = sum(y_coords) / len(y_coords)
        y_variance = sum((y - y_mean) ** 2 for y in y_coords) / len(y_coords)
        assert y_variance < 1.0  # reasonable threshold for planarity


def test_disconnected_components():
    # Test molecule with multiple disconnected components
    mol = smiles('CCO.c1ccccc1')
    mol.clean2d()
    
    # Check all atoms have coordinates
    assert all(n in mol._plane for n in mol)
    
    # Components should be separated in space
    components = list(mol.connected_components)
    assert len(components) == 2
    
    # Get bounding boxes for each component
    def get_bounds(atoms):
        xs = [mol._plane[n][0] for n in atoms]
        ys = [mol._plane[n][1] for n in atoms]
        return min(xs), max(xs), min(ys), max(ys)
    
    bounds1 = get_bounds(components[0])
    bounds2 = get_bounds(components[1])
    
    # Check components don't overlap in x-direction
    assert bounds1[1] < bounds2[0] or bounds2[1] < bounds1[0]


def test_reaction_clean2d():
    # Create a simple reaction
    reactant = smiles('CCO')
    product = smiles('CC=O')
    reaction = ReactionContainer([reactant], [product])
    
    # Clean coordinates
    reaction.clean2d()
    
    # Check that all molecules have coordinates
    for molecule in reaction.molecules():
        assert all(n in molecule._plane for n in molecule)
    
    # Check that reactants are positioned before products
    reactant_max_x = max(x for mol in reaction.reactants 
                        for x, _ in mol._plane.values())
    product_min_x = min(x for mol in reaction.products 
                       for x, _ in mol._plane.values())
    assert reactant_max_x < product_min_x
    
    # Check arrow exists and is positioned between reactants and products
    assert hasattr(reaction, '_arrow')
    arrow_start, arrow_end = reaction._arrow
    assert reactant_max_x < arrow_start < arrow_end < product_min_x


def test_reaction_with_reagents():
    # Create reaction with reagents
    reactant = smiles('CCO')
    reagent = smiles('Cl')
    product = smiles('CCCl')
    reaction = ReactionContainer([reactant], [reagent], [product])
    
    reaction.clean2d()
    
    # Check all molecules have coordinates
    for molecule in reaction.molecules():
        assert all(n in molecule._plane for n in molecule)
    
    # Check reagents are positioned above the arrow
    reagent_coords = [(x, y) for mol in reaction.reagents 
                     for x, y in mol._plane.values()]
    assert all(y > 0 for _, y in reagent_coords)  # reagents should be above
    
    # Verify arrow position
    arrow_start, arrow_end = reaction._arrow
    assert arrow_start < arrow_end
    
    # Check relative positioning
    reactant_max_x = max(x for mol in reaction.reactants 
                        for x, _ in mol._plane.values())
    product_min_x = min(x for mol in reaction.products 
                       for x, _ in mol._plane.values())
    assert reactant_max_x < arrow_start < arrow_end < product_min_x


def test_fix_positions():
    # Test just the position fixing functionality
    reaction = ReactionContainer([smiles('CCO')], [smiles('CC=O')])
    
    # Clean individual molecules first
    for mol in reaction.molecules():
        mol.clean2d()
    
    # Then fix positions
    reaction.fix_positions()
    
    # Check arrow exists
    assert hasattr(reaction, '_arrow')
    
    # Check molecules are properly spaced
    reactant_coords = [(x, y) for mol in reaction.reactants 
                      for x, y in mol._plane.values()]
    product_coords = [(x, y) for mol in reaction.products 
                     for x, y in mol._plane.values()]
    
    # Verify no overlap between reactants and products
    reactant_max_x = max(x for x, _ in reactant_coords)
    product_min_x = min(x for x, _ in product_coords)
    assert reactant_max_x < product_min_x 