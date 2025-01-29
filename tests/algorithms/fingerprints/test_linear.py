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
import numpy as np
import pytest


def test_linear_fingerprint_basic():
    # Test basic fingerprint generation
    mol = smiles('CCO')
    fp = mol.linear_fingerprint(min_radius=1, max_radius=2, length=1024)
    
    # Test array properties
    assert isinstance(fp, np.ndarray)
    assert fp.dtype == np.uint8
    assert fp.shape == (1024,)
    
    # Test binary nature
    assert set(np.unique(fp)).issubset({0, 1})


def test_linear_fingerprint_consistency():
    # Test that fingerprints are consistent for the same molecule
    mol = smiles('CCO')
    fp1 = mol.linear_fingerprint()
    fp2 = mol.linear_fingerprint()
    
    # Test exact equality of arrays
    np.testing.assert_array_equal(fp1, fp2)
    
    # Test different molecules give different fingerprints
    mol2 = smiles('CCC')
    fp3 = mol2.linear_fingerprint()
    assert not np.array_equal(fp1, fp3)


def test_linear_fingerprint_parameters():
    mol = smiles('CCO')
    
    # Test different radius parameters
    fp1 = mol.linear_fingerprint(min_radius=1, max_radius=2)
    fp2 = mol.linear_fingerprint(min_radius=1, max_radius=3)
    assert fp2.sum() >= fp1.sum()  # More radius should capture more features
    
    # Test different lengths
    fp3 = mol.linear_fingerprint(length=2048)
    assert fp3.shape == (2048,)
    assert isinstance(fp3, np.ndarray)
    assert fp3.dtype == np.uint8
    
    # Test number of active bits
    fp4 = mol.linear_fingerprint(number_active_bits=3)
    assert fp4.sum() >= fp1.sum()  # More active bits should set more bits


def test_linear_fingerprint_bit_pairs():
    # Test the number_bit_pairs parameter
    mol = smiles('CCCC')  # molecule with multiple similar fragments
    
    # Compare different number_bit_pairs settings
    fp1 = mol.linear_fingerprint(number_bit_pairs=1)
    fp2 = mol.linear_fingerprint(number_bit_pairs=2)
    fp3 = mol.linear_fingerprint(number_bit_pairs=4)
    
    # More bit pairs should potentially activate more bits
    assert fp1.sum() <= fp2.sum() <= fp3.sum()


def test_linear_fingerprint_complex_molecule():
    # Test with a more complex molecule
    mol = smiles('c1ccccc1CC(=O)O')
    fp = mol.linear_fingerprint()
    
    # Basic checks
    assert isinstance(fp, np.ndarray)
    assert fp.dtype == np.uint8
    
    # Should have reasonable number of bits set
    assert 0 < fp.sum() < len(fp)  # some bits should be set, but not all
    
    # Test with different parameters
    fp_large = mol.linear_fingerprint(max_radius=6, length=2048)
    assert fp_large.shape == (2048,)
    assert fp_large.sum() > 0


def test_linear_fingerprint_edge_cases():
    # Test single atom
    mol_single = smiles('C')
    fp_single = mol_single.linear_fingerprint()
    assert isinstance(fp_single, np.ndarray)
    assert fp_single.dtype == np.uint8
    assert fp_single.sum() > 0  # should have some bits set
    
    # Test disconnected components
    mol_disconnected = smiles('CC.CC')
    fp_disconnected = mol_disconnected.linear_fingerprint()
    assert isinstance(fp_disconnected, np.ndarray)
    assert fp_disconnected.dtype == np.uint8


def test_linear_fingerprint_arbitrary_length():
    # Test that non-power-of-2 lengths work but might have unexpected behavior
    mol = smiles('CCO')
    lengths = [100, 1000, 1500, 3000]
    
    for length in lengths:
        fp = mol.linear_fingerprint(length=length)
        assert isinstance(fp, np.ndarray)
        assert fp.dtype == np.uint8
        assert fp.shape == (length,)
        # The actual bits set might be fewer than expected due to masking
        assert 0 <= fp.sum() <= length


def test_linear_fingerprint_comparison():
    # Test fingerprint comparison between similar molecules
    mol1 = smiles('CCO')
    mol2 = smiles('CCC')
    mol3 = smiles('CCCO')
    
    fp1 = mol1.linear_fingerprint()
    fp2 = mol2.linear_fingerprint()
    fp3 = mol3.linear_fingerprint()
    
    # Calculate Tanimoto similarities
    def tanimoto(a, b):
        intersection = np.sum(np.logical_and(a, b))
        union = np.sum(np.logical_or(a, b))
        return intersection / union if union > 0 else 0.0
    
    # Similar molecules should have higher similarity
    sim12 = tanimoto(fp1, fp2)
    sim13 = tanimoto(fp1, fp3)
    sim23 = tanimoto(fp2, fp3)
    
    assert 0 <= sim12 <= 1
    assert 0 <= sim13 <= 1
    assert 0 <= sim23 <= 1 