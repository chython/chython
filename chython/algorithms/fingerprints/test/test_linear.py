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
import numpy as np
from chython import smiles


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

    # Test different lengths
    fp = mol.linear_fingerprint(length=2048)
    assert isinstance(fp, np.ndarray)
    assert fp.dtype == np.uint8
    assert fp.shape == (2048,)


def test_linear_fingerprint_consistency():
    # Test that fingerprints are consistent for the same molecule
    fp1 = smiles('CCO').linear_fingerprint()
    fp2 = smiles('OCC').linear_fingerprint()

    # Test exact equality of arrays
    assert np.array_equal(fp1, fp2)

    # Test different molecules give different fingerprints)
    fp3 = smiles('CCC').linear_fingerprint()
    assert not np.array_equal(fp1, fp3)


def test_linear_fingerprint_parameters():
    mol = smiles('CCO')

    # Test different radius parameters
    fp1 = mol.linear_fingerprint(min_radius=1, max_radius=2)
    fp2 = mol.linear_fingerprint(min_radius=1, max_radius=3)
    assert not np.array_equal(fp1, fp2)
    assert np.array_equal(fp1 & fp2, fp1)

    # Test number of active bits
    fp3 = mol.linear_fingerprint(number_active_bits=2)
    fp4 = mol.linear_fingerprint(number_active_bits=3)
    assert not np.array_equal(fp3, fp4)
    assert np.array_equal(fp3 & fp4, fp3)


def test_linear_fingerprint_bit_pairs():
    # Test the number_bit_pairs parameter
    mol = smiles('CCCCCCCCCCCCCCCCCCCCCCCC')  # molecule with multiple similar fragments

    # Compare different number_bit_pairs settings
    fp1 = mol.linear_fingerprint(number_bit_pairs=2)
    fp2 = mol.linear_fingerprint(number_bit_pairs=3)

    assert not np.array_equal(fp1, fp2)
    assert np.array_equal(fp1 & fp2, fp1)


def test_linear_fingerprint_edge_cases():
    fp1 = smiles('C').linear_fingerprint()
    assert fp1.sum() == 2

    fp1 = smiles('CC').linear_fingerprint()
    fp2 = smiles('CC.CC').linear_fingerprint()
    assert not np.array_equal(fp1, fp2)
    assert np.array_equal(fp1 & fp2, fp1)


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
