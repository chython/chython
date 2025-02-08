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
    mol = smiles('CO')  # methanol
    assert str(mol) in ('CO', 'OC')

    mol = smiles('c1ccccc1')  # benzene
    assert str(mol) == 'c1ccccc1'


def test_format_options():
    # Test different format options
    mol = smiles('C=1C=CC=CC=1')

    # Test asymmetric closures
    assert str(mol) == 'C=1C=CC=CC=1'
    assert format(mol, 'a') == 'C=1C=CC=CC1'

    assert format(mol, '!b') == 'C1CCCCC1'

    # Test disable stereo
    mol = smiles('C[C@H](O)CC')
    assert '@' in str(mol)
    assert '@' not in format(mol, '!s')

    mol = smiles('c1ccccc1')
    assert format(mol, 'A') == 'C:1:C:C:C:C:C:1'
    assert format(mol, 'Aa') == 'C:1:C:C:C:C:C1'
    assert format(mol, 'm') == '[cH:1]1[cH:2][cH:3][cH:4][cH:5][cH:6]1'
    assert format(mol, 'h') == '[cH]1[cH][cH][cH][cH][cH]1'

    assert format(mol, 'Ah') == '[CH]:1:[CH]:[CH]:[CH]:[CH]:[CH]:1'
    assert format(mol, 'Ah!b') == '[CH]1[CH][CH][CH][CH][CH]1'

    mol = smiles('[K+]')
    assert str(mol) == '[K+]'
    assert format(mol, '!z') == '[K]'

    mol = smiles('[CH3]')
    assert str(mol) == '[CH3] |^1:0|'
    assert format(mol, '!x') == '[CH3]'

    mol = smiles('CCO')
    assert len({format(mol, 'r') for _ in range(50)}) == 4


def test_smiles_comparison():
    # Test SMILES comparison functionality
    mol1 = smiles('CCO')
    mol2 = smiles('OCC')
    mol3 = smiles('CCC')

    assert mol1 == mol2  # same molecules
    assert mol1 != mol3  # different molecules
    assert hash(mol1) == hash(mol2)  # same hash for same molecules
