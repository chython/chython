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


def test_basic():
    # Test basic atom mapping in simple molecules
    mol1 = smiles('CC(=O)O')  # acetic acid
    mol2 = smiles('CC(O)=O')  # acetic acid

    assert mol1 <= mol2
    assert mol2 <= mol1

    mappings = list(mol1.get_mapping(mol2))
    assert len(mappings) == 1
    assert mappings[0] == {1: 1, 2: 2, 3: 4, 4: 3}
    assert not smiles('CC(O)O') <= mol1
    assert smiles('C[O-]') <= smiles('CC[O-]')
    assert not smiles('C[O-]') <= smiles('CCO')


def test_substructure_mapping():
    # Test mapping of a substructure
    mol = smiles('CCC(=O)OC')
    substructure = smiles('CC(=O)O')

    assert substructure < mol
    mappings = list(substructure.get_mapping(mol))
    assert len(mappings) == 1
    assert mappings[0] == {1: 2, 2: 3, 3: 4, 4: 5}


def test_multiple_mappings():
    # Test cases where multiple valid mappings exist
    mol = smiles('CC(=O)OC(=O)C')
    pattern = smiles('CC(=O)O')  # acetone pattern

    mappings = list(pattern.get_mapping(mol))
    assert len(mappings) == 2  # should find multiple matches
    assert {1: 1, 2: 2, 3: 3, 4: 4} in mappings
    assert {1: 7, 2: 5, 3: 6, 4: 4} in mappings
