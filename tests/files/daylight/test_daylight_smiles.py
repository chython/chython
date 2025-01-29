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
import pytest
from chython.exceptions import IncorrectSmiles
from chython.files.daylight.parser import parser
from chython.files.daylight.tokenize import smiles_tokenize


def test_daylight_smiles_basic():
    # Test basic SMILES
    result = parser(list(smiles_tokenize('CC')), True)
    assert len(result['atoms']) == 2
    assert len(result['bonds']) == 1
    
    result = parser(list(smiles_tokenize('O')), True)
    assert len(result['atoms']) == 1
    assert not result['bonds']


def test_daylight_smiles_empty():
    # Test empty SMILES
    tokens = list(smiles_tokenize(''))
    assert len(tokens) == 0  # empty string should produce empty token list
    
    # Empty token list should raise error when parsing
    with pytest.raises(IndexError):
        parser(tokens, True)


def test_daylight_smiles_invalid():
    # Test invalid SMILES
    with pytest.raises(IncorrectSmiles):
        parser(list(smiles_tokenize('C1CC')), True)  # unclosed cycle


def test_daylight_smiles_complex():
    # Test complex SMILES
    result = parser(list(smiles_tokenize('C1=CC=CC=C1')), True)
    assert len(result['atoms']) == 6
    assert len(result['bonds']) == 6
    
    result = parser(list(smiles_tokenize('C(=O)O')), True)
    assert len(result['atoms']) == 3
    assert len(result['bonds']) == 2


def test_daylight_smiles_charged():
    # Test charged species
    result = parser(list(smiles_tokenize('[NH4+]')), True)
    assert len(result['atoms']) == 1
    
    result = parser(list(smiles_tokenize('[OH-]')), True)
    assert len(result['atoms']) == 1


def test_daylight_smiles_aromatic():
    # Test aromatic SMILES
    result = parser(list(smiles_tokenize('c1ccccc1')), True)
    assert len(result['atoms']) == 6
    assert len(result['bonds']) == 6
    
    result = parser(list(smiles_tokenize('n1cccc1')), True)
    assert len(result['atoms']) == 5
    assert len(result['bonds']) == 5


def test_daylight_smiles_isotopes():
    # Test isotope labels
    result = parser(list(smiles_tokenize('[2H]O')), True)
    assert len(result['atoms']) == 2
    assert result['atoms'][0]['isotope'] == 2
    
    result = parser(list(smiles_tokenize('[13C]')), True)
    assert len(result['atoms']) == 1
    assert result['atoms'][0]['isotope'] == 13


def test_daylight_smiles_stereo():
    # Test stereochemistry
    result = parser(list(smiles_tokenize('C/C=C/C')), True)
    assert result['stereo_bonds']
    
    result = parser(list(smiles_tokenize('C/C=C\\C')), True)
    assert result['stereo_bonds'] 