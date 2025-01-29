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
from chython.exceptions import IncorrectSmiles, IncorrectSmarts
from chython.files.daylight.parser import parser
from chython.files.daylight.tokenize import smiles_tokenize, smarts_tokenize
from chython.files.daylight.smarts import smarts


def test_parser_basic():
    # Test basic parsing
    result = parser(list(smiles_tokenize('CC')), True)
    assert len(result['atoms']) == 2
    assert len(result['bonds']) == 1
    
    result = parser(list(smiles_tokenize('O')), True)
    assert len(result['atoms']) == 1
    assert not result['bonds']


def test_parser_cycles():
    # Test cycle notation
    result = parser(list(smiles_tokenize('C1CCC1')), True)
    assert len(result['atoms']) == 4
    assert len(result['bonds']) == 4
    
    # Test multiple cycles
    result = parser(list(smiles_tokenize('C1CC2CCC12')), True)
    assert len(result['atoms']) == 6
    assert len(result['bonds']) == 7


def test_parser_branches():
    # Test branched structures
    result = parser(list(smiles_tokenize('CC(C)C')), True)
    assert len(result['atoms']) == 4
    assert len(result['bonds']) == 3
    
    # Test nested branches
    result = parser(list(smiles_tokenize('CC(C(C)C)C')), True)
    assert len(result['atoms']) == 6
    assert len(result['bonds']) == 5


def test_parser_errors():
    # Test invalid structures
    with pytest.raises(IncorrectSmiles):
        parser(list(smiles_tokenize('C1CC')), True)  # unclosed cycle
    
    with pytest.raises(IncorrectSmiles):
        parser(list(smiles_tokenize('C((C)')), True)  # unmatched parentheses
    
    with pytest.raises(IncorrectSmiles):
        parser(list(smiles_tokenize('C1CC2')), True)  # unclosed cycles


def test_parser_stereo():
    # Test stereo bonds
    result = parser(list(smiles_tokenize('C/C=C/C')), True)
    assert result['stereo_bonds']
    
    result = parser(list(smiles_tokenize('C/C=C\\C')), True)
    assert result['stereo_bonds']


def test_parser_aromatic():
    # Test aromatic structures
    result = parser(list(smiles_tokenize('c1ccccc1')), True)
    assert len(result['atoms']) == 6
    assert len(result['bonds']) == 6
    
    result = parser(list(smiles_tokenize('n1ccccc1')), True)
    assert len(result['atoms']) == 6
    assert len(result['bonds']) == 6


def test_parser_smarts():
    # Test SMARTS parsing with supported primitives
    result = parser(list(smarts_tokenize('[D2]')), False)
    assert len(result['atoms']) == 1
    
    result = parser(list(smarts_tokenize('[h1]')), False)
    assert len(result['atoms']) == 1
    
    result = parser(list(smarts_tokenize('[!R]')), False)
    assert len(result['atoms']) == 1
    
    result = parser(list(smarts_tokenize('[r5]')), False)
    assert len(result['atoms']) == 1


def test_parser_smarts_errors():
    # Test invalid SMARTS
    with pytest.raises(IncorrectSmiles):
        parser(list(smarts_tokenize('[')), False)  # unclosed bracket
    
    with pytest.raises(IncorrectSmiles):
        parser(list(smarts_tokenize('[C')), False)  # unclosed bracket
    
    with pytest.raises(IncorrectSmiles):
        parser(list(smarts_tokenize('C1C')), False)  # unclosed ring 