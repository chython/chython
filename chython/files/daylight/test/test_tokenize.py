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
from chython.files.daylight.tokenize import smiles_tokenize, smarts_tokenize
from chython.containers import QueryBond
from chython.exceptions import IncorrectSmiles
from pytest import raises


def test_smiles_tokenize():
    assert smiles_tokenize('C') == [(0, {'element': 'C'})]
    assert smiles_tokenize('CC') == [(0, {'element': 'C'}), (0, {'element': 'C'})]
    assert smiles_tokenize('C=O') == [(0, {'element': 'C'}), (1, 2), (0, {'element': 'O'})]
    assert smiles_tokenize('C(O)N') == [(0, {'element': 'C'}), (2, None), (0, {'element': 'O'}),
                                        (3, None), (0, {'element': 'N'})]
    assert smiles_tokenize('C2CC2') == [(0, {'element': 'C'}), (6, 2), (0, {'element': 'C'}),
                                        (0, {'element': 'C'}), (6, 2)]


def test_smiles_tokenize_atom():
    assert smiles_tokenize('[NH4+]') == [(0, {'element': 'N', 'isotope': None, 'parsed_mapping': None, 'charge': 1,
                                              'implicit_hydrogens': 4, 'stereo': None})]
    assert smiles_tokenize('[14N]') == [(0, {'element': 'N', 'isotope': 14, 'parsed_mapping': None, 'charge': 0,
                                              'implicit_hydrogens': 0, 'stereo': None})]
    assert smiles_tokenize('[N@H]') == [(0, {'element': 'N', 'isotope': None, 'parsed_mapping': None, 'charge': 0,
                                             'implicit_hydrogens': 1, 'stereo': True})]
    assert smiles_tokenize('[N@@H--]') == [(0, {'element': 'N', 'isotope': None, 'parsed_mapping': None, 'charge': -2,
                                                'implicit_hydrogens': 1, 'stereo': False})]
    assert smiles_tokenize('[N@+3]') == [(0, {'element': 'N', 'isotope': None, 'parsed_mapping': None, 'charge': 3,
                                              'implicit_hydrogens': 0, 'stereo': True})]
    assert smiles_tokenize('[CH2:2]') == [(0, {'element': 'C', 'isotope': None, 'parsed_mapping': 2, 'charge': 0,
                                               'implicit_hydrogens': 2, 'stereo': None})]
    with raises(IncorrectSmiles):
        smiles_tokenize('[@N]')


def test_smarts_tokenize_atom():
    # Test basic SMARTS tokenization
    assert smarts_tokenize('[C]') == [(0, {'element': 'C'})]
    assert smarts_tokenize('[C,N]') == [(0, {'element': ['C', 'N']})]
    assert smarts_tokenize('[C+]') == [(0, {'charge': 1, 'element': 'C'})]
    assert smarts_tokenize('[#1]') == [(0, {'element': 1})]
    assert smarts_tokenize('[C;h1;@]') == [(0, {'element': 'C', 'implicit_hydrogens': [1], 'stereo': True})]
    assert smarts_tokenize('[O;z1,z2;x1]') == [(0, {'element': 'O', 'heteroatoms': [1], 'hybridization': [1, 2]})]
    assert smarts_tokenize('[Se;a;D1,D2;r4,r7:3]') == [(0, {'parsed_mapping': 3, 'element': 'Se', 'hybridization': 4, 'neighbors': [1, 2], 'ring_sizes': [4, 7]})]
    assert smarts_tokenize('[Cl;M]') == [(0, {'element': 'Cl', 'masked': True})]
    assert smarts_tokenize('[A:1]') == [(0, {'parsed_mapping': 1, 'element': 'A'})]
    assert smarts_tokenize('[M]') == [(0, {'element': 'M'})]


def test_smarts_tokenize_bonds():
    assert smarts_tokenize('[C][C]') == [(0, {'element': 'C'}), (0, {'element': 'C'})]
    assert smarts_tokenize('[C]-[C]') == [(0, {'element': 'C'}), (1, 1), (0, {'element': 'C'})]
    assert smarts_tokenize('[C]~[C]') == [(0, {'element': 'C'}), (1, 8), (0, {'element': 'C'})]
    assert smarts_tokenize('[C]!:[C]') == [(0, {'element': 'C'}), (10, [1, 2, 3]), (0, {'element': 'C'})]
    assert smarts_tokenize('[C]-,=[C]') == [(0, {'element': 'C'}), (10, [1, 2]), (0, {'element': 'C'})]
    assert smarts_tokenize('[C]-;@[C]') == [(0, {'element': 'C'}), (12, QueryBond(1, True)), (0, {'element': 'C'})]
    assert smarts_tokenize('[C]!-;!@[C]') == [(0, {'element': 'C'}), (12, QueryBond((2, 3, 4), False)),
                                              (0, {'element': 'C'})]
    assert smarts_tokenize('[C]-,=;!@[C]') == [(0, {'element': 'C'}), (12, QueryBond((1, 2), False)),
                                               (0, {'element': 'C'})]
