# -*- coding: utf-8 -*-
#
#  Copyright 2026 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from chython import smiles


pytest.importorskip('openclatura')


def test_iupac_ethanol():
    assert smiles('CCO').iupac == 'ethanol'


def test_iupac_benzene():
    assert smiles('c1ccccc1').iupac == 'benzene'


def test_iupac_acetic_acid():
    assert smiles('CC(=O)O').iupac == 'acetic acid'


def test_iupac_cached():
    mol = smiles('CCO')
    first = mol.iupac
    assert 'iupac' in mol.__dict__  # cached in instance dict
    assert mol.iupac is first  # same object, not recomputed
