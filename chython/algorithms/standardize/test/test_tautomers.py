# -*- coding: utf-8 -*-
#
#  Copyright 2025 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from pytest import mark


data = [
    ('CC1=NNC=C1', 'Cc1cc[nH]n1'),
    ('CC1=CC=NN1', 'Cc1cc[nH]n1'),
    ('CC1=NC=NN1', 'Cc1[nH]cnn1'),
    ('CC1=NN=CN1', 'Cc1[nH]cnn1'),
    ('CC1=NNN=N1', 'n1nn[nH]c1C'),
    ('CC1=NN=NN1', 'n1nn[nH]c1C'),
    ('CC1=CN=CN1', 'Cc1[nH]cnc1'),
    ('CC1=CNC=N1', 'Cc1[nH]cnc1'),
    ('CC1=CN=NN1', 'n1n[nH]cc1C'),
    ('CC1=CNN=N1', 'n1n[nH]cc1C'),
    ('CC1=NNN=C1', 'n1n[nH]cc1C'),
    ('CC1=NC2=NNC(C)=C2N1', 'Cc1n[nH]c2[nH]c(C)nc12'),
    ('CC1=NC2=C(N1)C(C)=NN2', 'Cc1n[nH]c2[nH]c(C)nc12'),
    ('CC1=NC2=C(C)NN=C2N1', 'Cc1n[nH]c2[nH]c(C)nc12'),
    ('COC(=N)NC', 'COC(N)=NC'),
    ('COC(N)=NC', 'COC(N)=NC'),
    ('CCN=C(N)NC', 'CCNC(=NC)N'),
    ('CCNC(=N)NC', 'CCNC(=NC)N'),
    ('CCNC(N)=NC', 'CCNC(=NC)N'),
    ('CNC(N)=NC(=N)NC', 'CNC(=N)NC(=N)NC'),
    ('CNC(=N)NC(=N)NC', 'CNC(=N)NC(=N)NC'),
    ('CNC(N)=NC(N)=NC', 'CNC(=N)NC(=N)NC'),
    ('CCN=CNC=NC', 'CCN=CN=CNC')
]


@mark.parametrize('raw,result', data)
def test_group(raw, result):
    tmp = smiles(raw)
    tmp.kekule()
    tmp.thiele()
    tmp.standardize_tautomers()
    assert tmp == smiles(result), f'{raw} > {tmp} != {result}'
