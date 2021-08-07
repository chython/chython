# -*- coding: utf-8 -*-
#
#  Copyright 2021 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from chython import SMILESRead, smiles, MoleculeContainer, CGRContainer
from chython.files._mdl import parse_error
from itertools import zip_longest
from io import StringIO


def test_bad_smiles_detection():
    data = ('[[C] C%[C] [C]] [] [(C)] C%(C) C(C%) C() C(C)) C((C) C%01 C01 S%% CBl CCr Cd C[H C(C [Zn=>=] [12] [+]' 
            '[C:1-] [C+5>0] [1+>-] [=*] [C+>-5] (=C) 1CC1 C#(C) C(C#) C=-C C.1C1 C1CC1C1CC1 C=1CC1 C#1CC-1 C1CC-1'
            'C(C)C( C1C CC= [->=]').replace(' ', '\n')

    with StringIO(data) as f, SMILESRead(f, ignore=False, store_log=True) as r:
        for s in r._data:
            assert isinstance(s, parse_error)


def test_invalid_ring_closure_ignoring():
    data = 'C1CC1C1CC1\nC=1CC1'
    with StringIO(data) as f, SMILESRead(f, ignore=True, store_log=True) as r:
        for t, v in zip_longest(r, ('C1CC1C2CC2', 'C1=CC1')):
            assert t == smiles(v)


def test_good_smiles():
    for t, v in zip('CN C-N C=N C#N C:N cn'.split(), (1, 1, 2, 3, 4, 4)):
        t = smiles(t)
        assert isinstance(t, MoleculeContainer)
        assert t.bond(1, 2).order == v

    t = smiles('C.N ')
    assert isinstance(t, MoleculeContainer)
    assert not t.has_bond(1, 2)

    for t, v1, v2 in zip('C[->.]N C[->=]N C[=>#]N C[#>:]N C[:>.]N [O->0][.>-][Na+>0]'.split(),
                         (1, 1, 2, 3, 4, None), (None, 2, 3, 4, None, 1)):
        t = smiles(t)
        assert isinstance(t, CGRContainer)
        assert t.bond(1, 2).order == v1
        assert t.bond(1, 2).p_order == v2

    for t in '[Fe] [13C] [C+4] [C:14] [C--:12] [C@@H] [C@H] [OH2] [OH3+]'.split():
        assert isinstance(smiles(t), MoleculeContainer)

    for t in 'C(C)C C(-C)C C(C)-C'.split():
        assert smiles(t).bonds_count == 2

    for t, v in zip('C1CC1 C-1CC-1 C%10CC%10 C-%11CC-%11 C12CC1C2 C-1-2CC-1C-2'.split(), (1, 1, 1, 1, 2, 2)):
        assert smiles(t).rings_count == v

    for t in ('[C+>-] [C++>--] [C-3>+2] [C0>+] [C-->0] [C*>^] [C^>*] [C+>-*>^] [C-3>+2^>*] '
              '[C0>-*>^] [C+2>0^>*] [n+>0]').split():
        assert isinstance(smiles(t), CGRContainer)
