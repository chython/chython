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
from chython import smiles, smarts, Reactor
from pytest import mark


data = [
    (('[B;D3;x2;z1:4]([O:5])([O:6])-[C;@@;h1:3]1([O;M][C;M]1)', '[Cl,Br,I;D1:1]-[C;a:2]'), ('[A;@:3]-[A:2]',),
     ('CC1O[C@@H]1B(O)O', 'Brc1ccccc1'), ('CC1O[C@H]1c1ccccc1',)),  # inverse stereo check
    (('[B;D3;x2;z1:4]([O:5])([O:6])-[C;@@;h1:3]1([O;M][C;M]1)', '[Cl,Br,I;D1:1]-[C;a:2]'), ('[A;@@:3]-[A:2]',),
     ('CC1O[C@@H]1B(O)O', 'Brc1ccccc1'), ('CC1O[C@@H]1c1ccccc1',)),  # keep stereo on RC
    (('[B;D3;x2;z1:4]([O:5])([O:6])-[C;@@;h1:3]1([O;M][C;M]1)', '[Cl,Br,I;D1:1]-[C;a:2]'), ('[A:3]-[A:2]',),
     ('CC1O[C@@H]1B(O)O', 'Brc1ccccc1'), ('CC1OC1c1ccccc1',)),  # drop stereo on RC
]


@mark.parametrize('patterns, products, source, result', data)
def test_transformer(patterns, products, source, result):
    for q, m in zip(patterns, source):
        assert smarts(q) <= smiles(m)

    reactor = Reactor([smarts(x) for x in patterns], [smarts(x) for x in products])
    out = {format(smiles(x), 'h') for x in result}
    assert {format(x, 'h') for x in next(reactor(*(smiles(x) for x in source))).products} == out
