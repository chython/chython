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
from chython import smiles, smarts, Transformer
from pytest import mark


data = [
    ('[C:1]Br', '[A:1][O;M]', 'C[C@H](OC)CBr', 'C[C@H](OC)CO'),  # keep stereo out of match
    ('[C:2][C:1]Br', '[A:2][A:1][O;M]', 'C[C@H](OC)CBr', 'C[C@H](OC)CO'),  # keep stereo inside match
    ('[C;M][C;@;h1:1]([O;M])[N;M]', '[A;@@:1]', 'CC[C@H](O)N', 'CC[C@@H](O)N'),  # inversion of stereo
    ('[C:1]Br', '[A:1][O;M]', 'C/C=C/CBr', 'C/C=C/CO'),  # keep stereo out of match
    ('[C:1]Br', '[A:1][O;M]', 'CC=[C@]=CCBr', 'CC=[C@]=CCO'),  # keep
    ('[C:1]Br', '[A:1][O;M]', 'CC=[C@]=CBr', 'CC=C=CO'),  # drop stereo on RC
    ('[C:1]Br', '[A:1][O;M]', 'C/C=C/Br', 'CC=CO'),  # drop stereo on RC
]


@mark.parametrize('pattern, replacement, source, result', data)
def test_transformer(pattern, replacement, source, result):
    transformer = Transformer(smarts(pattern), smarts(replacement))

    mol = smiles(source)
    out = {format(smiles(x), 'h') for x in ([result] if isinstance(result, str) else result)}
    assert {format(x, 'h') for x in transformer(mol)} == out
