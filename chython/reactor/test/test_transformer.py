# -*- coding: utf-8 -*-
#
#  Copyright 2025, 2026 Ramil Nugmanov <nougmanoff@protonmail.com>
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


# Extended stereo test cases.
#
# Extended stereo encodes epistemic state (what we know about the mixture), not geometry:
#   AND group (+N): racemic mixture — all centers in the group share relative config
#   OR group  (-N): pure compound, unknown absolute config at each center
#
# Rules implemented in _patcher:
#   1. Atoms outside RC: extended stereo always preserved
#   2. Template overrides stereo (swap/inversion): uncertainty is unchanged, preserve reactant's ext stereo
#   3. Stereo lost (center destroyed or retranslation fails): ext stereo dropped (masked by property)
#   4. MoleculeContainer template with extended stereo: template's label takes priority

ext_data = [
    # AND group outside RC preserved
    ('[C:1]Br', '[A:1][O;M]', 'C[C@H](CC)CBr |&1:1|', '[C@H](CO)(CC)C |&1:0|'),
    # OR group outside RC preserved
    ('[C:1]Br', '[A:1][O;M]', 'C[C@H](CC)CBr |o1:1|', '[C@H](CO)(CC)C |o1:0|'),
    # AND group in RC, stereo dropped
    ('[C:1]Br', '[A:1][O;M]', 'CC[C@@H](C)Br |&1:2|', 'CCC(C)O'),
    # AND group in RC, stereo inverted — swap doesn't resolve uncertainty
    ('[C;M][C;@;h1:1]([O;M])[N;M]', '[A;@@:1]', 'CC[C@H](O)N |&1:2|', 'O[C@@H](N)CC |&1:1|'),
    # OR group in RC, stereo inverted — swap doesn't resolve uncertainty
    ('[C;M][C;@;h1:1]([O;M])[N;M]', '[A;@@:1]', 'CC[C@H](O)N |o1:2|', 'O[C@@H](N)CC |o1:1|'),
    # two AND members: one loses stereo in RC, other outside RC keeps label
    ('[C:1]Br', '[A:1][N;M]', 'C[C@H](F)[C@H](O)Br |&1:1,3|', 'OC(N)[C@H](C)F |&1:3|'),
    # multiple independent AND/OR groups preserved
    ('[C:1]Br', '[A:1][O;M]', '[C@H](F)(O)[C@@H](Cl)CBr |&1:0,o1:3|',
     'O[C@@H]([C@@H](Cl)CO)F |&1:1,o1:2|'),
]


@mark.parametrize('pattern, replacement, source, expected', ext_data)
def test_extended_stereo(pattern, replacement, source, expected):
    t = Transformer(smarts(pattern), smarts(replacement))
    p = next(iter(t(smiles(source))))
    assert p == smiles(expected)
