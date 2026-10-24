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
"""Terminal signs over a spiro centre against identity, live group members, and `coupled_implied`: a bridgehead
whose sign the view fills from a signed partner."""
from chython.core import read_smiles


def canon(text):
    """`text` read and canonicalized: signs that state nothing go by the explicit rewrite."""
    m = read_smiles(text)
    m.canonicalize_stereo_groups()
    return m


SPIRO_TERMINALS = 'C[C@H]1CC2(C1)C[C@H](C)C2'       # 2,6-dimethylspiro[3.3]heptane, CID 21304360, centre unsigned
SPIRO_SIGNED = 'C[C@H]1C[C@]2(C1)C[C@H](C)C2'       # CID 21304360, centre signed
SPIRO_UNSIGNED = 'CC1CC2(C1)CC(C)C2'                # CID 21304360


def test_terminal_signs_over_an_unsigned_centre_are_the_unsigned_compound():
    m = canon(SPIRO_TERMINALS)
    assert m == canon(SPIRO_UNSIGNED)
    assert str(m) == str(canon(SPIRO_UNSIGNED))


def test_signing_the_spiro_centre_makes_the_terminals_live():
    m = read_smiles(SPIRO_SIGNED)
    assert all(m.stereo_is_live(a) for a in (2, 4, 7))


def test_mixed_groups_on_unsigned_spiro_terminals_are_the_racemate():
    # flipping one terminal inverts the molecule, so each terminal in AND states both enantiomers
    m = canon(SPIRO_TERMINALS + ' |&1:1,&2:6|')
    assert m != canon(SPIRO_UNSIGNED)
    assert m == read_smiles(str(m))
    assert m == canon(SPIRO_SIGNED + ' |&1:1,3,6|')


def test_every_live_group_member_is_live():
    # 3-methylpentan-2-ol, CID 11261; L-tartaric acid, CID 444305
    for text in ('C[C@H](O)[C@@H](C)CC |&1:1,&2:3|', '[C@@H]([C@H](C(=O)O)O)(C(=O)O)O |&1:0,1|'):
        m = read_smiles(text)
        assert m.live_stereo_groups(), text
        for members in m.live_stereo_groups().values():
            assert all(m.stereo_is_live(x) for x in members), text


def test_coupled_implied_names_the_filled_partner():
    half = read_smiles('C[C@@]12CCC(C1(C)C)CC2=O')          # camphor, CID 2537, C1 signed only
    units = {u['anchor']: u for u in half.stereo_units() if u['coupled']}
    assert units[2]['coupled_implied'] is False and units[5]['coupled_implied'] is True
    full = read_smiles('C[C@@]12CC[C@@H](C1(C)C)CC2=O')
    assert not any(u['coupled_implied'] for u in full.stereo_units())
    assert all(u['coupled_implied'] is False for u in read_smiles('CCO').stereo_units())
