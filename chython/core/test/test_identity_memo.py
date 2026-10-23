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
"""`identity_memo()`: identities reused across byte-equal arenas, answers unchanged."""
from chython.core import identity_memo, read_smiles
from chython.core import _core


DRAWINGS = (
    'C[C@H](O)CC',
    'C[C@@H](O)CC |&1:1|',
    'C[C@H](O)[C@@H](C)O |o1:1,3|',
    'O[C@H]1C[C@@H](O)C1',                               # cis-cyclobutane-1,3-diol, a mirror-symmetric pair
    'NC(=O)C1[C@H](C(N)=O)[C@H]1C(N)=O |&1:4,8|',        # cyclopropane-1,2,3-tricarboxamide, trans pair
    'C[C@@]12CC[C@@H](C1(C)C)CC2=O',                     # camphor
)


def _answers():
    out = []
    for text in DRAWINGS:
        m = read_smiles(text)
        c = m.copy()
        with c.edit() as e:
            e.set_parity(next(u['anchor'] for u in m.stereo_units() if u['parity']), 0)
        out.append((m.canonical_bytes, hash(m), str(m), m == c, str(c)))
    return out


def test_the_memo_changes_no_answer():
    plain = _answers()
    with identity_memo():
        assert _answers() == plain
        assert _answers() == plain                      # the second pass is answered from the memo


def test_the_memo_is_scoped_and_nests():
    assert _core._identity_memo_size() is None
    with identity_memo():
        read_smiles('C[C@H](O)CC').canonical_bytes
        filled = _core._identity_memo_size()
        assert filled
        with identity_memo():
            assert _core._identity_memo_size() == filled    # an inner scope shares the outer memo
        assert _core._identity_memo_size() == filled
    assert _core._identity_memo_size() is None              # leaving the outermost drops it


def test_the_memo_is_dropped_when_the_scope_raises():
    try:
        with identity_memo():
            raise KeyError
    except KeyError:
        pass
    assert _core._identity_memo_size() is None
