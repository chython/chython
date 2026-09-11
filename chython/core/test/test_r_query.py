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
"""A molecule carrying an R marker is refused as a query, and is a legitimate target.

`as_query()` is the one door: the comparison operators and `is_substructure` all route through it, so
each of the five spellings below refuses for the same reason.
"""
from pytest import raises
from chython.core import read_smiles


def _pair():
    """An R-bearing fragment and an ordinary molecule to compare it against."""
    return read_smiles('[R1]c1ccccc1'), read_smiles('Cc1ccccc1')


def test_the_le_operator_refuses_an_r():
    t, tol = _pair()
    with raises(ValueError, match='matches nothing'):
        t <= tol


def test_the_ge_operator_refuses_an_r():
    t, tol = _pair()
    with raises(ValueError, match='matches nothing'):
        tol >= t


def test_is_substructure_refuses_an_r():
    t, tol = _pair()
    with raises(ValueError, match='matches nothing'):
        t.is_substructure(tol)


def test_the_membership_test_refuses_an_r():
    t, tol = _pair()
    with raises(ValueError, match='matches nothing'):
        t in tol


def test_as_query_refuses_an_r():
    t, _ = _pair()
    with raises(ValueError, match='matches nothing'):
        t.as_query()


def test_an_r_bearing_molecule_is_a_legitimate_target():
    """The reverse direction is a question with an answer: `atom_admits` says no R is a carbon."""
    t, tol = _pair()
    assert (tol in t) is False


# The same refusal on a fragment the stickers enumerator built lives in
# `chython/reactions/test/test_stickers.py`: `sticky_fragments` is injected by `chython.reactions`, and
# `test_no_chython_two_imports.py` forbids this directory any import of this distribution but
# `chython.core`.
