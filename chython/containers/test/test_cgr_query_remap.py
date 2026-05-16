# -*- coding: utf-8 -*-
#
#  Copyright 2026 Tagir Akhmetshin <tagirshin@gmail.com>
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
"""Regression tests for :class:`QueryCGRContainer.remap`.

The override on ``QueryCGRContainer.remap`` forwarded a ``copy=`` kwarg to
``Graph.remap``, which does not accept it; every call (in-place or copying)
crashed with ``TypeError``.
"""
from chython.containers.cgr_query import QueryCGRContainer
from chython.periodictable import DynamicQueryElement


def _two_atom_qcgr() -> QueryCGRContainer:
    q = QueryCGRContainer()
    q.add_atom(DynamicQueryElement.from_symbol("C"), 1)
    q.add_atom(DynamicQueryElement.from_symbol("O"), 2)
    q.add_bond(1, 2, 1)
    return q


def test_remap_in_place_does_not_raise():
    q = _two_atom_qcgr()
    q.remap({1: 11, 2: 12})
    assert set(q._atoms) == {11, 12}
    assert set(q._bonds[11]) == {12}


def test_remap_copy_does_not_raise():
    q = _two_atom_qcgr()
    h = q.remap({1: 21, 2: 22}, copy=True)
    # original untouched
    assert set(q._atoms) == {1, 2}
    assert set(h._atoms) == {21, 22}
    assert set(h._bonds[21]) == {22}


def test_remap_preserves_side_state_dicts():
    # In-place remap must shuttle the per-atom side dicts to the new keys.
    q = _two_atom_qcgr()
    q._p_charges[1] = 1  # mark something distinctive
    q.remap({1: 99, 2: 98})
    assert 99 in q._p_charges and 1 not in q._p_charges
    assert q._p_charges[99] == 1
