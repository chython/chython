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
"""The resonance endpoint table: every row compiles, and every row still matches its own probe.

A pattern that matches nothing is invisible -- `fix_resonance` does not raise, it just stops
recognizing that endpoint.  The dialect makes it easy to write by accident: an unstated charge span
means neutral, not "any", so `[C,N,O;-]` written as `[C,N,O]` excludes every anion.  Hence the probe
column and `test_probe_matches`; the rest is structure.
"""
import pytest

from .. import RESONANCE_ROLES, resonance_rules, resonance_rules_by_role, resonance_table_text
from ...core import read_smiles


ROWS = resonance_rules()
IDS = [row.id for row in ROWS]


def test_table_is_not_empty():
    assert len(ROWS) > 10, 'the table lost rows; every one of them is a chython 2 opt-out'


def test_ids_unique():
    assert len(set(IDS)) == len(IDS)


def test_ids_are_namespaced():
    # a log record carries this id and nothing else that says where it came from
    assert all(row.id.startswith('resonance:') for row in ROWS)


def test_roles_in_vocabulary():
    assert {row.role for row in ROWS} <= set(RESONANCE_ROLES)


def test_every_role_used():
    """A role with no rows is either a dead branch in the pass or a table that lost its rows."""
    by_role = resonance_rules_by_role()
    assert set(by_role) == set(RESONANCE_ROLES)
    for role in RESONANCE_ROLES:
        assert by_role[role], f'no row plays the role {role}'


def test_grouping_is_a_partition():
    by_role = resonance_rules_by_role()
    assert sum(len(rows) for rows in by_role.values()) == len(ROWS)


@pytest.mark.parametrize('row', ROWS, ids=IDS)
def test_one_anchor(row):
    """The endpoint is the `:1` atom, and the loader is what enforces it; this pins the intent."""
    assert sum(n == 1 for n in row.query.map_numbers().values()) == 1
    assert row.anchor in row.query.query_numbers()


@pytest.mark.parametrize('row', ROWS, ids=IDS)
def test_probe_matches(row):
    """The row's own probe must match it, at the anchor.  The one test that catches a dead pattern."""
    molecule = read_smiles(row.probe)
    assert row.query.may_match(molecule), \
        f'{row.id}: the cheap screen already rejects its own probe {row.probe!r}'
    mappings = list(row.query.get_mapping(molecule))
    assert mappings, f'{row.id}: {row.smarts!r} matches nothing in its own probe {row.probe!r}'
    assert all(row.anchor in mapping for mapping in mappings)


@pytest.mark.parametrize('row', ROWS, ids=IDS)
def test_row_is_documented(row):
    """Every row carries a non-trivial comment saying what it claims and where it came from."""
    assert len(row.comment) > 20


def test_table_text_is_readable():
    """`resonance_table_text` is the `joinpath` site `chython/test/test_packaging.py` keys on."""
    text = resonance_table_text()
    assert text.startswith('#')
    header = next(line for line in text.splitlines() if not line.startswith('#'))
    assert header.split('\t') == ['id', 'role', 'smarts', 'probe', 'comment']
