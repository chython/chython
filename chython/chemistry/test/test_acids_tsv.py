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
"""The acid/base table: every row compiles, and every row still matches its own probe.

A pattern that matches nothing is invisible -- `neutralize` does not raise, it just stops recognizing
that site.  The dialect makes that easy to write by accident: an unstated charge span means neutral,
not "any", so `[O;D1]` written for a carboxylate matches no anion at all.  Hence the probe column and
`test_probe_matches`; the rest is structure, plus the one claim the table itself makes -- every row
names a CHARGED site, an `acid` holding an implicit hydrogen and a `base` able to take one.
"""
import pytest

from .. import ACID_ROLES, acids_rules, acids_rules_by_role, acids_table_text
from ...core import read_smiles


ROWS = acids_rules()
IDS = [row.id for row in ROWS]


def test_table_is_not_empty():
    assert len(ROWS) > 5, 'the table lost rows; each one is a site `neutralize` can no longer see'


def test_ids_unique():
    assert len(set(IDS)) == len(IDS)


def test_ids_are_namespaced():
    # a log record carries this id and nothing else that says where it came from
    assert all(row.id.startswith('acids:') for row in ROWS)


def test_roles_in_vocabulary():
    assert {row.role for row in ROWS} <= set(ACID_ROLES)


def test_every_role_used():
    """A role with no rows makes one half of the pass dead: a proton needs both ends."""
    by_role = acids_rules_by_role()
    assert set(by_role) == set(ACID_ROLES)
    for role in ACID_ROLES:
        assert by_role[role], f'no row plays the role {role}'


def test_grouping_is_a_partition():
    by_role = acids_rules_by_role()
    assert sum(len(rows) for rows in by_role.values()) == len(ROWS)


@pytest.mark.parametrize('row', ROWS, ids=IDS)
def test_one_anchor(row):
    """The site is the `:1` atom, and the loader is what enforces it; this pins the intent."""
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
def test_site_is_charged(row):
    """A neutral acid or base belongs to `enumerate_charged_forms`, not here.

    An `acid` row must land on a cation and a `base` row on an anion in its own probe, which is the
    table's header claim and the reason `neutralize` never creates charge.
    """
    molecule = read_smiles(row.probe)
    charges = {molecule.charge_of(mapping[row.anchor]) for mapping in row.query.get_mapping(molecule)}
    if row.role == 'acid':
        assert all(charge > 0 for charge in charges), f'{row.id}: an acid site must be a cation'
    else:
        assert all(charge < 0 for charge in charges), f'{row.id}: a base site must be an anion'


@pytest.mark.parametrize('row', [row for row in ROWS if row.role == 'acid'],
                         ids=[row.id for row in ROWS if row.role == 'acid'])
def test_acid_site_has_a_proton(row):
    """`h` reads IMPLICIT hydrogens, so an acid row that forgot it would match an aprotic cation."""
    molecule = read_smiles(row.probe)
    for mapping in row.query.get_mapping(molecule):
        assert molecule.implicit_h_of(mapping[row.anchor]), \
            f'{row.id}: matched an atom with no implicit hydrogen to give away'


@pytest.mark.parametrize('row', ROWS, ids=IDS)
def test_row_is_documented(row):
    """Every row carries a non-trivial comment saying what it claims."""
    assert len(row.comment) > 20


def test_table_text_is_readable():
    """`acids_table_text` is the `joinpath` site `chython/test/test_packaging.py` keys on."""
    text = acids_table_text()
    assert text.startswith('#')
    header = next(line for line in text.splitlines() if not line.startswith('#'))
    assert header.split('\t') == ['id', 'role', 'smarts', 'probe', 'comment']
