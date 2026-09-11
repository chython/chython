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
"""Every row of `functional.tsv` against its own `example` and `decoys`.

WHY THE TABLE CARRIES ITS OWN PROBES.  A group that is absent from a molecule and a group whose pattern
can never match anything are the same missing key in `mol.functional_groups()`, so a typo in a `z` or an
`x` is silent.  Two hundred rows of hand-written SMARTS is exactly the situation where silence is
expensive.
"""
import pytest
from .._tables import functional_rules, read_table
from ...core import read_smiles


GROUPS = tuple(functional_rules().values())


@pytest.mark.parametrize('group', GROUPS, ids=lambda group: group.name)
def test_every_row_matches_its_own_example(group):
    """The one gate a new row cannot pass by accident."""
    molecule = read_smiles(group.example)
    assert next(group.query.get_mapping(molecule), None), (
        f'{group.id} ({group.name}) does not match its own example {group.example!r}: {group.smarts}')


@pytest.mark.parametrize('group', [g for g in GROUPS if g.decoys], ids=lambda group: group.name)
def test_no_row_matches_its_decoys(group):
    """A decoy is the near miss this pattern has to keep rejecting -- the neighbouring group it would
    collapse into if a constraint were dropped, recorded next to the constraint that separates them."""
    for decoy in group.decoys:
        molecule = read_smiles(decoy)
        assert not next(group.query.get_mapping(molecule), None), (
            f'{group.id} ({group.name}) matches its own decoy {decoy!r}: {group.smarts}')


def test_every_row_has_an_example():
    """The column is not optional, which is what makes the gate above a ratchet rather than a sample."""
    missing = [row['name'] for row in read_table('functional.tsv') if not row['example']]
    assert not missing, f'functional.tsv rows with no example: {missing}'
