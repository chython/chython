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
"""Rotatable bond and H-bond donor/acceptor counts, over `tables/rotatable.tsv` and `tables/hbond.tsv`.
"""
from ._standardize import LogRecord
from ..core import recording
from ._tables import hbond_rules_by_role, rotatable_rules_by_role


def _mapped_bond(row, mapping):
    """The (low, high) stable-id pair the row's :1 and :2 matched, order-normalised."""
    a = mapping[row.numbers[1]]
    b = mapping[row.numbers[2]]
    return (a, b) if a < b else (b, a)


def rotatable_bonds_count(molecule) -> int:
    """Number of rotatable bonds, by the definition in `tables/rotatable.tsv`.

    A bond is counted once however many ways a pattern maps onto it: row 1 is symmetric in its two
    atoms and the matcher offers each bond in both directions, so the set is what makes this a count of
    bonds rather than of matches.  Charged and radical atoms count.
    """
    by_role = rotatable_rules_by_role()
    found = set()
    excluded = []
    for row in by_role['rotatable']:
        for mapping in row.query.get_mapping(molecule):
            found.add(_mapped_bond(row, mapping))
    for row in by_role['exclude']:
        for mapping in row.query.get_mapping(molecule):
            bond = _mapped_bond(row, mapping)
            if bond in found:
                found.discard(bond)
                excluded.append(LogRecord(row.id, bond, row.description))
    with recording(molecule, stage='rotatable') as log:
        log.extend(excluded)
    return len(found)


def hbond_atoms(molecule, role: str) -> frozenset:
    """Stable ids of the atoms `tables/hbond.tsv` types with `role`.

    Shared by the two counts and by `pharmacophore_invariants`, which is why it returns ids rather
    than a number.
    """
    rules = hbond_rules_by_role()
    if role not in rules:
        raise ValueError(f'role {role!r} is not one of {tuple(rules)}')
    out = set()
    for row in rules[role]:
        subject = row.numbers[1]          # the stable id of :1, from compile_smarts at load time
        for mapping in row.query.get_mapping(molecule):
            out.add(mapping[subject])
    return frozenset(out)


def hydrogen_bond_donors_count(molecule) -> int:
    """Count of hydrogen bond donor ATOMS, over `tables/hbond.tsv`."""
    return len(hbond_atoms(molecule, 'donor'))


def hydrogen_bond_acceptors_count(molecule) -> int:
    """Count of hydrogen bond acceptor ATOMS, over `tables/hbond.tsv`."""
    return len(hbond_atoms(molecule, 'acceptor'))
