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
"""Wildman-Crippen logP and molar refractivity over `tables/crippen.tsv`.

Wildman, Crippen, J. Chem. Inf. Comput. Sci. 1999, 39, 868.
"""
from ._standardize import LogRecord
from ._tables import crippen_rules_by_role, first_match
from ..core import recording


def crippen_contributions(molecule) -> dict:
    """Per-heavy-atom Wildman-Crippen contribution: `{n: (type, logp, mr)}`.

    Wildman, Crippen, J. Chem. Inf. Comput. Sci. 1999, 39, 868.  The value already includes the
    atom's hydrogens: an H row types the carrier and is multiplied by its total hydrogen count, so the
    answer does not depend on whether the caller made hydrogens explicit.  An explicit hydrogen atom
    gets its own zero-contribution entry, its contribution already being counted through the carrier.

    Reads only, never edits.  The rows must therefore constrain hydrogens with `H`, not `D`, or a type
    would depend on whether hydrogens are explicit -- see `crippen.tsv`'s header.
    """
    all_atom_ids = [a.n for a in molecule.atoms()]
    explicit_h_ids = frozenset(a.n for a in molecule.atoms() if a.element == 1)

    by_role = crippen_rules_by_role()
    heavy = first_match(by_role['heavy'], molecule)
    carriers = first_match(by_role['hydrogen'], molecule)

    out = {}
    lines = []
    for i in all_atom_ids:
        if i in explicit_h_ids:
            # not an untyped atom: already counted through its carrier's H row, so zero and no log line
            out[i] = ('-', 0.0, 0.0)
            continue
        row = heavy.get(i)
        if row is None:
            # no block catch-all matched: an element outside C/N/O/H/F/Cl/Br/I/P/S/metal
            lines.append(LogRecord('crippen:untyped', (i,),
                                   'no Wildman-Crippen type matches this atom; it contributes '
                                   'zero to logP and MR'))
            out[i] = ('-', 0.0, 0.0)
            continue
        logp, mr = row.logp, row.mr
        h_row = carriers.get(i)
        if h_row is not None:
            n = molecule.total_h_of(i)
            if n is None:
                lines.append(LogRecord('crippen:h-unknown', (i,),
                                       'the hydrogen count is unknown, so the hydrogen contribution '
                                       'is omitted; run kekule() then calc_implicit'))
            elif n:
                logp += h_row.logp * n
                mr += h_row.mr * n
        out[i] = (row.type, logp, mr)
    with recording(molecule, stage='crippen') as log:
        log.extend(lines)
    return out


def crippen_logp(molecule) -> float:
    """Wildman-Crippen atomic-contribution logP.

    Wildman, Crippen, J. Chem. Inf. Comput. Sci. 1999, 39, 868.
    """
    return sum(p[1] for p in crippen_contributions(molecule).values())


def crippen_mr(molecule) -> float:
    """Wildman-Crippen molar refractivity.

    Wildman, Crippen, J. Chem. Inf. Comput. Sci. 1999, 39, 868.  Four published types (N10, N12,
    Hal, Me2) have no MR value; they contribute zero and `tables/crippen.tsv` flags them, so an
    absent value never masquerades as a measured zero.
    """
    return sum(p[2] for p in crippen_contributions(molecule).values())
