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
"""Connectivity from a stored model: which atoms a geometry says are bonded.

The other half of `saturate()`.  A file stating coordinates and no bond -- an XYZ frame, a QM output --
gives a set of atoms and their positions, and the two questions that turns into a molecule are asked
by two explicit calls: `perceive_bonds()` says WHICH PAIRS are bonded, `saturate()` says at WHAT
ORDER.  Neither runs on read.

One threshold, `(r(n) + r(m)) * radius_multiplier` over `tables/covalent_radii.tsv`, and nothing else:
no element special case, no ring closure, no hydrogen rule.  A pair inside it is bonded and a pair
outside it is not, which is what makes the answer a function of the geometry the file stated.
"""
from ._tables import covalent_radii
from ..core import INFO, LogRecord, LOST, MoleculeContainer, recording


__all__ = ['perceive_bonds']


#: One id per outcome a caller filters on; see the table in `perceive_bonds`'s docstring.
_RULE = 'perceive:bonds'
_RULE_STATED = 'perceive:stated'
_RULE_NO_RADIUS = 'perceive:no-radius'


#: The 27 cells a pair can span, the cell being as wide as the longest bond the radii admit.
_NEIGHBOUR_CELLS = tuple((i, j, k) for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1))


def _pairs_within_reach(known: list, xyz: dict, reach: dict):
    """Every `(n, m)` pair, `n` before `m` in `known`, whose distance is inside its own threshold.

    A cell hash and not a pair loop, because the pair loop is quadratic and a 3400-atom model spends
    a second in it: bin the atoms on a grid as wide as the longest bond any pair of radii admits, and
    a bonded pair can then only be in the cell itself or one of its 26 neighbours.  The answer is the
    pair loop's exactly -- the same threshold decides every pair the grid brings together.
    """
    cell = 2. * max(reach.values())
    buckets: dict[tuple[int, int, int], list[int]] = {}
    for n in known:
        x, y, z = xyz[n]
        buckets.setdefault((int(x // cell), int(y // cell), int(z // cell)), []).append(n)

    position = {n: i for i, n in enumerate(known)}
    for (cx, cy, cz), members in buckets.items():
        neighbours = [m for i, j, k in _NEIGHBOUR_CELLS
                      for m in buckets.get((cx + i, cy + j, cz + k), ())]
        for n in members:
            nx, ny, nz = xyz[n]
            for m in neighbours:
                if position[m] <= position[n]:
                    continue          # each unordered pair is offered twice; take it once
                mx, my, mz = xyz[m]
                limit = reach[n] + reach[m]
                if (mx - nx) ** 2 + (my - ny) ** 2 + (mz - nz) ** 2 <= limit * limit:
                    yield n, m


def perceive_bonds(molecule: MoleculeContainer, *, model: int = 0,
                   radius_multiplier: float = 1.25) -> bool:
    """Add a single bond for every atom pair whose distance in `model` says they are bonded.

    Returns `True` when a bond was added.  Every bond is order 1, the order no valence rule has to
    justify; `saturate()` raises the ones a hydrogen count forces.  A bond the molecule already holds
    is left exactly as it is, order included -- perception adds connectivity and never revises it.

    ==========================  ==========  =====================================================
    rule                        severity    what it says
    ==========================  ==========  =====================================================
    ``perceive:bonds``          info        this many bonds were added, from which model
    ``perceive:stated``         info        this many pairs within reach were bonded already
    ``perceive:no-radius``      lost        these atoms state an element with no radius, so
                                            nothing is bonded to them
    ==========================  ==========  =====================================================

    `radius_multiplier` scales the sum of the two covalent radii.  The default answers every measured
    bond length and every measured nonbonded contact in `test_covalent_radii_tsv.py`, whose two corpora
    leave the window 1.2386 (fluorine's F-F bond, the longest bond relative to its radii) to 1.2961
    (cyclobutadiene's transannular carbons, the tightest contact) -- a 4.6% window, so the knob is not a
    free parameter.  A caller drawing a metal cluster or reading a stretched transition state moves it
    and reads the log.

    ONE THRESHOLD CANNOT ANSWER EVERY GEOMETRY, and the case that proves it is in that test file:
    bicyclo[1.1.1]pentane's bridgehead carbons are 1.845 A apart and not bonded, which no multiplier
    rejects while still reaching F2's bond.  A single distance rule is what this pass is; where a
    structure is strained enough for the two to overlap, the log is what a caller reads.

    A MODEL IS READ, NEVER GUESSED.  `model` indexes the conformer store, so a molecule carrying no
    geometry raises `IndexError` from the container rather than being handed an invented one.
    """
    conformer = molecule.conformer(model)          # IndexError names the model that is not there

    # Read everything first: the container answers no query once an edit session is open.  Sorted,
    # so the bonds are added in stable id order whatever the arena's own order is.
    radii = covalent_radii()
    numbers = sorted(molecule.atom_numbers)
    xyz = {n: conformer.xyz_of(n) for n in numbers}
    bonded = {n: set(molecule.neighbors_of(n)) for n in numbers}
    reach = {}
    unknown = []
    for n in numbers:
        element = molecule.element_of(n)
        if element in radii:
            reach[n] = radii[element] * radius_multiplier
        else:
            # An element with no row -- the R marker, or one past the survey the table covers.  A
            # radius is not invented for it, so it takes no bond and says so below.
            unknown.append(n)

    known = [n for n in numbers if n in reach]
    add = []
    stated = 0
    for n, m in _pairs_within_reach(known, xyz, reach):
        if m in bonded[n]:
            stated += 1
        else:
            add.append((n, m))
    add.sort()

    if add:
        # ALL-OR-NOTHING, and the edit scope is what enforces it: one refused bond discards the
        # journal, so the molecule is either the perceived one or the one that came in.
        with molecule.edit():
            for n, m in add:
                molecule.add_bond(n, m, 1)

    with recording(molecule, stage='perceive_bonds') as log:
        touched = tuple(sorted({n for pair in add for n in pair}))
        if add:
            log.append(LogRecord(_RULE, touched,
                                 f'{len(add)} bond(s) perceived from model {model} at '
                                 f'{radius_multiplier}x the covalent radii', INFO))
        else:
            log.append(LogRecord(_RULE, (),
                                 f'no bond perceived: no unbonded pair in model {model} lies within '
                                 f'{radius_multiplier}x the sum of its covalent radii', INFO))
        if stated:
            log.append(LogRecord(_RULE_STATED, (),
                                 f'{stated} pair(s) within reach are bonded already and are left as '
                                 'they are, order included', INFO))
        if unknown:
            log.append(LogRecord(_RULE_NO_RADIUS, tuple(unknown),
                                 f'{len(unknown)} atom(s) state an element the covalent radius table '
                                 f'has no row for, so nothing is bonded to them: '
                                 + ', '.join(f'{molecule.atom(n).atomic_symbol} at {n}'
                                             for n in unknown), LOST))
    return bool(add)
