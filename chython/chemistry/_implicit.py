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
"""Derive implicit hydrogen counts from the valence collection, and report what it rejects.

Owns no table: the valence rows live in `chython/core/valence_rules.tsv` and are compiled into the
extension.  Both functions delegate to `chython.core`, the one derivation every reader shares, so
`calc_implicit` and `check_valence` cannot answer differently about one atom.
"""
from ..core import MoleculeContainer
from ..core._core import valence_report


__all__ = ['calc_implicit', 'check_valence', 'environment_of']


def environment_of(molecule: MoleculeContainer, n: int) -> tuple[int, list[tuple[int, int]], int]:
    """`(order_sum, [(order, element), ...], aromatic_bonds)` for atom `n`.

    Explicit hydrogens are neighbours like any other and are already in `neighbors_of`; the
    molecule's own `explicit_h_of` is not added again, which would double-count them.

    The valence collection has no row for orders 4 and 8, so the policy on both lives here.  Order 8, a
    dative bond, contributes nothing: a donated lone pair is not a σ-bonding slot, and counting it would
    make every metal carbonyl a valence violation.  Order 4 is counted separately and returned as the
    third element, never as an environment entry, so a caller cannot accidentally sum it.
    """
    order_sum = 0
    aromatic = 0
    environment: list[tuple[int, int]] = []
    for m in molecule.neighbors_of(n):
        order = molecule.order_of(n, m)
        if order == 8:                     # dative: outside valence bookkeeping entirely
            continue
        elif order == 4:                   # aromatic: no row exists, do not pretend one does
            aromatic += 1
            continue
        order_sum += order
        environment.append((order, molecule.element_of(m)))
    return order_sum, environment, aromatic


def calc_implicit(molecule: MoleculeContainer, n: int) -> int | None:
    """Recompute and write atom `n`'s implicit hydrogen count.  Returns what was written.

    `None` when nothing local can derive the count; it is then written as `H_UNKNOWN`, never zero --
    an atom whose hydrogens nobody can derive is not an atom with no hydrogens.  That is also why this
    never raises: a metal the collection says nothing about must survive standardization.

    Two reasons for `None`: no valence row for this element in this charge and radical state, or the
    aromatic pnictogen whose count the ring decides (pyrrole versus pyridine).  `kekule()` resolves
    the second.  Dative and aromatic bonds are handled as `environment_of` describes.

    `MoleculeContainer.calc_implicit` is this, and holds the body; the function is kept because
    `check_valence` beside it is a function too, and the pair reads as one module.
    """
    return molecule.calc_implicit(n)


def check_valence(molecule: MoleculeContainer) -> list[tuple[int, str]]:
    """`[(atom, verdict)]` for every atom whose state the collection does not call `'valid'`.

    Two verdicts, and they are not the same claim.  `'violation'` means the collection describes this
    element in this charge and radical state and no row accepts what the molecule has -- a statement
    about the molecule.  `'unknown'` means no complete question could be put: nothing is described
    there, or the atom's aromatic class is the ring's to decide, or its count was never derived.
    Merging them is how a coverage hole gets mistaken for bad input, so they stay apart.

    An aromatic atom is a violation only when neither Kekule reading has a row, since a claim about
    the molecule must survive every form the ring could take.  Never raises and never edits.
    """
    return valence_report(molecule)
