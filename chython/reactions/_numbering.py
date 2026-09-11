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
"""Atom correspondence between two structures, and agreement between two mappings of one reaction.

The compare half of the reconstruction ladder, kept apart from the generate half.
"""

from ..core import MoleculeContainer, ReactionContainer


def fast_mapping(a: MoleculeContainer, b: MoleculeContainer) -> dict[int, int] | None:
    """`{a stable id: b stable id}` when `a` and `b` are the same structure, else `None`.

    Two canonical orderings composed through their shared positions, so it is stereo-aware: meso and
    (R,R) tartaric acid do not correspond.  Among automorphic atoms the choice is arbitrary but
    consistent, which is what `mapping_agrees()` scores modulo.
    """
    if a != b:
        return None
    ao = a.canonical_order()
    bo = b.canonical_order()
    back = {position: n for n, position in bo.items()}
    return {n: back[position] for n, position in ao.items()}


def mapping_agrees(produced: ReactionContainer, reference: ReactionContainer) -> tuple[int, int, int]:
    """`(agreed, disagreed, missing)` over the product atoms of two mappings of one reaction.

    Never routed through container equality: `__eq__` excludes map numbers, so it holds for every
    record and would report 100%.  Compared instead is, per product atom, which input atom it came
    from as `(input index, input stable id)`.  A disagreement is excused when the two input atoms lie
    in one automorphism orbit of the same input.  `missing` counts a product atom one side numbers and
    the other does not -- a template-created atom, or a partially mapped reference.
    """
    got = _pairs(produced)
    want = _pairs(reference)
    orbits = [m.automorphism_orbits() for m in _inputs(reference)]

    agreed = disagreed = missing = 0
    for key, target in want.items():
        source = got.get(key)
        if source is None:
            missing += 1
        elif source == target:
            agreed += 1
        elif source[0] == target[0] and _same_orbit(orbits[source[0]], source[1], target[1]):
            agreed += 1
        else:
            disagreed += 1
    missing += sum(1 for key in got if key not in want)
    return agreed, disagreed, missing


def _inputs(reaction):
    return [*reaction.reactants, *reaction.agents]


def _same_orbit(orbits, n, m):
    a = orbits.get(n)
    return a is not None and a == orbits.get(m)


def _pairs(reaction) -> dict[tuple[int, int], tuple[int, int]]:
    """`{product atom key: (input index, input stable id)}` for every numbered product atom.

    The key is `(product index, canonical position)` and not a stable id: a produced reaction's ids
    came from the patcher, so they are not comparable across the two sides while a position is.
    """
    sources = {}
    for i, molecule in enumerate(_inputs(reaction)):
        for n in molecule.atom_numbers:
            number = molecule.map_number_of(n)
            if number:
                sources[number] = (i, n)

    pairs = {}
    for j, product in enumerate(reaction.products):
        order = product.canonical_order()
        for n in product.atom_numbers:
            number = product.map_number_of(n)
            if not number:
                continue
            source = sources.get(number)
            if source is not None:
                pairs[j, order[n]] = source
    return pairs
