# -*- coding: utf-8 -*-
#
#  Copyright 2019-2026 Ramil Nugmanov <nougmanoff@protonmail.com>
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
"""Arranging a reaction's molecules left to right, as module functions over the core containers.

The molecules' planes are stored; the arrow and the signs are returned and never stored, because they
belong to one drawing at one style rather than to the reaction.
"""
from . import molecule as _molecule


Plane = dict[int, tuple[float, float]]


def layout2d(rxn, *, engine=None, force: bool = False) \
        -> tuple[list[Plane], tuple[float, float, float], list[tuple[float, float]]]:
    """Lay the reaction out and return `(planes, arrow, signs)` without touching anything.

    `planes` holds one `{n: (x, y)}` mapping per molecule, in `molecules()` order; `arrow` is
    `(x1, x2, y)`; `signs` is one `(x, y)` per `+` between two members of a side.  `clean2d()` is this
    plus the decision to keep the planes.  A molecule that already carries a layout keeps it and is only
    shifted into place; one that does not is laid out first by its own `layout2d`.

    :param engine: override the globally set 2d engine
    :param force: lay every molecule out again, whatever coordinates it already has
    """
    planes = []
    for m in rxn.molecules():
        planes.append(_molecule.layout2d(m, engine=engine, force=force))
    arrow, signs = _position(rxn, planes)
    return planes, arrow, signs


def clean2d(rxn, *, engine=None, force: bool = False) \
        -> tuple[tuple[float, float, float], list[tuple[float, float]]]:
    """Lay the reaction out, store each molecule's plane, and return the arrow and the signs.

    This always stores, unlike the molecule's `clean2d`: the arrangement is what the call produces, and
    a member that already had a layout still has to be moved onto the row.  `force` reaches only the
    members' own layouts.

    The arrow is not bit-idempotent across the first two calls, by one quantisation step: storing rounds
    a coordinate onto the arena's 1e-4 grid, so the arrow -- derived from the members' extents -- moves
    once by under 5e-5 and is exact from the second call onward.
    """
    planes, arrow, signs = layout2d(rxn, engine=engine, force=force)
    for m, plane in zip(rxn.molecules(), planes):
        _molecule._store_plane(m, plane)
    return arrow, signs


def _position(rxn, planes: list[Plane]) \
        -> tuple[tuple[float, float, float], list[tuple[float, float]]]:
    """Arrange the planes left to right and return the arrow span and the `+` sign positions.

    Shifts the given planes in place and never reads or writes the molecules.  The arrangement constants
    live here and nowhere else -- `ReactionStyle` sizes the arrow head and the `+` glyph and owns none of
    them.  The `3` is a minimum advance of `shift_x`, so the minimum arrow span is 2:
    `arrow_max = shift_x - 1` and the third unit is the clearance before the first product.  The row is
    the `y = 0` axis, since `_shift_plane_mean` centres every member's box on it.
    """
    # `planes` is in `molecules()` order and is sliced by position, not looked up per molecule: a
    # reaction may hold the same container object twice and each occurrence gets its own plane.
    reactants, agents, products = rxn.reactants, rxn.agents, rxn.products
    split = len(reactants) + len(agents)
    r_planes, g_planes, p_planes = planes[:len(reactants)], planes[len(reactants):split], planes[split:]

    shift_x = 0
    amount = len(reactants) - 1
    signs = []
    for m, plane in zip(reactants, r_planes):
        max_x = _molecule._shift_plane_mean(m, plane, shift_x)
        if amount:
            max_x += .2
            signs.append((max_x, 0.))
            amount -= 1
        shift_x = max_x + 1
    arrow_min = shift_x

    if agents:
        shift_x += .4
        for m, plane in zip(agents, g_planes):
            max_x = _molecule._shift_plane_min(m, plane, shift_x, .5)
            shift_x = max_x + 1
        shift_x += .4
        if shift_x - arrow_min < 3:
            shift_x = arrow_min + 3
    else:
        shift_x += 3
    arrow_max = shift_x - 1

    amount = len(products) - 1
    for m, plane in zip(products, p_planes):
        max_x = _molecule._shift_plane_mean(m, plane, shift_x)
        if amount:
            max_x += .2
            signs.append((max_x, 0.))
            amount -= 1
        shift_x = max_x + 1
    return (arrow_min, arrow_max, 0.), signs


__all__ = ['clean2d', 'layout2d']
