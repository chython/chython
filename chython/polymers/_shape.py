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
"""The shape key: a cheap, stereo-free screen for "could these atoms be that monomer?".

Read off the molecule itself -- no fragment is written or parsed -- as three sorted multisets:

| Part | One entry per |
|---|---|
| atoms | atom: element, charge, implicit hydrogens, aromatic |
| bonds | bond inside the unit: its two end atoms' entries and its order |
| attachments | cut bond leaving the unit: the entry of the atom inside |

A piece's atoms and bonds are computed once; a group of pieces is their union plus the cut bonds between
them, so trying a merge is merging sorted tuples.  Equal shapes are necessary for equal structures, not
sufficient: a match is a candidate, proven or rejected by the exact check on the units a covering
chooses.
"""

from typing import NamedTuple


__all__ = ['Shape', 'Parts', 'shape_of_core']


class Parts(NamedTuple):
    atoms: tuple
    bonds: tuple


class Shape:
    """Atom entries of one molecule, and the shape key of any set of its atoms."""
    __slots__ = ('mol', 'colour')

    def __init__(self, mol):
        self.mol = mol
        self.colour = {n: _colour(mol, n) for n in mol.atom_numbers}

    def piece(self, atoms):
        """The :class:`Parts` of one piece: its atoms and the bonds among them."""
        mol, colour = self.mol, self.colour
        inside = set(atoms)
        bonds = []
        for n in atoms:
            for m in mol.neighbors_of(n):
                if m in inside and n < m:
                    bonds.append(_bond(colour[n], colour[m], mol.order_of(n, m)))
        return Parts(tuple(sorted(colour[n] for n in atoms)), tuple(sorted(bonds)))

    def join(self, parts, internal, boundary):
        """The key of a union of pieces: `internal` are the cut bonds between them, `boundary` the cut
        bonds leaving the union as ``(inside atom, outside atom)``."""
        mol, colour = self.mol, self.colour
        atoms, bonds = [], []
        for p in parts:
            atoms.extend(p.atoms)
            bonds.extend(p.bonds)
        bonds.extend(_bond(colour[n], colour[m], mol.order_of(n, m)) for n, m in internal)
        return tuple(sorted(atoms)), tuple(sorted(bonds)), tuple(sorted(colour[n] for n, _ in boundary))


def shape_of_core(core):
    """The key of a library core: its atoms but the ``*`` markers, and one attachment per marker."""
    colour = {n: _colour(core, n) for n in core.atom_numbers if not core.atom(n).is_r}
    bonds, attachments = [], []
    for n in colour:
        for m in core.neighbors_of(n):
            if m in colour:
                if n < m:
                    bonds.append(_bond(colour[n], colour[m], core.order_of(n, m)))
            else:
                attachments.append(colour[n])
    return tuple(sorted(colour.values())), tuple(sorted(bonds)), tuple(sorted(attachments))


def _colour(mol, n):
    a = mol.atom(n)
    aromatic = any(mol.order_of(n, x) == 4 for x in mol.neighbors_of(n))
    return a.atomic_symbol, a.charge, mol.implicit_h_of(n) or 0, aromatic


def _bond(a, b, order):
    return (a, b, order) if a <= b else (b, a, order)
