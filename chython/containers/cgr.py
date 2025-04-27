# -*- coding: utf-8 -*-
#
#  Copyright 2017-2025 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from functools import cached_property
from typing import Dict, Iterator, Tuple, Optional, Collection
from .bonds import DynamicBond
from ..algorithms.fingerprints import FingerprintsCGR
from ..algorithms.isomorphism import Isomorphism
from ..algorithms.morgan import Morgan
from ..algorithms.rings import Rings
from ..algorithms.smiles import CGRSmiles
from ..periodictable import DynamicElement


class CGRContainer(CGRSmiles, Morgan, Rings, Isomorphism, FingerprintsCGR):
    __slots__ = ('_atoms', '_bonds', '__dict__')
    _atoms: Dict[int, DynamicElement]
    _bonds: Dict[int, Dict[int, DynamicBond]]

    def __init__(self):
        self._atoms = {}
        self._bonds = {}

    def atoms(self) -> Iterator[Tuple[int, DynamicElement]]:
        return iter(self._atoms.items())

    def bonds(self) -> Iterator[Tuple[int, int, DynamicBond]]:
        """
        Iterate other all bonds
        """
        seen = set()
        for n, m_bond in self._bonds.items():
            seen.add(n)
            for m, bond in m_bond.items():
                if m not in seen:
                    yield n, m, bond

    @cached_property
    def center_atoms(self) -> Tuple[int, ...]:
        """ Get list of atoms of reaction center (atoms with dynamic: bonds, charges, radicals).
        """
        center = {n for n, a in self._atoms.items() if a.is_dynamic}
        center.update(n for n, m_bond in self._bonds.items() if any(bond.is_dynamic for bond in m_bond.values()))
        return tuple(center)

    def substructure(self, atoms) -> 'CGRContainer':
        """
        Create substructure containing atoms from atoms list

        :param atoms: list of atoms numbers of substructure
        """
        atoms = set(atoms)
        sa = self._atoms
        sb = self._bonds

        sub = object.__new__(self.__class__)
        sub._atoms = {n: sa[n].copy() for n in atoms}

        sub._bonds = cb = {}
        for n in atoms:
            cb[n] = cbn = {}
            for m, bond in sb[n].items():
                if m in cb:  # bond partially exists. need back-connection.
                    cbn[m] = cb[m][n]
                elif m in atoms:
                    cbn[m] = bond.copy()
        return sub

    def augmented_substructure(self, atoms, deep: int = 1):
        atoms = set(atoms)
        bonds = self._bonds

        for _ in range(deep):
            n = {y for x in atoms for y in bonds[x]} | atoms
            if n == atoms:
                break
            atoms = n
        return self.substructure(atoms)

    def get_mapping(self, other: 'CGRContainer', /, *, automorphism_filter: bool = True,
                    searching_scope: Optional[Collection[int]] = None):
        """
        Get self to other CGR substructure mapping generator.

        :param other: CGR
        :param automorphism_filter: skip matches to the same atoms.
        :param searching_scope: substructure atoms list to localize isomorphism.
        """
        if isinstance(other, CGRContainer):
            return self._get_mapping(other, automorphism_filter=automorphism_filter, searching_scope=searching_scope)
        raise TypeError('CGRContainer expected')

    def __iter__(self):
        return iter(self._atoms)

    def __len__(self):
        return len(self._atoms)


__all__ = ['CGRContainer']
