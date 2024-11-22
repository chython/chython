# -*- coding: utf-8 -*-
#
#  Copyright 2018-2024 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from abc import ABC, abstractmethod
from functools import cached_property
from typing import Dict, Generic, Iterator, Optional, Tuple, TypeVar
from ..exceptions import AtomNotFound, MappingError, BondNotFound


Atom = TypeVar('Atom')
Bond = TypeVar('Bond')


class Graph(Generic[Atom, Bond], ABC):
    __slots__ = ('_atoms', '_bonds', '__dict__')
    __class_cache__ = {}

    _atoms: Dict[int, Atom]
    _bonds: Dict[int, Dict[int, Bond]]

    def __init__(self):
        self._atoms = {}
        self._bonds = {}

    def atom(self, n: int) -> Atom:
        return self._atoms[n]

    def has_atom(self, n: int) -> bool:
        return n in self._atoms

    def atoms(self) -> Iterator[Tuple[int, Atom]]:
        """
        iterate over all atoms
        """
        return iter(self._atoms.items())

    @property
    def atoms_count(self) -> int:
        return len(self._atoms)

    @property
    def atoms_numbers(self) -> Iterator[int]:
        return iter(self._atoms)

    def bond(self, n: int, m: int) -> Bond:
        try:
            return self._bonds[n][m]
        except KeyError as e:
            raise BondNotFound from e

    def has_bond(self, n: int, m: int) -> bool:
        try:
            self._bonds[n]  # check if atom exists
            return n in self._bonds[m]
        except KeyError:
            raise AtomNotFound

    def bonds(self) -> Iterator[Tuple[int, int, Bond]]:
        """
        iterate other all bonds
        """
        seen = set()
        for n, m_bond in self._bonds.items():
            seen.add(n)
            for m, bond in m_bond.items():
                if m not in seen:
                    yield n, m, bond

    @cached_property
    def bonds_count(self) -> int:
        return sum(len(x) for x in self._bonds.values()) // 2

    @abstractmethod
    def add_atom(self, atom: Atom, n: Optional[int] = None) -> int:
        """
        new atom addition
        """
        if n is None:
            n = max(self._atoms, default=0) + 1
        elif not isinstance(n, int):
            raise TypeError('mapping should be integer')
        elif n in self._atoms:
            raise MappingError('atom with same number exists')

        self._atoms[n] = atom
        self._bonds[n] = {}
        self.flush_cache()
        return n

    @abstractmethod
    def add_bond(self, n: int, m: int, bond: Bond):
        """
        Add bond.
        """
        if n == m:
            raise MappingError('atom loops impossible')
        if n not in self._bonds or m not in self._bonds:
            raise AtomNotFound('atoms not found')
        if n in self._bonds[m]:
            raise MappingError('atoms already bonded')

        self._bonds[n][m] = self._bonds[m][n] = bond
        self.flush_cache()

    def copy(self):
        """
        copy of graph
        """
        copy = object.__new__(self.__class__)
        copy._atoms = {n: atom.copy(full=True) for n, atom in self.atoms()}
        copy._bonds = cb = {}
        for n, m_bond in self._bonds.items():
            cb[n] = cbn = {}
            for m, bond in m_bond.items():
                if m in cb:  # bond partially exists. need back-connection.
                    cbn[m] = cb[m][n]
                else:
                    cbn[m] = bond.copy(full=True)
        return copy

    def remap(self, mapping: Dict[int, int]):
        """
        Change atom numbers

        :param mapping: mapping of old numbers to the new
        """
        if len(mapping) != len(set(mapping.values())) or \
                not (self._atoms.keys() - mapping.keys()).isdisjoint(mapping.values()):
            raise ValueError('mapping overlap')

        mg = mapping.get
        self._atoms = {mg(n, n): atom for n, atom in self.atoms()}
        self._bonds = {mg(n, n): {mg(m, m): bond for m, bond in m_bond.items()} for n, m_bond in self._bonds.items()}
        self.flush_cache()

    def union(self, other: 'Graph', *, remap: bool = False, copy: bool = True):
        """
        Merge Graphs into one.

        :param remap: if atoms has collisions then remap other graph atoms else raise exception.
        :param copy: keep original structure and return a new object
        """
        if self._atoms.keys() & other._atoms.keys():
            if not remap:
                raise MappingError('mapping of graphs is not disjoint')
            other = other.copy()
            other.remap({n: i for i, n in enumerate(other, start=max(self._atoms) + 1)})
        else:
            other = other.copy()  # make a copy
        u = self.copy() if copy else self
        u._atoms.update(other._atoms)
        u._bonds.update(other._bonds)
        if not copy:
            self.flush_cache()
        return u

    def flush_cache(self):
        self.__dict__.clear()

    def __copy__(self):
        return self.copy()

    def __or__(self, other):
        """
        G | H is union of graphs
        """
        return self.union(other, remap=True)

    def __ior__(self, other):
        """
        G =| H is union of graphs
        """
        return self.union(other, remap=True, copy=False)

    def __len__(self):
        return len(self._atoms)

    def __iter__(self) -> Iterator[int]:
        return iter(self._atoms)

    def __bool__(self):
        return bool(self._atoms)


__all__ = ['Graph']
