# -*- coding: utf-8 -*-
#
#  Copyright 2018-2021 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from ..algorithms.isomorphism import Isomorphism
from ..algorithms.morgan import Morgan
from ..algorithms.rings import Rings
from ..exceptions import AtomNotFound, MappingError


Atom = TypeVar('Atom')
Bond = TypeVar('Bond')


class Graph(Generic[Atom, Bond], Morgan, Rings, Isomorphism, ABC):
    __slots__ = ('_atoms', '_bonds', '_plane', '_charges', '_radicals', '__dict__', '__weakref__')
    __class_cache__ = {}

    _atoms: Dict[int, Atom]
    _bonds: Dict[int, Dict[int, Bond]]
    _charges: Dict[int, int]
    _radicals: Dict[int, bool]
    _plane: Dict[int, Tuple[float, float]]

    def __init__(self):
        self._atoms = {}
        self._bonds = {}
        self._charges = {}
        self._radicals = {}
        self._plane = {}

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
        return self._bonds[n][m]

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
    def add_atom(self, atom: Atom, _map: Optional[int] = None, *, charge: int = 0,
                 is_radical: bool = False, xy: Tuple[float, float] = (0., 0.)) -> int:
        """
        new atom addition
        """
        if _map is None:
            _map = max(self._atoms, default=0) + 1
        elif not isinstance(_map, int):
            raise TypeError('mapping should be integer')
        elif _map in self._atoms:
            raise MappingError('atom with same number exists')
        elif not isinstance(xy, tuple) or len(xy) != 2 or not isinstance(xy[0], float) or not isinstance(xy[1], float):
            raise TypeError('XY should be tuple with 2 float')
        elif not isinstance(is_radical, bool):
            raise TypeError('bool expected')
        elif not isinstance(charge, int):
            raise TypeError('formal charge should be int in range [-4, 4]')
        elif charge > 4 or charge < -4:
            raise ValueError('formal charge should be in range [-4, 4]')

        self._atoms[_map] = atom
        self._charges[_map] = charge
        self._radicals[_map] = is_radical
        self._plane[_map] = xy
        self._bonds[_map] = {}
        atom._attach_to_graph(self, _map)
        self.__dict__.clear()
        return _map

    @abstractmethod
    def add_bond(self, n: int, m: int, bond: Bond):
        """
        Add bond.
        """
        if n == m:
            raise ValueError('atom loops impossible')
        if n not in self._bonds or m not in self._bonds:
            raise AtomNotFound('atoms not found')
        if n in self._bonds[m]:
            raise ValueError('atoms already bonded')

        self._bonds[n][m] = self._bonds[m][n] = bond
        self.__dict__.clear()

    @abstractmethod
    def delete_atom(self, n: int):
        """
        Remove atom.
        """
        del self._atoms[n]
        del self._charges[n]
        del self._radicals[n]
        del self._plane[n]
        sb = self._bonds
        for m in sb.pop(n):
            del sb[m][n]
        self.__dict__.clear()

    @abstractmethod
    def delete_bond(self, n: int, m: int):
        """
        Remove bond.
        """
        del self._bonds[n][m]
        del self._bonds[m][n]
        self.__dict__.clear()

    @abstractmethod
    def copy(self):
        """
        copy of graph
        """
        copy = object.__new__(self.__class__)
        copy._charges = self._charges.copy()
        copy._radicals = self._radicals.copy()
        copy._plane = self._plane.copy()

        copy._bonds = cb = {}
        for n, m_bond in self._bonds.items():
            cb[n] = cbn = {}
            for m, bond in m_bond.items():
                if m in cb:  # bond partially exists. need back-connection.
                    cbn[m] = cb[m][n]
                else:
                    cbn[m] = bond.copy()

        copy._atoms = ca = {}
        for n, atom in self._atoms.items():
            atom = atom.copy()
            ca[n] = atom
            atom._attach_to_graph(copy, n)
        return copy

    @abstractmethod
    def union(self, other: 'Graph'):
        """
        Merge Graphs into one.
        """
        u = self.copy()
        u._charges.update(other._charges)
        u._radicals.update(other._radicals)
        u._plane.update(other._plane)

        ua = u._atoms
        for n, atom in other._atoms.items():
            ua[n] = atom = atom.copy()
            atom._attach_to_graph(u, n)

        ub = u._bonds
        for n, m_bond in other._bonds.items():
            ub[n] = ubn = {}
            for m, bond in m_bond.items():
                if m in ub:  # bond partially exists. need back-connection.
                    ubn[m] = ub[m][n]
                else:
                    ubn[m] = bond.copy()
        return u

    def flush_cache(self):
        self.__dict__.clear()

    def __copy__(self):
        return self.copy()

    def __or__(self, other):
        """
        G | H is union of graphs
        """
        return self.union(other)

    def __len__(self):
        return len(self._atoms)

    def __iter__(self) -> Iterator[int]:
        return iter(self._atoms)

    def __bool__(self):
        return bool(self._atoms)

    def __getstate__(self):
        state = {'atoms': self._atoms, 'bonds': self._bonds, 'plane': self._plane, 'charges': self._charges,
                 'radicals': self._radicals}
        import chython
        if chython.pickle_cache:
            state['cache'] = {k: v for k, v in self.__dict__.items() if k != '__cached_method___hash__'}
        return state

    def __setstate__(self, state):
        self._atoms = state['atoms']
        for n, a in state['atoms'].items():
            a._attach_to_graph(self, n)
        self._charges = state['charges']
        self._radicals = state['radicals']
        self._plane = state['plane']
        self._bonds = state['bonds']
        if 'cache' in state:
            self.__dict__.update(state['cache'])


__all__ = ['Graph']
