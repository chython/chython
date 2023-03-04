# -*- coding: utf-8 -*-
#
#  Copyright 2018-2023 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from ..algorithms.morgan import Morgan
from ..algorithms.rings import Rings
from ..exceptions import AtomNotFound, MappingError, BondNotFound


Atom = TypeVar('Atom')
Bond = TypeVar('Bond')


class Graph(Generic[Atom, Bond], Morgan, Rings, ABC):
    __slots__ = ('_atoms', '_bonds', '_charges', '_radicals', '_atoms_stereo', '_cis_trans_stereo', '_allenes_stereo',
                 '__dict__', '__weakref__')
    __class_cache__ = {}

    _atoms: Dict[int, Atom]
    _bonds: Dict[int, Dict[int, Bond]]
    _charges: Dict[int, int]
    _radicals: Dict[int, bool]
    _atoms_stereo: Dict[int, bool]
    _allenes_stereo: Dict[int, bool]
    _cis_trans_stereo: Dict[Tuple[int, int], bool]

    def __init__(self):
        self._atoms = {}
        self._bonds = {}
        self._charges = {}
        self._radicals = {}
        self._atoms_stereo = {}
        self._allenes_stereo = {}
        self._cis_trans_stereo = {}

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
    def add_atom(self, atom: Atom, n: Optional[int] = None, *, charge: int = 0, is_radical: bool = False) -> int:
        """
        new atom addition
        """
        if n is None:
            n = max(self._atoms, default=0) + 1
        elif not isinstance(n, int):
            raise TypeError('mapping should be integer')
        elif n in self._atoms:
            raise MappingError('atom with same number exists')
        elif not isinstance(is_radical, bool):
            raise TypeError('bool expected')
        elif not isinstance(charge, int):
            raise TypeError('formal charge should be int in range [-4, 4]')
        elif charge > 4 or charge < -4:
            raise ValueError('formal charge should be in range [-4, 4]')

        atom._attach_graph(self, n)
        self._atoms[n] = atom
        self._charges[n] = charge
        self._radicals[n] = is_radical
        self._bonds[n] = {}
        self.__dict__.clear()
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
        self.__dict__.clear()

    @abstractmethod
    def copy(self):
        """
        copy of graph
        """
        copy = object.__new__(self.__class__)
        copy._charges = self._charges.copy()
        copy._radicals = self._radicals.copy()

        copy._atoms = ca = {}
        for n, atom in self._atoms.items():
            atom = atom.copy()
            ca[n] = atom
            atom._attach_graph(copy, n)
        return copy

    @abstractmethod
    def remap(self, mapping: Dict[int, int], *, copy=False):
        """
        Change atom numbers

        :param mapping: mapping of old numbers to the new
        :param copy: keep original graph
        """
        if len(mapping) != len(set(mapping.values())) or \
                not (self._atoms.keys() - mapping.keys()).isdisjoint(mapping.values()):
            raise ValueError('mapping overlap')

        mg = mapping.get
        sc = self._charges
        sr = self._radicals

        if copy:
            h = self.__class__()
            ha = h._atoms
            hc = h._charges
            hr = h._radicals
            has = h._atoms_stereo
            hal = h._allenes_stereo
            hcs = h._cis_trans_stereo

            for n, atom in self._atoms.items():
                m = mg(n, n)
                atom = atom.copy()
                ha[m] = atom
                atom._attach_graph(h, m)
        else:
            ha = {}
            hc = {}
            hr = {}
            has = {}
            hal = {}
            hcs = {}

            for n, atom in self._atoms.items():
                m = mg(n, n)
                ha[m] = atom
                atom._change_map(m)  # change mapping number

        for n in self._atoms:
            m = mg(n, n)
            hc[m] = sc[n]
            hr[m] = sr[n]

        for n, stereo in self._atoms_stereo.items():
            has[mg(n, n)] = stereo
        for n, stereo in self._allenes_stereo.items():
            hal[mg(n, n)] = stereo
        for (n, m), stereo in self._cis_trans_stereo.items():
            hcs[(mg(n, n), mg(m, m))] = stereo

        if copy:
            return h  # noqa

        self._atoms = ha
        self._charges = hc
        self._radicals = hr
        self._atoms_stereo = has
        self._allenes_stereo = hal
        self._cis_trans_stereo = hcs
        self.flush_cache()
        return self

    @abstractmethod
    def union(self, other: 'Graph', *, remap: bool = False, copy: bool = True):
        """
        Merge Graphs into one.

        :param remap: if atoms has collisions then remap other graph atoms else raise exception.
        :param copy: keep original structure and return new object
        """
        if self._atoms.keys() & other._atoms.keys():
            if remap:
                other = other.remap({n: i for i, n in enumerate(other, start=max(self._atoms) + 1)}, copy=True)
            else:
                raise MappingError('mapping of graphs is not disjoint')

        u = self.copy() if copy else self
        u._charges.update(other._charges)
        u._radicals.update(other._radicals)

        ua = u._atoms
        for n, atom in other._atoms.items():
            ua[n] = atom = atom.copy()
            atom._attach_graph(u, n)

        u._atoms_stereo.update(other._atoms_stereo)
        u._allenes_stereo.update(other._allenes_stereo)
        u._cis_trans_stereo.update(other._cis_trans_stereo)
        return u, other

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

    def __getstate__(self):
        state = {'atoms': self._atoms, 'bonds': self._bonds, 'charges': self._charges,
                 'radicals': self._radicals}
        from chython import pickle_cache

        if pickle_cache:
            state['cache'] = {k: v for k, v in self.__dict__.items() if k != '__cached_method___hash__'}
        return state

    def __setstate__(self, state):
        self._atoms = state['atoms']
        for n, a in state['atoms'].items():
            a._attach_graph(self, n)
        self._charges = state['charges']
        self._radicals = state['radicals']
        self._bonds = state['bonds']
        if 'cache' in state:
            self.__dict__.update(state['cache'])


__all__ = ['Graph']
