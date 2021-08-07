# -*- coding: utf-8 -*-
#
#  Copyright 2020, 2021 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from typing import Optional, Tuple, Dict
from weakref import ref
from ...exceptions import IsConnectedAtom, IsNotConnectedAtom


class Core(ABC):
    __slots__ = ('__isotope', '_graph', '_map', '_backward')

    def __init__(self, isotope: Optional[int] = None):
        """
        Element object with specified isotope

        :param isotope: Isotope number of element
        """
        if isinstance(isotope, int):
            if isotope not in self.isotopes_distribution:
                raise ValueError(f'isotope number {isotope} impossible or not stable for {self.__class__.__name__}')
        elif isotope is not None:
            raise TypeError('integer isotope number required')
        self.__isotope = isotope

    def __repr__(self):
        if self.__isotope:
            return f'{self.__class__.__name__}({self.__isotope})'
        return f'{self.__class__.__name__}()'

    def __getstate__(self):
        return {'isotope': self.__isotope}

    def __setstate__(self, state):
        self.__isotope = state['isotope']

    def __int__(self):
        """
        Same as hash
        """
        return hash(self)

    @property
    @abstractmethod
    def atomic_number(self) -> int:
        """
        Element number
        """

    @property
    def isotope(self) -> Optional[int]:
        """
        Isotope number
        """
        return self.__isotope

    @property
    def atomic_mass(self) -> float:
        mass = self.isotopes_masses
        if self.__isotope is None:
            return sum(x * mass[i] for i, x in self.isotopes_distribution.items())
        return mass[self.__isotope]

    @property
    @abstractmethod
    def isotopes_distribution(self) -> Dict[int, float]:
        """
        Isotopes distribution in earth
        """

    @property
    @abstractmethod
    def isotopes_masses(self) -> Dict[int, float]:
        """
        Isotopes masses
        """

    @property
    @abstractmethod
    def atomic_radius(self) -> float:
        """
        Valence radius of atom
        """

    @property
    def charge(self) -> int:
        """
        Charge of atom
        """
        try:
            return self._graph()._charges[self._map]
        except AttributeError:
            raise IsNotConnectedAtom

    @property
    def is_radical(self) -> bool:
        """
        Radical state of atoms
        """
        try:
            return self._graph()._radicals[self._map]
        except AttributeError:
            raise IsNotConnectedAtom

    @property
    def x(self) -> float:
        """
        X coordinate of atom on 2D plane
        """
        try:
            return self._graph()._plane[self._map][0]
        except AttributeError:
            raise IsNotConnectedAtom

    @property
    def y(self) -> float:
        """
        Y coordinate of atom on 2D plane
        """
        try:
            return self._graph()._plane[self._map][1]
        except AttributeError:
            raise IsNotConnectedAtom

    @property
    def xy(self) -> Tuple[float, float]:
        """
        (X, Y) coordinates of atom on 2D plane
        """
        try:
            return self._graph()._plane[self._map]
        except AttributeError:
            raise IsNotConnectedAtom

    @property
    @abstractmethod
    def neighbors(self):
        """
        Neighbors count of atom
        """

    @property
    @abstractmethod
    def hybridization(self):
        """
        1 - if atom has zero or only single bonded neighbors, 2 - if has only one double bonded neighbor and any amount
        of single bonded, 3 - if has one triple bonded and any amount of double and single bonded neighbors or
        two double bonded and any amount of single bonded neighbors, 4 - if atom in aromatic ring.
        """

    @property
    def in_ring(self) -> bool:
        """
        Atom in any ring.
        """
        try:
            return self._map in self._graph().ring_atoms
        except AttributeError:
            raise IsNotConnectedAtom

    @property
    def ring_sizes(self) -> Tuple[int, ...]:
        """
        Atom rings sizes.
        """
        try:
            return self._graph().atoms_rings_sizes[self._map]
        except AttributeError:
            raise IsNotConnectedAtom
        except KeyError:
            return ()

    def copy(self) -> 'Core':
        """
        Detached from graph copy of element
        """
        copy = object.__new__(self.__class__)
        copy._Core__isotope = self.__isotope
        return copy

    def _attach_to_graph(self, graph, _map):
        try:
            self._graph
        except AttributeError:
            self._graph = ref(graph)
            self._map = _map
        else:
            raise IsConnectedAtom

    def _change_map(self, _map):
        try:
            self._graph
        except AttributeError:
            raise IsNotConnectedAtom
        else:
            self._map = _map


__all__ = ['Core']
