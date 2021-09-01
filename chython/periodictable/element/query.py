# -*- coding: utf-8 -*-
#
#  Copyright 2020, 2021 Ramil Nugmanov <nougmanoff@protonmail.com>
#  Copyright 2021 Dmitrij Zanadvornykh <zandmitrij@gmail.com>
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
from abc import ABC
from typing import Tuple, Type, List, Union
from .core import Core
from .element import Element
from ...exceptions import IsNotConnectedAtom


_inorganic = {'He', 'Ne', 'Ar', 'Kr', 'Xe', 'F', 'Cl', 'Br', 'I', 'C', 'N', 'O',
              'H', 'Si', 'P', 'S', 'Se', 'Ge', 'As', 'Sb', 'Te'}


class Query(Core, ABC):
    __slots__ = ()

    @Core.charge.setter
    def charge(self, charge):
        if not isinstance(charge, int):
            raise TypeError('formal charge should be int in range [-4, 4]')
        elif charge > 4 or charge < -4:
            raise ValueError('formal charge should be in range [-4, 4]')
        try:
            g = self._graph()
            g._charges[self._map] = charge
        except AttributeError:
            raise IsNotConnectedAtom
        else:
            g.flush_cache()

    @Core.is_radical.setter
    def is_radical(self, is_radical):
        if not isinstance(is_radical, bool):
            raise TypeError('bool expected')
        try:
            g = self._graph()
            g._radicals[self._map] = is_radical
        except AttributeError:
            raise IsNotConnectedAtom
        else:
            g.flush_cache()

    @property
    def neighbors(self) -> Tuple[int, ...]:
        try:
            return self._graph()._neighbors[self._map]
        except AttributeError:
            raise IsNotConnectedAtom

    @property
    def hybridization(self):
        try:
            return self._graph()._hybridizations[self._map]
        except AttributeError:
            raise IsNotConnectedAtom

    @neighbors.setter
    def neighbors(self, neighbors):
        try:
            g = self._graph()
            g._neighbors[self._map] = g._validate_neighbors(neighbors)
        except AttributeError:
            raise IsNotConnectedAtom
        else:
            g.flush_cache()

    @hybridization.setter
    def hybridization(self, hybridization):
        try:
            g = self._graph()
            g._hybridizations[self._map] = g._validate_hybridization(hybridization)
        except AttributeError:
            raise IsNotConnectedAtom
        else:
            g.flush_cache()

    @property
    def heteroatoms(self) -> Tuple[int, ...]:
        try:
            return self._graph()._heteroatoms[self._map]
        except AttributeError:
            raise IsNotConnectedAtom

    @heteroatoms.setter
    def heteroatoms(self, heteroatoms):
        try:
            g = self._graph()
            g._heteroatoms[self._map] = g._validate_neighbors(heteroatoms)
        except AttributeError:
            raise IsNotConnectedAtom
        else:
            g.flush_cache()

    @property
    def ring_sizes(self) -> Tuple[int, ...]:
        """
        Atom rings sizes.
        """
        try:
            return self._graph()._rings_sizes[self._map]
        except AttributeError:
            raise IsNotConnectedAtom
        except KeyError:
            return ()

    @ring_sizes.setter
    def ring_sizes(self, ring_sizes):
        try:
            g = self._graph()
            g._rings_sizes[self._map] = g._validate_rings(ring_sizes)
        except AttributeError:
            raise IsNotConnectedAtom
        else:
            g.flush_cache()

    @property
    def implicit_hydrogens(self) -> Tuple[int, ...]:
        try:
            return self._graph()._hydrogens[self._map]
        except AttributeError:
            raise IsNotConnectedAtom

    @implicit_hydrogens.setter
    def implicit_hydrogens(self, implicit_hydrogens):
        try:
            g = self._graph()
            g._hydrogens[self._map] = g._validate_neighbors(implicit_hydrogens)
        except AttributeError:
            raise IsNotConnectedAtom
        else:
            g.flush_cache()


class QueryElement(Query, ABC):
    __slots__ = ()

    @property
    def atomic_symbol(self) -> str:
        return self.__class__.__name__[5:]

    @classmethod
    def from_symbol(cls, symbol: str) -> Type[Union['QueryElement', 'AnyElement', 'AnyMetal']]:
        """
        get Element class by its symbol
        """
        if symbol == 'A':
            return AnyElement
        elif symbol == 'M':
            return AnyMetal
        try:
            element = next(x for x in QueryElement.__subclasses__() if x.__name__ == f'Query{symbol}')
        except StopIteration:
            raise ValueError(f'QueryElement with symbol "{symbol}" not found')
        return element

    @classmethod
    def from_atomic_number(cls, number: int) -> Type['QueryElement']:
        """
        get Element class by its number
        """
        try:
            element = next(x for x in QueryElement.__subclasses__() if x.atomic_number.fget(None) == number)
        except StopIteration:
            raise ValueError(f'QueryElement with number "{number}" not found')
        return element

    @classmethod
    def from_atom(cls, atom: Union['Element', 'Query']) -> 'Query':
        """
        get QueryElement or AnyElement object from Element object or copy of QueryElement or AnyElement
        """
        if isinstance(atom, Element):
            return cls.from_atomic_number(atom.atomic_number)(atom.isotope)
        elif not isinstance(atom, Query):
            raise TypeError('Element or Query expected')
        return atom.copy()

    def __eq__(self, other):
        """
        compare attached to molecules elements and query elements
        """
        if isinstance(other, Element):
            if self.atomic_number == other.atomic_number and self.charge == other.charge and \
                    self.is_radical == other.is_radical:
                if self.isotope and self.isotope != other.isotope:
                    return False
                if self.neighbors and other.neighbors not in self.neighbors:
                    return False
                if self.hybridization and other.hybridization not in self.hybridization:
                    return False
                if self.ring_sizes:
                    if self.ring_sizes[0]:
                        if set(self.ring_sizes).isdisjoint(other.ring_sizes):
                            return False
                    elif other.ring_sizes:  # not in ring expected
                        return False
                if self.implicit_hydrogens and other.implicit_hydrogens not in self.implicit_hydrogens:
                    return False
                if self.heteroatoms and other.heteroatoms not in self.heteroatoms:
                    return False
                return True
        elif isinstance(other, QueryElement) and self.atomic_number == other.atomic_number and \
                self.isotope == other.isotope and self.charge == other.charge and self.is_radical == other.is_radical \
                and self.neighbors == other.neighbors and self.hybridization == other.hybridization \
                and self.ring_sizes == other.ring_sizes and self.implicit_hydrogens == other.implicit_hydrogens \
                and self.heteroatoms == other.heteroatoms:
            # equal query element has equal query marks
            return True
        return False

    def __hash__(self):
        return hash((self.isotope or 0, self.atomic_number, self.charge, self.is_radical, self.neighbors,
                     self.hybridization, self.ring_sizes, self.implicit_hydrogens, self.heteroatoms))


class AnyElement(Query):
    __slots__ = ()

    def __init__(self, *args, **kwargs):
        super().__init__()

    @property
    def atomic_symbol(self) -> str:
        return 'A'

    @property
    def atomic_number(self) -> int:
        return 0

    def __eq__(self, other):
        """
        Compare attached to molecules elements and query elements
        """
        if isinstance(other, Element):
            if self.charge == other.charge and self.is_radical == other.is_radical:
                if self.neighbors and other.neighbors not in self.neighbors:
                    return False
                if self.hybridization and other.hybridization not in self.hybridization:
                    return False
                if self.ring_sizes:
                    if self.ring_sizes[0]:
                        if set(self.ring_sizes).isdisjoint(other.ring_sizes):
                            return False
                    elif other.ring_sizes:  # not in ring expected
                        return False
                if self.implicit_hydrogens and other.implicit_hydrogens not in self.implicit_hydrogens:
                    return False
                if self.heteroatoms and other.heteroatoms not in self.heteroatoms:
                    return False
                return True
        elif isinstance(other, AnyMetal):
            return False
        elif isinstance(other, Query) and self.charge == other.charge and self.is_radical == other.is_radical \
                and self.neighbors == other.neighbors and self.hybridization == other.hybridization \
                and self.ring_sizes == other.ring_sizes and self.implicit_hydrogens == other.implicit_hydrogens \
                and self.heteroatoms == other.heteroatoms:
            return True
        return False

    def __hash__(self):
        return hash((self.charge, self.is_radical, self.neighbors, self.hybridization, self.ring_sizes,
                     self.implicit_hydrogens, self.heteroatoms))


class AnyMetal(Query):
    """
    Charge and radical ignored any metal. Rings, hydrogens and heteroatoms count also ignored.

    Class designed for d-elements matching in standardization.
    """
    def __init__(self, *args, **kwargs):
        super().__init__()

    @property
    def atomic_symbol(self) -> str:
        return 'M'

    @property
    def atomic_number(self) -> int:
        return 0

    def __eq__(self, other):
        if isinstance(other, Element):
            if other.atomic_symbol not in _inorganic:
                if self.neighbors and other.neighbors not in self.neighbors:
                    return False
                if self.hybridization and other.hybridization not in self.hybridization:
                    return False
                return True
        elif isinstance(other, AnyMetal) and self.neighbors == other.neighbors \
                and self.hybridization == other.hybridization:
            return True
        return False

    def __hash__(self):
        return hash((self.neighbors, self.hybridization))


class ListElement(Query):
    __slots__ = ('_elements', '_numbers')

    def __init__(self, elements: List[str], *args, **kwargs):
        """
        Elements list
        """
        super().__init__()
        self._elements = tuple(elements)
        self._numbers = tuple(x.atomic_number.fget(None) for x in Element.__subclasses__() if x.__name__ in elements)

    @property
    def atomic_symbol(self) -> str:
        return ','.join(self._elements)

    @property
    def atomic_number(self) -> int:
        return 0

    def copy(self):
        copy = super().copy()
        copy._elements = self._elements
        copy._numbers = self._numbers
        return copy

    def __eq__(self, other):
        """
        Compare attached to molecules elements and query elements
        """
        if isinstance(other, Element):
            if other.atomic_number in self._numbers:
                if self.charge != other.charge or self.is_radical != other.is_radical:
                    return False
                if self.neighbors and other.neighbors not in self.neighbors:
                    return False
                if self.hybridization and other.hybridization not in self.hybridization:
                    return False
                if self.ring_sizes:
                    if self.ring_sizes[0]:
                        if set(self.ring_sizes).isdisjoint(other.ring_sizes):
                            return False
                    elif other.ring_sizes:  # not in ring expected
                        return False
                if self.implicit_hydrogens and other.implicit_hydrogens not in self.implicit_hydrogens:
                    return False
                if self.heteroatoms and other.heteroatoms not in self.heteroatoms:
                    return False
                return True
        elif isinstance(other, (AnyElement, AnyMetal)):
            return False
        elif isinstance(other, Query) and self.charge == other.charge and self.is_radical == other.is_radical \
                and self.neighbors == other.neighbors and self.hybridization == other.hybridization \
                and self.ring_sizes == other.ring_sizes and self.implicit_hydrogens == other.implicit_hydrogens \
                and self.heteroatoms == other.heteroatoms:
            if isinstance(other, ListElement):
                return self._numbers == other._numbers
            return other.atomic_number in self._numbers
        return False

    def __hash__(self):
        return hash((self._numbers, self.charge, self.is_radical, self.neighbors, self.hybridization,
                     self.ring_sizes, self.implicit_hydrogens, self.heteroatoms))

    def __getstate__(self):
        state = super().__getstate__()
        state['elements'] = self._elements
        return state

    def __setstate__(self, state):
        self._elements = state['elements']
        self._numbers = tuple(x.atomic_number.fget(None) for x in Element.__subclasses__()
                              if x.__name__ in state['elements'])
        super().__setstate__(state)

    def __repr__(self):
        return f'{self.__class__.__name__}([{",".join(self._elements)}])'


__all__ = ['Query', 'QueryElement', 'AnyElement', 'AnyMetal', 'ListElement']
