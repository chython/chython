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
from typing import Tuple, Dict, Type, List, Union
from .core import Core
from .element import Element
from ...exceptions import IsNotConnectedAtom


class Query(Core):
    __slots__ = ()

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
            g.flush_cache()
        except AttributeError:
            raise IsNotConnectedAtom

    @hybridization.setter
    def hybridization(self, hybridization):
        try:
            g = self._graph()
            g._hybridizations[self._map] = g._validate_hybridization(hybridization)
            g.flush_cache()
        except AttributeError:
            raise IsNotConnectedAtom

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
            g.flush_cache()
        except AttributeError:
            raise IsNotConnectedAtom

    @property
    def in_ring(self) -> bool:
        """
        Atom in any ring.
        """
        try:
            return bool(self._graph()._rings_sizes[self._map]) or self._map in self._graph().ring_atoms
        except AttributeError:
            raise IsNotConnectedAtom

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
            g.flush_cache()
        except AttributeError:
            raise IsNotConnectedAtom

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
            g.flush_cache()
        except AttributeError:
            raise IsNotConnectedAtom


class QueryElement(Query):
    __slots__ = ()

    @property
    def atomic_symbol(self) -> str:
        return self.__class__.__name__[5:]

    @classmethod
    def from_symbol(cls, symbol: str) -> Type[Union['QueryElement', 'AnyElement']]:
        """
        get Element class by its symbol
        """
        if symbol == 'A':
            return AnyElement
        try:
            element = next(x for x in QueryElement.__subclasses__() if x.__name__ == f'Query{symbol}')
        except StopIteration:
            raise ValueError(f'QueryElement with symbol "{symbol}" not found')
        return element

    @classmethod
    def from_atomic_number(cls, number: int) -> Type[Union['QueryElement', 'AnyElement']]:
        """
        get Element class by its number
        """
        if number == 0:
            return AnyElement
        try:
            element = next(x for x in QueryElement.__subclasses__() if x.atomic_number.fget(None) == number)
        except StopIteration:
            raise ValueError(f'QueryElement with number "{number}" not found')
        return element

    @classmethod
    def from_atom(cls, atom: Union['Element', 'QueryElement', 'AnyElement']) -> Union['QueryElement', 'AnyElement']:
        """
        get QueryElement or AnyElement object from Element object or copy of QueryElement or AnyElement
        """
        if isinstance(atom, Element):
            return cls.from_atomic_number(atom.atomic_number)(atom.isotope)
        elif not isinstance(atom, (QueryElement, AnyElement)):
            raise TypeError('Element, QueryElement or AnyElement expected')
        return atom.copy()

    @Core.charge.setter
    def charge(self, charge):
        try:
            g = self._graph()
            # remove stereo data
            if self._map in g._atoms_stereo:
                del g._atoms_stereo[self._map]
            if g._cis_trans_stereo:
                for nm, path in g._stereo_cis_trans_paths.items():
                    if self._map in path and nm in g._cis_trans_stereo:
                        del g._cis_trans_stereo[nm]
            if g._allenes_stereo:
                for c, path in g._stereo_allenes_paths.items():
                    if self._map in path and c in g._allenes_stereo:
                        del g._allenes_stereo[c]
            g._charges[self._map] = g._validate_charge(charge)
            g.flush_cache()
        except AttributeError:
            raise IsNotConnectedAtom

    @Core.is_radical.setter
    def is_radical(self, is_radical):
        try:
            g = self._graph()
            # remove stereo data
            if self._map in g._atoms_stereo:
                del g._atoms_stereo[self._map]
            if g._cis_trans_stereo:
                for nm, path in g._stereo_cis_trans_paths.items():
                    if self._map in path and nm in g._cis_trans_stereo:
                        del g._cis_trans_stereo[nm]
            if g._allenes_stereo:
                for c, path in g._stereo_allenes_paths.items():
                    if self._map in path and c in g._allenes_stereo:
                        del g._allenes_stereo[c]
            g._radicals[self._map] = g._validate_radical(is_radical)
            g.flush_cache()
        except AttributeError:
            raise IsNotConnectedAtom

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

    @property
    def isotopes_distribution(self) -> Dict[int, float]:
        return {}

    @property
    def isotopes_masses(self) -> Dict[int, float]:
        return {}

    @property
    def atomic_radius(self):
        return 0.5

    @Core.charge.setter
    def charge(self, charge):
        try:
            g = self._graph()
            g._charges[self._map] = g._validate_charge(charge)
            g.flush_cache()
        except AttributeError:
            raise IsNotConnectedAtom

    @Core.is_radical.setter
    def is_radical(self, is_radical):
        try:
            g = self._graph()
            g._radicals[self._map] = g._validate_radical(is_radical)
            g.flush_cache()
        except AttributeError:
            raise IsNotConnectedAtom

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
        elif isinstance(other, Query) and self.charge == other.charge and self.is_radical == other.is_radical \
                and self.neighbors == other.neighbors and self.hybridization == other.hybridization \
                and self.ring_sizes == other.ring_sizes and self.implicit_hydrogens == other.implicit_hydrogens \
                and self.heteroatoms == other.heteroatoms:
            return True
        return False

    def __hash__(self):
        return hash((self.charge, self.is_radical, self.neighbors, self.hybridization, self.ring_sizes,
                     self.implicit_hydrogens, self.heteroatoms))


class AnyMetal(AnyElement):
    """
    Charge and radical ignored any metal. Rings, hydrogens and heteroatoms count also ignored.

    Class designed for d-elements matching in standardization.
    """
    @property
    def atomic_symbol(self) -> str:
        return 'M'

    def __eq__(self, other):
        if isinstance(other, Element):
            if other.atomic_symbol not in {'He', 'Ne', 'Ar', 'Kr', 'Xe', 'F', 'Cl', 'Br', 'I', 'C', 'N', 'O', 'H', 'Si',
                                           'P', 'S', 'Se', 'Ge', 'As', 'Sb', 'Te'}:
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


class ListElement(AnyElement):
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
        elif isinstance(other, Query) and self.charge == other.charge and self.is_radical == other.is_radical \
                and self.neighbors == other.neighbors and self.hybridization == other.hybridization \
                and self.ring_sizes == other.ring_sizes and self.implicit_hydrogens == other.implicit_hydrogens \
                and self.heteroatoms == other.heteroatoms:
            if isinstance(other, ListElement):
                return self._numbers == other._numbers
            if isinstance(other, AnyElement):
                return True
            return other.atomic_number in self._numbers
        return False

    def copy(self) -> 'ListElement':
        """
        Detached from graph copy of element
        """
        copy = object.__new__(self.__class__)
        copy._Core__isotope = self._Core__isotope
        copy._elements = self._elements
        copy._numbers = self._numbers
        return copy

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


__all__ = ['QueryElement', 'AnyElement', 'AnyMetal', 'ListElement']
