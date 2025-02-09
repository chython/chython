# -*- coding: utf-8 -*-
#
#  Copyright 2020-2025 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from abc import ABC, abstractmethod
from functools import cached_property
from typing import Tuple, Type, List, Union, Optional
from .element import Element
from .groups import GroupXVIII


def _validate(value, prop):
    if value is None:
        return ()
    elif isinstance(value, int):
        if value < 0 or value > 14:
            raise ValueError(f'{prop} should be in range [0, 14]')
        return (value,)
    elif isinstance(value, (tuple, list)):
        if not all(isinstance(x, int) for x in value):
            raise TypeError(f'{prop} should be list or tuple of ints')
        if any(x < 0 or x > 14 for x in value):
            raise ValueError(f'{prop} should be in range [0, 14]')
        if len(set(value)) != len(value):
            raise ValueError(f'{prop} should be unique')
        return tuple(sorted(value))
    else:
        raise TypeError(f'{prop} should be int or list or tuple of ints')


class Query(ABC):
    __slots__ = ('_neighbors', '_hybridization', '_masked')

    def __init__(self, neighbors: Union[int, Tuple[int, ...], None] = None,
                 hybridization: Union[int, Tuple[int, ...], None] = None, masked: bool = False):
        self.neighbors = neighbors
        self.hybridization = hybridization
        self.masked = masked

    @property
    @abstractmethod
    def atomic_symbol(self) -> str:
        ...

    @property
    def neighbors(self) -> Tuple[int, ...]:
        return self._neighbors

    @neighbors.setter
    def neighbors(self, value):
        self._neighbors = _validate(value, 'neighbors')

    @property
    def hybridization(self) -> Tuple[int, ...]:
        return self._hybridization

    @hybridization.setter
    def hybridization(self, value):
        if value is None:
            self._hybridization = ()
        elif isinstance(value, int):
            if value < 1 or value > 4:
                raise ValueError('hybridization should be in range [1, 4]')
            self._hybridization = (value,)
        elif isinstance(value, (tuple, list)):
            if not all(isinstance(h, int) for h in value):
                raise TypeError('hybridizations should be list or tuple of ints')
            if any(h < 1 or h > 4 for h in value):
                raise ValueError('hybridizations should be in range [1, 4]')
            if len(set(value)) != len(value):
                raise ValueError('hybridizations should be unique')
            self._hybridization = tuple(sorted(value))
        else:
            raise TypeError('hybridization should be int or list or tuple of ints')

    @property
    def masked(self):
        return self._masked

    @masked.setter
    def masked(self, value):
        if not isinstance(value, bool):
            raise TypeError('masked should be bool')
        self._masked = value

    def copy(self, full=False):
        copy = object.__new__(self.__class__)
        copy._neighbors = self.neighbors
        copy._hybridization = self.hybridization

        copy._masked = self.masked if full else False
        return copy

    def __copy__(self):
        return self.copy()

    def __repr__(self):
        return f'{self.__class__.__name__}()'


class ExtendedQuery(Query, ABC):
    __slots__ = ('_charge', '_is_radical', '_heteroatoms', '_ring_sizes', '_implicit_hydrogens', '_stereo')

    def __init__(self, charge: int = 0, is_radical: bool = False, heteroatoms: Union[int, Tuple[int, ...], None] = None,
                 ring_sizes: Union[int, Tuple[int, ...], None] = None,
                 implicit_hydrogens: Union[int, Tuple[int, ...], None] = None, stereo: Optional[bool] = None, **kwargs):
        super().__init__(**kwargs)
        self.charge = charge
        self.is_radical = is_radical
        self.heteroatoms = heteroatoms
        self.ring_sizes = ring_sizes
        self.implicit_hydrogens = implicit_hydrogens
        self.stereo = stereo

    @property
    def charge(self) -> int:
        """
        Charge of atom
        """
        return self._charge

    @charge.setter
    def charge(self, value: int):
        if not isinstance(value, int):
            raise TypeError('formal charge should be int in range [-4, 4]')
        elif value > 4 or value < -4:
            raise ValueError('formal charge should be in range [-4, 4]')
        self._charge = value

    @property
    def is_radical(self) -> bool:
        """
        Radical state of atoms
        """
        return self._is_radical

    @is_radical.setter
    def is_radical(self, value: bool):
        if not isinstance(value, bool):
            raise TypeError('bool expected')
        self._is_radical = value

    @property
    def heteroatoms(self) -> Tuple[int, ...]:
        return self._heteroatoms

    @heteroatoms.setter
    def heteroatoms(self, value):
        self._heteroatoms = _validate(value, 'heteroatoms')

    @property
    def implicit_hydrogens(self) -> Tuple[int, ...]:
        return self._implicit_hydrogens

    @implicit_hydrogens.setter
    def implicit_hydrogens(self, value):
        self._implicit_hydrogens = _validate(value, 'implicit hydrogens')

    @property
    def ring_sizes(self) -> Tuple[int, ...]:
        """
        Atom rings sizes.
        """
        return self._ring_sizes

    @ring_sizes.setter
    def ring_sizes(self, value):
        if value is None:
            self._ring_sizes = ()
        elif isinstance(value, int):
            if value < 3 and value != 0:
                raise ValueError('rings should be greater or equal 3. ring equal to zero is no ring atom mark')
            self._ring_sizes = (value,)
        elif isinstance(value, (tuple, list)):
            if not all(isinstance(x, int) for x in value):
                raise TypeError('rings should be list or tuple of ints')
            if any(x < 3 for x in value):
                raise ValueError('rings should be greater or equal 3')
            if len(set(value)) != len(value):
                raise ValueError('rings should be unique')
            self._ring_sizes = tuple(sorted(value))
        else:
            raise TypeError('rings should be int or list or tuple of ints')

    @property
    def stereo(self):
        return self._stereo

    @stereo.setter
    def stereo(self, value: Optional[bool]):
        if value is not None and not isinstance(value, bool):
            raise TypeError('stereo should be bool')
        self._stereo = value

    def copy(self, full=False):
        copy = super().copy(full=full)
        copy._charge = self.charge
        copy._is_radical = self.is_radical
        copy._heteroatoms = self.heteroatoms
        copy._implicit_hydrogens = self.implicit_hydrogens
        copy._ring_sizes = self.ring_sizes

        copy._stereo = self.stereo if full else None
        return copy


class AnyMetal(Query):
    """
    Charge and radical ignored any metal. Rings, hydrogens and heteroatoms count also ignored.

    Class designed for d-elements matching in standardization.
    """
    __slots__ = ()

    @property
    def atomic_symbol(self) -> str:
        return 'M'

    def __eq__(self, other):
        if not isinstance(other, Element):
            return False
        if other.is_forming_single_bonds or isinstance(other, GroupXVIII):
            return False
        if self.neighbors and other.neighbors not in self.neighbors:
            return False
        if self.hybridization and other.hybridization not in self.hybridization:
            return False
        return True


class AnyElement(ExtendedQuery):
    __slots__ = ()

    @property
    def atomic_symbol(self) -> str:
        return 'A'

    def __eq__(self, other):
        """
        Compare attached to molecules elements and query elements
        """
        if not isinstance(other, Element):
            return False
        if self.charge != other.charge:
            return False
        if self.is_radical != other.is_radical:
            return False
        if self.neighbors and other.neighbors not in self.neighbors:
            return False
        if self.hybridization and other.hybridization not in self.hybridization:
            return False
        if self.ring_sizes:
            if self.ring_sizes[0]:
                if other.ring_sizes.isdisjoint(self.ring_sizes):
                    return False
            elif other.ring_sizes:  # not in ring expected
                return False
        if self.implicit_hydrogens and other.implicit_hydrogens not in self.implicit_hydrogens:
            return False
        if self.heteroatoms and other.heteroatoms not in self.heteroatoms:
            return False
        return True


class ListElement(ExtendedQuery):
    __slots__ = ('_elements', '__dict__')

    def __init__(self, elements: List[str], **kwargs):
        """
        Elements list
        """
        if not isinstance(elements, (list, tuple)) or not elements:
            raise ValueError('invalid elements list')
        tmp = []
        for x in elements:
            if isinstance(x, int):
                tmp.append(Element.from_atomic_number(x).__name__)
            elif isinstance(x, str):
                tmp.append(Element.from_symbol(x).__name__)
            else:
                raise ValueError(f'invalid element: {x}')
        super().__init__(**kwargs)
        self._elements = tuple(tmp)

    @property
    def atomic_symbol(self) -> str:
        return ','.join(self._elements)

    @cached_property
    def atomic_numbers(self):
        return tuple(x.atomic_number.fget(None) for x in Element.__subclasses__() if x.__name__ in self._elements)

    def copy(self, full=False):
        copy = super().copy(full=full)
        copy._elements = self._elements
        return copy

    def __eq__(self, other):
        """
        Compare attached to molecules elements and query elements
        """
        if not isinstance(other, Element):
            return False
        if other.atomic_number not in self.atomic_numbers:
            return False
        if self.charge != other.charge:
            return False
        if self.is_radical != other.is_radical:
            return False
        if self.neighbors and other.neighbors not in self.neighbors:
            return False
        if self.hybridization and other.hybridization not in self.hybridization:
            return False
        if self.ring_sizes:
            if self.ring_sizes[0]:
                if other.ring_sizes.isdisjoint(self.ring_sizes):
                    return False
            elif other.ring_sizes:  # not in ring expected
                return False
        if self.implicit_hydrogens and other.implicit_hydrogens not in self.implicit_hydrogens:
            return False
        if self.heteroatoms and other.heteroatoms not in self.heteroatoms:
            return False
        return True

    def __repr__(self):
        return f'{self.__class__.__name__}([{self.atomic_symbol}])'


class QueryElement(ExtendedQuery, ABC):
    __slots__ = ('_isotope',)

    def __init__(self, isotope: Optional[int] = None, **kwargs):
        super().__init__(**kwargs)
        self.isotope = isotope

    def __repr__(self):
        if self.isotope:
            return f'{self.__class__.__name__}({self.isotope})'
        return f'{self.__class__.__name__}()'

    @property
    def atomic_symbol(self) -> str:
        return self.__class__.__name__[5:]

    @property
    @abstractmethod
    def atomic_number(self) -> int:
        """
        Element number
        """

    @property
    def isotope(self):
        return self._isotope

    @isotope.setter
    def isotope(self, value: Optional[int]):
        if value is not None and not isinstance(value, int):
            raise TypeError('isotope must be an int')
        self._isotope = value

    @property
    @abstractmethod
    def mdl_isotope(self) -> int:
        ...

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
    def from_atom(cls, atom: 'Element', neighbors=False, hybridization=False, heteroatoms=False,
                  hydrogens=False, ring_sizes=False, stereo=False) -> 'QueryElement':
        """
        get QueryElement or AnyElement object from Element object or copy of QueryElement or AnyElement
        """
        if not isinstance(atom, Element):
            raise TypeError('Element or Query expected')

        # transfer true atomic props
        query = cls.from_atomic_number(atom.atomic_number)(atom.isotope)
        query._charge = atom.charge
        query._is_radical = atom.is_radical

        if neighbors:
            query._neighbors = (atom.neighbors,)
        if hybridization:
            query._hybridization = (atom.hybridization,)
        if heteroatoms:
            query._heteroatoms = (atom.heteroatoms,)
        if ring_sizes:
            query._ring_sizes = atom.ring_sizes
        if hydrogens and atom.implicit_hydrogens is not None:
            query._implicit_hydrogens = (atom.implicit_hydrogens,)
        if stereo:
            query._stereo = atom.stereo
        return query

    def copy(self, full=False):
        copy = super().copy(full=full)
        copy._isotope = self.isotope
        return copy

    def __eq__(self, other):
        """
        compare attached to molecules elements and query elements
        """
        if not isinstance(other, Element):
            return False
        if self.atomic_number != other.atomic_number:
            return False
        if self.charge != other.charge:
            return False
        if self.is_radical != other.is_radical:
            return False
        if self.isotope and self.isotope != other.isotope:
            return False
        if self.neighbors and other.neighbors not in self.neighbors:
            return False
        if self.hybridization and other.hybridization not in self.hybridization:
            return False
        if self.ring_sizes:
            if self.ring_sizes[0]:
                if other.ring_sizes.isdisjoint(self.ring_sizes):
                    return False
            elif other.ring_sizes:  # not in ring expected
                return False
        if self.implicit_hydrogens and other.implicit_hydrogens not in self.implicit_hydrogens:
            return False
        if self.heteroatoms and other.heteroatoms not in self.heteroatoms:
            return False
        return True


__all__ = ['Query', 'ExtendedQuery', 'QueryElement', 'AnyElement', 'AnyMetal', 'ListElement']
