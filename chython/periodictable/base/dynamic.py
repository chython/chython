# -*- coding: utf-8 -*-
#
#  Copyright 2020-2024 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from typing import Type, Optional
from .element import Element


class DynamicElement(ABC):
    __slots__ = ('_charge', '_is_radical', '_p_charge', '_p_is_radical', '_isotope')

    def __init__(self, isotope: Optional[int]):
        self._isotope = isotope
        self._charge = self._p_charge = 0
        self._is_radical = self._p_is_radical = False

    @property
    def isotope(self):
        return self._isotope

    @property
    def atomic_symbol(self) -> str:
        return self.__class__.__name__[7:]

    @property
    @abstractmethod
    def atomic_number(self) -> int:
        """
        Element number
        """

    @classmethod
    def from_symbol(cls, symbol: str) -> Type['DynamicElement']:
        """
        get DynamicElement class by its symbol
        """
        try:
            element = next(x for x in DynamicElement.__subclasses__() if x.__name__ == f'Dynamic{symbol}')
        except StopIteration:
            raise ValueError(f'DynamicElement with symbol "{symbol}" not found')
        return element

    @classmethod
    def from_atomic_number(cls, number: int) -> Type['DynamicElement']:
        """
        get DynamicElement class by its number
        """
        try:
            element = next(x for x in DynamicElement.__subclasses__() if x.atomic_number.fget(None) == number)
        except StopIteration:
            raise ValueError(f'DynamicElement with number "{number}" not found')
        return element

    @classmethod
    def from_atom(cls, atom: 'Element') -> 'DynamicElement':
        """
        get DynamicElement object from Element object
        """
        if not isinstance(atom, Element):
            raise TypeError('Element expected')
        dynamic = object.__new__(cls.from_atomic_number(atom.atomic_number))
        dynamic._isotope = atom.isotope
        dynamic._charge = dynamic._p_charge = atom.charge
        dynamic._is_radical = dynamic._p_is_radical = atom.is_radical
        return dynamic

    @classmethod
    def from_atoms(cls, atom1: 'Element', atom2: 'Element') -> 'DynamicElement':
        """
        get DynamicElement object from pair of Element objects
        """
        if not isinstance(atom1, Element) or not isinstance(atom2, Element):
            raise TypeError('Element expected')
        if atom1.atomic_number != atom2.atomic_number:
            raise ValueError('elements should be of the same type')
        if atom1.isotope != atom2.isotope:
            raise ValueError('elements should be of the same isotope')
        dynamic = object.__new__(cls.from_atomic_number(atom1.atomic_number))
        dynamic._isotope = atom1.isotope
        dynamic._charge = atom1.charge
        dynamic._p_charge = atom2.charge
        dynamic._is_radical = atom1.is_radical
        dynamic._p_is_radical = atom2.is_radical
        return dynamic

    @property
    def charge(self) -> int:
        return self._charge

    @property
    def is_radical(self) -> bool:
        return self._is_radical

    @property
    def p_charge(self) -> int:
        return self._p_charge

    @property
    def p_is_radical(self) -> bool:
        return self._p_is_radical

    def __eq__(self, other):
        """
        compare attached to molecules dynamic elements
        """
        return isinstance(other, DynamicElement) and self.atomic_number == other.atomic_number and \
            self.isotope == other.isotope and self.charge == other.charge and self.is_radical == other.is_radical and \
            self.p_charge == other.p_charge and self.p_is_radical == other.p_is_radical

    def __hash__(self):
        return hash((self.isotope or 0, self.atomic_number, self.charge, self.p_charge,
                     self.is_radical, self.p_is_radical))

    @property
    def is_dynamic(self) -> bool:
        """
        Atom has dynamic features
        """
        return self.charge != self.p_charge or self.is_radical != self.p_is_radical

    def copy(self):
        copy = object.__new__(self.__class__)
        copy._isotope = self.isotope
        copy._charge = self.charge
        copy._is_radical = self.is_radical
        copy._p_is_radical = self.p_is_radical
        copy._p_charge = self.p_charge
        return copy

    def __copy__(self):
        return self.copy()


__all__ = ['DynamicElement']
