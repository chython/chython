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
from typing import Type, Union, Optional
from .element import Element


class DynamicElement(ABC):
    __slots__ = ('_charge', '_is_radical', '_p_charge', '_p_is_radical', '_isotope')

    def __init__(self, isotope: Optional[int]):
        self._isotope = isotope

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
    def from_atom(cls, atom: Union['Element', 'DynamicElement']) -> 'DynamicElement':
        """
        get DynamicElement object from Element object or copy of DynamicElement object
        """
        if isinstance(atom, Element):
            return cls.from_atomic_number(atom.atomic_number)(atom.isotope)
        elif not isinstance(atom, DynamicElement):
            raise TypeError('Element or DynamicElement expected')
        return atom.copy()

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
