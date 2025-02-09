# -*- coding: utf-8 -*-
#
#  Copyright 2019-2025 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from typing import Optional, Tuple, Union, List, Set


class Bond:
    __slots__ = ('_order', '_in_ring', '_stereo')

    def __init__(self, order: int, *, stereo: Optional[bool] = None):
        if not isinstance(order, int):
            raise TypeError('invalid order value')
        elif order not in (1, 4, 2, 3, 8):
            raise ValueError('order should be from [1, 2, 3, 4, 8]')
        self._order = order
        self._stereo = stereo

    def __eq__(self, other):
        if isinstance(other, int):
            return self.order == other
        elif isinstance(other, Bond):
            return self.order == other.order
        return False

    def __repr__(self):
        return f'{self.__class__.__name__}({self.order})'

    def __int__(self):
        """
        Bond order.
        """
        return self.order

    def __hash__(self):
        """
        Bond order. Used in Morgan atoms ordering.
        """
        return self.order

    @property
    def order(self) -> int:
        return self._order

    @property
    def stereo(self) -> Optional[bool]:
        return self._stereo

    @property
    def in_ring(self) -> bool:
        return self._in_ring

    def copy(self, full=False, stereo=False) -> 'Bond':
        copy = object.__new__(self.__class__)
        copy._order = self.order
        if full:
            copy._stereo = self.stereo
            copy._in_ring = self.in_ring
        else:
            if stereo:
                copy._stereo = self.stereo
            else:
                copy._stereo = None
        return copy

    def __copy__(self):
        return self.copy()


class DynamicBond:
    __slots__ = ('_order', '_p_order')

    def __init__(self, order=None, p_order=None):
        if order is None:
            if not isinstance(p_order, int):
                raise TypeError('p_order should be int type')
        elif not isinstance(order, int):
            raise TypeError('order should be int type or None')
        elif p_order is not None and not isinstance(p_order, int):
            raise TypeError('p_order should be int type or None')

        if order not in (1, 4, 2, 3, None, 8) or p_order not in (1, 4, 2, 3, None, 8):
            raise ValueError('order or p_order should be from [1, 2, 3, 4, 8]')

        self._order = order
        self._p_order = p_order

    def __eq__(self, other):
        if isinstance(other, DynamicBond):
            return self.order == other.order and self.p_order == other.p_order
        return False

    def __repr__(self):
        return f'{self.__class__.__name__}({self.order}, {self.p_order})'

    def __int__(self):
        """
        Hash of bond orders.
        """
        return hash(self)

    def __hash__(self):
        """
        Hash of bond orders.
        """
        return hash((self.order or 0, self.p_order or 0))

    @property
    def is_dynamic(self) -> bool:
        """
        Bond has dynamic features
        """
        return self.order != self.p_order

    @property
    def order(self) -> Optional[int]:
        return self._order

    @property
    def p_order(self) -> Optional[int]:
        return self._p_order

    def copy(self) -> 'DynamicBond':
        copy = object.__new__(self.__class__)
        copy._order = self.order
        copy._p_order = self.p_order
        return copy

    def __copy__(self):
        return self.copy()

    @classmethod
    def from_bond(cls, bond: 'Bond') -> 'DynamicBond':
        if not isinstance(bond, Bond):
            raise TypeError('Bond expected')
        copy = object.__new__(cls)
        copy._order = copy._p_order = bond.order
        return copy


class QueryBond:
    __slots__ = ('_order', '_in_ring', '_stereo')

    def __init__(self, order: Union[int, List[int], Set[int], Tuple[int, ...]],
                 in_ring: Optional[bool] = None, stereo: Optional[bool] = None):
        if isinstance(order, (list, tuple, set)):
            if not all(isinstance(x, int) for x in order):
                raise TypeError('invalid order value')
            if any(x not in (1, 4, 2, 3, 8) for x in order):
                raise ValueError('order should be from [1, 2, 3, 4, 8]')
            order = tuple(sorted(set(order)))
        elif isinstance(order, int):
            if order not in (1, 4, 2, 3, 8):
                raise ValueError('order should be from [1, 2, 3, 4, 8]')
            order = (order,)
        else:
            raise TypeError('invalid order value')
        if in_ring is not None and not isinstance(in_ring, bool):
            raise TypeError('in_ring mark should be boolean or None')
        self._order = order
        self._in_ring = in_ring
        self.stereo = stereo

    def __eq__(self, other):
        if isinstance(other, Bond):
            if self.in_ring is not None:
                if self.in_ring != other.in_ring:
                    return False
            return other.order in self.order
        elif isinstance(other, QueryBond):
            return self.order == other.order and self.in_ring == other.in_ring
        elif isinstance(other, int):
            return other in self.order
        return False

    def __repr__(self):
        return f'{self.__class__.__name__}({self.order}, {self.in_ring})'

    def __int__(self):
        """
        Simple bond order or hash of sorted tuple of orders.
        """
        if len(self.order) == 1:
            return self.order[0]
        return hash(self.order)

    def __hash__(self):
        """
        Hash of orders and cycle mark. Used in Morgan atoms ordering.
        """
        return hash((self.order, self.in_ring))

    @property
    def order(self) -> Tuple[int, ...]:
        return self._order

    @property
    def in_ring(self) -> Optional[bool]:
        return self._in_ring

    @property
    def stereo(self) -> Optional[bool]:
        return self._stereo

    @stereo.setter
    def stereo(self, value):
        if value is not None and not isinstance(value, bool):
            raise TypeError('stereo mark should be boolean or None')
        self._stereo = value

    def copy(self, full=False) -> 'QueryBond':
        copy = object.__new__(self.__class__)
        copy._order = self.order
        if full:
            copy._in_ring = self.in_ring
            copy._stereo = self.stereo
        else:
            copy._in_ring = copy._stereo = None
        return copy

    def __copy__(self):
        return self.copy()

    @classmethod
    def from_bond(cls, bond: 'Bond', stereo=False, in_ring=False) -> 'QueryBond':
        if not isinstance(bond, Bond):
            raise TypeError('Bond expected')
        copy = object.__new__(cls)
        copy._order = (bond.order,)
        if in_ring:
            copy._in_ring = bond.in_ring
        else:
            copy._in_ring = None
        if stereo:
            copy._stereo = bond.stereo
        else:
            copy._stereo = None
        return copy


__all__ = ['Bond', 'DynamicBond', 'QueryBond']
