# -*- coding: utf-8 -*-
#
#  Copyright 2019-2022 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from weakref import ref
from ..exceptions import IsConnectedBond, IsNotConnectedBond


class Bond:
    __slots__ = ('__order', '__graph', '__n', '__m')

    def __init__(self, order: int):
        if not isinstance(order, int):
            raise TypeError('invalid order value')
        elif order not in (1, 4, 2, 3, 8):
            raise ValueError('order should be from [1, 2, 3, 4, 8]')
        self.__order = order

    def __eq__(self, other):
        if isinstance(other, Bond):
            return self.__order == other.order
        elif isinstance(other, int):
            return self.__order == other
        return False

    def __repr__(self):
        return f'{self.__class__.__name__}({self.__order})'

    def __int__(self):
        """
        Bond order.
        """
        return self.__order

    def __hash__(self):
        """
        Bond order. Used in Morgan atoms ordering.
        """
        return self.__order

    def __getstate__(self):
        return {'order': self.__order}

    def __setstate__(self, state):
        self.__order = state['order']

    @property
    def order(self) -> int:
        return self.__order

    @property
    def in_ring(self) -> bool:
        try:
            return self.__graph().is_ring_bond(self.__n, self.__m)
        except AttributeError:
            raise IsNotConnectedBond

    def copy(self) -> 'Bond':
        copy = object.__new__(self.__class__)
        copy._Bond__order = self.__order
        return copy

    @classmethod
    def from_bond(cls, bond):
        if isinstance(bond, cls):
            copy = object.__new__(cls)
            copy._Bond__order = bond.order
            return copy
        raise TypeError('Bond expected')

    def _attach_graph(self, graph, n, m):
        try:
            self.__graph
        except AttributeError:
            self.__graph = ref(graph)
            self.__n = n
            self.__m = m
        else:
            raise IsConnectedBond

    def _change_map(self, n, m):
        try:
            self.__graph
        except AttributeError:
            raise IsNotConnectedBond
        else:
            self.__n = n
            self.__m = m


class DynamicBond:
    __slots__ = ('__order', '__p_order')

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

        self.__order = order
        self.__p_order = p_order

    def __eq__(self, other):
        if isinstance(other, DynamicBond):
            return self.__order == other.order and self.__p_order == other.p_order
        return False

    def __repr__(self):
        return f'{self.__class__.__name__}({self.__order}, {self.__p_order})'

    def __int__(self):
        """
        Hash of bond orders.
        """
        return hash(self)

    def __hash__(self):
        """
        Hash of bond orders.
        """
        return hash((self.__order or 0, self.__p_order or 0))

    @property
    def is_dynamic(self) -> bool:
        """
        Bond has dynamic features
        """
        return self.__order != self.__p_order

    @property
    def order(self) -> Optional[int]:
        return self.__order

    @property
    def p_order(self) -> Optional[int]:
        return self.__p_order

    def copy(self) -> 'DynamicBond':
        copy = object.__new__(self.__class__)
        copy._DynamicBond__order = self.__order
        copy._DynamicBond__p_order = self.__p_order
        return copy

    @classmethod
    def from_bond(cls, bond):
        if isinstance(bond, Bond):
            copy = object.__new__(cls)
            copy._DynamicBond__order = copy._DynamicBond__p_order = bond.order
            return copy
        elif isinstance(bond, cls):
            copy = object.__new__(cls)
            copy._DynamicBond__order = bond.order
            copy._DynamicBond__p_order = bond.p_order
            return copy
        raise TypeError('DynamicBond expected')


class QueryBond:
    __slots__ = ('__order', '__in_ring')

    def __init__(self, order: Union[int, List[int], Set[int], Tuple[int, ...]], in_ring: Optional[bool] = None):
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
        self.__order = order
        self.__in_ring = in_ring

    def __eq__(self, other):
        if isinstance(other, Bond):
            if self.__in_ring is not None:
                if self.__in_ring != other.in_ring:
                    return False
            return other.order in self.__order
        elif isinstance(other, QueryBond):
            return self.__order == other.order and self.__in_ring == other.in_ring
        elif isinstance(other, int):
            return other in self.__order
        return False

    def __repr__(self):
        return f'{self.__class__.__name__}({self.__order}, {self.__in_ring})'

    def __int__(self):
        """
        Simple bond order or hash of sorted tuple of orders.
        """
        if len(self.__order) == 1:
            return self.__order[0]
        return hash(self.__order)

    def __hash__(self):
        """
        Hash of orders and cycle mark. Used in Morgan atoms ordering.
        """
        return hash((self.__order, self.__in_ring))

    @property
    def order(self) -> Tuple[int, ...]:
        return self.__order

    @property
    def in_ring(self) -> Optional[bool]:
        return self.__in_ring

    def copy(self) -> 'QueryBond':
        copy = object.__new__(self.__class__)
        copy._QueryBond__order = self.__order
        copy._QueryBond__in_ring = self.__in_ring
        return copy

    @classmethod
    def from_bond(cls, bond):
        if isinstance(bond, Bond):
            copy = object.__new__(cls)
            copy._QueryBond__order = (bond.order,)
            copy._QueryBond__in_ring = None
            return copy
        elif isinstance(bond, cls):
            copy = object.__new__(cls)
            copy._QueryBond__order = bond.order
            copy._QueryBond__in_ring = bond.in_ring
            return copy
        raise TypeError('QueryBond or Bond expected')


__all__ = ['Bond', 'DynamicBond', 'QueryBond']
