# -*- coding: utf-8 -*-
#
#  Copyright 2018-2025 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from typing import Tuple, Union
from .bonds import Bond, QueryBond
from .graph import Graph
from ..algorithms.isomorphism import QueryIsomorphism
from ..periodictable import Element, QueryElement
from ..periodictable.base import Query


class QueryContainer(Graph[Query, QueryBond], QueryIsomorphism):
    __slots__ = ('_smarts',)

    def __init__(self, smarts: str):
        super().__init__()
        self._smarts = smarts

    def __str__(self):
        return self._smarts

    def __repr__(self):
        return f'smarts({self._smarts})'

    def add_atom(self, atom: Union[Query, Element, int, str], *args, **kwargs):
        if not isinstance(atom, Query):
            # set only basic labels: charge, radical, isotope. use Query object directly for the full control.
            if isinstance(atom, Element):
                atom = QueryElement.from_atom(atom)
            elif isinstance(atom, str):
                atom = QueryElement.from_symbol(atom)()
            elif isinstance(atom, int):
                atom = QueryElement.from_atomic_number(atom)()
            else:
                raise TypeError('QueryElement object expected')
        return super().add_atom(atom, *args, **kwargs)

    def add_bond(self, n, m, bond: Union[QueryBond, Bond, int, Tuple[int, ...]]):
        if isinstance(bond, Bond):
            bond = QueryBond.from_bond(bond)
        elif not isinstance(bond, QueryBond):
            bond = QueryBond(bond)
        super().add_bond(n, m, bond)

    def union(self, other: 'QueryContainer', *, remap: bool = False, copy: bool = True) -> 'QueryContainer':
        if not isinstance(other, QueryContainer):
            raise TypeError('QueryContainer expected')
        return super().union(other, remap=remap, copy=copy)


__all__ = ['QueryContainer']
