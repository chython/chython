# -*- coding: utf-8 -*-
#
#  Copyright 2021 Aleksandr Sizov <murkyrussian@gmail.com>
#  Copyright 2021 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from lazy_object_proxy import Proxy
from ....containers import QueryContainer

_queries, _banned = [], []

q = QueryContainer()
q.add_atom('N', charge=1)
_queries.append(q)

q = QueryContainer()
q.add_atom('N', charge=0)
q.add_atom('A', hybridization=4)
q.add_bond(1, 2, 1)
_queries.append(q)

q = QueryContainer()
q.add_atom('N', charge=0, hybridization=1, heteroatoms=0)
q.add_atom('C', hybridization=(1, 2, 3))
q.add_bond(1, 2, 1)
_queries.append(q)

q = QueryContainer()
q.add_atom('N', charge=0)
q.add_atom('C')
q.add_atom('O')
q.add_bond(1, 2, 1)
q.add_bond(2, 3, 2)
_banned.append(q)

queries = Proxy(_queries)
banned = Proxy(_banned)


__all__ = ['queries', 'banned']
