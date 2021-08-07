# -*- coding: utf-8 -*-
#
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
from ...periodictable import ListElement


def _stripped_rules():
    from ...containers import QueryContainer
    rules = []

    # Ammonia, H-imine,guanidine,amidine. [H][N+]=,-,:
    q = QueryContainer()
    q.add_atom('N', charge=1, hydrogens=(1, 2, 3, 4))
    rules.append(q)
    return rules


def _rules():
    from ...containers import QueryContainer
    rules = _stripped_rules()

    # Phenoles, [H][O,S,Se]-Ar
    q = QueryContainer()
    q.add_atom(ListElement(['O', 'S', 'Se']), neighbors=1)
    q.add_atom(ListElement(['C', 'N']), hybridization=4)
    q.add_bond(1, 2, 1)
    rules.append(q)

    # Oxo-acids. [H][O,S,Se]-[C,N,P,S,Cl,Se,Br,I]=O
    q = QueryContainer()
    q.add_atom(ListElement(['O', 'S', 'Se']), neighbors=1)
    q.add_atom(ListElement(['C', 'N', 'P', 'S', 'Cl', 'Se', 'Br', 'I']))
    q.add_atom('O')
    q.add_bond(1, 2, 1)
    q.add_bond(2, 3, 2)
    rules.append(q)

    # Nitro acid. [H]O-[N+](=O)[O-]
    q = QueryContainer()
    q.add_atom('O', neighbors=1)
    q.add_atom('N', charge=1)
    q.add_atom('O', charge=-1)
    q.add_atom('O')
    q.add_bond(1, 2, 1)
    q.add_bond(2, 3, 1)
    q.add_bond(2, 4, 2)
    rules.append(q)

    # Halogen acids
    q = QueryContainer()
    q.add_atom(ListElement(['F', 'Cl', 'Br', 'I']), neighbors=0)
    rules.append(q)
    return rules


stripped_rules = Proxy(_stripped_rules)
rules = Proxy(_rules)


__all__ = ['stripped_rules', 'rules']
