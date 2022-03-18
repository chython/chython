# -*- coding: utf-8 -*-
#
#  Copyright 2021, 2022 Ramil Nugmanov <nougmanoff@protonmail.com>
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


def _xonyl_groups():
    from ...containers import QueryContainer
    rules = []

    q = QueryContainer()
    q.add_atom('N', neighbors=3, hybridization=2, charge=1)
    q.add_atom('O', neighbors=1, charge=-1)
    q.add_atom('O', neighbors=1)
    q.add_bond(1, 2, 1)
    q.add_bond(1, 3, 2)
    rules.append(q)

    q = QueryContainer()
    q.add_atom(ListElement(['C', 'P', 'S', 'Cl', 'As', 'Se', 'Br', 'Sb', 'Te', 'I', 'Bi', 'Po', 'At']))
    q.add_atom(ListElement(['S', 'O', 'Se', 'N']), neighbors=1)
    q.add_atom(ListElement(['S', 'O', 'Se', 'N']), neighbors=1)
    q.add_bond(1, 2, 1)
    q.add_bond(1, 3, 2)
    rules.append(q)

    q = QueryContainer()
    q.add_atom(ListElement(['C', 'P', 'S', 'Cl', 'As', 'Se', 'Br', 'Sb', 'Te', 'I', 'Bi', 'Po', 'At']))
    q.add_atom(ListElement(['S', 'O', 'Se', 'N']), neighbors=1, charge=-1)
    q.add_atom(ListElement(['S', 'O', 'Se', 'N']), neighbors=1)
    q.add_bond(1, 2, 1)
    q.add_bond(1, 3, 2)
    rules.append(q)
    return rules


def _substituents():
    """
    Rules for switchable functional groups remapping
    """
    from ...containers import QueryContainer
    rules = []

    # nitro addition
    #
    #      O [3]     [2] O -- *
    #     //            /:
    # * - N+   >>  * = N+
    #      \            \
    #       O- [2]       O- [3]
    #
    q = QueryContainer()
    q.add_atom('N', neighbors=3, charge=1, hybridization=2)
    q.add_atom('O', neighbors=1, charge=-1)
    q.add_atom('O', neighbors=2)
    q.add_bond(1, 2, 1)
    q.add_bond(1, 3, (1, 4))
    rules.append((q, ((2, 3), (3, 2))))

    #
    #      NH [3]          NH2 [3]
    #     //               /
    # * - X       >>  * - X
    #      \              :
    #       NH2 [2]       N [2]
    #
    q = QueryContainer()
    q.add_atom('C', neighbors=3)
    q.add_atom('N')
    q.add_atom('N', neighbors=2, hydrogens=0)
    q.add_bond(1, 2, 1)
    q.add_bond(1, 3, (2, 4))
    rules.append((q, ((2, 3), (3, 2))))

    #
    #      O [3]        [3] O - R
    #     //               /
    # * - X       >>  * - X
    #      \              \\
    #       OH [2]         O [2]
    #
    q = QueryContainer()
    q.add_atom(ListElement(['C', 'P', 'S', 'Cl', 'As', 'Se', 'Br', 'Sb', 'Te', 'I', 'Bi', 'Po', 'At']))
    q.add_atom(ListElement(['S', 'O', 'Se', 'N']), neighbors=2)
    q.add_atom(ListElement(['S', 'O', 'Se', 'N']), neighbors=1)
    q.add_bond(1, 2, 1)
    q.add_bond(1, 3, 2)
    rules.append((q, ((2, 3), (3, 2))))

    #
    #      O [3]            A
    #     //               /
    # * - X       >>  * - X
    #      \              \\
    #       OH [2]         O [2]
    #
    q = QueryContainer()
    q.add_atom(ListElement(['C', 'P', 'S', 'Cl', 'As', 'Se', 'Br', 'Sb', 'Te', 'I', 'Bi', 'Po', 'At']))
    q.add_atom('A')
    q.add_atom(ListElement(['S', 'O', 'Se', 'N']), neighbors=1)
    q.add_bond(1, 2, 1)
    q.add_bond(1, 3, 2)
    rules.append((q, ((3, 2),)))  # possible only: (3, 2)
    return rules


xonyl_groups = Proxy(_xonyl_groups)
substituents_groups = Proxy(_substituents)


__all__ = ['xonyl_groups', 'substituents_groups']
