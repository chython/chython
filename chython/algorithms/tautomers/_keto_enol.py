# -*- coding: utf-8 -*-
#
#  Copyright 2022 Ramil Nugmanov <nougmanoff@protonmail.com>
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


def _sugar_group():
    from ...containers import QueryContainer
    q = QueryContainer()
    q.add_atom(ListElement(['O', 'N']), hydrogens=(1, 2))  # enol
    q.add_atom(ListElement(['O', 'N']))  # ketone
    q.add_atom('C', hybridization=1)
    q.add_atom('C', hybridization=2)
    q.add_bond(1, 3, 1)
    q.add_bond(3, 4, 1)
    q.add_bond(2, 4, 2)
    return q


def _keto_rules():
    from ...containers import QueryContainer
    rules = []
    # first atom is H-acceptor
    # second is direction

    # C-C=[O,S,NH]
    q = QueryContainer()
    q.add_atom(ListElement(['O', 'S', 'N']), neighbors=1)
    q.add_atom('C', heteroatoms=1, hybridization=2, neighbors=(2, 3))
    q.add_bond(1, 2, 2)
    rules.append(q)

    # C-C=N-[C,N,O]
    q = QueryContainer()
    q.add_atom('N', neighbors=2)
    q.add_atom('C', heteroatoms=1, hybridization=2, neighbors=(2, 3))
    q.add_atom(ListElement(['C', 'N', 'O']))
    q.add_bond(1, 2, 2)
    q.add_bond(1, 3, 1)
    rules.append(q)

    # [C,H]-N=N-C
    q = QueryContainer()
    q.add_atom('N', neighbors=(1, 2), heteroatoms=1)
    q.add_atom('N')
    q.add_atom('C', heteroatoms=1)
    q.add_bond(1, 2, 2)
    q.add_bond(2, 3, 1)
    rules.append(q)

    # [S,O,NR;H]-C=N
    q = QueryContainer()
    q.add_atom('N')
    q.add_atom('C')
    q.add_atom(ListElement(['N', 'S', 'O']), hydrogens=1)
    q.add_bond(1, 2, 2)
    q.add_bond(2, 3, 1)
    rules.append(q)

    # [NH2]-C=N-R
    q = QueryContainer()
    q.add_atom('N', neighbors=2)
    q.add_atom('C')
    q.add_atom('N', neighbors=1)
    q.add_bond(1, 2, 2)
    q.add_bond(2, 3, 1)
    rules.append(q)

    # S=C-[N,O;H]
    q = QueryContainer()
    q.add_atom('S', neighbors=1)
    q.add_atom('C')
    q.add_atom(ListElement(['N', 'O']), hydrogens=(1, 2))
    q.add_bond(1, 2, 2)
    q.add_bond(2, 3, 1)
    rules.append(q)

    # O=C-[N,S;H]
    q = QueryContainer()
    q.add_atom('O')
    q.add_atom('C')
    q.add_atom(ListElement(['N', 'S']), hydrogens=(1, 2))
    q.add_bond(1, 2, 2)
    q.add_bond(2, 3, 1)
    rules.append(q)
    return rules


def _enol_rules():
    from ...containers import QueryContainer
    rules = []
    # first atom is H-donor
    # second is direction

    # C=C-[OH,SH,NH2]
    q = QueryContainer()
    q.add_atom(ListElement(['O', 'S', 'N']), neighbors=1)
    q.add_atom('C', heteroatoms=1, hybridization=2, neighbors=(2, 3))
    q.add_bond(1, 2, 1)
    rules.append(q)

    # C=C-[NH]-[C,N,O]
    q = QueryContainer()
    q.add_atom('N', neighbors=2)
    q.add_atom('C', heteroatoms=1, hybridization=2, neighbors=(2, 3))
    q.add_atom(ListElement(['C', 'N', 'O']))
    q.add_bond(1, 2, 1)
    q.add_bond(1, 3, 1)
    rules.append(q)

    # [C,H]-[NH]-N=C
    q = QueryContainer()
    q.add_atom('N', neighbors=(1, 2), hybridization=1, heteroatoms=1)
    q.add_atom('N')
    q.add_atom('C', heteroatoms=1)
    q.add_bond(1, 2, 1)
    q.add_bond(2, 3, 2)
    rules.append(q)

    # [S,O,NR;H]-C=N
    q = QueryContainer()
    q.add_atom(ListElement(['N', 'S', 'O']), hydrogens=1)
    q.add_atom('C')
    q.add_atom('N')
    q.add_bond(1, 2, 1)
    q.add_bond(2, 3, 2)
    rules.append(q)

    # [NH2]-C=N-R
    q = QueryContainer()
    q.add_atom('N', neighbors=1)
    q.add_atom('C')
    q.add_atom('N', neighbors=2)
    q.add_bond(1, 2, 1)
    q.add_bond(2, 3, 2)
    rules.append(q)

    # S=C-[N,O;H]
    q = QueryContainer()
    q.add_atom(ListElement(['N', 'O']), hydrogens=(1, 2))
    q.add_atom('C')
    q.add_atom('S', neighbors=1)
    q.add_bond(1, 2, 1)
    q.add_bond(2, 3, 2)
    rules.append(q)

    # O=C-[N,S;H]
    q = QueryContainer()
    q.add_atom(ListElement(['N', 'S']), hydrogens=(1, 2))
    q.add_atom('C')
    q.add_atom('O')
    q.add_bond(1, 2, 1)
    q.add_bond(2, 3, 2)
    rules.append(q)
    return rules


keto_rules = Proxy(_keto_rules)
enol_rules = Proxy(_enol_rules)
sugar_group = Proxy(_sugar_group)


__all__ = ['keto_rules', 'enol_rules', 'sugar_group']
