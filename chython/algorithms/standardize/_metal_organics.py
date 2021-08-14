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
from ...periodictable import ListElement, AnyMetal


def _rules():
    from ...containers import QueryContainer
    rules = []

    #
    #        R - N - C                 R - N - C
    #           /    ||                   /    ||
    #  A - M = C     || >>    A - M .. [C-]    ||
    #           \    ||                  \\    ||
    #        R - N - C               R - [N+]- C
    #
    q = QueryContainer()
    q.add_atom(AnyMetal())
    q.add_atom('C')
    q.add_atom('N', hybridization=1, neighbors=3, heteroatoms=0)
    q.add_atom('N', hybridization=1, neighbors=3, heteroatoms=0)
    q.add_atom('C', hybridization=2)
    q.add_atom('C', hybridization=2)
    q.add_bond(1, 2, 2)
    q.add_bond(2, 3, 1)
    q.add_bond(2, 4, 1)
    q.add_bond(3, 5, 1)
    q.add_bond(4, 6, 1)
    q.add_bond(5, 6, (1, 2))

    atom_fix = {2: (-1, None), 3: (1, None)}  # atom: (charge diff, new radical state or None)
    bonds_fix = ((1, 2, 8), (2, 3, 2))
    rules.append((q, atom_fix, bonds_fix))

    #
    # cyanide
    #
    q = QueryContainer()
    q.add_atom('C')
    q.add_atom('N')
    q.add_atom(AnyMetal())
    q.add_bond(1, 2, 3)
    q.add_bond(1, 3, 1)

    # NOTE: first fixed atom is metal. required for preventing errors in charge.
    atom_fix = {3: (1, None), 1: (-1, None)}
    bonds_fix = ((1, 3, 8),)
    rules.append((q, atom_fix, bonds_fix))

    q = QueryContainer()
    q.add_atom('C', charge=-1)
    q.add_atom('N')
    q.add_atom(AnyMetal())
    q.add_bond(1, 2, 3)
    q.add_bond(1, 3, 1)

    atom_fix = {}
    bonds_fix = ((1, 3, 8),)
    rules.append((q, atom_fix, bonds_fix))

    #
    # carbonyl
    #
    q = QueryContainer()
    q.add_atom('C')
    q.add_atom('O')
    q.add_atom(AnyMetal())
    q.add_bond(1, 2, 3)
    q.add_bond(1, 3, 1)

    atom_fix = {1: (-1, None), 2: (1, None)}
    bonds_fix = ((1, 3, 8),)
    rules.append((q, atom_fix, bonds_fix))

    q = QueryContainer()
    q.add_atom('C', is_radical=True)
    q.add_atom('O')
    q.add_atom(AnyMetal())
    q.add_bond(1, 2, 2)
    q.add_bond(1, 3, 1)

    atom_fix = {1: (-1, False), 2: (1, None)}
    bonds_fix = ((1, 3, 8), (1, 2, 3))
    rules.append((q, atom_fix, bonds_fix))

    q = QueryContainer()
    q.add_atom('C', neighbors=2)
    q.add_atom('O')
    q.add_atom(AnyMetal())
    q.add_bond(1, 2, 2)
    q.add_bond(1, 3, 1)

    atom_fix = {1: (-1, None), 2: (1, None)}
    bonds_fix = ((1, 3, 8), (1, 2, 3))
    rules.append((q, atom_fix, bonds_fix))

    q = QueryContainer()
    q.add_atom('C')
    q.add_atom('O')
    q.add_atom(AnyMetal())
    q.add_atom(AnyMetal())
    q.add_bond(1, 2, 2)
    q.add_bond(1, 3, 1)
    q.add_bond(1, 4, 1)

    atom_fix = {1: (-1, None), 2: (1, None)}
    bonds_fix = ((1, 3, 8), (1, 4, 8), (1, 2, 3))
    rules.append((q, atom_fix, bonds_fix))

    #
    # Ferrocene covalent uncharged
    #
    q = QueryContainer()
    q.add_atom(AnyMetal())
    q.add_atom('C')
    q.add_atom('C')
    q.add_atom('C')
    q.add_atom('C')
    q.add_atom('C')
    q.add_bond(1, 2, 1); q.add_bond(1, 3, 1); q.add_bond(1, 4, 1); q.add_bond(1, 5, 1); q.add_bond(1, 6, 1)
    q.add_bond(2, 3, 1); q.add_bond(3, 4, 1); q.add_bond(4, 5, 1); q.add_bond(5, 6, 1); q.add_bond(6, 2, 1)

    atom_fix = {1: (1, None), 2: (-1, None)}
    bonds_fix = ((1, 2, 8), (1, 3, 8), (1, 4, 8), (1, 5, 8), (1, 6, 8), (3, 4, 2), (5, 6, 2))
    rules.append((q, atom_fix, bonds_fix))

    #
    # Ferrocene covalent uncharged. invalid valence.
    #
    q = QueryContainer()
    q.add_atom(AnyMetal())
    q.add_atom('C')
    q.add_atom('C')
    q.add_atom('C')
    q.add_atom('C')
    q.add_atom('C')
    q.add_bond(1, 2, 1); q.add_bond(1, 3, 1); q.add_bond(1, 4, 1); q.add_bond(1, 5, 1); q.add_bond(1, 6, 1)
    q.add_bond(2, 3, 1); q.add_bond(3, 4, 2); q.add_bond(4, 5, 1); q.add_bond(5, 6, 2); q.add_bond(6, 2, 1)

    atom_fix = {1: (1, None), 2: (-1, None)}
    bonds_fix = ((1, 2, 8), (1, 3, 8), (1, 4, 8), (1, 5, 8), (1, 6, 8), (3, 4, 2), (5, 6, 2))
    rules.append((q, atom_fix, bonds_fix))

    #
    # Ferrocene covalent radical carbon
    #
    q = QueryContainer()
    q.add_atom(AnyMetal())
    q.add_atom('C', is_radical=True)
    q.add_atom('C', is_radical=True)
    q.add_atom('C', is_radical=True)
    q.add_atom('C', is_radical=True)
    q.add_atom('C', is_radical=True)
    q.add_bond(1, 2, 1); q.add_bond(1, 3, 1); q.add_bond(1, 4, 1); q.add_bond(1, 5, 1); q.add_bond(1, 6, 1)
    q.add_bond(2, 3, 1); q.add_bond(3, 4, 1); q.add_bond(4, 5, 1); q.add_bond(5, 6, 1); q.add_bond(6, 2, 1)

    atom_fix = {1: (1, None), 2: (-1, False), 3: (0, False), 4: (0, False), 5: (0, False), 6: (0, False)}
    bonds_fix = ((1, 2, 8), (1, 3, 8), (1, 4, 8), (1, 5, 8), (1, 6, 8), (3, 4, 2), (5, 6, 2))
    rules.append((q, atom_fix, bonds_fix))

    #
    # Ferrocene covalent charged
    #
    q = QueryContainer()
    q.add_atom(AnyMetal())
    q.add_atom('C', charge=-1)
    q.add_atom('C')
    q.add_atom('C')
    q.add_atom('C')
    q.add_atom('C')
    q.add_bond(1, 2, 1); q.add_bond(1, 3, 1); q.add_bond(1, 4, 1); q.add_bond(1, 5, 1); q.add_bond(1, 6, 1)
    q.add_bond(2, 3, 1); q.add_bond(3, 4, 2); q.add_bond(4, 5, 1); q.add_bond(5, 6, 2); q.add_bond(6, 2, 1)

    atom_fix = {}
    bonds_fix = ((1, 2, 8), (1, 3, 8), (1, 4, 8), (1, 5, 8), (1, 6, 8))
    rules.append((q, atom_fix, bonds_fix))

    #
    # Ferrocene coordinate radical carbon
    #
    q = QueryContainer()
    q.add_atom(AnyMetal())
    q.add_atom('C', is_radical=True)
    q.add_atom('C', is_radical=True)
    q.add_atom('C', is_radical=True)
    q.add_atom('C', is_radical=True)
    q.add_atom('C', is_radical=True)
    q.add_bond(1, 2, 8); q.add_bond(1, 3, 8); q.add_bond(1, 4, 8); q.add_bond(1, 5, 8); q.add_bond(1, 6, 8)
    q.add_bond(2, 3, 1); q.add_bond(3, 4, 1); q.add_bond(4, 5, 1); q.add_bond(5, 6, 1); q.add_bond(6, 2, 1)

    atom_fix = {1: (1, None), 2: (-1, False), 3: (0, False), 4: (0, False), 5: (0, False), 6: (0, False)}
    bonds_fix = ((1, 2, 8), (1, 3, 8), (1, 4, 8), (1, 5, 8), (1, 6, 8), (3, 4, 2), (5, 6, 2))
    rules.append((q, atom_fix, bonds_fix))

    #
    # Allyl complexes
    #
    q = QueryContainer()
    q.add_atom(AnyMetal())
    q.add_atom('C')
    q.add_atom('C')
    q.add_atom('C', is_radical=True)
    q.add_bond(1, 2, 8); q.add_bond(1, 3, 8); q.add_bond(1, 4, 8)
    q.add_bond(2, 3, 2)
    q.add_bond(3, 4, 1)

    atom_fix = {1: (1, None), 4: (-1, False)}
    bonds_fix = ()
    rules.append((q, atom_fix, bonds_fix))
    return rules


rules = Proxy(_rules)


__all__ = ['rules']
