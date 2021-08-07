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


def _rules():
    from ...containers import QueryContainer

    rules = []

    # Aromatic N-Oxide
    #
    #  : N :  >>  : [N+] :
    #    \\           \
    #     O           [O-]
    #
    q = QueryContainer()
    q.add_atom('N', neighbors=3, hybridization=4)
    q.add_atom('O', neighbors=1)
    q.add_bond(1, 2, 2)
    atom_fix = {1: {'_charges': 1}, 2: {'_charges': -1, '_hybridizations': 1}}
    bonds_fix = ((1, 2, 1),)
    rules.append((q, atom_fix, bonds_fix))

    # Aromatic N-Nitride?
    #
    #  : N :  >>  : [N+] :
    #    \\           \
    #     N           [N-]
    #
    q = QueryContainer()
    q.add_atom('N', neighbors=3, hybridization=4)
    q.add_atom('N', neighbors=(1, 2), hybridization=2)
    q.add_bond(1, 2, 2)
    atom_fix = {1: {'_charges': 1}, 2: {'_charges': -1, '_hybridizations': 1}}
    bonds_fix = ((1, 2, 1),)
    rules.append((q, atom_fix, bonds_fix))

    #
    # : [S+] : >> : S :
    #    |          \\
    #   [O-]         O
    #
    q = QueryContainer()
    q.add_atom('S', neighbors=3, hybridization=4, charge=1)
    q.add_atom('O', neighbors=1, charge=-1)
    q.add_bond(1, 2, 1)
    atom_fix = {1: {'_charges': 0}, 2: {'_charges': 0, '_hybridizations': 2}}
    bonds_fix = ((1, 2, 2),)
    rules.append((q, atom_fix, bonds_fix))

    #
    # [O-]-N:C:C:[N+]=O
    #
    q = QueryContainer()
    q.add_atom('O', neighbors=1, charge=-1)
    q.add_atom('N', neighbors=3)
    q.add_atom('C')
    q.add_atom('C')
    q.add_atom('N', neighbors=3, charge=1)
    q.add_atom('O', neighbors=1)
    q.add_bond(1, 2, 1)
    q.add_bond(2, 3, 4)
    q.add_bond(3, 4, 4)
    q.add_bond(4, 5, 4)
    q.add_bond(5, 6, 2)
    atom_fix = {2: {'_charges': 1}, 6: {'_charges': -1}}
    bonds_fix = ((5, 6, 1),)
    rules.append((q, atom_fix, bonds_fix))

    #
    # N : A : N - ?
    #  :     :
    #   C # C
    q = QueryContainer()
    q.add_atom('N', neighbors=2)
    q.add_atom('C', neighbors=2)
    q.add_atom('C', neighbors=2)
    q.add_atom('N', neighbors=(2, 3))
    q.add_atom(ListElement(['C', 'N']))
    q.add_bond(1, 2, 4)
    q.add_bond(2, 3, 3)
    q.add_bond(3, 4, 4)
    q.add_bond(4, 5, 4)
    q.add_bond(1, 5, 4)
    atom_fix = {}
    bonds_fix = ((2, 3, 4),)
    rules.append((q, atom_fix, bonds_fix))

    #
    # C:[N+]:[C-]
    #    \\
    #     O
    #
    q = QueryContainer()
    q.add_atom('N', neighbors=3, charge=1)
    q.add_atom('O', neighbors=1)
    q.add_atom('C', neighbors=(2, 3), charge=-1)
    q.add_atom('C', neighbors=(2, 3))
    q.add_bond(1, 2, 2)
    q.add_bond(1, 3, 4)
    q.add_bond(1, 4, 4)
    atom_fix = {2: {'_charges': -1, '_hybridizations': 1}, 3: {'_charges': 0}}
    bonds_fix = ((1, 2, 1),)
    rules.append((q, atom_fix, bonds_fix))

    #
    #  O=[N+] : C
    #     :     :
    #    O : N : C
    q = QueryContainer()
    q.add_atom('N', neighbors=3, charge=1)
    q.add_atom('O', neighbors=1)
    q.add_atom('O', neighbors=2)
    q.add_atom('C', neighbors=(2, 3))
    q.add_atom('C', neighbors=(2, 3))
    q.add_atom('N', neighbors=(2, 3))
    q.add_bond(1, 2, 2)
    q.add_bond(1, 3, 4)
    q.add_bond(1, 4, 4)
    q.add_bond(3, 6, 4)
    q.add_bond(4, 5, 4)
    q.add_bond(5, 6, 4)
    atom_fix = {1: {'_hybridizations': 2}, 3: {'_hybridizations': 1}, 4: {'_hybridizations': 2},
                5: {'_hybridizations': 2}, 6: {'_hybridizations': 1}}
    bonds_fix = ((1, 3, 1), (1, 4, 1), (3, 6, 1), (4, 5, 2), (5, 6, 1))
    rules.append((q, atom_fix, bonds_fix))
    return rules


rules = Proxy(_rules)


__all__ = ['rules']
