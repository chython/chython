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


def _rules():
    """
    Rules for switchable functional groups remapping
    """
    from ...containers import QueryContainer
    raw_rules = []

    # First atoms can be list of atoms or any atom, but other atoms should be same for reactants and products parts.

    # nitro-group switch charge
    #
    #      O
    #     //
    # * - N+
    #      \
    #       O-
    #
    atoms = ({'atom': 'N', 'neighbors': 3, 'hybridization': 2, 'charge': 1},
             {'atom': 'O', 'neighbors': 1, 'charge': -1}, {'atom': 'O', 'neighbors': 1})
    bonds = ((1, 2, 1), (1, 3, 2))
    _map = {2: 3, 3: 2}
    fix = {2: 3, 3: 2}
    raw_rules.append(((atoms, bonds), (atoms, bonds), _map, fix))

    # keep Xonyl-group in product
    atoms = ({'atom': ListElement(['C', 'N', 'P', 'S', 'Cl', 'As', 'Se', 'Br', 'Sb', 'Te', 'I', 'Bi', 'Po', 'At'])},
             {'atom': 'A'}, {'atom': 'O', 'neighbors': 1})
    bonds = ((1, 2, 1), (1, 3, 2))
    p_atoms = ({'atom': ListElement(['C', 'N', 'P', 'S', 'Cl', 'As', 'Se', 'Br', 'Sb', 'Te', 'I', 'Bi', 'Po', 'At'])},
               {'atom': 'O', 'neighbors': 1}, {'atom': 'O', 'neighbors': 1})
    p_bonds = ((1, 2, 1), (1, 3, 2))
    _map = {3: 2}
    fix = {3: 2, 2: 3}
    raw_rules.append(((atoms, bonds), (p_atoms, p_bonds), _map, fix))

    # oxo-acid salt switch charge
    #
    #      O
    #     //
    # * - X
    #      \
    #       O-
    #
    atoms = ({'atom': ListElement(['C', 'N', 'P', 'S', 'Cl', 'As', 'Se', 'Br', 'Sb', 'Te', 'I', 'Bi', 'Po', 'At'])},
             {'atom': 'O', 'neighbors': 1, 'charge': -1}, {'atom': 'O', 'neighbors': 1})
    bonds = ((1, 2, 1), (1, 3, 2))
    p_atoms = ({'atom': ListElement(['C', 'N', 'P', 'S', 'Cl', 'As', 'Se', 'Br', 'Sb', 'Te', 'I', 'Bi', 'Po', 'At'])},
               {'atom': 'O', 'neighbors': 1, 'charge': -1}, {'atom': 'O', 'neighbors': 1})
    p_bonds = ((1, 2, 1), (1, 3, 2))
    _map = {2: 3, 3: 2}
    fix = {2: 3, 3: 2}
    raw_rules.append(((atoms, bonds), (p_atoms, p_bonds), _map, fix))

    # oxo-acid switch hydrogen
    #
    #      O
    #     //
    # * - X
    #      \
    #       OH
    #
    atoms = ({'atom': ListElement(['C', 'N', 'P', 'S', 'Cl', 'As', 'Se', 'Br', 'Sb', 'Te', 'I', 'Bi', 'Po', 'At'])},
             {'atom': 'O', 'neighbors': 1}, {'atom': 'O', 'neighbors': 1})
    bonds = ((1, 2, 1), (1, 3, 2))
    p_atoms = ({'atom': ListElement(['C', 'N', 'P', 'S', 'Cl', 'As', 'Se', 'Br', 'Sb', 'Te', 'I', 'Bi', 'Po', 'At'])},
               {'atom': 'O', 'neighbors': 1}, {'atom': 'O', 'neighbors': 1})
    p_bonds = ((1, 2, 1), (1, 3, 2))
    _map = {2: 3, 3: 2}
    fix = {2: 3, 3: 2}
    raw_rules.append(((atoms, bonds), (p_atoms, p_bonds), _map, fix))

    # oxo-acid to salt keep Xonyl-group
    atoms = ({'atom': ListElement(['C', 'N', 'P', 'S', 'Cl', 'As', 'Se', 'Br', 'Sb', 'Te', 'I', 'Bi', 'Po', 'At'])},
             {'atom': 'O', 'neighbors': 1}, {'atom': 'O', 'neighbors': 1})
    bonds = ((1, 2, 1), (1, 3, 2))
    p_atoms = ({'atom': ListElement(['C', 'N', 'P', 'S', 'Cl', 'As', 'Se', 'Br', 'Sb', 'Te', 'I', 'Bi', 'Po', 'At'])},
               {'atom': 'O', 'neighbors': 1, 'charge': -1}, {'atom': 'O', 'neighbors': 1})
    p_bonds = ((1, 2, 1), (1, 3, 2))
    _map = {2: 3, 3: 2}
    fix = {2: 3, 3: 2}
    raw_rules.append(((atoms, bonds), (p_atoms, p_bonds), _map, fix))

    # oxo-acid salt to acid keep Xonyl-group
    atoms = ({'atom': ListElement(['C', 'N', 'P', 'S', 'Cl', 'As', 'Se', 'Br', 'Sb', 'Te', 'I', 'Bi', 'Po', 'At'])},
             {'atom': 'O', 'neighbors': 1, 'charge': -1}, {'atom': 'O', 'neighbors': 1})
    bonds = ((1, 2, 1), (1, 3, 2))
    p_atoms = ({'atom': ListElement(['C', 'N', 'P', 'S', 'Cl', 'As', 'Se', 'Br', 'Sb', 'Te', 'I', 'Bi', 'Po', 'At'])},
               {'atom': 'O', 'neighbors': 1}, {'atom': 'O', 'neighbors': 1})
    p_bonds = ((1, 2, 1), (1, 3, 2))
    _map = {2: 3, 3: 2}
    fix = {2: 3, 3: 2}
    raw_rules.append(((atoms, bonds), (p_atoms, p_bonds), _map, fix))

    # nitro addition
    #
    #      O [3]     [3] O -- *
    #     //            /
    # * - N+   >>  * = N+
    #      \            \
    #       O- [2]       O- [2]
    #
    atoms = ({'atom': 'N', 'neighbors': 3, 'charge': 1, 'hybridization': 2},
             {'atom': 'O', 'neighbors': 1, 'charge': -1}, {'atom': 'O', 'neighbors': 1})
    bonds = ((1, 2, 1), (1, 3, 2))
    p_atoms = ({'atom': 'N', 'neighbors': 3, 'charge': 1, 'hybridization': 2},
               {'atom': 'O', 'neighbors': 1, 'charge': -1}, {'atom': 'O', 'neighbors': 2})
    p_bonds = ((1, 2, 1), (1, 3, 1))
    _map = {2: 3, 3: 2}
    fix = {2: 3, 3: 2}
    raw_rules.append(((atoms, bonds), (p_atoms, p_bonds), _map, fix))

    # sulphate addition
    #
    #      O [3]        [2] O -- *
    #     //               /
    # * = S - *   >>  * = S - *
    #     |               \\
    #     O- [2]           O [3]
    #
    atoms = ({'atom': 'S', 'neighbors': 4, 'hybridization': 3}, {'atom': 'O', 'neighbors': 1, 'charge': -1},
             {'atom': 'O', 'neighbors': 1})
    bonds = ((1, 2, 1), (1, 3, 2))
    p_atoms = ({'atom': 'S', 'neighbors': 4, 'hybridization': 3}, {'atom': 'O', 'neighbors': 2},
               {'atom': 'O', 'neighbors': 1})
    p_bonds = ((1, 2, 1), (1, 3, 2))
    _map = {2: 3, 3: 2}
    fix = {2: 3, 3: 2}
    raw_rules.append(((atoms, bonds), (p_atoms, p_bonds), _map, fix))

    compiled_rules = []
    for (r_atoms, r_bonds), (p_atoms, p_bonds), _map, fix in raw_rules:
        r_q = QueryContainer()
        p_q = QueryContainer()
        for a in r_atoms:
            r_q.add_atom(**a)
        for n, m, b in r_bonds:
            r_q.add_bond(n, m, b)
        for a in p_atoms:
            p_q.add_atom(**a)
        for n, m, b in p_bonds:
            p_q.add_bond(n, m, b)
        compiled_rules.append((r_q, p_q, _map, fix))
    return compiled_rules


rules = Proxy(_rules)


__all__ = ['rules']
