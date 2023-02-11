# -*- coding: utf-8 -*-
#
#  Copyright 2021-2023 Ramil Nugmanov <nougmanoff@protonmail.com>
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


def _rules():
    from ... import smarts

    rules = []

    #
    #        R - N - C                 R - N - C
    #           /    ||                   /    ||
    #  A - M = C     || >>    A - M .. [C-]    ||
    #           \    ||                  \\    ||
    #        R - N - C               R - [N+]- C
    #
    q = smarts('[M:1]=[C:2]-1-[N;D3;x0;z1:3]-[C;z2:5]-,=[C;z2:6]-[N;D3;x0;z1:4]-1')
    atom_fix = {2: (-1, None), 3: (1, None)}  # atom: (charge diff, new radical state or None)
    bonds_fix = ((1, 2, 8), (2, 3, 2))
    rules.append((q, atom_fix, bonds_fix))

    #
    # cyanide
    #
    q = smarts('[M:3]-[C:1]#[N:2]')
    # NOTE: first fixed atom is metal. required for preventing errors in charge.
    # based on ordered dict nature.
    atom_fix = {3: (1, None), 1: (-1, None)}
    bonds_fix = ((1, 3, 8),)
    rules.append((q, atom_fix, bonds_fix))

    q = smarts('[M:3]-[C-:1]#[N:2]')
    atom_fix = {}
    bonds_fix = ((1, 3, 8),)
    rules.append((q, atom_fix, bonds_fix))

    #
    # cyanate/fulminate
    #
    q = smarts('[M:1]-[O,S;D2:2]-[C:3]#[N:4]')
    atom_fix = {1: (1, None), 4: (-1, None)}
    bonds_fix = ((1, 2, 8), (2, 3, 2), (3, 4, 2))
    rules.append((q, atom_fix, bonds_fix))

    q = smarts('[M:1]-[N:2]=[C:3]=[O,S;D1:4]')
    atom_fix = {1: (1, None), 2: (-1, None)}
    bonds_fix = ((1, 2, 8),)
    rules.append((q, atom_fix, bonds_fix))

    q = smarts('[M:1]-[O,S;D2:2]-[N+:3]#[C-:4]')
    atom_fix = {1: (1, None), 2: (-1, None)}
    bonds_fix = ((1, 2, 8),)
    rules.append((q, atom_fix, bonds_fix))

    q = smarts('[M:1]-[C:2]#[N+:3]-[O,S;D1-:4]')
    atom_fix = {1: (1, None), 2: (-1, None)}
    bonds_fix = ((1, 2, 8),)
    rules.append((q, atom_fix, bonds_fix))

    #
    # carbonyl
    #
    q = smarts('[M:3]-[C:1]#[O:2]')
    atom_fix = {1: (-1, None), 2: (1, None)}
    bonds_fix = ((1, 3, 8),)
    rules.append((q, atom_fix, bonds_fix))

    q = smarts('[M:3]-[C:1]=[O:2] |^1:1|')
    atom_fix = {1: (-1, False), 2: (1, None)}
    bonds_fix = ((1, 3, 8), (1, 2, 3))
    rules.append((q, atom_fix, bonds_fix))

    q = smarts('[M:3]-[C;D2:1]=[O:2]')
    atom_fix = {1: (-1, None), 2: (1, None)}
    bonds_fix = ((1, 3, 8), (1, 2, 3))
    rules.append((q, atom_fix, bonds_fix))

    q = smarts('[M:3]-[C:1](-[M:4])=[O:2]')
    atom_fix = {1: (-1, None), 2: (1, None)}
    bonds_fix = ((1, 3, 8), (1, 4, 8), (1, 2, 3))
    rules.append((q, atom_fix, bonds_fix))

    #
    # Ferrocene covalent uncharged
    #
    q = smarts('[M:1]-1-2-3-4-[C:2]-5-[C:3]-1-[C:4]-2-[C:5]-3-[C:6]-4-5')
    atom_fix = {1: (1, None), 2: (-1, None)}
    bonds_fix = ((1, 2, 8), (1, 3, 8), (1, 4, 8), (1, 5, 8), (1, 6, 8), (3, 4, 2), (5, 6, 2))
    rules.append((q, atom_fix, bonds_fix))

    #
    # Ferrocene covalent uncharged. invalid valence.
    #
    q = smarts('[M:1]-1-2-3-4-[C:2]-5-[C:3]-1=[C:4]-2-[C:5]-3=[C:6]-4-5')
    atom_fix = {1: (1, None), 2: (-1, None)}
    bonds_fix = ((1, 2, 8), (1, 3, 8), (1, 4, 8), (1, 5, 8), (1, 6, 8), (3, 4, 2), (5, 6, 2))
    rules.append((q, atom_fix, bonds_fix))

    #
    # Ferrocene covalent radical carbon
    #
    q = smarts('[M:1]-1-2-3-4-[C:2]-5-[C:3]-1=[C:4]-2-[C:5]-3=[C:6]-4-5 |^1:1|')
    atom_fix = {1: (1, None), 2: (-1, False), 3: (0, False), 4: (0, False), 5: (0, False), 6: (0, False)}
    bonds_fix = ((1, 2, 8), (1, 3, 8), (1, 4, 8), (1, 5, 8), (1, 6, 8), (3, 4, 2), (5, 6, 2))
    rules.append((q, atom_fix, bonds_fix))

    #
    # Ferrocene covalent charged
    #
    q = smarts('[M:1]-1-2-3-4-[C-:2]-5-[C:3]-1=[C:4]-2-[C:5]-3=[C:6]-4-5')
    atom_fix = {}
    bonds_fix = ((1, 2, 8), (1, 3, 8), (1, 4, 8), (1, 5, 8), (1, 6, 8))
    rules.append((q, atom_fix, bonds_fix))

    #
    # Ferrocene coordinate radical carbon
    #
    q = smarts('[M:1]~1~2~3~4~[C:2]-5-[C:3]~1=[C:4]~2-[C:5]~3=[C:6]~4-5 |^1:1|')
    atom_fix = {1: (1, None), 2: (-1, False), 3: (0, False), 4: (0, False), 5: (0, False), 6: (0, False)}
    bonds_fix = ((3, 4, 2), (5, 6, 2))
    rules.append((q, atom_fix, bonds_fix))

    #
    # Allyl complexes
    #
    q = smarts('[M:1]~1~2~[C;z2:2]=[C:3]~1-[C:4]~2 |^1:3|')
    atom_fix = {1: (1, None), 4: (-1, False)}
    bonds_fix = ()
    rules.append((q, atom_fix, bonds_fix))

    #
    # Phosphines
    #
    q = smarts('[M:1]-[P;D4;z1;+:2]-C')
    atom_fix = {1: (1, None), 2: (-1, None)}
    bonds_fix = ((1, 2, 8),)
    rules.append((q, atom_fix, bonds_fix))

    q = smarts('[M:1]-[P;D4;z1:2]-C')
    atom_fix = {}
    bonds_fix = ((1, 2, 8),)
    rules.append((q, atom_fix, bonds_fix))

    # Amines
    q = smarts('[M:1]-[N;z1;+:2]-C')
    atom_fix = {1: (1, None), 2: (-1, None)}
    bonds_fix = ((1, 2, 8),)
    rules.append((q, atom_fix, bonds_fix))

    q = smarts('[M:1]-[N;z1:2]-C')
    atom_fix = {}
    bonds_fix = ((1, 2, 8),)
    rules.append((q, atom_fix, bonds_fix))

    # ethers
    q = smarts('[M:1]-[O;D3;+:2](-C)-C')
    atom_fix = {1: (1, None), 2: (-1, None)}
    bonds_fix = ((1, 2, 8),)
    rules.append((q, atom_fix, bonds_fix))

    compiled_rules = []
    for q, atom_fix, bonds_fix in rules:
        any_atoms = [n for n, a in q.atoms() if a.atomic_symbol == 'A' and n not in atom_fix]
        any_atoms.extend(n for n, a in q.atoms() if a.atomic_symbol == 'M')
        compiled_rules.append((q, atom_fix, bonds_fix, any_atoms, False))
    return compiled_rules


rules = Proxy(_rules)


__all__ = ['rules']
