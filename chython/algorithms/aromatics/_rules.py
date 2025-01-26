# -*- coding: utf-8 -*-
#
#  Copyright 2021-2024 Ramil Nugmanov <nougmanoff@protonmail.com>
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

    # Aromatic N-Oxide
    #
    #  : N :  >>  : [N+] :
    #    \\           \
    #     O           [O-]
    #
    q = smarts('[N;a;D3:1]=[O;D1:2]')
    atom_fix = {1: 1, 2: -1}
    bonds_fix = ((1, 2, 1),)
    # query, atom fix, bond fix, allow multimatch
    rules.append((q, atom_fix, bonds_fix, False))

    # Aromatic N-Nitride?
    #
    #  : N :  >>  : [N+] :
    #    \\           \
    #     N           [N-]
    #
    q = smarts('[N;a;D3:1]=[N;D1,D2;z2:2]')
    atom_fix = {1: 1, 2: -1}
    bonds_fix = ((1, 2, 1),)
    rules.append((q, atom_fix, bonds_fix, False))

    #
    # : [S+] : >> : S :
    #    |          \\
    #   [O-]         O
    #
    q = smarts('[S;a;D3;+:1]-[O;D1;-:2]')
    atom_fix = {1: 0, 2: 0}
    bonds_fix = ((1, 2, 2),)
    rules.append((q, atom_fix, bonds_fix, False))

    #
    # [O-]-N:C:C:[N+]=O
    #
    q = smarts('[N;a;D3;+:1](=[O;D1:2]):[C:3]:[C:4]:[N;D3:5]-[O;D1;-:6]')
    atom_fix = {5: 1, 2: -1}
    bonds_fix = ((1, 2, 1),)
    rules.append((q, atom_fix, bonds_fix, False))

    #
    # N : A : N - ?
    #  :     :
    #   C # C
    q = smarts('[N;a;D2,D3;r5:3]:1:[C;D2;r5:1]#[C;D2;r5:2]:[N;D2;r5:5]:[C,N;r5:4]:1')
    atom_fix = {}
    bonds_fix = ((1, 2, 4),)
    rules.append((q, atom_fix, bonds_fix, False))

    #
    # C:[N+]:[C-]
    #    \\
    #     O
    #
    q = smarts('[N;a;D3;+:1](=[O;D1:2])(:[C;D2,D3;-:3]):[C;D2,D3:4]')
    atom_fix = {2: -1, 3: 0}
    bonds_fix = ((1, 2, 1),)
    rules.append((q, atom_fix, bonds_fix, False))

    #
    #  O=[N+] : C
    #     :     :
    #    O : N : C
    q = smarts('[N;a;D3;r5;+:1]:1(=[O;D1:2]):[O;D2;r5:3]:[N;D2,D3;r5:4]:[C;D2,D3;r5:5]:[C;D2,D3;r5:6]:1')
    atom_fix = {}
    bonds_fix = ((1, 3, 1), (1, 6, 1), (3, 4, 1), (4, 5, 1), (5, 6, 2))
    rules.append((q, atom_fix, bonds_fix, False))

    # bad complex representation
    #  :           :
    #   N - [M]  >> N ... [M]
    #  :           :
    q = smarts('[N;a;D3:1]-[M:2]')
    atom_fix = {}
    bonds_fix = ((1, 2, 8),)
    rules.append((q, atom_fix, bonds_fix, True))
    return rules


def _freaks():
    from ... import smarts

    rules = []

    q = smarts('[N,O,S;D2;r5;z1]1[A;r5]=,:[A;r5][A;r5]:[A;r5]1')
    rules.append(q)

    q = smarts('[N;D3;r5;z1]1[A;r5]=,:[A;r5][A;r5]:[A;r5]1')
    rules.append(q)
    return rules


rules = Proxy(_rules)
freak_rules = Proxy(_freaks)


__all__ = ['rules', 'freak_rules']
