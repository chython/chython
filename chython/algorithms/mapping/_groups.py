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


def _xonyl_groups():
    from ... import smarts

    rules = []
    # atom 1 - axis
    # atoms 2, 3 - swapping
    q = smarts('[N;D3;z2;+]([O;D1;-])=[O;D1]')
    rules.append(q)

    q = smarts('[N,O,S,Se;D1;z2:3]=[C,P,As,Sb,Bi,S,Se,Te,Po,Cl,Br,I,At:1][N,O,S,Se;D1:2]')
    rules.append(q)

    q = smarts('[N,O,S,Se;D1;z2:3]=[C,P,As,Sb,Bi,S,Se,Te,Po,Cl,Br,I,At:1][N,O,S,Se;D1;-:2]')
    rules.append(q)
    return rules


def _substituents():
    """
    Rules for switchable functional groups remapping
    """
    from ... import smarts

    rules = []

    # nitro addition
    #
    #      O [3]     [2] O -- *
    #     //            /:
    # * - N+   >>  * = N+
    #      \            \
    #       O- [2]       O- [3]
    #
    q = smarts('[N;D3;z2,z4;+]([O;D1;-])-,:[O;D2]')
    rules.append((q, ((2, 3), (3, 2))))

    #
    #      NH [3]          NH2 [3]
    #     //               /
    # * - X       >>  * - X
    #      \              :
    #       NH2 [2]       N [2]
    #
    q = smarts('[N;D2;h0;z2,z4:3]=,:[C;D3:1][N:2]')
    rules.append((q, ((2, 3), (3, 2))))

    #
    #      O [3]        [3] O - R
    #     //               /
    # * - X       >>  * - X
    #      \              \\
    #       OH [2]         O [2]
    #
    q = smarts('[N,O,S,Se;D1;z2:3]=[C,P,As,Sb,Bi,S,Se,Te,Po,Cl,Br,I,At:1][N,O,S,Se;D2:2]')
    rules.append((q, ((2, 3), (3, 2))))

    #
    #      O [3]            A
    #     //               /
    # * - X       >>  * - X
    #      \              \\
    #       OH [2]         O [2]
    #
    q = smarts('[N,O,S,Se;D1;z2:3]=[C,P,As,Sb,Bi,S,Se,Te,Po,Cl,Br,I,At:1][A:2]')
    rules.append((q, ((3, 2),)))  # possible only: (3, 2)
    return rules


xonyl_groups = Proxy(_xonyl_groups)
substituents_groups = Proxy(_substituents)


__all__ = ['xonyl_groups', 'substituents_groups']
