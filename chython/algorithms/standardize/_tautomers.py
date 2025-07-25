# -*- coding: utf-8 -*-
#
#  Copyright 2025 Ramil Nugmanov <nougmanoff@protonmail.com>
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
    """
    Rules working without Morgan.
    These rules are used for charge canonization for heterocycles
    """
    from ... import smarts

    rules = []
    # the order of patterns is important!
    # rules are divided into types:
    # * rules with known canonical hydrogen position - encoded as (smarts, 0),
    #   where atom 1 is desired H position, atom 2 is current H position
    # * rules with variable position - encoded as (smarts, 1), where H position must be selected between atom 1 or 2
    #   H must be on atom 1 in the pattern
    # * rules with variable position with out-of-selection hydrogen position - encoded as (smarts, 2),
    #   where H position must be moved from 3 to 1 or 2
    # * rules with variable position or hydrogen and double bond - encoded as (smarts, 3),
    #   where H position must be calculated globally

    # 1,2,4‐triazole
    # N1C=NC=N1>>N1C=NN=C1
    q = smarts('[N;a;r5;D2;h1;x1:2]:1:[N;r5;D2;h0;x1]:[C;r5]:[N;r5;D2;h0;x0:1]:[C;r5]:1')
    rules.append((q, 0))

    # tetrazole
    # N1N=CN=N1>>N1C=NN=N1
    q = smarts('[N;a;r5;D2;h0;x1:1]:1:[N;r5;D2;h1;x2:2]:[N;r5;D2;h0;x2]:[N;r5;D2;h0;x1]:[C;r5]:1')
    rules.append((q, 0))

    # pyrazole
    q = smarts('[N;a;r5;D2;h1;x1:1]:1:[N;r5;D2;h0;x1:2]:[C;r5]:[C;r5]:[C;r5]:1')
    rules.append((q, 1))

    # imidazole
    q = smarts('[N;a;r5;D2;h1;x0:1]:1:[C;r5]:[N;r5;D2;h0;x0:2]:[C;r5]:[C;r5]:1')
    rules.append((q, 1))

    # 1,2,3-triazole
    q = smarts('[N;a;r5;D2;h1;x1:1]:1:[N;r5;D2;h0;x2]:[N;r5;D2;h0;x1:2]:[C;r5]:[C;r5]:1')
    rules.append((q, 1))

    q = smarts('[N;a;r5;D2;h0;x1:1]:1:[N;r5;D2;h1;x2:3]:[N;r5;D2;h0;x1:2]:[C;r5]:[C;r5]:1')
    rules.append((q, 2))

    # amidine, guanidine, etc
    q = smarts('[N;h0,h1;z2:1]=[C;z2;x2,x3:3]-[N;h1,h2;z1:2]')
    rules.append((q, 3))
    return rules


rules = Proxy(_rules)


__all__ = ['rules']
