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


def _fixed_rules():
    """
    Rules working without Morgan.
    These rules are used for charge canonization for heterocycles
    """
    from ... import smarts

    rules = []

    # N1C=CN2[NH+]=CC=C12>>N1C=CC2=[NH+]C=CN12
    # first and second atoms are ryrrole-like
    # second atom will be charged
    # if fix is True, third atom will be uncharged else first
    q = smarts('[N;a;r5;+:1]:1:[N;r5:3]:2:[C;r5](:[N;r5:2]:[C;r5]:[C;r5]:2):[C;r5]:[C;r5]:1')
    rules.append((q, False))

    # N1C=C[N+]2=C1C=CN2>>N1C=CC2=[NH+]C=CN12
    q = smarts('[N;a;r5;+:3]:1:2:[N;r5:1]:[C;r5]:[C;r5]:[C;r5]:1:[N;r5:2]:[C;r5]:[C;r5]:2')
    rules.append((q, True))

    # N1C=C2C=CN[N+]2=C1>>N1C=CC2=C[NH+]=CN12
    q = smarts('[N;a;r5;+:3]:1:2:[N;r5:1]:[C;r5]:[C;r5]:[C;r5]:1:[C;r5]:[N;r5:2]:[C;r5]:2')
    rules.append((q, True))

    # N1C=C2NC=C[N+]2=C1>>N1C=CN2C=[NH+]C=C12
    q = smarts('[N;a;r5;+:3]:1:2:[C;r5](:[N;r5:1]:[C;r5]:[C;r5]:1):[C;r5]:[N;r5:2]:[C;r5]:2')
    rules.append((q, True))

    # N1C2=[NH+]C=CC2=CC=C1>>N1C=CC2=CC=C[NH+]=C12
    q = smarts('[N;a;r5;+:1]:1:[C;r6]:2:[N;r6:2]:[C;r6]:[C;r6]:[C;r6]:[C;r5]:2:[C;r5]:[C;r5]:1')
    rules.append((q, False))

    # N1C=C2C(=CC=[NH+]2)C=C1>>N1C=CC2=CC=[NH+]C=C12
    q = smarts('[N;a;r5;+:1]:1:[C;r6]:2:[C;r6]:[N;r6:2]:[C;r6]:[C;r6]:[C;r5]:2:[C;r5]:[C;r5]:1')
    rules.append((q, False))

    # C1=CC=2C(=[NH+]1)C=CNC=2>>N1C=CC2=C[NH+]=CC=C12
    q = smarts('[N;a;r5;+:1]:1:[C;r6]:2:[C;r6]:[C;r6]:[N;r6:2]:[C;r6]:[C;r5]:2:[C;r5]:[C;r5]:1')
    rules.append((q, False))

    # C=1C=[NH+]C=2C=1NC=CC=2>>N1C=CC2=[NH+]C=CC=C12
    q = smarts('[N;a;r5;+:1]:1:[C;r6]:2:[C;r6]:[C;r6]:[C;r6]:[N;r6:2]:[C;r5]:2:[C;r5]:[C;r5]:1')
    rules.append((q, False))

    # C1=2C=[NH+]C=C1NC=CC=2>>N1C=C2C=CC=[NH+]C2=C1
    q = smarts('[N;a;r5;+:1]:1:[C;r5]:[C;r6]:2:[N;r6:2]:[C;r6]:[C;r6]:[C;r6]:[C;r5]:2:[C;r5]:1')
    rules.append((q, False))

    # C1=2C=CNC=C1C=[NH+]C=2>>N1C=C2C=C[NH+]=CC2=C1
    q = smarts('[N;a;r5;+:1]:1:[C;r5]:[C;r6]:2:[C;r6]:[N;r6:2]:[C;r6]:[C;r6]:[C;r5]:2:[C;r5]:1')
    rules.append((q, False))
    return rules


def _morgan_rules():
    """
    Rules working with Morgan.
    These rules are for charge canonization
    """
    from ... import smarts

    rules = []

    # N1C=CC2=CC=[NH+]N12
    q = smarts('[N;a;r5;+:1]:1:[N;D3]:2:[N;r5:2]:[C;r5]:[C;r5]:[C;D3]:2:[C;r5]:[C;r5]:1')
    rules.append((q, False))

    # N1C=CC2=[N+]1NC=C2
    q = smarts('[N;a;D3;r5;+:3]:1:2:[N;r5:1]:[C;r5]:[C;r5]:[C;D3]:1:[C;r5]:[C;r5]:[N;r5:2]:2')
    rules.append((q, True))

    # N1C=CN2C=C[NH+]=C12
    q = smarts('[N;a;r5;+:1]:1:[C;D3]:2:[N;r5:2]:[C;r5]:[C;r5]:[N;D3]:2:[C;r5]:[C;r5]:1')
    rules.append((q, False))

    # N1C=C[N+]2=C1NC=C2
    q = smarts('[N;a;D3;r5;+:3]:1:2:[C;D3](:[N;r5:1]:[C;r5]:[C;r5]:1):[N;r5:2]:[C;r5]:[C;r5]:2')
    rules.append((q, True))

    # imidazole
    q = smarts('[N;a;r5;+:1]:1:[C;r5]:[N;r5:2]:[C;r5]:[C;r5]:1')
    rules.append((q, False))

    # pyrazole+
    q = smarts('[N;a;r5;+:1]:1:[N;r5:2]:[C;r5]:[C;r5]:[C;r5]:1')
    rules.append((q, False))
    return rules


fixed_rules = Proxy(_fixed_rules)
morgan_rules = Proxy(_morgan_rules)


__all__ = ['fixed_rules', 'morgan_rules']
