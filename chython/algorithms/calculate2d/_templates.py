# -*- coding: utf-8 -*-
#
#  Copyright 2024, 2025 Ramil Nugmanov <nougmanoff@protonmail.com>
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


def aligner(xy):
    x0, y0 = min(xy, key=lambda x: x[0])
    return [(round(x - x0, 4), round(y - y0, 4)) for x, y in xy]


def _rules():
    from ... import smarts

    rules = []

    # bicyclo[1.1.1]pentane
    q = smarts('[A;r4;D2;z1:2]-1-[A;D3,D4:1]-2-[A;r4;D2:5]-[A;D3,D4:3]-1-[A;r4;D2:4]-2')
    xy = {1: (0.0, 0.0), 2: (0.6674, -0.485), 3: (1.3348, 0.0), 4: (1.0799, 0.7846), 5: (0.2549, 0.7846)}
    sub = {1: (-0.825, 0.0), 3: (2.1598, 0.0)}
    rules.append((q, xy, sub))

    # adamantane
    q = smarts('[A;D3;r6:1]-1-2-[A;D2:2][A:6]-3-[A:7][A:8]([A;D2:3]-1)[A:9][A:10]([A;D2:4]-2)[A:5]-3')
    xy = {1: (0.0084, 0.8368), 2: (0.2398, 1.2541), 3: (0.2542, 0.4277), 4: (-0.4687, 0.8284), 5: (-0.7145, 1.2375),
          6: (0.0, 1.65), 7: (0.7144, 1.2375), 8: (0.7144, 0.4125), 9: (0.0, 0.0), 10: (-0.7145, 0.4125)}
    sub = {6: (0.0, 2.475), 7: (1.4289, 1.65), 8: (1.4289, 0.0),
           9: (0.0, -0.825), 10: (-1.429, 0.0), 5: (-1.429, 1.65)}
    rules.append((q, xy, sub))

    # cyclopropane
    q = smarts('[A;r3:1]-1-;@[A;r3:2]-;@[A;r3:3]-1')
    xy = {1: (0.825, 0.0), 2: (0.0, 0.0), 3: (0.4125, -0.7145)}
    rules.append((q, xy, {}))
    return rules


rules = Proxy(_rules)


__all__ = ['rules']
