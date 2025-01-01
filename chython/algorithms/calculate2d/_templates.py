# -*- coding: utf-8 -*-
#
#  Copyright 2024 Ramil Nugmanov <nougmanoff@protonmail.com>
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
    q = smarts('[A;r4;D2:2]-1-[A;D3,D4:1]-2-[A;r4;D2:5]-[A;D3,D4:3]-1-[A;r4;D2:4]-2')
    xy = [(0.0, 0.0), (0.6674, 0.485), (1.3348, 0.0), (1.0799, -0.7846), (0.2549, -0.7846)]
    rules.append((q, xy))

    # adamantane
    q = smarts('[A;D3;r6;z1:1]-1-2-[A:2][A:6]-3-[A:7][A:8]([A:3]-1)[A:9][A:10]([A:4]-2)[A:5]-3')
    xy = [(0.9254, -0.0085), (1.3673, -0.2636), (0.4835, -0.2636), (0.9254, 0.5018), (1.4053, 0.8113),
          (1.8737, 0.0), (1.4053, -0.8113), (0.4684, -0.8113), (0.0, 0.0), (0.4684, 0.8113)]
    rules.append((q, xy))
    return rules


rules = Proxy(_rules)


__all__ = ['rules']
