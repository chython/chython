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


def _rules():
    from ... import smarts

    rules = []


    #   C
    #  /  \
    # C    C
    # |    |
    # C    C
    #  \  /
    #   C
    q = smarts('[A;r6]!#;@1!#;@[A]!#;@[A]!#;@[A]!#;@[A]!#;@[A]1')
    xy = [(0, 2.31), (0, 0.77), (1.3337, 0), (2.6674, 0.77), (2.6674, 2.31), (1.3337, 3.08)]
    rules.append((q, xy))

    #
    #  C-,=C
    # ||   ||
    # C    C
    #  \  /
    #   C1
    q = smarts('[A;r5]-,:;@1-,:;@[A]!#;@[A]!#;@[A]!#;@[A]1')
    xy = [(1.2459, 0), (0, 0.9052), (0.4759, 2.3698), (2.0159, 2.3698), (2.4918, 0.9052)]
    rules.append((q, xy))

    #
    # C=,-C
    # \   /
    #   C
    #
    q = smarts('[A;r3]1-,=;@[A]-;@[A]1')
    xy = [(0, 0), (1.54, 0), (.77, -1.333)]
    rules.append((q, xy))

    #
    # C=,-C
    # |   |
    # C=,-C
    #
    q = smarts('[A;r4]-;@1-,=;@[A]-;@[A]-,=;@[A]1')
    xy = [(0, 0), (1.54, 0), (1.54, -1.54), (0, -1.54)]
    rules.append((q, xy))
    return rules


rules = Proxy(_rules)


__all__ = ['rules']
