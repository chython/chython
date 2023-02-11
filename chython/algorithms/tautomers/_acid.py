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


def _stripped_rules():
    from ... import smarts

    rules = []

    # Ammonia, H-imine,guanidine,amidine. [H][N+]=,-,:
    q = smarts('[N;h1,h2,h3,h4;+]')
    rules.append(q)
    return rules


def _rules():
    from ... import smarts

    rules = _stripped_rules()

    # Phenoles
    q = smarts('[O,S,Se;D1;z1][C,N;a]')
    rules.append(q)

    # Oxo-acids
    q = smarts('[O,S,Se;D1;z1][C,N,P,S,Se,Cl,Br,I]=O')
    rules.append(q)

    # Nitro acid
    q = smarts('[N;D3;z2;+]([O;D1:1])([O-])=O')
    rules.append(q)

    q = smarts('[F,Cl,Br,I;D0]')
    rules.append(q)
    return rules


stripped_rules = Proxy(_stripped_rules)
rules = Proxy(_rules)


__all__ = ['stripped_rules', 'rules']
