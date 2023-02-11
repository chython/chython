# -*- coding: utf-8 -*-
#
#  Copyright 2022, 2023 Ramil Nugmanov <nougmanoff@protonmail.com>
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


def _sugar_group():
    from ... import smarts

    return smarts('[N,O;D1,D2:1][C;z1][C;z2]=[N,O:2]')


def _keto_rules():
    from ... import smarts

    rules = []
    # first atom is H-acceptor
    # second is direction

    # C-C=[O,S,NH]
    q = smarts('[N,O,S;D1;z2]=[C;D2,D3;x1;z2]')
    rules.append(q)

    # C-C=N-[C,N,O]
    q = smarts('[N;D2;z2](=[C;D2,D3;x1;z2])[C,N,O]')
    rules.append(q)

    # [C,H]-N=N-C
    q = smarts('[N;D1,D2;x1;z2]=N[C;x1]')
    rules.append(q)

    # [S,O,NR;H]-C=N
    q = smarts('[N;z2]=C[N,O,S;h1]')
    rules.append(q)

    # [NH2]-C=N-R
    q = smarts('[N;D2;z2]=C[N;D1]')
    rules.append(q)

    # S=C-[N,O;H]
    q = smarts('[S;D1;z2]=C[N,O;h1,h2]')
    rules.append(q)

    # O=C-[N,S;H]
    q = smarts('[O;D1;z2]=C[N,S;h1,h2]')
    rules.append(q)
    return rules


def _enol_rules():
    from ... import smarts

    rules = []
    # first atom is H-donor
    # second is direction

    # C=C-[OH,SH,NH2]
    q = smarts('[N,O,S;D1;z1][C;D2,D3;x1;z2]')
    rules.append(q)

    # C=C-[NH]-[C,N,O]
    q = smarts('[N;D2;z1]([C;D2,D3;x1;z2])[C,N,O]')
    rules.append(q)

    # [C,H]-[NH]-N=C
    q = smarts('[N;D1,D2;x1;z1][N;z2]=[C;x1]')
    rules.append(q)

    # [S,O,NR;H]-C=N
    q = smarts('[N,O,S;h1;z1][C;z2]=N')
    rules.append(q)

    # [NH2]-C=N-R
    q = smarts('[N;D1;z1][C;z2]=[N;D2]')
    rules.append(q)

    # S=C-[N,O;H]
    q = smarts('[N,O;h1,h2;z1][C;z2]=[S;D1]')
    rules.append(q)

    # O=C-[N,S;H]
    q = smarts('[N,S;h1,h2;z1][C;z2]=O')
    rules.append(q)
    return rules


keto_rules = Proxy(_keto_rules)
enol_rules = Proxy(_enol_rules)
sugar_group = Proxy(_sugar_group)


__all__ = ['keto_rules', 'enol_rules', 'sugar_group']
