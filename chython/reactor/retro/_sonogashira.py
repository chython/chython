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

template = {
    'name': 'Sonogashira reaction',
    'description': 'Alkyne Ar-X CSP2-X Ac-X couplings',
    'templates': [
        {
            # Ar
            'product': '[C;D2;x0;z3:1](-;!@[C;a:2])#[C;D2;x0;M]',
            'reactants': [
                '[A:1]',
                '[A:2]-[Br;M]'
            ]
        },
        {
            # Ac
            'product': '[C;D2;x0;z3:1](-;!@[C;D3;x1;z2:2]=[O;M])#[C;D2;x0;M]',
            'reactants': [
                '[A:1]',
                '[A:2]-[Cl;M]'
            ]
        },
        {
            # CSP2
            'product': '[C;D2;x0;z3:1](-;!@[C;x0;z2:2])#[C;D2;x0;M]',
            'reactants': [
                '[A:1]',
                '[A:2]-[Br;M]'
            ]
        }
    ]
}
