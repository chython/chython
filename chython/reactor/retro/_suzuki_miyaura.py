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
    'name': 'Suzuki-Miyaura reaction',
    'description': 'Ar-Ar Ar-CSP2 CSP2-CSP2 couplings',
    'templates': [
        {
            'product': '[C;a;D3:1]-;!@[C;a:2]',
            'reactants': [
                '[A:1]-[B;M]([O;M])[O;M]',
                '[A:2]-[Br;M]'
            ]
        },
        {
            'product': '[C;a;D3:1]-;!@[C;z2:2]=[C;M]',
            'reactants': [
                '[A:1]-[B;M]([O;M])[O;M]',
                '[A:2]-[Br;M]'
            ]
        },
        {
            # reverse
            'product': '[C;a;D3:1]-;!@[C;z2:2]=[C;M]',
            'reactants': [
                '[A:1]-[Br;M]',
                '[A:2]-[B;M]([O;M])[O;M]'
            ]
        },
        {
            'product': '[C;D2,D3;z2:1](-;!@[C;z2:2]=[C;M])=[C;M]',
            'reactants': [
                '[A:1]-[B;M]([O;M])[O;M]',
                '[A:2]-[Br;M]'
            ]
        }
    ]
}
