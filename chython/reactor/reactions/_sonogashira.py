# -*- coding: utf-8 -*-
#
#  Copyright 2022-2024 Ramil Nugmanov <nougmanoff@protonmail.com>
#  Copyright 2023 Timur Gimadiev <timur.gimadiev@gmail.com>
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
    'description': 'Sonogashira reaction, C-C coupling reaction. It employs a palladium catalyst as well as copper'
                   'co-catalyst',
    'templates': [
        {
            'A': [
                # HC#C-R
                '[C;D1;x0;z3:1]#[C;D2;x0;M]'
            ],
            'B': [
                # Ar-Hal
                '[Cl,Br,I;D1:3]-[C;a:2]',
                # C=C-Hal
                '[Cl,Br,I;D1:3]-[C;x1;z2:2]=[C;x0;z2;M]',
                # R-C(=O)-Hal
                '[Cl,Br,I;D1:3]-[C;x2;z2:2]=[O;M]'
            ],
            'product': '[A:1]-[A:2]',
            'alerts': [],
            'ufe': {
                'A': '[A:1][At;M]',
                'B': 3
            }
        }
    ],
    'alerts': []
}
