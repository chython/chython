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
    'name': 'Macmillan',
    'description': 'Deoxygenative C-C coupling reaction',
    'templates': [
        {
            'A': [
                # Hal-Ar
                '[Cl,Br,I;D1:1]-[C;a:2]'
            ],
            'B': [
                # AlkCOOH
                '[C;z1:3]-[C;x2;z2:4](=[O:5])-[O;D1;x0;z1:6]',
                # Redox ester
                '[C;z1:3]-[C;x2;z2:4](=[O:5])-[O;D2;x1;z1:6]-[N;z1;x1;D3:7]1-[C;x2;z2:8](=[O:9])-[C;x0;z2,z4:10]!#[C;x0;z2,z4:11]-[C;x2;z2:12](=[O:13])-1'
            ],
            'product': '[A:2]-[A:3]',
            'alerts': [],
            'ufe': {
                'A': 1,
                'B': 3
            }
        }
    ],
    'alerts': []
}
