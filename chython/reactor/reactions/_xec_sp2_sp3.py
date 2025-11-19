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
    'name': 'XEC',
    'description': 'Cross-electrophile C-sp2-X C-sp3-X coupling reaction',
    'templates': [
        {
            'A': [
                # Hal-Ar
                '[Cl,Br,I;D1:1]-[C;a:2]',
                # Hal-pseudoaromatic, more specifically C5, C6 vinylic
                '[Cl,Br,I;D1:1]-[C;z2;r5,r6:2]'
            ],
            'B': [
                # sp3-C-X
                '[Cl,Br,I;D1:3]-[C;z1:4]'
            ],
            'product': '[A:2]-[A:4]',
            'alerts': [],
            'ufe': {
                'A': 1,
                'B': 3
            }
        }
    ],
    'alerts': []
}
