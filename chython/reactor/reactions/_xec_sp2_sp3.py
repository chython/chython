# -*- coding: utf-8 -*-
#
#  Copyright 2025 Kostia Chernichenko
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
                '[Cl,Br,I;D1:1]-[C;z2;r5,r6:2]',
                # Ar triflate
                '[C;a:2]-[O;D2;x1:1]-[S;x3;D4:10](=[O:11])(=[O:12])-[C;D4:13](-[F;D1:14])(-[F;D1:15])-[F;D1:16]'
            ],
            'B': [
                # sp3-C-X
                '[Cl,Br,I;D1:3]-[C;z1:4]'
            ],
            'product': '[A:2]-[A:4]',
            'alerts': [],
            'ufe': {
                'A': '[A:2][At;M]',
                'B': 3
            }
        }
    ],
    'alerts': []
}
