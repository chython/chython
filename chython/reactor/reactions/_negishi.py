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
    'name': 'Negishi coupling reaction',
    'description': 'Negishi C-C coupling reaction of sp2-C halides and organozinc reagents',
    'templates': [
        {
            'A': [
                # X-Ar
                '[Cl,Br,I;D1:1]-[C;a:2]',
                # X-C-sp2
                '[Cl,Br,I;D1:1]-[C;x1,x2;z2:2]=[C;x0,x1;z2;M]',
                # Ar triflate
                '[C;a:2]-[O;D2;x1:1]-[S;x3;D4:10](=[O:11])(=[O:12])-[C;D4;z1:13](-[F;D1:14])(-[F;D1:15])-[F;D1:16]',
                # Vinyl triflates
                '[C;x0,x1;z2;M]=[C;x1;z2:2]-[O;D2;x1:1]-[S;x3;D4:10](=[O:11])(=[O:12])-[C;D4;z1:13](-[F;D1:14])(-[F;D1:15])-[F;D1:16]',
            ],
            'B': [
                # C-Zn
                '[Zn]-[C:3]',
            ],
            'product': '[A:2]-[A:3]',
            'alerts': [],
            'ufe': {
                'A': 1,
                'B': '[A:3][At;M]'
            }
        }
    ],
    'alerts': []
}
