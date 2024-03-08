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
    'name': 'Amidation Reaction',
    'description': 'Amides formation from acids and amines',
    'templates': [
        {
            'A': [
                # [H,R]COOH
                '[O;x0;z2;M]=[C;x2:1][O;D1:2]'
            ],
            'B': [
                # Ar-NH2
                '[N;D1;x0;z1:3][C;a;M]',
                # Alk-NH2
                '[N;D1;x0;z1:3][C;z1;x1;M]',
                # Ar-NH-Ar
                '[N;D2;x0;z1:3]([C;a;M])[C;a;M]',
                # Alk-NH-Ar
                '[N;D2;x0;z1:3]([C;a;M])[C;z1;x1;M]',
                # Alk2NH
                '[N;D2;x0;z1:3]([C;z1;x1;M])[C;z1;x1;M]',
                # N1COCCC1
                '[N;D2;x0;z1;r5,r6,r7,r8:3]([C;z1;x2;M]-;@[O;M])[C;z1;x1;M]',
                # CNO[R,H]
                '[N;D2;x1;z1:3]([O;x1;z1;M])[C;z1;x1;M]',
                # C[NH]NAc
                '[N;D2;x1;z1:3]([N;D2;z1;x1;M][C;x2;z2;M]=[O;M])[C;z1;x1;M]'
            ],
            'product': '[A:1]-[A:3]',
            'alerts': [],
            'ufe': {
                'A': 2,  # use existing terminal atom
                'B': '[A:3][At;M]'  # add temporary terminal atom
            }
        }
    ],
    'alerts': ['[O;D1;x0;z1][C;z1;x1]', '[O;D1;z1][C,N;a]']  # global untolerant groups
}
