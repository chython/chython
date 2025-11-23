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
    'name': 'Mitsunobu reaction',
    'description': 'Mitsunobu coupling of aliphatic alcohols and various acidic nucleophiles such as phenols, carboxylic acids, thiols, etc.',
    'templates': [
        {
            'A': [
                # Chiral alcohol
                '[O;D1:1]-[C;@@;z1:2](-[C;M])-[C;M]',
            ],
            'B': [
                # Ar-OH
                '[O;D1;x0;z1:3][C;a;M]',
                # COOH
                '[O;D1;x0;z1:3]-[C;x2;z2;M](=[O;M])-[C;M]',
                # Thiol
                '[S;D1;x0;z1:3][C;M]',
                # Imides
                '[N;D2;x0;z1:3](-[C;x2,x3;z2;M](=[O;M]))-[C;x2,x3;z2;M](=[O;M])',
                # Sulfonamides
                '[N;D1,D2;z1:3]-[S;x3;M](=[O;M])(=[O;M])-[C;M]',
                # Azide anion
                '[N;-1:3]=[N;+1;M]=[N;-1;M]',
                # N-H heterocycles
                '[N;D2;a;h1:3](:[C,N,S;a;M]):[C,N;a;M]',
                '[H:4][N;D3:3](:[C,N,S;a;M]):[C,N;a;M]',
            ],
            'product': '[A;@@:2]-[A:3]',
            'alerts': [],
            'ufe': {
                'A': 1,
                'B': '[A:3][At;M]'
            }
        },
        {
            'A': [
                # OH-Alk
                '[O;D1:1]-[C;z1:2]',
            ],
            'B': [
                # Ar-OH
                '[O;D1;x0;z1:3][C;a;M]',
                # COOH
                '[O;D1;x0;z1:3]-[C;x2;z2;M](=[O;M])-[C;M]',
                # Thiol
                '[S;D1;x0;z1:3][C;M]',
                # Imides
                '[N;D2;x0;z1:3](-[C;x2,x3;z2;M](=[O;M]))-[C;x2,x3;z2;M](=[O;M])',
                # Sulfonamides
                '[N;D1,D2;z1:3]-[S;x3;M](=[O;M])(=[O;M])-[C;M]',
                # Azide anion
                '[N;-1:3]=[N;+1;M]=[N;-1;M]',
                # N-H heterocycles
                '[N;D2;a;h1:3](:[C,N,S;a;M]):[C,N;a;M]',
                '[H:4][N;D3:3](:[C,N,S;a;M]):[C,N;a;M]',
                # Thiophenol

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
