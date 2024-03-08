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
    'name': 'Fischer esterification',
    'description': 'Esters formation from alcohols and acids',
    'templates': [
        # reactants sets fully mixable
        {
            'A': [
                # C(=O)O
                '[O;D1;x0;z1:2]-[C;x2;z2:1]=[O;M]',
            ],
            'B': [
                # CO
                '[O;D1;x0;z1:3]-[C;x1;z1;M]'
            ],
            'product': '[A:1]-[A:3]',
            # condition-specific untolerant groups
            'alerts': [
                '[S;D1;x0;z1][C;x1;z1]',  # thiol
                '[O,S;D1;z1][A;a]'  # [thia]phenol
            ],
            'ufe': {
                'A': 2,
                'B': '[A:3][At;M]'
            }
        }
    ],
    'alerts': []  # global untolerant groups
}
