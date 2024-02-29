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
    'name': 'Amidation reaction',
    'description': 'Amide Coupling with Amines and Acids/Halo-Anhydrides',
    'templates': [
        {
            'product': '[N;D2,D3;z1:1]-;!@[C;x2;z2:2]=[O;M]',  # any SP3 nitrogen with carboxy
            'reactants': [
                '[A:1]',
                '[A:2]-[O;M]'
            ]
        },
        {
            'product': '[N;D2;z2;x0:1]-;!@[C;x2;z2:2]=[O;M]',  # C=N-C(=O)R
            'reactants': [
                '[A:1]',
                '[A:2]-[O;M]'
            ]
        }
    ]
}
