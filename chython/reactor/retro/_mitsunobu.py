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
    'name': 'Mitsunobu reaction',
    'description': 'Phenol-Alcohol Phenol-Acid Acid-Alcohol couplings',
    'templates': [
        {
            # Ph-O-Alk
            'product': '[O;D2;x0;z1:1](-;!@[C;a;M])[C;z1;x1:2]',
            'reactants': [
                '[A:1]',
                '[A:2]-[O;M]'
            ]
        },
        {
            # Ac-O-Alk
            'product': '[O;D2;x0;z1:1](-;!@[C;x2;z2;M]=[O;M])[C;z1;x1:2]',
            'reactants': [
                '[A:1]',
                '[A:2]-[O;M]'
            ]
        },
        {
            # Ph-O-Ac
            'product': '[O;D2;x0;z1:1](-;!@[C;D3;x2;z2:2]=[O;M])[C;a;M]',
            'reactants': [
                '[A:1]',
                '[A:2]-[O;M]'
            ]
        }
    ]
}
