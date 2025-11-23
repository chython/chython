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
    'name': 'Aliphatic C-N coupling reaction',
    'description': 'Nucleophilic substitution of aliphatic chlorides, bromides, iodides, sulfonates and sulfates with nitrogen nucleophiles such as amines, amides, hydrazines, hydrazones, N-H-aromatic heterocycles, etc.',
    'templates': [
        {
            'A': [
                # Hal-Alk
                '[Cl,Br,I;D1:1]-[C;z1:2]',
                # Alk-Sulfonate, i.e. mesylates, tosylates, triflates, nosylates, brosylates, etc.
                '[C;z1:2]-[O:1]-[S;D4:4](=[O:5])(=[O:6])-[C:7]',
                # Dimethyl sulfate, diethyl sulfate, etc.
                '[C;z1:2]-[O:1]-[S;D4:4](=[O:5])(=[O:6])-[O:7]-[C:8]',
            ],
            'B': [
                # Very generic nitrogen nucleophile, can be amine, amide, hydrazine, hydrazone, N-H-aromatic heterocycle, etc.
                # The only requirement is that nitrogen has at least one hydrogen attached to it
                '[N;h1,h2,h3:3]',
                '[N:3]-[H:9]'
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
