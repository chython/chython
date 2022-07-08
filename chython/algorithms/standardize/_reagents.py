# -*- coding: utf-8 -*-
#
#  Copyright 2022 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from lazy_object_proxy import Proxy


def _reagents():
    from ... import smiles

    tmp = ['N', 'O', 'F', 'Cl', 'Br', 'I', 'O=C=O',  # inorganic
           'C1=CC=CC=C1', 'C1CCCCC1', 'CC1=CC=CC=C1', 'CCCCCC', 'CCCCCCC',  # hydrocarbon
           'CO', 'CCO', 'CC(C)O',  # alcohol
           'OC=O', 'CC(=O)O',  # acid
           'CC(=O)OCC',  # ester
           'CCOCC', 'C1COCC1',  # ether
           ]

    rules = set()
    for x in tmp:
        x = smiles(x)
        x.thiele()
        rules.add(x)
    return rules


reagents_set = Proxy(_reagents)


__all__ = ['reagents_set']
