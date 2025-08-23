# -*- coding: utf-8 -*-
#
#  Copyright 2017-2025 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from typing import Union
from zlib import decompress
from .bonds import *
from .cgr import *
from .molecule import *
from .query import *
from .reaction import *


def unpach(data: bytes, /, *, compressed=True, skip_labels_calculation=False) -> Union[MoleculeContainer, ReactionContainer]:
    if compressed:
        data = decompress(data)
    try:
        return MoleculeContainer.unpack(data, compressed=False, skip_labels_calculation=skip_labels_calculation)
    except ValueError:
        pass
    # second try
    return ReactionContainer.unpack(data, compressed=False)


def from_rdkit(mol, /):
    return MoleculeContainer.from_rdkit(mol)


unpack = unpach
from_rdkit_molecule = from_rdkit


__all__ = [x for x in locals() if x.endswith('Container')]
__all__.extend(['Bond', 'QueryBond', 'unpack', 'unpach', 'from_rdkit', 'from_rdkit_molecule'])
