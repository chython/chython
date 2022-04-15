# -*- coding: utf-8 -*-
#
#  Copyright 2019-2022 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from importlib.util import find_spec
from .free_wilson import *
from .functional_groups import *
from .grid import *
from .retro import *


__all__ = ['functional_groups', 'fw_prepare_groups', 'fw_decomposition_tree',
           'grid_depict', 'GridDepict', 'retro_depict', 'RetroDepict']


if find_spec('rdkit'):
    from .rdkit import *
    __all__.extend(['from_rdkit_molecule', 'to_rdkit_molecule'])
