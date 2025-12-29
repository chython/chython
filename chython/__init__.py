# -*- coding: utf-8 -*-
#
#  Copyright 2014-2025 Ramil Nugmanov <nougmanoff@protonmail.com>
#  Copyright 2014-2019 Timur Madzhidov tmadzhidov@gmail.com features and API discussion
#  Copyright 2014-2019 Alexandre Varnek <varnek@unistra.fr> base idea of CGR approach
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
from typing import Literal
from .algorithms.depict import depict_settings
from .containers import *
from .files import *
from .reactor import *
from .utils import *


torch_device = 'cpu'  # AAM model device. Change before first `reset_mapping` call!
clean2d_engine: Literal['smilesdrawer', 'rdkit'] = 'smilesdrawer'
conformer_engine: Literal['rdkit', 'cdpkit'] = 'rdkit'


__all__ = []
