# -*- coding: utf-8 -*-
#
#  Copyright 2014-2021 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from .containers import *
from .files import *
from .reactor import *
from .utils import *


smiles = SMILESRead.create_parser(ignore=True, remap=False)
xyz = XYZRead.create_parser()


pickle_cache = False  # store cached attributes in pickle


__all__ = ['smiles', 'xyz', 'mdl_mol']

if 'INCHIRead' in locals():
    inchi = INCHIRead.create_parser()
    __all__.append('inchi')
