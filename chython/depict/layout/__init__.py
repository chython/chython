# -*- coding: utf-8 -*-
#
#  Copyright 2019-2026 Ramil Nugmanov <nougmanoff@protonmail.com>
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
"""2D layout, one module per container.

The two modules are exported and their contents are not: both define `clean2d` and `layout2d`, so a
star import would let one silently shadow the other.  Write `layout.molecule.clean2d(mol)`, or call the
method the core container carries.
"""
from . import molecule, reaction


__all__ = ['molecule', 'reaction']
