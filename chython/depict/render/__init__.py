# -*- coding: utf-8 -*-
#
#  Copyright 2026 Ramil Nugmanov <nougmanoff@protonmail.com>
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
"""Backends: one module per output format, each turning a `Scene` into bytes or a string.

A backend applies the two device facts -- the y-flip and the unit scale -- and nothing else; geometry
arrives finished, in molecule coordinates, y-up.  `Scene.to_svg()` and friends import these lazily, so
`scene.py` and `render/` do not cycle.
"""
from .svg import to_svg, to_svgz


__all__ = ['to_svg', 'to_svgz']
