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
"""Registers depiction onto the sealed core containers.

The core owns the eleven method names and their docstrings; this module supplies the bodies, because
`core` may not import a package above it and a `cdef class` cannot be extended from outside.
`_set_depict_fns` is all-or-nothing per group (the layout five, the drawing four, the 3D two).
"""
from ..core._core import _set_depict_fns
from . import figure as _figure
from .layout import molecule as _molecule, reaction as _reaction
from . import x3dom as _x3dom


def register():
    """Point the containers' eleven depiction methods at this package's implementations."""
    _set_depict_fns(clean2d=_molecule.clean2d, layout2d=_molecule.layout2d,
                    rescale2d=_molecule.rescale2d,
                    reaction_clean2d=_reaction.clean2d, reaction_layout2d=_reaction.layout2d,
                    depict=_figure.molecule_depict, scene=_figure.molecule_scene,
                    reaction_depict=_figure.reaction_depict, reaction_scene=_figure.reaction_scene,
                    depict3d=_x3dom.molecule_depict3d, view3d=_x3dom.molecule_view3d)


__all__ = ['register']
