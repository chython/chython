# -*- coding: utf-8 -*-
#
#  Copyright 2018-2026 Ramil Nugmanov <nougmanoff@protonmail.com>
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
"""2D layout, vector rendering and the overlay system.

A structure becomes a `Scene` -- a resolution-free tree of boxes, paths and text runs -- which a
backend turns into bytes; `render/svg.py` writes SVG and SVGZ, and `overlay.py` composes scalar data
onto the same scene.  `x3dom.py` is the 3D side and shares none of it -- a stored conformer straight to
X3DOM.  Importing this package runs `_hooks.register()`, which is what makes the eleven depiction
methods (`mol.clean2d()`, `mol.depict()`, `rxn.scene()`, `mol.view3d()`, ...) exist on the containers.
"""
from ._config import Clean2DEngine, cpk, get_clean2d_engine, set_clean2d_engine
from ._hooks import register as _register
from .colormap import NAMED_COLORMAPS, Colormap
from .field import ScalarField
from .figure import molecule_depict, molecule_scene, reaction_depict, reaction_scene
from .layout.molecule import clean2d, layout2d, rescale2d
from .layout.reaction import clean2d as reaction_clean2d, layout2d as reaction_layout2d
from .overlay import AtomField, AtomHalo, BondScale, Highlight, ValueLabels
from .scene import Box, Group, Path, Scene, Text, TextRun
from .style import DepictStyle, get_depict_style, set_depict_style
from .x3dom import JupyterWidget, molecule_depict3d, molecule_view3d


# `clean2d_engine` is deliberately absent: it is rebindable, and a copy taken at import time here
# would go stale on the first `set_clean2d_engine()`.  Read it through `get_clean2d_engine()`, or as
# `chython.clean2d_engine`.
__all__ = ['cpk', 'Clean2DEngine', 'get_clean2d_engine', 'set_clean2d_engine',
           'DepictStyle', 'get_depict_style', 'set_depict_style',
           'molecule_scene', 'molecule_depict', 'reaction_scene', 'reaction_depict',
           'clean2d', 'layout2d', 'rescale2d', 'reaction_clean2d', 'reaction_layout2d',
           'Scene', 'Group', 'Path', 'Text', 'TextRun', 'Box',
           'Highlight', 'AtomHalo', 'AtomField', 'BondScale', 'ValueLabels',
           'ScalarField', 'Colormap', 'NAMED_COLORMAPS',
           'molecule_depict3d', 'molecule_view3d', 'JupyterWidget']

_register()
