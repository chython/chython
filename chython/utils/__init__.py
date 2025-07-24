# -*- coding: utf-8 -*-
#
#  Copyright 2019-2025 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from .free_wilson import *
from .grid import *
from .retro import *
from ..containers.graph import Graph


_patched = False


def patch_pandas():
    """
    Fix pandas molecules representation.
    """
    # prevent recursive patching
    global _patched
    if _patched:
        return
    _patched = True

    from pandas.io.formats import printing
    from pandas.io.formats.printing import is_sequence

    def w(obj):
        if isinstance(obj, Graph):
            return False
        return is_sequence(obj)

    printing.is_sequence = w


__all__ = ['fw_prepare_groups', 'fw_decomposition_tree',
           'grid_depict', 'GridDepict', 'retro_depict', 'RetroDepict', 'patch_pandas']
