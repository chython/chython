# -*- coding: utf-8 -*-
#
#  Copyright 2021 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from itertools import zip_longest
from typing import List, Tuple, Union
from ..containers import MoleculeContainer, QueryContainer, CGRContainer

graphs = Union[MoleculeContainer, QueryContainer, CGRContainer]


def grid_depict(molecules: Union[List[graphs], Tuple[graphs, ...]], cols: int = 3):
    """
    Depict molecules grid.

    :param cols: number of molecules per row.
    """
    config = MoleculeContainer._render_config
    font_size = config['font_size']
    font125 = 1.25 * font_size

    planes = []
    render = []
    shift_y = 2.
    shift_x = 0.
    for ms in zip_longest(*[iter(molecules)] * cols):
        height = 0.
        for m in ms:
            if m is None:
                break
            min_y = min(y for x, y in m._plane.values())
            max_y = max(y for x, y in m._plane.values())
            h = max_y - min_y
            if height < h:  # get height of row
                height = h
            planes.append(m._plane.copy())

        # now align mols by height
        shift_y -= height / 2. + 2.
        max_x = 0.
        for m in ms:
            if m is None:
                break
            max_x = m._fix_plane_mean(max_x, shift_y) + 2.
            render.append(m.depict(embedding=True)[:5])
            if max_x > shift_x:  # get total width
                shift_x = max_x
        shift_y -= height / 2.

    # restore planes
    for p, m in zip(planes, molecules):
        m._plane = p

    width = shift_x + 3.0 * font_size
    height = -shift_y + 2.5 * font_size
    svg = [f'<svg width="{width:.2f}cm" height="{height:.2f}cm" '
           f'viewBox="{-font125:.2f} {-font125:.2f} {width:.2f} '
           f'{height:.2f}" xmlns="http://www.w3.org/2000/svg" version="1.1">']
    for atoms, bonds, define, masks, uid in render:
        svg.extend(MoleculeContainer._graph_svg(atoms, bonds, define, masks, uid, -font125, -font125, width, height))
    svg.append('</svg>')
    return '\n'.join(svg)


__all__ = ['grid_depict']
