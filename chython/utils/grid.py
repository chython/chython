# -*- coding: utf-8 -*-
#
#  Copyright 2021-2023 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from typing import List, Optional
from ..containers import MoleculeContainer
from ..algorithms.depict import _render_config, _graph_svg


def grid_depict(molecules: List[MoleculeContainer], labels: Optional[List[str]] = None, *, cols: int = 3,
                width=None, height=None, clean2d: bool = True):
    """
    Depict molecules grid.

    :param molecules: list of molecules
    :param labels: optional list of text labels
    :param cols: number of molecules per row
    :param width: set svg width param. by default auto-calculated.
    :param height: set svg height param. by default auto-calculated.
    :param clean2d: calculate coordinates if necessary.
    """
    font_size = _render_config['font_size']
    symbols_font_style = _render_config['symbols_font_style']
    font125 = 1.25 * font_size
    font75 = .75 * font_size

    planes = []
    render = []
    render_labels = []
    shift_y = 0.
    shift_x = 0.
    if labels is not None:
        assert len(molecules) == len(labels)
        labels = iter(labels)

    if clean2d:
        for m in molecules:
            if len(m) > 1:
                values = m._plane.values()
                min_x = min(x for x, _ in values)
                max_x = max(x for x, _ in values)
                min_y = min(y for _, y in values)
                max_y = max(y for _, y in values)
                if max_y - min_y < .01 and max_x - min_x < 0.01:
                    m.clean2d()

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

        max_x = 0.
        for m in ms:
            if m is None:
                break
            if labels is not None:
                render_labels.append(f'    <text x="{max_x:.2f}" y="{-shift_y:.2f}">{next(labels)}</text>')
                y = shift_y - height / 2. - font125  # blank
            else:
                y = shift_y - height / 2.
            max_x = m._fix_plane_mean(max_x, y) + 4. * font_size
            render.append(m.depict(_embedding=True)[:5])
            if max_x > shift_x:  # get total width
                shift_x = max_x
        shift_y -= height + 4. * font_size

    # restore planes
    for p, m in zip(planes, molecules):
        m._plane = p

    _width = shift_x - 1.5 * font_size
    _height = -shift_y - 1.5 * font_size
    if width is None:
        width = f'{_width:.2f}cm'
    if height is None:
        height = f'{_height:.2f}cm'
    svg = [f'<svg width="{width}" height="{height}" '
           f'viewBox="{-font125:.2f} {-font125:.2f} {_width:.2f} {_height:.2f}" '
           'xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" version="1.1">']
    for atoms, bonds, define, masks, uid in render:
        svg.extend(_graph_svg(atoms, bonds, define, masks, uid, -font125, -font125, _width, _height))
    svg.append(f'  <g font-size="{font75:.2f}" font-family="{symbols_font_style}">')
    svg.extend(render_labels)
    svg.append('  </g>')
    svg.append('</svg>')
    return '\n'.join(svg)


class GridDepict:
    """
    Grid depict for Jupyter notebooks.
    """
    def __init__(self, molecules: List[MoleculeContainer], labels: Optional[List[str]] = None, *, cols: int = 3):
        """
        :param molecules: list of molecules
        :param labels: optional list of text labels
        :param cols: number of molecules per row
        """
        self.molecules = molecules
        self.labels = labels
        self.cols = cols

    def _repr_svg_(self):
        return grid_depict(self.molecules, self.labels, cols=self.cols)


__all__ = ['grid_depict', 'GridDepict']
