# -*- coding: utf-8 -*-
#
#  Copyright 2021-2025 Ramil Nugmanov <nougmanoff@protonmail.com>
#  Copyright 2021 Alexander Sizov <murkyrussian@gmail.com>
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
from chython import MoleculeContainer
from typing import Tuple, List
from ..algorithms.depict import _render_config, _graph_svg


Tree = Tuple[MoleculeContainer, List['Tree']]


def retro_depict(tree: Tree, *, y_gap=3., x_gap=5., width=None, height=None, clean2d: bool = True) -> str:
    """
    Depict retrosynthetic tree.

    :param tree: Graph of molecules with recursive structure. Each node is tuple of molecule and list of child nodes.
        Child nodes list can be empty.
    :param y_gap: vertical gap between molecules.
    :param x_gap: horizontal gap between molecules.
    :param width: set svg width param. by default auto-calculated.
    :param height: set svg height param. by default auto-calculated.
    :param clean2d: calculate coordinates if necessary.
    """
    font_size = _render_config['font_size']
    font125 = 1.25 * font_size

    arrows = []
    columns = [[tree[0]]]
    current_layer = [(0, x) for x in tree[1]]
    while current_layer:
        next_layer = []
        column = []
        columns.append(column)
        for i, (j, (m, ms)) in enumerate(current_layer):
            column.append(m)
            arrows.append(j)  # i-th molecule in current column connected to j-th molecule in previous.
            next_layer.extend((i, x) for x in ms)
        current_layer = next_layer

    x_shift = 0.
    c_max_x = 0.
    c_max_y = 0.
    arrows = iter(arrows)
    render = []
    last_layer = []
    arrows_coords = []
    for column in columns:
        current_layer = []

        if clean2d:
            for m in column:
                if len(m) > 1:
                    min_x = min(a.x for _, a in m.atoms())
                    max_x = max(a.x for _, a in m.atoms())
                    min_y = min(a.y for _, a in m.atoms())
                    max_y = max(a.y for _, a in m.atoms())
                    if max_y - min_y < .01 and max_x - min_x < 0.01:
                        m.clean2d()

        heights = [max(a.y for _, a in m.atoms()) - min(a.y for _, a in m.atoms()) for m in column]
        y_shift = sum(heights) + y_gap * (len(heights) - 1)  # column height with gaps
        if y_shift > c_max_y:
            c_max_y = y_shift
        y_shift /= 2.  # center align

        for m, h in zip(column, heights):
            plane = [a.xy for _, a in m.atoms()]  # backup
            mx = m._fix_plane_min(x_shift, -y_shift)
            if mx > c_max_x:
                c_max_x = mx

            current_layer.append((mx + 1., y_shift - h / 2.))
            if x_shift:  # except first column
                arrows_coords.append((*last_layer[next(arrows)], x_shift - 1., y_shift - h / 2.))
                y_shift -= h + y_gap

            render.append(m.depict(_embedding=True)[:5])
            for (_, a), xy in zip(m.atoms(), plane):  # restore
                a.xy = xy

        x_shift = c_max_x + x_gap  # between columns gap
        last_layer = current_layer

    _width = c_max_x + 3.0 * font_size
    _height = c_max_y + 2.5 * font_size
    box_y = _height / 2.
    if width is None:
        width = f'{_width:.2f}cm'
    if height is None:
        height = f'{_height:.2f}cm'
    svg = [f'<svg width="{width}" height="{height}" '
           f'viewBox="{-font125:.2f} {-box_y:.2f} {_width:.2f} {_height:.2f}" '
           'xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" version="1.1">',
           '  <defs>\n    <marker id="arrow" markerWidth="10" markerHeight="10" '
           'refX="0" refY="3" orient="auto">\n      <path d="M0,0 L0,6 L9,3"/>\n    </marker>\n  </defs>']
    for atoms, bonds, define, masks, uid in render:
        svg.extend(_graph_svg(atoms, bonds, define, masks, uid, -font125, -box_y, _width, _height))

    svg.append('  <g fill="none" stroke="black" stroke-width=".04" marker-end="url(#arrow)">')
    for x1, y1, x2, y2 in arrows_coords:
        svg.append(f'    <line x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" y2="{y2:.2f}"/>')
    svg.append('  </g>')
    svg.append('</svg>')
    return '\n'.join(svg)


class RetroDepict:
    """
    Grid depict for Jupyter notebooks.
    """
    def __init__(self, tree: Tree):
        """
        :param tree: Graph of molecules with recursive structure.
            Each node is tuple of molecule and list of child nodes.
            Child nodes list can be empty.
        """
        self.tree = tree

    def _repr_svg_(self):
        return retro_depict(self.tree)


__all__ = ['retro_depict', 'RetroDepict']
