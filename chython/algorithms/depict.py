# -*- coding: utf-8 -*-
#
#  Copyright 2018-2025 Ramil Nugmanov <nougmanoff@protonmail.com>
#  Copyright 2019-2020 Dinar Batyrshin <batyrshin-dinar@mail.ru>
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
from asyncio import new_event_loop
from chython._cache import cached_method
from collections import defaultdict
from functools import partial
from math import atan2, sin, cos, hypot, pi
from os.path import join
from tempfile import TemporaryDirectory
from typing import TYPE_CHECKING, Union, Literal
from uuid import uuid4
from zlib import compress


if TYPE_CHECKING:
    from chython import ReactionContainer, MoleculeContainer

cpk = tuple('''
 #909090                                                                                         #D9FFFF
 #CC80FF #C2FF00                                         #FFB5B5 #101010 #3050F8 #FF0D0D #90E050 #B3E3F5
 #AB5CF2 #8AFF00                                         #BFA6A6 #F0C8A0 #FF8000 #C6C600 #1FF01F #80D1E3
 #8F40D4 #3DFF00 #E6E6E6 #BFC2C7 #A6A6AB #8A99C7 #9C7AC7 
                 #E06633 #F090A0 #50D050 #C88033 #7D80B0 #C28F8F #668F8F #BD80E3 #FFA100 #A62929 #5CB8D1
 #702EB0 #00FF00 #94FFFF #94E0E0 #73C2C9 #54B5B5 #3B9E9E 
                 #248F8F #0A7D8C #006985 #C0C0C0 #FFD98F #A67573 #668080 #9E63B5 #D47A00 #940094 #429EB0
 #57178F #00C900 #70D4FF 
                 #FFFFC7 #D9FFC7 #C7FFC7 #A3FFC7 #8FFFC7 #61FFC7 #45FFC7
                 #30FFC7 #1FFFC7 #00FF9C #00E675 #00D452 #00BF38 #00AB24 
                         #4DC2FF #4DA6FF #2194D6 #267DAB
                 #266696 #175487 #D0D0E0 #FFD123 #B8B8D0 #A6544D #575961 #9E4FB5 #AB5C00 #754F45 #428296
 #420066 #007D00 #70ABFA
                 #00BAFF #00A1FF #008FFF #0080FF #006BFF #545CF2 #785CE3
                 #8A4FE3 #A136D4 #B31FD4 #B31FBA #B30DA6 #BD0D87 #C70066
                         #CC0059 #D1004F #D90045 #E00038 
                 #E6002E #EB0026 #EB0026 #EB0026 #EB0026 #EB0026 #EB0026 #EB0026 #EB0026 #EB0026 #EB0026
'''.split())
_render_charge = {-4: '4-', -3: '3-', -2: '2-', -1: '-', 1: '+', 2: '2+', 3: '3+', 4: '4+'}
_render_config = {'carbon': False, 'dashes': (.2, .1), 'span_dy': .15, 'mapping': True, 'font_size': .5,
                  'span_size': .35, 'other_size': 0.3, 'monochrome': False, 'bond_color': 'black', 'bond_width': .04,
                  'other_color': 'black', 'bond_radius': .02, 'atom_radius': -.2, 'mapping_size': .25,
                  'atoms_colors': cpk, 'triple_space': .13, 'double_space': .06, 'mapping_color': '#0305A7',
                  'aromatic_space': .14, 'aromatic_dashes': (.15, .05), 'dx_m': .05, 'dy_m': .2,
                  'other_font_style': 'monospace', 'dx_ci': .05, 'dy_ci': 0.2, 'symbols_font_style': 'sans-serif',
                  'mapping_font_style': 'monospace', 'wedge_space': .08, 'arrow_color': 'black',
                  'query_color': '#5D8AA8', 'broken_color': 'red', 'formed_color': 'green',
                  'cgr_aromatic_space': .18, 'dx_nh': .075, 'dy_nh': .15,
                  'dynamic_bond_width': .08}

loop = browser = None


async def svg_render(s, t, width, height, scale):
    page = await browser.newPage()
    await page.setViewport({'deviceScaleFactor': scale, 'width': width, 'height': height})
    await page.goto(f'file://{s}')
    element = await page.querySelector('svg')
    await element.screenshot({'path': t, 'omitBackground': True})
    await page.close()


def svg2png(svg: str, width: int = 1000, height: int = 1000, scale: float = 1.):
    global loop, browser

    if loop is None:  # lazy browser launcher
        from pyppeteer import launch

        loop = new_event_loop()
        browser = loop.run_until_complete(launch())
    elif browser is None:
        raise ImportError('pyppeteer initialization failed')

    with TemporaryDirectory() as tmpdir:
        with open(s := join(tmpdir, 'input.svg'), 'w') as f:
            f.write(svg)

        loop.run_until_complete(svg_render(s, (t := join(tmpdir, 'output.png')), width, height, scale))

        with open(t, 'rb') as f:
            return f.read()


def _rotate_vector(x1, y1, x2, y2):
    """
    rotate x,y vector over x2-x1, y2-y1 angle
    """
    angle = atan2(y2, x2)
    cos_rad = cos(angle)
    sin_rad = sin(angle)
    return cos_rad * x1 - sin_rad * y1, sin_rad * x1 + cos_rad * y1


rotate_vector = _rotate_vector


def _optimal_label_direction(ax, ay, neighbor_xys):
    """Find the direction with the largest angular gap between bonds.

    Returns angle in radians pointing to the most open space around an atom.
    Coordinates should be in SVG space (y-inverted).
    """
    if not neighbor_xys:
        return pi / 4  # default: down-right in SVG
    angles = sorted(atan2(ny - ay, nx - ax) for nx, ny in neighbor_xys)
    n = len(angles)
    if n == 1:
        return angles[0] + pi  # opposite the single bond
    best_gap = 0
    best_angle = angles[0] + pi
    for i in range(n):
        if i == n - 1:
            gap = angles[0] + 2 * pi - angles[-1]
        else:
            gap = angles[i + 1] - angles[i]
        if gap > best_gap:
            best_gap = gap
            best_angle = angles[i] + gap / 2
    return best_angle


def _graph_svg(atoms, bonds, define, masks, uid, viewbox_x, viewbox_y, width, height):
    svg = [f'  <g id="{uid}-molecule">\n    <defs>']
    svg.extend(define)
    if bonds:
        if masks:
            svg.append(f'      <mask id="{uid}-mask">\n'
                       f'        <rect x="{viewbox_x:.2f}" y="{viewbox_y:.2f}" '
                       f'width="{width:.2f}" height="{height:.2f}" fill="white"/>')
            svg.extend(masks)
            svg.append('      </mask>\n    </defs>\n'
                       f'    <g fill="none" stroke="{_render_config["bond_color"]}" '
                       f'stroke-width="{_render_config["bond_width"]:.2f}" mask="url(#{uid}-mask)">')
            if len(bonds) == 1:  # SVG BUG adhoc
                svg.append(f'      <line x1="{viewbox_x:.2f}" y1="{viewbox_y:.2f}" '
                           f'x2="{viewbox_x + width:.2f}" y2="{viewbox_y:.2f}" stroke="none"/>')
        else:
            svg.append(f'    </defs>\n    <g fill="none" stroke="{_render_config["bond_color"]}" '
                       f'stroke-width="{_render_config["bond_width"]:.2f}">')
        svg.extend(bonds)
        svg.append('    </g>')
    else:
        svg.append('    </defs>')

    svg.extend(atoms)
    svg.append('  </g>')
    return svg


def _render_aromatic_bond(n_x, n_y, m_x, m_y, c_x, c_y):
    aromatic_space = _render_config['aromatic_space']
    dash3, dash4 = _render_config['aromatic_dashes']
    # n aligned xy
    mn_x, mn_y, cn_x, cn_y = m_x - n_x, m_y - n_y, c_x - n_x, c_y - n_y

    # nm reoriented xy
    mr_x, mr_y = hypot(mn_x, mn_y), 0
    cr_x, cr_y = _rotate_vector(cn_x, cn_y, mn_x, -mn_y)

    if cr_y and aromatic_space / cr_y < .65:
        if cr_y > 0:
            r_y = aromatic_space
        else:
            r_y = -aromatic_space
            cr_y = -cr_y

        ar_x = aromatic_space * cr_x / cr_y
        br_x = mr_x - aromatic_space * (mr_x - cr_x) / cr_y

        # backward reorienting
        an_x, an_y = _rotate_vector(ar_x, r_y, mn_x, mn_y)
        bn_x, bn_y = _rotate_vector(br_x, r_y, mn_x, mn_y)
        a_x, a_y = n_x + an_x, n_y + an_y
        b_x, b_y = n_x + bn_x, n_y + bn_y

        return f'      <line x1="{a_x:.2f}" y1="{-a_y:.2f}" x2="{b_x:.2f}" y2="{-b_y:.2f}" ' \
               f'stroke-dasharray="{dash3:.2f} {dash4:.2f}"/>'


def depict_settings(*, carbon: bool | None = None, aam: bool | None = None, monochrome: bool | None = None,
                    bond_color: str | None = None, aam_color: str | None = None, atoms_colors: tuple | None = None,
                    bond_width: float | None = None, wedge_space: float | None = None,
                    dashes: tuple[float, float] | None = None,
                    aromatic_dashes: tuple[float, float] | None = None,
                    dx_ci: float | None = None, dy_ci: float | None = None,
                    dx_m: float | None = None, dy_m: float | None = None,
                    span_dy: float | None = None, double_space: float | None = None,
                    triple_space: float | None = None, aromatic_space: float | None = None,
                    atom_radius: float | None = None, bond_radius: float | None = None,
                    font_size: float | None = None, other_size: float | None = None,
                    span_size: float | None = None, aam_size: float | None = None,
                    symbols_font_style: str | None = None, other_font_style: str | None = None,
                    other_color: str | None = None, arrow_color: str | None = None,
                    mapping_font_style: str | None = None,
                    query_color: str | None = None, broken_color: str | None = None,
                    formed_color: str | None = None, cgr_aromatic_space: float | None = None,
                    dx_nh: float | None = None, dy_nh: float | None = None):
    """
    Settings for depict of chemical structures.
    Pass parameters as keyword arguments to override current values. Unspecified parameters are left unchanged.

    :param carbon: if True, depict atom C
    :param font_size: font size. When changed, font-derived sizes (span_dy, span_size, other_size, mapping_size,
        dx_m, dy_m, dx_ci, dy_ci, dx_nh, dy_nh) are recalculated proportionally unless explicitly overridden
    :param aam_size: atom-to-atom mapping font size
    :param span_size: font size for hydrogen count
    :param other_size: isotope, radical, charges, neighbors and hybridization symbols size
    :param bond_width: bond width
    :param bond_color: color of bonds
    :param aam_color: atom-to-atom mapping color
    :param atoms_colors: atom colors where key is atomic number - 1, value is atom color (str)
    :param other_color: color for charges, radicals, isotopes
    :param arrow_color: color of reaction arrow and plus signs
    :param symbols_font_style: font style for atom symbols
    :param other_font_style: font style for charges, radicals, isotopes, hybridization and neighbors
    :param aam: if True, depict mapping
    :param monochrome: if True, colors of items in molecule not used
    :param dashes: first value is long of visible line, second is long of invisible line
    :param aromatic_space: space between simple and aromatic bonds
    :param triple_space: space between simple and triple bonds
    :param double_space: space between simple and double bonds
    :param aromatic_dashes: first value is long of visible line, second is long of invisible line
    :param atom_radius: radius of atoms spheres. if negative is multiplier to covalent radii
    :param bond_radius: radius of bonds spheres
    :param dx_ci: x-axis offset relative to the center of the atom symbol for radical, charges, isotope
    :param dy_ci: y-axis offset relative to the center of the atom symbol for radical, charges, isotope
    :param dx_m: x-axis offset relative to the center of the atom symbol for atom-to-atom mapping
    :param dy_m: y-axis offset relative to the center of the atom symbol for atom-to-atom mapping
    :param span_dy: y-axis offset relative to the center of the atom symbol for hydrogen count
    :param mapping_font_style: font style for mapping
    :param wedge_space: wedge bond width
    :param query_color: color for query atom/bond markers
    :param broken_color: color for broken bonds in CGR depiction
    :param formed_color: color for formed bonds in CGR depiction
    :param cgr_aromatic_space: space between aromatic bonds in CGR depiction
    :param dx_nh: x-axis offset for NH count display
    :param dy_nh: y-axis offset for NH count display
    """
    if carbon is not None: _render_config['carbon'] = carbon
    if dashes is not None: _render_config['dashes'] = dashes
    if aam is not None: _render_config['mapping'] = aam
    if monochrome is not None: _render_config['monochrome'] = monochrome
    if bond_color is not None: _render_config['bond_color'] = bond_color
    if bond_width is not None: _render_config['bond_width'] = bond_width
    if other_color is not None: _render_config['other_color'] = other_color
    if arrow_color is not None: _render_config['arrow_color'] = arrow_color
    if bond_radius is not None: _render_config['bond_radius'] = bond_radius
    if atom_radius is not None: _render_config['atom_radius'] = -atom_radius
    if atoms_colors is not None: _render_config['atoms_colors'] = atoms_colors
    if triple_space is not None: _render_config['triple_space'] = triple_space
    if double_space is not None: _render_config['double_space'] = double_space
    if aam_color is not None: _render_config['mapping_color'] = aam_color
    if aromatic_space is not None: _render_config['aromatic_space'] = aromatic_space
    if aromatic_dashes is not None: _render_config['aromatic_dashes'] = aromatic_dashes
    if other_font_style is not None: _render_config['other_font_style'] = other_font_style
    if symbols_font_style is not None: _render_config['symbols_font_style'] = symbols_font_style
    if mapping_font_style is not None: _render_config['mapping_font_style'] = mapping_font_style
    if wedge_space is not None: _render_config['wedge_space'] = wedge_space
    if query_color is not None: _render_config['query_color'] = query_color
    if broken_color is not None: _render_config['broken_color'] = broken_color
    if formed_color is not None: _render_config['formed_color'] = formed_color
    if cgr_aromatic_space is not None: _render_config['cgr_aromatic_space'] = cgr_aromatic_space

    # Handle font_size change: recalculate derived sizes proportionally
    if font_size is not None:
        _render_config['font_size'] = font_size
        if span_dy is None: _render_config['span_dy'] = 0.3 * font_size
        if span_size is None: _render_config['span_size'] = 0.7 * font_size
        if other_size is None: _render_config['other_size'] = 0.6 * font_size
        if aam_size is None: _render_config['mapping_size'] = 0.5 * font_size
        if dx_m is None: _render_config['dx_m'] = 0.1 * font_size
        if dy_m is None: _render_config['dy_m'] = 0.4 * font_size
        if dx_ci is None: _render_config['dx_ci'] = 0.1 * font_size
        if dy_ci is None: _render_config['dy_ci'] = 0.4 * font_size
        if dx_nh is None: _render_config['dx_nh'] = 0.15 * font_size
        if dy_nh is None: _render_config['dy_nh'] = 0.3 * font_size

    # Apply explicit overrides for font-derived values
    if span_dy is not None: _render_config['span_dy'] = span_dy
    if span_size is not None: _render_config['span_size'] = span_size
    if other_size is not None: _render_config['other_size'] = other_size
    if aam_size is not None: _render_config['mapping_size'] = aam_size
    if dx_m is not None: _render_config['dx_m'] = dx_m
    if dy_m is not None: _render_config['dy_m'] = dy_m
    if dx_ci is not None: _render_config['dx_ci'] = dx_ci
    if dy_ci is not None: _render_config['dy_ci'] = dy_ci
    if dx_nh is not None: _render_config['dx_nh'] = dx_nh
    if dy_nh is not None: _render_config['dy_nh'] = dy_nh


class DepictMolecule:
    __slots__ = ()

    def depict(self: Union['MoleculeContainer', 'DepictMolecule'], *, width=None, height=None, clean2d: bool = True,
               format: Literal['svg', 'png', 'svgz'] = 'svg', png_width=1000, png_height=1000, png_scale=1.,
               _embedding=False) -> Union[str, bytes]:
        """
        Depict molecule in SVG or PNG format.

        :param width: set svg width param. by default auto-calculated.
        :param height: set svg height param. by default auto-calculated.
        :param clean2d: calculate coordinates if necessary.
        :param format: output format - svg string, png bytes or gz compressed svg
        :param png_width, png_height: viewport size for PNG rendering
        :param png_scale: image scaling in PNG rendering
        """
        uid = str(uuid4())
        min_x = min(a.x for _, a in self.atoms())
        max_x = max(a.x for _, a in self.atoms())
        min_y = min(a.y for _, a in self.atoms())
        max_y = max(a.y for _, a in self.atoms())
        if clean2d and len(self) > 1 and max_y - min_y < .01 and max_x - min_x < 0.01:
            self.clean2d()
            min_x = min(a.x for _, a in self.atoms())
            max_x = max(a.x for _, a in self.atoms())
            min_y = min(a.y for _, a in self.atoms())
            max_y = max(a.y for _, a in self.atoms())

        bonds = self.__render_bonds()
        atoms, define, masks = self.__render_atoms(uid)
        if _embedding:
            return atoms, bonds, define, masks, uid, min_x, min_y, max_x, max_y

        font_size = _render_config['font_size']
        font125 = 1.25 * font_size
        _width = max_x - min_x + 4.0 * font_size
        _height = max_y - min_y + 2.5 * font_size
        viewbox_x = min_x - font125
        viewbox_y = -max_y - font125

        if width is None:
            width = f'{_width:.2f}cm'
        if height is None:
            height = f'{_height:.2f}cm'

        svg = [f'<svg width="{width}" height="{height}" '
               f'viewBox="{viewbox_x:.2f} {viewbox_y:.2f} {_width:.2f} {_height:.2f}" '
               'xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" version="1.1">']
        svg.extend(_graph_svg(atoms, bonds, define, masks, uid, viewbox_x, viewbox_y, _width, _height))
        svg.append('</svg>')
        svg = '\n'.join(svg)
        if format == 'svg':
            return svg
        elif format == 'png':
            return svg2png(svg, png_width, png_height, png_scale)
        elif format == 'svgz':
            return compress(svg.encode(), 9)
        raise ValueError(f'format must be svg, png or svgz, not {format}')

    @cached_method
    def _repr_svg_(self):
        return self.depict()

    def __render_bonds(self: Union['MoleculeContainer', 'DepictMolecule']):
        atoms = self._atoms
        svg = []
        double_space = _render_config['double_space']
        triple_space = _render_config['triple_space']
        wedge_space = _render_config['wedge_space']
        dash1, dash2 = _render_config['dashes']
        color = f' fill="{_render_config["bond_color"]}"'

        wedge = defaultdict(set)
        for n, m, s in self._wedge_map:
            wedge[n].add(m)
            wedge[m].add(n)

            nx, ny = atoms[n].xy
            mx, my = atoms[m].xy
            ny, my = -ny, -my
            dx, dy = _rotate_vector(0, wedge_space, mx - nx, ny - my)

            svg.append(f'      <path d="M{nx:.2f} {ny:.2f} L{mx + dx:.2f} {my + dy:.2f} '
                       f'L{mx - dx:.2f} {my - dy:.2f} Z"{s == 1 and color or ""}/>')

        for n, m, bond in self.bonds():
            if m in wedge[n]:
                continue
            nx, ny = atoms[n].xy
            mx, my = atoms[m].xy
            ny, my = -ny, -my
            if bond in (1, 4):
                svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"/>')
            elif bond == 2:
                dx, dy = _rotate_vector(0, double_space, mx - nx, ny - my)
                svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" x2="{mx + dx:.2f}" y2="{my - dy:.2f}"/>')
                svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" x2="{mx - dx:.2f}" y2="{my + dy:.2f}"/>')
            elif bond == 3:
                dx, dy = _rotate_vector(0, triple_space, mx - nx, ny - my)
                svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" x2="{mx + dx:.2f}" y2="{my - dy:.2f}"/>')
                svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"/>')
                svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" x2="{mx - dx:.2f}" y2="{my + dy:.2f}"/>')
            else:
                svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}" '
                           f'stroke-dasharray="{dash1:.2f} {dash2:.2f}"/>')

        for ring in self.aromatic_rings:
            cx = sum(atoms[n].x for n in ring) / len(ring)
            cy = sum(atoms[n].y for n in ring) / len(ring)

            for n, m in zip(ring, ring[1:]):
                nx, ny = atoms[n].xy
                mx, my = atoms[m].xy
                aromatic = _render_aromatic_bond(nx, ny, mx, my, cx, cy)
                if aromatic:
                    svg.append(aromatic)

            nx, ny = atoms[ring[-1]].xy
            mx, my = atoms[ring[0]].xy
            aromatic = _render_aromatic_bond(nx, ny, mx, my, cx, cy)
            if aromatic:
                svg.append(aromatic)
        return svg

    def __render_atoms(self: 'MoleculeContainer', uid):
        bonds = self._bonds

        carbon = _render_config['carbon']
        mapping = _render_config['mapping']
        span_size = _render_config['span_size']
        font_size = _render_config['font_size']
        monochrome = _render_config['monochrome']
        other_size = _render_config['other_size']
        atoms_colors = _render_config['atoms_colors']
        mapping_size = _render_config['mapping_size']
        dx_m, dy_m = _render_config['dx_m'], _render_config['dy_m']
        dx_ci, dy_ci = _render_config['dx_ci'], _render_config['dy_ci']
        symbols_font_style = _render_config['symbols_font_style']
        span_dy = _render_config['span_dy']
        other_font_style = _render_config['other_font_style']
        mapping_font_style = _render_config['mapping_font_style']

        if monochrome:
            map_fill = other_fill = 'black'
        else:
            map_fill = _render_config['mapping_color']
            other_fill = _render_config['other_color']

        font2 = .2 * font_size
        font3 = .3 * font_size
        font4 = .4 * font_size
        font5 = .5 * font_size
        font6 = .6 * font_size
        font7 = .7 * font_size
        font15 = .15 * font_size
        font25 = .25 * font_size
        stroke_width_s = font_size * .1
        stroke_width_o = other_size * .1
        stroke_width_m = mapping_size * .1

        svg = []
        maps = []
        symbols = []
        fill_zone = []
        others = []
        define = []
        mask = []

        for n, atom in self.atoms():
            x, y = atom.x, -atom.y
            symbol = atom.atomic_symbol
            if (symbol != 'C' or atom.charge or atom.is_radical or atom.isotope or carbon
                    or not bonds[n] or sum(b == 2 for b in bonds[n].values()) == 2):
                if atom.charge:
                    others.append(f'        <text x="{x:.2f}" y="{y:.2f}" dx="{dx_ci:.2f}" dy="-{dy_ci:.2f}">'
                                  f'{_render_charge[atom.charge]}{"↑" if atom.is_radical else ""}</text>')
                elif atom.is_radical:
                    others.append(f'        <text x="{x:.2f}" y="{y:.2f}" dx="{dx_ci:.2f}" dy="-{dy_ci:.2f}">↑</text>')
                if atom.isotope:
                    others.append(f'        <text x="{x:.2f}" y="{y:.2f}" dx="-{dx_ci:.2f}" dy="-{dy_ci:.2f}" '
                                  f'text-anchor="end">{atom.isotope}</text>')

                if len(symbol) > 1:
                    dx = font7
                    dx_mm = dx_m + font5
                    if symbol[-1] in ('l', 'i', 'r', 't'):
                        rx = font6
                        ax = font25
                    else:
                        rx = font7
                        ax = font15
                    fill_zone.append(f'          <ellipse cx="{x - ax:.2f}" cy="{y:.2f}" rx="{rx}" ry="{font4}"/>')
                else:
                    if symbol == 'I':
                        dx = font15
                        dx_mm = dx_m
                    else:
                        dx = font4
                        dx_mm = dx_m + font2
                    fill_zone.append(f'          <circle cx="{x:.2f}" cy="{y:.2f}" r="{font4:.2f}"/>')

                h = atom.implicit_hydrogens
                if h == 1:
                    h = 'H'
                elif h:
                    h = f'H<tspan  dy="{span_dy:.2f}" font-size="{span_size:.2f}">{h}</tspan>'
                else:
                    h = ''
                symbols.append(f'        <text id="{uid}-{n}" x="{x:.2f}" y="{y:.2f}" dx="-{dx:.2f}" dy="{font4:.2f}">'
                               f'{symbol}{h}</text>')

                svg.append(f'      <use xlink:href="#{uid}-{n}" '
                           f'fill="{"black" if monochrome else atoms_colors[atom.atomic_number - 1]}"/>')

                if mapping:
                    maps.append(f'        <text x="{x:.2f}" y="{y:.2f}" dx="-{dx_mm:.2f}" '
                                f'dy="{dy_m + font3:.2f}">{n}</text>')
            elif mapping:
                maps.append(f'        <text x="{x:.2f}" y="{y:.2f}" dx="-{dx_m:.2f}" dy="{dy_m:.2f}">{n}</text>')

        if svg:  # group atoms symbols
            if fill_zone:
                mask.append('        <g>')
                mask.extend(fill_zone)
                mask.append('        </g>')

            svg.insert(0, f'    <g font-size="{font_size:.2f}" font-family="{symbols_font_style}">')
            svg.append('    </g>')
            define.append(f'      <g id="{uid}-symbols" font-size="{font_size:.2f}" '
                          f'font-family="{symbols_font_style}">')
            define.extend(symbols)
            define.append('      </g>')
            mask.append('        <g stroke="black">\n          '
                        f'<use xlink:href="#{uid}-symbols" stroke-width="{stroke_width_s:.2f}"/>')

            if others:
                define.append(f'      <g id="{uid}-attrs" font-size="{other_size:.2f}" '
                              f'font-family="{other_font_style}">')
                define.extend(others)
                define.append('      </g>')
                svg.append(f'    <use xlink:href="#{uid}-attrs" fill="{other_fill}"/>')
                mask.append(f'          <use xlink:href="#{uid}-attrs" stroke-width="{stroke_width_o:.2f}"/>')
        if maps:
            if not svg:  # no atoms but maps
                mask.append('        <g stroke="black">')
            define.append(f'      <g id="{uid}-mapping" font-size="{mapping_size:.2f}" '
                          f'font-family="{mapping_font_style}" text-anchor="end">')
            define.extend(maps)
            define.append('      </g>')
            svg.append(f'    <use xlink:href="#{uid}-mapping" fill="{map_fill}"/>')
            mask.append(f'          <use xlink:href="#{uid}-mapping" stroke-width="{stroke_width_m:.2f}"/>\n'
                        '        </g>')
        elif svg:  # no maps but atoms
            mask.append('        </g>')
        return svg, define, mask


class DepictReaction:
    __slots__ = ()

    def depict(self: 'ReactionContainer', *, width=None, height=None, clean2d: bool = True,
               format: Literal['svg', 'png', 'svgz'] = 'svg',
               png_width=1000, png_height=1000, png_scale=1.) -> Union[str, bytes]:
        """
        Depict reaction in SVG format.

        :param width: set svg width param. by default auto-calculated.
        :param height: set svg height param. by default auto-calculated.
        :param clean2d: calculate coordinates if necessary.
        :param format: output format - svg string, png bytes or gz compressed svg
        :param png_width, png_height: viewport size for PNG rendering
        :param png_scale: image scaling in PNG rendering
        """
        arrow_color = _render_config['arrow_color']
        if not self._arrow:
            if clean2d:
                for m in self.molecules():
                    if len(m) > 1:
                        try:
                            min_x = min(a.x for _, a in m.atoms())
                            max_x = max(a.x for _, a in m.atoms())
                            min_y = min(a.y for _, a in m.atoms())
                            max_y = max(a.y for _, a in m.atoms())
                        except AttributeError:
                            m.clean2d()
                        else:
                            if max_y - min_y < .01 and max_x - min_x < 0.01:
                                m.clean2d()
            self.fix_positions()

        r_atoms = []
        r_bonds = []
        r_defines = []
        r_masks = []
        r_uids = []
        r_max_x = r_max_y = r_min_y = 0
        for m in self.molecules():
            atoms, bonds, define, masks, uid, min_x, min_y, max_x, max_y = m.depict(clean2d=False, _embedding=True)
            r_atoms.append(atoms)
            r_bonds.append(bonds)
            r_defines.append(define)
            r_masks.append(masks)
            r_uids.append(uid)
            if max_x > r_max_x:
                r_max_x = max_x
            if max_y > r_max_y:
                r_max_y = max_y
            if min_y < r_min_y:
                r_min_y = min_y

        font_size = _render_config['font_size']
        font125 = 1.25 * font_size
        _width = r_max_x + 4.0 * font_size
        _height = r_max_y - r_min_y + 2.5 * font_size
        viewbox_x = -font125
        viewbox_y = -r_max_y - font125

        if width is None:
            width = f'{_width:.2f}cm'
        if height is None:
            height = f'{_height:.2f}cm'

        svg = [f'<svg width="{width}" height="{height}" '
               f'viewBox="{viewbox_x:.2f} {viewbox_y:.2f} {_width:.2f} {_height:.2f}" '
               'xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" version="1.1">\n'
               '  <defs>\n    <marker id="arrow" markerWidth="10" markerHeight="10" '
               f'refX="0" refY="3" orient="auto">\n      <path d="M0,0 L0,6 L9,3" fill="{arrow_color}"/>\n'
               '    </marker>\n  </defs>\n'
               f'  <line x1="{self._arrow[0]:.2f}" y1="0" x2="{self._arrow[1]:.2f}" y2="0" '
               f'fill="none" stroke="{arrow_color}" stroke-width=".04" marker-end="url(#arrow)"/>']

        sings_plus = self._signs
        if sings_plus:
            svg.append(f'  <g fill="none" stroke="{arrow_color}" stroke-width=".04">')
            for x in sings_plus:
                svg.append(f'    <line x1="{x + .35:.2f}" y1="0" x2="{x + .65:.2f}" y2="0"/>')
                svg.append(f'    <line x1="{x + .5:.2f}" y1="0.15" x2="{x + .5:.2f}" y2="-0.15"/>')
            svg.append('  </g>')

        for atoms, bonds, define, masks, uid in zip(r_atoms, r_bonds, r_defines, r_masks, r_uids):
            svg.extend(_graph_svg(atoms, bonds, define, masks, uid, viewbox_x, viewbox_y, _width, _height))
        svg.append('</svg>')
        svg = '\n'.join(svg)
        if format == 'svg':
            return svg
        elif format == 'png':
            return svg2png(svg, png_width, png_height, png_scale)
        elif format == 'svgz':
            return compress(svg.encode(), 9)
        raise ValueError(f'format must be svg, png or svgz, not {format}')

    @cached_method
    def _repr_svg_(self):
        return self.depict()


class Depict:
    __slots__ = ()
    _extra_padding = 0.  # extra uniform padding for subclasses

    @cached_method
    def _repr_svg_(self):
        return self.depict()

    def depict(self, *, width=None, height=None, clean2d: bool = True,
               format: Literal['svg', 'png', 'svgz'] = 'svg', png_width=1000, png_height=1000, png_scale=1.,
               _embedding=False) -> Union[str, bytes]:
        """
        Depict structure in SVG, PNG or SVGZ format.

        :param width: set svg width param. by default auto-calculated.
        :param height: set svg height param. by default auto-calculated.
        :param clean2d: calculate coordinates if necessary.
        :param format: output format - svg string, png bytes or gz compressed svg
        :param png_width, png_height: viewport size for PNG rendering
        :param png_scale: image scaling in PNG rendering
        """
        if not self._plane or len(self._plane) < len(self._atoms):
            if clean2d:
                try:
                    self.clean2d()
                except (ImportError, Exception):
                    pass
            for n, atom in self._atoms.items():
                if n not in self._plane:
                    self._plane[n] = (0., 0.)
                    if hasattr(atom, 'xy'):
                        atom.xy = (0., 0.)
        values = self._plane.values()
        if not values:
            min_x = max_x = min_y = max_y = 0.
        else:
            min_x = min(x for x, _ in values)
            max_x = max(x for x, _ in values)
            min_y = min(y for _, y in values)
            max_y = max(y for _, y in values)
            if clean2d and len(self) > 1 and max_y - min_y < .01 and max_x - min_x < 0.01:
                try:
                    self.clean2d()
                except (ImportError, Exception):
                    pass
                else:
                    values = self._plane.values()
                    min_x = min(x for x, _ in values)
                    max_x = max(x for x, _ in values)
                    min_y = min(y for _, y in values)
                    max_y = max(y for _, y in values)

        config = _render_config
        bonds = self._render_bonds()
        atoms, masks = self._render_atoms()
        if _embedding:
            uid = str(uuid4())
            return atoms, bonds, [], self._masks_svg(masks), uid, min_x, min_y, max_x, max_y

        font_size = config['font_size']
        font125 = 1.25 * font_size
        ep = self._extra_padding
        _width = max_x - min_x + 3.0 * font_size + ep
        _height = max_y - min_y + 2.5 * font_size + ep
        viewbox_x = min_x - font125 - ep / 2
        viewbox_y = -max_y - font125 - ep / 2

        if width is None:
            width = f'{_width:.2f}cm'
        if height is None:
            height = f'{_height:.2f}cm'

        svg = [f'<svg width="{width}" height="{height}" '
               f'viewBox="{viewbox_x:.2f} {viewbox_y:.2f} {_width:.2f} '
               f'{_height:.2f}" xmlns="http://www.w3.org/2000/svg" version="1.1">']
        svg.extend(self._graph_svg(atoms, bonds, masks, viewbox_x, viewbox_y, _width, _height))
        svg.append('</svg>')
        svg = '\n'.join(svg)
        if format == 'svg':
            return svg
        elif format == 'png':
            return svg2png(svg, png_width, png_height, png_scale)
        elif format == 'svgz':
            return compress(svg.encode(), 9)
        raise ValueError(f'format must be svg, png or svgz, not {format}')

    @classmethod
    def _graph_svg(cls, atoms, bonds, masks, viewbox_x, viewbox_y, width, height):
        config = _render_config
        svg = ['  <g>']
        if bonds:
            if masks:
                uid = str(uuid4())
                svg.append(f'    <defs>\n      <mask id="mask-{uid}">\n'
                           f'        <rect x="{viewbox_x:.2f}" y="{viewbox_y:.2f}" '
                           f'width="{width:.2f}" height="{height:.2f}" fill="white"/>')
                svg.extend(cls._masks_svg(masks))
                svg.append('      </mask>\n    </defs>\n'
                           f'    <g fill="none" stroke="{config["bond_color"]}" '
                           f'stroke-width="{config["bond_width"]:.2f}"  mask="url(#mask-{uid})">')
                if len(bonds) == 1:  # SVG BUG adhoc
                    svg.append(f'      <line x1="{viewbox_x:.2f}" y1="{viewbox_y:.2f}" '
                               f'x2="{viewbox_x + width:.2f}" y2="{viewbox_y:.2f}" stroke="none"/>')
            else:
                svg.append(f'    <g fill="none" stroke="{config["bond_color"]}" '
                           f'stroke-width="{config["bond_width"]:.2f}">')
            svg.extend(bonds)
            svg.append('    </g>')

        if atoms:
            svg.append('    <g font-family="monospace">')
            svg.extend(atoms)
            svg.append('    </g>')
        svg.append('  </g>')
        return svg

    @classmethod
    def _masks_svg(cls, masks):
        config = _render_config

        font_size = config['font_size']
        other_size = config['other_size']
        mapping_size = config['mapping_size']
        other_font_style = config['other_font_style']
        symbols_font_style = config['symbols_font_style']

        svg = []
        stroke_width_s = font_size * .1
        stroke_width_o = other_size * .1
        stroke_width_m = mapping_size * .1

        if 'center' in masks:
            svg.append('        <g fill="black">')
            svg.extend(masks['center'])
            svg.append('        </g>')

        svg.append(f'        <g font-family="monospace" stroke="black">')
        if 'symbols' in masks:
            svg.append(f'          <g font-family="{symbols_font_style}" font-size="{font_size:.2f}" '
                       f'stroke-width="{stroke_width_s:.2f}">')
            svg.extend(masks['symbols'])
            svg.append('          </g>')
        if 'aam' in masks:
            svg.append(f'          <g font-size="{mapping_size:.2f}" stroke-width="{stroke_width_m:.2f}">')
            svg.extend(masks['aam'])
            svg.append('          </g>')
        if 'other' in masks:
            svg.append(f'          <g font-family="{other_font_style}" font-size="{other_size}" '
                       f'stroke-width="{stroke_width_o:.2f}">')
            svg.extend(masks['other'])
            svg.append('          </g>')
        if 'span' in masks:
            svg.append(f'          <g font-family="{symbols_font_style}" font-size="{font_size:.2f}" '
                       f'stroke-width="{stroke_width_o:.2f}">')
            svg.extend(masks['span'])
            svg.append('          </g>')
        svg.append('        </g>')
        return svg

    def _render_bonds(self):
        plane = self._plane
        config = _render_config

        broken = config['broken_color']
        formed = config['formed_color']
        dash1, dash2 = config['dashes']
        double_space = config['double_space']
        triple_space = config['triple_space']

        svg = []
        ar_bond_colors = defaultdict(dict)
        for n, m, bond in self.bonds():
            order, p_order = bond.order, bond.p_order
            nx, ny = plane[n]
            mx, my = plane[m]
            ny, my = -ny, -my
            rv = partial(rotate_vector, 0, x2=mx - nx, y2=ny - my)
            if order == 1:
                if p_order == 1:
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"/>')
                elif p_order == 4:
                    ar_bond_colors[n][m] = ar_bond_colors[m][n] = formed
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"/>')
                elif p_order == 2:
                    dx, dy = rv(double_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{formed}"/>')
                elif p_order == 3:
                    dx, dy = rv(triple_space)
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"'
                               f' stroke="{formed}"/>')
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{formed}"/>')
                elif p_order is None:
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"'
                               f' stroke="{broken}"/>')
                else:
                    dx, dy = rv(double_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" x2="{mx + dx:.2f}"'
                               f' y2="{my - dy:.2f}" stroke-dasharray="{dash1:.2f} {dash2:.2f}" stroke="{formed}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{broken}"/>')
            elif order == 4:
                if p_order == 4:
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"/>')
                elif p_order == 1:
                    ar_bond_colors[n][m] = ar_bond_colors[m][n] = broken
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"/>')
                elif p_order == 2:
                    ar_bond_colors[n][m] = ar_bond_colors[m][n] = broken
                    dx, dy = rv(double_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{formed}"/>')
                elif p_order == 3:
                    ar_bond_colors[n][m] = ar_bond_colors[m][n] = broken
                    dx, dy = rv(triple_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}"  stroke="{formed}"/>')
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{formed}"/>')
                elif p_order is None:
                    ar_bond_colors[n][m] = ar_bond_colors[m][n] = broken
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"'
                               f' stroke="{broken}"/>')
                else:
                    ar_bond_colors[n][m] = ar_bond_colors[m][n] = None
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"'
                               f' stroke="{broken}"/>')
            elif order == 2:
                if p_order == 2:
                    dx, dy = rv(double_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}"/>')
                elif p_order == 1:
                    dx, dy = rv(double_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{broken}"/>')
                elif p_order == 4:
                    ar_bond_colors[n][m] = ar_bond_colors[m][n] = formed
                    dx, dy = rv(double_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{broken}"/>')
                elif p_order == 3:
                    dx, dy = rv(triple_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}"/>')
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{formed}"/>')
                elif p_order is None:
                    dx, dy = rv(double_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}" stroke="{broken}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{broken}"/>')
                else:
                    dx, dy = rv(triple_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" x2="{mx + dx:.2f}"'
                               f' y2="{my - dy:.2f}" stroke-dasharray="{dash1:.2f} {dash2:.2f}" stroke="{formed}"/>')
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"'
                               f' stroke="{broken}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{broken}"/>')
            elif order == 3:
                if p_order == 3:
                    dx, dy = rv(triple_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}"/>')
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}"/>')
                elif p_order == 1:
                    dx, dy = rv(triple_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}"/>')
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"'
                               f' stroke="{broken}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" 'f'stroke="{broken}"/>')
                elif p_order == 4:
                    ar_bond_colors[n][m] = ar_bond_colors[m][n] = formed
                    dx, dy = rv(triple_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" x2="{mx + dx:.2f}" '
                               f'y2="{my - dy:.2f}" stroke="{broken}"/>')
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" x2="{mx - dx:.2f}" '
                               f'y2="{my + dy:.2f}" stroke="{broken}"/>')
                elif p_order == 2:
                    dx, dy = rv(triple_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}"/>')
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{broken}"/>')
                elif p_order is None:
                    dx, dy = rv(triple_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}" stroke="{broken}"/>')
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" '
                               f'x2="{mx:.2f}" y2="{my:.2f}" stroke="{broken}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{broken}"/>')
                else:
                    dx, dy = rv(double_space)
                    dx3 = 3 * dx
                    dy3 = 3 * dy
                    svg.append(f'      <line x1="{nx + dx3:.2f}" y1="{ny - dy3:.2f}" x2="{mx + dx3:.2f}" '
                               f'y2="{my - dy3:.2f}" stroke-dasharray="{dash1:.2f} {dash2:.2f}" stroke="{formed}"/>')
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}" stroke="{broken}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{broken}"/>')
                    svg.append(f'      <line x1="{nx - dx3:.2f}" y1="{ny + dy3:.2f}" x2="{mx - dx3:.2f}" '
                               f'y2="{my + dy3:.2f}" stroke="{broken}"/>')
            elif order is None:
                if p_order == 1:
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"'
                               f' stroke="{formed}"/>')
                elif p_order == 4:
                    ar_bond_colors[n][m] = ar_bond_colors[m][n] = formed
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"'
                               f' stroke="{formed}"/>')
                elif p_order == 2:
                    dx, dy = rv(double_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" x2="{mx + dx:.2f}" '
                               f'y2="{my - dy:.2f}" stroke="{formed}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" x2="{mx - dx:.2f}" '
                               f'y2="{my + dy:.2f}" stroke="{formed}"/>')
                elif p_order == 3:
                    dx, dy = rv(triple_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}" stroke="{formed}"/>')
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"'
                               f' stroke="{formed}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{formed}"/>')
                else:
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}" '
                               f'stroke-dasharray="{dash1:.2f} {dash2:.2f}" stroke="{formed}"/>')
            else:
                if p_order == 8:
                    svg.append(f'        <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}" '
                               f'stroke-dasharray="{dash1:.2f} {dash2:.2f}"/>')
                elif p_order == 1:
                    dx, dy = rv(double_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" x2="{mx + dx:.2f}"'
                               f' y2="{my - dy:.2f}" stroke-dasharray="{dash1:.2f} {dash2:.2f}" stroke="{broken}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{formed}"/>')
                elif p_order == 4:
                    ar_bond_colors[n][m] = ar_bond_colors[m][n] = None
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"'
                               f' stroke="{formed}"/>')
                elif p_order == 2:
                    dx, dy = rv(triple_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" x2="{mx + dx:.2f}"'
                               f' y2="{my - dy:.2f}" stroke-dasharray="{dash1:.2f} {dash2:.2f}" stroke="{broken}"/>')
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"'
                               f' stroke="{formed}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{formed}"/>')
                elif p_order == 3:
                    dx, dy = rv(double_space)
                    dx3 = 3 * dx
                    dy3 = 3 * dy
                    svg.append(f'      <line x1="{nx + dx3:.2f}" y1="{ny - dy3:.2f}" x2="{mx + dx3:.2f}" '
                               f'y2="{my - dy3:.2f}" stroke-dasharray="{dash1:.2f} {dash2:.2f}" stroke="{broken}"/>')
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}" stroke="{formed}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{formed}"/>')
                    svg.append(f'      <line x1="{nx - dx3:.2f}" y1="{ny + dy3:.2f}" '
                               f'x2="{mx - dx3:.2f}" y2="{my + dy3:.2f}" stroke="{formed}"/>')
                else:
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}" '
                               f'stroke-dasharray="{dash1:.2f} {dash2:.2f}" stroke="{broken}"/>')

        for ring in self.aromatic_rings:
            cx = sum(plane[x][0] for x in ring) / len(ring)
            cy = sum(plane[x][1] for x in ring) / len(ring)

            for n, m in zip(ring, ring[1:]):
                nx, ny = plane[n]
                mx, my = plane[m]
                aromatic = self.__render_aromatic_bond(nx, ny, mx, my, cx, cy, ar_bond_colors[n].get(m))
                if aromatic:
                    svg.append(aromatic)

            n, m = ring[-1], ring[0]
            nx, ny = plane[n]
            mx, my = plane[m]
            aromatic = self.__render_aromatic_bond(nx, ny, mx, my, cx, cy, ar_bond_colors[n].get(m))
            if aromatic:
                svg.append(aromatic)
        return svg

    def __render_aromatic_bond(self, n_x, n_y, m_x, m_y, c_x, c_y, color):
        config = _render_config

        dash1, dash2 = config['dashes']
        dash3, dash4 = config['aromatic_dashes']
        aromatic_space = config['aromatic_space']
        # n aligned xy
        mn_x, mn_y, cn_x, cn_y = m_x - n_x, m_y - n_y, c_x - n_x, c_y - n_y

        # nm reoriented xy
        mr_x, mr_y = hypot(mn_x, mn_y), 0
        cr_x, cr_y = rotate_vector(cn_x, cn_y, mn_x, -mn_y)

        if cr_y and aromatic_space / cr_y < .65:
            if cr_y > 0:
                r_y = aromatic_space
            else:
                r_y = -aromatic_space
                cr_y = -cr_y

            ar_x = aromatic_space * cr_x / cr_y
            br_x = mr_x - aromatic_space * (mr_x - cr_x) / cr_y

            # backward reorienting
            an_x, an_y = rotate_vector(ar_x, r_y, mn_x, mn_y)
            bn_x, bn_y = rotate_vector(br_x, r_y, mn_x, mn_y)
            if color:
                return f'      <line x1="{an_x + n_x:.2f}" y1="{-an_y - n_y:.2f}" x2="{bn_x + n_x:.2f}" ' \
                       f'y2="{-bn_y - n_y:.2f}" stroke-dasharray="{dash3:.2f} {dash4:.2f}" stroke="{color}"/>'
            elif color is None:
                dash3, dash4 = dash1, dash2
            return f'      <line x1="{an_x + n_x:.2f}" y1="{-an_y - n_y:.2f}"' \
                   f' x2="{bn_x + n_x:.2f}" y2="{-bn_y - n_y:.2f}" stroke-dasharray="{dash3:.2f} {dash4:.2f}"/>'

    def _render_atoms(self):
        bonds = self._bonds
        plane = self._plane
        charges = self._charges
        radicals = self._radicals
        p_charges = self._p_charges
        config = _render_config
        p_radicals = self._p_radicals

        carbon = config['carbon']
        font_size = config['font_size']
        other_size = config['other_size']
        monochrome = config['monochrome']
        atoms_colors = config['atoms_colors']
        dx_m, dy_m = config['dx_m'], config['dy_m']
        dx_ci, dy_ci = config['dx_ci'], config['dy_ci']
        symbols_font_style = config['symbols_font_style']
        font2 = .2 * font_size
        font4 = .4 * font_size
        font6 = .6 * font_size
        font7 = .7 * font_size
        font15 = .15 * font_size
        font25 = .25 * font_size
        din_charges = {m[0]: m[1] != n[1] for m, n in zip(charges.items(), p_charges.items())}
        din_radicals = {m[0]: m[1] != n[1] for m, n in zip(radicals.items(), p_radicals.items())}

        if monochrome:
            other_fill = 'black'
        else:
            other_fill = config['other_color']

        svg = []
        others = []
        mask = defaultdict(list)
        for n, atom in self._atoms.items():
            x, y = plane[n]
            y = -y
            symbol = atom.atomic_symbol
            if not bonds[n] or symbol != 'C' or carbon or atom.charge or atom.is_radical or atom.isotope \
                    or din_charges[n] or din_radicals[n]:

                if radicals[n]:
                    r = '↑' if p_radicals[n] else '↑↓'
                elif p_radicals[n]:
                    r = '↓↑'
                else:
                    r = ''

                if charges[n] != p_charges[n]:
                    t = _render_p_charge[charges[n]][p_charges[n]]
                    others.append(f'        <text x="{x:.2f}" y="{y:.2f}" dx="{dx_ci:.2f}" dy="-{dy_ci:.2f}" '
                                  f'font-size="{other_size:.2f}">{t}{r}</text>')
                    mask['other'].append(f'           <text x="{x:.2f}" y="{y:.2f}" dx="{dx_ci:.2f}" dy="-{dy_ci:.2f}">'
                                         f'{t}{r}</text>')
                if charges[n]:
                    t = _render_charge[charges[n]]
                    others.append(f'        <text x="{x:.2f}" y="{y:.2f}" dx="{dx_ci:.2f}" dy="-{dy_ci:.2f}" '
                                  f'font-size="{other_size:.2f}">{t}{r}</text>')
                    mask['other'].append(f'           <text x="{x:.2f}" y="{y:.2f}" dx="{dx_ci:.2f}" dy="-{dy_ci:.2f}">'
                                         f'{t}{r}</text>')
                elif r:
                    others.append(f'        <text x="{x:.2f}" y="{y:.2f}" dx="{dx_ci:.2f}" dy="-{dy_ci:.2f}" '
                                  f'font-size="{other_size:.2f}">{r}</text>')
                    mask['other'].append(f'           <text x="{x:.2f}" y="{y:.2f}" dx="{dx_ci:.2f}" dy="-{dy_ci:.2f}">'
                                         f'{r}</text>')

                if atom.isotope:
                    t = atom.isotope
                    others.append(f'        <text x="{x:.2f}" y="{y:.2f}" dx="-{dx_ci:.2f}" dy="-{dy_ci:.2f}" '
                                  f'font-size="{other_size:.2f}" text-anchor="end">{t}</text>')
                    mask['other'].append(f'           <text x="{x:.2f}" y="{y:.2f}" dx="-{dx_ci:.2f}"'
                                         f' dy="-{dy_ci:.2f}">{t}</text>')

                svg.append(f'      <g fill="{"black" if monochrome else atoms_colors[atom.atomic_number - 1]}" '
                           f'font-family="{symbols_font_style}">')
                if len(symbol) > 1:
                    dx = font7
                    if symbol[-1] in ('l', 'i', 'r', 't'):
                        rx = font6
                        ax = font25
                    else:
                        rx = font7
                        ax = font15
                    mask['center'].append(f'          <ellipse cx="{x - ax:.2f}" cy="{y:.2f}" rx="{rx}" ry="{font4}"/>')
                else:
                    if symbol == 'I':
                        dx = font15
                        dx_mm = dx_m
                    else:
                        dx = font4
                        dx_mm = dx_m + font2
                    mask['center'].append(f'          <circle cx="{x:.2f}" cy="{y:.2f}" r="{font4:.2f}"/>')
                svg.append(f'        <text x="{x:.2f}" y="{y:.2f}" dx="-{dx:.2f}" dy="{font4:.2f}" '
                           f'font-size="{font_size:.2f}">{symbol}</text>')
                mask['symbols'].append(f'           <text x="{x:.2f}" y="{y:.2f}" dx="-{dx:.2f}" '
                                       f'dy="{font4:.2f}">{symbol}</text>')
                svg.append('      </g>')

        if others:
            svg.append(f'      <g font-family="{config["other_font_style"]}" fill="{other_fill}" '
                       f'font-size="{other_size:.2f}">')
            svg.extend(others)
            svg.append('      </g>')
        return svg, mask


class DepictCGR(Depict):
    __slots__ = ()

    def _render_bonds(self):
        plane = self._plane
        config = _render_config

        broken = config['broken_color']
        formed = config['formed_color']
        dash1, dash2 = config['dashes']
        double_space = config['double_space']
        triple_space = config['triple_space']
        dw = config['dynamic_bond_width']
        sw = f' stroke-width="{dw:.2f}"'

        svg = []
        ar_bond_colors = defaultdict(dict)
        for n, m, bond in self.bonds():
            order, p_order = bond.order, bond.p_order
            nx, ny = plane[n]
            mx, my = plane[m]
            ny, my = -ny, -my
            rv = partial(rotate_vector, 0, x2=mx - nx, y2=ny - my)
            if order == 1:
                if p_order == 1:
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"/>')
                elif p_order == 4:
                    ar_bond_colors[n][m] = ar_bond_colors[m][n] = formed
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"/>')
                elif p_order == 2:
                    dx, dy = rv(double_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{formed}"{sw}/>')
                elif p_order == 3:
                    dx, dy = rv(triple_space)
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"'
                               f' stroke="{formed}"{sw}/>')
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{formed}"{sw}/>')
                elif p_order is None:
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"'
                               f' stroke="{broken}"{sw}/>')
                else:
                    dx, dy = rv(double_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" x2="{mx + dx:.2f}"'
                               f' y2="{my - dy:.2f}" stroke-dasharray="{dash1:.2f} {dash2:.2f}" stroke="{formed}"{sw}/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{broken}"{sw}/>')
            elif order == 4:
                if p_order == 4:
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"/>')
                elif p_order == 1:
                    ar_bond_colors[n][m] = ar_bond_colors[m][n] = broken
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"/>')
                elif p_order == 2:
                    ar_bond_colors[n][m] = ar_bond_colors[m][n] = broken
                    dx, dy = rv(double_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{formed}"{sw}/>')
                elif p_order == 3:
                    ar_bond_colors[n][m] = ar_bond_colors[m][n] = broken
                    dx, dy = rv(triple_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}"  stroke="{formed}"{sw}/>')
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{formed}"{sw}/>')
                elif p_order is None:
                    ar_bond_colors[n][m] = ar_bond_colors[m][n] = broken
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"'
                               f' stroke="{broken}"{sw}/>')
                else:
                    ar_bond_colors[n][m] = ar_bond_colors[m][n] = None
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"'
                               f' stroke="{broken}"{sw}/>')
            elif order == 2:
                if p_order == 2:
                    dx, dy = rv(double_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}"/>')
                elif p_order == 1:
                    dx, dy = rv(double_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{broken}"{sw}/>')
                elif p_order == 4:
                    ar_bond_colors[n][m] = ar_bond_colors[m][n] = formed
                    dx, dy = rv(double_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{broken}"{sw}/>')
                elif p_order == 3:
                    dx, dy = rv(triple_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}"/>')
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{formed}"{sw}/>')
                elif p_order is None:
                    dx, dy = rv(double_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}" stroke="{broken}"{sw}/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{broken}"{sw}/>')
                else:
                    dx, dy = rv(triple_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" x2="{mx + dx:.2f}"'
                               f' y2="{my - dy:.2f}" stroke-dasharray="{dash1:.2f} {dash2:.2f}" stroke="{formed}"{sw}/>')
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"'
                               f' stroke="{broken}"{sw}/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{broken}"{sw}/>')
            elif order == 3:
                if p_order == 3:
                    dx, dy = rv(triple_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}"/>')
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}"/>')
                elif p_order == 1:
                    dx, dy = rv(triple_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}"/>')
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"'
                               f' stroke="{broken}"{sw}/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" 'f'stroke="{broken}"{sw}/>')
                elif p_order == 4:
                    ar_bond_colors[n][m] = ar_bond_colors[m][n] = formed
                    dx, dy = rv(triple_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" x2="{mx + dx:.2f}" '
                               f'y2="{my - dy:.2f}" stroke="{broken}"{sw}/>')
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" x2="{mx - dx:.2f}" '
                               f'y2="{my + dy:.2f}" stroke="{broken}"{sw}/>')
                elif p_order == 2:
                    dx, dy = rv(triple_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}"/>')
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{broken}"{sw}/>')
                elif p_order is None:
                    dx, dy = rv(triple_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}" stroke="{broken}"{sw}/>')
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" '
                               f'x2="{mx:.2f}" y2="{my:.2f}" stroke="{broken}"{sw}/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{broken}"{sw}/>')
                else:
                    dx, dy = rv(double_space)
                    dx3 = 3 * dx
                    dy3 = 3 * dy
                    svg.append(f'      <line x1="{nx + dx3:.2f}" y1="{ny - dy3:.2f}" x2="{mx + dx3:.2f}" '
                               f'y2="{my - dy3:.2f}" stroke-dasharray="{dash1:.2f} {dash2:.2f}" stroke="{formed}"{sw}/>')
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}" stroke="{broken}"{sw}/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{broken}"{sw}/>')
                    svg.append(f'      <line x1="{nx - dx3:.2f}" y1="{ny + dy3:.2f}" x2="{mx - dx3:.2f}" '
                               f'y2="{my + dy3:.2f}" stroke="{broken}"{sw}/>')
            elif order is None:
                if p_order == 1:
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"'
                               f' stroke="{formed}"{sw}/>')
                elif p_order == 4:
                    ar_bond_colors[n][m] = ar_bond_colors[m][n] = formed
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"'
                               f' stroke="{formed}"{sw}/>')
                elif p_order == 2:
                    dx, dy = rv(double_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" x2="{mx + dx:.2f}" '
                               f'y2="{my - dy:.2f}" stroke="{formed}"{sw}/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" x2="{mx - dx:.2f}" '
                               f'y2="{my + dy:.2f}" stroke="{formed}"{sw}/>')
                elif p_order == 3:
                    dx, dy = rv(triple_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}" stroke="{formed}"{sw}/>')
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"'
                               f' stroke="{formed}"{sw}/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{formed}"{sw}/>')
                else:
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}" '
                               f'stroke-dasharray="{dash1:.2f} {dash2:.2f}" stroke="{formed}"{sw}/>')
            else:
                if p_order == 8:
                    svg.append(f'        <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}" '
                               f'stroke-dasharray="{dash1:.2f} {dash2:.2f}"/>')
                elif p_order == 1:
                    dx, dy = rv(double_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" x2="{mx + dx:.2f}"'
                               f' y2="{my - dy:.2f}" stroke-dasharray="{dash1:.2f} {dash2:.2f}" stroke="{broken}"{sw}/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{formed}"{sw}/>')
                elif p_order == 4:
                    ar_bond_colors[n][m] = ar_bond_colors[m][n] = None
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"'
                               f' stroke="{formed}"{sw}/>')
                elif p_order == 2:
                    dx, dy = rv(triple_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" x2="{mx + dx:.2f}"'
                               f' y2="{my - dy:.2f}" stroke-dasharray="{dash1:.2f} {dash2:.2f}" stroke="{broken}"{sw}/>')
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"'
                               f' stroke="{formed}"{sw}/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{formed}"{sw}/>')
                elif p_order == 3:
                    dx, dy = rv(double_space)
                    dx3 = 3 * dx
                    dy3 = 3 * dy
                    svg.append(f'      <line x1="{nx + dx3:.2f}" y1="{ny - dy3:.2f}" x2="{mx + dx3:.2f}" '
                               f'y2="{my - dy3:.2f}" stroke-dasharray="{dash1:.2f} {dash2:.2f}" stroke="{broken}"{sw}/>')
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}" stroke="{formed}"{sw}/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{formed}"{sw}/>')
                    svg.append(f'      <line x1="{nx - dx3:.2f}" y1="{ny + dy3:.2f}" '
                               f'x2="{mx - dx3:.2f}" y2="{my + dy3:.2f}" stroke="{formed}"{sw}/>')
                else:
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}" '
                               f'stroke-dasharray="{dash1:.2f} {dash2:.2f}" stroke="{broken}"{sw}/>')

        for ring in self.aromatic_rings:
            cx = sum(plane[x][0] for x in ring) / len(ring)
            cy = sum(plane[x][1] for x in ring) / len(ring)

            for n, m in zip(ring, ring[1:]):
                nx, ny = plane[n]
                mx, my = plane[m]
                aromatic = self.__render_aromatic_bond(nx, ny, mx, my, cx, cy, ar_bond_colors[n].get(m))
                if aromatic:
                    svg.append(aromatic)

            n, m = ring[-1], ring[0]
            nx, ny = plane[n]
            mx, my = plane[m]
            aromatic = self.__render_aromatic_bond(nx, ny, mx, my, cx, cy, ar_bond_colors[n].get(m))
            if aromatic:
                svg.append(aromatic)
        return svg

    def __render_aromatic_bond(self, n_x, n_y, m_x, m_y, c_x, c_y, color):
        config = _render_config

        dash1, dash2 = config['dashes']
        dash3, dash4 = config['aromatic_dashes']
        aromatic_space = config['cgr_aromatic_space']
        # n aligned xy
        mn_x, mn_y, cn_x, cn_y = m_x - n_x, m_y - n_y, c_x - n_x, c_y - n_y

        # nm reoriented xy
        mr_x, mr_y = hypot(mn_x, mn_y), 0
        cr_x, cr_y = rotate_vector(cn_x, cn_y, mn_x, -mn_y)

        if cr_y and aromatic_space / cr_y < .65:
            if cr_y > 0:
                r_y = aromatic_space
            else:
                r_y = -aromatic_space
                cr_y = -cr_y

            ar_x = aromatic_space * cr_x / cr_y
            br_x = mr_x - aromatic_space * (mr_x - cr_x) / cr_y

            # backward reorienting
            an_x, an_y = rotate_vector(ar_x, r_y, mn_x, mn_y)
            bn_x, bn_y = rotate_vector(br_x, r_y, mn_x, mn_y)
            if color:
                return f'      <line x1="{an_x + n_x:.2f}" y1="{-an_y - n_y:.2f}" x2="{bn_x + n_x:.2f}" ' \
                       f'y2="{-bn_y - n_y:.2f}" stroke-dasharray="{dash3:.2f} {dash4:.2f}" stroke="{color}"/>'
            elif color is None:
                dash3, dash4 = dash1, dash2
            return f'      <line x1="{an_x + n_x:.2f}" y1="{-an_y - n_y:.2f}"' \
                   f' x2="{bn_x + n_x:.2f}" y2="{-bn_y - n_y:.2f}" stroke-dasharray="{dash3:.2f} {dash4:.2f}"/>'

    def _render_atoms(self):
        bonds = self._bonds
        plane = self._plane
        charges = self._charges
        radicals = self._radicals
        p_charges = self._p_charges
        config = _render_config
        p_radicals = self._p_radicals

        carbon = config['carbon']
        font_size = config['font_size']
        other_size = config['other_size']
        monochrome = config['monochrome']
        atoms_colors = config['atoms_colors']
        dx_ci, dy_ci = config['dx_ci'], config['dy_ci']
        symbols_font_style = config['symbols_font_style']
        font4 = .4 * font_size
        font6 = .6 * font_size
        font7 = .7 * font_size
        font15 = .15 * font_size
        font25 = .25 * font_size
        din_charges = {m[0]: m[1] != n[1] for m, n in zip(charges.items(), p_charges.items())}
        din_radicals = {m[0]: m[1] != n[1] for m, n in zip(radicals.items(), p_radicals.items())}

        if monochrome:
            other_fill = 'black'
        else:
            other_fill = config['other_color']

        svg = []
        others = []
        mask = defaultdict(list)
        for n, atom in self._atoms.items():
            x, y = plane[n]
            y = -y
            symbol = atom.atomic_symbol
            if not bonds[n] or symbol != 'C' or carbon or atom.charge or atom.is_radical or atom.isotope \
                    or din_charges[n] or din_radicals[n]:

                if radicals[n]:
                    r = '↑' if p_radicals[n] else '↑↓'
                elif p_radicals[n]:
                    r = '↓↑'
                else:
                    r = ''

                if charges[n] != p_charges[n]:
                    t = _render_p_charge[charges[n]][p_charges[n]]
                    others.append(f'        <text x="{x:.2f}" y="{y:.2f}" dx="{dx_ci:.2f}" dy="-{dy_ci:.2f}" '
                                  f'font-size="{other_size:.2f}">{t}{r}</text>')
                    mask['other'].append(f'           <text x="{x:.2f}" y="{y:.2f}" dx="{dx_ci:.2f}" dy="-{dy_ci:.2f}">'
                                         f'{t}{r}</text>')
                if charges[n]:
                    t = _render_charge[charges[n]]
                    others.append(f'        <text x="{x:.2f}" y="{y:.2f}" dx="{dx_ci:.2f}" dy="-{dy_ci:.2f}" '
                                  f'font-size="{other_size:.2f}">{t}{r}</text>')
                    mask['other'].append(f'           <text x="{x:.2f}" y="{y:.2f}" dx="{dx_ci:.2f}" dy="-{dy_ci:.2f}">'
                                         f'{t}{r}</text>')
                elif r:
                    others.append(f'        <text x="{x:.2f}" y="{y:.2f}" dx="{dx_ci:.2f}" dy="-{dy_ci:.2f}" '
                                  f'font-size="{other_size:.2f}">{r}</text>')
                    mask['other'].append(f'           <text x="{x:.2f}" y="{y:.2f}" dx="{dx_ci:.2f}" dy="-{dy_ci:.2f}">'
                                         f'{r}</text>')

                if atom.isotope:
                    t = atom.isotope
                    others.append(f'        <text x="{x:.2f}" y="{y:.2f}" dx="-{dx_ci:.2f}" dy="-{dy_ci:.2f}" '
                                  f'font-size="{other_size:.2f}" text-anchor="end">{t}</text>')
                    mask['other'].append(f'           <text x="{x:.2f}" y="{y:.2f}" dx="-{dx_ci:.2f}"'
                                         f' dy="-{dy_ci:.2f}">{t}</text>')

                svg.append(f'      <g fill="{"black" if monochrome else atoms_colors[atom.atomic_number - 1]}" '
                           f'font-family="{symbols_font_style}">')
                if len(symbol) > 1:
                    dx = font7
                    if symbol[-1] in ('l', 'i', 'r', 't'):
                        rx = font6
                        ax = font25
                    else:
                        rx = font7
                        ax = font15
                    mask['center'].append(f'          <ellipse cx="{x - ax:.2f}" cy="{y:.2f}" rx="{rx}" ry="{font4}"/>')
                else:
                    if symbol == 'I':
                        dx = font15
                    else:
                        dx = font4
                    mask['center'].append(f'          <circle cx="{x:.2f}" cy="{y:.2f}" r="{font4:.2f}"/>')
                svg.append(f'        <text x="{x:.2f}" y="{y:.2f}" dx="-{dx:.2f}" dy="{font4:.2f}" '
                           f'font-size="{font_size:.2f}">{symbol}</text>')
                mask['symbols'].append(f'           <text x="{x:.2f}" y="{y:.2f}" dx="-{dx:.2f}" '
                                       f'dy="{font4:.2f}">{symbol}</text>')
                svg.append('      </g>')

        if others:
            svg.append(f'      <g font-family="{config["other_font_style"]}" fill="{other_fill}" '
                       f'font-size="{other_size:.2f}">')
            svg.extend(others)
            svg.append('      </g>')
        return svg, mask


class DepictQuery(Depict):
    __slots__ = ()
    _extra_padding = .5

    def _render_bonds(self):
        svg = []
        plane = self._plane
        config = _render_config

        dash1, dash2 = config['dashes']
        double_space = config['double_space']
        triple_space = config['triple_space']
        dash3, dash4 = config['aromatic_dashes']
        for n, m, bond in self.bonds():
            nx, ny = plane[n]
            mx, my = plane[m]
            ny, my = -ny, -my
            if bond == 1:
                svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"/>')
            elif bond == 4:
                dx, dy = rotate_vector(0, double_space, mx - nx, ny - my)
                svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" x2="{mx + dx:.2f}" y2="{my - dy:.2f}"/>')
                svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" x2="{mx - dx:.2f}" y2="{my + dy:.2f}" '
                           f'stroke-dasharray="{dash3:.2f} {dash4:.2f}"/>')
            elif bond == 2:
                dx, dy = rotate_vector(0, double_space, mx - nx, ny - my)
                svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" x2="{mx + dx:.2f}" y2="{my - dy:.2f}"/>')
                svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" x2="{mx - dx:.2f}" y2="{my + dy:.2f}"/>')
            elif bond == 3:
                dx, dy = rotate_vector(0, triple_space, mx - nx, ny - my)
                svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" x2="{mx + dx:.2f}" y2="{my - dy:.2f}"/>')
                svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"/>')
                svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" x2="{mx - dx:.2f}" y2="{my + dy:.2f}"/>')
            else:  # other query bonds
                svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}" '
                           f'stroke-dasharray="{dash1:.2f} {dash2:.2f}"/>')
        return svg

    def _render_atoms(self):
        plane = self._plane
        config = _render_config

        carbon = config['carbon']
        mapping = config['mapping']
        font_size = config['font_size']
        other_size = config['other_size']
        monochrome = config['monochrome']
        atoms_colors = config['atoms_colors']
        mapping_size = config['mapping_size']
        dx_m, dy_m = config['dx_m'], config['dy_m']
        other_font_style = config['other_font_style']
        dx_ci, dy_ci = config['dx_ci'], config['dy_ci']
        dx_nh, dy_nh = config['dx_nh'], config['dy_nh']
        symbols_font_style = config['symbols_font_style']

        font2 = .2 * font_size
        font3 = .3 * font_size
        font4 = .4 * font_size
        font5 = .5 * font_size
        font6 = .6 * font_size
        font7 = .7 * font_size
        font15 = .15 * font_size
        font25 = .25 * font_size

        if monochrome:
            map_fill = query_fill = other_fill = 'black'
        else:
            other_fill = config['other_color']
            map_fill = config['mapping_color']
            query_fill = config['query_color']

        svg = []
        maps = []
        others = []
        nghbrs = []
        hbrdztns = []
        mask = defaultdict(list)
        for n, atom in self._atoms.items():
            x, y = plane[n]
            y = -y
            neighbor_svgs = [(plane[m][0], -plane[m][1]) for m in self._bonds[n]]
            single = not self._bonds[n]
            symbol = atom.atomic_symbol
            if single or symbol != 'C' or carbon or atom.charge or atom.is_radical:
                if atom.charge:
                    t = _render_charge[atom.charge]
                    tt = "↑" if atom.is_radical else ""
                    others.append(f'        <text x="{x:.2f}" y="{y:.2f}" dx="{dx_ci:.2f}" dy="-{dy_ci:.2f}" '
                                  f'font-size="{other_size:.2f}">{t}{tt}</text>')
                    mask['other'].append(f'           <text x="{x:.2f}" y="{y:.2f}" dx="{dx_ci:.2f}" dy="-{dy_ci:.2f}">'
                                         f'{t}{tt}</text>')
                elif atom.is_radical:
                    others.append(f'        <text x="{x:.2f}" y="{y:.2f}" dx="{dx_ci:.2f}" dy="-{dy_ci:.2f}" '
                                  f'font-size="{other_size:.2f}">↑</text>')
                    mask['other'].append(f'           <text x="{x:.2f}" y="{y:.2f}" dx="{dx_ci:.2f}" dy="-{dy_ci:.2f}">'
                                         f'↑</text>')

                svg.append(f'      <g fill="{"black" if monochrome else atoms_colors[atom.atomic_number - 1]}" '
                           f'font-family="{symbols_font_style}">')
                if len(symbol) > 1:
                    dx = font7
                    dx_mm = dx_m + font5
                    dx_nhh = dx_nh + font5
                    if symbol[-1] in ('l', 'i', 'r', 't'):
                        rx = font6
                        ax = font25
                    else:
                        rx = font7
                        ax = font15
                    mask['center'].append(f'          <ellipse cx="{x - ax:.2f}" cy="{y:.2f}" rx="{rx}" ry="{font4}"/>')
                else:
                    if symbol == 'I':
                        dx = font15
                        dx_mm = dx_m
                        dx_nhh = dx_nh
                    else:
                        dx = font4
                        dx_mm = dx_m + font2
                        dx_nhh = dx_nh + font2
                    mask['center'].append(f'          <circle cx="{x:.2f}" cy="{y:.2f}" r="{font4:.2f}"/>')
                svg.append(f'        <text x="{x:.2f}" y="{y:.2f}" dx="-{dx:.2f}" dy="{font4:.2f}" '
                           f'font-size="{font_size:.2f}">{symbol}</text>')
                mask['symbols'].append(f'           <text x="{x:.2f}" y="{y:.2f}" dx="-{dx:.2f}" '
                                       f'dy="{font4:.2f}">{symbol}</text>')
                svg.append('      </g>')
                dx_nnh = dx_nhh

                if mapping:
                    maps.append(f'        <text x="{x:.2f}" y="{y:.2f}" dx="-{dx_mm:.2f}" dy="{dy_m + font3:.2f}" '
                                f'text-anchor="end">{n}</text>')
                    mask['aam'].append(f'           <text x="{x:.2f}" y="{y:.2f}" dx="-{dx_mm:.2f}" '
                                       f'dy="{dy_m + font3:.2f}" text-anchor="end">{n}</text>')
            elif mapping:
                dx_nnh = dx_nh
                maps.append(f'        <text x="{x:.2f}" y="{y:.2f}" dx="-{dx_m}" dy="{dy_m:.2f}" '
                            f'text-anchor="end">{n}</text>')
                mask['aam'].append(f'           <text x="{x:.2f}" y="{y:.2f}" dx="-{dx_m}" dy="{dy_m:.2f}" '
                                   f'text-anchor="end">{n}</text>')
            else:
                dx_nnh = dx_nh

            if atom.neighbors or atom.hybridization or getattr(atom, '_recursive_smarts', None):
                angle = _optimal_label_direction(x, y, neighbor_svgs)
                offset = hypot(dx_nnh, dy_nh)
                nh_dx = cos(angle) * offset
                nh_dy = sin(angle) * offset
                anchor = 'end' if nh_dx < 0 else 'start'
                level_sign = 1 if nh_dy >= 0 else -1

            if atom.neighbors:
                level = .6 * other_size
                nn = "".join(str(x) for x in atom.neighbors)
                nghbrs.append(f'        <text x="{x:.2f}" y="{y:.2f}" dx="{nh_dx:.2f}" dy="{nh_dy:.2f}" '
                              f'text-anchor="{anchor}">{nn}</text>')
            else:
                level = 0

            if atom.hybridization:
                hh = "".join(_render_hybridization[x] for x in atom.hybridization)
                h_dy = nh_dy + level * level_sign
                hbrdztns.append(f'        <text x="{x:.2f}" y="{y:.2f}" dx="{nh_dx:.2f}" dy="{h_dy:.2f}" '
                                f'text-anchor="{anchor}">{hh}</text>')
                level += .6 * other_size

            if getattr(atom, '_recursive_smarts', None):
                rc = f'$({len(atom._recursive_smarts)})'
                rc_dy = nh_dy + (level + .4 * other_size) * level_sign
                hbrdztns.append(f'        <text x="{x:.2f}" y="{y:.2f}" dx="{nh_dx:.2f}" dy="{rc_dy:.2f}" '
                                f'text-anchor="{anchor}">{rc}</text>')

        if nghbrs:
            svg.append(f'      <g fill="{query_fill}" font-family="{other_font_style}" font-size="{other_size:.2f}">')
            svg.extend(nghbrs)
            if hbrdztns:
                svg.extend(hbrdztns)
            svg.append('      </g>')
        elif hbrdztns:
            svg.append(f'      <g fill="{query_fill}" font-family="{other_font_style}" font-size="{other_size:.2f}">')
            svg.extend(hbrdztns)
            svg.append('    </g>')

        if others:
            svg.append(f'      <g fill="{other_fill}" font-family="{other_font_style}" font-size="{other_size:.2f}">')
            svg.extend(others)
            svg.append('      </g>')

        if mapping:
            svg.append(f'      <g fill="{map_fill}" font-size="{mapping_size:.2f}">')
            svg.extend(maps)
            svg.append('      </g>')

        return svg, mask


class DepictQueryCGR(Depict):
    __slots__ = ()
    _extra_padding = .5

    def _render_bonds(self):
        svg = []
        plane = self._plane
        config = _render_config

        dash1, dash2 = config['dashes']
        broken = config['broken_color']
        formed = config['formed_color']
        double_space = config['double_space']
        triple_space = config['triple_space']
        dash3, dash4 = config['aromatic_dashes']
        dw = config['dynamic_bond_width']
        sw = f' stroke-width="{dw:.2f}"'

        for n, m, bond in self.bonds():
            order, p_order = bond.order, bond.p_order
            nx, ny = plane[n]
            mx, my = plane[m]
            ny, my = -ny, -my
            rv = partial(rotate_vector, 0, x2=mx - nx, y2=ny - my)
            if order == 1:
                if p_order == 1:
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"/>')
                elif p_order == 4:
                    dx, dy = rv(double_space)
                    svg.append(
                        f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" x2="{mx + dx:.2f}" y2="{my - dy:.2f}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" x2="{mx - dx:.2f}"'
                               f' y2="{my + dy:.2f}" stroke-dasharray="{dash3:.2f} {dash4:.2f}" stroke="{formed}"{sw}/>')
                elif p_order == 2:
                    dx, dy = rv(double_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{formed}"{sw}/>')
                elif p_order == 3:
                    dx, dy = rv(triple_space)
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"'
                               f' stroke="{formed}"{sw}/>')
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{formed}"{sw}/>')
                elif p_order is None:
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"'
                               f' stroke="{broken}"{sw}/>')
                else:
                    dx, dy = rv(double_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}" stroke-dasharray="{dash1:.2f} '
                               f'{dash2:.2f}" stroke="{formed}"{sw}/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{broken}"{sw}/>')
            elif order == 4:
                if p_order == 4:
                    dx, dy = rv(double_space)
                    svg.append(
                        f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" x2="{mx + dx:.2f}" y2="{my - dy:.2f}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" x2="{mx - dx:.2f}"'
                               f' y2="{my + dy:.2f}" stroke-dasharray="{dash3:.2f} {dash4:.2f}"/>')
                elif p_order == 1:
                    dx, dy = rv(double_space)
                    svg.append(
                        f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" x2="{mx + dx:.2f}" y2="{my - dy:.2f}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" x2="{mx - dx:.2f}"'
                               f' y2="{my + dy:.2f}" stroke-dasharray="{dash3:.2f} {dash4:.2f}" stroke="{broken}"{sw}/>')
                elif p_order == 2:
                    dx, dy = rv(triple_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}"  stroke="{formed}"{sw}/>')
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" x2="{mx - dx:.2f}"'
                               f' y2="{my + dy:.2f}" stroke-dasharray="{dash3:.2f} {dash4:.2f}" stroke="{broken}"{sw}/>')
                elif p_order == 3:
                    dx, dy = rv(double_space)
                    dx3 = 3 * dx
                    dy3 = 3 * dy
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}" stroke="{formed}"{sw}/>')
                    svg.append(f'      <line x1="{nx - dx3:.2f}" y1="{ny + dy3:.2f}" x2="{mx - dx3:.2f}" '
                               f'y2="{my + dy3:.2f}" stroke="{formed}"{sw}/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}"/>')
                    svg.append(f'      <line x1="{nx + dx3:.2f}" y1="{ny - dy3:.2f}" x2="{mx + dx3:.2f}" '
                               f'y2="{my - dy3:.2f}" stroke-dasharray="{dash3:.2f} {dash4:.2f}" stroke="{broken}"{sw}/>')
                elif p_order is None:
                    dx, dy = rv(double_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" x2="{mx + dx:.2f}" '
                               f'y2="{my - dy:.2f}" stroke="{broken}"{sw}/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" x2="{mx - dx:.2f}"'
                               f' y2="{my + dy:.2f}" stroke-dasharray="{dash3:.2f} {dash4:.2f}" stroke="{broken}"{sw}/>')
                else:
                    dx, dy = rv(triple_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" x2="{mx + dx:.2f}"'
                               f' y2="{my - dy:.2f}" stroke-dasharray="{dash1:.2f} {dash2:.2f}" stroke="{formed}"{sw}/>')
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"'
                               f' stroke="{broken}"{sw}/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" x2="{mx - dx:.2f}"'
                               f' y2="{my + dy:.2f}" stroke-dasharray="{dash3:.2f} {dash4:.2f}" stroke="{broken}"{sw}/>')
            elif order == 2:
                if p_order == 2:
                    dx, dy = rv(double_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}"/>')
                elif p_order == 1:
                    dx, dy = rv(double_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{broken}"{sw}/>')
                elif p_order == 4:
                    dx, dy = rv(triple_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" x2="{mx + dx:.2f}"'
                               f' y2="{my - dy:.2f}" stroke="{broken}"{sw}/>')
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" x2="{mx - dx:.2f}"'
                               f' y2="{my + dy:.2f}" stroke-dasharray="{dash3:.2f} {dash4:.2f}" stroke="{formed}"{sw}/>')
                elif p_order == 3:
                    dx, dy = rv(triple_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}"/>')
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{formed}"{sw}/>')
                elif p_order is None:
                    dx, dy = rv(double_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}" stroke="{broken}"{sw}/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{broken}"{sw}/>')
                else:
                    dx, dy = rv(triple_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}" stroke-dasharray="{dash1:.2f} '
                               f'{dash2:.2f}" stroke="{formed}"{sw}/>')
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"'
                               f' stroke="{broken}"{sw}/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{broken}"{sw}/>')
            elif order == 3:
                if p_order == 3:
                    dx, dy = rv(triple_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}"/>')
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}"/>')
                elif p_order == 1:
                    dx, dy = rv(triple_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}"  stroke="{broken}"{sw}/>')
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" 'f'stroke="{broken}"{sw}/>')
                elif p_order == 4:
                    dx, dy = rv(double_space)
                    dx3 = 3 * dx
                    dy3 = 3 * dy
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}" stroke="{broken}"{sw}/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{broken}"{sw}/>')
                    svg.append(f'      <line x1="{nx - dx3:.2f}" y1="{ny + dy3:.2f}" '
                               f'x2="{mx - dx3:.2f}" y2="{my + dy3:.2f}"/>')
                    svg.append(f'      <line x1="{nx + dx3:.2f}" y1="{ny - dy3:.2f}" x2="{mx + dx3:.2f}" '
                               f'y2="{my - dy3:.2f}" stroke-dasharray="{dash3:.2f} {dash4:.2f}" stroke="{formed}"{sw}/>')
                elif p_order == 2:
                    dx, dy = rv(triple_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}"/>')
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{broken}"{sw}/>')
                elif p_order is None:
                    dx, dy = rv(triple_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}" stroke="{broken}"{sw}/>')
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" '
                               f'x2="{mx:.2f}" y2="{my:.2f}" stroke="{broken}"{sw}/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{broken}"{sw}/>')
                else:
                    dx, dy = rv(double_space)
                    dx3 = 3 * dx
                    dy3 = 3 * dy
                    svg.append(f'      <line x1="{nx + dx3:.2f}" y1="{ny - dy3:.2f}" x2="{mx + dx3:.2f}" '
                               f'y2="{my - dy3:.2f}" stroke-dasharray="{dash1:.2f} {dash2:.2f}" stroke="{formed}"{sw}/>')
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}" stroke="{broken}"{sw}/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{broken}"{sw}/>')
                    svg.append(f'      <line x1="{nx - dx3:.2f}" y1="{ny + dy3:.2f}" x2="{mx - dx3:.2f}" '
                               f'y2="{my + dy3:.2f}" stroke="{broken}"{sw}/>')
            elif order is None:
                if p_order == 1:
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"'
                               f' stroke="{formed}"{sw}/>')
                elif p_order == 4:
                    dx, dy = rv(double_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}" stroke="{formed}"{sw}/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" x2="{mx - dx:.2f}"'
                               f' y2="{my + dy:.2f}" stroke-dasharray="{dash3:.2f} {dash4:.2f}" stroke="{formed}"{sw}/>')
                elif p_order == 2:
                    dx, dy = rv(double_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" x2="{mx + dx:.2f}" '
                               f'y2="{my - dy:.2f}" stroke="{formed}"{sw}/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" x2="{mx - dx:.2f}" '
                               f'y2="{my + dy:.2f}" stroke="{formed}"{sw}/>')
                elif p_order == 3:
                    dx, dy = rv(triple_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}" stroke="{formed}"{sw}/>')
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"'
                               f' stroke="{formed}"{sw}/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{formed}"{sw}/>')
                else:
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}" '
                               f'stroke-dasharray="{dash1:.2f} {dash2:.2f}" stroke="{formed}"{sw}/>')
            else:
                if p_order == 8:
                    svg.append(f'        <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}" '
                               f'stroke-dasharray="{dash1:.2f} {dash2:.2f}"/>')
                elif p_order == 1:
                    dx, dy = rv(double_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" x2="{mx + dx:.2f}"'
                               f' y2="{my - dy:.2f}" stroke-dasharray="{dash1:.2f} {dash2:.2f}" stroke="{broken}"{sw}/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{formed}"{sw}/>')
                elif p_order == 4:
                    dx, dy = rv(triple_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" x2="{mx + dx:.2f}"'
                               f' y2="{my - dy:.2f}" stroke-dasharray="{dash1:.2f} {dash2:.2f}" stroke="{broken}"{sw}/>')
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"'
                               f' stroke="{formed}"{sw}/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" x2="{mx - dx:.2f}"'
                               f' y2="{my + dy:.2f}" stroke-dasharray="{dash3:.2f} {dash4:.2f}" stroke="{formed}"{sw}/>')
                elif p_order == 2:
                    dx, dy = rv(triple_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" x2="{mx + dx:.2f}"'
                               f' y2="{my - dy:.2f}" stroke-dasharray="{dash1:.2f} {dash2:.2f}" stroke="{broken}"{sw}/>')
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"'
                               f' stroke="{formed}"{sw}/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{formed}"{sw}/>')
                elif p_order == 3:
                    dx, dy = rv(double_space)
                    dx3 = 3 * dx
                    dy3 = 3 * dy
                    svg.append(f'      <line x1="{nx + dx3:.2f}" y1="{ny - dy3:.2f}" x2="{mx + dx3:.2f}" '
                               f'y2="{my - dy3:.2f}" stroke-dasharray="{dash1:.2f} {dash2:.2f}" stroke="{broken}"{sw}/>')
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}" stroke="{formed}"{sw}/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{formed}"{sw}/>')
                    svg.append(f'      <line x1="{nx - dx3:.2f}" y1="{ny + dy3:.2f}" '
                               f'x2="{mx - dx3:.2f}" y2="{my + dy3:.2f}" stroke="{formed}"{sw}/>')
                else:
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}" '
                               f'stroke-dasharray="{dash1:.2f} {dash2:.2f}" stroke="{broken}"{sw}/>')
        return svg

    def _render_atoms(self):
        plane = self._plane
        config = _render_config

        carbon = config['carbon']
        font_size = config['font_size']
        other_size = config['other_size']
        monochrome = config['monochrome']
        atoms_colors = config['atoms_colors']
        other_font_style = config['other_font_style']
        dx_ci, dy_ci = config['dx_ci'], config['dy_ci']
        dx_nh, dy_nh = config['dx_nh'], config['dy_nh']
        symbols_font_style = config['symbols_font_style']

        font2 = .2 * font_size
        font4 = .4 * font_size
        font5 = .5 * font_size
        font6 = .6 * font_size
        font7 = .7 * font_size
        font15 = .15 * font_size
        font25 = .25 * font_size

        if monochrome:
            query_fill = other_fill = 'black'
        else:
            other_fill = config['other_color']
            query_fill = config['query_color']

        svg = []
        others = []
        nghbrs = []
        hbrdztns = []
        mask = defaultdict(list)
        for n, atom in self._atoms.items():
            x, y = plane[n]
            y = -y
            neighbor_svgs = [(plane[m][0], -plane[m][1]) for m in self._bonds[n]]
            single = not self._bonds[n]
            symbol = atom.atomic_symbol
            if single or symbol != 'C' or carbon or atom.charge or atom.is_radical:
                if atom.charge:
                    t = _render_charge[atom.charge]
                    tt = "↑" if atom.is_radical else ""
                    others.append(f'        <text x="{x:.2f}" y="{y:.2f}" dx="{dx_ci:.2f}" dy="-{dy_ci:.2f}" '
                                  f'font-size="{other_size:.2f}">{t}{tt}</text>')
                    mask['other'].append(f'           <text x="{x:.2f}" y="{y:.2f}" dx="{dx_ci:.2f}" dy="-{dy_ci:.2f}">'
                                         f'{t}{tt}</text>')
                elif atom.is_radical:
                    others.append(f'        <text x="{x:.2f}" y="{y:.2f}" dx="{dx_ci:.2f}" dy="-{dy_ci:.2f}" '
                                  f'font-size="{other_size:.2f}">↑</text>')
                    mask['other'].append(f'           <text x="{x:.2f}" y="{y:.2f}" dx="{dx_ci:.2f}" dy="-{dy_ci:.2f}">'
                                         f'↑</text>')

                svg.append(f'      <g fill="{"black" if monochrome else atoms_colors[atom.atomic_number - 1]}" '
                           f'font-family="{symbols_font_style}">')
                if len(symbol) > 1:
                    dx = font7
                    dx_nhh = dx_nh + font5
                    if symbol[-1] in ('l', 'i', 'r', 't'):
                        rx = font6
                        ax = font25
                    else:
                        rx = font7
                        ax = font15
                    mask['center'].append(f'          <ellipse cx="{x - ax:.2f}" cy="{y:.2f}" rx="{rx}" ry="{font4}"/>')
                else:
                    if symbol == 'I':
                        dx = font15
                        dx_nhh = dx_nh
                    else:
                        dx = font4
                        dx_nhh = dx_nh + font2
                    mask['center'].append(f'          <circle cx="{x:.2f}" cy="{y:.2f}" r="{font4:.2f}"/>')
                svg.append(f'        <text x="{x:.2f}" y="{y:.2f}" dx="-{dx:.2f}" dy="{font4:.2f}" '
                           f'font-size="{font_size:.2f}">{symbol}</text>')
                mask['symbols'].append(f'           <text x="{x:.2f}" y="{y:.2f}" dx="-{dx:.2f}" '
                                       f'dy="{font4:.2f}">{symbol}</text>')
                svg.append('      </g>')
                dx_nnh = dx_nhh
            else:
                dx_nnh = dx_nh

            if atom.neighbors or atom.hybridization:
                angle = _optimal_label_direction(x, y, neighbor_svgs)
                offset = hypot(dx_nnh, dy_nh)
                nh_dx = cos(angle) * offset
                nh_dy = sin(angle) * offset
                anchor = 'end' if nh_dx < 0 else 'start'
                level_sign = 1 if nh_dy >= 0 else -1

            if atom.neighbors:
                level = .6 * other_size
                nn = ''.join(str(x) for x in atom.neighbors)
                if atom.p_neighbors:
                    pn = ''.join(str(x) for x in atom.p_neighbors)
                else:
                    pn = '0'
                nghbrs.append(f'        <text x="{x:.2f}" y="{y:.2f}" dx="{nh_dx:.2f}" dy="{nh_dy:.2f}" '
                              f'text-anchor="{anchor}">{nn}»{pn}</text>')
            else:
                level = 0

            if atom.hybridization:
                hh = ''.join(_render_hybridization[x] for x in atom.hybridization)
                if atom.p_hybridization:
                    ph = ''.join(_render_hybridization[x] for x in atom.p_hybridization)
                else:
                    ph = '0'
                h_dy = nh_dy + level * level_sign
                hbrdztns.append(f'        <text x="{x:.2f}" y="{y:.2f}" dx="{nh_dx:.2f}" dy="{h_dy:.2f}" '
                                f'text-anchor="{anchor}">{hh}»{ph}</text>')

        if nghbrs:
            svg.append(f'      <g fill="{query_fill}" font-family="{other_font_style}" font-size="{other_size:.2f}">')
            svg.extend(nghbrs)
            if hbrdztns:
                svg.extend(hbrdztns)
            svg.append('      </g>')
        elif hbrdztns:
            svg.append(f'      <g fill="{query_fill}" font-family="{other_font_style}" font-size="{other_size:.2f}">')
            svg.extend(hbrdztns)
            svg.append('      </g>')

        if others:
            svg.append(f'      <g fill="{other_fill}" font-family="{other_font_style}" font-size="{other_size:.2f}">')
            svg.extend(others)
            svg.append('      </g>')

        return svg, mask


_render_hybridization = {1: 's', 2: 'd', 3: 't', 4: 'a'}
_render_p_charge = {-4: {-3: '-4»-3', -2: '-4»-2', -1: '-4»-', 0: '-4»0', 1: '-4»+', 2: '-4»2', 3: '-4»3', 4: '-4»4'},
                    -3: {-4: '-3»-4', -2: '-3»-2', -1: '-3»-', 0: '-3»0', 1: '-3»+', 2: '-3»2', 3: '-3»3', 4: '-3»4'},
                    -2: {-4: '-2»-4', -3: '-2»-3', -1: '-2»-', 0: '-2»0', 1: '-2»+', 2: '-2»2', 3: '-2»3', 4: '-2»4'},
                    -1: {-4: '-»-4', -3: '-»-3', -2: '-»-2', 0: '-»0', 1: '-»+', 2: '-»2', 3: '-»3', 4: '-»4'},
                    0: {-4: '0»-4', -3: '0»-3', -2: '0»-2', -1: '0»-', 1: '0»+', 2: '0»2', 3: '0»3', 4: '0»4'},
                    1: {-4: '+»-4', -3: '+»-3', -2: '+»-2', -1: '+»-', 0: '+»0', 2: '+»2', 3: '+»3', 4: '+»4'},
                    2: {-4: '2»-4', -3: '2»-3', -2: '2»-2', -1: '2»-', 0: '2»0', 1: '2»+', 3: '2»3', 4: '2»4'},
                    3: {-4: '3»-4', -3: '3»-3', -2: '3»-2', -1: '3»-', 0: '3»0', 1: '3»+', 2: '3»2', 4: '3»4'},
                    4: {-4: '4»-4', -3: '4»-3', -2: '4»-2', -1: '4»-', 0: '4»0', 1: '4»+', 2: '4»2', 3: '4»3'}}


__all__ = ['DepictMolecule', 'DepictReaction', 'DepictCGR', 'DepictQuery', 'DepictQueryCGR', 'depict_settings']
