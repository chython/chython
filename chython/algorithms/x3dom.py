# -*- coding: utf-8 -*-
#
#  Copyright 2020-2024 Ramil Nugmanov <nougmanoff@protonmail.com>
#  Copyright 2020 Dinar Batyrshin <batyrshin-dinar@mail.ru>
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
from math import acos, sqrt
from typing import TYPE_CHECKING, Union
from .depict import _render_config


if TYPE_CHECKING:
    from chython import MoleculeContainer


def plane_normal(nmx, nmy, nmz, nox, noy, noz):
    # return normal to plane of two vectors nm and no
    # m <--- n
    #         \
    #          v
    #          o
    return nmy * noz - nmz * noy, nox * nmz - nmx * noz, nmx * noy - nmy * nox


def unit_vector(nmx, nmy, nmz):
    nmd = sqrt(nmx ** 2 + nmy ** 2 + nmz ** 2)
    return nmx / nmd, nmy / nmd, nmz / nmd


def get_angle(nx, ny, nz, mx, my, mz):
    ch = (nx * mx + ny * my + nz * mz) ** 2
    zn = (nx ** 2 + ny ** 2 + nz ** 2) * (mx ** 2 + my ** 2 + mz ** 2)
    if ch < .0001:
        return 1.
    elif zn < .0001:
        return .0
    else:
        return sqrt(1 - ch / zn)


def vector_normal(nmx, nmy, nmz):
    # return normal to vector nm
    if not -.0001 < nmx < .0001:
        return (- nmy - nmz) / nmx, 1, 1
    elif not -.0001 < nmy < .0001:
        return 1, (- nmx - nmz) / nmy, 1
    else:
        return 1, 1, (- nmx - nmy) / nmz


class JupyterWidget:
    def __init__(self, xml, width, height):
        self.xml = xml
        self.width = width
        self.height = height

    def _repr_html_(self):
        return ("<script type='text/javascript' src='https://www.x3dom.org/download/x3dom.js'></script>"
                "<link rel='stylesheet' type='text/css' href='https://www.x3dom.org/download/x3dom.css'>"
                f'<div style="width: {self.width}; height: {self.height}">{self.xml}</div>')

    def __html__(self):
        return self._repr_html_()


def _render_aromatic_bond(n_x, n_y, n_z, m_x, m_y, m_z, c_x, c_y, c_z):
    aromatic_space = _render_config['aromatic_space']

    # n aligned xyz
    nc_x, nc_y, nc_z = c_x - n_x, c_y - n_y, c_z - n_z
    mc_x, mc_y, mc_z = c_x - m_x, c_y - m_y, c_z - m_z

    nc_ln = sqrt(nc_x ** 2 + nc_y ** 2 + nc_z ** 2)
    mc_ln = sqrt(mc_x ** 2 + mc_y ** 2 + mc_z ** 2)
    sin1 = get_angle(m_x - n_x, m_y - n_y, m_z - n_z, nc_x, nc_y, nc_z)
    sin2 = get_angle(n_x - m_x, n_y - m_y, n_z - m_z, mc_x, mc_y, mc_z)

    if sin1 < .0001 or sin2 < .0001 or nc_ln < .0001 or mc_ln < .0001:
        return
    else:
        coef1 = aromatic_space / (nc_ln * sin1)
        coef2 = aromatic_space / (mc_ln * sin2)
        return nc_x * coef1, nc_y * coef1, nc_z * coef1, mc_x * coef2, mc_y * coef2, mc_z * coef2


def _render_dashes(nx, ny, nz, nmx, nmy, nmz, nm_ln, r_angle=None):
    bond_radius = _render_config['bond_radius']
    bond_color = _render_config['bond_color']

    if r_angle is None:
        dash1, dash2 = _render_config['aromatic_dashes']
        r_angle = acos(nmy / nm_ln)
    else:
        dash1, dash2 = _render_config['dashes']

    xml = []
    dashes_sum = dash1 + dash2
    if dashes_sum < .0001:
        raise ValueError('Dashes should be nonzero')

    d = dashes_sum / nm_ln
    dx, dy, dz = nmx * d, nmy * d, nmz * d
    b = int((nm_ln - dash1) // dashes_sum)
    t = (nm_ln - (b * dashes_sum)) / nm_ln
    nx, ny, nz = nx + nmx * t / 2, ny + nmy * t / 2, nz + nmz * t / 2
    for _ in range(b):
        xml.append(f"    <transform translation='{nx:.2f} {ny:.2f} {nz:.2f}' rotation='{nmz:.2f} 0 "
                   f"{-nmx:.2f} {r_angle:.2f}'>\n      <shape>\n        <appearance>\n"
                   f"          <material diffusecolor='{bond_color}'>\n          </material>\n"
                   f"       </appearance>\n        <cylinder radius='{bond_radius}' height='{dash1:.2f}'>\n"
                   "        </cylinder>\n      </shape>\n    </transform>\n")
        nx += dx
        ny += dy
        nz += dz
    xml.append(f"    <transform translation='{nx:.2f} {ny:.2f} {nz:.2f}' rotation='{nmz:.2f} 0 "
               f"{-nmx:.2f} {r_angle:.2f}'>\n      <shape>\n        <appearance>\n"
               f"          <material diffusecolor='{bond_color}'>\n          </material>\n"
               f"       </appearance>\n        <cylinder radius='{bond_radius}' height='{dash1:.2f}'>\n"
               "        </cylinder>\n      </shape>\n    </transform>\n")
    return xml


class X3domMolecule:
    __slots__ = ()

    def depict3d(self: Union['MoleculeContainer', 'X3domMolecule'], index: int = 0) -> str:
        """Get X3DOM XML string.

        :param index: index of conformer
        """
        if not hasattr(self, '_conformers'):
            raise ValueError('No conformers stored within structure')
        try:
            xyz = self._conformers[index]
        except IndexError:
            raise IndexError('Invalid conformer index')

        mx = sum(x for x, _, _ in xyz.values()) / len(xyz)
        my = sum(y for _, y, _ in xyz.values()) / len(xyz)
        mz = sum(z for _, _, z in xyz.values()) / len(xyz)
        xyz = {n: (x - mx, y - my, z - mz) for n, (x, y, z) in xyz.items()}
        atoms = self.__render_atoms(xyz)
        bonds = self.__render_bonds(xyz)
        return f'<x3d width=100% height=100%>\n  <scene>\n{atoms}{bonds}  </scene>\n</x3d>'

    def view3d(self, index: int = 0, width='600px', height='400px'):
        """
        Jupyter widget for 3D visualization.

        :param index: index of conformer
        :param width: widget width
        :param height: widget height
        """
        return JupyterWidget(self.depict3d(index), width, height)

    def __render_atoms(self: 'MoleculeContainer', xyz):
        font = _render_config['font_size']
        carbon = _render_config['carbon']
        radius = _render_config['atom_radius']
        colors = _render_config['atoms_colors']
        mapping_color = _render_config['mapping_color']

        if radius < 0:
            multiplier = -radius
            radius = 0
        elif not radius:
            multiplier = .2

        atoms = []
        if carbon:
            for n, a in self.atoms():
                r = radius or a.atomic_radius * multiplier
                fr = r * 0.71
                atoms.append(f"    <transform translation='{' '.join(format(x, '.2f') for x in xyz[n])}'>\n"
                             "      <billboard axisOfRotation='0 0 0'>\n        <group>\n"
                             f"          <transform translation='{fr:.2f} {fr:.2f} 0'>\n"
                             '            <shape>\n              <appearance>\n'
                             f"                <material diffuseColor='{mapping_color}'/>\n"
                             f"              </appearance>\n              <text string='{a.atomic_symbol}'>\n"
                             f"                <fontstyle family='sans' size='{font:.2f}' justify='first'/>\n"
                             "              </text>\n            </shape>\n          </transform>\n"
                             "          <shape>\n            <appearance>\n"
                             f"              <material diffuseColor='{colors[a.atomic_number - 1]}'/>\n"
                             f"            </appearance>\n            <sphere radius='{r:.2f}'/>\n"
                             "          </shape>\n        </group>\n      </billboard>\n    </transform>\n")
        else:
            for n, a in self.atoms():
                r = radius or a.atomic_radius * multiplier
                atoms.append(f"    <transform translation='{' '.join(format(x, '.2f') for x in xyz[n])}'>\n"
                             "      <shape>\n        <appearance>\n"
                             f"          <material diffuseColor='{colors[a.atomic_number - 1]}'/>\n"
                             f"        </appearance>\n        <sphere radius='{r:.2f}'/>\n"
                             "      </shape>\n    </transform>\n")
        return ''.join(atoms)

    def __render_bonds(self: 'MoleculeContainer', xyz):
        bonds = self._bonds

        bond_color = _render_config['bond_color']
        bond_radius = _render_config['bond_radius']
        double_space = _render_config['double_space']
        triple_space = _render_config['triple_space']
        r1 = triple_space * sqrt(3) / 3
        r2 = triple_space * sqrt(3) / 6

        xml = []
        lengths = {}
        doubles = {}
        half_triple = triple_space / 2
        for n, m, bond in self.bonds():
            nx, ny, nz = xyz[n]
            mx, my, mz = xyz[m]

            nmx, nmy, nmz = mx - nx, my - ny, mz - nz
            length = sqrt(nmx ** 2 + nmy ** 2 + nmz ** 2)
            if length < .001:
                continue

            rotation_angle = acos(nmy / length)
            lengths[(n, m)] = lengths[(m, n)] = (length, rotation_angle)
            x, y, z = nx + nmx / 2, ny + nmy / 2, nz + nmz / 2
            if bond in (1, 4):
                xml.append(f"    <transform translation='{x:.2f} {y:.2f} {z:.2f}' rotation='{nmz:.2f} 0 "
                           f"{-nmx:.2f} {rotation_angle:.2f}'>\n      <shape>\n        <appearance>\n"
                           f"          <material diffusecolor='{bond_color}'>\n          </material>\n"
                           f"       </appearance>\n        <cylinder radius='{bond_radius}' height='{length:.2f}'>\n"
                           "        </cylinder>\n      </shape>\n    </transform>\n")
            elif bond == 2:
                if n in doubles:
                    # normal for plane n m o
                    norm_x, norm_y, norm_z = plane_normal(nmx, nmy, nmz, *doubles[n])
                elif m in doubles:
                    # normal for plane n m o
                    norm_x, norm_y, norm_z = plane_normal(nmx, nmy, nmz, *doubles[m])
                else:
                    third = next((x for x in bonds[n] if x != m), None)
                    if third:
                        ox, oy, oz = xyz[third]
                        nox, noy, noz = ox - nx, oy - ny, oz - nz
                    else:
                        third = next((x for x in bonds[m] if x != n), None)
                        if third:
                            ox, oy, oz = xyz[third]
                            nox, noy, noz = ox - nx, oy - ny, oz - nz
                        else:
                            nox, noy, noz = vector_normal(nmx, nmy, nmz)

                    # normal for plane n m o
                    normx, normy, normz = unit_vector(*plane_normal(nmx, nmy, nmz, nox, noy, noz))

                    # normal for plane n m normal
                    norm_x, norm_y, norm_z = plane_normal(nmx, nmy, nmz, normx, normy, normz)

                doubles[n] = doubles[m] = (norm_x, norm_y, norm_z)
                norm_dist = sqrt(norm_x ** 2 + norm_y ** 2 + norm_z ** 2)

                if norm_dist < .0001:
                    coef = double_space * 10000
                else:
                    coef = double_space / norm_dist

                dx, dy, dz = norm_x * coef, norm_y * coef, norm_z * coef
                xml.append(
                    f"    <transform translation='{x + dx:.2f} {y + dy:.2f} {z + dz:.2f}' rotation='{nmz:.2f} 0 "
                    f"{-nmx:.2f} {rotation_angle:.2f}'>\n      <shape>\n        <appearance>\n"
                    f"          <material diffusecolor='{bond_color}'>\n          </material>\n"
                    f"       </appearance>\n        <cylinder radius='{bond_radius}' height='{length:.2f}'>\n"
                    "        </cylinder>\n      </shape>\n    </transform>\n")
                xml.append(
                    f"    <transform translation='{x - dx:.2f} {y - dy:.2f} {z - dz:.2f}' rotation='{nmz:.2f} 0 "
                    f"{-nmx:.2f} {rotation_angle:.2f}'>\n      <shape>\n        <appearance>\n"
                    f"          <material diffusecolor='{bond_color}'>\n          </material>\n"
                    f"       </appearance>\n        <cylinder radius='{bond_radius}' height='{length:.2f}'>\n"
                    "        </cylinder>\n      </shape>\n    </transform>\n")
            elif bond == 3:
                nox, noy, noz = vector_normal(nmx, nmy, nmz)

                # normal for plane n m o
                normx, normy, normz = unit_vector(*plane_normal(nmx, nmy, nmz, nox, noy, noz))
                vecrx, vecry, vecrz = normx * r1, normy * r1, normz * r1

                # normal for plane n m normal
                norm_x, norm_y, norm_z = unit_vector(*plane_normal(nmx, nmy, nmz, normx, normy, normz))
                vecx, vecy, vecz = norm_x * half_triple, norm_y * half_triple, norm_z * half_triple

                xml.append(f"    <transform translation='{x + vecrx:.2f} {y + vecry:.2f} {z + vecrz:.2f}'"
                           f" rotation='{nmz:.2f} 0 {-nmx:.2f} {rotation_angle:.2f}'>\n      <shape>\n"
                           f"        <appearance>\n          <material diffusecolor='{bond_color}'>\n"
                           f"          </material>\n       </appearance>\n        <cylinder radius='{bond_radius}'"
                           f" height='{length:.2f}'>\n        </cylinder>\n      </shape>\n    </transform>\n")

                xx, yy, zz = x - normx * r2, y - normy * r2, z - normz * r2
                xml.append(f"    <transform translation='{xx - vecx:.2f} {yy - vecy:.2f} {zz - vecz:.2f}'"
                           f" rotation='{nmz:.2f} 0 {-nmx:.2f} {rotation_angle:.2f}'>\n      <shape>\n"
                           f"        <appearance>\n          <material diffusecolor='{bond_color}'>\n"
                           f"          </material>\n       </appearance>\n        <cylinder radius='{bond_radius}'"
                           f" height='{length:.2f}'>\n        </cylinder>\n      </shape>\n    </transform>\n")
                xml.append(f"    <transform translation='{xx + vecx:.2f} {yy + vecy:.2f} {zz + vecz:.2f}'"
                           f" rotation='{nmz:.2f} 0 {-nmx:.2f} {rotation_angle:.2f}'>\n      <shape>\n"
                           f"        <appearance>\n          <material diffusecolor='{bond_color}'>\n"
                           f"          </material>\n       </appearance>\n        <cylinder radius='{bond_radius}'"
                           f" height='{length:.2f}'>\n        </cylinder>\n      </shape>\n    </transform>\n")
            else:
                xml.extend(_render_dashes(nx, ny, nz, nmx, nmy, nmz, length, r_angle=rotation_angle))

        for ring in self.aromatic_rings:
            cx = sum(xyz[n][0] for n in ring) / len(ring)
            cy = sum(xyz[n][1] for n in ring) / len(ring)
            cz = sum(xyz[n][2] for n in ring) / len(ring)

            for n, m in zip(ring, ring[1:]):
                nx, ny, nz = xyz[n]
                mx, my, mz = xyz[m]

                aromatic = _render_aromatic_bond(nx, ny, nz, mx, my, mz, cx, cy, cz)
                if aromatic:
                    veca_x, veca_y, veca_z, vecb_x, vecb_y, vecb_z = aromatic
                    ax, ay, az = nx + veca_x, ny + veca_y, nz + veca_z
                    abx, aby, abz = mx + vecb_x - ax, my + vecb_y - ay, mz + vecb_z - az
                    ab_ln = sqrt(abx ** 2 + aby ** 2 + abz ** 2)
                    if ab_ln < .0001:
                        continue
                    else:
                        xml.extend(_render_dashes(ax, ay, az, abx, aby, abz, ab_ln))

            i, j = ring[-1], ring[0]
            nx, ny, nz = xyz[i]
            mx, my, mz = xyz[j]
            aromatic = _render_aromatic_bond(nx, ny, nz, mx, my, mz, cx, cy, cz)
            if aromatic:
                veca_x, veca_y, veca_z, vecb_x, vecb_y, vecb_z = aromatic
                ax, ay, az = nx + veca_x, ny + veca_y, nz + veca_z
                abx, aby, abz = mx + vecb_x - ax, my + vecb_y - ay, mz + vecb_z - az
                ab_ln = sqrt(abx ** 2 + aby ** 2 + abz ** 2)
                if ab_ln < .0001:
                    continue
                else:
                    xml.extend(_render_dashes(ax, ay, az, abx, aby, abz, ab_ln))
        return ''.join(xml)


__all__ = ['X3domMolecule']
