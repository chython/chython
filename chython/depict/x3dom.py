# -*- coding: utf-8 -*-
#
#  Copyright 2020-2026 Ramil Nugmanov <nougmanoff@protonmail.com>
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
"""X3DOM output for a stored conformer: `mol.depict3d()` and `mol.view3d()`.

The bodies behind those two container methods, registered by `_hooks.register()`.  A model is READ --
`molecule.conformer(index)` -- and never computed, so a molecule with no geometry is refused here
rather than drawn against an invented one.  Rendering parameters are the module-level
`_X3DOM_DEFAULTS` below and are not part of `DepictStyle`, which is 2D throughout.

Spheres are sized from `Atom.atomic_radius`, the calculated radius; an R carries none, so the marker is
drawn as its own label in `R_COLOUR` the way the 2D side draws it.
"""
from math import acos, fsum, sqrt
from types import MappingProxyType
from ._config import R_COLOUR, cpk


#: This module's rendering parameters, frozen so a reader cannot mistake them for live settings.
#: `atoms_colors` references `_config.cpk`, the tree's one CPK palette, rather than copying 118 strings.
#: `atom_radius` NEGATIVE means a multiplier on `Atom.atomic_radius` rather than a fixed sphere size,
#: which is what makes an oxygen smaller than a carbon; a positive value is that size in angstroms.
_X3DOM_DEFAULTS = MappingProxyType({
    'carbon': False, 'dashes': (.2, .1), 'font_size': .5, 'bond_color': 'black', 'bond_radius': .02,
    'atom_radius': -.2, 'atoms_colors': cpk, 'triple_space': .13, 'double_space': .06,
    'mapping_color': '#0305A7', 'marker_color': R_COLOUR, 'aromatic_space': .14,
    'aromatic_dashes': (.15, .05)})


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
    aromatic_space = _X3DOM_DEFAULTS['aromatic_space']

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
    bond_radius = _X3DOM_DEFAULTS['bond_radius']
    bond_color = _X3DOM_DEFAULTS['bond_color']

    if r_angle is None:
        dash1, dash2 = _X3DOM_DEFAULTS['aromatic_dashes']
        r_angle = acos(nmy / nm_ln)
    else:
        dash1, dash2 = _X3DOM_DEFAULTS['dashes']

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


def _render_atoms(molecule, xyz):
    """Every atom as a sphere sized from its calculated radius; the R marker as its own label."""
    font = _X3DOM_DEFAULTS['font_size']
    labelled = _X3DOM_DEFAULTS['carbon']
    radius = _X3DOM_DEFAULTS['atom_radius']
    colors = _X3DOM_DEFAULTS['atoms_colors']
    label_color = _X3DOM_DEFAULTS['mapping_color']
    marker_color = _X3DOM_DEFAULTS['marker_color']

    # A negative `atom_radius` is a MULTIPLIER on the atom's own radius; a positive one is that size
    # for every atom.  See `_X3DOM_DEFAULTS`.
    if radius < 0:
        multiplier = -radius
        radius = 0.
    else:
        multiplier = .2

    atoms = []
    for atom in molecule.atoms():
        x, y, z = xyz[atom.n]
        if atom.element == 0:
            # No radius, so no sphere: a marker with a zero-radius sphere is a bond ending in nothing.
            atoms.append(_render_label(x, y, z, atom.atomic_symbol, marker_color, font))
            continue
        r = radius or atom.atomic_radius * multiplier
        colour = colors[atom.element - 1]
        atoms.append(_render_sphere(x, y, z, r, colour))
        if labelled:
            atoms.append(_render_label(x + r * .71, y + r * .71, z, atom.atomic_symbol,
                                       label_color, font))
    return ''.join(atoms)


def _render_sphere(x, y, z, r, colour):
    return (f"    <transform translation='{x:.2f} {y:.2f} {z:.2f}'>\n"
            "      <shape>\n        <appearance>\n"
            f"          <material diffuseColor='{colour}'/>\n"
            f"        </appearance>\n        <sphere radius='{r:.2f}'/>\n"
            "      </shape>\n    </transform>\n")


def _render_label(x, y, z, text, colour, font):
    """Text that turns to face the camera -- a billboard, so a label is readable from any angle."""
    return (f"    <transform translation='{x:.2f} {y:.2f} {z:.2f}'>\n"
            "      <billboard axisOfRotation='0 0 0'>\n        <shape>\n          <appearance>\n"
            f"            <material diffuseColor='{colour}'/>\n"
            f"          </appearance>\n          <text string='{text}'>\n"
            f"            <fontstyle family='sans' size='{font:.2f}' justify='middle'/>\n"
            "          </text>\n        </shape>\n      </billboard>\n    </transform>\n")


def _render_bonds(molecule, xyz):
    """Every bond as cylinders: one, two offset, three, or dashes for an order nothing pins."""
    bond_color = _X3DOM_DEFAULTS['bond_color']
    bond_radius = _X3DOM_DEFAULTS['bond_radius']
    double_space = _X3DOM_DEFAULTS['double_space']
    triple_space = _X3DOM_DEFAULTS['triple_space']
    r1 = triple_space * sqrt(3) / 3
    r2 = triple_space * sqrt(3) / 6

    xml = []
    doubles = {}
    half_triple = triple_space / 2
    for bond in molecule.bonds():
        n, m, order = bond.n, bond.m, int(bond)
        nx, ny, nz = xyz[n]
        mx, my, mz = xyz[m]

        nmx, nmy, nmz = mx - nx, my - ny, mz - nz
        length = sqrt(nmx ** 2 + nmy ** 2 + nmz ** 2)
        if length < .001:
            continue

        rotation_angle = acos(nmy / length)
        x, y, z = nx + nmx / 2, ny + nmy / 2, nz + nmz / 2
        if order in (1, 4):
            xml.append(f"    <transform translation='{x:.2f} {y:.2f} {z:.2f}' rotation='{nmz:.2f} 0 "
                       f"{-nmx:.2f} {rotation_angle:.2f}'>\n      <shape>\n        <appearance>\n"
                       f"          <material diffusecolor='{bond_color}'>\n          </material>\n"
                       f"       </appearance>\n        <cylinder radius='{bond_radius}' height='{length:.2f}'>\n"
                       "        </cylinder>\n      </shape>\n    </transform>\n")
        elif order == 2:
            if n in doubles:
                # normal for plane n m o
                norm_x, norm_y, norm_z = plane_normal(nmx, nmy, nmz, *doubles[n])
            elif m in doubles:
                # normal for plane n m o
                norm_x, norm_y, norm_z = plane_normal(nmx, nmy, nmz, *doubles[m])
            else:
                third = next((k for k in molecule.neighbors_of(n) if k != m), None)
                if third:
                    ox, oy, oz = xyz[third]
                    nox, noy, noz = ox - nx, oy - ny, oz - nz
                else:
                    third = next((k for k in molecule.neighbors_of(m) if k != n), None)
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
        elif order == 3:
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

    for ring in molecule.aromatic_rings:
        cx = fsum(xyz[n][0] for n in ring) / len(ring)
        cy = fsum(xyz[n][1] for n in ring) / len(ring)
        cz = fsum(xyz[n][2] for n in ring) / len(ring)

        for n, m in zip(ring, ring[1:]):
            nx, ny, nz = xyz[n]
            mx, my, mz = xyz[m]

            aromatic = _render_aromatic_bond(nx, ny, nz, mx, my, mz, cx, cy, cz)
            if aromatic:
                veca_x, veca_y, veca_z, vecb_x, vecb_y, vecb_z = aromatic
                ax, ay, az = nx + veca_x, ny + veca_y, nz + veca_z
                abx, aby, abz = mx + vecb_x - ax, my + vecb_y - ay, mz + vecb_z - az
                ab_ln = sqrt(abx ** 2 + aby ** 2 + abz ** 2)
                if ab_ln >= .0001:
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
            if ab_ln >= .0001:
                xml.extend(_render_dashes(ax, ay, az, abx, aby, abz, ab_ln))
    return ''.join(xml)


def molecule_depict3d(molecule, index=0):
    """Model `index` of `molecule` as an X3DOM document.  `MoleculeContainer.depict3d`'s body.

    The model is CENTRED on its own centroid, so a conformer read out of a crystal file is drawn at the
    origin rather than off-screen; nothing is stored, the way no drawing function stores anything.
    """
    if not molecule.has_3d:
        raise ValueError('no conformer stored on this molecule, and a conformer is read rather than '
                         'guessed: a 2D layout follows from the graph and a geometry does not.  '
                         '`chython.interop.conformers.generate_conformers` is the call that makes one')
    conformer = molecule.conformer(index)      # IndexError names the model that is not there

    xyz = {n: conformer.xyz_of(n) for n in molecule.atom_numbers}
    mx = fsum(x for x, _, _ in xyz.values()) / len(xyz)
    my = fsum(y for _, y, _ in xyz.values()) / len(xyz)
    mz = fsum(z for _, _, z in xyz.values()) / len(xyz)
    xyz = {n: (x - mx, y - my, z - mz) for n, (x, y, z) in xyz.items()}

    atoms = _render_atoms(molecule, xyz)
    bonds = _render_bonds(molecule, xyz)
    return f'<x3d width=100% height=100%>\n  <scene>\n{atoms}{bonds}  </scene>\n</x3d>'


def molecule_view3d(molecule, index=0, width='600px', height='400px'):
    """Model `index` in a Jupyter widget.  `MoleculeContainer.view3d`'s body, and `depict3d` in a div."""
    return JupyterWidget(molecule_depict3d(molecule, index), width, height)


__all__ = ['molecule_depict3d', 'molecule_view3d', 'JupyterWidget']
