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
"""Predefined ring-system shapes: pseudo-3D perspectives of bridged bicycles and cages.

| part | rule |
| --- | --- |
| skeleton | a ring system's atoms and bonds as an all-carbon, all-single-bond graph; a shape matches a system whose \
skeleton has its canonical form, whatever the elements and orders |
| model | `MODELS`: a skeleton's carbon frame from an MMFF-relaxed ETKDG conformer, mean bond 1, in its SMILES atom \
order, and its views best first, each seen orthographically, the model's z axis up |
| choice | per molecule, the drawing and symmetry whose reference substitutes most of the atoms the molecule \
substitutes, then leaving fewest of them inside its outline, drawn shapes first |
| drawn | `DRAWN`: a reference molecule and its ring atoms' points read from a drawing of it, mean bond 1; several \
may share a skeleton, and they replace the model's views |
| shape | `SHAPES`: points in `L`: the drawn skeletons, each other model seen along its first view |
| drawing | `shape(mol, atoms)`: the matched shape and its mirror image in plane coordinates, or None |

| shape | skeletons it draws |
| --- | --- |
| norbornane | camphor, norbornene |
| bicyclo[2.2.2]octane | quinuclidine, DABCO |
| bicyclo[3.2.1]octane | tropane |
| bicyclo[3.3.1]nonane | granatane |
| bicyclo[1.1.1]pentane, bicyclo[2.1.1]hexane, adamantane, twistane, cubane, prismane, tetrahedrane | |
"""
from functools import lru_cache
from itertools import islice
from math import cos, pi, sin

from ...core import MoleculeContainer, smiles
from .lattice import L


__all__ = ['DRAWN', 'MODELS', 'SHAPES', 'skeleton', 'shape']

EDGE = .08               # an atom this near its drawing's outline, in bonds, is on it
MAPPINGS = 48            # skeleton symmetries tried for a molecule's substituents


MODELS = {
    'bicyclo[1.1.1]pentane': ('C1C2CC1C2', (
            (0.0915, -0.9599, -0.2650),
            (0.7218, 0.4350, -0.5383),
            (-0.9923, -0.0799, 0.0950),
            (-0.7050, 0.7048, -0.0783),
            (-0.4622, -0.8643, 0.1983),
            (0.4735, 0.4423, 0.7617),
            (0.3396, -0.8640, 0.3717),
            (0.0208, 0.5696, 0.8217),
        ), [
            (-0.4535, 0.6337, -0.0706), (-0.0458, -0.1010, -0.6128), (0.7767, 0.0651, -0.0687),
            (0.0458, 0.1010, 0.6128), (-0.3232, -0.6988, 0.1394),
        ]),
    'bicyclo[2.1.1]hexane': ('C1CC2CC1C2', (
            (0.7218, 0.4350, -0.5383),
            (0.7050, -0.3532, 0.6150),
            (-0.3388, 0.6115, -0.7150),
            (-0.6013, -0.4429, -0.6650),
            (0.3337, 0.8678, 0.3683),
            (-0.2328, -0.2789, 0.9317),
            (0.6775, 0.7289, 0.0983),
            (-0.7656, 0.2279, 0.6017),
        ), [
            (-0.8145, 0.3577, 0.2226), (-0.7224, -0.4606, -0.3270), (0.2620, -0.5267, -0.3491),
            (0.5545, -0.3395, 0.5983), (0.1404, 0.5536, 0.3765), (0.5800, 0.4154, -0.5213),
        ]),
    'norbornane': ('C1CC2CCC1C2', (
            (0.5296, -0.6636, 0.5283),
            (0.8148, 0.5269, -0.2417),
            (0.0208, 0.5696, 0.8217),
            (0.3917, -0.6428, -0.6583),
        ), [
            (-0.8560, -0.5198, -0.1375), (-0.6698, 0.2858, -0.7011), (0.1136, 0.7066, -0.2395),
            (0.9227, 0.1353, -0.3901), (0.7365, -0.6703, 0.1736), (-0.1562, -0.4605, 0.5771),
            (-0.0907, 0.5229, 0.7174),
        ]),
    'bicyclo[2.2.2]octane': ('C1CC2CCC1CC2', (
            (0.9751, -0.0118, 0.2217),
            (0.1368, 0.9010, 0.4117),
            (0.1996, 0.4904, 0.8483),
            (-0.4622, -0.8643, 0.1983),
            (-0.6065, 0.0087, -0.7950),
            (0.8649, 0.4986, -0.0583),
            (-0.6380, -0.2401, 0.7317),
            (-0.3204, 0.3464, 0.8817),
        ), [
            (-0.8090, 0.4619, -0.5124), (-0.1028, 1.0511, -0.1228), (0.5983, 0.4991, 0.3300),
            (0.1633, -0.0432, 1.0497), (-0.5429, -0.6323, 0.6601), (-0.5983, -0.4991, -0.3300),
            (0.2926, -0.7133, -0.7321), (0.9988, -0.1242, -0.3425),
        ]),
    'bicyclo[3.2.1]octane': ('C1CC2CCC(C1)C2', (
            (0.7903, -0.3104, -0.5283),
            (-0.4614, 0.2368, 0.8550),
            (-0.7179, -0.1037, 0.6883),
            (-0.7098, -0.7031, 0.0417),
            (0.7448, -0.5427, 0.3883),
            (0.1424, -0.5617, -0.8150),
        ), [
            (1.1318, 0.0290, -0.1793), (0.6330, 0.8924, -0.1192), (-0.3600, 0.7659, -0.0520),
            (-0.7150, 0.1989, -0.7972), (-0.6391, -0.7455, -0.4712), (-0.2440, -0.6782, 0.4466),
            (0.7589, -0.6760, 0.4221), (-0.5656, 0.2135, 0.7501),
        ]),
    'bicyclo[3.3.1]nonane': ('C1CC2CCCC(C1)C2', (
            (-0.4248, 0.1120, 0.8983),
            (-0.2069, -0.0076, -0.9783),
            (-0.0898, 0.4937, 0.8650),
            (0.1979, 0.3047, -0.9317),
        ), [
            (0.9933, -0.5129, -0.5350), (0.4658, -1.1278, 0.0460), (-0.3485, -0.6921, 0.4390),
            (-1.0268, -0.4715, -0.2697), (-1.1749, 0.5055, -0.3958), (-0.3299, 1.0325, -0.4177),
            (0.3303, 0.7727, 0.2949), (1.1595, 0.3692, -0.1013), (-0.0689, 0.1244, 0.9396),
        ]),
    'adamantane': ('C1C2CC3CC1CC(C2)C3', (
            (-0.5535, 0.7750, -0.3050),
            (0.8837, 0.4510, 0.1250),
            (0.3042, -0.4462, -0.8417),
            (-0.2970, 0.9547, 0.0183),
            (-0.0902, 0.8584, -0.5050),
            (-0.9485, 0.0321, -0.3150),
            (-0.1699, -0.0301, 0.9850),
            (0.6125, -0.1174, -0.7817),
        ), [
            (-1.0692, -0.3320, -0.2737), (-0.8396, 0.4961, 0.2378), (-0.3493, 0.2417, 1.0714),
            (0.4883, -0.2530, 0.8397), (0.2512, -1.0769, 0.3248), (-0.2357, -0.8300, -0.5130),
            (0.3493, -0.2417, -1.0714), (0.5869, 0.5869, -0.5645), (-0.2512, 1.0769, -0.3248),
            (1.0692, 0.3320, 0.2737),
        ]),
    'cubane': ('C12C3C4C1C5C2C3C45', (
            (-0.1872, 0.5770, 0.7950),
            (0.8033, -0.0922, 0.5883),
            (-0.3591, -0.6410, 0.6783),
            (0.6712, -0.7289, -0.1350),
            (-0.7049, -0.3861, -0.5950),
            (0.2572, -0.9656, -0.0383),
            (-0.0736, 0.0675, 0.9950),
            (-0.7511, -0.4551, 0.4783),
        ), [
            (0.4733, 0.5917, -0.4194), (0.7381, 0.1380, 0.4315), (-0.1808, 0.2869, 0.7969),
            (-0.4456, 0.7406, -0.0541), (-0.7381, -0.1380, -0.4315), (0.1808, -0.2870, -0.7969),
            (0.4456, -0.7406, 0.0541), (-0.4733, -0.5917, 0.4194),
        ]),
    'prismane': ('C12C3C1C1C2C31', (
            (0.7122, -0.0225, 0.7017),
            (0.4119, -0.8633, -0.2917),
            (0.8257, 0.5113, 0.2383),
            (0.7218, 0.4350, -0.5383),
            (-0.4651, 0.7802, -0.4183),
            (0.5192, -0.2231, -0.8250),
            (-0.3433, -0.0888, -0.9350),
            (-0.4835, -0.7812, -0.3950),
        ), [
            (0.5496, -0.0418, 0.5286), (0.3259, -0.6499, -0.2336), (0.5679, 0.3047, -0.4097),
            (-0.3944, 0.5627, -0.3332), (-0.4127, 0.2162, 0.6051), (-0.6363, -0.3919, -0.1571),
        ]),
    'tetrahedrane': ('C12C3C1C23', (
            (-0.3844, -0.1226, 0.9150),
            (-0.8878, 0.3587, -0.2883),
            (0.3872, 0.9192, 0.0717),
            (-0.5980, 0.5750, -0.5583),
            (0.1463, -0.2580, -0.9550),
            (-0.7894, 0.0517, -0.6117),
            (0.7356, -0.6770, -0.0217),
            (-0.9998, 0.0023, -0.0183),
        ), [
            (0.5146, 0.2167, -0.2549), (-0.0392, -0.4321, 0.4342), (-0.2667, -0.2791, -0.4772),
            (-0.2086, 0.4945, 0.2979),
        ]),
}


def _inside(pts):
    """The indices of `pts` strictly inside their convex hull, farther than `EDGE` from it."""
    order = sorted(range(len(pts)), key=lambda i: pts[i])

    def half(seq):
        h = []
        for i in seq:
            while len(h) > 1 and (pts[h[-1]][0] - pts[h[-2]][0]) * (pts[i][1] - pts[h[-2]][1]) - \
                    (pts[h[-1]][1] - pts[h[-2]][1]) * (pts[i][0] - pts[h[-2]][0]) <= 1e-9:
                h.pop()
            h.append(i)
        return h
    hull = half(order)[:-1] + half(order[::-1])[:-1]
    out = set()
    for i, (x, y) in enumerate(pts):
        if all((pts[b][0] - pts[a][0]) * (y - pts[a][1]) - (pts[b][1] - pts[a][1]) * (x - pts[a][0]) >= EDGE *
               ((pts[b][0] - pts[a][0]) ** 2 + (pts[b][1] - pts[a][1]) ** 2) ** .5
               for a, b in zip(hull, hull[1:] + hull[:1])):
            out.add(i)
    return out


def _seen(xyz, view):
    """`xyz` seen orthographically along `view`, the model's z axis up where it is not the view."""
    n = sum(v * v for v in view) ** .5
    z = [v / n for v in view]
    up = (0., 0., 1.) if abs(z[2]) < .95 else (0., 1., 0.)
    x = (up[1] * z[2] - up[2] * z[1], up[2] * z[0] - up[0] * z[2], up[0] * z[1] - up[1] * z[0])
    m = sum(v * v for v in x) ** .5
    x = [v / m for v in x]
    y = (z[1] * x[2] - z[2] * x[1], z[2] * x[0] - z[0] * x[2], z[0] * x[1] - z[1] * x[0])
    return [(sum(a * b for a, b in zip(p, x)), sum(a * b for a, b in zip(p, y))) for p in xyz]


DRAWN = {                # name: (reference SMILES, its ring atoms' points in their SMILES order, the drawing read)
    'norbornane': ('C1CC2CCC1C2', [
        (0.66, -0.711), (1.176, 0.081), (0.276, 0.38), (-0.624, 0.189), (-1.164, -0.615),
        (-0.264, -0.423), (-0.06, 1.1),
    ], 'Wikimedia Commons NorbornaneNumbering.svg'),
    'camphor': ('CC1(C)C2CCC1(C)C(=O)C2', [
        (-0.196, 1.159), (0.254, 0.369), (-0.608, 0.122), (-1.086, -0.675), (-0.196, -0.433),
        (0.682, -0.664), (1.149, 0.122),
    ], 'Wikimedia Commons Alcanfor-2D.svg'),
    'bicyclo[1.1.1]pentane': ('C1C2CC1C2', [
        (-0.512, 0.845), (-0.185, 0.082), (-0.471, -0.858), (0.169, -0.109), (1.0, 0.041),
    ], 'Wikimedia Commons Bicyclo-1.1.1-pentane.svg'),
    'adamantane': ('C1C2CC3CC1CC(C2)C3', [
        (-1.264, -0.973), (-0.765, -0.106), (0.206, -0.383), (1.163, -0.113), (0.663, -0.973),
        (-0.3, -0.695), (-0.3, 0.296), (0.199, 1.163), (-0.765, 0.892), (1.163, 0.892),
    ], 'Wikimedia Commons Adamantane_acsv.svg'),
    'twistane': ('C1CC2C3CCC(C1C3)C2', [
        (0.254, 1.378), (-0.236, 1.378), (-0.506, 0.386), (0.513, -0.395), (0.271, -1.397),
        (-0.258, -1.397), (-0.522, -0.395), (0.513, 0.381), (1.494, 0.031), (-1.524, 0.031),
    ], 'Wikimedia Commons Twistane1.svg'),
    'bicyclo[3.2.1]octane': ('C1CC2CCC(C1)C2', [
        (1.267, -0.538), (0.648, -0.007), (-0.574, -0.007), (-1.223, -0.891), (-0.457, -0.435),
        (-0.059, 0.449), (1.12, 0.479), (-0.722, 0.95),
    ], 'Wikimedia Commons Tropanol_-_tropine.svg'),
    'bicyclo[2.2.2]octane': ('C1CC2CCC1CC2', [
        (-0.22, 0.707), (0.312, 1.419), (0.312, 0.337), (1.1, -0.09), (0.616, -0.907),
        (-0.22, -0.546), (-1.207, -0.831), (-0.694, -0.09),
    ], 'Wikimedia Commons Barrelene_structure.png'),
}


SHAPES = {name: (text, _seen(xyz, views[0])) for name, (text, views, xyz) in MODELS.items()}
SHAPES.update((name, (text, pts)) for name, (text, pts, _) in DRAWN.items() if len(pts) == len(smiles(text)))


def skeleton(mol, atoms):
    """`(skeleton, order)`: `atoms` of `mol` and their bonds as an all-carbon, all-single-bond molecule, `order[i]`
    the atom of `mol` the skeleton's atom `i + 1` stands for."""
    order = sorted(atoms)
    index = {}
    out = MoleculeContainer()
    with out.edit() as e:
        for x in order:
            index[x] = e.add_atom('C')
        for x in order:
            for y in mol.neighbors_of(x):
                if y in index and x < y:
                    e.add_bond(index[x], index[y], 1)
    return out, order


def _ring(mol):
    """The atoms of `mol`'s largest ring system, in its atom order."""
    system = max(({x for r in rs for x in r} for _, rs in _blocks(mol)), key=len)
    return [x for x in mol if x in system]


def _blocks(mol):
    from .rings import ring_blocks
    return ring_blocks(mol)


@lru_cache()
def _shapes():
    """Per skeleton, its drawings: `(name, skeleton, order, points, inside, bearing)`, points and the sets keyed by
    skeleton atom, `bearing` the atoms the reference molecule substitutes; drawn shapes first, in table order."""
    out = {}
    tables = [(name, text, [pts], True) for name, (text, pts, _) in DRAWN.items()]
    tables += [(name, text, [_seen(xyz, v) for v in views], False) for name, (text, views, xyz) in MODELS.items()]
    drawn = set()
    for name, text, seen, is_drawn in tables:
        mol = smiles(text)
        ring = _ring(mol)
        skel, order = skeleton(mol, set(ring))
        key = skel.canonical_bytes
        if not is_drawn and key in drawn:                # a drawn shape replaces the model's views
            continue
        if is_drawn:
            drawn.add(key)
        index = {x: i + 1 for i, x in enumerate(order)}
        bearing = {index[x] for x in ring if any(y not in index for y in mol.neighbors_of(x))}
        for p in seen:
            at = dict(zip(ring, p))
            inner = {index[ring[i]] for i in _inside(p)}
            out.setdefault(key, []).append((name, skel, {index[x]: at[x] for x in ring}, inner, bearing))
    return out


def shape(mol, atoms, rank=None):
    """The predefined shape of the ring system `atoms` and its mirror image, `[{atom: (x, y)}]` in plane coordinates,
    or None when no shape matches: the drawing and atom mapping whose reference substitutes most of the atoms `mol`
    substitutes, then leave fewest of them inside the outline, in table order, then by the canonical ranks the
    mapping gives the drawing's atoms, so the choice ignores atom order."""
    rank = mol.canonical_order() if rank is None else rank
    skel, order = skeleton(mol, atoms)
    found = _shapes().get(skel.canonical_bytes)
    if not found:
        return None
    mine = {i + 1 for i, x in enumerate(order) if any(y not in atoms for y in mol.neighbors_of(x))}
    best = None
    for k, (name, ref, pts, inner, bearing) in enumerate(found):
        for mapping in islice(ref.as_query().get_mapping(skel), MAPPINGS):
            subs = {i for i, x in mapping.items() if x in mine}
            score = (-len(subs & bearing), len(subs & inner), k,
                     tuple(rank[order[mapping[i] - 1]] for i in sorted(mapping)))
            if best is None or score < best[0]:
                best = (score, pts, mapping)
    if best is None:
        return None
    _, pts, mapping = best
    plane = {order[mapping[i] - 1]: (p[0] * L, p[1] * L) for i, p in pts.items()}
    return [plane, {x: (-p[0], p[1]) for x, p in plane.items()}]     # mirrored left to right: up stays up
