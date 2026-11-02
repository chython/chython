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
| choice | per molecule, the drawing and symmetry leaving fewest of the molecule's wedged stereocentres within \
`CROSS` of a bond crossing, where a wedge cannot say which bond is in front; then whose reference substitutes \
most of the atoms the molecule substitutes, then leaving fewest of them inside its outline, drawn shapes first |
| drawn | `DRAWN`: a reference molecule and its ring atoms' points read from a drawing of it, mean bond 1; several \
may share a skeleton, and they replace the model's views; a drawing in `BESIDE` draws only a cage with a ring fused \
to it, the bare skeleton keeping its model's views |
| shape | `SHAPES`: points in `L`: the drawn skeletons, each other model seen along its first view |
| depth | each drawing's atoms' z, toward the viewer: a model's view its own, a drawn skeleton's from its model \
aligned to the drawing by orthographic Procrustes over the skeleton's symmetries and the model's mirror image; read \
the way round that puts a one-atom bridge in front |
| drawing | `shape(mol, atoms)`: the matched shape and its mirror image in plane coordinates, or None |
| fused | `fused3d(mol, atoms, rings)`: a ring system of a shape and one more ring on one of its edges, a benzo- or \
pyrrolidino-fused cage: all of it in perspective, the ring a regular polygon in 3D \
seen as the cage is, the view the one whose shortest projected bond is longest; any other fused system is flat |

| shape | skeletons it draws |
| --- | --- |
| norbornane | camphor, norbornene |
| bicyclo[2.2.2]octane | quinuclidine, DABCO |
| bicyclo[3.2.1]octane | tropane |
| bicyclo[3.3.1]nonane | granatane; fused to a ring, pseudopelletierine |
| galanthamine | the benzofuro-benzazepine of galanthamine, lycoramine |
| bicyclo[1.1.1]pentane, bicyclo[2.1.1]hexane, adamantane, twistane, cubane, prismane, tetrahedrane | |
"""
from functools import lru_cache
from itertools import islice
from math import atan2, cos, pi, sin, tan

from ...core import MoleculeContainer, smiles
from .lattice import L


__all__ = ['DRAWN', 'MODELS', 'SHAPES', 'skeleton', 'shape', 'shape3d', 'fused3d']

EDGE = .08               # an atom this near its drawing's outline, in bonds, is on it
CROSS = .55              # an atom this near a crossing of two of its drawing's bonds, in bonds, is at it
SAMPLES = 96             # view directions a fused cage of a model skeleton is tried along
BESIDE = {'pseudopelletierine'}    # drawings kept for fused cages: the bare skeleton keeps its model's views
CLEAR = .6 * L           # a fused cage's view leaves every bond and unbonded pair at least this long
MAPPINGS = 48            # skeleton symmetries tried for a molecule's substituents
HANDED = .2              # a bridgehead's signed volume, in bonds cubed, once turned to its model's handedness


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
    'galanthamine': ('CN1CC[C@]23C=C[C@@H](O)C[C@@H]2Oc2c3c(C1)ccc2OC', [
        (2.143, 0.346), (1.243, -0.546), (0.512, 0.054), (-0.285, -0.509), (0.549, -0.97), (0.549, -1.958),
        (-0.285, -2.433), (-1.126, -1.958), (-1.126, -0.97), (-1.77, -0.019), (-1.126, 0.953), (-0.285, 0.471),
        (0.549, 0.953), (1.324, 0.273), (0.549, 1.941), (-0.285, 2.431), (-1.126, 1.941),
    ], 'Wikimedia Commons Galanthamine.svg'),
    'pseudopelletierine': ('CN1C2CCCC1CC(=O)C2', [
        (-0.031, 1.361), (-0.217, -0.103), (-1.308, -0.275), (-1.132, -0.896), (-0.727, 0.282),
        (0.385, 0.445), (1.175, 0.209), (1.175, -0.757), (0.681, -0.265),
    ], 'Wikimedia Commons Pseudopelletierine.svg'),
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
    """Per skeleton, its drawings: `(name, skeleton, points, inside, bearing, crossed, depth)`, points, depth and the
    sets keyed by skeleton atom, `bearing` the atoms the reference molecule substitutes, `crossed` those `_crossed`,
    `depth` each atom's z in bonds; drawn shapes first, in table order."""
    out = {}
    models = {}                                         # skeleton key -> (skeleton, {atom: xyz})
    for text, _, xyz in MODELS.values():
        mol = smiles(text)
        ring = _ring(mol)
        skel, order = skeleton(mol, set(ring))
        index = {x: i + 1 for i, x in enumerate(order)}
        models.setdefault(skel.canonical_bytes, (skel, {index[x]: p for x, p in zip(ring, xyz)}))
    tables = [(name, text, [pts], None, True) for name, (text, pts, _) in DRAWN.items()]
    tables += [(name, text, [_seen(xyz, v) for v in views], [_toward(xyz, v) for v in views], False)
               for name, (text, views, xyz) in MODELS.items()]
    drawn, beside = set(), set()
    for name, text, seen, depths, is_drawn in tables:
        mol = smiles(text)
        ring = _ring(mol)
        skel, order = skeleton(mol, set(ring))
        key = skel.canonical_bytes
        if not is_drawn and key in drawn and key not in beside:     # a drawn shape replaces the model's views
            continue
        if is_drawn:
            drawn.add(key)
            if name in BESIDE:
                beside.add(key)
        index = {x: i + 1 for i, x in enumerate(order)}
        bearing = {index[x] for x in ring if any(y not in index for y in mol.neighbors_of(x))}
        for k, p in enumerate(seen):
            at = dict(zip(ring, p))
            inner = {index[ring[i]] for i in _inside(p)}
            pts = {index[x]: at[x] for x in ring}
            if depths is None:
                depth = _aligned(skel, pts, *models[key]) if key in models else dict.fromkeys(pts, 0.)
            else:
                depth = {index[x]: z for x, z in zip(ring, depths[k])}
            out.setdefault(key, []).append((name, skel, pts, inner, bearing, _crossed(skel, pts),
                                            _bridge_front(skel, depth)))
    return out


def _bridge_front(skel, depth):
    """`depth`, or its negation, so that the one-atom bridges of `skel` lie toward the viewer: a projection reads
    either way round, and a drawing of a bridged bicycle shows its bridge passing in front, camphor's and
    7-oxanorbornane's alike."""
    heads = {x for x in depth if sum(1 for y in skel.neighbors_of(x)) == 3}
    bridges = [x for x in depth if x not in heads and len(skel.neighbors_of(x)) == 2 and
               all(y in heads for y in skel.neighbors_of(x))]
    if not bridges:
        return depth
    mean = sum(depth.values()) / len(depth)
    if sum(depth[x] for x in bridges) / len(bridges) < mean:
        return {x: -z for x, z in depth.items()}
    return depth


def _toward(xyz, view):
    """Each point of `xyz`'s z as `_seen` sees it along `view`: toward the viewer."""
    n = sum(v * v for v in view) ** .5
    return [sum(a * b / n for a, b in zip(p, view)) for p in xyz]


def _aligned(skel, pts, model, xyz):
    """`skel`'s atoms' z under the drawing `pts`: the model's points `xyz` of skeleton `model`, or their mirror image,
    turned and scaled so their orthographic projection comes nearest `pts`, over every mapping of `model` onto
    `skel`."""
    best = None
    for mirror in (1., -1.):
        for mapping in islice(model.as_query().get_mapping(skel), MAPPINGS):
            pairs = [((xyz[i][0], xyz[i][1], mirror * xyz[i][2]), pts[j]) for i, j in mapping.items()]
            fit = _procrustes([a for a, _ in pairs], [b for _, b in pairs])
            if fit is not None and (best is None or fit[0] < best[0]):
                best = (fit[0], {j: z for (_, j), z in zip(mapping.items(), fit[1])},
                        {j: a for a, (_, j) in zip((a for a, _ in pairs), mapping.items())})
    return _handed(skel, pts, best[1], best[2]) if best else dict.fromkeys(pts, 0.)


def _handed(skel, pts, z, model):
    """The depth `z` under the drawing `pts`, each bridgehead turned to the handedness the clearest one shares with
    the `model`'s points: a bridgehead drawn almost flat would otherwise read against its partner.  A turned one goes
    `HANDED` past flat, its neighbours kept."""
    def volume(p, a):
        (a1, a2, a3), (b1, b2, b3), (c1, c2, c3) = [[p[n][k] - p[a][k] for k in range(3)]
                                                    for n in sorted(skel.neighbors_of(a))]
        return a1 * (b2 * c3 - b3 * c2) - a2 * (b1 * c3 - b3 * c1) + a3 * (b1 * c2 - b2 * c1)
    heads = [x for x in pts if len(skel.neighbors_of(x)) == 3]
    if not heads:
        return z
    drawn = {x: (pts[x][0], pts[x][1], z[x]) for x in pts}
    ratio = {x: volume(drawn, x) / volume(model, x) for x in heads}
    clearest = max(heads, key=lambda x: abs(ratio[x]))
    want = 1 if ratio[clearest] > 0 else -1
    z = dict(z)
    for x in heads:
        if ratio[x] * want > 0:
            continue
        v0 = volume(drawn, x)
        drawn[x] = (*pts[x], z[x] + 1)
        slope = volume(drawn, x) - v0                   # the volume is linear in the bridgehead's own depth
        if abs(slope) < 1e-9:
            continue
        target = want * (1 if volume(model, x) > 0 else -1) * HANDED
        z[x] += (target - v0) / slope
        drawn[x] = (*pts[x], z[x])
    return z


def _procrustes(xs, ys):
    """`(residual, depths)`: the rotation and scale whose orthographic projection takes the 3D points `xs` nearest the
    2D points `ys`, and each point's z under it; None for a degenerate cloud."""
    n = len(xs)
    mx = [sum(p[k] for p in xs) / n for k in range(3)]
    my = [sum(p[k] for p in ys) / n for k in range(2)]
    xc = [[p[k] - mx[k] for k in range(3)] for p in xs]
    yc = [[p[k] - my[k] for k in range(2)] for p in ys]
    b = [[sum(y[r] * x[c] for x, y in zip(xc, yc)) for c in range(3)] for r in range(2)]       # Y^T X, 2 x 3
    sa, sb, sd = (sum(v * v for v in b[0]), sum(u * v for u, v in zip(*b)), sum(v * v for v in b[1]))
    det = sa * sd - sb * sb
    if det <= 1e-12:
        return None
    t = det ** .5
    k = (sa + sd + 2 * t) ** .5
    ra, rb, rd = (sa + t) / k, sb / k, (sd + t) / k                         # sqrt(B B^T)
    rdet = ra * rd - rb * rb
    ia, ib, id_ = rd / rdet, -rb / rdet, ra / rdet                          # its inverse
    rows = [[ia * b[0][c] + ib * b[1][c] for c in range(3)], [ib * b[0][c] + id_ * b[1][c] for c in range(3)]]
    r1, r2 = rows
    r3 = [r1[1] * r2[2] - r1[2] * r2[1], r1[2] * r2[0] - r1[0] * r2[2], r1[0] * r2[1] - r1[1] * r2[0]]
    scale = (ra + rd) / sum(v * v for x in xc for v in x)
    residual = sum((y[0] - scale * sum(a * c for a, c in zip(x, r1))) ** 2 +
                   (y[1] - scale * sum(a * c for a, c in zip(x, r2))) ** 2 for x, y in zip(xc, yc))
    return residual, [scale * sum(a * c for a, c in zip(x, r3)) for x in xc]


def _crossed(skel, pts):
    """The atoms of `skel` nearer than `CROSS` to a point where two of its bonds cross in the drawing `pts`."""
    bonds = [(b.n, b.m) for b in skel.bonds()]
    points = []
    for i, (a, b) in enumerate(bonds):
        for c, d in bonds[i + 1:]:
            if {a, b} & {c, d}:
                continue
            (x1, y1), (x2, y2), (x3, y3), (x4, y4) = pts[a], pts[b], pts[c], pts[d]
            den = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)
            if abs(den) < 1e-9:
                continue
            t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / den
            u = ((x1 - x3) * (y1 - y2) - (y1 - y3) * (x1 - x2)) / den
            if 0 < t < 1 and 0 < u < 1:
                points.append((x1 + t * (x2 - x1), y1 + t * (y2 - y1)))
    return {x for x, (px, py) in pts.items() if any(((px - cx) ** 2 + (py - cy) ** 2) ** .5 < CROSS
                                                    for cx, cy in points)}


def _wedged(mol):
    """The atoms of `mol` a drawing wedges: tetrahedral and allene stereocentres of a stated parity."""
    from ...core import SU_ALLENE, SU_TETRA
    return {u['anchor'] for u in mol.stereo_units()
            if u['kind'] in (SU_TETRA, SU_ALLENE) and mol.parity_of(u['anchor'])}


def shape(mol, atoms, rank=None):
    """The predefined shape of the ring system `atoms` and its mirror image, `[{atom: (x, y)}]` in plane coordinates,
    or None when no shape matches: the drawing and atom mapping leaving fewest wedged stereocentres of `mol` at a
    crossing, then whose reference substitutes most of the atoms `mol` substitutes, then leave fewest of them inside
    the outline, in table order, then by the canonical ranks the mapping gives the drawing's atoms, so the choice
    ignores atom order."""
    found = shape3d(mol, atoms, rank)
    return None if found is None else [p for p, _ in found]


def shape3d(mol, atoms, rank=None):
    """`shape`'s two drawings, each with its atoms' depth: `[({atom: (x, y)}, {atom: z})]`, z toward the viewer in
    plane units; the mirror image, turned left to right, keeps the depth."""
    ranked = _ranked(mol, atoms, rank)
    if not ranked:
        return None
    _, plane, z = ranked[0]
    return [(plane, z), ({x: (-p[0], p[1]) for x, p in plane.items()}, z)]     # mirrored left to right: up stays up


def _ranked(mol, atoms, rank=None, fused=False):
    """`[(score, {atom: (x, y)}, {atom: z})]` of every drawing and mapping of the shape of `atoms`, best first: of a
    `fused` cage its `BESIDE` drawing alone when it has one, of a bare one every other drawing."""
    rank = mol.canonical_order() if rank is None else rank
    skel, order = skeleton(mol, atoms)
    found = _shapes().get(skel.canonical_bytes)
    if not found:
        return []
    own = [f for f in found if f[0] in BESIDE]
    found = (own or found) if fused else [f for f in found if f[0] not in BESIDE]
    mine = {i + 1 for i, x in enumerate(order) if any(y not in atoms for y in mol.neighbors_of(x))}
    wedged = {i + 1 for i, x in enumerate(order) if x in _wedged(mol)}
    out = []
    for k, (name, ref, pts, inner, bearing, crossed, depth) in enumerate(found):
        for mapping in islice(ref.as_query().get_mapping(skel), MAPPINGS):
            subs = {i for i, x in mapping.items() if x in mine}
            score = (sum(mapping[i] in wedged for i in crossed), -len(subs & bearing), len(subs & inner), k,
                     tuple(rank[order[mapping[i] - 1]] for i in sorted(mapping)))
            plane = {order[mapping[i] - 1]: (p[0] * L, p[1] * L) for i, p in pts.items()}
            out.append((score, plane, {order[mapping[i] - 1]: depth[i] * L for i in pts}))
    out.sort(key=lambda t: t[0])
    return out


def fused3d(mol, atoms, rings, rank=None):
    """`(core, [({atom: (x, y)}, {atom: z})])`: the ring system `atoms` of `rings` as a shape and its mirror image,
    `core` the atoms the shape draws and the only ones with a depth: the system itself when a shape draws it whole,
    or all its rings but one fused on one edge.  That ring is a regular
    polygon in 3D, in the plane its fusion atoms' cage bonds make with their edge, seen as the cage is: of the
    cage's drawings and `SAMPLES` views of its model, the one leaving fewest wedged stereocentres at a crossing, then
    whose shortest projected bond and nearest two unbonded atoms are `CLEAR`, then fewest one-atom bridges inside
    its outline, a camphor's bridge standing out of it, then the clearest.  None otherwise."""
    whole = shape3d(mol, atoms, rank)
    if whole is not None:
        return frozenset(atoms), whole
    if len(rings) < 3:
        return None
    for i, ring in enumerate(rings):
        core = frozenset().union(*(r for k, r in enumerate(rings) if k != i))
        if len(core & set(ring)) != 2:
            continue
        bonds = [(b.n, b.m) for b in mol.bonds() if b.n in atoms and b.m in atoms]
        heads = {x for x in core if sum(y in core for y in mol.neighbors_of(x)) == 3}
        bridges = [x for x in core if x not in heads and sum(y in heads for y in mol.neighbors_of(x)) == 2]
        best = None
        for score, plane, z in _ranked(mol, core, rank, True) + _sampled(mol, core, rank):
            added = _ring3d(mol, plane, z, ring)
            if added is None:
                continue
            full = {**plane, **added}
            shortest = min(((full[a][0] - full[b][0]) ** 2 + (full[a][1] - full[b][1]) ** 2) ** .5 for a, b in bonds)
            pts = sorted(full)
            closest = min((((full[a][0] - full[b][0]) ** 2 + (full[a][1] - full[b][1]) ** 2) ** .5
                           for k, a in enumerate(pts) for b in pts[k + 1:] if b not in mol.neighbors_of(a)),
                          default=L)
            inner = {pts[k] for k in _inside([full[x] for x in pts])}
            clear = min(shortest, closest)
            key = (score[0], clear < CLEAR, sum(x in inner for x in bridges), -round(clear / L, 2), score[1:])
            if best is None or key < best[0]:
                best = key, full, z
        if best is not None:
            _, full, z = best
            return core, [(full, z), ({x: (-p[0], p[1]) for x, p in full.items()}, z)]
    return None


def _sampled(mol, atoms, rank=None):
    """`[(score, {atom: (x, y)}, {atom: z})]` of the model of the shape of `atoms` seen along `SAMPLES` directions
    spread over the sphere, each mapping onto `atoms`; a drawn shape, whose views are its drawing, offers none."""
    from math import acos, sqrt
    rank = mol.canonical_order() if rank is None else rank
    skel, order = skeleton(mol, atoms)
    key = skel.canonical_bytes
    if key not in _models() or any(not is_model for is_model in _drawn_keys().get(key, ())):
        return []
    model, xyz = _models()[key]
    wedged = {i + 1 for i, x in enumerate(order) if x in _wedged(mol)}
    out = []
    golden = pi * (3 - sqrt(5))
    for k in range(SAMPLES):
        h = 1 - 2 * (k + .5) / SAMPLES
        view = (sqrt(1 - h * h) * cos(golden * k), sqrt(1 - h * h) * sin(golden * k), h)
        ids = sorted(xyz)
        seen = dict(zip(ids, _seen([xyz[i] for i in ids], view)))
        depth = dict(zip(ids, _toward([xyz[i] for i in ids], view)))
        crossed = _crossed(model, seen)
        depth = _bridge_front(model, depth)
        for mapping in islice(model.as_query().get_mapping(skel), MAPPINGS):
            score = (sum(mapping[i] in wedged for i in crossed), 0, 0, len(_shapes()[key]),
                     tuple(rank[order[mapping[i] - 1]] for i in sorted(mapping)))
            out.append((score, {order[mapping[i] - 1]: (p[0] * L, p[1] * L) for i, p in seen.items()},
                        {order[mapping[i] - 1]: depth[i] * L for i in seen}))
    return out


@lru_cache()
def _models():
    """Skeleton key -> `(skeleton, {atom: xyz})` of each model."""
    out = {}
    for text, _, xyz in MODELS.values():
        mol = smiles(text)
        ring = _ring(mol)
        skel, order = skeleton(mol, set(ring))
        index = {x: i + 1 for i, x in enumerate(order)}
        out.setdefault(skel.canonical_bytes, (skel, {index[x]: p for x, p in zip(ring, xyz)}))
    return out


@lru_cache()
def _drawn_keys():
    """Skeleton key -> one `False` per drawn shape of it: a drawn skeleton keeps to its drawing."""
    out = {}
    for text, _, _ in DRAWN.values():
        mol = smiles(text)
        out.setdefault(skeleton(mol, set(_ring(mol)))[0].canonical_bytes, []).append(False)
    return out


def _ring3d(mol, plane, z, ring):
    """`{atom: (x, y)}` of the atoms `ring` adds to the cage `plane` at depth `z`: the regular polygon on the edge it
    shares, in 3D in the plane of that edge and the bonds leaving it into the cage, outward, projected as the cage is;
    None when the ring shares no single edge or an atom lands within half a bond of the cage."""
    n = len(ring)
    shared = [k for k, x in enumerate(ring) if x in plane]
    if len(shared) != 2 or (shared[1] - shared[0]) % n not in (1, n - 1):
        return None
    i, j = shared
    if (j - i) % n != 1:
        i, j = j, i                                     # ring[i] then ring[j] in the ring's order
    a, b = ring[i], ring[j]

    def point(x):
        return plane[x][0], plane[x][1], z[x]
    pa, pb = point(a), point(b)
    side = sum((u - v) ** 2 for u, v in zip(pa, pb)) ** .5
    e = tuple((v - u) / side for u, v in zip(pa, pb))
    out = [0., 0., 0.]
    for x, px in ((a, pa), (b, pb)):
        for y in mol.neighbors_of(x):
            if y in plane and y not in ring:
                for k in range(3):
                    out[k] += px[k] - point(y)[k]
    dot = sum(u * v for u, v in zip(out, e))
    w = [u - dot * v for u, v in zip(out, e)]
    norm = sum(u * u for u in w) ** .5
    if norm < 1e-9:
        return None
    w = [u / norm for u in w]
    apothem = side / (2 * tan(pi / n))
    r = side / (2 * sin(pi / n))
    centre = [(u + v) / 2 + apothem * c for u, v, c in zip(pa, pb, w)]
    tb = atan2(-apothem, side / 2)                      # b's angle about the centre, in the (e, w) frame
    added = {}
    for k in range(1, n - 1):
        t = tb + 2 * pi * k / n
        p3 = [c + r * (cos(t) * u + sin(t) * v) for c, u, v in zip(centre, e, w)]
        added[ring[(j + k) % n]] = (p3[0], p3[1])
    if any(((p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2) ** .5 < L / 2 for p in added.values() for q in plane.values()):
        return None
    return added


def _polygons(mol, plane, rest):
    """`plane` with each ring of `rest` added as a regular polygon on the one edge it shares with `plane`, on the side
    away from the drawing; None when a ring shares no single edge, or an atom it adds lands within half a bond of
    one drawn."""
    rest = list(rest)
    while rest:
        for ring in rest:
            shared = [i for i, x in enumerate(ring) if x in plane]
            if len(shared) == 2 and (shared[1] - shared[0]) % len(ring) in (1, len(ring) - 1):
                break
        else:
            return None
        rest.remove(ring)
        n = len(ring)
        i, j = shared
        if (j - i) % n != 1:
            i, j = j, i                                 # ring[i] then ring[j] in the ring's order
        (ax, ay), (bx, by) = plane[ring[i]], plane[ring[j]]
        side = ((ax - bx) ** 2 + (ay - by) ** 2) ** .5
        mx, my = (ax + bx) / 2, (ay + by) / 2
        nx, ny = (by - ay) / side, (ax - bx) / side     # a normal of the edge
        cx = sum(p[0] for p in plane.values()) / len(plane)
        cy = sum(p[1] for p in plane.values()) / len(plane)
        if (mx - cx) * nx + (my - cy) * ny < 0:
            nx, ny = -nx, -ny
        apothem = side / (2 * tan(pi / n))
        ox, oy = mx + nx * apothem, my + ny * apothem
        # from ring[j] onward round the centre, the turn's sense the one taking ring[i] to ring[j]
        t0 = atan2(by - oy, bx - ox)
        ti = atan2(ay - oy, ax - ox)
        step = 2 * pi / n if (t0 - ti) % (2 * pi) < pi else -2 * pi / n
        r = side / (2 * sin(pi / n))
        added = {}
        for k in range(1, n - 1):
            x = ring[(j + k) % n]
            added[x] = (ox + r * cos(t0 + k * step), oy + r * sin(t0 + k * step))
        if any(((p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2) ** .5 < L / 2 for p in added.values() for q in plane.values()):
            return None
        plane.update(added)
    return plane
