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
"""Authored drawings of ring blocks the lattice cannot embed without a conflict.

A template is a skeleton, matched against a block's graph alone: elements, charges and bond orders are ignored,
and the block must be the whole skeleton, no atom more or less. Coordinates are in bonds of length `L`, one per
skeleton atom in SMILES order, and drawn at that scale. Of several embeddings the one whose atoms, taken in skeleton
order, rank lowest is drawn. A handed template keeps its authored hand where the pose finds a mirror as wide.

| skeleton | drawing |
| --- | --- |
| morphinan | phenanthrene rings on the honeycomb; the ethylamine bridge from C13 to C9 crosses C8-C14 once |
| 4,5-epoxymorphinan | morphinan, the ether oxygen in the bay, both its bonds stretched |
| taxane | rings A and C regular, B an arc of 1.1 L chords; turned 30 degrees: C15's methyls fan into B on lattice |
| artemisinane | cyclohexane and lactone regular on the honeycomb; the peroxide bent inside the cycloheptane |
| porphine | a square: regular pyrroles on the diagonals, nitrogens inward, meso carbons on the axes at 120 degrees |
| phorbine | porphine; the isocyclic ring's two atoms where a regular pentagon on the beta-meso chord puts them |
| twistane | a hexagon of 1.5 L sides; the C2 bridge between opposite vertices bows, the other runs beside a side |
| crown | a crown-type ring: each repeating unit a hexagon arc, the heteroatoms the inward corners; `crown(n, period)` |
| strychnidine | handed: benzene on the left; benzene, lactam and cyclohexane regular on the honeycomb, the indoline a |
| | hexagon short of its top, the pyrrolidine above, the N-C4-C5 bridge down the right, the oxepane below it |
"""
from math import asin, atan2, cos, pi, radians, sin, sqrt
from ...core import smiles


__all__ = ['template_block', 'crown']

H = 3 ** .5 / 2
# morphine numbering: C1 C2 C3 C4 O C5 C6 C7 C8 C14 C9 N17 C16 C15 C13 C12 C11 C10
_EPOXY = ((-2 * H, -1.), (-2 * H, -2.), (-H, -2.5), (0., -2.), (H + .25, -1.5 - H / 2), (2 * H, -1.),
          (3 * H, -.5), (3 * H, .5), (2 * H, 1.), (H, .5), (0., 1.), (.25, 2.), (1., 1.25), (1.5, .25), (H, -.5),
          (0., -1.), (-H, -.5), (-H, .5))

# artemisinin numbering: C6 C7 C8 C8a C9 C10 O11 C12 C12a C5a C5 C4 C3 O13 O2 O1; hexagons fixed, the rest by
# least squares on unit bonds, regular ring angles and a 1.2 bond floor between atoms not bonded
_ARTEMISININ = ((-2 * H, 0.), (-2 * H, -1.), (-H, -1.5), (0., -1.), (H, -1.5), (2 * H, -1.), (2 * H, 0.), (H, .5),
                (0., 0.), (-H, .5), (-.988, 1.583), (-.469, 2.54), (.612, 2.473), (1.206, 1.542), (.16, 1.726),
                (-.023, .924))

# twistane numbering: C1 C2 C3 C4 C5 C6 C7 C8 C9 C10 of C1CC2CC3C1CC2CC3; hexagon C1 C2 C3 C8 C7 C6, bridges
# C6-C5-C4-C3 and C5-C10-C9-C8 inside, by least squares on unit inner bonds and a 0.8 bond floor between atoms
_TWISTANE = ((-3 * H / 2, .75), (-3 * H / 2, -.75), (0., -1.5), (-.4, -.48), (-.4, .56), (0., 1.5),
             (3 * H / 2, .75), (3 * H / 2, -.75), (.56, -.46), (.56, .42))

# strychnine numbering: C1 C2 N3 C4 C5 C6 C7 O8 C9 C10 C11 N13 C14 C15 C16 C17 C18 C19 C20 C21 C22 C23 C24 C25
_STRYCHNIDINE = ((1.73, 1.08), (2.09, 2.08), (3.11, 2.08), (5.3, 1.08), (5.3, -.76), (5.61, -1.68), (5.2, -2.57),
                 (4.33, -2.74), (4 * H, -2.), (3 * H, -2.5), (2 * H, -2.), (2 * H, -1.), (3 * H, -.5), (4 * H, -1.),
                 (5 * H, -.5), (5 * H, .5), (4 * H, 1.), (3 * H, .5), (H, .5), (0., 1.), (-H, .5), (-H, -.5),
                 (0., -1.), (H, -.5))


def _taxane():
    """Taxane in the atom order of its skeleton, taxol numbering C1-C14, C15, C20, O: ring A a regular hexagon with
    C15 on top; ring B's bonds from C1 to C11 chords of the circle whose remaining chord is A's C1-C11, 1.1 long but
    the unit C3-C8, room for C15's methyls in B; ring C a regular hexagon outside B on C3-C8, the oxetane a square
    outside C. Turned 30 degrees clockwise, so the lattice fans C15's methyls 90, 60, 90 degrees from its ring
    bonds."""
    chords = (1.1, 1.1, 1., 1.1, 1.1, 1.1)                  # C1-C2 C2-C3 C3-C8 C8-C9 C9-C10 C10-C11
    lo, hi = .55, 3.
    for _ in range(60):                                     # the radius closing the arc on a chord of sqrt 3
        r = (lo + hi) / 2
        g = pi - sum(asin(c / 2 / r) for c in chords)       # half the angle the C1-C11 chord subtends
        if g > 0 and 2 * r * sin(g) < 2 * H:
            lo = r
        else:
            hi = r
    oy = .5 + r * cos(g)
    a = -pi / 2 + g
    arc = [(r * cos(a), oy + r * sin(a))]
    for c in chords:                                        # C1 .. C11, counterclockwise
        a += 2 * asin(c / 2 / r)
        arc.append((r * cos(a), oy + r * sin(a)))

    def fuse(p, q, n, o):
        """Unit regular n-gon on edge p-q away from point o: its atoms from q's neighbour to p's."""
        ux, uy = q[0] - p[0], q[1] - p[1]
        left = (o[0] - p[0]) * -uy + (o[1] - p[1]) * ux < 0
        pts, x, y, a = [], q[0], q[1], atan2(uy, ux)
        for _ in range(n - 2):
            a += (1 if left else -1) * 2 * pi / n
            x, y = x + cos(a), y + sin(a)
            pts.append((x, y))
        return pts

    c1, c2, c3, c8, c9, c10, c11 = arc
    c7, c6, c5, c4 = fuse(c3, c8, 6, (0., oy))
    o, c20 = fuse(c4, c5, 4, c3)
    xy = (c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, (-H, -.5), (0., -1.), (H, -.5), (0., 1.), c20, o)
    k, j = cos(pi / 6), sin(pi / 6)
    return tuple((x * k + y * j, y * k - x * j) for x, y in xy)


def _porphine():
    """Porphine in the atom order of its skeleton: pyrrole k centred on the diagonal at 45 + 90 k degrees, its
    alpha carbons `r sin 72` off the axes on either side, so an alpha-meso-alpha angle is 120 degrees."""
    r = 1 / (2 * sin(pi / 5))                               # circumradius of the unit pentagon
    c = (sqrt(3) / 2 + r * sin(.35 * pi)) * sqrt(2)         # centre distance putting alpha at y = cos 30
    m = c * sqrt(.5) + r * cos(.35 * pi) + .5               # meso distance: alpha x plus a bond at -60 degrees
    rings, meso = [], []
    for k in range(4):
        t = pi / 4 + k * pi / 2
        o = (c * cos(t), c * sin(t))
        rings.append([(o[0] + r * cos(t + a), o[1] + r * sin(t + a))
                      for a in (-.6 * pi, -.2 * pi, .2 * pi, .6 * pi, pi)])     # alpha, beta, beta, alpha, N
        meso.append((m * cos(t + pi / 4), m * sin(t + pi / 4)))
    a, b, c_, d = rings
    return (a[1], a[2], a[3], meso[0], *b, meso[1], *c_, meso[2], *d, meso[3], a[0], a[4])


def _phorbine():
    """Porphine and the isocyclic ring's two atoms on the chord from pyrrole C's second beta carbon to its meso."""
    xy = _porphine()
    (px, py), (qx, qy), (ax, ay) = xy[12], xy[15], xy[13]
    ux, uy = qx - px, qy - py
    d = sqrt(ux * ux + uy * uy)
    ux, uy = ux / d, uy / d
    nx, ny = -uy, ux
    if (ax - px) * nx + (ay - py) * ny > 0:                 # away from the alpha carbon
        nx, ny = -nx, -ny
    r = d / (2 * sin(.4 * pi))
    mx, my = (px + qx) / 2, (py + qy) / 2
    h, w = r * (sin(.1 * pi) + sin(.3 * pi)), r * cos(.3 * pi)
    return (*xy, (mx - w * ux + h * nx, my - w * uy + h * ny), (mx + w * ux + h * nx, my + w * uy + h * ny))


TEMPLATES = (
    ('morphinan', 'C1CCCC2C1CC3NCCC24CCCCC34', tuple(_EPOXY[i] for i in (0, 1, 2, 3, 15, 16, 17, 10, 11, 12, 13, 14,
                                                                         5, 6, 7, 8, 9)), False),
    ('4,5-epoxymorphinan', 'C1CCC2OC3CCCC4C5NCCC34C2C1C5', _EPOXY, False),
    ('taxane', 'C12CC3C4C5CCC3CCC6CCC1.C26.C4O5', _taxane(), False),
    ('artemisinane', 'C1CCC2C(COC3C24C1CCC(O3)OO4)', _ARTEMISININ, False),
    ('porphine', 'C1=CC2=CC3=CC=C(N3)C=C4C=CC(=N4)C=C5C=CC(=N5)C=C1N2', _porphine(), False),
    ('phorbine', 'C1=CC2=CC3=CC=C(N3)C=C4C=C6C(=N4)C7=C5C=CC(=N5)C=C1N2.C6C7', _phorbine(), False),
    ('twistane', 'C1CC2CC3C1CC2CC3', _TWISTANE, False),
    ('strychnidine', 'C1CN2CC3=CCOC4CCN5C6C4C3CC2C61C7=CC=CC=C75', _STRYCHNIDINE, True),
)
_CACHE = {}


def _skeletons():
    """`(adjacency, coordinates, handed)` per template, by atom count; parsed on first use."""
    if not _CACHE:
        for name, s, xy, hand in TEMPLATES:
            mol = smiles(s)
            order = list(mol)
            idx = {x: i for i, x in enumerate(order)}
            adj = [sorted(idx[y] for y in mol.neighbors_of(x)) for x in order]
            assert len(adj) == len(xy), name
            _CACHE.setdefault(len(adj), []).append((adj, xy, hand))
    return _CACHE


def template_block(mol, atoms, rings, rank):
    """`({atom: axial point}, handed)` from the first template the block is, or None."""
    nb = {x: [y for y in mol.neighbors_of(x) if y in atoms] for x in atoms}
    bonds = sum(len(v) for v in nb.values()) // 2
    by_rank = sorted(atoms, key=rank.__getitem__)
    for adj, xy, hand in _skeletons().get(len(atoms), ()):
        if sum(len(v) for v in adj) // 2 != bonds:
            continue
        m = _match(adj, nb, by_rank)
        if m is not None:
            return {m[i]: (x / H, y - x / H / 2) for i, (x, y) in enumerate(xy)}, hand   # bonds are L: axial
    return None


def _match(adj, nb, by_rank):
    """Skeleton index to block atom, candidates tried in rank order so the first full match ranks lowest; equal
    bond counts make an edge-preserving injection an isomorphism."""
    n = len(adj)
    m = [None] * n
    used = set()

    def go(i):
        if i == n:
            return True
        back = [m[j] for j in adj[i] if j < i]
        for x in by_rank:
            if x in used or len(nb[x]) != len(adj[i]) or any(y not in nb[x] for y in back):
                continue
            m[i] = x
            used.add(x)
            if go(i + 1):
                return True
            used.discard(x)
        m[i] = None
        return False
    return m if go(0) else None


def crown(n, period):
    """`n` points of a crown in ring order, every bond one, the ring starting at a heteroatom: walking it, each
    carbon turns 60 degrees one way and each heteroatom the other way by what closes the ring, `60 (period - 1) -
    360 / units`; so 18-crown-6 is all 120 degrees, its oxygens the inward corners."""
    units = n // period
    hetero = radians(360 / units - 60 * (period - 1))
    pts, (x, y), h = [], (0., 0.), 0.
    for i in range(n):
        pts.append((x, y))
        x, y = x + cos(h), y + sin(h)
        h += hetero if (i + 1) % period == 0 else radians(60)
    cx, cy = sum(p[0] for p in pts) / n, sum(p[1] for p in pts) / n
    return [(p[0] - cx, p[1] - cy) for p in pts]
