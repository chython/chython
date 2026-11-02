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
"""The depth of a perspective cage: which face of it a stereocentre's substituent is on, which bond is in front.

A cage `shapes` draws is a projection, so its geometry can state a configuration a flat wedge would contradict:
with depth known, a substituent drawn on one face of its cage atom is one configuration and on the other face the
other.  The rear bond of each crossing is drawn broken, which is what makes the depth readable.

| part | rule |
| --- | --- |
| point | a cage atom in 3D: its plane point and its depth, `z` toward the viewer, `sign` choosing the reading |
| faces | of a cage atom of two cage bonds, its two tetrahedral directions, `TETRA` either side of the plane of those \
bonds about their outward bisector; of a bridgehead, of three, the one opposite their sum |
| reading | `reads`: the configuration drawn at a cage atom: each drawn substituent on the face whose projection lies \
nearest its bond, an implicit hydrogen on the face left, the parity in 3D by `core.wedge.tetrahedral_parity`; None \
when two faces project within `SPLIT` of each other or one shorter than `SHORT`, a ring seen edge or face on, or a \
drawn bond is not `READ` nearer its face than the other, between them |
| sign | `sign`: the depth under which the cage's bridgehead stereocentres read as stated; on a tie, with no one-atom \
bridge in front already, the one bringing the substituted atoms forward: a bridgehead has one face, so the cage alone \
decides it; a bridgehead coupled to one that reads as stated is shown with it, the \
cage stating both |
| crossing | a substituted stereocentre within `CROSS` of a crossing of two cage bonds is not read from its face, its \
drawn bond lost in the crossing, and takes a wedge with every stereocentre coupled to it |
| place | `place`: for each substituent of a cage stereocentre, the direction in the plane that reads as stated |
| found | `cages`: each shape drawing found in a plane, a fused system's core alone, rotated, scaled and possibly \
reflected, within `MATCH`; its depth sign the one under which most stereocentres read as stated; those it `shows`, \
and the rear bond of each crossing of two of its bonds with the crossing point |
"""
from itertools import islice
from math import atan2, cos, hypot, pi, radians, sin

from ...core.wedge import _Planar, tetrahedral_parity
from .lattice import L

TETRA = radians(54.75)  # half the tetrahedral angle: an exo bond's tilt off its atom's ring bond plane
SPLIT = radians(30)     # the least angle the two faces' projections leave apart to be told apart
SHORT = .3              # the shortest projection of a face, in bonds, a reader sees a direction in
READ = radians(15)      # a drawn substituent this much nearer its face than the other reads, else it is between
CROSS = .55             # a substituted stereocentre this near a crossing of the cage's bonds, in bonds, is unread
MATCH = .02 * L         # root mean square distance of a plane's cage atoms from a shape drawing found in it


def _unit(v):
    n = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]) ** .5
    return (v[0] / n, v[1] / n, v[2] / n) if n > 1e-12 else None


def _point(plane, z, x, sign):
    return plane[x][0], plane[x][1], sign * z[x]


def faces(mol, plane, z, a, sign=1):
    """The tetrahedral directions out of the cage at its atom `a`, unit 3D vectors: two at an atom of two cage bonds,
    one at a bridgehead; None for any other atom of `z` or a degenerate geometry."""
    ring = [y for y in mol.neighbors_of(a) if y in z]
    pa = _point(plane, z, a, sign)
    us = [_unit(tuple(q - p for p, q in zip(pa, _point(plane, z, y, sign)))) for y in ring]
    if None in us:
        return None
    if len(ring) == 3:
        e = _unit(tuple(-sum(u[k] for u in us) for k in range(3)))
        return None if e is None else [e]
    if len(ring) != 2:
        return None
    (ux, uy, uz), (vx, vy, vz) = us
    b = _unit((-ux - vx, -uy - vy, -uz - vz))
    n = _unit((uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx))
    if b is None or n is None:
        return None
    c, s = cos(TETRA), sin(TETRA)
    return [tuple(b[k] * c + n[k] * s for k in range(3)), tuple(b[k] * c - n[k] * s for k in range(3))]


def _readable(fs):
    if any(hypot(e[0], e[1]) < SHORT for e in fs):
        return False
    if len(fs) == 2:
        d = abs((atan2(fs[0][1], fs[0][0]) - atan2(fs[1][1], fs[1][0]) + pi) % (2 * pi) - pi)
        return d >= SPLIT
    return True


def _parity(mol, unit, plane, z, a, assigned, sign):
    """The parity of `unit` at cage atom `a` with each exo atom of `assigned` on its face, in 3D."""
    pa = _point(plane, z, a, sign)
    pts = {a: pa}
    for y in mol.neighbors_of(a):
        if y in z:
            pts[y] = _point(plane, z, y, sign)
    for y, e in assigned.items():
        pts[y] = tuple(p + L * v for p, v in zip(pa, e))
    return tetrahedral_parity(_Planar(mol, {x: p[:2] for x, p in pts.items()}), unit,
                              {x: p[2] for x, p in pts.items()})


def _assignments(fs, exo):
    """Each way of putting the atoms `exo` on distinct faces of `fs`."""
    if len(exo) == 1:
        return [{exo[0]: e} for e in fs]
    if len(exo) == 2 and len(fs) == 2:
        return [{exo[0]: fs[0], exo[1]: fs[1]}, {exo[0]: fs[1], exo[1]: fs[0]}]
    return []


def _units(mol):
    from ...core import SU_TETRA
    return {u['anchor']: u for u in mol.stereo_units() if u['kind'] == SU_TETRA and mol.parity_of(u['anchor'])}


def reads(mol, plane, z, a, sign=1, units=None):
    """The parity the drawing `plane` shows at the cage atom `a` under depth `sign * z`, each drawn substituent on the
    face nearest its bond; None when its faces cannot be told apart, or `a` has no substituent drawn and is no
    bridgehead."""
    unit = (units if units is not None else _units(mol)).get(a)
    fs = faces(mol, plane, z, a, sign)
    exo = [y for y in mol.neighbors_of(a) if y not in z]
    if unit is None or fs is None:
        return None
    if not exo:                                         # a bridgehead's hydrogen: the cage alone states it
        return _parity(mol, unit, plane, z, a, {}, sign) if len(fs) == 1 else None
    if not _readable(fs):
        return None
    ax, ay = plane[a]

    def miss(assigned):
        return sum(abs((atan2(plane[y][1] - ay, plane[y][0] - ax) - atan2(e[1], e[0]) + pi) % (2 * pi) - pi)
                   for y, e in assigned.items())
    options = _assignments(fs, exo)
    if not options:
        return None
    best = min(options, key=miss)
    for y, e in best.items():                           # each drawn bond plainly nearer its face than the other
        t = atan2(plane[y][1] - ay, plane[y][0] - ax)
        off = [abs((t - atan2(f[1], f[0]) + pi) % (2 * pi) - pi) for f in fs]
        own = abs((t - atan2(e[1], e[0]) + pi) % (2 * pi) - pi)
        if any(d - own < READ for f, d in zip(fs, off) if f is not e):
            return None
    return _parity(mol, unit, plane, z, a, best, sign)


def sign(mol, plane, z, units=None):
    """The depth sign under which the cage's bridgehead stereocentres read as stated; when none decides and the cage
    has no one-atom bridge, which the drawing already puts in front, the one bringing its substituted atoms forward;
    else `+1`."""
    units = _units(mol) if units is None else units
    votes = 0
    for a in z:
        if a in units and sum(y in z for y in mol.neighbors_of(a)) == 3:
            for s in (1, -1):
                if reads(mol, plane, z, a, s, units) == mol.parity_of(a):
                    votes += s
    if votes:
        return -1 if votes < 0 else 1
    heads = {x for x in z if sum(y in z for y in mol.neighbors_of(x)) == 3}
    if any(x not in heads and set(y for y in mol.neighbors_of(x) if y in z) <= heads for x in z):
        return 1                                        # a one-atom bridge: in front already
    forward = sum(z[x] for x in z if any(y not in z for y in mol.neighbors_of(x)))
    return -1 if forward < -1e-9 else 1


def place(mol, plane, z, a, sign=1, units=None):
    """`{substituent: angle}`: the plane direction for each drawn substituent of the cage stereocentre `a` that reads
    as its stated parity under depth `sign * z`; None when no placement does or its faces cannot be told apart."""
    unit = (units if units is not None else _units(mol)).get(a)
    fs = faces(mol, plane, z, a, sign)
    exo = sorted(y for y in mol.neighbors_of(a) if y not in z)
    if unit is None or fs is None or len(fs) != 2 or not exo or not _readable(fs):
        return None
    for assigned in _assignments(fs, exo):
        if _parity(mol, unit, plane, z, a, assigned, sign) == mol.parity_of(a):
            return {y: atan2(e[1], e[0]) for y, e in assigned.items()}
    return None


class Cage:
    """A shape drawing found in a plane: its atoms' depth, the stereocentres its geometry shows, the rear bonds."""
    __slots__ = ('atoms', 'z', 'shows', 'gaps')

    def __init__(self, atoms, z, shows, gaps):
        self.atoms = atoms
        self.z = z
        self.shows = shows          # stereocentres whose drawing reads as stated: no wedge needed
        self.gaps = gaps            # [((n, m), (x, y))]: the rear bond of a crossing, and where it is crossed


def cages(mol, plane):
    """Each cage `shapes` draws, found in `plane`; empty when none is."""
    from .rings import ring_blocks
    from .shapes import MAPPINGS, _shapes, fused3d, skeleton
    out = []
    units = None
    for block, rings in ring_blocks(mol):
        if not _shapes().get(skeleton(mol, block)[0].canonical_bytes) and not 3 <= len(rings):
            continue
        found = fused3d(mol, block, rings)
        if found is None:
            continue
        skel, order = skeleton(mol, found[0])
        found = _shapes().get(skel.canonical_bytes)
        if not found:
            continue
        best = None
        for _, ref, pts, _, _, _, depth in found:
            if not any(depth.values()):
                continue
            for mapping in islice(ref.as_query().get_mapping(skel), MAPPINGS):
                atoms = [order[mapping[i] - 1] for i in pts]
                for flip in (1., -1.):
                    fit = _similar([(flip * pts[i][0], pts[i][1]) for i in pts], [plane[x] for x in atoms])
                    if fit is not None and (best is None or fit[0] < best[0]):
                        best = fit[0], {x: flip * fit[1] * depth[i] for x, i in zip(atoms, pts)}
        if best is None or best[0] > MATCH:
            continue
        z = best[1]
        if units is None:
            units = _units(mol)
        s = sign(mol, plane, z, units)
        gaps = _rear(mol, plane, {a: s * v for a, v in z.items()})
        hidden = {a for a in z if any(y not in z for y in mol.neighbors_of(a))      # a drawn substituent's centre
                  for _, (x, y) in gaps if hypot(plane[a][0] - x, plane[a][1] - y) < CROSS * L}
        hidden |= {a for a in z if a in units and set(units[a]['coupled']) & hidden} | \
            {c for a in hidden if a in units for c in units[a]['coupled']}             # coupled: wedged together
        shows = {a for a in z if a in units and a not in hidden and
                 reads(mol, plane, z, a, s, units) == mol.parity_of(a)}
        shows |= {a for a in z if a in units and a not in hidden and set(units[a]['coupled']) & shows}  # its partner
        out.append(Cage(frozenset(z), {a: s * v for a, v in z.items()}, shows, gaps))
    return out


def _similar(src, dst):
    """`(rms, scale)` of the rotation and scale taking the points `src` nearest `dst`, translation free."""
    n = len(src)
    sx, sy = sum(p[0] for p in src) / n, sum(p[1] for p in src) / n
    dx, dy = sum(p[0] for p in dst) / n, sum(p[1] for p in dst) / n
    a = b = norm = 0.
    for (px, py), (qx, qy) in zip(src, dst):
        px, py, qx, qy = px - sx, py - sy, qx - dx, qy - dy
        a += px * qx + py * qy
        b += px * qy - py * qx
        norm += px * px + py * py
    if norm < 1e-12:
        return None
    scale = (a * a + b * b) ** .5 / norm
    c, s = a / (a * a + b * b) ** .5, b / (a * a + b * b) ** .5
    err = sum((qx - dx - scale * (c * (px - sx) - s * (py - sy))) ** 2 +
              (qy - dy - scale * (s * (px - sx) + c * (py - sy))) ** 2 for (px, py), (qx, qy) in zip(src, dst))
    return (err / n) ** .5, scale


def _rear(mol, plane, z):
    """`[((n, m), point)]`: of each two bonds of the cage `z` crossing in `plane`, the one behind, and where it is
    crossed."""
    bonds = sorted({(min(x, y), max(x, y)) for x in z for y in mol.neighbors_of(x) if y in z})
    out = []
    for i, (a, b) in enumerate(bonds):
        for c, d in bonds[i + 1:]:
            if {a, b} & {c, d}:
                continue
            (x1, y1), (x2, y2), (x3, y3), (x4, y4) = plane[a], plane[b], plane[c], plane[d]
            den = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)
            if abs(den) < 1e-12:
                continue
            t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / den
            u = ((x1 - x3) * (y1 - y2) - (y1 - y3) * (x1 - x2)) / den
            if not (0 < t < 1 and 0 < u < 1):
                continue
            front = z[a] + t * (z[b] - z[a]) >= z[c] + u * (z[d] - z[c])
            out.append(((c, d) if front else (a, b), (x1 + t * (x2 - x1), y1 + t * (y2 - y1))))
    return out


__all__ = ['Cage', 'cages', 'faces', 'place', 'reads', 'sign']
