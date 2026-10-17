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
"""The final orientation of a drawn component, and the packing of components.

| step | rule |
| --- | --- |
| upright | a mirror line off every lattice direction turns to vertical or horizontal, whichever draws wider; one on a
  diagonal lattice direction turns by a multiple of 60 degrees, or of 30 when acyclic or through a spiro atom, and
  only when the drawing is no narrower; a tie goes to a vertical bond, and a lone polygon of 7 or more off the
  lattice stands on an edge |
| level | a relaxed drawing `LEVEL` tall and over `TALL` times taller than wide turns by the multiple of 60 degrees
  drawing it widest |
| horizontal | a long drawing, spread `LONG` times more along its main axis than across and `SPAN` along it, turns by
  the multiple of 30 degrees laying that axis nearest horizontal when that at least halves the tilt; the pose's ring
  key picks it or its half turn |
| straight | a short acyclic drawing, under `SPAN` across, lays its one bond or the axis of its first linear atom
  horizontal |
| components | left to right by rank, `GAP` apart |
"""
from itertools import combinations
from math import atan2, cos, degrees, hypot, pi, radians, sin

from .lattice import L
from .quality import linear_atoms


__all__ = ['GAP', 'level', 'horizontal', 'lean', 'straight', 'upright']

GAP = 2 * L
LEVEL = 8 * L           # height from which a drawing taller than wide is levelled
TALL = 1.25             # height to width ratio over which it is
LONG = 1.25             # spread along the main axis over the spread across it from which a drawing is long
SPAN = 2 * L            # and its extent along that axis from which it is; a spiro pair of rings from `SPIRO`
SPIRO = 1.5 * L


def level(pts):
    """Turn `pts` about the origin by the multiple of 60 degrees drawing it widest, when it is `LEVEL` tall and over
    `TALL` times taller than wide."""
    def box(k):
        c, s = cos(radians(60 * k)), sin(radians(60 * k))
        xs = [x * c - y * s for x, y in pts.values()]
        ys = [x * s + y * c for x, y in pts.values()]
        return max(xs) - min(xs), max(ys) - min(ys)

    w, h = box(0)
    if h < LEVEL or h < TALL * w:
        return
    k = max(range(6), key=lambda k: (round(box(k)[0] - box(k)[1], 3), -k))
    c, s = cos(radians(60 * k)), sin(radians(60 * k))
    for a, (x, y) in pts.items():
        pts[a] = (x * c - y * s, x * s + y * c)


def horizontal(pts, rings, span=SPAN):
    """Turn `pts` about its centroid by the multiple of 30 degrees bringing its main axis, the direction it spreads
    most along, nearest horizontal, when it is long - spread `LONG` times more along that axis than across it and
    `SPAN` from end to end - and the turn at least halves the tilt. Of the turn and the half turn after it, the one
    the pose's `_lean` of `rings` puts first wins."""
    n = len(pts)
    if n < 3:
        return
    cx = sum(p[0] for p in pts.values()) / n
    cy = sum(p[1] for p in pts.values()) / n
    xx = sum((p[0] - cx) ** 2 for p in pts.values())
    yy = sum((p[1] - cy) ** 2 for p in pts.values())
    xy = sum((p[0] - cx) * (p[1] - cy) for p in pts.values())
    mid, half = (xx + yy) / 2, hypot((xx - yy) / 2, xy)
    if mid - half <= 0 or mid + half < LONG * LONG * (mid - half):
        return
    t = atan2(2 * xy, xx - yy) / 2
    k = round(t / radians(30))
    if not k or 2 * abs(t - radians(30 * k)) > abs(t):
        return
    u, v = cos(t), sin(t)
    along = [(p[0] - cx) * u + (p[1] - cy) * v for p in pts.values()]
    if max(along) - min(along) < span:
        return
    c, s = cos(radians(-30 * k)), sin(radians(-30 * k))
    turned = {a: ((x - cx) * c - (y - cy) * s, (x - cx) * s + (y - cy) * c) for a, (x, y) in pts.items()}
    best = None
    for f in (1, -1):
        q = {a: (f * x, f * y) for a, (x, y) in turned.items()}
        key = _lean(q, rings)
        if best is None or key < best[0]:
            best = key, q
    for a, (x, y) in best[1].items():
        pts[a] = (cx + x, cy + y)


def straight(mol, pts, rank):
    """Turn a short acyclic `pts`, under `SPAN` across, about its centroid so its bond, when it has one, or the axis
    through the linear atom of least rank lies horizontal: an azide, a nitrile, nitric oxide."""
    if len(pts) == 2:
        a, b = pts
    else:
        axis = sorted(linear_atoms(mol) & set(pts), key=rank.__getitem__)
        if not axis:
            return
        a, b = mol.neighbors_of(axis[0])
    if max(hypot(p[0] - q[0], p[1] - q[1]) for p, q in combinations(pts.values(), 2)) >= SPAN:
        return
    t = atan2(pts[b][1] - pts[a][1], pts[b][0] - pts[a][0])
    t = (t + pi / 2) % pi - pi / 2
    if abs(t) < 1e-9:
        return
    n = len(pts)
    cx = sum(p[0] for p in pts.values()) / n
    cy = sum(p[1] for p in pts.values()) / n
    c, s = cos(-t), sin(-t)
    for x, (px, py) in pts.items():
        pts[x] = (cx + (px - cx) * c - (py - cy) * s, cy + (px - cx) * s + (py - cy) * c)


def lean(pts, rings):
    """Of `pts`, its two mirror images about its centroid and its half turn, the one the pose's `_lean` of `rings` puts
    first, `pts` itself on a tie: a steroid with ring D to the upper right."""
    if not rings:
        return
    n = len(pts)
    cx = sum(p[0] for p in pts.values()) / n
    cy = sum(p[1] for p in pts.values()) / n
    best = None
    for k, (fx, fy) in enumerate(((1, 1), (1, -1), (-1, 1), (-1, -1))):
        q = {a: (fx * (x - cx), fy * (y - cy)) for a, (x, y) in pts.items()}
        key = (_lean(q, rings), k)
        if best is None or key < best[0]:
            best = key, q
    for a, (x, y) in best[1].items():
        pts[a] = (cx + x, cy + y)


def _lean(pts, rings):
    """The pose's ring key, least first: fewest centres of `rings` to the upper right of their mean negated, then
    the smaller rings' pull that way negated."""
    if not rings:
        return 0, 0
    cs = [(sum(pts[x][0] for x in r) / len(r), sum(pts[x][1] for x in r) / len(r)) for r in rings]
    mx = sum(c[0] for c in cs) / len(cs)
    my = sum(c[1] for c in cs) / len(cs)
    return (-sum(c[0] - mx > .01 and c[1] - my > .01 for c in cs),
            -round(sum((c[0] - mx + c[1] - my) / len(r) for c, r in zip(cs, rings)), 3))


def _lattice(p, q):
    """The bond pq lies on a multiple of 30 degrees."""
    t = degrees(atan2(q[1] - p[1], q[0] - p[0]))
    return min(t % 30, -t % 30) < 2


def upright(mol, atoms, pts, orbit):
    """Turn `pts` about its centroid so a mirror line through it, mapping each atom onto one of its orbit and bonds
    onto bonds, lies vertical or horizontal; untouched when a mirror line already does. A mirror on a diagonal lattice
    direction turns only by a lattice turn, 30 degrees when `atoms` is acyclic or the mirror passes a spiro atom and
    60 otherwise, and only when that draws no narrower. Among turns drawing equally wide, one giving a vertical bond
    wins; a lone ring of 7 or more with no bond on the lattice turns freely, so an octagon stands on an edge."""
    if orbit is None or len(atoms) < 3:
        return
    n = len(pts)
    cx = sum(p[0] for p in pts.values()) / n
    cy = sum(p[1] for p in pts.values()) / n
    by = {}
    for x in pts:
        by.setdefault(orbit[x], []).append(x)
    a = max(min(by.values(), key=len), key=lambda x: (hypot(pts[x][0] - cx, pts[x][1] - cy), -x))
    ax, ay = pts[a][0] - cx, pts[a][1] - cy
    if hypot(ax, ay) < .1 * L:
        return
    bonds = {(b.n, b.m) for b in mol.bonds() if b.n in atoms}
    bonds |= {(m, n) for n, m in bonds}
    tol = .1 * L
    axes = []
    for b in by[orbit[a]]:
        if b == a:
            t = atan2(ay, ax)
        else:
            t = atan2(pts[b][0] - pts[a][0], pts[a][1] - pts[b][1])     # the perpendicular of ab
        u, v = cos(t), sin(t)
        image = {}
        for x, (px, py) in pts.items():
            dx, dy = px - cx, py - cy
            d = 2 * (dx * u + dy * v)
            rx, ry = cx + d * u - dx, cy + d * v - dy
            y = next((y for y in by[orbit[x]] if abs(pts[y][0] - rx) < tol and abs(pts[y][1] - ry) < tol), None)
            if y is None:
                break
            image[x] = y
        else:
            if all((image[n], image[m]) in bonds for n, m in bonds):
                axes.append(t)
    if not axes:
        return
    if len(atoms) > 6 and len(bonds) == 2 * len(atoms) and not any(_lattice(pts[n], pts[m]) for n, m in bonds):
        on = []                                             # a lone polygon off the lattice turns onto an edge
    else:
        on = [t for t in axes if min(degrees(t) % 30, -degrees(t) % 30) < 2]
        if any(min(degrees(t) % 90, -degrees(t) % 90) < 2 for t in on):
            return
    best = None
    if on:                                                  # a lattice turn, and only one no narrower
        xs = [p[0] for p in pts.values()]
        ys = [p[1] for p in pts.values()]
        best = (-round(max(xs) - min(xs) - (max(ys) - min(ys)), 3) + 1e-3, False, 2), None, None
    ringed = any(mol.in_ring_of(x) for x in atoms)
    spiro = [x for x in atoms if len(mol.neighbors_of(x)) == 4 and mol.in_ring_of(x) and
             all(len([z for z in mol.neighbors_of(y) if mol.in_ring_of(z)]) == 2 for y in mol.neighbors_of(x))]
    for t in on or axes:
        free = not ringed or any(abs((pts[x][0] - cx) * sin(t) - (pts[x][1] - cy) * cos(t)) < tol for x in spiro)
        for k, goal in enumerate((pi / 2, 0.)):
            if on and abs(sin((6 if free else 3) * (goal - t))) > .05:
                continue
            c, s = cos(goal - t), sin(goal - t)
            xs = [c * (p[0] - cx) - s * (p[1] - cy) for p in pts.values()]
            ys = [s * (p[0] - cx) + c * (p[1] - cy) for p in pts.values()]
            flat = any(abs((c * (pts[n][0] - pts[m][0]) - s * (pts[n][1] - pts[m][1]))) < 1e-3 * L for n, m in bonds)
            key = (-round(max(xs) - min(xs) - (max(ys) - min(ys)), 3), not flat, k)
            if best is None or key < best[0]:
                best = key, c, s
    _, c, s = best
    if c is None:
        return
    for x, (px, py) in pts.items():
        dx, dy = px - cx, py - cy
        pts[x] = (cx + c * dx - s * dy, cy + s * dx + c * dy)
