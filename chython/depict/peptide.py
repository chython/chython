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
"""Peptide cross-links drawn as routed lines, derived from `(mol, plane)` alone.

A link is routed only when the backbone lies on the 120-degree lattice and the link is longer than
`ROUTE` -- any other plane keeps it a plain bond, so a stored peptide plane draws the same way again.

| Kind | Drawing |
| --- | --- |
| `bracket` | both ends on one side of the backbone: a leg from each, a bar `BAR` past the atoms it spans |
| `head_to_tail` | from each end's exit down to a bar `BAR` below every atom |
| `teleport` | a dashed stub ending in a numbered badge at each end |

Corners are circular arcs, one cubic each: a sharp bend would read as a carbon.  A bracket leg runs straight
along its end's free direction -- a ring's exterior bisector, a chain's 120-degree slot -- into the bar, and
rounds at `BEND` whatever the angle; a leg too shallow, long or blocked is a stub of `JOG` then a vertical.
"""
from collections import deque
from math import acos, atan2, cos, degrees, hypot, pi, radians, sin, tan
from typing import NamedTuple

from ..core.monomers import monomers
from .bonds import ray_box_exit, segment_hits_box
from .scene import Path, Text, TextRun, circle, curve, line, move


__all__ = ['Link', 'link_paths', 'routed_links', 'rounded']

L = .825
ROUTE = 1.5 * L         # a shorter link stays a bond
BAR = .6                # bar clearance past the outermost atom
NEST = .45              # an enclosing bracket's step outward
STUB = .55              # teleport stub
BADGE = .22             # teleport badge radius
BADGE_TEXT = .3
CORNER = .25
BEND = .5               # a bracket's corner radius
LEG = 1.5 * L           # the longest horizontal run of a straight leg
JOG = .5 * L            # a leg's stub before it turns vertical
CLOSE = .3              # a straight leg's clearance from an atom
DASHES = (.08, .06)


class Link(NamedTuple):
    """One routed cross-link: atoms `a`, `b`, its `kind`, and `side` +1 above / -1 below for a bracket."""
    a: int
    b: int
    kind: str
    side: int = 0


def _on_lattice(plane, backbone):
    """Every backbone bond of length `L`, or `2 L` for a bridge, at 30 deg + a multiple of 60, to the rounding
    of a stored plane."""
    for a, b in zip(backbone, backbone[1:]):
        dx, dy = plane[b][0] - plane[a][0], plane[b][1] - plane[a][1]
        d = degrees(atan2(dy, dx)) - 30
        if abs(d - 60 * round(d / 60)) > .05 or min(abs(hypot(dx, dy) - s) for s in (L, 2 * L)) > .01:
            return False
    return True


def _host(mol, x, idx, links):
    """The first backbone atom BFS reaches from `x`, never across a cross-link."""
    seen, queue = {x}, deque([x])
    while queue:
        y = queue.popleft()
        if y in idx:
            return y
        for z in mol.neighbors_of(y):
            if z not in seen and frozenset((y, z)) not in links:
                seen.add(z)
                queue.append(z)
    return None


def routed_links(mol, plane) -> tuple:
    """The cross-links of `mol` to route in `plane`, or `()` when the plane is not a peptide layout: each
    component on its own."""
    components = mol.connected_components
    if len(components) == 1:
        return _routed(mol, plane)
    return tuple(k for c in components if len(c) > 1 for k in _routed(mol.substructure(c), plane))


def _routed(mol, plane):
    pep = monomers(mol)
    if pep is None or not _on_lattice(plane, pep.backbone):
        return ()
    idx = {a: q for q, a in enumerate(pep.backbone)}
    keys = frozenset(frozenset(x) for x in pep.crosslinks)
    out, spans = [], []
    for q, (a, b) in enumerate(pep.crosslinks):
        if hypot(plane[a][0] - plane[b][0], plane[a][1] - plane[b][1]) <= ROUTE:
            continue
        if pep.head_to_tail and q == len(pep.crosslinks) - 1:
            out.append(Link(a, b, 'head_to_tail'))
            continue
        ha, hb = _host(mol, a, idx, keys), _host(mol, b, idx, keys)
        if ha is None or hb is None:
            out.append(Link(a, b, 'teleport'))
            continue
        sa = (plane[a][1] > plane[ha][1]) - (plane[a][1] < plane[ha][1])
        sb = (plane[b][1] > plane[hb][1]) - (plane[b][1] < plane[hb][1])
        ka, kb = sorted((idx[ha], idx[hb]))
        spans.append((kb - ka, ka, kb, a, b, sa if sa == sb else 0))
    accepted = []
    for _, ka, kb, a, b, s in sorted(spans):
        if s and not any(s == s2 and (k1 < ka < k2 < kb or ka < k1 < kb < k2) for k1, k2, s2 in accepted):
            accepted.append((ka, kb, s))
            out.append(Link(a, b, 'bracket', s))
        else:
            out.append(Link(a, b, 'teleport'))
    return tuple(out)


def rounded(points, r=CORNER):
    """A polyline with every corner a circular-arc cubic of radius `r`, clipped to half the shorter side."""
    pts = [p for i, p in enumerate(points)
           if not i or hypot(p[0] - points[i - 1][0], p[1] - points[i - 1][1]) > 1e-9]
    out = [move(*pts[0])]
    for i in range(1, len(pts) - 1):
        p, c, q = pts[i - 1], pts[i], pts[i + 1]
        u = _norm((c[0] - p[0], c[1] - p[1]))
        v = _norm((q[0] - c[0], q[1] - c[1]))
        th = acos(max(-1., min(1., u[0] * v[0] + u[1] * v[1])))    # the turning angle
        if th < 1e-3:
            out.append(line(*c))
            continue
        t = min(r * tan(th / 2), hypot(c[0] - p[0], c[1] - p[1]) / 2, hypot(q[0] - c[0], q[1] - c[1]) / 2)
        k = 4 / 3 * tan(th / 4) * t / tan(th / 2)                     # cubic handle of an arc of angle th
        a = (c[0] - u[0] * t, c[1] - u[1] * t)
        b = (c[0] + v[0] * t, c[1] + v[1] * t)
        out.append(line(*a))
        out.append(curve(a[0] + u[0] * k, a[1] + u[1] * k, b[0] - v[0] * k, b[1] - v[1] * k, *b))
    out.append(line(*pts[-1]))
    return tuple(out)


def _norm(v):
    n = hypot(*v) or 1.
    return v[0] / n, v[1] / n


def _exit(mol, plane, backbone, x, links):
    """Where a head-to-tail line leaves `x`: N straight down, C at 120 degrees to its bond, the lower side."""
    p = plane[x]
    if x == backbone[0]:
        d = 270.
    elif x == backbone[-1]:
        q = plane[backbone[-2]]
        last = degrees(atan2(p[1] - q[1], p[0] - q[0]))
        d = min((last + 60, last - 60), key=lambda t: sin(radians(t)))
    else:
        dx, dy = _stub_direction(mol, plane, x, links)
        d = degrees(atan2(dy, dx))
    return p[0] + L * cos(radians(d)), p[1] + L * sin(radians(d))


def _headings(mol, plane, x, s, toward, links):
    """The free directions of `x`: those steep toward side `s` first, leaning toward the x of `toward` before
    leaning away, so a bracket narrows to its bar; then by steepness.  Those pointing away from `s` dropped
    unless none is left."""
    p = plane[x]
    steep = sin(radians(25))
    ds = sorted(_free_directions(mol, plane, x, links),
                key=lambda d: (s * d[1] <= steep, d[0] * (toward - p[0]) < -1e-6, -round(s * d[1], 6)))
    return [d for d in ds if s * d[1] > 1e-6] or ds[:1]


def _leg(mol, plane, x, us, level, ink):
    """Points after `x` to the bar at height `level`: straight along a direction of `us` when that is steep,
    short and clear, else a stub of `JOG` along one and then vertical, the one crossing least, the order of
    `us` breaking a tie.  `JOG` is half a bond past `x`'s label, so the vertical runs between the lattice's
    columns."""
    p = plane[x]
    for u in us:
        if abs(u[1]) > sin(radians(25)):
            t = (level - p[1]) / u[1]
            q = (p[0] + u[0] * t, p[1] + u[1] * t)
            if t > 0 and abs(q[0] - p[0]) <= LEG and not _hits(mol, plane, x, p, q, ink):
                return [q]
    best = None
    for u in us:
        j = JOG + (ray_box_exit(p, u, ink[x]) if x in ink else 0.)
        e = (p[0] + u[0] * j, p[1] + u[1] * j)
        n = _hits(mol, plane, x, p, e, ink) + _hits(mol, plane, x, e, (e[0], level), ink)
        if best is None or n < best[0]:
            best = n, e
    e = best[1]
    return [e, (e[0], level)]


def _hits(mol, plane, x, p, q, ink):
    """How many label boxes, bonds and atoms, other than `x`'s, the segment `p`-`q` crosses or passes within
    `CLOSE` of."""
    n = sum(y != x and segment_hits_box(p, q, bx) for y, bx in ink.items())
    dx, dy = q[0] - p[0], q[1] - p[1]
    d = dx * dx + dy * dy or 1.
    for y, (cx, cy) in plane.items():
        if y != x:
            t = max(0., min(1., ((cx - p[0]) * dx + (cy - p[1]) * dy) / d))
            n += hypot(p[0] + t * dx - cx, p[1] + t * dy - cy) < CLOSE
    return n + sum(x not in (b.n, b.m) and _cross(p, q, plane[b.n], plane[b.m]) for b in mol.bonds())


def _cross(p, q, a, b):
    """Do the segments `p`-`q` and `a`-`b` properly intersect?"""
    def side(o, u, v):
        return (u[0] - o[0]) * (v[1] - o[1]) - (u[1] - o[1]) * (v[0] - o[0])
    return side(p, q, a) * side(p, q, b) < 0 and side(a, b, p) * side(a, b, q) < 0


def _free_directions(mol, plane, x, links):
    """Straight on from a single vertical bond, then 120 degrees either side of it; else the bisector of each gap
    between drawn bonds."""
    p = plane[x]
    angs = sorted(atan2(plane[y][1] - p[1], plane[y][0] - p[0]) for y in mol.neighbors_of(x)
                  if frozenset((x, y)) not in links)
    if not angs:
        return [(0., 1.), (0., -1.)]
    if len(angs) == 1:
        ts = [angs[0] + radians(120), angs[0] - radians(120)]
        if abs(cos(angs[0])) < 1e-6:
            ts.insert(0, angs[0] + pi)
    else:
        ts = [angs[q] + ((angs[(q + 1) % len(angs)] - angs[q]) % (2 * pi)) / 2 for q in range(len(angs))]
    return [(cos(t), sin(t)) for t in ts]


def _stub_direction(mol, plane, x, links):
    """120 degrees from a single drawn bond, the side with more room; else the widest gap's bisector."""
    p = plane[x]
    angs = sorted(atan2(plane[y][1] - p[1], plane[y][0] - p[0]) for y in mol.neighbors_of(x)
                  if frozenset((x, y)) not in links)
    if not angs:
        return 0., 1.
    if len(angs) == 1:
        candidates = [angs[0] + radians(120), angs[0] - radians(120)]
    else:
        gap, lo = max(((angs[(q + 1) % len(angs)] - angs[q]) % (2 * pi), angs[q]) for q in range(len(angs)))
        candidates = [lo + gap / 2]

    def room(t):
        e = (p[0] + cos(t) * .9, p[1] + sin(t) * .9)
        return min((hypot(e[0] - q[0], e[1] - q[1]) for y, q in plane.items() if y != x), default=9.)
    t = max(candidates, key=room)
    return cos(t), sin(t)


def link_paths(mol, plane, links, boxes, style) -> list:
    """`Path` and `Text` nodes for `links`, the line ends trimmed at the label boxes like bonds: each
    component on its own, the teleport numbers running on."""
    if not links:
        return []
    components = mol.connected_components
    if len(components) == 1:
        return _paths(mol, plane, links, boxes, style, 0)[0]
    out, number = [], 0
    for c in components:
        atoms = set(c)
        own = [k for k in links if k.a in atoms]
        if own:
            sub = {x: plane[x] for x in c}
            nodes, number = _paths(mol.substructure(c), sub, own, {x: b for x, b in boxes.items() if x in atoms},
                                   style, number)
            out.extend(nodes)
    return out


def _paths(mol, plane, links, boxes, style, number):
    """`link_paths` of one component, teleports numbered after `number`: `(nodes, last number)`."""
    bond = style.bond
    colour, width = bond.colour, bond.width
    background = style.page.background or '#ffffff'
    keys = frozenset(frozenset((k.a, k.b)) for k in links)
    ink = {x: lb.box for x, lb in boxes.items() if lb.text is not None}
    ys = [y for _, y in plane.values()]
    backbone = monomers(mol).backbone
    out, levels = [], []

    def clip(pts, x0, x1):
        pts = list(pts)
        for i, j, x in ((0, 1, x0), (-1, -2, x1)):
            if x in ink:
                p, q = pts[i], pts[j]
                u = _norm((q[0] - p[0], q[1] - p[1]))
                t = ray_box_exit(p, u, ink[x]) + bond.trim
                pts[i] = (p[0] + u[0] * t, p[1] + u[1] * t)
        return pts

    def draw(pts, dashes=None):
        out.append(Path([rounded(pts)], stroke=colour, width=width, dashes=dashes, cap=bond.cap,
                        join='round'))

    teleports = []
    for k in sorted(links, key=lambda k: abs(plane[k.a][0] - plane[k.b][0])):
        a, b = (k.a, k.b) if plane[k.a][0] <= plane[k.b][0] else (k.b, k.a)
        (xa, ya), (xb, yb) = plane[a], plane[b]
        if k.kind == 'head_to_tail':
            lo = min(ys) - BAR
            (ua, va), (ub, vb) = _exit(mol, plane, backbone, a, keys), _exit(mol, plane, backbone, b, keys)
            draw(clip([(xa, ya), (ua, va), (ua, lo), (ub, lo), (ub, vb), (xb, yb)], a, b))
        elif k.kind == 'bracket':
            s = k.side
            ua, ub = _headings(mol, plane, a, s, xb, keys), _headings(mol, plane, b, s, xa, keys)
            stubs = [(xa + ua[0][0] * L, ya + ua[0][1] * L), (xb + ub[0][0] * L, yb + ub[0][1] * L)]
            la, lb = stubs[:1], stubs[1:]
            x0, x1 = min(xa, la[0][0]), max(xb, lb[0][0])
            for _ in range(2):                            # a leg widens the span it clears
                level = max([s * y for x, y in plane.values() if x0 - .1 <= x <= x1 + .1] +
                            [s * y for x, y in stubs + la[:-1] + lb[:-1]]) + BAR
                for lo, hi, s2, l2 in levels:
                    if s2 == s and x0 <= lo + 1e-6 and hi <= x1 + 1e-6:
                        level = max(level, l2 + NEST)
                la, lb = _leg(mol, plane, a, ua, s * level, ink), _leg(mol, plane, b, ub, s * level, ink)
                x0, x1 = min(x0, la[-1][0]), max(x1, lb[-1][0])
            levels.append((x0, x1, s, level))
            pts = clip([(xa, ya), *la, *lb[::-1], (xb, yb)], a, b)
            out.append(Path([rounded(pts, BEND)], stroke=colour, width=width, cap=bond.cap, join='round'))
        else:
            teleports.append((min(xa, xb), a, b))
    for _, a, b in sorted(teleports):
        number += 1
        for q in (a, b):
            x, y = plane[q]
            d = _stub_direction(mol, plane, q, keys)
            ex, ey = x + d[0] * STUB, y + d[1] * STUB
            cx, cy = ex + d[0] * BADGE, ey + d[1] * BADGE
            draw(clip([(x, y), (ex, ey)], q, None), DASHES)
            out.append(Path([circle(cx, cy, BADGE)], stroke=colour, width=width * .75, fill=background))
            out.append(Text([TextRun(str(number), size=BADGE_TEXT)], x=cx, y=cy - BADGE_TEXT / 3, anchor='middle',
                            fill=colour))
    return out, number
