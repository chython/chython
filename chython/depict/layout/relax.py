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
"""The ring-shape post-steps run on a projected plane after the controller's search, and the plane predicates its
score and edits share.

| pass | rule |
| --- | --- |
| `regularize` | a ring of 4, 5, 7 or 8 atoms, no chord through it, becomes its regular polygon |
| | the fused bond, or the atom carrying the largest branch, stays; a moved atom's branch turns with its exit |
| | unfused and refused so: the unit bond carrying the largest branches stays, then the second-heaviest atom |
| | fused on more bonds: the largest fused piece holds, each other turns rigidly with its bond or its atom's exit |
| | fused along a path of three or more atoms, a peri ring: the path holds, the rest an arc of unit bonds |
| | an atom's two branches go onto the thirds of its outside gap, turned off them by up to 30 degrees for room |
| | a macrocycle atom's bond back into the macrocycle stretches; its one other branch bisects the widest gap |
| `wheels` | an atom of four ring bonds, one face a ring of 6 to 8 atoms, at least one pentagon and the rest pentagons \
or gaps, each pentagon fused at it: that ring holds and the rest of the turn splits evenly, a gap narrowing to 70 \
degrees or its own before a pentagon stays squashed; each pentagon's far atoms an arc of unit bonds, a branch on its \
atom's open gap |
| `fair` | a spike of a ring of 9 or more atoms, a turn of 110 degrees or sharper: the run of atoms turning its way, \
between the nearest turns the other way or fused atoms, becomes a circular arc of equal bonds |
| `seat` | a small ring fused into a ring of 9 or more atoms by one path: it turns and shifts rigidly, the two \
macrocycle bonds held at L, to the pose that evens the angles at the four atoms around those bonds |

A step is kept only when no defect count grows on the atoms it moved: `ovl61`, `crs` and `aob` of the scorecard. No
step may break a cis/trans configuration or bend an atom of `linear_atoms`.
"""
from collections import deque
from math import atan2, cos, hypot, pi, radians, sin, sqrt

from ...core.wedge import cis_trans_frame
from .lattice import L
from .quality import _inside, linear_atoms


__all__ = ['regularize', 'wheels', 'fair', 'seat', 'cross_atoms', 'count_defects', 'defect_atoms']

CLOSE = .61 * L         # the scorecard's `ovl61`
ALIGNS = 0, .5, 1       # fractions of the way a regularized ring's one-bond branch turns onto its exit
SHORT = .45 * L         # a bonded pair closer than this is close: a planar cage's .53 L spoke is not
ON_BOND = .25 * L       # the scorecard's `aob`
SPIKE = 110 * pi / 180     # a turn `fair` smooths: a ring's interior angle of 70 degrees or less
STRAIGHT = pi / 36         # a turn the other way below this still joins a `fair` run
FANS = (0., pi / 12, -pi / 12, pi / 6, -pi / 6, None)  # a two-branch ring atom's turns off the thirds; None: rigid
ROOM = .7 * L           # the least distance a repair leaves between a moved atom and one it is not bonded to
SEAT_TURNS = tuple(k * pi / 180 for k in range(-30, 31))    # `seat`'s rigid turns of a fused ring
SEAT_OFF = .05 * L                      # a macrocycle bond into a fused ring this far off L unseats it
SEAT_GAIN = .05                         # the least drop in `seat`'s angle cost, in squared radians


def _cross(p, q, a, b):
    def side(o, u, v):
        return (u[0] - o[0]) * (v[1] - o[1]) - (u[1] - o[1]) * (v[0] - o[0])
    return side(p, q, a) * side(p, q, b) < -1e-12 and side(a, b, p) * side(a, b, q) < -1e-12     # collinear is not


def _seg_dist(x, p, q):
    dx, dy = q[0] - p[0], q[1] - p[1]
    d = dx * dx + dy * dy
    t = 0. if not d else max(0., min(1., ((x[0] - p[0]) * dx + (x[1] - p[1]) * dy) / d))
    return hypot(p[0] + t * dx - x[0], p[1] + t * dy - x[1])


def count_defects(plane, atoms, bonds, moved=None):
    """`(close, crossing, on_bond)` among `atoms` and `bonds`; with `moved`, only those touching a moved atom. A
    bonded pair is close only within `SHORT`."""
    if moved is None:
        moved = set(atoms)
    link = {(a, b) for a, b in bonds} | {(b, a) for a, b in bonds}
    close = cross = on = 0
    mv = [(x, plane[x]) for x in atoms if x in moved]
    for i, (a, (xa, ya)) in enumerate(mv):                 # moved against unmoved, then each moved pair once
        for b in atoms:
            if b not in moved:
                xb, yb = plane[b]
                if abs(xa - xb) < CLOSE and abs(ya - yb) < CLOSE:
                    close += hypot(xa - xb, ya - yb) < (SHORT if (a, b) in link else CLOSE)
        for b, (xb, yb) in mv[i + 1:]:
            if abs(xa - xb) < CLOSE and abs(ya - yb) < CLOSE:
                close += hypot(xa - xb, ya - yb) < (SHORT if (a, b) in link else CLOSE)
    box = {}                                                # a bond's bounding box, widened by `ON_BOND`
    for a, b in bonds:
        (xa, ya), (xb, yb) = plane[a], plane[b]
        box[a, b] = min(xa, xb) - ON_BOND, max(xa, xb) + ON_BOND, min(ya, yb) - ON_BOND, max(ya, yb) + ON_BOND
    touched = [(n, m) for n, m in bonds if n in moved or m in moved]
    for a, b in bonds:
        x0, x1, y0, y1 = box[a, b]
        for c, e in touched:
            u0, u1, v0, v1 = box[c, e]
            if u0 <= x1 and x0 <= u1 and v0 <= y1 and y0 <= v1 and (c, e) != (a, b) and \
                    len({a, b, c, e}) == 4 and _cross(plane[a], plane[b], plane[c], plane[e]):
                cross += 1
    every = [(x, plane[x]) for x in atoms]
    for a, b in bonds:
        x0, x1, y0, y1 = box[a, b]
        for x, p in every if a in moved or b in moved else mv:     # an untouched bond needs a moved atom
            if x0 <= p[0] <= x1 and y0 <= p[1] <= y1 and x != a and x != b and \
                    _seg_dist(p, plane[a], plane[b]) < ON_BOND:
                on += 1
    return close, cross, on


def defect_atoms(plane, atoms, bonds):
    """Atoms in a close pair, a bonded one within `SHORT`, on a crossing bond, or on a bond, with the ends of that
    bond."""
    out = set()
    link = {(a, b) for a, b in bonds} | {(b, a) for a, b in bonds}
    every = [(x, plane[x]) for x in atoms]
    for i, (a, (xa, ya)) in enumerate(every):
        for b, (xb, yb) in every[i + 1:]:
            if abs(xa - xb) < CLOSE and abs(ya - yb) < CLOSE and \
                    hypot(xa - xb, ya - yb) < (SHORT if (a, b) in link else CLOSE):
                out.update((a, b))
    box = []
    for a, b in bonds:
        (xa, ya), (xb, yb) = plane[a], plane[b]
        box.append((a, b, min(xa, xb) - ON_BOND, max(xa, xb) + ON_BOND, min(ya, yb) - ON_BOND, max(ya, yb) + ON_BOND))
    for i, (a, b, x0, x1, y0, y1) in enumerate(box):
        for c, e, u0, u1, v0, v1 in box[i + 1:]:
            if u0 <= x1 and x0 <= u1 and v0 <= y1 and y0 <= v1 and len({a, b, c, e}) == 4 and \
                    _cross(plane[a], plane[b], plane[c], plane[e]):
                out.update((a, b, c, e))
        for x, p in every:
            if x0 <= p[0] <= x1 and y0 <= p[1] <= y1 and x != a and x != b and \
                    _seg_dist(p, plane[a], plane[b]) < ON_BOND:
                out.update((a, b, x))
    return out


def _bonds(mol, plane, rank):
    """The bonds within `plane`, each from its lower-ranked atom, in rank order: a float sum or a moved side is
    then a function of the molecule, not of its atom ids."""
    out = [(b.n, b.m) if rank[b.n] < rank[b.m] else (b.m, b.n) for b in mol.bonds() if b.n in plane]
    return sorted(out, key=lambda b: (rank[b[0]], rank[b[1]]))


def _frames(mol, plane):
    """`[(unit, frame, frame atoms)]` of the cis/trans units with a stated parity whose frame lies in `plane`."""
    out = []
    for u in mol.stereo_units():
        if u['kind'] != 1 or not u['parity']:
            continue
        f = cis_trans_frame(mol, u)
        if f is not None and all(x in plane for x in f):
            ends = set(f) | set(mol.neighbors_of(f[1])) | set(mol.neighbors_of(f[2]))
            out.append((u, f, ends))
    return out


def _ct_broken(mol, plane, frames, moved):
    """Units among `frames` touching `moved` that `plane` draws against their stated parity: the named pair on
    the wrong sides of the bond, an end's other substituent on its named one's side, or any substituent on the
    bond's line, which reads the same whichever substituent the parity names."""
    bad = 0
    for u, (near, anchor, partner, far), ends in frames:
        if not ends & moved:
            continue
        (ax, ay), (bx, by) = plane[anchor], plane[partner]
        dx, dy = bx - ax, by - ay
        eps = 1e-6 * (dx * dx + dy * dy)

        def side(x):
            return dx * (plane[x][1] - ay) - dy * (plane[x][0] - ax)
        sn, sf = side(near), side(far)
        ok = abs(sn) > eps and abs(sf) > eps and (sn * sf > 0) == (u['parity'] == 2)
        for end, other, ref in ((anchor, partner, sn), (partner, anchor, sf)):
            for x in mol.neighbors_of(end):
                if x != other and x != near and x != far and x in plane:
                    t = side(x)
                    ok = ok and abs(t) > eps and t * ref < 0
        bad += not ok
    return bad


def _bent(mol, plane, linear, moved):
    """Atoms of `linear` next to or in `moved` that `plane` draws more than 5 degrees off straight."""
    bad = 0
    for x in linear:
        a, b = mol.neighbors_of(x)
        if x not in moved and a not in moved and b not in moved or a not in plane or b not in plane:
            continue
        (px, py), (ax, ay), (bx, by) = plane[x], plane[a], plane[b]
        u = hypot(ax - px, ay - py) * hypot(bx - px, by - py)
        bad += not u or ((ax - px) * (bx - px) + (ay - py) * (by - py)) / u > -cos(pi / 36)
    return bad


def _near(planes, atoms, bonds, moved, margin=2 * L):
    """The `atoms` and `bonds` within `margin` of the moved atoms, or of a bond holding one, in any of `planes`."""
    ends = set(moved)
    for n, m in bonds:
        if n in moved or m in moved:
            ends.add(n)
            ends.add(m)
    xs = [p[x][0] for p in planes for x in ends]
    ys = [p[x][1] for p in planes for x in ends]
    x0, x1, y0, y1 = min(xs) - margin, max(xs) + margin, min(ys) - margin, max(ys) + margin
    if len(planes) == 1:
        p = planes[0]
        near = [a for a in atoms if x0 <= p[a][0] <= x1 and y0 <= p[a][1] <= y1]
        keep = []
        for n, m in bonds:
            (ax, ay), (bx, by) = p[n], p[m]
            if (ax <= x1 or bx <= x1) and (ax >= x0 or bx >= x0) and (ay <= y1 or by <= y1) and \
                    (ay >= y0 or by >= y0):
                keep.append((n, m))
        return near, keep
    near = [a for a in atoms if any(x0 <= p[a][0] <= x1 and y0 <= p[a][1] <= y1 for p in planes)]
    keep = []
    for n, m in bonds:
        bx = [p[z][0] for p in planes for z in (n, m)]
        by = [p[z][1] for p in planes for z in (n, m)]
        if min(bx) <= x1 and max(bx) >= x0 and min(by) <= y1 and max(by) >= y0:
            keep.append((n, m))
    return near, keep


def _rotate(p, o, c, s):
    x, y = p[0] - o[0], p[1] - o[1]
    return o[0] + c * x - s * y, o[1] + s * x + c * y


def _widest(o, pts):
    """Bisector angle of the widest gap between the directions from `o` to `pts`."""
    angs = sorted(atan2(p[1] - o[1], p[0] - o[0]) for p in pts)
    gaps = [((angs[(i + 1) % len(angs)] - a) % (2 * pi) or 2 * pi, a) for i, a in enumerate(angs)]
    g, a = max(gaps)
    return a + g / 2


def _exit(plane, ring, i):
    """Outward bisector angle of the ring at its `i`-th atom: away from the ring's inside at a reflex atom too, and
    the chord's outer normal at a straight one."""
    n = len(ring)
    p, a, b = plane[ring[i]], plane[ring[i - 1]], plane[ring[(i + 1) % n]]
    ux, uy = p[0] - a[0], p[1] - a[1]
    vx, vy = p[0] - b[0], p[1] - b[1]
    du, dv = hypot(ux, uy) or 1., hypot(vx, vy) or 1.
    sx, sy = ux / du + vx / dv, uy / du + vy / dv
    wind = sum(plane[ring[k - 1]][0] * plane[ring[k]][1] - plane[ring[k]][0] * plane[ring[k - 1]][1]
               for k in range(n))
    nx, ny = (b[1] - a[1], a[0] - b[0]) if wind > 0 else (a[1] - b[1], b[0] - a[0])
    if hypot(sx, sy) < .1:
        return atan2(ny, nx)
    if sx * nx + sy * ny < 0:
        return atan2(-sy, -sx)
    return atan2(sy, sx)


def _branch(mol, x, ring, elastic=False):
    """Atoms reached from `x` off the ring, or None when that reaches the ring again. An `elastic` atom keeps the
    neighbours that do not reach it, and the bonds to the rest stretch."""
    stop = set(ring)
    out = set()
    for y in mol.neighbors_of(x):
        if y in stop:
            continue
        seen = {x, y}
        queue = deque([y])
        while queue:
            a = queue.popleft()
            for b in mol.neighbors_of(a):
                if b in stop and b != x:
                    if not elastic:
                        return None
                    seen = None
                    break
                if b not in seen:
                    seen.add(b)
                    queue.append(b)
            if seen is None:
                break
        if seen is not None:
            seen.discard(x)
            out |= seen
    return out


def _polygon(plane, ring, fixed, anchor):
    """Regular positions `{atom: point}` of `ring`; `fixed` is a pair of ring neighbours or None."""
    n = len(ring)
    wind = 0.
    for i in range(n):
        (x1, y1), (x2, y2) = plane[ring[i]], plane[ring[(i + 1) % n]]
        wind += x1 * y2 - x2 * y1
    turn = 2 * pi / n * (1 if wind > 0 else -1)
    if fixed is not None:
        i = ring.index(fixed[0])
        if ring[(i + 1) % n] != fixed[1]:
            i = ring.index(fixed[1])
        a, b = plane[ring[i]], plane[ring[(i + 1) % n]]
        side = hypot(b[0] - a[0], b[1] - a[1])
        ang = atan2(b[1] - a[1], b[0] - a[0])
        out = {ring[i]: a, ring[(i + 1) % n]: b}
        p = b
        for k in range(2, n):
            ang += turn
            p = (p[0] + side * cos(ang), p[1] + side * sin(ang))
            out[ring[(i + k) % n]] = p
        return out
    # no fused bond: the anchor stays, the centre on its inward bisector, the side one bond
    radius = L / (2 * sin(pi / n))
    i = ring.index(anchor)
    e = _exit(plane, ring, i)
    a = plane[anchor]
    c = (a[0] - radius * cos(e), a[1] - radius * sin(e))
    return {ring[(i + k) % n]: (c[0] + radius * cos(e + k * turn), c[1] + radius * sin(e + k * turn))
            for k in range(n)}


def _pieces(mol, plane, ring):
    """`[(ring atoms touched, atoms)]` for each connected piece of `plane` off `ring`."""
    stop = set(ring)
    seen = set()
    out = []
    for x in ring:
        for y in mol.neighbors_of(x):
            if y in stop or y in seen or y not in plane:
                continue
            piece, touch = {y}, set()
            queue = deque([y])
            while queue:
                a = queue.popleft()
                for b in mol.neighbors_of(a):
                    if b in stop:
                        touch.add(b)
                    elif b not in piece and b in plane:
                        piece.add(b)
                        queue.append(b)
            seen |= piece
            out.append((touch, piece))
    return out


def _length(plane, a, b):
    return hypot(plane[a][0] - plane[b][0], plane[a][1] - plane[b][1])


def _regular(plane, ring):
    """Whether `ring` is its regular polygon of side L."""
    n = len(ring)
    return all(abs(_length(plane, a, b) - L) < 1e-6 for a, b in zip(ring, ring[1:] + ring[:1])) and \
        all(abs(t - 2 * pi / n) < 1e-6 for t in _turns(plane, ring))


def _hub(mol, plane, ring, rank):
    """`{atom: point}` of `ring` regular and of each piece off it turned rigidly with it: the largest piece fused on
    a ring bond holds, a piece fused on another bond follows that bond, a branch its atom's exit. None for a chord,
    a piece touching the ring anywhere else, or no fused piece."""
    n = len(ring)
    pos = {x: i for i, x in enumerate(ring)}
    for x in ring:
        if any(y in pos and (pos[y] - pos[x]) % n not in (1, n - 1) for y in mol.neighbors_of(x)):
            return None
    pieces = _pieces(mol, plane, ring)
    for touch, _ in pieces:
        if len(touch) > 2 or len(touch) == 2 and (pos[min(touch, key=pos.get)] - pos[max(touch, key=pos.get)]) % n \
                not in (1, n - 1):
            return None
    fused = [(t, p) for t, p in pieces if len(t) == 2]
    if not fused:
        return None
    # a bond of one L first: the lattice draws a pentagon's fifth bond long
    touch, big = max(fused, key=lambda e: (abs(_length(plane, *e[0]) - L) < 1e-6, len(e[1]),
                                           -min(rank[z] for z in e[1])))
    a, b = touch
    target = _polygon(plane, ring, (a, b) if ring[(pos[a] + 1) % n] == b else (b, a), None)
    out = dict(target)
    for touch, piece in pieces:
        if piece is big:
            continue
        if len(touch) == 2:                                 # about the bond's midpoint: an atom order is no choice
            a, b = touch
            (o0, o1), (t0, t1) = (plane[a], plane[b]), (target[a], target[b])
            turn = atan2(t1[1] - t0[1], t1[0] - t0[0]) - atan2(o1[1] - o0[1], o1[0] - o0[0])
            o0, t0 = ((o0[0] + o1[0]) / 2, (o0[1] + o1[1]) / 2), ((t0[0] + t1[0]) / 2, (t0[1] + t1[1]) / 2)
        else:
            a, = touch
            o0, t0 = plane[a], target[a]
            turn = _exit(target, ring, pos[a]) - _exit(plane, ring, pos[a])
        c, s = cos(turn), sin(turn)
        for y in piece:
            q = _rotate(plane[y], o0, c, s)
            out[y] = (q[0] - o0[0] + t0[0], q[1] - o0[1] + t0[1])
    return out


def _peri(plane, ring, shared):
    """`{atom: point}` of `ring` with its `shared` atoms, one path along it, held and the rest an arc of unit bonds
    on their side; None when the free atoms are not one run or the arc cannot span the path's ends."""
    n = len(ring)
    free = [i for i, x in enumerate(ring) if x not in shared]
    if not free:
        return None
    i0 = next((i for i in free if ring[i - 1] in shared), None)
    if i0 is None or any(ring[(i0 + k) % n] in shared for k in range(len(free))):
        return None
    run = [ring[(i0 + k) % n] for k in range(len(free))]
    p, q = plane[ring[i0 - 1]], plane[ring[(i0 + len(free)) % n]]
    cx = sum(plane[x][0] for x in run) / len(run)
    cy = sum(plane[x][1] for x in run) / len(run)
    best = None
    for left in (True, False):
        pts = _arc(p, q, len(run) + 1, L, left)
        if pts is None:
            continue
        d = sum(hypot(u - cx, v - cy) for u, v in pts)
        if best is None or d < best[0]:
            best = d, pts
    if best is None:
        return None
    out = {x: plane[x] for x in ring}
    out.update(zip(run, best[1]))
    return out


def _thirds(mol, plane, ring, j, roots, shift):
    """`{atom: point}` of the two branches at `roots` of ring atom `ring[j]`, each turned rigidly about it onto the
    thirds of the gap outside its ring bonds, in their order, then both by `shift`; empty when the branches meet."""
    x = ring[j]
    px, py = plane[x]
    sides = []
    for y in roots:
        seen = {x, y}
        queue = deque([y])
        while queue:
            for b in mol.neighbors_of(queue.popleft()):
                if b not in seen:
                    seen.add(b)
                    queue.append(b)
        seen.discard(x)
        sides.append(seen)
    if sides[0] & sides[1]:
        return {}

    def angle(z):
        return atan2(plane[z][1] - py, plane[z][0] - px)
    a, b = angle(ring[j - 1]), angle(ring[(j + 1) % len(ring)])
    if (b - a) % (2 * pi) < pi:
        a, b = b, a                                         # the outside, the wider arc, runs from a to b
    g = (b - a) % (2 * pi)
    out = {}
    for k, (y, side) in enumerate(sorted(zip(roots, sides), key=lambda e: (angle(e[0]) - a) % (2 * pi)), 1):
        d = a + k * g / 3 + shift - angle(y)
        c, s = cos(d), sin(d)
        for z in side:
            out[z] = _rotate(plane[z], (px, py), c, s)
    return out


def regularize(mol, plane, rings, rank):
    """Regularize in place each ring of `rings` (atom tuples in cyclic order) that allows it; returns the count. A
    ring left irregular tries once more after the rest: a dioxolane drawn on a furanose's long bond is regular once
    the furanose pulled that bond to L."""
    done = _regularize(mol, plane, rings, rank)
    again = [r for r in rings if len(r) != 6 and not _regular(plane, r)]
    if again and done:
        done += _regularize(mol, plane, again, rank, rings)
    return done


def _regularize(mol, plane, todo, rank, rings=None):
    """`regularize`'s one round over the rings of `todo`, fused as in `rings` (`todo` itself when None)."""
    if rings is None:
        rings = todo
    atoms = sorted(plane, key=rank.__getitem__)
    bonds = _bonds(mol, plane, rank)
    frames = _frames(mol, plane)
    linear = linear_atoms(mol) & set(plane)
    edges = {}                                              # ring bond -> the small rings holding it
    large = set()                                           # atoms of a macrocycle: a branch may close through it
    for r in rings:
        if len(r) > 8:
            large.update(r)
            continue
        for a, b in zip(r, r[1:] + r[:1]):
            edges.setdefault(frozenset((a, b)), set()).add(r)
    done = 0
    # a hexagon last and only when not regular: the lattice draws it so, save one fused to a pentagon's long bond
    for ring in sorted(todo, key=lambda r: (len(r) == 6, sorted(rank[x] for x in r))):
        n = len(ring)
        if n not in (4, 5, 6, 7, 8) or n == 6 and _regular(plane, ring):
            continue
        # an atom of a bond another ring holds; a spiro atom is not, its other ring moving as its branch
        shared = [x for i, x in enumerate(ring) if len(edges[frozenset((ring[i - 1], x))]) > 1 or
                  len(edges[frozenset((x, ring[(i + 1) % n]))]) > 1]
        fixed = None
        if len(shared) == 2:
            i = ring.index(shared[0])
            if shared[1] in (ring[i - 1], ring[(i + 1) % n]):
                fixed = tuple(shared)
        peri = None
        if len(shared) > 2:                                 # one path of fused bonds, each ortho: else a hub
            rs = set(ring)
            held = [others for a, b in zip(ring, ring[1:] + ring[:1])
                    if len(others := [r for r in edges[frozenset((a, b))] if r is not ring])]
            if len(held) == len(shared) - 1 and all(len(rs.intersection(r)) == 2 for o in held for r in o):
                peri = _peri(plane, ring, shared)
        if len(shared) > 1 and fixed is None and peri is None:  # fused on more than one bond: each fused piece turns
            if (hub := _hub(mol, plane, ring, rank)) is None:
                continue
            branches = None
        else:
            hub = None
            branches = {x: _branch(mol, x, ring, x in large) for x in ring if x not in shared}
            if any(v is None for v in branches.values()):
                continue
        if branches is None:
            poses = [(fixed, None)]
        elif len(shared) == 1:
            poses = [(fixed, shared[0])]
        elif shared:
            poses = [(fixed, None)]
        else:                                               # the heaviest atom held, then the heaviest unit bond
            heavy = sorted(ring, key=lambda x: (-len(branches[x]), rank[x]))
            pairs = sorted(((a, b) for a, b in zip(ring, ring[1:] + ring[:1]) if abs(_length(plane, a, b) - L) < 1e-6),
                           key=lambda e: (-len(branches[e[0]]) - len(branches[e[1]]), heavy.index(e[0])))
            poses = [(None, heavy[0])] + [(e, None) for e in pairs[:1]] + [(None, x) for x in heavy[1:2]]
        start = dict(plane)
        ctx = atoms, bonds, frames, linear
        # each branch turned with its exit, then a one-bond branch turned toward it: an N-alkyl the lattice squeezed
        if any(_posed(mol, plane, ring, shared, fixed, anchor, peri or hub, branches, large, start, ctx, align)
               for align in ALIGNS if branches or not align for fixed, anchor in poses):
            done += 1
    return done


def _posed(mol, plane, ring, shared, fixed, anchor, peri, branches, large, old, ctx, align=0):
    """Set `ring` regular in `plane` from `old`, holding the bond `fixed` or the atom `anchor`, each branch turned with
    its atom's exit, a one-bond branch that `align` fraction of its way onto that exit, or to the points `peri` with no
    `branches`; True when that grows no defect, cis/trans break or bend, else `plane` is set back to `old`."""
    atoms, bonds, frames, linear = ctx
    fans = []                                               # ring atoms carrying two branches
    if branches is None:
        plane.update(peri)
    else:
        target = peri or _polygon(old, ring, fixed, anchor)
        held = set(shared) | set(fixed or ()) | {anchor}
        for j, x in enumerate(ring):
            if x in held:
                continue
            e0 = _exit(old, ring, j)
            plane[x] = target[x]
            if branches[x]:
                e1 = _exit(target, ring, j)
                kept = [y for y in mol.neighbors_of(x) if y in branches[x]]
                if x in large and len(kept) == 1:           # the branch bisects the widest gap left by the chain
                    e0 = atan2(old[kept[0]][1] - old[x][1], old[kept[0]][0] - old[x][0])
                    e1 = _widest(target[x], [target.get(z, plane[z]) for z in mol.neighbors_of(x) if z != kept[0]])
                elif align and len(kept) == 1:
                    e = atan2(old[kept[0]][1] - old[x][1], old[kept[0]][0] - old[x][0]) - e0
                    e0 += align * atan2(sin(e), cos(e))
                elif x not in large and len(kept) == 2:
                    fans.append((j, kept))
                c, s = cos(e1 - e0), sin(e1 - e0)
                o, t = old[x], target[x]
                for y in branches[x]:
                    q = _rotate(old[y], o, c, s)
                    plane[y] = (q[0] - o[0] + t[0], q[1] - o[1] + t[1])
    rigid = dict(plane)
    # two branches on the thirds, then turned off them toward a clear side; last as the lattice drew them
    for shift in (FANS if fans else (None,)):
        plane.update(rigid)
        if shift is not None:
            for j, kept in fans:
                plane.update(_thirds(mol, rigid, ring, j, kept, shift))
        moved = {x for x in plane if plane[x] != old[x]}
        if not moved:
            return False
        near, local = _near((old, plane), atoms, bonds, moved)
        if all(a <= b for a, b in zip(count_defects(plane, near, local, moved),
                                      count_defects(old, near, local, moved))) and \
                _ct_broken(mol, plane, frames, moved) <= _ct_broken(mol, old, frames, moved) and \
                _bent(mol, plane, linear, moved) <= _bent(mol, old, linear, moved):
            return True
    plane.update(old)
    return False


def wheels(mol, plane, rings, rank):
    """Share out in place the angle at an atom of four ring bonds, one face a ring of 6 to 8 atoms and the rest
    pentagons or gaps: that ring holds, its two bonds stay, and the rest of the turn splits evenly over the other
    faces, a pentagon taking no more than its 108 degrees. Each pentagon's two far atoms become an arc of unit bonds
    between its moved neighbours, a branch turning with its atom's exit; kept only when no defect count grows. Returns
    the atoms so drawn."""
    small = [r for r in rings if len(r) <= 8]
    atoms = sorted(plane, key=rank.__getitem__)
    bonds = _bonds(mol, plane, rank)
    frames = _frames(mol, plane)
    linear = linear_atoms(mol) & set(plane)
    done = 0
    for x in atoms:
        nb = mol.neighbors_of(x)
        if len(nb) != 4 or not all(y in plane for y in nb):
            continue
        px, py = plane[x]
        order = sorted(nb, key=lambda y: atan2(plane[y][1] - py, plane[y][0] - px))
        faces = []
        for i, u in enumerate(order):
            v = order[(i + 1) % 4]
            faces.append((u, v, next((r for r in small if _face(r, u, x, v)), None)))
        held = [k for k, (_, _, r) in enumerate(faces) if r is not None and len(r) != 5]
        pent = [r for _, _, r in faces if r is not None and len(r) == 5]
        if len(held) != 1 or not pent or any(not any(r and y in r for _, _, r in faces) for y in nb) or \
                any(all(sum(r is not None and y in r for _, _, r in faces) < 2 for y in nb if y in q) for q in pent):
            continue                                        # a spiro pentagon turns in `spiros`
        faces = faces[held[0]:] + faces[:held[0]]
        system = {z for _, _, r in faces if r for z in r}
        movers = {v for _, v, _ in faces[1:3]} | {z for r in pent for z in r if z != x and z not in nb}
        if any(z in r for r in small if not any(r is f for _, _, f in faces) for z in movers):
            continue                                        # a moved atom held by a ring off the wheel
        u0, v0, hex_ = faces[0]
        turn = (atan2(plane[v0][1] - py, plane[v0][0] - px) - atan2(plane[u0][1] - py, plane[u0][0] - px)) % (2 * pi)
        rest = 2 * pi - turn
        gaps = sum(r is None for _, _, r in faces[1:])
        now = [(atan2(plane[v][1] - py, plane[v][0] - px) - atan2(plane[u][1] - py, plane[u][0] - px)) % (2 * pi)
               for u, v, _ in faces[1:]]
        cap = 3 * pi / 5
        wants = []
        for g in (rest / 3, radians(70), *(w for w, (_, _, r) in zip(now, faces[1:]) if r is None)):
            if not gaps:                                    # a gap narrows, down to its own, before a pentagon
                p_ = rest / 3
            else:
                p_ = min(cap, (rest - g * gaps) / len(pent))
                g = (rest - p_ * len(pent)) / gaps
            want = [p_ if r else g for _, _, r in faces[1:]]
            if want not in wants and any(abs(a - b) >= radians(5) for a, b in zip(want, now)):
                wants.append(want)
        for want in wants:
            new = {}
            ang = atan2(plane[v0][1] - py, plane[v0][0] - px)
            for (u, v, _), w in zip(faces[1:3], want):
                ang += w
                d = hypot(plane[v][0] - px, plane[v][1] - py) or L
                new[v] = (px + d * cos(ang), py + d * sin(ang))
            ok = True
            for u, v, r in faces[1:]:
                if r is None or len(r) != 5:
                    continue
                i = r.index(x)
                step = 1 if r[(i + 1) % 5] == v else -1
                far = [r[(i + 2 * step) % 5], r[(i + 3 * step) % 5]]  # from v's side to u's
                a, b = new.get(v, plane[v]), new.get(u, plane[u])
                best = None
                for left in (True, False):
                    pts = _arc(a, b, 3, L, left)
                    if pts is not None:
                        d = sum(hypot(q[0] - px, q[1] - py) for q in pts)
                        if best is None or d > best[0]:
                            best = d, pts
                if best is None:
                    ok = False
                    break
                new.update(zip(far, best[1]))
            if not ok:
                continue
            old = dict(plane)
            plane.update(new)
            for z in new:                                       # a branch turns with its atom's exit
                ring_nb = [y for y in mol.neighbors_of(z) if y in system]
                rest_nb = [y for y in mol.neighbors_of(z) if y not in system]
                if not rest_nb:
                    continue
                side = set()
                for y in rest_nb:
                    part = _side(mol, z, y)
                    if part is None or part & system:
                        side = None
                        break
                    side |= part
                if side is None:
                    ok = False
                    break
                o = old[z]                                      # the branches' own heading: a straight atom has two
                e0 = atan2(sum((old[y][1] - o[1]) / (_length(old, z, y) or L) for y in rest_nb),
                           sum((old[y][0] - o[0]) / (_length(old, z, y) or L) for y in rest_nb))
                t = plane[z]                                    # the widest gap no ring of the wheel fills
                angs = sorted((atan2(plane[y][1] - t[1], plane[y][0] - t[0]), y) for y in ring_nb)
                gaps = [((angs[(i + 1) % len(angs)][0] - a) % (2 * pi) or 2 * pi, a) for i, (a, y) in enumerate(angs)
                        if len(angs) == 1 or not any(_face(r, y, z, angs[(i + 1) % len(angs)][1]) for r in small)]
                g, a = max(gaps) if gaps else (0., _widest(t, [plane[y] for y in ring_nb]))
                e1 = a + g / 2
                c, s_ = cos(e1 - e0), sin(e1 - e0)
                o, t = old[z], plane[z]
                for y in side:
                    q = _rotate(old[y], o, c, s_)
                    plane[y] = (q[0] - o[0] + t[0], q[1] - o[1] + t[1])
            moved = {z for z in plane if plane[z] != old[z]}
            if ok and moved:
                near, local = _near((old, plane), atoms, bonds, moved)
                if all(a <= b for a, b in zip(count_defects(plane, near, local, moved),
                                              count_defects(old, near, local, moved))) and \
                        _ct_broken(mol, plane, frames, moved) <= _ct_broken(mol, old, frames, moved) and \
                        _bent(mol, plane, linear, moved) <= _bent(mol, old, linear, moved):
                    done += 1
                    break
            plane.update(old)
    return done


def _face(ring, u, x, v):
    """Whether `ring` runs u-x-v."""
    n = len(ring)
    if x not in ring:
        return False
    i = ring.index(x)
    return {ring[i - 1], ring[(i + 1) % n]} == {u, v}


def _turns(plane, ring):
    """Signed turn in radians at each atom of `ring`, positive toward the ring's inside."""
    n = len(ring)
    wind = sum(plane[ring[i]][0] * plane[ring[(i + 1) % n]][1] - plane[ring[(i + 1) % n]][0] * plane[ring[i]][1]
               for i in range(n))
    out = []
    for i in range(n):
        (ax, ay), (bx, by), (cx, cy) = plane[ring[i - 1]], plane[ring[i]], plane[ring[(i + 1) % n]]
        ux, uy, vx, vy = bx - ax, by - ay, cx - bx, cy - by
        t = atan2(ux * vy - uy * vx, ux * vx + uy * vy)
        out.append(t if wind > 0 else -t)
    return out


def _arc(p, q, k, side, left):
    """`k - 1` points of the circular arc of `k` equal chords of length `side` from `p` to `q`, bulging to the left
    of `p -> q` when `left`; None when the chords cannot span it."""
    span = hypot(q[0] - p[0], q[1] - p[1])
    if not 0 < span < k * side * .999:
        return None
    lo, hi = 0., 2 * pi / k                             # sin(k f / 2) / sin(f / 2) falls from k to 0 on this range
    for _ in range(60):
        f = (lo + hi) / 2
        if sin(k * f / 2) / sin(f / 2) > span / side:
            lo = f
        else:
            hi = f
    f = (lo + hi) / 2
    s = -1 if left else 1                               # a leftward bulge turns right along the chain
    ang = atan2(q[1] - p[1], q[0] - p[0]) - s * (k - 1) * f / 2
    out, x = [], p
    for _ in range(k - 1):
        x = (x[0] + side * cos(ang), x[1] + side * sin(ang))
        out.append(x)
        ang += s * f
    return out


def fair(mol, plane, rings, rank):
    """Fair in place each spike of a ring of `rings` of 9 or more atoms; returns the count."""
    atoms = sorted(plane, key=rank.__getitem__)
    bonds = _bonds(mol, plane, rank)
    frames = _frames(mol, plane)
    linear = linear_atoms(mol) & set(plane)
    member = {}
    for r in rings:
        for x in r:
            member[x] = member.get(x, 0) + 1
    done = 0
    for ring in sorted(rings, key=lambda r: sorted(rank[x] for x in r)):
        n = len(ring)
        if n < 9:
            continue
        turns = _turns(plane, ring)
        spikes = [i for i in range(n) if abs(turns[i]) >= SPIKE and member[ring[i]] == 1]
        for i in sorted(spikes, key=lambda i: rank[ring[i]]):
            turns = _turns(plane, ring)
            if abs(turns[i]) < SPIKE:
                continue
            sign = 1 if turns[i] > 0 else -1

            def free(j):
                return sign * turns[j % n] > -STRAIGHT and member[ring[j % n]] == 1
            a, b = i, i
            while free(a - 1) and b - a < n:
                a -= 1
            while free(b + 1) and b - a < n:
                b += 1
            if b - a >= n - 2:
                continue
            run = [ring[j % n] for j in range(a, b + 1)]
            branches = {x: _branch(mol, x, ring) for x in run}
            if any(v is None for v in branches.values()):
                continue
            p, q = plane[ring[(a - 1) % n]], plane[ring[(b + 1) % n]]
            chain = [ring[j % n] for j in range(a - 1, b + 2)]
            side = sum(hypot(plane[x][0] - plane[y][0], plane[x][1] - plane[y][1])
                       for x, y in zip(chain, chain[1:])) / (len(chain) - 1)
            tip = plane[ring[i]]                    # the arc bulges to the spike's side of its chord
            left = (q[0] - p[0]) * (tip[1] - p[1]) - (q[1] - p[1]) * (tip[0] - p[0]) > 0
            target = _arc(p, q, len(chain) - 1, side, left)
            if target is None:
                continue
            old = dict(plane)
            plane.update(zip(run, target))
            for j, x in zip(range(a, b + 1), run):
                if branches[x]:
                    e1 = _exit(plane, ring, j % n)
                    e0 = _exit(old, ring, j % n)
                    c, s = cos(e1 - e0), sin(e1 - e0)
                    o, t = old[x], plane[x]
                    for y in branches[x]:
                        r = _rotate(old[y], o, c, s)
                        plane[y] = (r[0] - o[0] + t[0], r[1] - o[1] + t[1])
            moved = {x for x in plane if plane[x] != old[x]}
            if not moved:                           # the spike is already the arc
                continue
            near, local = _near((old, plane), atoms, bonds, moved)
            before = count_defects(old, near, local, moved)
            after = count_defects(plane, near, local, moved)
            if any(u > v for u, v in zip(after, before)) or \
                    _ct_broken(mol, plane, frames, moved) > _ct_broken(mol, old, frames, moved) or \
                    _bent(mol, plane, linear, moved) > _bent(mol, old, linear, moved):
                plane.update(old)
            else:
                done += 1
    return done


def _side(mol, a, b):
    """Atoms reached from `b` without the bond `a`-`b`, or None when that reaches `a`: the bond is in a ring."""
    seen = {b}
    stack = [b]
    while stack:
        x = stack.pop()
        for y in mol.neighbors_of(x):
            if y == a:
                if x != b:
                    return None
                continue
            if y not in seen:
                seen.add(y)
                stack.append(y)
    return seen


def _moved(plane, part, pivot, other, move):
    """`{atom: point}` of `part`: mirrored across the line `pivot`-`other` (angle None) or turned by the angle
    about `pivot`; then pushed along the bond to `stretch` times its length, or to L for a stretch of 0."""
    angle, stretch = move
    o, t = plane[pivot], plane[other]
    if angle is None:
        dx, dy = t[0] - o[0], t[1] - o[1]
        d = dx * dx + dy * dy or 1.
        out = {}
        for x in part:
            px, py = plane[x][0] - o[0], plane[x][1] - o[1]
            k = 2 * (px * dx + py * dy) / d
            out[x] = (o[0] + k * dx - px, o[1] + k * dy - py)
    else:
        c, s = cos(angle), sin(angle)
        out = {x: _rotate(plane[x], o, c, s) for x in part}
    if stretch != 1:                                        # `part` holds `other`, never `pivot`
        a = out[other]
        k = stretch - 1 if stretch else L / (hypot(a[0] - o[0], a[1] - o[1]) or L) - 1
        dx, dy = (a[0] - o[0]) * k, (a[1] - o[1]) * k
        out = {x: (p[0] + dx, p[1] + dy) for x, p in out.items()}
    return out


def cross_atoms(mol):
    """Acyclic atoms of four neighbours, drawn as a right-angled cross, a tetrahedral centre included."""
    return {x for x in mol if len(mol.neighbors_of(x)) == 4 and not mol.in_ring_of(x)}


def _kept(mol, plane, move, part, ctx):
    """The least room of `part` in `move`, or None when the move crowds it within `ROOM` (or closer than it was),
    grows a defect count, breaks a cis/trans or linear atom or puts an atom inside a ring."""
    atoms, bonds, frames, linear, rings = ctx
    room = _room(mol, move, atoms, part)
    if room < min(_room(mol, plane, atoms, part), ROOM) or \
            _ct_broken(mol, move, frames, part) > _ct_broken(mol, plane, frames, part) or \
            _bent(mol, move, linear, part) > _bent(mol, plane, linear, part) or \
            _inr(move, rings, part) > _inr(plane, rings, part):
        return None
    near, local = _near((plane, move), atoms, bonds, part)
    if any(a > b for a, b in zip(count_defects(move, near, local, part), count_defects(plane, near, local, part))):
        return None
    return room


def _context(mol, plane, rank):
    atoms = sorted(plane, key=rank.__getitem__)
    return atoms, _bonds(mol, plane, rank), _frames(mol, plane), linear_atoms(mol) & set(plane), \
        _cycles(mol, plane, rank)


def _room(mol, plane, atoms, moved):
    """The least distance from a `moved` atom to one of `atoms` it is not bonded to."""
    least = float('inf')
    pts = [(y, plane[y]) for y in atoms]
    for x in moved:
        px, py = plane[x]
        nb = mol.neighbors_of(x)
        for y, (qx, qy) in pts:
            d = (px - qx) ** 2 + (py - qy) ** 2
            if d < least and y != x and y not in nb:
                least = d
    return sqrt(least)


def _cycles(mol, plane, rank, size=8):
    """Every simple cycle of 3 to `size` atoms within `plane`, each from its lowest-ranked atom: unlike an SSSR,
    a function of the molecule, not of its atom ids."""
    ring = [x for x in plane if mol.in_ring_of(x)]
    out = []
    for s in sorted(ring, key=rank.__getitem__):
        stack = [(s, [s])]
        while stack:
            x, path = stack.pop()
            for y in mol.neighbors_of(x):
                if y == s and len(path) >= 3 and rank[path[1]] < rank[path[-1]]:
                    out.append(tuple(path))
                elif y in plane and rank[y] > rank[s] and y not in path and len(path) < size and mol.in_ring_of(y):
                    stack.append((y, path + [y]))
    return out


def _inr(plane, rings, part=None):
    """Atoms inside a ring polygon on the other side of `part`: the scorecard's `inr` across the moved bond; with no
    `part`, every atom inside a ring it is not in."""
    n = 0
    for r in rings:
        poly = [plane[x] for x in r]
        x0, x1 = min(p[0] for p in poly), max(p[0] for p in poly)
        y0, y1 = min(p[1] for p in poly), max(p[1] for p in poly)
        ring, moved = set(r), part is not None and r[0] in part
        n += sum(_inside(q, poly) for x, q in plane.items() if x0 <= q[0] <= x1 and y0 <= q[1] <= y1 and
                 (x not in ring if part is None else (x in part) != moved))     # the ring's box first
    return n


def _evenness(mol, plane, xs):
    """Sum over `xs` of the squared differences between each gap between two bonds and its atom's even share."""
    cost = 0.
    for x in xs:
        px, py = plane[x]
        a = sorted(atan2(plane[y][1] - py, plane[y][0] - px) for y in mol.neighbors_of(x))
        if len(a) < 2:
            continue
        even = 2 * pi / len(a) if len(a) > 2 else 2 * pi / 3
        gaps = [v - u for u, v in zip(a, a[1:])] + [2 * pi - a[-1] + a[0]]
        cost += sum((g - even) ** 2 for g in gaps if g < pi) if len(a) == 2 else sum((g - even) ** 2 for g in gaps)
    return cost


def _circles(p, q, r):
    """The two points at distance `r` from both `p` and `q`, or none."""
    dx, dy = q[0] - p[0], q[1] - p[1]
    d = hypot(dx, dy)
    if d < 1e-9 or d > 2 * r:
        return ()
    h = (r * r - d * d / 4) ** .5
    mx, my = p[0] + dx / 2, p[1] + dy / 2
    return (mx - dy / d * h, my + dx / d * h), (mx + dy / d * h, my - dx / d * h)


def seat(mol, plane, rings, rank):
    """Turn and shift rigidly, in place, each ring of `rings` of 3 to 8 atoms other than a hexagon, fused into a ring
    of 9 or more by one path, with all it carries, where a macrocycle bond into it is `SEAT_OFF` off L, so both come
    to L and the angles at the four atoms around them come nearest even: a dihydrofuran pinched into a cembranoid
    ring. A hexagon already sits on the lattice. `rings` are canonical, as `ring_blocks` gives them: an SSSR's start
    and direction follow atom ids. Returns the rings seated."""
    rings = [tuple(r) for r in rings if all(x in plane for x in r)]
    macros = [r for r in rings if len(r) > 8]
    if not macros:
        return 0
    ctx = None
    done = 0
    for small in sorted((r for r in rings if len(r) <= 8 and len(r) != 6), key=lambda r: sorted(rank[x] for x in r)):
        rs = set(small)
        for macro in macros:
            n = len(macro)
            ends = [macro[i] for i in range(n) if macro[i] in rs and macro[i - 1] not in rs] + \
                [macro[i] for i in range(n) if macro[i] in rs and macro[(i + 1) % n] not in rs]
            if len(ends) != 2 or len(rs & set(macro)) < 2:
                continue
            i, j = macro.index(ends[0]), macro.index(ends[1])
            a, pa, b, pb = ends[0], macro[i - 1], ends[1], macro[(j + 1) % n]
            if pa in rs or pb in rs:
                continue
            if abs(_length(plane, a, pa) - L) < SEAT_OFF and abs(_length(plane, b, pb) - L) < SEAT_OFF:
                continue
            part, todo = set(small), list(small)
            while todo:
                for y in mol.neighbors_of(todo.pop()):
                    if y not in part and y not in (pa, pb):
                        part.add(y)
                        todo.append(y)
            if pa in part or pb in part:
                continue
            if ctx is None:
                ctx = _context(mol, plane, rank)
            around = (a, b, pa, pb)
            cost = _evenness(mol, plane, around)
            ax, ay = plane[a]
            best = None
            for t in SEAT_TURNS:
                c, s = cos(t), sin(t)
                turned = {x: _rotate(plane[x], plane[a], c, s) for x in part}
                vx, vy = turned[b][0] - ax, turned[b][1] - ay
                for q in _circles(plane[pa], (plane[pb][0] - vx, plane[pb][1] - vy), L):
                    move = dict(plane)
                    for x in part:
                        move[x] = (turned[x][0] - ax + q[0], turned[x][1] - ay + q[1])
                    if hypot(q[0] - ax, q[1] - ay) > L:
                        continue
                    k = _evenness(mol, move, around)
                    if k < cost - SEAT_GAIN and (best is None or k < best[0] - 1e-9) and \
                            _kept(mol, plane, move, part, ctx) is not None:
                        best = k, move
            if best is not None:
                plane.update(best[1])
                done += 1
    return done
