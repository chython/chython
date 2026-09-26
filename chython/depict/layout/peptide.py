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
"""Peptide layout: the backbone on the 120-degree lattice, side chains alternating up and down.

| Stage | Rule |
| --- | --- |
| backbone | one bond of `L` per vertex; first bond at 30 deg, vertex k turns -60 odd, +60 even |
| steps | inverted vertices come in even runs, so the chain keeps its axis; one ending the chain may be odd |
| backbone ring | a five-ring on N-C-alpha is its lattice hexagon less the top or bottom vertex: a horizontal chord |
| bridge | the amide between two backbone rings is `BRIDGE` long, along its lattice direction |
| side chains | acyclic atoms walked natively at 120 deg; ring systems are engine tiles, cached by form |
| freedoms | reflection, torsion flips at depth 2/3, then a stretched root bond -- never a rotation |
| stretch removal | pair relaxation, then steps near a stretched anchor; never at the cost of a collision or a bracket |

`peptide_layout` takes a `tile(sub)` callable returning a plane at mean bond `L`, and names no engine.
The plane it returns is complete; cross-links are drawn by `depict.peptide` from the plane alone.
"""
from collections import deque
from math import atan2, cos, degrees, hypot, pi, radians, sin

from ...core.monomers import _bfs_path
from ...core.wedge import cis_trans_frame, cis_trans_parity


__all__ = ['peptide_layout']

L = .825
CLOSE = .8 * L          # atom to atom, or to a reserved slot
ON_BOND = .45 * L       # atom to a bond it does not belong to
STRETCHES = (0., .3, .6, 1., 1.5, 2., 3.)
BUDGET = 40             # re-layouts spent on removing stretch
BRIDGE = 2 * L          # the amide between two backbone rings
FLAT = 2.5              # side solver: cost per half-bond of backbone height; an inversion costs 1
WINDOW = 12             # side solver: states kept above the cheapest


# --- lattice ----------------------------------------------------------------------------------------------
def _walk(start, dirs, bridges=()):
    """`start`, then one bond per direction in degrees: `L`, or `BRIDGE` for a bond index in `bridges`."""
    out = [start]
    for k, d in enumerate(dirs):
        x, y = out[-1]
        s = BRIDGE if k in bridges else L
        out.append((x + s * cos(radians(d)), y + s * sin(radians(d))))
    return out


def _backbone_dirs(n, inv):
    """Bond directions of an n-atom chain: 30, -30, 30, ...; a vertex in `inv` turns the other way."""
    dirs, d = [], 30
    for k in range(1, n):
        dirs.append(d)
        if k < n - 1:
            t = -60 if k % 2 else 60
            d += -t if k in inv else t
    return dirs


def _default_side(k):
    """+1 up, -1 down: where vertex k's substituents point in the plain zigzag."""
    return 1 if k % 2 else -1


def _side(k, inv):
    return -_default_side(k) if k in inv else _default_side(k)


def _solve(n, pins, differ, same=(), pin_cost=10, bridges=()):
    """Inverted vertices: a union of even runs inside 1..n-2, the last one odd if it ends at n-2, by dynamic
    programming.

    Cost 1 per inverted vertex, `pin_cost` per violated pin `{k: must_invert}`, 100 per violated pair of
    `differ` `(k, k+1)` that must have `inv[k] != inv[k+1]`, 5 per violated pair of `same` that should have
    them equal, `FLAT` per half-bond of height, max less min, so a chain buys a step back down; a bond index
    in `bridges` rises twice.  States more than `WINDOW` above the cheapest are dropped.  States: run 0 not
    inverted, 1 odd, 2 even; bond direction; height above the lowest vertex; the height so far.
    """
    ks = range(1, n - 1)
    if not ks:
        return set()
    best = {(0, 30, 1, 1): (FLAT, None)}                  # path as (t, parent), newest first
    differ, same = set(differ), set(same)
    for k in ks:
        new, rise, want = {}, 2 if k in bridges else 1, pins.get(k)
        turn = -60 if k % 2 else 60
        split, keep = (k - 1, k) in differ, (k - 1, k) in same
        for (s, d, y, h), (c, path) in best.items():
            for t in (0, 1, 2):
                if s == 1 and t != 2 or s != 1 and t == 2:      # an odd run must continue, once
                    continue
                cc = c + (t != 0)
                if want is not None and want != (t != 0):
                    cc += pin_cost
                if split and (s != 0) == (t != 0):
                    cc += 100
                if keep and (s != 0) != (t != 0):
                    cc += 5
                dd = (d + (-turn if t else turn)) % 360
                yy = y + rise * _RISE[dd]
                hh = max(h, yy) - min(0, yy)
                cc += FLAT * (hh - h)
                key = (t, dd, max(yy, 0), hh)
                if key not in new or cc < new[key][0]:
                    new[key] = (cc, (t, path))
        low = min(c for c, _ in new.values())
        best = {key: v for key, v in new.items() if v[0] <= low + WINDOW}
    path = min(best.values(), key=lambda t: t[0])[1]  # an odd run may end the chain: it bends one bond
    out = set()
    for k in reversed(ks):
        if path[0]:
            out.add(k)
        path = path[1]
    return out


_RISE = {30: 1, 90: 2, 150: 1, 210: -1, 270: -2, 330: -1}   # half-bonds of height per bond direction


# --- geometry ---------------------------------------------------------------------------------------------
def _norm(v):
    d = hypot(*v) or 1.
    return v[0] / d, v[1] / d


def _rot(p, a):
    c, s = cos(a), sin(a)
    return p[0] * c - p[1] * s, p[0] * s + p[1] * c


def _seg_dist(p, a, b):
    dx, dy = b[0] - a[0], b[1] - a[1]
    t = max(0., min(1., ((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / (dx * dx + dy * dy or 1.)))
    return hypot(p[0] - a[0] - t * dx, p[1] - a[1] - t * dy)


def _orient(a, b, c):
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


def _cross(p1, p2, p3, p4):
    """A proper crossing.  Collinear lattice bonds give orientations of float noise, hence the tolerance;
    collinear overlap is `ON_BOND`'s to count."""
    return _orient(p1, p2, p3) * _orient(p1, p2, p4) < -1e-9 and _orient(p3, p4, p1) * _orient(p3, p4, p2) < -1e-9


def _reflect(pts, sub, p, q):
    """Reflect the atoms `sub` of `pts` about the line p-q, in place."""
    ex, ey = _norm((q[0] - p[0], q[1] - p[1]))
    for x in sub:
        dx, dy = pts[x][0] - p[0], pts[x][1] - p[1]
        k = dx * ex + dy * ey
        pts[x] = (p[0] + 2 * k * ex - dx, p[1] + 2 * k * ey - dy)


# --- placed atoms and the collision index -----------------------------------------------------------------
class _Scene:
    """Placed atoms, reserved slots and a uniform grid of cell `L` over atoms and bonds.

    Cross-links are never obstacles: they are routed at draw time.  `bad()` counts, for a candidate
    `{atom: point}`: atoms within `CLOSE` of a placed atom or a slot, new bonds crossing placed ones, and
    atoms within `ON_BOND` of a foreign bond, both ways.  `_check`, when given, receives the grid count and
    the brute-force count for every candidate.
    """
    __slots__ = ('mol', 'pos', 'virtual', 'slots', 'groups', 'links', '_atoms', '_bonds', '_check')

    def __init__(self, mol, links, check=None):
        self.mol = mol
        self.pos = {}
        self.virtual = {}
        self.slots = {}
        self.groups = []                             # placed `_Group`s, for relaxation
        self.links = links                           # frozensets of the cross-link bonds
        self._atoms = {}                             # cell -> atoms
        self._bonds = {}                             # cell -> bond frozensets
        self._check = check

    # grid
    @staticmethod
    def _cell(p):
        return int(p[0] // L), int(p[1] // L)

    def _bond_cells(self, p, q):
        (i0, j0), (i1, j1) = self._cell((min(p[0], q[0]), min(p[1], q[1]))), \
            self._cell((max(p[0], q[0]), max(p[1], q[1])))
        return [(i, j) for i in range(i0, i1 + 1) for j in range(j0, j1 + 1)]

    def _placed_bonds(self, x):
        return [frozenset((x, y)) for y in self.mol.neighbors_of(x)
                if y in self.pos and frozenset((x, y)) not in self.links]

    def commit(self, pts):
        for x, p in pts.items():
            if x in self.pos:
                self.drop((x,))
            self.pos[x] = p
            self._atoms.setdefault(self._cell(p), set()).add(x)
            for b in self._placed_bonds(x):
                u, v = b
                for c in self._bond_cells(self.pos[u], self.pos[v]):
                    self._bonds.setdefault(c, set()).add(b)

    def drop(self, atoms):
        for x in atoms:
            for b in self._placed_bonds(x):
                u, v = b
                for c in self._bond_cells(self.pos[u], self.pos[v]):
                    self._bonds[c].discard(b)
            self._atoms[self._cell(self.pos[x])].discard(x)
            del self.pos[x]

    # collisions
    def bad(self, pts, skip):
        ends = list(pts.values())                    # a new bond to a placed atom reaches past the points
        ends.extend(self.pos[y] for x in pts for y in self.mol.neighbors_of(x) if y in self.pos and y not in pts)
        xs = [p[0] for p in ends]
        ys = [p[1] for p in ends]
        (i0, j0), (i1, j1) = self._cell((min(xs) - L, min(ys) - L)), self._cell((max(xs) + L, max(ys) + L))
        atoms, bonds = set(), set()
        for i in range(i0, i1 + 1):
            for j in range(j0, j1 + 1):
                atoms.update(self._atoms.get((i, j), ()))
                bonds.update(self._bonds.get((i, j), ()))
        n = self._count(pts, skip, atoms, bonds)
        if self._check is not None:
            self._check(n, self._brute(pts, skip))
        return n

    def _brute(self, pts, skip):
        pos = self.pos
        bonds = {b for x in pos for b in self._placed_bonds(x)}
        return self._count(pts, skip, set(pos), bonds)

    def _count(self, pts, skip, atoms, bonds):
        mol, pos = self.mol, self.pos
        atoms = [(y, pos[y]) for y in atoms if y not in skip and y not in pts]
        old = [tuple(b) for b in bonds if not b & pts.keys()]
        new = [(x, y) for x in pts for y in mol.neighbors_of(x) if frozenset((x, y)) not in self.links and
               (y in pts and x < y or y in pos and y not in pts)]
        n = 0
        for p in pts.values():
            for q in self.slots.values():
                if hypot(p[0] - q[0], p[1] - q[1]) < CLOSE:
                    n += 1
            for _, q in atoms:
                if hypot(p[0] - q[0], p[1] - q[1]) < CLOSE:
                    n += 1
        for a, b in new:
            pa, pb = pts[a], pts[b] if b in pts else pos[b]
            for c, d in old:
                if len({a, b, c, d}) == 4 and _cross(pa, pb, pos[c], pos[d]):
                    n += 1
            for y, q in atoms:
                if y != a and y != b and _seg_dist(q, pa, pb) < ON_BOND:
                    n += 1
        for x, p in pts.items():
            for c, d in old:
                if x != c and x != d and _seg_dist(p, pos[c], pos[d]) < ON_BOND:
                    n += 1
        return n

    def nb(self, a):
        out = [self.pos[y] for y in self.mol.neighbors_of(a)
               if y in self.pos and frozenset((a, y)) not in self.links]
        if a in self.virtual:
            out.append(self.virtual[a])
        return out

    def outward(self, a):
        p = self.pos[a]
        s = [0., 0.]
        for q in self.nb(a):
            d = _norm((q[0] - p[0], q[1] - p[1]))
            s[0] -= d[0]
            s[1] -= d[1]
        return _norm(s)

    def free_directions(self, a, g):
        """`g` directions sharing the largest angular gap between `a`'s placed neighbours evenly."""
        p = self.pos[a]
        angs = sorted(atan2(q[1] - p[1], q[0] - p[0]) for q in self.nb(a))
        if not angs:
            return [(cos(2 * pi * q / g), sin(2 * pi * q / g)) for q in range(g)]
        if len(angs) == 1:
            if g == 1:                               # a chain end continues the zigzag at 120 degrees
                d = angs[0] + radians(120)
                return [(cos(d), sin(d))]
            lo, gap = angs[0], 2 * pi
        else:
            gap, lo = max(((angs[(q + 1) % len(angs)] - angs[q]) % (2 * pi), angs[q]) for q in range(len(angs)))
        return [(cos(lo + gap * (q + 1) / (g + 1)), sin(lo + gap * (q + 1) / (g + 1))) for q in range(g)]


# --- group tiles: native acyclic walk, engine ring systems ------------------------------------------------
def _ring_systems(mol, members):
    """`{atom: frozenset}`: rings inside one residue, those sharing an atom merged into one system.

    A ring closed through a cross-link or the chain is no tile: its link is routed at draw time.
    """
    systems = []
    for ring in mol.sssr:
        ring = set(ring)
        if not any(ring <= m for m in members):
            continue
        for s in [s for s in systems if s & ring]:
            systems.remove(s)
            ring |= s
        systems.append(ring)
    return {x: frozenset(s) for s in systems for x in s}


class _Tiles:
    """Ring-system tiles with a one-atom shell of ghost neighbours, cached by canonical form."""
    __slots__ = ('mol', 'tile', 'cache')

    def __init__(self, mol, tile):
        self.mol = mol
        self.tile = tile
        self.cache = {}

    def __call__(self, system):
        mol = self.mol
        shell = {y for x in system for y in mol.neighbors_of(x) if y not in system}
        sub = mol.substructure(sorted(system | shell))
        key = sub.canonical_bytes
        order = sub.canonical_order()
        if key not in self.cache:
            plane = self.tile(sub)
            self.cache[key] = {order[x]: plane[x] for x in sub}
        tile = self.cache[key]
        return {x: tile[order[x]] for x in sub}


def _subtree(mol, x, within, blocked):
    """Atoms reachable from `x` inside `within` without passing `blocked`."""
    seen, stack = {x}, [x]
    while stack:
        y = stack.pop()
        for z in mol.neighbors_of(y):
            if z in within and z not in seen and z not in blocked:
                seen.add(z)
                stack.append(z)
    return seen


def _kid_directions(k, d, s, linear):
    """`(direction, sign)` per child of an atom reached along `d`, whose next turn is `s` * 60 degrees."""
    if linear and k == 1:
        return [(d, s)]
    if k == 1:
        return [(d + 60 * s, -s)]
    if k == 2:
        return [(d + 60 * s, -s), (d - 60 * s, s)]
    if k == 3:
        return [(d, s), (d + 90, -1), (d - 90, 1)]
    return [(d + 180 + 360 * (q + 1) / (k + 1), 1) for q in range(k)]


def _group_tile(mol, group, a, rings, tiles, cis_trans):
    """`group` hung on `a` in a local frame: `a` at the origin, the slot axis along +x.

    Acyclic atoms walk the lattice at 120 degrees (straight through an sp atom); a ring system is its
    tile, mapped rigidly onto the bond that enters it.  Returns `(points, spiro)`, `spiro` when the
    anchor itself belongs to a ring system of the group.
    """
    pts, order, stack = {a: (0., 0.)}, {a: 0}, []

    def enter_ring(system, t, x0, frame):
        """Place `system` from tile `t` so that `frame(t point) -> local point`; queue its exits."""
        for y in system:
            pts[y] = frame(t[y])
            order[y] = len(order)
        for z in sorted(system):
            for w in mol.neighbors_of(z):
                if w in group and w not in system and w not in pts:
                    q = frame(t[w])
                    stack.append((z, w, degrees(atan2(q[1] - pts[z][1], q[0] - pts[z][0])), 1))

    ra = rings.get(a)
    spiro = ra is not None and bool(ra & group)
    if spiro:
        t = tiles(ra)
        cx = sum(t[y][0] for y in ra) / len(ra) - t[a][0]
        cy = sum(t[y][1] for y in ra) / len(ra) - t[a][1]
        ang = -atan2(cy, cx)
        ta = t[a]
        enter_ring(ra - {a}, t, a, lambda p: _rot((p[0] - ta[0], p[1] - ta[1]), ang))
    else:
        root = min(y for y in mol.neighbors_of(a) if y in group)
        stack.append((a, root, 0., 1))
    while stack:
        p, x, d, s = stack.pop()
        if x in pts:
            continue
        px = (pts[p][0] + L * cos(radians(d)), pts[p][1] + L * sin(radians(d)))
        system = rings.get(x)
        if system is not None:
            t = tiles(system)
            ang = radians(d) - atan2(t[x][1] - t[p][1], t[x][0] - t[p][0])
            tx = t[x]

            def frame(q, tx=tx, ang=ang, px=px):
                v = _rot((q[0] - tx[0], q[1] - tx[1]), ang)
                return px[0] + v[0], px[1] + v[1]
            enter_ring(system, t, x, frame)
            continue
        pts[x] = px
        order[x] = len(order)
        kids = [y for y in mol.neighbors_of(x) if y in group and y not in pts]
        if not kids:
            continue
        size = {y: len(_subtree(mol, y, group, pts.keys() | {x})) for y in kids}
        kids.sort(key=lambda y: (-size[y], y))
        for y, (dy, sy) in zip(kids, _kid_directions(len(kids), d, s, mol.hybridization_of(x) == 3)):
            stack.append((x, y, dy, sy))

    # a stated cis/trans bond reads back from the walk, or its far side is reflected about the bond
    for unit in cis_trans:
        frame = cis_trans_frame(mol, unit)
        if frame is None or not all(y in pts for y in frame):
            continue
        near, u, v, far = frame
        if cis_trans_parity(mol, unit, plane=pts) != unit['parity']:
            if order[u] > order[v]:
                u, v = v, u
            _reflect(pts, _subtree(mol, v, group, {u}) - {v}, pts[u], pts[v])
    return pts, spiro


def _torsions(mol, group, a, tile):
    """Single acyclic bonds `(u, v)` at depth 2 and 3 from `a`, with the subtree beyond `v`.

    Reflecting that subtree about u-v keeps every angle and every cis/trans bond.
    """
    root = [y for y in mol.neighbors_of(a) if y in group]
    depth, par, queue = {r: 1 for r in root}, {r: a for r in root}, deque(root)
    while queue:
        x = queue.popleft()
        for y in mol.neighbors_of(x):
            if y in group and y not in depth:
                depth[y], par[y] = depth[x] + 1, x
                queue.append(y)
    out = []
    for v, d in sorted(depth.items()):
        if d not in (2, 3):
            continue
        u = par[v]
        if mol.order_of(u, v) != 1 or mol.bond_in_ring(u, v):
            continue
        sub = _subtree(mol, v, group, {u})
        if len(sub) > 1:
            out.append((u, v, sub - {v}))
    return out


# --- placing one group ------------------------------------------------------------------------------------
class _Group:
    __slots__ = ('atoms', 'tile', 'spiro', 'torsions', 'anchor', 'direction', 'stretch')

    def __init__(self, atoms, tile, spiro, torsions, anchor, direction):
        self.atoms = atoms
        self.tile = tile
        self.spiro = spiro
        self.torsions = torsions
        self.anchor = anchor
        self.direction = direction
        self.stretch = 0.


def _conform(sc, g, stretches=STRETCHES):
    """`(bad, stretch, points)` for group `g` on its slot; commits nothing.

    Freedoms, cheapest first: reflection, one or two torsion flips, then a longer root bond.  The first
    collision-free candidate wins, otherwise the one with fewest collisions.
    """
    ax, ay = sc.pos[g.anchor]
    ang = atan2(g.direction[1], g.direction[0])
    combos = [()] + [(f,) for f in g.torsions] + \
        [(f, h) for i, f in enumerate(g.torsions) for h in g.torsions[i + 1:]]
    if g.spiro:
        stretches = (0.,)
    best = None
    for s in stretches:
        for refl in (1, -1):
            base = {}
            for x in g.atoms:
                tx, ty = g.tile[x]
                dx, dy = _rot((tx + s * L, ty * refl), ang)
                base[x] = (ax + dx, ay + dy)
            for combo in combos:
                pts = dict(base)
                for u, v, sub in combo:
                    _reflect(pts, sub, pts[u] if u in pts else sc.pos[u], pts[v])
                b = sc.bad(pts, {g.anchor})
                if best is None or b < best[0]:
                    best = (b, s, pts)
                if not b:
                    return best
    return best


def _place(sc, g):
    b, s, pts = _conform(sc, g)
    sc.commit(pts)
    g.stretch = s
    sc.groups.append(g)
    return b


def _relax(sc):
    """A stretched group re-placed at zero stretch together with one placed group anchored within 3 L."""
    for g in sc.groups:
        if not g.stretch:
            continue
        p = sc.pos[g.anchor]
        for h in sc.groups:
            q = sc.pos[h.anchor]
            if h is g or h.anchor in g.atoms or g.anchor in h.atoms or hypot(q[0] - p[0], q[1] - p[1]) > 3 * L:
                continue
            keep = {x: sc.pos[x] for x in g.atoms | h.atoms}
            sc.drop(keep)
            b1, _, p1 = _conform(sc, g, (0.,))
            if not b1:
                sc.commit(p1)
                b2, s2, p2 = _conform(sc, h, (0., h.stretch) if h.stretch else (0.,))
                if not b2:
                    sc.commit(p2)
                    g.stretch, h.stretch = 0., s2
                    break
                sc.drop(p1)
            sc.commit(keep)


def _groups(mol, atoms, placed):
    """Connected unplaced parts of `atoms`, each with its placed neighbours inside `atoms`."""
    seen, out = set(), []
    for s in sorted(atoms):
        if s in seen or s in placed:
            continue
        comp = _subtree(mol, s, atoms, placed)
        seen |= comp
        att = sorted({y for x in comp for y in mol.neighbors_of(x) if y in placed and y in atoms})
        out.append((comp, att))
    return out


def _ring_path(mol, a, b, comp):
    """Shortest a..b path through `comp`, never the direct a-b bond."""
    prev, queue = {a: None}, deque([a])
    while queue:
        x = queue.popleft()
        for y in mol.neighbors_of(x):
            if y == b and x != a:
                path = [b, x]
                while prev[path[-1]] is not None:
                    path.append(prev[path[-1]])
                return path[::-1]
            if y in comp and y not in prev:
                prev[y] = x
                queue.append(y)
    return None


def _polygon(sc, path, out):
    """A regular polygon on the placed edge path[0]-path[-1], on the side of unit vector `out`."""
    r = len(path)
    p0, p1 = sc.pos[path[0]], sc.pos[path[-1]]
    ex, ey = _norm((p1[0] - p0[0], p1[1] - p0[1]))
    nx, ny = -ey, ex
    if nx * out[0] + ny * out[1] < 0:
        nx, ny = -nx, -ny
    radius = L / (2 * sin(pi / r))
    h = radius * cos(pi / r)
    c = ((p0[0] + p1[0]) / 2 + nx * h, (p0[1] + p1[1]) / 2 + ny * h)
    t0 = atan2(p0[1] - c[1], p0[0] - c[0])
    t1 = atan2(p1[1] - c[1], p1[0] - c[0])
    d = (t0 - t1 + pi) % (2 * pi) - pi                # the step from p1 to p0, continued past p0
    sc.commit({x: (c[0] + radius * cos(t0 + i * d), c[1] + radius * sin(t0 + i * d))
               for i, x in enumerate(path[1:-1], 1)})


def _open_hexagon(sc, path, out):
    """A five-ring on the placed edge path[0]-path[-1] as the lattice hexagon on the side of `out`, less
    its top or bottom vertex: the chord is horizontal and every other ring angle is 120 degrees."""
    p0, p1 = sc.pos[path[0]], sc.pos[path[-1]]
    ex, ey = _norm((p1[0] - p0[0], p1[1] - p0[1]))
    nx, ny = -ey, ex
    if nx * out[0] + ny * out[1] < 0:
        nx, ny = -nx, -ny
    h = L * 3 ** .5 / 2
    c = ((p0[0] + p1[0]) / 2 + nx * h, (p0[1] + p1[1]) / 2 + ny * h)
    t0 = atan2(p0[1] - c[1], p0[0] - c[0])
    t1 = atan2(p1[1] - c[1], p1[0] - c[0])
    d = (t0 - t1 + pi) % (2 * pi) - pi                # the step from p1 to p0, continued past p0
    cell = [(c[0] + L * cos(t0 + i * d), c[1] + L * sin(t0 + i * d)) for i in range(1, 5)]
    gap = max(range(4), key=lambda i: abs(cell[i][1] - c[1]))          # the vertex straight above or below
    sc.commit(dict(zip(path[1:-1], (p for i, p in enumerate(cell) if i != gap))))


def _place_groups(sc, ctx, order, skip=frozenset()):
    """One pass: every unplaced part of each residue in `order` hung on its placed attachment.

    Returns `(collisions, progress)`.  A part closing a ring on two bonded placed atoms is a polygon on
    that edge; a part with any other two attachments waits for the engine.
    """
    mol, members = sc.mol, ctx.peptide.members
    bad, by_atom, closing = 0, {}, []
    for i in order:
        for comp, att in _groups(mol, members[i], sc.pos.keys()):
            if comp & skip or not att:
                continue
            if len(att) == 1:
                by_atom.setdefault(att[0], []).append(comp)
            elif len(att) == 2 and att[1] in mol.neighbors_of(att[0]):
                closing.append((comp, att))
    progress = False
    for comp, (a, b) in closing:
        path = _ring_path(mol, a, b, comp)
        if path:
            _polygon(sc, path, sc.outward(a))
            progress = True
    jobs = []
    for a, comps in by_atom.items():
        dirs = sc.free_directions(a, len(comps))
        comps.sort(key=len, reverse=True)
        slots = sorted(range(len(dirs)), key=lambda q: abs(q - (len(dirs) - 1) / 2))
        for comp, q in zip(comps, slots):
            jobs.append((len(comp), min(comp), comp, a, dirs[q]))
    jobs.sort(key=lambda t: t[:2])
    for q, (_, _, comp, a, d) in enumerate(jobs):
        sc.slots[q] = (sc.pos[a][0] + L * d[0], sc.pos[a][1] + L * d[1])
    for q, (_, _, comp, a, d) in enumerate(jobs):
        del sc.slots[q]
        tile, spiro = _group_tile(mol, comp, a, ctx.rings, ctx.tiles, ctx.cis_trans)
        bad += _place(sc, _Group(comp, tile, spiro, () if spiro else _torsions(mol, comp, a, tile), a, d))
    return bad, progress or bool(jobs)


def _snap(deg):
    return (round((deg - 30) / 60) * 60 + 30) % 360


def _branch_walk(sc, atoms, a, right, depth):
    """Drop from `a` along the lattice as an armchair until past `depth`, then run parallel as a zigzag."""
    o = sc.outward(a)
    down = o[1] < 0
    v = 270 if down else 90
    h = 1 if right else -1
    sgn = -1 if down else 1
    diag = (v - 60 * h * sgn) % 360                   # the armchair's drifting diagonal
    horiz = (diag - 60 * h * sgn) % 360               # its partner in the horizontal zigzag
    cur = _snap(degrees(atan2(o[1], o[0])))
    x, y = sc.pos[a]
    phase = 'drop'
    for q, at in enumerate(atoms):
        if q:
            if phase == 'drop':
                if (y < depth if down else y > depth) and cur == diag:
                    phase, cur = 'run', horiz
                else:
                    cur = diag if cur == v else v
            else:
                cur = diag if cur == horiz else horiz
        x, y = x + L * cos(radians(cur)), y + L * sin(radians(cur))
        sc.commit({at: (x, y)})


# --- the layout -------------------------------------------------------------------------------------------
class _Context:
    """What survives across re-layouts: the segmentation, ring systems, the tile cache, stated cis/trans, and
    each layout by its inverted vertices."""
    __slots__ = ('mol', 'peptide', 'rings', 'tiles', 'cis_trans', 'check', 'done')

    def __init__(self, mol, peptide, tile, check):
        self.mol = mol
        self.peptide = peptide
        self.rings = _ring_systems(mol, peptide.members)
        self.tiles = _Tiles(mol, tile)
        self.cis_trans = [u for u in mol.stereo_units() if u['kind'] == 1 and u['parity']]
        self.check = check
        self.done = {}


def peptide_layout(mol, peptide, tile, *, rounds=3, _check=None):
    """`{atom: (x, y)}` for every atom of `mol`, the backbone of `peptide` on the lattice, or None.

    None when a part bonded to the placed atoms fits no slot: a ring the backbone path runs through
    (4-aminobenzoic acid, isonipecotic acid) has its rest attached at two non-adjacent vertices.

    :param peptide: `monomers(mol)`, not None.
    :param tile: `tile(sub) -> {atom: (x, y)}` at mean bond `L`, for ring systems and leftover atoms.
    :param rounds: passes over the stretched anchors; each takes the first pin that improves on it.
    """
    ctx = _Context(mol, peptide, tile, _check)
    best = _layout(ctx, {})
    if best[0] is None:
        return None
    extra, hopeless, budget = {}, set(), BUDGET
    n = len(peptide.backbone)
    for _ in range(rounds):
        progress = False
        for v in sorted(best[2] - hopeless):
            if v not in best[2]:                      # a step near another anchor freed it
                continue
            for step in _steps(v, n):
                if step.keys() & extra.keys():
                    continue
                if not budget:
                    return _flatten(ctx, best, extra)
                budget -= 1
                res = _layout(ctx, {**extra, **step})
                if res[1] < best[1]:
                    best, extra, progress = res, {**extra, **step}, True
                    break
            else:
                hopeless.add(v)
        if not progress or not best[1][2]:
            break
    return _flatten(ctx, best, extra)


def _flatten(ctx, best, extra):
    """The plane of `best`, or of a lower backbone found by dropping `same` pairs one at a time, taken when it
    is no worse.  Only pairs the solve without any `same` breaks are tried."""
    bb = ctx.peptide.backbone
    ys = [best[0][x][1] for x in bb]
    low, same = max(ys) - min(ys), best[3]
    free = _layout(ctx, extra, same, low)[4]
    relax = set()
    for p in sorted(q for q in same if (q[0] in free) != (q[1] in free)):
        res = _layout(ctx, extra, relax | {p}, low)
        if res[0] is not None and res[1] <= best[1]:
            ys = [res[0][x][1] for x in bb]
            best, low = res, max(ys) - min(ys)
            relax.add(p)
    return best[0]


def _steps(v, n):
    """Pins to try near a stretched vertex `v`, nearest first: a step at j, or j kept uninverted."""
    out = []
    for j in sorted(range(max(1, v - 3), min(n - 2, v + 3) + 1), key=lambda j: (abs(j - v), j)):
        if j + 1 <= n - 2:
            out.append({j: True, j + 1: True})
        out.append({j: False})
        if j == n - 2:
            out.append({j: True})
    return out


def _host(mol, x, idx, links):
    """Backbone index of the first backbone atom BFS reaches from `x`, never across a cross-link."""
    seen, queue = {x}, deque([x])
    while queue:
        y = queue.popleft()
        if y in idx:
            return idx[y]
        for z in mol.neighbors_of(y):
            if z not in seen and frozenset((y, z)) not in links:
                seen.add(z)
                queue.append(z)
    return None


def _turn_leaves(mol, peptide, plane, idx, links):
    """At a three-bond link end, trade a terminal neighbour's slot for the link's when that points the link
    more toward its side of the backbone and the slot is clear: the leg then leaves straight for the bar.  A
    link end hanging as a leaf takes its chain's other slot when that leans it toward its partner on the same
    side, so the bracket's two legs mirror each other."""
    for q, (a, b) in enumerate(peptide.crosslinks):
        if peptide.head_to_tail and q == len(peptide.crosslinks) - 1:
            continue
        for x, o in ((a, b), (b, a)):
            _lean_in(mol, peptide, plane, idx, links, x, o)
        for x in (a, b):
            k = _host(mol, x, idx, links)
            if k is None or x in idx:
                continue
            side = plane[x][1] - plane[peptide.backbone[k]][1]
            ns = [y for y in mol.neighbors_of(x) if frozenset((x, y)) not in links]
            if not side or len(ns) != 2:
                continue
            px, py = plane[x]
            us = [_norm((plane[y][0] - px, plane[y][1] - py)) for y in ns]
            free = _norm((-us[0][0] - us[1][0], -us[0][1] - us[1][1]))
            for y, u in zip(ns, us):
                if mol.degree_of(y) != 1 or side * u[1] <= side * free[1] + .1:
                    continue
                r = hypot(plane[y][0] - px, plane[y][1] - py)
                t = (px + free[0] * r, py + free[1] * r)
                if all(hypot(t[0] - c[0], t[1] - c[1]) >= CLOSE for z, c in plane.items() if z not in (x, y)):
                    plane[y] = t
                    break


def _lean_in(mol, peptide, plane, idx, links, x, o):
    """Reflect the leaf link end `x` through its neighbour's bond axis when that brings it nearer `o` in x, keeps
    it on its side of the backbone and is clear."""
    ns = [y for y in mol.neighbors_of(x) if frozenset((x, y)) not in links]
    k, ko = _host(mol, x, idx, links), _host(mol, o, idx, links)
    if len(ns) != 1 or x in idx or k is None or ko is None:
        return
    y = ns[0]
    zs = [z for z in mol.neighbors_of(y) if z != x]
    if len(zs) != 1:
        return
    side = plane[x][1] - plane[peptide.backbone[k]][1]
    if side * (plane[o][1] - plane[peptide.backbone[ko]][1]) <= 0:
        return
    (px, py), (cx, cy), (zx, zy) = plane[x], plane[y], plane[zs[0]]
    ax, ay = _norm((cx - zx, cy - zy))
    vx, vy = px - cx, py - cy
    d = vx * ax + vy * ay
    t = (cx + 2 * d * ax - vx, cy + 2 * d * ay - vy)
    if abs(t[0] - plane[o][0]) >= abs(px - plane[o][0]) - .1 or side * (t[1] - cy) <= 0:
        return
    if all(hypot(t[0] - c[0], t[1] - c[1]) >= CLOSE for z, c in plane.items() if z != x):
        plane[x] = t


def _link_sides(mol, peptide, idx, links):
    """Pins `{vertex: must_invert}`: links shortest span first, each on the first side crossing no link."""
    n = len(peptide.backbone)
    spans = []
    for q, (a, b) in enumerate(peptide.crosslinks):
        if peptide.head_to_tail and q == len(peptide.crosslinks) - 1:
            continue
        ka, kb = _host(mol, a, idx, links), _host(mol, b, idx, links)
        if ka is not None and kb is not None:
            spans.append(tuple(sorted((ka, kb))))
    spans.sort(key=lambda t: t[1] - t[0])
    chosen, pins = [], {}
    for ka, kb in spans:
        for s in (_default_side(ka), -_default_side(ka)):
            if any(s == s2 and (k1 < ka < k2 < kb or ka < k1 < kb < k2) for k1, k2, s2 in chosen):
                continue
            if any(kx in pins and pins[kx] != (s != _default_side(kx)) for kx in (ka, kb)):
                continue
            chosen.append((ka, kb, s))
            for kx in (ka, kb):
                if 1 <= kx <= n - 2:
                    pins[kx] = s != _default_side(kx)
            break
    return pins


def _layout(ctx, extra, relax=(), below=None):
    """One layout with the extra pins `extra` and the `same` pairs in `relax` dropped:
    `(plane, (collisions, lost links, stretch), near, same, inverted)`, the plane None when atoms bonded to it
    are left over.  `(None, None, None, same, inverted)` when the backbone is not lower than `below`."""
    mol, peptide = ctx.mol, ctx.peptide
    bb, members = peptide.backbone, peptide.members
    n = len(bb)
    idx = {a: q for q, a in enumerate(bb)}
    links = frozenset(frozenset(x) for x in peptide.crosslinks)

    # a backbone ring keeps N and C-alpha on one side, and the vertices either side of it on the other
    differ, same, ring_jobs, pairs = set(), set(), [], []
    for i in peptide.main:
        for comp, att in _groups(mol, members[i], idx):
            ks = sorted(idx[a] for a in att)
            if len(ks) == 2 and ks[1] - ks[0] == 1:
                ring_jobs.append((comp, bb[ks[0]], bb[ks[1]]))
                pairs.append(ks)
                if ks[0] >= 1 and ks[1] <= n - 2:
                    differ.add((ks[0], ks[1]))
                same.update(p for p in ((ks[0] - 1, ks[0]), (ks[1], ks[1] + 1)) if p[0] >= 1 and p[1] <= n - 2)
    pins = _link_sides(mol, peptide, idx, links)
    for k, v in extra.items():
        pins.setdefault(k, v)
    starts = {k0 for k0, _ in pairs}
    bridges = {k1 + 1 for _, k1 in pairs if k1 + 2 in starts}
    same -= {p for k in bridges for p in ((k - 1, k), (k, k + 1))}     # a bridge clears its carbonyl
    inv = frozenset(_solve(n, pins, differ, same - set(relax), bridges=bridges))
    if below is not None:
        ys = [y for _, y in _walk((0., 0.), _backbone_dirs(n, inv), bridges)]
        if max(ys) - min(ys) > below - 1e-6:
            return None, None, None, same, inv
    if inv not in ctx.done:
        ctx.done[inv] = _draw(ctx, inv, bridges, ring_jobs, idx, links) + (same, inv)
    return ctx.done[inv]


def _draw(ctx, inv, bridges, ring_jobs, idx, links):
    """The layout of `_layout` for the inverted vertices `inv`, less its `same`."""
    mol, peptide = ctx.mol, ctx.peptide
    bb, members, n = peptide.backbone, peptide.members, len(peptide.backbone)
    dirs = _backbone_dirs(n, inv)

    sc = _Scene(mol, links, ctx.check)
    sc.commit(dict(zip(bb, _walk((0., 0.), dirs, bridges))))
    if peptide.head_to_tail and dirs:
        last = dirs[-1]
        for e, d in ((bb[0], 270), (bb[-1], min((last + 60, last - 60), key=lambda t: sin(radians(t))))):
            sc.virtual[e] = (sc.pos[e][0] + L * cos(radians(d)), sc.pos[e][1] + L * sin(radians(d)))
    for comp, a, b in ring_jobs:
        path = _ring_path(mol, a, b, comp)
        if path and len(path) == 5:
            _open_hexagon(sc, path, sc.outward(b if idx[a] == 0 else a))
        elif path:
            _polygon(sc, path, sc.outward(b if idx[a] == 0 else a))

    # branch spurs are reserved while the main chain's side chains go down
    walks, reserved = [], set()
    for br in peptide.branches:
        free = {x for x in members[br.residue] if x not in sc.pos}
        spur = None
        for a in sorted(x for x in members[br.residue] if x in sc.pos):
            p = _bfs_path(mol, a, br.host, free)
            if p and (spur is None or len(p) < len(spur)):
                spur = p
        if spur is None or spur[0] not in idx:
            continue
        atoms = spur[1:] + list(br.path)
        walks.append((spur[0], atoms))
        reserved.update(atoms)
    bad, _ = _place_groups(sc, ctx, peptide.main, reserved)
    for a, atoms in walks:
        ys = [p[1] for p in sc.pos.values()]
        down = sc.outward(a)[1] < 0
        k = idx[a]
        _branch_walk(sc, atoms, a, n - 1 - k > k, min(ys) - 1.2 * L if down else max(ys) + 1.2 * L)

    order = list(peptide.main) + [j for br in peptide.branches for j in br.chain]
    order += [j for j in range(len(members)) if j not in set(order)]
    for _ in range(6):
        b, progress = _place_groups(sc, ctx, order)
        bad += b
        if not progress:
            break
    _relax(sc)
    plane = dict(sc.pos)
    _turn_leaves(mol, peptide, plane, idx, links)
    rest = [x for x in mol if x not in plane]
    if any(y in plane for x in rest for y in mol.neighbors_of(x)):
        return None, (bad, 0, 0.), set()
    if rest:
        extra_plane = ctx.tiles.tile(mol.substructure(rest))
        mx = max(x for x, _ in plane.values()) + 2
        plane.update({x: (px + mx, py) for x, (px, py) in extra_plane.items()})
    stretch = sum(g.stretch for g in sc.groups)
    near = {idx[g.anchor] for g in sc.groups if g.stretch and g.anchor in idx}
    return plane, (bad, _lost(mol, peptide, plane, idx, links), stretch), near


def _lost(mol, peptide, plane, idx, links):
    """Links whose two ends leave the backbone on different sides, or level with it: no bracket takes them."""
    lost = 0
    for q, (a, b) in enumerate(peptide.crosslinks):
        if peptide.head_to_tail and q == len(peptide.crosslinks) - 1:
            continue
        ka, kb = _host(mol, a, idx, links), _host(mol, b, idx, links)
        if ka is None or kb is None:
            continue
        sa, sb = (plane[x][1] - plane[peptide.backbone[k]][1] for x, k in ((a, ka), (b, kb)))
        lost += not (sa * sb > 0 and min(abs(sa), abs(sb)) > 1e-3)
    return lost
