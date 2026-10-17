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
"""The triangular lattice T of spacing `L`: honeycomb vertices plus hexagon centres.

| term | meaning |
| --- | --- |
| point | axial integer `(q, r)`; plane `x = L q cos 30`, `y = L (q / 2 + r)` |
| `DIRS[k]` | unit step at `30 + 60 k` degrees |
| sublattice | `(q - r) % 3`: 0 and 1 are the honeycomb, 2 the hexagon centres |
| honeycomb step | from sublattice 0 along `DIRS[0, 2, 4]`, from 1 along `DIRS[1, 3, 5]` |
| `Board` | atoms at points and bonds as segments; a placement that would add a defect is refused |

A hexagon centred at `c` has vertices `c + DIRS[k]` in order: vertical sides, a vertex at top.
"""
from math import hypot, sqrt


__all__ = ['L', 'DIRS', 'Board', 'add', 'sub', 'scale', 'project', 'sublattice', 'honeycomb_dirs', 'direction',
           'rotate', 'reflect', 'apply', 'SYMMETRIES']

L = .825
C30 = sqrt(3) / 2
DIRS = ((1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1))
CLOSE = .8 * L          # atom to atom
ON_BOND = .55 * L       # atom to a bond it is not in: refuses the flanks of a chord, at .5 L
CELL = 2 * L            # bucket side of the spatial hash


def add(p, q):
    return p[0] + q[0], p[1] + q[1]


def sub(p, q):
    return p[0] - q[0], p[1] - q[1]


def scale(p, k):
    return p[0] * k, p[1] * k


def project(p):
    return L * C30 * p[0], L * (p[0] / 2 + p[1])


def sublattice(p):
    return (p[0] - p[1]) % 3


def honeycomb_dirs(p):
    """The three honeycomb step indices of `p`, empty for a hexagon centre."""
    s = sublattice(p)
    return (0, 2, 4) if s == 0 else (1, 3, 5) if s == 1 else ()


def direction(v):
    """`(k, n)` when `v` is `n` unit steps along `DIRS[k]`, else None."""
    q, r = v
    if q and not r:
        return (0, q) if q > 0 else (3, -q)
    if r and not q:
        return (1, r) if r > 0 else (4, -r)
    if q and q == -r:
        return (5, q) if q > 0 else (2, -q)
    return None


def rotate(p, k):
    """`p` turned by `60 k` degrees about the origin."""
    q, r = p
    for _ in range(k % 6):
        q, r = -r, q + r
    return q, r


def reflect(p):
    """`p` mirrored in the horizontal axis: `DIRS[k]` goes to `DIRS[5 - k]`."""
    q, r = p
    return q, -q - r


# the 12 lattice symmetries fixing the origin, as (mirror, turn)
SYMMETRIES = tuple((m, k) for m in (False, True) for k in range(6))


def apply(p, g):
    m, k = g
    return rotate(reflect(p) if m else p, k)


def _orient(a, b, c):
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


def _seg_dist(x, p, q):
    dx, dy = q[0] - p[0], q[1] - p[1]
    d = dx * dx + dy * dy
    t = 0. if not d else max(0., min(1., ((x[0] - p[0]) * dx + (x[1] - p[1]) * dy) / d))
    return hypot(p[0] + t * dx - x[0], p[1] + t * dy - x[1])


def _crosses(p, q, a, b):
    """Proper crossing of segments `pq` and `ab`, exact on integer axial points."""
    return (_orient(p, q, a) * _orient(p, q, b) < 0) and (_orient(a, b, p) * _orient(a, b, q) < 0)


class Board:
    """Atoms at lattice points and bonds between them; `bad_*` count the defects a placement would add."""
    __slots__ = ('at', 'pos', 'bonds', '_xy', '_bxy', '_abuck', '_bbuck')

    def __init__(self):
        self.at = {}            # point -> atoms there; more than one only on a forced placement
        self.pos = {}           # atom -> point
        self.bonds = {}         # (n, m) sorted -> (p, q)
        self._xy = {}           # atom -> projected point
        self._bxy = {}          # bond -> projected ends
        self._abuck = {}        # cell -> atoms; each atom in exactly one
        self._bbuck = {}        # cell -> bonds; a bond in each cell its widened box meets

    @staticmethod
    def _key(xy):
        return int(xy[0] // CELL), int(xy[1] // CELL)

    def _bond_keys(self, p, q):
        (x1, y1), (x2, y2) = project(p), project(q)
        m = ON_BOND
        i1, j1 = self._key((min(x1, x2) - m, min(y1, y2) - m))
        i2, j2 = self._key((max(x1, x2) + m, max(y1, y2) + m))
        return [(i, j) for i in range(i1, i2 + 1) for j in range(j1, j2 + 1)]

    def near_atoms(self, xy, pad=1):
        i, j = self._key(xy)
        for di in range(-pad, pad + 1):
            for dj in range(-pad, pad + 1):
                yield from self._abuck.get((i + di, j + dj), ())

    def near_bonds(self, keys):
        seen = set()
        for k in keys:
            for b in self._bbuck.get(k, ()):
                if b not in seen:
                    seen.add(b)
                    yield b

    # --- tests --------------------------------------------------------------------------------------------

    def bad_atom(self, p):
        """Defects an atom at `p` would add against the board: a taken point, a bond through it."""
        if p in self.at:
            return 1
        xy = project(p)
        bad = 0
        for b in self._bbuck.get(self._key(xy), ()):
            u, v = self._bxy[b]
            if _seg_dist(xy, u, v) < ON_BOND:
                bad += 1
        return bad

    def bad_bond(self, p, q):
        """Defects a bond `pq` would add: an atom on it, a proper crossing; its own end points excepted."""
        bad = 0
        xp, xq = project(p), project(q)
        x0, x1 = min(xp[0], xq[0]) - ON_BOND, max(xp[0], xq[0]) + ON_BOND
        y0, y1 = min(xp[1], xq[1]) - ON_BOND, max(xp[1], xq[1]) + ON_BOND
        keys = self._bond_keys(p, q)
        pos, xys = self.pos, self._xy
        for k in keys:
            for x in self._abuck.get(k, ()):
                xy = xys[x]
                if x0 <= xy[0] <= x1 and y0 <= xy[1] <= y1 and pos[x] != p and pos[x] != q and \
                        _seg_dist(xy, xp, xq) < ON_BOND:
                    bad += 1
        q0, q1 = (p[0], q[0]) if p[0] < q[0] else (q[0], p[0])     # axial box: a proper crossing overlaps it
        r0, r1 = (p[1], q[1]) if p[1] < q[1] else (q[1], p[1])
        for b in self.near_bonds(keys):
            u, v = self.bonds[b]
            if (u[0] > q0 or v[0] > q0) and (u[0] < q1 or v[0] < q1) and (u[1] > r0 or v[1] > r0) and \
                    (u[1] < r1 or v[1] < r1) and u != p and u != q and v != p and v != q and _crosses(p, q, u, v):
                bad += 1
        return bad

    # --- edits --------------------------------------------------------------------------------------------

    def put(self, x, p):
        self.at.setdefault(p, []).append(x)
        self.pos[x] = p
        self._xy[x] = xy = project(p)
        self._abuck.setdefault(self._key(xy), []).append(x)

    def take(self, x):
        p = self.pos.pop(x)
        here = self.at[p]
        here.remove(x)
        if not here:
            del self.at[p]
        self._abuck[self._key(self._xy.pop(x))].remove(x)

    def link(self, n, m):
        b = (n, m) if n < m else (m, n)
        if b in self.bonds:
            return
        p, q = self.pos[n], self.pos[m]
        self.bonds[b] = (p, q)
        self._bxy[b] = project(p), project(q)
        for k in self._bond_keys(p, q):
            self._bbuck.setdefault(k, []).append(b)

    def unlink(self, n, m):
        b = (n, m) if n < m else (m, n)
        p, q = self.bonds.pop(b)
        del self._bxy[b]
        for k in self._bond_keys(p, q):
            self._bbuck[k].remove(b)
