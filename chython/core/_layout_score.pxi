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
# LAYOUT SCORE: `LayoutScore`, the seven terms of a 2D layout's score, and their change under a move.
#
# A geometry kernel for `chython.depict.layout.score`: it reads flat arrays and knows nothing of a molecule. Atoms
# are positions 0..n-1 in the caller's canonical order, coordinates `double[2n]`, and every threshold arrives in
# `thresholds`, so the layout's constants stay in `depict/layout/`. The arithmetic follows the Python oracle line for
# line: `count_defects` (pairs once, each crossing counted from both of its bonds), `_crowded`, `_bent`,
# `_straight` and `_polygon_rmsd`, over the region `relax._near` keeps.

from libc.math cimport atan2, cos, fabs, hypot, sin, fmod


# --- thresholds, in the order `LayoutScore` takes them ------------------------------------------------
DEF T_L = 0
DEF T_CLOSE = 1
DEF T_SHORT = 2
DEF T_ON_BOND = 3
DEF T_CROWD = 4
DEF T_STRETCHED = 5
DEF T_BENT = 6          # degrees
DEF T_STRAIGHT = 7      # cosine
DEF T_MARGIN = 8
DEF T_COUNT = 9
DEF TERMS = 7
DEF PI = 3.141592653589793


cdef inline double _len2(const double *xy, int32_t n, int32_t m) noexcept nogil:
    return hypot(xy[2 * n] - xy[2 * m], xy[2 * n + 1] - xy[2 * m + 1])


cdef inline double _side(double ox, double oy, double ux, double uy, double vx, double vy) noexcept nogil:
    return (ux - ox) * (vy - oy) - (uy - oy) * (vx - ox)


cdef inline bint _crossing(const double *xy, int32_t a, int32_t b, int32_t c, int32_t e) noexcept nogil:
    cdef double px = xy[2 * a], py = xy[2 * a + 1], qx = xy[2 * b], qy = xy[2 * b + 1]
    cdef double ax = xy[2 * c], ay = xy[2 * c + 1], bx = xy[2 * e], by = xy[2 * e + 1]
    return _side(px, py, qx, qy, ax, ay) * _side(px, py, qx, qy, bx, by) < -1e-12 and \
        _side(ax, ay, bx, by, px, py) * _side(ax, ay, bx, by, qx, qy) < -1e-12


cdef inline double _seg(double x, double y, double px, double py, double qx, double qy) noexcept nogil:
    cdef double dx = qx - px, dy = qy - py, d = dx * dx + dy * dy, t = 0.
    if d:
        t = ((x - px) * dx + (y - py) * dy) / d
        t = 0. if t < 0. else (1. if t > 1. else t)
    return hypot(px + t * dx - x, py + t * dy - y)


cdef class LayoutScore:
    """The score terms of one component's planes: `defects`, `closure`, `crowd`, `stretched`, `bent`, `ring` and
    `straight`, as `chython.depict.layout.score` defines them.

    `bonds`, `closures` and `drawn`: flat endpoint pairs. `adj` and `near2`: CSR offsets and members, `near2` every
    atom within two bonds, itself included. `bent` and `straight`: per atom, 1 where the term may hold. `rings`: CSR
    of the rings of at most 8 atoms, each in ring order. `linear` is not needed: the flags carry it."""
    cdef int32_t n, nb
    cdef int32_t[::1] bonds, closures, drawn, adj_off, adj, near2_off, near2, ring_off, ring
    cdef uint8_t[::1] bent, straight
    cdef double t[T_COUNT]

    def __init__(self, int32_t n, bonds, closures, drawn, adj_off, adj, near2_off, near2, bent, straight, ring_off,
                 ring, thresholds):
        cdef int k
        self.n = n
        self.bonds, self.closures, self.drawn = bonds, closures, drawn
        self.nb = len(self.bonds) // 2
        self.adj_off, self.adj, self.near2_off, self.near2 = adj_off, adj, near2_off, near2
        self.bent, self.straight = bent, straight
        self.ring_off, self.ring = ring_off, ring
        if len(thresholds) != T_COUNT:
            raise ValueError(f'{T_COUNT} thresholds, not {len(thresholds)}')
        for k in range(T_COUNT):
            self.t[k] = thresholds[k]

    cdef bint _bonded(self, int32_t a, int32_t b) noexcept:
        cdef int32_t k
        for k in range(self.adj_off[a], self.adj_off[a + 1]):
            if self.adj[k] == b:
                return True
        return False

    cdef bint _within2(self, int32_t a, int32_t b) noexcept:
        cdef int32_t k
        for k in range(self.near2_off[a], self.near2_off[a + 1]):
            if self.near2[k] == b:
                return True
        return False

    cdef double _bent_at(self, const double *xy, int32_t x) noexcept:
        cdef int32_t lo = self.adj_off[x], hi = self.adj_off[x + 1], d = hi - lo, i, j
        cdef double angs[3]
        cdef double g, dev, a, w
        if not self.bent[x]:
            return 0.
        for i in range(d):
            a = atan2(xy[2 * self.adj[lo + i] + 1] - xy[2 * x + 1], xy[2 * self.adj[lo + i]] - xy[2 * x]) * 180. / PI
            a = fmod(a, 360.)
            if a < 0.:
                a += 360.
            angs[i] = a
        for i in range(1, d):                               # sort the two or three angles
            w = angs[i]
            j = i - 1
            while j >= 0 and angs[j] > w:
                angs[j + 1] = angs[j]
                j -= 1
            angs[j + 1] = w
        dev = 0.
        if d == 2:
            g = fmod(angs[1] - angs[0], 360.)
            if g < 0.:
                g += 360.
            w = fmod(angs[0] - angs[1], 360.)
            if w < 0.:
                w += 360.
            dev = fabs((g if g < w else w) - 120.)
        else:
            for i in range(3):
                g = fmod(angs[(i + 1) % 3] - angs[i], 360.)
                if g < 0.:
                    g += 360.
                if fabs(g - 120.) > dev:
                    dev = fabs(g - 120.)
        return 1. if dev > self.t[T_BENT] else 0.

    cdef double _straight_at(self, const double *xy, int32_t x) noexcept:
        cdef int32_t a, b
        cdef double px, py, u
        if not self.straight[x]:
            return 0.
        a, b = self.adj[self.adj_off[x]], self.adj[self.adj_off[x] + 1]
        px, py = xy[2 * x], xy[2 * x + 1]
        u = hypot(xy[2 * a] - px, xy[2 * a + 1] - py) * hypot(xy[2 * b] - px, xy[2 * b + 1] - py)
        if not u:
            return 0.
        return 1. if ((xy[2 * a] - px) * (xy[2 * b] - px) + (xy[2 * a + 1] - py) * (xy[2 * b + 1] - py)) / u < \
            self.t[T_STRAIGHT] else 0.

    cdef double _polygon(self, const double *xy, int32_t lo, int32_t hi) noexcept:
        cdef int32_t m = hi - lo, k, w
        cdef double cx = 0., cy = 0., best = -1., a, b, s, h, c, tt, r, qx, qy, px, py
        for k in range(lo, hi):
            cx += xy[2 * self.ring[k]]
            cy += xy[2 * self.ring[k] + 1]
        cx /= m
        cy /= m
        s = self.t[T_L] / (2 * sin(PI / m))
        for w in (1, -1):
            a = b = 0.
            for k in range(m):
                px, py = xy[2 * self.ring[lo + k]] - cx, xy[2 * self.ring[lo + k] + 1] - cy
                qx, qy = cos(2 * PI * w * k / m), sin(2 * PI * w * k / m)
                a += px * qx + py * qy
                b += py * qx - px * qy
            h = hypot(a, b)
            if h:
                c, tt = a / h, b / h
            else:
                c, tt = 1., 0.
            r = 0.
            for k in range(m):
                px, py = xy[2 * self.ring[lo + k]] - cx, xy[2 * self.ring[lo + k] + 1] - cy
                qx, qy = cos(2 * PI * w * k / m), sin(2 * PI * w * k / m)
                r += (px - s * (c * qx - tt * qy)) ** 2 + (py - s * (tt * qx + c * qy)) ** 2
            if best < 0. or r < best:
                best = r
        return (best / m) ** .5 / self.t[T_L]

    cdef void _parts(self, const double *xy, const uint8_t *part, const int32_t *atoms, int32_t na,
                     const int32_t *bonds, int32_t nbk, double *out) noexcept:
        """The seven terms of `xy` over the atoms and bonds given, those of `part` where the term is per move."""
        cdef double close = self.t[T_CLOSE], short = self.t[T_SHORT], onb = self.t[T_ON_BOND], x0, x1, y0, y1
        cdef double crowd_d = self.t[T_CROWD], u0, u1, v0, v1, px, py
        cdef int32_t i, j, a, b, c, e, k, x, nring
        cdef double defects = 0., closure = 0., crowd = 0., stretched = 0., bent = 0., ring = 0., straight = 0.
        for i in range(na):                                 # close pairs, each once
            a = atoms[i]
            for j in range(i + 1, na):
                b = atoms[j]
                if fabs(xy[2 * a] - xy[2 * b]) < close and fabs(xy[2 * a + 1] - xy[2 * b + 1]) < close:
                    if _len2(xy, a, b) < (short if self._bonded(a, b) else close):
                        defects += 1
        for i in range(nbk):                                # crossings, from both bonds; atoms on bonds
            a, b = bonds[2 * i], bonds[2 * i + 1]
            x0 = (xy[2 * a] if xy[2 * a] < xy[2 * b] else xy[2 * b]) - onb
            x1 = (xy[2 * a] if xy[2 * a] > xy[2 * b] else xy[2 * b]) + onb
            y0 = (xy[2 * a + 1] if xy[2 * a + 1] < xy[2 * b + 1] else xy[2 * b + 1]) - onb
            y1 = (xy[2 * a + 1] if xy[2 * a + 1] > xy[2 * b + 1] else xy[2 * b + 1]) + onb
            for j in range(nbk):
                if j == i:
                    continue
                c, e = bonds[2 * j], bonds[2 * j + 1]
                if a == c or a == e or b == c or b == e:
                    continue
                u0 = (xy[2 * c] if xy[2 * c] < xy[2 * e] else xy[2 * e]) - onb
                u1 = (xy[2 * c] if xy[2 * c] > xy[2 * e] else xy[2 * e]) + onb
                v0 = (xy[2 * c + 1] if xy[2 * c + 1] < xy[2 * e + 1] else xy[2 * e + 1]) - onb
                v1 = (xy[2 * c + 1] if xy[2 * c + 1] > xy[2 * e + 1] else xy[2 * e + 1]) + onb
                if u0 <= x1 and x0 <= u1 and v0 <= y1 and y0 <= v1 and _crossing(xy, a, b, c, e):
                    defects += 1
            for k in range(na):
                x = atoms[k]
                if x == a or x == b:
                    continue
                px, py = xy[2 * x], xy[2 * x + 1]
                if x0 <= px <= x1 and y0 <= py <= y1 and \
                        _seg(px, py, xy[2 * a], xy[2 * a + 1], xy[2 * b], xy[2 * b + 1]) < onb:
                    defects += 1
        for i in range(len(self.closures) // 2):
            a, b = self.closures[2 * i], self.closures[2 * i + 1]
            if part[a] or part[b]:
                closure += fabs(_len2(xy, a, b) - self.t[T_L])
        for i in range(self.n):                             # crowded pairs with an atom of `part`, each once
            if not part[i]:
                continue
            for k in range(na):
                b = atoms[k]
                if self._within2(i, b) or (part[b] and b <= i):
                    continue
                if fabs(xy[2 * i] - xy[2 * b]) < crowd_d and fabs(xy[2 * i + 1] - xy[2 * b + 1]) < crowd_d and \
                        _len2(xy, i, b) < crowd_d:
                    crowd += 1
        for i in range(len(self.drawn) // 2):
            a, b = self.drawn[2 * i], self.drawn[2 * i + 1]
            if (part[a] or part[b]) and _len2(xy, a, b) > self.t[T_STRETCHED]:
                stretched += 1
        for x in range(self.n):                             # part and its neighbours: bent and straight atoms
            if not part[x]:
                for k in range(self.adj_off[x], self.adj_off[x + 1]):
                    if part[self.adj[k]]:
                        break
                else:
                    continue
            bent += self._bent_at(xy, x)
            straight += self._straight_at(xy, x)
        nring = len(self.ring_off) - 1
        for i in range(nring):
            for k in range(self.ring_off[i], self.ring_off[i + 1]):
                if part[self.ring[k]]:
                    ring += self._polygon(xy, self.ring_off[i], self.ring_off[i + 1])
                    break
        out[0], out[1], out[2], out[3], out[4], out[5], out[6] = \
            defects, closure, crowd, stretched, bent, ring, straight

    def terms(self, double[::1] xy):
        """The seven terms of the whole component drawn as `xy`."""
        cdef int32_t i
        cdef double out[TERMS]
        cdef uint8_t *part = <uint8_t *> malloc(self.n + 4 * self.n + 1)
        if part == NULL:
            raise MemoryError()
        cdef int32_t *atoms = <int32_t *> (part + self.n + (4 - self.n % 4) % 4)
        try:
            for i in range(self.n):
                part[i] = 1
                atoms[i] = i
            self._parts(&xy[0], part, atoms, self.n, &self.bonds[0] if self.nb else NULL, self.nb, out)
        finally:
            free(part)
        return int(out[0]), out[1], int(out[2]), int(out[3]), int(out[4]), out[5], int(out[6])

    def delta(self, double[::1] xy, int32_t[::1] moved, double[::1] to):
        """The change of `terms` when the atoms `moved` go to `to` (`double[2 len(moved)]`), over the region
        `relax._near` keeps: atoms and bonds within the margin of a moved atom or of a bond holding one, in either
        plane."""
        cdef int32_t n = self.n, i, k, a, b, na = 0, nbk = 0
        cdef double x0, x1, y0, y1, m = self.t[T_MARGIN], ax, ay, lo, hi
        cdef double before[TERMS]
        cdef double after[TERMS]
        cdef size_t scratch = 2 * n * sizeof(double) + n * sizeof(int32_t) + 2 * self.nb * sizeof(int32_t) + 2 * n
        cdef char *block = <char *> malloc(scratch + 16)
        if block == NULL:
            raise MemoryError()
        cdef double *new = <double *> block
        cdef int32_t *atoms = <int32_t *> (block + 2 * n * sizeof(double))
        cdef int32_t *bonds = atoms + n
        cdef uint8_t *part = <uint8_t *> (bonds + 2 * self.nb)
        cdef uint8_t *ends = part + n
        try:
            memcpy(new, &xy[0], 2 * n * sizeof(double))
            memset(part, 0, 2 * n)
            for i in range(len(moved)):
                part[moved[i]] = 1
                new[2 * moved[i]], new[2 * moved[i] + 1] = to[2 * i], to[2 * i + 1]
            for i in range(n):
                ends[i] = part[i]
            for i in range(self.nb):
                a, b = self.bonds[2 * i], self.bonds[2 * i + 1]
                if part[a] or part[b]:
                    ends[a] = ends[b] = 1
            x0 = y0 = 1e300
            x1 = y1 = -1e300
            for i in range(n):
                if ends[i]:
                    for k in range(2):
                        ax = xy[2 * i] if k == 0 else new[2 * i]
                        ay = xy[2 * i + 1] if k == 0 else new[2 * i + 1]
                        x0, x1 = (ax if ax < x0 else x0), (ax if ax > x1 else x1)
                        y0, y1 = (ay if ay < y0 else y0), (ay if ay > y1 else y1)
            x0, x1, y0, y1 = x0 - m, x1 + m, y0 - m, y1 + m
            for i in range(n):
                if (x0 <= xy[2 * i] <= x1 and y0 <= xy[2 * i + 1] <= y1) or \
                        (x0 <= new[2 * i] <= x1 and y0 <= new[2 * i + 1] <= y1):
                    atoms[na] = i
                    na += 1
            for i in range(self.nb):
                a, b = self.bonds[2 * i], self.bonds[2 * i + 1]
                lo = min(min(xy[2 * a], xy[2 * b]), min(new[2 * a], new[2 * b]))
                hi = max(max(xy[2 * a], xy[2 * b]), max(new[2 * a], new[2 * b]))
                ax = min(min(xy[2 * a + 1], xy[2 * b + 1]), min(new[2 * a + 1], new[2 * b + 1]))
                ay = max(max(xy[2 * a + 1], xy[2 * b + 1]), max(new[2 * a + 1], new[2 * b + 1]))
                if lo <= x1 and hi >= x0 and ax <= y1 and ay >= y0:
                    bonds[2 * nbk], bonds[2 * nbk + 1] = a, b
                    nbk += 1
            self._parts(&xy[0], part, atoms, na, bonds, nbk, before)
            self._parts(new, part, atoms, na, bonds, nbk, after)
        finally:
            free(block)
        return (int(after[0] - before[0]), after[1] - before[1], int(after[2] - before[2]),
                int(after[3] - before[3]), int(after[4] - before[4]), after[5] - before[5],
                int(after[6] - before[6]))
