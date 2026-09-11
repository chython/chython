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
"""Scalar field interpolation and smooth isoline extraction.

Two rules the rest follows from: an atom with NO value contributes NOTHING and is not zero, so an absent
atom cannot pull the field toward a neutral midpoint; and a traced vertex is Newton-refined onto the true
isoline before it is fitted, since marching squares places it by linear interpolation with O(h²) error.
"""
from dataclasses import dataclass, field
from math import exp, hypot, sqrt

from .scene import Box, move, curve, close


__all__ = ['ScalarField', 'Grid', 'contour_levels', 'convex_hull', 'in_polygon',
           'isolines', 'refine', 'sample', 'to_cubics', 'trim_asymptote']


class Grid:
    """Sampled scalar field on a regular grid.

    `min_x`/`min_y` is the world-coordinate origin of cell (0, 0), `spacing` the cell size in x and y,
    `nx`/`ny` the column and row counts.  `z` is row-major -- `z[j * nx + i]` is cell (i, j) -- and None
    there means "field undefined here".
    """
    __slots__ = ('min_x', 'min_y', 'spacing', 'nx', 'ny', 'z')

    def __init__(self, min_x: float, min_y: float, spacing: float, nx: int, ny: int,
                 z: list[float | None]):
        self.min_x = min_x
        self.min_y = min_y
        self.spacing = spacing
        self.nx = nx
        self.ny = ny
        self.z = z

    def __iter__(self):
        yield self.min_x
        yield self.min_y
        yield self.spacing
        yield self.nx
        yield self.ny
        yield self.z


@dataclass(frozen=True, slots=True)
class ScalarField:
    """Shepard / Gaussian-kernel scalar field over named atom positions.

    :param values: `{n: value}`.  Only atoms listed here contribute to the field.
    :param plane: `{n: (x, y)}` for ALL atoms in the plane -- a superset of `values`.
    :param sigma: Gaussian width, in molecule units.
    :param cutoff: beyond this distance from every atom in `values`, `at()` returns None.
    :param hull: optional counter-clockwise convex hull polygon; points outside it return None.
    """
    values: dict[int, float]
    plane: dict[int, tuple[float, float]]
    sigma: float = field(kw_only=True)
    cutoff: float = field(kw_only=True)
    hull: tuple[tuple[float, float], ...] | None = None

    def __post_init__(self):
        missing = set(self.values) - set(self.plane)
        if missing:
            first = next(iter(sorted(missing)))
            raise ValueError(
                f'value given for atom id {first} which is not in the plane')

    def at(self, x: float, y: float) -> float | None:
        """Gaussian-weighted field value at (x, y), or None if outside all cutoffs / hull.

            f(p) = Σᵢ vᵢ · exp( −‖p − pᵢ‖² / 2σ² )

        Deliberately UN-normalized (no ÷ Σ wᵢ): normalising would make a single atom's field constant
        everywhere (f = wᵢ·vᵢ / wᵢ = vᵢ), producing no isoline at any level below vᵢ and rendering a
        decaying property as a flat disk.  The un-normalized sum decays to zero far from all atoms.
        """
        if self.hull is not None and self.hull:
            if not in_polygon(x, y, self.hull):
                return None

        sigma2 = self.sigma * self.sigma
        cutoff2 = self.cutoff * self.cutoff
        in_range = False
        wv_sum = 0.
        for atom_id, v in self.values.items():
            ax, ay = self.plane[atom_id]
            dx = x - ax
            dy = y - ay
            d2 = dx * dx + dy * dy
            if d2 > cutoff2:
                continue
            in_range = True
            w = exp(-d2 / (2. * sigma2))
            wv_sum += w * v

        if not in_range:
            return None
        # The guard is "no atom was within cutoff", not a threshold on the weight sum.
        return wv_sum

    def gradient(self, x: float, y: float) -> tuple[float, float]:
        """Analytic gradient of the un-normalized Gaussian field at (x, y).

            ∂f/∂xⱼ = Σᵢ vᵢ · (∂wᵢ/∂xⱼ),  ∂wᵢ/∂xⱼ = −wᵢ·(xⱼ−aᵢⱼ)/σ²
        """
        sigma2 = self.sigma * self.sigma
        cutoff2 = self.cutoff * self.cutoff
        in_range = False
        gx = 0.
        gy = 0.

        for atom_id, v in self.values.items():
            ax, ay = self.plane[atom_id]
            ddx = x - ax
            ddy = y - ay
            d2 = ddx * ddx + ddy * ddy
            if d2 > cutoff2:
                continue
            in_range = True
            w = exp(-d2 / (2. * sigma2))
            coeff = -w / sigma2
            gx += coeff * ddx * v
            gy += coeff * ddy * v

        if not in_range:
            return 0., 0.
        return gx, gy

    def bounds(self, pad: float) -> Box:
        """Axis-aligned bounding box of all atoms in `values`, expanded by `pad`.

        :raises ValueError: `pad` would invert the box.  An inverted box does not propagate harmlessly:
            `sample` clamps it to a 2×2 grid whose corners lie on the wrong sides of each other, which
            marching squares reads as case 15, so the mistake would surface far away as "no contours".
        """
        xs = [self.plane[i][0] for i in self.values]
        ys = [self.plane[i][1] for i in self.values]
        min_x, max_x = min(xs) - pad, max(xs) + pad
        min_y, max_y = min(ys) - pad, max(ys) + pad
        if min_x > max_x or min_y > max_y:
            raise ValueError(
                f'pad {pad} inverts the bounding box of these atoms, which spans '
                f'{max(xs) - min(xs):.4g} × {max(ys) - min(ys):.4g}: a negative pad may not exceed '
                f'half the smaller side')
        return Box(min_x, min_y, max_x, max_y)


def sample(field: ScalarField, box: Box, spacing: float) -> Grid:
    """Sample `field` on a regular grid that covers `box`.

    The grid is row-major: z[j * nx + i] is the sample at column i, row j.
    None is stored where `field.at()` returns None.
    """
    nx = max(2, int((box.max_x - box.min_x) / spacing) + 1)
    ny = max(2, int((box.max_y - box.min_y) / spacing) + 1)
    z: list[float | None] = []
    for j in range(ny):
        y = box.min_y + j * spacing
        for i in range(nx):
            x = box.min_x + i * spacing
            z.append(field.at(x, y))
    return Grid(box.min_x, box.min_y, spacing, nx, ny, z)


def _ms_edge_point(i: int, j: int, spacing: float, min_x: float, min_y: float,
                   v00: float, v10: float, v01: float, v11: float,
                   level: float, edge: int) -> tuple[float, float]:
    """World coordinate of a level crossing on edge `edge` of cell (i, j).

    Cell corners: BL=(i,j), BR=(i+1,j), TR=(i+1,j+1), TL=(i,j+1).
    Edges: 0=bottom (BL→BR), 1=right (BR→TR), 2=top (TL→TR), 3=left (BL→TL).
    """
    x0 = min_x + i * spacing
    y0 = min_y + j * spacing
    s = spacing
    if edge == 0:   # bottom: BL → BR
        t = (level - v00) / (v10 - v00) if v10 != v00 else 0.5
        return x0 + t * s, y0
    elif edge == 1:  # right: BR → TR
        t = (level - v10) / (v11 - v10) if v11 != v10 else 0.5
        return x0 + s, y0 + t * s
    elif edge == 2:  # top: TL → TR  (note: traced left-to-right so TL first)
        t = (level - v01) / (v11 - v01) if v11 != v01 else 0.5
        return x0 + t * s, y0 + s
    else:           # left: BL → TL
        t = (level - v00) / (v01 - v00) if v01 != v00 else 0.5
        return x0, y0 + t * s


def _build_ms_table():
    table = []
    for case in range(16):
        bl = bool(case & 1)
        br = bool(case & 2)
        tr = bool(case & 4)
        tl = bool(case & 8)
        crossings = []
        if bl != br:
            crossings.append(0)
        if br != tr:
            crossings.append(1)
        if tr != tl:
            crossings.append(2)
        if tl != bl:
            crossings.append(3)
        if len(crossings) == 0:
            table.append([])
        elif len(crossings) == 2:
            table.append([(crossings[0], crossings[1])])
        elif len(crossings) == 4:
            table.append(None)  # saddle: cases 5 and 10
        else:
            table.append([])  # degenerate
    return table


_MSTABLE = _build_ms_table()


def _saddle_segments(case: int, centre: float, level: float) -> list[tuple[int, int]]:
    """Resolve saddle ambiguity by centre value.

    Each edge pair that shares exactly one corner cuts that corner off from the rest of the cell:
    {0,1} cuts off BR, {1,2} TR, {2,3} TL, {0,3} BL.  The centre says which diagonal is connected through
    the middle, so the other diagonal's two corners are isolated and the two segments must cut those off.
    Case 5 is BL and TR above the level, case 10 BR and TL, which is why they answer with opposite pairs.
    """
    if case == 5:
        if centre >= level:
            return [(0, 1), (2, 3)]
        else:
            return [(0, 3), (1, 2)]
    else:  # case 10
        if centre >= level:
            return [(0, 3), (1, 2)]
        else:
            return [(0, 1), (2, 3)]


def _round_pt(x: float, y: float) -> tuple[int, int]:
    """Round a world point to an integer key at 1e-9 precision."""
    return round(x * 1e9), round(y * 1e9)


def isolines(grid: Grid, level: float) -> list[list[tuple[float, float]]]:
    """Trace isolines for `level` from `grid` using marching squares.

    Returns a list of polylines.  Closed polylines repeat their first point as the last; open ones
    (running off the grid edge / cutoff boundary) do not.  A cell with a None corner is skipped — the
    field is undefined there, and emitting a segment would invent a wall at the cutoff boundary.  That
    skip is observable: a level at or near zero grazes the None ring and comes back as many open
    fragments of that ragged edge rather than as one ring.
    """
    nx, ny = grid.nx, grid.ny
    s = grid.spacing
    mx, my = grid.min_x, grid.min_y

    # Collect raw segments, then chain them into polylines through endpoint adjacency.
    raw_segments: list[tuple[tuple[float, float], tuple[float, float]]] = []

    for j in range(ny - 1):
        for i in range(nx - 1):
            v00 = grid.z[j * nx + i]
            v10 = grid.z[j * nx + (i + 1)]
            v11 = grid.z[(j + 1) * nx + (i + 1)]
            v01 = grid.z[(j + 1) * nx + i]

            if v00 is None or v10 is None or v11 is None or v01 is None:
                continue

            bl = v00 >= level
            br = v10 >= level
            tr = v11 >= level
            tl = v01 >= level

            case = (int(bl)) | (int(br) << 1) | (int(tr) << 2) | (int(tl) << 3)

            if case == 0 or case == 15:
                continue

            if case == 5 or case == 10:
                centre = (v00 + v10 + v11 + v01) * 0.25
                pairs = _saddle_segments(case, centre, level)
            else:
                pairs = _MSTABLE[case]
                if not pairs:
                    continue

            for edge_a, edge_b in pairs:
                pa = _ms_edge_point(i, j, s, mx, my, v00, v10, v01, v11, level, edge_a)
                pb = _ms_edge_point(i, j, s, mx, my, v00, v10, v01, v11, level, edge_b)
                raw_segments.append((pa, pb))

    if not raw_segments:
        return []

    # key -> [(neighbour_key, my_xy, neighbour_xy)]; edge identity is the frozenset of the two keys
    from collections import defaultdict
    adj: dict[tuple[int, int],
              list[tuple[tuple[int, int], tuple[float, float], tuple[float, float]]]] = defaultdict(list)
    for pa, pb in raw_segments:
        ka = _round_pt(*pa)
        kb = _round_pt(*pb)
        adj[ka].append((kb, pa, pb))
        adj[kb].append((ka, pb, pa))

    used_edges: set = set()   # frozenset of two keys
    chains: list[list[tuple[float, float]]] = []

    def _walk(start_k, start_xy):
        chain = [start_xy]
        cur_k = start_k
        while True:
            found = False
            for (nk, my_xy, their_xy) in adj[cur_k]:
                edge = frozenset((cur_k, nk))
                if edge not in used_edges:
                    used_edges.add(edge)
                    chain.append(their_xy)
                    cur_k = nk
                    found = True
                    break
            if not found:
                break
        return chain, cur_k

    # free ends (degree 1) first, so no chain is started from its middle
    all_keys = set(adj.keys())
    free_ends = [k for k in all_keys if len(adj[k]) == 1]
    loop_starts = [k for k in all_keys if len(adj[k]) == 2]

    processed_starts: set = set()

    # open chains
    for start_k in free_ends:
        if start_k in processed_starts:
            continue
        _, start_xy, _ = adj[start_k][0]
        chain, end_k = _walk(start_k, start_xy)
        if len(chain) > 1:
            chains.append(chain)
        processed_starts.add(start_k)

    # closed loops
    for start_k in loop_starts:
        any_unused = any(frozenset((start_k, nk)) not in used_edges
                         for (nk, _, _) in adj[start_k])
        if not any_unused:
            continue
        _, start_xy, _ = adj[start_k][0]
        chain, end_k = _walk(start_k, start_xy)
        if len(chain) > 1:
            chain.append(chain[0])
            chains.append(chain)

    return chains


#: How far below a chain's OWN strongest gradient a vertex may sit and still be on a contour.  A ratio, so
#: it is free of the field's units and of the values' magnitude.  Measured on a phenol charge field at the
#: shipped defaults: |grad f| varies by a factor of 3.6 along a genuine contour and reaches 7e-17 out on
#: the asymptote, so one millionth sits in that gap with ten orders of clearance on the noise side.
_GRADIENT_COLLAPSE = 1e-6


def trim_asymptote(field: ScalarField,
                   polyline: list[tuple[float, float]]) -> list[list[tuple[float, float]]]:
    """The runs of `polyline` that lie on a real contour: the asymptotic tail is dropped, not the chain.

    Where the sum has decayed to nothing the sign of `f` is floating-point residue rather than field, so
    marching squares traces a ragged curve that is not a contour.  One chain can be both -- a nodal line
    where two Gaussians cancel (|grad f| ~ 1) near the atoms, asymptote (~1e-17) far out -- so collapse is
    measured against the chain's own strongest gradient and never as a threshold on `f`, small on either.

    :return: the maximal runs of surviving vertices, each of at least two, in the chain's own order.
        Nothing dropped gives the whole chain as one run, so a caller's "is this closed" test still holds.
    """
    if len(polyline) < 2:
        # One vertex bounds nothing and strokes nothing; `to_cubics` would emit a bare move.
        return []
    magnitudes = [hypot(*field.gradient(x, y)) for x, y in polyline]
    strongest = max(magnitudes)
    if strongest <= 0.:
        # The field does not vary anywhere on this chain: the degenerate case of the asymptote, and the
        # one the ratio cannot express, so it is answered here rather than as 0/0.
        return []

    floor = strongest * _GRADIENT_COLLAPSE
    if all(g >= floor for g in magnitudes):
        return [list(polyline)]

    runs: list[list[tuple[float, float]]] = []
    run: list[tuple[float, float]] = []
    for point, g in zip(polyline, magnitudes):
        if g >= floor:
            run.append(point)
            continue
        if len(run) >= 2:
            runs.append(run)
        run = []
    if len(run) >= 2:
        runs.append(run)
    return runs


def refine(field: ScalarField, polyline: list[tuple[float, float]],
           level: float, steps: int, *, spacing: float) -> list[tuple[float, float]]:
    """Newton-refine each vertex of `polyline` onto the true isoline at `level`.

    Each step: p ← p - (f(p) - level) * ∇f / |∇f|².  Vertices where |∇f| < 1e-9 or where `at()` returns
    None are left unchanged.

    `spacing` is the grid spacing the polyline was traced on, and it is the divergence guard: a
    marching-squares vertex is within one cell of the true isoline by construction, so a step longer than
    one cell is not refinement -- the vertex keeps its traced position and the remaining steps with it.
    """
    limit2 = spacing * spacing
    result = []
    for x, y in polyline:
        rx, ry = x, y
        for _ in range(steps):
            fval = field.at(rx, ry)
            if fval is None:
                break
            gx, gy = field.gradient(rx, ry)
            g2 = gx * gx + gy * gy
            if g2 < 1e-18:  # |∇f| < 1e-9
                break
            step = (fval - level) / g2
            dx = -step * gx
            dy = -step * gy
            if dx * dx + dy * dy > limit2:
                rx, ry = x, y   # diverging: keep the marching-squares point
                break
            rx += dx
            ry += dy
        result.append((rx, ry))
    return result


def to_cubics(field: ScalarField, polyline: list[tuple[float, float]],
              level: float, *, closed: bool) -> tuple:
    """Fit one cubic Bezier per span.  Tangents come from the field gradient rotated 90°.

    At each vertex the tangent direction is (-gy, gx) normalised, with its sign chosen
    to point in the direction of travel (chord to next vertex).  Each span uses handles
    at chord_length / 3 along the tangent from each endpoint.

    Returns a tuple of scene segments: (M, ...), (C, ...), ..., optionally (Z,).
    """
    pts = list(polyline)
    if not pts:
        return ()

    # For a closed polyline the last point repeats the first — drop the duplicate
    if closed and len(pts) >= 2 and pts[0] == pts[-1]:
        pts = pts[:-1]

    n = len(pts)
    if n < 2:
        # Single point: emit just a move
        return (move(*pts[0]),)

    def _tangent_at(idx: int, pts: list, closed: bool) -> tuple[float, float]:
        """Unit tangent at pts[idx], perpendicular to ∇f, signed to follow travel direction."""
        x, y = pts[idx]
        gx, gy = field.gradient(x, y)
        g2 = gx * gx + gy * gy
        if g2 < 1e-18:
            # Gradient vanishes: fall back to chord direction
            if closed:
                next_idx = (idx + 1) % n
                prev_idx = (idx - 1) % n
            else:
                next_idx = min(idx + 1, n - 1)
                prev_idx = max(idx - 1, 0)
            nx_ = pts[next_idx][0] - pts[prev_idx][0]
            ny_ = pts[next_idx][1] - pts[prev_idx][1]
            nn = sqrt(nx_ * nx_ + ny_ * ny_)
            if nn < 1e-18:
                return 1., 0.
            return nx_ / nn, ny_ / nn
        # Isoline tangent: gradient rotated 90°
        g = sqrt(g2)
        tx, ty = -gy / g, gx / g
        if closed:
            next_idx = (idx + 1) % n
        else:
            next_idx = min(idx + 1, n - 1)
        cx = pts[next_idx][0] - x
        cy = pts[next_idx][1] - y
        if tx * cx + ty * cy < 0.:
            tx, ty = -tx, -ty
        return tx, ty

    tangents = [_tangent_at(i, pts, closed) for i in range(n)]

    segs = [move(*pts[0])]

    spans = n if closed else n - 1
    for k in range(spans):
        i0 = k
        i1 = (k + 1) % n if closed else k + 1
        x0, y0 = pts[i0]
        x1, y1 = pts[i1]
        chord = sqrt((x1 - x0) ** 2 + (y1 - y0) ** 2)
        if chord < 1e-9:
            # A marching-squares ring occasionally produces a near-zero-length closing segment when the
            # trace starts next to its wrap point.  Skip it — Z already closes back to the M point.
            continue
        h = chord / 3.
        t0x, t0y = tangents[i0]   # forward tangent at i0: cp1 = pts[i0] + h * t0
        t1x, t1y = tangents[i1]   # forward tangent at i1: cp2 = pts[i1] - h * t1
        # G1 continuity: cp2 of THIS span and cp1 of the NEXT are `pts[i1] ∓ h*t1`, exactly collinear and
        # opposite, so no conditional flip is needed.
        cp1x = x0 + h * t0x
        cp1y = y0 + h * t0y
        cp2x = x1 - h * t1x    # backward handle: always -tangent
        cp2y = y1 - h * t1y
        segs.append(curve(cp1x, cp1y, cp2x, cp2y, x1, y1))

    if closed:
        segs.append(close())

    return tuple(segs)


def contour_levels(count: int, vmin: float, vmax: float) -> list[float]:
    """Return `count` evenly-spaced levels strictly inside [vmin, vmax].

    Formula: `vmin + (i + 1) * (vmax - vmin) / (count + 1)` for i in range(count).  Interior is not
    cosmetic: a level at or just above the sampled minimum runs along the cutoff boundary, where
    `isolines`' None-corner skip chops it into open fragments.  At the one production call site
    (`overlay._field_bands`) [vmin, vmax] is the range the field was SAMPLED over and not the colormap's
    domain -- a contour can only exist between the field's own extremes.
    """
    if count <= 0:
        return []
    span = vmax - vmin
    step = span / (count + 1)
    return [vmin + (i + 1) * step for i in range(count)]


def convex_hull(points: list[tuple[float, float]], pad: float) -> tuple[tuple[float, float], ...]:
    """Monotone-chain convex hull of `points`, each edge offset outward by `pad`.

    Degenerate cases: 0 points → empty tuple; 1 point → axis-aligned square of half-width `pad`;
    2 points or a collinear set → rectangle of half-width `pad` around the segment.
    """
    pts = list(points)
    if not pts:
        return ()

    if len(pts) == 1:
        x, y = pts[0]
        return ((x - pad, y - pad), (x + pad, y - pad),
                (x + pad, y + pad), (x - pad, y + pad))

    # Andrew's monotone chain
    pts_sorted = sorted(set(pts))

    def cross(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    lower = []
    for p in pts_sorted:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)

    upper = []
    for p in reversed(pts_sorted):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)

    hull_pts = lower[:-1] + upper[:-1]

    if len(hull_pts) < 3:
        # Collinear: a rectangle around the segment between the two extreme points
        p0, p1 = pts_sorted[0], pts_sorted[-1]
        dx = p1[0] - p0[0]
        dy = p1[1] - p0[1]
        length = sqrt(dx * dx + dy * dy)
        if length < 1e-12:
            # All same point
            x, y = p0
            return ((x - pad, y - pad), (x + pad, y - pad),
                    (x + pad, y + pad), (x - pad, y + pad))
        # Unit perpendicular
        px = -dy / length * pad
        py = dx / length * pad
        # Also extend endpoints by pad along segment direction
        ux = dx / length * pad
        uy = dy / length * pad
        return (
            (p0[0] - ux - px, p0[1] - uy - py),
            (p1[0] + ux - px, p1[1] + uy - py),
            (p1[0] + ux + px, p1[1] + uy + py),
            (p0[0] - ux + px, p0[1] - uy + py),
        )

    if pad == 0.:
        return tuple(hull_pts)

    # Offset each edge outward by `pad`, then recompute vertices as edge intersections.  "Outward" for a
    # counter-clockwise hull means to the right of the edge direction.
    n = len(hull_pts)
    offset_edges = []
    for k in range(n):
        p0 = hull_pts[k]
        p1 = hull_pts[(k + 1) % n]
        dx = p1[0] - p0[0]
        dy = p1[1] - p0[1]
        length = sqrt(dx * dx + dy * dy)
        if length < 1e-12:
            continue
        # Outward normal: (dy, -dx) / length
        nx_ = dy / length
        ny_ = -dx / length
        a = (p0[0] + pad * nx_, p0[1] + pad * ny_)
        b = (p1[0] + pad * nx_, p1[1] + pad * ny_)
        offset_edges.append((a, b))

    if not offset_edges:
        return tuple(hull_pts)

    # Recompute vertices as intersections of consecutive offset edges
    new_verts = []
    m = len(offset_edges)
    for k in range(m):
        a1, b1 = offset_edges[k]
        a2, b2 = offset_edges[(k + 1) % m]
        # Intersect lines a1+t*(b1-a1) and a2+s*(b2-a2)
        d1x = b1[0] - a1[0]
        d1y = b1[1] - a1[1]
        d2x = b2[0] - a2[0]
        d2y = b2[1] - a2[1]
        denom = d1x * d2y - d1y * d2x
        if abs(denom) < 1e-12:
            # Parallel edges (very short edge): use midpoint of endpoint pair
            new_verts.append(((b1[0] + a2[0]) * 0.5, (b1[1] + a2[1]) * 0.5))
        else:
            t = ((a2[0] - a1[0]) * d2y - (a2[1] - a1[1]) * d2x) / denom
            x = a1[0] + t * d1x
            y = a1[1] + t * d1y
            new_verts.append((x, y))

    return tuple(new_verts)


def in_polygon(x: float, y: float, polygon: tuple[tuple[float, float], ...]) -> bool:
    """Crossing-number test: True if (x, y) is inside `polygon`.

    Uses the half-open edge convention `(y1 <= y) != (y2 <= y)` so a vertex exactly
    on the ray is counted once rather than twice.
    """
    n = len(polygon)
    if n < 3:
        return False
    crossings = 0
    for k in range(n):
        x1, y1 = polygon[k]
        x2, y2 = polygon[(k + 1) % n]
        if (y1 <= y) != (y2 <= y):
            x_cross = x1 + (y - y1) * (x2 - x1) / (y2 - y1)
            if x < x_cross:
                crossings += 1
    return crossings % 2 == 1
