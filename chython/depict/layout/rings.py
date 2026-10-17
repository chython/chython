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
"""Ring blocks on the lattice, by ear decomposition.

| step | rule |
| --- | --- |
| block | SSSR rings sharing two or more atoms, merged; a spiro atom joins two blocks |
| first ring | the ring with most fused neighbours, six first; drawn as a fixed shape |
| small shapes | 3 lattice triangle, 4 rhombus, 5 hexagon less its top vertex, 6 hexagon, 7 and 8 bulged hexagons |
| large shapes | the boundary of the roundest polyhex of that perimeter, bulged or cut by one for an odd size |
| ear | the unplaced run of the ring with most placed atoms; a best-first walk of exactly its bond count |
| moves | unit step, chord of `sqrt 3`, straight 2 and 3 steps; the longer, the dearer |
| cost | turn away from 120 degrees, a hexagon-centre point, a small ring inside a placed ring |
| ranking | conflicts, broken stated cis/trans bonds, exits with no free unit step, cost |
| first ring | turned on its shape to honour its cis/trans bonds, then, alone in its block, to put exits outward |
| | an atom of two exits also keeps both ring neighbours off 240-degree corners |
| bridge | one atom across a hexagon's para atoms: the rest on the lattice, the atom `LIFT` off the centre |
| round | a macrocycle drawn as a ribbon, isoperimetric ratio under `ROUND`, is walked again on each polyhex of up
  to three rows, roundest first; the `TURNS` cheapest clean turns of each go to the engine |
| regrid | a macrocycle whose small-ring unit walked off its own drawing: each unit is placed as drawn alone |
| | at the lattice symmetry and shift nearest the walk; chain atoms over `KEEP` bonds away keep their points |
| | the rest rewalk; kept when conflicts, broken bonds and crowding do not grow |

Every choice is ordered by `rank`, the canonical position of an atom, so the embedding ignores atom order.
"""
from functools import lru_cache
from heapq import heappop, heappush
from math import atan2, cos, degrees, fsum, hypot, pi, sin

from ...core.wedge import SU_CIS_TRANS, cis_trans_frame
from .lattice import C30, DIRS, L, ON_BOND, SYMMETRIES, Board, _crosses, _seg_dist, add, apply, project, sub, sublattice


__all__ = ['ring_blocks', 'bridge_block', 'embed_block', 'macro_block', 'planar_block', 'round_block', 'roundness',
           'sphere_block', 'cage_block']

MOVES = tuple([(d, 0.) for d in DIRS] +
              [(add(DIRS[k], DIRS[(k + 1) % 6]), 2.) for k in range(6)] +
              [((2 * d[0], 2 * d[1]), 3.) for d in DIRS] +
              [((3 * d[0], 3 * d[1]), 6.) for d in DIRS])
TURN = 1.5          # per 60 degrees away from 120
CENTRE = 1.       # a point on the block's hexagon-centre sublattice
INSIDE = 8.       # a point of a small ring's ear inside a placed ring
EXPANSIONS = 40000   # per ear route
GREED = 4.           # a long ear's second route: priced per move left, so it dives for any closing walk
RETRIES = 160000     # route expansions over all retries of one block
TUTTE_REDRAWN = 3    # best planar drawings redrawn with length weights
TUTTE_SPREAD = .5    # least distance of two atoms of a planar drawing, in L: adamantane outlined at L
BEND = (.5, .4, .3, .2)  # offsets tried for a bent atom, in median bonds
CLEAR = .3           # least distance of a bent atom from a bond it is not in, in median bonds
SPHERE = 30          # fewest atoms of a cage drawn as a projected sphere; a smaller one stays planar
BRIDGE = .02         # a saturated macrocycle bridge's give against a unit exit's 1: it holds its 120 degrees
EXIT = .8            # least distance of a macrocycle's substituent exit from another exit or block atom, in bonds
EYE = 3.             # the perspective eye's distance from a sphere's centre, in radii
SPHERE_STEPS = 400   # subspace iterations of a sphere's coordinates
CAGE = 12            # fewest atoms of a dense cage drawn as a projected model
DENSE = .4           # least rings per atom of a dense cage: dodecahedrane .55, adamantane .3
KNOT = .25           # least rings per atom of a cage whose flat drawing stretches a bond to `STRAIN` times another
STRAIN = 2.5         # bond length ratio of a flat drawing that a `KNOT`-dense cage trades for a model
CAGE_STEPS = 300     # stress-majorization iterations of a dense cage's coordinates
VIEWS = 400          # view directions tried for a dense cage, spread evenly over the sphere
KEEP = 5             # bonds from a regridded unit beyond which a macrocycle chain atom keeps its walked point
SNAP = .02           # squared residual per atom, in squared L, below which a walked unit is not bent
ROUND = .35          # isoperimetric ratio of a macrocycle polygon below which it is a ribbon
TURNS = 8            # cheapest clean turns of one polyhex a ribbon macrocycle is redrawn with
LIFT = .2            # a one-atom bridge's offset from the centre of the hexagon it spans, in bonds


def ring_blocks(mol):
    """`[(atoms, rings)]`: frozenset of atoms and the SSSR rings of each block, the blocks in order of their sorted
    ranks, since two blocks on one spiro atom share their lowest rank; each ring in cyclic order from its lowest rank
    towards its lower-ranked neighbour."""
    groups = []
    for r in mol.sssr:
        s, rs = set(r), [tuple(r)]
        while True:
            hit = [g for g in groups if len(g[0] & s) >= 2]
            if not hit:
                break
            for g in hit:
                groups.remove(g)
                s |= g[0]
                rs.extend(g[1])
        groups.append((s, rs))
    rank = mol.canonical_order()
    groups.sort(key=lambda g: sorted(rank[x] for x in g[0]))
    return [(frozenset(s), [_rotated(r, rank) for r in _canonical_basis(mol, s, rs, rank)]) for s, rs in groups]


def _rotated(ring, rank):
    """`ring` from its lowest-ranked atom towards the lower-ranked of that atom's two ring neighbours."""
    i = min(range(len(ring)), key=lambda i: rank[ring[i]])
    seq = tuple(ring[i:]) + tuple(ring[:i])
    return seq[:1] + seq[1:][::-1] if rank[seq[-1]] < rank[seq[1]] else seq


def _canonical_basis(mol, atoms, rings, rank, budget=50000):
    """A minimum cycle basis of the block chosen in rank order, so that a cage's rings ignore atom order; the
    given `rings` when the cycle enumeration exceeds `budget` steps."""
    if len(rings) < 2:
        return rings
    top = max(map(len, rings))
    cycles = _cycles(mol, atoms, rank, top, budget, whole=True)
    if cycles is None:
        return rings
    cycles.sort(key=lambda c: (len(c), sorted(rank[x] for x in c)))
    index = {}
    basis = []                                              # GF(2) row echelon: {pivot: bitset}
    out = []
    for c in cycles:
        v = 0
        for i in range(len(c)):
            e = frozenset((c[i], c[i - 1]))
            v ^= 1 << index.setdefault(e, len(index))
        for piv, row in basis:
            if v >> piv & 1:
                v ^= row
        if v:
            basis.append((v.bit_length() - 1, v))
            out.append(c)
            if len(out) == len(rings):
                return out
    return rings


# --- shapes -------------------------------------------------------------------------------------------------

def _hexagon(c=(0, 0)):
    return [add(c, d) for d in DIRS]


def _boundary(cells):
    """The boundary cycle of a polyhex as points in order, or None when it is not one simple cycle."""
    edges = {}
    for c in cells:
        vs = _hexagon(c)
        for k in range(6):
            e = frozenset((vs[k], vs[(k + 1) % 6]))
            edges[e] = edges.get(e, 0) + 1
    out = [tuple(e) for e, n in edges.items() if n == 1]
    adj = {}
    for u, v in out:
        adj.setdefault(u, []).append(v)
        adj.setdefault(v, []).append(u)
    if any(len(v) != 2 for v in adj.values()):
        return None
    start = min(adj)
    cycle, prev, cur = [start], None, start
    while True:
        a, b = adj[cur]
        nxt = b if a == prev else a
        if nxt == start:
            break
        cycle.append(nxt)
        prev, cur = cur, nxt
    return cycle if len(cycle) == len(adj) else None


def _cell_neighbours(c):
    # hexagon centres adjacent to `c` across an edge: two unit steps summed
    return [add(c, add(DIRS[k], DIRS[(k + 1) % 6])) for k in range(6)]


@lru_cache(maxsize=None)
def _polyhexes():
    """{perimeter: cells} for the roundest simply connected polyhex found of each even perimeter."""
    best = {}

    def keep(cells):
        b = _boundary(cells)
        if b is None:
            return
        n = len(b)
        if n not in best or len(cells) > len(best[n]):
            best[n] = tuple(cells)
    # spiral prefixes: cells by distance from the origin cell
    grid = {(0, 0)}
    frontier = [(0, 0)]
    for _ in range(4):
        nxt = []
        for c in frontier:
            for d in _cell_neighbours(c):
                if d not in grid:
                    grid.add(d)
                    nxt.append(d)
        frontier = nxt
    order = sorted(grid, key=lambda c: (round(hypot(*project(c)), 6), atan2(project(c)[1], project(c)[0])))
    for k in range(1, len(order) + 1):
        keep(order[:k])
    # ribbons of one and two rows
    step = add(DIRS[5], DIRS[0])                            # horizontal, to the right-hand neighbour cell
    down = add(DIRS[4], DIRS[5])
    for a in range(1, 16):
        row = [(step[0] * i, step[1] * i) for i in range(a)]
        keep(row)
        for b in range(1, a + 1):
            keep(row + [add(down, (step[0] * i, step[1] * i)) for i in range(b)])
    return best


def _shape(n):
    """The points of a first ring of size `n`, in cyclic order."""
    h = _hexagon()
    if n == 3:
        return [(0, 0), DIRS[0], DIRS[1]]
    if n == 4:
        return [(0, 0), DIRS[0], add(DIRS[0], DIRS[1]), DIRS[1]]
    if n == 5:
        return [h[2], h[3], h[4], h[5], h[0]]             # top vertex dropped: the chord h0-h2 is horizontal
    if n == 6:
        return h
    if n == 7:
        return _bulge(h, [(5, 0)])
    if n == 8:
        return _bulge(h, [(5, 0), (2, 3)])
    table = _polyhexes()
    if n in table:
        return _boundary(table[n])
    if n - 1 in table:
        return _bulge(_boundary(table[n - 1]), None)
    if n + 1 in table:
        b = _boundary(table[n + 1])
        return b[1:]                                        # a chord over the dropped corner
    return None


def _bulge(cycle, edges):
    """`cycle` with each edge `(i, i + 1)` replaced by a two-step detour through the cell centre beyond it;
    `None` picks the rightmost edge."""
    n = len(cycle)
    if edges is None:
        i = max(range(n), key=lambda i: (project(cycle[i])[0] + project(cycle[(i + 1) % n])[0],
                                         project(cycle[i])[1]))
        edges = [(i, (i + 1) % n)]
    cx = sum(project(p)[0] for p in cycle) / n
    cy = sum(project(p)[1] for p in cycle) / n
    out = list(cycle)
    for i, j in sorted(edges, reverse=True):
        u, v = cycle[i], cycle[j]
        cands = [add(u, d) for d in DIRS if sub(add(u, d), v) in DIRS]
        c = max(cands, key=lambda p: (project(p)[0] - cx) ** 2 + (project(p)[1] - cy) ** 2)
        out.insert(i + 1, c)
    return out


# --- ears ---------------------------------------------------------------------------------------------------

def _angle(u, v):
    a = degrees(atan2(*reversed(project(u)))) - degrees(atan2(*reversed(project(v))))
    return abs((a + 180) % 360 - 180)


@lru_cache(maxsize=None)
def _turn(u, v):
    """Cost of the angle between bonds `u` and `v` leaving one atom."""
    return abs(_angle(u, v) - 120) / 60 * TURN


def _gap(board, mol, x, v, placed):
    """Cost of a new bond leaving `x` along `v`, against the bonds `x` already has."""
    p = board.pos[x]
    worst = 0.
    for y in mol.neighbors_of(x):
        if y in placed and y in board.pos and _bonded(board, x, y):
            a = _angle(v, sub(board.pos[y], p))
            worst = max(worst, max(0., 120 - a) / 60 * TURN)
    return worst


def _bonded(board, x, y):
    return ((x, y) if x < y else (y, x)) in board.bonds


def _inside(pt, polys):
    x, y = project(pt)
    for poly in polys:
        inside = False
        n = len(poly)
        for i in range(n):
            (x1, y1), (x2, y2) = poly[i], poly[(i + 1) % n]
            if (y1 > y) != (y2 > y) and x < (x2 - x1) * (y - y1) / (y2 - y1) + x1:
                inside = not inside
        if inside:
            return True
    return False


STRETCH = 2 / (3 ** .5 - 1) / L     # least price per unit of length beyond a unit step: the chord's
EXACT = 4                           # moves left up to which `_bound` is exact rather than a length bound


@lru_cache(maxsize=None)
def _least(r, d, u):
    """Least price and turn cost of exactly `r` moves summing to the step `d` after the move `u`, obstacles and the
    turn at the far end aside; inf when none does."""
    if not r:
        return 0. if d == (0, 0) else float('inf')
    if hypot(*project(d)) > 3 * r * L + 1e-9:
        return float('inf')
    back = (-u[0], -u[1])
    return min(price + _turn(back, v) + _least(r - 1, sub(d, v), v) for v, price in MOVES)


def _bound(r, d, u):
    """A lower bound on the cost of the last `r` moves of a walk still `d` from its end after the move `u`: exact
    for a short tail."""
    if r <= EXACT:
        return _least(r, d, u)
    reach = hypot(*project(d))
    return float('inf') if reach > 3 * r * L else max(0., reach - r * L) * STRETCH


REACH = (3 * L + ON_BOND) ** 2      # squared: beyond it from a move's start, a point is off the move by ON_BOND
CROSS = (6 * L) ** 2                # squared: beyond it between their starts, two moves cannot cross


def _near(v):
    """The lattice points other than its ends within `ON_BOND` of the move `v`, as offsets from its start: the points
    inside a straight move, and the two flanking a diagonal one. Every other lattice point is `sqrt(3) / 2` of a bond
    away or farther."""
    t = project(v)
    return tuple((q, r) for q in range(-4, 5) for r in range(-4, 5)
                 if (q, r) not in ((0, 0), v) and _seg_dist(project((q, r)), (0., 0.), t) < ON_BOND)


NEAR = tuple(_near(v) for v, _ in MOVES)


def _near_path(xq, xprev):
    """Whether `xq` is within `ON_BOND` of a move of the projected walk `xprev`."""
    px, py = xq
    for (sx, sy), t in zip(xprev, xprev[1:]):
        if (sx - px) ** 2 + (sy - py) ** 2 <= REACH and _seg_dist(xq, (sx, sy), t) < ON_BOND:
            return True
    return False


def _near_bond(xp, xq, xs):
    """Whether a projected point of `xs` is within `ON_BOND` of the move from `xp` to `xq`."""
    px, py = xp
    for x in xs:
        if (x[0] - px) ** 2 + (x[1] - py) ** 2 <= REACH and _seg_dist(x, xp, xq) < ON_BOND:
            return True
    return False


def _crosses_path(p, q, xp, prev, xprev):
    """Whether the move from `p` to `q` properly crosses a move of the walk `prev`, projected `xprev`."""
    px, py = xp
    for j in range(len(prev) - 1):
        s, t = prev[j], prev[j + 1]
        if s == p or t == p:
            continue
        sx, sy = xprev[j]
        if (sx - px) ** 2 + (sy - py) ** 2 <= CROSS and _crosses(p, q, s, t):
            return True
    return False


def _route(board, mol, a, b, m, centre, polys, placed, spent, greed=0.):
    """`(points, cost)` of the `m` new atoms of an ear from `a` to `b`, the cheapest walk of `m + 1` moves, or
    None; best-first by cost plus `_bound` of the moves left, plus `greed` per move left. Its expansions are added to
    `spent[0]`, and it stops when that reaches `spent[1]`."""
    pa, pb = board.pos[a], board.pos[b]
    heap = [(0., 0, 0., 0, pa, 0, None, (), (project(pa),), ())]
    # an ear from a lattice point walks the lattice, where `NEAR` answers ON_BOND exactly; `inner` holds the lattice
    # points near the walk's moves
    grid = all(type(v) is int for v in pa)
    seen = set()
    count = 0
    limit = min(EXPANSIONS, spent[1] - spent[0])
    while heap and count < limit:
        _, _, cost, _, p, k, last, path, xpath, inner = heappop(heap)
        key = (p, k, last)
        if key in seen:
            continue
        seen.add(key)
        if k == m + 1:
            spent[0] += count
            return list(path[:-1]), cost
        left = m - k                                        # moves left after the next one
        xp = xpath[-1]
        for i, (v, price) in enumerate(MOVES):
            q = add(p, v)
            xq = project(q)
            if left == 0:
                if q != pb:
                    continue
            else:
                if q == pb or q in path:
                    continue
                h = _bound(left, sub(pb, q), v)
                if h == float('inf'):
                    continue
                if board.bad_atom(q):
                    continue
                if (q == pa or q in inner) if grid else _near_path(xq, xpath):
                    continue
            if board.bad_bond(p, q):
                continue
            if any(add(p, o) in path for o in NEAR[i]) if grid else _near_bond(xp, xq, xpath[1:-1]):
                continue
            prev = (pa,) + path
            if _crosses_path(p, q, xp, prev, xpath):
                continue
            c = cost + price
            if k == 0:
                c += _gap(board, mol, a, v, placed)
            else:
                c += _turn(sub(prev[-2], p), v)
            if left == 0:
                c += _gap(board, mol, b, sub(p, q), placed)
                h = 0.
            else:
                if centre is not None and sublattice(q) == centre:
                    c += CENTRE
                if polys and _inside(q, polys):
                    c += INSIDE
            count += 1
            # a tie goes deepest first
            heappush(heap, (c + h + greed * left, -k, c, count, q, k + 1, i, path + (q,), xpath + (xq,),
                            inner + tuple(add(p, o) for o in NEAR[i]) if grid and NEAR[i] else inner))
    spent[0] += count
    return None


def embed_block(mol, atoms, rings, rank, k=1):
    """`({atom: point}, conflicts)` for one block: a lattice embedding and the defects it could not avoid, a bond out
    of the block with no free unit step counted among them. With `k` over 1, a list of up to `k` distinct such
    results, best first, the search going on until `k` are clean or its budget is spent.

    The first cycle is the ring of most fused neighbours, or the longest simple cycle when the block holds a ring
    over 8: a long first cycle leaves short ears, whose routes are cheap. A block that conflicts retries from every
    longest simple cycle, then every SSSR ring, and keeps the fewest conflicts, then the fewest broken cis/trans
    bonds, then the fewest crowded exits, then the cheapest walk; last come the simple cycles of up to 8 atoms outside
    the basis, since the basis may hold only the cycle that puts a crowded atom inside."""
    fused = {r: sum(len(set(r) & set(s)) >= 2 for s in rings if s is not r) for r in rings}

    def key(r):
        return -fused[r], len(r) != 6, len(r), sorted(rank[x] for x in r)
    longest = [c for c in _long_cycles(mol, atoms, rank) if _shape(len(c)) is not None]
    first = longest[0] if longest and max(map(len, rings)) > 8 else min(rings, key=key)
    ct = _ring_ct(mol, atoms, rank)
    exo = {x: n for x in atoms if (n := sum(y not in atoms for y in mol.neighbors_of(x)))}
    best = _embed(mol, atoms, rings, rank, first, [0, EXPANSIONS * len(rings)], ct, exo)
    found = [] if best is None else [best]
    if best is None or any(best[1:4]) or k > 1:
        seen = {frozenset(first)}
        spent = [0, RETRIES]
        cands = longest + sorted(rings, key=key) + sorted(_short_cycles(mol, atoms, rank), key=lambda c: (
            len(c) != 6, len(c), sorted(rank[x] for x in c)))
        for c in cands:
            if frozenset(c) in seen or _shape(len(c)) is None:
                continue
            seen.add(frozenset(c))
            if spent[0] >= RETRIES:
                break
            r = _embed(mol, atoms, rings, rank, c, spent, ct, exo)
            if r is not None:
                found.append(r)
                if best is None or r[1:] < best[1:]:
                    best = r
            if best is not None and not any(best[1:4]) and sum(not any(f[1:4]) for f in found) >= k:
                break
    if k > 1:
        out, seen = [], set()
        for f in sorted(found, key=lambda f: f[1:]):
            if max(map(len, rings)) > 8:
                r = _regrid(mol, atoms, rings, rank, f[0], ct, exo)
                if r is not None and r[1:4] <= f[1:4]:
                    f = r
            key = tuple(sorted(f[0].items()))
            if key not in seen:
                seen.add(key)
                out.append((f[0], f[1] + f[3]))
            if len(out) == k:
                break
        return out
    if best is not None and max(map(len, rings)) > 8:
        r = _regrid(mol, atoms, rings, rank, best[0], ct, exo)
        if r is not None and r[1:4] <= best[1:4]:
            best = r
    return None if best is None else (best[0], best[1] + best[3])


def roundness(poly):
    """The isoperimetric ratio `4 pi area / perimeter^2` of a polygon of cartesian points: 1 for a circle."""
    n = len(poly)
    area = abs(fsum(poly[i - 1][0] * poly[i][1] - poly[i][0] * poly[i - 1][1] for i in range(n))) / 2
    per = fsum(hypot(poly[i][0] - poly[i - 1][0], poly[i][1] - poly[i - 1][1]) for i in range(n))
    return 4 * pi * area / per ** 2


@lru_cache(maxsize=None)
def _round_shapes(n):
    """`(roundness, cycle)` of `n` lattice points rounder than `ROUND`, roundest first: boundaries of polyhexes of up to
    three rows, each row starting under the one above or a cell either side, of perimeter `n` or bulged from `n - 1`."""
    step = add(DIRS[5], DIRS[0])
    down = add(DIRS[4], DIRS[5])
    found = {}
    for a in range(1, (n + 2) // 4 + 1):                  # a row of `a` cells bounds the perimeter by 4a + 2
        for b in range(a + 1):
            for c in range(b + 1):
                for o in (-1, 0, 1) if b else (0,):
                    for u in (-1, 0, 1) if c else (0,):
                        cells = [(step[0] * i, step[1] * i) for i in range(a)]
                        cells += [add(down, (step[0] * (i + o), step[1] * (i + o))) for i in range(b)]
                        cells += [add(add(down, down), (step[0] * (i + o + u), step[1] * (i + o + u)))
                                  for i in range(c)]
                        cycle = _boundary(cells)
                        if cycle is None or len(cycle) not in (n, n - 1):
                            continue
                        if len(cycle) == n - 1:
                            cycle = _bulge(cycle, None)
                        k = roundness([project(p) for p in cycle])
                        if k > ROUND:
                            found.setdefault(tuple(sorted(cycle)), (k, tuple(cycle)))
    return sorted(found.values(), key=lambda t: (-t[0], t[1]))


def round_block(mol, atoms, rings, rank):
    """Per polyhex of `_round_shapes`, `(roundness, [{atom: point}])`: the `TURNS` cheapest clean embeddings of the
    block with its largest ring laid first on it, every turn and mirror tried, the other rings walked as ears."""
    macro = max(rings, key=lambda r: (len(r), [-rank[x] for x in r]))
    ct = _ring_ct(mol, atoms, rank)
    exo = {x: n for x in atoms if (n := sum(y not in atoms for y in mol.neighbors_of(x)))}
    n = len(macro)
    turns = [macro[k:] + macro[:k] for k in range(n)]
    turns += [t[::-1] for t in turns]
    todo = [r for r in rings if r is not macro]
    for k, shape in _round_shapes(n):
        found = []
        for j, t in enumerate(turns):
            if _ct_broken(dict(zip(t, shape)), ct):
                continue
            board = Board()
            for x, p in zip(t, shape):
                board.put(x, p)
            for i in range(n):
                board.link(t[i], t[(i + 1) % n])
            placed = set(t)
            small = [[project(p) for p in shape]]
            conflicts, total = _ears(mol, board, placed, list(todo), rank, sublattice((0, 0)), small,
                                     [0, EXPANSIONS * len(rings)])
            r = _finish(mol, atoms, board, placed, conflicts, total, ct, exo)
            if r is not None and not any(r[1:4]):
                found.append((r[4], j, r[0]))
        found.sort(key=lambda f: f[:2])
        yield k, [f[2] for f in found[:TURNS]]


def _regrid(mol, atoms, rings, rank, pos, ct, exo):
    """`_embed`'s answer with each unit of small rings, merged by shared atom, in its own lattice drawing, or None when
    no unit of `pos` is bent or one finds no free placement: a local honeycomb turned by the lattice symmetry and
    shifted by the step nearest its walked points, the ears of the rings over 8 then walked between the units."""
    groups = _groups([r for r in rings if len(r) <= 8], rank)
    units = []
    bent = False
    for u, rs in sorted(groups, key=lambda g: sorted(rank[x] for x in g[0])):
        r = embed_block(mol, u, rs, rank)
        if r is None or r[1]:
            return None
        fits = _fits(r[0], pos, sorted(u, key=rank.__getitem__))
        bent |= fits[0][0] >= SNAP * len(u)
        units.append((u, rs, fits))
    if not bent:
        return None
    board = Board()
    placed = set()
    small = []
    for u, rs, fits in units:
        bonds = [(b.n, b.m) for b in mol.bonds() if b.n in u and b.m in u]
        new = next((t for _, t in fits if not any(board.bad_atom(p) for p in t.values()) and
                    not any(board.bad_bond(t[a], t[b]) for a, b in bonds)), None)
        if new is None:
            return None
        for x, p in new.items():
            board.put(x, p)
            placed.add(x)
        for a, b in bonds:
            board.link(a, b)
        small += [[project(new[x]) for x in r] for r in rs]
    # a chain atom more than KEEP bonds from every unit keeps its walked point; the ears rewalk the rest
    far = {x: 0 for x in placed}
    front = list(placed)
    while front:
        nxt = []
        for x in front:
            for y in mol.neighbors_of(x):
                if y in atoms and y not in far:
                    far[y] = far[x] + 1
                    nxt.append(y)
        front = nxt
    for x in sorted(atoms, key=rank.__getitem__):
        if x not in placed and far.get(x, 0) > KEEP and not board.bad_atom(pos[x]):
            board.put(x, pos[x])
            placed.add(x)
    for b in mol.bonds():
        n, m = b.n, b.m
        if n in placed and m in placed and not _bonded(board, n, m) and not board.bad_bond(pos[n], pos[m]):
            board.link(n, m)
    todo = [r for r in rings if len(r) > 8]
    conflicts, total = _ears(mol, board, placed, todo, rank, None, small, [0, 2 * EXPANSIONS * len(todo)])
    return _finish(mol, atoms, board, placed, conflicts, total, ct, exo)


def _fits(ideal, pos, xs):
    """`[(squared residual, ideal moved onto pos)]` over `xs`, cheapest first: each lattice symmetry, shifted by the
    integer step nearest the centroids' offset or one lattice step from it."""
    out = []
    n = len(xs)
    for j, g in enumerate(SYMMETRIES):
        moved = {x: apply(ideal[x], g) for x in xs}
        dq = round(fsum(pos[x][0] - moved[x][0] for x in xs) / n)
        dr = round(fsum(pos[x][1] - moved[x][1] for x in xs) / n)
        for k, d in enumerate(((0, 0),) + DIRS):
            t = {x: (q + dq + d[0], r + dr + d[1]) for x, (q, r) in moved.items()}
            err = fsum(hypot(*(sub(project(t[x]), project(pos[x])))) ** 2 for x in xs)
            out.append((err, j, k, t))
    out.sort(key=lambda e: e[:3])
    return [(e[0], e[3]) for e in out]


def planar_block(mol, atoms, rings, rank):
    """`{atom: axial point}` off the lattice with no crossing, or None: the barycentric (Tutte) drawing inside a
    peripheral cycle, an induced cycle whose removal leaves the block connected, as a regular polygon. Of the
    candidate cycles the one holding most atoms bonded out of the block wins, then the one whose drawing keeps its
    closest pair widest against the mean bond. The winner is redrawn with each bond's weight scaled by its length
    over the median inner bond, which any positive weights keep crossing-free. Its median bond is `L` unless two
    atoms come closer than `TUTTE_SPREAD` L."""
    order = sorted(atoms, key=rank.__getitem__)
    nb = {x: sorted((y for y in mol.neighbors_of(x) if y in atoms), key=rank.__getitem__) for x in order}
    bonds = [(x, y) for x in order for y in nb[x] if rank[x] < rank[y]]
    exo = {x for x in atoms if any(y not in atoms for y in mol.neighbors_of(x))}
    seen = set()
    good = []
    cands = sorted(list(rings) + _long_cycles(mol, atoms, rank, keep=12),
                   key=lambda c: (-len(c), sorted(rank[x] for x in c)))
    for c in cands:
        k = frozenset(c)
        if k in seen:
            continue
        seen.add(k)
        if not _peripheral(nb, atoms, c):
            continue
        pts = _tutte(nb, atoms, c, rank)
        if _clean(pts, atoms, bonds):
            good.append(((len(exo.intersection(c)), round(_spread(pts, bonds), 9)), len(good), pts, c))
    best = None
    for score, _, pts, c in sorted(good, key=lambda g: (g[0], -g[1]), reverse=True)[:TUTTE_REDRAWN]:
        pts = _reweigh(nb, atoms, bonds, c, pts, rank)
        if _clean(pts, atoms, bonds):
            score = score[0], round(_spread(pts, bonds), 9)
            if best is None or score > best[0]:
                best = score, pts
    if best is None:
        return None
    return _to_axial(_even(nb, atoms, bonds, best[1], rank), bonds, TUTTE_SPREAD)


def bridge_block(mol, atoms, rings, rank):
    """`{atom: axial point}` of a block bridged by one atom across the para atoms of a hexagon, or None: the block
    less that atom embedded clean on the lattice, the atom inside its hexagon, `LIFT` bonds off the centre towards
    the half whose atoms bond fewer atoms outside the hexagon. The bridge atom bonds nothing but its two bridgeheads;
    of several, the lowest-ranked whose remainder embeds clean is drawn."""
    for p in sorted(atoms, key=rank.__getitem__):
        if len(mol.neighbors_of(p)) != 2:
            continue
        a, b = mol.neighbors_of(p)
        rest = atoms - {p}
        through = [r for r in rings if p in r]
        if len(through) > 2:
            continue
        edges = set()
        for r in through:
            edges ^= {frozenset(e) for e in zip(r, r[1:] + r[:1])}
        edges -= {frozenset((p, a)), frozenset((p, b))}
        six = _cycle_of(edges) if len(through) == 2 else next((r for r in rings if p not in r and a in r and b in r),
                                                             None)
        if six is None or len(six) != 6 or {a, b} - set(six) or (six.index(a) - six.index(b)) % 6 != 3:
            continue
        rs = [r for r in rings if p not in r]
        if len(through) == 2:
            rs.append(_rotated(six, rank))
        r = embed_block(mol, rest, rs, rank)
        if r is None or r[1]:
            continue
        pts = {x: project(q) for x, q in r[0].items()}
        cx, cy = fsum(pts[x][0] for x in six) / 6, fsum(pts[x][1] for x in six) / 6
        if any(abs(hypot(pts[x][0] - cx, pts[x][1] - cy) - L) > 1e-6 for x in six):
            continue
        (ax, ay), (bx, by) = pts[a], pts[b]
        nx, ny = (ay - by) / 2 / L, (bx - ax) / 2 / L       # unit normal of a -> b

        def crowd(sign):
            half = [x for x in six if x not in (a, b) and sign * ((pts[x][0] - cx) * nx + (pts[x][1] - cy) * ny) > 0]
            out = sum(y not in six and y != p for x in half for y in mol.neighbors_of(x))
            return out, min(map(rank.__getitem__, half))
        k = 1 if crowd(1) < crowd(-1) else -1
        px, py = cx + k * LIFT * L * nx, cy + k * LIFT * L * ny
        q = px / L / C30
        return {**r[0], p: (q, py / L - q / 2)}
    return None


def _cycle_of(edges):
    """The atoms of `edges` in cyclic order, or None when they are not one simple cycle."""
    nb = {}
    for e in edges:
        x, y = e
        nb.setdefault(x, []).append(y)
        nb.setdefault(y, []).append(x)
    if not nb or any(len(v) != 2 for v in nb.values()):
        return None
    start = min(nb)
    ring, prev = [start], None
    while True:
        x = ring[-1]
        y = nb[x][0] if nb[x][0] != prev else nb[x][1]
        if y == start:
            break
        prev = x
        ring.append(y)
    return tuple(ring) if len(ring) == len(nb) else None


def _spectral(order, adj, shift):
    """Three coordinates per atom: the eigenvectors of `shift + A` after the constant one, by subspace iteration."""
    n = len(order)
    one = [n ** -.5] * n

    def orth(vs):
        out = []
        for v in vs:
            for u in [one] + out:
                d = fsum(a * b for a, b in zip(u, v))
                v = [a - d * b for a, b in zip(v, u)]
            k = fsum(a * a for a in v) ** -.5
            out.append([a * k for a in v])
        return out
    vs = orth([[(i * 7919 + k * 104729) % 1009 / 1009 - .5 for i in range(n)] for k in range(3)])
    for _ in range(SPHERE_STEPS):
        vs = orth([[shift * v[i] + fsum(v[j] for j in adj[i]) for i in range(n)] for v in vs])
    return [[vs[0][i], vs[1][i], vs[2][i]] for i in range(n)]


def _seen(p, x, y, z):
    """`p` seen in perspective from `EYE` radii down `z`, `x` right and `y` up."""
    r = max(fsum(a * a for a in q) for q in p) ** .5
    out = []
    for q in p:
        u, v, w = (fsum(c * e for c, e in zip(q, axis)) for axis in (x, y, z))
        f = EYE * r / (EYE * r - w)
        out.append((u * f, v * f))
    return out


def _cross3(a, b):
    return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]


def sphere_block(mol, atoms, rings, rank):
    """`{atom: axial point}` of a closed cage, or None: a block of at least `SPHERE` atoms, each with three bonds in
    it, a fullerene among them. Its topological coordinates, the three Laplacian eigenvectors after the constant one,
    are seen in perspective from `EYE` radii down the axis of its lowest-ranked smallest ring, that ring's
    lowest-ranked atom up. Front and back overlap, as in a model held to the eye."""
    order = sorted(atoms, key=rank.__getitem__)
    if len(order) < SPHERE:
        return None
    idx = {x: i for i, x in enumerate(order)}
    adj = [[idx[y] for y in mol.neighbors_of(x) if y in atoms] for x in order]
    if any(len(a) != 3 for a in adj):
        return None
    # 3 + A: its top eigenvector is the constant one and the next three span the coordinates
    p = _spectral(order, adj, 3)
    ring = min(rings, key=lambda r: (len(r), sorted(rank[x] for x in r)))
    z = _unit([fsum(p[idx[x]][k] for x in ring) for k in range(3)])
    top = p[idx[min(ring, key=rank.__getitem__)]]
    d = fsum(a * b for a, b in zip(top, z))
    y = _unit([a - d * b for a, b in zip(top, z)])
    pts = dict(zip(order, _seen(p, _cross3(y, z), y, z)))
    return _to_axial(pts, [(a, b) for a in order for b in mol.neighbors_of(a) if b in atoms and idx[a] < idx[b]], 0)


def cage_block(mol, atoms, rings, rank, dense=DENSE):
    """`{atom: axial point}` of a dense cage, or None: a block of at least `CAGE` atoms and `dense` rings per atom.
    Spectral coordinates are relaxed by stress majorization toward graph distances, 1 for a bond and 1.633 for a
    tetrahedral 1,3 pair, then seen in perspective from `EYE` radii down whichever of `VIEWS` directions has the
    fewest crossings, crowded pairs, atoms on bonds and uneven bonds. Edges cross, as in a model held to the eye."""
    order = sorted(atoms, key=rank.__getitem__)
    n = len(order)
    if n < CAGE or len(rings) < dense * n:
        return None
    idx = {x: i for i, x in enumerate(order)}
    adj = [[idx[y] for y in mol.neighbors_of(x) if y in atoms] for x in order]
    hop = []
    for s in range(n):
        d = [-1] * n
        d[s] = 0
        queue = [s]
        for u in queue:
            for v in adj[u]:
                if d[v] < 0:
                    d[v] = d[u] + 1
                    queue.append(v)
        hop.append([1.633 + (k - 2) * .9 if k > 1 else 1. for k in d])
    p = [[3 * a for a in q] for q in _spectral(order, adj, max(map(len, adj)))]
    for _ in range(CAGE_STEPS):                             # SMACOF, weights the inverse squared target
        q = []
        for i in range(n):
            acc, ws = [0., 0., 0.], 0.
            for j in range(n):
                if i != j:
                    t = hop[i][j]
                    w = t ** -2
                    dx = [p[i][k] - p[j][k] for k in range(3)]
                    d = fsum(a * a for a in dx) ** .5 or 1e-9
                    for k in range(3):
                        acc[k] += w * (p[j][k] + t * dx[k] / d)
                    ws += w
            q.append([a / ws for a in acc])
        p = q
    c = [fsum(q[k] for q in p) / n for k in range(3)]
    p = [[q[k] - c[k] for k in range(3)] for q in p]
    bonds = [(i, j) for i in range(n) for j in adj[i] if i < j]
    best = None
    for s in range(VIEWS):
        z = 1 - 2 * (s + .5) / VIEWS
        r = (1 - z * z) ** .5
        a = s * pi * (3 - 5 ** .5)
        z = [r * cos(a), r * sin(a), z]
        x = _unit(_cross3([1, 0, 0] if abs(z[0]) < .9 else [0, 1, 0], z))
        seen = _seen(p, x, _cross3(z, x), z)
        crossed = sum(1 for k, (a, b) in enumerate(bonds) for c, d in bonds[:k]
                      if len({a, b, c, d}) == 4 and _crosses(seen[a], seen[b], seen[c], seen[d]))
        ls = [hypot(seen[a][0] - seen[b][0], seen[a][1] - seen[b][1]) for a, b in bonds]
        m = sorted(ls)[len(ls) // 2]
        close = sum(1 for i in range(n) for j in range(i)
                    if hypot(seen[i][0] - seen[j][0], seen[i][1] - seen[j][1]) < .5 * m)
        close += sum(1 for i in range(n) for a, b in bonds
                     if i != a and i != b and _seg_dist(seen[i], seen[a], seen[b]) < .25 * m)
        score = crossed + 3 * close + 5 * (max(ls) / min(ls) - 1)
        if best is None or score < best[0]:
            best = score, seen
    return _to_axial(dict(zip(order, best[1])), [(order[a], order[b]) for a, b in bonds], 0)


def _groups(small, rank):
    """`[(atoms, rings)]`: the units of a macrocycle block, its `small` rings sharing an atom, merged."""
    groups = []
    for r in sorted(small, key=lambda r: sorted(rank[x] for x in r)):
        hit = [g for g in groups if g[0] & set(r)]
        groups = [g for g in groups if g not in hit] + [(set(r).union(*(g[0] for g in hit)),
                                                         [r] + [q for g in hit for q in g[1]])]
    return groups


def _unit(v):
    k = fsum(a * a for a in v) ** -.5
    return [a * k for a in v]


def macro_block(mol, atoms, rings, rank):
    """`{atom: axial point}` off the lattice, or None: a block of one ring over 8 atoms whose smaller rings form
    units, each met by the macrocycle along one run of two or more atoms. Each unit keeps its own lattice drawing,
    rigid, its bulk outside; the macrocycle is the closed polygon of its linker atoms and each unit's chord, every
    corner turned as near its ideal as closure allows: 120 degrees at a linker, a radial exit at a unit's end. A
    saturated linker between two units holds its 120 degrees, the units' exits giving way, unless that crowds a
    unit's substituent exits. The
    linkers are at most two per unit: a macrocycle mostly of chain stays on the lattice."""
    big = [r for r in rings if len(r) > 8]
    small = [r for r in rings if len(r) <= 8]
    if len(big) != 1 or not small:
        return None
    groups = _groups(small, rank)
    ring = list(big[0])
    if set(ring).union(*(g[0] for g in groups)) != set(atoms):
        return None
    owner = {x: k for k, (u, _) in enumerate(groups) for x in u}
    free = [i for i, x in enumerate(ring) if x not in owner]
    if not free or len(free) > 2 * len(groups):
        return None
    i = min(free, key=lambda i: rank[ring[i]])              # from the lowest-ranked linker atom
    ring = ring[i:] + ring[:i]
    if rank[ring[-1]] < rank[ring[1]]:
        ring = ring[:1] + ring[1:][::-1]
    n = len(ring)
    verts, sides, want, units, bridges = [], [], [], [], []
    i = 0
    while i < n:
        k = owner.get(ring[i])
        if k is None:                                       # a linker atom
            x = ring[i]
            verts.append(x)
            sides.append(L)
            if owner.get(ring[i - 1]) is not None and owner.get(ring[(i + 1) % n]) is not None and all(
                    mol.bond(x, y).order == 1 for y in mol.neighbors_of(x)):
                bridges.append(len(want))                   # a saturated bridge of two units
            want.append(60.)
            i += 1
            continue
        j = i
        while j + 1 < n and owner.get(ring[j + 1]) == k:
            j += 1
        if j == i or any(g == k for g, _ in units):
            return None                                     # a spiro contact, or a unit met twice
        e, x = ring[i], ring[j]
        u, rs = groups[k]
        if len(rs) == 1:                                    # a lone ring: its regular polygon
            q = rs[0]
            rad = L / 2 / sin(pi / len(q))
            pts = {a: (rad * cos(2 * pi * t / len(q)), rad * sin(2 * pi * t / len(q))) for t, a in enumerate(q)}
        else:
            r = embed_block(mol, u, rs, rank)
            if r is None or r[1]:
                return None
            pts = {a: project(p) for a, p in r[0].items()}
        pts = _chord_frame(pts, e, x, [a for a in u if a not in ring[i:j + 1]])
        if pts is None:
            return None
        units.append((k, pts))
        c = pts[x][0]
        # a radial exit points away from the atom's unit bonds; the walk turns onto the chord and off it
        te, tx = (_radial(mol, pts, a, u) for a in (e, x))
        verts += [e, x]
        sides += [c, L]
        want += [(-(te + 180) + 180) % 360 - 180, (tx + 180) % 360 - 180]
        i = j + 1
    bonds = [(b.n, b.m) for b in mol.bonds() if b.n in atoms and b.m in atoms]
    if bridges:                                             # stiff bridges, while the unit exits stay clear
        out = _macro_draw(verts, sides, want, [BRIDGE if i in bridges else 1. for i in range(len(want))], units)
        if out is not None and _exits_clear(mol, out, atoms, bonds):
            return _macro_done(out, atoms, bonds)
    return _macro_done(_macro_draw(verts, sides, want, None, units), atoms, bonds)


def _macro_draw(verts, sides, want, give, units):
    """The macrocycle polygon closed nearest `want`, each unit set on its chord; None when it does not close."""
    turns = _close(sides, want, give)
    if turns is None or any(not -150 < t < 150 for t in turns):
        return None
    out = {}
    phi, px, py = 0., 0., 0.
    for v, side, t in zip(verts, sides, [0.] + turns[1:]):
        phi += t
        out[v] = (px, py)
        px, py = px + side * cos(phi * pi / 180), py + side * sin(phi * pi / 180)
    frame = {}
    phi = 0.
    for v, t in zip(verts, [0.] + turns[1:]):
        phi += t
        frame[v] = phi * pi / 180
    for _, pts in units:
        e = next(a for a in pts if a in frame and pts[a] == (0., 0.))
        (ox, oy), co, si = out[e], cos(frame[e]), sin(frame[e])
        for a, (qx, qy) in pts.items():
            out[a] = (ox + qx * co - qy * si, oy + qx * si + qy * co)
    return out


def _exits_clear(mol, out, atoms, bonds):
    """Whether each exit of a block atom of two block bonds, a bond length along their outward bisector, keeps
    `EXIT` bond lengths from every other exit and block atom."""
    ends = []
    for x in atoms:
        nb = [y for y in mol.neighbors_of(x) if y in atoms]
        if len(nb) != 2 or len(mol.neighbors_of(x)) == 2:
            continue
        (px, py), (ax, ay), (bx, by) = out[x], out[nb[0]], out[nb[1]]
        wx, wy = 2 * px - ax - bx, 2 * py - ay - by
        d = hypot(wx, wy)
        if d < 1e-9:
            return False
        ends.append((x, px + wx / d * L, py + wy / d * L))
    near = EXIT * L
    return all(hypot(qx - rx, qy - ry) >= near for x, qx, qy in ends for y, (rx, ry) in out.items() if y != x) and \
        all(hypot(a[1] - b[1], a[2] - b[2]) >= near for i, a in enumerate(ends) for b in ends[i + 1:])


def _macro_done(out, atoms, bonds):
    """`out` in axial units, or None when absent, unclean or crowded."""
    xs = sorted(atoms)
    if out is None or len(out) != len(atoms) or not _clean(out, atoms, bonds) or min(
            hypot(out[a][0] - out[b][0], out[a][1] - out[b][1]) for i, a in enumerate(xs) for b in xs[i + 1:]) < .6 * L:
        return None
    return _to_axial(out, bonds)


def _chord_frame(pts, e, x, bulk):
    """`pts` moved so `e` is the origin and `x` on +x, the `bulk` atoms below; None with no bulk."""
    if not bulk:
        return None
    (ex, ey), (xx, xy) = pts[e], pts[x]
    c = hypot(xx - ex, xy - ey)
    ca, sa = (xx - ex) / c, (xy - ey) / c
    out = {a: ((px - ex) * ca + (py - ey) * sa, (py - ey) * ca - (px - ex) * sa) for a, (px, py) in pts.items()}
    out[e] = (0., 0.)
    if sum(out[a][1] for a in bulk) > 0:
        out = {a: (px, -py) for a, (px, py) in out.items()}
    return out


def _radial(mol, pts, a, unit):
    """The direction, degrees, pointing from `a` away from its bonds inside `unit`."""
    x, y = pts[a]
    dx = dy = 0.
    for b in mol.neighbors_of(a):
        if b in unit:
            vx, vy = pts[b][0] - x, pts[b][1] - y
            d = hypot(vx, vy)
            dx -= vx / d
            dy -= vy / d
    return degrees(atan2(dy, dx))


def _close(sides, want, give=None, steps=50):
    """Turns, degrees, nearest `want` in least squares that close the polygon of `sides` counterclockwise; the
    turn at the first vertex is implied. A corner's `give` is its inverse stiffness, 1 by default. None when no
    closure is found."""
    m = len(sides)
    give = give or [1.] * m
    tau = [w * pi / 180 for w in want]
    goal = list(tau)
    for _ in range(steps):
        phi = [0.] * m
        for i in range(1, m):
            phi[i] = phi[i - 1] + tau[i]
        g = [sum(tau) - 2 * pi, sum(s * cos(f) for s, f in zip(sides, phi)),
             sum(s * sin(f) for s, f in zip(sides, phi))]
        jac = [[1.] * m, [0.] * m, [0.] * m]
        tx = ty = 0.
        for i in range(m - 1, 0, -1):
            tx += -sides[i] * sin(phi[i])
            ty += sides[i] * cos(phi[i])
            jac[1][i], jac[2][i] = tx, ty
        d = [gl - t for gl, t in zip(goal, tau)]
        rhs = [gi + sum(jr[i] * d[i] for i in range(m)) for gi, jr in zip(g, jac)]
        a = [[sum(give[i] * ja[i] * jb[i] for i in range(m)) for jb in jac] for ja in jac]
        lam = _solve3(a, rhs)
        if lam is None:
            return None
        step = [d[i] - give[i] * sum(jac[r][i] * lam[r] for r in range(3)) for i in range(m)]
        tau = [t + s for t, s in zip(tau, step)]
        if max(map(abs, step)) < 1e-12:
            break
    phi = 0.
    x = y = 0.
    for i in range(m):
        if i:
            phi += tau[i]
        x += sides[i] * cos(phi)
        y += sides[i] * sin(phi)
    if hypot(x, y) > 1e-6 or abs(sum(tau) - 2 * pi) > 1e-6:
        return None
    return [t * 180 / pi for t in tau]


def _solve3(a, b):
    """`a x = b` for a 3 by 3 `a`, Cramer's rule; None when singular."""
    def det(m):
        return (m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
                m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]))
    d = det(a)
    if abs(d) < 1e-12:
        return None
    return [det([[b[r] if c == k else a[r][c] for c in range(3)] for r in range(3)]) / d for k in range(3)]


def _even(nb, atoms, bonds, pts, rank, steps=200):
    """`pts` moved toward one bond length and apart where atoms crowd; a step that leaves `_clean`, or brings
    the closest pair nearer against the mean bond, is halved."""
    order = sorted(atoms, key=rank.__getitem__)
    size = sorted(hypot(pts[a][0] - pts[b][0], pts[a][1] - pts[b][1]) for a, b in bonds)[len(bonds) // 2]
    pts = _bend(nb, order, bonds, pts, size)
    rate = .1
    for _ in range(steps):
        force = {x: [0., 0.] for x in order}
        for a, b in bonds:
            dx, dy = pts[b][0] - pts[a][0], pts[b][1] - pts[a][1]
            d = hypot(dx, dy) or 1e-9
            f = (d - size) / d
            force[a][0] += f * dx
            force[a][1] += f * dy
            force[b][0] -= f * dx
            force[b][1] -= f * dy
        for i, a in enumerate(order):
            for b in order[i + 1:]:
                if b in nb[a]:
                    continue
                dx, dy = pts[b][0] - pts[a][0], pts[b][1] - pts[a][1]
                d = hypot(dx, dy) or 1e-9
                if d < 1.5 * size:
                    f = (1.5 * size - d) / d * .5
                    force[a][0] -= f * dx
                    force[a][1] -= f * dy
                    force[b][0] += f * dx
                    force[b][1] += f * dy
        new = {x: (pts[x][0] + rate * force[x][0], pts[x][1] + rate * force[x][1]) for x in order}
        if _clean(new, atoms, bonds) and _spread(new, bonds) >= _spread(pts, bonds):
            pts = new
            rate = min(.2, rate * 1.5)
        else:
            rate /= 2
            if rate < 1e-4:
                break
    return pts


def _bend(nb, order, bonds, pts, size):
    """`pts` with each atom of two neighbours on their line moved off it by the offset of `BEND` that keeps the drawing
    clean and leaves it farthest from any other bond: a straight run of the barycentric drawing reads as one bond.
    The move turns counterclockwise about the outer polygon's centre, the origin, so a cage's bridges twist one way;
    an atom near the centre moves away from its nearest other atom. Each line is read before any move."""
    src, pts = pts, dict(pts)
    for x in order:
        if len(nb[x]) != 2:
            continue
        (ax, ay), (bx, by), (px, py) = src[nb[x][0]], src[nb[x][1]], src[x]
        dx, dy = bx - ax, by - ay
        d = hypot(dx, dy) or 1e-9
        if abs(dx * (py - ay) - dy * (px - ax)) / d > .05 * size:
            continue
        near = min((y for y in order if y != x and y not in nb[x]),
                   key=lambda y: hypot(pts[y][0] - px, pts[y][1] - py), default=None)
        if near is None:
            continue
        nx, ny = -dy / d, dx / d
        turn = px * ny - py * nx                            # counterclockwise about the outer polygon's centre
        if abs(turn) > .1 * size:
            if turn < 0:
                nx, ny = -nx, -ny
        elif nx * (pts[near][0] - px) + ny * (pts[near][1] - py) > 0:
            nx, ny = -nx, -ny
        best = None
        for t in BEND:
            new = {**pts, x: (px + nx * t * size, py + ny * t * size)}
            gap = min((_seg_dist(new[x], new[a], new[b]) for a, b in bonds if x != a and x != b), default=size)
            if gap >= CLEAR * size and _clean(new, pts.keys(), bonds) and (best is None or gap > best[0]):
                best = gap, new
        if best is not None:
            pts = best[1]
    return pts


def _spread(pts, bonds):
    """The closest pair of atoms over the mean bond: the figure `_even` raises."""
    xs = list(pts)
    gap = min(hypot(pts[a][0] - pts[b][0], pts[a][1] - pts[b][1]) for i, a in enumerate(xs) for b in xs[i + 1:])
    return gap * len(bonds) / sum(hypot(pts[a][0] - pts[b][0], pts[a][1] - pts[b][1]) for a, b in bonds)


def _peripheral(nb, atoms, cycle):
    ring = set(cycle)
    if sum(y in ring for x in cycle for y in nb[x]) != 2 * len(cycle):   # a chord
        return False
    rest = [x for x in atoms if x not in ring]
    if not rest:
        return True
    reach = {rest[0]}
    stack = [rest[0]]
    while stack:
        for y in nb[stack.pop()]:
            if y not in ring and y not in reach:
                reach.add(y)
                stack.append(y)
    return len(reach) == len(rest)


def _tutte(nb, atoms, cycle, rank):
    n = len(cycle)
    i = min(range(n), key=lambda i: rank[cycle[i]])
    seq = tuple(cycle[i:]) + tuple(cycle[:i])
    if rank[seq[-1]] < rank[seq[1]]:
        seq = seq[:1] + seq[1:][::-1]
    pts = {x: (cos(2 * pi * j / n + pi / 2), sin(2 * pi * j / n + pi / 2)) for j, x in enumerate(seq)}
    for x in atoms:
        pts.setdefault(x, (0., 0.))
    return _relax_inner(nb, atoms, set(seq), pts, rank, None)


def _relax_inner(nb, atoms, outer, pts, rank, weight, tol=1e-12, steps=2000):
    """Gauss-Seidel sweeps moving each atom off `outer` to the `weight`-ed mean of its neighbours."""
    pts = dict(pts)
    inner = sorted((x for x in atoms if x not in outer), key=rank.__getitem__)
    for _ in range(steps):
        worst = 0.
        for x in inner:
            if weight is None:
                w = [1.] * len(nb[x])
            else:
                w = [weight[x, y] for y in nb[x]]
            t = sum(w)
            p = (sum(k * pts[y][0] for k, y in zip(w, nb[x])) / t, sum(k * pts[y][1] for k, y in zip(w, nb[x])) / t)
            worst = max(worst, abs(p[0] - pts[x][0]) + abs(p[1] - pts[x][1]))
            pts[x] = p
        if worst < tol:
            break
    return pts


def _reweigh(nb, atoms, bonds, cycle, pts, rank, rounds=30):
    """The barycentric drawing `pts` inside `cycle` redrawn `rounds` times, each bond's weight multiplied by the
    square root of its length over the median inner bond: a nested ring grows toward one bond length."""
    outer = set(cycle)
    weight = {}
    for a, b in bonds:
        weight[a, b] = weight[b, a] = 1.
    for _ in range(rounds):
        ls = {(a, b): hypot(pts[a][0] - pts[b][0], pts[a][1] - pts[b][1]) for a, b in bonds}
        inner = sorted(v for (a, b), v in ls.items() if a not in outer or b not in outer)
        if not inner:
            break
        size = inner[len(inner) // 2]
        for (a, b), v in ls.items():
            weight[a, b] = weight[b, a] = weight[a, b] * (v / size) ** .5
        pts = _relax_inner(nb, atoms, outer, pts, rank, weight, 1e-10, 3000)
    return pts


def _clean(pts, atoms, bonds):
    xs = sorted(atoms)
    for i, a in enumerate(xs):
        for b in xs[i + 1:]:
            if hypot(pts[a][0] - pts[b][0], pts[a][1] - pts[b][1]) < 1e-3:
                return False
    size = min(hypot(pts[a][0] - pts[b][0], pts[a][1] - pts[b][1]) for a, b in bonds)
    for j, (a, b) in enumerate(bonds):
        for c, d in bonds[j + 1:]:
            if len({a, b, c, d}) == 4 and _fcross(pts[a], pts[b], pts[c], pts[d]):
                return False
        for x in atoms:
            if x != a and x != b and _seg_dist(pts[x], pts[a], pts[b]) < .25 * size:
                return False
    return True


def _fcross(p, q, a, b):
    def side(o, u, v):
        return (u[0] - o[0]) * (v[1] - o[1]) - (u[1] - o[1]) * (v[0] - o[0])
    return side(p, q, a) * side(p, q, b) < -1e-12 and side(a, b, p) * side(a, b, q) < -1e-12


def _to_axial(pts, bonds, spread=.8):
    """Scaled so the median bond is `L`, and no two atoms closer than `spread` L; then in axial units."""
    lengths = sorted(hypot(pts[a][0] - pts[b][0], pts[a][1] - pts[b][1]) for a, b in bonds)
    xs = list(pts)
    k = 1 / lengths[len(lengths) // 2]                      # in units of L
    if spread:
        gap = min(hypot(pts[a][0] - pts[b][0], pts[a][1] - pts[b][1]) for i, a in enumerate(xs) for b in xs[i + 1:])
        k = max(k, spread / gap)
    out = {}
    for x, (px, py) in pts.items():
        q = px * k / C30
        out[x] = (q, py * k - q / 2)
    return out


def _long_cycles(mol, atoms, rank, keep=3, budget=20000):
    """The `keep` longest simple cycles of the block found within `budget` steps, rank-ordered."""
    out = sorted(_cycles(mol, atoms, rank, len(atoms), budget),
                 key=lambda c: (-len(c), sorted(rank[x] for x in c)))
    return out[:keep]


def _short_cycles(mol, atoms, rank, top=8, budget=20000):
    """The simple cycles of the block of at most `top` atoms found within `budget` steps."""
    return _cycles(mol, atoms, rank, top, budget)


def _cycles(mol, atoms, rank, top, budget, whole=False):
    """Simple cycles of at most `top` atoms, each once, from its lowest-ranked atom, found within `budget` steps;
    None when `whole` asks for all of them and the budget runs out."""
    order = sorted(atoms, key=rank.__getitem__)
    nb = {x: sorted((y for y in mol.neighbors_of(x) if y in atoms), key=rank.__getitem__) for x in atoms}
    found = {}
    steps = 0
    for s in order:
        low = rank[s]
        stack = [(s, [s])]
        while stack and steps < budget:
            x, path = stack.pop()
            steps += 1
            for y in nb[x]:
                if y == s and len(path) > 2:
                    c = frozenset(path)
                    if c not in found or len(path) > len(found[c]):
                        found[c] = tuple(path)
                elif rank[y] > low and y not in path and len(path) < top:
                    stack.append((y, path + [y]))
        if stack and whole:
            return None
    return list(found.values())


def _embed(mol, atoms, rings, rank, first, spent, ct, exo):
    """`(points, conflicts, broken, crowded, cost)` of the ear walk from `first`, or None; `spent` is `_route`'s
    budget, `broken` the frames of `ct` the points break, `crowded` the bonds out of the block, `exo` counting them
    per atom, with no free unit step to take."""
    board = Board()
    placed = set()
    shape = _shape(len(first))
    # the ring starts from its lowest-ranked atom, in the direction of its lower-ranked neighbour, then turns on
    # its shape to the first rotation, mirrored second, that breaks fewest of the cis/trans bonds it holds, then
    # lays fewest long edges of its shape on a bond fused to another ring: a pentagon's gap, then fewest fused bonds
    # beside a 240-degree corner, which leaves the fused ring no room outside; alone in its block,
    # then leaves fewest exits without an outward step, fewest two-exit atoms beside a 240-degree corner, fewest
    # exits on one, and, over 8 atoms, fewest exits off a 120-degree corner: an odd ring's bulge
    first = tuple(first)
    i = min(range(len(first)), key=lambda i: rank[first[i]])
    seq = first[i:] + first[:i]
    if rank[seq[-1]] < rank[seq[1]]:
        seq = seq[:1] + seq[1:][::-1]
    own = [f for f in ct if all(x in seq for x in f[:4])]
    alone = len(rings) == 1 and any(x in exo for x in seq)
    n = len(seq)
    gaps = [j for j in range(n) if sub(shape[(j + 1) % n], shape[j]) not in DIRS]
    fused = {frozenset(e) for r in rings if set(r) != set(seq) for e in zip(r, r[1:] + r[:1])}
    if own or alone or gaps and fused:
        turns = [seq[k:] + seq[:k] for k in range(n)]
        turns += [t[::-1] for t in turns]
        out = _outward(shape)

        def fit(t):
            edges = [frozenset((t[j], t[(j + 1) % n])) in fused for j in range(n)]
            broken = _ct_broken(dict(zip(t, shape)), own), sum(edges[j] for j in gaps), \
                sum(e and min(out[j], out[(j + 1) % n]) < 3 for j, e in enumerate(edges))
            if not alone:
                return broken
            crowd = sum(max(0, exo.get(x, 0) - k) for x, k in zip(t, out))
            pocket = sum(exo.get(x, 0) > 1 and (out[j - 1] < 2 or out[(j + 1) % n] < 2) for j, x in enumerate(t))
            return broken, crowd, pocket, sum(x in exo for x, k in zip(t, out) if k < 2), \
                n > 8 and sum(x in exo for x, k in zip(t, out) if k != 3)
        seq = min(turns, key=fit)
    for x, p in zip(seq, shape):
        board.put(x, p)
        placed.add(x)
    for j in range(len(seq)):
        board.link(seq[j], seq[(j + 1) % len(seq)])
    centre = sublattice((0, 0)) if len(first) >= 5 else None
    small = [[project(board.pos[x]) for x in seq]]           # placed rings' polygons, for the inside cost
    todo = [r for r in rings if set(r) != set(first)]
    conflicts, total = _ears(mol, board, placed, todo, rank, centre, small, spent)
    return _finish(mol, atoms, board, placed, conflicts, total, ct, exo)


def _ears(mol, board, placed, todo, rank, centre, small, spent):
    """`(conflicts, cost)` of walking the ears of the `todo` rings onto `board`, each the ring with most placed
    atoms; `small` collects the placed rings' polygons."""
    conflicts = 0
    total = 0.
    while todo:
        cand = [r for r in todo if sum(x in placed for x in r) >= 2]
        if not cand:
            break
        r = min(cand, key=lambda r: (len(r) > 8, -sum(x in placed for x in r), len(r) != 6,
                                     sorted(rank[x] for x in r)))
        runs = _runs(r, placed)
        if not runs:
            todo.remove(r)
            small.append([project(board.pos[x]) for x in r])
            continue
        a, mid, b = min(runs, key=lambda t: sorted(rank[x] for x in t[1]))
        if rank[b] < rank[a]:                               # walk from the lower-ranked end: order-free
            a, mid, b = b, mid[::-1], a
        polys = small if len(r) <= 8 else ()
        found = _route(board, mol, a, b, len(mid), centre, polys, placed, spent)
        if found is None and len(mid) > EXACT:
            found = _route(board, mol, a, b, len(mid), centre, polys, placed, spent, GREED)
        if found is None:
            pts = _straight(board, a, b, len(mid))
            conflicts += 1
        else:
            pts, cost = found
            total += cost
        prev = a
        for x, p in zip(mid, pts):
            board.put(x, p)
            placed.add(x)
            board.link(prev, x)
            prev = x
        board.link(prev, b)
    return conflicts, total


def _finish(mol, atoms, board, placed, conflicts, total, ct, exo):
    """`_embed`'s answer from a walked `board`, the bonds no ear drew linked and counted when bad; None when an atom
    is unplaced."""
    for x in atoms:
        if x not in placed:                                 # not reached through any ring
            return None
    for bond in mol.bonds():
        n, m = bond.n, bond.m
        if n in atoms and m in atoms and not _bonded(board, n, m):
            conflicts += board.bad_bond(board.pos[n], board.pos[m])
            board.link(n, m)
    return dict(board.pos), conflicts, _ct_broken(board.pos, ct), _crowded(board, exo), total


def _outward(shape):
    """Per point of a cyclic `shape`, the lattice directions leaving it to the outside: 3 at a 120-degree corner,
    1 at a 240-degree one."""
    poly = [project(p) for p in shape]
    n = len(shape)
    out = []
    for i, p in enumerate(shape):
        ring = (sub(shape[i - 1], p), sub(shape[(i + 1) % n], p))
        out.append(sum(d not in ring and not _inside((p[0] + d[0] / 2, p[1] + d[1] / 2), (poly,)) for d in DIRS))
    return out


def _crowded(board, exo):
    """Bonds out of the block, `exo` counting them per atom, beyond the free unit steps of their atom."""
    bad = 0
    for x, need in exo.items():
        p = board.pos[x]
        free = sum(not board.bad_atom(q) and not board.bad_bond(p, q) for q in (add(p, d) for d in DIRS))
        bad += max(0, need - free)
    return bad


def _ring_ct(mol, atoms, rank):
    """`[(a, anchor, partner, b, same)]`: each stated cis/trans bond of the block with a block neighbour on both
    ends, `a` and `b` the lowest-ranked of those, `same` when they share a side."""
    out = []
    for u in mol.stereo_units():
        if u['kind'] != SU_CIS_TRANS or not u['parity'] or u['anchor'] not in atoms:
            continue
        f = cis_trans_frame(mol, u)
        if f is None:
            continue
        near, anchor, partner, far = f
        if partner not in atoms or partner not in mol.neighbors_of(anchor):
            continue
        a = [y for y in mol.neighbors_of(anchor) if y in atoms and y != partner]
        b = [y for y in mol.neighbors_of(partner) if y in atoms and y != anchor]
        if a and b:
            a, b = min(a, key=rank.__getitem__), min(b, key=rank.__getitem__)
            out.append((a, anchor, partner, b, (u['parity'] == 2) ^ (a != near) ^ (b != far)))
    return out


def _ct_broken(pos, frames):
    """Frames whose two substituents in `pos` lie on the wrong sides of their bond, or on its line."""
    bad = 0
    for a, anchor, partner, b, same in frames:
        if a not in pos or b not in pos:
            continue
        pa, pb = project(pos[anchor]), project(pos[partner])
        dx, dy = pb[0] - pa[0], pb[1] - pa[1]
        sa, sb = (dx * (q[1] - pa[1]) - dy * (q[0] - pa[0]) for q in (project(pos[a]), project(pos[b])))
        if abs(sa) < 1e-6 or abs(sb) < 1e-6 or (sa * sb > 0) != same:
            bad += 1
    return bad


def _runs(ring, placed):
    """`[(a, [new atoms], b)]`: maximal unplaced runs of `ring` between two placed atoms, in ring order."""
    n = len(ring)
    start = next((i for i in range(n) if ring[i] in placed), None)
    if start is None:
        return []
    seq = ring[start:] + ring[:start]
    out = []
    i = 0
    while i < n:
        if seq[i] in placed:
            j = i + 1
            while j < n and seq[j] not in placed:
                j += 1
            if j > i + 1:
                out.append((seq[i], list(seq[i + 1:j]), seq[j % n]))
            i = j
        else:
            i += 1
    return out


def _straight(board, a, b, m):
    """A fallback ear: `m` points on the straight line from `a` to `b`, not necessarily on the lattice."""
    pa, pb = board.pos[a], board.pos[b]
    return [(pa[0] + (pb[0] - pa[0]) * (i + 1) / (m + 1), pa[1] + (pb[1] - pa[1]) * (i + 1) / (m + 1))
            for i in range(m)]
