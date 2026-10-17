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
"""The controller's edits of a plane: a ring system redrawn, or an acyclic side turned, mirrored or stretched.

| kind | move | offered |
| --- | --- | --- |
| `embed` | a ring system drawn as one of its `variants`; each piece not holding the root turns rigidly onto its \
anchor's exit | each undrawn system, in rank order |
| `flip` | the smaller side of an acyclic bond mirrored across the bond | each frontier bond |
| `rotate` | the side turned about the pivot onto a target of `targets` | up to `TARGETS` per frontier bond |
| `stretch` | the bond lengthened or shortened by `STEP`, between `L` and `LONGEST` | up to 2 per frontier bond |
| `pivot` | the smaller half at a spiro atom turned so its ring bonds' bisector is opposite the other half's | each \
frontier spiro atom |
| `even` | the terminal ligands of an acyclic atom of 5 or more neighbours, one held at most, set at `L` on its \
`360 / n` slots | each frontier such atom |
| `bend` | the side turned `BEND_TURN` either way at a pivot of two neighbours | only with `violate` |

| target of a rotation | order |
| --- | --- |
| 120 degrees off each other bond of the pivot | 1 |
| the bisector of each gap between the pivot's other bonds that is not a ring's corner | 2 |
| the widest gap between ring bonds outside every ring, split evenly among the pivot's acyclic bonds | 3 |
| `360 / n` either way from the side's bond, at an acyclic pivot of n of 4 or more | 4 |

Targets within `CLEAR` of another bond of the pivot are dropped, and so are those within `SAME` of one listed
before or of the side's bond; the `TARGETS` nearest the side's bond are kept. A move that breaks a stated
cis/trans or bends a linear atom is `masked`.
"""
from math import atan2, cos, hypot, pi, radians, sin
from typing import NamedTuple

from .lattice import L, project
from .quality import _inside
from .relax import _bent, _ct_broken, _moved, _side
from .score import frontier
from .rings import bridge_block, cage_block, embed_block, macro_block, planar_block, sphere_block
from .templates import template_block


__all__ = ['Move', 'TARGET_TYPES', 'targets', 'targets_typed', 'outside_angles', 'candidates', 'arc_moves', 'apply',
           'masked', 'broken', 'move_key', 'variants', 'embed']

TARGETS = 4
SAME = radians(5)
CLEAR = radians(30)
STEP = L / 2
LONGEST = 3 * L
BEND_TURN = pi / 6
ARCS = 8                # arcs offered per ring per step
ARC_TURN = pi / 6       # an arc's turn about its first end
KINDS = {'embed': 0, 'flip': 1, 'arc_flip': 2, 'pivot': 3, 'even': 4, 'rotate': 5, 'arc_rotate': 6, 'stretch': 7,
         'arc_stretch': 8, 'bend': 9, 'pin': 10, 'turn': 11}


class Move(NamedTuple):
    kind: str               # one of `KINDS`
    pivot: int              # the atom held; an embed's block index
    other: int              # the side's atom on the bond; an embed's variant index
    arg: float              # a turn in radians, a stretch's length; 0 for a flip or an embed
    part: frozenset         # the atoms moved; empty for an embed


def move_key(scope, move):
    """The order of two equally scored moves: by kind, then by the ranks of their atoms."""
    if move.kind == 'embed':
        return 0, move.pivot, move.other, 0.
    return KINDS[move.kind], scope.rank[move.pivot], scope.rank[move.other], round(move.arg, 9)


def _dir(plane, x, y):
    return atan2(plane[y][1] - plane[x][1], plane[y][0] - plane[x][0]) % (2 * pi)


def _off(a, b):
    d = (a - b) % (2 * pi)
    return min(d, 2 * pi - d)


def _gaps(dirs):
    """`[(start, width, i, j)]` counterclockwise between consecutive sorted directions, `i`, `j` their indices."""
    order = sorted(range(len(dirs)), key=dirs.__getitem__)
    out = []
    for k, i in enumerate(order):
        j = order[(k + 1) % len(order)]
        out.append((dirs[i], (dirs[j] - dirs[i]) % (2 * pi) or 2 * pi, i, j))
    return out


def _corners(ring):
    return {(ring[i], frozenset((ring[i - 1], ring[(i + 1) % len(ring)]))) for i in range(len(ring))}


def outside_angles(scope, plane, x, k, ring=None):
    """`k` directions splitting evenly the widest gap at `x` between its bonds to `ring` (default: its ring bonds)
    that is no ring's corner, nor inside the ring of its two bonds; empty when `x` has fewer than two such bonds."""
    nb = [y for y in scope.adj[x] if (y in ring if ring is not None else frozenset((x, y)) in scope.ring_bonds)]
    if len(nb) < 2 or not k:
        return []
    dirs = [_dir(plane, x, y) for y in nb]
    polys = [[plane[z] for z in r] for r in scope.rings if (x, frozenset(nb)) in _corners(r) and
             all(z in plane for z in r)] if len(nb) == 2 else []
    best = None
    for s, w, i, j in _gaps(dirs):
        if len(nb) > 2 and (x, frozenset((nb[i], nb[j]))) in scope.corners:
            continue
        a = s + w / 2
        probe = plane[x][0] + .5 * L * cos(a), plane[x][1] + .5 * L * sin(a)
        if any(_inside(probe, p) for p in polys):           # a reflex corner's widest gap is the ring's inside
            continue
        if best is None or w > best[1] + 1e-9:
            best = s, w
    if best is None:
        return []
    s, w = best
    return [(s + w * t / (k + 1)) % (2 * pi) for t in range(1, k + 1)]


TARGET_TYPES = ('lattice', 'bisector', 'split', 'turn')


def targets_typed(scope, plane, pivot, other):
    """`[(direction, type)]` from `pivot` the side through `other` may turn onto, nearest its bond first; `type`
    names its row of the rotation-target table, an index into `TARGET_TYPES`."""
    cur = _dir(plane, pivot, other)
    rest = [y for y in scope.adj[pivot] if y != other]
    held = [_dir(plane, pivot, y) for y in rest]
    listed = []
    for h in held:
        listed += [((h + 2 * pi / 3) % (2 * pi), 0), ((h - 2 * pi / 3) % (2 * pi), 0)]
    if len(held) > 1:
        for s, w, i, j in _gaps(held):
            if (pivot, frozenset((rest[i], rest[j]))) not in scope.corners:
                listed.append(((s + w / 2) % (2 * pi), 1))
    exits = sum(frozenset((pivot, y)) not in scope.ring_bonds for y in scope.adj[pivot])
    listed += [(a, 2) for a in outside_angles(scope, plane, pivot, exits)]
    n = len(scope.adj[pivot])
    if n >= 4 and not any(frozenset((pivot, y)) in scope.ring_bonds for y in rest):
        listed += [((cur + 2 * pi / n) % (2 * pi), 3), ((cur - 2 * pi / n) % (2 * pi), 3)]
    out = []
    for a, t in listed:
        if all(_off(a, h) >= CLEAR for h in held) and _off(a, cur) >= SAME and all(_off(a, b) >= SAME for b, _ in out):
            out.append((a, t))
    out.sort(key=lambda p: round(_off(p[0], cur), 9))
    return out[:TARGETS]


def targets(scope, plane, pivot, other):
    """Absolute directions from `pivot` the side through `other` may turn onto, nearest its bond first."""
    return [a for a, _ in targets_typed(scope, plane, pivot, other)]


def _sided(scope, n, m):
    """`(pivot, other, part)` of the acyclic bond `n`-`m`: the smaller side moves."""
    if (n, m) not in scope.sides:
        scope.sides[n, m] = frozenset(_side(scope.mol, n, m))
    side = scope.sides[n, m]
    if 2 * len(side) <= len(scope.atoms):
        return n, m, side
    return m, n, frozenset(scope.atoms) - side


def _length(plane, n, m):
    return hypot(plane[n][0] - plane[m][0], plane[n][1] - plane[m][1])


def candidates(scope, plane, front, violate=False):
    """The primitive moves on acyclic bonds with an atom in `front`, in rank order of their bonds, then the `pivot` and
    `even` moves of atoms in `front`, in rank order; `violate` offers only bends."""
    out = []
    for n, m in scope.bonds:
        if (n, m) not in scope.acyclic or n not in front and m not in front:
            continue
        pivot, other, part = _sided(scope, n, m)
        if violate:
            if len(scope.adj[pivot]) == 2 and pivot not in scope.linear:
                out += [Move('bend', pivot, other, BEND_TURN, part), Move('bend', pivot, other, -BEND_TURN, part)]
            continue
        out.append(Move('flip', pivot, other, 0., part))
        cur = _dir(plane, pivot, other)
        for t in targets(scope, plane, pivot, other):
            out.append(Move('rotate', pivot, other, (t - cur + pi) % (2 * pi) - pi, part))
        d = _length(plane, pivot, other)
        for s in (d + STEP, d - STEP):
            if L - 1e-9 <= s <= LONGEST + 1e-9:
                out.append(Move('stretch', pivot, other, s, part))
    if not violate:
        for x in scope.atoms:
            if x in front:
                out += _atom_moves(scope, plane, x)
    return out


def _spiro_halves(scope, x):
    """`[(ring neighbours, half)]` of a spiro atom `x`: each of its two rings' neighbour pair with the atoms reached
    from it without `x`; empty for any other atom."""
    pairs = [ys for z, ys in scope.corners if z == x]
    if len(scope.adj[x]) != 4 or len(pairs) != 2 or pairs[0] & pairs[1]:
        return []
    out = []
    for ys in pairs:
        start = min(ys, key=scope.rank.__getitem__)
        seen, stack = {start}, [start]
        while stack:
            y = stack.pop()
            for z in scope.adj[y]:
                if z != x and z not in seen:
                    seen.add(z)
                    stack.append(z)
        if seen & pairs[0] and seen & pairs[1]:                 # one ring system through both pairs: not spiro
            return []
        out.append((ys, frozenset(seen)))
    return out


def _bisector(plane, x, ys):
    """The direction halving the narrower angle at `x` between its bonds to the two atoms of `ys`."""
    u, v = (_dir(plane, x, y) for y in ys)
    d = (v - u) % (2 * pi)
    return (u + d / 2) % (2 * pi) if d < pi else (v + (2 * pi - d) / 2) % (2 * pi)


def _atom_moves(scope, plane, x):
    """The `pivot` move of a spiro atom's smaller half and the `even` move of an acyclic atom's terminal ligands."""
    out = []
    halves = _spiro_halves(scope, x)
    if halves:
        (ya, a), (yb, b) = sorted(halves, key=lambda h: (len(h[1]), min(scope.rank[z] for z in h[1])))
        turn = (_bisector(plane, x, yb) + pi - _bisector(plane, x, ya) + pi) % (2 * pi) - pi
        if abs(turn) >= SAME:
            out.append(Move('pivot', x, min(ya, key=scope.rank.__getitem__), turn, a))
    nb = scope.adj[x]
    tips = [y for y in nb if len(scope.adj[y]) == 1]
    if len(nb) >= 5 and len(tips) >= 2 and len(nb) - len(tips) <= 1 and \
            not any(frozenset((x, y)) in scope.ring_bonds for y in nb):
        out.append(Move('even', x, x, 0., frozenset(tips)))
    return out


def _even(scope, plane, x, tips):
    """Terminal ligands `tips` of `x` at `L` on the `360 / n` slots from its held bond, or from its lowest-ranked
    ligand when it holds none, each keeping its turn order."""
    nb = scope.adj[x]
    n = len(nb)
    held = [y for y in nb if y not in tips]
    first = held[0] if held else min(tips, key=scope.rank.__getitem__)
    r = _dir(plane, x, first)
    order = sorted((y for y in tips if y != first), key=lambda y: (round((_dir(plane, x, y) - r) % (2 * pi), 9),
                                                                  scope.rank[y]))
    px, py = plane[x]
    out = {} if held else {first: (px + L * cos(r), py + L * sin(r))}
    for k, y in enumerate(order, 1):
        a = r + 2 * pi * k / n
        out[y] = (px + L * cos(a), py + L * sin(a))
    return out


def _arc_carry(scope, ring, inner, i, j):
    """The atoms an arc `inner` of `ring` between `i` and `j` moves: itself and all it reaches without `i`, `j` or the
    rest of `ring`; None when that reaches the rest of `ring`, so the arc cannot move rigidly."""
    rest = set(ring) - set(inner) - {i, j}
    seen, stack = set(inner), list(inner)
    while stack:
        x = stack.pop()
        for y in scope.adj[x]:
            if y in seen or y == i or y == j:
                continue
            if y in rest:
                return None
            seen.add(y)
            stack.append(y)
    return frozenset(seen)


def arc_moves(scope, plane, front):
    """The arc moves of each ring with an arc end in `front`: up to `ARCS` arcs of at most half the ring, nearest
    the trouble first, arcs as near taken in turn by length, shortest first, then by the ranks of their ends; each a
    crankshaft `arc_flip` across its chord, an `arc_rotate` of `ARC_TURN` either way about its first end, and an
    `arc_stretch` of `STEP` either way along its first bond."""
    rank = scope.rank
    hops = {x: 0 for x in frontier(scope, plane, 0)}        # bonds from the nearest trouble atom
    todo = list(hops)
    while todo:
        nxt = []
        for x in todo:
            for y in scope.adj[x]:
                if y not in hops:
                    hops[y] = hops[x] + 1
                    nxt.append(y)
        todo = nxt
    far = len(scope.atoms)
    out = []
    for ring in scope.rings:
        n = len(ring)
        arcs = []
        for k in range(1, n // 2 + 1):
            for s in range(n):
                inner = tuple(ring[(s + t) % n] for t in range(k))
                i, j = ring[s - 1], ring[(s + k) % n]
                if i not in front and j not in front:
                    continue
                part = _arc_carry(scope, ring, inner, i, j)
                if part is not None:
                    near = min(hops.get(x, far) for x in (i, j, *inner))
                    arcs.append((near, k, sorted((rank[i], rank[j])), i, j, inner, part))
        arcs.sort(key=lambda a: a[:3])
        turn = {}                                           # an arc's turn among those as near and as long
        for a in arcs:
            turn[id(a)] = turn.setdefault((a[0], a[1]), 0)
            turn[a[0], a[1]] += 1
        arcs.sort(key=lambda a: (a[0], turn[id(a)], a[1], a[2]))
        seen = set()
        for *_, i, j, inner, part in arcs:
            if part in seen:
                continue
            seen.add(part)
            if len(seen) > ARCS:
                break
            a = inner[0]
            out.append(Move('arc_flip', i, j, 0., part))
            out += [Move('arc_rotate', i, a, t, part) for t in (ARC_TURN, -ARC_TURN)]
            out += [Move('arc_stretch', i, a, t, part) for t in (STEP, -STEP)]
    return out


def apply(scope, plane, move):
    """`{atom: point}` of the atoms `move` puts elsewhere."""
    if move.kind == 'even':
        return _even(scope, plane, move.pivot, move.part)
    if move.kind == 'arc_flip':
        return _moved(plane, move.part, move.pivot, move.other, (None, 1.))
    if move.kind == 'arc_rotate':
        return _moved(plane, move.part, move.pivot, move.other, (move.arg, 1.))
    if move.kind == 'arc_stretch':
        (ix, iy), (ax, ay) = plane[move.pivot], plane[move.other]
        d = _length(plane, move.pivot, move.other) or 1.
        dx, dy = (ax - ix) / d * move.arg, (ay - iy) / d * move.arg
        return {x: (plane[x][0] + dx, plane[x][1] + dy) for x in move.part}
    if move.kind == 'pivot':
        c, s = cos(move.arg), sin(move.arg)
        ox, oy = plane[move.pivot]
        return {x: (ox + c * (plane[x][0] - ox) - s * (plane[x][1] - oy),
                    oy + s * (plane[x][0] - ox) + c * (plane[x][1] - oy)) for x in move.part}
    if move.kind == 'flip':
        return _moved(plane, move.part, move.pivot, move.other, (None, 1.))
    if move.kind in ('rotate', 'bend'):
        return _moved(plane, move.part, move.pivot, move.other, (move.arg, 1.))
    if move.kind == 'stretch':
        return _moved(plane, move.part, move.pivot, move.other, (0., move.arg / _length(plane, move.pivot, move.other)))
    raise ValueError(f'no primitive move {move.kind}')


def masked(scope, plane, moved, longest=LONGEST):
    """True when `moved` breaks a stated cis/trans, bends a linear atom that `plane` drew straight, or leaves a bond
    between a moved and an unmoved atom shorter than `L` (or than it was) or longer than `longest`."""
    new = {**plane, **moved}
    part = set(moved)
    mol = scope.mol
    for x in part:
        for y in scope.adj[x]:
            if y not in part:
                d = _length(new, x, y)
                if d > longest + 1e-9 or d < min(L, _length(plane, x, y)) - 1e-9:
                    return True
    return broken(scope, plane, moved) > _ct_broken(mol, plane, scope.frames, part) or \
        _bent(mol, new, scope.linear, part) > _bent(mol, plane, scope.linear, part)


def broken(scope, plane, moved):
    """The stated cis/trans units touching `moved` that `plane` with `moved` applied draws against their parity."""
    return _ct_broken(scope.mol, {**plane, **moved}, scope.frames, set(moved))


# --- ring systems ------------------------------------------------------------------------------------------

VARIANTS = 4            # embeddings offered per ring system, mirrors included
UNEVEN = 1.5            # a clean variant's bond spread over the evenest's, beyond which it is dropped


def _block_drawings(mol, atoms, rings, rank):
    """`(axial points, clean)` of each block function that draws the block, in the order the spec ranks them. A
    block holding a ring of more than 8 atoms gets its polygon drawing (`planar_block`) only when no other function
    draws it: a macrocycle is a lattice meander, never a circle, if any meander exists."""
    macro = max(map(len, rings)) > 8
    found = False
    t = template_block(mol, atoms, rings, rank)
    if t is not None:
        found = True
        yield t[0], True
    for pts, bad in embed_block(mol, atoms, rings, rank, k=VARIANTS):
        found = True
        yield pts, not bad
    for f in ((macro_block, bridge_block, sphere_block, cage_block) if macro else
              (bridge_block, planar_block, sphere_block, cage_block, macro_block)):
        p = f(mol, atoms, rings, rank)
        if p is not None:
            found = True
            yield p, True
    if macro and not found:
        p = planar_block(mol, atoms, rings, rank)
        if p is not None:
            yield p, True


def _signature(order, pts):
    """Equal for two drawings of one block that differ by a turn and a shift only."""
    a, b = order[0], order[1]
    (ax, ay), (bx, by) = pts[a], pts[b]
    dist = tuple(round(hypot(pts[x][0] - ax, pts[x][1] - ay), 3) for x in order)
    side = tuple(round((bx - ax) * (pts[x][1] - ay) - (by - ay) * (pts[x][0] - ax), 3) for x in order)
    return dist, side


def variants(scope, block):
    """Up to `VARIANTS` drawings of the ring system `block = (atoms, rings)` in plane coordinates, each followed by its
    mirror image, distinct up to a turn and a shift: the clean drawings, fewest broken cis/trans units of the system
    first, then evenest bonds, one whose bond spread exceeds by `UNEVEN` the evenest's of as many broken units
    dropped; the drawings with conflicts only when none is clean."""
    atoms, rings = block
    order = sorted(atoms, key=scope.rank.__getitem__)
    bonds = [(n, m) for n, m in scope.bonds if n in atoms and m in atoms]
    clean, dirty, seen = [], [], set()
    for pts, ok in _block_drawings(scope.mol, atoms, rings, scope.rank):
        plane = {x: project(p) for x, p in pts.items()}
        for v in (plane, {x: (p[0], -p[1]) for x, p in plane.items()}):
            k = _signature(order, v)
            if k not in seen:
                seen.add(k)
                (clean if ok else dirty).append(v)
        if len(clean) >= 2 * VARIANTS:
            break
    frames = [f for f in scope.frames if set(f[1]) <= atoms]
    keyed = [(_ct_broken(scope.mol, v, frames, atoms), _spread(v, bonds), v) for v in clean]
    least = {}
    for ct, d, _ in keyed:
        least[ct] = min(least.get(ct, d), d)
    keyed = [t for t in keyed if t[1] <= least[t[0]] * UNEVEN]
    clean = [v for _, _, v in sorted(keyed, key=lambda t: (t[0], round(t[1], 2)))]
    return (clean or dirty)[:VARIANTS]


def _spread(plane, bonds):
    """The longest bond of `plane` over its shortest."""
    lengths = [hypot(plane[n][0] - plane[m][0], plane[n][1] - plane[m][1]) for n, m in bonds]
    return max(lengths) / (min(lengths) or 1.) if lengths else 1.


def _pieces(scope, atoms):
    """`[(anchor, piece)]`: the parts of the component left without `atoms`, each with the one atom of `atoms` it
    hangs on, by the rank of its lowest atom."""
    left = set(scope.atoms) - atoms
    out = []
    while left:
        start = min(left, key=scope.rank.__getitem__)
        piece, stack = {start}, [start]
        while stack:
            x = stack.pop()
            for y in scope.adj[x]:
                if y in left and y not in piece:
                    piece.add(y)
                    stack.append(y)
        left -= piece
        anchors = sorted({y for x in piece for y in scope.adj[x] if y in atoms}, key=scope.rank.__getitem__)
        out.append((anchors[0], frozenset(piece)))
    return out


def _toward(plane, a, piece, adj):
    """Direction from `a` to the mean of its neighbours in `piece`."""
    ys = [y for y in adj[a] if y in piece]
    mx = sum(plane[y][0] for y in ys) / len(ys)
    my = sum(plane[y][1] for y in ys) / len(ys)
    return atan2(my - plane[a][1], mx - plane[a][0])


def _turn(p, o, a):
    c, s = cos(a), sin(a)
    x, y = p[0] - o[0], p[1] - o[1]
    return c * x - s * y, s * x + c * y


def embed(scope, plane, atoms, variant, root):
    """`{atom: point}` drawing the ring system `atoms` as `variant`: the piece holding `root` stays, the system turns so
    its bond to that piece leaves on its exit direction, and every other piece turns rigidly onto its anchor's exit."""
    pieces = _pieces(scope, atoms)
    exits = {}                                              # anchor -> its pieces, each with an exit in `variant`
    for a, piece in pieces:
        exits.setdefault(a, []).append(piece)
    out = {a: dict(zip(ps, outside_angles(scope, variant, a, len(ps), atoms))) for a, ps in exits.items()}
    held = next(((a, p) for a, p in pieces if root in p), None)
    if held is None:
        o, phi = root, 0.
    else:
        o, piece = held
        t = out[o].get(piece)
        phi = 0. if t is None else _toward(plane, o, piece, scope.adj) - t
    new = {}
    for x in atoms:
        dx, dy = _turn(variant[x], variant[o], phi)
        new[x] = plane[o][0] + dx, plane[o][1] + dy
    for a, piece in pieces:
        if held is not None and piece is held[1]:
            continue
        t = out[a].get(piece)
        psi = 0. if t is None else t + phi - _toward(plane, a, piece, scope.adj)
        for x in piece:
            dx, dy = _turn(plane[x], plane[a], psi)
            new[x] = new[a][0] + dx, new[a][1] + dy
    return new
