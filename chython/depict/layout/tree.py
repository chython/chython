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
"""The tree embedding: the controller's initial plane, every tree bond at `L`.

| rule | value |
| --- | --- |
| root | the atom of least eccentricity, then of least rank |
| tree | depth first from the root, neighbours in rank order, those continuing the ring of the bond in first; a \
bond back to a placed atom is a closure |
| order | an atom's children placed together, by branch size, largest first, then rank |
| degree 1 to 3 | bonds 120 degrees apart; a chain turns alternately unless the side rule turns it, through a root \
of two unbranched acyclic branches too when no cis/trans is stated: diethyl ether is one zigzag |
| linear atom | its two bonds 180 degrees apart |
| ring atom of two ring bonds, degree 4 or more | ring bonds 120 degrees apart; the rest split the outside 240 evenly |
| metal on an oxo chain | of four or more acyclic bonds, two through an oxygen to another metal: those two on \
one line, the others split over its two sides, evenly within each |
| any other degree n of 4 or more | bonds `360 / n` apart; an acyclic atom of two terminal double-bonded oxygens, \
else of exactly two terminal atoms, puts those two on one line and its other two bonds on the other; a cross of two \
branches and two terminal atoms puts \
the branches on one line, the parent side counted, when it is a carbon of two ring-holding branches (bisphenol A), \
a carbon of two like halogens (a CF2 chain), a sulfonyl or selenonyl, a chain phosphate or a dimethylsilyl; a \
silicon puts two like carbon branches of `TWIN` or more atoms on one line (a di-tert-butylsilyl) |
| side | an atom's children as placed or mirrored about its parent bond: the stated cis/trans kept first, then \
fewest children within `NEAR` of a placed atom or bonds crossing a placed one, a ring atom's children on their \
better side counted too, then least distance from the children to their placed ring-mates, then, below an \
acyclic parent, the heaviest acyclic child trans to the grandparent (a root's heaviest other child) when its \
unbranched run is at most `ZIGZAG` atoms; \
a ring path so curls round its ring and fused hexagons tile the honeycomb |

A closure bond keeps whatever length the tree gives it.
"""
from collections import deque
from math import ceil, cos, floor, pi, sin
from typing import NamedTuple

from .lattice import L
from .relax import _ct_broken


__all__ = ['Grid', 'Tree', 'tree_embedding']

NEAR = .75 * L          # a child this close to a placed atom counts against its side
ZIGZAG = 12             # the longest unbranched run kept trans to its grandparent; a longer one folds to pack


class Tree(NamedTuple):
    plane: dict             # atom -> (x, y)
    parent: dict            # atom -> its tree parent; the root maps to None
    root: int
    closures: list          # (n, m) bonds outside the tree, n of lower rank


def _eccentricity(adj, x):
    seen = {x: 0}
    queue = deque([x])
    while queue:
        y = queue.popleft()
        for z in adj[y]:
            if z not in seen:
                seen[z] = seen[y] + 1
                queue.append(z)
    return max(seen.values())


def _oxo_chain(scope, x):
    """The two neighbours of a metal `x` of four or more acyclic bonds that each bridge to another metal, as the
    oxygens of a polyoxometalate chain do, or None."""
    mol = scope.mol
    if len(scope.adj[x]) < 4 or not mol.atom(x).is_metal or \
            any(frozenset((x, y)) in scope.ring_bonds for y in scope.adj[x]):
        return None
    chain = [y for y in scope.adj[x] if len(scope.adj[y]) == 2 and all(mol.atom(z).is_metal for z in scope.adj[y])]
    return chain if len(chain) == 2 else None


def _chain_slots(n):
    """`(angle, chain)` slots of a metal on an oxo chain: the two chain bonds on one line, the `n - 2` others split
    over its two sides, each side's evenly."""
    up, down = (n - 1) // 2, (n - 2) // 2
    return [(0., True), (pi, True)] + [(pi * (k + 1) / (up + 1), False) for k in range(up)] + \
        [(pi + pi * (k + 1) / (down + 1), False) for k in range(down)]


def _slots(scope, x, n):
    """Relative angles of the `n` bonds of `x`, as `(angle, ring)`: the first ring slot or slot 0 is the reference."""
    ring = [y for y in scope.adj[x] if frozenset((x, y)) in scope.ring_bonds]
    if n == 2 and x in scope.linear:
        return [(0., False), (pi, False)]
    if n <= 3:
        return [(k * 2 * pi / 3, bool(ring)) for k in range(n)]
    if len(ring) == 2:
        rest = n - 2
        return [(0., True), (2 * pi / 3, True)] + \
            [(2 * pi / 3 + 4 * pi / 3 * k / (rest + 1), False) for k in range(1, rest + 1)]
    return [(2 * pi * k / n, False) for k in range(n)]


def tree_embedding(scope):
    """The `Tree` of `scope`'s component: `scope.closures` is left to the caller, which builds the `Scope` again with
    `Tree.closures`."""
    adj, rank = scope.adj, scope.rank
    atoms = scope.atoms
    root = min(atoms, key=lambda x: (_eccentricity(adj, x), rank[x]))
    parent = {root: None}
    children = {x: [] for x in atoms}
    closures = []
    mates = {x: set() for x in atoms}
    for r in scope.rings:
        if r[0] in mates:
            for x in r:
                mates[x].update(r)
    stack = [root]
    seen = set()
    while stack:                                            # iterative: a long chain is no recursion
        x = stack.pop()
        if x in seen:
            continue
        seen.add(x)
        p = parent[x]
        ahead = sorted(adj[x], key=lambda y: (not (p is not None and y in mates[x] and y in mates[p] and
                                                   _share(scope.rings, (p, x, y))), rank[y]))
        for y in reversed(ahead):                           # the ring of the bond in walks on first
            if y not in seen:
                parent[y] = x
                stack.append(y)
    for x in atoms:
        if parent[x] is not None:
            children[parent[x]].append(x)
    tree = {frozenset((x, p)) for x, p in parent.items() if p is not None}
    for n, m in scope.bonds:
        if frozenset((n, m)) not in tree:
            closures.append((n, m))
    size, ringed = {}, {}                               # subtree atoms, and of them ring atoms
    in_ring = {x for b in scope.ring_bonds for x in b}
    down = {}                                           # atoms of the unbranched acyclic run from an atom down
    for x in reversed(_order(root, children)):
        size[x] = 1 + sum(size[y] for y in children[x])
        ringed[x] = (x in in_ring) + sum(ringed[y] for y in children[x])
        down[x] = 0 if x in in_ring else 1 + (down[children[x][0]] if len(children[x]) == 1 else 0)
    run = {}                                            # atoms of the whole unbranched acyclic run through an atom
    for x in _order(root, children):
        p = parent[x]
        run[x] = run[p] if p is not None and p not in in_ring and x not in in_ring and len(children[p]) == 1 and \
            len(adj[p]) == 2 else down[x]

    def bare(y):                                            # an unbranched acyclic chain from `y` on
        while True:
            if y in in_ring or len(children[y]) > 1:
                return False
            if not children[y]:
                return True
            y = children[y][0]

    plane = {root: (0., 0.)}
    back = {root: pi}                                       # direction from an atom to its parent
    turn = {root: 1}
    frames = scope.frames
    grid = Grid()
    grid.add(root, plane[root])
    queue = deque([root])
    while queue:
        x = queue.popleft()
        kids = sorted(children[x], key=lambda y: (-size[y], rank[y]))
        if not kids:
            continue
        angles = _place(scope, x, kids, back[x], turn[x], parent[x],
                        [(len(scope.atoms) - size[x], len(in_ring) - ringed[x])] + [(size[y], ringed[y]) for y in kids])
        mirrored = [2 * back[x] - a for a in angles]
        near_atoms = grid.atoms(plane[x], 2 * L + NEAR)     # every atom a child or a grandchild can come near
        near_bonds = grid.bonds(plane[x], 3 * L)
        best = None
        for k, side in enumerate((angles, mirrored)):
            _put(plane, x, kids, side)
            here = [f for f in frames if f[2] <= plane.keys()]  # every substituent placed: `near` follows input order
            ct = _ct_broken(scope.mol, plane, here, set(kids) | {x}) if here else 0
            near = sum(1 for y in kids for z in near_atoms if z != x and _dist(plane, y, z) < NEAR) + \
                sum(1 for y in kids for a, b in near_bonds if x not in (a, b) and _crosses(plane, x, y, a, b))
            pull = sum(_dist(plane, y, z) for y in kids for z in mates[y] if z in plane and z != x and z not in kids)
            ahead = sum(_ahead(scope, plane, near_atoms, near_bonds, x, kids, y, a + pi,
                               -turn[x] if j == 0 else turn[x],
                               sorted(children[y], key=lambda g: (-size[g], rank[g])))
                        for j, (y, a) in enumerate(zip(kids, side)) if mates[y])   # a ring path looks ahead
            cis = 0                                     # the heaviest branch trans to the grandparent: a zigzag
            g = parent[x]
            if g is not None and g not in in_ring and not mates[kids[0]] and run[kids[0]] <= ZIGZAG:
                ref = parent[g] if parent[g] is not None else next(          # through a root: its heaviest other
                    (y for y in sorted(children[g], key=lambda y: (-size[y], rank[y])) if y != x), None)   # child
                if ref is not None:
                    (gx, gy), (px, py), (qx, qy), (kx, ky) = plane[ref], plane[g], plane[x], plane[kids[0]]
                    cis = ((qx - px) * (gy - py) - (qy - py) * (gx - px)) * \
                        ((qx - px) * (ky - py) - (qy - py) * (kx - px)) > 0
            score = (ct, near + ahead, round(pull, 9), cis, k)
            if best is None or score < best[0]:
                best = (score, side)
        angles = best[1]
        _put(plane, x, kids, angles)
        for y in kids:
            grid.add(y, plane[y])
            grid.link(x, y, plane[x], plane[y])
        chain = parent[x] is None and len(kids) == 2 and not frames and all(map(bare, kids))   # one zigzag through
        for k, (y, a) in enumerate(zip(kids, angles)):
            back[y] = a + pi
            turn[y] = (1, -1)[k] * turn[x] if chain else -turn[x] if k == 0 else turn[x]
            queue.append(y)
    return Tree(plane, parent, root, closures)


def _order(root, children):
    out, stack = [], [root]
    while stack:
        x = stack.pop()
        out.append(x)
        stack.extend(children[x])
    return out


def _ahead(scope, plane, atoms, bonds, x, kids, y, back, turn, children):
    """The fewest troubles `y`'s children meet on either side: atoms within `NEAR` and bonds crossing, against the
    placed `atoms` and `bonds` near `x` and `x`'s `kids` just placed."""
    if not children:
        return 0
    angles = _place(scope, y, children, back, turn, x)
    best = None
    for side in (angles, [2 * back - a for a in angles]):
        spot = [(plane[y][0] + L * cos(a), plane[y][1] + L * sin(a)) for a in side]
        n = 0
        for q in spot:
            n += sum(1 for z in atoms if ((q[0] - plane[z][0]) ** 2 + (q[1] - plane[z][1]) ** 2) ** .5 < NEAR)
            n += sum(1 for z in kids if z != y and ((q[0] - plane[z][0]) ** 2 + (q[1] - plane[z][1]) ** 2) ** .5 <
                     NEAR)
            plane['?'] = q
            n += sum(1 for a, b in bonds if y not in (a, b) and _crosses(plane, y, '?', a, b))
            n += sum(1 for z in kids if z != y and _crosses(plane, y, '?', x, z))
        del plane['?']
        best = n if best is None else min(best, n)
    return best


class Grid:
    """Atoms and bonds by square cell of side `L`: the ones near a point, without a scan of all."""
    __slots__ = ('_atoms', '_bonds')

    def __init__(self):
        self._atoms, self._bonds = {}, {}

    @staticmethod
    def _cell(p):
        return floor(p[0] / L), floor(p[1] / L)

    def add(self, x, p):
        self._atoms.setdefault(self._cell(p), []).append(x)

    def link(self, a, b, pa, pb):
        """Hold the bond `a`-`b`, drawn from `pa` to `pb`, in the cells of both ends."""
        for c in {self._cell(pa), self._cell(pb)}:
            self._bonds.setdefault(c, []).append((a, b))

    def _around(self, table, p, r):
        cx, cy = floor(p[0] / L), floor(p[1] / L)
        k = ceil(r / L)
        out = []
        get = table.get
        for i in range(cx - k, cx + k + 1):
            for j in range(cy - k, cy + k + 1):
                found = get((i, j))
                if found:
                    out += found
        return out

    def atoms(self, p, r):
        """The atoms in cells within `r` of `p`: a superset of those within `r`."""
        return self._around(self._atoms, p, r)

    def bonds(self, p, r):
        """The bonds with an end in a cell within `r` of `p`, each once."""
        return set(self._around(self._bonds, p, r))


def _share(rings, atoms):
    return any(all(x in r for x in atoms) for r in rings)


def _crosses(plane, a, b, c, d):
    """True when bonds `a`-`b` and `c`-`d`, sharing no atom, cross properly."""
    if a == c or a == d or b == c or b == d:
        return False
    (ax, ay), (bx, by), (cx, cy), (dx, dy) = plane[a], plane[b], plane[c], plane[d]
    if max(ax, bx) < min(cx, dx) or max(cx, dx) < min(ax, bx) or max(ay, by) < min(cy, dy) or \
            max(cy, dy) < min(ay, by):                  # boxes apart: no crossing, at four comparisons
        return False
    ux, uy, vx, vy = bx - ax, by - ay, dx - cx, dy - cy
    if (ux * (cy - ay) - uy * (cx - ax)) * (ux * (dy - ay) - uy * (dx - ax)) >= -1e-12:
        return False
    return (vx * (ay - cy) - vy * (ax - cx)) * (vx * (by - cy) - vy * (bx - cx)) < -1e-12


def _dist(plane, a, b):
    return ((plane[a][0] - plane[b][0]) ** 2 + (plane[a][1] - plane[b][1]) ** 2) ** .5


def _ring_cross(sizes):
    """True for a cross of two branches holding a ring and two single terminal atoms, as in bisphenol A."""
    sizes = [(s, r) for s, r in sizes if s]             # a root has no parent side
    return sorted(r > 0 for _, r in sizes) == [False, False, True, True] and \
        sorted(s for s, r in sizes if not r)[-1] == 1


HALOGENS = ('F', 'Cl', 'Br', 'I')


def _dioxo(scope, x):
    """The two terminal oxygens doubly bonded to `x`, in rank order, or none: a sulfonyl, a selenonyl."""
    mol = scope.mol
    out = [y for y in scope.adj[x] if len(scope.adj[y]) == 1 and mol.atom(y).atomic_symbol == 'O' and
           mol.bond(x, y).order == 2]
    return out if len(out) == 2 else []


def _ends(scope, x):
    """The two terminal neighbours of `x`, in rank order, when it has exactly two: a gem-dimethyl, a CF2."""
    out = [y for y in scope.adj[x] if len(scope.adj[y]) == 1]
    return out if len(out) == 2 else []


def _pair_cross(scope, x, sizes):
    """True for a carbon of two branches and two terminal atoms of one halogen, a CF2 or CCl2 of a chain, a sulfur,
    selenium, tellurium or phosphorus of two branches and two terminal oxygens, a sulfonyl, a selenonyl or a phosphate
    of a chain, or a silicon of four carbons, two of them terminal, a dimethylsilyl of a carbon chain."""
    mol = scope.mol
    ends = [mol.atom(y).atomic_symbol for y in scope.adj[x] if len(scope.adj[y]) == 1]
    if len(ends) != 2 or ends[0] != ends[1] or sorted(s > 1 for s, _ in sizes if s) != [False, False, True, True]:
        return False
    here = mol.atom(x).atomic_symbol
    return here == 'C' and ends[0] in HALOGENS or here in ('S', 'Se', 'Te', 'P') and ends[0] == 'O' or \
        here == 'Si' and all(mol.atom(y).atomic_symbol == 'C' for y in scope.adj[x])


TWIN = 4                # the least atoms of each of two like branches a cross puts on one line


def _twins(scope, kids, sizes):
    """`kids` with two carbon branches of one size of `TWIN` or more atoms, the largest such pair, on a cross's
    opposite free slots, the first and the third: a di-tert-butyl silicon's two tert-butyls on one line."""
    by = {}
    for y, (k, r) in zip(kids, sizes):
        if k >= TWIN and not r and scope.mol.atom(y).atomic_symbol == 'C':
            by.setdefault(k, []).append(y)
    pairs = [(k, v) for k, v in by.items() if len(v) == 2]
    if not pairs:
        return kids
    a, b = max(pairs)[1]
    rest = [y for y in kids if y not in (a, b)]
    return [a, rest[0], b] + rest[1:]


def _cross(kids, sizes, parented):
    """`kids` reordered onto a cross's free slots, the first of which is the parent's neighbour, then its opposite,
    then the other neighbour: the two heaviest branches, the parent side counted, on one line."""
    up, sz = sizes[0], dict(zip(kids, sizes[1:]))
    if len(set(sz.values())) == 1:                  # equal branches: no pair to choose, the slots stand
        return kids
    heavy = sorted(kids, key=lambda y: -sz[y])
    if not parented:                                    # slots 0 to 3 free: heavy at 0 and 2, light at 1 and 3
        return [heavy[0], heavy[2], heavy[1]] + heavy[3:]
    if up >= sz[heavy[1]]:                              # the heaviest child opposite the parent
        return [heavy[1], heavy[0], heavy[2]]
    return [heavy[0], heavy[2], heavy[1]]              # the two heaviest children opposite each other


def _put(plane, x, kids, angles):
    px, py = plane[x]
    for y, a in zip(kids, angles):
        plane[y] = (px + L * cos(a), py + L * sin(a))


def _place(scope, x, kids, back, turn, parent, sizes=None):
    """Absolute angles of `kids` from `x`: the slot set of `x` turned so the parent's slot points `back`. With
    `sizes`, `(atoms, ring atoms)` of the parent side and of each of `kids`, a carbon cross of two ring-holding branches
    and two terminal atoms puts the ring branches on one line and the terminal atoms on the other."""
    n = len(scope.adj[x])
    chain = _oxo_chain(scope, x)
    if chain is not None:                                   # chain slots stand where ring slots would
        slots = _chain_slots(n)
        ring_of = {y: y in chain for y in scope.adj[x]}
    else:
        slots = _slots(scope, x, n)
        ring_of = {y: frozenset((x, y)) in scope.ring_bonds for y in scope.adj[x]}
    given = kids
    oxo = _dioxo(scope, x) or _ends(scope, x) if n == 4 and len(kids) >= 3 and not any(ring_of.values()) else []
    if oxo and parent not in oxo:                           # O=X=O or two end atoms on one line, the rest on the other
        rest = [y for y in kids if y not in oxo]
        kids = [oxo[0], rest[0], oxo[1]] if parent is not None else [rest[0], oxo[0], rest[1], oxo[1]]
    elif sizes is not None and n == 4 and len(kids) >= 3 and not any(ring_of.values()) and \
            (scope.mol.atom(x).atomic_symbol == 'C' and _ring_cross(sizes) or _pair_cross(scope, x, sizes)):
        kids = _cross(kids, [s for s, _ in sizes], parent is not None)
    elif sizes is not None and n == 4 and len(kids) >= 3 and not any(ring_of.values()) and \
            scope.mol.atom(x).atomic_symbol == 'Si':
        kids = _twins(scope, kids, sizes[1:])
    free = list(range(n))
    if parent is not None:                                  # the parent takes its kind's first slot
        k = next((k for k in free if slots[k][1] == ring_of[parent]), free[0])
        free.remove(k)
        base = back - slots[k][0]
    else:
        base = 0.
    others = [y for y in scope.adj[x] if y != parent and y not in kids]     # closure partners keep their slots
    angles = {}
    for y in kids + others:                                 # a ring neighbour to a ring slot when one is free
        k = next((k for k in free if slots[k][1] == ring_of[y]), free[0])
        free.remove(k)
        angles[y] = base + slots[k][0]
    out = [angles[y] for y in given]
    if n <= 3 and parent is not None and turn < 0:          # a chain turns alternately: mirror about the parent
        out = [2 * back - a for a in out]
    return out
