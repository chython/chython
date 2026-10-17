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
"""What a policy reads of an `Env`: atom, pair, global and step features, each invariant to translation, rotation and
atom order. Atoms are listed in canonical rank order, and a pair or a step names atoms by their position in that list.

| atom feature | value |
| --- | --- |
| element | one-hot over `ELEMENTS`, the last class every other element |
| degree | neighbours / 6 |
| ring | in a ring; its smallest ring / 8; its ring systems / 2 |
| flags | linear atom; in a stated cis/trans frame; labelled (not a neutral carbon) |
| reach | atoms within two bonds / 12 |
| strain | the largest deviation of its bond gaps from even, / 180 degrees |
| room | its nearest unbonded atom, / L, at most 3 |
| state | in the frontier; drawn (acyclic, or in a drawn ring system) |

| pair feature | value |
| --- | --- |
| hops | bonds between the two, at most `HOPS`, an integer |
| distance | in the plane, / L |
| bond | 0 none, 1 to 3 its order, 4 aromatic, an integer; ring bond; it crosses another bond |

Pairs are every pair of the component up to `SPARSE` atoms; above it, each atom's pairs within 3 bonds and its
`NEAREST` nearest in the plane.

| global feature | value |
| --- | --- |
| budget | steps taken / cap; ring systems undrawn / systems |
| state | broken cis/trans units; open closures; the `TERMS`, each / atoms |

| step feature | value |
| --- | --- |
| kind | index into `KINDS`, `stop` the last |
| target | a rotation's index into `TARGET_TYPES`; `len(TARGET_TYPES)` for an arc move; one more for any other |
| values | its `delta_terms`, each / atoms; broken units it adds; atoms it moves / atoms; its argument: \
a turn / pi, a stretch / L, a variant's place / `VARIANTS` |
| atoms | the positions of its pivot and other atom (-1 for none) and of every atom it moves |
"""
from heapq import nsmallest
from math import atan2, degrees, hypot, pi

from .edits import KINDS, TARGET_TYPES, VARIANTS, _dir, targets_typed
from .lattice import L
from .relax import _cross
from .score import TERMS, frontier, open_closures


__all__ = ['ELEMENTS', 'atom_features', 'pair_features', 'global_features', 'step_features', 'order']

ELEMENTS = ('C', 'N', 'O', 'S', 'P', 'F', 'Cl', 'Br', 'I', 'B', 'Si', 'Se')
HOPS = 8
SPARSE = 400
NEAREST = 32
ATOM = len(ELEMENTS) + 1 + 12
PAIR = 5
GLOBAL = 4 + len(TERMS)
STEP = 3 + len(TERMS)


def order(env):
    """The component's atoms in canonical rank order, and each atom's position."""
    atoms = env.scope.atoms
    return atoms, {x: k for k, x in enumerate(atoms)}


def _hops(env):
    """Bonds between every two atoms, at most `HOPS`, from a breadth-first search per atom; cached on the scope's
    sides cache under the key `None`."""
    cache = env.scope.sides
    if None not in cache:
        adj = env.scope.adj
        out = {}
        for s in env.scope.atoms:
            dist, frontier_ = {s: 0}, [s]
            for d in range(1, HOPS):
                nxt = []
                for x in frontier_:
                    for y in adj[x]:
                        if y not in dist:
                            dist[y] = d
                            nxt.append(y)
                frontier_ = nxt
            out[s] = dist
        cache[None] = out
    return cache[None]


def _strain(env, x):
    nb = env.scope.adj[x]
    if len(nb) < 2:
        return 0.
    p = env.plane[x]
    angs = sorted(degrees(atan2(env.plane[y][1] - p[1], env.plane[y][0] - p[0])) % 360 for y in nb)
    gaps = [(angs[(i + 1) % len(angs)] - angs[i]) % 360 for i in range(len(angs))]
    even = 180. if x in env.scope.linear else 360. / len(nb)
    return max(abs(g - even) for g in gaps) / 180.


def atom_features(env):
    """One vector of `ATOM` floats per atom, in `order`."""
    scope, plane, mol = env.scope, env.plane, env.mol
    atoms, _ = order(env)
    front = frontier(scope, plane)
    framed = {x for f in scope.frames for x in f[1]}
    drawn = {x for (a, _), d in zip(env.blocks, env.drawn) if d >= 0 for x in a}
    cyclic = {x for (a, _) in env.blocks for x in a}
    small = {}
    for r in scope.rings:
        for x in r:
            small[x] = min(small.get(x, 99), len(r))
    systems = {}
    for a, _ in env.blocks:
        for x in a:
            systems[x] = systems.get(x, 0) + 1
    hops = _hops(env)
    out = []
    for x in atoms:
        a = mol.atom(x)
        e = [0.] * (len(ELEMENTS) + 1)
        e[ELEMENTS.index(a.atomic_symbol) if a.atomic_symbol in ELEMENTS else len(ELEMENTS)] = 1.
        px, py = plane[x]
        room = min((hypot(plane[y][0] - px, plane[y][1] - py) for y in atoms
                    if y != x and y not in scope.adj[x]), default=3 * L) / L
        out.append(e + [
            len(scope.adj[x]) / 6, float(x in small), small.get(x, 0) / 8 if x in small else 0., systems.get(x, 0) / 2,
            float(x in scope.linear), float(x in framed), float(a.atomic_symbol != 'C' or a.charge != 0),
            sum(1 for d in hops[x].values() if 0 < d <= 2) / 12, _strain(env, x), min(room, 3.),
            float(x in front), float(x in drawn or x not in cyclic)])
    return out


def _crossing(env):
    bonds = env.scope.bonds
    plane = env.plane
    out = set()
    for k, (a, b) in enumerate(bonds):
        for c, d in bonds[k + 1:]:
            if len({a, b, c, d}) == 4 and _cross(plane[a], plane[b], plane[c], plane[d]):
                out.add(frozenset((a, b)))
                out.add(frozenset((c, d)))
    return out


def pair_features(env):
    """`(pairs, vectors)`: position pairs `(i, j)`, `i < j`, and one vector of `PAIR` numbers each."""
    scope, plane, mol = env.scope, env.plane, env.mol
    atoms, pos = order(env)
    hops = _hops(env)
    if len(atoms) <= SPARSE:
        pairs = [(i, j) for i in range(len(atoms)) for j in range(i + 1, len(atoms))]
    else:
        keep = set()
        for x in atoms:
            px, py = plane[x]
            keep |= {tuple(sorted((pos[x], pos[y]))) for y, d in hops[x].items() if 0 < d <= 3}
            near = nsmallest(NEAREST + 1, atoms, key=lambda y: (hypot(plane[y][0] - px, plane[y][1] - py), pos[y]))
            keep |= {tuple(sorted((pos[x], pos[y]))) for y in near if y != x}
        pairs = sorted(keep)
    crossing = _crossing(env)
    out = []
    for i, j in pairs:
        x, y = atoms[i], atoms[j]
        b = mol.bond(x, y) if y in scope.adj[x] else None
        order_ = 0 if b is None else (4 if b.order == 4 else min(int(b.order), 3))
        out.append((min(hops[x].get(y, HOPS), HOPS), hypot(plane[x][0] - plane[y][0], plane[x][1] - plane[y][1]) / L,
                    order_, float(b is not None and frozenset((x, y)) in scope.ring_bonds),
                    float(b is not None and frozenset((x, y)) in crossing)))
    return pairs, out


def global_features(env):
    """`GLOBAL` floats of the state."""
    n = len(env.scope.atoms)
    systems = len(env.blocks) or 1
    return [env.steps / env.cap, sum(d == -1 for d in env.drawn) / systems, float(env.ct),
            len(open_closures(env.scope, env.plane)), *(t / n for t in env.terms)]


def step_features(env, step):
    """`(kind, target, values, atoms)` of a legal step: `values` `STEP` floats, `atoms` positions with the pivot and
    other atom first."""
    m = step.move
    _, pos = order(env)
    n = len(pos)
    kind = KINDS.get(m.kind, len(KINDS))
    if m.kind == 'rotate' and not getattr(env, 'primitive', False):    # a primitive turn says its angle, in `arg`
        a = (_dir(env.plane, m.pivot, m.other) + m.arg) % (2 * pi)
        typed = targets_typed(env.scope, env.plane, m.pivot, m.other)
        target = min(typed, key=lambda t: min((t[0] - a) % (2 * pi), (a - t[0]) % (2 * pi)))[1] if typed else \
            len(TARGET_TYPES) + 1
    elif m.kind.startswith('arc_'):
        target = len(TARGET_TYPES)
    else:
        target = len(TARGET_TYPES) + 1
    if m.kind in ('stretch',):
        arg = m.arg / L
    elif m.kind == 'arc_stretch':
        arg = m.arg / L
    elif m.kind in ('rotate', 'bend', 'arc_rotate', 'pivot'):
        arg = m.arg / pi
    elif m.kind == 'embed':
        arg = (m.other + 1) / VARIANTS                  # variants come in canonical order, each before its mirror
    else:
        arg = 0.
    moved = sorted(pos[x] for x in step.moved)
    values = [*(t / n for t in step.terms), float(step.ct - env.ct), len(moved) / n, arg]
    site = [pos.get(m.pivot, -1) if m.kind != 'embed' else -1, pos.get(m.other, -1) if m.kind != 'embed' else -1]
    return kind, target, values, site + moved
