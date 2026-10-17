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
"""The controller's score of a plane, compared as a tuple: lower is better.

| place | term |
| --- | --- |
| 0 | defects: the `count_defects` total of close pairs, crossings and atoms on bonds |
| 1 | closure: the sum over closure bonds of their distance from `L` |
| 2 | conventions: a crowded pair 1, a stretched bond 1 (a closure excepted), a bent acyclic atom or a straight ring \
atom `BEND`, a ring of at most 8 atoms `RING` per L of RMSD from its regular polygon, a main axis off horizontal \
`TILT` per 90 degrees; a pair crowds nearer than `CROWD`, neither bonded nor sharing a neighbour |

A `delta` scores only the atoms a move touches and their surroundings within `relax._near`; added to the parent's
score it gives the child's. `key` rounds to `1e-9`, so float noise never orders two planes.

| frontier | atoms |
| --- | --- |
| trouble | `defect_atoms`, both atoms of a crowded pair, both ends of an open closure, the atoms of a broken stated \
cis/trans |
| frontier | trouble and every atom within 2 bonds of it |
"""
from array import array
from math import atan2, cos, degrees, fsum, hypot, radians

from .lattice import L
from .quality import _polygon_rmsd, linear_atoms
from .relax import CLOSE, ON_BOND, SHORT, _bonds, _ct_broken, _frames, _near, _side, count_defects, defect_atoms
from .rings import ring_blocks


BRIDGED = ('P', 'S', 'Si')   # an oxygen between two of these or two metals is drawn straight: a polyphosphate, a
                             # polyoxometalate chain


def bridged(mol):
    """The oxygens of two neighbours, both of `BRIDGED` or both metals, and the silicons of two silicon neighbours,
    a polysilane; a `Scope` draws the acyclic ones straight."""
    out = set()
    for x in mol:
        ns = mol.neighbors_of(x)
        if len(ns) != 2:
            continue
        here = mol.atom(x).atomic_symbol
        if here == 'O' and (all(mol.atom(y).atomic_symbol in BRIDGED for y in ns) or
                            all(mol.atom(y).is_metal for y in ns)) or \
                here == 'Si' and all(mol.atom(y).atomic_symbol == 'Si' for y in ns):
            out.add(x)
    return out


__all__ = ['Scope', 'TERMS', 'terms', 'delta_terms', 'combine', 'full', 'delta', 'add', 'key', 'better', 'frontier',
           'open_closures']

CROWD = 1.2 * L         # two unbonded atoms nearer than this crowd
STRETCHED = 1.1 * L     # a bond, a closure excepted, longer than this is stretched
OFF = .1 * L            # a closure bond further than this from L is open
BEND = 2.5              # a bent atom weighs more than two stretched bonds
BENT = 5.               # degrees off 120 from which an acyclic atom of two or three neighbours is bent
RING = 4.               # a ring of at most 8 atoms, per L of RMSD from its regular polygon
STRAIGHT = 170.         # degrees from which a ring atom of two neighbours, not linear by its bonds, is straight
TILT = 2.               # a drawing's main axis at 90 degrees from horizontal; at 45, half
LONG = 1.5              # spread along the main axis over across it from which a drawing has one
FLAT = 15.              # degrees of tilt that still count as horizontal: half a 30-degree turn
RADIUS = 2              # bonds from the trouble the frontier reaches
TERMS = ('defects', 'closure', 'crowd', 'stretched', 'bent', 'ring', 'straight', 'tilt')


class Scope:
    """What a component's score and edits read and never change."""
    __slots__ = ('mol', 'atoms', 'rank', 'bonds', 'acyclic', 'ring_bonds', 'corners', 'closures', 'frames', 'linear',
                 'adj', 'sides', 'within2', 'rings', 'drawn_bonds', 'kernel', 'positions', 'orient')

    def __init__(self, mol, atoms, rank, closures):
        self.mol = mol
        self.rank = rank
        self.atoms = sorted(atoms, key=rank.__getitem__)
        whole = {x: None for x in atoms}
        self.bonds = _bonds(mol, whole, rank)
        self.acyclic = {b for b in self.bonds if _side(mol, *b) is not None}
        self.rings = [r for _, rs in ring_blocks(mol) for r in rs]   # canonical: an SSSR of a bridged system is not
        self.ring_bonds = {frozenset((r[i - 1], r[i])) for r in self.rings for i in range(len(r))}
        self.corners = {(r[i], frozenset((r[i - 1], r[(i + 1) % len(r)]))) for r in self.rings for i in range(len(r))}
        self.closures = list(closures)
        open_ = {frozenset(b) for b in self.closures}
        self.drawn_bonds = [b for b in self.bonds if frozenset(b) not in open_]    # every bond but a closure
        self.frames = _frames(mol, whole)
        self.linear = (linear_atoms(mol) | {x for x in bridged(mol) if not any(
            frozenset((x, y)) in self.ring_bonds for y in mol.neighbors_of(x))}) & set(atoms)
        self.adj = {x: sorted(mol.neighbors_of(x), key=rank.__getitem__) for x in atoms}
        self.sides = {}                                     # (n, m) -> atoms reached from m without the bond
        self.within2 = {x: {z for y in self.adj[x] for z in (y, *self.adj[y])} for x in atoms}
        self.kernel = self.positions = None             # the compiled score, built on first use
        self.orient = False                             # `tilt` counts: a drawing whose orientation is its own


def _length(plane, n, m):
    return hypot(plane[n][0] - plane[m][0], plane[n][1] - plane[m][1])


def _bent(scope, plane, x):
    nb = scope.adj[x]
    if len(nb) not in (2, 3) or x in scope.linear or any(frozenset((x, y)) in scope.ring_bonds for y in nb):
        return False
    px, py = plane[x]
    angs = sorted(degrees(atan2(plane[y][1] - py, plane[y][0] - px)) % 360 for y in nb)
    gaps = [(angs[(i + 1) % len(angs)] - angs[i]) % 360 for i in range(len(angs))]
    dev = abs(min(gaps) - 120) if len(nb) == 2 else max(abs(g - 120) for g in gaps)
    return dev > BENT


def _straight(scope, plane, x):
    """A ring atom of two neighbours, not linear by its bonds, drawn on its neighbours' line."""
    nb = scope.adj[x]
    if len(nb) != 2 or x in scope.linear or not any(frozenset((x, y)) in scope.ring_bonds for y in nb):
        return False
    (px, py), (ax, ay), (bx, by) = plane[x], plane[nb[0]], plane[nb[1]]
    u = hypot(ax - px, ay - py) * hypot(bx - px, by - py)
    return bool(u) and ((ax - px) * (bx - px) + (ay - py) * (by - py)) / u < cos(radians(STRAIGHT))


def _crowded(scope, plane, part, atoms):
    """Pairs nearer than `CROWD` with an atom in `part`, each once, two atoms bonded or sharing a neighbour excepted:
    the ligands of one centre sit at its slots, `L` apart at 60 degrees."""
    rank = scope.rank
    for a in part:
        xa, ya = plane[a]
        skip = scope.within2[a]
        for b in atoms:
            if b in skip or b in part and rank[b] <= rank[a]:
                continue
            xb, yb = plane[b]
            if abs(xa - xb) < CROWD and abs(ya - yb) < CROWD and hypot(xa - xb, ya - yb) < CROWD:
                yield a, b


def _parts(scope, plane, part, atoms, bonds):
    defects = sum(count_defects(plane, atoms, bonds))       # every pair of the subset: a pair far from `part` cancels
    closure = fsum(abs(_length(plane, n, m) - L) for n, m in scope.closures if n in part or m in part)
    crowd = sum(1 for _ in _crowded(scope, plane, part, atoms))
    stretched = sum(_length(plane, n, m) > STRETCHED for n, m in scope.drawn_bonds if n in part or m in part)
    touched = part | {y for x in part for y in scope.adj[x]}
    bent = sum(_bent(scope, plane, x) for x in touched)
    ring = fsum(_polygon_rmsd([plane[x] for x in r]) for r in scope.rings if len(r) <= 8 and not part.isdisjoint(r))
    straight = sum(_straight(scope, plane, x) for x in touched)
    return defects, closure, crowd, stretched, bent, ring, straight


def combine(terms):
    """The score of the `TERMS`: defects, closure, conventions."""
    defects, closure, crowd, stretched, bent, ring, straight, tilt = terms
    return defects, closure, crowd + stretched + BEND * (bent + straight) + RING * ring + TILT * tilt


def tilt(plane):
    """The main axis's angle from horizontal beyond `FLAT`, as a fraction of the rest of 90 degrees, of a drawing
    spread `LONG` times more along it than across; 0 for a rounder one."""
    pts = list(plane.values())
    n = len(pts)
    if n < 3:
        return 0.
    cx, cy = fsum(p[0] for p in pts) / n, fsum(p[1] for p in pts) / n
    sxx = fsum((p[0] - cx) ** 2 for p in pts)
    syy = fsum((p[1] - cy) ** 2 for p in pts)
    sxy = fsum((p[0] - cx) * (p[1] - cy) for p in pts)
    half, det = (sxx + syy) / 2, sxx * syy - sxy * sxy
    spread = max(0., half * half - det) ** .5
    if half + spread <= LONG * max(half - spread, 1e-12):
        return 0.
    a = abs(degrees(.5 * atan2(2 * sxy, sxx - syy))) % 180
    return max(0., min(a, 180 - a) - FLAT) / (90 - FLAT)


def _tilt(scope, plane):
    return tilt(plane) if scope.orient else 0.


def _oracle_terms(scope, plane):
    """`terms` in Python: the reference the kernel is tested against."""
    return _parts(scope, plane, set(scope.atoms), scope.atoms, scope.bonds) + (_tilt(scope, plane),)


def _oracle_delta(scope, plane, moved):
    """`delta_terms` in Python: the reference the kernel is tested against."""
    part = set(moved)
    new = {**plane, **moved}
    atoms, bonds = _near((plane, new), scope.atoms, scope.bonds, part)
    a = _parts(scope, plane, part, atoms, bonds) + (_tilt(scope, plane),)
    b = _parts(scope, new, part, atoms, bonds) + (_tilt(scope, new),)
    return tuple(y - x for x, y in zip(a, b))


def _kernel(scope):
    """The scope's compiled `LayoutScore`, built once: atoms are positions in `scope.atoms`."""
    if scope.kernel is None:
        from ...core._core import LayoutScore
        pos = {x: k for k, x in enumerate(scope.atoms)}
        n = len(scope.atoms)

        def flat(pairs):
            return array('i', [pos[a] for p in pairs for a in p])

        def csr(lists):
            off, out = array('i', [0]), array('i')
            for xs in lists:
                out.extend(xs)
                off.append(len(out))
            return off, out
        adj_off, adj = csr([[pos[y] for y in scope.adj[x]] for x in scope.atoms])
        near_off, near = csr([[pos[y] for y in scope.within2[x]] for x in scope.atoms])
        ring_off, ring = csr([[pos[x] for x in r] for r in scope.rings if len(r) <= 8 and all(x in pos for x in r)])
        bent = bytearray(n)
        straight = bytearray(n)
        for x in scope.atoms:
            nb = scope.adj[x]
            in_ring = any(frozenset((x, y)) in scope.ring_bonds for y in nb)
            bent[pos[x]] = len(nb) in (2, 3) and x not in scope.linear and not in_ring
            straight[pos[x]] = len(nb) == 2 and x not in scope.linear and in_ring
        scope.kernel = LayoutScore(n, flat(scope.bonds), flat(scope.closures), flat(scope.drawn_bonds), adj_off, adj,
                                   near_off, near, memoryview(bent), memoryview(straight), ring_off, ring,
                                   (L, CLOSE, SHORT, ON_BOND, CROWD, STRETCHED, BENT, cos(radians(STRAIGHT)), 2 * L))
        scope.positions = pos
    return scope.kernel


def _xy(scope, plane):
    return array('d', [c for x in scope.atoms for c in plane[x]])


def kernel_terms(scope, plane):
    """`terms` by the compiled kernel, `tilt` by Python."""
    return _kernel(scope).terms(_xy(scope, plane)) + (_tilt(scope, plane),)


def kernel_delta(scope, plane, moved):
    """`delta_terms` by the compiled kernel."""
    k = _kernel(scope)
    pos = scope.positions
    d = k.delta(_xy(scope, plane), array('i', [pos[x] for x in moved]),
                array('d', [c for p in moved.values() for c in p]))
    return d + ((tilt({**plane, **moved}) - tilt(plane)) if scope.orient else 0.,)


terms = kernel_terms
delta_terms = kernel_delta


def full(scope, plane):
    """The score of `plane` over the whole component."""
    return combine(terms(scope, plane))


def delta(scope, plane, moved):
    """The change of `full` when the points of `moved` replace their atoms' points in `plane`."""
    return combine(delta_terms(scope, plane, moved))


def add(score, d):
    return tuple(x + y for x, y in zip(score, d))


def key(score):
    return score[0], round(score[1], 9), round(score[2], 9)


def better(a, b):
    """`a` scores strictly lower than `b`."""
    return key(a) < key(b)


def open_closures(scope, plane):
    return [(n, m) for n, m in scope.closures if abs(_length(plane, n, m) - L) > OFF]


def frontier(scope, plane, radius=RADIUS):
    """Atoms within `radius` bonds of a defect, a crowded pair, an open closure or a stated cis/trans drawn against its
    parity: the trouble itself at 0."""
    near = set(defect_atoms(plane, scope.atoms, scope.bonds))
    for f in scope.frames:
        if f[2] <= plane.keys() and _ct_broken(scope.mol, plane, [f], f[2]):
            near |= f[2]
    for a, b in _crowded(scope, plane, set(scope.atoms), scope.atoms):
        near.update((a, b))
    for n, m in open_closures(scope, plane):
        near.update((n, m))
    for _ in range(radius):
        near |= {y for x in near for y in scope.adj[x]}
    return near
