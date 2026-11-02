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
"""The cascade layout: templated ring systems on a honeycomb tree, then a deterministic fit, then the pose.

| stage | rule |
| --- | --- |
| start | `primitive.Initial`: the tree, each ring system's `system_drawing`, substituent exits, boundary bonds at `L` |
| fit | while a defect, a broken stated cis/trans or an atom inside a ring is left, up to `FIT` moves: of the flips, \
anchored rotations and stretches of tree bonds whose subtree holds some atoms of the trouble but not all, none \
putting an atom inside a ring, a substituent `in_ring_angle` of its atom, two bonds of an acyclic atom or, but in \
the fit allowing atoms off even, two substituents of a ring atom within `PINCH`, the move repairing most cis/trans \
units, then a stretch along a linear atom of no metal, then removing most defects, then bending fewest substituents \
within `EXIT` of their atom's ring bonds, then squeezing fewest pairs of a ring atom's two substituents (see \
`_squeezed`), a bend or squeeze kept only for a defect it clears, then leaving fewest atoms off even (see `_off`), \
then taking most atoms out of rings, then not stretching a bond to a terminal atom; ties to the kind in `KIND` \
order, then the smaller move, then a stretch of a bond between two sp3 atoms, then the step order; no move breaks a \
stated cis/trans or a linear atom, and no stretch leaves a bond past `STRETCHED` |
| unstick | a defect, an atom off even or a long bond left: a flip, bend or stretch taken at a bond joining two \
massive fragments, or a cheap one near the trouble, when the repair after it leaves a better scene (see `_unstick`) |
| settle | a ring substituent inside a ring, bent onto it or off its bisector put on the bisector, the rest repaired \
round it (see `_settle`) |
| unload | nothing wrong left: flips opening the most locally overloaded corners (see `_unload`) |
| nudge | a defect left: an atom of it moved by one of `NUDGES` when that removes a defect |
| pose | `level`, `horizontal`, `lean` of the ortho-fused systems, `straight` when acyclic, `upright`; a component \
whose largest system is a shape keeps the shape's orientation |
| components | left to right by rank, `pose.GAP` apart |

Covered and frozen: a molecule whose every ring system is a shape, one ring of up to `TEMPLATE` atoms, or rings of
up to `TEMPLATE` atoms sharing edges only is `covered`; the drawing of each covered corpus row is held by its
`drawing_hash` in `GOLDEN`, and a test fails on any change until `bless` writes the file again:

    python -m chython.depict.layout.cascade bless

The ring systems: each drawn once, by kind, its drawing and mirror image in plane coordinates.

| kind | drawing |
| --- | --- |
| shape | `shapes.fused3d`: a predefined cage, rings fused to it on one edge each as regular polygons |
| circle | a crown `cycles` calls a `Circle`: `templates.crown`, its edge-fused rings regular polygons outward; a ring \
of `MEDIUM` atoms with smaller rings fused on its edges: the regular polygon, those outward (`_medium`) |
| ring | one ring of up to `TEMPLATE` atoms: its regular polygon |
| system | rings of up to `TEMPLATE` atoms: the first clean `edits._block_drawings` drawing, else the first; an \
ortho-fused system's the most regular of the first `REGULAR` clean ones, its rings evened by `relax.regularize` |
| | an ortho-fused macrocycle's: its rings of up to `TEMPLATE` atoms evened by `relax.regularize` too, unless that \
adds a defect |
| units | a macrocycle holding a ring of up to `TEMPLATE` atoms: the first clean drawing, if it has no defect of \
its own |
| tree | any other macrocycle: None, so it rides the tree |

A system or units drawing is a function of the system's key alone: its all-carbon skeleton with one stub per bond
leaving it. It is computed once on the key renumbered in its canonical order, that order deciding every tie, and
cached;
a molecule gets it through the skeleton mapping whose canonical ranks come first, so no drawing depends on atom
order or on which molecule was drawn before.
"""
from argparse import ArgumentParser
from hashlib import sha1
from itertools import islice
from math import atan2, cos, hypot, pi, sin, tan
from pathlib import Path

from ...core import MoleculeContainer
from .lattice import L, project
from .relax import count_defects
from .shapes import fused3d, shape
from .tree import Grid


__all__ = ['TEMPLATE', 'FIT', 'STRETCHED', 'GOLDEN', 'system_drawing', 'covered', 'fit', 'cascade_layout',
           'drawing_hash', 'bless', 'in_ring_angle']

TEMPLATE = 9            # the largest ring a polygon or a system drawing covers
MAPPINGS = 48           # skeleton symmetries tried when mapping a cached drawing onto a molecule
CACHE = 4096
REGULAR = 8             # clean system drawings compared for the most regular
LONGEST_BOND = 2 * L    # the longest bond a nudge of a bridged system's drawing leaves
SPACE = .7 * L          # two atoms nearer than this, neither bonded nor geminal, are spaced by a flip or stretch
FIT = 32                # fitting moves at most per component
STRETCHED = 3 * L       # the longest a fitting stretch leaves a bond
EXIT = pi / 2 - 1e-6    # the least angle a fitting move leaves between a substituent and its atom's ring bonds
PINCH = pi / 3 - 1e-6   # the least angle a fitting move leaves between two bonds of an acyclic atom, or two
                        # substituents of a ring atom
LARGE = 300             # heavy atoms of a component past which smiles-drawer draws it, not the fit
CRAMP = pi / 8          # a substituent's bond this near a ring bond of its atom is folded onto the ring, unless
                        # the atom has two of each: then `PINCH`
EVEN = pi / 36          # how far below `360 / n` a gap of an acyclic atom of n bonds is off even
UNSTICK = 8             # flips and bends the fit tries with a repair after each, when stuck
AFTER = 8               # repair moves after such a move
MASSIVE = 6             # atoms on each side of a bridge the stuck fit moves whole
REACH = 2               # steps a kick's stretch at a bridge takes at most
WIDE = 12               # kicks at bridges tried per round
POLISH = 4              # rounds of kicks once no defect is left
CROWDED = .99 * L       # unbonded atoms nearer than the lattice's 1 L overload a corner
LOOSE = 1.5 * L         # and nearer than this, a tie between two flips decided
UNLOAD = .01            # the least overload, in L squared, a flip of `_unload` removes
SETTLE = 8              # substituents a component's settle tries
SMALL_SETTLE = 40       # the most atoms a settled substituent holds
NUDGES = tuple((r * cos(k * pi / 6), r * sin(k * pi / 6)) for r in (L / 4, L / 2) for k in range(12))
KIND = {'flip': 0, 'rotate': 1, 'stretch': 2}
GOLDEN = Path(__file__).resolve().parents[3] / 'test' / 'layout_golden.tsv'
DIGITS = 4              # decimals of `L` a drawing hash keeps: platform noise below them is not a change
MEDIUM = (10, 12)       # a fused ring this size is a regular polygon, off the lattice


def system_drawing(mol, atoms, rings, rank):
    """`(kind, [drawing, mirror image])` of the ring system `atoms` with `rings`, each drawing `{atom: (x, y)}`, or
    None to leave it to the tree; `kind` one of the table's."""
    found = fused3d(mol, atoms, rings, rank)
    if found is not None:
        return 'shape', [p for p, _ in found[1]]
    found = _circle(mol, atoms, rings, rank)
    if found is None:
        found = _medium(mol, atoms, rings)
    if found is None:
        found = _peri(mol, atoms, rings)
    if found is not None:
        return 'circle', [found, {x: (-p[0], p[1]) for x, p in found.items()}]
    if len(rings) == 1:
        ring = rings[0]
        if len(ring) > TEMPLATE:
            return None
        r = L / (2 * sin(pi / len(ring)))
        plane = {x: (r * cos(2 * pi * i / len(ring)), r * sin(2 * pi * i / len(ring))) for i, x in enumerate(ring)}
        return 'ring', [plane, {x: (-p[0], p[1]) for x, p in plane.items()}]
    if min(map(len, rings)) > TEMPLATE:
        return None
    key, order = _key(mol, atoms)
    drawn = _drawing(key)
    if drawn is None:
        return None
    ref, pts = drawn
    best = None
    for mapping in islice(ref.as_query().get_mapping(key), MAPPINGS):
        mine = {order[mapping[i] - 1]: p for i, p in pts.items()}
        tie = tuple(rank[x] for _, x in sorted((i, order[mapping[i] - 1]) for i in pts))
        if best is None or tie < best[0]:
            best = (tie, mine)
    if best is None:
        return None
    plane = best[1]
    return ('system' if max(map(len, rings)) <= TEMPLATE else 'units'), [plane, {x: (-p[0], p[1]) for x, p in
                                                                           plane.items()}]


def _circle(mol, atoms, rings, rank):
    """The crown drawing of a block whose main cycle `cycles` calls a `Circle`, `{atom: (x, y)}`, or None: the
    cycle on `templates.crown`; each ring fused on one of its edges the regular polygon on that edge, away from the
    centre; None when an atom is left or a stated
    cis/trans unit would break, its substituents where their exits put them."""
    from ...core.chains import MACRO, Circle, cycles
    from .templates import crown
    if max(map(len, rings)) < MACRO:
        return None
    component = next(c for c in mol.connected_components if atoms <= set(c))
    found = cycles(mol.substructure(component))
    if not isinstance(found, Circle) or not set(found.ring) <= atoms:
        return None
    plane = {x: (p[0] * L, p[1] * L) for x, p in zip(found.ring, crown(len(found.ring), found.period))}
    for r in rings:
        if set(r) <= set(plane):
            continue
        shared = [x for x in r if x in plane]
        if len(shared) != 2 or shared[1] not in mol.neighbors_of(shared[0]):
            return None
        r = list(r)
        i = r.index(shared[0])
        path = r[i:] + r[:i]                            # from one shared atom round to the other
        if path[1] == shared[1]:
            path = [path[0]] + path[1:][::-1]
        (ax, ay), (bx, by) = plane[path[0]], plane[path[-1]]
        mx, my = (ax + bx) / 2, (ay + by) / 2
        k = len(r)
        h = L / (2 * tan(pi / k))
        nx, ny = -(by - ay) / L, (bx - ax) / L          # the edge's normal, turned outward
        if nx * mx + ny * my < 0:
            nx, ny = -nx, -ny
        cx, cy = mx + nx * h, my + ny * h
        rad = L / (2 * sin(pi / k))
        t0 = atan2(ay - cy, ax - cx)
        t1 = atan2(by - cy, bx - cx)
        step = min((2 * pi / k, -2 * pi / k),           # the way round that ends on the other shared atom
                   key=lambda d: abs((t0 + d * (k - 1) - t1 + pi) % (2 * pi) - pi))
        for j, x in enumerate(path[1:-1], 1):
            plane[x] = (cx + rad * cos(t0 + step * j), cy + rad * sin(t0 + step * j))
    if set(plane) != set(atoms):
        return None
    from .relax import _ct_broken, _frames
    probe = dict(plane)                                 # each substituent where its exit will put it: the bisector
    for x in found.ring:
        ring_nb = [y for y in mol.neighbors_of(x) if y in plane]
        ux = sum(plane[x][0] - plane[y][0] for y in ring_nb)
        uy = sum(plane[x][1] - plane[y][1] for y in ring_nb)
        d = hypot(ux, uy) or 1.
        for y in mol.neighbors_of(x):
            if y not in plane:
                probe[y] = (plane[x][0] + ux / d * L, plane[x][1] + uy / d * L)
    frames = _frames(mol, probe)
    return None if _ct_broken(mol, probe, frames, set(probe)) else plane


def _medium(mol, atoms, rings):
    """A block of a ring of `MEDIUM` atoms and smaller rings fused to it on one edge each, `{atom: (x, y)}`, or
    None: the medium ring a regular polygon, every other ring one outward on its edge, off the lattice; None when a
    ring shares no single edge, lands on the drawing, or a stated cis/trans unit would break."""
    from .relax import _ct_broken, _frames
    from .shapes import _polygons
    big = max(rings, key=len)
    if not MEDIUM[0] <= len(big) <= MEDIUM[1] or len(rings) < 2 or \
            sum(len(r) >= MEDIUM[0] for r in rings) > 1:
        return None
    n = len(big)
    r = L / (2 * sin(pi / n))
    plane = {x: (r * cos(2 * pi * i / n), r * sin(2 * pi * i / n)) for i, x in enumerate(big)}
    plane = _polygons(mol, plane, [q for q in rings if q is not big])
    if plane is None or set(plane) != set(atoms):
        return None
    probe = dict(plane)                                 # each substituent on its atom's outward bisector
    for x in atoms:
        ring_nb = [y for y in mol.neighbors_of(x) if y in plane]
        ux = sum(plane[x][0] - plane[y][0] for y in ring_nb)
        uy = sum(plane[x][1] - plane[y][1] for y in ring_nb)
        d = hypot(ux, uy) or 1.
        for y in mol.neighbors_of(x):
            if y not in plane:
                probe[y] = (plane[x][0] + ux / d * L, plane[x][1] + uy / d * L)
    return None if _ct_broken(mol, probe, _frames(mol, probe), set(probe)) else plane


def _peri(mol, atoms, rings):
    """A block whose largest ring, of 7 to `TEMPLATE` atoms, closes round an ortho-fused set of smaller rings at a
    peri atom, `{atom: (x, y)}`, or None: the smaller rings regular polygons each fused on one edge, the largest
    ring's free atoms an arc of equal bonds from one placed atom to the other, outward -- an indole's pyrrole stays
    regular and the azepine bridging it takes the strain.  None when the rest is no such arc or it lands on the
    drawing."""
    from .shapes import _polygons
    big = max(rings, key=len)
    small = [r for r in rings if r is not big]
    if not 7 <= len(big) <= TEMPLATE or not small or any(len(r) >= len(big) for r in small) or \
            not any(sum(x in r for r in rings) == 3 for x in big):
        return None
    first = small[0]
    n = len(first)
    rad = L / (2 * sin(pi / n))
    plane = {x: (rad * cos(2 * pi * i / n), rad * sin(2 * pi * i / n)) for i, x in enumerate(first)}
    plane = _polygons(mol, plane, small[1:])
    if plane is None:
        return None
    k = len(big)
    free = [i for i, x in enumerate(big) if x not in plane]
    if not free:
        return None
    start = next(i for i in range(k) if big[i] in plane and big[(i + 1) % k] not in plane)
    path = [big[start]]
    i = (start + 1) % k
    while big[i] not in plane:
        path.append(big[i])
        i = (i + 1) % k
    path.append(big[i])
    if len(path) - 2 != len(free):
        return None
    (ax, ay), (bx, by) = plane[path[0]], plane[path[-1]]
    c = hypot(bx - ax, by - ay)
    m = len(path) - 1                                   # bonds of the arc, each `L`
    if c >= m * L or c < 1e-9:
        return None
    lo, hi = 1e-6, 2 * pi / m - 1e-6                    # the arc's turn per bond: chord `c` of `m` chords `L`
    for _ in range(60):
        phi = (lo + hi) / 2
        if L * sin(m * phi / 2) / sin(phi / 2) > c:
            lo = phi
        else:
            hi = phi
    r = L / (2 * sin(phi / 2))
    mx, my = (ax + bx) / 2, (ay + by) / 2
    gx = sum(p[0] for p in plane.values()) / len(plane)
    gy = sum(p[1] for p in plane.values()) / len(plane)
    nx, ny = -(by - ay) / c, (bx - ax) / c
    if (mx - gx) * nx + (my - gy) * ny < 0:
        nx, ny = -nx, -ny                               # bulging away from what is drawn
    h = r * cos(m * phi / 2)
    cx, cy = mx - nx * h, my - ny * h
    ta = atan2(ay - cy, ax - cx)
    step = min((phi, -phi), key=lambda d: abs((ta + m * d - atan2(by - cy, bx - cx) + pi) % (2 * pi) - pi))
    added = {x: (cx + r * cos(ta + step * j), cy + r * sin(ta + step * j)) for j, x in enumerate(path[1:-1], 1)}
    if any(hypot(q[0] - p[0], q[1] - p[1]) < L / 2 for q in added.values() for p in plane.values()):
        return None
    plane.update(added)
    return plane if set(plane) == set(atoms) else None


def _key(mol, atoms):
    """`(key molecule, order)`: `atoms` and their bonds as an all-carbon, all-single-bond graph with one stub atom
    per bond leaving them; `order[i]` the atom of `mol` the key's atom `i + 1` stands for, None for a stub."""
    inside = sorted(atoms)
    out = MoleculeContainer()
    order = []
    with out.edit() as e:
        index = {}
        for x in inside:
            index[x] = e.add_atom('C')
            order.append(x)
        for x in inside:
            for y in mol.neighbors_of(x):
                if y in index:
                    if x < y:
                        e.add_bond(index[x], index[y], 1)
                else:
                    e.add_bond(index[x], e.add_atom('C'), 1)
                    order.append(None)
    return out, order


_DRAWINGS = {}           # key bytes -> drawing; cleared past `CACHE` keys


def _drawing(key):
    """`(key molecule, {key atom: (x, y)})` of the system the key stands for, or None; computed once per key."""
    k = key.canonical_bytes
    if k in _DRAWINGS:
        return _DRAWINGS[k]
    from .edits import _block_drawings                   # the block functions need the full layout package
    from .rings import ring_blocks
    key = _canonical_copy(key)                          # numbered by canonical rank: one drawing per key bytes
    rank = key.canonical_order()
    blocks = [b for b in ring_blocks(key) if len(b[0]) > 1]
    out = None
    CLEANS = REGULAR
    if blocks:
        atoms, rings = max(blocks, key=lambda b: len(b[0]))
        ortho = all(sum(x in r for r in rings) <= 2 for x in atoms)     # no peri or bridged atom
        first, best = None, None
        for pts, clean in _block_drawings(key, atoms, rings, rank):
            plane = {x: project(p) for x, p in pts.items()}
            if first is None:
                first = plane
            if clean and not ortho:
                best = (0., plane)
                break
            if clean:                                   # of the clean drawings the most regular, the first on a tie
                bad = _irregularity(plane, rings)
                if best is None or bad < best[0] - 1e-9:
                    best = (bad, plane)
                if bad < 1e-6:
                    break
                CLEANS -= 1
                if not CLEANS:
                    break
        if best is not None:
            first = best[1]
        if first is not None and not ortho and max(map(len, rings)) <= TEMPLATE:
            _untangle(key, first, atoms)                # a bridged small-ring system's own defects, nudged out once
        if first is not None and ortho:
            from .relax import regularize                   # the key has no stated geometry to keep: even its rings
            ring = key.substructure(atoms)               # the drawing holds the ring atoms alone
            kept = dict(first)
            regularize(ring, first, [tuple(r) for r in rings], {x: rank[x] for x in ring})
            if max(map(len, rings)) > TEMPLATE and _defects(key, first, atoms) > _defects(key, kept, atoms):
                first.clear()                           # a macrocycle evened onto itself: drawn as it was
                first.update(kept)
        if first is not None and (max(map(len, rings)) <= TEMPLATE or not _defects(key, first, atoms)):
            out = (key, first)
    if len(_DRAWINGS) >= CACHE:
        _DRAWINGS.clear()
    _DRAWINGS[k] = out
    return out


def _defects(key, plane, atoms):
    """The defects of the ring system `atoms` of `key` drawn as `plane`."""
    return sum(count_defects(plane, list(atoms), [(n, m) for n in atoms for m in key.neighbors_of(n) if m in atoms
                                                  and n < m]))


def _untangle(key, plane, atoms):
    """While `plane`'s drawing of `atoms` has atoms close or on a bond, up to `FIT` times: nudge in place the atom
    of a defect by the move of `NUDGES` leaving fewest of those, then fewest crossings, ties to the shorter move and
    the lower rank, until none leaves fewer."""
    rank = key.canonical_order()
    bonds = [(n, m) for n in atoms for m in key.neighbors_of(n) if m in atoms and n < m]
    order = sorted(atoms, key=rank.__getitem__)
    def worst(plane):                                   # atoms close or on a bond first, then crossings
        close, cross, on = count_defects(plane, order, bonds)
        return close + on, cross
    now = worst(plane)
    for _ in range(FIT):
        if not now[0]:
            return
        best = None
        from .relax import defect_atoms
        trouble = defect_atoms(plane, order, bonds)
        for x in (x for x in order if x in trouble):
            old = plane[x]
            for k, (dx, dy) in enumerate(NUDGES):
                plane[x] = (old[0] + dx, old[1] + dy)
                if any(hypot(plane[x][0] - plane[y][0], plane[x][1] - plane[y][1]) > LONGEST_BOND
                       for y in key.neighbors_of(x) if y in atoms):
                    continue
                d = worst(plane)
                if d < now and (best is None or (d, k) < best[0][:2]):
                    best = ((d, k, rank[x]), x, plane[x])
            plane[x] = old
        if best is None:
            return
        (now, _, _), x, q = best
        plane[x] = q


def _irregularity(plane, rings):
    """The worst ring of `plane`: its largest bond's departure from `L`, in `L`s, plus its largest angle's from the
    regular polygon's, in radians."""
    worst = 0.
    for r in rings:
        n = len(r)
        ideal = pi - 2 * pi / n
        bond = max(abs(hypot(plane[a][0] - plane[b][0], plane[a][1] - plane[b][1]) / L - 1)
                   for a, b in zip(r, r[1:] + r[:1]))
        angle = 0.
        for i in range(n):
            a, x, b = plane[r[i - 1]], plane[r[i]], plane[r[(i + 1) % n]]
            t = abs((atan2(a[1] - x[1], a[0] - x[0]) - atan2(b[1] - x[1], b[0] - x[0]) + pi) % (2 * pi) - pi)
            angle = max(angle, abs(t - ideal))
        worst = max(worst, bond + angle)
    return worst


def _canonical_copy(mol):
    """`mol` with its atoms numbered in canonical rank order, so equal canonical bytes give one numbering."""
    rank = mol.canonical_order()
    order = sorted(mol, key=rank.__getitem__)
    out = MoleculeContainer()
    with out.edit() as e:
        index = {x: e.add_atom(mol.atom(x).atomic_symbol) for x in order}
        for n, m in sorted(((n, m) for n in mol for m in mol.neighbors_of(n) if rank[n] < rank[m]),
                           key=lambda b: (rank[b[0]], rank[b[1]])):
            e.add_bond(index[n], index[m], mol.bond(n, m).order)
    return out


def fit(initial):
    """The plane of `initial`'s component after the fit: `{atom: (x, y)}`. The fit runs by flips and bends alone and
    with stretches allowed, neither leaving an atom off even, an acyclic atom's gaps at `360 / n` or a ring atom's lone
    substituent on its bisector, the better by `_judge` winning, the flips on a tie; when the winner leaves a defect or
    a broken cis/trans, the fit allowing atoms off even is run too, and the one leaving less wins; when it leaves a bond
    `_cramped`, that fit is run too and wins when `_judge` weighs it better."""
    flips, left = _fit(initial, True, False)
    judged = _judge(initial, flips, left)
    if judged == ((0, 0), 0, 0, 0, 0, False, 0.):       # nothing a stretch could better
        return flips
    plane, rest = _fit(initial, True)
    stretched = _judge(initial, plane, rest)
    if judged <= stretched:
        plane, rest = flips, left
    else:
        judged = stretched
    if rest > (0, 0):
        other, more = _fit(initial, False)
        if more < rest:
            return other
    elif judged[2]:                                     # a bond cramped onto its ring: the fit off even may free it
        other, more = _fit(initial, False)
        if _judge(initial, other, more) < judged:
            return other
    return plane


def _judge(initial, plane, left):
    """`(defects and broken cis/trans, atoms inside rings, atoms `_cramped`, atoms off even, ring atoms `_squeezed`, an
    unbonded pair nearer than `SPACE`, stretch in `L`s)` of a fitted `plane`: how `fit` weighs a drawing by flips
    against one by stretches."""
    scope = initial.scope
    rings = [r for r in scope.rings if r[0] in plane]
    grid = Grid()
    for x in scope.atoms:
        grid.add(x, plane[x])
    close = any(y not in scope.within2[x] and hypot(plane[x][0] - plane[y][0], plane[x][1] - plane[y][1]) < SPACE
                for x in scope.atoms for y in grid.atoms(plane[x], SPACE) if y != x)
    stretch = sum(max(hypot(plane[p][0] - plane[x][0], plane[p][1] - plane[x][1]) - L, 0.)
                  for p, x, _ in initial.edges)
    every = set(scope.atoms)
    return (left, len(_inside_rings(plane, rings)), _cramped(scope, plane), _uneven(scope, plane, every),
            _squeezed(scope, plane, every), close, round(stretch / L, 2))


def _fit(initial, even, stretch=True):
    """`(plane, (defects, broken cis/trans))` after the fit; `even` forbids a move leaving an atom off even, no
    `stretch` any stretch of the repair."""
    env = initial.env()
    scope, mol = initial.scope, initial.mol
    orient, scope.orient = scope.orient, False          # the fit reads defects; the pose orients afterwards
    rings = [r for r in scope.rings if r[0] in env.plane]   # this component's own
    _repair(env, scope, mol, rings, even, stretch=stretch, gem=stretch)
    _unstick(env, scope, mol, rings, even, stretch)
    if even:
        _settle(env, initial, rings)
    _nudge(env, scope, mol)
    if even and env.defects <= 0 and not env.ct:
        _space(env, scope, mol)
        _unload(env, scope, mol, rings)
    scope.orient = orient
    return dict(env.plane), (max(env.defects, 0), env.ct)


def _repair(env, scope, mol, rings, even, held=frozenset(), banned=frozenset(), budget=FIT, stretch=True, gem=False):
    """The fit's moves, up to `budget`, while a defect, a broken stated cis/trans or an atom inside a ring is left; no
    move shifts an atom of `held` or turns a bond of `banned`, nor, unless `stretch`, stretches one, nor, with `gem`,
    turns a ring atom's two substituents off even: the fit that stretches lengthens their bonds instead."""
    from .relax import _ct_broken, defect_atoms
    for _ in range(budget):
        inside = _inside_rings(env.plane, rings)
        if env.defects <= 0 and not env.ct and not inside:
            break
        trouble = defect_atoms(env.plane, scope.atoms, scope.bonds) | inside
        if env.ct:                                      # a broken stated cis/trans: its unit's atoms are trouble too
            trouble |= {x for f in scope.frames if _ct_broken(mol, env.plane, [f], set(f[2])) for x in f[2]}
        best = _best(env, scope, mol, rings, even, held, trouble,
                     (m for m in env._moves() if (m.pivot, m.other) not in banned and (stretch or m.kind != 'stretch')),
                     bool(inside), gem)
        if best is None:
            break
        _, moved, d, ct = best
        env.plane.update(moved)
        env.terms = tuple(a + b for a, b in zip(env.terms, d))
        env.ct += ct
        env._legal = None


def _bridges(initial):
    """The tree bonds on no ring joining two fragments of `MASSIVE` or more atoms each, and the tree bonds beside them:
    where a flip or stretch moves a whole fragment clear of another."""
    scope, n = initial.scope, len(initial.scope.atoms)
    core = {(p, x) for p, x, part in initial.edges if frozenset((p, x)) not in scope.ring_bonds and
            min(len(part), n - len(part)) >= MASSIVE}
    ends = {y for e in core for y in e}
    return {(p, x) for p, x, _ in initial.edges if (p, x) in core or p in ends or x in ends}


def _unstick(env, scope, mol, rings, even, stretch=True):
    """While the drawing leaves a defect, an atom off even or a long bond, up to `FIT` times: of the flips, anchored
    rotations and, when `stretch`, lengthening stretches up to `REACH` steps of a tree bond, none putting an atom off
    even when `even`, the `WIDE` at `_bridges` adding fewest defects, and while a defect is left the `UNSTICK` adding
    fewest of the flips and rotations whose subtree holds some atoms of the trouble but not all, each taken and up to
    `AFTER` repair moves run after it with that bond fixed; the one best by `_scene`, kept when better than now, at
    most `POLISH` times once no defect is left: a fragment moved clear only once the rest is repaired round it."""
    from .edits import masked
    from .primitive import LONGEST
    from .relax import defect_atoms
    from .score import delta_terms
    every = set(scope.atoms)
    bridges = _bridges(env.initial)
    for step in range(FIT):
        now = _scene(env, scope, every, rings)
        if now[:5] == (0, 0, 0, 0, 0.) or not now[0] and step >= POLISH:
            return
        trouble = defect_atoms(env.plane, scope.atoms, scope.bonds) if env.defects > 0 else set()
        cheap, wide = [], []
        for m in env._moves():
            at = (m.pivot, m.other) in bridges
            if m.kind == 'stretch':
                if not stretch or not at or m.arg <= 0:
                    continue
                kicks = [m._replace(arg=k * m.arg) for k in range(1, REACH + 1)]
            elif m.kind in ('flip', 'rotate') and (at or m.part & trouble and not trouble <= m.part):
                kicks = [m]
            else:
                continue
            for k in kicks:
                moved = env._moved(k)
                cut = {k.pivot, k.other}
                if masked(scope, env.plane, moved, LONGEST) or even and \
                        _uneven(scope, {**env.plane, **moved}, cut, stretch) > _uneven(scope, env.plane, cut, stretch):
                    continue
                d = delta_terms(scope, env.plane, moved)
                (wide if at else cheap).append(((d[0], KIND[k.kind], len(k.part), env.tie(k)), k, moved, d))
        tries = sorted(cheap, key=lambda o: o[0])[:UNSTICK] + sorted(wide, key=lambda o: o[0])[:WIDE]
        best = None
        saved = (dict(env.plane), env.terms, env.ct)
        for _, m, moved, d in tries:
            env.plane.update(moved)
            env.terms = tuple(a + b for a, b in zip(env.terms, d))
            env._legal = None
            _repair(env, scope, mol, rings, even, banned=frozenset({(m.pivot, m.other)}), budget=AFTER,
                    stretch=stretch, gem=stretch)
            got = _scene(env, scope, every, rings)
            if got < now and (best is None or got < best[0]):
                best = (got, dict(env.plane), env.terms, env.ct)
            env.plane, env.terms, env.ct = dict(saved[0]), saved[1], saved[2]
            env._legal = None
        if best is None:
            return
        _, env.plane, env.terms, env.ct = best
        env._legal = None


def _scene(env, scope, every, rings):
    """`(defects, broken cis/trans, atoms inside rings, atoms off even, stretch, -nearest unbonded pair)` of `env`'s
    drawing, distances in `L`s to 2 decimals: the order `_unstick` compares results by."""
    plane = env.plane
    grid = Grid()
    for x in scope.atoms:
        grid.add(x, plane[x])
    near = min((hypot(plane[x][0] - plane[y][0], plane[x][1] - plane[y][1]) for x in scope.atoms
                for y in grid.atoms(plane[x], L) if x < y and y not in scope.within2[x]), default=L)
    stretch = sum(max(hypot(plane[p][0] - plane[x][0], plane[p][1] - plane[x][1]) - L, 0.)
                  for p, x, _ in env.initial.edges)
    return (max(env.defects, 0), env.ct, len(_inside_rings(plane, rings)), _uneven(scope, plane, every),
            round(stretch / L, 2), -round(min(near, L) / L, 2))


def _best(env, scope, mol, rings, even, held, trouble, moves, inside=True, gem=False):
    """`(key, moved, delta terms, delta cis/trans)` of the fit's best of `moves` on `trouble`, or None."""
    from .edits import masked
    from .primitive import LONGEST
    from .relax import _ct_broken, _inr
    from .score import delta_terms
    best = None
    for m in moves:
        if m.kind == 'turn' or not m.part & trouble or trouble <= m.part:
            continue
        if not m.part & held:
            moved = env._moved(m)
            if m.kind == 'stretch':
                moved = _even_line(env, scope, mol, m, moved)
        elif held < m.part:                             # the same picture: the rest moved back
            moved = _inverse(env, m)
        else:
            continue
        if m.kind == 'stretch':
            (px, py), (qx, qy) = moved.get(m.pivot, env.plane[m.pivot]), moved.get(m.other, env.plane[m.other])
            if ((px - qx) ** 2 + (py - qy) ** 2) ** .5 > STRETCHED + 1e-9:
                continue
        if masked(scope, env.plane, moved, LONGEST):
            continue
        part = set(moved)
        d = delta_terms(scope, env.plane, moved)
        ct = _ct_broken(mol, {**env.plane, **moved}, scope.frames, set(moved)) - \
            _ct_broken(mol, env.plane, scope.frames, set(moved))
        if ct > 0 or (not ct and d[0] > 0) or (ct, d[0]) == (0, 0) and not inside:
            continue                                    # nothing inside a ring: a move must repair to count
        inr = _inr({**env.plane, **moved}, rings, part) - _inr(env.plane, rings, part)
        if inr > 0 or (ct, d[0], inr) >= (0, 0, 0):
            continue                                    # never an atom into a ring
        if _pinched(scope, {**env.plane, **moved}, part, even) > _pinched(scope, env.plane, part, even):
            continue                                    # never two bonds of an atom squeezed together
        bent = _bent_exits(scope, {**env.plane, **moved}, part) - _bent_exits(scope, env.plane, part)
        if bent > 0 and (ct, d[0]) >= (0, 0):
            continue                                    # a substituent bent onto its ring only to clear a defect
        if _cramped(scope, {**env.plane, **moved}, part) > _cramped(scope, env.plane, part):
            continue                                    # never a substituent folded onto its ring bond
        squeezed = _squeezed(scope, {**env.plane, **moved}, part) - _squeezed(scope, env.plane, part)
        if squeezed > 0 and (ct, d[0]) >= (0, 0):
            continue
        shaped = env.initial.shaped
        if _inward(scope, {**env.plane, **moved}, part, shaped) > _inward(scope, env.plane, part, shaped):
            continue                                    # never a substituent across its own ring
        uneven = _uneven(scope, {**env.plane, **moved}, part, gem) - _uneven(scope, env.plane, part, gem) if even else 0
        if uneven > 0:
            continue                                    # the even fit never bends: it stretches or flips instead
        sp3 = m.kind == 'stretch' and mol.hybridization_of(m.pivot) == 1 and mol.hybridization_of(m.other) == 1
        end = m.kind == 'stretch' and len(scope.adj[m.other]) == 1         # a terminal bond keeps its length
        line = m.kind == 'stretch' and not end and (m.pivot in scope.linear or m.other in scope.linear) and \
            not mol.atom(m.pivot).is_metal and not mol.atom(m.other).is_metal
        key = (ct, not line, d[0], bent, squeezed, uneven, inr, end, KIND[m.kind], len(moved), not sp3, env.tie(m))
        if best is None or key < best[0]:
            best = (key, moved, d, ct)
    return best


def _even_line(env, scope, mol, m, moved):
    """`moved` of the stretch `m` with its linear end, of two non-terminal neighbours, no metal among the three, moved
    half as far: the stretch shared by both bonds of the line, Si-Si-Si kept symmetric."""
    for x in (m.pivot, m.other):
        if x in scope.linear and not mol.atom(x).is_metal and all(len(scope.adj[y]) > 1 and not mol.atom(y).is_metal
                                                                  for y in scope.adj[x]):
            (px, py), (qx, qy) = env.plane[m.pivot], env.plane[m.other]
            d = hypot(qx - px, qy - py) or 1.
            dx, dy = (qx - px) / d * m.arg / 2, (qy - py) / d * m.arg / 2
            return {**moved, x: (env.plane[x][0] + dx, env.plane[x][1] + dy)}
    return moved


def _inverse(env, m):
    """`{atom: point}` of the atoms of `env`'s plane outside `m.part` moved by the inverse of `m`: the drawing `m`
    makes, up to a rigid motion, with `m.part` left in place."""
    from .primitive import _axis
    plane = env.plane
    (px, py), (qx, qy) = plane[m.pivot], plane[m.other]
    rest = [x for x in plane if x not in m.part]
    if m.kind == 'flip':
        dx, dy = _axis(px, py, qx, qy, m.arg)
        d = dx * dx + dy * dy or 1.
        out = {}
        for x in rest:
            ux, uy = plane[x][0] - px, plane[x][1] - py
            k = 2 * (ux * dx + uy * dy) / d
            out[x] = (px + k * dx - ux, py + k * dy - uy)
        return out
    if m.kind == 'rotate':
        t = -env._turn(m)
        c, s = cos(t), sin(t)
        return {x: (px + c * (plane[x][0] - px) - s * (plane[x][1] - py),
                    py + s * (plane[x][0] - px) + c * (plane[x][1] - py)) for x in rest}
    d = hypot(qx - px, qy - py) or 1.                   # stretch
    dx, dy = (px - qx) / d * m.arg, (py - qy) / d * m.arg
    return {x: (plane[x][0] + dx, plane[x][1] + dy) for x in rest}


def _settle(env, initial, rings):
    """Each ring substituent of up to `SMALL_SETTLE` atoms, not of a shape's projection, inside a ring, bent onto its
    ring or off its bisector, in rank order, up to `SETTLE`: the flip or anchored rotation of its tree bond putting it
    on the bisector outside every ring, then the repair with that bond never turned again; kept when that leaves no more
    defects or broken cis/trans and fewer atoms inside rings, bent or off even, the first such move of the fewest
    defects added."""
    from .relax import _ct_broken, _inr
    from .score import delta_terms
    scope, mol = initial.scope, initial.mol
    every = set(scope.atoms)

    def worse(plane):
        return len(_inside_rings(plane, rings)) + _bent_exits(scope, plane, every) + _uneven(scope, plane, every)
    targets = [(p, x, part) for p, x, part in initial.edges if len(part) <= SMALL_SETTLE and len(scope.adj[p]) == 3
               and p not in initial.projected
               and sum(frozenset((p, y)) in scope.ring_bonds for y in scope.adj[p]) == 2]
    tries = 0
    for p, x, part in sorted(targets, key=lambda e: scope.rank[e[1]]):
        if tries >= SETTLE:
            break
        if not (_off(scope, env.plane, p) or _bent_exits(scope, env.plane, part) or
                part & _inside_rings(env.plane, rings)):
            continue
        tries += 1
        now = (max(env.defects, 0), env.ct)
        bad = worse(env.plane)
        options = []
        for m in env._moves():
            if m.pivot != p or m.other != x or m.kind not in ('flip', 'rotate'):
                continue
            moved = env._moved(m)
            new = {**env.plane, **moved}
            if _off(scope, new, p) or _bent_exits(scope, new, part) or \
                    _inr(new, rings, set(moved)) > _inr(env.plane, rings, set(moved)):
                continue
            d = delta_terms(scope, env.plane, moved)
            options.append(((d[0], env.tie(m)), moved, d))
        for _, moved, d in sorted(options, key=lambda o: o[0]):
            saved = (dict(env.plane), env.terms, env.ct)
            env.ct += _ct_broken(mol, {**env.plane, **moved}, scope.frames, set(moved)) - \
                _ct_broken(mol, env.plane, scope.frames, set(moved))
            env.plane.update(moved)
            env.terms = tuple(a + b for a, b in zip(env.terms, d))
            env._legal = None
            _repair(env, scope, mol, rings, True, banned=frozenset({(p, x)}))
            if (max(env.defects, 0), env.ct) <= now and worse(env.plane) < bad:
                break
            env.plane, env.terms, env.ct = saved
            env._legal = None


def _overload(scope, plane, moved, grid):
    """`(tight, loose)` local overload between the atoms of `moved` and the rest, `grid` holding every atom where
    `plane` has it: over each such pair, neither bonded nor geminal, the squared shortfall below `CROWDED` and below
    `LOOSE`, in `L`s."""
    tight = loose = 0.
    for x, (px, py) in moved.items():
        for y in grid.atoms((px, py), LOOSE):
            if y in moved or y in scope.within2[x]:
                continue
            d = hypot(px - plane[y][0], py - plane[y][1])
            if d < LOOSE:
                loose += ((LOOSE - d) / L) ** 2
                if d < CROWDED:
                    tight += ((CROWDED - d) / L) ** 2
    return tight, loose


def _unload(env, scope, mol, rings):
    """Up to `FIT` times: the flip of a tree bond lowering most the tight `_overload` of its subtree against the rest,
    by `UNLOAD` or more, then the loose, adding no defect, broken cis/trans, atom off even, substituent
    `in_ring_angle` or atom inside a ring: a crowded corner opened by a flip where nothing is wrong."""
    from .edits import masked
    from .primitive import LONGEST
    from .relax import _ct_broken, _inr
    from .score import delta_terms
    shaped = env.initial.shaped
    for _ in range(FIT):
        plane = env.plane
        grid = Grid()
        for x in scope.atoms:
            grid.add(x, plane[x])
        if not any(y not in scope.within2[x] and hypot(plane[x][0] - plane[y][0], plane[x][1] - plane[y][1]) < CROWDED
                   for x in scope.atoms for y in grid.atoms(plane[x], CROWDED) if y != x):
            return                                      # no corner overloaded
        best = None
        for m in env._moves():
            if m.kind != 'flip':
                continue
            moved = env._moved(m)
            (t1, l1), (t0, l0) = _overload(scope, plane, moved, grid), \
                _overload(scope, plane, {x: plane[x] for x in moved}, grid)
            key = (round(t1 - t0, 6), round(l1 - l0, 6), len(m.part), env.tie(m))
            if t1 - t0 > -UNLOAD or best is not None and key >= best[0]:
                continue
            new = {**plane, **moved}
            part, cut = set(moved), {m.pivot, m.other}
            if masked(scope, plane, moved, LONGEST) or \
                    _uneven(scope, new, cut) > _uneven(scope, plane, cut) or \
                    _inward(scope, new, cut, shaped) > _inward(scope, plane, cut, shaped) or \
                    delta_terms(scope, plane, moved)[0] > 0 or \
                    _ct_broken(mol, new, scope.frames, part) > _ct_broken(mol, plane, scope.frames, part) or \
                    _inr(new, rings, part) > _inr(plane, rings, part):
                continue
            best = (key, moved)
        if best is None:
            return
        env.plane.update(best[1])
        env._legal = None


def _space(env, scope, mol):
    """While two atoms, neither bonded nor geminal, are nearer than `SPACE`, up to `FIT` times: the flip or stretch
    parting most such pairs that adds no defect, breaks no stated cis/trans and leaves no atom off even."""
    from .edits import masked
    from .primitive import LONGEST
    from .relax import _ct_broken
    from .score import delta_terms
    def crowded(plane):
        g = Grid()
        for x in scope.atoms:
            g.add(x, plane[x])
        return {(x, y) for x in scope.atoms for y in g.atoms(plane[x], SPACE) if x < y and y not in scope.within2[x]
                and hypot(plane[x][0] - plane[y][0], plane[x][1] - plane[y][1]) < SPACE}
    for _ in range(FIT):
        pairs = crowded(env.plane)
        if not pairs:
            return
        trouble = {x for p in pairs for x in p}
        best = None
        for m in env._moves():
            if m.kind not in ('flip', 'stretch') or not m.part & trouble or trouble <= m.part:
                continue
            moved = env._moved(m)
            (px, py), (qx, qy) = env.plane[m.pivot], moved[m.other]
            if m.kind == 'stretch' and hypot(px - qx, py - qy) > STRETCHED + 1e-9 or \
                    masked(scope, env.plane, moved, LONGEST):
                continue
            new = {**env.plane, **moved}
            part = set(moved)
            if delta_terms(scope, env.plane, moved)[0] > 0 or \
                    _uneven(scope, new, part) > _uneven(scope, env.plane, part) or \
                    _ct_broken(mol, new, scope.frames, part) > _ct_broken(mol, env.plane, scope.frames, part):
                continue
            n = len(crowded(new))
            key = (n, KIND[m.kind], len(m.part), env.tie(m))
            if n < len(pairs) and (best is None or key < best[0]):
                best = (key, moved)
        if best is None:
            return
        env.plane.update(best[1])
        env._legal = None


def _nudge(env, scope, mol):
    """While a defect is left, up to `FIT` times: move the atom of the trouble, with its terminal neighbours, by
    whichever of `NUDGES`, or an end atom round its neighbour by whichever of 24 directions, removes most defects and
    breaks no stated cis/trans: a crowded cage atom pulled clear, but never an atom of a cage whose depth places a
    stereocentre's substituent, nor that substituent."""
    from .edits import masked
    from .primitive import LONGEST
    from .relax import _ct_broken, defect_atoms
    from .score import delta_terms
    for _ in range(FIT):
        if env.defects <= 0:
            return
        best = None
        faced = {v for bond in env.initial.faced for v in bond}
        faced |= {x for b in env.initial.rigid if b & faced for x in b}
        for x in sorted(defect_atoms(env.plane, scope.atoms, scope.bonds), key=scope.rank.__getitem__):
            if x in faced:                              # a cage its depth reads, and its faced substituents, stay
                continue
            part = [x] + [y for y in scope.adj[x] if len(scope.adj[y]) == 1]
            tries = [(hypot(dx, dy), {y: (env.plane[y][0] + dx, env.plane[y][1] + dy) for y in part})
                     for dx, dy in NUDGES]
            if len(scope.adj[x]) == 1:                  # an end atom also swings round its neighbour
                (ax, ay), = [env.plane[y] for y in scope.adj[x]]
                tries += [(L, {x: (ax + L * cos(k * pi / 12), ay + L * sin(k * pi / 12))}) for k in range(24)]
            for k, (far, moved) in enumerate(tries):
                if masked(scope, env.plane, moved, LONGEST):
                    continue
                d = delta_terms(scope, env.plane, moved)
                if d[0] >= 0 or _ct_broken(mol, {**env.plane, **moved}, scope.frames, set(moved)) > \
                        _ct_broken(mol, env.plane, scope.frames, set(moved)):
                    continue
                key = (d[0], far, scope.rank[x], k)
                if best is None or key < best[0]:
                    best = (key, moved, d)
        if best is None:
            return
        _, moved, d = best
        env.plane.update(moved)
        env.terms = tuple(a + b for a, b in zip(env.terms, d))
        env._legal = None


def _bent_exits(scope, plane, part):
    """Bonds from a ring atom of two ring bonds and one substituent to that substituent, an end in `part`, nearer
    than `EXIT` to one of the atom's ring bonds."""
    n = 0
    for x, y in scope.bonds:
        if frozenset((x, y)) in scope.ring_bonds or x not in part and y not in part:
            continue
        for r, s in ((x, y), (y, x)):
            ring = [z for z in scope.adj[r] if frozenset((r, z)) in scope.ring_bonds]
            if len(ring) != 2 or len(scope.adj[r]) != 3:       # a ring atom of one substituent
                continue
            u = atan2(plane[s][1] - plane[r][1], plane[s][0] - plane[r][0])
            n += sum(abs((u - atan2(plane[z][1] - plane[r][1], plane[z][0] - plane[r][0]) + pi) % (2 * pi) - pi) <
                     EXIT for z in ring)
    return n


def _squeezed(scope, plane, part):
    """Ring atoms of two ring bonds and two substituents, in `part` or bonded to it, whose substituents are nearer
    than a third of the exterior angle less `EVEN`: off the even slots, which split it in three."""
    n = 0
    for x in {z for y in part for z in (y, *scope.adj[y])}:
        ring = [y for y in scope.adj[x] if frozenset((x, y)) in scope.ring_bonds]
        exo = [y for y in scope.adj[x] if frozenset((x, y)) not in scope.ring_bonds]
        if len(ring) != 2 or len(exo) != 2:
            continue
        px, py = plane[x]
        u, v = (atan2(plane[y][1] - py, plane[y][0] - px) for y in ring)
        a, b = (atan2(plane[y][1] - py, plane[y][0] - px) for y in exo)
        inner = abs((u - v + pi) % (2 * pi) - pi)
        n += abs((a - b + pi) % (2 * pi) - pi) < (2 * pi - inner) / 3 - EVEN
    return n


def _inward(scope, plane, part, shaped):
    """Bonds from a ring atom, not in `shaped`, to a substituent, an end in `part`, `in_ring_angle` of the atom."""
    n = 0
    for x, y in scope.bonds:
        if frozenset((x, y)) in scope.ring_bonds or x not in part and y not in part:
            continue
        for r, s in ((x, y), (y, x)):
            if r not in shaped and in_ring_angle(scope, plane, r, atan2(plane[s][1] - plane[r][1],
                                                                        plane[s][0] - plane[r][0])):
                n += 1
    return n


def in_ring_angle(scope, plane, a, t):
    """True when the direction `t` from `a` lies inside the angle, below 180 degrees, of `a`'s two bonds in a ring
    of up to `TEMPLATE` atoms."""
    px, py = plane[a]
    for r in scope.rings:
        if a not in r or len(r) > TEMPLATE:
            continue
        i = r.index(a)
        u, v = (atan2(plane[z][1] - py, plane[z][0] - px) for z in (r[i - 1], r[(i + 1) % len(r)]))
        inner = (v - u) % (2 * pi)
        if inner > pi:                                  # the smaller angle, from u anticlockwise
            u, inner = v, 2 * pi - inner
        if 1e-6 < (t - u) % (2 * pi) < inner - 1e-6:
            return True
    return False


def _pinched(scope, plane, part, ring=True):
    """Atoms, in `part` or bonded to it, two of whose bonds are nearer than `PINCH`: of an acyclic atom of three or
    more bonds any two, of a ring atom, when `ring`, two substituents."""
    n = 0
    for x in {z for y in part for z in (y, *scope.adj[y])}:
        if len(scope.adj[x]) < 3:
            continue
        exo = [y for y in scope.adj[x] if frozenset((x, y)) not in scope.ring_bonds]
        if len(exo) < 2 or not ring and len(exo) < len(scope.adj[x]):
            continue
        px, py = plane[x]
        t = sorted(atan2(plane[y][1] - py, plane[y][0] - px) for y in exo)
        if len(exo) < len(scope.adj[x]):                # a ring atom: the gaps between its substituents alone
            gap = min(v - u for u, v in zip(t, t[1:]))
            gap = min(gap, 2 * pi - t[-1] + t[0]) if len(exo) > 2 else min(gap, 2 * pi - gap)
        else:
            gap = min(min(v - u for u, v in zip(t, t[1:])), 2 * pi - t[-1] + t[0])
        n += gap < PINCH
    return n


def _cramped(scope, plane, part=None):
    """Ring atoms, in `part` or bonded to it, or of the whole drawing, with a substituent's bond folded onto a ring
    bond, which neither `_pinched` nor `_squeezed`, both between substituents, counts: nearer than `PINCH` at an atom
    of two ring bonds and two substituents, whose exterior angle has a slot for each, nearer than `CRAMP` elsewhere."""
    n = 0
    for x in scope.atoms if part is None else {z for y in part for z in (y, *scope.adj[y])}:
        if len(scope.adj[x]) < 3:
            continue
        px, py = plane[x]
        t = sorted((atan2(plane[y][1] - py, plane[y][0] - px), frozenset((x, y)) in scope.ring_bonds)
                   for y in scope.adj[x])
        least = PINCH if len(t) == 4 and sum(r for _, r in t) == 2 else CRAMP
        n += any((v - u) % (2 * pi) < least and r != q for (u, r), (v, q) in zip(t, t[1:] + t[:1]))
    return n


def _uneven(scope, plane, part, gem=False):
    """Atoms, in `part` or bonded to it, `_off` even."""
    return sum(_off(scope, plane, x, gem) for x in {z for y in part for z in (y, *scope.adj[y])})


def _off(scope, plane, x, gem=False):
    """True when `x` is off even: an acyclic atom of two bonds, not linear, `EVEN` or more off 120 degrees; of n >= 3
    bonds a gap of which is `EVEN` or more below `360 / n`; a ring atom of one substituent that substituent `EVEN` or
    more off its bisector, and with `gem` of two a gap of its exterior angle twice `EVEN` off a third of it; never a
    metal, nor an atom of a shape's projection."""
    if x in scope.projected:                            # a 3D view's angles are the view's, not off even
        return False
    k = len(scope.adj[x])
    ring = sum(frozenset((x, y)) in scope.ring_bonds for y in scope.adj[x])
    if k < 2 or ring and (ring, k) not in ((2, 3), (2, 4) if gem else ((2, 3),)) or scope.mol.atom(x).is_metal:
        return False                                    # a metal keeps its own geometry
    px, py = plane[x]
    t = sorted(atan2(plane[y][1] - py, plane[y][0] - px) for y in scope.adj[x])
    gaps = [v - u for u, v in zip(t, t[1:])] + [2 * pi - t[-1] + t[0]]
    if ring and k == 4:                                 # two substituents: the exterior angle in three equal gaps
        bonds = sorted((atan2(plane[y][1] - py, plane[y][0] - px), frozenset((x, y)) in scope.ring_bonds)
                       for y in scope.adj[x])
        split = [((v - u) % (2 * pi), r and q) for (u, r), (v, q) in zip(bonds, bonds[1:] + bonds[:1])]
        inner = [g for g, both in split if both]
        if len(inner) != 1:                             # a substituent between the ring bonds: inside its ring
            return False
        third = (2 * pi - inner[0]) / 3
        return any(abs(g - third) > 2 * EVEN for g, both in split if not both)
    if ring:                                            # the two gaps beside the substituent equal
        y = next(y for y in scope.adj[x] if frozenset((x, y)) not in scope.ring_bonds)
        u = atan2(plane[y][1] - py, plane[y][0] - px)
        sides = [abs((u - atan2(plane[z][1] - py, plane[z][0] - px) + pi) % (2 * pi) - pi) for z in scope.adj[x]
                 if z != y]
        return abs(sides[0] - sides[1]) > 2 * EVEN
    if k == 2:                                          # a chain atom at 120, linear exempt
        return x not in scope.linear and abs(min(gaps) - 2 * pi / 3) > EVEN
    return min(gaps) < 2 * pi / k - EVEN


def _inside_rings(plane, rings):
    """The atoms inside a ring polygon they are not in, with that ring's atoms."""
    from .relax import _inside
    out = set()
    for r in rings:
        poly = [plane[x] for x in r]
        ring = set(r)
        found = {x for x in plane if x not in ring and _inside(plane[x], poly)}
        if found:
            out |= found | ring
    return out


def cascade_layout(mol, log=None):
    """`{atom: (x, y)}` for every atom of `mol`, mean bond near `L`."""
    from ...core import AutomorphismBudgetExceeded
    from .pose import GAP, SPAN, SPIRO, horizontal, lean, level, straight, upright
    from .primitive import Initial
    from .rings import ring_blocks
    if not len(mol):
        return {}
    rank = mol.canonical_order()
    blocks = ring_blocks(mol)
    try:
        orbit = mol.automorphism_orbits()
    except AutomorphismBudgetExceeded:
        orbit = None
    plane = {}
    x0 = 0.
    for c in sorted((set(c) for c in mol.connected_components), key=lambda c: min(rank[x] for x in c)):
        large = sum(mol.atom(x).atomic_symbol != 'H' for x in c) > LARGE
        if large:
            pts = _large(mol, c)
        else:
            initial = Initial(mol, c, rank)
            pts = fit(initial)
        if not large and not initial.upright:
            rings = [r for b, rs in blocks if b <= c for r in rs]
            level(pts)
            spiro = any(sum(x in r for r in rings) > 1 and sum(1 for y in mol.neighbors_of(x) if
                                                               any(x in r and y in r for r in rings)) == 4 for x in c)
            horizontal(pts, rings, SPIRO if spiro else SPAN)
            lean(pts, [r for b, rs in blocks if b <= c and all(sum(x in q for q in rs) <= 2 for x in b) for r in rs])
            if not rings:
                straight(mol, pts, rank)
            upright(mol, c, pts, orbit)
        xs = [p[0] for p in pts.values()]
        ys = [p[1] for p in pts.values()]
        dx, dy = x0 - min(xs), -(min(ys) + max(ys)) / 2
        for a, (x, y) in pts.items():
            plane[a] = (x + dx, y + dy)
        x0 = max(xs) + dx + GAP
    return plane


def _large(mol, atoms):
    """`{atom: (x, y)}` of a component over `LARGE` heavy atoms, drawn by smiles-drawer and scaled to a mean bond of
    `L`: the fit's repair scores every move over the whole drawing at each step, so its cost grows with the square of
    the size and a 500-atom oligonucleotide takes many minutes."""
    from .molecule import _engine_layout
    sub = mol.substructure(atoms)
    plane = _engine_layout(sub, 'smilesdrawer')
    bonds = [hypot(plane[b.n][0] - plane[b.m][0], plane[b.n][1] - plane[b.m][1]) for b in sub.bonds()]
    k = L * len(bonds) / (sum(bonds) or 1.)
    return {x: (k * plane[x][0], k * plane[x][1]) for x in atoms}


def covered(mol):
    """True when every ring system of `mol` is in the cascade's frozen domain, a crown among them."""
    from .rings import ring_blocks
    for atoms, rings in ring_blocks(mol):
        if max(map(len, rings)) > TEMPLATE and _circle(mol, atoms, rings, mol.canonical_order()) is not None:
            continue                                    # a crown
        if len(rings) == 1:
            if len(rings[0]) > TEMPLATE:
                return False
        elif shape(mol, atoms) is None and (max(map(len, rings)) > TEMPLATE or any(
                len(set(a) & set(b)) not in (0, 2) for i, a in enumerate(rings) for b in rings[:i])):
            return False
    return True


def drawing_hash(mol, plane):
    """A digest of `plane` by canonical rank, `DIGITS` decimals of `L`, centred: equal for the same drawing in any
    atom order."""
    rank = mol.canonical_order()
    cx = sum(p[0] for p in plane.values()) / len(plane)
    cy = sum(p[1] for p in plane.values()) / len(plane)
    pts = sorted((rank[x], round((p[0] - cx) / L, DIGITS) + 0., round((p[1] - cy) / L, DIGITS) + 0.)
                 for x, p in plane.items())
    return sha1(repr(pts).encode()).hexdigest()[:16]


def golden(corpus=None):
    """`{id: hash}` of the cascade's drawing of every covered molecule of the corpus."""
    from ...core import smiles
    from .molecule import layout2d
    from .scorecard import CORPUS, read_corpus
    out = {}
    for row in read_corpus(corpus or CORPUS):
        if '>' in row['smiles']:
            continue
        mol = smiles(row['smiles'])
        if covered(mol):
            out[row['id']] = drawing_hash(mol, layout2d(mol, engine='chython', force=True, chain=False))
    return out


def bless(path=GOLDEN, corpus=None):
    """Write `GOLDEN`: the cascade's drawings of the covered corpus, frozen."""
    rows = golden(corpus)
    Path(path).write_text('id\thash\n' + ''.join(f'{i}\t{h}\n' for i, h in sorted(rows.items())))
    return len(rows)


def main(argv=None):
    ap = ArgumentParser(prog='python -m chython.depict.layout.cascade')
    ap.add_argument('verb', choices=('bless',))
    ap.add_argument('--out', default=str(GOLDEN))
    a = ap.parse_args(argv)
    print(f'{bless(a.out)} covered drawings frozen in {a.out}')


if __name__ == '__main__':
    main()
