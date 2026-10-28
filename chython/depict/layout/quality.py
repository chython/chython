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
"""Layout quality, one number per metric, for any engine's plane.

| key | counts |
| --- | --- |
| `ovl` | atom pairs closer than `0.30 L` |
| `ovl61` | atom pairs closer than `0.61 L`, where labels collide |
| `crs` | non-adjacent bond pairs properly crossing |
| `aob` | atoms within `0.25 L` of a bond they are not in |
| `lbl` | bonds through the label box of an atom they are not in |
| `lin` | atoms of `linear_atoms` whose two bonds are more than 5 degrees off straight (in `shape`) |
| `inr` | atoms inside a ring of at most 8 of another ring block, a cage's own atoms exempt (in `shape`) |

The metrics of stability and convention follow below the defects. Every function takes a plane already
rescaled by `rescaled()`.
"""
from itertools import combinations, islice
from math import atan2, cos, degrees, fsum, hypot, pi, sin, sqrt

from ...core import smiles
from ...core.wedge import cis_trans_parity
from ..bonds import segment_hits_box
from ..label import labels
from ..style import get_depict_style
from .rings import ring_blocks


L = .825


def linear_atoms(mol):
    """Atoms drawn straight: two neighbours and a triple bond, or a carbon or nitrogen of two double bonds."""
    out = set()
    for x in mol:
        if len(mol.neighbors_of(x)) != 2:
            continue
        h = mol.hybridization_of(x)
        if h == 3 or h == 5 and mol.atom(x).atomic_symbol in ('C', 'N'):
            out.add(x)
    return out


def rescaled(mol, plane):
    """A copy of `plane` at mean bond `L`, in any unit; unchanged when `mol` has no bond of nonzero length."""
    bonds = [hypot(plane[b.n][0] - plane[b.m][0], plane[b.n][1] - plane[b.m][1]) for b in mol.bonds()]
    mean = fsum(bonds) / len(bonds) if bonds else 0.
    if mean < 1e-9:
        return dict(plane)
    return {x: (p[0] * L / mean, p[1] * L / mean) for x, p in plane.items()}


def _cross(p, q, a, b):
    def side(o, u, v):
        return (u[0] - o[0]) * (v[1] - o[1]) - (u[1] - o[1]) * (v[0] - o[0])
    return side(p, q, a) * side(p, q, b) < -1e-12 and side(a, b, p) * side(a, b, q) < -1e-12     # collinear is not


def _seg_dist(x, p, q):
    dx, dy = q[0] - p[0], q[1] - p[1]
    d = dx * dx + dy * dy
    t = 0. if not d else max(0., min(1., ((x[0] - p[0]) * dx + (x[1] - p[1]) * dy) / d))
    return hypot(p[0] + t * dx - x[0], p[1] + t * dy - x[1])


def defects(mol, plane, style=None, routed=()):
    """Counts of the geometric defects of the module table; a plane too collapsed to label scores `lbl` 0.

    `routed` holds `(n, m)` bonds drawn as routed paths rather than segments (`depict.peptide`), left out."""
    atoms = list(mol)
    routed = {frozenset(x) for x in routed}
    bonds = [(b.n, b.m) for b in mol.bonds() if frozenset((b.n, b.m)) not in routed]
    out = {'ovl': 0, 'ovl61': 0, 'crs': 0, 'aob': 0, 'lbl': 0}
    for a, b in combinations(atoms, 2):
        d = hypot(plane[a][0] - plane[b][0], plane[a][1] - plane[b][1])
        out['ovl'] += d < .30 * L
        out['ovl61'] += d < .61 * L
    for (a, b), (c, e) in combinations(bonds, 2):
        if len({a, b, c, e}) == 4 and _cross(plane[a], plane[b], plane[c], plane[e]):
            out['crs'] += 1
    for a, b in bonds:
        for x in atoms:
            if x != a and x != b and _seg_dist(plane[x], plane[a], plane[b]) < .25 * L:
                out['aob'] += 1
    if bonds and len(atoms) > 2:
        try:
            boxes = labels(mol, plane, style or get_depict_style())
        except ValueError:                                  # a degenerate span: already counted by `ovl`
            return out
        for a, b in bonds:
            for x, lb in boxes.items():
                if x != a and x != b and lb.text is not None and segment_hits_box(plane[a], plane[b], lb.box):
                    out['lbl'] += 1
    return out


def _polygon_rmsd(pts):
    """RMSD of `pts`, in ring order, from the regular polygon of side `L`: rigid fit, both windings. The radius is
    fixed, so a collapsed ring is far from its polygon rather than a perfect polygon of radius 0."""
    n = len(pts)
    cx, cy = fsum(p[0] for p in pts) / n, fsum(p[1] for p in pts) / n
    ps = [(p[0] - cx, p[1] - cy) for p in pts]
    best = None
    for w in (1, -1):
        qs = [(cos(2 * pi * w * k / n), sin(2 * pi * w * k / n)) for k in range(n)]
        a = fsum(p[0] * q[0] + p[1] * q[1] for p, q in zip(ps, qs))
        b = fsum(p[1] * q[0] - p[0] * q[1] for p, q in zip(ps, qs))
        s = L / (2 * sin(pi / n))
        h = hypot(a, b)
        c, t = (a / h, b / h) if h else (1., 0.)
        r = fsum((p[0] - s * (c * q[0] - t * q[1])) ** 2 + (p[1] - s * (t * q[0] + c * q[1])) ** 2
                 for p, q in zip(ps, qs))
        best = r if best is None else min(best, r)
    return sqrt(best / n) / L


def shape(mol, plane):
    lengths = [hypot(plane[b.n][0] - plane[b.m][0], plane[b.n][1] - plane[b.m][1]) for b in mol.bonds()]
    out = {'cv': 0., 'sqz': 0, 'ang': 0., 'rng': 0., 'axis': None, 'ct': 0, 'lin': 0, 'inr': 0}
    linear = linear_atoms(mol)
    if lengths:
        mean = fsum(lengths) / len(lengths)
        if mean:
            out['cv'] = sqrt(fsum((x - mean) ** 2 for x in lengths) / len(lengths)) / mean
    for x in mol:
        nb = list(mol.neighbors_of(x))
        if len(nb) < 2:
            continue
        angs = sorted(degrees(atan2(plane[y][1] - plane[x][1], plane[y][0] - plane[x][0])) % 360 for y in nb)
        gaps = [(angs[(q + 1) % len(angs)] - angs[q]) % 360 for q in range(len(angs))]
        out['sqz'] += min(gaps) < 30
        if len(nb) == 2 and x in linear:
            out['lin'] += abs(min(gaps) - 180) > 5
        if not mol.in_ring_of(x) and x not in linear and len(nb) in (2, 3):
            dev = abs(min(gaps) - 120) if len(nb) == 2 else max(abs(g - 120) for g in gaps)
            out['ang'] = max(out['ang'], dev)
    rings = [r for r in mol.sssr if 3 <= len(r) <= 8]
    if rings:
        out['rng'] = fsum(_polygon_rmsd([plane[x] for x in r]) for r in rings) / len(rings)
        blocks = [b for b, _ in ring_blocks(mol)]
        for r in rings:
            own = next(b for b in blocks if b.issuperset(r))
            poly = [plane[x] for x in r]
            out['inr'] += sum(_inside(plane[x], poly) for x in mol if x not in own)
    pts = list(plane.values())
    if len(pts) > 1:
        cx, cy = fsum(p[0] for p in pts) / len(pts), fsum(p[1] for p in pts) / len(pts)
        sxx = fsum((p[0] - cx) ** 2 for p in pts)
        syy = fsum((p[1] - cy) ** 2 for p in pts)
        sxy = fsum((p[0] - cx) * (p[1] - cy) for p in pts)
        tr, det = sxx + syy, sxx * syy - sxy * sxy
        disc = sqrt(max(0., tr * tr / 4 - det))
        lo, hi = tr / 2 - disc, tr / 2 + disc
        if hi > 1.5 * max(lo, 1e-12):
            a = abs(degrees(.5 * atan2(2 * sxy, sxx - syy))) % 180
            out['axis'] = min(a, 180 - a)
    for u in mol.stereo_units():
        if u['kind'] == 1 and u['parity'] and cis_trans_parity(mol, u, plane=plane) != u['parity']:
            out['ct'] += 1
    return out


def _inside(p, poly):
    """True when `p` is inside the polygon `poly`, by the even-odd rule."""
    c = False
    for (x1, y1), (x2, y2) in zip(poly[-1:] + poly[:-1], poly):
        if (y1 > p[1]) != (y2 > p[1]) and p[0] < (x2 - x1) * (p[1] - y1) / (y2 - y1) + x1:
            c = not c
    return c


def perm(mol, lay, k=10):
    """Atoms of `k` random-order re-reads of `mol` whose `lay` position has no same-element atom of the
    reference's within `1e-3 L`, each component's centroid aligned: zero for a layout blind to atom order.
    Where components sit on the page is not scored.

    The reference is `mol` re-read from its own SMILES, so it carries what every re-read carries: a parity
    stored on a double bond the stereo perception calls labile is no part of the molecule's identity and is
    not written, so `mol` as parsed would be compared against re-reads holding less than it does."""
    def centred(m, p):
        p = rescaled(m, p)
        out = []
        for c in m.connected_components:
            cx, cy = fsum(p[x][0] for x in c) / len(c), fsum(p[x][1] for x in c) / len(c)
            out.extend((m.element_of(x), p[x][0] - cx, p[x][1] - cy) for x in c)
        return out
    mol = smiles(format(mol))
    ref = centred(mol, lay(mol))
    miss = 0
    for _ in range(k):
        other = smiles(format(mol, 'r'))
        free = list(ref)
        for e, x, y in centred(other, lay(other)):
            j = next((j for j, (f, u, v) in enumerate(free) if f == e and hypot(u - x, v - y) < 1e-3 * L), None)
            if j is None:
                miss += 1
            else:
                free.pop(j)
    return miss


def shift_rmsd(pairs):
    n = len(pairs)
    dx = fsum(q[0] - p[0] for p, q in pairs) / n
    dy = fsum(q[1] - p[1] for p, q in pairs) / n
    return sqrt(fsum((q[0] - p[0] - dx) ** 2 + (q[1] - p[1] - dy) ** 2 for p, q in pairs) / n) / L


def _automorphisms(mol, cap):
    return islice(mol.as_query().get_mapping(mol), cap)


def _ring_systems(mol):
    systems = []
    for ring in mol.sssr:
        ring = set(ring)
        for s in [s for s in systems if s & ring]:
            systems.remove(s)
            ring |= s
        systems.append(ring)
    return systems


def pose(mol, plane, ref, cap=1000):
    """Distance, in `L`, of each ring system of `plane` from the reviewed `ref`, the best automorphism,
    translation only; an acyclic molecule is one unit."""
    units = sorted(sorted(s) for s in _ring_systems(mol)) or [list(mol)]
    best = [float('inf')] * len(units)
    for m in _automorphisms(mol, cap):
        for k, u in enumerate(units):
            best[k] = min(best[k], shift_rmsd([(ref[x], plane[m[x]]) for x in u]))
    return best


def _environment(mol, x):
    """What a mapped atom must keep to count as unchanged: its element and each mapped neighbour's bond."""
    env = []
    for b in mol.bonds():
        if x in (b.n, b.m):
            y = b.m if b.n == x else b.n
            env.append((mol.map_number_of(y), int(b.order)))
    return mol.element_of(x), tuple(sorted(env))


def rxn(reaction, planes):
    """RMSD, in `L`, of the mapped atoms whose bonds do not change, each reactant against each product;
    None when no pair of molecules shares two such atoms, or a map number repeats on one side."""
    mols = list(reaction.molecules())
    nr, np_ = len(reaction.reactants), len(reaction.products)
    sides = [(mols[:nr], planes[:nr]), (mols[len(mols) - np_:], planes[len(mols) - np_:])]
    index = []
    for ms, ps in sides:
        seen = {}
        for m, p in zip(ms, ps):
            for x in m:
                k = m.map_number_of(x)
                if not k:
                    continue
                if k in seen:
                    return None
                seen[k] = (id(m), _environment(m, x), p[x])
        index.append(seen)
    groups = {}
    for k, (ma, ea, pa) in index[0].items():
        if k in index[1] and index[1][k][1] == ea:
            groups.setdefault((ma, index[1][k][0]), []).append((pa, index[1][k][2]))
    pairs = [g for g in groups.values() if len(g) >= 2]
    if not pairs:
        return None
    return sqrt(fsum(shift_rmsd(g) ** 2 * len(g) for g in pairs) / sum(len(g) for g in pairs))


def ser(core, members, cap=50):
    """RMSD, in `L`, of `core` in each member against the first member, the best embedding each; None when
    fewer than two members hold the core."""
    q = core.as_query()
    ref = None
    out = []
    for mol, plane in members:
        embeds = list(islice(q.get_mapping(mol), cap))
        if not embeds:
            continue
        if ref is None:
            ref = {c: plane[x] for c, x in embeds[0].items()}
            continue
        out.append(min(shift_rmsd([(ref[c], plane[x]) for c, x in e.items()]) for e in embeds))
    return fsum(out) / len(out) if out else None
