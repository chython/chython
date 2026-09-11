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
"""Bond geometry: trimmed analytically against the labels, chained into paths, one shape per order.

Nothing here paints over anything -- a bond stops at the label's box, so no `<mask>` is needed and all
three backends can express it.  Every length and colour comes from `DepictStyle`; the four module
constants below are construction limits rather than taste.  Coordinates are y-up, in molecule units,
unflipped: the flip belongs to the backend.
"""
from math import fsum, hypot

from ..core import LogRecord
from .scene import Box, EMPTY_BOX, Path, circle, polyline
from .style import DepictStyle


__all__ = ['bond_paths', 'chains', 'has_ink', 'inner_line', 'ray_box_exit', 'segment_hits_box', 'trim',
           'trim_ink', 'trim_per_end']


Point = tuple[float, float]
Segment = tuple[Point, Point]

# How far off a bond's own line the ring centroid may sit before its inner line is refused: the point at
# which `offset / |perpendicular|` stops being a projection and starts being an extrapolation.
_COLLINEAR_LIMIT = .65

# "The same ring corner" versus "a gap".  Two inner lines meeting at a corner agree to about 1e-5 on a
# `clean2d()` layout; 1e-3 is three orders above that and far below the smallest real break, `bond.trim`.
_MERGE_TOLERANCE = 1e-3

# The shortest dash a pattern may be compensated down to -- a renderability limit, not a look.  `Path`
# refuses a non-positive dash length, and round-cap compensation can drive a fine pattern negative.
_MIN_DASH = .01

# Below this a bond has no derivable direction (two atoms on one point, which a layout can produce).
# Far below any visible distance, because the answer it guards is `None`: a bond dropped from the picture.
_MIN_LENGTH = 1e-9


def ray_box_exit(origin: Point, direction: Point, box: Box) -> float:
    """How far along `direction` a ray from `origin` travels before it leaves `box`.  0 if it starts out.

    The slab method.  `direction` need not be normalized; the result is in units of `direction`.
    """
    x, y = origin
    dx, dy = direction
    if box.min_x > box.max_x:            # the empty box
        return 0.
    if not (box.min_x <= x <= box.max_x and box.min_y <= y <= box.max_y):
        return 0.                        # the ray starts outside: nothing to skip
    far = float('inf')
    if dx:
        far = min(far, ((box.max_x if dx > 0. else box.min_x) - x) / dx)
    if dy:
        far = min(far, ((box.max_y if dy > 0. else box.min_y) - y) / dy)
    return 0. if far == float('inf') else max(far, 0.)


def trim(p: Point, q: Point, box_p: Box, box_q: Box, clearance: float) -> Segment | None:
    """The visible part of the segment `p`->`q`, or None when the labels leave nothing.

    Both ends move: the start out of `box_p`, the end back out of `box_q`, and each by a further
    `clearance` so the stroke does not touch the ink.  `None` rather than a zero-length segment, because
    a zero-length stroke with a round cap is a dot and a reader takes a dot for a radical.  `clearance`
    is spent unconditionally, even at an unlabelled end; :func:`trim_ink` is the ink-conditional form.
    """
    return trim_per_end(p, q, box_p, box_q, clearance, clearance)


def trim_per_end(p: Point, q: Point, box_p: Box, box_q: Box, clearance_p: float,
                 clearance_q: float) -> Segment | None:
    """`trim` with the two clearances stated apart.  One geometry for three callers that differ.

    A chain's interior vertex takes no clearance (the stroke goes through it), a ring's inner line takes
    it only where the vertex carries ink, and a wedge takes none at its narrow end.  Public because
    `wedge.py` is the third caller.
    """
    dx = q[0] - p[0]
    dy = q[1] - p[1]
    length = hypot(dx, dy)
    if length < _MIN_LENGTH:
        return None
    ux, uy = dx / length, dy / length
    head = ray_box_exit(p, (ux, uy), box_p) + clearance_p
    tail = ray_box_exit(q, (-ux, -uy), box_q) + clearance_q
    if head + tail >= length:
        return None
    return (p[0] + ux * head, p[1] + uy * head), (q[0] - ux * tail, q[1] - uy * tail)


def segment_hits_box(p: Point, q: Point, box: Box) -> bool:
    """Does the segment `p`-`q` touch `box`?  Liang-Barsky, so a diagonal miss is a miss.

    A bounding-box test would call every long diagonal bond a hit on every annotation in its corner.

    Here rather than beside its callers because both of them ask the same question about the same thing:
    `label.py` asks it of the bond axes to choose where a map number goes, and `figure.py` asks it of the
    paths actually drawn to decide whether that number needs a plate under it.  Two copies would answer
    differently the first time one of them learned about curves.
    """
    x, y = p
    dx = q[0] - x
    dy = q[1] - y
    if (max(x, q[0]) < box.min_x or min(x, q[0]) > box.max_x
            or max(y, q[1]) < box.min_y or min(y, q[1]) > box.max_y):
        return False                     # cheap rejection first: most bonds are nowhere near
    t0, t1 = 0., 1.
    for delta, low, near, far in ((dx, x, box.min_x, box.max_x), (dy, y, box.min_y, box.max_y)):
        if abs(delta) < _MIN_LENGTH:
            if low < near or low > far:
                return False
            continue
        a = (near - low) / delta
        b = (far - low) / delta
        if a > b:
            a, b = b, a
        if a > t0:
            t0 = a
        if b < t1:
            t1 = b
        if t0 > t1:
            return False
    return True


def inner_line(p: Point, q: Point, centroid: Point, offset: float) -> Segment | None:
    """The line `offset` inside the bond `p`->`q`, on the centroid's side, as its two points.

    Serves both the ring double bond's second line and the `'dashed-inner'` arc.  Each end is shortened
    along the vertex bisector to where the inset line crosses the vertex-to-centroid line, which on a
    regular ring is where the neighbouring inner line arrives -- so a run closes with no spur and no
    adjacency lookup.  The amount depends on the vertex angle, so no constant suits a five- and a
    six-ring at once.

    `None` when the centroid is nearly collinear with the bond (the inset explodes) or when a skewed ring
    projects the crossing outside the bond's own footprint, where the line would read as another ring.
    """
    dx, dy = q[0] - p[0], q[1] - p[1]
    length = hypot(dx, dy)
    if length < _MIN_LENGTH:
        return None
    ux, uy = dx / length, dy / length
    # the centroid in the bond's own frame: +x along p->q, +y to its left
    cx, cy = centroid[0] - p[0], centroid[1] - p[1]
    along = cx * ux + cy * uy
    across = -cx * uy + cy * ux
    if not across or offset / abs(across) >= _COLLINEAR_LIMIT:
        return None
    side = offset if across > 0. else -offset
    across = abs(across)
    head = offset * along / across                          # where p's bisector crosses the inset line
    tail = length - offset * (length - along) / across      # and q's
    head = min(max(head, 0.), length)                       # a skewed ring projects them outside the
    tail = min(max(tail, 0.), length)                       # bond: clamp to its own footprint
    if tail <= head:
        return None
    return ((p[0] + head * ux - side * uy, p[1] + head * uy + side * ux),
            (p[0] + tail * ux - side * uy, p[1] + tail * uy + side * ux))


def _key(n: int, m: int) -> tuple[int, int]:
    """A bond's key: the two stable ids, low first.  `skip`, `widths` and `colours` are keyed this way."""
    return (n, m) if n < m else (m, n)


def _adjacency(mol) -> dict[int, dict[int, int]]:
    """`{n: {neighbour: order}}` from ONE pass over the bonds."""
    graph = {atom.n: {} for atom in mol.atoms()}
    for bond in mol.bonds():
        graph[bond.n][bond.m] = bond.order
        graph[bond.m][bond.n] = bond.order
    return graph


def _check_boxes(mol, boxes):
    """`labels()` returns one `Label` per atom; anything less would fail as a bare `KeyError` in a walk."""
    for sid in mol:
        if sid not in boxes:
            raise ValueError('boxes has no entry for atom %d: pass the whole dict labels() returned, '
                             'which has an entry for every atom including the unlabelled ones' % sid)


def _check_plane(mol, plane):
    """The same guard for the coordinates."""
    for sid in mol:
        if sid not in plane:
            raise ValueError('plane has no coordinates for atom %d: pass the whole mapping '
                             'coordinates() returned, which has an entry for every atom' % sid)


def chains(mol, boxes, plain=None) -> list[tuple[int, ...]]:
    """Runs of atoms whose bonds can be drawn as ONE stroke, longest first.

    `plain` is the set of bond keys (low-first) the notation draws as exactly one plain stroke along the
    bond axis, at this width and colour; `None` means the order-1 bonds.  `bond_paths` computes it,
    because it is the notation's business and not the order's (under `'circle'` and `'dashed-inner'` an
    order-4 ring bond draws an order-1 bond's shape).

    A run passes through an atom only when it is a bare vertex of degree 2 with two plain bonds.  Every
    plain bond appears in exactly one run; a closed run repeats its first atom last, so a caller closes
    the path rather than stacking two caps where a join belongs.
    """
    _check_boxes(mol, boxes)
    graph = _adjacency(mol)
    if plain is None:
        plain = {_key(n, m) for n, neighbours in graph.items()
                 for m, order in neighbours.items() if order == 1}
    passable = {sid for sid, neighbours in graph.items()
                if len(neighbours) == 2 and boxes[sid].text is None
                and all(_key(sid, other) in plain for other in neighbours)}
    remaining = set(plain)
    out = []
    # every run that has an end starts at one: an atom a stroke cannot pass through
    for start in sorted(sid for sid in graph if sid not in passable):
        for first in sorted(graph[start]):
            if _key(start, first) not in remaining:
                continue
            remaining.discard(_key(start, first))
            run = [start, first]
            previous, current = start, first
            while current in passable:
                following = next(x for x in graph[current] if x != previous)
                step = _key(current, following)
                if step not in remaining:
                    break
                remaining.discard(step)
                run.append(following)
                previous, current = current, following
            out.append(tuple(run))
    # what is left is a cycle of bare degree-2 vertices -- a ring with no label and no branch on it
    while remaining:
        start, first = min(remaining)
        remaining.discard((start, first))
        run = [start, first]
        previous, current = start, first
        while current != start:
            following = next(x for x in graph[current] if x != previous)
            remaining.discard(_key(current, following))
            run.append(following)
            previous, current = current, following
        out.append(tuple(run))
    out.sort(key=len, reverse=True)
    return out


def bond_paths(mol, plane, boxes, style: DepictStyle, *, skip=frozenset(), widths=None, colours=None,
               log=None) -> list[Path]:
    """Every non-stereo bond, as `Path` objects.  Stereo bonds are `wedge.py`'s.

    Order by order:
      1  a line, chained with its neighbours where it can be
      2  two lines, `bond.spacing` apart -- offset INWARD in a ring, straddling the axis otherwise
      3  three lines, the middle one on the axis
      4  aromatic: `bond.aromatic` chooses kekule (default), circle or dashed-inner
      8  dative: one line, DASHED (`bond.dative_dashes`), with no head -- the container does not record
         which atom donates, so a head would be a claim picked from atom order

    A bond whose two labels leave nothing to draw is skipped and reported through `log` -- see `trim`.
    `skip` names bonds somebody else draws (a wedge draws its own).  `widths` and `colours` override the
    style per bond, keyed low-first; a bond differing from its chain neighbour breaks the chain, since
    one path cannot taper or change colour half way along.
    """
    _check_boxes(mol, boxes)
    _check_plane(mol, plane)
    bond_style = style.bond
    skipped = frozenset(_key(*key) for key in skip)
    per_bond_width = {} if widths is None else {_key(*k): float(v) for k, v in widths.items()}
    per_bond_colour = {} if colours is None else {_key(*k): v for k, v in colours.items()}
    graph = _adjacency(mol)
    orders = {_key(n, m): order for n, neighbours in graph.items() for m, order in neighbours.items()}
    centroids = _bond_centroids(mol, plane)

    def paint(key):
        return per_bond_width.get(key, bond_style.width), per_bond_colour.get(key, bond_style.colour)

    # which bonds alternate is asked of the core, on a throwaway copy -- see `_kekule_doubles`
    doubles = frozenset()
    unresolved = frozenset()
    if bond_style.aromatic == 'kekule':
        doubles, unresolved = _kekule_doubles(mol, orders, log)

    # which bonds are one plain stroke, not which order they are: an order-4 bond the notation gave no
    # second line draws what an order-1 bond draws.  `doubles` carries the whole difference between the
    # three aromatic notations, so one expression covers all of them.
    plain = {key for key, order in orders.items()
             if order == 1 or (order == 4 and key not in doubles)}

    out = []
    for run in chains(mol, boxes, plain):
        for ids, closed, attrs in _split(run, paint, skipped):
            width, colour = attrs
            for points, is_closed in _strokes(ids, closed, plane, boxes, bond_style.trim, log):
                out.append(_path([polyline(points, closed=is_closed)], width, colour, bond_style))

    for key in sorted(orders):
        order = orders[key]
        if key in plain or key in skipped:
            continue
        n, m = key
        p, q = plane[n], plane[m]
        box_p, box_q = boxes[n].box, boxes[m].box
        width, colour = paint(key)
        axis = trim_ink(p, q, box_p, box_q, bond_style.trim, log=log, atoms=key)
        if axis is None:
            _crowded(log, n, m)
            continue
        dashes = None
        if order == 2 or (order == 4 and key in doubles):
            lines = _double_lines(axis, p, q, centroids.get(key), box_p, box_q, bond_style, log, n, m)
        elif order == 3:
            lines = _triple_lines(axis, bond_style.triple_spacing)
        elif order == 8:
            lines = [axis]
            dashes = _dash_pattern(bond_style, bond_style.dative_dashes)
        else:
            # any order this module does not model, drawn as the plain line the arena guarantees is a
            # contact.  An order-4 bond never lands here: it is in `doubles` or in `plain`.
            lines = [axis]
        out.append(_path([polyline(pair) for pair in lines], width, colour, bond_style, dashes=dashes))

    if bond_style.aromatic == 'circle':
        for ring in mol.aromatic_rings:
            ornament = _circle_path(ring, plane, bond_style)
            if ornament is not None:
                out.append(ornament)
    elif bond_style.aromatic == 'dashed-inner':
        for ring in mol.aromatic_rings:
            ornament = _dashed_path(ring, plane, boxes, bond_style)
            if ornament is not None:
                out.append(ornament)
    elif unresolved:
        # 'kekule' asked for, but an aromatic system has no Kekule form: its rings fall back to the
        # circle, which commits to nothing beyond "aromatic".  Single lines would read as a cyclopentane.
        for ring in mol.aromatic_rings:
            if all(sid in unresolved for sid in ring):
                ornament = _circle_path(ring, plane, bond_style)
                if ornament is not None:
                    out.append(ornament)
    return out


def _path(subpaths, width: float, colour: str, bond_style, dashes=None) -> Path:
    """The one place a bond's `Path` is built, so cap, join and miter limit come off the style once."""
    return Path(subpaths, stroke=colour, width=width, dashes=dashes, cap=bond_style.cap,
                join=bond_style.join, miter_limit=bond_style.miter_limit)


def _crowded(log, n: int, m: int):
    if log is not None:
        log.append(LogRecord('depict:crowded', (n, m),
                             'the labels on both atoms leave no room for the bond; it was not drawn'))


def _no_room_inside(log, n: int, m: int):
    if log is not None:
        log.append(LogRecord('depict:crowded-inner', (n, m),
                             'the labels leave no room for the inner line of this multiple bond; '
                             'only the bond axis was drawn'))


def _tight(log, atoms):
    """Rung 2 of `trim_ink` fired: drawn closer to the ink than the style asked for.

    A distinct id from `depict:crowded`, which means the bond was lost altogether.
    """
    if log is not None:
        log.append(LogRecord('depict:tight', tuple(atoms),
                             'the labels leave less room than bond.trim asks for; the line was drawn '
                             'up to the label box with no clearance rather than dropped'))


def has_ink(box: Box) -> bool:
    """Does this box hold a glyph?  An unlabelled atom's box is a point and the empty box is inverted."""
    return box.max_x > box.min_x and box.max_y > box.min_y


def trim_ink(a: Point, b: Point, box_p: Box, box_q: Box, clearance: float, *, log=None,
             atoms=()) -> Segment | None:
    """A segment trimmed only where the vertex actually carries ink -- the rule for every drawn line here.

    `bond.trim` is clearance from a label's ink box, so spending it at a bare vertex would punch a hole in
    the drawing (14.5% of a default bond, at every ring corner and branch point).  Public `trim()` keeps
    its unconditional semantics for callers that want it.

    Two rungs, because clearance is a preference and not a constraint.  Rung 1 is the ink-conditional
    clearance; when that leaves nothing, rung 2 spends none and tries again -- still stopping at the
    label's box, so the stroke never touches a glyph -- and only then does the refusal stand.  Otherwise
    hydrogen peroxide at the `acs` preset draws as two unbonded oxygens: a bond that exists must never be
    drawn as nothing.  Rung 2 puts a `depict:tight` record in `log`; `log` and `atoms` are both optional.
    """
    clearance_p = clearance if has_ink(box_p) else 0.
    clearance_q = clearance if has_ink(box_q) else 0.
    segment = trim_per_end(a, b, box_p, box_q, clearance_p, clearance_q)
    if segment is not None or not (clearance_p or clearance_q):
        return segment      # rung 1 answered, or spent nothing so rung 2 asks the same question
    segment = trim_per_end(a, b, box_p, box_q, 0., 0.)
    if segment is not None:
        _tight(log, atoms)
    return segment


def _split(run, paint, skipped):
    """One chain, cut where a caller's override or a skip says the stroke cannot continue.

    Yields `(ids, closed, attrs)`.  A closed run stays closed only when every one of its bonds is drawn
    and all of them share a width and a colour; otherwise it is rotated so a cut falls at the ends and
    comes back as one or more open runs.
    """
    closed = run[0] == run[-1]
    ids = list(run[:-1]) if closed else list(run)
    count = len(ids)
    if closed:
        edges = [(ids[i], ids[(i + 1) % count]) for i in range(count)]
    else:
        edges = [(ids[i], ids[i + 1]) for i in range(count - 1)]
    attrs = [None if _key(*edge) in skipped else paint(_key(*edge)) for edge in edges]
    if closed and all(attr is not None and attr == attrs[0] for attr in attrs):
        return [(tuple(ids), True, attrs[0])]
    if closed:   # rotate so index 0 begins a run rather than continuing one
        start = next(i for i in range(count)
                     if attrs[i] is None or attrs[i - 1] is None or attrs[i] != attrs[i - 1])
        edges = edges[start:] + edges[:start]
        attrs = attrs[start:] + attrs[:start]
    out = []
    current = []
    current_attrs = None
    for edge, attr in zip(edges, attrs):
        if attr is None:
            if current:
                out.append((tuple(current), False, current_attrs))
                current = []
            continue
        if current and attr == current_attrs:
            current.append(edge[1])
        else:
            if current:
                out.append((tuple(current), False, current_attrs))
            current = [edge[0], edge[1]]
            current_attrs = attr
    if current:
        out.append((tuple(current), False, current_attrs))
    return out


def _strokes(ids, closed, plane, boxes, clearance, log):
    """The points of one stroke: trimmed at its two ends only, and through every vertex between them.

    An interior vertex is bare by construction, so the stroke passes exactly through the atom point and
    the corner is a join rather than two butt caps.  The two ends go through `trim_ink`, so a run ending
    on a bare vertex is not shortened there either.  A terminal bond the labels swallow is dropped from
    the run and reported.

    A run may come back to its start and still not close: the walk stops at a labelled atom, so a ring
    interrupted only by one label arrives as a cycle whose joining atom carries a glyph.  Closing it would
    draw the stroke through that glyph, so it becomes one open run trimmed at both ends against that box.
    """
    if closed and not has_ink(boxes[ids[0]].box):
        return [([plane[sid] for sid in ids], True)]
    ids = ([*ids, ids[0]] if closed else list(ids))
    while len(ids) > 2:
        if trim_ink(plane[ids[0]], plane[ids[1]], boxes[ids[0]].box, EMPTY_BOX, clearance) is not None:
            break
        _crowded(log, ids[0], ids[1])
        ids = ids[1:]
    while len(ids) > 2:
        if trim_ink(plane[ids[-2]], plane[ids[-1]], EMPTY_BOX, boxes[ids[-1]].box,
                    clearance) is not None:
            break
        _crowded(log, ids[-2], ids[-1])
        ids = ids[:-1]
    # the two probes above carry no `log`: the same question is asked again below to build the geometry,
    # and the surviving call is the one that discloses.
    if len(ids) == 2:
        segment = trim_ink(plane[ids[0]], plane[ids[1]], boxes[ids[0]].box, boxes[ids[1]].box,
                           clearance, log=log, atoms=_key(ids[0], ids[1]))
        if segment is None:
            _crowded(log, ids[0], ids[1])
            return []
        return [([segment[0], segment[1]], False)]
    head = trim_ink(plane[ids[0]], plane[ids[1]], boxes[ids[0]].box, EMPTY_BOX, clearance, log=log,
                    atoms=_key(ids[0], ids[1]))
    tail = trim_ink(plane[ids[-2]], plane[ids[-1]], EMPTY_BOX, boxes[ids[-1]].box, clearance, log=log,
                    atoms=_key(ids[-2], ids[-1]))
    points = [head[0]] + [plane[sid] for sid in ids[1:-1]] + [tail[1]]
    return [(points, False)]


def _bond_centroids(mol, plane) -> dict[tuple[int, int], Point]:
    """Each ring bond's inward direction, as the centroid of the smallest ring that holds it.

    Smallest because that is the ring a chemist reads the bond as belonging to: a fused bond's inner line
    goes inside its own six-ring, not inside the ten-membered perimeter.
    """
    out = {}
    for ring in sorted(mol.rings, key=len):
        centroid = (fsum(plane[sid][0] for sid in ring) / len(ring),
                    fsum(plane[sid][1] for sid in ring) / len(ring))
        for a, b in zip(ring, ring[1:] + ring[:1]):
            out.setdefault(_key(a, b), centroid)
    return out


def _kekule_doubles(mol, orders, log=None) -> tuple[frozenset, frozenset]:
    """Which aromatic bonds get the second line -- asked of `kekule()`, on a throwaway copy.

    Returns `(doubles, unresolved)`: the bond keys that alternate, and the stable ids of every atom in an
    aromatic system with no Kekule form at all.

    The copy is O(1) and preserves stable ids, so its orders are keyed like `orders`; `kekule()` repairs,
    and on a throwaway that repair is discarded with the copy.  `kekule_copy()` is the wrong door -- it
    raises on an unresolved system, taking a whole picture down over one bad ring.  An unresolved system is
    still rewritten with the best matching found, so its bonds are excluded here by their atoms rather than
    trusted; `bond_paths` draws those rings as circles and `log` says so.
    """
    copy = mol.copy()
    result = copy.kekule()
    unresolved = frozenset(sid for system in result.unresolved for sid in system)
    doubles = frozenset(key for key, order in orders.items()
                        if order == 4 and key[0] not in unresolved and key[1] not in unresolved
                        and copy.order_of(*key) == 2)
    if log is not None:
        for system in result.unresolved:
            log.append(LogRecord('depict:no-kekule', tuple(system),
                                 'this aromatic system has no Kekule form; no alternating second '
                                 'lines were drawn and its rings carry the aromatic circle instead'))
    return doubles, unresolved


def _double_lines(axis, p, q, centroid, box_p, box_q, bond_style, log, n, m):
    """A double bond's two lines: the axis and its partner, inside the ring when there is one."""
    if centroid is not None:
        inner = inner_line(p, q, centroid, bond_style.spacing)
        if inner is not None:
            trimmed = trim_ink(inner[0], inner[1], box_p, box_q, bond_style.trim, log=log,
                               atoms=_key(n, m))
            if trimmed is not None:
                return [axis, trimmed]
            _no_room_inside(log, n, m)
            return [axis]
        # the ring is there but its centroid is unusable (nearly collinear, or skewed past the bond):
        # straddling states no side, where an inset on a guessed side would state the wrong one
    return _straddle(axis, bond_style.spacing)


def _straddle(axis, spacing):
    """Two lines `spacing` apart, symmetric about the axis.  For a bond with no ring to lean into."""
    (ax, ay), (bx, by) = axis
    length = hypot(bx - ax, by - ay)
    ux, uy = (bx - ax) / length, (by - ay) / length
    dx, dy = -uy * spacing / 2., ux * spacing / 2.
    return [((ax + dx, ay + dy), (bx + dx, by + dy)), ((ax - dx, ay - dy), (bx - dx, by - dy))]


def _triple_lines(axis, spacing):
    """Three lines: one on the axis, one either side at `bond.triple_spacing`."""
    (ax, ay), (bx, by) = axis
    length = hypot(bx - ax, by - ay)
    ux, uy = (bx - ax) / length, (by - ay) / length
    dx, dy = -uy * spacing, ux * spacing
    return [((ax + dx, ay + dy), (bx + dx, by + dy)), axis,
            ((ax - dx, ay - dy), (bx - dx, by - dy))]


def _circle_path(ring, plane, bond_style) -> Path | None:
    """The aromatic circle: `bond.aromatic_inset` inside the ring's inscribed radius.

    The inscribed radius and not the vertex radius, because the inset is a gap from the BONDS -- a
    circle placed by the vertices would touch the bonds of any ring that is not equilateral.
    """
    centroid = (fsum(plane[sid][0] for sid in ring) / len(ring),
                fsum(plane[sid][1] for sid in ring) / len(ring))
    apothem = min(hypot((plane[a][0] + plane[b][0]) / 2. - centroid[0],
                        (plane[a][1] + plane[b][1]) / 2. - centroid[1])
                  for a, b in zip(ring, ring[1:] + ring[:1]))
    radius = apothem - bond_style.aromatic_inset
    if radius <= 0.:
        return None   # the inset swallowed the ring: no circle rather than a dot in the middle of it
    return _path([circle(centroid[0], centroid[1], radius)], bond_style.width, bond_style.colour,
                 bond_style)


def _dashed_path(ring, plane, boxes, bond_style) -> Path | None:
    """One ring's dashed inner arc: the inset lines, merged at the corners they share.

    A join and an SVG dash phase apply only within one subpath, so a ring emitted as six two-point
    subpaths anchors a dash to every corner and stacks two round caps where a join belongs.  Consecutive
    lines therefore merge into one subpath, closed when the walk comes back round.
    """
    centroid = (fsum(plane[sid][0] for sid in ring) / len(ring),
                fsum(plane[sid][1] for sid in ring) / len(ring))
    segments = []
    for a, b in zip(ring, ring[1:] + ring[:1]):
        segment = inner_line(plane[a], plane[b], centroid, bond_style.aromatic_dash_inset)
        if segment is not None:
            segment = trim_ink(segment[0], segment[1], boxes[a].box, boxes[b].box, bond_style.trim)
        segments.append(segment)
    runs = _merge_cycle(segments)
    if not runs:
        return None
    return _path([polyline(points, closed=closed) for points, closed in runs], bond_style.width,
                 bond_style.colour, bond_style,
                 dashes=_dash_pattern(bond_style, bond_style.aromatic_dashes))


def _merge_cycle(segments) -> list[tuple[list[Point], bool]]:
    """A ring's worth of inset lines, some of them missing, joined into as few runs as possible.

    `breaks` is empty exactly when nothing is missing and every corner met, which is the closed ring.
    """
    count = len(segments)
    breaks = [i for i in range(count)
              if segments[i] is None or segments[i - 1] is None
              or not _near(segments[i - 1][1], segments[i][0])]
    if not breaks:
        return [([segment[0] for segment in segments], True)]
    order = breaks[0]
    out = []
    current = []
    for step in range(count):
        segment = segments[(order + step) % count]
        if segment is None:
            if current:
                out.append((current, False))
                current = []
            continue
        if current and _near(current[-1], segment[0]):
            current.append(segment[1])
        else:
            if current:
                out.append((current, False))
            current = [segment[0], segment[1]]
    if current:
        out.append((current, False))
    return out


def _near(a: Point, b: Point) -> bool:
    return hypot(a[0] - b[0], a[1] - b[1]) < _MERGE_TOLERANCE


def _dash_pattern(bond_style, pattern) -> tuple[float, ...]:
    """A dash pattern compensated for a round cap.

    A round cap extends every dash by half the stroke width at each end, so a .15 dash renders .19 long:
    shrink the painted lengths by the width and grow the gaps by it.  A butt cap needs no compensation.
    `_MIN_DASH` is the floor because `Path` refuses a non-positive dash length.

    `pattern` is passed in because there are two notations -- `bond.aromatic_dashes` for the inner arc,
    `bond.dative_dashes` for the dative bond -- and the compensation is a property of the cap.
    """
    if bond_style.cap != 'round':
        return tuple(float(value) for value in pattern)
    return tuple(max(value - bond_style.width, _MIN_DASH) if not index % 2
                 else value + bond_style.width for index, value in enumerate(pattern))
