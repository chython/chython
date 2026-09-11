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
"""Bond geometry: trimmed analytically, chained into paths, one shape per order.

A line stops where the label's box starts rather than hiding under a mask, and a run of unlabelled
atoms is ONE path with a miter join rather than one `<line>` per bond.
"""
from math import hypot, isclose
from pytest import approx, raises
from chython import smiles
from chython.depict.bonds import (bond_paths, chains, inner_line, ray_box_exit, segment_hits_box, trim,
                                  trim_ink, trim_per_end, _kekule_doubles)
from chython.depict.label import labels
from chython.depict.scene import Box
from chython.depict.style import DepictStyle


# A fused/heteroaromatic corpus for the two geometry invariants below.  The counts the tests assert are
# taken over THIS list -- widen it and they move.
RINGS = ('c1ccccc1', 'c1ccncc1', 'c1ccoc1', 'c1cc[nH]c1', 'c1ccc2ccccc2c1', 'c1ccc2cc3ccccc3cc2c1',
         'c1ccc(-c2ccccc2)cc1', 'c1ccc2[nH]ccc2c1', 'c1ccc2ncccc2c1', 'Cc1ccccc1', 'Oc1ccccc1')

# `bond.aromatic` defaults to 'dashed-inner', so every test ABOUT the alternating lines asks for them.
KEKULE = DepictStyle().tuned(**{'bond.aromatic': 'kekule'})


def _setup(text, style=None):
    style = style or DepictStyle()
    mol = smiles(text)
    mol.clean2d()
    plane = mol.coordinates()
    return mol, plane, labels(mol, plane, style), style


def _doubles(mol):
    """The aromatic bonds the `'kekule'` notation gives a second line, as low-first keys."""
    orders = {(b.n, b.m) if b.n < b.m else (b.m, b.n): b.order for b in mol.bonds()}
    doubles, _ = _kekule_doubles(mol, orders)
    return doubles


def _points(subpath):
    """The endpoint of every segment in a subpath, as `(x, y)`; the argumentless `('Z',)` has none."""
    return [(s[-2], s[-1]) for s in subpath if len(s) > 1]


def _centroid(plane, ring):
    return (sum(plane[n][0] for n in ring) / len(ring), sum(plane[n][1] for n in ring) / len(ring))


def _inside(point, polygon):
    """Strict point-in-polygon by the crossing number.  `polygon` is a sequence of `(x, y)`."""
    x, y = point
    inside = False
    for i in range(len(polygon)):
        ax, ay = polygon[i - 1]
        bx, by = polygon[i]
        if (ay > y) != (by > y) and x < ax + (y - ay) * (bx - ax) / (by - ay):
            inside = not inside
    return inside


def _forced_plus_inner_line(p, q, centroid, offset):
    """THE CONTROL for `inner_line`: the same geometry with the side forced to +1.

    A dropped sign puts every inner line on the fixed side of its bond, which for half a ring is outside
    it.  This function IS that bug, so the assertions below can be shown to discriminate.
    """
    dx, dy = q[0] - p[0], q[1] - p[1]
    length = hypot(dx, dy)
    ux, uy = dx / length, dy / length
    cx, cy = centroid[0] - p[0], centroid[1] - p[1]
    cr_x, cr_y = cx * ux + cy * uy, -cx * uy + cy * ux
    if not cr_y or offset / abs(cr_y) >= .65:
        return None
    cr_y = abs(cr_y)                                   # the dropped sign: always the +y side
    a_x = offset * cr_x / cr_y
    b_x = length - offset * (length - cr_x) / cr_y
    a_x = min(max(a_x, 0.), length)
    b_x = min(max(b_x, 0.), length)
    if b_x <= a_x:
        return None
    return ((p[0] + a_x * ux - offset * uy, p[1] + a_x * uy + offset * ux),
            (p[0] + b_x * ux - offset * uy, p[1] + b_x * uy + offset * ux))


def test_a_ray_leaves_a_box_at_the_face_it_hits():
    """the whole trim, in one function: where does the line cross the label's box"""
    box = Box(-1., -1., 1., 1.)
    assert ray_box_exit((0., 0.), (1., 0.), box) == approx(1.)
    assert ray_box_exit((0., 0.), (0., 1.), box) == approx(1.)
    assert ray_box_exit((0., 0.), (1., 1.), box) == approx(1.), 'a diagonal leaves at the corner'


def test_a_ray_leaves_a_wide_box_by_its_nearest_face():
    """the slab method, not the corner: a wide label is left through its top, not its side"""
    assert ray_box_exit((0., 0.), (0., 1.), Box(-1., -.25, 1., .25)) == approx(.25)
    assert ray_box_exit((0., 0.), (1., 0.), Box(-1., -.25, 1., .25)) == approx(1.)
    assert ray_box_exit((0., 0.), (-1., 0.), Box(-.3, -.25, 1., .25)) == approx(.3), 'the near face'


def test_a_ray_from_outside_a_box_does_not_move_the_end():
    assert ray_box_exit((5., 5.), (1., 0.), Box(-1., -1., 1., 1.)) == approx(0.)


def test_a_ray_from_a_degenerate_box_travels_nothing():
    """an unlabelled atom's box is a point, so its bond starts exactly at the vertex"""
    assert ray_box_exit((0., 0.), (1., 0.), Box(0., 0., 0., 0.)) == approx(0.)


def test_a_ray_from_the_empty_box_travels_nothing():
    """`EMPTY_BOX` is inverted, and a min > max box must not read as "everywhere\""""
    from chython.depict.scene import EMPTY_BOX

    assert ray_box_exit((0., 0.), (1., 0.), EMPTY_BOX) == approx(0.)


def test_ray_box_exit_is_in_units_of_the_direction():
    """the docstring's promise: an unnormalized direction scales the answer"""
    box = Box(-1., -1., 1., 1.)
    assert ray_box_exit((0., 0.), (2., 0.), box) == approx(.5)


def test_a_segment_crossing_a_box_hits_it_and_one_stopping_short_does_not():
    box = Box(0., 0., 1., 1.)
    assert segment_hits_box((-1., .5), (2., .5), box)
    assert segment_hits_box((.4, .4), (.6, .6), box), 'a segment wholly inside is a hit'
    assert not segment_hits_box((-1., .5), (-.5, .5), box), 'stops before the box'
    assert not segment_hits_box((-1., 2.), (2., 2.), box), 'passes above it'


def test_a_segment_that_only_a_bounding_box_test_would_call_a_hit_is_a_miss():
    """the reason this is Liang-Barsky and not `Box.union`: a long diagonal bond past a number's corner

    Both endpoints are outside the box and its x and y spans both overlap the box's, so the two-interval
    test says yes -- but the line passes the corner without touching.
    """
    box = Box(0., 0., 1., 1.)
    assert not segment_hits_box((-1., .5), (.5, -1.), box)
    assert segment_hits_box((-1., 1.), (.5, -.5), box), 'the same diagonal moved onto the corner'


def test_a_segment_grazing_an_edge_counts_as_touching():
    """`touch`, as the docstring says: a bond lying exactly along the padded box's edge is worth plating"""
    box = Box(0., 0., 1., 1.)
    assert segment_hits_box((-1., 1.), (2., 1.), box)
    assert segment_hits_box((0., -1.), (0., 2.), box)


def test_a_degenerate_segment_is_a_point_test():
    """a curve contributes its endpoints, and a closed curve can hand in the same point twice"""
    box = Box(0., 0., 1., 1.)
    assert segment_hits_box((.5, .5), (.5, .5), box)
    assert not segment_hits_box((2., 2.), (2., 2.), box)


def test_nothing_hits_the_empty_box():
    """an annotation with no ink has no box, and must not plate the whole figure"""
    from chython.depict.scene import EMPTY_BOX

    assert not segment_hits_box((-1., 0.), (1., 0.), EMPTY_BOX)


def test_trim_shortens_both_ends_by_the_boxes_and_the_clearance():
    """A BOX IS ABSOLUTE, in molecule coordinates: the far box is written where the far atom is"""
    start, end = trim((0., 0.), (2., 0.), Box(-.2, -.2, .2, .2), Box(1.8, -.2, 2.2, .2), .05)
    assert start == approx((.25, 0.))
    assert end == approx((1.75, 0.))


def test_trim_leaves_an_unlabelled_end_alone():
    point = Box(0., 0., 0., 0.)
    start, end = trim((0., 0.), (1., 0.), point, Box(.8, -.2, 1.2, .2), 0.)
    assert start == approx((0., 0.))
    assert end == approx((.8, 0.))


def test_public_trim_spends_its_clearance_even_where_there_is_no_ink():
    """`trim()`'s contract is UNCONDITIONAL; the `has_ink` gate belongs to this module's own drawing

    A NON-ZERO clearance against two bare ends is the only shape that tells gated from ungated apart:
    at clearance `0.` the two give the same answer.
    """
    point = Box(0., 0., 0., 0.)
    start, end = trim((0., 0.), (1., 0.), point, point, .05)
    assert start == approx((.05, 0.)), 'public trim() gated its clearance: `trim_ink` leaked downward'
    assert end == approx((.95, 0.)), 'and it spends the same at the far end'


def test_trim_refuses_when_the_boxes_leave_nothing():
    """two crowded labels: a zero-length stroke with a round cap is a dot, which reads as a radical"""
    assert trim((0., 0.), (.6, 0.), Box(-.3, -.3, .3, .3), Box(.3, -.3, .9, .3), 0.) is None
    # and the control: the same bond with the two boxes just clear of each other still draws
    assert trim((0., 0.), (.6, 0.), Box(-.2, -.2, .2, .2), Box(.4, -.2, .8, .2), 0.) is not None


def test_the_clearance_is_dropped_before_the_bond_is():
    """rung 2: clearance is a PREFERENCE, so it is spent only if it fits

    The boxes below leave .1 of real space over a 1.0 bond and .06 of clearance at each end asks for
    .12, so rung 1 refuses and rung 2 draws from box edge to box edge.
    """
    p, q = (0., 0.), (1., 0.)
    box_p, box_q = Box(-.45, -.2, .45, .2), Box(.55, -.2, 1.45, .2)
    assert trim_per_end(p, q, box_p, box_q, .06, .06) is None, 'the premise: rung 1 refuses'
    log = []
    start, end = trim_ink(p, q, box_p, box_q, .06, log=log, atoms=(7, 9))
    assert start == approx((.45, 0.)), 'rung 2 stops at the label box, not short of it and not in it'
    assert end == approx((.55, 0.))
    assert [(r.rule, r.atoms) for r in log] == [('depict:tight', (7, 9))], 'drawn tighter, undisclosed'


def test_a_genuinely_overlapping_pair_of_labels_still_refuses():
    """rung 2 is not a licence to draw through a glyph: overlapping boxes are still no room"""
    log = []
    assert trim_ink((0., 0.), (.6, 0.), Box(-.4, -.3, .4, .3), Box(.2, -.3, 1., .3), .06,
                    log=log, atoms=(1, 2)) is None
    assert log == [], 'nothing was drawn, so nothing was drawn tightly'


def test_a_roomy_bond_never_reaches_the_second_rung():
    """the control: on an ordinary bond rung 1 answers, to the last bit, and nothing is logged"""
    p, q = (0., 0.), (1., 0.)
    box_p, box_q = Box(-.2, -.15, .2, .15), Box(.8, -.15, 1.2, .15)
    log = []
    assert trim_ink(p, q, box_p, box_q, .06, log=log, atoms=(1, 2)) == \
        trim_per_end(p, q, box_p, box_q, .06, .06)
    assert log == []


def test_hydrogen_peroxide_gets_a_bond_at_the_acs_preset():
    """`OO` at `acs`: .0084 of real space over a .825 bond, and 2 x .07 of clearance asked for

    Measured: the two label boxes leave head + tail = .81664 of a .82503 bond, so the ink does not
    overlap and the clearance alone would drop the only bond hydrogen peroxide has.
    """
    style = DepictStyle.preset('acs')
    mol = smiles('OO')
    mol.clean2d()
    plane = mol.coordinates()
    log = []
    paths = bond_paths(mol, plane, labels(mol, plane, style), style, log=log)
    assert sum(len(sub) - 1 for p in paths for sub in p.subpaths) == 1, 'the bond was not drawn'
    assert [(r.rule, r.atoms) for r in log] == [('depict:tight', (1, 2))]


def test_two_overlapping_labels_still_lose_their_bond_and_say_crowded():
    """`[NH3]~[BH3]` at `acs`, which rung 2 CANNOT save, and must not pretend to

    Measured at that preset: length .82503, head .42753, tail .47530 -- the ink boxes overlap before any
    clearance, so nothing is drawn and the id is `depict:crowded` rather than `depict:tight`.
    """
    style = DepictStyle.preset('acs')
    mol = smiles('[NH3]~[BH3]')
    mol.clean2d()
    plane = mol.coordinates()
    log = []
    paths = bond_paths(mol, plane, labels(mol, plane, style), style, log=log)
    assert sum(len(sub) - 1 for p in paths for sub in p.subpaths) == 0, 'drawn through a glyph'
    assert [(r.rule, r.atoms) for r in log] == [('depict:crowded', (1, 2))]


def test_an_ordinary_picture_never_reaches_the_second_rung():
    """the whole-picture form of the guard above: no `depict:tight` means rung 1 answered everywhere

    Every rung-2 firing that reaches a stroke is logged, so the absence of the record pins the geometry
    without transcribing a coordinate that goes stale when a font metric moves.
    """
    for text in ('CCO', 'c1ccccc1', 'CC(=O)Oc1ccccc1C(=O)O', 'OS(=O)(=O)O', 'C[N+](C)(C)C'):
        style = DepictStyle.preset('acs')
        mol = smiles(text)
        mol.clean2d()
        plane = mol.coordinates()
        log = []
        bond_paths(mol, plane, labels(mol, plane, style), style, log=log)
        assert not [r for r in log if r.rule == 'depict:tight'], f'{text}: rung 2 became the policy'


def test_trim_refuses_a_zero_length_bond():
    """two atoms on one point: a direction cannot be derived, so there is nothing to draw"""
    point = Box(0., 0., 0., 0.)
    assert trim((1., 1.), (1., 1.), point, point, 0.) is None


def test_trim_is_symmetric_in_its_arguments():
    a, b = Box(-.2, -.2, .2, .2), Box(.8, -.3, 1.2, .3)
    forward = trim((0., 0.), (1., 0.), a, b, .02)
    backward = trim((1., 0.), (0., 0.), b, a, .02)
    assert forward[0] == approx(backward[1])
    assert forward[1] == approx(backward[0])


def test_trim_moves_along_the_bond_and_not_along_an_axis():
    """a diagonal bond: both ends move ALONG it, so the trimmed segment is still collinear"""
    box = Box(-.2, -.2, .2, .2)
    start, end = trim((0., 0.), (3., 4.), box, Box(2.8, 3.8, 3.2, 4.2), .1)
    assert (end[0] - start[0]) * 4. == approx((end[1] - start[1]) * 3.), 'off the bond axis'
    assert hypot(start[0], start[1]) > .2, 'the start did not clear the box'


def test_no_output_contains_a_mask_or_a_white_knockout():
    """the trim is the mechanism; nothing is painted over anything

    Structural rather than a search for `#ffffff`, which nothing here can emit anyway: what a
    reintroduced mask changes is the stroke LENGTH, so the drawn end must fall short of the oxygen's
    centre by at least the half-width of its own label box.
    """
    mol, plane, boxes, style = _setup('CCO')
    paths = bond_paths(mol, plane, boxes, style)
    assert all(getattr(p, 'fill', None) is None and p.stroke != '#ffffff' for p in paths)

    oxygen = next(n for n in plane if boxes[n].text is not None)
    ox, oy = plane[oxygen]
    reach = min(hypot(x - ox, y - oy) for p in paths for sub in p.subpaths
                for x, y in _points(sub))
    box = boxes[oxygen].box
    assert reach >= (box.max_x - box.min_x) / 2, \
        'a full-length line under a mask would reach the atom centre; a trimmed one stops at the box'


def test_no_bond_is_trimmed_at_a_bare_vertex():
    """`bond.trim` is clearance from a label's INK, so a bare vertex takes none -- at either end

    "Nearest ink to this atom" is satisfied by ANY stroke, so only two of the five discriminate:
    `C#CC#C` for the bond AXIS, which only a triple bond exposes, and `CCCCCC` for `_strokes`' head/tail
    pair -- both 0 here and .06 with the trim ungated.  NOT `C=C=C` or a terminal `=CH2`, which read .09
    at head: that is `_straddle`'s `spacing / 2` and not a trim at all.  The positive control at the end
    is the mirror defect -- a stroke ending inside a glyph -- which a gap-only test cannot see.
    """
    for text in ('CC(C)C', 'c1ccccc1', 'CC=CC', 'C#CC#C', 'CCCCCC'):
        mol, plane, boxes, style = _setup(text)
        paths = bond_paths(mol, plane, boxes, style)
        drawn = [point for p in paths for sub in p.subpaths for point in _points(sub)]
        bare = [n for n in plane if boxes[n].text is None]
        assert bare, text
        for n in bare:
            x, y = plane[n]
            assert min(hypot(px - x, py - y) for px, py in drawn) < 1e-9, \
                f'{text}: the drawing stops short of bare atom {n}'

    # the control: a LABELLED atom is still cleared, by its own box and by `trim` on top of it
    mol, plane, boxes, style = _setup('CCO')
    paths = bond_paths(mol, plane, boxes, style)
    oxygen = next(n for n in plane if boxes[n].text is not None)
    ox, oy = plane[oxygen]
    reach = min(hypot(x - ox, y - oy) for p in paths for sub in p.subpaths for x, y in _points(sub))
    assert reach >= boxes[oxygen].box.width / 2 + style.bond.trim, \
        'a bond ends inside the oxygen glyph: the clearance was gated where the ink is'


def test_a_chain_of_unlabelled_atoms_is_one_path():
    """hexane: five bonds, ONE path, four miter joins -- the disconnected-bond fix"""
    mol, plane, boxes, style = _setup('CCCCCC')
    paths = bond_paths(mol, plane, boxes, style)
    assert len(paths) == 1
    assert len(paths[0].subpaths) == 1
    assert len(paths[0].subpaths[0]) == 6, 'M plus five L'
    assert paths[0].join == 'miter'


def test_a_chain_passes_exactly_through_its_interior_vertices():
    """the join is AT the vertex: a trimmed interior end would reopen the notch the chain closes"""
    mol, plane, boxes, style = _setup('CCCCCC')
    interior = _points(bond_paths(mol, plane, boxes, style)[0].subpaths[0])[1:-1]
    assert len(interior) == 4
    for x, y in interior:
        assert any(hypot(px - x, py - y) < 1e-9 for px, py in plane.values()), \
            'an interior chain point is not an atom point'


def test_a_label_breaks_a_chain():
    """the oxygen's box interrupts the stroke, so the two halves cannot be one continuous path"""
    mol, plane, boxes, style = _setup('CCCOCCC')
    paths = bond_paths(mol, plane, boxes, style)
    assert sum(len(p.subpaths) for p in paths) >= 2
    assert all(len(sub) >= 2 for p in paths for sub in p.subpaths)


def test_a_branch_point_breaks_a_chain():
    """three bonds cannot be one stroke through a degree-3 atom without drawing a bond twice"""
    mol, plane, boxes, style = _setup('CC(C)C')
    paths = bond_paths(mol, plane, boxes, style)
    total_segments = sum(len(sub) - 1 for p in paths for sub in p.subpaths)
    assert total_segments == 3, f'three bonds, {total_segments} segments drawn'


def test_a_ring_is_chained_and_closed():
    mol, plane, boxes, style = _setup('C1CCCCC1')
    paths = bond_paths(mol, plane, boxes, style)
    assert sum(len(sub) - 1 for p in paths for sub in p.subpaths) == 6
    assert any(sub[-1][0] == 'Z' for p in paths for sub in p.subpaths), \
        'a closed ring closes the path rather than repeating its first point'


def test_a_closed_ring_does_not_repeat_its_first_point():
    """`Z` is the closure; an M...L back to the start would stack two caps where a join belongs"""
    mol, plane, boxes, style = _setup('C1CCCCC1')
    sub, = (s for p in bond_paths(mol, plane, boxes, style) for s in p.subpaths)
    points = _points(sub)
    assert len(points) == 6, 'six vertices, no repeat'
    assert hypot(points[0][0] - points[-1][0], points[0][1] - points[-1][1]) > 1e-6


def test_every_bond_is_drawn_exactly_once():
    """chaining must not double a bond or drop one -- and a dropped bond is invisible in a picture

    This asks for the kekule notation, so an aromatic ring adds one inner line per double of the
    alternating pattern, `len(ring) // 2` in an isolated ring.  Tetrahydropyran is the one member that
    sees the WRAP-AROUND bond: a ring interrupted by exactly one label reaches `_strokes` as a cycle
    whose first and last atom are the same, and without the repeat it draws 5 of its 6 bonds while every
    other assertion here stays true.
    """
    for text in ('CCCCCC', 'CC(C)C', 'C1CCCCC1', 'C1CCOCC1', 'c1ccccc1C(=O)O',
                 'CC(=O)OC1=CC=CC=C1C(=O)O'):
        mol, plane, boxes, style = _setup(text, KEKULE)
        drawn = sum(len(sub) - 1 for p in bond_paths(mol, plane, boxes, style) for sub in p.subpaths)
        expected = sum(1 for _ in mol.bonds())
        # a double bond draws two lines, a triple three
        expected += sum(bond.order - 1 for bond in mol.bonds() if bond.order in (2, 3))
        expected += sum(len(ring) // 2 for ring in mol.aromatic_rings)
        assert drawn == expected, text

    # under `'circle'` the whole aromatic perimeter is a single stroke, so every heteroaromatic ring has
    # the shape above.  The circle is four cubics and no bond, so it is counted out by its `C` commands.
    style = DepictStyle().tuned(**{'bond.aromatic': 'circle'})
    mol, plane, boxes, style = _setup('c1ccncc1', style)
    drawn = sum(len(sub) - 1 for p in bond_paths(mol, plane, boxes, style) for sub in p.subpaths
                if not any(segment[0] == 'C' for segment in sub))
    assert drawn == 6, 'pyridine under circle lost a ring bond: the perimeter did not wrap around'


def test_a_double_bond_is_two_lines_in_one_path():
    mol, plane, boxes, style = _setup('C=C')
    paths = bond_paths(mol, plane, boxes, style)
    assert len(paths) == 1
    assert len(paths[0].subpaths) == 2, 'ONE path: one element, one paint, one join'


def test_the_two_lines_of_a_double_bond_are_the_style_spacing_apart():
    """the acyclic straddle AND the ring inset: the two lines are `spacing` apart in both notations

    A ring double takes `inner_line` at `spacing` on the centroid's side rather than `_straddle` at
    `spacing / 2` either side of the axis, so measuring only ethene cannot tell the two apart.  For the
    inset the quantity measured is the PERPENDICULAR distance from the axis.
    """
    style = DepictStyle().tuned(**{'bond.spacing': .2})
    mol, plane, boxes, style = _setup('C=C', style)
    (first, second), = (p.subpaths for p in bond_paths(mol, plane, boxes, style))
    separation = hypot(first[0][1] - second[0][1], first[0][2] - second[0][2])
    assert separation == approx(.2, abs=1e-9)

    mol, plane, boxes, style = _setup('C1=CC=CC=C1', style)
    seen = 0
    for path in bond_paths(mol, plane, boxes, style):
        if len(path.subpaths) != 2:
            continue
        seen += 1
        axis, inner = path.subpaths
        ax, ay, bx, by = axis[0][1], axis[0][2], axis[1][1], axis[1][2]
        length = hypot(bx - ax, by - ay)
        for point in _points(inner):
            across = ((point[0] - ax) * (by - ay) - (point[1] - ay) * (bx - ax)) / length
            assert abs(across) == approx(.2, abs=1e-9), 'the inner line is not `spacing` off the axis'
    assert seen == 3, f'three ring doubles in kekule benzene, found {seen}'


def test_the_two_lines_of_a_double_bond_are_the_same_length():
    """one trim on the axis, applied to both lines: unequal lengths are what per-line trimming gives"""
    mol, plane, boxes, style = _setup('CC=O')
    for path in bond_paths(mol, plane, boxes, style):
        if len(path.subpaths) != 2:
            continue
        lengths = [hypot(s[1][1] - s[0][1], s[1][2] - s[0][2]) for s in path.subpaths]
        assert lengths[0] == approx(lengths[1], abs=1e-9)


def test_a_double_bond_in_a_ring_puts_its_second_line_inside():
    """outside would read as a different ring; every chemist's eye expects the inner line"""
    mol, plane, boxes, style = _setup('C1=CC=CC=C1')
    centre_x = sum(x for x, _ in plane.values()) / len(plane)
    centre_y = sum(y for _, y in plane.values()) / len(plane)
    seen = 0
    for path in bond_paths(mol, plane, boxes, style):
        if len(path.subpaths) != 2:
            continue
        seen += 1
        first, second = path.subpaths
        d1 = hypot((first[0][1] + first[1][1]) / 2 - centre_x, (first[0][2] + first[1][2]) / 2 - centre_y)
        d2 = hypot((second[0][1] + second[1][1]) / 2 - centre_x,
                   (second[0][2] + second[1][2]) / 2 - centre_y)
        assert min(d1, d2) < max(d1, d2), 'the two lines are at different radii'
    assert seen == 3, f'kekule benzene has three ring doubles, found {seen}'


def test_a_terminal_double_bond_is_drawn_symmetrically():
    """a carbonyl's two lines straddle the bond axis; there is no ring to pick a side from"""
    mol, plane, boxes, style = _setup('CC=O')
    seen = 0
    for path in bond_paths(mol, plane, boxes, style):
        if len(path.subpaths) != 2:
            continue
        seen += 1
        first, second = path.subpaths
        # both offsets from the axis, equal and opposite
        assert isclose(hypot(first[0][1] - second[0][1], first[0][2] - second[0][2]),
                       style.bond.spacing, abs_tol=1e-9)
        # equal and opposite, not both on one side: the mean of the two lines is ON the bond axis
        (n, m), = [(b.n, b.m) for b in mol.bonds() if b.order == 2]
        (px, py), (qx, qy) = plane[n], plane[m]
        mid_x = (first[0][1] + second[0][1]) / 2 - px
        mid_y = (first[0][2] + second[0][2]) / 2 - py
        assert mid_x * (qy - py) - mid_y * (qx - px) == approx(0., abs=1e-9), 'the mean is off the axis'
    assert seen == 1, 'the one double bond'


def test_a_triple_bond_is_three_lines_with_the_centre_on_the_axis():
    mol, plane, boxes, style = _setup('CC#N')
    triple = [p for p in bond_paths(mol, plane, boxes, style) if len(p.subpaths) == 3]
    assert len(triple) == 1
    (n, m), = [(b.n, b.m) for b in mol.bonds() if b.order == 3]
    (px, py), (qx, qy) = plane[n], plane[m]
    on_axis = [s for s in triple[0].subpaths
               if abs((s[0][1] - px) * (qy - py) - (s[0][2] - py) * (qx - px)) < 1e-9]
    assert len(on_axis) == 1, 'exactly one of the three lines is the axis itself'
    offsets = sorted(round((s[0][1] - px) * (qy - py) - (s[0][2] - py) * (qx - px), 9)
                     / hypot(qx - px, qy - py) for s in triple[0].subpaths)
    assert offsets[0] == approx(-offsets[2], abs=1e-9), 'the outer lines straddle the axis'
    assert offsets[2] - offsets[0] == approx(2. * style.bond.triple_spacing, abs=1e-9)


def test_an_aromatic_ring_is_dashed_inner_by_default():
    """the reference look: a solid perimeter with a dashed line inside it, and no alternating second lines

    Benzene is TWO paths under the default -- one closed six-point perimeter run and one closed dashed
    arc -- and the count is asserted so a regression to per-bond subpaths, which anchors a dash to every
    corner, fails rather than merely looking worse.
    """
    mol, plane, boxes, style = _setup('c1ccccc1')
    paths = bond_paths(mol, plane, boxes, style)
    assert len(paths) == 2, f'one perimeter, one dashed arc -- got {len(paths)}'
    dashed, = [p for p in paths if p.dashes]
    # The notation is (.15, .05); what is drawn is that pattern compensated for the DEFAULT ROUND CAP,
    # which lengthens each dash by half the stroke width at each end.  See the compensation test below.
    painted, gap = style.bond.aromatic_dashes
    assert dashed.dashes == approx((painted - style.bond.width, gap + style.bond.width))
    assert len(dashed.subpaths) == 1 and dashed.subpaths[0][-1] == ('Z',), 'one closed run'
    assert not any(len(p.subpaths) == 2 for p in paths), 'no alternating second lines'
    assert not any(seg[0] == 'C' for p in paths for sub in p.subpaths for seg in sub), 'not the circle'


def test_the_kekule_notation_is_still_available():
    """the control for the test above: the alternating lines are a notation, not a deleted feature"""
    mol, plane, boxes, style = _setup('c1ccccc1', KEKULE)
    paths = bond_paths(mol, plane, boxes, style)
    assert sum(1 for p in paths if len(p.subpaths) == 2) == 3, 'three alternating lines'
    assert not any(p.dashes for p in paths)


def test_a_kekulized_molecule_draws_no_inner_ring_under_the_default():
    """CORRECT, and pinned so nobody "fixes" it: `aromatic_rings` filters on order 4

    A Kekule structure has no order-4 bond, so it has no aromatic ring to ornament -- and it needs
    none, because the alternating lines it already carries say the same thing.  The aromatic spelling of
    the same compound is beside it so the assertion cannot be passing for want of a ring.
    """
    mol, plane, boxes, style = _setup('C1=CC=CC=C1')
    assert mol.aromatic_rings == [], 'premise: a kekulized ring is not an aromatic ring'
    paths = bond_paths(mol, plane, boxes, style)
    assert not any(p.dashes for p in paths), 'an inner line on a Kekule ring says nothing new'
    assert sum(1 for p in paths if len(p.subpaths) == 2) == 3, 'and the three doubles are drawn'

    mol, plane, boxes, style = _setup('c1ccccc1')
    assert any(p.dashes for p in bond_paths(mol, plane, boxes, style)), \
        'the control: the aromatic spelling of the same ring DOES get the arc'


def test_kekule_drawing_does_not_mutate_the_molecule():
    """a picture is not a repair pass: the aromatic bonds are still aromatic afterwards"""
    mol, plane, boxes, style = _setup('c1ccccc1C(=O)O', KEKULE)
    before = [(b.n, b.m, b.order) for b in mol.bonds()]
    bond_paths(mol, plane, boxes, style)
    assert [(b.n, b.m, b.order) for b in mol.bonds()] == before


def test_kekule_drawing_leaks_no_repair_back_into_the_molecule():
    """the alternation is `kekule()`'s answer, taken on a COPY -- and `kekule()` REPAIRS

    It may move a charge or derive a hydrogen count, so the snapshot is every order, charge, implicit
    hydrogen count and radical flag.  Imidazole and pyrrole are where a repair is available to leak.
    """
    for text in ('c1cnc[nH]1', 'c1cc[nH]c1', 'Cn1cnc2c1c(=O)n(C)c(=O)n2C'):
        mol, plane, boxes, style = _setup(text, KEKULE)
        before = (sorted((b.n, b.m, b.order) for b in mol.bonds()),
                  {a.n: (a.charge, a.implicit_h, a.is_radical) for a in mol.atoms()})
        bond_paths(mol, plane, boxes, style)
        after = (sorted((b.n, b.m, b.order) for b in mol.bonds()),
                 {a.n: (a.charge, a.implicit_h, a.is_radical) for a in mol.atoms()})
        assert after == before, f'{text}: drawing it changed it'


def test_a_heteroaromatic_never_doubles_a_bond_its_heteroatom_cannot_carry():
    """no second line on a heteroatom that cannot carry one

    Each of these four has exactly ONE Kekule form, so `expected` is chemistry and not an arbitrary
    choice among forms.  The `wrong` column is the bond a maximal matching draws -- graph theory with no
    element, charge or valence in it -- asserted absent so a regression to it fails by name.
    """
    cases = ((('c1ccoc1'), {(1, 5), (2, 3)}, {(3, 4)}),                       # furan, O4
             (('c1cc[nH]c1'), {(1, 5), (2, 3)}, {(3, 4)}),                    # pyrrole, N4-H
             (('c1cnc[nH]1'), {(1, 2), (3, 4)}, {(1, 5)}),                    # imidazole, N5-H
             ('Cn1cnc2c1c(=O)n(C)c(=O)n2C',                                   # caffeine
              {(3, 4), (5, 6)}, {(2, 6), (4, 5), (7, 9), (11, 13)}))
    for text, expected, wrong in cases:
        doubles = _doubles(smiles(text))
        assert doubles == expected, f'{text}: {sorted(doubles)}, expected {sorted(expected)}'
        assert not (doubles & wrong), f'{text}: the greedy matching\'s bond is back'


def test_the_drawn_furan_puts_no_second_line_on_its_oxygen():
    """the same invariant through the DRAWN OUTPUT, not through the helper, in case the call site errs"""
    mol, plane, boxes, style = _setup('c1ccoc1', KEKULE)
    oxygen, = [a.n for a in mol.atoms() if a.atomic_symbol == 'O']
    for path in bond_paths(mol, plane, boxes, style):
        if len(path.subpaths) != 2:
            continue
        for x, y in _points(path.subpaths[0]):
            near = min(plane, key=lambda n: hypot(plane[n][0] - x, plane[n][1] - y))
            assert near != oxygen, 'a double line ends on furan\'s oxygen'


def test_kekule_keeps_the_counts_the_matching_got_right():
    """benzene three, naphthalene five: the COUNT, not which bonds

    Benzene has two Kekule forms and naphthalene three, all correct, so pinning a chosen set would pin
    an arbitrary choice; the count and the no-atom-twice property are what a chemist checks.
    """
    for text, count in (('c1ccccc1', 3), ('c1ccc2ccccc2c1', 5)):
        doubles = _doubles(smiles(text))
        assert len(doubles) == count, f'{text}: {len(doubles)} doubles, expected {count}'
        seen = [n for key in doubles for n in key]
        assert len(seen) == len(set(seen)), f'{text}: an atom carries two double lines'


def test_the_kekule_alternation_is_deterministic():
    """two pictures of one molecule must not differ: `doubles` is the same set both times"""
    mol = smiles('c1ccc2[nH]ccc2c1')
    assert _doubles(mol) == _doubles(mol)


def test_an_aromatic_system_with_no_kekule_form_gets_the_circle_and_a_log_line():
    """the fallback, and its premise verified rather than assumed

    A circle says "aromatic", which is the whole of what is known; five plain single lines would say
    cyclopentane.  `log` is the caller's only signal that the picture is weaker than the default.
    """
    mol = smiles('c1cccc1')
    assert mol.copy().kekule().unresolved, 'the premise: this system really has no Kekule form'

    mol, plane, boxes, style = _setup('c1cccc1', KEKULE)
    log = []
    paths = bond_paths(mol, plane, boxes, style, log=log)
    assert any(any(seg[0] == 'C' for seg in sub) for p in paths for sub in p.subpaths), \
        'the aromatic circle is cubics, and it is not there'
    assert not any(len(p.subpaths) == 2 and all(len(s) == 2 for s in p.subpaths) for p in paths), \
        'an unresolved system gets no second lines at all'
    records = [r for r in log if r.rule == 'depict:no-kekule']
    assert len(records) == 1, f'one record per unresolved system, found {len(records)}'
    assert set(records[0].atoms) == set(mol), 'the record names the atoms of the system'


def test_an_ordinary_aromatic_molecule_logs_no_kekule_failure():
    """the control: a log asserted only for presence is satisfied by one that fires on everything"""
    for text in ('c1ccccc1', 'c1ccoc1', 'c1ccc2ccccc2c1', 'Cn1cnc2c1c(=O)n(C)c(=O)n2C'):
        mol, plane, boxes, style = _setup(text, KEKULE)
        log = []
        bond_paths(mol, plane, boxes, style, log=log)
        assert not [r for r in log if r.rule == 'depict:no-kekule'], f'{text}: a false failure'


def test_kekule_never_gives_one_atom_two_double_lines():
    """the alternation is a matching: two inner lines meeting at one atom is a valence a chemist reads

    Naphthalene catches a per-ring alternation with no memory: its two rings share a bond, so a second
    ring alternating from scratch can double up on the shared atoms.
    """
    for text in ('c1ccc2ccccc2c1', 'c1ccc2cc3ccccc3cc2c1', 'c1ccc2[nH]ccc2c1'):
        mol, plane, boxes, style = _setup(text, KEKULE)
        seen = {}
        for path in bond_paths(mol, plane, boxes, style):
            if len(path.subpaths) != 2:
                continue
            # the axis subpath's two ends are the two atoms of the bond it belongs to
            for x, y in _points(path.subpaths[0]):
                near = min(plane, key=lambda n: hypot(plane[n][0] - x, plane[n][1] - y))
                seen[near] = seen.get(near, 0) + 1
        assert seen and max(seen.values()) == 1, f'{text}: an atom carries two double lines'


def test_an_aromatic_ring_can_be_drawn_as_a_circle():
    style = DepictStyle().tuned(**{'bond.aromatic': 'circle'})
    mol, plane, boxes, style = _setup('c1ccccc1', style)
    paths = bond_paths(mol, plane, boxes, style)
    assert any(any(seg[0] == 'C' for seg in sub) for p in paths for sub in p.subpaths), \
        'the inner circle is cubics'
    assert not any(len(p.subpaths) == 2 and all(len(s) == 2 for s in p.subpaths) for p in paths), \
        'a circle replaces the alternating lines, it does not join them'


def test_an_aromatic_perimeter_is_chained_under_circle():
    """benzene under `circle`: the six ring bonds are ONE closed stroke, plus the circle

    The chain predicate is "the notation draws exactly one plain stroke along this axis", not "order 1",
    and under `circle` that holds of every aromatic bond.  The count is asserted so a regression to
    per-bond paths -- two butt caps meeting at 120 degrees at every corner -- fails rather than looking
    merely worse.
    """
    style = DepictStyle().tuned(**{'bond.aromatic': 'circle'})
    mol, plane, boxes, style = _setup('c1ccccc1', style)
    paths = bond_paths(mol, plane, boxes, style)
    assert len(paths) == 2, 'one perimeter, one circle'
    perimeter, = [p for p in paths if not any(seg[0] == 'C' for sub in p.subpaths for seg in sub)]
    assert len(perimeter.subpaths) == 1
    assert len(_points(perimeter.subpaths[0])) == 6
    assert perimeter.subpaths[0][-1] == ('Z',), 'an unbroken ring closes'


def test_a_fused_aromatic_perimeter_is_chained_under_circle_and_dashed_inner():
    """naphthalene: five paths -- three chained runs and two ornaments

    Three runs and not two closed rings: the fusion carbons are degree-3 branch points, so each ring
    contributes an open five-bond run and the fusion bond is a path of its own.
    """
    for mode in ('circle', 'dashed-inner'):
        style = DepictStyle().tuned(**{'bond.aromatic': mode})
        mol, plane, boxes, style = _setup('c1ccc2ccccc2c1', style)
        paths = bond_paths(mol, plane, boxes, style)
        assert len(paths) == 5, f'{mode}: three chained runs and two ornaments, got {len(paths)}'
        strokes = [p for p in paths if not p.dashes
                   and not any(seg[0] == 'C' for sub in p.subpaths for seg in sub)]
        assert sorted(len(_points(sub)) for p in strokes for sub in p.subpaths) == [2, 6, 6]
        drawn = sum(len(sub) - 1 for p in strokes for sub in p.subpaths)
        assert drawn == 11, f'{mode}: eleven bonds, {drawn} segments'


def test_a_hetero_ring_perimeter_does_not_close_through_its_label():
    """pyridine under `circle`, and tetrahydropyran: the walk comes back, the stroke still may not close

    Closing a cycle whose joining atom carries a glyph draws the stroke straight through it, so it is one
    OPEN run trimmed at both ends against the same box.
    """
    for text, mode in (('c1ccncc1', 'circle'), ('C1CCOCC1', 'kekule')):
        style = DepictStyle().tuned(**{'bond.aromatic': mode})
        mol, plane, boxes, style = _setup(text, style)
        paths = bond_paths(mol, plane, boxes, style)
        strokes = [p for p in paths if not any(seg[0] == 'C' for sub in p.subpaths for seg in sub)]
        assert not any(seg == ('Z',) for p in strokes for sub in p.subpaths for seg in sub), \
            f'{text}: the perimeter closed through the heteroatom label'
        hetero = next(n for n in plane if boxes[n].text is not None)
        hx, hy = plane[hetero]
        reach = min(hypot(x - hx, y - hy) for p in strokes for sub in p.subpaths
                    for x, y in _points(sub))
        assert reach >= boxes[hetero].box.width / 2, f'{text}: a stroke ends inside the glyph'


def test_kekule_does_not_chain_the_bonds_it_gave_a_second_line():
    """the discriminating case for the generalized predicate: widen it to "order 4" and this fails

    A predicate of "order 4 always" swallows the doubled bonds' axes into the perimeter run and their
    second lines vanish -- toluene loses three lines and still looks like a ring.
    """
    mol, plane, boxes, style = _setup('Cc1ccccc1', KEKULE)
    paths = bond_paths(mol, plane, boxes, style)
    doubles = [p for p in paths if len(p.subpaths) == 2]
    assert len(doubles) == 3, f'three ring doubles, found {len(doubles)}'
    for path in doubles:
        assert all(len(_points(sub)) == 2 for sub in path.subpaths), 'a doubled bond was chained'
    drawn = sum(len(sub) - 1 for p in paths for sub in p.subpaths)
    assert drawn == 10, f'seven bonds and three second lines, {drawn} segments drawn'


def test_the_aromatic_circle_sits_inside_the_ring_by_the_style_inset():
    """`aromatic_inset` is the gap from the bonds, so the radius is the apothem less the inset"""
    style = DepictStyle().tuned(**{'bond.aromatic': 'circle', 'bond.aromatic_inset': .3})
    mol, plane, boxes, style = _setup('c1ccccc1', style)
    circle, = [p for p in bond_paths(mol, plane, boxes, style)
               if any(seg[0] == 'C' for sub in p.subpaths for seg in sub)]
    ring, = mol.aromatic_rings
    cx, cy = _centroid(plane, ring)
    radii = {round(hypot(x - cx, y - cy), 6) for sub in circle.subpaths for x, y in _points(sub)}
    apothem = min(hypot((plane[a][0] + plane[b][0]) / 2 - cx, (plane[a][1] + plane[b][1]) / 2 - cy)
                  for a, b in zip(ring, ring[1:] + ring[:1]))
    assert len(radii) == 1, 'a circle has one radius'
    assert radii.pop() == approx(apothem - .3, abs=1e-6)   # the set rounds to six places


def test_an_aromatic_ring_can_be_drawn_with_a_dashed_inner_arc():
    style = DepictStyle().tuned(**{'bond.aromatic': 'dashed-inner'})
    mol, plane, boxes, style = _setup('c1ccccc1', style)
    assert any(p.dashes for p in bond_paths(mol, plane, boxes, style))


def test_the_dashed_inner_arc_reads_its_own_inset_and_not_the_circles():
    """TWO insets: the dashed arc reads its own, not the circle's

    Each field is tuned ALONE, nothing else in the suite telling them apart, and the molecule and its
    layout are built ONCE so a second `clean2d()` cannot pose as a moved arc.
    """
    mol = smiles('c1ccccc1')
    mol.clean2d()
    plane = mol.coordinates()
    ring, = mol.aromatic_rings
    cx, cy = _centroid(plane, ring)

    def arc_radius(style):
        dashed, = [p for p in bond_paths(mol, plane, labels(mol, plane, style), style) if p.dashes]
        first = dashed.subpaths[0]
        return hypot((first[0][1] + first[1][1]) / 2 - cx, (first[0][2] + first[1][2]) / 2 - cy)

    apothem = min(hypot((plane[a][0] + plane[b][0]) / 2 - cx, (plane[a][1] + plane[b][1]) / 2 - cy)
                  for a, b in zip(ring, ring[1:] + ring[:1]))
    base = DepictStyle().tuned(**{'bond.aromatic': 'dashed-inner'})
    # 1e-4 on the two that compare against the APOTHEM, and nothing tighter is available: the apothem is
    # a `min()` over a `clean2d()` hexagon whose six differ by 3.655e-05, so each residual is 3.5e-05 --
    # the irregularity itself, ~3x below this tolerance and 2600x below the .26 the third assertion needs.
    assert arc_radius(base) == approx(apothem - base.bond.aromatic_dash_inset, abs=1e-4)
    assert arc_radius(base.tuned(**{'bond.aromatic_dash_inset': .3})) == approx(apothem - .3, abs=1e-4)
    # the control needs no tolerance: one layout answers both sides, so the claim is that the number does
    # not change at all rather than by less than the hexagon's own error.
    assert arc_radius(base.tuned(**{'bond.aromatic_inset': .4})) == approx(arc_radius(base), abs=1e-12), \
        "the dashed arc moved with the CIRCLE's inset: the two fields are crossed"


def test_the_dashed_inner_ring_of_benzene_is_one_closed_subpath():
    """SUBPATHS, not elements: SVG restarts the dash phase and the join at every M

    Six two-point subpaths would anchor a dash to every corner and stack two caps where a join belongs,
    while `len(paths) == 1` reports that as correct.
    """
    style = DepictStyle().tuned(**{'bond.aromatic': 'dashed-inner'})
    mol, plane, boxes, style = _setup('c1ccccc1', style)
    dashed, = [p for p in bond_paths(mol, plane, boxes, style) if p.dashes]
    assert len(dashed.subpaths) == 1, 'one ring, one continuous run'
    assert len(_points(dashed.subpaths[0])) == 6, 'six corners'
    assert dashed.subpaths[0][-1] == ('Z',), 'an unbroken ring closes'


def test_a_label_breaks_the_dashed_inner_ring_open():
    """pyridine: the trim at N moves both ends beside it, so the run cannot close through the glyph"""
    style = DepictStyle().tuned(**{'bond.aromatic': 'dashed-inner'})
    mol, plane, boxes, style = _setup('c1ccncc1', style)
    dashed, = [p for p in bond_paths(mol, plane, boxes, style) if p.dashes]
    assert not any(seg == ('Z',) for sub in dashed.subpaths for seg in sub), \
        'the ring closed through the nitrogen label'
    assert len(dashed.subpaths) == 1, 'still ONE open run all the way round, not six'
    assert len(_points(dashed.subpaths[0])) == 7, 'six lines end to end, open: one point more'


def test_an_inner_line_is_inside_its_own_ring_and_the_forced_side_is_not():
    """the invariant, over a fused and heteroaromatic corpus, with its control beside it

    `inner_line` returns the two POINTS and picks the side itself so a caller cannot drop the sign.  The
    assertion is point-in-polygon against the bond's own ring; the control is the same geometry with the
    side forced positive, and were it to score as well the assertion would be measuring nothing.
    """
    style = DepictStyle()
    total = inside = control_inside = 0
    for text in RINGS:
        mol, plane, boxes, style = _setup(text, style)
        for ring in mol.aromatic_rings:
            centroid = _centroid(plane, ring)
            polygon = [plane[n] for n in ring]
            for a, b in zip(ring, ring[1:] + ring[:1]):
                total += 1
                line = inner_line(plane[a], plane[b], centroid, style.bond.aromatic_dash_inset)
                assert line is not None, f'{text}: a regular ring bond has an inner line'
                if all(_inside(point, polygon) for point in line):
                    inside += 1
                control = _forced_plus_inner_line(plane[a], plane[b], centroid,
                                                  style.bond.aromatic_dash_inset)
                if control is not None and all(_inside(point, polygon) for point in control):
                    control_inside += 1
    assert total > 60, f'the corpus shrank: {total} ring bonds'
    assert inside == total, f'{inside}/{total} inner lines are inside their ring'
    # .6 and not merely "< total": the control is the dropped sign, so roughly half of every ring's bonds
    # must come out on the wrong side.  `< total` would pass at 98/99.  Measured here: 30/99.
    assert control_inside < total * .6, \
        f'the control scored {control_inside}/{total}: too close to the real thing to be measuring it'


def test_a_ring_double_bond_draws_its_inner_line_inside_the_ring():
    """the same invariant one level up: `bond_paths` must not undo the side `inner_line` chose"""
    style = KEKULE
    total = inside = 0
    for text in RINGS:
        mol, plane, boxes, style = _setup(text, style)
        rings = [(r, _centroid(plane, r), [plane[n] for n in r]) for r in mol.aromatic_rings]
        for path in bond_paths(mol, plane, boxes, style):
            if len(path.subpaths) != 2:
                continue
            total += 1
            axis, inner = path.subpaths
            mid = ((inner[0][1] + inner[1][1]) / 2, (inner[0][2] + inner[1][2]) / 2)
            axis_mid = ((axis[0][1] + axis[1][1]) / 2, (axis[0][2] + axis[1][2]) / 2)
            # the ring this bond belongs to: the one whose polygon the bond's own midpoint sits on
            best = min(rings, key=lambda r: hypot(r[1][0] - axis_mid[0], r[1][1] - axis_mid[1]))
            if _inside(mid, best[2]):
                inside += 1
    assert total > 25, f'too few ring doubles to measure: {total}'
    assert inside == total, f'{inside}/{total} ring doubles put their second line inside the ring'


def test_the_inner_lines_of_a_ring_meet_at_its_corners():
    """the shortening comes off the vertex bisector, so consecutive inner lines close with no spur

    A fixed shrink cannot do this -- the right amount depends on the vertex angle, so a constant that
    closes a hexagon leaves a gap in a five-ring, which is why furan is in the list.  And the amount is
    NOT `bond.trim`, which is clearance from a label's ink.
    """
    style = DepictStyle()
    worst = 0.
    for text in ('c1ccccc1', 'c1ccncc1', 'c1ccc2ccccc2c1', 'c1ccoc1'):
        mol, plane, boxes, style = _setup(text, style)
        for ring in mol.aromatic_rings:
            centroid = _centroid(plane, ring)
            lines = [inner_line(plane[a], plane[b], centroid, style.bond.aromatic_dash_inset)
                     for a, b in zip(ring, ring[1:] + ring[:1])]
            for first, second in zip(lines, lines[1:] + lines[:1]):
                worst = max(worst, hypot(first[1][0] - second[0][0], first[1][1] - second[0][1]))
    assert worst < 1e-3, f'worst corner gap {worst:g}: the runs cannot merge into one subpath'


def test_an_inner_line_refuses_a_centroid_collinear_with_the_bond():
    """the inset explodes there, and a line through the wrong ring is worse than no line"""
    assert inner_line((0., 0.), (1., 0.), (.5, 0.), .14) is None, 'exactly collinear'
    assert inner_line((0., 0.), (1., 0.), (.5, .2), .14) is None, 'offset/|perp| >= .65'
    assert inner_line((0., 0.), (1., 0.), (.5, .25), .14) is not None, 'and just inside the limit'


def test_an_inner_line_refuses_when_the_bisectors_leave_nothing():
    """a skewed ring projects both bisector crossings past the bond's own footprint

    Both land behind the start, so the clamp leaves the two ends coincident; without it the "inner" line
    would sit beside a bond it does not belong to.
    """
    assert inner_line((0., 0.), (1., 0.), (-2., .3), .14) is None


def test_an_inner_line_is_clamped_to_the_bonds_footprint():
    """the clamp, where it still leaves something: no end may overhang the bond it insets"""
    line = inner_line((0., 0.), (1., 0.), (-.4, .35), .14)
    assert line is not None
    assert line[0][0] == approx(0.), 'the near end was clamped to the start'
    assert 0. <= line[1][0] <= 1.


def test_an_inner_line_stays_on_its_bond_however_distorted_the_ring():
    """the clamp, fuzzed: whatever the layout does, an inner line lies ON the bond it insets

    In bond-axis coordinates the invariant is exact for every sample; an axis-aligned bounding box holds
    for 17% of them and point-in-polygon for 44%, so both would pin noise.  Refusing is a legal answer
    for a nearly collinear centroid, hence the closing assertion that most bonds are still drawn.
    """
    from math import isfinite
    from random import Random

    rnd = Random(7)
    inset = DepictStyle().bond.aromatic_dash_inset
    total = drawn = 0
    for text in ('c1ccccc1', 'c1ccc2ccccc2c1', 'c1cc[nH]c1', 'c1ccncc1'):
        mol = smiles(text)
        mol.clean2d()
        base = {a.n: (a.x, a.y) for a in mol.atoms()}
        rings = [tuple(r) for r in mol.aromatic_rings]
        for _ in range(500):
            plane = {n: (x + rnd.uniform(-1.6, 1.6), y + rnd.uniform(-1.6, 1.6))
                     for n, (x, y) in base.items()}
            for ring in rings:
                centre = (sum(plane[n][0] for n in ring) / len(ring),
                          sum(plane[n][1] for n in ring) / len(ring))
                for a, b in zip(ring, ring[1:] + ring[:1]):
                    total += 1
                    line = inner_line(plane[a], plane[b], centre, inset)
                    if line is None:
                        continue
                    drawn += 1
                    (ax, ay), (bx, by) = plane[a], plane[b]
                    dx, dy = bx - ax, by - ay
                    d2 = dx * dx + dy * dy
                    for px, py in line:
                        assert isfinite(px) and isfinite(py)
                        t = ((px - ax) * dx + (py - ay) * dy) / d2
                        assert -1e-9 <= t <= 1. + 1e-9, f'{text}: inner line overhangs its bond at t={t}'
    assert drawn > total * .6, f'only {drawn}/{total} drawn: the test would pass by refusing everything'


def test_a_dative_bond_is_dashed_headless_and_distinguishable_from_a_single_bond():
    """the ammonia-boron trifluoride adduct, written NEUTRAL: dashed, headless, not a single bond

    `~` is the dative bond in SMILES and without it this is an ordinary single bond, so the order-8 guard
    is an assertion of its own.  No head: the container stores order 8 and not the arrow's direction, so
    a head at either end would be picked out of atom order.  The B-F bonds are NOT dashed, which is what
    makes the notation legible.
    """
    mol = smiles('[NH3]~[B](F)(F)F')
    mol.clean2d()
    plane = mol.coordinates()
    style = DepictStyle()
    assert any(b.order == 8 for b in mol.bonds()), 'the SMILES has to carry a dative bond'
    paths = bond_paths(mol, plane, labels(mol, plane, style), style)
    assert all(p.fill is None for p in paths), 'no head: the direction is not a fact to draw'
    dashed = [p for p in paths if p.dashes]
    assert len(dashed) == 1, 'one dative bond, one dashed line'
    painted, gap = style.bond.dative_dashes                  # compensated for the default round cap
    assert dashed[0].dashes == approx((painted - style.bond.width, gap + style.bond.width))
    assert len(dashed[0].subpaths) == 1 and len(_points(dashed[0].subpaths[0])) == 2, 'a plain line'
    solid = [p for p in paths if not p.dashes]
    assert len(solid) == 3, 'and the three B-F single bonds are solid, so the reader can tell them apart'


def test_a_crowded_bond_is_skipped_and_logged():
    """two labels whose boxes overlap: no bond is drawn, and the caller is told which"""
    mol = smiles('OS(=O)(=O)O')
    mol.clean2d()
    plane = mol.coordinates()
    style = DepictStyle().tuned(**{'label.size': 2.5, 'label.pad': .5})   # absurd, on purpose
    log = []
    bond_paths(mol, plane, labels(mol, plane, style), style, log=log)
    assert log, 'a bond was dropped and nothing said so'
    assert all(record.rule.startswith('depict:') for record in log)
    assert len(log[0].atoms) == 2


def test_a_roomy_picture_logs_nothing():
    """the control for the record above: an ordinary molecule must not log at all"""
    mol, plane, boxes, style = _setup('CCO')
    log = []
    bond_paths(mol, plane, boxes, style, log=log)
    assert log == []


def test_the_line_width_comes_from_the_style():
    mol, plane, boxes, style = _setup('CCO', DepictStyle().tuned(**{'bond.width': .077}))
    assert all(p.width == approx(.077) for p in bond_paths(mol, plane, boxes, style) if p.stroke)


def test_bond_colour_comes_from_the_style():
    mol, plane, boxes, style = _setup('CCO', DepictStyle().tuned(**{'bond.colour': '#883300'}))
    assert all(p.stroke == '#883300' for p in bond_paths(mol, plane, boxes, style) if p.stroke)


def test_every_bond_is_round_capped_by_default():
    """a bond is drawn as SEVERAL paths, and their ends have to meet

    Ethanol is one chain run and one trimmed stub either side of the O; a ring adds an inner line, a
    double bond a second line, a stereo centre a wedge.  Each stops at its own end, and two butt caps
    arriving at one point from two angles show the notch between them, so every end is the same disc.
    The default is asserted rather than only the plumbing, `bond.cap` being what a figure is drawn with.
    """
    style = DepictStyle()
    assert style.bond.cap == 'round', 'the default itself, not merely that it reaches the path'
    for text in ('CCO', 'c1ccccc1O', 'CC(=O)Nc1ccccc1'):
        mol, plane, boxes, style = _setup(text)
        paths = bond_paths(mol, plane, boxes, style)
        assert paths, text
        assert all(p.cap == 'round' for p in paths if p.stroke), text


def test_the_cap_join_and_miter_limit_come_from_the_style():
    style = DepictStyle().tuned(**{'bond.cap': 'round', 'bond.join': 'bevel', 'bond.miter_limit': 2.})
    mol, plane, boxes, style = _setup('CCCCCC', style)
    path, = bond_paths(mol, plane, boxes, style)
    assert (path.cap, path.join, path.miter_limit) == ('round', 'bevel', 2.)


def test_a_skipped_bond_is_not_drawn():
    """wedges draw their own bond, so `skip` has to remove it from here"""
    mol, plane, boxes, style = _setup('CCCCCC')
    keys = [(b.n, b.m) if b.n < b.m else (b.m, b.n) for b in mol.bonds()]
    paths = bond_paths(mol, plane, boxes, style, skip={keys[2]})
    drawn = sum(len(sub) - 1 for p in paths for sub in p.subpaths)
    assert drawn == 4, f'five bonds less one skipped, {drawn} drawn'
    assert len(paths) == 2, 'and the chain is cut in two where the bond went'


def test_a_per_bond_width_overrides_the_style_and_breaks_the_chain():
    """a single path cannot taper, so a bond drawn wider is a path of its own"""
    mol, plane, boxes, style = _setup('CCCCCC')
    keys = [(b.n, b.m) if b.n < b.m else (b.m, b.n) for b in mol.bonds()]
    paths = bond_paths(mol, plane, boxes, style, widths={keys[2]: .2})
    assert sum(len(sub) - 1 for p in paths for sub in p.subpaths) == 5
    assert sorted(p.width for p in paths) == [approx(style.bond.width), approx(style.bond.width),
                                              approx(.2)]


def test_a_per_bond_colour_overrides_the_style_and_breaks_the_chain():
    mol, plane, boxes, style = _setup('CCCCCC')
    keys = [(b.n, b.m) if b.n < b.m else (b.m, b.n) for b in mol.bonds()]
    paths = bond_paths(mol, plane, boxes, style, colours={keys[0]: '#ff0000'})
    assert sum(len(sub) - 1 for p in paths for sub in p.subpaths) == 5
    assert sorted(p.stroke for p in paths) == ['#000000', '#ff0000']


def test_an_override_is_read_low_first_whichever_way_the_caller_wrote_it():
    """the key convention, stated once: `(low, high)`, and the reverse is not a second bond"""
    mol, plane, boxes, style = _setup('CCCCCC')
    keys = [(b.n, b.m) if b.n < b.m else (b.m, b.n) for b in mol.bonds()]
    reversed_key = (keys[2][1], keys[2][0])
    paths = bond_paths(mol, plane, boxes, style, widths={reversed_key: .2})
    assert any(p.width == approx(.2) for p in paths), 'a reversed key found no bond'


def test_chains_never_include_a_labelled_atom_in_the_middle():
    mol, plane, boxes, style = _setup('CCOCC')
    oxygen = [a.n for a in mol.atoms() if a.atomic_symbol == 'O'][0]
    found = False
    for chain in chains(mol, boxes):
        assert oxygen not in chain[1:-1], 'a label interrupts the stroke'
        found = found or oxygen in chain
    assert found, 'the oxygen still ends two chains'


def test_chains_cover_every_single_bond_exactly_once():
    """the chaining invariant at its own level, so a failure names the walker and not the drawing"""
    for text in ('CCCCCC', 'CC(C)C', 'C1CCCCC1', 'CCOCC', 'CC(=O)OC1=CC=CC=C1C(=O)O', 'CCC.CCC'):
        mol, plane, boxes, style = _setup(text)
        singles = {(b.n, b.m) if b.n < b.m else (b.m, b.n) for b in mol.bonds() if b.order == 1}
        walked = []
        for chain in chains(mol, boxes):
            assert len(chain) >= 2, 'a one-atom chain draws nothing'
            for a, b in zip(chain, chain[1:]):
                walked.append((a, b) if a < b else (b, a))
        assert sorted(walked) == sorted(singles), text


def test_chains_pass_through_a_bare_degree_two_vertex_and_stop_at_a_branch():
    mol, plane, boxes, style = _setup('CC(C)CC')
    branch = [a.n for a in mol.atoms() if a.degree == 3][0]
    for chain in chains(mol, boxes):
        assert branch not in chain[1:-1], 'a branch point cannot be passed through'
    assert max(len(c) for c in chains(mol, boxes)) == 3, 'the two-bond run through the bare CH2'


def test_a_chain_does_not_pass_through_a_multiple_bond():
    """a double bond is its own shape, so the run stops at it rather than drawing it as a line"""
    mol, plane, boxes, style = _setup('CCC=CCC')
    doubled = {a for b in mol.bonds() if b.order == 2 for a in (b.n, b.m)}
    for chain in chains(mol, boxes):
        for a, b in zip(chain, chain[1:]):
            assert not (a in doubled and b in doubled), 'the double bond was chained as a line'
        assert not (set(chain[1:-1]) & doubled), 'a run passed through a double-bonded atom'


def test_chains_of_a_disconnected_molecule_do_not_bridge_the_components():
    mol, plane, boxes, style = _setup('CCC.CCC')
    components = [set(c) for c in mol.connected_components]
    assert len(components) == 2
    for chain in chains(mol, boxes):
        assert any(set(chain) <= component for component in components)


def test_a_closed_chain_reports_its_closure():
    """cyclohexane: the walk came back, and it says so by repeating the first atom last"""
    mol, plane, boxes, style = _setup('C1CCCCC1')
    chain, = chains(mol, boxes)
    assert chain[0] == chain[-1]
    assert len(set(chain)) == 6


def test_the_dash_pattern_is_compensated_for_a_round_cap():
    """a round cap extends every dash by half the stroke width at each end

    With a butt cap there is nothing to compensate, so the configured lengths go through unchanged: one
    number in the style, two honest renderings of it.  Round is the default, so the butt case is the one
    that has to say so.
    """
    style = DepictStyle().tuned(**{'bond.aromatic': 'dashed-inner', 'bond.width': .04,
                                   'bond.aromatic_dashes': (.15, .05)})
    mol, plane, boxes, style = _setup('c1ccccc1', style)
    round_capped, = [p.dashes for p in bond_paths(mol, plane, boxes, style) if p.dashes]
    assert round_capped == approx((.15 - .04, .05 + .04))

    style = style.tuned(**{'bond.cap': 'butt'})
    butt, = [p.dashes for p in bond_paths(mol, plane, boxes, style) if p.dashes]
    assert butt == approx((.15, .05))


def test_a_compensated_dash_never_goes_non_positive():
    """`Path` refuses a non-positive dash length, so the floor is load-bearing"""
    style = DepictStyle().tuned(**{'bond.aromatic': 'dashed-inner', 'bond.cap': 'round',
                                   'bond.width': .3, 'bond.aromatic_dashes': (.15, .05)})
    mol, plane, boxes, style = _setup('c1ccccc1', style)
    dashes, = [p.dashes for p in bond_paths(mol, plane, boxes, style) if p.dashes]
    assert dashes[0] > 0.


def test_bond_paths_needs_a_box_for_every_atom():
    """`labels()` returns one per atom; a caller passing a filtered dict gets told, not a KeyError"""
    mol, plane, boxes, style = _setup('CCO')
    with raises(ValueError):
        bond_paths(mol, plane, {n: b for n, b in boxes.items() if n != 1}, style)


def test_bond_paths_needs_a_point_for_every_atom():
    """the same guard on the OTHER mapping, which would otherwise be a KeyError from inside a walk"""
    mol, plane, boxes, style = _setup('CCO')
    with raises(ValueError, match='coordinates'):
        bond_paths(mol, {n: p for n, p in plane.items() if n != 1}, boxes, style)
