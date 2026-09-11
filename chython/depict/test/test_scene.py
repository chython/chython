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
"""The scene IR: three primitives, molecule coordinates, y-up.

Three backends draw the same picture and two of them (PDF, PostScript) have no <mask>, no CSS and no
reusable <defs>, so geometry is objects rather than SVG: a test asserts on a curve, not on a string.
"""
from math import isclose
from pytest import approx, mark, raises
from chython.depict.scene import (BLACK, EMPTY_BOX, Box, Group, Path, Scene, Text, TextRun, circle,
                                  close, curve, ellipse, line, move, polyline, rounded_box, rgb, to_hex)


def test_a_path_holds_subpaths_and_is_immutable():
    p = Path([[move(0., 0.), line(1., 1.)]], stroke=BLACK, width=.06)
    assert p.subpaths == (((('M', 0., 0.), ('L', 1., 1.))),)
    assert p.stroke == '#000000'
    with raises(AttributeError):
        p.width = .1


def test_a_stroked_path_needs_a_width():
    with raises(ValueError, match='width'):
        Path([[move(0., 0.), line(1., 0.)]], stroke=BLACK)


def test_a_path_that_neither_fills_nor_strokes_is_refused():
    """an invisible path is a bug in the caller, not a thing to serialize"""
    with raises(ValueError, match='fill or stroke'):
        Path([[move(0., 0.), line(1., 0.)]])


def test_an_empty_subpath_is_refused():
    with raises(ValueError, match='empty'):
        Path([[]], fill=BLACK)


def test_a_subpath_must_start_with_a_move():
    with raises(ValueError, match='M'):
        Path([[line(1., 1.)]], fill=BLACK)


def test_bounds_of_a_line_include_the_stroke_width():
    p = Path([[move(0., 0.), line(1., 0.)]], stroke=BLACK, width=.1)
    assert p.bounds == approx((-.05, -.05, 1.05, .05))


def test_bounds_of_a_curve_use_the_control_hull():
    """conservative on purpose: a hull never clips, and a tight bezier bound is arithmetic no viewBox
    needs -- the cost of being generous is whitespace, the cost of being wrong is a cropped picture"""
    p = Path([[move(0., 0.), curve(0., 2., 1., 2., 1., 0.)]], fill=BLACK)
    assert p.bounds == approx((0., 0., 1., 2.))


def test_a_fill_only_path_is_not_inflated():
    p = Path([[move(0., 0.), line(1., 1.), close()]], fill=BLACK)
    assert p.bounds == approx((0., 0., 1., 1.))


def test_circle_is_four_cubics_and_closes():
    segments = circle(1., 2., .5)
    assert segments[0][0] == 'M'
    assert [s[0] for s in segments] == ['M', 'C', 'C', 'C', 'C', 'Z']
    box = Path([segments], fill=BLACK).bounds
    assert box == approx((.5, 1.5, 1.5, 2.5)), 'the control hull of a k=0.5523 circle is its bbox'


def test_circle_passes_through_its_four_cardinal_points():
    segments = circle(0., 0., 1.)
    ends = [(s[-2], s[-1]) for s in segments if s[0] in 'MC']
    for point in ((1., 0.), (0., 1.), (-1., 0.), (0., -1.)):
        assert any(isclose(x, point[0], abs_tol=1e-12) and isclose(y, point[1], abs_tol=1e-12)
                   for x, y in ends), point


def test_a_circle_is_an_ellipse_with_one_radius():
    """ONE arc approximation in the package: `circle` delegates, so the two can never disagree"""
    assert circle(1., 2., .5) == ellipse(1., 2., .5, .5)


def test_an_ellipse_passes_through_its_four_cardinal_points():
    """KAPPA is per axis -- the quarter arcs are independent in x and y -- so a wide ellipse is not a
    scaled circle's control points with a scaled radius."""
    segments = ellipse(0., 0., 2., .5)
    ends = [(s[-2], s[-1]) for s in segments if s[0] in 'MC']
    for point in ((2., 0.), (0., .5), (-2., 0.), (0., -.5)):
        assert any(isclose(x, point[0], abs_tol=1e-12) and isclose(y, point[1], abs_tol=1e-12)
                   for x, y in ends), point
    assert Path([segments], fill=BLACK).bounds == approx((-2., -.5, 2., .5))


def test_a_rounded_box_is_four_lines_and_four_corner_arcs_within_its_box():
    segments = rounded_box(Box(0., 0., 4., 2.), .5)
    assert [s[0] for s in segments] == ['M', 'L', 'C', 'L', 'C', 'L', 'C', 'L', 'C', 'Z']
    assert Path([segments], fill=BLACK).bounds == approx((0., 0., 4., 2.)), \
        'the corner control points are ON the box, so a plate never spills past what it was measured from'


def test_a_rounded_boxs_radius_is_clamped_to_the_stadium():
    """the caller is a plate behind a one-digit number, where the radius asked for exceeds the box

    The shape wanted at that limit is the stadium -- fully round ends -- and not an exception for having
    asked for more than the box has room for, so the radius is clamped and the box is still respected.
    """
    stadium = rounded_box(Box(0., 0., 4., 2.), 1.)
    assert rounded_box(Box(0., 0., 4., 2.), 50.) == stadium
    assert Path([stadium], fill=BLACK).bounds == approx((0., 0., 4., 2.))


def test_a_rounded_box_with_no_radius_is_the_rectangle():
    """no arcs at all, rather than four degenerate cubics a renderer has to collapse"""
    assert rounded_box(Box(0., 0., 4., 2.), 0.) == polyline(((4., 0.), (4., 2.), (0., 2.), (0., 0.)),
                                                            closed=True)


def test_polyline_open_and_closed():
    assert polyline([(0., 0.), (1., 0.)]) == (('M', 0., 0.), ('L', 1., 0.))
    assert polyline([(0., 0.), (1., 0.)], closed=True) == (('M', 0., 0.), ('L', 1., 0.), ('Z',))


def test_text_carries_runs_so_a_label_is_one_object():
    """CH3 is one anchored label with three runs, not three <text> elements a caller has to place

    A run carries its own dy/dx, so a subscript is a property of the run and the label is placed once.
    """
    t = Text([TextRun('CH', size=.4), TextRun('3', size=.28, dy=-.12)], x=1., y=2., anchor='middle')
    assert len(t.runs) == 2
    assert t.anchor == 'middle'
    assert t.runs[1].dy == approx(-.12)


def test_text_refuses_an_unknown_anchor():
    with raises(ValueError, match='anchor'):
        Text([TextRun('C', size=.4)], x=0., y=0., anchor='centre')


def test_text_refuses_no_runs():
    with raises(ValueError, match='run'):
        Text([], x=0., y=0.)


def test_group_opacity_is_the_overlap_mechanism():
    """group opacity, not per-child alpha: two overlapping bands inside one group composite ONCE

    This is what removes the polygon boolean union from the contour code -- filled bands drawn in
    painter's order inside a 0.35-opacity group look like one translucent shape.
    """
    g = Group([Path([[move(0., 0.), line(1., 0.)]], stroke=BLACK, width=.1)], opacity=.35)
    assert g.opacity == approx(.35)
    assert g.bounds == approx((-.05, -.05, 1.05, .05))


def test_group_opacity_is_range_checked():
    with raises(ValueError, match='opacity'):
        Group([], opacity=1.5)


def test_an_empty_group_has_empty_bounds_and_does_not_poison_a_union():
    scene = Scene([Group([]), Path([[move(0., 0.), line(1., 1.)]], fill=BLACK)])
    assert scene.bounds == approx((0., 0., 1., 1.))


def test_scene_bounds_is_the_union_of_its_children():
    scene = Scene([Path([[move(0., 0.), line(1., 0.)]], fill=BLACK),
                   Path([[move(2., -1.), line(2., 3.)]], fill=BLACK)])
    assert scene.bounds == approx((0., -1., 2., 3.))


def test_scene_bounds_can_be_stated_and_then_wins():
    """a caller that wants a fixed frame -- a grid cell, an animation -- states it and is obeyed"""
    scene = Scene([Path([[move(0., 0.), line(1., 0.)]], fill=BLACK)], bounds=Box(-1., -1., 5., 5.))
    assert scene.bounds == approx((-1., -1., 5., 5.))


def test_an_empty_scene_has_a_degenerate_box_and_not_a_crash():
    assert Scene([]).bounds == approx((0., 0., 0., 0.))


def test_frame_inflates_computed_bounds_but_not_a_stated_one():
    """backends call `scene.frame(margin)` to get the render box; the asymmetry is the contract"""
    computed = Scene([Path([[move(0., 0.), line(1., 0.)]], fill=BLACK)])
    assert computed.frame(.5) == approx((-.5, -.5, 1.5, .5)), 'computed bounds are inflated by margin'

    stated = Scene([Path([[move(0., 0.), line(1., 0.)]], fill=BLACK)], bounds=Box(0., 0., 2., 1.))
    assert stated.frame(.5) == approx((0., 0., 2., 1.)), 'stated bounds are returned as given'


def test_box_helpers():
    box = Box(0., 0., 2., 1.)
    assert (box.width, box.height) == approx((2., 1.))
    assert box.inflate(.5) == approx((-.5, -.5, 2.5, 1.5))
    assert box.union(Box(-1., 0., 1., 1.)) == approx((-1., 0., 2., 1.))
    assert box.translated(1., -2.) == approx((1., -2., 3., -1.))
    assert box.contains(Box(.5, .1, 1.5, .9))
    assert not box.contains(Box(.5, .1, 2.5, .9))


def test_box_of_an_empty_iterable_is_the_empty_box_and_not_the_origin():
    """the bug this exists to prevent: seeding a union with Box(0, 0, 0, 0) puts the origin in every box"""
    assert Box.of([]) is EMPTY_BOX
    assert Box.of([Box(3., 3., 4., 4.)]) == approx((3., 3., 4., 4.))
    assert Box.of([Box(3., 3., 4., 4.), Box(-1., 0., 0., 1.)]) == approx((-1., 0., 4., 4.))


def test_the_empty_box_is_contained_by_anything_so_an_empty_group_never_fails_a_check():
    assert Box(0., 0., 1., 1.).contains(EMPTY_BOX)


def test_translating_a_path_moves_every_coordinate_and_changes_nothing_else():
    # The control points are at y=-1, so the untranslated hull is (0., -1., 1., 0.), pinned by the last
    # assertion below; a curve bulging the other way satisfies neither it nor the translated bounds.
    p = Path([[move(0., 0.), curve(0., -1., 1., -1., 1., 0.), close()]], fill=BLACK, even_odd=True)
    moved = p.translated(2., -1.)
    assert moved.bounds == approx((2., -2., 3., -1.))
    assert moved.fill == p.fill and moved.even_odd == p.even_odd
    assert p.bounds == approx((0., -1., 1., 0.)), 'and the original is untouched'


def test_translating_a_text_moves_its_anchor():
    t = Text([TextRun('N')], x=1., y=2.)
    assert t.translated(-1., 1.).x == approx(0.)
    assert t.translated(-1., 1.).y == approx(3.)


def test_translating_a_group_moves_its_children_and_its_clip():
    clip = Path([[move(0., 0.), line(1., 0.), line(1., 1.), close()]], fill=BLACK)
    g = Group([Path([[move(0., 0.), line(1., 1.)]], stroke=BLACK, width=.1)], clip=clip)
    moved = g.translated(5., 5.)
    assert moved.children[0].bounds.min_x == approx(4.95)
    assert moved.clip.bounds == approx((5., 5., 6., 6.))


def test_colours_are_hex_strings_and_the_two_names_resolve():
    assert rgb(255, 0, 128) == '#ff0080'
    assert to_hex('black') == '#000000'
    assert to_hex('white') == '#ffffff'
    assert to_hex('#ABCDEF') == '#abcdef'
    assert to_hex(None) is None


def test_a_three_digit_colour_expands_to_six():
    """`#abc` is the CSS short form and means `#aabbcc` -- doubling each digit, not padding with zeros

    The wrong expansion (`#0a0b0c`) is legal, dark, and silently not the colour the caller asked for.
    """
    assert to_hex('#abc') == '#aabbcc'
    assert to_hex('#FFF') == '#ffffff'
    assert to_hex('#000') == '#000000'


def test_an_unparsable_colour_is_refused_where_it_is_written():
    """not at draw time, three modules later, in a backend that cannot say which path it came from"""
    with raises(ValueError, match='colour'):
        Path([[move(0., 0.), line(1., 0.)]], fill='cornflowerblue')


def test_dashes_are_a_tuple_of_positive_lengths():
    p = Path([[move(0., 0.), line(1., 0.)]], stroke=BLACK, width=.05, dashes=(.1, .05))
    assert p.dashes == (.1, .05)
    with raises(ValueError, match='dash'):
        Path([[move(0., 0.), line(1., 0.)]], stroke=BLACK, width=.05, dashes=(.1, -.05))


def test_to_svg_is_reachable_from_the_scene():
    """serialization hangs off the scene, and scene.py must not import render/ at module level or the two
    would cycle.  All three doors, each being its own lazy import -- `_repr_svg_` is the one Jupyter calls
    unasked, where a `NameError` surfaces only as a cell that prints nothing."""
    scene = Scene([Path([[move(0., 0.), line(1., 0.)]], stroke=BLACK, width=.1)])
    assert scene.to_svg().startswith('<svg')
    assert scene.to_svgz()[:2] == b'\x1f\x8b'
    assert scene._repr_svg_().startswith('<svg')


def test_the_print_backends_say_so_until_they_exist():
    """and each message names the release AND the thing to do instead

    Matching on `PDF`/`EPS` alone is satisfied by "lands with the print formats", which leaves a caller
    waiting; the redirect to `to_svg()` is therefore asserted too.
    """
    for call, fmt in ((Scene([]).to_pdf, 'PDF'), (Scene([]).to_eps, 'EPS')):
        with raises(NotImplementedError) as e:
            call()
        assert fmt in str(e.value)
        assert '3.1' in str(e.value), 'the message must name the release, not an unnamed future'
        assert 'to_svg()' in str(e.value), 'and must name what to do instead'


# Every construction this IR refuses, and the fragment of the message naming the reason.  An unreached
# `raise` is indistinguishable from one whose condition is inverted or whose field name is misspelled.
_REFUSALS = [
    (lambda: rgb(256, 0, 0), 'out of range', 'a channel above 255'),
    (lambda: rgb(0, -1, 0), 'out of range', 'a negative channel'),
    (lambda: to_hex(0x00ff00), 'string', 'an int is not a colour'),
    (lambda: to_hex('#gg0000'), 'unparsable', 'seven characters that are not hex'),
    (lambda: to_hex('#xyz'), 'unparsable', 'four characters that are not hex'),
    (lambda: polyline([]), 'at least one point', 'a polyline through nothing'),
    (lambda: circle(0., 0., 0.), 'radius', 'a circle with no radius'),
    (lambda: ellipse(0., 0., 1., 0.), 'radii', 'an ellipse flat in one axis'),
    (lambda: rounded_box(Box(0., 0., 1., 1.), -.1), 'negative', 'a negative corner radius'),
    (lambda: rounded_box(EMPTY_BOX, .1), 'empty box', 'the corners of nothing'),
    (lambda: Path([[move(0., 0.), ('Q', 1., 1.)]], fill=BLACK), 'unknown segment', 'a quadratic'),
    (lambda: Path([], fill=BLACK), 'at least one subpath', 'a path with no subpaths'),
    (lambda: Path([[move(0., 0.), line(1., 0.)]], stroke=BLACK, width=0.), 'width', 'a zero-width stroke'),
    (lambda: Path([[move(0., 0.), line(1., 0.)]], fill=BLACK, cap='flat'), 'cap', 'an unknown cap'),
    (lambda: Path([[move(0., 0.), line(1., 0.)]], fill=BLACK, join='mitre'), 'join', 'an unknown join'),
    (lambda: Path([[move(0., 0.), line(1., 0.)]], stroke=BLACK, width=.05, miter_limit=.5),
     'miter limit', 'a miter limit below the ratio of 1 it is'),
    (lambda: TextRun(''), 'draws nothing', 'a run with no text'),
    (lambda: TextRun('C', size=0.), 'size', 'a run with no size'),
    (lambda: TextRun('C', weight='light'), 'weight', 'a weight no base-14 font has'),
    (lambda: TextRun('C', style='oblique'), 'style', 'a style no base-14 font has'),
]


@mark.parametrize('build,match', [(b, m) for b, m, _ in _REFUSALS], ids=[i for _, _, i in _REFUSALS])
def test_the_ir_refuses_it_where_it_is_written(build, match):
    with raises(ValueError, match=match):
        build()


def test_scene_bounds_includes_a_text_child():
    """the commonest real scene is a label over a bond, and its bounds must include the label"""
    scene = Scene([Path([[move(0., 0.), line(1., 0.)]], stroke=BLACK, width=.1),
                   Text([TextRun('CH', size=.4)], x=2., y=2.)])
    assert scene.bounds.max_x > 1.05
