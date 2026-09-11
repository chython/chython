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
"""The legend, and the two things it must never do: lie about resolution, or overlap the structure.

Discrete swatches rather than a gradient, one per band, so a reader can match them one to one.  Placement
is measured: the tick labels are typeset through `metrics`, the bar's box INCLUDES them, and
`place_colorbar` puts that whole box outside the content box.
"""
from pytest import approx, raises
from chython.depict.colorbar import colorbar, legend_side, place_colorbar
from chython.depict.colormap import Colormap
from chython.depict.field import contour_levels
from chython.depict.overlay import Swatch, tiled_swatches
from chython.depict.scene import Box, Path, Text
from chython.depict.style import DepictStyle


def _swatches_of(nodes):
    return [n for n in nodes if isinstance(n, Path) and n.fill is not None and n.stroke is None]


def _rules_of(nodes):
    """A swatch with no interval of values behind it draws as a rule: a stroke, no fill."""
    return [n for n in nodes if isinstance(n, Path) and n.fill is None and n.stroke is not None]


def _adjacent(pairs, cmap):
    """`(value, colour)` pairs as swatches, each covering up to the next one.

    For the tests whose subject is the COLOURS and their order rather than the intervals.
    """
    ordered = sorted(pairs)
    return [Swatch(v, c, v, ordered[i + 1][0] if i + 1 < len(ordered) else cmap.vmax, True)
            for i, (v, c) in enumerate(ordered)]


def _bar(cmap=None, levels=5, swatches=None, side='right', length=4.):
    """`colorbar` takes the `Swatch`es the figure DREW and derives nothing itself, so the helper supplies
    them: `levels=` for the no-bands case (a plain scale tiled edge to edge), `swatches=` otherwise."""
    style = DepictStyle()
    cmap = cmap if cmap is not None else Colormap.named('coolwarm').fitted([-.4, .6])
    if swatches is None:
        swatches = tiled_swatches(cmap, contour_levels(levels, cmap.vmin, cmap.vmax))
    return colorbar(cmap, swatches=swatches, style=style, side=side, length=length)


def test_a_bar_has_one_swatch_per_level():
    nodes, _ = _bar(levels=7)
    assert len(_swatches_of(nodes)) == 7


def test_each_swatch_is_the_colour_it_was_handed_and_not_one_recomputed():
    """the caller passes the colours the bands were FILLED with and the bar draws exactly those.  The four
    below are not `contour_levels` of anything, so no re-derivation from the domain could produce them."""
    cmap = Colormap.named('coolwarm').fitted([-.4, .6])
    given = [(-.35, '#5671d8'), (-.25, '#7493e7'), (-.15, '#94b1ee'), (-.05, '#b8c9f1')]
    nodes, _ = _bar(cmap=cmap, swatches=_adjacent(given, cmap))
    assert [s.fill for s in _swatches_of(nodes)] == [c for _, c in given]


def test_the_swatches_run_low_to_high_whatever_order_they_arrive_in():
    """`bands_of` is in PAINTER'S order -- outermost band first, which for a diverging field is the level
    nearest zero -- and a bar is read low to high, so the bar sorts"""
    cmap = Colormap.named('coolwarm').fitted([-.4, .6])
    given = [(-.05, '#b8c9f1'), (-.35, '#5671d8'), (-.15, '#94b1ee'), (-.25, '#7493e7')]
    nodes, _ = _bar(cmap=cmap, swatches=_adjacent(given, cmap), side='right', length=4.)
    swatches = _swatches_of(nodes)
    by_height = [s.fill for s in sorted(swatches, key=lambda s: s.bounds.min_y)]
    assert by_height == ['#5671d8', '#7493e7', '#94b1ee', '#b8c9f1']


def test_a_plain_scale_tiles_with_no_gap_and_no_overlap():
    """`tiled_swatches` is the no-bands case -- a halo's bar -- where nothing is there for a gap to mean,
    so the strip is continuous.  A FIELD's bar is allowed gaps: see below."""
    nodes, _ = _bar(levels=5, side='right', length=5.)
    swatches = _swatches_of(nodes)
    spans = sorted((s.bounds.min_y, s.bounds.max_y) for s in swatches)
    for (_, top), (bottom, _) in zip(spans, spans[1:]):
        assert top == approx(bottom), 'a gap or an overlap between two swatches reads as a defect'


def test_a_swatch_sits_at_its_own_VALUES_along_the_domain_and_not_in_the_nth_equal_block():
    """THE AXIS IS THE MAP'S DOMAIN: four bands crowded into the bottom third draw four swatches in the
    bottom third, not four equal blocks, so two figures on one `domain=` stay comparable"""
    cmap = Colormap.named('coolwarm').fitted([-.4, .6])
    assert (cmap.vmin, cmap.vmax) == (-.6, .6), 'a diverging map symmetrizes: the field is at one end'
    given = [(-.35, '#5671d8'), (-.25, '#7493e7'), (-.15, '#94b1ee'), (-.05, '#b8c9f1')]
    nodes, _ = _bar(cmap=cmap, swatches=_adjacent(given, cmap), side='right', length=6.)
    swatches = sorted(_swatches_of(nodes), key=lambda s: s.bounds.min_y)

    for swatch, (value, _) in zip(swatches, given):
        assert swatch.bounds.min_y == approx(cmap.normalised(value) * 6.), 'not at its ordinal place'
    # the last one runs to the top of the domain, because on the page it covers everything past its level
    assert swatches[-1].bounds.max_y == approx(6.)
    # and the four together do NOT fill the bar: -0.6..-0.35 is empty, and it is empty in the picture too
    assert swatches[0].bounds.min_y > .1 * 6.


def test_a_level_with_no_band_leaves_a_gap_where_it_belongs():
    """a field that traced no closed contour at one level must not have that level's stretch of the bar
    filled in by its neighbours -- a bands-range axis closes the gap up and shows a continuous strip"""
    cmap = Colormap.named('coolwarm').fitted([-.4, .4])
    drawn = [Swatch(-.3, '#5671d8', -.4, -.2, True), Swatch(.1, '#f4c6ad', .1, .4, True)]
    nodes, _ = _bar(cmap=cmap, swatches=drawn, side='right', length=8.)
    spans = sorted((s.bounds.min_y, s.bounds.max_y) for s in _swatches_of(nodes))
    assert len(spans) == 2
    assert spans[0][1] < spans[1][0], 'the two swatches must not touch'
    gap = (spans[0][1], spans[1][0])
    assert gap == approx((cmap.normalised(-.2) * 8., cmap.normalised(.1) * 8.))


def test_a_swatch_with_no_interval_is_a_rule_and_not_a_block_of_colour():
    """an open contour encloses no region, so no interval of values wears its colour and the bar marks the
    one value it does mean -- as for every level in `fill=False` mode"""
    cmap = Colormap.named('coolwarm').fitted([-.4, .4])
    nodes, _ = _bar(cmap=cmap, swatches=[Swatch(-.2, '#5671d8', -.2, -.2, False)], side='right', length=8.)
    assert _swatches_of(nodes) == [], 'a block of colour would claim an interval the field never filled'
    rule = _rules_of(nodes)
    assert len(rule) == 1 and rule[0].stroke == '#5671d8'
    box = rule[0].bounds
    # the bounds of a stroke include half its width, so the VALUE is the centre line
    assert (box.min_y + box.max_y) / 2. == approx(cmap.normalised(-.2) * 8.), 'and at its own value'
    assert box.height < DepictStyle().page.legend_breadth, 'across the strip, not along it'
    assert box.width > DepictStyle().page.legend_breadth


def test_a_vertical_bar_runs_the_length_it_was_given():
    _, box = _bar(side='right', length=4.)
    assert box.max_y - box.min_y >= 4.


def test_a_horizontal_bar_runs_along_x_instead():
    nodes, box = _bar(side='bottom', length=6.)
    assert box.max_x - box.min_x >= 6.
    assert box.max_y - box.min_y < 6.


def test_the_ticks_are_min_max_and_zero():
    """three ticks, and the exact strings -- the format is the deliverable, not an accident"""
    from chython.depict.label import MINUS

    cmap = Colormap.named('coolwarm').fitted([-.4, .6])
    assert (cmap.vmin, cmap.vmax) == (-.6, .6), 'a diverging map symmetrizes; see the note below'
    nodes, _ = _bar(cmap=cmap)
    written = sorted(''.join(r.text for r in n.runs) for n in nodes if isinstance(n, Text))
    assert written == ['+0.00', '+0.60', MINUS + '0.60'], written


def test_a_map_that_does_not_span_zero_gets_two_ticks():
    # viridis([.2, .9]) stays [.2, .9]: a sequential map does not symmetrize
    """no zero tick where zero is off the scale -- a tick outside the bar is worse than none"""
    nodes, _ = _bar(cmap=Colormap.named('viridis').fitted([.2, .9]))
    written = [n for n in nodes if isinstance(n, Text)]
    assert len(written) == 2


def test_a_negative_tick_is_written_with_a_typographic_minus():
    from chython.depict.label import MINUS

    nodes, _ = _bar(cmap=Colormap.named('coolwarm').fitted([-.4, .6]))
    written = ''.join(r.text for n in nodes if isinstance(n, Text) for r in n.runs)
    assert MINUS in written
    assert '-' not in written


def test_a_tick_label_never_sits_on_the_strip():
    """Both sides, because the two get the tick from different geometry.

    A `Text`'s y is its BASELINE, so a tick placed at the gap alone has its glyphs back inside the strip
    -- digits on the colours they are labelling, which is the one thing a legend may not do.
    """
    for side in ('right', 'bottom'):
        nodes, _ = _bar(side=side)
        strip = Box.of([n.bounds for n in nodes if isinstance(n, Path)])
        for tick in (n for n in nodes if isinstance(n, Text)):
            ink = tick.bounds
            clear = (ink.max_x <= strip.min_x or ink.min_x >= strip.max_x
                     or ink.max_y <= strip.min_y or ink.min_y >= strip.max_y)
            assert clear, f'{side} tick {ink} overlaps the swatch strip {strip}'


def test_a_vertical_bars_tick_is_centred_on_its_own_value():
    """Not baselined on it: the label would read half a cap height above the colour it names."""
    cmap = Colormap.named('viridis').fitted([.2, .9])
    nodes, _ = _bar(cmap=cmap, side='right', length=4.)
    ticks = sorted((n.bounds for n in nodes if isinstance(n, Text)),
                   key=lambda ink: ink.min_y)
    assert len(ticks) == 2
    assert (ticks[0].min_y + ticks[0].max_y) / 2. == approx(0.)      # vmin, at the strip's foot
    assert (ticks[1].min_y + ticks[1].max_y) / 2. == approx(4.)      # vmax, at its head


def test_the_box_includes_the_tick_labels_and_not_just_the_strip():
    nodes, box = _bar(cmap=Colormap.named('coolwarm').fitted([-.4, .6]))
    strip = Box.of([n.bounds for n in nodes if isinstance(n, Path)])
    assert box.max_x > strip.max_x, 'the labels are to the right of a vertical strip and count'
    # No contains loop: the box is the union of those bounds, so it contains every one by construction.


def test_placement_puts_the_bar_clear_of_the_content_on_the_right():
    style = DepictStyle()
    nodes, box = _bar(side='right')
    content = Box(0., 0., 6., 4.)
    placed = place_colorbar(nodes, box, content, style, 'right')
    assert placed.bounds.min_x >= content.max_x


def test_placement_puts_the_bar_below_the_content_on_the_bottom():
    style = DepictStyle()
    nodes, box = _bar(side='bottom')
    content = Box(0., 0., 6., 4.)
    placed = place_colorbar(nodes, box, content, style, 'bottom')
    assert placed.bounds.max_y <= content.min_y


def test_auto_puts_a_tall_molecule_beside_and_a_wide_one_below():
    from chython.depict.overlay import AtomField

    style = DepictStyle()
    overlays = (AtomField({1: .1, 2: -.1}),)
    assert legend_side(style, overlays, Box(0., 0., 2., 9.)) == 'right'
    assert legend_side(style, overlays, Box(0., 0., 9., 2.)) == 'bottom'


def test_auto_shows_nothing_when_no_overlay_carries_a_scale():
    from chython.depict.overlay import Highlight

    style = DepictStyle()
    assert legend_side(style, (Highlight(atoms=[1]),), Box(0., 0., 2., 9.)) is None
    assert legend_side(style, (), Box(0., 0., 2., 9.)) is None


def test_none_suppresses_the_bar_even_when_an_overlay_has_a_scale():
    from chython.depict.overlay import AtomField

    style = DepictStyle().tuned(**{'page.legend': 'none'})
    assert legend_side(style, (AtomField({1: .1}),), Box(0., 0., 2., 9.)) is None


def test_an_explicit_side_overrides_the_shape_heuristic():
    from chython.depict.overlay import AtomField

    style = DepictStyle().tuned(**{'page.legend': 'bottom'})
    assert legend_side(style, (AtomField({1: .1}),), Box(0., 0., 2., 9.)) == 'bottom'


def test_two_overlays_with_different_scales_are_refused_rather_than_drawn_wrong():
    """one bar cannot label two domains, and drawing the first silently is the bug"""
    from chython.depict.overlay import AtomField, BondScale

    style = DepictStyle()
    overlays = (AtomField({1: -.4, 2: .6}), BondScale({(1, 2): 900.}, encode='color'))
    with raises(ValueError, match='two different scales'):
        legend_side(style, overlays, Box(0., 0., 6., 4.))


def test_no_swatches_gives_no_bar_rather_than_an_empty_frame():
    nodes, box = _bar(swatches=[])
    assert nodes == []
    assert (box.min_x, box.min_y, box.max_x, box.max_y) == (0., 0., 0., 0.)
