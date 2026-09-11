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
"""The SVG backend: one y-flip, one scale, and no <mask> anywhere.

Three properties pinned here: the y-flip is PER COORDINATE, not `scale(S, -S)`, which would mirror the
glyphs and need a counter-transform on every label; `<mask>` never appears, bond trimming being analytic;
and the output is DETERMINISTIC, so two renders of one molecule are one string.
"""
from gzip import decompress
from re import findall, search
from pytest import approx, raises
from chython.depict.metrics import text_box
from chython.depict.render.svg import format_number, to_svg, to_svgz
from chython.depict.scene import BLACK, Box, Group, Path, Scene, Text, TextRun, circle, close, curve, line, move
from chython.depict.style import DepictStyle, PageStyle


def _one_line(y0=0., y1=1.):
    return Scene([Path([[move(0., y0), line(1., y1)]], stroke=BLACK, width=.1)])


def test_the_document_opens_and_closes_with_the_process_default_style():
    """exercises get_depict_style() -- at least one test must cover the default path"""
    svg = to_svg(_one_line())
    assert svg.startswith('<svg')
    assert svg.rstrip().endswith('</svg>')
    assert 'xmlns="http://www.w3.org/2000/svg"' in svg


def test_no_mask_and_no_uuid():
    svg = to_svg(_one_line(), style=DepictStyle())
    assert 'mask' not in svg
    assert 'clip-path' not in svg, 'nothing in this scene is clipped, so nothing should be emitted'
    assert not findall(r'[0-9a-f]{8}-[0-9a-f]{4}-', svg), 'no generated id: the output must be stable'


def test_the_output_is_byte_identical_across_renders():
    assert to_svg(_one_line(), style=DepictStyle()) == to_svg(_one_line(), style=DepictStyle())


def test_y_is_negated_per_coordinate_and_there_is_no_negative_scale():
    """a point at y=+1 in molecule space must come out ABOVE one at y=0, without mirroring anything"""
    svg = to_svg(_one_line(0., 1.), style=DepictStyle())
    d = search(r'\sd="([^"]+)"', svg).group(1)
    assert d.startswith('M0 0'), d
    assert '-1' in d, 'the endpoint at y=+1 must be emitted as y=-1'
    assert 'scale(' not in svg, 'a negative scale would mirror the glyphs'
    assert 'transform' not in svg


def test_the_viewbox_is_the_scene_bounds_plus_the_margin():
    style = DepictStyle().tuned(**{'page.margin': .5})
    svg = to_svg(Scene([Path([[move(0., 0.), line(2., 1.)]], fill=BLACK)]), style=style)
    numbers = [float(v) for v in search(r'viewBox="([^"]+)"', svg).group(1).split()]
    assert numbers == approx([-.5, -1.5, 3., 2.]), 'min_x, -max_y, width, height'


def test_the_physical_size_comes_from_the_page_style():
    scene = Scene([Path([[move(0., 0.), line(2., 0.)]], fill=BLACK)])
    svg = to_svg(scene, style=DepictStyle().tuned(**{'page.scale_mm': 10., 'page.margin': 0.}))
    assert search(r'width="([^"]+)"', svg).group(1) == '20mm'

    wide = DepictStyle(page=PageStyle(width_mm=83., scale_mm=None, margin=0.))
    svg = to_svg(scene, style=wide)
    assert search(r'width="([^"]+)"', svg).group(1) == '83mm'
    # `-max_y` is -0. here and `format_number` normalizes it, so `"0 -0 2 0"` cannot be emitted
    assert 'viewBox="0 0 2 0"' in svg


def test_a_stated_scene_bounds_wins_over_the_content():
    svg = to_svg(Scene([Path([[move(0., 0.), line(1., 0.)]], fill=BLACK)], bounds=Box(-2., -2., 4., 4.)),
                 style=DepictStyle().tuned(**{'page.margin': 0.}))
    assert 'viewBox="-2 -4 6 6"' in svg


def test_a_stated_scene_bounds_is_not_inflated_by_margin():
    """stating bounds fixes the frame; a margin around it would make it something other than what was stated"""
    svg = to_svg(Scene([Path([[move(0., 0.), line(1., 0.)]], fill=BLACK)], bounds=Box(-2., -2., 4., 4.)),
                 style=DepictStyle().tuned(**{'page.margin': .5}))
    assert 'viewBox="-2 -4 6 6"' in svg, 'stated bounds survive a non-zero margin unchanged'


def test_a_stroked_path_carries_its_paint_attributes():
    svg = to_svg(Scene([Path([[move(0., 0.), line(1., 0.)]], stroke='#123456', width=.05,
                             cap='round', join='round')]), style=DepictStyle())
    assert 'stroke="#123456"' in svg
    assert 'stroke-width="0.05"' in svg
    assert 'stroke-linecap="round"' in svg
    assert 'stroke-linejoin="round"' in svg
    assert 'fill="none"' in svg, 'an unfilled path must say so or SVG fills it black'


def test_default_cap_and_join_are_not_emitted():
    """SVG's defaults are butt and miter; writing them again is bytes in every figure for nothing"""
    svg = to_svg(_one_line(), style=DepictStyle())
    assert 'stroke-linecap' not in svg
    assert 'stroke-linejoin' not in svg


def test_a_miter_limit_is_emitted_only_where_a_miter_join_is_in_force():
    """`stroke-miterlimit` has no meaning for a round or bevel join, so it is written only with miter

    A long thin wedge needs it: at a sharp angle the default limit of 4 turns the spike into a bevel.
    """
    corner = [[move(0., 0.), line(1., 0.), line(1., 1.)]]
    mitred = Path(corner, stroke=BLACK, width=.05, miter_limit=8.)
    assert 'stroke-miterlimit="8"' in to_svg(Scene([mitred]), style=DepictStyle())

    rounded = Path(corner, stroke=BLACK, width=.05, join='round', miter_limit=8.)
    assert 'stroke-miterlimit' not in to_svg(Scene([rounded]), style=DepictStyle())


def test_a_runs_dx_is_emitted_and_is_not_negated():
    """only y is flipped: a dx is an advance along x, which the scene and SVG already agree about"""
    svg = to_svg(Scene([Text([TextRun('C', size=.4), TextRun('H', size=.4, dx=.1)], x=0., y=0.)]),
                 style=DepictStyle())
    assert 'dx="0.1"' in svg
    assert 'dx="-0.1"' not in svg


def test_dashes_become_a_dasharray():
    svg = to_svg(Scene([Path([[move(0., 0.), line(1., 0.)]], stroke=BLACK, width=.02,
                             dashes=(.09, .07))]), style=DepictStyle())
    assert 'stroke-dasharray="0.09 0.07"' in svg


def test_a_curve_is_emitted_as_one_C_command():
    svg = to_svg(Scene([Path([[move(0., 0.), curve(0., 1., 1., 1., 1., 0.)]], fill=BLACK)]),
                 style=DepictStyle())
    d = search(r'\sd="([^"]+)"', svg).group(1)
    assert d == 'M0 0C0 -1 1 -1 1 0', d


def test_several_subpaths_are_one_path_element():
    scene = Scene([Path([[move(0., 0.), line(1., 0.)], [move(0., .2), line(1., .2)]],
                        stroke=BLACK, width=.04)])
    svg = to_svg(scene, style=DepictStyle())
    assert svg.count('<path') == 1, 'a double bond is ONE path: one element, one join, one paint'
    assert svg.count('M') == 2


def test_a_closed_subpath_ends_with_Z():
    svg = to_svg(Scene([Path([[move(0., 0.), line(1., 0.), line(1., 1.), close()]], fill=BLACK)]),
                 style=DepictStyle())
    assert search(r'\sd="([^"]+)"', svg).group(1).endswith('Z')


def test_even_odd_is_emitted_when_asked_and_not_otherwise():
    hole = Path([circle(0., 0., 1.), circle(0., 0., .5)], fill=BLACK, even_odd=True)
    assert 'fill-rule="evenodd"' in to_svg(Scene([hole]), style=DepictStyle())
    assert 'fill-rule' not in to_svg(Scene([Path([circle(0., 0., 1.)], fill=BLACK)]), style=DepictStyle())


def test_a_label_is_one_text_element_with_a_tspan_per_run():
    scene = Scene([Text([TextRun('CH', size=.4), TextRun('3', size=.28, dy=-.12)],
                        x=1., y=2., anchor='middle', fill=BLACK)])
    svg = to_svg(scene, style=DepictStyle())
    assert svg.count('<text') == 1
    assert svg.count('<tspan') == 2
    assert 'text-anchor="middle"' in svg
    assert 'y="-2"' in svg, 'the baseline y is negated like every other coordinate'


def test_a_subscript_run_carries_a_positive_dy_in_svg_space():
    """the run's dy is NEGATIVE in the y-up scene and must come out POSITIVE in SVG"""
    svg = to_svg(Scene([Text([TextRun('C', size=.4), TextRun('3', size=.28, dy=-.12)], x=0., y=0.)]),
                 style=DepictStyle())
    assert 'dy="0.12"' in svg


def test_two_shifted_runs_emit_the_difference_because_svg_dy_is_cumulative():
    """`TextRun.dy` is absolute from the label's baseline; `<tspan dy>` is relative to the pen

    TWO shifted runs distinguish the readings -- a label with one is identical under both.  `NH2+` is the
    shape an atom label is built from, and emitting each run's own `dy` puts its charge 0.112 below where
    `text_box` measured it, so the knock-out rectangle covers white space and the glyph sits outside it.
    """
    label = Text([TextRun('N', size=.4), TextRun('H', size=.4),
                  TextRun('2', size=.28, dy=-.112), TextRun('+', size=.28, dy=.16)], x=0., y=0.)
    svg = to_svg(Scene([label]), style=DepictStyle())

    shifts = [float(v) for v in findall(r'dy="([^"]+)"', svg)]
    assert shifts == approx([.112, -.272]), 'the deltas between consecutive runs, not the shifts'

    # summed the way a viewer sums them, they put the pen at the superscript's own absolute dy...
    assert -sum(shifts) == approx(.16)
    # ...which is the height `text_box` measured that run at, so mask and glyph agree
    assert text_box(label).max_y == approx(text_box(Text([TextRun('+', size=.28, dy=.16)])).max_y)


def test_the_font_family_names_a_family_svg_can_resolve():
    """PS_NAME carries `Times-Roman` which is a PostScript font name; CSS matches family names, so
    SVG_FAMILY carries the real names that viewers resolve"""
    from xml.etree.ElementTree import fromstring
    svg = to_svg(Scene([Text([TextRun('C', size=.4, family='times')], x=0., y=0.)]),
                 style=DepictStyle(), standalone=True)
    root = fromstring(svg)
    ns = 'http://www.w3.org/2000/svg'
    tspan = root.find(f'.//{{{ns}}}tspan')
    assert tspan is not None
    assert tspan.get('font-family') == '"Times New Roman",Times,serif'


def test_bold_and_italic_runs_say_so():
    svg = to_svg(Scene([Text([TextRun('R', size=.4, style='italic', weight='bold')], x=0., y=0.)]),
                 style=DepictStyle())
    assert 'font-style="italic"' in svg
    assert 'font-weight="bold"' in svg


def test_text_is_xml_escaped():
    svg = to_svg(Scene([Text([TextRun('<&>', size=.4)], x=0., y=0.)]), style=DepictStyle())
    assert '&lt;&amp;&gt;' in svg
    assert '<&>' not in svg


def test_a_group_with_opacity_becomes_a_g_with_group_opacity():
    scene = Scene([Group([Path([[move(0., 0.), line(1., 0.)]], fill=BLACK)], opacity=.35)])
    svg = to_svg(scene, style=DepictStyle())
    assert 'opacity="0.35"' in svg
    assert 'fill-opacity' not in svg, 'GROUP opacity, so overlapping children composite once'


def test_a_group_without_opacity_or_clip_emits_no_wrapper():
    """structure the output does not need is bytes in every figure and a node in every DOM"""
    svg = to_svg(Scene([Group([Path([[move(0., 0.), line(1., 0.)]], fill=BLACK)])]), style=DepictStyle())
    assert '<g' not in svg
    assert '<path' in svg


def test_a_clip_becomes_a_clippath_with_a_deterministic_id():
    scene = Scene([Group([Path([[move(0., 0.), line(1., 0.)]], fill=BLACK)],
                         clip=Path([circle(0., 0., 1.)], fill=BLACK))])
    svg = to_svg(scene, style=DepictStyle())
    ids = findall(r'id="([^"]+)"', svg)
    assert ids == ['clip0'], 'sequential, so two renders of one scene are the same string'
    assert 'clip-path="url(#clip0)"' in svg


def test_a_clip_with_even_odd_emits_clip_rule():
    """a donut clip -- two subpaths with even_odd -- clips as a hole; without clip-rule it is solid"""
    clip = Path([circle(0., 0., 2.), circle(0., 0., .5)], fill=BLACK, even_odd=True)
    scene = Scene([Group([Path([[move(0., 0.), line(1., 0.)]], fill=BLACK)], clip=clip)])
    svg = to_svg(scene, style=DepictStyle())
    assert 'clip-rule="evenodd"' in svg


def test_a_clip_without_even_odd_does_not_emit_clip_rule():
    """SVG's default clip-rule is nonzero; writing it is bytes in every figure for nothing"""
    clip = Path([circle(0., 0., 1.)], fill=BLACK)
    scene = Scene([Group([Path([[move(0., 0.), line(1., 0.)]], fill=BLACK)], clip=clip)])
    svg = to_svg(scene, style=DepictStyle())
    assert 'clip-rule' not in svg


def test_numbers_are_short_and_have_no_exponent():
    assert format_number(0.) == '0'
    assert format_number(1.) == '1'
    assert format_number(-0.) == '0', 'negative zero would make two renders differ for nothing'
    assert format_number(.10000000001) == '0.1'
    assert format_number(1e-9) == '0'
    assert format_number(1234.5) == '1234.5'
    assert 'e' not in format_number(1e-7)


def test_format_number_four_decimal_places_is_the_quantum():
    """two later backends inherit this function; pin the rounding quantum so it cannot drift"""
    assert format_number(1 / 3) == '0.3333'
    assert format_number(2.00001) == '2', 'rounds to integer when the fractional part vanishes at 4 dp'


def test_format_number_rejects_non_finite():
    """inf and nan in path data are not valid SVG; an early ValueError beats a viewer that chokes"""
    with raises(ValueError):
        format_number(float('inf'))
    with raises(ValueError):
        format_number(float('nan'))


def test_an_unknown_primitive_raises():
    """_render's else branch.  A stated bounds keeps scene.frame() off the unknown child's .bounds."""
    with raises(TypeError):
        to_svg(Scene([object()], bounds=Box(0., 0., 1., 1.)))  # type: ignore[list-item]


def test_standalone_adds_the_xml_declaration_and_the_default_does_not():
    assert to_svg(_one_line(), style=DepictStyle()).startswith('<svg')
    assert to_svg(_one_line(), style=DepictStyle(), standalone=True).startswith('<?xml')


def test_a_background_is_a_rect_and_transparent_is_nothing():
    style_with_bg = DepictStyle().tuned(**{'page.background': '#ffffff'})
    assert '<rect' in to_svg(_one_line(), style=style_with_bg)
    assert '<rect' not in to_svg(_one_line(), style=DepictStyle())


def test_svgz_is_gzip_of_the_standalone_document():
    style = DepictStyle()
    data = to_svgz(_one_line(), style=style)
    assert data[:2] == b'\x1f\x8b'
    assert decompress(data).decode() == to_svg(_one_line(), style=style, standalone=True)


def test_svgz_is_reproducible():
    """gzip writes an mtime by default, so two renders of one figure would differ byte for byte"""
    style = DepictStyle()
    assert to_svgz(_one_line(), style=style) == to_svgz(_one_line(), style=style)


def test_an_empty_scene_renders_an_empty_document_framed_on_the_origin():
    """with no children `Scene.bounds` collapses to the ORIGIN, not to the inverted empty box

    The alternative is an `inf` in the viewBox, which `format_number` refuses.  The frame is the margin
    around a point -- .35 each way, so .7 by .7 -- and not a zero-extent viewBox.
    """
    svg = to_svg(Scene([]), style=DepictStyle())
    assert svg.startswith('<svg') and svg.endswith('</svg>')
    assert 'viewBox="-0.35 -0.35 0.7 0.7"' in svg, 'the default .35 margin around the origin'
    assert '<path' not in svg and '<text' not in svg, 'nothing to draw, so nothing is drawn'


def test_a_zero_width_scene_with_a_stated_page_width_does_not_divide_by_zero():
    """`width_mm / width` has no answer for a scene with no width, so the scale falls back to 0

    The document then states 0mm and a zero-extent viewBox, whose rendering the SVG spec disables -- not
    a crash, and not a minimum extent nothing could justify.
    """
    style = DepictStyle(page=PageStyle(width_mm=83., scale_mm=None, margin=0.))
    svg = to_svg(Scene([Path([[move(1., 0.), line(1., 2.)]], fill=BLACK)]), style=style)
    assert 'width="0mm"' in svg and 'height="0mm"' in svg
    assert 'viewBox="1 -2 0 2"' in svg, 'the geometry is still there; only the physical scale collapsed'


def test_the_document_parses_as_xml():
    """the cheap total check: every attribute quoted, every element closed, every entity escaped"""
    from xml.etree.ElementTree import fromstring

    scene = Scene([Group([Path([circle(0., 0., 1.)], fill=BLACK, even_odd=True),
                          Text([TextRun('<C>', size=.4), TextRun('3', size=.3, dy=-.1)],
                               x=0., y=0., anchor='middle')], opacity=.5,
                         clip=Path([circle(0., 0., 2.)], fill=BLACK))])
    fromstring(to_svg(scene, standalone=True))
