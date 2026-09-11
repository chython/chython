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
"""The five overlay kinds, and the layer each belongs in.

An overlay returns (under, over) and never sorts: fields and halos go UNDER the structure, value labels
and bond-width scaling are part of it or go OVER.  Overlapping highlights are ONE group with group
opacity, never a boolean union -- two 50%-opaque discs drawn independently make a 75%-opaque lens where
they meet, reading as a third highlighted region.  Filled contour bands nest by the same mechanism.
"""
from math import hypot

from pytest import approx, raises
from chython import smiles
from chython.depict.colormap import Colormap
from chython.depict.field import convex_hull, in_polygon
from chython.depict.overlay import (AtomField, AtomHalo, BondScale, Highlight, ValueLabels,
                                    bands_of, render_overlays, scale_of)
from chython.depict.label import labels
from chython.depict.scene import Group, Text
from chython.depict.style import DepictStyle


def _phenol():
    """a public compound with a heteroatom, a ring and a labelled position"""
    mol = smiles('Oc1ccccc1')
    mol.clean2d()
    style = DepictStyle()
    return mol, mol.coordinates(), labels(mol, mol.coordinates(), style), style


def _charges(mol):
    """synthetic per-atom scalars, shaped like a charge map: negative on O, alternating on the ring"""
    return {a.n: (-.45 if a.atomic_symbol == 'O' else .05 * (-1) ** a.n)
            for a in mol.atoms()}


def test_a_highlight_is_a_disc_under_the_structure():
    mol, plane, boxes, style = _phenol()
    under, over = Highlight(atoms=[1, 2]).render(mol, plane, boxes, style)
    assert over == []
    assert len(under) == 1 and isinstance(under[0], Group), 'one group, so opacity composites once'
    assert len(under[0].children) == 2


def test_a_highlight_group_carries_the_opacity_and_the_children_do_not():
    """the whole point: two overlapping discs must not darken where they meet"""
    mol, plane, boxes, style = _phenol()
    under, _ = Highlight(atoms=[1, 2], color='#ffd54f').render(mol, plane, boxes, style)
    # 0.45 is the style.highlight.opacity default.  Path has no fill_opacity field, so per-child opacity
    # is structurally impossible and group opacity is the only compositor available.
    assert under[0].opacity == approx(0.45)


def test_highlighting_a_bond_covers_the_span_between_its_atoms():
    mol, plane, boxes, style = _phenol()
    under, _ = Highlight(bonds=[(2, 3)]).render(mol, plane, boxes, style)
    box = under[0].children[0].bounds
    for n in (2, 3):
        # BOTH axes: on a vertical or horizontal bond the two atoms share one coordinate, so an x-only
        # test is satisfied by a disc drawn on one endpoint alone
        assert box.min_x - 1e-9 <= plane[n][0] <= box.max_x + 1e-9
        assert box.min_y - 1e-9 <= plane[n][1] <= box.max_y + 1e-9


def test_highlighting_an_atom_the_molecule_does_not_have_is_refused():
    mol, plane, boxes, style = _phenol()
    with raises(ValueError, match='999'):
        Highlight(atoms=[999]).render(mol, plane, boxes, style)


def test_highlighting_a_pair_that_is_not_bonded_is_refused():
    mol, plane, boxes, style = _phenol()
    with raises(ValueError, match='not bonded'):
        Highlight(bonds=[(1, 4)]).render(mol, plane, boxes, style)


def test_a_halo_colours_each_named_atom_by_its_value():
    mol, plane, boxes, style = _phenol()
    values = _charges(mol)
    oxygen = [a.n for a in mol.atoms() if a.atomic_symbol == 'O'][0]
    under, _ = AtomHalo(values).render(mol, plane, boxes, style)
    discs = {id(child): child for child in under[0].children}
    assert len(discs) == len(values)
    fills = {child.fill for child in under[0].children}
    assert len(fills) > 1, 'different values, different colours'


def test_a_halo_appears_only_on_the_atoms_named():
    mol, plane, boxes, style = _phenol()
    oxygen = [a.n for a in mol.atoms() if a.atomic_symbol == 'O'][0]
    under, _ = AtomHalo({oxygen: -.45}).render(mol, plane, boxes, style)
    assert len(under[0].children) == 1


def test_a_halo_on_an_empty_mapping_draws_nothing_and_does_not_divide_by_zero():
    mol, plane, boxes, style = _phenol()
    assert AtomHalo({}).render(mol, plane, boxes, style) == ([], [])


def test_a_field_produces_nested_filled_bands_outermost_first():
    """painter's order IS the nesting: each band is painted over the ones outside it, so the visible
    colour anywhere is the innermost band covering that point -- and no union is computed"""
    mol, plane, boxes, style = _phenol()
    under, over = AtomField(_charges(mol)).render(mol, plane, boxes, style)
    assert over == []
    # ONE LEVEL IS TWO SIBLING PATHS by default: the fill, then its boundary stroke.  Order the FILLS
    # only -- over all children the areas come in equal pairs and a non-strict sort passes by accident.
    fills = [p for p in under[0].children if p.fill is not None]
    assert len(fills) >= 3
    areas = [(p.bounds.max_x - p.bounds.min_x) * (p.bounds.max_y - p.bounds.min_y) for p in fills]
    assert areas == sorted(areas, reverse=True), 'outermost band painted first'
    assert all(a > b for a, b in zip(areas, areas[1:])), 'and strictly, so equal pairs cannot pass it'


def test_a_field_band_is_filled_with_no_stroke_and_is_made_of_cubics():
    """the fill Path carries no stroke of its own -- the boundary is its own sibling Path, so a printer
    that drops strokes still gets the bands and one that drops fills still gets the contours"""
    mol, plane, boxes, style = _phenol()
    under, _ = AtomField(_charges(mol)).render(mol, plane, boxes, style)
    band = under[0].children[0]
    assert band.fill is not None and band.stroke is None
    assert any(seg[0] == 'C' for sub in band.subpaths for seg in sub)
    outline = under[0].children[1]
    assert outline.stroke is not None and outline.fill is None, 'the boundary is the next sibling'


def test_an_unfilled_field_is_stroked_isolines_instead():
    mol, plane, boxes, style = _phenol()
    under, _ = AtomField(_charges(mol), fill=False).render(mol, plane, boxes, style)
    line = under[0].children[0]
    assert line.stroke is not None and line.fill is None


def test_a_field_with_explicit_levels_draws_exactly_those():
    """the three levels are MEASURED to lie inside this field and to trace one isoline each -- see the
    note below on why the obvious [-.2, 0., .2] does not"""
    mol, plane, boxes, style = _phenol()
    under, _ = AtomField(_charges(mol), levels=[-.35, -.25, -.15], fill=False).render(mol, plane, boxes,
                                                                                      style)
    assert len(under[0].children) == 3


def test_a_level_the_field_never_reaches_draws_nothing_and_is_not_an_error():
    """the field is a weighted SUM, so its range is NOT the values' range -- a caller asking for a level
    above the maximum is making an ordinary mistake, not an invalid request"""
    mol, plane, boxes, style = _phenol()
    under, _ = AtomField(_charges(mol), levels=[-.25, 5.], fill=False).render(mol, plane, boxes, style)
    assert len(under[0].children) == 1


def _end_points(path):
    """The on-curve points of a Path: a move and a line state theirs, a cubic's is its last pair."""
    out = []
    for sub in path.subpaths:
        for seg in sub:
            if seg[0] in ('M', 'L'):
                out.append((seg[1], seg[2]))
            elif seg[0] == 'C':
                out.append((seg[5], seg[6]))
    return out


def test_a_field_is_unclipped_by_default_because_its_bands_lie_outside_its_own_atoms():
    """A BAND IS A RING AROUND THE ATOMS, so a hull drawn through them cuts it.

    Measured on this phenol charge field: 173 of the 432 band vertices fall outside the padded convex hull
    of the valued atoms, across 5 of the 9 bands, sitting 0.570--1.216 from the nearest valued atom against
    a pad of 0.55.  The field limits itself at `contour.cutoff` instead, a boundary from the chemistry.
    """
    mol, plane, boxes, style = _phenol()
    overlay = AtomField(_charges(mol))
    assert overlay.clip is None, 'the DEFAULT is what this test is about'
    under, _ = overlay.render(mol, plane, boxes, style)
    assert under[0].clip is None, 'and no clip path reaches the group'

    hull = convex_hull([plane[a] for a in overlay.values], style.field.contour.pad)
    fills = [p for p in under[0].children if p.fill is not None]
    outside = [(p, v) for p in fills for v in _end_points(p) if not in_polygon(v[0], v[1], hull)]
    assert outside, 'no vertex outside the hull would mean the hull had nothing to cut'
    assert len({id(p) for p, _ in outside}) > 1, 'and it is not one stray outermost ring'
    excess = max(min(hypot(v[0] - plane[a][0], v[1] - plane[a][1]) for a in overlay.values)
                 for _, v in outside)
    assert excess > 2. * style.field.contour.pad, 'nor a rounding-width excursion past the pad'


def test_asking_for_a_hull_clip_still_cuts_the_drawing():
    """'hull' and 'box' stay available: the clip is a renderer clip PATH, so the geometry is untouched and
    the cut shows only in the bounds, the Group intersecting its children's box with the clip's."""
    mol, plane, boxes, style = _phenol()
    # clip=None spelled out: this test is about 'hull' and must not read the default's value.
    free = AtomField(_charges(mol), clip=None).render(mol, plane, boxes, style)[0][0]
    clipped = AtomField(_charges(mol), clip='hull').render(mol, plane, boxes, style)[0][0]

    assert clipped.clip is not None, 'a clip was asked for'
    assert [_end_points(a) for a in clipped.children] == [_end_points(a) for a in free.children], \
        'the same geometry: a clip path cuts at render time and never edits a contour'
    for edge in ('min_x', 'min_y'):
        assert getattr(clipped.bounds, edge) > getattr(free.bounds, edge), f'{edge} was not cut'
    assert clipped.bounds.max_x < free.bounds.max_x, 'max_x was not cut'


def test_the_zero_level_of_a_mixed_sign_field_is_stroked_open_and_never_filled():
    """Level 0 of a mixed-sign field is a nodal line, open at both ends and enclosing no region, so it is
    stroked even in the default `fill=True` mode.  Measured: the two chains trace 90 and 104 vertices and
    draw 60 and 68, so trimming the asymptotic tail does not turn either into a band."""
    mol, plane, boxes, style = _phenol()
    log = []
    under, _ = AtomField(_charges(mol), levels=[0.]).render(mol, plane, boxes, style, log=log)
    assert under, 'a level whose only chains are open must still draw: the nodal line is chemistry'
    drawn = under[0].children
    assert drawn
    for path in drawn:
        assert path.fill is None and path.stroke is not None, 'stroked, not filled'
        for sub in path.subpaths:
            assert sub[-1][0] != 'Z', 'and left open: an open chain closed is an invented arc'

    assert log, 'a level that drew no band and dropped vertices must not do so silently'
    message = ' '.join(r.message for r in log)
    assert 'open contour' in message and 'asymptotic' in message
    # the traced chains reach 3.9 from the nearest valued atom, where the field is 4e-18 and the sign is
    # floating-point noise; the drawn ones stop at 2.41
    reach = max(min(hypot(x - plane[a][0], y - plane[a][1]) for a in _charges(mol))
                for path in drawn for x, y in _end_points(path))
    assert reach < 2.5, reach


def test_each_band_carries_the_interval_of_values_its_colour_covers():
    """A band's colour shows between its own level and the next band's level OUTWARD, so the level is an
    END of the interval and never its middle.  This phenol field runs -0.4388..+0.0340 inside a domain of
    ±0.45: the innermost band runs down to the domain edge, the outermost stops at -0.0132."""
    mol, plane, boxes, style = _phenol()
    overlay = AtomField(_charges(mol))
    swatches = sorted(bands_of(overlay, mol, plane, style))
    cmap = scale_of(overlay)
    assert len(swatches) == 9 and all(s.filled for s in swatches)

    for swatch in swatches:
        assert swatch.lo < swatch.hi
        assert swatch.level in (swatch.lo, swatch.hi), 'a level is an end of its interval, not its middle'
    for low, high in zip(swatches, swatches[1:]):
        assert low.hi == approx(high.lo), 'and the intervals meet, so no value falls between two bands'
    assert swatches[0].lo == approx(cmap.vmin)
    assert swatches[-1].hi < cmap.vmax, 'the part of the domain this field never reaches stays empty'


def test_a_level_that_traces_nothing_says_so_in_the_log_rather_than_vanishing():
    """`levels=[..., 5.]` is an ordinary caller mistake and not an error -- but the swatch count and the
    band count then disagree, and a reader has no way to tell which level went missing."""
    mol, plane, boxes, style = _phenol()
    log = []
    AtomField(_charges(mol), levels=[-.25, 5.], fill=False).render(mol, plane, boxes, style, log=log)
    assert log
    message = ' '.join(r.message for r in log)
    assert 'never reaches' in message and '+5' in message
    assert '-0.25:' in message, 'every level is accounted for, not only the missing one'


def test_a_plane_in_the_wrong_units_is_reported_rather_than_drawn_as_discs():
    """sigma is a multiple of the mean bond length and cutoff is absolute, so a plane in picometres puts
    the effective sigma outside the cutoff and every atom draws as a hard-edged disc.  One log line."""
    mol, plane, boxes, style = _phenol()
    picometres = {n: (x * 100., y * 100.) for n, (x, y) in plane.items()}
    log = []
    AtomField(_charges(mol)).render(mol, picometres, boxes, style, log=log)
    assert log, 'a field truncated inside its own sigma said nothing'
    message = ' '.join(r.message for r in log)
    assert 'sigma' in message.lower() and str(int(style.field.contour.cutoff)) in message

    quiet = []
    AtomField(_charges(mol)).render(mol, plane, boxes, style, log=quiet)
    assert not quiet, 'and an ordinary clean2d layout must not warn'


def test_bond_scale_widens_a_bond_and_leaves_the_others_alone():
    """this one does not draw: it returns per-bond width overrides that bonds.py consumes"""
    mol, plane, boxes, style = _phenol()
    widths = BondScale({(2, 3): 1.}).bond_widths(mol, style)
    assert widths[(2, 3)] > style.bond.width
    assert (3, 4) not in widths


def test_bond_scale_maps_the_data_range_onto_the_width_range():
    mol, plane, boxes, style = _phenol()
    scale = BondScale({(2, 3): 0., (3, 4): 1.}, width_range=(.02, .10))
    widths = scale.bond_widths(mol, style)
    assert widths[(2, 3)] == approx(.02)
    assert widths[(3, 4)] == approx(.10)


def test_bond_scale_keys_are_order_independent():
    mol, plane, boxes, style = _phenol()
    assert BondScale({(3, 2): 1.}).bond_widths(mol, style).keys() == {(2, 3)}


def test_bond_scale_can_also_colour():
    mol, plane, boxes, style = _phenol()
    scale = BondScale({(2, 3): 0., (3, 4): 1.}, encode='color', colormap=Colormap.named('coolwarm'))
    colours = scale.bond_colours(mol, style)
    assert colours[(2, 3)] != colours[(3, 4)]


def test_value_labels_go_over_the_structure_and_read_the_number():
    mol, plane, boxes, style = _phenol()
    oxygen = [a.n for a in mol.atoms() if a.atomic_symbol == 'O'][0]
    under, over = ValueLabels({oxygen: -.4512}).render(mol, plane, boxes, style)
    assert under == []
    text = over[0] if isinstance(over[0], Text) else over[0].children[0]
    assert ''.join(run.text for run in text.runs) == '−0.45', 'two places, typographic minus'


def test_the_caller_s_format_string_is_the_one_used():
    """`{:+.3f}` on .5 is '+0.500': the leading plus is KEPT and only the ASCII minus is swapped for the
    typographic one, a hyphen beside a numeral being the wrong glyph where '+' is not"""
    mol, plane, boxes, style = _phenol()
    over = ValueLabels({1: .5}, fmt='{:+.3f}').render(mol, plane, boxes, style)[1]
    text = over[0] if isinstance(over[0], Text) else over[0].children[0]
    assert ''.join(run.text for run in text.runs) == '+0.500'


def test_a_value_label_is_offset_clear_of_the_atom_label():
    mol, plane, boxes, style = _phenol()
    oxygen = [a.n for a in mol.atoms() if a.atomic_symbol == 'O'][0]
    over = ValueLabels({oxygen: -.45}).render(mol, plane, boxes, style)[1]
    text = over[0] if isinstance(over[0], Text) else over[0].children[0]
    assert not _overlaps(text.bounds, boxes[oxygen].box)


def test_a_bond_keyed_value_label_sits_at_the_bond_midpoint():
    mol, plane, boxes, style = _phenol()
    over = ValueLabels({(2, 3): 1.42}).render(mol, plane, boxes, style)[1]
    text = over[0] if isinstance(over[0], Text) else over[0].children[0]
    centre = ((text.bounds.min_x + text.bounds.max_x) / 2, (text.bounds.min_y + text.bounds.max_y) / 2)
    mid = ((plane[2][0] + plane[3][0]) / 2, (plane[2][1] + plane[3][1]) / 2)
    # "nearer the midpoint than either endpoint" is scale-free where a raw tolerance is not: the layout's
    # bond is .825, so `< .5` on x alone also passes for a label parked on an endpoint.
    assert hypot(centre[0] - mid[0], centre[1] - mid[1]) < min(
        hypot(centre[0] - plane[n][0], centre[1] - plane[n][1]) for n in (2, 3))


def test_render_overlays_concatenates_in_the_order_given():
    """two fields, and the caller's order decides which is on top -- the code does not guess"""
    mol, plane, boxes, style = _phenol()
    # told apart by child count -- a halo on two atoms, a highlight on one.  Identity cannot see order:
    # every call builds fresh Groups, so `reversed_under[0] is not under[0]` passes with order ignored.
    first = AtomHalo({1: 1., 2: -1.})
    second = Highlight(atoms=[3])
    under, over = render_overlays([first, second], mol, plane, boxes, style)
    assert [len(g.children) for g in under] == [2, 1]
    reversed_under, _ = render_overlays([second, first], mol, plane, boxes, style)
    assert [len(g.children) for g in reversed_under] == [1, 2]


def test_no_overlay_draws_anything_for_an_empty_list():
    mol, plane, boxes, style = _phenol()
    assert render_overlays([], mol, plane, boxes, style) == ([], [])


def test_highlight_label_goes_to_the_over_list_not_under():
    """a label rendered under bonds and atom symbols is occluded at 83 mm; it must be in over"""
    mol, plane, boxes, style = _phenol()
    under, over = Highlight(atoms=[1], label='A').render(mol, plane, boxes, style)
    assert any(isinstance(n, Text) for n in over), 'label must be in the over list'
    assert not any(isinstance(n, Text) for n in under), 'and not in under'


def test_scale_of_highlight_returns_none():
    assert scale_of(Highlight(atoms=[1])) is None


def test_scale_of_value_labels_returns_none():
    assert scale_of(ValueLabels({1: 0.5})) is None


def test_scale_of_bond_scale_width_only_returns_none():
    assert scale_of(BondScale({(2, 3): 0.5}, encode='width')) is None


def test_scale_of_bond_scale_color_returns_colormap():
    result = scale_of(BondScale({(2, 3): 0., (3, 4): 1.}, encode='color'))
    assert isinstance(result, Colormap)


def test_scale_of_atom_halo_returns_colormap():
    result = scale_of(AtomHalo({1: 0., 2: 1.}))
    assert isinstance(result, Colormap)


def test_scale_of_atom_field_returns_colormap():
    result = scale_of(AtomField({1: -0.4, 2: 0.4}))
    assert isinstance(result, Colormap)


def test_scale_of_and_render_fit_the_same_domain():
    """the bar is labelled by `scale_of` and the picture coloured by `render`, both through one `_fitted`
    so they cannot drift; the case that catches a drift is `domain=`"""
    values = {1: -.4, 2: .1}
    plain = scale_of(AtomField(values))
    assert (plain.vmin, plain.vmax) == (-.4, .4), 'a diverging map symmetrizes what it was fitted to'
    stated = scale_of(AtomField(values, domain=(-1., .5)))
    assert (stated.vmin, stated.vmax) == (-1., .5), 'and an explicit domain is used unchanged'
    # BondScale reaches the same helper through a different attribute path
    assert scale_of(BondScale(dict(zip([(2, 3), (3, 4)], values.values())), encode='color',
                              domain=(-1., .5))).vmin == -1.


# the style fields below were literals here while the field they duplicate went unread.  Each test says
# the same thing: TUNE THE FIELD, SEE THE PICTURE CHANGE.

def test_a_value_label_takes_its_format_from_the_style():
    mol, plane, boxes, style = _phenol()
    tuned = style.tuned(**{'field.value_format': '{:.4f}'})
    over = ValueLabels({1: .5}).render(mol, plane, boxes, tuned)[1]
    assert ''.join(r.text for r in over[0].runs) == '0.5000'
    # and the caller's own fmt still wins over the style's
    over = ValueLabels({1: .5}, fmt='{:+.1f}').render(mol, plane, boxes, tuned)[1]
    assert ''.join(r.text for r in over[0].runs) == '+0.5'


def test_a_value_label_takes_its_colour_and_size_from_the_style():
    mol, plane, boxes, style = _phenol()
    tuned = style.tuned(**{'field.value_colour': '#ff0000', 'field.value_scale': .9})
    text = ValueLabels({1: .5}).render(mol, plane, boxes, tuned)[1][0]
    assert text.fill == '#ff0000'
    assert text.runs[0].size == approx(tuned.label.size * .9)
    plain = ValueLabels({1: .5}).render(mol, plane, boxes, style)[1][0]
    assert plain.runs[0].size == approx(style.label.size * style.field.value_scale)


def test_bond_widths_span_the_style_s_declared_width_range():
    """`bond_min_width` / `bond_max_width` are the declaration; `style.bond.width * .5` and `* 2.5` would
    be a second spelling of one default"""
    mol, plane, boxes, style = _phenol()
    tuned = style.tuned(**{'field.bond_min_width': .05, 'field.bond_max_width': .5})
    widths = BondScale({(2, 3): 0., (3, 4): 1.}).bond_widths(mol, tuned)
    assert widths[(2, 3)] == approx(.05)
    assert widths[(3, 4)] == approx(.5)
    plain = BondScale({(2, 3): 0., (3, 4): 1.}).bond_widths(mol, style)
    assert plain[(2, 3)] == approx(style.field.bond_min_width)
    assert plain[(3, 4)] == approx(style.field.bond_max_width)


def test_a_field_group_takes_its_opacity_from_the_style():
    mol, plane, boxes, style = _phenol()
    tuned = style.tuned(**{'field.contour.fill_opacity': .2})
    under, _ = AtomField(_charges(mol)).render(mol, plane, boxes, tuned)
    assert under[0].opacity == approx(.2)
    # the overlay's own value still overrides the style's
    under, _ = AtomField(_charges(mol), opacity=.9).render(mol, plane, boxes, tuned)
    assert under[0].opacity == approx(.9)


def _overlaps(a, b):
    return a.min_x < b.max_x and b.min_x < a.max_x and a.min_y < b.max_y and b.min_y < a.max_y
