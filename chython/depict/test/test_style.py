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
"""The style tree: every rendering constant, in one place, immutable.

Two properties a flat process-wide settings dict cannot give: a caller may hold two styles at once (a
paper's figure and its SI), and a typo'd key is refused rather than a silent no-op.  Both tested here.
"""
from dataclasses import FrozenInstanceError, replace
from pytest import approx, raises
from chython.depict.style import (PRESETS, AtomStyle, BondStyle, DepictStyle, HighlightStyle,
                                  LabelStyle, PageStyle, get_depict_style, set_depict_style)


def test_the_default_style_is_complete_and_frozen():
    style = DepictStyle()
    assert style.bond.width > 0.
    assert style.label.size > 0.
    assert style.page.margin >= 0.
    with raises(FrozenInstanceError):
        style.bond.width = .1


def test_nested_replace_is_the_plain_dataclass_road():
    """`tuned()` is sugar, not a new mechanism -- `replace` has to keep working or the tree is a DSL"""
    style = DepictStyle()
    other = replace(style, bond=replace(style.bond, width=.09))
    assert other.bond.width == approx(.09)
    assert style.bond.width != approx(.09), 'the original is untouched'


def test_tuned_takes_dotted_keys_and_returns_a_new_style():
    style = DepictStyle()
    other = style.tuned(**{'bond.width': .055, 'label.size': .45})
    assert (other.bond.width, other.label.size) == approx((.055, .45))
    assert other is not style
    assert style.bond.width != approx(.055)


def test_tuned_leaves_every_untouched_field_alone():
    style = DepictStyle()
    other = style.tuned(**{'bond.width': .055})
    assert other.bond.spacing == approx(style.bond.spacing)
    assert other.atom == style.atom
    assert other.page == style.page


def test_an_unknown_dotted_key_is_refused_and_the_message_names_the_field():
    """a misspelled key must be refused, not accepted as a silent no-op"""
    with raises(KeyError, match='bond.widht'):
        DepictStyle().tuned(**{'bond.widht': .055})
    with raises(KeyError, match='bonds'):
        DepictStyle().tuned(**{'bonds.width': .055})


def test_an_undotted_key_is_refused_with_a_message_that_says_how_to_spell_it():
    with raises(KeyError, match='dotted'):
        DepictStyle().tuned(width=.055)


def test_tuned_reaches_the_second_level_where_the_tree_has_one():
    style = DepictStyle().tuned(**{'field.contour.refine': 4})
    assert style.field.contour.refine == 4
    assert DepictStyle().field.contour.refine == 2, 'and the original is untouched'


def test_the_bond_defaults_that_carry_a_measured_number_are_pinned():
    """the four measured numbers behind the reference look, asserted because every drawing test TUNES
    them first and so nothing else reads the defaults at all

    TWO INSETS, deliberately: `aromatic_inset` was chosen against the circle, and the dashed inner line
    sits closer to its bond.  `test_bonds.py` tunes each of the two alone.
    """
    bond = BondStyle()
    assert bond.aromatic_inset == approx(.22), "the CIRCLE's inset from the ring bonds"
    assert bond.aromatic_dash_inset == approx(.14), "V2's aromatic_space: the dashed arc sits closer"
    assert bond.aromatic_dashes == approx((.15, .05))
    assert bond.dative_dashes == approx((.2, .1)), "V2's dative pattern"
    assert not hasattr(bond, 'dashes'), \
        'the field is named for its notation: a bare `dashes` beside `aromatic_dashes` says nothing'
    assert not hasattr(bond, 'dative_head'), \
        'order 8 draws headless: the direction is not a fact a container carries, so the knob read nothing'


def test_the_default_aromatic_notation_is_the_dashed_inner_ring():
    """one field default, and the presets inherit it because each builds a fresh `BondStyle`

    The process default is asserted too: `_DEFAULT_STYLE` is built from `DepictStyle()` at import, so a
    preset-only assertion would pass on a `mol.depict()` that still drew alternating lines.
    """
    assert BondStyle().aromatic == 'dashed-inner'
    assert DepictStyle().bond.aromatic == 'dashed-inner'
    assert get_depict_style().bond.aromatic == 'dashed-inner'
    for name in PRESETS:
        assert DepictStyle.preset(name).bond.aromatic == 'dashed-inner', name


def test_the_other_two_aromatic_notations_are_still_reachable():
    for mode in ('kekule', 'circle'):
        assert DepictStyle().tuned(**{'bond.aromatic': mode}).bond.aromatic == mode
    with raises(ValueError, match='dashed-inner'):
        BondStyle(aromatic='dashed')


def test_the_default_style_writes_the_map_numbers_a_structure_carries():
    assert AtomStyle().map_numbers is True
    assert DepictStyle().atom.map_numbers is True
    assert get_depict_style().atom.map_numbers is True
    for name in PRESETS:
        assert DepictStyle.preset(name).atom.map_numbers is True, name


def test_enhanced_stereo_groups_are_drawn_by_default():
    """withholding `&1` draws a single enantiomer where the file said the centre is racemic"""
    assert AtomStyle().stereo_groups is True
    assert DepictStyle().atom.stereo_groups is True
    for name in PRESETS:
        assert DepictStyle.preset(name).atom.stereo_groups is True, name
    assert DepictStyle().tuned(**{'atom.stereo_groups': False}).atom.stereo_groups is False


def test_the_two_annotation_rows_are_pinned():
    """the numbers that keep a descriptor off a map number, in fractions of `label.size`"""
    label = LabelStyle()
    assert label.annotation_rise == approx(.40), 'the stereo row, above the baseline'
    assert label.annotation_drop == approx(.40), 'the map row, below it'
    with raises(ValueError, match='annotation rise'):
        LabelStyle(annotation_rise=-.1)
    with raises(ValueError, match='annotation drop'):
        LabelStyle(annotation_drop=-.1)


def test_a_bad_value_is_refused_at_construction_and_not_at_draw_time():
    with raises(ValueError, match='width'):
        BondStyle(width=-.05)
    with raises(ValueError, match='aromatic inset'):
        BondStyle(aromatic_inset=-1.)
    with raises(ValueError, match='aromatic dash inset'):
        BondStyle(aromatic_dash_inset=-1.)
    with raises(ValueError, match='family'):
        LabelStyle(family='comic')
    with raises(ValueError, match='colour'):
        AtomStyle(carbon_colour='cornflowerblue')
    with raises(ValueError, match='margin'):
        PageStyle(margin=-1.)


def test_colours_are_normalized_on_the_way_in():
    assert AtomStyle(carbon_colour='#ABC').carbon_colour == '#aabbcc'
    assert AtomStyle(carbon_colour='black').carbon_colour == '#000000'


def test_every_preset_is_a_complete_style():
    for name in PRESETS:
        style = DepictStyle.preset(name)
        assert isinstance(style, DepictStyle)
        assert style.bond.width > 0. and style.label.size > 0.


def test_the_acs_preset_differs_from_the_default_and_says_how():
    acs = DepictStyle.preset('acs')
    assert acs != DepictStyle()
    assert acs.label.family == 'helvetica'
    assert acs.page.width_mm == approx(83.), 'ACS single-column is 83 mm'


def test_a_preset_is_tunable_which_is_the_documented_idiom():
    style = DepictStyle.preset('acs').tuned(**{'bond.width': .055, 'label.size': .45})
    assert (style.bond.width, style.label.size) == approx((.055, .45))
    assert style.page.width_mm == approx(83.), 'and the rest of the preset survives'


def test_an_unknown_preset_is_refused_and_the_message_lists_the_real_ones():
    with raises(ValueError, match='nature'):
        DepictStyle.preset('nature')


def test_the_process_default_round_trips():
    original = get_depict_style()
    try:
        set_depict_style(DepictStyle.preset('acs'))
        assert get_depict_style() == DepictStyle.preset('acs')
    finally:
        set_depict_style(original)
    assert get_depict_style() == original


def test_setting_a_non_style_as_the_default_is_refused():
    with raises(TypeError, match='DepictStyle'):
        set_depict_style({'bond': {'width': .05}})


def test_a_style_is_hashable_so_a_scene_can_be_cached_against_it():
    assert hash(DepictStyle()) == hash(DepictStyle())
    assert hash(DepictStyle().tuned(**{'bond.width': .09})) != hash(DepictStyle())


def test_two_styles_coexist():
    """a paper's figure and its SI, in one process"""
    thin = DepictStyle().tuned(**{'bond.width': .04})
    thick = DepictStyle().tuned(**{'bond.width': .09})
    assert (thin.bond.width, thick.bond.width) == approx((.04, .09))


def test_a_bad_intermediate_segment_is_refused_and_the_message_names_the_dotted_path():
    """three-segment key, wrong middle: `field.contoru.refine` names `field.contoru` in the error"""
    with raises(KeyError, match='field.contoru'):
        DepictStyle().tuned(**{'field.contoru.refine': 4})


def test_a_negative_highlight_outline_width_is_refused():
    with raises(ValueError, match='outline width'):
        HighlightStyle(outline_width=-.01)
