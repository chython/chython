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
"""Shipped glyph metrics: the depictor measures before it places.

The table is DATA on purpose -- there is no measurement API a library can call in SVG, and three
backends measuring separately would disagree.  These tests pin two different things: that the shipped
tables are complete and sane (always run), and that they still match the AFMs they came from (skipped
when matplotlib is absent, the `needs_inchi` shape).
"""
from pytest import approx, mark, raises, importorskip
from chython.depict.metrics import (FAMILIES, METRICS_DIRECTORY, PDF_BASE_FONT, advance,
                                    font_metrics, text_box)
from chython.depict.scene import Text, TextRun


ELEMENT_LETTERS = set('ABCDEFGHIKLMNOPRSTUVWXYZabcdefghiklmnoprstuvwxyz')
LABEL_EXTRAS = set('0123456789+-.()[]*−')


@mark.parametrize('family', sorted(FAMILIES))
def test_every_character_a_label_can_contain_is_in_the_table(family):
    """no fallback, no guessed width: a glyph the depictor can emit must be measurable, a missing one
    surfacing as a mis-centred label in a picture"""
    table = font_metrics(family)
    missing = sorted((ELEMENT_LETTERS | LABEL_EXTRAS) - set(table))
    assert not missing, f'{family} cannot measure {missing}'


@mark.parametrize('family', sorted(FAMILIES))
def test_widths_are_positive_and_plausible(family):
    for char, glyph in font_metrics(family).items():
        # Upper bound is 2000, not 1000: some glyphs legitimately exceed 1 em (Helvetica @ is 1015).
        assert 0 < glyph.wx <= 2000, (char, glyph)
        assert glyph.llx <= glyph.urx and glyph.lly <= glyph.ury, (char, glyph)


def test_helvetica_capital_c_is_its_documented_width():
    """one hard number, so a regenerated table that silently shifted is caught"""
    glyph = font_metrics('helvetica')['C']
    assert (glyph.wx, glyph.llx, glyph.lly, glyph.urx, glyph.ury) == (722, 44, -19, 681, 737)


def test_the_two_families_differ():
    assert font_metrics('helvetica')['C'].wx != font_metrics('times')['C'].wx


def test_an_unknown_family_is_refused_by_name():
    with raises(ValueError, match='family'):
        font_metrics('comic')


def test_the_family_list_agrees_with_the_pdf_font_names():
    """a family the metrics know and PDF cannot name would be a runtime failure in the PDF backend"""
    assert set(PDF_BASE_FONT) == set(FAMILIES)


def test_advance_scales_with_size_and_sums_over_characters():
    table = font_metrics('helvetica')
    expected = (table['C'].wx + table['l'].wx) / 1000. * .4
    assert advance('Cl', 'helvetica', .4) == approx(expected)
    assert advance('Cl', 'helvetica', .8) == approx(2. * expected)


def test_advance_of_nothing_is_nothing():
    assert advance('', 'helvetica', .4) == 0.


def test_advance_refuses_a_character_it_cannot_measure():
    """silently substituting a width is how a label ends up off-centre with nothing to blame"""
    with raises(KeyError, match='中'):
        advance('中', 'helvetica', .4)


def test_a_start_anchored_box_begins_at_the_anchor():
    metrics = font_metrics('helvetica')
    box = text_box(Text([TextRun('CH', size=.4)], x=1., y=2., anchor='start'))
    assert box.min_x == approx(1. + metrics['C'].llx / 1000. * .4)
    # Pinned, not bounded: the pen advances C.wx before H, so H's ink ends at (C.wx + H.urx)*scale.  An
    # inequality also passes with the pen advance between glyphs missing.
    assert box.max_x == approx(1. + (metrics['C'].wx + metrics['H'].urx) / 1000. * .4)


def test_a_middle_anchored_box_straddles_the_anchor():
    box = text_box(Text([TextRun('CH', size=.4)], x=0., y=0., anchor='middle'))
    assert box.min_x < 0. < box.max_x
    assert abs(abs(box.min_x) - abs(box.max_x)) < .02, 'roughly symmetric about the anchor'


def test_an_end_anchored_box_ends_at_the_anchor():
    # H's right ink edge is at x - advance('H') + H.urx_scaled = x + (H.urx - H.wx) * scale.  The full
    # advance of 'CH' would wrongly subtract C.wx_scaled too: the pen before H is there, not at zero.
    box = text_box(Text([TextRun('CH', size=.4)], x=1., y=0., anchor='end'))
    assert box.max_x == approx(1. - (advance('H', 'helvetica', .4)
                                     - font_metrics('helvetica')['H'].urx / 1000. * .4))


def test_the_box_is_the_ink_and_not_the_em():
    """the label knock-out has to be the size of the INK, or a symbol sits in a hole too big for it"""
    box = text_box(Text([TextRun('c', size=1.)], x=0., y=0.))
    assert box.height < .8, 'a lowercase c has no ascender and no descender; its box must say so'


def test_a_subscript_run_lowers_the_box_and_extends_it():
    """CH3 -- the 3's dy is negative because the scene is y-up, and the box must follow it down"""
    plain = text_box(Text([TextRun('CH', size=.4)], x=0., y=0.))
    with_sub = text_box(Text([TextRun('CH', size=.4), TextRun('3', size=.28, dy=-.12)], x=0., y=0.))
    assert with_sub.max_x > plain.max_x, 'the subscript advances the label'
    assert with_sub.min_y < plain.min_y, 'and hangs below it'


def test_runs_advance_from_where_the_previous_run_ended():
    two_runs = text_box(Text([TextRun('C', size=.4), TextRun('H', size=.4)], x=0., y=0.))
    one_run = text_box(Text([TextRun('CH', size=.4)], x=0., y=0.))
    assert two_runs.max_x == approx(one_run.max_x), 'splitting a label into runs must not move it'


def test_a_runs_dx_shifts_only_that_run():
    shifted = text_box(Text([TextRun('C', size=.4), TextRun('H', size=.4, dx=.1)], x=0., y=0.))
    plain = text_box(Text([TextRun('C', size=.4), TextRun('H', size=.4)], x=0., y=0.))
    assert shifted.max_x == approx(plain.max_x + .1)


@mark.parametrize('family', sorted(FAMILIES))
def test_the_shipped_table_still_matches_the_afm_it_came_from(family, tmp_path):
    """the drift gate: regenerate into a temporary directory and diff, as the core's TSV tests do

    Skipped when matplotlib is absent -- it is the AFM SOURCE, a developer dependency, and the shipped
    tables are complete without it.  `scripts/` skips for a stronger reason: no wheel installs it, so a
    generator that cannot be reached is an absent gate rather than a failing one.
    """
    importorskip('matplotlib', reason='matplotlib supplies the source AFMs')
    generator = importorskip('scripts.gen_font_metrics', reason='`scripts/` is not installed with the '
                                                                'package; run from a source checkout')

    generator.compile_tables(tmp_path)
    shipped = METRICS_DIRECTORY.joinpath(f'{family}.tsv').read_text(encoding='utf-8')
    assert (tmp_path / f'{family}.tsv').read_text(encoding='utf-8') == shipped, \
        'regenerate with `python scripts/gen_font_metrics.py compile`'


@mark.parametrize('family', sorted(FAMILIES))
def test_the_adobe_paragraph_survives_the_round_trip(family):
    """Adobe's permission paragraph is a LICENCE CONDITION on its own wording ("this paragraph is not
    modified"), and the header wraps it, so unwrapping must give back every character.

    A reflow that dropped or joined a word breaches the licence while every other test stays green.
    Skipped like the drift gate above, the wording it compares against living in `scripts/`.
    """
    generator = importorskip('scripts.gen_font_metrics', reason='`scripts/` is not installed with the '
                             'package; run from a source checkout')
    ADOBE_COPYRIGHT, ADOBE_PARAGRAPH = generator.ADOBE_COPYRIGHT, generator.ADOBE_PARAGRAPH

    header = []
    for line in METRICS_DIRECTORY.joinpath(f'{family}.tsv').read_text(encoding='utf-8').splitlines():
        if not line.startswith('#'):
            break
        header.append(line[1:].strip())

    assert ADOBE_COPYRIGHT in header, 'the copyright notice must be retained verbatim'
    joined = ' '.join(x for x in header if x)
    assert ADOBE_PARAGRAPH in joined, 'the permission paragraph must round-trip character for character'
    assert 'MODIFICATION NOTICE' in joined, 'a derived table must prominently note that it is derived'
