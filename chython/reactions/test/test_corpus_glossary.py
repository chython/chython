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
"""`docs/glossary.rst` is the two corpora, and this is the gate that keeps it from becoming a copy.

    functional.tsv + protective.tsv  ->  the tables in docs/glossary.rst
    `-- the authority                    `-- generated; the drift test below is what lets a reader
                                             trust a name they read there

The same shape `test_element_tables.py::test_the_compiled_tables_are_the_two_tsvs` has, for the same
reason: a page transcribed by hand states last month's corpus and says nothing when a row is added.
"""
from pytest import skip

from .gen_corpus_glossary import PAGE, compile_glossary, rst_text
from .._tables import functional_rules, protective_rules


def _page():
    if not PAGE.is_file():
        skip('no docs/ beside this package -- an installed copy, not a checkout')
    # `encoding='utf-8'`: the page holds an em dash and a prime, so the locale's codec decides whether
    # the generated text is found in it -- on cp1252 both decode to something else and the drift test
    # fails with the page unchanged.
    return PAGE.read_text(encoding='utf-8')


def test_the_page_is_the_two_tables():
    """Regenerating changes nothing, or the page has drifted from the corpora."""
    assert compile_glossary() in _page(), \
        'docs/glossary.rst has drifted from the corpora; run gen_corpus_glossary.py'


def test_every_row_of_both_corpora_is_listed():
    """A glossary that omits a row is worse than none: a reader concludes the group is absent."""
    text = _page()
    for name in functional_rules():
        assert '``%s``' % name in text, name
    for name in protective_rules():
        assert '``%s``' % name in text, name


def test_every_row_is_listed_by_id_as_well_as_by_name():
    """The id is what a consumer stores, so the page a reader looks a name up in resolves it to one."""
    text = _page()
    for rule in functional_rules().values():
        assert '``%s``' % rule.id in text, rule.id
    for rule in protective_rules().values():
        assert '``%s``' % rule.id in text, rule.id


def test_the_page_states_the_counts_it_lists():
    """Both counts are generated, so neither can be the number of rows there used to be."""
    text = _page()
    assert '%d functional groups' % len(functional_rules()) in text
    assert '%d protecting groups' % len(protective_rules()) in text


def test_a_markdown_code_span_in_a_description_becomes_an_rst_literal():
    """35 descriptions spell a name or a pattern in single backticks, which is a title reference in rst
    and a literal in Markdown.  Rewritten rather than escaped, so the column stays readable in the TSV."""
    assert rst_text('see `primary_amide`') == 'see ``primary_amide``'


def test_a_description_rst_would_read_as_markup_is_refused():
    """The generator refuses rather than emitting markup the TSV did not mean.

    A bare `*` or `_` outside a code span is emphasis and a reference in rst, and there is none in the
    corpora today -- so the moment one is written, the generator says so instead of rendering it.
    """
    from pytest import raises

    with raises(ValueError, match='reads as rst markup'):
        rst_text('a trailing reference_')
    with raises(ValueError, match='reads as rst markup'):
        rst_text('an *emphasis*')
