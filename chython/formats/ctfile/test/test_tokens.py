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
"""The V3000 physical-line layer, tested below the grammar.

Join first, tokenize second: a quoted V3000 value may be split across a continuation, so a tokenizer
run on physical lines sees half a string.  Every declared count and width in CTfile is a hint, so a
disagreement is reported and the values kept, and no byte is discarded.
"""

from .._tokens import V30_PREFIX, emit_v30, join_continuations, parse_list, quote_value, tokenize


def _phys(*bodies):
    """Physical lines with the prefix attached, the way they appear in a file."""
    return [V30_PREFIX + b for b in bodies]


# joining, first

def test_quoted_string_split_across_continuation():
    """``FIELDDISP`` is a fixed-layout quoted string that exceeds 80 columns, so real files break it
    mid-string; joining before tokenizing makes it one value."""
    joined = join_continuations(_phys('FIELDDISP="    0.0000    0.00-', '00    DR    ALL"'))
    assert joined == ['FIELDDISP="    0.0000    0.0000    DR    ALL"']
    tokens = tokenize(joined[0])
    assert tokens == ['FIELDDISP="    0.0000    0.0000    DR    ALL"'], tokens


def test_the_join_is_byte_exact_and_does_not_strip_the_continued_part():
    """The spaces on either side of a break are inside the quotes and are part of the value, so
    ``rstrip()``ing a part corrupts a layout string."""
    assert join_continuations(_phys('A="x   -', '   y"')) == ['A="x      y"']


def test_a_line_without_the_prefix_keeps_its_first_seven_characters():
    """Slicing at column 7 unconditionally eats seven bytes of a prefixless line, so it is kept."""
    log = []
    assert join_continuations(['BEGIN ATOM'], log) == ['BEGIN ATOM']
    assert any('without M  V30 prefix' in x for x in log), log


def test_the_prefix_is_accepted_one_space_short():
    """``M  V30`` with no trailing space is what some writers emit for an empty continuation, and
    treating it as prefixless would log a false alarm on every such line."""
    assert join_continuations(['M  V30 BEGIN ATOM', 'M  V30']) == ['BEGIN ATOM', '']


def test_line_endings_are_removed_but_nothing_else_is():
    assert join_continuations([V30_PREFIX + 'COUNTS 1 0 0 0 0\r\n']) == ['COUNTS 1 0 0 0 0']


def test_a_dangling_continuation_yields_what_there_is():
    """A trailing ``-`` with no following line is out of grammar; dropping the logical line would
    lose a whole block."""
    log = []
    assert join_continuations(_phys('ATOMS=(2 1 -'), log) == ['ATOMS=(2 1 ']
    assert any('dangling' in x for x in log), log


def test_a_hyphen_that_is_not_at_the_end_is_not_a_continuation():
    """Negative coordinates end in digits, but a value may legitimately contain a hyphen."""
    assert join_continuations(_phys('1 C -1.5 -2.5 0 0')) == ['1 C -1.5 -2.5 0 0']


# tokenizing, second

def test_positional_tokens_and_key_value_pairs_come_back_as_written():
    assert tokenize('1 C 0 0 0 0 CHG=1 MASS=15') == \
        ['1', 'C', '0', '0', '0', '0', 'CHG=1', 'MASS=15']


def test_a_quoted_value_keeps_its_spaces_and_its_quotes():
    """The quotes stay in the token: stripping them makes ``"1"`` and ``1`` indistinguishable, and
    for ``SEQID`` those are not the same thing.  The caller that knows the keyword's type unquotes."""
    assert tokenize('FIELDNAME="Molecular Weight"') == ['FIELDNAME="Molecular Weight"']


def test_a_doubled_quote_inside_a_quoted_value_does_not_end_it():
    assert tokenize('F="say ""hi"" now" G=2') == ['F="say ""hi"" now"', 'G=2']


def test_a_parenthesised_list_is_one_token():
    assert tokenize('ATOMS=(3 1 2 5) FIELDNAME=X') == ['ATOMS=(3 1 2 5)', 'FIELDNAME=X']


def test_a_closing_paren_inside_quotes_does_not_close_the_list():
    """A ``FIELDDATA`` list carrying a name with a parenthesis in it -- ``(R)-`` -- is the real case;
    counting parentheses without knowing about quotes truncates it."""
    assert tokenize('X=(1 "a) b") Y=2') == ['X=(1 "a) b")', 'Y=2']


def test_an_unquoted_value_with_spaces_ends_where_the_next_keyword_begins():
    """The spec forbids this and writers produce it anyway.  Splitting it into several positional
    tokens shifts every positional value after it, so an atom line's element reads from the wrong
    column."""
    log = []
    tokens = tokenize('FIELDDATA=Molecular Weight: 375,40 FIELDNAME=X', log)
    assert tokens == ['FIELDDATA=Molecular Weight: 375,40', 'FIELDNAME=X'], tokens
    assert any('unquoted value with spaces' in x for x in log), log


def test_an_unquoted_value_with_spaces_and_no_following_keyword_runs_to_end_of_line():
    log = []
    assert tokenize('FIELDDATA=Molecular Weight: 375,40', log) == \
        ['FIELDDATA=Molecular Weight: 375,40']
    assert any('unquoted value with spaces' in x for x in log), log


def test_a_bare_value_followed_by_another_keyword_is_not_reported_as_spaceful():
    """`LABEL=Boc CSTATE=(...)` is conforming: the value ends at the space, the look-ahead stops at
    that same space, and nothing is absorbed.  The report is for a value that grew past its first
    bare run, so it names `FIELDDATA` and not every non-final keyword in the file."""
    log = []
    assert tokenize('ATOMS=(2 1 2) LABEL=Boc CSTATE=(4 3 0 0 1) ESTATE=E', log) == \
        ['ATOMS=(2 1 2)', 'LABEL=Boc', 'CSTATE=(4 3 0 0 1)', 'ESTATE=E']
    assert not log, log


def test_a_spaceful_value_taken_to_end_of_line_does_not_absorb_trailing_layout_spaces():
    """Inside quotes trailing spaces are data; an unquoted value cannot express significant edge
    whitespace, so the end-of-line branch stops at the last non-space."""
    assert tokenize('FIELDDATA=a b   ') == ['FIELDDATA=a b']


def test_only_the_keyword_alphabet_can_end_a_spaceful_value():
    """The look-ahead accepts only the CTfile keyword alphabet -- upper case, digits, underscore,
    then ``=`` -- so a colon, a lower-case word or an ``=``-free capital cannot cut a value short."""
    assert tokenize('FIELDDATA=a: b Weight ALL x FIELDNAME=X') == \
        ['FIELDDATA=a: b Weight ALL x', 'FIELDNAME=X']


def test_a_lower_case_key_is_still_a_key_when_it_is_written_as_one():
    """The ``KEY=`` scan does not filter by case -- only the look-ahead bounding a spaceful value
    does.  A writer emitting ``fieldname=`` is nonconforming but unambiguous."""
    assert tokenize('fieldname=X') == ['fieldname=X']


def test_an_unterminated_quote_is_reported_and_the_rest_of_the_line_is_the_value():
    log = []
    assert tokenize('F="never closed', log) == ['F="never closed']
    assert any('unterminated quoted' in x for x in log), log


def test_an_unterminated_list_is_reported_and_the_rest_of_the_line_is_the_value():
    log = []
    assert tokenize('ATOMS=(3 1 2', log) == ['ATOMS=(3 1 2']
    assert any('unterminated parenthesised' in x for x in log), log


def test_leading_and_repeated_spaces_produce_no_empty_tokens():
    assert tokenize('   1    C   0') == ['1', 'C', '0']


def test_an_empty_line_tokenizes_to_nothing():
    assert tokenize('') == []


def test_a_bare_quoted_positional_token_is_handled():
    """``SDT``-style positional strings exist, and a positional token is scanned by a different
    branch than a ``KEY=`` value."""
    assert tokenize('1 "a b" 2') == ['1', '"a b"', '2']


# list values

def test_a_list_is_parsed_from_either_the_raw_token_or_the_parenthesised_part():
    assert parse_list('ATOMS=(3 1 2 5)') == ['1', '2', '5']
    assert parse_list('(3 1 2 5)') == ['1', '2', '5']


def test_the_declared_list_count_is_checked_and_not_trusted():
    """The count is a writer's claim: truncating to it drops S-group atoms, padding invents refs."""
    log = []
    assert parse_list('(9 1 2)', log) == ['1', '2']
    assert any('count 9 disagrees with 2' in x for x in log), log


def test_a_list_without_a_leading_count_yields_all_of_its_items():
    log = []
    assert parse_list('(a b c)', log) == ['a', 'b', 'c']
    assert any('without a leading count' in x for x in log), log


def test_an_empty_list_is_empty_and_not_an_error():
    assert parse_list('ATOMS=()') == []


# writing back

def test_a_value_is_quoted_only_when_the_grammar_requires_it():
    assert quote_value('X') == 'X'
    assert quote_value('a b') == '"a b"'
    assert quote_value('') == '""', 'an empty value has no bare spelling at all'
    assert quote_value('a"b') == '"a""b"'


def test_a_value_ending_in_a_hyphen_is_quoted_because_the_hyphen_is_the_continuation_marker():
    """`LABEL=NH3+Cl-` bare puts a data hyphen last on its physical line, and the reader joins the
    next line onto the value: the S-group block's own `END SGROUP` becomes part of the label.
    Asserted through the joiner, which is the reader that has to get the value back.
    """
    assert quote_value('NH3+Cl-') == '"NH3+Cl-"'
    assert join_continuations(emit_v30('LABEL=' + quote_value('NH3+Cl-'))
                              + [V30_PREFIX + 'END SGROUP']) == ['LABEL="NH3+Cl-"', 'END SGROUP']


def test_a_well_formed_list_is_emitted_bare_even_though_it_holds_spaces():
    """Quoting turns a list into a string, and a reader looking for ``CSTATE=(4 ...)`` needs a list."""
    assert quote_value('(4 1 2 3)') == '(4 1 2 3)'


def test_a_short_line_is_emitted_as_one_physical_line():
    assert emit_v30('COUNTS 1 0 0 0 0') == [V30_PREFIX + 'COUNTS 1 0 0 0 0']


def test_every_emitted_physical_line_fits_the_eighty_column_limit():
    lines = emit_v30('FIELDDATA=' + 'x' * 300)
    assert all(len(x) <= 80 for x in lines), [len(x) for x in lines]
    assert all(x.startswith(V30_PREFIX) for x in lines)


def test_wrapping_round_trips_byte_exactly_through_the_joiner():
    """The wrap position is a free choice: the reader reassembles it identically.  Asserted over runs
    of significant spaces, where an off-by-one in the break would show."""
    for content in ('A=1 B=2',
                    'FIELDDISP="' + ' ' * 60 + 'DR    ALL  0       0"',
                    'FIELDDATA=' + 'y' * 400,
                    'X=' + 'a b ' * 40):
        assert join_continuations(emit_v30(content)) == [content], content


def test_the_writer_does_not_break_inside_a_quoted_string():
    """Strict out, permissive in: the reader tolerates a mid-string break, the writer never emits one."""
    content = 'FIELDNAME=A FIELDDISP="' + 'q' * 90 + '"'
    lines = emit_v30(content)
    # the break must land before the opening quote, so the first line carries no quote at all
    assert lines[0].count('"') == 0, lines


def test_a_single_unwrappable_token_still_produces_valid_lines():
    """A 4000-character ``FIELDDATA`` is one token with no break point, so the break falls where it
    must -- refusing would refuse a real data field."""
    lines = emit_v30('FIELDDATA=' + 'z' * 4000)
    assert all(len(x) <= 80 for x in lines)
    assert len(lines) > 50
    assert join_continuations(lines) == ['FIELDDATA=' + 'z' * 4000]


def test_the_two_operations_compose_in_the_documented_order():
    """End to end: emit, join, tokenize returns the tokens that went in.  Emit, tokenize, join does not."""
    tokens_in = ['1', 'DAT', '0', 'ATOMS=(1 4)', 'FIELDNAME="Molecular Weight"',
                 'FIELDDISP="    0.0000    0.0000    DR    ALL  0       0"',
                 'FIELDDATA=375,40']
    lines = emit_v30(' '.join(tokens_in))
    assert len(lines) > 1, 'the fixture must actually wrap or it tests nothing'
    assert tokenize(join_continuations(lines)[0]) == tokens_in
