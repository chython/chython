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
"""RDfile framing and metadata.

``test/MR.rdf`` is the fixture that matters: four records, two molecule and two reaction, holding
three V2000 CTABs and three V3000 ones.  A reader that decides the version per file -- or that
handles only one record kind -- reads it wrong, which is why both are per-record here.
"""

from pytest import mark

from chython.formats.ctfile import split_rdf_records


# --- Splitter: record boundaries

def test_mr_rdf_has_four_records_of_both_kinds(root):
    with (root / 'test' / 'MR.rdf').open(encoding='utf8') as f:
        records = list(split_rdf_records(f))
    assert [tag for tag, _ in records] == ['$MFMT', '$MFMT', '$RFMT', '$RFMT']


def test_datum_value_with_trailing_plus_is_data_not_continuation():
    """A trailing '+' in a $DATUM value is data, not a continuation marker: the format has no such rule,
    CTfile p.46 using column position rather than a marker character.  Treating '+' as one makes the
    ``$MFMT`` after a value like '95%+' invisible and destroys the following record.
    """
    lines = ['$MFMT',
             'mol', '', '',
             '  0  0  0     0  0            999 V2000',
             'M  END',
             '$DTYPE yield',
             '$DATUM 95%+',
             '$MFMT',
             'mol2', '', '',
             '  0  0  0     0  0            999 V2000',
             'M  END']
    records = list(split_rdf_records(lines))
    assert len(records) == 2, [t for t, _ in records]
    # Also assert the value is intact: '95%+' is stored as-is, not stripped to '95%'
    from chython.formats.ctfile import parse_rdf_record
    first = parse_rdf_record(records[0][0], records[0][1], [])
    assert first.meta['yield'] == '95%+'


def test_rfmt_inside_a_datum_value_does_not_split_the_record():
    """A $RFMT line after an 80-char $DATUM line is value text, not a delimiter: the positional
    continuation rule (CTfile p.46) says a line of >= 80 characters continues onto the next one whatever
    that line contains.
    """
    long_datum = '$DATUM ' + 'x' * 73   # exactly 80 characters
    assert len(long_datum) == 80
    lines = ['$RDFILE 1',
             '$DATM    01/02/17 17:17',
             '$MFMT',
             'one', '', '',
             '  0  0  0     0  0            999 V2000',
             'M  END',
             '$DTYPE comment',
             long_datum,
             '$RFMT is not a record here',
             '$MFMT',
             'two', '', '',
             '  0  0  0     0  0            999 V2000',
             'M  END']
    records = list(split_rdf_records(lines))
    assert len(records) == 2, [t for t, _ in records]
    assert records[0][1][0] == 'one'
    assert records[1][1][0] == 'two'


def test_positional_continuation_suppression_is_logged():
    """When a record tag is suppressed inside a positional continuation, a log line is produced."""
    long_datum = '$DATUM ' + 'x' * 73
    assert len(long_datum) == 80
    log = []
    lines = ['$MFMT',
             'mol', '', '',
             '  0  0  0     0  0            999 V2000',
             'M  END',
             '$DTYPE k',
             long_datum,
             '$RFMT suppressed line',
             '$MFMT',
             'mol2', '', '',
             '  0  0  0     0  0            999 V2000',
             'M  END']
    records = list(split_rdf_records(lines, log))
    assert len(records) == 2
    assert any('suppressed' in e for e in log), log


def test_registry_reference_on_record_tag_is_logged():
    """A registry number on the ``$MFMT``/``$RFMT`` tag line is logged as unsupported, not stored."""
    log = []
    lines = ['$MFMT ext-reg-123',
             'one', '', '',
             '  0  0  0     0  0            999 V2000',
             'M  END']
    records = list(split_rdf_records(lines, log))
    assert len(records) == 1
    assert any('unsupported' in entry and 'registry' in entry for entry in log)


def test_line_before_first_record_is_logged():
    """A non-header line appearing before the first $MFMT/$RFMT tag is logged."""
    log = []
    lines = ['some garbage line',
             '$MFMT',
             'mol', '', '',
             '  0  0  0     0  0            999 V2000',
             'M  END']
    list(split_rdf_records(lines, log))
    assert len(log) == 1
    assert 'before the first record' in log[0]


# --- The region before the first record tag arms no positional continuation

def test_long_dtype_before_the_first_record_does_not_swallow_it():
    """A stray 80-column ``$DTYPE`` in the header must not consume the first record: no data field is open
    before any record tag, so nothing there can start a wrapped logical line, and treating it as one arms
    a continuation that eats the ``$MFMT`` behind it.
    """
    counts = '  0  0  0     0  0            999 V2000'
    stray = '$DTYPE junk'.ljust(80, 'x')
    assert len(stray.rstrip()) == 80   # long enough to arm, if anything here could arm
    lines = ['$RDFILE 1', stray,
             '$MFMT', 'mol', '', '', counts, 'M  END',
             '$MFMT', 'mol2', '', '', counts, 'M  END']
    records = list(split_rdf_records(lines, []))
    assert len(records) == 2, [t for t, _ in records]
    assert records[0][1][0] == 'mol'
    assert records[1][1][0] == 'mol2'


def test_stray_dtype_before_the_first_record_is_logged_exactly_once():
    """The stray header line is reported and the record behind it is not, so the log cannot be non-empty
    and false -- five ``before the first record`` lines, four about the first record's own tag and body,
    claim the file lacks a record it plainly has.
    """
    counts = '  0  0  0     0  0            999 V2000'
    stray = '$DTYPE junk'.ljust(80, 'x')
    log = []
    lines = ['$RDFILE 1', stray,
             '$MFMT', 'mol', '', '', counts, 'M  END']
    list(split_rdf_records(lines, log))
    assert len(log) == 1, log
    assert 'before the first record' in log[0]
    assert '$DTYPE' in log[0]
    assert not any('$MFMT' in entry or 'mol' in entry for entry in log), log


# --- One arming predicate, read by the splitter and by the field parser alike

def test_advance_reports_the_continuation_the_shared_predicate_exposes():
    """``_DataFieldState.advance`` must read ``continuation_open`` rather than recompute the conjunction:
    sharing the state is not sharing the decision, and a predicate written twice moves in one place only.
    """
    from chython.formats.ctfile._rdf import _DataFieldState

    state = _DataFieldState()
    for line in ['some ctab line', '$DTYPE k', '$DATUM ' + 'x' * 73, '$MFMT absorbed', 'tail',
                 '$MEREG 1', '$DATUM ' + 'y' * 73]:
        exposed = state.continuation_open
        was_continuation, _ = state.advance(line)
        assert was_continuation == exposed, line


# --- A registry cross-reference is the same construct wherever it appears

def test_registry_cross_reference_is_unsupported_in_both_positions():
    """``$MEREG`` in the structure body and in the data-field tail read alike: the ``unsupported: `` prefix
    means *the file is fine and we are the limitation*, which is true of a registry tag in either position.
    """
    from chython.formats.ctfile import parse_rdf_fields, parse_rdf_record

    registry = '$MEREG 42'
    body_log = []
    parse_rdf_record('$MFMT', [registry, 'one', '', '',
                               '  0  0  0     0  0            999 V2000', 'M  END'], body_log)
    tail_log = []
    parse_rdf_fields(['$DTYPE a', '$DATUM v', registry], tail_log)

    body = [e for e in body_log if 'MEREG' in e]
    tail = [e for e in tail_log if 'MEREG' in e]
    assert len(body) == 1, body_log
    assert len(tail) == 1, tail_log
    assert str(body[0]).startswith('unsupported: '), body[0]
    assert str(tail[0]).startswith('unsupported: '), tail[0]
    assert body[0] == tail[0], (body[0], tail[0])


# --- A $DATM timestamp: stored where the format puts it, reported where it does not

def test_the_header_timestamp_is_stored_and_the_tail_one_is_reported():
    """One position is valid and the other is not, so the two do not read alike.  ``$DATM`` on line 2 is
    the spec's file timestamp, kept verbatim in the splitter's ``header`` dict with nothing to log; in a
    ``$DTYPE``/``$DATUM`` tail it is reported *unprefixed*, since ``unsupported: `` would claim the file
    is fine and chython the limitation.
    """
    from chython.formats.ctfile import parse_rdf_fields

    stamp = '$DATM    01/02/17 17:17'
    header_log = []
    header = {}
    list(split_rdf_records(['$RDFILE 1', stamp,
                            '$MFMT', 'mol', '', '',
                            '  0  0  0     0  0            999 V2000', 'M  END'],
                           header_log, header=header))
    tail_log = []
    parse_rdf_fields(['$DTYPE a', '$DATUM v', stamp], tail_log)

    assert header == {'date': '01/02/17 17:17'}
    assert header_log == [], header_log

    tail = [e for e in tail_log if 'DATM' in e]
    assert len(tail) == 1, tail_log
    assert not str(tail[0]).startswith('unsupported: '), tail[0]
    assert 'unrecognised' not in tail[0], tail[0]


def test_a_second_header_timestamp_keeps_the_later_one_and_says_so():
    """One file, two ``$DATM`` lines: one of them is not the file's timestamp.

    Keeping the first and staying silent would leave the file looking well formed, so the later value
    wins -- the same last-one-wins rule a repeated data-field name gets -- and the displacement is a
    log line.
    """
    log = []
    header = {}
    list(split_rdf_records(['$RDFILE 1', '$DATM    01/02/17 17:17', '$DATM    02/03/18 18:18',
                            '$MFMT', 'mol', '', '',
                            '  0  0  0     0  0            999 V2000', 'M  END'], log, header=header))
    assert header == {'date': '02/03/18 18:18'}
    assert len(log) == 1, log


def test_a_well_formed_file_logs_nothing_at_all(root):
    """Every real RDfile has a ``$DATM`` on line 2, and it costs no log line: a predicate that is constant
    over a whole format carries no information about that format, and the construct has somewhere to go.
    """
    log = []
    header = {}
    with (root / 'test' / 'MR.rdf').open(encoding='utf8') as f:
        list(split_rdf_records(f, log, header=header))
    assert log == [], log
    assert header == {'date': '01/02/17 17:17'}


# --- Field parser: $DTYPE/$DATUM

@mark.parametrize('line', ['$DATUM MUTADAT', '$DATUMMUTADAT'])
def test_datum_value_starting_with_dollar_letters_is_not_eaten(line):
    """The value must be taken by prefix, not by ``lstrip('$DATUM')``, which strips a character *set*: the
    unspaced ``'$DATUMMUTADAT'`` is eaten whole by that, while ``'$DATUM MUTADAT'`` halts at the space and
    so discriminates nothing.  Both spellings are here for that contrast.
    """
    from chython.formats.ctfile import parse_rdf_fields

    fields = parse_rdf_fields(['$DTYPE registry', line], [])
    assert fields['registry'] == 'MUTADAT'


def test_rdf_fields_are_a_dict():
    """``{name: value}`` in file order, the shape ``meta`` holds."""
    from chython.formats.ctfile import parse_rdf_fields

    lines = ['$DTYPE TEMP', '$DATUM 100', '$DTYPE SOLVENT', '$DATUM water']
    assert parse_rdf_fields(lines) == {'TEMP': '100', 'SOLVENT': 'water'}


def test_a_repeated_dtype_merges_and_says_so():
    """One value per name is what a mapping holds, so the two value lines join and the collision is
    reported -- the same rule ``parse_data_fields`` applies to a repeated SDF field name."""
    from chython.formats.ctfile import parse_rdf_fields

    log = []
    assert parse_rdf_fields(['$DTYPE K', '$DATUM a', '$DTYPE K', '$DATUM b'], log) == {'K': 'a\nb'}
    assert any('appears twice' in x for x in log), log


def test_dtype_with_no_datum_yields_empty_value():
    """A ``$DTYPE`` followed immediately by another ``$DTYPE`` must not produce a field whose value
    is the string ``'None'``."""
    from chython.formats.ctfile import parse_rdf_fields

    fields = parse_rdf_fields(['$DTYPE a', '$DTYPE b', '$DATUM v'], [])
    assert fields['a'] == ''
    assert fields['b'] == 'v'


def test_datum_without_dtype_is_logged():
    """A ``$DATUM`` with no preceding ``$DTYPE`` produces exactly one log message."""
    from chython.formats.ctfile import parse_rdf_fields

    log = []
    parse_rdf_fields(['$DATUM orphan'], log)
    assert len(log) == 1
    assert '$DATUM' in log[0] and '$DTYPE' in log[0]


def test_second_datum_under_one_dtype_keeps_the_first_value_and_logs_the_other():
    """Two ``$DATUM`` lines under one ``$DTYPE``: the first value is kept, the second reported.  A field's
    value is the ``$DATUM`` that follows the ``$DTYPE`` naming it, so the second has no name of its own;
    silent replacement would leave the caller unable to learn the reader chose between two values.
    """
    from chython.formats.ctfile import parse_rdf_fields

    log = []
    fields = parse_rdf_fields(['$DTYPE a', '$DATUM v', '$DATUM w'], log)
    assert fields == {'a': 'v'}, fields
    assert len(log) == 1, log
    assert '$DATUM w' in log[0], log[0]
    assert not str(log[0]).startswith('unsupported: '), log[0]


def test_datum_permissive_continuation_appends_with_newline():
    """A non-``$`` line after a short ``$DATUM`` is a permissive continuation joined with ``\\n``."""
    from chython.formats.ctfile import parse_rdf_fields

    fields = parse_rdf_fields(['$DTYPE x',
                               '$DATUM first line',
                               'second line'], [])
    assert fields['x'] == 'first line\nsecond line'


def test_positional_continuation_value_is_complete():
    """An 80-char ``$DATUM`` line continues onto the next; the resulting value contains both parts
    concatenated with **no separator** -- which is what tells a positional continuation from the
    permissive one above, that joins with a newline."""
    from chython.formats.ctfile import parse_rdf_fields

    suffix = 'A' * 73
    datum = '$DATUM ' + suffix   # exactly 80 characters
    assert len(datum) == 80

    fields = parse_rdf_fields(['$DTYPE k', datum, 'tail'], [])
    assert fields['k'] == suffix + 'tail'


def test_positional_continuation_chain():
    """Three chained 80-char lines join into one value with no separator anywhere in it."""
    from chython.formats.ctfile import parse_rdf_fields

    suffix1 = 'A' * 73
    cont1 = 'B' * 80
    cont2 = 'C' * 5
    datum = '$DATUM ' + suffix1   # 80 chars
    assert len(datum) == 80
    assert len(cont1) == 80

    fields = parse_rdf_fields(['$DTYPE chain', datum, cont1, cont2], [])
    assert fields['chain'] == suffix1 + cont1 + cont2


def test_rfmt_inside_positional_continuation_is_in_value():
    """A ``$RFMT`` line that the splitter kept inside a positional continuation is concatenated into
    the field value, not silently dropped."""
    from chython.formats.ctfile import parse_rdf_fields

    suffix = 'x' * 73
    datum = '$DATUM ' + suffix   # 80 chars
    continuation_text = '$RFMT is not a delimiter here'
    fields = parse_rdf_fields(['$DTYPE comment', datum, continuation_text], [])
    assert fields['comment'] == suffix + continuation_text


# --- parse_rdf_record: container types and metadata

def test_mr_rdf_records_carry_the_right_container_and_meta(root):
    """The tag decides the container, and both kinds hold their ``$DTYPE`` pairs on ``meta``."""
    from chython.core import MoleculeContainer
    from chython.core.reaction import ReactionContainer
    from chython.formats.ctfile import parse_rdf_record

    with (root / 'test' / 'MR.rdf').open(encoding='utf8') as f:
        records = list(split_rdf_records(f))
    parsed = [parse_rdf_record(tag, lines, []) for tag, lines in records]

    assert [isinstance(r, ReactionContainer) for r in parsed] == [False, False, True, True]
    assert all(isinstance(r, MoleculeContainer) for r in parsed[:2])
    assert [r.meta['CdId'] for r in parsed] == ['MOL V2000', 'MOL V3000', 'RXN V2000', 'RXN V3000']


def test_the_ctab_versions_in_mr_rdf_are_genuinely_mixed(root):
    """Three V2000 CTABs and three V3000 -- the reason the sniff is per-CTAB and not per-file.

    The version is framing and no container holds it, so it comes out of the ``header=`` dict.
    """
    from chython.core.reaction import ReactionContainer
    from chython.formats.ctfile import V2000_STAMP, V3000_STAMP, parse_rdf_record

    with (root / 'test' / 'MR.rdf').open(encoding='utf8') as f:
        records = list(split_rdf_records(f))
    versions = []
    parsed = []
    for tag, lines in records:
        header = {}
        parsed.append(parse_rdf_record(tag, lines, [], header=header))
        versions.append(header['version'])
    assert versions == [V2000_STAMP, V3000_STAMP, V2000_STAMP, V3000_STAMP]
    assert sum(len(list(r.molecules())) for r in parsed
               if isinstance(r, ReactionContainer)) == 4


def test_reaction_record_header_fields(root):
    """The name line is ``reaction.title``; ``program`` and ``comment`` have no container home and
    come out of ``header=``."""
    from chython.core.reaction import ReactionContainer
    from chython.formats.ctfile import parse_rdf_record

    with (root / 'test' / 'MR.rdf').open(encoding='utf8') as f:
        records = list(split_rdf_records(f))
    # records[2] is the first $RFMT, whose $RXN name line is 'title3'
    header = {}
    reaction = parse_rdf_record(records[2][0], records[2][1], [], header=header)
    assert isinstance(reaction, ReactionContainer)
    assert reaction.title == 'title3'
    assert header['program'] == ''
    assert header['comment'] == ''


def test_mireg_line_in_record_body_is_logged_and_filtered():
    """A ``$MIREG`` line inside the record body is removed from the body AND logged as unsupported."""
    from chython.formats.ctfile import parse_rdf_record

    log = []
    lines = ['$MIREG 456',
             'one', '', '',
             '  0  0  0     0  0            999 V2000',
             'M  END']
    molecule = parse_rdf_record('$MFMT', lines, log)
    assert molecule.atom_count == 0
    assert any('unsupported' in e and 'registry' in e for e in log), log


# --- Column-80 padding does not arm positional continuation

def test_space_padded_datum_does_not_arm_continuation():
    """A $DATUM padded to column 80 with trailing spaces must not arm positional continuation.

    Content length is what the spec counts; vendor tools sometimes pad to column 80 with spaces.
    ``'$DATUM 95%'.ljust(80)`` is 80 characters but only 10 are content.
    """
    padded = '$DATUM 95%'.ljust(80)
    assert len(padded) == 80
    assert len(padded.rstrip()) == 10   # confirm it is padding, not content
    lines = ['$MFMT',
             'mol', '', '',
             '  0  0  0     0  0            999 V2000',
             'M  END',
             '$DTYPE yield',
             padded,
             '$MFMT',
             'mol2', '', '',
             '  0  0  0     0  0            999 V2000',
             'M  END']
    records = list(split_rdf_records(lines))
    assert len(records) == 2, records


# --- A wrapped $DTYPE name continues via positional continuation, not discards the $DATUM

def test_long_dtype_continues_name_not_discards_datum():
    """A $DTYPE line of >= 80 characters continues the field name onto the next line.  Without that, the
    $DATUM following is absorbed by a continuation branch that does nothing while ``current_lines`` is
    ``None``, and the field is emitted with an empty value and no log line.
    """
    from chython.formats.ctfile import parse_rdf_fields

    # $DTYPE padded to exactly 80 chars of content (no trailing spaces)
    long_dtype = ('$DTYPE name').ljust(80, 'x')   # 80 chars, no spaces
    assert len(long_dtype.rstrip()) == 80
    name_tail = '_suffix'   # continuation of the name
    fields = parse_rdf_fields([long_dtype, name_tail, '$DATUM value'], [])
    # The name is the concatenation of both lines, and the value is not empty: it was not lost.
    assert fields == {long_dtype[len('$DTYPE'):].strip() + name_tail: 'value'}


def test_long_dtype_absorbs_dollar_line_and_logs():
    """When the continuation of a long $DTYPE name starts with $, it is logged."""
    from chython.formats.ctfile import parse_rdf_fields

    long_dtype = ('$DTYPE name').ljust(80, 'x')
    assert len(long_dtype.rstrip()) == 80
    log = []
    fields = parse_rdf_fields([long_dtype, '$DATUM value'], log)
    # The $DATUM was absorbed into the name (positional continuation), logged, and emitted with
    # value == '' (no further $DATUM followed)
    assert len(log) == 1
    assert 'absorbed' in log[0]


# --- An absorbed keyword is reported whichever target swallowed it

def test_long_datum_absorbs_dollar_line_and_logs():
    """A keyword absorbed into a wrapped ``$DATUM`` value is reported, as one in a name is.  Composed with
    the rule that a second ``$DATUM`` cannot displace a stored value, a silent value branch is the
    difference between a record reporting two malformed lines and one reporting neither.
    """
    from chython.formats.ctfile import parse_rdf_fields

    long_datum = '$DATUM ' + 'v' * 73   # exactly 80 content characters
    assert len(long_datum.rstrip()) == 80
    log = []
    fields = parse_rdf_fields(['$DTYPE a', long_datum, '$DTYPE b'], log)
    assert len(fields) == 1, fields
    assert 'absorbed' in log[0], log
    assert '$DTYPE b' in log[0], log[0]


def test_the_absorbed_keyword_report_is_one_rule_for_both_targets():
    """The name branch and the value branch produce the same message, differing only in the target, so a
    change to the wording or to the "starts with a dollar" test moves both.
    """
    from chython.formats.ctfile import parse_rdf_fields

    absorbed = '$DTYPE b'
    name_log = []
    parse_rdf_fields([('$DTYPE a').ljust(80, 'x'), absorbed], name_log)
    value_log = []
    parse_rdf_fields(['$DTYPE a', '$DATUM ' + 'v' * 73, absorbed], value_log)

    assert len(name_log) == 1, name_log
    assert len(value_log) == 1, value_log
    assert str(name_log[0]).replace('$DTYPE name', '<target>') == \
           str(value_log[0]).replace('$DATUM value', '<target>'), (name_log[0], value_log[0])


def test_a_keyword_absorbed_into_a_value_leaves_no_silent_loss():
    """The composed case: a wrapped value swallows a ``$DTYPE`` and a second ``$DATUM`` follows.  Two
    fields' worth of content is malformed, so the reader keeps what it can -- ``a``'s wrapped value --
    and reports both the swallowed keyword and the displaced value rather than neither.
    """
    from chython.formats.ctfile import parse_rdf_fields

    log = []
    fields = parse_rdf_fields(['$DTYPE a', '$DATUM ' + 'v' * 73, '$DTYPE b', '$DATUM w'], log)
    assert fields == {'a': 'v' * 73 + '$DTYPE b'}, fields
    assert len(log) == 2, log
    assert 'absorbed' in log[0], log[0]
    assert '$DATUM w' in log[1], log[1]


# --- The 79/80 continuation boundary

def test_positional_continuation_boundary_79_is_not_armed():
    """A $DATUM line of exactly 79 content characters does not arm positional continuation."""
    datum_79 = ('$DATUM x').ljust(79, 'x')
    assert len(datum_79) == 79
    assert len(datum_79.rstrip()) == 79   # no trailing whitespace
    lines = ['$MFMT', 'mol', '', '',
             '  0  0  0     0  0            999 V2000', 'M  END',
             '$DTYPE k', datum_79,
             '$MFMT', 'mol2', '', '',
             '  0  0  0     0  0            999 V2000', 'M  END']
    records = list(split_rdf_records(lines))
    assert len(records) == 2


def test_positional_continuation_boundary_80_is_armed():
    """A $DATUM line of exactly 80 content characters DOES arm positional continuation."""
    datum_80 = ('$DATUM x').ljust(80, 'x')
    assert len(datum_80) == 80
    assert len(datum_80.rstrip()) == 80
    # The $MFMT that follows the 80-char line is a suppressed continuation, not a new record.
    lines = ['$MFMT', 'mol', '', '',
             '  0  0  0     0  0            999 V2000', 'M  END',
             '$DTYPE k', datum_80,
             '$MFMT suppressed',
             'extra', '', '',
             '  0  0  0     0  0            999 V2000', 'M  END']
    records = list(split_rdf_records(lines))
    assert len(records) == 1


# --- CRLF line endings

def test_crlf_terminated_lines_parse_correctly():
    """Files from Windows tooling have CRLF endings; no stray ``\\r`` should reach field values."""
    lines_crlf = ['$MFMT\r\n',
                  'mol\r\n',
                  '\r\n',
                  '\r\n',
                  '  0  0  0     0  0            999 V2000\r\n',
                  'M  END\r\n',
                  '$DTYPE name\r\n',
                  '$DATUM value\r\n',
                  '$MFMT\r\n',
                  'mol2\r\n',
                  '\r\n',
                  '\r\n',
                  '  0  0  0     0  0            999 V2000\r\n',
                  'M  END\r\n']
    records = list(split_rdf_records(lines_crlf))
    assert len(records) == 2
    # No stray \r in any record-body line
    for _, body in records:
        assert all('\r' not in ln for ln in body), body
    # No stray \r in the parsed field value
    from chython.formats.ctfile import parse_rdf_record
    first = parse_rdf_record(records[0][0], records[0][1], [])
    assert '\r' not in first.meta.get('name', ''), first.meta


# --- parse_rxn_record emits each version message exactly once

def test_parse_rxn_record_version_message_appears_exactly_once():
    """``parse_rxn_record`` sniffs the version into a throwaway log so ``parse_rxn``'s internal sniff
    is the one authoritative call.  A version-disagreement message must appear exactly once."""
    from chython.formats.ctfile import parse_rxn_record

    # $RXN without V3000 in the tag, but M  V30 lines in the body: sniff_rxn_version logs exactly
    # one message about this disagreement.  A double sniff would produce two.
    body = ['$RXN', 'title', '', '',
            'M  V30 COUNTS 0 0',
            'M  END']
    log = []
    parse_rxn_record(body, [], log)
    # The specific message sniff_rxn_version appends on this disagreement:
    version_msgs = [m for m in log if 'read as V3000' in m]
    assert len(version_msgs) == 1, version_msgs


# --- CTAB body lines do not arm positional continuation

def test_long_v30_line_in_ctab_body_does_not_suppress_next_record():
    """A V3000 body line of 82 characters must not arm positional continuation.

    ``M  V30`` lines can legitimately exceed column 80.  The continuation rule applies only to
    ``$DTYPE``/``$DATUM`` logical lines, not to CTAB body lines.  The long line is the immediate
    predecessor of ``$MFMT`` (no intervening line that would reset ``prev_long``).
    """
    long_v30 = 'M  V30 ' + 'x' * 75   # 7 + 75 = 82 content chars; immediately before $MFMT
    assert len(long_v30.rstrip()) == 82
    lines = ['$MFMT', 'mol', '', '',
             '  0  0  0     0  0            999 V2000',
             long_v30,     # last line of record 1's body, immediately followed by $MFMT
             '$MFMT',
             'mol2', '', '',
             '  0  0  0     0  0            999 V2000',
             'M  END']
    records = list(split_rdf_records(lines))
    assert len(records) == 2


# --- Unrecognised $-led lines end the data-field logical line without data loss

def test_unrecognised_dollar_line_is_logged_and_does_not_arm_continuation():
    """An 80-character ``$``-led line that is not ``$DTYPE``/``$DATUM`` must be logged and must not arm
    positional continuation -- otherwise it is dropped silently and the ``$DTYPE`` after it is absorbed
    into the previous value.  The 79-character variant witnesses the logging half on its own.
    """
    from chython.formats.ctfile import parse_rdf_fields

    long_dollar = '$100 for the reagent'.ljust(80, '.')
    assert len(long_dollar) == 80
    log = []
    fields = parse_rdf_fields(
        ['$DTYPE note', '$DATUM see below', long_dollar, '$DTYPE yield', '$DATUM 95'], log
    )
    assert fields == {'note': 'see below', 'yield': '95'}, fields
    assert len(log) == 1, log
    assert 'unrecognised' in log[0], log[0]


def test_short_unrecognised_dollar_line_is_logged_and_does_not_arm_continuation():
    """A 79-character ``$``-led non-keyword line is logged and the following ``$DTYPE`` is unharmed.  Below
    the arming threshold, so this proves "logged the drop" only; the 80-char test above covers the other
    half.
    """
    from chython.formats.ctfile import parse_rdf_fields

    short_dollar = '$100 for the reagent'.ljust(79, '.')
    assert len(short_dollar) == 79
    log = []
    fields = parse_rdf_fields(
        ['$DTYPE note', '$DATUM see below', short_dollar, '$DTYPE yield', '$DATUM 95'], log
    )
    assert fields == {'note': 'see below', 'yield': '95'}, fields
    assert len(log) == 1, log
    assert 'unrecognised' in log[0], log[0]


def test_unrecognised_dollar_line_inside_dtype_does_not_mangle_name():
    """A registry-style ``$``-led line inside a ``$DTYPE`` block does not absorb into the name: an 80-column
    ``$MEREG`` after ``$DTYPE a`` otherwise reaches the wrapped-name branch, giving ``'a$DTYPE b'`` and
    losing ``b``'s field entirely.
    """
    from chython.formats.ctfile import parse_rdf_fields

    mangler = '$MEREG 12'.ljust(80, 'x')
    assert len(mangler.rstrip()) == 80
    log = []
    fields = parse_rdf_fields(['$DTYPE a', mangler, '$DTYPE b', '$DATUM w'], log)
    # Field 'a' has no value (no $DATUM); field 'b' has value 'w'.  `'a$DTYPE b'` as a key is the
    # failure this is here to catch, so the names are asserted and not just the values.
    assert fields == {'a': '', 'b': 'w'}, fields


def test_splitter_and_parser_agree_on_unrecognised_dollar_line():
    """``split_rdf_records`` and ``parse_rdf_fields`` share the logical-line rule through
    ``_DataFieldState``: over an 80-char ``$``-led non-keyword line the splitter must frame one record
    and the parser must yield two fields.  Independent arming logic in the two fails one or the other.
    """
    from chython.formats.ctfile import parse_rdf_fields

    long_dollar = '$100 for the reagent'.ljust(80, '.')
    assert len(long_dollar) == 80

    # Build a record body as if the splitter already framed it.  The splitter sees no record tags
    # inside (the $-led line is not $MFMT/$RFMT), so framing is unambiguous.  The parser then
    # receives exactly these lines.
    record_lines = ['$DTYPE note', '$DATUM see below', long_dollar, '$DTYPE yield', '$DATUM 95']

    # Splitter: feed lines with a surrounding $MFMT wrapper -- the $-led line must not split it.
    mol_lines = (['$MFMT', 'mol', '', '',
                  '  0  0  0     0  0            999 V2000', 'M  END']
                 + record_lines)
    records = list(split_rdf_records(mol_lines))
    assert len(records) == 1, f'splitter produced {len(records)} records, expected 1'

    # Parser: the same record_lines must yield two correct fields.
    log = []
    fields = parse_rdf_fields(record_lines, log)
    assert fields == {'note': 'see below', 'yield': '95'}, \
        f'parser produced {len(fields)} fields, expected 2'


# --- A line with no open $DATUM value to join is dropped, but never silently

def test_line_between_dtype_and_datum_is_logged():
    """Text between a ``$DTYPE`` and its ``$DATUM`` has nowhere to go -- the name is on the ``$DTYPE`` line
    and no value is open -- so the line is dropped and reported, since a silent drop leaves the caller
    unable to learn the record held a line the reader could not place.
    """
    from chython.formats.ctfile import parse_rdf_fields

    log = []
    fields = parse_rdf_fields(['$DTYPE a', 'stray text', '$DATUM v'], log)
    assert fields == {'a': 'v'}, fields
    assert len(log) == 1, log
    assert 'stray text' in log[0], log[0]
    assert not str(log[0]).startswith('unsupported: '), log[0]


def test_permissive_continuation_of_an_orphan_datum_is_logged():
    """A non-``$`` line after an orphan ``$DATUM`` is a second lost line and gets a second message: the
    orphan message covers one logical line, and a permissive continuation is a different physical line with
    content of its own.
    """
    from chython.formats.ctfile import parse_rdf_fields

    log = []
    fields = parse_rdf_fields(['$DATUM v', 'more text'], log)
    assert fields == {}
    assert len(log) == 2, log
    assert 'with no preceding' in log[0], log[0]
    assert 'more text' in log[1], log[1]


def test_text_after_a_closed_logical_line_does_not_rejoin_the_previous_value():
    """Once a ``$``-led non-keyword line has ended the logical line, plain text does not resume it.  The
    field parser must read that off the shared state rather than re-derive it, or the text is appended to a
    value the file closed two lines ago -- data invented, and the splitter disagreeing.
    """
    from chython.formats.ctfile import parse_rdf_fields

    log = []
    fields = parse_rdf_fields(['$DTYPE k', '$DATUM v', '$NOTAKEYWORD', 'orphaned text'], log)
    assert fields == {'k': 'v'}, fields
    assert len(log) == 2, log
    assert 'unrecognised field keyword' in log[0], log[0]
    assert 'orphaned text' in log[1], log[1]


# --- The stream surface: RDFRead and RDFWrite

def test_rdfread_iterates_molecules_and_reactions(root):
    """An RDfile interleaves the two kinds by design, so iteration has to yield both."""
    from chython.core import MoleculeContainer
    from chython.core.reaction import ReactionContainer
    from chython.formats.ctfile import RDFRead

    with RDFRead(root / 'test' / 'MR.rdf') as f:
        objects = list(f)
        assert not f.failed, [x.error for x in f.failed]
    assert [type(x) is ReactionContainer for x in objects] == [False, False, True, True]
    assert isinstance(objects[0], MoleculeContainer)


def test_rdfread_meta_follows_the_current_record(root):
    from chython.formats.ctfile import RDFRead

    with RDFRead(root / 'test' / 'MR.rdf') as f:
        seen = []
        for _ in f:
            seen.append(f.meta.get('CdId'))
    assert seen == ['MOL V2000', 'MOL V3000', 'RXN V2000', 'RXN V3000']


#: This file builds its records inline, so these two do too.  `$DATM` on line 2 is where the format
#: puts a timestamp -- see `test_the_header_timestamp_is_stored_and_the_tail_one_is_reported`.
_METHANE = ['methane', '  test', '', '  1  0  0  0  0  0            999 V2000',
            '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0', 'M  END']
_RDF_HEADER = ['$RDFILE 1', '$DATM    01/02/17 17:17']


def test_a_molecule_record_carries_its_dtype_pairs(tmp_path):
    """A molecule record's metadata is the molecule's, not a wrapper's."""
    from chython.formats.ctfile import RDFRead

    path = tmp_path / 'a.rdf'
    path.write_text('\n'.join(_RDF_HEADER + ['$MFMT'] + _METHANE
                              + ['$DTYPE K', '$DATUM v']) + '\n')
    with RDFRead(path) as f:
        mol = next(iter(f))
        assert f.meta == mol.meta, 'the reader delegates rather than keeping a second copy'
    assert mol.meta == {'K': 'v'}


def test_a_molecule_record_round_trips_its_metadata(tmp_path):
    """``write(mol)`` with nothing else said writes the molecule's own ``$DTYPE`` pairs."""
    from chython.formats.ctfile import RDFRead, RDFWrite
    from .._sdf import parse_record

    mol = parse_record(_METHANE)
    mol.meta['K'] = 'v'
    path = tmp_path / 'b.rdf'
    with RDFWrite(path) as w:
        w.write(mol)
    with RDFRead(path) as f:
        assert next(iter(f)).meta == {'K': 'v'}


def test_all_four_rdf_fixtures_read(root):
    from chython.formats.ctfile import RDFRead

    for name in ('MR.rdf', 'ions.rdf', 'standardize.rdf', 'reaction_centerslist.rdf'):
        with RDFRead(root / 'test' / name) as f:
            objects = list(f)
        assert objects, name
        assert not f.failed, (name, [x.error for x in f.failed])


def test_read_record_raises_stopiteration_at_end_of_file(root):
    """The same end-of-file convention SDFRead has.  A ``None`` return would be a second one."""
    from pytest import raises

    from chython.formats.ctfile import RDFRead

    with RDFRead(root / 'test' / 'MR.rdf') as f:
        for _ in range(4):
            f.read_record()
        with raises(StopIteration):
            f.read_record()


def test_the_log_belongs_to_the_current_record(tmp_path):
    """``f.log`` is per-record, so a clean record must not inherit the previous record's damage.  The damage
    is a reacting-centre code in the bond block, written by hand because no fixture in ``test/`` carries one
    -- all are 0 in columns 19-21, ``reaction_centerslist.rdf`` being named for a computed CGR list instead.
    """
    from chython.formats.ctfile import RDFRead

    def record(centre):
        return ['$MFMT', 'ethane', '', '',
                '  2  1  0  0  0  0            999 V2000',
                '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
                '    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
                f'  1  2  1  0  0  0{centre:3d}',
                'M  END']

    path = tmp_path / 'centres.rdf'
    path.write_text('\n'.join(['$RDFILE 1', *record(0), *record(1)]) + '\n', encoding='utf8')

    with RDFRead(path) as f:
        f.read_record()
        first = list(f.log)
        f.read_record()
        second = list(f.log)
    assert not any('reacting-centre' in x for x in first), first
    assert any('reacting-centre' in x for x in second), second


def test_rdfwrite_round_trips_both_kinds(root, tmp_path):
    """The container alone is enough to write the record back, both kinds, metadata included."""
    from chython.core.reaction import ReactionContainer
    from chython.formats.ctfile import RDFRead, RDFWrite

    with RDFRead(root / 'test' / 'MR.rdf') as f:
        records = [f.read_record() for _ in range(4)]

    path = tmp_path / 'out.rdf'
    with RDFWrite(path) as out:
        for record in records:
            out.write(record)

    with RDFRead(path) as f:
        again = [f.read_record() for _ in range(4)]
        assert not f.failed, [x.error for x in f.failed]
    assert [isinstance(r, ReactionContainer) for r in again] == \
           [isinstance(r, ReactionContainer) for r in records]
    assert [r.meta.get('CdId') for r in again] == [r.meta.get('CdId') for r in records]


def test_the_writer_emits_its_own_header_once_and_not_at_all_when_appending(tmp_path):
    """``$RDFILE 1`` opens a file and must not turn up in the middle of one.

    The flag is the writer's own, not ``self._file.tell()``: a caller may hand this class a pipe, and
    a pipe cannot be asked where it is.
    """
    from chython.formats.ctfile import RDFWrite, mol

    molecule = mol('\n'.join(['x', '', '', '  1  0  0  0  0  0            999 V2000',
                              '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
                              'M  END']))
    path = tmp_path / 'header.rdf'
    with RDFWrite(path) as out:
        out.write(molecule)
        out.write(molecule)
    with RDFWrite(path, append=True) as out:
        out.write(molecule)

    lines = path.read_text(encoding='utf8').split('\n')
    assert [x for x in lines if x.startswith('$RDFILE')] == ['$RDFILE 1']
    assert len([x for x in lines if x.startswith('$DATM')]) == 1
    assert len([x for x in lines if x.startswith('$MFMT')]) == 3


def test_write_accepts_an_iterable_of_pairs_as_well_as_a_mapping(root, tmp_path):
    """RETARGETED: ``record.fields`` was ``[DataField, ...]`` and is gone, but the writer still takes
    any iterable of pairs -- a caller that built its fields in order has nothing to convert.
    """
    from chython.formats.ctfile import RDFRead, RDFWrite

    with RDFRead(root / 'test' / 'MR.rdf') as f:
        molecule = f.read_record()
    path = tmp_path / 'fields.rdf'
    with RDFWrite(path) as out:
        out.write(molecule, meta=list(molecule.meta.items()))
    with RDFRead(path) as f:
        assert f.read_record().meta == molecule.meta


def test_a_long_meta_value_is_wrapped_and_survives_the_next_record(tmp_path):
    """The writer side of the 80-column continuation rule: no physical line it emits reaches 80 characters
    unless it is a real continuation, since an unwrapped 200-character ``$DATUM`` opens a continuation on
    re-read and swallows the ``$MFMT`` after it.
    """
    from chython.formats.ctfile import RDFRead, RDFWrite, mol

    molecule = mol('\n'.join(['x', '', '', '  1  0  0  0  0  0            999 V2000',
                              '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
                              'M  END']))
    long = 'A' * 200
    path = tmp_path / 'long.rdf'
    with RDFWrite(path) as out:
        out.write(molecule, meta={'note': long})
        out.write(molecule, meta={'note': 'short'})

    with RDFRead(path) as f:
        records = [f.read_record() for _ in range(2)]
        assert not f.failed, [x.error for x in f.failed]
    assert records[0].meta['note'] == long
    assert records[1].meta['note'] == 'short'


def test_a_long_meta_name_is_wrapped_too(tmp_path):
    """The reader continues a wrapped ``$DTYPE`` name, so one rule covers every logical line."""
    from chython.formats.ctfile import RDFRead, RDFWrite, mol

    molecule = mol('\n'.join(['x', '', '', '  1  0  0  0  0  0            999 V2000',
                              '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
                              'M  END']))
    name = 'N' * 150
    path = tmp_path / 'name.rdf'
    with RDFWrite(path) as out:
        out.write(molecule, meta={name: 'v'})
        out.write(molecule, meta={'after': 'w'})

    with RDFRead(path) as f:
        records = [f.read_record() for _ in range(2)]
        assert not f.failed, [x.error for x in f.failed]
    assert records[0].meta[name] == 'v'
    assert records[1].meta['after'] == 'w'


def test_a_value_exactly_filling_the_datum_column_round_trips(tmp_path):
    """A value whose ``$DATUM`` line reaches exactly 80 columns round-trips byte for byte.  Such a line
    arms positional continuation, so the writer appends an empty physical line: the reader absorbs it as
    the continuation, adds nothing, disarms, and sees the next keyword normally.  Three lengths --
    ``n=73``, ``n=153``, ``n=233`` -- because each wrap boundary can break for its own reason.
    """
    from io import StringIO

    from chython.formats.ctfile import RDFRead, RDFWrite, mol

    mol_text = '\n'.join(['x', '', '', '  1  0  0  0  0  0            999 V2000',
                          '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
                          'M  END'])
    molecule = mol(mol_text)

    for n in (73, 153, 233):
        v = 'A' * n
        assert len('$DATUM ' + v) % 80 == 0, 'not a whole-multiple length'
        buf = StringIO()
        with RDFWrite(buf) as out:
            log = out.write(molecule, meta={'A': v, 'B': 'sentinel'})
        assert not any('multiple of 80' in x for x in log), (n, log)
        with RDFRead(StringIO(buf.getvalue())) as f:
            got = dict(f.read_record().meta)
        assert got.get('A') == v, (n, repr(got.get('A'))[:40])
        assert got.get('B') == 'sentinel', (n, got)


def test_a_value_with_trailing_space_at_80_columns_does_not_arm_empty_line(tmp_path):
    """A $DATUM line of 80 raw characters ending in spaces does not arm continuation: the rule measures
    content without trailing whitespace, so ``_continues`` is ``False`` and the empty-line append gate
    mirrors it.  The value reads back stripped and the following field survives.
    """
    from io import StringIO

    from chython.formats.ctfile import RDFRead, RDFWrite, mol

    mol_text = '\n'.join(['x', '', '', '  1  0  0  0  0  0            999 V2000',
                          '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
                          'M  END'])
    molecule = mol(mol_text)

    # '$DATUM ' + 'A'*72 + ' ' is 80 raw, 79 rstripped -> does not arm
    v = 'A' * 72 + ' '
    assert len('$DATUM ' + v) == 80
    assert len(('$DATUM ' + v).rstrip()) == 79
    buf = StringIO()
    with RDFWrite(buf) as out:
        log = out.write(molecule, meta={'key': v, 'sentinel': 'after'})
    # No log about wrapping at a space for the final chunk
    assert not any('wraps at a space' in x for x in log), log
    with RDFRead(StringIO(buf.getvalue())) as f:
        got = dict(f.read_record().meta)
    # The value should be the rstripped version (trailing space removed by reader)
    assert got.get('key') == 'A' * 72, repr(got.get('key'))
    # The sentinel field must survive -- not swallowed by the empty line
    assert got.get('sentinel') == 'after', got


def test_values_of_all_length_classes_round_trip(tmp_path):
    """Values from 1 to 400 characters all survive write/read for ``$DATUM`` and ``$DTYPE`` -- one chunk
    (< 73), whole multiples of 80, and everything else above 80.  The following sentinel field is checked
    too, the silent failure being the sentinel swallowed.
    """
    from io import StringIO

    from chython.formats.ctfile import RDFRead, RDFWrite, mol

    mol_text = '\n'.join(['x', '', '', '  1  0  0  0  0  0            999 V2000',
                          '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
                          'M  END'])
    molecule = mol(mol_text)

    # $DATUM sweep: value of length n
    bad_datum = []
    for n in range(1, 401):
        v = 'y' * n
        buf = StringIO()
        with RDFWrite(buf) as out:
            out.write(molecule, meta={'A': v, 'B': 'sentinel'})
        with RDFRead(StringIO(buf.getvalue())) as f:
            got = dict(f.read_record().meta)
        if got.get('A') != v or got.get('B') != 'sentinel':
            bad_datum.append(n)

    # $DTYPE sweep: name of length n
    bad_dtype = []
    for n in range(1, 401):
        name = 'N' * n
        buf = StringIO()
        with RDFWrite(buf) as out:
            out.write(molecule, meta={name: 'value', 'after': 'sentinel'})
        with RDFRead(StringIO(buf.getvalue())) as f:
            got = dict(f.read_record().meta)
        if got.get(name) != 'value' or got.get('after') != 'sentinel':
            bad_dtype.append(n)

    # Multiline values at several boundary lengths: a genuine newline and a wrapped long line
    # use the same rejoin mechanism and must not be confused.
    bad_multiline = []
    for n in (72, 73, 152, 153, 232, 233):
        v = 'y' * n + '\nmore text'
        buf = StringIO()
        with RDFWrite(buf) as out:
            out.write(molecule, meta={'A': v, 'B': 'sentinel'})
        with RDFRead(StringIO(buf.getvalue())) as f:
            got = dict(f.read_record().meta)
        if got.get('A') != v or got.get('B') != 'sentinel':
            bad_multiline.append(n)

    assert not bad_datum, f'$DATUM lengths that failed: {bad_datum}'
    assert not bad_dtype, f'$DTYPE lengths that failed: {bad_dtype}'
    assert not bad_multiline, f'multiline lengths that failed: {bad_multiline}'


def test_a_multi_line_meta_value_keeps_its_lines(tmp_path):
    from chython.formats.ctfile import RDFRead, RDFWrite, mol

    molecule = mol('\n'.join(['x', '', '', '  1  0  0  0  0  0            999 V2000',
                              '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
                              'M  END']))
    path = tmp_path / 'multi.rdf'
    with RDFWrite(path) as out:
        out.write(molecule, meta={'note': 'one\ntwo\nthree'})
    with RDFRead(path) as f:
        assert f.read_record().meta['note'] == 'one\ntwo\nthree'


def test_the_v3000_writer_changes_the_ctab_and_not_the_rdfile_framing(root, tmp_path):
    """``$RDFILE``/``$MFMT``/``$DTYPE`` are not versioned; only the CTAB the record holds is."""
    from chython.formats.ctfile import ERDFWrite, RDFRead

    with RDFRead(root / 'test' / 'MR.rdf') as f:
        molecule = f.read_record()
    path = tmp_path / 'v3000.rdf'
    with ERDFWrite(path) as out:
        out.write(molecule)

    text = path.read_text(encoding='utf8')
    assert 'M  V30 BEGIN CTAB' in text
    assert text.startswith('$RDFILE 1\n')
    with RDFRead(path) as f:
        again = f.read_record()
        assert f.version == 'V3000'
    assert again.meta == molecule.meta


def test_a_reaction_keeps_its_metadata_after_the_record_is_gone(root):
    """The container is the durable home, so iteration is not a lossy path for a reaction:
    ``ReactionContainer.meta`` is where an RDfile's ``$DTYPE``/``$DATUM`` pairs live, and a caller who keeps
    the container and drops the record must not lose them.
    """
    from chython.core.reaction import ReactionContainer
    from chython.formats.ctfile import RDFRead

    with RDFRead(root / 'test' / 'MR.rdf') as f:
        reactions = [x for x in f if isinstance(x, ReactionContainer)]
    assert [r.meta.get('CdId') for r in reactions] == ['RXN V2000', 'RXN V3000']


def test_a_repeated_dtype_merges_on_the_container_and_is_reported(tmp_path):
    """INVERTED: this asserted a repeated ``$DTYPE`` kept both values on ``record.fields`` and
    collapsed last-wins on the container.  There is one storage now, so there is no lossless second
    place to keep them -- the values merge and the collision is a log line on the container, which is
    what leaves the caller able to learn the file said the name twice.
    """
    from chython.formats.ctfile import RDFRead

    path = tmp_path / 'dup.rdf'
    path.write_text('\n'.join(['$RDFILE 1', '$RFMT', '$RXN', 'name', '', '', '  0  0',
                               '$DTYPE k', '$DATUM first',
                               '$DTYPE k', '$DATUM second']) + '\n', encoding='utf8')
    with RDFRead(path) as f:
        reaction = f.read_record()
    assert reaction.meta == {'k': 'first\nsecond'}
    assert any('appears twice' in x for x in reaction.log), reaction.log


def test_pach_refuses_a_reaction_carrying_metadata_and_takes_the_waiver(root):
    """``pach`` has no metadata field and refuses rather than dropping silently, so read-then-pach raises
    until the caller waives it.  The entry point for a reaction is ``reaction_pach_dump``; the bare
    ``pach_dump`` beside it takes a *molecule*.  ``title`` is waived alongside ``meta`` because an
    RDfile's reaction carries a name line too, and that refusal is a different one.
    """
    from pytest import raises

    from chython.core.reaction import ReactionContainer, reaction_pach_dump
    from chython.formats.ctfile import RDFRead

    with RDFRead(root / 'test' / 'MR.rdf') as f:
        reaction = next(x for x in f if isinstance(x, ReactionContainer))
    assert reaction.meta
    with raises(ValueError, match='metadata key'):
        reaction_pach_dump(reaction, drop=['title'])
    assert reaction_pach_dump(reaction, drop=['meta', 'title'])


def test_the_file_timestamp_is_on_the_reader_and_not_in_any_record(root):
    """One ``$DATM`` describes the whole file, so it is reader state and not record metadata, kept verbatim:
    the payload looks like ``01/02/17 17:17`` and whether that is 2017 or 1917 is not the reader's business.
    """
    from chython.formats.ctfile import RDFRead

    with RDFRead(root / 'test' / 'MR.rdf') as f:
        record = f.read_record()
        assert f.date == '01/02/17 17:17'
        assert not any('$DATM' in x for x in f.log)
        assert 'date' not in record.meta


def test_a_timestamp_outside_the_header_is_reported_as_a_broken_file(tmp_path):
    """The other half: a valid position stores, an invalid one reports -- and unprefixed."""
    from chython.formats.ctfile import parse_rdf_fields

    log = []
    parse_rdf_fields(['$DTYPE a', '$DATUM v', '$DATM    01/02/17 17:17'], log)
    assert any('$DATM' in x for x in log), log
    assert not any(str(x).startswith('unsupported: ') for x in log), log


def test_an_unparsable_record_costs_the_record_and_not_the_file(tmp_path):
    """``failed`` is the third option: not a raise that loses the file, not a silent skip."""
    from chython.formats.ctfile import RDFRead

    good = ['$MFMT', 'ethane', '', '',
            '  1  0  0  0  0  0            999 V2000',
            '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
            'M  END']
    bad = ['$RFMT', '$RXN', 'broken', '', '']   # a $RXN with no counts line at all
    path = tmp_path / 'mixed.rdf'
    path.write_text('\n'.join(['$RDFILE 1', *bad, *good]) + '\n', encoding='utf8')

    with RDFRead(path) as f:
        objects = list(f)
    assert len(objects) == 1
    assert len(f.failed) == 1
    assert f.failed[0].position == 0
    assert f.failed[0].lines[0] == '$RXN'


def test_sdfread_survives_an_rxn_block_and_still_reads_its_data_fields(tmp_path):
    """An SD file is not supposed to hold one, and input is garbage by default, so we read it.  Iteration is
    tested beside ``read_record()`` because a reaction record answers different questions from a molecule one.
    """
    from chython.formats.ctfile import SDFRead
    from chython.core.reaction import ReactionContainer

    component = ['one', '', '', '  1  0  0  0  0  0            999 V2000',
                 '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
                 'M  END']
    text = '\n'.join(['$RXN', 'title', '', '', '  1  1',
                      '$MOL', *component, '$MOL', *component,
                      '>  <ID>', '7', '', '$$$$']) + '\n'
    path = tmp_path / 'in.sdf'
    path.write_text(text, encoding='utf8')

    with SDFRead(path) as f:
        reaction = f.read_record()
        assert not f.failed, [x.error for x in f.failed]
        assert isinstance(reaction, ReactionContainer)
        assert any('$RXN' in x for x in reaction.log), reaction.log
        # the SDF data fields after the last M END are still read, and nothing was logged as stray
        assert reaction.meta == {'ID': '7'}
        assert not any('outside any field' in x for x in reaction.log), reaction.log
        # the properties a caller reaches for do not raise on a reaction record
        assert f.sgroups is None and f.unknown_hydrogens == ()

    with SDFRead(path) as f:
        assert [type(x) for x in f] == [ReactionContainer]


def test_framing_damage_reaches_the_caller_and_is_not_a_record_field(tmp_path):
    """The splitter's log has to go somewhere, and a record's own log is not that somewhere: a registry
    reference on a tag line, a stray line before the first record and a suppressed delimiter are all
    decisions about where records begin and end.  Handing the splitter no log loses every one of them.
    """
    from chython.formats.ctfile import RDFRead

    path = tmp_path / 'framing.rdf'
    path.write_text('\n'.join(['$RDFILE 1', 'stray header line',
                               '$MFMT reg-99', 'ethane', '', '',
                               '  1  0  0  0  0  0            999 V2000',
                               '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
                               'M  END']) + '\n', encoding='utf8')
    with RDFRead(path) as f:
        molecule = f.read_record()
    assert any('before the first record' in x for x in f.file_log), f.file_log
    assert any('registry' in x for x in f.file_log), f.file_log
    # and none of it is attributed to the record, which did not cause any of it
    assert not any('registry' in x for x in molecule.log), molecule.log


# --- Fix 2 write half: meta=None falls back to the reaction container's own meta

def test_rdfwrite_uses_reaction_meta_when_none_is_passed(root, tmp_path):
    """``meta=None`` on a reaction write means "use the container's own metadata".

    ``meta={}`` is still how a caller writes a record with no metadata regardless of what the
    container holds.  The asymmetry with molecules (which have no ``.meta``) is a storage fact
    documented in the ``write`` docstring.
    """
    from chython.core.reaction import ReactionContainer
    from chython.formats.ctfile import RDFRead, RDFWrite

    with RDFRead(root / 'test' / 'MR.rdf') as f:
        reactions = [x for x in f if isinstance(x, ReactionContainer)]
    assert all(r.meta for r in reactions), 'fixture reactions must carry metadata'

    path = tmp_path / 'reactions.rdf'
    with RDFWrite(path) as out:
        for rxn in reactions:
            out.write(rxn)          # meta=None -- should use rxn.meta

    with RDFRead(path) as f:
        back = [x for x in f if isinstance(x, ReactionContainer)]
    assert not f.failed, [x.error for x in f.failed]
    assert [r.meta.get('CdId') for r in back] == [r.meta.get('CdId') for r in reactions]


# --- Fix 3: record tag is written only after the emitter succeeds

def test_write_refusal_leaves_no_tag_in_the_file(tmp_path):
    """A writer that refuses must leave the file as it found it: no ``$MFMT``/``$RFMT`` tag with no record
    body under it.  The refusal here is a coordinate too wide for the V2000 atom line's 10-character
    column -- about the file format, not the chemistry -- but the framing is what is under test, so any
    refusal will do.
    """
    from io import StringIO

    from pytest import raises

    from chython.formats.ctfile import RDFWrite
    from chython.formats.ctfile._errors import MalformedCtfile
    from chython.formats.ctfile._sdf import parse_record

    mol = parse_record(['methane', '', '',
                        '  1  0  0  0  0  0            999 V2000',
                        '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
                        'M  END'], [])
    with mol.edit():
        mol.set_xy(next(iter(mol.atom_numbers)), 123456., 0.)

    buf = StringIO()
    w = RDFWrite(buf)
    with raises(MalformedCtfile, match='10-character column'):
        w.write(mol)
    w.close()

    text = buf.getvalue()
    # The file header may be present (a file with a header and no records is a valid empty RDfile),
    # but no record tag must appear under it.
    assert '$MFMT' not in text, repr(text[:80])
    assert '$RFMT' not in text, repr(text[:80])


# --- Fix 4: $RXN block under $MFMT is rescued in RDFRead, not filed as a failure

def test_mfmt_holding_a_rxn_block_is_rescued_and_logged(tmp_path):
    """A ``$MFMT`` record whose body is a ``$RXN`` block is read as a reaction, not a failure.

    The file is wrong (``$RFMT`` is the correct tag), but reading it beats refusing it.  The log
    line is unprefixed -- the file is broken, not a construct chython declines to model.
    """
    from chython.core.reaction import ReactionContainer
    from chython.formats.ctfile import RDFRead

    component = ['mol', '', '', '  1  0  0  0  0  0            999 V2000',
                 '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
                 'M  END']
    text = '\n'.join(['$RDFILE 1', '$MFMT', '$RXN', 'title', '', '', '  1  1',
                      '$MOL', *component, '$MOL', *component]) + '\n'
    path = tmp_path / 'rxn_in_mfmt.rdf'
    path.write_text(text, encoding='utf8')

    with RDFRead(path) as f:
        reaction = f.read_record()
        assert not f.failed, [x.error for x in f.failed]
    assert isinstance(reaction, ReactionContainer)
    assert any('$RXN' in x for x in reaction.log), reaction.log
    assert not any(str(x).startswith('unsupported: ') for x in reaction.log), reaction.log


# --- Fix 5: data-field cut is after the last M  END, not at the first >

def test_component_title_starting_with_gt_does_not_cut_the_reaction(tmp_path):
    """A ``$MOL`` component whose title starts with ``>`` must not be mistaken for the data delimiter: the
    SDF delimiter is a column-0 ``>`` and so is a title like ``> product name``, so the data block is found
    after the last ``M  END`` rather than from the top of the record.
    """
    from chython.core.reaction import ReactionContainer
    from chython.formats.ctfile import SDFRead

    component_with_gt = ['> product name', '', '', '  1  0  0  0  0  0            999 V2000',
                         '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
                         'M  END']
    reactant = ['reactant', '', '', '  1  0  0  0  0  0            999 V2000',
                '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
                'M  END']
    text = '\n'.join(['$RXN', 'title', '', '', '  1  1',
                      '$MOL', *reactant, '$MOL', *component_with_gt,
                      '>  <ID>', '7', '', '$$$$']) + '\n'
    path = tmp_path / 'gt_title.sdf'
    path.write_text(text, encoding='utf8')

    with SDFRead(path) as f:
        reaction = f.read_record()
        assert not f.failed, [x.error for x in f.failed]
    assert isinstance(reaction, ReactionContainer)
    assert len(list(reaction.reactants)) == 1
    assert len(list(reaction.products)) == 1
    assert reaction.meta == {'ID': '7'}


# --- Fix 6: append=True on a non-existent path still writes the file header

def test_append_to_a_fresh_path_writes_the_header(tmp_path):
    """Opening a new file with ``append=True`` must still write the ``$RDFILE 1`` header: a fresh file has
    none, and our own reader accepts headerless files while a third-party one may not.
    """
    from chython.formats.ctfile import RDFRead, RDFWrite, mol

    molecule = mol('\n'.join(['x', '', '', '  1  0  0  0  0  0            999 V2000',
                              '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
                              'M  END']))
    path = tmp_path / 'fresh.rdf'
    with RDFWrite(path, append=True) as out:
        out.write(molecule, meta={'a': 'b'})
    text = path.read_text(encoding='utf8')
    assert text.startswith('$RDFILE 1\n'), repr(text[:40])

    # A real append to an existing file must NOT add a second header.
    with RDFWrite(path, append=True) as out:
        out.write(molecule, meta={'c': 'd'})
    text = path.read_text(encoding='utf8')
    assert text.count('$RDFILE 1') == 1, repr(text[:80])
    with RDFRead(path) as f:
        objects = list(f)
    assert len(objects) == 2


# --- Fix 7: leading/trailing whitespace on name or value is logged as unsupported

def test_write_logs_whitespace_on_name_and_value(tmp_path):
    """Whitespace on a ``$DTYPE`` name or ``$DATUM`` value that the reader will strip is logged, or the
    caller's data is truncated silently.  ``unsupported: `` because the format cannot carry the construct.
    """
    from chython.formats.ctfile import RDFRead, RDFWrite, mol

    molecule = mol('\n'.join(['x', '', '', '  1  0  0  0  0  0            999 V2000',
                              '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
                              'M  END']))
    path = tmp_path / 'ws.rdf'
    with RDFWrite(path) as out:
        log_name = out.write(molecule, meta={' N ': 'v'})
        log_val = out.write(molecule, meta={'A': '  padded  '})
    # Both cases produce an unsupported: log line
    assert any(str(x).startswith('unsupported: ') and 'name' in x for x in log_name), log_name
    assert any(str(x).startswith('unsupported: ') and 'value' in x for x in log_val), log_val
    # The value reads back stripped, so the caller knows from the log what happened
    with RDFRead(path) as f:
        r1 = f.read_record()
        r2 = f.read_record()
    assert r1.meta.get('N') == 'v', r1.meta
    assert r2.meta.get('A') == 'padded', r2.meta


def _meta_round_trip(value, tmp_path, name='K'):
    """``(write log, value read back)`` for one metadata field through a real file."""
    from chython.formats.ctfile import RDFRead, RDFWrite, mol

    molecule = mol('\n'.join(['x', '', '', '  1  0  0  0  0  0            999 V2000',
                              '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
                              'M  END']))
    path = tmp_path / 'ws.rdf'
    with RDFWrite(path) as out:
        log = out.write(molecule, meta={name: value})
    with RDFRead(path) as f:
        record = f.read_record()
    return log, record.meta.get(name)


def test_write_does_not_claim_a_loss_a_wrapped_value_does_not_suffer(tmp_path):
    """A long value's trailing whitespace survives, so nothing is reported.  The reader strips only the
    first physical chunk of a logical line and concatenates positional continuations raw, so a value long
    enough to wrap carries its tail on a continuation and comes back byte for byte -- and `unsupported: `
    is the prefix a caller screens on, so a false hit there costs more than most wrong log lines.
    """
    value = 'z' * 73 + '   '
    log, read_back = _meta_round_trip(value, tmp_path)

    assert read_back == value, 'a wrapped value keeps its trailing whitespace'
    assert not [x for x in log if 'whitespace' in x], log


def test_write_reports_edge_whitespace_exactly_at_the_73_column_boundary(tmp_path):
    """``$DATUM `` occupies 7 of the 80 columns the continuation rule measures, so 73 characters is the
    longest value entirely inside the chunk the reader strips -- the last width that loses.  Both sides of
    the boundary are asserted, since one side alone passes for an off-by-one.
    """
    lost, lost_back = _meta_round_trip('z' * 70 + '   ', tmp_path)
    kept, kept_back = _meta_round_trip('z' * 73 + '   ', tmp_path)

    assert lost_back == 'z' * 70, 'a 73-column value is stripped whole'
    assert any(str(x).startswith('unsupported: ') and 'trailing whitespace' in x for x in lost), lost

    assert kept_back == 'z' * 73 + '   '
    assert not [x for x in kept if 'whitespace' in x], kept


def test_a_permissive_continuation_keeps_both_edges_and_says_nothing(tmp_path):
    """A later line of a multi-line value is appended raw, so neither edge is touched: the reader stores a
    permissive continuation as its own element and does not strip it, so whitespace on it is carried and
    reporting it would be the same false claim as the wrapped case above.
    """
    value = 'first\n  second  '
    log, read_back = _meta_round_trip(value, tmp_path)

    assert read_back == value
    assert not [x for x in log if 'whitespace' in x], log


def test_a_wrap_landing_inside_whitespace_is_still_reported(tmp_path):
    """A value past 73 characters ending in whitespace is not always safe: with the 80-column cut inside
    that whitespace the first chunk stops short of 80, arms no continuation, and the remainder is rejoined
    with a ``\\n`` in it.  A different mechanism owns this -- `_wrap_logical_line` predicts the newline --
    so suppressing the `unsupported: ` line above does not suppress it.
    """
    value = 'z' * 71 + '   '   # 74 characters, so `$DATUM ` + value cuts at 80 inside the spaces
    log, read_back = _meta_round_trip(value, tmp_path)

    assert read_back != value, 'this case really does lose the value'
    assert any('wraps at a space' in x and 'newline' in x for x in log), log
    assert not [x for x in log if str(x).startswith('unsupported: ')], 'a mangling is not a missing feature'


# --- chython 2 differential: side assignment and component order

def _drain(reader):
    """Every record a reader will give up; `read_record()` signals end of file by raising.  Written as a
    loop rather than a comprehension because PEP 479 would turn that `StopIteration` into a `RuntimeError`
    exactly when the oracle and this reader disagree on the record count.
    """
    out = []
    while True:
        try:
            out.append(reader.read_record())
        except StopIteration:
            return out


def test_v2_agrees_on_the_sides_of_every_rdf_fixture(root, oracle_session):
    """chython 2 reading the same four files, asked for *counts* rather than molecules: side assignment and
    component order are what this stream can get wrong in a way no hand-written test catches.

    `agents` is compared and no fixture exercises it -- all eleven reactions have zero agents and RXN
    V2000's counts line has no agent field -- so that comparison is `0 == 0` everywhere.  The final
    assertion names files the oracle could not read, since `tolerant=True` must not silently shrink
    coverage.
    """
    from chython.core.reaction import ReactionContainer
    from chython.formats.ctfile import RDFRead

    names = ('MR.rdf', 'ions.rdf', 'standardize.rdf', 'reaction_centerslist.rdf')
    expected = oracle_session.read_rdf({n: root / 'test' / n for n in names}, tolerant=True)

    unread = []
    for name in names:
        reference = expected.get(name)
        if reference is None:
            unread.append(name)
            continue      # V2 cannot read the file; collect and report below
        with RDFRead(root / 'test' / name) as f:
            ours = [{'reactants': len(r.reactants),
                     'products': len(r.products),
                     'agents': len(r.agents),
                     'atoms': [m.atom_count for m in r.molecules()]}
                    for r in _drain(f) if isinstance(r, ReactionContainer)]
        theirs = [x for x in reference if x is not None]
        assert ours == theirs, name
    # every fixture must have been readable by V2; a missing one is a broken oracle or a moved file
    assert not unread, f'chython 2 could not read {unread}'


def test_a_value_split_before_a_dollar_round_trips_and_the_reader_says_it_absorbed_a_keyword():
    """The one round trip that is exact *and* logged.  The writer chunks at the arming length, so a boundary
    can fall just before a ``$`` inside a value; the reader rejoins it byte for byte and still reports
    absorbing something keyword-shaped, since it cannot know who wrote the file."""
    from chython.formats.ctfile._rdf import _wrap_logical_line, parse_rdf_fields

    value = 'a' * 73 + '$xyz'
    writer_log = []
    physical = _wrap_logical_line('$DATUM ' + value, writer_log, 'a $DATUM value')
    assert [len(x) for x in physical] == [80, 4]
    assert physical[1].startswith('$')
    assert not writer_log        # the writer had nothing to decide; the chunking is what it is

    reader_log = []
    fields = parse_rdf_fields(['$DTYPE NAME', *physical], reader_log)
    assert fields['NAME'] == value
    assert any('absorbed a keyword' in x for x in reader_log), reader_log
