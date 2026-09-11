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
"""Framing, version sniffing and data fields.

Version is sniffed per record, not per file: ``test/implicit.sdf`` mixes V2000 and V3000 records, and
a per-file choice reads the odd one out as a molecule with zero atoms and no error.
"""

from pytest import raises

from chython.core import MoleculeContainer
from .._errors import MalformedCtfile
from .._sdf import (RECORD_SEPARATOR, UNPARSED_KEY, emit_record, parse_data_fields, parse_record,
                    sniff_version, split_records)
from .._sgroup import UNSUPPORTED
from .._v2000 import V2000_STAMP
from .._v3000 import V3000_STAMP


_MINIMAL = ['methane', '  test', '', '  1  0  0  0  0  0            999 V2000',
            '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0', 'M  END']


# framing

def test_records_are_split_on_the_separator():
    lines = _MINIMAL + [RECORD_SEPARATOR] + _MINIMAL + [RECORD_SEPARATOR]
    assert len(list(split_records(lines))) == 2


def test_a_final_record_without_a_separator_is_still_a_record():
    """Single-record exports omit the separator constantly.  Dropping the record loses the file."""
    assert len(list(split_records(_MINIMAL))) == 1


def test_the_separator_is_matched_as_a_prefix():
    """Real writers pad it, and a file edited across platforms leaves a carriage return on it."""
    for spelling in ('$$$$', '$$$$   ', '$$$$\r', '$$$$ 1'):
        assert len(list(split_records(_MINIMAL + [spelling] + _MINIMAL))) == 2, spelling


def test_trailing_blank_lines_are_not_a_record():
    """A file ending in a separator plus newlines would otherwise yield an empty final record."""
    assert len(list(split_records(_MINIMAL + [RECORD_SEPARATOR, '', '  ', '']))) == 1


def test_carriage_returns_are_stripped_from_every_line():
    """A stray ``\\r`` left by a CRLF file breaks every integer parse in a fixed-column field."""
    records = list(split_records([x + '\r' for x in _MINIMAL]))
    assert not any(x.endswith('\r') for x in records[0]), records[0]


def test_separator_text_inside_a_data_value_does_not_split_the_record():
    """The test is a prefix, not a substring: a value containing the characters is not a separator."""
    lines = _MINIMAL + ['>  <NOTE>', 'money$$$$money', '', RECORD_SEPARATOR]
    records = list(split_records(lines))
    assert len(records) == 1
    assert 'money$$$$money' in records[0]


# sniffing

def test_the_stamp_is_read_from_the_counts_line():
    assert sniff_version(_MINIMAL, []) == V2000_STAMP


def test_a_v3000_body_outranks_a_v2000_stamp():
    """Shipped writers produce this.  The body is the thing that can actually be parsed, so it wins."""
    log = []
    lines = list(_MINIMAL)
    lines.insert(-1, 'M  V30 BEGIN CTAB')  # inside the record, ahead of `M  END`
    assert sniff_version(lines, log) == V3000_STAMP
    assert any('stamped' in x for x in log), log


def test_a_stamp_one_column_off_is_still_read():
    """The single most common malformation in a hand-edited file.  Refusing it reads nothing."""
    lines = list(_MINIMAL)
    lines[3] = '  1  0  0  0  0  0            999  V3000'
    assert sniff_version(lines, []) == V3000_STAMP


def test_no_stamp_at_all_reads_as_v2000_and_says_so():
    """V3000 cannot express a CTAB without its own keyword lines, so their absence settles it."""
    log = []
    lines = list(_MINIMAL)
    lines[3] = '  1  0  0  0  0  0            999'
    assert sniff_version(lines, log) == V2000_STAMP
    assert any('no version stamp' in x for x in log), log


def test_the_m_v30_search_stops_at_m_end():
    """A data field mentioning ``M  V30`` after the record ends must not change the version."""
    lines = _MINIMAL + ['>  <COMMENT>', 'M  V30 BEGIN CTAB', '']
    assert sniff_version(lines, []) == V2000_STAMP


def test_a_mixed_version_file_is_sniffed_per_record(root):
    """The regression this module exists for, on the real file that motivated it."""
    path = root / 'test' / 'implicit.sdf'
    if not path.exists():
        return
    with path.open(encoding='utf8', errors='replace') as f:
        records = list(split_records(f))
    versions = {sniff_version(r, []) for r in records}
    assert V3000_STAMP in versions and V2000_STAMP in versions, (
        f'implicit.sdf should mix versions, got {versions}; if the file changed, this test is no '
        f'longer evidence and a synthetic mixed file should replace it')


# data fields

def test_data_fields_are_a_dict():
    lines = _MINIMAL + ['>  <NAME>', 'ethanol', '', '>  <MP>', '-114', '']
    assert parse_data_fields(lines) == {'NAME': 'ethanol', 'MP': '-114'}


def test_a_multi_line_value_keeps_its_lines():
    """An IUPAC name wrapped at 80 columns is one value on three lines; the line split is data."""
    assert parse_data_fields(_MINIMAL + ['>  <NOTE>', 'first', 'second', '']) == \
        {'NOTE': 'first\nsecond'}


def test_a_repeated_name_merges_and_says_so():
    """INVERTED: this asserted two `DataField`s survived a repeated name.  One value per name is what a
    mapping can spell, so the values merge, and the merge is reported rather than left to be
    discovered."""
    log = []
    lines = _MINIMAL + ['>  <K>', 'a', '', '>  <K>', 'b', '']
    assert parse_data_fields(lines, log) == {'K': 'a\nb'}
    assert any('appears twice' in x for x in log), log


def test_a_field_stated_with_no_value_is_not_the_same_as_a_field_never_stated():
    """``>  <TAG>`` then a blank line is a field stated with an empty value; an absent field is a
    different fact, and both must survive a round trip.  Do not "simplify" the blank-line branch in
    :func:`parse_data_fields` into ``if not line.strip(): break``."""
    empty = parse_data_fields(_MINIMAL + ['>  <EMPTY>', '', '>  <FULL>', 'x', ''])
    assert empty == {'EMPTY': '', 'FULL': 'x'}

    absent = parse_data_fields(_MINIMAL + ['>  <FULL>', 'x', ''])
    assert absent == {'FULL': 'x'}, 'no field conjured for the gap'

    # and the writer must not turn the stated-but-empty one into the absent one
    lines, _ = emit_record(parse_record(_MINIMAL, []), meta=empty)
    assert parse_data_fields(lines) == {'EMPTY': '', 'FULL': 'x'}


def test_a_field_number_in_the_header_is_reported():
    """INVERTED: nothing keeps the verbatim header.  There is no model for a field or a registry
    number, so neither is written back, and a header carrying one says `unsupported: `."""
    log = []
    lines = _MINIMAL + ['> 25 <MELTING.POINT> DT12 42', '180', '']
    assert parse_data_fields(lines, log) == {'MELTING.POINT': '180'}
    assert any(str(x).startswith(UNSUPPORTED) for x in log), log


def test_a_header_with_no_angle_brackets_keeps_the_field_under_its_header_text():
    log = []
    assert parse_data_fields(_MINIMAL + ['> DT1', 'value', ''], log) == {'DT1': 'value'}
    assert any('no <name>' in x for x in log), log


def test_two_blank_lines_between_fields_do_not_end_the_block():
    fields = parse_data_fields(_MINIMAL + ['>  <A>', '1', '', '', '>  <B>', '2', ''])
    assert list(fields) == ['A', 'B']


def test_data_before_m_end_is_not_a_field():
    """The block starts after ``M  END``.  A property line is not a data field."""
    assert parse_data_fields(['>  <A>', '1', ''] + _MINIMAL) == {}


def test_a_line_before_any_field_is_stored_under_the_unparsed_key():
    """The bucket exists because the input posture is store and log, never drop."""
    log = []
    lines = _MINIMAL + ['stray', '>  <K>', 'a', '']
    assert parse_data_fields(lines, log) == {UNPARSED_KEY: 'stray', 'K': 'a'}
    assert any('outside any field' in x for x in log), log

# whole records


def test_parse_record_returns_a_molecule_carrying_its_meta():
    mol = parse_record(_MINIMAL + ['>  <NAME>', 'methane', '', '$$$$'])
    assert isinstance(mol, MoleculeContainer) and mol.meta == {'NAME': 'methane'}


def test_parse_record_fills_a_header_dict():
    header = {}
    molecule = parse_record(_MINIMAL, header=header)
    assert set(header) == {'version', 'program', 'comment', 'sgroups', 'unknown_hydrogens'}
    assert molecule.title == 'methane', 'the name line is the molecule\'s, not the header dict\'s'


def test_parse_record_dispatches_on_the_sniffed_version():
    header = {}
    mol = parse_record(_MINIMAL + ['>  <ID>', '7', ''], [], header=header)
    assert [mol.element_of(s) for s in mol.atom_numbers] == [6]
    assert mol.meta == {'ID': '7'}
    assert not header['unknown_hydrogens']
    assert header['version'] == 'V2000'
    assert mol.title == 'methane'


def test_the_molecules_meta_is_one_dict_and_not_a_view_rebuilt_per_read():
    """INVERTED: this asserted the `FieldsView` object was identity-stable over a list of
    `DataField`s.  There is no second storage to view now -- `mol.meta` IS the dict -- and a caller
    holding it across an edit still has the live one."""
    mol = parse_record(_MINIMAL + ['>  <ID>', '7', ''])
    held = mol.meta
    assert mol.meta is held
    held['NEW'] = 'added'
    assert mol.meta == {'ID': '7', 'NEW': 'added'}


def test_a_write_to_meta_is_what_the_writer_writes():
    """INVERTED: this asserted a write to the view reached the `fields` list underneath."""
    mol = parse_record(_MINIMAL + ['>  <ID>', '7', ''])
    mol.meta['ID'] = 'x'
    assert parse_record(emit_record(mol)[0]).meta == {'ID': 'x'}


def test_emit_record_writes_the_molecule_s_own_meta():
    mol = parse_record(_MINIMAL + ['>  <NAME>', 'methane', '', '$$$$'])
    assert '>  <NAME>' in emit_record(mol)[0]


def test_emit_record_writes_none_when_told_none():
    mol = parse_record(_MINIMAL + ['>  <NAME>', 'methane', '', '$$$$'])
    assert '>  <NAME>' not in emit_record(mol, meta={})[0]


def test_the_unparsed_bucket_is_not_written_back():
    """It holds lines that were not a field, so re-emitting it would invent one."""
    mol = parse_record(_MINIMAL + ['stray', '$$$$'])
    lines, log = emit_record(mol)
    assert f'>  <{UNPARSED_KEY}>' not in lines and any('not written back' in str(x) for x in log)


def test_mol_round_trips_its_data_fields():
    """The done-when, on one line."""
    text = '\n'.join(_MINIMAL + ['>  <NAME>', 'methane', '', '$$$$'])
    assert parse_record(emit_record(parse_record(text.split('\n')))[0]).meta == {'NAME': 'methane'}


def test_a_value_with_a_newline_is_written_as_two_lines_and_read_back_as_one_value():
    mol = parse_record(_MINIMAL)
    mol.meta['NOTE'] = 'line one\nline two'
    lines, _ = emit_record(mol)
    assert lines[lines.index('>  <NOTE>') + 1:lines.index('>  <NOTE>') + 3] == ['line one',
                                                                                'line two']
    assert parse_record(lines).meta == {'NOTE': 'line one\nline two'}


def test_a_record_round_trips_through_emit_and_parse():
    header = {}
    mol = parse_record(_MINIMAL + ['>  <ID>', '7', ''], [], header=header)
    lines, _ = emit_record(mol, header['sgroups'], title='methane')
    assert lines[-1] == RECORD_SEPARATOR
    mol2 = parse_record(lines, [])
    assert [mol2.element_of(s) for s in mol2.atom_numbers] == [6]
    assert mol2.implicit_h_of(next(iter(mol2.atom_numbers))) == 4
    assert mol2.meta == {'ID': '7'}


def test_a_field_value_is_closed_by_a_blank_line_on_write():
    """Without it the next ``>`` line is read as more of the previous value."""
    mol = parse_record(_MINIMAL)
    lines, _ = emit_record(mol, meta={'A': '1', 'B': '2'})
    assert parse_record(lines, []).meta == {'A': '1', 'B': '2'}


def test_a_sequence_of_pairs_is_accepted_where_a_mapping_is_expected():
    """A caller that built its fields as an ordered list of pairs need not go through a dict first."""
    mol = parse_record(_MINIMAL)
    lines, _ = emit_record(mol, meta=[('A', 1)])
    assert parse_record(lines, []).meta == {'A': '1'}


def test_an_unknown_version_is_refused_by_name():
    mol = parse_record(_MINIMAL, [])
    with raises(MalformedCtfile, match='V2000'):
        emit_record(mol, version='V4000')
