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
"""XYZ reader tests.

A trailing ``[mutant: ...]`` line names the implementation line whose deletion or
alteration makes the test fail.
"""

from pathlib import Path

import pytest

from ...core import MoleculeContainer, smiles
from ..xyz import (XYZAtom, XYZFrame, _looks_like_atom_line, _looks_like_truncated_atom_line,
                   _normalize_element, _parse_count, xyz, xyz_conformers)


_DATA = Path(__file__).resolve().parent.parent.parent.parent / 'test'


# fixtures


@pytest.fixture
def ch_xyz():
    return (_DATA / 'ch.xyz').read_text(encoding='utf-8')


@pytest.fixture
def truncated_xyz():
    return (_DATA / 'xyz_truncated.xyz').read_text(encoding='utf-8')


@pytest.fixture
def extended_xyz():
    return (_DATA / 'xyz_extended.xyz').read_text(encoding='utf-8')


@pytest.fixture
def two_frame_xyz():
    return (_DATA / 'xyz_two_frame.xyz').read_text(encoding='utf-8')


# well-formed parsing


def test_ch_xyz_parses_all_frames(ch_xyz):
    """ch.xyz contains eight XYZ blocks concatenated; xyz() must return all eight.

    [mutant: replace `frames.append(frame)` with `return [frame]` after the first]
    """
    frames = xyz(ch_xyz)
    assert len(frames) == 8, f'expected 8 frames, got {len(frames)}'


def test_ch_xyz_first_frame_atom_count(ch_xyz):
    """The first frame of ch.xyz has 30 atoms, as its count line states.

    [mutant: remove the inner atom-parsing loop so atoms is always empty]
    """
    frames = xyz(ch_xyz)
    assert frames[0].stated_count == 30
    assert len(frames[0]) == 30


def test_ch_xyz_last_frame_atom_count(ch_xyz):
    """The last frame has 27 atoms and a 'charge=1' annotation in the comment line.

    [mutant: mistake the comment line for an atom line]
    """
    frames = xyz(ch_xyz)
    last = frames[-1]
    assert last.stated_count == 27
    assert len(last) == 27
    assert 'charge=1' in last.title


def test_all_frames_have_correct_atom_counts(ch_xyz):
    """Every frame in ch.xyz has stated_count == len(atoms) and an empty log: a non-empty log on
    a well-formed file means the reader invented a problem.

    [mutant: always log a spurious 'record: count' message]
    """
    frames = xyz(ch_xyz)
    for i, frame in enumerate(frames):
        assert frame.stated_count == len(frame), (
            f'frame {i}: stated {frame.stated_count}, got {len(frame)}')
        assert frame.log == [], f'frame {i}: unexpected log entries: {frame.log}'


def test_element_symbols_are_valid(ch_xyz):
    """Every atom element in ch.xyz is a recognized element symbol; the corpus holds only
    C, N, O, H, Na, all in proper case.

    [mutant: skip step 1 (direct match) in _normalize_element]
    """
    from chython.core._core import element_symbols
    valid = frozenset(element_symbols()[1:])
    frames = xyz(ch_xyz)
    for i, frame in enumerate(frames):
        for j, atom in enumerate(frame.atoms):
            assert atom.element in valid, (
                f'frame {i} atom {j}: {atom.element!r} is not a valid element symbol')


def test_coordinates_are_stored(ch_xyz):
    """x, y and z are all stored in XYZAtom: the container cannot hold z but the intermediate
    must.

    [mutant: store only x and y, set z to 0.0 always]
    """
    frames = xyz(ch_xyz)
    # ch.xyz has non-zero z throughout; a lost z reads as 0.0
    any_nonzero_z = any(abs(a.z) > 0.001 for f in frames for a in f.atoms)
    assert any_nonzero_z, 'all z values are zero -- z coordinates were not stored'


def test_comment_line_preserved(ch_xyz):
    """The comment line is stored verbatim in XYZFrame.title, where the second callable looks for
    charge/radical annotations.

    [mutant: always set frame.title to empty string]
    """
    frames = xyz(ch_xyz)
    assert frames[-1].title.strip() == 'charge=1'


# multi-frame


def test_two_frame_file_returns_two_frames(two_frame_xyz):
    """A two-frame XYZ file produces exactly two XYZFrame objects.

    [mutant: break after the first frame is appended]
    """
    frames = xyz(two_frame_xyz)
    assert len(frames) == 2


def test_two_frame_counts(two_frame_xyz):
    """First frame has 2 atoms; second frame has 3 atoms.

    [mutant: off-by-one in atom-count loop]
    """
    frames = xyz(two_frame_xyz)
    assert len(frames[0]) == 2
    assert len(frames[1]) == 3


def test_two_frame_second_frame_elements(two_frame_xyz):
    """Second frame contains O, H, H -- read from the second frame, not the first.

    [mutant: re-read atoms from the start of the string for each frame]
    """
    frames = xyz(two_frame_xyz)
    elems = [a.element for a in frames[1].atoms]
    assert elems == ['O', 'H', 'H']


# truncated file


def test_truncated_file_returns_frame(truncated_xyz):
    """A file whose atom block is shorter than its count line states still yields a frame: we keep
    what we have rather than discarding it.

    [mutant: return [] when atoms_found != count]
    """
    frames = xyz(truncated_xyz)
    assert len(frames) == 1


def test_truncated_file_logs_record_prefix(truncated_xyz):
    """The count discrepancy is logged with the 'record:' prefix, not silently: a reader that kept
    2 atoms and logged nothing passes a count check while hiding the missing data.

    [mutant: remove the count-mismatch log append]
    """
    frames = xyz(truncated_xyz)
    log_entries = frames[0].log
    assert any(str(e).startswith('record:') for e in log_entries), (
        f'expected a record: entry; got: {log_entries}')


def test_truncated_file_atom_count_matches_reality(truncated_xyz):
    """xyz_truncated.xyz states 5 atoms but has 2; we keep the 2 that are there.

    [mutant: fill missing atoms with dummy entries instead of stopping at EOF]
    """
    frames = xyz(truncated_xyz)
    assert frames[0].stated_count == 5
    assert len(frames[0]) == 2


def test_truncated_file_caller_log_gets_record_entry(truncated_xyz):
    """The caller-supplied log also receives the record-prefix entry; the global log and frame.log
    must agree, or aggregated log processing misses the issue.

    [mutant: append only to frame.log and not to the caller's log]
    """
    global_log: list = []
    frames = xyz(truncated_xyz, log=global_log)
    assert any(str(e).startswith('record:') for e in global_log), global_log


# extended XYZ


def test_extended_xyz_returns_correct_atom_count(extended_xyz):
    """An extended-XYZ comment line does not prevent parsing the atom block.

    [mutant: return [] when the comment looks like extended XYZ]
    """
    frames = xyz(extended_xyz)
    assert len(frames) == 1
    assert len(frames[0]) == 3


def test_extended_xyz_logs_unsupported_for_lattice(extended_xyz):
    """The Lattice field is logged as 'unsupported:' because we do not model it -- that prefix is
    how a caller asks 'does this record use features we do not model?'.

    [mutant: use 'atom:' prefix instead of 'unsupported:' for Lattice]
    """
    frames = xyz(extended_xyz)
    unsupported = [e for e in frames[0].log if str(e).startswith('unsupported:')]
    keys = {e for e in unsupported if 'Lattice' in e}
    assert keys, f'expected unsupported: entry for Lattice; log: {frames[0].log}'


def test_extended_xyz_logs_unsupported_for_properties(extended_xyz):
    """Properties field is logged as unsupported: because we do not model it.

    [mutant: omit Properties from the unsupported logging loop]
    """
    frames = xyz(extended_xyz)
    unsupported = [e for e in frames[0].log if str(e).startswith('unsupported:')]
    keys = {e for e in unsupported if 'Properties' in e}
    assert keys, f'expected unsupported: entry for Properties; log: {frames[0].log}'


def test_extended_xyz_logs_unsupported_for_energy(extended_xyz):
    """The energy field is also an unmodelled extended-XYZ field.

    [mutant: only log Lattice and Properties, skip other keys]
    """
    frames = xyz(extended_xyz)
    unsupported = [e for e in frames[0].log if str(e).startswith('unsupported:')]
    keys = {e for e in unsupported if 'energy' in e}
    assert keys, f'expected unsupported: entry for energy; log: {frames[0].log}'


def test_extended_xyz_atom_elements_correct(extended_xyz):
    """Atom elements are read correctly even when the comment carries extended metadata.

    [mutant: confuse the comment line's key=value tokens with atom tokens]
    """
    frames = xyz(extended_xyz)
    elems = [a.element for a in frames[0].atoms]
    assert elems == ['C', 'C', 'N']


# element normalization


def test_wrong_case_symbol_is_corrected():
    """'CL' is corrected to 'Cl' and logged with the 'atom:' prefix; storing 'CL' raw would fail
    later when the second callable calls add_atom('CL').

    [mutant: skip the case-normalization step in _normalize_element]
    """
    text = '1\n\nCL 0.0 0.0 0.0\n'
    log: list = []
    frames = xyz(text, log=log)
    assert len(frames) == 1
    assert frames[0].atoms[0].element == 'Cl'
    assert any(str(e).startswith('atom:') and 'CL' in e for e in log), log


def test_wrong_case_fe_corrected():
    """'FE' corrected to 'Fe' -- iron, not a typo.

    [mutant: only normalize two-char symbols starting with a vowel]
    """
    text = '1\n\nFE 0.0 0.0 0.0\n'
    frames = xyz(text)
    assert frames[0].atoms[0].element == 'Fe'
    assert any('FE' in e for e in frames[0].log)


def test_atomic_number_symbol_normalized():
    """A bare integer is treated as an atomic number; some codes write '6' for carbon.

    [mutant: skip the isdigit branch in _normalize_element]
    """
    text = '1\n\n6 0.0 0.0 0.0\n'
    log: list = []
    frames = xyz(text, log=log)
    assert frames[0].atoms[0].element == 'C'
    assert any(str(e).startswith('atom:') and 'atomic number' in e for e in log), log


def test_trailing_digit_stripped():
    """'C1' has its trailing digit stripped to 'C'; Tinker-style XYZ writes 'C1', 'N3'.

    [mutant: skip the trailing-character-stripping step]
    """
    text = '1\n\nC1 0.0 0.0 0.0\n'
    log: list = []
    frames = xyz(text, log=log)
    assert frames[0].atoms[0].element == 'C'
    assert any(str(e).startswith('atom:') and 'C1' in e for e in log), log


def test_unknown_element_stored_and_logged():
    """An unknown symbol like 'X' or 'Du' is stored as-is and logged with 'atom:', so the second
    callable sees the original token; dropping the atom would shift every later atom's index.

    [mutant: replace unknown symbols with 'C' silently]
    """
    text = '2\n\nX 0.0 0.0 0.0\nDu 1.0 0.0 0.0\n'
    log: list = []
    frames = xyz(text, log=log)
    assert len(frames[0].atoms) == 2
    assert frames[0].atoms[0].element == 'X'
    assert frames[0].atoms[1].element == 'Du'
    assert any('atom:' in e and 'X' in e for e in log), log
    assert any('atom:' in e and 'Du' in e for e in log), log


def test_deuterium_stored_as_hydrogen_with_isotope():
    """'D' is deuterium, which many quantum-chemistry codes write in XYZ output: stored as 'H'
    with isotope 2, since the container cannot hold it as 'D'.

    [mutant: skip the 'D' special case and store 'D' raw]
    """
    text = '1\n\nD 0.0 0.0 0.0\n'
    log: list = []
    frames = xyz(text, log=log)
    atom = frames[0].atoms[0]
    assert atom.element == 'H'
    assert atom.isotope == 2
    assert any('atom:' in e and 'D' in e for e in log), log


def test_tritium_stored_as_hydrogen_with_isotope():
    """'T' is tritium: stored as 'H' with isotope 3.

    [mutant: handle D but not T]
    """
    text = '1\n\nT 0.0 0.0 0.0\n'
    frames = xyz(text)
    atom = frames[0].atoms[0]
    assert atom.element == 'H'
    assert atom.isotope == 3


def test_unknown_element_atom_count_not_affected():
    """An unrecognized symbol still counts as an atom: dropping it would renumber the atoms after
    it and break any index-based reference.

    [mutant: skip appending the atom when normalization returns an unknown symbol]
    """
    text = '3\n\nC 0.0 0.0 0.0\nX 1.0 0.0 0.0\nN 2.0 0.0 0.0\n'
    frames = xyz(text)
    assert len(frames[0].atoms) == 3


# line-ending tolerance


def test_crlf_line_endings(two_frame_xyz):
    """CRLF line endings in an XYZ file parse identically to LF.

    [mutant: split on '\\n' only instead of using splitlines()]
    """
    crlf = two_frame_xyz.replace('\n', '\r\n')
    frames = xyz(crlf)
    assert len(frames) == 2
    assert len(frames[0]) == 2
    assert len(frames[1]) == 3


def test_missing_trailing_newline():
    """An XYZ file with no final newline still parses its last atom.

    [mutant: require splitlines() to return a trailing empty element]
    """
    text = '2\n\nC 0.0 0.0 0.0\nN 1.0 0.0 0.0'  # no trailing newline
    frames = xyz(text)
    assert len(frames) == 1
    assert len(frames[0]) == 2


def test_blank_lines_between_frames():
    """Blank lines between XYZ frames cost no frames; some trajectory writers insert one as a
    separator.

    [mutant: treat a blank line as a count-line parse failure and abort]
    """
    text = '1\nframe 1\nC 0.0 0.0 0.0\n\n\n2\nframe 2\nC 0.0 0.0 0.0\nN 1.0 0.0 0.0\n'
    frames = xyz(text)
    assert len(frames) == 2
    assert len(frames[0]) == 1
    assert len(frames[1]) == 2


def test_empty_string_returns_no_frames():
    """Empty input yields an empty list, not an exception.

    [mutant: raise ValueError on empty input]
    """
    assert xyz('') == []


def test_whitespace_only_string_returns_no_frames():
    """All-whitespace input contains no frames.

    [mutant: try to parse '   ' as a count line]
    """
    assert xyz('   \n\n\t\n') == []


# coordinate storage


def test_z_coordinates_are_nonzero_when_present():
    """z from a 3D file is stored in XYZAtom.z, not silently zeroed: the container keeps only x
    and y, but the intermediate must hold all three.

    [mutant: always set z = 0.0 in XYZAtom.__init__]
    """
    text = '1\n\nC 1.0 2.0 3.0\n'
    frames = xyz(text)
    assert frames[0].atoms[0].z == pytest.approx(3.0)


def test_xy_coordinates_stored():
    """x and y are round-tripped without modification.

    [mutant: swap x and y when storing]
    """
    text = '1\n\nC 1.23 -4.56 7.89\n'
    frames = xyz(text)
    atom = frames[0].atoms[0]
    assert atom.x == pytest.approx(1.23)
    assert atom.y == pytest.approx(-4.56)
    assert atom.z == pytest.approx(7.89)


# log aggregation


def test_caller_log_receives_all_frame_logs(truncated_xyz, two_frame_xyz):
    """A caller-supplied log receives entries from every frame, not just the last: without
    aggregation a caller processing a trajectory sees only the last frame's issues.

    [mutant: clear log between frames or only append at the end]
    """
    combined = truncated_xyz + two_frame_xyz
    global_log: list = []
    frames = xyz(combined, log=global_log)
    assert len(frames) == 3  # 1 truncated + 2 from two_frame
    assert any('record:' in e for e in global_log), global_log

    for frame in frames:
        for entry in frame.log:
            assert entry in global_log, (
                f'frame entry {entry!r} missing from global log')


# _normalize_element unit tests


def test_normalize_direct_match():
    """A symbol already in proper case needs no normalization.

    [mutant: always apply case correction even when the symbol matches]
    """
    sym, iso, msg = _normalize_element('C')
    assert sym == 'C'
    assert iso == 0
    assert msg is None


def test_normalize_two_char_direct():
    """Two-char symbol in proper case: no correction needed.

    [mutant: capitalize all input before the direct-match check]
    """
    sym, iso, msg = _normalize_element('Fe')
    assert sym == 'Fe'
    assert iso == 0
    assert msg is None


def test_normalize_all_upper():
    """All-uppercase symbols are corrected, with an 'atom:' log entry.

    [mutant: fall through to 'unrecognized' when case normalization is skipped]
    """
    sym, iso, msg = _normalize_element('CL')
    assert sym == 'Cl'
    assert iso == 0
    assert msg is not None and str(msg).startswith('atom:')


def test_normalize_atomic_number():
    """Bare integer maps to the element at that atomic number.

    [mutant: return raw string for isdigit inputs]
    """
    sym, iso, msg = _normalize_element('8')
    assert sym == 'O'
    assert iso == 0
    assert msg is not None and 'atomic number' in msg


def test_normalize_out_of_range_atomic_number():
    """An out-of-range atomic number is stored raw and logged.

    [mutant: wrap out-of-range numbers modulo 118]
    """
    sym, iso, msg = _normalize_element('200')
    assert sym == '200'   # stored raw
    assert msg is not None and 'out of range' in msg


def test_normalize_deuterium():
    """'D' → 'H', isotope 2.

    [mutant: treat D as unknown]
    """
    sym, iso, msg = _normalize_element('D')
    assert sym == 'H'
    assert iso == 2
    assert msg is not None


def test_normalize_tritium():
    """'T' → 'H', isotope 3.

    [mutant: treat T as unknown]
    """
    sym, iso, msg = _normalize_element('T')
    assert sym == 'H'
    assert iso == 3
    assert msg is not None


def test_normalize_trailing_digit():
    """'N3' has trailing digit stripped to 'N'.

    [mutant: accept 'N3' as-is without stripping]
    """
    sym, iso, msg = _normalize_element('N3')
    assert sym == 'N'
    assert msg is not None and 'N3' in msg


def test_normalize_unknown():
    """Completely unknown symbol is stored raw.

    [mutant: replace unknown symbols with 'C']
    """
    sym, iso, msg = _normalize_element('Xx')
    assert sym == 'Xx'
    assert msg is not None and 'unrecognized' in msg


# understated count


def test_understated_count_is_logged():
    """A count smaller than the actual atom-line count logs a 'record:' entry.  The reader keeps
    only the stated N atoms, so multi-frame stays synchronised, but the surplus must be reported
    or an atom vanishes with no indication.

    [mutant: remove the 'elif i < n' surplus-peek block so surplus lines are silent]
    """
    # count=2, but three atom lines follow before the next frame or EOF
    text = '2\nunderstated header\nC 0.0 0.0 0.0\nN 1.0 0.0 0.0\nO 2.0 0.0 0.0\n'
    log: list = []
    frames = xyz(text, log=log)
    assert len(frames) == 1
    assert len(frames[0]) == 2
    record_entries = [e for e in frames[0].log if str(e).startswith('record:') and 'surplus' in e]
    assert record_entries, f'expected a record: surplus entry; got: {frames[0].log}'
    assert '1' in record_entries[0], f'expected surplus count in message: {record_entries[0]!r}'


def test_second_frame_resyncs_after_understated_count():
    """When the first frame's count is understated, frame N+1 is still found and parsed: the outer
    loop scans for a line that parses as a positive integer count, and surplus atom lines do not.

    This does not discriminate whether the surplus-peek advances past the surplus lines; that is
    pinned by ``test_a_surplus_line_is_reported_once``.

    [mutant: stop the outer loop on the first non-integer line after a frame]
    """
    text = (
        '2\nfirst frame (understated)\n'
        'C 0.0 0.0 0.0\nN 1.0 0.0 0.0\n'
        'O 2.0 0.0 0.0\n'          # surplus -- NOT in the first frame
        '3\nsecond frame\n'
        'C 0.0 0.0 0.0\nH 1.0 0.0 0.0\nH -1.0 0.0 0.0\n'
    )
    frames = xyz(text)
    assert len(frames) == 2, f'expected 2 frames; got {len(frames)}'
    assert [a.element for a in frames[1].atoms] == ['C', 'H', 'H']
    assert frames[1].log == [], f'second frame has unexpected log: {frames[1].log}'


# a count line that is not one


@pytest.mark.parametrize('header', ['3 atoms', '3.0', ' 3 # water', 'natoms=3', '-3'])
def test_a_header_that_states_no_count_is_reported_not_dropped(header):
    """Three well-formed atom lines behind an unreadable header are still three atoms lost.  No
    frame is produced -- a reconstructed count is an invented frame -- but the block is named with
    how many atom lines went unread.

    [mutant: `continue` on a bad count line without appending a log entry]
    """
    log: list = []
    frames = xyz(header + '\nc\nC 0 0 0\nN 1 0 0\nO 2 0 0\n', log=log)
    assert frames == []
    record_entries = [e for e in log if str(e).startswith('record:')]
    assert record_entries, f'a whole block vanished with an empty log: {log}'
    assert '3 atom line(s)' in record_entries[0], record_entries[0]


def test_a_byte_order_mark_costs_its_frame_and_says_so():
    """A Windows-written trajectory glues a UTF-8 BOM to the first count line.  That frame cannot
    be read, and returning the second frame as if it were the whole file is the silent
    first-frame-only mode this reader exists to prevent.

    [mutant: strip the BOM silently, or skip the block with no log entry]
    """
    log: list = []
    frames = xyz('﻿2\nfirst\nC 0 0 0\nN 1 0 0\n2\nsecond\nO 0 0 0\nS 1 0 0\n', log=log)
    assert [f.title for f in frames] == ['second']
    assert any(str(e).startswith('record:') and '2 atom line(s)' in e for e in log), log


def test_a_stray_block_does_not_swallow_the_frame_that_follows_it():
    """The skip stops at the next count line rather than consuming to end of file, or one bad header
    would cost every frame after it too.

    [mutant: advance to the end of the input instead of breaking on the next count line]
    """
    log: list = []
    frames = xyz('two\nc\nC 0 0 0\n2\nreal\nO 0 0 0\nS 1 0 0\n', log=log)
    assert [f.title for f in frames] == ['real']
    assert [a.element for a in frames[0].atoms] == ['O', 'S']


# a truncated surplus atom line


def test_a_truncated_atom_line_past_the_count_is_reported_like_one_inside_it():
    """A truncated atom line one past the stated count -- the shape a killed job leaves -- is
    reported the same way as one inside the count; the damage does not depend on whether it fits.

    [mutant: count only complete atom lines in the surplus peek]
    """
    log: list = []
    frames = xyz('1\nc\nC 0.0 0.0 0.0\nN 1.0 0.0\n', log=log)
    assert len(frames) == 1
    assert len(frames[0]) == 1            # the surplus atom is not rescued
    assert any(str(e).startswith('record:') and 'surplus' in e for e in frames[0].log), frames[0].log
    assert any(str(e).startswith('atom:') and 'N 1.0 0.0' in e for e in frames[0].log), frames[0].log
    assert frames[0].log == log           # every frame entry reaches the caller


def test_a_surplus_line_is_reported_once():
    """Reported lines are consumed.  Left in the stream they meet the outer loop, which reports
    them again as 'this is not a count line' -- one thing, two names.

    [mutant: peek without advancing, so the outer loop meets the surplus lines again]
    """
    log: list = []
    xyz('1\nc\nC 0.0 0.0 0.0\nN 1.0 0.0 0.0\nO 2.0 0.0 0.0\n', log=log)
    assert len(log) == 1, log


# the line discriminants


def test_extra_columns_after_the_coordinates_are_an_atom_line():
    """Trailing columns are ignored by the atom parser, so the discriminant ignores them too: it
    is tokens 1..3 that decide, not the last three tokens.

    [mutant: check the last three tokens instead of tokens 1..3]
    """
    assert _looks_like_atom_line('C 0.0 0.0 0.0 junk')
    assert not _looks_like_atom_line('C 0.0 junk 0.0')


def test_a_count_line_is_never_mistaken_for_a_truncated_atom_line():
    """Both surplus shapes stay distinguishable from the next frame's header, which is what makes
    the one-line lookahead unambiguous rather than heuristic.

    [mutant: accept a single-token line as a truncated atom line]
    """
    assert not _looks_like_truncated_atom_line('12')
    assert _looks_like_truncated_atom_line('N 1.0 0.0')
    assert not _looks_like_truncated_atom_line('N 1.0 x')


def test_one_definition_of_a_count_line():
    """A shape rejected as a frame header must not be accepted as the end of a frame, or a frame
    ends on a line the reader then refuses to start a frame with.

    [mutant: a second int() test somewhere in the reader]
    """
    assert _parse_count('3') == 3
    assert _parse_count('0') == 0
    assert _parse_count('-3') is None
    assert _parse_count('3.0') is None
    assert _parse_count('3 atoms') is None


# the documented contract


def test_a_malformed_atom_line_is_the_third_reason_the_counts_differ():
    """`stated_count` is what the file said and is never revised to match what was found; an atom
    line that could not be read is an `atom:` entry, not a `record:` one.

    [mutant: overwrite stated_count with len(atoms), or report a record: line for a malformed atom]
    """
    frame = xyz('3\nc\nC 0 0 0\nBOGUS\nN 1 0 0\n')[0]
    assert frame.stated_count == 3
    assert len(frame) == 2
    assert [str(e).split(':')[0] for e in frame.log] == ['atom']


def test_an_atom_is_hashable_because_it_compares_by_value():
    """Both classes are public through the facade, and a value `__eq__` without `__hash__` makes
    `set(frame.atoms)` raise on a class whose whole purpose is to be collected.

    [mutant: delete __hash__, or hash by identity while __eq__ compares by value]
    """
    one = XYZAtom('C', 0, 1.0, 2.0, 3.0)
    same = XYZAtom('C', 0, 1.0, 2.0, 3.0)
    other = XYZAtom('N', 0, 1.0, 2.0, 3.0)
    assert one == same and hash(one) == hash(same)
    assert len({one, same, other}) == 2
    assert len(set(xyz('2\nc\nC 0 0 0\nC 0 0 0\n')[0].atoms)) == 1
    assert isinstance(hash(XYZFrame()), int)      # the sibling stays hashable by identity


@pytest.mark.parametrize('bad', [None, b'1\n\nC 0 0 0\n', 42, ['1', '', 'C 0 0 0']])
def test_a_non_string_argument_is_refused_at_the_call_site(bad):
    """A non-string argument is refused at the call site, not several lines into the parse from
    inside the element normalizer.

    [mutant: drop the isinstance check and let splitlines() decide]
    """
    with pytest.raises(TypeError):
        xyz(bad)


# the frame's own log, and why it is the frame's


def test_a_frames_damage_reaches_the_frame_with_nothing_passed_in():
    """`XYZFrame.log` is the destination, and the caller's `log=` list is the second copy.

    The frame is not a `MoleculeContainer` -- XYZ states no bond, so there is nothing to put a
    container's log on until `chython.chemistry.saturate()` is run on a molecule some other pass built
    -- which is why the record object keeps a `log` of its own where a reader that returns a container
    would write to `mol.log`.

    [mutant: append only to the caller's list in the atom-line branches]
    """
    frames = xyz('2\ntwo damaged lines\nCL 0.0 0.0 0.0\nX 1.0 0.0 0.0\n')   # no log= anywhere
    assert len(frames) == 1
    assert not isinstance(frames[0], MoleculeContainer)
    assert [str(x.rule) for x in frames[0].log] == ['xyz:symbol-case-corrected',
                                                    'xyz:symbol-unrecognized'], frames[0].log


# frames onto a molecule the caller brings


#: Two frames of one water, the second nudged along x.  The topology is the caller's; XYZ states none.
TRAJECTORY = """3
step 0
O  0.000  0.000  0.000
H  0.757  0.586  0.000
H -0.757  0.586  0.000
3
step 1
O  0.010  0.000  0.000
H  0.767  0.586  0.000
H -0.747  0.586  0.000
"""


def _water():
    """Water with both hydrogens explicit, in the trajectory's O, H, H order."""
    return smiles('O([H])[H]')


def test_frames_become_conformers_of_a_molecule_the_caller_brings():
    """Two frames, two models, each carrying its ordinal in the file.

    [mutant: store the first frame only]
    """
    molecule = _water()
    assert xyz_conformers(molecule, xyz(TRAJECTORY)) == 2
    assert len(molecule.conformers) == 2
    assert [c.ext_index for c in molecule.conformers] == [0, 1]
    assert molecule.conformer(1).coordinates[0] == (0.01, 0.0, 0.0)


def test_frames_append_to_the_models_a_molecule_already_carries():
    """A second call keeps the first call's models and adds after them.

    [mutant: drop the existing models before storing]
    """
    molecule = _water()
    assert xyz_conformers(molecule, xyz(TRAJECTORY)) == 2
    assert xyz_conformers(molecule, xyz(TRAJECTORY)) == 2
    assert len(molecule.conformers) == 4
    assert [c.ext_index for c in molecule.conformers] == [0, 1, 0, 1]


def test_a_frame_of_the_wrong_length_is_logged_and_skipped():
    """All-or-nothing per frame: the short frame stores nothing and the other two still land."""
    molecule = _water()
    text = TRAJECTORY + '2\nshort\nO 0. 0. 0.\nH 1. 0. 0.\n'
    log: list = []
    assert xyz_conformers(molecule, xyz(text), log=log) == 2
    assert len(molecule.conformers) == 2
    assert any('atom(s) where the molecule holds' in x for x in log), log
    assert any('atom(s) where the molecule holds' in x for x in molecule.log), molecule.log


def test_a_transposed_frame_is_refused_on_the_element_sequence():
    """The positional match is the whole contract, so a frame that only agrees on the count is not a
    frame of this molecule.

    [mutant: check the count and not the element sequence]
    """
    molecule = _water()
    swapped = '3\nswapped\nH  0.757  0.586  0.000\nO  0.  0.  0.\nH -0.757  0.586  0.000\n'
    log: list = []
    assert xyz_conformers(molecule, xyz(swapped), log=log) == 0
    assert not molecule.has_3d
    assert any('position 0' in x for x in log), log
