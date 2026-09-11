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
"""The file surface: what a caller gets from a path, a buffer, a bad record, an early exit -- plus
that ``from chython import SDFRead`` resolves to this reader, which no lower test can state.
"""

from io import StringIO
from pathlib import Path

from pytest import raises

from .._facade import mol
from .._sdf import RECORD_SEPARATOR, parse_record
from .._stream import ESDFWrite, FailedRecord, SDFRead, SDFWrite
from .._v3000 import V3000_STAMP
from ....core import read_smiles


_METHANE = ['methane', '  test', '', '  1  0  0  0  0  0            999 V2000',
            '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0', 'M  END']
_BROKEN = ['broken', '', '', '  not a counts line', 'M  END']
#: Two carbons joined by bond type 4.  Nonsense as chemistry; what is under test is that the bond
#: type reaches the file, which needs no ring.
_AROMATIC = ['aromatic', '', '', '  2  1  0  0  0  0            999 V2000',
             '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
             '    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
             '  1  2  4  0  0  0  0', 'M  END']


def _sdf(*records):
    """One SDF text from record line-lists, each closed by the separator."""
    out = []
    for record in records:
        out.extend(record)
        out.append(RECORD_SEPARATOR)
    return '\n'.join(out) + '\n'


# registration

def test_the_package_root_names_resolve_to_this_reader():
    """``chython.SDFRead`` is what a user reaches for, and until it is this class the reader cannot be
    used.  ``formats.__all__`` is asserted too, so the export list and the resolved object cannot
    drift apart."""
    import chython
    from chython import formats

    assert chython.SDFRead is SDFRead
    assert chython.SDFWrite is SDFWrite
    assert chython.ESDFWrite is ESDFWrite
    assert 'SDFRead' in formats.__all__


# what a file is

def test_a_buffer_is_read():
    with SDFRead(StringIO(_sdf(_METHANE, _METHANE))) as f:
        assert len(f.read()) == 2


def test_a_str_path_and_a_path_object_read_the_same(tmp_path):
    path = tmp_path / 'x.sdf'
    path.write_text(_sdf(_METHANE, _METHANE))
    with SDFRead(str(path)) as f:
        by_str = len(f.read())
    with SDFRead(Path(path)) as f:
        assert len(f.read()) == by_str == 2


def test_anything_else_is_refused_by_type():
    """Refused at construction, not as a confusing parse error two hundred records later."""
    with raises(TypeError, match='invalid file'):
        SDFRead(42)


def test_the_callers_buffer_is_not_closed_for_them():
    """Closing someone else's handle is a bug that surfaces in *their* next read, far from here."""
    buffer = StringIO(_sdf(_METHANE))
    with SDFRead(buffer) as f:
        f.read()
    assert not buffer.closed
    f.close(force=True)
    assert buffer.closed


def test_a_file_the_reader_opened_is_closed_on_exit(tmp_path):
    path = tmp_path / 'x.sdf'
    path.write_text(_sdf(_METHANE))
    with SDFRead(str(path)) as f:
        f.read()
    assert f._file.closed


# reading

def test_iteration_yields_one_molecule_per_record():
    with SDFRead(StringIO(_sdf(_METHANE, _METHANE, _METHANE))) as f:
        assert [m.atom_count for m in f] == [1, 1, 1]


def test_read_takes_a_count_and_leaves_the_rest():
    with SDFRead(StringIO(_sdf(_METHANE, _METHANE, _METHANE))) as f:
        assert len(f.read(2)) == 2
        assert len(f.read()) == 1


def test_tell_counts_records_not_lines():
    with SDFRead(StringIO(_sdf(_METHANE, _METHANE))) as f:
        assert f.tell() == -1
        f.read_structure()
        assert f.tell() == 0
        f.read_structure()
        assert f.tell() == 1


def test_the_side_data_follows_the_current_record():
    """The reader holds `meta`/`title`, so they must move with the record: stale side data from the
    previous one is worse than none."""
    a = _METHANE + ['>  <ID>', 'first', '']
    b = list(_METHANE)
    b[0] = 'ethane-ish'
    b += ['>  <ID>', 'second', '']
    with SDFRead(StringIO(_sdf(a, b))) as f:
        f.read_structure()
        assert f.meta == {'ID': 'first'}
        assert f.title == 'methane'
        f.read_structure()
        assert f.meta == {'ID': 'second'}
        assert f.title == 'ethane-ish'


def test_before_the_first_record_the_side_data_is_empty_rather_than_an_error():
    with SDFRead(StringIO(_sdf(_METHANE))) as f:
        assert f.record is None
        assert f.meta == {} and f.title == '' and f.log == [] and f.unknown_hydrogens == ()


def test_read_record_hands_over_the_container_and_the_reader_keeps_the_framing():
    """INVERTED: `read_record` returned a wrapper carrying both.  The container holds the meta and the
    title itself now, and only what no container holds stays on the reader."""
    with SDFRead(StringIO(_sdf(_METHANE + ['>  <ID>', '7', '']))) as f:
        mol = f.read_record()
        assert f.version == 'V2000'
    assert mol.atom_count == 1
    assert mol.meta == {'ID': '7'}
    assert mol.title == 'methane'


# a bad record

def test_a_bad_record_is_collected_and_the_good_ones_keep_coming():
    """Why this class exists rather than a loop over ``parse_record``: skipping silently loses a
    tenth of a pipeline's input with nothing to show, and raising loses the other nine tenths."""
    with SDFRead(StringIO(_sdf(_METHANE, _BROKEN, _METHANE))) as f:
        mols = f.read()
    assert len(mols) == 2
    assert len(f.failed) == 1
    failed, = f.failed
    assert isinstance(failed, FailedRecord)
    assert failed.position == 1
    assert 'not a counts line' in failed.text
    assert 'counts line' in str(failed.error)


def test_a_failed_record_keeps_its_place_in_the_file():
    """Its position is the record index, so a report can say *which* record to go and look at."""
    with SDFRead(StringIO(_sdf(_METHANE, _BROKEN, _METHANE, _BROKEN))) as f:
        f.read()
    assert [x.position for x in f.failed] == [1, 3]


def test_a_record_that_only_needed_recovering_is_not_a_failure():
    """A repaired-and-logged malformation must not land in `failed`, or "damaged" and "unreadable"
    stop being different."""
    record = list(_METHANE)
    record[3] = '  1  0  0  0  0  0            999'
    with SDFRead(StringIO(_sdf(record))) as f:
        f.read_structure()
        assert not f.failed
        assert any('no version stamp' in x for x in f.log), f.log


# writing

def test_write_then_read_returns_the_molecule_and_its_meta(tmp_path):
    path = tmp_path / 'out.sdf'
    with SDFRead(StringIO(_sdf(_METHANE))) as f:
        mol = f.read_record()
    with SDFWrite(str(path)) as w:
        w.write(mol, meta={'ID': '7'}, title='methane')
    with SDFRead(str(path)) as f:
        again = f.read_record()
    assert again.atom_count == 1
    assert again.implicit_h_of(next(iter(again.atom_numbers))) == 4
    assert again.meta == {'ID': '7'}
    assert again.title == 'methane'


def test_append_adds_a_record_instead_of_truncating(tmp_path):
    path = tmp_path / 'out.sdf'
    m = mol(_sdf(_METHANE))
    with SDFWrite(str(path)) as w:
        w.write(m)
    with SDFWrite(str(path), append=True) as w:
        w.write(m)
    with SDFRead(str(path)) as f:
        assert len(f.read()) == 2


def test_the_unknown_set_survives_a_write_with_nothing_asked_of_the_caller(corpus):
    """An unknown count round trips through the file surface with nothing asked of the caller: the
    writer must not state ``IMPL_H0`` and turn "nobody said" into "there are none"."""
    checked = 0
    for name, records in corpus.items():
        for record in records:
            with SDFRead(StringIO('\n'.join(record) + '\n')) as f:
                mol = f.read_record()
                unknown, sgroups = f.unknown_hydrogens, f.sgroups
            if not unknown or mol.aromatic_bond_count:
                continue  # aromatic records cannot be written at all; that is test_v2000's business
            buffer = StringIO()
            with SDFWrite(buffer) as w:
                w.write(mol, sgroups=sgroups)
            buffer.seek(0)
            with SDFRead(buffer) as f:
                again = f.read_record()
                assert set(f.unknown_hydrogens) == set(unknown), name
            assert again.unknown_h_count == mol.unknown_h_count, name
            checked += 1
    assert checked, 'no writable record in the corpus has an unknown count; this test proved nothing'


def test_esdfwrite_writes_a_v3000_record(tmp_path):
    path = tmp_path / 'out.sdf'
    m = mol(_sdf(_METHANE))
    with ESDFWrite(str(path)) as w:
        w.write(m, title='methane')
    assert V3000_STAMP in path.read_text(encoding='utf-8')
    with SDFRead(str(path)) as f:
        record = f.read_record()
        assert f.version == V3000_STAMP
    assert record.atom_count == 1


def test_both_writers_write_an_aromatic_bond_as_bond_type_four():
    """Bond type 4 is the CTfile column for an aromatic bond, so both writers emit it -- V2000 as
    ``  1  2  4`` and V3000 (``ESDFWrite``) as ``M  V30 1 4 1 2``.  Tested as a pair: one version
    changing without the other is this package's signature failure."""
    m = mol('\n'.join(_AROMATIC))
    assert m.aromatic_bond_count == 1

    buf = StringIO()
    with SDFWrite(buf) as w:
        w.write(m)
    assert '  1  2  4  0  0  0  0' in buf.getvalue().splitlines(), buf.getvalue()

    buf = StringIO()
    with ESDFWrite(buf) as w:
        w.write(m)
    assert 'M  V30 1 4 1 2' in buf.getvalue().splitlines(), buf.getvalue()


def test_the_writer_returns_the_log_rather_than_swallowing_it():
    """A writer that had to decide something returns the log rather than swallowing it."""
    with SDFWrite(StringIO()) as w:
        assert isinstance(w.write(mol('\n'.join(_METHANE))), list)


# the corpus

def test_the_corpus_reads_through_the_stream_with_nothing_failed(corpus, root):
    """Through the file rather than a record list: the only way framing, sniffing and the stream are
    checked end to end."""
    total = 0
    for name in corpus:
        path = root / 'test' / name
        with SDFRead(str(path)) as f:
            mols = f.read()
        assert not f.failed, f'{name}: {f.failed}'
        assert len(mols) == len(corpus[name]), name
        total += len(mols)
    assert total == 512, total


# a failure that is not the file's

def test_a_parser_bug_lands_in_failed_like_any_other_bad_record(monkeypatch):
    """Any exception from one record is filed rather than raised, so the tail of the file survives a
    parser bug; the exception object is kept for whoever needs to tell the two cases apart."""
    from .. import _stream

    def boom(lines, log, **kwargs):
        raise RuntimeError('parser tripped')

    monkeypatch.setattr(_stream, 'parse_record', boom)
    with SDFRead(StringIO(_sdf(_METHANE, _METHANE))) as f:
        assert f.read() == []
    assert [type(x.error).__name__ for x in f.failed] == ['RuntimeError', 'RuntimeError']
    assert [x.position for x in f.failed] == [0, 1]


def test_the_writer_declares_no_slot_it_never_assigns():
    """The version is the class-level `_stamp` the subclass overrides, not an instance slot."""
    assert SDFWrite.__slots__ == ()
    assert not hasattr(SDFWrite(StringIO()), '_version')
    assert SDFWrite._stamp != ESDFWrite._stamp


def test_the_reader_yields_molecules_carrying_their_meta(tmp_path):
    path = tmp_path / 'a.sdf'
    path.write_text(_sdf(_METHANE + ['>  <NAME>', 'methane', '']))
    with SDFRead(path) as f:
        mol = next(iter(f))
        assert f.meta == mol.meta
    assert mol.meta == {'NAME': 'methane'}


def test_the_readers_title_is_the_molecules_title(tmp_path):
    """One name line, one spelling, one type -- the reader delegates rather than keeping a copy."""
    path = tmp_path / 'a.sdf'
    path.write_text(_sdf(_METHANE))
    with SDFRead(path) as f:
        mol = next(iter(f))
        assert f.title == mol.title == 'methane'


def test_a_write_round_trips_meta_with_no_keyword(tmp_path):
    """The done-when: out on a write with nothing passed."""
    mol = parse_record(_METHANE + ['>  <NAME>', 'methane', ''])
    path = tmp_path / 'b.sdf'
    with SDFWrite(path) as w:
        w.write(mol)
    with SDFRead(path) as f:
        assert next(iter(f)).meta == {'NAME': 'methane'}


def test_the_writer_and_the_facade_state_the_same_default():
    """Both spell "I did not say" as `None`; a caller switching between them must not change meaning."""
    m = read_smiles('CCO')
    m.meta['ACTIVITY'] = '5.0'
    buf = StringIO()
    with SDFWrite(buf) as f:
        f.write(m)
    assert ('>  <ACTIVITY>' in buf.getvalue()) is ('>  <ACTIVITY>' in mol(m))


def test_the_writer_writes_no_data_fields_when_told_none():
    m = read_smiles('CCO')
    m.meta['ACTIVITY'] = '5.0'
    buf = StringIO()
    with SDFWrite(buf) as f:
        f.write(m, meta={})
    assert 'ACTIVITY' not in buf.getvalue()


def test_esdfwrite_writes_the_molecules_own_data_fields_too():
    """The V3000 writer is the same writer with another stamp; the default must not differ."""
    m = read_smiles('CCO')
    m.meta['ACTIVITY'] = '5.0'
    buf = StringIO()
    with ESDFWrite(buf) as f:
        f.write(m)
    assert '>  <ACTIVITY>' in buf.getvalue()
