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
"""The file surface: :class:`SDFRead`, :class:`SDFWrite`, :class:`ESDFWrite`.

Streams wrapped around :func:`~._sdf.parse_record` and :func:`~._sdf.emit_record`; no chemistry and no
column offsets here.  A reader never raises for a bad record and never skips it silently: it lands in
:attr:`SDFRead.failed`.  The title, the data fields and the log are the container's own; what is left
on the stream is the record's framing -- version, program, comment, S-group store, unknown hydrogens.
"""

from pathlib import Path

from ._rxn import parse_rxn_record
from ._sdf import emit_record, parse_data_fields, parse_record, split_records
from ._v2000 import V2000_STAMP
from ._v3000 import V3000_STAMP
from ...core import LogRecord, REPAIRED


__all__ = ['SDFRead', 'SDFWrite', 'ESDFWrite', 'FailedRecord']


class FailedRecord:
    """A record the reader could not parse: its ``lines``, the ``error``, and its ``position``.

    Kept rather than raised and rather than dropped -- see :attr:`SDFRead.failed`.
    """
    __slots__ = ('position', 'lines', 'error')

    def __init__(self, position, lines, error):
        self.position = position
        self.lines = list(lines)
        self.error = error

    @property
    def text(self):
        return '\n'.join(self.lines)

    def __repr__(self):
        return f'FailedRecord({self.position}, {type(self.error).__name__}: {self.error!s:.60})'


class _FileBacked:
    """Path-or-buffer handling and the context manager, shared by the reader and the writers.

    Accepts ``str``, :class:`pathlib.Path` or any object with the right method, and remembers which it
    was: a file this class opened is closed on exit, a buffer the caller passed only with
    ``close(force=True)``.
    """
    __slots__ = ('_file', '_is_buffer')

    def _open(self, file, mode, probe):
        # `errors='surrogateescape'` in BOTH directions, and it is what makes the round trip byte for
        # byte: a name line is not required to be UTF-8, a reader may not refuse a record for it, and a
        # writer must put the byte back.  A buffer the caller opened is the caller's policy.
        #
        # `encoding='utf-8'` states the codec that pair works over: the locale's would decide which byte
        # becomes which character, so one file would read as a different title on a cp1252 host and
        # `surrogateescape` would hand back a byte the file never held.  `newline='\n'` when writing is
        # the same rule in the other direction -- what a record's bytes are may not be a property of the
        # host, and text mode writes CRLF on Windows.
        opened = {'encoding': 'utf-8', 'errors': 'surrogateescape'}
        if 'r' not in mode:
            opened['newline'] = '\n'
        if isinstance(file, str):
            self._file = open(file, mode, **opened)
            self._is_buffer = False
        elif isinstance(file, Path):
            self._file = file.open(mode, **opened)
            self._is_buffer = False
        elif hasattr(file, probe):
            self._file = file
            self._is_buffer = True
        else:
            raise TypeError(f'invalid file: expected a path or an object with .{probe}(), '
                            f'got {type(file).__name__}')

    def close(self, force=False):
        """Close the file.  A buffer the caller opened is closed only with `force`."""
        if not self._is_buffer or force:
            self._file.close()

    def __enter__(self):
        return self

    def __exit__(self, _type, value, traceback):
        self.close()


class SDFRead(_FileBacked):
    """MDL SDF/molfile reader.  Iterate it for molecules; ``with`` works; a path or a buffer works.

    ::

        with SDFRead('input.sdf') as f:
            for mol in f:
                ...
            if f.failed:
                print(f'{len(f.failed)} unparsable record(s)')

    :attr:`meta`, :attr:`log` and :attr:`title` are the container's own and are reachable from it;
    :attr:`sgroups`, :attr:`version`, :attr:`program`, :attr:`comment` and :attr:`unknown_hydrogens`
    are the record's framing, which no container holds, so they are reader state for the record most
    recently returned.
    """
    __slots__ = ('_records', '_position', '_record', '_header', 'failed', 'ignore_stereo')

    def __init__(self, file, *, ignore_stereo=False):
        """:param ignore_stereo: skip the stereo step; the constitution is read either way."""
        self._open(file, 'r', 'read')
        self._records = split_records(self._file)
        self._position = -1
        self._record = None
        #: The current record's framing facts; see the properties below.
        self._header = {}
        #: ``[FailedRecord, ...]`` -- every record that could not be parsed, in file order.
        self.failed = []
        self.ignore_stereo = ignore_stereo

    @property
    def record(self):
        """The container most recently read, or ``None``.  A `MoleculeContainer`, or a
        `ReactionContainer` for a record that held a ``$RXN``."""
        return self._record

    @property
    def meta(self):
        """``{name: value}`` of the current record's data fields -- the container's own `meta`."""
        return self._record.meta if self._record is not None else {}

    @property
    def sgroups(self):
        """The current record's :class:`~._sgroup.SGroupStore`.  ``mol.sgroups`` is the same records as
        dicts; this is the parsed store, for a caller re-emitting them."""
        return self._header.get('sgroups')

    @property
    def log(self):
        """Every recovery made while reading the current record -- the container's own `log`."""
        return self._record.log if self._record is not None else []

    @property
    def title(self):
        """The record's first line -- ``mol.title``, and the same ``str``."""
        return self._record.title if self._record is not None else ''

    @property
    def version(self):
        """The CTAB version actually read, which is not always the one stamped."""
        return self._header.get('version')

    @property
    def program(self):
        """Line 2 of the molfile: the program that wrote it."""
        return self._header.get('program', '')

    @property
    def comment(self):
        """Line 3 of the molfile."""
        return self._header.get('comment', '')

    @property
    def unknown_hydrogens(self):
        """Stable ids of the current record's atoms whose implicit hydrogen count nobody stated."""
        return self._header.get('unknown_hydrogens', ())

    def tell(self):
        """The index of the record most recently read; ``-1`` before the first."""
        return self._position

    def read_record(self):
        """The next record's container.  Raises ``StopIteration`` at end of file.

        Unparsable records are recorded in :attr:`failed` and skipped, so this returns the next record
        that *is* parsable rather than propagating the failure.  Ask :attr:`failed` afterwards.

        A record whose first line is ``$RXN`` comes back as a `ReactionContainer` instead -- an SD file
        should not hold one, and reading it beats refusing it.
        """
        for lines in self._records:
            self._position += 1
            self._header = {}
            try:
                if lines and lines[0].startswith('$RXN'):
                    self._record = self._read_embedded_reaction(lines)
                    return self._record
                self._record = parse_record(lines, [], ignore_stereo=self.ignore_stereo,
                                            header=self._header)
            except Exception as e:  # noqa: one record's failure, of any kind, must not end the file
                # A parser bug and a malformed record are filed identically; the exception is kept
                # for whoever needs to tell them apart.
                self._record = None
                self._header = {}
                self.failed.append(FailedRecord(self._position, lines, e))
                continue
            return self._record
        raise StopIteration

    def _read_embedded_reaction(self, lines):
        """A record holding a ``$RXN`` block, read as a `ReactionContainer`.

        Read rather than filed as a :class:`FailedRecord`, with a log line, and through the same
        :func:`~._rxn.parse_rxn_record` an RDfile's ``$RFMT`` goes through.
        """
        own = [LogRecord('sdf:record-holds-rxn', (),
                         'record holds a $RXN block, which an SD file should not; read as a reaction',
                         REPAIRED)]
        # `parse_data_fields` skips to the FIRST `M  END` and a reaction has one per component, so the
        # data-field block can only begin after the LAST one: cut at the first column-0 `>` after it.
        # A `$MOL` component's title line (line 2 of a molfile) is arbitrary text and may start with
        # `>`, so scanning from the top would cut on a component titled `> product name`.
        last_mend = max((i for i, x in enumerate(lines) if x.startswith('M  END')), default=-1)
        cut = next((i for i, x in enumerate(lines)
                    if i > last_mend and x.startswith('>')), len(lines))
        fields = parse_data_fields(['M  END', *lines[cut:]], own) if cut < len(lines) else {}
        reaction = parse_rxn_record(lines[:cut], fields, [], ignore_stereo=self.ignore_stereo,
                                    header=self._header)
        # `SDFRead.log` IS `record.log`, so what this method observed about the record has nowhere else
        # to go: the rescue and the data-field lines belong to the reaction it produced.
        reaction.log.absorb('read', own)
        return reaction

    def read_structure(self):
        """The next molecule -- or the reaction of a record that held a ``$RXN``, which an SD file
        should not.  Same call as :meth:`read_record`, kept as the spelling a streaming caller reaches
        for."""
        return self.read_record()

    def read(self, amount=None):
        """The whole file as a list of molecules, or the next `amount` of them."""
        out = []
        while amount is None or len(out) < amount:
            try:
                out.append(self.read_structure())
            except StopIteration:
                break
        return out

    def __iter__(self):
        return self

    def __next__(self):
        return self.read_structure()


class SDFWrite(_FileBacked):
    """MDL V2000 SDF writer.  ``write(mol)`` per record; ``with`` works; a path or a buffer works.

    A molecule holding an aromatic bond is refused rather than kekulised on the caller's behalf; the
    message names ``kekule()``.  Use :class:`ESDFWrite` for V3000.
    """
    __slots__ = ()
    _stamp = V2000_STAMP

    def __init__(self, file, *, append=False):
        self._open(file, 'a' if append else 'w', 'write')

    def write(self, mol, *, sgroups=None, meta=None, title=None):
        """Write one record.  Returns the log of anything the writer had to decide.

        :param sgroups: an :class:`~._sgroup.SGroupStore`, or ``None`` for the molecule's own
        :param title: a replacement name line, or ``None`` for the molecule's own
        :param meta: ``{name: value}``, or any iterable of pairs -- the SDF data fields

        ``title`` and ``sgroups`` default to ``None``, meaning "what the molecule holds" rather than
        "empty", so ``write(mol)`` in a read loop keeps the title, the S-groups and the atom aliases.
        ``sgroups=SGroupStore()`` still writes none: "I did not say" and "I said none" differ.

        ``meta`` defaults the same way: ``None`` is the molecule's own data fields, so ``write(mol)`` in
        a read loop round-trips them.  ``meta={}`` writes none.
        """
        lines, log = emit_record(mol, sgroups, meta, version=self._stamp, title=title,
                                 separator=True)
        self._file.write('\n'.join(lines))
        self._file.write('\n')
        return log


class ESDFWrite(SDFWrite):
    """MDL V3000 SDF writer.  Same surface as :class:`SDFWrite`, extended CTAB.

    V3000 is what to reach for when a structure does not fit V2000's fixed columns -- more than 999
    atoms or bonds, a coordinate outside the 10-character column, or an AND/OR stereo group, which
    V2000 cannot spell at all.  A large charge is not one of them: V2000 writes 0 in the ccc column
    and the truth in `M  CHG`.  `needs_v3000()` is that list as a predicate.

    It does not buy an aromatic bond: order ``4`` in a structure record states a query, so
    ``emit_v3000`` refuses it exactly as V2000 does and points at ``kekule()``.
    """
    __slots__ = ()
    _stamp = V3000_STAMP
