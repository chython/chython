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
"""RDfile: records of molecules OR reactions, with ``$DTYPE``/``$DATUM`` metadata.

Framing is state-aware where SD framing is not: a ``$DATUM`` value is free text that may span
physical lines, so a value line beginning ``$RFMT`` is a file somebody really has.  The continuation
rule (CTfile specification p.46) applies to ``$DTYPE``/``$DATUM`` logical lines only -- V3000's
``M  V30`` lines can legitimately exceed column 80:

1. a physical line whose **content** (trailing whitespace excluded) is ``>= 80`` characters
   **positionally continues** onto the next, concatenated with no separator; a trailing ``+`` is
   data, not a marker;
2. a non-``$`` line following a ``$DATUM`` is a **permissive continuation**, joined with a newline,
   because vendors wrap short lines without 80-character fill;
3. a ``$``-led line ends the value unless positional continuation is open, in which case it is value
   text and its absorption is logged;
4. any other ``$``-led line ends the logical line and is dropped with one log line;
5. a line kept inside the logical line with no open ``$DATUM`` value to join is malformed: dropped
   and logged;
6. a second ``$DATUM`` under one ``$DTYPE`` names nothing, so the **first** value is kept.

``$MIREG``/``$MEREG``/``$RIREG``/``$REREG`` carry ``unsupported: ``.  ``$DATM`` does not: on line 2
it is the header timestamp the specification puts there, stored verbatim on :attr:`RDFRead.date`,
while in a data-field tail it is a header keyword where the format has no room for it.
"""

from collections.abc import Mapping
from pathlib import Path
from time import strftime

from ._rxn import emit_rxn, parse_rxn_record
from ._sdf import emit_record, parse_record
from ._stream import _FileBacked, FailedRecord
from ._v2000 import V2000_STAMP
from ._v3000 import V3000_STAMP
from ...core import LogRecord, LOST, REPAIRED, ReactionContainer


__all__ = ['RDF_HEADER', 'split_rdf_records', 'parse_rdf_fields', 'parse_rdf_record',
           'RDFRead', 'RDFWrite', 'ERDFWrite']


RDF_HEADER = '$RDFILE 1'
_MOL_RECORD = '$MFMT'
_RXN_RECORD = '$RFMT'
_RECORD_TAGS = (_MOL_RECORD, _RXN_RECORD)
_DTYPE = '$DTYPE'
_DATUM = '$DATUM'
#: The file-level timestamp.  On line 2 of the header it is stored on the reader; a vendor also emits
#: one in the ``$DTYPE``/``$DATUM`` tail, where the format has no room for it and it is reported.
_DATM = '$DATM'
#: ``$MIREG``/``$MEREG``/``$RIREG``/``$REREG`` may follow the record tag on its own line, and a
#: vendor also puts one in the ``$DTYPE``/``$DATUM`` tail.  Both positions read the same way.
_REGISTRY_TAGS = ('$MIREG', '$MEREG', '$RIREG', '$REREG')
#: Physical content length at which the spec mandates positional continuation.
_CONTINUATION_LENGTH = 80
#: Columns a ``$DTYPE ``/``$DATUM `` keyword and its separating space occupy on the first physical line
#: of a logical line.  Both keywords are six characters, so one constant serves both.
_KEYWORD_WIDTH = len(_DATUM) + 1


def _registry_message(line):
    """The one wording for a registry cross-reference, wherever in the record it was found.

    An ``unsupported: `` line: the format has the construct and this layer has nothing to store it in,
    which position within the record does not change.
    """
    return f'unsupported: {line[:6].strip()} registry cross-reference not stored: {line!r:.60}'


def _datm_payload(line):
    """The raw text after ``$DATM``, verbatim: ``'01/02/17 17:17'``.

    Unparsed on purpose -- whether ``01/02/17`` is January the second or the first of February, and
    whether the year is 2017 or 1917, is a guess this layer has no basis for.
    """
    return line[len(_DATM):].strip()


def _log_absorbed_keyword(log, line, target):
    """Report a ``$``-led line that a positional continuation swallowed, whatever swallowed it.

    *target* names what absorbed it -- a wrapped ``$DTYPE`` name or an open ``$DATUM`` value.  One
    wording for both, because it means the same thing in each: the file's fill to column 80 was
    accidental and a line that looks like a keyword was read as text.  A continuation not starting
    with ``$`` is ordinary wrapped content and says nothing worth a line.
    """
    if line.startswith('$'):
        log.append(LogRecord('rdf:absorbed-keyword', (),
                             f'absorbed a keyword as a positional continuation of a long {target}: '
                             f'{line!r:.40}'))


def _continues(line):
    """``True`` when this physical line mandates positional continuation onto the next.

    Length is measured without trailing whitespace: a vendor that space-pads a short value to column
    80 arms no continuation, while genuine fill carries content to column 80 without padding.
    """
    return len(line.rstrip()) >= _CONTINUATION_LENGTH


class _DataFieldState:
    """Shared tracker for whether a ``$DTYPE``/``$DATUM`` logical line is open.

    Positional continuation is armed only while such a line is open, and both ``split_rdf_records``
    and ``parse_rdf_fields`` read that rule from here rather than each deriving it.  The logical line
    is entered by a ``$DTYPE`` or ``$DATUM`` (or a continuation of one) and left by any ``$``-led line
    that is neither and is not a positional continuation; a non-``$`` line after a ``$DATUM`` is a
    permissive continuation and stays inside it.  An orphan ``$DATUM`` opens the line too, so its
    positional continuation is consumed without a second message.

    **Reset between records** (call ``reset()``).
    """

    __slots__ = ('in_field', 'prev_long')

    def __init__(self):
        self.in_field = False   # inside a $DTYPE/$DATUM logical line
        self.prev_long = False  # previous line was >= _CONTINUATION_LENGTH content characters

    def reset(self):
        self.in_field = False
        self.prev_long = False

    @property
    def continuation_open(self):
        """``True`` when the next line is the positional continuation of the previous one.

        *The* arming predicate, written once: the previous data-field line's content reached column 80
        **and** a ``$DTYPE``/``$DATUM`` logical line is open.  Both callers read it rather than
        re-deriving the conjunction.
        """
        return self.prev_long and self.in_field

    def advance(self, line):
        """Update state for *line* (already stripped of ``\r``/``\n``).

        Returns ``(was_continuation, was_permissive)`` -- the two ways *line* belongs to the logical
        line that was open **before** this call, at most one of them ``True``:

        - ``was_continuation``: *line* is the positional-continuation target of the previous long
          line, whatever it starts with.  This is ``continuation_open`` read before the update.
        - ``was_permissive``: *line* is a non-``$``-led line inside an open logical line, which a
          vendor wrote by wrapping a value without 80-character fill.

        Recognising ``$DTYPE`` and ``$DATUM`` is the caller's job; everything else about where a line
        belongs is decided here.
        """
        was_continuation = self.continuation_open
        was_permissive = False
        long = _continues(line)

        if was_continuation:
            # Positional continuation: the open logical line absorbs this line, whatever it is.
            self.prev_long = long
            # in_field stays True
        elif line.startswith(_DTYPE) or line.startswith(_DATUM):
            self.in_field = True
            self.prev_long = long
        elif self.in_field and not line.startswith('$'):
            # Permissive continuation of a $DATUM value: non-$-led line keeps the field open.
            was_permissive = True
            self.prev_long = long
        else:
            # A $-led line that is not $DTYPE/$DATUM, or any line when not in a field.
            self.in_field = False
            self.prev_long = False

        return was_continuation, was_permissive


def split_rdf_records(lines, log=None, *, header=None):
    """Yield ``(tag, [record lines])`` for each ``$MFMT``/``$RFMT`` record.  The tag is excluded.

    File-level ``$RDFILE``/``$DATM`` lines before the first record are consumed and not yielded.  Pass
    a dict as *header* to keep what they said: the timestamp arrives as ``header['date']``, verbatim
    and unparsed, and :attr:`RDFRead.date` is that key.

    Delimiter recognition is suppressed **only** while a positional continuation of a
    ``$DTYPE``/``$DATUM`` logical line is open, so a long ``M  V30`` line does not swallow the
    following ``$MFMT``; every suppression is logged.  The region before the first record tag has no
    logical line to continue into, so a long ``$DTYPE`` there arms nothing.
    """
    log = [] if log is None else log
    tag = None
    current = []
    state = _DataFieldState()
    for raw in lines:
        line = raw.rstrip('\r\n')
        continuation_open = state.continuation_open
        if not continuation_open and line.startswith(_RECORD_TAGS):
            if tag is not None:
                yield tag, current
            tag = line[:5]
            current = []
            state.reset()
            # A registry number on the tag line is metadata we do not model; say so once per record.
            rest = line[5:].strip()
            if rest:
                log.append(LogRecord('rdf:registry-reference', (),
                                     f'unsupported: {tag} registry reference {rest!r:.30} not stored',
                                     LOST))
            continue
        if tag is None:
            # The header region holds no logical line, so a stray $DTYPE here arms nothing and the
            # next $MFMT stays a delimiter.
            state.reset()
            if line.startswith(_DATM):
                # Where the specification puts the timestamp, so it is stored rather than reported.
                # A second $DATM in one header is malformed; the last one wins and says so.
                if header is not None:
                    if 'date' in header:
                        log.append(LogRecord('rdf:duplicate-datm-header', (),
                                             f'a second {_DATM} in the file header; the later timestamp is '
                                             f'kept: {line!r:.60}', REPAIRED))
                    header['date'] = _datm_payload(line)
            elif not (line.startswith(RDF_HEADER) or not line.strip()):
                log.append(LogRecord('rdf:pre-record-line', (),
                                     f'line before the first record, ignored: {line!r:.60}', LOST))
            continue
        # Inside a record body: check for a suppressed record tag
        if continuation_open and line.startswith(_RECORD_TAGS):
            log.append(LogRecord('rdf:record-tag-suppressed', (),
                                 f'record tag suppressed inside a positional continuation: {line!r:.40}'))
        current.append(line)
        state.advance(line)
    if tag is not None:
        yield tag, current


def parse_rdf_fields(lines, log=None):
    """``{name: value}`` from a record's ``$DTYPE``/``$DATUM`` tail, in file order.

    A repeated name merges its value lines and is reported, exactly as ``parse_data_fields`` does it:
    one value per name is what a mapping holds.  The ``$DATUM`` read is a **prefix slice, not** ``lstrip``, which
    strips a character set and would eat ``$DATUMMUTADAT`` whole; a trailing ``+`` is data, the format
    having no marker-character continuation.

    Continuation state is ``_DataFieldState``'s, shared with ``split_rdf_records``: a positional
    continuation is concatenated onto the ``$DTYPE`` name while no ``$DATUM`` has arrived, otherwise
    into the open value, and a ``$``-led line absorbed either way is logged.  A ``$``-led line that is
    *not* absorbed ends the logical line and is dropped with one line, as is a line with no open value
    to join -- a registry reference prefixed ``unsupported: ``, the rest unprefixed, a broken file not
    being a construct chython declines to model.
    """
    log = [] if log is None else log
    fields = {}
    name = None
    current_lines = None
    state = _DataFieldState()

    def flush():
        if name in fields:
            log.append(LogRecord('rdf:duplicate-dtype', (),
                                 f'{_DTYPE} {name!r:.30} appears twice and one value per name is what `meta` '
                                 f'holds; the values are merged', REPAIRED))
        fields.setdefault(name, []).extend(current_lines if current_lines is not None else [])

    for line in lines:
        was_continuation, was_permissive = state.advance(line)

        if was_continuation:
            # Positional continuation: which target absorbs it decides the wording and nothing else.
            if current_lines is None and name is not None:
                # Continuing a $DTYPE name line.
                _log_absorbed_keyword(log, line, f'{_DTYPE} name')
                name = name + line
            elif current_lines is not None:
                _log_absorbed_keyword(log, line, f'{_DATUM} value')
                current_lines[-1] = current_lines[-1] + line
            # else: orphan $DATUM continuation — consumed with no extra log line (see docstring).
        elif line.startswith(_DTYPE):
            if name is not None:
                flush()
            name = line[len(_DTYPE):].strip()
            current_lines = None
        elif line.startswith(_DATUM):
            if name is None:
                log.append(LogRecord('rdf:orphan-datum', (),
                                     f'{_DATUM} with no preceding {_DTYPE}, ignored: {line!r:.60}', LOST))
            elif current_lines is not None:
                # A second $DATUM under one $DTYPE names nothing, so the field keeps the value the
                # $DTYPE named and this one goes -- reported, never silently.
                log.append(LogRecord('rdf:duplicate-datum', (),
                                     f'a second {_DATUM} under one {_DTYPE}: the first value is kept and '
                                     f'this one dropped: {line!r:.60}', LOST))
            else:
                # Prefix slice, not `lstrip`, which strips a character set and would eat the leading
                # letters of a value like `MUTADAT`.
                current_lines = [line[len(_DATUM):].strip()]
        elif line.startswith(_DATM):
            # A header keyword in a position the format has no room for, and NOT `unsupported: `: the
            # header's own timestamp is stored, so nothing is missing here and this is a report about
            # the file rather than about a gap in chython.
            log.append(LogRecord('rdf:misplaced-datm', (),
                                 f'{_DATM} outside the file header, ignored: {line!r:.60}', LOST))
        elif was_permissive and current_lines is not None:
            # Permissive continuation: vendor wrapped without 80-character fill.
            current_lines.append(line)
        elif line.startswith(_REGISTRY_TAGS):
            # A registry cross-reference reads the same here as in the structure body.
            log.append(LogRecord('rdf:registry-reference', (), _registry_message(line), LOST))
        elif line.startswith('$'):
            # A $-led line that is no keyword we know and is not a positional continuation.
            # The data-field logical line has ended; log and drop the line.
            log.append(LogRecord('rdf:unrecognised-keyword', (),
                                 f'unrecognised field keyword, ignored: {line!r:.60}', LOST))
        else:
            # No open $DATUM value to join: between a $DTYPE and its $DATUM, after an orphan $DATUM,
            # or a CTAB body line leaking in.  Nowhere to put it, so it goes -- but never silently.
            log.append(LogRecord('rdf:line-outside-datum', (),
                                 f'line outside any $DATUM value, ignored: {line!r:.60}', LOST))

    if name is not None:
        flush()
    return {k: '\n'.join(v) for k, v in fields.items()}


def parse_rdf_record(tag, lines, log=None, *, ignore_stereo=False, header=None):
    """A `MoleculeContainer` for ``$MFMT``, a `ReactionContainer` for ``$RFMT``.

    The structure runs to the first ``$DTYPE``; everything from there is metadata.  A record with no
    ``$DTYPE`` is all structure, the common case in a reaction-only RDfile.

    A registry cross-reference in the structure body is removed from it and logged as unsupported, with
    the same wording ``parse_rdf_fields`` gives one in the metadata tail.  A molecule record's
    ``$DTYPE`` pairs land on ``mol.meta``, the same place a reaction record's do.  A ``$MFMT`` record
    whose body starts with ``$RXN`` is rescued and logged unprefixed, the file being wrong.

    `header`, when given a dict, is filled by whichever parser ran -- see :func:`~._sdf.parse_record`
    and :func:`~._rxn.parse_rxn_record`.

    The registry references and the ``$DTYPE``/``$DATUM`` lines are about the record, so they land on
    whichever container the record produced as well as on `log`.
    """
    log = [] if log is None else log
    # This function's own lines, held aside because the container that will carry them is built below.
    own = []
    cut = next((i for i, x in enumerate(lines) if x.startswith(_DTYPE)), len(lines))
    body = []
    for x in lines[:cut]:
        if x.startswith(_REGISTRY_TAGS):
            own.append(LogRecord('rdf:registry-reference', (), _registry_message(x), LOST))
        else:
            body.append(x)
    fields = parse_rdf_fields(lines[cut:], own)

    # A $MFMT slot that holds a $RXN block is a malformed RDfile; read it as a reaction rather than
    # filing a FailedRecord.  The log line is how a caller learns their input mixed the two kinds.
    rescued = tag != _RXN_RECORD and bool(body) and body[0].startswith('$RXN')
    if rescued:
        own.append(LogRecord('rdf:mfmt-contains-rxn', (),
                             '$MFMT record holds a $RXN block, which an RDfile $MFMT should not; '
                             'read as a reaction', REPAIRED))
    # Before the structure parser writes, so the caller's list stays in file order.
    log.extend(own)

    if tag == _RXN_RECORD or rescued:
        container = parse_rxn_record(body, fields, log, ignore_stereo=ignore_stereo, header=header)
    else:
        container = parse_record(body, log, ignore_stereo=ignore_stereo, header=header)
        if fields:
            # RDfile molecule records carry $DTYPE, not `>  <name>` fields, so they arrive separately
            # and land in the same place.
            container.meta.update(fields)
    container.log.absorb('read', own)
    return container


class RDFRead(_FileBacked):
    """MDL RDfile reader.  Iterate it for molecules **and** reactions; ``with`` works; path or buffer.

    ::

        with RDFRead('input.rdf') as f:
            for structure in f:            # a MoleculeContainer or a ReactionContainer
                ...
            if f.failed:
                print(f'{len(f.failed)} unparsable record(s)')

    An RDfile interleaves ``$MFMT`` and ``$RFMT`` by design, so iteration yields whichever the record
    held; ``isinstance`` tells them apart.

    :attr:`meta`, :attr:`log` and :attr:`title` are the container's own and are reachable from it;
    :attr:`version`, :attr:`program` and :attr:`comment` are the record's framing, which no container
    holds, so they are reader state for the record most recently returned.

    :attr:`date` (the file's own ``$DATM``, never copied into a container's metadata) and
    :attr:`file_log` (framing damage found between records) belong to the FILE, so they are reader
    state that accumulates as records are consumed.
    """
    __slots__ = ('_records', '_position', '_record', '_header', '_record_header', 'failed',
                 'file_log', 'ignore_stereo')

    def __init__(self, file, *, ignore_stereo=False):
        """:param ignore_stereo: skip the stereo step; the constitution is read either way."""
        self._open(file, 'r', 'read')
        # Filled by the splitter as it lazily consumes the first lines, so `date` is None until the
        # first record is asked for.
        self._header = {}
        #: The current record's framing facts -- version, program, comment.  Separate from
        #: `_header`, which is the FILE's and is where `date` comes from.
        self._record_header = {}
        #: Framing damage, in file order: a stray line before the first record, a registry reference
        #: on a tag line, a record tag swallowed by a positional continuation.  Separate from a
        #: record's own log because a splitter decision is about where records begin and end, so
        #: attributing one to a single record would be a guess.
        self.file_log = []
        self._records = split_rdf_records(self._file, self.file_log, header=self._header)
        self._position = -1
        self._record = None
        #: ``[FailedRecord, ...]`` -- every record that could not be parsed, in file order.
        self.failed = []
        self.ignore_stereo = ignore_stereo

    # ------------------------------------------------------------------ the current record's context

    @property
    def record(self):
        """The container most recently read, or ``None`` -- a molecule for ``$MFMT``, a reaction for
        ``$RFMT``."""
        return self._record

    @property
    def meta(self):
        """``{name: value}`` of the current record's ``$DTYPE``/``$DATUM`` pairs -- the container's
        own `meta`."""
        return self._record.meta if self._record is not None else {}

    @property
    def log(self):
        """Every recovery made while reading the current record -- the container's own `log`."""
        return self._record.log if self._record is not None else []

    @property
    def title(self):
        """The record's name line -- ``container.title``, and the same ``str``."""
        return self._record.title if self._record is not None else ''

    @property
    def version(self):
        """The CTAB version actually read, which is not always the one stamped."""
        return self._record_header.get('version')

    @property
    def program(self):
        """Line 2 of the record header: the program that wrote it."""
        return self._record_header.get('program', '')

    @property
    def comment(self):
        """Line 3 of the record header."""
        return self._record_header.get('comment', '')

    @property
    def date(self):
        """The file header's ``$DATM`` payload, **verbatim and unparsed**, or ``None``.

        It looks like ``'01/02/17 17:17'``; neither the century nor the day/month order is decidable
        from the string.  ``None`` until the first record has been read, the header being consumed
        lazily with it, and ``None`` afterwards for a file carrying no timestamp.
        """
        return self._header.get('date')

    def tell(self):
        """The index of the record most recently read; ``-1`` before the first."""
        return self._position

    # ----------------------------------------------------------------------------------- reading

    def read_record(self):
        """The next record.  Raises ``StopIteration`` at end of file, as :class:`~._stream.SDFRead`
        does.

        Unparsable records are recorded in :attr:`failed` and skipped, so this returns the next record
        that *is* parsable rather than propagating the failure.
        """
        for tag, lines in self._records:
            self._position += 1
            self._record_header = {}
            try:
                self._record = parse_rdf_record(tag, lines, [], ignore_stereo=self.ignore_stereo,
                                                header=self._record_header)
            except Exception as e:
                # One unreadable record never costs the file.  A parser bug of ours is filed the same
                # way a malformed record is: a caller cannot act on the difference.
                self._record = None
                self._record_header = {}
                self.failed.append(FailedRecord(self._position, lines, e))
                continue
            return self._record
        raise StopIteration

    def read_structure(self):
        """The next molecule or reaction.  Same call as :meth:`read_record`, kept as the spelling a
        streaming caller reaches for."""
        return self.read_record()

    def read(self, amount=None):
        """The whole file as a list of molecules and reactions, or the next `amount` of them."""
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


def _wrap_logical_line(text, log, what):
    """*text* as physical lines the reader rejoins, byte for byte.

    Not cosmetic wrapping.  The reader concatenates a physical line whose content reaches column 80
    onto the next with no separator (specification p.46), so a logical line written long in one
    physical line arms that rule against whatever follows the record -- an unwrapped 200-character
    ``$DATUM`` swallows the next ``$DTYPE`` or ``$MFMT`` -- while one wrapped SHORT is rejoined with a
    ``\n`` between the pieces, a short non-``$`` line being a permissive continuation.  So the chunks
    are exactly 80 characters, the arming length, and only the last is short.

    A final chunk whose content reaches column 80 arms a continuation into the line after the record;
    an empty physical line appended after it absorbs that arming, adds nothing to the value and arms
    nothing further.  A chunk that does not arm must not get one, and arming is decided by content
    length (``_continues`` ignores trailing whitespace), so a space-padded 80-column chunk needs none.

    Two reader messages are expected on this path and are not damage: a chunk boundary falling just
    before a ``$`` inside a value makes the reader report absorbing a keyword, though the value comes
    back byte for byte.  The one case that does not round-trip is a non-final chunk 80 characters long
    ending in whitespace: it arms nothing, so the next line is a permissive continuation and the value
    reads back with a ``\n`` in it -- logged, wrong where it can be seen.
    """
    chunks = []
    rest = text
    while len(rest) > _CONTINUATION_LENGTH:
        chunks.append(rest[:_CONTINUATION_LENGTH])
        rest = rest[_CONTINUATION_LENGTH:]
    if len(rest.rstrip()) >= _CONTINUATION_LENGTH:
        # The last chunk arms continuation into the next line; the empty line absorbs that arming
        # without adding to the value.
        chunks.append(rest)
        chunks.append('')
    else:
        chunks.append(rest)
    for chunk in chunks[:-1]:
        if len(chunk.rstrip()) < _CONTINUATION_LENGTH:
            log.append(LogRecord('rdf:wrap-at-space', (),
                                 f'{what} wraps at a space, and a space-padded line continues nothing; the '
                                 f'value will read back with a newline in it', REPAIRED))
            break
    return chunks


def _log_stripped_whitespace(log, text, what, *, first):
    """Report the edge whitespace of one physical line the reader will eat -- and only that.

    ``parse_rdf_fields`` reads the first physical chunk as ``line[len(keyword):].strip()``, while every
    positional continuation after it is concatenated raw and a permissive continuation is appended raw
    as a new element.  So leading whitespace on a keyword-carrying line always goes, and trailing
    whitespace goes only when the whole text fits inside that first chunk --
    ``len(text) <= _CONTINUATION_LENGTH - _KEYWORD_WIDTH``, 73 characters.  A longer text carries its
    tail on an untouched continuation and round-trips exactly, so there is nothing to report;
    reporting it anyway would be a false ``unsupported: ``, which claims the format is the limitation.
    """
    if not first:
        # Appended raw as its own element; neither edge is touched.
        return
    if text[:1].isspace():
        log.append(LogRecord('rdf:leading-whitespace', (),
                             f'unsupported: {what} has leading whitespace that the format strips on read: '
                             f'{text!r:.40}', LOST))
    if text[-1:].isspace() and len(text) <= _CONTINUATION_LENGTH - _KEYWORD_WIDTH:
        log.append(LogRecord('rdf:trailing-whitespace', (),
                             f'unsupported: {what} has trailing whitespace that the format strips on read: '
                             f'{text!r:.40}', LOST))


class RDFWrite(_FileBacked):
    """MDL V2000 RDfile writer.  ``write(molecule)`` or ``write(reaction)``; path or buffer; ``with``.

    The record tag follows what it is handed -- ``$MFMT`` for a molecule, ``$RFMT`` for a reaction.
    Metadata is written as ``$DTYPE``/``$DATUM`` pairs and never as an SDF's ``>  <name>``.

    Use :class:`ERDFWrite` for V3000 CTABs.
    """
    __slots__ = ('_started',)
    _stamp = V2000_STAMP

    def __init__(self, file, *, append=False):
        """:param append: add to an existing RDfile, whose header is already in it."""
        if append and isinstance(file, (str, Path)):
            # A new or empty file needs a header even in append mode, and tell() is no use here:
            # _FileBacked accepts any object with a .write(), and a pipe has no tell().
            try:
                already_started = Path(file).stat().st_size > 0
            except FileNotFoundError:
                already_started = False
        else:
            already_started = append
        self._open(file, 'a' if append else 'w', 'write')
        self._started = already_started

    def write(self, data, *, meta=None, title=None):
        """Write one record -- a molecule as ``$MFMT``, a reaction as ``$RFMT``.  Returns the log.

        :param meta: ``{name: value}``, or any iterable of pairs -- the record's ``$DTYPE``/``$DATUM``
            pairs.  ``None`` means whatever the container holds, a molecule's now as well as a
            reaction's; pass ``meta={}`` to write a record with none.
        :param title: a replacement name line, or ``None`` for the structure's own.
        """
        log = []
        if meta is None:
            meta = data.meta
        if not self._started:
            self._file.write(f'{RDF_HEADER}\n{_DATM}    {strftime("%m/%d/%y %H:%M")}\n')
            self._started = True
        # Build the record lines BEFORE writing the tag.  A refusal from the emitter must leave the
        # file as it found it; a tag with no record body underneath it is not a valid empty RDfile.
        if isinstance(data, ReactionContainer):
            lines, log = emit_rxn(data, version=self._stamp, title=title, log=log)
            self._file.write(f'{_RXN_RECORD}\n')
        else:
            # `separator=False`: `$$$$` is SD framing.  An RDfile record is framed by the next tag and
            # its metadata is the $DTYPE/$DATUM tail below.
            # `meta={}`: an RDfile record's metadata is the $DTYPE/$DATUM tail written below, never an
            # SDF's `>  <name>` block, so the molfile emitter must write none of it.
            lines, log = emit_record(data, meta={}, version=self._stamp, title=title,
                                     separator=False, log=log)
            self._file.write(f'{_MOL_RECORD}\n')
        self._file.write('\n'.join(lines))
        self._file.write('\n')

        # `Mapping` and not `dict`: a caller may hand back any mapping, and an iterable of pairs is
        # accepted for one that built its fields in order.
        for name, value in (meta.items() if isinstance(meta, Mapping) else (meta or ())):
            value = str(value)
            # Which edge survives depends on where the wrap puts it, so the name and every value part
            # ask `_log_stripped_whitespace` rather than deciding here.
            _log_stripped_whitespace(log, name, f'{_DTYPE} name', first=True)
            self._file.write('\n'.join(_wrap_logical_line(f'{_DTYPE} {name}', log,
                                                          f'the $DTYPE name {name!r:.30}')) + '\n')
            # A newline inside a value is a new physical line: that is how the reader stores a
            # multi-line value, and each of those lines is then wrapped by the same rule.
            for i, part in enumerate(value.split('\n')):
                # Only the first physical line carries the keyword; the rest are permissive
                # continuations.
                what = f'the $DATUM value of {name!r:.30}'
                _log_stripped_whitespace(log, part, what, first=not i)
                first = f'{_DATUM} {part}' if not i else part
                if i and part.startswith('$'):
                    # The format has no escape for this, so say what will happen: the reader ends a
                    # value at a $-led line.
                    log.append(LogRecord('rdf:dollar-in-value', (),
                                         f'{what} has a line beginning with $, which the reader will take for '
                                         f'a keyword and not for value text: {part!r:.40}', LOST))
                self._file.write('\n'.join(_wrap_logical_line(first, log, what)) + '\n')
        return log


class ERDFWrite(RDFWrite):
    """MDL V3000 RDfile writer.  Same surface as :class:`RDFWrite`, extended CTABs.

    The RDfile framing is identical -- ``$RDFILE``, ``$MFMT``/``$RFMT``, ``$DTYPE``/``$DATUM`` are not
    versioned -- so only the stamp handed to the CTAB emitters changes.
    """
    __slots__ = ()
    _stamp = V3000_STAMP
