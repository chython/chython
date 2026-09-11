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
"""SDF and molfile framing: record separation, per-record version sniffing, data fields.

``$$$$`` is matched as a line prefix (writers pad it, and a data value may contain the characters)
and a missing final separator is normal.  The version stamp is sniffed per record, because one
``.sdf`` may mix versions, and a record body carrying ``M  V30`` outranks a wrong stamp.
"""

from collections.abc import Mapping

from ...core import LogRecord, LOST, REPAIRED
from ._errors import MalformedCtfile
from ._sgroup import UNSUPPORTED
from ._v2000 import V2000_STAMP, emit_v2000, parse_v2000
from ._v3000 import V3000_STAMP, emit_v3000, parse_v3000


__all__ = ['split_records', 'sniff_version', 'parse_data_fields', 'parse_record', 'emit_record',
           'RECORD_SEPARATOR', 'UNPARSED_KEY']


#: The record separator, matched as a line prefix.
RECORD_SEPARATOR = '$$$$'

# The spec puts the version stamp at columns 34-39 of the fourth line.  A stamp one column off is a
# common malformation, so the search covers the whole line and this span is only what a writer fills.
_STAMP_SPAN = (33, 39)


def split_records(lines):
    """Yield each record's lines, separator excluded.  `lines` is any iterable of strings.

    Line endings are stripped, including a lone ``\\r`` from a file written on one platform and
    edited on another.  A trailing chunk with no separator is yielded if it holds anything at all,
    and blank padding at the end of a file is not a record.
    """
    current = []
    for line in lines:
        line = line.rstrip('\n').rstrip('\r')
        if line.startswith(RECORD_SEPARATOR):
            yield current
            current = []
        else:
            current.append(line)
    if any(line.strip() for line in current):
        yield current


def sniff_version(lines, log=None):
    """``'V2000'`` or ``'V3000'`` for one record's lines.

    The body outranks the stamp: a ``M  V30`` line means V3000 whatever the counts line claims, and
    an absent stamp with no ``M  V30`` means V2000, since V3000 cannot express a CTAB without its own
    keyword lines.
    """
    log = [] if log is None else log
    stamped = ''
    if len(lines) > 3:
        line = lines[3]
        for candidate in (V3000_STAMP, V2000_STAMP):
            # The declared span first, then anywhere on the line, so a stamp written one column off
            # is still read rather than silently becoming a V2000 record.
            if line[_STAMP_SPAN[0]:_STAMP_SPAN[1]].strip() == candidate or candidate in line:
                stamped = candidate
                break

    body = ''
    for line in lines[4:]:
        if line.startswith('M  V30'):
            body = V3000_STAMP
            break
        if line.startswith('M  END'):
            break

    if body and stamped and body != stamped:
        log.append(LogRecord('sdf:stamp-body-mismatch', (),
                             f'counts line is stamped {stamped} but the record body is {body}; '
                             f'read as {body}', REPAIRED))
    version = body or stamped or V2000_STAMP
    if not stamped and version == V2000_STAMP:
        log.append(LogRecord('sdf:no-stamp', (),
                             'no version stamp on the counts line; read as V2000', REPAIRED))
    return version


#: Where a line the reader could not attribute to any field goes.  The key is stable, so a caller may
#: grep for it.  `emit_record` does NOT write it back -- see there.
UNPARSED_KEY = 'chython_unparsed_metadata'


def parse_data_fields(lines, log=None):
    """Parse the data-field block that follows ``M  END``.  Returns ``{name: value}``.

    A field is a ``> <NAME>`` header, value lines and a blank line; the name is taken from between the
    first ``<`` and the next ``>``.  `lines` may be the whole record; everything up to and including
    ``M  END`` is skipped.  Value lines are kept verbatim and joined with ``'\n'``.

    A REPEATED NAME MERGES, its value lines appended to the earlier ones.  One value per name is what a
    mapping can spell, so the merge is reported rather than left to be discovered.  The verbatim
    ``> ...`` header is dropped for the same reason: the spec also allows a field number and an external
    registry number there, neither is modelled, and a header carrying one is reported ``unsupported: ``.
    """
    log = [] if log is None else log
    fields = {}
    current = None
    started = False
    for line in lines:
        if not started:
            if line.startswith('M  END'):
                started = True
            continue
        if line.startswith(RECORD_SEPARATOR):
            # The separator ends the record, so it is not data.  It matters now that an unattributable
            # line is STORED rather than logged and dropped: a caller handing `parse_record` the lines
            # with the `$$$$` still on them would otherwise get it back as a metadata value.
            break
        if line.startswith('>'):
            name = _field_name(line)
            if name is None:
                log.append(LogRecord('sdf:no-field-name', (),
                                     f'data field header with no <name>: {line!r:.60}; '
                                     f'field kept under its header text'))
                name = line[1:].strip()
            else:
                start = line.find('<')
                rest = (line[1:start] + line[line.find('>', start) + 1:]).strip()
                if rest:
                    log.append(LogRecord('sdf:field-number-not-modelled', (),
                                         f'{UNSUPPORTED}data field {name!r:.30}: the header line carries '
                                         f'{rest!r:.30} beside the name -- a field or registry number is not '
                                         f'modelled and is not written back', LOST))
            if name in fields:
                log.append(LogRecord('sdf:duplicate-field', (),
                                     f'data field {name!r:.30} appears twice and one value per name is what '
                                     f'`meta` holds; the values are merged'))
            current = name
            fields.setdefault(name, [])
        elif not line.strip():
            # A blank line closes the value.  It does not end the block: the next `>` opens the next
            # field, and files with two blank lines between fields are ordinary.
            current = None
        elif current is not None:
            fields[current].append(line)
        else:
            log.append(LogRecord('sdf:unparsed-data', (),
                                 f'data outside any field, stored under {UNPARSED_KEY!r}: {line!r:.60}'))
            fields.setdefault(UNPARSED_KEY, []).append(line)
    return {k: '\n'.join(v) for k, v in fields.items()}


def _field_name(line):
    """The name between the first ``<`` and the next ``>``, or ``None``."""
    start = line.find('<')
    if start < 0:
        return None
    end = line.find('>', start + 1)
    if end < 0:
        return None
    return line[start + 1:end]


def parse_record(lines, log=None, *, ignore_stereo=False, header=None):
    """Parse one record's lines into a `MoleculeContainer`.

    The version is sniffed rather than assumed, and the data fields are parsed whether or not the
    record came from an SDF -- a molfile simply has none.  The fields land on ``mol.meta``, the name
    line on ``mol.title``, and everything the reader recovered on ``mol.log`` as well as on `log`.

    `header`, when given a dict, is filled with what the FILE said and no container holds:
    ``version``, ``program``, ``comment``, ``sgroups`` (the :class:`~._sgroup.SGroupStore`) and
    ``unknown_hydrogens``.  The name line is not among them: it is ``mol.title``, the same ``str``.
    An out-parameter, so the ordinary call still returns one object -- the shape
    :func:`~._rdf.split_rdf_records` already uses.
    """
    log = [] if log is None else log
    # THIS RECORD'S OWN LIST, not the caller's, and it is what the version sniffer, the CTAB parser and
    # the data-field parser all write to.  `parse_v2000`/`parse_v3000` alias it onto `ctab.log`, so
    # `build` returns it plus the build's own lines -- one list, in file order, holding exactly this
    # record.  A caller reusing one `log=` across records therefore gets record 3's lines appended
    # rather than record 1's folded onto record 3's molecule.
    own = []
    version = sniff_version(lines, own)
    ctab = parse_v3000(lines, own) if version == V3000_STAMP else parse_v2000(lines, own)
    ctab.meta = parse_data_fields(lines, own)
    mol, store, build_log = ctab.build(ignore_stereo=ignore_stereo)
    log.extend(build_log)
    if header is not None:
        header['version'] = version
        header['program'] = ctab.program
        header['comment'] = ctab.comment
        header['sgroups'] = store
        header['unknown_hydrogens'] = ctab.unknown_hydrogens
    return mol


def emit_record(mol, sgroups=None, meta=None, *, version=V2000_STAMP, title=None, program='',
                comment='', separator=True, log=None):
    """Render one record: molfile lines, then data fields, then ``$$$$``.

    `version` chooses the emitter; V2000 is the default and refuses rather than truncating when a
    structure does not fit its fixed columns, naming V3000 as the fix.  `title`, `sgroups` and `meta`
    default to ``None``, meaning "whatever the molecule holds" (see :func:`~._sgroup.resolve_output`);
    passing any of them explicitly overrides it, and ``meta={}`` writes no data fields at all.
    """
    log = [] if log is None else log
    if version == V3000_STAMP:
        lines, log = emit_v3000(mol, sgroups, title=title, program=program, comment=comment,
                                log=log)
    elif version == V2000_STAMP:
        lines, log = emit_v2000(mol, sgroups, title=title, program=program, comment=comment,
                                log=log)
    else:
        raise MalformedCtfile(f'unknown CTfile version {version!r}; expected {V2000_STAMP} or '
                              f'{V3000_STAMP}')
    if meta is None:
        meta = mol.meta
    for name, value in (meta.items() if isinstance(meta, Mapping) else meta):
        if name == UNPARSED_KEY:
            log.append(LogRecord('sdf:unparsed-not-written', (),
                                 f'the {UNPARSED_KEY} field is not written back: it holds lines the reader could '
                                 f'not attribute to any field, not a field of its own', LOST))
            continue
        lines.append(f'>  <{name}>')
        lines.extend(str(value).split('\n'))
        lines.append('')  # the blank line is what closes a value; without it the next field is data
    if separator:
        lines.append(RECORD_SEPARATOR)
    return lines, log
