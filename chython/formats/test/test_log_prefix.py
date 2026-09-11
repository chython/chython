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
"""The one log-line convention this tree enforces: a line reporting an unmodelled construct carries
the `unsupported: ` prefix, so `any(str(x).startswith('unsupported') for x in log)` answers "did this file
state something we do not model?".  Every other line is prose for a human.  Deliberately one-way --
this must never be widened to require a prefix on every line.
"""

from ast import Call, Constant, JoinedStr, parse, unparse, walk
from pathlib import Path

from chython.formats.ctfile._sgroup import UNSUPPORTED


#: Words that mean "we did not model this", as opposed to "the file was broken and we coped".
#: See _WHITELIST for the ones that are repairs.
_MARKERS = ('is not modelled', 're-emitted', 'is a query field', 'ignored')

#: Log-emitting source lines that contain a marker word but report a repair, not a gap in chython.
#: Each entry is `(file, substring of the message)`.  Adding one claims the file was malformed and we
#: recovered; adding wrongly makes `unsupported` under-report, so say why.
_WHITELIST = (
    # A wedge on an atom that cannot carry one is the drawing being wrong.
    ('wedge.py', 'wedge drawn on a non-stereogenic centre'),
    # A bond stereo value outside the legal set is a malformed field.
    ('_v2000.py', 'is not 0, 1, 3, 4 or 6'),
    # Truncated or unreadable lines: malformed input.
    ('_v3000.py', 'atom line with unreadable index'),
    ('_v3000.py', 'atom index'),
    ('_v3000.py', 'sgroup line too short'),
    ('_v3000.py', 'unrecognised V3000 line'),
    ('_v3000.py', 'positional field'),
    # A bond CFG value outside the legal V3000 wedge-code set is a malformed field.
    ('_v3000.py', 'not a V3000 wedge code, ignored'),
    ('_sdf.py', 'data outside any field'),
    # An IMPL_H S-group naming several atoms, or carrying no datum, contradicts its own extension.
    ('_hydrogens.py', 'expected one of each, ignored'),
    # A stated valence below the bonds already drawn cannot be a total valence: a malformed field.
    ('_hydrogens.py', 'cannot be a total valence'),
    # A line appearing before the first $MFMT/$RFMT record tag is malformed RDfile structure.
    ('_rdf.py', 'line before the first record'),
    # A $DATUM line with no preceding $DTYPE: the file is broken, not a feature gap.
    ('_rdf.py', '{_DATUM} with no preceding {_DTYPE}'),
    # An unrecognised $-led line between data fields: broken file.  A construct the format really has
    # and chython declines to model -- $MIREG/$MEREG/$RIREG/$REREG -- gets the prefixed message.
    ('_rdf.py', 'unrecognised field keyword, ignored'),
    # A $DATM in a $DTYPE/$DATUM tail is the header timestamp where the format has no room for it.
    # The header's own $DATM is stored on `RDFRead.date` and logs nothing.
    ('_rdf.py', 'outside the file header, ignored'),
    # A plain line with no open $DATUM value to join it to: malformed RDfile structure, not a gap.
    ('_rdf.py', 'line outside any $DATUM value'),
    # CML's <atomParity> carries a signed number; a non-numeric value is a malformed field.
    ('_cml.py', 'atomParity value'),
    # CML's schema fixes atomRefs4 at four, so three or five is the document contradicting itself.
    ('_cml.py', 'atoms, not four'),
    # A reference to an atom id that does not exist in the document is malformed input.
    ('_cml.py', 'atomParity references unknown atom'),
    # The parity atom may stand in for its own implicit hydrogen once; naming it twice is no geometry.
    ('_cml.py', 'phantom direction twice'),
    # A file claiming MDL's dictionary and writing a term outside it is malformed -- the MDL wedge
    # vocabulary itself is read.  A `convention` naming some *other* dictionary carries the prefix.
    ('_cml.py', 'bondStereo conventionValue'),
    # <bondStereo/> empty with no convention states nothing; W, H, C, T and convention="MDL" are read.
    ('_cml.py', 'empty bondStereo'),
)


#: Modules outside `formats/` that write into a writer's log and so are in this ratchet's scope.
#: `core/wedge.py` emits the coordinate-free cis/trans loss for four writers, each handing it the
#: caller's log; scope follows who writes the line, not which directory it sits in.
_EXTERNAL = ('core/wedge.py',)


def _sources():
    root = Path(__file__).resolve().parent.parent          # chython/formats/
    files = [p for p in root.rglob('*.py') if 'test' not in p.parts]
    for name in _EXTERNAL:
        path = root.parent / name                          # chython/
        assert path.is_file(), f'{name} is named in _EXTERNAL and does not exist'
        files.append(path)
    return sorted(files)


def _whitelisted(path, text):
    return any(path.name == name and probe in text for name, probe in _WHITELIST)


def _messages(source):
    """Every message a `LogRecord` in this source is built with, as `(is_fstring, text)`.

    AST and not a regex, because the message is no longer the first literal in the call -- the rule id
    is -- and a multi-line `LogRecord(...)` is ordinary.  Still not a second implementation of the thing
    being checked: this reads argument positions, never meaning.
    """
    for node in walk(parse(source)):
        if not isinstance(node, Call) or getattr(node.func, 'id', None) != 'LogRecord':
            continue
        arg = next((k.value for k in node.keywords if k.arg == 'message'), None)
        if arg is None and len(node.args) > 2:
            arg = node.args[2]
        if isinstance(arg, Constant) and isinstance(arg.value, str):
            yield False, arg.value
        elif isinstance(arg, JoinedStr):
            # Rebuilt with each placeholder as its own source text, so `f'{UNSUPPORTED}...'` still tests
            # as prefixed and a whitelist probe like `'{_DATUM} with no preceding {_DTYPE}'` still matches.
            yield True, ''.join(v.value if isinstance(v, Constant) else '{' + unparse(v.value) + '}'
                                for v in arg.values)


def test_unmodelled_constructs_are_prefixed():
    """Every log line reporting something chython does not model starts with `unsupported: `."""
    offenders = []
    for path in _sources():
        source = path.read_text(encoding='utf8')
        for is_fstring, message in _messages(source):
            if not any(m in message for m in _MARKERS):
                continue
            if _whitelisted(path, message):
                continue
            # Accept the literal prefix, or `f'{UNSUPPORTED}...'` as the two sgroup sites in
            # _v3000.py write it.  The `f` is required: without it the braces are literal characters.
            if not (message.startswith(UNSUPPORTED) or (is_fstring and message.startswith('{UNSUPPORTED}'))):
                offenders.append(f'{path.name}: {message[:70]}')
    assert not offenders, 'unmodelled-construct log lines without the prefix:\n' + '\n'.join(offenders)


def test_the_question_is_answerable_on_a_real_file():
    """The bar from the spec: the prefix answers a real question on a real record."""
    from chython.formats.ctfile import parse_record

    # A V3000 record with a LINKNODE line: nothing in chython models one.
    lines = ['linknode', '', '',
             '  0  0  0     0  0            999 V3000',
             'M  V30 BEGIN CTAB',
             'M  V30 COUNTS 2 1 0 0 0',
             'M  V30 BEGIN ATOM',
             'M  V30 1 C 0 0 0 0',
             'M  V30 2 C 1 0 0 0',
             'M  V30 END ATOM',
             'M  V30 BEGIN BOND',
             'M  V30 1 1 1 2',
             'M  V30 END BOND',
             'M  V30 LINKNODE 1 2 1 1 2 1 3',
             'M  V30 END CTAB',
             'M  END']
    log = []
    parse_record(lines, log)
    assert any(str(x).startswith(UNSUPPORTED) for x in log), log


def test_sgroup_marker_survives_merge():
    """The prefix is still at position 0 after `merge_log` prepends the sgroup location.

    The source-level test above cannot see runtime composition: `'sgroup 1 SUP: unsupported: ...'`
    would pass it and break every caller doing `startswith('unsupported: ')`.
    """
    from chython.formats.ctfile import parse_v3000

    # A V3000 record containing a SUP S-group with a SAP= keyword.  SAP holds atom indices and is
    # in _INDEX_VALUED -- it is dropped and logged as an unsupported construct.
    lines = ['sap_test', '', '',
             '  0  0  0     0  0            999 V3000',
             'M  V30 BEGIN CTAB',
             'M  V30 COUNTS 2 1 1 0 0',
             'M  V30 BEGIN ATOM',
             'M  V30 1 C 0 0 0 0',
             'M  V30 2 C 1 0 0 0',
             'M  V30 END ATOM',
             'M  V30 BEGIN BOND',
             'M  V30 1 1 1 2',
             'M  V30 END BOND',
             'M  V30 BEGIN SGROUP',
             'M  V30 1 SUP 0 ATOMS=(2 1 2) SAP=(3 1 2 1)',
             'M  V30 END SGROUP',
             'M  V30 END CTAB',
             'M  END']
    log = []
    parse_v3000(lines, log)
    sap_lines = [x for x in log if 'SAP' in x]
    assert sap_lines, f'expected a SAP log line, got: {log}'
    for line in sap_lines:
        assert str(line).startswith(UNSUPPORTED), (
            f'marker buried mid-string -- startswith check would miss it: {line!r}')
