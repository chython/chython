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
"""The V3000 physical-line layer: continuation joining, tokenizing, and their inverses.

Join physical lines first -- a line ending in ``-`` continues, and the next line's ``M  V30 ``
prefix is stripped -- then tokenize the logical line into positional values and ``KEY=value`` pairs.
The other order breaks on a quoted string split across a continuation, which real files contain.
"""
from ...core import LogRecord


__all__ = ['join_continuations', 'tokenize', 'parse_list', 'quote_value', 'emit_v30', 'V30_PREFIX']


V30_PREFIX = 'M  V30 '
_PREFIX_LEN = len(V30_PREFIX)
# CTfile: a V3000 physical line is at most 80 characters including the prefix.
_LINE_LIMIT = 80


def join_continuations(lines, log=None):
    """Join V3000 physical lines into logical lines.

    `lines` is an iterable of strings, each expected to carry the ``M  V30 `` prefix; returns the
    logical lines with the prefix removed.

    The concatenation is byte-exact -- a continued part keeps its trailing whitespace, because a
    quoted value may be split mid-string and the spaces around the break are inside the quotes
    (Pipeline Pilot's ``FIELDDISP`` strings depend on it).  A prefixless line is passed through
    intact and reported, never sliced at column 7, which would eat seven characters of data.
    """
    out = []
    parts = []
    for line in lines:
        line = line.rstrip('\r\n')
        if line.startswith(V30_PREFIX):
            body = line[_PREFIX_LEN:]
        elif line.startswith('M  V30'):
            # prefix present but short one space -- some writers emit `M  V30` with a single
            # trailing space, or none at all on an empty continuation.
            body = line[6:]
            if body.startswith(' '):
                body = body[1:]
        else:
            if log is not None:
                log.append(LogRecord('ctab:missing-v30-prefix', (),
                                     f'V3000 line without M  V30 prefix: {line!r:.60}'))
            body = line
        if body.endswith('-'):
            parts.append(body[:-1])
        elif parts:
            parts.append(body)
            out.append(''.join(parts))
            parts = []
        else:
            out.append(body)
    if parts:
        # A trailing `-` with nothing after it.  The spec has no such thing; take what we have.
        if log is not None:
            log.append(LogRecord('ctab:dangling-continuation', (), 'dangling continuation'))
        out.append(''.join(parts))
    return out


def tokenize(line, log=None):
    """Split one logical V3000 line into a list of tokens.

    A token is either a positional value or a ``KEY=value`` pair, returned as the raw string
    including the ``KEY=``.  Values may be:

    * double-quoted, with ``""`` for a literal quote;
    * a parenthesised list ``(N v1 ... vN)``, which may itself contain quoted strings;
    * a bare run of non-space characters;
    * a bare run *containing* spaces -- which the spec forbids and writers nonetheless produce.
      Such a value is consumed up to the next ``KEY=`` or end of line and reported.

    A scanner rather than a regex because of the last case: deciding where an unquoted value ends
    requires looking ahead for a following ``KEY=``, which an alternation over the whole line cannot
    do -- it splits the value into several positional tokens instead.
    """
    tokens = []
    n = len(line)
    i = 0
    while i < n:
        if line[i] == ' ':
            i += 1
            continue
        start = i
        # scan a key
        key_end = -1
        while i < n and line[i] not in ' ="(':
            i += 1
        if i < n and line[i] == '=':
            key_end = i
            i += 1
        if key_end < 0:
            # positional token: it may still be quoted or parenthesised
            if i < n and line[i] == '"':
                i = _skip_quoted(line, i, log)
            elif i < n and line[i] == '(':
                i = _skip_parens(line, i, log)
            else:
                while i < n and line[i] != ' ':
                    i += 1
            tokens.append(line[start:i])
            continue
        # KEY= ... : the value
        if i < n and line[i] == '"':
            i = _skip_quoted(line, i, log)
        elif i < n and line[i] == '(':
            i = _skip_parens(line, i, log)
        else:
            j = i
            while j < n and line[j] != ' ':
                j += 1
            # An unquoted value with a space in it (a spec violation, reported once): the value ends
            # at the next `KEY=` on this line, or at end of line.  In every V3000 block line the
            # positional values come *before* the keywords -- ATOM is `index type x y z aamap
            # KEY=...`, SGROUP is `index type parent KEY=...` -- so a bare token after a `KEY=` pair
            # is the rest of an unquoted value, not a positional one a consumer would read as an
            # atom index.  `FIELDDATA=Molecular Weight: 375,40` is such a file.
            if j < n:
                k = _next_key(line, j)
                # `len(line.rstrip())` and not `n`: an unquoted value has no way to express
                # significant edge whitespace, so trailing spaces on the line are layout, not data.
                end = k if k >= 0 else len(line.rstrip())
                # `end == j` is the conforming line: the bare run ended and the next `KEY=` -- or
                # the line's trailing layout -- begins at that very space, so no space was taken
                # into the value and there is nothing to report.  Reporting on `j < n` instead
                # names every non-final key on every keyworded line.
                if end > j and log is not None:
                    log.append(LogRecord('ctab:unquoted-spaces', (),
                                         f'unquoted value with spaces for {line[start:key_end]}'))
                j = end
            i = j
        tokens.append(line[start:i])
    return tokens


def _skip_quoted(line, i, log):
    """Return the index just past the double-quoted run starting at `line[i] == '"'`."""
    n = len(line)
    i += 1
    while i < n:
        if line[i] == '"':
            if i + 1 < n and line[i + 1] == '"':  # `""` is a literal quote
                i += 2
                continue
            return i + 1
        i += 1
    if log is not None:
        log.append(LogRecord('ctab:unterminated-quote', (), 'unterminated quoted value'))
    return n


def _skip_parens(line, i, log):
    """Return the index just past the parenthesised run starting at `line[i] == '('`.

    Quoted strings inside the parentheses are skipped whole, so a `)` inside quotes does not
    close the list.  Nesting is not a CTfile feature and is not supported.
    """
    n = len(line)
    i += 1
    while i < n:
        if line[i] == '"':
            i = _skip_quoted(line, i, log)
            continue
        if line[i] == ')':
            return i + 1
        i += 1
    if log is not None:
        log.append(LogRecord('ctab:unterminated-list', (), 'unterminated parenthesised list'))
    return n


def _next_key(line, i):
    """Index of the start of the next `KEY=` token at or after `i`, or ``-1`` when there is none.

    Used only to bound an unquoted value that contains spaces.  A `KEY=` must be preceded by a
    space and must consist of upper-case letters, digits and underscores -- the CTfile keyword
    alphabet -- so a value like `Molecular Weight: 375,40` cannot be mistaken for one.  The
    sentinel must not be a valid position: the caller distinguishes "stop before the key" from
    "take the rest of the line".
    """
    n = len(line)
    j = i
    while j < n:
        if line[j] != ' ':
            j += 1
            continue
        k = j + 1
        while k < n and (line[k].isupper() or line[k].isdigit() or line[k] == '_'):
            k += 1
        if k > j + 1 and k < n and line[k] == '=':
            return j
        j = k if k > j else j + 1
    return -1


def parse_list(value, log=None):
    """Parse a ``(N v1 ... vN)`` list value into a list of strings.

    `value` may be the raw ``KEY=(...)`` token or just the ``(...)`` part.  The declared count is
    checked but **not** trusted: when it disagrees with how many values are present the values win
    and the disagreement is reported, because a wrong count is a writer bug while the values are
    the data.
    """
    if '=' in value and not value.startswith('('):
        value = value.split('=', 1)[1]
    value = value.strip()
    if value.startswith('('):
        value = value[1:]
    if value.endswith(')'):
        value = value[:-1]
    items = value.split()
    if not items:
        return []
    try:
        declared = int(items[0])
    except ValueError:
        if log is not None:
            log.append(LogRecord('ctab:list-no-count', (),
                                 f'list without a leading count: {value!r:.40}'))
        return items
    rest = items[1:]
    if declared != len(rest) and log is not None:
        log.append(LogRecord('ctab:list-count-mismatch', (),
                             f'list count {declared} disagrees with {len(rest)} values'))
    return rest


def quote_value(value):
    """Quote a V3000 value if the grammar requires it.

    Quoting is required when the value is empty, holds a space, a parenthesis or a quote, or **ends
    with a hyphen**; a literal quote is doubled.  A value needing no quotes is emitted bare, as
    reference writers do.  One exception: a well-formed ``(...)`` list is emitted bare despite its
    spaces and parentheses, since quoting would turn it into a string and a reader looking for
    ``CSTATE=(4 ...)`` wants a list.

    The trailing hyphen is the continuation marker, so a line-final data hyphen is read as one: bare
    ``LABEL=NH3+Cl-`` ends its physical line and the reader joins ``END SGROUP`` onto the label.
    Quoting puts a ``"`` last instead.  `emit_v30` wrapping is safe without this -- it appends its
    own marker, giving ``--``, and the reader strips one -- but the last physical line has no marker
    to hide behind.
    """
    if value == '':
        return '""'
    if value.startswith('(') and value.endswith(')') and '"' not in value:
        return value
    if value.endswith('-') or any(c in value for c in ' ()"'):
        return '"' + value.replace('"', '""') + '"'
    return value


def emit_v30(content):
    """Render one logical line as the physical ``M  V30 `` lines it needs, as a list of strings.

    Wraps at 80 columns with a trailing ``-`` continuation marker.  Unlike the reader, the writer
    never breaks inside a quoted string: it backtracks to the last break point outside quotes.  A
    single token longer than the available width -- a 4000-character ``FIELDDATA`` -- cannot be
    wrapped, so there the break falls where it must; the reader rejoins byte-exactly either way.
    """
    room = _LINE_LIMIT - _PREFIX_LEN
    if len(content) <= room:
        return [V30_PREFIX + content]
    out = []
    pos = 0
    n = len(content)
    while n - pos > room:
        cut = _break_at(content, pos, pos + room - 1)
        out.append(V30_PREFIX + content[pos:cut] + '-')
        pos = cut
    out.append(V30_PREFIX + content[pos:])
    return out


def _break_at(content, start, limit):
    """Choose a break position in `content[start:]` at or before `limit`.

    Three preferences, in order: a token boundary (the character after a space outside quotes), any
    position outside quotes, then `limit`.  The first keeps ``KEY=value`` pairs intact, the second
    keeps quoted strings intact, the third is the unwrappable-single-token case.
    """
    last_token = -1
    last_safe = -1
    in_quote = False
    i = start
    while i <= limit and i < len(content):
        c = content[i]
        if c == '"':
            if in_quote and i + 1 < len(content) and content[i + 1] == '"':
                i += 2
                continue
            in_quote = not in_quote
        elif not in_quote:
            last_safe = i
            if c == ' ':
                last_token = i + 1
        i += 1
    if start < last_token <= limit:
        return last_token
    if last_safe > start:
        return last_safe
    return limit
