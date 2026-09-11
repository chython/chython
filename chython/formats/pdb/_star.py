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
"""The STAR/CIF grammar and no vocabulary: ``data_`` blocks, ``loop_`` tables, quoted values,
``;``-delimited text fields and the two nulls.  CIF spells *inapplicable* ``.`` and *unknown* ``?``,
kept apart as the :data:`INAPPLICABLE` and :data:`UNKNOWN` singletons (:func:`is_null` for either).
:func:`parse_star` streams: it consumes lines and yields each block as the next ``data_`` closes it.
Nothing raises -- every grammar-level malformation is stored as best it can be and logged.
"""
from collections.abc import Iterable, Iterator
from pathlib import Path

from ...core import LogRecord, LOST, REPAIRED


__all__ = ['INAPPLICABLE', 'UNKNOWN', 'StarBlock', 'StarLoop', 'is_null', 'parse_star']


class _Null:
    """A CIF null.  Two instances exist and they are not equal to each other or to anything else."""
    __slots__ = ('_spelling',)

    def __init__(self, spelling: str):
        self._spelling = spelling

    def __repr__(self):
        return f'<CIF {self._spelling}>'

    def __str__(self):
        return self._spelling

    def __bool__(self):
        # Falsy, so `value or default` reads naturally at a call site that does not care which null
        # it got.  A caller that does care compares identity.
        return False


#: ``.`` -- the item does not apply to this row.  A polymer sequence number on a water oxygen.
INAPPLICABLE = _Null('.')
#: ``?`` -- the item applies and its value is not known.
UNKNOWN = _Null('?')


def is_null(value) -> bool:
    """True for either CIF null.  False for ``''``, ``0`` and ``0.0``, which are stated values."""
    return value is INAPPLICABLE or value is UNKNOWN


class StarLoop:
    """One ``loop_``: its header tags, lowercased, and its rows.

    Rows are lists of values in header order.  A row is exactly as long as ``tags`` -- a short final
    row is padded with :data:`UNKNOWN` and logged rather than dropped, because the atoms in the rows
    before it are not at fault.
    """
    __slots__ = ('tags', 'rows', '_index')

    def __init__(self, tags: list[str], rows: list[list]):
        self.tags = tags
        self.rows = rows
        # First occurrence wins: a duplicated header tag is logged by the parser, and a lookup has
        # to resolve to something.
        self._index = {}
        for i, tag in enumerate(tags):
            self._index.setdefault(tag, i)

    def has(self, tag: str) -> bool:
        return tag.lower() in self._index

    def column(self, tag: str) -> int:
        """Column position of *tag*, or -1 when the loop does not carry it."""
        return self._index.get(tag.lower(), -1)

    def value(self, row: list, tag: str, default=None):
        """The value of *tag* in *row*, or *default* when the loop has no such column."""
        i = self._index.get(tag.lower(), -1)
        return default if i < 0 else row[i]

    def __len__(self):
        return len(self.rows)

    def __repr__(self):
        return f'StarLoop({len(self.tags)} tags, {len(self.rows)} rows)'


class StarBlock:
    """One ``data_`` block: its scalar items and its loops.

    ``items`` maps a lowercased tag to its value.  ``loops`` is every ``loop_`` in the block, in file
    order.  :meth:`loop` finds the one carrying a given tag, which is how a vocabulary layer asks for
    ``_atom_site`` without caring where in the file it sat.
    """
    __slots__ = ('name', 'items', 'loops')

    def __init__(self, name: str | None):
        self.name = name
        self.items: dict = {}
        self.loops: list[StarLoop] = []

    def loop(self, tag: str) -> StarLoop | None:
        """The first loop carrying *tag*, or ``None``."""
        tag = tag.lower()
        for loop in self.loops:
            if loop.has(tag):
                return loop
        return None

    def get(self, tag: str, default=None):
        """The scalar item *tag*, or *default*.  Does not look inside loops."""
        return self.items.get(tag.lower(), default)

    def categories(self) -> set:
        """Every ``_category`` name present, from scalar items and from loop headers alike.

        A vocabulary layer uses this to name the categories it did not read.
        """
        found = set()
        for tag in self.items:
            found.add(tag.split('.', 1)[0])
        for loop in self.loops:
            for tag in loop.tags:
                found.add(tag.split('.', 1)[0])
        return found

    def __repr__(self):
        return f'StarBlock({self.name!r}, {len(self.items)} items, {len(self.loops)} loops)'


# --------------------------------------------------------------------------- line sources

def _iter_lines(source: str | Path | Iterable[str]) -> Iterator[str]:
    """Lines from *source*.

    A ``str`` is always the file's text and never a path; a :class:`~pathlib.Path` is always a path,
    read line by line so a 25 MB entry is never held whole.  Anything else is an iterable of lines.
    """
    if isinstance(source, Path):
        with source.open(encoding='utf8', errors='replace') as f:
            yield from f
    elif isinstance(source, str):
        yield from source.splitlines()
    else:
        yield from source


# --------------------------------------------------------------------------- lexer

#: Token kinds.  `bare` is the only kind a keyword or a tag can arrive as: a value written `'loop_'`
#: or in a `;` field is data, and treating it as syntax is the classic CIF reader bug.
_BARE = 'bare'
_QUOTED = 'quoted'
_TEXT = 'text'


def _tokens(lines: Iterable[str], log: list) -> Iterator[tuple]:
    """Yield ``(kind, value, lineno)`` for every token in *lines*.

    One left-to-right pass per line, no backtracking and no recursion.
    """
    # Manual counter, not `enumerate`: a `;` field consumes lines from the same iterator, so an
    # `enumerate` counter would be wrong on every line after the first text field.
    source = iter(lines)
    lineno = 0
    for raw in source:
        lineno += 1
        # CRLF and a bare CR both appear in files written on other platforms.
        line = raw.rstrip('\n').rstrip('\r')

        if lineno == 1 and line.startswith('#\\#CIF_2.0'):
            # CIF 2.0 announces itself with this exact magic and adds triple-quoted strings, lists and
            # tables.  Read on with the CIF 1.1 grammar, right for every construct the two share.
            log.append(LogRecord('star:cif2-grammar', (),
                                 'unsupported: CIF 2.0 syntax declared by the leading magic comment; '
                                 'the file is read with the CIF 1.1 grammar and any CIF 2.0-only '
                                 'value will not parse',
                                 LOST))
            continue

        i = 0
        if line[:1] == ';':
            # A `;` opens a text field only in column 1.  Collect until a line that also opens in
            # column 1 with `;`.
            open_text_line = lineno
            collected = [line[1:]]
            terminator = None
            for raw2 in source:
                lineno += 1
                line2 = raw2.rstrip('\n').rstrip('\r')
                if line2[:1] == ';':
                    terminator = line2
                    break
                collected.append(line2)
            yield _TEXT, '\n'.join(collected), open_text_line
            if terminator is None:
                log.append(LogRecord('star:unterminated-text-field', (),
                                     f'record: multi-line text field opened on line {open_text_line} '
                                     f'is never closed; its value runs to the end of the file'))
                return
            # Whatever follows the closing `;` on its line is lexed normally: the format allows only
            # whitespace there, but a writer that put a value there still gets read.
            line = terminator
            i = 1

        n = len(line)
        while i < n:
            ch = line[i]
            if ch in ' \t':
                i += 1
                continue
            if ch == '#':
                break  # comment to end of line
            if ch == "'" or ch == '"':
                # The quote closes only when the next character is whitespace or end of line, so
                # `O5'` and `can't` inside a quoted value need no escaping and get none.
                j = i + 1
                while j < n:
                    if line[j] == ch and (j + 1 == n or line[j + 1] in ' \t'):
                        break
                    j += 1
                if j >= n:
                    log.append(LogRecord('star:unterminated-quote', (),
                                         f'record: quoted value opened on line {lineno} is not '
                                         f'closed; the rest of the line is taken as the value',
                                         REPAIRED))
                    yield _QUOTED, line[i + 1:], lineno
                    break
                yield _QUOTED, line[i + 1:j], lineno
                i = j + 1
                continue
            j = i
            while j < n and line[j] not in ' \t':
                j += 1
            yield _BARE, line[i:j], lineno
            i = j


# --------------------------------------------------------------------------- parser

def _value(kind: str, token: str):
    """A token as a value: the two nulls only when written bare and alone."""
    if kind is _BARE:
        if token == '.':
            return INAPPLICABLE
        if token == '?':
            return UNKNOWN
    return token


def parse_star(source: str | Path | Iterable[str], *, log: list | None = None) \
        -> Iterator[StarBlock]:
    """Yield every ``data_`` block in *source*.

    *log* is the caller's list; grammar-level damage is appended to it.  Nothing here raises on a
    malformed file: a block always comes back, holding as much as the text supported.
    """
    log = [] if log is None else log
    lines = _iter_lines(source)

    block: StarBlock | None = None
    pending_tag: str | None = None
    pending_line = 0
    loop_tags: list[str] | None = None
    loop_values: list | None = None
    in_loop_header = False
    save_frame: str | None = None
    # In a file that is not CIF every word is an orphan value, so they are counted into one line.
    orphans = 0
    orphan_first: str | None = None
    saw_data_block = False

    def close_loop():
        """Turn the flat value list into rows and attach it to the block."""
        nonlocal loop_tags, loop_values, in_loop_header
        if loop_tags is None:
            return
        tags, values = loop_tags, loop_values
        loop_tags = loop_values = None
        in_loop_header = False
        if not tags:
            if values:
                log.append(LogRecord('star:loop-no-headers', (),
                                     f'record: loop_ with no header tags carries {len(values)} '
                                     f'value(s); they have no name and are not stored',
                                     LOST))
            return
        width = len(tags)
        rows = [values[i:i + width] for i in range(0, len(values), width)]
        if rows and len(rows[-1]) < width:
            short = rows[-1]
            log.append(LogRecord('star:short-loop-row', (),
                                 f'record: loop_ over {tags[0].split(".", 1)[0]} states {width} '
                                 f'column(s); its last row holds {len(short)}, padded to width with '
                                 f'unknown',
                                 REPAIRED))
            short.extend([UNKNOWN] * (width - len(short)))
        elif not rows:
            log.append(LogRecord('star:empty-loop', (),
                                 f'record: loop_ over {tags[0].split(".", 1)[0]} states {width} '
                                 f'column(s) and holds no rows'))
        if block is not None:
            block.loops.append(StarLoop(tags, rows))

    def ensure_block():
        """A block to put values in, for a file whose first item precedes its ``data_``."""
        nonlocal block
        if block is None:
            log.append(LogRecord('star:pre-block-item', (),
                                 'record: an item appears before any data_ block; it is stored in '
                                 'an unnamed block',
                                 REPAIRED))
            block = StarBlock(None)
        return block

    for kind, token, lineno in _tokens(lines, log):
        lowered = token.lower() if kind is _BARE else ''

        # ---- save frames.  Dictionaries use them; a data file does not.  Skip the frame's
        # contents rather than mixing dictionary definitions into the block's items.
        if save_frame is not None:
            if lowered == 'save_':
                save_frame = None
            continue
        if kind is _BARE and lowered.startswith('save_') and len(token) > 5:
            close_loop()
            pending_tag = None
            save_frame = token[5:]
            log.append(LogRecord('star:save-frame', (),
                                 f'unsupported: STAR save frame {save_frame!r} on line {lineno} '
                                 f'is skipped; nothing in this grammar layer models a frame',
                                 LOST))
            continue

        if kind is _BARE and lowered.startswith('data_'):
            close_loop()
            if pending_tag is not None:
                log.append(LogRecord('star:item-no-value', (),
                                     f'record: item {pending_tag} on line {pending_line} has no '
                                     f'value',
                                     LOST))
                pending_tag = None
            if block is not None:
                yield block
            block = StarBlock(token[5:])
            saw_data_block = True
            continue

        if kind is _BARE and lowered == 'global_':
            close_loop()
            pending_tag = None
            log.append(LogRecord('star:global-block', (),
                                 f'unsupported: STAR global_ block on line {lineno} is skipped; '
                                 f'its values apply to every block and nothing here models that',
                                 LOST))
            continue

        if kind is _BARE and lowered == 'loop_':
            close_loop()
            if pending_tag is not None:
                log.append(LogRecord('star:item-no-value', (),
                                     f'record: item {pending_tag} on line {pending_line} has no '
                                     f'value',
                                     LOST))
                pending_tag = None
            ensure_block()
            loop_tags = []
            loop_values = []
            in_loop_header = True
            continue

        if kind is _BARE and lowered == 'stop_':
            # STAR's explicit loop terminator.  mmCIF never writes one; reading it costs a branch.
            close_loop()
            continue

        is_tag = kind is _BARE and token[:1] == '_'

        if in_loop_header:
            if is_tag:
                if lowered in loop_tags:
                    log.append(LogRecord('star:duplicate-loop-tag', (),
                                         f'record: loop_ header repeats {lowered} on line {lineno}; '
                                         f'the first column of that name is the one read'))
                loop_tags.append(lowered)
                continue
            in_loop_header = False  # first value ends the header; fall through and store it

        if loop_tags is not None:
            if is_tag:
                close_loop()
                # fall through to the scalar-item branch below
            else:
                loop_values.append(_value(kind, token))
                continue

        if is_tag:
            if pending_tag is not None:
                log.append(LogRecord('star:item-no-value', (),
                                     f'record: item {pending_tag} on line {pending_line} has no '
                                     f'value',
                                     LOST))
            pending_tag = lowered
            pending_line = lineno
            continue

        if pending_tag is None:
            orphans += 1
            if orphan_first is None:
                orphan_first = (f'record: value {str(token)[:20]!r} on line {lineno} belongs to no '
                                f'item')
            continue
        target = ensure_block()
        if pending_tag in target.items:
            log.append(LogRecord('star:duplicate-item', (),
                                 f'record: duplicate item {pending_tag} on line {lineno}; the first '
                                 f'value is the one kept'))
        else:
            target.items[pending_tag] = _value(kind, token)
        pending_tag = None

    close_loop()
    if pending_tag is not None:
        log.append(LogRecord('star:item-no-value', (),
                             f'record: item {pending_tag} on line {pending_line} has no value',
                             LOST))
    if orphan_first is not None:
        log.append(LogRecord('star:orphan-value', (),
                             orphan_first if orphans == 1
                             else f'{orphan_first} (and {orphans - 1} more value(s))',
                             LOST))
    if not saw_data_block and block is None:
        # No `data_` anywhere means no block, whatever the text is; nothing scores how CIF-like it
        # looked.  A file whose items merely precede its first `data_` has a block and its own line.
        log.append(LogRecord('star:no-data-block', (),
                             'record: no data_ block is stated anywhere in the text, so there is '
                             'no CIF here; no block is read',
                             LOST))
    if block is not None:
        yield block
