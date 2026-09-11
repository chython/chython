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
"""Every append to a log emits a `LogRecord`.  The count of the sites that do not only goes down.

`by_severity`, `by_stage`, `repaired()` and `lost()` are what `core/_log.py` exists to make possible, and
a bare sentence arrives with `severity=INFO` and `stage=''` -- invisible to all four.  518 sites across
33 files, so the conversion was ratcheted per file; the allowance is now empty and stays empty.

A COUNT AND NOT A PER-SITE ALLOW-LIST, on purpose: a substring allow-list of that many entries would be
longer than the diff it guards, and a count cannot be satisfied by moving a line.
"""
from pathlib import Path
from re import compile


ROOT = Path(__file__).resolve().parent.parent          # chython/

#: `<anything>log.append(` with an argument that is not a record.  `_MC_RECORD` and `mc_record` are the
#: `.pxi` spellings of `LogRecord`, which cannot be imported at module level there.  Tuple and list
#: literals are excluded: `out.append((r, code))` in `core/wedge.py` is a return value, not a log.  So are
#: `_as_bytes` and `blob_bytes`, which fill an S-group's own `log` list of byte blobs -- a different `log`.
_NOT_A_RECORD = r'\s*(?!LogRecord\b|_MC_RECORD\b|mc_record\(|_as_bytes\(|blob_bytes\()(?![(\[])'
_LOG_APPEND = compile(r'\b\w*log\.append\(' + _NOT_A_RECORD)
#: The same, for the modules that alias the log as `out`.
_OUT_APPEND = compile(r'\bout\.append\(' + _NOT_A_RECORD)
_ALIAS = compile(r'out = \[\] if log is None else log')

#: Sites still appending a non-`LogRecord`, per file.  Lower the number, or delete the entry at zero.
#: Empty, and `test_the_ratchet_is_closed` is what keeps it that way.
_REMAINING = {}


def _sources():
    return sorted(p for p in ROOT.rglob('*')
                  if p.suffix in ('.py', '.pxi', '.pyx') and 'test' not in p.parts)


def _count(path):
    source = path.read_text(encoding='utf8')
    n = len(_LOG_APPEND.findall(source))
    if _ALIAS.search(source):
        n += len(_OUT_APPEND.findall(source))
    return n


def test_no_file_appends_more_bare_records_than_it_is_allowed():
    over = []
    for path in _sources():
        name = str(path.relative_to(ROOT))
        found, allowed = _count(path), _REMAINING.get(name, 0)
        if found > allowed:
            over.append(f'{name}: {found} bare append(s), {allowed} allowed')
    assert not over, 'the ratchet only turns one way:\n' + '\n'.join(over)


def test_no_allowance_is_stale():
    """A number that is too high is a ratchet that has stopped ratcheting."""
    stale = []
    for name, allowed in sorted(_REMAINING.items()):
        found = _count(ROOT / name)
        if found < allowed:
            stale.append(f'{name}: {found} bare append(s), {allowed} still allowed -- lower it')
    assert not stale, '\n'.join(stale)


def test_the_ratchet_is_closed():
    """Nothing may re-enter through the allow-list: the conversion is finished."""
    assert _REMAINING == {}, f'still allowed: {sorted(_REMAINING)}'
