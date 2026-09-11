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
"""The mechanical half of the code standard, enforced instead of re-swept.

No linter is installed, and `pycodestyle` would not read the `.pxi` files if one were: these checks
cover every tracked text file under `chython/` and `docs/`, Cython included.  What is NOT here is the
prose half -- crisp comments, no chython-2 memoir -- because the words that would spell it out are
also legitimate: `_pach.pxi` says "v2" 52 times about a WIRE FORMAT version.  CLAUDE.md states that
half and the diff is where it is checked.

Each deliberate exception is a named entry with a reason, so widening one is a visible edit.
"""
from collections import Counter
from pathlib import Path
from re import compile as re_compile
from subprocess import run


ROOT = Path(__file__).resolve().parent.parent.parent

#: Lines over 120 columns are generated data, and wrapping them would change what they generate.
WIDE_BY_CONSTRUCTION = {
    'chython/core/_smiles_read.pxi': 'the element perfect-hash table',
    'chython/core/test/v3_fixtures.py': 'base64 arena fixtures',
    'chython/core/test/v4_fixtures.py': 'arena version-4 oracle bytes',
}

#: `# --- label ---...` block rules, in the two spellings the tree uses (with and without a closing `#`).
SECTION_RULE = re_compile(r'^(#+ --- .*?) (-+)( #)?$')

#: The characters docutils accepts as a section adornment.  A run of one of them under a line of text is
#: an underline, and it titles a section only if it covers that text.
ADORNMENT = frozenset('=-~^"\'`:.+*#_<>')


def _tracked(*suffixes):
    """Tracked text files under `chython/` and `docs/`, which are the only source directories.

    `git ls-files` and not a glob: an untracked scratch file is not held to the standard, and a file
    `.gitignore` covers is not one a sweep would have reached either.
    """
    out = []
    for path in run(['git', 'ls-files', 'chython', 'docs'], cwd=ROOT,
                    capture_output=True, text=True).stdout.split('\n'):
        if not path or not path.endswith(suffixes):
            continue
        try:
            out.append((path, (ROOT / path).read_text(encoding='utf-8')))
        except (UnicodeDecodeError, FileNotFoundError):
            continue
    assert out, 'the git ls-files read found nothing; it has stopped working'
    return out


TEXT = ('.py', '.pxi', '.pyx', '.pxd', '.tsv', '.md', '.rst', '.txt', '.js', '.yml', '.cfg')
CODE = ('.py', '.pxi', '.pyx', '.pxd')


def test_every_file_ends_in_exactly_one_newline():
    """A missing one makes the next append land on the last line; a second one is a blank line in a diff.

    Empty files are exempt because an empty file has no line to terminate.
    """
    missing = [p for p, s in _tracked(*TEXT) if s and not s.endswith('\n')]
    extra = [p for p, s in _tracked(*TEXT) if s.endswith('\n\n')]
    assert not missing, 'no final newline:\n  %s' % '\n  '.join(missing)
    assert not extra, 'blank line at end of file:\n  %s' % '\n  '.join(extra)


def test_no_carriage_returns():
    """Data files under the repo root are fixed-column formats where CRLF is correct; source is not."""
    offenders = [p for p, s in _tracked(*TEXT) if '\r' in s]
    assert not offenders, 'CRLF line endings:\n  %s' % '\n  '.join(offenders)


def test_no_trailing_whitespace():
    """Invisible, and it makes a whitespace-only diff hunk out of an unrelated edit.

    A TSV is checked for trailing SPACES only: a trailing tab there is a required empty column, so
    `residues.tsv`'s rows end in tabs by design.
    """
    offenders = []
    for path, src in _tracked(*TEXT):
        strip = ' \t' if not path.endswith('.tsv') else ' '
        for n, line in enumerate(src.split('\n'), 1):
            if line.strip() and line != line.rstrip(strip):
                offenders.append('%s:%d' % (path, n))
    assert not offenders, 'trailing whitespace:\n  %s' % '\n  '.join(offenders[:40])


def test_blank_lines_are_empty_and_never_run_to_three():
    """Two blank lines separate top-level definitions; three separate nothing, and a blank line that
    holds spaces is one no editor shows."""
    padded, runs = [], []
    for path, src in _tracked(*TEXT):
        for n, line in enumerate(src.split('\n'), 1):
            if line and not line.strip():
                padded.append('%s:%d' % (path, n))
        if '\n\n\n\n' in src:
            runs.append(path)
    assert not padded, 'whitespace-only line:\n  %s' % '\n  '.join(padded[:40])
    assert not runs, 'three or more consecutive blank lines:\n  %s' % '\n  '.join(runs)


def test_code_stays_within_120_columns():
    """The convention `pycodestyle --max-line-length=120` would enforce, extended to Cython.

    `pycodestyle` reads only `.py`, which left the `.pxi` half of the tree unlinted.
    """
    offenders = []
    for path, src in _tracked(*CODE):
        if path in WIDE_BY_CONSTRUCTION:
            continue
        offenders += ['%s:%d (%d columns)' % (path, n, len(line))
                      for n, line in enumerate(src.split('\n'), 1) if len(line) > 120]
    assert not offenders, (
        'over 120 columns:\n  %s\n'
        'Wrap it, or -- only if the line is generated data -- name the file in WIDE_BY_CONSTRUCTION '
        'with the reason.' % '\n  '.join(offenders[:40]))


def test_the_named_wide_files_are_still_wide():
    """The other direction: an exception whose reason has expired is one to delete, not to carry."""
    for path, reason in WIDE_BY_CONSTRUCTION.items():
        src = (ROOT / path).read_text(encoding='utf-8')
        assert any(len(line) > 120 for line in src.split('\n')), \
            '%s is within 120 columns now; drop it from WIDE_BY_CONSTRUCTION (%s)' % (path, reason)


def test_section_rule_comments_line_up():
    """`# --- label -----` block rules are padded to one column per file, and drift breaks the column.

    Only a CLUSTER is judged: where a file's rules span more than three columns they are labels of
    different widths rather than one padded block, and a file with a single rule has no column to keep.
    An off-by-one in a block of otherwise equal rules is the drift this catches.
    """
    offenders = []
    for path, src in _tracked(*CODE, '.tsv'):
        widths = Counter(len(line) for line in src.split('\n') if SECTION_RULE.match(line))
        if len(widths) < 2 or max(widths) - min(widths) > 3:
            continue
        column, _ = widths.most_common(1)[0]
        offenders.append('%s: rules at %s, expected all %d'
                         % (path, sorted(widths), column))
    assert not offenders, 'section rules out of column:\n  %s' % '\n  '.join(offenders)


def test_rst_section_underlines_cover_their_titles():
    """A short underline is a build warning, and it is one nobody reading the page can see.

    An inline literal counts by its source width, which is where the off-by-one comes from: a heading
    spelled ``foo()`` is 9 columns to the rule and 5 to a reader. Two is the shortest adornment docutils
    recognises; a run holding a space is a table border and not an underline.
    """
    offenders = []
    for path, src in _tracked('.rst'):
        lines = src.split('\n')
        for n, line in enumerate(lines[1:], 2):    # `n` is 1-based, so the title above is `lines[n - 2]`
            rule = line.rstrip()
            if len(rule) < 2 or len(set(rule)) != 1 or rule[0] not in ADORNMENT:
                continue
            title = lines[n - 2].rstrip()
            if not title or (len(set(title)) == 1 and title[0] in ADORNMENT):
                continue                           # a transition, or the overline of an overlined title
            if len(rule) < len(title):
                offenders.append('%s:%d (%d columns under a %d-column title)'
                                 % (path, n, len(rule), len(title)))
    assert not offenders, 'section underline shorter than its title:\n  %s' % '\n  '.join(offenders)
