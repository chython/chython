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
"""Every Python sample in ``docs/`` is executed, and no sample may opt out by being unexecutable.

This is not a style gate.  It is the comparator between two representations of one API -- the one the
library has and the one the documentation claims -- and the class of defect it catches (absent methods,
absent facade names, changed signatures) is invisible to a Sphinx build by construction, because Sphinx
renders a code block without reading it.

**The gate is keyed on the property, not on a directory.**  Three assertions, and the last two are what
keep the first honest:

* every ``testcode::`` block executes without raising;
* every ``testoutput::`` block equals what its sample printed;
* **no ``code-block:: python`` survives anywhere in ``docs/``.**

Without the third, the gate would be trivially defeatable and would defeat itself the first time
someone documented a new feature: a sample written as ``code-block`` renders identically, runs never,
and reports nothing.  "Which directive did you use" is therefore not left to an author's memory.
``code-block:: bash`` is untouched -- this gate is about Python -- and so is anything under
``docs/_build``, which is output.

Without the second, a page could run every sample and still print a number no reader would ever see: a
``testoutput`` is compared by ``sphinx.ext.doctest`` and by nothing a test run invokes, so an expected
output is only as true as the last person to read it.  Comparison here is exact after trailing
whitespace and surrounding blank lines go, since a sample whose output needs a wildcard to match is a
sample stating something it does not know.  Where the two comparators could then disagree -- an rst body
may not carry a trailing space, which ``sphinx -b doctest`` nevertheless demands -- the sample is what
must change, so a printed line ending in whitespace fails here rather than passing here and failing
there.

State is shared **down a page** and never across pages, which is how the pages are actually written: a
molecule parsed in the first sample is used by the fifth.  A page is therefore all-or-nothing, and a
failure names the file and the line the block starts on so the sample is one click away.

A sample that cannot run in a bare checkout does not get an exemption here; it gets a ``testsetup::``
block that makes it runnable.  Six samples read files (``molecules.sdf``, ``reactions.rdf``) that a
``testsetup`` now writes into a temporary directory, which is strictly better than the alternative of
letting them not run: the sample is true, and the file it reads is one the reader can see being made.
The one genuine exception is a sample needing a Java JAR, and it carries ``:skipif:`` with the reason
in the directive rather than silence.
"""

from pathlib import Path
from pytest import fail, mark, skip


#: ``.. <directive>::`` at the start of a line, with its indentation and argument.
_DIRECTIVE = '.. '

#: The directives whose bodies this gate executes, in the order a page's blocks must run.
#: ``testsetup`` first within a page is Sphinx's rule too, and it is what lets a sample that reads a
#: file be honest about where the file came from.
_EXECUTED = ('testsetup', 'testcode')

#: The directive holding what the preceding executed block must have printed.
_COMPARED = 'testoutput'

#: The directive that must NOT appear with a Python argument.  See the module docstring: a sample
#: written this way renders identically to an executed one and is silently never run.
_UNGATED = 'code-block'

#: Languages that mean "this is Python" to Sphinx.
_PYTHON = frozenset(('python', 'python3', 'py'))


def _doc_root():
    """``docs/`` beside the repository's ``chython/``, or ``None`` in an installed package."""
    for parent in Path(__file__).resolve().parents:
        if (parent / 'chython').is_dir() and (parent / 'docs').is_dir():
            return parent / 'docs'
    return None


def _pages():
    root = _doc_root()
    return sorted(root.glob('*.rst')) if root is not None else []


def _blocks(path):
    """``[(directive, language, line_number, source)]`` for one page, in file order.

    The body of a directive is the indented run that follows it, with option lines (``:name: value``)
    and leading blanks skipped.  Dedented to the body's own first-line indent so it compiles.

    `encoding='utf-8'`: the pages are UTF-8 and `docs/depiction.rst` holds a `⁻`, which the locale codec
    Windows hands `read_text` cannot decode.
    """
    lines = path.read_text(encoding='utf-8').splitlines()
    out = []
    i = 0
    while i < len(lines):
        stripped = lines[i].lstrip()
        if not stripped.startswith(_DIRECTIVE) or '::' not in stripped:
            i += 1
            continue
        head, _, argument = stripped[len(_DIRECTIVE):].partition('::')
        directive = head.strip()
        if directive not in _EXECUTED and directive not in (_COMPARED, _UNGATED):
            i += 1
            continue
        outer = len(lines[i]) - len(stripped)
        start = i + 1
        options = {}
        j = i + 1
        while j < len(lines):
            option = lines[j].strip()
            if not option:
                j += 1
            elif option.startswith(':') and option.count(':') >= 2:
                key, _, value = option[1:].partition(':')
                options[key.strip()] = value.strip()
                j += 1
            else:
                break
        body, base = [], None
        while j < len(lines):
            line = lines[j]
            if not line.strip():
                body.append('')
                j += 1
                continue
            indent = len(line) - len(line.lstrip())
            if base is None:
                if indent <= outer:
                    break
                base = indent
            elif indent < base:
                break
            body.append(line[base:])
            j += 1
        out.append((directive, argument.strip().lower(), start, '\n'.join(body).rstrip(), options))
        i = j
    return out


def test_doc_has_no_ungated_python_sample():
    """The assertion that keeps the executing one honest -- see the module docstring.

    A ``code-block:: python`` renders exactly like a ``testcode::`` and runs never.  Leaving the choice
    to an author's memory is how 112 unexecuted samples accumulated in the first place, so the choice
    is not left to memory.

    [mutant: changing a page's ``testcode::`` back to ``code-block:: python`` -- this test fails and
    names the page and line, where the executing test below simply stops covering that sample and stays
    green, which is the whole point of having both.]
    """
    pages = _pages()
    if not pages:
        skip('no docs/ beside this package -- an installed copy, not a checkout')

    ungated = [f'{p.name}:{line}' for p in pages
               for directive, language, line, _, _ in _blocks(p)
               if directive == _UNGATED and language in _PYTHON]
    assert not ungated, ('these Python samples are never executed; write them as `.. testcode::` so '
                         f'this suite runs them: {", ".join(ungated)}')


def _comparable(text):
    """One printed or expected block, with trailing whitespace and surrounding blanks gone.

    Neither difference is visible to a reader of the rendered page, and an rst body carries the second
    by construction -- the directive's blank line and the dedent leave them behind.
    """
    return '\n'.join(line.rstrip() for line in text.strip().splitlines())


@mark.parametrize('page', [p.name for p in _pages()] or ['<no docs/>'])
def test_every_documented_sample_runs(page):
    """Each page's executed blocks, in order, in one namespace and a temporary directory.

    One test per **page** rather than per block, because the pages share state down their length on
    purpose -- a molecule parsed in the first sample is the subject of the fifth.  Splitting per block
    would either re-run every predecessor or report a cascade of `NameError`s naming the wrong sample.
    The failure message names the block that raised, so the granularity of the *report* is still the
    block.

    A `testoutput::` is compared against what the block before it printed.  Its own `:skipif:` mirrors
    that block's, which is how a page states a sample needing an optional dependency; a `testoutput`
    whose sample did not run and which carries no such condition is the page contradicting itself and
    fails as one.

    [mutant: changing one digit of any `testoutput` body -- this test fails and names the page, the
    line and both texts.  Before the comparison existed a wrong expected output was caught by nothing a
    test run invokes.]
    """
    root = _doc_root()
    if root is None:
        skip('no docs/ beside this package -- an installed copy, not a checkout')

    from contextlib import redirect_stdout
    from io import StringIO
    from os import chdir, getcwd
    from tempfile import TemporaryDirectory

    path = root / page
    blocks = [b for b in _blocks(path) if b[0] in _EXECUTED or b[0] == _COMPARED]
    if not blocks:
        skip(f'{page} documents no Python sample')

    namespace = {'__name__': f'doc_{page.replace(".", "_")}'}
    printed = None  # what the last executed block printed, or None when it was skipped
    was = getcwd()
    with TemporaryDirectory() as scratch:
        try:
            chdir(scratch)
            for directive, language, line, source, options in blocks:
                if 'skipif' in options:
                    try:
                        if eval(options['skipif'], dict(namespace)):
                            if directive != _COMPARED:
                                printed = None
                            continue
                    except Exception as e:  # a broken condition must not read as a skip
                        fail(f'{page}:{line} has an unevaluable :skipif: -- {e!r}')
                if directive == _COMPARED:
                    if printed is None:
                        fail(f'{page}:{line} states the output of a sample that did not run')
                    if _comparable(printed) != _comparable(source):
                        fail(f'{page}:{line} states an output the sample does not print.\n\n'
                             f'expected:\n{source}\n\nprinted:\n{printed.rstrip()}')
                    # `strip('\n')` and not `strip()`: the trailing space this looks for is usually on
                    # the LAST line, which a full strip would remove before the check could see it.
                    ragged = [n for n, text in enumerate(printed.strip('\n').splitlines(), 1)
                              if text != text.rstrip()]
                    if ragged:
                        fail(f'{page}:{line} cannot state what the sample prints: line(s) '
                             f'{", ".join(str(n) for n in ragged)} end in whitespace, which an rst body '
                             'may not carry and `sphinx -b doctest` compares exactly.  Make the sample '
                             'print no trailing space -- a slice that ends mid-gap is the usual cause.')
                    continue
                captured = StringIO()
                try:
                    with redirect_stdout(captured):
                        exec(compile(source, f'{page}:{line}', 'exec'), namespace)
                except BaseException as e:
                    fail(f'{page}:{line} ({directive}) raised {type(e).__name__}: {e}\n\n{source}')
                printed = captured.getvalue()
        finally:
            chdir(was)
