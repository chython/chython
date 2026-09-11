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
"""Every file the installed package opens at runtime is declared as package data.

`pyproject.toml` sets `include-package-data = false`, which is deliberate: the wheel then holds exactly
what `[tool.setuptools.package-data]` names, instead of whatever `MANIFEST.in` happened to sweep into the
sdist. The cost of that choice is that adding a runtime data file and forgetting to name it produces a
wheel that imports fine and fails only when that file is first opened -- on a user's machine, in a code
path a test suite run from a source checkout can never reach, because in a checkout the file is simply
there.

This test is that missing gate. It reads the resource names out of the source and checks each one against
the declarations, so the failure lands on the developer who added the file rather than on the person who
installed the wheel.

WHY THE TOML IS PARSED BY REGEX. `tomllib` arrived in 3.11 and this package supports 3.10, so on the
oldest supported interpreter a `tomllib`-based version of this test would skip -- which is the same as not
having it, since 3.10 is a version people build wheels on. The block being read is small, hand-maintained
and flat, so a regex is enough; if it ever stops being any of those, this test failing to find a
declaration it can see with its own eyes is the signal to switch to a real parser.
"""
from pathlib import Path
from pytest import skip
from re import DOTALL, finditer, search


ROOT = Path(__file__).resolve().parent.parent.parent
PACKAGE = ROOT / 'chython'


def _declared():
    """package name -> list of declared filename patterns, from pyproject.toml."""
    # `encoding='utf-8'` here and on the source sweep below: the tree is UTF-8 and `read_text` without
    # it asks the locale, which is cp1252 on the Windows runner -- where `depict/field.py`'s `∇` raises.
    text = (ROOT / 'pyproject.toml').read_text(encoding='utf-8')
    block = search(r'^\[tool\.setuptools\.package-data\]\n(.*?)(?=^\[|\Z)', text, DOTALL | 8)
    assert block, 'pyproject.toml has no [tool.setuptools.package-data] section'

    out = {}
    for line in block.group(1).splitlines():
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        key, _, value = line.partition('=')
        package = key.strip().strip('\'"')
        out[package] = [m.group(1) for m in finditer(r'[\'"]([^\'"]+)[\'"]', value)]
    return out


def _resource_reads():
    """(package, filename or None) for every `files(...).joinpath(...)` in the source.

    This is the one shape in the tree that reads a file shipped INSIDE the package, as opposed to a file
    the user names -- `open(file)` in a reader takes a path from the caller and is not package data.

    The name is None when the argument is not a literal. `libinchi/wrapper.py` picks its filename from
    `sys.platform` before joining it, and a scanner that only matched quoted strings would skip that read
    entirely -- reporting a clean sweep while covering half the resources in the tree. So a computed name
    is reported rather than dropped, and checked in the weaker form its shape allows: the package must
    declare SOMETHING. That catches the failure this whole test exists for -- a package that ships no data
    at all -- and leaves the exact filenames to the reading code, which is the only thing that knows them.
    """
    found = set()
    for path in PACKAGE.rglob('*.py'):
        if '/test/' in path.as_posix():
            continue
        source = path.read_text(encoding='utf-8')
        if 'joinpath' not in source:
            continue
        package = '.'.join(path.relative_to(ROOT).parent.parts)
        for m in finditer(r'joinpath\(\s*([^)]*?)\s*\)', source):
            argument = m.group(1)
            literal = search(r'^[\'"]([^\'"]+)[\'"]$', argument)
            found.add((package, literal.group(1) if literal else None))
    return found


def _covered(name, patterns):
    """A declaration covers a name if it matches literally or as a glob."""
    from fnmatch import fnmatch
    return any(fnmatch(name, p) for p in patterns)


def test_every_runtime_resource_is_declared_as_package_data():
    declared = _declared()
    missing = []
    for package, name in sorted(_resource_reads(), key=lambda x: (x[0], x[1] or '')):
        patterns = declared.get(package, [])
        if name is None:
            if not patterns:
                missing.append('%s opens a resource whose name it computes, and declares no data at all'
                               % package)
        elif not _covered(name, patterns):
            missing.append('%s opens %r, which no package-data entry covers' % (package, name))
    assert not missing, \
        'these files would be absent from a wheel:\n  ' + '\n  '.join(missing) + \
        '\nadd them to [tool.setuptools.package-data] in pyproject.toml'


def test_every_chemistry_table_is_in_package_data():
    """Every `*.tsv` under `chython/chemistry/tables/` must be declared in package-data.

    `read_table` opens files with `joinpath(f'tables/{name}')` -- an f-string the literal scanner in
    `test_every_runtime_resource_is_declared_as_package_data` cannot see.  That test falls back to the
    weaker check ("the package declares SOMETHING"), so deleting a single table's entry is invisible to
    it.  This check enumerates the directory directly and compares against the declared patterns, so
    the failure lands on the commit that removed the entry rather than on a user running from a wheel.
    """
    tables_dir = PACKAGE / 'chemistry' / 'tables'
    patterns = _declared().get('chython.chemistry', [])
    missing = sorted(
        f.name for f in tables_dir.glob('*.tsv')
        if not _covered(f'tables/{f.name}', patterns)
    )
    assert not missing, (
        'these chemistry tables would be absent from a wheel: %s\n'
        "add them as 'tables/<name>.tsv' to [tool.setuptools.package-data] in pyproject.toml"
        % missing
    )


def test_every_declaration_matches_a_file_that_exists():
    """The other direction: a declaration naming nothing is either a typo or a leftover.

    MATCHED AGAINST THE PATH RELATIVE TO THE PACKAGE, not the bare filename, because that is what
    setuptools matches a package-data pattern against. `chython.chemistry` declares
    `tables/resonance.tsv`, and a bare-filename comparison could neither see the file -- a
    non-recursive `iterdir()` never enters `tables/` -- nor tell `tables/resonance.tsv` from a
    `resonance.tsv` in the package root, which are different declarations shipping different things.
    Relative paths make both exact, and make this test stricter than the filename version was rather
    than looser.

    ONE DELIBERATE EXCEPTION. The three `libinchi.*` names are one library under three platform
    spellings, and only the host platform's is ever present -- the other two matching nothing is the
    intended outcome, not a stale entry. So a declaration is judged by whether its PACKAGE has a match,
    not whether every pattern does.
    """
    empty = []
    for package, patterns in _declared().items():
        directory = ROOT.joinpath(*package.split('.'))
        assert directory.is_dir(), '%s is declared but is not a directory' % package
        names = [f.relative_to(directory).as_posix() for f in directory.rglob('*') if f.is_file()]
        if not any(_covered(name, patterns) for name in names):
            empty.append('%s declares %s and none of them exist' % (package, patterns))
    assert not empty, '\n'.join(empty)


def _staged_lib_dirs():
    """All build/lib* directories present under the repo root."""
    build = ROOT / 'build'
    if not build.is_dir():
        return []
    return list(build.glob('lib*'))


def _staged_files():
    """Every staged file, as (staging directory, path relative to it)."""
    return [(lib, staged.relative_to(lib))
            for lib in _staged_lib_dirs() for staged in sorted(lib.rglob('*')) if staged.is_file()]


def test_staging_holds_nothing_the_source_tree_lost():
    """Every staged file must still have a counterpart in the source tree.

    A stale staged file cannot be caught by any import-level test: it is invisible until someone
    builds a wheel, which is why the gate has to read the directory.  `bdist_wheel` zips everything it
    finds under `build/lib*`, and that directory is never emptied, so a deleted module, a renamed data
    file or a moved package leaves its old copy there to be packaged.  Renaming a package is the sharp
    case: both copies land in one wheel, so the old import path keeps working out of it.
    `prune_stale_staging()` in `setup.py` is what clears them; this test verifies it did.

    Extensions are exempt and are `test_v2_boundary.py`'s business: `_core.cpython-312-darwin.so` has
    no counterpart in the source tree by construction, and neither has `libinchi.*`, which is built
    into `build/inchi/` and staged from there.

    No build directory means a clean checkout -- there are no staged files to be wrong, so the test
    passes vacuously.  This is intentional and noted here because a silent pass for a structural
    reason is exactly the failure mode the rest of this file guards against; the negative control
    below is what keeps the two apart.
    """
    offenders = [str(lib.relative_to(ROOT) / rel) for lib, rel in _staged_files()
                 if not _built_artefact(rel) and not (ROOT / rel).exists()]

    assert not offenders, (
        'these files are staged but no longer exist in the source tree, and a wheel built now would '
        'ship them:\n  ' + '\n  '.join(sorted(offenders))
        + '\nrun `python setup.py build_ext --inplace` to prune them'
    )


def _built_artefact(rel):
    """True for the staged files that legitimately have no counterpart in the source tree."""
    return rel.name.startswith('libinchi.') or rel.name.endswith('.pyd') or '.cpython-' in rel.name


def test_the_staging_gate_can_fail():
    """Negative control for test_staging_holds_nothing_the_source_tree_lost.

    Without it that test passes vacuously the day `_staged_files()` stops finding anything -- the same
    silent-green failure mode `test_the_gate_can_fail` guards against for the resource tests.  It is
    skipped rather than failed on a clean checkout, because "nothing is staged" is a legitimate state
    of the tree and only the scanner going quiet while files are there is a defect.

    A STAGING DIRECTORY HOLDING ONLY BUILT ARTEFACTS IS THAT SAME LEGITIMATE STATE, and this used to
    fail on it.  `python setup.py build_ext --inplace` stages the extension and `libinchi` and runs no
    `build_py`, so in a fresh worktree -- where nobody has built a wheel or an sdist -- `build/lib*`
    holds exactly two files, both of them correctly exempt.  Asserting a `.py` among them reported "the
    scanner has stopped working" about a scanner that had just found every file there was, which is a
    false alarm on a state every worktree passes through.  The claim below is therefore made on the
    files the scanner returned: it found some, and among the ordinary modules -- if any are staged at
    all -- at least one is a `.py`.
    """
    staged = _staged_files()
    if not staged:
        skip('no build/lib* staging directory; nothing for the scanner to find')

    ordinary = [rel for _, rel in staged if not _built_artefact(rel)]
    if not ordinary:
        skip('only built artefacts are staged; `build_py` has never run here, which is what a '
             'worktree built with `build_ext --inplace` looks like')

    # the scanner sees real files, and the exemption is narrow enough to leave ordinary modules in
    assert any(rel.suffix == '.py' for rel in ordinary), \
        'the staging scanner found no Python files at all; it has stopped working'
    assert not _built_artefact(Path('chython/core/_structure.py')), 'the exemption is too broad'
    assert _built_artefact(Path('chython/core/_core.cpython-312-darwin.so')), \
        'extensions must stay exempt; they have no source counterpart by construction'


def test_the_gate_can_fail():
    """A negative control, because both tests above pass vacuously if the scanners find nothing.

    Without this, deleting the body of `_resource_reads` would leave a green suite.
    """
    reads = _resource_reads()
    assert reads, 'the resource scanner found nothing; it has stopped working'
    assert _declared(), 'the pyproject scanner found nothing; it has stopped working'

    # both shapes, because the computed-name branch is the one that was missing at first and the one a
    # future simplification would drop again
    assert any(name is not None for _, name in reads), 'the scanner no longer sees literal names'
    assert any(name is None for _, name in reads), 'the scanner no longer sees computed names'

    # a name nobody declares must be reported as missing
    assert not _covered('definitely-not-shipped.tsv', _declared()['chython.core'])
    # and one that is declared must be accepted
    assert _covered('valence_rules.tsv', _declared()['chython.core'])

    # A DECLARATION IN A SUBDIRECTORY IS MATCHED AS A PATH. Both halves, because the mistake here is
    # silent in both directions: a bare filename must NOT satisfy a `tables/` declaration (that is the
    # wheel-ships-nothing failure), and the relative path must.
    chemistry = _declared()['chython.chemistry']
    assert not _covered('resonance.tsv', chemistry), \
        'a bare filename satisfies a `tables/` declaration; the two are different files to setuptools'
    assert _covered('tables/resonance.tsv', chemistry)
