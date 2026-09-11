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
"""Every import of chython 2 from outside chython 2, enumerated -- and there are none.

WHAT THIS GATE HOLDS SHUT is the direction: nothing outside chython 2 imports chython 2, and chython
builds exactly one extension.

WHY A GATE AND NOT A DOCUMENT.  A list in prose goes stale in a week, and it goes stale *silently in
the wrong direction*: a new chython 2 import nobody decided to add reads as an unchanged document.

THE V2 SOURCE IS NOT IN THE WORKING TREE -- `algorithms/`, `containers/`, `files/` and `reactor/` are
deleted, unported modules included, because a directory kept "as reference" is a directory a
half-ported module can reach into for a valence table or a mixin.  What reads them is git
(`git show 5e39eb5:<path>`) and, for what was and was not ported,
`docs/superpowers/research/2026-09-03-v2-ported-audit.md`.  That does not retire this scan: restoring
one of those directories out of git to get at a table is the cheapest possible way to re-create the
dependency, and the scan below is the thing that makes it fail.

The property is a ratchet, and it is TWO assertions, not one:

* an edge that is not in `ALLOWED` fails -- you cannot add a dependency on chython 2 by accident,
  only by writing a line in this file that says who you are and why;
* an entry in `ALLOWED` that no longer matches any edge ALSO fails -- when you cut a dependency you
  must delete its line, so the ledger shrinks with the code instead of describing a tree that has
  moved on.

`ALLOWED` is empty, so the first assertion is the whole gate and the second only fires if someone
declares an edge and then removes it without tidying up.  Both are kept: an empty `ALLOWED` is a state
this file can leave and come back to.

WHAT COUNTS AS AN EDGE.  A module outside the chython 2 packages naming one of them in an `import`
statement, at module level or inside a function -- a deferred import is still an import, so a scanner
reading only the top of each file reports a tree cleaner than the one that exists.

WHAT DOES NOT COUNT.  Mentioning a V2 name in a string, a comment or a docstring: prose is not a
dependency.  Neither is a subprocess: `chython/core/test/oracle.py` runs an out-of-tree chython 2.24
under `-I` to keep the differential tests alive, and it is not on this ledger because spawning an
interpreter is not an import, which is what lets those tests run with V2 off the import graph.

ONE PACKAGE AND ONE EXTENSION.  The import ledger holds the Python half; the extension ledger further
down holds the compiled half.  `test_one_chython_one_so` at the bottom is the conjunction.
"""
from ast import Call, Constant, Import, ImportFrom, Name, parse, walk
from pathlib import Path
from re import compile as re_compile


ROOT = Path(__file__).resolve().parent.parent.parent
PACKAGE = ROOT / 'chython'

# THE FORBIDDEN NAMES, and none of them is a directory.  They all stay listed on purpose.  An import of
# a package that is not there raises `ModuleNotFoundError`, so the gate looks redundant; what it
# actually catches is somebody restoring a directory out of git to get at a table, which is the cheapest
# possible way to re-create the dependency this file exists to prevent.
V2_PACKAGES = frozenset({'algorithms', 'containers', 'files', 'periodictable', 'reactor', 'utils'})

# (importing module, chython 2 package) -> why it is still there, and whose job it is.
#
# EMPTY.  A reason field is not decoration: it says whether the next person is looking at a real feature
# gap, a test that needs an out-of-process oracle, or a line somebody forgot, and those three want
# completely different work.  A test oracle in particular does NOT belong here -- it reaches chython 2
# through `chython/core/test/oracle.py`, which spawns an out-of-tree 2.24 under `-I`, pins the version
# and asserts the child's `chython.__file__` lies outside this repository, so an oracle cannot silently
# become a test of this tree against itself.
#
# Adding an entry back is legitimate -- write the pair and the reason -- but it is a decision somebody
# makes in this file, which is the whole point.
ALLOWED = {}


def _imports_into(packages):
    """{(importing module, subpackage)} for every edge into `packages` from outside them.

    Parametrised on the subpackage set for one reason, and it is not reuse: with `ALLOWED` empty, every
    assertion that only reads `_v2_imports()` passes just as happily when the scanner is BROKEN as when
    the tree is clean, and a walk that finds no files at all is indistinguishable from success.  Pointing
    the same scan at `chython.core` -- a package the whole library imports and always will -- gives a
    control that cannot go quiet, which is what `test_the_scan_of_the_real_tree_still_works` uses.
    """
    found = set()
    for path in sorted(PACKAGE.rglob('*.py')):
        parts = path.relative_to(ROOT).with_suffix('').parts
        is_init = parts[-1] == '__init__'
        module = '.'.join(parts[:-1] if is_init else parts)
        # the package a relative import counts up from: a package's __init__ IS its package
        package = module if is_init else '.'.join(parts[:-1])
        if module.split('.')[1:2] and module.split('.')[1] in packages:
            continue  # chython 2 importing itself is not an edge to cut, it dies wholesale

        # `encoding='utf-8'` on every source read in this file: the tree is UTF-8, and `read_text`
        # without it asks the locale -- cp1252 on the Windows runner, where `depict/field.py`'s `∇`
        # raises and a scanner that cannot read a file reports a clean tree instead of a failure.
        for node in walk(parse(path.read_text(encoding='utf-8'))):
            for target in _targets(node, package):
                named = target.split('.')
                if len(named) > 1 and named[0] == 'chython' and named[1] in packages:
                    found.add((module, named[1]))
    return found


def _v2_imports():
    """{(importing module, chython 2 package)} for every edge into chython 2 from outside it."""
    return _imports_into(V2_PACKAGES)


def _targets(node, package):
    """Absolute module names one import statement refers to.

    RELATIVE IMPORTS ARE THE WHOLE DIFFICULTY, and getting them wrong is silent: an off-by-one in the
    level arithmetic resolves `from ..algorithms import x` to `chython.containers.algorithms`, which
    matches no V2 package, so the edge vanishes and the scan reports a clean tree.  `level` counts up
    from the importing module's PACKAGE, so level 1 is that package and each further level strips one
    more part -- and `test_the_scanner_resolves_relative_imports` pins a real four-level case.
    """
    if isinstance(node, Import):
        return [alias.name for alias in node.names]
    if not isinstance(node, ImportFrom):
        return []
    if node.level:
        parts = package.split('.')
        up = node.level - 1
        base = parts[:len(parts) - up] if up <= len(parts) else ['chython']
        prefix = '.'.join(base)
        module = f'{prefix}.{node.module}' if node.module else prefix
    else:
        module = node.module or ''
    # the module itself, and each name after `import`, since any of them may be a submodule
    return [module] + [f'{module}.{alias.name}' for alias in node.names]


def test_no_new_dependency_on_chython_2():
    """An import of chython 2 from outside it must be declared here first.

    With `ALLOWED` empty this is the gate: "nothing imports chython 2" is a property that has to keep
    being true rather than one that was achieved, and the V2 source is one `git checkout` away from
    being importable again.
    """
    undeclared = sorted(_v2_imports() - set(ALLOWED))
    assert not undeclared, (
        'these modules import chython 2 and are not in ALLOWED:\n  '
        + '\n  '.join('%s -> chython.%s' % edge for edge in undeclared)
        + '\n\nchython 2 is being deleted. If this dependency is genuinely needed for now, add it to '
          'ALLOWED in this file with a reason saying who removes it and when. If it is not, use the '
          'chython 3 equivalent instead.')


def test_the_ledger_has_no_stale_entries():
    """THE RATCHET. Cutting a dependency means deleting its line, so the ledger shrinks with the code.

    Without this, `ALLOWED` would only ever grow and a tree with two edges left would look exactly like
    one with twenty.  With it, the size of this dict is the remaining distance.

    It is trivially true on an empty `ALLOWED` and is kept anyway, because an empty ledger is a state
    this file can leave: the next temporary edge somebody declares is one nobody would remember to
    delete if the ratchet had been taken out on the grounds that it had nothing left to measure.
    """
    stale = sorted(set(ALLOWED) - _v2_imports())
    assert not stale, (
        'these ALLOWED entries no longer match any import -- the dependency is gone:\n  '
        + '\n  '.join('%s -> chython.%s' % edge for edge in stale)
        + '\n\nDelete them. That is the point of this test: the ledger is only useful if it shrinks '
          'when the code does.')


def test_the_scanner_resolves_relative_imports():
    """A negative control on the one thing whose failure mode is a clean-looking report.

    Both ledger tests pass vacuously if `_v2_imports` returns nothing, and pass *misleadingly* if it
    silently drops relative imports -- so the resolver is exercised directly, on sources written here.

    THE CASES ARE WRITTEN HERE AND NOT PINNED TO A REAL EDGE.  Any real edge is one somebody is paid to
    delete, so pinning one makes the control fail the day the tree improves -- the one reason a control
    must never fail.

    THE LEVEL MUST BE >= 2 or the case proves nothing.  With `level == 1` the off-by-one gives the same
    answer as the correct arithmetic, so a one-dot import passes either way; only two or more dots
    separate them.  Both spellings of the error are covered below: counting from the parent when the
    importer is a package `__init__`, and counting from the module itself when it is not.
    """
    cases = [
        # (importing module, is it a package __init__, source, what it must resolve to)
        ('chython.formats.ctfile.test.conftest', False,
         'from ....files.SDFrw import SDFRead', {'files'}),
        ('chython.depict.grid', False, 'from ..containers import MoleculeContainer', {'containers'}),
        # a package's __init__ counts up from ITS OWN package, not its parent -- the off-by-one here
        # makes every edge out of a subpackage invisible
        ('chython.formats', True, 'from .ctfile import mol', set()),
        ('chython.formats', True, 'from ..files import xyz', {'files'}),
        ('chython.core', True, 'from ..periodictable import C', {'periodictable'}),
        # absolute, and a submodule named after `import`
        ('chython.anything', False, 'from chython import containers', {'containers'}),
        ('chython.anything', False, 'import chython.reactor', {'reactor'}),
        # not edges: chython 3 packages, and a V2 name that is only a string
        ('chython.anything', False, 'from ..core import MoleculeContainer', set()),
        ('chython.anything', False, 'x = "from chython.containers import MoleculeContainer"', set()),
    ]
    for module, is_init, source, expected in cases:
        package = module if is_init else module.rsplit('.', 1)[0]
        found = set()
        for node in walk(parse(source)):
            for target in _targets(node, package):
                named = target.split('.')
                if len(named) > 1 and named[0] == 'chython' and named[1] in V2_PACKAGES:
                    found.add(named[1])
        assert found == expected, \
            'in %s (%s): %r resolved to %s, expected %s' % (
                module, 'package' if is_init else 'module', source, sorted(found), sorted(expected))


def test_the_scan_of_the_real_tree_still_works():
    """The control on the scanner itself, and it cannot be a comparison against `ALLOWED`.

    `bool(_v2_imports()) == bool(ALLOWED)` holds for free with both sides empty, and every other
    assertion in this file reads `_v2_imports()` expecting nothing -- so a scanner that walked no files,
    failed to parse, or resolved every import to the wrong name would turn the entire module green.

    So the scan is run against `chython.core` instead.  Every layer above the core imports it, by
    construction and for as long as there is a library, so a scan that comes back empty for that prefix
    is a broken scan and nothing else.  Both spellings are checked, because they fail separately: the
    walk finding files at all, and `_targets` resolving what it finds.
    """
    edges = _imports_into(frozenset({'core'}))
    assert edges, 'the scanner found no imports of chython.core; the walk or the parse is broken'

    # Most in-tree imports of the core are RELATIVE (`from ..core import ...`), and they are spread over
    # every layer above it.  If `_targets` had lost its relative-import arithmetic the survivors would be
    # the handful of absolute ones, clustered in one or two places -- so breadth is what separates a
    # working resolver from a scan that found only the easy half.
    layers = {module.split('.')[1] for module, _ in edges if module.count('.') > 1}
    assert len(layers) >= 3, (
        'only %s import chython.core; the relative-import arithmetic in _targets has probably '
        'regressed, since every layer above the core imports it' % (sorted(layers) or 'nothing'))


def test_a_package_importing_itself_is_not_counted():
    """chython 2's packages were deleted together, so an edge among them was never work: counting one
    would bury the edges that matter under hundreds that resolve themselves when a directory goes.

    Asserted on `chython.core`, not on chython 2: with no V2 edges left, the V2 spelling of this test
    passed whether the skip worked or not.  The core imports itself constantly -- every `.pxi` layer's
    Python-side neighbour does -- so if the skip were broken this list would be long.
    """
    # A slice and not `[1]`, because the facade is `chython` -- no second part at all -- and it can be
    # one of the importers, so indexing raises rather than answering.
    assert not [module for module, _ in _imports_into(frozenset({'core'}))
                if module.split('.')[1:2] == ['core']]


# ---------------------------------------------------------------------------------------------------
# The compiled half.  One package, one extension: that is the shape chython ships in.
#
# A separate ratchet, because an extension is not deleted by cutting an import -- `setup.py` names it
# explicitly, and dropping the name is the edit.  A `.pyx` that leaves this ledger is read out of git.
TARGET_EXTENSION = 'chython.core._core'

# the interpreter tag every extension filename carries between the module name and the suffix:
# `_core.cpython-310-darwin.so`, `_core.cp312-win_amd64.pyd`, `_core.pypy310-pp73-darwin.so`
_TAGGED = re_compile(r'\.(?:cpython|cp|pypy)-?\d')

EXTENSIONS = {
    TARGET_EXTENSION:
        'the whole core is one translation unit, so every cdef call between its .pxi layers is a '
        'static call the C compiler can inline',
}


def _declared_extensions():
    """Extension module names from `setup.py`, read rather than executed.

    `setup.py` cannot be imported to ask it: it runs `cythonize`, compiles libinchi, and prunes the
    staging directory as import side effects.  The list is a literal, so `ast` reads it exactly and
    for free.
    """
    tree = parse((ROOT / 'setup.py').read_text(encoding='utf-8'))
    names = set()
    for node in walk(tree):
        if isinstance(node, Call) and isinstance(node.func, Name) and node.func.id == 'Extension' \
                and node.args and isinstance(node.args[0], Constant):
            names.add(node.args[0].value)
    return names


def test_the_extension_ledger_matches_setup_py():
    """Adding or removing a compiled module is a declaration, in both directions.

    The removing direction is the ratchet again, and it matters more here than for imports: an
    extension that leaves the `Extension` list does NOT leave the build, because `build/lib*` is never
    emptied and `bdist_wheel` zips whatever is in it, so a deleted `.pyx` can still ship its stale `.so`
    in the next wheel.  `setup.py:prune_stale_staging()` is what empties the staging directory; this
    test is what notices the removal happened at all.
    """
    declared, actual = set(EXTENSIONS), _declared_extensions()
    assert not actual - declared, (
        'setup.py builds extensions this ledger does not declare: %s\n'
        'chython 3 ships one extension. Adding a second needs a reason written here.'
        % sorted(actual - declared))
    assert not declared - actual, (
        'this ledger declares extensions setup.py no longer builds: %s\n'
        'Delete them here too -- and check that no stale .so is left in build/lib* or in the source '
        'tree, because neither is cleaned by removing the Extension entry.' % sorted(declared - actual))


def test_no_stale_extension_in_the_source_tree():
    """An in-place `.so` whose module is no longer built is a trap that answers imports.

    In-place builds put the artefact next to its source, and nothing removes it when the source goes.
    A stale one keeps an import working in the developer's checkout long after the module has been
    deleted -- so the tests pass locally and fail for everyone else.  Only module
    identity is checked, not the ABI tag: several tags for a module that IS still built just means the
    developer builds for several interpreters, which is their business.

    AN EXTENSION IS A SHARED LIBRARY CARRYING AN INTERPRETER TAG, which is `prune_stale_staging()`'s
    rule in `setup.py` and is here for the same reason it is there: `core/libinchi.so` is a shared
    library the in-place build stages beside the extension on Linux, so a rule reading every `.so` as a
    module calls the supported state stale -- on Linux only, since the same file is `libinchi.dylib` on
    macOS.  Both halves of the rule are load-bearing: the tag alone matches every `.pyc` under
    `__pycache__`, and the suffix alone matches libinchi.
    """
    built = _declared_extensions()
    stale = sorted(str(p.relative_to(ROOT)) for p in PACKAGE.rglob('*')
                   if p.suffix in ('.so', '.pyd', '.dll') and _TAGGED.search(p.name)
                   and '.'.join((*p.relative_to(ROOT).parts[:-1], p.name.split('.')[0])) not in built)
    assert not stale, (
        'these compiled modules are in the source tree but built from nothing:\n  ' + '\n  '.join(stale)
        + '\n\nDelete them. They will answer an import that should fail.')


def test_one_chython_one_so():
    """The criterion, stated as a test rather than as an aspiration.

    It asserts a conjunction of the two ledgers above, which is the part neither of them says on its
    own: an empty `ALLOWED` with extra extensions still building is not one package and one extension,
    and neither is one extension with a live import edge.  Both, together, is the shape chython ships in.
    """
    remaining = sorted(set(EXTENSIONS) - {TARGET_EXTENSION})
    v2_imports = sorted(_v2_imports())
    assert not remaining and not v2_imports, (
        'chython 3 is not yet one package and one extension.\n'
        '  extensions still built beside %s: %d\n    %s\n'
        '  imports of chython 2 from outside it: %d\n'
        'Both lists reach zero together, and then chython 2 can be deleted outright.'
        % (TARGET_EXTENSION, len(remaining), '\n    '.join(remaining) or '-', len(v2_imports)))
