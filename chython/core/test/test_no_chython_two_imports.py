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
"""The core depends on nothing above it, and its tests do not reach into chython 2.

THE RULE.  The dependency direction is one-way:

    core <- chemistry <- {featurize, reactions} <- {formats, depict, interop} <- the facade

so `chython.core` may import `chython.core` and the standard library, and nothing else from this
distribution.  In particular NOTHING IMPORTS `chython` ITSELF.  That is not a style preference:
`chython/__init__.py` pulls in about 130 modules, so a single `from chython import smiles` anywhere
under `chython/core/` couples the core's own test suite to every layer above it.  Compiling the core
to a `.so` does not help -- the package `__init__` runs first either way.

WHY THIS IS A TEST RATHER THAN A NOTE IN A DOCUMENT.  Reaching upward for a cheap second opinion is
easy to do one file at a time, and the cost is invisible until the suite cannot run without the thing
it is replacing.  A note does not stop that; a test does.

A SUBPROCESS IS NOT AN IMPORT, and the distinction is the whole point.  `oracle.py` compares against
chython 2, permanently -- it runs it in ANOTHER INTERPRETER against an INSTALLED
copy.  Nothing enters this process's `sys.modules`, the working tree is not on the oracle's
`sys.path`, and deleting chython 2 from this repository changes nothing about it.  So:

    subprocess.run([oracle_interpreter, '-I', '-c', 'import chython.periodictable ...'])   OK
    import chython.periodictable                                                           NOT OK

Do not "simplify" the first into the second.  The string `chython.periodictable` appearing inside a
subprocess payload is fine and this test is written to allow it -- it looks at import STATEMENTS,
not at text.
"""
import ast
import pathlib
import re


CORE = pathlib.Path(__file__).resolve().parent.parent
ROOT = CORE.parent

# everything in the distribution that the core may not name.  `chython.core` is the one allowed
# `chython.*` prefix; the bare facade is the worst of the lot and is checked for separately
FORBIDDEN = ('periodictable', 'containers', 'files', 'algorithms', 'reactor',
             'exceptions', 'reactor', 'core.deprecated')

# EMPTY, and it is meant to stay that way.  DO NOT ADD TO THIS LIST -- an entry is an instance of the
# blocker this test exists to have removed, and the second assertion below fails if one appears.
KNOWN_DEBT = set()


def sources():
    """Every file under `chython/core/`, tests included, that could carry an import."""
    out = []
    for path in sorted(CORE.rglob('*')):
        if path.suffix in ('.py', '.pyx', '.pxi', '.pxd') and path.is_file():
            out.append(path)
    return out


def python_imports(path):
    """Dotted module names imported by a `.py` file, from its AST -- not from its text.

    An AST and not a grep because a subprocess payload, a docstring and a comment all mention the
    forbidden names legitimately, and a grep would either fail on those or be taught exceptions
    until it stopped meaning anything.
    """
    # `encoding='utf-8'` here and in `cython_imports` below: `core/` is UTF-8 -- `_molecule_container.pxi`
    # writes `4.3×10⁹` about an int32 overflow -- and `read_text` without it asks the locale, which on the
    # Windows runner is cp1252 and raises, turning a scan that found nothing forbidden into an error.
    tree = ast.parse(path.read_text(encoding='utf-8'), filename=str(path))
    out = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            out.extend(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            if node.level:
                # a relative import: resolve it against this file's package so that `from ...`
                # reaching out of the core is caught rather than skipped
                package = path.relative_to(ROOT.parent).parent.parts
                if node.level - 1:
                    package = package[:-(node.level - 1)]
                base = '.'.join(package)
                if node.module:
                    out.append(f'{base}.{node.module}')
                else:
                    # `from .. import files` names its target in the alias, not in `module`, and
                    # that is the shape a refactor reaching sideways actually takes
                    out.extend(f'{base}.{alias.name}' for alias in node.names)
            elif node.module:
                out.append(node.module)
    return out


# a Cython `import`/`cimport` statement at the start of a line; `.pxi` files are textually included
# so they are not parseable on their own and there is no AST to ask
CYTHON_IMPORT = re.compile(r'^\s*(?:from\s+([\w.]+)\s+c?import\b|c?import\s+([\w.]+))', re.M)


def cython_imports(path):
    out = []
    for a, b in CYTHON_IMPORT.findall(path.read_text(encoding='utf-8')):
        out.append(a or b)
    return out


def offenders():
    """(file, module) for every import of this distribution that is not `chython.core`."""
    found = []
    for path in sources():
        reader = python_imports if path.suffix == '.py' else cython_imports
        for module in reader(path):
            head = module.split('.')
            if head[0] != 'chython':
                continue
            if len(head) == 1 or head[1] != 'core':
                found.append((path.name, module))
    return found


def test_the_core_imports_nothing_else_from_this_distribution():
    """The core's suite runs against `chython.core` and the standard library and nothing else.

    A failure names the file and the module.  The fix is never to add the file to `KNOWN_DEBT`: it
    is either to state the answer directly, or to reach chython 2 through `oracle.ask`, which
    runs it in another interpreter and is therefore not an import.
    """
    bad = [(f, m) for f, m in offenders() if f not in KNOWN_DEBT]
    assert not bad, (
        'chython/core/ must not import anything from this distribution except chython.core -- '
        'importing the facade or any chython 2 package loads ~130 modules and makes the core\'s '
        f'suite die with chython 2:\n' + '\n'.join(f'  {f}: {m}' for f, m in sorted(bad)))


def test_the_known_debt_has_not_grown():
    """The allow-list is a ratchet: it may shrink, and adding to it must be a visible decision."""
    assert {f for f, _ in offenders()} <= KNOWN_DEBT, \
        'a new file reached into chython 2; remove the import rather than widening KNOWN_DEBT'
    # and every entry still earns its place: a stale name would quietly re-open the hole
    for name in KNOWN_DEBT:
        assert (CORE / 'test' / name).is_file(), f'{name} is gone; drop it from KNOWN_DEBT'


def test_the_forbidden_names_are_the_packages_that_actually_exist():
    """The list is only a guard while at least one name it holds is a module that exists.

    A list matching nothing makes the guard above a no-op that passes because nothing can match it,
    so this test says so loudly instead.
    """
    present = {name for name in FORBIDDEN
               if (ROOT / name).is_dir() or (ROOT / f'{name}.py').is_file()}
    assert present, 'none of the forbidden packages exist any more: chython 2 is gone, so this ' \
                    'file and oracle.py have done their job and can be deleted together'


def test_the_scan_actually_reaches_the_files_it_claims_to():
    # the guard is evidence only if it is looking at something: pin the shape of the walk so that a
    # renamed suffix or a moved directory shows up as a failure here and not as a silent pass
    paths = sources()
    names = {p.name for p in paths}
    assert len(paths) > 60, len(paths)
    assert '_valence.pxi' in names and '_core.pyx' in names
    assert 'test_valence.py' in names and 'oracle.py' in names
    assert sum(1 for p in paths if p.suffix == '.pxi') > 15


def test_the_scan_would_catch_an_offender(tmp_path):
    """Plant one and check the reader sees it -- in both flavours of import, and relative too."""
    py = tmp_path / 'planted.py'
    py.write_text('from chython import smiles\n'
                  'import chython.periodictable\n'
                  'from chython.files.daylight.smiles import smiles as s\n'
                  'from chython.core._core import read_smiles\n'
                  'import json\n')
    found = python_imports(py)
    assert found == ['chython', 'chython.periodictable', 'chython.files.daylight.smiles',
                     'chython.core._core', 'json']

    pxi = tmp_path / 'planted.pxi'
    pxi.write_text('# import chython.periodictable in a comment is not an import\n'
                   'from chython.containers import MoleculeContainer\n'
                   'cimport chython.algorithms\n')
    assert cython_imports(pxi) == ['chython.containers', 'chython.algorithms']

    # a subprocess payload naming a forbidden module is NOT an import, and must not be flagged
    ok = tmp_path / 'subprocess_user.py'
    ok.write_text('from chython.core.test.oracle import ask\n'
                  'ask("from chython.periodictable.base.element import _elements_map")\n')
    assert python_imports(ok) == ['chython.core.test.oracle']


def test_a_relative_import_out_of_the_core_is_caught(tmp_path):
    # `from ..files import x` inside chython/core/test/ is the shape a refactor produces, and a
    # scan that only looked at absolute names would miss every one of them
    package = tmp_path / 'chython' / 'core' / 'test'
    package.mkdir(parents=True)
    planted = package / 'planted.py'
    planted.write_text('from . import gen_valence_rules\n'      # chython.core.test
                       'from .. import _core\n'                 # chython.core
                       'from ...files import SDFrw\n')          # chython.files -- the offender
    # `python_imports` resolves against the path, so build one that looks like the real tree
    global ROOT
    saved, ROOT = ROOT, tmp_path / 'chython'
    try:
        found = python_imports(planted)
    finally:
        ROOT = saved
    assert found == ['chython.core.test.gen_valence_rules', 'chython.core._core',
                     'chython.files']
    assert [m for m in found if m.split('.')[:2] != ['chython', 'core']] == ['chython.files']
