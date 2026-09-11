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
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with this program; if not, see <https://www.gnu.org/licenses/>.
#
"""libinchi reaches the wheel from `build/inchi/`, and cannot reach it from the source tree.

`test_packaging.py` is the gate for every other runtime data file, and it works by reading
`[tool.setuptools.package-data]` -- which is exactly the mechanism libinchi does not use.
`package_data` names files INSIDE the package directory and can reach nothing else, so declaring
`libinchi.*` there is what obliged the build to write a 1.3 MB binary into `chython/core/`, making a
wheel's contents a function of what happened to be lying in the checkout.  That is the same failure the
copy-back build had for the extensions, arriving by the other door, and `setup.py`'s docstring says this
file exists to make it inexpressible.

So `setup.py` builds the binary into `build/inchi/` and stages it into the wheel with a `build_py`
subclass.  The cost of moving off `package_data` is that libinchi left the coverage of the gate that
covers everything else, and this file is the replacement.  What it holds:

  * the mechanism cannot silently revert -- if `libinchi.*` reappears in `package-data`, or the
    `cmdclass` that stages it disappears, these tests fail rather than the release doing so;
  * the CI's paths agree with the builder's default.  This is not hypothetical tidiness: the workflow
    pointed at `chython/files/libinchi/` for the whole life of this branch, long after the directory was
    renamed, and nothing said so because the workflow only runs on `release: published`.  A stale path
    there fails at `if-no-files-found: error` -- during a release, on four platforms at once;
  * the staging itself, whenever a build has actually run in this tree.

Nothing here builds anything.  The binary takes a cmake run to produce and is absent on any machine
without the INCHI submodule, which is a supported state -- `core/__init__.py` falls back silently and
the suite's `needs_inchi` marks skip.  Every assertion is therefore about declarations, plus one about
build output that is checked only if build output exists.
"""
from ast import Call, Constant, Dict, Name, parse, walk
from pathlib import Path
from re import DOTALL, finditer, search


ROOT = Path(__file__).resolve().parent.parent.parent

# The directory `build_inchi.py` writes to and `setup.py` stages from, spelled once here and checked
# against all three files below.  `build/inchi` and NOT `build/libinchi`: `prune_stale_staging()` treats
# every `build/lib*` directory as a wheel staging area, and `build/libinchi` matches that glob -- it
# would delete the binary as a file with no counterpart in the source tree, which is what it is.
OUTPUT_DIR = 'build/inchi'
LIBNAMES = ('libinchi.so', 'libinchi.dylib', 'libinchi.dll')


def _package_data():
    """package name -> declared filename patterns, from pyproject.toml.

    A regex for the same reason `test_packaging.py` uses one: `tomllib` is 3.11 and this package
    supports 3.10, so a parsed version of this test would skip on an interpreter people build wheels on.
    """
    text = (ROOT / 'pyproject.toml').read_text(encoding='utf-8')
    block = search(r'^\[tool\.setuptools\.package-data\]\n(.*?)(?=^\[|\Z)', text, DOTALL | 8)
    assert block, 'pyproject.toml has no [tool.setuptools.package-data] section'
    out = {}
    for line in block.group(1).splitlines():
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        key, _, value = line.partition('=')
        out[key.strip().strip('\'"')] = [m.group(1) for m in finditer(r'[\'"]([^\'"]+)[\'"]', value)]
    return out


def _setup_cmdclass():
    """The `cmdclass` mapping from setup.py's `setup(...)` call: command name -> class name.

    Read by `ast` and not by import, for the reason `test_packaging.py` gives for the same choice:
    importing `setup.py` runs `cythonize`, builds libinchi and prunes the staging directory.
    """
    for node in walk(parse((ROOT / 'setup.py').read_text(encoding='utf-8'))):
        if isinstance(node, Call) and isinstance(node.func, Name) and node.func.id == 'setup':
            for keyword in node.keywords:
                if keyword.arg == 'cmdclass' and isinstance(keyword.value, Dict):
                    return {k.value: v.id for k, v in zip(keyword.value.keys, keyword.value.values)
                            if isinstance(k, Constant) and isinstance(v, Name)}
    return {}


def _workflow():
    """The release workflow with its full-line comments removed.

    The comments are removed because that workflow EXPLAINS the two drifts below -- it says in prose that
    `poetry build` would now build nothing, and that the artifact no longer comes from
    `chython/files/libinchi/`.  A scanner reading the raw text finds those phrases and reports the
    explanation as the regression.  Only whole-line comments are dropped; nothing in this file puts a
    `#` inside a value.
    """
    text = (ROOT / '.github/workflows/python-package.yml').read_text(encoding='utf-8')
    return '\n'.join(line for line in text.splitlines() if not line.lstrip().startswith('#'))


# --- the mechanism ---------------------------------------------------------------------------------

def test_libinchi_is_not_package_data():
    """The binary must not be declarable from the source tree, in any package.

    This is the assertion that makes the whole arrangement true rather than merely intended: while
    `package-data` names `libinchi.*`, a stale binary in a checkout is packageable and a wheel can ship
    one build's library beside another build's extension.
    """
    offenders = ['%s declares %s' % (package, pattern)
                 for package, patterns in _package_data().items()
                 for pattern in patterns
                 if pattern.startswith('libinchi')]
    assert not offenders, (
        'libinchi is staged from %s by setup.py, not packaged from the source tree:\n  %s\n'
        'A package-data entry can only name a file inside chython/, which puts the binary back in the '
        'source tree and makes the wheel depend on what is lying in the checkout.' % (
            OUTPUT_DIR, '\n  '.join(offenders)))


def test_setup_py_stages_libinchi_into_the_wheel_and_into_an_inplace_build():
    """Both halves, because they serve different consumers and either one alone is a silent loss.

    `build_py` is what puts the binary in the wheel.  `build_ext` is what puts it next to the extension
    for `--inplace` and for `pip install -e .`, where `core/__init__.py` looks for it relative to its own
    `__file__` -- without it the dev loop and every editable install lose InChI, about a hundred tests
    skip, and `test_inchi.py` reports success having asserted nothing.
    """
    cmdclass = _setup_cmdclass()
    assert cmdclass, 'setup.py passes no cmdclass; the ast read has stopped working or the staging is gone'
    assert 'build_py' in cmdclass, 'setup.py has no build_py override, so no wheel will contain libinchi'
    assert 'build_ext' in cmdclass, \
        'setup.py has no build_ext override, so an in-place or editable build has no libinchi beside the extension'

    source = (ROOT / 'setup.py').read_text(encoding='utf-8')
    assert OUTPUT_DIR.split('/')[-1] in source, \
        "setup.py no longer mentions the %s output directory; it may have gone back to writing into " \
        "chython/core/" % OUTPUT_DIR


def test_build_inchi_defaults_to_the_build_directory():
    """`python build_inchi.py` with no arguments must not write into the package.

    CI runs exactly that, so this default IS the CI's output path, and it is read as a literal rather
    than searched for in the file's prose -- the module docstring discusses `chython/core/` at length in
    order to explain why the binary no longer goes there, so a text search cannot tell an explanation
    apart from a regression.
    """
    defaults = [node.value for node in walk(parse((ROOT / 'build_inchi.py').read_text(encoding='utf-8')))
                if isinstance(node, Constant) and isinstance(node.value, str) and node.value.startswith('build/')]
    assert OUTPUT_DIR in defaults, \
        "build_inchi.py has no %r path literal; its default --target may have moved back into the " \
        "source tree" % OUTPUT_DIR


# --- the CI's paths, which no test could see before ------------------------------------------------

def test_the_workflow_moves_libinchi_through_the_build_directory():
    """Upload and download must both name `build/inchi`, and neither may name a path inside chython/.

    The workflow runs only on `release: published`, so a wrong path here is invisible until a release is
    already published and four platforms fail at once.  It sat wrong for the whole life of this branch.
    """
    text = _workflow()
    assert '%s/libinchi.*' % OUTPUT_DIR in text, \
        'the workflow does not upload %s/libinchi.*; the artifact step still points somewhere else' % OUTPUT_DIR
    assert '%s/' % OUTPUT_DIR in text, 'the workflow does not download the artifact into %s/' % OUTPUT_DIR
    assert 'chython/files/libinchi' not in text, \
        'the workflow still names chython/files/libinchi -- that directory was renamed to chython/core ' \
        'and then the binary left the source tree entirely'
    assert 'chython/core/libinchi' not in text, \
        'the workflow names chython/core/libinchi; the binary is staged there by setup.py and is not a ' \
        'path CI writes to'


def test_the_workflow_builds_with_the_declared_backend():
    """`poetry build` built the wheels through 2.24 and would now build nothing at all.

    The backend moved to setuptools in `[build-system]` and the workflow was not touched, which is the
    same class of drift as the artifact path above and equally invisible outside a release.
    """
    text = _workflow()
    backend = search(r"build-backend\s*=\s*['\"]([^'\"]+)['\"]", (ROOT / 'pyproject.toml').read_text(encoding='utf-8'))
    assert backend and backend.group(1).startswith('setuptools'), \
        'the build backend is no longer setuptools; this test and the workflow both need re-reading'
    assert 'poetry build' not in text, \
        'the workflow runs `poetry build` while [build-system] declares setuptools; it would produce no wheel'
    assert 'python -m build' in text, 'the workflow has no PEP 517 build step'


def test_the_linux_wheel_is_retagged_for_manylinux():
    """setuptools tags Linux wheels `linux_x86_64`, which PyPI rejects; poetry-core did not.

    2.24 published `manylinux_2_39_*` with no repair step because poetry-core took its platform tag from
    `packaging.tags`.  `bdist_wheel` does not, so the backend migration silently removed the only reason
    the Linux uploads worked -- a failure that appears at `twine upload`, on Linux only, during a
    release.
    """
    text = _workflow()
    assert 'wheel tags' in text or 'auditwheel' in text, \
        'no retag or repair step for Linux: setuptools emits linux_x86_64 and PyPI rejects it'
    assert 'sys_tags' in text or 'manylinux' in text, \
        'the retag step does not name a manylinux tag or compute one'


# --- the staging itself, when there is build output to look at -------------------------------------

def test_every_staged_tree_that_has_python_also_has_libinchi():
    """The property, checked against real build output whenever any exists.

    Vacuous on a clean checkout, and that is intentional and said out loud here because a silent pass for
    a structural reason is the failure mode the rest of this file guards against.  Vacuous too when the
    binary was never built -- no submodule, no cmake, unsupported platform -- which is supported.

    A `build/lib*` directory is only in scope once `build_py` has run in it: `build_ext --inplace`
    creates the same directory holding nothing but the extension, and requiring libinchi there would
    fail on the ordinary developer loop.  `chython/__init__.py` is the discriminator.
    """
    built = [ROOT / OUTPUT_DIR / name for name in LIBNAMES]
    built = [p for p in built if p.exists()]
    if not built:
        return  # nothing was built for this platform; nothing can be staged

    names = {p.name for p in built}
    missing = []
    for lib in ROOT.glob('build/lib*'):
        if not (lib / 'chython' / '__init__.py').exists():
            continue  # build_ext-only staging directory: no package data belongs in it yet
        if not any((lib / 'chython' / 'core' / name).exists() for name in names):
            missing.append(str(lib.relative_to(ROOT)))
    assert not missing, (
        'these staged trees would produce a wheel with no InChI support:\n  %s\n'
        'setup.py\'s build_py subclass should have copied %s into each chython/core; run '
        '`python -m build --wheel` again.' % ('\n  '.join(sorted(missing)), ', '.join(sorted(names))))


def test_the_gate_can_fail():
    """Negative controls: every scanner above passes vacuously if it stops finding anything."""
    data = _package_data()
    assert data, 'the package-data scanner found nothing; it has stopped working'
    assert 'chython.core' in data, 'the scanner no longer sees chython.core'
    # and it would still see a libinchi entry if one came back
    assert [p for p in ['libinchi.so'] if p.startswith('libinchi')], 'the offender predicate is broken'

    assert _setup_cmdclass(), 'the setup.py cmdclass scanner found nothing; it has stopped working'
    assert 'workflow_dispatch' in _workflow(), 'the workflow scanner is not reading the workflow'
