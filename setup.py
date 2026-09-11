# -*- coding: utf-8 -*-
#
#  Copyright 2023-2026 Ramil Nugmanov <nougmanoff@protonmail.com>
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
"""The compiled half of the build.  Metadata lives in `pyproject.toml`.

This file exists because the extension list is COMPUTED -- the platform's optimisation flags and
Cython's directives -- and `pyproject.toml` cannot express a computation.  What is deprecated is
`python setup.py <command>` as a way to build or publish; supplying `ext_modules` to a PEP 517
backend is supported.

Nothing is copied back into the source tree and `pyproject.toml` globs no build output: the backend
places what it built, so a wheel's contents cannot become a function of what lies in the checkout.

DEVELOPER LOOP: `python setup.py build_ext --inplace`, the one place the legacy command line is still
the answer -- an in-place build is not a distribution, so PEP 517 has no equivalent.
"""
import sys
from pathlib import Path
from shutil import copyfile
from sysconfig import get_platform

from Cython.Build import cythonize
from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext
from setuptools.command.build_py import build_py


# The two ways this file runs differ only here.  `python setup.py build_ext --inplace` runs it as a
# script, so `__file__` exists and its directory is `sys.path[0]`.  A PEP 517 backend `exec`s the
# source with neither, and `import build_inchi` then raises `ModuleNotFoundError` though the file is
# right there.  The working directory is the source tree either way, which makes the fallback sound.
ROOT = Path(__file__).resolve().parent if '__file__' in globals() else Path.cwd().resolve()
sys.path.insert(0, str(ROOT))

from build_inchi import build as build_libinchi, libname_for_platform  # noqa: E402  -- needs ROOT


# `-fsigned-char`, and not only on aarch64: plain `char` is UNSIGNED on Linux ARM and signed
# everywhere else this builds -- x86-64 Linux, both macOS arches, MSVC -- so without it one wheel in
# the matrix reads a negative `char` as its 256-complement and nothing in CI would show it.  The flag
# is what makes the arm wheel the same program as the one the tests ran against; it is a no-op on
# x86-64, which is why it is not conditional within Linux.  MSVC needs nothing (its default is signed;
# `/J` is the opposite switch), and the InChI build gets the same flag from `build_inchi.py`.
#
# `-g0` ON LINUX ONLY, because distutils compiles an extension with CPython's own `CFLAGS` and those
# carry `-g`: measured on 3.0's cp312 wheel, `.debug*` was 17.50 MB of a 20.82 MB `_core.so` whose
# `.text` is 2.52 MB.  These flags are appended after CPython's and gcc takes the last of `-g`/`-g0`,
# which is what makes a flag here able to cancel one from there.  `.symtab` is not debug info and
# survives, so a C-level backtrace keeps its function names.  Nothing to do on the other two: clang
# leaves DWARF in the `.o` files and emits a `.dSYM` only when asked, and MSVC writes a separate `.pdb`.
platform = get_platform()
if platform == 'win-amd64':
    extra_compile_args = ['/O2']
elif platform.startswith('linux'):
    extra_compile_args = ['-O3', '-g0', '-fsigned-char']
else:
    extra_compile_args = ['-fsigned-char']

# libinchi is a BUILD OUTPUT: built into `build/inchi/`, put in the wheel by `StageLibinchi` below, and
# written into the source tree only by an in-place build.
#
# Unconditional, because `build()` returns early when the binary is already there -- free on every
# later build, and on every CI wheel job, which downloads one binary per OS instead of building it --
# and warns rather than raising when cmake or the InChI sources are absent, so an sdist without the
# submodule still yields a working wheel minus InChI.
LIBINCHI_DIR = ROOT / 'build' / 'inchi'
libname = libname_for_platform()
if libname:
    build_libinchi(LIBINCHI_DIR / libname)

# ONE EXTENSION, ratcheted by `chython/test/test_v2_boundary.py`: a second one needs a reason written
# there.  `prune_stale_staging()` below is what clears the artefact an extension leaves in `build/lib*`
# once it stops being built.
extensions = [
    # the whole core is one translation unit: _core.pyx includes the .pxi layers, so every
    # cdef call between them is a static call the C compiler can inline
    Extension('chython.core._core',
              ['chython/core/_core.pyx'],
              extra_compile_args=extra_compile_args),
]


def prune_stale_staging():
    """Delete anything staged under `build/lib*` that the source tree no longer has.

    The staging directory is never emptied, so `bdist_wheel` packages what previous builds left there.
    A deleted `.pyx` keeps shipping its `.so`, and a renamed package ships under both names: after
    `files/ctfile/` became `formats/ctfile/`, `import chython.files.ctfile` worked out of a wheel and
    served pre-move code no test covers.  Nothing fails -- the file is simply in the directory being
    zipped.  Pure Python strands the more easily of the two, moving a package being cheap and not
    looking like touching the build.

    Three authorities, one per failure mode:

    * an extension is claimed by the `Extension` list above and nothing else -- its staged name
      (`_core.cpython-310-darwin.so`) matches no filename in the source, so it cannot be reconciled
      against the tree.  Only ABI-tagged names are considered, which keeps `core/libinchi.so` -- a
      shared library on Linux, legitimately not an extension -- out of this rule without a path
      exclusion that could drift;
    * `libinchi.*` is claimed by `StageLibinchi` below, the one staged file with no source counterpart
      by design.  Named explicitly rather than left to the ordering that would re-stage it anyway,
      since a rule that is only harmless by accident is not one to rely on;
    * everything else is claimed by a counterpart at the same path in the source tree -- one rule for
      a deleted module, a moved package and a renamed data file alike.
    """
    wanted = {ext.name for ext in extensions}
    for lib in ROOT.glob('build/lib*'):
        if lib == LIBINCHI_DIR:  # cannot happen while the directory is named `inchi`, and is fatal if it does
            continue
        for staged in sorted(p for p in lib.rglob('*') if p.is_file()):
            rel = staged.relative_to(lib)
            if staged.name.startswith('libinchi.'):
                continue  # staged by StageLibinchi from build/inchi/; see the docstring
            if staged.name.endswith(('.pyd',)) or '.cpython-' in staged.name:
                module = '.'.join((*rel.parts[:-1], rel.name.split('.')[0]))
                if module not in wanted:
                    staged.unlink()
                    print('pruned stale extension %s (no longer built from source)' % rel)
            elif not (ROOT / rel).exists():
                staged.unlink()
                print('pruned stale %s (no longer in the source tree)' % rel)

        # directories the pruning emptied would otherwise stay as empty packages in the wheel
        for directory in sorted((p for p in lib.rglob('*') if p.is_dir()), reverse=True):
            if not any(directory.iterdir()):
                directory.rmdir()
                print('pruned emptied directory %s' % directory.relative_to(lib))


def stage_libinchi(into: Path) -> bool:
    """Copy the built libinchi into *into* -- a `chython/core` directory, staged or in the source tree.

    Returns False when there is nothing to copy: an unsupported platform, or no cmake, or no INCHI
    submodule.  A supported outcome and not a failure -- the build is then one without InChI support,
    which `core/__init__.py`'s silent fallback and the `needs_inchi` skips are for.
    """
    if not libname:
        return False
    source = LIBINCHI_DIR / libname
    if not source.exists():
        print('libinchi absent from %s; the build will have no InChI support' % LIBINCHI_DIR)
        return False
    into.mkdir(parents=True, exist_ok=True)
    copyfile(source, into / libname)
    # printed as given, neither made relative to ROOT: `build_py` hands over a RELATIVE `build_lib`, and
    # `Path.relative_to` raises rather than coping when one side is relative and the other absolute
    print('staged %s -> %s' % (source, into / libname))
    return True


class StageLibinchi(build_py):
    """Put libinchi in the wheel.

    Not `package_data`, which names files inside the SOURCE tree -- declaring `libinchi.*` there is what
    required the binary to sit in `chython/core/`.  Naming it here means a stale binary left in a
    checkout by an in-place build is no longer packageable: the wheel gets this build's, or none.
    """
    def run(self):
        super().run()
        stage_libinchi(Path(self.build_lib) / 'chython' / 'core')


class BuildExtInplaceStagesLibinchi(build_ext):
    """The developer loop, and the editable install, which `build_py` above does not serve.

    `build_ext --inplace` and `pip install -e .` both write the extension into the source tree, where
    `core/__init__.py` looks for libinchi beside it.  Without this both silently lose InChI: roughly a
    hundred tests skip and `test_inchi.py` reports success having asserted nothing.

    So the binary does land in `chython/core/` for an in-place build -- the exception the module
    docstring draws -- and only that way, which is what stops a distribution build picking it up.
    """
    def run(self):
        super().run()
        if self.inplace:
            stage_libinchi(ROOT / 'chython' / 'core')


prune_stale_staging()

setup(ext_modules=cythonize(extensions, language_level=3,
                            compiler_directives={'freethreading_compatible': True}),
      cmdclass={'build_py': StageLibinchi, 'build_ext': BuildExtInplaceStagesLibinchi})
