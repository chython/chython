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
from Cython.Build import build_ext, cythonize
from pathlib import Path
from setuptools import Extension
from setuptools.dist import Distribution
from shutil import copyfile, which
from subprocess import run
from sysconfig import get_platform
from tempfile import TemporaryDirectory
from warnings import warn


platform = get_platform()
if platform == 'win-amd64':
    libname = 'libinchi.dll'
    extra_compile_args = ['/O2']
elif platform.startswith('linux'):
    libname = 'libinchi.so'
    extra_compile_args = ['-O3']
elif platform.startswith('macosx'):
    libname = 'libinchi.dylib'
    extra_compile_args = []
else:
    libname = None
    extra_compile_args = []


def build_libinchi(target: Path):
    """
    Compile the bundled IUPAC InChI source into a shared library at `target`.

    The InChI C source lives in `INCHI/INCHI-1-SRC/` and ships a CMake project
    for the `libinchi` target. We build it out-of-source and copy the produced
    binary (named `libinchi.{dll,so,dylib}` on the respective platforms) into
    the package. If the source tree or CMake is unavailable the build is skipped
    with a warning so that a wheel is still produced (InChI features disabled).
    """
    source = Path('INCHI/INCHI-1-SRC/INCHI-1-SRC/INCHI_API/libinchi/src')
    if not source.is_dir():
        warn(f'InChI source not found at {source}; skipping libinchi build')
        return
    if which('cmake') is None:
        warn('cmake not found; skipping libinchi build')
        return

    with TemporaryDirectory() as tmp:
        run(['cmake', '-S', str(source), '-B', tmp, '-DCMAKE_BUILD_TYPE=Release'], check=True)
        run(['cmake', '--build', tmp, '--config', 'Release', '--target', 'libinchi'], check=True)
        # CMake places the shared library under lib/ (Unix) or bin/ (Windows DLL);
        # multi-config generators add a Release/ level. rglob finds it at any depth.
        for produced in Path(tmp).rglob(libname):
            copyfile(produced, target)
            break
        else:
            warn(f'libinchi build produced no {libname}; skipping')


if libname:
    build_libinchi(Path('chython/files/libinchi') / libname)

extensions = [
    Extension('chython.algorithms._isomorphism',
              ['chython/algorithms/_isomorphism.pyx'],
              extra_compile_args=extra_compile_args),
    Extension('chython.containers._pack_v2',
              ['chython/containers/_pack_v2.pyx'],
              extra_compile_args=extra_compile_args),
    Extension('chython.containers._unpack_v0v2',
              ['chython/containers/_unpack_v0v2.pyx'],
              extra_compile_args=extra_compile_args),
    Extension('chython.files._xyz',
              ['chython/files/_xyz.pyx'],
              extra_compile_args=extra_compile_args),
    Extension('chython.algorithms._rings',
              ['chython/algorithms/_rings.pyx'],
              extra_compile_args=extra_compile_args)
]

ext_modules = cythonize(extensions, language_level=3,
                        compiler_directives={'freethreading_compatible': True})
cmd = build_ext(Distribution({'ext_modules': ext_modules}))
cmd.ensure_finalized()
cmd.run()

for output in cmd.get_outputs():
    output = Path(output)
    copyfile(output, output.relative_to(cmd.build_lib))
