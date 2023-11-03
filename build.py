# -*- coding: utf-8 -*-
#
#  Copyright 2023 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from shutil import copyfile
from sysconfig import get_platform


platform = get_platform()
if platform == 'win-amd64':
    libname = 'libinchi.dll'
    extra_compile_args = ['/O2']
elif platform == 'linux-x86_64':
    libname = 'libinchi.so'
    extra_compile_args = ['-O3']
elif platform.startswith('macosx') and platform.endswith('x86_64'):
    libname = 'libinchi.dynlib'
    extra_compile_args = []
else:
    libname = None
    extra_compile_args = []

if libname:
    copyfile(Path('INCHI') / libname, Path('chython/files/libinchi') / libname)

extensions = [
    Extension('chython.algorithms._isomorphism',
              ['chython/algorithms/_isomorphism.pyx'],
              extra_compile_args=extra_compile_args),
    Extension('chython.containers._pack',
              ['chython/containers/_pack.pyx'],
              extra_compile_args=extra_compile_args),
    Extension('chython.containers._unpack',
              ['chython/containers/_unpack.pyx'],
              extra_compile_args=extra_compile_args),
    Extension('chython.containers._cpack',
              ['chython/containers/_cpack.pyx'],
              extra_compile_args=extra_compile_args),
    Extension('chython.files._xyz',
              ['chython/files/_xyz.pyx'],
              extra_compile_args=extra_compile_args)
]

ext_modules = cythonize(extensions, language_level=3)
cmd = build_ext(Distribution({'ext_modules': ext_modules}))
cmd.ensure_finalized()
cmd.run()

for output in cmd.get_outputs():
    output = Path(output)
    copyfile(output, output.relative_to(cmd.build_lib))
