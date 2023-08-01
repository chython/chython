#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  Copyright 2014-2023 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from distutils.command.sdist import sdist
from distutils.command.build import build
from distutils.util import get_platform
from importlib.util import find_spec
from pathlib import Path
from setuptools import setup, Extension, find_packages


class _sdist(sdist):
    def finalize_options(self):
        super().finalize_options()
        self.distribution.data_files.append(('lib', ['INCHI/libinchi.so',
                                                     'INCHI/libinchi.dll', 'INCHI/libinchi.dynlib']))


cmd_class = {'sdist': _sdist}


if find_spec('wheel'):
    from wheel.bdist_wheel import bdist_wheel

    class _bdist_wheel(bdist_wheel):
        def finalize_options(self):
            super().finalize_options()
            self.root_is_pure = False
            platform = get_platform()
            if platform == 'win-amd64':
                self.distribution.data_files.append(('lib', ['INCHI/libinchi.dll']))
            elif platform == 'linux-x86_64':
                self.distribution.data_files.append(('lib', ['INCHI/libinchi.so']))
            elif platform == 'macosx-11-x86_64':
                self.distribution.data_files.append(('lib', ['INCHI/libinchi.dynlib']))

    cmd_class['bdist_wheel'] = _bdist_wheel


if find_spec('cython'):
    class _build(build):
        def finalize_options(self):
            super().finalize_options()
            from Cython.Build import cythonize
            self.distribution.ext_modules = cythonize(self.distribution.ext_modules, language_level=3)

    cmd_class['build'] = _build


setup(
    name='chython',
    version='1.66',
    packages=find_packages(),
    url='https://github.com/chython/chython',
    license='LGPLv3',
    author='Dr. Ramil Nugmanov',
    author_email='nougmanoff@protonmail.com',
    python_requires='>=3.8',
    cmdclass=cmd_class,
    ext_modules=[Extension('chython.algorithms._isomorphism', ['chython/algorithms/_isomorphism.pyx'],
                           extra_compile_args=['-O3']),
                 Extension('chython.containers._unpack', ['chython/containers/_unpack.pyx'],
                           extra_compile_args=['-O3']),
                 Extension('chython.containers._pack', ['chython/containers/_pack.pyx'], extra_compile_args=['-O3']),
                 Extension('chython.containers._cpack', ['chython/containers/_cpack.pyx'], extra_compile_args=['-O3']),
                 Extension('chython.files._xyz', ['chython/files/_xyz.pyx'], extra_compile_args=['-O3'])],
    setup_requires=['wheel', 'cython'],
    install_requires=['CachedMethods>=0.1.4,<0.2', 'lazy-object-proxy>=1.6', 'lxml>=4.1', 'py-mini-racer>=0.4.0',
                      'numpy>=1.18'],
    extras_require={'pytest': ['pytest'], 'mapping': ['chytorch-rxnmap>=1.4']},
    package_data={'chython.algorithms': ['_isomorphism.pyx'], 'chython.algorithms.calculate2d': ['clean2d.js'],
                  'chython.containers': ['_pack.pyx', '_unpack.pyx', '_cpack.pyx'], 'chython.files': ['_xyz.pyx']},
    data_files=[],
    zip_safe=False,
    long_description=(Path(__file__).parent / 'README.rst').read_text('utf8'),
    classifiers=['Environment :: Plugins',
                 'Intended Audience :: Science/Research',
                 'License :: OSI Approved :: GNU Lesser General Public License v3 or later (LGPLv3+)',
                 'Operating System :: OS Independent',
                 'Programming Language :: Python',
                 'Programming Language :: Python :: 3 :: Only',
                 'Programming Language :: Python :: 3.8',
                 'Topic :: Scientific/Engineering',
                 'Topic :: Scientific/Engineering :: Chemistry',
                 'Topic :: Scientific/Engineering :: Information Analysis',
                 'Topic :: Software Development',
                 'Topic :: Software Development :: Libraries',
                 'Topic :: Software Development :: Libraries :: Python Modules']
)
