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
"""
Build libinchi from the bundled INCHI submodule and copy it into
chython/files/libinchi/.

Run standalone:
    python build_inchi.py [--target <path>]

If --target is not given the binary is written to chython/files/libinchi/.
"""
import argparse
import sys
from pathlib import Path
from shutil import copyfile, which
from subprocess import run
from sysconfig import get_platform
from tempfile import TemporaryDirectory
from warnings import warn


def libname_for_platform() -> str | None:
    p = get_platform()
    if p == 'win-amd64':
        return 'libinchi.dll'
    if p.startswith('linux'):
        return 'libinchi.so'
    if p.startswith('macosx'):
        return 'libinchi.dylib'
    return None


def build(target: Path) -> bool:
    """
    Build libinchi and write the binary to *target*.
    Returns True on success, False if skipped (source or cmake missing).
    """
    if target.exists():
        return True

    libname = libname_for_platform()
    if libname is None:
        warn(f'Unsupported platform {get_platform()}; skipping libinchi build')
        return False

    source = Path('INCHI/INCHI-1-SRC/INCHI_API/libinchi/src')
    if not source.is_dir():
        warn(f'InChI source not found at {source}; skipping libinchi build')
        return False
    if which('cmake') is None:
        warn('cmake not found; skipping libinchi build')
        return False

    target.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory() as tmp:
        run(['cmake', '-S', str(source), '-B', tmp, '-DCMAKE_BUILD_TYPE=Release'], check=True)
        run(['cmake', '--build', tmp, '--config', 'Release', '--target', 'libinchi'], check=True)
        for produced in Path(tmp).rglob(libname):
            copyfile(produced, target)
            print(f'libinchi written to {target}')
            return True
        warn(f'libinchi build produced no {libname}; skipping')
        return False


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Build libinchi shared library')
    parser.add_argument('--target', type=Path,
                        default=Path('chython/files/libinchi') / (libname_for_platform() or 'libinchi'),
                        help='Destination path for the compiled binary')
    args = parser.parse_args()
    sys.exit(0 if build(args.target) else 1)
