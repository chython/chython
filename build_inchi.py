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
Build libinchi from the bundled INCHI submodule.

Run standalone:
    python build_inchi.py [--target <path>]

If --target is not given the binary is written to `build/inchi/`, which is a BUILD OUTPUT DIRECTORY
and not part of the source tree.  `setup.py` stages it from there into the wheel; an in-place build
additionally copies it next to the extension, because `core/__init__.py` loads it relative to its own
`__file__` and `core` is its only consumer (the bridge is `_inchi.pxi`, inside that same extension).

WHY NOT chython/core/ DIRECTLY, which is where this used to write.  A build that writes a binary into
the source tree makes the wheel's contents a function of what happens to be lying in the checkout --
the same failure `setup.py`'s docstring describes for the extensions and the copy-back build it
replaced.  Writing to `build/` costs one staging step and makes the source tree's cleanliness a
property of the build rather than of `.gitignore`.

`build/inchi/` and NOT `build/libinchi/`: `prune_stale_staging()` in `setup.py` treats every
`build/lib*` directory as a wheel staging area, and `build/libinchi` matches that glob -- it would
delete this binary as a file with no counterpart in the source tree, which is exactly what it is.
"""
import argparse
import sys
from pathlib import Path
from shutil import copyfile, which
from subprocess import run
from sysconfig import get_platform
from tempfile import TemporaryDirectory
from warnings import warn


# Resolved from this file and not from the working directory.  `setup.py` imports this module, and a
# build backend is entitled to run from anywhere; a relative `INCHI/...` silently degrades into the
# "source not found, skipping" warning path, which produces a wheel with no InChI and no error.
ROOT = Path(__file__).resolve().parent


def cmake_args() -> list[str]:
    """Everything the bundled ``CMakeLists.txt`` does not decide and a wheel cannot be left to guess.

    `core/__init__.py` loads this binary with `ctypes` INTO THE RUNNING PROCESS, so its architectures,
    its minimum OS version and its `char` signedness are the interpreter's and the extension's, not
    the build machine's defaults.

    * ``-fsigned-char`` -- plain `char` is unsigned on Linux ARM.  `setup.py` compiles the extension
      with the flag for the same reason, and InChI's own gcc option set lists it; its public types are
      `S_CHAR`, spelled `signed char`, so the flag is about the library's internal bare `char`.
    * ``CMAKE_SHARED_LINKER_FLAGS=-Wl,-s`` on Linux -- ``INCHI_API/libinchi/src/CMakeLists.txt`` gives
      gcc-like compilers ``-g;-O1`` through ``target_compile_options``, which lands after
      ``CMAKE_C_FLAGS`` and after the ``Release`` config's own flags, so neither ``CMAKE_BUILD_TYPE``
      nor a ``-g0`` above can cancel it.  A link-time strip can: measured on 3.0's Linux wheel,
      ``.debug*`` was 3.05 MB of a 4.36 MB ``libinchi.so`` whose ``.text`` is 1.00 MB.  Linux only --
      ld64 deprecates ``-s`` and the linked Mach-O carries no DWARF to begin with.
    * ``CMAKE_OSX_ARCHITECTURES`` -- a universal2 interpreter needs both slices, or InChI is absent on
      the arch the dylib lacks while `import chython` still succeeds and every InChI test skips.
    * ``CMAKE_OSX_DEPLOYMENT_TARGET`` -- ``wheel``'s ``calculate_macosx_platform_tag`` raises the
      wheel's tag to cover every binary inside it, so a dylib stamped with the build machine's macOS
      version tags the whole wheel with that version: built on macOS 26 the wheel installed on macOS
      26 and nowhere earlier, whatever the extension had been compiled for.  With the flag the tag is
      `sysconfig`'s own -- `macosx_10_9_universal2` -- and cmake raises the arm64 slice to 11.0 itself.

    `sysconfig.get_platform()` is the source throughout -- `linux-aarch64`, `macosx-10.9-universal2`,
    `macosx-11.0-arm64` -- being what the interpreter says about itself.
    """
    platform = get_platform()
    if platform.startswith('win'):
        return []                                       # MSVC's `char` is signed and the rest is mac
    args = ['-DCMAKE_C_FLAGS=-fsigned-char']
    if platform.startswith('linux'):
        args.append('-DCMAKE_SHARED_LINKER_FLAGS=-Wl,-s')
    parts = platform.split('-')
    if platform.startswith('macosx') and len(parts) == 3:
        _, target, arch = parts
        archs = 'arm64;x86_64' if arch in ('universal2', 'fat64', 'intel') else arch
        args += [f'-DCMAKE_OSX_DEPLOYMENT_TARGET={target}', f'-DCMAKE_OSX_ARCHITECTURES={archs}']
    return args


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

    source = ROOT / 'INCHI/INCHI-1-SRC/INCHI_API/libinchi/src'
    if not source.is_dir():
        warn(f'InChI source not found at {source}; skipping libinchi build')
        return False
    if which('cmake') is None:
        warn('cmake not found; skipping libinchi build')
        return False

    extra = cmake_args()

    target.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory() as tmp:
        run(['cmake', '-S', str(source), '-B', tmp, '-DCMAKE_BUILD_TYPE=Release', *extra], check=True)
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
                        default=ROOT / 'build/inchi' / (libname_for_platform() or 'libinchi'),
                        help='Destination path for the compiled binary')
    args = parser.parse_args()
    sys.exit(0 if build(args.target) else 1)
