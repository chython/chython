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
"""Every wheel in the matrix is the same program, and the release ships a source distribution.

`test_libinchi_staging.py` holds where the InChI binary comes from and that it reaches the wheel.  This
file holds the two properties that are invisible in a single build and only fail on one row of the
release matrix:

  * **`char` signedness.**  Plain `char` is UNSIGNED on Linux ARM and signed on x86-64 Linux, both macOS
    arches and MSVC.  The core reads `signed char` fields out of InChI's structs, so without
    `-fsigned-char` the aarch64 wheel reads a chloride's -1 charge as 255 -- in a wheel nothing in CI
    runs, since the tests run on the arch that was already right.
  * **the macOS floor and slices.**  `build_inchi.py` derives `CMAKE_OSX_DEPLOYMENT_TARGET` and
    `CMAKE_OSX_ARCHITECTURES` from `sysconfig.get_platform()`, because `wheel`'s
    `calculate_macosx_platform_tag` RAISES a wheel's tag to cover every binary inside it: a dylib
    stamped with the build machine's macOS version tags the whole wheel with that version, and a dylib
    with one slice makes InChI absent on the other arch while `import chython` still succeeds and every
    InChI test skips.

Plus the sdist, which is neither -- it is what `pip install chython` falls back to for every platform
the matrix does not cover.

`cmake_args()` is exercised rather than scanned: it is a pure function of `sysconfig.get_platform()` and
the flags are what the assertion is about.  The workflow is read as text, by the regex
`test_libinchi_staging.py` explains -- PyYAML is not in `[dependency-groups] dev`, and a test that skips
without it does not run.  Nothing here builds anything.
"""
from ast import Assign, Constant, List, Name, parse, walk
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
from re import M, finditer
from struct import unpack_from


ROOT = Path(__file__).resolve().parent.parent.parent

# The one flag that makes a negative `char` negative on every platform this ships to.  gcc and clang
# both take it; MSVC's default is already signed and its switch (`/J`) is the opposite one.
SIGNED_CHAR = '-fsigned-char'

# The two that keep DWARF out of the Linux wheel: `-g0` cancels the `-g` distutils inherits from
# CPython's `CFLAGS`, and the link-time strip is the only thing that reaches the `-g` InChI's own
# CMakeLists gives gcc.  Neither has a macOS or MSVC counterpart -- those toolchains put debug info in a
# `.dSYM` and a `.pdb`, outside the wheel either way.
NO_DEBUG_INFO = '-g0'
STRIP = '-Wl,-s'


def _build_inchi():
    """The builder module, loaded from the repo root by path.

    Not importable by name -- it is a build script beside `setup.py` and not part of the package -- and
    safe to load: it defines two functions and a path, and its command line is behind `__main__`.
    """
    spec = spec_from_file_location('_build_inchi_under_test', ROOT / 'build_inchi.py')
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _extra_compile_args():
    """Every list assigned to `extra_compile_args` in setup.py, by `ast`.

    One per platform branch, and read rather than imported for the reason `test_packaging.py` gives:
    importing `setup.py` runs `cythonize`, builds libinchi and prunes the staging directory.
    """
    out = []
    for node in walk(parse((ROOT / 'setup.py').read_text(encoding='utf-8'))):
        if isinstance(node, Assign) and isinstance(node.value, List) \
                and any(isinstance(t, Name) and t.id == 'extra_compile_args' for t in node.targets):
            out.append([e.value for e in node.value.elts if isinstance(e, Constant)])
    return out


def _jobs():
    """job name -> its block, from the release workflow."""
    text = (ROOT / '.github/workflows/python-package.yml').read_text(encoding='utf-8')
    body = text.split('\njobs:\n', 1)[1]
    bounds = [(m.start(), m.group(1)) for m in finditer(r'^  ([A-Za-z_][\w-]*):$', body, M)]
    bounds.append((len(body), None))
    return {name: body[start:bounds[i + 1][0]] for i, (start, name) in enumerate(bounds[:-1]) if name}


# --- char signedness -------------------------------------------------------------------------------

def test_the_extension_is_compiled_with_signed_char_wherever_the_compiler_is_not_msvc():
    """Both non-Windows branches, and not the aarch64 one alone.

    The flag is a no-op on x86-64 and on macOS, whose `char` is signed anyway; applying it there is what
    makes the arm wheel the same program as the one the tests ran against, rather than a build that
    differs from every other row by one flag nobody reads.
    """
    branches = _extra_compile_args()
    assert len(branches) >= 2, 'setup.py no longer assigns extra_compile_args per platform; the ast read is stale'
    offenders = [flags for flags in branches
                 if not any(f.startswith('/') for f in flags) and SIGNED_CHAR not in flags]
    assert not offenders, (
        'these setup.py compile-flag branches do not force signed `char`: %s\n'
        'Plain `char` is unsigned on Linux ARM, so the aarch64 wheel reads every negative `signed char` '
        'from InChI as its 256-complement -- a chloride charge of -1 becomes 255.  Nothing in CI shows '
        'it: the tests run on the arch that was already right.' % offenders)


def test_libinchi_is_compiled_with_signed_char_on_every_platform_that_takes_the_flag():
    """The library and the extension must agree, being one process: `core/__init__.py` loads the dylib
    with `ctypes`, so a bare `char` inside InChI that the extension reads back as signed is one program
    compiled two ways."""
    module = _build_inchi()
    for platform in ('linux-aarch64', 'linux-x86_64', 'macosx-10.9-universal2', 'macosx-11.0-arm64'):
        module.get_platform = lambda p=platform: p
        flags = ' '.join(module.cmake_args())
        assert SIGNED_CHAR in flags, \
            'build_inchi.py passes no %s on %s, so libinchi and the extension disagree about the sign ' \
            'of a bare `char`' % (SIGNED_CHAR, platform)


def test_the_windows_build_asks_cmake_for_nothing():
    """MSVC's `char` is signed and the deployment target and architectures are macOS concepts, so the
    Windows branch has nothing to add -- stated as a test because an empty return looks like an
    oversight."""
    module = _build_inchi()
    module.get_platform = lambda: 'win-amd64'
    assert module.cmake_args() == [], 'build_inchi.py now passes cmake flags on Windows; they need a reason'


# --- the macOS floor and slices --------------------------------------------------------------------

def test_the_macos_dylib_takes_its_floor_and_slices_from_the_interpreter():
    """`macosx-10.9-universal2` must produce both slices and a 10.9 floor, not the build machine's.

    Without the floor the wheel is tagged for whatever macOS built it -- `macosx_26_0_universal2`,
    installable nowhere earlier, whatever the extension itself was compiled for.  Without the
    architectures the dylib is the build machine's arch alone, and a universal2 wheel installs on an
    Intel Mac with no InChI and no error.
    """
    module = _build_inchi()
    module.get_platform = lambda: 'macosx-10.9-universal2'
    flags = module.cmake_args()
    assert '-DCMAKE_OSX_DEPLOYMENT_TARGET=10.9' in flags, \
        'no deployment target: `wheel` raises the wheel tag to the macOS version of the machine that ' \
        'built the dylib, and the release installs on that version and nowhere earlier -- got %s' % flags
    assert '-DCMAKE_OSX_ARCHITECTURES=arm64;x86_64' in flags, \
        'a universal2 interpreter needs both slices, or InChI is silently absent on one arch -- got %s' % flags


def test_a_single_arch_macos_interpreter_gets_that_arch_and_nothing_else():
    """What CI actually builds: `setup-python` is one arch per runner, so the flags must not widen it to
    universal2 -- an x86_64 slice on an arm64 runner needs an SDK that may not be there."""
    module = _build_inchi()
    module.get_platform = lambda: 'macosx-11.0-arm64'
    flags = module.cmake_args()
    assert '-DCMAKE_OSX_ARCHITECTURES=arm64' in flags, flags
    assert '-DCMAKE_OSX_DEPLOYMENT_TARGET=11.0' in flags, flags


def test_the_linux_build_is_told_nothing_about_architectures():
    """The macOS flags are macOS-only: cmake warns on an unused `CMAKE_OSX_*`, and a Linux cross-build is
    not something this file expresses."""
    module = _build_inchi()
    module.get_platform = lambda: 'linux-aarch64'
    assert module.cmake_args() == ['-DCMAKE_C_FLAGS=%s' % SIGNED_CHAR, '-DCMAKE_SHARED_LINKER_FLAGS=%s' % STRIP], \
        module.cmake_args()


def test_the_libinchi_job_pins_the_interpreter_it_reads_the_platform_from():
    """One dylib per OS serves five interpreters, so on macOS this job's python sets the floor for all
    five wheels.  Left to the runner's default it is the runner image's macOS version, which raises every
    mac wheel's tag -- so the job pins the oldest interpreter in the wheel matrix, whose floor cannot
    raise any of them."""
    job = _jobs().get('libinchi')
    assert job, 'the release workflow has no libinchi job, or the job scanner has stopped working'
    assert 'actions/setup-python' in job, (
        'the libinchi job builds with the runner image\'s default python.  `build_inchi.py` reads '
        '`sysconfig.get_platform()`, so on macOS that interpreter decides the minimum OS version of '
        'every wheel carrying this dylib.')


# --- debug info stays out of the wheel -------------------------------------------------------------

def test_the_linux_extension_asks_for_no_debug_info():
    """distutils compiles an extension with CPython's own `CFLAGS`, and those carry `-g`.

    Measured on 3.0's cp312 Linux wheel: `.debug*` was 17.50 MB of a 20.82 MB `_core.so` whose `.text` is
    2.52 MB, against a 2.18 MB `.pyd` for the same code on Windows.  These flags are appended after
    CPython's and gcc takes the last of `-g`/`-g0`, which is what lets a flag here cancel one from there.
    """
    linux = [flags for flags in _extra_compile_args() if '-O3' in flags]
    assert len(linux) == 1, \
        'setup.py no longer has exactly one -O3 branch, so the Linux branch is not identifiable here'
    assert NO_DEBUG_INFO in linux[0], (
        'the Linux extension is compiled without %s, so the `-g` in CPython\'s CFLAGS stands and the '
        'wheel ships DWARF -- 84%% of `_core.so` when this was last measured' % NO_DEBUG_INFO)


def test_libinchi_is_stripped_at_link_time_on_linux_and_nowhere_else():
    """`CMAKE_BUILD_TYPE=Release` does not remove InChI's `-g`, because InChI's own CMakeLists adds it.

    `INCHI_API/libinchi/src/CMakeLists.txt` gives gcc-like compilers `-g;-O1` through
    `target_compile_options`, which lands after both `CMAKE_C_FLAGS` and the config's flags -- so a `-g0`
    passed in cannot win and only a link-time strip reaches it.  3.05 MB of a 4.36 MB `libinchi.so`.
    """
    module = _build_inchi()
    for platform, wanted in (('linux-x86_64', True), ('linux-aarch64', True),
                             ('macosx-10.9-universal2', False), ('macosx-11.0-arm64', False),
                             ('win-amd64', False)):
        module.get_platform = lambda p=platform: p
        stripping = any(STRIP in flag for flag in module.cmake_args())
        assert stripping is wanted, (
            '%s: cmake_args() %s the link-time strip.  Without it on Linux libinchi carries 3 MB of '
            'DWARF into every wheel; with it on macOS, ld64 deprecates `-s` and the linked Mach-O has no '
            'DWARF to remove.' % (platform, 'lacks' if wanted else 'should not pass'))


def _elf_debug_bytes(data):
    """(total, `.debug*` bytes) for an ELF image, or None when it is not one.

    `struct` and not `readelf`: the property is worth checking on any Linux row of the test matrix, and a
    binutils dependency would make it skip instead.
    """
    if data[:4] != b'\x7fELF' or data[4] != 2:            # ELF64 only; nothing here builds 32-bit
        return None
    shoff, = unpack_from('<Q', data, 0x28)
    shentsize, shnum, shstrndx = unpack_from('<HHH', data, 0x3a)
    names = unpack_from('<Q', data, shoff + shstrndx * shentsize + 24)[0]
    debug = 0
    for i in range(shnum):
        offset, = unpack_from('<I', data, shoff + i * shentsize)
        name = data[names + offset:data.index(b'\0', names + offset)]
        if name.startswith(b'.debug'):
            debug += unpack_from('<Q', data, shoff + i * shentsize + 32)[0]
    return len(data), debug


def test_a_built_linux_binary_carries_no_dwarf():
    """The property the two flags above are for, against build output whenever there is some.

    Vacuous off Linux and on a checkout nobody has built -- said out loud, because a silent pass for a
    structural reason is what the declarations above are the guard against.  This is the assertion that
    would survive a gcc that ignored `-g0` or a cmake that dropped the linker flag.
    """
    offenders = []
    for path in [*(ROOT / 'chython' / 'core').glob('_core*.so'), ROOT / 'build' / 'inchi' / 'libinchi.so']:
        if not path.exists():
            continue
        measured = _elf_debug_bytes(path.read_bytes())
        if measured and measured[1]:
            offenders.append('%s: %.2f MB of %.2f MB is .debug*'
                             % (path.name, measured[1] / 1048576, measured[0] / 1048576))
    assert not offenders, (
        'these binaries would ship their debug info inside the wheel:\n  %s\n'
        'The extension takes `%s` from setup.py and libinchi takes `%s` from build_inchi.py; a build '
        'that predates either needs redoing before it is measured.' % ('\n  '.join(offenders),
                                                                      NO_DEBUG_INFO, STRIP))


# --- the source distribution -----------------------------------------------------------------------

def test_the_release_publishes_a_source_distribution():
    """Wheels cover four platforms and five interpreters; everything else -- musl, FreeBSD, macOS
    x86_64, a CPython newer than this matrix -- gets `pip install chython` answering "no matching
    distribution" unless there is an sdist to fall back to."""
    jobs = _jobs()
    building = [name for name, job in jobs.items() if 'build --sdist' in job]
    assert building, (
        'no job runs `python -m build --sdist`: the release is wheels only, so PyPI has no source '
        'distribution and any platform outside the wheel matrix cannot install chython at all')
    assert len(building) == 1, \
        'more than one job builds the sdist (%s); it is a property of the source and every job would ' \
        'produce the same file' % ', '.join(sorted(building))
    assert 'twine upload' in jobs[building[0]], \
        'the %s job builds an sdist and does not upload it' % building[0]


def test_the_gate_can_fail():
    """Negative controls: every scanner here passes vacuously if it stops finding anything."""
    assert _extra_compile_args(), 'the setup.py compile-flag scanner found nothing'
    jobs = _jobs()
    assert 'binary' in jobs and 'libinchi' in jobs, \
        'the workflow job scanner is not reading the jobs: %s' % sorted(jobs)
    assert 'strategy' in jobs['binary'], 'the job scanner returns a truncated block'
    assert hasattr(_build_inchi(), 'cmake_args'), 'build_inchi.py has no cmake_args(); the flags moved somewhere else'
