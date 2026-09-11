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
"""`chython.formats` imports without dragging the whole library in behind it: nothing under
`formats/` imports `chython` itself.  A violation breaks no test -- the import still works -- it only
reverses the dependency arrows.  Compiling to a `.so` does not help: `import chython.core` still runs
`chython/__init__.py` on the way to the submodule.  Each test runs in a subprocess, since by
collection time some other test module has already imported the facade.
"""
from subprocess import run
from sys import executable


def _probe(source):
    """Run `source` in a fresh interpreter; return it, asserting a clean exit."""
    done = run([executable, '-c', source], capture_output=True, text=True)
    assert done.returncode == 0, done.stderr
    return done.stdout


def test_importing_formats_does_not_execute_the_facade():
    """The rule itself: `import chython.formats` must not run `chython/__init__.py`.

    A shim package is planted in `sys.modules` first, so creating the empty parent -- harmless and
    unavoidable -- is not mistaken for executing the facade, whose body would bind `clean2d_engine`.
    """
    out = _probe('''
import sys, types
# stand in a package object for `chython` so the submodule import does not have to create one by
# executing the real `__init__`. If `chython.formats` reaches back through the facade, the import
# system will run the real body and `clean2d_engine` will appear.
shim = types.ModuleType('chython')
shim.__path__ = [__import__('os').path.join(__import__('os').getcwd(), 'chython')]
shim._sentinel = True
sys.modules['chython'] = shim

import chython.formats

print('facade_executed', hasattr(sys.modules['chython'], 'clean2d_engine'))
print('sentinel_survived', getattr(sys.modules['chython'], '_sentinel', False))
''')
    assert 'facade_executed False' in out, \
        'importing chython.formats executed chython/__init__.py -- something in formats/ imports the ' \
        'facade (`from chython import ...`, or a `..` that resolves to the package root)'
    assert 'sentinel_survived True' in out


def test_formats_does_not_import_the_chython_two_containers():
    """The same rule one layer finer: `formats` is built on the core, not on chython 2.

    The listed packages are deleted, so importing one raises before this scan sees it.  The prefixes
    stay to catch a module restored out of git -- the cheapest way to re-create the dependency.
    """
    out = _probe('''
import sys, types
shim = types.ModuleType('chython')
shim.__path__ = [__import__('os').path.join(__import__('os').getcwd(), 'chython')]
sys.modules['chython'] = shim

import chython.formats

leaked = sorted(m for m in sys.modules
                if m.startswith(('chython.containers', 'chython.algorithms', 'chython.reactor',
                                 'chython.periodictable', 'chython.files')))
print('LEAKED', leaked)
''')
    assert 'LEAKED []' in out, \
        'chython.formats pulled in chython 2 packages: %s' % out.strip()


def test_the_probe_can_fail():
    """A negative control: importing the facade *does* trip both detectors.

    Without it, a broken `_probe` leaves both tests above green while measuring nothing.  The leak
    half watches `chython.core`, which the facade imports in its first statement and always will.
    """
    out = _probe('''
import sys
import chython
print('facade_executed', hasattr(sys.modules['chython'], 'clean2d_engine'))
leaked = [m for m in sys.modules if m.startswith('chython.core')]
print('LEAKED', sorted(leaked)[:1])
''')
    assert 'facade_executed True' in out, 'the detector no longer sees the facade being executed'
    assert 'LEAKED []' not in out, 'the sys.modules prefix scan no longer sees an imported subpackage'
