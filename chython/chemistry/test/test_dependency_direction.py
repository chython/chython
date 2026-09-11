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
"""`core <- chemistry`: nothing in that chain imports the facade, and the tables do not import the
passes.  Checked twice -- statically, so a failure names the offending line, and at runtime with the
facade replaced by an empty stub.  `import chython.chemistry` alone would NOT prove the runtime half:
Python runs the parent package's `__init__` before the child's, however clean the child is.
"""
import ast
from pathlib import Path
from subprocess import run
from sys import executable

import pytest


ROOT = Path(__file__).resolve().parent.parent.parent                             # chython/
PACKAGE = 'chemistry'

# The layer below plus this package itself.
ALLOWED = ('chython.core', 'chython.chemistry')

# The modules that read the tables, and must not come to depend on the modules that apply them, so
# that a row stays reviewable without reading a pass.
KNOWLEDGE = ('_tables.py', '_smarts.py', '_residues.py')
PASSES = ('chython.chemistry._standardize', 'chython.chemistry._resonance',
          'chython.chemistry._implicit', 'chython.chemistry._counts',
          'chython.chemistry._crippen', 'chython.chemistry._maccs',
          'chython.chemistry._pharmacophore', 'chython.chemistry._qed',
          'chython.chemistry._tpsa')


def _imports(path: Path):
    """`(lineno, dotted target)` for every import in one file, absolute or relative.

    Relative imports resolve against the file's own package, so `from ..core import x` and
    `from chython.core import x` are the same fact.  The last path component is dropped
    unconditionally, module or `__init__.py` alike: level 1 means "my package", which for
    `chemistry/_tables.py` is `chython.chemistry` and not `chython.chemistry._tables`.  A target one
    component too long still passes the layer rule, so `test_the_loader_rule_can_fail` pins the
    resolution itself.
    """
    tree = ast.parse(path.read_text(encoding='utf-8'))
    parts = path.relative_to(ROOT.parent).with_suffix('').parts[:-1]

    out = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            out.extend((node.lineno, alias.name) for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            if not node.level:
                out.append((node.lineno, node.module or ''))
            else:
                # level 1 is this package, level 2 its parent, and so on
                base = parts[:len(parts) - node.level + 1]
                out.append((node.lineno, '.'.join(base + ((node.module,) if node.module else ()))))
    return out


def test_nothing_in_the_package_imports_above_itself():
    """Every in-library import resolves to `chython.core` or to this package.  Statically.

    `test/` is included deliberately: a test reaching for the facade to build a fixture breaks
    isolation as effectively as production code would.
    """
    offences = []
    for path in sorted((ROOT / PACKAGE).rglob('*.py')):
        for lineno, target in _imports(path):
            if not target.startswith('chython'):
                continue                                           # stdlib or a third party
            if target == 'chython' or not target.startswith(ALLOWED):
                offences.append(f'{path.relative_to(ROOT.parent)}:{lineno}: {target}')
    assert not offences, (
        f'chython.{PACKAGE} imports outside `core <- chemistry`:\n  '
        + '\n  '.join(offences) + '\n\nThe dependency direction is what keeps the '
        '`lazy_object_proxy` wrappers out of chython/__init__.py.  Move the code down a '
        'layer or pass the value in; do not widen ALLOWED without a ruling.')


@pytest.mark.parametrize('module', KNOWLEDGE)
def test_the_table_loaders_do_not_import_the_passes(module):
    """The knowledge half stays readable without the passes: a loader imports no pass."""
    offences = [f'{module}:{lineno}: {target}'
                for lineno, target in _imports(ROOT / PACKAGE / module)
                if target in PASSES]
    assert not offences, (
        f'chython/chemistry/{module} reads the tables and must not import the code that applies '
        'them:\n  ' + '\n  '.join(offences) + '\n\nThe direction is loader -> pass and never back. '
        'If a pass has a value the loader needs, that value is knowledge and belongs in a table '
        'column or in the loader.')


def test_the_loader_rule_can_fail():
    """Negative control: pins the RESOLVED names, not just "something matched".

    A typo in `PASSES` or an off-by-one in `_imports` would leave the rule above green forever.
    Asserting the exact targets fails if the resolver drifts by one component either way.
    """
    targets = {target for _, target in _imports(ROOT / PACKAGE / '_resonance.py')}
    expected = {'chython.chemistry._implicit', 'chython.chemistry._standardize'}
    assert expected <= targets, (
        '`_resonance.py` imports `_implicit` and `_standardize`, and the scanner must resolve both '
        f'to their real dotted names.  It reported:\n  ' + '\n  '.join(sorted(targets))
        + f'\n\nmissing: {sorted(expected - targets)}')
    assert expected <= set(PASSES), 'PASSES no longer names the passes `_resonance.py` imports'


# Run with `chython` replaced by an empty package, so the facade is unreachable rather than unused.
_SCRIPT = """
import sys, types

stub = types.ModuleType('chython')
stub.__path__ = ['__PACKAGE_ROOT__']
sys.modules['chython'] = stub

import chython.chemistry
from chython.core import read_smiles

# and prove it does something, not just that it imports: the pass is registered onto the core
# container by chython.chemistry, and that hook is the whole interface between the layers
molecule = read_smiles('CN(=O)=O')
assert molecule.standardize(), 'the pass did not fire on a pentavalent nitro group'
assert molecule.smiles == 'C[N+]([O-])=O', molecule.smiles

leaked = sorted(m for m in sys.modules if m.startswith('chython.') and not
                m.startswith(('chython.core', 'chython.chemistry')))
sys.stdout.write('LEAKED\\t%s\\n' % ','.join(leaked))
sys.stdout.write('PROXY\\t%s\\n' % ('lazy_object_proxy' in sys.modules))
sys.stdout.write('FACADE\\t%s\\n' % (sys.modules['chython'] is stub))
"""


def test_the_package_works_with_the_facade_never_executed():
    """The claim, executed: standardize a molecule in an interpreter where `chython` is empty."""
    # substitution rather than `%`, because the script formats its own output with `%s`
    script = _SCRIPT.replace('__PACKAGE_ROOT__', str(ROOT))
    result = run([executable, '-c', script], capture_output=True, text=True, cwd=str(ROOT.parent))
    assert result.returncode == 0, (
        'chython.chemistry cannot be used without the facade:\n' + result.stderr)

    reported = dict(line.split('\t') for line in result.stdout.splitlines() if '\t' in line)
    assert reported['FACADE'] == 'True', 'something replaced the stub with the real facade'
    assert reported['LEAKED'] == '', (
        f"importing this package pulled in {reported['LEAKED']}.  The static test above "
        'should have caught it; if it did not the import is dynamic, and a dynamic import of the '
        'facade is the same dependency wearing a hat')
    assert reported['PROXY'] == 'False', (
        'lazy_object_proxy was imported, so something on this path still needs the facade to be '
        'lazy.  That library is what the layout exists to keep deleted')
