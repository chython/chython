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
"""`core <- reactions`, nothing here imports the facade, and the corpus does not import the enumerators.

Sideways is forbidden too: `chython.chemistry` is a sibling and not a layer below, which is why
`_tables.py` carries its own `read_table` rather than importing the one next door.
"""
import ast
from pathlib import Path
from subprocess import run
from sys import executable

import pytest


ROOT = Path(__file__).resolve().parent.parent.parent                             # chython/
PACKAGE = 'reactions'

# The layer below, and this package itself.  `chython.chemistry` is deliberately NOT here: it is a
# sibling.
ALLOWED = ('chython.core', 'chython.reactions')

# The corpus half, which must not come to depend on the half that decides which templates to try.
KNOWLEDGE = ('_tables.py',)
DRIVERS = ('chython.reactions._enumerate',)


def _imports(path: Path):
    """`(lineno, dotted target)` for every import in one file, absolute or relative.

    The last path component is dropped unconditionally: level 1 means "my package", so
    `reactions/_tables.py` resolves to `chython.reactions`.  Getting that wrong leaves the layer rule
    green and breaks the corpus rule, which needs an exact name -- `test_the_corpus_rule_can_fail` pins
    the resolution.
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
                base = parts[:len(parts) - node.level + 1]
                out.append((node.lineno, '.'.join(base + ((node.module,) if node.module else ()))))
    return out


def test_nothing_in_the_package_imports_above_or_beside_itself():
    """Every in-library import resolves to `chython.core` or to this package.  Statically.

    `test/` is included: a fixture reaching for the facade makes the package untestable in isolation
    just as effectively as production code would.
    """
    offences = []
    for path in sorted((ROOT / PACKAGE).rglob('*.py')):
        for lineno, target in _imports(path):
            # `chython.` WITH THE DOT: `chython_rxnmap` is a separate distribution -- the model weights,
            # which `attention/_session.py` imports -- and a prefix test without the boundary reads it as
            # this library.
            if target != 'chython' and not target.startswith('chython.'):
                continue                                           # stdlib or a third party
            if target == 'chython' or not target.startswith(ALLOWED):
                offences.append(f'{path.relative_to(ROOT.parent)}:{lineno}: {target}')
    assert not offences, (
        f'chython.{PACKAGE} imports outside `core <- reactions`:\n  '
        + '\n  '.join(offences) + '\n\nThe dependency direction is what keeps `lazy_object_proxy` '
        'out of the facade, and `chython.chemistry` is a sibling rather than a layer below. Move the '
        'code down a layer or pass the value in; do not widen ALLOWED without a ruling.')


@pytest.mark.parametrize('module', KNOWLEDGE)
def test_the_corpus_loader_does_not_import_the_enumerators(module):
    """The knowledge half stays readable without the code that drives it."""
    offences = [f'{module}:{lineno}: {target}'
                for lineno, target in _imports(ROOT / PACKAGE / module)
                if target in DRIVERS]
    assert not offences, (
        f'chython/reactions/{module} composes the corpus and must not import the code that '
        'enumerates it:\n  ' + '\n  '.join(offences) + '\n\nThe direction is loader -> driver and '
        'never back.  If the driver has a value the loader needs, that value is knowledge and belongs '
        'in a table column.')


def test_the_corpus_rule_can_fail():
    """Negative control, pinning the RESOLVED NAMES and not merely "something matched".

    Without it a typo in `DRIVERS` or an off-by-one in `_imports` leaves the rule above green forever.
    """
    targets = {target for _, target in _imports(ROOT / PACKAGE / '_enumerate.py')}
    expected = {'chython.reactions._tables', 'chython.core'}
    assert expected <= targets, (
        '`_enumerate.py` imports `_tables` and `chython.core`, and the scanner must resolve both to '
        'their real dotted names.  It reported:\n  ' + '\n  '.join(sorted(targets))
        + f'\n\nmissing: {sorted(expected - targets)}')
    # and the direction the rule screens for is the reverse of the one just proven to exist
    assert 'chython.reactions._enumerate' in DRIVERS, 'DRIVERS no longer names the driver module'


# Run with `chython` replaced by an empty package: the facade is unreachable, not merely unused.
_SCRIPT = """
import sys, types

stub = types.ModuleType('chython')
stub.__path__ = ['__PACKAGE_ROOT__']
sys.modules['chython'] = stub

import chython.reactions
from chython.core import read_smiles

# and prove it does something, not just that it imports: the injection hook is the whole interface
acid = read_smiles('CC(=O)O')
amine = read_smiles('CCN')
names = sorted(r.name for r in acid @ amine)
assert names == ['amidation'], names
assert 'carboxylic_acid' in acid.functional_groups(), acid.functional_groups()

leaked = sorted(m for m in sys.modules if m.startswith('chython.') and not
                m.startswith(('chython.core', 'chython.reactions')))
sys.stdout.write('LEAKED\\t%s\\n' % ','.join(leaked))
sys.stdout.write('PROXY\\t%s\\n' % ('lazy_object_proxy' in sys.modules))
sys.stdout.write('FACADE\\t%s\\n' % (sys.modules['chython'] is stub))
"""


def test_the_package_works_with_the_facade_never_executed():
    """The claim, executed: enumerate a reaction in an interpreter where `chython` is empty.

    `mol @ mol` specifically, because `__matmul__` resolves through the type's slot: it is compiled into
    the core and its body arrives by injection, so it is the piece likeliest to need the facade.
    """
    # substitution rather than `%`, because the script formats its own output with `%s`
    script = _SCRIPT.replace('__PACKAGE_ROOT__', str(ROOT))
    result = run([executable, '-c', script], capture_output=True, text=True, cwd=str(ROOT.parent))
    assert result.returncode == 0, (
        'chython.reactions cannot be used without the facade:\n' + result.stderr)

    reported = dict(line.split('\t') for line in result.stdout.splitlines() if '\t' in line)
    assert reported['FACADE'] == 'True', 'something replaced the stub with the real facade'
    assert reported['LEAKED'] == '', (
        f"importing this package pulled in {reported['LEAKED']}.  The static test above should have "
        'caught it; if it did not the import is dynamic, and a dynamic import of the facade is the '
        'same dependency wearing a hat')
    assert reported['PROXY'] == 'False', (
        'lazy_object_proxy was imported, so something on this path still needs the facade to be lazy. '
        'That library is what the layout exists to keep deleted, and this corpus is where every one '
        'of chython 2 s Proxy objects lived')
