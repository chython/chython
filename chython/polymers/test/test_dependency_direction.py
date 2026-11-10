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
"""`core <- polymers`, and nothing here imports the facade or a sibling package."""
import ast
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent.parent                             # chython/
PACKAGE = 'polymers'

# The layer below, and this package itself.
ALLOWED = ('chython.core', 'chython.polymers')


def _imports(path: Path):
    """`(lineno, dotted target)` for every import in one file, absolute or relative."""
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
    """Every in-library import resolves to `chython.core` or to this package, `test/` included."""
    offences = []
    for path in sorted((ROOT / PACKAGE).rglob('*.py')):
        for lineno, target in _imports(path):
            if target != 'chython' and not target.startswith('chython.'):
                continue
            if target == 'chython' or not target.startswith(ALLOWED):
                offences.append(f'{path.relative_to(ROOT.parent)}:{lineno}: {target}')
    assert not offences, (f'chython.{PACKAGE} imports outside `core <- polymers`:\n  ' + '\n  '.join(offences))
