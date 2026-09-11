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
"""The descriptor tables load on first use, and none of them loads on import.

`tables/maccs.tsv` compiles 158 queries and `tables/qed_alerts.tsv` 64, which is the cost `import
chython` would pay on a path that never asks for a descriptor.  Asserted in a subprocess: the caches are
module-level dicts, so once any test in the session has touched one, an in-process check would read a
warm cache and pass whatever the import does.
"""
from pathlib import Path
from subprocess import run
from sys import executable


#: The repo root, so the subprocess imports this checkout and not an installed copy.
ROOT = Path(__file__).resolve().parents[3]

#: One line per cache, printed as `NAME<tab>True|False`.  The import must leave every one of them cold.
_SCRIPT = """
import sys
sys.path.insert(0, '__PACKAGE_ROOT__')
import chython.chemistry
from chython.chemistry import _tables as t

for name in ('_MACCS_CACHE', '_MACCS_CORPUS_CACHE', '_QED_ALERTS_CACHE'):
    sys.stdout.write('%s\\t%s\\n' % (name, bool(getattr(t, name))))
sys.stdout.write('NUMPY\\t%s\\n' % ('numpy' in sys.modules))
"""


def _probe(root):
    # substitution rather than `%`, because the script formats its own output with `%s`
    script = _SCRIPT.replace('__PACKAGE_ROOT__', str(root))
    result = run([executable, '-c', script], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    return dict(line.split('\t') for line in result.stdout.splitlines())


def test_no_descriptor_table_loads_on_import():
    state = _probe(ROOT)
    cold = [name for name, loaded in state.items() if name != 'NUMPY' and loaded == 'False']
    assert len(cold) == 3, state
    assert 'True' not in [state[n] for n in state if n != 'NUMPY'], state


def test_importing_chemistry_does_not_import_numpy():
    """numpy is the `ml` extra and `chython.chemistry` is in the base install.

    `_maccs.py` imports numpy inside `maccs_keys`, not at module level, precisely so that this holds --
    and a module-level import there is invisible until something measures it, since a dev environment
    has numpy either way.
    """
    state = _probe(ROOT)
    assert state['NUMPY'] == 'False', 'importing chython.chemistry pulled numpy in'


def test_the_tables_load_on_first_use_and_stay_loaded():
    from chython.chemistry._tables import maccs_corpus, maccs_rules, qed_alerts

    assert len(maccs_rules()) == 166
    assert maccs_rules() is maccs_rules()          # the second call is the cache
    assert len(maccs_corpus()) == 330              # 165 keys with a definition, set and unset
    assert maccs_corpus() is maccs_corpus()
    assert qed_alerts() is qed_alerts()
