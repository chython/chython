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
"""numpy is optional, and the library is usable without it.

The motivating deployment is a serverless function, where the bundle is uploaded on every deploy and
numpy is by a wide margin the largest dependency chython can pull in.  Nothing on the representation
path needs it -- parse, standardize, kekulize, canonicalize, depict, react, write -- so nothing on the
representation path may import it, and the array-answering surface has to fail at the call with an
error that names the extra rather than at `import chython` with a traceback from inside a featurizer.

WHY THESE TESTS RUN IN SUBPROCESSES.  Both questions are about the import graph, and this process has
already answered them the wrong way round: pytest has imported chython, and numpy is in `sys.modules`
because some other test asked for a fingerprint.  Neither question can be asked in-process, so each
one gets a fresh interpreter -- the same reason `core/test/oracle.py` spawns rather than imports.

THIS FILE IS A RATCHET.  A single module-level `from numpy import ...` anywhere `chython.chemistry`
imports eagerly -- a featurizer, say -- drags numpy into the façade for every caller.  A measurement
written into a comment decays silently; this one fails the build.
"""
from pathlib import Path
from re import DOTALL, MULTILINE, search
from subprocess import run
from sys import executable


ROOT = Path(__file__).resolve().parent.parent.parent

#: Refuse numpy to everything downstream, from inside the interpreter under test.
#:
#: A `sys.meta_path` finder rather than uninstalling numpy or scrubbing `sys.path`: the point is to
#: reproduce a machine where numpy was never installed, while still importing the chython under test
#: from this checkout.  Raising from `find_spec` -- rather than returning None -- is what makes the
#: failure look like absence instead of falling through to the real finders.
_BLOCK_NUMPY = '''
import sys


class _NoNumpy:
    def find_spec(self, name, path=None, target=None):
        if name == 'numpy' or name.startswith('numpy.'):
            raise ImportError("No module named 'numpy'")
        return None


sys.meta_path.insert(0, _NoNumpy())
for _name in [n for n in sys.modules if n == 'numpy' or n.startswith('numpy.')]:
    del sys.modules[_name]
'''


def _python(script):
    """Run `script` in a fresh interpreter with this checkout importable.  Returns the CompletedProcess.

    No `-I`: unlike `core/test/oracle.py`, the whole point here is to import the tree under test, and `-I`
    drops the working directory from `sys.path`.  `cwd=ROOT` is what puts this checkout first.
    """
    return run([executable, '-c', script], cwd=ROOT, capture_output=True, text=True, timeout=300)


def test_importing_chython_does_not_import_numpy():
    """The façade must not pay for, or require, an optional dependency."""
    out = _python('import sys\n'
                  'import chython\n'
                  "assert 'numpy' not in sys.modules, sorted(n for n in sys.modules "
                  "if n.startswith('numpy'))\n"
                  "print('clean')")
    assert out.returncode == 0, f'`import chython` imported numpy or failed:\n{out.stderr}'
    assert 'clean' in out.stdout


def test_the_representation_path_works_with_numpy_absent():
    """Parse, repair, canonicalize and write, on an interpreter where numpy cannot be imported.

    This is the deployment being bought, so it is asserted end to end rather than by importing the
    modules and trusting them.  Aspirin covers a ring, an aromatic system to kekulize and re-perceive,
    a carboxylic acid for standardization, and stereo-free CTfile output.

    NO CANONICAL STRING IS PINNED HERE.  This test owns one question -- does the path run without
    numpy -- and a literal expected SMILES would make it fail for the unrelated reason that the
    canonical order changed, which `core/test/` already covers and covers better.  What is asserted
    instead is the property that cannot hold by accident: canonicalizing is idempotent, and reparsing
    the output reproduces it.  Both walk the whole pipeline; neither cares what it spells.
    """
    out = _python(_BLOCK_NUMPY + '''
from chython import smiles, mol, smarts

m = smiles('CC(=O)Oc1ccccc1C(=O)O')
assert m.atom_count == 13
m.kekule()
m.standardize()
m.thiele()
m.canonicalize()
once = str(m)
m.canonicalize()
assert str(m) == once, (once, str(m))
again = smiles(once)
again.canonicalize()
assert str(again) == once, (once, str(again))

# the writers, and the reader that has to come back through the arena
text = mol(m)
assert 'V2000' in text or 'V3000' in text
assert again == m, 'reparsing the canonical form gave a different molecule'

# substructure matching, which is the other half of the library people deploy
assert smarts('[O;D1]-[C;D3](=O)-[C;a]') < m, 'the aryl carboxylic acid did not match'
assert len(m.split()) == 1

# the descriptors that are NOT distance-based stay available, which is the boundary `chython[ml]`
# draws inside one surface: these four read tables and counts, so they never reach the binder.
assert m.tpsa > 0
assert isinstance(m.crippen_logp, float)
assert m.rings_count == 1
assert m.bertz_ct > 0

# and QED with them: its eight inputs are those tables and counts, so the score is numpy-free even
# though `maccs_keys` -- which shares the aromatic ring count with it -- is not
assert 0. < m.qed < 1.

import sys
assert 'numpy' not in sys.modules, sorted(n for n in sys.modules if n.startswith('numpy'))
print('representation path clean')
''')
    assert out.returncode == 0, f'the representation path needed numpy:\n{out.stderr}'
    assert 'representation path clean' in out.stdout


def test_the_array_surface_names_the_extra_when_numpy_is_absent():
    """All twenty-five numpy-backed entry points raise ImportError naming `chython[ml]`.

    THE LIST IS EXHAUSTIVE ON PURPOSE, and it is the reason this test earns its runtime.  A caller who
    hits one of these on a minimal install has to be told which extra to install, and the ones most
    likely to be missed are the ones that do not look like array methods: `morgan_bit_set`,
    `morgan_hash_set`, `morgan_hash_counts` and their `linear_` counterparts answer a plain set or
    dict, and still need numpy because every fingerprint spelling builds the same invariant vector
    first; `maccs_bit_set` answers a frozenset over the `uint8[167]` vector `maccs_keys` fills.
    Likewise `wiener_index` and friends answer a number, and need numpy because `distance_matrix` is
    the only shortest-path code in the core and they all read its output.
    `rxn.modeling_view` answers dicts, but it is a dict assembly over the same invariant arrays and
    needs numpy for exactly the same reason.  Guessing which of these was safe is exactly the mistake
    the enumeration prevents.

    Each name is called separately, and a `lambda` rather than a `getattr` because five of them are
    properties: `m.wiener_index` raises on attribute access, so the call has to be deferred by
    something that also defers an attribute read.
    """
    out = _python(_BLOCK_NUMPY + '''
from chython import smiles
from chython.chemistry import pharmacophore_invariants

m = smiles('CC(=O)Oc1ccccc1C(=O)O')
rxn = smiles('[CH3:1][OH:2]>>[CH3:1][NH2:3]')
calls = {
    # the folded fingerprints -- these do look like array methods
    'morgan_fingerprint': lambda: m.morgan_fingerprint(),
    'morgan_count_vector': lambda: m.morgan_count_vector(),
    'linear_fingerprint': lambda: m.linear_fingerprint(),
    'linear_count_vector': lambda: m.linear_count_vector(),
    # ...and the unfolded ones, which answer a set or a dict and need numpy anyway
    'morgan_hash_counts': lambda: m.morgan_hash_counts(),
    'morgan_hash_set': lambda: m.morgan_hash_set(),
    'morgan_bit_set': lambda: m.morgan_bit_set(),
    'linear_hash_counts': lambda: m.linear_hash_counts(),
    'linear_hash_set': lambda: m.linear_hash_set(),
    'linear_bit_set': lambda: m.linear_bit_set(),
    # the vector everything above is built on, and the two matrices
    'atom_invariants': lambda: m.atom_invariants(),
    'adjacency_matrix': lambda: m.adjacency_matrix(),
    'distance_matrix': lambda: m.distance_matrix(),
    # the descriptors derived from the distance matrix -- four of these five are properties
    'eccentricities': lambda: m.eccentricities(),
    'wiener_index': lambda: m.wiener_index,
    'graph_radius': lambda: m.graph_radius,
    'graph_diameter': lambda: m.graph_diameter,
    'balaban_j': lambda: m.balaban_j,
    # and the ones in `chemistry`, which reach the core's message through `require_numpy`
    'pharmacophore_invariants': lambda: pharmacophore_invariants(m),
    'maccs_keys': lambda: m.maccs_keys(),
    # `maccs_bit_set` answers a frozenset and is here for `morgan_bit_set`'s reason: it reads the
    # `uint8[167]` vector off `maccs_keys` and only then picks the set bits out of it
    'maccs_bit_set': lambda: m.maccs_bit_set(),
    # the ML views -- mol.transition_view and rxn.transition_view are distinct containers
    'm.state_view': lambda: m.state_view(),
    'm.transition_view': lambda: m.transition_view(),
    'rxn.transition_view': lambda: rxn.transition_view(),
    # rxn.modeling_view answers dicts, but it is a dict assembly over the same arrays
    'rxn.modeling_view': lambda: rxn.modeling_view(),
}
for name, call in calls.items():
    try:
        call()
    except ImportError as e:
        assert 'chython[ml]' in str(e), f'{name} does not name the extra: {e}'
    else:
        raise AssertionError(f'{name} answered without numpy; it is not on the numpy path any more')
print('all', len(calls), 'named the extra')
''')
    assert out.returncode == 0, f'an array entry point did not name the extra:\n{out.stderr}'
    assert 'named the extra' in out.stdout


def test_numpy_is_declared_as_an_extra_and_not_as_a_runtime_dependency():
    """A packaging test, so the install shape cannot drift back from the code shape.

    Regex and not `tomllib` for `test_packaging.py`'s reason: `tomllib` is 3.11 and this package
    supports 3.10, so a parsed version of this test would skip on the oldest interpreter people build
    wheels on -- which is the same as not having it.
    """
    text = (ROOT / 'pyproject.toml').read_text(encoding='utf-8')

    runtime = search(r'^dependencies = \[(.*?)^\]', text, DOTALL | MULTILINE)
    assert runtime, 'pyproject.toml has no [project] dependencies array'
    requirements = [line.strip().strip(",'\"") for line in runtime.group(1).splitlines()
                    if line.strip() and not line.strip().startswith('#')]
    assert not [r for r in requirements if r.startswith('numpy')], \
        f'numpy is back in the runtime dependencies: {requirements}'

    extras = search(r'^\[project\.optional-dependencies\]\n(.*?)(?=^\[|\Z)', text,
                    DOTALL | MULTILINE)
    assert extras, 'pyproject.toml has no [project.optional-dependencies] section'
    ml = search(r"^ml = \[(.*?)\]", extras.group(1), MULTILINE)
    assert ml, 'there is no `ml` extra for the array surface to point at'
    assert 'numpy' in ml.group(1), f'the `ml` extra does not provide numpy: {ml.group(1)}'
