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
"""chytorch's compiled `_unpack`, consulted out of process, as the differential for `state_view`.

    CHYTORCH_PATH=/path/to/chytorch pytest chython/core/test/test_ml_unpack_differential.py

WHY A SUBPROCESS.  The child must not import chython: an in-process comparison against a library that
also imports chython proves nothing about which implementation produced a number.  It receives pach v2
records as base64 and returns arrays as JSON, so its only chython dependency is the byte format.

`-I` IS LOAD-BEARING.  Without it the child inherits `sys.path` and `PYTHONPATH` from the parent, the
parent's chython is importable, and a differential silently becomes a tautology.  The extension is
loaded by file path, not via a package import, so no `__init__` runs and neither `chython` nor
`chytorch` ever reaches `sys.modules`.

NO VERSION PIN.  chytorch declares no release version, so `verify()` asserts what the comparison
actually needs -- that the child loaded chytorch's `_unpack` extension and that neither `chython` nor
`chytorch` is in its `sys.modules` -- and `test_the_harness_can_disagree` proves the channel can
report a mismatch.
"""
import json
from base64 import b64encode
from os import environ
from pathlib import Path
from subprocess import PIPE, run as _run
from sys import executable

from pytest import mark, skip


__all__ = ['ENV_VAR', 'probe', 'require', 'requires_oracle', 'unpack', 'verify']

ENV_VAR = 'CHYTORCH_PATH'
ISOLATION = '-I'

# The extension is loaded by file path so that no chytorch `__init__.py` runs and the child's
# `sys.modules` remains free of both `chython` and `chytorch`. {so_glob!r} is the only placeholder.
_PREAMBLE = """
import importlib.util as _iu, glob as _glob, json, sys


def _emit(payload):
    sys.stdout.write('@@' + json.dumps(payload) + '@@')


_so_candidates = _glob.glob({so_glob!r})
if not _so_candidates:
    raise RuntimeError('_unpack extension not found; set CHYTORCH_PATH to the checkout root')
_spec = _iu.spec_from_file_location('_unpack', _so_candidates[0])
_mod = _iu.module_from_spec(_spec)
_spec.loader.exec_module(_mod)
unpack_graph = _mod.unpack_graph
_UNPACK_PATH = _so_candidates[0]
"""


def path():
    """Where chytorch lives.  Raises ``KeyError`` when ``CHYTORCH_PATH`` is unset."""
    return Path(environ[ENV_VAR])


def _so_glob():
    """Glob pattern for the compiled `_unpack` extension in the chytorch checkout."""
    return str(path() / 'chytorch' / 'utils' / 'data' / 'molecule' / '_unpack*.so')


def run(source, payload=None):
    """Run `source` in a child that has chytorch's `_unpack` and no chython; return what it emitted."""
    script = _PREAMBLE.format(so_glob=_so_glob()) + source
    argv = [executable, ISOLATION, '-c', script]
    result = _run(argv, stdout=PIPE, stderr=PIPE, input=json.dumps(payload or {}), text=True)
    if result.returncode:
        raise RuntimeError(f'chytorch oracle exited {result.returncode}:\n{result.stderr}')
    head, _, rest = result.stdout.partition('@@')
    body, _, _ = rest.partition('@@')
    if not body:
        raise RuntimeError(f'chytorch oracle emitted nothing:\n{result.stdout}\n{result.stderr}')
    return json.loads(body)


def probe():
    """The child's `_unpack` path, and whether `chython` or `chytorch` reached the child."""
    return run("""
_emit({'unpack': _UNPACK_PATH, 'chython': 'chython' in sys.modules,
       'chytorch': 'chytorch' in sys.modules})
""")


def verify():
    """Assert the two properties the differential rests on, or raise saying which failed."""
    info = probe()
    if 'chytorch' not in info['unpack']:
        raise RuntimeError(f"the child loaded an _unpack that is not chytorch's: {info['unpack']}")
    if info['chython']:
        raise RuntimeError('chython reached the oracle child; the differential is a tautology')
    if info['chytorch']:
        raise RuntimeError('chytorch reached the oracle child via package import; '
                           'the extension must be loaded by file path, not by package')
    return info


def available():
    """True when `verify()` passes; False on any exception. Called once at import time."""
    try:
        verify()
    except Exception:
        return False
    return True


def require():
    """Skip the calling test when chytorch is not reachable.

    Set ``CHYTORCH_PATH`` to a chytorch checkout to enable the differential.
    """
    try:
        verify()
    except Exception as e:
        skip(f'chytorch oracle unavailable ({e}); '
             f'set {ENV_VAR} to a chytorch checkout to enable the differential')


requires_oracle = mark.skipif(not available(),
                              reason=(f'chytorch oracle unavailable; '
                                      f'set {ENV_VAR} to a chytorch checkout to enable the differential'))


def unpack(records, max_neighbors=14, max_distance=10):
    """`_unpack` over pach v2 records: one dict of `atoms`, `neighbors`, `distances` each.

    Uses `unpack_graph(data, max_neighbors, max_distance)` -- no CLS row, no symmetric attention,
    cross-component value is 1 (matching `disconnected=1` in the encoding).
    """
    payload = {'records': [b64encode(r).decode('ascii') for r in records],
               'max_neighbors': max_neighbors, 'max_distance': max_distance}
    return run("""
from base64 import b64decode
payload = json.loads(sys.stdin.read())
out = []
for encoded in payload['records']:
    atoms, neighbors, distances = unpack_graph(b64decode(encoded),
                                               payload['max_neighbors'], payload['max_distance'])
    out.append({'atoms': atoms.tolist(), 'neighbors': neighbors.tolist(),
                'distances': distances.tolist()})
_emit(out)
""", payload)
