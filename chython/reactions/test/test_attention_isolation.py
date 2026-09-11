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
"""The model is not loaded until the mapper is called, and the library works without it installed.

WHAT IS BEING PROTECTED.  The weights are an 84 MB file inside a separate distribution and ONNX Runtime
is tens of MB of shared library.  `chython.reactions` imports `attention/` at its own init to register the
container method, so a single module-level `from numpy import ...` or `import onnxruntime` in
`attention/__init__.py` puts all of it behind `import chython` -- for every caller, including the
serverless deployment `chython/test/test_optional_numpy.py` exists to keep small.

The three heavy names are checked together because they arrive together: `_encode.py` needs numpy,
`_session.py` needs the runtime and the weights, and both are imported inside the function body.

WHY SUBPROCESSES.  Every question here is about `sys.modules`, and this interpreter has already answered
them the wrong way round -- pytest imported chython, and another test in this package imported numpy to
build a matrix.  Each question gets a fresh interpreter, for the reason `test_optional_numpy.py` gives.
"""
from pathlib import Path
from subprocess import run
from sys import executable

from pytest import mark, raises, skip

from ...core import ReactionContainer, read_smiles as smiles
from ..attention import attention_available


ROOT = Path(__file__).resolve().parent.parent.parent.parent

#: The names that may not be in `sys.modules` after an import of this library.
HEAVY = ('numpy', 'onnxruntime', 'chython_rxnmap')


def _python(script):
    """Run `script` in a fresh interpreter with this checkout importable.

    No `-I`, for the reason `chython/test/test_optional_numpy.py` states: the tree under test is the one
    in the working directory, and `-I` drops it from `sys.path`.
    """
    return run([executable, '-c', script], cwd=str(ROOT), capture_output=True, text=True, timeout=300)


_LEAK_CHECK = '''
import sys

import %s

leaked = sorted(n for n in sys.modules if n.split('.')[0] in %r)
sys.stdout.write('LEAKED\\t%%s\\n' %% ','.join(leaked))
'''


@mark.parametrize('module', ['chython', 'chython.reactions', 'chython.reactions.attention'])
def test_the_model_is_absent_from_sys_modules_after_an_import(module):
    """Importing the façade, the package, or the mapper's own package loads none of the three.

    `chython.reactions.attention` is in the list deliberately: it is the module that would most
    reasonably import its own dependencies at the top, and the one place the rule is easiest to break by
    moving a function-level import up for readability.
    """
    out = _python(_LEAK_CHECK % (module, HEAVY))
    assert out.returncode == 0, f'`import {module}` failed:\n{out.stderr}'
    reported = dict(line.split('\t') for line in out.stdout.splitlines() if '\t' in line)
    assert reported['LEAKED'] == '', (
        f"`import {module}` loaded {reported['LEAKED']}.  The weights are 84 MB and the runtime is a "
        'shared library of comparable size; both are imported inside `attention_mapping`, and moving '
        'either import to module level makes every caller of `import chython` pay for them.')


def test_asking_whether_the_mapper_is_available_imports_nothing():
    """`attention_available()` answers from `importlib.util.find_spec`, which does not execute a module.

    A caller branches on this before deciding whether to map, so the question itself must be free --
    otherwise the answer `False` is the only cheap one and the answer `True` costs what it was asked to
    avoid.
    """
    out = _python('''
import sys

from chython import attention_available

answer = attention_available()
leaked = sorted(n for n in sys.modules if n.split('.')[0] in %r)
sys.stdout.write('ANSWER\\t%%s\\nLEAKED\\t%%s\\n' %% (answer, ','.join(leaked)))
''' % (HEAVY,))
    assert out.returncode == 0, f'`attention_available()` raised:\n{out.stderr}'
    reported = dict(line.split('\t') for line in out.stdout.splitlines() if '\t' in line)
    assert reported['ANSWER'] in ('True', 'False')
    assert reported['LEAKED'] == '', (
        f"`attention_available()` imported {reported['LEAKED']}.  It is a `find_spec` pair and must "
        'stay one: the question exists so a caller can avoid the cost, not pay it to ask.')


def test_the_facade_still_maps_when_the_model_is_installed():
    """The counterpart to the ratchet above: the lazy import resolves when called.

    Without this the leak test passes trivially on an installation where the import would fail anyway.
    """
    if not attention_available():
        skip('needs `chython[mapping]`: onnxruntime and chython_rxnmap')
    out = _python('''
import sys

from chython import smiles

rxn = smiles('CC(=O)O.CCN>>CC(=O)NCC.O')
result = rxn.attention_mapping()
sys.stdout.write('CHANGED\\t%s\\nRUNTIME\\t%s\\n'
                 % (bool(result), 'onnxruntime' in sys.modules))
''')
    assert out.returncode == 0, f'the mapper failed on a façade install:\n{out.stderr}'
    reported = dict(line.split('\t') for line in out.stdout.splitlines() if '\t' in line)
    assert reported['CHANGED'] == 'True'
    assert reported['RUNTIME'] == 'True', 'the mapper ran without importing the runtime it needs'


def test_the_mapper_names_the_extra_when_the_runtime_is_absent():
    """With `onnxruntime` unimportable the call raises ImportError naming `chython[mapping]`.

    NAMING THE EXTRA IS THE WHOLE POINT.  A caller on a minimal install gets one chance to learn what to
    install, and a bare `No module named 'onnxruntime'` from four frames down does not say which package
    of chython's asked for it.
    """
    out = _python('''
import sys


class _NoRuntime:
    def find_spec(self, name, path=None, target=None):
        if name.split('.')[0] in ('onnxruntime', 'chython_rxnmap'):
            raise ImportError("No module named %r" % name)
        return None


sys.meta_path.insert(0, _NoRuntime())

from chython import attention_available, smiles

assert not attention_available(), 'find_spec was supposed to be blocked'

rxn = smiles('CC=O.O>>CC(O)O')
try:
    rxn.attention_mapping()
except ImportError as e:
    sys.stdout.write('MESSAGE\\t%s\\n' % e)
else:
    sys.stdout.write('MESSAGE\\tno error\\n')
''')
    assert out.returncode == 0, f'the subprocess failed for an unrelated reason:\n{out.stderr}'
    reported = dict(line.split('\t') for line in out.stdout.splitlines() if '\t' in line)
    assert 'chython[mapping]' in reported['MESSAGE'], (
        f"the mapper reported {reported['MESSAGE']!r}, which does not name the extra to install")


def test_the_encoder_is_importable_on_its_own():
    """`_encode` and `_assign` need numpy and neither needs the runtime or the weights.

    The split is what lets the encoder differential and the assignment units run on a machine with no
    model; if the two files ever merge, those two suites become model-gated and the coverage claim they
    carry stops being checkable in CI.
    """
    out = _python('''
import sys


class _NoRuntime:
    def find_spec(self, name, path=None, target=None):
        if name.split('.')[0] in ('onnxruntime', 'chython_rxnmap'):
            raise ImportError("No module named %r" % name)
        return None


sys.meta_path.insert(0, _NoRuntime())

from chython.core import read_smiles
from chython.reactions.attention._assign import greedy_mapping, side_adjacency
from chython.reactions.attention._encode import encode_reaction

encoded = encode_reaction([read_smiles('CCO')], [read_smiles('CC=O')])
sys.stdout.write('TOKENS\\t%d\\n' % encoded.atoms.shape[0])
''')
    assert out.returncode == 0, f'the encoder needs the runtime to import:\n{out.stderr}'
    reported = dict(line.split('\t') for line in out.stdout.splitlines() if '\t' in line)
    assert reported['TOKENS'] == '9', '1 rxn_cls, then 1 mol_cls + 3 atoms on each of the two sides'


def test_the_container_method_exists_without_the_model():
    """The method is the core's, and only its body comes from `chython.reactions`.

    So the name resolves on any install; what an install without the extra changes is what happens when
    it is called, which the test above pins.
    """
    rxn = ReactionContainer([smiles('CCO')], [smiles('CC=O')])
    assert callable(rxn.attention_mapping)
    with raises(TypeError):
        rxn.attention_mapping(1.75)                    # keyword-only, so a positional is a TypeError
