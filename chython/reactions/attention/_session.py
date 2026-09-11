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
"""The ONNX Runtime session over the `chython-rxnmap` weights.

STDLIB ONLY AT MODULE LEVEL.  `chython.reactions` imports this module at its own init so that
`attention_available()` answers without a heavy import; the runtime and the weights are reached inside
`get_session`, which is the first call that needs either.

The artefact is fixed and carries no schema: one file, one public name, `chython_rxnmap.model_path`.
Its properties are constants in `_encode.py` rather than something negotiated here, because a different
model is a different distribution with its own algorithm, not a second version this code adapts to.
"""
from functools import cache
from importlib.util import find_spec
from os import cpu_count


#: What ONNX Runtime is asked for when the caller names no thread count.  chython 2's value.
THREAD_CEILING = 8


def attention_available() -> bool:
    """Whether the runtime and the weights of `chython[mapping]` are installed.  Loads neither.

    `find_spec` locates a top-level distribution without executing it, so this stays cheap enough to
    call in a loop and cheap enough for a documentation sample's `:skipif:`.  It answers about the
    INSTALLATION and not about a reaction, which is why it is a function here and not a container
    method.

    A PROBE NEVER RAISES.  `find_spec` propagates whatever a finder raises, and a finder that fails --
    a broken egg-link, a `sys.meta_path` entry of someone else's -- is an installation this cannot use;
    `False` is the answer, and the ImportError with the extra's name in it belongs to `get_session`,
    which is where the caller who ignored this went next.
    """
    try:
        return find_spec('onnxruntime') is not None and find_spec('chython_rxnmap') is not None
    except (ImportError, ValueError):
        return False


def default_threads() -> int:
    """The intra-op thread count when the caller names none."""
    return min(cpu_count() or 4, THREAD_CEILING)


@cache
def get_session(threads: int):
    """The loaded model, one session per thread count.

    CACHED ON `threads`, SO A SECOND VALUE COSTS A SECOND LOADED MODEL -- 80.8 MiB of weights read from
    disk again and held for the process's life.  A caller varying the count per reaction pays for it
    once per distinct value, which is the honest price of letting the count be an argument at all.

    Inter-op parallelism is 1: the graph is a single chain of eight encoder layers, so there is nothing
    for a second op thread to run.
    """
    try:
        import onnxruntime as ort
    except ImportError:
        raise ImportError('attention mapping needs ONNX Runtime, which is an extra: '
                          '`pip install chython[mapping]`') from None
    try:
        from chython_rxnmap import model_path
    except ImportError:
        raise ImportError('attention mapping needs the model weights, which are their own '
                          'distribution because they are 80 MiB: `pip install chython[mapping]`') from None

    options = ort.SessionOptions()
    options.inter_op_num_threads = 1
    options.intra_op_num_threads = threads
    options.graph_optimization_level = ort.GraphOptimizationLevel.ORT_ENABLE_ALL
    return ort.InferenceSession(model_path, options, providers=['CPUExecutionProvider'])


__all__ = ['THREAD_CEILING', 'attention_available', 'default_threads', 'get_session']
