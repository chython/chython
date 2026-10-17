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
"""The policy network and its tensors. Skipped without torch."""
import pytest

torch = pytest.importorskip('torch')

from chython.core import smiles  # noqa: E402
from chython.depict.layout.env import Env  # noqa: E402
from chython.depict.layout.greedy import best_step  # noqa: E402
from chython.depict.layout.learn.policy import Policy  # noqa: E402
from chython.depict.layout.learn.tensors import permuted, state, steps  # noqa: E402
from chython.depict.layout.rings import ring_blocks  # noqa: E402


CHOLESTEROL = 'CC(C)CCCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C'


def _env(text, n=2):
    mol = smiles(text)
    env = Env(mol, set(mol), mol.canonical_order(), ring_blocks(mol), everywhere=True)
    for _ in range(n):
        s = best_step(env, env.legal())
        if s.move.kind == 'stop':
            break
        env.step(s)
    return env


def test_one_logit_per_legal_step_and_a_value():
    env = _env(CHOLESTEROL)
    net = Policy()
    logits, value = net(state(env), steps(env, env.legal()))
    assert logits.shape == (len(env.legal()),) and value.shape == ()
    assert torch.isfinite(logits).all() and torch.isfinite(value)


def test_the_policy_ignores_how_atoms_are_listed():
    torch.manual_seed(0)
    env = _env(CHOLESTEROL)
    net = Policy().eval()
    s, a = state(env), steps(env, env.legal())
    order = torch.randperm(s['atoms'].shape[0])
    ps, pa = permuted(s, a, order)
    l1, v1 = net(s, a)
    l2, v2 = net(ps, pa)
    assert torch.allclose(l1, l2, atol=1e-5) and torch.allclose(v1, v2, atol=1e-5)


def test_stop_alone_still_has_a_logit():
    mol = smiles('CCCCCC')
    env = Env(mol, set(mol), mol.canonical_order(), [])           # no trouble: only `stop`
    logits, _ = Policy()(state(env), steps(env, env.legal()))
    assert logits.shape == (1,)


def test_a_large_record_runs_on_sparse_pairs():
    mol = smiles('C' * 401)
    env = Env(mol, set(mol), mol.canonical_order(), [])          # the frontier's few steps: the pairs are the point
    s = state(env)
    assert s['pairs'].shape[0] < 40 * len(mol)
    logits, _ = Policy()(s, steps(env, env.legal()))
    assert torch.isfinite(logits).all()


@pytest.mark.skipif(not torch.backends.mps.is_available(), reason='no MPS')
def test_mps_agrees_with_cpu():
    torch.manual_seed(0)
    env = _env(CHOLESTEROL)
    net = Policy().eval()
    s, a = state(env), steps(env, env.legal())
    l1, v1 = net(s, a)
    dev = torch.device('mps')
    l2, v2 = net.to(dev)({k: v.to(dev) for k, v in s.items()}, {k: v.to(dev) for k, v in a.items()})
    assert torch.allclose(l1, l2.cpu(), atol=1e-4) and torch.allclose(v1, v2.cpu(), atol=1e-4)
