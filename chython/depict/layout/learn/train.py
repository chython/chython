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
"""N2: the policy trained toward liked drawings by outcome, then refined on its own routes.

| stage | rule |
| --- | --- |
| examples | each demonstration that reached its goal within `SAME`, replayed state by state: its legal steps, \
their `target_set`, and the cost to go, `LAMBDA` per step left plus the final distance |
| imitation | loss −log of the policy's probability on the target set, plus the value's squared error to the cost \
to go |
| refinement | routes sampled from the policy from each training start, each costing `cost` of where it ends: one \
outcome by two routes costs the same, and the shorter route less; policy gradient on the value's advantage |
| split | a source is held out when the first byte of its SHA-1 is under `HELD` × 256 |
| chooser | `chooser(model)` is a `choose(env, steps)`: the step of highest logit, ties by `env.tie` |

    python -m chython.depict.layout.learn.train --demos RUN --out DIR [--epochs N] [--rollouts N] [--seed S]
"""
import json
import random
from argparse import ArgumentParser
from hashlib import sha1
from pathlib import Path

import torch

from ....core import smiles as read_smiles
from ..env import Env
from ..rings import ring_blocks
from .demos import images
from .policy import Policy
from .replay import key, read
from .targets import SAME, Goal, cost, target_set
from .tensors import state, steps


__all__ = ['HELD', 'goal_of', 'held_out', 'split', 'examples', 'imitate', 'refine', 'load', 'chooser', 'main']

HELD = .1
PATIENCE = 30            # epochs without a validation gain before training stops
GAIN = 1e-3             # the least validation gain that counts


def held_out(source):
    return sha1(source.encode()).digest()[0] < HELD * 256


def goal_of(mol, plane, atoms):
    """The `Goal` of the component `atoms` of `mol` drawn as `plane`."""
    rank = mol.canonical_order()
    order = sorted(atoms, key=rank.__getitem__)
    maps = images(mol)
    return Goal(order, [[plane[m[x]] for x in order] for m in maps])


def _component(mol, ranks):
    rank = mol.canonical_order()
    return {x for x in mol if rank[x] in set(ranks)}


def examples(trajectory, plane):
    """The training examples of one demonstration `trajectory` toward the liked `plane`."""
    mol = read_smiles(trajectory.smiles)
    atoms = _component(mol, trajectory.meta['component'])
    goal = goal_of(mol, plane, atoms)
    rank = mol.canonical_order()
    env = Env(mol, atoms, rank, [b for b in ring_blocks(mol) if b[0] <= atoms], everywhere=True)
    final = trajectory.meta['distance']
    out = []
    keys = list(trajectory.keys) + [('stop',)]
    for i, k in enumerate(keys):
        legal = env.legal()
        mask, _ = target_set(env, legal, goal)
        out.append((state(env), steps(env, legal), torch.tensor(mask), cost(final, len(keys) - 1 - i)))
        found = [s for s in legal if key(env, s) == k]
        if len(found) != 1:
            break
        env.step(found[0])
    return out


def split(sources):
    """`(train, valid)` of the sources not held out: valid when the second byte of the SHA-1 is under `HELD` × 256."""
    train, valid = [], []
    for x in sources:
        (valid if sha1(x.encode()).digest()[1] < HELD * 256 else train).append(x)
    return train, valid


def _loss(model, example, device):
    s, a, mask, c = example
    logits, value = model({k: v.to(device) for k, v in s.items()}, {k: v.to(device) for k, v in a.items()})
    mask = mask.to(device)
    return torch.logsumexp(logits, 0) - torch.logsumexp(logits[mask], 0) + .5 * (value - c) ** 2


def imitate(model, data, epochs, lr=1e-3, seed=0, device='cpu', valid=None, patience=PATIENCE):
    """Train `model` on `data`, `(state, steps, mask, cost)` examples, until `epochs` or until the loss on `valid`
    (default `data`) has not improved by `GAIN` for `patience` epochs; the learning rate halves after `patience` // 3
    stalled epochs, and the best weights are kept. Returns `{'train', 'valid', 'best', 'epoch'}`."""
    torch.manual_seed(seed)
    model.to(device).train()
    opt = torch.optim.Adam(model.parameters(), lr=lr)
    sched = torch.optim.lr_scheduler.ReduceLROnPlateau(opt, factor=.5, patience=max(1, patience // 3))
    valid = data if valid is None else valid
    order = list(range(len(data)))
    rng = random.Random(seed)
    history = {'train': [], 'valid': [], 'best': None, 'epoch': 0}
    best_state, stalled = None, 0
    for epoch in range(epochs):
        model.train()
        rng.shuffle(order)
        total = 0.
        for i in order:
            loss = _loss(model, data[i], device)
            opt.zero_grad()
            loss.backward()
            opt.step()
            total += float(loss)
        model.eval()
        with torch.no_grad():
            v = sum(float(_loss(model, e, device)) for e in valid) / max(1, len(valid))
        history['train'].append(total / max(1, len(data)))
        history['valid'].append(v)
        sched.step(v)
        if history['best'] is None or v < history['best'] - GAIN:
            history['best'], history['epoch'], stalled = v, epoch, 0
            best_state = {k: t.detach().clone() for k, t in model.state_dict().items()}
        else:
            stalled += 1
            if stalled >= patience:
                break
    if best_state is not None:
        model.load_state_dict(best_state)
        history['best'] = min(history['valid'])
    return history


def refine(model, starts, rollouts, lr=3e-4, seed=0, device='cpu', entropy=.01):
    """Policy-gradient refinement: from each `(env factory, goal)` of `starts`, `rollouts` sampled routes, each
    rewarded −`cost` of its end; returns the mean reward per pass."""
    torch.manual_seed(seed)
    model.to(device).train()
    opt = torch.optim.Adam(model.parameters(), lr=lr)
    means = []
    for _ in range(rollouts):
        rewards = []
        for make, goal in starts:
            env = make()
            logps, values, ents = [], [], []
            while not env.done:
                legal = env.legal()
                s = {k: v.to(device) for k, v in state(env).items()}
                a = {k: v.to(device) for k, v in steps(env, legal).items()}
                logits, value = model(s, a)
                p = torch.distributions.Categorical(logits=logits)
                i = p.sample()
                logps.append(p.log_prob(i))
                ents.append(p.entropy())
                values.append(value)
                env.step(legal[int(i)])
            c = cost(goal.distance(env), env.steps, env.terms[0])     # one outcome, one cost, whatever the route
            rewards.append(-c)
            values = torch.stack(values)
            advantage = (values - c).detach()               # better than the value expected: reinforce
            loss = -(torch.stack(logps) * advantage).sum() + .5 * ((values - c) ** 2).sum() \
                - entropy * torch.stack(ents).sum()
            opt.zero_grad()
            loss.backward()
            opt.step()
        means.append(sum(rewards) / max(1, len(rewards)))
    return means


def load(directory):
    """The policy in `directory`: `policy.pt`, sized by `manifest.json` when there is one, whose `kind` says which
    environment it drives: `primitive`, or the edit environment by default."""
    directory = Path(directory).expanduser()
    meta = json.loads((directory / 'manifest.json').read_text()) if (directory / 'manifest.json').exists() else {}
    model = Policy(**meta.get('size', {}))
    model.load_state_dict(torch.load(directory / 'policy.pt'))
    model.environment = meta.get('kind')
    return model


def chooser(model, device='cpu'):
    """A `choose(env, steps)` taking the step of highest logit, ties by `env.tie`."""
    model.to(device).eval()

    def choose(env, legal):
        with torch.no_grad():
            logits, _ = model({k: v.to(device) for k, v in state(env).items()},
                              {k: v.to(device) for k, v in steps(env, legal).items()})
        top = float(logits.max())
        return min((s for s, x in zip(legal, logits.tolist()) if x >= top - 1e-6), key=lambda s: env.tie(s.move))
    choose.everywhere = True                                 # a policy edits where nothing is wrong, too
    choose.primitive = getattr(model, 'environment', None) == 'primitive'
    return choose


def _dataset(demos, poses, rebuild=False):
    """`{source: examples}` of the demonstrations in `demos` that reached their goal, not held out; cached in
    `demos/examples.pt`."""
    path = Path(demos) / 'examples.pt'
    if path.exists() and not rebuild:
        return torch.load(path, weights_only=False)
    out = {}
    for t in read(demos):
        if t.meta['distance'] >= SAME or held_out(t.source) or t.source not in poses:
            continue
        mol = read_smiles(t.smiles)
        out.setdefault(t.source, []).extend(examples(t, dict(zip(mol, map(tuple, poses[t.source])))))
    torch.save(out, path)
    return out


def main(argv=None):
    from ..scorecard import POSES, _read_poses
    ap = ArgumentParser(prog='python -m chython.depict.layout.learn.train')
    ap.add_argument('--demos', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--poses', default=str(POSES))
    ap.add_argument('--epochs', type=int, default=500)
    ap.add_argument('--patience', type=int, default=PATIENCE)
    ap.add_argument('--width', type=int, default=32)
    ap.add_argument('--layers', type=int, default=3)
    ap.add_argument('--heads', type=int, default=4)
    ap.add_argument('--lr', type=float, default=1e-3)
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--device', default='cpu')
    ap.add_argument('--rebuild', action='store_true')
    a = ap.parse_args(argv)
    data = _dataset(a.demos, _read_poses(a.poses), a.rebuild)
    train, valid = split(sorted(data))
    model = Policy(a.width, a.heads, a.layers)
    history = imitate(model, [e for x in train for e in data[x]], a.epochs, a.lr, a.seed, a.device,
                      [e for x in valid for e in data[x]], a.patience)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    torch.save(model.state_dict(), out / 'policy.pt')
    params = sum(p.numel() for p in model.parameters())
    (out / 'manifest.json').write_text(json.dumps({'demos': a.demos, 'size': model.size, 'parameters': params,
                                                    'train': len(train), 'valid': len(valid), 'lr': a.lr,
                                                    'seed': a.seed, 'history': history}, indent=1))
    print(f"{model.size} {params} params: {len(train)}+{len(valid)} sources, best valid {history['best']:.4f} at "
          f"epoch {history['epoch']} of {len(history['train'])}")

if __name__ == '__main__':
    main()
