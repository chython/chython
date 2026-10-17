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
"""The judge: a network mimicking the user's taste, the probability-like score that a finished drawing is one the user
would accept.

| part | rule |
| --- | --- |
| network | the policy's `Encoder` and a scalar head on its global token |
| view | a drawing of one component read as a `PrimitiveEnv` at that plane: the policy's own features |
| pairs | within one component: a liked drawing over each negative |
| negatives | the user's disliked drawings, the policy's own drawings, and `perturbed` liked drawings: a few \
random primitive steps away, so a step off the liked drawing is worse |
| loss | Bradley–Terry: softplus of the negative's score minus the positive's |
| training | to a plateau of the validation loss, as `train.imitate` does; accuracy is the share of pairs \
ordered right |
"""
import random

import torch
from torch import nn

from ..relax import _ct_broken
from ..score import terms
from .policy import HEADS, LAYERS, WIDTH, Encoder
from .tensors import state


__all__ = ['Judge', 'view', 'perturbed', 'pairs_of', 'train_judge', 'judge_score']

PATIENCE = 30
GAIN = 1e-3


class Judge(nn.Module):
    def __init__(self, width=WIDTH, heads=HEADS, layers=LAYERS):
        super().__init__()
        self.size = {'width': width, 'heads': heads, 'layers': layers}
        self.encoder = Encoder(width, heads, layers)
        self.head = nn.Sequential(nn.Linear(width, width), nn.GELU(), nn.Linear(width, 1))

    def forward(self, s):
        x = self.encoder(s)
        return self.head(x[-1]).squeeze(-1)


def view(initial, plane):
    """`initial`'s environment drawn as `plane`."""
    env = initial.env()
    env.plane = {x: tuple(plane[x]) for x in env.plane}
    env.terms = terms(env.scope, env.plane)
    env.ct = _ct_broken(env.mol, env.plane, env.scope.frames, set(env.plane))
    return env


def perturbed(initial, plane, steps=3, seed=0):
    """`plane` after `steps` random primitive steps, whole-drawing turns excepted: the judge sees no orientation."""
    rng = random.Random(seed)
    env = view(initial, plane)
    for _ in range(steps):
        legal = [s for s in env.legal() if s.move.kind not in ('stop', 'turn')]
        if not legal:
            break
        env.step(rng.choice(legal))
    return dict(env.plane)


def pairs_of(initial, liked, negatives):
    """`[(liked state, negative state)]` tensors for one component."""
    good = state(view(initial, liked))
    return [(good, state(view(initial, n))) for n in negatives]


def _scores(model, pair, device):
    a, b = pair
    return model({k: v.to(device) for k, v in a.items()}), model({k: v.to(device) for k, v in b.items()})


def train_judge(model, data, epochs, lr=1e-3, seed=0, device='cpu', valid=None, patience=PATIENCE):
    """Train on `data`, `(liked, negative)` state pairs, until `epochs` or a plateau of the loss on `valid`; the best
    weights kept. Returns `{'train', 'valid', 'accuracy', 'best', 'epoch'}`."""
    torch.manual_seed(seed)
    model.to(device).train()
    opt = torch.optim.Adam(model.parameters(), lr=lr)
    sched = torch.optim.lr_scheduler.ReduceLROnPlateau(opt, factor=.5, patience=max(1, patience // 3))
    valid = data if valid is None else valid
    order = list(range(len(data)))
    rng = random.Random(seed)
    history = {'train': [], 'valid': [], 'accuracy': [], 'best': None, 'epoch': 0}
    best_state, stalled = None, 0
    for epoch in range(epochs):
        model.train()
        rng.shuffle(order)
        total = 0.
        for i in order:
            good, bad = _scores(model, data[i], device)
            loss = nn.functional.softplus(bad - good)
            opt.zero_grad()
            loss.backward()
            opt.step()
            total += float(loss)
        model.eval()
        with torch.no_grad():
            v = right = 0.
            for p in valid:
                good, bad = _scores(model, p, device)
                v += float(nn.functional.softplus(bad - good))
                right += float(good > bad)
        n = max(1, len(valid))
        history['train'].append(total / max(1, len(data)))
        history['valid'].append(v / n)
        history['accuracy'].append(right / n)
        sched.step(v / n)
        if history['best'] is None or v / n < history['best'] - GAIN:
            history['best'], history['epoch'], stalled = v / n, epoch, 0
            best_state = {k: t.detach().clone() for k, t in model.state_dict().items()}
        else:
            stalled += 1
            if stalled >= patience:
                break
    if best_state is not None:
        model.load_state_dict(best_state)
    return history


def judge_score(model, initial, plane, device='cpu'):
    """The judge's score of `initial`'s component drawn as `plane`."""
    model.to(device).eval()
    with torch.no_grad():
        return float(model({k: v.to(device) for k, v in state(view(initial, plane)).items()}))
