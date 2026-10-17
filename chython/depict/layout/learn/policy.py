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
"""The policy: a graph transformer over atoms with legal steps as query tokens.

| part | design |
| --- | --- |
| atoms | a linear embedding of the atom features; one global token embedding the global features |
| attention | `LAYERS` pre-norm layers, width `WIDTH`, `HEADS` heads; each pair adds a per-head bias from its hops, \
an RBF of its distance and its bond; a pair not listed (sparse states) is masked; the global token sees and is seen by \
all |
| steps | kind, target type and values embedded to `WIDTH`; each step attends to its member atoms and the global token |
| heads | a logit per step; the value of the state from the global token |
"""
import torch
from torch import nn

from ..edits import KINDS, TARGET_TYPES
from ..features import ATOM, GLOBAL, HOPS, STEP


__all__ = ['Encoder', 'Policy', 'WIDTH', 'HEADS', 'LAYERS']

WIDTH = 32
HEADS = 4
LAYERS = 3
RBF = 16


class _Layer(nn.Module):
    def __init__(self, d, heads):
        super().__init__()
        self.heads = heads
        self.norm1, self.norm2 = nn.LayerNorm(d), nn.LayerNorm(d)
        self.qkv = nn.Linear(d, 3 * d)
        self.out = nn.Linear(d, d)
        self.ffn = nn.Sequential(nn.Linear(d, 2 * d), nn.GELU(), nn.Linear(2 * d, d))

    def forward(self, x, bias):
        n, d = x.shape
        h = self.heads
        q, k, v = self.qkv(self.norm1(x)).reshape(n, 3, h, d // h).permute(1, 2, 0, 3)
        att = (q @ k.transpose(-1, -2)) / (d // h) ** .5 + bias
        x = x + self.out((att.softmax(-1) @ v).permute(1, 0, 2).reshape(n, d))
        return x + self.ffn(self.norm2(x))


class Encoder(nn.Module):
    """Atom tokens and a global token through `layers` attention layers with pair biases: `(n + 1) × width`."""
    def __init__(self, width=WIDTH, heads=HEADS, layers=LAYERS):
        super().__init__()
        self.heads = heads
        self.atom = nn.Linear(ATOM, width)
        self.glob = nn.Linear(GLOBAL, width)
        self.hops = nn.Embedding(HOPS + 1, heads)
        self.register_buffer('centres', torch.linspace(0., 4., RBF))
        self.rbf = nn.Linear(RBF, heads)
        self.bond = nn.Embedding(5, heads)
        self.flags = nn.Linear(2, heads)
        self.layers = nn.ModuleList(_Layer(width, heads) for _ in range(layers))
        self.norm = nn.LayerNorm(width)

    def _bias(self, s, n):
        h = self.heads
        bias = torch.full((h, n + 1, n + 1), float('-inf'), device=s['atoms'].device)
        idx = torch.arange(n + 1, device=bias.device)
        bias[:, idx, idx] = 0.
        bias[:, n, :] = 0.
        bias[:, :, n] = 0.
        if len(s['pairs']):
            rbf = torch.exp(-((s['distance'][:, None] - self.centres) / .25) ** 2)
            b = (self.hops(s['hops']) + self.rbf(rbf) + self.bond(s['bond']) + self.flags(s['flags'])).T
            i, j = s['pairs'][:, 0], s['pairs'][:, 1]
            bias[:, i, j] = b
            bias[:, j, i] = b
        return bias

    def forward(self, s):
        n = s['atoms'].shape[0]
        x = torch.cat([self.atom(s['atoms']), self.glob(s['glob'])[None]])
        bias = self._bias(s, n)
        for layer in self.layers:
            x = layer(x, bias)
        return self.norm(x)


class Policy(nn.Module):
    def __init__(self, width=WIDTH, heads=HEADS, layers=LAYERS):
        super().__init__()
        d = width
        self.size = {'width': width, 'heads': heads, 'layers': layers}
        self.encoder = Encoder(width, heads, layers)
        self.kind = nn.Embedding(len(KINDS) + 1, d)
        self.target = nn.Embedding(len(TARGET_TYPES) + 2, d)
        self.values = nn.Linear(STEP, d)
        self.query, self.key, self.val = nn.Linear(d, d), nn.Linear(d, d), nn.Linear(d, d)
        self.logit = nn.Sequential(nn.Linear(d, d), nn.GELU(), nn.Linear(d, 1))
        self.value = nn.Sequential(nn.Linear(d, d), nn.GELU(), nn.Linear(d, 1))

    def forward(self, s, a):
        """`(logits, value)`: one logit per step of `a`, and the state's value."""
        n = s['atoms'].shape[0]
        x = self.encoder(s)
        t = self.kind(a['kind']) + self.target(a['target']) + self.values(a['values'])
        members = a['members']
        allowed = torch.zeros(len(t), n + 1, dtype=torch.bool, device=x.device)
        rows = torch.arange(len(t), device=x.device)[:, None].expand_as(members)
        keep = members >= 0
        allowed[rows[keep], members[keep]] = True
        allowed[:, n] = True                                # the global token always
        att = (self.query(t) @ self.key(x).T) / x.shape[1] ** .5
        att = att.masked_fill(~allowed, float('-inf')).softmax(-1)
        h = t + att @ self.val(x)
        return self.logit(h).squeeze(-1), self.value(x[n]).squeeze(-1)
