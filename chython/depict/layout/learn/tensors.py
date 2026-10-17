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
"""Torch tensors of an environment's state and legal steps, from `layout.features`.

| tensor | shape | content |
| --- | --- | --- |
| `atoms` | n × `ATOM` | atom features, atoms in canonical rank order |
| `pairs` | P × 2 | position pairs `i < j` |
| `hops`, `bond` | P | integers: bonds between, bond kind |
| `distance` | P | in `L` |
| `flags` | P × 2 | ring bond; crosses another bond |
| `glob` | `GLOBAL` | global features |
| `kind`, `target` | S | per legal step |
| `values` | S × `STEP` | per legal step |
| `members` | S × A | atom positions a step reads, its pivot and other atom first; -1 pads |
"""
import torch

from ..features import atom_features, global_features, pair_features, step_features


__all__ = ['state', 'steps', 'permuted']


def state(env):
    """The tensors of `env`'s state."""
    pairs, feats = pair_features(env)
    return {'atoms': torch.tensor(atom_features(env), dtype=torch.float32),
            'pairs': torch.tensor(pairs, dtype=torch.long).reshape(-1, 2),
            'hops': torch.tensor([f[0] for f in feats], dtype=torch.long),
            'distance': torch.tensor([f[1] for f in feats], dtype=torch.float32),
            'bond': torch.tensor([f[2] for f in feats], dtype=torch.long),
            'flags': torch.tensor([f[3:] for f in feats], dtype=torch.float32).reshape(-1, 2),
            'glob': torch.tensor(global_features(env), dtype=torch.float32)}


def steps(env, legal):
    """The tensors of `legal`, steps of `env`."""
    rows = [step_features(env, s) for s in legal]
    width = max([1] + [len(r[3]) for r in rows])
    members = torch.full((len(rows), width), -1, dtype=torch.long)
    for i, r in enumerate(rows):
        if r[3]:
            members[i, :len(r[3])] = torch.tensor(r[3], dtype=torch.long)
    return {'kind': torch.tensor([r[0] for r in rows], dtype=torch.long),
            'target': torch.tensor([r[1] for r in rows], dtype=torch.long),
            'values': torch.tensor([r[2] for r in rows], dtype=torch.float32),
            'members': members}


def permuted(s, a, order):
    """`s` and `a` with atoms listed as `order` (new position -> old position): for testing equivariance."""
    inv = torch.empty_like(order)
    inv[order] = torch.arange(len(order))
    pairs = inv[s['pairs']]
    pairs, _ = pairs.sort(dim=1)
    members = a['members'].clone()
    keep = members >= 0
    members[keep] = inv[members[keep]]
    return {**s, 'atoms': s['atoms'][order], 'pairs': pairs}, {**a, 'members': members}
