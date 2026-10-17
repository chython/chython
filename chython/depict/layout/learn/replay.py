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
"""Replay: a state is its molecule and the keys of the steps taken, so a trajectory is a few hundred bytes and any of
its states rebuilds exactly.

| part | rule |
| --- | --- |
| key | `(kind, system, variant)` for an embed, `('stop',)`, `('pin', pins, let-go pairs)`, else `(kind, rank of \
pivot, rank of other, argument to 1e-9, sorted ranks of the atoms moved)`: a function of the molecule, not of ids |
| replay | each key must name exactly one legal step; anything else raises `ReplayMismatch` |
| store | append-only JSONL, one file per writing process, `replay-<pid>.jsonl` in the run directory |
| internal | a trajectory of an internal molecule holds its hash as `source` and no SMILES: it is re-read on rebuild |
"""
import json
import os
from pathlib import Path
from typing import NamedTuple

from ....core import smiles as read_smiles
from ..env import Env
from ..rings import ring_blocks


__all__ = ['Trajectory', 'ReplayMismatch', 'key', 'replay', 'write', 'read', 'rebuild']


class ReplayMismatch(ValueError):
    """Raised when a key names no legal step, or more than one."""
    def __init__(self, index, key, found):
        super().__init__(f'step {index}: key {key!r} names {found} legal steps, not one')
        self.index = index
        self.key = key


class Trajectory(NamedTuple):
    source: str             # a public id, or `h:` and a hash for an internal molecule
    smiles: str             # None for an internal molecule
    keys: list              # the step keys, in order
    meta: dict              # `everywhere` and whatever the writer records: seeds, final distance


def key(env, step):
    """The replay key of a legal `step` of `env`."""
    m = step.move
    if m.kind == 'stop':
        return ('stop',)
    if m.kind == 'embed':
        return ('embed', m.pivot, m.other)
    if m.kind == 'pin':
        return ('pin', tuple(sorted(step.state['pins'].items())), tuple(sorted(step.state['relax'])))
    rank = env.scope.rank
    return (m.kind, rank[m.pivot], rank[m.other], round(m.arg, 9), tuple(sorted(rank[x] for x in m.part)))


def replay(mol, keys, everywhere=False, component=None):
    """The `Env` of `mol` after the steps `keys` name: of the component whose atoms have the canonical ranks
    `component`, or of every atom as one."""
    rank = mol.canonical_order()
    atoms = set(mol) if component is None else {x for x in mol if rank[x] in set(component)}
    env = Env(mol, atoms, rank, [b for b in ring_blocks(mol) if b[0] <= atoms], everywhere=everywhere)
    for i, k in enumerate(keys):
        found = [s for s in env.legal() if key(env, s) == k]
        if len(found) != 1:
            raise ReplayMismatch(i, k, len(found))
        env.step(found[0])
    return env


def _tuples(x):
    return tuple(_tuples(y) for y in x) if isinstance(x, list) else x


def write(run_dir, trajectories):
    """Append `trajectories` to this process's file in `run_dir`."""
    path = Path(run_dir) / f'replay-{os.getpid()}.jsonl'
    with open(path, 'a') as f:
        for t in trajectories:
            f.write(json.dumps({'source': t.source, 'smiles': t.smiles, 'keys': t.keys, 'meta': t.meta}) + '\n')


def read(run_dir):
    """Every trajectory in `run_dir`, file by file in name order, each file in the order written."""
    out = []
    for path in sorted(Path(run_dir).glob('replay-*.jsonl')):
        for line in path.read_text().splitlines():
            r = json.loads(line)
            out.append(Trajectory(r['source'], r['smiles'], [_tuples(k) for k in r['keys']], r['meta']))
    return out


def rebuild(trajectory, smiles_of):
    """The `Env` before each step of `trajectory` and after its last; `smiles_of(source)` gives an internal
    molecule's SMILES, read again at the source."""
    text = trajectory.smiles if trajectory.smiles is not None else smiles_of(trajectory.source)
    mol = read_smiles(text)
    everywhere = trajectory.meta.get('everywhere', False)
    component = trajectory.meta.get('component')
    return [replay(mol, trajectory.keys[:i], everywhere, component) for i in range(len(trajectory.keys) + 1)]
