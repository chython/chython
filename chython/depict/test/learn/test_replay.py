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
"""Replay keys, exact replay and the replay store."""
import pytest

from chython.core import smiles
from chython.depict.layout.env import Env
from chython.depict.layout.greedy import best_step
from chython.depict.layout.learn.replay import ReplayMismatch, Trajectory, key, read, rebuild, replay, write
from chython.depict.layout.rings import ring_blocks


CHOLESTEROL = 'CC(C)CCCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C'


def _run(mol, steps=10, everywhere=True):
    env = Env(mol, set(mol), mol.canonical_order(), ring_blocks(mol), everywhere=everywhere)
    keys = []
    while not env.done and len(keys) < steps:
        s = best_step(env, env.legal())
        keys.append(key(env, s))
        env.step(s)
    return env, keys


def _canonical(env):
    rank = env.scope.rank
    return sorted((rank[x], p) for x, p in env.plane.items())


def test_a_run_replays_bit_for_bit_from_any_atom_order():
    mol = smiles(CHOLESTEROL)
    env, keys = _run(mol)
    ref = _canonical(env)
    for _ in range(12):
        other = smiles(format(mol, 'r'))
        assert _canonical(replay(other, keys, everywhere=True)) == ref


def test_an_edited_key_does_not_replay():
    mol = smiles(CHOLESTEROL)
    _, keys = _run(mol)
    assert len(keys) > 1
    bad = list(keys)
    bad[0] = (bad[0][0], 10 ** 6) + tuple(bad[0][2:])
    with pytest.raises(ReplayMismatch) as info:
        replay(mol, bad, everywhere=True)
    assert info.value.index == 0


def test_the_store_round_trips_and_keeps_internal_structures_out(tmp_path):
    mol = smiles(CHOLESTEROL)
    _, keys = _run(mol)
    public = Trajectory('pubchem:5997', CHOLESTEROL, keys, {'everywhere': True})
    internal = Trajectory('h:0123456789ab', None, keys, {'everywhere': True})
    write(tmp_path, [public, internal])
    back = read(tmp_path)
    assert back == [public, internal]
    text = ''.join(p.read_text() for p in tmp_path.iterdir())
    assert CHOLESTEROL in text and text.count('"smiles": null') == 1
    envs = rebuild(back[1], lambda source: CHOLESTEROL)
    assert len(envs) == len(keys) + 1 and _canonical(envs[-1]) == _canonical(replay(mol, keys, everywhere=True))
