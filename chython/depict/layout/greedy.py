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
"""Greedy driving of an environment, and the model-free chooser.

| part | rule |
| --- | --- |
| `drive` | until `done`: `choose(env, env.legal())` names one legal step, which is taken |
| `best_step` | the model-free chooser: the step of least rank, ties by the environment's `tie`; `stop` once no step \
betters the state |
"""
__all__ = ['drive', 'best_step', 'peptide_greedy']


def drive(env, choose):
    """Take the steps `choose` picks until the environment is done; returns it."""
    while not env.done:
        env.step(choose(env, env.legal()))
    return env


def best_step(env, steps):
    """The step of least rank, ties by `env.tie`; `stop` when none betters the state. Embedding is never stopped."""
    stop = next((s for s in steps if s.move.kind == 'stop'), None)
    moves = [s for s in steps if s.move.kind != 'stop']
    if not moves:
        return stop
    best = min(moves, key=lambda s: (s.rank, env.tie(s.move)))
    if stop is not None and best.rank >= env.rank:
        return stop
    return best


def peptide_greedy(choose=best_step):
    """A `peptide_layout` chooser driving a `PeptideEnv` greedily with `choose`."""
    from .env import PeptideEnv

    def run(ctx):
        env = drive(PeptideEnv(ctx), choose)
        return None if env.plane is None else env.inv
    return run
