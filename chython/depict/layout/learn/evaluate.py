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
"""Held-out evaluation: the policy against the model-free chooser on liked drawings no training saw.

| measure | per component |
| --- | --- |
| distance | the aligned distance, in `L`, of the post-stepped drawing to the liked one |
| steps | steps taken before `stop` |
| ct | stated cis/trans units the final plane draws against their parity |
| defects | the scorecard's `ovl`, `crs` and `aob` of the final plane, rescaled |

    python -m chython.depict.layout.learn.evaluate --model DIR [--poses P] [--corpus C]
"""
from argparse import ArgumentParser
from statistics import median

from ....core import smiles as read_smiles
from ..env import Env
from ..greedy import best_step, drive
from ..quality import defects, rescaled, shape
from ..rings import ring_blocks
from .demos import post
from .train import chooser, goal_of, held_out, load


__all__ = ['compare', 'summary', 'main']


def compare(mol, plane, choose):
    """`[(distance, steps, ct, defects)]` per component of `mol`, each driven greedily by `choose` from its tree
    embedding; `choose` None is the model-free chooser on the frontier."""
    rank = mol.canonical_order()
    out = []
    for c in sorted((set(c) for c in mol.connected_components), key=lambda c: min(rank[x] for x in c)):
        env = Env(mol, c, rank, [b for b in ring_blocks(mol) if b[0] <= c], everywhere=choose is not None)
        drive(env, choose or best_step)
        p = post(env)
        goal = goal_of(mol, plane, c)
        d = goal.distance(env)
        sub = mol.substructure(sorted(c))
        q = rescaled(sub, {x: p[x] for x in c})
        bad = defects(sub, q)
        out.append((d, env.steps, shape(sub, q)['ct'], bad['ovl'] + bad['crs'] + bad['aob']))
    return out


def main(argv=None):
    from ..scorecard import CORPUS, POSES, _read_poses, read_corpus
    ap = ArgumentParser(prog='python -m chython.depict.layout.learn.evaluate')
    ap.add_argument('--model', required=True)
    ap.add_argument('--poses', default=str(POSES))
    ap.add_argument('--corpus', default=str(CORPUS))
    a = ap.parse_args(argv)
    choose = chooser(load(a.model))
    rows = {r['id']: r for r in read_corpus(a.corpus)}
    poses = _read_poses(a.poses)
    mine, base = [], []
    for source in sorted(poses):
        if not held_out(source) or source not in rows or '>' in rows[source]['smiles']:
            continue
        mol = read_smiles(rows[source]['smiles'])
        plane = dict(zip(mol, map(tuple, poses[source])))
        mine += compare(mol, plane, choose)
        base += compare(mol, plane, None)
        print(source, [round(r[0], 3) for r in mine[-1:]], [round(r[0], 3) for r in base[-1:]], flush=True)
    for name, rs in (('policy', mine), ('model-free', base)):
        print(summary(name, rs))


def summary(name, rows):
    """One line of a `compare` run: components, median distance, within `SAME`, median steps, ct, defect rows."""
    return (f'{name}: {len(rows)} components, median distance {median(r[0] for r in rows):.3f} L, within 0.15 L '
            f'{sum(r[0] < .15 for r in rows)}, median steps {median(r[1] for r in rows)}, ct {sum(r[2] for r in rows)},'
            f' defect rows {sum(r[3] > 0 for r in rows)}')


if __name__ == '__main__':
    main()
