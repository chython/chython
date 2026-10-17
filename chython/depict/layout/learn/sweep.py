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
"""A size sweep: each policy size trained to its validation plateau on one demonstration run, then evaluated on the
held-out poses against the model-free chooser.

    python -m chython.depict.layout.learn.sweep --demos RUN --out DIR [--sizes 32x3x4 64x3x4 ...]

| part | rule |
| --- | --- |
| size | `WIDTHxLAYERSxHEADS` |
| run | `train.main` per size into `DIR/<size>`, one after another, never two at once |
| report | `DIR/sweep.md`: per size its parameters, best validation loss, epochs, and the held-out `summary` |
"""
from argparse import ArgumentParser
from pathlib import Path

from ....core import smiles as read_smiles
from .evaluate import compare, summary
from .train import chooser, held_out, load, main as train_main


__all__ = ['SIZES', 'main']

SIZES = ('32x3x4', '64x3x4', '64x4x8', '128x4x8')


def _held(poses, rows):
    out = []
    for source in sorted(poses):
        if held_out(source) and source in rows and '>' not in rows[source]['smiles']:
            mol = read_smiles(rows[source]['smiles'])
            out.append((mol, dict(zip(mol, map(tuple, poses[source])))))
    return out


def main(argv=None):
    import json
    from ..scorecard import CORPUS, POSES, _read_poses, read_corpus
    ap = ArgumentParser(prog='python -m chython.depict.layout.learn.sweep')
    ap.add_argument('--demos', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--sizes', nargs='+', default=list(SIZES))
    ap.add_argument('--seed', type=int, default=0)
    a = ap.parse_args(argv)
    out = Path(a.out).expanduser()
    out.mkdir(parents=True, exist_ok=True)
    held = _held(_read_poses(POSES), {r['id']: r for r in read_corpus(CORPUS)})
    base = [r for mol, plane in held for r in compare(mol, plane, None)]
    lines = ['| size | parameters | best valid | epochs | held out |', '| --- | --- | --- | --- | --- |']
    for size in a.sizes:
        w, l, h = map(int, size.split('x'))
        where = out / size
        train_main(['--demos', a.demos, '--out', str(where), '--width', str(w), '--layers', str(l), '--heads', str(h),
                    '--seed', str(a.seed)])
        meta = json.loads((where / 'manifest.json').read_text())
        choose = chooser(load(where))
        mine = [r for mol, plane in held for r in compare(mol, plane, choose)]
        lines.append(f"| {size} | {meta['parameters']} | {meta['history']['best']:.4f} | "
                     f"{len(meta['history']['train'])} | {summary('policy', mine)} |")
        (out / 'sweep.md').write_text('\n'.join(lines + ['', summary('model-free', base)]) + '\n')
        print(lines[-1], flush=True)


if __name__ == '__main__':
    main()
