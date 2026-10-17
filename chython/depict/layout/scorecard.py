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
"""The clean2d scorecard: every engine over `test/layout_corpus.tsv`, one subprocess each.

    python -m chython.depict.layout.scorecard run [--engines E ...] [--corpus P] [--out DIR]
    python -m chython.depict.layout.scorecard report [--out DIR]
    python -m chython.depict.layout.scorecard pose ID ENGINE [--out DIR]

A worker appends one JSON line per row and flushes it, so a crash loses at most the row it died on.  The
runner records that row as failed and restarts the worker after it.  An engine whose toolkit is not
installed is recorded absent, which is not a failure.

`peptide` is the backbone layout over the default engine; `chython` is `layout2d` as it ships with the chython
engine, the backbone layout included; `greedy` is the controller driven greedily by its model-free chooser, the
baseline a trained policy replaces; `policy:DIR` is the trained policy in `DIR/policy.pt` driving the same path.
"""
import json
import os
import subprocess
import sys
import time
from math import isfinite
from pathlib import Path
from statistics import fmean, median, quantiles

from ...core import ReactionContainer, smiles
from ..peptide import routed_links
from .molecule import layout2d
from .quality import defects, perm, pose, rescaled, rxn, ser, shape


ENGINES = ('smilesdrawer', 'rdkit', 'cdk', 'obabel', 'indigo', 'peptide', 'chython', 'greedy')
ROOT = Path(__file__).resolve().parents[3]
CORPUS = ROOT / 'test' / 'layout_corpus.tsv'
THREADS = 2              # torch threads of a policy engine


def read_corpus(path):
    lines = Path(path).read_text().splitlines()
    keys = lines[0].split('\t')
    return [dict(zip(keys, s.split('\t'))) for s in lines[1:] if s]


def layout_fn(engine):
    if engine == 'peptide':
        return lambda mol: layout2d(mol, force=True)
    if engine == 'chython':
        return lambda mol: layout2d(mol, engine=engine, force=True, peptide=True)
    if engine == 'greedy':
        from .greedy import best_step
        return lambda mol: layout2d(mol, engine='chython', force=True, peptide=True, choose=best_step)
    if engine.startswith('policy:'):
        import torch
        from .learn.train import chooser, load
        torch.set_num_threads(THREADS)                      # a forward pass this small gains nothing from more
        choose = chooser(load(engine[7:]))
        return lambda mol: layout2d(mol, engine='chython', force=True, peptide=True, choose=choose)
    return lambda mol: layout2d(mol, engine=engine, force=True, peptide=False)


def _valid(mol, plane):
    return set(plane) == set(mol) and all(isfinite(v) for p in plane.values() for v in p)


def score_row(row, lay):
    rec = {'id': row['id'], 'stratum': row['stratum'], 'ok': False, 'error': None, 't_us': None,
           'plane': None, 'metrics': None}
    try:
        obj = smiles(row['smiles'])
        mols = list(obj.molecules()) if isinstance(obj, ReactionContainer) else [obj]
        t = time.perf_counter()
        planes = [lay(m) for m in mols]
        rec['t_us'] = (time.perf_counter() - t) * 1e6
        if not all(_valid(m, p) for m, p in zip(mols, planes)):
            rec['error'] = 'plane misses an atom or holds a non-finite coordinate'
            return rec
        links = [[(k.a, k.b) for k in routed_links(m, p)] for m, p in zip(mols, planes)]
        planes = [p if k else rescaled(m, p) for m, p, k in zip(mols, planes, links)]   # routed: the lattice's L
        rec['plane'] = [[list(p[x]) for x in m] for m, p in zip(mols, planes)]
        if isinstance(obj, ReactionContainer):
            rec['metrics'] = {'rxn': rxn(obj, planes)}
        else:
            (mol,), (plane,), (routed,) = mols, planes, links
            rec['metrics'] = {**defects(mol, plane, routed=routed), **shape(mol, plane), 'perm': perm(mol, lay)}
        rec['ok'] = True
    except Exception as e:                                  # an engine's failure is a result, not a crash
        rec['error'] = f'{type(e).__name__}: {e}'
    return rec


def worker(engine, corpus, out, after=None, lay=None):
    rows = read_corpus(corpus)
    if after is not None:
        rows = rows[[r['id'] for r in rows].index(after) + 1:]
    lay = lay or layout_fn(engine)
    with open(out, 'a') as f:
        try:
            lay(smiles('CC'))                               # warm up, and learn whether the toolkit is here
        except Exception as e:                              # whatever a missing toolkit raises: a JVM, a jar
            f.write(json.dumps({'absent': f'{type(e).__name__}: {e}'}) + '\n')
            return
        for row in rows:
            f.write(json.dumps(score_row(row, lay)) + '\n')
            f.flush()


def _rss_mb(pid):
    r = subprocess.run(['ps', '-o', 'rss=', '-p', str(pid)], capture_output=True, text=True)
    return int(r.stdout.strip() or 0) / 1024


def _last_id(path):
    if not path.exists():
        return None
    ids = [json.loads(s).get('id') for s in path.read_text().splitlines() if s]
    return ids[-1] if ids else None


def run(engines, corpus, outdir, rss_mb=4096, stall_s=120, cmd=None):
    """One engine at a time: a subprocess killed for a crash, an RSS above `rss_mb` or `stall_s` without a
    line records the row after the last written one as failed, then resumes past it."""
    cmd = cmd or [sys.executable, '-m', 'chython.depict.layout.scorecard', 'worker']
    env = dict(os.environ)                                  # the worker imports this tree, not an installed copy
    env['PYTHONPATH'] = os.pathsep.join(filter(None, (str(ROOT), env.get('PYTHONPATH'))))
    ids = [r['id'] for r in read_corpus(corpus)]
    rows = {r['id']: r for r in read_corpus(corpus)}
    for engine in engines:
        out = Path(outdir) / f'{engine}.jsonl'
        out.unlink(missing_ok=True)
        after = None
        while True:
            args = [*cmd, engine, str(corpus), str(out)] + ([after] if after else [])
            proc = subprocess.Popen(args, env=env)
            size, seen = -1, time.monotonic()
            reason = None
            while proc.poll() is None:
                time.sleep(.5)
                now = out.stat().st_size if out.exists() else 0
                if now != size:
                    size, seen = now, time.monotonic()
                if _rss_mb(proc.pid) > rss_mb:
                    reason = f'RSS above {rss_mb} MB'
                elif time.monotonic() - seen > stall_s:
                    reason = f'no row in {stall_s} s'
                if reason:
                    proc.kill()
                    proc.wait()
                    break
            if proc.returncode == 0 and reason is None:
                break
            last = _last_id(out)
            k = ids.index(last) + 1 if last in ids else 0
            if k >= len(ids):
                break
            dead = rows[ids[k]]
            with open(out, 'a') as f:
                f.write(json.dumps({'id': dead['id'], 'stratum': dead['stratum'], 'ok': False, 't_us': None,
                                    'error': reason or f'worker exited {proc.returncode}', 'plane': None,
                                    'metrics': None}) + '\n')
            if k + 1 >= len(ids):
                break
            after = dead['id']


POSES = CORPUS.with_name('layout_poses.tsv')
COUNTS = ('ovl', 'ovl61', 'crs', 'aob', 'lbl', 'sqz', 'ct', 'inr', 'perm')
MEANS = ('cv', 'rng', 'ang', 'axis', 'rxn')


def _load(outdir):
    out = {}
    for p in sorted(Path(outdir).glob('*.jsonl')):
        recs = [json.loads(s) for s in p.read_text().splitlines() if s]
        out[p.stem] = 'absent' if recs and 'absent' in recs[0] else {r['id']: r for r in recs}
    return out


def _read_poses(poses):
    if poses is None or not Path(poses).exists():
        return {}
    return {i: json.loads(p) for i, p in (s.split('\t') for s in Path(poses).read_text().splitlines() if s)}


def _mol(row):
    return smiles(row['smiles'])


def summarize(outdir, corpus, poses=None):
    rows = read_corpus(corpus)
    refs = _read_poses(poses)
    table = {}
    for engine, recs in _load(outdir).items():
        if recs == 'absent':
            table[engine] = 'absent'
            continue
        table[engine] = per = {}
        for stratum in dict.fromkeys(r['stratum'] for r in rows):
            mine = [r for r in rows if r['stratum'] == stratum and r['id'] in recs]
            ok = [recs[r['id']] for r in mine if recs[r['id']]['ok']]
            ts = sorted(r['t_us'] for r in ok if r['t_us'] is not None)
            s = {'n': len(mine), 'failed': len(mine) - len(ok)}
            for k in COUNTS:
                vs = [r['metrics'][k] for r in ok if k in r['metrics']]
                s[k] = 100 * sum(v > 0 for v in vs) / len(vs) if vs else None
            for k in MEANS:
                vs = [r['metrics'][k] for r in ok if r['metrics'].get(k) is not None]
                s[k] = fmean(vs) if vs else None
            s['t50'] = median(ts) if ts else None
            s['t95'] = quantiles(ts, n=20)[-1] if len(ts) > 1 else (ts[0] if ts else None)
            series = {}
            for r in mine:
                if r['stratum'] == 'ser' and recs[r['id']]['ok']:
                    series.setdefault(r['id'].rsplit(':', 1)[0], []).append(r)
            vs = []
            for members in series.values():
                ms = [(_mol(r), recs[r['id']]['plane'][0]) for r in members]
                v = ser(smiles(members[0]['core']), [(m, dict(zip(m, map(tuple, p)))) for m, p in ms])
                if v is not None:
                    vs.append(v)
            s['ser'] = fmean(vs) if vs else None
            reviewed = [r for r in mine if r['id'] in refs and recs[r['id']]['ok']]
            miss = units = 0                                # ring systems off their reviewed pose
            for r in reviewed:
                m = _mol(r)
                plane = dict(zip(m, map(tuple, recs[r['id']]['plane'][0])))
                d = pose(m, plane, dict(zip(m, map(tuple, refs[r['id']]))))
                miss += sum(v > .15 for v in d)
                units += len(d)
            s['pose'] = f'{miss}/{units}' if reviewed else None
            per[stratum] = s
    return table


KEYS = ('n', 'failed', *COUNTS, *MEANS, 'ser', 'pose', 't50', 't95')


def _cell(v):
    return '' if v is None else v if isinstance(v, str) else f'{v:.3g}'


def report(outdir, corpus, poses=None):
    outdir = Path(outdir)
    table = summarize(outdir, corpus, poses)
    md = ['| engine | stratum | ' + ' | '.join(KEYS) + ' |', '| --- ' * (len(KEYS) + 2) + '|']
    tsv = ['engine\tstratum\t' + '\t'.join(KEYS)]
    for engine, per in table.items():
        if per == 'absent':
            md.append(f'| {engine} | absent |' + ' |' * len(KEYS))
            continue
        for stratum, s in per.items():
            md.append(f'| {engine} | {stratum} | ' + ' | '.join(_cell(s[k]) for k in KEYS) + ' |')
            tsv.append(f'{engine}\t{stratum}\t' + '\t'.join(_cell(s[k]) for k in KEYS))
    (outdir / 'scorecard.md').write_text('\n'.join(md) + '\n')
    (outdir / 'scorecard.tsv').write_text('\n'.join(tsv) + '\n')
    _gallery(outdir, corpus)


def _gallery(outdir, corpus):
    from ..figure import molecule_scene                     # figure imports layout: resolved at call time
    recs = {e: r for e, r in _load(outdir).items() if r != 'absent'}
    rows = read_corpus(corpus)
    for stratum in dict.fromkeys(r['stratum'] for r in rows):
        html = ['<html><body><table border="1"><tr><th>id</th>' + ''.join(f'<th>{e}</th>' for e in recs)
                + '</tr>']
        for row in (r for r in rows if r['stratum'] == stratum):
            obj = smiles(row['smiles'])
            mols = list(obj.molecules()) if isinstance(obj, ReactionContainer) else [obj]
            cells = []
            for e, rs in recs.items():
                rec = rs.get(row['id'])
                if not rec or not rec['ok']:
                    cells.append(f'<td>{rec["error"] if rec else "missing"}</td>')
                    continue
                svgs = [molecule_scene(m, plane=dict(zip(m, map(tuple, p)))).to_svg()
                        for m, p in zip(mols, rec['plane'])]
                cells.append('<td>' + ''.join(svgs) + '</td>')
            html.append(f'<tr><td>{row["id"]}</td>' + ''.join(cells) + '</tr>')
        html.append('</table></body></html>')
        (outdir / f'gallery-{stratum}.html').write_text('\n'.join(html))


def set_pose(outdir, corpus, poses, row_id, engine):
    rec = _load(outdir)[engine][row_id]
    if not rec['ok']:
        raise ValueError(f'{engine} has no plane for {row_id}: {rec["error"]}')
    refs = _read_poses(poses)
    refs[row_id] = rec['plane'][0]
    Path(poses).write_text(''.join(f'{i}\t{json.dumps(p)}\n' for i, p in sorted(refs.items())))


def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser(prog='python -m chython.depict.layout.scorecard')
    sub = ap.add_subparsers(dest='cmd', required=True)
    w = sub.add_parser('worker')
    w.add_argument('engine')
    w.add_argument('corpus')
    w.add_argument('out')
    w.add_argument('after', nargs='?')
    r = sub.add_parser('run')
    r.add_argument('--engines', nargs='+', default=list(ENGINES))
    r.add_argument('--corpus', default=str(CORPUS))
    r.add_argument('--out', default='scorecard')
    for name in ('report', 'pose'):
        x = sub.add_parser(name)
        if name == 'pose':
            x.add_argument('id')
            x.add_argument('engine')
        x.add_argument('--corpus', default=str(CORPUS))
        x.add_argument('--out', default='scorecard')
    a = ap.parse_args(argv)
    if a.cmd == 'worker':
        worker(a.engine, a.corpus, a.out, a.after)
    elif a.cmd == 'run':
        Path(a.out).mkdir(exist_ok=True)
        run(a.engines, a.corpus, a.out)
    elif a.cmd == 'report':
        report(a.out, a.corpus, POSES)
    else:
        set_pose(a.out, a.corpus, POSES, a.id, a.engine)


if __name__ == '__main__':
    main()
