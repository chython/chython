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
"""AlphaZero rounds for the layout controller: a tester picks molecules, the policy plays them by tree search, the judge
scores the drawings, and all three learn; the run resumes where it stopped.

    python -m chython.depict.layout.learn.alphazero --run DIR [--rounds N] [--workers N] ...

| stage | rule |
| --- | --- |
| pool | the corpus' molecules, largest component each, up to `max_atoms` |
| tester | `Tester` over `descriptors` predicts a molecule's outcome; `molecules` picks favour the low, a `UNIFORM` \
share at random, and `liked` liked molecules not held out so the judge has its pairs |
| demos | once per run, `demos.primitive_demonstration` toward each liked pose not held out, in `workers` processes; \
`demos.jsonl` |
| bootstrap | with no policy yet, `round-0000/policy.pt` learns the demonstrations alone |
| self-play | `mcts.episode` per pick in `workers` processes, reward `mcts.outcome`: the judge's score, the metrics, \
the steps, and the distance to the liked drawing where the user liked one |
| judge | liked poses over disliked drawings, `perturbed` liked drawings and this round's drawings of the \
same molecules |
| policy | visit shares and outcomes of the last `window` rounds' records, `keep` per episode and `examples` in all, \
plus up to `examples` / 2 demonstration records, to a plateau of the validation policy loss; value loss \
weight `VALUE` |
| eval | the held-out liked poses drawn greedily by the policy: distance to the liked drawing, judge score, penalty |
| state | `round-NNNN/` holds `policy.pt`, `judge.pt`, `tester.pt`, `records.jsonl` and `metrics.json`; \
`state.json` the last round done |
"""
import json
import random
from argparse import ArgumentParser
from dataclasses import asdict, dataclass, fields
from math import exp
from pathlib import Path
from statistics import median
from time import perf_counter

import torch
from torch import nn

from ....core import smiles as read_smiles
from ..greedy import drive
from ..primitive import Initial
from ..rings import ring_blocks
from .demos import aligned, images, primitive_demonstration
from .judge import Judge, judge_score, pairs_of, perturbed, train_judge
from .mcts import HARD, Search, episode, outcome, penalty, rebuild, step_key
from .policy import Policy
from .tensors import state, steps
from .train import chooser, held_out, split


__all__ = ['Config', 'DESCRIPTORS', 'descriptors', 'Tester', 'rounds', 'main']

DESCRIPTORS = 12
UNIFORM = .25
THREADS = 2
VALUE = 1.


@dataclass
class Config:
    corpus: str = ''
    poses: str = ''
    dislikes: str = ''
    molecules: int = 200
    sims: int = 32
    workers: int = 3
    width: int = 64
    layers: int = 4
    heads: int = 8
    epochs: int = 1000
    patience: int = 30
    window: int = 8
    max_atoms: int = 60
    perturbations: int = 4
    keep: int = 12
    liked: int = 100
    examples: int = 8000
    seed: int = 0


def descriptors(mol):
    """`DESCRIPTORS` numbers of a molecule, invariant to its atom order."""
    n = len(mol)
    bonds = sum(1 for _ in mol.bonds())
    rings = [len(r) for r in mol.sssr]
    degrees = [len(mol.neighbors_of(x)) for x in mol]
    hetero = sum(mol.atom(x).atomic_symbol not in ('C', 'H') for x in mol)
    units = list(mol.stereo_units())
    return [n / 50, bonds / 50, len(rings) / 6, max(rings, default=0) / 12, len(ring_blocks(mol)) / 4, hetero / 10,
            sum(u['kind'] == 1 for u in units) / 4, sum(u['kind'] != 1 for u in units) / 6, max(degrees, default=0) / 6,
            sum(d >= 3 for d in degrees) / 20, (bonds - n + len(mol.connected_components)) / 6,
            sum(d == 1 for d in degrees) / 20]


class Tester(nn.Module):
    """Predicts the outcome a molecule will get; picks the ones it expects to go worst."""
    def __init__(self, width=32):
        super().__init__()
        self.net = nn.Sequential(nn.Linear(DESCRIPTORS, width), nn.GELU(), nn.Linear(width, width), nn.GELU(),
                                 nn.Linear(width, 1))

    def forward(self, x):
        return self.net(x).squeeze(-1)

    def fit(self, data, epochs=200, lr=3e-3):
        """Fit to `(descriptors, outcome)` pairs; returns the last mean squared error."""
        if not data:
            return None
        x = torch.tensor([d for d, _ in data], dtype=torch.float32)
        y = torch.tensor([o for _, o in data], dtype=torch.float32)
        opt = torch.optim.Adam(self.parameters(), lr=lr)
        loss = None
        for _ in range(epochs):
            loss = ((self(x) - y) ** 2).mean()
            opt.zero_grad()
            loss.backward()
            opt.step()
        return float(loss)

    def pick(self, pool, count, seed=0, temperature=1.):
        """`count` sources of `pool`, `[(source, mol)]`, low predicted outcomes favoured, a `UNIFORM` share at
        random."""
        if not pool:
            return []
        with torch.no_grad():
            pred = self(torch.tensor([descriptors(m) for _, m in pool], dtype=torch.float32)).tolist()
        low = min(pred)
        weight = [exp(-(p - low) / temperature) for p in pred]
        total = sum(weight)
        weight = [(1 - UNIFORM) * w / total + UNIFORM / len(pool) for w in weight]
        rng = random.Random(seed)
        return [pool[i][0] for i in rng.choices(range(len(pool)), weights=weight, k=count)]


def _largest(mol):
    return max((set(c) for c in mol.connected_components), key=len)


def _corpus(cfg):
    from ..scorecard import read_corpus
    out = {}
    for r in read_corpus(cfg.corpus):
        if '>' in r['smiles']:
            continue
        mol = read_smiles(r['smiles'])
        if len(_largest(mol)) <= cfg.max_atoms:
            out[r['id']] = r['smiles']
    return out


def _liked(cfg, texts):
    from ..scorecard import _read_poses
    out = {}
    for source, plane in _read_poses(cfg.poses).items():
        if source in texts:
            mol = read_smiles(texts[source])
            if len(plane) == len(mol):
                full = dict(zip(mol, map(tuple, plane)))
                out[source] = {x: full[x] for x in _largest(mol)}
    return out


def _dislikes(cfg, texts):
    from ..review import read_drawings
    out = {}
    for source, kept in read_drawings(cfg.dislikes).items():
        if source in texts:
            mol = read_smiles(texts[source])
            comp = _largest(mol)
            out[source] = [{x: tuple(p) for x, p in zip(mol, r['plane']) if x in comp} for r in kept
                           if len(r['plane']) == len(mol)]
    return out


_MODELS = {}


def _load(path, kind, size):
    key = (str(path), kind)
    stamp = Path(path).stat().st_mtime if Path(path).exists() else None
    if key not in _MODELS or _MODELS[key][0] != stamp:
        model = (Policy if kind == 'policy' else Judge)(**size)
        if stamp is not None:
            model.load_state_dict(torch.load(path))
        _MODELS[key] = (stamp, model.eval())
    return _MODELS[key][1]


def _goals(mol, initial, liked):
    """The liked drawing of `initial`'s component under each automorphism, point lists in scope order."""
    atoms = initial.scope.atoms
    return [[liked[m[x]] for x in atoms] for m in images(mol) if all(m[x] in liked for x in atoms)] or \
        [[liked[x] for x in atoms]]


def _distance(env, goals):
    mine = [env.plane[x] for x in env.scope.atoms]
    return min(aligned(mine, g) for g in goals)


def _reward(judge, initial, goals=None):
    def reward(env):
        return outcome(env, judge_score(judge, initial, env.plane) if judge is not None else 0.,
                       _distance(env, goals) if goals else None)
    return reward


def _play(job):
    source, text, policy, judge, size, sims, seed, liked = job
    torch.set_num_threads(1)
    t = perf_counter()
    try:
        mol = read_smiles(text)
        initial = Initial(mol, _largest(mol))
        jm = _load(judge, 'judge', size) if Path(judge).exists() else None
        goals = _goals(mol, initial, {x: tuple(p) for x, p in zip(mol, liked) if p is not None}) if liked else None
        search = Search(_load(policy, 'policy', size), _reward(jm, initial, goals), sims=sims)
        records, final = episode(source, initial, search, seed=seed)
        return {'source': source, 'records': records, 'final': [list(final.plane[x]) for x in initial.scope.atoms],
                'outcome': records[-1]['outcome'] if records else None, 'penalty': penalty(final),
                'judge': judge_score(jm, initial, final.plane) if jm is not None else 0.,
                'distance': _distance(final, goals) if goals else None, 'steps': final.steps,
                'seconds': perf_counter() - t}
    except Exception as e:                                  # one molecule's failure is a result, not a crash
        return {'source': source, 'error': f'{type(e).__name__}: {e}', 'seconds': perf_counter() - t}


def _policy_examples(records, initials):
    out = []
    for r in records:
        initial = initials.get(r['source'])
        if initial is None:
            continue
        env = rebuild(initial, r)
        legal = env.legal()
        if not legal:
            continue
        share = torch.tensor([r['visits'].get(repr(step_key(env, s)), 0.) for s in legal])
        if float(share.sum()) <= 0:
            continue
        out.append((r['source'], state(env), steps(env, legal), share / share.sum(), r['outcome']))
    return out


def _train_policy(model, data, cfg, device='cpu'):
    sources = sorted({d[0] for d in data})
    train_s, valid_s = split(sources)
    train = [d for d in data if d[0] in set(train_s)] or data
    valid = [d for d in data if d[0] in set(valid_s)] or data
    opt = torch.optim.Adam(model.parameters(), lr=1e-3)
    sched = torch.optim.lr_scheduler.ReduceLROnPlateau(opt, factor=.5, patience=max(1, cfg.patience // 3))
    rng = random.Random(cfg.seed)
    history = {'train': [], 'valid': [], 'best': None, 'epoch': 0}
    best, stalled = None, 0

    def loss_of(d):
        _, s, a, pi, z = d
        logits, value = model({k: v.to(device) for k, v in s.items()}, {k: v.to(device) for k, v in a.items()})
        return -(pi.to(device) * torch.log_softmax(logits, 0)).sum(), (value - z) ** 2
    history['value'] = []
    for epoch in range(cfg.epochs):
        model.train()
        rng.shuffle(train)
        total = 0.
        for d in train:
            ce, sq = loss_of(d)
            loss = ce + VALUE * sq
            opt.zero_grad()
            loss.backward()
            opt.step()
            total += float(ce)
        model.eval()
        with torch.no_grad():
            parts = [loss_of(d) for d in valid]
        v = sum(float(ce) for ce, _ in parts) / len(valid)       # the plateau is the policy's: the value is noisier
        history['train'].append(total / len(train))
        history['valid'].append(v)
        history['value'].append(sum(float(sq) for _, sq in parts) / len(valid))
        sched.step(v)
        if history['best'] is None or v < history['best'] - 1e-3:
            history['best'], history['epoch'], stalled = v, epoch, 0
            best = {k: t.detach().clone() for k, t in model.state_dict().items()}
        else:
            stalled += 1
            if stalled >= cfg.patience:
                break
    if best is not None:
        model.load_state_dict(best)
    return history


def _evaluate(policy, judge, liked, texts):
    choose = chooser(policy)
    choose.primitive = True
    rows = []
    for source in sorted(liked):
        if not held_out(source):
            continue
        mol = read_smiles(texts[source])
        initial = Initial(mol, set(liked[source]))
        env = drive(initial.env(), choose)
        rows.append({'source': source, 'distance': _distance(env, _goals(mol, initial, liked[source])),
                     'judge': judge_score(judge, initial, env.plane) if judge is not None else 0.,
                     'penalty': penalty(env),
                     'defects': env.defects, 'steps': env.steps})
    if not rows:
        return {'count': 0}
    return {'count': len(rows), 'distance': median(r['distance'] for r in rows),
            'within': sum(r['distance'] < .15 for r in rows), 'judge': median(r['judge'] for r in rows),
            'penalty': median(r['penalty'] for r in rows), 'defect_rows': sum(r['defects'] > 0 for r in rows),
            'steps': median(r['steps'] for r in rows), 'rows': rows}


def _full(plane, text):
    """`plane`, `{atom: (x, y)}` of the largest component, as a point list over every atom of `text`, holes as
    None; a job carries this rather than atom ids."""
    return [list(plane[x]) if x in plane else None for x in read_smiles(text)]


def _map(fn, jobs, workers):
    if not workers:
        return [fn(j) for j in jobs]
    from concurrent.futures import ProcessPoolExecutor
    from multiprocessing import get_context
    with ProcessPoolExecutor(workers, mp_context=get_context('spawn')) as ex:
        return list(ex.map(fn, jobs, chunksize=1))


def _demonstrate(job):
    source, text, liked = job
    torch.set_num_threads(1)
    t = perf_counter()
    try:
        mol = read_smiles(text)
        initial = Initial(mol, _largest(mol))
        goals = _goals(mol, initial, {x: tuple(p) for x, p in zip(mol, liked) if p is not None})
        distance, records = primitive_demonstration(initial, goals, source=source)
        env = rebuild(initial, records[-1])
        z = outcome(env, 0., distance)
        for rec in records:
            rec['outcome'] = z
        return {'source': source, 'distance': distance, 'records': records, 'seconds': perf_counter() - t}
    except Exception as e:                                  # one molecule's failure is a result, not a crash
        return {'source': source, 'error': f'{type(e).__name__}: {e}', 'seconds': perf_counter() - t}


def demonstrations(run, cfg, texts, liked):
    """The run's demonstration records, `demos.jsonl`, computed on first call: one route per liked pose not held
    out."""
    path = run / 'demos.jsonl'
    if not path.exists():
        jobs = [(s, texts[s], _full(liked[s], texts[s])) for s in sorted(liked) if not held_out(s)]
        results = _map(_demonstrate, jobs, cfg.workers)
        with open(path.with_suffix('.tmp'), 'w') as f:
            for x in results:
                if 'error' not in x:
                    for rec in x['records']:
                        f.write(json.dumps(rec) + '\n')
        summary = [{k: v for k, v in x.items() if k != 'records'} for x in results]
        (run / 'demos.json').write_text(json.dumps(summary))
        path.with_suffix('.tmp').rename(path)
    return [json.loads(l) for l in path.read_text().splitlines() if l]


def _demo_records(run, cfg, r):
    path = run / 'demos.jsonl'
    if not path.exists():
        return []
    recs = [json.loads(l) for l in path.read_text().splitlines() if l]
    if len(recs) > cfg.examples // 2:
        recs = random.Random(cfg.seed * 43 + r).sample(recs, cfg.examples // 2)
    return recs


def bootstrap(run, cfg, texts, liked):
    """`round-0000/policy.pt`: a policy learnt from the demonstrations alone, once."""
    here = run / 'round-0000'
    if (here / 'policy.pt').exists():
        return None
    here.mkdir(parents=True, exist_ok=True)
    size = {'width': cfg.width, 'heads': cfg.heads, 'layers': cfg.layers}
    records = demonstrations(run, cfg, texts, liked)
    initials = {}
    for rec in records:
        if rec['source'] not in initials:
            mol = read_smiles(texts[rec['source']])
            initials[rec['source']] = Initial(mol, _largest(mol))
    policy = Policy(**size)
    data = _policy_examples(records, initials)
    history = _train_policy(policy, data, cfg) if data else {'best': None}
    torch.save(policy.state_dict(), here / 'policy.pt')
    (here / 'manifest.json').write_text(json.dumps({'size': size, 'kind': 'primitive'}))
    metrics = {'round': 0, 'config': asdict(cfg), 'policy': {'examples': len(data), 'best': history.get('best'),
                                                             'epochs': len(history.get('train', [])),
                                                             'train': history.get('train', []),
                                                             'valid': history.get('valid', []),
                                                             'value': history.get('value', [])},
               'eval': _evaluate(policy, None, liked, texts)}
    (here / 'metrics.json').write_text(json.dumps(metrics))
    return metrics


def _round(run, cfg, r, texts, liked, dislikes):
    t0 = perf_counter()
    size = {'width': cfg.width, 'heads': cfg.heads, 'layers': cfg.layers}
    prev = run / f'round-{r - 1:04d}'
    here = run / f'round-{r:04d}'
    here.mkdir(parents=True, exist_ok=True)
    pool_texts = texts
    pool = [(s, read_smiles(t)) for s, t in pool_texts.items()]
    tester = Tester()
    if (prev / 'tester.pt').exists():
        tester.load_state_dict(torch.load(prev / 'tester.pt'))
    picks = tester.pick(pool, cfg.molecules, seed=cfg.seed * 100003 + r)
    shown = sorted(s for s in liked if not held_out(s) and s not in set(picks))
    shown = random.Random(cfg.seed * 41 + r).sample(shown, min(cfg.liked, len(shown)))
    picks = sorted(set(picks)) + sorted(shown)
    policy_path, judge_path = prev / 'policy.pt', prev / 'judge.pt'
    jobs = [(s, pool_texts[s], policy_path, judge_path, size, cfg.sims, cfg.seed * 7919 + r * 131 + i,
             _full(liked[s], pool_texts[s]) if s in liked else None) for i, s in enumerate(picks)]
    results = _map(_play, jobs, cfg.workers)
    good = [x for x in results if 'error' not in x]
    rng = random.Random(cfg.seed * 31 + r)
    with open(here / 'records.jsonl', 'w') as f:
        for x in good:
            recs = x['records']
            kept = recs if len(recs) <= cfg.keep else sorted(rng.sample(recs[:-1], cfg.keep - 1),
                                                             key=lambda q: q['steps']) + recs[-1:]
            for rec in kept:
                f.write(json.dumps(rec) + '\n')
    selfplay = {'played': len(results), 'errors': len(results) - len(good),
                'outcome': median([x['outcome'] for x in good if x['outcome'] is not None] or [0.]),
                'penalty': median([x['penalty'] for x in good] or [0.]),
                'judge': median([x['judge'] for x in good] or [0.]),
                'steps': median([x['steps'] for x in good] or [0]),
                'distance': median([x['distance'] for x in good if x['distance'] is not None] or [0.]),
                'defect_rows': sum(x['penalty'] >= HARD for x in good), 'seconds': perf_counter() - t0}
    initials = {}
    for x in good:
        mol = read_smiles(pool_texts[x['source']])
        initials[x['source']] = Initial(mol, _largest(mol))
    # the judge: liked over disliked, perturbed and this round's own drawings of the same molecules
    judge = Judge(**size)
    if judge_path.exists():
        judge.load_state_dict(torch.load(judge_path))
    pairs, by = [], {}
    for x in good:
        by.setdefault(x['source'], []).append(x['final'])
    for source, plane in liked.items():
        if held_out(source):
            continue
        mol = read_smiles(texts[source])
        initial = Initial(mol, set(plane))
        own = [dict(zip(initial.scope.atoms, map(tuple, f))) for f in by.get(source, [])]
        negatives = dislikes.get(source, []) + own + \
            [perturbed(initial, plane, steps=3, seed=cfg.seed + r * 17 + k) for k in range(cfg.perturbations)]
        negatives = [n for n in negatives if set(n) == set(plane) and n != plane]
        pairs += [(source, p) for p in pairs_of(initial, plane, negatives)]
    train_s, valid_s = split(sorted({s for s, _ in pairs}))
    jt = [p for s, p in pairs if s in set(train_s)] or [p for _, p in pairs]
    jv = [p for s, p in pairs if s in set(valid_s)] or jt
    judge_history = train_judge(judge, jt, cfg.epochs, valid=jv, patience=cfg.patience, seed=cfg.seed + r) \
        if jt else {'best': None, 'accuracy': []}
    torch.save(judge.state_dict(), here / 'judge.pt')
    # the policy: visit shares and outcomes of the last `window` rounds
    records = []
    for k in range(max(1, r - cfg.window + 1), r + 1):
        path = run / f'round-{k:04d}' / 'records.jsonl'
        if path.exists():
            records += [json.loads(l) for l in path.read_text().splitlines() if l]
    for rec in records:
        if rec['source'] not in initials and rec['source'] in pool_texts:
            mol = read_smiles(pool_texts[rec['source']])
            initials[rec['source']] = Initial(mol, _largest(mol))
    policy = Policy(**size)
    if policy_path.exists():
        policy.load_state_dict(torch.load(policy_path))
    if len(records) > cfg.examples:
        records = random.Random(cfg.seed * 37 + r).sample(records, cfg.examples)
    shown = _demo_records(run, cfg, r)
    for rec in shown:
        if rec['source'] not in initials:
            mol = read_smiles(texts[rec['source']])
            initials[rec['source']] = Initial(mol, _largest(mol))
    data = _policy_examples(records + shown, initials)
    policy_history = _train_policy(policy, data, cfg) if data else {'best': None}
    torch.save(policy.state_dict(), here / 'policy.pt')
    (here / 'manifest.json').write_text(json.dumps({'size': size, 'kind': 'primitive'}))
    # the tester: descriptors to outcomes of this round
    fit = tester.fit([(descriptors(read_smiles(pool_texts[x['source']])), x['outcome']) for x in good
                      if x['outcome'] is not None])
    torch.save(tester.state_dict(), here / 'tester.pt')
    evaluation = _evaluate(policy, judge, liked, texts)
    metrics = {'round': r, 'config': asdict(cfg), 'selfplay': selfplay,
               'policy': {'examples': len(data), 'best': policy_history.get('best'),
                          'epochs': len(policy_history.get('train', [])), 'train': policy_history.get('train', []),
                          'valid': policy_history.get('valid', []), 'value': policy_history.get('value', [])},
               'judge': {'pairs': len(pairs), 'best': judge_history.get('best'),
                         'accuracy': (judge_history.get('accuracy') or [None])[-1],
                         'train': judge_history.get('train', []), 'valid': judge_history.get('valid', [])},
               'tester': {'mse': fit}, 'eval': evaluation, 'seconds': perf_counter() - t0}
    (here / 'metrics.json').write_text(json.dumps(metrics))
    (run / 'state.json').write_text(json.dumps({'round': r}))
    return metrics


def rounds(run, cfg, count):
    """Play `count` rounds after the last one `run` holds."""
    run = Path(run).expanduser()
    run.mkdir(parents=True, exist_ok=True)
    (run / 'config.json').write_text(json.dumps(asdict(cfg), indent=1))
    done = json.loads((run / 'state.json').read_text())['round'] if (run / 'state.json').exists() else 0
    texts = _corpus(cfg)
    liked = _liked(cfg, texts)
    out = []
    if not done:
        first = bootstrap(run, cfg, texts, liked)
        if first is not None:
            out.append(first)
    for r in range(done + 1, done + count + 1):
        dislikes = _dislikes(cfg, texts)                    # read every round: the user's votes keep arriving
        liked = _liked(cfg, texts)
        out.append(_round(run, cfg, r, texts, liked, dislikes))
    return out


def main(argv=None):
    from ..scorecard import CORPUS, POSES
    ap = ArgumentParser(prog='python -m chython.depict.layout.learn.alphazero')
    ap.add_argument('--run', required=True)
    ap.add_argument('--rounds', type=int, default=1000000)
    for f in fields(Config):
        if f.name not in ('corpus', 'poses', 'dislikes'):
            ap.add_argument(f'--{f.name.replace("_", "-")}', type=type(f.default), default=f.default)
    a = ap.parse_args(argv)
    torch.set_num_threads(THREADS)                         # a shared machine: never every core
    rest = {f.name: getattr(a, f.name) for f in fields(Config) if f.name not in ('corpus', 'poses', 'dislikes')}
    cfg = Config(corpus=str(CORPUS), poses=str(POSES), dislikes=str(POSES.with_name('layout_dislikes.tsv')), **rest)
    for m in rounds(a.run, cfg, a.rounds):
        e, play = m['eval'], m.get('selfplay', {})
        print(f"round {m['round']}: played {play.get('played')}, outcome {play.get('outcome')}, "
              f"judge acc {m.get('judge', {}).get('accuracy')}, policy best {m['policy']['best']}, "
              f"eval distance {e.get('distance')}, within {e.get('within')}, defect rows {e.get('defect_rows')}",
              flush=True)


if __name__ == '__main__':
    main()
