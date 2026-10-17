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
"""A review page over `test/layout_corpus.tsv`: like, dislike and comment on each drawing.

    python -m chython.depict.layout.review [--port 8765] [--engine E] [--baseline DIR[:E]] [--corpus P]

| Piece | Rule |
| --- | --- |
| drawing | the current tree's `--engine` (`chython` by default, `greedy`, or `policy:DIR` for a trained policy), in a \
worker subprocess restarted when a depict source changes |
| verdict | `test/layout_review.tsv`: id, verdict, comment, commit, plane hash, time, engine; written on every vote |
| like | also writes the plane to `test/layout_poses.tsv`, the goal training aims at; single molecules only; on the \
`chython` engine, approves the drawing: it goes to `test/layout_approved.tsv`, the evaluation set the engine is \
held to by `test_every_approved_drawing_is_kept`; the engine never reads it |
| dislike | removes the pose only when the disliked drawing is the pose: a policy falling short keeps its goal; \
the disliked drawing itself goes to `test/layout_dislikes.tsv`, a negative for the judge; a dislike or a cleared \
verdict withdraws an approval |
| page | opens on the unreviewed rows: a molecule the engine draws as approved, or as a full analog's approved \
drawing, needs no second look; an analog drawn otherwise stays unreviewed |
| changed | the drawing's plane hash differs from the one the verdict was given on |
| baseline | `DIR/E.jsonl` of a scorecard run (`E` the engine after a colon, `chython` by default), drawn from its \
stored plane |

Keys: `j`/`k` move, `l` like, `d` dislike, `c` comment, `n` atom numbers, `b` baseline.
"""
import json
import os
import subprocess
import sys
import threading
import time
from argparse import ArgumentParser
from hashlib import sha1
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from select import select
from urllib.parse import parse_qs, urlparse
from .approved import LIBRARY, approve, held, withdraw
from .scorecard import CORPUS, COUNTS, POSES, ROOT, _load, _read_poses, layout_fn, read_corpus, score_row


REVIEW = CORPUS.with_name('layout_review.tsv')
FIELDS = ('id', 'verdict', 'comment', 'commit', 'hash', 'time', 'engine')
SOURCES = ROOT / 'chython' / 'depict'
STALL_S = 120                                               # a worker silent this long on one row is killed


# --- verdicts -------------------------------------------------------------------------------------------------

def _escape(s):
    return s.replace('\\', '\\\\').replace('\t', '\\t').replace('\n', '\\n').replace('\r', '')


def _unescape(s):
    out, it = [], iter(s)
    for c in it:
        if c == '\\':
            c = next(it, '')
            out.append({'t': '\t', 'n': '\n'}.get(c, c))
        else:
            out.append(c)
    return ''.join(out)


def read_review(path):
    if not Path(path).exists():
        return {}
    lines = Path(path).read_text().splitlines()
    out = {}
    for s in lines[1:]:
        if s:
            r = dict(zip(FIELDS, map(_unescape, s.split('\t'))))
            r['engine'] = r.get('engine') or 'chython'      # a row from before the engine column judged chython
            out[r['id']] = r
    return out


def write_review(path, review):
    body = ''.join('\t'.join(_escape(r.get(k, '')) for k in FIELDS) + '\n' for _, r in sorted(review.items()))
    tmp = Path(f'{path}.tmp')
    tmp.write_text('\t'.join(FIELDS) + '\n' + body)
    tmp.replace(path)


def read_drawings(path):
    """`{id: [{plane, engine, hash, time}]}` of the drawings kept in `path`."""
    if not Path(path).exists():
        return {}
    return {i: json.loads(r) for i, r in (s.split('\t', 1) for s in Path(path).read_text().splitlines() if s)}


def write_drawings(path, drawings):
    tmp = Path(f'{path}.tmp')
    tmp.write_text(''.join(f'{i}\t{json.dumps(r)}\n' for i, r in sorted(drawings.items())))
    tmp.replace(path)


def write_poses(path, poses):
    tmp = Path(f'{path}.tmp')
    tmp.write_text(''.join(f'{i}\t{json.dumps(p)}\n' for i, p in sorted(poses.items())))
    tmp.replace(path)


def plane_hash(plane):
    return sha1(json.dumps([[[round(v, 3) for v in xy] for xy in p] for p in plane]).encode()).hexdigest()[:12]


def _commit():
    try:
        head = subprocess.run(['git', 'rev-parse', '--short', 'HEAD'], cwd=ROOT, capture_output=True,
                              text=True).stdout.strip()
        dirty = subprocess.run(['git', 'status', '--porcelain', '--', str(SOURCES)], cwd=ROOT,
                               capture_output=True, text=True).stdout.strip()
    except OSError:
        return ''
    return f'{head}+' if dirty else head


# --- render worker --------------------------------------------------------------------------------------------

def render(req, lay=None):
    """One corpus row as SVG, plain and atom-numbered: laid out by `lay`, or drawn from `req['plane']`."""
    from ...core import ReactionContainer, smiles
    from ..figure import molecule_scene                     # figure imports layout: resolved at call time
    from ..overlay import ValueLabels

    if req.get('plane') is None:
        rec = score_row(req, lay or layout_fn('chython'))
    else:
        rec = {'id': req['id'], 'ok': True, 'error': None, 't_us': None, 'plane': req['plane'],
               'metrics': req.get('metrics')}
    out = {k: rec[k] for k in ('id', 'ok', 'error', 't_us', 'metrics')}
    if not rec['ok']:
        return out
    obj = smiles(req['smiles'])
    mols = list(obj.molecules()) if isinstance(obj, ReactionContainer) else [obj]
    planes = [dict(zip(m, map(tuple, p))) for m, p in zip(mols, rec['plane'])]
    out['svg'] = [molecule_scene(m, plane=p).to_svg() for m, p in zip(mols, planes)]
    out['num'] = [molecule_scene(m, plane=p, overlays=[ValueLabels({n: n for n in m}, fmt='{:d}',
                                                                   color='#c00000')]).to_svg()
                  for m, p in zip(mols, planes)]
    out['plane'] = rec['plane']
    out['hash'] = plane_hash(rec['plane'])
    out['approved'] = all(held(m, p) for m, p in zip(mols, planes))
    return out


def worker(engine='chython'):
    lay = layout_fn(engine)
    for line in sys.stdin:
        req = json.loads(line)
        try:
            out = render(req, lay)
        except Exception as e:                              # a drawing's failure is a result, not a crash
            out = {'id': req['id'], 'ok': False, 'error': f'{type(e).__name__}: {e}'}
        sys.stdout.write(json.dumps(out) + '\n')
        sys.stdout.flush()


def source_stamp(root=SOURCES):
    """Latest modification time over the depict sources and the compiled core."""
    paths = [p for p in Path(root).rglob('*.py') if 'test' not in p.parts]
    paths += (Path(root).parent / 'core').glob('_core*.so')
    return max((p.stat().st_mtime for p in paths), default=0.)


class Renderer:
    """A worker subprocess serving `render`, replaced when `source_stamp` moves; results cached per stamp."""
    def __init__(self, stall_s=STALL_S, engine='chython'):
        self.stall_s = stall_s
        self.engine = engine
        self.lock = threading.Lock()
        self.proc = None
        self.stamp = None
        self.cache = {}

    def _spawn(self):
        self.close()
        env = {**os.environ, 'PYTHONPATH': os.pathsep.join(filter(None, (str(ROOT), os.environ.get('PYTHONPATH'))))}
        self.proc = subprocess.Popen([sys.executable, '-m', 'chython.depict.layout.review', 'worker', '--engine',
                                      self.engine], cwd=ROOT,
                                     env=env, stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True, bufsize=1)

    def close(self):
        if self.proc is not None:
            self.proc.kill()
            self.proc.wait()
            self.proc = None

    def fresh(self):
        """Restart on a source change, dropping the cache; respawn a worker killed on a stall, keeping it."""
        stamp = source_stamp()
        if stamp != self.stamp:
            self.stamp = stamp
            self.cache.clear()
            self._spawn()
        elif self.proc is None or self.proc.poll() is not None:
            self._spawn()

    def __call__(self, row, plane=None, metrics=None):
        key = (row['id'], plane is not None)
        if key in self.cache and self.stamp == source_stamp():     # a drawing done never waits behind prefill
            return self.cache[key]
        with self.lock:
            self.fresh()
            if key in self.cache:
                return self.cache[key]
            req = {'id': row['id'], 'stratum': row['stratum'], 'smiles': row['smiles'], 'plane': plane,
                   'metrics': metrics}
            self.proc.stdin.write(json.dumps(req) + '\n')
            self.proc.stdin.flush()
            ready, _, _ = select([self.proc.stdout], [], [], self.stall_s)
            line = self.proc.stdout.readline() if ready else ''
            if not line:
                self.close()
                out = {'id': row['id'], 'ok': False, 'error': f'no drawing in {self.stall_s} s'}
            else:
                out = json.loads(line)
            self.cache[key] = out
            return out

    def cached(self, row_id):
        return self.cache.get((row_id, False))


# --- server ---------------------------------------------------------------------------------------------------

class Review:
    def __init__(self, corpus=CORPUS, review=REVIEW, poses=POSES, baseline=None, stall_s=STALL_S, engine='chython',
                 approved=None):
        self.rows = {r['id']: r for r in read_corpus(corpus)}
        self.engine = engine
        self.review_path, self.poses_path = Path(review), Path(poses)
        if approved is None:                                # the library for the real review, a neighbour for any other
            approved = LIBRARY if Path(review).resolve() == REVIEW.resolve() else Path(review).with_name('approved.tsv')
        self.approved_path = Path(approved)
        self.dislikes_path = self.poses_path.with_name(self.poses_path.name.replace('poses', 'dislikes'))
        self.review = read_review(review)
        self.base = {}
        if baseline is not None:
            where, _, name = str(baseline).rpartition(':') if ':' in str(baseline) else (baseline, '', 'chython')
            recs = _load(where).get(name or 'chython')
            self.base = recs if isinstance(recs, dict) else {}
        self.render = Renderer(stall_s, engine)
        self.lock = threading.Lock()

    def listing(self):
        out = []
        for i, row in self.rows.items():
            live, base = self.render.cached(i), self.base.get(i)
            rec = live or base
            v = self.review.get(i, {})
            out.append({'id': i, 'stratum': row['stratum'], 'live': live is not None,
                        'approved': bool(live and live.get('approved')),
                        'metrics': rec and rec.get('metrics'), 'error': rec and rec.get('error'),
                        'verdict': v.get('verdict', ''), 'comment': v.get('comment', ''), 'hash': v.get('hash', ''),
                        'changed': bool(v.get('hash') and live and live.get('hash') and v['hash'] != live['hash'])})
        return {'rows': out, 'baseline': bool(self.base), 'counts': COUNTS}

    def draw(self, row_id, base=False):
        row = self.rows[row_id]
        if base:
            rec = self.base.get(row_id)
            if not rec or not rec['ok']:
                return {'id': row_id, 'ok': False, 'error': rec['error'] if rec else 'not in the baseline'}
            return self.render(row, rec['plane'], rec['metrics'])
        return self.render(row)

    def vote(self, row_id, verdict=None, comment=None):
        if row_id not in self.rows:
            raise KeyError(row_id)
        live = self.render.cached(row_id) or self.draw(row_id)     # the drawing on screen: never queue behind prefill
        with self.lock:
            v = self.review.setdefault(row_id, {'id': row_id, 'verdict': '', 'comment': ''})
            if verdict is not None:
                v['verdict'] = verdict
                v['hash'] = live.get('hash', '')
                v['commit'] = _commit()
                v['engine'] = self.engine
                self._pose(row_id, verdict, live)
                self._approve(row_id, verdict, live)
            if comment is not None:
                v['comment'] = comment
            v['time'] = time.strftime('%Y-%m-%dT%H:%M:%S')
            if not v['verdict'] and not v['comment']:
                del self.review[row_id]
            write_review(self.review_path, self.review)
            return {**v, 'changed': False}

    def _approve(self, row_id, verdict, live):
        """A like on the `chython` engine approves the drawing on screen; any other verdict withdraws it."""
        plane = live.get('plane')
        if verdict == 'like' and self.engine == 'chython' and plane and len(plane) == 1:
            approve(row_id, self.rows[row_id]['smiles'], plane[0], self.approved_path)
        else:
            withdraw(row_id, self.approved_path)

    def _pose(self, row_id, verdict, live):
        plane = live.get('plane')
        if not plane or len(plane) != 1:                    # a reaction's sides are not posed
            return
        poses = _read_poses(self.poses_path)
        if verdict == 'like':
            poses[row_id] = plane[0]
        elif row_id in poses and plane_hash([poses[row_id]]) == plane_hash(plane):
            del poses[row_id]
        else:
            poses = None
        if poses is not None:
            write_poses(self.poses_path, poses)
        if verdict == 'dislike':                            # a disliked drawing is a negative the judge learns from
            kept = read_drawings(self.dislikes_path)
            mine = kept.setdefault(row_id, [])
            h = plane_hash(plane)
            if all(r['hash'] != h for r in mine):
                mine.append({'plane': plane[0], 'engine': self.engine, 'hash': h,
                             'time': time.strftime('%Y-%m-%dT%H:%M:%S')})
                write_drawings(self.dislikes_path, kept)

    def prefill(self, stop):
        """Draw every row in corpus order while the page is idle, so the defect filters see live counts."""
        while not stop.is_set():
            stamp = self.render.stamp
            for i, row in self.rows.items():
                if stop.is_set() or self.render.stamp != stamp:
                    break
                if self.render.cached(i) is None:
                    self.render(row)
            else:
                stop.wait(2.)
                with self.render.lock:
                    self.render.fresh()


def handler(app):
    class Handler(BaseHTTPRequestHandler):
        def log_message(self, *args):
            pass

        def _send(self, body, kind='application/json', code=200):
            data = body.encode() if isinstance(body, str) else json.dumps(body).encode()
            self.send_response(code)
            self.send_header('Content-Type', kind)
            self.send_header('Content-Length', str(len(data)))
            self.send_header('Cache-Control', 'no-store')
            self.end_headers()
            self.wfile.write(data)

        def do_GET(self):
            url = urlparse(self.path)
            q = {k: v[0] for k, v in parse_qs(url.query).items()}
            if url.path == '/':
                self._send(PAGE, 'text/html; charset=utf-8')
            elif url.path == '/list':
                self._send(app.listing())
            elif url.path == '/draw' and q.get('id') in app.rows:
                self._send(app.draw(q['id'], q.get('base') == '1'))
            else:
                self._send({'error': 'not found'}, code=404)

        def do_POST(self):
            if urlparse(self.path).path != '/vote':
                return self._send({'error': 'not found'}, code=404)
            req = json.loads(self.rfile.read(int(self.headers.get('Content-Length', 0))) or b'{}')
            if req.get('verdict') not in (None, '', 'like', 'dislike'):
                return self._send({'error': 'verdict is like, dislike or empty'}, code=400)
            try:
                self._send(app.vote(req.get('id'), req.get('verdict'), req.get('comment')))
            except KeyError:
                self._send({'error': 'no such row'}, code=404)
    return Handler


def serve(app, port=8765, prefill=True):
    server = ThreadingHTTPServer(('127.0.0.1', port), handler(app))
    server.daemon_threads = True
    stop = threading.Event()
    if prefill:
        threading.Thread(target=app.prefill, args=(stop,), daemon=True).start()
    server.stop = stop
    return server


PAGE = r'''<!doctype html>
<html><head><meta charset="utf-8"><title>clean2d review</title><style>
body{font:13px system-ui,sans-serif;margin:0;background:#f4f4f4}
#bar{position:sticky;top:0;z-index:2;background:#fff;border-bottom:1px solid #ccc;padding:6px 10px;
  display:flex;gap:10px;align-items:center;flex-wrap:wrap}
#grid{display:grid;grid-template-columns:repeat(auto-fill,minmax(var(--w,340px),1fr));gap:8px;padding:8px}
.card{background:#fff;border:2px solid #ddd;border-radius:6px;padding:6px;display:flex;flex-direction:column}
.card.focus{border-color:#1a73e8;box-shadow:0 0 0 2px #1a73e855}
.card.like{background:#effaef}.card.dislike{background:#fdeeee}
.head{display:flex;justify-content:space-between;gap:6px;align-items:baseline}
.head a{font-weight:600;color:#222;text-decoration:none}
.bad{color:#b00020;font-weight:600}.chg{color:#e37400;font-weight:600}.muted{color:#888}
.pic{min-height:180px;display:flex;flex-wrap:wrap;align-items:center;justify-content:center;gap:8px}
.pic svg{max-width:100%;height:auto;max-height:var(--h,320px);width:auto}
.pic.big svg{max-height:none}
.row{display:flex;gap:4px;align-items:center}
.row button{border:1px solid #bbb;background:#fff;border-radius:4px;cursor:pointer;font-size:15px}
.row button.on{background:#1a73e8;color:#fff;border-color:#1a73e8}
textarea{flex:1;font:12px system-ui;resize:vertical;min-height:22px;height:22px}
</style></head><body>
<div id="bar">
 <b>clean2d review</b>
 <select id="stratum"><option value="">all strata</option></select>
 <select id="defect"><option value="">any drawing</option><option value="*">any defect</option></select>
 <select id="verdict"><option value="">any verdict</option><option value="none" selected>unreviewed</option>
  <option value="like">liked</option><option value="dislike">disliked</option>
  <option value="comment">commented</option><option value="changed">changed since review</option></select>
 <input id="q" placeholder="id filter" size="14">
 <label><input type="checkbox" id="num"> numbers (n)</label>
 <label><input type="checkbox" id="base"> baseline (b)</label>
 <label>size <input type="range" id="size" min="220" max="900" value="340"></label>
 <span id="stat" class="muted"></span>
 <span class="muted">j/k move · l like · d dislike · c comment · Esc leave</span>
</div>
<div id="grid"></div>
<script>
const $ = s => document.querySelector(s);
let rows = [], shown = [], focus = 0, counts = [];
const drawn = new Map();

async function list() {
  const r = await (await fetch('/list')).json();
  counts = r.counts;
  if (!r.baseline) $('#base').disabled = true;
  if (!rows.length) {
    for (const s of [...new Set(r.rows.map(x => x.stratum))]) $('#stratum').add(new Option(s, s));
    for (const c of counts) $('#defect').add(new Option(c, c));
  }
  rows = r.rows;
  const live = rows.filter(x => x.live).length;
  $('#stat').textContent = `${live}/${rows.length} drawn by this tree`;
  return r;
}
function defects(x) {
  if (!x.metrics) return [];
  return counts.filter(c => x.metrics[c]).map(c => `${c} ${x.metrics[c]}`);
}
function keep(x) {
  const s = $('#stratum').value, d = $('#defect').value, v = $('#verdict').value, q = $('#q').value;
  if (s && x.stratum !== s) return false;
  if (q && !x.id.includes(q)) return false;
  if (d === '*' && !defects(x).length && !x.error) return false;
  if (d && d !== '*' && !(x.metrics && x.metrics[d])) return false;
  if (v === 'none' && (x.verdict || x.approved)) return false;
  if ((v === 'like' || v === 'dislike') && x.verdict !== v) return false;
  if (v === 'comment' && !x.comment) return false;
  if (v === 'changed' && !x.changed) return false;
  return true;
}
function link(id) {
  const m = id.match(/^pubchem:(\d+)$/);
  return m ? `https://pubchem.ncbi.nlm.nih.gov/compound/${m[1]}` : null;
}
function card(x, i) {
  const el = document.createElement('div');
  el.className = 'card ' + x.verdict; el.dataset.id = x.id; el.dataset.i = i;
  const href = link(x.id);
  el.innerHTML = `<div class="head"><span>${href ? `<a href="${href}" target="_blank">${x.id}</a>` : `<b>${x.id}</b>`}
    <span class="muted">${x.stratum}</span></span><span class="info"></span></div>
    <div class="pic muted">…</div>
    <div class="row"><button data-v="like" title="l">👍</button><button data-v="dislike" title="d">👎</button>
    <textarea placeholder="comment (c)">${x.comment.replace(/</g, '&lt;')}</textarea></div>`;
  el.querySelectorAll('button').forEach(b => {
    if (b.dataset.v === x.verdict) b.classList.add('on');
    b.onclick = () => vote(x, x.verdict === b.dataset.v ? '' : b.dataset.v);
  });
  const ta = el.querySelector('textarea');
  ta.onchange = () => comment(x, ta.value);
  ta.onkeydown = e => { if (e.key === 'Enter' && !e.shiftKey) { e.preventDefault(); ta.blur(); }
                        if (e.key === 'Escape') ta.blur(); };
  el.onclick = e => { if (e.target.tagName !== 'TEXTAREA') setFocus(i, false); };
  info(el, x);
  return el;
}
function info(el, x) {
  const d = defects(x);
  el.querySelector('.info').innerHTML = (x.changed ? '<span class="chg">changed </span>' : '') +
    (x.error ? `<span class="bad">${x.error}</span>` : d.length ? `<span class="bad">${d.join(' · ')}</span>` :
     x.metrics ? '<span class="muted">clean</span>' : '');
}
async function paint(el) {
  const x = rows.find(r => r.id === el.dataset.id), base = $('#base').checked, key = x.id + (base ? '#b' : '');
  let r = drawn.get(key);
  if (!r) {
    r = await (await fetch(`/draw?id=${encodeURIComponent(x.id)}${base ? '&base=1' : ''}`)).json();
    drawn.set(key, r);
    if (!base) { x.metrics = r.metrics; x.error = r.error; x.live = true;
                 x.changed = !!(x.hash && r.hash && x.hash !== r.hash); }
  }
  const pic = el.querySelector('.pic');
  pic.classList.remove('muted');
  pic.innerHTML = r.ok ? ($('#num').checked ? r.num : r.svg).join('') : `<span class="bad">${r.error}</span>`;
  if (!base) info(el, x);
  el.dataset.painted = key + ($('#num').checked ? '#n' : '');
}
const seen = new IntersectionObserver(es => es.forEach(e => { if (e.isIntersecting) paint(e.target); }),
                                      {rootMargin: '400px'});
function build() {
  seen.disconnect();
  shown = rows.filter(keep);
  const grid = $('#grid'); grid.innerHTML = '';
  shown.forEach((x, i) => { const el = card(x, i); grid.appendChild(el); seen.observe(el); });
  focus = Math.min(focus, Math.max(shown.length - 1, 0));
  setFocus(focus, false);
}
function repaint() { document.querySelectorAll('.card[data-painted]').forEach(paint); }
function setFocus(i, scroll = true) {
  const cards = document.querySelectorAll('.card');
  if (!cards.length) return;
  focus = Math.max(0, Math.min(i, cards.length - 1));
  cards.forEach(c => c.classList.remove('focus'));
  cards[focus].classList.add('focus');
  if (scroll) cards[focus].scrollIntoView({block: 'nearest', behavior: 'smooth'});
}
async function post(body) {
  return (await fetch('/vote', {method: 'POST', body: JSON.stringify(body)})).json();
}
async function vote(x, v) {
  const r = await post({id: x.id, verdict: v});
  x.verdict = r.verdict || ''; x.hash = r.hash; x.changed = false;
  const el = document.querySelector(`.card[data-id="${CSS.escape(x.id)}"]`);
  el.className = 'card focus ' + x.verdict;
  el.querySelectorAll('button').forEach(b => b.classList.toggle('on', b.dataset.v === x.verdict));
  info(el, x);
}
async function comment(x, text) { const r = await post({id: x.id, comment: text}); x.comment = r.comment || ''; }
document.addEventListener('keydown', e => {
  if (['TEXTAREA', 'INPUT', 'SELECT'].includes(e.target.tagName) && e.target.type !== 'checkbox') return;
  const x = shown[focus], el = document.querySelectorAll('.card')[focus];
  if (e.key === 'j') setFocus(focus + 1);
  else if (e.key === 'k') setFocus(focus - 1);
  else if (e.key === 'l' && x) vote(x, x.verdict === 'like' ? '' : 'like');
  else if (e.key === 'd' && x) vote(x, x.verdict === 'dislike' ? '' : 'dislike');
  else if (e.key === 'c' && el) { e.preventDefault(); el.querySelector('textarea').focus(); }
  else if (e.key === 'n') { $('#num').checked = !$('#num').checked; repaint(); }
  else if (e.key === 'b' && !$('#base').disabled) { $('#base').checked = !$('#base').checked; repaint(); }
});
['#stratum', '#defect', '#verdict'].forEach(s => $(s).onchange = async () => { await list(); build(); });
$('#q').oninput = build;
$('#num').onchange = repaint; $('#base').onchange = repaint;
$('#size').oninput = e => { document.body.style.setProperty('--w', e.target.value + 'px');
                            document.body.style.setProperty('--h', e.target.value + 'px'); };
setInterval(list, 5000);
list().then(build);
</script></body></html>
'''


def main(argv=None):
    ap = ArgumentParser(prog='python -m chython.depict.layout.review')
    ap.add_argument('mode', nargs='?', default='serve', choices=('serve', 'worker'))
    ap.add_argument('--port', type=int, default=8765)
    ap.add_argument('--corpus', default=CORPUS)
    ap.add_argument('--review', default=REVIEW)
    ap.add_argument('--poses', default=POSES)
    ap.add_argument('--baseline', help='a scorecard --out DIR, `:E` naming the engine whose jsonl to draw')
    ap.add_argument('--engine', default='chython', help='chython, greedy, or policy:DIR of a trained policy')
    a = ap.parse_args(argv)
    if a.mode == 'worker':
        return worker(a.engine)
    app = Review(a.corpus, a.review, a.poses, a.baseline, engine=a.engine)
    server = serve(app, a.port)
    print(f'http://localhost:{server.server_address[1]}', flush=True)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.stop.set()
        app.render.close()


if __name__ == '__main__':
    main()
