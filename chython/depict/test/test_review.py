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
import json
import sys
import threading
from urllib.request import Request, urlopen
from chython.depict.layout.review import Review, read_review, render, serve, write_review
from chython.depict.layout.scorecard import _read_poses
from pytest import mark

# the review page's HTTP server: a dev tool whose sockets close early on Windows (WinError 10038)
served = mark.skipif(sys.platform == 'win32', reason='the review server is not supported on Windows')


def test_a_comment_survives_tabs_newlines_and_backslashes(tmp_path):
    rows = {'a': {'id': 'a', 'verdict': 'like', 'comment': 'x\ty\nz \\t', 'commit': 'abc', 'hash': 'h',
                  'time': 't', 'engine': 'greedy'}}
    write_review(tmp_path / 'r.tsv', rows)
    assert read_review(tmp_path / 'r.tsv') == rows


def test_render_draws_both_sides_of_a_reaction_and_hashes_the_plane():
    out = render({'id': 'r', 'stratum': 'rxn', 'smiles': 'CCO.CC(=O)O>>CCOC(C)=O'})
    assert out['ok'] and len(out['svg']) == len(out['num']) == 3
    again = render({'id': 'r', 'smiles': 'CCO.CC(=O)O>>CCOC(C)=O', 'plane': out['plane']})
    assert again['hash'] == out['hash']


@served
def test_a_like_is_stored_and_posed_and_a_dislike_unposes(tmp_path):
    corpus = tmp_path / 'corpus.tsv'
    corpus.write_text('id\tstratum\tsmiles\nhard:phenol\thard\tOc1ccccc1\nrxn:1\trxn\tCC>>CC\n')
    app = Review(corpus, tmp_path / 'review.tsv', tmp_path / 'poses.tsv')
    server = serve(app, 0, prefill=False)
    threading.Thread(target=server.serve_forever, daemon=True).start()
    url = f'http://127.0.0.1:{server.server_address[1]}'

    def post(body):
        return json.load(urlopen(Request(f'{url}/vote', data=json.dumps(body).encode(), method='POST')))
    try:
        drawn = json.load(urlopen(f'{url}/draw?id=hard:phenol'))
        assert drawn['ok'] and drawn['metrics']['ovl61'] == 0 and '<svg' in drawn['svg'][0]
        assert post({'id': 'hard:phenol', 'verdict': 'like', 'comment': 'fine'})['hash'] == drawn['hash']
        assert post({'id': 'rxn:1', 'verdict': 'like'})['verdict'] == 'like'
        assert read_review(tmp_path / 'review.tsv')['hard:phenol']['comment'] == 'fine'
        assert list(_read_poses(tmp_path / 'poses.tsv')) == ['hard:phenol']     # a reaction is not posed
        listed = {r['id']: r for r in json.load(urlopen(f'{url}/list'))['rows']}
        assert listed['hard:phenol']['verdict'] == 'like' and not listed['hard:phenol']['changed']
        post({'id': 'hard:phenol', 'verdict': 'dislike'})
        assert _read_poses(tmp_path / 'poses.tsv') == {}
    finally:
        server.shutdown()
        app.render.close()


def test_an_old_review_file_reads_as_the_chython_engine_and_keeps_every_field(tmp_path):
    path = tmp_path / 'review.tsv'
    path.write_text('id\tverdict\tcomment\tcommit\thash\ttime\nhard:a\tlike\tok\tabc\t123\tT1\n')
    review = read_review(path)
    assert review['hard:a']['engine'] == 'chython'
    write_review(path, review)
    again = read_review(path)
    assert again == review and again['hard:a']['comment'] == 'ok' and again['hard:a']['hash'] == '123'


@served
def test_a_dislike_of_another_drawing_keeps_the_liked_pose(tmp_path):
    corpus = tmp_path / 'corpus.tsv'
    corpus.write_text('id\tstratum\tsmiles\nhard:phenol\thard\tOc1ccccc1\n')
    poses = tmp_path / 'poses.tsv'
    poses.write_text('hard:phenol\t' + json.dumps([[0., 0.]] * 7) + '\n')     # a liked drawing, not this one
    app = Review(corpus, tmp_path / 'review.tsv', poses, engine='greedy')
    try:
        v = app.vote('hard:phenol', 'dislike')
        assert v['engine'] == 'greedy' and read_review(tmp_path / 'review.tsv')['hard:phenol']['engine'] == 'greedy'
        assert list(_read_poses(poses)) == ['hard:phenol']                      # the goal stays
        app.vote('hard:phenol', 'like')
        assert _read_poses(poses)['hard:phenol'] != [[0., 0.]] * 7              # a like makes it the goal
        app.vote('hard:phenol', 'dislike')
        assert _read_poses(poses) == {}                                         # disliking the goal itself unposes
    finally:
        app.render.close()


@served
def test_a_baseline_names_its_engine(tmp_path):
    from chython.depict.layout.scorecard import score_row, layout_fn
    corpus = tmp_path / 'corpus.tsv'
    corpus.write_text('id\tstratum\tsmiles\nhard:phenol\thard\tOc1ccccc1\n')
    run = tmp_path / 'run'
    run.mkdir()
    rec = score_row({'id': 'hard:phenol', 'stratum': 'hard', 'smiles': 'Oc1ccccc1', 'core': ''}, layout_fn('greedy'))
    (run / 'greedy.jsonl').write_text(json.dumps(rec) + '\n')
    app = Review(corpus, tmp_path / 'review.tsv', tmp_path / 'poses.tsv', baseline=f'{run}:greedy')
    try:
        assert app.base and app.draw('hard:phenol', base=True)['ok']
    finally:
        app.render.close()


@served
def test_a_vote_on_a_drawn_row_does_not_wait_for_the_renderer(tmp_path):
    corpus = tmp_path / 'corpus.tsv'
    corpus.write_text('id\tstratum\tsmiles\nhard:phenol\thard\tOc1ccccc1\n')
    app = Review(corpus, tmp_path / 'review.tsv', tmp_path / 'poses.tsv')
    try:
        drawn = app.draw('hard:phenol')
        app.render.lock.acquire()                       # the renderer busy, as prefill keeps it
        try:
            out = []
            t = threading.Thread(target=lambda: out.append(app.vote('hard:phenol', 'like')))
            t.start()
            t.join(5)
            assert out and out[0]['hash'] == drawn['hash']
        finally:
            app.render.lock.release()
    finally:
        app.render.close()


@served
def test_a_drawn_row_is_served_while_the_renderer_is_busy(tmp_path):
    corpus = tmp_path / 'corpus.tsv'
    corpus.write_text('id\tstratum\tsmiles\nhard:phenol\thard\tOc1ccccc1\n')
    app = Review(corpus, tmp_path / 'review.tsv', tmp_path / 'poses.tsv')
    try:
        drawn = app.draw('hard:phenol')
        app.render.lock.acquire()                       # the renderer busy, as prefill keeps it
        try:
            out = []
            t = threading.Thread(target=lambda: out.append(app.draw('hard:phenol')))
            t.start()
            t.join(5)
            assert out and out[0]['hash'] == drawn['hash']
        finally:
            app.render.lock.release()
    finally:
        app.render.close()


def test_a_worker_killed_on_a_stall_keeps_every_drawing_done():
    from chython.depict.layout.review import Renderer
    r = Renderer()
    try:
        r.fresh()
        r.cache[('hard:a', False)] = {'id': 'hard:a', 'ok': True}
        r.close()                                       # as a stalled row kills the worker
        r.fresh()
        assert r.proc is not None and r.cache == {('hard:a', False): {'id': 'hard:a', 'ok': True}}
    finally:
        r.close()


@served
def test_a_dislike_keeps_the_disliked_drawing(tmp_path):
    from chython.depict.layout.review import read_drawings
    corpus = tmp_path / 'corpus.tsv'
    corpus.write_text('id\tstratum\tsmiles\nhard:phenol\thard\tOc1ccccc1\n')
    app = Review(corpus, tmp_path / 'review.tsv', tmp_path / 'poses.tsv', engine='greedy')
    try:
        drawn = app.draw('hard:phenol')
        app.vote('hard:phenol', 'dislike')
        kept = read_drawings(tmp_path / 'dislikes.tsv')
        assert kept['hard:phenol'][0]['plane'] == drawn['plane'][0] and kept['hard:phenol'][0]['engine'] == 'greedy'
        app.vote('hard:phenol', 'dislike')                  # the same drawing is not kept twice
        assert len(read_drawings(tmp_path / 'dislikes.tsv')['hard:phenol']) == 1
    finally:
        app.render.close()


@served
def test_a_like_approves_the_drawing_and_a_dislike_withdraws_it(tmp_path):
    from chython.depict.layout.approved import read_library
    corpus = tmp_path / 'corpus.tsv'
    corpus.write_text('id\tstratum\tsmiles\nhard:phenol\thard\tOc1ccccc1\n')
    app = Review(corpus, tmp_path / 'review.tsv', tmp_path / 'poses.tsv')
    try:
        drawn = app.draw('hard:phenol')
        app.vote('hard:phenol', 'like')
        assert read_library(tmp_path / 'approved.tsv') == {'hard:phenol': ('Oc1ccccc1', [[round(v, 6) for v in p]
                                                                                       for p in drawn['plane'][0]])}
        app.vote('hard:phenol', 'dislike')
        assert read_library(tmp_path / 'approved.tsv') == {}
    finally:
        app.render.close()
