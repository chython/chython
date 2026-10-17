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
"""The scorecard records every row whatever an engine does with it: a plane, an error, or a crash."""
import json
import sys
from math import nan

import pytest

from chython.depict.layout.scorecard import (ENGINES, layout_fn, read_corpus, report, run, score_row, set_pose,
                                             summarize, worker)


def _corpus(tmp_path, rows):
    p = tmp_path / 'corpus.tsv'
    p.write_text('id\tstratum\tsmiles\tcore\n' + ''.join('\t'.join(r) + '\n' for r in rows))
    return p


def _line(mol):
    return {x: (.825 * k, 0.) for k, x in enumerate(mol)}


def test_a_plane_is_scored():
    rec = score_row({'id': 'hard:propane', 'stratum': 'hard', 'smiles': 'CCC', 'core': ''}, _line)
    assert rec['ok'] and rec['metrics']['ovl'] == 0 and len(rec['plane'][0]) == 3 and rec['t_us'] >= 0


def test_a_raising_engine_is_recorded():
    def boom(mol):
        raise RuntimeError('no layout')
    rec = score_row({'id': 'x', 'stratum': 'hard', 'smiles': 'CC', 'core': ''}, boom)
    assert not rec['ok'] and 'no layout' in rec['error']


def test_a_nan_or_short_plane_fails_the_row():
    nan_lay = (lambda m: {x: (nan, 0.) for x in m})
    short = (lambda m: {next(iter(m)): (0., 0.)})
    for lay in (nan_lay, short):
        rec = score_row({'id': 'x', 'stratum': 'hard', 'smiles': 'CCO', 'core': ''}, lay)
        assert not rec['ok'] and rec['metrics'] is None


def test_a_reaction_row_scores_alignment():
    row = {'id': 'rxn:1', 'stratum': 'rxn', 'smiles': '[CH3:1][CH2:2][OH:3]>>[CH3:1][CH2:2][O:3]C', 'core': ''}
    rec = score_row(row, _line)
    assert rec['ok'] and rec['metrics']['rxn'] is not None


def test_a_crashing_interpreter_is_resumed(tmp_path):
    corpus = _corpus(tmp_path, [('a', 'hard', 'CC', ''), ('b', 'hard', 'CCC', ''), ('c', 'hard', 'CCCC', '')])
    fake = tmp_path / 'fake_worker.py'
    fake.write_text(
        'import os, sys\n'
        'from chython.depict.layout.scorecard import worker\n'
        'def lay(m):\n'
        '    if len(m) == 3:\n'
        '        os._exit(3)\n'
        '    return {x: (.825 * k, 0.) for k, x in enumerate(m)}\n'
        'engine, corpus, out = sys.argv[1:4]\n'
        'worker(engine, corpus, out, after=sys.argv[4] if len(sys.argv) > 4 else None, lay=lay)\n')
    run(['fake'], corpus, tmp_path, cmd=[sys.executable, str(fake)])
    recs = [json.loads(s) for s in (tmp_path / 'fake.jsonl').read_text().splitlines()]
    assert [r['id'] for r in recs] == ['a', 'b', 'c']
    assert [r['ok'] for r in recs] == [True, False, True] and 'exited' in recs[1]['error']


def test_read_corpus_keeps_columns(tmp_path):
    assert read_corpus(_corpus(tmp_path, [('a', 'hard', 'CC', '')])) == \
        [{'id': 'a', 'stratum': 'hard', 'smiles': 'CC', 'core': ''}]


def _records(tmp_path, engine, rows, lay):
    with open(tmp_path / f'{engine}.jsonl', 'w') as f:
        for r in rows:
            f.write(json.dumps(score_row(r, lay)) + '\n')


SERIES = [('ser:ph:1', 'ser', 'Cc1ccccc1', 'c1ccccc1'), ('ser:ph:2', 'ser', 'Oc1ccccc1', 'c1ccccc1'),
          ('hard:propane', 'hard', 'CCC', '')]


def test_summary_per_engine_and_stratum(tmp_path):
    corpus = _corpus(tmp_path, SERIES)
    rows = read_corpus(corpus)
    _records(tmp_path, 'line', rows, _line)
    (tmp_path / 'gone.jsonl').write_text(json.dumps({'absent': 'no toolkit'}) + '\n')
    s = summarize(tmp_path, corpus)
    assert s['gone'] == 'absent'
    assert s['line']['hard']['n'] == 1 and s['line']['hard']['failed'] == 0
    assert s['line']['ser']['ser'] == pytest.approx(0, abs=1e-9)


def test_report_writes_table_and_gallery(tmp_path):
    corpus = _corpus(tmp_path, SERIES)
    _records(tmp_path, 'line', read_corpus(corpus), _line)
    report(tmp_path, corpus)
    md = (tmp_path / 'scorecard.md').read_text()
    assert '| engine |' in md and 'line' in md
    assert '<svg' in (tmp_path / 'gallery-hard.html').read_text()


def test_a_reviewed_pose_scores_zero_on_its_engine(tmp_path):
    corpus = _corpus(tmp_path, SERIES)
    _records(tmp_path, 'line', read_corpus(corpus), _line)
    poses = tmp_path / 'poses.tsv'
    set_pose(tmp_path, corpus, poses, 'hard:propane', 'line')
    assert summarize(tmp_path, corpus, poses)['line']['hard']['pose'] == '0/1'


def test_a_toolkit_failing_its_warm_up_is_absent(tmp_path):
    corpus = _corpus(tmp_path, [('a', 'hard', 'CC', '')])

    def lay(mol):
        raise ValueError('no JVM')
    worker('stub', corpus, tmp_path / 'stub.jsonl', lay=lay)
    assert [json.loads(s) for s in (tmp_path / 'stub.jsonl').read_text().splitlines()] == [
        {'absent': 'ValueError: no JVM'}]


def test_a_routed_cross_link_is_not_scored_as_a_straight_bond():
    # gramicidin S, PubChem 73357: head-to-tail, the closing amide routed below the backbone
    s = 'N%99C(C(C)C)C(=O)NC(CCCN)C(=O)NC(CC(C)C)C(=O)NC(Cc1ccccc1)C(=O)N%10CCCC%10C(=O)NC(C(C)C)C(=O)' \
        'NC(CCCN)C(=O)NC(CC(C)C)C(=O)NC(Cc1ccccc1)C(=O)N%11CCCC%11C(=O)%99'
    rec = score_row({'id': 'peptide:gramicidin S', 'stratum': 'peptide', 'smiles': s, 'core': ''},
                    layout_fn('chython'))
    assert rec['ok'] and (rec['metrics']['crs'], rec['metrics']['aob'], rec['metrics']['lbl']) == (0, 0, 0)


def test_a_routed_plane_is_stored_at_its_lattice_scale():
    """a long routed link would shrink a mean-bond rescale, and a redrawn plane off the lattice routes nothing"""
    from chython import smiles
    from chython.depict.peptide import routed_links
    s = 'N%99C(C(C)C)C(=O)NC(CCCN)C(=O)NC(CC(C)C)C(=O)NC(Cc1ccccc1)C(=O)N%10CCCC%10C(=O)NC(C(C)C)C(=O)' \
        'NC(CCCN)C(=O)NC(CC(C)C)C(=O)NC(Cc1ccccc1)C(=O)N%11CCCC%11C(=O)%99'
    rec = score_row({'id': 'peptide:gramicidin S', 'stratum': 'peptide', 'smiles': s, 'core': ''},
                    layout_fn('chython'))
    mol = smiles(s)
    assert [k.kind for k in routed_links(mol, dict(zip(mol, map(tuple, rec['plane'][0]))))] == ['head_to_tail']


def test_the_chython_engine_is_the_cascade():
    from chython.core import smiles
    from chython.depict.layout.cascade import cascade_layout
    assert 'controller' not in ENGINES
    mol = smiles('CC(=O)Oc1ccccc1C(=O)O')
    a, b = layout_fn('chython')(mol), cascade_layout(mol)

    def centred(p):
        cx, cy = sum(x for x, _ in p.values()) / len(p), sum(y for _, y in p.values()) / len(p)
        return {k: (round(x - cx, 6), round(y - cy, 6)) for k, (x, y) in p.items()}
    assert centred(a) == centred(b)


def test_the_greedy_engine_drives_the_controller_with_the_model_free_chooser():
    assert 'greedy' in ENGINES
    rec = score_row({'id': 'hard:aspirin', 'stratum': 'hard', 'smiles': 'CC(=O)Oc1ccccc1C(=O)O', 'core': ''},
                    layout_fn('greedy'))
    assert rec['ok'] and rec['metrics']['ct'] == 0 and rec['metrics']['perm'] == 0
