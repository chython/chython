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
"""A dead stereo-group member -- on a unit unconfigured or not stereogenic -- is written by every CTAB and
MRV writer as stored; after `canonicalize_stereo_groups()` the text is the text of the molecule without it."""
from pytest import mark

from chython.core import INFO, read_smiles

from .._mrv import write_mrv
from ...ctfile import mol


PAIRS = [
    ('CC(O)CC |o1:1|', 'CC(O)CC'),
    ('CCCC |a:1|', 'CCCC'),
    ('C[C@H](C)C |&1:1|', 'C[C@H](C)C'),
    ('CC=C=CC |&1:2|', 'CC=C=CC'),
    ('C[C@H](O)CC.CC(O)CC |&1:1,&2:6|', 'C[C@H](O)CC.CC(O)CC |&1:1|'),
    ('C[C@H](O)CC.CCCC |a:1,o1:6|', 'C[C@H](O)CC.CCCC |a:1|'),
]

WRITERS = [
    ('v2000', lambda m, log: mol(m, version=2000, log=log)),
    ('v3000', lambda m, log: mol(m, version=3000, log=log)),
    ('mrv', lambda m, log: write_mrv(m, log=log)),
]


def _canonical(text):
    m = read_smiles(text)
    m.canonicalize_stereo_groups()
    return m


@mark.parametrize('dead, clean', PAIRS)
@mark.parametrize('name, write', WRITERS)
def test_a_dead_member_is_written_until_the_rewrite_drops_it(dead, clean, name, write):
    a, b = read_smiles(dead), read_smiles(clean)
    assert a != b
    log = []
    write(a, log)
    assert not [r for r in log if r.rule.endswith(':stereo-group-dead')], log
    a, b = _canonical(dead), _canonical(clean)
    assert a == b
    log_a, log_b = [], []
    assert write(a, log_a) == write(b, log_b), name
    assert [r.rule for r in log_a] == [r.rule for r in log_b]


def test_the_default_version_is_moved_by_a_dead_member_until_the_rewrite():
    log = []
    assert 'V3000' in mol(read_smiles('CC(O)CC |o1:1|'), log=log).split('\n')[3]
    log = []
    text = mol(_canonical('CC(O)CC |o1:1|'), log=log)
    assert 'V2000' in text.split('\n')[3]
    assert not any(r.rule == 'v2000:enhanced-stereo-not-written' for r in log), log


def test_the_v2000_chiral_flag_reads_the_stored_groups():
    # an OR group on a configured centre clears the flag; a dead OR member beside an ABS one does too, as stored,
    # until the rewrite drops it
    live = mol(read_smiles('C[C@H](O)CC |o1:1|'), version=2000)
    beside = mol(read_smiles('C[C@H](O)CC.CC(O)CC |o1:6|'), version=2000)
    bare = mol(read_smiles('C[C@H](O)CC.CC(O)CC'), version=2000)
    assert live.split('\n')[3][12:15] == '  0'
    assert beside != bare
    assert mol(_canonical('C[C@H](O)CC.CC(O)CC |o1:6|'), version=2000) == bare
