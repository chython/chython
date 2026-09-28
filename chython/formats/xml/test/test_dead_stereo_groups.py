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
"""A dead stereo-group member -- on a unit unconfigured or not stereogenic -- is written by no CTAB or
MRV writer: the text is the text of the molecule without it, and each dropped member is one INFO line."""
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


@mark.parametrize('dead, clean', PAIRS)
@mark.parametrize('name, write', WRITERS)
def test_a_dead_member_writes_the_text_of_no_member(dead, clean, name, write):
    a, b = read_smiles(dead), read_smiles(clean)
    assert a == b
    log_a, log_b = [], []
    assert write(a, log_a) == write(b, log_b), name
    dropped = [r for r in log_a if r.rule == f'{name}:stereo-group-dead']
    assert dropped and all(r.severity == INFO for r in dropped), log_a
    assert [r for r in log_a if r.rule != f'{name}:stereo-group-dead'] == log_b


def test_the_default_version_is_not_moved_by_a_dead_member():
    log = []
    text = mol(read_smiles('CC(O)CC |o1:1|'), log=log)
    assert 'V2000' in text.split('\n')[3]
    assert not any(r.rule == 'v2000:enhanced-stereo-not-written' for r in log), log


def test_the_v2000_chiral_flag_reads_live_groups_only():
    # an OR group on a configured centre clears the flag; a dead OR member beside an ABS one does not
    live = mol(read_smiles('C[C@H](O)CC |o1:1|'), version=2000)
    beside = mol(read_smiles('C[C@H](O)CC.CC(O)CC |o1:6|'), version=2000)
    bare = mol(read_smiles('C[C@H](O)CC.CC(O)CC'), version=2000)
    assert live.split('\n')[3][12:15] == '  0'
    assert beside == bare
