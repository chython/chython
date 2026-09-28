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
"""DEAD STEREO-GROUP MEMBERS: stored, ignored by identity, and never written.

A member is dead when its unit is unconfigured or not stereogenic in any combination of AND phases. A whole
group is dead when inverting its live members restates the molecule in every phase of the other AND groups:
a meso pair under any kind. The reader keeps the byte; `str()` and every format spec write the molecule as if
it were absent.
"""
from pytest import mark

from chython.core import INFO, read_smiles, write_smiles


# (with dead members, the same record with the dead members removed)
PAIRS = [
    ('CC(O)CC |o1:1|', 'CC(O)CC'),                                  # unconfigured centre
    ('CCCC |a:1|', 'CCCC'),                                         # no unit at all
    ('C[C@H](C)C |&1:1|', 'C[C@H](C)C'),                            # configured, not stereogenic
    ('CC=C=CC |&1:2|', 'CC=C=CC'),                                  # unconfigured allene
    ('CC=CC |&1:1|', 'CC=CC'),                                      # unconfigured double bond
    ('C[C@H](O)CC.CC(O)CC |&1:1,&2:6|', 'C[C@H](O)CC.CC(O)CC |&1:1|'),     # one collection emptied
    ('C[C@H](O)CC.CC(O)CC |&1:1,6|', 'C[C@H](O)CC.CC(O)CC |&1:1|'),        # one member of two dropped
    ('C[C@H](O)CC.CCCC |o1:1,a:6|', 'C[C@H](O)CC.CCCC |o1:1|'),
]

SPECS = ['', '!e', '!x', '!s', 'i', 'A', 'm', 'h', 'i!e']


@mark.parametrize('dead, clean', PAIRS)
def test_str_is_a_function_of_identity(dead, clean):
    a, b = read_smiles(dead), read_smiles(clean)
    assert a == b and hash(a) == hash(b)
    assert a.stereo_groups() != b.stereo_groups(), 'the reader stores the dead byte'
    for spec in SPECS:
        assert format(a, spec) == format(b, spec), spec
    assert str(a) == str(b)


@mark.parametrize('text, written', [
    ('CC(O)CC |o1:1|', 'C(C)C(O)C'),
    ('CCCC |a:1|', 'C(C)CC'),
    ('C[C@H](C)C |&1:1|', 'C(C)(C)C'),
    ('C[C@H](O)C[C@H](O)C |a:1,o1:4|', 'O[C@@H](C)C[C@@H](C)O |a:1,o1:4|'),
])
def test_the_written_string(text, written):
    assert str(read_smiles(text)) == written


def test_a_pseudoasymmetric_member_between_two_and_collections_is_live():
    """Pentane-2,3,4-triol: C3 is stereogenic in some combination of the two AND phases, so it stays."""
    m = read_smiles('C[C@H](O)[C@H](O)[C@H](O)C |&1:1,&2:5,o1:3|')
    assert m.dead_stereo_groups() == {}
    assert (2, 1) in m.live_stereo_groups()
    assert 'o1:' in str(m)


def test_each_dead_member_is_one_info_record_where_the_tail_states_collections():
    m = read_smiles('C[C@H](O)CC.CC(O)CC.CCCC |&1:1,6,o1:11|')
    log = []
    write_smiles(m, '', log=log)
    dead = [r for r in log if r.rule == 'smiles:stereo-group-dead']
    assert len(dead) == 2 and all(r.severity == INFO for r in dead), log
    for spec in ('!e', '!x', '!s'):
        log = []
        write_smiles(m, spec, log=log)
        assert not any(r.rule == 'smiles:stereo-group-dead' for r in log), (spec, log)


def test_live_dead_and_stored_memberships():
    m = read_smiles('C[C@H](O)CC.CC(O)CC |&1:1,6|')
    stored = m.stereo_groups()
    assert m.canonical_stereo_groups() == {(3, 1): [2, 7]}
    assert m.live_stereo_groups() == {(3, 1): [2]}
    assert m.dead_stereo_groups() == {(3, 1): [7]}
    assert m.stereo_is_live(2) and not m.stereo_is_live(7)
    allene = read_smiles('CC=C=CC |&1:2|')
    assert not allene.stereo_is_live((2, 4)) and not allene.stereo_is_live(3)
    configured = read_smiles('C[CH]=[C@]=[CH]C')
    assert configured.stereo_is_live((2, 4)) and configured.stereo_is_live(3)
    assert m.stereo_groups() == stored, 'asking changes nothing stored'


def test_copy_and_pach_keep_the_dead_byte():
    m = read_smiles('CC(O)CC |o1:1|')
    assert m.copy().stereo_groups() == m.stereo_groups() == {(2, 1): [2]}
    back = type(m).unpack(m.pack())
    assert back.stereo_groups() == m.stereo_groups()
    assert str(back) == str(m) == 'C(C)C(O)C'


# --- a group whose inversion restates the molecule ------------------------------------------------------------

#: meso-pentane-2,4-diol in two atom orders, the pair of centres at 1 and 4 in both
MESO_DIOL = ['C[C@H](O)C[C@H](O)C', 'O[C@@H](C)C[C@@H](C)O']


@mark.parametrize('text', MESO_DIOL)
@mark.parametrize('tail', ['|&1:1,4|', '|o1:1,4|', '|a:1,4|'])
def test_a_group_over_a_meso_pair_is_dead(text, tail):
    plain = read_smiles(MESO_DIOL[0])
    m = read_smiles(f'{text} {tail}')
    assert m == plain and hash(m) == hash(plain)
    assert str(m) == str(plain) == 'O[C@@H](C)C[C@@H](C)O'
    assert m.live_stereo_groups() == {}
    assert sorted(m.dead_stereo_groups().popitem()[1]) == [2, 5]


@mark.parametrize('tail, other', [
    ('|&1:1,&2:4|', '|&1:1,4|'),       # a split meso pair: neither group restates the molecule alone
    ('|o1:1,o2:4|', '|o1:1,4|'),
    ('|&1:1|', ''),                    # C4 left absolute
    ('|a:1|', ''),
    ('|o1:1|', ''),
])
def test_a_group_over_part_of_a_meso_pair_is_live(tail, other):
    first, second = (read_smiles(f'{t} {tail}') for t in MESO_DIOL)
    assert first == second and str(first) == str(second)
    assert first.dead_stereo_groups() == {}
    assert first != read_smiles(f'{MESO_DIOL[0]} {other}')


def test_a_group_over_a_meso_pair_and_a_chiral_centre_is_live():
    """meso-Pentane-2,4-diol with butan-2-amine: inverting the group gives the other amine enantiomer."""
    first = read_smiles('C[C@H](O)C[C@H](O)C.C[C@H](N)CC |&1:1,4,8|')
    second = read_smiles('C[C@@H](N)CC.O[C@@H](C)C[C@@H](C)O |&1:1,6,9|')
    assert first == second and str(first) == str(second)
    assert first.dead_stereo_groups() == {}
    assert first != read_smiles('C[C@H](O)C[C@H](O)C.C[C@H](N)CC')
    assert first != read_smiles('C[C@H](O)C[C@H](O)C.C[C@H](N)CC |&1:1,4|')


def test_only_the_meso_group_of_two_dies():
    """meso-Pentane-2,4-diol and racemic ephedrine: the ephedrine group keeps canonical id 1 in every order."""
    plain = read_smiles('C[C@H](O)C[C@H](O)C.C[C@H](NC)[C@@H](O)c1ccccc1 |&1:8,11|')
    forms = [
        read_smiles('C[C@H](O)C[C@H](O)C.C[C@H](NC)[C@@H](O)c1ccccc1 |&1:1,4,&2:8,11|'),
        read_smiles('C[C@H](NC)[C@@H](O)c1ccccc1.C[C@H](O)C[C@H](O)C |&1:13,16,&2:1,4|'),
        read_smiles('C[C@H](NC)[C@@H](O)c1ccccc1.C[C@H](O)C[C@H](O)C |&2:13,16,&1:1,4|'),
    ]
    for m in forms:
        assert m == plain and hash(m) == hash(plain) and str(m) == str(plain)
        assert list(m.live_stereo_groups()) == [(3, 1)]
        assert len(m.dead_stereo_groups()) == 1
    assert str(plain) == 'c1cccc(c1)[C@@H](O)[C@H](NC)C.C[C@@H](O)C[C@H](C)O |&1:6,8|'


@mark.parametrize('text', ['C[C@H](O)[C@H](O)[C@H](O)C', 'C[C@H](O)[C@@H](O)[C@H](O)C'])
@mark.parametrize('kind', ['&1', 'o1', 'a'])
def test_a_group_over_every_centre_of_a_meso_triol_is_dead(text, kind):
    """Both meso pentane-2,3,4-triols: the three centres inverted together are the same compound."""
    m = read_smiles(f'{text} |{kind}:1,3,5|')
    assert m == read_smiles(text) and str(m) == str(read_smiles(text))
    assert m.live_stereo_groups() == {}


def test_a_group_is_dead_only_in_every_phase_of_the_other_and_groups():
    """Pentane-2,3,4-triol, C2 and C4 drawn unlike: inverting o1:3 restates the like phase only, so it is live."""
    for text in ('C[C@@H](O)[C@H](O)[C@@H](O)C', 'C[C@H](O)[C@H](O)[C@@H](O)C'):
        m = read_smiles(f'{text} |&1:1,&2:5,o1:3|')
        assert m.dead_stereo_groups() == {}, text
        back = read_smiles(str(m))
        assert back == m and str(back) == str(m), text
    back = read_smiles(str(m))
    assert back == m and str(back) == str(m)
