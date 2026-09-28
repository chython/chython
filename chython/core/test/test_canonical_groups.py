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
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with this program; if not, see <https://www.gnu.org/licenses/>.
#
"""Enhanced stereo in identity and in the canonical string: kind and membership count, ids do not.

| surface                                  | reads groups as                                  |
| ---------------------------------------- | ------------------------------------------------- |
| `==`, `hash`, `canonical_bytes`          | ABS/OR/AND kind and partition, configured anchors |
| `str()` (canonical SMILES)               | the same partition seeds the atom order           |
| `format(mol, '!e')`                      | nothing: equal to the ungrouped string            |

An AND group is a mixture of the drawn configuration and its full inverse, so inverting every member
of one AND group is the same statement: identity and the canonical string read each AND group in a
canonical phase.  An OR group is one isolated isomer whose absolute configuration is unknown, and its
signs tell such isomers apart; OR and ABS have no phase.  An explicit ABS label states each member
absolute on its own, so it counts against an unlabelled sign; a label on an unconfigured or
non-stereogenic unit states nothing.
"""
from re import findall

from pytest import mark

from chython.core._core import read_reaction_smiles, read_smiles


def _orders(m, count=20):
    return {str(read_smiles(format(m, 'r'))) for _ in range(count)}


@mark.parametrize('text', [
    # (1R,2R)-cyclohexane-1,2-diamine, CAS 20439-47-8: one of its two equivalent centres in an OR group
    'N[C@@H]1CCCC[C@H]1N |o1:1|',
    # D-mannitol, CAS 69-65-8: C-2 and C-3 in one OR group, which the twofold axis maps onto C-5 and C-4
    'OC[C@@H](O)[C@@H](O)[C@H](O)[C@H](O)CO |o1:2,4|',
    # three copies of butan-2-ol, the partition spanning components
    'C[C@H](O)CC.C[C@H](O)CC.C[C@H](O)CC |&1:1,6,&2:11|',
    # a salt, the group on the organic component
    'C[C@H](O)CC.[Na+].C[C@@H](O)CC |&1:1,o1:8|',
])
def test_every_atom_order_writes_one_string(text):
    m = read_smiles(text)
    out = _orders(m)
    assert len(out) == 1, out
    back = read_smiles(out.pop())
    assert back == m and back.canonical_bytes == m.canonical_bytes and hash(back) == hash(m)


@mark.parametrize('a,b,equal', [
    ('C[C@H](O)CC |&1:1|', 'C[C@H](O)CC', False),                    # racemate vs single enantiomer
    ('C[C@H](O)CC |&1:1|', 'C[C@H](O)CC |o1:1|', False),             # AND vs OR
    ('C[C@H](O)CC |a:1|', 'C[C@H](O)CC', False),                     # an explicit ABS label counts
    ('C[C@H](O)CC.Cl |&1:1|', 'C[C@H](O)CC.Cl', False),              # through the component split
    ('C[C@H](O)[C@H](C)CC |o1:1,o2:3|', 'C[C@H](O)[C@H](C)CC |o2:1,o1:3|', True),   # ids are labels
    ('C[C@H](O)[C@H](C)CC |o1:1,o2:3|', 'C[C@H](O)[C@H](C)CC |o1:1,3|', False),     # partition counts
    ('C[C@H](O)C[C@H](N)C |&1:1,&2:4|', 'C[C@H](O)C[C@H](N)C |&2:1,&1:4|', True),
    ('C[C@H](O)C[C@H](N)C |&1:1,o1:4|', 'C[C@H](O)C[C@H](N)C |&3:1,o7:4|', True),
    ('C[C@H](O)C[C@H](N)C |&1:1,o1:4|', 'C[C@H](O)C[C@H](N)C |o1:1,&1:4|', False),  # kind is not an id
    # D-mannitol, CAS 69-65-8: swapped ids across the twofold axis, and a partition that crosses it
    ('OC[C@@H](O)[C@@H](O)[C@H](O)[C@H](O)CO |&1:2,4,&2:6,8|',
     'OC[C@@H](O)[C@@H](O)[C@H](O)[C@H](O)CO |&2:2,4,&1:6,8|', True),
    ('OC[C@@H](O)[C@@H](O)[C@H](O)[C@H](O)CO |o1:2,&1:4,&2:6,o2:8|',
     'OC[C@@H](O)[C@@H](O)[C@H](O)[C@H](O)CO |o2:2,&2:4,&1:6,o1:8|', True),
    ('OC[C@@H](O)[C@@H](O)[C@H](O)[C@H](O)CO |&1:2,4,&2:6,8|',
     'OC[C@@H](O)[C@@H](O)[C@H](O)[C@H](O)CO |&1:2,8,&2:4,6|', False),
    ('C[C@H](O)CC.C[C@H](O)CC.C[C@H](O)CC |&1:1,6,&2:11|',
     'C[C@H](O)CC.C[C@H](O)CC.C[C@H](O)CC |&1:1,&2:6,11|', True),   # one partition up to symmetry
])
def test_identity_reads_kind_and_partition(a, b, equal):
    ma, mb = read_smiles(a), read_smiles(b)
    assert (ma == mb) is equal
    assert (ma.canonical_bytes == mb.canonical_bytes) is equal
    if equal:
        assert hash(ma) == hash(mb) and format(ma, '!e') == format(mb, '!e')


@mark.parametrize('a,b,equal', [
    ('C[C@H](O)CC |a:1|', 'CC[C@H](C)O |a:2|', True),                # one compound in two atom orders
    ('C[C@H](O)CC |a:1|', 'C[C@@H](O)CC |a:1|', False),              # the sign inverted
    ('C[C@H](O)[C@H](C)CC |a:1,&1:3|', 'C[C@H](O)[C@H](C)CC |&1:3|', False),   # ABS beside a group
    ('CC(O)CC |a:1|', 'CC(O)CC', True),                              # unconfigured: states nothing
    ('C[C@H](C)CC |a:1|', 'CC(C)CC', True),                          # non-stereogenic: states nothing
    # two copies of butan-2-ol: which copy carries the label is a label, how many carry it is not
    ('C[C@H](O)CC.C[C@H](O)CC |a:1,6|', 'C[C@H](O)CC.C[C@H](O)CC |a:6,1|', True),
    ('C[C@H](O)CC.C[C@H](O)CC |a:1|', 'C[C@H](O)CC.C[C@H](O)CC |a:6|', True),
    ('C[C@H](O)CC.C[C@H](O)CC |a:1|', 'C[C@H](O)CC.C[C@H](O)CC |a:1,6|', False),
])
def test_an_explicit_abs_label_is_identity(a, b, equal):
    # the writer states every stored label, so an empty one is compared on the extension-free string
    ma, mb = read_smiles(a), read_smiles(b)
    assert (ma == mb) is equal
    assert (ma.canonical_bytes == mb.canonical_bytes) is equal
    if equal:
        assert hash(ma) == hash(mb) and format(ma, '!e') == format(mb, '!e')
    else:
        assert str(ma) != str(mb)


@mark.parametrize('text', ['C[C@H](O)CC |a:1|', 'C[C@H](O)[C@H](C)CC |a:1,&1:3|',
                           'C[C@H](O)CC.C[C@H](O)CC |a:1|'])
def test_an_abs_label_writes_one_string_in_every_atom_order(text):
    m = read_smiles(text)
    out = _orders(m, 30)
    assert out == {str(m)}, out
    assert ' |a:' in str(m) and (1, 0) in m.canonical_stereo_groups()
    back = read_smiles(str(m))
    assert back == m and hash(back) == hash(m)


@mark.parametrize('text', [
    'N[C@@H]1CCCC[C@H]1N |o1:1|',
    'OC[C@@H](O)[C@@H](O)[C@H](O)[C@H](O)CO |o1:2,4|',
])
def test_the_extension_free_string_ignores_groups(text):
    grouped, plain = read_smiles(text), read_smiles(text.split(' ')[0])
    assert format(grouped, '!e') == format(plain, '!e')
    assert {format(read_smiles(format(grouped, 'r')), '!e') for _ in range(20)} == {format(plain, '!e')}


@mark.parametrize('a,b,equal', [
    ('C[C@H](O)CC |&1:1|', 'C[C@@H](O)CC |&1:1|', True),             # one member: no parity stated
    ('C[C@H](O)CC |o1:1|', 'C[C@@H](O)CC |o1:1|', False),            # OR: the two isolated isomers
    ('C[C@H](O)CC |&1:1|', 'C[C@@H](O)CC |o1:1|', False),            # a mixture is neither isomer
    ('C[C@H](O)CC |a:1|', 'C[C@@H](O)CC |a:1|', False),              # ABS is absolute
    ('C[C@H](O)[C@H](C)CC |&1:1,3|', 'C[C@@H](O)[C@@H](C)CC |&1:1,3|', True),   # the whole group inverted
    ('C[C@H](O)[C@H](C)CC |&1:1,3|', 'C[C@@H](O)[C@H](C)CC |&1:1,3|', False),   # one member inverted
    ('C[C@H](O)[C@H](C)CC |o1:1,3|', 'C[C@@H](O)[C@@H](C)CC |o1:1,3|', False),  # OR has no phase
    ('C[C@H](O)[C@H](C)CC |o1:1,3|', 'C[C@@H](O)[C@H](C)CC |o1:1,3|', False),
    ('C[C@H](O)[C@H](C)CC |a:1,3|', 'C[C@@H](O)[C@@H](C)CC |a:1,3|', False),   # an ABS pair is not a group
    ('C[C@H](O)[C@H](C)CC', 'C[C@@H](O)[C@@H](C)CC', False),
    ('C[C@H](O)[C@H](C)CC |&1:1,3|', 'C[C@H](O)[C@H](C)CC |o1:1,3|', False),    # AND pair vs OR pair
    ('C[C@H](O)[C@H](C)CC |&1:1,3|', 'C[C@@H](O)[C@@H](C)CC |o1:1,3|', False),
    # an AND pair beside an OR pair: only the AND pair's phase moves
    ('C[C@H](O)[C@H](C)C[C@H](O)[C@H](C)CC |&1:1,3,o1:6,8|',
     'C[C@@H](O)[C@@H](C)C[C@H](O)[C@H](C)CC |&1:1,3,o1:6,8|', True),
    ('C[C@H](O)[C@H](C)C[C@H](O)[C@H](C)CC |&1:1,3,o1:6,8|',
     'C[C@H](O)[C@H](C)C[C@@H](O)[C@@H](C)CC |&1:1,3,o1:6,8|', False),
    # two groups, each in its own phase: inverting one of them restates the molecule
    ('C[C@H](O)[C@H](C)C[C@H](O)[C@H](C)CC |&1:1,3,&2:6,8|',
     'C[C@H](O)[C@H](C)C[C@@H](O)[C@@H](C)CC |&1:1,3,&2:6,8|', True),
    ('C[C@H](O)[C@H](C)C[C@H](O)[C@H](C)CC |&1:1,3,&2:6,8|',
     'C[C@H](O)[C@H](C)C[C@@H](O)[C@H](C)CC |&1:1,3,&2:6,8|', False),
    # ABS beside a group: the group's phase moves, the ABS centre does not
    ('C[C@H](O)[C@H](N)[C@H](F)C |a:1,&1:3,5|', 'C[C@H](O)[C@@H](N)[C@@H](F)C |a:1,&1:3,5|', True),
    ('C[C@H](O)[C@H](N)[C@H](F)C |a:1,&1:3,5|', 'C[C@@H](O)[C@H](N)[C@H](F)C |a:1,&1:3,5|', False),
    # meso pentane-2,4-diol against the racemic pair, and each against its own inversion
    ('C[C@H](O)C[C@@H](O)C |&1:1,4|', 'C[C@@H](O)C[C@H](O)C |&1:1,4|', True),
    ('C[C@H](O)C[C@H](O)C |&1:1,4|', 'C[C@@H](O)C[C@@H](O)C |&1:1,4|', True),
    ('C[C@H](O)C[C@@H](O)C |&1:1,4|', 'C[C@H](O)C[C@H](O)C |&1:1,4|', False),
])
def test_a_group_is_read_in_its_canonical_phase(a, b, equal):
    ma, mb = read_smiles(a), read_smiles(b)
    assert (ma == mb) is equal
    assert (ma.canonical_bytes == mb.canonical_bytes) is equal
    assert (str(ma) == str(mb)) is equal
    if equal:
        assert hash(ma) == hash(mb)


@mark.parametrize('text', [
    'C[C@@H](O)CC |&1:1|',
    'C[C@@H](O)[C@H](C)CC |o1:1,3|',
    'C[C@@H](O)CC |o1:1|',
    'C[C@H](O)[C@H](C)C[C@H](O)[C@@H](C)CC |&1:1,3,o1:6,8|',
    'C[C@H](O)[C@H](C)C[C@@H](O)[C@@H](C)CC |&1:1,3,&2:6,8|',
    'C[C@H](O)[C@@H](N)[C@@H](F)C |a:1,&1:3,5|',
    'C[C@@H](O)C[C@H](O)C |&1:1,4|',
    # pentane-2,3,4-triol: C3 is pseudo-asymmetric, and stereogenic only in some phases of the groups
    'C[C@@H](O)[C@H](O)[C@@H](O)C |&1:1,&2:5,o1:3|',
])
def test_the_phase_is_one_string_in_every_atom_order(text):
    m = read_smiles(text)
    out = {str(read_smiles(format(m, 'r'))) for _ in range(20)}
    assert out == {str(m)}, out
    back = read_smiles(str(m))
    assert back == m and str(back) == str(m)


@mark.parametrize('tail,equal', [
    ('&1:1,&2:5', True),        # C2 and C4 free: like or unlike, and C3 either way
    ('&1:1,5', True),           # unlike, stated relatively: inverting the group turns 3r into 3s
    ('&1:1,3,5', False),        # all three relative to each other: the two meso forms
    ('a:1,5', False),           # C2 and C4 absolute, so C3 is fixed against them
    ('o1:1,5', False),          # one isolated isomer: OR has no phase, so C3 is fixed as under ABS
])
def test_a_pseudo_asymmetric_centre_is_read_against_its_neighbours_phases(tail, equal):
    # C3 of pentane-2,3,4-triol is stereogenic only when C2 and C4 differ, and a group's phase can
    # decide whether they do: its parity is kept wherever some phase makes it stereogenic
    a = read_smiles('C[C@@H](O)[C@H](O)[C@@H](O)C |%s|' % tail)
    b = read_smiles('C[C@@H](O)[C@@H](O)[C@@H](O)C |%s|' % tail)
    assert len(findall(r'\[C@@?H\]', str(a).split(' ')[0])) == 3, 'C3 carries a sign'
    assert (a == b) is equal
    assert (str(a) == str(b)) is equal


def test_the_stored_signs_are_untouched():
    # the phase is a reading: the arena keeps what was read, and `!e` writes it
    m = read_smiles('C[C@H](O)CC |&1:1|')
    before = m.to_bytes()
    assert str(m) == 'C(C)[C@H](O)C |&1:2|'
    assert m.to_bytes() == before
    assert format(m, '!e') == 'C(C)[C@@H](O)C'
    # an OR singleton writes the sign it was read with
    assert str(read_smiles('C[C@H](O)CC |o1:1|')) == 'C(C)[C@@H](O)C |o1:2|'
    assert str(read_smiles('C[C@@H](O)CC |o1:1|')) == 'C(C)[C@H](O)C |o1:2|'


def test_a_group_joined_across_the_arrow_keeps_retention_apart_from_inversion():
    # one map number in an AND group on both sides joins the two into one collection, which states the
    # product's configuration relative to the reactant's: those molecules keep their stored signs
    left = '[CH3:1][C@H:2]([OH:3])[CH2:4][CH3:5]>>'
    retention = read_reaction_smiles(left + '[CH3:1][C@H:2]([NH2:3])[CH2:4][CH3:5] |&1:1,6|')
    inversion = read_reaction_smiles(left + '[CH3:1][C@@H:2]([NH2:3])[CH2:4][CH3:5] |&1:1,6|')
    assert format(retention, 'm') != format(inversion, 'm')
    assert format(retention, 'm').endswith('|&1:2,7|')
