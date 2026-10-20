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

A member is dead when its unit is unconfigured or not stereogenic in any combination of the context: the
other AND and OR groups and the unconfigured stereogenic units. A whole group is dead when inverting its live
members restates the molecule in every such combination: a meso pair under any kind. A labelling is rewritten
to the least one stating the same mixture, and an ABS or unlabelled parity whose inversion leaves the mixture
unchanged is cleared. The reader keeps the byte; `str()` and every format spec write the molecule as if it
were absent.
"""
import re

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


#: (template, centre slots): pentane-2,3,4-triol and 2,3,4-trihydroxyglutaric acid, both symmetric chains
TRIOL = ('C[C%sH](O)[C%sH](O)[C%sH](O)C', ('1', '3', '5'))
GLUTARIC = ('OC(=O)[C%sH](O)[C%sH](O)[C%sH](O)C(=O)O', ('3', '5', '7'))
FLIP = {'@': '@@', '@@': '@', '': ''}


def both_orders(compound, signs, tail=''):
    """`compound` with `signs` ('' for open) written from either end; `tail` names the forward slots."""
    template, slots = compound
    swap = dict(zip(slots, reversed(slots)))
    backward = re.sub(r'(?<=[:,])\d+', lambda x: swap[x.group()], tail)
    return (f'{template % signs} {tail}'.strip(),
            f'{template % tuple(FLIP[x] for x in reversed(signs))} {backward}'.strip())


ANY_TRIOL = [('@', '@', '@'), ('@', '@@', '@'), ('@', '@', '@@'), ('@@', '@', '@')]


@mark.parametrize('signs', ANY_TRIOL)
def test_an_or_member_between_two_and_collections_states_nothing(signs):
    """Pentane-2,3,4-triol, C2 and C4 in two AND groups: each OR choice at C3 is all four diastereomers."""
    plain = read_smiles(both_orders(TRIOL, (signs[0], '', signs[2]), '|&1:1,&2:5|')[0])
    for text in both_orders(TRIOL, signs, '|&1:1,&2:5,o1:3|'):
        m = read_smiles(text)
        assert m == plain and hash(m) == hash(plain) and str(m) == str(plain), text
        assert list(m.dead_stereo_groups()) == [(2, 1)], text
        back = read_smiles(str(m))
        assert back == m and str(back) == str(m), text


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
def test_a_group_over_a_meso_pair_reads_as_absolute(text, tail):
    """A label over a configuration that is its own mirror image states the drawn diastereomer and no
    enantiomer, as ABS does; unlabelled is unknown, so it is not the unlabelled drawing."""
    plain = read_smiles(MESO_DIOL[0])
    absolute = read_smiles(f'{MESO_DIOL[0]} |a:1,4|')
    m = read_smiles(f'{text} {tail}')
    assert m == absolute and hash(m) == hash(absolute) and m != plain
    assert str(m) == 'O[C@@H](C)C[C@@H](C)O |a:1,4|'
    assert m.live_stereo_groups() == {(1, 0): [2, 5]} and m.dead_stereo_groups() == {}


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


def test_the_meso_group_of_two_reads_as_absolute():
    """meso-Pentane-2,4-diol and racemic pseudoephedrine: the meso group reads as ABS and the pseudoephedrine
    group keeps canonical id 1 in every order."""
    plain = read_smiles('C[C@H](O)C[C@H](O)C.C[C@H](NC)[C@@H](O)c1ccccc1 |a:1,4,&1:8,11|')
    forms = [
        read_smiles('C[C@H](O)C[C@H](O)C.C[C@H](NC)[C@@H](O)c1ccccc1 |&1:1,4,&2:8,11|'),
        read_smiles('C[C@H](NC)[C@@H](O)c1ccccc1.C[C@H](O)C[C@H](O)C |&1:13,16,&2:1,4|'),
        read_smiles('C[C@H](NC)[C@@H](O)c1ccccc1.C[C@H](O)C[C@H](O)C |&2:13,16,&1:1,4|'),
    ]
    for m in forms:
        assert m == plain and hash(m) == hash(plain) and str(m) == str(plain)
        assert sorted(m.live_stereo_groups()) == [(1, 0), (3, 1)]
        assert m.dead_stereo_groups() == {}
    assert str(plain) == 'c1cccc(c1)[C@H](O)[C@@H](NC)C.C[C@@H](O)C[C@H](C)O |a:13,16,&1:6,8|'


@mark.parametrize('text', ['C[C@H](O)[C@H](O)[C@H](O)C', 'C[C@H](O)[C@@H](O)[C@H](O)C'])
@mark.parametrize('kind', ['&1', 'o1', 'a'])
def test_a_group_over_every_centre_of_a_meso_triol_reads_as_absolute(text, kind):
    """Both meso pentane-2,3,4-triols: the three centres inverted together are the same compound, so the group
    states the diastereomer as ABS does, and not the unlabelled drawing."""
    m = read_smiles(f'{text} |{kind}:1,3,5|')
    assert m == read_smiles(f'{text} |a:1,3,5|') and m != read_smiles(text)
    assert list(m.live_stereo_groups()) == [(1, 0)]


#: cyclopropane-1,2,3-tricarboxylic acid and its triamide, C1 open and C2, C3 trans: the cis,cis,trans isomer
TRANS_PAIR = ['OC(=O)C1[C@H](C(=O)O)[C@H]1C(=O)O', 'NC(=O)C1[C@H](C(N)=O)[C@H]1C(N)=O']


@mark.parametrize('text', TRANS_PAIR, ids=['acid', 'amide'])
@mark.parametrize('kind', ['&1', 'o1', 'a'])
def test_a_group_over_a_pseudo_asymmetric_trans_pair_reads_as_absolute(text, kind):
    """Each member is refused in every combination, and the pair pins the cis,cis,trans isomer, which is its
    own mirror image: the group states it as ABS does, and `validate_stereo` keeps both signs."""
    flat = read_smiles(text.replace('@', ''))
    absolute = read_smiles(f'{text} |a:4,8|')
    m = read_smiles(f'{text} |{kind}:4,8|')
    assert m == absolute and hash(m) == hash(absolute) and str(m) == str(absolute)
    assert m != read_smiles(text) and m != flat
    assert re.search(r' \|a:\d+,\d+\|$', str(m)) and m.stereo_is_live(5) and m.stereo_is_live(9)
    assert m.copy().validate_stereo() == []
    c = m.copy()
    c.canonicalize_stereo_groups()
    assert c == absolute and str(c) == str(absolute)


@mark.parametrize('text', [t.replace('@H]1C', '@@H]1C') for t in TRANS_PAIR], ids=['acid', 'amide'])
@mark.parametrize('tail', ['', ' |a:4,8|', ' |o1:4,8|'])
def test_a_cis_pair_beside_an_open_pseudo_asymmetric_centre_states_nothing(text, tail):
    """C2, C3 cis with C1 open: the completions are all-cis and cis,cis,trans, those of the unsigned acid, so
    the pair clears whole, as tropine's bridgeheads do with C3 open."""
    flat = read_smiles(text.replace('@', ''))
    m = read_smiles(text + tail)
    assert m == flat and hash(m) == hash(flat) and str(m) == str(flat) and '@' not in str(m)
    assert m.canonical_bytes == flat.canonical_bytes


@mark.parametrize('text', [t.replace('@H]1C', '@@H]1C') for t in TRANS_PAIR], ids=['acid', 'amide'])
def test_a_cis_pair_in_and_beside_an_open_centre_is_a_mixture(text):
    """In AND the members' mirror leaves C1, so the two phases are all-cis and cis,cis,trans: the mixture,
    as tropine's bridgeheads in AND with C3 open, and not the unsigned acid."""
    m = read_smiles(f'{text} |&1:4,8|')
    assert m != read_smiles(text.replace('@', ''))
    assert str(read_smiles(str(m))) == str(m)


@mark.parametrize('signs', ANY_TRIOL)
def test_an_absolute_centre_between_two_and_collections_states_nothing(signs):
    """`|a:3,&1:1,&2:5|` is `|&1:1,&2:5|` with C3 open, and its string reads back to itself."""
    plain = read_smiles(both_orders(TRIOL, (signs[0], '', signs[2]), '|&1:1,&2:5|')[0])
    for text in both_orders(TRIOL, signs, '|a:3,&1:1,&2:5|'):
        m = read_smiles(text)
        assert m == plain and str(m) == str(plain), text
        back = read_smiles(str(m))
        assert back == m and str(back) == str(m), text


# --- a plain parity whose inversion restates the completions --------------------------------------------------

@mark.parametrize('compound', [TRIOL, GLUTARIC], ids=['triol', 'trihydroxyglutaric'])
@mark.parametrize('signs', [('@', '@', ''), ('@@', '@', '')])
def test_a_centre_beside_an_open_one_is_kept(compound, signs):
    """One end and the middle configured, the other end open: the middle inverted is another mixture."""
    other = (signs[0], FLIP[signs[1]], '')
    first = [read_smiles(t) for t in both_orders(compound, signs)]
    second = [read_smiles(t) for t in both_orders(compound, other)]
    assert first[0] == first[1] and str(first[0]) == str(first[1])
    assert second[0] == second[1] and str(second[0]) == str(second[1])
    assert first[0] != second[0] and str(first[0]) != str(second[0])
    assert len(re.findall(r'@@?', str(first[0]))) == 2
    for m in first + second:
        assert str(read_smiles(str(m))) == str(m)


@mark.parametrize('compound', [TRIOL, GLUTARIC], ids=['triol', 'trihydroxyglutaric'])
@mark.parametrize('sign', ['@', '@@'])
def test_a_middle_centre_between_two_open_ones_states_nothing(compound, sign):
    """Both ends open: the completions of C2 and C4 are every diastereomer, whichever sign C3 has."""
    plain = read_smiles(both_orders(compound, ('', '', ''))[0])
    for text in both_orders(compound, ('', sign, '')):
        m = read_smiles(text)
        assert m == plain and hash(m) == hash(plain) and str(m) == str(plain), text
        assert '@' not in str(m)


@mark.parametrize('tail', ['|&1:1,6,11,16|', '|a:1,6,11,16|', ''])
def test_both_enantiomers_drawn_apart_read_back_under_any_label(tail):
    """trans-2-Methylcyclohexane-1-carboxylic acid and its enantiomer as two components: inverting all four
    centres swaps the components, so an AND over all of them states the pair as ABS does, and either reads
    back to itself; the unlabelled pair is another identity."""
    pair = 'C[C@H]1CCCC[C@@H]1C(=O)O.C[C@@H]1CCCC[C@H]1C(=O)O'
    m = read_smiles(f'{pair} {tail}'.strip())
    back = read_smiles(str(m))
    assert back == m and str(back) == str(m)
    assert (m == read_smiles(pair)) == (not tail)
    if tail:
        assert m == read_smiles(f'{pair} |a:1,6,11,16|')


@mark.parametrize('text', ['CC(O)CC |&1:1|', 'C[C@H](O)CC |w:1.0,&1:1|'])
def test_validate_stereo_heals_an_unsigned_and_singleton(text):
    """Butan-2-ol with no sign, or a wavy one, alone in an AND group is the racemate whatever its sign: the
    centre takes one, the group is kept, and the repair is logged."""
    m = read_smiles(text)
    assert m.validate_stereo() == []
    assert m == read_smiles('C[C@H](O)CC |&1:1|')
    assert [r.severity for r in m.log if r.rule == 'stereo:racemate-unsigned'] == ['repaired']


@mark.parametrize('text', ['C[C@H](O)C(N)C |&1:1,3|', 'CC(O)CC |o1:1|'])
def test_validate_stereo_leaves_an_unsigned_member_that_is_not_an_and_singleton(text):
    """Beside a signed member the unsigned centre's relation is unknown; an OR names one isomer, and with no
    sign not which."""
    m = read_smiles(text)
    m.validate_stereo()
    assert not [r for r in m.log if r.rule == 'stereo:racemate-unsigned']
