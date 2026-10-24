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
"""`MoleculeContainer.canonicalize_stereo_groups()`: the explicit rewrite of the stored collections to the
canonical labelling.  `str()` and identity follow the stored collections, so they move with it.

| property                   | pinned by                                                          |
| -------------------------- | ------------------------------------------------------------------ |
| identity follows the store | `==`, `hash`, `canonical_bytes`, `str()` move exactly when it does  |
| the string's labelling     | reading `str()` back gives an arena the method leaves alone         |
| wedges and descriptors     | an inverted unit's wedges and `R`/`S` are inverted with it          |
| convergence                | two drawings differing in ids or AND phase give one `to_bytes()`    |
| idempotence                | the second call answers False                                      |
| pach                       | `pack`/`unpack` and `to_bytes`/`from_bytes` keep the labelling      |
"""
from pytest import mark

from chython.core._core import MoleculeContainer, read_smiles


# (name, drawing, the same compound drawn with other ids or the other AND phase)
PAIRS = [
    ('butan-2-ol', 'C[C@H](O)CC |&1:1|', 'C[C@@H](O)CC |&2:1|'),
    ('pseudoephedrine', 'C[C@H](NC)[C@@H](O)c1ccccc1 |&1:1,4|', 'C[C@@H](NC)[C@H](O)c1ccccc1 |&3:1,4|'),
    ('pentane-2,4-diol', 'C[C@H](O)C[C@@H](O)C |&1:1,&2:4|', 'C[C@H](O)C[C@@H](O)C |&2:1,&1:4|'),
    ('pentane-2,3,4-triol', 'C[C@@H](O)[C@H](O)[C@@H](O)C |&1:1,&2:5,o1:3|',
     'C[C@@H](O)[C@@H](O)[C@@H](O)C |&2:1,&1:5,o2:3|'),
    ('hexane-2,3,4,5-tetrol', 'C[C@@H](O)[C@H](O)[C@H](O)[C@@H](C)O |o1:1,3,&1:5,7|',
     'C[C@@H](O)[C@H](O)[C@@H](O)[C@H](C)O |o2:1,3,&4:5,7|'),
    ('D-mannitol', 'OC[C@@H](O)[C@@H](O)[C@H](O)[C@H](O)CO |&1:2,4,&2:6,8|',
     'OC[C@H](O)[C@H](O)[C@H](O)[C@H](O)CO |&2:2,4,&1:6,8|'),
]


def _identity(m):
    return m == m.copy(), hash(m), m.canonical_bytes, str(m)


def _muconic(first, second):
    """Muconic acid with its two C=C axes in AND collections `first` and `second`."""
    m = read_smiles('OC(=O)/C=C/C=C/C(=O)O')
    with m.edit() as e:
        for pair, group in zip(sorted(m.chiral_bonds()), (first, second)):
            e.set_stereo_group(pair, 3, group)
    return m


@mark.parametrize('name, one, other', PAIRS)
def test_identity_and_the_string_move_exactly_when_the_store_does(name, one, other):
    for text in (one, other):
        m = read_smiles(text)
        before = _identity(m)
        changed = m.canonicalize_stereo_groups()
        assert (_identity(m) == before) is not changed, name


@mark.parametrize('name, one, other', PAIRS)
def test_the_stored_labelling_is_the_one_the_string_states(name, one, other):
    for text in (one, other):
        m = read_smiles(text)
        m.canonicalize_stereo_groups()
        back = read_smiles(str(m))
        # the method reads in slot order, which for a read-back is the order `str()` wrote in
        assert not back.canonicalize_stereo_groups(), (name, str(m))
        assert str(back) == str(m) and back == m
        assert sorted((k, len(v)) for k, v in back.stereo_groups().items()) == \
               sorted((k, len(v)) for k, v in m.stereo_groups().items())


@mark.parametrize('name, one, other', PAIRS)
def test_two_labellings_of_one_compound_converge(name, one, other):
    a = read_smiles(one)
    b = read_smiles(other)
    assert (a == b) == (str(a) == str(b)), name        # identity follows the stored labelling
    a.canonicalize_stereo_groups()
    b.canonicalize_stereo_groups()
    assert a == b and hash(a) == hash(b) and str(a) == str(b), name


@mark.parametrize('name, one, other', PAIRS)
def test_idempotent(name, one, other):
    for text in (one, other):
        m = read_smiles(text)
        m.canonicalize_stereo_groups()
        data = m.to_bytes()
        assert not m.canonicalize_stereo_groups(), name
        assert m.to_bytes() == data


def test_swapped_ids_become_canonical_ids():
    m = read_smiles('C[C@@H](O)CC |&2:1|')
    assert m.stereo_groups() == {(3, 2): [2]} and str(m) == 'C(C)[C@H](O)C |&2:2|'
    assert m.canonicalize_stereo_groups()
    assert m.stereo_groups() == {(3, 1): [2]}
    assert str(m) == 'C(C)[C@H](O)C |&1:2|'


def test_an_and_collection_takes_its_canonical_phase():
    one = read_smiles('C[C@H](O)CC |&1:1|')
    other = read_smiles('C[C@@H](O)CC |&1:1|')
    assert one != other                     # the stored phases differ until the explicit rewrite
    changed = {one.canonicalize_stereo_groups(), other.canonicalize_stereo_groups()}
    assert changed == {True, False}         # exactly one of the two phases is the canonical one
    assert one.parity_of(2) == other.parity_of(2) and one == other


def test_or_and_abs_keep_their_signs():
    m = read_smiles('C[C@H](O)[C@@H](C)[C@H](O)C |o1:1,a:6|')
    signs = {n: m.parity_of(n) for n in (2, 7)}
    m.canonicalize_stereo_groups()
    assert {n: m.parity_of(n) for n in (2, 7)} == signs


def test_a_dead_label_is_dropped():
    m = read_smiles('CC(O)CC |o1:1|')
    assert m.stereo_groups() == {(2, 1): [2]}
    assert m.canonicalize_stereo_groups()
    assert m.stereo_groups() == {}
    assert [r.rule for r in m.log.by_stage('stereo-groups')] == ['stereo-groups:dead']


def test_the_meso_centre_of_the_triol_is_dropped():
    m = read_smiles('C[C@@H](O)[C@H](O)[C@@H](O)C |&1:1,&2:5,o1:3|')
    assert m.canonicalize_stereo_groups()
    rules = [(r.rule, r.atoms) for r in m.log.by_stage('stereo-groups')]
    assert ('stereo-groups:dead', (4,)) in rules
    assert ('stereo-groups:not-stereogenic', (4,)) in rules
    assert not m.parity_of(4)
    assert set(m.stereo_groups()) == {(3, 1), (3, 2)}


def test_the_triol_unlabelled_end_moves_the_and_group_to_one_terminal():
    # an unlabelled sign reads as ABS, OR or AND: with it, one AND terminal states what C2 unlabelled
    # beside a C3-C4 AND pair states
    text = 'C[C@H](O)[C@H](O)[C@@H](O)C |&1:3,5|'
    m = read_smiles(text)
    twin = read_smiles('C[C@H](O)[C@@H](O)[C@H](O)C |&1:3,5|')     # &1 in its other phase
    assert m != twin and str(m) != str(twin)                        # two stored phases, two records
    assert read_smiles(str(m)) == m
    assert m.live_stereo_groups() == {(3, 1): [4, 6]}   # the labels that pair with the stored signs
    assert m.canonicalize_stereo_groups() and twin.canonicalize_stereo_groups()
    assert m == twin and hash(m) == hash(twin) and str(m) == str(twin)
    assert m.stereo_groups() == {(3, 1): [2]}           # one AND terminal states what the pair did
    assert not m.canonicalize_stereo_groups()


@mark.parametrize('name, one, other', PAIRS)
def test_a_canonicalized_read_back_writes_the_same_arena(name, one, other):
    for text in (one, other):
        m = read_smiles(text)
        m.canonicalize_stereo_groups()
        back = read_smiles(str(m))
        again = read_smiles(str(back))
        assert format(again, 'i') == format(back, 'i'), name
        assert again.stereo_groups() == back.stereo_groups()


def test_one_record_per_changed_atom_at_its_stage():
    m = read_smiles('C[C@H](O)C[C@@H](O)C |&1:1,&2:4|')
    before = len(m.log)
    assert m.canonicalize_stereo_groups()
    records = m.log[before:]
    assert records and all(r.stage == 'stereo-groups' and r.severity == 'info' for r in records)
    assert all(len(r.atoms) == 1 for r in records)
    assert all(r.rule in ('stereo-groups:relabelled', 'stereo-groups:phase') for r in records)
    count = len(m.log)
    assert not m.canonicalize_stereo_groups()
    assert len(m.log) == count              # nothing changed, nothing recorded


def test_nothing_to_do_is_a_pure_read():
    m = read_smiles('C[C@H](O)CC')
    copy = m.copy()
    assert not m.canonicalize_stereo_groups()
    assert m.shares_arena_with(copy)
    assert not len(m.log.by_stage('stereo-groups'))


def test_a_shared_arena_is_not_written():
    m = read_smiles('C[C@@H](O)CC |&2:1|')
    copy = m.copy()
    m.canonicalize_stereo_groups()
    assert copy.stereo_groups() == {(3, 2): [2]}


def test_axis_collections_are_renumbered():
    one = _muconic(1, 2)
    other = _muconic(2, 1)
    assert one == other                     # the two axes are exchanged by the molecule's symmetry
    one.canonicalize_stereo_groups()
    other.canonicalize_stereo_groups()
    assert one == other and str(one) == str(other)
    assert not one.canonicalize_stereo_groups()


CASES = [PAIRS[2][2], PAIRS[4][2], PAIRS[3][2]]


def _stereo_state(m):
    return (m.stereo_groups(), m.bond_stereo_groups(), {n: m.parity_of(n) for n in m.atom_numbers},
            str(m), format(m, 'i'))


@mark.parametrize('canonical', [False, True])
@mark.parametrize('text', CASES)
def test_pach_and_arena_bytes_keep_the_labelling(text, canonical):
    m = read_smiles(text)
    if canonical:
        m.canonicalize_stereo_groups()
    state = _stereo_state(m)
    assert _stereo_state(MoleculeContainer.unpack(m.pack())) == state
    assert _stereo_state(MoleculeContainer.from_bytes(m.to_bytes())) == state


@mark.parametrize('canonical', [False, True])
def test_pach_keeps_an_axis_labelling(canonical):
    m = _muconic(2, 1)
    if canonical:
        m.canonicalize_stereo_groups()
    state = _stereo_state(m)
    assert _stereo_state(MoleculeContainer.unpack(m.pack())) == state
    assert _stereo_state(MoleculeContainer.from_bytes(m.to_bytes())) == state


def test_the_drawn_labelling_survives_pach_uncanonicalized():
    m = read_smiles('C[C@H](O)C[C@@H](O)C |&2:1,&1:4|')
    assert m.stereo_groups() == {(3, 2): [2], (3, 1): [5]}
    assert MoleculeContainer.unpack(m.pack()).stereo_groups() == {(3, 2): [2], (3, 1): [5]}
    assert format(MoleculeContainer.unpack(m.pack()), 'i') == format(m, 'i') == 'C[C@H](O)C[C@@H](O)C |&1:4,&2:1|'


def test_random_and_input_order_writes_state_the_stored_ids():
    m = read_smiles('C[C@@H](O)CC |&2:1|')
    assert format(m, 'i') == 'C[C@@H](O)CC |&2:1|'
    assert '|&2:' in format(m, 'r')
    assert str(m) == 'C(C)[C@H](O)C |&2:2|'     # the canonical write states the stored id too
