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
"""The five S-group invariants of RULES.md section 1.4, each as the test that section demands.

Written against the STORAGE and not against a file format: every fixture is built with `add_atom`, so
a failure here is the arena's and never a parser's.  The MDL epic owns the reader that produces these
records and the writer that consumes them; what is asserted here is only what storage promises, which
is that a record goes in, comes back, and stays correct across an edit that moves every atom index.
"""

import pytest

from chython.core import MoleculeContainer


# Deliberately not valid UTF-8, and deliberately containing a NUL: files in the wild hold bytes like
# these, the blob stores them exactly, and the length is carried separately from the bytes precisely so
# a NUL is not a terminator.  See `blob_rec_t`.
RAW = b'\xff\xfe\x80 caf\xe9 latin-1, a NUL:\x00 and a newline:\n'


def _chain(n=4, element='C'):
    """An n-atom chain, and its stable ids."""
    mol = MoleculeContainer()
    with mol.edit() as e:
        ids = [e.add_atom(element) for _ in range(n)]
        for a, b in zip(ids, ids[1:]):
            e.add_bond(a, b)
    return mol, ids


def _round_trip(mol):
    return MoleculeContainer.from_bytes(mol.to_bytes())


# ------------------------------------------------------------------------------------------------
# Invariant 1: a record with zero atom references SURVIVES.
# ------------------------------------------------------------------------------------------------

def test_a_record_that_loses_every_atom_is_still_there_and_still_carries_its_payload():
    """RULES section 1.4 invariant 1, and the empty-versus-absent distinction of section 6.3.

    An emptied DAT record still asserts that a field was attached to something.  Dropping it would be
    a fidelity loss with no diagnostic anywhere -- the file simply comes back one record short -- which
    is why this is the invariant that bites first.
    """
    mol, ids = _chain(3)
    mol.set_sgroups([{'type': b'DAT', 'name': b'BATCH', 'atoms': (ids[0],),
                      'data': [b'lot-42'], 'index': 1}])
    with mol.edit() as e:
        e.delete_atom(ids[0])

    assert len(mol.sgroups) == 1, 'the record was dropped when its last atom died'
    rec = mol.sgroups[0]
    assert rec['atoms'] == ()
    # THE PAYLOAD IS THE POINT.  A surviving record with its data thrown away would pass a count
    # assertion and still have lost everything the record was for.
    assert rec['type'] == b'DAT'
    assert rec['name'] == b'BATCH'
    assert rec['data'] == (b'lot-42',)
    assert rec['index'] == 1
    assert _round_trip(mol).sgroups == mol.sgroups


def test_an_emptied_record_is_not_confused_with_a_record_that_was_always_empty():
    """The two are the same VALUE and must be, which is what makes the log the only difference.

    Asserted rather than left implicit because it is the reason `sgroup_log` exists at all: storage
    cannot distinguish them, so if the loss is not reported at the moment it happens it is not
    recoverable afterwards from the bytes.
    """
    mol, ids = _chain(3)
    mol.set_sgroups([{'type': b'DAT', 'atoms': (ids[0],)}])
    with mol.edit() as e:
        e.delete_atom(ids[0])

    born_empty, _ = _chain(2)
    born_empty.set_sgroups([{'type': b'DAT', 'atoms': ()}])

    assert mol.sgroups == born_empty.sgroups
    assert mol.sgroup_log and not born_empty.sgroup_log


# ------------------------------------------------------------------------------------------------
# Invariant 2: references are remapped, and a dead one is dropped AND REPORTED -- never zeroed.
# ------------------------------------------------------------------------------------------------

def test_surviving_references_are_remapped_and_never_silently_become_atom_zero():
    """RULES section 1.4 invariant 2.  Id 0 must not become reachable by omission.

    The fixture deletes the FIRST atom on purpose: that is the deletion that shifts every remaining
    index down by one, so a carry that forgot to remap would leave references that are still in range
    and still plausible -- pointing one atom off.  A test that deleted the last atom would pass against
    a carry that did nothing at all.
    """
    mol, ids = _chain(5)
    mol.set_sgroups([{'type': b'SRU', 'atoms': tuple(ids[1:]), 'index': 1}])
    with mol.edit() as e:
        e.delete_atom(ids[0])

    assert mol.sgroups[0]['atoms'] == tuple(ids[1:]), 'references did not follow their atoms'
    assert not mol.sgroup_log, 'no reference was lost, so nothing should have been reported'


def test_a_bond_pair_dies_with_either_endpoint_and_the_loss_is_reported():
    mol, ids = _chain(4)
    mol.set_sgroups([{'type': b'SRU', 'atoms': tuple(ids),
                      'bonds': ((ids[0], ids[1]), (ids[2], ids[3])), 'index': 1}])
    with mol.edit() as e:
        e.delete_atom(ids[0])

    rec = mol.sgroups[0]
    assert rec['bonds'] == ((ids[2], ids[3]),), 'the surviving pair was lost or was not remapped'
    assert rec['atoms'] == tuple(ids[1:])
    assert len(mol.sgroup_log) == 1


def test_the_loss_report_counts_records_and_not_references():
    """One record losing four things is one event; four records losing one thing each is four.

    Pinned because the number is the whole value of the report -- a count of references would say "4"
    for both, and a reader trying to find out how much of its file survived cannot use that.
    """
    mol, ids = _chain(5)
    mol.set_sgroups([{'type': b'DAT', 'atoms': tuple(ids[:4]), 'index': 1},
                     {'type': b'DAT', 'atoms': (ids[4],), 'index': 2}])
    with mol.edit() as e:
        for n in ids[:4]:
            e.delete_atom(n)

    assert len(mol.sgroup_log) == 1
    assert mol.sgroup_log[0].startswith('1 sgroup')


def test_a_cstate_whose_bond_dies_keeps_its_vector_tail_by_becoming_unresolved():
    """The rule the DERIVED tail run forces, and the reason it is a better rule anyway.

    A CSTATE's tail is a blob handle whose position is computed from `cstates_len`, and the blob is
    copied byte for byte by every carry.  Compacting a dead pair out would shorten the run while every
    tail stayed, so each surviving pair would read the tail of the pair before it -- tails re-paired
    with the wrong bonds, silently.  So a dead CSTATE is demoted to "unresolved", which is a state the
    model already has: it is what a CSTATE whose bond index did not resolve on read looks like.

    Two pairs, and the FIRST one dies, because that is the ordering under which a compaction would
    shift the survivor's tail rather than leave it where it was.
    """
    mol, ids = _chain(5)
    mol.set_sgroups([{'type': b'SUP', 'atoms': tuple(ids),
                      'cstates': [((ids[0], ids[1]), b'  1.0  0.0  0.0'),
                                  ((ids[2], ids[3]), b'  0.0  1.0  0.0')],
                      'index': 1}])
    with mol.edit() as e:
        e.delete_atom(ids[0])

    cstates = mol.sgroups[0]['cstates']
    assert len(cstates) == 2, 'a CSTATE was compacted out, which moves the tails'
    assert cstates[0] == (None, b'  1.0  0.0  0.0'), 'the dead pair lost its tail or kept its atoms'
    assert cstates[1] == ((ids[2], ids[3]), b'  0.0  1.0  0.0'), 'the surviving tail was re-paired'
    assert len(mol.sgroup_log) == 1, 'the demotion is still a loss and must be reported'


def test_an_unresolved_cstate_survives_a_carry_untouched():
    """(NO_REF, NO_REF) must not be fed to the index map.

    `newidx` is an array subscripted by atom index; the sentinel is 0xFFFFFFFF, so indexing with it
    reads four gigabytes past the end.  The check has to come BEFORE the subscript, and this is the
    fixture that reaches it.
    """
    mol, ids = _chain(3)
    mol.set_sgroups([{'type': b'SUP', 'atoms': tuple(ids),
                      'cstates': [(None, b'unparsed vector text')], 'index': 1}])
    with mol.edit() as e:
        e.delete_atom(ids[0])

    assert mol.sgroups[0]['cstates'] == ((None, b'unparsed vector text'),)


# ------------------------------------------------------------------------------------------------
# Invariant 3: record order is file order, and values keep their order within a key.
# ------------------------------------------------------------------------------------------------

def test_record_order_is_preserved_across_a_round_trip_and_an_edit():
    """RULES section 1.4 invariant 3.

    The records are given DESCENDING Sgroup numbers so that any sort -- by index, by type, by anything
    -- produces a different answer from the order they went in.  A fixture in ascending order would
    pass against a writer that sorted.
    """
    mol, ids = _chain(4)
    mol.set_sgroups([{'type': b'DAT', 'name': b'third', 'atoms': (ids[3],), 'index': 30},
                     {'type': b'SUP', 'name': b'second', 'atoms': (ids[2],), 'index': 20},
                     {'type': b'SRU', 'name': b'first', 'atoms': (ids[1],), 'index': 10}])
    expect = (b'third', b'second', b'first')
    assert tuple(r['name'] for r in mol.sgroups) == expect
    assert tuple(r['name'] for r in _round_trip(mol).sgroups) == expect
    with mol.edit() as e:
        e.delete_atom(ids[0])
    assert tuple(r['name'] for r in mol.sgroups) == expect


def test_a_repeated_field_keyword_keeps_every_value_and_their_order():
    """`fields` is a sequence of pairs and NOT a mapping, which is what this asserts.

    A dict keyed on the keyword would silently keep the last value of a repeated one.  Real files
    repeat keywords, so the pair sequence is the model and file order is the only order there is.

    NOT AN END-TO-END GUARANTEE, and the missing half is the reader.  No file can currently reach this
    property: the CTfile reader collapses repeated keywords into a mapping at parse time, so today the
    only caller that can exercise it is one building arena records directly, as this test does.  What is
    asserted here is that STORAGE does not lose the order -- so that when the reader stops collapsing,
    nothing on this side has to change.  Read as an end-to-end promise it would be overclaiming.
    """
    mol, ids = _chain(2)
    mol.set_sgroups([{'type': b'GEN', 'atoms': (ids[0],),
                      'fields': ((b'SMT', b'one'), (b'NATREPLACE', b'x/y'), (b'SMT', b'two'),
                                 (b'SMT', b'three'))}])
    assert _round_trip(mol).sgroups[0]['fields'] == ((b'SMT', b'one'), (b'NATREPLACE', b'x/y'),
                                                     (b'SMT', b'two'), (b'SMT', b'three'))


def test_data_lines_keep_their_order():
    mol, ids = _chain(2)
    lines = [b'line %d' % i for i in range(12)]
    mol.set_sgroups([{'type': b'DAT', 'atoms': (ids[0],), 'data': lines}])
    assert _round_trip(mol).sgroups[0]['data'] == tuple(lines)


# ------------------------------------------------------------------------------------------------
# Invariant 4: bytes in, bytes out.
# ------------------------------------------------------------------------------------------------

@pytest.mark.parametrize('field', ['type', 'subtype', 'name', 'disp_tail'])
def test_every_fixed_string_slot_round_trips_undecodable_bytes(field):
    """RULES section 1.4 invariant 4, applied to the slots it does NOT name.

    The rule is stated about `data`, and `data` is the field where an undecodable byte is most likely.
    It is not the field where a decode is most likely to be ADDED, though -- `type` and `name` look like
    identifiers, and a convenience `str` accessor on one of them is exactly the change this asserts
    against.  So all four fixed slots are held to the same standard, by the same bytes.
    """
    mol, ids = _chain(2)
    mol.set_sgroups([{'type': b'DAT', 'atoms': (ids[0],), field: RAW}])
    assert _round_trip(mol).sgroups[0][field] == RAW
    with pytest.raises(UnicodeDecodeError):
        RAW.decode('utf8')


def test_data_and_field_values_round_trip_undecodable_bytes():
    mol, ids = _chain(2)
    mol.set_sgroups([{'type': b'DAT', 'atoms': (ids[0],), 'data': [RAW, b'', RAW + RAW],
                      'fields': ((RAW, RAW),), 'log': [RAW]}])
    rec = _round_trip(mol).sgroups[0]
    assert rec['data'] == (RAW, b'', RAW + RAW), 'an empty datum between two full ones is a datum'
    assert rec['fields'] == ((RAW, RAW),)
    assert rec['log'] == (RAW,)


def test_the_title_round_trips_undecodable_bytes():
    """An SDF name line is not required to be UTF-8 either, and it is handle 0 of the same blob.

    THE PROMISE IS THE ROUND TRIP AND NOT THE TYPE: `title` is `str`, decoded with `surrogateescape`,
    and re-encoding with the same handler gives back the byte the blob holds.
    """
    mol, _ = _chain(2)
    mol.set_title(RAW)
    assert _round_trip(mol).title.encode('utf8', 'surrogateescape') == RAW


def test_a_str_title_comes_back_as_the_same_str():
    """`str` in, `str` out, and the blob still stores bytes.

    `title` decodes with `surrogateescape`, which is what lets a byte no codec accepts survive the
    decode, so the accessor is symmetric with `set_title`; `test_title.py` pins the round trip.
    """
    mol, _ = _chain(2)
    mol.set_title('caf\xe9')
    assert mol.title == 'caf\xe9'
    assert isinstance(mol.title, str)


def test_a_title_and_no_sgroups_is_a_molecule_that_still_carries_a_blob():
    """The title is handle 0, so it is the one blob user that needs no records at all.

    Worth its own test because the validation of the blob is deliberately independent of the validation
    of the records: sizing the blob check to its only caller is how an untrusted segment goes unwalked.
    """
    mol, _ = _chain(3)
    mol.set_title(b'aspirin')
    assert mol.sgroups == ()
    assert _round_trip(mol).title == 'aspirin'
    assert MoleculeContainer().title == '', 'a molecule with no blob has no title, not an error'


def test_a_title_survives_every_edit():
    mol, ids = _chain(4)
    mol.set_title(b'kept')
    with mol.edit() as e:
        e.delete_atom(ids[0])
        e.add_atom('N')
    assert mol.title == 'kept'


# ------------------------------------------------------------------------------------------------
# Invariant 5: 0xFFFF is a sentinel for three fields, so their real maximum is 0xFFFE.
# ------------------------------------------------------------------------------------------------

def test_the_greatest_real_sgroup_number_is_distinguishable_from_no_number():
    """RULES section 1.4 invariant 5 and section 6.3.

    0xFFFE stored and read back as 0xFFFE is the assertion; a field that saturated or that treated its
    own maximum as "none" would fail it.  Both numbered records are also given a numbered PARENT, so
    the sentinel is exercised in the field where a wrong answer is a broken hierarchy rather than a
    wrong label.
    """
    mol, ids = _chain(2)
    mol.set_sgroups([{'type': b'A', 'atoms': (ids[0],), 'index': 0xFFFE, 'ext_index': 0xFFFE},
                     {'type': b'B', 'atoms': (ids[1],), 'index': 0, 'parent': 0xFFFE},
                     {'type': b'C', 'atoms': ()}])
    got = _round_trip(mol).sgroups
    assert (got[0]['index'], got[0]['ext_index'], got[0]['parent']) == (0xFFFE, 0xFFFE, 0xFFFF)
    assert (got[1]['index'], got[1]['parent']) == (0, 0xFFFE)
    assert got[2]['index'] == 0xFFFF, 'an unset number is the sentinel and not 0'


def test_sgroup_number_zero_is_a_real_number():
    """V2000 `M  STY` writes the number in a three-character field where `  0` is representable.

    Nothing in CTfile forbids it, so using 0 as "none" would force a renumber-and-log on read and break
    the one thing this storage promises unconditionally -- that `index` round-trips verbatim.
    """
    mol, ids = _chain(2)
    mol.set_sgroups([{'type': b'DAT', 'atoms': (ids[0],), 'index': 0}])
    assert _round_trip(mol).sgroups[0]['index'] == 0


def test_a_number_past_the_sentinel_is_refused():
    mol, ids = _chain(2)
    for field in ('index', 'ext_index', 'parent'):
        with pytest.raises(ValueError, match='0..65534'):
            mol.set_sgroups([{'type': b'DAT', 'atoms': (ids[0],), field: 0x10000}])


def test_a_dangling_parent_is_refused_whatever_the_record_order():
    """Refused at the boundary, and in BOTH orders, because file order is not hierarchy order.

    A record may name a parent declared after it, so the check is a second pass.  A one-pass check
    would accept the forward reference and reject the backward one, which is a rule about file layout
    masquerading as a rule about hierarchy.
    """
    mol, ids = _chain(2)
    with pytest.raises(ValueError, match='parent 7'):
        mol.set_sgroups([{'type': b'DAT', 'atoms': (ids[0],), 'parent': 7}])

    # forward reference: the parent is declared second and must be accepted
    mol.set_sgroups([{'type': b'DAT', 'atoms': (ids[0],), 'index': 1, 'parent': 2},
                     {'type': b'SUP', 'atoms': (ids[1],), 'index': 2}])
    assert mol.sgroups[0]['parent'] == 2


# ------------------------------------------------------------------------------------------------
# The FIELDDISP anchor: xy_t is exactly F10.4, which is why it is the right type and not merely one
# that fits.
# ------------------------------------------------------------------------------------------------

def test_a_fielddisp_anchor_round_trips_all_four_decimals_exactly():
    """A float32 would lose the fourth decimal on a five-digit coordinate, silently and only there.

    So the values tested are not small: `-98765.4321` is the case that separates an exact fixed point
    from an approximate float, and `0.0001` is the case that separates it from an integer.
    """
    mol, ids = _chain(2)
    for x, y in ((49.5979, -3.8125), (-98765.4321, 98765.4321), (0.0001, -0.0001),
                 (0.0, 0.0), (214748.0, -214748.0)):
        mol.set_sgroups([{'type': b'DAT', 'atoms': (ids[0],), 'disp': (x, y),
                          'disp_tail': b'    DA    ALL  1       5'}])
        got = _round_trip(mol).sgroups[0]
        assert got['disp'] == pytest.approx((x, y), abs=1e-9)
        assert got['disp_tail'] == b'    DA    ALL  1       5'


def test_an_absent_anchor_is_distinguishable_from_the_origin():
    """(0, 0) is a legal anchor, which is the whole reason SGROUP_FLAG_DISP exists.

    Without the flag, "no FIELDDISP" and "a FIELDDISP at the origin" would be the same stored bytes and
    the writer would have to invent one of them.
    """
    mol, ids = _chain(2)
    mol.set_sgroups([{'type': b'DAT', 'atoms': (ids[0],), 'disp': (0.0, 0.0)},
                     {'type': b'DAT', 'atoms': (ids[1],)}])
    got = _round_trip(mol).sgroups
    assert got[0]['disp'] == (0.0, 0.0)
    assert got[1]['disp'] is None


def test_an_anchor_outside_the_fixed_point_range_is_refused():
    mol, ids = _chain(2)
    with pytest.raises(ValueError, match='214748'):
        mol.set_sgroups([{'type': b'DAT', 'atoms': (ids[0],), 'disp': (300000.0, 0.0)}])


# ------------------------------------------------------------------------------------------------
# Atom aliases, which share the storage and are a separate view.
# ------------------------------------------------------------------------------------------------

def test_an_alias_follows_its_atom_across_an_edit_and_is_not_an_sgroup():
    mol, ids = _chain(4)
    mol.set_aliases({ids[1]: 'Ph', ids[3]: b'OMe'})
    mol.set_sgroups([{'type': b'DAT', 'atoms': (ids[2],), 'name': b'not an alias'}])

    assert mol.aliases == {ids[1]: b'Ph', ids[3]: b'OMe'}
    assert len(mol.sgroups) == 1, 'aliases leaked into the S-group view'
    assert mol.sgroups[0]['name'] == b'not an alias'

    with mol.edit() as e:
        e.delete_atom(ids[0])
    assert mol.aliases == {ids[1]: b'Ph', ids[3]: b'OMe'}, 'aliases did not follow their atoms'
    assert _round_trip(mol).aliases == mol.aliases


def test_an_alias_disappears_from_the_view_when_its_atom_dies():
    """It leaves the VIEW and stays in storage, which is invariant 1 applied to an alias too.

    Stated as a test because it is the one place aliases behave unlike a dict: a label with no atom has
    nothing to label, so it is not in `aliases` -- but the record is still there and the loss is still
    reported, exactly as for every other record.
    """
    mol, ids = _chain(3)
    mol.set_aliases({ids[0]: b'Ph'})
    with mol.edit() as e:
        e.delete_atom(ids[0])
    assert mol.aliases == {}
    assert mol.sgroup_log, 'the alias was dropped from the view with no diagnostic'
    assert mol.sgroup_log == ('1 alias(es) lost the atom they label',)


def test_a_lost_alias_and_a_lost_sgroup_reference_are_two_different_log_lines():
    """One count for both would send a reader to the wrong half of the file.

    An alias is stored as a record like any other, and storage has no reason to care which it is.  A
    writer does: it emits an alias as a display label on one atom and an S-group as its own block, so
    "something lost a reference" is not actionable.  The two are counted apart and reported apart.

    THE ORDER IS ASSERTED BECAUSE IT IS PART OF THE CONTRACT.  Aliases come first because a V2000
    record puts its `A`/`V` alias lines in the atom-adjacent part and the `M  ST*` properties block
    after them, so a reader diffing this log against a file walks both in the same direction.
    """
    mol, ids = _chain(4)
    mol.set_sgroups([{'type': b'DAT', 'name': b'BATCH', 'atoms': (ids[0], ids[1]),
                      'data': [b'lot-42']}])
    mol.set_aliases({ids[0]: b'Ph', ids[1]: b'OMe'})
    with mol.edit() as e:
        e.delete_atom(ids[0])
        e.delete_atom(ids[1])
    assert mol.sgroup_log == ('2 alias(es) lost the atom they label',
                              '1 sgroup record(s) lost a reference to a deleted atom')


def test_setting_one_view_leaves_the_other_alone():
    """Three independent views over one segment, so each setter must be surgical.

    The failure this catches is the easy one to write: a setter that rebuilds the segment from its own
    argument and drops whatever the other two views held.
    """
    mol, ids = _chain(3)
    mol.set_title(b'title')
    mol.set_sgroups([{'type': b'DAT', 'atoms': (ids[0],), 'name': b'rec'}])
    mol.set_aliases({ids[1]: b'Ph'})
    assert (mol.title, len(mol.sgroups), mol.aliases) == ('title', 1, {ids[1]: b'Ph'})

    mol.set_title(b'retitled')
    assert (len(mol.sgroups), mol.aliases) == (1, {ids[1]: b'Ph'})
    mol.set_sgroups([{'type': b'SUP', 'atoms': (ids[2],)}])
    assert (mol.title, mol.aliases) == ('retitled', {ids[1]: b'Ph'})
    mol.set_aliases({})
    assert (mol.title, len(mol.sgroups)) == ('retitled', 1)
    assert mol.sgroups[0]['type'] == b'SUP'


def test_an_alias_labels_exactly_one_atom():
    mol, ids = _chain(3)
    with pytest.raises(ValueError, match='exactly one atom'):
        mol.set_sgroups([{'type': b'', 'atoms': (ids[0], ids[1]), '_alias': True}])


# ------------------------------------------------------------------------------------------------
# The boundary: what the container refuses, and what it never invents.
# ------------------------------------------------------------------------------------------------

def test_an_unknown_key_is_an_error_and_not_ignored():
    """A misspelled key that silently did nothing is a lost reference set with no diagnostic."""
    mol, ids = _chain(2)
    with pytest.raises(ValueError, match='patom'):
        mol.set_sgroups([{'type': b'DAT', 'atoms': (ids[0],), 'patom': (ids[1],)}])


def test_what_sgroups_returns_is_what_set_sgroups_accepts():
    """The read shape and the write shape are ONE shape, which is a property and not a coincidence.

    A read-only key -- `_alias`, here -- is how the two quietly stop being one, because the obvious
    round trip then fails on the caller's own output.
    """
    mol, ids = _chain(4)
    mol.set_sgroups([{'type': b'DAT', 'name': b'x', 'atoms': (ids[0], ids[1]),
                      'patoms': (ids[0],), 'bonds': ((ids[0], ids[1]),),
                      'cstates': [((ids[2], ids[3]), b'tail'), (None, b'raw')],
                      'data': [b'd'], 'fields': ((b'k', b'v'),), 'log': [b'l'],
                      'index': 3, 'ext_index': 4, 'parent': 0xFFFF,
                      'disp': (1.5, -2.5), 'disp_tail': b'tail'}])
    mol.set_aliases({ids[3]: b'Ph'})
    before = mol.sgroups
    mol.set_sgroups(list(before))
    assert mol.sgroups == before
    assert mol.aliases == {ids[3]: b'Ph'}, 'a round trip through set_sgroups ate the aliases'


def test_a_reference_to_an_atom_this_molecule_does_not_have_is_refused():
    """Refused with a KeyError from the id map, and that is the right shape: it is a caller bug.

    NOT the same thing as a reference that DIES, which is data and is reported.  A caller naming an atom
    that never existed has made a mistake, and quietly dropping it would hide it.
    """
    mol, ids = _chain(2)
    with pytest.raises(KeyError):
        mol.set_sgroups([{'type': b'DAT', 'atoms': (9999,)}])


def test_patoms_are_carried_and_remapped_like_atoms():
    """The parent-atom subset of a MUL group, which is the one reference list with no pairs.

    Easy to leave out of a carry loop, because nothing else in the record depends on it -- so it gets
    its own assertion rather than riding on the atoms one.
    """
    mol, ids = _chain(5)
    mol.set_sgroups([{'type': b'MUL', 'atoms': tuple(ids[1:]), 'patoms': (ids[1], ids[2]),
                      'index': 1}])
    with mol.edit() as e:
        e.delete_atom(ids[0])
    assert mol.sgroups[0]['patoms'] == (ids[1], ids[2])
    with mol.edit() as e:
        e.delete_atom(ids[1])
    assert mol.sgroups[0]['patoms'] == (ids[2],)


def test_sgroups_survive_a_remap():
    """`remap` relabels stable ids without touching an index, so the references must move with them."""
    mol, ids = _chain(3)
    mol.set_title(b'kept')
    mol.set_sgroups([{'type': b'DAT', 'atoms': (ids[0], ids[2]),
                      'bonds': ((ids[0], ids[1]),), 'index': 1}])
    mol.set_aliases({ids[0]: b'Ph'})
    mol.remap({ids[0]: 100, ids[1]: 200, ids[2]: 300})
    rec = mol.sgroups[0]
    assert rec['atoms'] == (100, 300)
    assert rec['bonds'] == ((100, 200),)
    assert mol.aliases == {100: b'Ph'}
    assert mol.title == 'kept'


def test_to_bytes_is_stable_across_reads_of_a_molecule_carrying_sgroups():
    """`to_bytes()` IS a molecule identity in v4, so a read must not move a byte of it.

    The S-group segments are persistent, so they are inside that identity -- which means every read
    that builds a derived cache has to leave them alone, and this is the fixture that would notice.
    """
    mol, ids = _chain(6)
    mol.set_title(b'identity')
    mol.set_sgroups([{'type': b'SRU', 'atoms': tuple(ids), 'bonds': ((ids[0], ids[1]),),
                      'data': [RAW], 'index': 1}])
    mol.set_aliases({ids[0]: b'Ph'})
    before = mol.to_bytes()
    mol.rings
    mol.atoms_order
    mol.connected_components
    str(mol)
    assert mol.to_bytes() == before
    assert _round_trip(mol).to_bytes() == before


def test_a_molecule_with_no_sgroups_pays_nothing_for_the_feature():
    """The three segments are absent, not empty, on a molecule that has none of them.

    Trailing empty table entries are not written, so a plain molecule's header is unchanged by this
    release -- which is worth an assertion because the alternative costs 24 bytes on every molecule
    ever stored, and would only ever be noticed as a size regression.
    """
    plain, _ = _chain(6)
    plain_len = len(plain.to_bytes())
    titled, _ = _chain(6)
    titled.set_title(b'')
    assert plain.sgroups == () and plain.aliases == {} and plain.title == ''
    assert len(titled.to_bytes()) > plain_len, 'an empty title cost nothing, so it stored nothing'
