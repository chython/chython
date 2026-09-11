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
"""`to_bytes()` must be a molecule identity: no READ may change the serialised bytes.

The v3 arena could not promise this. Building a derived segment wrote its table entry and bumped
`total_len`, and both of those live inside the persistent prefix that `to_bytes()` returns -- so the
v3 docstring had to tell callers to hash `bytes[128:]` instead. The v4 arena keeps derived segments
out of the buffer entirely, which is what makes the whole slice stable.

There are two independent ways the promise can break, and this module tests both because only the
first one is fixed by moving the derived table:

  1. A read writes the HEADER -- a segment table entry, `total_len`, or a flag bit. This is the v3
     defect, and `test_derived_build_does_not_touch_the_header` is the direct measurement of it.
  2. A read writes a PERSISTENT SEGMENT. This is subtler and survives the format change untouched,
     because two pieces of derived information are deliberately stored inside persistent segments:
     `HE_IN_RING` in `halfedge_t.flags` (SEG_CSR_EDGE) and the in-ring / ring-count fields in
     `atom_t` (SEG_ATOMS). Those bytes are only safe because ring perception is EAGER -- it runs
     inside `rebuild_derived` at edit-exit and again on ingest, never lazily on first read. Nothing
     in the format enforces that, so `test_to_bytes_stable_under_every_read` is what stops a future
     change from making ring perception lazy and silently breaking serialised identity again.

Both are measured twice: once over the seven v3 fixtures, and once over five AROMATIC molecules at
the bottom of the file. The second set exists because the fixtures are v3 bytes and v3 refused order
4, so no fixture can carry `HE_AROMATIC` -- and `HE_AROMATIC` sits in the same persistent flags byte
as `HE_IN_RING`, which means hazard 2 has a second occupant that the fixture corpus is structurally
unable to see.
"""
from base64 import b64decode
from pytest import mark, raises

from chython.core import MoleculeContainer

from .v3_fixtures import V3_FIXTURES


# Every public read that can build a derived segment, or that reaches one that does. Named
# individually in three groups rather than swept from dir(): a sweep cannot tell a no-argument read
# from one taking a stable id (every Cython method is the same type), and it would drag in mutators.
NO_ARG_READS = (
    # the lazy stereo unit table (SEG_STEREO_UNIT)
    'stereo_units', 'stereogenic_units', 'chiral_atoms', 'chiral_bonds', 'validate_stereo',
    # canonical labelling, which reaches component labels and the feature words
    'canonical_order', 'automorphism_orbits', 'is_asymmetric',
    # the lazy component labels (SEG_COMPONENT_LABEL)
    'component_labels',
    # stereo groups and their canonical form
    'stereo_groups', 'canonical_stereo_groups', 'canonical_stereo_group_ambiguities', 'wedges',
)

PROPERTY_READS = (
    'stereo_truncated', 'atoms_order', 'atoms_order_classes', '_union_feature_words',
    'connected_components', 'connected_components_count',
    # ring perception -- writes HE_IN_RING and the atom ring fields, both PERSISTENT
    'rings', 'sssr', 'rings_count',
    'has_stereo_groups', 'has_coordinates', 'element_counts',
)

ATOM_READS = (
    'is_chiral', 'features_of', 'edge_words_of', 'unit_of', 'in_ring_of', 'ring_sizes_of',
    'ring_sizes_word_of', 'ring_count_of', 'macrocycle_of', 'hybridization_of',
    'heteroatoms_of', 'degree_of', 'total_h_of', 'parity_of', 'stereo_of', 'stereo_group_of',
)

BOND_READS = ('order_of', 'bond_in_ring', 'shares_ring', 'wedge_of')

READ_COUNT = len(NO_ARG_READS) + len(PROPERTY_READS) + len(ATOM_READS) + len(BOND_READS)


def _exercise(mol):
    """Call every derived-building read. Returns how many distinct reads actually ran."""
    n = 0
    for name in NO_ARG_READS:
        getattr(mol, name)()
        n += 1
    for name in PROPERTY_READS:
        getattr(mol, name)
        n += 1
    sids = list(mol.atom_numbers)
    for name in ATOM_READS:
        method = getattr(mol, name)
        for sid in sids:
            method(sid)
        n += 1
    for name in BOND_READS:
        method = getattr(mol, name)
        for a in sids:
            for b in mol.neighbors_of(a):
                method(a, b)
        n += 1
    return n


def _load(fixture, warmth):
    return MoleculeContainer.from_bytes(V3_FIXTURES[fixture][warmth])


def _payload(data):
    """The bytes after the header -- the first persistent SEGMENT onwards.

    Version-aware on purpose, so that the persistent-segment test measures the same region on a v3
    buffer and on a later one and can be run against each. All are 24 fixed bytes plus
    `seg_count` entries of 8; v3 always spends 13 for a 128-byte header, while a later buffer
    spends one past its highest used segment and is usually shorter. Reading `seg_count` out of the
    buffer is therefore load-bearing and not tidiness -- and the fixtures here are the one
    population whose header IS 128 bytes after normalisation, since a v3 buffer keeps its table so
    that no payload has to move. v3's `seg_count` field is v4's; only the four bytes at 20..23
    differ in meaning (v3's `total_len`, v4's `seg_count` plus a reserved half-word). Versions 5 and
    6 have version 4's header shape exactly -- their changes are a new segment id and a narrower
    conformer record, neither of them a header field -- so one arm serves all three.
    """
    version = int.from_bytes(data[4:6], 'little')
    if version == 3:
        return data[128:]
    elif version in (4, 5, 6):
        seg_count = int.from_bytes(data[20:22], 'little')
        return data[24 + 8 * seg_count:]
    raise AssertionError('unknown arena version %d' % version)


ALL = sorted(V3_FIXTURES)


@mark.parametrize('fixture', ALL)
@mark.parametrize('warmth', ('cold', 'warm'))
def test_to_bytes_stable_under_every_read(fixture, warmth):
    """The headline promise: reads do not change the serialised bytes."""
    mol = _load(fixture, warmth)
    before = mol.to_bytes()
    assert _exercise(mol) == READ_COUNT, 'the read set did not run'
    assert mol.to_bytes() == before


@mark.parametrize('fixture', ALL)
def test_the_reads_really_build_derived_segments(fixture):
    """Anti-vacuity for the test above (ruling F102).

    Without this, a build in which every read above was a no-op -- or in which the derived caches
    were computed eagerly and never appended -- would pass `test_to_bytes_stable_under_every_read`
    for the wrong reason. So assert that the reads genuinely grew derived state, and that they grew
    it OUTSIDE the persistent prefix.
    """
    mol = _load(fixture, 'cold')
    persistent = mol.persistent_len
    grown_before = mol.total_len
    _exercise(mol)
    assert mol.total_len > grown_before, 'no derived segment was built, so the test above is vacuous'
    assert mol.persistent_len == persistent, 'a read resized the persistent prefix'
    assert len(mol.to_bytes()) == persistent, 'to_bytes() is not exactly the persistent prefix'


@mark.parametrize('fixture', ALL)
def test_derived_build_does_not_touch_the_header(fixture):
    """The v3 defect, measured directly rather than argued.

    In v3 the first lazy read wrote a segment table entry and `total_len`, both inside the first 128
    bytes. This asserts the header is untouched, which is a strictly stronger statement than the
    whole-buffer comparison above -- if some persistent segment ALSO changes, this test still
    localises the damage to the header or not.
    """
    mol = _load(fixture, 'cold')
    before = mol.to_bytes()
    _exercise(mol)
    after = mol.to_bytes()
    assert after[:32] == before[:32], 'a read wrote the fixed header'
    assert after == before


@mark.parametrize('fixture', ALL)
def test_a_read_does_not_change_a_persistent_segment(fixture):
    """Ring perception writes HE_IN_RING and the atom ring fields, both PERSISTENT.

    That is only safe while ring perception stays eager. This is the test that fails the day it is
    made lazy, which is the point of writing it separately from the whole-buffer comparison: the
    failure message says which region moved.
    """
    mol = _load(fixture, 'cold')
    before = _payload(mol.to_bytes())
    # ring perception first and alone, so that a failure here cannot be blamed on stereo
    mol.rings_count
    mol.sssr
    for sid in mol.atom_numbers:
        mol.in_ring_of(sid)
        mol.ring_sizes_of(sid)
    assert _payload(mol.to_bytes()) == before, 'ring perception wrote a persistent segment'
    # and now the stereo half, which writes the parity into SEG_PARITY on a MUTATING path only
    mol.stereo_units()
    mol.validate_stereo()
    assert _payload(mol.to_bytes()) == before, 'stereo perception wrote a persistent segment'
    # and everything else
    _exercise(mol)
    assert _payload(mol.to_bytes()) == before, 'some read wrote a persistent segment'


@mark.parametrize('fixture', ALL)
def test_repeated_reads_are_idempotent(fixture):
    """A second pass over the same reads must not differ from the first.

    A cache that rebuilds itself on every access would pass the single-pass tests while doing the
    work every time; this catches the version of that bug that also perturbs bytes.
    """
    mol = _load(fixture, 'cold')
    _exercise(mol)
    once = mol.to_bytes()
    grown = mol.total_len
    _exercise(mol)
    assert mol.to_bytes() == once
    assert mol.total_len == grown, 'a derived segment was rebuilt and re-appended'


@mark.parametrize('fixture', ALL)
def test_identity_survives_a_round_trip_through_a_warmed_molecule(fixture):
    """Serialising a read-warmed molecule must give the same bytes as serialising a cold one.

    This is the property a dedup key needs, and the one v3 could not offer: two Structures holding
    the same molecule must serialise identically regardless of what has been asked of them. It is
    NOT implied by per-molecule stability -- two molecules could each be stable at different bytes.
    """
    cold = _load(fixture, 'cold')
    warm = _load(fixture, 'cold')
    _exercise(warm)
    assert warm.to_bytes() == cold.to_bytes()
    # and the same across the copy that shares the arena, plus a clone that does not
    assert warm.copy().to_bytes() == cold.to_bytes()


def test_the_harness_can_fail():
    """Ruling F102: show the comparison catching a molecule that really did change.

    Every assertion above is of the form "these bytes did not move". A harness that could not
    observe bytes moving would pass all of them while measuring nothing, so make it observe one.
    A MUTATION is expected to change the bytes -- that is the control.
    """
    mol = _load(ALL[0], 'cold')
    before = mol.to_bytes()
    with mol.edit() as m:
        m.set_charge(mol.atom_numbers[0], 1)
    assert mol.to_bytes() != before, 'the comparison cannot see a change, so it proves nothing'


def test_to_bytes_is_not_the_whole_allocation():
    """`to_bytes()` returns the persistent prefix, and derived bytes are not in it.

    Stated as a test because the size relation is the whole design: if `to_bytes()` ever returned
    `total_len` bytes again, every identity test above would still pass on a cold molecule and fail
    only on a warm one, which is exactly the v3 failure mode.
    """
    mol = _load('a_ring', 'cold')
    _exercise(mol)
    assert mol.total_len > mol.persistent_len
    assert len(mol.to_bytes()) == mol.persistent_len
    assert bytes(mol.persistent_view) == mol.to_bytes()


# ── the same promise on molecules the v3 fixtures cannot contain ─────────────────────────────────
#
# Every fixture above is v3 bytes, and v3 refused order 4, so nothing in `V3_FIXTURES` carries an
# aromatic bond. That left the identity promise UNMEASURED on exactly the molecules arena v4 added,
# and the risk is specific rather than theoretical: `HE_AROMATIC` lives in the halfedge flags field,
# which is PERSISTENT and which `perceive_rings` also writes -- `HE_IN_RING` is the neighbouring bit.
# A ring pass that rebuilt that byte from its own perception instead of OR-ing into it would clear
# the aromatic bit on first read, silently, and only on a molecule no fixture holds.

AROMATIC_CASES = {
    # benzene: the flags byte takes HE_AROMATIC from the input and HE_IN_RING from perception, so
    # every half-edge here is a place the two writers could collide
    'benzene': ('C' * 6, [(i, (i + 1) % 6, 4) for i in range(6)]),
    # naphthalene: the fusion bonds carry two ring memberships as well as the aromatic flag
    'naphthalene': ('C' * 10, [(0, 1, 4), (1, 2, 4), (2, 3, 4), (3, 4, 4), (4, 5, 4), (5, 0, 4),
                               (4, 6, 4), (6, 7, 4), (7, 8, 4), (8, 9, 4), (9, 5, 4)]),
    # a mixed record: one ring aromatic, one Kekule, joined by an acyclic single bond. The case a
    # three-state flag could not describe, and where a normalising read would show up
    'mixed biaryl': ('C' * 12, [(i, (i + 1) % 6, 4) for i in range(6)] +
                     [(6 + i, 6 + (i + 1) % 6, 2 if i % 2 == 0 else 1) for i in range(6)] +
                     [(0, 6, 1)]),
    # an aromatic bond in NO ring: chemical nonsense, and therefore the sharpest test of whether
    # perception rewrites what the caller stored
    'acyclic aromatic': ('CC', [(0, 1, 4)]),
    # aromatic thiophene: the heteroatom whose electron budget the stereo pass reads
    'thiophene': ('SCCCC', [(0, 1, 4), (1, 2, 4), (2, 3, 4), (3, 4, 4), (4, 0, 4)]),
}


def _aromatic(case):
    atoms, bonds = AROMATIC_CASES[case]
    mol = MoleculeContainer()
    with mol.edit():
        sids = [mol.add_atom(e) for e in atoms]
        for i, j, order in bonds:
            mol.add_bond(sids[i], sids[j], order)
    return mol


AROMATIC = sorted(AROMATIC_CASES)


@mark.parametrize('case', AROMATIC)
def test_to_bytes_stable_under_every_read_on_an_aromatic_molecule(case):
    """The headline promise, extended to the molecules arena v4 made storable."""
    mol = _aromatic(case)
    before = mol.to_bytes()
    aromatic_before = mol.aromatic_bond_count
    assert aromatic_before, 'the case carries no aromatic bond, so it measures nothing new'
    assert _exercise(mol) == READ_COUNT, 'the read set did not run'
    assert mol.to_bytes() == before
    assert mol.aromatic_bond_count == aromatic_before


@mark.parametrize('case', AROMATIC)
def test_the_reads_build_derived_state_on_an_aromatic_molecule_too(case):
    """Anti-vacuity: the reads must actually grow something here as well (ruling F102)."""
    mol = _aromatic(case)
    persistent = mol.persistent_len
    grown = mol.total_len
    _exercise(mol)
    assert mol.total_len > grown, 'no derived segment was built, so the test above is vacuous'
    assert mol.persistent_len == persistent
    assert len(mol.to_bytes()) == persistent


@mark.parametrize('case', AROMATIC)
def test_aromatic_identity_survives_serialisation_and_a_second_warming(case):
    """Cold bytes, warm bytes and a round trip through both must be one byte string.

    `from_bytes` clears the derived segments and re-derives, so this is where a flags byte rebuilt
    from perception rather than read from the buffer would part company with the original -- and it
    is also where the `(order == 4) == HE_AROMATIC` agreement check would refuse a buffer a writer had
    corrupted, which is the loud failure rather than the quiet one.
    """
    cold = _aromatic(case)
    raw = cold.to_bytes()
    _exercise(cold)
    assert cold.to_bytes() == raw
    warm = MoleculeContainer.from_bytes(raw)
    _exercise(warm)
    assert warm.to_bytes() == raw
    assert MoleculeContainer.from_bytes(warm.to_bytes()).to_bytes() == raw
    assert warm.copy().to_bytes() == raw


def test_the_aromatic_harness_can_fail():
    """Ruling F102, for this section: show the comparison seeing an aromatic bond change.

    The control is `kekule()`, a MUTATION whose whole job is to change these bytes. If the comparison
    could not observe it, the three tests above would prove nothing.
    """
    mol = _aromatic('benzene')
    before = mol.to_bytes()
    assert mol.kekule().changed
    assert mol.to_bytes() != before
    assert mol.aromatic_bond_count == 0
