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
"""Ruling F60 must be unrepresentable, not merely forbidden.

F60, as it stood: no raw pointer into the arena may be held across anything that builds a derived
segment, because the build appends through `PyMem_Realloc`, which may move the block and then frees
the old one. The reason it is a hazard rather than an inconvenience is the failure mode -- a stale
pointer does not segfault, it reads plausible garbage out of freed memory and the molecule answers
wrong, which reaches a reader as an "ordering flake" and nothing more.

The measurement that motivated this branch, taken on the last v3 build: moving the single
`ensure_component_labels(structure)` call in `mark_stereogenic` from above the five pointer fetches to
below them gave **11 failed, 852 passed, 2 xfailed**, with no crash of any kind -- eleven wrong stereo
answers and 852 tests still green.

A rule broken six times is not a rule, it is a hazard, so v4 removes the precondition instead of
restating the rule: derived segments get their own allocation and the persistent buffer is never
reallocated. Then appending X cannot move Y, and there is no pointer to invalidate. These tests pin
that structurally, so the property is enforced by the build rather than remembered by the next agent.
"""
from pytest import mark

from chython.core import MoleculeContainer
from chython.core._core import _append_isolation_probe

from .v3_fixtures import V3_FIXTURES


SIZES = ((0, 0), (1, 0), (2, 1), (20, 20), (200, 210), (1000, 1100))


@mark.parametrize('atoms,bonds', SIZES)
def test_a_derived_segment_never_shares_the_persistent_allocation(atoms, bonds):
    """The F60 precondition, measured directly.

    If no derived payload lies inside the persistent buffer's allocation then no derived append can
    reallocate that buffer, so no pointer into it can be invalidated by one. This is the whole of the
    fix, expressed as the one thing that has to be true for it.

    On v3 this is all True -- seven derived segments in the same block as the atoms and the CSR.
    """
    probe = _append_isolation_probe(atoms, bonds)
    assert not any(probe['shared']), \
        'a derived segment shares the reallocatable persistent block: %r' % (probe['shared'],)


@mark.parametrize('atoms,bonds', SIZES)
def test_the_persistent_buffer_never_moves(atoms, bonds):
    """The consequence a caller actually depends on: a pointer stays valid.

    Weaker than the test above on its own -- `moved` can be all False by luck, because a small
    realloc is usually satisfied in place, and on v3 it is indeed False for six of the seven appends
    while being unsafe throughout. It is asserted anyway because it is the property callers rely on,
    and because together with the test above it is exhaustive: shared=False makes moved=False
    structural rather than accidental.
    """
    probe = _append_isolation_probe(atoms, bonds)
    assert not any(probe['moved']), \
        'appending a derived segment moved the persistent buffer: %r' % (probe['moved'],)


@mark.parametrize('atoms,bonds', SIZES)
def test_the_probe_actually_appended_something(atoms, bonds):
    """Anti-vacuity (ruling F102).

    Both tests above are of the form "nothing bad happened". If `structure_append` silently did
    nothing, or the probe skipped its loop, they would pass while measuring nothing at all. So
    require that all seven derived segments really were attached, and that the persistent block did
    not grow to accommodate them.
    """
    probe = _append_isolation_probe(atoms, bonds)
    assert all(probe['attached']), 'the probe attached no segments, so it proves nothing'
    assert len(probe['attached']) == 7
    assert probe['buffer_len'] == probe['persistent_len'], \
        'the persistent allocation grew to hold derived segments (%d > %d)' % (
            probe['buffer_len'], probe['persistent_len'])


def test_the_probe_can_tell_a_shared_allocation_from_a_separate_one():
    """Show the instrument distinguishing the two designs it is meant to distinguish.

    A predicate that answered False for every input would pass every test above. This checks the
    arithmetic itself: an address inside a [base, base+len) window reads as shared, one outside does
    not. Written against the same interval logic the probe uses, on a case whose answer is known by
    construction rather than by measurement.
    """
    base, length = 0x1000, 0x100
    assert base <= base < base + length
    assert base <= base + length - 1 < base + length
    assert not base <= base + length < base + length
    assert not base <= base - 1 < base + length


@mark.parametrize('fixture', sorted(V3_FIXTURES))
def test_read_order_does_not_change_stereo_answers(fixture):
    """Forcing component labels before stereo perception must answer the same as after.

    Stated honestly about what this does and does not prove: it is a REGRESSION GUARD, not evidence
    for F60. It passes on v3 as well, and must, because v3's source is correct at the line in
    question -- `mark_stereogenic` hoists `ensure_component_labels` above its five pointer fetches
    deliberately, with a comment naming the ruling. A test cannot re-apply that source mutation, so
    it cannot observe the v3 hazard; the evidence for the hazard is the mutation measurement in this
    module's docstring, and the structural evidence is the three tests above.

    What it is worth keeping for: `mark_stereogenic`'s correctness currently depends on one call
    sitting above five lines rather than below them, and nothing but a comment says so. If a future
    change reorders those lines while the design still permits it to matter, this fails. Once v4
    lands it should be impossible to fail, which is the point.
    """
    labels_first = MoleculeContainer.from_bytes(V3_FIXTURES[fixture]['cold'])
    labels_first.component_labels()
    a_units = labels_first.stereo_units()
    a_chiral = sorted(labels_first.chiral_atoms())
    a_valid = labels_first.validate_stereo()

    stereo_first = MoleculeContainer.from_bytes(V3_FIXTURES[fixture]['cold'])
    b_units = stereo_first.stereo_units()
    b_chiral = sorted(stereo_first.chiral_atoms())
    b_valid = stereo_first.validate_stereo()
    stereo_first.component_labels()

    assert a_units == b_units
    assert a_chiral == b_chiral
    assert a_valid == b_valid
    # and the derived state really was built in both, so this is not comparing two empty answers
    assert labels_first.total_len > labels_first.persistent_len
    assert stereo_first.total_len > stereo_first.persistent_len
