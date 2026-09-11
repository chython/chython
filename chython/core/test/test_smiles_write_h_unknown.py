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
"""H_UNKNOWN on the way OUT: what a string can say about a count nobody stated.

The arena's side of the sentinel is `test_h_unknown.py`.  This file is the writer's, and it exists
because SMILES HAS NO SPELLING FOR "UNSTATED" ANYWHERE:

* inside brackets an absent H term means ZERO -- `[13C]` is a hydrogen-free carbon, exactly as
  `[13CH0]` is -- so a bracketed atom cannot decline to answer;
* a bare `C` does mean "the reader derives it", which is the closest thing to the truth, and it is
  available only when no other property forces the bracket.

So the writer omits the term in both cases and REPORTS the atom in `smw_traversal`'s `unknown_h`,
whichever spelling it got.  The first assertion of the file is the one that matters most: the number
15 never reaches the string.  `[CH15]O` is not a hydrogen count anyone can read -- it is the raw
nibble printed as a number.

The second half is stereo, where the sentinel is not merely unspellable but ACTIVE: a sign is a
statement about four POSITIONS in the written order, an implicit hydrogen occupies the position it is
written in, and an unknown count does not say whether there is one.  Those signs are refused and
reported in `lost`.  Refused only where the count can actually move a position, which is why
`halomethane` (four heavy neighbours, no room for a hydrogen at all) keeps its `@` and
`bromochlorofluoromethane` loses it.
"""
from itertools import permutations

from pytest import mark

from chython.core import H_UNKNOWN, MoleculeContainer
from chython.core._core import smw_traversal, write_smiles


# ------------------------------------------------------------------------------------------------
# FIXTURE PLUMBING, the same shape as the other writer test files': atoms as (element, implicit_h),
# bonds as (i, j, order), and `order` a creation order over the atom list so a claim can be swept.
def build(atoms, bonds, order=None):
    m = MoleculeContainer()
    sids = {}
    for j in (range(len(atoms)) if order is None else order):
        element, hydrogens = atoms[j]
        sids[j] = m.add_atom(element, implicit_h=hydrogens)
    for a, b, o in bonds:
        m.add_bond(sids[a], sids[b], o)
    return m, sids


def configure(m, sid, frame, want):
    """Store the parity that makes `translate_stereo(sid, frame)` answer `want`.

    A STORED PARITY IS NOT A CONFIGURATION -- it is a configuration relative to the atom's refs order,
    which the CREATION order decides.  So a sweep calling `set_parity(sid, 2)` in every creation order
    is sweeping different MOLECULES, and the moment perception stopped refusing `halomethane` that
    sweep started reporting two strings for what looked like one input.  Asking the production
    `translate_stereo` which value lands the wanted arrangement of atom identities is the technique
    `test_smiles_write_stereo.py` uses, for the same reason: deriving it here would reimplement the
    arithmetic under test.
    """
    for parity in (1, 2):
        m.set_parity(sid, parity)
        if m.translate_stereo(sid, frame) == want:
            return parity
    raise AssertionError('neither parity gives %r in frame %r' % (want, frame))


def sign_in(smiles):
    if '@@' in smiles:
        return 2
    if '@' in smiles:
        return 1
    return 0


# Methanol whose carbon states nothing.  Nothing else forces a bracket, so the bare spelling is
# reachable and the count comes back from the valence model.
METHANOL = ([(6, H_UNKNOWN), (8, 1)], [(0, 1, 1)])
# The same carbon with a second reason for brackets, one per row: (property, value, spec, string).
# Every one of them ends up stating ZERO hydrogens, and every one is a loss.
FORCED = [('isotope', 'isotope', 13, '', '[13C]O'),
          ('charge', 'charge', 1, '', '[N+]O'),
          ('radical', 'radical', True, '', '[C]O |^1:0|'),
          ('map number', 'map', 7, 'm', '[C:7]O')]

# CHFClBr with the hydrogen unstated: three heavy neighbours, so the frame has one unnamed direction
# and the arena will not say whether it is a hydrogen.
BROMOCHLOROFLUOROMETHANE = ([(6, H_UNKNOWN), (9, 0), (17, 0), (35, 0)],
                            [(0, 1, 1), (0, 2, 1), (0, 3, 1)])
# CFClBrI with the count unstated: FOUR heavy neighbours leave no room for a hydrogen, so the
# sentinel says nothing the frame needed and the sign survives.
HALOMETHANE = ([(6, H_UNKNOWN), (9, 0), (17, 0), (35, 0), (53, 0)],
               [(0, 1, 1), (0, 2, 1), (0, 3, 1), (0, 4, 1)])
# 1,3-difluoroallene, one terminal's hydrogen unstated.  The axial sign is over the TERMINALS'
# directions, so this is the same defect one atom further from the anchor.
DIFLUOROALLENE = ([(6, 1), (6, 0), (6, H_UNKNOWN), (9, 0), (9, 0)],
                  [(0, 1, 2), (1, 2, 2), (0, 3, 1), (2, 4, 1)])
DIFLUOROALLENE_KNOWN = ([(6, 1), (6, 0), (6, 1), (9, 0), (9, 0)],
                        [(0, 1, 2), (1, 2, 2), (0, 3, 1), (2, 4, 1)])
# A tetrasubstituted allene: both terminals carry two heavy substituents, so no terminal needs a
# hydrogen position and an unknown count on one of them changes nothing.
TETRAHALOALLENE = ([(6, 0), (6, 0), (6, H_UNKNOWN), (35, 0), (9, 0), (17, 0), (53, 0)],
                   [(0, 1, 2), (1, 2, 2), (0, 3, 1), (0, 4, 1), (2, 5, 1), (2, 6, 1)])


# ------------------------------------------------------------------------------------------------
# THE NUMBER THAT MUST NOT APPEAR.
def test_the_nibble_is_never_printed_as_a_count():
    """`[CH15]O` was the output, and it is the reason this file exists.

    Not "wrong by one" -- 15 is the sentinel's bit pattern, and printing it invents a pentadecavalent
    carbon out of a flag.  Asserting the absence of the substring as well as the whole string, because
    the substring is what would survive a partial fix somewhere else in the token.
    """
    m, _ = build(*METHANOL)
    assert write_smiles(m) == 'CO'
    assert 'H15' not in write_smiles(m)
    assert '15' not in write_smiles(m)


def test_no_spec_prints_it_either():
    """Every accepted spec, since the H term is written from more than one branch.

    `h` is in the list on purpose: its contract is "state the count explicitly", and the one atom
    here has no count to state.
    """
    m, _ = build(*METHANOL)
    for spec in ('', 'h', 'A', 'm', '!s', '!b', '!z', '!x', 'a', 'i', 'hA', '!bh'):
        assert '15' not in write_smiles(m, spec), spec


# ------------------------------------------------------------------------------------------------
# THE TWO SPELLINGS.
def test_a_bare_atom_leaves_the_count_to_the_reader():
    """The good case, and the only one that loses nothing a reader would notice: `CO`.

    A bare symbol means "derive the count from the valence model", which is the same thing the
    molecule says -- so the string's answer is the best available one.  It is still in `unknown_h`,
    because "derive it" and "nobody knows" are not the same claim and a caller comparing the two
    molecules will find a carbon with three hydrogens where this one has none.
    """
    m, sids = build(*METHANOL)
    assert write_smiles(m) == 'CO'
    assert smw_traversal(m)['unknown_h'] == (sids[0],)


def test_forcing_hydrogens_cannot_state_what_nobody_knows():
    """`h` brackets every atom and states its count -- except this one.

    Bracketing it would state a ZERO, which is the one number the molecule rules out saying.  So `h`
    leaves the atom bare, and that is a DELIBERATE hole in the key's contract rather than a bug in
    it: the alternative is a string that lies.
    """
    m, sids = build(*METHANOL)
    assert write_smiles(m, 'h') == 'C[OH]'
    assert smw_traversal(m, 'h')['unknown_h'] == (sids[0],)


@mark.parametrize('name,prop,value,spec,expected', FORCED)
def test_another_property_forces_the_bracket_and_the_string_then_says_zero(name, prop, value,
                                                                           spec, expected):
    """The unavoidable loss, one row per reason: the H term is omitted and a bracket reads that as 0.

    Nothing else can be done in this notation -- there is no `[13CH?]` -- so the whole of the
    writer's obligation is to REPORT it, which is the second assertion.  The charge row uses nitrogen
    because a carbocation with an unstated hydrogen count is not a molecule anyone holds; the reason
    under test is the bracket, and every row reaches it by its own door.
    """
    m, sids = build([(7 if prop == 'charge' else 6, H_UNKNOWN), (8, 1)], [(0, 1, 1)])
    if prop == 'isotope':
        m.set_isotope(sids[0], value)
    elif prop == 'charge':
        m.set_charge(sids[0], value)
    elif prop == 'radical':
        m.set_radical(sids[0], value)
    else:
        m.set_map_number(sids[0], value)
    assert write_smiles(m, spec) == expected
    assert smw_traversal(m, spec)['unknown_h'] == (sids[0],)


def test_the_report_is_in_emission_order_and_counts_the_same_atoms_as_the_arena():
    """`unknown_h` is a report about atoms, so it is ordered like `order` and agrees with
    `unknown_h_count`.

    The arena's counter is the independent side of this: two mechanisms count the same sentinel, one
    over slots and one over the emission sequence, and a writer that lost an atom on the way would
    disagree with it.
    """
    m, sids = build([(6, H_UNKNOWN), (8, 1), (6, H_UNKNOWN)],
                    [(0, 1, 1), (1, 2, 1)])
    probe = smw_traversal(m)
    assert m.unknown_h_count == 2
    assert set(probe['unknown_h']) == {sids[0], sids[2]}
    assert len(probe['unknown_h']) == 2
    positions = [probe['order'].index(sid) for sid in probe['unknown_h']]
    assert positions == sorted(positions)


def test_an_atom_with_a_stated_count_is_not_reported():
    """The control.  Zero is a count."""
    m, _ = build([(6, 3), (8, 1)], [(0, 1, 1)])
    assert write_smiles(m) == 'CO'
    assert smw_traversal(m)['unknown_h'] == ()
    assert m.unknown_h_count == 0


# ------------------------------------------------------------------------------------------------
# STEREO.  Where the sentinel is not just unspellable but moves a written position.
def test_a_centre_whose_hydrogen_position_is_unknown_loses_its_sign():
    """No `@`, and the atom named in `lost`.

    The sign would be a claim about four positions in the written order, and the arena will not say
    whether one of them exists.  Writing one anyway would be the exact failure ruling F26 exists for,
    reached by a different road: a sign that is right for some inputs and silently wrong for others.

    WHO REFUSES, measured rather than assumed: perception does, before the writer is asked --
    `stereo_units()` is empty here although `parity_of` still answers 2.  So the writer's own guard in
    `smw_sign_of` is not what produces this, and the report is: a STATED parity with no unit under it
    is a loss, whatever declined to name the unit.  That rule is why the string does not go out silent
    on the arena's behalf.
    """
    m, sids = build(*BROMOCHLOROFLUOROMETHANE)
    m.set_parity(sids[0], 2)
    assert m.parity_of(sids[0]) == 2 and m.stereo_units() == []
    s = write_smiles(m)
    assert sign_in(s) == 0
    probe = smw_traversal(m)
    assert probe['lost'] == (sids[0],)
    assert probe['unknown_h'] == (sids[0],)


def test_the_same_centre_keeps_its_sign_once_the_count_is_stated():
    """The control for the test above, and the proof that the refusal is the sentinel's doing and not
    the fixture's: one field changes and the `@` comes back."""
    m, sids = build([(6, 1), (9, 0), (17, 0), (35, 0)], [(0, 1, 1), (0, 2, 1), (0, 3, 1)])
    m.set_parity(sids[0], 2)
    assert sign_in(write_smiles(m)) != 0
    probe = smw_traversal(m)
    assert probe['lost'] == ()
    assert probe['unknown_h'] == ()


def test_a_fully_substituted_centre_KEEPS_its_sign_now_that_perception_narrowed():
    """FOUR heavy neighbours leave no room for a hydrogen, so the count could not have moved a
    position -- and as of arena 6835438 perception agrees and emits the unit.

    This assertion is the reversal of what this file measured on 2026-09-02, and the reversal is the
    point.  It read `== 0` and `lost == (sid,)` then, recorded as a MEASUREMENT of somebody else's
    refusal rather than as this file's rule, with the prediction written into the docstring: "if
    perception ever narrows its refusal, this test changes to `!= 0` and `lost == ()` and nothing in
    the writer has to move".  Perception narrowed, the test changed, and NOTHING IN THE WRITER MOVED --
    which is what the writer's own guard (`smw_h_frame_unknown`) was kept unreachable for.  The atom is
    still in `unknown_h`: the count is still unstated, it just never mattered to a position.
    """
    m, sids = build(*HALOMETHANE)
    m.set_parity(sids[0], 2)
    assert len(m.stereo_units()) == 1
    assert sign_in(write_smiles(m)) != 0
    probe = smw_traversal(m)
    assert probe['lost'] == ()
    assert probe['unknown_h'] == (sids[0],)


def test_an_axial_sign_is_refused_when_a_TERMINAL_says_nothing():
    """The allene case: the four directions belong to the TERMINALS, so an unknown count two bonds
    from the anchor is what removes the sign.

    The anchor itself has a perfectly known count here (a chain carbon has no hydrogens), which is
    why reading the sentinel at the signed atom alone would miss this entirely.
    """
    m, sids = build(*DIFLUOROALLENE)
    m.set_parity(sids[1], 2)
    assert sign_in(write_smiles(m)) == 0
    probe = smw_traversal(m)
    assert probe['lost'] == (sids[1],)
    assert probe['unknown_h'] == (sids[2],)   # the loss is at the centre, the sentinel at a terminal


def test_the_axial_control_writes_a_sign():
    m, sids = build(*DIFLUOROALLENE_KNOWN)
    m.set_parity(sids[1], 2)
    assert sign_in(write_smiles(m)) != 0
    assert smw_traversal(m)['lost'] == ()


def test_a_tetrasubstituted_axis_KEEPS_its_sign_now_too():
    """Both terminals full, so no terminal needed a hydrogen position, and the axis survives -- the
    same reversal as `halomethane`, one bond further out.

    Worth keeping as its own row rather than folding into that one: the axial refusal lives in
    `smw_allene_order`, a different piece of the writer from the tetrahedral path, and the two were
    reported and narrowed separately.  Which atom carries the sentinel is asserted because it is NOT
    the signed atom -- the loss would have been at the centre, the unstated count is at a terminal."""
    m, sids = build(*TETRAHALOALLENE)
    m.set_parity(sids[1], 2)
    assert len(m.stereo_units()) == 1
    assert sign_in(write_smiles(m)) != 0
    probe = smw_traversal(m)
    assert probe['lost'] == ()
    assert probe['unknown_h'] == (sids[2],)


def test_the_refusal_is_reported_under_suppressed_bond_tokens_too():
    """`!b` drops every bond token, and a tetrahedral sign normally survives that.  This one does
    not, so it has to be reported on that path as well -- the `!b` branch of `smw_directions` is a
    separate piece of code and would have been a separate hole."""
    m, sids = build(*BROMOCHLOROFLUOROMETHANE)
    m.set_parity(sids[0], 2)
    assert sign_in(write_smiles(m, '!b')) == 0
    assert smw_traversal(m, '!b')['lost'] == (sids[0],)


# ------------------------------------------------------------------------------------------------
# INVARIANCE.  A refusal that reads uninitialised scratch is a class of bug one measurement cannot
# see, so every creation order is swept.
def sweep(fixture, anchor, frame, want, spec=''):
    """Every creation order of `fixture`, one CONFIGURATION, the set of strings."""
    atoms, bonds = fixture
    strings = set()
    for order in permutations(range(len(atoms))):
        m, sids = build(atoms, bonds, order)
        if frame is None:
            m.set_parity(sids[anchor], 2)
        else:
            configure(m, sids[anchor], tuple(sids[f] for f in frame), want)
        strings.add(write_smiles(m, spec))
    return strings


# (name, fixture, anchor, frame, want).  `frame is None` is the refused-sign case, where the stored
# value cannot reach the string and so needs no configuring -- and a sweep that went from one string
# to two THERE would mean the refusal itself had stopped being invariant.
SWEEPS = [('bromochlorofluoromethane', BROMOCHLOROFLUOROMETHANE, 0, None, 0),
          ('halomethane', HALOMETHANE, 0, (1, 2, 3, 4), 1)]


@mark.parametrize('name,fixture,anchor,frame,want', SWEEPS)
def test_the_string_is_the_same_from_every_creation_order(name, fixture, anchor, frame, want):
    """Full factorial over the creation orders, with a CONFIGURATION stated.

    The point is not only canonical invariance, which the other files sweep too: before the guard
    landed, the unmatched frame walked off the end of the initialised part of a four-entry array, and
    the parity it computed was whatever the stack held.  A single-order assertion cannot see that;
    24 and 120 orders that all agree can.
    """
    assert len(sweep(fixture, anchor, frame, want)) == 1


@mark.parametrize('name,fixture,anchor,frame,want,n', [SWEEPS[0] + (12,), SWEEPS[1] + (48,)])
def test_the_sweep_can_fail(name, fixture, anchor, frame, want, n):
    """Ruling F102: the sweeps above are worth nothing unless they are shown able to fail.

    Stored-slot order (`i`) is the writer with its canonicalisation switched off, and it produces a
    different string for most creation orders of the same molecule.  Both numbers are recorded, per
    fixture, so that a change which quietly collapsed either would be visible -- and the halomethane
    row is here because it is the one whose sign now reaches the string, which makes it the row where
    a broken configuration would hide.
    """
    assert len(sweep(fixture, anchor, frame, want, 'i')) == n
