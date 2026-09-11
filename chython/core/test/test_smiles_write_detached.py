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
"""DETACHED SMILES: a fragment whose cut bonds are RING BONDS, so fragments re-join by concatenation.

The one idea being tested is that `%12` is not a new syntax.  It is the notation's own way of saying
"this bond's other end is somewhere else in the text", so `CC%12` and `N%12` become one molecule by
being written next to each other, and no reader needs to be told anything.  Everything in this file
follows from that: the attachment ids come out of the ring-closure pool (so an internal closure may
never take one), an atom's token is computed over the WHOLE molecule (because the fragment is only ever
read after a join), and a tetrahedral sign survives while a cis/trans one cannot.

WHY THE SURFACE IS WIDER THAN A TEXT JOIN.  chython 2's `sticky_smiles` is the name a porting reader
looks for; it joins TEXT, so it needs the two attachment atoms to land at the two ENDS of the string
and carries exactly two of them, in one component.  Joining BONDS puts no condition on where the
attachment atom sits -- which is why the tests below attach to an atom in the middle of a ring, to an
aromatic atom, to a double bond, and to a salt.

Public compounds throughout, and the two graph-theoretic fixtures (`k2n_methyl`) are there because a
molecule needs eleven simultaneously open rings before the attachment ids and the closure numbers can
collide at all, and no small drug-like molecule has eleven.
"""
from itertools import permutations

from pytest import mark, raises

from chython.core import MoleculeContainer
from chython.core._core import DetachedSmiles, detached_smiles, smw_traversal, write_smiles


# ------------------------------------------------------------------------------------------------
# FIXTURE PLUMBING, the same shape as the other writer test files'.
def build(atoms, bonds, order=None):
    m = MoleculeContainer()
    sids = {}
    for j in (range(len(atoms)) if order is None else order):
        element, hydrogens = atoms[j]
        sids[j] = m.add_atom(element, implicit_h=hydrogens)
    for a, b, o in bonds:
        m.add_bond(sids[a], sids[b], o)
    return m, sids


def sign_in(smiles):
    if '@@' in smiles:
        return 2
    if '@' in smiles:
        return 1
    return 0


def configure(m, sid, frame, want):
    """Store the parity that makes `translate_stereo(sid, frame)` answer `want`.

    A STORED PARITY IS NOT A CONFIGURATION -- it is a configuration relative to the atom's refs order,
    which the CREATION order decides.  So a sweep that called `set_parity(sid, 2)` in every creation
    order would be sweeping different molecules and would rightly produce two strings.  Asking the
    production `translate_stereo` which of the two values lands the wanted arrangement of atom
    identities is the same technique `test_smiles_write_stereo.py` uses, and for the same reason:
    computing it here would reimplement the arithmetic under test.
    """
    for parity in (1, 2):
        m.set_parity(sid, parity)
        if m.translate_stereo(sid, frame) == want:
            return parity
    raise AssertionError('neither parity gives %r in frame %r' % (want, frame))


def _frame(sids, frame):
    return tuple(None if f is None else sids[f] for f in frame)


def k2n_methyl(k):
    """Two atoms bridged by `k` others, plus a methyl: `k - 1` rings, all open at once, one cut point.

    The rings are what force the closure numbers past 9 -- every bridge after the first closes at the
    same atom -- and the methyl is the only acyclic bond in it, so it is the only thing that can be
    cut.  Not a compound anyone has in a bottle; it is here for the id space and nothing else.
    """
    m = MoleculeContainer()
    a = m.add_atom(6, implicit_h=0)
    b = m.add_atom(6, implicit_h=0)
    for _ in range(k):
        x = m.add_atom(6, implicit_h=0)
        m.add_bond(a, x, 1)
        m.add_bond(b, x, 1)
    me = m.add_atom(6, implicit_h=3)
    m.add_bond(a, me, 1)
    return m, a, me


# Ethanol.  The smallest molecule with a cuttable bond, and the one the atom-token rule is visible on.
ETHANOL = ([(6, 3), (6, 2), (8, 1)], [(0, 1, 1), (1, 2, 1)])
# Sodium acetate: a real SALT, so the dropped side's growth rule has an unrelated component to leave
# alone.  Growing the retained set from the keeps instead would take the sodium with the acetate half.
SODIUM_ACETATE = ([(6, 3), (6, 0), (8, 0), (8, 0), (11, 0)], [(0, 1, 1), (1, 2, 2), (1, 3, 1)])
# Benzene with STORED aromatic bonds (order 4), for the ring refusal.
BENZENE = ([(6, 1)] * 6, [(0, 1, 4), (1, 2, 4), (2, 3, 4), (3, 4, 4), (4, 5, 4), (5, 0, 4)])
# Toluene, stored aromatic: the attachment sits on an AROMATIC atom.
TOLUENE = ([(6, 0)] + [(6, 1)] * 5 + [(6, 3)],
           [(0, 1, 4), (1, 2, 4), (2, 3, 4), (3, 4, 4), (4, 5, 4), (5, 0, 4), (0, 6, 1)])
# 1-fluoroethan-1-amine: one tetrahedral centre with three heavy neighbours and an implicit hydrogen,
# so every cut moves a written position and the sign has to move with it.
FLUOROETHYLAMINE = ([(6, 3), (6, 1), (7, 2), (9, 0)], [(1, 0, 1), (1, 2, 1), (1, 3, 1)])
# (Z)-1,2-difluoroethene: the cis/trans unit, whose configuration is two tokens on two bonds.
DIFLUOROETHENE = ([(9, 0), (6, 1), (6, 1), (9, 0)], [(0, 1, 1), (1, 2, 2), (2, 3, 1)])
# 1,3-difluoroallene: the axial unit, whose four directions belong to the TERMINALS.
DIFLUOROALLENE = ([(6, 1), (6, 0), (6, 1), (9, 0), (9, 0)],
                  [(0, 1, 2), (1, 2, 2), (0, 3, 1), (2, 4, 1)])
# Propene, for a cut across a DOUBLE bond.
PROPENE = ([(6, 2), (6, 1), (6, 3)], [(0, 1, 2), (1, 2, 1)])


# ------------------------------------------------------------------------------------------------
# THE CUT.
def test_a_cut_gives_two_halves_that_partition_the_molecule():
    """Ethanol as a methyl and a hydroxymethyl, each naming the same attachment.

    The partition is the assertion, not the strings: every atom is in exactly one fragment, which is
    what makes a join total.  A rule that decided the dropped side by reachability from the KEEP would
    also satisfy this on a connected molecule, which is why the salt test below exists.
    """
    m, sids = build(*ETHANOL)
    keep_methyl = detached_smiles(m, {10: (sids[0], sids[1])})
    keep_rest = detached_smiles(m, {10: (sids[1], sids[0])})
    assert keep_methyl.text == 'C%10'
    assert keep_rest.text == 'C%10O'
    assert set(keep_methyl.order) | set(keep_rest.order) == set(sids.values())
    assert not set(keep_methyl.order) & set(keep_rest.order)
    assert keep_methyl.open_ids == keep_rest.open_ids == (10,)


def test_the_cut_is_an_ORDERED_pair():
    """`(keep, drop)` and not a bond, because there is nothing in `C-C` to say which half is wanted.

    The two orderings give two different fragments of two different sizes, so the order is load-bearing
    rather than a convention that could have been either way.
    """
    m, sids = build(*ETHANOL)
    assert detached_smiles(m, {10: (sids[0], sids[1])}).atom_count == 1
    assert detached_smiles(m, {10: (sids[1], sids[0])}).atom_count == 2


def test_the_dropped_side_grows_from_the_DROP_so_a_salt_keeps_its_counter_ion():
    """Sodium acetate cut at the C-C: the sodium is in neither the retained nor the dropped path, and
    it stays.

    THE MEASURED CONSEQUENCE OF A DESIGN CHOICE.  `dropped` is the closure of the named drop atoms
    under the uncut bonds; the retained set is everything else.  The other formulation -- retained is
    what the keeps reach -- looks equivalent and is not: an unrelated component reaches no keep, so it
    would vanish from the fragment without anybody naming it, and a caller cutting an ester in a
    hydrochloride salt would silently lose the HCl.
    """
    m, sids = build(*SODIUM_ACETATE)
    m.set_charge(sids[3], -1)
    m.set_charge(sids[4], 1)
    f = detached_smiles(m, {10: (sids[0], sids[1])})
    assert f.text == 'C%10.[Na+]'
    assert set(f.order) == {sids[0], sids[4]}


def test_no_cuts_is_the_whole_molecule():
    """The degenerate case is not special-cased, and a fragment with nothing open is a molecule.

    It is the base of a join: a complete component can be joined to a fragment as itself.
    """
    m, _ = build(*ETHANOL)
    f = detached_smiles(m, {})
    assert f.text == write_smiles(m)
    assert f.open_ids == ()
    assert str(f) == write_smiles(m)


def test_an_empty_molecule_is_an_empty_fragment():
    f = detached_smiles(MoleculeContainer(), {})
    assert f.text == '' and f.order == () and f.open_ids == ()


# ------------------------------------------------------------------------------------------------
# THE REFUSALS.  Every one of them names atoms, because the caller's next move is to edit the cut list
# and a message they cannot act on is a message that sends them to read this source.
def test_a_ring_bond_is_refused_and_the_message_names_the_path_round():
    """One id cannot carry both ends of a ring opening.

    Not a limitation of the notation -- `C1CCCCC1` opens a ring with one number -- but of a CUT: the
    two ends would both be attachments, and two occurrences of `%10` in one fragment is a ring closure
    inside it, which re-forms the bond the caller asked to break.  Silently.  So it is refused, and the
    path is in the message because the fix is a second cut somewhere along it.
    """
    m, sids = build(*BENZENE)
    with raises(ValueError) as e:
        detached_smiles(m, {10: (sids[0], sids[1])})
    assert 'ring' in str(e.value)
    for sid in sids.values():
        assert str(sid) in str(e.value)          # the whole cycle, so a second cut can be chosen


def test_two_cuts_open_a_ring():
    """The corollary, and the reason the refusal above is not a dead end: cut TWO bonds and the ring
    opens, with the retained atom carrying two attachments.

    Benzene cut at both bonds of one carbon is that carbon plus a five-atom chain -- and the five-atom
    chain is the fragment with two open ids, which is exactly how a linker is written.
    """
    m, sids = build(*BENZENE)
    f = detached_smiles(m, {10: (sids[1], sids[0]), 11: (sids[5], sids[0])})
    assert f.open_ids == (10, 11)
    assert f.atom_count == 5
    assert '%10' in f.text and '%11' in f.text
    other = detached_smiles(m, {10: (sids[0], sids[1]), 11: (sids[0], sids[5])})
    assert other.atom_count == 1
    assert other.text.count('%1') == 2           # both ids on the one retained atom


def test_atoms_that_are_not_bonded_are_refused():
    m, sids = build(*ETHANOL)
    with raises(ValueError, match='not bonded'):
        detached_smiles(m, {10: (sids[0], sids[2])})


def test_one_atom_cannot_be_both_sides():
    m, sids = build(*ETHANOL)
    with raises(ValueError, match='both'):
        detached_smiles(m, {10: (sids[0], sids[0])})


def test_a_bond_named_by_two_cuts_is_refused():
    """Two ids on one bond would write `%10%11` at the retained atom and leave both dangling."""
    m, sids = build(*ETHANOL)
    with raises(ValueError, match='two cuts'):
        detached_smiles(m, {10: (sids[0], sids[1]), 11: (sids[0], sids[1])})


def test_an_atom_cannot_be_kept_by_one_cut_and_dropped_by_another():
    """A contradiction the caller has to resolve, and it is checked BEFORE the walk so that the walk's
    answer cannot be blamed for it."""
    m, sids = build(*ETHANOL)
    with raises(ValueError, match='dropped side'):
        detached_smiles(m, {10: (sids[0], sids[1]), 11: (sids[1], sids[2])})


@mark.parametrize('bad', [0, 1, 9, 100, 255])
def test_an_id_outside_ten_to_ninety_nine_is_refused(bad):
    """The floor is a MEASUREMENT, not a preference: `%05` is rejected by RDKit 2026.03.4 and by
    chython 2, so a fixed-width low spelling is not available, and a BARE digit is indistinguishable
    from an ordinary ring closure -- which is the one thing an attachment must never be mistaken for.
    """
    m, sids = build(*ETHANOL)
    with raises(ValueError, match='outside'):
        detached_smiles(m, {bad: (sids[0], sids[1])})


def test_an_unknown_atom_is_refused():
    m, sids = build(*ETHANOL)
    with raises(KeyError):
        detached_smiles(m, {10: (sids[0], 999)})


def test_a_cut_must_be_a_pair():
    m, sids = build(*ETHANOL)
    with raises(TypeError):
        detached_smiles(m, {10: (sids[0], sids[1], sids[2])})
    with raises(TypeError):
        detached_smiles(m, [(10, sids[0], sids[1])])
    with raises(TypeError):
        detached_smiles(m, {'10': (sids[0], sids[1])})


# ------------------------------------------------------------------------------------------------
# THE ATTACHMENT IS A RING BOND, and shares the ring bonds' number space.
def test_the_attachment_is_always_percent_two_digits():
    m, sids = build(*ETHANOL)
    for i in (10, 11, 42, 99):
        text = detached_smiles(m, {i: (sids[0], sids[1])}).text
        assert text == 'C%%%d' % i
        assert '%0' not in text


def test_the_attachment_is_written_in_the_ring_bond_position():
    """`atom ringbond* branch*` is the grammar, and an attachment is one of the ringbonds.

    Visible in the string: the `%10` comes before the branch, not after it and not inside it.  This is
    the same fact `smw_direction_order` relies on when it says a tetrahedral sign survives -- the
    attachment holds a POSITION in the neighbour order, and it holds the one a closure would.
    """
    m, sids = build(*FLUOROETHYLAMINE)
    f = detached_smiles(m, {10: (sids[1], sids[2])})
    assert f.text == 'C%10(C)F'
    assert f.text.index('%10') < f.text.index('(')


def test_an_internal_closure_never_takes_an_attachment_number():
    """Eleven open rings and a cut using id 10: the closures step over 10 and take 12 instead.

    `%12` and `12` are ONE ring bond to every reader, so a closure reusing an attachment's number
    would be paired with the attachment on a join -- silently, giving a valid molecule that is not the
    one asked for.  Withheld for the whole string and never released, because a fragment's numbers must
    mean the same thing at every position in it: the join happens at the string level and knows nothing
    about where a number was in scope.
    """
    m, a, me = k2n_methyl(12)
    f = detached_smiles(m, {10: (a, me)})
    assert f.open_ids == (10,)
    assert 10 not in f.closure_ids
    assert max(f.closure_ids) > 10               # it needed eleven numbers and found an eleventh
    assert '%10' in f.text


def test_reserve_withholds_another_fragments_ids():
    """`reserve` is how the fragments of one join agree: every attachment id in the whole set is
    withheld from every fragment's own closures.

    Needed because a number in scope anywhere in the joined text is in scope everywhere in it -- the
    reader has no notion of "this part of the string".
    """
    m, a, me = k2n_methyl(11)
    plain = detached_smiles(m, {11: (a, me)})
    reserved = detached_smiles(m, {11: (a, me)}, reserve=(10, 12))
    assert 10 in plain.closure_ids
    assert 10 not in reserved.closure_ids and 12 not in reserved.closure_ids
    assert len(reserved.closure_ids) == len(plain.closure_ids)


def test_both_fragments_write_the_bond_token():
    """A cut across a double bond gives `C=%10` on BOTH sides.

    Measured rather than reasoned: RDKit, Indigo, OpenBabel and chython 2 all accept a ring bond whose
    order is stated at both ends and all four reject a clash, so writing the token at one end only
    would make the join order-dependent and buy nothing.
    """
    m, sids = build(*PROPENE)
    left = detached_smiles(m, {10: (sids[0], sids[1])})
    right = detached_smiles(m, {10: (sids[1], sids[0])})
    assert left.text == 'C=%10'
    assert right.text.count('=%10') == 1
    assert str(DetachedSmiles.join(left, right)) == 'C=%10.C=%10C'


def test_the_attachments_report_names_the_bond():
    """`smw_traversal` exposes the cut as the traversal sees it, so a test can check the boundary
    without reading it back out of the string."""
    m, sids = build(*ETHANOL)
    probe = smw_traversal(m, cuts={10: (sids[0], sids[1])})
    assert probe['attachments'] == ((sids[0], sids[1], 10),)
    assert probe['order'] == (sids[0],)
    assert smw_traversal(m)['attachments'] == ()


# ------------------------------------------------------------------------------------------------
# THE ATOM TOKEN COMES FROM THE WHOLE MOLECULE, because a fragment is only ever read after a join.
def test_the_atom_token_is_written_from_the_whole_molecule():
    """`C%10` and not `[CH3]%10`, and the difference is the whole rule.

    The methyl carbon has three stored hydrogens and, in the MOLECULE, one heavy neighbour -- which is
    what the valence model derives 3 from, so the bare spelling is faithful and no bracket is needed.
    Count the neighbours in the FRAGMENT instead and there are none, the model derives 4, the stored 3
    disagrees, and the writer brackets the atom to state it.  Both spellings describe the same fragment
    correctly; only one describes the JOINED molecule correctly, and the joined molecule is the only
    thing a fragment is ever read as.
    """
    m, sids = build(*ETHANOL)
    assert detached_smiles(m, {10: (sids[0], sids[1])}).text == 'C%10'
    # And when something else forces the bracket, the count inside it is still the molecule's.
    assert detached_smiles(m, {10: (sids[0], sids[1])}, 'h').text == '[CH3]%10'


def test_an_aromatic_atom_can_carry_an_attachment():
    """Toluene's ring, cut from its methyl: the ring stays aromatic and the attachment sits on an
    aromatic carbon.

    An atom inside a ring is never at an END of the string, so a text join has nothing to cut; the
    `%10` token goes at the atom's own written position instead.
    """
    m, sids = build(*TOLUENE)
    ring = detached_smiles(m, {10: (sids[0], sids[6])})
    assert ring.text == 'c1c%10cccc1'
    assert str(DetachedSmiles.join(ring, detached_smiles(m, {10: (sids[6], sids[0])}))) == \
        'c1c%10cccc1.C%10'


# ------------------------------------------------------------------------------------------------
# STEREO.  A tetrahedral sign survives a cut; a cis/trans configuration cannot.
@mark.parametrize('drop', [0, 2, 3])
def test_a_tetrahedral_sign_survives_a_cut_and_is_the_written_orders(drop):
    """The anti-drift test, with a cut: the sign in the fragment is `translate_stereo` of the order the
    fragment writes -- and the attachment is one of the four positions in that order.

    Two independent implementations of the same arithmetic on one agreed input, as in
    `test_smiles_write_stereo.py`: the C path through the writer, and `translate_stereo` through the
    container.  What the cut adds is that the frame CHANGES -- the dropped neighbour moves out of the
    branch list and into the ring-bond group -- so a writer that copied the whole molecule's sign
    across would fail here for two of the three cuts.
    """
    m, sids = build(*FLUOROETHYLAMINE)
    for parity in (1, 2):
        m.set_parity(sids[1], parity)
        cuts = {10: (sids[1], sids[drop])}
        written = smw_traversal(m, cuts=cuts)['directions'][sids[1]]
        expected = m.translate_stereo(sids[1], written)
        assert expected in (1, 2)
        assert sids[drop] in written              # the dropped atom still holds a position
        # parity 2 (odd) is `@`: the same anchor `test_smiles_write_stereo.py` states.
        assert sign_in(detached_smiles(m, cuts).text) == (1 if expected == 2 else 2)


def test_the_sign_is_not_the_whole_molecules_sign():
    """The control for the test above, and the shape of the bug it would catch: cutting the amine off
    1-fluoroethanamine flips the character, because the nitrogen moves from the third written position
    to the first ring-bond one and one exchange is odd."""
    m, sids = build(*FLUOROETHYLAMINE)
    m.set_parity(sids[1], 2)
    assert sign_in(write_smiles(m)) != sign_in(detached_smiles(m, {10: (sids[1], sids[2])}).text)


def test_a_cis_trans_unit_is_refused_when_a_reference_is_dropped():
    """(Z)-1,2-difluoroethene cut at a C-F: no `/` anywhere, and the anchor is in `lost`.

    A cis/trans configuration is TWO tokens on TWO bonds and its content is their relation, so one of
    them landing in another string leaves nothing readable.  Unlike a tetrahedral sign, which is four
    positions at ONE atom and therefore entirely inside whichever fragment holds that atom, this cannot
    be rescued by giving the attachment a position -- there is no direction to put on a bond the string
    does not contain.
    """
    m, sids = build(*DIFLUOROETHENE)
    m.set_parity(sids[1], 2)
    assert '/' in write_smiles(m)
    f = detached_smiles(m, {10: (sids[2], sids[3])})
    assert '/' not in f.text and '\\' not in f.text
    assert f.lost == (sids[1],)


def test_a_cis_trans_unit_is_refused_when_the_partner_terminal_is_dropped():
    """The other way a cut splits the unit: the double bond itself is cut."""
    m, sids = build(*DIFLUOROETHENE)
    m.set_parity(sids[1], 2)
    f = detached_smiles(m, {10: (sids[1], sids[2])})
    assert '/' not in f.text and '\\' not in f.text
    assert f.lost == (sids[1],)


def test_a_cis_trans_unit_the_cut_does_not_touch_survives():
    """The control.  (Z)-hex-3-ene cut at a terminal methyl keeps the configuration, so the refusals
    above are about the CUT and not about the presence of a cut."""
    atoms = [(6, 3), (6, 2), (6, 1), (6, 1), (6, 2), (6, 3)]
    bonds = [(0, 1, 1), (1, 2, 1), (2, 3, 2), (3, 4, 1), (4, 5, 1)]
    m, sids = build(atoms, bonds)
    m.set_parity(sids[2], 2)
    f = detached_smiles(m, {10: (sids[1], sids[0])})
    assert '/' in f.text or '\\' in f.text
    assert f.lost == ()


def test_an_axial_sign_survives_a_dropped_substituent():
    """1,3-difluoroallene cut at a terminal C-F: the sign stays, because the attachment holds the
    position the fluorine held.

    The axial frame is the two TERMINALS' directions, and an attachment is a direction like any other.
    """
    m, sids = build(*DIFLUOROALLENE)
    m.set_parity(sids[1], 2)
    f = detached_smiles(m, {10: (sids[2], sids[4])})
    assert sign_in(f.text) != 0
    assert f.lost == ()


def test_an_axial_sign_is_refused_when_a_TERMINAL_is_dropped():
    """Cutting the axis itself: the sign is a claim about the relation between the two ends and one end
    is gone."""
    m, sids = build(*DIFLUOROALLENE)
    m.set_parity(sids[1], 2)
    f = detached_smiles(m, {10: (sids[1], sids[2])})
    assert sign_in(f.text) == 0
    assert f.lost == (sids[1],)


def test_a_dropped_atoms_own_configuration_is_not_reported_as_lost():
    """`lost` is about this string's coverage of what it claims to describe, and it does not claim to
    describe the dropped part.

    So the amine-side fragment of a chiral molecule reports nothing, although a centre's configuration
    is certainly absent from it -- the centre is absent from it.
    """
    m, sids = build(*FLUOROETHYLAMINE)
    m.set_parity(sids[1], 2)
    f = detached_smiles(m, {10: (sids[2], sids[1])})
    assert f.order == (sids[2],)
    assert f.lost == ()


# ------------------------------------------------------------------------------------------------
# JOIN.
def test_two_halves_join_back():
    """Concatenation with a `.`, and the ids are closed.

    `.` and not nothing: each fragment stays its own component of the TEXT, and the ring bond is what
    makes it one molecule.  That is also what keeps every fragment's atom tokens valid -- a fragment's
    first atom is still a component leader after the join, so nothing it wrote is read differently.
    """
    m, sids = build(*ETHANOL)
    j = DetachedSmiles.join(detached_smiles(m, {10: (sids[0], sids[1])}),
                            detached_smiles(m, {10: (sids[1], sids[0])}))
    assert str(j) == 'C%10.C%10O'
    assert j.open_ids == ()
    assert j.closure_ids == (10,)                # the join CLOSED it; it is a ring bond now
    assert j.order == (sids[0], sids[1], sids[2])
    assert j.atom_count == 3


def test_a_partial_join_leaves_an_id_open():
    """An id in one fragment stays open, so a molecule can be assembled in stages and the result is
    still a `DetachedSmiles` rather than a string."""
    m, sids = build(*BENZENE)
    middle = detached_smiles(m, {10: (sids[1], sids[0]), 11: (sids[5], sids[0])})
    m2, sids2 = build(*ETHANOL)
    cap = detached_smiles(m2, {10: (sids2[0], sids2[1])})
    j = DetachedSmiles.join(middle, cap)
    assert j.open_ids == (11,)
    assert 10 in j.closure_ids
    assert '%11' in j.text


def test_an_id_in_three_fragments_is_refused():
    m, sids = build(*ETHANOL)
    f = detached_smiles(m, {10: (sids[0], sids[1])})
    with raises(ValueError, match='more than two'):
        DetachedSmiles.join(f, f, f)


def test_a_closure_that_collides_with_an_attachment_is_refused_and_names_reserve():
    """The collision `reserve` exists for, and the message says so.

    Refused rather than renumbered: the texts are already written, and renumbering one would mean
    editing it -- which is the text surgery this whole design exists to avoid.
    """
    ringy, a, me = k2n_methyl(12)
    big = detached_smiles(ringy, {11: (a, me)})
    assert 10 in big.closure_ids
    m, sids = build(*ETHANOL)
    small = detached_smiles(m, {10: (sids[0], sids[1])})
    with raises(ValueError, match='reserve'):
        DetachedSmiles.join(big, small)
    fixed = detached_smiles(ringy, {11: (a, me)}, reserve=(10,))
    j = DetachedSmiles.join(fixed, small)         # and with the reservation it goes through
    assert j.open_ids == (10, 11)


def test_join_needs_detached_smiles():
    m, sids = build(*ETHANOL)
    with raises(ValueError):
        DetachedSmiles.join()
    with raises(TypeError):
        DetachedSmiles.join(detached_smiles(m, {}), 'C%10')


def test_a_fragment_with_an_open_id_is_not_a_molecule():
    """Documented and asserted: `str()` of an open fragment is a dangling ring bond, which every reader
    rejects.  That is the point -- a fragment cannot be mistaken for a molecule, so nothing downstream
    can accidentally treat one as a SMILES."""
    m, sids = build(*ETHANOL)
    f = detached_smiles(m, {10: (sids[0], sids[1])})
    assert '%10' in str(f) and f.open_ids == (10,)
    assert str(DetachedSmiles.join(f, detached_smiles(m, {10: (sids[1], sids[0])}))).count('%10') == 2


# ------------------------------------------------------------------------------------------------
# THE CXSMILES TAIL, which is the one part of a fragment that cannot simply be concatenated.
def test_the_tail_is_reindexed_on_a_join():
    """A radical index counts atoms from the START of the whole string, so the second fragment's
    indices shift by the first's atom count.

    Kept as STRUCTURE and formatted once, rather than edited as text: an index list is easy to shift
    and hard to shift correctly by regular expression, and the writer would then be parsing the format
    it also writes.
    """
    m, sids = build([(6, 3), (6, 2), (8, 0)], [(0, 1, 1), (1, 2, 1)])
    m.set_radical(sids[2], True)
    left = detached_smiles(m, {10: (sids[0], sids[1])})
    right = detached_smiles(m, {10: (sids[1], sids[0])})
    assert str(left) == 'C%10'                   # no radical in this half, so no tail at all
    assert str(right) == 'C%10[O] |^1:1|'
    assert str(DetachedSmiles.join(left, right)) == 'C%10.C%10[O] |^1:2|'


def test_two_fragments_and_groups_are_renumbered():
    """Two fragments' `&1` are two DIFFERENT groups, and keeping the number would claim their atoms
    invert together."""
    m, sids = build(*FLUOROETHYLAMINE)
    m.set_parity(sids[1], 2)
    m.set_stereo_group(sids[1], 3, 1)
    f = detached_smiles(m, {10: (sids[1], sids[2])})
    g = detached_smiles(m, {11: (sids[1], sids[2])})
    assert str(f).endswith('|&1:0|')
    assert str(DetachedSmiles.join(f, g)).endswith('|&1:0,&2:3|')


def test_the_labels_are_concatenated_on_a_join_and_not_shifted():
    """`$...$` is POSITIONAL -- one entry per atom -- so a join concatenates the lists, and a fragment
    carrying no label still owes the field a blank for each of its own atoms."""
    m, sids = build([(6, 3), (6, 2), (8, 1)], [(0, 1, 1), (1, 2, 1)])
    m.set_aliases({sids[0]: b'Me', sids[2]: b'OH'})
    left = detached_smiles(m, {10: (sids[0], sids[1])})
    right = detached_smiles(m, {10: (sids[1], sids[0])})
    assert str(left) == 'C%10 |$Me$|'
    assert str(right) == 'C%10O |$;OH$|'
    assert str(DetachedSmiles.join(left, right)) == 'C%10.C%10O |$Me;;OH$|'


def test_no_tail_is_written_under_the_no_cxsmiles_key():
    m, sids = build([(6, 3), (6, 2), (8, 0)], [(0, 1, 1), (1, 2, 1)])
    m.set_radical(sids[2], True)
    f = detached_smiles(m, {10: (sids[1], sids[0])}, '!x')
    assert str(f) == 'C%10[O]'
    assert f.tail == ([], [], {}, {}, [])


# ------------------------------------------------------------------------------------------------
# INVARIANCE.  A fragment is written by the same canonical machinery as a molecule, so the claim is the
# same: the same molecule and the same cut give the same text from every creation order.
@mark.parametrize('name,fixture,keep,drop', [('ethanol', ETHANOL, 1, 0),
                                             ('fluoroethylamine', FLUOROETHYLAMINE, 1, 2)])
def test_a_fragment_is_the_same_from_every_creation_order(name, fixture, keep, drop):
    atoms, bonds = fixture
    strings = set()
    for order in permutations(range(len(atoms))):
        m, sids = build(atoms, bonds, list(order))
        if name == 'fluoroethylamine':
            configure(m, sids[1], _frame(sids, (0, 2, 3, None)), 2)
        strings.add(detached_smiles(m, {10: (sids[keep], sids[drop])}).text)
    assert len(strings) == 1


def test_the_sweep_can_fail():
    """Ruling F102: the sweep above is worth nothing unless it is shown able to fail.

    Stored-slot order (`i`) is the writer with its canonicalisation switched off, and the number of
    distinct strings is recorded so that a change which quietly collapsed it would be visible.
    """
    atoms, bonds = FLUOROETHYLAMINE
    strings = set()
    for order in permutations(range(len(atoms))):
        m, sids = build(atoms, bonds, list(order))
        configure(m, sids[1], _frame(sids, (0, 2, 3, None)), 2)
        strings.add(detached_smiles(m, {10: (sids[1], sids[2])}, 'i').text)
    assert len(strings) == 4


# ------------------------------------------------------------------------------------------------
# NOT A SPEC KEY, and therefore not a cache key.
def test_there_is_no_format_key_for_a_detached_string():
    """A detached fragment is a FUNCTION CALL, deliberately.

    `format(mol, ...)` keys are cacheable identities of a molecule; a fragment is an identity of a
    molecule AND a cut list, and its text is not canonical for anything smaller than the pair.  Giving
    it a letter would put it in reach of any code that caches on a spec string, and the first such
    cache would be wrong for every cut but one.
    """
    m, _ = build(*ETHANOL)
    with raises(ValueError, match='unknown format key'):
        write_smiles(m, 'd')
