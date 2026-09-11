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
"""`@` and `@@`: the fixtures ruling F26 exists for.

Ruling F26's whole cost was that no fixture swept creation orders with a CONFIGURATION set, so a
writer emitting the stored parity byte looked correct.  Every sweep here does, and every one of them
is shown able to fail: the same molecule written in stored-slot order gives twelve, sixteen and
eleven distinct strings where the canonical writer gives one.

Three claims are checked against an oracle outside chython, because two of them are conventions and a
convention cannot be verified against the code that implements it:

* the implicit hydrogen occupies the position it is WRITTEN in -- RDKit 2026.03.4;
* every spelling the writer produces from every creation order is ONE molecule -- RDKit;
* core parity 2 (odd) in the ruling-F26 refs frame is `@` -- measured by the MDL epic through
  chython 2, which is the only external anchor that value has.  The core defines `even` and `odd` and
  nothing else, so this is a convention SHARED WITH THE READER rather than a fact about the arena.

The RDKit tests skip where RDKit is absent; the rest do not depend on it.
"""
from itertools import permutations
from math import factorial
from random import Random

from pytest import importorskip, mark

from chython.core import MoleculeContainer
from chython.core._core import smw_stereo_seed_labels, smw_traversal, write_smiles


# ------------------------------------------------------------------------------------------------
# FIXTURE PLUMBING.  Same shape as test_smiles_write.py's, plus the sid map -- a configuration has to
# be stated in terms of ATOM IDENTITIES to mean the same thing under two different creation orders,
# and the ids are what identify them.
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
    """Store the parity that makes `translate_stereo(sid, frame)` answer `want`; return it.

    THE POINT OF SEARCHING RATHER THAN COMPUTING.  A test that wants "this configuration" under many
    creation orders has to convert a frame of atom identities into a stored parity, and the stored
    parity is relative to the refs order, which the creation order decides.  Computing it here would
    mean reimplementing `translate_parity` in the test -- the same arithmetic the writer uses, so a
    sign error would cancel and the test would pass on a broken writer.  Trying both values and
    asking the PRODUCTION function which one lands where we want has no such blind spot: there are
    only two, and `translate_stereo` is not the code under test.
    """
    for parity in (1, 2):
        m.set_parity(sid, parity)
        if m.translate_stereo(sid, frame) == want:
            return parity
    raise AssertionError('neither parity gives %r in frame %r' % (want, frame))


def sign_in(smiles):
    """The one stereo sign in a single-centre string: 2 for `@@`, 1 for `@`, 0 for none."""
    if '@@' in smiles:
        return 2
    if '@' in smiles:
        return 1
    return 0


# CHFClBr -- one centre, three heavy neighbours and an IMPLICIT hydrogen, so it is the fixture the
# positional-hydrogen rule lives or dies on.
BROMOCHLOROFLUOROMETHANE = ([(6, 1), (9, 0), (17, 0), (35, 0)],
                            [(0, 1, 1), (0, 2, 1), (0, 3, 1)])
# CFClBrI -- one centre, FOUR heavy neighbours and no unnamed direction at all: the path with no
# positional rule to get wrong, which is what makes it the control for the fixture above.
HALOMETHANE = ([(6, 0), (9, 0), (17, 0), (35, 0), (53, 0)],
               [(0, 1, 1), (0, 2, 1), (0, 3, 1), (0, 4, 1)])
# Ethyl methyl sulfoxide -- one centre whose fourth direction is a LONE PAIR, which has no
# positional rule (measured below) and sits last.
SULFOXIDE = ([(16, 0), (6, 3), (8, 0), (6, 2), (6, 3)],
             [(0, 1, 1), (0, 2, 2), (0, 3, 1), (3, 4, 1)])
# trans-1,2-dichlorocyclopropane -- TWO centres, and each one's direction list contains a RING
# CLOSURE, so the closure's place in the written order is under test and not just a branch's.
DICHLOROCYCLOPROPANE = ([(6, 1), (6, 1), (6, 2), (17, 0), (17, 0)],
                        [(0, 1, 1), (1, 2, 1), (2, 0, 1), (0, 3, 1), (1, 4, 1)])
# (Z)-but-2-ene, for the one thing this file asserts is NOT written yet.
BUTENE = ([(6, 3), (6, 1), (6, 1), (6, 3)], [(0, 1, 1), (1, 2, 2), (2, 3, 1)])


SINGLE_CENTRE = [('bromochlorofluoromethane', BROMOCHLOROFLUOROMETHANE, 0, (1, 2, 3, None)),
                 ('halomethane', HALOMETHANE, 0, (1, 2, 3, 4)),
                 ('sulfoxide', SULFOXIDE, 0, (1, 2, 3, None))]


def _frame(sids, frame):
    return tuple(None if f is None else sids[f] for f in frame)


# ------------------------------------------------------------------------------------------------
# THE ANCHOR.
def test_the_anchor_string_itself_comes_back_out():
    """The anchor with no arithmetic in the way: the written order IS the refs order.

    The MDL epic measured the parity convention on `F[C@](Cl)(Br)I` -- chython 2's parser reads that
    string into the bool this core stores as parity 2.  Build the same molecule with fluorine at slot
    0 and the carbon at slot 1 and two things line up exactly: the refs are `(F, Cl, Br, I)`, heavy
    neighbours in ascending slot order, and stored-slot output starts at slot 0, so the string's
    written order is `(F, Cl, Br, I)` too.  The permutation between them is the identity, so parity 2
    must produce that string CHARACTER FOR CHARACTER -- no inversion count to get wrong, and nothing
    left between the measurement and the assertion.

    Every other test in this file rests on this one; `test_parity_two_in_the_refs_frame_is_the_...`
    below is the same claim with a permutation in the middle, which is what makes it a check of the
    frame arithmetic rather than of the convention.
    """
    atoms, bonds = HALOMETHANE
    m, sids = build(atoms, bonds, order=[1, 0, 2, 3, 4])
    assert [u['refs'] for u in m.stereo_units() if u['anchor'] == sids[0]] == \
        [(sids[1], sids[2], sids[3], sids[4])]
    m.set_parity(sids[0], 2)
    assert write_smiles(m, 'i') == 'F[C@](Cl)(Br)I'
    m.set_parity(sids[0], 1)
    assert write_smiles(m, 'i') == 'F[C@@](Cl)(Br)I'


def test_parity_two_in_the_refs_frame_is_the_molecule_the_anchor_predicts():
    """The one test that ties the core's `even`/`odd` to a configuration in the world.

    The chain, and every link is external to this file:

    1. the MDL epic measured that a negative signed volume is anticlockwise, is SMILES `@`, and is
       core parity 2 -- through chython 2, whose parser reads `F[C@](Cl)(Br)I` into that same bool;
    2. so CHFClBr with parity 2 and refs `(F, Cl, Br, H)` is `@` READ IN THAT ORDER;
    3. `F[C@?H](Cl)Br` writes the same four directions as `(F, H, Cl, Br)`, which is
       `(0, 3, 1, 2)` of the refs -- two inversions, EVEN -- so the sign does not change and the
       molecule is `F[C@H](Cl)Br`;
    4. the writer produces `[C@@H](F)(Cl)Br`, whose written order is `(H, F, Cl, Br)`: one more
       transposition, so the sign flips, which is why the string says `@@` and not `@`.

    RDKit closes the loop: those two strings must be one molecule.  Delete any link and the test
    fails -- an inverted anchor fails at 4, a hydrogen written last fails at 4, and a wrong frame in
    `smw_direction_order` fails at 3.
    """
    atoms, bonds = BROMOCHLOROFLUOROMETHANE
    m, sids = build(atoms, bonds)
    assert configure(m, sids[0], _frame(sids, (1, 2, 3, None)), 2) == 2
    assert write_smiles(m) == '[C@@H](F)(Cl)Br'

    chem = importorskip('rdkit.Chem')
    assert chem.MolToSmiles(chem.MolFromSmiles('[C@@H](F)(Cl)Br')) == \
        chem.MolToSmiles(chem.MolFromSmiles('F[C@H](Cl)Br'))


def test_the_two_parities_are_the_two_spellings():
    """Nothing is lost between them: parity 1 and parity 2 differ in the sign and in nothing else."""
    atoms, bonds = BROMOCHLOROFLUOROMETHANE
    m, sids = build(atoms, bonds)
    m.set_parity(sids[0], 1)
    one = write_smiles(m)
    m.set_parity(sids[0], 2)
    two = write_smiles(m)
    assert one == '[C@H](F)(Cl)Br'
    assert two == '[C@@H](F)(Cl)Br'
    assert one.replace('@H', '@@H') == two


# ------------------------------------------------------------------------------------------------
# THE POSITIONAL HYDROGEN RULE.
def test_implicit_hydrogen_takes_the_position_it_is_written_in():
    """One configuration, two spellings, opposite signs -- because the hydrogen moved.

    In stored-slot order the start atom is slot 0, so creating the carbon first puts it at the head
    of the string with its hydrogen FIRST in the written order, and creating fluorine first gives the
    carbon a parent and puts the hydrogen SECOND.  That is one transposition, so the same
    configuration must be spelled with opposite signs.  A writer that put the hydrogen last
    unconditionally would emit the same sign twice and be wrong for exactly one of the two.

    The measurement this encodes (RDKit 2026.03.4, 2026-09-02): `[C@H](F)(Cl)Br` and
    `F[C@@H](Cl)Br` are one molecule, and `[C@H](F)(Cl)Br` and `F[C@H](Cl)Br` are two.
    """
    atoms, bonds = BROMOCHLOROFLUOROMETHANE
    leading, lsids = build(atoms, bonds, order=[0, 1, 2, 3])
    parented, psids = build(atoms, bonds, order=[1, 0, 2, 3])
    configure(leading, lsids[0], _frame(lsids, (1, 2, 3, None)), 2)
    configure(parented, psids[0], _frame(psids, (1, 2, 3, None)), 2)
    first = write_smiles(leading, 'i')
    second = write_smiles(parented, 'i')
    assert first == '[C@@H](F)(Cl)Br'
    assert second == 'F[C@H](Cl)Br'
    assert sign_in(first) != sign_in(second)

    chem = importorskip('rdkit.Chem')
    assert chem.MolToSmiles(chem.MolFromSmiles(first)) == \
        chem.MolToSmiles(chem.MolFromSmiles(second))


def test_the_lone_pair_has_no_positional_rule_and_sits_last():
    """A sulfoxide's fourth direction does NOT move when the sulfur leads its component.

    Measured the same day: `[S@](=O)(C)CC` and `O=[S@](C)CC` are one molecule to RDKit, so unlike a
    hydrogen the lone pair keeps its place whether or not there is a preceding atom.  Last rather
    than first is the remaining choice; last is what chython 2 does (its frame is the three named
    substituents with the fourth direction fixed at the end) and it makes the lone pair's position in
    the written order equal to its position in `refs`, so the permutation is over the named
    directions alone.

    The test states it as: one configuration, sulfur leading and sulfur parented, SAME sign.  Which
    is the opposite of the hydrogen test above, and that contrast is the whole content.
    """
    atoms, bonds = SULFOXIDE
    leading, lsids = build(atoms, bonds, order=[0, 1, 2, 3, 4])
    parented, psids = build(atoms, bonds, order=[1, 0, 2, 3, 4])
    configure(leading, lsids[0], _frame(lsids, (1, 2, 3, None)), 2)
    configure(parented, psids[0], _frame(psids, (1, 2, 3, None)), 2)
    first = write_smiles(leading, 'i')
    second = write_smiles(parented, 'i')
    assert first == '[S@](C)(=O)CC'
    assert second == 'C[S@](=O)CC'
    assert sign_in(first) == sign_in(second)

    chem = importorskip('rdkit.Chem')
    assert chem.MolToSmiles(chem.MolFromSmiles(first)) == \
        chem.MolToSmiles(chem.MolFromSmiles(second))


# ------------------------------------------------------------------------------------------------
# THE ANTI-DRIFT TEST.
@mark.parametrize('name,fixture,centre,frame', SINGLE_CENTRE)
def test_the_sign_in_the_string_is_translate_stereo_of_the_written_order(name, fixture, centre,
                                                                         frame):
    """The C sign path against `MoleculeContainer.translate_stereo`, on the writer's OWN order.

    `smw_traversal`'s `directions` key is `smw_direction_order`'s output -- the very list
    `smw_sign_of` translates the parity into -- so this compares two independent implementations of
    the same arithmetic on one agreed input: the C `translate_parity` reached through the writer, and
    the Python-facing `translate_stereo` reached through the container.  They share
    `translate_parity` and nothing else; the writer's frame construction, the perm matching and the
    parity-to-sign mapping are all only on one side.

    Run over both parities and every creation order, canonical and stored, so the frame changes under
    the test rather than being fixed by the fixture.
    """
    atoms, bonds = fixture
    for order in permutations(range(len(atoms))):
        m, sids = build(atoms, bonds, order=list(order))
        for parity in (1, 2):
            m.set_parity(sids[centre], parity)
            for spec in ('', 'i'):
                written = smw_traversal(m, spec)['directions'][sids[centre]]
                expected = m.translate_stereo(sids[centre], written)
                assert expected in (1, 2), (name, order, parity)
                # parity 2 (odd) is `@` -- the anchor, and the only place this file states it.
                assert sign_in(write_smiles(m, spec)) == (1 if expected == 2 else 2), \
                    (name, order, parity, spec, written)


# ------------------------------------------------------------------------------------------------
# THE SWEEPS.  Ruling F26.
SWEEPS = [('bromochlorofluoromethane', BROMOCHLOROFLUOROMETHANE, [(0, (1, 2, 3, None))], 12),
          ('halomethane', HALOMETHANE, [(0, (1, 2, 3, 4))], 48),
          ('sulfoxide', SULFOXIDE, [(0, (1, 2, 3, None))], 16),
          ('dichlorocyclopropane', DICHLOROCYCLOPROPANE,
           [(0, (1, 2, 3, None)), (1, (0, 2, 4, None))], 11)]


def _sweep(fixture, centres, spec):
    """Every creation order, one fixed configuration, the set of strings produced."""
    atoms, bonds = fixture
    seen = set()
    orders = list(permutations(range(len(atoms))))
    for order in orders:
        m, sids = build(atoms, bonds, order=list(order))
        for centre, frame in centres:
            configure(m, sids[centre], _frame(sids, frame), 2)
        seen.add(write_smiles(m, spec))
    return seen, len(orders)


@mark.parametrize('name,fixture,centres,stored_count', SWEEPS)
def test_canonical_stereo_output_does_not_depend_on_the_creation_order(name, fixture, centres,
                                                                       stored_count):
    """One configuration, every creation order, ONE string.  The fixture ruling F26 asked for.

    Exhaustive rather than sampled: these molecules are four and five atoms, so 24 and 120 orders is
    the whole group and there is nothing left to sample.
    """
    seen, count = _sweep(fixture, centres, '')
    assert len(seen) == 1, (name, count, sorted(seen)[:4])
    assert count == factorial(len(fixture[0]))


@mark.parametrize('name,fixture,centres,stored_count', SWEEPS)
def test_stored_order_stereo_is_not_creation_order_invariant(name, fixture, centres, stored_count):
    """The could-have-failed evidence for the test above (ruling F102), as a NUMBER.

    Stored order is the creation order by definition, so the same sweep must produce many strings --
    twelve, forty-eight, sixteen and eleven, and the counts are asserted rather than just `> 1` so that a
    change which quietly collapses them is a failure and not a silent weakening.  If this ever
    reports one string, the test above has stopped measuring anything.
    """
    seen, _ = _sweep(fixture, centres, 'i')
    assert len(seen) == stored_count, (name, sorted(seen))


@mark.parametrize('name,fixture,centres,stored_count', SWEEPS)
def test_every_stored_order_spelling_is_the_same_molecule_to_rdkit(name, fixture, centres,
                                                                   stored_count):
    """The other half of the sweep, and the one that needs an oracle.

    That the canonical path gives one string says the writer is CONSISTENT.  That all forty-eight (or
    sixteen, or eleven) stored-order spellings are one molecule to RDKit says it is RIGHT: a frame
    error would make some of them enantiomers of the others, which no amount of internal agreement
    could detect.  This is the test that would have caught the defect ruling F26 was written about.
    """
    chem = importorskip('rdkit.Chem')
    seen, _ = _sweep(fixture, centres, 'i')
    assert len(seen) == stored_count, name
    canonical = {chem.MolToSmiles(chem.MolFromSmiles(s)) for s in seen}
    assert len(canonical) == 1, (name, sorted(canonical))


# ------------------------------------------------------------------------------------------------
# THE STEREO SEED: a constitutional symmetry that the CONFIGURATION breaks.
#
# 2,3-dichlorobutane is the tetrahedral case of the defect the cis/trans suite's diene shows on a
# double bond.  Its constitution is symmetric end to end -- swap C1 with C4, C2 with C3 and the two
# chlorines, and the graph maps onto itself -- so the stereo-blind refinement leaves the two centres
# in one class and the extremal search has a tie to break.  Give the two centres OPPOSITE
# configurations (the meso compound) and that swap is no longer a symmetry of the molecule, but the
# order still does not know it: whichever centre the creation order put first came out first.  The
# seed is what tells it, and this block reads the mechanism directly instead of inferring it from a
# string.
DICHLOROBUTANE = ([(6, 3), (6, 1), (6, 1), (6, 3), (17, 0), (17, 0)],
                  [(0, 1, 1), (1, 2, 1), (2, 3, 1), (1, 4, 1), (2, 5, 1)])
# The two frames are MIRROR IMAGES of each other under that swap -- (C1, C3, Cl) at one centre and
# (C4, C2, Cl) at the other -- so equal `want` values mean the swap carries one configuration onto
# the other, and the pair (2, 1) is the meso compound while (2, 2) and (1, 1) are the enantiomers.
DICHLOROBUTANE_FRAMES = [(1, (0, 2, 4, None)), (2, (3, 1, 5, None))]


def _dichlorobutane(order, wants):
    m, sids = build(*DICHLOROBUTANE, order=order)
    for (centre, frame), want in zip(DICHLOROBUTANE_FRAMES, wants):
        configure(m, sids[centre], _frame(sids, frame), want)
    return m, sids


def _dichlorobutane_sweep(wants, spec=''):
    seen = set()
    for order in permutations(range(6)):
        m, _ = _dichlorobutane(list(order), wants)
        seen.add(write_smiles(m, spec))
    return seen


def test_a_meso_compound_gives_one_string():
    """720 creation orders, one string -- the tetrahedral half of what the stereo seed fixes."""
    assert len(_dichlorobutane_sweep((2, 1))) == 1, sorted(_dichlorobutane_sweep((2, 1)))


def test_the_meso_stored_order_spellings_are_one_molecule_to_rdkit():
    """The could-have-failed number (forty) and the oracle, in one test because they share the sweep.

    Forty stored-order spellings of one compound: a frame error here would show up as some of them
    being the chiral diastereomer rather than the meso one, which is a difference RDKit sees and
    internal agreement cannot.
    """
    chem = importorskip('rdkit.Chem')
    seen = _dichlorobutane_sweep((2, 1), 'i')
    assert len(seen) == 40, sorted(seen)[:4]
    assert len({chem.MolToSmiles(chem.MolFromSmiles(s)) for s in seen}) == 1, sorted(seen)[:4]


def test_the_four_dichlorobutane_configurations_are_three_compounds():
    """Two enantiomers and one meso: four combinations, THREE strings, and the collapse is asserted.

    A seed that over-separated would give these four strings, each one perfectly invariant over
    creation orders, so every sweep in this file would still pass.  The meso compound is the same
    substance whichever centre is called R, so its two spellings have to be one string.
    """
    chem = importorskip('rdkit.Chem')
    ours = {}
    for wants in ((2, 2), (1, 1), (2, 1), (1, 2)):
        seen = _dichlorobutane_sweep(wants)
        assert len(seen) == 1, (wants, sorted(seen))
        ours[wants] = seen.pop()
    assert len(set(ours.values())) == 3, ours
    assert ours[(2, 1)] == ours[(1, 2)], ours
    assert ours[(2, 2)] != ours[(1, 1)], ours          # enantiomers: two compounds, two strings
    theirs = {wants: chem.MolToSmiles(chem.MolFromSmiles(s)) for wants, s in ours.items()}
    assert len(set(theirs.values())) == 3, theirs
    assert theirs[(2, 1)] == theirs[(1, 2)], theirs


def test_the_canonical_order_pins_the_MESO_halves_by_ITSELF_and_leaves_the_chiral_pair_tied():
    """The mechanism, measured through the public order rather than through a string (ruling F102).

    `canonical_order()` is NOT STEREO-BLIND -- the search's leaf certificate carries a parity tail and
    its orbit prune refines by parity -- so the unseeded order pins the meso halves on its own, over
    all 720 creation orders, and the seed agrees with it rather than rescuing it.  A stereo-blind order
    puts EITHER centre first depending on the creation order, the two centres being one refinement
    class with the tie falling to slot order, and then only the seed can pin them.  Both are asserted,
    because the seed is a public argument and has to be sigma-equivariant whether or not anything
    depends on it; `smw_stereo_seed_labels` is subsumed for this case.

    The chiral diastereomers are the control and their answer is the interesting one: the ambiguity
    stays, and must, because there the swap IS an automorphism of the configured molecule, so the two
    labellings describe the same string and nothing needs breaking.  A search that pinned this case
    too would be inventing an asymmetry -- ruling F95's sigma-equivariance requirement failing in the
    direction that still looks like success -- so the new prune is required to leave it alone, and this
    is where that is checked.  Both enantiomers are swept, not just one, because a prune that broke
    the tie in a handedness-DEPENDENT way would pass on a single row.

    The last assertion is the one that says the pin tracks CONFIGURATION and not numbering: the meso
    compound presented as `(2, 1)` and as `(1, 2)` is the same substance with the fixture's two centres
    exchanged, so a correct order must reach the opposite answer for the two -- the same absolute
    labelling, read through a fixture whose own numbering flipped.  Equal answers there would mean the
    order was still keying on slots.
    """
    answers = {}
    for wants, values in (((2, 1), 1), ((1, 2), 1), ((2, 2), 2), ((1, 1), 2)):
        blind = set()
        seeded = set()
        for order in permutations(range(6)):
            m, sids = _dichlorobutane(list(order), wants)
            positions = m.canonical_order()
            blind.add(positions[sids[1]] < positions[sids[2]])
            positions = m.canonical_order(smw_stereo_seed_labels(m))
            seeded.add(positions[sids[1]] < positions[sids[2]])
        assert len(blind) == values, (wants, blind)
        assert len(seeded) == values, (wants, seeded)
        assert blind == seeded, (wants, blind, seeded)
        answers[wants] = blind
    assert answers[(2, 1)] != answers[(1, 2)], answers


def test_the_seed_splits_the_meso_centres_and_leaves_the_chiral_ones_tied():
    """The labels themselves: all six distinct for meso, palindromic for the enantiomer.

    Two claims one string cannot make. `atoms_order` is the same three classes in both cases, so the
    difference is entirely the seed's, and the palindrome in the second case is the automorphism still
    standing where it should.
    """
    m, sids = _dichlorobutane(None, (2, 1))
    classes = [m.atoms_order[sids[j]] for j in range(6)]
    assert classes == [2, 1, 1, 2, 3, 3], classes
    labels = smw_stereo_seed_labels(m)
    assert len({labels[sids[j]] for j in range(6)}) == 6, labels

    m, sids = _dichlorobutane(None, (2, 2))
    assert [m.atoms_order[sids[j]] for j in range(6)] == classes
    labels = [smw_stereo_seed_labels(m)[sids[j]] for j in range(6)]
    # The swap sigma written out rather than a reversal: sigma is (C1 C4)(C2 C3)(Cl Cl), which is not
    # index reversal for this fixture's bond list, and a reversal check would pass for the wrong reason.
    for j, k in ((0, 3), (1, 2), (4, 5)):
        assert labels[j] == labels[k], (j, k, labels)
    assert len(set(labels)) == 3, labels


# ------------------------------------------------------------------------------------------------
# WHAT IS AND IS NOT WRITTEN.
def test_no_sign_under_the_no_stereo_key():
    atoms, bonds = BROMOCHLOROFLUOROMETHANE
    m, sids = build(atoms, bonds)
    m.set_parity(sids[0], 2)
    assert write_smiles(m, '!s') == 'C(F)(Cl)Br'
    assert smw_traversal(m, '!s')['directions'] == {}


def test_a_stated_parity_is_written_even_where_it_is_not_stereogenic():
    """Fidelity: the writer spells what the molecule HOLDS, and does not audit it.

    Dichlorofluoromethane's carbon anchors a unit -- four directions, two of them the same chlorine
    twice over -- and `stereogenic_units()` correctly refuses it, so a parity stated on it means
    nothing.  It is still written.  Dropping it would be a silent edit of the input, and the caller
    who wants to know already has `stereo_rejections`; the reader stores what it reads for the same
    reason, so the round trip is stable in both directions.

    Overturning this is one condition in `smw_sign_of` -- it is a decision, not an accident.
    """
    m, sids = build([(6, 1), (9, 0), (17, 0), (17, 0)], [(0, 1, 1), (0, 2, 1), (0, 3, 1)])
    assert [u['anchor'] for u in m.stereogenic_units()] == []
    assert any(u['anchor'] == sids[0] for u in m.stereo_units())
    m.set_parity(sids[0], 2)
    assert sign_in(write_smiles(m)) != 0


def test_cis_trans_configuration_is_written_as_a_direction_and_never_as_a_sign():
    """A cis/trans unit's configuration is a property of a BOND, so no `@` may appear for it.

    `smw_sign_of`'s tetrahedral path must not be reached for kind 1; reached, it puts a sign on an atom
    where the sign means nothing.  The direction itself is covered in full by
    test_smiles_write_cis_trans.py, including the convention it rests on.
    """
    atoms, bonds = BUTENE
    m, sids = build(atoms, bonds)
    anchors = [u['anchor'] for u in m.stereogenic_units() if u['kind'] != 0]
    assert anchors, 'but-2-ene must have a stereogenic double bond'
    m.set_parity(anchors[0], 2)
    written = write_smiles(m)
    assert '@' not in written, written
    assert written.count('/') + written.count('\\') == 2, written


def test_the_sweep_helper_does_not_silently_accept_a_frame_it_cannot_state():
    """`configure` propagates rather than guessing -- checked, since every sweep in this file trusts it.

    For a LEGAL frame both parities are reachable, so `configure` can only fail by
    `translate_stereo` refusing, and refusing is what a frame naming a non-neighbour must do.  If
    that ever became silent, every sweep above would still pass while measuring an unstated
    configuration.
    """
    atoms, bonds = BROMOCHLOROFLUOROMETHANE
    m, sids = build(atoms, bonds)
    stranger = m.add_atom('C')
    try:
        configure(m, sids[0], (sids[1], sids[2], stranger, None), 2)
    except (AssertionError, KeyError, ValueError):
        return
    raise AssertionError('configure accepted a frame naming a non-neighbour')


def test_a_random_creation_order_sample_agrees_with_the_exhaustive_sweep():
    """The sampling helper the constitution file uses, applied to stereo on a bigger molecule.

    Six atoms is 720 orders, which is still exhaustive; the point is that the sweep shape scales past
    the fixtures above without the answer changing.
    """
    atoms = [(6, 1), (6, 1), (6, 2), (17, 0), (17, 0), (6, 3)]
    bonds = [(0, 1, 1), (1, 2, 1), (2, 0, 1), (0, 3, 1), (1, 4, 1), (2, 5, 1)]
    rng = Random(20260902)
    seen = set()
    for _ in range(200):
        order = list(range(len(atoms)))
        rng.shuffle(order)
        m, sids = build(atoms, bonds, order=order)
        configure(m, sids[0], _frame(sids, (1, 2, 3, None)), 2)
        configure(m, sids[1], _frame(sids, (0, 2, 4, None)), 2)
        seen.add(write_smiles(m))
    assert len(seen) == 1, sorted(seen)[:4]
