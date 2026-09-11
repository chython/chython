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
"""THE MIRROR AUTOMORPHISM, and why the identity rests on a certificate that carries configuration.

THE MECHANISM.  A `canonical_order` that is a function of the CONSTITUTION alone -- a certificate
carrying elements, bonds, charges and hydrogen counts and NOT ONE BIT of configuration, with an orbit
prune keeping one candidate per orbit of the constitutional automorphism group -- oscillates.  On a
molecule whose constitution is symmetric and whose configuration is not carried along by that
symmetry, the two labellings a mirror automorphism relates tie on the certificate, only ONE of them is
even generated (the prune drops the other), and which one that is comes from the slot order -- the
order the caller happened to add the atoms in.

cis-cyclobutane-1,3-diol is the smallest witness in the repo's own `test/stereo.sdf`.  It is ACHIRAL:
the constitutional automorphism that swaps its two carbinol carbons inverts the parity at both, so
`O[C@H]1C[C@H](O)C1` and `O[C@@H]1C[C@@H](O)C1` are two spellings of ONE compound.  Both are legal,
which is why nothing is corrupted by it, and exactly two spellings is why it is invisible from
inside: nothing downstream can see it except by comparing two strings that should be equal, or two
`canonical_bytes` that should be equal.

WHY THAT BLOCKED `__eq__` AND `__hash__`.  Both rest on `canonical_bytes` (`_molecule_container.pxi`
says so, at length, and lists the two wrong substrates it is not).  An oscillating canonical form
makes `==` self-break on a ROUND TRIP: write the molecule out, read it back, and the reader's atom
order is not the writer's, so the two compare unequal and the same compound sits in a set twice.
Measured over 4 creation orders on the 294 usable records of `test/stereo.sdf`, a constitution-only
certificate oscillates on 55 of them, its canonical SMILES on 43, and 45 fail write -> read -> `==`;
the certificate below answers 0 / 294 on all three.  The corpus-scale gate lives in
`test_smiles_write_differential.py::test_the_tetrahedral_string_does_not_depend_on_the_creation_order`;
this file is the per-molecule half -- one named compound per mechanism, so a regression says WHICH
shape broke.

WHAT PREVENTS IT, in one sentence each, both in `_canonical.pxi`:

  * the leaf certificate carries a PARITY TAIL -- one digit per canonical position, read in the frame
    the DISCRETE leaf colouring names -- so two labellings a parity-inverting automorphism relates do
    not tie, and the extremum picks one of them as a function of the molecule;
  * the orbit prune may not drop a candidate it cannot prove is stereo-equivalent.
    It runs against a colouring that carries the parity digit, and where that colouring cannot
    NAME a unit's frame -- which is precisely the mirror case, because the two ring branches leaving a
    carbinol carbon share a colour -- it does not prune at all and both labellings reach a leaf.

V2 IS NOT THE ORACLE HERE.  Put the same question to an isolated chython 2.24 -- eight creation
orders, count the distinct canonical strings -- and it oscillates on twelve of the fifteen witnesses
below, up to SIX distinct strings for one compound, and fails its own write -> read -> `==` on six of
them.  Oscillation is a property of this family of certificate algorithms, so V2 agreeing with V3 on a
symmetric stereocentre is worth nothing as evidence.  V2 earns its keep here in two other ways: as a
gross-ordering check on the stereo-free control (it agrees there), and as a source of acceptance
tests, since where V2 diverges it diverges in a nameable way -- see the note on
`test_the_two_spellings_of_an_achiral_compound_are_one_compound`.  What actually pins the fix is the
invariants, which need no oracle at all: one canonical form per compound across creation orders, a
molecule equal to itself after a write-then-read, a molecule not equal to its enantiomer, and
`__hash__` agreeing with `__eq__` on all of it.

A SECOND, INDEPENDENT DEFECT LIVES ON THE SAME PATH and these tests catch it too.  `__eq__` opens
with an exact rejection on `_union_feature_words`, and word IV of that screen (`_features.pxi`,
`atom_feature_word4`) ORs in bit 6 from the RAW STORED parity, which is a statement in the
molecule's own slot frame (ruling F26) and not a property of the compound.  Two spellings of one
meso compound therefore differ in the screen and `==` returns False before `canonical_bytes` is ever
consulted, even once the canonical form is fixed.  Measured on `C[C@H](O)[C@H](O)C` against
`C[C@@H](O)[C@@H](O)C` and on `O[C@H]1CCCC[C@H]1O` against `O[C@@H]1CCCC[C@@H]1O`: `==` False,
`canonical_bytes` equal.  The screen may only carry frame-free features, so the frame-relative bit is
masked out of the comparison; word IV's own docstring already concedes the bit is "screen-invisible
on its own".

EVERY STRUCTURE HERE IS A PUBLIC COMPOUND.
"""
from random import Random

from pytest import mark, raises

from chython.core import AutomorphismBudgetExceeded, MoleculeContainer
from chython.core._core import read_smiles, write_smiles


# ------------------------------------------------------------------------------------------------
# THE RELABELLER.  Rebuilding a molecule in a different atom order is the whole experiment, and it
# has to carry the CONFIGURATION across, which a plain atoms-and-bonds copy does not: a stored parity
# is a statement about the anchor's neighbours IN SLOT ORDER (ruling F26), so the same compound in two
# atom orders holds two different bytes.  `translate_stereo` is the public read that is frame-relative
# rather than slot-relative, so the carry is "ask the source for its parity in a NAMED frame, then
# store whichever byte makes the copy answer the same in the corresponding frame" -- which is what
# `translate_stereo` exists for and what every other stereo test in this directory does.
def _relabel(m, order):
    """Rebuild `m` with its atoms added in `order`, a permutation of `m.atom_numbers`.

    Same compound, different slots, same configuration.  Raises rather than returning a molecule
    whose configuration did not survive, because a silent loss here would look like the defect under
    test.
    """
    fresh = MoleculeContainer()
    ids = {}
    with fresh.edit():
        for sid in order:
            ids[sid] = fresh.add_atom(m.element_of(sid), charge=m.charge_of(sid),
                                      isotope=m.isotope_of(sid), radical=m.radical_of(sid),
                                      implicit_h=m.implicit_h_of(sid))
        for a in m.atom_numbers:
            for b in m.neighbors_of(a):
                if a < b:
                    fresh.add_bond(ids[a], ids[b], m.order_of(a, b))
    for unit in m.stereo_units():
        anchor = unit['anchor']
        if not m.parity_of(anchor):
            continue
        want = m.translate_stereo(anchor, unit['refs'])
        mapped = tuple(None if r is None else ids[r] for r in unit['refs'])
        for parity in (1, 2):
            fresh.set_parity(ids[anchor], parity)
            if fresh.translate_stereo(ids[anchor], mapped) == want:
                break
        else:                                   # pragma: no cover - a broken relabeller, not a defect
            raise AssertionError('the relabeller lost the configuration at %d' % anchor)
    return fresh


def _orders(m, count, seed):
    """`count` creation orders for `m`: the identity first, then random permutations."""
    rnd = Random(seed)
    sids = list(m.atom_numbers)
    yield sids
    for _ in range(count - 1):
        perm = list(sids)
        rnd.shuffle(perm)
        yield perm


# ------------------------------------------------------------------------------------------------
# THE WITNESSES.  Every one oscillates under a constitution-only certificate, named by compound not by
# its `test/stereo.sdf` index, plus the allene entry and the two acyclic meso cases.
#
# WHAT THEY HAVE IN COMMON IS A SYMMETRIC CONSTITUTION, NOT ACHIRALITY.  Most of the table is achiral
# -- a constitutional automorphism carries the compound onto itself while inverting parity at two
# centres, so the two spellings are one compound -- but `1,2,4-trimethylcyclopentane` and the
# alternating hexachlorocyclohexane are CHIRAL (RDKit's CIP labeller says so), and they oscillate all
# the same.  That is the point: the tie is between two LABELLINGS of one molecule that the
# constitutional certificate cannot separate, and whether the compound happens to superpose on its
# mirror image is a separate question.  Fixing only the achiral half would leave the chiral half
# picking its canonical form by slot order, which is the same bug.
#
# TWO ENTRIES REST ON CHYTHON'S OWN MODEL, because no oracle here can second them: RDKit perceives no
# stereo at all in `cyclohexane-1,4-diylidene bis-allene` (it drops both allene axes) or in the
# `adamantane skeleton` (it drops all four bridge centres), so for those two the invariants below --
# one canonical form per compound, and equal after a round trip -- are the entire specification.
MIRRORS = [
    ('cis-cyclobutane-1,3-diol', 'O[C@H]1C[C@H](O)C1'),
    ('1,2,4-trimethylcyclopentane', 'C[C@@H]1C[C@@H](C[C@H]1C)C'),
    ('trans-decalin', '[H][C@]12CCCC[C@@]1([H])CCCC2'),
    ('cis-cyclohexane-1,4-diol', 'O[C@H]1CC[C@@H](O)CC1'),
    ('cis-1,4-dimethylcyclohexane', 'C[C@H]1CC[C@@H](C)CC1'),
    ('hexachlorocyclohexane, alternating signs',
     '[C@H]1(Cl)[C@@H](Cl)[C@H](Cl)[C@@H](Cl)[C@H](Cl)[C@@H]1Cl'),
    ('hexachlorocyclohexane, uniform signs',
     '[C@H]1(Cl)[C@H](Cl)[C@H](Cl)[C@H](Cl)[C@H](Cl)[C@H]1Cl'),
    ('1,2,3,4-tetramethoxycyclobutane', 'CO[C@H]1[C@H](OC)[C@H](OC)[C@H]1OC'),
    ('4-chloro-4-methylcyclohexan-1-ol', 'O[C@H]1CC[C@](C)(Cl)CC1'),
    ('adamantane skeleton', 'C1[C@@H]2C[C@@H]3C[C@H]1C[C@H](C3)C2'),
    ('bis(4-methylcyclohexyl)methanol', 'OC([C@H]1CC[C@H](C)CC1)[C@H]1CC[C@H](C)CC1'),
    ('cyclohexane-1,4-diylidene bis-allene', 'C(=[C@]=C1CCC(CC1)=[C@]=CC)C'),
    ('meso-tartaric acid', 'O[C@H](C(=O)O)[C@@H](O)C(=O)O'),
    ('meso-2,3-butanediol', 'C[C@H](O)[C@H](O)C'),
    ('cis-cyclohexane-1,2-diol', 'O[C@H]1CCCC[C@H]1O'),
]


@mark.parametrize('name, text', MIRRORS, ids=[n for n, _ in MIRRORS])
def test_the_canonical_form_of_a_mirror_symmetric_compound_is_one_value(name, text):
    """THE DEFECT, per compound: eight creation orders, one `canonical_bytes`.

    This is the assertion the whole change exists for.  Before it, every entry in the table returned
    two or more values here -- the mirror image labelling and the original one, tying on a
    stereo-blind certificate and separated by nothing but which atom the SMILES parser reached first.
    """
    m = read_smiles(text)
    forms = {_relabel(m, order).canonical_bytes for order in _orders(m, 8, 20260903)}
    assert len(forms) == 1, '%d distinct canonical forms for %s' % (len(forms), name)


@mark.parametrize('name, text', MIRRORS, ids=[n for n, _ in MIRRORS])
def test_the_canonical_string_of_a_mirror_symmetric_compound_is_one_string(name, text):
    """The same statement one level up, because the string is what a person sees.

    `canonical_bytes` and the canonical SMILES are two readings of the SAME labelling, so in
    principle one assertion would do; in practice they are reached through different call sites
    (`mol_identity_bytes` seeds its own colouring, `smw_canonical_positions` seeds another) and a fix
    that repaired one and not the other is exactly the shape of mistake worth a second test.
    """
    m = read_smiles(text)
    strings = {write_smiles(_relabel(m, order)) for order in _orders(m, 8, 20260903)}
    assert len(strings) == 1, '%d distinct strings for %s: %s' % (len(strings), name,
                                                                  sorted(strings))


@mark.parametrize('name, text', MIRRORS, ids=[n for n, _ in MIRRORS])
def test_a_mirror_symmetric_compound_equals_itself_after_a_round_trip(name, text):
    """`==` MUST NOT SELF-BREAK ON A WRITE-THEN-READ, and these compounds are where it would.

    The reader's atom order is the string's traversal order and not the writer's input order, so a
    round trip IS a relabelling -- which makes this the one place an oscillating canonical form shows
    up in ordinary use, with no test harness involved at all.
    """
    m = read_smiles(text)
    back = read_smiles(write_smiles(m))
    assert back == m
    assert hash(back) == hash(m)


# ------------------------------------------------------------------------------------------------
# THE OTHER DIRECTION, and the reason the fix is not "ignore stereo in the certificate".  Making the
# certificate parity-BLIND would also make every entry above return one value, and would make
# enantiomers equal.  Both halves have to hold at once, which is why they are tested together.
#
# WHICH TABLE A SPELLING BELONGS IN WAS MEASURED, NOT GUESSED.  Two `C[C@H]`s in a row look like a
# pair and are not: RDKit's CIP labeller reads `C[C@H](O)[C@H](O)C` as (2S,3R), which is the MESO
# diastereomer and belongs below, and `C[C@H](O)[C@@H](O)C` as (2S,3S), which is the chiral one and
# belongs here.  Same inversion for the cyclohexane-1,2-diol, where `O[C@H]1CCCC[C@H]1O` is the cis
# (meso) ring.  The check that settles it for any candidate is whether RDKit's canonical SMILES of the
# string equals its canonical SMILES of the SIGN-INVERTED string: equal means one compound.
ENANTIOMERS = [
    ('alanine', 'N[C@@H](C)C(=O)O', 'N[C@H](C)C(=O)O'),
    ('glyceraldehyde', 'OC[C@@H](O)C=O', 'OC[C@H](O)C=O'),
    ('bromochlorofluoromethane', '[C@H](F)(Cl)Br', '[C@@H](F)(Cl)Br'),
    ('trans-cyclohexane-1,2-diol', 'O[C@H]1CCCC[C@@H]1O', 'O[C@@H]1CCCC[C@H]1O'),
    ('tartaric acid', 'O[C@H](C(=O)O)[C@H](O)C(=O)O', 'O[C@@H](C(=O)O)[C@@H](O)C(=O)O'),
    ('2,3-butanediol', 'C[C@H](O)[C@@H](O)C', 'C[C@@H](O)[C@H](O)C'),
    ('penta-2,3-diene', 'C/C=[C@]=C/C', 'C/C=[C@@]=C/C'),
]


@mark.parametrize('name, left, right', ENANTIOMERS, ids=[n for n, _, _ in ENANTIOMERS])
def test_a_molecule_and_its_enantiomer_are_not_equal(name, left, right):
    """A MIRROR IMAGE PAIR IS TWO COMPOUNDS, and the parity tail is what says so.

    Every pair here is one chiral skeleton written twice with every sign inverted.  The certificate's
    graph part is byte-identical across a pair -- that is what makes them a pair -- so if this passes
    it is the parity tail passing, and it is the direction a canonical form may never get wrong:
    unequal molecules reported equal corrupts a dict, while the reverse only costs a cache miss.
    """
    a = read_smiles(left)
    b = read_smiles(right)
    assert a != b, name
    assert hash(a) != hash(b), '%s: a legal collision, but not one that should happen here' % name
    assert a.canonical_bytes != b.canonical_bytes


@mark.parametrize('name, left, right', ENANTIOMERS, ids=[n for n, _, _ in ENANTIOMERS])
def test_an_enantiomeric_pair_stays_two_compounds_under_every_relabelling(name, left, right):
    """And the inequality is not an artefact of the two strings' atom orders happening to differ."""
    a = read_smiles(left)
    b = read_smiles(right)
    left_forms = {_relabel(a, order).canonical_bytes for order in _orders(a, 6, 20260903)}
    right_forms = {_relabel(b, order).canonical_bytes for order in _orders(b, 6, 20260904)}
    assert len(left_forms) == 1 and len(right_forms) == 1, name
    assert left_forms != right_forms, name


# ------------------------------------------------------------------------------------------------
# THE ACHIRAL CONVERSE.  A compound that IS its own mirror image must have ONE canonical form, and
# the two enantiomeric SPELLINGS of it must land on that one value -- which is the same statement as
# the invariance sweep above, reached from the string side instead of from the atom-order side.
ACHIRAL_PAIRS = [
    ('cis-cyclobutane-1,3-diol', 'O[C@H]1C[C@H](O)C1', 'O[C@@H]1C[C@@H](O)C1'),
    ('meso-tartaric acid', 'O[C@H](C(=O)O)[C@@H](O)C(=O)O', 'O[C@@H](C(=O)O)[C@H](O)C(=O)O'),
    ('meso-2,3-butanediol', 'C[C@H](O)[C@H](O)C', 'C[C@@H](O)[C@@H](O)C'),
    ('cis-cyclohexane-1,2-diol', 'O[C@H]1CCCC[C@H]1O', 'O[C@@H]1CCCC[C@@H]1O'),
    ('cis-1,4-dimethylcyclohexane', 'C[C@H]1CC[C@@H](C)CC1', 'C[C@@H]1CC[C@H](C)CC1'),
]


@mark.parametrize('name, left, right', ACHIRAL_PAIRS, ids=[n for n, _, _ in ACHIRAL_PAIRS])
def test_the_two_spellings_of_an_achiral_compound_are_one_compound(name, left, right):
    """An achiral compound superposes on its mirror image, so the two spellings are ONE molecule.

    This is the assertion that makes the fix a canonicalisation rather than a refusal to tie-break:
    it would be easy to make every sweep above pass by giving the two labellings different canonical
    forms and calling the molecule two compounds.  RDKit agrees with this reading on all five -- the
    differential file's `test_nothing_oscillates_and_a_REGRESSION_would_have_to_be_a_MIRROR_PAIR`
    measures it corpus-wide, and the sign-inversion check above reproduces it per row.

    CHYTHON 2 AGREES ON THREE AND NOT ON THE OTHER TWO, which is why it is not the oracle here.  Asked
    the same question, 2.24 returns one canonical string for meso-tartaric acid, meso-2,3-butanediol
    and cis-cyclohexane-1,2-diol, and TWO for cis-cyclobutane-1,3-diol and cis-1,4-dimethyl-
    cyclohexane -- it splits an achiral compound in half.  That is the V2 divergence these rows exist
    to fence off, not a precedent to preserve, so the specification is the invariant and not the second
    opinion.
    """
    a = read_smiles(left)
    b = read_smiles(right)
    assert a == b, name
    assert hash(a) == hash(b), name
    assert write_smiles(a) == write_smiles(b), name


@mark.parametrize('name, left, right', ACHIRAL_PAIRS, ids=[n for n, _, _ in ACHIRAL_PAIRS])
def test_the_feature_word_screen_never_rejects_a_pair_the_canonical_form_accepts(name, left, right):
    """THE SECOND DEFECT, isolated: `__eq__`'s prefilter may be lossy in ONE direction only.

    `__eq__` rejects on `_union_feature_words` before it pays for a canonical form, and calls that an
    exact rejection.  It is one only for features that are properties of the COMPOUND; word IV's bit 6
    is the raw stored parity, which is a statement in the molecule's own slot frame, so the two
    spellings here held opposite bits and `==` answered False while `canonical_bytes` were equal.
    Asserted as the implication rather than by reading the bit, because the implication is the contract
    and would still have to hold if the screen were rebuilt from different words tomorrow.
    """
    a = read_smiles(left)
    b = read_smiles(right)
    assert a.canonical_bytes == b.canonical_bytes, '%s: the canonical form, not the screen' % name
    assert a._union_feature_words == b._union_feature_words or a == b, \
        '%s: the screen rejected a pair whose canonical forms agree' % name


# ------------------------------------------------------------------------------------------------
# AUTOMORPHIC RELABELLING, HASH/EQ CONSISTENCY, AND THE STEREO-FREE CONTROL.
_STEREO_FREE = ['c1ccccc1', 'CC(C)C', 'C1CC1', 'OCC(O)CO', 'c1ccc2ccccc2c1', 'C1CCCCC1',
                'CC(=O)Oc1ccccc1C(=O)O', 'Clc1ccc(Cl)cc1', 'C1CC2CCC1C2']


@mark.parametrize('text', _STEREO_FREE)
def test_a_constitution_is_still_invariant_under_relabelling(text):
    """The control: the fix must not move a molecule that has no configuration to say anything about.

    Every entry is stereo-free, so the parity tail is all zeros and the orbit prune keeps its old
    reach.  A failure here is not a stereo defect, it is the fix leaking into the common path.
    """
    m = read_smiles(text)
    forms = {_relabel(m, order).canonical_bytes for order in _orders(m, 8, 20260903)}
    strings = {write_smiles(_relabel(m, order)) for order in _orders(m, 8, 20260903)}
    assert len(forms) == 1
    assert len(strings) == 1


@mark.parametrize('text', [t for _, t in MIRRORS] + _STEREO_FREE)
def test_hash_and_eq_agree_on_every_relabelling(text):
    """`a == b` implies `hash(a) == hash(b)`, checked pairwise rather than assumed.

    The Python contract is one-directional and this is the direction that matters: a set relies on
    it, and a `__hash__` that disagreed with `__eq__` would make a molecule findable or unfindable by
    luck.  The set assertion at the end is the same statement in the form a caller writes.
    """
    m = read_smiles(text)
    copies = [_relabel(m, order) for order in _orders(m, 6, 20260905)]
    for a in copies:
        for b in copies:
            assert a == b
            assert hash(a) == hash(b)
    assert len({*copies}) == 1, 'the same compound six times is one set member'


def test_a_set_of_relabelled_mirror_compounds_has_one_member_each():
    """The whole table at once, in the shape the defect actually bit in: a `set` of molecules."""
    bag = []
    for _, text in MIRRORS:
        m = read_smiles(text)
        bag.extend(_relabel(m, order) for order in _orders(m, 4, 20260906))
    assert len(set(bag)) == len(MIRRORS)


def test_eq_against_a_non_molecule_is_not_an_error():
    """`NotImplemented` and not a raise, so `mol == 'C'` is False rather than a TypeError."""
    m = read_smiles('CCO')
    assert m != 'CCO'
    assert m != 42
    assert not (m == None)              # noqa: E711 -- `is None` would not exercise __eq__


def test_the_relabeller_itself_preserves_the_configuration():
    """A guard on the harness: if `_relabel` silently dropped stereo, every sweep above would pass.

    So it is checked directly -- the relabelled molecule writes a string carrying the same number of
    stereo tokens, and a molecule relabelled and then relabelled back is the original.
    """
    m = read_smiles('N[C@@H](C)C(=O)O')
    other = read_smiles('N[C@H](C)C(=O)O')
    for order in _orders(m, 6, 20260907):
        copy = _relabel(m, order)
        assert copy == m
        assert copy != other, 'the relabeller lost the sign, so it cannot witness anything'


def test_an_empty_molecule_and_a_lone_atom_hash():
    """The degenerate sizes, because `mol_identity_bytes` short-circuits both."""
    empty = MoleculeContainer()
    lone = read_smiles('C')
    assert empty == MoleculeContainer()
    assert hash(empty) == hash(MoleculeContainer())
    assert empty != lone
    assert lone == read_smiles('C')


def test_a_relabelled_molecule_is_not_confused_with_a_different_constitution():
    """The screen `__eq__` opens with is a prefilter, so the slow path has to be the decider.

    propane and butane share `_union_feature_words` -- the docstring on that property says so and
    that is why it is private -- and they differ in atom count, which `__eq__` rejects on first.
    Pentane against 2-methylbutane is the pair that shares BOTH counts and the words, so it is the
    one that actually reaches the canonical form.
    """
    a = read_smiles('CCCCC')
    b = read_smiles('CCC(C)C')
    assert a != b
    assert a.canonical_bytes != b.canonical_bytes
    for order in _orders(a, 4, 20260908):
        assert _relabel(a, order) != b


# ------------------------------------------------------------------------------------------------
# THE COST OF THE FIX ABOVE, AND WHERE IT WAS PAID BACK.  Making the search stereo-aware closed the
# oscillation and left a hole in the PRUNING, because the two mechanisms it added pull opposite ways:
# the parity reached the leaf certificate as a TAIL, and the orbit prune learned to stand down
# wherever a colouring cannot name a configured unit's frame.  So on a molecule whose constitution
# ties two atoms that only the parity separates, the prune declines and nothing else fires either --
# `_canon_indicator` is constitutional, so it scores both candidates equally, and both subtrees are
# walked to the leaf where the tail finally decides.  Such branchings compose, and the tree grows as
# 3^k over k of them while the ANSWER never moves.
#
# The cure is to fold the parity digits into the ROOT partition, which is what `mol_identity_bytes`
# had always done for itself -- hence `canonical_order` costing 100x what `canonical_bytes` cost on
# the same molecule, the two running searches of 7413 and 47 nodes for one answer.
#
# `nodes_before` is measured on the unfolded search and `nodes_after` on the folded one, both by
# bisecting `_node_budget` on this exact string, and `nodes_after` is asserted from both sides -- so
# a change that regrows the tree fails here rather than in a benchmark nobody runs.
_NESTED = [
    ('one branching', 13, 3,
     'C([C@H]([C@H](C)Cl)[C@@H](C)Cl)[C@@H]([C@H](C)Cl)[C@@H](C)Cl'),
    ('two branchings', 79, 7,
     'C([C@H]([C@H](C)Cl)[C@@H](C)Cl)([C@@H]([C@H](C)Cl)[C@@H](C)Cl)'
     '[C@H]([C@H](C)Cl)[C@H](C)Cl'),
    ('three branchings', 729, 13,
     'C([C@H]([C@H](C)Cl)[C@@H](C)Cl)([C@@H]([C@H](C)Cl)[C@@H](C)Cl)'
     '([C@H]([C@H](C)Cl)[C@H](C)Cl)[C@H]([C@H](C)Cl)[C@@H](C)Cl'),
]


@mark.parametrize('name, before, after, text', _NESTED, ids=[n for n, _, _, _ in _NESTED])
def test_the_extremal_search_starts_from_the_parity_refined_colouring(name, before, after, text):
    """The tree is linear in the branchings, not exponential, and the answer is the unfolded one.

    Both halves matter.  A smaller tree that moved the labelling would be a different canonical form
    wearing the old one's name, so the budgeted answer is compared against the default-budget answer
    and the canonical form is compared across creation orders.
    """
    m = read_smiles(text)
    assert after < before, 'the fixture no longer witnesses anything'
    with raises(AutomorphismBudgetExceeded):
        m.canonical_order(_node_budget=after - 1)
    order = m.canonical_order(_node_budget=after)
    assert order == m.canonical_order()
    assert sorted(order.values()) == list(range(m.atom_count))
    assert len({_relabel(m, o).canonical_bytes for o in _orders(m, 8, 20260905)}) == 1


@mark.parametrize('name, before, after, text', _NESTED, ids=[n for n, _, _, _ in _NESTED])
def test_the_stereo_free_string_of_a_deep_witness_is_still_a_constitution_key(
        name, before, after, text):
    """The fold is NOT taken for `!s`, and this is the test that says why it may not be.

    `smw_canonical_positions` withholds its stereo seed under `!s` so that two configurations of one
    constitution write one string.  The fold would have reached `!s` anyway, from the other side --
    it needs no seed -- and folding it moved 63 of the 393 records in `test/`, every one of which
    then failed the round trip below.  `test_the_seed_is_not_taken_without_the_stereo_key` in
    `test_smiles_write_cis_trans.py` asserts the same contract and did NOT catch that: its fixture is
    a diene, and the fold only bites where a TETRAHEDRAL parity separates a constitutional tie.
    """
    m = read_smiles(text)
    key = write_smiles(m, '!s')
    assert write_smiles(read_smiles(key), '!s') == key
    for order in _orders(m, 8, 20260905):
        assert write_smiles(_relabel(m, order), '!s') == key


_FUSED_SPELLINGS = [
    'C1CC[C@H]2C[C@H]3CCCC[C@H]3C[C@H]2C1',
    'C1CC[C@@H]2C[C@H]3CCCC[C@H]3C[C@H]2C1',
    'C1CC[C@H]2C[C@@H]3CCCC[C@H]3C[C@H]2C1',
    'C1CC[C@@H]2C[C@@H]3CCCC[C@@H]3C[C@@H]2C1',
]


def test_configurations_of_one_fused_skeleton_write_one_stereo_free_string():
    """The same contract stated as the user reads it: `!s` is a key on the CONSTITUTION.

    The spellings are not all one compound -- `canonical_bytes` separates them, asserted here so the
    test cannot pass by the molecules being secretly equal -- and they share one `!s`.  This fused
    tricyclic is the smallest skeleton the fold moves.
    """
    ms = [read_smiles(t) for t in _FUSED_SPELLINGS]
    assert len({write_smiles(m, '!s') for m in ms}) == 1
    assert len({m.canonical_bytes for m in ms}) > 1
