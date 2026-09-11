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
"""The patcher: a template applied to molecules, on the container's own edit session.

`read_smirks` works out what a template MEANS; this is what it DOES.  Five groups of tests, and the
last four each pin a place where a reactor can quietly answer for something the template did not say:

**Building.**  Intermolecular, intramolecular, deletion by absence, the deletion closure, and `M`.
The closure is the one piece of V2's reactor carried forward as behaviour rather than as a rewrite:
an unmapped fragment hanging off a deleted atom goes with it, so a template that means to keep the
alkyl of an ester has to map the alkyl.  `M` is how a template says "matched, and not mine to
delete".

**Composition.**  Which molecules come out.  A reaction's reactants are the INPUTS THE MATCH TOUCHED
and nothing else, and its products are every component of the patched graph that holds a touched or
created atom -- which is what keeps a counter-ion from evaporating when the reaction centre is the
anion it pairs with.  Both halves are asserted, because either one alone is satisfiable by a wrong
rule.

**Hydrogens.**  Recomputed for the atoms the patch wrote and for the surviving neighbours of a
deletion, and for nobody else.  An atom that merely sat inside the match keeps its stored count
EXACTLY, wrong or not -- the input is the input.  Where the valence collection has no answer (an
aromatic bond in the reaction centre, no row for the state) the atom gets `H_UNKNOWN` and never a
guessed zero, and `kekule()` plus `chython.chemistry.calc_implicit` is the caller's repair.

**Stereo (N5, N7, N8, N10).**  A single substitution at a configured tetrahedral centre keeps its
configuration, and the patcher does not lift a finger to make that happen: the arena re-bases the
parity positionally on apply, because a parity is stored against CSR ascending-neighbour order and
not against insertion order.  So V2's neighbour-SET comparison -- which cannot see a substitution at
all, both sides having the same count -- is not ported.  What this layer owes is the reporting: when
the arena drops a parity the drop is a log record naming the anchor (N7), the enhanced-stereo group
goes with it in the same operation (N8), and no candidate rule of the patcher's own ever invents a
parity the input did not carry (N10).

**Identity and survival (N9, N11).**  Duplicate embeddings are deduped on STRUCTURE -- never on a
formatted SMILES string, whose canonical form can oscillate on a symmetric stereocentre -- and a
candidate that fails takes itself out of the enumeration rather than the enumeration with it.
"""
from pytest import raises
from .._core import read_smiles as smiles, read_smirks


def products_of(template, *molecules, **kwargs):
    """Every reaction the template gives, as lists of product SMILES.

    A convenience for the build tests only.  Nothing that asserts about stereo goes through it: N9
    says a formatted string is not the identity of a stereo-bearing product, and a test that reads
    one is a test that can pass for the wrong reason.
    """
    return [[str(p) for p in r.products] for r in template(*molecules, **kwargs)]


# --- building -----------------------------------------------------------------------------------

def test_intermolecular_amide():
    """Two inputs, one bond made and one broken.  The plainest thing a template can be."""
    t = read_smirks('[C;z2:1][Cl;D1].[N;D1;h2:2]>>[C:1][N:2]')
    log = []
    reactions = list(t(smiles('CC(=O)Cl'), smiles('CCN'), log=log))
    assert len(reactions) == 1
    r = reactions[0]
    assert [str(m) for m in r.reactants] == ['C(C)(=O)Cl', 'C(C)N']
    assert [str(m) for m in r.products] == ['C(=O)(NCC)C']
    assert log == []


def test_intramolecular_cyclization():
    """One input, a ring closed inside it.

    V2's reactor could not express this: its patcher worked over a list of molecules and a template
    fragment could not reach across two atoms of ONE of them without the caller pre-uniting the
    inputs itself.  Here the sides are ordinary SMARTS and the ring bond is an ordinary product
    bond, so there is nothing to special-case -- the test exists to prove there is nothing.
    """
    t = read_smirks('[C;D1;h3:1][C:2][C:3][C:4][Br;D1]>>[C:1]1[C:2][C:3][C:4]1')
    assert products_of(t, smiles('CCCCBr')) == [['C1CCC1']]


def test_deletion_is_by_absence():
    """A reactant atom the product side does not mention is deleted.  No `:100` convention.

    Ester hydrolysis with the alkyl MAPPED, which is what a template must do when it means to keep
    the fragment -- see the closure test for what happens when it does not.
    """
    t = read_smirks('[C:1](=[O:2])[O:3][C:4]>>[C:1](=[O:2])[O;D1:3].[C:4][O;D1:5]')
    assert products_of(t, smiles('CC(=O)OCC')) == [['C(C)(=O)O', 'C(C)O']]


def test_deletion_closure_takes_the_unmapped_fragment():
    """An unmapped fragment hanging off a deleted atom goes with it.  V2's `_get_deleted`, restated.

    The same hydrolysis with the alkyl UNMAPPED: the ester oxygen is deleted, and the ethyl behind
    it now reaches no surviving matched atom, so it leaves too.  This is V2's behaviour and it is
    kept: a leaving group is usually several atoms and a template naming only its attachment point
    means to lose all of it.

    One deliberate difference from V2: it starts the walk at a neighbour of the deleted atom without
    testing that neighbour for membership in the doomed set, so a chain of two deleted atoms lets the
    walk cross the first one.  Here nothing crosses a deleted atom.
    """
    t = read_smirks('[C:1][O][C]>>[C:1][O;D1:9]')
    assert products_of(t, smiles('CC(=O)OCC')) == [['C(C)(=O)O']]


def test_mask_exempts_an_atom_from_deletion():
    """`M` says "matched, and not mine to delete".  The contrast is the whole test.

    Both templates match ethyl acetate's ester oxygen and the carbon behind it and mention neither
    in the product; the masked one keeps the methyl-and-oxygen intact, the bare one deletes them.
    """
    masked = read_smirks('[C:1][O:2][C;M]>>[C:1][O:2]')
    assert products_of(masked, smiles('CC(=O)OCC')) == [['C(C)OC(C)=O']]

    bare = read_smirks('[C:1][O:2][C]>>[C:1][O:2]')
    assert products_of(bare, smiles('CC(=O)OCC')) == [['C(C)(=O)O']]


def test_created_atom_is_born_and_bonded():
    """A product-side atom with no reactant partner is created, with the properties the string states.

    EXPLICIT-ONLY, as V2 was: an unstated charge on a created atom is zero rather than "whatever
    the neighbour had", which is why the created oxygen here is neutral without saying so.
    """
    t = read_smirks('[C:1][Br;D1]>>[C:1][O;D1;-:2]')
    r = next(iter(t(smiles('CCBr'))))
    p = r.products[0]
    created = [n for n in p.atom_numbers if p.map_number_of(n) == 0 and p.element_of(n) == 8]
    assert len(created) == 1
    assert p.charge_of(created[0]) == -1


# --- composition --------------------------------------------------------------------------------

def test_untouched_input_is_not_a_reactant():
    """An input the match never reached is not part of the reaction at all.

    Not a filter over the products -- the reactants too.  A caller who hands in a whole reaction
    mixture gets back the reaction that happened, not the mixture with an arrow in it.
    """
    t = read_smirks('[N;D1;h2:1]>>[N;D1;h2:1]')
    r = next(iter(t(smiles('CCN'), smiles('c1ccccc1'))))
    assert [str(m) for m in r.reactants] == ['C(C)N']
    assert [str(m) for m in r.products] == ['C(C)N']


def test_counter_ion_survives_as_its_own_product():
    """A component of a TOUCHED input that the match did not reach still comes out.

    The rule is per-input, not per-component: sodium acetate's reaction centre is the anion, and the
    sodium is a separate connected component of the same input that no template atom matched.  It is
    a product, on its own, because a protonation that silently discards the counter-ion is a
    protonation that does not balance -- and because the caller wrote the salt down on purpose.
    """
    t = read_smirks('[C;z2:1](=[O;D1:2])[O;D1;-:3]>>[C:1](=[O:2])[O;D1:3]')
    assert products_of(t, smiles('CC(=O)[O-].[Na+]')) == [['C(C)(=O)O', '[Na+]']]


def test_reactants_are_snapshots_not_the_patched_graph():
    """The reactant side is the input as it arrived, and the input itself is never mutated.

    The patcher works on a copy, and it is worth an assertion because the copy is what makes an
    enumeration over several embeddings independent.
    """
    m = smiles('CCBr')
    t = read_smirks('[C:1][Br;D1]>>[C:1][O;D1:2]')
    r = next(iter(t(m)))
    assert str(r.reactants[0]) == 'C(C)Br'
    assert str(m) == 'C(C)Br'


def test_the_patcher_numbers_its_products_from_one():
    # THE REACTOR IMPOSES ITS MAPPING.  `[C:1][Br;D1]>>[C:1][O;D1;-:2]` both deletes (the bromine
    # pairs with nothing) and creates (the oxide has no reactant partner), so the product carries one
    # atom of each class beside the two carried ones -- and the created one stays 0.
    t = read_smirks('[C:1][Br;D1]>>[C:1][O;D1;-:2]')
    p = next(iter(t(smiles('CCBr')))).products[0]
    assert {n: p.map_number_of(n) for n in p.atom_numbers} == {1: 1, 2: 2, 4: 0}


def test_an_inputs_own_map_numbers_are_replaced_and_not_carried():
    # Imposed, not preserved: preserving cannot be made 1-1 across two inputs, so it is not attempted
    # on one either.  The input's 7/8/9 do not survive.
    t = read_smirks('[C:1][Br;D1]>>[C:1][O;D1;-:2]')
    p = next(iter(t(smiles('[CH3:7][CH2:8][Br:9]')))).products[0]
    assert sorted(p.map_number_of(n) for n in p.atom_numbers) == [0, 1, 2]


def test_call_refuses_a_non_molecule_eagerly():
    """The refusal comes from the CALL, not from the first `next()`.

    A generator that validates its arguments lazily reports a caller's type error from somewhere
    inside a `for` loop, several frames from the mistake.  Both the wrong type and the empty call
    are checked before the generator exists.

    A collection handed in whole gets its own message, because a `MoleculeContainer` iterates over
    atom ids: `template(*molecules)` would otherwise reach the matcher as a list of ints.
    """
    t = read_smirks('[C:1][Br;D1]>>[C:1][O;D1:2]')
    with raises(TypeError):
        t('CCBr')
    with raises(TypeError):
        t(smiles('CCBr'), 'CCBr')
    with raises(ValueError):
        t()
    with raises(TypeError, match=r'unpack the list at the call'):
        t([smiles('CCBr')])


# --- hydrogens ----------------------------------------------------------------------------------

def test_hydrogens_are_recomputed_where_the_patch_wrote():
    """A deletion's surviving neighbour gets its count re-derived.

    Losing the bromine leaves a carbon that had one hydrogen with two, and nothing else in the
    molecule changes.
    """
    t = read_smirks('[C;z1:1][Br;D1]>>[C:1]')
    p = next(iter(t(smiles('CC(N)Br')))).products[0]
    assert p.implicit_h_of(2) == 2
    assert str(p) == 'C(C)N'


def test_matched_but_unwritten_atom_keeps_its_stored_count():
    """An atom that merely sat inside the match is not recomputed, even when its count is wrong.

    `[CH1]CCl` states one hydrogen on a carbon whose valence rules would give three.  The template
    matches that carbon -- it is `:1`, it is paired, it is part of the reaction centre by anybody's
    reading -- and changes nothing about it.  So it comes out with one hydrogen still.

    This is the input-is-garbage rule pointed at the patcher: a reaction is not a repair pass, and
    an atom's stored count is a fact about the input that only an explicit repair may overwrite.
    Recomputing the whole match instead would be easier to write and would silently launder every
    hand-drawn hydrogen count in the reaction centre.
    """
    m = smiles('[CH1]CCl')
    assert m.implicit_h_of(1) == 1
    t = read_smirks('[C:1][C:2][Cl;D1]>>[C:1][C:2][O;D1:3]')
    p = next(iter(t(m))).products[0]
    assert p.implicit_h_of(1) == 1
    assert p.implicit_h_of(2) == 2


def test_an_aromatic_reaction_centre_is_derived_LIKE_ANYTHING_ELSE():
    """THIS TEST USED TO ASSERT `None` HERE, and the `None` was the patcher's own limitation.

    `smk_hydrogens` walked the CSR itself and returned `H_UNKNOWN` on sight of any order-4 bond,
    justified by "the valence collection has no row for an atom holding aromatic bonds".  True, and
    beside the point: the collection is not asked about order 4, it is asked about the bond-order sum
    the atom's aromatic CLASS implies, and `arom_classify_atom` supplies the class.  A mono-substituted
    aromatic carbon must take one ring double bond -- every Kekule form of the ring agrees, which is
    why the answer is available before anybody kekulises -- so its sum is 4 and neutral carbon at sum 4
    has no hydrogens.

    The patcher now delegates to `_hydrogens.pxi`, the one derivation every reader shares, and this
    product's counts are ATOM FOR ATOM the counts the same molecule gets read straight from
    `c1ccccc1O` -- asserted below, because "0 instead of None" on its own could as easily be a
    regression as a fix.
    """
    t = read_smirks('[C;a:1][Br;D1]>>[C;a:1][O;D1:2]')
    p = next(iter(t(smiles('c1ccccc1Br')))).products[0]
    assert p.implicit_h_of(6) == 0, 'the substituted ring carbon takes a ring double bond in every form'
    assert p.implicit_h_of(1) == 1
    assert p.unknown_h_count == 0

    reference = smiles('c1ccccc1O')
    assert ([p.implicit_h_of(n) for n in p.atom_numbers]
            == [reference.implicit_h_of(n) for n in reference.atom_numbers])


def test_the_ambiguous_aromatic_atom_is_STILL_unknown_after_a_patch():
    """The honest half, and the one the derivation cannot close: pyrrole versus pyridine.

    N-demethylating N-methylpyrrole is the case, and it has to be a patch that changes a bond AT the
    nitrogen: only `changed` atoms are recomputed, so substituting a ring carbon leaves the azole
    nitrogen's stated count alone -- correctly, the input said the number.  Take its substituent away
    and the nitrogen is two-coordinate and aromatic with nothing left stating its hydrogens, which is
    the pyrrole-versus-pyridine choice exactly: with a hydrogen it donates its lone pair and takes no
    ring double bond, without one it must take one, and only the RING decides.

    So `H_UNKNOWN` survives here, and it is now the NARROW case rather than every aromatic atom the
    patch touched.  `kekule()` closes it, in one call and not two; see
    `chython/chemistry/test/test_reaction_hydrogen_repair.py`.
    """
    t = read_smirks('[N;a;D3:1]-[C;D1;z1]>>[N;a;D2:1]')
    p = next(iter(t(smiles('Cn1cccc1')))).products[0]
    unknown = [n for n in p.atom_numbers if p.implicit_h_of(n) is None]
    assert [p.element_of(n) for n in unknown] == [7], 'only the pnictogen, and only because of the ring'
    assert p.unknown_h_count == 1


def test_created_atom_hydrogens_are_derived():
    """A created atom is born UNKNOWN and then derived like anything else the patch wrote."""
    t = read_smirks('[C:1][Br;D1]>>[C:1][O;D1:2]')
    p = next(iter(t(smiles('CCBr')))).products[0]
    created = [n for n in p.atom_numbers if p.element_of(n) == 8]
    assert len(created) == 1
    assert p.implicit_h_of(created[0]) == 1
    assert p.unknown_h_count == 0


# --- stereo: N5, N7, N8, N10 --------------------------------------------------------------------

def test_single_substitution_re_bases_the_configuration_for_keep_to_hold():
    """N5.  One arm replaced at a configured centre: the arena re-bases, `@=` keeps what it re-based.

    The arena stores a parity against CSR ascending-neighbour order and re-bases it on every apply,
    so replacing bromine with a freshly-numbered oxygen -- which lands at a DIFFERENT position in
    that order -- is re-based positionally rather than dropped.  The patcher contributes nothing to
    that and must not: a neighbour-SET comparison, V2's guard, answers "unchanged" for exactly this
    case, so it agrees here by accident and disagrees elsewhere.

    What the patcher DOES contribute is the reaction-centre drop, which is why the template says `@=`:
    the substitution happens at the centre, so a template that states nothing gets nothing.  The
    re-base is what this test is about, and `@=` is how it becomes observable -- an implementation that
    dropped the sign instead of re-basing it, or re-based it wrongly, fails on the tag.

    Asserted on canonical bytes rather than on a SMILES string, per N9, and against BOTH tags: an
    equality test that only checks the expected form passes for any implementation that emits a
    constant.
    """
    t = read_smirks('[C;z1:1][Br;D1]>>[C@=:1][O;D1:2]')
    log = []
    p = next(iter(t(smiles('C[C@H](N)Br'), log=log))).products[0]
    assert p.canonical_bytes == smiles('C[C@H](N)O').canonical_bytes
    assert p.canonical_bytes != smiles('C[C@@H](N)O').canonical_bytes
    assert log == []


def test_configuration_away_from_the_centre_is_untouched():
    """N7's negative control.  A reaction that changes nothing at a labelled centre logs nothing.

    The template rewrites a methyl-oxygen bond two bonds away from the stereocentre; the centre
    keeps its stored parity and the log stays empty.  Without this control, N7's requirement is
    satisfiable by a patcher that reports every centre in every molecule it touches.
    """
    t = read_smirks('[C;D1;h3:1][O;D2:2]>>[C:1][O:2]')
    log = []
    p = next(iter(t(smiles('C[C@H](N)OC'), log=log))).products[0]
    assert p.parity_of(2) == 2
    assert log == []


def test_dropped_parity_is_logged():
    """N7.  The drop is a log record and not a silent `else`.

    Removing an arm rather than replacing one leaves the arena nothing to re-base the sign against,
    so the parity goes -- correctly, because three-coordinate carbon with one hydrogen is not a
    stereocentre.  The drop is the right answer; silence about it is not, because a caller cannot then
    tell "this reaction is not stereospecific" from "chython does not know".  The record names the
    anchor and the rule, so both the atom and the template that reached it are recoverable.
    """
    t = read_smirks('[C;z1:1][Br;D1]>>[C:1]')
    log = []
    p = next(iter(t(smiles('C[C@H](N)Br'), log=log))).products[0]
    assert p.parity_of(2) == 0
    assert len(log) == 1
    assert log[0].rule == 'smirks:[C;z1:1][Br;D1]>>[C:1]'
    assert log[0].atoms == (2,)


def test_dropped_parity_clears_its_group():
    """N8.  The parity and the enhanced-stereo group go in one operation, never one without the other.

    A group id left dangling behind a cleared parity is invisible to a caller and rejoins a later pass
    that sets a parity at that atom to a group the atom may not belong to.  Here the group is gone from
    `stereo_groups()`, which is the only place the answer could hide.
    """
    t = read_smirks('[C;z1:1][Br;D1]>>[C:1]')
    m = smiles('C[C@H](N)Br |&1:1|')
    assert m.stereo_groups()
    log = []
    p = next(iter(t(m, log=log))).products[0]
    assert p.parity_of(2) == 0
    assert p.stereo_groups() == {}
    assert len(log) == 1


def test_patcher_invents_no_configuration():
    """N10.  The patcher has no candidate rule of its own, so it cannot label an unlabelled centre.

    The product here IS a tetrahedral centre by every candidate list's reckoning -- four different
    substituents -- and the input said nothing about its configuration.  It comes out saying nothing
    about it: parity zero, no group, empty log.  A candidate list with no symmetry reduction reports
    isopropane's CH and a CF3 carbon, and a reactor that consults one is a reactor that can assert a
    configuration nobody stated.
    """
    t = read_smirks('[C:1][Br;D1]>>[C:1][O;D1:2]')
    log = []
    p = next(iter(t(smiles('CC(N)Br'), log=log))).products[0]
    assert all(p.parity_of(n) == 0 for n in p.atom_numbers)
    assert p.stereo_groups() == {}
    assert log == []


# --- identity and survival: N9, N11 -------------------------------------------------------------

def test_duplicate_embeddings_are_deduped_on_structure():
    """N9.  Two embeddings that give the same product give one reaction.

    Dimethyl ether matches this template both ways round, and both ways give methanol from the same
    atoms.  A dedupe keyed on `str(reaction)`, as V2's is, is a live hazard the moment a template emits
    stereo: canonical output can oscillate on a symmetric stereocentre, so the key can disagree between
    runs.  The key here is the touched inputs plus the products' canonical bytes.

    With the automorphism filter off the count is unchanged, which is the point: the filter reduces
    the SEARCH, and the dedupe is what makes the ANSWER a set.
    """
    t = read_smirks('[C;D1;h3:1][O;D2:2][C;D1;h3]>>[C:1][O;D1:2]')
    assert products_of(t, smiles('COC')) == [['CO']]
    t = read_smirks('[C;D1;h3:1][O;D2:2][C;D1;h3]>>[C:1][O;D1:2]')
    assert len(list(t(smiles('COC'), automorphism_filter=False))) == 1


def test_distinct_sites_are_not_deduped():
    """The other half of N9: two sites that give different molecules are two reactions.

    1,2-dibromopropane has two C-Br bonds and they are not equivalent, so substituting one is not
    substituting the other and both answers are wanted.  A dedupe that keyed on the touched INPUTS
    alone -- which is half of the key -- would collapse them, so the products' canonical bytes are
    the other half.

    Its symmetric cousin is the control below.
    """
    t = read_smirks('[C:1][Br;D1]>>[C:1][O;D1:2]')
    assert products_of(t, smiles('BrCC(Br)C')) == [['C(CO)(C)Br'], ['C(CBr)(C)O']]


def test_equivalent_sites_give_one_answer():
    """1,4-dibromobutane's two bromines ARE equivalent, and the reaction is reported once.

    Both halves of the machinery agree here and that is worth pinning: the automorphism filter never
    offers the second embedding, and if it did the structural key would reject it -- the product is
    the same molecule either way.  With the filter off the count is still one, which says the dedupe
    is doing its own work and not riding on the search.
    """
    t = read_smirks('[C:1][Br;D1]>>[C:1][O;D1:2]')
    assert products_of(t, smiles('BrCCCCBr')) == [['C(CBr)CCO']]
    assert len(list(t(smiles('BrCCCCBr'), automorphism_filter=False))) == 1


def test_a_pathological_template_cannot_kill_the_enumeration():
    """N11.  Every candidate is guarded end to end -- the patch AND the dedupe key.

    A guard around the patcher alone, with the dedupe outside it, lets a bad key escape, propagate out
    of the generator and take the tail of the enumeration with it.  Here the whole candidate is inside
    one guard and a failure becomes a log line and a skipped candidate.

    THE GUARD CANNOT BE TRIPPED BY ANYTHING IN THIS TEST, and that is the finding rather than a gap:
    the arena stores what it is told.  Eighteen fluorines on one carbon, a nonexistent uranium
    isotope and a `+8` carbon are all written and none of them raises -- an unsatisfiable valence
    yields `H_UNKNOWN` and a formal charge at the edge of the storable range is stored.  So what is
    asserted is the invariant the guard is there to protect: each of these enumerates to completion,
    with every candidate accounted for.  A future template language that can state something the
    arena refuses is what will exercise the `except` branch, and it will find it already written.
    """
    cases = [
        # a valence no collection has a row for
        ('[C;D0:1]>>[C:1]' + ''.join('(-[F;D1:%d])' % n for n in range(2, 20)), 'C'),
        # an isotope that does not exist
        ('[C;D0:1]>>[C:1][300U;D1:2]', 'C'),
        # the far end of the storable charge range
        ('[C;D0:1]>>[C;+8:1]', 'C'),
        # every matched atom deleted but one
        ('[C:1][O:2]>>[C:1]', 'CO'),
    ]
    for pattern, molecule in cases:
        log = []
        reactions = list(read_smirks(pattern)(smiles(molecule), log=log))
        assert len(reactions) == 1, pattern
        assert reactions[0].products, pattern
        assert log == [], pattern


def test_the_reactants_carry_the_same_numbers_as_the_products():
    # 1-1 across the arrow, and contiguous from 1 over BOTH sides at once -- the acid takes 1-3 and the
    # alcohol 4-6, so nothing collides.  The acid's leaving OH is on one side only and stays 0.
    t = read_smirks('[C:1]([O;D1;h1:2])=[O:3].[C;z1:4][O;D1;h1:5]>>[C:1](=[O:3])[O:5][C:4]')
    r = next(iter(t(smiles('CC(=O)O'), smiles('CCO'))))
    acid, alcohol = r.reactants
    product = r.products[0]
    assert sorted(acid.map_number_of(n) for n in acid.atom_numbers) == [0, 1, 2, 3]
    assert sorted(alcohol.map_number_of(n) for n in alcohol.atom_numbers) == [4, 5, 6]
    assert sorted(product.map_number_of(n) for n in product.atom_numbers) == [1, 2, 3, 4, 5, 6]
    source = {}
    for molecule in r.reactants:
        for n in molecule.atom_numbers:
            if molecule.map_number_of(n):
                source[molecule.map_number_of(n)] = molecule.element_of(n)
    for n in product.atom_numbers:              # a pair is one atom, so one element
        assert source[product.map_number_of(n)] == product.element_of(n)


def test_the_second_input_is_numbered_at_its_own_atom_numbers():
    # The union renumbers everything after the first, so a reactant snapshot and the products live in
    # two different atom-number spaces and the numbering has to cross back.  Pinned with an alcohol
    # whose own numbers are 11-13, which are neither 1..N nor the 5-7 the union gives it.
    t = read_smirks('[C:1]([O;D1;h1:2])=[O:3].[C;z1:4][O;D1;h1:5]>>[C:1](=[O:3])[O:5][C:4]')
    alcohol = smiles('CCO')
    alcohol.remap({1: 11, 2: 12, 3: 13})
    r = next(iter(t(smiles('CC(=O)O'), alcohol)))
    out = r.reactants[1]
    assert list(out.atom_numbers) == [11, 12, 13]
    assert {n: out.map_number_of(n) for n in out.atom_numbers} == {11: 4, 12: 5, 13: 6}


def test_two_inputs_numbered_from_one_do_not_collide():
    # `union` carries `map_number`, so without the patcher imposing its own numbering two inputs each
    # numbered 1..N arrive in one product with three numbers used twice and a driver cannot tell which
    # atom a number means.  Imposing it is what makes 1-1 guaranteeable rather than merely usual.
    t = read_smirks('[C:1]([O;D1;h1:2])=[O:3].[C;z1:4][O;D1;h1:5]>>[C:1](=[O:3])[O:5][C:4]')
    acid, alcohol = smiles('CC(=O)O'), smiles('CCO')
    for molecule in (acid, alcohol):
        for i, n in enumerate(molecule.atom_numbers, 1):
            molecule.set_map_number(n, i)
    view = next(iter(t(acid, alcohol))).modeling_view()
    assert view.collisions == {'reactants': (), 'products': ()}
    assert view.unmapped == {'reactants': 1, 'products': 0}


def test_a_substrate_past_the_map_number_ceiling_still_reacts():
    # `atom_t.map_number` is 16 bits with a declared ceiling of 9999, so a bigger substrate cannot be
    # mapped 1-1 at all.  The reaction is still produced and the overflow comes back at 0 with a log
    # line -- refusing here would make an imposed mapping cost a peptide-sized input its reaction.
    t = read_smirks('[C;D1;h3:1][C:2]>>[O;D1:3][C:1][C:2]')
    log = []
    r = next(iter(t(smiles('C' * 10050), log=log)))
    p = r.products[0]
    numbers = sorted(n for n in (p.map_number_of(i) for i in p.atom_numbers) if n)
    assert numbers == list(range(1, 10000))                  # capped, and still 1-1 over what it reached
    assert sum(1 for i in p.atom_numbers if p.map_number_of(i) == 0) == 10050 - 9999 + 1
    assert any('ceiling on a map number' in record.message for record in log)
