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
"""The product-side post-filter (N4): the half of a product atom that TESTS instead of building.

`D`, `h`, `H`, `x`, `z`, `r`, `R`, `M`-the-metal-test and `@`-the-ring-bond are read against the
PATCHED molecule.  Nothing about that is a compromise: the graph a product primitive describes is the
graph it is read in, which is why a cyclization states its ring size as a product-side `r5` and needs
no ring-size argument, no relational vocabulary and no second kernel.  The reader has already refused
anything that neither builds nor tests, so a template cannot carry a token with no effect.

**A rejection is a negative result, not an error.**  A candidate whose product fails a check is not
yielded and the enumeration continues -- the same shape as a reactant side that does not match, moved
later in the pipeline.  It does get one log line, which the matcher's silence does not, because a
template that matched and was then rejected is otherwise indistinguishable from one that never
matched: the line names the primitive in the author's own notation so they can find it in their
string.

**The bits are the matcher's bits.**  Every check compiles through `prim_apply`, the same primitive
compiler the reactant side's boxes are built with, so three semantics that are easy to reimplement
wrongly come across for free and are pinned here as consequences rather than as separate rules:

* `D`, `x` and `z` exclude a dative bond where `degree_of()` and `heteroatoms_of()` do not;
* an implicit count the patch could not derive is `H_UNKNOWN`, which answers NEITHER `h0` nor `!h0`
  -- an atom whose hydrogen count nobody knows cannot answer a question about it in either
  direction;
* a ring-size demand is a multi-hot test, so `r5,r6` is a disjunction and `r5;r6` would be a
  conjunction over two spans.
"""
from pytest import raises
from .._core import read_smiles as smiles, read_smirks, IncorrectSmirks


# one bond formed, one leaving group lost, and a slot for the check under test.  Deliberately the
# dullest reaction in the file: what is being measured is the filter, not the patch.
SUBSTITUTION = '[C:1][Br;D1]>>[C;%s:1][I;D1:2]'
# 4-bromobutan-1-ol closing to tetrahydrofuran -- N4's motivating case, where the ring the primitive
# asks about is one the patch created and no reactant-side primitive could have seen
CYCLIZATION = '[O;D1;h1:1][C:2][C:3][C:4][C:5][Br;D1]>>[O;%s:1]1[C:2][C:3][C:4][C:5]1'


def outcome(template, molecule):
    """`(product SMILES or None, log messages)` for a single-site template.

    Returns the product string rather than the container because every fixture here is
    stereo-free -- N9's objection to keying on a formatted SMILES is about stereo, and there is
    none in this file.
    """
    log = []
    reactions = list(read_smirks(template)(smiles(molecule), log=log))
    assert len(reactions) < 2, 'the fixtures are single-site by construction'
    messages = [record.message for record in log]
    if not reactions:
        return None, messages
    assert len(reactions[0].products) == 1
    return str(reactions[0].products[0]), messages


# --- the nine primitives ------------------------------------------------------------------------

def test_degree_holds_and_fails():
    assert outcome(SUBSTITUTION % 'D2', 'CCBr')[0] == 'C(C)I'
    assert outcome(SUBSTITUTION % 'D3', 'CCBr')[0] is None


def test_implicit_hydrogens_are_read_after_the_recompute():
    # the reactant carbon has two hydrogens and keeps them, but the count this asks about is the one
    # the patch wrote -- the filter runs after the hydrogen pass for exactly this reason
    assert outcome(SUBSTITUTION % 'h2', 'CCBr')[0] == 'C(C)I'
    assert outcome(SUBSTITUTION % 'h3', 'CCBr')[0] is None


def test_total_hydrogens():
    assert outcome(SUBSTITUTION % 'H2', 'CCBr')[0] == 'C(C)I'
    assert outcome(SUBSTITUTION % 'H1', 'CCBr')[0] is None


def test_heteroatoms():
    assert outcome(SUBSTITUTION % 'x1', 'CCBr')[0] == 'C(C)I'
    assert outcome(SUBSTITUTION % 'x0', 'CCBr')[0] is None


def test_hybridization():
    assert outcome(CYCLIZATION % 'z1', 'OCCCCBr')[0] == 'O1CCCC1'
    assert outcome(CYCLIZATION % 'z2', 'OCCCCBr')[0] is None


def test_ring_size_is_read_on_the_ring_the_patch_created():
    """N4's whole point.  No reactant-side primitive can say this: the ring does not exist yet."""
    assert outcome(CYCLIZATION % 'r5', 'OCCCCBr')[0] == 'O1CCCC1'
    assert outcome(CYCLIZATION % 'r6', 'OCCCCBr')[0] is None


def test_ring_count():
    assert outcome(CYCLIZATION % 'R1', 'OCCCCBr')[0] == 'O1CCCC1'
    assert outcome(CYCLIZATION % 'R0', 'OCCCCBr')[0] is None


def test_the_metal_test_on_an_inherited_element():
    """`M` in the element position is the metal test and a CHECK; `M` as a modifier is the mask.

    Two meanings for one letter, told apart by position, and only the mask is refused on a product
    side (there is nothing to protect where nothing is matched).  A product atom writing `M` states
    no element, so it has to pair and inherit one -- which makes this the narrowest test there is:
    build the element from the reactant, then assert what class it belongs to.
    """
    assert outcome('[Na;D0;*:1]>>[M:1]', '[Na+]')[0] == '[Na]'
    assert outcome('[C:1][Br;D1]>>[M:1][I;D1:2]', 'CCBr')[0] is None
    assert outcome('[C:1][Br;D1]>>[!M:1][I;D1:2]', 'CCBr')[0] == 'C(C)I'


def test_a_ring_bond_is_checked_on_the_bond_and_not_on_its_atoms():
    """The half-edge word, never the atom's aggregate, which ORs every incident bond together."""
    assert outcome(CYCLIZATION.replace('[O;%s:1]1', '[O:1]1-;@'), 'OCCCCBr')[0] == 'O1CCCC1'
    assert outcome(CYCLIZATION.replace('[O;%s:1]1', '[O:1]1-;!@'), 'OCCCCBr')[0] is None


# --- what the compiler gives for free -----------------------------------------------------------

def test_a_check_on_a_created_atom():
    assert outcome('[C:1][Br;D1]>>[C:1][I;D1:2]', 'CCBr')[0] == 'C(C)I'
    assert outcome('[C:1][Br;D1]>>[C:1][I;D2:2]', 'CCBr')[0] is None


def test_an_unknown_hydrogen_count_answers_NEITHER_direction():
    """A xenon the valence collection has no row for gets `H_UNKNOWN`, not a guessed zero.

    So `h0` fails -- and `!h0` fails as well, since a negated demand is satisfied by data and never by
    the absence of it.  Worth checking on a ported template: chython 2 compares `h` against `None`, so
    `None != 0` is True and `!h0` there matches an atom whose count nobody stated.
    """
    assert outcome('[C:1][Br;D1]>>[C:1][Xe;h0:2]', 'CCBr')[0] is None
    assert outcome('[C:1][Br;D1]>>[C:1][Xe;!h0:2]', 'CCBr')[0] is None


def test_an_alternative_is_a_disjunction_and_a_clause_a_conjunction():
    assert outcome(CYCLIZATION % 'r5,r6', 'OCCCCBr')[0] == 'O1CCCC1'
    assert outcome(CYCLIZATION % 'r6,r7', 'OCCCCBr')[0] is None
    assert outcome(CYCLIZATION % 'r5;R1', 'OCCCCBr')[0] == 'O1CCCC1'
    assert outcome(CYCLIZATION % 'r5;R2', 'OCCCCBr')[0] is None


def test_a_negation_holds_where_the_positive_form_does_not():
    assert outcome(CYCLIZATION % '!r6', 'OCCCCBr')[0] == 'O1CCCC1'
    assert outcome(CYCLIZATION % '!r5', 'OCCCCBr')[0] is None


def test_a_check_does_not_bring_the_neutral_charge_default_with_it():
    """The box is built from the check primitives and NOTHING else.

    A query atom's unstated defaults are supplied at seal -- neutral charge above all -- and a check
    box must not have them: the charge is the build half's business, and a filled box would reject
    every product atom the template deliberately charged.
    """
    assert outcome('[C;D1;h3:1][O;D1;h1]>>[C;D1:1][O;D1;-:2]', 'CO')[0] == 'C[O-]'


# --- the report ---------------------------------------------------------------------------------

def test_a_rejection_names_the_primitive_in_the_authors_own_notation():
    product, messages = outcome(CYCLIZATION % 'r6', 'OCCCCBr')
    assert product is None
    assert len(messages) == 1
    assert '`r6`' in messages[0]
    assert 'map number 1' in messages[0]


def test_a_rejection_on_a_bond_names_both_atoms():
    product, messages = outcome(CYCLIZATION.replace('[O;%s:1]1', '[O:1]1-;!@'), 'OCCCCBr')
    assert product is None
    assert len(messages) == 1
    assert '`!@`' in messages[0]
    assert 'between atoms' in messages[0]


def test_one_rejected_candidate_does_not_take_the_others_with_it():
    """4-bromobutan-2-ol's two bromides, one primary and one secondary, against a demand for `h2`."""
    log = []
    template = read_smirks('[C:1][Br;D1]>>[C;h2:1][I;D1:2]')
    reactions = list(template(smiles('CC(Br)CCBr'), log=log))
    assert len(reactions) == 1
    assert str(reactions[0].products[0]) == 'C(CC(Br)C)I'
    assert len(log) == 1
    assert '`h2`' in log[0].message


# --- the refusals the reader keeps making -------------------------------------------------------

def test_the_any_charge_wildcard_is_neither_role_and_says_so():
    """`*` withdraws the reactant side's neutral default.  A patch has no default to withdraw."""
    with raises(IncorrectSmirks) as exc:
        read_smirks('[Na;D0;*:1]>>[M;*:1]')
    assert 'the any-charge wildcard' in str(exc.value)
    assert 'State the charge outright' in str(exc.value)


def test_the_mask_is_still_refused_on_a_product_side():
    with raises(IncorrectSmirks) as exc:
        read_smirks('[Na;D0;*:1]>>[Na;M:1]')
    assert 'protects a MATCHED atom' in str(exc.value)
