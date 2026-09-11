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
"""The SMIRKS reader: the arrow, the two sides, and what only a two-sided string can say.

`read_smirks` applies nothing to a molecule -- there is no patcher yet -- so every test here is about
what the READER worked out.  Three kinds:

**Structure.**  Which atoms pair by map number, which are deleted, which are created.  Deletion is by
ABSENCE, which is the whole of chython 2's `:100` / `:200` leaving-group convention replaced by
nothing at all.

**Refusals.**  The three-part `reactants>agents>products` form, a single `>`, more than one arrow,
whitespace inside a side, a map number naming two atoms of one side.  A message about one side carries
that side's own byte offsets, and the tests assert them: the two sides go through one lexer, and one
lexer cannot have a shared offset space.

**The log.**  The product side is explicit-only -- an unstated charge is zero, not "whatever it
matched" -- and that is silent by construction.  The mapped-pair lint is what makes it visible, and it
is a log line rather than a refusal because neutralizing a cation is a legitimate thing to mean.

**The classification.**  Every product primitive builds something or checks something, and one that
does neither is refused rather than ignored, which is what keeps a dead product-side form out of a
corpus.

Templates here name public compounds and public reactions only.
"""
from pytest import mark, raises
from chython.core import IncorrectSmarts, IncorrectSmirks, ReactionTemplate, read_smirks


# The primitive kinds the classification reports, spelled here rather than left as bare integers in
# the assertions.  They are `_query_boxes.pxi`'s own PRIM_* / BPRIM_* constants, which are internal to
# the extension -- a test asserting on them is asserting on the core's numbering on purpose, because
# that numbering is what the patcher will switch on.
ELEMENT, ANY, METAL, ISOTOPE, CHARGE, RADICAL = 1, 2, 3, 4, 5, 6
DEGREE, IMPLICIT_H, TOTAL_H, HETEROATOMS, HYBRIDIZATION, RING_SIZE, RING_COUNT = 7, 8, 9, 10, 11, 12, 13
STEREO, NO_ISOTOPE = 14, 15
BOND_ORDER, BOND_AROMATIC, BOND_RING = 20, 21, 22


# ----------------------------------------------------------------------------------------------
# THE ARROW AND THE SIDES
# ----------------------------------------------------------------------------------------------
def test_the_two_sides_are_read_into_two_different_kinds_of_object():
    """The structural change the notation exists for: the reactant side is a query and the product
    side is not.  A sealed `QueryContainer` on one side, a patch spec on the other.

    With both sides `QueryContainer`s, as in chython 2, a stereo override guarded by
    `isinstance(ra, Element)` cannot execute -- a `QueryElement` is not an `Element`.  The product side
    is not a query here, so there is nothing to route around."""
    t = read_smirks('[C;D1;x1:1][Br;D1]>>[C;D1:1][O;D1;h1]')
    assert isinstance(t, ReactionTemplate)
    assert t.reactants.atom_count == 2
    assert t.product_atom_count == 2
    assert t.product_bonds == ((1, 2),)
    # the reactant side seals; nothing on the template hands out a sealed product side, because a
    # patch has no boxes to compile
    assert t.reactants.atom_count_sealed() == 2


def test_whitespace_may_surround_the_arrow():
    """`A >> B` is how a human writes it, so the reader takes it.  It cannot use `read_smarts`'s rule
    that the string ends at the first space, because a SMIRKS is not one token."""
    spaced = read_smirks('[C;D1;x1:1][Br;D1] >> [C;D1:1][O;D1;h1]')
    tight = read_smirks('[C;D1;x1:1][Br;D1]>>[C;D1:1][O;D1;h1]')
    assert spaced.mapped_pairs == tight.mapped_pairs
    assert spaced.deleted_atoms == tight.deleted_atoms


def test_the_template_remembers_the_string_it_came_from():
    t = read_smirks(' [C:1]>>[C:1] ')
    assert t.smirks == '[C:1]>>[C:1]'
    assert repr(t) == "read_smirks('[C:1]>>[C:1]')"


def test_a_template_cannot_be_built_any_other_way():
    """Owner decision: one construction entry point.  A second one would be a second notation."""
    with raises(TypeError, match='read_smirks'):
        ReactionTemplate()


def test_bytes_read_the_same_as_str():
    assert read_smirks(b'[C:1]>>[C:1]').smirks == read_smirks('[C:1]>>[C:1]').smirks


def test_a_non_string_is_a_type_error():
    with raises(TypeError, match='str or bytes'):
        read_smirks(42)


# ----------------------------------------------------------------------------------------------
# MAP NUMBERS: PAIRING, DELETION BY ABSENCE, CREATION
# ----------------------------------------------------------------------------------------------
def test_a_map_number_on_both_sides_pairs_the_two_atoms():
    """Finkelstein: the carbon carries through, the bromide leaves, the iodide arrives."""
    t = read_smirks('[C;D1;x1:1][Br;D1:2]>>[C;D1:1][I;D1]')
    assert t.reactant_map_numbers == {1: 1, 2: 2}
    assert t.product_map_numbers == {1: 1}
    assert t.mapped_pairs == {1: (1, 1)}


def test_a_reactant_map_number_absent_from_the_product_side_is_a_deletion():
    """The standard SMIRKS rule, and it is enough on its own.  chython 2 needed the `:100` / `:200`
    numbering convention to tell its patcher which matched atoms to remove, because its two sides were
    separate objects with no way to express absence.  Ported templates may keep those numbers as
    documentation; nothing reads them."""
    t = read_smirks('[C;D1;x1:1][Br;D1:100]>>[C;D1:1][O;D1;h1]')
    assert t.deleted_atoms == {2}
    assert 100 not in t.mapped_pairs


def test_an_unmapped_reactant_atom_is_deleted_too():
    """It pairs with nothing, so it falls to the same rule as a number the product side dropped."""
    t = read_smirks('[C;D1;x1:1][Br;D1]>>[C;D1:1][O;D1;h1]')
    assert t.deleted_atoms == {2}


def test_a_masked_reactant_atom_is_never_deleted():
    """`M` survives from chython 2 with its meaning intact: an atom the pattern must match and must
    never remove, for a template naming something purely as context."""
    plain = read_smirks('[C:1][N;D1]>>[C:1]')
    assert plain.deleted_atoms == {2}
    masked = read_smirks('[C:1][N;D1;M]>>[C:1]')
    assert masked.deleted_atoms == frozenset()
    assert masked.reactants.masked_atoms() == {2}


def test_a_product_atom_with_no_reactant_partner_is_created():
    t = read_smirks('[C;D1;x1:1][Br;D1]>>[C;D1:1][O;D1;h1]')
    assert t.created_atoms == {2}


def test_a_product_only_map_number_pairs_with_nothing_and_says_so():
    """A number written only on the product side names an atom that is still created; the number
    itself does nothing, and a template author who thought otherwise gets told."""
    log = []
    t = read_smirks('[C:1][Br;D1]>>[C:1][O;D1;h1:7]', log)
    assert t.created_atoms == {2}
    assert t.mapped_pairs == {1: (1, 1)}
    assert [str(x) for x in log] == ['map number 7 is on the product side only, so it pairs with nothing; that atom '
                                     'is created']


def test_a_map_number_on_two_atoms_of_one_side_is_refused():
    """Ambiguity, not chemistry: there is no answer to "which of these two pairs".  Both sides get the
    check, and the message names which one."""
    with raises(IncorrectSmirks, match='map number 1 is on two atoms of the reactant side'):
        read_smirks('[C:1][C:1]>>[C:1]')
    with raises(IncorrectSmirks, match='map number 1 is on two atoms of the product side'):
        read_smirks('[C:1]>>[C:1][C:1]')


# ----------------------------------------------------------------------------------------------
# COMPONENT GROUPING ACROSS THE ARROW (N-G1, N-G2)
# ----------------------------------------------------------------------------------------------
def test_the_reactant_side_carries_its_component_groups_through():
    """The intramolecular case the whole notation change was asked for.  Grouping is the SMARTS
    reader's operator, so a SMIRKS gets it for free -- and chython 2's pattern-fusing constructor,
    which had no notion of groups at all, does not come across."""
    intra = read_smirks('([C;D1;x1:1][Br;D1].[O;D1;h1:2])>>([C;D1:1][O:2])')
    assert intra.reactants.component_groups() == ((frozenset({1, 2}), 0), (frozenset({3}), 0))

    inter = read_smirks('([C;D1;x1:1][Br;D1]).([O;D1;h1:2])>>([C;D1:1][O:2])')
    assert inter.reactants.component_groups() == ((frozenset({1, 2}), 0), (frozenset({3}), 1))

    loose = read_smirks('[C;D1;x1:1][Br;D1].[O;D1;h1:2]>>[C;D1:1][O:2]')
    assert loose.reactants.component_groups() == ((frozenset({1, 2}), None), (frozenset({3}), None))


# ----------------------------------------------------------------------------------------------
# REFUSALS
# ----------------------------------------------------------------------------------------------
def test_the_three_part_form_is_refused_and_says_why():
    """N-G7.  Not an oversight: an agent is matched and never patched, so it is a third semantics for
    the same lexer.  Refusing is reversible; a half-implemented agent side is not."""
    with raises(IncorrectSmirks, match='matched and never patched'):
        read_smirks('[C:1][Br;D1]>[K+].[I-]>[C:1][I]')


@mark.parametrize('pattern,fragment', [
    ('[C:1][C:1]', 'no `>>`'),
    ('[C:1]>[C:1]', 'a single `>` at position 5 is not the SMIRKS arrow'),
    ('[C:1]>>[C:1]>>[C:1]', '4 `>` characters'),
    ('[C:1]>>[C:1] [C:1]', 'the product side holds whitespace at position 5'),
    ('[C:1] [C:1]>>[C:1]', 'the reactant side holds whitespace at position 5'),
])
def test_refusals(pattern, fragment):
    with raises(IncorrectSmirks) as e:
        read_smirks(pattern)
    assert fragment in str(e.value)


def test_a_side_that_is_not_a_smarts_names_the_side_and_keeps_its_own_offsets():
    """Two sides through one lexer buys the shared dialect and costs a shared offset space.  So the
    message says which side, quotes it, and the position inside it is the lexer's own."""
    with raises(IncorrectSmirks) as e:
        read_smirks('[C:1]>>[C:1]1')
    assert 'the product side `[C:1]1` is not a readable SMARTS' in str(e.value)
    assert 'ring bond 1 opens at position 5 and never closes' in str(e.value)


def test_an_empty_side_is_refused_through_the_lexer():
    with raises(IncorrectSmirks, match='no atoms in the string'):
        read_smirks('>>[C:1]')


def test_it_is_catchable_as_a_smarts_error():
    """`IncorrectSmirks` subclasses `IncorrectSmarts`, which subclasses `IncorrectSmiles`, so a
    pipeline reading several notations catches one exception and not three."""
    with raises(IncorrectSmarts):
        read_smirks('[C:1]')


def test_a_non_ascii_string_is_refused():
    with raises(IncorrectSmirks, match='non-ASCII'):
        read_smirks('[C:1]>>[C:1]—')


def test_a_contradictory_reactant_side_is_refused_at_read_time():
    """The reactant side seals, so a term that can never match is an error here rather than a silent
    failure to match at the first use.  The product side does not seal -- it is a patch."""
    with raises(IncorrectSmirks, match='the reactant side does not compile'):
        read_smirks('[C;D1;D2:1]>>[C:1]')


# ----------------------------------------------------------------------------------------------
# THE MAPPED-PAIR LINT (N1, N-G3).  The price of explicit-only product semantics.
# ----------------------------------------------------------------------------------------------
def test_the_lint_fires_when_the_product_side_drops_a_charge():
    """N-G3.  A log line and NOT a refusal: neutralizing a cation is a legitimate thing for a template
    to mean, and the reader's job is to make the silent case visible rather than to decide it."""
    log = []
    read_smirks('[N;+:1][C;D1]>>[N:1]', log)
    assert [str(x) for x in log] == ['map number 1 states a charge on the reactant side and none on the product side, '
                                     'so the product atom is neutral']


def test_the_lint_is_silent_when_the_product_side_restates_the_charge():
    """The negative control N-G3 asks for.  Quaternization of an amine, charge stated on both sides."""
    log = []
    read_smirks('[N;+:1][C;D1:2]>>[N;+:1][C;D1:2]', log)
    assert log == []


def test_the_lint_covers_isotope_and_radical_too():
    log = []
    read_smirks('[13C:1][Br;D1]>>[C:1][O;D1;h1]', log)
    assert [str(x) for x in log] == ['map number 1 states an isotope on the reactant side and none on the product '
                                     'side, so the product atom has no mass number']


def test_a_negated_property_is_not_a_statement_the_lint_reports():
    """`[C;!+]` says something about charge, but what it says is satisfied by the neutral atom an
    unstated product charge produces.  There is no surprise in it, so there is no line."""
    log = []
    read_smirks('[C;!+:1][Br;D1]>>[C:1][O;D1;h1]', log)
    assert log == []


def test_the_lint_only_looks_at_paired_atoms():
    """An atom being deleted has no product side to disagree with, and a created atom has no reactant
    side.  Neither can produce a line."""
    log = []
    read_smirks('[C:1][N;+;D4]>>[C:1][O;-;D1]', log)
    assert log == []


# ----------------------------------------------------------------------------------------------
# THE EXTENSION TAIL.  One tail, reactant atoms then product atoms (N-G10).
# ----------------------------------------------------------------------------------------------
def test_a_product_side_stereo_group_is_recorded_as_a_directive():
    """`racemize` is spelled as an AND group on the product side, which is a positive statement about
    a mixture and not a dropped label.  This is the routing test -- what the field MEANS is
    `test_smirks_stereo.py`'s -- and the index space is the reaction's: two reactant atoms, so the
    product carbon is atom 2.

    No sign on the product atom, and that is the point of the pair: a group needs no frame, so unlike
    a sign it says nothing about the anchor's directions and is accepted on an atom whose degree the
    string does not state."""
    log = []
    t = read_smirks('[C;@:1][Br;D1]>>[C:1][O;D1;h1] |&1:2|', log)
    assert t.product_stereo_groups == {1: (3, 1)}       # (STEREO_AND, group 1)
    assert log == []


def test_the_reactant_side_still_declines_a_stereo_group(capsys):
    """N-G10.  A query cannot test a stereo group, which is what `read_smarts` says about the same
    field; only the product half of the tail is a directive.  One field may name both sides, and each
    index is routed on its own."""
    log = []
    t = read_smirks('[C;@:1][Br;D1]>>[C:1][O;D1;h1] |&1:0,2|', log)
    assert t.product_stereo_groups == {1: (3, 1)}
    assert [str(x) for x in log] == ['the extension field &1:0,2 names 1 reactant atom(s); a query cannot test a '
                                     'stereo group, so that part was not applied']
    assert capsys.readouterr().out == ''


def test_a_radical_field_applies_to_the_reactant_side_and_records_on_the_product_side():
    """The radical is the one field the reactant side CAN test, so it is applied there exactly as
    `read_smarts` applies it -- and then the mapped-pair lint notices that the product side says
    nothing about it, which is the two halves of this reader agreeing."""
    log = []
    t = read_smirks('[C:1]>>[C:1] |^1:0|', log)
    assert t.product_radicals == frozenset()
    assert [str(x) for x in log] == ['map number 1 states a radical on the reactant side and none on the product '
                                     'side, so the product atom is not a radical']

    log = []
    t = read_smirks('[C:1]>>[C:1] |^1:1|', log)
    assert t.product_radicals == {1}
    assert log == []


def test_an_index_past_the_reaction_is_dropped_with_a_count_of_both_sides():
    log = []
    read_smirks('[C:1]>>[C:1] |&1:9|', log)
    assert [str(x) for x in log] == ['the extension field names atom 9, but the reaction has 2 atom(s) (1 reactant, '
                                     '1 product); the mark was dropped']


def test_an_unterminated_tail_is_reported_and_ignored():
    log = []
    t = read_smirks('[C:1]>>[C:1] |&1:1', log)
    assert t.product_stereo_groups == {}
    assert [str(x) for x in log] == ['the extension block after the SMIRKS is not terminated and was ignored: |&1:1']


def test_a_field_this_reader_cannot_apply_is_named():
    log = []
    read_smirks('[C:1]>>[C:1] |c:0|', log)
    assert [str(x) for x in log] == ['the extension field c:0 says nothing this reader can apply to either side']


def test_omitting_the_log_discards_it_rather_than_failing():
    """Every line above is a line a caller may not want; none of them is load-bearing."""
    assert read_smirks('[N;+:1]>>[N:1]').mapped_pairs == {1: (1, 1)}


# ----------------------------------------------------------------------------------------------
# PRODUCT-PRIMITIVE CLASSIFICATION (N4, N-G6).  Build, check, or refused -- never ignored.
# ----------------------------------------------------------------------------------------------
def test_a_product_atom_is_split_into_what_it_builds_and_what_it_checks():
    """The same bracket says both kinds of thing, and the reader separates them once, here, so
    neither the patcher nor the post-filter has to know the other exists."""
    t = read_smirks('[C;D1;x1:1][Br;D1]>>[C;D1:1][O;D1;h1]')
    assert t.product_atom_build == {1: ((ELEMENT, 6),), 2: ((ELEMENT, 8),)}
    assert t.product_atom_check == {1: ((((DEGREE, 1, 0),),),),
                                    2: ((((DEGREE, 1, 0),),), (((IMPLICIT_H, 1, 0),),))}


def test_every_product_atom_and_bond_has_an_entry_in_both_dictionaries():
    """An empty tuple and a missing key are the same fact, and one of the two spellings makes the
    patcher write `.get(...)` at every site.  So every key is present."""
    t = read_smirks('[C:1][O:2]>>[C:1][O:2]')
    assert t.product_atom_build == {1: ((ELEMENT, 6),), 2: ((ELEMENT, 8),)}
    assert t.product_atom_check == {1: (), 2: ()}
    assert t.product_bond_build == {(1, 2): ()}
    assert t.product_bond_check == {(1, 2): ()}


def test_an_untokenised_product_bond_states_nothing_and_a_written_one_states_its_order():
    """A bond with no expression is a single, and the patcher supplies that from an empty build list
    exactly as `query_seal` supplies it for a query.  Two lowercase atoms are the one case where the
    lexer writes the order itself -- aromatic, for the reason it is in SMILES."""
    assert read_smirks('[C:1][O:2]>>[C:1]=[O:2]').product_bond_build == {(1, 2): ((BOND_ORDER, 2),)}
    assert read_smirks('[C:1][C:2]>>[c:1][c:2]').product_bond_build == {(1, 2): ((BOND_AROMATIC, 0),)}
    assert read_smirks('[C:1][C:2]>>[c:1]:[c:2]').product_bond_build == {(1, 2): ((BOND_AROMATIC, 0),)}


def test_a_ring_bond_on_a_product_is_a_check_and_not_a_thing_to_build():
    """`@` on a bond asks whether the result is cyclic.  The patcher cannot make a bond cyclic by
    setting a flag -- the ring is a consequence of the atoms it joined -- so it post-filters."""
    t = read_smirks('[C:1][O:2]>>[C:1]-;@[O:2]')
    assert t.product_bond_build == {(1, 2): ((BOND_ORDER, 1),)}
    assert t.product_bond_check == {(1, 2): ((((BOND_RING, 0, 0),),),)}


def test_a_product_side_r_is_how_a_cyclization_states_its_ring_size():
    """Owner decision: the ring size of a cyclization is a product-side `r` primitive, not an
    argument to anything.  It classifies as a check, which is what makes that decision implementable:
    the patcher closes the ring and the post-filter rejects the sizes the template did not mean."""
    t = read_smirks('([C;D1;x1:1][Br;D1].[O;D1;h1:2])>>([C;r5:1][O;r5:2])')
    assert t.product_atom_check == {1: ((((RING_SIZE, 5, 0),),),), 2: ((((RING_SIZE, 5, 0),),),)}
    assert t.product_atom_build == {1: ((ELEMENT, 6),), 2: ((ELEMENT, 8),)}


def test_a_negated_check_is_still_a_check():
    """Negation only defeats a BUILD primitive.  "not in a five-ring" is a perfectly good question to
    ask of a result, and the negation travels with the primitive for the post-filter to apply."""
    t = read_smirks('[C:1]>>[C;!r5:1]')
    assert t.product_atom_check == {1: ((((RING_SIZE, 5, 1),),),)}


def test_a_disjunction_of_checks_keeps_its_shape():
    """`,` between two checks is answerable -- either alternative satisfies the clause -- so the
    clause survives as a clause rather than being flattened into a conjunction."""
    t = read_smirks('[C:1]>>[C;r5,r6:1]')
    assert t.product_atom_check == {1: ((((RING_SIZE, 5, 0),), ((RING_SIZE, 6, 0),)),)}


def test_the_metal_test_is_a_check_and_the_element_still_has_to_come_from_somewhere():
    """`[M]` as the first primitive is the metal test, which reads perfectly well on a product: the
    result must turn out to be a metal.  It names no element, though, so the atom needs a partner to
    inherit one from -- and a metal salt template has one."""
    t = read_smirks('[M:1][O;D1;h1:2]>>[M:1][O;-:2]')
    assert t.product_atom_check == {1: ((((METAL, 0, 0),),),), 2: ()}
    assert t.product_inherited_elements == {1}


def test_an_unstated_product_element_is_inherited_from_the_partner():
    """The element is the one field with no default: there is no neutral element the way there is a
    neutral charge.  So `[C:1]` states carbon, `[A:1]` says "whatever it matched" out loud, and a
    bracket with neither means the same as `[A:1]`."""
    assert read_smirks('[C:1]>>[A:1]').product_inherited_elements == {1}
    assert read_smirks('[C:1]>>[C:1]').product_inherited_elements == frozenset()
    # transmutation: the product states a different element from the reactant, which is a build
    t = read_smirks('[C:1][Br;D1]>>[Si:1][Cl;D1]')
    assert t.product_atom_build == {1: ((ELEMENT, 14),), 2: ((ELEMENT, 17),)}
    assert t.product_inherited_elements == frozenset()


@mark.parametrize('pattern,fragment', [
    # a build primitive offered as one of several alternatives: nothing to build
    ('[C:1]>>[C,N:1]', 'offers an element as one of several alternatives'),
    ('[C:1]>>[C;+,+2:1]', 'offers a charge as one of several alternatives'),
    # `~` is spelled as an OR of five orders, so it falls to the same rule -- and it SHOULD, because
    # "any bond" is not a bond a patcher can make
    ('[C:1]>>[C:1]~[O;D1;h1]', 'the product bond between atoms 1 and 2 offers a bond order as one '
                               'of several alternatives'),
    ('[C:1]>>[C:1]-,=[O;D1]', 'offers a bond order as one of several alternatives'),
    # a negated build primitive.  Not every build primitive can even be written negated -- the lexer
    # already refuses `!@` on its own ground, that a configuration has an opposite to write instead --
    # so these are the two spellings that reach this layer.
    ('[C:1]>>[C;!+:1]', 'negates a charge; a patch states what to build'),
    ('[C:1]>>[C;!C:1]', 'negates an element'),
    # one field, two statements, in either spelling.  The isotope has no second spelling to collide
    # with: it is a prefix, and a bracket has one prefix position.
    ('[C:1]>>[C;+;+2:1]', 'states the charge twice'),
    ('[C:1]>>[C;+&+2:1]', 'states the charge twice'),
    ('[C:1]>>[C:1]=;#[O]', 'states the bond order twice'),
    # `M` in the masked sense
    ('[C:1]>>[C;M:1]', 'protects a MATCHED atom from deletion'),
    # no element and no partner to inherit one from
    ('[C:1]>>[C:1][D1]', 'product atom 2 states no element and pairs with no reactant atom'),
    ('[C:1]>>[C:1][A]', 'product atom 2 writes `A`'),
    ('[C:1]>>[C:1][A:7]', 'product atom 2 (map number 7) writes `A`'),
])
def test_a_product_primitive_that_cannot_be_placed_is_refused_at_read_time(pattern, fragment):
    """N-G6.  Each of these parses as a SMARTS and would be silently ignored by a patcher that trusted
    its input.  The reader refuses instead, so the dead form never reaches a corpus."""
    with raises(IncorrectSmirks) as e:
        read_smirks(pattern)
    assert fragment in str(e.value)


@mark.parametrize('pattern', [
    '[C;D1;x1:1][Br;D1]>>[C;D1:1][I;D1]',                     # Finkelstein
    '[N;D1;z1:1][C;z1:2]>>[N;+;D2:1][C:2]',                    # quaternization, charge on the product
    '[C;z2:1]=[C;z2:2]>>[C;z1:1][C;z1:2]',                     # hydrogenation
    '[c:1]:[c:2]>>[c:1]:[c:2]',                                # aromatic bonds carried through
    '[C:1][O;D1;h1:2]>>[C:1][O;-:2]',                          # deprotonation
    '[13C:1][Br;D1]>>[13C:1][O;D1;h1]',                        # isotope restated on the product
    # inversion, and the whole of its spelling: `@~` is the unit's other state, so no arm of the anchor
    # has to be named and the reactant side need not sign anything
    '[C:1][Br;D1]>>[C@~:1][I;D1]',
    '[C;@:1][Br;D1]>>[C@~:1][I;D1]',                           # and the same, narrowed to a configured one
    '[C:1][Br;D1]>>[C;&1:1][I;D1]',                            # racemisation, in the bracket
    '[C:1][Br;D1]>>[C;o2:1][I;D1]',                            # and the OR kind of the same thing
    '([C;D1;x1:1][Br;D1].[O;D1;h1:2])>>([C:1][O:2])',          # intramolecular etherification
    '[C:1][Br;D1]>>[C:1][O;D1;h1] |^1:0|',                     # a tail on the reactant side
])
def test_every_product_form_a_template_needs_is_accepted(pattern):
    """The negative control N-G6 asks for.  A classification that refuses too much is as bad as one
    that refuses nothing, and these are the shapes the corpus is made of."""
    assert isinstance(read_smirks(pattern), ReactionTemplate)
