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
"""What a product side can say about a configuration, and what a reactant side can.

Everything a template says about configuration it says in SMARTS, and there is no keyword argument
anywhere that says it instead.  Every statement is on the PRODUCT side:

| the product side says                     | what happens                                                    |
|-------------------------------------------|-----------------------------------------------------------------|
| nothing, at the reaction centre           | dropped -- the template did not say what the reaction did to it |
| nothing, anywhere else                    | whatever the arena re-based is carried                          |
| `@=`                                      | the configuration comes through unchanged, whatever the kind    |
| `@~`                                      | the unit's OTHER state, whatever the kind                       |
| `&<n>` / `o<n>`, no sign                  | racemised: the centre is CONFIGURED and grouped                 |
| a sign in a group with other signed atoms | a DRAWN member of a correlated set                              |
| a sign with none of the above             | refused at `read_smirks` -- nothing to be relative to           |
| `/` and `\\` on both ends of a double bond | the geometry, stated OUTRIGHT: E or Z, whatever came in          |
| `/` or `\\` on one end only                | refused at `read_smirks` -- half a geometry is not one           |

So the SN2 inversion a template exists for is one string with no arms named at all and no sign on the
reactant side:

    [C;z1:1][Br;D1] >> [C@~:1][I;D1:2]

**`@=` AND `@~` ARE RELATIVE AND TAKE NO FRAME.**  The arena re-based the parity into the molecule's own
frame already, and a flip of a value in a frame is a flip in every frame.  **A REACTANT-SIDE SIGN IS
MATCHING SELECTIVITY AND NOTHING ELSE** -- it narrows the template to a substrate that arrives
configured, and no product statement reads it.

**ONE TOKEN FOR EVERY KIND.**  All four units chython models -- tetrahedral, cis/trans, allene,
atropisomer -- hold exactly two states, so "the other one" is well defined for each and `@~` covers a
geometry as readily as a parity.  Which also means the token addresses whatever unit the ATOM anchors: a
template meaning a tetrahedral centre and nothing else narrows its reactant side (`z1`, `D3`) rather
than relying on the token.

**THE REACTION CENTRE'S CONFIGURATION IS DROPPED BY DEFAULT.**  A unit the patch wrote part of -- the
anchor of a tetrahedral centre, either terminal of a double bond -- loses its configuration unless the
template stated one of the things above.  The arena would have carried it: replacing one arm at a centre
keeps the positional frame, so the re-based sign is readable and says "retained", which is a claim about
the course of the reaction that the template never made.  A unit the patch did not touch is not the
reaction centre and nothing here reaches it.

**ABSOLUTE CONFIGURATION HAS NO SPELLING**, and that is chemistry rather than economy.  A configuration
cannot appear where nothing chiral acted, and in a substrate-controlled diastereoselection an absolute
product sign is actively wrong -- it would turn the enantiomeric substrate into the same absolute
product.  A reaction that creates a centre either racemises it (`&<n>`) or states a configuration
relative to one the substrate already had.

**A GEOMETRY IS THE ONE EXCEPTION, and it is not one.**  `/` and `\\` state a cis/trans geometry
outright, because an alkene's two faces are not enantiomeric: the enantiomeric substrate does not give
the enantiomeric product, so nothing chiral has to have acted for a Wittig to make its alkene E.  On a
reactant side the same pair is the matching selectivity a sign is -- it narrows the template to an alkene
that arrives with the stated geometry -- so the two sides read it with different code and refuse the same
malformed strings.

**Relative configuration is that second thing, and it is a correlated group** rather than a new token:
an AND group already means "as drawn, or all members flipped", which is a fixed relative configuration
presented as a racemate.  So a diastereoselective template draws its centres' signs in one group and the
substrate decides the rest -- resolved against a member that arrives configured, and emitted as a real
group when none does.  One template, and the answer follows the substrate.

**A group is a directive too, and only where stated.**  The product side is explicit-only for enhanced
stereo exactly as it is for every other property: a template that states a group gets that group,
renumbered to an id free in the patched molecule, and a template that states none leaves whatever the
input carried alone.  A racemic centre is CONFIGURED and grouped -- never unconfigured -- so a group on a
centre the patch left unconfigured writes a parity too.

**Refusals are at read time.**  A sign with nothing to be relative to, or whose anchor cannot be a
tetrahedral centre in the product, is not a template that misbehaves at application time: `read_smirks`
refuses the string.  The failures that survive to application are the honest ones -- the template is well
formed and the MOLECULE cannot hold what it asks for -- and those are log records, never refusals,
because the input is the input.
"""
from pytest import mark, raises
from .._core import read_smiles as smiles, read_smarts, read_smirks, IncorrectSmarts, IncorrectSmirks


# `C[C@H](N)Br` is 1-bromoethan-1-amine, a configured centre with three heavy neighbours and one
# implicit hydrogen.  `SUBSTITUTION` names none of its arms on purpose: neither `@=` nor `@~` needs a
# frame, so the canonical retention and inversion templates are this short.
SUBSTITUTION = '[C:1][Br;D1]'
SELECTIVE = '[C;@:1][Br;D1]'
SUBSTRATE = 'C[C@H](N)Br'


def product_of(template, molecule, **kwargs):
    """The single product of the single reaction, as a container.

    A container and not a string: N9's rule is that a formatted SMILES is not the identity of a
    stereo-bearing product, so everything below compares `canonical_bytes` or reads `parity_of`.
    """
    reactions = list(read_smirks(template)(smiles(molecule), **kwargs))
    assert len(reactions) == 1, 'the fixtures are single-site by construction'
    assert len(reactions[0].products) == 1
    return reactions[0].products[0]


# --- `@=` and `@~`, the two relative statements ---------------------------------------------------

def test_keep_retains_and_asks_nothing_of_the_reactant_side():
    """`@=` states the course, and an UNSTATED product side does not.

    Both are asserted in one test because this is the reaction centre: retention is a statement, and a
    template that makes none gets no configuration rather than the one the arena could still read.
    Nothing on the reactant side has to be signed, so the same template also serves an unconfigured
    substrate -- which is the second half of the test.
    """
    stated = product_of(SUBSTITUTION + '>>[C@=:1][I;D1:2]', SUBSTRATE)
    unstated = product_of(SUBSTITUTION + '>>[C:1][I;D1:2]', SUBSTRATE)

    assert stated.canonical_bytes == smiles('C[C@H](N)I').canonical_bytes
    assert unstated.canonical_bytes == smiles('C[CH](N)I').canonical_bytes
    assert stated.canonical_bytes != unstated.canonical_bytes

    assert product_of(SUBSTITUTION + '>>[C@=:1][I;D1:2]',
                      'C[CH](N)Br').canonical_bytes == smiles('C[CH](N)I').canonical_bytes


def test_invert_gives_the_other_configuration():
    """The Walden inversion, and the whole of the notation for it: one token, and no arms.

    The reactant side names one direction -- the leaving bromide -- and says nothing about the other
    three.  That is the template the design exists to make writable: the reaction does not care what the
    other three are, so the string does not say.
    """
    product = product_of(SUBSTITUTION + '>>[C@~:1][I;D1:2]', SUBSTRATE)

    assert product.canonical_bytes == smiles('C[C@@H](N)I').canonical_bytes
    assert product.canonical_bytes != smiles('C[C@H](N)I').canonical_bytes


def test_invert_is_silent_on_a_substrate_that_arrives_unconfigured():
    """`@~` is a conditional statement, so an unconfigured substrate is not a failure to report.

    "Whatever came in comes out the other way" is about a configuration the substrate has; one that
    arrives without it has nothing for the token to be about, and inventing one is the thing this whole
    language refuses.  Silent rather than logged because the alternative is a record on every
    unconfigured substrate of every corpus row that states inversion -- which is most of them.
    """
    log = []
    product = product_of(SUBSTITUTION + '>>[C@~:1][I;D1:2]', 'C[CH](N)Br', log=log)

    assert product.parity_of(2) == 0
    assert log == []


def test_one_template_covers_three_and_four_coordinate_centres_alike():
    """The measurable payoff of taking no frame: ONE inversion template, every centre shape.

    A statement drawn against a frame needs three or four directions named, so a template written for a
    centre with an implicit hydrogen could not match one with four heavy neighbours, and one whose arms
    were spelled `[C:2]` could not match an anionic arm -- `[A]` is neutral, and `[*]` is refused on a
    product side.  With no arms to name, none of that can go wrong: the three substrates here are a
    3-heavy centre, a 4-heavy centre and one carrying an alkoxide.
    """
    template = SUBSTITUTION + '>>[C@~:1][I;D1:2]'
    for substrate, inverted in (('C[C@H](N)Br', 'C[C@@H](N)I'),
                                ('C[C@](N)(O)Br', 'C[C@@](N)(O)I'),
                                ('C[C@H]([O-])Br', 'C[C@@H]([O-])I')):
        assert product_of(template, substrate).canonical_bytes == smiles(inverted).canonical_bytes


def test_keep_and_invert_are_reported_as_atom_sets():
    """What the reader hands the patcher: two sets of product atoms, and no values.

    Neither token names a sign, so there is nothing per atom to report -- which is the introspection
    shape of "relative, and no frame".
    """
    t = read_smirks(SUBSTITUTION + '>>[C@~:1][I;D1:2]')
    assert t.product_stereo_invert == frozenset({1}) and t.product_stereo_keep == frozenset()
    assert t.product_stereo_correlated == {}

    t = read_smirks(SUBSTITUTION + '>>[C@=:1][I;D1:2]')
    assert t.product_stereo_keep == frozenset({1}) and t.product_stereo_invert == frozenset()


# --- a reactant-side sign is selectivity ----------------------------------------------------------

def test_a_reactant_sign_narrows_the_match_and_states_nothing_about_the_product():
    """`[C;@]` is "configured", and that is the whole of its content on a reactant side.

    Two claims in one test, because they are one fact: the sign selects the substrates that arrive
    configured, and the product it yields is whatever the product side said -- the same product the
    unsigned template gives on the same substrate.
    """
    assert list(read_smirks(SELECTIVE + '>>[C@~:1][I;D1:2]')(smiles('C[CH](N)Br'))) == []

    selective = product_of(SELECTIVE + '>>[C@~:1][I;D1:2]', SUBSTRATE)
    plain = product_of(SUBSTITUTION + '>>[C@~:1][I;D1:2]', SUBSTRATE)
    assert selective.canonical_bytes == plain.canonical_bytes


def test_the_reactant_signs_own_character_says_nothing():
    """`@` and `@@` are the same reactant-side query, and with no arms named they have to be.

    With fewer than three directions named there is no frame, so the sign is unenforceable as a value
    and the matcher widens it to "configured, either sign".  `[C;@,@@]` is therefore accepted too and
    means the same thing written twice -- not a contradiction, since neither character is being compared
    with anything.
    """
    one = product_of('[C;@:1][Br;D1]>>[C@~:1][I;D1:2]', SUBSTRATE)
    other = product_of('[C;@@:1][Br;D1]>>[C@~:1][I;D1:2]', SUBSTRATE)
    both = product_of('[C;@,@@:1][Br;D1]>>[C@~:1][I;D1:2]', SUBSTRATE)

    assert one.canonical_bytes == other.canonical_bytes == both.canonical_bytes
    assert one.canonical_bytes == smiles('C[C@@H](N)I').canonical_bytes


def test_naming_the_arms_narrows_the_match_further():
    """With three directions named the reactant sign is a VALUE again -- and the product side is not.

    An anchor naming three or four directions has a frame, so its sign is enforced against it and the
    query selects one enantiomer.  Which is why reordering the arms swaps which character matches, and
    why exactly one of the two orders below matches this substrate.  The PRODUCT side is unaffected: in
    both matching spellings the inversion is `@~`, with no frame consulted.

    So arms are for narrowing the match, and a template that does not want to narrow simply omits them.
    """
    plain = product_of(SUBSTITUTION + '>>[C@~:1][I;D1:2]', SUBSTRATE)
    one = '[C;@%s:1]([C:2])([N:3])[Br;D1]>>[C@~:1]([C:2])([N:3])[I;D1:4]'
    other = '[C;@%s:1]([N:3])([C:2])[Br;D1]>>[C@~:1]([N:3])([C:2])[I;D1:4]'

    for template in (one, other):
        matching = [t for t in (template % '', template % '@')
                    if list(read_smirks(t)(smiles(SUBSTRATE)))]
        assert len(matching) == 1, 'a framed sign selects one enantiomer'
        assert product_of(matching[0], SUBSTRATE).canonical_bytes == plain.canonical_bytes


# --- every kind, one token ------------------------------------------------------------------------

def test_a_double_bond_at_the_reaction_centre_is_dropped_kept_or_turned_over():
    """The kind with no sign spelling, and all three answers for it in one place.

    Substituting at a vinyl carbon writes one of the unit's own two atoms, so the geometry goes unless
    the template speaks.  `@=` holds it and `@~` gives the other geometry -- the same two tokens as at a
    tetrahedral centre, because a cis/trans unit is two-state in exactly the same sense.
    """
    dropped = product_of(SUBSTITUTION + '>>[C:1][O;D1;H1:2]', 'C/C=C/Br')
    kept = product_of(SUBSTITUTION + '>>[C@=:1][O;D1;H1:2]', 'C/C=C/Br')
    turned = product_of(SUBSTITUTION + '>>[C@~:1][O;D1;H1:2]', 'C/C=C/Br')

    assert dropped.canonical_bytes == smiles('CC=CO').canonical_bytes
    assert kept.canonical_bytes == smiles('C/C=C/O').canonical_bytes
    assert turned.canonical_bytes == smiles('C/C=C\\O').canonical_bytes


def test_either_terminal_addresses_the_bond_kind():
    """Which terminal ANCHORS a cis/trans unit is a fact about slot order, not about chemistry.

    A template addresses the atom it means and the unit is found either way -- for `@=`, which spares it
    from the drop, and for `@~`, which has to locate the anchor to flip the parity stored there.
    """
    near_keep = product_of(SUBSTITUTION + '>>[C@=:1][O;D1;H1:2]', 'C/C=C/Br')
    far_keep = product_of('[C:3]=[C:1][Br;D1]>>[C@=:3]=[C:1][O;D1;H1:2]', 'C/C=C/Br')
    assert far_keep.canonical_bytes == near_keep.canonical_bytes

    near_turn = product_of(SUBSTITUTION + '>>[C@~:1][O;D1;H1:2]', 'C/C=C/Br')
    far_turn = product_of('[C:3]=[C:1][Br;D1]>>[C@~:3]=[C:1][O;D1;H1:2]', 'C/C=C/Br')
    assert far_turn.canonical_bytes == near_turn.canonical_bytes
    assert far_turn.canonical_bytes != far_keep.canonical_bytes


def test_the_token_addresses_whatever_unit_the_atom_anchors():
    """One token for every kind cuts both ways, and this is the edge of it.

    `[C;@]` asks only that the atom be CONFIGURED, and the bromine-bearing carbon of this alkene is: it
    anchors the cis/trans unit.  So a template written for an SN2 and applied to a vinyl halide turns
    that GEOMETRY over -- the token said "the other state of the unit here", and that is the unit here.
    A template meaning a tetrahedral centre and nothing else says so on its reactant side, which is why
    every corpus row that inverts carries `z1` or a degree.
    """
    turned = product_of(SELECTIVE + '>>[C@~:1][I;D1:2]', 'C/C(Br)=C/C')
    assert turned.canonical_bytes == smiles('C/C(I)=C\\C').canonical_bytes

    narrowed = '[C;@;z1:1][Br;D1]>>[C@~:1][I;D1:2]'
    assert list(read_smirks(narrowed)(smiles('C/C(Br)=C/C'))) == []
    assert product_of(narrowed, SUBSTRATE).canonical_bytes == smiles('C[C@@H](N)I').canonical_bytes


# --- enhanced stereo, in the bracket --------------------------------------------------------------

def test_a_group_is_written_in_the_bracket():
    """`&<n>` AND and `o<n>` OR, on the atom, which is the spelling templates use.

    CXSMILES puts the same two kinds in a `|...|` tail addressed by zero-based index over the atoms as
    written -- workable for a molecule serialized once, and miserable for a reaction, where one tail
    indexes both sides end to end.  In the bracket it names one atom on one side and needs no counting.
    """
    racemic = product_of(SELECTIVE + '>>[C;&1:1][I;D1:2]', SUBSTRATE)
    either = product_of(SELECTIVE + '>>[C;o1:1][I;D1:2]', SUBSTRATE)

    assert racemic.stereo_groups() == {(3, 1): [2]}
    assert either.stereo_groups() == {(2, 1): [2]}
    assert racemic.parity_of(2) and either.parity_of(2), 'a grouped centre is a configured one'


def test_a_group_on_a_created_centre_configures_it():
    """RACEMISATION IS CONFIGURED-AND-GROUPED, and this is the test that says so.

    An unconfigured centre and a racemic one are different facts, and only one of them can be
    depicted, canonicalised or told apart from a CH2 -- so a racemate is not the absence of a
    configuration.  The substrate here has no configuration for the patch to carry, so an
    implementation that only ever CLEARED a parity would leave parity 0 with a group beside it: a
    mixture of one unconfigured thing.
    """
    product = product_of('[C;D3;z1:1]([C:2])([N:3])[Br;D1]>>[C;&1:1]([C:2])([N:3])[O;D1:5]',
                         'CC(N)Br')

    assert product.parity_of(2) != 0
    assert product.stereo_groups() == {(3, 1): [2]}


def test_a_group_needs_a_stereogenic_unit_and_says_so_when_there_is_none():
    """A group on a CH2 is a false claim about a mixture, so no group is written and a line is logged.

    The atom the template groups here has two hydrogens after the patch, so nothing about it can be
    one of a set of configurations.  Logged and skipped rather than refused: whether the product atom
    is stereogenic depends on the MOLECULE, which the string cannot see.
    """
    log = []
    product = product_of('[C;D2;z1:1]([C:2])[Br;D1]>>[C;&1:1]([C:2])[H]', 'CCBr', log=log)

    assert product.stereo_groups() == {}
    assert len(log) == 1 and 'no stereogenic unit' in log[0].message


@mark.parametrize('smirks', [
    '[C:1]=[C:2][Br;D1]>>[C:1]=[C;&1:2][O;D1;h1:3]',      # the anchor
    '[C:1]=[C:2][Br;D1]>>[C;&1:1]=[C:2][O;D1;h1:3]',      # and the far terminal, which anchors nothing
])
def test_a_group_on_a_double_bond_is_an_e_z_mixture(smirks):
    """A GROUP MEANS "both of this unit's two states", and a cis/trans unit has two of them.

    A vinyl substitution with no geometric control gives an E/Z mixture, which is the same statement a
    racemate is and wants the same notation.  The group moves onto the unit's anchor whichever terminal
    the template named it on, because the parity and the group byte are one statement about one unit
    and the matcher reads both off the anchor -- so the two spellings give one answer.
    """
    product = product_of(smirks, 'C/C=C/Br')

    assert format(product, 'x') == 'C/C=C/O |&1:1|'
    assert product.stereo_groups() == {(3, 1): [2]}


def test_a_group_fabricates_the_geometry_the_substrate_never_drew():
    """The N3 rule holds for every kind: a grouped unit is a configured one.

    The substrate's double bond arrives unconfigured, so there is nothing to carry -- and a group
    written beside no geometry would say "not known, and a mixture", which are two different claims.
    So the patcher writes one of the two states and groups it, exactly as it does for a created centre.
    """
    product = product_of('[C:1]=[C:2][Br;D1]>>[C:1]=[C;&1:2][O;D1;h1:3]', 'CC=CBr')

    assert product.parity_of(2) != 0
    assert product.stereo_groups() == {(3, 1): [2]}


def test_an_allene_is_grouped_on_its_centre_and_only_there():
    """An allene is named on ONE atom (`stereo_unit_partner`), so its group has one place to go.

    Not an inconsistency with the double bond's two spellings: which atoms name a unit is a fact about
    the kind, and the same rule sends `@~` to the same atom.  A group on a terminal reaches nothing and
    is logged away with everything else that names no stereogenic unit.
    """
    axis = '[C:1]=[C:2]=[C:3][Br;D1]>>[C:1]=[C;&1:2]=[C:3][O;D1;h1:4]'
    end = '[C:1]=[C:2]=[C:3][Br;D1]>>[C;&1:1]=[C:2]=[C:3][O;D1;h1:4]'
    log = []

    product = product_of(axis, 'CC(F)=[C]=C(F)Br')
    assert product.parity_of(4) != 0
    assert product.stereo_groups() == {(3, 1): [4]}

    assert product_of(end, 'CC(F)=[C]=C(F)Br', log=log).stereo_groups() == {}
    assert len(log) == 1 and 'no stereogenic unit' in log[0].message


def test_one_unit_takes_one_group():
    """Both terminals grouped is one statement made twice, and the two brackets may not even agree.

    Resolved rather than refused, for the reason every other group case is: which atoms share a unit is
    a fact about the patched molecule.  The lower id wins and the collision is reported, so a template
    that meant two mixtures learns it got one.
    """
    log = []
    product = product_of('[C:1]=[C:2][Br;D1]>>[C;&1:1]=[C;&2:2][O;D1;h1:3]', 'CC=CBr', log=log)

    assert len(product.stereo_groups()) == 1
    assert len(log) == 1 and 'one unit takes one group' in log[0].message


def test_a_template_group_id_is_renumbered_around_the_input():
    """Template group ids are template-local, so an id the input already uses is not reused.

    The substrate carries `&1` on the very atom the template groups as its own `&1`, and the two are
    unrelated statements that happen to have collided on a number.  The template's group gets a free
    id instead of joining a mixture the atom may not belong to -- the same reasoning as N8's rule that
    a cleared parity takes its group with it.
    """
    product = product_of(SELECTIVE + '>>[C;&1:1][I;D1:2]', SUBSTRATE + ' |&1:1|')

    assert product.stereo_groups() == {(3, 2): [2]}


def test_a_group_the_template_does_not_state_is_carried():
    """Explicit-only cuts both ways: no group on the product side is no statement about the group.

    The template keeps the configuration and states no group, so the input's AND membership survives --
    the patcher has not been told the centre left the mixture and does not guess that it did.
    """
    product = product_of(SUBSTITUTION + '>>[C@=:1][I;D1:2]', SUBSTRATE + ' |&1:1|')

    assert product.stereo_groups() == {(3, 1): [2]}


def test_the_absolute_group_still_has_only_a_tail_spelling():
    """`|a:idx|` survives, and has no bracket form, because `a` in a bracket is the aromatic flag.

    An absolute group is the CXSMILES default and states nothing a template needs to say twice, so
    losing the bracket spelling costs nothing.  It is asserted here so the tail's routing is not
    silently dropped along with the notation that replaced its other two kinds.
    """
    product = product_of(SUBSTITUTION + '>>[C@=:1][I;D1:2] |a:2|', SUBSTRATE)

    assert product.stereo_groups() == {(1, 0): [2]}


def test_a_group_is_not_something_a_query_can_test():
    """The one bracket token that a plain query and a SMIRKS REACTANT side both refuse.

    A group says a configuration is one of a set -- a fact about a molecule and its mixture, not a
    property of an atom -- so there is nothing in a target to compare it with.  Both doors are shut by
    one line, in the seal, which is the only thing `read_smarts` and a reactant side have in common.
    """
    with raises(IncorrectSmarts):
        read_smarts('[C;&1]')
    with raises(IncorrectSmirks):
        read_smirks('[C;&1:1][Br;D1]>>[C:1][I;D1:2]')


def test_the_bracket_group_does_not_shadow_the_high_and_or_aromatic_oxygen():
    """One character of lookahead is all the group token costs, and this is the boundary it draws.

    `&` followed by a digit is a group and `&` followed by anything else is Daylight's high AND; `o`
    followed by a digit is a group and a bare `o` is aromatic oxygen.  No primitive name is a digit,
    so the group token collides with nothing -- the only cost is that an aromatic oxygen in a group
    has to be spelled `[o;o1]`.
    """
    assert read_smarts('[C&D1]').is_substructure(smiles('CC'))
    assert read_smarts('[o]').is_substructure(smiles('c1ccoc1'))
    assert read_smarts('[o;D2]').is_substructure(smiles('c1ccoc1'))
    # and the group spellings all reach the seal's refusal rather than the lexer's "names no primitive"
    for text in ('[C;&1]', '[C&1]', '[C;o1]', '[o;o1]'):
        with raises(IncorrectSmarts, match='enhanced-stereo group'):
            read_smarts(text)


def test_a_bracket_group_and_a_tail_group_on_one_atom_are_refused():
    """One atom, one group.  Two statements of it is a contradiction, not a precedence question.

    The bracket groups are collected before the tail is read, so the tail is the side that can see the
    collision -- which is why the refusal lives there and names the field.
    """
    with raises(IncorrectSmirks, match='one atom, one group'):
        read_smirks(SELECTIVE + '>>[C;&1:1][I;D1:2] |&2:2|')


# --- against what the arena did (N7, N8) ----------------------------------------------------------

def test_a_configuration_the_patch_destroyed_is_reported_and_not_invented():
    """When the patch outruns `rebase_parity` there is nothing left for `@~` to turn over.

    Replacing TWO of the anchor's four directions is past what the arena can re-base a sign against, so
    the parity is dropped and reported (N7).  `@~` cannot rescue it: "the other state" is a statement
    ABOUT a configuration, and the patcher may not invent one to make the template come true.  So the
    drop record stands alone -- one line, not two, because the token adds nothing where it has nothing
    to be about.
    """
    log = []
    product = product_of('[C:1]([C:2])([N;D1:3])[Br;D1]>>[C@~:1]([C:2])([O;D1:6])[S;D1:7]',
                         SUBSTRATE, log=log)

    assert product.parity_of(2) == 0
    assert len(log) == 1 and log[0].atoms == (2,)
    assert 'past what the arena could re-base the sign against' in log[0].message


def test_a_template_group_is_not_cleared_with_a_dropped_parity():
    """N8 clears the group of a parity the arena dropped -- and not of one the template restates.

    Same patch, same substrate carrying `&1`, three outcomes: unstated drops the parity and clears the
    group with it; a group alone re-configures the centre and gives it a fresh id; and the id is fresh
    rather than the input's 1, because the two statements are unrelated.  Asserted together because
    the clear and the write happen in one edit scope and an ordering bug between them is invisible
    from any one of the three.
    """
    swap_two = '[C:1]([C:2])([N;D1:3])[Br;D1]>>[C%s:1]([C:2])([O;D1:6])[S;D1:7]'
    substrate = SUBSTRATE + ' |&1:1|'

    assert product_of(swap_two % '', substrate).stereo_groups() == {}
    grouped = product_of(swap_two % ';&1', substrate)
    assert grouped.stereo_groups() == {(3, 2): [2]}
    assert grouped.parity_of(2) != 0, 'the group re-configured what the arena dropped'


# --- refusals, at the string ----------------------------------------------------------------------

def test_a_sign_with_nothing_to_be_relative_to_raises():
    """A bare product sign could only be an absolute setting, so it is not a template.

    Absolute configuration has no spelling, and this is not a gap to be worked around: a reaction that
    creates a centre out of an achiral one either racemises it or states a configuration relative to one
    the substrate already had.  Two spellings of the same mistake -- a sign on a mapped atom, and a sign
    on a product atom that pairs with no reactant atom at all.  Signing the reactant side does not help,
    since a reactant sign is a query and no product statement reads it.
    """
    with raises(IncorrectSmirks, match='racemise'):
        read_smirks(SUBSTITUTION + '>>[C;@:1][I;D1:2]')
    with raises(IncorrectSmirks, match='relative to'):
        read_smirks(SUBSTITUTION + '>>[C:1][C;@;D1:2]')
    with raises(IncorrectSmirks, match='racemise'):
        read_smirks(SELECTIVE + '>>[C;@@:1][I;D1:2]')


def test_a_sign_on_a_non_tetrahedral_anchor_raises():
    """A sign states a tetrahedral configuration and there is no other reading of it.

    A cis/trans, allene or atropisomer configuration has no sign spelling -- `@=` carries one through a
    patch, `@~` turns it over and `/` `\\` state a geometry outright -- so a sign on an anchor whose
    product bonds are not all single is a template that cannot mean anything, and it fails at
    `read_smirks` rather than at some later application against some later molecule.  Both a double bond
    and an aromatic one are checked: the question is the order the patch BUILDS, not the element.
    """
    with raises(IncorrectSmirks, match='TETRAHEDRAL'):
        read_smirks('[C:1]=[C:2]>>[C;@:1]=[C:2]')
    with raises(IncorrectSmirks, match='TETRAHEDRAL'):
        read_smirks('[C:1]([C:2])([N:3])[Br;D1]>>[C;@:1](:[C:2])([N:3])[O;D1:5]')


def test_a_created_atom_can_neither_keep_nor_invert():
    """`@=` and `@~` on an atom that pairs with nothing: there is no configuration for either to be about."""
    with raises(IncorrectSmirks, match='pairs with no reactant atom'):
        read_smirks(SUBSTITUTION + '>>[C:1][C@=:2]')
    with raises(IncorrectSmirks, match='pairs with no reactant atom'):
        read_smirks(SUBSTITUTION + '>>[C:1][C@~:2]')


def test_keep_invert_and_a_sign_are_one_field():
    """All three state the configuration, so stating two of them is stating it twice."""
    for product in ('[C;@;@=:1][I;D1:2]', '[C;@;@~:1][I;D1:2]', '[C;@=;@~:1][I;D1:2]'):
        with raises(IncorrectSmirks, match='states the stereo sign twice'):
            read_smirks(SELECTIVE + '>>' + product)


def test_keep_and_invert_cannot_be_part_of_a_query():
    """A query has no reactant to have had a configuration, so the seal refuses both tokens.

    The lexer reads them on either side of the arrow, because one lexer serves both -- the refusal is at
    the seal, which is the only thing the product side never reaches.  Same shape as `#0`.
    """
    for token in ('@=', '@~'):
        with raises(IncorrectSmirks, match='cannot be part of a query'):
            read_smirks('[C%s:1][Br;D1]>>[C:1][I;D1:2]' % token)
        with raises(IncorrectSmarts, match='cannot be part of a query'):
            read_smarts('[C%s]' % token)


# --- relative configuration, which is a correlated group (N12) ------------------------------------

# Directed epoxidation of an allylic alcohol: the oxygen is delivered syn to the hydroxyl, so what the
# reaction fixes is the configuration of the three centres RELATIVE to each other and not any one of
# them absolutely.  Every signed member names three or four of its directions, which this class of
# template gets for free -- the substituents are already in the pattern, because the selectivity is
# what depends on them.  The reactant side signs nothing at all: that is what lets one template take
# both a configured and an unconfigured substrate.
EPOXIDATION = ('[C:1]([O;D1;h1:2])([C:6])[C:3](-[C:7])=[C:4]-[C:8]'
               '>>[C;&1;@:1]([O:2])([C:6])[C;&1;@:3]1([C:7])[C;&1;@@:4]([C:8])[O:5]1')
# pent-3-en-2-ol, and the diastereomer the template draws from its (R) enantiomer
ALLYLIC = 'C[C@H](O)C(C)=CC'
SYN = 'O1[C@@H](C)[C@@]1([C@H](C)O)C'


def test_a_correlated_group_resolves_against_a_carried_configuration():
    """The substrate knows its own absolute configuration, so the product does too -- and gets no group.

    Case 3 of N3: some member arrives configured, so the drawn set is mirrored as a whole if that is
    what agreeing with it takes, and what comes out is a single diastereomer of known absolute
    configuration.  No group, because there is no mixture left to describe.
    """
    product = product_of(EPOXIDATION, ALLYLIC)

    assert product.canonical_bytes == smiles(SYN).canonical_bytes
    assert product.stereo_groups() == {}


def test_the_same_template_on_an_achiral_substrate_gives_the_racemate_of_one_diastereomer():
    """Case 2: nothing decided which enantiomer, so the honest answer says so -- and says which pair.

    All three centres in one AND group, which is exactly "as drawn, or all three flipped": the relative
    configuration is fixed and the absolute one is not.  That is the product of epoxidising a racemic
    allylic alcohol, and no absolute-setting notation could have expressed it.

    `canonical_bytes` does not carry group membership, so the diastereomer is checked against `SYN` and
    the mixture against `stereo_groups()` -- one assertion each, for two separate claims.
    """
    product = product_of(EPOXIDATION, 'CC(O)C(C)=CC')

    assert product.canonical_bytes == smiles(SYN).canonical_bytes
    assert len(product.stereo_groups()) == 1
    assert sorted(next(iter(product.stereo_groups().values()))) == [2, 4, 6]
    assert next(iter(product.stereo_groups()))[0] == 3, 'AND, not OR: one diastereomer, both hands'


def test_the_enantiomeric_substrate_gives_the_ENANTIOMERIC_product():
    """The measurement that makes absolute setting indefensible, run as a test.

    Feed the mirror-image alcohol and every sign in the product flips: the relative configuration the
    template states is preserved and the absolute one follows the substrate, which is what a directed
    epoxidation actually does.  A template that had spelled its product's configuration absolutely
    would have returned `SYN` here as well -- the same absolute product from both enantiomers, which is
    not a reaction.
    """
    product = product_of(EPOXIDATION, 'C[C@@H](O)C(C)=CC')

    assert product.canonical_bytes == smiles('O1[C@H](C)[C@]1([C@@H](C)O)C').canonical_bytes
    assert product.canonical_bytes != smiles(SYN).canonical_bytes
    assert product.stereo_groups() == {}


def test_the_drawn_set_is_reported_as_signs_and_frames():
    """What the reader hands the patcher: the sign, and the arm order it is drawn against.

    This is the one reading where a frame is computed, and where it belongs -- a statement about several
    centres at once has to be drawn somewhere.  Ascending product-atom order, `None` in the slot of an
    implicit hydrogen, which is the shape `translate_stereo` takes.  `product_stereo_keep` and
    `product_stereo_invert` are empty: the three statements are disjoint by construction.
    """
    t = read_smirks(EPOXIDATION)

    assert t.product_stereo_keep == frozenset() and t.product_stereo_invert == frozenset()
    assert t.product_stereo_correlated == {1: (1, (2, 3, 4, None)),
                                           4: (1, (1, 5, 6, 8)),
                                           6: (2, (4, 7, 8, None))}
    assert t.product_stereo_groups == {1: (3, 1), 4: (3, 1), 6: (3, 1)}


def test_the_OR_kind_correlates_the_same_way():
    """`o<n>` groups its members as "one of the two, nobody knows which", and correlates identically.

    The kind decides what the group CLAIMS about the mixture, not whether the signs inside it are
    relative to each other -- so a template that would rather say "a single enantiomer, set by a
    reagent this pattern does not name" writes `o` and gets the same relative configuration.
    """
    t = read_smirks(EPOXIDATION.replace('&1', 'o1'))
    assert len(t.product_stereo_correlated) == 3

    product = product_of(EPOXIDATION.replace('&1', 'o1'), 'CC(O)C(C)=CC')
    assert product.canonical_bytes == smiles(SYN).canonical_bytes
    assert next(iter(product.stereo_groups()))[0] == 2, 'OR'


def test_a_correlated_member_states_three_or_four_directions():
    """A drawn configuration needs a frame, and a frame is three named directions or four.

    Not a widening candidate the way a reactant side's frameless sign is: there the sign is a query and
    the value is unenforceable, so widening it costs a template nothing.  Here the sign IS the value, so
    a member with two directions has stated something with no content, and the refusal is at the string.
    """
    with raises(IncorrectSmirks, match='names 2 direction'):
        read_smirks('[C:1]([O;D1;h1:2])([C:6])[C:3](-[C:7])=[C:4]-[C:8]'
                    '>>[C;&1;@:1]([O:2])[C;&1;@:3]1([C:7])[C;&1;@@:4]([C:8])[O:5]1')


def test_a_reactant_sign_beside_a_correlated_product_sign_is_only_selectivity():
    """Nothing to reconcile: one is a query, the other is a drawn configuration.

    Signing the epoxidation's carbinol on the reactant side narrows the template to a substrate that
    arrives configured, and the drawn set still resolves against whatever that configuration turns out
    to be.  Which is why the same product comes back as from the unsigned template -- the sign selected
    a substrate, it did not state a course.
    """
    signed = EPOXIDATION.replace('[C:1]([O;D1;h1:2])', '[C;@:1]([O;D1;h1:2])', 1)

    assert list(read_smirks(signed)(smiles('CC(O)C(C)=CC'))) == []
    assert product_of(signed, ALLYLIC).canonical_bytes == smiles(SYN).canonical_bytes


def test_signed_atoms_in_DIFFERENT_groups_are_not_correlated():
    """Correlation is group membership and nothing looser -- two groups are two statements.

    Split the epoxidation's three members across `&1` and `&2` and none of the groups has two signed
    members, so no sign has anything to be relative to and each falls through to the absolute reading,
    which does not exist.  Refused with the message that names all three alternatives.
    """
    with raises(IncorrectSmirks, match='relative to'):
        read_smirks('[C:1]([O;D1;h1:2])([C:6])[C:3](-[C:7])=[C:4]-[C:8]'
                    '>>[C;&1;@:1]([O:2])([C:6])[C;&2;@:3]1([C:7])[C;&1;@@:4]([C:8])[O:5]1')


def test_a_member_the_molecule_cannot_configure_is_logged_and_the_rest_still_apply():
    """A drawn member that lands on a non-centre is skipped, and does not take the group with it.

    2-methylbut-3-en-2-ol epoxidises to a product whose second epoxide carbon carries two identical
    methyls, so that centre is not stereogenic and no configuration of it exists to write.  Input is
    input: the member is logged and skipped, the members that CAN be configured still get the relative
    configuration they were drawn with, and the group covers those.
    """
    log = []
    product = product_of(EPOXIDATION, 'CC(O)C(C)=C(C)C', log=log)

    assert any('no tetrahedral centre' in record.message or
               'directions the patched molecule does not have' in record.message for record in log)
    assert len(product.stereo_groups()) == 1
    assert sorted(next(iter(product.stereo_groups().values()))) == [2, 4]

# --- a drawn geometry, `/` and `\` ----------------------------------------------------------------

# Acetaldehyde plus ethyl bromide onto but-2-ene, and hexa-2,4-diene onto itself: one template per
# geometry, differing in one character.  Neither reactant side states a geometry -- the first has none to
# state and the second deliberately declines to -- so the product's E or Z comes from the drawing alone.
OLEFINATION = '[C;h3:1][C;h1:2]=[O;D1:3].[C;h3:4][C;h2:5][Br;D1]>>[C:1]/[C:2]=[C:5]%s[C:4]'
ISOMERISATION = '[C;h3:1][C;h1:2]=[C;h1:3][C;h3:4]>>[C:1]/[C:2]=[C:3]%s[C:4]'
DIENE = ('[C;h3:1][C;h1:2]=[C;h1:3][C;h1:4]=[C;h1:5][C;h3:6]'
         '>>[C:1]/[C:2]=[C:3]/[C:4]=[C:5]/[C:6]')


def olefination_product(product_side, **kwargs):
    """The single product of the olefination, whose two reactants make `product_of` unusable."""
    reactions = list(read_smirks(OLEFINATION % product_side)(smiles('CC=O'), smiles('CCBr'), **kwargs))
    assert len(reactions) == 1, 'the fixture is single-site by construction'
    assert len(reactions[0].products) == 1
    return reactions[0].products[0]


def test_a_geometry_is_drawn_on_a_double_bond_the_reaction_creates():
    """The statement absolute configuration does not get, and the reason it is not the same statement.

    An olefination makes its alkene E or Z by mechanism, and neither answer needs a chiral influence to
    have acted: the enantiomeric substrate does not give the enantiomeric product, because an alkene's
    two faces are not enantiomeric.  So the template says which, in one character, and the two spellings
    give the two molecules.
    """
    assert olefination_product('/').canonical_bytes == smiles('C/C=C/C').canonical_bytes
    assert olefination_product('\\').canonical_bytes == smiles('C/C=C\\C').canonical_bytes


def test_a_drawn_geometry_overrides_whatever_arrived():
    """Absolute means absolute: three substrates, one answer.

    `@~` is conditional -- it turns over what came in and is silent where nothing did -- and a drawing is
    not.  E, Z and an unconfigured double bond all come out E, which is what makes this the statement an
    isomerisation is written with.
    """
    for substrate in ('C/C=C\\C', 'C/C=C/C', 'CC=CC'):
        assert product_of(ISOMERISATION % '/', substrate).canonical_bytes == \
            smiles('C/C=C/C').canonical_bytes


def test_the_drawing_is_reported_as_two_terminals_two_substituents_and_one_bit():
    """What the reader hands the patcher: the chain's two ends, one marked arm each, and `trans`.

    Keyed by the terminals rather than by the directed bonds, because the terminals are what a unit is
    anchored at -- and reduced to a bit at read time, so nothing downstream has to know which end of
    which bond a `/` was written from.
    """
    assert read_smirks(ISOMERISATION % '/').product_stereo_geometry == {(2, 3): (1, 4, True)}
    assert read_smirks(ISOMERISATION % '\\').product_stereo_geometry == {(2, 3): (1, 4, False)}
    assert read_smirks(SUBSTITUTION + '>>[C@~:1][I;D1:2]').product_stereo_geometry == {}


def test_a_direction_is_read_from_the_atom_it_is_written_from():
    """`[C:2](/[C:1])=` and `[C:1]/[C:2]=` are opposite statements, and that is the notation.

    A direction names a side relative to the atom written FIRST, so moving the substituent into a branch
    turns the statement over without changing a character of it.  Which is why nothing normalises the
    pair: `smk_journal_directions` is the one journal walk whose key is not low-first.
    """
    chain = ISOMERISATION % '/'
    branch = '[C;h3:1][C;h1:2]=[C;h1:3][C;h3:4]>>[C:2](\\[C:1])=[C:3]/[C:4]'
    turned = '[C;h3:1][C;h1:2]=[C;h1:3][C;h3:4]>>[C:2](/[C:1])=[C:3]/[C:4]'

    assert product_of(branch, 'CC=CC').canonical_bytes == product_of(chain, 'CC=CC').canonical_bytes
    assert product_of(turned, 'CC=CC').canonical_bytes == smiles('C/C=C\\C').canonical_bytes


def test_one_direction_between_two_double_bonds_states_both_geometries():
    """`C/C=C/C=C/C` is three directions doing four jobs, and the middle one does two of them.

    Which is why the statement is read from the CHAIN TERMINALS and not from the directed bonds: the
    single bond between the two alkenes is a substituent of one terminal of each, so a walk over the
    directions would find the second geometry unstated and refuse the string.
    """
    template = read_smirks(DIENE)
    assert template.product_stereo_geometry == {(2, 3): (1, 4, True), (4, 5): (3, 6, True)}
    assert product_of(DIENE, 'CC=CC=CC').canonical_bytes == smiles('C/C=C/C=C/C').canonical_bytes


def test_a_chain_the_molecule_holds_no_cis_trans_unit_for_is_logged_and_skipped():
    """An ODD number of double bonds is an allene, whose configuration has no `/` spelling.

    The reader does not check the kind -- how many double bonds a chain has is a fact about the string,
    but which unit the patched molecule holds is a fact about the molecule -- so this is a log record and
    the product comes out unconfigured.  Input is input, on the template's side of the arrow too.
    """
    log = []
    template = ('[C;h3:1][C;h1:2]=[C:3]=[C;h1:4][C;h3:5]'
                '>>[C:1]/[C:2]=[C:3]=[C:4]/[C:5]')
    product = product_of(template, 'CC=C=CC', log=log)

    assert any('holds no cis/trans unit' in record.message for record in log)
    assert product.parity_of(2) == 0 and product.parity_of(4) == 0


# --- the geometry a query asks for ----------------------------------------------------------------

# Crotyl bromide's displacement, selective for the E isomer: the alkene is not the reaction centre, so
# the product side states no geometry and the arena carries the one the reactant side demanded.
GEOMETRY_SELECTIVE = '[C:1]/[C:2]=[C:3]/[C:4][Br;D1]>>[C:1][C:2]=[C:3][C:4][I;D1:5]'


def test_a_query_geometry_matches_that_geometry_and_no_other():
    """The reactant-side reading of the same pair: not a drawing but a demand.

    E and Z are two molecules and the query separates them, which is what makes a template selective for
    one alkene of a mixture.  Written in the target's own spelling or in another one -- the demand is
    about the unit the target holds, and a SMILES string is not that unit (N9).
    """
    e, z = smiles('C/C=C/C'), smiles('C/C=C\\C')

    assert read_smarts('C/C=C/C').is_substructure(e)
    assert not read_smarts('C/C=C/C').is_substructure(z)
    assert read_smarts('C/C=C\\C').is_substructure(z)
    assert not read_smarts('C/C=C\\C').is_substructure(e)
    assert read_smarts('[CH3]/[CH]=[CH]/[CH3]').is_substructure(e)


def test_an_unconfigured_double_bond_answers_no_geometry_query():
    """Ruling F54 on this kind too: a target that never said is not a target that said either."""
    for pattern in ('C/C=C/C', 'C/C=C\\C'):
        assert not read_smarts(pattern).is_substructure(smiles('CC=CC'))


def test_a_query_direction_is_read_from_the_atom_it_is_written_from():
    """One notation on both sides of the arrow: `C(/C)=C/C` demands what `C/C=C/C` denies."""
    assert read_smarts('C(/C)=C/C').is_substructure(smiles('C/C=C\\C'))
    assert read_smarts('C(\\C)=C/C').is_substructure(smiles('C/C=C/C'))


def test_either_terminal_of_the_target_unit_may_anchor_it():
    """Which end anchors the target's unit is a fact about the TARGET's slot order, so both are read.

    But-2-ene's own symmetry says so twice: the query maps onto it in both orientations, so a reading
    that only tried the query's first terminal would answer one of them by accident.
    """
    assert len(list(read_smarts('C/C=C/C').get_mapping(smiles('C/C=C/C')))) == 2


def test_a_geometry_mixture_answers_either_query():
    """An E/Z mixture is one molecule holding both states, so both queries hit it -- ruling F86.

    The group's decision is made from the same SF_* mask the tetrahedral half returns, which is why a
    geometry needed no machinery of its own for it.
    """
    mixture = smiles('C/C=C/C |&1:1,2|')

    assert read_smarts('C/C=C/C').is_substructure(mixture)
    assert read_smarts('C/C=C\\C').is_substructure(mixture)


def test_a_query_over_a_diene_states_both_geometries():
    """Read from the terminals here as well, so the shared direction makes `C/C=C/C=C/C` two demands."""
    assert read_smarts('C/C=C/C=C/C').is_substructure(smiles('C/C=C/C=C/C'))
    assert not read_smarts('C/C=C/C=C/C').is_substructure(smiles('C/C=C/C=C\\C'))
    assert read_smarts('C/C=C\\C=C/C').is_substructure(smiles('C/C=C\\C=C/C'))


def test_a_reactant_geometry_is_selectivity_and_the_arena_carries_the_geometry():
    """A geometry beside the arrow is matching selectivity, the same as a reactant-side sign is.

    E-crotyl bromide is displaced and the Z isomer is not a site at all.  Nothing on the product side
    mentions the alkene, and it is not the reaction centre either, so the E survives the patch.
    """
    template = read_smirks(GEOMETRY_SELECTIVE)

    assert product_of(GEOMETRY_SELECTIVE, 'C/C=C/CBr').canonical_bytes == \
        smiles('C/C=C/CI').canonical_bytes
    assert not list(template(smiles('C/C=C\\CBr')))
    assert not list(template(smiles('CC=CCBr')))


# --- a geometry's own refusals, at the string -----------------------------------------------------

def test_half_a_geometry_is_not_a_geometry():
    """One end marked says which side one substituent is on, which does not say which geometry it is."""
    with raises(IncorrectSmirks, match='a direction on one end only'):
        read_smirks('[C;h3:1][C;h1:2]=[C;h1:3][C;h3:4]>>[C:1]/[C:2]=[C:3][C:4]')


def test_both_substituents_of_one_terminal_on_one_side_is_refused():
    """No geometry puts a terminal's two substituents on the same side of the bond.

    Legal and usual to mark both -- `[C:2](/[F:5])(\\[C:1])=` says one thing twice -- so what is refused
    is the CONTRADICTION and not the redundancy.
    """
    both = '[C:1][C:2]([F:5])=[C:3][Br:6]>>[C:2](%s[F:5])(%s[C:1])=[C:3]/[Br:6]'
    read_smirks(both % ('/', '\\'))
    with raises(IncorrectSmirks, match='same side'):
        read_smirks(both % ('/', '/'))


def test_a_drawn_geometry_and_keep_or_invert_are_two_answers_to_one_question():
    """One states the geometry outright and the other takes the reactant's, so a string may not do both."""
    for token in ('@=', '@~'):
        with raises(IncorrectSmirks, match='two answers to one question'):
            read_smirks('[C;h3:1][C;h1:2]=[C;h1:3][C;h3:4]'
                        '>>[C:1]/[C%s:2]=[C:3]/[C:4]' % token)


def test_a_direction_that_names_no_double_bond_is_refused():
    """Dead surface, refused as N4 refuses its own: a `/` says nothing on its own.

    A single bond between two saturated carbons has no side to be on, and an aromatic bond is not a chain
    of double bonds either -- the question is the order the patch BUILDS.
    """
    with raises(IncorrectSmirks, match='name no chain of double bonds'):
        read_smirks('[C:1][C:2]>>[C:1]/[C:2]')
    with raises(IncorrectSmirks, match='name no chain of double bonds'):
        read_smirks('[C:1][C:2]=[C:3][O:4]>>[C:1]/[C:2]:[C:3]/[O:4]')


def test_the_same_refusals_hold_on_the_reactant_side():
    """Two readers, one vocabulary: the query seal walks the journal with the CSR, the patcher without.

    So each refusal above is asserted at the other door too -- the sides are read by different code, and
    a rule that held on only one of them would be a rule the dialect does not have.
    """
    with raises(IncorrectSmarts, match='a direction on one end only'):
        read_smarts('C/C=CC')
    with raises(IncorrectSmarts, match='same side'):
        read_smarts('C/C(\\C)=CC')
    with raises(IncorrectSmarts, match='no chain of double bonds'):
        read_smarts('C/CC')
    with raises(IncorrectSmirks, match='a direction on one end only'):
        read_smirks('[C:1]/[C:2]=[C:3][C:4]>>[C:1][C:2]=[C:3][C:4]')


def test_a_direction_combines_with_no_other_bond_token():
    """It carries the single bond itself, so `-` is already implied and an alternative has no reading.

    `!/` is the one that looks meaningful and is not: a bond has two sides, so "not this one" names the
    other one, and there is a token for it.
    """
    with raises(IncorrectSmarts, match='name every side but one'):
        read_smarts('[C]!/[C]')
    with raises(IncorrectSmarts, match='combines with nothing'):
        read_smarts('[C]/,\\[C]')
    with raises(IncorrectSmarts, match='combines with nothing'):
        read_smarts('[C]/;-[C]')


def test_a_direction_on_a_ring_closure_has_two_readings():
    """The label's two ends state the bond from opposite atoms, so a side there names no one side.

    Refused at either end rather than picked, since a template stating a geometry across a ring closure
    can restate the same bond in the chain.
    """
    with raises(IncorrectSmarts, match='two readings on a closure'):
        read_smarts('C/1=C/C1')
    with raises(IncorrectSmarts, match='two readings on a closure'):
        read_smarts('C=1CC/1')
