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
import pytest

from itertools import permutations

from chython.core import MoleculeContainer, QueryContainer


# the encoding MoleculeContainer.set_stereo_group validates:
# 0 unspecified, 1 abs, 2 or, 3 and -- note OR is 2 and AND is 3, not the other way round
SG_UNSPECIFIED, SG_ABS, SG_OR, SG_AND = 0, 1, 2, 3


# The carbon states implicit_h=1 because the core never derives one: without it the centre has
# three directions, anchors no unit, and every assertion below is about nothing.
def _chiral_target(parity, *, group=None, fourth=None):
    m = MoleculeContainer()
    with m.edit():
        elements = ('C', 'F', 'Cl', 'Br') if fourth is None else ('C', 'F', 'Cl', 'Br', fourth)
        sids = [m.add_atom(e, implicit_h=1 if k == 0 and fourth is None else 0)
                for k, e in enumerate(elements)]
        for s in sids[1:]:
            m.add_bond(sids[0], s, 1)
        m.set_parity(sids[0], parity)
        if group is not None:
            m.set_stereo_group(sids[0], group[0], group[1])
    return m, sids


# QueryContainer.add_atom() takes no element and add_bond() takes no order: both are primitives.
# There is no QueryContainer.set_stereo -- the primitive's name is 'stereo' (PRIM_NAMES).
# Follow chython/core/test/test_isomorphism.py's one_atom(), which is the established idiom.
def _chiral_query(parity):
    q = QueryContainer()
    qids = [q.add_atom() for _ in range(4)]
    for sid, number in zip(qids, (6, 9, 17, 35)):
        q.atom_primitive(sid, 'element', number)
    for s in qids[1:]:
        q.add_bond(qids[0], s)
        q.bond_primitive(qids[0], s, 'bond_order', 1)
    # a second primitive on the same atom needs an explicit operator between them, exactly as
    # test_isomorphism.py's one_atom() helper does
    q.atom_operator(qids[0], 'and_low')
    q.atom_primitive(qids[0], 'stereo', parity)
    return q, qids


def _chain_query(parity):
    """The same centre with one arm grown to a three-carbon chain: C(F)(Cl)CCC, stereo on atom 0."""
    q = QueryContainer()
    qids = [q.add_atom() for _ in range(6)]
    for sid, number in zip(qids, (6, 9, 17, 6, 6, 6)):
        q.atom_primitive(sid, 'element', number)
    for i, j in ((0, 1), (0, 2), (0, 3), (3, 4), (4, 5)):
        q.add_bond(qids[i], qids[j])
        q.bond_primitive(qids[i], qids[j], 'bond_order', 1)
    q.atom_operator(qids[0], 'and_low')
    q.atom_primitive(qids[0], 'stereo', parity)
    return q


def test_matching_parity_matches():
    """Also the readiness test: the stereo atom is first in the plan, so a kernel that tested the
    primitive when the anchor was bound would read three unmapped directions and refuse. Only an
    expected MATCH can tell readiness from luck -- an opposite-parity case is refused by a broken
    kernel and a correct one alike, which is why `test_opposite_parity_does_not_match` below proves
    nothing about when the check runs.
    """
    m, _ = _chiral_target(1)
    q, _ = _chiral_query(1)
    assert q.is_substructure(m)


def test_opposite_parity_does_not_match():
    m, _ = _chiral_target(1)
    q, _ = _chiral_query(2)
    assert not q.is_substructure(m)


def test_a_query_without_a_stereo_primitive_matches_either_parity():
    q = QueryContainer()
    qids = [q.add_atom() for _ in range(4)]
    for sid, number in zip(qids, (6, 9, 17, 35)):
        q.atom_primitive(sid, 'element', number)
    for s in qids[1:]:
        q.add_bond(qids[0], s)
        q.bond_primitive(qids[0], s, 'bond_order', 1)
    assert q.is_substructure(_chiral_target(1)[0])
    assert q.is_substructure(_chiral_target(2)[0])


def test_an_unconfigured_target_does_not_match_a_stereo_query():
    # a stereogenic centre with no stated parity: the unit exists, parity_of is 0
    m = MoleculeContainer()
    with m.edit():
        sids = [m.add_atom(e, implicit_h=1 if k == 0 else 0)
                for k, e in enumerate(('C', 'F', 'Cl', 'Br'))]
        for s in sids[1:]:
            m.add_bond(sids[0], s, 1)
    assert m.unit_of(sids[0]) is not None and m.parity_of(sids[0]) == 0
    q, _ = _chiral_query(1)
    assert not q.is_substructure(m)


def test_a_query_that_leaves_a_named_direction_unaccounted_refuses():
    # Ruling F67: the target centre has four heavy neighbours, the query names three, so one named
    # target direction is the image of nothing and the parity statement is about another molecule.
    # BOTH query signs are tried against both target parities: each sign is satisfiable by one
    # enantiomer, so a single sign would pin only half of the refusal.
    for sign in (1, 2):
        q, _ = _chiral_query(sign)
        for parity in (1, 2):
            m, sids = _chiral_target(parity, fourth='I')
            assert m.unit_of(sids[0]) is not None, 'the target centre really is a unit'
            assert not q.is_substructure(m), \
                f'under-specified query (sign {sign}) must not match parity {parity}'


def test_and_group_matches_either_parity():
    q, _ = _chiral_query(1)
    for parity in (1, 2):
        m, _ = _chiral_target(parity, group=(SG_AND, 1))
        assert q.is_substructure(m), f'AND must match parity {parity}'


def test_abs_group_matches_only_its_own_parity():
    q, _ = _chiral_query(1)
    assert q.is_substructure(_chiral_target(1, group=(SG_ABS, 0))[0])
    assert not q.is_substructure(_chiral_target(2, group=(SG_ABS, 0))[0])


def test_the_check_runs_at_the_last_direction_not_the_anchor():
    """Readiness at a depth greater than the anchor's: one arm is a three-carbon chain, so the
    anchor binds early while the fourth direction is reached much later. A kernel that tested the
    primitive at the anchor's own position -- or at any fixed offset from it -- refuses this match.
    """
    m = MoleculeContainer()
    with m.edit():
        sids = [m.add_atom(e, implicit_h=h) for e, h in
                (('C', 1), ('F', 0), ('Cl', 0), ('C', 2), ('C', 2), ('C', 3))]
        for i, j in ((0, 1), (0, 2), (0, 3), (3, 4), (4, 5)):
            m.add_bond(sids[i], sids[j], 1)
        m.set_parity(sids[0], 1)
    assert m.unit_of(sids[0]) is not None

    assert _chain_query(1).is_substructure(m)
    # the opposite primitive on the same shape must refuse, so the match above is not luck
    assert not _chain_query(2).is_substructure(m)


def test_unspecified_group_behaves_as_abs():
    q, _ = _chiral_query(1)
    assert q.is_substructure(_chiral_target(1, group=(SG_UNSPECIFIED, 0))[0])
    assert not q.is_substructure(_chiral_target(2, group=(SG_UNSPECIFIED, 0))[0])


def _query_in_element_order(numbers, parity, orders=None):
    """The `_chiral_query` shape with the neighbour ELEMENTS in the caller's creation order.

    The query's own ruling-F26 order is the order its atoms were created in, so this helper is how a
    test states a frame that differs from the target's.  `orders` gives the bond order per neighbour
    when the default single bond will not do; `parity` 0 leaves the stereo primitive off entirely,
    which is the control every refusal test below needs -- without it a refusal could just as well be
    the graph failing to match.
    """
    q = QueryContainer()
    qids = [q.add_atom() for _ in numbers]
    for sid, number in zip(qids, numbers):
        q.atom_primitive(sid, 'element', number)
    if orders is None:
        orders = [1] * (len(numbers) - 1)
    for s, order in zip(qids[1:], orders):
        q.add_bond(qids[0], s)
        q.bond_primitive(qids[0], s, 'bond_order', order)
    if parity:
        q.atom_operator(qids[0], 'and_low')
        q.atom_primitive(qids[0], 'stereo', parity)
    return q, qids


def test_the_primitive_is_read_in_the_querys_own_frame():
    """A sign is a statement about an ORDER, and each side has its own.

    The target's order is ruling F26 over the molecule -- here F, Cl, Br, then the implicit
    hydrogen's unnamed direction -- and the query's is the order its atoms were created in.  This
    query lists the same three neighbours BACKWARDS (Br, Cl, F), an odd permutation, so the sign that
    describes the very same three-dimensional arrangement is the opposite one: `@@` must match a
    target whose stored parity is `@`, and `@` must not.

    This is also why feature word IV's stereo bit cannot screen a stereo primitive and why
    `SPAN_COVERED[3]` cannot help either: bit 6 holds the parity in the MOLECULE's frame, the
    primitive's value is in the QUERY's, and no per-atom mask can compare them.
    """
    m, _ = _chiral_target(1)
    assert _query_in_element_order((6, 35, 17, 9), 2)[0].is_substructure(m)
    assert not _query_in_element_order((6, 35, 17, 9), 1)[0].is_substructure(m)
    # the same target read in its own order wants the sign it stores, so the flip above is the
    # frame's doing and not a global inversion
    assert _query_in_element_order((6, 9, 17, 35), 1)[0].is_substructure(m)


def test_a_stereo_primitive_on_a_non_tetrahedral_target_unit_refuses():
    """F77 case 1: the matched atom anchors a unit, but not an atom-kind one.

    `FC(Cl)=C(Br)I` anchors its cis/trans unit on the F/Cl carbon, and that carbon is the only atom
    this query's stereo atom can map to -- three neighbours, two of them F and Cl, the third reached
    by a double bond.  A tetrahedral sign says nothing about a double bond's geometry, so the kernel
    refuses rather than translating a parity through a frame of a different kind.  The parity is set
    so the refusal cannot be blamed on ruling F54's unconfigured case.

    ONLY THE PROBE PINS THE `kind == SU_TETRA` CLAUSE ITSELF, and no container-level case can: a
    cis/trans unit's refs are the substituents of BOTH alkene carbons (F, Cl, Br, I here), so the
    anchor's own third neighbour -- the far alkene carbon -- is never among them and the frame
    accounting refuses first (ruling F67).  Deleting the kind clause leaves this test green; deleting
    it breaks `test_the_frame_check_refuses_a_short_tetrahedral_record`.  This test's job is that the
    two refusals compose into the right container-level answer.
    """
    m = MoleculeContainer()
    with m.edit():
        sids = [m.add_atom(e) for e in ('C', 'F', 'Cl', 'C', 'Br', 'I')]
        for i, j, order in ((0, 1, 1), (0, 2, 1), (0, 3, 2), (3, 4, 1), (3, 5, 1)):
            m.add_bond(sids[i], sids[j], order)
        m.set_parity(sids[0], 1)
    unit = m.unit_of(sids[0])
    assert unit is not None and unit['kind'] != 0, 'the anchor really carries a non-atom kind'
    assert unit['parity'] == 1, 'and it really is configured'

    assert _query_in_element_order((6, 9, 17, 6), 0, orders=(1, 1, 2))[0].is_substructure(m), \
        'the graph itself matches, so a refusal below is the stereo primitive talking'
    for parity in (1, 2):
        assert not _query_in_element_order((6, 9, 17, 6), parity, orders=(1, 1, 2))[0] \
            .is_substructure(m), f'a tetrahedral sign must not read a cis/trans frame ({parity})'


def test_a_stereo_query_atom_with_too_few_neighbours_asks_only_that_it_be_CONFIGURED():
    """Fewer than three named directions is not a frame, so the SIGN is dropped and nothing else is.

    Two named neighbours leave two of the four slots to guess at, and there is no permutation to read
    a sign against -- so the value is unenforceable.  What remains enforceable is the half that never
    needed a frame: bit 7 says a configuration is there at all, and the box demands it whichever sign
    was written.  Hence `@` and `@@` are interchangeable HERE and both match a configured centre,
    while an unconfigured one matches neither.

    Refusing instead (F77 case 2) makes every such pattern match nothing -- including the SMIRKS
    spelling `[C;@:1][Br;D1]`, "a configured centre losing its bromide", which is how a template says
    it does not care what the other three directions are.  Two named directions
    rather than one, so the widening is not resting on a degenerate frame either way.
    """
    m, _ = _chiral_target(1)
    plain, _ = _chiral_target(0)
    assert _query_in_element_order((6, 9, 17), 0)[0].is_substructure(m)
    for parity in (1, 2):
        assert _query_in_element_order((6, 9, 17), parity)[0].is_substructure(m), \
            f'a configured centre satisfies "configured" whichever sign asked ({parity})'
        assert not _query_in_element_order((6, 9, 17), parity)[0].is_substructure(plain), \
            f'and an unconfigured one satisfies neither ({parity})'


def test_the_frame_check_refuses_a_short_tetrahedral_record():
    """F77 case 3, at the probe: a SU_TETRA record with fewer than four directions.

    No molecule can reach this branch -- perception has exactly one atom-kind emit site and it passes
    the literal 4 -- so the refusal is pinned here, where the two arguments can be chosen, plus a
    check on real molecules that the literal really is what comes out.  Kind 0 is SU_TETRA and kind 1
    is SU_CIS_TRANS, the same encoding `stereo_units()['kind']` reports.
    """
    from chython.core._core import _stereo_frame_probe

    assert _stereo_frame_probe(0, 4)
    assert not _stereo_frame_probe(0, 3), 'a short atom-kind record is refused'
    for kind in (1, 2, 3):
        assert not _stereo_frame_probe(kind, 4), \
            f'a full-width record of kind {kind} is not a tetrahedron'

    seen = 0
    for m in (_chiral_target(1)[0], _chiral_target(1, fourth='I')[0]):
        for unit in m.stereo_units():
            if unit['kind'] == 0:
                seen += 1
                assert unit['n_refs'] == 4, 'perception emits atom kinds at full width only'
    assert seen == 2, 'both molecules really did produce an atom-kind unit to check'


def test_an_and_group_does_not_excuse_an_under_specified_frame():
    """The F67 refusal runs BEFORE the group is consulted, and has to.

    An AND group says the target's sign is arbitrary within its group, so the sign is not compared --
    but a query that leaves one named target direction unaccounted for is not describing this
    molecule's centre at all, and no group can make it.  Both refusals are the same code path in a
    different order, so this is the test that keeps the group branch from swallowing F67.
    """
    for sign in (1, 2):
        q, _ = _chiral_query(sign)
        for parity in (1, 2):
            m, sids = _chiral_target(parity, group=(SG_AND, 1), fourth='I')
            assert m.unit_of(sids[0])['n_refs'] == 4
            assert not q.is_substructure(m), \
                f'AND must not rescue an under-specified frame (sign {sign}, parity {parity})'


def test_a_parity_validate_stereo_would_report_still_matches():
    """Ruling F76: the kernel reads the stated parity, not the stereogenicity marks.

    `CH3-CH(Cl)-CH3` anchors a tetrahedral unit -- four directions, two methyls, a chlorine and an
    implicit hydrogen -- that no molecule information depends on, so `mark_stereogenic` leaves it
    unmarked and `validate_stereo` reports the stored sign as unjustified.  The kernel matches it
    anyway: it goes through the unmarked door, translates the parity, and answers.  BOTH signs match
    here, and that is not a weakness of the test but the reason the unit is unmarked -- the two
    methyls are interchangeable, so one embedding translates to `@` and the other to `@@`.

    `validate_stereo()` CLEARS what it reports, which is why it is called last; the match after it is
    the proof that the stated parity, and nothing else, was what the primitive was reading.
    """
    m = MoleculeContainer()
    with m.edit():
        sids = [m.add_atom('C', implicit_h=1), m.add_atom('C', implicit_h=3),
                m.add_atom('C', implicit_h=3), m.add_atom('Cl')]
        for s in sids[1:]:
            m.add_bond(sids[0], s, 1)
        m.set_parity(sids[0], 1)
    unit = m.unit_of(sids[0])
    assert unit['parity'] == 1 and not unit['stereogenic'], 'configured, and not stereogenic'

    for parity in (1, 2):
        assert _query_in_element_order((6, 6, 6, 17), parity)[0].is_substructure(m), \
            f'an unjustified parity is still matchable ({parity})'

    assert m.validate_stereo() == [sids[0]], 'and validate_stereo does report it'
    assert not _query_in_element_order((6, 6, 6, 17), 1)[0].is_substructure(m), \
        'once cleared there is no parity to read, so the stereo query stops matching'


# ---------------------------------------------------------------------------
# Ruling F87 (the sign is per box) and ruling F88 (a drawn hydrogen)
# ---------------------------------------------------------------------------

def _query_with_tokens(numbers, tokens, orders=None):
    """`_query_in_element_order`, but the stereo atom's token stream is written out by the caller.

    `tokens` is a list of `('op', name)` and `(primitive_name, value)` pairs, appended to atom 0
    after its element -- which is how a test writes a disjunction like `[C@,N]`, the shape ruling F87
    is about.  Atom 0's element comes first, so a token list starting with an operator continues it.
    """
    q = QueryContainer()
    qids = [q.add_atom() for _ in numbers]
    for sid, number in zip(qids, numbers):
        q.atom_primitive(sid, 'element', number)
    if orders is None:
        orders = [1] * (len(numbers) - 1)
    for s, order in zip(qids[1:], orders):
        q.add_bond(qids[0], s)
        q.bond_primitive(qids[0], s, 'bond_order', order)
    for name, value in tokens:
        if name == 'op':
            q.atom_operator(qids[0], value)
        else:
            q.atom_primitive(qids[0], name, value)
    return q, qids


def _nitrogen_target():
    """N(F)(Cl)Br: three heavy neighbours on a nitrogen, no parity anywhere."""
    m = MoleculeContainer()
    with m.edit():
        sids = [m.add_atom(e) for e in ('N', 'F', 'Cl', 'Br')]
        for s in sids[1:]:
            m.add_bond(sids[0], s, 1)
    return m, sids


def test_a_sign_on_one_arm_of_a_disjunction_does_not_bind_the_other():
    """Ruling F87, case 1: `[C@,N]` is "an @-configured carbon, or any nitrogen".

    The nitrogen arm names no configuration, so it must match a nitrogen that states none.  With the
    sign hoisted to the atom this returned False and a real embedding was dropped in silence.  The
    control is the same query without the primitive, which shows the graph itself matches.
    """
    m, _ = _nitrogen_target()
    q, _ = _query_with_tokens((6, 9, 17, 35),
                              [('op', 'and_high'), ('stereo', 1), ('op', 'or'), ('element', 7)])
    control, _ = _query_with_tokens((6, 9, 17, 35), [('op', 'or'), ('element', 7)])
    assert control.is_substructure(m), 'the graph matches, so a refusal is the primitive talking'
    assert q.is_substructure(m), '[C@,N] must match a nitrogen through its sign-free arm'
    # box_counts is in plan order, and this query's root is not the stereo atom -- '[C,N]' allows
    # two elements, so no element bucket can seed it
    assert sorted(q.box_counts()) == [1, 1, 1, 2], 'and both arms are still there as two boxes'

    # the carbon arm keeps discriminating, so this is not the sign going unchecked
    assert q.is_substructure(_chiral_target(1)[0])
    assert not q.is_substructure(_chiral_target(2)[0])


def test_a_sign_paired_with_a_plain_primitive_leaves_that_primitive_alone():
    """Ruling F87, case 2, and the merge block that makes it work.

    `[C;@,D3]` matched through its D3 arm states nothing about configuration, so a carbon with no
    parity at all satisfies it.  The two boxes differ in one span only -- degree -- so a merge into one
    box would carry the `@` demand along with it (F87); `box_counts` pins that they stay two.
    """
    bare, sids = _chiral_target(1)
    with bare.edit():
        bare.set_parity(sids[0], 0)
    assert bare.parity_of(sids[0]) == 0, 'the target states no configuration'

    q, _ = _query_with_tokens((6, 9, 17, 35),
                              [('op', 'and_low'), ('stereo', 1), ('op', 'or'), ('degree', 3)])
    control, _ = _query_with_tokens((6, 9, 17, 35),
                                    [('op', 'and_low'), ('charge', 0), ('op', 'or'), ('degree', 3)])
    assert control.is_substructure(bare)
    assert q.is_substructure(bare), '[C;@,D3] must match through D3'
    assert sorted(q.box_counts()) == [1, 1, 1, 2], \
        'the sign difference blocks the merge (ruling F87)'


def test_a_centre_cannot_satisfy_both_signs_at_once():
    """Ruling F87, case 3: an ANDed pair of signs is unsatisfiable, not a construction error.

    `[C;@;@@]` seals, and refuses every target.  `[C;@,N;@@]` is the same contradiction reached
    through a disjunction -- its carbon arm ANDs both signs -- and an atom-level union of the two signs
    reads as "either sign will do" and MATCHES a parity-2 carbon.  The control replaces the first sign
    with a charge, leaving a single-sign box that does match.
    """
    both, _ = _query_with_tokens((6, 9, 17, 35),
                                 [('op', 'and_low'), ('stereo', 1),
                                  ('op', 'and_low'), ('stereo', 2)])
    assert both.box_counts() == [1, 1, 1, 1], 'it seals: an unsatisfiable query is constructible'
    for parity in (1, 2):
        assert not both.is_substructure(_chiral_target(parity)[0]), \
            f'no configuration satisfies both signs (parity {parity})'

    disjunction, _ = _query_with_tokens((6, 9, 17, 35),
                                        [('op', 'and_low'), ('stereo', 1), ('op', 'or'),
                                         ('element', 7), ('op', 'and_low'), ('stereo', 2)])
    control, _ = _query_with_tokens((6, 9, 17, 35),
                                    [('op', 'and_low'), ('charge', 0), ('op', 'or'),
                                     ('element', 7), ('op', 'and_low'), ('stereo', 2)])
    assert control.is_substructure(_chiral_target(2)[0]), 'the shape matches with one sign'
    assert not disjunction.is_substructure(_chiral_target(2)[0]), \
        '[C;@,N;@@] on a carbon ANDs the two signs and cannot match'


def _explicit_hydrogen_target(parity, hydrogens=1):
    """C(F)(Cl)(Br) with the fourth direction DRAWN as an H atom, and no implicit hydrogens.

    `hydrogens=2` drops the bromine for a second drawn H, which is the two-hydrogen centre ruling
    F88 must not open a door to.

    "AND NO IMPLICIT HYDROGENS" IS SAID OUT LOUD, in `implicit_h=0`, because the fixture's whole point
    is that the fourth direction is the DRAWN atom and not an implicit count.  Left to `add_atom`'s
    default the count is `H_UNKNOWN`, an unknown implicit count makes `total_h` unknown too, and the
    `H1` primitive stops matching a centre that plainly has one hydrogen.
    """
    m = MoleculeContainer()
    with m.edit():
        elements = ['C', 'F', 'Cl'] + (['Br', 'H'] if hydrogens == 1 else ['H', 'H'])
        sids = [m.add_atom(e, implicit_h=0) for e in elements]
        for s in sids[1:]:
            m.add_bond(sids[0], s, 1)
        m.set_parity(sids[0], parity)
    return m, sids


def test_a_drawn_hydrogen_at_the_centre_answers_like_an_implicit_one():
    """Ruling F88: whether a hydrogen is drawn is an input-representation choice.

    `stereo_units`' reference order puts hydrogen directions after heavy ones precisely so the tuple
    does not move when a hydrogen starts or stops being drawn; matching must not read it either, and
    MDL input carries explicit H at stereocentres routinely.  Both representations are asserted
    against both signs, so this pins the agreement and not just one half of it.
    """
    for parity, matching in ((1, 1), (2, 2)):
        drawn, dsids = _explicit_hydrogen_target(parity)
        implicit, isids = _chiral_target(parity)
        assert drawn.unit_of(dsids[0])['refs'][3] == dsids[4], 'the H really is a named direction'
        assert implicit.unit_of(isids[0])['refs'][3] is None, 'and here it really is unnamed'
        for sign in (1, 2):
            q, _ = _query_in_element_order((6, 9, 17, 35), sign)
            expected = sign == matching
            assert q.is_substructure(drawn) is expected, \
                f'drawn hydrogen, parity {parity}, sign {sign}'
            assert q.is_substructure(implicit) is expected, \
                f'implicit hydrogen, parity {parity}, sign {sign}'


def test_a_centre_with_two_drawn_hydrogens_is_not_a_frame_the_query_can_read():
    """F88 does not open a door to a non-stereocentre, and it has to close it itself.

    Perception does NOT refuse this centre a unit -- it records what the input stated, so
    `C([H])([H])(F)Cl` with a parity emits a full four-ref tetrahedral record -- and the kernel cannot
    consult `stereogenic` to notice, because that is the marked door's product and ruling F76 sends
    matching through the unmarked one.  So the refusal comes from the frame: two of the four
    directions are hydrogens, the query's unnamed direction could be either, and perception's
    tie-break between them is not an answer.  The control names the same three neighbours with no
    sign demanded and MATCHES -- so the graph is not what refuses, the primitive is.
    """
    m, sids = _explicit_hydrogen_target(1, hydrogens=2)
    unit = m.unit_of(sids[0])
    assert unit is not None and unit['n_refs'] == 4 and unit['parity'] == 1
    assert not unit['stereogenic'], 'two hydrogens: the centre is not stereogenic'

    control, _ = _query_in_element_order((6, 9, 17, 1), 0)
    assert control.is_substructure(m), 'the graph matches, so a refusal is the primitive talking'
    for sign in (1, 2):
        q, _ = _query_in_element_order((6, 9, 17, 1), sign)
        assert not q.is_substructure(m), f'an ambiguous hydrogen pairing must refuse (sign {sign})'


def test_the_implicit_hydrogen_primitive_decides_the_drawn_case_on_its_own():
    """The `h` primitive is an ordinary box screen and is independent of ruling F88.

    `h1` counts IMPLICIT hydrogens, so it separates the two representations -- that is its job, and
    F88 does not touch it.  `H1` counts total hydrogens and so accepts both.  Stated as a test
    because "a drawn hydrogen must not change the answer" could otherwise be over-read into the
    primitives that exist precisely to ask about it.
    """
    drawn, _ = _explicit_hydrogen_target(1)
    implicit, _ = _chiral_target(1)
    q_implicit, _ = _query_with_tokens((6, 9, 17, 35),
                                       [('op', 'and_low'), ('implicit_h', 1),
                                        ('op', 'and_low'), ('stereo', 1)])
    q_total, _ = _query_with_tokens((6, 9, 17, 35),
                                    [('op', 'and_low'), ('total_h', 1),
                                     ('op', 'and_low'), ('stereo', 1)])
    assert q_implicit.is_substructure(implicit)
    assert not q_implicit.is_substructure(drawn), 'h counts implicit hydrogens only'
    assert q_total.is_substructure(implicit)
    assert q_total.is_substructure(drawn), 'H counts both, and the sign still matches'


def test_a_query_without_a_stereo_primitive_does_not_build_the_unit_table():
    """Ruling F76's cost half: only QFLAG_HAS_STEREO opens the door in matcher_init.

    SEG_STEREO_UNIT is the one lazily built derived segment, and building it is what a stereo query
    pays for.  A query that names no configuration must not pay, and `total_len` is what notices:
    the arena grows by the table when the door is opened and not otherwise.  Without this test the
    guard could be deleted and every other test would still pass.
    """
    plain, _ = _query_in_element_order((6, 9, 17, 35), 0)
    stereo, _ = _query_in_element_order((6, 9, 17, 35), 1)

    m, _ = _chiral_target(1)
    before = m.total_len
    assert plain.is_substructure(m)
    assert m.total_len == before, 'a query with no stereo primitive builds no unit table'

    m2, _ = _chiral_target(1)
    assert m2.total_len == before
    assert stereo.is_substructure(m2)
    assert m2.total_len > before, 'and a stereo query does build one'


_TWO_CENTRE_BONDS = ((0, 1), (1, 2), (2, 3), (1, 4), (2, 5))
_TWO_CENTRE_ELEMENTS = (6, 6, 6, 6, 17, 9)      # C C C C Cl F -- 2-chloro-3-fluorobutane


def _two_centre_target(p1, p2, *, kind=SG_OR, group=1, groups=None):
    """2-chloro-3-fluorobutane: two independent stereocentres, each with a hydrogen direction.

    `groups` overrides `group` per centre as (g1, g2) when the two centres belong to different
    groups. Creation order matches `_two_centre_query`'s, so the two frames agree and the parity
    values below are directly comparable.
    """
    g1, g2 = (group, group) if groups is None else groups
    m = MoleculeContainer()
    with m.edit():
        sids = [m.add_atom(e, implicit_h=h) for e, h in
                (('C', 3), ('C', 1), ('C', 1), ('C', 3), ('Cl', 0), ('F', 0))]
        for i, j in _TWO_CENTRE_BONDS:
            m.add_bond(sids[i], sids[j], 1)
        m.set_parity(sids[1], p1)
        m.set_parity(sids[2], p2)
        m.set_stereo_group(sids[1], kind, g1)
        m.set_stereo_group(sids[2], kind, g2)
    assert m.unit_of(sids[1]) is not None and m.unit_of(sids[2]) is not None
    return m, sids


def _two_centre_query(p1, p2):
    q = QueryContainer()
    qids = [q.add_atom() for _ in _TWO_CENTRE_ELEMENTS]
    for sid, number in zip(qids, _TWO_CENTRE_ELEMENTS):
        q.atom_primitive(sid, 'element', number)
    for i, j in _TWO_CENTRE_BONDS:
        q.add_bond(qids[i], qids[j])
        q.bond_primitive(qids[i], qids[j], 'bond_order', 1)
    for sid, parity in ((qids[1], p1), (qids[2], p2)):
        q.atom_operator(sid, 'and_low')
        q.atom_primitive(sid, 'stereo', parity)
    return q, qids


def test_or_group_matches_the_stated_combination():
    m, _ = _two_centre_target(1, 2)
    q, _ = _two_centre_query(1, 2)
    assert q.is_substructure(m)


def test_or_group_matches_the_fully_inverted_combination():
    m, _ = _two_centre_target(1, 2)
    q, _ = _two_centre_query(2, 1)
    assert q.is_substructure(m), 'both units flip together, so the mirror image matches'


def test_or_group_rejects_a_partial_flip():
    m, _ = _two_centre_target(1, 2)
    q, _ = _two_centre_query(1, 1)
    assert not q.is_substructure(m), 'one flipped and one not is the combination OR forbids'


def test_and_group_accepts_the_partial_flip_that_or_rejects():
    m, _ = _two_centre_target(1, 2, kind=SG_AND)
    q, _ = _two_centre_query(1, 1)
    assert q.is_substructure(m), 'AND is per-unit, so the units are independent'


def test_independent_or_groups_decide_independently():
    m, _ = _two_centre_target(1, 2, groups=(1, 2))
    q, _ = _two_centre_query(1, 1)
    assert q.is_substructure(m), 'separate groups may flip separately'


def test_group_ids_are_renumbered_canonically():
    # the same molecule with the same members, differing only in the opaque input id
    a, _ = _two_centre_target(1, 2, group=1)
    b, _ = _two_centre_target(1, 2, group=63)
    groups = a.canonical_stereo_groups()
    assert groups, 'there is a group to renumber'
    assert groups == b.canonical_stereo_groups(), 'opaque input ids must not survive'
    # comparing sorted KEYS alone would pass on any two one-group molecules: compare memberships


def test_renumbering_is_a_bijection_on_the_ids_present():
    # ruling F68: groups partition the atoms, so nothing ever merges -- the count and the
    # partition are invariant and only the labels move
    m, sids = _two_centre_target(1, 2, groups=(7, 40))
    groups = m.canonical_stereo_groups()
    assert len(groups) == 2
    # 0 is not a group id -- set_stereo_group requires 1..63 for OR and AND -- so the canonical
    # ids are dense from 1, and the renumbering's output stays a legal input
    assert sorted(k[1] for k in groups) == [1, 2]
    assert sorted(v for members in groups.values() for v in members) == sorted(sids[1:3])


def test_canonical_groups_are_invariant_under_a_group_id_permutation():
    a, _ = _two_centre_target(1, 2, groups=(1, 2))
    b, _ = _two_centre_target(1, 2, groups=(2, 1))
    ga, gb = a.canonical_stereo_groups(), b.canonical_stereo_groups()
    assert len(ga) == 2, 'both groups are present, so the comparison is about something'
    assert ga == gb, 'a canonical id may not depend on which opaque id the caller picked'


def test_canonical_groups_are_invariant_under_a_wide_id_relabelling():
    a, _ = _two_centre_target(1, 2, groups=(1, 2))
    for lo, hi in ((3, 4), (7, 40), (62, 63), (40, 7)):
        b, _ = _two_centre_target(1, 2, groups=(lo, hi))
        assert a.canonical_stereo_groups() == b.canonical_stereo_groups(), f'ids {lo},{hi} moved the view'


# --- the decision variable: three cases the eight tests above do not separate -------------------
# Ruling F86 puts the OR decision at the complete mapping rather than in the DFS frame, and the
# tests above all have exactly one embedding and adjacent group members, so they would pass under
# a decision that never resets, never survives a distance, and never lets the search continue.

_SPACED_BONDS = ((0, 1), (1, 2), (1, 3), (3, 4), (4, 5), (5, 6), (5, 7))
# 2-chloro-5-fluorohexane: the two centres are four bonds apart, so between the position where the
# first group member is read and the position where the second is, the DFS crosses two CH2 atoms
_SPACED_ATOMS = (('C', 3), ('C', 1), ('Cl', 0), ('C', 2), ('C', 2), ('C', 1), ('F', 0), ('C', 3))
_SPACED_NUMBERS = (6, 6, 17, 6, 6, 6, 9, 6)


def _spaced_target(p1, p2, *, kind=SG_OR, groups=(1, 1)):
    m = MoleculeContainer()
    with m.edit():
        sids = [m.add_atom(e, implicit_h=h) for e, h in _SPACED_ATOMS]
        for i, j in _SPACED_BONDS:
            m.add_bond(sids[i], sids[j], 1)
        m.set_parity(sids[1], p1)
        m.set_parity(sids[5], p2)
        m.set_stereo_group(sids[1], kind, groups[0])
        m.set_stereo_group(sids[5], kind, groups[1])
    assert m.unit_of(sids[1]) is not None and m.unit_of(sids[5]) is not None
    return m, sids


def _spaced_query(p1, p2):
    q = QueryContainer()
    qids = [q.add_atom() for _ in _SPACED_NUMBERS]
    for sid, number in zip(qids, _SPACED_NUMBERS):
        q.atom_primitive(sid, 'element', number)
    for i, j in _SPACED_BONDS:
        q.add_bond(qids[i], qids[j])
        q.bond_primitive(qids[i], qids[j], 'bond_order', 1)
    for sid, parity in ((qids[1], p1), (qids[5], p2)):
        q.atom_operator(sid, 'and_low')
        q.atom_primitive(sid, 'stereo', parity)
    return q, qids


def test_an_or_group_decision_survives_the_distance_between_its_members():
    """The flip is the group's, not the frame's, so it has to outlive the positions in between.

    The eight tests above put the two members on adjacent atoms; here they are four bonds apart and
    only the fully inverted combination matches, so the choice made when the first member is read
    has to still be the choice when the second one is, two CH2 atoms later.  The partial flips are
    asserted in the same fixture: without them "matches" could be a decision that forgot.
    """
    m, _ = _spaced_target(1, 2)
    q, _ = _spaced_query(2, 1)
    assert q.is_substructure(m), 'both members flip together across the whole chain'
    assert _spaced_query(1, 2)[0].is_substructure(m), 'and the stated combination still matches'
    for p1, p2 in ((1, 1), (2, 2)):
        assert not _spaced_query(p1, p2)[0].is_substructure(m), \
            f'a partial flip ({p1}, {p2}) is what OR forbids however far apart the members are'


_FOUR_BONDS = ((0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (1, 6), (2, 7), (3, 8), (4, 9))
# 2-chloro-3-fluoro-4-bromo-5-iodohexane: four centres, each with one hydrogen direction
_FOUR_ATOMS = (('C', 3), ('C', 1), ('C', 1), ('C', 1), ('C', 1), ('C', 3),
               ('Cl', 0), ('F', 0), ('Br', 0), ('I', 0))
_FOUR_NUMBERS = (6, 6, 6, 6, 6, 6, 17, 9, 35, 53)
_FOUR_CENTRES = (1, 2, 3, 4)


def _four_centre_target(parities, groups=(1, 1, 2, 2), *, kind=SG_OR):
    m = MoleculeContainer()
    with m.edit():
        sids = [m.add_atom(e, implicit_h=h) for e, h in _FOUR_ATOMS]
        for i, j in _FOUR_BONDS:
            m.add_bond(sids[i], sids[j], 1)
        for k, parity, group in zip(_FOUR_CENTRES, parities, groups):
            m.set_parity(sids[k], parity)
            m.set_stereo_group(sids[k], kind, group)
    for k in _FOUR_CENTRES:
        assert m.unit_of(sids[k]) is not None, f'centre {k} anchors no unit'
    return m, sids


def _four_centre_query(parities):
    q = QueryContainer()
    qids = [q.add_atom() for _ in _FOUR_NUMBERS]
    for sid, number in zip(qids, _FOUR_NUMBERS):
        q.atom_primitive(sid, 'element', number)
    for i, j in _FOUR_BONDS:
        q.add_bond(qids[i], qids[j])
        q.bond_primitive(qids[i], qids[j], 'bond_order', 1)
    for k, parity in zip(_FOUR_CENTRES, parities):
        q.atom_operator(qids[k], 'and_low')
        q.atom_primitive(qids[k], 'stereo', parity)
    return q, qids


def test_two_or_groups_need_all_four_flip_combinations():
    """Two groups of two: each of the four assignments is the only one that answers some query.

    A decision that flipped both groups together, or that flipped only the group it met first,
    would answer three of these four correctly.  The mixed queries are here for the other half of
    the claim: the freedom is per group, and inside a group there is none.
    """
    m, _ = _four_centre_target((1, 1, 1, 1))
    for parities, why in (((1, 1, 1, 1), 'plain, plain'), ((1, 1, 2, 2), 'plain, flipped'),
                          ((2, 2, 1, 1), 'flipped, plain'), ((2, 2, 2, 2), 'flipped, flipped')):
        assert _four_centre_query(parities)[0].is_substructure(m), f'{why} is a legal assignment'
    for parities in ((2, 1, 1, 1), (1, 2, 1, 1), (1, 1, 2, 1), (1, 1, 1, 2), (2, 2, 2, 1)):
        assert not _four_centre_query(parities)[0].is_substructure(m), \
            f'{parities} splits a group, and no per-group choice can satisfy that'


def _two_chain_target(order):
    """Two disjoint 2-chloro-3-fluorobutanes in one container, each its own OR group.

    `order` is a pair of (parities, group) in creation order, which is also the order the kernel's
    root candidate loop walks -- so putting the chain that cannot satisfy the query first is how
    this fixture reaches a failing group decision before a passing one.
    """
    m = MoleculeContainer()
    chains = []
    with m.edit():
        for parities, group in order:
            sids = [m.add_atom(e, implicit_h=h) for e, h in
                    (('C', 3), ('C', 1), ('C', 1), ('C', 3), ('Cl', 0), ('F', 0))]
            for i, j in _TWO_CENTRE_BONDS:
                m.add_bond(sids[i], sids[j], 1)
            m.set_parity(sids[1], parities[0])
            m.set_parity(sids[2], parities[1])
            m.set_stereo_group(sids[1], SG_OR, group)
            m.set_stereo_group(sids[2], SG_OR, group)
            chains.append(sids)
    for sids in chains:
        assert m.unit_of(sids[1]) is not None and m.unit_of(sids[2]) is not None
    return m, chains


def test_a_failed_group_decision_does_not_end_the_search():
    """One site refuses the group demand and another satisfies it: the answer is yes.

    Both creation orders are asserted because which chain the root loop reaches first is the
    kernel's business: whichever it is, one of the two calls below has to walk past a complete
    mapping that the group decision refused and keep going.  `count` pins the other half -- exactly
    one of the two sites is an answer, so a decision that passed everything would say two.

    `get_mapping` is asserted as well as `count`, and for a reason of its own: it is the entry point
    that returns to Python between solutions and so re-borrows the arena through matcher_reseat.  In
    the order where the refusing site comes second, the refusal is made after a resume, so the state
    the resume rebuilds has to include the flag that turns the group decision on.
    """
    for first, second in ((((1, 2), 1), ((1, 1), 2)), (((1, 1), 1), ((1, 2), 2))):
        m, chains = _two_chain_target((first, second))
        q, _ = _two_centre_query(1, 1)
        assert q.is_substructure(m), f'the satisfiable site is found with {first} first'
        assert q.count(m) == 1, 'and only one of the two sites satisfies the group demand'
        assert len(list(q.get_mapping(m))) == 1, 'and the resuming entry point agrees'


def test_a_group_decision_does_not_leak_between_two_matches():
    """The decision is per solution, so two mappings out of one matcher must not share it.

    The container holds one site that needs the flip and one that does not, in separate groups, and
    a single query matches both.  A decision variable that survived from the first solution to the
    second would refuse the second (its group was already committed the other way) and `count`
    would read 1.  Repeating the call pins the same thing across matcher instances.
    """
    m, chains = _two_chain_target((((1, 1), 1), ((2, 2), 2)))
    q, _ = _two_centre_query(1, 1)
    assert q.count(m) == 2, 'one site plain, one site flipped, both in the same walk'
    assert q.count(m) == 2, 'and the second walk sees the same two'
    mapped = [sorted(mp.values()) for mp in q.get_mapping(m)]
    assert len(mapped) == 2 and mapped[0] != mapped[1], 'two distinct sites, not one twice'


# --- ruling F79: the canonical view is a read ---------------------------------------------------

def test_the_canonical_group_view_writes_nothing():
    """Ruling F79: no renumbering of the segment, so a shared arena cannot notice the call.

    `copy()` shares the arena, and `total_len` is what a segment append would move.  The group
    bytes are asserted afterwards through `stereo_groups()`, which reads the segment the caller
    could see renumbered: the opaque ids 7 and 40 have to still be 7 and 40.
    """
    m, sids = _two_centre_target(1, 2, groups=(7, 40))
    other = m.copy()
    assert m.shares_arena_with(other)
    before, generation = m.total_len, m.generation
    raw = m.stereo_groups()
    assert sorted(raw) == [(SG_OR, 7), (SG_OR, 40)], 'the opaque ids going in'

    view = m.canonical_stereo_groups()
    assert sorted(k[1] for k in view) == [1, 2], 'the view renumbered'
    assert m.total_len == before, 'the arena did not grow'
    assert m.generation == generation, 'and nothing rebound'
    assert m.shares_arena_with(other), 'still the same arena'
    assert m.stereo_groups() == raw, 'the stored ids are untouched'
    assert other.stereo_groups() == raw, 'and the sharer sees them unchanged too'
    assert other.canonical_stereo_groups() == view, 'the sharer reads the same view'


# --- ruling F89: the membership fixpoint, and the tie it cannot break ---------------------------
#
# Every fixture from here to the end of this section is a polychlorocycloalkane: ring carbons 0..n-1,
# each with its own chlorine and one implicit hydrogen (the core derives none, so the fixture states
# it).  Every ring centre is then a stereo unit, the constitution alone separates none of them, and
# the group partition is the only thing left to canonicalise on -- which is what makes this skeleton
# the witness for both halves of F89.  n = 4 is 1,2,3,4-tetrachlorocyclobutane, n = 6 the lindane
# skeleton, n = 8 octachlorocyclooctane.
#
# THE ENCODING AXIS, which every sweep below varies and which decides how a fixture may be written.
# Ruling F26 stores a unit's refs in CSR-slot ASCENDING order, so the parity BYTE is a function of the
# order the atoms were created in: one molecule written two ways is two byte patterns, and a fixture
# that holds the bytes fixed while permuting the creation order sweeps DIFFERENT MOLECULES.  So a
# fixture here names its stereochemistry in the frame the RING names -- (the next ring atom, the
# previous one, the chlorine, the hydrogen) -- and `_ring` sets whatever bytes reproduce that reading
# in whatever order the atoms went in.  The stored bytes really do differ between two encodings of one
# molecule, which `test_one_molecule_in_two_atom_orders_reads_the_same` asserts rather than assumes.
#
# The ring frame is also where these tests' GROUND TRUTH comes from, with no reference to the
# implementation.  A rotation carries the frame onto itself, so a rotation is a symmetry of the
# annotated molecule only if it PRESERVES every frame parity; a reflection reverses next/prev, one
# transposition, so it is a symmetry only if it INVERTS every one.  Two consequences are used
# repeatedly below.  A reflection through an ATOM fixes that atom and would need p == 3 - p there, so
# it is never a symmetry of a fully substituted ring: only reflections through bond midpoints can be.
# And whichever map survives must carry the group partition onto itself as a bijection, or it
# exchanges no group ids at all.


def _ring(order, frame, groups, kind=SG_OR):
    """A polychlorocycloalkane in one creation `order`, carrying `frame` and `groups`.

    `order` is the creation order of all 2n atoms -- carbons 0..n-1, then their chlorines n..2n-1;
    `frame[k]` is carbon k's parity IN THE RING FRAME and `groups[k]` its stereo group number.  The
    molecule is a function of `frame` and `groups` alone: `order` chooses only how it is encoded.
    """
    n = len(frame)
    m = MoleculeContainer()
    with m.edit():
        sids = [0] * (2 * n)
        for slot in order:
            sids[slot] = m.add_atom('C', implicit_h=1) if slot < n else m.add_atom('Cl')
        for k in range(n):
            m.add_bond(sids[k], sids[(k + 1) % n], 1)
            m.add_bond(sids[k], sids[k + n], 1)
        for k in range(n):
            m.set_parity(sids[k], 1)
            m.set_stereo_group(sids[k], kind, groups[k])
    for k in range(n):
        assert m.unit_of(sids[k]) is not None, 'every ring centre is a stereo unit'
    # which bytes to flip is decided BEFORE the edit opens, because translate_stereo refuses a
    # container with pending edits; flipping one centre's byte cannot change another centre's frame
    flip = [(k, 3 - m.parity_of(sids[k])) for k in range(n) if _ring_parity(m, sids, k) != frame[k]]
    with m.edit():
        for k, parity in flip:
            m.set_parity(sids[k], parity)
    assert tuple(_ring_parity(m, sids, k) for k in range(n)) == tuple(frame), 're-encoding failed'
    return m, sids


def _ring_parity(m, sids, k):
    """Carbon k's parity in the ring frame, which asks the same question in every encoding."""
    n = len(sids) // 2
    return m.translate_stereo(sids[k], (sids[(k + 1) % n], sids[(k - 1) % n], sids[k + n], None))


def _ring_view(m, sids):
    """(kind, id) -> the RING POSITIONS of the members, so two encodings can be compared at all."""
    slot = {s: k for k, s in enumerate(sids)}
    return {key: tuple(sorted(slot[a] for a in members))
            for key, members in m.canonical_stereo_groups().items()}


def _ring_orders(n, exhaustive=False):
    """Creation orders that vary the encoding and nothing else.

    Every dihedral relabelling of the ring, because that is where ruling F26's byte moves the most,
    and orders that are NOT ring symmetries so that a sweep cannot pass by handling only the
    symmetric ones: every permutation of the carbons where that is affordable, of the first four
    otherwise.  Each carbon order is taken with the chlorines created before it, after it, and
    interleaved with it, since the byte depends on the slot distance between a centre and its refs.
    """
    heads = [tuple((start + step * k) % n for k in range(n))
             for start in range(n) for step in (1, -1)]
    heads += ([tuple(p) for p in permutations(range(n))] if exhaustive else
              [tuple(p) + tuple(range(4, n)) for p in permutations(range(4))])
    chlorines = tuple(range(n, 2 * n))
    return [o for head in heads
            for o in (head + chlorines, chlorines + head,
                      tuple(x for pair in zip(head, chlorines) for x in pair))]


def _group_shapes(m, sids):
    """(kind, id) -> what the groups ARE, with no reference to any labelling.

    Member sets cannot be compared across two creation orders as stable ids (they are the
    relabelling) and must not be compared as `canonical_order()` positions either: that order is
    stereo-blind, so on a symmetric skeleton it is only canonical up to a symmetry that moves group
    members between groups -- measured on the 3+1 OR fixture below, whose view every reading certifies
    pinned while its join with `canonical_order()` positions takes twelve distinct values over the 96
    encodings, which is why the view is not joined to that order at all.  A group's size and its
    members' RING-FRAME parities are labelling-free, and they tell these fixtures' groups apart.  The
    frame is the load-bearing word: `parity_of` returns the stored byte, which ruling F26 makes a
    function of the creation order, so a shape built on it would compare two encodings of one
    molecule on the one thing that legitimately differs between them.
    """
    slot = {s: k for k, s in enumerate(sids)}
    return {key: (len(members), tuple(sorted(_ring_parity(m, sids, slot[a]) for a in members)))
            for key, members in m.canonical_stereo_groups().items()}


def test_one_molecule_in_two_atom_orders_reads_the_same():
    """Ruling F95: the reading is a function of the molecule, and the stored parity byte is not.

    1,2,3,4-tetrachlorocyclobutane, every ring carbon the same parity in the ring frame, OR groups
    {C1}, {C2} and {C3,C4}.  No rotation carries that partition onto itself and no reflection can
    preserve an all-equal frame at all, so every id is pinned and the whole id -> members map is a
    property of the molecule -- there is nothing here for an ambiguity class to hide.

    The test asserts the premise as well as the conclusion: the three encodings named below store
    three DIFFERENT byte patterns, so an implementation that seeds its refinement on the byte cannot
    pass.  Measured, and this is the defect the ruling names: with the byte in the seed the two
    encodings (0,1,2,3) and (1,2,0,3) of one molecule report `ambiguities() == ()` -- every id pinned
    -- while disagreeing about which group is which id.
    """
    reference = None
    bytes_seen = set()
    for order in _ring_orders(4, exhaustive=True):
        m, sids = _ring(order, (1, 1, 1, 1), (1, 2, 3, 3))
        bytes_seen.add(tuple(m.parity_of(sids[k]) for k in range(4)))
        view = _ring_view(m, sids)
        assert m.canonical_stereo_group_ambiguities() == (), \
            f'creation order {order}: this molecule has no symmetry that exchanges its groups'
        if reference is None:
            reference = view
            assert sorted(view.items()) == [((SG_OR, 1), (1,)), ((SG_OR, 2), (0,)),
                                            ((SG_OR, 3), (2, 3))], 'three groups, as built'
        else:
            assert view == reference, f'creation order {order} moved a pinned id'
    assert len(bytes_seen) > 1, \
        'the sweep must vary the ENCODING: one byte pattern would make this test about nothing'


def test_the_two_encodings_that_forced_ruling_f95_agree():
    """The witness itself, kept as a fixture: two creation orders, one molecule, one reading.

    1,2,3,4-tetrachlorocyclobutane, OR groups {C1,C2} and {C3,C4}, ring-frame parities (1, 2, 2, 2).
    Created in the order (C1, C2, C3, C4) it stores the bytes (1, 1, 1, 2); created in the order
    (C2, C3, C1, C4) it stores (1, 2, 1, 1) -- ruling F26 puts a unit's refs in slot ascending order,
    so the SAME stereochemistry has two byte patterns, and this test asserts that premise so it cannot
    pass by comparing a molecule with itself.

    At the byte-seeded seed this replaced, the two encodings reported `ambiguities() == ()` -- every id
    pinned -- while disagreeing about which group was id 1.  Two readings both claiming to be pinned and
    contradicting each other is what no ambiguity class can excuse, and it is the whole of ruling F95.

    What the fix costs is visible here too, and is asserted rather than hidden: no dihedral map of this
    ring exchanges the two groups while preserving the frame parities (a rotation by two would need
    parity 1 at C3, a reflection would have to invert all four), so the molecule is really pinned, and
    the reading reports one ambiguity class covering both groups anyway.  On a ring a centre's next and
    previous carbons carry the same colour until something separates them, so its parity has no frame
    the colouring can name and contributes nothing -- rule 3's licensed direction, saying less than is
    known rather than more.
    """
    a, asids = _ring((0, 1, 2, 3, 4, 5, 6, 7), (1, 2, 2, 2), (1, 1, 2, 2))
    b, bsids = _ring((1, 2, 0, 3, 4, 5, 6, 7), (1, 2, 2, 2), (1, 1, 2, 2))
    assert tuple(a.parity_of(asids[k]) for k in range(4)) == (1, 1, 1, 2)
    assert tuple(b.parity_of(bsids[k]) for k in range(4)) == (1, 2, 1, 1), \
        'the two encodings must differ in the stored bytes, or this test compares nothing'
    assert _ring_view(a, asids) == _ring_view(b, bsids), 'one molecule, one id -> members map'
    assert a.canonical_stereo_group_ambiguities() == b.canonical_stereo_group_ambiguities()
    assert a.canonical_stereo_group_ambiguities() == (frozenset({(SG_OR, 1), (SG_OR, 2)}),), \
        'reported tied though the molecule is pinned: the coarse side of rule 3, and it is agreed on'


def test_the_canonical_view_survives_a_relabelling_and_a_swap_of_the_stored_ids():
    """Two axes at once, because either alone admits a rule F80 forbids.

    Groups {C1,C2,C3} and {C4} on the tetrachlorocyclobutane: same kind, same element, and on this
    skeleton the same constitutional class, so nothing but the membership itself distinguishes them.
    The fixture is built in every encoding AND with the two stored ids exchanged, and the
    (kind, id) -> shape map must not move.

    Measured, both mutations: rank the groups by their ascending STORED byte and the id swap flips
    which group is id 1, so the map moves (that rule survived the fixture this test replaced, which
    varied the creation order only).  Take the membership feedback out of the seed -- ruling F80 step
    3 -- and the four ring carbons stay in one refinement class, the singleton lands on an arbitrary
    canonical position, and the creation order decides whether it is id 1 or id 2.
    """
    reference = None
    for groups in ((1, 1, 1, 2), (2, 2, 2, 1)):
        for order in _ring_orders(4, exhaustive=True):
            m, sids = _ring(order, (1, 2, 1, 2), groups)
            shapes = _group_shapes(m, sids)
            assert sorted(shapes.values()) == [(1, (2,)), (3, (1, 1, 2))], 'the two groups, as built'
            assert m.canonical_stereo_group_ambiguities() == (), \
                'sizes 3 and 1: no symmetry can exchange them, so both ids are pinned'
            if reference is None:
                reference = shapes
            else:
                assert shapes == reference, f'stored ids {groups}, creation order {order}, moved it'


def test_two_interchangeable_groups_are_reported_as_one_ambiguity_class():
    """Ruling F89's second half: a tie the fixpoint cannot break is exposed, not broken.

    Groups {C1,C2} and {C3,C4} on the same ring.  The rotation by two carries one onto the other and
    carries every RING-FRAME parity to an equal one under both patterns below, so it is an
    automorphism of the parity-annotated molecule and the two assignments of ids 1 and 2 describe the
    SAME mixture.  No refinement over the membership can separate them -- both groups have two
    members, one of each class, in every round -- so the view reports the pair as one ambiguity class
    instead of choosing on the creation order.

    What the tie costs is stated here as a class and not as a measurement, because on THIS fixture the
    two groups are indistinguishable in every labelling-free observable -- same size, same frame
    parities, same ring positions up to the rotation -- so the raw id -> members direction is measured
    to take exactly one value over the 96 encodings this sweep visits, and the fixture cannot witness
    the movement the class admits.  The witness for that lives in the off-ring AND fixture at the end
    of this file, where the two members ARE distinguishable and the raw direction takes two values.
    Here the assertion is that the report exposes the tie at all rather than picking a side.
    """
    for parities in ((1, 2, 1, 2), (1, 1, 1, 1)):
        reference = None
        for order in _ring_orders(4, exhaustive=True):
            m, sids = _ring(order, parities, (1, 1, 2, 2))
            shapes = _group_shapes(m, sids)
            ambiguities = m.canonical_stereo_group_ambiguities()
            assert len(ambiguities) == 1, f'{parities}: one class of interchangeable groups'
            assert ambiguities[0] == frozenset(shapes), \
                'and it covers both keys, so no id in this view may be compared on its own'
            assert len(shapes) == 2, 'two groups still, because F89 forbids merging them'
            collapsed = (sorted(shapes.values()), sorted(ambiguities[0]))
            if reference is None:
                reference = collapsed
            else:
                assert collapsed == reference, \
                    f'{parities}, creation order {order}: the collapsed reading moved anyway'


def test_the_view_does_not_merge_two_groups_it_cannot_tell_apart():
    """An ambiguity is not an equivalence: merging would change the mixture.

    Two OR groups of one member each describe four stereoisomers; one OR group of two members
    describes two.  So even where no invariant rule can say which group is which, the count and the
    memberships must survive -- what degrades is only the id -> members direction.
    """
    m, sids = _ring((0, 1, 2, 3, 4, 5, 6, 7), (1, 1, 1, 1), (1, 1, 2, 2))
    view = m.canonical_stereo_groups()
    assert len(view) == 2, 'two groups in, two groups out'
    assert sorted(len(v) for v in view.values()) == [2, 2], 'and the memberships are untouched'
    assert sorted(a for members in view.values() for a in members) == sorted(sids[:4])
    assert m.stereo_groups() == {(SG_OR, 1): [sids[0], sids[1]], (SG_OR, 2): [sids[2], sids[3]]}, \
        'the stored partition is what it was'


def test_the_second_fixpoint_round_separates_what_the_first_cannot():
    """Ruling F80 step 3 is a FIXPOINT and not one pass, and this is what the extra pass buys.

    Hexachlorocyclohexane with ring-frame parities (1, 1, 1, 2, 2, 2) and OR groups {C1}, {C2},
    {C3,C4}, {C5}, {C6}: four groups of one and one of two.  The molecule's only symmetry is the
    reflection through the midpoints of the C1-C6 and C3-C4 bonds -- it inverts every frame parity,
    which is what that parity pattern demands of a reflection -- and it exchanges C1 with C6 and C2
    with C5 while fixing the pair.  So the truth is two ambiguity classes of two groups each.

    Round 1 cannot say that much: all four singleton groups carry the same key (one member, one
    refinement class, because the constitution puts all six carbons in one class), so all four are
    tied together.  Round 1's labels colour the ring singleton/singleton/pair/pair/singleton/singleton,
    refining that colouring separates {C1,C6} from {C2,C5}, and round 2 hands the two orbits different
    keys.  Measured with the feedback loop capped: at one round the four singletons come back as ONE
    class, at two rounds the answer below is complete.
    """
    reference = None
    for order in _ring_orders(6):
        m, sids = _ring(order, (1, 1, 1, 2, 2, 2), (1, 2, 3, 3, 4, 5))
        view = _ring_view(m, sids)
        by_slot = {members: key for key, members in view.items()}
        assert sorted(by_slot) == [(0,), (1,), (2, 3), (4,), (5,)], 'five groups, as built'
        assert m.canonical_stereo_group_ambiguities() == \
            (frozenset({by_slot[(1,)], by_slot[(4,)]}), frozenset({by_slot[(0,)], by_slot[(5,)]})), \
            f'creation order {order}: the reflection exchanges C2 with C5 and C1 with C6'
        assert by_slot[(2, 3)] == (SG_OR, 5), \
            'the pair is the only group of two, so its id is pinned whatever is tied around it'
        collapsed = (sorted(view.values()), m.canonical_stereo_group_ambiguities())
        if reference is None:
            reference = collapsed
        else:
            assert collapsed == reference, f'creation order {order} moved the reading'


def test_the_fixpoint_iterates_past_the_first_feedback_round():
    """The loop runs to a fixpoint and not for a fixed number of rounds, and this fixture needs two.

    The same ring with every frame parity EQUAL and OR groups {C1}, {C2}, {C4}, {C3,C5,C6}.  Equal
    parities forbid every reflection (a reflection would have to turn each 1 into a 2) and the triple
    is carried onto itself by no rotation but the identity, so the molecule has no symmetry at all and
    all four ids are pinned.

    Round 1 ties the three singletons: each carries the key (one member, one class).  Its labels
    colour the ring A A B A B B, whose only automorphism is the identity, so refining THAT colouring
    separates all six ring positions -- C4 first, with two triple neighbours where C1 and C2 have one
    each -- and round 2 gives the three singleton groups three different keys.  Measured with the
    feedback capped at one round: the three come back as one ambiguity class, which is weaker than the
    molecule allows; at two rounds the reading below is complete.
    """
    reference = None
    for order in _ring_orders(6):
        m, sids = _ring(order, (1, 1, 1, 1, 1, 1), (1, 2, 3, 4, 3, 3))
        view = _ring_view(m, sids)
        assert m.canonical_stereo_group_ambiguities() == (), \
            f'creation order {order}: nothing here is symmetric, so nothing is ambiguous'
        if reference is None:
            reference = view
            assert sorted(view.items()) == [((SG_OR, 1), (1,)), ((SG_OR, 2), (3,)),
                                            ((SG_OR, 3), (0,)), ((SG_OR, 4), (2, 4, 5))], \
                'the three singletons rank before the triple, each separated by its neighbours'
        else:
            assert view == reference, f'creation order {order} moved a pinned id'


def test_a_tie_does_not_move_the_id_of_a_group_it_does_not_involve():
    """An ambiguity has to stay inside its own class, and ranking on position alone does not.

    Octachlorocyclooctane with ALTERNATING ring-frame parities and OR groups {C1,C3,C4}, {C2,C6} and
    {C5,C7,C8}.  The rotation by four preserves an alternating pattern and carries the two triples
    onto each other while fixing the pair, so that pair of triples is a genuine F89 tie.  Nothing else
    survives: the odd rotations invert the parities, the four bond-midpoint reflections that do
    preserve them all move {C2,C6}, and a reflection through an atom is never a symmetry.  {C2,C6} is
    therefore pinned -- it is the only group of two -- and its id must be a function of the molecule.

    Measured: rank the groups by their smallest canonical position alone and the tied pair takes ids
    {1, 2} on some encodings and {1, 3} on others, because the extremal search chooses which of the
    two comes first; {C2,C6} follows it between 2 and 3, so a caller comparing the ONE key F89
    promises is pinned reads two different answers.  Ranking on the fixpoint label first, and on the
    position only inside a label, confines the choice to the tied block.
    """
    reference = None
    frame = tuple(1 + k % 2 for k in range(8))
    for order in _ring_orders(8):
        m, sids = _ring(order, frame, (1, 2, 1, 1, 3, 2, 3, 3))
        view = _ring_view(m, sids)
        by_slot = {members: key for key, members in view.items()}
        assert sorted(by_slot) == [(0, 2, 3), (1, 5), (4, 6, 7)], 'three groups, as built'
        assert m.canonical_stereo_group_ambiguities() == \
            (frozenset({by_slot[(0, 2, 3)], by_slot[(4, 6, 7)]}),), \
            'the rotation by four exchanges the two triples and nothing else'
        if reference is None:
            reference = by_slot[(1, 5)]
            assert reference == (SG_OR, 1), \
                'the label orders groups by member count first, so the pair leads'
        else:
            assert by_slot[(1, 5)] == reference, \
                f'creation order {order} moved the id of the group that is not tied'


def test_ambiguity_classes_are_ordered_by_their_smallest_canonical_id():
    """Ruling F92: the tuple `canonical_stereo_group_ambiguities()` returns compares with `==`.

    Its order has to come from the ids it reports and not from the stored bytes the caller happened to
    use, and this molecule makes that testable three times over.  Octachlorocyclooctane, alternating
    frame parities again, OR groups {C1,C5} and six singletons: the rotation by four is the only
    symmetry (every reflection it allows moves the pair, every odd rotation inverts the parities), so
    it sorts the singletons into THREE orbits of two -- {C2,C6}, {C3,C7} and {C4,C8} -- beside the
    pair, which no symmetry can exchange with a group of one.  All 5040 relabellings of the seven
    stored numbers and every encoding must return the same three frozensets in the same order.

    Measured: numbering the classes in ascending stored-byte order, the rule this replaced, returns the
    same three frozensets in a different order on 4200 of the 5040 relabellings, while
    `canonical_stereo_groups()` stays identical key for key on all 5039 comparisons -- so `==` on the
    tuple was False for two forms of one molecule whose group view was character for character the
    same, and no fixture saw it because none built more than one class.
    """
    frame = tuple(1 + k % 2 for k in range(8))
    partition = (1, 2, 3, 4, 1, 5, 6, 7)

    def check(m, sids, subject):
        view = _ring_view(m, sids)
        by_slot = {members: key for key, members in view.items()}
        assert len(view) == 7 and len(by_slot) == 7, f'{subject}: seven groups in, seven groups out'
        assert m.canonical_stereo_group_ambiguities() == \
            (frozenset({by_slot[(3,)], by_slot[(7,)]}), frozenset({by_slot[(2,)], by_slot[(6,)]}),
             frozenset({by_slot[(1,)], by_slot[(5,)]})), \
            f'{subject}: three classes, ordered by the smallest id each of them holds'
        assert by_slot[(0, 4)] == (SG_OR, 7), \
            f'{subject}: the pair is pinned by its member count, whatever is tied around it'

    for perm in permutations(range(1, 8)):      # every relabelling of the seven stored numbers
        groups = tuple(perm[g - 1] for g in partition)
        check(*_ring(range(16), frame, groups), f'stored ids {groups}')
    for order in _ring_orders(8):               # and every encoding of one labelling
        check(*_ring(order, frame, partition), f'creation order {order}')


def test_the_fixpoint_runs_past_a_third_round_when_the_molecule_needs_it():
    """The loop is not two rounds with a loop around it: this molecule needs four.

    Hexachlorocyclohexane with frame parities (1, 1, 1, 1, 1, 2) and OR groups {C1}, {C2}, {C4}, {C5}
    and {C3,C6}.  The molecule has no symmetry whatsoever -- the single odd parity leaves no rotation
    but the identity, and a reflection would have to invert five equal parities into five of the other
    value -- so all five ids are pinned.

    Measured with the feedback capped: at three rounds `canonical_stereo_group_ambiguities()` reports
    {C1},{C4} as one class and {C2},{C5} as another, two ambiguities the molecule does not have; the
    full fixpoint pins all five.  The round histogram over all 12,992 (frame parity pattern, group
    partition) fixtures of this ring is 7,472 settled by round 1, 1,152 by round 2, 4,248 by round 3
    and 120 needing a fourth, against the hard bound of n rounds the loop carries.

    The subjects differ only in the stored numbers and the encoding, and the assertion is the full
    id -> members map, so a round lost anywhere in the loop shows up here.
    """
    reference = [None]
    partition = (1, 2, 3, 4, 5, 3)

    def check(m, sids, subject):
        view = _ring_view(m, sids)
        assert m.canonical_stereo_group_ambiguities() == (), \
            f'{subject}: the fourth round pins all five groups'
        if reference[0] is None:
            reference[0] = view
            assert sorted(view.items()) == [((SG_OR, 1), (0,)), ((SG_OR, 2), (3,)),
                                            ((SG_OR, 3), (1,)), ((SG_OR, 4), (4,)),
                                            ((SG_OR, 5), (2, 5))], \
                'four singletons, then the pair, each separated by its surroundings'
        else:
            assert view == reference[0], f'{subject} moved a pinned id'

    for perm in permutations(range(1, 6)):      # every relabelling of the five stored numbers
        groups = tuple(perm[g - 1] for g in partition)
        check(*_ring(tuple(range(12)), (1, 1, 1, 1, 1, 2), groups), f'stored ids {groups}')
    for order in _ring_orders(6):               # and every encoding of one labelling
        check(*_ring(order, (1, 1, 1, 1, 1, 2), partition), f'creation order {order}')


def test_a_forged_unspecified_byte_reads_the_same_in_both_views():
    """The canonical view renumbers OR and AND and passes every other kind through unchanged.

    `set_stereo_group` forces the number to 0 for kinds 0 and 1, so a nonzero number under kind 0
    can only arrive from a forged buffer -- `from_bytes` checks the segment's LENGTH and not its
    bytes.  When it does, `canonical_stereo_groups()` must report the same key `stereo_groups()`
    does: the view's job is to renumber the two numbered kinds, and reporting (0, 0) for a stored
    (0, 5) would be the view inventing a normalisation of a byte it does not own (ruling F79).
    """
    m = MoleculeContainer()
    with m.edit():
        sids = [m.add_atom('C', implicit_h=1), m.add_atom('C', implicit_h=3),
                m.add_atom('Cl'), m.add_atom('F'), m.add_atom('Br')]
        for s in sids[1:]:
            m.add_bond(sids[0], s, 1)
        m.set_parity(sids[0], 1)
        m.set_stereo_group(sids[0], SG_OR, 37)          # 0x80 | 37 == 165
    buffer = bytearray(m.to_bytes())
    assert buffer.count(165) == 1, 'forged buffer: the group byte has to be the only 165 in it'
    buffer[buffer.index(165)] = 5                       # kind 0, number 5: unreachable via the API
    forged = MoleculeContainer.from_bytes(bytes(buffer))
    assert forged.stereo_groups() == {(SG_UNSPECIFIED, 5): [sids[0]]}, 'the byte survived the trip'
    assert forged.canonical_stereo_groups() == {(SG_UNSPECIFIED, 5): [sids[0]]}, \
        'so the canonical view has to agree with it'
    assert forged.canonical_stereo_group_ambiguities() == (), 'and nothing is ambiguous here'
    assert m.canonical_stereo_groups() == {(SG_OR, 1): [sids[0]]}, \
        'while a real OR group IS renumbered, so the pass-through is not just inaction'


# --- the word IV screen -------------------------------------------------------------------------

def test_a_stereo_query_is_screened_out_before_the_unit_table_is_built():
    """The bit is 'a parity is configured', not 'the parity is odd'.

    `query_may_match` runs before matcher_init opens the stereo door, so a query that demands a
    configured centre and a target that has none never pays for the unit table -- `total_len` is
    what notices, exactly as the door test above uses it.  The second half is why the bit is not
    the stereo VALUE bit: parity 1 is EVEN, so the value bit is clear on this target, and a screen
    built on that bit would refuse a molecule that matches.
    """
    bare, _ = _chiral_target(0)
    q, _ = _chiral_query(1)
    before = bare.total_len
    assert not q.is_substructure(bare), 'nothing to match: no centre is configured'
    assert bare.total_len == before, 'and the signature screen said so before the door opened'

    configured, _ = _chiral_target(1)
    before = configured.total_len
    assert q.is_substructure(configured), 'an EVEN parity is configured and must not be screened'
    assert configured.total_len > before, 'so this one does open the door'


# 2,3-dichlorobutane: C1 a methyl, C2 and C3 the stereocentres, C4 a methyl, and a chlorine on each
# centre.  The two centres are ONE constitutional class, so nothing about them differs except their
# stereochemistry -- and the frame that decides what "differs" means is (the methyl, the other centre,
# the chlorine, the implicit hydrogen), which asks the same question at both of them and in every
# encoding, where the stored byte does not (ruling F26, argued at the ring fixtures above).
#
# The swap of the two centres, (C1 C4)(C2 C3)(Cl Cl), carries that frame onto that frame IN ORDER, so
# it preserves both parities exactly when the two are EQUAL.  Equal frame parities are therefore the
# C2-symmetric diastereomer, whose two centres are homotopic and whose two stereo groups nothing can
# tell apart; opposite parities are the meso diastereomer, whose centres are R and S and are
# separable by any rule that reads the parity at all.
_BUTANE_ATOMS = (('C', 3), ('C', 1), ('C', 1), ('C', 3), ('Cl', 0), ('Cl', 0))
_BUTANE_BONDS = ((0, 1), (1, 2), (2, 3), (1, 4), (2, 5))
_BUTANE_FRAMES = ((0, 2, 4), (3, 1, 5))     # C2's refs then C3's, both (methyl, centre, chlorine)


def _butane(order, frame, groups=(1, 2), kind=SG_OR):
    """2,3-dichlorobutane with one stereo group per centre, in one creation `order`.

    `frame` is each centre's parity in the frame above, so the molecule is a function of `frame`
    alone; `order` chooses only how it is encoded, and `groups` only how it is numbered.
    """
    m = MoleculeContainer()
    with m.edit():
        sids = [0] * 6
        for slot in order:
            element, h = _BUTANE_ATOMS[slot]
            sids[slot] = m.add_atom(element, implicit_h=h)
        for i, j in _BUTANE_BONDS:
            m.add_bond(sids[i], sids[j], 1)
        for centre, number in zip((1, 2), groups):
            m.set_parity(sids[centre], 1)
            m.set_stereo_group(sids[centre], kind, number)
    assert m.unit_of(sids[1]) is not None and m.unit_of(sids[2]) is not None
    flip = [(c, 3 - m.parity_of(sids[c])) for c in (1, 2)
            if _butane_parity(m, sids, c) != frame[c - 1]]
    with m.edit():
        for c, parity in flip:
            m.set_parity(sids[c], parity)
    assert tuple(_butane_parity(m, sids, c) for c in (1, 2)) == tuple(frame), 're-encoding failed'
    return m, sids


def _butane_parity(m, sids, centre):
    """Centre C2 (1) or C3 (2) read in the frame (methyl, the other centre, chlorine, hydrogen)."""
    return m.translate_stereo(sids[centre],
                              tuple(sids[r] for r in _BUTANE_FRAMES[centre - 1]) + (None,))


def test_the_canonical_view_ignores_the_stored_id_where_nothing_else_separates_the_groups():
    """Ruling F80 step 1: the seed carries the group KIND and never the stored group id.

    The C2-symmetric diastereomer, one OR group per centre.  The two centres are interchangeable --
    one constitutional class, and the swap preserves both frame parities because they are equal -- so
    the stored ids are the ONLY thing that tells the two groups apart, and a rule that read them would
    have nothing to contradict it.  That is what makes this fixture, and not the id-permutation tests
    above, the one that measures the seed: on the meso form the parity separates the centres and the
    id in the seed would change nothing.  Measured -- with the stored id folded into the seed this test
    fails and every other test in the core suite passes.

    A view built on the stored ids says group 1 is the first centre in one molecule and the second in
    the other, which is a distinction the caller's choice of ids invented.
    """
    a, sids = _butane(tuple(range(6)), (1, 1), groups=(1, 2))
    b, _ = _butane(tuple(range(6)), (1, 1), groups=(2, 1))
    va, vb = a.canonical_stereo_groups(), b.canonical_stereo_groups()
    assert sorted(va) == [(SG_OR, 1), (SG_OR, 2)], 'two OR groups, densely numbered'
    assert sorted(v for members in va.values() for v in members) == sorted(sids[1:3])
    assert va == vb, 'the stored id is the only difference between the two molecules'


def test_the_canonical_view_pairs_a_group_with_the_same_parity_under_every_relabelling():
    """Ruling F80 step 1, the other half: the seed carries the PARITY, so the order it feeds is
    label-invariant on a molecule whose constitution alone cannot separate its centres.

    The meso diastereomer, whose two centres are one constitutional class but opposite in the frame.
    An unseeded canonical order ties them and breaks the tie by arena slot -- measured: over the 720
    encodings below `canonical_order` puts the frame-parity-1 centre at position 1 in some and at
    position 5 in others.  The parity in the seed splits the class, and then which group is numbered 1
    stops depending on the caller's atom order.

    The observable has to say what SITS in each group and not merely which positions the groups
    occupy: ids are handed out BY ascending member position, so a comparison of positions alone reads
    back {1: [smaller], 2: [larger]} whatever the order did, and would pass against an unseeded call.
    The member's parity is the thing that moves, and it is read in the FRAME: pairing an id with the
    stored byte would compare two encodings on the one thing that legitimately differs between them.
    Measured -- drop the parity term from the seed entirely, the coarse landing ruling F95 allows, and
    this test fails on 360 of the 720 encodings.
    """
    def view_by_parity(order):
        m, sids = _butane(order, (1, 2))
        return {key: sorted(_butane_parity(m, sids, 1 if v == sids[1] else 2) for v in members)
                for key, members in m.canonical_stereo_groups().items()}

    reference = view_by_parity(tuple(range(6)))
    assert sorted(reference) == [(SG_OR, 1), (SG_OR, 2)], 'two groups to tell apart'
    assert sorted(v for members in reference.values() for v in members) == [1, 2], 'and two parities'
    for order in permutations(range(6)):
        assert view_by_parity(order) == reference, f'creation order {order} moved the pairing'


def test_an_and_pair_off_a_ring_is_reported_as_one_ambiguity_class_in_every_encoding():
    """Ruling F89's second half again, on the two axes the ring fixtures never reach: AND, and no ring.

    The C2-symmetric diastereomer with one AND group per centre.  The swap of the two centres
    preserves both frame parities, so it is an automorphism of the annotated molecule and the two
    groups are interchangeable: the view must report one class covering both keys and must still
    report TWO groups, because an AND pair of one member each is not one AND group of two.

    Both halves of the fixture are deliberate.  The kind is a term of the seed in its own right, and
    every other ambiguity fixture in this file is OR.  And a unit's directions come from a methyl and a
    chlorine here rather than from two ring neighbours, so the parity term of the seed is read in a
    frame no ring supplies.  The second loop is what stops the first from passing vacuously: give the
    same skeleton OPPOSITE frame parities and the swap can no longer preserve them, so the tie has to
    disappear and every id has to pin -- if it did not, the first loop would be measuring the skeleton
    rather than the stereochemistry.

    This is also the file's one measurement of what a surviving tie actually costs, because it is the
    only ambiguity fixture whose tied members are told apart by something other than the labelling:
    the raw id -> WHICH CENTRE direction takes two values over the 720 encodings -- both assignments
    occur -- while the collapsed reading asserted below takes one.  The opposite-parity half takes one
    value in the raw direction too, which is the same statement as its empty ambiguity tuple.
    """
    reference = None
    for order in permutations(range(6)):
        m, sids = _butane(order, (1, 1), kind=SG_AND)
        view = m.canonical_stereo_groups()
        ambiguities = m.canonical_stereo_group_ambiguities()
        assert sorted(view) == [(SG_AND, 1), (SG_AND, 2)], 'two AND groups, densely numbered'
        assert sorted(len(v) for v in view.values()) == [1, 1], 'and neither of them absorbed the other'
        assert ambiguities == (frozenset(view),), \
            f'creation order {order}: the swap of the two centres is an automorphism, so they tie'
        collapsed = (sorted(sorted(_butane_parity(m, sids, 1 if v == sids[1] else 2) for v in members)
                            for members in view.values()), sorted(ambiguities[0]))
        if reference is None:
            reference = collapsed
        else:
            assert collapsed == reference, f'creation order {order} moved the collapsed reading'

    pinned = None
    for order in permutations(range(6)):
        m, sids = _butane(order, (1, 2), kind=SG_AND)
        assert m.canonical_stereo_group_ambiguities() == (), \
            f'creation order {order}: opposite parities leave the swap no way to preserve them'
        view = {key: sorted(_butane_parity(m, sids, 1 if v == sids[1] else 2) for v in members)
                for key, members in m.canonical_stereo_groups().items()}
        if pinned is None:
            pinned = view
        else:
            assert view == pinned, f'creation order {order} moved a pinned id'


# --- ruling F95 on the three NON-TETRAHEDRAL kinds ----------------------------------------------
# Every fixture above this line is tetrahedral, and the seed's parity term reaches the bond kinds by
# different arithmetic: ruling F56 makes the wholesale exchange of the two direction PAIRS even, so
# `_frame_free_parity_code` sorts each pair on its own and flips one bit per reversed pair, where
# SU_TETRA sorts all four slots and translates the permutation through the parity table.  That branch
# was never measured on the encoding axis, and the encoding axis is where ruling F95's defect lived.
#
# The skeleton is TWO IDENTICAL DISCONNECTED COMPONENTS with one stereo group each, which is what puts
# the ground truth beyond the implementation's reach.  The component swap is a constitutional
# automorphism; each component's stereochemistry is stated in the same LOGICAL frame, so the swap
# carries that frame onto that frame in order and therefore preserves both parities exactly when the
# two are EQUAL.  Equal frame parities are two interchangeable groups; opposite parities are two groups
# that any rule reading the parity at all must pin, and that half is what stops the first from passing
# on the skeleton alone.
#
# Each kind is carried by a molecule with FOUR named directions and not two, and the obvious smaller
# fixtures are the reason.  Hold one frame parity fixed on 2-butene (SU_CIS_TRANS) or penta-2,3-diene
# (SU_ALLENE) and the stored byte is 1 in EVERY creation order -- measured over all 24 and all 120 of
# them.  Each of their direction pairs is (a carbon, an unnamed direction), and an unnamed direction's
# key sits below every atom's, so the pair sorts the same way always and there is no encoding axis to
# sweep.  Give each pair a second named direction and the sort order becomes a function of the creation
# order: 1,2-dichloro-1,2-difluoroethene stores both bytes over its 720 orders and
# 1,3-dichloro-1,3-difluoroallene over its 5,040.  Over the orders swept below each fixture stores
# THREE of the four byte PATTERNS, and two of the three occur in the tied and the pinned half alike
# (asserted, not assumed) -- so the stored byte does not determine even whether the molecule is
# symmetric.
SU_CIS_TRANS, SU_ALLENE, SU_ATROPISOMER = 1, 2, 3

# (atoms as (element, implicit_h); bonds as (i, j, order); the atoms that may anchor the unit; the
#  unit's four directions in the order the frame names them; the atom carrying the stereo group; the
#  kind `stereo_units()` must report; the name).  Hydrogen counts are stated because the core derives
#  none.  Two atoms may anchor a cis/trans or an atropisomer unit, because ruling F45 lets a bond kind
#  anchor at either end and leaves the choice to slot order.
_DIFLUOROETHENE = ((('C', 0), ('C', 0), ('F', 0), ('Cl', 0), ('F', 0), ('Cl', 0)),
                   ((0, 1, 2), (0, 2, 1), (0, 3, 1), (1, 4, 1), (1, 5, 1)),
                   (0, 1), (2, 3, 4, 5), 0, SU_CIS_TRANS, '1,2-dichloro-1,2-difluoroethene')
_DIFLUOROALLENE = ((('C', 0), ('C', 0), ('C', 0), ('F', 0), ('Cl', 0), ('F', 0), ('Cl', 0)),
                   ((0, 1, 2), (1, 2, 2), (0, 3, 1), (0, 4, 1), (2, 5, 1), (2, 6, 1)),
                   (1,), (3, 4, 5, 6), 1, SU_ALLENE, '1,3-dichloro-1,3-difluoroallene')
_HALOBIPHENYL = ((('C', 0), ('C', 0), ('C', 1), ('C', 1), ('C', 1), ('C', 1),
                  ('C', 0), ('C', 0), ('C', 1), ('C', 1), ('C', 1), ('C', 1), ('Cl', 0), ('F', 0)),
                 ((0, 1, 2), (1, 2, 1), (2, 3, 2), (3, 4, 1), (4, 5, 2), (5, 0, 1),
                  (6, 7, 2), (7, 8, 1), (8, 9, 2), (9, 10, 1), (10, 11, 2), (11, 6, 1),
                  (0, 6, 1), (1, 12, 1), (7, 13, 1)),
                 (0, 6), (1, 5, 7, 11), 0, SU_ATROPISOMER, "2-chloro-2'-fluorobiphenyl")
_KIND_SPECS = (_DIFLUOROETHENE, _DIFLUOROALLENE, _HALOBIPHENYL)


def _two_components(spec, order, parities, kind=SG_OR):
    """`spec` doubled, written in creation `order`, one stereo group per component.

    `parities[c]` is component c's parity in the frame `spec` names, so the molecule is a function of
    `parities` alone and `order` chooses only how it is encoded.  Returns the molecule, the logical
    slot -> atom id table over both components, and the two anchors in component order.
    """
    atoms, bonds, owners, _, group_atom, expected, name = spec
    n = len(atoms)
    m = MoleculeContainer()
    sids = [0] * (2 * n)
    with m.edit():
        for slot in order:
            element, h = atoms[slot % n]
            sids[slot] = m.add_atom(element, implicit_h=h)
        for shift in (0, n):
            for i, j, o in bonds:
                m.add_bond(sids[i + shift], sids[j + shift], o)
        for number, shift in enumerate((0, n), 1):
            m.set_stereo_group(sids[group_atom + shift], kind, number)
    # which end anchors is a choice ruling F45 leaves to slot order, so the anchor is DISCOVERED from
    # perception rather than assumed; the frame is stated on the refs, which no such choice touches
    slot = {s: k for k, s in enumerate(sids)}
    anchors = [None, None]
    units = m.stereo_units()
    assert len(units) == 2, f'{name}: {len(units)} units perceived, so the fixture is not two of one'
    for unit in units:
        k = slot[unit['anchor']]
        assert k % n in owners, f'{name}: unit anchored at logical {k}, outside {owners}'
        assert anchors[k // n] is None, f'{name}: both units landed in one component'
        assert unit['kind'] == expected, f"{name}: kind {unit['kind']} perceived, not {expected}"
        anchors[k // n] = unit['anchor']
    with m.edit():
        for a in anchors:
            m.set_parity(a, 1)
    # the flip list is built BEFORE the edit opens: translate_stereo refuses pending edits, and one
    # component's byte cannot change the other's frame -- they are not even connected
    flip = [(a, 3 - m.parity_of(a)) for a, want, shift in zip(anchors, parities, (0, n))
            if _component_parity(m, sids, spec, a, shift) != want]
    with m.edit():
        for a, parity in flip:
            m.set_parity(a, parity)
    assert _component_parities(m, sids, spec, anchors) == tuple(parities), 're-encoding failed'
    return m, sids, anchors


def _component_parity(m, sids, spec, anchor, shift):
    """One component's parity in the frame `spec` names, which asks the same question in every
    encoding -- unlike `parity_of`, whose byte ruling F26 makes a function of the creation order."""
    return m.translate_stereo(anchor, tuple(sids[r + shift] for r in spec[3]))


def _component_parities(m, sids, spec, anchors):
    n = len(spec[0])
    return tuple(_component_parity(m, sids, spec, anchors[c], c * n) for c in (0, 1))


def _kind_view(m, sids, spec, anchors):
    """(kind, id) -> the FRAME PARITIES of the members, the labelling-free observable.

    Member ids cannot be compared across two creation orders -- they ARE the relabelling -- and the
    stored bytes must not be compared either, since they are the one thing that legitimately differs
    between two encodings of one molecule.  What each group HOLDS is a fact about the molecule.
    """
    n = len(spec[0])
    parities = _component_parities(m, sids, spec, anchors)
    slot = {s: k for k, s in enumerate(sids)}
    return {key: sorted(parities[slot[a] // n] for a in members)
            for key, members in m.canonical_stereo_groups().items()}


def _kind_orders(n):
    """Creation orders that vary the encoding and nothing else: every cyclic rotation of the n atom
    slots and every rotation reversed, 2n orders -- 24 for the ethene, 28 for the allene, 56 for the
    biphenyl.  A rotation moves the slot distance between a unit's anchor and each of its refs through
    every value, which is what ruling F26's byte is a function of, and a reversal inverts every such
    distance as well as swapping the two components' turn to be created.
    """
    rotations = [tuple((start + k) % n for k in range(n)) for start in range(n)]
    return rotations + [o[::-1] for o in rotations]


@pytest.mark.parametrize('spec', _KIND_SPECS, ids=[s[6] for s in _KIND_SPECS])
def test_a_non_tetrahedral_kind_reads_the_same_in_every_encoding(spec):
    """Ruling F95 for SU_CIS_TRANS, SU_ALLENE and SU_ATROPISOMER: the reading is a function of the
    molecule, and the stored parity byte is not.

    Two copies of the fixture, one stereo group each, in both numbered kinds and in both halves of the
    ground truth above.  Opposite frame parities: the component swap cannot preserve them, so the two
    groups are distinguishable and every id must be pinned AND must sit on the same parity in every
    encoding.  Equal frame parities: the swap is an automorphism of the annotated molecule, so the two
    groups are interchangeable and the report must be one ambiguity class covering both keys -- while
    still reporting TWO groups of one member each, because two tied groups are not one group of two.

    The premise is asserted rather than assumed, twice over.  Each half stores three distinct byte
    patterns over the orders swept, so an implementation seeding its refinement on the byte cannot pass
    by luck; and the two halves SHARE byte patterns -- (1, 1) and (1, 2) occur in both -- so no rule
    that reads the byte can even tell the symmetric molecule from the pinned one here.

    Measured against the mutant that gives the bond kinds their stored byte back -- one line in
    `_frame_free_parity_code`, `return parity + 1` in place of the pair-sorted code -- both halves fail
    on all three fixtures, and NOTHING ELSE IN THE CORE SUITE DOES: 3 failed, 857 passed.  The pinned
    half moves its id -> parity pairing on 4 of 24, 4 of 28 and 16 of 56 encodings and reports a tie
    where the molecule is pinned on 8, 8 and 32 of them; the tied half reports every id PINNED on the
    same 8, 8 and 32, which is the unsafe direction rule 3 forbids -- two encodings of one molecule both
    claiming to be pinned and naming different ids, which is the whole of ruling F95.
    """
    shared = None
    for kind in (SG_OR, SG_AND):
        halves = {}
        for parities, truth in (((1, 2), 'pinned'), ((1, 1), 'tied')):
            reference = None
            bytes_seen = set()
            for order in _kind_orders(2 * len(spec[0])):
                m, sids, anchors = _two_components(spec, order, parities, kind)
                bytes_seen.add(tuple(m.parity_of(a) for a in anchors))
                view = _kind_view(m, sids, spec, anchors)
                ambiguities = m.canonical_stereo_group_ambiguities()
                assert sorted(view) == [(kind, 1), (kind, 2)], 'two groups, densely numbered'
                assert sorted(len(v) for v in view.values()) == [1, 1], \
                    'and neither of them absorbed the other'
                if truth == 'pinned':
                    assert ambiguities == (), \
                        f'order {order}: opposite parities leave the swap nothing to preserve'
                    assert sorted(view.values()) == [[1], [2]], 'one parity each, and opposite'
                else:
                    assert ambiguities == (frozenset(view),), \
                        f'order {order}: the component swap is an automorphism, so the groups tie'
                if reference is None:
                    reference = view
                else:
                    assert view == reference, f'order {order} moved the id -> parity pairing'
            assert len(bytes_seen) == 3, \
                f'the sweep must vary the ENCODING, and this half stored {sorted(bytes_seen)}'
            halves[truth] = bytes_seen
        common = halves['pinned'] & halves['tied']
        assert common == {(1, 1), (1, 2)}, \
            f'one byte pattern must serve both a tied and a pinned molecule, not {sorted(common)}'
        assert shared is None or shared == common, 'and the same ones under either group kind'
        shared = common
