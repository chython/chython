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
"""Order 4 is a STORED order, and this file is the evidence rather than the argument.

The arena stores order 4 rather than making a parser kekulise on the way in, on the input-fidelity
invariant -- a molecule must come back out as the file wrote it, and a file that drew an aromatic
ring said something a Kekule ring does not. Two spellings of one chemical object are therefore both
storable, they are DIFFERENT objects to every derived answer, and `kekule()` / `thiele()` are the only
crossings between them.

What that decision costs is an audit of every reader of a bond order, and the outcome of the audit is
the interesting part: **nothing in the core refuses an aromatic molecule.** Hybridization answers 4.
The feature words separate order 4 from a dative bond. The canonical bond word already folded the
flag in. The stereo surface reads degree, hydrogen count and ring membership and never an order --
an axis rule that read a ring double bond instead would lose every aromatic biaryl -- and takes one
arithmetic correction to an electron budget. The tests below are organised by that audit:
representation, then serialisation, then each consumer asked whether it can answer.
"""
from itertools import permutations

from pytest import raises

from chython.core import MoleculeContainer, QueryContainer


def ring(elements, orders):
    m = MoleculeContainer()
    ids = [m.add_atom(e) for e in elements]
    for i, o in enumerate(orders):
        m.add_bond(ids[i], ids[(i + 1) % len(ids)], o)
    return m, ids


def benzene(order=4):
    return ring([6] * 6, [order] * 6)


def kekule_benzene():
    return ring([6] * 6, [2, 1, 2, 1, 2, 1])


def kinds(m):
    """The multiset of stereo-unit kinds, as a sorted list.

    `stereo_units()` yields DICTS (`kind`, `anchor`, `refs`, ...), not tuples, so every comparison
    here goes through a named key.
    """
    return sorted(u['kind'] for u in m.stereo_units())


def anchors(m, kind):
    return {u['anchor'] for u in m.stereo_units() if u['kind'] == kind}


def canonical_string(m):
    """A labelling-independent encoding of a molecule, keyed by canonical POSITION.

    The same idiom `test_canonical.py` uses, and the reason it is not written inline: `canonical_order()`
    returns a dict from stable id to position, so iterating it yields stable ids in slot order and
    indexing it with a position is a lookup by the wrong key. Both mistakes give a plausible tuple that
    varies with creation order, which is exactly the failure an invariance test is trying to detect --
    so it must not be able to manufacture one.
    """
    position = m.canonical_order()
    elements = [None] * len(position)
    for sid, p in position.items():
        elements[p] = m.element_of(sid)
    edges = sorted((min(position[a], position[b]), max(position[a], position[b]), m.order_of(a, b))
                   for a in m.atom_numbers for b in m.neighbors_of(a))
    return repr((elements, edges))


def biphenyl(order=4, substituted=False):
    """Two six-rings joined by one acyclic single bond, each ring written with `order`.

    `substituted` puts a methyl ortho to each pivot -- a 2,2'-disubstituted biaryl, which is what
    an atropisomer axis needs distinguishable ortho pairs for.
    """
    m = MoleculeContainer()
    a = [m.add_atom(6) for _ in range(12)]
    for base in (0, 6):
        for i in range(6):
            m.add_bond(a[base + i], a[base + (i + 1) % 6], order)
    m.add_bond(a[0], a[6], 1)
    if substituted:
        for pivot_ortho in (1, 7):
            m.add_bond(a[pivot_ortho], m.add_atom(6), 1)
    # THE COUNTS ARE STATED, and for this fixture that is not optional.  `add_atom` stores H_UNKNOWN
    # when nothing is said, and axis detection reads the hydrogen count -- a pivot with three heavy
    # neighbours and an unknown hydrogen may or may not have a fourth direction, and refusing is the
    # only safe answer.  Leaving them unsaid does not weaken this fixture, it EMPTIES it: every unit
    # list comes back `[]` and four tests here stop asserting on an axis.  Degree decides: a methyl
    # carries 3, an aromatic CH 1, a substituted ring carbon 0.
    for s in m.atom_numbers:
        m.set_hydrogens(s, 3 if m.degree_of(s) == 1 else 3 - m.degree_of(s))
    return m, a


# ── representation ──────────────────────────────────────────────────────────────────────


def test_hybridization_of_an_aromatic_ring_carbon_is_four_and_its_kekule_twins_is_two():
    """The headline claim, on both spellings, because the point is that they DIFFER.

    `derive_scalars` tests for an aromatic bond before it counts double bonds, so an aromatic ring
    carbon reports 4 -- V2's `hybridization == 4` -- while the same carbon written Kekule reports 2
    from its one double bond. Nothing derives the first from the second: aromaticity is stored, not
    perceived, and a molecule that did not say it was aromatic is not told that it is.
    """
    arom, ids = benzene()
    for s in ids:
        assert arom.hybridization_of(s) == 4

    kek, kids = kekule_benzene()
    for s in kids:
        assert kek.hybridization_of(s) == 2


def test_hybridization_four_needs_only_one_aromatic_bond_and_outranks_the_double_bonds():
    """Which branch wins when an atom has both, stated as a test rather than left to the reader.

    An atom carrying one aromatic bond and one exocyclic double bond is sp2 either way, so the
    interesting case is the one where the two branches disagree: two cumulated double bonds report
    5 (chython 3's own value for an allene or a sulfone centre) and an aromatic bond reports 4. The
    aromatic test runs first, so 4 wins -- a delocalised atom's geometry is settled by the ring it
    is part of, and 5 would claim two localised pi systems it does not have.
    """
    m = MoleculeContainer()
    c = [m.add_atom(6) for _ in range(4)]
    m.add_bond(c[0], c[1], 4)       # one aromatic bond, in no ring at all
    m.add_bond(c[0], c[2], 2)
    m.add_bond(c[0], c[3], 2)       # ... plus two double bonds: 5 without the aromatic branch
    assert m.hybridization_of(c[0]) == 4
    assert m.hybridization_of(c[2]) == 2


def test_the_aromatic_bond_count_is_exact_and_is_kekule_is_exactly_its_zero():
    arom, _ = benzene()
    assert arom.aromatic_bond_count == 6
    assert not arom.is_kekule

    kek, _ = kekule_benzene()
    assert kek.aromatic_bond_count == 0
    assert kek.is_kekule

    empty = MoleculeContainer()
    assert empty.aromatic_bond_count == 0
    assert empty.is_kekule, 'a molecule with no bonds is trivially Kekule'


def test_a_molecule_may_hold_one_aromatic_ring_beside_one_alternating_ring():
    """The case that rules out a three-state KEKULE/AROMATIC/MIXED flag.

    Two rings joined by a single bond, one written aromatic and one written Kekule -- two
    differently-drawn inputs in one record, which a reaction file produces routinely. A count
    describes this honestly; an enum would have to call it "mixed" and then somebody would have to
    decide whether the alternating ring OUGHT to be aromatic, which is a perception question no
    stored state can answer.
    """
    m, a = biphenyl(order=4)
    for i in range(6):                       # rewrite the second ring as Kekule
        m.set_order(a[6 + i], a[6 + (i + 1) % 6], 2 if i % 2 == 0 else 1)
    assert m.aromatic_bond_count == 6
    assert not m.is_kekule
    assert {m.hybridization_of(s) for s in a[:6]} == {4}
    assert {m.hybridization_of(s) for s in a[6:]} == {2}


def test_the_count_follows_every_operation_that_can_change_a_bond():
    """Recomputed from the surviving edges, never incremented -- which is what makes it exact.

    Each step below is a different path into the arena: an in-place edit, a scoped compaction, a
    clone, a renumbering and a deletion. A counter maintained incrementally would survive most of
    them and be wrong on one; there is no path that reports a bond the graph does not have.
    """
    m, ids = benzene()
    m.set_order(ids[0], ids[1], 1)                      # in place
    assert m.aromatic_bond_count == 5
    with m.edit():                                      # through a scope
        m.set_order(ids[1], ids[2], 2)
    assert m.aromatic_bond_count == 4
    assert m.copy().aromatic_bond_count == 4            # clone carries it
    other = m.copy()
    other.remap({s: s + 100 for s in other.atom_numbers})
    assert other.aromatic_bond_count == 4               # remap does not touch orders
    m.delete_atom(ids[3])                               # loses the two bonds at that atom
    assert m.aromatic_bond_count == 2
    m.delete_bond(ids[4], ids[5])
    assert m.aromatic_bond_count == 1


# ── serialisation ───────────────────────────────────────────────────────────────────────


def test_an_aromatic_molecule_round_trips_byte_identically():
    arom, _ = benzene()
    raw = arom.to_bytes()
    back = MoleculeContainer.from_bytes(raw)
    assert back.to_bytes() == raw
    assert back.aromatic_bond_count == 6
    assert [back.order_of(a, b) for a, b in ((1, 2), (2, 3))] == [4, 4]
    # `pack()`/`unpack`, the legacy pach record, keeps order 4 as well -- the aromatic flag is a bond
    # fact and has to survive every serialisation.
    assert not arom.is_kekule and back.is_kekule == arom.is_kekule


def test_a_buffer_whose_order_and_aromatic_flag_disagree_is_refused():
    """The pair is ONE fact written twice, so `from_bytes` requires them to agree.

    Two places read it -- an order switch reads `order == 4`, a topology switch reads the flag --
    and a buffer that set one without the other would answer one question aromatic and the other
    Kekule. Both directions are forged here because they fail for different reasons: an order-4
    half-edge without the flag loses its feature bit, a flagged order-1 half-edge gains one.
    """
    arom = benzene()[0].to_bytes()
    off, flag_bit = _halfedge_offsets(arom)[0], 2
    stripped = bytearray(arom)
    stripped[off + 6] &= ~flag_bit                       # order 4, flag cleared
    with raises(ValueError, match='HE_AROMATIC'):
        MoleculeContainer.from_bytes(bytes(stripped))

    kek = kekule_benzene()[0].to_bytes()
    flagged = bytearray(kek)
    flagged[_halfedge_offsets(kek)[0] + 6] |= flag_bit   # order 1 or 2, flag set
    with raises(ValueError, match='HE_AROMATIC'):
        MoleculeContainer.from_bytes(bytes(flagged))


def test_a_reserved_halfedge_flag_bit_is_refused():
    """Every bit above the defined mask is reserved, and reserved means rejected rather than ignored.

    The flags field is 16 bits. Five are defined -- in_ring, aromatic, and three of CIP code -- and
    silently masking the rest would make a future flag unversioned: a v4 reader would accept a v5
    buffer and answer as though the flag were absent, which is exactly the failure the version field
    exists to prevent.

    THE MASK MOVED ONCE ALREADY, when the CIP code took bits 2-4, and this test is why that was safe
    to do: it is written against the mask rather than against a bit number, so widening the defined
    region moves which bits it probes and does not weaken what it proves.
    """
    raw = benzene()[0].to_bytes()
    for bit in (0x20, 0x40, 0x80):
        forged = bytearray(raw)
        forged[_halfedge_offsets(raw)[0] + 6] |= bit
        with raises(ValueError, match='reserved'):
            MoleculeContainer.from_bytes(bytes(forged))
    forged = bytearray(raw)
    forged[_halfedge_offsets(raw)[0] + 7] |= 1           # the high half of the same field
    with raises(ValueError, match='reserved'):
        MoleculeContainer.from_bytes(bytes(forged))


def test_a_halfedge_cip_code_outside_the_domain_is_refused():
    """The field holds three bits and four descriptors are defined, so 5, 6 and 7 are reachable.

    A width check would pass all three. They are refused by a DOMAIN check, for the same reason the
    implicit-hydrogen nibble is bounded at 14 rather than at 15: the width of a field is not the
    domain of what it holds, and a code with no name would surface later as an index error in Python
    rather than as a rejected buffer here.
    """
    raw = benzene()[0].to_bytes()
    for code in (5, 6, 7):
        forged = bytearray(raw)
        forged[_halfedge_offsets(raw)[0] + 6] |= code << 2
        with raises(ValueError, match='CIP code'):
            MoleculeContainer.from_bytes(bytes(forged))


def _halfedge_offsets(data):
    """Byte offset of each half-edge record, read out of the segment table by hand."""
    from struct import unpack_from
    off, length = unpack_from('<II', data, 24 + 8 * 2)   # table entry 2 is the CSR edge array
    return [off + 8 * i for i in range(length // 8)]


# ── substructure search: the consumer a gate would have broken ──────────────────────────


def _bond_query(order_name, value=0, negated=False):
    """A two-carbon query whose bond carries one primitive."""
    q = QueryContainer()
    a, b = q.add_atom(), q.add_atom()
    for sid in (a, b):
        q.atom_primitive(sid, 'element', 6)
    q.add_bond(a, b)
    q.bond_primitive(a, b, order_name, value, negated=negated)
    return q


def test_an_aromatic_bond_is_found_by_an_aromatic_query_and_by_nothing_else():
    """The whole reason the feature words were fixed instead of gated.

    Gating would have made substructure search refuse on aromatic molecules, which is most
    molecules. The fix costs no bit: order 4 shares the order-8 bit with a dative bond and is
    separated from it by the aromatic topology bit, which fires exactly when the order is 4.
    """
    arom, _ = benzene()
    assert _bond_query('bond_aromatic').count(arom) == 12, 'six bonds, both directions'
    for order in (1, 2, 3, 8):
        assert _bond_query('bond_order', order).count(arom) == 0, \
            'order %d must not match an aromatic bond' % order


def test_a_coordination_bond_query_does_not_match_an_aromatic_bond_and_the_converse():
    """The one place the shared bit could have gone wrong, in both directions.

    Order 4 and order 8 occupy the same order bit, so a demand for one must forbid the other's
    topology bit. `bond_order 8` additionally forbids the aromatic bit; `bond_aromatic` demands it.
    """
    arom, _ = benzene()
    dative = ring([6] * 6, [8] * 6)[0]
    assert _bond_query('bond_order', 8).count(dative) == 12
    assert _bond_query('bond_order', 8).count(arom) == 0
    assert _bond_query('bond_aromatic').count(dative) == 0
    assert _bond_query('bond_aromatic').count(arom) == 12


def test_an_aromatic_query_does_not_match_a_kekule_molecule_or_the_reverse():
    """Stated as a test because it is a DECISION and users will meet it.

    The two spellings are different graphs and neither matches the other. The alternative -- an order
    that matches both -- would make `bond_order 1` ambiguous on every aromatic ring. A caller who
    wants one answer for both must normalise first.
    """
    arom, _ = benzene()
    kek, _ = kekule_benzene()
    assert _bond_query('bond_aromatic').count(kek) == 0
    assert _bond_query('bond_order', 2).count(arom) == 0
    # three double bonds, each matched from both ends: the count is over EMBEDDINGS, not bonds
    assert _bond_query('bond_order', 2).count(kek) == 6
    assert _bond_query('bond_order', 1).count(kek) == 6


def test_a_ring_bond_query_matches_either_spelling_and_an_acyclic_aromatic_bond_as_well():
    """`bond_ring` is about topology, so both spellings of a real ring answer it.

    AND SO DOES AN AROMATIC BOND IN NO RING AT ALL, which is a consequence of the layout rather than
    an accident, and is recorded here because it is the one answer somebody will file as a bug. The
    topology triple is a one-hot span of three bits over a full word, and order 4 has no order bit of
    its own, so `W0_BIT_RING_AROM` has to double as the aromatic bond's ORDER value -- it is set from
    the aromatic flag, before the ring test, precisely so that an aromatic bond outside a ring cannot
    fall through to `W0_BIT_NOT_RING` and become indistinguishable from a dative bond. `bond_ring`
    accepts that bit, so it accepts the bond.

    The alternative needs a spare bit, and word 0 is sixty-four of sixty-four (as is word 2). The
    trade is: an aromatic bond drawn outside any perceived ring answers "ring" -- a structure no
    correct file contains, since aromaticity is a ring property -- against every genuine aromatic
    bond in every real molecule answering "coordination bond" instead. This is the cheaper wrong
    answer, and it is only reachable through the arena's own API on a structure a chemist would
    reject.
    """
    for m, _ in (benzene(), kekule_benzene()):
        assert _bond_query('bond_ring').count(m) == 12
    chain = MoleculeContainer()
    x, y = chain.add_atom(6), chain.add_atom(6)
    chain.add_bond(x, y, 4)
    assert not chain.bond_in_ring(x, y), 'ring PERCEPTION is not fooled: there is no ring here'
    assert _bond_query('bond_ring').count(chain) == 2, 'but the screening bit says ring, by design'
    assert _bond_query('bond_aromatic').count(chain) == 2
    # a single bond in no ring is still refused, so `bond_ring` has not become a tautology
    single = MoleculeContainer()
    p, r = single.add_atom(6), single.add_atom(6)
    single.add_bond(p, r, 1)
    assert _bond_query('bond_ring').count(single) == 0


def test_a_negated_coordination_order_is_refused_rather_than_answered_wrongly():
    """The one demand the layout cannot state, and it says so.

    "not order 8" is {1, 2, 3, aromatic}, which over the shared bit is a disjunction, and a box is
    a conjunction of forbidden bits. Forbidding the order-8 bit alone would silently drop every
    aromatic bond. A construction error is the honest answer, and the caller can spell the same
    demand as an OR term.

    The error surfaces at SEAL, not at `bond_primitive`, because the container only journals tokens
    and the box compiler is what discovers the contradiction. So the refusal is provoked by a match,
    which is where a caller would meet it.
    """
    arom, _ = benzene()
    q = _bond_query('bond_order', 8, negated=True)
    with raises(ValueError, match='coordination'):
        q.count(arom)
    for order in (1, 2, 3):
        # the same demand over an order with a bit of its own is expressible, and answers
        assert _bond_query('bond_order', order, negated=True).count(arom) == 12


# ── stereo: the surface that was nearly gated ───────────────────────────────────────────


def test_an_aromatic_biaryl_keeps_its_atropisomer_axis():
    """A refusal here would have been the branch's worst mistake, so this is its guard.

    Axis detection reads degree, hydrogen count and ring membership -- never a bond order -- so an
    aromatic-written biaryl passes exactly as its Kekule twin does. The unit is reported on the
    pivot bond of a 2,2'-disubstituted biaryl in both spellings.
    """
    arom, ids = biphenyl(order=4, substituted=True)
    assert 3 in kinds(arom), 'SU_ATROPISOMER (kind 3) must be perceived on an aromatic biaryl'
    assert anchors(arom, 3) <= {ids[0], ids[6]}, 'and on the pivot bond, not somewhere else'


def test_the_aromatic_biaryl_and_its_kekule_twin_report_THE_SAME_UNIT_SET():
    """Equality, not containment -- and the reason is a cut that was already there.

    The premise this test was first written on was wrong and is worth recording, because it is the
    plausible one: a Kekule benzene has three order-2 bonds, so it looks as though it must offer
    three extra cis/trans candidates per ring that the aromatic spelling cannot. It does not.
    `SU_MIN_STEREO_RING = 8` refuses a cis/trans unit on any ring smaller than eight members, so
    benzene's double bonds are cut before an order is ever compared, and both spellings report the
    axis alone.

    That makes ruling F100 -- the collision between a ring double bond and an axis -- unreachable
    from either spelling at benzene size, which is why storing order 4 FIXES it rather than merely
    avoiding it: `_is_chain_bond` excludes an aromatic bond outright, so above the ring-size cut the
    aromatic spelling has no candidate to collide even where the Kekule one would.
    """
    arom, _ = biphenyl(order=4, substituted=True)
    kek, a = biphenyl(order=1, substituted=True)
    for base in (0, 6):                                  # a real alternating ring, not all-order-2
        for i in (0, 2, 4):
            kek.set_order(a[base + i], a[base + (i + 1) % 6], 2)
    assert kek.aromatic_bond_count == 0 and arom.aromatic_bond_count == 12
    # EQUALITY OF THE TWO SPELLINGS IS THE CLAIM, and it is asserted on its own line so that a
    # failure says which half of the sentence broke.
    assert kinds(arom) == kinds(kek)
    # 1 is SU_CIS_TRANS, and its absence is the thing the docstring argues for: neither spelling
    # offers a ring cis/trans candidate at benzene size.  The two methyl carbons carry three
    # hydrogens each, reach four directions and are therefore TETRAHEDRAL CANDIDATES (kind 0,
    # `stereogenic=False`) in both spellings alike; a methyl is not a stereocentre, and no claim here
    # depends on it not being a candidate.
    assert 1 not in kinds(arom)
    assert [u['kind'] for u in arom.stereo_units() if u['stereogenic']] == [3], \
        'the axis is the only STEREOGENIC unit, in both spellings'
    assert [u['kind'] for u in kek.stereo_units() if u['stereogenic']] == [3]


def test_a_benzylic_stereocentre_is_perceived_through_an_aromatic_bond():
    """An aromatic bond is a DIRECTION like any other: order 4 is not refused by the tetra walk.

    1-phenylethanol's carbinol carbon reaches four directions as phenyl, methyl, hydroxyl and one
    hydrogen, and one of those directions arrives over an aromatic bond. Refusing order 4 there -- the
    way order 3 and order 8 are refused -- would have lost the commonest stereocentre in medicinal
    chemistry.

    THE HYDROGEN IS STATED, not derived: the core does not compute implicit hydrogens (there is no
    `calc_implicit` here; that is the standardisation layer's job), so a carbon left at the arena's
    default reaches three directions and is not a unit in EITHER spelling. `set_hydrogens` is how a
    parser will say what the file said, and it is how this test says it.
    """
    m, a = ring([6] * 6, [4] * 6)
    c = m.add_atom(6)
    m.add_bond(a[0], c, 1)
    m.add_bond(c, m.add_atom(8), 1)
    m.add_bond(c, m.add_atom(6), 1)
    m.set_hydrogens(c, 1)
    assert c in anchors(m, 0), 'the carbinol carbon must be a tetrahedral unit'
    # and the aromatic bond is what carries one of its four directions
    assert m.order_of(c, a[0]) == 1 and m.hybridization_of(a[0]) == 4


def test_an_aromatic_sulfur_keeps_the_lone_pair_its_kekule_twin_has():
    """The one real defect the audit found, and the test that would have caught it.

    `_sulfur_lone_pair` spends one electron per bonding electron, so charging 4 for an aromatic bond
    put thiophene's sulfur 5 electrons over budget per ring bond and denied it a pair. The rule is
    one sigma electron per aromatic bond plus ONE delocalised pi for the atom, which reproduces the
    Kekule budget exactly: 2 + 1 == 1 + 1 + 1 here, and 3 + 1 == 1 + 1 + 2 at a ring fusion.

    Measured through a shape whose verdict the budget decides: a sulfonium-like S with three
    directions plus the pair is a tetrahedral unit, and it must be one in both spellings.
    """
    for orders in ([4, 4, 4, 4, 4], [1, 2, 1, 2, 1]):
        m = MoleculeContainer()
        # The sulfur is the ANCHOR, so its count is stated: perception refuses an anchor whose
        # hydrogen count is unknown, and `add_atom` stores H_UNKNOWN for an unsaid one.  A ring
        # sulfur with two ring bonds and an exocyclic methyl carries no hydrogen, and saying so is
        # what leaves the lone-pair budget as the only thing this test is measuring.
        a = [m.add_atom(16, implicit_h=0)] + [m.add_atom(6) for _ in range(4)]
        for i, o in enumerate(orders):
            m.add_bond(a[i], a[(i + 1) % 5], o)
        # substituents that make the ring positions distinguishable
        m.add_bond(a[1], m.add_atom(8), 1)
        m.add_bond(a[2], m.add_atom(7), 1)
        assert m.hybridization_of(a[0]) == (4 if orders[0] == 4 else 1)
        exo = m.add_atom(6)
        m.add_bond(a[0], exo, 1)          # S with three sigma directions; the pair completes four
        assert a[0] in anchors(m, 0), 'the sulfur must reach four directions in both spellings'


def test_every_stereo_read_answers_on_an_aromatic_molecule():
    """There is no refusal on this surface, and this is the test that says so.

    Five methods were gated on `is_kekule` for part of an afternoon: `stereo_units`, `chiral_atoms`,
    `chiral_bonds`, `is_chiral` and `validate_stereo`. The gate was wrong -- it would have refused
    most drug-like molecules to protect against a defect that was not there -- and its removal is
    what this asserts. `validate_stereo` is in the list on purpose: it MUTATES, so it exercises the
    clone-and-clear path on an aromatic buffer as well as the read paths.
    """
    m, a = biphenyl(order=4, substituted=True)
    assert m.stereo_units()
    m.chiral_atoms()
    m.chiral_bonds()
    m.is_chiral(a[0])
    m.stereogenic_units()
    assert m.validate_stereo() == []
    assert m.aromatic_bond_count == 12, 'and nothing normalised the bonds behind the caller'


# ── canonical form ──────────────────────────────────────────────────────────────────────


def test_the_canonical_order_of_an_aromatic_molecule_is_creation_order_independent():
    """The completeness question the invariant refinement has to answer for a NEW bond order.

    An aromatic ring is more symmetric than its Kekule twin -- benzene's six bonds are all alike,
    where the Kekule form alternates -- so if the bond invariant failed to separate anything, the
    refinement would stall and the canonical order would depend on the order the atoms were created
    in. Each molecule below is built in every permutation of its ring, and the number of DISTINCT
    canonical orders (as element/order sequences, which is what is invariant under renumbering) must
    be one.
    """
    for elements in ([6] * 6, [7] + [6] * 5, [7, 6, 6, 7, 6, 6]):
        seen = set()
        for perm in permutations(range(len(elements))):
            m = MoleculeContainer()
            ids = {}
            for position in perm:
                ids[position] = m.add_atom(elements[position])
            for i in range(len(elements)):
                m.add_bond(ids[i], ids[(i + 1) % len(elements)], 4)
            seen.add(canonical_string(m))
        assert len(seen) == 1, '%d distinct canonical forms over %d creation orders' % (
            len(seen), len(list(permutations(range(len(elements))))))


def test_an_aromatic_ring_and_its_kekule_twin_have_different_canonical_forms():
    """The consequence a caller has to know about, stated where it can be found.

    Same atoms, same connectivity, different stored bonds -- so different feature words, different
    Morgan colours and a different canonical order. Two spellings of benzene are not equal and do
    not hash alike. This is not a defect to be papered over with a normalising comparison; it is
    what storing the input faithfully means, and the fix for a caller who wants them equal is to
    kekulise or aromatise both first.
    """
    arom, _ = benzene()
    kek, _ = kekule_benzene()
    assert canonical_string(arom) != canonical_string(kek)
    assert {arom.hybridization_of(s) for s in arom.atom_numbers} == {4}
    assert {kek.hybridization_of(s) for s in kek.atom_numbers} == {2}


def test_naphthalene_written_aromatic_canonicalises_to_one_order_from_many_starts():
    """A fused system, because a single ring cannot show a wrong answer that only fusion produces.

    Naphthalene's two ring-fusion carbons carry three aromatic bonds each and are the only atoms
    that do; every other position is one of two orbits. Built from each of its ten atoms as the
    first-created one, the canonical order must be the same molecule every time.
    """
    edges = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0),
             (4, 6), (6, 7), (7, 8), (8, 9), (9, 5)]
    seen = set()
    for start in range(10):
        m = MoleculeContainer()
        ids = {}
        for k in range(10):
            ids[(start + k) % 10] = m.add_atom(6)
        for i, j in edges:
            m.add_bond(ids[i], ids[j], 4)
        seen.add(canonical_string(m))
    assert len(seen) == 1, '%d distinct canonical embeddings over 10 creation orders' % len(seen)
    # the fusion carbons are the only three-aromatic-bond atoms, so the shape really is what it says
    assert sorted(len(m.neighbors_of(s)) for s in m.atom_numbers) == [2] * 8 + [3] * 2


# ── what does NOT change ────────────────────────────────────────────────────────────────


def test_no_read_on_the_public_surface_changes_the_stored_bonds():
    """The invariant `aromatic_bond_count` exists to make testable, swept over the whole surface.

    Every public no-argument read and every per-atom read is called on an aromatic molecule, and the
    count must be untouched afterwards. Nothing sanitises, normalises, kekulises or aromatises
    behind the caller's back.

    THE MUTATORS ARE EXCLUDED BY NAME, and only the mutators: `dir()` does not distinguish a read
    from a write, so calling everything would call `delete_atom` (which of course changes the count)
    and prove nothing. The list is a prefix rule over the API's own naming convention plus the two
    scoped operations, so a NEW read cannot accidentally land in it -- a new read has to be named
    like a read to be skipped, and no read in this class is.

    `kekule()` IS THE ENTRY THIS MECHANISM WAS BUILT FOR, and it earned its place the intended way:
    the test was written with an empty exception list on a branch that had no `kekule()`, the
    kekulisation epic merged, and this assertion failed naming it. Nobody had to remember to come back
    and add it. `thiele()` will arrive the same way. The floor below is a named number rather than
    `> 0` so that a refactor which made most of the surface unreachable would fail here instead of
    passing vacuously.

    THE ALTERNATIVE SPELLINGS ARE SWEPT LIKE ANY OTHER READ, and not skipped by name: whether
    `atoms_count` mutates is exactly as much this test's business as whether `atom_count` does.
    """
    changes_bonds = {'edit', 'remap', 'copy_to',
                     'kekule',      # arrived by merge; caught by this test, not by memory
                     'thiele'}      # not yet implemented; listed so its arrival is not a surprise
    m, a = biphenyl(order=4, substituted=True)
    before = m.aromatic_bond_count
    answered = 0
    for name in sorted(dir(m)):
        if name.startswith('_') or name in changes_bonds or \
                name.startswith(('set_', 'add_', 'delete_')):
            continue
        attribute = getattr(type(m), name, None)
        try:
            if isinstance(attribute, property):
                getattr(m, name)
            elif callable(attribute):
                try:
                    attribute(m)                       # no-argument methods
                except TypeError:
                    attribute(m, a[0])                 # per-atom reads
            else:
                continue
        except Exception:
            continue                                   # a read that refuses is not a mutation
        answered += 1
        assert m.aromatic_bond_count == before, '%s changed the stored bonds' % name
    assert answered >= 40, 'only %d members answered, so this swept almost nothing' % answered
