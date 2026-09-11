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
# The primitive→box→DNF term compiler.
#
# ---------------------------------------------------------------------------
# Primitive constants and wbox_t
# ---------------------------------------------------------------------------

DEF Q_BOX_MAX_ANY = 8
DEF Q_ATOM_MAX_BOXES = 64

cdef enum:
    PRIM_ELEMENT      = 1    # value = atomic number 1..118
    PRIM_ANY          = 2    # [A]: any atom, no element constraint
    PRIM_METAL        = 3    # [M]
    PRIM_ISOTOPE      = 4    # value = ABSOLUTE mass number
    PRIM_CHARGE       = 5    # value = -4..8; saturates at the span edges
    PRIM_RADICAL      = 6    # value ignored
    PRIM_DEGREE       = 7
    PRIM_IMPLICIT_H   = 8
    PRIM_TOTAL_H      = 9
    PRIM_HETEROATOMS  = 10
    PRIM_HYBRIDIZATION = 11  # value 1..6 (see validator at :263)
    PRIM_RING_SIZE    = 12   # value 3..; multi-hot span
    PRIM_RING_COUNT   = 13
    PRIM_STEREO       = 14   # value 1 = '@', 2 = '@@'; contributes NO box bits (see prim_apply)
    PRIM_NO_ISOTOPE   = 15   # value ignored: demand the 'isotope absent' bit
    PRIM_ANY_CHARGE   = 16   # '*': value ignored; TOUCHES the charge span and forbids nothing
    PRIM_R_MARKER     = 17   # '#0': value ignored; a patch BUILD, and refused by prim_apply below
    PRIM_STEREO_KEEP  = 18   # '@=': value ignored; a patch BUILD, and refused by prim_apply below
    PRIM_STEREO_INVERT = 19  # '@~': value ignored; a patch BUILD, and refused by prim_apply below
    BPRIM_ORDER       = 20   # value 1, 2, 3, 8
    BPRIM_AROMATIC    = 21
    BPRIM_RING        = 22


cdef enum:
    # The stereo sign, as a bitmask, so that a box can carry "no demand" (0), one sign, or -- when
    # two stereo primitives are ANDed into the same box -- the contradiction 3.  Ruling F87: the
    # sign lives PER BOX, because a query atom's boxes are a disjunction and each disjunct states
    # its own configuration.  '[C@,N]' has an @-demanding box and a sign-free one; hoisting the
    # sign to the atom made the nitrogen box demand @ as well and dropped a real embedding.
    QSIGN_CW    = 1          # '@'  -- odd parity in the query's own reference frame
    QSIGN_CCW   = 2          # '@@' -- even parity in the query's own reference frame
    QSIGN_BOTH  = 3          # '[C;@;@@]': ANDed contradiction, refused at match time, not at seal
    QSIGN_FREE  = 4          # kernel-side only: some admitting box demands no sign at all


cdef struct wbox_t:                      # a box under construction, before the arena exists
    uint64_t neg[4]
    uint64_t touched[4]
    uint64_t any_mask[Q_BOX_MAX_ANY]
    uint32_t any_word[Q_BOX_MAX_ANY]
    uint32_t any_count
    uint32_t element_set_size            # how many elements this box still allows
    uint32_t element_single              # the element, when element_set_size == 1
    uint8_t sign                         # QSIGN_* bitmask; 0 = this disjunct says nothing (F87)


# Metal element masks for word 0 and word 1.
# Metal = any element that is not is_forming_single_bonds and not GroupXVIII.
# Non-metal elements (25 total): H(1) He(2) B(5) C(6) N(7) O(8) F(9) Ne(10)
#   Si(14) P(15) S(16) Cl(17) Ar(18) Ge(32) As(33) Se(34) Br(35) Kr(36)
#   Sb(51) Te(52) I(53) Xe(54) At(85) Rn(86) Og(118)
# 93 metals total: 34 light (e<=56) contribute 34 bits + 1 heavy-marker bit in word 0;
#   59 heavy (57<=e<=118, excluding 85/86/118) contribute 59 bits in word 1.
# Cross-checked against V2 isomorphism.py AnyMetal mask v1 = 0x0060707ffc1fff87 (word 0 matches).
DEF METAL_W0_ALLOWED  = 0x0060707ffc1fff87   # allowed metal bits in word 0 (element span)
DEF METAL_W1_ALLOWED  = 0x1fffffffcfffffff   # allowed metal bits in word 1 (heavy span)
DEF METAL_W0_FORBIDDEN = 0x019f8f8003e00078  # W0_ELEMENT_SPAN & ~METAL_W0_ALLOWED
DEF METAL_W1_FORBIDDEN = 0x2000000030000000  # W1_ELEMENT_SPAN & ~METAL_W1_ALLOWED


# ---------------------------------------------------------------------------
# Span table — single source for every one-hot feature span
# ---------------------------------------------------------------------------

cdef enum:
    SPAN_TOPOLOGY    = 0
    SPAN_ORDER       = 1
    SPAN_RADICAL     = 2
    SPAN_HETEROATOMS = 3
    SPAN_DEGREE      = 4
    SPAN_IMPLICIT_H  = 5
    SPAN_EXPLICIT_H  = 6
    SPAN_TOTAL_H     = 7
    SPAN_CHARGE      = 8
    SPAN_ISOTOPE     = 9
    SPAN_HYBRIDIZATION = 10
    SPAN_RING_COUNT  = 11
    SPAN_PARITY_KNOWN = 12

DEF SPAN_COUNT = 13

cdef extern from *:
    """
    /* Every one-hot feature span, as (word, mask) pairs.  A box that forbids all of a span it
       touched can never match.  The element span is absent on purpose: it straddles words 0 and
       1, so box_unsatisfiable tests it by inspecting the neg bits directly (see W0_LIGHT_ELEMENT_SPAN
       in _features.pxi).  ring_sizes is absent because it is multi-hot -- forbidding all of it
       means 'in no ring', which is satisfiable.
       Index names: SPAN_TOPOLOGY=0, SPAN_ORDER=1, SPAN_RADICAL=2, SPAN_HETEROATOMS=3,
       SPAN_DEGREE=4, SPAN_IMPLICIT_H=5, SPAN_EXPLICIT_H=6, SPAN_TOTAL_H=7, SPAN_CHARGE=8,
       SPAN_ISOTOPE=9, SPAN_HYBRIDIZATION=10, SPAN_RING_COUNT=11, SPAN_PARITY_KNOWN=12. */
    static const unsigned int SPAN_WORD[13] = {0, 0, 1, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3};
    static const unsigned long long SPAN_MASK[13] = {
        0x4600000000000000ULL,   /* word 0 bits 57/58/62: bond topology triple  */
        0xB800000000000000ULL,   /* word 0 bits 59/60/61/63: bond order         */
        0xC000000000000000ULL,   /* word 1 bits 62-63: radical / not radical    */
        0x00000000000001FFULL,   /* word 2 bits 0-8: heteroatom count           */
        0x000000000001FE00ULL,   /* word 2 bits 9-16: degree                    */
        0x00000000003E0000ULL,   /* word 2 bits 17-21: implicit hydrogens       */
        0x0000000007C00000ULL,   /* word 2 bits 22-26: explicit hydrogens       */
        0x00000001F8000000ULL,   /* word 2 bits 27-32: total hydrogens          */
        0x00003FFE00000000ULL,   /* word 2 bits 33-45: charge, biased by 4      */
        0xFFFFC00000000000ULL,   /* word 2 bits 46-63: isotope delta, 63 = none */
        0x000000000000003FULL,   /* word 3 bits 0-5: hybridization              */
        0x00FF800000000000ULL,   /* word 3 bits 47-55: ring count               */
        0x0000000000000180ULL};  /* word 3 bits 7-8: parity configured / not    */
    /* Bits that boxes_merge's span loop actually compares, per word.  Words 0-2 are fully covered
       (element + topology + order; element + radical; the seven word-2 count spans).  Word 3 is
       not: bits 6-46 (stereo, ring_sizes) and 56-63 (aromatic ring count) are outside every
       SPAN_MASK entry, so two boxes differing only there would otherwise look identical to the
       span loop.  Derived as: OR of SPAN_MASK entries for each word, plus the element span.
       Word 0: W0_ELEMENT_SPAN | SPAN_TOPOLOGY | SPAN_ORDER = 0xFFFFFFFFFFFFFFFF
       Word 1: W1_ELEMENT_SPAN | SPAN_RADICAL = 0xFFFFFFFFFFFFFFFF
       Word 2: all seven word-2 spans = 0xFFFFFFFFFFFFFFFF
       Word 3: SPAN_HYBRIDIZATION | SPAN_RING_COUNT = 0x00FF80000000003F
       Word 3 bit 6, the stereo VALUE bit, is outside every SPAN_MASK entry and outside SPAN_COVERED,
       and for ruling F78's reason: prim_apply puts no bit in a box for the sign, because a stored
       parity is in the molecule's frame and a query's is in its own, so there is nothing for the
       span loop to compare.  Bits 7-8 are the parity-configured span (configured / not configured),
       and they are TWO bits so that a box can state the demand at all: over a one-bit span
       wbox_forbid_one_hot has no rest of the span to forbid, and box_unsatisfiable would read
       "forbids the whole span" off any box that merely excluded the bit.

       They are in SPAN_MASK and deliberately NOT here.  Membership of SPAN_COVERED only LOOSENS
       boxes_merge: same_uncovered compares the bits OUTSIDE this mask for equality, so putting a span
       in stops that check from seeing it.  What protects the parity demand across a merge is the
       `bi.sign != bj.sign` guard, which predates the span: box.sign is written only by the
       PRIM_STEREO branch that sets bit 8 in the same breath, so sign != 0 if and only if the demand
       is present, and two boxes that could differ in it differ in sign first and are refused a merge
       before same_uncovered is reached.  Leaving the span uncovered therefore costs nothing and keeps
       the fallback strict rather than the guard wide.  One consequence to state rather than leave for
       a reader to rediscover: a span in SPAN_MASK but not here can never be the SINGLE differing span
       boxes_merge widens, because differing there fails same_uncovered first.  Bits 7-8's SPAN_MASK
       entry is thus inert for merging while still load-bearing for wbox_forbid_one_hot and
       box_unsatisfiable, which is the strict direction and not a hazard. */
    static const unsigned long long SPAN_COVERED[4] = {
        0xFFFFFFFFFFFFFFFFULL, 0xFFFFFFFFFFFFFFFFULL,
        0xFFFFFFFFFFFFFFFFULL, 0x00FF80000000003FULL};
    """
    const uint32_t SPAN_WORD[SPAN_COUNT]
    const uint64_t SPAN_MASK[SPAN_COUNT]
    const uint64_t SPAN_COVERED[4]


cdef inline int wbox_forbid_one_hot(wbox_t *box, uint32_t word, uint64_t span, uint64_t bit,
                                    bint negated) noexcept nogil:
    """value in {bit} over a one-hot span, or its complement when negated."""
    box.touched[word] |= span
    if negated:
        box.neg[word] |= bit
    else:
        box.neg[word] |= span & ~bit
    return 0


cdef int wbox_require_any(wbox_t *box, uint32_t word, uint64_t span, uint64_t bit) except -1:
    """At least one of `bit` present -- the only shape a multi-hot span admits."""
    # Entries are never merged, not even two on the same word: [c;r5;r6] is two independent
    # demands, and OR-ing them into one mask would turn the AND into an OR.
    box.touched[word] |= span
    if box.any_count >= Q_BOX_MAX_ANY:
        raise ValueError('too many positive multi-hot constraints on one query atom; the cap '
                         'is %d' % Q_BOX_MAX_ANY)
    box.any_mask[box.any_count] = bit
    box.any_word[box.any_count] = word
    box.any_count += 1
    return 0


cdef int prim_apply(wbox_t *box, uint32_t kind, int32_t value, bint negated) except -1:
    """Apply one SMARTS primitive to a wbox_t under construction.

    ORs the primitive's forbidden bits into box.neg, its span into box.touched, appends
    an any entry for a positive test on a multi-hot span, and maintains element_set_size /
    element_single for the screen and for isotope resolution.
    """
    cdef uint32_t bit, e
    cdef int32_t delta
    cdef uint64_t ring_span

    ring_span = <uint64_t> 0x7FFFFFC00000    # word 3 bits 22-46

    if kind == PRIM_ELEMENT:
        e = <uint32_t> value
        if e < 1 or e > 118:
            raise ValueError('element number %d is out of range 1..118' % value)
        if not negated:
            # Update element bookkeeping before touching (touched == 0 means no prior element)
            if not box.touched[0] & W0_ELEMENT_SPAN:
                # First non-negated element primitive
                box.element_set_size = 1
                box.element_single = e
            elif box.element_set_size == 1:
                if box.element_single != e:
                    box.element_set_size = 0   # unsatisfiable: two different elements ANDed
            # If already 0 (unsatisfiable after prior conflict), leave it
        box.touched[0] |= W0_ELEMENT_SPAN
        box.touched[1] |= W1_ELEMENT_SPAN
        if negated:
            # Forbid exactly this element's bit(s), touch nothing extra
            if e <= 56:
                box.neg[0] |= <uint64_t> 1 << (57 - e)
            else:
                # Word 1 alone carries the identity, so forbidding its bit excludes exactly this
                # element.  Word 0 bit 0 is the heavy *marker*, shared by every element above 56 —
                # forbidding it here would make [!U] reject thorium and every other heavy element.
                # Light elements leave word 1's element span empty, so they pass untouched.
                box.neg[1] |= <uint64_t> 1 << (e - 57)
        else:
            # Forbid every other element: span minus this element's bit
            if e <= 56:
                box.neg[0] |= W0_ELEMENT_SPAN & ~(<uint64_t> 1 << (57 - e))
                box.neg[1] |= W1_ELEMENT_SPAN            # all heavy-element bits forbidden
            else:
                box.neg[0] |= W0_ELEMENT_SPAN & ~(<uint64_t> 1)   # forbid all light elements
                box.neg[1] |= W1_ELEMENT_SPAN & ~(<uint64_t> 1 << (e - 57))

    elif kind == PRIM_ANY:
        pass   # touches nothing; any atom is allowed

    elif kind == PRIM_METAL:
        if negated:
            # Forbid metal elements (keep only non-metals)
            box.neg[0] |= <uint64_t> METAL_W0_ALLOWED & W0_ELEMENT_SPAN
            box.neg[1] |= <uint64_t> METAL_W1_ALLOWED & W1_ELEMENT_SPAN
        else:
            # Forbid non-metal elements
            box.neg[0] |= <uint64_t> METAL_W0_FORBIDDEN
            box.neg[1] |= <uint64_t> METAL_W1_FORBIDDEN
        box.touched[0] |= W0_ELEMENT_SPAN
        box.touched[1] |= W1_ELEMENT_SPAN

    elif kind == PRIM_ISOTOPE:
        if box.element_set_size != 1:
            raise ValueError('an isotope primitive needs a settled element in the same box')
        delta = value - <int32_t> MDL_ISOTOPE[box.element_single]
        bit = 46 + _bit_of(delta, -8, 8)
        wbox_forbid_one_hot(box, 2, SPAN_MASK[SPAN_ISOTOPE], <uint64_t> 1 << bit, negated)

    elif kind == PRIM_NO_ISOTOPE:
        wbox_forbid_one_hot(box, 2, SPAN_MASK[SPAN_ISOTOPE],
                            <uint64_t> 1 << 63, negated)

    elif kind == PRIM_CHARGE:
        bit = 33 + _bit_of(value, -4, 8)
        wbox_forbid_one_hot(box, 2, SPAN_MASK[SPAN_CHARGE], <uint64_t> 1 << bit, negated)

    elif kind == PRIM_ANY_CHARGE:
        # "any charge" is not a demand, it is the WITHDRAWAL of one: `box_fill_defaults` neutralises
        # a charge span nobody touched, so touching it and forbidding nothing is the whole primitive.
        # One box and zero bits, which is why this is a primitive rather than a thirteen-way OR over
        # the span -- the OR is what the chemistry layer's `_free_charge_and_radical` builds today,
        # and it costs thirteen boxes per wildcard atom to say the same thing.
        #
        # It composes with an explicit charge in the same box the way every other box bit composes:
        # a box is a conjunction of forbidden bits, `*` forbids none, so `[C;*;+2]` is still +2.
        if negated:
            raise ValueError('`!*` is not a demand: `*` withdraws the charge default, and there is '
                             'nothing to negate.  Spell the charge you mean')
        box.touched[2] |= SPAN_MASK[SPAN_CHARGE]

    elif kind == PRIM_R_MARKER:
        # `#0` is the R marker, and an R matches nothing -- `MoleculeContainer.as_query` refuses a
        # molecule carrying one for the same reason.  It is a SMIRKS product-side build spelling, and
        # a product side never seals, so this branch is reached only by a query that holds one.  The
        # refusal is here rather than in the lexer because one lexer serves both sides of an arrow.
        raise ValueError('`#0` is the R marker, which matches nothing, so it cannot be part of a '
                         'query; it states an attachment point to BUILD on a SMIRKS product side')

    elif kind == PRIM_STEREO_KEEP:
        # `@=` says the configuration here is the one the reactant had, which is a statement about a
        # CHANGE.  A query changes nothing, so it has nothing to keep -- and "is configured, either
        # way" is not what the token means and is not spellable by design (see `PRIM_STEREO` below).
        # Reached only by a query, for the reason the R marker's branch gives.
        raise ValueError('`@=` carries a configuration through a patch, so it cannot be part of a '
                         'query; it states on a SMIRKS product side that the reactant\'s '
                         'configuration survives the reaction')

    elif kind == PRIM_STEREO_INVERT:
        # Same case as `@=` one branch up: a statement about a CHANGE, and a query changes nothing.
        raise ValueError('`@~` inverts a configuration through a patch, so it cannot be part of a '
                         'query; it states on a SMIRKS product side that the reactant\'s '
                         'configuration comes out as the other one')

    elif kind == PRIM_RADICAL:
        # word 1 bits 62-63: bit 62 = not-radical, bit 63 = radical
        if negated:
            box.touched[1] |= SPAN_MASK[SPAN_RADICAL]
            box.neg[1] |= <uint64_t> 1 << 63
        else:
            box.touched[1] |= SPAN_MASK[SPAN_RADICAL]
            box.neg[1] |= <uint64_t> 1 << 62

    elif kind == PRIM_DEGREE:
        bit = 9 + _bit_of(value, 0, 7)
        wbox_forbid_one_hot(box, 2, SPAN_MASK[SPAN_DEGREE], <uint64_t> 1 << bit, negated)

    elif kind == PRIM_IMPLICIT_H:
        bit = 17 + _bit_of(value, 0, 4)
        wbox_forbid_one_hot(box, 2, SPAN_MASK[SPAN_IMPLICIT_H], <uint64_t> 1 << bit, negated)

    elif kind == PRIM_TOTAL_H:
        bit = 27 + _bit_of(value, 0, 5)
        wbox_forbid_one_hot(box, 2, SPAN_MASK[SPAN_TOTAL_H], <uint64_t> 1 << bit, negated)

    elif kind == PRIM_HETEROATOMS:
        bit = _bit_of(value, 0, 8)
        wbox_forbid_one_hot(box, 2, SPAN_MASK[SPAN_HETEROATOMS], <uint64_t> 1 << bit, negated)

    elif kind == PRIM_HYBRIDIZATION:
        if value < 1 or value > 6:
            raise ValueError('hybridization value %d is out of range 1..6' % value)
        bit = <uint32_t> (value - 1)
        wbox_forbid_one_hot(box, 3, SPAN_MASK[SPAN_HYBRIDIZATION], <uint64_t> 1 << bit, negated)

    elif kind == PRIM_RING_SIZE:
        if value < 3:
            raise ValueError('ring size %d is below the minimum of 3' % value)
        if value <= 24:
            bit = 22 + <uint32_t> value
        elif value <= 32:
            bit = 22       # bucket: sizes 25-32
        elif value <= 48:
            bit = 23       # bucket: sizes 33-48
        else:
            bit = 24       # bucket: sizes > 48
        if negated:
            box.touched[3] |= ring_span
            box.neg[3] |= <uint64_t> 1 << bit
        else:
            wbox_require_any(box, 3, ring_span, <uint64_t> 1 << bit)

    elif kind == PRIM_RING_COUNT:
        bit = 47 + _bit_of(value, 0, 8)
        wbox_forbid_one_hot(box, 3, SPAN_MASK[SPAN_RING_COUNT], <uint64_t> 1 << bit, negated)

    elif kind == PRIM_STEREO:
        # Validated here and nowhere else -- this is the only gate every stereo token passes
        # through, whatever built it.
        if value != 1 and value != 2:
            raise ValueError("stereo value %d is not 1 ('@') or 2 ('@@')" % value)
        if negated:
            # '[C;!@]' has no meaning to give: the negation of "even in the query's frame" is
            # "odd OR not configured at all", and the second half is what an absent primitive
            # already says.  Refusing is a construction error, not a silent widening.
            raise ValueError('a stereo primitive cannot be negated')
        # NO BOX BIT FOR THE SIGN, deliberately.  Feature word IV bit 6 carries the parity as stored
        # in the MOLECULE's ruling-F26 frame, while this value is a statement in the QUERY's frame,
        # and the two frames differ by whatever permutation the embedding turns out to be.  A box
        # demanding bit 6 for '@@' would therefore reject a target whose stored parity is even but
        # whose parity re-expressed in the query's frame is odd, which is a true match:
        # test_the_primitive_is_read_in_the_querys_own_frame is exactly that molecule.
        # So the sign cannot be screened by feature bits at all; it rides the box as a side channel
        # (ruling F87) and the kernel reads it off whichever box admitted the candidate.  Bit 6
        # stays outside SPAN_COVERED[3] as a consequence: nothing may put it in a box, so there is
        # nothing for the span loop to compare.
        box.sign |= <uint8_t> (QSIGN_CW if value == 1 else QSIGN_CCW)
        # What CAN be screened is that a configuration is there at all: either sign demands a
        # configured centre (ruling F54 -- an unset parity is not a wedge that happens to point the
        # other way), and "configured" is frame-free.  This is the whole of the box's contribution
        # and it is what puts bit 7 into the query signature, so query_may_match refuses a target
        # with no configured centre at all before matcher_init builds the unit table for it.
        wbox_forbid_one_hot(box, 3, SPAN_MASK[SPAN_PARITY_KNOWN], <uint64_t> 1 << 7, False)

    elif kind == BPRIM_ORDER:
        # Order 4 (aromatic) is NOT one of the values here, and the reason is the word-0 budget:
        # all sixty-four bits are spent (see `w0_bond_bits`), so a stored order 4 shares
        # W0_BIT_ORDER8 with a dative bond and is told apart from it by W0_BIT_RING_AROM, which
        # fires exactly when the order is 4.  `bond_aromatic` is how a query demands one.
        if value == 1:
            bit = W0_BIT_ORDER1
        elif value == 2:
            bit = W0_BIT_ORDER2
        elif value == 3:
            bit = W0_BIT_ORDER3
        elif value == 8:
            bit = W0_BIT_ORDER8
            if negated:
                # "NOT a coordination bond" is the one demand this layout cannot state, and it is
                # refused rather than approximated.  The exact set is {1, 2, 3, aromatic}, which
                # over these bits is "W0_BIT_ORDER8 clear OR W0_BIT_RING_AROM set" -- a disjunction,
                # and a box is a conjunction of forbidden bits.  Forbidding W0_BIT_ORDER8 alone
                # would also reject every aromatic bond, silently dropping true embeddings, so the
                # construction error is preferable: the caller can spell the same thing as an OR
                # term over four positive primitives ('-,=,#,:'), which is what SMARTS writes
                # anyway.  Unreachable from the SMARTS surface, which has no bond negation beyond
                # '!@'.
                raise ValueError('a negated coordination-bond order cannot be expressed; write the '
                                 'orders you do want as an OR term instead')
            # Positive: forbid the rest of the order span AND the aromatic bit, or an aromatic bond
            # would satisfy a demand for a coordination bond.
            box.neg[0] |= <uint64_t> 1 << W0_BIT_RING_AROM
            box.touched[0] |= W0_TOPOLOGY_SPAN
        else:
            raise ValueError('bond order value %d is not valid (use 1, 2, 3, or 8)' % value)
        wbox_forbid_one_hot(box, 0, <uint64_t> W0_ORDER_SPAN, <uint64_t> 1 << bit, negated)

    elif kind == BPRIM_AROMATIC:
        # Aromatic bond: topology bit W0_BIT_RING_AROM (57), which `w0_bond_bits` sets iff the
        # stored order is 4.  The other topology bits (W0_BIT_NOT_RING=58 and
        # W0_BIT_RING_PLAIN=62) are forbidden, and forbidding them is also what excludes orders
        # 1/2/3 -- every such bond takes one of those two.
        wbox_forbid_one_hot(box, 0, W0_TOPOLOGY_SPAN,
                            <uint64_t> 1 << W0_BIT_RING_AROM, negated)

    elif kind == BPRIM_RING:
        # Ring bond: either plain ring (bit 62) or aromatic ring (bit 57). Not-ring (bit 58)
        # is forbidden. Negated means acyclic: forbid both ring bits.
        if negated:
            # Not in a ring: forbid bits 57 and 62
            box.touched[0] |= W0_TOPOLOGY_SPAN
            box.neg[0] |= (<uint64_t> 1 << W0_BIT_RING_AROM) | (<uint64_t> 1 << W0_BIT_RING_PLAIN)
        else:
            # In a ring: forbid bit 58 (not-ring)
            box.touched[0] |= W0_TOPOLOGY_SPAN
            box.neg[0] |= <uint64_t> 1 << W0_BIT_NOT_RING

    else:
        raise ValueError('unknown primitive kind %d' % kind)
    return 0


with cython.warn.undeclared(False):
    # bare so Python can import it, guarded so warn.undeclared stays quiet
    PRIM_NAMES = {
        'element': PRIM_ELEMENT, 'any': PRIM_ANY, 'metal': PRIM_METAL,
        'isotope': PRIM_ISOTOPE, 'charge': PRIM_CHARGE, 'radical': PRIM_RADICAL,
        'degree': PRIM_DEGREE, 'implicit_h': PRIM_IMPLICIT_H, 'total_h': PRIM_TOTAL_H,
        'heteroatoms': PRIM_HETEROATOMS, 'hybridization': PRIM_HYBRIDIZATION,
        'ring_size': PRIM_RING_SIZE, 'ring_count': PRIM_RING_COUNT,
        'stereo': PRIM_STEREO, 'no_isotope': PRIM_NO_ISOTOPE, 'any_charge': PRIM_ANY_CHARGE,
        'r_marker': PRIM_R_MARKER, 'stereo_keep': PRIM_STEREO_KEEP,
        'stereo_invert': PRIM_STEREO_INVERT,
        'bond_order': BPRIM_ORDER, 'bond_aromatic': BPRIM_AROMATIC, 'bond_ring': BPRIM_RING,
    }
    # The only name -> opcode table.  Every probe and QueryContainer goes through it, so an
    # unknown operator name raises KeyError exactly as an unknown primitive already does.
    OPC_NAMES = {'or': OPC_OR, 'and_low': OPC_AND_LOW, 'and_high': OPC_AND_HIGH}


def _prim_probe(str name, int32_t value, bint negated):
    return _prim_probe_seeded(name, value, negated, 0)


def _prim_probe_seeded(str name, int32_t value, bint negated, int element):
    cdef wbox_t box
    cdef uint32_t i
    cdef list any_out = []
    memset(&box, 0, sizeof(wbox_t))
    if element:
        prim_apply(&box, PRIM_ELEMENT, element, False)
        memset(&box.neg, 0, sizeof(box.neg))       # keep only the settled element, not its mask
        memset(&box.touched, 0, sizeof(box.touched))
    prim_apply(&box, PRIM_NAMES[name], value, negated)
    for i in range(box.any_count):
        any_out.append((box.any_word[i], box.any_mask[i]))
    return {'neg': (box.neg[0], box.neg[1], box.neg[2], box.neg[3]),
            'touched': (box.touched[0], box.touched[1], box.touched[2], box.touched[3]),
            'any': tuple(any_out), 'sign': box.sign}


def _query_alloc_probe(uint32_t atom_count, uint32_t bond_count, uint32_t box_count,
                       uint32_t bond_box_count, uint32_t any_count, uint32_t closure_count,
                       uint32_t component_count, uint32_t automorphism_count=0,
                       uint32_t stereo_count=0):
    cdef Query q = query_alloc(atom_count, bond_count, box_count, bond_box_count, any_count,
                               closure_count, component_count, automorphism_count, 0, stereo_count)
    cdef qclosure_t *closures = query_closures(q)
    cdef uint32_t *demand = query_element_demand(q)
    cdef int seg
    cdef uint32_t i
    cdef list segments = []
    cdef list demand_out = []
    for seg in range(QSEG_COUNT):
        segments.append((q.header.segments[seg].offset, q.header.segments[seg].length))
    for i in range(120):
        demand_out.append(demand[i])
    return {'magic': q.header.magic, 'version': q.header.version,
            'atom_count': q.header.atom_count, 'total_len': q.header.total_len,
            'automorphism_count': q.header.automorphism_count,
            'stereo_count': q.header.stereo_count,
            'segments': segments, 'element_demand': demand_out,
            'closure_reads': [closures[0].to_index, closures[0].bond_index]}


# ---------------------------------------------------------------------------
# Logic compiler
# ---------------------------------------------------------------------------

cdef enum:
    OPC_PRIM     = 0   # a primitive: (kind, value, negated)
    OPC_AND_LOW  = 1   # ';'
    OPC_OR       = 2   # ','
    OPC_AND_HIGH = 3   # '&' or an implicit juxtaposition


cdef struct qtoken_t:
    uint32_t opcode
    uint32_t kind
    int32_t value
    bint negated


cdef struct wterm_t:                 # one atom's disjunction under construction
    wbox_t boxes[Q_ATOM_MAX_BOXES]
    uint32_t count


cdef int box_fill_defaults(wbox_t *box) except -1:
    """A default constrains a span that no primitive in THIS box touched.

    Charge -> neutral, radical -> not a radical.  Everything else stays free: a query that says
    nothing about degree matches any degree.  Called after the cross-product, per final box, so
    that `[C;+]` -- whose single box touches both spans -- keeps its charge.

    NOTE on bond boxes (Watch item): compile_term applies this function unconditionally, so bond
    boxes also receive neutral-charge and not-radical demands in words 1 and 2.  This is harmless
    because `_fold_bond_into_atom` folds bond boxes by ORing neg[0] alone, and the kernel tests bonds
    against an edge-feature word (word 0 only) -- words 1-3 of a bond box are never consulted.  Do not
    "fix" these defaults without verifying every consumer of bond boxes.
    """
    if not (box.touched[2] & SPAN_MASK[SPAN_CHARGE]):
        box.neg[2] |= SPAN_MASK[SPAN_CHARGE] & ~(<uint64_t> 1 << 37)
        box.touched[2] |= SPAN_MASK[SPAN_CHARGE]
    if not (box.touched[1] & SPAN_MASK[SPAN_RADICAL]):
        box.neg[1] |= <uint64_t> 1 << 63
        box.touched[1] |= SPAN_MASK[SPAN_RADICAL]
    return 0


cdef bint box_unsatisfiable(wbox_t *box) noexcept nogil:
    """A box that forbids every bit of a span it touched can never match.

    sign == QSIGN_BOTH ('[C;@;@@]') is NOT pruned here, on purpose.  Pruning it would delete the
    only box of '[C;@;@@]' and the @@-half of '[C;@,N;@@]', and compile_term_unmerged turns an
    empty term into ValueError -- which ruling F87 forbids for a merely unsatisfiable stereo query,
    for the same reason F77 case 2 gives.  The kernel refuses such a box at match time instead.
    """
    cdef uint32_t i
    cdef uint64_t light
    # Element span straddles words 0 and 1; test by inspecting neg bits directly.
    # light = bits 1-56, the light-element identity bits (bit 0 is the heavy-element flag).
    # No element matches if: all light bits forbidden AND (heavy flag forbidden OR all heavy
    # identity bits forbidden).
    light = <uint64_t> W0_LIGHT_ELEMENT_SPAN
    if (box.neg[0] & light) == light and (
            box.neg[0] & <uint64_t> 1 or
            (box.neg[1] & <uint64_t> W1_ELEMENT_SPAN) == <uint64_t> W1_ELEMENT_SPAN):
        return True
    for i in range(SPAN_COUNT):
        if box.touched[SPAN_WORD[i]] & SPAN_MASK[i] and \
                box.neg[SPAN_WORD[i]] & SPAN_MASK[i] == SPAN_MASK[i]:
            return True
    for i in range(box.any_count):
        # a positive multi-hot demand for a bit the same box forbids
        if box.neg[box.any_word[i]] & box.any_mask[i]:
            return True
    return False


cdef int boxes_merge(wterm_t *term) except -1:
    """Collapse boxes that differ in exactly one SPAN_MASK span into one box.

    Two boxes are eligible only when: (a) their any-lists are identical, (b) their element words
    are identical (Ruling 2), (c) their stereo signs are identical (ruling F87 -- the sign is not a
    feature bit, so no span loop would notice it, and merging '[C;@,D3]' into one box is exactly how
    the @ demand leaked onto the D3 disjunct), (d) their bits outside SPAN_COVERED are identical
    (otherwise a negated ring_size bit lives in the uncovered region and would be silently lost),
    and (e) at most one SPAN_MASK span differs.  diff_cnt == 0 is a sound dedup of identical boxes.
    Repeat until stable (capped at Q_ATOM_MAX_BOXES passes).
    """
    cdef wbox_t *bi
    cdef wbox_t *bj
    cdef uint32_t pass_cnt, i, j, k, m, diff_cnt, diff_idx, diff_word, new_count
    cdef uint64_t diff_mask
    cdef bint changed, same_any, same_uncovered
    cdef uint8_t merged[Q_ATOM_MAX_BOXES]

    memset(merged, 0, Q_ATOM_MAX_BOXES)
    pass_cnt = 0
    changed = True
    while changed and pass_cnt < Q_ATOM_MAX_BOXES:
        changed = False
        pass_cnt += 1
        for i in range(term.count):
            if merged[i]:
                continue
            bi = &term.boxes[i]
            for j in range(i + 1, term.count):
                if merged[j]:
                    continue
                bj = &term.boxes[j]
                # Same any list (count and entries in order)
                if bi.any_count != bj.any_count:
                    continue
                same_any = True
                for k in range(bi.any_count):
                    if bi.any_word[k] != bj.any_word[k] or bi.any_mask[k] != bj.any_mask[k]:
                        same_any = False
                        break
                if not same_any:
                    continue
                # Element words: any difference blocks the merge (Ruling 2)
                if (bi.neg[0] & <uint64_t> W0_ELEMENT_SPAN !=
                        bj.neg[0] & <uint64_t> W0_ELEMENT_SPAN):
                    continue
                if (bi.neg[1] & <uint64_t> W1_ELEMENT_SPAN !=
                        bj.neg[1] & <uint64_t> W1_ELEMENT_SPAN):
                    continue
                # Stereo sign: any difference blocks the merge (ruling F87)
                if bi.sign != bj.sign:
                    continue
                # Uncovered bits (ring_sizes, stereo, aromatic ring count in word 3) must be
                # identical; otherwise a negated ring_size demand would be silently dropped.
                same_uncovered = True
                for k in range(4):
                    if (bi.neg[k] & ~SPAN_COVERED[k]) != (bj.neg[k] & ~SPAN_COVERED[k]):
                        same_uncovered = False
                        break
                if not same_uncovered:
                    continue
                # Count SPAN_MASK spans where neg differs
                diff_cnt = 0
                diff_idx = 0
                for k in range(SPAN_COUNT):
                    m = SPAN_WORD[k]
                    if (bi.neg[m] & SPAN_MASK[k]) != (bj.neg[m] & SPAN_MASK[k]):
                        diff_cnt += 1
                        diff_idx = k
                if diff_cnt > 1:
                    continue
                if diff_cnt == 1:
                    # i absorbs j -- OR the allowed halves of the differing span.  Widening
                    # never makes a box unsatisfiable.  element_set_size is inherited from i;
                    # sound because the element words are identical.
                    diff_word = SPAN_WORD[diff_idx]
                    diff_mask = SPAN_MASK[diff_idx]
                    bi.neg[diff_word] = ((bi.neg[diff_word] & ~diff_mask) |
                                         (bi.neg[diff_word] & bj.neg[diff_word] & diff_mask))
                # diff_cnt == 0 means every covered and uncovered bit matches, so absorbing j
                # is a plain dedup and the touched union below is the whole of it.
                for k in range(4):
                    bi.touched[k] |= bj.touched[k]
                merged[j] = 1
                changed = True
    # Compact: remove absorbed boxes
    new_count = 0
    for i in range(term.count):
        if not merged[i]:
            if new_count != i:
                term.boxes[new_count] = term.boxes[i]
            new_count += 1
    term.count = new_count
    return 0


cdef int _cross_term(wterm_t *result, wterm_t *left, wterm_t *right) except -1:
    """Cross-product: result = left × right.  When left.count == 0, result is a copy of right.

    result must not alias left or right.  No additional wterm_t is allocated by the caller;
    three are already live in compile_term_unmerged's frame (out, alt, tmp).
    """
    cdef wbox_t *res
    cdef wbox_t *lb
    cdef wbox_t *rb
    cdef uint32_t i, j, k, ai, n_boxes, lac, rac
    if left.count == 0:
        result[0] = right[0]
        return 0
    n_boxes = left.count * right.count
    if n_boxes > Q_ATOM_MAX_BOXES:
        raise ValueError('box count would exceed the cap of %d' % Q_ATOM_MAX_BOXES)
    memset(result, 0, sizeof(wterm_t))
    for i in range(left.count):
        lb = &left.boxes[i]
        for j in range(right.count):
            rb = &right.boxes[j]
            res = &result.boxes[i * right.count + j]
            for k in range(4):
                res.neg[k] = lb.neg[k] | rb.neg[k]
                res.touched[k] = lb.touched[k] | rb.touched[k]
            # ANDing two disjuncts ANDs their stereo demands; CW|CCW = QSIGN_BOTH, which no
            # target satisfies and which the kernel -- not this function -- refuses.
            res.sign = lb.sign | rb.sign
            lac = lb.any_count
            rac = rb.any_count
            if lac + rac > Q_BOX_MAX_ANY:
                raise ValueError(
                    'too many positive multi-hot constraints on one query atom; '
                    'the cap is %d' % Q_BOX_MAX_ANY)
            res.any_count = lac + rac
            for ai in range(lac):
                res.any_mask[ai] = lb.any_mask[ai]
                res.any_word[ai] = lb.any_word[ai]
            for ai in range(rac):
                res.any_mask[lac + ai] = rb.any_mask[ai]
                res.any_word[lac + ai] = rb.any_word[ai]
            # element bookkeeping: intersect
            # (element_set_size is what PRIM_ISOTOPE's delta resolution needs)
            if lb.touched[0] & W0_ELEMENT_SPAN and rb.touched[0] & W0_ELEMENT_SPAN:
                if (lb.element_set_size == 1 and rb.element_set_size == 1 and
                        lb.element_single == rb.element_single):
                    res.element_set_size = 1
                    res.element_single = lb.element_single
                else:
                    res.element_set_size = 0
                    res.element_single = 0
            elif lb.touched[0] & W0_ELEMENT_SPAN:
                res.element_set_size = lb.element_set_size
                res.element_single = lb.element_single
            elif rb.touched[0] & W0_ELEMENT_SPAN:
                res.element_set_size = rb.element_set_size
                res.element_single = rb.element_single
            else:
                res.element_set_size = 0
                res.element_single = 0
    result.count = left.count * right.count
    return 0


cdef int compile_term_unmerged(wterm_t *out, qtoken_t *tokens, uint32_t count) except -1:
    """Tokens -> cross-product -> box_fill_defaults -> prune unsatisfiable.

    Raises ValueError for malformed streams (trailing op, two ops in a row, empty stream) and
    for wholly unsatisfiable terms.  Does NOT call boxes_merge.
    """
    cdef wbox_t cur
    cdef wterm_t alt, tmp
    cdef qtoken_t *tok
    cdef uint32_t i, j
    cdef bint expect_prim

    if count == 0:
        raise ValueError('malformed query term')

    memset(&cur, 0, sizeof(wbox_t))
    memset(&alt, 0, sizeof(wterm_t))
    memset(out, 0, sizeof(wterm_t))
    expect_prim = True

    for i in range(count):
        tok = &tokens[i]
        if tok.opcode == OPC_PRIM:
            if not expect_prim:
                raise ValueError('malformed query term')
            prim_apply(&cur, tok.kind, tok.value, tok.negated)
            expect_prim = False

        elif tok.opcode == OPC_AND_HIGH:
            if expect_prim:
                raise ValueError('malformed query term')
            expect_prim = True

        elif tok.opcode == OPC_OR:
            if expect_prim:
                raise ValueError('malformed query term')
            # push cur into alt
            if alt.count >= Q_ATOM_MAX_BOXES:
                raise ValueError('too many alternatives; the cap is %d' % Q_ATOM_MAX_BOXES)
            alt.boxes[alt.count] = cur
            alt.count += 1
            memset(&cur, 0, sizeof(wbox_t))
            expect_prim = True

        elif tok.opcode == OPC_AND_LOW:
            if expect_prim:
                raise ValueError('malformed query term')
            # push cur into alt, then cross(out, alt), clear alt and cur
            if alt.count >= Q_ATOM_MAX_BOXES:
                raise ValueError('too many alternatives; the cap is %d' % Q_ATOM_MAX_BOXES)
            alt.boxes[alt.count] = cur
            alt.count += 1
            _cross_term(&tmp, out, &alt)
            out[0] = tmp
            memset(&alt, 0, sizeof(wterm_t))
            memset(&cur, 0, sizeof(wbox_t))
            expect_prim = True

    if expect_prim:
        raise ValueError('malformed query term')

    # End of stream: push cur into alt, then cross(out, alt)
    if alt.count >= Q_ATOM_MAX_BOXES:
        raise ValueError('too many alternatives; the cap is %d' % Q_ATOM_MAX_BOXES)
    alt.boxes[alt.count] = cur
    alt.count += 1
    _cross_term(&tmp, out, &alt)
    out[0] = tmp

    # Apply defaults per box, then prune unsatisfiable boxes
    j = 0
    for i in range(out.count):
        box_fill_defaults(&out.boxes[i])
        if not box_unsatisfiable(&out.boxes[i]):
            if j != i:
                out.boxes[j] = out.boxes[i]
            j += 1
    out.count = j

    if out.count == 0:
        raise ValueError('this query term can never match anything')
    return 0


cdef int compile_term(wterm_t *out, qtoken_t *tokens, uint32_t count) except -1:
    """Full pipeline: compile_term_unmerged then boxes_merge."""
    compile_term_unmerged(out, tokens, count)
    boxes_merge(out)
    return 0


def _compile_probe(list ops, bint merge=True):
    """Python test probe: compile a token list, return a list of box dicts."""
    cdef qtoken_t *tokens
    cdef qtoken_t *tok
    cdef wterm_t term
    cdef wbox_t *box
    cdef uint32_t i, j
    cdef list out = [], any_out
    cdef tuple op
    tokens = <qtoken_t *> PyMem_Malloc((len(ops) + 1) * sizeof(qtoken_t))
    if tokens is NULL:
        raise MemoryError()
    try:
        for i in range(len(ops)):
            op = ops[i]
            tok = &tokens[i]
            if op[0] == 'prim':
                tok.opcode = OPC_PRIM
                tok.kind = PRIM_NAMES[op[1]]
                tok.value = op[2]
                tok.negated = op[3]
            else:
                tok.opcode = OPC_NAMES[op[0]]
                tok.kind = 0
                tok.value = 0
                tok.negated = False
        if merge:
            compile_term(&term, tokens, <uint32_t> len(ops))
        else:
            compile_term_unmerged(&term, tokens, <uint32_t> len(ops))
    finally:
        PyMem_Free(tokens)
    for i in range(term.count):
        box = &term.boxes[i]
        any_out = []
        for j in range(box.any_count):
            any_out.append((box.any_word[j], box.any_mask[j]))
        out.append({'neg': (box.neg[0], box.neg[1], box.neg[2], box.neg[3]),
                    'any': tuple(any_out), 'sign': box.sign})
    return out
