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
# The stereo unit table: which atoms and bonds could carry a configuration, and what their
# directions are.  This fragment holds the derived segment and the candidate rules.
# Whether a candidate is stereogenic -- whether flipping it yields a different molecule -- is a
# different question, answered against the automorphism group in `_canonical.pxi`, and it is
# deliberately not asked here.  A candidate is a place a parity CAN be written; nothing in this
# file reads or writes one.
#
# WHY THE TABLE IS DERIVED AND NOT PERSISTENT. A parity is meaningful only relative to an ordered
# list of directions, and this design fixes that list as a function of the anchor's CSR neighbour
# order alone (the exact order is ruling F26, below).  The CSR is persistent, so the list is
# recoverable from the graph at any time and storing it would be storing a second copy of something
# already on disk -- one that a bond edit could leave disagreeing with the first.  Tetrahedral
# stereo therefore costs one byte per atom in SEG_PARITY, and this table is a cache in front of the
# graph.
#
# THE ANCHOR NO-COLLISION INVARIANT, WHICH THREE THINGS DEPEND ON. Every unit is keyed by one
# anchor atom -- the centre for atom kinds, the lower-indexed terminal for bond kinds -- and no
# atom may anchor two units.  That is not an accident of which kinds this file emits; it is what
# the storage assumes:
#
#   1. A parity is keyed by anchor slot in SEG_PARITY, so a second unit on the same anchor would
#      overwrite the first one's configuration.
#   2. `stereo_unit_of` is a lookup BY ANCHOR.  With two units per anchor it would return an
#      arbitrary one of them, and the automorphism filter's per-atom questions would get an
#      arbitrary answer.
#   3. The perception scratch is bounded at one record per atom, so a colliding kind would also
#      be a buffer overrun rather than a wrong answer.
#
# It does NOT hold by itself, and the hybridization case analysis that spec 3.3 first offered for it
# -- sp3 centre, sp2 cis/trans terminal, sp allene centre, an atropisomer pivot whose degree is fully
# spent -- was wrong in two of its three bullets.  Two ordinary molecules break it:
#
#   * a SULFUR cumulene terminal.  `PhS(CH3)=NCH3` reaches four ATOM directions as two sigma, one pi
#     and a lone pair, so it is already SU_TETRA, while its two sigma positions also read as a
#     cis/trans terminal's pair.  Ruling F43 refuses that terminal, on the terminal's own merits: a
#     pyramidal sulfur has no plane for two in-plane positions to lie in.  See `_terminal_pair`.
#   * an ATROPISOMER pivot in a ring of eight or more.  A Kekule aryl pivot always carries a ring
#     double bond, so it is always a cis/trans terminal too; below ring size 8 the small-ring cut
#     drops that chain, and at 8 and up it does not.  Ruling F45 resolves it by ANCHOR CHOICE -- a
#     bond kind may anchor at either end -- rather than by refusing the axis, because refusing a pivot
#     that has a ring double bond would kill every biaryl.  See pass 3.
#
# What holds the invariant up is therefore three local rules, each argued where it is written, plus
# one gate.  `_stereo_emit` is that gate: it refuses a second record on a taken anchor, once, on the
# only path that can violate it, and there is no repair -- a parity is keyed by anchor slot in
# `SEG_PARITY`, so a second record on one anchor would still overwrite the first.  Hence a raise and
# not a fallback.  A later task widening a kind's domain must do what ruling F45 did: ask
# `_anchor_taken` and choose, or refuse locally with a reason.  Amide rotamers are the named
# example -- their pivot is an sp2 carbon that can also be a cis/trans terminal.
#
# WHAT A DIRECTION IS. Four of them, and a candidate needs exactly four:
#
#   * each sigma-bonded neighbour is one direction, its atom slot
#   * each pi-bonded neighbour is ONE direction, not two -- the second lobe of a double bond is
#     not somewhere a substituent can sit
#   * an explicit hydrogen is one direction, NAMED by its atom slot like any other neighbour
#   * each implicit hydrogen is one direction with no atom of its own
#   * a lone pair is one direction with no atom of its own, on sulfur only (spec 4.2, D4)
#   * a triple or an order-8 bond refuses the atom outright: see `_perceive_stereo_units`
#
# WHERE A HYDROGEN DIRECTION SORTS, AND WHY THAT IS THE WHOLE OF THE RULE (ruling F26). A `refs`
# entry is SU_NO_REF if and only if the direction has NO ATOM -- an implicit hydrogen or the lone
# pair.  The order is: heavy slots in CSR ascending order, then hydrogen slots in ascending order,
# then the unnamed directions, implicit hydrogens before the lone pair.
#
# That order is applied per DIRECTION LIST, and a bond kind has two of them (spec 3.2): an atom kind
# packs four directions in one F26 order, a bond kind packs two F26-ordered pairs with the anchor's
# pair first.  So a bond kind's `refs` can hold SU_NO_REF in the middle and need not ascend across
# the pair boundary; `stereo_unit_t` says the same thing at the declaration `translate_stereo`
# reads through.
#
# The property being bought is that a stored parity survives a change of representation.
# Explicitness is a drawing choice and not chemistry -- parsers, standardisation and the depiction
# layer all add and drop explicit hydrogens -- so the list a parity is measured against must not
# be re-based when one appears.  It is the hydrogen direction's POSITION that has to be fixed for
# that, not its identity: CHFClBr gives (F, Cl, Br, none) implied and (F, Cl, Br, H) drawn, the
# three heavy directions in the same three places and the hydrogen in the same place, so the same
# parity still names the same configuration.  The obvious rule -- a hydrogen is an ordinary neighbour
# sorted by its own slot -- loses exactly this, because an H's slot depends on when it was added,
# so it can land anywhere in the list and every explicitation silently re-bases every stored parity
# in the molecule.
#
# Erasing the hydrogen's identity instead would buy the same invariance and cost more: it makes
# `[2H]C([H])(Cl)Br` -- a centre whose only distinguishing feature is which of two hydrogens is
# which -- a record with two identical `none` directions, so it could be called a candidate but
# never assigned a configuration.  Naming them costs nothing and makes it expressible.  The one
# residual: a centre with one drawn and one implied hydrogen orders the two by slot once the second
# is drawn, so its parity depends on the new atom's slot.  In practice the new slot is the higher
# one and nothing moves, and heavy directions are already ordered by slot, so this adds no
# dependence the design did not already have.  There is deliberately no machinery for it.
#
# WHY THE LONE PAIR IS COUNTED FROM AN ELECTRON BUDGET. Sulfur is the only element whose lone pair
# counts, but "sulfur" is not enough to know whether there IS one: a sulfoxide has a pair and
# reaches four directions with two sigma and one pi, while a sulfone has spent it and reaches four
# with two sigma and two pi.  `_sulfur_lone_pair` therefore counts electrons -- six, less the
# charge, less one per bonding electron spent -- which decides both cases with no element table and
# no bond pattern.  Two consequences worth naming: the ylide drawing `[O-][S+](C)C` and the
# double-bond drawing `O=S(C)C` get the SAME verdict, because both spend the same electrons; and a
# sulfonium `[S+](C)(C)C` gets its pair, which is what makes it a centre.
#
# At most ONE lone-pair direction, even where the budget allows two.  Two lone pairs on one atom
# are the same direction twice: nothing can ever tell them apart, so a record that needs both of
# them to reach four -- dimethyl sulfide, two sigma and two pairs -- is not a candidate, and
# counting them both would invent a centre that no later stage could remove, since nothing after
# this point re-examines a rejected atom.
#
# HOW MUCH DISTINGUISHABILITY IS TESTED HERE. Almost none, on purpose.  "Four DISTINGUISHABLE
# directions" is in the end a question about the automorphism group, and the automorphism filter
# can only ever REMOVE units.  So the asymmetry that matters is: admitting a candidate that turns
# out not to be stereogenic costs one filtered record, while refusing one loses it permanently.
# The local test therefore rejects only what it can be certain of -- two pi directions whose
# neighbours are TERMINAL and carry the same atom record, which is the symmetric sulfone -- and
# lets everything else through.  Note the terminal requirement: two nitrogens double-bonded to the
# same sulfur have identical atom records however different their substituents are, so comparing
# records alone would drop a real sulfodiimide centre.  The comparator is `_atom_colour_equal`
# from `_canonical.pxi`, deliberately the same one the automorphism group uses, so this cheap test
# can never disagree with the exact one that follows it.


cdef enum:
    SU_TETRA = 0
    SU_CIS_TRANS = 1
    SU_ALLENE = 2
    SU_ATROPISOMER = 3
    SU_HELICAL = 4          # RESERVED: never produced; see the spec's non-goals

# The same four kinds on the PYTHON surface, because `stereo_units()` hands out `{'kind': <int>}`
# and until now there was no name for the integer it hands out.  Both toolkit converters in
# `chython/chimera/` independently grew their own by-value copy of this enum -- each commented as
# avoiding an import of the extension that `chimera/_iupac.py` makes routinely -- and a by-value
# copy of a domain is a domain that drifts silently, since nothing compares the two.  Published
# from the enum itself so there is one statement of it, beside the field it names (RULES.md 6).
#
# SU_HELICAL is published with the rest even though nothing produces it: a consumer switching on
# `kind` needs to be able to SAY "not this one", and leaving the reserved value nameless is how it
# ends up written as a bare 4 in somebody's else-branch.
#
# THROUGH `globals()` for the same reason `H_UNKNOWN` is: these are C constants, and the obvious
# `SU_TETRA = SU_TETRA` at module scope is a Python name shadowing the C name it means to publish.
# A string key is the one place the C name survives to be read.
globals()['SU_TETRA'] = SU_TETRA
globals()['SU_CIS_TRANS'] = SU_CIS_TRANS
globals()['SU_ALLENE'] = SU_ALLENE
globals()['SU_ATROPISOMER'] = SU_ATROPISOMER
globals()['SU_HELICAL'] = SU_HELICAL


# `stereo_unit_t.spare` is one byte shared between two owners, and the split is the reason nothing
# may ever assign it wholesale:
#   * the LOW nibble is flag bits -- SU_STEREOGENIC = 1.  Values 2, 4, 8 are free for a future task
#     that needs them; a writer must set them with `|=` so that a plain `=` cannot erase the mask above.
#     See the comment on SU_STEREOGENIC about the future value that §4.6 once reserved here.
#   * the HIGH nibble is a 4-BIT MASK OF WHICH `refs` SLOTS HOLD AN UNNAMED DIRECTION, written here,
#     read as `spare >> SU_UNNAMED_SHIFT`.  Bit i set means slot i is a direction with no atom of its
#     own; bit i clear with `refs[i] == SU_NO_REF` means slot i is not a direction at all.  Four
#     slots, four bits, so it cannot overflow the nibble.
#
# WHY A MASK AND NOT A COUNT (ruling F41).  For an atom kind the unnamed directions occupy the tail,
# so a count locates them and the two spellings are equivalent -- `popcount(mask)` is exactly that
# count.  For a BOND kind they sit at per-pair positions and a count is
# ambiguous precisely where the decision lives:
#
#   CH2=CHCH3     two unnamed -> NOT stereogenic: both implicit hydrogens are on one terminal
#   CH3CH=CHCH3   two unnamed -> stereogenic:     one on each terminal
#
# The automorphism filter cannot tell those apart from `2`, and it can from `0b1011` against
# `0b1010`.  As it happens
# the first of those two is refused HERE -- `_terminal_pair` rejects a terminal all of whose
# directions are unnamed, for the same reason it is not stereogenic -- so today the ambiguity is a
# statement about what the field can express rather than about a record that exists.  The oxime below
# is a record that does exist, and the next case is the reason not to bet on the refusal holding: a
# kind whose direction list is not a pair, added later, would need the slots.  The largest
# real population of E/Z units is worse still: `CH3CH=NOH` gives refs `(CH3, None, O, None)`, where
# slot 1 is the carbon's implicit hydrogen -- a real direction -- and slot 3 is the nitrogen's lone
# pair, which `_terminal_pair` deliberately does not count, so it is a slot with no direction at all.
# A count says `1` for both and cannot say which slot.  The mask says `0b0010`.
#
# This is not a format change: SEG_STEREO_UNIT is derived and rebuilt from the graph on every read.
cdef enum:
    SU_FLAG_MASK = 0x0F     # `spare & SU_FLAG_MASK` is the flag nibble
    SU_UNNAMED_SHIFT = 4    # `spare >> SU_UNNAMED_SHIFT` is the slot mask; it is the top nibble
    SU_UNNAMED_MASK = 0x0F  # ...and this is its width, for a reader that masks rather than shifts


cdef enum:
    # The flag nibble's bits.  `mark_stereogenic` sets SU_STEREOGENIC with `|=` -- it is the ONE
    # sanctioned exception to "`_stereo_emit` is the only writer of `spare`", and only for this
    # nibble: a plain assignment would erase the unnamed-direction mask sharing the byte.
    SU_STEREOGENIC = 1      # the unit really is stereogenic: no automorphism is a witness for it
    # Values 2, 4, 8 are free.  §4.6 once allocated value 2 for a geometric-realizability mark
    # (SU_UNREALIZABLE) on small-ring cis/trans units.  That stage proved provably empty on a
    # Kekulé-only arena: `mark_unrealizable` would call `_terminals_share_small_ring` on exactly
    # the units that already survived `_terminals_share_small_ring` in perception's pass 2, so the
    # intersection is always empty.  The full argument is in `_terminals_share_small_ring`'s
    # comment and at the refusal site.  When the arena can distinguish aromatic bonds from Kekulé
    # order 2, value 2 is available for that purpose with the full §4.6 machinery.


# A typed global rather than an enum member, for the same reason as `CANON_NO_SLOT` in
# `_canonical.pxi`: 0xFFFFFFFF is not representable in the `int` a C enumerator has to fit, and
# the comparisons against it happen inside `nogil` where a DEF would be a Python int.
#
# NOTE ON SU_NO_REF IN PERMUTATION MATCHING.  The matching in `translate_stereo` must NOT
# rely on SU_NO_REF == 0xFFFFFFFF sorting last under a numeric comparison.  Ruling F26 puts
# NAMED hydrogens after the heavy atoms by rule, not by value: a drawn hydrogen's slot can be
# any uint32_t and may well be LOWER than a heavy-atom slot added later, so a numeric sort
# gives the wrong reference order the moment an explicit H is involved.  Every SU_NO_REF match
# uses slot EQUALITY and position ORDER inside the array.
cdef uint32_t SU_NO_REF = 0xFFFFFFFF


# ------------------------------------------------------------------------------------------------
# PARITY TRANSLATION.  The permutation table and helpers for `translate_stereo`.
#
# PERM_PARITY_4 is indexed by `_perm_index4`, which encodes a 4-permutation as an integer in
# 0..23 using the factorial number system (Lehmer code / factoriadic):
#
#   index = l0*6 + l1*2 + l2*1
#
# where l_i is the rank of perm[i] among the elements of {0,1,2,3} not yet used by
# positions 0..i-1.  l3 is always 0 and contributes nothing.  The encoding gives a
# bijection between the 24 permutations of {0,1,2,3} and the integers 0..23, and parity
# is determined by the number of inversions modulo 2.
#
# WHY N==4 IS THE ONLY CASE.  Perception's four `_stereo_emit` call sites (TETRA, ALLENE,
# CIS_TRANS, ATROPISOMER) all pass n_refs=4, so u.n_refs is always 4 for every unit that
# reaches `translate_stereo`.  The one site that passes anything else is
# `_stereo_anchor_collision_probe`, whose second emit is intentionally refused and whose unit
# is never translated.  The tetrahedral candidate rule requires exactly four directions, so nothing in
# the design produces a unit with fewer.  The n==3 and n==2 tables are therefore both dead; any
# later kind with a different list length should re-derive its own table rather than try to
# extend this one.
#
# The slice assignment below runs at Python import time (it is a module-scope statement), which
# is a one-time cost.  The contents are a literal, so no computation happens -- just a memcpy
# into the static C array.
#
#   index  permutation  inversions  parity
#       0  (0,1,2,3)        0         0
#       1  (0,1,3,2)        1         1
#       2  (0,2,1,3)        1         1
#       3  (0,2,3,1)        2         0
#       4  (0,3,1,2)        2         0
#       5  (0,3,2,1)        3         1
#       6  (1,0,2,3)        1         1
#       7  (1,0,3,2)        2         0
#       8  (1,2,0,3)        2         0
#       9  (1,2,3,0)        3         1
#      10  (1,3,0,2)        3         1
#      11  (1,3,2,0)        4         0
#      12  (2,0,1,3)        2         0
#      13  (2,0,3,1)        3         1
#      14  (2,1,0,3)        3         1
#      15  (2,1,3,0)        4         0
#      16  (2,3,0,1)        4         0
#      17  (2,3,1,0)        5         1
#      18  (3,0,1,2)        3         1
#      19  (3,0,2,1)        4         0
#      20  (3,1,0,2)        4         0
#      21  (3,1,2,0)        5         1
#      22  (3,2,0,1)        5         1
#      23  (3,2,1,0)        6         0
# NOTE: PERM_PARITY_4, _perm_index4, and permutation_parity_of serve SU_TETRA (atom kind) only.
# Bond-kind units (SU_BOND and similar) use a direct pair-XOR computation in translate_stereo
# and never consult this table.
cdef uint8_t PERM_PARITY_4[24]
PERM_PARITY_4[:] = [0, 1, 1, 0, 0, 1,
                    1, 0, 0, 1, 1, 0,
                    0, 1, 1, 0, 0, 1,
                    1, 0, 0, 1, 1, 0]


cdef inline uint32_t _popcount4(uint32_t x) noexcept nogil:
    """Hamming weight of x in 0..15 (4-bit popcount, branchless parallel).

    Used by `_perm_index4` to count how many values less than perm[i] are still
    available -- i.e. not yet consumed by earlier positions.
    """
    x = x - ((x >> 1) & 0x5u)
    x = (x & 0x3u) + ((x >> 2) & 0x3u)
    return x


cdef inline uint32_t _perm_index4(uint32_t *perm) noexcept nogil:
    """Factorial-base index of a permutation of {0,1,2,3}, in 0..23.

    index = l0*6 + l1*2 + l2*1, where l_i is the rank of perm[i] among
    the elements of {0,1,2,3} not consumed by positions 0..i-1.
    l3 is always 0 and contributes nothing.
    """
    cdef uint32_t used = 0
    cdef uint32_t k, mask
    cdef uint32_t index
    # position 0: rank * 3! = rank * 6
    k = perm[0]
    index = k * 6u          # rank of k in {0,1,2,3} is exactly k (none used yet)
    used = 1u << k
    # position 1: rank * 2! = rank * 2
    k = perm[1]
    mask = (1u << k) - 1u
    index += _popcount4(mask & ~used) * 2u
    used |= 1u << k
    # position 2: rank * 1! = rank * 1
    k = perm[2]
    mask = (1u << k) - 1u
    index += _popcount4(mask & ~used)
    # position 3: always 0, nothing to add
    return index


cdef inline uint8_t permutation_parity_of(uint32_t *perm) noexcept nogil:
    """Parity of the pre-built permutation `perm[0:4]` of {0,1,2,3}: 0 even, 1 odd.

    `perm` must be a valid permutation of {0,1,2,3}; it is the caller's responsibility to
    build it correctly (see `translate_stereo`).  The precondition is hard: an invalid perm
    index that `_perm_index4` computes out-of-range for PERM_PARITY_4 is UB.
    `perm` must be zero-initialised before the caller fills it, so an unmatched slot (which
    should never happen if the caller validated the order) reads PERM_PARITY_4[0] = 0
    rather than crashing; the caller's validation is the contract, the init is the seatbelt.

    n is always 4; see the PERM_PARITY_4 comment above for why.
    """
    return PERM_PARITY_4[_perm_index4(perm)]


cdef uint8_t translate_parity(uint8_t parity, uint32_t *perm) noexcept nogil:
    """The caller's parity re-expressed in the permuted direction order.

    `parity` is the stored three-state value (0 unset, 1 even, 2 odd); `perm[0:4]` is the
    pre-built permutation from `translate_stereo`.  Returns 0 when parity is 0 (unset).

    The parity field is passed as a value, not read from the unit record.  The unit record's
    `parity` field is always 0 (perception never writes it; the true value lives in the
    anchor atom's `SEG_PARITY` byte).  `_stereo_emit` is the only place in the codebase
    that writes the unit record -- passing parity as a value rather than reading the unit
    keeps that invariant intact.

    Stored parity XOR permutation parity = result parity, where XOR is applied to the
    zero-based value (1→0 even, 2→1 odd) and then re-biased to 1/2.
    """
    if parity == 0:
        return 0
    cdef uint8_t pp = permutation_parity_of(perm)
    # (parity - 1) is 0 for even, 1 for odd; XOR with pp, then bias back to 1/2
    return <uint8_t> (((parity - 1) ^ pp) + 1)


# 24 bytes, the same width as `atom_t`, and packed for the same reason every record in this arena
# is: the segment's bytes are read back through this declaration and nothing else, so the layout
# has to be the declaration rather than whatever padding the platform would have chosen.
# `refs` is always four entries whatever `n_refs` says; the unused tail is SU_NO_REF, so a reader
# that ignores `n_refs` sees directions with no atom rather than stale slots.
cdef packed struct stereo_unit_t:
    uint8_t kind            # SU_TETRA .. SU_ATROPISOMER
    uint8_t parity          # RESERVED, always 0.  The parity is a byte per atom in SEG_PARITY;
                            # this field cannot be dropped: the struct is packed to 24 bytes
                            # (1+1+1+1+4+16) and deleting one puts `anchor`, a `uint32_t`, at an
                            # odd offset in every record -- SU_COUNT_HEADER is four words precisely
                            # to keep every record 8-aligned.
    uint8_t n_refs
    uint8_t spare           # low nibble: SU_STEREOGENIC and three free bits, set with |= only.
                            # high nibble: mask of which `refs` slots hold an unnamed direction.
                            # See SU_UNNAMED_SHIFT and ruling F41.
    uint32_t anchor         # atom slot
    uint32_t refs[4]        # the anchor's directions, per kind (spec 3.2):
                            #   ATOM kinds (SU_TETRA) pack four directions in ONE ruling-F26 order:
                            #     heavy slots CSR ascending, then explicit hydrogens ascending, then
                            #     SU_NO_REF for each direction with no atom of its own.
                            #   BOND kinds (SU_CIS_TRANS, SU_ALLENE, SU_ATROPISOMER) pack TWO
                            #     F26-ordered PAIRS, the anchor's end first.  So SU_NO_REF can appear
                            #     in the MIDDLE -- an implicit-hydrogen but-2-ene is (C, None, C,
                            #     None) -- and the four entries need not ascend across the pair
                            #     boundary, since the anchor's pair leads whatever its slots are.


# The segment is [count][truncated][marked][reserved][records...]: the count lives in the segment's
# own first bytes so that `structure_stereo_unit_count` needs no second segment, and a FOUR-word
# header rather than three keeps every record 8-aligned exactly as the original two-word one did.
# This is also why an EMPTY table is sixteen bytes and not zero -- `structure_append` of length 0
# leaves the segment absent, `structure_has` keeps reporting it missing, and perception would
# silently re-run on every read.
#
# The second word was the padding and now carries `mark_stereogenic`'s truncation answer, which is
# the ONE fact about this table that is not derivable from the table itself: every undecided unit of
# a record whose symmetry search ran out of budget is marked stereogenic, and "everything is marked"
# is indistinguishable from a decision unless the table says so.  It lives here rather than in a C
# field on `Structure` because the flag has to survive the same way the table does -- the table is
# built once and read many times, and every one of those reads is entitled to know (ruling F62; it
# reaches Python as `MoleculeContainer.stereo_truncated`).
#
# The THIRD word says whether `mark_stereogenic` has run on this table at all (ruling F70).  It has
# to be a word rather than an inference from the segment's presence, because the segment can now be
# built WITHOUT the marking pass: `ensure_stereo_units_unmarked` serves the callers that read pure
# constitution -- kind, refs, anchor and the unnamed mask -- which are the journal apply, which paid
# two full budgeted stereogenicity searches per edit before the split, and the isomorphism kernel,
# which by ruling F76 asks only what the target states.  Every reader that wants a mark goes
# through `ensure_stereo_units`, which reads this word and marks if it is 0, so "the segment is
# present" and "the marks are decided" are two separate facts and only the second one gates the
# search.  The FOURTH word is reserved and written 0.
#
# Not a format change: the segment is DERIVED, rebuilt from the graph on every read, and import
# clears its table entry outright (the derived-segment list in `_molecule_arena.pxi`).  Both header
# offsets and the record base are computed from this symbol, so widening it moves nothing a
# serialised buffer can see and STRUCT_VERSION does not change.
DEF SU_COUNT_HEADER = 16

cdef str SU_ANCHOR_COLLISION_MSG = (
    'two stereo units claim one anchor atom; a parity is keyed by anchor slot, so the second would overwrite the first')


cdef inline uint32_t structure_stereo_unit_count(Structure structure) noexcept nogil:
    """How many records the table holds. Zero when the table has not been built: the absent
    segment reads off the zero page, which is what makes every loop below safe without a
    `structure_has` guard of its own."""
    return (<uint32_t *> structure.segment(SEG_STEREO_UNIT))[0]


cdef inline bint structure_stereo_truncated(Structure structure) noexcept nogil:
    """Was this table's stereogenicity marking taken conservatively because a search was truncated?

    False for an absent segment, which reads off the zero page -- the same reason
    `structure_stereo_unit_count` needs no `structure_has` guard.
    """
    return (<uint32_t *> structure.segment(SEG_STEREO_UNIT))[1] != 0


cdef inline stereo_unit_t *structure_stereo_units(Structure structure) noexcept nogil:
    """The record array, past the count word. Never dereference further than
    `structure_stereo_unit_count` says: for an absent segment this points into the zero page,
    which is only ZERO_PAGE_SIZE long, and a real table is longer than that at 171 units."""
    return <stereo_unit_t *> (<char *> structure.segment(SEG_STEREO_UNIT) + SU_COUNT_HEADER)


cdef inline int structure_invalidate_stereo_units(Structure structure) except -1:
    """Forget the table so the next reader derives it again.

    The block is RETIRED, not freed (`structure_retire`): it is held one deep and released with the
    arena, so a pointer taken before this call still reads the previous, correct table instead of
    freed memory.  The block is tracked, so the invalidate does not grow the serialised buffer.

    WHAT GOES STALE AND WHAT CANNOT.  The table's SU_STEREOGENIC marks are computed by
    `mark_stereogenic` from the stored parities of OTHER units (`_stereo_consistent`), so a table
    built while a parity was still set can carry a mark the new parities do not justify.  The refs
    themselves are pure constitution and never go stale, and the record's own `parity` field cannot
    go stale either: perception always writes 0 there and `_unit_dict` reads the live value from
    SEG_PARITY instead.  So the marks are the whole of what this discards.

    WHO CALLS IT (ruling F74).  `validate_stereo`, which clones the arena to clear the parity
    it cannot justify -- and `structure_clone` copies every segment, DERIVED SEGMENTS INCLUDED,
    so the clone inherits marks computed against parities the clone has just cleared.  The clone
    copies the derived caches rather than dropping them, and only two of the seven have an `ensure_*`
    guard, so a dropped cache would leave its unguarded readers answering from the zero page -- which
    is what keeps this call necessary.  The journal apply does NOT call it: after ruling F70 the apply
    builds its table through `ensure_stereo_units_unmarked` and the rebuild between harvest and replay
    has dropped the segment anyway, so there is nothing there to invalidate (measured).
    """
    return structure_retire(structure, SEG_STEREO_UNIT)


cdef inline stereo_unit_t *stereo_unit_of(Structure structure, uint32_t slot) noexcept nogil:
    """The unit anchored at atom `slot`, or NULL. Call `ensure_stereo_units` first.

    A linear scan: the table has one record per stereo unit, not per atom, and molecules with
    enough units for this to matter do not exist.  It is a FUNCTION of the slot only because of
    the anchor no-collision invariant above -- with two records per anchor it would return
    whichever came first, which is why that invariant is asserted rather than assumed.
    """
    cdef uint32_t k
    cdef uint32_t count = structure_stereo_unit_count(structure)
    cdef stereo_unit_t *units = structure_stereo_units(structure)
    for k in range(count):
        if units[k].anchor == slot:
            return &units[k]
    return NULL


cdef inline int _stereo_emit(stereo_unit_t *out, Py_ssize_t *count, uint8_t *anchored,
                             uint8_t kind, uint32_t anchor, uint32_t *refs,
                             uint8_t n_refs, uint8_t unnamed) noexcept nogil:
    """Append one record, refusing a second unit on an anchor that already has one.

    Returns 0, or -1 when the anchor is taken -- the caller turns that into a raise.  Every kind
    must come through here; the refusal is the only enforcement of the invariant the storage
    layout depends on, and a kind that appended a record by hand would bypass it.

    `unnamed` is the mask of `refs` slots holding a direction with no atom of its own -- bit i for
    slot i, ruling F41.  It goes in the high nibble of `spare`, leaving the low nibble to the flag
    nibble (of which `SU_STEREOGENIC` is bit 0); it is zero here, so this is the one place `spare`
    may be assigned rather than or-ed.  A consumer that wants a scalar count takes its popcount.
    """
    cdef stereo_unit_t *u
    cdef int k
    if anchored[anchor]:
        return -1
    anchored[anchor] = 1
    u = out + count[0]
    count[0] = count[0] + 1
    u.kind = kind
    u.parity = 0          # always 0: the anchor's SEG_PARITY byte is the parity
    u.n_refs = n_refs
    u.spare = <uint8_t> ((unnamed & SU_UNNAMED_MASK) << SU_UNNAMED_SHIFT)
    u.anchor = anchor
    for k in range(4):
        u.refs[k] = refs[k]
    return 0


cdef inline bint _anchor_taken(uint8_t *anchored, uint32_t slot) noexcept nogil:
    """Has some record already claimed atom `slot` as its anchor?

    The ONE place the no-collision invariant is reasoned about outside `_stereo_emit`, which is the
    single gate that enforces it (ruling F42 deleted the second, post-build scan: it was unreachable
    while `_stereo_emit` is the only emitter, and deleting it changed no test).  A bond kind may
    anchor at either end, so it asks this and relocates rather than colliding -- see pass 3.
    """
    return anchored[slot] != 0


cdef inline int _bond_electrons(halfedge_t *e) noexcept nogil:
    """How many of THIS atom's electrons the bond spends, for `_sulfur_lone_pair`'s budget.

    The order, except for an aromatic bond, which spends ONE -- its sigma -- and leaves the pi to
    the per-atom term the two callers add once (`+ 1` when the atom has any aromatic bond).  Both
    halves are needed and the split is what makes the answer Kekule-independent:

        an aromatic CH in benzene    2 sigma + 1 delocalised pi = 3
        its Kekule twin              order 1 + order 2          = 3
        a ring-fusion carbon         3 sigma + 1                = 4
        its Kekule twin              1 + 1 + 2                  = 4
        thiophene's sulfur           2 sigma + 1                = 3, so 6 - 3 >= 2, pair kept
        its Kekule twin              1 + 1                      = 2, pair kept

    The delocalised pi is ONE electron per atom however many aromatic bonds it has, because that is
    what an atom contributes to a ring's pi system, and adding it per bond would charge a fusion
    carbon three.  Without the correction an aromatic ring atom's `spent` comes out five too high --
    a bare `spent += e.order` charges 4 per aromatic bond -- and `_sulfur_lone_pair` denies a pair to
    sulfurs that have one.  This is the only raw order read in the file: every other test here is on
    degree, hydrogen count, ring membership or `_is_chain_bond`.
    """
    return 1 if e.flags & HE_AROMATIC else <int> e.order


cdef inline bint _sulfur_lone_pair(atom_t *a, int spent) noexcept nogil:
    """Does this sulfur still have a lone pair, given the bonding electrons it has `spent`?

    Six valence electrons, less the formal charge, less one for each electron in a bond to a
    neighbour or to a hydrogen -- so a sigma bond and a hydrogen cost one each and a double bond
    costs two.  Two electrons left is a pair.  Counting electrons rather than matching a bond
    pattern is what makes the sulfoxide / sulfone split fall out of the rule, and what makes the
    S=O and [S+]-[O-] drawings of one sulfoxide agree.

    Radicals are DELIBERATELY outside the budget: `at_radical` is not consulted, so an odd
    electron count reports a full pair for what is really a pair plus an unpaired electron, or an
    unpaired electron alone.  `[S.](C)(C)C` is therefore admitted as a candidate -- which is the
    right answer, a sulfuranyl radical is pyramidal, but it is the right answer for the wrong
    reason and a later task that wants radicals treated exactly must change this line, not add a
    case around it.
    """
    return (6 - <int> a.charge - spent) >= 2


cdef inline bint _pi_directions_indistinguishable(atom_t *p, atom_t *q) noexcept nogil:
    """Are these two pi-bonded neighbours the same direction twice -- O=S=O?

    Only when both are TERMINAL, so that the atom record is the whole of what hangs off the
    direction, and the records match.  A non-terminal pair is left to the automorphism group; see
    the fragment comment on why refusing too much here is the expensive mistake.

    Named for the case it was written for, but the test is direction-kind agnostic -- "both ends
    are single atoms and those atoms are identical" -- so the cumulene and atropisomer rules below
    use it for their sigma pairs too.  Keeping one comparator is the point: it is
    `_atom_colour_equal`, the same one the automorphism group uses, so no cheap local verdict here
    can contradict the automorphism filter's exact one.
    """
    return p.degree == 1 and q.degree == 1 and _atom_colour_equal(p, q)


# ------------------------------------------------------------------------------------------------
# CUMULENES.  A maximal chain of consecutive double bonds; odd atom count is axial (an allene),
# even is cis/trans-like.  Three things about the walk are not obvious.
#
# WHY AN AROMATIC RING IS THE TRAP HERE, AND WHAT `HE_AROMATIC` NOW SETTLES (ruling F31, revised).
# The arena stores order 4 with HE_AROMATIC set, so `_is_chain_bond` -- written for this and until
# now excluding nothing -- has become the real discriminator: an aromatic-written ring contributes
# NO chain bonds, offers no cis/trans candidates, and cannot collide with an atropisomer unit on a
# pivot.  Ruling F100's order-dependence (79 of 300 creation orders) was a Kekule artifact of exactly
# this trap and is gone for an aromatic-written molecule.
#
# The trap remains for a KEKULE-written one, and its observable shape is not the one it reads like.
# Kekulisation ALTERNATES, so benzene's order-2 half-edges are three separate two-atom chains rather
# than one six-atom chain -- a walk keyed on `order == 2` does not run round the ring.  What it does
# instead is offer all three of those chains as cis/trans candidates, and every terminal passes the
# local test (one ring neighbour and one hydrogen are plainly different), so Kekule benzene reports
# three units, toluene four, and a Kekule biaryl reports one on every pivot -- which then collides
# with the atropisomer unit on the same atom.  chython 2 offers the same candidates (spec 10.6: six
# spurious ones on VS055).  Both spellings are storable by design (a file's bonds are stored as the file
# drew them), so both behaviours are live and the small-ring cut below is still load-bearing.
#
# WHY THE RING RULE IS A CANDIDATE RULE HERE AND NOT ONLY A REALIZABILITY FLAG, AND WHY IT IS FOR
# EVEN CHAINS ONLY (ruling F44).  An EVEN chain whose two terminals share a ring smaller than
# SU_MIN_STEREO_RING is not a candidate: the ring path holds the terminals' substituents cis, so
# there is no second configuration for a parity to name.  That threshold is the spec's 4.6
# realizability number, and the honest reading is that on a Kekule-only arena it is ALSO the only
# available spelling of "this double bond is an aromatic ring's, not a stereogenic one".
#
# The argument is about CIS/TRANS and does not transfer to an axial unit, which is why the cut sits
# below the odd/even split.  An allene's two terminals are perpendicular; "cis" is not defined for
# them and the two configurations are enantiomers that no ring path can equate.  Applying it to odd
# chains cost 1,2-cyclohexadiene and 1,2-cycloheptadiene their axial candidates -- strained but
# chirally distinct, with a literature on enantioselective trapping -- while 1,2-cyclononadiene, the
# textbook resolved cyclic allene, survived: wrong on three of four canonical members of the exact
# axis this epic is judged on.  Three consequences, all deliberate:
#
#   * rings of 8 and up are admitted for cis/trans too, so trans-cyclooctene and the macrocyclic
#     cumulenes the epic pins as a gate survive.  Excluding every ring double bond instead would
#     have been simpler and would have lost them permanently, since this is the only stage that can
#     ADMIT a unit.
#   * a later realizability stage must not apply the same cut twice, and it will find nothing left
#     to suppress for small rings.  Spec 4.6 says so.
#   * a small-ring ALLENE is emitted here and may well be unrealizable on geometric grounds -- a
#     three-ring cannot hold one.  `SU_UNREALIZABLE` (value 2 in the flag nibble) is the mechanism
#     for that question; this build does not take that stage, and the flag nibble's own enum comment
#     on value 2 records why.
#
# WHY THE WALK TERMINATES.  It starts only from an atom with exactly one chain bond and refuses to
# step onto an atom with more than two, so every walked vertex has chain-degree at most two: the
# component of a chain-degree-1 vertex is then a simple path, and the walk cannot revisit.  A
# forged double-bond ring has no degree-1 vertex at all and is never entered; a forged T-junction
# is refused at the step onto it.  Nothing here relies on a step counter.


cdef enum:
    # The file's ONLY ring-size number: below this a shared ring holds a cis/trans unit's
    # substituents cis, so the unit is not admitted.  Editing this line moves the perception
    # boundary AND the boundary that `test_the_small_ring_cut_admits_at_threshold_and_refuses_below`
    # directly measures; the two move together because they share this one constant.
    #
    # WHY THIS IS A CANDIDATE RULE AND NOT A LATER MARKING STAGE (ruling F63, spec §4.6).  §4.6
    # specified a separate realizability mark (`SU_UNREALIZABLE`) set after stereogenicity is
    # decided.  That stage is provably empty on a Kekulé-only arena: it would call
    # `_terminals_share_small_ring` on exactly the units that already survived
    # `_terminals_share_small_ring` here -- same predicate, same inputs, empty intersection.  The
    # reason the candidate rule is the ONLY place the test can live: on a Kekulé arena an aromatic
    # ring bond is order 2 and indistinguishable from a small-ring alkene's, so admitting the unit
    # then marking it unrealizable would also admit every benzene ring as a cis/trans candidate,
    # flooding `stereo_units()` with junk on every aromatic molecule.  The refusal here is therefore
    # correct and permanent on a Kekulé arena; it is NOT conflating "not stereogenic" with "not
    # realizable" -- a small-ring alkene IS genuinely stereogenic in the graph-theoretic sense, and
    # a later marking stage could record that, but on a Kekulé arena it cannot do so without also
    # admitting the aromatic ring bonds.
    #
    # WHEN TO REVISIT: THE PRECONDITION HAS NOW BEEN MET, AND THE WORK IS DELIBERATELY NOT TAKEN.
    # The arena stores order 4 with HE_AROMATIC live, so `_is_chain_bond` really does exclude an
    # aromatic ring bond and the candidate rule could be relaxed to admit non-aromatic small-ring
    # alkenes, with a separate `SU_UNREALIZABLE` mark (value 2 in the flag nibble, still free)
    # carrying what this cut currently conflates.  It is not done here because it CHANGES ANSWERS --
    # cyclohexene would gain a unit -- and a perception change of that size is not something to land
    # beside a storage change.  Note also that the cut cannot simply be deleted even then: a
    # Kekulé-written aromatic ring is still storable and still reaches this test.
    # `test_cyclohexene_double_bond_is_excluded_by_small_ring_cut` and the boundary test name the
    # distinction so that it survives to that future task.
    SU_MIN_STEREO_RING = 8


cdef inline bint _is_chain_bond(halfedge_t *e) noexcept nogil:
    """Is this half-edge a link in a cumulene chain?

    Kekule order 2 and NOT flagged aromatic, and both halves are now live: a bond written aromatic is
    stored as order 4 with HE_AROMATIC, so the flag test excludes it and the order test would have
    too.  Keeping both is not redundancy -- the pair is an invariant `structure_from_bytes` enforces,
    and this predicate reads the half of it that states the INTENT ("not a localised double bond")
    rather than the half that states the encoding.

    A DECISION, not a side effect: an aromatic-written cis/trans unit therefore does not exist.  A
    stereocentre spelled across a bond the file drew as aromatic names a configuration in a
    delocalised system, which is either an aromatic ring's (where there is no configuration to name)
    or a mis-drawn double bond.  Neither is a unit this pass should invent; whoever wants the second
    one back must kekulise first, which is what `kekule()` is for.
    """
    return e.order == 2 and not (e.flags & HE_AROMATIC)


cdef inline int _chain_degree(uint32_t *ptr, halfedge_t *edges, uint32_t i) noexcept nogil:
    """How many chain bonds atom `i` has: 1 makes it a terminal, 2 an interior, 0 neither."""
    cdef uint32_t k
    cdef int n = 0
    for k in range(ptr[i], ptr[i + 1]):
        if _is_chain_bond(&edges[k]):
            n += 1
    return n


cdef inline uint32_t _chain_next(uint32_t *ptr, halfedge_t *edges, uint32_t cur,
                                 uint32_t prev) noexcept nogil:
    """The next atom along the chain from `cur`, arriving from `prev`; SU_NO_REF at the end.
    Pass SU_NO_REF as `prev` for the first step: no atom slot can equal it."""
    cdef uint32_t k
    for k in range(ptr[cur], ptr[cur + 1]):
        if _is_chain_bond(&edges[k]) and edges[k].to != prev:
            return edges[k].to
    return SU_NO_REF


cdef inline bint _cumulene_walk(atom_t *atoms, uint32_t *ptr, halfedge_t *edges, uint32_t t,
                                uint32_t *other, uint32_t *other_prev,
                                uint32_t *n_atoms) noexcept nogil:
    """Walk the chain from terminal `t` to its far end.  False when the chain is not one.

    Reports the far terminal, the atom the walk arrived at it from (that terminal's own chain
    neighbour, which its direction pair has to exclude) and the chain's ATOM count, which is what
    decides axial against cis/trans.
    """
    cdef uint32_t prev = SU_NO_REF
    cdef uint32_t cur = t
    cdef uint32_t nxt
    cdef uint32_t count = 1
    cdef int cd
    while True:
        nxt = _chain_next(ptr, edges, cur, prev)
        if nxt == SU_NO_REF:
            break
        cd = _chain_degree(ptr, edges, nxt)
        if cd > 2:
            return False        # a branched double-bond subgraph is not a chain; also what
                                # bounds this loop -- see the fragment comment
        if cd == 2 and atoms[nxt].degree != 2:
            return False        # a cumulene INTERIOR is sp: exactly two neighbours, no more
        prev = cur
        cur = nxt
        count += 1
    other[0] = cur
    other_prev[0] = prev
    n_atoms[0] = count
    return True


cdef inline uint32_t _chain_nth(uint32_t *ptr, halfedge_t *edges, uint32_t t,
                                uint32_t steps) noexcept nogil:
    """The atom `steps` chain bonds along from terminal `t`.  Used for the axial anchor, which is
    the chain's centre atom; the walk that found the chain has already proved the path exists."""
    cdef uint32_t prev = SU_NO_REF
    cdef uint32_t cur = t
    cdef uint32_t nxt
    cdef uint32_t left = steps
    while left:
        nxt = _chain_next(ptr, edges, cur, prev)
        prev = cur
        cur = nxt
        left -= 1
    return cur


cdef inline uint32_t _ring_prototype_size(Structure structure, uint32_t words, uint32_t word,
                                         uint64_t mask) noexcept nogil:
    """The ring size of the prototype whose bit is `mask` in word `word` of the ring bitmap.

    A prototype's size IS the number of atoms carrying its bit: `_fill_descriptors` sets the bit on
    every atom of the cycle and on no other, so counting them recovers `psize` exactly with no new
    storage and no new segment field.  One pass over atoms, run only for a prototype the two
    terminals actually share -- which is at most a handful of prototypes on any real molecule.
    """
    cdef uint64_t *bits = structure_ring_bits(structure)
    cdef uint32_t i
    cdef uint32_t size = 0
    for i in range(structure.header.atom_count):
        if bits[<size_t> i * words + word] & mask:
            size += 1
    return size


cdef inline bint _terminals_share_small_ring(Structure structure, atom_t *atoms,
                                             uint32_t a, uint32_t b) noexcept nogil:
    """Do the chain's two terminals sit together in a ring smaller than SU_MIN_STEREO_RING?

    The size test and the identity test are ONE test, deliberately: a shared prototype is found and
    then that prototype's own size is measured.  Asking them independently -- "is either atom on some
    small ring" and separately "do they share some ring" -- refuses a twelve-ring double bond with a
    cyclopropane fused at each terminal, where the cyclopropanes constrain nothing about the
    twelve-ring.  Losing a candidate is the expensive
    direction (see the fragment header), so the exact test is worth its one pass over atoms.
    `structure_shares_ring` is the cheap early-out that keeps that pass off the common case, and it
    is also what keeps bicyclohexylidene -- a double bond joining two ring atoms without being in a
    ring itself -- out of the exclusion.

    In the other direction the rule still over-emits: an eighteen-membered alternating macrocycle
    reports nine cis/trans units, none of which is realizable as two configurations.  That half is
    harmless -- stereogenicity marking can only remove a candidate, never add one, and after ruling
    F63 there is no later realizability stage at all -- and stays.
    """
    cdef uint32_t words, k, bit
    cdef uint64_t *bits
    cdef uint64_t shared
    if not at_in_ring(&atoms[a]) or not at_in_ring(&atoms[b]):
        return False
    if not structure_shares_ring(structure, a, b):
        return False
    words = structure_ring_words(structure)
    bits = structure_ring_bits(structure)
    for k in range(words):
        shared = bits[<size_t> a * words + k] & bits[<size_t> b * words + k]
        while shared:
            bit = <uint32_t> _lo_bit64(shared)
            if _ring_prototype_size(structure, words, k,
                                    <uint64_t> 1 << bit) < SU_MIN_STEREO_RING:
                return True
            shared &= shared - 1
    return False


cdef inline bint _terminal_pair(atom_t *atoms, uint32_t *ptr, halfedge_t *edges, uint32_t t,
                                uint32_t chain_nb, uint32_t *out, int *unnamed) noexcept nogil:
    """Fill `out[0:2]` with terminal `t`'s two directions besides the chain and `unnamed` with the
    pair's 2-bit unnamed-slot mask.  False when `t` cannot be a cumulene terminal.

    Ruling F26's order within the pair: the heavy slot first, then an explicit hydrogen's slot,
    then SU_NO_REF.  A pair is NEVER sorted by raw value -- SU_NO_REF being 0xFFFFFFFF makes that
    look right until a drawn hydrogen has the lower slot, and then it silently re-bases the parity
    of every explicitated cumulene.

    `unnamed` is the mask of pair slots holding a direction with no atom of its own -- bit 0 for
    `out[0]`, bit 1 for `out[1]`, shifted by the caller for the far terminal.  An empty slot is NOT
    an unnamed direction: there is no direction there at all, and ruling F41 is that the record has
    to tell those two apart, because `CH3CH=NOH` puts a real implicit hydrogen in one slot and a
    nothing in another and a scalar count gives both the same answer.

    A lone pair is not one of the two directions -- a cumulene terminal's directions are its two
    IN-PLANE SIGMA positions, and spec 4.2's lone-pair direction completes an ATOM's four, not a
    terminal's two.  Ruling F43: that same argument DISQUALIFIES the terminal outright, and it is a
    statement about the terminal and not a dodge of the anchor collision that follows from it.  A
    pyramidal sulfur has no plane for two in-plane positions to lie in: `R2S=NR` and `R2S=CR2` hold
    their configuration AT SULFUR -- a sulfimide or ylide stereocentre, in this epic's scope and
    already perceived as SU_TETRA -- and a putative E/Z across the `S=N` names no information that
    pyramidal unit does not already carry.  So the cis/trans unit is not relocated; it does not exist.

    Refusing on `element == 16` is airtight rather than a hack, and the arithmetic is why: a terminal
    allows `n_sigma + n_h + n_imp <= 2` while an ATOM needs four directions, so the fourth is
    reachable only through a lone pair, and only sulfur's counts (spec 4.2, D4).  The electron budget
    is passed the same `spent` the tetrahedral walk passes, so the two agree by construction.
    """
    cdef uint32_t k
    cdef halfedge_t *e
    cdef uint32_t heavy[2]
    cdef uint32_t hydro[2]
    cdef int n_heavy = 0
    cdef int n_h = 0
    cdef int n_imp
    cdef int spent = 0
    cdef bint any_aromatic = False
    cdef int mask = 0
    for k in range(ptr[t], ptr[t + 1]):
        e = &edges[k]
        # every bond spends electrons, the chain bond included, so this runs before the skip below
        spent += _bond_electrons(e)
        if e.flags & HE_AROMATIC:
            any_aromatic = True
        if e.to == chain_nb and _is_chain_bond(e):
            continue                      # the chain bond is not one of the two directions
        if e.order == 3 or e.order == 8:
            # A terminal carrying a triple bond is not sp2, and order 8 carries no geometry --
            # the same refusal the tetrahedral walk makes, for the same two reasons.
            return False
        # BOTH bounds checks are load-bearing, exactly as in `_perceive_stereo_units`: these are
        # two-element STACK arrays and the count that rejects an over-substituted terminal is only
        # known after this loop has finished writing.  A forged carbon with twenty-four fluorines
        # overruns `heavy` by twenty-two words unless the write is guarded.  The counters keep
        # incrementing past two; it is the count that must stay honest, not the array.
        if atoms[e.to].element == 1:
            if n_h < 2:
                hydro[n_h] = e.to
            n_h += 1
        else:
            if n_heavy < 2:
                heavy[n_heavy] = e.to
            n_heavy += 1
    n_imp = at_implicit_h(&atoms[t])
    if n_imp == H_UNKNOWN:
        # A TERMINAL WHOSE HYDROGEN COUNT NOBODY RECORDED, refused only where the count could have
        # changed the answer.  When it can, the refusal is about chemistry before it is about
        # arithmetic: the missing number is exactly the one that decides whether this terminal has an
        # isomer at all.  One hydrogen and one heavy substituent is a genuine E/Z pair; two hydrogens
        # is `=CH2` and has none (the `n_heavy + n_h == 0` case below).  A record that does not say
        # which cannot be assigned a configuration, and guessing would put a parity on half the
        # ethenes in a badly written file.
        #
        # TWO NAMED DIRECTIONS ALREADY FILL THE TERMINAL, so there the missing number cannot matter --
        # a terminal has two in-plane positions and both are spoken for.  Refusing anyway would lose
        # the stated geometry of a fully substituted double bond or allene terminus, as in the
        # four-direction case in `_perceive_stereo_units`, so the sentinel reads as zero here for the
        # same reason: any other value is refused by `> 2` below whatever it is.
        #
        # The sentinel is never added to `spent`, which is why this test sits above the arithmetic:
        # 15 in the electron budget makes the sulfur lone-pair test read a corrupted total.
        if n_heavy + n_h < 2:
            return False
        n_imp = 0
    spent += n_imp
    if any_aromatic:
        spent += 1                        # the one delocalised pi electron; see `_bond_electrons`
    if atoms[t].element == 16 and _sulfur_lone_pair(&atoms[t], spent):
        # Ruling F43, argued in full above: a pyramidal sulfur has no in-plane pair of positions.
        # This is also the one shape where the anchor no-collision invariant would otherwise be
        # false -- such a sulfur reaches four directions as two sigma, one pi and the pair, so
        # SU_TETRA is already anchored on it, and `PhS(CH3)=NCH3` raised out of
        # `stereo_units()` whenever the sulfur happened to hold the lower of the two terminal slots.
        return False
    if n_heavy + n_h + n_imp > 2:
        return False                      # more than two in-plane directions is not a terminal
    if n_heavy + n_h == 0:
        # No NAMED direction at all, so nothing a parity could be measured against: either the
        # terminal is a =CH2, whose two implicit hydrogens are both protium and therefore one
        # direction twice over -- ethene, propene, isobutene, none of which has a cis/trans isomer
        # -- or it is a forged =C with no directions whatsoever.  This is the same rejection the
        # tetrahedral walk makes for two sulfur lone pairs, and it is safe in the strong sense:
        # a pair of implicit hydrogens can never become distinguishable later, so nothing is
        # lost permanently.
        return False
    if n_heavy == 2 and _pi_directions_indistinguishable(&atoms[heavy[0]], &atoms[heavy[1]]):
        return False                      # (CH3)2C= : one direction twice over
    out[0] = SU_NO_REF
    out[1] = SU_NO_REF
    # `n_heavy + n_h <= 2` past the count test above, so neither fill can leave the window.
    for k in range(<uint32_t> n_heavy):
        out[k] = heavy[k]
    for k in range(<uint32_t> n_h):
        out[n_heavy + k] = hydro[k]
    # The implicit hydrogens take the slots after every named one, so their bits are the tail of the
    # pair.  At most one of them is reachable: two would need `n_heavy + n_h == 0`, refused above.
    for k in range(<uint32_t> n_imp):
        mask |= 1 << (n_heavy + n_h + <int> k)
    unnamed[0] = mask
    return True


# ------------------------------------------------------------------------------------------------
# ATROPISOMERS.  Spec 4.5's rule, and the ONLY heuristic in the design -- everything else here is a
# graph property.  It lives in `_is_atropisomer_axis` alone so that the barrier it stands in for can
# be retuned without touching perception.
#
# The rule: the bond is single and not in a ring, both ends are ring atoms, each end's ring carries
# at least one ortho substituent, and the two ends' ortho pairs are distinguishable.  One condition
# is not in the spec's sentence and is load-bearing anyway: NEITHER PIVOT MAY CARRY A HYDROGEN
# DIRECTION.  Spec 3.3's case analysis argues the anchor cannot collide because a pivot's degree 3
# is "fully consumed by two ring bonds plus the pivot -- no room for an exocyclic double bond".
# That rules out a cumulene collision and not a tetrahedral one: a saturated pivot -- bicyclohexyl,
# not biphenyl -- has two ring bonds, the pivot bond and an implicit hydrogen, which is four
# directions and therefore already a tetrahedral unit on that same atom.  Perception would raise on
# an ordinary molecule.  Requiring a hydrogen-free pivot is also the right chemistry: a saturated
# C-C bond rotates however crowded its ortho positions are, and it is the aryl ends' rigidity that
# makes the biaryl axis a configuration rather than a conformation.
#
# The other half of the collision -- the pivot's own ring double bond, which a Kekule aryl pivot
# ALWAYS has -- is not handled here at all.  Below ring size 8 the cumulene small-ring cut happens to
# drop that chain, which is why plain biphenyl never collided; at 8 and up it does not, and an
# ortho-substituted biaryl of eight-membered rings would raise.  Ruling F45 puts the resolution in
# pass 3's anchor choice instead, so this end of the rule does not depend on the cut's threshold.
# Refusing a pivot with a ring double bond is wrong twice over: it kills every biaryl, and it makes a
# candidate rule depend on a realizability number.
#
# An explicit hydrogen is NOT an ortho substituent.  A hydrogen is not what hinders rotation, and
# explicitness is a drawing choice -- counting a drawn one would make a molecule an atropisomer and
# its implicit-hydrogen twin not one, which is the representation dependence ruling F26 exists to
# keep out of the record.


cdef inline bint _is_ring_fusion_atom(atom_t *atoms, uint32_t *ptr, halfedge_t *edges,
                                      uint32_t o) noexcept nogil:
    """Is `o` a ring-fusion atom -- degree three or more with every one of its bonds in a ring?

    Ruling F46: such an ortho neighbour hinders rotation about a biaryl axis the way a substituent
    does, because the ring fused there occupies the ortho position.  1,1'-binaphthyl is the case that
    forces the rule: its hindrance is entirely the PERI hydrogen on C8, so nothing hangs off C8a as an
    exocyclic substituent and the plain ortho test perceived nothing at all -- while BINOL, the same
    scaffold plus two hydroxyls, was already a candidate.  Binaphthyl and BINAP are most of why
    `kind = 3` exists, so half the kind was invisible.

    A separate predicate rather than another clause inside the ortho loop, because it is a different
    argument about a different atom: the loop asks what hangs OFF the ring, this asks what the ring
    is fused TO.
    """
    cdef uint32_t k
    if atoms[o].degree < 3:
        return False
    for k in range(ptr[o], ptr[o + 1]):
        if not (edges[k].flags & HE_IN_RING):
            return False
    return True


cdef inline bint _atropisomer_end(atom_t *atoms, uint32_t *ptr, halfedge_t *edges, uint32_t pivot,
                                  uint32_t partner, uint32_t *out) noexcept nogil:
    """One end of a candidate axis: fill `out[0:2]` with the pivot's two ring directions in CSR
    ascending order and answer whether the end qualifies.

    No bounds guard on `out`, and that is a precondition rather than an omission: the degree test
    below runs BEFORE the loop, so exactly two of the pivot's three half-edges are not the partner
    bond.  The tetrahedral perception's guards are needed because there the count is only
    known afterwards.

    THE HYDROGEN TEST IS A TRUTHINESS TEST AND IS CORRECT ON THE SENTINEL, which is why this one
    needed no code.  It reads "the pivot carries no hydrogen direction", and H_UNKNOWN (15) is
    truthy, so an axis whose pivot has no recorded hydrogen count is refused.  That is the strict
    direction and the one we want: an unrecorded count might be a hydrogen, and a hydrogen-bearing
    pivot is a saturated centre that rotates (see the ATROPISOMERS block above), so admitting it on
    the strength of a missing number would perceive an axis through a single bond that turns freely.
    """
    cdef uint32_t k, m, o
    cdef halfedge_t *e
    cdef uint32_t ring[2]
    cdef int n_ring = 0
    cdef bint hindered = False
    if atoms[pivot].degree != 3 or at_implicit_h(&atoms[pivot]) or at_explicit_h(&atoms[pivot]):
        return False
    for k in range(ptr[pivot], ptr[pivot + 1]):
        e = &edges[k]
        if e.to == partner:
            continue
        if not (e.flags & HE_IN_RING):
            return False        # a third acyclic bond: this is a branch point, not a biaryl pivot
        ring[n_ring] = e.to
        n_ring += 1
    for m in range(2):
        o = ring[m]
        if _is_ring_fusion_atom(atoms, ptr, edges, o):
            hindered = True         # ruling F46: the fused ring is itself the ortho substituent
            break
        for k in range(ptr[o], ptr[o + 1]):
            e = &edges[k]
            if e.to == pivot or (e.flags & HE_IN_RING):
                continue        # inside the ring, so not a substituent hanging off it
            if atoms[e.to].element == 1:
                continue        # a drawn hydrogen is not what hinders rotation
            hindered = True
            break
        if hindered:
            break
    if not hindered:
        return False
    # The rule's distinguishability clause, spelled with the SHARED comparator so that no local
    # verdict here can contradict the automorphism filter's exact one.  It is UNREACHABLE BY
    # CONSTRUCTION, not merely unreached by today's molecules: the comparator requires `degree == 1`
    # on both arguments and both of these are ring atoms, whose degree is at least two, so it can
    # never fire for any comparator that keeps that terminal precondition -- and the precondition is
    # load-bearing where the comparator is used on pi directions.  The line is therefore dead until
    # the automorphism filter either widens the comparator to non-terminal pairs or replaces this
    # call with the automorphism test, which is what actually decides it: the ortho pair of a
    # 2,6-disubstituted end is symmetric, and the evidence is an automorphism swapping ring[0] with
    # ring[1].  It is kept, rather than deleted, because it is the rule as spec 4.5 words it and
    # because deleting it would hide the obligation.
    # `test_a_symmetric_ortho_pair_is_still_a_candidate_here` pins the current answer so that a
    # later automorphism edit is a visible change.
    if _pi_directions_indistinguishable(&atoms[ring[0]], &atoms[ring[1]]):
        return False
    out[0] = ring[0]
    out[1] = ring[1]
    return True


cdef bint _is_atropisomer_axis(Structure structure, halfedge_t *e, uint32_t a,
                               uint32_t b) noexcept nogil:
    """Is the bond `a`-`b`, whose half-edge is `e`, an atropisomer axis?  Spec 4.5's heuristic, in
    one place.

    Symmetric in `a` and `b`, so a caller may pass the bond either way round; the anchor rule is the
    caller's, not this predicate's.  The half-edge is an argument rather than something this looks up
    with `csr_find_at`, because the only caller is a sweep that already holds it and this predicate
    is invoked on every bond in the molecule.

    KNOWN OVER-ADMISSION.  The ortho test below counts a ring-FUSION neighbour as hindering
    (ruling F46, which is what makes 1,1'-binaphthyl a candidate at all), and a fusion atom is not
    always a peri position -- some fused biaryls are admitted whose real rotational barrier is low.
    Over-admission is perception's safe direction: `mark_stereogenic` can remove a candidate;
    nothing after this point re-examines a refused one.
    """
    cdef atom_t *atoms = structure.atoms()
    cdef uint32_t *ptr = csr_ptr(structure)
    cdef halfedge_t *edges = csr_edges(structure)
    cdef uint32_t refs_a[2]
    cdef uint32_t refs_b[2]
    if e.order != 1 or (e.flags & HE_IN_RING):
        return False
    if not at_in_ring(&atoms[a]) or not at_in_ring(&atoms[b]):
        return False
    return (_atropisomer_end(atoms, ptr, edges, a, b, refs_a)
            and _atropisomer_end(atoms, ptr, edges, b, a, refs_b))


cdef Py_ssize_t _perceive_stereo_units(Structure structure, stereo_unit_t *out,
                                       uint8_t *anchored) noexcept nogil:
    """Fill `out` with the molecule's stereo units and return how many there are; -1 on a
    collision (see `_stereo_emit`).

    `out` must have room for one record per atom and `anchored` must be a zeroed byte per atom.

    Three passes, one per kind that has its own walk: tetrahedral over atoms, cumulene over
    double-bond chains, atropisomer over bonds.  THE RESULTING UNIT ORDER IS NOT PROMISED and
    nothing may come to depend on it -- the tetrahedral records happen to be in ascending anchor
    order because that pass is a single ascending sweep, and `stereo_unit_of` is a linear scan
    precisely so that the cumulene and atropisomer passes need not preserve it.
    """
    cdef atom_t *atoms = structure.atoms()
    cdef uint32_t *ptr = csr_ptr(structure)
    cdef halfedge_t *edges = csr_edges(structure)
    cdef atom_t *a
    cdef halfedge_t *e
    cdef uint32_t n = structure.header.atom_count
    cdef uint32_t i, k, pi_a, pi_b
    # zeroed rather than left to the walk that fills them: Cython cannot see that `_cumulene_walk`
    # writes all three whenever it returns True, and an uninitialised read is the kind of thing a
    # `nogil` warning is worth keeping clean for
    cdef uint32_t far = 0, far_prev = 0, chain_len = 0
    cdef uint32_t refs[4]
    cdef uint32_t h_refs[4]
    cdef int n_ref, n_h, n_unnamed, spent, pi_count, unnamed_mask, near_mask, far_mask
    cdef uint32_t pivot, other_pivot
    cdef bint refused, any_aromatic
    cdef Py_ssize_t count = 0

    for i in range(n):
        a = &atoms[i]
        n_ref = 0
        n_h = 0
        n_unnamed = 0
        spent = 0
        pi_count = 0
        pi_a = 0
        pi_b = 0
        refused = False
        any_aromatic = False
        for k in range(ptr[i], ptr[i + 1]):
            e = &edges[k]
            # A triple bond, and chython's order 8 for everything else, refuse the atom outright.
            # For a triple bond the valence is already spent, so no legal record reaches four
            # directions past one -- the early-out only ever fires on a forged one, and refusing
            # that is also the answer we want.  Order 8 carries no geometry to build a direction
            # on: it stands in for dative, ionic and unspecified bonds alike.
            if e.order == 3 or e.order == 8:
                refused = True
                break
            # ORDER 4 IS NOT REFUSED HERE and must not be: an aromatic bond is a direction like any
            # other, and a benzylic centre reaches four directions through one.  What it changes is
            # the electron budget, which `_bond_electrons` carries.
            spent += _bond_electrons(e)
            if e.flags & HE_AROMATIC:
                any_aromatic = True
            # BOTH bounds checks below are load-bearing: `refs` and `h_refs` are four-element
            # STACK arrays, while the count that rejects a hypercoordinate record is only computed
            # after this loop has finished writing into them.  So a phosphorus with seven fluorines
            # overruns `refs` by three, and a carbon drawn with seven hydrogens overruns `h_refs`,
            # unless the write is guarded here.  The counters keep incrementing past four; it is
            # the count that must stay honest, not the array.
            if atoms[e.to].element == 1:
                # an explicit hydrogen is a NAMED direction (ruling F26), collected separately so
                # that it can sort after every heavy one.  It is deliberately kept out of the pi
                # bookkeeping: nothing hangs off a hydrogen, so it can never be the repeated pi
                # direction the rejection below looks for.
                if n_h < 4:
                    h_refs[n_h] = e.to
                n_h += 1
            else:
                if e.order == 2:
                    pi_count += 1
                    if pi_count == 1:
                        pi_a = e.to
                    elif pi_count == 2:
                        pi_b = e.to
                if n_ref < 4:
                    refs[n_ref] = e.to
                n_ref += 1
        if refused:
            continue

        if at_implicit_h_unknown(a):
            # AN UNRECORDED HYDROGEN COUNT REFUSES THE UNIT ONLY WHERE IT COULD HAVE CHANGED THE
            # ANSWER, which is where the named directions leave room for a hydrogen.  The two lines
            # below are arithmetic on a count and the sentinel is not one, so it may not simply fall
            # through: `spent` would go 15 over and read as a fully consumed valence, and `n_unnamed`
            # would go to 15.
            #
            # The decision is that a stereocentre needs to know its directions, and an unrecorded
            # hydrogen count usually leaves the NUMBER of them unknown: three heavy neighbours and one
            # hydrogen is a centre, three heavy neighbours and nothing is not, and the missing number
            # is precisely which.  Perceiving a unit there would let `set_parity` write a
            # configuration onto an atom that may have no fourth direction to measure it against.
            #
            # AT FOUR NAMED DIRECTIONS THE MISSING NUMBER CANNOT MATTER, and refusing anyway is not
            # caution -- `CFClBrI` with the sentinel on its carbon would perceive no unit, so a stated
            # parity on a FULLY SUBSTITUTED centre would be silently unwritable, which is exactly the
            # shape a query format hands over.  Four named directions leave no room: any
            # value the count could have taken puts the total past four and would be refused below, so
            # the frame is fully determined and the only self-consistent reading of the sentinel here
            # is zero.  This is the same argument as the electron budget's -- what is unknowable is
            # allowed not to be known, provided nothing downstream reads it as a measurement -- and it
            # is why the count is not added rather than defaulted: a zero would be indistinguishable
            # from a derivation, and `implicit_h_of` must still answer None for this atom.
            if n_ref + n_h < 4:
                continue
        else:
            spent += at_implicit_h(a)
            n_unnamed += at_implicit_h(a)
        if any_aromatic:
            spent += 1                    # the one delocalised pi electron; see `_bond_electrons`
        if a.element == 16 and _sulfur_lone_pair(a, spent):
            n_unnamed += 1        # at most one; two pairs are the same direction twice

        if n_ref + n_h + n_unnamed != 4:
            continue
        # `pi_count == 2` exactly, and that is deliberate rather than an oversight: `pi_a` and
        # `pi_b` hold only the FIRST two pi neighbours, so a three-pi record cannot be compared
        # pairwise here without a loop.  A forged S(=N)(=N)(=N)C is therefore admitted.  No real
        # record reaches four directions with three pi bonds -- the valence is spent -- and the
        # automorphism test removes it in any case, so the loop would be code that never runs.
        if pi_count == 2 and _pi_directions_indistinguishable(&atoms[pi_a], &atoms[pi_b]):
            continue

        # Heavy slots (already in `refs[0:n_ref]`, CSR ascending), then hydrogen slots ascending,
        # then the unnamed directions.  The direction count above is exactly four, so `n_ref + n_h`
        # is at most four and neither fill can leave the array.
        for k in range(<uint32_t> n_h):
            refs[n_ref + k] = h_refs[k]
        for k in range(<uint32_t> (n_ref + n_h), 4):
            refs[k] = SU_NO_REF
        # An atom kind's unnamed directions are its TAIL, so the mask is a run of `n_unnamed` bits
        # starting where the named ones stopped, and its popcount is that count.  Ruling F41 is about
        # bond kinds; this spelling is the same information.
        unnamed_mask = ((1 << n_unnamed) - 1) << (n_ref + n_h)
        if _stereo_emit(out, &count, anchored, SU_TETRA, i, refs, 4, <uint8_t> unnamed_mask):
            return -1

    # Pass 2: cumulenes.  Each chain is walked from whichever terminal has the lower slot, so
    # `other <= i` is how the second terminal's visit is dropped rather than a visited array.
    for i in range(n):
        if _chain_degree(ptr, edges, i) != 1:
            continue
        if not _cumulene_walk(atoms, ptr, edges, i, &far, &far_prev, &chain_len):
            continue
        if far <= i:
            continue
        near_mask = 0
        far_mask = 0
        if not _terminal_pair(atoms, ptr, edges, i, _chain_next(ptr, edges, i, SU_NO_REF),
                              &refs[0], &near_mask):
            continue
        if not _terminal_pair(atoms, ptr, edges, far, far_prev, &refs[2], &far_mask):
            continue
        # Two F26-ordered pairs, so the far terminal's two slot bits shift into the record's slots
        # 2 and 3.  This is the whole reason `spare` carries a mask and not a count (ruling F41).
        n_unnamed = near_mask | (far_mask << 2)
        # Odd atom count -> axial, anchored on the centre atom; even -> cis/trans-like, anchored on
        # the lower-indexed terminal.  `chain_len >> 1` is the centre's distance from either
        # terminal exactly because the count is odd.
        if chain_len & 1:
            if _stereo_emit(out, &count, anchored, SU_ALLENE,
                            _chain_nth(ptr, edges, i, chain_len >> 1), refs, 4,
                            <uint8_t> n_unnamed):
                return -1
        # The small-ring cut applies to cis/trans ONLY -- ruling F44, argued in the fragment comment:
        # its "the ring path holds them cis" premise is not a statement about an axial unit, whose
        # terminals are perpendicular and whose two configurations are enantiomers.
        #
        # THIS REFUSAL SUPPRESSES UNITS THAT ARE GENUINELY STEREOGENIC (ruling F63, §4.6).  A small-
        # ring alkene -- cyclohexene is the textbook case -- IS stereogenic in the graph-theoretic
        # sense: no automorphism of the constitution exchanges its two configurations, so it would
        # survive `mark_stereogenic`.  Geometry makes it inaccessible, not symmetry.  The two
        # questions ("is this stereogenic?" and "is this realizable?") have the same answer here, but
        # for different reasons, and conflating them is the defect the spec calls out.
        #
        # ON A KEKULÉ ARENA THIS IS THE ONLY PLACE THE REFUSAL CAN LIVE.  §4.6 specified a separate
        # post-stereogenicity marking stage that would set SU_UNREALIZABLE on survivors.  That stage
        # is provably empty here: it would call `_terminals_share_small_ring` on the units that
        # already passed `_terminals_share_small_ring` at this line -- same predicate, same arena,
        # empty intersection.  Measured: `mark_unrealizable` with a raise-on-mark probe ran against
        # the full test suite (1656 tests), the automorphism corpora (ring4 k=1..6, ring6 k=1..4), and
        # an aromatic corpus (benzene, naphthalene, pyridine, indole, biphenyl, fused polycyclic)
        # and never fired.  The root cause is also structural: on a Kekulé arena an aromatic ring
        # bond is order 2 and indistinguishable from a small-ring alkene, so admitting the unit
        # (to mark it later) would flood `stereo_units()` with one junk candidate per benzene ring.
        # `test_cyclohexene_double_bond_is_excluded_by_small_ring_cut` documents the boundary from
        # the outside; `test_the_small_ring_cut_admits_at_threshold_and_refuses_below` pins the value.
        # When HE_AROMATIC carries a live flag, `_is_chain_bond` already excludes aromatic bonds,
        # so the unit can then be admitted here and a SU_UNREALIZABLE mark (value 2, currently
        # free in the flag nibble) becomes meaningful.  That is the condition under which to revisit.
        elif _terminals_share_small_ring(structure, atoms, i, far):
            continue
        elif _stereo_emit(out, &count, anchored, SU_CIS_TRANS, i, refs, 4, <uint8_t> n_unnamed):
            return -1

    # Pass 3: atropisomers, over bonds taken from their lower-indexed end, so each axis is seen once.
    for i in range(n):
        for k in range(ptr[i], ptr[i + 1]):
            if edges[k].to <= i:
                continue
            if not _is_atropisomer_axis(structure, &edges[k], i, edges[k].to):
                continue
            # ANCHOR CHOICE (ruling F45).  A bond kind may anchor at EITHER end, so a taken lower
            # pivot is not a collision -- it is a reason to anchor at the other one.  This is not
            # hypothetical: an ortho-substituted biaryl of eight-membered rings has a Kekule ring
            # double bond on each pivot that the small-ring cut does not reach (it stops below 8), so
            # the pivot is a cis/trans terminal too, and before this the molecule raised out of
            # `stereo_units()`.  Refusing such a pivot instead would kill EVERY biaryl, because a
            # Kekule aromatic pivot always carries a ring double bond.
            #
            # Relocating costs no meaning.  Parity is stored against `refs`, whose two pairs follow
            # the anchor, so the record says the same thing about the same molecule either way; only
            # which atom is keyed in SEG_PARITY changes.  That makes the choice depend on slot order,
            # i.e. on input order -- exactly as "the lower-indexed terminal" already does for every
            # bond kind -- and there is deliberately no attempt to canonicalise it.
            pivot = i
            other_pivot = edges[k].to
            if _anchor_taken(anchored, pivot):
                pivot = edges[k].to
                other_pivot = i
                if _anchor_taken(anchored, pivot):
                    # Both pivots claimed.  This IS reachable, and rarely: an ortho-substituted
                    # bi(cyclooctatetraenyl) whose two rings each alternate from the pivot gives each
                    # pivot the LOWER terminal of its own ring double bond, so the cumulene pass
                    # anchors at both pivots and the axis has nowhere left to go
                    # (test_a_biaryl_of_eight_rings_can_lose_its_axis_to_both_pivots).  It is a
                    # REFUSAL rather than an assertion because losing one candidate is a thing
                    # perception is allowed to do and raising on a valid molecule is not.
                    #
                    # The principled fix is one more step of the same relocation: a cis/trans unit may
                    # also anchor at EITHER terminal, so the colliding one moves to its own other
                    # terminal and frees a pivot.  Deferred deliberately -- it turns anchor choice into
                    # a cascade (the relocated unit may collide in turn), and the molecules that need
                    # it are eight-ring biaryls, not chemistry anyone draws.
                    continue
            # Both ends already qualified inside the predicate; these two calls are what collects
            # the refs.  Kept separate so the heuristic stays in one function.  The anchor's pair
            # leads, which is what makes the relocation above meaning-preserving.
            _atropisomer_end(atoms, ptr, edges, pivot, other_pivot, &refs[0])
            _atropisomer_end(atoms, ptr, edges, other_pivot, pivot, &refs[2])
            if _stereo_emit(out, &count, anchored, SU_ATROPISOMER, pivot, refs, 4, 0):
                return -1
    return count


# ---------------------------------------------------------------------------------------------
# STEREOGENICITY.  Perception above emits CANDIDATES; everything below decides which of them a
# molecule can actually hold two configurations of.  The predicate is exact and it is one
# sentence: a candidate unit U is stereogenic unless some automorphism of the CONSTITUTION is a
# WITNESS against it -- STABILIZES U SETWISE, acts ODDLY on U's four directions, and is
# STEREO-CONSISTENT at every other unit (it may permute the others, but only in a way their
# stored parities allow).  Such an automorphism carries U's two configurations onto each other,
# so they are the same molecule and U names nothing.
#
# SETWISE, NOT POINTWISE, AND THAT IS RULING F61.  For an ATOM kind the two are the same thing:
# the unit is named on one atom, and nothing but that atom can play its part.  For a BOND kind the
# unit is named on TWO atoms and the anchor is merely whichever of them is keyed in SEG_PARITY -- a
# function of slot order, not of chemistry -- so an automorphism that EXCHANGES the two terminals
# still carries the unit onto itself and is as much a witness as one that fixes them.  Reading the
# predicate as "fixes the anchor" hides every terminal-exchanging witness and therefore OVER-MARKS,
# and it does so on a symmetric macrocycle rather than on something exotic: see
# `test_a_terminal_exchanging_witness_unmarks_a_macrocyclic_double_bond`, where the quarter rotation
# of a 20-ring exchanges an alkene's ends, induces the 4-cycle [2, 3, 1, 0] on the four directions
# and so leaves the bond with one stereoisomer, not two.  This is why phase 4 below runs TWO pinned
# searches for a bond kind: anchor onto anchor, and anchor onto the other terminal.
#
# "Acts oddly" is the whole content.  An automorphism that stabilizes the unit and acts EVENLY is a
# relabelling of the same configuration and says nothing either way, which is why the enumeration
# below filters PERM_PARITY_4 to the twelve odd permutations and searches only for those.
#
# Stereo-consistency is a two-colouring problem, and it is solved as one: each unit is a variable
# over {even, odd}, an automorphism relates variable V to variable sigma(V) by the parity of the
# permutation it induces on V's directions, a stored parity pins a variable to a value, and the
# question "is there an assignment satisfying all of it" is bipartiteness.  Union-find with a
# parity bit answers it in one pass.  A SELF-LOOP with odd parity -- the automorphism fixes some
# other unit's anchor and acts oddly there -- is an immediate contradiction, and that single clause
# is what makes both centres of 1,4-dimethylcyclohexane stereogenic: the arm swap is odd at each
# of them, so it cannot be a witness against either.
# ---------------------------------------------------------------------------------------------


cdef enum:
    SG_UNDECIDED = 0
    SG_YES = 1              # stereogenic
    SG_NO = 2               # refused: some witness exists, or the chemistry rules it out


cdef inline bint _is_group_15_or_16(uint8_t element) noexcept nogil:
    """N/P/As/Sb/Bi and O/S/Se/Te/Po -- the elements whose stereocentres invert through their lone
    pair fast enough that a hydrogen on them names nothing at room temperature.

    Written as a membership test and NOT as `element != 6`, which is the shape it is easy to reach
    for and which is wrong: `O[SiH](CCC)C` is a silicon centre with a hydrogen and it is a genuine,
    resolvable stereocentre, as are the germanium and tin analogues.  Silicon has no lone pair to
    invert through.  Neither has a quaternary ammonium `[N+](C)(C)(C)CC` -- but that has no hydrogen
    either, so the hydrogen half of the test excludes it without needing to know why.
    """
    return (element == 7 or element == 15 or element == 33 or element == 51 or element == 83
            or element == 8 or element == 16 or element == 34 or element == 52 or element == 84)


cdef inline bint _anchor_is_protic(atom_t *atoms, stereo_unit_t *u) noexcept nogil:
    """Does the anchor carry a hydrogen it can invert through?

    Both hydrogen counts, because a drawn hydrogen and an implied one are the same chemistry:
    `at_explicit_h` is derived from the CSR by `derive_scalars`, so this is the same answer as
    walking the adjacency for `element == 1` and it is O(1).  A question about the ANCHOR only --
    a bond kind's far terminal has its own atom record and its hydrogens are not this unit's
    inversion path.

    H_UNKNOWN IS A LIVE INPUT HERE and the unit decides how it reads, which is why this takes the
    unit rather than the anchor id.  A fully substituted anchor is perceived with the sentinel on it
    (see the two H_UNKNOWN cuts above, both narrowed to "only where the count could have changed the
    frame"), so the sentinel does reach this function, and reading it as a raw count would make
    `at_implicit_h(...) != 0` true and deny the centre -- an R4N+ whose count nobody recorded would
    lose its configuration to a comparison against a value that is not a count.

    The unit answers it exactly: NO UNNAMED DIRECTION MEANS NO IMPLICIT HYDROGEN, whatever the nibble
    says, because every one of the four directions is named and an implicit hydrogen is by definition
    unnamed.  Where the unit does have an unnamed slot the sentinel reads as PROTIC, the strict
    direction: "might carry a hydrogen it can invert through" and "does" get the same answer, and that
    answer denies the centre rather than granting it one.  That combination is unreachable today --
    perception admits the sentinel only at a full frame -- and it is spelled out anyway, because the
    safe reading must not depend on the other cut staying exactly as narrow as it is now.
    """
    if not _is_group_15_or_16(atoms[u.anchor].element):
        return False
    if at_explicit_h(&atoms[u.anchor]) != 0:
        return True
    if at_implicit_h_unknown(&atoms[u.anchor]):
        return ((u.spare >> SU_UNNAMED_SHIFT) & SU_UNNAMED_MASK) != 0
    return at_implicit_h(&atoms[u.anchor]) != 0


cdef inline uint32_t _direction_key(stereo_unit_t *u, uint32_t *colour, uint32_t slot) noexcept nogil:
    """A comparable label for direction slot `slot`: two slots with the same key may be
    interchangeable, two with different keys certainly are not.

    Three sources, and the two sentinels have to be told apart (ruling F41): a NAMED direction is
    keyed by its atom's colour, an UNNAMED one (an implicit hydrogen, a lone pair) by a shared
    sentinel, and an EMPTY slot -- no direction there at all -- by a different sentinel.  An empty
    slot can only correspond to an empty slot, and an unnamed direction to an unnamed direction.
    """
    if u.refs[slot] != SU_NO_REF:
        return colour[u.refs[slot]] + 2
    if (u.spare >> SU_UNNAMED_SHIFT) & (1 << slot):
        return 1            # a direction with no atom of its own
    return 0                # no direction at all


cdef inline bint _directions_separated(stereo_unit_t *u, uint32_t *colour) noexcept nogil:
    """Are the slots of every one of this unit's direction lists pairwise distinct under `colour`?

    If they are, NO non-identity permutation of a direction list is available at all, so no ODD one
    is either, and the unit is stereogenic without a search.  Sound for any colouring an
    automorphism must preserve -- refinement classes (the search itself prunes on them) and exact
    orbits both qualify, which is why shortcuts 2 and 3 are one function called twice.

    Per direction LIST rather than across all four, which is stronger and is what the geometry
    says: an atom kind has one list of four, a bond kind two of two, and the anchor pin already
    forbids a cis/trans or atropisomer unit's two pairs from trading places.  Across-all-four would
    lose 2-butene, whose two methyls share an orbit while each terminal's own pair is separated.
    """
    cdef uint32_t width = 4 if u.kind == SU_TETRA else 2
    cdef uint32_t lists = 1 if u.kind == SU_TETRA else 2
    cdef uint32_t base, i, j, l
    for l in range(lists):
        base = l * width
        for i in range(base, base + width):
            for j in range(i + 1, base + width):
                if _direction_key(u, colour, i) == _direction_key(u, colour, j):
                    return False
    return True


cdef inline bint _list_has_two_unnamed(stereo_unit_t *u) noexcept nogil:
    """Shortcut 0a: some direction list holds two directions with no atom of their own.

    Two unnamed directions in one list are indistinguishable BY CONSTRUCTION -- there is nothing to
    tell them apart with, since neither has an atom, an isotope or a substituent -- so swapping them
    is an odd action that is always available and the unit is never stereogenic.  This is a
    statement about the record, not about the graph, so it is tested here and not by the search:
    the automorphism group of the constitution cannot see a lone pair at all.

    Per direction list and per kind: an atom kind's four slots are one list, a bond kind's are two
    pairs, and `CH2=C=CH2` fails on each pair separately while `CH3-CH=C=CH-CH3` fails on neither.
    """
    cdef uint32_t mask = (u.spare >> SU_UNNAMED_SHIFT) & SU_UNNAMED_MASK
    if u.kind == SU_TETRA:
        return _popcount4(mask) >= 2
    return _popcount4(mask & 0b0011) >= 2 or _popcount4(mask & 0b1100) >= 2


cdef inline bint _permutation_expressible(stereo_unit_t *u, uint32_t r, bint axis_pinned,
                                          bint crossed) noexcept nogil:
    """Can the pinned search be ASKED for `PERM_ODD_4[r]`, and does the answer then mean it?

    Only named slots carry a pin, so a permutation is expressible only when the named pins alone
    determine the action on the nameless ones.  Where they do not, the search is not merely wasted --
    it is unsound, because some automorphism satisfies the partial pins and gets credited with a
    parity that is not the one it really has.

    `axis_pinned` says the caller has pinned BOTH atoms a bond kind is named on, and `crossed` says
    which way round: false for anchor-onto-anchor, true for anchor-onto-the-other-terminal.  Three
    things then have to hold.

    PAIRS MAP ONTO PAIRS (bond kinds).  An automorphism that stabilizes a bond kind's axis fixes or
    exchanges the two terminals wholesale, so a permutation splitting one pair across both describes
    no automorphism at all.  `[2, 1, 0, 3]` on penta-2,3-diene is the case: it pins only "the two
    methyls trade places", which the allene's end exchange satisfies -- and the end exchange is the
    pair exchange, which is EVEN, while this row is odd.  Searched, it silently unmarks the allene.

    THE ROW MUST CROSS THE PAIRS EXACTLY WHEN THE AXIS PINS DO (`axis_pinned`).  With the terminals
    pinned onto themselves, a direction of the anchor can only land on a direction of the anchor, so
    a crossing row asks for something the pins forbid and the search would answer "no automorphism"
    for a reason that has nothing to do with the molecule; with the terminals pinned across, only a
    crossing row is coherent.  This is what splits the four rows a bond kind can express into two
    per pin set, and it is only sound because the caller really does pin both terminals -- SU_ALLENE
    is anchored at the chain CENTRE, has no partner to pin, and passes `axis_pinned` false so that
    its one search keeps both halves (its C2 axis performs the pair exchange itself).

    A NAMELESS SLOT MAY ONLY TAKE A NAMELESS SLOT'S PLACE, AND ONLY ONE OF ITS OWN FLAVOUR.  Sigma
    maps atoms to atoms, so it can never carry an implicit hydrogen onto a named neighbour; and a
    PINNED slot -- a lone pair, mask bit clear -- is not a direction that moves at all (ruling F55),
    so it may not trade with a real unnamed direction either.  For an atom kind that leaves nothing
    to move: its four slots are one list, at most one of them is nameless once 0a has had its say,
    and a lone permutation of one element is the identity.  Written as the general test anyway, so
    the atom and bond cases are one rule rather than two branches.

    WHAT THAT CLAUSE DOES *NOT* SAY IS THAT A NAMELESS SLOT KEEPS ITS INDEX.  It does not: the
    wholesale pair exchange carries the nameless slot of pair 0 onto the nameless slot of pair 1, so
    slot 1 legally becomes slot 3 (ruling F55, and the same rule `translate_stereo` obeys).  What is
    invariant is STRUCTURAL -- a nameless slot corresponds to a nameless slot at the same OFFSET
    WITHIN ITS PAIR -- and the flavour test above already enforces exactly that, without a second
    check: ruling F26 orders each pair named-first, then unnamed directions, then empty slots, so two
    pairs admit a flavour-preserving correspondence only when their flavours agree slot by slot, and
    then the offsets agree too.  An explicit offset test would be redundant on every record obeying
    F26 and WRONG on one that does not (a pair written nameless-first has a legal correspondence at a
    different offset), which is why the invariant is argued here rather than asserted below.  It is
    also inert on the rows this loop actually sees: demanding the absolute index instead breaks no
    test and moves no verdict on 822 fuzzed records or either macrocycle family, because when both
    pairs carry one nameless slot the only offset-preserving crossing correspondence is [2, 3, 0, 1],
    which is EVEN and so never enumerated, and the two ODD crossing rows, [2, 3, 1, 0] and
    [3, 2, 0, 1], each send a named slot onto a nameless one and die below regardless.

    All three conditions REFUSE rows, so the effect of getting one too strong is over-marking (a
    witness that exists is not looked for) and never a lost configuration.

    NONE OF THIS IS INSURANCE AGAINST A CASE THAT CANNOT HAPPEN.  Replacing this whole function with
    `return True` fails seven tests, and the reason a nameless slot reaches phase 4 at all is that
    `_directions_separated` is ALL-OR-NOTHING across a unit's lists: it refuses the unit when ANY
    list has a repeated key, so a list carrying a nameless slot rides along whenever some other list
    is unseparated (`CH3-CH=CCl2`), and for an atom kind the nameless slot sits in the same single
    list as the repeat (2,3,4-trichloropentane's middle carbon).  With the gate gone, row
    [0, 1, 3, 2] there writes `pin[refs[2]] = refs[3] == SU_NO_REF`, which IS `CANON_NO_SLOT`, the
    identity satisfies what is left of the pins, and the middle centre disappears.  `_pin_named_row`
    is the second lock on that same door.
    """
    cdef uint32_t i, j
    cdef uint32_t mask = (u.spare >> SU_UNNAMED_SHIFT) & SU_UNNAMED_MASK
    if u.kind != SU_TETRA:
        if (PERM_ODD_4[r][0] < 2) != (PERM_ODD_4[r][1] < 2):
            return False
        if (PERM_ODD_4[r][2] < 2) != (PERM_ODD_4[r][3] < 2):
            return False
        if axis_pinned and (PERM_ODD_4[r][0] >= 2) != (crossed != 0):
            return False
    for i in range(4):
        if u.refs[i] != SU_NO_REF:
            continue
        j = PERM_ODD_4[r][i]
        if u.refs[j] != SU_NO_REF:
            return False
        if ((mask >> i) & 1) != ((mask >> j) & 1):
            return False
    return True


cdef inline bint _pin_named_row(stereo_unit_t *u, uint32_t r, uint32_t *pin) noexcept nogil:
    """Pin every NAMED direction of `u` onto the slot `PERM_ODD_4[r]` sends it to.  False -- and no
    pin written for that slot onward -- when the row would need a named direction to land on a
    nameless one.

    SU_NO_REF AND CANON_NO_SLOT ARE THE SAME VALUE, 0xFFFFFFFF.  So `pin[u.refs[i]] = u.refs[j]`
    with a nameless `j` does not pin that direction to nothing, it silently UNPINS it: the search
    then leaves the slot free, finds an automorphism that ignores it entirely, and credits that
    automorphism with this row's odd parity.  Finding 2 of the automorphism filter review is that
    mechanism seen from the other side -- with the gate above removed, row [0, 1, 3, 2] on
    2,3,4-trichloropentane's middle carbon writes exactly this pin, the identity satisfies what is
    left, and the centre
    vanishes.  `_permutation_expressible` refuses every such row before we get here (the nameless
    slots map among themselves, and a bijection of a finite set that maps a subset into itself maps
    the complement into the complement), so this returning False is the second lock on the one door;
    it is written because the collision is invisible at the assignment and costs one comparison.  And
    what it protects is worse than one relaxed slot: `pin` has already been consumed by
    `_canon_search_order`, which keys on the pinned slot SET, so writing `CANON_NO_SLOT` back in here
    would desynchronize the search order from the pin array rather than merely free a slot.
    Measured, not just argued: deleting the two lines moves no verdict on 822 fuzzed records, the 22
    named cases or either macrocycle family, and breaks no test.
    """
    cdef uint32_t i
    for i in range(4):
        if u.refs[i] == SU_NO_REF:
            continue
        if u.refs[PERM_ODD_4[r][i]] == SU_NO_REF:
            return False
        pin[u.refs[i]] = u.refs[PERM_ODD_4[r][i]]
    return True


cdef inline uint32_t _pf_find(uint32_t *parent, uint8_t *par, uint32_t x,
                              uint8_t *acc) noexcept nogil:
    """Root of `x`, with `acc[0]` receiving the parity accumulated along the way.

    No path compression: the set has one node per stereo unit plus one, the walk is over a forest
    that at most `count` unions ever built, and compression under a parity label has to fold the
    labels as it relinks -- a correct but fiddly loop guarding a cost that no molecule has.
    """
    cdef uint8_t p = 0
    while parent[x] != x:
        p ^= par[x]
        x = parent[x]
    acc[0] = p
    return x


cdef inline bint _pf_union(uint32_t *parent, uint8_t *par, uint32_t x, uint32_t y,
                           uint8_t rel) noexcept nogil:
    """Record `parity(x) XOR parity(y) == rel`; False when that contradicts what is already known.

    False is the whole point of the structure.  Two units related odd-ly to each other and both
    pinned even is a contradiction, and so is the degenerate case x == y with rel == 1 -- the
    self-loop clause, which lands here as `px ^ py == 0 != 1`.
    """
    cdef uint8_t px = 0, py = 0
    cdef uint32_t rx = _pf_find(parent, par, x, &px)
    cdef uint32_t ry = _pf_find(parent, par, y, &py)
    if rx == ry:
        return (px ^ py) == rel
    parent[rx] = ry
    par[rx] = px ^ py ^ rel
    return True


# The twelve ODD permutations of four direction slots, lexicographically -- the complete set of
# actions that can be a witness, and the only ones the enumeration below searches for.  This is
# PERM_PARITY_4 filtered, not a second source of truth: `_permutation_parity_probe` exposes the
# table's parity computation to a test that regenerates the filter in Python and compares.
cdef uint8_t PERM_ODD_4[12][4]
PERM_ODD_4[:] = [[0, 1, 3, 2], [0, 2, 1, 3], [0, 3, 2, 1], [1, 0, 2, 3],
                 [1, 2, 3, 0], [1, 3, 0, 2], [2, 0, 3, 1], [2, 1, 0, 3],
                 [2, 3, 1, 0], [3, 0, 1, 2], [3, 1, 2, 0], [3, 2, 0, 1]]


cdef uint32_t stereo_unit_partner(Structure structure, stereo_unit_t *u) noexcept nogil:
    """The other atom a BOND kind is named on -- the far cis/trans terminal, the other biaryl pivot
    -- or SU_NO_REF for a kind that is named on one atom.

    An atom kind has no partner and neither has SU_ALLENE, whose anchor is the chain's centre and
    whose name is that one atom (spec 3.2); `chiral_bonds` keys on this and `chiral_atoms` gets the
    rest, which is why the allene lands with the atoms.
    """
    cdef atom_t *atoms = structure.atoms()
    cdef uint32_t *ptr = csr_ptr(structure)
    cdef halfedge_t *edges = csr_edges(structure)
    cdef uint32_t far = 0, far_prev = 0, chain_len = 0
    cdef uint32_t k, j
    if u.kind == SU_CIS_TRANS:
        # The chain, not just the neighbour: butatriene's terminals are three bonds apart and the
        # unit is named on the two ENDS.
        if _cumulene_walk(atoms, ptr, edges, u.anchor, &far, &far_prev, &chain_len):
            return far
        return SU_NO_REF
    if u.kind == SU_ATROPISOMER:
        # The pivot adjacent to both of the far end's recorded ring neighbours.  Reading the axis off
        # as "the anchor's one bond that is not in a ring" is what this replaces, and that is wrong
        # on 2-methylbiphenyl, where the methyl is an acyclic bond on the pivot too.
        for k in range(ptr[u.anchor], ptr[u.anchor + 1]):
            j = edges[k].to
            if (csr_find_at(ptr, edges, j, u.refs[2]) is not NULL
                    and csr_find_at(ptr, edges, j, u.refs[3]) is not NULL):
                return j
        return SU_NO_REF
    return SU_NO_REF


cdef inline int _induced_permutation(uint32_t *sigma, stereo_unit_t *v, stereo_unit_t *w,
                                    bint swapped, bint swap_known, uint32_t *perm) noexcept nogil:
    """Fill `perm[0:4]` with the permutation an automorphism induces on `v`'s direction slots when
    it carries `v` onto `w`: `perm[i]` is the slot of `w` that `v`'s slot `i` lands in.  -1 when the
    automorphism does not carry `v` onto `w` after all.

    `swapped` says whether `v`'s anchor mapped onto `w`'s PARTNER rather than onto its anchor, and
    is consulted only when a direction list holds no named direction to locate itself by -- a
    cumulene terminal bearing nothing but an implicit hydrogen.  `swap_known` is false for
    SU_ALLENE, whose anchor is the chain centre and therefore says nothing about which terminal went
    where; a nameless list on an allene is refused rather than guessed at, which can only refuse a
    witness and so only over-mark.

    THE PARITY OF THIS PERMUTATION IS READ WITH `permutation_parity_of` FOR EVERY KIND, including
    the bond kinds that `translate_stereo` deliberately does NOT use it for.  That is
    sound, and the reason is ruling F56: the permutations reachable here preserve the pair PARTITION
    (each of `v`'s lists maps wholly into one of `w`'s, by construction below), and on those
    `permutation_parity_of` equals the pair decomposition `swap0 XOR swap1`, because the pair
    exchange (0 2)(1 3) is a product of two transpositions and contributes nothing to either.  The
    exchange is reachable: an allene's C2 axis performs it.
    """
    cdef uint32_t vmask = (v.spare >> SU_UNNAMED_SHIFT) & SU_UNNAMED_MASK
    cdef uint32_t wmask = (w.spare >> SU_UNNAMED_SHIFT) & SU_UNNAMED_MASK
    cdef uint32_t width = 4 if v.kind == SU_TETRA else 2
    cdef uint32_t lists = 1 if v.kind == SU_TETRA else 2
    cdef uint32_t used = 0
    cdef uint32_t base, wbase, i, j, l, img
    cdef bint found
    if v.kind != w.kind:
        return -1
    for i in range(4):
        perm[i] = 0
    for l in range(lists):
        base = l * width
        # Which of `w`'s lists this one lands in.  The first NAMED slot settles it by lookup, and
        # ruling F26 puts the named slots first, so this is the first slot for almost every record.
        wbase = 4
        for i in range(base, base + width):
            if v.refs[i] == SU_NO_REF:
                continue
            img = sigma[v.refs[i]]
            for j in range(4):
                if w.refs[j] == img:
                    wbase = (j // width) * width
                    break
            break
        if wbase == 4:
            if not swap_known:
                return -1
            wbase = (l ^ (<uint32_t> 1 if swapped else <uint32_t> 0)) * width
        # ONE loop over the list's slots, named and sentinel alike.  A named slot is matched by its
        # image's slot, a sentinel by its FLAVOUR -- an unnamed direction may only correspond to an
        # unnamed direction and an empty slot to an empty slot (ruling F41).  Greedy is unambiguous
        # because 0a has already refused every unit with two unnamed directions in one list, so each
        # list holds at most one of each flavour.
        for i in range(base, base + width):
            found = False
            for j in range(wbase, wbase + width):
                if used & (<uint32_t> 1 << j):
                    continue
                if v.refs[i] != SU_NO_REF:
                    found = w.refs[j] == sigma[v.refs[i]]
                else:
                    found = (w.refs[j] == SU_NO_REF
                             and ((vmask >> i) & 1) == ((wmask >> j) & 1))
                if found:
                    perm[i] = j
                    used |= <uint32_t> 1 << j
                    break
            if not found:
                return -1
    return 0


cdef inline int _stereo_consistent(stereo_unit_t *units, uint32_t count,
                                   uint32_t *anchor_of, uint32_t *partner_of, uint32_t *sigma,
                                   uint32_t self_index, uint32_t *parent, uint8_t *par,
                                   uint8_t *constrains, uint8_t *parities) noexcept nogil:
    """Can `sigma` be a witness as far as the units OTHER than `self_index` are concerned?

    1 yes, 0 no.  Every such unit contributes one constraint -- `parity(V) XOR parity(sigma(V))` is
    the parity of the induced permutation -- and every configured one is pinned to its stored value
    against a ground node.  Union-find with a parity bit decides the whole system at once.

    `self_index` is EXEMPT: the automorphism acting oddly on it is the hypothesis under test, and
    including it would contradict itself by construction (a self-loop with odd parity).  Units with
    `constrains[k] == 0` are exempt too -- 0a refused them, so they have no configuration to be
    consistent about and pinning one would invent a constraint out of a lone pair.

    `par` is the union-find's parity bits, indexed by UNIT.  `parities` is SEG_PARITY, indexed by
    ATOM SLOT.  The two must not be confused: `par` records path parity in the union-find forest
    and is modified throughout, while `parities` is the molecule's stored three-state parity read
    once per unit.
    """
    cdef uint32_t ground = count
    cdef uint32_t k, img, w_index, cand
    cdef uint32_t perm[4]
    cdef stereo_unit_t *v
    cdef stereo_unit_t *w
    cdef uint8_t p
    cdef bint swapped
    for k in range(count + 1):
        parent[k] = k
        par[k] = 0
    for k in range(count):
        if k == self_index or not constrains[k]:
            continue
        v = &units[k]
        p = parities[v.anchor] if parities is not NULL else 0
        if p:
            if not _pf_union(parent, par, k, ground, 1 if p == 2 else 0):
                return 0
        # WHERE THE IMAGE UNIT IS ANCHORED IS NOT sigma OF WHERE THIS ONE IS.  A bond kind anchors
        # at its lower-indexed end, which is a function of input order and not of the constitution,
        # so an automorphism exchanging a cis/trans unit's terminals carries this record onto a
        # record anchored at the OTHER one.  Both places are therefore looked up.  Refusing to look
        # (treating a missing unit at sigma[anchor] as an inconsistency) is safe but over-marks, and
        # it over-marks on a symmetric macrocycle rather than on something exotic.
        w_index = SU_NO_REF
        swapped = False
        img = sigma[v.anchor]
        cand = anchor_of[img]
        if cand != SU_NO_REF and units[cand].kind == v.kind:
            w_index = cand
        else:
            cand = partner_of[img]
            if cand != SU_NO_REF and units[cand].kind == v.kind:
                w_index = cand
                swapped = True
        if w_index == SU_NO_REF:
            return 0
        w = &units[w_index]
        if _induced_permutation(sigma, v, w, swapped, v.kind != SU_ALLENE, perm):
            return 0
        if not _pf_union(parent, par, k, w_index, permutation_parity_of(perm)):
            return 0
    return 1


cdef int mark_stereogenic(Structure structure) except -1:
    """Set SU_STEREOGENIC on every candidate unit that really is one.  Returns 1 when some unit's
    answer had to be taken conservatively because a search was truncated, 0 otherwise.

    Five decisions in increasing cost, and every one of them can only be reached by a unit the
    cheaper ones did not settle:

      0. CHEMISTRY, which the automorphism group cannot see.  Two directions with no atom of their
         own in one list are indistinguishable, and a hydrogen on a group-15/16 anchor inverts.
      1. A TRIVIAL GROUP.  No automorphism at all, so no witness: every survivor is stereogenic,
         and this is where nearly every real molecule is decided.
      2. SEPARATED REFINEMENT CLASSES.  An automorphism preserves the refinement, so a unit whose
         direction lists are pointwise distinct under it admits no permutation, odd or even.
      3. SEPARATED ORBITS.  The same test against the exact partition, which is coarser and so
         decides units that the refinement left open.
      4. THE SEARCH.  For each of the twelve odd permutations, enumerate the automorphisms that
         STABILIZE THE UNIT and realise that permutation, and ask each whether the OTHER units'
         stored parities can live with it.  One that can is a witness and the unit is refused.
         Stabilizing takes two pinned searches for a kind named on two atoms -- terminals fixed and
         terminals exchanged (ruling F61) -- and one for a kind named on a single atom.
      5. TRUNCATION.  A search that ran out of budget proves nothing, so the unit is marked
         stereogenic conservatively: dropping it here is irreversible, since `mark_stereogenic` is
         not called again without rebuilding the table.

    THE BRIEF'S PARITY-SEEDED FIXED-POINT LOOP IS DELIBERATELY NOT HERE, and the report says so at
    length.  Re-refining with `rank[i] * 3 + parity(i)` as the seed asks for `parity(V) == parity(W)`
    where the consistency clause asks for `parity(V) == parity(W) XOR induced`, so it is strictly
    stricter than the predicate and would over-mark; and the predicate is idempotent, so there is
    nothing for a second round to find.  The group is computed ONCE, from the constitution.
    """
    cdef uint32_t count = structure_stereo_unit_count(structure)
    if count == 0:
        return 0
    cdef uint32_t n = structure.header.atom_count
    # HOISTED DELIBERATELY, ABOVE THE FIVE POINTER FETCHES BELOW (ruling F60).
    # `ensure_component_labels` calls `structure_append`, which reallocs the arena, so it is the one
    # thing in this function that can move every pointer into it.  Phase 4 needs the labels to pin
    # the anchor's complement, and `ensure_stereo_units` re-fetches after `mark_stereogenic` for
    # exactly this reason, so appending here is legal -- appending BELOW this line would not be.
    #
    # THE SYMPTOM, SO THAT WHOEVER MOVES IT RECOGNISES THEIR OWN MUTATION: nothing crashes.  The
    # verdicts land in freed memory and the reads that follow return plausible garbage -- record
    # dependent, so most fixtures still pass -- e.g. `validate_stereo()` reporting [2, 4, 6] where [4]
    # is correct.  Moved to just below the `units` fetch it took out 13 tests in
    # `test_stereo_perception.py`; moved a line or two differently, 9.  The count is placement
    # dependent and the shape is not: silent wrong stereo answers, never a segfault.
    ensure_component_labels(structure)
    # Read-only from here on -- neither compute_atoms_order nor mol_automorphisms appends to the
    # arena -- so these pointers stay valid for the whole function (ruling F60).  A call that
    # could append would have to be hoisted above them or the pointers re-fetched after it.
    cdef atom_t *atoms = structure.atoms()
    cdef uint32_t *ptr = csr_ptr(structure)
    cdef halfedge_t *edges = csr_edges(structure)
    cdef uint32_t *comp = structure_component_labels(structure)
    cdef stereo_unit_t *units = structure_stereo_units(structure)
    cdef stereo_unit_t *u
    cdef uint8_t *parities = NULL
    if structure_has(structure, SEG_PARITY):
        parities = structure_parities(structure)
    cdef uint32_t *cls = NULL
    cdef uint32_t *orbits = NULL
    cdef uint32_t *anchor_of = NULL
    cdef uint32_t *partner_of = NULL
    cdef uint32_t *partner_atom = NULL
    cdef uint32_t *pin = NULL
    cdef uint32_t *order = NULL
    cdef uint32_t *anchor = NULL
    cdef uint32_t *sigma = NULL
    cdef uint32_t *cursor = NULL
    cdef uint8_t *taken = NULL
    cdef uint8_t *verdict = NULL
    cdef uint8_t *constrains = NULL
    cdef uint32_t *parent = NULL
    cdef uint8_t *par = NULL
    cdef uint32_t k, i, j, r, s, sets, part, budget, home = 0, undecided = 0, flags = 0
    cdef bint multi = False
    cdef uint64_t spent = 0
    cdef Py_ssize_t classes
    cdef pinned_search_t st
    cdef bint witness, cut, truncated = False
    try:
        cls = <uint32_t *> PyMem_Malloc(<size_t> n * sizeof(uint32_t))
        orbits = <uint32_t *> PyMem_Malloc(<size_t> n * sizeof(uint32_t))
        anchor_of = <uint32_t *> PyMem_Malloc(<size_t> n * sizeof(uint32_t))
        partner_of = <uint32_t *> PyMem_Malloc(<size_t> n * sizeof(uint32_t))
        partner_atom = <uint32_t *> PyMem_Malloc(<size_t> count * sizeof(uint32_t))
        pin = <uint32_t *> PyMem_Malloc(<size_t> n * sizeof(uint32_t))
        order = <uint32_t *> PyMem_Malloc(<size_t> n * sizeof(uint32_t))
        anchor = <uint32_t *> PyMem_Malloc(<size_t> n * sizeof(uint32_t))
        sigma = <uint32_t *> PyMem_Malloc(<size_t> n * sizeof(uint32_t))
        cursor = <uint32_t *> PyMem_Malloc(<size_t> n * sizeof(uint32_t))
        taken = <uint8_t *> PyMem_Malloc(<size_t> n)
        verdict = <uint8_t *> PyMem_Malloc(<size_t> count)
        constrains = <uint8_t *> PyMem_Malloc(<size_t> count)
        parent = <uint32_t *> PyMem_Malloc(<size_t> (count + 1) * sizeof(uint32_t))
        par = <uint8_t *> PyMem_Malloc(<size_t> (count + 1))
        if (cls is NULL or orbits is NULL or anchor_of is NULL or partner_of is NULL
                or partner_atom is NULL or pin is NULL
                or order is NULL or anchor is NULL or sigma is NULL or cursor is NULL
                or taken is NULL or verdict is NULL or constrains is NULL or parent is NULL
                or par is NULL):
            raise MemoryError('stereogenicity scratch allocation failed')

        # --- 0. chemistry, one loop, kind-independent above the two kind-aware helpers ---
        for i in range(n):
            anchor_of[i] = SU_NO_REF
            partner_of[i] = SU_NO_REF
            pin[i] = CANON_NO_SLOT
        for k in range(count):
            u = &units[k]
            anchor_of[u.anchor] = k
            j = stereo_unit_partner(structure, u)
            partner_atom[k] = j     # the SECOND atom this unit is named on, for phase 4's pin set B
            if j != SU_NO_REF:
                partner_of[j] = k
            constrains[k] = 1
            verdict[k] = SG_UNDECIDED
            if _list_has_two_unnamed(u):
                verdict[k] = SG_NO
                constrains[k] = 0   # nothing to be consistent ABOUT: not a parity carrier at all
            elif _anchor_is_protic(atoms, u):
                verdict[k] = SG_NO  # still a parity carrier as far as the GRAPH is concerned, so it
                                    # keeps constraining -- the refusal is chemistry, not symmetry
            else:
                undecided += 1

        # --- 1./2. the refinement, then the group ---
        if undecided:
            with nogil:
                classes = compute_atoms_order(structure, cls, NULL)
            if classes < 0:
                raise MemoryError('stereogenicity refinement failed to allocate')
            for k in range(count):
                if verdict[k] == SG_UNDECIDED and _directions_separated(&units[k], cls):
                    verdict[k] = SG_YES
                    undecided -= 1
        if undecided:
            mol_automorphisms(structure, NULL, orbits, &flags)
            if flags & (CANON_ASYMMETRIC | CANON_BUDGET_EXCEEDED):
                # ASYMMETRIC: no automorphism exists, so no witness can.  BUDGET_EXCEEDED: the
                # orbits may be finer than the truth, so neither shortcut 3 nor a search pruning on
                # them may conclude anything -- both land on the same conservative answer.
                if flags & CANON_BUDGET_EXCEEDED:
                    truncated = True
                for k in range(count):
                    if verdict[k] == SG_UNDECIDED:
                        verdict[k] = SG_YES
                undecided = 0
            else:
                for k in range(count):
                    if verdict[k] == SG_UNDECIDED and _directions_separated(&units[k], orbits):
                        verdict[k] = SG_YES
                        undecided -= 1

        # --- 4./5. the pinned odd-stabilizer enumeration ---
        #
        # RESTRICTED TO THE ANCHOR'S OWN CONNECTED COMPONENT, and the verdicts are IDENTICAL -- not
        # merely conservative -- in both directions.  Both halves of that are needed, because a
        # restriction that only failed to LOSE witnesses could still gain one.
        #
        # NO WITNESS IS LOST.  Let sigma be a witness for unit U: an automorphism realising an odd
        # permutation on U's four directions, stabilizing U setwise, and satisfying
        # `_stereo_consistent`.  Components are blocks of the automorphism group, so sigma maps U's
        # component C onto C.  Define sigma' = sigma on C and the identity everywhere else.  A
        # disjoint union's automorphisms compose componentwise, so sigma' is an automorphism of the
        # whole record; it induces the same odd permutation on U, because every one of U's directions
        # lies in C (see the next paragraph); and it is stereo-consistent, because every unit inside C
        # is checked exactly as it was, while every unit V outside C has sigma'(V) = V with the
        # IDENTITY induced permutation, so V's clause reads `parity(V) == parity(V) XOR 0`.  sigma' is
        # therefore a witness, and it is one the restricted search can find.
        #
        # THAT LAST STEP IS "UNCHANGED", NOT "TRIVIALLY TRUE".  `_induced_permutation` can refuse a unit
        # under the identity itself -- an SU_ALLENE whose four slots are all SU_NO_REF would, since
        # `wbase == 4` with `swap_known` false -- and such a unit's clause fails whether it sits inside C
        # or outside it.  Perception emits no such unit (13 allene / cumulene / ketene / ketenimine /
        # CO2 / CS2 shapes probed, none), and the refusal is sigma-INDEPENDENT, so it cannot separate the
        # restricted search from the unrestricted one.  What the argument needs is that V's clause reads
        # the same before and after, which it does.
        #
        # NO UNIT STRADDLES TWO COMPONENTS, which is what makes "U's component" well defined.  A
        # direction is a CSR neighbour of one of the atoms the unit is named on -- of the anchor for
        # SU_TETRA, of the two chain terminals for SU_ALLENE/SU_CIS_TRANS, of the two pivots for
        # SU_ATROPISOMER (`_perceive_stereo_units`, all three passes) -- and the two named atoms of a
        # bond kind are joined by the chain or by the axis bond, so they share a component too.  Every
        # atom this loop pins to itself for U is therefore in C, and every atom it pins as complement
        # is outside it; the two sets do not meet, which is why the reset below can be two loops.
        #
        # NO WITNESS IS GAINED.  A sigma' the restricted search returns is a genuine automorphism of
        # the whole record satisfying the whole predicate -- the pins only ADD constraints -- so a unit
        # it refuses is genuinely not stereogenic.
        #
        # WHAT IT BUYS.  A witness stabilizes its unit setwise, hence stabilizes C, so every partial
        # map that sends C somewhere else is pure budget burn: repeating one symmetric fragment as k
        # separate components multiplies the group by k! without changing any single unit's answer,
        # which is how five copies of 1,3,5,7-tetramethylcyclooctane -- 60 atoms -- exhaust the budget
        # without these pins and get flagged instead of decided.  Salts, solvates and one side of a
        # reaction are that shape.  Measured on k disjoint 1,3,5,7-tetramethylcyclooctanes, without
        # these pins -> with them: k=5 (60 atoms) 28.1 ms -> 0.0 ms, k=6 (72 atoms) 37.6 ms -> 0.1 ms,
        # k=160 (1,920 atoms) 1754 ms -> 34 ms, and `stereo_truncated` false on all three where it is
        # true without them.
        # CONNECTED records are the control and they must not move: 1,800 / 3,600 / 5,760 atoms ran at
        # 1.045 / 1.084 / 1.050 of the unrestricted build before the `multi` guard below, and at
        # 1.004 / 1.021 / 1.005 with it.
        #
        # THE FLAG IS NOT GONE FROM DISJOINT COPIES -- THE THRESHOLD MOVED, from ~60 atoms to ~5,760,
        # and the ordering against connected records INVERTED.  480 disjoint copies (5,760 atoms, 1,920
        # units) report `stereo_truncated is True` again, while a CONNECTED 5,760-atom methylated
        # macrocycle decides in 130 ms and a 7,680-atom one in 229 ms.  The mechanism is this
        # restriction's own cost: every unit's search walks a pinned prefix of length ~n, so `spent`
        # accrues Theta(n) per unit and the whole-call budget goes as Theta(n * unit_count), i.e. ~n^2
        # against a fixed CANON_MAX_NODES_CALL.  Measured on k disjoint copies, k / atoms / ms /
        # truncated: 20/240/0.5/F, 40/480/1.7/F, 80/960/6.6/F, 160/1920/34/F, 320/3840/178/F,
        # 480/5760/462/TRUE, 640/7680/894/TRUE.  `marked` STAYS EXACTLY CORRECT across the threshold
        # (4*k at every k, including 480 and 640), so this is a soundness-preserving false alarm and not
        # a wrong answer -- a caller reading the flag on a 6,000-atom salt is told "unproven" about
        # marks that are in fact right.  The per-component unit list named below is the fix; it is
        # deliberately not done here.
        #
        # THE COMPONENT-WISE SKIP IN `_stereo_consistent` WAS MEASURED AND DECLINED.  With the
        # complement pinned to the identity, every unit outside C has an identity induced permutation
        # and its clause is provably trivial, so skipping those units is exactly equivalent -- but it
        # is worth 19% at 480 atoms, 5% at 1,920 and 3% at 3,840, i.e. the gain SHRINKS with size,
        # because what dominates is the O(n) pinned prefix this loop writes and not the O(count) scan.
        # It would also couple the two changes: a later reader who removed the pins would have to
        # remove the skip.  Not worth it for 3%.  What WOULD pay, if this ever needs to be faster, is a
        # per-component unit list, which shrinks the prefix and the scan together.
        if undecided:
            # ONE O(n) SCAN INSTEAD OF TWO PER UNIT.  A single-component record is the overwhelmingly
            # common shape and the loops below write nothing on it, but they still WALK n atoms twice
            # for every undecided unit, which measured 1.5-6% of the whole search on connected records
            # from 600 to 7,680 atoms.  `label_components` knows the component count and
            # `ensure_component_labels` discards it, so this recovers it in one pass.  `comp[0]` is
            # safe: `count > 0` was checked at the top, and a unit implies an atom.
            for i in range(n):
                if comp[i] != comp[0]:
                    multi = True
                    break
            st.n = n
            st.cls = cls
            st.pin = pin
            st.order = order
            st.anchor = anchor
            st.sigma = sigma
            st.cursor = cursor
            st.taken = taken
            with nogil:
                for k in range(count):
                    if verdict[k] != SG_UNDECIDED:
                        continue
                    if spent >= CANON_MAX_NODES_CALL:
                        truncated = True
                        verdict[k] = SG_YES
                        continue
                    u = &units[k]
                    part = partner_atom[k]
                    # BOTH ATOMS THE UNIT IS NAMED ON ARE PINNED, and every NAMED direction.  The
                    # pinned SLOT SET is the same for all twelve permutations and for both pin sets
                    # below -- only the images change -- and `_canon_search_order`'s prefix is a
                    # function of the slot set alone, so it is built once per unit.
                    #
                    # Pinning the PARTNER is what makes "stabilizes the unit setwise" the thing that
                    # is searched for.  It costs no witness: a sigma that maps U's directions onto
                    # U's directions must map U's axis onto itself, since a direction is a neighbour
                    # of one of the two named atoms, so the only two possibilities are the two pin
                    # sets below.  It buys pruning, and it makes the crossing rows of set B mean
                    # what they say.
                    pin[u.anchor] = u.anchor
                    if part != SU_NO_REF:
                        pin[part] = part
                    for i in range(4):
                        if u.refs[i] != SU_NO_REF:
                            pin[u.refs[i]] = u.refs[i]
                    # EVERY ATOM OUTSIDE THE ANCHOR'S COMPONENT IS PINNED TO ITSELF, which is the
                    # restriction argued above.  These pins are constant for the whole unit -- both pin
                    # sets and all twelve rows -- so the pinned SLOT SET is still a function of the unit
                    # alone and `_canon_search_order`'s prefix is still built once per unit.  A pinned
                    # depth draws from a one-element candidate list, so the complement costs one budget
                    # node each, once per search, in exchange for the k! the search does not walk.
                    #
                    # THIS PIN IS WHAT THE SPEEDUP IS ATTRIBUTABLE TO, and EXACTLY TWO TESTS FAIL IF IT
                    # GOES -- measured with this loop and its reset both deleted:
                    # `test_k_identical_components_mark_k_times_one_copy_and_do_not_truncate` and
                    # `test_an_apply_that_re_bases_a_parity_leaves_the_marking_to_the_next_reader`, in
                    # both cases on the `stereo_truncated is False` assertion.  The k-copies timings
                    # return to 28.1 / 37.6 / 1754 ms for k = 5 / 6 / 160, and NO VERDICT MOVES -- the
                    # pins buy time, nothing else.  The verdicts are pinned by other tests:
                    # `test_two_units_in_different_components_are_decided_the_same_in_both_orders` and
                    # `test_a_multi_component_record_answers_component_by_component`, neither of which
                    # notices the pins going (their records are far too small to reach the budget) and
                    # both of which notice the RESET going -- see below.
                    #
                    # THE `if multi` PAIR MUST STAY A PAIR, and `home` is why.  The reset below reads
                    # the `home` this line writes, so the two loops are guarded by the SAME condition
                    # and nothing between them may `continue` past the reset -- every early exit in
                    # this unit's body is a `break` out of the row loops, above.  Insert a `continue`
                    # between the pin and the reset and the complement pins leak into the next
                    # component's search; that is a wrong ANSWER, not a slowdown.
                    if multi:
                        home = comp[u.anchor]
                        for i in range(n):
                            if comp[i] != home:
                                pin[i] = i
                    _canon_search_order(n, ptr, edges, pin, order, anchor, taken)
                    witness = False
                    cut = False
                    # SET A is the identity on the axis, SET B exchanges the two terminals.  A unit
                    # is stereogenic only when NEITHER finds a witness (ruling F61); a kind named on
                    # one atom -- SU_TETRA, and SU_ALLENE, whose anchor is the chain centre -- has
                    # no set B and no partner to pin.
                    #
                    # SET B IS STRUCTURALLY UNCOVERABLE FOR SU_ATROPISOMER, and that is a fact about
                    # the kind rather than a gap in the fixtures -- the same standing as the "acts
                    # evenly" row that `test_an_automorphism_that_acts_evenly_is_not_a_witness`
                    # documents.  `_is_atropisomer_axis` requires a single bond with both ends in a
                    # ring and the axis itself NOT `HE_IN_RING`, so the axis is never a ring bond and
                    # the constitution's automorphism group factors over the two rings joined by it.
                    # A sigma that exchanges the terminals therefore acts as a product of a ring
                    # exchange and per-ring flips, and any odd action it has on the four directions is
                    # already delivered by a single-ring flip that FIXES both terminals -- which is a
                    # set A witness.  Measured: no biaryl changes verdict between `sets = 1` and
                    # `sets = 2`.  Set B is live and load-bearing for SU_CIS_TRANS, which is where the
                    # macrocycle test exercises it.
                    sets = 2 if part != SU_NO_REF else 1
                    for s in range(sets):
                        if s:
                            pin[u.anchor] = part
                            pin[part] = u.anchor
                        for r in range(12):
                            # Not every odd permutation can be asked of a search that pins atoms:
                            # `_permutation_expressible` is where that is decided, and the rows it
                            # refuses are the ones whose pins would be satisfied by an automorphism
                            # of a different parity -- or, now, the ones that contradict this pin
                            # set's axis images outright.
                            if not _permutation_expressible(u, r, part != SU_NO_REF, s != 0):
                                continue
                            if not _pin_named_row(u, r, pin):
                                continue    # a named direction onto a nameless one: see the helper
                            budget = CANON_MAX_NODES_SEARCH
                            mol_find_pinned_begin(&st)
                            while mol_find_pinned_next(ptr, edges, atoms, &st, &budget):
                                # A candidate that fails the consistency test proves NOTHING -- only
                                # an exhausted enumeration is a decision, which is why this loop
                                # resumes the search rather than restarting it.
                                if _stereo_consistent(units, count, anchor_of, partner_of,
                                                      sigma, k, parent, par, constrains, parities):
                                    witness = True
                                    break
                            spent += <uint64_t> CANON_MAX_NODES_SEARCH - budget
                            if st.truncated:
                                truncated = True
                                cut = True
                            if witness or cut:
                                break
                        if witness or cut:
                            break
                    pin[u.anchor] = CANON_NO_SLOT
                    if part != SU_NO_REF:
                        pin[part] = CANON_NO_SLOT
                    for i in range(4):
                        if u.refs[i] != SU_NO_REF:
                            pin[u.refs[i]] = CANON_NO_SLOT
                    # AND THE COMPLEMENT PINS GO WITH THEM.  `pin` is reused across units, and the next
                    # unit's component is a different set, so a leftover complement pin silently
                    # OVER-restricts its search -- the next component ends up pinned to the identity,
                    # the identity is even, no witness is found, and units are wrongly MARKED.  Never a
                    # crash: measured with this loop deleted, a record holding 1,4-dimethylcyclohexane
                    # and methylcyclohexane as two components reports THREE chiral atoms where two is
                    # correct, with `stereo_truncated` false, and 18 of 149 salt-shaped records
                    # over-mark for 24 spurious marks in total.  Bond kinds leak the same way: a
                    # spurious cis/trans mark appears on an ethylidenecyclohexane sitting behind any
                    # component that has a unit.
                    #
                    # `test_two_units_in_different_components_are_decided_the_same_in_both_orders` (5 of
                    # its 8 records disagree) and `test_a_multi_component_record_answers_component_by_component`
                    # (4 of its 8) are the tests that notice.  BOTH NEED A MARKED COMPONENT BEFORE A
                    # REFUSED ONE in at least one case, because the leak's first victim is the second
                    # component processed -- a record whose refusals all come first passes even with
                    # this loop gone.  The converse does NOT hold: two of the surviving cases do carry a
                    # mark before a refusal, because sensitivity depends on the component and not only
                    # on its position.  So those tests carry each case in BOTH orders rather than
                    # trusting a rule about ordering to tell them which half is sensitive.
                    if multi:
                        for i in range(n):
                            if comp[i] != home:
                                pin[i] = CANON_NO_SLOT
                    verdict[k] = SG_NO if witness else SG_YES

        for k in range(count):
            if verdict[k] == SG_YES:
                # READ-MODIFY-WRITE ON THE FLAG NIBBLE, never an assignment.  `spare` also carries
                # perception's unnamed-direction mask in its high nibble, and `spare = SU_STEREOGENIC`
                # would erase it -- silently, and only for the units that ARE stereogenic, which is
                # the subset every later task reads.  This is the one sanctioned exception to
                # `_stereo_emit` being the only writer of the field.
                units[k].spare |= SU_STEREOGENIC
    finally:
        PyMem_Free(cls)
        PyMem_Free(orbits)
        PyMem_Free(anchor_of)
        PyMem_Free(partner_of)
        PyMem_Free(partner_atom)
        PyMem_Free(pin)
        PyMem_Free(order)
        PyMem_Free(anchor)
        PyMem_Free(sigma)
        PyMem_Free(cursor)
        PyMem_Free(taken)
        PyMem_Free(verdict)
        PyMem_Free(constrains)
        PyMem_Free(parent)
        PyMem_Free(par)
    return 1 if truncated else 0


cdef int ensure_stereo_units_unmarked(Structure structure) except -1:
    """Fill SEG_STEREO_UNIT with PERCEPTION ONLY, leaving the marks undecided; idempotent.

    Ruling F70.  The table this leaves behind is pure constitution -- `kind`, `anchor`, `refs`,
    `n_refs` and `spare`'s unnamed-direction nibble -- and its header says so: word `[2]` is 0, which
    is how `ensure_stereo_units` knows the marking pass still owes it a run.  Word `[1]` is 0 too,
    because a truncation answer only exists once a search has been asked for one.

    FOR CALLERS THAT READ CONSTITUTION AND NOTHING ELSE.  Three of them exist.

    The journal apply is the first: `_harvest_parities` reads `kind`/`refs`/`anchor` and
    `rebase_parity` reads `kind`/`refs`/`spare`'s high nibble; neither has any use for
    SU_STEREOGENIC, and calling the marking variant made every apply on a molecule carrying one
    parity bit pay two full budgeted stereogenicity searches -- 3449 ms against 0.3 ms for twenty
    trivial `add_atom` calls on six copies of 1,3,5,7-tetramethylcyclooctane, measured.

    The isomorphism kernel is the second, by ruling F76: a stereo primitive asks whether the target
    STATES the configuration it names, and whether that statement is justified is validate_stereo's
    question, not the kernel's -- so matching reads `kind`/`refs`/`n_refs` and the anchor atom's
    parity byte in SEG_PARITY, and never SU_STEREOGENIC or the marks.  Routing it through the marking variant
    would make every stereo query pay a witness search per target, and would drag ruling F62's
    truncation policy into the answer.

    `canonical_stereo_group_ids` is the third, by ruling F95: it re-bases each configured parity onto
    the frame its own colouring names, which needs `kind`, `refs` and the unnamed nibble and nothing
    about whether the unit is stereogenic -- an unmarked table is not merely sufficient there, it is
    the honest input, since a group id must not move because a symmetry search ran out of budget.

    Any OTHER caller must stay on `ensure_stereo_units`: the opt-out is deliberately the narrow
    direction, because an unmarked variant named at a handful of call sites cannot silently unmark a
    reader, whereas unmarked-by-default would have broken every caller not found by grep.

    Not nogil, and the order of what it does is load-bearing: it reallocates the arena through
    `structure_append` -> `PyMem_Realloc`, so every caller must invoke it *before* taking any
    pointer into the arena buffer, and it must itself finish perceiving before it appends.  That
    is why perception writes into a scratch block rather than into the segment: the record count
    is not known until the walk is over, the segment cannot be sized before that, and a walk that
    held arena pointers across the append would be walking freed memory.

    The scratch is one record per atom, which is exact rather than generous -- see the anchor
    no-collision invariant in the fragment comment, which `_stereo_emit` asserts and which this
    bound is the third consumer of.
    """
    if structure_has(structure, SEG_STEREO_UNIT):
        return 0
    cdef uint32_t n = structure.header.atom_count
    cdef size_t slots = <size_t> n if n else 1
    cdef Py_ssize_t count = 0
    cdef stereo_unit_t *scratch = <stereo_unit_t *> PyMem_Malloc(slots * sizeof(stereo_unit_t))
    cdef uint8_t *anchored = <uint8_t *> PyMem_Malloc(slots)
    if scratch is NULL or anchored is NULL:
        PyMem_Free(scratch)
        PyMem_Free(anchored)
        raise MemoryError('stereo unit scratch allocation failed')
    try:
        memset(anchored, 0, slots)
        with nogil:
            count = _perceive_stereo_units(structure, scratch, anchored)
        if count < 0:
            # The walk runs in `nogil` and cannot raise, so it reports the refusal as a negative
            # count and this is where it becomes an exception.
            raise RuntimeError(SU_ANCHOR_COLLISION_MSG)
        # There is deliberately NO second scan over the finished table (ruling F42).  `_stereo_emit`
        # is the single gate: it is the only emitter, so a post-build scan could not fire on any path
        # perception takes, and an enforcement that cannot be observed is worse than one that is live
        # and covered.  A later kind that fills a record by hand would bypass both, so no kind may;
        # `_anchor_taken` is the one place outside the gate where the invariant is reasoned about.
        structure_append(structure, SEG_STEREO_UNIT,
                         SU_COUNT_HEADER + <size_t> count * sizeof(stereo_unit_t))
        # every pointer taken AFTER the append that may have moved the buffer
        (<uint32_t *> structure.segment(SEG_STEREO_UNIT))[0] = <uint32_t> count
        (<uint32_t *> structure.segment(SEG_STEREO_UNIT))[1] = 0    # truncated: no search has run
        (<uint32_t *> structure.segment(SEG_STEREO_UNIT))[2] = 0    # marked: `ensure_stereo_units`
        (<uint32_t *> structure.segment(SEG_STEREO_UNIT))[3] = 0    # reserved
        if count:
            memcpy(structure_stereo_units(structure), scratch,
                   <size_t> count * sizeof(stereo_unit_t))
    finally:
        PyMem_Free(scratch)
        PyMem_Free(anchored)
    return 0


cdef int ensure_stereo_units(Structure structure) except -1:
    """Fill SEG_STEREO_UNIT and DECIDE its SU_STEREOGENIC marks if that has not happened yet.

    The invariant this carries, and every reader that wants the marks is entitled to it:
    `stereo_units()`'s `stereogenic` key and `MoleculeContainer.stereo_truncated` are never absent
    and never undecided.  There is no window in which a unit's record says "candidate" and means
    "not asked yet" -- it either comes through here or it comes through
    `ensure_stereo_units_unmarked`, whose callers read constitution only: the journal apply
    (ruling F70) and the isomorphism kernel (ruling F76).

    Idempotent on both halves, and the two halves are separately idempotent: the segment is built at
    most once and the marking pass runs at most once, gated on header word `[2]` rather than on the
    segment's presence, because after F70 the segment can exist with the marks still owed.

    DOES NOT RAISE WHEN THE STEREOGENICITY SEARCH RAN OUT OF BUDGET; it returns the conservative
    table and records that fact in the header word (ruling F62).  The direction of the error is what
    makes that sound: a unit is unmarked only by a witness FOUND, so a search that stops early can
    leave a unit MARKED whose witness it never reached -- it can never take away a mark the unit
    deserves.  The table is therefore an over-approximation of "could this be a stereocentre", which
    is the safe side of that question, and the marking is exactly right on every record measured to
    truncate so far.  With the witness search restricted to the anchor's own component, two shapes
    reach the flag: a CONNECTED record with a large group -- cyclo[CH(CH3)CH2]6 with every hydrogen
    explicit, 54 atoms, `marked == 6`, which is
    `test_a_connected_record_can_still_exhaust_the_budget` -- and disjoint copies once the record is
    big enough to spend the whole-call budget on pinned prefixes, ~5,760 atoms and up (480 copies of
    1,3,5,7-tetramethylcyclooctane: `marked == 4 * copies` exactly, the value the decisive single-copy
    answer predicts).

    The word is PER MOLECULE and says nothing about WHICH units were affected, so a caller cannot
    narrow it -- it reads the flag as `MoleculeContainer.stereo_truncated` and decides what an
    unproven mark is worth to it.

    Why this differs from `automorphism_orbits`, which does raise on the same exhaustion: that
    answers an exactness question, where a truncated search returns SOME labelling in place of THE
    labelling and there is no safe degraded answer; this answers a soundness question, where the
    degraded answer is the conservative one.  Exactness there, soundness here, deliberately.

    Not nogil and it reallocates the arena, for the same reason the unmarked half does; ruling F60
    applies unchanged, and the segment pointer below is re-fetched after the build rather than held
    across it.
    """
    ensure_stereo_units_unmarked(structure)
    # RE-FETCHED AFTER the call that may have appended a segment and moved the buffer (ruling F60).
    if (<uint32_t *> structure.segment(SEG_STEREO_UNIT))[2] != 0:
        return 0
    # The marking must run AFTER the append -- it writes into the segment, and perception's scratch
    # is already gone by then.  The truncation answer is RECORDED, not raised: see the docstring,
    # ruling F62, and `test_a_connected_record_can_still_exhaust_the_budget`.
    cdef int truncated = mark_stereogenic(structure)
    # ...and re-fetched again, because `mark_stereogenic` is not required to leave the arena in place.
    (<uint32_t *> structure.segment(SEG_STEREO_UNIT))[1] = 1 if truncated else 0
    (<uint32_t *> structure.segment(SEG_STEREO_UNIT))[2] = 1
    return 0


# ------------------------------------------------------------------------------------------------
# PARITY RE-BASING.  A stored parity is a statement about ONE list of directions in ONE order, and
# both halves of that are properties of the CSR row it was written against.  An edit rebuilds the
# arena, so every stored parity has to be either re-expressed in the new order or dropped -- and
# because `_apply` is the only writer of the arena, doing it there is what makes a stale sign
# structurally impossible instead of a thing every future call site has to remember.
#
# WHAT IS AND IS NOT A RE-BASE.  Slots never move relative to one another: `add_atom` appends,
# deletion compacts monotonically, and `remap` rewrites stable ids in place (ruling F64), so the
# CSR order of any two surviving neighbours is invariant across an apply.  What CHANGES a row is
# a direction entering or leaving it, and the interesting case is the direction that has no atom of
# its own: bond an atom to the anchor and state its implicit hydrogen away in the SAME edit and the
# unit stays alive with the same four directions, one of which is now named -- and named atoms sort
# by slot, so the new one can land in FRONT of an old one and shift it.  That is a permutation of
# the anchor's directions with the unit intact, and it is the whole reason this function exists.
#
# THE RULE, in one sentence: match each old direction to the new position holding the same thing,
# and the answer is the stored bit XOR the parity of that permutation -- the same ARITHMETIC
# `translate_stereo` would do if it were handed the old ref order.  Only the arithmetic is shared:
# the two functions are not the same question, and reading them as one is what produced ruling F72's
# wrong keep.  `translate_stereo` interprets a frame the CALLER ordered, so re-listing the two pairs
# the other way round is a legitimate re-ordering of one geometry and must come out even; this
# function compares two CANONICAL frames on the same anchor, where the anchor's own pair leads in
# both, so a named atom that moves between pairs crossed the double bond rather than being re-listed.
# See `rebase_parity`'s docstring, which states the asymmetry and the guard that enforces it.
#
# WHAT "THE SAME THING" MEANS, and it is POSITION rather than IDENTITY.  A direction with no atom of
# its own -- an implicit hydrogen, a lone pair -- matches whatever occupies its position afterwards,
# and an atom that stops being a direction leaves its position to whatever takes it.  That is not a
# leniency, it is ruling F26's promise: the ref tuple's shape does not change when a hydrogen starts
# or stops being drawn, so neither may the sign, and a rule keyed on the direction SET would clear
# the parity of every molecule whose hydrogens get explicitated -- which is most of them, and both
# directions of it (`test_explicitating_the_hydrogen_leaves_the_sign_alone`,
# `test_implicitating_the_hydrogen_leaves_the_sign_alone`).  ONE substituted direction is the same
# case with a heavier atom in the hydrogen's role and is treated the same way; TWO at once is not a
# frame any more, and the sign dies rather than be guessed at.
#
# THE CORRESPONDENCE IS THREE-WAY TYPED (ruling F69), and the law is `_direction_key`'s own, written
# there for the stereogenicity search and binding here for the same reason: *an empty slot can only
# correspond to an empty slot, and an unnamed direction to an unnamed direction*.  A ref slot is
#
#   NAMED    -- `refs[i] != SU_NO_REF`; its identity is the atom, and it matches by atom;
#   UNNAMED  -- `SU_NO_REF` with bit i SET in `spare >> SU_UNNAMED_SHIFT`: a direction with no atom
#               of its own, an implicit hydrogen or a sulfur lone pair;
#   EMPTY    -- `SU_NO_REF` with bit i CLEAR: not a direction at all.  `CH3CH=NOH` is the record
#               that exists -- refs `(CH3, None, O, None)`, mask `0b0010`, where slot 1 is the
#               carbon's implicit hydrogen and slot 3 is the oxime nitrogen's lone pair, which
#               `_terminal_pair` deliberately does not count as a direction.
#
# THIS LAW IS DEFENCE IN DEPTH TODAY, and saying so is the honest way to carry it.  No PERCEIVED
# direction list holds both an UNNAMED and an EMPTY slot, and that is structural rather than a survey: a
# SU_TETRA record's four slots are all real directions -- neighbour, implicit hydrogen, or the one sulfur
# lone pair -- so EMPTY cannot occur there at all, and EMPTY otherwise occurs only as the second slot of
# a bond-kind pair, where `_terminal_pair` refuses a pair whose BOTH slots are nameless.  With at most
# one nameless slot per list the pairing is forced whichever way the mask reads, so on reachable frames
# the mask confirms rather than decides.  It is still carried and still enforced, for two reasons: it is
# what keeps a FORGED frame honest (`test_rebase_obeys_the_empty_versus_unnamed_law`), and it becomes
# decisive the day a kind with two nameless slots in one list is perceived -- which is a kind away, not a
# rewrite away.
#
# So the permutation is built in FOUR PASSES -- named by identity, then empty to empty, then unnamed
# to unnamed, all three in ascending order, and only then whatever is left on each side pairs up
# ascending.  THE LEFTOVERS ARE THE SUBSTITUTIONS, and they are what no arithmetic can guess through:
# an old position with nothing of its own to correspond to, faced with a new position that holds
# something the old frame never mentioned.
#
# THE LEFTOVER BUDGET IS ONE OLD AND ONE NEW POSITION PER DIRECTION LIST, NOT PER UNIT.  A direction
# list is all four slots for SU_TETRA and each of the two pairs for the bond kinds -- the same
# per-list rather than per-unit accounting `_directions_separated` and `_list_has_two_unnamed` make,
# for the same geometric reason.  The granularity is load-bearing in both directions:
#
#   * explicitating BOTH implicit hydrogens of a cis/trans double bond is one substitution in EACH
#     pair.  A unit-wide budget of one would drop that sign, and ruling F26 forbids it.
#   * TWO arrivals into ONE terminal -- `(C5, unnamed)` becoming `(Cl, Br)` -- is a guess, and a
#     unit-wide count of vanished NAMED directions sees only one and keeps it.  It also keeps a
#     tetrahedral `(F, Cl, Br, unnamed)` turning into `(Cl, Br, I, At)`, where the two arrivals may
#     take the two free positions either way round and the two ways disagree.
#   * one arrival replacing one departure anywhere is re-based positionally, exactly as before.
#
# `RB_GONE` -- an old NAMED direction whose atom the edit deleted -- is simply a leftover on the old
# side and counts against that same budget.
#
# AND EACH LIST CORRESPONDS TO ITSELF (ruling F75).  For a bond kind, pair 0 is the anchor's own end in
# both frames and the anchor is the same atom, so an old direction found in the OTHER pair crossed the
# double bond rather than being re-listed: the frame was rebuilt, not re-ordered, and the sign dies.
# `translate_stereo` accepts that same exchange as even and is right to, because the frame IT reads is a
# caller's ordering; the divergence is a decision, and `rebase_parity`'s docstring is where it is argued.
cdef int RB_DROP = -1       # ...and the answer when the old frame is not in the new table at all
cdef uint32_t RB_GONE = 0xFFFFFFFEu    # an old direction whose ATOM the edit deleted; see below


cdef int rebase_parity(Structure structure, uint32_t anchor_slot, uint8_t old_kind,
                       uint8_t old_parity, uint32_t *old_refs,
                       uint8_t old_unnamed_mask) noexcept nogil:
    """The stored parity re-expressed against the CURRENT unit at `anchor_slot`, or `RB_DROP`.

    `old_refs[0:4]` is the pre-edit ref list translated into CURRENT slots, with `SU_NO_REF` for
    every pre-edit slot that named no atom and `RB_GONE` for one whose atom the edit deleted.
    `old_unnamed_mask` is that same pre-edit record's `spare >> SU_UNNAMED_SHIFT`, and it is not
    optional: it is the only thing that tells an UNNAMED old slot (a real direction with no atom of its
    own) from an EMPTY one (no direction there at all).  On a PERCEIVED frame that distinction never
    changes the answer -- see the fragment comment for the structural reason -- and it is required
    anyway, because the law it enforces is the one thing standing between a forged frame and a computed
    sign, and because it becomes decisive the day a kind with two nameless slots in one direction list
    is perceived.  `old_parity` is the pre-edit stored value, 0
    for none, 1 even, 2 odd; `old_kind` is the pre-edit kind.  The table must already be built --
    `ensure_stereo_units_unmarked` is enough, since only constitution is read here; take
    `structure`'s pointers after that call, not before (ruling F60).

    Returns 0 for `old_parity == 0` -- nothing configured, nothing to re-base, the same guard
    `translate_parity` opens with -- 1 or 2 for a sign that still means something, and `RB_DROP`
    when it does not:

      * the anchor anchors no unit -- the frame is gone (a fourth heavy neighbour on a tetrahedral
        centre, a deleted bond, a chain that grew and moved its anchor);
      * the kind changed under the same anchor -- a sign about an axis is not a sign about a
        centre.  Structurally redundant today, since a kind change also changes the directions and
        the leftover budget would drop it anyway, and kept because the pair arithmetic below is only
        meaningful for the kind that was measured;
      * MORE THAN ONE leftover position on either side of ONE DIRECTION LIST.  A leftover is a
        position whose content the typed correspondence could not match: an old direction that is
        not a direction of the new unit, and a new direction the old frame never mentioned.  One
        such pair is a substitution and is re-based positionally, which is ruling F26's promise about
        hydrogens with a heavier atom in the hydrogen's place; two in one list leaves the two
        arrivals interchangeable and the two ways round disagree, so the sign dies rather than be
        guessed at.  PER LIST, not per unit -- see the fragment comment for why both halves of that
        are needed.
      * a bond kind ANY of whose named directions is now in the OTHER pair, INCLUDING the wholesale
        exchange (ruling F75).  Pair 0 is the anchor's own end in both frames -- ruling F26,
        `_stereo_emit`, and the cumulene and atropisomer packing all fix it -- and this function looks
        the new unit up by the SAME anchor slot the old frame was harvested from, so a named atom that
        was in pair 0 and is now in pair 1 was not re-listed: it is bonded to the other end of the
        double bond, all four bonds having been broken and remade.  The old sign constrains nothing
        about that, so it dies.

        AND THIS IS WHERE `translate_stereo` DIVERGES, deliberately: it accepts the wholesale exchange
        as even (`test_dichlorobut_2_ene_pair_exchange_does_not_flip`) and must, because ITS frame
        comes from a CALLER, who may legitimately list the two pairs either way round -- there the
        exchange is a re-ordering of one geometry.  Here both frames are PERCEIVED, canonical and
        keyed to one anchor, so the pair order is not anybody's to choose and an exchange can only be
        a migration.  The two functions ask different questions and the divergence is the answer to
        both, not an inconsistency to tidy up (ruling F75 reverses F72, which held the opposite; the
        cost of the F72 reading was a reachable wrong keep on `C(Cl)(Br)=C(F)(I)`).

        THIS GUARD'S PREMISE IS THE ANCHOR.  An AXIS-KEYED harvest -- the follow-up that would let a
        Kekule shift or `thiele` keep a configuration whose anchor moved -- would compare frames
        anchored at DIFFERENT atoms, and there old pair 0 legitimately can correspond to new pair 1.
        Whoever writes that has to revisit this, not extend it.
    """
    cdef stereo_unit_t *u
    cdef uint32_t perm[4]
    cdef bint matched[4]        # old position i has been given a new position
    cdef uint8_t new_unnamed_mask
    cdef uint32_t old_kinds[4]  # 0 EMPTY, 1 UNNAMED, 2 NAMED-or-GONE
    cdef uint32_t new_kinds[4]
    cdef uint32_t i, j, used = 0
    cdef uint32_t width, lists, l, base, old_left, new_left
    cdef uint32_t pp
    if old_parity == 0:
        return 0
    u = stereo_unit_of(structure, anchor_slot)
    if u is NULL:
        return RB_DROP
    if u.kind != old_kind:
        return RB_DROP
    new_unnamed_mask = <uint8_t> (u.spare >> SU_UNNAMED_SHIFT)
    perm[0] = 0; perm[1] = 0; perm[2] = 0; perm[3] = 0
    # THE THREE-WAY CLASSIFICATION of both sides, ruling F69 and `_direction_key`'s law.  `RB_GONE`
    # is typed with the NAMED slots it came from: it was a named direction, and losing its atom makes
    # it a leftover rather than an unnamed direction that could pair with one.
    for i in range(4):
        matched[i] = False
        if old_refs[i] != SU_NO_REF:
            old_kinds[i] = 2
        elif (old_unnamed_mask >> i) & 1:
            old_kinds[i] = 1
        else:
            old_kinds[i] = 0
        if u.refs[i] != SU_NO_REF:
            new_kinds[i] = 2
        elif (new_unnamed_mask >> i) & 1:
            new_kinds[i] = 1
        else:
            new_kinds[i] = 0
    # PASS 1 -- NAMED to NAMED by identity, and named FIRST so that the typed passes below take what
    # is left rather than race them.  Deliberately across ALL FOUR new slots rather than within the old
    # slot's own list, and for a bond kind that is what DETECTS a crossing instead of laundering it: an
    # atom that moved to the other pair is found there and refused just below.  Searching only its own
    # list would leave it unmatched, hand its position to whatever arrived, and re-base a sign across a
    # bond that was broken and remade.
    for i in range(4):
        if old_refs[i] == SU_NO_REF or old_refs[i] == RB_GONE:
            continue
        for j in range(4):
            if u.refs[j] == old_refs[i] and not (used & (1u << j)):
                perm[i] = j
                used |= 1u << j
                matched[i] = True
                break
    width = 4 if old_kind == SU_TETRA else 2
    lists = 1 if old_kind == SU_TETRA else 2
    # THE LIST CORRESPONDENCE IS THE IDENTITY, and it is PINNED rather than read off (ruling F75).
    # Every pass after this one matches an old slot only against a new slot of its OWN list, which is
    # what keeps a cis/trans terminal's implicit hydrogen from being handed the OTHER terminal's
    # implicit hydrogen: both are unnamed, both are unmatched, and an ascending sweep across all four
    # slots would pair them and then produce a permutation the pair guard refuses -- measured, on
    # `test_substituting_the_implicit_hydrogen_rebases_a_cis_trans_sign`.
    #
    # For SU_TETRA there is one list and it corresponds to itself.  For a bond kind list l corresponds
    # to list l and to nothing else, because pair 0 is the ANCHOR'S OWN end of the frame in both the
    # old record and the new one (ruling F26 and `_stereo_emit`; cumulenes walk from the lower-slot
    # terminal and atropisomers pack the anchor's pair first), the anchor is the same atom, and slots
    # never move relative to one another across an apply -- the invariant this fragment opens with.  So
    # a named direction that turns up in the OTHER pair did not get re-listed, it MIGRATED across the
    # double bond, and that is refused here rather than re-based.  `translate_stereo` DOES accept that
    # exchange, and must, because its frame is the CALLER'S ordering rather than a perceived one; the
    # divergence between the two is the decision, argued in the docstring.
    if lists == 2:
        for i in range(4):
            if matched[i] and (perm[i] >> 1) != (i >> 1):
                return RB_DROP
    # PASSES 2 AND 3 -- EMPTY to EMPTY, then UNNAMED to UNNAMED, each in ascending order and each
    # within its own list.  An empty slot can only correspond to an empty slot and an unnamed direction
    # to an unnamed direction, so the oxime frame `(CH3, unnamed, O, empty)` hands its empty slot to the
    # new empty slot and its unnamed slot to whatever unnamed slot survives in its own pair -- never to
    # the hydrogen that got drawn.
    for l in range(2):          # l == 0 is the EMPTY class, l == 1 the UNNAMED one
        for i in range(4):
            if matched[i] or old_kinds[i] != l:
                continue
            base = 0 if lists == 1 else (i >> 1) * 2u
            for j in range(base, base + width):
                if new_kinds[j] == old_kinds[i] and not (used & (1u << j)):
                    perm[i] = j
                    used |= 1u << j
                    matched[i] = True
                    break
    # THE LEFTOVER BUDGET, counted PER DIRECTION LIST before pass 4 assigns anything: one old and one
    # new position at most.  Every pass above matches one old slot against a new slot of the
    # corresponding list, so the two counts agree list by list; both are checked anyway, because the
    # cost is four comparisons and the alternative is trusting that agreement.
    for l in range(lists):
        base = l * width
        old_left = 0
        new_left = 0
        for i in range(base, base + width):
            if not matched[i]:
                old_left += 1
            if not (used & (1u << i)):
                new_left += 1
        if old_left > 1 or new_left > 1:
            return RB_DROP
    # PASS 4 -- the leftovers pair up in ascending order within the corresponding list.  THESE ARE THE
    # SUBSTITUTIONS: an old position with nothing of its own to correspond to takes the position of
    # whatever replaced it, which is the same correspondence `translate_stereo`'s phase 3 makes
    # between a caller's Nones and the stored unnamed slots, so the two directions of the arithmetic
    # agree by construction.
    for i in range(4):
        if not matched[i]:
            base = 0 if lists == 1 else (i >> 1) * 2u
            for j in range(base, base + width):
                if not (used & (1u << j)):
                    perm[i] = j
                    used |= 1u << j
                    matched[i] = True
                    break
            if not matched[i]:
                # its own list is full while another has room.  The budget above counts both sides of
                # every list, so I have no case that reaches this; it returns rather than leaving
                # `perm[i]` at its initialised 0, because a refusal is the only safe thing to do with a
                # permutation that was never completed.
                return RB_DROP
    if old_kind == SU_TETRA:
        pp = permutation_parity_of(perm)
    else:
        # BOND KINDS (rulings F56 and F75): each within-pair transposition is odd, so the answer is the
        # XOR of the two, and the pairs themselves cannot have traded places -- the crossing refusal
        # above has already sent any frame that looks like an exchange to RB_DROP.  This guard is
        # defence in depth over that one: two conditions on one invariant, so a future edit to either
        # cannot quietly compute a sign for a permutation that leaves its pair.
        if perm[0] > 1u or perm[1] > 1u or perm[2] < 2u or perm[3] < 2u:
            return RB_DROP
        pp = (1u if perm[0] != 0u else 0u) ^ (1u if perm[2] != 2u else 0u)
    return <int> (((old_parity - 1u) ^ pp) + 1u)


# ------------------------------------------------------------------------------------------------
# DEFERRED VALIDATION.  Which stated parities the finished molecule cannot justify.
#
# ASKED ON DEMAND AND NOT INSIDE `add_atom_stereo`.  Judging one configuration at a time WHILE the
# molecule is being built stores two of 2,3,4-trichloropentane's three configurations and raises
# nothing: the middle centre is stereogenic only because its two arms are enantiomeric, and when the
# middle sign is offered the second arm does not exist yet.  Ordering is the answer, not a better
# predicate.  A parity is accepted into SEG_PARITY UNCONDITIONALLY at journal-append time
# (`set_parity` range-checks 0..2 and nothing else), and the judgement happens once, on demand, here.
#
# WHAT COUNTS AS UNJUSTIFIED, and there are exactly two shapes of it:
#   * the atom anchors no unit at all -- a lone carbon mid-edit, a CH2Cl2 carbon with two directions;
#   * the atom anchors a unit that `mark_stereogenic` did not mark, so flipping the sign would give
#     the same molecule back and the sign is not a fact about it.
#
# TRUNCATION CANNOT CAUSE A WRONG ANSWER HERE, and the reader's first instinct runs the other way, so
# it is worth stating: `mark_stereogenic`'s decision 5 marks every unit of a record whose search ran
# out of budget, so a parity sitting on a truncated candidate is on a MARKED unit and is justified.
# The over-approximation is on the safe side of this question too -- a conservative mark keeps input,
# it never discards it (ruling F62; `test_a_truncated_record_still_validates_clean`).
cdef int collect_stereo_rejections(Structure structure,
                                   uint32_t **slots_out, uint32_t *count_out) except -1:
    """The atom slots carrying a parity that no stereogenic unit backs, ascending by SLOT.

    Writes a `PyMem_Malloc`ed array of slots to `slots_out[0]` and its length to `count_out[0]`; the
    CALLER FREES the array.  Nothing is allocated when there is nothing to report -- `slots_out[0]`
    stays NULL and `count_out[0]` stays 0 -- which is what keeps the common molecule's validation a
    pure read, and it is why this counts before it allocates rather than growing a buffer.

    ASCENDING BY SLOT, which is not ascending by stable id: `remap` rewrites ids in place and slot
    order is fixed at build time, so the two orders agree only until somebody remaps.  The container
    sorts the stable ids it maps these to (`validate_stereo`), and that is where the ORDER of the
    reported list is decided.

    Calls the MARKED `ensure_stereo_units` (ruling F70): the whole question is whether a stated
    parity is justified, and "justified" is the SU_STEREOGENIC mark.  That call REALLOCATES THE ARENA
    (ruling F60), and `structure_parity_at` inlines the absent-segment guard, so no pointer into the
    arena is needed here.

    This function only REPORTS.  It does not clear, and it must not: the bits live in an arena that
    `copy()` shares outright, so clearing one is a clone-and-rebind that only the container can do
    (ruling F65).
    """
    cdef uint32_t n
    cdef uint32_t i, k, count = 0
    cdef stereo_unit_t *u
    cdef uint32_t *out
    slots_out[0] = NULL
    count_out[0] = 0
    # BEFORE any pointer into the arena: this builds the table and decides its marks, and both
    # halves can move the buffer.
    ensure_stereo_units(structure)
    n = structure.header.atom_count
    if n == 0:
        return 0
    for i in range(n):
        if structure_parity_at(structure, i):
            u = stereo_unit_of(structure, i)
            if u is NULL or not (u.spare & SU_STEREOGENIC):
                count += 1
    if count == 0:
        return 0
    out = <uint32_t *> PyMem_Malloc(<size_t> count * sizeof(uint32_t))
    if out is NULL:
        raise MemoryError('stereo rejection list allocation failed')
    k = 0
    for i in range(n):
        if structure_parity_at(structure, i):
            u = stereo_unit_of(structure, i)
            if u is NULL or not (u.spare & SU_STEREOGENIC):
                out[k] = i
                k += 1
    slots_out[0] = out
    count_out[0] = count
    return 0


def _permutation_parity_probe(p):
    """`permutation_parity_of` on a Python 4-tuple, so a test can regenerate PERM_ODD_4's filter."""
    cdef uint32_t perm[4]
    cdef int i
    for i in range(4):
        perm[i] = <uint32_t> p[i]
    return permutation_parity_of(perm)


def _odd_permutation_table():
    """PERM_ODD_4 as a tuple of tuples, for the test that proves it is exactly the odd half."""
    cdef int r, i
    cdef list rows = []
    cdef list row
    for r in range(12):
        row = []
        for i in range(4):
            row.append(PERM_ODD_4[r][i])
        rows.append(tuple(row))
    return tuple(rows)


def _stereo_unit_record_size():
    return sizeof(stereo_unit_t)


def _stereo_anchor_collision_probe():
    """Emit twice on one anchor through `_stereo_emit`, so a test can watch the gate refuse.

    `_stereo_emit` is the SINGLE enforcement of the anchor no-collision invariant (ruling F42), and no
    molecule reaches it: every kind that could collide either refuses locally or relocates its anchor,
    which is the point.  So the live gate is only observable through a probe, and this one exercises
    the gate itself rather than a copy of its logic.  Raises the same RuntimeError
    `ensure_stereo_units` does, and AssertionError if the refusal failed to fire.
    """
    cdef stereo_unit_t scratch[2]
    cdef uint8_t anchored = 0
    cdef uint32_t refs[4]
    cdef Py_ssize_t count = 0
    cdef int k
    for k in range(4):
        refs[k] = SU_NO_REF
    if _stereo_emit(scratch, &count, &anchored, SU_TETRA, 0, refs, 4, 0x0F):
        raise AssertionError('the first unit on a free anchor was refused')
    if _stereo_emit(scratch, &count, &anchored, SU_CIS_TRANS, 0, refs, 2, 0) == 0:
        raise AssertionError('a second unit on a taken anchor was accepted')
    raise RuntimeError(SU_ANCHOR_COLLISION_MSG)


# --- canonical stereo group ids (ruling F79) ----------------------------------------------------
#
# The stored group id is OPAQUE: it comes from the input file, it is 1..63 because that is what fits
# beside the kind in one byte, and MDL's numbering carries no meaning beyond "these atoms belong
# together".  Two molecules that are the same molecule with the same groups can therefore differ in
# every stored id, which makes the stored ids unusable in anything derived -- a signature, a hash, a
# comparison.  What is well defined is the PARTITION: which atoms share a group, and of which kind.
#
# So the canonical id is a VIEW and not a renumbering: nothing is written back to the segment.  A
# renumbering would have to clone the arena (ruling F65 forbids writing through a shared one), bump
# the generation, and detach every copy() that shared it -- all to store a number that is a function
# of what is already stored.  RULING F79: compute it, return it, leave the segment alone.
#
# The order is ruling F80's: groups of one kind are ordered by where their members sit in the
# molecule's CANONICAL atom order, which is what makes the answer independent of both the stored ids
# and the order the atoms were added in.


# RULING F95 -- THE STORED PARITY BYTE MAY NOT REACH A SEED, BECAUSE IT IS FRAME-RELATIVE.
#
# A parity says "these four directions, IN THIS ORDER, turn this way", and the order is ruling F26's:
# the anchor's heavy neighbours in CSR-SLOT ASCENDING order, then its explicit hydrogens, then the
# directions with no atom of their own.  Slot order is the order the atoms were added in, so ONE
# MOLECULE WRITTEN TWICE HAS TWO BYTE PATTERNS -- which is the entire reason `translate_stereo`
# exists.  A seed that reads the raw byte therefore colours the ENCODING and not the molecule.
#
# Measured, on 1,2,3,4-tetrachlorocyclobutane with OR groups {C0,C1} and {C2,C3}: creation order
# (0,1,2,3) stores (1,1,1,2) and creation order (1,2,0,3) stores (1,2,1,1), and those are the SAME
# molecule -- every centre reads the same parity in the chemical frame (next ring atom, previous ring
# atom, its chlorine, its hydrogen), which is `translate_stereo`'s own definition of a frame change.
# With the byte in the seed the two encodings exchanged which group got canonical id 1 while BOTH
# reported every id PINNED.  No ambiguity class can excuse that: both readings claim to be pinned and
# they disagree.  AT HEAD THAT SWEEP MOVES NOTHING: over the C4/C5/C6 rings, every ring-frame parity
# pattern (2^n) x every set partition of the centres (Bell(n)), each re-encoded 24 ways for C4 and C5
# and 12 for C6 and each order compared against the first, the view moved on 0 of 5,520 / 38,272 /
# 142,912 pairs.  The before-numbers those sweeps produced are deliberately not quoted: the
# re-encodings are drawn pseudorandomly, so how many pairs move depends on WHICH orders were drawn --
# three independent draws give three different counts -- while the zeros above do not depend on the
# draw and every draw tried agrees on them.  The order set that is fully specified is the test file's:
# `_ring_orders`, every permutation of the ring carbons crossed with three placements of the
# chlorines, where `test_one_molecule_in_two_atom_orders_reads_the_same` asserts the zero.
#
# THE FIX IS TO RE-BASE THE PARITY ONTO A FRAME THE COLOURING NAMES, which is the trick ruling F89
# already licenses for the group membership: the term becomes a function of derived data only.
# `_direction_key` keys a direction by its atom's COLOUR (with distinct sentinels for a direction
# with no atom and for no direction at all), and corresponding atoms of two encodings carry the same
# colour -- compute_atoms_order numbers its classes by ascending invariant key, not by slot.  So
# "the directions in ascending colour key" names the SAME frame in every encoding, and the parity
# read in it is a fact about the molecule.  Where two of a unit's directions share a key the frame is
# genuinely unnamed -- swapping those two directions flips the parity, and nothing here can say which
# way round they go -- so the unit contributes NO parity that round.  That costs nothing at the
# fixpoint: refinement only splits classes, so a pair separated in one round stays separated in every
# later one, and a pair that is never separated is one whose exchange the ambiguity report already
# owns (ruling F89).
#
# Nothing is written back.  `translate_parity` and the direction keys are pure arithmetic on a value
# the caller passes in, so this whole path is the read ruling F65 requires -- no parity is stored in
# any frame but the one the arena already holds.
#
# THE STANDING REQUIREMENT FOR ANY FUTURE SEED TERM, and the one line that would have prevented three
# rounds of this: A SEED TERM MUST BE sigma-EQUIVARIANT, NOT MERELY ENCODING-INVARIANT.  It has to
# take equal values at u and at sigma(u) for EVERY automorphism sigma of the annotated molecule, not
# just the same value each time this build reads one fixed encoding.  That is exactly how the stored
# byte failed: it is perfectly stable for a fixed atom order -- read it twice, get it twice -- and
# meaningless across atom orders, because sigma carries a unit's directions to the corresponding
# directions of its image while the SLOTS those directions occupy are the caller's.  A term computed
# from derived data (a refinement class, a fixpoint label, a parity re-based onto a frame those name)
# is equivariant because sigma preserves the derived data itself; a term computed from anything the
# arena stores per slot is not, however invariant it looks in one encoding.  Both tests are cheap to
# run and only the first one is the ruling: a sweep over re-encodings of one molecule catches an
# encoding-variant term, and a sweep over molecules WITH SYMMETRY catches a term that is
# encoding-invariant but not equivariant.


cdef inline uint32_t _frame_free_parity_code(stereo_unit_t *u, uint8_t parity,
                                            uint32_t *colour) noexcept nogil:
    """`parity` re-expressed in the frame `colour` names, as a seed digit: 2 even, 3 odd, 1 unnamed.

    1 means "this unit carries a configured parity whose frame this colouring cannot name" -- two of
    its directions share a colour key.  It is still information, and invariant information: WHETHER a
    parity is configured does not depend on the atom order, only its value does.

    Two shapes, because the two geometries admit different reorderings, and both are read off the
    same key:

    * SU_TETRA has one list of four directions and any permutation of it is available, so the frame
      is "the four slots in ascending key" and the answer is the stored parity translated by the
      permutation that gets there (`translate_parity`, the same arithmetic `translate_stereo` does).
    * the bond kinds have two ordered pairs, and ruling F56 is that the wholesale pair exchange is
      EVEN and contributes nothing while each within-pair swap is odd.  So the frame is "each pair in
      ascending key", the two pairs need no canonical order at all, and the answer is the stored
      parity XOR one bit per pair that had to be reversed.

    A pair whose second slot holds an unnamed direction or a pinned non-direction always sorts
    reversed (those keys are 1 and 0, below every named atom's), so such a unit's code is a fixed
    flip of its stored value.  That is not a mistake to correct: the frame is a DEFINITION, and any
    definition that is the same in every encoding serves.  Ruling F55's within-pair pin is about
    whether a caller's requested order is legal, and no order is being requested here.
    """
    cdef uint32_t perm[4]
    cdef uint32_t key[4]
    cdef uint32_t i, j, best, bestkey, l, pp = 0
    if not _directions_separated(u, colour):
        return 1                            # per direction LIST, which is exactly the gate needed
    if u.kind == SU_TETRA:
        for i in range(4):
            key[i] = _direction_key(u, colour, i)
        for i in range(4):                  # perm[i] = the slot holding the i-th smallest key
            best = 0
            bestkey = 0xFFFFFFFF
            for j in range(4):
                if key[j] < bestkey:
                    bestkey = key[j]
                    best = j
            perm[i] = best
            key[best] = 0xFFFFFFFF          # a real key is colour + 2 <= n + 2, never this
        return translate_parity(parity, perm) + 1
    for l in range(2):
        if _direction_key(u, colour, 2 * l) > _direction_key(u, colour, 2 * l + 1):
            pp ^= 1
    return (((parity - 1) ^ pp) + 1) + 1


cdef void _frame_free_parity_seed(Structure structure, stereo_unit_t *units, uint32_t nunits,
                                  uint32_t *partner, uint32_t *colour, uint32_t *par,
                                  uint32_t n) noexcept nogil:
    """Fill `par[0:n]` with each atom's parity code in the frame `colour` names; 0 where there is none.

    The code lands on the atom or atoms the unit is NAMED ON, and for a bond kind that is both ends.
    Ruling F45 is explicit that a bond kind may anchor at either end and that the choice follows slot
    order -- "there is deliberately no attempt to canonicalise it" -- so which terminal is keyed in
    SEG_PARITY is itself frame-relative, and a digit written to the anchor alone would smuggle the
    encoding back in through the atom it landed on.  Writing it to both ends is symmetric and so
    invariant.  SU_ALLENE needs nothing extra: its anchor is the chain's centre, which is the same
    atom in every encoding, and `stereo_unit_partner` answers SU_NO_REF for it.

    Combined with `max` rather than assignment, so that an atom reached twice does not depend on
    which unit the table lists first (a list order that is slot order, hence frame-relative).  Two
    units can reach one atom: an atropisomer pivot bearing a ring double bond is a cis/trans terminal
    too, which is the collision ruling F45's relocation exists for.  No claim is made here that it
    cannot happen -- `max` makes the answer order-free either way.

    A configured parity on an atom that anchors NO unit is ignored: a parity is a statement about a
    frame, and an atom with no unit has no frame for it to be a statement about.
    """
    cdef uint32_t i, code, a
    for i in range(n):
        par[i] = 0
    for i in range(nunits):
        a = units[i].anchor
        if a >= n or not structure_parity_at(structure, a):
            continue
        code = _frame_free_parity_code(&units[i], structure_parity_at(structure, a), colour)
        if code > par[a]:
            par[a] = code
        if partner[i] != SU_NO_REF and code > par[partner[i]]:
            par[partner[i]] = code


# --- the canonical search's stereo seam (`_canonical.pxi`'s two hooks) --------------------------
#
# WHY THE CANONICAL SEARCH NEEDS THESE AND NOT JUST THE CERTIFICATE.  `mol_identity_bytes` below reads
# the parity digits AFTER `mol_canonical_order` has chosen a labelling, which is sound only if the
# labelling itself is a function of the molecule.  Without these hooks it is not, on any molecule whose
# constitutional symmetry inverts a parity: the search's orbit prune calls two such labellings
# interchangeable, keeps the one with the lower slot, and the digits then differ between two encodings
# of one compound.  These two functions are how the search sees a configuration.
#
# THEY LIVE HERE AND NOT THERE because `_canonical.pxi` is included first and must not know what a
# parity is -- the same reason `mol_certificate_words` takes its stereo term as an opaque array.  The
# pointers are installed at the bottom of this fragment, at import, once.


cdef int _canon_stereo_prepare(Structure structure) except -1:
    """`_canon_prepare_hook`: build the unit table before the search takes any arena pointer.

    The UNMARKED variant, for the reason `mol_identity_bytes` and `canonical_stereo_group_ids` give:
    what the digits read is constitution plus the anchor's parity byte in SEG_PARITY, never a
    stereogenicity mark, so making every canonical order pay a budgeted witness search per unit
    would buy nothing -- and
    would let a canonical labelling move because a symmetry search ran out of budget, which is the
    worst kind of dependency to introduce into a hash.
    """
    return ensure_stereo_units_unmarked(structure)


cdef bint _canon_stereo_digits(Structure structure, uint32_t *colour, uint32_t *digits_out,
                               uint32_t *scratch, uint8_t *unnamed_out) noexcept nogil:
    """`_canon_stereo_hook`: fill `digits_out[0:n]` with ruling F95's frame-free parity codes.

    Returns whether `colour` NAMES every configured unit's frame -- which is what the caller gates its
    orbit prune on, so the answer must err towards False and never towards True.  Two ways it can be
    False, and the second is the loose one:

      * some unit's code came back 1, "a parity is configured and this colouring cannot read its
        value" -- two of the unit's directions share a colour key.  This is the mirror case: the two
        ring branches leaving a carbinol carbon of cis-cyclobutane-1,3-diol share a colour in every
        colouring the constitution admits, so nothing coarser than a discrete leaf can name that frame.
      * two units landed a digit on ONE atom.  `_frame_free_parity_seed` combines them with `max` and
        neither code is recoverable afterwards, so the digit does not determine the configuration and
        the caller may not prune on a colouring refined by it.  The collision needs an atropisomer
        pivot that is also a cis/trans terminal (the case ruling F45's relocation exists for); refusing
        to prune there costs nodes on a molecule that has one and nothing on any other.

    `scratch` is 2n uint32_t: the per-unit partners `_frame_free_parity_seed` wants, then a per-atom
    count of how many units reached it.  Nothing is allocated -- this runs once per tree node.

    `unnamed_out`, when the caller supplies it, receives WHICH ATOMS a False answer is about: 1 on an
    atom whose own digit is unreadable or collided, 2 on the rest of such a unit's frame -- its anchor,
    its partner and its directions.  A symmetry that fixes every marked atom carries each unnamed unit
    onto itself with every direction fixed, so it cannot invert the parity the digits cannot read, and
    the caller may prune with it after all.  The two values keep the marking free of slot order: only a
    1 makes a unit contribute its frame, so a 2 written by one unit can never recruit the next.

    THE SINGLE DEFINITION IS `_frame_free_parity_seed`, called rather than reimplemented, because the
    digits the SEARCH maximises and the digits `mol_identity_bytes` PUBLISHES have to be the same
    function of the same colouring.  Two copies that agreed today would be a labelling chosen under
    one definition and read out under another the first time either moved.
    """
    cdef uint32_t n = structure.header.atom_count
    cdef uint32_t nunits = structure_stereo_unit_count(structure)
    cdef stereo_unit_t *units
    cdef uint32_t *partner = scratch
    cdef uint32_t *reached = scratch + n
    cdef uint32_t i, j, a, p, r
    cdef bint named = True
    for i in range(n):
        digits_out[i] = 0
        reached[i] = 0
    if unnamed_out is not NULL:
        memset(unnamed_out, 0, <size_t> n)
    if not nunits:
        return True
    if nunits > n:
        # Unreachable by the anchor no-collision invariant (`ensure_stereo_units_unmarked` sizes its
        # scratch at one record per atom on the strength of it).  Checked because the alternative to a
        # conservative answer here is a write past the end of a borrowed buffer.
        return False
    units = structure_stereo_units(structure)
    for i in range(nunits):
        partner[i] = stereo_unit_partner(structure, &units[i])
        a = units[i].anchor
        if a >= n or not structure_parity_at(structure, a):
            continue
        reached[a] += 1
        p = partner[i]
        if p != SU_NO_REF and p < n:
            reached[p] += 1
    _frame_free_parity_seed(structure, units, nunits, partner, colour, digits_out, n)
    for i in range(n):
        if reached[i] > 1 or (reached[i] and digits_out[i] == 1):
            named = False
            if unnamed_out is NULL:
                break
            unnamed_out[i] = 1
    if named or unnamed_out is NULL:
        return named
    for i in range(nunits):
        a = units[i].anchor
        p = partner[i]
        if a >= n:
            continue
        if unnamed_out[a] != 1 and not (p != SU_NO_REF and p < n and unnamed_out[p] == 1):
            continue
        if not unnamed_out[a]:
            unnamed_out[a] = 2
        if p != SU_NO_REF and p < n and not unnamed_out[p]:
            unnamed_out[p] = 2
        for j in range(4):
            r = units[i].refs[j]
            if r != SU_NO_REF and r < n and not unnamed_out[r]:
                unnamed_out[r] = 2
    return False


# INSTALLED AT IMPORT, ONCE, and never cleared.  Written here rather than passed as a parameter at the
# four `mol_canonical_order` call sites so that there is exactly one canonical order in the process: a
# hook a caller could forget would give `mol_identity_bytes` and the SMILES writer two different
# labellings of one molecule, and they would then disagree about which molecules are equal.
_canon_stereo_hook = _canon_stereo_digits
_canon_prepare_hook = _canon_stereo_prepare


cdef inline bint _sg_key_less(uint32_t b1, uint32_t b2, uint32_t *count, uint32_t *off,
                              uint32_t *memb) noexcept nogil:
    """Is group byte b1's membership key below b2's?  Key = (member count, member classes ascending).

    A total order on the keys and nothing else: the byte values themselves are never compared, so two
    groups with the same key stay tied and one canonical id cannot depend on which id the caller
    happened to store.
    """
    cdef uint32_t i
    if count[b1] != count[b2]:
        return count[b1] < count[b2]
    for i in range(count[b1]):
        if memb[off[b1] + i] != memb[off[b2] + i]:
            return memb[off[b1] + i] < memb[off[b2] + i]
    return False


cdef int canonical_stereo_group_ids(Structure structure, uint8_t *ids_out,
                                    uint8_t *amb_out) except -1:
    """Canonical id per SEG_STEREO_GROUPS byte value, into a caller-owned 256-byte array.

    `ids_out[b]` receives the canonical group id of the group whose stored byte is `b`, and 0 for a
    byte no atom of this molecule carries.  Indexing by the whole byte rather than by the group
    number keeps the kinds apart: OR 1 and AND 1 are different groups that share a group number, and
    they get their own canonical ids.

    `amb_out`, when not NULL, receives the AMBIGUITY CLASSES: a 1-based class number per stored byte,
    0 for a byte whose canonical id is pinned.  Two bytes sharing a class number are groups this
    molecule's own symmetry can exchange, so which of their ids each one got is arbitrary -- see
    ruling F89 below.  The ids inside one class are contiguous, so the id BLOCK a class owns is
    invariant even though its contents are not, and the CLASS NUMBERS ascend with the smallest
    canonical id of each class (ruling F92, argued at the ambiguity pass) -- never with the stored
    byte, which is the caller's and not this molecule's.
    Only the two numbered kinds are ever reported: kinds 0 and 1 report a stored number instead.

    Reads the arena and writes nothing to it (ruling F79).  Every canonicalisation call below is a
    read as well, so a caller holding a `copy()` that shares this arena cannot observe the call.

    RULING F80, in three steps, and RULING F89, which is why step 3 is a fixpoint and why the
    residual is exposed instead of broken:

    1. The canonical order is seeded with the ENCODING-INVARIANT part of the stereo state -- whether a
       parity is configured, its value RE-BASED ONTO A FRAME THE COLOURING NAMES (ruling F95, argued
       above the helpers: the stored byte is relative to ruling F26's slot order and so is a fact
       about the caller's atom order, not about the molecule), and the group's KIND -- and never with
       the stored group id or the stored parity byte, both of which are things being replaced.  A seed
       containing the group id would make the order depend on the ids and the ids depend on the order;
       a seed containing the raw byte made two encodings of one molecule answer differently while both
       claimed every id was pinned.

       The seed FOLDS IN the constitutional class rather than replacing it: compute_atoms_order takes
       `seed` in place of `_atom_invariant`, not beside it, so a seed of stereo state alone would
       start the refinement with carbon and chlorine in one class and hand the extremal search a much
       coarser tree to break ties in.  `cls[v] * 16` leaves the sixteen stereo codes room underneath:
       four parity codes (none, frame unnamed, even, odd) times four kinds, which is the arithmetic at
       the seed itself.  Every term of the seed must also be sigma-EQUIVARIANT and not merely
       encoding-invariant -- stated as a standing requirement in the ruling F95 block above.

    2. Groups of one kind are ordered by their MEMBERSHIP KEY -- member count, then the members'
       refinement classes ascending -- and only inside one key by their smallest member position.
       The key leads so that each tied block of groups takes a CONTIGUOUS run of ids, which is what
       `canonical_stereo_group_ambiguities` promises its caller: a key outside every class may be
       compared on its own, and a class may be renumbered inside itself without stepping over an id
       that belongs to a group the caller must not touch.

       Measured on octachlorocyclooctane with alternating ring-frame parities and OR groups
       (1, 2, 1, 1, 3, 2, 3, 3) -- the interchangeable triples {C1,C3,C4} / {C5,C7,C8} beside the
       pinned pair {C2,C6} -- over the 120 encodings its test sweeps.  Ranking on position alone
       gives ids {1, 3} to the tied triples and 2 to the PINNED pair, on every encoding: an answer
       that is invariant (with a frame-free seed a pinned group's position cannot move -- an
       automorphism of the fixpoint colouring maps a group only onto a group of its own label) but
       whose one class is not a run of ids.  With the key leading the pair leads at id 1 on its member
       count and the tie is confined to {2, 3}.  So this ordering earns its place on the contiguity
       alone -- position ranking is invariant here and merely gives an answer a caller cannot use.

       Two distinct bytes have disjoint, non-empty member sets, so no two groups share a smallest
       position and the order within a key is total.  That does NOT settle canonicity: it says nothing
       about whether the positions being compared are the same positions on the next read of the same
       molecule.  The order is canonical only up to the
       automorphism group of whatever colouring the seed produced, and a step-1 seed knows each
       atom's KIND but not WHICH ATOMS SHARE A GROUP.  1,2,3,4-tetrachlorocyclobutane with ring-frame
       parities (1, 2, 1, 2) and OR groups (1, 1, 1, 2) is the witness: without step 3 the four ring
       carbons stay in one refinement class, the singleton group lands on whichever canonical position
       the extremal search hands it, and the creation order decides whether it is id 1 or id 2 -- two
       creation orders, one canonical graph, one partition, ids swapped, and nothing ambiguous about
       the molecule (a group of three and a group of one cannot be exchanged by anything).

    3. So the membership itself is fed back, to a fixpoint.  Each round re-refines with a colour per
       atom of (its class this round, a label for the MULTISET OF CLASSES its co-members carry).
       That label is a function of the partition and never of a stored id or a slot, so every round
       is label-invariant; refinement is monotone, so the class count only grows and the loop ends in
       at most n rounds.  A group alone in its label at the fixpoint is pinned: its key differs from
       every other key of its kind, so step 2 ranks it without ever reaching the position tie-break.
       More than two rounds are sometimes needed and are measured to be, by replacing the bound below
       with a literal cap: at one round five of this file's stereo-group fixtures fail, at two rounds
       three, at three rounds exactly one -- always
       test_the_fixpoint_runs_past_a_third_round_when_the_molecule_needs_it, a C6 ring with four
       singleton OR groups and one pair, whose singletons are told apart only by the FOURTH round and
       where a cap of three invents two ambiguity classes the molecule does not have -- and at the full
       bound none.  Over ALL 12,992 (ring-frame parity pattern, group partition) fixtures of that C6
       ring, 7,472 settle after one round, 1,152 need a second, 4,248 a third and 120 a fourth; the C8
       ring with its all-equal and its alternating frame runs 1,748 / 4,560 / 1,956 / 16 over its 8,280
       fixtures.  So a fourth round is not exotic, and n is the only bound the loop can honestly carry.

       RULING F89: when two groups still share a label after the fixpoint -- equivalently, when a
       refinement class still spans both of them, see the proof at the ambiguity pass -- no invariant
       rule can separate them, and this reports that instead of breaking it on creation order.  In the
       cyclobutane above, given OR groups (1, 1, 2, 2) instead, the ring's rotation by two carries
       {C1,C2} onto {C3,C4} and each RING-FRAME parity onto an equal one -- a rotation preserves that
       frame, so it preserves the parities read in it -- so it is an automorphism of the
       parity-annotated molecule and either assignment describes the same mixture: the id -> members
       direction is genuinely ambiguous while the SET of member sets is not.  The two are NOT merged --
       two OR groups describe four stereoisomers where one describes two -- and nothing raises,
       following ruling F62's shape: expose the degraded guarantee, do not refuse the answer.  A ring
       is where that guarantee degrades most often, and not by accident: a centre's next and previous
       ring atoms carry the same colour until something separates them, so its parity has no frame the
       colouring can name and folds in nothing.  Sharing a label is necessary for the exchange
       and not quite sufficient (colour refinement is incomplete), so a rare pinned pair may be
       reported ambiguous; the error is towards saying less than is known, never more.

    4. ABS and unspecified are not numbered: their canonical id is their stored group number.  ABS is
       the single absolute bucket rather than one group among several -- set_stereo_group forces its
       number to 0 -- and kind 0 is the absence of a group.  A forged buffer carrying a nonzero
       number for either is passed through rather than folded to 0, so the view reports what is
       stored instead of inventing a partition, and it reports it for both kinds alike.
    """
    memset(ids_out, 0, 256)
    if amb_out is not NULL:
        memset(amb_out, 0, 256)
    cdef uint32_t n = structure.header.atom_count
    if n == 0 or not structure_has(structure, SEG_STEREO_GROUPS):
        return 0

    cdef uint8_t *sg = structure_stereo_groups(structure)
    cdef stereo_unit_t *units
    cdef uint32_t nunits
    cdef uint32_t i, j, best, kind, c, nlab, nsel, span, nmemb
    cdef int b
    cdef uint32_t min_pos[256]
    cdef uint32_t count[256]
    cdef uint32_t off[257]
    cdef uint32_t cursor[256]
    cdef uint32_t label[256]
    cdef uint32_t sel[256]
    cdef bint any_group = False

    for b in range(256):
        min_pos[b] = 0xFFFFFFFF
    for i in range(n):
        if sg[i]:
            any_group = True
            break
    if not any_group:
        return 0

    # The unit table, for the direction frames ruling F95's parity term is re-based onto.  The
    # UNMARKED variant, which is the third caller of it and of the same shape as the isomorphism
    # kernel (ruling F76): what is read below is constitution -- `kind`, `anchor`, `refs` and
    # `spare`'s unnamed nibble -- plus the anchor atom's parity byte in SEG_PARITY, and never a
    # stereogenicity mark.  Routing this through `ensure_stereo_units` would make every canonical group id pay a
    # budgeted witness search per unit for an answer it does not read.
    #
    # It appends a segment, so it REALLOCATES the arena: it runs before any pointer into the buffer is
    # taken, and `sg` -- taken above for the any_group scan -- is re-borrowed after it (ruling F60).
    ensure_stereo_units_unmarked(structure)
    sg = structure_stereo_groups(structure)
    units = structure_stereo_units(structure)
    nunits = structure_stereo_unit_count(structure)

    # cls, two seed buffers to alternate between, the round's classes, the order, the members
    # grouped by byte, this round's parity codes, and each unit's partner atom.
    cdef uint32_t *cls = <uint32_t *> PyMem_Malloc(<size_t> (8 * <size_t> n + 1)
                                                   * sizeof(uint32_t))
    if cls is NULL:
        raise MemoryError()
    cdef uint32_t *sa = cls + n
    cdef uint32_t *sb = sa + n
    cdef uint32_t *cur = sb + n
    cdef uint32_t *order = cur + n
    cdef uint32_t *memb = order + n         # each group's member classes, grouped by stored byte
    cdef uint32_t *par = memb + n           # ruling F95's parity code per atom, 0..3
    cdef uint32_t *partner = par + n        # per UNIT, and nunits <= n by the anchor invariant
    cdef uint32_t *seed_in = sa
    cdef uint32_t *seed_out = sb
    cdef uint32_t *swap
    cdef uint32_t flags = 0
    cdef Py_ssize_t classes = 0, prev = -1      # the loop below always assigns before reading
    cdef Py_ssize_t rounds = 0                  # against the bound proved at the loop
    try:
        if compute_atoms_order(structure, cls, NULL) < 0:
            raise MemoryError()
        for i in range(nunits):
            partner[i] = stereo_unit_partner(structure, &units[i])
        # Round 0's colouring is `cls` itself -- the stereo-blind refinement classes, which are
        # invariant by construction -- so its parity term is already frame-free and there is no
        # reason to leave it out.  Measured, because ruling F95 asks for no parity here: dropping the
        # term from round 0 and folding it only from round 1 returns the same group view AND the same
        # ambiguity classes on every one of the 127,152 C4-C7 polychlorocycloalkanes (each ring-frame
        # parity pattern by each group partition) and the same all-zero encoding-sweep counters -- so
        # the fixpoint recovers on round 1 exactly what round 0 folds in, and neither spelling is
        # coarser.  It is kept because it is free, and because "the seed is the colouring plus the
        # parity read in that colouring" is one rule for every round rather than one plus an exception.
        _frame_free_parity_seed(structure, units, nunits, partner, cls, par, n)
        for i in range(n):
            # 0..15: (F95 parity code 0..3) * 4 + kind 0..3.  The stored group NUMBER is deliberately
            # absent, and so is the stored parity BYTE (ruling F95, argued above the helpers).
            seed_in[i] = (cls[i] * 4 + par[i]) * 4 + sg_kind(sg[i])

        # Step 3's fixpoint.  TERMINATION, stated because the loop has no other guard: the quantity
        # that grows is `classes`, the number of refinement classes of the round's partition.  Each
        # round seeds compute_atoms_order with (this round's class, label, parity code), which REFINES
        # this round's partition -- the class is the leading digit, so the other two can only split a
        # class further, never merge two -- and compute_atoms_order refines its seed, so classes never
        # falls even though a parity code may appear in a later round than it was absent in; the break
        # fires the first round it fails to rise.  It is bounded by n, so the round count is bounded
        # by n as well: round k cannot be reached unless classes rose on each of the k - 1 rounds
        # before it, so classes >= k, and classes <= n.  The cap below turns that bound into a
        # branch instead of trusting it.
        while True:
            classes = compute_atoms_order(structure, cur, seed_in)
            if classes < 0:
                raise MemoryError()
            if classes == prev:
                break                       # refinement learned nothing new: seed_in is the finest
            prev = classes
            # Re-borrowed per round rather than held across the call (ruling F60).  Nothing in
            # compute_atoms_order appends a segment today; the rule is about what a reader may assume.
            sg = structure_stereo_groups(structure)
            units = structure_stereo_units(structure)

            memset(count, 0, sizeof(count))
            for i in range(n):
                if sg[i]:
                    count[sg[i]] += 1
            off[0] = 0
            for b in range(256):
                off[b + 1] = off[b] + count[b]
                cursor[b] = off[b]
            for i in range(n):
                if sg[i]:
                    memb[cursor[sg[i]]] = cur[i]
                    cursor[sg[i]] += 1
            for b in range(256):            # each group's classes ascending, so equal multisets
                for i in range(off[b] + 1, off[b + 1]):     # compare element by element
                    c = memb[i]
                    j = i
                    while j > off[b] and memb[j - 1] > c:
                        memb[j] = memb[j - 1]
                        j -= 1
                    memb[j] = c
            # Labels are handed out in ASCENDING KEY ORDER, key = (member count, member classes
            # ascending), and never in stored-byte order.  Both consumers of `seed` document that
            # they read its equality classes only, but the class numbering compute_atoms_order builds
            # out of them is ordered BY VALUE and the extremal search sees that order, so a label
            # taken from the stored byte puts the caller's ids back into the answer through the side
            # door.  Measured: with first-seen-byte labels, exchanging the two stored ids of
            # test_the_canonical_view_survives_a_relabelling_and_a_swap_of_the_stored_ids swaps the
            # canonical ids.  The key is invariant because a refinement class number is (its own
            # docstring) a function of the keys and not of any slot.
            nsel = 0
            for b in range(256):
                label[b] = 0
                if count[b]:
                    j = nsel                # insertion sort of the present bytes by key
                    while j > 0 and _sg_key_less(<uint32_t> b, sel[j - 1], count, off, memb):
                        sel[j] = sel[j - 1]
                        j -= 1
                    sel[j] = <uint32_t> b
                    nsel += 1
            nlab = 0
            for i in range(nsel):
                if i and not _sg_key_less(sel[i - 1], sel[i], count, off, memb):
                    label[sel[i]] = nlab    # equal keys: one label, so the groups stay tied
                else:
                    nlab += 1
                    label[sel[i]] = nlab
            # (class, co-member multiset, the parity read in THIS round's colouring) as one integer.
            # A label is handed out per PRESENT BYTE, of any kind -- step 4 makes a forged kind-0 byte
            # reportable, so all 255 nonzero bytes can be present -- hence nlab <= 255 and span <=
            # 256; with the parity code's four values under it, `(cur[i] * span + label) * 4` stays
            # inside uint32 for any molecule under 2**32 / 1024 = 4.19 M atoms.
            span = nlab + 1
            _frame_free_parity_seed(structure, units, nunits, partner, cur, par, n)
            for i in range(n):
                seed_out[i] = (cur[i] * span + label[sg[i]]) * 4 + par[i]
            swap = seed_in
            seed_in = seed_out
            seed_out = swap
            rounds += 1
            if rounds >= <Py_ssize_t> n:
                # The hard bound, and the PROOF that taking the current assignment answers rather
                # than degrades -- the epic's rule is that a bound may not be written without one.
                # By the paragraph at the loop, round k is reachable only if classes rose on each
                # round before it, so classes >= k; a round numbered n therefore has classes == n
                # exactly and the partition is DISCRETE, one atom per class.  Refining a discrete
                # partition returns it unchanged, so the next round would find classes == prev and
                # break at the top -- before recomputing anything.  The state this break leaves,
                # `seed_in` just swapped in and `label` from this round, is bit for bit the state
                # that break would leave, so on the reachable path the bound is not a policy at all.
                # And if compute_atoms_order ever stopped refining its seed, this would fire with the
                # partition still coarse: more groups sharing a label, hence MORE groups reported as
                # tied by amb_out, and the assignment is still label-invariant because a label is
                # still a function of the partition.  A coarser answer that names its own imprecision
                # is ruling F62's shape.  A hang is not.
                break

        # Raises AutomorphismBudgetExceeded on a truncated search and writes nothing: a labelling
        # from a truncated extremal search is a different labelling, not an approximate one, and a
        # group id built on it would be a wrong answer that reads like a right one.
        mol_canonical_order(structure, seed_in, order, &flags, True)
        # Re-borrowed after the calls above rather than reused across them (ruling F60).  None of
        # them canonicalises through the arena -- each allocates its own scratch and appends no
        # segment -- but the rule is about what a reader may assume, not about what today's callee
        # happens to do.
        sg = structure_stereo_groups(structure)
        for i in range(n):
            if sg[i] and order[i] < min_pos[sg[i]]:
                min_pos[sg[i]] = order[i]
    finally:
        PyMem_Free(cls)

    # Kinds are numbered independently, each densely from 1, by ASCENDING (fixpoint label, smallest
    # member position) -- step 2's ordering, whose measurement lives in the docstring.  The label
    # leads so that a tied block of groups takes a CONTIGUOUS run of ids: ranked on position alone,
    # the octachlorocyclooctane fixture returns its two interchangeable triples as ids {1, 3} with the
    # PINNED pair between them at 2, on all 120 encodings, and then a class is not a run and a caller
    # renumbering one has to step over an id it does not own.
    # A selection sort over at most 64 candidate bytes per kind, twice: this runs once per call, not
    # per atom, and the array it sorts is bounded by the encoding rather than by the molecule.
    cdef uint8_t nxt
    for kind in range(2, 4):
        nxt = 1
        while True:
            best = 0xFFFFFFFF
            c = 0xFFFFFFFF
            b = -1
            for i in range(64):
                j = (kind << 6) | i
                if count[j] and (label[j] < c or (label[j] == c and min_pos[j] < best)):
                    c = label[j]
                    best = min_pos[j]
                    b = <int> j
            if b < 0:
                break
            ids_out[b] = nxt
            count[b] = 0                # taken: the next pass must not find it again.  From here on
            nxt += 1                    # PRESENCE is min_pos[b] != 0xFFFFFFFF, not count[b]
    for b in range(256):                # step 4: kinds 0 and 1 keep the stored number, as stored
        if (b >> 6) < 2 and min_pos[b] != 0xFFFFFFFF:
            ids_out[b] = <uint8_t> (b & 0x3f)

    if amb_out is NULL:
        return 0
    # F89: at a fixpoint, two groups share a refinement class if and only if they share a LABEL.
    # (=>) a class holding a member of each gives both atoms the colour (class, own label), and a
    # fixpoint cannot leave two atoms of different colour in one class, so the labels are equal;
    # (<=) equal labels mean equal member-class multisets, so every class of one is a class of the
    # other.  So the label is the whole ambiguity relation and no union-find over the classes is
    # needed: one over these classes computes exactly this partition.
    #
    # Only the two NUMBERED kinds can be ambiguous: kinds 0 and 1 report the stored number, which is
    # nothing this function chose.
    #
    # RULING F92, and the reason the outer loop walks CANONICAL IDS rather than stored bytes: classes
    # are numbered by ascending (kind, smallest canonical id in the class), so the numbering -- and
    # the tuple order a caller builds from it -- is a function of the molecule.  That key is exactly
    # invariant, not merely usually so, and the argument is one line: the ambiguity is by
    # construction confined WITHIN a class (each class is one label block, and step 2 hands a label
    # block a contiguous run of ids whose position is decided by the label's key), so the MINIMUM id
    # over a whole class is unchanged by any permutation the ambiguity permits, and permuting inside
    # one class is the only freedom there is.  Classes are disjoint sets of ids, so no two share a
    # minimum and the order is total.
    #
    # DO NOT WALK `for i in range(64)` HERE.  That numbers the classes in ascending STORED BYTE order,
    # so relabelling the stored numbers returns the same frozensets in a DIFFERENT order -- the caller's
    # ids leaking back into a public answer, the same failure as the Major one level up.  The witness is
    # octachlorocyclooctane with alternating ring-frame parities, one OR pair and six singletons forming
    # three tied pairs; `canonical_stereo_groups()` stays identical key for key across its relabellings.
    cdef uint32_t nclass = 0
    cdef uint32_t byid[65]                  # canonical id -> stored byte, rebuilt for each kind
    for kind in range(2, 4):
        for i in range(65):
            byid[i] = 0xFFFFFFFF
        for i in range(64):
            j = (kind << 6) | i
            if min_pos[j] != 0xFFFFFFFF:
                byid[ids_out[j]] = j        # ids_out is a bijection onto 1..k for these two kinds
        for i in range(1, 65):
            j = byid[i]
            if j == 0xFFFFFFFF or amb_out[j]:
                continue                # absent, or already spoken for by a smaller id of its class
            nmemb = 0
            for b in range(kind << 6, (kind << 6) + 64):
                if b != <int> j and min_pos[b] != 0xFFFFFFFF and label[b] == label[j]:
                    nmemb += 1          # the whole kind: byte order does not track id order
            if not nmemb:
                continue                # alone in its label: its id is pinned
            nclass += 1
            amb_out[j] = <uint8_t> nclass
            for b in range(kind << 6, (kind << 6) + 64):
                if b != <int> j and min_pos[b] != 0xFFFFFFFF and label[b] == label[j]:
                    amb_out[b] = <uint8_t> nclass
    return 0


cdef bytes mol_identity_bytes(Structure structure):
    """The molecule's canonical form as bytes: what `==` compares and what `hash()` hashes.

    Two parts, in this order, and both read off CANONICAL POSITIONS so that neither mentions the
    caller's atom order:

      * `mol_certificate_words`' graph string -- per position, the atom's own invariant word
        (element, isotope, charge, radical, implicit hydrogen count, ring membership) and its bonds
        to higher positions with their orders and aromatic bits;
      * one parity digit per position, ruling F95's frame-free code, so that two molecules
        differing only in a configured parity do not compare equal.

    WHY NOT `signature`, WHICH IS THE OBVIOUS CANDIDATE AND IS WRONG.  `signature` is the OR of
    every atom's four feature words -- a molecule-level screen.  Propane, butane and pentane all
    return `(866942928268820480, 4611686018427387904, 9223372174432275457, 72198331526283521)`,
    because every field the OR keeps (elements present, degrees present, charges, hybridizations,
    bond orders) is identical across the three; `signature`-based equality reports
    `smiles('CCC') == smiles('CCCC')` as equal.  A screen answers "may these match" and is built to be
    permissive; equality needs the opposite bias.

    WHY NOT A CANONICAL SMILES STRING.  A string is only as sound as the labelling behind it, and a
    labelling that takes the FIRST discrete leaf rather than the extremal one is a function of the
    input order and not of the molecule: sixty relabelings of one cubane skeleton produce forty
    distinct strings that way (measured in `_canonical.pxi`'s fragment comment).  The words below come
    from the extremal search, which is what makes a hash sound.

    WHY THE PARITY DIGIT NEEDS A STEREO-AWARE LABELLING AND NOT A REFINED STEREO-BLIND ONE.  The digit
    is read in the frame the CANONICAL POSITIONS name, so positions pinned only up to the automorphism
    group of a STEREO-BLIND colouring are not enough: where a molecule's own symmetry exchanges two
    stereo units carrying different parities, which one takes which position is arbitrary and the digit
    strings of two encodings of one molecule can differ -- a FALSE NEGATIVE.

    THE MEMBERSHIP FIXPOINT IS NOT THE CURE, and that is the trap in this whole area.  A fixpoint --
    `canonical_stereo_group_ids` above, seed, refine, fold back, repeat, of which the seed below is
    round zero -- can only split a tie some colouring CAN name, and in the mirror case there is none:
    on cis-cyclobutane-1,3-diol the two ring branches leaving a carbinol carbon are interchangeable in
    EVERY colouring the constitution admits, so every round returns the same classes and the anchor's
    frame stays unnamed however many rounds are run.  The tie is not one the colouring failed to
    split, it is one the colouring cannot see, and no amount of refining a stereo-blind invariant
    reaches it.

    WHAT SETTLES IT is the extremal search itself being stereo-aware, in `_canonical.pxi`: the leaf
    certificate carries a PARITY TAIL, and the orbit prune feeds `mol_automorphisms` the parity-refined
    colouring `cls * 4 + digit` instead of the constitutional one, so a symmetry only collapses two
    candidates when it preserves configuration too -- and where a configured unit's frame cannot be
    named at all the prune stands down rather than guessing.  The digits come from
    `_canon_stereo_digits` just below, which calls the same `_frame_free_parity_seed` used here, so the
    search's reading of a parity and this function's reading of it cannot drift apart.
    `test_canonical_mirror.py` pins the invariants directly.

    ONE STRUCTURAL GAP IS LEFT, in the same direction: a false negative, never a false positive.
    `_frame_free_parity_seed` combines digits with `max`, so when two units land a digit on one atom
    the digits merge irrecoverably; `_canon_stereo_digits` detects that and declines to prune, which
    keeps the search honest, but the tail it then compares is computed from merged digits.  No witness
    for it has been found.  The asymmetry of the failure mode is what makes the gap tolerable: equal
    molecules reported unequal costs a cache miss, unequal molecules reported equal would corrupt a
    dict, and that direction stays closed because a digit is a function of the position frame and the
    frame is a function of the structure.

    Raises `AutomorphismBudgetExceeded` through `mol_canonical_order` on a truncated extremal search,
    and there is no degraded answer for the reason stated there.  An empty molecule hashes as `b''`.
    """
    cdef uint32_t n = structure.header.atom_count
    cdef uint32_t *units_scratch = NULL
    cdef uint64_t *cert = NULL
    cdef stereo_unit_t *units
    cdef uint32_t *cls
    cdef uint32_t *order
    cdef uint32_t *par
    cdef uint32_t *partner
    cdef uint32_t *seed
    cdef uint32_t nunits, i
    cdef uint32_t flags = 0
    cdef size_t cert_len
    if n == 0:
        return b''

    # The UNMARKED unit table, as in `canonical_stereo_group_ids` and for the same reason: what is
    # read here is constitution plus the anchor's parity byte in SEG_PARITY, never a stereogenicity mark, so
    # making an equality test pay a budgeted witness search per unit would buy nothing.  It appends
    # a segment and so REALLOCATES the arena (ruling F60) -- hence before any pointer is taken.
    ensure_stereo_units_unmarked(structure)
    units = structure_stereo_units(structure)
    nunits = structure_stereo_unit_count(structure)

    cert_len = mol_certificate_len(structure, True)
    units_scratch = <uint32_t *> PyMem_Malloc(<size_t> 5 * n * sizeof(uint32_t))
    if units_scratch is NULL:
        raise MemoryError('canonical identity scratch allocation failed')
    cls = units_scratch
    order = cls + n
    par = order + n
    partner = par + n           # per UNIT, and nunits <= n by the anchor invariant
    seed = partner + n
    try:
        cert = <uint64_t *> PyMem_Malloc(cert_len * sizeof(uint64_t))
        if cert is NULL:
            raise MemoryError('canonical identity certificate allocation failed')
        if compute_atoms_order(structure, cls, NULL) < 0:
            raise MemoryError('atom order refinement failed to allocate')
        for i in range(nunits):
            partner[i] = stereo_unit_partner(structure, &units[i])
        # Round 0 of `canonical_stereo_group_ids`' construction: the stereo-blind classes, plus the
        # parity read in the frame those classes name.  Folded INTO the class rather than beside it
        # (`cls[i] * 4 + par[i]`) because compute_atoms_order takes a seed in place of the atom
        # invariant, so a seed of parity alone would start the refinement with carbon and chlorine
        # in one class.
        # The seed is not what carries the identity: the search refines by parity at every node, so it
        # would reach a stereo-distinguishing labelling from a bare colouring.  It stays here
        # because it is free at this point (`cls` is computed either way), it is sigma-equivariant
        # by ruling F95, and starting the refinement already split costs the extremal search tree
        # nodes it would otherwise have to branch through.  It is an optimisation, not a mechanism.
        _frame_free_parity_seed(structure, units, nunits, partner, cls, par, n)
        for i in range(n):
            seed[i] = cls[i] * 4 + par[i]
        mol_canonical_order(structure, seed, order, &flags, True)
        # The digits that go INTO the string are re-read in the POSITION frame, which is discrete --
        # so a centre whose parity had no frame the coarse colouring could name still contributes a
        # real digit here rather than the "frame unnamed" code.  This is the step that makes
        # enantiomers of an asymmetric skeleton compare unequal.
        for i in range(n):
            cls[i] = order[i] + 1
        _frame_free_parity_seed(structure, units, nunits, partner, cls, par, n)
        mol_certificate_words(structure, order, par, cert)
        return <bytes> (<char *> cert)[:cert_len * sizeof(uint64_t)]
    finally:
        PyMem_Free(cert)
        PyMem_Free(units_scratch)
