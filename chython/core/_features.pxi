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
# The feature-word encoding: four u64 per atom, one u64 per half-edge, and the element
# index.  These are derived segments -- `structure_from_bytes` clears them and
# `rebuild_derived` fills them.
#
# Shared vocabulary rather than molecule-only, which is why the file is not named for a side.
# The molecule side writes these words (`fill_features`, `fill_edge_words`) and the query side
# reads the same layout to build its forbidden masks, so `_bit_of` and the W0_/W1_ span
# constants are the contract between the two.


cdef inline uint32_t element_bucket_begin(Structure structure, uint32_t element) noexcept nogil:
    return structure_element_index(structure)[element]


cdef inline uint32_t element_bucket_end(Structure structure, uint32_t element) noexcept nogil:
    return structure_element_index(structure)[element + 1]


cdef inline uint32_t _bit_of(int32_t value, int32_t lo, int32_t hi) noexcept nogil:
    """Saturating value -> bit offset within a one-bit-per-value field."""
    if value <= lo:
        return 0
    elif value >= hi:
        return <uint32_t> (hi - lo)
    return <uint32_t> (value - lo)


# Word 0 named spans and bit positions.  Every bit in a span is mutually exclusive with
# every other bit in the same span (exactly one fires per bond per traversal step).  The
# topology triple is the "ring-arom / not-ring / ring-plain" classification: aromatic implies
# in-ring, so the three states are disjoint.  The query primitive compiler reads these
# constants too.
DEF W0_ELEMENT_SPAN   = 0x01FFFFFFFFFFFFFF   # bits 0-56
DEF W0_TOPOLOGY_SPAN  = 0x4600000000000000   # bits 57, 58, 62
DEF W0_ORDER_SPAN     = 0xB800000000000000   # bits 59, 60, 61, 63
DEF W0_BIT_RING_AROM  = 57
DEF W0_BIT_NOT_RING   = 58
DEF W0_BIT_RING_PLAIN = 62
DEF W0_BIT_ORDER1     = 59
DEF W0_BIT_ORDER2     = 60
DEF W0_BIT_ORDER3     = 61
DEF W0_BIT_ORDER8     = 63
DEF W1_ELEMENT_SPAN        = 0x3FFFFFFFFFFFFFFF   # word 1 bits 0-61: heavy-element identity span
# Word IV's stereo bit, and the mask that drops it.
#
# NAMED BECAUSE A SECOND READER OF THE WORD NEEDED IT.  `atom_feature_word4` sets bit 6 from the
# SEG_PARITY byte (0 none, 1 even, 2 odd), and its own docstring explains at length why no query box
# may demand that bit: the stored value is a statement in the molecule's ruling-F26 slot frame, a
# query's is in its own, and the two differ by the embedding's permutation.  The consequence nobody
# had written down is that the bit is not frame-free for a MOLECULE-to-MOLECULE comparison either.
# Two spellings of one meso compound -- `C[C@H](O)[C@H](O)C` and `C[C@@H](O)[C@@H](O)C` -- hold
# opposite stored parities at both centres and therefore differ in the union row, while their
# canonical forms are equal, so `__eq__`'s union-word screen rejected a pair that IS equal.  A screen
# is allowed to be lossy in one direction only; a false NEGATIVE is a wrong answer.  Anything
# comparing union rows across two molecules masks with W4_FRAME_FREE_MASK; the per-atom words handed
# to the isomorphism kernel are untouched, because there the frame is the embedding's and the kernel
# handles it.
DEF W4_BIT_STEREO          = 6
DEF W4_FRAME_FREE_MASK     = 0xFFFFFFFFFFFFFFBF   # word IV minus bit 6 (the stored parity)
DEF W0_LIGHT_ELEMENT_SPAN  = 0x01FFFFFFFFFFFFFE   # word 0 bits 1-56: the light elements only
# (bit 0 of word 0 is the heavy-element flag, not an element identity bit; excluded here so
#  that forbidding every bit in W0_LIGHT_ELEMENT_SPAN + the flag together means no element matches)


# The two halves of word 0, factored out because `fill_features` and `fill_edge_words` both
# build them and must agree bit for bit -- the kernel ANDs a query box's neg[0] against either
# one.  Keeping them as functions makes that agreement structural instead of a promise in a
# docstring.

cdef inline uint64_t w0_element_bits(uint32_t element) noexcept nogil:
    """Word 0's element span: one-hot for the light elements, the shared flag for the heavy ones."""
    if element > 56:
        return 1
    elif element:
        return <uint64_t> 1 << (57 - element)
    return 0


cdef inline uint64_t w0_bond_bits(halfedge_t *e) noexcept nogil:
    """Word 0's topology triple and order span for one half-edge."""
    cdef uint64_t w
    # AROMATIC IS TESTED FIRST, and the order matters. Word 0 has no free bit -- bits 0-56 are the
    # element span, 57/58/62 the topology triple, 59/60/61/63 the order span, sixty-four of sixty-four
    # -- so a stored order 4 cannot have an order bit of its own, and word II is full as well (bits
    # 0-61 heavy-element identity, 62/63 the radical pair). What makes that survivable rather than a
    # screening hole is that order 4 and HE_AROMATIC are set together by construction
    # (`_emit_half` writes the flag from the order, `structure_from_bytes` rejects a buffer where they
    # disagree), so W0_BIT_RING_AROM fires EXACTLY when the order is 4 and is the aromatic bond's own
    # value. That is why this branch precedes the ring test: an aromatic bond outside a perceived ring
    # would otherwise take W0_BIT_NOT_RING and become indistinguishable from a dative bond, which is
    # the one wrong answer the layout could still produce.
    #
    # An order-4 bond therefore shares W0_BIT_ORDER8 with a dative bond. The pair separates them, and
    # the obligation that falls out of it is on the query side and is discharged in `_query_boxes.pxi`:
    # a BPRIM_ORDER demand for 8 forbids W0_BIT_RING_AROM as well as the other order bits, so
    # "coordination bond" cannot admit an aromatic one. BPRIM_AROMATIC already demands this bit and
    # BPRIM_RING already accepts it, so those two needed nothing.
    if e.flags & HE_AROMATIC:
        w = <uint64_t> 1 << W0_BIT_RING_AROM
    elif not e.flags & HE_IN_RING:
        w = <uint64_t> 1 << W0_BIT_NOT_RING
    else:
        w = <uint64_t> 1 << W0_BIT_RING_PLAIN
    if e.order == 1:
        w |= <uint64_t> 1 << W0_BIT_ORDER1
    elif e.order == 2:
        w |= <uint64_t> 1 << W0_BIT_ORDER2
    elif e.order == 3:
        w |= <uint64_t> 1 << W0_BIT_ORDER3
    else:                                  # dative (order 8) or aromatic (order 4)
        w |= <uint64_t> 1 << W0_BIT_ORDER8
    return w


cdef inline uint64_t atom_feature_word4(atom_t *a, uint8_t parity) noexcept nogil:
    """Feature word IV for ONE atom: hybridization, the stereo bit, ring sizes and ring counts.

    THE SINGLE DEFINITION (ruling F78).  `fill_features` calls this rather than inlining the
    arithmetic, so `refresh_parity_features` -- which every parity writer that runs after
    `rebuild_derived` calls -- cannot drift from the derivation whose output it has to reproduce.
    A `to_bytes`/`from_bytes` round trip re-derives from scratch and is the oracle: what a writer
    leaves behind must be byte-identical to it, which it is only while there is one formula.

    `parity` is the three-state SEG_PARITY byte for this atom (0 none, 1 even, 2 odd).  An ODD
    parity is screened as bit 6, and whether one is configured at all as the TWO-BIT span in bits
    7 and 8: bit 7 "a parity is configured", bit 8 "none is".  Two bits and not one because the
    span is what a query box can forbid: `wbox_forbid_one_hot` states a demand by forbidding the
    rest of its span, and over a one-bit span there is no rest, so a box could not say "configured"
    at all.

    Bit 6 remains screen-invisible on its own: a box may not demand it, because the stored value is in
    the molecule's ruling-F26 frame and a query's is in its own, and the two differ by the embedding's
    permutation (see prim_apply's PRIM_STEREO branch).  "Configured" carries no frame, so it screens.
    """
    cdef uint64_t w4 = <uint64_t> 1 << (at_hybridization(a) - 1)
    if parity == 2:
        w4 |= <uint64_t> 1 << W4_BIT_STEREO
    if parity:
        w4 |= <uint64_t> 1 << 7
    else:
        w4 |= <uint64_t> 1 << 8
    w4 |= <uint64_t> a.ring_sizes << 22
    w4 |= <uint64_t> 1 << (47 + _bit_of(at_ring_count(a), 0, 8))
    # STILL THE CONSTANT `1 << 56`, and now deliberately so rather than for want of aromatic bonds.
    # `perceive_rings` passes 0 for the aromatic ring count on every atom and only
    # `_atom_field_probe` ever writes a non-zero one, so no corpus can distinguish this term from
    # any other formula over that count.  Storing order 4 did NOT make the count derivable: the
    # ring bitmap is VERTEX-scoped (`_fill_descriptors`), and "every bond of this prototype is
    # order 4" is a question about its EDGES.  The vertex-only approximation -- every vertex of the
    # prototype has two aromatic bonds -- is wrong on biphenylene, whose central four-ring has four
    # such vertices and two single bonds of its own, so it would report an aromatic cyclobutadiene.
    # Deriving it properly means carrying an edge set per prototype, which costs memory on every
    # molecule to feed a term no query can screen on (it is outside every SPAN_MASK entry).  Left
    # constant until something needs the number; whoever needs it must widen the prototype, not
    # this formula.
    w4 |= <uint64_t> 1 << (56 + _bit_of(at_aromatic_ring_count(a), 0, 7))
    return w4


cdef void refresh_parity_features(Structure structure) noexcept nogil:
    """Re-derive feature word IV after a parity was written outside `rebuild_derived` (F78).

    Word IV is the only feature word a parity reaches, so word I..III and their union entries
    are left alone; word IV is recomputed for EVERY atom and the union row rebuilt from those,
    because the union is an OR and a bit that has gone cannot be un-ORed out of it.  O(atoms), paid
    only by a writer that actually changed a bit.

    A no-op when the segment is absent (nothing to keep current) or when its length does not match
    the atom count, which is the only way the per-atom writes below could run off the end.  On return
    word IV and the union row state what SEG_PARITY holds, whatever wrote it: `structure_clone` copies
    derived segments verbatim, so a clone whose parities were rewritten needs this call to agree with
    its own segment.
    """
    cdef uint32_t n = structure.header.atom_count
    cdef uint64_t *feat
    cdef atom_t *atoms
    cdef uint8_t *par
    cdef uint64_t union4 = 0
    cdef uint64_t w4
    cdef uint32_t i, base
    if structure_seg_len(structure, SEG_FEATURES) != 32 * (1 + <size_t> n):
        return
    feat = structure_features(structure)
    atoms = structure.atoms()
    par = NULL
    if structure_has(structure, SEG_PARITY):
        par = structure_parities(structure)
    for i in range(n):
        base = 4 + 4 * i
        w4 = atom_feature_word4(&atoms[i], par[i] if par is not NULL else 0)
        feat[base + 3] = w4
        union4 |= w4
    feat[3] = union4


cdef int fill_features(Structure structure) except -1:
    cdef uint32_t n = structure.header.atom_count
    structure_append(structure, SEG_FEATURES, 32 * (1 + <size_t> n))

    cdef uint64_t *feat = structure_features(structure)
    cdef atom_t *atoms = structure.atoms()
    cdef uint32_t *ptr = csr_ptr(structure)
    cdef halfedge_t *edges = csr_edges(structure)
    cdef atom_t *a
    cdef uint8_t *par = NULL
    cdef uint64_t w1, w2, w3, w4
    cdef uint32_t i, k, element, ih, eh, th, base
    cdef uint32_t degree, heteroatoms
    cdef uint8_t nb_element
    cdef int32_t delta
    if structure_has(structure, SEG_PARITY):
        par = structure_parities(structure)

    with nogil:
        for i in range(n):
            a = &atoms[i]
            element = a.element
            # --- word I: element and incident bonds ---
            #
            # The two COUNTS the `D` and `x` primitives screen on are accumulated here rather than
            # read off `a.degree` / `a.heteroatoms`, and they are NOT the same numbers.  A dative
            # bond is a coordination contact, not a substituent: `[Fe]~N(C)(C)C` is trimethylamine
            # donating its lone pair, and the nitrogen has three substituents in every sense a rule
            # cares about -- which is why `derive_scalars` already gives it `z1`, sp3, counting no
            # dative bond towards hybridization.  `D` and `x` now agree with `z`.
            #
            # `a.degree` and `a.heteroatoms` KEEP counting the contact, because they are structural:
            # `rebuild_derived` derives degree as the CSR row length, `_pach.pxi` writes it as that
            # row length, and `_stereo.pxi` reads it as connectivity.  Two different facts, so two
            # counts -- neither is a stale copy of the other, and each is derived in one place only.
            w1 = w0_element_bits(element)
            degree = 0
            heteroatoms = 0
            for k in range(ptr[i], ptr[i + 1]):
                w1 |= w0_bond_bits(&edges[k])
                if edges[k].order == 8:
                    continue
                degree += 1
                nb_element = atoms[edges[k].to].element
                if element_is_heteroatom(nb_element):
                    heteroatoms += 1

            # --- word II: heavy element and radical ---
            w2 = 0
            if element > 56:
                w2 = <uint64_t> 1 << (element - 57)
            if at_radical(a):
                w2 |= <uint64_t> 1 << 63
            else:
                w2 |= <uint64_t> 1 << 62

            # --- word III: counts, charge, isotope ---
            ih = at_implicit_h(a)
            eh = at_explicit_h(a)
            w3 = <uint64_t> 1 << _bit_of(<int32_t> heteroatoms, 0, 8)
            w3 |= <uint64_t> 1 << (9 + _bit_of(<int32_t> degree, 0, 7))
            w3 |= <uint64_t> 1 << (22 + _bit_of(<int32_t> eh, 0, 4))
            if ih == H_UNKNOWN:
                # EVERY bit of both spans, which in this encoding means "no h or H demand can be
                # satisfied here" and NOT "any of them can".  The kernel's whole atom test is
                # `f[2] & box.neg[2]` (`_box_admits`) and a box states `h2` by forbidding the
                # rest of its span, so an atom carrying the full span is refused by any box that
                # touched it at all -- positive `h2` and negated `!h2` alike.  That is the intended
                # reading: an atom whose implicit hydrogen count nobody recorded cannot answer a
                # question about its implicit hydrogen count, in either direction.
                #
                # DIVERGES FROM chython 2 on the negated form, deliberately.  There `h` compared
                # against `Element.implicit_hydrogens`, so `None != 2` was True and `[C;!h2]` matched
                # an atom with no count at all -- a match granted by the absence of data rather than
                # by the data.  Here it does not match.  The positive form agrees with V2 (`None == 2`
                # was False there too), and the total-H span follows the implicit one because a total
                # is a sum and a sum with an unknown term is unknown.
                #
                # The explicit span above stays EXACT: an explicit hydrogen is an atom someone drew,
                # so its count is known even when the implicit one is not.
                #
                # SPAN_MASK is `_query_boxes.pxi`'s, a LATER fragment, and reaching forward to it is
                # deliberate: those masks are the layout contract between this writer and the query
                # side, so a second spelling here is exactly the drift the SPAN_MASK comment warns
                # about.  Sound because every verbatim extern block is emitted ahead of every
                # function body -- see the include-order section of RULES.md, which records this as
                # the one function-body forward reference in the core.
                w3 |= SPAN_MASK[SPAN_IMPLICIT_H] | SPAN_MASK[SPAN_TOTAL_H]
            else:
                th = ih + eh
                w3 |= <uint64_t> 1 << (17 + _bit_of(<int32_t> ih, 0, 4))
                w3 |= <uint64_t> 1 << (27 + _bit_of(<int32_t> th, 0, 5))
            w3 |= <uint64_t> 1 << (33 + _bit_of(a.charge, -4, 8))
            if a.isotope:
                delta = <int32_t> a.isotope - <int32_t> MDL_ISOTOPE[element]
                w3 |= <uint64_t> 1 << (46 + _bit_of(delta, -8, 8))
            else:
                w3 |= <uint64_t> 1 << 63

            # --- word IV: hybridization, stereo, rings ---
            # In `atom_feature_word4` rather than here, because a parity writer that runs after
            # this pass has to reproduce it exactly (ruling F78) and two copies of the formula
            # would drift.
            w4 = atom_feature_word4(a, par[i] if par is not NULL else 0)

            base = 4 + 4 * i
            feat[base] = w1
            feat[base + 1] = w2
            feat[base + 2] = w3
            feat[base + 3] = w4
            feat[0] |= w1
            feat[1] |= w2
            feat[2] |= w3
            feat[3] |= w4
    return 0


cdef int fill_edge_words(Structure structure) except -1:
    """One u64 per half-edge: the TARGET atom's element bits plus THIS bond's topology and
    order bits, in feature word 0's layout.

    The kernel's inner loop tests a candidate's element and the bond reaching it in a single
    AND against a query box's neg[0], without touching the candidate's feature words at all.
    Word 0's element span is one-hot and so are the topology triple and the order span, so a
    forbidden mask over this word is exact.
    """
    cdef size_t half_edges = <size_t> 2 * structure.header.bond_count
    # Allocate at least one word even when there are no bonds so that the buffer is
    # a real allocation rather than a pointer into the shared read-only zero page.
    # Consequence: structure_has(s, SEG_EDGE_WORD) returns True for a single-atom
    # molecule, so it is NOT a bond-existence test.
    structure_append(structure, SEG_EDGE_WORD, 8 * (half_edges if half_edges else 1))

    cdef uint64_t *words = structure_edge_words(structure)
    cdef atom_t *atoms = structure.atoms()
    cdef halfedge_t *edges = csr_edges(structure)
    cdef halfedge_t *e
    cdef uint32_t k
    with nogil:
        for k in range(half_edges):
            e = &edges[k]
            words[k] = w0_element_bits(atoms[e.to].element) | w0_bond_bits(e)
    return 0


# Fields that survive substructure embedding. chython compares molecule atoms by atomic
# number, isotope, charge and radical only (periodictable/base/element.py:402), and bonds
# by order, so those values must appear in any superstructure. An atom's LOCAL
# environment need not: degree, hydrogen counts, heteroatom count, hybridization and the
# ring descriptors all change when an atom gains neighbours, and ring_sizes is built from
# relevant cycles, which subgraph embedding does not preserve. Screening on any of them
# rejects true substructures -- CCC really is a substructure of CC(C)C, yet propane's
# middle carbon has degree 2 and isobutane has no degree-2 atom at all.
cdef extern from *:
    """
    static const unsigned long long SIG_MASK[4] = {
        0xB9FFFFFFFFFFFFFFULL,   /* word I: every bit but 57, 58 and 62 -- the bond topology
                                    triple. Ring membership and aromaticity both change under
                                    embedding, so neither may screen. Order bit 63 stays. */
        0xFFFFFFFFFFFFFFFFULL,   /* word II: heavy elements and radical, both exact */
        0xFFFFFFFE00000000ULL,   /* word III: bits 33-63 only -- charge and isotope */
        0x0000000000000000ULL};  /* word IV: hybridization, stereo, rings -- none survive */
    """
    const uint64_t SIG_MASK[4]


cdef bint sig_contains(uint64_t *big, uint64_t *small) noexcept nogil:
    cdef uint32_t k
    cdef uint64_t want
    for k in range(4):
        want = small[k] & SIG_MASK[k]
        if (big[k] & want) != want:
            return False
    return True


cdef int fill_element_index(Structure structure) except -1:
    cdef uint32_t n = structure.header.atom_count
    structure_append(structure, SEG_ELEMENT_INDEX, (120 + <size_t> n) * sizeof(uint32_t))

    cdef uint32_t *offset = structure_element_index(structure)
    cdef uint32_t *idx = offset + 120
    cdef atom_t *atoms = structure.atoms()
    cdef uint32_t cursor[120]
    cdef uint32_t i, e
    with nogil:
        # count element e into slot e + 1, so an inclusive prefix sum turns
        # the array directly into bucket starts
        for i in range(n):
            offset[atoms[i].element + 1] += 1
        for e in range(1, 120):
            offset[e] += offset[e - 1]
        # offset[e] is now the start of bucket e, and offset[119] == n
        for e in range(120):
            cursor[e] = offset[e]
        for i in range(n):
            e = atoms[i].element
            idx[cursor[e]] = i
            cursor[e] += 1
    return 0


cdef int rebuild_derived(Structure structure) except -1:
    # THIS CALL INVALIDATES EVERY ARENA POINTER (Ruling F60).  It appends six derived
    # segments, each through structure_append -> PyMem_Realloc, which is free to MOVE the
    # buffer; the old block is then freed.  So no caller may hold an atom_t*, halfedge_t*,
    # uint32_t* or any other pointer into `structure` across this call — re-fetch afterwards.
    # The same applies to ensure_stereo_units / ensure_component_labels and to anything else
    # that builds a lazy segment.  Reading a stale pointer here does not crash: it returns
    # plausible garbage out of freed memory, which is how a corrupt `_numbers` list once
    # went undetected for two fix rounds and was found only by bisecting an "ordering flake".
    #
    # Re-derive degree from the CSR — it may be forged in a packed buffer.
    # Fetch atoms/ptr before any structure_append calls below could realloc.
    cdef atom_t *atoms = structure.atoms()
    cdef uint32_t *ptr = csr_ptr(structure)
    cdef uint32_t n = structure.header.atom_count
    cdef uint32_t i, d
    cdef int rc
    for i in range(n):
        d = ptr[i + 1] - ptr[i]
        atoms[i].degree = <uint8_t> (d if d < 255 else 255)
    with nogil:
        derive_scalars(structure)
        rc = mark_bridges(structure)
    if rc:
        raise MemoryError('bridge detection scratch allocation failed')
    # SEG_RELEVANT_RINGS holds a minimum cycle basis, not the relevant set: the relevant set is
    # exponential (2**20 on a 20-benzene cyclophane) while the prototypes it comes from are
    # polynomial, so the basis is what a caller can be handed as cycles
    perceive_rings(structure)
    fill_features(structure)
    fill_edge_words(structure)
    fill_element_index(structure)
    return 0
