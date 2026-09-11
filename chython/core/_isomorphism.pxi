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
# The subgraph isomorphism kernel: a DFS over query positions, one AND per feature word.
#
# matcher_next is a RESUMABLE STATE MACHINE, not a recursive function, and that is the whole
# shape of this file.  A caller wants embeddings one at a time -- get_mapping is a generator,
# is_substructure stops at the first hit -- so the search has to suspend with its stack intact
# and resume where it left off.  C recursion cannot suspend: the frames are gone the moment the
# function returns.  So the recursion is reified: `depth` is the position being filled,
# `candidate[depth]` is that level's cursor into its candidate source, and `mapping` and `used`
# are the partial assignment.  One call runs the loop until `depth` reaches atom_count (an
# embedding, return True) or underflows 0 (exhausted, return False).  The next call re-enters at
# exactly the state the previous one left, which is why `candidate[depth]` always points PAST
# the atom that level is currently using.
#
# Positions are already in DFS order -- query_seal laid the arena out that way -- so there is no
# ordering work here at all: position p's parent is qatom_t.back and it is always < p.
#
# Component groups (Task 13) are matched.  A query carrying a component group constrains which
# molecule component each query component lands in: same group → same molecule component, different
# groups → different molecule components, no group → unconstrained.


# Module-level cdef rather than DEF: 0xFFFFFFFF does not fit a C int, so a DEF of it is a Python
# object and comparing depth against it inside `nogil` would need the GIL.  Same shape as
# Q_NO_SLOT in _query_seal.pxi.
#
# MATCH_UNSET is written into m.mapping on every backtrack.  closures_admit (Task 11) reads
# m.mapping for closure target positions, but query_seal assigns each closure to the HIGHER of
# its two endpoints (_query_seal.pxi) and positions are filled in ascending DFS order, so
# the target index is strictly less than the owning position and is always mapped when
# closures_admit runs.  Tasks 12-13 rely on the same monotonicity.
#
# MATCH_UNSET is deliberately NOT 0xFFFFFFFF, which is what SU_NO_REF (_stereo.pxi) is.  Values
# read out of m.mapping are compared against a unit record's refs in stereo_admits, where
# SU_NO_REF means "this direction is not an atom"; sharing one bit pattern would let an unmapped
# position pair with that slot instead of being caught.  The readiness ordering makes that
# unreachable today, so this is defence in depth and not a fix -- but a value nobody compares
# against costs nothing to keep distinct.
cdef uint32_t MATCH_UNSET = 0xFFFFFFFD   # mapping[position]: this position is not filled
cdef uint32_t MATCH_DONE = 0xFFFFFFFF    # depth: the search is over; every later call is False


cdef enum:
    CAND_ELEMENT_BUCKET = 0   # scan the target's index for one element
    CAND_ALL_ATOMS = 1        # scan 0 .. n-1: a root that admits more than one element
    CAND_NEIGHBOURS = 2       # scan the CSR half-edges of mapping[qatom_t.back]


cdef struct matcher_t:
    # Every pointer here is borrowed from a Query and a Structure the caller keeps alive for the
    # whole search; nothing is owned except the nine arrays matcher_init allocates.  Because the
    # arena pointers are cached, no arena append may happen during a search -- structure_append
    # reallocates, and a cached pointer would dangle.  A query arena is immutable after seal for
    # exactly this reason, which is why there is no query_append to worry about.  A search that
    # RETURNS TO PYTHON between solutions cannot promise that much, so the two generators that do
    # call matcher_reseat on every resume; a search that runs to completion under nogil needs
    # nothing.
    qatom_t *atoms
    qbox_t *boxes
    qany_t *anys
    # Closure data: closures_admit (Task 11) reads these for positions that own closures.
    qclosure_t *closures
    qbond_t *bonds
    qbox_t *bond_boxes
    # Component group support (Task 13).  group_count is 0 when QFLAG_HAS_GROUP is clear, so
    # group_admits and the anchor helpers short-circuit immediately for ungrouped queries.
    qcomp_t *components
    uint32_t component_count
    int32_t *group_anchor          # one entry per group, -1 = unclaimed; owned (malloc'd)
    uint32_t group_count           # highest group number + 1, or 0 if no component groups
    uint32_t *component_of_position  # position -> component index; owned (malloc'd)
    uint32_t *labels               # owned (malloc'd): molecule component labels, computed fresh
                                   # at init rather than borrowed from the arena: no arena append
                                   # may happen during a search (the no-append rule above)
    # The automorphism group query_seal computed, as permutations of positions, the identity
    # excluded.  Read only at a complete solution, and only when the caller asked for the filter.
    uint32_t *automorphisms
    uint32_t automorphism_count
    bint automorphism_filter
    # Stereo matching (Task 10).  stereo_count is 0 unless the query carries QFLAG_HAS_STEREO, and
    # every cost below -- the unit table included -- is paid only then.  `units`, `satoms`,
    # `parities` and `stereo_groups` point into the MOLECULE arena and so are re-borrowed by
    # matcher_reseat; `stereo` points into the immutable query arena and is not.
    qstereo_t *stereo
    uint32_t stereo_count
    stereo_unit_t *units
    uint32_t unit_count
    uint8_t *sign_state            # owned (malloc'd, only when stereo_count): per position, which
                                   # configurations the boxes that admitted its current candidate
                                   # accept (ruling F87).  Written by matcher_next, where the edge
                                   # word is live; read by stereo_admits, which runs at a LATER
                                   # position than the anchor whenever the frame closes late.
    atom_t *satoms                 # the molecule's atoms, read for whether a direction is a
                                   # hydrogen (ruling F88)
    uint8_t *parities              # SEG_PARITY, or NULL when the molecule configures none: the
                                   # zero page covers 4096 bytes and this array is indexed by atom
    uint8_t *stereo_groups         # SEG_STEREO_GROUPS, or NULL when the molecule has none: the
                                   # zero page covers 4096 bytes and this array is indexed by atom
    bint has_or_group              # some atom of the MOLECULE carries an OR group (ruling F86): the
                                   # group decision at a complete mapping is skipped unless it does,
                                   # so an ABS/AND-only target pays one test of a register per
                                   # solution.  Recomputed by matcher_reseat with the rest of the
                                   # borrowed stereo state, and only when stereo_count is nonzero.
    uint64_t *features             # structure_features(structure) + 4, i.e. past the union row
    uint64_t *edge_words
    uint32_t *csr_begin
    halfedge_t *edges
    uint32_t *element_slots        # structure_element_index(structure) + 120: the bucket contents
    uint32_t *mapping              # position -> molecule atom index; MATCH_UNSET when unfilled,
                                   # though nothing reads that yet -- see the note above
    uint8_t *used                  # molecule atom index -> taken
    uint32_t *candidate            # per-position cursor into its candidate source
    uint8_t *candidate_kind        # one CAND_* per position, decided once at init
    uint32_t *bucket_begin         # per-position element bucket bounds, resolved once at init
    uint32_t *bucket_end
    uint32_t depth
    uint32_t atom_count
    uint32_t structure_atom_count


cdef int _root_element(qbox_t *boxes, uint32_t box_count) noexcept nogil:
    """The single element EVERY box of a root position demands, or -1 when there is not one.

    qatom_t carries no element field, so the DFS seed has to come back out of the box masks.  A
    root's neg[0] holds element bits only -- nothing folded a bond into it -- so the light and
    heavy element sets a box leaves open are exact, and _element_from_neg (_query_seal.pxi) decodes
    them.  This deliberately shares _term_exact_element's strictness: `[C,N]` has two boxes
    allowing one element each, and no single bucket covers it, so it must fall back to scanning
    every atom.  Only the precondition differs -- that one reads a wbox_t with its `touched`
    masks, this one reads the sealed arena, which keeps no `touched` and needs none.
    """
    cdef uint32_t i
    cdef int result = -1
    cdef int e

    if box_count == 0:
        return -1
    for i in range(box_count):
        e = _element_from_neg(boxes[i].neg[0], boxes[i].neg[1])
        if e < 0:
            return -1
        if result < 0:
            result = e
        elif result != e:
            return -1
    return result


cdef inline bint _box_admits(qbox_t *box, qany_t *anys, uint64_t edge_word,
                             uint64_t *f) noexcept nogil:
    """One box against one candidate: four ANDs, then any positive multi-hot demands.

    Factored out of atom_admits so that _stereo_sign_mask can ask the same question of the same
    boxes (ruling F87 puts the stereo sign ON the box, so the two must agree about which boxes
    admitted a candidate, and a second hand-written copy of this test would drift).
    """
    cdef uint32_t k
    if edge_word & box.neg[0]:
        return False
    if f[1] & box.neg[1] or f[2] & box.neg[2] or f[3] & box.neg[3]:
        return False
    for k in range(box.any_count):
        if not (f[anys[box.any_begin + k].word] & anys[box.any_begin + k].mask):
            return False
    return True


cdef inline bint atom_admits(matcher_t *m, uint64_t edge_word,
                             uint32_t position, uint32_t candidate) noexcept nogil:
    """Does molecule atom `candidate`, reached over `edge_word`, satisfy query position?

    The disjunction over boxes short-circuits on the first box that admits, and each box costs
    four ANDs and a branch.  `edge_word` is word 0's comparison partner and which word the caller
    passes is the whole correctness question:

      * a ROOT position's neg[0] holds element bits only, so the aggregate per-atom word 0
        (features[4 * candidate]) is correct -- its element bits are exact, and the root box
        forbids none of its topology or order bits;
      * every NON-ROOT position's neg[0] holds element bits PLUS the bond query_seal folded into
        it, so it must be given the incident HALF-EDGE word from edge_words[k].  The aggregate
        word 0 ORs every incident bond's order bit together, so a folded box forbidding `single`
        would reject acetone's carbonyl carbon and query O=C-C would find nothing.

    Word 0 is never read off the candidate's own row: f[0] is deliberately absent below.  No atom
    primitive writes word 0 outside the element span (an atom's ring membership is word 3), so the
    element bits carried by `edge_word` settle everything word 0 can say.
    """
    # An R (element 0) is a marker, not an atom a query can name.  The refusal is here, at the
    # one gate every query primitive passes through, rather than per primitive: `[A]`, `[#6]` and a
    # bare `c` all arrive as an `atom_admits` call, and a marker admits none of them.
    if m.satoms[candidate].element == 0:
        return False
    cdef qatom_t *qa = m.atoms + position
    cdef qbox_t *boxes = m.boxes + qa.box_begin
    cdef uint64_t *f = m.features + 4 * candidate
    cdef uint32_t b
    for b in range(qa.box_count):
        if _box_admits(&boxes[b], m.anys, edge_word, f):
            return True
    return False


cdef inline uint8_t _stereo_sign_mask(matcher_t *m, uint64_t edge_word,
                                     uint32_t position, uint32_t candidate) noexcept nogil:
    """Which configurations the boxes that ADMIT this candidate are willing to accept (ruling F87).

    A query atom's boxes are a disjunction, and each box carries its own sign, so the answer is a
    union over the admitting boxes only -- the boxes that rejected the candidate on element, degree
    or anything else have no say about its configuration.  atom_admits stops at the first admitting
    box; this one has to see them all, which is why it is a separate pass and not a by-product.

    The three contributions, and why the encoding needs a fourth bit:

      * a box with no stereo primitive accepts either configuration, which is QSIGN_FREE -- and it
        dominates: '[C;@,D3]' matched through its D3 disjunct states nothing about stereo at all,
        so no refusal of any kind may follow.
      * a box naming one sign contributes that sign.
      * a box that ANDed both signs ('[C;@;@@]') accepts nothing and so contributes NOTHING.  This
        is the reason the mask cannot be a plain OR of the box signs: QSIGN_CW | QSIGN_CCW read off
        one box means "no configuration", read off two boxes means "either", and only skipping the
        contradictory box keeps the two apart.

    A zero result therefore means every admitting box was contradictory, and the caller refuses.
    Only positions whose qatom_t carries a QATOM_STEREO_* flag are ever asked (the flag is set iff
    some box has a sign, which is also exactly when query_seal emits an SU_TETRA qstereo_t record --
    a geometry record's terminal carries no flag and needs no mask, its demand being in the record),
    so the early return costs one test on the ordinary stereo-query atom and its value is never read.
    """
    cdef qatom_t *qa = m.atoms + position
    cdef qbox_t *boxes = m.boxes + qa.box_begin
    cdef uint64_t *f = m.features + 4 * candidate
    cdef uint32_t b
    cdef uint8_t mask = 0
    if not (qa.flags & (QATOM_STEREO_CW | QATOM_STEREO_CCW)):
        return 0
    for b in range(qa.box_count):
        if not _box_admits(&boxes[b], m.anys, edge_word, f):
            continue
        if boxes[b].sign == QSIGN_BOTH:
            continue
        elif boxes[b].sign:
            mask |= boxes[b].sign
        else:
            mask |= <uint8_t> QSIGN_FREE
    return mask


cdef inline bint closures_admit(matcher_t *m, uint32_t position,
                                uint32_t candidate) noexcept nogil:
    """Do all closures owned by this position accept the molecule bond they require?

    Called immediately after atom_admits succeeds, before marking `used`.  Returns False if
    any closure bond is absent from the molecule or no box in its qbond_t admits the bond.

    Only neg[0] of each bond box is consulted -- see _fold_bond_into_atom (_query_seal.pxi)
    for why words 1-3 of a bond box are never used: compile_term runs box_fill_defaults over
    bond terms, which writes a not-a-radical (word 1) and neutral-charge (word 2) demand that
    no bond feature could ever satisfy.  The element half of the edge_word AND is harmlessly
    zero: a closure's boxes touch only the bond spans in word 0.  A future multi-hot bond
    primitive would need its own any-entry check here, because any entries are tested against
    the candidate's feature words -- not the half-edge word -- and are not consulted here.

    Invariant: closures[qa.closure_begin + c].to_index < position for all c.  query_seal
    assigns each closure to the HIGHER of its two endpoints (_query_seal.pxi); positions
    are filled in ascending order, so the target is always already mapped.  No MATCH_UNSET
    guard is needed and none is added -- it would be dead code.
    """
    cdef qatom_t *qa = m.atoms + position
    cdef halfedge_t *he
    cdef uint64_t edge_word
    cdef qbond_t *qb
    cdef uint32_t c, other, b
    cdef bint ok
    for c in range(qa.closure_count):
        other = m.mapping[m.closures[qa.closure_begin + c].to_index]
        he = csr_find_at(m.csr_begin, m.edges, candidate, other)
        if he == NULL:
            return False
        # csr_find_at returns a pointer into the CSR edge array; subtracting the base gives the
        # index, which is shared with the edge_words array (both are parallel to the CSR edges).
        edge_word = m.edge_words[he - m.edges]
        qb = m.bonds + m.closures[qa.closure_begin + c].bond_index
        ok = False
        for b in range(qb.box_count):
            if not (edge_word & m.bond_boxes[qb.box_begin + b].neg[0]):
                ok = True
                break
        if not ok:
            return False
    return True


cdef inline bint _stereo_frame_ok(uint8_t kind, uint8_t n_refs) noexcept nogil:
    """Can this target unit record be compared with a query's tetrahedral frame? (ruling F77)

    Two refusals, not two assertions.  A sign about an AXIS is not a sign about a CENTRE, so a
    matched atom that anchors SU_CIS_TRANS, SU_ALLENE or SU_ATROPISOMER refuses -- translate_stereo
    refuses the analogous case for the same reason.  And `n_refs` is checked rather than assumed:
    the gather below reads four slots, and a record that named fewer directions would have it
    reading padding as a direction.  Perception emits SU_TETRA with n_refs = 4 at its single emit
    site, so the second half is a guard over a value no molecule currently produces; it is pinned
    through _stereo_frame_probe rather than through a molecule.
    """
    return kind == SU_TETRA and n_refs == 4


cdef inline bint _scan_or_group(uint8_t *groups, uint32_t n) noexcept nogil:
    """Does any atom of the molecule carry an OR stereo group?

    One pass over one byte per atom, at match init and at every resume that re-borrows the arena
    (matcher_reseat), and only for a query that carries a stereo primitive.  Re-scanning on resume
    rather than caching across it is deliberate: a `for hit in q.get_mapping(mol):` body may edit the
    molecule, and the flag has to describe the arena the next candidate will be read from.
    """
    cdef uint32_t i
    for i in range(n):
        if sg_kind(groups[i]) == 2:
            return True
    return False


cdef enum:
    # What one stereo record says about the mapping in front of it.  A mask, not a verdict: an OR
    # group needs to know which of the two configurations the record would have accepted, because
    # the choice between them belongs to the group and not to the record (ruling F86).
    SF_REFUSE = 0            # no configuration of this centre satisfies the record
    SF_PLAIN = 1             # the target's stored configuration satisfies it
    SF_FLIPPED = 2           # the target's MIRROR configuration would satisfy it
    SF_NOTHING = 0x80        # the record states nothing here, so it constrains no group


cdef uint8_t _stereo_geometry_flips(matcher_t *m, qstereo_t *qs, uint8_t *group_out) noexcept nogil:
    """The SU_CIS_TRANS half of _stereo_record_flips: which states of one drawn geometry satisfy it.

    A `/` and `\\` pair states the geometry OUTRIGHT, so the demand is in the record (`spare`'s high
    byte) rather than in a box: there is no `m.sign_state` to read and no QSIGN_FREE disjunct to make
    the statement conditional.  What is shared with the tetrahedral half is everything after that --
    the target unit is read in the QUERY's frame, and the answer is a mask so that an E/Z mixture's
    group makes the choice (ruling F86) exactly as a racemate's does.
    """
    cdef stereo_unit_t *u = NULL
    cdef uint32_t anchor, other_term, k, want
    cdef uint32_t perm[4]
    cdef uint8_t parity, translated, demand, out

    group_out[0] = 0
    anchor = m.mapping[qs.position]
    other_term = m.mapping[qs.refs[2]]
    # The unit is anchored at ONE of the two terminals and which one is a fact about the target's slot
    # order, so both are tried.  A record naming an allene or an atropisomer here refuses: a sign about
    # an axis is not a statement about which side of a double bond a substituent is on.
    for k in range(m.unit_count):
        if (m.units[k].anchor == anchor or m.units[k].anchor == other_term) \
                and m.units[k].kind == SU_CIS_TRANS and m.units[k].n_refs == 4:
            u = m.units + k
            break
    if u is NULL:
        return SF_REFUSE
    # The unit's OWN anchor, not the query's: which terminal anchors is the target's fact, and the
    # parity is stored against that slot.
    parity = m.parities[u.anchor] if m.parities is not NULL else 0
    if not parity:
        return SF_REFUSE          # ruling F54 again: unconfigured is not "either geometry"
    for k in range(2):
        # the anchor's own pair leads `refs` (ruling F26), so which terminal anchors decides which of
        # the two marked substituents the frame starts with
        if (u.anchor == anchor) == (k == 0):
            want = m.mapping[qs.refs[0]]
        else:
            want = m.mapping[qs.refs[1]]
        if u.refs[2 * k] == want:
            perm[2 * k] = 2 * k
            perm[2 * k + 1] = 2 * k + 1
        elif u.refs[2 * k + 1] == want:
            perm[2 * k] = 2 * k + 1
            perm[2 * k + 1] = 2 * k
        else:
            # the query marks a substituent the unit does not carry on that terminal, so there is no
            # frame to read the geometry in
            return SF_REFUSE
    if m.stereo_groups is not NULL:
        group_out[0] = m.stereo_groups[u.anchor]
    translated = <uint8_t> translate_parity(parity, perm)
    demand = <uint8_t> (qs.spare >> 8)
    out = 0
    if translated == demand:
        out |= SF_PLAIN
    if (3 ^ translated) == demand:
        out |= SF_FLIPPED
    return out


cdef uint8_t _stereo_record_flips(matcher_t *m, qstereo_t *qs, uint8_t *group_out) noexcept nogil:
    """Which configurations of one query stereo record's anchor satisfy it, as an SF_* mask.

    Called once per record per candidate acceptance by stereo_admits, and again per record at a
    complete mapping by stereo_groups_admit -- the function is a pure read of `m.mapping`,
    `m.sign_state` and the arena, so the two callers agree by construction rather than by comment.

    `group_out` receives the anchor's SEG_STEREO_GROUPS byte, or 0 when the molecule has no such
    segment or the record states nothing; both callers need the group and neither should have to
    recover the anchor to get it.

    The comparison is made in the QUERY's frame: the record's `refs` are the anchor's query
    neighbours in the query's own ruling-F26 order, their images are looked up in the target unit's
    ref order, and translate_parity re-expresses the target's stored value under the permutation
    between the two.  The stored value comes from SEG_PARITY at the anchor's slot -- the unit
    record's parity field is always 0 and _stereo_emit is its only writer.

    WHICH sign is demanded is not read from the record (ruling F87): `qs.sign` is the union over the
    anchor's boxes and is reporting only.  The demand is `m.sign_state[qs.position]`, written by
    matcher_next when that position's candidate was accepted, because only there is the edge word
    live that says which of the anchor's boxes admitted it.  SF_FLIPPED is decided from the same
    mask, against the opposite parity (`3 ^ translated`, since the values are 1 and 2): a box that
    demands both signs accepts either choice, and one that demands neither refused above.
    """
    cdef stereo_unit_t *u
    cdef uint32_t i, j, k, anchor, want, refs_used, hslot, hcount
    cdef uint32_t perm[4]
    cdef uint8_t parity, sign_mask, translated, out
    cdef bint found

    group_out[0] = 0
    if (qs.spare & 0xFF) == SU_CIS_TRANS:
        # A GEOMETRY, not a centre: its own reading, and it must come before the sign_state load --
        # nothing wrote that byte, the record's demand being absolute and in the record.
        return _stereo_geometry_flips(m, qs, group_out)
    sign_mask = m.sign_state[qs.position]
    if sign_mask & <uint8_t> QSIGN_FREE:
        # Some box that admitted the anchor names no configuration, so this disjunct of the
        # query says nothing about stereo and nothing below may refuse on its behalf -- not the
        # frame accounting, not F77, not F67.  '[C@,N]' matching a nitrogen is this line.
        return SF_NOTHING
    if sign_mask == 0:
        # Every admitting box was self-contradictory ('[C;@;@@]').  Refused here rather than
        # pruned at compile time: ruling F87 keeps an unsatisfiable stereo query constructible.
        return SF_REFUSE
    # More than four query neighbours cannot be a tetrahedron, so no configuration satisfies the
    # record.  Refused at MATCH time, not at seal time: a query that cannot be satisfied is not a
    # construction error, and nothing in the stereo epic raises.
    if qs.n_refs > 4:
        return SF_REFUSE
    # FEWER THAN THREE NAMES NO FRAME, WHICH IS NOT THE SAME AS NAMING AN IMPOSSIBLE ONE.  Ruling F77
    # case 2 refused here, and the refusal made every such pattern match nothing at all -- `[C;@]`,
    # and the SMIRKS spelling `[C;@:1][Br;D1]` that says "a configured centre losing its bromide"
    # without wanting to enumerate the other three directions.  There is no permutation to read a
    # sign against, so the sign's VALUE is unenforceable and the honest reading is the part that is
    # still enforceable: the box already demands a configured centre, frame-free, via bit 7.  So this
    # widens to "configured, either sign", and `@` and `@@` are interchangeable in that position.
    if qs.n_refs < 3:
        return SF_NOTHING
    anchor = m.mapping[qs.position]
    u = NULL
    # A linear scan of the molecule's unit table, inside the DFS candidate loop.  The table is
    # one record per configured centre -- single digits for most molecules -- so an index would
    # be a per-target allocation to save a handful of loads.  What would justify one: a stereo
    # query with many stereo positions run against targets whose unit count is in the hundreds
    # (a peptide or a polysaccharide), where this becomes O(positions * units) per candidate.
    for k in range(m.unit_count):
        if m.units[k].anchor == anchor:
            u = m.units + k
            break
    if u is NULL:
        return SF_REFUSE          # the target names no stereo unit here, so it states nothing
    if not _stereo_frame_ok(u.kind, u.n_refs):
        return SF_REFUSE
    parity = m.parities[anchor] if m.parities is not NULL else 0
    if not parity:
        return SF_REFUSE          # ruling F54: no parity configured is not "no wedge drawn"
    # The permutation from the query's direction order to the unit's.  Its DIRECTION does not
    # matter -- a permutation and its inverse share a parity -- but its being a permutation of
    # 0..3 does, which is what the total, injective matching below establishes.
    refs_used = 0
    for i in range(4):
        perm[i] = 0
    for i in range(4):
        if i < qs.n_refs:
            want = m.mapping[qs.refs[i]]
            found = False
            for j in range(4):
                if refs_used & (<uint32_t> 1 << j):
                    continue
                if u.refs[j] == want:
                    perm[i] = j
                    refs_used |= <uint32_t> 1 << j
                    found = True
                    break
            if not found:
                # A named query direction whose image is not one of the target's four: the query
                # describes a neighbour the unit does not have, so there is no permutation to
                # take a parity through.  `want` cannot be MATCH_UNSET here (readiness is the
                # last of the refs to be mapped), and since Minor 4 the two sentinels differ
                # anyway, so an unmapped ref could not be mistaken for the target's unnamed slot.
                return SF_REFUSE
        else:
            # The query's UNNAMED direction, and ruling F88: it may pair with a hydrogen,
            # whether the target drew that hydrogen or left it implicit.  Whether a hydrogen is
            # drawn is an input-representation choice and stereo_units' reference order is built
            # so it does not move the tuple; matching must not read it either.
            #
            # The pairing is by hydrogen-ness of the whole ref list rather than by "take the one
            # slot left over", and that is deliberate: a centre with TWO hydrogen directions has
            # two directions the query cannot tell apart, and pairing with either would answer
            # from perception's arbitrary tie-break between them.  Perception does NOT refuse
            # such a centre a unit -- measured, `C([H])([H])(F)Cl` with a stored parity emits
            # SU_TETRA with refs (F, Cl, H, H) -- and `stereogenic` cannot be consulted here
            # (it is the marked door's product and ruling F76 sends the kernel through the
            # unmarked one), so the refusal has to be made from the frame itself.
            #
            # RULING F67 survives as the hcount == 0 arm: C(F)(Cl)(Br)I against a query naming
            # three neighbours pairs its unnamed direction with nothing, rather than with the
            # iodine.  Measured, with that refusal removed: `@` matches the parity-1 target and
            # `@@` the parity-2 one, so both signs become satisfiable and the query gets a
            # confident answer about a molecule it never described.  On success refs_used is
            # 0xF: every direction accounted for.
            #
            # What is counted is every direction that is not a named heavy atom: SU_NO_REF (an
            # implicit hydrogen or a lone pair) and any explicit H neighbour.
            hcount = 0
            hslot = 4
            for j in range(4):
                # SU_NO_REF first: it is not an atom index and must not be dereferenced.
                if u.refs[j] == SU_NO_REF or m.satoms[u.refs[j]].element == 1:
                    hcount += 1
                    hslot = j
            if hcount != 1 or refs_used & (<uint32_t> 1 << hslot):
                return SF_REFUSE
            perm[i] = hslot
            refs_used |= <uint32_t> 1 << hslot
    # Among directions that are NOT hydrogens the flavour is still not distinguished: a
    # three-neighbour query matches a centre whose fourth direction is a lone pair the same way
    # it matches one whose fourth is a hydrogen, unless the query says `h`.  The `h` primitive is
    # an ordinary box screen and decides that case on its own, before this function runs at all.
    if m.stereo_groups is not NULL:
        group_out[0] = m.stereo_groups[anchor]
    # sign_mask holds QSIGN_CW, QSIGN_CCW or both by here (QSIGN_FREE and 0 returned above), and
    # translate_parity returns 1 or 2, so this is a bit test.  Both bits set means two different
    # boxes admitted the candidate demanding opposite configurations, which either one satisfies.
    translated = <uint8_t> translate_parity(parity, perm)
    out = 0
    if sign_mask & translated:
        out |= SF_PLAIN
    if sign_mask & (3 ^ translated):
        out |= SF_FLIPPED
    return out


cdef bint stereo_admits(matcher_t *m, uint32_t position) noexcept nogil:
    """Do the stereo primitives whose frame is complete at `position` accept the mapping so far?

    Called after m.mapping[position] is written, for every position, and short-circuited by
    `m.stereo_count` at the call site so an ordinary query pays one test of a register.  A record's
    `readiness` is the position at which the last of its directions becomes mapped, computed at seal
    time (query_seal, QSEG_STEREO), which is why there is no "am I ready yet" test here.

    This is the per-unit half of the decision, and the group kind decides how much of it is made
    here:

    * unspecified (0) and ABS (1): the stored configuration is the only one the target has, so a
      record that does not accept it refuses the candidate now.
    * AND (3): the racemate is present, so either configuration is there to be matched and only
      SF_REFUSE -- the frame accounting -- can refuse.  F67 is about whether the query describes the
      target's frame, and an AND group does not make an under-specified frame acceptable.
    * OR (2): the sign is the GROUP's to choose, not this centre's, so nothing is decided here beyond
      SF_REFUSE; stereo_groups_admit makes the choice at the complete mapping (ruling F86).  This arm
      still prunes: a record whose frame does not fit refuses the candidate here, inside the DFS, and
      only the question of which configuration was taken is deferred.
    """
    cdef uint32_t r
    # group_byte is initialised here only to satisfy the control-flow analysis: it is an out-parameter
    # of _stereo_record_flips, which writes it before any early return, and Cython cannot see that.
    cdef uint8_t flips, kind, group_byte = 0

    for r in range(m.stereo_count):
        if m.stereo[r].readiness != position:
            continue
        flips = _stereo_record_flips(m, m.stereo + r, &group_byte)
        if flips == SF_REFUSE:
            return False
        if flips & SF_NOTHING:
            continue
        kind = sg_kind(group_byte)
        if kind == 2 or kind == 3:
            continue
        if not (flips & SF_PLAIN):
            return False
    return True


cdef bint stereo_groups_admit(matcher_t *m) noexcept nogil:
    """Can every OR group be given ONE configuration that satisfies all of its matched units?

    Ruling F86: the OR decision is made here, at a complete mapping, rather than as a per-group
    variable in the DFS frame.  Every record is re-read -- _stereo_record_flips is a pure function of
    the mapping, so the answers are the ones stereo_admits already saw -- and each group accumulates
    the intersection of its members' acceptable configurations.  A group survives if that intersection
    is non-empty; the mapping survives if every group it touched does.

    Cost: O(records) per complete mapping, against O(2 ** groups) for the alternative of searching
    the assignments in the frame.  The two are not bounds on the same quantity -- the assignment
    search multiplies the DFS, this multiplies the solutions -- and the measurement is in the task
    report.

    `seen` is the set of group numbers met, and is load-bearing: `plain` and `flipped` start at zero,
    so without it an untouched group would be indistinguishable from a group whose intersection came
    out empty.  Group numbers are 1..63 (set_stereo_group's range) and 0 is what ABS and unspecified
    carry, so one uint64_t covers every group a molecule can have and bit 0 is never an OR member's.

    A record that states nothing (SF_NOTHING) constrains no group, including one whose other members
    do constrain it: '[C@,N]' matching the nitrogen is a disjunct that made no claim, and a group must
    not be cornered by a claim that was not made.  A record that refuses outright cannot appear here
    -- stereo_admits refused that candidate at the record's readiness position and the mapping has not
    changed since -- but is handled anyway, and in the safe direction: it clears both of its group's
    bits, so the group, and with it the mapping, fails.
    """
    cdef uint64_t seen = 0, plain = 0, flipped = 0, bit
    cdef uint32_t r
    cdef uint8_t flips, group_byte = 0      # written by the callee; see stereo_admits

    for r in range(m.stereo_count):
        flips = _stereo_record_flips(m, m.stereo + r, &group_byte)
        if flips & SF_NOTHING:
            continue
        if sg_kind(group_byte) != 2:
            continue
        bit = <uint64_t> 1 << sg_group(group_byte)
        if not (seen & bit):
            seen |= bit
            plain |= bit
            flipped |= bit
        if not (flips & SF_PLAIN):
            plain &= ~bit
        if not (flips & SF_FLIPPED):
            flipped &= ~bit
    return not (seen & ~(plain | flipped))


cdef inline bint group_admits(matcher_t *m, uint32_t position, uint32_t candidate) noexcept nogil:
    """Does placing `candidate` at `position` respect the component-group constraints?

    This function is side-effect-free: it only reads m.group_anchor and never writes it.
    The anchor is set (and cleared on backtrack) by _set_group_anchor/_maybe_clear_anchor in
    matcher_next, after all checks for a position have passed.

    group_count == 0 is the fast path: ungrouped queries never call this function (the caller
    guards with `if m.group_count`), but the guard below also makes it safe if called directly,
    and prevents a NULL dereference on component_of_position.
    """
    if m.group_count == 0:
        return True
    cdef qcomp_t *comps = m.components
    cdef uint32_t comp = m.component_of_position[position]
    cdef int32_t group = comps[comp].group
    cdef uint32_t label, g
    if group < 0 or comps[comp].begin != position:
        return True                     # unconstrained, or not this component's root
    label = m.labels[candidate]
    if m.group_anchor[group] >= 0:
        return <uint32_t> m.group_anchor[group] == label
    for g in range(m.group_count):
        if g != <uint32_t> group and m.group_anchor[g] == <int32_t> label:
            return False                # another group already owns this molecule component
    return True


cdef inline void _set_group_anchor(matcher_t *m, uint32_t position, uint32_t candidate) noexcept nogil:
    """If position is an unclaimed group root, record which molecule component it landed in."""
    cdef qcomp_t *comps = m.components
    cdef uint32_t comp = m.component_of_position[position]
    cdef int32_t group = comps[comp].group
    if group < 0 or comps[comp].begin != position:
        return
    if m.group_anchor[group] < 0:
        m.group_anchor[group] = <int32_t> m.labels[candidate]


cdef inline void _maybe_clear_anchor(matcher_t *m, uint32_t position) noexcept nogil:
    """On backtrack past position: clear group_anchor if position is the lowest-position root of
    its group.

    The DFS invariant guarantees that when this function runs, all positions 0..position are
    still mapped (backtrack unwinds in reverse order, so nothing below position has been cleared
    yet).  Therefore every `begin < position` root IS currently mapped -- the `!= MATCH_UNSET`
    test would be vacuously true and is not written.  The scan simply asks: is there any component
    root with the same group number and a lower begin?  If yes, that root set the anchor and will
    clear it when the DFS reaches it; do not touch it.  If no, this position IS the lowest root,
    it set the anchor, and it must clear it.
    """
    cdef qcomp_t *comps = m.components
    cdef uint32_t comp = m.component_of_position[position]
    cdef int32_t group = comps[comp].group
    cdef uint32_t ci
    if group < 0 or comps[comp].begin != position:
        return
    # Clear only when this is the lowest-position root for this group.
    for ci in range(m.component_count):
        if m.components[ci].group == group and m.components[ci].begin < position:
            return  # a lower root exists; it set the anchor
    m.group_anchor[group] = -1


cdef inline uint32_t _candidate_begin(matcher_t *m, uint32_t position) noexcept nogil:
    """The first cursor value for a position.  For CAND_NEIGHBOURS the parent must be mapped."""
    if m.candidate_kind[position] == CAND_ELEMENT_BUCKET:
        return m.bucket_begin[position]
    if m.candidate_kind[position] == CAND_ALL_ATOMS:
        return 0
    return m.csr_begin[m.mapping[m.atoms[position].back]]


cdef inline uint32_t _candidate_end(matcher_t *m, uint32_t position) noexcept nogil:
    """One past the last cursor value for a position."""
    if m.candidate_kind[position] == CAND_ELEMENT_BUCKET:
        return m.bucket_end[position]
    if m.candidate_kind[position] == CAND_ALL_ATOMS:
        return m.structure_atom_count
    return m.csr_begin[m.mapping[m.atoms[position].back] + 1]


cdef bint query_may_match(Query query, Structure structure) noexcept nogil:
    """A sound lower bound: False means no embedding exists, True means maybe.

    Checks, in increasing cost order:
      1. Atom count: a query with more atoms than the target cannot embed.
      2. Element demand: QSEG_DEMAND_LIST holds (element, count) pairs built at seal from
         the QSEG_ELEMENT_DEMAND histogram.  Each pair names an element whose every atom in the
         query demands exactly that element.  The molecule must have at least count atoms of that
         element in its element index; if not, no embedding exists.  Iterating the compact list
         (O(demanded elements), typically 1-3 entries) is strictly cheaper than scanning the full
         118-slot histogram used before Task 14 fix 1.
      3. Signature: every bit in the query's four demand words must be present in the
         molecule's union row (structure_features, word 0).  A query bit is placed there
         only when every box of some atom requires it, and a required bit implies the
         matching molecule atom actually has it -- so the union row carries it too.
         All four words are used here (unlike sig_contains, which masks out ring and
         hybridization bits that change under embedding): a query demand bit is a property
         the match atom genuinely has, not a context-sensitive descriptor.
    """
    cdef uint32_t *demand_list = query_demand_list(query)
    cdef uint64_t *union_words = structure_features(structure)
    cdef uint32_t dl_len = query.header.segments[QSEG_DEMAND_LIST].length // (2 * sizeof(uint32_t))
    cdef uint32_t dl_pos, elem, cnt, w_idx
    if query.header.atom_count > structure.header.atom_count:
        return False
    for dl_pos in range(dl_len):
        elem = demand_list[dl_pos * 2]
        cnt  = demand_list[dl_pos * 2 + 1]
        if cnt > element_bucket_end(structure, elem) - element_bucket_begin(structure, elem):
            return False
    for w_idx in range(4):
        if (union_words[w_idx] & query.header.signature[w_idx]) != query.header.signature[w_idx]:
            return False
    return True


cdef void matcher_free(matcher_t *m) noexcept nogil:
    """Release the search state.  Safe on a zeroed struct, so an early-out init needs no undo.

    Also marks the matcher exhausted.  Without that store `m.depth` keeps whatever the search left
    behind, and a matcher_next on a freed struct would index candidate_kind through NULL instead of
    returning False -- no current caller does that, but "safe on a zeroed struct" invites the
    broader reading, and one store makes it true for Tasks 11-14 as well.
    """
    free(m.mapping)
    free(m.used)
    free(m.candidate)
    free(m.candidate_kind)
    free(m.bucket_begin)
    free(m.bucket_end)
    free(m.group_anchor)
    free(m.component_of_position)
    free(m.labels)
    free(m.sign_state)
    m.sign_state = NULL
    m.mapping = NULL
    m.used = NULL
    m.candidate = NULL
    m.candidate_kind = NULL
    m.bucket_begin = NULL
    m.bucket_end = NULL
    m.group_anchor = NULL
    m.component_of_position = NULL
    m.labels = NULL
    m.depth = MATCH_DONE


cdef int matcher_init(matcher_t *m, Query query, Structure structure,
                      bint automorphism_filter=False) except -1:
    """Resolve every pointer, decide each position's candidate source, and arm depth 0.

    Stores no reference to either object: the caller must keep both alive for the whole search.
    A matcher that cannot possibly match -- an empty target, or a query with more atoms than the
    target has -- is armed as already exhausted, with nothing allocated; matcher_next returns
    False and matcher_free is still safe to call.

    `automorphism_filter` asks for one embedding per automorphism orbit.  The group is already in
    the arena -- query_seal computed it -- so this costs nothing at init beyond two loads.
    """
    cdef qatom_t *atoms
    cdef qbox_t *boxes
    cdef uint32_t n, sn, i, e, gc

    cdef int element

    memset(m, 0, sizeof(matcher_t))
    m.depth = MATCH_DONE          # every early return below leaves an exhausted matcher

    n = query.header.atom_count
    sn = structure.header.atom_count
    m.atom_count = n
    m.structure_atom_count = sn
    # No caller can reach the n == 0 arm: sealed() raises ValueError('an empty query matches
    # nothing') and every entry point seals before calling this.  Kept so the kernel is total on
    # its own terms rather than on its callers' -- do not go hunting for the caller that needs it.
    if n == 0 or sn == 0 or n > sn:
        return 0

    # Pre-search screen: reject without allocating when the molecule provably cannot hold the
    # query.  This is a performance gate only -- removing it changes no count and no result
    # (ruling 5 / Task 14).  It must sit before any allocation so that matcher_free on the
    # zero-initialised struct is still safe on a False return.
    if not query_may_match(query, structure):
        return 0

    # THE UNMARKED DOOR (ruling F76), and BEFORE the first arena pointer (ruling F60).
    # Unmarked because a stereo primitive asks whether the target STATES the configuration the
    # query names; whether that statement is justified is validate_stereo's question, not the
    # kernel's.  Routing matching through ensure_stereo_units would make every stereo query pay a
    # stereogenicity witness search per target -- the cost ruling F70 removed from the journal
    # apply -- and would drag ruling F62's truncation policy into the kernel; through this door
    # truncation cannot affect a match at all.  Built once per match attempt, never in the DFS,
    # and only for a query that carries a stereo primitive: this call can realloc the arena, so
    # every pointer below is taken after it.
    if query.header.flags & QFLAG_HAS_STEREO:
        ensure_stereo_units_unmarked(structure)
        m.stereo = query_stereo(query)
        m.stereo_count = query.header.stereo_count
        m.units = structure_stereo_units(structure)
        m.unit_count = structure_stereo_unit_count(structure)
        if structure_has(structure, SEG_PARITY):
            m.parities = structure_parities(structure)
        if structure_has(structure, SEG_STEREO_GROUPS):
            m.stereo_groups = structure_stereo_groups(structure)
            m.has_or_group = _scan_or_group(m.stereo_groups, structure.header.atom_count)

    m.atoms = query.atoms()
    m.satoms = structure.atoms()
    m.boxes = query_boxes(query)
    m.anys = query_any(query)
    m.closures = query_closures(query)
    m.bonds = query_bonds(query)
    m.bond_boxes = query_bond_boxes(query)
    m.components = query_components(query)
    m.component_count = query.header.component_count
    m.automorphisms = query_automorphisms(query)
    m.automorphism_count = query.header.automorphism_count
    m.automorphism_filter = automorphism_filter
    # fill_features writes a union row first, so atom i's row starts at 4 + 4 * i; offsetting the
    # base by 4 makes m.features + 4 * i atom i's row.
    m.features = structure_features(structure) + 4
    m.edge_words = structure_edge_words(structure)
    m.csr_begin = csr_ptr(structure)
    m.edges = csr_edges(structure)
    m.element_slots = structure_element_index(structure) + 120

    # Determine group_count: one past the highest group number, or 0 if no groups.
    gc = 0
    if query.header.flags & QFLAG_HAS_GROUP:
        for i in range(m.component_count):
            if m.components[i].group >= 0:
                e = <uint32_t> m.components[i].group + 1
                if e > gc:
                    gc = e
    m.group_count = gc

    # libc rather than PyMem so that matcher_free is honestly nogil: a caller running the whole
    # search with the GIL released still has to release the state at the end.
    m.mapping = <uint32_t *> malloc(n * sizeof(uint32_t))
    m.candidate = <uint32_t *> malloc(n * sizeof(uint32_t))
    m.candidate_kind = <uint8_t *> malloc(n)
    m.bucket_begin = <uint32_t *> malloc(n * sizeof(uint32_t))
    m.bucket_end = <uint32_t *> malloc(n * sizeof(uint32_t))
    m.used = <uint8_t *> malloc(sn)
    if m.stereo_count:
        # Ruling F87's per-position sign demand.  Allocated only for a stereo query, like every
        # other stereo cost in this struct; left NULL otherwise, and nothing reads it then because
        # every read is behind the same `m.stereo_count` guard as this allocation.
        m.sign_state = <uint8_t *> malloc(n)
    if gc > 0:
        m.group_anchor = <int32_t *> malloc(gc * sizeof(int32_t))
        m.component_of_position = <uint32_t *> malloc(n * sizeof(uint32_t))
        m.labels = <uint32_t *> malloc(sn * sizeof(uint32_t))
    if (m.mapping is NULL or m.candidate is NULL or m.candidate_kind is NULL or
            m.bucket_begin is NULL or m.bucket_end is NULL or m.used is NULL or
            (m.stereo_count and m.sign_state is NULL) or
            (gc > 0 and (m.group_anchor is NULL or m.component_of_position is NULL or
                         m.labels is NULL))):
        matcher_free(m)
        raise MemoryError('matcher allocation failed')

    memset(m.used, 0, sn)
    if m.stereo_count:
        # Every entry is written at the position's accept before it is read, so this zeroing is
        # hygiene rather than correctness -- but a zero reads as "no admitting box accepts any
        # configuration", which refuses, and that is the safe direction for a bug to fail in.
        memset(m.sign_state, 0, n)

    if gc > 0:
        # Initialise all anchors to unclaimed and fill component_of_position from the ranges.
        for i in range(gc):
            m.group_anchor[i] = -1
        for i in range(m.component_count):
            for e in range(m.components[i].begin, m.components[i].end):
                m.component_of_position[e] = i
        # label_components uses a local DFS stack and never appends to the arena, so the
        # structure pointers cached above remain valid throughout the search.
        if label_components(structure, m.labels) < 0:
            matcher_free(m)
            raise MemoryError('component label allocation failed')

    atoms = m.atoms
    boxes = m.boxes
    for i in range(n):
        m.mapping[i] = MATCH_UNSET
        m.candidate[i] = 0
        m.bucket_begin[i] = 0
        m.bucket_end[i] = 0
        if not atoms[i].flags & QATOM_ROOT:
            m.candidate_kind[i] = CAND_NEIGHBOURS
            continue
        element = _root_element(boxes + atoms[i].box_begin, atoms[i].box_count)
        if element < 0:
            m.candidate_kind[i] = CAND_ALL_ATOMS
        else:
            m.candidate_kind[i] = CAND_ELEMENT_BUCKET
            e = <uint32_t> element
            m.bucket_begin[i] = element_bucket_begin(structure, e)
            m.bucket_end[i] = element_bucket_end(structure, e)

    m.depth = 0
    m.candidate[0] = _candidate_begin(m, 0)
    return 0


cdef inline void matcher_reseat(matcher_t *m, Structure structure) noexcept nogil:
    """Re-borrow the arena pointers, for a search that has let Python run.

    The struct comment says no arena append may happen during a search, and a search that yields
    to Python cannot enforce that: the body of a `for hit in q.get_mapping(mol):` loop may call
    stereo_units() or component_labels() on the same molecule, and structure_append reallocates
    through PyMem_Realloc.  Rather than forbid that from a docstring, the two generators re-seat
    on every resume, which removes the precondition instead of documenting it.

    Only pointers into the MOLECULE arena are listed.  The query arena is immutable after seal
    (there is no query_append), and everything else in the struct is malloc'd and owned.  A
    caller running the whole search under nogil, like count(), appends nothing and needs none of
    this.  Every segment read here is built eagerly by rebuild_derived, so re-seating can only
    ever produce the same addresses or the moved ones -- never a segment that was absent at init
    and is present now, which would change the search mid-flight.

    The stereo pointers are the one group that is NOT built by rebuild_derived: SEG_STEREO_UNIT is
    lazy.  matcher_init builds it eagerly for exactly this reason, so it is present at every resume
    and re-seating can only find it moved -- but the count is re-read with the pointer, because a
    resume that rebuilt the table would otherwise pair a new base with a stale count.  The query
    side is left alone: there is no query_append.

    SEG_PARITY is persistent, laid out once, so unlike the lazy SEG_STEREO_UNIT it can only ever
    be found moved and never found newly present -- a NULL m.parities stays NULL for the whole
    search.
    """
    m.features = structure_features(structure) + 4
    m.edge_words = structure_edge_words(structure)
    m.csr_begin = csr_ptr(structure)
    m.edges = csr_edges(structure)
    m.element_slots = structure_element_index(structure) + 120
    m.satoms = structure.atoms()
    if m.stereo_count:
        m.units = structure_stereo_units(structure)
        m.unit_count = structure_stereo_unit_count(structure)
        m.parities = NULL
        if structure_has(structure, SEG_PARITY):
            m.parities = structure_parities(structure)
        m.stereo_groups = NULL
        m.has_or_group = False
        if structure_has(structure, SEG_STEREO_GROUPS):
            m.stereo_groups = structure_stereo_groups(structure)
            m.has_or_group = _scan_or_group(m.stereo_groups, structure.header.atom_count)


cdef inline bint mapping_is_canonical(matcher_t *m) noexcept nogil:
    """Is this embedding the lexicographically smallest member of its automorphism orbit?

    For each stored automorphism sigma, compare the permuted mapping against the current one
    position by position.  The first position where they differ decides: if the permuted image is
    smaller, some other embedding of the same orbit is smaller than this one and will be (or was)
    reported in its place, so this one is a duplicate.  Equal all the way through means sigma fixes
    this embedding, which is not a reason to reject it.

    Comparing only against the STORED rows -- a subset of the group when QFLAG_PARTIAL_AUTOMORPHISM
    is set -- keeps this safe in the only direction that matters.  A smaller subset admits more
    embeddings, never fewer: the true orbit minimum is still minimal against a subset, so it is
    never rejected, and the filter degrades to reporting some duplicates rather than losing hits.
    """
    cdef uint32_t row, i, n = m.atom_count
    cdef uint32_t *sigma
    for row in range(m.automorphism_count):
        sigma = m.automorphisms + row * n
        for i in range(n):
            if m.mapping[sigma[i]] < m.mapping[i]:
                return False
            elif m.mapping[sigma[i]] > m.mapping[i]:
                break
    return True


cdef bint matcher_next(matcher_t *m) noexcept nogil:
    """Fill m.mapping with the next embedding and return True, or return False when exhausted."""
    cdef uint32_t position, cursor, stop, c
    cdef uint64_t edge_word
    cdef uint8_t kind
    cdef bint filled

    if m.depth == MATCH_DONE:
        return False
    if m.depth == m.atom_count:
        # Resuming from the embedding the previous call reported: give back its last position and
        # carry on from that position's cursor, which already points past the atom it used.
        m.depth -= 1
        if m.group_count:
            _maybe_clear_anchor(m, m.depth)
        m.used[m.mapping[m.depth]] = 0
        m.mapping[m.depth] = MATCH_UNSET

    while True:
        position = m.depth
        kind = m.candidate_kind[position]
        cursor = m.candidate[position]
        stop = _candidate_end(m, position)
        filled = False
        while cursor < stop:
            if kind == CAND_ELEMENT_BUCKET:
                c = m.element_slots[cursor]
                edge_word = m.features[4 * c]
            elif kind == CAND_ALL_ATOMS:
                c = cursor
                edge_word = m.features[4 * c]
            else:
                c = m.edges[cursor].to
                edge_word = m.edge_words[cursor]
            cursor += 1
            # one byte load before four: injectivity is the cheaper of the two tests
            if m.used[c]:
                continue
            if atom_admits(m, edge_word, position, c):
                if not closures_admit(m, position, c):
                    continue
                if m.group_count and not group_admits(m, position, c):
                    continue
                # Stereo reads mapping[position], so the write comes first and is taken back if the
                # check refuses; mapping is scratch until `filled` is set, so nothing else can see
                # it.  The group anchor is set only afterwards -- it is the one side effect here,
                # and a refused candidate must not leave it claimed.
                m.mapping[position] = c
                if m.stereo_count:
                    # The sign demand has to be recorded HERE, not in stereo_admits: it depends on
                    # which of this position's boxes admitted `c`, and `edge_word` -- the half-edge
                    # word that decides that -- is live only at this call site (ruling F87).  A
                    # record whose frame closes at a later position reads it back then.
                    m.sign_state[position] = _stereo_sign_mask(m, edge_word, position, c)
                    if not stereo_admits(m, position):
                        m.mapping[position] = MATCH_UNSET
                        continue
                if m.group_count:
                    _set_group_anchor(m, position, c)
                m.used[c] = 1
                filled = True
                break
        m.candidate[position] = cursor
        if filled:
            m.depth += 1
            if m.depth == m.atom_count:
                # The OR group decision (ruling F86), before the automorphism filter and before the
                # yield: a mapping the groups refuse is not an embedding, so it must not be counted,
                # returned, or offered to mapping_is_canonical as a representative.  `has_or_group`
                # is false for every target without an OR group, which is nearly all of them.
                if ((not m.has_or_group or stereo_groups_admit(m))
                        and (not m.automorphism_filter or mapping_is_canonical(m))):
                    return True
                # An automorphic duplicate, or a group demand no single choice per group can meet.
                # Give the last position back and keep searching, which is exactly what the
                # resume-from-embedding path at the top of this function does -- the difference is
                # only that the caller never saw this one.
                m.depth -= 1
                if m.group_count:
                    _maybe_clear_anchor(m, m.depth)
                m.used[m.mapping[m.depth]] = 0
                m.mapping[m.depth] = MATCH_UNSET
                continue
            # The child's source may depend on this position's atom (CAND_NEIGHBOURS), so its
            # cursor is armed here rather than at init.
            m.candidate[m.depth] = _candidate_begin(m, m.depth)
        else:
            if position == 0:
                m.depth = MATCH_DONE
                return False
            m.depth -= 1
            if m.group_count:
                _maybe_clear_anchor(m, m.depth)
            m.used[m.mapping[m.depth]] = 0
            m.mapping[m.depth] = MATCH_UNSET


def _stereo_frame_probe(int kind, int n_refs):
    """`_stereo_frame_ok` on a forged (kind, n_refs) pair, so its two refusals can be pinned.

    Ruling F77 cases 1 and 3.  Case 1 is reachable from a molecule -- an sp2 carbon with three
    neighbours anchors SU_CIS_TRANS and a stereo query can match it -- and is pinned that way too.
    Case 3 is not: `_perceive_stereo_units` has one SU_TETRA emit site and it passes the literal 4,
    so no molecule can produce a tetrahedral record with another value.  Per ruling F75 this probe
    pins the guard's behaviour and says nothing about whether the guard could be relaxed.
    """
    return _stereo_frame_ok(<uint8_t> kind, <uint8_t> n_refs)
