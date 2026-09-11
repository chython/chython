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
# The journal and the seal: query_seal and the automorphism group.
#
# ---------------------------------------------------------------------------
# The journal and the seal — Task 8
# ---------------------------------------------------------------------------
# A query is built by appending operations to a journal, exactly as MoleculeContainer does,
# and sealed once into an arena.  journal_t lives in _molecule_container.pxi, which is
# included after this layer, so the query gets its own record type: the molecule journal
# records values, this one records tokens.

cdef enum:
    QOP_ADD_ATOM = 1        # a = stable id
    QOP_ATOM_TOKEN = 2      # a = stable id, then opcode/kind/value/negated
    QOP_ADD_BOND = 3        # a, b = stable ids
    QOP_BOND_TOKEN = 4      # a, b = stable ids, then opcode/kind/value/negated
    QOP_SET_GROUP = 5       # a = stable id, v = group number (-1 = ungrouped)
    QOP_SET_MASKED = 6      # a = stable id
    QOP_SET_MAP = 7         # a = stable id, v = map number
    QOP_SET_STEREO_GROUP = 8  # a = stable id, kind = STEREO_OR/STEREO_AND, v = group number
    QOP_SET_BOND_DIRECTION = 9  # a, b = stable ids, v = SMI_DIR_UP/SMI_DIR_DOWN, oriented a -> b


cdef struct qop_t:          # 28 bytes
    uint32_t op
    uint32_t a
    uint32_t b
    uint32_t opcode
    uint32_t kind
    int32_t value
    uint32_t negated


cdef uint32_t Q_NO_SLOT = 0xFFFFFFFF    # a slot that does not exist


cdef inline void _token_from_op(qtoken_t *dst, qop_t *src) noexcept nogil:
    """Copy one journal token's payload into the box compiler's scratch.

    Atom tokens and bond tokens are gathered by two separate loops in `query_seal` that differ
    only in which op code and which slot they select.  The copy itself must not differ, or a
    primitive would compile one way on an atom and another way on a bond.
    """
    dst.opcode = src.opcode
    dst.kind = src.kind
    dst.value = src.value
    dst.negated = src.negated != 0


cdef inline size_t _alloc_at_least(uint32_t want) noexcept nogil:
    """A malloc count of at least one: PyMem_Malloc(0) may legally return NULL."""
    if want:
        return want
    return 1


cdef inline uint32_t _term_any_total(wterm_t *term) noexcept nogil:
    """How many QSEG_ANY records one term's boxes will occupy.

    The sizing pass must agree with `_emit_boxes` about this number or the arena is short; both
    the atom half and the bond half of the sizing pass ask the same question.
    """
    cdef uint32_t i, want = 0
    for i in range(term.count):
        want += term.boxes[i].any_count
    return want


cdef inline uint32_t _emit_boxes(qbox_t *out_boxes, uint32_t box_cursor, qany_t *out_any,
                                uint32_t *any_cursor, wterm_t *term) noexcept nogil:
    """Write one compiled term's boxes into the arena; returns the advanced box cursor.

    `query_seal` emits boxes twice, once for atoms and once for bonds, into two different
    segments -- and the records are identical, down to clearing `spare`.  Written out by hand
    the two copies were twelve lines each of `out_boxes[box_cursor].neg[0] = ...`, which is both
    how the duplication got there and why it would have drifted.  `any_cursor` is by pointer
    because atoms and bonds draw their any-words from one shared QSEG_ANY run.

    `sign` (ruling F87) is copied like any other field; a bond term can never carry one, because
    PRIM_STEREO is not a BPRIM_* and no bond token stream can reach that branch of prim_apply.
    """
    cdef uint32_t k, u
    cdef qbox_t *box
    cdef qany_t *any
    cdef wbox_t *src
    for k in range(term.count):
        box = &out_boxes[box_cursor]
        src = &term.boxes[k]
        for u in range(4):
            box.neg[u] = src.neg[u]
        box.any_begin = any_cursor[0]
        box.any_count = <uint16_t> src.any_count
        box.sign = src.sign
        box.spare = 0
        for u in range(src.any_count):
            any = &out_any[any_cursor[0]]
            any.mask = src.any_mask[u]
            any.word = src.any_word[u]
            any.spare = 0
            any_cursor[0] += 1
        box_cursor += 1
    return box_cursor


cdef uint32_t _find_bond_slot(uint32_t *bond_a, uint32_t *bond_b, uint32_t bond_count,
                              uint32_t *slot_of, uint32_t next_id,
                              uint32_t n, uint32_t m) noexcept nogil:
    """The bond slot for a pair of stable ids, or Q_NO_SLOT.

    A linear scan rather than an open-addressed table keyed by (min, max): a query with more
    than a few dozen bonds does not exist in practice, and seal runs once per query.
    """
    cdef uint32_t k, s, t
    if n >= next_id or m >= next_id:
        return Q_NO_SLOT
    s = slot_of[n]
    t = slot_of[m]
    if s == Q_NO_SLOT or t == Q_NO_SLOT:
        return Q_NO_SLOT
    for k in range(bond_count):
        if (bond_a[k] == s and bond_b[k] == t) or (bond_a[k] == t and bond_b[k] == s):
            return k
    return Q_NO_SLOT


cdef inline uint32_t _single_bit_index(uint64_t w) noexcept nogil:
    """The index of the only set bit in w.  Seal-time only, so a loop is fast enough."""
    cdef uint32_t i
    for i in range(64):
        if w >> i & <uint64_t> 1:
            return i
    return 64


cdef uint32_t _term_element_card(wterm_t *term) noexcept nogil:
    """How many distinct elements at least one box of this disjunction allows.

    Q_NO_SLOT (0xFFFFFFFF) when any box leaves the element span untouched -- an atom that
    accepts every element must never win a component root, because the DFS would then have to
    seed from every atom of the target instead of from one element bucket.
    """
    cdef uint32_t i, total
    cdef uint64_t light = 0, heavy = 0
    cdef bint heavy_open = False
    cdef wbox_t *box

    for i in range(term.count):
        box = &term.boxes[i]
        if not box.touched[0] & <uint64_t> W0_ELEMENT_SPAN:
            return Q_NO_SLOT
        light |= ~box.neg[0] & <uint64_t> W0_LIGHT_ELEMENT_SPAN
        # Word 0 bit 0 is the heavy-element marker shared by every element above 56: a box that
        # forbids it allows no heavy element at all, whatever its word-1 bits say.
        if not box.neg[0] & <uint64_t> 1:
            heavy_open = True
            heavy |= ~box.neg[1] & <uint64_t> W1_ELEMENT_SPAN
    total = _popcount64(light)
    if heavy_open:
        total += _popcount64(heavy)
    return total


cdef inline int _element_from_neg(uint64_t neg0, uint64_t neg1) noexcept nogil:
    """The single element one box's forbidden masks leave open, or -1 when it is not exactly one.

    The `57 - bit` / `57 + bit` decoding was starting to multiply -- _term_exact_element and
    _root_element (_isomorphism.pxi) had it verbatim -- so it lives here once.  Callers keep their
    own preconditions: _term_exact_element checks wbox_t.touched first, _root_element relies on a
    sealed root's neg[0] holding element bits only.  Neither precondition changes the arithmetic.

    Word 0 bit 0 is the heavy-element marker shared by every element above 56: a box that forbids
    it allows no heavy element at all, whatever its word-1 bits say.  An untouched element span
    needs no special case -- it leaves all 56 light bits open, and 56 != 1.
    """
    cdef uint64_t light = ~neg0 & <uint64_t> W0_LIGHT_ELEMENT_SPAN
    cdef uint64_t heavy

    if neg0 & <uint64_t> 1:
        heavy = 0
    else:
        heavy = ~neg1 & <uint64_t> W1_ELEMENT_SPAN
    if _popcount64(light) + _popcount64(heavy) != 1:
        return -1
    if light:
        return <int> (57 - _single_bit_index(light))
    return <int> (57 + _single_bit_index(heavy))


cdef int _term_exact_element(wterm_t *term) noexcept nogil:
    """The one element that EVERY box of this disjunction demands, or -1.

    This is the QSEG_ELEMENT_DEMAND rule, and it is deliberately stricter than
    _term_element_card: `[C,N]` has two boxes allowing one element each, but no molecule atom
    is guaranteed to be carbon, so it must contribute to no histogram slot.  Keeping the
    histogram a sound lower bound is what makes Task 14's screen safe.
    """
    cdef uint32_t i
    cdef int result = -1
    cdef int e
    cdef wbox_t *box

    if term.count == 0:
        return -1
    for i in range(term.count):
        box = &term.boxes[i]
        if not box.touched[0] & <uint64_t> W0_ELEMENT_SPAN:
            return -1
        e = _element_from_neg(box.neg[0], box.neg[1])
        if e < 0:
            return -1
        if result < 0:
            result = e
        elif result != e:
            return -1
    return result


cdef bint _term_is_metal_only(wterm_t *term) noexcept nogil:
    """Does EVERY box of this disjunction constrain the element to exactly the 93 metals?

    That is `[M]`, and the QATOM_METAL_ELEMENT rule.  PRIM_METAL is the one primitive that forbids
    exactly METAL_W*_FORBIDDEN, so the masks are compared rather than counted: a hand-written list
    of 93 elements that happened to have the metals' cardinality is not the metal wildcard, and
    `[M,C]` -- metals in one box, carbon in the other -- is not one either, because the atom it
    describes may be a carbon.

    Every box, not some box, for the same reason _term_exact_element demands every box: an atom's
    boxes are a disjunction, so a guarantee has to hold in all of them.
    """
    cdef uint32_t i
    cdef wbox_t *box

    if term.count == 0:
        return False
    for i in range(term.count):
        box = &term.boxes[i]
        if not box.touched[0] & <uint64_t> W0_ELEMENT_SPAN:
            return False
        if (box.neg[0] & <uint64_t> W0_ELEMENT_SPAN) != <uint64_t> METAL_W0_FORBIDDEN:
            return False
        if (box.neg[1] & <uint64_t> W1_ELEMENT_SPAN) != <uint64_t> METAL_W1_FORBIDDEN:
            return False
    return True


cdef int _fold_bond_into_atom(wterm_t *out, wterm_t *atom, wterm_t *bond) except -1:
    """Cross an atom's disjunction with its tree bond's, folding word 0 only.

    The folded neg[0] is tested against the incident HALF-EDGE word that fill_edge_words builds
    (_features.pxi), not against the candidate atom's own feature word 0.  A half-edge word
    carries the target atom's element bits plus exactly one topology bit and exactly one order
    bit -- the bond the DFS arrived by -- all in feature word 0's layout.  Every span in it is
    one-hot, so a forbidden-bit test over it is EXACT, and folding costs the matcher nothing:
    one AND settles the atom and the bond that reaches it together.

    It is NOT valid against the aggregate per-atom word 0 that fill_features writes.  That word
    ORs the topology and order bits of EVERY incident bond together, so acetone's carbonyl carbon
    carries the single-bond bit from its two methyls as well as its own double-bond bit, and a
    folded box that forbids single would reject it.  Task 10 may pass features[4 * candidate] for
    a ROOT position, whose box holds no folded bond demand; for every other position it must pass
    the half-edge word or query O=C-C stops matching CC(=O)C.

    Words 1-3 come from the atom box untouched, and that is load-bearing: compile_term runs
    box_fill_defaults over bond terms too, so every bond box carries a neutral-charge (word 2)
    and not-a-radical (word 1) demand that no bond feature could ever satisfy.  ORing those in
    would make every folded box unsatisfiable.  Task 11 checks closures against an edge word,
    which is word 0 as well, so words 1-3 of a bond box are never consulted anywhere.  Do not
    "fix" the defaults in box_fill_defaults without revisiting both.
    """
    cdef uint32_t i, j, k, box_idx, ai, aac, bac
    cdef wbox_t *ab
    cdef wbox_t *bb
    cdef wbox_t *ob

    if atom.count * bond.count > Q_ATOM_MAX_BOXES:
        raise ValueError('box count would exceed the cap of %d' % Q_ATOM_MAX_BOXES)
    memset(out, 0, sizeof(wterm_t))
    for i in range(atom.count):
        ab = &atom.boxes[i]
        aac = ab.any_count
        for j in range(bond.count):
            bb = &bond.boxes[j]
            bac = bb.any_count
            box_idx = i * bond.count + j
            ob = &out.boxes[box_idx]
            ob[0] = ab[0]
            ob.neg[0] |= bb.neg[0]
            ob.touched[0] |= bb.touched[0]
            # Bond boxes never carry any-entries today (no bond primitive is multi-hot), so this
            # concatenation is a copy.  It is written in full so a future multi-hot bond primitive
            # does not silently lose its demands -- but note that an any entry is checked against
            # the candidate's own feature word, f[word], not against the half-edge word the folded
            # neg[0] is tested against.  A multi-hot bond primitive therefore cannot simply be
            # folded here: it would need its own check on the closure/tree edge instead.
            if aac + bac > Q_BOX_MAX_ANY:
                raise ValueError('too many positive multi-hot constraints on one query atom; '
                                 'the cap is %d' % Q_BOX_MAX_ANY)
            for ai in range(bac):
                k = aac + ai
                ob.any_mask[k] = bb.any_mask[ai]
                ob.any_word[k] = bb.any_word[ai]
            ob.any_count = aac + bac
    out.count = atom.count * bond.count
    # Prune boxes the fold made impossible, then merge what the fold made mergeable.
    k = 0
    for i in range(out.count):
        if not box_unsatisfiable(&out.boxes[i]):
            if k != i:
                out.boxes[k] = out.boxes[i]
            k += 1
    out.count = k
    if out.count == 0:
        raise ValueError('this bond and the atom it leads to can never match together')
    boxes_merge(out)
    return 0


# ---------------------------------------------------------------------------
# The automorphism group — Task 12
# ---------------------------------------------------------------------------
# The group is computed INSIDE query_seal, from scratch that exists there and nowhere else: the
# UNFOLDED atom terms, the per-bond terms, and the query's CSR.  It cannot be recovered from a
# sealed arena, and that is not a stylistic preference:
#
#   * a non-root position's boxes are its atom term crossed with its tree bond's term (pass 6
#     folds them), so two identical atoms land in different classes the moment one of them roots
#     a component.  A partition over sealed boxes admits only the permutations that preserve the
#     DFS tree, which for cyclopropane collapses S3 to a single swap;
#   * qbond_t carries no endpoints, and a tree bond's box_count is 0 because its boxes moved into
#     the atom, so the arena records the query's adjacency but not its per-edge constraints.  A
#     pass running after seal cannot tell a single bond from a double one on a tree edge.
#
# There is deliberately NO query_append.  The arena is immutable after seal and has to stay that
# way: a suspended get_mapping generator caches raw arena pointers in its matcher_t, so a
# PyMem_Realloc driven by a second, filtered search on the same QueryContainer would move the
# buffer out from under it -- a use-after-free reachable from pure Python.

DEF Q_AUTOMORPHISM_MAX_ROWS = 1024      # twelve interchangeable atoms have 12! symmetries
DEF Q_AUTOMORPHISM_MAX_NODES = 100000   # backtracking budget, counted in candidates CONSIDERED
# One word for the box count, then per box: neg[4], the any count, and two words per any entry.
DEF Q_TERM_HASH_WORDS = 1 + Q_ATOM_MAX_BOXES * (5 + 2 * Q_BOX_MAX_ANY)


cdef bint _wterm_equal(wterm_t *a, wterm_t *b) noexcept nogil:
    """Byte-exact equality of two compiled disjunctions.

    This is the acceptance test behind every stored automorphism, and it compares BYTES on
    purpose: a hash here would let a collision manufacture a symmetry that is not one, and the
    filter would then drop legitimate matches.  Only the constraints are compared -- neg and the
    any list -- not wbox_t's construction bookkeeping (touched, element_set_size), which says
    nothing about what the term matches.

    Boxes are compared in order, so two terms holding the same boxes permuted compare unequal.
    That costs symmetry the filter would otherwise exploit, which makes it over-report; it can
    never invent a symmetry.

    The two callers are not equally dependent on this.  For ATOM terms the check is a collision
    guard by construction: the refinement seed hashes the whole atom term, so two slots that
    survive into one class already agree byte for byte and only an _xxh64 collision could bring an
    unequal pair here.  For BOND terms it is load-bearing -- the seed covers a bond term only
    through the neighbour multiset, which 1-WL can collapse (an even cycle with alternating bond
    labels leaves every slot in one class), so this is the only place the difference is seen.
    That asymmetry is why weakening the seed kills atom-term tests while the verification half is
    reachable only through bond terms; see the two alternating-ring tests in test_isomorphism.py.
    """
    cdef uint32_t i, k
    cdef wbox_t *ba
    cdef wbox_t *bb
    if a.count != b.count:
        return False
    for i in range(a.count):
        ba = &a.boxes[i]
        bb = &b.boxes[i]
        for k in range(4):
            if ba.neg[k] != bb.neg[k]:
                return False
        if ba.any_count != bb.any_count:
            return False
        for k in range(ba.any_count):
            if ba.any_mask[k] != bb.any_mask[k] or ba.any_word[k] != bb.any_word[k]:
                return False
    return True


cdef uint64_t _wterm_hash(wterm_t *term, uint64_t *scratch, uint64_t seed) noexcept nogil:
    """A 64-bit digest of everything _wterm_equal compares.  `scratch` holds Q_TERM_HASH_WORDS."""
    cdef uint32_t i, k, fill = 0
    cdef wbox_t *box
    scratch[fill] = <uint64_t> term.count
    fill += 1
    for i in range(term.count):
        box = &term.boxes[i]
        for k in range(4):
            scratch[fill] = box.neg[k]
            fill += 1
        scratch[fill] = <uint64_t> box.any_count
        fill += 1
        for k in range(box.any_count):
            scratch[fill] = box.any_mask[k]
            fill += 1
            scratch[fill] = <uint64_t> box.any_word[k]
            fill += 1
    return _xxh64(scratch, fill, seed)


cdef inline uint32_t _bond_between(uint32_t *csr_head, uint32_t *csr_nbr, uint32_t *csr_bnd,
                                  uint32_t s, uint32_t t) noexcept nogil:
    """The bond slot joining two query slots, or Q_NO_SLOT.  A CSR scan: query degrees are tiny."""
    cdef uint32_t k
    for k in range(csr_head[s], csr_head[s + 1]):
        if csr_nbr[k] == t:
            return csr_bnd[k]
    return Q_NO_SLOT


cdef int _compute_automorphisms(uint32_t atom_count, uint32_t bond_count,
                               wterm_t *atom_terms, wterm_t *bond_terms,
                               uint32_t *csr_head, uint32_t *csr_nbr, uint32_t *csr_bnd,
                               uint32_t *comp_of, int32_t *comp_group, uint32_t *pos_of,
                               uint32_t **rows_out, uint32_t *row_count_out,
                               uint32_t *flags_out) except -1:
    """The query's automorphism group, as permutations of DFS positions.

    THE VERIFICATION IS THE AUTHORITY; THE PARTITION IS ONLY A PRUNER.  A candidate permutation
    is accepted only when, checked exactly:

      * every slot's UNFOLDED atom term is byte-equal to its image's, and its component's group
        number equals its image's (a permutation that moves a slot from a group-0 component into a
        group-1 one is not a symmetry of the constraint system component groups impose), and
      * the edge set maps onto itself: for each pair of slots a bond exists on one side exactly
        when it exists on the other, and the two bond terms are byte-equal.

    The refinement below only narrows which candidates are worth trying.  A partition that is too
    coarse costs candidates the verification then rejects; a partition that is too fine loses
    automorphisms, which makes the filter filter LESS -- over-reporting matches, never dropping a
    real one.  Both directions are safe.  A wrong verification is not, which is why _xxh64 appears
    in the refinement and never in the acceptance test.

    On success rows_out[0] is NULL with row_count_out[0] == 0 for the trivial group, or a PyMem
    block of row_count_out[0] * atom_count uint32_t that the caller owns and must free.
    flags_out[0] receives QFLAG_ASYMMETRIC when the group is trivial, or
    QFLAG_PARTIAL_AUTOMORPHISM when either cap was hit.  A partial group filters less and
    therefore over-reports; it can never drop a match, because the lexicographically smallest
    member of an orbit is still smallest when tested against a SUBSET of the group.

    THESE ROWS MUST NEVER BE USED TO DERIVE ORBITS.  The enumeration is lexicographic in slot
    order, so a truncated one holds permutations of the LAST slots only and its union claims the
    first slots are fixed: six waters in one record give 5 orbits where there are 2.  That
    under-reports symmetry, which for stereo means inventing stereocentres.  `mapping_is_canonical`
    is the sole consumer and it needs any subgroup, which is why the truncation is safe for it and
    for nothing else.  The molecule side asks its own question per unresolved pair instead --
    see `mol_automorphisms` in _canonical.pxi.
    """
    cdef uint64_t *scratch = NULL
    cdef uint64_t *key = NULL
    cdef uint64_t *nb = NULL
    cdef uint64_t *bond_hash = NULL
    cdef uint32_t *cls = NULL
    cdef uint32_t *nxt = NULL
    cdef uint32_t *idx = NULL
    cdef uint32_t *tmp = NULL
    cdef uint32_t *sigma = NULL
    cdef uint32_t *cursor = NULL
    cdef uint32_t *rows = NULL
    cdef uint8_t *taken = NULL
    cdef uint64_t pair[2]
    cdef uint32_t s, t, u, k, d, cand, b1, b2, maxdeg = 0, depth = 0, rounds = 0
    cdef uint32_t nodes = 0, row_count = 0
    cdef Py_ssize_t classes = 0, newclasses = 0
    cdef bint ok = True, found = False, identity = True, overflow = False

    rows_out[0] = NULL
    row_count_out[0] = 0
    flags_out[0] = 0
    if atom_count < 2:
        flags_out[0] = QFLAG_ASYMMETRIC
        return 0
    for s in range(atom_count):
        d = csr_head[s + 1] - csr_head[s]
        if d > maxdeg:
            maxdeg = d

    try:
        scratch = <uint64_t *> PyMem_Malloc(Q_TERM_HASH_WORDS * sizeof(uint64_t))
        key = <uint64_t *> PyMem_Malloc(atom_count * sizeof(uint64_t))
        nb = <uint64_t *> PyMem_Malloc(_alloc_at_least(maxdeg) * sizeof(uint64_t))
        bond_hash = <uint64_t *> PyMem_Malloc(_alloc_at_least(bond_count) * sizeof(uint64_t))
        cls = <uint32_t *> PyMem_Malloc(atom_count * sizeof(uint32_t))
        nxt = <uint32_t *> PyMem_Malloc(atom_count * sizeof(uint32_t))
        idx = <uint32_t *> PyMem_Malloc(atom_count * sizeof(uint32_t))
        tmp = <uint32_t *> PyMem_Malloc(atom_count * sizeof(uint32_t))
        if (scratch is NULL or key is NULL or nb is NULL or bond_hash is NULL or cls is NULL or
                nxt is NULL or idx is NULL or tmp is NULL):
            raise MemoryError('automorphism scratch allocation failed')

        # Round 0: the atom term and the owning component's group.  Rounds after that fold in the
        # neighbours' classes and the joining bond's term, sorted so the CSR's incidental
        # neighbour order cannot leak in.  Equal inputs hash equal, so two slots an automorphism
        # relates can never be driven into different classes, whatever the hash does.
        for s in range(atom_count):
            key[s] = _wterm_hash(&atom_terms[s], scratch,
                                 <uint64_t> <int64_t> comp_group[comp_of[s]])
        for k in range(bond_count):
            bond_hash[k] = _wterm_hash(&bond_terms[k], scratch, 0)
        classes = _classify(key, idx, tmp, atom_count, cls)
        while classes < <Py_ssize_t> atom_count and rounds < atom_count:
            rounds += 1
            for s in range(atom_count):
                d = 0
                for k in range(csr_head[s], csr_head[s + 1]):
                    pair[0] = <uint64_t> cls[csr_nbr[k]]
                    pair[1] = bond_hash[csr_bnd[k]]
                    nb[d] = _xxh64(pair, 2, 0)
                    d += 1
                _sort_words(nb, d)
                key[s] = _xxh64(nb, d, <uint64_t> cls[s])
            newclasses = _classify(key, idx, tmp, atom_count, nxt)
            if newclasses <= classes:
                # A fixed point, or a hash collision merged two classes.  Keep the previous
                # partition, which is at least as fine, and stop: neither direction is unsafe --
                # see the docstring -- and refining further cannot help once it stops splitting.
                break
            classes = newclasses
            memcpy(cls, nxt, <size_t> atom_count * sizeof(uint32_t))

        if classes == <Py_ssize_t> atom_count:
            # Every class a singleton.  An automorphism preserves the refinement, so the only one
            # left is the identity.  This is the common case and it costs nothing beyond the
            # refinement -- no enumeration and no row buffer.
            flags_out[0] = QFLAG_ASYMMETRIC
            return 0

        rows = <uint32_t *> PyMem_Malloc(
            <size_t> Q_AUTOMORPHISM_MAX_ROWS * atom_count * sizeof(uint32_t))
        sigma = <uint32_t *> PyMem_Malloc(atom_count * sizeof(uint32_t))
        cursor = <uint32_t *> PyMem_Malloc(atom_count * sizeof(uint32_t))
        taken = <uint8_t *> PyMem_Malloc(atom_count * sizeof(uint8_t))
        if rows is NULL or sigma is NULL or cursor is NULL or taken is NULL:
            raise MemoryError('automorphism scratch allocation failed')
        memset(taken, 0, atom_count)
        for s in range(atom_count):
            sigma[s] = Q_NO_SLOT
        cursor[0] = 0

        # Backtracking over slots in slot order: sigma[0 .. depth-1] is a partial injection whose
        # verification already holds on every pair inside it, so a complete assignment needs no
        # further check.
        while True:
            found = False
            t = cursor[depth]
            while t < atom_count:
                cand = t
                t += 1
                if taken[cand] or cls[cand] != cls[depth]:
                    continue
                nodes += 1
                if nodes > Q_AUTOMORPHISM_MAX_NODES:
                    overflow = True
                    break
                if comp_group[comp_of[cand]] != comp_group[comp_of[depth]]:
                    continue
                if not _wterm_equal(&atom_terms[depth], &atom_terms[cand]):
                    continue
                ok = True
                for u in range(depth):
                    b1 = _bond_between(csr_head, csr_nbr, csr_bnd, depth, u)
                    b2 = _bond_between(csr_head, csr_nbr, csr_bnd, cand, sigma[u])
                    if b1 == Q_NO_SLOT:
                        if b2 != Q_NO_SLOT:
                            ok = False
                            break
                    elif b2 == Q_NO_SLOT:
                        ok = False
                        break
                    elif not _wterm_equal(&bond_terms[b1], &bond_terms[b2]):
                        ok = False
                        break
                if not ok:
                    continue
                sigma[depth] = cand
                taken[cand] = 1
                found = True
                break
            cursor[depth] = t
            if overflow:
                break
            if found:
                depth += 1
                if depth < atom_count:
                    cursor[depth] = 0
                    continue
                identity = True
                for s in range(atom_count):
                    if sigma[s] != s:
                        identity = False
                        break
                if not identity:
                    if row_count == Q_AUTOMORPHISM_MAX_ROWS:
                        overflow = True
                    else:
                        # Rows are permutations of POSITIONS: mapping_is_canonical indexes
                        # m.mapping by position, and slot order is not position order.
                        for s in range(atom_count):
                            rows[row_count * atom_count + pos_of[s]] = pos_of[sigma[s]]
                        row_count += 1
                depth -= 1
                taken[sigma[depth]] = 0
                sigma[depth] = Q_NO_SLOT
                if overflow:
                    break
            elif depth == 0:
                break
            else:
                depth -= 1
                taken[sigma[depth]] = 0
                sigma[depth] = Q_NO_SLOT

        if overflow:
            flags_out[0] = QFLAG_PARTIAL_AUTOMORPHISM
        if row_count == 0:
            PyMem_Free(rows)
            rows = NULL
            if not overflow:
                # The refinement left a class with more than one member, but the verification
                # rejected every non-identity candidate: the group really is trivial.
                flags_out[0] = QFLAG_ASYMMETRIC
        row_count_out[0] = row_count
        rows_out[0] = rows
        rows = NULL          # ownership handed to the caller; the finally must not free it
    finally:
        PyMem_Free(scratch)
        PyMem_Free(key)
        PyMem_Free(nb)
        PyMem_Free(bond_hash)
        PyMem_Free(cls)
        PyMem_Free(nxt)
        PyMem_Free(idx)
        PyMem_Free(tmp)
        PyMem_Free(sigma)
        PyMem_Free(cursor)
        PyMem_Free(taken)
        PyMem_Free(rows)
    return 0


# ------------------------------------------------------------------------------------------------
# THE GEOMETRY A QUERY ASKS FOR: `/` AND `\` READ BACK OFF THE JOURNAL
# ------------------------------------------------------------------------------------------------
#
# A direction compiles to no box, because which side of a double bond a substituent sits on is not a
# property of the bond it is written on.  What it is half of -- a statement about the double bond's two
# ends -- is what QSEG_STEREO carries: one record per chain of double bonds, naming the two terminals,
# the substituent a direction marked on each, and whether the two stand on opposite sides.  The kernel
# reads the target's cis/trans unit in that frame, the way it reads a centre's in the query's F26 order.
#
# READ FROM THE CHAIN'S TERMINALS, not from the marked bonds, for the reason the SMIRKS product side is
# (`smk_directions`): one single bond between two chains marks a side for both of them, so `C/C=C/C=C/C`
# is three marks making two geometries.
#
# A CHAIN BOND IS A BOND WHOSE EXPRESSION IS EXACTLY `=`.  `[C]=,#[C]` pins no double bond, so a
# direction beside one falls to the "names no chain" refusal rather than being guessed at.

cdef struct qgeom_t:            # one geometry a query states, in query SLOTS
    uint32_t term_a
    uint32_t term_b
    uint32_t mark_a             # the substituent of term_a that a direction named
    uint32_t mark_b
    uint8_t trans               # the two marks stand on opposite sides


cdef int _seal_geometries(qop_t *ops, uint32_t op_count, uint32_t *slot_of, uint32_t next_id,
                          uint32_t *stable_of, uint32_t *bond_a, uint32_t *bond_b,
                          uint32_t bond_count, uint32_t atom_count, uint32_t *csr_head,
                          uint32_t *csr_nbr, uint32_t *csr_bnd, qgeom_t *out,
                          uint32_t *count_out) except -1:
    """Group the journal's directions into geometries; `out` takes one record per geometry.

    Called after the CSR is built and before QSEG_STEREO is sized, so `out` has room for one record
    per bond.  Every refusal here is a ValueError, which is what the readers turn into their own
    exception -- a direction that states nothing is a defect in the string, not in the target.
    """
    cdef uint8_t *dbl = NULL         # per bond: its expression is exactly `=`
    cdef uint8_t *bdir = NULL        # per bond: SMI_DIR_UP / SMI_DIR_DOWN, 0 for no direction
    cdef uint8_t *used = NULL        # per bond: a geometry spent this direction
    cdef uint32_t *bfrom = NULL      # per bond: the slot the direction was written FROM
    cdef uint8_t *seen = NULL        # per atom: a chain already walked through it
    cdef uint32_t i, k, b, s, t, arms, other, prev, cur, hops, tokens, count = 0
    cdef uint32_t term[2]
    cdef uint32_t mark[2]
    cdef int side[2]
    cdef int d
    cdef bint plain, reached
    cdef qop_t *qop

    count_out[0] = 0
    if not bond_count:
        return 0
    try:
        dbl = <uint8_t *> PyMem_Malloc(_alloc_at_least(bond_count) * sizeof(uint8_t))
        bdir = <uint8_t *> PyMem_Malloc(_alloc_at_least(bond_count) * sizeof(uint8_t))
        used = <uint8_t *> PyMem_Malloc(_alloc_at_least(bond_count) * sizeof(uint8_t))
        bfrom = <uint32_t *> PyMem_Malloc(_alloc_at_least(bond_count) * sizeof(uint32_t))
        seen = <uint8_t *> PyMem_Malloc(_alloc_at_least(atom_count) * sizeof(uint8_t))
        if dbl is NULL or bdir is NULL or used is NULL or bfrom is NULL or seen is NULL:
            raise MemoryError('query seal scratch allocation failed')
        memset(dbl, 0, _alloc_at_least(bond_count) * sizeof(uint8_t))
        memset(bdir, 0, _alloc_at_least(bond_count) * sizeof(uint8_t))
        memset(used, 0, _alloc_at_least(bond_count) * sizeof(uint8_t))
        memset(seen, 0, _alloc_at_least(atom_count) * sizeof(uint8_t))

        for i in range(op_count):
            qop = &ops[i]
            if qop.op != QOP_SET_BOND_DIRECTION:
                continue
            b = _find_bond_slot(bond_a, bond_b, bond_count, slot_of, next_id, qop.a, qop.b)
            if bdir[b] and (bdir[b] != <uint8_t> qop.value or bfrom[b] != slot_of[qop.a]):
                raise ValueError('bond %d-%d carries two directions, which are two statements about '
                                 'one side' % (qop.a, qop.b))
            bdir[b] = <uint8_t> qop.value
            bfrom[b] = slot_of[qop.a]

        for b in range(bond_count):
            tokens = 0
            plain = False
            for i in range(op_count):
                qop = &ops[i]
                if qop.op != QOP_BOND_TOKEN or _find_bond_slot(bond_a, bond_b, bond_count, slot_of,
                                                               next_id, qop.a, qop.b) != b:
                    continue
                tokens += 1
                plain = (qop.opcode == OPC_PRIM and qop.kind == BPRIM_ORDER and qop.value == 2
                         and not qop.negated)
            if tokens == 1 and plain:
                dbl[b] = 1

        for s in range(atom_count):
            if seen[s]:
                continue
            arms = 0
            for k in range(csr_head[s], csr_head[s + 1]):
                if dbl[csr_bnd[k]]:
                    arms += 1
            if arms != 1:
                continue
            # Walk to the far terminal.  An all-double ring has no atom with one arm and is never
            # entered; an atom with three is no chain end this can order, and the walk gives up on it
            # -- its directions then fall to the "names no chain" refusal below.
            cur = s
            prev = Q_NO_SLOT
            hops = 0
            reached = False
            while hops <= bond_count:
                arms = 0
                other = Q_NO_SLOT
                for k in range(csr_head[cur], csr_head[cur + 1]):
                    if dbl[csr_bnd[k]] and csr_nbr[k] != prev:
                        arms += 1
                        other = csr_nbr[k]
                if not arms:
                    reached = True
                    break
                if arms > 1:
                    break
                prev = cur
                cur = other
                hops += 1
            if not reached:
                continue
            term[0] = s
            term[1] = cur
            seen[s] = 1
            seen[cur] = 1

            for i in range(2):
                t = term[i]
                mark[i] = Q_NO_SLOT
                side[i] = 0
                for k in range(csr_head[t], csr_head[t + 1]):
                    b = csr_bnd[k]
                    if dbl[b] or not bdir[b]:
                        continue
                    # the same statement read from the other end is upside down
                    d = bdir[b] if bfrom[b] == t else 3 - <int> bdir[b]
                    used[b] = 1
                    if mark[i] == Q_NO_SLOT:
                        mark[i] = csr_nbr[k]
                        side[i] = d
                    elif d == side[i]:
                        raise ValueError('atom %d puts both of its substituents on the same side of '
                                         'the double bond it terminates, and no geometry does that'
                                         % stable_of[t])
            if mark[0] == Q_NO_SLOT and mark[1] == Q_NO_SLOT:
                continue
            if mark[0] == Q_NO_SLOT or mark[1] == Q_NO_SLOT:
                raise ValueError('the double bond between atoms %d and %d carries a direction on one '
                                 'end only; a geometry is a statement about both, so a `/` or `\\` is '
                                 'needed on a substituent of each'
                                 % (stable_of[term[0]], stable_of[term[1]]))
            # opposite directions, each read from its own terminal, means opposite sides
            if term[0] < term[1]:
                out[count].term_a = term[0]
                out[count].term_b = term[1]
                out[count].mark_a = mark[0]
                out[count].mark_b = mark[1]
            else:
                out[count].term_a = term[1]
                out[count].term_b = term[0]
                out[count].mark_a = mark[1]
                out[count].mark_b = mark[0]
            out[count].trans = 1 if side[0] != side[1] else 0
            count += 1

        for b in range(bond_count):
            if bdir[b] and not used[b]:
                raise ValueError('bond %d-%d carries a `/` or `\\` that names no chain of double '
                                 'bonds, and a direction states nothing on its own: it says which '
                                 'side of a geometry a substituent is on, so there has to be a '
                                 'geometry for it to be part of'
                                 % (stable_of[bond_a[b]], stable_of[bond_b[b]]))
        count_out[0] = count
    finally:
        PyMem_Free(dbl)
        PyMem_Free(bdir)
        PyMem_Free(used)
        PyMem_Free(bfrom)
        PyMem_Free(seen)
    return 0


cdef Query query_seal(qop_t *ops, uint32_t op_count, uint32_t next_id,
                      uint32_t **position_to_n):
    """Seal a query journal into an arena laid out as a DFS plan.

    `next_id` is one past the highest stable id the caller handed out, so the id-to-slot table
    is a plain array.  On success `position_to_n` receives a freshly PyMem_Malloc'ed
    uint32_t[atom_count] that the caller owns and must free: the arena is indexed by DFS
    position and deliberately holds no caller identifiers, but Task 10 reports mappings in the
    caller's namespace.
    """
    cdef uint32_t *slot_of = NULL
    cdef uint32_t *stable_of = NULL
    cdef uint8_t *masked_of = NULL
    cdef uint16_t *map_of = NULL
    cdef int32_t *group_of = NULL
    cdef uint32_t *bond_a = NULL
    cdef uint32_t *bond_b = NULL
    cdef wterm_t *atom_terms = NULL
    cdef wterm_t *bond_terms = NULL
    cdef qtoken_t *scratch = NULL
    cdef qtoken_t *tok
    cdef uint32_t *csr_head = NULL
    cdef uint32_t *csr_cur = NULL
    cdef uint32_t *csr_nbr = NULL
    cdef uint32_t *csr_bnd = NULL
    cdef uint32_t *comp_of = NULL
    cdef uint32_t *order = NULL
    cdef uint32_t *back = NULL
    cdef uint32_t *back_bond = NULL
    cdef uint32_t *pos_of = NULL
    cdef uint32_t *stack = NULL
    cdef uint8_t *visited = NULL
    cdef uint32_t *card = NULL
    cdef uint32_t *rarity = NULL
    # Per slot, the element-wildcard flags QATOM_ANY_ELEMENT / QATOM_METAL_ELEMENT.  Computed with
    # `card` because both read the PRE-FOLD term, and the fold overwrites atom_terms in place.
    cdef uint8_t *wild_of = NULL
    cdef uint8_t *is_tree = NULL
    cdef uint32_t *clo_owner = NULL
    cdef uint32_t *clo_to = NULL
    cdef uint32_t *clo_bond = NULL
    cdef int32_t *comp_group_of = NULL
    cdef uint32_t *auto_rows = NULL
    # Per slot, the signs its stereo primitives named: bit 0 for '@' (1), bit 1 for '@@' (2).
    cdef uint8_t *stereo_of = NULL
    cdef uint32_t demand[120]

    cdef uint32_t atom_count = 0, bond_count = 0, comp_count = 0, closure_count = 0
    cdef uint32_t i, k, fill, s, t, u, v, n, m, deg, next_pos, sp, root, best, best_deg
    # `comp` walks components and `n_nbrs` counts one slot's neighbours. Both were spelled `n`
    # before the n/m sweep gave that letter to the atom; they MUST stay declared, because Cython
    # answers an undeclared local with an inferred Python object and only a build warning.
    cdef uint32_t comp, n_nbrs
    cdef uint32_t best_bond, bslot, box_total, any_total, bond_box_total, pa, pb, flags
    cdef uint32_t auto_count = 0, auto_flags = 0
    cdef int32_t comp_group, exact
    cdef uint32_t box_cursor, any_cursor, bond_box_cursor, clo_cursor
    cdef wterm_t folded
    cdef wterm_t *wt
    cdef wbox_t *wbx
    cdef Query q = None
    cdef qop_t *qop
    cdef qatom_t *out_atoms
    cdef qatom_t *qa
    cdef qbox_t *out_boxes
    cdef qany_t *out_any
    cdef qbond_t *out_bonds
    cdef qbond_t *qb
    cdef qbox_t *out_bond_boxes
    cdef qclosure_t *out_closures
    cdef qcomp_t *out_comps
    cdef qcomp_t *qc
    cdef uint32_t *out_demand
    cdef uint32_t *out_demand_list
    cdef qstereo_t *out_stereo
    cdef qstereo_t *qs
    cdef qgeom_t *geoms = NULL
    cdef uint32_t geom_count = 0
    cdef uint32_t stereo_count = 0, stereo_cursor, ready, swap
    cdef uint32_t nbrs[4]
    cdef uint32_t demand_list_count = 0, dl_cursor
    cdef uint32_t *ids
    cdef object exc
    # Task 14: signature words (screen demand) — computed post-fold, written to header
    cdef uint64_t sig[4]
    cdef uint64_t atom_sig[4]
    cdef uint64_t span_contrib, allowed_bits, M
    cdef uint32_t si, w
    cdef int e_common, e_box

    memset(demand, 0, sizeof(demand))
    try:
        # --- pass 1: dense slots for atoms and bonds -----------------------------------
        # Two walks so a bond may legally precede its atoms in the journal.
        for i in range(op_count):
            qop = &ops[i]
            if qop.op == QOP_ADD_ATOM:
                atom_count += 1
            elif qop.op == QOP_ADD_BOND:
                bond_count += 1

        slot_of = <uint32_t *> PyMem_Malloc(_alloc_at_least(next_id) * sizeof(uint32_t))
        stable_of = <uint32_t *> PyMem_Malloc(_alloc_at_least(atom_count) * sizeof(uint32_t))
        masked_of = <uint8_t *> PyMem_Malloc(_alloc_at_least(atom_count) * sizeof(uint8_t))
        map_of = <uint16_t *> PyMem_Malloc(_alloc_at_least(atom_count) * sizeof(uint16_t))
        group_of = <int32_t *> PyMem_Malloc(_alloc_at_least(atom_count) * sizeof(int32_t))
        bond_a = <uint32_t *> PyMem_Malloc(_alloc_at_least(bond_count) * sizeof(uint32_t))
        bond_b = <uint32_t *> PyMem_Malloc(_alloc_at_least(bond_count) * sizeof(uint32_t))
        if (slot_of is NULL or stable_of is NULL or masked_of is NULL or map_of is NULL or
                group_of is NULL or bond_a is NULL or bond_b is NULL):
            raise MemoryError('query seal scratch allocation failed')
        for i in range(next_id):
            slot_of[i] = Q_NO_SLOT
        for i in range(atom_count):
            masked_of[i] = 0
            map_of[i] = 0
            group_of[i] = -1

        fill = 0
        for i in range(op_count):
            qop = &ops[i]
            if qop.op == QOP_ADD_ATOM:
                n = qop.a
                if n >= next_id:
                    raise ValueError('atom id %d is past next_id %d' % (n, next_id))
                if slot_of[n] != Q_NO_SLOT:
                    raise ValueError('atom %d was added twice' % n)
                slot_of[n] = fill
                stable_of[fill] = n
                fill += 1

        fill = 0
        for i in range(op_count):
            qop = &ops[i]
            if qop.op != QOP_ADD_BOND:
                continue
            n = qop.a
            m = qop.b
            if n == m:
                raise ValueError('bond %d-%d is a self loop' % (n, m))
            if n >= next_id or slot_of[n] == Q_NO_SLOT:
                raise ValueError('bond %d-%d references unknown atom %d' % (n, m, n))
            if m >= next_id or slot_of[m] == Q_NO_SLOT:
                raise ValueError('bond %d-%d references unknown atom %d' % (n, m, m))
            s = slot_of[n]
            t = slot_of[m]
            # A linear scan, not a hash table: a query with more than a few dozen bonds does not
            # exist in practice, and seal runs once per query.
            for k in range(fill):
                if (bond_a[k] == s and bond_b[k] == t) or (bond_a[k] == t and bond_b[k] == s):
                    raise ValueError('duplicate bond %d-%d' % (n, m))
            bond_a[fill] = s
            bond_b[fill] = t
            fill += 1

        # atom-scoped journal entries, and the validation Ruling 4 asks for
        for i in range(op_count):
            qop = &ops[i]
            if (qop.op == QOP_ATOM_TOKEN or qop.op == QOP_SET_GROUP or
                    qop.op == QOP_SET_MASKED or qop.op == QOP_SET_MAP or
                    qop.op == QOP_SET_STEREO_GROUP):
                n = qop.a
                if n >= next_id or slot_of[n] == Q_NO_SLOT:
                    raise ValueError('unknown atom %d' % n)
                s = slot_of[n]
                if qop.op == QOP_SET_STEREO_GROUP:
                    # THE ONE JOURNAL OP A QUERY CANNOT COMPILE, and it is refused here so that both
                    # doors are shut by one line: `read_smarts` seals, and so does a SMIRKS reactant
                    # side, so `[C;&1]` fails at the string on either.  An enhanced-stereo group says
                    # a configuration is one of a SET -- a fact about a molecule and its mixture, not
                    # a property of an atom -- and there is nothing in a target to compare it with.
                    # It is journalled rather than rejected in the lexer because the SMIRKS product
                    # side is the one place it means something, and that side never seals.
                    raise ValueError('atom %d carries an enhanced-stereo group; a group states that '
                                     'a configuration is one of a set, which is a fact about a '
                                     'molecule rather than something a query can test' % n)
                elif qop.op == QOP_SET_GROUP:
                    group_of[s] = qop.value
                elif qop.op == QOP_SET_MASKED:
                    masked_of[s] = 1
                elif qop.op == QOP_SET_MAP:
                    if qop.value < 0 or qop.value > MAP_NUMBER_MAX:
                        raise ValueError('map number %d is out of range 0..%d'
                                         % (qop.value, MAP_NUMBER_MAX))
                    map_of[s] = <uint16_t> qop.value
            elif qop.op == QOP_BOND_TOKEN or qop.op == QOP_SET_BOND_DIRECTION:
                # A token for a pair that never got a QOP_ADD_BOND would otherwise vanish
                # silently, dropping a constraint the caller wrote.
                if _find_bond_slot(bond_a, bond_b, bond_count, slot_of, next_id,
                                   qop.a, qop.b) == Q_NO_SLOT:
                    raise ValueError('unknown bond %d-%d' % (qop.a, qop.b))
                # A direction is validated here as a bond reference and read as a geometry later, by
                # `_seal_geometries`: on its own it names half a statement, and the whole one needs
                # the CSR to find the chain of double bonds the two halves are about.

        # --- pass 2: compile every atom and every bond ---------------------------------
        # wterm_t is ~3.6 KB, far too large for a stack array of n.
        atom_terms = <wterm_t *> PyMem_Malloc(_alloc_at_least(atom_count) * sizeof(wterm_t))
        bond_terms = <wterm_t *> PyMem_Malloc(_alloc_at_least(bond_count) * sizeof(wterm_t))
        scratch = <qtoken_t *> PyMem_Malloc((op_count + 1) * sizeof(qtoken_t))
        stereo_of = <uint8_t *> PyMem_Malloc(_alloc_at_least(atom_count) * sizeof(uint8_t))
        if atom_terms is NULL or bond_terms is NULL or scratch is NULL or stereo_of is NULL:
            raise MemoryError('query seal scratch allocation failed')
        memset(stereo_of, 0, _alloc_at_least(atom_count) * sizeof(uint8_t))

        for s in range(atom_count):
            n = stable_of[s]
            fill = 0
            for i in range(op_count):
                qop = &ops[i]
                if qop.op == QOP_ATOM_TOKEN and slot_of[qop.a] == s:
                    _token_from_op(&scratch[fill], qop)
                    fill += 1
            if fill == 0:
                # Even [A] emits PRIM_ANY, so a token-free atom means an unclosed bracket.
                raise ValueError('atom %d has no primitives' % n)
            try:
                compile_term(&atom_terms[s], scratch, fill)
            except ValueError as exc:
                # ValueError only: every raise reachable from prim_apply is one, the stereo
                # validation included.
                raise ValueError('atom %d: %s' % (n, exc))
            # Ruling F87: the sign is decided PER BOX, and the kernel reads it off whichever box
            # admitted the candidate (matcher_next -> _stereo_sign_mask).  What is collected here
            # is only the union over the atom's boxes, and it is used for exactly two things: the
            # QATOM_STEREO_* flags, which tell the kernel "this atom has at least one box with a
            # sign, so compute the mask", and the qstereo_t.sign field, which is reporting only.
            # It must be read off the compiled boxes rather than the raw tokens: a token scan
            # cannot tell '[C;@,N]' (one box demands @, one demands nothing) from '[C;@;@@]' (one
            # box demands both), and that conflation is what dropped embeddings before F87.
            for i in range(atom_terms[s].count):
                stereo_of[s] |= atom_terms[s].boxes[i].sign

        for bslot in range(bond_count):
            fill = 0
            for i in range(op_count):
                qop = &ops[i]
                if qop.op == QOP_BOND_TOKEN and \
                        _find_bond_slot(bond_a, bond_b, bond_count, slot_of, next_id,
                                        qop.a, qop.b) == bslot:
                    _token_from_op(&scratch[fill], qop)
                    fill += 1
            if fill == 0:
                # An implicit bond in SMARTS means a single bond.
                tok = &scratch[0]
                tok.opcode = OPC_PRIM
                tok.kind = BPRIM_ORDER
                tok.value = 1
                tok.negated = False
                fill = 1
            try:
                compile_term(&bond_terms[bslot], scratch, fill)
            except ValueError as exc:
                raise ValueError('bond %d-%d: %s' % (stable_of[bond_a[bslot]],
                                                     stable_of[bond_b[bslot]], exc))

        # --- element demand histogram and the per-atom element cardinality -------------
        card = <uint32_t *> PyMem_Malloc(_alloc_at_least(atom_count) * sizeof(uint32_t))
        rarity = <uint32_t *> PyMem_Malloc(_alloc_at_least(atom_count) * sizeof(uint32_t))
        wild_of = <uint8_t *> PyMem_Malloc(_alloc_at_least(atom_count) * sizeof(uint8_t))
        if card is NULL or rarity is NULL or wild_of is NULL:
            raise MemoryError('query seal scratch allocation failed')
        for s in range(atom_count):
            card[s] = _term_element_card(&atom_terms[s])
            # The element-wildcard witness, derived here rather than journalled: Q_NO_SLOT is
            # _term_element_card's own answer for "some box names no element at all", so the two
            # facts come off one walk of the pre-fold term.
            if card[s] == Q_NO_SLOT:
                wild_of[s] = QATOM_ANY_ELEMENT
            elif _term_is_metal_only(&atom_terms[s]):
                wild_of[s] = QATOM_METAL_ELEMENT
            else:
                wild_of[s] = 0
            exact = _term_exact_element(&atom_terms[s])
            if exact > 0 and exact < 120:
                demand[exact] += 1
                rarity[s] = <uint32_t> exact         # the element, for now
            else:
                rarity[s] = Q_NO_SLOT
        # The histogram is built here, PRE-fold, while the boxes the arena ends up holding are
        # post-fold.  That is sound only because of boxes_merge's element-words guard: the fold
        # itself never touches element bits, pruning only removes boxes, and boxes_merge refuses
        # to merge two boxes whose element words differ, so no post-fold box can allow an element
        # its pre-fold ancestor forbade.  Drop that guard and this histogram silently over-counts.
        #
        # Second pass, once the histogram is complete: rarity becomes how many query atoms
        # demand this atom's element.  It breaks ties between atoms of equal cardinality, so
        # C-C-O roots at the oxygen rather than at the middle carbon.
        for s in range(atom_count):
            if rarity[s] != Q_NO_SLOT:
                rarity[s] = demand[rarity[s]]

        # --- pass 3: CSR over the query graph ------------------------------------------
        # Count degrees, prefix-sum, scatter.  csr_cur is the per-row write cursor so csr_head
        # keeps the classic ptr form the DFS reads degrees from.
        csr_head = <uint32_t *> PyMem_Malloc((atom_count + 1) * sizeof(uint32_t))
        csr_cur = <uint32_t *> PyMem_Malloc(_alloc_at_least(atom_count) * sizeof(uint32_t))
        csr_nbr = <uint32_t *> PyMem_Malloc(_alloc_at_least(2 * bond_count) * sizeof(uint32_t))
        csr_bnd = <uint32_t *> PyMem_Malloc(_alloc_at_least(2 * bond_count) * sizeof(uint32_t))
        if csr_head is NULL or csr_cur is NULL or csr_nbr is NULL or csr_bnd is NULL:
            raise MemoryError('query seal scratch allocation failed')
        for s in range(atom_count + 1):
            csr_head[s] = 0
        for bslot in range(bond_count):
            csr_head[bond_a[bslot] + 1] += 1
            csr_head[bond_b[bslot] + 1] += 1
        for s in range(atom_count):
            csr_head[s + 1] += csr_head[s]
            csr_cur[s] = 0
        for bslot in range(bond_count):
            s = bond_a[bslot]
            t = bond_b[bslot]
            csr_nbr[csr_head[s] + csr_cur[s]] = t
            csr_bnd[csr_head[s] + csr_cur[s]] = bslot
            csr_cur[s] += 1
            csr_nbr[csr_head[t] + csr_cur[t]] = s
            csr_bnd[csr_head[t] + csr_cur[t]] = bslot
            csr_cur[t] += 1

        # --- the geometries `/` and `\` state, which need the CSR and nothing later ----
        geoms = <qgeom_t *> PyMem_Malloc(_alloc_at_least(bond_count) * sizeof(qgeom_t))
        if geoms is NULL:
            raise MemoryError('query seal scratch allocation failed')
        _seal_geometries(ops, op_count, slot_of, next_id, stable_of, bond_a, bond_b, bond_count,
                         atom_count, csr_head, csr_nbr, csr_bnd, geoms, &geom_count)

        # --- pass 4: components --------------------------------------------------------
        comp_of = <uint32_t *> PyMem_Malloc(_alloc_at_least(atom_count) * sizeof(uint32_t))
        stack = <uint32_t *> PyMem_Malloc(_alloc_at_least(atom_count) * sizeof(uint32_t))
        visited = <uint8_t *> PyMem_Malloc(_alloc_at_least(atom_count) * sizeof(uint8_t))
        if comp_of is NULL or stack is NULL or visited is NULL:
            raise MemoryError('query seal scratch allocation failed')
        for s in range(atom_count):
            comp_of[s] = Q_NO_SLOT
        for s in range(atom_count):
            if comp_of[s] != Q_NO_SLOT:
                continue
            comp_of[s] = comp_count
            sp = 0
            stack[sp] = s
            sp += 1
            while sp:
                sp -= 1
                v = stack[sp]
                for k in range(csr_head[v], csr_head[v + 1]):
                    u = csr_nbr[k]
                    if comp_of[u] == Q_NO_SLOT:
                        comp_of[u] = comp_count
                        stack[sp] = u
                        sp += 1
            comp_count += 1

        # --- pass 5: the DFS order -----------------------------------------------------
        order = <uint32_t *> PyMem_Malloc(_alloc_at_least(atom_count) * sizeof(uint32_t))
        back = <uint32_t *> PyMem_Malloc(_alloc_at_least(atom_count) * sizeof(uint32_t))
        back_bond = <uint32_t *> PyMem_Malloc(_alloc_at_least(atom_count) * sizeof(uint32_t))
        pos_of = <uint32_t *> PyMem_Malloc(_alloc_at_least(atom_count) * sizeof(uint32_t))
        is_tree = <uint8_t *> PyMem_Malloc(_alloc_at_least(bond_count) * sizeof(uint8_t))
        if order is NULL or back is NULL or back_bond is NULL or pos_of is NULL or is_tree is NULL:
            raise MemoryError('query seal scratch allocation failed')
        for s in range(atom_count):
            visited[s] = 0
            pos_of[s] = Q_NO_SLOT
        for bslot in range(bond_count):
            is_tree[bslot] = 0

        next_pos = 0
        for comp in range(comp_count):
            # The root is the atom the element index can seed from most cheaply: fewest elements
            # allowed, then -- among atoms equally constrained -- the element the fewest other
            # query atoms demand, then the highest degree, then the lowest slot.
            root = Q_NO_SLOT
            for s in range(atom_count):
                if comp_of[s] != comp:
                    continue
                if root == Q_NO_SLOT:
                    root = s
                    continue
                if card[s] != card[root]:
                    if card[s] < card[root]:
                        root = s
                    continue
                if rarity[s] != rarity[root]:
                    if rarity[s] < rarity[root]:
                        root = s
                    continue
                if csr_head[s + 1] - csr_head[s] > csr_head[root + 1] - csr_head[root]:
                    root = s
            # root cannot still be Q_NO_SLOT: pass 4 labelled every slot and bumped comp_count
            # once per fill, so every label in 0..comp_count-1 owns at least its seed slot, and
            # the loop above takes the first match unconditionally.  Under boundscheck=False a
            # stale Q_NO_SLOT here would be an out-of-bounds write, not an exception.
            visited[root] = 1
            order[next_pos] = root
            back[next_pos] = next_pos
            back_bond[next_pos] = Q_NO_SLOT
            pos_of[root] = next_pos
            next_pos += 1
            sp = 0
            stack[sp] = root
            sp += 1
            while sp:
                v = stack[sp - 1]
                best = Q_NO_SLOT
                best_deg = 0
                best_bond = Q_NO_SLOT
                for k in range(csr_head[v], csr_head[v + 1]):
                    u = csr_nbr[k]
                    if visited[u]:
                        continue
                    deg = csr_head[u + 1] - csr_head[u]
                    if best == Q_NO_SLOT or deg > best_deg or (deg == best_deg and u < best):
                        best = u
                        best_deg = deg
                        best_bond = csr_bnd[k]
                if best == Q_NO_SLOT:
                    sp -= 1
                    continue
                visited[best] = 1
                order[next_pos] = best
                back[next_pos] = pos_of[v]
                back_bond[next_pos] = best_bond
                pos_of[best] = next_pos
                is_tree[best_bond] = 1
                next_pos += 1
                stack[sp] = best
                sp += 1

        # --- group agreement, per component --------------------------------------------
        # A group is per-atom in the journal but per-component in the arena.  This runs BEFORE the
        # fold because _compute_automorphisms needs the per-component group number and has to run
        # before the fold itself -- see the block above that function for why.
        flags = 0
        comp_group_of = <int32_t *> PyMem_Malloc(_alloc_at_least(comp_count) * sizeof(int32_t))
        if comp_group_of is NULL:
            raise MemoryError('query seal scratch allocation failed')
        for comp in range(comp_count):
            comp_group = -1
            for s in range(atom_count):
                if comp_of[s] != comp or group_of[s] < 0:
                    continue
                if comp_group < 0:
                    comp_group = group_of[s]
                elif comp_group != group_of[s]:
                    raise ValueError('component spans two groups')
            comp_group_of[comp] = comp_group
            if comp_group >= 0:
                flags |= QFLAG_HAS_GROUP
        for s in range(atom_count):
            if masked_of[s]:
                flags |= QFLAG_HAS_MASKED

        # --- the automorphism group ----------------------------------------------------
        # Last use of the UNFOLDED atom terms: the fold below overwrites them in place.
        _compute_automorphisms(atom_count, bond_count, atom_terms, bond_terms,
                               csr_head, csr_nbr, csr_bnd, comp_of, comp_group_of, pos_of,
                               &auto_rows, &auto_count, &auto_flags)
        for s in range(atom_count):
            if stereo_of[s]:
                stereo_count += 1
        stereo_count += geom_count
        if stereo_count and auto_count:
            # A STEREO PRIMITIVE IS INVISIBLE TO THE AUTOMORPHISM SEARCH.  Its sign is not in the
            # boxes _wterm_equal compares, by necessity (prim_apply's PRIM_STEREO branch), so a
            # permutation that exchanges two neighbours the query cannot tell apart is accepted as
            # a symmetry even though it inverts the frame the primitive is a statement about --
            # an odd permutation of a centre's directions is exactly a change of sign.
            # mapping_is_canonical DROPS embeddings, so an over-large group can lose a real match,
            # while a group that is too small only over-reports duplicates (its docstring says so).
            # Discarding the rows is the safe direction.  The group is reported PARTIAL, not
            # ASYMMETRIC: it is being under-reported, not proven trivial.
            #
            # TEACHING _wterm_equal ABOUT THE SIGN CANNOT REPLACE THIS, and the reason is worth
            # keeping so nobody tries: in `[C@](F)(F)Cl` the two fluorines carry no sign at all, so
            # their terms compare equal under any sign-aware test, and yet exchanging them is the
            # odd permutation that inverts the CENTRE's frame.  The permutation is illegal because
            # of what it does to a third atom's directions, which is not a property of any pair of
            # terms.  The real fix is to reject candidate permutations that act oddly on a stereo
            # frame -- inside _compute_automorphisms, where the whole permutation is in hand.  That
            # is a performance recovery (it restores the filter for stereo queries), not a
            # correctness fix, and it is parked for the whole-branch review.
            #
            # A GEOMETRY RECORD IS IN THE SAME POSITION, which is why geom_count is in the count
            # above: a direction is in no box either, and exchanging a terminal's two substituents
            # turns the geometry over.
            auto_count = 0
            auto_flags = QFLAG_PARTIAL_AUTOMORPHISM
        flags |= auto_flags

        # --- pass 6: fold tree bonds, close the rest -----------------------------------
        clo_owner = <uint32_t *> PyMem_Malloc(_alloc_at_least(bond_count) * sizeof(uint32_t))
        clo_to = <uint32_t *> PyMem_Malloc(_alloc_at_least(bond_count) * sizeof(uint32_t))
        clo_bond = <uint32_t *> PyMem_Malloc(_alloc_at_least(bond_count) * sizeof(uint32_t))
        if clo_owner is NULL or clo_to is NULL or clo_bond is NULL:
            raise MemoryError('query seal scratch allocation failed')

        for i in range(atom_count):
            if back_bond[i] == Q_NO_SLOT:
                continue
            s = order[i]
            bslot = back_bond[i]
            try:
                _fold_bond_into_atom(&folded, &atom_terms[s], &bond_terms[bslot])
            except ValueError as exc:
                raise ValueError('bond %d-%d: %s' % (stable_of[bond_a[bslot]],
                                                     stable_of[bond_b[bslot]], exc))
            atom_terms[s] = folded

        for bslot in range(bond_count):
            if is_tree[bslot]:
                continue
            pa = pos_of[bond_a[bslot]]
            pb = pos_of[bond_b[bslot]]
            if pa > pb:
                clo_owner[closure_count] = pa
                clo_to[closure_count] = pb
            else:
                clo_owner[closure_count] = pb
                clo_to[closure_count] = pa
            clo_bond[closure_count] = bslot
            closure_count += 1

        # --- signature: the screen's demand words (Task 14) ----------------------------
        # Walk each atom's (possibly folded) disjunction and collect bits that EVERY box
        # of that atom requires.  A bit is "required" by one box when the box forbids every
        # other bit in its span, i.e. (M & ~neg[w]) is exactly one bit.  Per-atom results
        # are OR'd into the four query-signature words.
        #
        # Element bits use _element_from_neg rather than the general span loop:
        # the element span straddles words 0 and 1, so the "single remaining bit" check
        # would have to span two words; _element_from_neg already handles that correctly.
        #
        # Non-element spans: SPAN_WORD/SPAN_MASK enumerates every one-hot span in words 0-3.
        # The topology and order bits from a folded tree bond appear in neg[0] of a non-root
        # atom and are therefore sound: fill_features ORs every incident bond's bits into the
        # atom's aggregate word 0, so a match's atom carries the bit and the union row too.
        memset(sig, 0, 4 * sizeof(uint64_t))
        for s in range(atom_count):
            wt = &atom_terms[s]
            memset(atom_sig, 0, 4 * sizeof(uint64_t))

            # --- element contribution ---------------------------------------------------
            # Every box must decode via _element_from_neg to the same definite element.
            # [C,N] has two boxes that each allow one element but disagree, so it contributes
            # nothing -- and that is correct: no single atom is guaranteed to be carbon.
            e_common = -2   # -2 = not yet seen, -1 = no universal requirement
            for k in range(wt.count):
                wbx = &wt.boxes[k]
                e_box = _element_from_neg(wbx.neg[0], wbx.neg[1])
                if e_box <= 0:
                    e_common = -1
                    break
                if e_common == -2:
                    e_common = e_box
                elif e_common != e_box:
                    e_common = -1
                    break
            if e_common > 0:
                if e_common <= 56:
                    atom_sig[0] |= <uint64_t> 1 << (57 - <uint32_t> e_common)
                else:
                    atom_sig[0] |= <uint64_t> 1               # heavy-element marker (bit 0)
                    atom_sig[1] |= <uint64_t> 1 << (<uint32_t> e_common - 57)

            # --- non-element one-hot span contributions ---------------------------------
            # For each span, check whether every box pins the span to the same single bit.
            for si in range(SPAN_COUNT):
                w = SPAN_WORD[si]
                M = SPAN_MASK[si]
                span_contrib = 0
                for k in range(wt.count):
                    allowed_bits = M & ~wt.boxes[k].neg[w]
                    # A box pins the span to exactly one bit iff allowed_bits has exactly
                    # one bit set.  A zero means the box is unsatisfiable (already pruned),
                    # and multiple bits mean the box does not pin this span.
                    if not allowed_bits or (allowed_bits & (allowed_bits - 1)):
                        span_contrib = 0
                        break
                    if k == 0:
                        span_contrib = allowed_bits   # first box's required bit
                    elif span_contrib != allowed_bits:
                        span_contrib = 0              # boxes disagree on this span
                        break
                atom_sig[w] |= span_contrib

            # OR this atom's contribution into the query signature.
            for k in range(4):
                sig[k] |= atom_sig[k]

        # --- pass 7: emit ---------------------------------------------------------------
        box_total = 0
        any_total = 0
        for s in range(atom_count):
            box_total += atom_terms[s].count
            any_total += _term_any_total(&atom_terms[s])
        bond_box_total = 0
        for bslot in range(bond_count):
            if not is_tree[bslot]:
                bond_box_total += bond_terms[bslot].count
                any_total += _term_any_total(&bond_terms[bslot])

        # Count demanded elements for the compact demand list (QSEG_DEMAND_LIST).
        demand_list_count = 0
        for i in range(1, 119):
            if demand[i]:
                demand_list_count += 1

        q = query_alloc(atom_count, bond_count, box_total, bond_box_total, any_total,
                        closure_count, comp_count, auto_count, demand_list_count, stereo_count)
        # Take arena pointers only after the last allocation.  query_alloc is the only one, and it
        # stays that way: the arena is immutable after seal.
        if auto_count:
            memcpy(query_automorphisms(q), auto_rows,
                   <size_t> auto_count * atom_count * sizeof(uint32_t))
        out_atoms = q.atoms()
        out_boxes = query_boxes(q)
        out_any = query_any(q)
        out_bonds = query_bonds(q)
        out_bond_boxes = query_bond_boxes(q)
        out_closures = query_closures(q)
        out_comps = query_components(q)
        out_demand = query_element_demand(q)
        out_demand_list = query_demand_list(q)
        out_stereo = query_stereo(q)

        box_cursor = 0
        any_cursor = 0
        clo_cursor = 0
        for i in range(atom_count):
            s = order[i]
            qa = &out_atoms[i]
            qa.box_begin = box_cursor
            qa.box_count = <uint16_t> atom_terms[s].count
            qa.back = back[i]
            qa.flags = wild_of[s]
            if masked_of[s]:
                qa.flags |= QATOM_MASKED
            if back_bond[i] == Q_NO_SLOT:
                qa.flags |= QATOM_ROOT
            # Both bits set means the atom named both signs; the frame and the sign itself live in
            # this position's QSEG_STEREO record, and these bits are the cheap per-atom witness that
            # one exists.
            if stereo_of[s] & 1:
                qa.flags |= QATOM_STEREO_CW
            if stereo_of[s] & 2:
                qa.flags |= QATOM_STEREO_CCW
            qa.map_number = map_of[s]
            qa.spare = 0            # still reserved headroom: the readiness position went into
                                    # QSEG_STEREO, where the frame it belongs to already is
            box_cursor = _emit_boxes(out_boxes, box_cursor, out_any, &any_cursor, &atom_terms[s])
            # closures owned by this position, in bond-slot order
            qa.closure_begin = clo_cursor
            qa.closure_count = 0
            for k in range(closure_count):
                if clo_owner[k] == i:
                    out_closures[clo_cursor].to_index = clo_to[k]
                    out_closures[clo_cursor].bond_index = clo_bond[k]
                    clo_cursor += 1
                    qa.closure_count += 1

        # --- QSEG_STEREO: the frame each stereo primitive is a statement about ------------
        stereo_cursor = 0
        for i in range(atom_count):
            s = order[i]
            if not stereo_of[s]:
                continue
            qs = &out_stereo[stereo_cursor]
            stereo_cursor += 1
            qs.position = i
            qs.sign = stereo_of[s]
            qs.spare = 0
            deg = csr_head[s + 1] - csr_head[s]
            qs.n_refs = <uint8_t> (deg if deg < 255 else 255)
            for k in range(4):
                qs.refs[k] = Q_NO_SLOT
            ready = i
            if deg <= 4:
                # THE QUERY'S OWN RULING-F26 ORDER: named neighbours by ascending query slot --
                # which is the order the caller created them in -- and then the unnamed direction,
                # which needs no entry.  The CSR is built in bond-slot order, so the ascent is
                # imposed here rather than inherited.  This order is what makes '@' a statement in
                # the QUERY's frame and therefore independent of the target: the kernel re-expresses
                # it in the target unit's order (translate_parity) instead of comparing raw values.
                # Every hardcoded 1 / 2 in a stereo test depends on this being the creation order.
                n_nbrs = 0
                for k in range(csr_head[s], csr_head[s + 1]):
                    nbrs[n_nbrs] = csr_nbr[k]
                    n_nbrs += 1
                for k in range(n_nbrs):
                    for u in range(k + 1, n_nbrs):
                        if nbrs[u] < nbrs[k]:
                            swap = nbrs[k]
                            nbrs[k] = nbrs[u]
                            nbrs[u] = swap
                for k in range(n_nbrs):
                    qs.refs[k] = pos_of[nbrs[k]]
                    if pos_of[nbrs[k]] > ready:
                        ready = pos_of[nbrs[k]]
            # The DFS position at which the last of this frame's directions becomes mapped.  The
            # anchor's own position is in the maximum because a lone stereo atom in a one-atom
            # component still has to be tested somewhere, and because nothing guarantees the
            # neighbours come later: the plan is ordered by rarity, not by adjacency.
            qs.readiness = ready

        # A GEOMETRY IS A RECORD ABOUT TWO ATOMS, and the fields say so by holding four positions
        # rather than an anchor and its directions: `position` is one terminal, `refs` is that
        # terminal's marked substituent, the other terminal's, and the other terminal.  `spare` is the
        # low byte SU_CIS_TRANS -- SU_TETRA being 0, which is what every record above writes -- and the
        # high byte the parity demanded in the frame `(marked, other, marked, other)`, 1 for trans.
        # `sign` stays 0: nothing here came off a box, so there is no per-box sign to report.
        for k in range(geom_count):
            qs = &out_stereo[stereo_cursor]
            stereo_cursor += 1
            qs.position = pos_of[geoms[k].term_a]
            qs.sign = 0
            qs.n_refs = 3
            qs.spare = <uint16_t> (SU_CIS_TRANS | ((2 - <int> geoms[k].trans) << 8))
            qs.refs[0] = pos_of[geoms[k].mark_a]
            qs.refs[1] = pos_of[geoms[k].mark_b]
            qs.refs[2] = pos_of[geoms[k].term_b]
            qs.refs[3] = Q_NO_SLOT
            ready = qs.position
            for u in range(3):
                if qs.refs[u] > ready:
                    ready = qs.refs[u]
            qs.readiness = ready

        bond_box_cursor = 0
        for bslot in range(bond_count):
            qb = &out_bonds[bslot]
            qb.flags = 0
            qb.box_begin = bond_box_cursor
            if is_tree[bslot]:
                # A tree bond's boxes now live in the atom it leads to.  The record stays: it is
                # the home of the stereo tri-state, and Task 16's differential test needs to see
                # that the bond exists at all.
                qb.box_count = 0
                continue
            qb.box_count = <uint16_t> bond_terms[bslot].count
            bond_box_cursor = _emit_boxes(out_bond_boxes, bond_box_cursor, out_any, &any_cursor,
                                         &bond_terms[bslot])

        # Components, in the order their positions were emitted.  The DFS finishes one component
        # before starting the next, so each component's positions are a contiguous run.
        i = 0
        for comp in range(comp_count):
            qc = &out_comps[comp]
            qc.begin = i
            while i < atom_count and comp_of[order[i]] == comp:
                i += 1
            qc.end = i
            qc.group = comp_group_of[comp]

        for i in range(120):
            out_demand[i] = demand[i]

        # Write the compact demand list (QSEG_DEMAND_LIST): one (element, count) pair per
        # demanded element.  query_may_match iterates this list directly, skipping the 118-slot
        # scan over the full histogram.
        dl_cursor = 0
        for i in range(1, 119):
            if demand[i]:
                out_demand_list[dl_cursor] = i           # element number
                out_demand_list[dl_cursor + 1] = demand[i]  # required count
                dl_cursor += 2

        # Write the pre-computed signature (Task 14).
        for k in range(4):
            q.header.signature[k] = sig[k]

        # QFLAG_HAS_STEREO is live as of Task 10: it is set exactly when QSEG_STEREO is non-empty,
        # and matcher_init reads it to decide whether to build the target's stereo unit table at
        # all.  A query with no stereo primitive must pay nothing for stereo, which is what that
        # equivalence buys.  The scan stays a scan over the atom flags rather than a
        # `if stereo_count` so that the flag and the per-atom bits cannot disagree.
        # There is deliberately no bond half: qbond_t.flags has no stereo bits defined, and a
        # direction is not a demand on a bond anyway -- the geometry it is half of is a record in
        # QSEG_STEREO, counted here so the kernel builds the target's unit table for it.
        for i in range(atom_count):
            if out_atoms[i].flags & (QATOM_STEREO_CW | QATOM_STEREO_CCW):
                flags |= QFLAG_HAS_STEREO
        if geom_count:
            flags |= QFLAG_HAS_STEREO
        q.header.flags = flags

        ids = <uint32_t *> PyMem_Malloc(_alloc_at_least(atom_count) * sizeof(uint32_t))
        if ids is NULL:
            raise MemoryError('query seal allocation failed')
        for i in range(atom_count):
            ids[i] = stable_of[order[i]]
        position_to_n[0] = ids
    finally:
        PyMem_Free(slot_of)
        PyMem_Free(stable_of)
        PyMem_Free(masked_of)
        PyMem_Free(map_of)
        PyMem_Free(group_of)
        PyMem_Free(bond_a)
        PyMem_Free(bond_b)
        PyMem_Free(atom_terms)
        PyMem_Free(bond_terms)
        PyMem_Free(scratch)
        PyMem_Free(csr_head)
        PyMem_Free(csr_cur)
        PyMem_Free(csr_nbr)
        PyMem_Free(csr_bnd)
        PyMem_Free(comp_of)
        PyMem_Free(order)
        PyMem_Free(back)
        PyMem_Free(back_bond)
        PyMem_Free(pos_of)
        PyMem_Free(stack)
        PyMem_Free(visited)
        PyMem_Free(card)
        PyMem_Free(rarity)
        PyMem_Free(wild_of)
        PyMem_Free(is_tree)
        PyMem_Free(clo_owner)
        PyMem_Free(clo_to)
        PyMem_Free(clo_bond)
        PyMem_Free(comp_group_of)
        PyMem_Free(auto_rows)
        PyMem_Free(stereo_of)
        PyMem_Free(geoms)
    return q


def _seal_probe(list ops):
    """Python test probe: seal a journal of op tuples and read the arena back as plain data.

    Nothing here needs to be fast; clarity wins, since this is the only window into seal the
    tests get.
    """
    cdef qop_t *buf = NULL
    cdef uint32_t *ids = NULL
    cdef Query q = None
    cdef uint32_t op_count = <uint32_t> len(ops)
    cdef uint32_t next_id = 1
    cdef uint32_t i, k, p, ai
    cdef qop_t *rec
    cdef qatom_t *atoms
    cdef qatom_t *qa
    cdef qbox_t *boxes
    cdef qbox_t *bx
    cdef qany_t *any_records
    cdef qany_t *an
    cdef qbond_t *bonds
    cdef qbond_t *qb
    cdef qbox_t *bond_boxes
    cdef qclosure_t *closures
    cdef qclosure_t *cl
    cdef qcomp_t *comps
    cdef qcomp_t *qc
    cdef uint32_t *demand
    cdef tuple op
    cdef dict slot_of_number = {}
    cdef list order_out = [], back_out = [], roots_out = [], comps_out = [], closures_out = []
    cdef list boxes_out = [], bond_box_counts = [], bond_boxes_out = []
    cdef list masked_out = [], map_out = [], demand_out = [], autos_out = [], stereo_out = []
    cdef list wildcards_out = []
    cdef list box_list, any_list
    cdef uint32_t *autos
    cdef qstereo_t *stereo
    cdef qstereo_t *qs

    for op in ops:
        if op[0] == 'bond' or op[0] == 'btoken' or op[0] == 'bop':
            if <uint32_t> op[1] >= next_id:
                next_id = <uint32_t> op[1] + 1
            if <uint32_t> op[2] >= next_id:
                next_id = <uint32_t> op[2] + 1
        else:
            if <uint32_t> op[1] >= next_id:
                next_id = <uint32_t> op[1] + 1

    buf = <qop_t *> PyMem_Malloc(_alloc_at_least(op_count) * sizeof(qop_t))
    if buf is NULL:
        raise MemoryError()
    try:
        for i in range(op_count):
            op = ops[i]
            rec = &buf[i]
            memset(rec, 0, sizeof(qop_t))
            if op[0] == 'atom':
                rec.op = QOP_ADD_ATOM
                rec.a = op[1]
                # slot i is the i-th ('atom', n) tuple: query_seal hands back position ->
                # stable id, and every test reads slots, so the probe derives them here.
                slot_of_number[op[1]] = len(slot_of_number)
            elif op[0] == 'token':
                rec.op = QOP_ATOM_TOKEN
                rec.a = op[1]
                rec.opcode = OPC_PRIM
                rec.kind = PRIM_NAMES[op[2]]
                rec.value = op[3]
                rec.negated = 1 if op[4] else 0
            elif op[0] == 'op':
                rec.op = QOP_ATOM_TOKEN
                rec.a = op[1]
                rec.opcode = OPC_NAMES[op[2]]
            elif op[0] == 'bond':
                rec.op = QOP_ADD_BOND
                rec.a = op[1]
                rec.b = op[2]
            elif op[0] == 'btoken':
                rec.op = QOP_BOND_TOKEN
                rec.a = op[1]
                rec.b = op[2]
                rec.opcode = OPC_PRIM
                rec.kind = PRIM_NAMES[op[3]]
                rec.value = op[4]
                rec.negated = 1 if op[5] else 0
            elif op[0] == 'bop':
                rec.op = QOP_BOND_TOKEN
                rec.a = op[1]
                rec.b = op[2]
                rec.opcode = OPC_NAMES[op[3]]
            elif op[0] == 'group':
                rec.op = QOP_SET_GROUP
                rec.a = op[1]
                rec.value = op[2]
            elif op[0] == 'masked':
                rec.op = QOP_SET_MASKED
                rec.a = op[1]
            elif op[0] == 'map':
                rec.op = QOP_SET_MAP
                rec.a = op[1]
                rec.value = op[2]
            else:
                raise ValueError('unknown journal op %r' % (op[0],))
        q = query_seal(buf, op_count, next_id, &ids)
    finally:
        PyMem_Free(buf)

    atoms = q.atoms()
    boxes = query_boxes(q)
    any_records = query_any(q)
    bonds = query_bonds(q)
    bond_boxes = query_bond_boxes(q)
    closures = query_closures(q)
    comps = query_components(q)
    demand = query_element_demand(q)
    try:
        for p in range(q.header.atom_count):
            qa = &atoms[p]
            order_out.append(slot_of_number[ids[p]])
            back_out.append(qa.back)
            masked_out.append(bool(qa.flags & QATOM_MASKED))
            if qa.flags & QATOM_ANY_ELEMENT:
                wildcards_out.append('any')
            elif qa.flags & QATOM_METAL_ELEMENT:
                wildcards_out.append('metal')
            else:
                wildcards_out.append(None)
            map_out.append(qa.map_number)
            if qa.flags & QATOM_ROOT:
                roots_out.append(p)
            box_list = []
            for k in range(qa.box_count):
                bx = &boxes[qa.box_begin + k]
                any_list = []
                for ai in range(bx.any_count):
                    an = &any_records[bx.any_begin + ai]
                    any_list.append((an.word, an.mask))
                box_list.append({'neg': (bx.neg[0], bx.neg[1], bx.neg[2], bx.neg[3]),
                                 'any': tuple(any_list), 'sign': bx.sign})
            boxes_out.append(box_list)
            for k in range(qa.closure_count):
                cl = &closures[qa.closure_begin + k]
                closures_out.append((cl.to_index, cl.bond_index))
    finally:
        PyMem_Free(ids)

    for i in range(q.header.bond_count):
        qb = &bonds[i]
        bond_box_counts.append(qb.box_count)
        box_list = []
        for k in range(qb.box_count):
            bx = &bond_boxes[qb.box_begin + k]
            box_list.append({'neg': (bx.neg[0], bx.neg[1], bx.neg[2], bx.neg[3])})
        bond_boxes_out.append(box_list)
    for i in range(q.header.component_count):
        qc = &comps[i]
        comps_out.append((qc.begin, qc.end, qc.group))
    for i in range(120):
        demand_out.append(demand[i])
    autos = query_automorphisms(q)
    for i in range(q.header.automorphism_count):
        box_list = []
        for k in range(q.header.atom_count):
            box_list.append(autos[i * q.header.atom_count + k])
        autos_out.append(tuple(box_list))
    stereo = query_stereo(q)
    for i in range(q.header.stereo_count):
        qs = &stereo[i]
        stereo_out.append({'position': qs.position, 'readiness': qs.readiness,
                           'refs': (qs.refs[0], qs.refs[1], qs.refs[2], qs.refs[3]),
                           'sign': qs.sign, 'n_refs': qs.n_refs,
                           'kind': qs.spare & 0xFF, 'demand': qs.spare >> 8})

    return {'stereo': stereo_out,
            'atom_count': q.header.atom_count, 'bond_count': q.header.bond_count,
            'component_count': q.header.component_count,
            'closure_count': <uint32_t> len(closures_out), 'flags': q.header.flags,
            'order': order_out, 'back': back_out, 'roots': roots_out, 'components': comps_out,
            'closures': closures_out, 'boxes': boxes_out,
            'bond_box_counts': bond_box_counts, 'bond_boxes': bond_boxes_out,
            'masked': masked_out, 'map_numbers': map_out, 'element_demand': demand_out,
            'wildcards': wildcards_out,
            'automorphism_count': q.header.automorphism_count, 'automorphisms': autos_out}
