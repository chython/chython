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
# CIP DESCRIPTOR ASSIGNMENT.  `_stereo.pxi` says which sites can carry a configuration and stores the
# parity; this fragment says which letter the configuration is called.
#
# RULES.md 1.5 governs and is not restated here.  Three of its consequences shape every function below:
# a stored descriptor is never read as an input, the aromatic form is ranked as it stands rather than
# kekulised, and a site this fragment cannot decide is refused rather than guessed.
#
# Phase 1 decides rules 1a, 1b and 2 on SU_TETRA.  A site those three tie is `cip:undecided`.

# A ring system is the connected component of the ring-bond subgraph.  It bounds the visited-mask half
# of a digraph node's key: outside a ring system a branch cannot depend on the path that reached it, so
# a wider mask would be untested width.  Union-find over ring bonds, the same shape as the ring-system
# count in `_descriptors.pxi`.
cdef inline uint32_t _cip_uf_find(uint32_t *parent, uint32_t i) noexcept nogil:
    while parent[i] != i:
        parent[i] = parent[parent[i]]
        i = parent[i]
    return i


cdef int cip_ring_systems(Structure structure, uint32_t *out) except -1:
    """1-based ring-system id per atom slot, 0 for an atom on no ring.  Returns the system count.

    `out` is caller-supplied and atom_count wide.  `HE_IN_RING` is what makes this a scan rather than a
    perception: `mark_bridges` set the bit on seal, so a bridge -- the bond between bicyclohexyl's two
    rings -- fuses nothing.
    """
    cdef uint32_t n = structure.header.atom_count
    cdef uint32_t *ptr = csr_ptr(structure)
    cdef halfedge_t *edges = csr_edges(structure)
    cdef uint32_t i, k, ra, rb, systems = 0
    cdef uint32_t *parent = <uint32_t *> malloc(<size_t> (n if n else 1) * sizeof(uint32_t))
    if parent is NULL:
        raise MemoryError()
    try:
        for i in range(n):
            parent[i] = i
            out[i] = 0
        for i in range(n):
            for k in range(ptr[i], ptr[i + 1]):
                if edges[k].flags & HE_IN_RING:
                    ra = _cip_uf_find(parent, i)
                    rb = _cip_uf_find(parent, edges[k].to)
                    if ra != rb:
                        parent[ra] = rb
        # dense 1-based ids, assigned in first-appearance order of each root
        for i in range(n):
            for k in range(ptr[i], ptr[i + 1]):
                if edges[k].flags & HE_IN_RING:
                    ra = _cip_uf_find(parent, i)
                    if out[ra] == 0:
                        systems += 1
                        out[ra] = systems
                    out[i] = out[ra]
                    break
    finally:
        free(parent)
    return <int> systems


# THE RANK WORD.  One uint64_t per digraph node, compared as an integer so a rule-1 step is a word
# compare rather than a struct compare.  Higher word ranks higher.  The field widths are the domain and
# are stated once, here (RULES.md 6.1).
DEF CIP_Z_SHIFT     = 48     # bits 63-48: atomic number x CIP_Z_SCALE
DEF CIP_Z_SCALE     = 2      # a phase-2 mancude duplicate carries the MEAN of its positions, and x2
                             # makes every mean this fragment can meet an integer in the same field
DEF CIP_REAL_BIT    = 47     # bit 47: 1 real atom, 0 duplicate -- rule 1b, above the back count
DEF CIP_BACK_SHIFT  = 32     # bits 46-32: spheres from the duplicate UP TO the atom it duplicates
DEF CIP_BACK_MAX    = 0x7fff
DEF CIP_MASS_SCALE  = 1000   # bits 31-0: mass in millidaltons -- rule 2


cdef inline uint64_t cip_word_of(uint32_t z, uint32_t isotope, bint duplicate,
                                 uint32_t back) noexcept nogil:
    """The rule-1a/1b/2 key of one digraph node.  See the DEF block above for the layout.

    RULE 1B IS ENCODED RELATIVE, NOT ABSOLUTE, and that choice is what keeps the digraph shared.  A
    duplicate ranks by the root distance of the atom it duplicates -- nearer the root ranks higher -- and
    `back` is how many spheres UP the path that atom sits from the duplicate itself: 0 for the duplicate
    a double bond adds (it duplicates its own neighbour), 2 for the duplicate of the parent across a
    multiple bond, and the length of the jump for a ring closure.  Two nodes are only ever compared at
    equal sphere (the comparator advances both frontiers in lockstep), where a larger `back` is exactly a
    smaller absolute distance, so the ordering is the absolute one.  Storing the absolute distance
    instead would make one cyclohexyl branch intern differently at every depth it hangs from.

    A REAL atom passes back 0, and its bit 47 puts it above every duplicate of the same element
    regardless of the field.  `word == 0` is reserved for the phantom the comparator pads a short child
    group with, which is why z 0 never reaches here.

    `element_mass` takes isotope 0 to mean the abundance-weighted average, which is what an atom with no
    stated isotope carries, so rule 2 needs no special case for it.  It answers 0 for an element with no
    measured mass in `isotopes.tsv`; rule 2 is then blind between two such atoms of one element, which
    is the correct outcome -- there is no mass to rank them by -- and the rule 1a/1b fields above are
    unaffected.
    """
    cdef uint32_t b = back if back < CIP_BACK_MAX else CIP_BACK_MAX
    cdef uint64_t w = (<uint64_t> z * CIP_Z_SCALE) << CIP_Z_SHIFT
    if not duplicate:
        w |= (<uint64_t> 1) << CIP_REAL_BIT
    w |= (<uint64_t> b) << CIP_BACK_SHIFT
    w |= <uint64_t> <uint32_t> round(element_mass(z, isotope) * CIP_MASS_SCALE)
    return w


cdef inline uint64_t cip_rank_word(atom_t *a, bint duplicate, uint32_t back) noexcept nogil:
    return cip_word_of(a.element, a.isotope, duplicate, back)


# THE SHARED DIGRAPH.  CIP's hierarchical digraph is a tree over simple paths, which is what makes it
# exponential and per-centre.  Two things make it one DAG per molecule instead: nodes are hash-consed on
# their FULL key, so identical branches are one node; and the key is depth-invariant (see `cip_word_of`),
# so a branch hanging further out is still the same node.
DEF CIP_NO_ID = 0xffffffff
DEF CIP_NODE_INIT = 256      # node pool, child pool and intern table all start from this and double
DEF CIP_DEPTH_MAX = 1024     # a branch deeper than this refuses rather than recursing on.  Depth is
                             # bounded by the longest simple path, so this is a stack guard and not a
                             # chemistry bound: no shipped corpus reaches it
DEF CIP_KIDS_MAX = 8         # children of one node: real neighbours, the duplicates their orders add,
                             # the parent's duplicates and the implicit hydrogens.  A sulfate sulfur
                             # entered from a double-bonded oxygen reaches 5; above 8 the site refuses
DEF CIP_NODE_UNDECIDABLE = 1 # order 4 or 8 in this node's bonds, an unknown hydrogen count, or the
                             # same anywhere below it
DEF CIP_MEMO_INIT = 1024     # comparison verdicts; a steroid's centres re-ask the same pairs
DEF CIP_MEMO_ABSENT = 2      # not a verdict: the three verdicts are -1, 0 and 1

# All bits set, and a `cdef` rather than a `DEF` because a 64-bit literal of this width is a Python int
# to Cython, which cannot be compared under nogil.  No node pair can spell it: both ids would have to be
# CIP_NO_ID, and CIP_NO_ID is never a node.
cdef uint64_t CIP_MEMO_EMPTY = <uint64_t> -1


cdef struct cip_node_t:
    uint64_t word            # cip_word_of: rules 1a, 1b, 2
    uint32_t first_child     # index into ctx.child; equals ctx.n_child when n_children is 0
    uint32_t n_children
    uint8_t flags            # CIP_NODE_UNDECIDABLE.  PART OF THE INTERN KEY
    uint8_t ordered          # ctx.corder[first_child ...] holds the CIP priority order


# The context is C pointers only, never the `Structure` object: a struct holding an extension type would
# need refcounting this fragment has no way to do, and `pinned_search_t` next door is the precedent.  One
# context per molecule, shared by every unit in it -- that sharing is the whole point.
cdef struct cip_ctx_t:
    uint32_t n               # atom count
    uint32_t *ptr            # csr_ptr
    halfedge_t *edges        # csr_edges
    atom_t *atoms
    uint32_t *ring_system    # per atom, cip_ring_systems
    uint8_t *on_path         # per atom: 1 while the walk holds it, for expansion clause 4
    uint32_t *path_depth     # per atom: the sphere it sits at on the current path, for rule 1b's back
    cip_node_t *node         # the interned node pool
    uint32_t *child          # child-id runs in the canonical intern order, indexed by node.first_child
    uint32_t *corder         # the same runs in CIP priority order, filled lazily
    uint32_t *bucket         # open-addressed intern table, CIP_NO_ID for empty
    uint32_t mask            # bucket count - 1; always a power of two minus one
    uint32_t n_nodes
    uint32_t cap_nodes
    uint32_t n_child
    uint32_t cap_child
    uint64_t *memo_key       # comparison verdicts, keyed (a << 32) | b
    int8_t *memo_val
    uint32_t memo_mask
    uint32_t n_memo


cdef void cip_ctx_free(cip_ctx_t *ctx) noexcept nogil:
    """Safe on a partly built context, so `cip_ctx_init` has one failure path rather than six."""
    free(ctx.ring_system)
    free(ctx.on_path)
    free(ctx.path_depth)
    free(ctx.node)
    free(ctx.child)
    free(ctx.corder)
    free(ctx.bucket)
    free(ctx.memo_key)
    free(ctx.memo_val)
    ctx.memo_key = NULL
    ctx.memo_val = NULL
    ctx.ring_system = NULL
    ctx.on_path = NULL
    ctx.path_depth = NULL
    ctx.node = NULL
    ctx.child = NULL
    ctx.corder = NULL
    ctx.bucket = NULL


cdef int cip_ctx_init(cip_ctx_t *ctx, Structure structure) except -1:
    """One context per molecule: the pointers the walk reads, and the pools it fills."""
    cdef uint32_t i
    cdef uint32_t n = structure.header.atom_count
    cdef size_t span = <size_t> (n if n else 1)
    ctx.n = n
    ctx.ptr = csr_ptr(structure)
    ctx.edges = csr_edges(structure)
    ctx.atoms = structure.atoms()
    ctx.n_nodes = 0
    ctx.n_child = 0
    ctx.cap_nodes = CIP_NODE_INIT
    ctx.cap_child = CIP_NODE_INIT * 4
    ctx.mask = CIP_NODE_INIT * 2 - 1
    ctx.ring_system = <uint32_t *> malloc(span * sizeof(uint32_t))
    ctx.on_path = <uint8_t *> calloc(span, sizeof(uint8_t))
    ctx.path_depth = <uint32_t *> calloc(span, sizeof(uint32_t))
    ctx.node = <cip_node_t *> malloc(<size_t> ctx.cap_nodes * sizeof(cip_node_t))
    ctx.child = <uint32_t *> malloc(<size_t> ctx.cap_child * sizeof(uint32_t))
    ctx.corder = <uint32_t *> malloc(<size_t> ctx.cap_child * sizeof(uint32_t))
    ctx.bucket = <uint32_t *> malloc(<size_t> (ctx.mask + 1) * sizeof(uint32_t))
    ctx.n_memo = 0
    ctx.memo_mask = CIP_MEMO_INIT - 1
    ctx.memo_key = <uint64_t *> malloc(<size_t> CIP_MEMO_INIT * sizeof(uint64_t))
    ctx.memo_val = <int8_t *> malloc(<size_t> CIP_MEMO_INIT * sizeof(int8_t))
    if (ctx.ring_system is NULL or ctx.on_path is NULL or ctx.path_depth is NULL or ctx.node is NULL
            or ctx.child is NULL or ctx.corder is NULL or ctx.bucket is NULL
            or ctx.memo_key is NULL or ctx.memo_val is NULL):
        cip_ctx_free(ctx)
        raise MemoryError()
    for i in range(ctx.mask + 1):
        ctx.bucket[i] = CIP_NO_ID
    for i in range(CIP_MEMO_INIT):
        ctx.memo_key[i] = CIP_MEMO_EMPTY
    cip_ring_systems(structure, ctx.ring_system)
    return 0


# A hash only picks a bucket; every candidate is compared on the FULL key, so a collision costs one probe
# and can never hand back a wrong id.  That is what lets the comparator treat `id_a == id_b` as a proof of
# structural identity rather than as a gamble on a digest.
cdef inline uint64_t _cip_hash(uint64_t word, uint8_t flags, uint32_t *kids,
                               uint32_t n_kids) noexcept nogil:
    """FNV-1a over the whole key."""
    cdef uint64_t h = <uint64_t> 0xcbf29ce484222325
    cdef uint64_t prime = <uint64_t> 0x100000001b3
    cdef uint32_t i, b
    for b in range(8):
        h = (h ^ ((word >> (8 * b)) & 0xff)) * prime
    h = (h ^ <uint64_t> flags) * prime
    for i in range(n_kids):
        h = (h ^ <uint64_t> kids[i]) * prime
    return h


cdef int _cip_reserve(cip_ctx_t *ctx, uint32_t n_kids) except -1:
    """Room for one more node and its children, and a bucket table under 70% load.

    CALLED BEFORE THE PROBE, not after it: a rehash moves every bucket, so a slot computed first and
    used after a growth would write the new node into the wrong place.
    """
    cdef uint32_t i, id_, slot, cap
    cdef cip_node_t *node
    cdef uint32_t *child
    cdef uint32_t *bucket
    if ctx.n_nodes + 1 > ctx.cap_nodes:
        cap = ctx.cap_nodes * 2
        node = <cip_node_t *> realloc(ctx.node, <size_t> cap * sizeof(cip_node_t))
        if node is NULL:
            raise MemoryError()
        ctx.node = node
        ctx.cap_nodes = cap
    if ctx.n_child + n_kids > ctx.cap_child:
        cap = ctx.cap_child * 2
        while cap < ctx.n_child + n_kids:
            cap *= 2
        child = <uint32_t *> realloc(ctx.child, <size_t> cap * sizeof(uint32_t))
        if child is NULL:
            raise MemoryError()
        ctx.child = child
        child = <uint32_t *> realloc(ctx.corder, <size_t> cap * sizeof(uint32_t))
        if child is NULL:
            raise MemoryError()
        ctx.corder = child
        ctx.cap_child = cap
    if (ctx.n_nodes + 1) * 10 > (ctx.mask + 1) * 7:
        cap = (ctx.mask + 1) * 2
        bucket = <uint32_t *> malloc(<size_t> cap * sizeof(uint32_t))
        if bucket is NULL:
            raise MemoryError()
        for i in range(cap):
            bucket[i] = CIP_NO_ID
        free(ctx.bucket)
        ctx.bucket = bucket
        ctx.mask = cap - 1
        for id_ in range(ctx.n_nodes):
            slot = <uint32_t> (_cip_hash(ctx.node[id_].word, ctx.node[id_].flags,
                                         &ctx.child[ctx.node[id_].first_child],
                                         ctx.node[id_].n_children) & ctx.mask)
            while ctx.bucket[slot] != CIP_NO_ID:
                slot = (slot + 1) & ctx.mask
            ctx.bucket[slot] = id_
    return 0


cdef uint32_t cip_intern(cip_ctx_t *ctx, uint64_t word, uint8_t flags, uint32_t *kids,
                         uint32_t n_kids) except? CIP_NO_ID:
    """The id of the node `(word, flags, kids)`, created once and returned to every later request.

    `kids` must already be in the canonical order `_cip_sort_ids_desc` fixes, or two spellings of one
    node would intern twice.
    """
    cdef uint64_t h
    cdef uint32_t slot, id_, i
    cdef bint same
    _cip_reserve(ctx, n_kids)
    h = _cip_hash(word, flags, kids, n_kids)
    slot = <uint32_t> (h & ctx.mask)
    while True:
        id_ = ctx.bucket[slot]
        if id_ == CIP_NO_ID:
            break
        if (ctx.node[id_].word == word and ctx.node[id_].flags == flags and
                ctx.node[id_].n_children == n_kids):
            same = True
            for i in range(n_kids):
                if ctx.child[ctx.node[id_].first_child + i] != kids[i]:
                    same = False
                    break
            if same:
                return id_
        slot = (slot + 1) & ctx.mask
    id_ = ctx.n_nodes
    ctx.node[id_].word = word
    ctx.node[id_].flags = flags
    ctx.node[id_].ordered = 0
    ctx.node[id_].first_child = ctx.n_child
    ctx.node[id_].n_children = n_kids
    for i in range(n_kids):
        ctx.child[ctx.n_child + i] = kids[i]
    ctx.n_child += n_kids
    ctx.n_nodes += 1
    ctx.bucket[slot] = id_
    return id_


cdef inline uint32_t cip_intern_dup(cip_ctx_t *ctx, uint32_t atom,
                                    uint32_t back) except? CIP_NO_ID:
    """A duplicate of `atom`, childless, `back` spheres below the atom it stands for -- see
    `cip_word_of` for why that count and not an absolute distance."""
    return cip_intern(ctx, cip_word_of(ctx.atoms[atom].element, ctx.atoms[atom].isotope, True, back),
                      0, NULL, 0)


cdef inline uint32_t cip_intern_hydrogen(cip_ctx_t *ctx) except? CIP_NO_ID:
    """An implicit hydrogen: a REAL childless atom, since the count in the arena stands for atoms the
    molecule has and not for phantoms."""
    return cip_intern(ctx, cip_word_of(1, 0, False, 0), 0, NULL, 0)


cdef inline void _cip_sort_ids_desc(cip_ctx_t *ctx, uint32_t *ids, uint32_t n) noexcept nogil:
    """Insertion sort of `ids` by node word descending, then by id ascending.  n <= CIP_KIDS_MAX.

    THIS IS THE INTERN ORDER, NOT THE CIP ORDER.  Its only job is to make the key canonical, so that two
    spellings of one child multiset hash the same; because an id is itself a function of the subtree,
    `(word, id)` is a function of the multiset and nothing else.  The CIP priority order is computed
    separately, into `ctx.corder` -- reordering `ctx.child` in place would invalidate every bucket
    already written.
    """
    cdef uint32_t i, j, v
    for i in range(1, n):
        v = ids[i]
        j = i
        while j > 0 and (ctx.node[ids[j - 1]].word < ctx.node[v].word or
                         (ctx.node[ids[j - 1]].word == ctx.node[v].word and ids[j - 1] > v)):
            ids[j] = ids[j - 1]
            j -= 1
        ids[j] = v


cdef uint32_t cip_subtree(cip_ctx_t *ctx, uint32_t atom, uint32_t parent, uint32_t in_order,
                          uint32_t depth) except? CIP_NO_ID:
    """The interned id of the branch rooted at `atom`, entered from `parent` across `in_order`.

    Equal id means identical truncated digraph.  `depth` is the sphere `atom` sits at, 1 for a direction
    of an anchor.  Five expansion clauses, and the node's children are all of them together:

    1. every incident bond except the one back to `parent`: one REAL child, plus `order - 1` duplicates
       of the far atom.  That is how rule 1a sees a double bond;
    2. the bond back to `parent`: `in_order - 1` duplicates OF THE PARENT and no real child.  A carbonyl
       oxygen entered from its carbon has a duplicate carbon among its children; omitting this is the
       classic way to rank an aldehyde below an alcohol;
    3. `at_implicit_h` childless real hydrogens.  They are a count in the arena and atoms nowhere, so
       nothing else would ever put them in the digraph -- and a CH3 / CH2R comparison needs them;
    4. a neighbour already on the walk's path within the same ring system: a childless duplicate, which
       is what terminates a ring.  THE ANCHOR IS ON THE PATH (see `cip_branch`), so a ring closing back
       onto it closes here rather than walking a second lap;
    5. CIP_NODE_UNDECIDABLE when a bond this node expands has order 4 or 8, when the implicit hydrogen
       count is H_UNKNOWN, or when any child carries the bit.  Order 4 needs the mancude mean phase 1
       does not compute, order 8 has no CIP precedence, and an unknown count is not a count.  The bit is
       part of the intern key, so a node carrying it never shares an id with one that does not.

    A NODE THE BIT REACHES IS CHILDLESS: the flag is decided from this node's own bonds and hydrogen
    count BEFORE any child is expanded, and a node carrying it returns at once.  That is not a shortcut
    with a cost -- the bit propagates to the root of the direction, `cip_rank_directions` refuses the
    whole site the moment a direction's root carries it, and no comparison whose answer is used can
    therefore contain one.  What it buys is the pathological case: without it, a stereocentre hanging off
    a large fused arene expands that arene's whole path-tree to learn what its first aromatic bond
    already said (2002 us on the worst record of `test/stereo.sdf`, 20 us with it).
    """
    if depth > CIP_DEPTH_MAX:
        raise OverflowError('CIP digraph deeper than %d spheres' % CIP_DEPTH_MAX)
    # the duplicate and hydrogen counts below are counted DOWN in a `while` rather than walked by a
    # `for ... in range`, because a range variable a loop body never reads is an unused entry, and this
    # tree treats a Cython warning as a build failure
    cdef uint32_t k, target, order, kid, back, reps
    cdef uint32_t kids[CIP_KIDS_MAX]
    cdef uint32_t n_kids = 0
    cdef uint8_t flags = 0
    cdef atom_t *a = &ctx.atoms[atom]
    if at_implicit_h_unknown(a):
        flags |= CIP_NODE_UNDECIDABLE
    if in_order == 4 or in_order == 8:
        flags |= CIP_NODE_UNDECIDABLE
    for k in range(ctx.ptr[atom], ctx.ptr[atom + 1]):
        if ctx.edges[k].order == 4 or ctx.edges[k].order == 8:
            flags |= CIP_NODE_UNDECIDABLE
            break
    if flags & CIP_NODE_UNDECIDABLE:
        return cip_intern(ctx, cip_rank_word(a, False, 0), flags, kids, 0)
    ctx.on_path[atom] = 1
    ctx.path_depth[atom] = depth
    try:
        for k in range(ctx.ptr[atom], ctx.ptr[atom + 1]):
            target = ctx.edges[k].to
            order = ctx.edges[k].order      # 1, 2 or 3 here: 4 and 8 returned above
            if ctx.on_path[target] and (target == parent or
                                        ctx.ring_system[target] == ctx.ring_system[atom]):
                # clauses 2 and 4.  Every duplicate here stands for an atom the path already holds, so
                # they all carry the same jump back up to it; the parent alone contributes no closure
                # duplicate, because it is this node's parent rather than a ring closure.
                back = depth + 1 - ctx.path_depth[target]
                reps = order        # 1, 2 or 3
                if target == parent:
                    reps -= 1
                while reps:
                    if n_kids >= CIP_KIDS_MAX:
                        raise OverflowError('CIP node with more than %d children' % CIP_KIDS_MAX)
                    kids[n_kids] = cip_intern_dup(ctx, target, back)
                    n_kids += 1
                    reps -= 1
                continue
            kid = cip_subtree(ctx, target, atom, order, depth + 1)                      # clause 1
            flags |= ctx.node[kid].flags & CIP_NODE_UNDECIDABLE                         # clause 5
            if n_kids >= CIP_KIDS_MAX:
                raise OverflowError('CIP node with more than %d children' % CIP_KIDS_MAX)
            kids[n_kids] = kid
            n_kids += 1
            reps = order - 1                                                            # clause 1
            while reps:
                if n_kids >= CIP_KIDS_MAX:
                    raise OverflowError('CIP node with more than %d children' % CIP_KIDS_MAX)
                kids[n_kids] = cip_intern_dup(ctx, target, 0)
                n_kids += 1
                reps -= 1
        if not at_implicit_h_unknown(a):
            reps = at_implicit_h(a)                                                     # clause 3
            while reps:
                if n_kids >= CIP_KIDS_MAX:
                    raise OverflowError('CIP node with more than %d children' % CIP_KIDS_MAX)
                kids[n_kids] = cip_intern_hydrogen(ctx)
                n_kids += 1
                reps -= 1
        _cip_sort_ids_desc(ctx, kids, n_kids)
        return cip_intern(ctx, cip_rank_word(a, False, 0), flags, kids, n_kids)
    finally:
        ctx.on_path[atom] = 0


cdef uint32_t cip_branch(cip_ctx_t *ctx, uint32_t atom, uint32_t parent,
                         uint32_t in_order) except? CIP_NO_ID:
    """One direction of a stereo unit, sphere 1.  The ANCHOR goes on the path first: a ring closing back
    onto it has to close on a duplicate, and left off the path the anchor would be expanded again and the
    branch would carry a whole extra lap of real atoms."""
    ctx.on_path[parent] = 1
    ctx.path_depth[parent] = 0
    try:
        return cip_subtree(ctx, atom, parent, in_order, 1)
    finally:
        ctx.on_path[parent] = 0


cdef inline int _cip_memo_get(cip_ctx_t *ctx, uint64_t key) noexcept nogil:
    """The recorded verdict for one node pair, or CIP_MEMO_ABSENT."""
    cdef uint32_t slot = <uint32_t> (((key * <uint64_t> 0x9e3779b97f4a7c15) >> 32) & ctx.memo_mask)
    while ctx.memo_key[slot] != CIP_MEMO_EMPTY:
        if ctx.memo_key[slot] == key:
            return ctx.memo_val[slot]
        slot = (slot + 1) & ctx.memo_mask
    return CIP_MEMO_ABSENT


cdef int _cip_memo_put(cip_ctx_t *ctx, uint64_t key, int verdict) except -1:
    """Record one verdict, growing the table past 70% load.  `n_memo` counts occupancy."""
    cdef uint32_t slot, i, cap
    cdef uint64_t *keys
    cdef int8_t *vals
    if (ctx.n_memo + 1) * 10 > (ctx.memo_mask + 1) * 7:
        cap = (ctx.memo_mask + 1) * 2
        keys = <uint64_t *> malloc(<size_t> cap * sizeof(uint64_t))
        vals = <int8_t *> malloc(<size_t> cap * sizeof(int8_t))
        if keys is NULL or vals is NULL:
            free(keys)
            free(vals)
            raise MemoryError()
        for i in range(cap):
            keys[i] = CIP_MEMO_EMPTY
        for i in range(ctx.memo_mask + 1):
            if ctx.memo_key[i] != CIP_MEMO_EMPTY:
                slot = <uint32_t> (((ctx.memo_key[i] * <uint64_t> 0x9e3779b97f4a7c15) >> 32)
                                   & (cap - 1))
                while keys[slot] != CIP_MEMO_EMPTY:
                    slot = (slot + 1) & (cap - 1)
                keys[slot] = ctx.memo_key[i]
                vals[slot] = ctx.memo_val[i]
        free(ctx.memo_key)
        free(ctx.memo_val)
        ctx.memo_key = keys
        ctx.memo_val = vals
        ctx.memo_mask = cap - 1
    slot = <uint32_t> (((key * <uint64_t> 0x9e3779b97f4a7c15) >> 32) & ctx.memo_mask)
    while ctx.memo_key[slot] != CIP_MEMO_EMPTY:
        if ctx.memo_key[slot] == key:
            return 0
        slot = (slot + 1) & ctx.memo_mask
    ctx.memo_key[slot] = key
    ctx.memo_val[slot] = <int8_t> verdict
    ctx.n_memo += 1
    return 0


cdef int _cip_frontier_grow(uint32_t **fa, uint32_t **fb, uint32_t **ga, uint32_t **gb,
                            uint32_t *cap, uint32_t want) except -1:
    """Grow ALL FOUR frontier buffers together, so one capacity describes every one of them.

    Growing only the next-sphere pair would be an overflow one swap later: after a swap that pair is the
    current frontier and the smaller pair is the one being filled.
    """
    cdef uint32_t c = cap[0]
    cdef uint32_t i
    cdef uint32_t *p
    cdef uint32_t **buffers[4]
    while c < want:
        c *= 2
    buffers[0] = fa
    buffers[1] = fb
    buffers[2] = ga
    buffers[3] = gb
    for i in range(4):
        p = <uint32_t *> realloc(buffers[i][0], <size_t> c * sizeof(uint32_t))
        if p is NULL:
            raise MemoryError()
        buffers[i][0] = p
    cap[0] = c
    return 0


cdef inline void _cip_frontier_swap(uint32_t **fa, uint32_t **fb,
                                    uint32_t **ga, uint32_t **gb) noexcept nogil:
    """The next sphere becomes the current one; the old current buffers are reused for the next."""
    cdef uint32_t *t = fa[0]
    fa[0] = ga[0]
    ga[0] = t
    t = fb[0]
    fb[0] = gb[0]
    gb[0] = t


cdef object _cip_node_key(cip_ctx_t *ctx, uint32_t id_):
    """The node unfolded into nested tuples: `(word, flags, children)`.  A TEST DOOR ONLY.

    An interned id is a counter in its own pool, so ids from two contexts say nothing about each other,
    while this unfolding is a function of the branch alone and compares across molecules.  It is
    exponential in the sharing it undoes, which is why nothing but a test calls it.
    """
    cdef uint32_t i
    cdef list kids = []
    for i in range(ctx.node[id_].n_children):
        kids.append(_cip_node_key(ctx, ctx.child[ctx.node[id_].first_child + i]))
    return (ctx.node[id_].word, ctx.node[id_].flags, tuple(kids))


# THE COMPARATOR.  Rule 1 is compared SPHERE BY SPHERE, not depth first, and the two genuinely differ: a
# difference at sphere 2 of a low-priority branch outranks a difference at sphere 3 of a high-priority
# one, so a lexicographic walk over the same digraph would answer a different question.  Three properties
# make the sphere-wise form affordable:
#
#   a == b -> 0 at once       interning is on the full key, so one id is one truncated digraph.  This is
#                             an exact tie certificate rather than a digest gamble
#   every verdict memoized    a steroid's centres re-ask the same pairs
#   a short child group       pads with word 0, which IS CIP's phantom of atomic number 0, so no padding
#                             node is materialised and phantoms cannot multiply down the spheres
#
# `cip_children_ordered` and `cip_compare_nodes` are mutually recursive and it terminates: ordering a
# node's children compares subtrees one sphere shallower than the node itself.
#
# NOTHING MAY INTERN WHILE A COMPARISON IS IN FLIGHT.  `cip_children_ordered` hands back a pointer into
# `ctx.corder`, and `_cip_reserve` can move that array; `cip_rank_directions` therefore interns all four
# directions before it sorts any of them.
cdef uint32_t *cip_children_ordered(cip_ctx_t *ctx, uint32_t id_) except NULL:
    """The node's children in descending CIP priority, cached on the node.

    The cache is `ctx.corder`, parallel to `ctx.child` and never a reordering of it: `ctx.child` is the
    intern key and every bucket already written depends on its order.  Ties fall back to the id, so the
    order is total and deterministic even where CIP is silent.
    """
    cdef uint32_t n = ctx.node[id_].n_children
    cdef uint32_t base = ctx.node[id_].first_child
    cdef uint32_t i, j, v
    cdef int cmp_
    if ctx.node[id_].ordered:
        return &ctx.corder[base]
    for i in range(n):
        ctx.corder[base + i] = ctx.child[base + i]
    for i in range(1, n):
        v = ctx.corder[base + i]
        j = i
        while j > 0:
            cmp_ = cip_compare_nodes(ctx, ctx.corder[base + j - 1], v)
            if cmp_ < 0 or (cmp_ == 0 and ctx.corder[base + j - 1] < v):
                break
            ctx.corder[base + j] = ctx.corder[base + j - 1]
            j -= 1
        ctx.corder[base + j] = v
    ctx.node[id_].ordered = 1
    return &ctx.corder[base]


cdef int cip_compare_nodes(cip_ctx_t *ctx, uint32_t a, uint32_t b) except -2:
    """-1 when subtree `a` outranks `b`, 1 when `b` outranks `a`, 0 when rules 1a/1b/2 tie.

    `fa[i]` and `fb[i]` are the paired nodes of the current sphere.  They are paired because the previous
    sphere compared EQUAL at every position -- the only way this loop continues -- so the two frontiers
    have equal length and their entries correspond one to one.
    """
    cdef uint64_t key
    cdef uint32_t *fa = NULL
    cdef uint32_t *fb = NULL
    cdef uint32_t *ga = NULL
    cdef uint32_t *gb = NULL
    cdef uint32_t *ka
    cdef uint32_t *kb
    cdef uint32_t cap = CIP_KIDS_MAX
    cdef uint32_t la = 1
    cdef uint32_t lb = 0
    cdef uint32_t i, j, ca, cb, kmax
    cdef uint64_t wa, wb
    cdef int verdict = 0
    cdef bint settled = False
    if a == b:
        return 0
    key = (<uint64_t> a << 32) | <uint64_t> b
    verdict = _cip_memo_get(ctx, key)
    if verdict != CIP_MEMO_ABSENT:
        return verdict
    verdict = 0
    if ctx.node[a].word != ctx.node[b].word:
        verdict = -1 if ctx.node[a].word > ctx.node[b].word else 1
        _cip_memo_put(ctx, key, verdict)
        return verdict
    fa = <uint32_t *> malloc(<size_t> cap * sizeof(uint32_t))
    fb = <uint32_t *> malloc(<size_t> cap * sizeof(uint32_t))
    ga = <uint32_t *> malloc(<size_t> cap * sizeof(uint32_t))
    gb = <uint32_t *> malloc(<size_t> cap * sizeof(uint32_t))
    if fa is NULL or fb is NULL or ga is NULL or gb is NULL:
        free(fa)
        free(fb)
        free(ga)
        free(gb)
        raise MemoryError()
    try:
        fa[0] = a
        fb[0] = b
        while la:
            lb = 0
            for i in range(la):
                ka = cip_children_ordered(ctx, fa[i])
                kb = cip_children_ordered(ctx, fb[i])
                ca = ctx.node[fa[i]].n_children
                cb = ctx.node[fb[i]].n_children
                kmax = ca if ca > cb else cb
                for j in range(kmax):
                    # a child group shorter than its partner pads with word 0: CIP's phantom
                    wa = ctx.node[ka[j]].word if j < ca else 0
                    wb = ctx.node[kb[j]].word if j < cb else 0
                    if wa != wb:
                        verdict = -1 if wa > wb else 1
                        settled = True
                        break
                if settled:
                    break
                # the groups tied, so ca == cb and the next sphere stays paired
                if lb + ca > cap:
                    _cip_frontier_grow(&fa, &fb, &ga, &gb, &cap, lb + ca)
                    ka = cip_children_ordered(ctx, fa[i])
                    kb = cip_children_ordered(ctx, fb[i])
                for j in range(ca):
                    ga[lb + j] = ka[j]
                    gb[lb + j] = kb[j]
                lb += ca
            if settled:
                break
            _cip_frontier_swap(&fa, &fb, &ga, &gb)
            la = lb
        _cip_memo_put(ctx, key, verdict)
        return verdict
    finally:
        free(fa)
        free(fb)
        free(ga)
        free(gb)


cdef inline int _cip_cmp_direction(cip_ctx_t *ctx, uint32_t x, uint32_t y) except -2:
    """As `cip_compare_nodes`, with CIP_NO_ID -- a direction with no atom of its own -- lowest.

    An SU_NO_REF direction is the site's implicit hydrogen or its lone pair.  Both rank below every named
    ligand: the lone pair is CIP's phantom, and hydrogen is the lightest element there is, so one branch
    serves both.
    """
    if x == y:
        return 0
    if x == CIP_NO_ID:
        return 1
    if y == CIP_NO_ID:
        return -1
    return cip_compare_nodes(ctx, x, y)


cdef int cip_rank_directions(cip_ctx_t *ctx, stereo_unit_t *u, uint32_t *order_out) except -1:
    """`order_out[0:4]`: the unit's refs POSITIONS in descending CIP priority.

    Returns 0 on a full ranking and 1 when rules 1a/1b/2 do not settle the unit -- a tie between two
    directions, or a ligand carrying CIP_NODE_UNDECIDABLE.  The caller refuses on 1; it never breaks a
    tie itself.
    """
    cdef uint32_t ids[4]
    cdef uint32_t i, j, v
    cdef int cmp_
    cdef halfedge_t *e
    # EVERY SUBTREE FIRST: interning may move ctx.corder, and the sort below holds pointers into it
    for i in range(4):
        order_out[i] = i
        ids[i] = CIP_NO_ID
        if u.refs[i] != SU_NO_REF:
            e = csr_find_at(ctx.ptr, ctx.edges, u.anchor, u.refs[i])
            if e is NULL:
                raise KeyError('a stereo unit names a direction that is not a bond')
            ids[i] = cip_branch(ctx, u.refs[i], u.anchor, e.order)
            if ctx.node[ids[i]].flags & CIP_NODE_UNDECIDABLE:
                return 1
    for i in range(1, 4):
        v = order_out[i]
        j = i
        while j > 0:
            cmp_ = _cip_cmp_direction(ctx, ids[order_out[j - 1]], ids[v])
            if cmp_ <= 0:
                break
            order_out[j] = order_out[j - 1]
            j -= 1
        order_out[j] = v
    for i in range(3):
        if _cip_cmp_direction(ctx, ids[order_out[i]], ids[order_out[i + 1]]) == 0:
            return 1
    return 0


# --- the letter a ranking implies ------------------------------------------------------------------
# The codes are ATOM_CIP_CODES' indices, which _molecule_container.pxi defines as
# (None, 'R', 'S', 'r', 's', 'M', 'P', 'm', 'p').  Named here rather than spelled as 1 and 2 at the one
# call site, because the pair that matters is (code, parity) and the derivation below ties them.
DEF CIP_CODE_R = 1
DEF CIP_CODE_S = 2

# WHICH LETTER AN EVEN PARITY MEANS, in `translate_parity`'s three-state encoding (1 even, 2 odd) and in
# DESCENDING PRIORITY order.  Derived, not chosen, and derived twice:
#
#   * the SMILES reader takes `@` to be parity 2, odd, in the SMILES neighbour order; OpenSMILES says `@`
#     means the last three neighbours run anticlockwise viewed FROM the first, while CIP views a centre
#     with its LOWEST-priority ligand pointing away and calls 1->2->3 clockwise R -- the opposite
#     viewpoint, which reverses the sense.  So for a neighbour order (lowest, r1, r2, r3) `@` is R, and
#     re-expressing (lowest, r1, r2, r3) as (r1, r2, r3, lowest) is a 4-cycle, which is odd: odd becomes
#     even, so R is even;
#   * `[C@H](F)(Cl)Br` is (S)-bromochlorofluoromethane, and
#     `test_assign_labels_bromochlorofluoromethane_as_written` is that calibration.
#
# If the calibration fails, re-derive.  Do not flip this constant without writing down the derivation
# that flips with it.
DEF CIP_R_PARITY = 1
