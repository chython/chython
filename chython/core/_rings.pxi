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
# Ring perception: bridge detection and Vismara relevant-cycle prototypes.
#
# Everything that touches the graph runs in C on malloc'd arrays with fixed-width bitsets. The
# only Python-level code here is the exception raising in `perceive_rings`, which happens once
# per structure and never inside a loop over atoms or bonds.
#
# The predecessor of this file did the prototype half with Python ints, lists, dicts and
# frozensets. On a 650-atom / 937-bond honeycomb flake that spent 73% of its time in the object
# machinery alone -- bignum arithmetic 19%, int boxing 18%, list operations 15%, allocator and
# gc 11%, generic compare/hash 10% -- and allocated about 2.9 million Python objects per
# derivation, 422500 of them the empty predecessor lists of one BFS per root.
#


# Guards against the pathological end of the input space. Both raise rather than degrade:
# a partial ring segment would make bond_in_ring, ring_sizes_of and the isomorphism feature
# words quietly disagree with the graph, and every consumer of those trusts them.
cdef Py_ssize_t PROTOTYPE_LIMIT = 200000
cdef time_t DEADLINE_SECONDS = 60

cdef Py_ssize_t HASH_INIT = 1024
cdef Py_ssize_t POOL_INIT = 64


ctypedef struct ring_ctx_t:
    uint32_t n
    Py_ssize_t ne                 # ring bonds, i.e. edges surviving bridge detection
    Py_ssize_t wv                 # uint64 words in a vertex bitset
    Py_ssize_t we                 # uint64 words in an edge bitset
    time_t deadline

    # CSR over the ring subgraph only. reid[k] is the canonical id of the edge in slot k, the
    # same id in both of its two slots, so an edge bitset needs no lookup table.
    uint32_t *rptr                # n + 1
    uint32_t *radj                # 2 * ne
    uint32_t *reid                # 2 * ne
    uint32_t *eu                  # ne
    uint32_t *ev                  # ne

    # per-root BFS scratch, reused for every root
    int32_t *dv                   # n, -1 = unreached
    uint32_t *queue               # n
    uint32_t *pred_cnt            # n
    uint32_t *pred_val            # 2 * ne, vertex v's slots start at rptr[v]
    uint64_t *vset                # n * wv
    uint32_t *cycbuf              # n + 2
    uint32_t *cycbuf2             # n + 2
    uint64_t *tmprow              # we

    # prototype pool
    Py_ssize_t pcount
    Py_ssize_t pcap
    uint32_t *psize               # pcap, cycle length
    uint64_t *pvbits              # pcap * wv, vertices on some cycle of the prototype
    uint64_t *pebits              # pcap * we, edges of the representative
    Py_ssize_t *pcyc_ofs          # pcap, offset into pcyc
    uint32_t *pcyc                # cyc_cap, representative vertex sequences
    Py_ssize_t cyc_len
    Py_ssize_t cyc_cap

    # open-addressed dedupe over pebits rows
    Py_ssize_t hcap
    Py_ssize_t *hslot             # hcap, -1 = empty, else a prototype index

    # GF(2) elimination table indexed by pivot bit
    uint64_t *basis               # ne * we
    uint8_t *basis_used           # ne

    Py_ssize_t *order             # pcap, prototype indices in (size, cycle) order
    Py_ssize_t *order_tmp         # pcap, merge sort scratch
    Py_ssize_t *keep              # pcap, the relevant prototypes
    Py_ssize_t keep_count
    Py_ssize_t *basis_idx         # pcap, the minimum cycle basis
    Py_ssize_t basis_count


cdef inline int _hi_bit64(uint64_t w) noexcept nogil:
    # index of the highest set bit; w must be nonzero. Hand-rolled rather than
    # __builtin_clzll for the same reason _popcount64 is: MSVC has neither.
    cdef int r = 0
    if w >> 32:
        w >>= 32
        r += 32
    if w >> 16:
        w >>= 16
        r += 16
    if w >> 8:
        w >>= 8
        r += 8
    if w >> 4:
        w >>= 4
        r += 4
    if w >> 2:
        w >>= 2
        r += 2
    if w >> 1:
        r += 1
    return r


cdef inline int _lo_bit64(uint64_t w) noexcept nogil:
    return _hi_bit64(w & (~w + <uint64_t> 1))


cdef inline uint32_t _popcount64(uint64_t w) noexcept nogil:
    cdef uint32_t c = 0
    while w:
        w &= w - 1
        c += 1
    return c


cdef inline Py_ssize_t _row_hi(uint64_t *row, Py_ssize_t we) noexcept nogil:
    cdef Py_ssize_t w = we - 1
    while w >= 0:
        if row[w]:
            return w * 64 + _hi_bit64(row[w])
        w -= 1
    return -1


cdef void _free_ctx(ring_ctx_t *ctx) noexcept nogil:
    free(ctx.rptr); free(ctx.radj); free(ctx.reid); free(ctx.eu); free(ctx.ev)
    free(ctx.dv); free(ctx.queue); free(ctx.pred_cnt); free(ctx.pred_val)
    free(ctx.vset); free(ctx.cycbuf); free(ctx.cycbuf2); free(ctx.tmprow)
    free(ctx.psize); free(ctx.pvbits); free(ctx.pebits); free(ctx.pcyc_ofs); free(ctx.pcyc)
    free(ctx.hslot); free(ctx.basis); free(ctx.basis_used)
    free(ctx.order); free(ctx.order_tmp); free(ctx.keep); free(ctx.basis_idx)
    memset(ctx, 0, sizeof(ring_ctx_t))


cdef int _prepare_ctx(ring_ctx_t *ctx, Structure structure) noexcept nogil:
    """Build the ring-subgraph CSR and every fixed-size scratch buffer.

    Leaves ctx.ne == 0 when bridge detection found no ring bond, which the caller reads as
    "nothing to perceive" -- the descriptors have already been cleared by then.
    """
    cdef uint32_t n = structure.header.atom_count
    cdef uint32_t *ptr = csr_ptr(structure)
    cdef halfedge_t *edges = csr_edges(structure)
    cdef uint32_t v, u, k, t, nid, half
    cdef Py_ssize_t i, ne, wv, we

    ctx.n = n
    ctx.rptr = <uint32_t *> calloc(<size_t> n + 1, sizeof(uint32_t))
    ctx.dv = <int32_t *> malloc(<size_t> n * sizeof(int32_t))
    ctx.queue = <uint32_t *> malloc(<size_t> n * sizeof(uint32_t))
    ctx.pred_cnt = <uint32_t *> malloc(<size_t> n * sizeof(uint32_t))
    ctx.cycbuf = <uint32_t *> malloc((<size_t> n + 2) * sizeof(uint32_t))
    ctx.cycbuf2 = <uint32_t *> malloc((<size_t> n + 2) * sizeof(uint32_t))
    if (ctx.rptr is NULL or ctx.dv is NULL or ctx.queue is NULL or ctx.pred_cnt is NULL
            or ctx.cycbuf is NULL or ctx.cycbuf2 is NULL):
        return -1

    # pass 1: ring degrees into rptr[v + 1], then prefix sums in place
    for v in range(n):
        for k in range(ptr[v], ptr[v + 1]):
            if edges[k].flags & HE_IN_RING:
                ctx.rptr[v + 1] += 1
    for v in range(n):
        ctx.rptr[v + 1] += ctx.rptr[v]
    half = ctx.rptr[n]
    if half == 0:
        ctx.ne = 0
        return 0
    ne = <Py_ssize_t> (half // 2)
    ctx.ne = ne
    wv = (<Py_ssize_t> n + 63) // 64
    we = (ne + 63) // 64
    ctx.wv = wv
    ctx.we = we

    ctx.radj = <uint32_t *> malloc(<size_t> half * sizeof(uint32_t))
    ctx.reid = <uint32_t *> malloc(<size_t> half * sizeof(uint32_t))
    ctx.eu = <uint32_t *> malloc(<size_t> ne * sizeof(uint32_t))
    ctx.ev = <uint32_t *> malloc(<size_t> ne * sizeof(uint32_t))
    ctx.pred_val = <uint32_t *> malloc(<size_t> half * sizeof(uint32_t))
    ctx.vset = <uint64_t *> malloc(<size_t> n * <size_t> wv * sizeof(uint64_t))
    ctx.tmprow = <uint64_t *> malloc(<size_t> we * sizeof(uint64_t))
    ctx.basis = <uint64_t *> calloc(<size_t> ne * <size_t> we, sizeof(uint64_t))
    ctx.basis_used = <uint8_t *> calloc(<size_t> ne, sizeof(uint8_t))
    ctx.hslot = <Py_ssize_t *> malloc(<size_t> HASH_INIT * sizeof(Py_ssize_t))
    if (ctx.radj is NULL or ctx.reid is NULL or ctx.eu is NULL or ctx.ev is NULL
            or ctx.pred_val is NULL or ctx.vset is NULL or ctx.tmprow is NULL
            or ctx.basis is NULL or ctx.basis_used is NULL or ctx.hslot is NULL):
        return -1
    ctx.hcap = HASH_INIT
    for i in range(HASH_INIT):
        ctx.hslot[i] = -1

    # pass 2: neighbours, in the original CSR order so edge ids come out canonical
    for v in range(n):
        ctx.pred_cnt[v] = 0
    for v in range(n):
        for k in range(ptr[v], ptr[v + 1]):
            if edges[k].flags & HE_IN_RING:
                ctx.radj[ctx.rptr[v] + ctx.pred_cnt[v]] = edges[k].to
                ctx.pred_cnt[v] += 1

    # pass 3: dense edge ids ascending, backfilling the twin slot
    nid = 0
    for v in range(n):
        for k in range(ctx.rptr[v], ctx.rptr[v + 1]):
            u = ctx.radj[k]
            if u > v:
                ctx.reid[k] = nid
                ctx.eu[nid] = v
                ctx.ev[nid] = u
                for t in range(ctx.rptr[u], ctx.rptr[u + 1]):
                    if ctx.radj[t] == v:
                        ctx.reid[t] = nid
                        break
                nid += 1
    if <Py_ssize_t> nid != ne:
        return -1                 # a half-edge without its twin: the CSR is not symmetric
    return 0


cdef int _grow_pool(ring_ctx_t *ctx) noexcept nogil:
    cdef Py_ssize_t cap = ctx.pcap * 2 if ctx.pcap else POOL_INIT
    cdef void *p
    p = realloc(ctx.psize, <size_t> cap * sizeof(uint32_t))
    if p is NULL:
        return -1
    ctx.psize = <uint32_t *> p
    p = realloc(ctx.pvbits, <size_t> cap * <size_t> ctx.wv * sizeof(uint64_t))
    if p is NULL:
        return -1
    ctx.pvbits = <uint64_t *> p
    p = realloc(ctx.pebits, <size_t> cap * <size_t> ctx.we * sizeof(uint64_t))
    if p is NULL:
        return -1
    ctx.pebits = <uint64_t *> p
    p = realloc(ctx.pcyc_ofs, <size_t> cap * sizeof(Py_ssize_t))
    if p is NULL:
        return -1
    ctx.pcyc_ofs = <Py_ssize_t *> p
    p = realloc(ctx.order, <size_t> cap * sizeof(Py_ssize_t))
    if p is NULL:
        return -1
    ctx.order = <Py_ssize_t *> p
    p = realloc(ctx.order_tmp, <size_t> cap * sizeof(Py_ssize_t))
    if p is NULL:
        return -1
    ctx.order_tmp = <Py_ssize_t *> p
    p = realloc(ctx.keep, <size_t> cap * sizeof(Py_ssize_t))
    if p is NULL:
        return -1
    ctx.keep = <Py_ssize_t *> p
    p = realloc(ctx.basis_idx, <size_t> cap * sizeof(Py_ssize_t))
    if p is NULL:
        return -1
    ctx.basis_idx = <Py_ssize_t *> p
    ctx.pcap = cap
    return 0


cdef int _grow_cyc(ring_ctx_t *ctx, Py_ssize_t need) noexcept nogil:
    cdef Py_ssize_t cap = ctx.cyc_cap
    cdef void *p
    if ctx.cyc_len + need <= cap:
        return 0
    if cap == 0:
        cap = 1024
    while cap < ctx.cyc_len + need:
        cap *= 2
    p = realloc(ctx.pcyc, <size_t> cap * sizeof(uint32_t))
    if p is NULL:
        return -1
    ctx.pcyc = <uint32_t *> p
    ctx.cyc_cap = cap
    return 0


cdef inline uint64_t _row_hash(uint64_t *row, Py_ssize_t we) noexcept nogil:
    cdef uint64_t h = <uint64_t> 14695981039346656037UL
    cdef Py_ssize_t w
    for w in range(we):
        h ^= row[w]
        h *= <uint64_t> 1099511628211UL
    h ^= h >> 33
    h *= <uint64_t> 0xff51afd7ed558ccd
    h ^= h >> 33
    return h


cdef Py_ssize_t _seen_lookup(ring_ctx_t *ctx, uint64_t *row) noexcept nogil:
    """Index of the prototype whose edge set is `row`, or -1.

    Replaces a frozenset of edge ids per candidate plus a set of those frozensets. Two
    prototypes with the same edge set are the same cycle, so the row is the whole key.
    """
    cdef Py_ssize_t mask = ctx.hcap - 1
    cdef Py_ssize_t pos = <Py_ssize_t> (_row_hash(row, ctx.we) & <uint64_t> mask)
    cdef Py_ssize_t idx
    cdef size_t nbytes = <size_t> ctx.we * sizeof(uint64_t)
    while True:
        idx = ctx.hslot[pos]
        if idx < 0:
            return -1
        if memcmp(ctx.pebits + <size_t> idx * <size_t> ctx.we, row, nbytes) == 0:
            return idx
        pos = (pos + 1) & mask


cdef inline void _hash_put(ring_ctx_t *ctx, Py_ssize_t idx) noexcept nogil:
    cdef Py_ssize_t mask = ctx.hcap - 1
    cdef Py_ssize_t pos = <Py_ssize_t> (
        _row_hash(ctx.pebits + <size_t> idx * <size_t> ctx.we, ctx.we) & <uint64_t> mask)
    while ctx.hslot[pos] >= 0:
        pos = (pos + 1) & mask
    ctx.hslot[pos] = idx


cdef int _seen_place(ring_ctx_t *ctx, Py_ssize_t idx) noexcept nogil:
    cdef Py_ssize_t cap, i
    cdef Py_ssize_t *slots
    if (idx + 1) * 2 >= ctx.hcap:
        cap = ctx.hcap * 2
        slots = <Py_ssize_t *> malloc(<size_t> cap * sizeof(Py_ssize_t))
        if slots is NULL:
            return -1
        free(ctx.hslot)
        ctx.hslot = slots
        ctx.hcap = cap
        for i in range(cap):
            ctx.hslot[i] = -1
        for i in range(idx):
            _hash_put(ctx, i)
    _hash_put(ctx, idx)
    return 0


cdef inline int _cmp_proto(ring_ctx_t *ctx, Py_ssize_t a, Py_ssize_t b) noexcept nogil:
    # Size first, then the edge bitset by memcmp. Distinct prototypes have distinct edge sets,
    # so this is a total order -- which is all the size-class grouping in _filter_relevant
    # needs. It deliberately says nothing about vertex order inside a cycle: a canonical
    # rotation costs a rotate and a reversed compare per prototype and no caller reads it.
    cdef size_t nbytes = <size_t> ctx.we * sizeof(uint64_t)
    if ctx.psize[a] != ctx.psize[b]:
        return -1 if ctx.psize[a] < ctx.psize[b] else 1
    return memcmp(ctx.pebits + <size_t> a * <size_t> ctx.we,
                  ctx.pebits + <size_t> b * <size_t> ctx.we, nbytes)


cdef void _sort_order(ring_ctx_t *ctx) noexcept nogil:
    """Bottom-up merge sort of prototype indices by (size, edge set)."""
    cdef Py_ssize_t nn = ctx.pcount
    cdef Py_ssize_t width = 1
    cdef Py_ssize_t i, l, m, r, a, b, o
    cdef Py_ssize_t *src = ctx.order
    cdef Py_ssize_t *dst = ctx.order_tmp
    cdef Py_ssize_t *swap
    for i in range(nn):
        ctx.order[i] = i
    while width < nn:
        i = 0
        while i < nn:
            l = i
            m = i + width
            r = i + 2 * width
            if m > nn:
                m = nn
            if r > nn:
                r = nn
            a = l
            b = m
            o = l
            while a < m and b < r:
                if _cmp_proto(ctx, src[b], src[a]) < 0:
                    dst[o] = src[b]
                    b += 1
                else:
                    dst[o] = src[a]
                    a += 1
                o += 1
            while a < m:
                dst[o] = src[a]
                a += 1
                o += 1
            while b < r:
                dst[o] = src[b]
                b += 1
                o += 1
            i = r
        swap = src
        src = dst
        dst = swap
        width *= 2
    if src is not ctx.order:
        memcpy(ctx.order, src, <size_t> nn * sizeof(Py_ssize_t))


cdef inline Py_ssize_t _edge_id(ring_ctx_t *ctx, uint32_t a, uint32_t b) noexcept nogil:
    cdef uint32_t k
    for k in range(ctx.rptr[a], ctx.rptr[a + 1]):
        if ctx.radj[k] == b:
            return <Py_ssize_t> ctx.reid[k]
    return -1


cdef inline bint _meets_only_root(ring_ctx_t *ctx, uint32_t y, uint32_t z,
                                  uint32_t root) noexcept nogil:
    """vset[y] & vset[z] == {root}, the Vismara prototype admission test."""
    cdef Py_ssize_t w
    cdef Py_ssize_t rw = <Py_ssize_t> (root >> 6)
    cdef uint64_t got, want
    cdef uint64_t *ry = ctx.vset + <size_t> y * <size_t> ctx.wv
    cdef uint64_t *rz = ctx.vset + <size_t> z * <size_t> ctx.wv
    for w in range(ctx.wv):
        got = ry[w] & rz[w]
        want = (<uint64_t> 1 << (root & 63)) if w == rw else <uint64_t> 0
        if got != want:
            return False
    return True


cdef Py_ssize_t _representative(ring_ctx_t *ctx, uint32_t root, uint32_t y, uint32_t z,
                                int64_t apex) noexcept nogil:
    """One cycle of the prototype (root, y, z, apex) into cycbuf; returns its length.

    Following the first predecessor of each arm is enough. The admission test makes every
    shortest root..y path vertex-disjoint from every shortest root..z path apart from the root,
    and the apex, when present, is strictly farther from the root than either arm, so it cannot
    lie on one. The walk is therefore always a simple cycle, and every cycle the prototype
    generates has the same length and the same relevance -- which is why one representative
    decides the prototype and expanding over all shortest-path combinations buys nothing -- and
    costs without bound: a 150-atom macrocyclic aryl sulfone with 20 para-substituted benzenes has
    one 100-membered relevant cycle per arc choice, 2**20 of them, measured at 36 s and a 423 MB
    arena where the prototype is 21 records.
    """
    cdef Py_ssize_t la = 0
    cdef Py_ssize_t lb = 0
    cdef Py_ssize_t i
    cdef uint32_t cur = y
    ctx.cycbuf[la] = y
    la += 1
    while cur != root:
        cur = ctx.pred_val[ctx.rptr[cur]]
        ctx.cycbuf[la] = cur
        la += 1
    cur = z
    ctx.cycbuf2[lb] = z
    lb += 1
    while cur != root:
        cur = ctx.pred_val[ctx.rptr[cur]]
        ctx.cycbuf2[lb] = cur
        lb += 1
    # cycbuf runs y..root and cycbuf2 runs z..root; walk y -> root -> z, then close via apex
    for i in range(lb - 2, -1, -1):
        ctx.cycbuf[la] = ctx.cycbuf2[i]
        la += 1
    if apex >= 0:
        ctx.cycbuf[la] = <uint32_t> apex
        la += 1
    return la


cdef int _emit(ring_ctx_t *ctx, uint32_t root, uint32_t y, uint32_t z,
               int64_t apex) noexcept nogil:
    """Record the prototype unless its representative is degenerate or already seen."""
    cdef Py_ssize_t L = _representative(ctx, root, y, z, apex)
    cdef Py_ssize_t we = ctx.we
    cdef Py_ssize_t wv = ctx.wv
    cdef Py_ssize_t i, j, eid, idx
    cdef uint64_t *row = ctx.tmprow
    cdef uint64_t *pv
    memset(row, 0, <size_t> we * sizeof(uint64_t))
    for i in range(L):
        # wrap by branch, not by %: signed modulo drags in Cython's sign-correcting helper,
        # and this runs once per edge of every candidate cycle
        j = i + 1
        if j == L:
            j = 0
        eid = _edge_id(ctx, ctx.cycbuf[i], ctx.cycbuf[j])
        if eid < 0:
            return 0                                  # closing edge is not a ring bond
        if row[eid >> 6] >> (eid & 63) & 1:
            return 0                                  # an edge twice: not a simple cycle
        row[eid >> 6] |= <uint64_t> 1 << (eid & 63)
    if _seen_lookup(ctx, row) >= 0:
        return 0
    if ctx.pcount == ctx.pcap and _grow_pool(ctx):
        return -1
    if _grow_cyc(ctx, L):
        return -1
    idx = ctx.pcount
    ctx.psize[idx] = <uint32_t> L
    memcpy(ctx.pebits + <size_t> idx * <size_t> we, row, <size_t> we * sizeof(uint64_t))
    pv = ctx.pvbits + <size_t> idx * <size_t> wv
    for i in range(wv):
        pv[i] = (ctx.vset[<size_t> y * <size_t> wv + i]
                 | ctx.vset[<size_t> z * <size_t> wv + i])
    if apex >= 0:
        pv[apex >> 6] |= <uint64_t> 1 << (apex & 63)
    ctx.pcyc_ofs[idx] = ctx.cyc_len
    memcpy(ctx.pcyc + ctx.cyc_len, ctx.cycbuf, <size_t> L * sizeof(uint32_t))
    ctx.cyc_len += L
    if _seen_place(ctx, idx):
        return -1
    ctx.pcount = idx + 1
    return 0


cdef int _build_prototypes(ring_ctx_t *ctx) noexcept nogil:
    """Every relevant-cycle prototype candidate, deduplicated.

    A prototype is Vismara's (root, y, z, apex): a root, two vertices equidistant from it, and
    either an apex above both (even cycles) or the ring bond y-z (odd cycles). It stands for all
    the cycles obtained by choosing a shortest root..y path and a shortest root..z path -- all
    of one length, all relevant or none. There are at most mu*|E| prototypes, and that
    polynomial bound is on the prototypes, never on the cycles they generate.

    Candidates from a single BFS *tree* per root miss relevant cycles whenever shortest paths
    tie (K3,3 yields 8 of 9), so the search records every predecessor and works over the
    shortest-path DAG. Restricting each search to {w : w >= root} finds every cycle exactly
    once, from its lowest-numbered vertex, and collapses expansion work to the prototype count
    -- without it a 98-atom graphene flake generates 13567 candidates instead of 36.
    """
    cdef uint32_t n = ctx.n
    cdef Py_ssize_t wv = ctx.wv
    cdef uint32_t root, v, cur, nb, y, z, k, p
    cdef Py_ssize_t head, tail, si, pi, pj, w, i
    cdef int32_t d
    cdef uint64_t *acc
    cdef int rc

    for root in range(n):
        if ctx.rptr[root] == ctx.rptr[root + 1]:
            continue                                  # not on any ring bond
        if time(NULL) > ctx.deadline:
            return -3

        memset(ctx.dv, 0xff, <size_t> n * sizeof(int32_t))
        memset(ctx.pred_cnt, 0, <size_t> n * sizeof(uint32_t))
        ctx.dv[root] = 0
        ctx.queue[0] = root
        head = 0
        tail = 1
        while head < tail:
            cur = ctx.queue[head]
            head += 1
            d = ctx.dv[cur] + 1
            for k in range(ctx.rptr[cur], ctx.rptr[cur + 1]):
                nb = ctx.radj[k]
                if nb < root:
                    continue
                if ctx.dv[nb] < 0:
                    ctx.dv[nb] = d
                    ctx.pred_val[ctx.rptr[nb] + ctx.pred_cnt[nb]] = cur
                    ctx.pred_cnt[nb] += 1
                    ctx.queue[tail] = nb
                    tail += 1
                elif ctx.dv[nb] == d:
                    ctx.pred_val[ctx.rptr[nb] + ctx.pred_cnt[nb]] = cur
                    ctx.pred_cnt[nb] += 1

        # vset[x] = every vertex on some shortest root..x path. The queue is in nondecreasing
        # distance order, so one forward pass suffices. Rows of unreached vertices stay stale
        # from the previous root and are never read: every use is guarded by dv >= 1.
        for si in range(tail):
            cur = ctx.queue[si]
            acc = ctx.vset + <size_t> cur * <size_t> wv
            memset(acc, 0, <size_t> wv * sizeof(uint64_t))
            acc[cur >> 6] |= <uint64_t> 1 << (cur & 63)
            for pi in range(ctx.pred_cnt[cur]):
                p = ctx.pred_val[ctx.rptr[cur] + pi]
                for w in range(wv):
                    acc[w] |= ctx.vset[<size_t> p * <size_t> wv + w]

        # odd prototypes: a ring bond whose endpoints are equidistant from the root
        for i in range(ctx.ne):
            y = ctx.eu[i]
            z = ctx.ev[i]
            if ctx.dv[y] < 1 or ctx.dv[z] < 1 or ctx.dv[y] != ctx.dv[z]:
                continue
            if not _meets_only_root(ctx, y, z, root):
                continue
            rc = _emit(ctx, root, y, z, -1)
            if rc:
                return rc

        # even prototypes: an apex with two distinct DAG predecessors
        for v in range(n):
            if ctx.dv[v] < 2:
                continue
            for pi in range(ctx.pred_cnt[v]):
                y = ctx.pred_val[ctx.rptr[v] + pi]
                for pj in range(pi + 1, ctx.pred_cnt[v]):
                    z = ctx.pred_val[ctx.rptr[v] + pj]
                    if not _meets_only_root(ctx, y, z, root):
                        continue
                    rc = _emit(ctx, root, y, z, <int64_t> v)
                    if rc:
                        return rc

        if ctx.pcount > PROTOTYPE_LIMIT:
            return -2
    return 0


cdef Py_ssize_t _reduce_row(ring_ctx_t *ctx, uint64_t *row) noexcept nogil:
    """Reduce `row` against the elimination table in place; pivot bit, or -1 if it vanished."""
    cdef Py_ssize_t we = ctx.we
    cdef Py_ssize_t h = _row_hi(row, we)
    cdef Py_ssize_t w
    cdef uint64_t *br
    while h >= 0:
        if not ctx.basis_used[h]:
            return h
        br = ctx.basis + <size_t> h * <size_t> we
        for w in range(we):
            row[w] ^= br[w]
        h = _row_hi(row, we)
    return -1


cdef inline void _basis_store(ring_ctx_t *ctx, Py_ssize_t pivot, uint64_t *row) noexcept nogil:
    memcpy(ctx.basis + <size_t> pivot * <size_t> ctx.we, row,
           <size_t> ctx.we * sizeof(uint64_t))
    ctx.basis_used[pivot] = 1


cdef int _filter_relevant(ring_ctx_t *ctx) noexcept nogil:
    """Keep prototypes whose representative is not spanned by strictly shorter ones.

    One representative per prototype is sufficient and is what makes the filter polynomial:
    relevance is uniform over a prototype, so the representative's verdict is the family's.
    """
    cdef Py_ssize_t total = ctx.pcount
    cdef Py_ssize_t we = ctx.we
    cdef Py_ssize_t i = 0
    cdef Py_ssize_t j, k, idx, pivot
    cdef uint32_t size
    _sort_order(ctx)
    memset(ctx.basis_used, 0, <size_t> ctx.ne * sizeof(uint8_t))
    ctx.keep_count = 0
    while i < total:
        if time(NULL) > ctx.deadline:
            return -3
        size = ctx.psize[ctx.order[i]]
        j = i
        while j < total and ctx.psize[ctx.order[j]] == size:
            j += 1
        # classify this size class against the basis of strictly smaller cycles
        for k in range(i, j):
            idx = ctx.order[k]
            memcpy(ctx.tmprow, ctx.pebits + <size_t> idx * <size_t> we,
                   <size_t> we * sizeof(uint64_t))
            if _reduce_row(ctx, ctx.tmprow) >= 0:
                ctx.keep[ctx.keep_count] = idx
                ctx.keep_count += 1
        # then admit the whole class
        for k in range(i, j):
            idx = ctx.order[k]
            memcpy(ctx.tmprow, ctx.pebits + <size_t> idx * <size_t> we,
                   <size_t> we * sizeof(uint64_t))
            pivot = _reduce_row(ctx, ctx.tmprow)
            if pivot >= 0:
                _basis_store(ctx, pivot, ctx.tmprow)
        i = j
    return 0


cdef int _select_basis(ring_ctx_t *ctx) noexcept nogil:
    """A minimum cycle basis, one representative each, greedily in nondecreasing size.

    Greedy insertion over relevant-cycle representatives reaches full rank, so the result has
    exactly `mu = |E| - |V| + components` members without mu ever being computed: the rank of
    the relevant cycles is mu by definition, and no representative that raises the rank is
    skipped. Taking the shortest independent ones first is what makes the basis minimum.
    """
    cdef Py_ssize_t we = ctx.we
    cdef Py_ssize_t i, idx, pivot
    memset(ctx.basis_used, 0, <size_t> ctx.ne * sizeof(uint8_t))
    ctx.basis_count = 0
    for i in range(ctx.keep_count):
        idx = ctx.keep[i]
        memcpy(ctx.tmprow, ctx.pebits + <size_t> idx * <size_t> we,
               <size_t> we * sizeof(uint64_t))
        pivot = _reduce_row(ctx, ctx.tmprow)
        if pivot >= 0:
            _basis_store(ctx, pivot, ctx.tmprow)
            ctx.basis_idx[ctx.basis_count] = idx
            ctx.basis_count += 1
    return 0


cdef int _write_ring_segment(Structure structure, ring_ctx_t *ctx) except -1:
    """SEG_RELEVANT_RINGS as [count][offsets...][sentinel][indices...]."""
    cdef uint32_t count = <uint32_t> ctx.basis_count
    cdef uint32_t total = 0
    cdef uint32_t *out
    cdef uint32_t i, pos
    cdef Py_ssize_t j, ofs, idx
    if count == 0:
        return 0
    for j in range(ctx.basis_count):
        total += ctx.psize[ctx.basis_idx[j]]
    structure_append(structure, SEG_RELEVANT_RINGS,
                     <size_t> (1 + count + 1 + total) * sizeof(uint32_t))
    out = structure_rings(structure)
    out[0] = count
    pos = 0
    for i in range(count):
        out[1 + i] = pos
        pos += ctx.psize[ctx.basis_idx[i]]
    out[1 + count] = pos
    pos = 2 + count
    for i in range(count):
        idx = ctx.basis_idx[i]
        ofs = ctx.pcyc_ofs[idx]
        for j in range(ctx.psize[idx]):
            out[pos] = ctx.pcyc[ofs + j]
            pos += 1
    return 0


cdef int _fill_descriptors(Structure structure, ring_ctx_t *ctx) except -1:
    """Per-atom ring sizes, ring bitmap and ring counts, one bit per relevant-cycle prototype.

    The bitmap is prototype-scoped, not basis-scoped: `shares_ring(a, b)` asks whether some
    relevant prototype covers both atoms. That is exact whenever a prototype generates a single
    cycle -- every fused, spiro and cage system in practice -- and conservative only where two
    atoms sit on arcs of one prototype that no single cycle of it uses together, as in a
    macrocyclic cyclophane. Scoping it to the stored basis instead would be worse: which of a
    cage's equivalent faces the basis drops is an artefact of the greedy order, so cubane atoms
    that plainly share a face would answer False.

    Ring sizes and counts likewise derive from the prototypes rather than the basis. The basis
    is one cycle short of the full relevant set on nearly every polycycle -- C60 has 32 faces at
    circuit rank 31 -- so reading descriptors off the basis would lose a real ring.
    """
    cdef uint32_t n = ctx.n
    cdef Py_ssize_t count = ctx.keep_count
    cdef uint32_t words = <uint32_t> ((count + 63) // 64)
    cdef Py_ssize_t wv = ctx.wv
    cdef Py_ssize_t r, w
    cdef atom_t *atoms
    cdef uint64_t *bits
    cdef uint64_t *pv
    cdef uint64_t word
    cdef uint32_t v, k, size, total
    if words == 0 or n == 0:
        return 0
    structure_append(structure, SEG_RING_BITS,
                     <size_t> n * <size_t> words * sizeof(uint64_t))
    bits = structure_ring_bits(structure)
    atoms = structure.atoms()          # structure_append may have moved the buffer
    with nogil:
        for r in range(count):
            size = ctx.psize[ctx.keep[r]]
            pv = ctx.pvbits + <size_t> ctx.keep[r] * <size_t> wv
            for w in range(wv):
                word = pv[w]
                while word:
                    v = <uint32_t> (w * 64 + _lo_bit64(word))
                    word &= word - 1
                    at_add_ring_size(&atoms[v], size)
                    bits[<size_t> v * <size_t> words + (r >> 6)] |= <uint64_t> 1 << (r & 63)
        for v in range(n):
            total = 0
            for k in range(words):
                total += _popcount64(bits[<size_t> v * <size_t> words + k])
            if total > 255:
                total = 255
            at_set_ring_counts(&atoms[v], <uint8_t> total, 0)
    return 0


cdef int _raise_rc(int rc, uint32_t n) except -1:
    if rc == 0:
        return 0
    if rc == -1:
        raise MemoryError('ring perception scratch allocation failed')
    if rc == -2:
        raise ValueError('ring perception exceeded the %d relevant-cycle prototype limit on a '
                         '%d-atom structure' % (int(PROTOTYPE_LIMIT), int(n)))
    raise ValueError('ring perception exceeded its %d s deadline on a %d-atom structure'
                     % (int(DEADLINE_SECONDS), int(n)))


cdef int perceive_rings(Structure structure) except -1:
    """Fill every ring descriptor: per-atom sizes and counts, the bitmap, the cycle basis.

    The clearing pass runs unconditionally and before every early return: `_apply` memcpys each
    atom_t field forward, so a field a perception pass may not write has to be cleared here or
    it survives from the pre-edit structure.
    """
    cdef uint32_t n = structure.header.atom_count
    cdef atom_t *atoms = structure.atoms()
    cdef atom_t *a
    cdef uint32_t v
    cdef int rc = 0
    cdef ring_ctx_t ctx
    with nogil:
        for v in range(n):
            a = &atoms[v]
            a.ring_sizes = 0
            a.ring_counts = 0
    if n == 0 or structure.header.bond_count == 0:
        return 0

    memset(&ctx, 0, sizeof(ring_ctx_t))
    try:
        with nogil:
            ctx.deadline = time(NULL) + <time_t> DEADLINE_SECONDS
            rc = _prepare_ctx(&ctx, structure)
            if rc == 0 and ctx.ne:
                rc = _build_prototypes(&ctx)
                if rc == 0:
                    rc = _filter_relevant(&ctx)
                if rc == 0:
                    rc = _select_basis(&ctx)
        _raise_rc(rc, n)
        if ctx.keep_count:
            _write_ring_segment(structure, &ctx)
            _fill_descriptors(structure, &ctx)
    finally:
        _free_ctx(&ctx)
    return 0


cdef int mark_bridges(Structure structure) noexcept nogil:
    """Flag every half-edge that lies on a cycle, and every atom that carries one.

    Order-8 bonds are not edges here. A dative bond is a coordination arrow, not a ring
    closure, and admitting it destroys the ring it appears to create: ferrocene's Fe-Cp
    fragment is a wheel, whose five triangles weigh 15 against four triangles plus the Cp
    five-ring at 17, so a minimum cycle basis drops the Cp ring -- and since the ring is in no
    minimum basis it is not a relevant cycle either, so it vanishes from ring_sizes_of as
    well. No downstream algorithm can recover it. Excluding the bond here is what keeps it.

    This is the only place the exclusion is needed: `_prepare_ctx` reads only half-edges with
    HE_IN_RING set, so prototypes, the basis and every per-atom descriptor inherit it.
    """
    cdef uint32_t n = structure.header.atom_count
    if n == 0:
        return 0
    cdef uint32_t *ptr = csr_ptr(structure)
    cdef halfedge_t *edges = csr_edges(structure)
    cdef atom_t *atoms = structure.atoms()

    cdef halfedge_t *e
    cdef uint32_t h
    for h in range(2 * structure.header.bond_count):
        e = &edges[h]
        e.flags = <uint16_t> (e.flags & ~HE_IN_RING)

    cdef uint32_t *disc = <uint32_t *> malloc(n * sizeof(uint32_t))
    cdef uint32_t *low = <uint32_t *> malloc(n * sizeof(uint32_t))
    # explicit DFS stack: vertex, next half-edge to examine, parent half-edge index
    cdef uint32_t *st_v = <uint32_t *> malloc(n * sizeof(uint32_t))
    cdef uint32_t *st_k = <uint32_t *> malloc(n * sizeof(uint32_t))
    cdef uint32_t *st_p = <uint32_t *> malloc(n * sizeof(uint32_t))
    if disc is NULL or low is NULL or st_v is NULL or st_k is NULL or st_p is NULL:
        free(disc); free(low); free(st_v); free(st_k); free(st_p)
        return -1

    cdef uint32_t NONE = 0xffffffff
    cdef bint any_ring
    cdef uint32_t i, root, timer = 0, top, v, k, child, parent, twin
    for i in range(n):
        disc[i] = NONE
        low[i] = NONE

    for root in range(n):
        if disc[root] != NONE:
            continue
        st_v[0] = root
        st_k[0] = ptr[root]
        st_p[0] = NONE
        disc[root] = timer
        low[root] = timer
        timer += 1
        top = 0
        while True:
            v = st_v[top]
            k = st_k[top]
            if k < ptr[v + 1]:
                st_k[top] = k + 1
                if k == st_p[top]:
                    continue                       # never walk back up the tree edge
                e = &edges[k]
                if e.order == 8:
                    continue                       # dative bonds are not ring edges
                child = e.to
                if disc[child] == NONE:            # tree edge: descend
                    twin = ptr[child]
                    while twin < ptr[child + 1] and edges[twin].to != v:
                        twin += 1
                    top += 1
                    st_v[top] = child
                    st_k[top] = ptr[child]
                    st_p[top] = twin
                    disc[child] = timer
                    low[child] = timer
                    timer += 1
                else:                              # back edge: on a cycle by definition
                    _mark_ring(edges, ptr, v, k, child)
                    if disc[child] < low[v]:
                        low[v] = disc[child]
            else:
                if top == 0:
                    break
                top -= 1
                parent = st_v[top]
                if low[v] < low[parent]:
                    low[parent] = low[v]
                if low[v] <= disc[parent]:         # tree edge parent-v is not a bridge
                    twin = st_p[top + 1]           # the child -> parent half-edge
                    edges[twin].flags |= HE_IN_RING
                    _mark_twin(edges, ptr, parent, v)

    for v in range(n):
        any_ring = False
        for k in range(ptr[v], ptr[v + 1]):
            if edges[k].flags & HE_IN_RING:
                any_ring = True
                break
        at_set_in_ring(&atoms[v], any_ring)

    free(disc); free(low); free(st_v); free(st_k); free(st_p)
    return 0


cdef inline void _mark_twin(halfedge_t *edges, uint32_t *ptr, uint32_t frm,
                            uint32_t to) noexcept nogil:
    cdef uint32_t k
    for k in range(ptr[frm], ptr[frm + 1]):
        if edges[k].to == to:
            edges[k].flags |= HE_IN_RING
            return


cdef inline void _mark_ring(halfedge_t *edges, uint32_t *ptr, uint32_t v, uint32_t k,
                            uint32_t child) noexcept nogil:
    edges[k].flags |= HE_IN_RING
    _mark_twin(edges, ptr, child, v)
