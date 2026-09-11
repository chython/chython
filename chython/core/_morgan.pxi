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
# Atom ordering: refinement of the graph to an equitable partition.
#
# One key per atom, one sort, one hash. Round k hashes each atom's class together with its
# neighbours' classes and bond orders, sorted so the neighbour order the CSR happens to hold
# cannot leak in; atoms are then grouped by hash into new classes. Classes only ever split,
# because an atom's own class seeds its hash, so the count is monotone and the fixed point
# arrives in at most n rounds.
#
# Hashing is safe here in a way that is worth stating: equal inputs always hash equal, so two
# symmetry-equivalent atoms can never be driven into different classes. A 64-bit collision can
# only merge two classes that a previous round had separated -- it makes the partition coarser,
# never wrong, and coarser is the harmless direction, since a coarse partition just claims less.
# The loop watches the class count anyway and stops on any drop, keeping the finer partition, so
# a collision costs discrimination on that one molecule and nothing else. Keeping the degraded
# result and logging a note about it is the alternative, and it is a worse answer.
#
# The starting class packs the atom invariant into one uint64: element, isotope, charge, radical,
# implicit hydrogens, ring membership. Ring membership is largely re-derivable by refinement and is
# kept because it belongs to the invariant, not for reach.
# Everything is fixed-constant integer arithmetic, so the result is reproducible across runs,
# platforms and Python versions.


# XXH64's round function, tail mixer and avalanche, specialised to whole 8-byte lanes. This is
# xxhash's mixing, not a substitute for it; only the byte-tail cases it can never reach are gone.
DEF XXH_P1 = 0x9E3779B185EBCA87
DEF XXH_P2 = 0xC2B2AE3D27D4EB4F
DEF XXH_P3 = 0x165667B19E3779F9
DEF XXH_P4 = 0x85EBCA77C2B2AE63
DEF XXH_P5 = 0x27D4EB2F165667C5


cdef inline uint64_t _rotl64(uint64_t x, int r) noexcept nogil:
    return (x << r) | (x >> (64 - r))


cdef inline uint64_t _xxh64(uint64_t *words, uint32_t count, uint64_t seed) noexcept nogil:
    cdef uint64_t h = seed + <uint64_t> XXH_P5 + 8 * <uint64_t> count
    cdef uint64_t lane
    cdef uint32_t i
    for i in range(count):
        lane = _rotl64(words[i] * <uint64_t> XXH_P2, 31) * <uint64_t> XXH_P1
        h = _rotl64(h ^ lane, 27) * <uint64_t> XXH_P1 + <uint64_t> XXH_P4
    h ^= h >> 33
    h *= <uint64_t> XXH_P2
    h ^= h >> 29
    h *= <uint64_t> XXH_P3
    h ^= h >> 32
    return h


cdef inline uint64_t _atom_invariant(atom_t *a) noexcept nogil:
    """The atom's own starting class, packed into one word. Charge is biased by 4 to keep the
    field unsigned over the -4 .. +8 range atom_t allows.

    The implicit nibble goes in RAW, sentinel included, and needs no H_UNKNOWN branch: this is a
    packed key that is only ever compared for equality, and "unknown" is a perfectly good fourteenth
    value to be equal or unequal to.  Two atoms with no recorded count land in one class, an atom
    with no count and an atom with three land in different ones, and both are the right answer -- a
    refinement class is a claim about indistinguishability, and a missing number distinguishes.
    """
    cdef uint64_t v = <uint64_t> a.element
    v = (v << 16) | <uint64_t> a.isotope
    v = (v << 4) | <uint64_t> (<int> a.charge + 4)
    v = (v << 1) | (<uint64_t> 1 if at_radical(a) else <uint64_t> 0)
    v = (v << 4) | <uint64_t> at_implicit_h(a)
    v = (v << 1) | (<uint64_t> 1 if at_in_ring(a) else <uint64_t> 0)
    # The R index is part of the intrinsic record, so it belongs in the word the identity bytes emit:
    # `R1` and `R2` are different markers, and a fragment dedup key that collapsed them would merge
    # two attachment patterns.  Zero for every real element, which leaves their invariants unchanged.
    v = (v << 8) | <uint64_t> at_r_index(a)
    return v


cdef inline void _sort_words(uint64_t *w, uint32_t count) noexcept nogil:
    """Insertion sort. Degrees are single digits; nothing cleverer pays off."""
    cdef uint64_t x
    cdef uint32_t i, j
    for i in range(1, count):
        x = w[i]
        j = i
        while j and w[j - 1] > x:
            w[j] = w[j - 1]
            j -= 1
        w[j] = x


cdef Py_ssize_t _classify(uint64_t *key, uint32_t *idx, uint32_t *tmp, uint32_t n,
                          uint32_t *out) noexcept nogil:
    """Group atoms by equal key, writing dense 1-based classes into `out`; return the count.

    Bottom-up merge sort over indices, so equal keys land adjacent and the grouping is one pass.
    """
    cdef Py_ssize_t total = n
    cdef Py_ssize_t width = 1, i, mid, right, li, ri, k, classes
    cdef uint32_t *src = idx
    cdef uint32_t *dst = tmp
    cdef uint32_t *swap
    for i in range(total):
        idx[i] = <uint32_t> i
    while width < total:
        i = 0
        while i < total:
            mid = i + width
            if mid > total:
                mid = total
            right = i + 2 * width
            if right > total:
                right = total
            li = i
            ri = mid
            k = i
            while li < mid and ri < right:
                if key[src[ri]] < key[src[li]]:
                    dst[k] = src[ri]
                    ri += 1
                else:
                    dst[k] = src[li]
                    li += 1
                k += 1
            while li < mid:
                dst[k] = src[li]
                li += 1
                k += 1
            while ri < right:
                dst[k] = src[ri]
                ri += 1
                k += 1
            i += 2 * width
        swap = src
        src = dst
        dst = swap
        width *= 2
    classes = 0
    for i in range(total):
        if i == 0 or key[src[i]] != key[src[i - 1]]:
            classes += 1
        out[src[i]] = <uint32_t> classes
    return classes


cdef Py_ssize_t compute_atoms_order(Structure structure, uint32_t *rank,
                                    uint32_t *seed) noexcept nogil:
    """Write each atom's 1-based class into `rank`; return the class count, -1 on failure.

    Atoms sharing a class are the ones refinement cannot tell apart. Ranks are dense and ordered
    by the round-0 invariant, but rounds after that order by hash, so only the partition carries
    meaning -- not which class got which number.

    `seed` is NULL for the plain case, where the starting classes come from the atom records.
    Pass per-atom starting labels instead to refine a partition the caller already has: stereo
    needs exactly that, re-refining as stereocentre differentiation feeds it new distinctions.
    Labels need not be dense or 1-based; only their equality classes are read.
    """
    cdef uint32_t n = structure.header.atom_count
    if n == 0:
        return 0

    cdef uint32_t *ptr = csr_ptr(structure)
    cdef halfedge_t *edges = csr_edges(structure)
    cdef atom_t *atoms = structure.atoms()
    cdef uint32_t v, u, k, d, maxdeg = 0
    cdef Py_ssize_t classes, newclasses, rounds

    for v in range(n):
        d = ptr[v + 1] - ptr[v]
        if d > maxdeg:
            maxdeg = d

    # one allocation carved into six spans: u64 first so the u32 spans stay aligned
    cdef char *block = <char *> malloc(<size_t> (n + maxdeg) * sizeof(uint64_t)
                                       + <size_t> 4 * n * sizeof(uint32_t))
    if block is NULL:
        return -1
    cdef uint64_t *key = <uint64_t *> block
    cdef uint64_t *nb = key + n
    cdef uint32_t *cur = <uint32_t *> (nb + maxdeg)
    cdef uint32_t *nxt = cur + n
    cdef uint32_t *idx = nxt + n
    cdef uint32_t *tmp = idx + n

    if seed is NULL:
        for v in range(n):
            key[v] = _atom_invariant(&atoms[v])
    else:
        for v in range(n):
            key[v] = <uint64_t> seed[v]
    classes = _classify(key, idx, tmp, n, cur)

    rounds = 0
    while classes < <Py_ssize_t> n and rounds < <Py_ssize_t> n:
        rounds += 1
        for v in range(n):
            d = 0
            for k in range(ptr[v], ptr[v + 1]):
                u = edges[k].to
                nb[d] = (<uint64_t> cur[u] << 8) | <uint64_t> edges[k].order
                d += 1
            _sort_words(nb, d)
            key[v] = _xxh64(nb, d, <uint64_t> cur[v])
        newclasses = _classify(key, idx, tmp, n, nxt)
        if newclasses <= classes:  # fixed point, or a collision merged classes: keep the finer
            break
        classes = newclasses
        memcpy(cur, nxt, <size_t> n * sizeof(uint32_t))

    memcpy(rank, cur, <size_t> n * sizeof(uint32_t))
    free(block)
    return classes
