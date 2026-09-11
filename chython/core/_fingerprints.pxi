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
# FINGERPRINTS, as the product of two independent things: a per-atom LABEL and a FRAGMENT
# ENUMERATOR.  The enumerators here know nothing about chemistry -- they walk the CSR, read a
# `uint32[n]` vector of labels, and return the unfolded truth as `dict[int, int]`, a 64-bit fragment
# hash to how many times the molecule carries it.  Folding is one shared step layered on that dict.
#
# WHY THE UNFOLDED DICT IS THE ONLY THING AN ENUMERATOR RETURNS.  Every folded spelling is derived
# from it in three lines, so a folded and an unfolded answer cannot disagree.  A folded-only count
# approximation -- one extra bit lit per repeat -- has no unfolded spelling that says the same thing.
#
# WHY THE LABEL VECTOR IS A PUBLIC ARGUMENT.  `invariants=` is the whole expansion story: FCFP is
# Morgan over pharmacophore labels and costs no code here, SiRMS is one new enumerator over labels
# that already exist.  A caller with its own atom typing passes its own array, with no subclassing
# and no injection.
#
# INCLUDED AFTER `_molecule_container.pxi` and `_molecule_topology.pxi`: the wrappers in the class
# call forward into these `cdef` functions, which RULES.md §7.1 says is free, and `_numpy_load` and
# `_NP_*` are bound in the topology fragment.


# --- the featurization atom invariant -------------------------------------------------------------
#
# NOT `_morgan.pxi::_atom_invariant`, and the reason is content rather than taste.  That one has no
# degree term, spends sixteen bits on the isotope, keeps the implicit and explicit hydrogen nibbles
# apart, and yields a packed key rather than a mixed hash -- every one of which is right for a
# canonical refinement class and wrong here, where the two hydrogen counts must SUM (`C` and
# `[H]C([H])([H])[H]` are one molecule) and the result must fold.  The names are close enough to
# confuse: `fp_atom_invariant` serves featurization, `_atom_invariant` serves canonical refinement,
# and neither may be swapped for the other.

DEF FP_H_MAX = H_IMPLICIT_MAX + H_EXPLICIT_MAX  # total hydrogen count; five bits hold it with two to spare

# One seed per hashed domain, so a fragment of one kind can never collide with a fragment of
# another at equal content.  ASCII, for no reason beyond being greppable.
DEF FP_SEED_ATOM  = 0x63686641   # 'chfA' -- one atom's label
DEF FP_SEED_SHELL = 0x63686653   # 'chfS' -- one Morgan shell
DEF FP_SEED_PATH  = 0x63686650   # 'chfP' -- one linear path


cdef inline uint32_t fp_atom_invariant(atom_t *a) noexcept nogil:
    """One atom's default label: element, isotope, charge, radical, total H, degree, in-ring flag.

    Packed into a word and then mixed down to 32 bits, because a label is fed to a fragment hash and
    a packed key's low bits are its in-ring flag -- fine for equality, useless under a shift.

    AN UNKNOWN HYDROGEN COUNT READS AS ZERO, and this is the one place in the tree where the sentinel
    collapses rather than propagating.  RULES.md §6.4 -- if "I mean zero" and "nobody knows" produce
    one value the encoding is wrong -- governs a STORED field, where the two facts must stay tellable
    apart.  This is not one.  A fingerprint is a lossy screen whose purpose is to bring near-matches
    together, and the case that decides it is retrieval: an MDL record of a Suzuki palladium catalyst
    whose hydrogens nobody could derive has to screen against the curated form of the same catalyst,
    where they are zero.  Splitting them would break the screen at the job it is wanted for, and the
    sentinel almost always becomes zero once the record is repaired anyway.  `_atom_invariant` keeps
    its sentinel for the opposite reason -- a refinement class is a claim about indistinguishability,
    and a missing number distinguishes.  Two answers, two questions, and neither is a stale copy.

    So `FP_H_MAX` is `H_IMPLICIT_MAX + H_EXPLICIT_MAX` and the five-bit field holds nothing else.

    DEGREE IS THE HEAVY-ATOM DEGREE, `a.degree - at_explicit_h(a)`, and the subtraction is what makes
    the hydrogen claim above true.  `read_smiles('C')` gives the carbon `implicit_h=4, degree=0`;
    `read_smiles('[H]C([H])([H])[H]')` gives it `explicit_h=4, degree=4`.  The totals already agree,
    so leaving the raw CSR row length in would split one molecule into two labels over nothing but
    how it was written -- and `explicit_h` is by construction the count of hydrogen neighbours, so it
    is exactly the term to remove.  A DATIVE BOND STAYS COUNTED: a coordination contact is a real
    structural difference, not a difference in spelling, which is where this parts company with the
    SMARTS `D` primitive.  A full byte, so no bound to declare and no saturation to get wrong.
    """
    cdef uint64_t h
    cdef uint64_t v
    if at_implicit_h_unknown(a):
        h = <uint64_t> at_explicit_h(a)          # the sentinel is not a count; read it as no implicit
    else:
        h = <uint64_t> at_implicit_h(a) + <uint64_t> at_explicit_h(a)
    v = <uint64_t> a.element
    v = (v << 16) | <uint64_t> a.isotope
    v = (v << 4) | <uint64_t> (<int> a.charge - CHARGE_MIN)
    v = (v << 1) | (<uint64_t> 1 if at_radical(a) else <uint64_t> 0)
    v = (v << 5) | h
    v = (v << 8) | <uint64_t> (a.degree - at_explicit_h(a))   # heavy-atom degree; see above
    v = (v << 1) | (<uint64_t> 1 if at_in_ring(a) else <uint64_t> 0)
    return <uint32_t> _xxh64(&v, 1, FP_SEED_ATOM)


cdef object fp_atom_invariants(Structure structure):
    """`(n,)` uint32 of the default label, one per atom in the molecule's own atom order."""
    _numpy_load()
    cdef uint32_t n_atoms = structure.header.atom_count
    cdef object out = _NP_EMPTY(n_atoms, dtype='uint32')
    if n_atoms == 0:
        return out
    cdef uint32_t[::1] labels = out
    cdef atom_t *atoms = structure.atoms()
    cdef uint32_t n
    with nogil:
        for n in range(n_atoms):
            labels[n] = fp_atom_invariant(&atoms[n])
    return out


cdef object _fp_labels(Structure structure, object invariants):
    """The label vector an enumerator reads: the default one, or the caller's after validation.

    RETURNS THE ARRAY AND NOT A POINTER.  The caller keeps it alive for the length of the walk; a
    `uint32_t *` handed out from here would outlive the buffer it points at the moment the temporary
    was dropped.

    A NON-`uint32` VECTOR IS REFUSED RATHER THAN CONVERTED.  `astype` on a float array truncates and
    on a negative int wraps, and both are a caller's mistake worth hearing about -- so the message
    says what to pass instead.
    """
    if invariants is None:
        return fp_atom_invariants(structure)
    _numpy_load()
    cdef uint32_t n_atoms = structure.header.atom_count
    cdef object out = _NP_ASCONTIGUOUS(invariants)
    if out.dtype.name != 'uint32':
        raise ValueError(f'invariants must be a uint32 vector, got {out.dtype.name}; '
                         f'`atom_invariants()` returns the right dtype, and a caller building its '
                         f"own scheme should say dtype='uint32'")
    if out.ndim != 1:
        raise ValueError(f'invariants must be a one-dimensional vector, '
                         f'got {out.ndim} dimensions')
    if out.shape[0] != n_atoms:
        raise ValueError(f'invariants needs one entry per atom, in this molecule\'s atom order: '
                         f'expected {n_atoms}, got {out.shape[0]}')
    return out


# --- argument validation --------------------------------------------------------------------------

cdef int _fp_check_radii(int min_radius, int max_radius) except -1:
    """Both radius bounds, refused rather than clamped.

    An unvalidated `min_radius=0` behaves as 2 and drops every singleton fragment -- a wrong
    fingerprint with no way to notice.
    """
    if min_radius < 1:
        raise ValueError(f'min_radius must be at least 1, got {min_radius}: radius 1 is the atom '
                         f'by itself and there is no smaller fragment')
    if max_radius < min_radius:
        raise ValueError(f'max_radius must be at least min_radius, got max_radius={max_radius} '
                         f'below min_radius={min_radius}')
    return 0


# --- morgan: circular fragments -------------------------------------------------------------------

cdef dict fp_morgan_counts(Structure structure, uint32_t min_radius, uint32_t max_radius,
                           object invariants):
    """Circular fragments of radius `min_radius..max_radius`, as hash to how many atoms carry it.

    RADIUS COUNTS SHELLS AND STARTS AT 1: radius 1 is the bare atom label, radius 2 the atom plus its
    bonded neighbours, and radius `r` costs `r - 1` expansions.

    One shell's identifier is the hash of the centre's previous identifier followed by its
    neighbours', SORTED, so the result cannot depend on CSR order and therefore cannot depend on the
    order the molecule was written in.  The bond order is mixed into each neighbour's word rather
    than stored beside it: one word per neighbour is what lets `_sort_words` canonicalise the
    multiset in a single pass, and a mix that can collide is no worse than the hash it feeds.

    EVERY SHELL IS HASHED, RADIUS 1 INCLUDED.  Returning the 32-bit label unmixed at radius 1 would
    leave a fingerprint's high bits constant there, and section 3.2's `number_active_bits * log2
    (length) <= 64` rule is a promise about all 64.
    """
    cdef object arr = _fp_labels(structure, invariants)   # validated even for an empty molecule
    cdef dict counts = {}
    cdef uint32_t n_atoms = structure.header.atom_count
    cdef const uint32_t[::1] labels
    cdef const uint32_t *lab
    cdef atom_t *atoms
    cdef uint32_t *ptr
    cdef halfedge_t *edges
    cdef halfedge_t *e
    cdef uint32_t n, k, radius, deg, width = 0
    cdef uint64_t *swap
    cdef uint64_t *block
    cdef uint64_t *cur
    cdef uint64_t *nxt
    cdef uint64_t *words

    if n_atoms == 0:
        return counts
    labels = arr
    lab = &labels[0]
    atoms = structure.atoms()
    ptr = csr_ptr(structure)
    edges = csr_edges(structure)
    # THE RAW CSR ROW LENGTH, not the heavy-atom degree `fp_atom_invariant` uses.  This sizes the
    # word buffer the walk writes into, and the walk visits every half-edge including one to a
    # hydrogen atom -- subtracting `at_explicit_h` here would undersize the buffer and corrupt memory.
    for n in range(n_atoms):
        if atoms[n].degree > width:
            width = atoms[n].degree

    # one block, three regions (RULES.md §5.2): this shell's identifiers, the next shell's, and the
    # word buffer one atom's hash is built in -- `width + 1` for the centre plus its neighbours.
    block = <uint64_t *> PyMem_Malloc(<size_t> (2 * n_atoms + width + 1) * sizeof(uint64_t))
    if block is NULL:
        raise MemoryError('morgan fingerprint scratch allocation failed')
    cur = block
    nxt = block + n_atoms
    words = block + 2 * n_atoms
    try:
        with nogil:
            for n in range(n_atoms):
                words[0] = <uint64_t> lab[n]
                cur[n] = _xxh64(words, 1, FP_SEED_SHELL)
        for radius in range(1, max_radius + 1):
            if radius > 1:
                with nogil:
                    for n in range(n_atoms):
                        words[0] = cur[n]
                        deg = 0
                        for k in range(ptr[n], ptr[n + 1]):
                            e = edges + k
                            deg += 1
                            words[deg] = cur[e.to] ^ (<uint64_t> e.order * <uint64_t> XXH_P1)
                        _sort_words(words + 1, deg)
                        nxt[n] = _xxh64(words, deg + 1, FP_SEED_SHELL)
                swap = cur
                cur = nxt
                nxt = swap
            if radius >= min_radius:
                # the GIL is back for exactly this: a dict update, once per atom per kept shell
                for n in range(n_atoms):
                    counts[cur[n]] = counts.get(cur[n], 0) + 1
    finally:
        PyMem_Free(block)
    return counts


# --- the folder, shared by every family -----------------------------------------------------------
#
# A hash contributes the low `log2(length)` bits as its first position, then the next `log2(length)`
# for a second active bit, and so on.  All three spellings walk the same positions, which is the
# point: `bit_set`, `fingerprint().nonzero()` and `count_vector().nonzero()` are the same answer in
# three shapes, and there is one place where that could ever stop being true.

cdef inline uint32_t _fp_log2(uint32_t length) noexcept nogil:
    """`log2` of a power of two, by shifting.  The caller has already checked that it is one."""
    cdef uint32_t width = 0
    while (<uint32_t> 1 << width) < length:
        width += 1
    return width


cdef inline uint32_t _fp_position(uint64_t h, uint32_t i, uint32_t width, uint64_t mask) noexcept nogil:
    """Bit position for the i-th active bit of hash `h`, given a window of `width` bits and `mask`.

    ONE PLACE FOR ONE EXPRESSION.  All three fold functions read the same `log2(length)`-bit window
    of `h` at offset `i * width`; extracting it here means the three spellings of the same answer
    are structurally identical rather than agreeing only by test.
    """
    return <uint32_t> ((h >> (i * width)) & mask)


cdef int _fp_check_folding(int length, int number_active_bits) except -1:
    """The folding arguments, refused rather than silently truncated.

    THE ACTIVE-BIT BOUND IS ABOUT THE HASH'S WIDTH.  A fragment hash is 64 bits, so `n` slices of
    `log2(length)` bits each need `n * log2(length) <= 64`; past that the folder has shifted the hash
    away entirely and every further bit is the constant zero.  `length=1024, number_active_bits=7`
    is the first pair over the bound: accepting it returns a fingerprint whose last bits are all
    zero, wasting a tenth of the vector and biasing every similarity computed from it.
    """
    if length < 2 or (length & (length - 1)):
        raise ValueError(f'length must be a power of two and at least 2, got {length}')
    if number_active_bits < 1:
        raise ValueError(f'number_active_bits must be at least 1, got {number_active_bits}')
    cdef uint32_t width = _fp_log2(<uint32_t> length)
    if (<uint32_t> number_active_bits) * width > 64:
        raise ValueError(
            f'number_active_bits={number_active_bits} needs '
            f'{number_active_bits * <int> width} bits of hash at length={length}, and a fragment '
            f'hash is 64 bits wide: every bit past the {64 // <int> width}th would be a constant')
    return 0


cdef set fp_fold_bit_set(dict counts, uint32_t length, uint32_t active):
    """The folded bit POSITIONS, `0 <= p < length`.  A count is read as presence and nothing more."""
    cdef uint32_t width = _fp_log2(length)
    cdef uint64_t mask = <uint64_t> length - 1
    cdef uint64_t h
    cdef uint32_t i
    cdef set out = set()
    for h in counts:
        for i in range(active):
            out.add(_fp_position(h, i, width, mask))
    return out


cdef object fp_fold_binary(dict counts, uint32_t length, uint32_t active):
    """`(length,)` uint8 of 0 and 1: one where some fragment folded onto that position."""
    _numpy_load()
    cdef object out = _NP_ZEROS(length, dtype='uint8')
    cdef uint8_t[::1] bits = out
    cdef uint32_t width = _fp_log2(length)
    cdef uint64_t mask = <uint64_t> length - 1
    cdef uint64_t h
    cdef uint32_t i
    for h in counts:
        for i in range(active):
            bits[_fp_position(h, i, width, mask)] = 1
    return out


cdef object fp_fold_counted(dict counts, uint32_t length, uint32_t active):
    """`(length,)` uint32: every fragment adds its FULL count to every position it activates.

    A COLLISION ADDS, IT DOES NOT REPLACE, which is what makes the vector's total equal the unfolded
    total times `active` -- the property that says folding lost no weight, only resolution.

    SATURATES AT `uint32` MAX RATHER THAN WRAPPING.  A wrap would turn the largest count in the
    molecule into the smallest, and a saturated maximum is at least monotone.  Reaching it needs four
    billion copies of one fragment, so this is a correctness floor and not a live case -- the branch
    is deliberately untested; see the note next to `test_the_count_vector_carries_at_least_as_much_
    weight_as_the_binary_one` in `test_fingerprints.py`.
    """
    _numpy_load()
    cdef object out = _NP_ZEROS(length, dtype='uint32')
    cdef uint32_t[::1] vec = out
    cdef uint32_t width = _fp_log2(length)
    cdef uint64_t mask = <uint64_t> length - 1
    cdef uint64_t h, c, room
    cdef uint32_t i, p
    for h, c in counts.items():
        for i in range(active):
            p = _fp_position(h, i, width, mask)
            room = <uint64_t> 0xffffffff - <uint64_t> vec[p]
            vec[p] += <uint32_t> (c if c < room else room)
    return out


# --- linear: simple paths -------------------------------------------------------------------------

cdef inline uint64_t _fp_path_hash(uint32_t *path, uint8_t *orders, uint32_t depth,
                                   const uint32_t *lab, uint64_t *words) noexcept nogil:
    """Hash one simple path, read from whichever end gives the larger word sequence.

    THE END IS CHOSEN BY LABELS, NEVER BY ATOM INDEX.  Two paths through different atoms carrying the
    same labels and orders must land on one hash -- that is the whole point of a fragment key -- and
    an index tie-break would split them.

    The forward sequence is written first and reversed in place when the other end wins.  `count` is
    odd and labels sit on the even positions, so a reversal keeps labels and bond orders in their own
    slots and no second buffer is needed.
    """
    cdef uint32_t count = 2 * depth + 1
    cdef uint32_t i, m
    cdef uint64_t t
    for i in range(depth + 1):
        words[2 * i] = <uint64_t> lab[path[i]]
    for i in range(1, depth + 1):
        words[2 * i - 1] = <uint64_t> orders[i]
    for i in range(count // 2):
        if words[i] != words[count - 1 - i]:
            if words[i] < words[count - 1 - i]:
                for m in range(count // 2):
                    t = words[m]
                    words[m] = words[count - 1 - m]
                    words[count - 1 - m] = t
            break
    return _xxh64(words, count, FP_SEED_PATH)


cdef dict fp_linear_counts(Structure structure, uint32_t min_radius, uint32_t max_radius,
                           object invariants):
    """Every simple path of `min_radius..max_radius` ATOMS, as hash to how many paths carry it.

    THE RADII COUNT ATOMS, NOT BONDS: length 1 is a lone atom label, length 2 is a bond, length `r`
    spans `r - 1` bonds.  So `sum(counts.values())` at length 2 is exactly `bond_count`.  `min_radius`
    defaults to 1, so this surface takes the same arguments as the Morgan family.

    EACH UNDIRECTED PATH IS COUNTED ONCE.  A depth-first walk from every atom finds each path from
    both ends, so one direction is chosen -- `path[0] < path[depth]` on the dense index -- and the
    other dropped.  Deduplicating afterwards through a set of tuples is the same answer and a great
    deal more allocation.

    THE WALK IS EXPONENTIAL IN `max_radius` on a densely fused ring system, and there is no cap: a
    silent truncation would be a wrong fingerprint with nothing to notice it by.  The default of 4 is
    cheap everywhere; a caller asking for 20 on a steroid is asking for the walk it gets.

    THE GIL IS HELD FOR THE WALK.  A path is counted as it is found and the counter is a dict; the
    alternative -- a growable hash buffer drained afterwards -- buys nothing at these sizes and adds
    a resize path to get wrong.
    """
    cdef object arr = _fp_labels(structure, invariants)   # validated even for an empty molecule
    cdef dict counts = {}
    cdef uint32_t n_atoms = structure.header.atom_count
    cdef const uint32_t[::1] labels
    cdef const uint32_t *lab
    cdef uint32_t *ptr
    cdef halfedge_t *edges
    cdef halfedge_t *e
    cdef uint32_t start, depth, k, to, walk_depth
    cdef uint64_t key
    cdef size_t size_words, size_path, size_cursor, size_orders
    cdef char *block
    cdef uint64_t *words
    cdef uint32_t *path
    cdef uint32_t *cursor
    cdef uint8_t *orders
    cdef uint8_t *on_path

    if n_atoms == 0:
        return counts
    if max_radius == 0:
        # chosen: raise (not clamp to empty) -- max_radius == 0 has no meaning and silence would
        # hide a direct caller that forgot _fp_check_radii.  fp_linear_counts holds the GIL so
        # raising is clean; the guard also prevents 2*max_radius-1 underflowing the scratch size.
        raise ValueError('max_radius must be at least 1 in fp_linear_counts')
    labels = arr
    lab = &labels[0]
    ptr = csr_ptr(structure)
    edges = csr_edges(structure)
    # A simple path visits distinct atoms, so no path can exceed n_atoms atoms.  Clamp the scratch
    # to the reachable bound: linear_hash_counts(1, 2**31-1) on ethanol must not reserve 34 GB.
    walk_depth = max_radius if max_radius < n_atoms else n_atoms

    # one block, five regions (RULES.md §5.2).  THE ORDER IS THE ALIGNMENT: the uint64 words come
    # first at offset zero, the uint32 stack arrays next at a multiple of eight, and the two byte
    # arrays last where alignment cannot be violated.  walk_depth is the clamped maximum path
    # length; max_radius is kept in the loop condition below for the case where it is smaller.
    size_words = <size_t> (2 * walk_depth - 1) * sizeof(uint64_t)
    size_path = <size_t> walk_depth * sizeof(uint32_t)
    size_cursor = <size_t> walk_depth * sizeof(uint32_t)
    size_orders = <size_t> walk_depth * sizeof(uint8_t)
    block = <char *> PyMem_Malloc(size_words + size_path + size_cursor + size_orders
                                  + <size_t> n_atoms * sizeof(uint8_t))
    if block is NULL:
        raise MemoryError('linear fingerprint scratch allocation failed')
    words = <uint64_t *> block
    path = <uint32_t *> (block + size_words)
    cursor = <uint32_t *> (block + size_words + size_path)
    orders = <uint8_t *> (block + size_words + size_path + size_cursor)
    on_path = <uint8_t *> (block + size_words + size_path + size_cursor + size_orders)
    try:
        for k in range(n_atoms):                # zero once; backtracking unstamps as it goes
            on_path[k] = 0
        for start in range(n_atoms):
            depth = 0
            path[0] = start
            orders[0] = 0                       # no bond precedes the first atom
            cursor[0] = ptr[start]
            on_path[start] = 1
            if min_radius <= 1:
                key = _fp_path_hash(path, orders, 0, lab, words)
                counts[key] = counts.get(key, 0) + 1
            while True:
                if depth + 1 < walk_depth and cursor[depth] < ptr[path[depth] + 1]:
                    k = cursor[depth]
                    cursor[depth] += 1
                    e = edges + k
                    to = e.to
                    if on_path[to]:
                        continue
                    depth += 1
                    path[depth] = to
                    orders[depth] = e.order
                    cursor[depth] = ptr[to]
                    on_path[to] = 1
                    # one direction per undirected path; the ends of a simple path never coincide
                    if depth + 1 >= min_radius and path[0] < path[depth]:
                        key = _fp_path_hash(path, orders, depth, lab, words)
                        counts[key] = counts.get(key, 0) + 1
                elif depth == 0:
                    on_path[start] = 0          # only stamp remaining; deeper atoms backtracked
                    break
                else:
                    on_path[path[depth]] = 0
                    depth -= 1
    finally:
        PyMem_Free(block)
    return counts
