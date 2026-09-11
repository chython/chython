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
# ML VIEWS: `TensorEncoding`, `mol_state_view`, `mol_transition_view`, `reaction_transition_view`.
# One `TensorEncoding` instance is built once and reused alongside a dataset; passing these as keyword
# arguments per structure would cost a kwargs dict on every molecule, which is real against a ~5 µs
# view budget.  AFTER the topology fragment, because `mol_state_view` calls `csr_bfs_all` from there.


# --- domains --------------------------------------------------------------------------------------
#
# Declared here and nowhere else (RULES.md §6.1).  A width is not a bound: `max_neighbors` fits a
# uint8 because the arena's degree is one, and 255 is the largest value it can express, not a chemical
# claim.  `unknown_h` tops out at 15 because that is H_UNKNOWN itself, the only 4-bit value that is a
# sentinel rather than a count.

DEF ML_H_MAX = 15                # inclusive; 15 IS H_UNKNOWN, so `unknown_h=15` is a pass-through
DEF ML_DEGREE_MAX = 255          # the arena's degree is a uint8
DEF ML_TOKEN_MIN = -2147483648   # int32 minimum; a token outside this range cannot round-trip
DEF ML_TOKEN_MAX = 2147483647    # int32 maximum


# --- vocabulary probe table -----------------------------------------------------------------------
#
# The key packs the five state numbers into one uint32: element 7 bits | h_before 4 | n_before 8 |
# h_after 4 | n_after 8 = 31 bits.  Degree gets 8 because the arena's is a uint8 and a metal centre
# can exceed 15; hydrogens get exactly 4, since the fourth bit pattern 15 is H_UNKNOWN itself.
#
# Open addressing, linear probing, power-of-two capacity at a load factor under 0.5.  1024 slots is
# 8 KB and L1-resident, so a lookup is typically one probe.  Key 0 marks an empty slot; a vocabulary
# entry for key 0 cannot be stored, so an atom whose packed key equals 0 (a lone R marker: element 0,
# h 0, degree 0) always reads `unknown`.  The vocabulary refuses that key explicitly.

DEF ML_ELEMENT_MAX = 127         # 7 bits
DEF ML_KEY_EMPTY = 0

cdef struct ml_vocab_t:
    uint32_t *keys
    int32_t *values
    uint32_t mask               # capacity - 1; capacity is a power of two


cdef inline uint32_t ML_VOCAB_KEY(uint32_t element, uint32_t hb, uint32_t nb,
                                  uint32_t ha, uint32_t na) noexcept nogil:
    return (element << 24) | (hb << 20) | (nb << 12) | (ha << 8) | na


cdef class TensorEncoding:
    """How a view's integers are laid out: shifts, clamps, padding and an optional vocabulary.

    BUILT ONCE AND REUSED.  Passing these as keyword arguments per structure costs a kwargs dict on
    every molecule, which is real against a ~5 µs view.  A caller keeps one instance beside its
    dataset.

    EVERY DEFAULT IS IDENTITY EXCEPT `unknown_h`.  `TensorEncoding()` yields physical values: atomic
    numbers, hydrogen counts, heavy degrees, and bond-count distances with -1 for a pair that has no
    path.  `unknown_h` defaults to 0 because a trained vocabulary keyed on `implicit_hydrogens or 0`
    misses every atom whose count is unstated if the key carries 15; `unknown_h=15` recovers the
    sentinel, and `ReactionContainer.modeling_view()` uses exactly that.

    ZERO IS OFF for `max_distance`, `max_neighbors`, `width` and `pad_diagonal`.  The first three have
    no meaningful zero: a clamp to zero flattens the column it clamps, and a width of zero admits no
    atom.  `pad_diagonal` is off for a different reason: a requested zero is indistinguishable from
    `pad` in the default `pad=0` configuration; the knob exists to place one non-masked cell on each
    fully padded row, so that an attention softmax over that row has at least one finite entry and does
    not produce NaN.  `pad_diagonal=1` is the canonical non-zero spelling.

    ORDER OF OPERATIONS, which is what a wrong reading gets backwards: the vocabulary key is built
    from physical unclamped values with `unknown_h` already applied, THEN the clamps apply, THEN the
    shifts.  Building the key post-clamp yields the unknown token for every atom in a structure that
    trips a clamp.

    TOKENS ARE NEVER SHIFTED.  The vocabulary owns its own id space.
    """
    cdef readonly int32_t element_shift, hydrogen_shift, neighbor_shift, distance_shift
    cdef readonly int32_t disconnected, unknown_h, max_distance, max_neighbors
    cdef readonly int32_t width, pad, pad_diagonal, unknown
    cdef readonly object vocabulary
    cdef ml_vocab_t _table

    def __init__(self, int32_t element_shift=0, int32_t hydrogen_shift=0, int32_t neighbor_shift=0,
                 int32_t distance_shift=0, int32_t disconnected=-1, int32_t unknown_h=0,
                 int32_t max_distance=0, int32_t max_neighbors=0, int32_t width=0, int32_t pad=0,
                 int32_t pad_diagonal=0, object vocabulary=None, int32_t unknown=-1):
        if unknown_h < 0 or unknown_h > ML_H_MAX:
            raise ValueError(f'unknown_h must be 0..{ML_H_MAX}, got {unknown_h}')
        if max_distance < 0:
            raise ValueError(f'max_distance must be non-negative, got {max_distance}')
        if max_neighbors < 0 or max_neighbors > ML_DEGREE_MAX:
            raise ValueError(f'max_neighbors must be 0..{ML_DEGREE_MAX}, got {max_neighbors}')
        if width < 0:
            raise ValueError(f'width must be non-negative, got {width}')

        self.element_shift = element_shift
        self.hydrogen_shift = hydrogen_shift
        self.neighbor_shift = neighbor_shift
        self.distance_shift = distance_shift
        self.disconnected = disconnected
        self.unknown_h = unknown_h
        self.max_distance = max_distance
        self.max_neighbors = max_neighbors
        self.width = width
        self.pad = pad
        self.pad_diagonal = pad_diagonal
        self.unknown = unknown
        self._table.keys = NULL
        self._table.values = NULL
        self._table.mask = 0
        self.vocabulary = None if vocabulary is None else dict(vocabulary)
        if self.vocabulary is not None:
            self._compile_vocabulary()

    cdef int _compile_vocabulary(self) except -1:
        """One open-addressed table, built once, read in the same C pass as the columns.

        The domains are checked here because a key that does not fit its field collides with another
        key, and a collision in a token table is a wrong training example with no symptom.
        """
        cdef uint32_t capacity = 8
        while capacity < 2 * <uint32_t> len(self.vocabulary):
            capacity <<= 1
        cdef uint32_t *keys = <uint32_t *> PyMem_Malloc(<size_t> capacity * sizeof(uint32_t))
        if keys is NULL:
            raise MemoryError('vocabulary table allocation failed')
        cdef int32_t *values = <int32_t *> PyMem_Malloc(<size_t> capacity * sizeof(int32_t))
        if values is NULL:
            PyMem_Free(keys)
            raise MemoryError('vocabulary table allocation failed')

        cdef uint32_t i, key, slot
        cdef object raw, token, element, hb, nb, ha, na
        for i in range(capacity):
            keys[i] = ML_KEY_EMPTY
        try:
            for raw, token in self.vocabulary.items():
                if len(raw) != 5:
                    raise ValueError('a vocabulary key is five numbers '
                                     '(element, h_before, n_before, h_after, n_after), '
                                     f'got {raw!r}')
                element, hb, nb, ha, na = raw
                if not 0 <= element <= ML_ELEMENT_MAX:
                    raise ValueError(f'vocabulary key element must be 0..{ML_ELEMENT_MAX}, '
                                     f'got {element} in {raw!r}')
                if not 0 <= hb <= ML_H_MAX or not 0 <= ha <= ML_H_MAX:
                    raise ValueError(f'vocabulary key hydrogen count must be 0..{ML_H_MAX}, '
                                     f'got {raw!r}')
                if not 0 <= nb <= ML_DEGREE_MAX or not 0 <= na <= ML_DEGREE_MAX:
                    raise ValueError(f'vocabulary key neighbor count must be 0..{ML_DEGREE_MAX}, '
                                     f'got {raw!r}')
                if not ML_TOKEN_MIN <= token <= ML_TOKEN_MAX:
                    raise ValueError(f'vocabulary token must fit int32, got {token} for {raw!r}')
                key = ML_VOCAB_KEY(element, hb, nb, ha, na)
                if key == ML_KEY_EMPTY:
                    raise ValueError('vocabulary key (0, 0, 0, 0, 0) packs to the empty-slot marker; '
                                     'the entry cannot be stored and that atom state always reads unknown')
                slot = key & (capacity - 1)
                while keys[slot] != ML_KEY_EMPTY:
                    if keys[slot] == key:
                        break
                    slot = (slot + 1) & (capacity - 1)
                keys[slot] = key
                values[slot] = <int32_t> token
        except BaseException:
            PyMem_Free(keys)
            PyMem_Free(values)
            raise
        self._table.keys = keys
        self._table.values = values
        self._table.mask = capacity - 1
        return 0

    def __dealloc__(self):
        PyMem_Free(self._table.keys)
        PyMem_Free(self._table.values)

    def __reduce__(self):
        return _rebuild_tensor_encoding, (self.element_shift, self.hydrogen_shift,
                                          self.neighbor_shift, self.distance_shift,
                                          self.disconnected, self.unknown_h, self.max_distance,
                                          self.max_neighbors, self.width, self.pad,
                                          self.pad_diagonal, self.vocabulary, self.unknown)

    def __repr__(self):
        """Only the knobs that were set: twelve zeros hide the one value that is not zero."""
        cdef object args, name, default, value
        args = []
        for name, default in (('element_shift', 0), ('hydrogen_shift', 0), ('neighbor_shift', 0),
                              ('distance_shift', 0), ('disconnected', -1), ('unknown_h', 0),
                              ('max_distance', 0), ('max_neighbors', 0), ('width', 0), ('pad', 0),
                              ('pad_diagonal', 0), ('unknown', -1)):
            value = getattr(self, name)
            if value != default:
                args.append(f'{name}={value}')
        if self.vocabulary is not None:
            args.append(f'vocabulary=<{len(self.vocabulary)} entries>')
        return f'TensorEncoding({", ".join(args)})'


def _rebuild_tensor_encoding(element_shift, hydrogen_shift, neighbor_shift, distance_shift,
                              disconnected, unknown_h, max_distance, max_neighbors, width, pad,
                              pad_diagonal, vocabulary, unknown):
    """The unpickler for `TensorEncoding`; a `cdef class` has no keyword `__init__` to call directly."""
    return TensorEncoding(element_shift, hydrogen_shift, neighbor_shift, distance_shift,
                          disconnected, unknown_h, max_distance, max_neighbors, width, pad,
                          pad_diagonal, vocabulary, unknown)


# --- view 1: atom state ---------------------------------------------------------------------------

cdef class StateView:
    """Per-atom columns and the pairwise distance block for one molecule.

    Five arrays and no methods: a framework wraps them, and anything it wants to do to them it does
    faster in its own tensor library.  `tokens` is None when the encoding carries no vocabulary.
    """
    cdef readonly object elements, hydrogens, neighbors, distances, tokens

    def __repr__(self):
        return f'StateView({self.elements.shape[0]} atoms' \
               f'{"" if self.tokens is None else ", tokens"})'


cdef inline int32_t _ml_clamp(int32_t value, int32_t limit, int32_t shift) noexcept nogil:
    """A clamp then a shift, in that order.  `limit == 0` is off."""
    if limit and value > limit:
        value = limit
    return value + shift


cdef inline int32_t _ml_token(const uint32_t *keys, const int32_t *values, uint32_t mask,
                              int32_t unknown, uint32_t key) noexcept nogil:
    """The token for a packed state key, or `unknown` on a miss.

    `keys`, `values` and `mask` come from a compiled `TensorEncoding._table`.  Linear probing over
    a table whose load factor is under 0.5, so the walk terminates on an empty slot; a full table
    cannot occur because the capacity is sized to at least twice the entry count.
    """
    cdef uint32_t slot = key & mask
    while keys[slot] != ML_KEY_EMPTY:
        if keys[slot] == key:
            return values[slot]
        slot = (slot + 1) & mask
    return unknown


def mol_state_view(MoleculeContainer molecule not None, TensorEncoding encoding=None):
    """Element, implicit hydrogen count, heavy degree and pairwise distance, as int32 arrays.

    THE CONTAINER IS THE INPUT AND THE PACKED BYTES ARE NOT.  A bytes-to-arrays kernel would need
    its own hydrogen and degree derivation, and one derivation shared by every reader is a standing
    rule; it would also forfeit H_UNKNOWN, R markers and every repair pass.

    `encoding=None` means physical values with an unstated hydrogen count reported as 0.
    """
    require_numpy()
    if encoding is None:
        encoding = TensorEncoding()

    cdef Structure structure = molecule._structure
    cdef uint32_t n_atoms = structure.header.atom_count
    cdef uint32_t width = n_atoms
    if encoding.width:
        if n_atoms > <uint32_t> encoding.width:
            raise ValueError(f'{n_atoms} atoms do not fit width={encoding.width}; '
                             'a truncated structure is a wrong training example')
        width = <uint32_t> encoding.width

    cdef object elements = _NP_EMPTY(width, dtype='int32')
    cdef object hydrogens = _NP_EMPTY(width, dtype='int32')
    cdef object neighbors = _NP_EMPTY(width, dtype='int32')
    cdef object distances = _NP_EMPTY((width, width), dtype='int32')
    cdef int32_t[::1] el_out = elements
    cdef int32_t[::1] h_out = hydrogens
    cdef int32_t[::1] n_out = neighbors
    cdef int32_t[:, ::1] d_out = distances
    cdef bint has_vocabulary = encoding.vocabulary is not None
    cdef object tokens = _NP_EMPTY(width, dtype='int32') if has_vocabulary else None
    cdef int32_t[::1] t_out = tokens if has_vocabulary else None

    cdef atom_t *atoms = structure.atoms()
    cdef uint32_t *ptr = csr_ptr(structure)
    cdef halfedge_t *edges = csr_edges(structure)
    cdef uint32_t n_half = ptr[n_atoms] if n_atoms else 0

    # one block, two regions (RULES.md §5.2): flattened half-edge targets, then the BFS queue.
    # `align8()` is not applied because both regions are `uint32_t`: the second region starts at
    # `block + n_half`, which is 4-byte aligned — exactly what a `uint32_t *` requires.  A block
    # holding mixed-width regions must use `align8()` between them.
    cdef uint32_t *block = NULL
    if n_atoms:
        block = <uint32_t *> PyMem_Malloc(<size_t> (n_half + n_atoms) * sizeof(uint32_t))
        if block is NULL:
            raise MemoryError('state view scratch allocation failed')
    cdef uint32_t *to = block
    cdef uint32_t *queue = block + n_half
    cdef uint32_t i, j, k, vocab_key
    cdef int32_t h, d, vocab_key_h, vocab_unknown
    cdef int32_t pad = encoding.pad
    cdef int32_t pad_diagonal = encoding.pad_diagonal
    cdef int32_t unknown_h = encoding.unknown_h
    cdef int32_t disconnected = encoding.disconnected
    cdef int32_t distance_shift = encoding.distance_shift
    cdef int32_t max_distance = encoding.max_distance
    cdef int32_t element_shift = encoding.element_shift
    cdef int32_t hydrogen_shift = encoding.hydrogen_shift
    cdef int32_t neighbor_shift = encoding.neighbor_shift
    cdef int32_t max_neighbors = encoding.max_neighbors
    cdef uint32_t *vocab_keys = encoding._table.keys
    cdef int32_t *vocab_values = encoding._table.values
    cdef uint32_t vocab_mask = encoding._table.mask
    vocab_unknown = encoding.unknown
    cdef atom_t *a
    try:
        with nogil:
            for i in range(n_atoms):
                # the hydrogen nibble is read once and branched on; H_UNKNOWN is 15 and is a
                # sentinel rather than a count, so it never reaches the shift as a number.
                a = &atoms[i]
                h = a.hydrogens & 0x0f
                el_out[i] = <int32_t> a.element + element_shift
                h_out[i] = (unknown_h if h == H_UNKNOWN else h) + hydrogen_shift
                n_out[i] = _ml_clamp(<int32_t> a.degree, max_neighbors, neighbor_shift)
                if has_vocabulary:
                    # key from physical unclamped values, unknown_h applied; clamps and shifts come
                    # after (ORDER OF OPERATIONS in TensorEncoding docstring).
                    vocab_key_h = unknown_h if h == H_UNKNOWN else h
                    vocab_key = ML_VOCAB_KEY(<uint32_t> a.element, <uint32_t> vocab_key_h,
                                            <uint32_t> a.degree, <uint32_t> vocab_key_h,
                                            <uint32_t> a.degree)
                    t_out[i] = _ml_token(vocab_keys, vocab_values, vocab_mask,
                                         vocab_unknown, vocab_key)
            for i in range(n_atoms, width):
                el_out[i] = pad
                h_out[i] = pad
                n_out[i] = pad
                if has_vocabulary:
                    t_out[i] = pad

            if n_atoms:
                for k in range(n_half):
                    to[k] = edges[k].to
                csr_bfs_all(ptr, to, n_atoms, &d_out[0, 0], width, queue)
                # in place over the block the BFS filled: a physical -1 becomes `disconnected`
                # verbatim, and every real distance is clamped then shifted.  A shifted sentinel
                # collides with a real distance, so the two paths never meet.
                for i in range(n_atoms):
                    for j in range(n_atoms):
                        d = d_out[i, j]
                        if d < 0:
                            d_out[i, j] = disconnected
                        else:
                            d_out[i, j] = _ml_clamp(d, max_distance, distance_shift)
                    for j in range(n_atoms, width):
                        d_out[i, j] = pad
            for i in range(n_atoms, width):
                for j in range(width):
                    d_out[i, j] = pad
                if pad_diagonal:
                    d_out[i, i] = pad_diagonal
    finally:
        PyMem_Free(block)

    cdef StateView view = StateView.__new__(StateView)
    view.elements = elements
    view.hydrogens = hydrogens
    view.neighbors = neighbors
    view.distances = distances
    view.tokens = tokens
    return view


# --- view 2: transition ---------------------------------------------------------------------------

cdef class TransitionView:
    """Per-atom state on each side of a transformation, over the union graph.

    A MOLECULE IS THE `before == after` CASE, so one layout and one vocabulary serve molecules and
    reactions.  `bond_before` / `bond_after` are the reaction centre: (1, 0) a broken bond, (0, 1) a
    formed one, (1, 2) an order change; 0 is a bond absent on that side, which is the whole of what a
    dynamic bond carried.

    `unmapped` and `collisions` REPORT AND DO NOT REFUSE.  A badly mapped record still produces a
    view, and the pipeline boundary decides what to do about it.
    """
    cdef readonly object elements, h_before, n_before, h_after, n_after, map_numbers
    cdef readonly object distances, bonds, bond_before, bond_after, tokens
    cdef readonly object unmapped, collisions

    def __repr__(self):
        return f'TransitionView({self.elements.shape[0]} atoms, {self.bonds.shape[0]} bonds' \
               f'{"" if self.tokens is None else ", tokens"})'


cdef struct ml_union_t:
    # One side-resolved union graph, struct-of-arrays, owned by the caller's single block.
    # `n_atoms` union atoms and `n_bonds` union bonds.  `bond_i[k] < bond_j[k]`;
    # `order_before` / `order_after` are the per-side orders, 0 where the bond is absent.
    uint32_t n_atoms
    uint32_t n_bonds
    uint8_t *element
    uint8_t *h_before
    uint8_t *n_before
    uint8_t *h_after
    uint8_t *n_after
    uint16_t *map_number
    uint32_t *bond_i
    uint32_t *bond_j
    uint8_t *order_before
    uint8_t *order_after


cdef object _ml_fill_transition_arrays(ml_union_t *union_graph, TensorEncoding encoding,
                                       object unmapped, object collisions):
    """Apply the encoding to a built union and return the `TransitionView`.

    THE ONE PLACE THE ENCODING IS APPLIED to a transition.  Molecules and reactions differ in how
    the union is built and in nothing after it, so a second copy of the shift-and-clamp arithmetic
    here would be a second convention nobody declared.

    Allocates the numpy arrays and the BFS scratch; the union's own memory belongs to the caller.
    """
    cdef uint32_t n_atoms = union_graph.n_atoms
    cdef uint32_t n_bonds = union_graph.n_bonds
    cdef uint32_t width = n_atoms
    if encoding.width:
        if n_atoms > <uint32_t> encoding.width:
            raise ValueError(f'{n_atoms} atoms do not fit width={encoding.width}; '
                             'a truncated structure is a wrong training example')
        width = <uint32_t> encoding.width

    cdef object elements = _NP_EMPTY(width, dtype='int32')
    cdef object h_before = _NP_EMPTY(width, dtype='int32')
    cdef object n_before = _NP_EMPTY(width, dtype='int32')
    cdef object h_after = _NP_EMPTY(width, dtype='int32')
    cdef object n_after = _NP_EMPTY(width, dtype='int32')
    cdef object map_numbers = _NP_EMPTY(width, dtype='int32')
    cdef object distances = _NP_EMPTY((width, width), dtype='int32')
    cdef object bonds = _NP_EMPTY((n_bonds, 2), dtype='int32')
    cdef object bond_before = _NP_EMPTY(n_bonds, dtype='int32')
    cdef object bond_after = _NP_EMPTY(n_bonds, dtype='int32')
    cdef bint has_vocabulary = encoding.vocabulary is not None
    cdef object tokens = _NP_EMPTY(width, dtype='int32') if has_vocabulary else None

    cdef int32_t[::1] el_out = elements
    cdef int32_t[::1] hb_out = h_before
    cdef int32_t[::1] nb_out = n_before
    cdef int32_t[::1] ha_out = h_after
    cdef int32_t[::1] na_out = n_after
    cdef int32_t[::1] mn_out = map_numbers
    cdef int32_t[:, ::1] d_out = distances
    cdef int32_t[:, ::1] b_out = bonds
    cdef int32_t[::1] ob_out = bond_before
    cdef int32_t[::1] oa_out = bond_after
    cdef int32_t[::1] t_out = tokens if has_vocabulary else None

    # hoist every encoding field used in the nogil kernel (one convention, one copy — see `mol_state_view`)
    cdef int32_t element_shift = encoding.element_shift
    cdef int32_t hydrogen_shift = encoding.hydrogen_shift
    cdef int32_t neighbor_shift = encoding.neighbor_shift
    cdef int32_t distance_shift = encoding.distance_shift
    cdef int32_t disconnected = encoding.disconnected
    cdef int32_t unknown_h = encoding.unknown_h
    cdef int32_t max_distance = encoding.max_distance
    cdef int32_t max_neighbors = encoding.max_neighbors
    cdef int32_t pad = encoding.pad
    cdef int32_t pad_diagonal = encoding.pad_diagonal
    cdef uint32_t *vocab_keys = encoding._table.keys
    cdef int32_t *vocab_values = encoding._table.values
    cdef uint32_t vocab_mask = encoding._table.mask
    cdef int32_t vocab_unknown = encoding.unknown

    # one block, three regions (RULES.md §5.2): union CSR built by counting sort, its flattened
    # targets, and the BFS queue.  Every bond contributes two half-edges; `n_atoms + 1` is the ptr
    # array.  All three regions are uint32_t, so no align8() between them is needed.
    cdef size_t n_half = <size_t> 2 * n_bonds
    cdef uint32_t *block = <uint32_t *> PyMem_Malloc(
        (<size_t> n_atoms + 1 + n_half + n_atoms) * sizeof(uint32_t))
    if block is NULL:
        raise MemoryError('transition view scratch allocation failed')
    cdef uint32_t *csr_p = block
    cdef uint32_t *csr_t = block + n_atoms + 1
    cdef uint32_t *queue = csr_t + n_half
    cdef uint32_t i, j, k, n, m
    cdef int32_t d, hb, ha
    try:
        with nogil:
            for i in range(n_atoms):
                # `unknown_h` is applied before the shift and before the key, exactly as in
                # `mol_state_view`; H_UNKNOWN is a sentinel and never reaches arithmetic as a count.
                hb = unknown_h if union_graph.h_before[i] == H_UNKNOWN \
                    else <int32_t> union_graph.h_before[i]
                ha = unknown_h if union_graph.h_after[i] == H_UNKNOWN \
                    else <int32_t> union_graph.h_after[i]
                el_out[i] = <int32_t> union_graph.element[i] + element_shift
                hb_out[i] = hb + hydrogen_shift
                ha_out[i] = ha + hydrogen_shift
                nb_out[i] = _ml_clamp(<int32_t> union_graph.n_before[i], max_neighbors,
                                      neighbor_shift)
                na_out[i] = _ml_clamp(<int32_t> union_graph.n_after[i], max_neighbors,
                                      neighbor_shift)
                mn_out[i] = <int32_t> union_graph.map_number[i]
                if has_vocabulary:
                    # key from physical unclamped values with unknown_h applied; clamps and shifts
                    # come after (ORDER OF OPERATIONS in TensorEncoding docstring).
                    t_out[i] = _ml_token(vocab_keys, vocab_values, vocab_mask, vocab_unknown,
                                         ML_VOCAB_KEY(union_graph.element[i], <uint32_t> hb,
                                                      union_graph.n_before[i], <uint32_t> ha,
                                                      union_graph.n_after[i]))
            for i in range(n_atoms, width):
                el_out[i] = pad
                hb_out[i] = pad
                ha_out[i] = pad
                nb_out[i] = pad
                na_out[i] = pad
                mn_out[i] = pad
                if has_vocabulary:
                    t_out[i] = pad

            # bond orders are chemical values and are never shifted: 0 means absent on that side,
            # and a shift would put a real order where the absence marker is.
            for k in range(n_bonds):
                b_out[k, 0] = <int32_t> union_graph.bond_i[k]
                b_out[k, 1] = <int32_t> union_graph.bond_j[k]
                ob_out[k] = <int32_t> union_graph.order_before[k]
                oa_out[k] = <int32_t> union_graph.order_after[k]

            # union CSR by counting sort, O(V+E): clear, degree tally, prefix sum, then place.
            # `csr_p` doubles as a running cursor during the place pass and is restored by the
            # shift-down that follows.  Without the restore, `csr_p[i]` would be `csr_p[i+1]` and
            # `csr_bfs_all` would walk every edge one position off, producing a plausible-looking
            # but wrong distance matrix with no test to catch it.
            for i in range(n_atoms + 1):
                csr_p[i] = 0
            for k in range(n_bonds):
                csr_p[union_graph.bond_i[k] + 1] += 1
                csr_p[union_graph.bond_j[k] + 1] += 1
            for i in range(n_atoms):
                csr_p[i + 1] += csr_p[i]
            for k in range(n_bonds):
                n = union_graph.bond_i[k]
                m = union_graph.bond_j[k]
                csr_t[csr_p[n]] = m
                csr_p[n] += 1
                csr_t[csr_p[m]] = n
                csr_p[m] += 1
            for i in range(n_atoms, 0, -1):
                csr_p[i] = csr_p[i - 1]
            csr_p[0] = 0

            if n_atoms:
                csr_bfs_all(csr_p, csr_t, n_atoms, &d_out[0, 0], width, queue)
                for i in range(n_atoms):
                    for j in range(n_atoms):
                        d = d_out[i, j]
                        if d < 0:
                            d_out[i, j] = disconnected
                        else:
                            d_out[i, j] = _ml_clamp(d, max_distance, distance_shift)
                    for j in range(n_atoms, width):
                        d_out[i, j] = pad
            for i in range(n_atoms, width):
                for j in range(width):
                    d_out[i, j] = pad
                if pad_diagonal:
                    d_out[i, i] = pad_diagonal
    finally:
        PyMem_Free(block)

    cdef TransitionView tv = TransitionView.__new__(TransitionView)
    tv.elements = elements
    tv.h_before = h_before
    tv.n_before = n_before
    tv.h_after = h_after
    tv.n_after = n_after
    tv.map_numbers = map_numbers
    tv.distances = distances
    tv.bonds = bonds
    tv.bond_before = bond_before
    tv.bond_after = bond_after
    tv.tokens = tokens
    tv.unmapped = unmapped
    tv.collisions = collisions
    return tv


def mol_transition_view(MoleculeContainer molecule not None, TensorEncoding encoding=None):
    """`transition_view` for a molecule: the same state on both sides, over its own graph."""
    require_numpy()
    if encoding is None:
        encoding = TensorEncoding()

    cdef Structure structure = molecule._structure
    cdef uint32_t n_atoms = structure.header.atom_count
    cdef atom_t *atoms = structure.atoms()
    cdef uint32_t *mol_ptr = csr_ptr(structure)
    cdef halfedge_t *edges = csr_edges(structure)
    cdef uint32_t n_half = mol_ptr[n_atoms] if n_atoms else 0
    cdef uint32_t n_bonds = n_half // 2

    # one block, ten regions (RULES.md §5.2).  The four leading uint8 columns share no alignment
    # padding between them (alignment 1 each); the fifth uint8 column is padded to align the
    # uint16_t map_number region that follows.  bond_i and bond_j are packed into one aligned
    # block (2*n_bonds uint32_t = 8*n_bonds bytes, always 8-byte aligned); order_before and
    # order_after are packed into one aligned block (2*n_bonds uint8_t, padded to 8 bytes).
    cdef size_t need = (<size_t> 4 * n_atoms
                        + align8(<size_t> n_atoms)                       # n_after, padded to align map_number
                        + align8(<size_t> n_atoms * sizeof(uint16_t))    # map_number
                        + align8(<size_t> 2 * n_bonds * sizeof(uint32_t))  # bond_i + bond_j
                        + align8(<size_t> 2 * n_bonds))                  # order_before + order_after
    cdef char *block = <char *> PyMem_Malloc(need if need else 8)
    if block is NULL:
        raise MemoryError('transition union allocation failed')
    cdef ml_union_t ug
    cdef char *cursor = block
    ug.n_atoms = n_atoms
    ug.n_bonds = n_bonds
    ug.element = <uint8_t *> cursor; cursor += n_atoms
    ug.h_before = <uint8_t *> cursor; cursor += n_atoms
    ug.n_before = <uint8_t *> cursor; cursor += n_atoms
    ug.h_after = <uint8_t *> cursor; cursor += n_atoms
    ug.n_after = <uint8_t *> cursor; cursor += align8(<size_t> n_atoms)
    ug.map_number = <uint16_t *> cursor
    cursor += align8(<size_t> n_atoms * sizeof(uint16_t))
    ug.bond_i = <uint32_t *> cursor                   # bond_i and bond_j share one aligned block
    ug.bond_j = (<uint32_t *> cursor) + n_bonds
    cursor += align8(<size_t> 2 * n_bonds * sizeof(uint32_t))
    ug.order_before = <uint8_t *> cursor              # order_before and order_after share one block
    ug.order_after = <uint8_t *> cursor + n_bonds

    cdef uint32_t i, k, w, fill = 0
    cdef uint8_t h, eord
    cdef atom_t *a
    try:
        with nogil:
            for i in range(n_atoms):
                a = &atoms[i]
                h = a.hydrogens & 0x0f
                ug.element[i] = a.element
                ug.h_before[i] = h
                ug.h_after[i] = h
                ug.n_before[i] = a.degree
                ug.n_after[i] = a.degree
                ug.map_number[i] = a.map_number
            # each undirected bond once, taken from the half-edge whose source index is smaller;
            # this normalises the pair to `i < j` without a sort.
            for i in range(n_atoms):
                for k in range(mol_ptr[i], mol_ptr[i + 1]):
                    w = edges[k].to
                    eord = edges[k].order
                    if i < w:
                        ug.bond_i[fill] = i
                        ug.bond_j[fill] = w
                        ug.order_before[fill] = eord
                        ug.order_after[fill] = eord
                        fill += 1
            # `n_bonds` was sized from `n_half // 2`; write the count actually filled in case a
            # stored graph has a duplicate half-edge (self-loop or arena error).
            ug.n_bonds = fill
        return _ml_fill_transition_arrays(&ug, encoding,
                                          {'reactants': 0, 'products': 0},
                                          {'reactants': (), 'products': ()})
    finally:
        PyMem_Free(block)


# --- view 3: reaction transition ------------------------------------------------------------------

DEF ML_NO_INDEX = 0xffffffff


# `uidx` is `map_number -> union index`, valid only for touched slots, and it is what PAIRS the two
# sides -- an unmapped atom therefore takes a union row without an entry in it.  `rseen` / `pseen` mark
# which side already claimed a map number.  `urow` is indexed by running ordinal over all atoms of both
# sides (left first, then right) and holds the union row the atom took, or `ML_NO_INDEX` for a rejected
# collision; the bond passes read it for both endpoints, which is why they need no map number of their
# own.  `head` / `link` are a per-atom singly linked list over bond slots, keyed on the lower union index
# of the pair.  `side` codes: 1 = reactant-only, 2 = product-only, 3 = both.  The uidx table is max map
# number + 1, at most 10000 slots -- initialized slot by slot, not memset, because a three-atom record
# must not pay for all.
cdef struct ml_merge_t:
    uint32_t *uidx
    uint8_t *rseen
    uint8_t *pseen
    uint8_t *side              # 1 = reactant-only, 2 = product-only, 3 = both
    uint32_t *urow             # union row per atom ordinal; ML_NO_INDEX = a rejected collision
    uint32_t *head
    uint32_t *link


cdef inline void _ml_union_row(ml_union_t *g, ml_merge_t *merge, uint32_t row, atom_t *at,
                               uint32_t mn, uint8_t side) noexcept nogil:
    """One union row from one atom, on the side that atom came from.

    `h_after` is written equal to `h_before`, which a paired product atom then overwrites.  An atom on
    one side only keeps it -- the convention `ReactionModelingView.states` states -- and an unmapped
    atom is one of those, the record having said nothing that would pair it.
    """
    merge.side[row] = side
    merge.head[row] = ML_NO_INDEX
    g.element[row] = at.element
    g.h_before[row] = at.hydrogens & 0x0f
    g.h_after[row] = at.hydrogens & 0x0f
    g.map_number[row] = <uint16_t> mn


cdef void _ml_union_degrees(ml_union_t *union_graph, ml_merge_t *merge) noexcept nogil:
    """Per-side heavy degree over the union, by the rule `ReactionModelingView.states` states.

    An atom present on both sides counts, on each side, the bonds that side has.  A REACTANT-ONLY ATOM
    counts, after, only its reactant-only neighbours: it left with the fragment it belongs to, and the
    bond to the atom that stayed is not a bond it still has.  A product-only atom is the mirror.

    One pass over the union bonds, both endpoints of each, so no scratch and no second pass.
    """
    cdef uint32_t k, n, m, n_atoms = union_graph.n_atoms
    for n in range(n_atoms):
        union_graph.n_before[n] = 0
        union_graph.n_after[n] = 0
    for k in range(union_graph.n_bonds):
        n = union_graph.bond_i[k]
        m = union_graph.bond_j[k]
        if union_graph.order_before[k]:
            # `side == 1` is reactant-only: after the reaction it keeps only the neighbours that
            # left with it
            if merge.side[n] != 2:
                union_graph.n_before[n] += 1
                if merge.side[n] == 1 and merge.side[m] == 1:
                    union_graph.n_after[n] += 1
            if merge.side[m] != 2:
                union_graph.n_before[m] += 1
                if merge.side[n] == 1 and merge.side[m] == 1:
                    union_graph.n_after[m] += 1
        if union_graph.order_after[k]:
            if merge.side[n] != 1:
                union_graph.n_after[n] += 1
                if merge.side[n] == 2 and merge.side[m] == 2:
                    union_graph.n_before[n] += 1
            if merge.side[m] != 1:
                union_graph.n_after[m] += 1
                if merge.side[n] == 2 and merge.side[m] == 2:
                    union_graph.n_before[m] += 1


def reaction_transition_view(reactants, products, TensorEncoding encoding=None):
    """The union of a mapped reaction's two sides as int32 arrays.

    AGENTS ARE EXCLUDED BY NOT BEING PASSED.  The caller hands over the two sides it wants unioned, so
    there is no side-kind test inside the kernel to keep in step with the container's notion of one.

    AN UNMAPPED ATOM IS PLACED ON THE SIDE IT CAME FROM, as reactant-only or product-only, and counted
    in `unmapped`.  Left out, it takes the degree of every neighbour that stayed down with it: an aryl
    bromide whose bromine carries no number reads two heavy neighbours before the coupling and three
    after, so every aryl halide of one ring gives the same transition state.  The price is that an atom
    the record leaves bare on BOTH sides is two rows rather than one, nothing in the record pairing
    them, and `unmapped` is what a consumer that cannot accept that reads.

    A BADLY MAPPED RECORD IS REPORTED.  A map number claimed twice on one side is listed in `collisions`
    and the first claim keeps the row; the second contributes nothing, since a union cannot hold two
    atoms at one key.  Refusals live at the answer boundary and this is not one.
    """
    require_numpy()
    if encoding is None:
        encoding = TensorEncoding()

    cdef list left = list(reactants)
    cdef list right = list(products)
    cdef MoleculeContainer mol
    for mol in left + right:
        mol._require_clean()

    cdef Structure structure
    cdef atom_t *atoms
    cdef uint32_t *ptr
    cdef halfedge_t *edges

    # capacities: every atom on both sides could be mapped and distinct, and every bond could be one.
    # `n_mol` and `i` are hoisted before the capacity loop so both `for i in range(n_mol)` uses
    # below are typed rather than implicitly declared Python objects.
    cdef uint32_t cap_atoms = 0, cap_bonds = 0, max_mn = 0
    cdef uint32_t n_mol, i
    cdef atom_t *at
    for mol in left + right:
        structure = mol._structure
        atoms = structure.atoms()
        ptr = csr_ptr(structure)
        n_mol = structure.header.atom_count
        cap_atoms += n_mol
        cap_bonds += ptr[n_mol] // 2 if n_mol else 0
        for i in range(n_mol):
            if atoms[i].map_number > max_mn:
                max_mn = atoms[i].map_number

    # `slots` = max map number + 1, at most 10000: the touch pass initializes only touched slots,
    # so a record with one map number of 9999 costs two uint32_t writes, not 10000 memset bytes.
    cdef uint32_t slots = max_mn + 1
    cdef size_t need = (align8(<size_t> slots * sizeof(uint32_t))            # uidx
                        + align8(<size_t> slots)                             # rseen
                        + align8(<size_t> slots)                             # pseen
                        + align8(<size_t> cap_atoms)                         # side
                        + align8(<size_t> cap_atoms * sizeof(uint32_t))      # urow
                        + align8(<size_t> cap_atoms * sizeof(uint32_t))      # head
                        + align8(<size_t> cap_bonds * sizeof(uint32_t))      # link
                        + <size_t> 4 * cap_atoms                             # element + h_before + n_before + h_after
                        + align8(<size_t> cap_atoms)                         # n_after, padded to align map_number
                        + align8(<size_t> cap_atoms * sizeof(uint16_t))      # map_number
                        + 2 * align8(<size_t> cap_bonds * sizeof(uint32_t))  # bond_i + bond_j
                        + 2 * align8(<size_t> cap_bonds))                    # order_before + order_after
    cdef char *block = <char *> PyMem_Malloc(need if need else 8)
    if block is NULL:
        raise MemoryError('reaction union allocation failed')

    cdef ml_merge_t merge
    cdef ml_union_t union_graph
    cdef char *cursor = block
    merge.uidx = <uint32_t *> cursor; cursor += align8(<size_t> slots * sizeof(uint32_t))
    merge.rseen = <uint8_t *> cursor; cursor += align8(<size_t> slots)
    merge.pseen = <uint8_t *> cursor; cursor += align8(<size_t> slots)
    merge.side = <uint8_t *> cursor; cursor += align8(<size_t> cap_atoms)
    merge.urow = <uint32_t *> cursor; cursor += align8(<size_t> cap_atoms * sizeof(uint32_t))
    merge.head = <uint32_t *> cursor; cursor += align8(<size_t> cap_atoms * sizeof(uint32_t))
    merge.link = <uint32_t *> cursor; cursor += align8(<size_t> cap_bonds * sizeof(uint32_t))
    union_graph.element = <uint8_t *> cursor; cursor += cap_atoms
    union_graph.h_before = <uint8_t *> cursor; cursor += cap_atoms
    union_graph.n_before = <uint8_t *> cursor; cursor += cap_atoms
    union_graph.h_after = <uint8_t *> cursor; cursor += cap_atoms
    union_graph.n_after = <uint8_t *> cursor; cursor += align8(<size_t> cap_atoms)
    union_graph.map_number = <uint16_t *> cursor
    cursor += align8(<size_t> cap_atoms * sizeof(uint16_t))
    union_graph.bond_i = <uint32_t *> cursor
    cursor += align8(<size_t> cap_bonds * sizeof(uint32_t))
    union_graph.bond_j = <uint32_t *> cursor
    cursor += align8(<size_t> cap_bonds * sizeof(uint32_t))
    union_graph.order_before = <uint8_t *> cursor; cursor += align8(<size_t> cap_bonds)
    union_graph.order_after = <uint8_t *> cursor

    cdef uint32_t n_union = 0, n_bonds = 0
    cdef uint32_t unmapped_r = 0, unmapped_p = 0
    cdef list collisions_r = [], collisions_p = []
    cdef uint32_t k, mn, n, m, e, base, right_base
    cdef uint32_t edge_to
    cdef uint8_t order
    try:
        # touch pass: initialize exactly the slots this record will read.  A blanket memset of up
        # to 10000 slots on a three-atom record would cost more than the whole view.
        for mol in left + right:
            structure = mol._structure
            atoms = structure.atoms()
            n_mol = structure.header.atom_count
            with nogil:
                for i in range(n_mol):
                    mn = atoms[i].map_number
                    if mn:
                        merge.uidx[mn] = ML_NO_INDEX
                        merge.rseen[mn] = 0
                        merge.pseen[mn] = 0

        # reactant atoms, in container order, take the low union indices
        base = 0
        for mol in left:
            structure = mol._structure
            atoms = structure.atoms()
            n_mol = structure.header.atom_count
            for i in range(n_mol):
                at = atoms + i
                mn = at.map_number
                if mn:
                    if merge.rseen[mn]:
                        collisions_r.append(mn)
                        merge.urow[base + i] = ML_NO_INDEX
                        continue
                    merge.rseen[mn] = 1
                    merge.uidx[mn] = n_union
                else:
                    unmapped_r += 1
                _ml_union_row(&union_graph, &merge, n_union, at, mn, 1)
                merge.urow[base + i] = n_union
                n_union += 1
            base += n_mol
        right_base = base  # starting ordinal for product atoms in the urow array

        # product atoms: a known map number reads its after state onto the existing row, a new one
        # appends -- which is what keeps a component contiguous in the union order
        for mol in right:
            structure = mol._structure
            atoms = structure.atoms()
            n_mol = structure.header.atom_count
            for i in range(n_mol):
                at = atoms + i
                mn = at.map_number
                if mn:
                    if merge.pseen[mn]:
                        collisions_p.append(mn)
                        merge.urow[base + i] = ML_NO_INDEX
                        continue
                    merge.pseen[mn] = 1
                    if merge.rseen[mn]:
                        n = merge.uidx[mn]
                        merge.side[n] = 3
                        union_graph.h_after[n] = at.hydrogens & 0x0f
                        merge.urow[base + i] = n
                        continue
                    merge.uidx[mn] = n_union
                else:
                    unmapped_p += 1
                _ml_union_row(&union_graph, &merge, n_union, at, mn, 2)
                merge.urow[base + i] = n_union
                n_union += 1
            base += n_mol

        # reactant bonds, both endpoints holding a union row: appended in reactant order, `n < m`, so
        # the `n >= m` guard is what visits each bond once.
        base = 0
        for mol in left:
            structure = mol._structure
            ptr = csr_ptr(structure)
            edges = csr_edges(structure)
            n_mol = structure.header.atom_count
            for i in range(n_mol):
                n = merge.urow[base + i]
                if n == ML_NO_INDEX:
                    continue
                for k in range(ptr[i], ptr[i + 1]):
                    edge_to = edges[k].to
                    m = merge.urow[base + edge_to]
                    if m == ML_NO_INDEX or n >= m:
                        continue
                    union_graph.bond_i[n_bonds] = n
                    union_graph.bond_j[n_bonds] = m
                    union_graph.order_before[n_bonds] = edges[k].order
                    union_graph.order_after[n_bonds] = 0
                    merge.link[n_bonds] = merge.head[n]
                    merge.head[n] = n_bonds
                    n_bonds += 1
            base += n_mol

        # product bonds: merged onto a reactant bond when the pair already exists, appended otherwise
        base = right_base
        for mol in right:
            structure = mol._structure
            ptr = csr_ptr(structure)
            edges = csr_edges(structure)
            n_mol = structure.header.atom_count
            for i in range(n_mol):
                n = merge.urow[base + i]
                if n == ML_NO_INDEX:
                    continue
                for k in range(ptr[i], ptr[i + 1]):
                    edge_to = edges[k].to
                    m = merge.urow[base + edge_to]
                    if m == ML_NO_INDEX or n >= m:
                        continue
                    order = edges[k].order
                    e = merge.head[n]
                    while e != ML_NO_INDEX:
                        if union_graph.bond_j[e] == m:
                            break
                        e = merge.link[e]
                    if e != ML_NO_INDEX:
                        union_graph.order_after[e] = order
                        continue
                    union_graph.bond_i[n_bonds] = n
                    union_graph.bond_j[n_bonds] = m
                    union_graph.order_before[n_bonds] = 0
                    union_graph.order_after[n_bonds] = order
                    merge.link[n_bonds] = merge.head[n]
                    merge.head[n] = n_bonds
                    n_bonds += 1
            base += n_mol

        union_graph.n_atoms = n_union
        union_graph.n_bonds = n_bonds
        with nogil:
            _ml_union_degrees(&union_graph, &merge)
        return _ml_fill_transition_arrays(
            &union_graph, encoding,
            {'reactants': unmapped_r, 'products': unmapped_p},
            {'reactants': tuple(sorted(set(collisions_r))),
             'products': tuple(sorted(set(collisions_p)))})
    finally:
        PyMem_Free(block)
