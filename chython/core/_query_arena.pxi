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
# Query arena layout.  A query atom's constraints are expressed as a sequence of box records
# (qbox_t).  Each box holds four "forbidden" 64-bit masks: if the candidate atom's feature
# word AND the corresponding neg word is non-zero, the candidate is rejected.  The neg
# (forbidden) encoding, rather than an allowed-set encoding, makes the test exact for
# one-hot spans: since exactly one bit per span fires per atom, a forbidden-bit test and a
# not-in-allowed-set test are equivalent.  One-hot spans cover element, charge-range,
# hybridization and similar descriptors; multi-bit spans (ring-size bitmaps) use separate
# qany_t records for OR conditions.


DEF QUERY_MAGIC = 0x43485951        # 'CHYQ'
DEF QUERY_VERSION = 1


cdef enum:
    QSEG_ATOMS = 0
    QSEG_BOXES = 1
    QSEG_ANY = 2
    QSEG_BONDS = 3
    QSEG_BOND_BOXES = 4
    QSEG_CLOSURES = 5
    QSEG_COMPONENTS = 6
    QSEG_ELEMENT_DEMAND = 7
    # QSEG_DEMAND_LIST holds the demanded-element pairs used by query_may_match: one
    # (element, count) uint32_t pair per demanded element.  Built at seal from the
    # QSEG_ELEMENT_DEMAND histogram; query_may_match iterates this list (O(entries)) rather than
    # scanning the full 118-slot histogram (O(118)) -- the compact form is both faster and more
    # expressive of intent.  A zero-length list means no element is demanded.
    QSEG_DEMAND_LIST = 8
    # QSEG_AUTOMORPHISM holds header.automorphism_count rows of atom_count uint32_t: one
    # permutation of DFS positions per row, the identity excluded.  query_alloc sizes it exactly,
    # like every other segment -- query_seal computes the group before it emits, so nothing is ever
    # appended to a sealed arena and no consumer needs a count-based guard: a zero-length segment
    # routes to the zero page and header.automorphism_count is 0 in exactly that case, so every
    # loop over it is empty.  It stays after the count-derived segments because its size depends on
    # the query's symmetry rather than on its atom and bond counts.
    QSEG_AUTOMORPHISM = 9
    # QSEG_STEREO holds header.stereo_count qstereo_t records: one per query atom carrying a stereo
    # primitive, and one per geometry a `/` and `\` pair states -- the frame the statement is about,
    # plus the DFS position at which that frame is complete.  Sized by the query's stereo content, so
    # it sits with QSEG_AUTOMORPHISM after the count-derived segments.  Zero-length when the query
    # states no configuration at all, which is also exactly when QFLAG_HAS_STEREO is clear.
    QSEG_STEREO = 10
    QSEG_COUNT = 11


cdef enum:
    QFLAG_HAS_GROUP = 1
    QFLAG_HAS_STEREO = 2
    QFLAG_HAS_MASKED = 4
    QFLAG_ASYMMETRIC = 8                # Task 12: the automorphism group is trivial
    QFLAG_PARTIAL_AUTOMORPHISM = 16     # Task 12: the group search hit its node or its row cap


cdef enum:
    QATOM_MASKED = 1        # qatom_t.flags bit 0
    QATOM_STEREO_CW = 2     # bit 1: some box of this atom demands '@'
    QATOM_STEREO_CCW = 4    # bit 2: some box of this atom demands '@@'.  Both bits set means the
                            # atom's boxes disagree, or one box ANDed the two -- which box applies
                            # is a per-candidate question the kernel answers (ruling F87)
    QATOM_ROOT = 8          # bit 3: a component root -- seed from element buckets
    QATOM_ANY_ELEMENT = 16  # bit 4: some box of this atom leaves the element span untouched, so the
                            # atom constrains no element at all -- `[A]`, `[*]`, or a bracket that
                            # only counts, like `[D2]`
    QATOM_METAL_ELEMENT = 32  # bit 5: EVERY box constrains the element to exactly the 93 metals --
                            # `[M]`.  Exclusive with QATOM_ANY_ELEMENT by construction: an untouched
                            # span is not the metal mask.
                            #
                            # BOTH ARE DERIVED AT SEAL FROM THE COMPILED TERM, not journalled the way
                            # QATOM_MASKED is.  A mask is a caller's declaration and nothing else can
                            # tell you about it; "does this atom name an element" is a question the
                            # boxes already answer, so a second, writable copy of the answer could
                            # only drift from them.  They exist because a rule table needs the
                            # wildcard-versus-named distinction: a wildcard is shared CONTEXT, so two
                            # matches of one rule may overlap there, while a named atom identifies the
                            # site being repaired.


cdef packed struct qatom_t:      # 24 bytes; spare is reserved headroom for the stereo epic
    uint32_t box_begin
    uint16_t box_count
    uint16_t flags
    uint32_t back
    uint32_t closure_begin
    uint16_t closure_count
    uint16_t map_number
    uint32_t spare


cdef packed struct qbox_t:       # 40 bytes
    uint64_t neg[4]
    uint32_t any_begin
    uint16_t any_count
    uint8_t sign                 # QSIGN_* bitmask, ruling F87; 0 on every bond box
    uint8_t spare


cdef packed struct qany_t:       # 16 bytes
    uint64_t mask
    uint32_t word
    uint32_t spare


cdef packed struct qbond_t:      # 8 bytes
    uint32_t box_begin
    uint16_t box_count
    uint16_t flags


cdef packed struct qclosure_t:   # 8 bytes
    uint32_t to_index
    uint32_t bond_index


cdef packed struct qstereo_t:    # 28 bytes: one stereo statement a query makes
    # TWO KINDS OF RECORD SHARE THE FIELDS, and `spare`'s low byte says which -- SU_TETRA for a
    # centre's sign and SU_CIS_TRANS for the geometry a `/` and `\` pair states.  The comments below
    # describe the tetrahedral reading; for a geometry `position` is one terminal of the chain of
    # double bonds, `refs` is (that terminal's marked substituent, the other terminal's, the other
    # terminal, Q_NO_SLOT) with `n_refs` 3, `sign` is 0, and `spare`'s high byte is the parity
    # demanded in the frame `(marked, other, marked, other)` -- 1 for trans.  See _seal_geometries.
    uint32_t position       # the DFS position of the atom the primitive sits on
    uint32_t readiness      # the DFS position at which the last of `refs` becomes mapped, so the
                            # greatest of `position` and the mapped `refs`; the kernel tests the
                            # primitive there and needs no "am I ready yet" test of its own
    uint32_t refs[4]        # the atom's query neighbours as DFS positions, in the QUERY's ruling-F26
                            # order (ascending query slot, i.e. creation order), Q_NO_SLOT-padded
    uint8_t sign            # the UNION of the atom's per-box signs (QSIGN_*), for reporting only.
                            # Ruling F87 moved the decision onto the box: an atom's boxes are a
                            # disjunction and each names its own configuration, so the kernel reads
                            # the box that admitted the candidate, never this field.  It is kept
                            # because it says at a glance what an atom's stereo content is, and
                            # because it is what makes the record's presence auditable from Python.
    uint8_t n_refs          # named query neighbours, 0..4.  Anything but 3 or 4 cannot describe a
                            # tetrahedral frame and the kernel refuses the match (ruling F77 case 2)
    uint16_t spare


cdef packed struct qcomp_t:      # 12 bytes
    uint32_t begin
    uint32_t end
    int32_t group


cdef packed struct QueryHeader:
    uint32_t magic
    uint32_t version
    uint32_t flags
    uint32_t atom_count
    uint32_t bond_count
    uint32_t box_count
    uint32_t component_count
    uint32_t automorphism_count      # rows in QSEG_AUTOMORPHISM, the identity excluded
    uint32_t stereo_count            # qstereo_t records in QSEG_STEREO
    uint32_t total_len
    uint32_t reserved[2]
    uint64_t signature[4]            # the screen's demand words; Task 14 fills it
    segment_t segments[QSEG_COUNT]


cdef class Query:
    cdef char *buffer
    cdef QueryHeader *header
    cdef size_t total_len

    def __cinit__(self):
        self.buffer = NULL
        self.header = NULL
        self.total_len = 0

    def __dealloc__(self):
        PyMem_Free(self.buffer)
        self.buffer = NULL

    cdef inline void *segment(self, int seg) noexcept nogil:
        if self.header.segments[seg].length == 0:
            return <void *> _zero_page
        return <void *> (self.buffer + self.header.segments[seg].offset)

    cdef inline qatom_t *atoms(self) noexcept nogil:
        return <qatom_t *> self.segment(QSEG_ATOMS)


cdef inline qbox_t *query_boxes(Query q) noexcept nogil:
    return <qbox_t *> q.segment(QSEG_BOXES)


cdef inline qany_t *query_any(Query q) noexcept nogil:
    return <qany_t *> q.segment(QSEG_ANY)


cdef inline qbond_t *query_bonds(Query q) noexcept nogil:
    return <qbond_t *> q.segment(QSEG_BONDS)


cdef inline qbox_t *query_bond_boxes(Query q) noexcept nogil:
    return <qbox_t *> q.segment(QSEG_BOND_BOXES)


cdef inline qclosure_t *query_closures(Query q) noexcept nogil:
    return <qclosure_t *> q.segment(QSEG_CLOSURES)


cdef inline qcomp_t *query_components(Query q) noexcept nogil:
    return <qcomp_t *> q.segment(QSEG_COMPONENTS)


cdef inline uint32_t *query_element_demand(Query q) noexcept nogil:
    return <uint32_t *> q.segment(QSEG_ELEMENT_DEMAND)


cdef inline uint32_t *query_demand_list(Query q) noexcept nogil:
    """Pointer to the compact (element, count) demand list; length from the segment header."""
    return <uint32_t *> q.segment(QSEG_DEMAND_LIST)


cdef inline uint32_t *query_automorphisms(Query q) noexcept nogil:
    return <uint32_t *> q.segment(QSEG_AUTOMORPHISM)


cdef inline qstereo_t *query_stereo(Query q) noexcept nogil:
    return <qstereo_t *> q.segment(QSEG_STEREO)


cdef Query query_alloc(uint32_t atom_count, uint32_t bond_count, uint32_t box_count,
                       uint32_t bond_box_count, uint32_t any_count, uint32_t closure_count,
                       uint32_t component_count, uint32_t automorphism_count,
                       uint32_t demand_list_len=0, uint32_t stereo_count=0):
    cdef size_t offset = sizeof(QueryHeader)
    cdef size_t atoms_len = align8(atom_count * sizeof(qatom_t))
    cdef size_t boxes_len = align8(box_count * sizeof(qbox_t))
    cdef size_t any_len = align8(any_count * sizeof(qany_t))
    cdef size_t bonds_len = align8(bond_count * sizeof(qbond_t))
    cdef size_t bond_boxes_len = align8(bond_box_count * sizeof(qbox_t))
    cdef size_t closures_len = align8(closure_count * sizeof(qclosure_t))
    cdef size_t components_len = align8(component_count * sizeof(qcomp_t))
    cdef size_t demand_len = 120 * sizeof(uint32_t)   # always present: 480 bytes, 8-aligned
    # demand_list: demand_list_len (element, count) pairs, each 2 x uint32_t = 8 bytes.  8 bytes
    # is already 8-aligned, so align8 is a no-op here; the expression makes the invariant explicit.
    cdef size_t demand_list_size = align8(<size_t> demand_list_len * 2 * sizeof(uint32_t))
    cdef size_t autos_len = align8(<size_t> automorphism_count * atom_count * sizeof(uint32_t))
    cdef size_t stereo_len = align8(<size_t> stereo_count * sizeof(qstereo_t))
    cdef size_t total = (offset + atoms_len + boxes_len + any_len + bonds_len +
                         bond_boxes_len + closures_len + components_len + demand_len +
                         demand_list_size + autos_len + stereo_len)

    cdef Query q = Query.__new__(Query)
    q.buffer = <char *> PyMem_Malloc(total)
    if q.buffer is NULL:
        raise MemoryError('query allocation failed')
    memset(q.buffer, 0, total)
    q.total_len = total
    q.header = <QueryHeader *> q.buffer
    q.header.magic = QUERY_MAGIC
    q.header.version = QUERY_VERSION
    q.header.atom_count = atom_count
    q.header.bond_count = bond_count
    q.header.box_count = box_count
    q.header.component_count = component_count
    q.header.automorphism_count = automorphism_count
    q.header.stereo_count = stereo_count
    q.header.total_len = <uint32_t> total

    # Lay segments out in index order.  Zero-count segments get offset = current end and
    # length = 0; Query.segment() routes them to the zero page.  QSEG_ELEMENT_DEMAND is
    # always present.  QSEG_AUTOMORPHISM is zero-length exactly when the group is trivial.
    q.header.segments[QSEG_ATOMS].offset = <uint32_t> offset
    q.header.segments[QSEG_ATOMS].length = <uint32_t> atoms_len
    offset += atoms_len

    q.header.segments[QSEG_BOXES].offset = <uint32_t> offset
    q.header.segments[QSEG_BOXES].length = <uint32_t> boxes_len
    offset += boxes_len

    q.header.segments[QSEG_ANY].offset = <uint32_t> offset
    q.header.segments[QSEG_ANY].length = <uint32_t> any_len
    offset += any_len

    q.header.segments[QSEG_BONDS].offset = <uint32_t> offset
    q.header.segments[QSEG_BONDS].length = <uint32_t> bonds_len
    offset += bonds_len

    q.header.segments[QSEG_BOND_BOXES].offset = <uint32_t> offset
    q.header.segments[QSEG_BOND_BOXES].length = <uint32_t> bond_boxes_len
    offset += bond_boxes_len

    q.header.segments[QSEG_CLOSURES].offset = <uint32_t> offset
    q.header.segments[QSEG_CLOSURES].length = <uint32_t> closures_len
    offset += closures_len

    q.header.segments[QSEG_COMPONENTS].offset = <uint32_t> offset
    q.header.segments[QSEG_COMPONENTS].length = <uint32_t> components_len
    offset += components_len

    q.header.segments[QSEG_ELEMENT_DEMAND].offset = <uint32_t> offset
    q.header.segments[QSEG_ELEMENT_DEMAND].length = <uint32_t> demand_len
    offset += demand_len

    q.header.segments[QSEG_DEMAND_LIST].offset = <uint32_t> offset
    q.header.segments[QSEG_DEMAND_LIST].length = <uint32_t> demand_list_size
    offset += demand_list_size

    q.header.segments[QSEG_AUTOMORPHISM].offset = <uint32_t> offset
    q.header.segments[QSEG_AUTOMORPHISM].length = <uint32_t> autos_len
    offset += autos_len

    q.header.segments[QSEG_STEREO].offset = <uint32_t> offset
    q.header.segments[QSEG_STEREO].length = <uint32_t> stereo_len

    return q


def _query_header_size():
    return sizeof(QueryHeader)


def _query_record_sizes():
    return {'qatom_t': sizeof(qatom_t), 'qbox_t': sizeof(qbox_t), 'qany_t': sizeof(qany_t),
            'qbond_t': sizeof(qbond_t), 'qclosure_t': sizeof(qclosure_t),
            'qcomp_t': sizeof(qcomp_t), 'qstereo_t': sizeof(qstereo_t)}


def _query_segment_ids():
    return {'QSEG_ATOMS': <int> QSEG_ATOMS, 'QSEG_BOXES': <int> QSEG_BOXES,
            'QSEG_ANY': <int> QSEG_ANY, 'QSEG_BONDS': <int> QSEG_BONDS,
            'QSEG_BOND_BOXES': <int> QSEG_BOND_BOXES, 'QSEG_CLOSURES': <int> QSEG_CLOSURES,
            'QSEG_COMPONENTS': <int> QSEG_COMPONENTS,
            'QSEG_ELEMENT_DEMAND': <int> QSEG_ELEMENT_DEMAND,
            'QSEG_DEMAND_LIST': <int> QSEG_DEMAND_LIST,
            'QSEG_AUTOMORPHISM': <int> QSEG_AUTOMORPHISM,
            'QSEG_STEREO': <int> QSEG_STEREO,
            'QSEG_COUNT': <int> QSEG_COUNT}
