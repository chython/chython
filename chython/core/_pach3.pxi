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
# ==================================================================================================
# pach VERSIONS 3 AND 4.  One design and one reader; version 3 carries display coordinates and
# version 4 does not.  `docs/pach.rst` is the wire layout a third party implements against and this
# comment is the same layouts stated where the code is.
#
# Header, 12 bytes:
#    0      version   u8    3 | 4
#    1      flags     u8    bit 0 = map block present; bits 1-7 reserved, must be 0
#    2-3    atoms     u16
#    4-5    bonds     u16
#    6-7    stereo    u16
#    8-9    sgroups   u16   enhanced-stereo entries
#    10-11  reserved  u16   must be 0
#
# Every count is in the header and every stride is constant, so a record's length is arithmetic and
# `pach_record_length` is O(1).  The header is the one place slack is deliberate: a future block needs
# a COUNT, and spare bits inside a fixed-stride record cannot hold one.
# ==================================================================================================

DEF PACH3_VERSION_XY = 3
DEF PACH3_VERSION_FLAT = 4
DEF PACH3_HEADER_LEN = 12
DEF PACH3_ATOM_XY_LEN = 9
DEF PACH3_ATOM_FLAT_LEN = 3
DEF PACH3_BOND_LEN = 5
DEF PACH3_STEREO_LEN = 9
DEF PACH3_SGROUP_LEN = 3
DEF PACH3_MAP_LEN = 2
DEF PACH3_FLAG_MAP = 0x01
# The isotope field is 6 bits spelling `MDL_ISOTOPE[z] - 32 + value` for value 1..63, and 0 for unset:
# mass numbers 31 below the element's MDL reference to 31 above it.  Measured maximum shift over
# chython's 436 nuclides is 8, so the field cannot fill.
DEF PACH3_ISOTOPE_BIAS = 32
# int24 two's complement at XY_SCALE: +/-838.8607 Angstrom, exact rather than rounded, so a decode and
# re-encode of a coordinate-bearing record is byte-stable.
DEF PACH3_XY_LIMIT = 8388607


cdef inline uint32_t _p3_u16(const unsigned char *p) noexcept nogil:
    return p[0] | (<uint32_t> p[1] << 8)


cdef inline uint32_t _p3_present(Py_ssize_t at, Py_ssize_t length, uint32_t declared,
                                 Py_ssize_t stride) noexcept nogil:
    """How many of a block's declared records are in the buffer.

    THE HEADER'S COUNTS DEFINE THE LAYOUT AND THE BUFFER'S LENGTH DECIDES HOW MUCH OF IT ARRIVED.  A
    block starts where the declared counts put it and stops at whichever comes first, its own count or
    the buffer's end -- so a missing byte costs the record it falls in and none of the records behind
    it, and a count larger than the block holds reads what follows the block as its own records.
    """
    cdef Py_ssize_t have
    if at >= length:
        return 0
    have = (length - at) // stride
    return declared if <Py_ssize_t> declared <= have else <uint32_t> have


cdef inline void _p3_put_u16(unsigned char *p, uint32_t v) noexcept nogil:
    p[0] = <unsigned char> (v & 0xff)
    p[1] = <unsigned char> ((v >> 8) & 0xff)


cdef inline int32_t _p3_i24(const unsigned char *p) noexcept nogil:
    cdef uint32_t v = p[0] | (<uint32_t> p[1] << 8) | (<uint32_t> p[2] << 16)
    if v & 0x800000:
        return <int32_t> v - 0x1000000
    return <int32_t> v


cdef inline void _p3_put_i24(unsigned char *p, int32_t v) noexcept nogil:
    cdef uint32_t u = <uint32_t> v & 0xffffff
    p[0] = <unsigned char> (u & 0xff)
    p[1] = <unsigned char> ((u >> 8) & 0xff)
    p[2] = <unsigned char> ((u >> 16) & 0xff)


cdef inline Py_ssize_t _pach3_size(uint32_t atoms, uint32_t bonds, uint32_t stereo, uint32_t sgroups,
                                  bint want_xy, bint want_map) noexcept nogil:
    """The length arithmetic the header states, from the counts rather than from a buffer."""
    cdef Py_ssize_t out = PACH3_HEADER_LEN \
        + <Py_ssize_t> atoms * (PACH3_ATOM_XY_LEN if want_xy else PACH3_ATOM_FLAT_LEN) \
        + <Py_ssize_t> bonds * PACH3_BOND_LEN \
        + <Py_ssize_t> stereo * PACH3_STEREO_LEN \
        + <Py_ssize_t> sgroups * PACH3_SGROUP_LEN
    if want_map:
        out += <Py_ssize_t> atoms * PACH3_MAP_LEN
    return out


cdef Py_ssize_t _pach3_length(const unsigned char *data, Py_ssize_t length) noexcept nogil:
    """Byte length of the version 3 or 4 record at `data`, or -1 when the header is not all there.

    Reads the header and nothing else: a length that had to walk the record would make a caller
    stepping through a concatenated store quadratic in the store.
    """
    if length < PACH3_HEADER_LEN:
        return -1
    return _pach3_size(_p3_u16(data + 2), _p3_u16(data + 4), _p3_u16(data + 6), _p3_u16(data + 8),
                       data[0] == PACH3_VERSION_XY, (data[1] & PACH3_FLAG_MAP) != 0)


cdef int _pach3_refuse_losses(MoleculeContainer mol, uint32_t drop_mask) except -1:
    """Everything the arena holds that versions 3 and 4 have no field for, refused by name.

    A conformer set is one of them: the coordinate block is 2D display geometry, so a 3D conformer is
    a loss this asks about rather than a coordinate it could write. Map numbers, wedges, stereo groups
    and stereo each have a block of their own and are not asked about here.
    """
    cdef Structure structure = mol._structure
    cdef atom_t *atoms = structure.atoms()
    cdef halfedge_t *edges = csr_edges(structure)
    cdef uint32_t i, k
    if not (drop_mask & PACH_DROP_CIP):
        for i in range(structure.header.atom_count):
            if atoms[i].reserved & ATOM_CIP_MASK:
                raise ValueError('atom %d carries a cip descriptor and the pach format has no field '
                                 'for one; pass drop=[\'cip\'] to write the record without it'
                                 % atoms[i].n)
        for k in range(2 * structure.header.bond_count):
            if edges[k].flags & HE_CIP_MASK:
                raise ValueError('a bond carries a cip descriptor and the pach format has no field '
                                 'for one; pass drop=[\'cip\'] to write the record without it')
    if not (drop_mask & PACH_DROP_SGROUPS) and structure_sgroup_count(structure):
        raise ValueError('this molecule carries %d sgroups and the pach format has no field for them; '
                         'pass drop=[\'sgroups\'] to write the record without them'
                         % structure_sgroup_count(structure))
    if not (drop_mask & PACH_DROP_CONFORMERS) and structure_has(structure, SEG_CONFORMERS):
        raise ValueError('this molecule carries 3D conformers and the pach format holds 2D display '
                         'coordinates only; pass drop=[\'conformers\'] to write the record without '
                         'them, or to_bytes() to keep them')
    if not (drop_mask & PACH_DROP_TITLE) and len(blob_bytes(structure, SEG_OPAQUE_BLOB, 0)):
        raise ValueError('this molecule carries a title and the pach format has no text of any kind; '
                         'pass drop=[\'title\'] to write the record without it')
    # `_meta` and not `meta`, so asking the question does not create the dict it is asking about
    if not (drop_mask & PACH_DROP_META) and mol._meta:
        raise ValueError('this molecule carries %d metadata key(s) and the pach format has no field '
                         'for any of them; pass drop=[\'meta\'] to write the record without them'
                         % len(mol._meta))
    return 0


cdef int _pach3_put_atom(atom_t *a, unsigned char *out) except -1:
    """The three bytes both versions share, at `out[0:3]`.

        byte 0   bit 7 clear: bits 6-0 are the atomic number.  Bit 7 set: bits 6-0 are the R index.
        byte 1   isotope 6 | radical 1 | h_pinned 1
        byte 2   implicit H 4 | charge + 4, 4

    No field crosses a byte boundary, and the hydrogen nibble is the arena's own: 0..14 and 15 for
    H_UNKNOWN, with no translation either way.
    """
    cdef int32_t shift, charge
    if a.element == 0:
        if a.isotope:
            raise ValueError('atom %d is an R marker carrying isotope %d, and the isotope field is a '
                             'shift from an element\'s reference mass' % (a.n, a.isotope))
        out[0] = <unsigned char> (0x80 | (at_r_index(a) & 0x7f))
    else:
        out[0] = <unsigned char> a.element
    out[1] = 0
    if a.isotope:
        shift = <int32_t> a.isotope - <int32_t> MDL_ISOTOPE[a.element] + PACH3_ISOTOPE_BIAS
        if shift < 1 or shift > 63:
            raise ValueError('atom %d carries isotope %d and the pach isotope field reaches 31 mass '
                             'numbers either side of element %d\'s reference %d'
                             % (a.n, a.isotope, a.element, MDL_ISOTOPE[a.element]))
        out[1] = <unsigned char> shift
    if at_radical(a):
        out[1] |= 0x40
    if at_h_pinned(a):
        out[1] |= 0x80
    charge = <int32_t> a.charge + 4
    if charge < 0 or charge > 15:
        raise ValueError('atom %d carries charge %d and the pach charge field holds -4 to +11'
                         % (a.n, a.charge))
    out[2] = <unsigned char> ((a.hydrogens & 0x0f) | (<unsigned char> charge << 4))
    return 0


cdef inline void _pach3_unit_frame(stereo_unit_t *u, uint32_t *want) noexcept nogil:
    """The unit's four directions in the ORDER THE RECORD USES, which is stated here and nowhere else.

        SU_TETRA     the three named directions in refs order, then the implied one -- refs[3] when
                     all four are named, and the unnamed direction otherwise
        bond kinds   per list, the named direction and then the list's other one, the anchor's list
                     first.  `SU_NO_REF` sits in the MIDDLE of the four slots for a bond kind -- an
                     implicit-hydrogen but-2-ene is (C, None, C, None) -- which is why this is built
                     per list rather than by compacting the four.

    Ruling F41 puts at most one unnamed direction in a list, so "the named one" is `refs[base]` unless
    that is `SU_NO_REF`.  F26 orders each list's unnamed direction last, so the frame EQUALS `refs` on
    every unit perception emits today and `smi_perm_of` answers the identity --
    `test_pach3.py:test_the_frame_is_refs_order_on_every_configured_unit` is that claim's harness.  The
    frame is built anyway because then a revised F26 cannot change what a stored record means.  The
    decoder inverts exactly this, and `translate_parity` is an XOR by the permutation's parity and
    therefore its own inverse, which makes the pair a fixed point.
    """
    cdef uint32_t j, base, r
    cdef uint32_t k = 0
    if u.kind == SU_TETRA:
        want[0] = SU_NO_REF; want[1] = SU_NO_REF; want[2] = SU_NO_REF; want[3] = SU_NO_REF
        for j in range(4):
            r = u.refs[j]
            if r != SU_NO_REF:
                want[k] = r
                k += 1
    else:
        for base in range(0, 4, 2):
            r = u.refs[base]
            if r == SU_NO_REF:
                want[base] = u.refs[base + 1]
                want[base + 1] = SU_NO_REF
            else:
                want[base] = r
                want[base + 1] = u.refs[base + 1]


cdef int _pach3_stereo_record(Structure structure, stereo_unit_t *u, uint8_t parity,
                              unsigned char *out) except -1:
    """One nine-byte stereo record, or nothing written and 0 returned when the unit is unconfigured.

        bytes 0-1  slot0    bytes 2-3  slot1    bytes 4-5  slot2    bytes 6-7  slot3
        byte 8     kind 3 | sign 1 | reserved 4

        kind 0 SU_TETRA         centre  d0      d1      d2
        kind 1 SU_CIS_TRANS     end A   a0      end B   b0
        kind 2 SU_ALLENE        end A   a0      end B   b0     (the chain centre is not stored)
        kind 3 SU_ATROPISOMER   pivot A oA0     pivot B oB0

    A slot is an ATOM INDEX in the record's own atom block, which is arena slot order.  A
    configuration states the owner of each direction list and all but one of its directions; the last
    is implied by identity -- the owner's direction not already named -- so there are no masks and no
    sentinels and every slot holds a real atom.

    `sign` is one bit because a record exists only for a configured unit: 0 is the even parity and 1
    the odd one, and the unset third state needs no spelling.  The reserved nibble is where a
    three-state parity would go without touching the stride.  The sign is the parity in the frame
    `_pach3_unit_frame` builds, which is what lets ruling F26 be revised without invalidating a stored
    record: nothing in the record depends on the direction ORDERING, only on identity.
    """
    cdef atom_t *atoms = structure.atoms()
    cdef uint32_t *ptr = csr_ptr(structure)
    cdef halfedge_t *edges = csr_edges(structure)
    cdef uint8_t kind = u.kind
    cdef uint32_t anchor = u.anchor
    cdef uint32_t owner_a, owner_b
    cdef uint32_t terms[2]
    cdef uint32_t inwards[2]
    cdef uint32_t want[4]
    cdef uint32_t perm[4]
    if parity == 0:
        return 0
    if u.n_refs != 4:
        # UNREACHABLE THROUGH PERCEPTION and kept anyway: `_stereo.pxi`'s "WHY N==4 IS THE ONLY CASE"
        # states that all four `_stereo_emit` sites pass `n_refs=4` and that the one site which does
        # not, `_stereo_anchor_collision_probe`, has its unit refused before translation.  That is an
        # invariant of another file over a table this one only reads, so the guard states it here.
        raise ValueError('the stereo unit anchored at atom %d orders %d reference direction(s) and a '
                         'parity is a fact about four; pass drop=[\'stereo\'] to write the record '
                         'without it' % (atoms[anchor].n, u.n_refs))
    _pach3_unit_frame(u, want)
    if kind == SU_TETRA:
        if want[2] == SU_NO_REF:
            raise ValueError('the tetrahedral centre at atom %d names %d direction(s) with an atom of '
                             'their own and the pach record states three; pass drop=[\'stereo\'] to '
                             'write the record without it'
                             % (atoms[anchor].n, 1 if want[1] == SU_NO_REF else 2))
        _p3_put_u16(out, anchor)
        _p3_put_u16(out + 2, want[0])
        _p3_put_u16(out + 4, want[1])
        _p3_put_u16(out + 6, want[2])
    else:
        # A list with TWO unnamed directions has no second slot to name and, per ruling F41
        # (`_list_has_two_unnamed`), is not stereogenic -- so a configured unit cannot present one.
        if want[0] == SU_NO_REF or want[2] == SU_NO_REF:
            raise ValueError('the stereo unit anchored at atom %d has a direction list with no named '
                             'direction, so the pach record has nothing to state it against; pass '
                             'drop=[\'stereo\'] to write the record without it' % atoms[anchor].n)
        if kind == SU_ALLENE:
            # The anchor is the chain's CENTRE and the record names the two ENDS, so the ends come
            # from the chain and which end leads comes from which one owns the first named direction.
            if not _pach_allene_ends(ptr, edges, anchor, terms, inwards):
                raise ValueError('atom %d anchors an allene whose chain this molecule does not hold; '
                                 'pass drop=[\'stereo\'] to write the record without it'
                                 % atoms[anchor].n)
            if csr_find_at(ptr, edges, terms[0], want[0]) is not NULL:
                owner_a = terms[0]
                owner_b = terms[1]
            else:
                owner_a = terms[1]
                owner_b = terms[0]
        else:
            owner_a = anchor
            owner_b = stereo_unit_partner(structure, u)
            if owner_b == SU_NO_REF:
                raise ValueError('the stereo unit anchored at atom %d names two atoms and the second '
                                 'is not in this molecule; pass drop=[\'stereo\'] to write the record '
                                 'without it' % atoms[anchor].n)
        _p3_put_u16(out, owner_a)
        _p3_put_u16(out + 2, want[0])
        _p3_put_u16(out + 4, owner_b)
        _p3_put_u16(out + 6, want[2])
    smi_perm_of(u, want, perm)
    out[8] = <unsigned char> (kind | ((translate_parity(parity, perm) - 1) << 3))
    return 1


cdef int _pach3_stereo_block(Structure structure, unsigned char *out, uint32_t *count) except -1:
    """Every configured unit as a nine-byte record, into room for `unit_count` of them.

    Called only when `stereo` was not dropped, so a unit this format cannot state is a refusal here
    rather than a skip -- `drop=['stereo']` is what a caller who wants the record anyway passes.  The
    table holds unconfigured units too, which is why the header's count is what was WRITTEN.  One unit
    per anchor, so that count cannot exceed the atom count and the header's u16 field cannot fill.
    """
    cdef stereo_unit_t *units = structure_stereo_units(structure)
    cdef stereo_unit_t *u
    cdef uint32_t total = structure_stereo_unit_count(structure)
    cdef uint32_t i
    cdef uint32_t written = 0
    for i in range(total):
        u = units + i
        written += <uint32_t> _pach3_stereo_record(structure, u,
                                                  structure_parity_at(structure, u.anchor),
                                                  out + written * PACH3_STEREO_LEN)
    count[0] = written
    return 0


cdef bytes _pach3_encode(MoleculeContainer mol, uint32_t drop_mask):
    """One version 3 or version 4 pach record, uncompressed.

    Version 3 when the molecule has coordinates and `coordinates` was not dropped, else version 4.
    ATOM ORDER IS ARENA SLOT ORDER and stable ids are not written: bonds and stereo slots address
    positions in the atom block, so a reader reproduces the order by reading it, and the record has no
    12-bit id ceiling to run into.
    """
    cdef Structure structure = mol._structure
    cdef atom_t *atoms
    cdef uint32_t *ptr
    cdef halfedge_t *edges
    cdef halfedge_t *rev
    cdef xy_t *xy
    cdef uint8_t *groups
    cdef uint32_t n, nb, i, j, k, a1, a2
    cdef uint32_t unit_count, stereo_count, sgroup_count
    cdef uint8_t sg_value
    cdef bint want_xy
    cdef bint want_map
    cdef bint drop_wedges
    cdef unsigned char wedge
    cdef uint32_t wedged
    cdef Py_ssize_t alloc, at
    cdef unsigned char *buf
    cdef bytes out

    mol._require_clean()
    _pach3_refuse_losses(mol, drop_mask)
    # DERIVED and unmarked, and it can reallocate, so every pointer below is taken after it.  The
    # decoder asks the same question, and the two directions have to ask the same one or a stored
    # configuration means two things.
    ensure_stereo_units_unmarked(structure)
    structure = mol._structure
    atoms = structure.atoms()
    ptr = csr_ptr(structure)
    edges = csr_edges(structure)
    n = structure.header.atom_count
    nb = structure.header.bond_count
    if n > 65535:
        raise ValueError('this molecule has %d atoms and the pach atom count is a 16 bit field; use '
                         'to_bytes(), which is lossless' % n)
    if nb > 65535:
        raise ValueError('this molecule has %d bonds and the pach bond count is a 16 bit field; use '
                         'to_bytes(), which is lossless' % nb)
    want_xy = structure_has(structure, SEG_XY) and not (drop_mask & PACH_DROP_COORDINATES)
    xy = NULL
    if want_xy:
        xy = structure_xy(structure)
    unit_count = 0
    if not (drop_mask & PACH_DROP_STEREO):
        unit_count = structure_stereo_unit_count(structure)
    # EITHER NAME DROPS THIS BLOCK.  `drop=['stereo']` drops it along with the configurations, a caller
    # asking for a record without stereo meaning without stereo; `drop=['stereo_groups']` drops it alone
    # and is the name the version 0 and 2 writers' own refusal tells a caller to pass.
    # `drop=['sgroups']` is a different field -- `_pach3_refuse_losses` asks about the CTfile S-group
    # records, which have no block here at all.
    groups = NULL
    sgroup_count = 0
    if structure_has(structure, SEG_STEREO_GROUPS) \
            and not (drop_mask & (PACH_DROP_STEREO | PACH_DROP_STEREO_GROUPS)):
        groups = structure_stereo_groups(structure)
        for i in range(n):
            if groups[i]:
                sgroup_count += 1
        if not sgroup_count:
            groups = NULL

    want_map = False
    if not (drop_mask & PACH_DROP_MAP_NUMBER):
        for i in range(n):
            if atoms[i].map_number:
                want_map = True
                break
    # ONE ALLOCATION, at the stereo block's upper bound: the table holds unconfigured units, so how
    # many records there are is known only once they are written.  The stereo block is written in
    # place at its own offset and the record is the prefix `buf[:at]`, which is what a scratch buffer
    # and a memcpy would otherwise be for.
    alloc = _pach3_size(n, nb, unit_count, sgroup_count, want_xy, want_map)
    buf = <unsigned char *> PyMem_Malloc(alloc)
    if buf is NULL:
        raise MemoryError('pach record allocation failed')
    try:
        memset(buf, 0, alloc)
        buf[0] = PACH3_VERSION_XY if want_xy else PACH3_VERSION_FLAT
        buf[1] = PACH3_FLAG_MAP if want_map else 0
        _p3_put_u16(buf + 2, n)
        _p3_put_u16(buf + 4, nb)
        _p3_put_u16(buf + 8, sgroup_count)
        at = PACH3_HEADER_LEN
        for i in range(n):
            _pach3_put_atom(&atoms[i], buf + at)
            if want_xy:
                if xy[i].x < -PACH3_XY_LIMIT or xy[i].x > PACH3_XY_LIMIT \
                        or xy[i].y < -PACH3_XY_LIMIT or xy[i].y > PACH3_XY_LIMIT:
                    raise ValueError('atom %d sits at (%r, %r) and the pach coordinate field reaches '
                                     '+/-838.8607; pass drop=[\'coordinates\'] to write the record '
                                     'without a drawing'
                                     % (atoms[i].n, xy_read_x(&xy[i]), xy_read_y(&xy[i])))
                _p3_put_i24(buf + at + 3, xy[i].x)
                _p3_put_i24(buf + at + 6, xy[i].y)
            at += PACH3_ATOM_XY_LEN if want_xy else PACH3_ATOM_FLAT_LEN
        drop_wedges = (drop_mask & PACH_DROP_WEDGES) != 0 or not want_xy
        # A version 4 record has a wedge nibble and no drawing for it to mean anything against, so the
        # nibble is written 0 and the loss is on the molecule's log rather than raised: dropping the
        # coordinates is what the caller asked for and the wedge went with them.
        wedged = 0
        if not want_xy and not (drop_mask & PACH_DROP_WEDGES):
            for k in range(2 * nb):
                if edges[k].wedge:
                    wedged += 1
            if wedged:
                mol._log_event('pach:wedge-lost', 'pach',
                               'version 4 carries no coordinates, so %d wedge(s) were not written'
                               % wedged, mc_lost())
        for i in range(n):
            for k in range(ptr[i], ptr[i + 1]):
                j = edges[k].to
                if j < i:
                    continue                       # the other half-edge wrote this bond
                a1 = i
                a2 = j
                wedge = 0
                if not drop_wedges:
                    wedge = edges[k].wedge
                    if wedge == 0:
                        # the wedge's narrow end is j, and the pair is written narrow end first
                        rev = csr_find_at(ptr, edges, j, i)
                        if rev is not NULL and rev.wedge:
                            a1 = j
                            a2 = i
                            wedge = rev.wedge
                _p3_put_u16(buf + at, a1)
                _p3_put_u16(buf + at + 2, a2)
                buf[at + 4] = <unsigned char> (edges[k].order | (wedge << 4))
                at += PACH3_BOND_LEN
        stereo_count = 0
        if unit_count:
            _pach3_stereo_block(structure, buf + at, &stereo_count)
            _p3_put_u16(buf + 6, stereo_count)
            at += <Py_ssize_t> stereo_count * PACH3_STEREO_LEN
        # ENHANCED STEREO, three bytes an entry and only for an atom that carries one: `atom u16` and
        # the arena's own packed group byte, kind in bits 7-6 and group in bits 5-0.  Its own block
        # rather than a field in a stereo record because `set_stereo_group` accepts ANY atom, so an
        # atom owning no unit has no record to hold it.
        if groups is not NULL:
            for i in range(n):
                sg_value = groups[i]
                if sg_value:
                    _p3_put_u16(buf + at, i)
                    buf[at + 2] = sg_value
                    at += PACH3_SGROUP_LEN
        if want_map:
            for i in range(n):
                _p3_put_u16(buf + at, atoms[i].map_number)
                at += PACH3_MAP_LEN
        # `at` and not `alloc`: the record ends where the last block ended, and an unconfigured unit
        # left room the header does not declare.  `pach_record_length` reads the counts, so the two
        # have to agree -- `test_the_declared_length_is_the_buffer_length` is where they are compared.
        out = <bytes> buf[:at]
    finally:
        PyMem_Free(buf)
    return out


cdef struct pach3_atom_t:
    uint8_t element            # 0 for an R marker
    uint8_t r_index
    int8_t charge
    uint8_t radical
    uint8_t pinned
    uint8_t hydrogens          # the arena's own nibble, H_UNKNOWN included
    uint16_t isotope           # absolute mass number, 0 for unset
    int32_t x
    int32_t y


cdef MoleculeContainer _pach3_build(pach3_atom_t *pa, uint32_t atoms_count, edge_edit_t *edits,
                                   uint32_t bonds_count, bint want_xy, bint want_parity,
                                   const uint8_t *groups, const uint16_t *maps):
    """Lay out the arena and wrap it in a container.

    STABLE IDS ARE NOT STORED, so the numbers are 1..n in record order.  `pack()` renumbers and the
    arena's `n` is a container's private label; that is what removes the 12-bit id ceiling rather than
    widening a field.  `rebuild_derived` REALLOCATES, which is why nothing here reads `atoms` after
    it -- the numbers come from the loop index, not from the buffer (ruling F60).

    `want_parity` is the header's DECLARED stereo count, not the number of configurations that turn out
    to resolve: the caller has the header before it has the graph a record resolves against, and the
    price of asking the earlier question is one byte per atom on a record whose every configuration is
    then dropped.

    `groups` is one packed group byte per atom in atom order, or NULL when no entry survived reading --
    the enhanced-stereo block names an atom rather than a unit, so the caller resolves it against the
    atom block alone and hands the finished row over.

    `maps` is a map number per atom in atom order, or NULL when all atoms are unmapped -- written
    beside `atoms[i].n` as a persistent field no derived segment reads.
    """
    cdef Structure structure
    cdef atom_t *atoms
    cdef pach3_atom_t *src
    cdef atom_t *dst
    cdef xy_t *xy
    cdef xy_t *xyp
    cdef uint32_t i, seg_mask = 0
    cdef int rc
    cdef list numbers = []
    cdef dict index_of = {}
    cdef MoleculeContainer mol
    if want_xy:
        seg_mask = SEG_MASK_XY
    if want_parity:
        # `_pach3_apply_stereo` runs after this function seals the arena, so the segment it writes into
        # is named here or nowhere (the persistent block is laid out once).
        seg_mask |= SEG_MASK_PARITY
    if groups is not NULL:
        seg_mask |= SEG_MASK_STEREO
    structure = structure_alloc_full(atoms_count, bonds_count, False, seg_mask, NULL)
    atoms = structure.atoms()
    for i in range(atoms_count):
        src = pa + i
        dst = atoms + i
        dst.element = src.element
        dst.charge = src.charge
        dst.isotope = src.isotope
        dst.n = i + 1
        if maps is not NULL:
            dst.map_number = maps[i]
        # the explicit nibble is derived and `rebuild_derived` fills it from the CSR
        at_set_h(dst, src.hydrogens, 0)
        if src.element == 0:
            at_set_r_index(dst, src.r_index)
        if src.radical:
            at_set_radical(dst, True)
        if src.pinned:
            at_set_h_pinned(dst, True)
    with nogil:
        rc = csr_build(structure, edits, bonds_count)
    if rc:
        raise MemoryError('csr scratch allocation failed')
    if want_xy:
        xy = structure_xy(structure)
        for i in range(atoms_count):
            src = pa + i
            xyp = xy + i
            xyp.x = src.x
            xyp.y = src.y
    if groups is not NULL:
        # `sg_len` is `align8(atom_count)`, so one byte an atom is inside the segment
        memcpy(structure_stereo_groups(structure), groups, atoms_count)
    rebuild_derived(structure)
    for i in range(atoms_count):
        numbers.append(i + 1)
        index_of[i + 1] = i
    mol = MoleculeContainer.__new__(MoleculeContainer)
    mol._structure = structure
    mol._numbers = numbers
    mol._index_of = index_of
    mol._next_id = atoms_count + 1
    mol._first_pending = atoms_count + 1
    return mol


cdef inline bint _pach3_seen(edge_edit_t *edits, uint32_t count, uint32_t a1,
                             uint32_t a2) noexcept nogil:
    """Whether this pair is already in the kept-bond array.  Linear, and deliberately so: a duplicate
    is a damage report and the array is the record's own bond count."""
    cdef uint32_t k
    cdef edge_edit_t *e
    for k in range(count):
        e = edits + k
        if (e.src == a1 and e.dst == a2) or (e.src == a2 and e.dst == a1):
            return True
    return False


cdef int _pach3_apply_wedges(MoleculeContainer mol, edge_edit_t *edits, uint8_t *wedges,
                             uint32_t count) except -1:
    """Each record's wedge onto the half-edge LEAVING its `a1`, which is the narrow end.

    After the build, not before: `csr_build` clears both halves' wedge, and `rebuild_derived`
    reallocates the persistent buffer, so the pointers have to be taken here.
    """
    cdef Structure structure = mol._structure
    cdef uint32_t *ptr = csr_ptr(structure)
    cdef halfedge_t *edges = csr_edges(structure)
    cdef halfedge_t *he
    cdef edge_edit_t *e
    cdef uint32_t k
    for k in range(count):
        if not wedges[k]:
            continue
        e = edits + k
        he = csr_find_at(ptr, edges, e.src, e.dst)
        if he is not NULL:
            he.wedge = wedges[k]
    return 0


# THE PERMUTATION IS BUILT FROM refs INDICES AND NEVER FROM DIRECTION VALUES, which is the one place
# the decoder does not mirror the writer's spelling.  `_pach3_unit_frame` hands `smi_perm_of` a frame
# whose direction lists are still in `refs` order, and there its in-order pairing of the anonymous
# `SU_NO_REF` positions is exact.  A record read from its FAR owner (the fallback below) presents the
# two lists swapped, and then the two anonymous positions swap with them -- so pairing them in order
# adds a transposition and inverts the parity of every unit whose both lists carry an unnamed
# direction, which is every `C/C=C/C`.  An index is not anonymous, so this asks the question that has
# one answer: which `refs` slot is the record's i-th direction.

cdef inline bint _p3_pair_index(stereo_unit_t *u, uint32_t base, uint32_t named,
                                uint32_t *index) noexcept nogil:
    """The `refs` index at which the list `refs[base:base+2]` holds `named`.

    False when the list does not hold `named` at all, which is a record disagreeing with the graph it
    was decoded into.  The list's other index is `2 * base + 1 - index[0]`, so naming one names both.
    The two entries are hoisted rather than addressed because `refs` lives in a packed struct
    (RULES.md 2.3).
    """
    cdef uint32_t r0 = u.refs[base]
    cdef uint32_t r1 = u.refs[base + 1]
    if r0 == named:
        index[0] = base
        return True
    if r1 == named:
        index[0] = base + 1
        return True
    return False


cdef inline bint _p3_tetra_indices(stereo_unit_t *u, const uint32_t *named3,
                                   uint32_t *perm) noexcept nogil:
    """`perm[0:3]` the `refs` index of each direction the record names, `perm[3]` the one left over.

    False when the three the record names are not three of the unit's four.  `used` is a match mask, so
    a name is consumed by one slot only -- an unnamed direction is `SU_NO_REF` in `refs` and the record
    never states one, but the mask is what keeps that a property of the loop rather than of the input.
    The four entries are hoisted rather than addressed because `refs` lives in a packed struct
    (RULES.md 2.3), and the inner loop reads all four of them.
    """
    cdef uint32_t j, k
    cdef uint32_t used = 0
    cdef bint found
    cdef uint32_t refs[4]
    for j in range(4):
        refs[j] = u.refs[j]
    perm[0] = 0; perm[1] = 0; perm[2] = 0; perm[3] = 0
    for k in range(3):
        found = False
        for j in range(4):
            if refs[j] == named3[k] and not (used & (1u << j)):
                perm[k] = j
                used |= 1u << j
                found = True
                break
        if not found:
            return False
    for j in range(4):
        if not (used & (1u << j)):
            perm[3] = j
            return True
    return False


cdef inline stereo_unit_t *_p3_allene_unit(Structure structure, uint32_t a, uint32_t b,
                                           uint32_t *base_a) noexcept nogil:
    """The allene unit holding `a` in one direction list and `b` in the other, with which list is `a`'s.

    An allene's anchor is the chain CENTRE, which the record does not store, so this is the one kind
    looked up by refs membership rather than by anchor.  The chain is not walked either way.

    The four entries are hoisted rather than addressed because `refs` lives in a packed struct
    (RULES.md 2.3), and each candidate unit reads all four of them twice.
    """
    cdef stereo_unit_t *units = structure_stereo_units(structure)
    cdef stereo_unit_t *u
    cdef uint32_t total = structure_stereo_unit_count(structure)
    cdef uint32_t i, r0, r1, r2, r3
    for i in range(total):
        u = units + i
        if u.kind != SU_ALLENE:
            continue
        r0 = u.refs[0]
        r1 = u.refs[1]
        r2 = u.refs[2]
        r3 = u.refs[3]
        if (r0 == a or r1 == a) and (r2 == b or r3 == b):
            base_a[0] = 0
            return u
        if (r2 == a or r3 == a) and (r0 == b or r1 == b):
            base_a[0] = 2
            return u
    return NULL


cdef int _pach3_apply_stereo(MoleculeContainer mol, const unsigned char *data, Py_ssize_t at,
                             uint32_t count, uint32_t n, list problems) except -1:
    """The stereo block onto a built molecule: one configuration per nine bytes, each dropped alone.

    `_pach3_stereo_record` inverted.  The unit is found by IDENTITY and never by walking a chain -- by
    anchor for a centre, by anchor at either owner slot for the two axis kinds, by refs membership for
    an allene, whose anchor is the chain centre and is not stored.  The record's direction frame is then
    rebuilt as a permutation of that unit's `refs` -- `_pach3_unit_frame` is the writer's half and the
    two are read together -- and `translate_parity` runs back.  It is an XOR by the permutation's
    parity, so the writer's call and this one are the same call and the pair is a fixed point rather
    than a near miss.

    A record that does not resolve costs its own configuration and nothing else -- a loop over a
    million stored records must not stop at one of them.
    """
    cdef Structure structure
    cdef atom_t *atoms
    cdef uint32_t *ptr
    cdef halfedge_t *edges
    cdef stereo_unit_t *u
    cdef uint32_t k, slot0, slot1, slot2, slot3, kind, sign, reserved, base_a, base_b
    cdef uint32_t named[3]
    cdef uint32_t perm[4]
    cdef bint wrote = False
    cdef object err
    if count == 0:
        return 0
    structure = mol._structure
    # DERIVING THE UNIT TABLE CAN FAIL and must not become this function's exception: the core refuses
    # a graph in which two units claim one anchor, and a damaged record can decode to such a graph.
    try:
        ensure_stereo_units_unmarked(structure)
    except Exception as err:
        problems.append('this record\'s graph has no usable stereo unit table (%s), so the %d '
                        'configuration(s) it states were dropped' % (err, count))
        return 0
    atoms = structure.atoms()                      # after, not before: the ensure may reallocate
    ptr = csr_ptr(structure)
    edges = csr_edges(structure)
    for k in range(count):
        slot0 = _p3_u16(data + at)
        slot1 = _p3_u16(data + at + 2)
        slot2 = _p3_u16(data + at + 4)
        slot3 = _p3_u16(data + at + 6)
        kind = data[at + 8] & 0x07
        sign = (data[at + 8] >> 3) & 1
        reserved = data[at + 8] >> 4
        at += PACH3_STEREO_LEN
        if reserved:
            problems.append('stereo record %d byte 8 bits 4-7 are reserved and must be 0; %d is '
                            'ignored' % (k, reserved))
        if slot0 >= n or slot1 >= n or slot2 >= n or slot3 >= n:
            problems.append('stereo record %d names atom index %d and this record has %d atom(s); the '
                            'configuration was dropped'
                            % (k, max(max(slot0, slot1), max(slot2, slot3)), n))
            continue
        if kind > SU_ATROPISOMER:
            problems.append('stereo record %d states kind %d and pach has four, 0 tetrahedral, 1 '
                            'cis/trans, 2 allene and 3 atropisomer; the configuration was dropped'
                            % (k, kind))
            continue
        base_a = 0
        if kind == SU_TETRA:
            u = stereo_unit_of(structure, slot0)
            if u is not NULL and u.kind != SU_TETRA:
                u = NULL
        elif kind == SU_ALLENE:
            u = _p3_allene_unit(structure, slot1, slot3, &base_a)
            # Slot0 and slot2 are the terminals, which the refs search did not read, so they are an
            # adjacency cross-check: one CSR lookup per end and no chain walk.
            if u is not NULL and (csr_find_at(ptr, edges, slot0, slot1) is NULL
                                  or csr_find_at(ptr, edges, slot2, slot3) is NULL):
                u = NULL
        else:
            # Ruling F45 relocates an axis' anchor to the other pivot, so the unit is looked for at
            # BOTH ends and the record stays readable either way round.  When it is the far owner that
            # anchors, `refs[0:2]` is slot2's list and slot1's is `refs[2:4]`.
            u = stereo_unit_of(structure, slot0)
            if u is NULL or u.kind != kind or stereo_unit_partner(structure, u) != slot2:
                u = stereo_unit_of(structure, slot2)
                base_a = 2
                if u is NULL or u.kind != kind or stereo_unit_partner(structure, u) != slot0:
                    u = NULL
        if u is NULL:
            if kind == SU_TETRA:
                problems.append('stereo record %d states a tetrahedral centre at atom %d, which '
                                'anchors no stereo unit of that kind in this molecule; the '
                                'configuration was dropped' % (k, atoms[slot0].n))
            else:
                problems.append('stereo record %d states %s over atoms %d and %d, which anchor no '
                                'stereo unit of that kind in this molecule; the configuration was '
                                'dropped'
                                % (k, smi_kind_name(<uint8_t> kind), atoms[slot0].n, atoms[slot2].n))
            continue
        if u.n_refs != 4:
            # The encoder's twin, and unreachable for the same reason: `_stereo.pxi`'s "WHY N==4 IS THE
            # ONLY CASE" states that every perception emit site passes `n_refs=4`.  A decode path must
            # not assume an invariant another file maintains over a table rebuilt from a forged
            # record's own graph, so the guard stays and reports rather than raising.
            problems.append('stereo record %d resolves to a unit ordering %d reference direction(s) '
                            'and a parity is a fact about four; the configuration was dropped'
                            % (k, u.n_refs))
            continue
        if structure_parity_at(structure, u.anchor):
            problems.append('stereo record %d resolves to the unit anchored at atom %d, which an '
                            'earlier record already configured; the later one was dropped'
                            % (k, atoms[u.anchor].n))
            continue
        # THE RECORD'S FRAME, as a permutation of `refs`: slot1 and slot3 are the named directions and
        # the implied ones are whichever slots those leave.
        if kind == SU_TETRA:
            named[0] = slot1
            named[1] = slot2
            named[2] = slot3
            if not _p3_tetra_indices(u, named, perm):
                problems.append('stereo record %d names three directions of the centre at atom %d and '
                                'at least one is not one of its four; the configuration was dropped'
                                % (k, atoms[slot0].n))
                continue
        else:
            base_b = 2 - base_a
            if not _p3_pair_index(u, base_a, slot1, &perm[0]) \
                    or not _p3_pair_index(u, base_b, slot3, &perm[2]):
                problems.append('stereo record %d names directions %d and %d, which are not the ones '
                                'the unit over atoms %d and %d orders; the configuration was dropped'
                                % (k, atoms[slot1].n, atoms[slot3].n, atoms[slot0].n, atoms[slot2].n))
                continue
            perm[1] = 2 * base_a + 1 - perm[0]
            perm[3] = 2 * base_b + 1 - perm[2]
        structure_set_parity(structure, u.anchor, translate_parity(<uint8_t> (sign + 1), perm))
        wrote = True
    if wrote:
        refresh_parity_features(structure)
    return 0


cdef tuple _pach3_decode(const unsigned char *data, Py_ssize_t length):
    """One version 3 or version 4 record.  `(MoleculeContainer or None, problems)`.

    MemoryError is the one exception that escapes, and it is not a statement about the bytes.  The
    graph derivation's own refusals are caught at the build -- `_pach_derivation_lost` -- because they
    ARE a statement about the bytes: a well-formed record can state a graph `perceive_rings` refuses.

    THE ATOM BLOCK IS THE ONE ALL-OR-NOTHING PART.  Bonds, stereo and the rest address positions in
    it, so a half-read atom block would make every later reference mean an atom that is not there;
    everything after it is read as far as the buffer goes and the shortfall is reported.
    """
    cdef list problems = []
    cdef unsigned char version, flags
    cdef uint32_t n, i, iso_field, chg, r_index
    cdef bint want_xy
    cdef Py_ssize_t stride, at
    cdef pach3_atom_t *pa
    cdef pach3_atom_t *dst
    cdef MoleculeContainer mol
    cdef edge_edit_t *edits
    cdef uint8_t *wedges
    cdef edge_edit_t *e
    cdef uint32_t nb, nb_declared, ns_declared, nsg_declared, nm_declared, k, a1, a2, order, wedge
    cdef uint32_t nb_have, ns_have, nsg_have, nm_have
    cdef uint32_t slot, value, kind, group
    cdef Py_ssize_t bonds_at, stereo_at, groups_at, map_at
    cdef size_t pa_len, edits_len, wedges_len, sg_len, maps_len
    cdef uint8_t *sg
    cdef uint16_t *maps = NULL
    cdef bint sg_any = False
    cdef bint map_any = False
    cdef unsigned char *block = NULL
    cdef object err

    if length < PACH3_HEADER_LEN:
        problems.append('a version 3 or 4 pach record is at least a 12 byte header and this buffer is '
                        '%d byte(s)' % length)
        return (None, problems)
    version = data[0]
    flags = data[1]
    want_xy = version == PACH3_VERSION_XY
    stride = PACH3_ATOM_XY_LEN if want_xy else PACH3_ATOM_FLAT_LEN
    n = _p3_u16(data + 2)
    nb_declared = _p3_u16(data + 4)
    ns_declared = _p3_u16(data + 6)
    nsg_declared = _p3_u16(data + 8)
    if flags & ~<unsigned char> PACH3_FLAG_MAP:
        problems.append('header flags is 0x%02x and only bit 0 is defined; the rest are ignored'
                        % flags)
    if _p3_u16(data + 10):
        problems.append('header bytes 10-11 are reserved and must be 0; %d is ignored'
                        % _p3_u16(data + 10))
    if PACH3_HEADER_LEN + <Py_ssize_t> n * stride > length:
        problems.append('the header declares %d atoms and the atom block does not fit in %d byte(s)'
                        % (n, length))
        return (None, problems)
    if not n:
        # In block order, and the map block needs no sentence of its own: its count is `n`.
        if nb_declared:
            problems.append('the header declares %d bond(s) and no atoms for them to name; none were '
                            'read' % nb_declared)
        if ns_declared:
            problems.append('the header declares %d stereo configuration(s) and no atoms for them to '
                            'name; none were read' % ns_declared)
        if nsg_declared:
            problems.append('the header declares %d enhanced-stereo entr(ies) and no atoms for them to '
                            'name; none were read' % nsg_declared)
        try:
            return (_pach3_build(NULL, 0, NULL, 0, want_xy, False, NULL, NULL), problems)
        except ValueError as err:
            _pach_derivation_lost(problems, err)
            return (None, problems)
    # One block carries every scratch region so there is one allocation, one NULL check and one free.
    # THE BOND REGIONS ARE SIZED FROM WHAT THE BUFFER HOLDS, not from the declared count: the loop below
    # ranges over `nb_have`, and a declared 65535 in a 15 byte buffer would otherwise allocate 851,960
    # bytes of the two of them.  Every region that can be empty gets zero bytes and a NULL pointer, so a
    # stray write faults instead of landing in the region behind it (RULES.md 5.3).
    bonds_at = PACH3_HEADER_LEN + <Py_ssize_t> n * stride
    nb_have = _p3_present(bonds_at, length, nb_declared, PACH3_BOND_LEN)
    pa_len = align8(n * sizeof(pach3_atom_t))
    edits_len = align8(nb_have * sizeof(edge_edit_t))
    wedges_len = align8(nb_have * sizeof(uint8_t))
    sg_len = align8(n * sizeof(uint8_t)) if nsg_declared else 0
    maps_len = align8(n * sizeof(uint16_t)) if flags & PACH3_FLAG_MAP else 0
    block = <unsigned char *> PyMem_Malloc(pa_len + edits_len + wedges_len + sg_len + maps_len)
    if block is NULL:
        raise MemoryError('pach record scratch allocation failed')
    pa = <pach3_atom_t *> block
    edits = <edge_edit_t *> (block + pa_len) if edits_len else NULL
    wedges = <uint8_t *> (block + pa_len + edits_len) if wedges_len else NULL
    sg = <uint8_t *> (block + pa_len + edits_len + wedges_len) if sg_len else NULL
    maps = <uint16_t *> (block + pa_len + edits_len + wedges_len + sg_len) if maps_len else NULL
    try:
        memset(block, 0, pa_len + edits_len + wedges_len + sg_len + maps_len)
        at = PACH3_HEADER_LEN
        for i in range(n):
            dst = pa + i
            # ONE RULE FOR THE ELEMENT BYTE: 1..118 with bit 7 clear is the atomic number, 0x80|0..99
            # is the R index, and everything else is read as a bare R marker -- element 0 is the one
            # code the arena has that is not a claim about which element this is.
            if data[at] & 0x80:
                r_index = data[at] & 0x7f
                if r_index > R_INDEX_MAX:
                    problems.append('atom %d states R index %d and the maximum is %d; read as a bare '
                                    'R' % (i, r_index, R_INDEX_MAX))
                    r_index = 0
                dst.r_index = <uint8_t> r_index
            elif data[at] == 0 or data[at] > 118:
                problems.append('atom %d\'s element byte is %d, which names neither an element nor an '
                                'R index; read as a bare R marker' % (i, data[at]))
            else:
                dst.element = data[at]
            iso_field = data[at + 1] & 0x3f
            if iso_field:
                if dst.element == 0:
                    problems.append('atom %d is an R marker carrying isotope field %d; ignored'
                                    % (i, iso_field))
                else:
                    dst.isotope = <uint16_t> (<int32_t> MDL_ISOTOPE[dst.element]
                                              - PACH3_ISOTOPE_BIAS + <int32_t> iso_field)
            dst.radical = 1 if data[at + 1] & 0x40 else 0
            dst.pinned = 1 if data[at + 1] & 0x80 else 0
            dst.hydrogens = data[at + 2] & 0x0f
            chg = data[at + 2] >> 4
            if <int32_t> chg - 4 > CHARGE_MAX:
                problems.append('atom %d states charge %d and the arena holds %d to %d; clamped'
                                % (i, <int32_t> chg - 4, CHARGE_MIN, CHARGE_MAX))
                dst.charge = CHARGE_MAX
            else:
                dst.charge = <int8_t> (<int32_t> chg - 4)
            if want_xy:
                dst.x = _p3_i24(data + at + 3)
                dst.y = _p3_i24(data + at + 6)
            at += stride
        # ONE CLIP RULE FOR EVERY COUNTED BLOCK.  Each block starts where the DECLARED counts put it and
        # reads its own `*_have`, so a shortfall costs the records that did not arrive and leaves the
        # offsets of every block behind it where the header states them.  `_pach3_length` answers from
        # the same arithmetic.  A `*_declared` is what the header says and nothing reassigns one, which
        # is what these three offsets rest on.
        stereo_at = bonds_at + <Py_ssize_t> nb_declared * PACH3_BOND_LEN
        groups_at = stereo_at + <Py_ssize_t> ns_declared * PACH3_STEREO_LEN
        map_at = groups_at + <Py_ssize_t> nsg_declared * PACH3_SGROUP_LEN
        nm_declared = 0
        if flags & PACH3_FLAG_MAP:
            nm_declared = n                            # the map block states one number per atom
        if nb_have < nb_declared:
            problems.append('the header declares %d bond(s) and the buffer holds %d; the rest of the '
                            'block was not read' % (nb_declared, nb_have))
        nb = 0
        at = bonds_at
        for k in range(nb_have):
            a1 = _p3_u16(data + at)
            a2 = _p3_u16(data + at + 2)
            order = data[at + 4] & 0x0f
            wedge = data[at + 4] >> 4
            at += PACH3_BOND_LEN
            if a1 >= n or a2 >= n:
                problems.append('bond %d names atom index %d and this record has %d atom(s); the '
                                'bond is dropped' % (k, a2 if a2 >= n else a1, n))
                continue
            if a1 == a2:
                problems.append('bond %d joins atom %d to itself; the bond is dropped' % (k, a1))
                continue
            if _pach3_seen(edits, nb, a1, a2):
                problems.append('bond %d repeats the pair (%d, %d); the repeat is dropped'
                                % (k, a1, a2))
                continue
            # An unreadable order is still a connection, and the arena has no "order unknown" the way
            # it has H_UNKNOWN, so the bond is stored single and the substitution is reported.
            if order != 1 and order != 2 and order != 3 and order != 4 and order != 8:
                problems.append('bond %d states order %d, which is not one of 1, 2, 3, 4 and 8; read '
                                'as single' % (k, order))
                order = 1
            if wedge:
                if not want_xy:
                    problems.append('bond %d carries wedge %d in a version 4 record, which has no '
                                    'coordinates for it to mean anything against; read as none'
                                    % (k, wedge))
                    wedge = 0
                elif wedge > 3:
                    problems.append('bond %d states wedge %d and the codes are 0 to 3; read as none'
                                    % (k, wedge))
                    wedge = 0
            e = edits + nb
            e.src = a1
            e.dst = a2
            e.order = <uint8_t> order
            wedges[nb] = <uint8_t> wedge
            nb += 1
        # THE GROUP AND MAP ROWS ARE READ BEFORE THE BUILD AND THE STEREO BLOCK APPLIED AFTER IT: the
        # persistent segments are laid out once, so the build has to be told that a parity segment is
        # wanted and handed both finished rows.
        ns_have = _p3_present(stereo_at, length, ns_declared, PACH3_STEREO_LEN)
        if ns_have < ns_declared:
            problems.append('the header declares %d stereo configuration(s) and the buffer holds %d; '
                            'the rest of the block was not read' % (ns_declared, ns_have))
        nsg_have = _p3_present(groups_at, length, nsg_declared, PACH3_SGROUP_LEN)
        if nsg_have < nsg_declared:
            problems.append('the header declares %d enhanced-stereo entr(ies) and the buffer holds %d; '
                            'the rest of the block was not read' % (nsg_declared, nsg_have))
        at = groups_at
        for k in range(nsg_have):
            slot = _p3_u16(data + at)
            value = data[at + 2]
            at += PACH3_SGROUP_LEN
            if slot >= n:
                problems.append('enhanced-stereo entry %d names atom index %d and this record has '
                                '%d atom(s); the entry was dropped' % (k, slot, n))
                continue
            kind = sg_kind(<uint8_t> value)
            group = sg_group(<uint8_t> value)
            # A ZERO BYTE STATES THE DEFAULT, which is what an absent entry states, so it is not
            # damage on its own -- only a group index with no kind to own it is.
            if kind == 0:
                if group:
                    problems.append('enhanced-stereo entry %d states no kind and group index %d; '
                                    'the entry was dropped' % (k, group))
                continue
            if kind == 1:
                if group:
                    problems.append('enhanced-stereo entry %d is abs and carries group index %d, '
                                    'which only or and and take; read as abs with none'
                                    % (k, group))
                    value = sg_pack(1, 0)
            elif group == 0:
                problems.append('enhanced-stereo entry %d is %s and states no group index, and 1 '
                                'to 63 is what one takes; the entry was dropped'
                                % (k, 'or' if kind == 2 else 'and'))
                continue
            if sg[slot]:
                problems.append('enhanced-stereo entry %d repeats atom index %d; the repeat was '
                                'dropped' % (k, slot))
                continue
            sg[slot] = <uint8_t> value
            sg_any = True
        nm_have = _p3_present(map_at, length, nm_declared, PACH3_MAP_LEN)
        if nm_have < nm_declared:
            problems.append('the header declares %d map number(s) and the buffer holds %d; the rest of '
                            'the block was not read' % (nm_declared, nm_have))
        at = map_at
        for i in range(nm_have):
            value = _p3_u16(data + at)
            at += PACH3_MAP_LEN
            if value > MAP_NUMBER_MAX:
                # THE ONLY GUARD ON THIS PATH.  A pach build allocates its own arena, so the
                # `structure_from_bytes` validator never sees these atoms.
                problems.append('atom %d states map number %d and the arena holds 0 to %d; read as none'
                                % (i, value, MAP_NUMBER_MAX))
                continue
            maps[i] = <uint16_t> value
            map_any = True
        try:
            mol = _pach3_build(pa, n, edits, nb, want_xy, _p3_u16(data + 6) != 0,
                               sg if sg_any else NULL, maps if map_any else NULL)
        except ValueError as err:
            _pach_derivation_lost(problems, err)
            return (None, problems)
        _pach3_apply_wedges(mol, edits, wedges, nb)
        _pach3_apply_stereo(mol, data, stereo_at, ns_have, n, problems)
        # RETURNED INSIDE THE `try`, which the `finally` below still covers: the build's own refusal
        # arm leaves `mol` unassigned, and a return after the `finally` would read a name that path
        # never wrote.
        return (mol, problems)
    finally:
        PyMem_Free(block)


def pach_load(data, *, compressed=None):
    """Read one pach record.  `(MoleculeContainer or None, problems)`.  MemoryError is the one
    exception that escapes, and it is not a statement about the bytes.  A graph the derivation refuses
    -- `perceive_rings` has a prototype limit and a deadline, and a well-formed record can state a
    graph that reaches either -- is `None` and a sentence, not the raise `edit()`'s seal would give.

    `problems` is a list of sentences about what the record got wrong.  A molecule and a non-empty
    list together is the normal outcome for a damaged record: the decoder stores what it can read and
    says what it could not, because a loop over forty thousand stored records must not be stopped by
    one of them.  `None` means nothing at all could be built and the list says why.

    `compressed` defaults to sniffing, and the sniff is exact rather than heuristic: a raw record's
    first byte is its version, one of 0, 2, 3 and 4, and a zlib header's low nibble is its compression
    method, always 8, so neither value can be the other.  `True` and `False` state it instead, which a
    caller who would rather hear that its store is not what it thought can pass.

    Trailing bytes after the record are ignored, so a caller walking a concatenated stream can hand the
    rest of the buffer over and use `pach_record_length` to advance.
    """
    cdef list problems = []
    cdef const unsigned char[::1] view
    cdef bytes raw = bytes(data)
    cdef bint looks_raw = len(raw) > 0 and raw[0] in (0, 2, 3, 4)
    if not len(raw):
        problems.append('the buffer is empty; a pach record is at least a 4 byte header')
        return (None, problems)
    if compressed is True and looks_raw:
        problems.append('compressed=True was stated and the buffer begins with a pach version byte, '
                        'so it is a raw record')
        return (None, problems)
    # `compressed=False` hands the bytes straight to the decoder even when they do not look like a
    # record, because the decoder's own report of WHAT is wrong with them is more use to a caller than
    # this function's report that they do not begin with a known version byte.
    if compressed is not False and not looks_raw:
        raw = _pach_decompress(raw, problems)
        if raw is None:
            return (None, problems)
    if not len(raw):
        problems.append('the buffer is empty; a pach record is at least a 4 byte header')
        return (None, problems)
    view = raw
    if view[0] == PACH3_VERSION_XY or view[0] == PACH3_VERSION_FLAT:
        return _pach3_decode(&view[0], view.shape[0])
    return _pach_decode(&view[0], view.shape[0])


def pach_dump(MoleculeContainer mol not None, *, bint compressed=True, drop=None, version=None):
    """Write one pach record.  `version` selects the layout.

    `None` is version 3 when the molecule has coordinates and version 4 when it does not, which is the
    whole of the choice a caller normally makes.  `2` writes the legacy record, which loses
    atropisomers, every stereo group, every wedge and every map number and refuses rather than losing
    them quietly.  `3` and `4` state the third-generation layout outright, and `4` writes no coordinate
    block even for a molecule that has one.

    AN EXPLICIT VERSION IS A REQUEST THE ENCODER ANSWERS, and two things answer `3` with version 4: a
    molecule with no drawing, because a coordinate block of zeros would state a position nothing
    recorded, and a caller who passed `drop=['coordinates']`, because a waiver named at the door wins
    over the version asked for beside it.

    Raises `ValueError` naming any field the arena holds and the chosen version cannot carry.  `drop`
    waives those refusals: an iterable of field names, or `'*'` for all of them.  The names are
    `map_number`, `title`, `sgroups`, `cip`, `wedges`, `stereo_groups`, `stereo`, `meta`,
    `coordinates` and `conformers`; `conformers` is asked only by versions 3 and 4.  An unrecognised
    name is refused rather than ignored.
    """
    cdef uint32_t mask = 0
    cdef object name
    if drop is None:
        pass
    elif drop == '*':
        mask = PACH_DROP_ALL
    else:
        for name in drop:
            if name not in _PACH_DROP_NAMES:
                raise ValueError('%r is not a droppable field; the drop names are %s'
                                 % (name, ', '.join(sorted(_PACH_DROP_NAMES))))
            mask |= <uint32_t> <int> _PACH_DROP_NAMES[name]
    cdef bytes raw
    # `None` AND `3` ARE ONE ARM.  Both ask for the coordinates the molecule has, and `drop=` is the
    # waiver that wins tree-wide: stripping `PACH_DROP_COORDINATES` back out of the caller's mask here
    # would validate the name and then discard it, and the record would come back version 3 with the
    # drawing the caller waived.
    if version is None or version == PACH3_VERSION_XY:
        raw = _pach3_encode(mol, mask)
    elif version == 2:
        raw = _pach_encode(mol, mask & 0xff)
    elif version == PACH3_VERSION_FLAT:
        raw = _pach3_encode(mol, mask | PACH_DROP_COORDINATES)
    else:
        raise ValueError('%r is not a writable pach version; they are 2, 3, 4 and None for 3-or-4 '
                         'by whether the molecule has coordinates' % (version,))
    if compressed:
        return zlib.compress(raw, 9)
    return raw


def pach_record_length(data, *, compressed=None):
    """How many bytes the first pach record in `data` occupies.

    The length is not stored: it is a function of the header counts and the per-version strides,
    which is why a caller walking a stream of concatenated records needs this rather than arithmetic
    of its own.  Raises `ValueError` -- it returns a number and has no way to say "unknown", so it
    is an answer boundary like `unpack` and not a loop-safe door like `pach_load`.
    """
    cdef list problems = []
    cdef bytes raw = bytes(data)
    cdef const unsigned char[::1] view
    cdef uint32_t atoms_count, ct_count, i, deg_sum = 0
    cdef unsigned char version
    # The same exact sniff as `pach_load`: the known version bytes are 0, 2, 3 and 4, and a zlib
    # header's low nibble is its compression method, always 8, so neither value can be the other.
    if compressed is not False and not (len(raw) and raw[0] in (0, 2, 3, 4)):
        raw = _pach_decompress(raw, problems)
        if raw is None:
            raise ValueError(problems[0])
    view = raw
    if view.shape[0] < 4:
        raise ValueError('a pach record is at least a 4 byte header and this buffer is %d byte(s)'
                         % view.shape[0])
    version = view[0]
    if version == PACH3_VERSION_XY or version == PACH3_VERSION_FLAT:
        if view.shape[0] < PACH3_HEADER_LEN:
            raise ValueError('a version %d pach record is at least a 12 byte header and this buffer '
                             'is %d byte(s)' % (version, view.shape[0]))
        return _pach3_length(&view[0], view.shape[0])
    if version != 0 and version != 2:
        raise ValueError('byte 0 is %d, which is not a pach version; the molecule versions are 0, 2, '
                         '3 and 4' % version)
    atoms_count = (view[1] << 4) | (view[2] >> 4)
    ct_count = ((view[2] & 0x0f) << 8) | view[3]
    if 4 + <Py_ssize_t> 9 * atoms_count > view.shape[0]:
        raise ValueError('the header declares %d atoms and the buffer is only %d bytes, so the '
                         'record\'s length cannot be computed' % (atoms_count, view.shape[0]))
    for i in range(atoms_count):
        deg_sum += view[4 + 9 * <Py_ssize_t> i + 1] & 0x0f
    return (4 + 9 * <Py_ssize_t> atoms_count + 3 * <Py_ssize_t> (deg_sum // 2)
            + _pach_order_block_len(deg_sum // 2, version) + 4 * <Py_ssize_t> ct_count)
