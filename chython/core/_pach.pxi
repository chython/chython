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
# THE PACH FORMAT: chython's stored-structure format from 1.1 to 2.24, read and written here.
# ==================================================================================================
#
# WHY THIS FRAGMENT EXISTS.  Stored data outlives the code that wrote it.  There are databases whose
# keys are `MoleculeContainer.pack()` output and nothing else, and a V3 that cannot read them is a V3
# nobody can migrate to.  THE SPECIFICATION is chython 2's own writer and reader, read out of git --
# `git show 5e39eb5:chython/containers/_pack_v2.pyx` and `:chython/containers/_unpack_v0v2.pyx` -- plus
# the three corpora in `core/test/`.
#
# TWO VERSIONS EXIST AND BOTH ARE READ.  Byte 0 is the version.
#
#   version 0   chython 1.1 through 1.44.  Not an inherited format: chython's own first one.
#   version 2   chython 1.45 through 2.24.
#
# The two differ in ONE block -- how bond orders are packed -- and are byte-identical everywhere else.
# There is no version 1 in any release.
#
# THE LAYOUT, big-endian, bit-packed, no alignment anywhere:
#
#   header, 4 bytes      byte 0        format version
#                        12 bits       atom count            (data[1] << 4 | data[2] >> 4)
#                        12 bits       cis/trans entry count ((data[2] & 15) << 8 | data[3])
#
#   atom block, 9 bytes per atom, in the writer's own atom order:
#                        12 bits       atom number, 1..4095 -- the stable id, NOT an index
#                        4 bits        neighbour count, 0..15
#                        4 bits        stereo nibble (below)
#                        5 bits        isotope shift (0 = unset; see PACH_ISOTOPE_BIAS)
#                        7 bits        atomic number, 1..118
#                        16 + 16 bits  x, y as float16
#                        3 bits        implicit hydrogen count, 7 = "not stated"
#                        4 bits        formal charge + 4
#                        1 bit         radical flag
#
#   connection table     the atoms' neighbour lists, concatenated in atom-block order, each entry a
#                        12-bit atom NUMBER, two entries to three bytes.  Length 3 * bond count, and
#                        the bond count is half the sum of the neighbour counts -- the table is
#                        bidirected and carries every bond twice.
#
#   bond order block     one 3-bit value per bond, `order - 1`, in the order the connection table is
#                        CONSUMED: a bond's order appears where the table first names it, which is at
#                        the lower-positioned of its two atoms.  version 2 packs the values as a flat
#                        MSB-first bitstream, ceil(3 * bonds / 8) bytes; version 0 packs five values
#                        into two bytes with one pad bit at the TOP, ceil(bonds / 5) * 2 bytes.
#
#   cis/trans block      4 bytes per entry: two 12-bit atom numbers naming the two TERMINALS of a
#                        cumulene chain, 7 pad bits, 1 sign bit.
#
# THE STEREO NIBBLE is two 2-bit fields, `tetrahedron | allene`, and the writer chose between them on
# the atom's neighbour COUNT: an atom with exactly two neighbours got the allene field (0b0010 /
# 0b0011), anything else the tetrahedron field (0b1000 / 0b1100).  The second bit is the sign.  The
# 2.24 reader collapses all four with a catch-all `else: True`, so the distinction the writer made is
# not recoverable from its own answers, and this decoder does not pretend otherwise: the sign is read
# and WHICH UNIT IT BELONGS TO IS DECIDED BY THE GRAPH.  That is not a loss -- an allene centre has
# two neighbours and a tetrahedral centre has three or four, so the graph knows.
#
# ==================================================================================================
# `pack` AND `to_bytes` ARE TWO FORMATS ON PURPOSE, AND HERE IS THE DIVISION.
# ==================================================================================================
#
#   `pack`/`unpack`  (this fragment)  THE LEGACY WIRE FORMAT.  Its only job is compatibility with
#       bytes chython 1.x and 2.x wrote and with readers that expect them.  It is small -- 9 bytes an
#       atom, and zlib on top -- and it is LOSSY: see the losses list below.  Nothing new should be
#       stored in it, and nothing in it may ever change, because a format whose spelling drifts is
#       not a format a stored key can be read with.
#
#   `to_bytes`/`from_bytes`  (`_molecule_arena.pxi`)  THE ARENA, verbatim.  It is the molecule's own
#       memory, so it is lossless by construction -- S-groups, CIP codes, wedges, enhanced stereo
#       groups, map numbers, the title, 32-bit atom counts, absolute isotopes -- and it is the
#       identity `__reduce__` pickles through.  It is bigger per atom and compresses less well.
#
# So: `to_bytes` for anything this release writes, `unpack` for anything an earlier one wrote.  The
# two are told apart WITHOUT A FLAG, by the first byte: a pach record's is its version, 0 or 2, and
# the arena's is part of `CHY3`.  `MoleculeContainer.unpack` dispatches on it, which is why a caller
# holding a column of stored keys of unknown vintage can feed all of them to one door.
#
# ==================================================================================================
# WHAT THE ARENA HOLDS AND PACH CANNOT.  Every one of these makes the WRITER REFUSE, by name, unless
# the caller passes `drop=`.  A serialiser that quietly throws data away is worse than one that
# cannot serialise, because the loss is discovered by whoever reads the record back years later.
# ==================================================================================================
#
#   map_number      no field.  Atom-to-atom mapping is the point of a reaction record.
#   title           no field, and pach has no text of any kind.
#   meta            no field: record metadata, an SDF data field or an RDfile DTYPE/DATUM pair.
#   sgroups         no field: superatoms, data S-groups, polymer brackets, their labels.
#   cip             no field: atom and bond CIP descriptors.
#   wedges          no field: which bond was drawn as a wedge, and from which end.
#   stereo_groups   no field: enhanced stereo (AND/OR/ABS) group membership.
#   stereo          only three unit kinds have a slot -- tetrahedral, allene, cis/trans.  An
#                   atropisomer's or a helical unit's parity has nowhere to go.
#
# THE ONE FIELD THAT IS DEGRADED RATHER THAN REFUSED IS THE COORDINATE, and it is the exception because
# refusing it would refuse every drawn molecule there is.  pach carries x and y as float16 -- about
# three significant decimal digits, |x| < 65504 -- and NOT AS A PRESENCE FLAG: it writes 0 for an atom
# with no coordinate, so "this record has no drawing" and "every atom sits at the origin" are the same
# four bytes.  The decoder resolves that in favour of the first reading whenever every coordinate is
# zero BYTE, because a molecule read from a string and packed is the common case and it has no drawing
# all.  Both halves of this are losses; neither is `drop`-able, and both are stated here instead.
#
# And these the writer refuses outright, with no `drop` spelling, because the record would be a lie
# rather than a subset: an atom number above 4095, an atom with more than 15 neighbours, an implicit
# hydrogen count above 6 (the field is 3 bits with 7 taken by the sentinel), an isotope more than 15
# mass units from its element's MDL reference, a formal charge outside -4..+8.
#
# TWO THINGS THE FORMAT CARRIES THAT THE ARENA CANNOT HOLD EXACTLY, reported rather than hidden:
#
#   * a float16 coordinate.  The arena stores a display coordinate as a x10000 fixed-point int32, so
#     the float16 `1.2001953125` comes back as `1.2002`.  Half of 1e-4 is the worst error.
#   * an implicit hydrogen count of 7..14 in a v0/v2 record cannot occur (the field caps at 6), but
#     the 3-bit field's value 7 is the "not stated" sentinel and the arena's is 15 -- the two
#     spellings are different and neither is a count.  `implicit_h_of` answers None for both.
#
# ==================================================================================================
# INPUT BY DEFAULT IS GARBAGE, AND THE DECODER NEVER RAISES.
# ==================================================================================================
#
# `pach_load` returns `(molecule_or_None, problems)` where `problems` is a list of sentences.  A
# truncated, bit-flipped or version-mismatched record must not take down a loop over forty thousand
# of them, so nothing in the decode path raises: it stores what it can and says what was wrong.
# `MoleculeContainer.unpack` is the ANSWER BOUNDARY and raises `ValueError` when the decoder could
# not build a molecule at all -- a caller who asked for a molecule and cannot have one is told so.
# Both doors, one decoder.
#
# NO REPAIR RUNS HERE, IN EITHER DIRECTION.  A record written aromatic is stored aromatic and a
# record written Kekule is stored Kekule; the decoder does not aromatise, the writer does not
# kekulise, and neither calls the standardiser.  IO is not a mutator of representation, and a writer
# that repairs its input is the same violation as a reader that normalises.
#
# ==================================================================================================
# WHAT THIS DECODER READS AROUND.  Each is a statement about bytes that already exist in stored data,
# and none is reproduced by the V3 writer.  Every one has a test in `core/test/test_pach.py`.
# ==================================================================================================
#
#   * the v0 order block's first value is read `a >> 4`, four bits, with no mask.  A set pad bit --
#     bit 7, which the format does not define -- therefore yields an order of 9..16, and 2.24 built a
#     Bond with it.  Masked to three bits here.
#   * the connection table is trusted to be symmetric: an entry whose partner does not name it back
#     raises KeyError out of the C extension.  Dropped with a note here.
#   * the implicit hydrogen count is written `<unsigned char> h << 5`, so a count of 8 becomes 0 and
#     a count of 7 becomes the "not stated" sentinel.  The V3 writer refuses instead.
#   * the charge field holds `charge + 4` in four bits and admits +5..+11, which no chython 2 Element
#     could hold.  Read, reported and clamped here.
#   * the isotope shift is written `isotope - common_isotope` with no range check, so an exotic
#     isotope wraps into a neighbouring element's mass or into the "unset" spelling.  Reported here.
#   * `pack(check=True)` refuses a molecule with no bonds, although the format has no such
#     restriction -- the order block is simply zero bytes long, so a lone sodium cation is written here.
#   * the cis/trans sign is written into bit 0 of a byte whose other seven bits are documented as
#     padding, and read back as `if d:` -- the whole byte.  Read the same way here, deliberately: a
#     record whose pad bits are set was read as a sign by the writer's own reader, and a decoder that
#     disagreed with it would report a different molecule than the one that was stored.
# ==================================================================================================

# The 12-bit atom number field, and the reason a V3 molecule is not always packable: arena stable ids
# are 32-bit and never reused, so a molecule that has been edited enough has ids this cannot hold.
DEF PACH_MAX_NUMBER = 4095
DEF PACH_MAX_DEGREE = 15
# 3 bits with 7 spoken for.  The arena's own cap is 14 and its sentinel is 15; neither maps.
DEF PACH_MAX_IMPLICIT_H = 6
DEF PACH_H_UNKNOWN = 7
# The 5-bit isotope field spells `MDL_ISOTOPE[z] - 16 + shift` for shift 1..31, and 0 for "unset" --
# so mass numbers from 15 below the element's MDL reference to 15 above it.  A CONSTANT AND NOT A
# TABLE: chython 2's `common_isotopes[z]` is `MDL_ISOTOPE[z] - 16` for every z in 1..118, verified
# element by element against `_elements.pxi`, so one isotope table serves both.
DEF PACH_ISOTOPE_BIAS = 16

# WHICH ARENA PARITY chython 2's SIGN BIT MEANS, per unit kind.  A sign is not a fact about an atom,
# it is a fact about an atom RELATIVE TO AN ORDER of its neighbours, and the two codebases order
# neighbours differently (chython 2: `_bonds` insertion order with hydrogens moved last; the arena:
# ruling F26).  `_pach_*_frame` below rebuilds chython 2's order as arena slots and `translate_parity`
# re-expresses the parity in the arena's -- that part is arithmetic and settles itself.
#
# What arithmetic cannot settle is the one bit of CONVENTION left over: whether chython 2's True is
# the even parity or the odd one in its own frame.  Both codebases' sign algebra is the parity of a
# permutation (chython 2's `_tetrahedron_translate` and `_alkene_translate` are permutation-parity
# tables and nothing else -- checked entry by entry), so the leftover is exactly one bit per unit
# kind, and it was measured rather than reasoned: `test_every_decoded_parity_says_what_the_smiles_
# reader_says` re-reads each corpus record's own source string with the V3 SMILES reader and compares
# atom by atom, over both parities, on hundreds of centres.  These three constants are that
# measurement.  They are the parity a sign of TRUE means; False means the other one.
# The measurement, on 49 tetrahedral, 16 cis/trans and 3 allene centres: with all three set to 2 every
# tetrahedral and every allene centre agreed and EVERY cis/trans centre disagreed -- all of them, not
# most, which is the signature of a convention bit rather than of a broken frame.  So cis/trans is 1
# and the other two are 2.  chython 2 spelled a tetrahedron's sign against a frame of four neighbours
# of one atom and a cumulene's against two neighbours of each end, and nothing ever required the two
# spellings to agree about which arrangement is the even one.
#
# THE ALLENE BIT MOVES WITH THE SMILES READER'S AXIAL FRAME, since the reader is what it is measured
# against: all three corpus allenes state their configuration on a chain end whose second direction is
# an implicit hydrogen, so the position that reader gives an unwritten hydrogen (`smi_written_pair`)
# decides this constant.  CDK 2.12 reads all three the way the current reader does.
DEF PACH_TETRA_TRUE = 2
DEF PACH_ALLENE_TRUE = 2
DEF PACH_CIS_TRANS_TRUE = 1


cdef object zlib
import zlib


cdef struct pach_atom_t:
    uint16_t number
    uint8_t degree
    uint8_t element
    int8_t charge
    uint8_t radical
    uint8_t hydrogens          # 0..6, or PACH_H_UNKNOWN
    uint16_t isotope           # absolute mass number, 0 for unset
    uint8_t sign               # 0 none, 1 False, 2 True
    int32_t x                  # already in the arena's x10000 fixed point
    int32_t y


# --------------------------------------------------------------------------------------------------
# float16.  Two functions, deliberately not each other's inverse.
# --------------------------------------------------------------------------------------------------

cdef inline double _pach_scale2(double x, int e) noexcept nogil:
    """`x * 2 ** e` by repeated doubling. Exact in binary floating point over this fragment's range,
    and it costs no libc dependency the rest of the core does not already have."""
    cdef int k
    if e > 0:
        for k in range(e):
            x *= 2.
    else:
        for k in range(-e):
            x *= .5
    return x


cdef inline double _pach_f16_decode(unsigned char a, unsigned char b) noexcept nogil:
    """chython 2's `double_from_bytes`, value for value.

    An exponent field of 31 is infinity or NaN in IEEE 754 and this reads it as an ordinary exponent
    of 16, so a record carrying 0x7c00 decodes to 131072.0 rather than inf.  That is what chython 2
    answered for those bytes, no stored record written by its own writer can contain them (the writer
    refuses anything with |x| >= 2**16 and writes a hard zero instead), and reproducing it keeps a
    bit-flipped record decoding to the molecule the previous reader saw rather than to a NaN that
    would poison every coordinate arithmetic downstream.
    """
    cdef int e = (a >> 2) & 0x1f
    cdef double x = (((a & 0x03) << 8) | b) / 1024.
    if e:
        x += 1.
        e -= 15
    else:
        e = -14
    x = _pach_scale2(x, e)
    if a >> 7:
        return -x
    return x


cdef inline void _pach_f16_encode(double x, unsigned char *p) noexcept nogil:
    """The NEAREST float16 to `x`, big-endian, or a hard zero when the value has no float16 at all.

    chython 2 TRUNCATES here (`bits = <unsigned short> f | ...`), losing up to a full unit in the last
    place downwards on every coordinate whose mantissa did not fit.  This ROUNDS, and what compatibility
    demands is only that chython 2 can read what this writes, which it can, because a correctly rounded
    float16 is a float16.  The price is that a decode/encode round trip of a coordinate-bearing record
    is not byte-identical, which is measured and reported rather than papered over.

    The zero cases are chython 2's, and they are the reason there is no infinity in a pach record:
    anything that would overflow the exponent field is written as 0 rather than as inf.  Ties round
    away from zero rather than to even; the field has 3 decimal digits and a tie is a coordinate whose
    last bit was never information.
    """
    cdef uint32_t bits
    cdef uint32_t sign = 0
    cdef int e
    cdef double m
    if x != x or x == 0.:                       # NaN has no spelling here either
        p[0] = 0
        p[1] = 0
        return
    if x < 0.:
        sign = 0x8000
        x = -x
    if x >= 65520.:                             # rounds to the exponent field's overflow
        p[0] = 0
        p[1] = 0
        return
    if x < 6.103515625e-05:                     # 2 ** -14: subnormal, or nothing
        # The IEEE encoding is monotone in the integer, so a subnormal that rounds up to 1024 becomes
        # the smallest normal with the right bits by arithmetic rather than by a branch.
        bits = <uint32_t> round(_pach_scale2(x, 24))
        p[0] = <unsigned char> ((sign | bits) >> 8)
        p[1] = <unsigned char> (sign | bits)
        return
    # Normalise INTO [1, 2) in both directions.  `x` is already known to be at least 2 ** -14, so the
    # downward loop cannot run past the smallest normal exponent.
    e = 0
    m = x
    while m >= 2.:
        m *= .5
        e += 1
    while m < 1.:
        m *= 2.
        e -= 1
    # Monotone again: a mantissa that rounds up to 2048 carries into the exponent field correctly.
    bits = <uint32_t> ((e + 15) << 10) + <uint32_t> round(_pach_scale2(m, 10)) - 1024u
    p[0] = <unsigned char> ((sign | bits) >> 8)
    p[1] = <unsigned char> (sign | bits)


# --------------------------------------------------------------------------------------------------
# Reading the blocks.  Every accessor takes the buffer's length and reports "not there" rather than
# reading past it, which is what makes the decoder total over arbitrary bytes.
# --------------------------------------------------------------------------------------------------

cdef inline int _pach_number_at(const unsigned char *data, Py_ssize_t length, Py_ssize_t base,
                                uint32_t k) noexcept nogil:
    """The k-th 12-bit number of a table starting at `base`, or -1 when the buffer ends first."""
    cdef Py_ssize_t at = base + 3 * (k >> 1)
    if k & 1:
        if at + 3 > length:
            return -1
        return ((data[at + 1] & 0x0f) << 8) | data[at + 2]
    if at + 2 > length:
        return -1
    return (data[at] << 4) | (data[at + 1] >> 4)


cdef inline int _pach_order_at(const unsigned char *data, Py_ssize_t length, Py_ssize_t base,
                               uint32_t k, unsigned char version) noexcept nogil:
    """The k-th 3-bit bond order value, or -1 when the buffer ends first.

    Version 2's block is a flat MSB-first bitstream of 3-bit values; chython 2 unpacks it with an
    eight-state machine over the `3 3 2 | 1 3 3 1 | 2 3 3` pattern, which is the same thing written
    out longhand.  Version 0's is five values to a 16-bit group with ONE PAD BIT AT THE TOP, so the
    j-th value of a group is `word >> (12 - 3 * j)`; chython 2 reads the first one as `a >> 4` with no
    mask, which lets that pad bit -- a bit the format does not define -- become part of an order.
    """
    cdef Py_ssize_t at
    cdef uint32_t bit, group, word
    if version == 2:
        bit = 3 * k
        at = base + (bit >> 3)
        bit &= 7
        if bit <= 5:
            if at + 1 > length:
                return -1
            return (data[at] >> (5 - bit)) & 7
        if at + 2 > length:
            return -1
        return (((data[at] << 8) | data[at + 1]) >> (13 - bit)) & 7
    group = k // 5
    at = base + 2 * <Py_ssize_t> group
    if at + 2 > length:
        return -1
    word = (data[at] << 8) | data[at + 1]
    return (word >> (12 - 3 * (k - 5 * group))) & 7


cdef inline Py_ssize_t _pach_order_block_len(uint32_t bonds, unsigned char version) noexcept nogil:
    """Bytes the order block occupies. Zero bonds is zero bytes in both versions."""
    if version == 2:
        return (3 * <Py_ssize_t> bonds + 7) // 8
    return ((<Py_ssize_t> bonds + 4) // 5) * 2


# --------------------------------------------------------------------------------------------------
# THE DECODER, constitution and coordinates.  One function, because every block's interpretation
# depends on the ones before it -- the bond count is a function of the atom block, the order block's
# position is a function of the bond count -- and splitting it would mean handing the whole parse state
# across the seam.
# --------------------------------------------------------------------------------------------------

cdef int _pach_derivation_lost(list problems, object err) except -1:
    """A record whose graph the derivation refuses, reported rather than raised.

    Every builder ends in `rebuild_derived`, which runs `perceive_rings` -- an answer boundary with a
    relevant-cycle prototype limit and a wall-clock deadline (`_rings.pxi:_raise_rc`).  A well-formed
    record can state a graph that trips either: a 150-atom complete graph is 56 KB and reaches the
    prototype limit.  The decode path is not an answer boundary, so `pach_load` and
    `reaction_pach_load` answer `None` and this sentence where `edit()`'s seal raises.

    Shared by both decoders, so the four call sites cannot word one loss two ways.
    """
    problems.append('this record\'s graph could not be derived (%s), so no molecule was read' % err)
    return 0


cdef tuple _pach_decode(const unsigned char *data, Py_ssize_t length):
    """`(MoleculeContainer or None, problems)` for one pach record.

    NOTHING HERE RAISES ON THE RECORD'S CONTENT.  A truncated, bit-flipped or forged record produces
    a molecule holding what could be read and a list of sentences saying what could not, or `None` and
    a list saying why nothing at all could be built.  MemoryError is the one exception that escapes,
    and it is not a statement about the bytes.  The graph derivation's own refusals are caught at the
    build -- `_pach_derivation_lost` -- because they are a statement about the bytes.
    """
    cdef list problems = []
    cdef unsigned char version, nibble, hcr
    cdef uint32_t atoms_count, declared, ct_count, i, j, k, n_edits = 0, bonds_count, deg_sum = 0
    cdef uint32_t other, pos
    cdef int value, order, mass
    cdef Py_ssize_t table_at, order_at, ct_at, atom_at
    cdef bint want_xy = False
    cdef bint want_parity = False
    cdef bint short_table = False
    cdef bint short_orders = False
    cdef bint dup, symmetric
    cdef pach_atom_t *pa = NULL
    cdef int32_t *slot_of = NULL
    cdef uint32_t *nbr = NULL
    cdef uint32_t *off = NULL
    cdef uint8_t *seen = NULL
    cdef edge_edit_t *edits = NULL
    cdef MoleculeContainer mol = None
    cdef object err

    if length < 4:
        problems.append('a pach record is at least a 4 byte header and this buffer is %d byte(s)'
                        % length)
        return (None, problems)
    version = data[0]
    if version != 0 and version != 2:
        problems.append('byte 0 is %d, which is not a pach version; only 0 (chython 1.1 to 1.44) and '
                        '2 (1.45 to 2.24) were ever written' % version)
        return (None, problems)
    atoms_count = (data[1] << 4) | (data[2] >> 4)
    ct_count = ((data[2] & 0x0f) << 8) | data[3]

    if atoms_count == 0:
        # Not an error.  chython 2 packed an empty molecule as a bare header and an empty molecule is
        # a molecule; refusing it here would make a legitimate stored record unreadable.
        if ct_count:
            problems.append('the header declares %d cis/trans entries for a record with no atoms; '
                            'they can name nothing and were dropped' % ct_count)
        try:
            return (_pach_build(NULL, 0, NULL, 0, False, False), problems)
        except ValueError as err:
            _pach_derivation_lost(problems, err)
            return (None, problems)

    pa = <pach_atom_t *> PyMem_Malloc(atoms_count * sizeof(pach_atom_t))
    slot_of = <int32_t *> PyMem_Malloc(4096 * sizeof(int32_t))
    off = <uint32_t *> PyMem_Malloc((atoms_count + 1) * sizeof(uint32_t))
    seen = <uint8_t *> PyMem_Malloc(atoms_count * sizeof(uint8_t))
    if pa is NULL or slot_of is NULL or off is NULL or seen is NULL:
        PyMem_Free(pa)
        PyMem_Free(slot_of)
        PyMem_Free(off)
        PyMem_Free(seen)
        raise MemoryError('pach decode scratch allocation failed')
    try:
        for i in range(4096):
            slot_of[i] = -1
        # PyMem_Malloc does not zero, and `seen` is read for a neighbour the walk has not reached yet:
        # garbage there makes a bond's order be consumed at the wrong end, which silently permutes
        # every order after it.
        memset(seen, 0, atoms_count * sizeof(uint8_t))

        # ---- the atom block.  A truncated one is NOT fatal: the atoms that are wholly present are
        # real atoms and the caller is entitled to them, so the count is lowered and said out loud.
        if 4 + <Py_ssize_t> 9 * atoms_count > length:
            declared = atoms_count
            atoms_count = <uint32_t> ((length - 4) // 9)
            problems.append('the header declares %d atoms and the buffer holds %d whole atom '
                            'block(s); the record is truncated' % (declared, atoms_count))
            if atoms_count == 0:
                return (None, problems)

        for i in range(atoms_count):
            atom_at = 4 + 9 * <Py_ssize_t> i
            pa[i].number = <uint16_t> ((data[atom_at] << 4) | (data[atom_at + 1] >> 4))
            pa[i].degree = data[atom_at + 1] & 0x0f
            deg_sum += pa[i].degree

            # chython 2's own reader: 0 is "no configuration", 0b0010 and 0b1000 are False, and
            # EVERYTHING ELSE is True -- a catch-all that swallows the tetrahedron/allene distinction
            # its writer made in the same nibble.  Read the same way, because that is the molecule
            # the previous reader reported for these bytes.  Which unit the sign belongs to is then
            # decided by the graph, and the graph knows: an allene centre has two neighbours.
            nibble = data[atom_at + 2] >> 4
            if nibble == 0:
                pa[i].sign = 0
            elif nibble == 0b0010 or nibble == 0b1000:
                pa[i].sign = 1
            else:
                pa[i].sign = 2

            pa[i].element = data[atom_at + 3] & 0x7f
            if pa[i].element < 1 or pa[i].element > 118:
                problems.append('atom %d states atomic number %d, which is not an element; nothing '
                                'can be built from this record' % (pa[i].number, pa[i].element))
                return (None, problems)
            if pa[i].number == 0:
                problems.append('atom number 0 appears in the atom block; a stable id is 1..4095 and '
                                '0 is how "no atom" is spelled, so nothing can be built from this '
                                'record')
                return (None, problems)
            if slot_of[pa[i].number] >= 0:
                problems.append('a duplicate atom number, %d, appears in the atom block; the '
                                'connection table cannot say which of the two it means, so nothing '
                                'can be built from this record' % pa[i].number)
                return (None, problems)
            slot_of[pa[i].number] = <int32_t> i

            k = ((data[atom_at + 2] & 0x0f) << 1) | (data[atom_at + 3] >> 7)
            if k:
                mass = <int> MDL_ISOTOPE[pa[i].element] - PACH_ISOTOPE_BIAS + <int> k
                if mass < 1:
                    problems.append('atom %d states isotope shift %d, which is mass number %d for '
                                    'element %d; chython 2 wrote this field with no range check, so '
                                    'the isotope is left unset rather than invented'
                                    % (pa[i].number, k, mass, pa[i].element))
                    pa[i].isotope = 0
                else:
                    pa[i].isotope = <uint16_t> mass
            else:
                pa[i].isotope = 0

            pa[i].x = <int32_t> round(_pach_f16_decode(data[atom_at + 4],
                                                       data[atom_at + 5]) * XY_SCALE)
            pa[i].y = <int32_t> round(_pach_f16_decode(data[atom_at + 6],
                                                       data[atom_at + 7]) * XY_SCALE)
            # PRESENCE IS DECIDED ON THE STORED BYTES AND NOT ON THE SCALED RESULT.  A coordinate of
            # 1e-5 is a real coordinate that the arena's fixed point rounds to zero, and judging
            # presence after the rounding would turn a drawing whose atoms happen to sit very close to
            # the origin into a molecule with no drawing at all -- a second loss on top of the
            # quantisation, and one caused by this decoder rather than by either format.
            if data[atom_at + 4] or data[atom_at + 5] or data[atom_at + 6] or data[atom_at + 7]:
                want_xy = True

            hcr = data[atom_at + 8]
            k = hcr >> 5
            if k == PACH_H_UNKNOWN:
                # The two codebases spell "nobody said" differently -- 7 in three bits there, 15 in
                # four here -- and neither is a count.  Translating the spelling is not normalising
                # the molecule; writing 7 into a field whose 7 means seven hydrogens would be.
                pa[i].hydrogens = H_UNKNOWN
            else:
                pa[i].hydrogens = <uint8_t> k
            value = <int> ((hcr >> 1) & 0x0f) - 4
            if value > CHARGE_MAX:
                problems.append('atom %d states formal charge %+d; the field holds charge+4 in four '
                                'bits so it admits up to +11, which no chython 2 element could hold '
                                'and none ever wrote, and it is clamped to %+d'
                                % (pa[i].number, value, CHARGE_MAX))
                value = CHARGE_MAX
            pa[i].charge = <int8_t> value
            pa[i].radical = hcr & 0x01

        # ---- the connection table.  The bond count is not stored anywhere: it is half the sum of the
        # neighbour counts, which is why an odd sum is a corrupt record rather than a rounding.
        if deg_sum & 1:
            problems.append('the neighbour counts sum to %d, an odd number; a bidirected table '
                            'carries every bond twice, so at least one entry is missing' % deg_sum)
        bonds_count = deg_sum // 2
        table_at = 4 + 9 * <Py_ssize_t> atoms_count
        order_at = table_at + 3 * <Py_ssize_t> bonds_count
        ct_at = order_at + _pach_order_block_len(bonds_count, version)

        off[0] = 0
        for i in range(atoms_count):
            off[i + 1] = off[i] + pa[i].degree
        # +1 so that a zero-bond molecule still gets a live pointer rather than NULL, which the
        # allocation check below would read as a failure.
        nbr = <uint32_t *> PyMem_Malloc((off[atoms_count] + 1) * sizeof(uint32_t))
        # ONE SLOT PER HALF-EDGE AND NOT PER BOND, deliberately.  A well-formed table names every bond
        # twice, so `bonds_count` slots are enough for it -- but `bonds_count` is computed FROM THE
        # DAMAGED HEADER, and a table where atom 0 names nine neighbours and none of them names it back
        # has a degree sum of 9, a `bonds_count` of 4, and up to 9 pairs to record.  Sizing the buffer
        # by the number of entries that can possibly be examined is a bound that holds whatever the
        # bytes say, which is the only kind of bound worth having here: the version of this line that
        # trusted `bonds_count` corrupted the heap on a mutated record and crashed thousands of records
        # later, in a decoder whose entire purpose is that it cannot be brought down by one bad record.
        edits = <edge_edit_t *> PyMem_Malloc((off[atoms_count] + 1) * sizeof(edge_edit_t))
        if nbr is NULL or edits is NULL:
            raise MemoryError('pach decode scratch allocation failed')

        for k in range(off[atoms_count]):
            value = _pach_number_at(data, length, table_at, k)
            if value < 0:
                nbr[k] = SU_NO_REF
                short_table = True
            elif slot_of[value] < 0:
                nbr[k] = SU_NO_REF
                problems.append('the connection table names atom %d, which the atom block does not '
                                'declare; the entry was dropped' % value)
            else:
                nbr[k] = <uint32_t> slot_of[value]
        if short_table:
            problems.append('the connection table needs %d bytes and the buffer ends before they do; '
                            'the entries past the end were dropped' % (3 * bonds_count))

        # A BOND'S ORDER IS SPENT WHERE THE TABLE FIRST NAMES THE PAIR, so this loop has to walk the
        # table in exactly chython 2's order -- atoms in block order, each atom's neighbours in table
        # order -- and consume one order value per fresh pair.  Everything that is wrong with a pair
        # is reported AFTER its order has been consumed, because dropping a bond must not shift the
        # orders of every bond after it.
        pos = 0
        for i in range(atoms_count):
            seen[i] = 1
            for j in range(off[i], off[i + 1]):
                other = nbr[j]
                if other != SU_NO_REF and other != i and seen[other]:
                    # the partner was walked already, so this half-bond's order was spent there --
                    # unless the partner's own list does not name this atom back, which is the
                    # asymmetric table chython 2 raised KeyError out of the C extension on.
                    symmetric = False
                    for k in range(off[other], off[other + 1]):
                        if nbr[k] == i:
                            symmetric = True
                            break
                    if not symmetric:
                        problems.append('atom %d names atom %d as a neighbour and is not named back; '
                                        'the half-bond was dropped'
                                        % (pa[i].number, pa[other].number))
                    continue
                order = _pach_order_at(data, length, order_at, pos, version)
                pos += 1
                if order < 0:
                    short_orders = True
                    order = 0
                if other == SU_NO_REF:
                    continue
                if other == i:
                    problems.append('atom %d names itself as a neighbour; the arena holds no self '
                                    'bonds and the entry was dropped' % pa[i].number)
                    continue
                dup = False
                for k in range(off[i], j):
                    if nbr[k] == other:
                        dup = True
                        break
                if dup:
                    problems.append('atom %d names atom %d twice; the second bond was dropped'
                                    % (pa[i].number, pa[other].number))
                    continue
                # RECIPROCITY IS REQUIRED IN THIS DIRECTION TOO.  The branch above catches the
                # asymmetric pair whose partner was walked first; this catches the one whose partner
                # comes later, and without it a table saying "0 is bonded to 5" while atom 5 says
                # nothing produced a bond that only one side of the record ever declared.  Both
                # directions drop it and report it once, which is the same answer read from either end.
                symmetric = False
                for k in range(off[other], off[other + 1]):
                    if nbr[k] == i:
                        symmetric = True
                        break
                if not symmetric:
                    problems.append('atom %d names atom %d as a neighbour and is not named back; the '
                                    'half-bond was dropped' % (pa[i].number, pa[other].number))
                    continue
                order += 1
                if order == 5 or order == 6 or order == 7:
                    problems.append('the bond between atoms %d and %d states order %d; three bits '
                                    'admit 1..8, the arena holds 1, 2, 3, 4 and 8, and it is stored '
                                    'as 8 -- present, order unspecified'
                                    % (pa[i].number, pa[other].number, order))
                    order = 8
                edits[n_edits].src = i
                edits[n_edits].dst = other
                edits[n_edits].order = <uint8_t> order
                n_edits += 1
        if short_orders:
            problems.append('the bond order block ends before the last bond; the orders it does not '
                            'reach were read as single')

        # A v2 record's configurations are per-atom signs plus a cis/trans block, and both are applied
        # AFTER the seal -- so whether this record states a configuration at all has to be answered
        # before `_pach_build` lays the arena out.  A false positive costs one byte per atom on a record
        # whose every sign turns out unusable; a false negative would be a refusal from
        # `structure_set_parity`, so the scan is over `sign` as stored and not over what survives.
        want_parity = ct_count != 0
        if not want_parity:
            for i in range(atoms_count):
                if pa[i].sign:
                    want_parity = True
                    break
        try:
            mol = _pach_build(pa, atoms_count, edits, n_edits, want_xy, want_parity)
        except ValueError as err:
            _pach_derivation_lost(problems, err)
            return (None, problems)
        _pach_apply_stereo(mol, problems, pa, atoms_count, nbr, off,
                           data, length, ct_at, ct_count, slot_of)
    finally:
        PyMem_Free(pa)
        PyMem_Free(slot_of)
        PyMem_Free(off)
        PyMem_Free(seen)
        PyMem_Free(nbr)
        PyMem_Free(edits)
    return (mol, problems)


cdef MoleculeContainer _pach_build(pach_atom_t *pa, uint32_t atoms_count,
                                   edge_edit_t *edits, uint32_t bonds_count, bint want_xy,
                                   bint want_parity):
    """Lay out the arena and wrap it in a container.

    The sequence is `from_bytes`'s and its comments are the authority for the order; the one thing
    worth noting is that `rebuild_derived` appends the derived segments and therefore REALLOCATES,
    so an `atom_t *` taken before it is a pointer into a freed block afterwards.

    NO COORDINATE SEGMENT is allocated when every coordinate is zero.  pach cannot tell "this record
    has no drawing" from "every atom sits at the origin" -- it writes 0 for both -- and of the two
    readings the first is the one that describes real stored data, since a molecule read from a SMILES
    string and packed has no coordinates at all.

    NO PARITY SEGMENT is allocated for a record that states no configuration, on the same principle:
    absent is unset, and the stored corpus is mostly stereo-free.
    """
    cdef Structure structure
    cdef atom_t *atoms
    cdef xy_t *xy
    cdef uint32_t i, seg_mask = 0
    cdef int rc
    if want_xy:
        seg_mask = SEG_MASK_XY
    if want_parity:
        # `_pach_apply_stereo` runs after this function seals the arena, so the segment it writes into
        # is named here or nowhere (the persistent block is laid out once).
        seg_mask |= SEG_MASK_PARITY
    structure = structure_alloc_full(atoms_count, bonds_count, False, seg_mask, NULL)
    atoms = structure.atoms()
    for i in range(atoms_count):
        atoms[i].element = pa[i].element
        atoms[i].charge = pa[i].charge
        atoms[i].isotope = pa[i].isotope
        atoms[i].n = pa[i].number
        # The explicit nibble is DERIVED and `rebuild_derived` fills it from the CSR; only the implicit
        # count is stored, and the sentinel goes through unchanged.
        at_set_h(&atoms[i], pa[i].hydrogens, 0)
        if pa[i].radical:
            at_set_radical(&atoms[i], True)
    with nogil:
        rc = csr_build(structure, edits, bonds_count)
    if rc:
        raise MemoryError('csr scratch allocation failed')
    if want_xy:
        xy = structure_xy(structure)
        for i in range(atoms_count):
            xy[i].x = pa[i].x
            xy[i].y = pa[i].y
    rebuild_derived(structure)
    atoms = structure.atoms()                      # ruling F60: the rebuild moved the buffer

    cdef uint32_t n
    cdef uint32_t high = 0
    cdef list numbers = []
    cdef dict index_of = {}
    for i in range(atoms_count):
        n = atoms[i].n
        numbers.append(n)
        index_of[n] = i
        if n > high:
            high = n
    cdef MoleculeContainer mol = MoleculeContainer.__new__(MoleculeContainer)
    mol._structure = structure
    mol._numbers = numbers
    mol._index_of = index_of
    mol._next_id = high + 1
    mol._first_pending = high + 1
    return mol


# --------------------------------------------------------------------------------------------------
# STEREO.  A sign is not a fact about an atom: it is a fact about an atom RELATIVE TO AN ORDER of its
# neighbours, and chython 2 and the arena order neighbours differently.  So the whole of the
# translation is
#
#     rebuild chython 2's order as arena slots  ->  `smi_perm_of`  ->  `translate_parity`
#
# and `translate_parity` is XOR with the permutation's parity, hence its own inverse: the encoder runs
# the identical three steps in the other direction.  There is no new permutation algebra here, which
# is the point -- chython 2's `_tetrahedron_translate` and `_alkene_translate` are permutation-parity
# tables and nothing else, so the two codebases' sign arithmetic is already the same function and only
# the frames and one convention bit per unit kind differ.
#
# THE FRAMES, from chython 2.24's `algorithms/stereo.py`, which is the only place they are written
# down:
#
#   tetrahedron  `stereogenic_tetrahedrons[n] = tuple(x for x in bonds[n] if atoms[x] != H)`, so the
#                non-hydrogen neighbours in `_bonds` insertion order -- which IS the pach connection
#                table's order for that atom, since the table is written by iterating `_bonds[n]`.
#                `_translate_tetrahedron_sign` appends the explicit hydrogen LAST when the frame has
#                only three heavy neighbours ("hydrogen always last in order").
#   cumulene     `stereogenic_cumulenes[path] = (nn[0], mn[0], sn, sm)` where `nn`/`mn` are each
#                terminal's neighbours excluding the chain, hydrogens and order-8 bonds, and
#                `sn`/`sm` are the second one or None.  So positions 0 and 2 belong to one terminal
#                and 1 and 3 to the other, and `_translate_cis_trans_sign`'s
#                `n2 is None and atoms[nn] == H` admits an explicit hydrogen into a None slot --
#                exactly the arena's own pair order (ruling F26: heavy, then a drawn hydrogen, then
#                nothing), which is why the two frames differ by a permutation and not by content.
#   allene       the same tuple, keyed by the chain's centre.
#
# WHICH TERMINAL LEADS DOES NOT MATTER -- exchanging the two ends of a cumulene frame is two
# transpositions and therefore even -- but WHICH PAIR OF THE ARENA'S `refs` each end is does, whenever
# both ends have an unnamed slot: `smi_perm_of` pairs unnamed positions off in order, so an end with an
# unnamed slot must be laid out against the `refs` pair that owns the matching unnamed slot.
# `_pach_cumulene_want` settles that by identity, on `refs[0]`, rather than by assuming which terminal
# the arena anchored the unit on.
# --------------------------------------------------------------------------------------------------

cdef bint _pach_tetra_frame(atom_t *atoms, uint32_t *nbr, uint32_t *off, uint32_t t,
                            uint32_t *want) noexcept nogil:
    """chython 2's tetrahedral frame for atom `t` as four arena slots. False when it has none."""
    cdef uint32_t k, x
    cdef uint32_t heavy = 0
    cdef uint32_t hydro = SU_NO_REF
    cdef uint32_t n_h = 0
    want[0] = SU_NO_REF
    want[1] = SU_NO_REF
    want[2] = SU_NO_REF
    want[3] = SU_NO_REF
    for k in range(off[t], off[t + 1]):
        x = nbr[k]
        if x == SU_NO_REF:
            continue
        if atoms[x].element == 1:
            n_h += 1
            if hydro == SU_NO_REF:
                hydro = x
        elif heavy < 4:
            want[heavy] = x
            heavy += 1
        else:
            return False                  # five heavy neighbours order no frame of four
    if n_h > 1:
        # `[C@](H)(H)(F)Cl`: two directions with no atom of their own cannot be told apart, so an
        # order over them is not an order and a sign against it names nothing.
        return False
    if heavy == 4:
        return n_h == 0
    if heavy != 3:
        return False
    want[3] = hydro                       # SU_NO_REF when the fourth direction has no atom at all
    return True


cdef bint _pach_end_pair(atom_t *atoms, uint32_t *ptr, halfedge_t *edges, uint32_t *nbr,
                         uint32_t *off, uint32_t t, uint32_t inward,
                         uint32_t *out) noexcept nogil:
    """chython 2's two directions of cumulene terminal `t`, whose chain neighbour is `inward`."""
    cdef uint32_t k, j, x
    cdef uint32_t n_heavy = 0
    cdef uint32_t hydro = SU_NO_REF
    cdef uint8_t order
    out[0] = SU_NO_REF
    out[1] = SU_NO_REF
    for k in range(off[t], off[t + 1]):
        x = nbr[k]
        if x == SU_NO_REF or x == inward:
            continue
        order = 0
        for j in range(ptr[t], ptr[t + 1]):
            if edges[j].to == x:
                order = edges[j].order
                break
        if order == 8:
            continue                      # chython 2's `b != 8`: no geometry on an unspecified bond
        if atoms[x].element == 1:
            if hydro == SU_NO_REF:
                hydro = x
        elif n_heavy < 2:
            out[n_heavy] = x
            n_heavy += 1
        else:
            return False
    if n_heavy == 0:
        return False                      # chython 2 requires a heavy substituent on both ends
    if n_heavy == 1:
        out[1] = hydro
    return True


cdef bint _pach_cumulene_want(stereo_unit_t *u, uint32_t *pair_a,
                              uint32_t *pair_b, uint32_t *want) noexcept nogil:
    """Interleave two end pairs into chython 2's `(nn0, mn0, sn, sm)` frame.

    The end holding `refs[0]` leads, so that the unnamed slots of `want` and of `refs` are in the same
    order and `smi_perm_of` pairs them correctly.  See the block comment above.
    """
    cdef uint32_t lead0, lead1, far0, far1
    if pair_a[0] == u.refs[0]:
        lead0 = pair_a[0]; lead1 = pair_a[1]; far0 = pair_b[0]; far1 = pair_b[1]
    elif pair_b[0] == u.refs[0]:
        lead0 = pair_b[0]; lead1 = pair_b[1]; far0 = pair_a[0]; far1 = pair_a[1]
    else:
        return False
    want[0] = lead0
    want[1] = far0
    want[2] = lead1
    want[3] = far1
    return True


cdef bint _pach_allene_ends(uint32_t *ptr, halfedge_t *edges, uint32_t centre,
                            uint32_t *terms, uint32_t *inwards) noexcept nogil:
    """The two terminals of the cumulene chain centred on `centre`, and their inward neighbours."""
    cdef uint32_t k, cur, prev, nxt
    cdef uint32_t n = 0
    for k in range(ptr[centre], ptr[centre + 1]):
        if not _is_chain_bond(&edges[k]):
            continue
        if n == 2:
            return False
        prev = centre
        cur = edges[k].to
        while True:
            nxt = _chain_next(ptr, edges, cur, prev)
            if nxt == SU_NO_REF:
                break
            prev = cur
            cur = nxt
        terms[n] = cur
        inwards[n] = prev
        n += 1
    return n == 2


cdef int _pach_apply_stereo(MoleculeContainer mol, list problems, pach_atom_t *pa,
                            uint32_t atoms_count, uint32_t *nbr, uint32_t *off,
                            const unsigned char *data, Py_ssize_t length, Py_ssize_t ct_at,
                            uint32_t ct_count, int32_t *slot_of) except -1:
    """Write the record's configurations into SEG_PARITY.

    UNMARKED, for `smi_stereo`'s reason (ruling F70): every question asked here is about constitution
    -- what kind of unit is anchored here, what its reference directions are -- and none of them is
    "is this stereogenic?".  A record that states a configuration is storing what the record said.
    """
    cdef Structure structure = mol._structure
    cdef uint32_t i, k, kind_true
    cdef uint32_t sn, sm, anchor, partner, inward, far, far_prev, n_chain
    cdef uint32_t want[4]
    cdef uint32_t perm[4]
    cdef uint32_t pair_a[2]
    cdef uint32_t pair_b[2]
    cdef uint32_t terms[2]
    cdef uint32_t inwards[2]
    cdef uint8_t parity
    cdef int value
    cdef bint wrote = False
    cdef bint any_sign = False
    cdef atom_t *atoms
    cdef uint32_t *ptr
    cdef halfedge_t *edges
    cdef stereo_unit_t *u

    for i in range(atoms_count):
        if pa[i].sign:
            any_sign = True
            break
    if not any_sign and not ct_count:
        return 0

    # DERIVING THE UNIT TABLE CAN FAIL, and on a damaged record it must not become this function's
    # exception.  The core refuses a graph in which two stereo units claim one anchor atom -- a parity
    # is keyed by anchor slot in SEG_PARITY, so the second unit would overwrite the first -- and a
    # mutated record can decode to such a graph: one mutant in forty-one thousand did.  It is not a
    # defect of the decoder, and the same graph built through `add_atom`/`add_bond` refuses the same
    # way; but a caller looping over a store gets the graph and the sentence rather than a traceback,
    # and can then decide.
    try:
        ensure_stereo_units_unmarked(structure)
    except Exception as err:
        problems.append('this record\'s graph has no usable stereo unit table (%s), so the %d '
                        'configuration(s) it states were dropped'
                        % (err, ct_count + (1 if any_sign else 0)))
        return 0
    atoms = structure.atoms()                      # after, not before: the ensure may reallocate
    ptr = csr_ptr(structure)
    edges = csr_edges(structure)

    # ---- the atom nibbles: tetrahedral centres, and allene centres, whose sign chython 2 also kept
    # on the atom.  Which of the two a given atom is is a question about the graph, and the graph is
    # asked rather than the nibble's two-bit field, which chython 2's reader collapses (see the header).
    for i in range(atoms_count):
        if not pa[i].sign:
            continue
        u = stereo_unit_of(structure, i)
        if u is NULL:
            problems.append('atom %d carries a configuration and anchors no stereo unit in this '
                            'molecule; there is nowhere to put it and it was dropped' % pa[i].number)
            continue
        if u.n_refs != 4:
            problems.append('atom %d carries a configuration over %d reference direction(s); four is '
                            'what an ordered frame needs, so it was dropped'
                            % (pa[i].number, u.n_refs))
            continue
        if u.kind == SU_TETRA:
            kind_true = PACH_TETRA_TRUE
            if not _pach_tetra_frame(atoms, nbr, off, i, want):
                problems.append('atom %d carries a configuration and its neighbours do not form a '
                                'frame chython 2 could have measured it against; it was dropped'
                                % pa[i].number)
                continue
        elif u.kind == SU_ALLENE:
            kind_true = PACH_ALLENE_TRUE
            if not _pach_allene_ends(ptr, edges, i, terms, inwards):
                problems.append('atom %d anchors an allene whose chain this record does not describe; '
                                'the configuration was dropped' % pa[i].number)
                continue
            if not _pach_end_pair(atoms, ptr, edges, nbr, off, terms[0], inwards[0], pair_a) \
                    or not _pach_end_pair(atoms, ptr, edges, nbr, off, terms[1], inwards[1], pair_b) \
                    or not _pach_cumulene_want(u, pair_a, pair_b, want):
                problems.append('atom %d anchors an allene whose terminals do not form a frame '
                                'chython 2 could have measured its sign against; it was dropped'
                                % pa[i].number)
                continue
        else:
            problems.append('atom %d carries a configuration and anchors %s, which pach has no field '
                            'for and chython 2 never wrote; it was dropped'
                            % (pa[i].number, smi_kind_name(u.kind)))
            continue
        smi_perm_of(u, want, perm)
        parity = translate_parity(<uint8_t> (kind_true if pa[i].sign == 2 else 3 - kind_true), perm)
        structure_set_parity(structure, i, parity)
        wrote = True

    # ---- the cis/trans block.  Its entries name the two TERMINALS of an even-length cumulene chain
    # and chython 2 kept the sign on the chain's central BOND; the arena keeps it on the terminal it
    # anchored the unit at, which is why this loop looks the unit up from either end.
    for k in range(ct_count):
        # THE ENTRY'S STRIDE IS FOUR BYTES AND ITS NUMBER PAIR OCCUPIES THREE, so the pair is addressed
        # as the two numbers of ITS OWN triple and not as numbers 2k and 2k+1 of one long table.  The
        # difference is invisible in a record with a single entry and shifts every later entry by one
        # byte in a record with two, which is exactly the shape of bug a fixture corpus catches.
        value = _pach_number_at(data, length, ct_at + 4 * <Py_ssize_t> k, 0)
        if value < 0:
            problems.append('the cis/trans block declares %d entries and the buffer ends after %d; '
                            'the rest were dropped' % (ct_count, k))
            break
        sn = <uint32_t> value
        value = _pach_number_at(data, length, ct_at + 4 * <Py_ssize_t> k, 1)
        if value < 0:
            problems.append('the cis/trans block declares %d entries and the buffer ends after %d; '
                            'the rest were dropped' % (ct_count, k))
            break
        sm = <uint32_t> value
        if ct_at + 4 * <Py_ssize_t> k + 4 > length:
            problems.append('the cis/trans entry naming atoms %d and %d has no sign byte; it was '
                            'dropped' % (sn, sm))
            break
        # THE WHOLE BYTE, not bit 0.  chython 2 wrote the sign into bit 0 of a byte it documents as
        # padding and read it back as `if d:`, so a record whose pad bits are set was a True to the
        # writer's own reader.  Reading only bit 0 here would report a different molecule than the one
        # that was stored.
        parity = 2 if data[ct_at + 4 * <Py_ssize_t> k + 3] else 1
        if sn > PACH_MAX_NUMBER or sm > PACH_MAX_NUMBER \
                or slot_of[sn] < 0 or slot_of[sm] < 0:
            problems.append('a cis/trans entry names atoms %d and %d and the atom block declares at '
                            'least one of them nowhere; it was dropped' % (sn, sm))
            continue
        anchor = <uint32_t> slot_of[sn]
        partner = <uint32_t> slot_of[sm]
        u = stereo_unit_of(structure, anchor)
        if u is NULL or u.kind != SU_CIS_TRANS or stereo_unit_partner(structure, u) != partner:
            anchor = <uint32_t> slot_of[sm]
            partner = <uint32_t> slot_of[sn]
            u = stereo_unit_of(structure, anchor)
            if u is NULL or u.kind != SU_CIS_TRANS or stereo_unit_partner(structure, u) != partner:
                problems.append('a cis/trans entry names atoms %d and %d, which anchor no cis/trans '
                                'unit in this molecule; it was dropped' % (sn, sm))
                continue
        if u.n_refs != 4:
            problems.append('the cis/trans unit of atoms %d and %d orders %d reference direction(s); '
                            'four is what a sign needs and it was dropped' % (sn, sm, u.n_refs))
            continue
        inward = _chain_next(ptr, edges, anchor, SU_NO_REF)
        if inward == SU_NO_REF or not _cumulene_walk(atoms, ptr, edges, anchor,
                                                     &far, &far_prev, &n_chain):
            problems.append('a cis/trans entry names atoms %d and %d and the chain between them is '
                            'not one this molecule holds; it was dropped' % (sn, sm))
            continue
        if not _pach_end_pair(atoms, ptr, edges, nbr, off, anchor, inward, pair_a) \
                or not _pach_end_pair(atoms, ptr, edges, nbr, off, far, far_prev, pair_b) \
                or not _pach_cumulene_want(u, pair_a, pair_b, want):
            problems.append('a cis/trans entry names atoms %d and %d, whose terminals do not form a '
                            'frame chython 2 could have measured its sign against; it was dropped'
                            % (sn, sm))
            continue
        smi_perm_of(u, want, perm)
        parity = translate_parity(<uint8_t> (PACH_CIS_TRANS_TRUE if parity == 2
                                             else 3 - PACH_CIS_TRANS_TRUE), perm)
        structure_set_parity(structure, anchor, parity)
        wrote = True

    if wrote:
        refresh_parity_features(structure)
    return 0


# --------------------------------------------------------------------------------------------------
# THE ENCODER.  Writing the bit fields is the easy half; the hard half is that a writer must not lose
# anything silently, so every field the arena holds and pach has no slot for is a NAMED REFUSAL that
# the caller can waive one name at a time.
# --------------------------------------------------------------------------------------------------

DEF PACH_DROP_MAP_NUMBER = 0x01
DEF PACH_DROP_TITLE = 0x02
DEF PACH_DROP_SGROUPS = 0x04
DEF PACH_DROP_CIP = 0x08
DEF PACH_DROP_WEDGES = 0x10
DEF PACH_DROP_STEREO_GROUPS = 0x20
DEF PACH_DROP_STEREO = 0x40
DEF PACH_DROP_META = 0x80
# Not a loss the LEGACY versions can be asked to take -- 0 and 2 always carry a coordinate field --
# so `_pach_refuse_losses` never reads this bit.  It selects version 4 over version 3.
DEF PACH_DROP_COORDINATES = 0x100
# Read by `_pach3_refuse_losses` only: version 2 has no conformer field either, and its refusal set is
# frozen with the format.
DEF PACH_DROP_CONFORMERS = 0x200
DEF PACH_DROP_ALL = 0x3ff

cdef dict _PACH_DROP_NAMES = {'map_number': PACH_DROP_MAP_NUMBER, 'title': PACH_DROP_TITLE,
                              'sgroups': PACH_DROP_SGROUPS, 'cip': PACH_DROP_CIP,
                              'wedges': PACH_DROP_WEDGES,
                              'stereo_groups': PACH_DROP_STEREO_GROUPS,
                              'stereo': PACH_DROP_STEREO, 'meta': PACH_DROP_META,
                              'coordinates': PACH_DROP_COORDINATES,
                              'conformers': PACH_DROP_CONFORMERS}


cdef inline void _pach_put_number(unsigned char *p, uint32_t k, uint32_t v) noexcept nogil:
    """The k-th 12-bit number of a table, two to three bytes. The buffer starts zeroed and the even
    index of a triple is always written before its odd one, which is why the low nibble is an OR."""
    cdef uint32_t at = 3 * (k >> 1)
    if k & 1:
        p[at + 1] |= <unsigned char> (v >> 8)
        p[at + 2] = <unsigned char> v
    else:
        p[at] = <unsigned char> (v >> 4)
        p[at + 1] = <unsigned char> (v << 4)


cdef inline void _pach_put_order(unsigned char *p, uint32_t k, uint8_t v) noexcept nogil:
    """The k-th 3-bit order value of a version 2 block: a flat MSB-first bitstream."""
    cdef uint32_t bit = 3 * k
    cdef uint32_t at = bit >> 3
    bit &= 7
    if bit <= 5:
        p[at] |= <unsigned char> (v << (5 - bit))
    else:
        p[at] |= <unsigned char> (v >> (bit - 5))
        p[at + 1] |= <unsigned char> (v << (13 - bit))


cdef bint _pach_chain_middle(uint32_t *ptr, halfedge_t *edges, uint32_t a, uint32_t b,
                             uint32_t *ma, uint32_t *mb) noexcept nogil:
    """The CENTRAL BOND of the even cumulene chain from terminal `a` to terminal `b`.

    chython 2 kept a cis/trans sign on that bond and emitted the record's cis/trans entry where the
    connection table first names it, so a writer that wants its blocks in chython 2's order has to know
    which bond it is.  64 chain atoms is a bound and not a limit -- the longest cumulene anybody has
    made is a few dozen -- and a chain past it is reported by the caller rather than truncated.
    """
    cdef uint32_t path[64]
    cdef uint32_t cnt = 1
    cdef uint32_t prev = SU_NO_REF
    cdef uint32_t cur = a
    cdef uint32_t nxt
    path[0] = a
    while True:
        nxt = _chain_next(ptr, edges, cur, prev)
        if nxt == SU_NO_REF:
            break
        if cnt >= 64:
            return False
        path[cnt] = nxt
        cnt += 1
        prev = cur
        cur = nxt
    if cur != b or cnt < 2 or (cnt & 1):
        return False
    ma[0] = path[cnt // 2 - 1]
    mb[0] = path[cnt // 2]
    return True


cdef int _pach_refuse_losses(MoleculeContainer mol, uint32_t drop_mask) except -1:
    """Everything the arena holds that pach has no field for, refused by name.

    Not a validation pass: a validator answers "is this molecule legal", and every molecule here is.
    This answers "would writing it lose something", which is a question about the FORMAT, and the only
    honest answers are "no", "yes and here is what" and "yes and you said you did not mind".
    """
    cdef Structure structure = mol._structure
    cdef atom_t *atoms = structure.atoms()
    cdef halfedge_t *edges = csr_edges(structure)
    cdef uint32_t n = structure.header.atom_count
    cdef uint32_t i, k
    if not (drop_mask & PACH_DROP_MAP_NUMBER):
        for i in range(n):
            if atoms[i].map_number:
                raise ValueError('atom %d carries map_number %d and the pach format has no field for '
                                 'one; pass drop=[\'map_number\'] to write the record without it, or '
                                 'to_bytes() to keep it'
                                 % (atoms[i].n, atoms[i].map_number))
    if not (drop_mask & PACH_DROP_CIP):
        for i in range(n):
            if atoms[i].reserved & ATOM_CIP_MASK:
                raise ValueError('atom %d carries a cip descriptor and the pach format has no field '
                                 'for one; pass drop=[\'cip\'] to write the record without it'
                                 % atoms[i].n)
        for k in range(2 * structure.header.bond_count):
            if edges[k].flags & HE_CIP_MASK:
                raise ValueError('a bond carries a cip descriptor and the pach format has no field '
                                 'for one; pass drop=[\'cip\'] to write the record without it')
    if not (drop_mask & PACH_DROP_WEDGES):
        for k in range(2 * structure.header.bond_count):
            if edges[k].wedge:
                raise ValueError('a bond carries a wedge and the pach format has no field for one; '
                                 'pass drop=[\'wedges\'] to write the record without it')
    if not (drop_mask & PACH_DROP_STEREO_GROUPS) and structure_has(structure, SEG_STEREO_GROUPS):
        raise ValueError('this molecule carries enhanced stereo_groups and the pach format has no '
                         'field for them; pass drop=[\'stereo_groups\'] to write the record without '
                         'them')
    if not (drop_mask & PACH_DROP_SGROUPS) and structure_sgroup_count(structure):
        raise ValueError('this molecule carries %d sgroups and the pach format has no field for them; '
                         'pass drop=[\'sgroups\'] to write the record without them'
                         % structure_sgroup_count(structure))
    if not (drop_mask & PACH_DROP_TITLE) and len(blob_bytes(structure, SEG_OPAQUE_BLOB, 0)):
        raise ValueError('this molecule carries a title and the pach format has no text of any kind; '
                         'pass drop=[\'title\'] to write the record without it')
    # `_meta` and not `meta`, so asking the question does not create the dict it is asking about
    if not (drop_mask & PACH_DROP_META) and mol._meta:
        raise ValueError('this molecule carries %d metadata key(s) and the pach format has no field '
                         'for any of them; pass drop=[\'meta\'] to write the record without them'
                         % len(mol._meta))
    return 0


cdef bytes _pach_encode(MoleculeContainer mol, uint32_t drop_mask):
    """One version 2 pach record, uncompressed.

    ONLY VERSION 2 IS EVER WRITTEN.  Version 0 differs in the bond order block alone and every reader
    that understands 0 understands 2; a writer with a version switch would be offering a choice whose
    only possible use is to make a record harder to read.

    THE CONNECTION TABLE IS EMITTED IN CSR ORDER, which is ascending by neighbour slot, where chython 2
    emitted its `_bonds` insertion order.  That is the first of the two reasons a re-encoded record is
    not byte-identical to the original, and it is not fixable: the arena does not store an insertion
    order, so there is nothing to reproduce.  It costs nothing, because the order block is consumed in
    connection-table order and this writer's own table is what its own reader walks.
    """
    cdef Structure structure = mol._structure
    cdef uint32_t n, nb, i, j, k, e, deg, num, ct_count = 0, pos = 0, npos = 0, written = 0
    cdef uint32_t partner, inward, far, far_prev, n_chain
    cdef Py_ssize_t size, table_at, order_at, ct_at, atom_at
    cdef int shift, charge
    cdef uint8_t h, parity, v2, kind_true, nibble
    cdef uint32_t want[4]
    cdef uint32_t perm[4]
    cdef uint32_t pair_a[2]
    cdef uint32_t pair_b[2]
    cdef uint32_t terms[2]
    cdef uint32_t inwards[2]
    cdef atom_t *atoms
    cdef uint32_t *ptr
    cdef halfedge_t *edges
    cdef xy_t *xy = NULL
    cdef stereo_unit_t *u
    cdef unsigned char *buf = NULL
    cdef uint8_t *nib = NULL
    cdef uint8_t *seen = NULL
    cdef uint32_t *nbr = NULL
    cdef uint32_t *ct_n = NULL
    cdef uint32_t *ct_m = NULL
    cdef uint32_t *ct_a = NULL
    cdef uint32_t *ct_b = NULL
    cdef uint8_t *ct_s = NULL
    cdef bytes out

    mol._require_clean()
    _pach_refuse_losses(mol, drop_mask)
    # The unit table is DERIVED and unmarked, which is what the decoder built its frames against; the
    # two directions have to ask the same question or the sign they exchange means two things.  It can
    # reallocate, so every pointer below is taken after it.
    ensure_stereo_units_unmarked(structure)
    structure = mol._structure
    atoms = structure.atoms()
    ptr = csr_ptr(structure)
    edges = csr_edges(structure)
    n = structure.header.atom_count
    nb = structure.header.bond_count
    if structure_has(structure, SEG_XY):
        xy = structure_xy(structure)

    # ---- per-atom structural refusals: the R-marker check (a field the pach format has no room for)
    # comes first, then field-width bounds.  A missing field writes an unreadable record; a wrapped
    # field decodes to a different molecule.  Both are worse than not writing.
    for i in range(n):
        if atoms[i].element == 0:
            raise ValueError('atom %d is an R marker and the pach format has no field for one: its '
                             'atom block is a fixed 9 bytes and byte 3 holds the atomic number with '
                             'no room for an index. Use `to_bytes`, which is lossless, or write '
                             'SMILES or a CTfile' % atoms[i].n)
        if atoms[i].n > PACH_MAX_NUMBER:
            raise ValueError('atom number %d does not fit the pach format\'s 12 bit atom number '
                             'field, whose greatest value is 4095; arena stable ids are 32 bit and '
                             'never reused, so an edited molecule can outgrow the format. remap() to '
                             'small numbers, or use to_bytes()' % atoms[i].n)
        deg = ptr[i + 1] - ptr[i]
        if deg > PACH_MAX_DEGREE:
            raise ValueError('atom %d has %d neighbours and the pach format\'s neighbour count field '
                             'is 4 bits, so 15 is the most it can state'
                             % (atoms[i].n, deg))
        h = at_implicit_h(&atoms[i])
        if h != H_UNKNOWN and h > PACH_MAX_IMPLICIT_H:
            raise ValueError('atom %d states %d implicit hydrogens and the pach format\'s field is 3 '
                             'bits with 7 reserved for "not stated", so 6 is the most it can hold. '
                             'chython 2 wrote this field as `<unsigned char> count << 5` and silently '
                             'turned 8 into 0; this writer refuses instead'
                             % (atoms[i].n, h))
        charge = atoms[i].charge
        if charge < -4 or charge > 11:
            raise ValueError('atom %d states formal charge %+d and the pach format holds charge+4 in '
                             '4 bits, so -4..+11 is its whole range'
                             % (atoms[i].n, charge))
        if atoms[i].isotope:
            shift = <int> atoms[i].isotope - <int> MDL_ISOTOPE[atoms[i].element] + PACH_ISOTOPE_BIAS
            if shift < 1 or shift > 31:
                raise ValueError('atom %d states isotope %d and the pach format stores an isotope as '
                                 'a 5 bit offset from its element\'s MDL reference mass %d, so only '
                                 '%d..%d can be written'
                                 % (atoms[i].n, atoms[i].isotope,
                                    MDL_ISOTOPE[atoms[i].element],
                                    <int> MDL_ISOTOPE[atoms[i].element] - 15,
                                    <int> MDL_ISOTOPE[atoms[i].element] + 15))

    nib = <uint8_t *> PyMem_Malloc((n + 1) * sizeof(uint8_t))
    seen = <uint8_t *> PyMem_Malloc((n + 1) * sizeof(uint8_t))
    nbr = <uint32_t *> PyMem_Malloc((2 * nb + 1) * sizeof(uint32_t))
    ct_n = <uint32_t *> PyMem_Malloc((n + 1) * sizeof(uint32_t))
    ct_m = <uint32_t *> PyMem_Malloc((n + 1) * sizeof(uint32_t))
    ct_a = <uint32_t *> PyMem_Malloc((n + 1) * sizeof(uint32_t))
    ct_b = <uint32_t *> PyMem_Malloc((n + 1) * sizeof(uint32_t))
    ct_s = <uint8_t *> PyMem_Malloc((n + 1) * sizeof(uint8_t))
    if nib is NULL or seen is NULL or nbr is NULL or ct_n is NULL or ct_m is NULL \
            or ct_a is NULL or ct_b is NULL or ct_s is NULL:
        PyMem_Free(nib); PyMem_Free(seen); PyMem_Free(nbr)
        PyMem_Free(ct_n); PyMem_Free(ct_m); PyMem_Free(ct_a); PyMem_Free(ct_b); PyMem_Free(ct_s)
        raise MemoryError('pach encode scratch allocation failed')
    try:
        memset(nib, 0, (n + 1) * sizeof(uint8_t))
        memset(seen, 0, (n + 1) * sizeof(uint8_t))
        # THE FRAME BUILDERS TAKE A NEIGHBOUR LIST AND ITS OFFSETS, and the CSR is exactly that, so the
        # writer hands them the same functions the reader used with `ptr` as the offset array.  One
        # frame source for both directions is what makes the round trip a fixed point rather than a
        # coincidence.
        for k in range(2 * nb):
            nbr[k] = edges[k].to

        # ---- the configurations.  Atom kinds become a nibble in the atom block, cis/trans becomes an
        # entry in the trailing block, and everything else is a loss with a name.
        for i in range(n):
            parity = structure_parity_at(structure, i)
            if not parity:
                continue
            u = stereo_unit_of(structure, i)
            if u is NULL or u.n_refs != 4:
                if drop_mask & PACH_DROP_STEREO:
                    continue
                raise ValueError('atom %d carries a configuration that no stereo unit of this '
                                 'molecule can express, so the pach format has nowhere to put it; '
                                 'pass drop=[\'stereo\'] to write the record without it'
                                 % atoms[i].n)
            if u.kind == SU_TETRA or u.kind == SU_ALLENE:
                if u.kind == SU_TETRA:
                    kind_true = PACH_TETRA_TRUE
                    if not _pach_tetra_frame(atoms, nbr, ptr, i, want):
                        if drop_mask & PACH_DROP_STEREO:
                            continue
                        raise ValueError('atom %d carries a configuration and its neighbours do not '
                                         'form a frame the pach format can state it against; pass '
                                         'drop=[\'stereo\'] to write the record without it'
                                         % atoms[i].n)
                else:
                    kind_true = PACH_ALLENE_TRUE
                    if not _pach_allene_ends(ptr, edges, i, terms, inwards) \
                            or not _pach_end_pair(atoms, ptr, edges, nbr, ptr,
                                                  terms[0], inwards[0], pair_a) \
                            or not _pach_end_pair(atoms, ptr, edges, nbr, ptr,
                                                  terms[1], inwards[1], pair_b) \
                            or not _pach_cumulene_want(u, pair_a, pair_b, want):
                        if drop_mask & PACH_DROP_STEREO:
                            continue
                        raise ValueError('atom %d anchors an allene whose terminals do not form a '
                                         'frame the pach format can state a sign against; pass '
                                         'drop=[\'stereo\'] to write the record without it'
                                         % atoms[i].n)
                smi_perm_of(u, want, perm)
                v2 = translate_parity(parity, perm)
                # THE FIELD IS CHOSEN BY THE NEIGHBOUR COUNT, not by the unit kind, because that is
                # what chython 2's writer did: two neighbours got the allene pair of bits and anything
                # else the tetrahedron pair.  Its reader collapses all four values anyway, so the choice
                # only matters for byte identity -- and reproducing it costs one comparison.
                deg = ptr[i + 1] - ptr[i]
                if deg == 2:
                    nibble = 0x30 if v2 == kind_true else 0x20
                else:
                    nibble = 0xc0 if v2 == kind_true else 0x80
                nib[i] = nibble
            elif u.kind == SU_CIS_TRANS:
                kind_true = PACH_CIS_TRANS_TRUE
                partner = stereo_unit_partner(structure, u)
                inward = _chain_next(ptr, edges, i, SU_NO_REF)
                if partner == SU_NO_REF or inward == SU_NO_REF \
                        or not _cumulene_walk(atoms, ptr, edges, i, &far, &far_prev, &n_chain) \
                        or far != partner \
                        or not _pach_chain_middle(ptr, edges, i, partner,
                                                  &ct_a[ct_count], &ct_b[ct_count]) \
                        or not _pach_end_pair(atoms, ptr, edges, nbr, ptr, i, inward, pair_a) \
                        or not _pach_end_pair(atoms, ptr, edges, nbr, ptr, partner, far_prev,
                                              pair_b) \
                        or not _pach_cumulene_want(u, pair_a, pair_b, want):
                    if drop_mask & PACH_DROP_STEREO:
                        continue
                    raise ValueError('the cis/trans configuration anchored at atom %d has no frame '
                                     'the pach format can state it against; pass drop=[\'stereo\'] '
                                     'to write the record without it' % atoms[i].n)
                smi_perm_of(u, want, perm)
                v2 = translate_parity(parity, perm)
                ct_n[ct_count] = atoms[i].n
                ct_m[ct_count] = atoms[partner].n
                ct_s[ct_count] = 1 if v2 == kind_true else 0
                ct_count += 1
            elif drop_mask & PACH_DROP_STEREO:
                continue
            else:
                raise ValueError('atom %d carries a configuration of a kind the pach format has no '
                                 'field for -- it has three, tetrahedral, allene and cis/trans; pass '
                                 'drop=[\'stereo\'] to write the record without it'
                                 % atoms[i].n)

        table_at = 4 + 9 * <Py_ssize_t> n
        order_at = table_at + 3 * <Py_ssize_t> nb
        ct_at = order_at + _pach_order_block_len(nb, 2)
        size = ct_at + 4 * <Py_ssize_t> ct_count
        buf = <unsigned char *> PyMem_Malloc(size)
        if buf is NULL:
            raise MemoryError('pach encode buffer allocation failed')
        try:
            memset(buf, 0, size)
            buf[0] = 2
            buf[1] = <unsigned char> (n >> 4)
            buf[2] = <unsigned char> ((n << 4) | (ct_count >> 8))
            buf[3] = <unsigned char> ct_count

            for i in range(n):
                atom_at = 4 + 9 * <Py_ssize_t> i
                num = atoms[i].n
                deg = ptr[i + 1] - ptr[i]
                buf[atom_at] = <unsigned char> (num >> 4)
                buf[atom_at + 1] = <unsigned char> ((num << 4) | deg)
                shift = 0
                if atoms[i].isotope:
                    shift = (<int> atoms[i].isotope - <int> MDL_ISOTOPE[atoms[i].element]
                             + PACH_ISOTOPE_BIAS)
                buf[atom_at + 2] = <unsigned char> (nib[i] | (shift >> 1))
                buf[atom_at + 3] = <unsigned char> (((shift & 1) << 7) | atoms[i].element)
                if xy is not NULL:
                    _pach_f16_encode(xy_read_x(xy + i), buf + atom_at + 4)
                    _pach_f16_encode(xy_read_y(xy + i), buf + atom_at + 6)
                h = at_implicit_h(&atoms[i])
                if h == H_UNKNOWN:
                    h = PACH_H_UNKNOWN
                buf[atom_at + 8] = <unsigned char> ((h << 5) | ((atoms[i].charge + 4) << 1)
                                                    | (1 if at_radical(&atoms[i]) else 0))

            # ---- one walk writes the connection table, the order block and the cis/trans block, and
            # that is not an optimisation: an order is spent where the table FIRST names a pair, so the
            # two blocks are the same traversal seen twice, and chython 2 emitted the cis/trans entry
            # in the same place.  Writing them in three passes would mean stating the rule three times.
            for i in range(n):
                seen[i] = 1
                for k in range(ptr[i], ptr[i + 1]):
                    j = edges[k].to
                    _pach_put_number(buf + table_at, npos, atoms[j].n)
                    npos += 1
                    if seen[j]:
                        continue
                    _pach_put_order(buf + order_at, pos, <uint8_t> (edges[k].order - 1))
                    pos += 1
                    for e in range(ct_count):
                        if (ct_a[e] == i and ct_b[e] == j) or (ct_a[e] == j and ct_b[e] == i):
                            _pach_put_number(buf + ct_at + 4 * <Py_ssize_t> written, 0, ct_n[e])
                            _pach_put_number(buf + ct_at + 4 * <Py_ssize_t> written, 1, ct_m[e])
                            buf[ct_at + 4 * <Py_ssize_t> written + 3] = ct_s[e]
                            written += 1
                            break
            out = <bytes> buf[:size]
        finally:
            PyMem_Free(buf)
    finally:
        PyMem_Free(nib); PyMem_Free(seen); PyMem_Free(nbr)
        PyMem_Free(ct_n); PyMem_Free(ct_m); PyMem_Free(ct_a); PyMem_Free(ct_b); PyMem_Free(ct_s)
    return out


# --------------------------------------------------------------------------------------------------
# THE PYTHON DOORS.  Three functions and two container methods, and the split between them is the whole
# of the garbage policy: `pach_load` is the loop-safe door and never raises on a record's content;
# `MoleculeContainer.unpack` is the ANSWER boundary and raises when a caller who asked for a molecule
# cannot have one.  One decoder underneath both.
# --------------------------------------------------------------------------------------------------

cdef bytes _pach_decompress(object data, list problems):
    """`data` through zlib, or None with a sentence saying why not."""
    try:
        return zlib.decompress(data)
    except Exception as err:
        problems.append('the buffer is not zlib compressed data: %s' % err)
        return None
