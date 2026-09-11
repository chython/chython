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
# The molecule arena: `Structure`, one allocation carved into numbered segments.
#
# A segment's layout is part of the serialised format, so a packed buffer and a live arena hold
# the same bytes and `structure_from_bytes` is a copy rather than a parse.  This file holds the
# storage and its accessors only; the derived encoding the matcher screens on is in
# `_features.pxi`, and the element tables it reads are in `_elements.pxi`.


# Segments come in two blocks, and only the first one is a format.
#
# PERSISTENT segments are serialised. Their ids are part of the on-disk format and are FROZEN
# FOREVER: a packed buffer's segment table is indexed by these numbers, so a new persistent segment
# APPENDS to this block and never reorders it. Ids 0-4 are also v3's ids, unchanged, which is what
# lets this build read a v3 buffer without moving a byte.
#
# DERIVED segments are rebuilt on demand and never serialised. Their ids are COMPILE-TIME ONLY --
# they name a slot in the `Structure` object, not an entry in any buffer -- so renumbering them is
# free. That is exactly what makes adding a persistent segment cheap, and it is why v3's table ran
# out at twelve of thirteen: seven of those twelve were derived caches occupying slots in a
# serialised table that never carried them.
cdef enum:
    SEG_ATOMS = 0
    SEG_CSR_PTR = 1
    SEG_CSR_EDGE = 2
    SEG_XY = 3
    SEG_STEREO_GROUPS = 4
    SEG_SGROUP_RECORD = 5
    SEG_SGROUP_INDEX = 6
    # OPAQUE, AND THE NAME HAS TO SAY SO.  This segment holds S-group payload the arena must never
    # interpret: `data` fields whose bytes may not be valid UTF-8, and whose fidelity requirement is
    # precisely that NOBODY DECODES THEM -- a decode at store time loses an undecodable byte for good.
    # A name mentioning text invites exactly the convenience decode the invariant forbids.  A SEGMENT
    # NAME IS FREE TO CHANGE WHILE ONLY THE CORE READS IT; the NUMBER is not, and that is CHECKED
    # rather than asserted: `to_bytes()` is a molecule identity, so a renumbering would reprice every
    # stored key -- `test_arena_v3_compat.py` reads frozen v3 bytes and
    # `test_segment_ids_are_dense_and_unique` pins this id at 7.
    SEG_OPAQUE_BLOB = 7
    # THE THIRD COORDINATE, AND IT IS NOT IN `xy_t`.  A conformer is a whole set of positions for
    # every atom, so a molecule may hold several -- a multi-record SDF, a PDB `MODEL` stack, extXYZ
    # frames, a generator's output.  A one-z-per-atom layout (`SEG_Z`, not adopted) cannot hold two of
    # them, so a reader facing a file that states three drops two at the input boundary.  Layout is in
    # `_conformers.pxi`; the design is
    # `docs/superpowers/specs/2026-09-05-chython3-seg-conformers-design.md`.
    #
    # APPENDED rather than inserted, which is the whole reason this cost nothing: no stored buffer
    # moves a byte, because every reader sizes the table from the buffer's own `seg_count` and a
    # molecule with no conformers writes neither the segment nor its table entry.
    SEG_CONFORMERS = 8
    # THE STATED PARITY, one byte per atom: 0 none, 1 even, 2 odd.  PERSISTENT because it is the only
    # copy -- a parity is what the input said, and `SEG_STEREO_UNIT` beside it is a derived table of
    # FRAMES a parity is read against.  One byte per atom rather than a record per configuration for
    # the reason `structure_alloc_full` gives: a record count is a function of perception, which needs
    # the built CSR, so it is not knowable when the persistent block is laid out.
    #
    # ABSENT IS UNSET, so a molecule with no stereo writes no segment and pays no byte, exactly as
    # SEG_XY and SEG_STEREO_GROUPS do.  A write into an absent segment would land on the shared zero
    # page: `structure_set_parity` refuses it rather than checking for it at every read.
    #
    # Two molecules differing only in a stated parity are different bytes and so different molecules.
    SEG_PARITY = 9
    SEG_PERSISTENT_COUNT = 10

    SEG_RING_BITS = 10
    SEG_RELEVANT_RINGS = 11
    SEG_FEATURES = 12
    SEG_ELEMENT_INDEX = 13
    SEG_EDGE_WORD = 14
    SEG_COMPONENT_LABEL = 15
    SEG_STEREO_UNIT = 16
    SEG_COUNT = 17

    # Table entries physically present in a buffer this build writes. Only the persistent block is
    # meaningful; entries SEG_PERSISTENT_COUNT..SEG_TABLE_MAX-1 are written zero and reserved for
    # persistent segments a later release adds. Growth is a `seg_count` bump plus, once this many
    # are spent, a larger SEG_TABLE_MAX -- and a larger SEG_TABLE_MAX does NOT invalidate old
    # buffers, because every reader takes the table size from the buffer's own `seg_count`.
    #
    # Thirteen is not arbitrary: 24 fixed bytes + 13 * 8 = 128, which is v3's header size, so v4
    # segment payloads begin at exactly the offset v3's did.  `SEG_CONFORMERS` and `SEG_PARITY` have
    # spent two of the five spare entries and three remain.
    SEG_TABLE_MAX = 13

    SEG_MASK_XY = 1
    SEG_MASK_STEREO = 2
    SEG_MASK_PARITY = 4

    STRUCT_MAGIC = 0x43485933
    STRUCT_VERSION = 6
    STRUCT_VERSION_V5 = 5
    STRUCT_VERSION_V4 = 4
    STRUCT_VERSION_V3 = 3
    V3_HEADER_LEN = 128
    V3_SEG_COUNT = 13
    V3_PERSISTENT_COUNT = 5
    V4_PERSISTENT_COUNT = 9
    FLAG_WIDE_INDEX = 1
    FLAG_TOPOLOGY_DIRTY = 2

    HE_IN_RING = 1
    HE_AROMATIC = 2

    # 0xFFFF is "none" in both of sgroup_t's u16 reference fields. `ext_index` must round-trip a
    # vendor file's Sgroup number VERBATIM including zero -- nothing in the CTfile spec forbids
    # numbering from 0, and V2000 `M  STY` writes it in a field where `  0` is representable -- so
    # zero cannot double as a sentinel without making "child of group 0" and "no parent" the same
    # bits. Spending one value out of 65536 is the cheaper trade, and it also gives a record built
    # programmatically (by a reactor, a standardiser, a writer) a way to say it has no file number
    # rather than inventing one.
    SGROUP_NO_INDEX = 0xFFFF
    # AND THE BOUND IS 0xFFFE, NAMED BEFORE ANY VALIDATOR EXISTS TO GET IT WRONG.  A u16 holds
    # 0..0xFFFF; the greatest REAL index is one less, because the top value is spent above.  This is
    # the same shape as H_UNKNOWN against H_IMPLICIT_MAX (see the hydrogen block below), and that one
    # was diagnosed only AFTER a reader had derived its bound from the nibble's width and admitted the
    # sentinel as a count -- on three write paths, the worst of them a valence clamp that put a
    # computed number onto the value meaning "nobody could compute it".  The family has now cost three
    # commits across two epics, so this time the bound is written first: any validator for `index`,
    # `ext_index` or `parent` cites SGROUP_INDEX_MAX, and NEVER the field's width.  Asserted as the
    # RELATION rather than as two numbers, because two literals are two things that can drift.
    SGROUP_INDEX_MAX = 0xFFFE

    # THE SAME ARGUMENT ONE FIELD WIDER.  `conformer_t.ext_index` round-trips a file's own model
    # number verbatim, and a PDB `MODEL 0` is representable -- nothing in the spec forbids numbering
    # from zero -- so zero cannot double as "this conformer has no file number".  The top u32 value is
    # spent instead, exactly as SGROUP_NO_INDEX spends the top u16.
    CONF_NO_INDEX = 0xFFFFFFFF
    # AND THIS IS A DIFFERENT KIND OF NUMBER, WHICH IS WHY IT IS NOT `CONF_NO_INDEX - 1`.  It bounds
    # how many models one molecule may hold.  The count's field is a u32 and a u32 admits four billion
    # models; `count * atom_count * sizeof(xyz_t)` overruns the 4 GiB buffer limit long before that,
    # so THE WIDTH IS NOT THE BOUND (§6.3) and deriving one from it would be the H_UNKNOWN family's
    # mistake in a new field.  65535 is large enough that no chemistry file reaches it and small
    # enough that the worst case stays inside the buffer limit: 65535 models of a 660-atom chain is
    # 519 MB, of which the record table is 262 kB.  Exported, so a format module refusing a
    # trajectory cites this rather than guessing.
    CONF_MAX_MODELS = 0xFFFF
    # THE WIDEST INDEX A FILE MAY STATE, one below the sentinel -- the same shape SGROUP_INDEX_MAX has
    # above.  A validator needs the bound named before it can cite it, and a file stating 0xFFFFFFFF
    # would otherwise be stored as "no number" rather than refused.
    CONF_EXT_INDEX_MAX = 0xFFFFFFFE
    # THE CONFORMER RECORD BEFORE THIS VERSION, in bytes.  A version-4 or version-5 buffer's record is
    # four words wide where this build's is one, and `structure_from_bytes` needs the old stride twice:
    # to check the segment's declared length against the buffer's own version, and to walk the table
    # while re-laying it.  Named because a literal 16 at either site would be a second copy of a layout
    # the struct no longer states.
    CONFORMER_RECORD_V5 = 16

    # AND THE OTHER THREE u16 COUNTS HAVE NO SENTINEL, SO THEIR MAXIMUM IS THE FULL WIDTH.  This is
    # the same family read the other way, and it is written down for the same reason: the next reader
    # of this block has just been told three times that a u16 field's real bound is 0xFFFE, and
    # `data_len`, `fields_len` and `log_len` do not spend a value, so 0xFFFF is a legal count for all
    # three.  A validator that "helpfully" capped them at SGROUP_INDEX_MAX would refuse one legal
    # record in 65536 for a symmetry that does not exist.  A record exceeding this is REFUSED rather
    # than truncated -- see `sgroup_check_lists`; truncating is the silent loss the whole segment exists
    # to prevent, and refusing at the boundary is the standing posture for input we cannot store.
    SGROUP_LIST_MAX = 0xFFFF

    # `sgroup_t.flags`.  Two bits, and the second one is the interesting one.
    #
    # SGROUP_FLAG_DISP -- `disp` holds a FIELDDISP anchor.  Needed because (0, 0) is a legal anchor
    # and "no anchor" must not be spelled the same way; the empty-versus-stated-zero rule (section 6.3
    # of RULES.md) applied to a coordinate pair.
    #
    # SGROUP_FLAG_ALIAS -- the record is NOT a CTfile S-group but a V2000 `A  <n>` / MRV mrvAlias atom
    # display label: `atoms_len == 1` and blob handle `strings_off` is the label.  It rides in this
    # segment ON PURPOSE.  An alias is stable-id-keyed annotation whose only storage requirement is
    # exactly the S-group's -- survive a remap, be dropped AND REPORTED when its atom dies -- so
    # giving it its own segment would duplicate `structure_carry_sgroups` to gain nothing but a
    # second thing to forget.  `sgroups` filters these out and `aliases` materialises them, so no
    # caller sees the sharing.
    SGROUP_FLAG_DISP = 1
    SGROUP_FLAG_ALIAS = 2
    SGROUP_FLAG_DEFINED = 3


# Field domains for atom_t — one declaration each; every validator must use these constants.
# Placed above the struct so the field comments can reference them backward.
#
# THESE ARE STORAGE BOUNDS AND NOT CHEMICAL DOMAINS.  A validator asks "can the field hold this",
# never "is this chemistry"; the arena stores garbage faithfully by design.  Named here because
# CHARGE_MIN/CHARGE_MAX sit close enough to a future valence table to be mistaken for its domain,
# and they are not: the valence rule collection is 1036 rows over 118 elements and their charges run
# -4..+4 with nothing above +4.  A table keyed on charge must therefore answer "no rule" for +5..+8
# rather than a confident zero -- an absent row is a GAP, not a violation -- and any range a comment
# names must say WHICH range it is.
DEF CHARGE_MIN     = -4      # atom_t.charge: what the FIELD may hold, and nothing more
DEF CHARGE_MAX     =  8      # NOT the domain of any chemical rule -- see below
DEF ISOTOPE_MAX    = 65535   # atom_t.isotope  (uint16 max; 0 = unset)
DEF MAP_NUMBER_MAX = 9999    # atom_t.map_number
DEF H_NIBBLE_MAX   = 15      # atom_t.hydrogens: the 4-bit NIBBLE'S WIDTH, not a count's bound
DEF H_EXPLICIT_MAX = 15      # atom_t.hydrogens high nibble: a count, and it may reach the width
DEF H_IMPLICIT_MAX = 14      # atom_t.hydrogens low nibble: a count STOPS AT 14, because ...
DEF H_UNKNOWN      = 15      # ... 15 is not a count there, it is "nobody knows"
# THREE CONSTANTS FOR TWO NIBBLES, AND NOT ONE FOR ALL OF IT.  H_NIBBLE_MAX is a fact about the
# LAYOUT -- four bits hold 0..15 -- and it is not the bound of anything a caller may state.  The
# explicit nibble's count bound happens to equal it and is named separately anyway; the implicit
# nibble's does not, because H_UNKNOWN takes the top value.  A validator that bounds an implicit
# count by H_NIBBLE_MAX therefore accepts the sentinel AS A COUNT, which destroys the only property
# the third state has: that it stays distinguishable from a real answer.  Same defect as a valence
# table returning a confident 0 where it has no row.  So no bound is spelled `H_NIBBLE_MAX` outside
# this block, and where a bound is named it says WHICH range it bounds.
#
# THE IMPLICIT NIBBLE HAS A SENTINEL AND THE EXPLICIT ONE DOES NOT, which is not an asymmetry for its
# own sake.  An explicit hydrogen is an atom someone drew: it is there or it is not, and the count is
# a fact about the record.  An implicit count is a DERIVED number, and there are records for which no
# derivation exists -- an aromatic atom in a ring that will not kekulise, a `[B-]` at twelve bonds, a
# copper at an odd coordination number.  Storing 0 for those and publishing the exception on the side
# is a storage gap wearing a parameter's clothes.
#
# HOW COMMON IS IT: RARE, AND THE JUSTIFICATION DELIBERATELY DOES NOT REST ON THE COUNT.  The reader
# measurements disagree by two orders of magnitude -- "no valence rule" on public NCI 5K has been
# reported as both 155 atoms of 82,157 and 4 atoms of 5,012 molecules, because the two counts measure
# different things: a table miss is not the same event as an atom that ends up with no answer, and a
# reader may miss the table and still derive a count another way.  THE ZERO IS NOW MEASURED, not
# quoted: the SMILES reader epic swept 5,536 real molecules (4,990 NCI plus a 546-molecule SDF) and
# found 0 atoms stored with the sentinel, WITH A NEGATIVE CONTROL FIRST -- the probe `N(F)(F)(F)F`
# comes back 1, so the instrument was shown able to answer non-zero before its zero was believed
# (RULES.md section 6.3).  The transferable warning: THE PROBE MOLECULE MUST BE UNBRACKETED.  A bracket
# STATES its hydrogen count, so `[N](F)(F)(F)F` exercises the path that keeps the caller's 0 and never
# reaches the derivation that has no answer -- testing the wrong door and getting a zero that means
# nothing.  An unbracketed atom with no valence rule stores H_UNKNOWN and logs it.
#
# So the population is empty today and that is a fact rather than a report.  It is still not an
# argument against the state.  A representation must be able to hold "no answer" for the same reason a
# validator must not bound a count by a nibble's width: the population being empty today is a fact
# about today's corpus, and the alternative is fabricating a 0 that no later consumer can tell from
# a measurement.  So no comment or test in this file quotes a population as the reason -- if one
# needs a number, it names the corpus and the reader that produced it, and both change.
#
# CORRESPONDENCE WITH chython 2: `Element.implicit_hydrogens` returns `None` for exactly this
# population, and that None is the oracle this sentinel reproduces.  So `implicit_h_of` answers None
# too, and `int(None)` raising in a caller's arithmetic is the point: a fabricated 0 would not raise,
# it would just make the mass, the formula and the valence quietly wrong.  15 rather than a separate
# flag byte because the nibble already had a value nobody could reach -- no atom carries fifteen
# implicit hydrogens -- so the state costs no space and cannot be lost by a writer that forgets it.


# CIP DESCRIPTOR STORAGE.  Storage and nothing else: no descriptor is computed here, or anywhere in
# this release.  A stored descriptor got here because an input STATED it, and that is the whole reason
# the arena holds one.
#
# THE VALUES ARE DEFINED BY POSITION IN `ATOM_CIP_CODES` / `BOND_CIP_CODES` (`_molecule_container.pxi`),
# not by a constant per descriptor, and that is deliberate.  Two of the eight atom spellings and two of
# the four bond spellings are the SAME LETTERS -- M and P are axial descriptors on both sides -- so a
# per-value constant family would need either two names for one letter or one name shared across two
# fields of different widths.  Both are how a reverse mapping ends up disagreeing with a forward one.
# One ordered table per field, index IS the stored code, and no second spelling exists to drift.
#
# CASE IS SIGNIFICANT AND NOTHING HERE UPPERCASES.  Lowercase r/s are the pseudo-asymmetric
# descriptors of CIP's auxiliary rules -- a different determination about a different kind of centre --
# so 'r' is not a spelling of 'R' and a `.upper()` anywhere on this path destroys a distinction the
# input made.  The display forms of some file formats have no case, and that uppercasing belongs at
# those write boundaries, never here.
DEF ATOM_CIP_MAX = 8         # codes 1..8 = R S r s M P m p; 0 = no descriptor
DEF BOND_CIP_MAX = 4         # codes 1..4 = E Z M P;         0 = no descriptor

# `atom_t.reserved` STOPS BEING "MUST STAY ZERO" HERE, so it needs the discipline `halfedge_t.flags`
# already has: a mask of what is DEFINED, with everything above it rejected by `structure_from_bytes`.
# See the note on HE_FLAG_DEFINED for what happens to a reserved region nobody validates -- the top
# fourteen half-edge flag bits were dead payload no reader could see, found by a byte-flip sweep.  The
# low nibble (bits 0-3) is the atom's CIP code, bits 4-11 the R index, bits 12-31 the CGR / Query
# hook, required to be zero.
DEF ATOM_CIP_MASK        = 0x0000000f
# Bits 4-11 of `reserved`: the R index of an R atom (element 0), 0 when it has none.  An index without
# element 0 is refused on load -- the two fields are one fact and a mismatch is a corrupt record.
DEF ATOM_R_INDEX_MASK  = 0x00000ff0
DEF ATOM_R_INDEX_SHIFT = 4
DEF R_INDEX_MAX = 99         # atom_t.reserved bits 4-11 hold 0..255; the DOMAIN is two decimal digits,
                             # so an index fits a three-character CTfile symbol column and a two-
                             # character depiction label.  A stored 100-255 is representable and
                             # illegal -- `structure_from_bytes` rejects it.
DEF ATOM_RESERVED_DEFINED = 0x00000fff

# The bond's code in `halfedge_t.flags`, above HE_IN_RING and HE_AROMATIC.  Three bits for four
# values, so code 5 is reachable in the field and refused by the validator -- the width is not the
# domain, exactly as for the hydrogen nibbles above.
DEF HE_CIP_SHIFT = 2
DEF HE_CIP_MASK  = 0x001c
DEF ATOM_FLAGS_RESERVED = 0x82   # bits 1 and 7


cdef packed struct atom_t:
    uint8_t element        # atomic number 1-118
    int8_t charge          # CHARGE_MIN..CHARGE_MAX
    uint8_t hydrogens      # low nibble implicit 0..H_IMPLICIT_MAX or H_UNKNOWN, high nibble explicit 0..H_EXPLICIT_MAX
    uint8_t flags          # 0 radical, 1 RESERVED, 2 in_ring, 3-5 hybridization, 6 h_pinned, 7 RESERVED
                           # `ATOM_FLAGS_RESERVED` bits are reserved; every writer leaves them 0, and
                           # a version-5 buffer that sets one is refused.  In version 3 and version 4
                           # these bits ARE the atom's parity -- `MoleculeContainer.from_bytes` adopts
                           # them.  Widening this byte is a layout change: `to_bytes()` is a molecule
                           # identity, so every stored key is repriced.
    uint16_t isotope       # 0..ISOTOPE_MAX (absolute mass number; 0 = unset)
    uint16_t map_number    # 0..MAP_NUMBER_MAX (AAM; 0 = unset)
    uint32_t n             # the atom's number; never reused
    uint8_t degree         # derived, inlined
    uint8_t heteroatoms    # derived, inlined
    uint32_t ring_sizes    # derived, inlined: sizes 3-24 plus 3 buckets
    uint16_t ring_counts   # low byte total, high byte aromatic
    uint32_t reserved      # bits 0-3 CIP code (ATOM_CIP_MASK), bits 4-11 R index (ATOM_R_INDEX_MASK),
                           # bits 12-31 CGR / Query hook, must stay zero.  NOT a free word any more:
                           # see ATOM_RESERVED_DEFINED, which is what `structure_from_bytes` rejects a
                           # stray bit against.


cdef packed struct halfedge_t:
    uint32_t to
    uint8_t order          # 1, 2, 3, 4 (aromatic), 8
    uint8_t wedge          # 0 none, 1 up, 2 down, 3 either
    uint16_t flags         # bit 0 in_ring, bit 1 aromatic, bits 2-4 CIP code, bits 5-15 spare


cdef struct edge_edit_t:
    uint32_t src
    uint32_t dst
    uint8_t order


# Bond orders a caller may state, AND the orders the arena stores.  One set for both.
#
# INPUT FIDELITY IS THE INVARIANT.  A source that says aromatic is stored aromatic; a source that
# says Kekule is stored Kekule; nothing in this layer normalises either way.  Order 4 is stored as
# order 4 with HE_AROMATIC set, and a molecule may hold one type-4 ring and one alternating ring at
# the same time -- that is a faithful record of two differently-written inputs, not a defect.
# Kekulising on apply would be a silent normalisation of the caller's representation, and only
# `kekule()` and `thiele()` may change which representation a molecule holds.
#
# NO SURFACE HERE NEEDS A GATE ON A STORED AROMATIC BOND.  Ring perception reads topology and never
# an order; hybridization answers 4 because 4 is the domain table's reserved value for aromatic; the
# canonical bond word folds HE_AROMATIC in as a distinct value; and stereo takes a five-line
# electron-budget correction (`_stereo.pxi`, `spent`) instead.  So `is_kekule` and
# `aromatic_bond_count` are the whole surface, and that is the point rather than a shortfall: they
# let a caller SEE which representation it holds.  No refusal helper lives at the bottom of this
# file -- an unused refusal reads as a policy the core has, and this one is not one.
cdef frozenset ALLOWED_ORDERS = frozenset((1, 2, 3, 4, 8))
cdef str ALLOWED_ORDERS_MSG = 'order must be 1, 2, 3, 4 (aromatic) or 8'

# The half-edge flag bits that are DEFINED.  Everything above them is reserved and must be zero;
# `from_bytes` rejects a buffer that sets one, which is what keeps a future bit free to mean
# something.  Without this check the top fourteen bits of every half-edge were dead payload: a
# byte-flip sweep over the frozen v3 records found them to be the only non-padding bytes in a
# persistent segment that no reader could see.
DEF HE_FLAG_DEFINED = 0x1f         # HE_IN_RING | HE_AROMATIC | HE_CIP_MASK

# "No atom index" in SEG_SGROUP_INDEX.  Only a CSTATE pair uses it: a CSTATE whose bond index did not
# resolve on read keeps its whole value as opaque text, and the pair slots have to say that rather
# than name atom 0.  Atom indices are dense from 0, so the top value is the only free one, and by the
# rule in the S-group block above the greatest storable index is therefore 0xFFFFFFFE -- not a real
# constraint (the 4 GiB buffer limit bites four billion atoms earlier) but named because the family's
# whole lesson is that an unnamed bound gets re-derived from the field's width.
#
# NEITHER AN ENUM MEMBER NOR A `DEF`, AND BOTH REASONS ARE MEASURED.  This one value cost two failed
# spellings before the third, and the failures are opposite kinds:
#
# * In the S-group `cdef enum` beside SGROUP_NO_INDEX it does not just widen itself.  0xFFFFFFFF does
#   not fit an `int`, so the C compiler retypes THE WHOLE ENUM as `unsigned int`, and every comparison
#   anywhere in the translation unit between an `int` local and ANY member of that enum becomes a
#   signed/unsigned mismatch.  It produced five `-Wsign-compare` warnings in
#   `_structure_resolve_persistent` and `structure_clone` -- code this change never touched, comparing
#   `int seg` against SEG_* ids that happen to share the enum.  One member's value is a property of
#   every other member, which is a coupling no comment inside the enum would make visible.
# * As a `DEF` it stops being a C constant at all: Cython types an out-of-`int`-range literal as a
#   Python object, so `src[k] == SGROUP_NO_REF` inside `_sgroup_carry_pairs` needs the GIL and the
#   `nogil` carry will not compile.  Four errors, all on the one line.
#
# A module-level `cdef uint32_t` is neither, and it is already the house spelling for exactly this
# value -- MATCH_UNSET/MATCH_DONE (_isomorphism.pxi), AROM_NO_ATOM (_kekule.pxi) and SU_NO_REF
# (_stereo.pxi) all carry the same note.  SHARING SU_NO_REF'S BIT PATTERN IS DELIBERATE AND SAFE
# HERE, where MATCH_UNSET's is not: nothing ever compares an S-group index slot against a stereo
# unit's `refs`, so the two sentinels have no meeting point to be confused at.
cdef uint32_t SGROUP_NO_REF = 0xFFFFFFFF


# XY_SCALE: the arena stores display coordinates as a fixed-point int32 (molecule units multiplied
# by XY_SCALE).  Every write multiplies by XY_SCALE; every read divides by it.  Declared here,
# beside xy_t, because §6.1 of RULES.md requires "every field's domain declared exactly once" and
# this file owns the struct and its read accessors.
#
# A site that wants the raw integer rather than the float — e.g. a serialiser passing the stored
# value through unchanged — still cites XY_SCALE and carries a one-line comment saying so.
DEF XY_SCALE = 10000

cdef packed struct xy_t:
    int32_t x
    int32_t y


# A CONFORMER'S POSITION, AND `xy_t` IS NOT ITS PREFIX BY ACCIDENT.  Same fixed-point grid, same
# XY_SCALE, so the x,y a 3D file states reach `SEG_XY` and `SEG_CONFORMERS` as the SAME bits and a
# projection taken from a conformer is byte-identical to the layout the file produced.
#
# A SEPARATE STRUCT AND NOT A WIDENED `xy_t`, which is the crystals design's Q1 and the one decision
# here that a later author is most likely to want to undo.  Widening is the smaller edit and costs
# every purely 2D molecule four bytes per atom -- and `to_bytes()` is a molecule identity, so that
# reprices every stored key for a molecule that has no third coordinate at all.  A segment costs
# nothing to a molecule that does not use it.
#
# RULES.md §1.4 applies and is the trap: the slot's range is +-214748.0, which is WIDER than any
# CTfile `F10.4` field's, so `set_xyz` range-checks the SLOT and nothing else.  A value that will not
# fit in ten columns is the writer's problem, at write time, with a log line -- not a refusal here.
cdef packed struct xyz_t:
    int32_t x
    int32_t y
    int32_t z


# WHAT A CONFORMER KNOWS ABOUT ITSELF: the number the file gave it, and nothing else.
#
# NO PROPERTY IS STORED, because singling one out is arbitrary -- extXYZ's comment line, an SDF
# property and a MOL2 remark each state a different set and none of them states a unit.  A value
# belongs to whatever names it, and this record names coordinates.
#
# What is NOT here, each for a reason rather than a schedule: a provenance string and a conformer
# name (both free text, and the arena's one string segment exists under the invariant that nobody
# decodes it), and per-conformer velocities, forces or charges (no chython consumer; a segment is
# cheaper to append later than to shrink).
cdef packed struct conformer_t:
    uint32_t ext_index      # the file's own MODEL/frame number, VERBATIM; CONF_NO_INDEX = none


# The segment's own header, so that `count` travels with the payload rather than being derived from
# the segment's length.  Both would work; this way `structure_from_bytes` can CHECK one against the
# other, and a stated length disagreeing with the derived one is a corrupt buffer rather than a
# silently different molecule.
cdef packed struct conformer_hdr_t:
    uint32_t count          # models present, 1..CONF_MAX_MODELS; 0 is not written at all
    uint32_t reserved       # must stay zero


cdef packed struct segment_t:
    uint32_t offset
    uint32_t length


# ONE LEVEL DOWN FROM `segment_t`, AND DELIBERATELY THE SAME TWO FIELDS.  A segment is already
# (offset, length), so byte-level variability was never the missing piece -- RECORD FRAMING was: N
# byte strings of differing size, addressed in O(1) by their position in the list (a `handle`).  The
# blob's on-disk layout is
#
#     off 0        uint32 count
#     off 4        uint32 payload_off        == 8 + count * 8, so always 8-aligned
#     off 8        blob_rec_t rec[count]     off is relative to payload_off and 8-ALIGNED
#     payload_off  the bytes, each record 8-aligned, inter-record padding zero
#
# `off` is 8-aligned so a record can be read as `uint32 *` (S-group reference runs) with no unaligned
# load; `len` is stored EXACT because the alignment padding would otherwise destroy it, and storing
# both is smaller and plainer than a `count + 1` offset array plus a per-record length prefix.
cdef packed struct blob_rec_t:
    uint32_t off
    uint32_t len


# One CTfile S-group, or one atom alias -- see SGROUP_FLAG_ALIAS.  48 bytes, every field naturally
# aligned within an 8-aligned stride, which is what lets `disp` and the four u32 counts be read
# without an unaligned load.
#
# WHY THIS IS TWO SEGMENTS AND NOT ONE.  Everything a molecule EDIT can invalidate is an atom
# reference, and every atom reference in a record lives in one contiguous run of `SEG_SGROUP_INDEX`
# starting at `refs_off`, in the fixed order atoms | patoms | bonds | cstates.  So `_apply`'s carry
# is one loop over one array with one map, and the opaque half (`SEG_OPAQUE_BLOB`) is copied byte for
# byte by everything, forever.  The split is not structured-versus-unstructured; it is
# INVALIDATABLE-versus-OPAQUE, and that is why `type` -- a string, structured by any normal reading --
# sits on the opaque side.
#
# EVERY `*_len` IS A COUNT OF SLOTS OR HANDLES, never of the things they encode.  `bonds_len` and
# `cstates_len` are EVEN: two slots per endpoint pair.  `fields_len` is EVEN: two handles per
# (keyword, value).  Counting pairs instead would put the same number in two units in one struct,
# which is how a doubling or a halving gets written.
#
# THERE IS NO `cstate_tails_len`, AND ITS ABSENCE IS THE DESIGN.  A CSTATE is a bond reference plus a
# vector tail that this library has no opinion about, so each pair needs one string beside it --
# `cstates_len // 2` of them, a number this record already states.  Giving it a field would let the
# two counts disagree, and a disagreement there does not corrupt bytes, it RE-PAIRS TAILS WITH THE
# WRONG BONDS, which is silent and is the exact failure the structured CSTATE model exists to prevent.
# So the tails are a run whose length is DERIVED, and the only way to change how many there are is to
# change how many pairs there are.
cdef packed struct sgroup_t:
    uint32_t refs_off       # first slot of this record's run in SEG_SGROUP_INDEX
    uint32_t strings_off    # first handle of this record's run in SEG_OPAQUE_BLOB
    xy_t disp               # FIELDDISP anchor; MEANINGFUL ONLY under SGROUP_FLAG_DISP
    uint32_t atoms_len      # atom indices
    uint32_t patoms_len     # atom indices -- the parent-atom subset of a MUL group
    uint32_t bonds_len      # atom indices, two per endpoint pair
    uint32_t cstates_len    # atom indices, two per pair; (NO_REF, NO_REF) = never resolved
    uint16_t index          # the file's Sgroup number, VERBATIM; SGROUP_NO_INDEX = unnumbered
    uint16_t ext_index      # V2000's external number; SGROUP_NO_INDEX = none
    uint16_t parent         # another record's `index`; SGROUP_NO_INDEX = none
    uint16_t flags          # SGROUP_FLAG_*
    uint16_t data_len       # FIELDDATA handles
    uint16_t fields_len     # unmodelled keywords, two handles per (keyword, value)
    uint16_t log_len        # reader diagnostics carried with the record
    uint16_t spare          # must stay zero; from_bytes rejects non-zero


cdef packed struct StructureHeader:
    uint32_t magic
    uint16_t version
    uint16_t flags
    uint32_t atom_count
    uint32_t bond_count
    uint32_t persistent_len   # == the whole buffer; the buffer is never grown
    uint16_t seg_count        # table entries present in THIS buffer; the table size is data
    uint16_t reserved0        # must stay zero
    segment_t segments[SEG_TABLE_MAX]
    # THE STRUCT IS 128 BYTES AND THE HEADER USUALLY IS NOT.  `segments` is sized SEG_TABLE_MAX so
    # that C can name every entry, but `header_len` is 24 + 8 * seg_count, and the writer stops the
    # table one past the highest entry the molecule uses: 48 bytes for a molecule with no
    # coordinates, no stereo groups and no S-groups, which is most of them.  Anything that reads
    # `segments[i]` for `i >= seg_count` is reading the FIRST SEGMENT'S PAYLOAD -- atom records --
    # and anything that WRITES there corrupts them.  `_structure_resolve_persistent`,
    # `structure_from_bytes` and `_alloc_probe` are the three readers, and each bounds itself; a
    # fourth must do the same.  Trailing entries only: ids are positional, so an interior empty
    # (a text blob with no coordinates) still occupies its slot.
    #
    # 24 fixed bytes + 13 * 8 = 128 is therefore the MAXIMUM header, reached only by a molecule that
    # uses the last persistent segment -- and by every v3 buffer, whose thirteen entries are why the
    # two versions' payloads line up at all.
    #
    # Two things about this layout are load-bearing, and both are why v4 costs a version byte and
    # nothing else:
    #
    # 1. THE TABLE STARTS AT OFFSET 24 IN v3 AND IN v4, and persistent ids 0-4 mean the same thing
    #    in both. v4 bought `seg_count` by DELETING v3's `total_len` -- the field that made a read
    #    change the serialised bytes -- so the four bytes came from the defect rather than from the
    #    table. With seg_count == 13 the header is 128 bytes, exactly v3's, so a v3 buffer's segment
    #    payloads are already where a v4 reader looks for them and nothing is relocated on ingest.
    #    What a v3 buffer does NOT agree about is bytes 20-23, which hold its `total_len`, and table
    #    entries 5-11, which name its DERIVED segments -- and under v4 ids 5-7 are persistent
    #    S-group segments. So the v3 read path must supply seg_count itself and zero every entry at
    #    or above SEG_PERSISTENT_COUNT; a warm v3 buffer otherwise looks like it carries S-groups at
    #    meaningless offsets. `structure_from_bytes` does both, and the compatibility suite runs on
    #    frozen warm v3 bytes precisely to keep it doing them.
    #
    # 2. GROWTH IS A `seg_count` BUMP. Only persistent segments are in this table at all -- derived
    #    caches live in the `Structure` object, where their ids are compile-time and renumbering them
    #    is free -- so the ten persistent slots in use leave three spare before SEG_TABLE_MAX must
    #    rise, and raising it invalidates nothing, because a reader sizes the table from the buffer's
    #    own `seg_count` and not from this constant. A reader accepts any `seg_count`; an entry it
    #    does not know must be EMPTY, and a non-empty unknown segment is an error naming its index,
    #    so an old build handed a molecule carrying data it cannot model says so instead of quietly
    #    writing the molecule back without it.
    #
    # `to_bytes()` IS a molecule identity: `persistent_len` equals the buffer length, no read writes
    # inside it, and a hash or a dedup key is taken over the whole slice. Nothing has to be skipped at
    # the front of the buffer.


# A segment whose element stride times its count can exceed ZERO_PAGE_SIZE must be
# guarded with structure_has() before reading; the zero page only covers 4096 bytes.
DEF ZERO_PAGE_SIZE = 4096

cdef char _zero_page[ZERO_PAGE_SIZE]


# SEG_COUNT - SEG_PERSISTENT_COUNT, spelled as a literal because a `DEF` is evaluated by the Cython
# compiler and cannot read an enumerator.  So it is a duplicate of two numbers thirty lines above,
# which is exactly the kind of thing that goes stale when a segment is added -- hence the module-init
# check next to `globals().update(_segment_ids())`, which refuses to import a build where the three
# numbers disagree.  A wrong value here would silently under-size `_derived` and `_retired` and the
# arena would write one pointer past the end of both.
DEF SEG_DERIVED_COUNT = 7


cdef class Structure:
    """One arena: a persistent buffer that never moves, plus derived caches beside it.

    THE PERSISTENT BUFFER IS ALLOCATED ONCE AND NEVER REALLOCATED. Every derived segment gets its
    own allocation, recorded here rather than in the serialised header. That single property is what
    retires ruling F60: appending segment X cannot move segment Y, because they are different
    allocations, so a pointer taken into the arena stays valid across `ensure_stereo_units`,
    `ensure_component_labels`, `rebuild_derived` and anything else that builds a cache. There is no
    ordering rule left to remember, and no invalidation comment left to honour.

    `_seg_base` is the resolved pointer for every segment, persistent and derived alike, and it is
    `_zero_page` for a segment that is absent. So `segment()` is one aligned load with no branch,
    and a loop over an absent segment needs no `structure_has` guard -- the property the ring bitmap
    and the stereo unit table already depended on, now uniform and cheaper.
    """
    cdef char *buffer                              # the persistent allocation; persistent_len bytes
    cdef StructureHeader *header
    cdef size_t buffer_len                         # bytes at `buffer`; always == persistent_len
    cdef size_t total_len                          # buffer_len + every derived allocation
    cdef bint owns
    cdef void *_seg_base[SEG_COUNT]
    cdef uint32_t _seg_len[SEG_COUNT]
    # Derived allocations this object owns and must free. `_retired` is one deep per derived segment:
    # when a cache is replaced, the old block is retired rather than freed, so a pointer taken before
    # the replacement reads the PREVIOUS CORRECT table instead of freed memory. That is byte for byte
    # the semantics v3 already had -- it left the old table stranded inside the buffer and said so --
    # except that here the block is tracked and released at __dealloc__ rather than leaked.
    cdef void *_derived[SEG_DERIVED_COUNT]
    cdef void *_retired[SEG_DERIVED_COUNT]
    # How many bonds are stored with order 4.  DERIVED, and therefore here and not in the header:
    # a header field would make `to_bytes()` depend on it, which is the defect this format version
    # exists to fix.  Off-header it is also a uint32 with no ceiling, where the only spare header
    # room was two bytes and would have capped silently at 65535 aromatic bonds.
    #
    # NOT MAINTAINED -- RECOMPUTED.  Every mutation goes through `_apply`, which rebuilds the CSR
    # from scratch, so there is no incremental update path that could get out of step with the bond
    # orders.  `csr_build` and `structure_from_bytes` each set it inside a loop they already run over
    # every half-edge, so it costs nothing and cannot disagree with what it counts.  That is why
    # there is no stored "representation" flag: a flag can contradict the bonds, a count cannot.
    cdef uint32_t aromatic_bond_count

    def __dealloc__(self):
        cdef int i
        for i in range(SEG_DERIVED_COUNT):
            if self._derived[i] is not NULL:
                PyMem_Free(self._derived[i])
                self._derived[i] = NULL
            if self._retired[i] is not NULL:
                PyMem_Free(self._retired[i])
                self._retired[i] = NULL
        if self.owns and self.buffer is not NULL:
            PyMem_Free(self.buffer)
            self.buffer = NULL

    cdef inline void *segment(self, int seg) noexcept nogil:
        return self._seg_base[seg]

    cdef inline atom_t *atoms(self) noexcept nogil:
        return <atom_t *> self._seg_base[SEG_ATOMS]


cdef inline void _structure_blank(Structure structure) noexcept nogil:
    """Point every segment at the zero page and own nothing. Call before any layout."""
    cdef int seg
    for seg in range(SEG_COUNT):
        structure._seg_base[seg] = <void *> _zero_page
        structure._seg_len[seg] = 0
    for seg in range(SEG_DERIVED_COUNT):
        structure._derived[seg] = NULL
        structure._retired[seg] = NULL


cdef inline void _structure_resolve_persistent(Structure structure) noexcept nogil:
    """Resolve the persistent table into `_seg_base` / `_seg_len`.

    The header's table is the truth on disk; these arrays are the truth in memory. Resolving once
    here is what turns `segment()` into a single load, and it is also the only place that reads an
    offset out of the packed header -- RULES.md 2.3 forbids binding a pointer into a packed struct,
    and hoisting the two values into locals is how this obeys it.

    STOPS AT `seg_count`, WHICH IS USUALLY LESS THAN SEG_PERSISTENT_COUNT.  Trailing empty entries
    are not written, so the bytes at the missing entries' positions are the first segment's PAYLOAD:
    reading them would resolve `_seg_base[SEG_XY]` to a pointer computed out of atom records.  The
    entries past the table keep `_structure_blank`'s zero page and a zero length, which is exactly
    what an absent segment reads as anyway -- so no caller can tell "absent" from "not in the table",
    and none should be able to.
    """
    cdef int seg
    cdef int present = structure.header.seg_count
    cdef uint32_t offset, length
    if present > SEG_PERSISTENT_COUNT:
        present = SEG_PERSISTENT_COUNT
    for seg in range(present):
        offset = structure.header.segments[seg].offset
        length = structure.header.segments[seg].length
        structure._seg_len[seg] = length
        if length == 0:
            structure._seg_base[seg] = <void *> _zero_page
        else:
            structure._seg_base[seg] = <void *> (structure.buffer + offset)
    for seg in range(present, SEG_PERSISTENT_COUNT):
        structure._seg_len[seg] = 0
        structure._seg_base[seg] = <void *> _zero_page


cdef inline xy_t *structure_xy(Structure structure) noexcept nogil:
    return <xy_t *> structure.segment(SEG_XY)


cdef inline double xy_read_x(xy_t *p) noexcept nogil:
    """ONE implementation of the x-coordinate read — divide the stored fixed-point integer by
    XY_SCALE.  Six callers, cited BY NAME because a line number goes stale on the next edit above
    it and a name does not: `MoleculeContainer.union`, `.xy_of` and `.coordinates` in
    `_molecule_container.pxi`, `Atom.x` and `Atom.y` in `_molecule_views.pxi`, and `_pach_encode`
    in `_pach.pxi`.  All go through here; the literal 10000 does not appear elsewhere."""
    return p.x / <double> XY_SCALE


cdef inline double xy_read_y(xy_t *p) noexcept nogil:
    """ONE implementation of the y-coordinate read.  Same six callers as xy_read_x."""
    return p.y / <double> XY_SCALE


# ------------------------------------------------------------------------------------------------
# SEG_CONFORMERS: 3D geometry, and the one place its layout arithmetic lives.
#
# The payload is one contiguous run:
#
#     offset 0                    conformer_hdr_t        count, reserved
#     offset 8                    conformer_t[count]     4 bytes each
#     offset 8 + 4 * count        xyz_t[count][atoms]    12 bytes each, MODEL-MAJOR
#
# MODEL-MAJOR, so one model's positions are contiguous.  That is the access pattern every consumer
# has -- write a model, project a model, superpose a model -- and atom-major would stride every one
# of them.  It also makes a single model's coordinates one `memcpy`.
#
# `count` TRAVELS WITH THE PAYLOAD even though the segment's length would yield it.  Both work; this
# way `structure_from_bytes` can check one against the other, so a stated length that disagrees with
# the derived one is a corrupt buffer rather than a silently different molecule.  An empty segment is
# not written at all -- a zero-length segment is ABSENT (`structure_has` says so), so "count == 0"
# and "no segment" are the same state and there is no third one to keep consistent.

cdef inline size_t conformer_seg_len_for(uint32_t models, uint32_t atom_count,
                                         size_t record) noexcept nogil:
    """`conformer_seg_len` at a stated record width, for the one caller reading a buffer whose record
    is not this build's: `structure_from_bytes`, checking a version-4 or version-5 length and then
    walking its table.  `record` is `sizeof(conformer_t)` for every other caller."""
    if models == 0:
        return 0
    return align8(sizeof(conformer_hdr_t) + <size_t> models * record
                  + <size_t> models * atom_count * sizeof(xyz_t))


cdef inline size_t conformer_seg_len(uint32_t models, uint32_t atom_count) noexcept nogil:
    """The segment's byte length for `models` models of an `atom_count`-atom molecule. ONE COPY.

    Every other site -- the allocator, the copy, the ingest validator, the container -- calls this
    rather than repeating the arithmetic, because §6.1 makes a layout a field whose domain is
    declared exactly once and four copies of a formula are four things that can disagree about one
    buffer.  Returns 0 for 0 models, which is what makes "no conformers" cost nothing.

    `size_t` and not `uint32_t`: on a 64-bit build the product cannot overflow before the caller's
    own 4 GiB check sees it, and returning the narrow type would wrap first and pass that check.
    """
    return conformer_seg_len_for(models, atom_count, sizeof(conformer_t))


cdef inline uint32_t structure_conformer_count(Structure structure) noexcept nogil:
    """How many models this molecule holds. 0 when the segment is absent, via the zero page.

    Reads through `segment()` rather than the header table, so an absent segment resolves to the zero
    page and the count reads 0 with no branch -- the same trick `structure_xy` relies on.
    """
    return (<conformer_hdr_t *> structure.segment(SEG_CONFORMERS)).count


cdef inline conformer_t *structure_conformer_records(Structure structure) noexcept nogil:
    """The per-model record table. Only valid when the count is non-zero."""
    return <conformer_t *> (<char *> structure.segment(SEG_CONFORMERS) + sizeof(conformer_hdr_t))


cdef inline xyz_t *structure_conformer_xyz(Structure structure, uint32_t model) noexcept nogil:
    """Model `model`'s positions, atom-indexed. Only valid when `model < count`.

    THE CALLER BOUNDS `model` AND THIS DOES NOT, which is the house pattern for an arena accessor
    (`csr_edges`, `structure_xy` and `atoms()` all bound nothing either) -- a check here would run on
    every atom of every model and answer a question the caller has already answered.  The two entry
    points that take a model from outside the core, `set_xyz` and `xyz_of`, both range-check it.
    """
    cdef char *base = <char *> structure.segment(SEG_CONFORMERS)
    cdef uint32_t count = (<conformer_hdr_t *> base).count
    return <xyz_t *> (base + sizeof(conformer_hdr_t) + <size_t> count * sizeof(conformer_t)
                      + <size_t> model * structure.header.atom_count * sizeof(xyz_t))


# THREE READERS AND NOT A CAST TO `xy_t *`.  `xyz_t`'s first two fields ARE `xy_t`'s, deliberately, so
# reusing `xy_read_x` on an `xyz_t *` would compile and would give the right answer today -- and it
# would be type punning between two packed structs, which is the shape §2.3 exists to forbid and which
# stops being right the moment either struct gains a field.  The scale is still declared exactly once:
# all five readers divide by the DEF, and the literal 10000 appears in none of them.

cdef inline double xyz_read_x(xyz_t *p) noexcept nogil:
    return p.x / <double> XY_SCALE


cdef inline double xyz_read_y(xyz_t *p) noexcept nogil:
    return p.y / <double> XY_SCALE


cdef inline double xyz_read_z(xyz_t *p) noexcept nogil:
    return p.z / <double> XY_SCALE


cdef inline uint8_t *structure_stereo_groups(Structure structure) noexcept nogil:
    return <uint8_t *> structure.segment(SEG_STEREO_GROUPS)


#: The highest OR/AND group id one `SEG_STEREO_GROUPS` byte holds beside its two kind bits.  The
#: domain of a group id, declared here because `sg_pack` is what bounds it; `set_stereo_group` is the
#: one place that refuses, and both readers renumber a file's larger id rather than lose the group.
DEF STEREO_GROUP_MAX = 0x3f


cdef inline uint8_t sg_kind(uint8_t v) noexcept nogil:
    return v >> 6


cdef inline uint8_t sg_group(uint8_t v) noexcept nogil:
    return v & 0x3f


cdef inline uint8_t sg_pack(uint8_t kind, uint8_t group) noexcept nogil:
    return ((kind & 0x03) << 6) | (group & 0x3f)


cdef inline uint8_t *structure_parities(Structure structure) noexcept nogil:
    """The stated parities, one byte per atom: 0 none, 1 even, 2 odd.

    GUARD WITH `structure_has` BEFORE INDEXING THIS.  An absent segment resolves to the zero page,
    which is ZERO_PAGE_SIZE bytes and shorter than a 5000-atom molecule's row -- the same obligation
    `structure_stereo_groups` carries.  `structure_parity_at` is the guarded single read; this is the
    pointer a loop hoists once after establishing the segment is there.

    `structure_parity_at`, `structure_set_parity` and `structure_clear_parities` test
    `_seg_len[SEG_PARITY]` directly rather than calling `structure_has` because `structure_has` is
    declared further down in the same translation unit and is not visible at their call sites -- the
    inlined test is the identical check.
    """
    return <uint8_t *> structure.segment(SEG_PARITY)


cdef inline uint8_t structure_parity_at(Structure structure, uint32_t slot) noexcept nogil:
    """One atom's three-state parity: 0 none configured, 1 even, 2 odd.

    A wedge does not configure a parity (Ruling F54), so a wedge-drawn centre answers 0 here exactly
    as an undrawn one does -- `wedge_of` and `Bond.wedge` are what answer that question.
    """
    if structure._seg_len[SEG_PARITY] == 0:
        return 0
    return (<uint8_t *> structure._seg_base[SEG_PARITY])[slot]


cdef inline int structure_set_parity(Structure structure, uint32_t slot, uint8_t value) except -1:
    """State one atom's parity.  REFUSES AN ABSENT SEGMENT rather than growing one.

    The persistent block is laid out once and never reallocated, so a writer names the segment first:
    `SEG_MASK_PARITY` at `structure_alloc_full` for a reader building its own arena, `OP_WANT_PARITY`
    for one writing after a seal.  A writer that did neither would otherwise store into the shared
    zero page and corrupt every absent segment in the process, so this is a defect raising rather
    than input being rejected.
    """
    if structure._seg_len[SEG_PARITY] == 0:
        raise RuntimeError('this arena has no parity segment; SEG_MASK_PARITY or OP_WANT_PARITY names '
                           'one before a parity can be stated')
    (<uint8_t *> structure._seg_base[SEG_PARITY])[slot] = value
    return 0


cdef inline void structure_clear_parities(Structure structure) noexcept nogil:
    """Every stated parity gone.  A no-op when the segment is absent -- there is nothing to clear."""
    if structure._seg_len[SEG_PARITY]:
        memset(structure._seg_base[SEG_PARITY], 0, structure._seg_len[SEG_PARITY])


cdef inline uint32_t *csr_ptr(Structure structure) noexcept nogil:
    return <uint32_t *> structure.segment(SEG_CSR_PTR)


cdef inline halfedge_t *csr_edges(Structure structure) noexcept nogil:
    return <halfedge_t *> structure.segment(SEG_CSR_EDGE)


cdef inline halfedge_t *csr_find_at(uint32_t *ptr, halfedge_t *edges,
                                    uint32_t i, uint32_t j) noexcept nogil:
    """Linear scan of atom i's adjacency list for a half-edge leading to j.

    Returns a pointer into the edges array, or NULL when no such half-edge exists.  The
    pointer arithmetic `he - edges` gives the half-edge's index, which is parallel to the
    edge_words array -- callers that need the edge_word use this rather than a separate
    lookup.  Extracted from csr_find so the isomorphism kernel can call it with the cached
    csr_begin and edges pointers from matcher_t without needing a Structure object.
    """
    cdef uint32_t k
    for k in range(ptr[i], ptr[i + 1]):
        if edges[k].to == j:
            return &edges[k]
    return NULL


cdef inline halfedge_t *csr_find(Structure structure, uint32_t i, uint32_t j) noexcept nogil:
    return csr_find_at(csr_ptr(structure), csr_edges(structure), i, j)


cdef inline uint8_t at_implicit_h(atom_t *a) noexcept nogil:
    """The implicit nibble RAW, so H_UNKNOWN (15) comes back as 15 and not as 0.

    DELIBERATELY NOT MASKED TO A SAFE ZERO.  A reader that has not been taught about the sentinel then
    computes with 15, which is a visibly absurd hydrogen count that shows up in the first test; the
    same reader given a silent 0 computes a plausible wrong answer nobody sees.  Every arithmetic
    caller must ask `at_implicit_h_unknown` first -- `_features.pxi`, `_stereo.pxi` and `float()` do.
    """
    return a.hydrogens & 0x0f


cdef inline bint at_implicit_h_unknown(atom_t *a) noexcept nogil:
    """Is this atom's implicit hydrogen count unknown?  The one question to ask before arithmetic."""
    return (a.hydrogens & 0x0f) == H_UNKNOWN


cdef inline uint8_t at_explicit_h(atom_t *a) noexcept nogil:
    return a.hydrogens >> 4


cdef inline void at_set_h(atom_t *a, uint8_t implicit, uint8_t explicit) noexcept nogil:
    a.hydrogens = (implicit & 0x0f) | ((explicit & 0x0f) << 4)


# EVERY WRITER BELOW IS A BYTE-WIDE READ-MODIFY-WRITE, and that is an invariant about callers and not
# just an implementation note.  `|=`, `&=` and the hybridization store alike load `flags`, alter some
# bits and store the whole byte back, because four independent facts share one byte
# (radical, in_ring, hybridization, h_pinned).  Under `freethreading_compatible` that makes the bits a
# SINGLE memory location in the C11 sense: two threads
# setting two DIFFERENT flags on the same atom race and one update is lost, where two separate bytes
# would have been safe without any synchronisation.  The race is between any two flags, not only the
# pair that happens to appear in one expression.
#
# So: `flags` is written under the structure's own mutation discipline -- one writer per structure, and
# a journal applied by `_apply` -- and NEVER concurrently per-field.  Nothing here takes a lock, and
# adding per-flag atomics would cost the packing its whole advantage (measured: packed `atom_t` is 28%
# faster on a bandwidth-bound scan than one byte per flag, because a composite predicate is one load
# plus ALU work instead of three loads).  The discipline is the cheaper half of that trade, which is
# why it is stated rather than enforced -- an unstated concurrency invariant is an unexported domain,
# and RULES.md 6.3 says what happens to those at a boundary.
cdef inline bint at_radical(atom_t *a) noexcept nogil:
    return (a.flags & 0x01) != 0


cdef inline void at_set_radical(atom_t *a, bint value) noexcept nogil:
    if value:
        a.flags |= 0x01
    else:
        a.flags &= ~(<uint8_t> 0x01)


cdef inline bint at_in_ring(atom_t *a) noexcept nogil:
    return (a.flags & 0x04) != 0


cdef inline void at_set_in_ring(atom_t *a, bint value) noexcept nogil:
    if value:
        a.flags |= 0x04
    else:
        a.flags &= ~(<uint8_t> 0x04)


cdef inline bint at_h_pinned(atom_t *a) noexcept nogil:
    return (a.flags & 0x40) != 0


cdef inline void at_set_h_pinned(atom_t *a, bint value) noexcept nogil:
    if value:
        a.flags |= 0x40
    else:
        a.flags &= ~(<uint8_t> 0x40)


cdef inline uint8_t at_hybridization(atom_t *a) noexcept nogil:
    return (a.flags >> 3) & 0x07


cdef inline void at_set_hybridization(atom_t *a, uint8_t z) noexcept nogil:
    a.flags = (a.flags & 0xc7) | ((z & 0x07) << 3)


cdef inline uint8_t at_cip(atom_t *a) noexcept nogil:
    """The atom's stored CIP code, 0 for none.  See the ATOM_CIP_MAX block for what the codes are.

    IN `reserved` AND NOT IN `flags`, and the reason is a layout fact rather than a preference: a
    new flag in `flags` is a layout change -- `to_bytes()` is a molecule identity, so every stored key
    is repriced.  `reserved` is already in those bytes as zero, so an atom with no descriptor
    serialises to exactly the same bytes it did before this field existed.  Nothing was repriced.
    """
    return <uint8_t> (a.reserved & ATOM_CIP_MASK)


cdef inline void at_set_cip(atom_t *a, uint8_t code) noexcept nogil:
    """Read-modify-write of the low nibble (bits 0-3), preserving bits 4-11 (R index) and bits 12-31
    (CGR / Query hook).

    Not a whole-word store even though bits 4-31 may be zero today, because a whole-word store is the
    assumption that makes the first user of those bits lose its value to a descriptor write.
    """
    a.reserved = (a.reserved & ~(<uint32_t> ATOM_CIP_MASK)) | (code & ATOM_CIP_MASK)


cdef inline uint8_t at_r_index(atom_t *a) noexcept nogil:
    return <uint8_t> ((a.reserved & ATOM_R_INDEX_MASK) >> ATOM_R_INDEX_SHIFT)


cdef inline void at_set_r_index(atom_t *a, uint8_t index) noexcept nogil:
    a.reserved = (a.reserved & ~(<uint32_t> ATOM_R_INDEX_MASK)) | \
                 ((<uint32_t> index << ATOM_R_INDEX_SHIFT) & ATOM_R_INDEX_MASK)


cdef inline uint8_t he_cip(halfedge_t *e) noexcept nogil:
    return <uint8_t> ((e.flags & HE_CIP_MASK) >> HE_CIP_SHIFT)


cdef inline void he_set_cip(halfedge_t *e, uint8_t code) noexcept nogil:
    """ONE HALF-EDGE.  A bond is two, and a CIP descriptor is not directional, so both must be set or
    the answer depends on which end you ask from.

    Deliberately NOT given a both-halves signature here: this layer has no way to reach the twin
    without a CSR lookup, and a helper that took a Structure would let a caller set one half and think
    it had set the bond.  The single place that writes a bond descriptor is the re-apply loop in
    `_apply`, which calls this twice from one site -- the same shape as `_emit_half` for order, and for
    the same reason stated there: the two halves must agree, and writing them out at two call sites is
    how they drift apart.
    """
    e.flags = <uint16_t> ((e.flags & ~(<uint16_t> HE_CIP_MASK))
                          | ((code << HE_CIP_SHIFT) & HE_CIP_MASK))


cdef inline uint8_t at_ring_count(atom_t *a) noexcept nogil:
    return a.ring_counts & 0xff


cdef inline uint8_t at_aromatic_ring_count(atom_t *a) noexcept nogil:
    return a.ring_counts >> 8


cdef inline void at_set_ring_counts(atom_t *a, uint8_t total, uint8_t aromatic) noexcept nogil:
    a.ring_counts = total | (<uint16_t> aromatic << 8)


cdef inline uint32_t *structure_rings(Structure structure) noexcept nogil:
    return <uint32_t *> structure.segment(SEG_RELEVANT_RINGS)


cdef inline uint64_t *structure_ring_bits(Structure structure) noexcept nogil:
    return <uint64_t *> structure.segment(SEG_RING_BITS)


cdef inline uint32_t structure_ring_words(Structure structure) noexcept nogil:
    # Derived from the bitmap's own length, not from SEG_RELEVANT_RINGS: the bitmap carries one
    # bit per relevant-cycle prototype while SEG_RELEVANT_RINGS carries a minimum cycle basis,
    # and those two counts differ on almost every polycycle.
    cdef uint32_t n = structure.header.atom_count
    cdef uint32_t bits = structure_seg_len(structure, SEG_RING_BITS)
    if n == 0 or bits == 0:
        return 0
    return <uint32_t> (bits // (<size_t> n * sizeof(uint64_t)))


cdef inline bint structure_shares_ring(Structure structure, uint32_t a,
                                       uint32_t b) noexcept nogil:
    """Do atoms `a` and `b` sit on some common relevant-cycle prototype?

    `structure_ring_words() == 0` subsumes a `structure_has(SEG_RING_BITS)` guard of its own: an
    absent segment has length zero and the word count is derived from that length.

    Named on the structure rather than the atom pair because the bitmap is a segment, not an
    `atom_t` field.  Callers: `MoleculeContainer.shares_ring` and the stereo layer's small-ring cut.
    """
    cdef uint32_t words = structure_ring_words(structure)
    cdef uint64_t *bits
    cdef uint32_t k
    if words == 0:
        return False
    bits = structure_ring_bits(structure)
    for k in range(words):
        if bits[<size_t> a * words + k] & bits[<size_t> b * words + k]:
            return True
    return False


cdef inline void at_add_ring_size(atom_t *a, uint32_t size) noexcept nogil:
    if size < 3:
        return
    elif size <= 24:
        a.ring_sizes |= <uint32_t> 1 << size
    elif size <= 32:
        a.ring_sizes |= 1
    elif size <= 48:
        a.ring_sizes |= 2
    else:
        a.ring_sizes |= 4


cdef inline uint64_t *structure_features(Structure structure) noexcept nogil:
    return <uint64_t *> structure.segment(SEG_FEATURES)


cdef inline uint64_t *structure_edge_words(Structure structure) noexcept nogil:
    return <uint64_t *> structure.segment(SEG_EDGE_WORD)


cdef inline uint32_t *structure_component_labels(Structure structure) noexcept nogil:
    return <uint32_t *> structure.segment(SEG_COMPONENT_LABEL)


# SEG_STEREO_UNIT's accessors are in `_stereo.pxi`, beside the `stereo_unit_t` they return and the
# perception that fills them: unlike every segment above it, its records are not a flat array of a
# built-in type, and its length is a function of the molecule's chemistry rather than of its atom
# or bond count.


cdef inline uint32_t *structure_element_index(Structure structure) noexcept nogil:
    return <uint32_t *> structure.segment(SEG_ELEMENT_INDEX)


cdef inline size_t align8(size_t n) noexcept nogil:
    return (n + 7) & ~(<size_t> 7)


cdef inline size_t structure_header_len(uint16_t seg_count) noexcept nogil:
    """Byte offset of the first segment payload. The table size is data, not a constant."""
    return 24 + <size_t> 8 * seg_count


cdef Structure structure_alloc_full(uint32_t atom_count, uint32_t bond_count, bint wide,
                            uint32_t seg_mask, const uint32_t *var_len = NULL,
                            uint32_t conf_models = 0):
    """Lay out the persistent block ONCE. Nothing ever grows it again.

    `var_len`, when given, is three lengths in SEG_SGROUP_RECORD, SEG_SGROUP_INDEX, SEG_OPAQUE_BLOB
    order -- the variable-length persistent segments, whose sizes are a function of the data rather
    than of the atom and bond counts. A caller sizes them first (`blob_size_for`,
    `structure_sgroup_var_len`) and passes the answer here, which is the same two-phase shape
    `csr_build` uses. There is deliberately no way to grow one afterwards: growth would mean a
    realloc, and a persistent buffer that can be reallocated is the whole of ruling F60.

    `conf_models` IS A COUNT AND NOT A LENGTH, unlike every entry in `var_len`, and that asymmetry is
    deliberate.  SEG_CONFORMERS' length is a function of the model count AND `atom_count`, both of
    which this function already has, so passing the length would put the same arithmetic in every
    caller and give §6.1 four places to disagree about one layout.  `conformer_seg_len` below is the
    one place it is computed; a caller states how many models it has and nothing else.
    """
    cdef size_t atoms_len = align8(atom_count * sizeof(atom_t))
    cdef size_t ptr_len = align8((<size_t> atom_count + 1) * sizeof(uint32_t))
    cdef size_t edge_len = align8(<size_t> 2 * bond_count * sizeof(halfedge_t))
    cdef size_t xy_len = 0
    if seg_mask & SEG_MASK_XY:
        xy_len = align8(atom_count * sizeof(xy_t))
    cdef size_t sg_len = 0
    if seg_mask & SEG_MASK_STEREO:
        sg_len = align8(atom_count * sizeof(uint8_t))
    cdef size_t par_len = 0
    if seg_mask & SEG_MASK_PARITY:
        par_len = align8(atom_count)
    cdef size_t rec_len = 0, idx_len = 0, blob_len = 0
    if var_len is not NULL:
        rec_len = align8(var_len[0])
        idx_len = align8(var_len[1])
        blob_len = align8(var_len[2])
    cdef size_t conf_len = conformer_seg_len(conf_models, atom_count)

    # TRAILING EMPTY TABLE ENTRIES ARE NOT WRITTEN.  `seg_count` is one past the HIGHEST entry this
    # molecule actually uses, so the header is 24 + 8 * that and not a constant 128.  Measured on a
    # 13-atom aspirin with no coordinates, no stereo groups and no S-groups: 48 bytes of header
    # instead of 128, on a 704-byte molecule -- the 80 bytes saved are more than a tenth of the
    # record.
    #
    # Legal under the rule this header's own comment states: the table size is DATA, every reader
    # sizes it from the buffer's `seg_count`, and an entry a reader does not find is empty.  What is
    # NOT legal is dropping an INTERIOR empty -- ids are positional, so a molecule with a text blob
    # and no coordinates still writes an empty entry 3.  Trailing only.
    #
    # The first three are always written even when empty (a single atom has no bonds, so
    # SEG_CSR_EDGE's length is 0), because every reader dereferences those three unconditionally --
    # `structure_from_bytes` checks their sizes before it validates anything else.  So the floor is
    # SEG_CSR_EDGE + 1, and `structure_from_bytes` enforces it on ingest rather than trusting it.
    cdef uint16_t seg_count = SEG_CSR_EDGE + 1
    if par_len:
        seg_count = SEG_PARITY + 1
    elif conf_len:
        seg_count = SEG_CONFORMERS + 1
    elif blob_len:
        seg_count = SEG_OPAQUE_BLOB + 1
    elif idx_len:
        seg_count = SEG_SGROUP_INDEX + 1
    elif rec_len:
        seg_count = SEG_SGROUP_RECORD + 1
    elif sg_len:
        seg_count = SEG_STEREO_GROUPS + 1
    elif xy_len:
        seg_count = SEG_XY + 1
    cdef size_t offset = structure_header_len(seg_count)
    cdef size_t total = (offset + atoms_len + ptr_len + edge_len + xy_len + sg_len
                         + rec_len + idx_len + blob_len + conf_len + par_len)
    if total > 0xFFFFFFFF:
        raise OverflowError('structure exceeds the 4 GiB buffer limit')
    cdef Structure structure = Structure.__new__(Structure)
    structure.buffer = <char *> PyMem_Malloc(total)
    if structure.buffer is NULL:
        raise MemoryError('structure allocation failed')
    memset(structure.buffer, 0, total)
    _structure_blank(structure)
    structure.owns = True
    structure.buffer_len = total
    structure.total_len = total
    structure.header = <StructureHeader *> structure.buffer
    structure.header.magic = STRUCT_MAGIC
    structure.header.version = STRUCT_VERSION
    structure.header.flags = FLAG_WIDE_INDEX if wide else 0
    structure.header.atom_count = atom_count
    structure.header.bond_count = bond_count
    structure.header.seg_count = seg_count

    structure.header.segments[SEG_ATOMS].offset = <uint32_t> offset
    structure.header.segments[SEG_ATOMS].length = <uint32_t> atoms_len
    offset += atoms_len
    structure.header.segments[SEG_CSR_PTR].offset = <uint32_t> offset
    structure.header.segments[SEG_CSR_PTR].length = <uint32_t> ptr_len
    offset += ptr_len
    structure.header.segments[SEG_CSR_EDGE].offset = <uint32_t> offset
    structure.header.segments[SEG_CSR_EDGE].length = <uint32_t> edge_len
    offset += edge_len
    if xy_len:
        structure.header.segments[SEG_XY].offset = <uint32_t> offset
        structure.header.segments[SEG_XY].length = <uint32_t> xy_len
        offset += xy_len
    if sg_len:
        structure.header.segments[SEG_STEREO_GROUPS].offset = <uint32_t> offset
        structure.header.segments[SEG_STEREO_GROUPS].length = <uint32_t> sg_len
        offset += sg_len
    if rec_len:
        structure.header.segments[SEG_SGROUP_RECORD].offset = <uint32_t> offset
        structure.header.segments[SEG_SGROUP_RECORD].length = <uint32_t> rec_len
        offset += rec_len
    if idx_len:
        structure.header.segments[SEG_SGROUP_INDEX].offset = <uint32_t> offset
        structure.header.segments[SEG_SGROUP_INDEX].length = <uint32_t> idx_len
        offset += idx_len
    if blob_len:
        structure.header.segments[SEG_OPAQUE_BLOB].offset = <uint32_t> offset
        structure.header.segments[SEG_OPAQUE_BLOB].length = <uint32_t> blob_len
        offset += blob_len
    if conf_len:
        structure.header.segments[SEG_CONFORMERS].offset = <uint32_t> offset
        structure.header.segments[SEG_CONFORMERS].length = <uint32_t> conf_len
        offset += conf_len
    if par_len:
        structure.header.segments[SEG_PARITY].offset = <uint32_t> offset
        structure.header.segments[SEG_PARITY].length = <uint32_t> par_len
        offset += par_len
    structure.header.persistent_len = <uint32_t> offset
    _structure_resolve_persistent(structure)
    # AFTER the resolve, because the count is written THROUGH the resolved base pointer rather than
    # at a hand-computed offset -- the one arithmetic here is `conformer_seg_len`'s and nothing
    # recomputes it.  The whole buffer was memset to zero, so `reserved` and every coordinate are
    # already zero, and zero is the origin, which is what an unplaced atom reads.  `ext_index` needs
    # a write: zero is a REAL model number (a PDB `MODEL 0` is representable), so the sentinel has to
    # be stamped in rather than left as the memset's zero, or every fresh conformer would claim to
    # have been read from a file that numbered it 0.
    cdef conformer_t *rec
    cdef uint32_t model
    if conf_len:
        (<conformer_hdr_t *> structure.segment(SEG_CONFORMERS)).count = conf_models
        rec = structure_conformer_records(structure)
        for model in range(conf_models):
            rec[model].ext_index = <uint32_t> CONF_NO_INDEX
    return structure


cdef Structure structure_alloc(uint32_t atom_count, uint32_t bond_count, bint wide):
    return structure_alloc_full(atom_count, bond_count, wide, 0, NULL)


# ------------------------------------------------------------------------------------------------
# THE BLOB: an ordered list of byte strings in one segment, addressed by handle.
#
# BYTES IN, BYTES OUT, NO IMPLICIT DECODE.  MDL FIELDDATA is byte-oriented and an undecodable byte
# from another vendor's file has to survive the round trip; choosing an encoding is a policy decision
# for a layer above storage, and making it here would lose the byte permanently.  So this primitive
# knows nothing about text, and `blob_bytes` returns `bytes` rather than `str`.
#
# AN EMPTY BLOB IS 8 BYTES, NOT 0.  A zero-length segment is ABSENT -- `structure_has` keeps saying so
# and a builder that checks it would re-run forever -- so a blob that exists with no records still
# writes its two header words.  SEG_STEREO_UNIT learned this the same way.
#
# AN ABSENT BLOB COSTS NOTHING AND NEEDS NO GUARD.  `structure.segment()` answers the zero page for a
# segment that is not there, so `blob_count` reads 0 out of it and every loop over a blob is safe with
# no `structure_has` in front, exactly as the ring bitmap and the stereo unit table already rely on.
# ------------------------------------------------------------------------------------------------

cdef inline void *structure_blob(Structure structure) noexcept nogil:
    """The blob segment's base, for the one operation that treats it as bytes: copying all of it."""
    return structure.segment(SEG_OPAQUE_BLOB)


cdef inline uint32_t blob_count(Structure structure, int seg) noexcept nogil:
    """Number of records in a blob segment; 0 when the segment is absent."""
    return (<const uint32_t *> structure.segment(seg))[0]


cdef inline const uint8_t *blob_at(Structure structure, int seg, uint32_t handle,
                                   uint32_t *out_len) noexcept nogil:
    """Bytes of one record, or NULL when `handle` is past the end.

    `out_len` receives the EXACT length, which is why it is stored: the payload is 8-aligned and the
    padding would otherwise be indistinguishable from content.
    """
    cdef const uint32_t *head = <const uint32_t *> structure.segment(seg)
    if handle >= head[0]:
        out_len[0] = 0
        return NULL
    cdef const blob_rec_t *rec = (<const blob_rec_t *> (<const char *> head + 8)) + handle
    out_len[0] = rec.len
    return <const uint8_t *> head + head[1] + rec.off


cdef inline const uint32_t *blob_u32_at(Structure structure, int seg, uint32_t handle,
                                        uint32_t *out_n) noexcept nogil:
    """The same record as a `uint32` array. Sound only because payload records are 8-aligned.

    `out_n` receives the number of WHOLE words; a record whose length is not a multiple of four is a
    caller error and the tail is not reported.
    """
    # INITIALISED BECAUSE CYTHON CANNOT SEE AN OUT-PARAMETER BEING WRITTEN.  `blob_at` assigns
    # `nbytes` through the pointer on every path, but taking the address is not an assignment as far
    # as the control-flow analysis is concerned, so it emits "might be referenced before assignment"
    # -- and a Cython warning is not decoration in this file (RULES section 9.7: the same analysis
    # silently turns an undeclared local into a Python object).  The 0 is unreachable, not a default.
    cdef uint32_t nbytes = 0
    cdef const uint8_t *p = blob_at(structure, seg, handle, &nbytes)
    out_n[0] = nbytes >> 2
    return <const uint32_t *> p


cdef bytes blob_bytes(Structure structure, int seg, uint32_t handle):
    """A copy of one record, for Python. Empty bytes for a handle past the end.

    An absent record and an empty one are the same answer here ON PURPOSE, and it is the one place in
    this family where that is right: a handle past the end is a caller bug, not data, and the callers
    are all `for handle in range(...)` loops over counts the record itself supplies.  `blob_at`
    returns NULL for the one and a valid pointer for the other, so C code can still tell them apart.
    """
    cdef uint32_t n = 0                       # see `blob_u32_at` -- out-parameter, not a default
    cdef const uint8_t *p = blob_at(structure, seg, handle, &n)
    if p is NULL:
        return b''
    return <bytes> (<const char *> p)[:n]


cdef size_t blob_size_for(list items) except? 0:
    """Bytes a blob holding `items` needs, header and alignment padding included.

    SIZE FIRST, ALLOCATE, THEN FILL -- the same two-phase shape as `csr_build`, and for the same
    reason: `structure_put_blob` writes into a segment that was sized by this function, and the
    alternative (a blob that grows) would need a realloc of the persistent buffer, which is the whole
    of ruling F60.  `structure_put_blob` re-checks the size it was given rather than trusting the
    caller to have called this.

    EVERY ITEM IS TYPE-CHECKED, here by `<bytes?>` and in `structure_put_blob` by assignment to a
    `cdef bytes`.  An unchecked cast of a `str` reads its length out of the unicode header and its
    payload from the wrong offset, so a caller that forgot to encode one item stores sixteen bytes of
    CPython object header followed by a truncated string, and stores it silently.  A TypeError is the
    whole of the difference.
    """
    cdef size_t total = 8 + align8(8 * <size_t> len(items))
    cdef object item
    for item in items:
        total += align8(len(<bytes?> item))
    return total


cdef int structure_put_blob(Structure structure, int seg, list items) except -1:
    """Write `items` into an already-allocated blob segment.

    Raises when the segment is not exactly `blob_size_for(items)` long, which catches the two-phase
    contract being broken in either direction -- a caller who sized for different items, and a caller
    who forgot to size at all and got a segment that is absent.
    """
    cdef size_t need = blob_size_for(items)
    cdef uint32_t seg_len = structure_seg_len(structure, seg)
    if seg_len != need:
        raise ValueError('blob segment %d is %d bytes, these records need %d'
                         % (seg, int(seg_len), int(need)))
    cdef char *base = <char *> structure.segment(seg)
    cdef uint32_t *head = <uint32_t *> base
    cdef uint32_t count = <uint32_t> len(items)
    cdef uint32_t payload_off = <uint32_t> (8 + align8(8 * <size_t> count))
    head[0] = count
    head[1] = payload_off
    cdef blob_rec_t *rec = <blob_rec_t *> (base + 8)
    cdef uint32_t cursor = 0
    cdef uint32_t i
    cdef bytes item
    cdef Py_ssize_t n
    for i in range(count):
        item = items[i]                             # the assignment is the check; see `blob_size_for`
        n = len(item)
        rec[i].off = cursor
        rec[i].len = <uint32_t> n
        if n:
            memcpy(base + payload_off + cursor, <const char *> item, <size_t> n)
        cursor += <uint32_t> align8(<size_t> n)
    return 0


# ------------------------------------------------------------------------------------------------
# S-GROUP STORAGE.  The MDL epic owns the model; this owns the bytes.
# ------------------------------------------------------------------------------------------------

cdef inline sgroup_t *structure_sgroups(Structure structure) noexcept nogil:
    return <sgroup_t *> structure.segment(SEG_SGROUP_RECORD)


cdef inline uint32_t structure_sgroup_count(Structure structure) noexcept nogil:
    """Records in the segment. Zero for an absent segment, because the zero page has zero length."""
    return structure_seg_len(structure, SEG_SGROUP_RECORD) // <uint32_t> sizeof(sgroup_t)


cdef inline uint32_t *structure_sgroup_index(Structure structure) noexcept nogil:
    return <uint32_t *> structure.segment(SEG_SGROUP_INDEX)


cdef inline uint32_t sgroup_refs_len(sgroup_t *rec) noexcept nogil:
    """Slots in this record's reference run -- the four sub-lists are contiguous, in field order."""
    return rec.atoms_len + rec.patoms_len + rec.bonds_len + rec.cstates_len


cdef inline uint32_t sgroup_strings_len(sgroup_t *rec) noexcept nogil:
    """Handles in this record's blob run: four fixed, three counted lists, then the CSTATE tails.

    The fixed four are type, subtype, name, disp_tail -- ALWAYS PRESENT, empty bytes when unset,
    because a positional run cannot skip a slot and still be positional.  An empty `type` is a record
    whose file did not name one, which is a thing that happens.

    The CSTATE tails come LAST so that the three lists whose lengths are FIELDS stay contiguous and the
    one whose length is DERIVED is the remainder -- which means a mis-set `data_len` runs off the end of
    the record's own run and is caught, instead of quietly eating the tails.
    """
    return 4 + rec.data_len + rec.fields_len + rec.log_len + (rec.cstates_len >> 1)


cdef void structure_sgroup_var_len(Structure src, uint32_t *out) noexcept nogil:
    """The three `var_len` lengths needed to hold `src`'s S-groups unchanged.

    SIZING IS TRIVIAL, AND THAT IS THE DESIGN'S POINT.  Record count and order are preserved across
    every edit -- an emptied record STAYS PRESENT AND EMPTY (invariant 1) -- so the record segment
    never shrinks; the blob is copied byte for byte; and only the index segment can lose slots, which
    it does by compacting into a shorter run.  Sizing at `src`'s lengths is therefore always enough
    and never wrong, and the slack the compaction leaves in the index segment is zeroed and
    unreferenced.  No allocation guesswork anywhere in the carry path.
    """
    out[0] = structure_seg_len(src, SEG_SGROUP_RECORD)
    out[1] = structure_seg_len(src, SEG_SGROUP_INDEX)
    out[2] = structure_seg_len(src, SEG_OPAQUE_BLOB)


cdef int structure_carry_sgroups(Structure dst, Structure src, const int32_t *newidx,
                                 uint32_t *alias_lost) except -1:
    """Carry `src`'s S-groups into `dst`, remapping atom indices through `newidx`.

    `newidx[i]` is the new index of src atom `i`, or -1 if the atom is gone -- exactly the array
    `_apply` already computes.  Returns the number of records that LOST at least one reference, so the
    caller logs precisely that many events rather than guessing.  A record that loses everything is
    still written: invariant 1, and it is the reason this function returns a loss count instead of
    dropping records.

    ALIASES ARE COUNTED SEPARATELY, into `alias_lost`, and the split is a request from the format
    epic rather than a distinction storage cares about.  An alias is a record here like any other, but
    a writer emits it as a display label on one atom and an S-group as its own block, so one
    undifferentiated count sends a reader looking in the wrong half of the file.  `alias_lost` is
    ADDED TO and not assigned, so a caller can accumulate across calls; initialise it.

    An atom slot dies iff its atom dies.  A bond or CSTATE pair dies iff EITHER endpoint dies -- and
    the two are then treated differently, which `_sgroup_carry_pairs` explains: a bond pair is removed,
    a CSTATE pair is demoted to "unresolved" so that it keeps its vector tail.

    WHAT THIS DOES NOT MAINTAIN, and the docstring says it because the comment in the caller will be
    read later: a live pair is not a live BOND.  `delete_bond` with both endpoints alive leaves a pair
    naming two live atoms with no bond between them, and an atom-index map cannot see that.  The check
    belongs where the knowledge is -- a writer resolves each pair against the live bond table and
    drops-and-reports what does not resolve, because emitting it would produce a file whose SBL names
    a bond that is not in the bond block.  This layer keeps REFERENTIAL integrity of atom indices, not
    CHEMICAL integrity of bonds.
    """
    cdef uint32_t n = structure_sgroup_count(src)
    if n == 0:
        return 0
    cdef sgroup_t *s = structure_sgroups(src)
    cdef sgroup_t *d = structure_sgroups(dst)
    cdef uint32_t *si = structure_sgroup_index(src)
    cdef uint32_t *di = structure_sgroup_index(dst)
    cdef uint32_t r, k, run, cursor = 0, kept, pair_lost
    cdef uint32_t lost_records = 0
    cdef bint lost_here
    cdef int32_t a
    if structure_seg_len(dst, SEG_SGROUP_RECORD) < structure_seg_len(src, SEG_SGROUP_RECORD):
        raise ValueError('destination sgroup record segment is smaller than the source')
    if structure_seg_len(dst, SEG_SGROUP_INDEX) < structure_seg_len(src, SEG_SGROUP_INDEX):
        raise ValueError('destination sgroup index segment is smaller than the source')

    for r in range(n):
        memcpy(&d[r], &s[r], sizeof(sgroup_t))
        d[r].refs_off = cursor
        lost_here = False
        # atoms, then patoms: one slot each, dropped when the atom is gone.
        run = s[r].atoms_len
        kept = 0
        for k in range(run):
            a = newidx[si[s[r].refs_off + k]]
            if a >= 0:
                di[cursor + kept] = <uint32_t> a
                kept += 1
        d[r].atoms_len = kept
        lost_here |= kept != run
        cursor += kept
        run = s[r].patoms_len
        kept = 0
        for k in range(run):
            a = newidx[si[s[r].refs_off + s[r].atoms_len + k]]
            if a >= 0:
                di[cursor + kept] = <uint32_t> a
                kept += 1
        d[r].patoms_len = kept
        lost_here |= kept != run
        cursor += kept
        # bonds, then cstates: two slots each, and BOTH must live.  A CSTATE that never resolved on
        # read carries (NO_REF, NO_REF) and must survive that test, so the sentinel is checked before
        # `newidx` is indexed with it -- reading newidx[0xFFFFFFFF] is the bug this ordering avoids.
        # The `demote` flag differs between the two calls; `_sgroup_carry_pairs` says why.
        pair_lost = 0
        cursor += _sgroup_carry_pairs(si + s[r].refs_off + s[r].atoms_len + s[r].patoms_len,
                                      s[r].bonds_len, di + cursor, newidx, &d[r].bonds_len,
                                      False, &pair_lost)
        cursor += _sgroup_carry_pairs(si + s[r].refs_off + s[r].atoms_len + s[r].patoms_len
                                      + s[r].bonds_len,
                                      s[r].cstates_len, di + cursor, newidx, &d[r].cstates_len,
                                      True, &pair_lost)
        lost_here |= pair_lost != 0
        if lost_here:
            if s[r].flags & SGROUP_FLAG_ALIAS:
                alias_lost[0] += 1
            else:
                lost_records += 1
    return lost_records


cdef inline uint32_t _sgroup_carry_pairs(uint32_t *src, uint32_t run, uint32_t *dst,
                                         const int32_t *newidx, uint32_t *out_len,
                                         bint demote, uint32_t *out_lost) noexcept nogil:
    """Carry `run` slots of endpoint pairs; returns the slots written and sets `out_len` to the same.

    Two return channels for one number because the caller needs it as a cursor step AND as a field on
    the record, and computing it twice is how the two stop agreeing.  `out_lost` is incremented per
    dead pair and is a SEPARATE channel from the length for the reason `demote` exists.

    `demote` IS THE DIFFERENCE BETWEEN A BOND LIST AND A CSTATE LIST, and it is not a preference:

    * A bond pair is only a pair.  A dead one is COMPACTED OUT and the list gets shorter.
    * A CSTATE pair owns a vector tail in the blob, and the tails are a run whose length is DERIVED
      from `cstates_len` (see `sgroup_strings_len`).  Compacting a dead CSTATE out would shorten that
      run while the blob -- copied byte for byte by every carry -- keeps all its tails, and every
      surviving pair would then read the tail of the pair before it.  So a dead CSTATE is DEMOTED to
      (NO_REF, NO_REF) instead: the count does not move, the tail stays with its own pair, and the
      resulting state is one the model already has and already round-trips, because it is
      indistinguishable from a CSTATE whose bond index did not resolve when the file was read.

    That is the whole reason the loss count is reported out of band here.  For cstates the length is
    invariant across the carry, so a caller comparing lengths would see no loss and report none.
    """
    cdef uint32_t k, kept = 0
    cdef int32_t a, b
    for k in range(0, run, 2):
        if src[k] == SGROUP_NO_REF or src[k + 1] == SGROUP_NO_REF:
            # Never resolved to a bond on read; it is text and stays text, so it carries unchanged.
            dst[kept] = SGROUP_NO_REF
            dst[kept + 1] = SGROUP_NO_REF
            kept += 2
            continue
        a = newidx[src[k]]
        b = newidx[src[k + 1]]
        if a >= 0 and b >= 0:
            dst[kept] = <uint32_t> a
            dst[kept + 1] = <uint32_t> b
            kept += 2
        else:
            out_lost[0] += 1
            if demote:
                dst[kept] = SGROUP_NO_REF
                dst[kept + 1] = SGROUP_NO_REF
                kept += 2
    out_len[0] = kept
    return kept


cdef Structure structure_respan(Structure src, const uint32_t *var_len):
    """A fresh persistent block: `src`'s chemistry verbatim, S-group segments at NEW sizes.

    THE ONE OPERATION `structure_append` CANNOT DO.  A persistent segment never grows in place (F60),
    and the S-group segments are the only persistent ones whose size is a function of the DATA rather
    than of the atom and bond counts -- so replacing a molecule's S-groups means replacing its buffer.
    That is what this is, and it is why setting S-groups is a whole-set operation with no incremental
    form: an incremental one would need exactly the realloc F60 forbids.

    The atom, CSR and coordinate payloads are copied BYTE FOR BYTE, so atom indices, bond order and
    every derived answer are unchanged by construction -- which is why the caller re-derives rather
    than this function: the caller knows it is about to fill the new segments, and `rebuild_derived`
    reads none of them.

    The S-group segments arrive ZEROED and are the caller's to fill.  A zeroed record segment is not a
    valid one (`sgroup_strings_len` of an all-zero record is 4 handles it does not own), so a caller
    that allocates and does not fill leaves a molecule whose `to_bytes` its own `from_bytes` refuses --
    caught by the round trip rather than passed on.
    """
    cdef uint32_t seg_mask = 0
    if structure_has(src, SEG_XY):
        seg_mask |= SEG_MASK_XY
    if structure_has(src, SEG_STEREO_GROUPS):
        seg_mask |= SEG_MASK_STEREO
    if structure_has(src, SEG_PARITY):
        seg_mask |= SEG_MASK_PARITY
    # CONFORMERS CARRY THROUGH A RESPAN, and forgetting them here would have been a silent geometry
    # loss on a path that has nothing to do with geometry: setting S-groups on a molecule read from a
    # 3D file would have returned it flat.  Passed as the model count rather than a mask because the
    # segment's size depends on it, so unlike SEG_XY it cannot be re-derived from the atom count.
    cdef uint32_t conf_models = structure_conformer_count(src)
    cdef Structure out = structure_alloc_full(src.header.atom_count, src.header.bond_count,
                                              (src.header.flags & FLAG_WIDE_INDEX) != 0,
                                              seg_mask, var_len, conf_models)
    # `flags` carries FLAG_WIDE_INDEX, which alloc_full has already set from the argument above; every
    # other flag is a property of the chemistry and so is carried.  Copied as a whole word rather than
    # bit by bit, because a flag added later must travel with the molecule by default -- the failure
    # mode of the alternative is a new flag silently cleared by an S-group edit.
    out.header.flags = src.header.flags
    # Payload copies, by segment, using each destination's OWN length as the bound.  alloc_full sizes
    # these from the same atom and bond counts, so the lengths are equal; taking the destination's is
    # what makes that an assumption the copy cannot outrun.
    memcpy(out.atoms(), src.atoms(), structure_seg_len(out, SEG_ATOMS))
    memcpy(csr_ptr(out), csr_ptr(src), structure_seg_len(out, SEG_CSR_PTR))
    memcpy(csr_edges(out), csr_edges(src), structure_seg_len(out, SEG_CSR_EDGE))
    if seg_mask & SEG_MASK_XY:
        memcpy(structure_xy(out), structure_xy(src), structure_seg_len(out, SEG_XY))
    if seg_mask & SEG_MASK_STEREO:
        memcpy(structure_stereo_groups(out), structure_stereo_groups(src),
               structure_seg_len(out, SEG_STEREO_GROUPS))
    if conf_models:
        # THE WHOLE SEGMENT IN ONE MEMCPY, header included, because the atom count and the model count
        # are both unchanged -- so the destination's layout is the source's byte for byte and the
        # records travel with the coordinates.  The destination's own length is the bound, for the
        # reason stated above.
        memcpy(out.segment(SEG_CONFORMERS), src.segment(SEG_CONFORMERS),
               structure_seg_len(out, SEG_CONFORMERS))
    if seg_mask & SEG_MASK_PARITY:
        memcpy(structure_parities(out), structure_parities(src),
               structure_seg_len(out, SEG_PARITY))
    out.aromatic_bond_count = src.aromatic_bond_count
    return out


cdef Structure structure_with_parity(Structure src):
    """`src`'s persistent block, byte for byte, plus an EMPTY parity segment.

    For one caller: ingesting a version-3 or version-4 buffer, whose parities are in the atom flags and
    which therefore arrives with no segment to put them in.  The persistent block is laid out once and
    `structure_from_bytes` copies its buffer verbatim, so the segment cannot be added to the buffer that
    came in -- this builds the one it would have had.

    NOT `structure_respan`, which zeroes the S-group segments for a caller that is about to refill
    them.  Here every payload including the S-groups is carried, so the caller re-derives and fills
    nothing but the parities.

    Derived segments are not copied.  The caller runs `rebuild_derived` afterwards -- ingest does that
    anyway, so adoption costs one allocation and no re-derivation.
    """
    cdef uint32_t seg_mask = SEG_MASK_PARITY
    cdef uint32_t sg_var[3]
    cdef uint32_t *var_len = NULL
    if structure_has(src, SEG_XY):
        seg_mask |= SEG_MASK_XY
    if structure_has(src, SEG_STEREO_GROUPS):
        seg_mask |= SEG_MASK_STEREO
    if structure_has(src, SEG_OPAQUE_BLOB):
        structure_sgroup_var_len(src, sg_var)
        var_len = sg_var
    cdef uint32_t conf_models = structure_conformer_count(src)
    cdef Structure out = structure_alloc_full(src.header.atom_count, src.header.bond_count,
                                              (src.header.flags & FLAG_WIDE_INDEX) != 0,
                                              seg_mask, var_len, conf_models)
    out.header.flags = src.header.flags
    memcpy(out.atoms(), src.atoms(), structure_seg_len(out, SEG_ATOMS))
    memcpy(csr_ptr(out), csr_ptr(src), structure_seg_len(out, SEG_CSR_PTR))
    memcpy(csr_edges(out), csr_edges(src), structure_seg_len(out, SEG_CSR_EDGE))
    if seg_mask & SEG_MASK_XY:
        memcpy(structure_xy(out), structure_xy(src), structure_seg_len(out, SEG_XY))
    if seg_mask & SEG_MASK_STEREO:
        memcpy(structure_stereo_groups(out), structure_stereo_groups(src),
               structure_seg_len(out, SEG_STEREO_GROUPS))
    if var_len is not NULL:
        # THE TWO S-GROUP SEGMENTS AND THE BLOB, at the sizes `structure_sgroup_var_len` read off `src`
        # -- so every destination length equals its source's and the destination's is the bound.
        memcpy(out.segment(SEG_SGROUP_RECORD), src.segment(SEG_SGROUP_RECORD),
               structure_seg_len(out, SEG_SGROUP_RECORD))
        memcpy(out.segment(SEG_SGROUP_INDEX), src.segment(SEG_SGROUP_INDEX),
               structure_seg_len(out, SEG_SGROUP_INDEX))
        memcpy(structure_blob(out), structure_blob(src), structure_seg_len(out, SEG_OPAQUE_BLOB))
    if conf_models:
        memcpy(out.segment(SEG_CONFORMERS), src.segment(SEG_CONFORMERS),
               structure_seg_len(out, SEG_CONFORMERS))
    out.aromatic_bond_count = src.aromatic_bond_count
    return out


cdef Structure structure_migrate_conformers(Structure src):
    """`src`'s payload with SEG_CONFORMERS re-laid at this build's record width.

    For one caller: ingesting a version-4 or version-5 buffer that carries conformers.  Its record is
    `CONFORMER_RECORD_V5` bytes and this build's is `sizeof(conformer_t)`, so `src`'s own accessors
    cannot read its coordinates -- both bases are computed here from the old stride.  Each model's
    `ext_index` is the old record's first word, which is what that field held there too; the three
    words after it are dropped.

    Derived segments are not copied.  The caller runs `rebuild_derived` afterwards -- ingest does that
    anyway, so the migration costs one allocation and no re-derivation.
    """
    cdef uint32_t models = structure_conformer_count(src)
    cdef uint32_t atoms = src.header.atom_count
    cdef uint32_t seg_mask = 0
    cdef uint32_t sg_var[3]
    cdef uint32_t *var_len = NULL
    if structure_has(src, SEG_XY):
        seg_mask |= SEG_MASK_XY
    if structure_has(src, SEG_STEREO_GROUPS):
        seg_mask |= SEG_MASK_STEREO
    if structure_has(src, SEG_PARITY):
        seg_mask |= SEG_MASK_PARITY
    if structure_has(src, SEG_OPAQUE_BLOB):
        structure_sgroup_var_len(src, sg_var)
        var_len = sg_var
    cdef Structure out = structure_alloc_full(atoms, src.header.bond_count,
                                              (src.header.flags & FLAG_WIDE_INDEX) != 0,
                                              seg_mask, var_len, models)
    out.header.flags = src.header.flags
    memcpy(out.atoms(), src.atoms(), structure_seg_len(out, SEG_ATOMS))
    memcpy(csr_ptr(out), csr_ptr(src), structure_seg_len(out, SEG_CSR_PTR))
    memcpy(csr_edges(out), csr_edges(src), structure_seg_len(out, SEG_CSR_EDGE))
    if seg_mask & SEG_MASK_XY:
        memcpy(structure_xy(out), structure_xy(src), structure_seg_len(out, SEG_XY))
    if seg_mask & SEG_MASK_STEREO:
        memcpy(structure_stereo_groups(out), structure_stereo_groups(src),
               structure_seg_len(out, SEG_STEREO_GROUPS))
    if var_len is not NULL:
        memcpy(out.segment(SEG_SGROUP_RECORD), src.segment(SEG_SGROUP_RECORD),
               structure_seg_len(out, SEG_SGROUP_RECORD))
        memcpy(out.segment(SEG_SGROUP_INDEX), src.segment(SEG_SGROUP_INDEX),
               structure_seg_len(out, SEG_SGROUP_INDEX))
        memcpy(structure_blob(out), structure_blob(src), structure_seg_len(out, SEG_OPAQUE_BLOB))
    if seg_mask & SEG_MASK_PARITY:
        memcpy(structure_parities(out), structure_parities(src),
               structure_seg_len(out, SEG_PARITY))
    # THE OLD SEGMENT WALKED AT THE OLD STRIDE.  `src`'s length was checked as an equality against
    # `conformer_seg_len_for(models, atoms, CONFORMER_RECORD_V5)`, so both bases below are inside it.
    cdef const char *old_rec = <const char *> src.segment(SEG_CONFORMERS) + sizeof(conformer_hdr_t)
    cdef const xyz_t *old_xyz = <const xyz_t *> (old_rec + <size_t> models * CONFORMER_RECORD_V5)
    cdef conformer_t *rec = structure_conformer_records(out)
    cdef uint32_t model
    for model in range(models):
        rec[model].ext_index = (<const uint32_t *> (old_rec
                                                    + <size_t> model * CONFORMER_RECORD_V5))[0]
    if models and atoms:
        memcpy(structure_conformer_xyz(out, 0), old_xyz,
               <size_t> models * atoms * sizeof(xyz_t))
    out.aromatic_bond_count = src.aromatic_bond_count
    return out


cdef Structure structure_component_graph(Structure src, const uint32_t *slots, uint32_t m,
                                         uint32_t bonds, uint32_t *local):
    """One connected component of `src` as a structure of its own: `m` atoms, `bonds` bonds.

    For the canonical order, which decomposes a multi-component record and canonicalises each
    component alone (`_canon_order_split`).  `slots` is the component's slots ASCENDING and `local` is
    the caller's scratch of `src.header.atom_count` words, filled here with the inverse map.

    ASCENDING IS A PRECONDITION AND NOT A CONVENIENCE.  A monotone renumbering carries every
    slot-order frame onto itself, so the parity bytes stay valid without touching them -- the same
    reason `split()` preserves stereo where `substructure()` cannot.  Fed a permuted `slots` this would
    silently invert configurations.

    Carried: the atom records byte for byte (so the derived scalars, `in_ring` included, arrive already
    right -- they are component-local facts), the CSR with `to` remapped, the header flags, the
    parities.  NOT carried: coordinates, S-groups, stereo groups, conformers, and every derived
    segment.  Nothing the canonical search reads consults them; the caller runs `rebuild_derived`,
    which is what gives the component its own rings and features.
    """
    cdef uint32_t seg_mask = SEG_MASK_PARITY if structure_has(src, SEG_PARITY) else 0
    cdef Structure out = structure_alloc_full(m, bonds, (src.header.flags & FLAG_WIDE_INDEX) != 0,
                                             seg_mask)
    cdef atom_t *sa = src.atoms()
    cdef atom_t *da = out.atoms()
    cdef uint32_t *sptr = csr_ptr(src)
    cdef halfedge_t *se = csr_edges(src)
    cdef uint32_t *dptr = csr_ptr(out)
    cdef halfedge_t *de = csr_edges(out)
    cdef uint8_t *sp
    cdef uint8_t *dp
    cdef uint32_t i, j, s, fill = 0, arom = 0
    out.header.flags = src.header.flags
    for i in range(m):
        local[slots[i]] = i
    dptr[0] = 0
    for i in range(m):
        s = slots[i]
        da[i] = sa[s]
        for j in range(sptr[s], sptr[s + 1]):
            de[fill] = se[j]
            de[fill].to = local[se[j].to]
            if se[j].order == 4:
                arom += 1
            fill += 1
        dptr[i + 1] = fill
    # Recounted from the orders just copied, in the loop that copied them -- the rule every path that
    # WALKS the half-edges follows, and here the walk is unavoidable anyway.  `arom` counted half-edges.
    out.aromatic_bond_count = arom // 2
    if seg_mask:
        sp = structure_parities(src)
        dp = structure_parities(out)
        for i in range(m):
            dp[i] = sp[slots[i]]
    return out


cdef Structure structure_clone(Structure src):
    """Copy the persistent buffer verbatim, and every derived cache that has been built.

    The derived caches are copied rather than dropped or re-derived, which is not an optimisation but
    a correctness requirement: only two of the seven have an `ensure_*` guard in front of them, and
    the other five are filled by `rebuild_derived` at edit-exit and on ingest. A clone that dropped
    them would leave `structure_features`, the edge words, the element index and the ring tables
    reading the zero page, and their readers -- which do not guard -- would answer from zeros. That
    is the same class of silent wrong answer as F60, arrived at from the other direction.

    Segment layout is a function of the atom and bond counts alone, so no re-derivation is needed
    anyway: `remap()` relabels stable ids, and nothing derived depends on a stable id.
    """
    cdef Structure out = Structure.__new__(Structure)
    cdef int seg, slot
    cdef uint32_t length
    out.buffer = <char *> PyMem_Malloc(src.buffer_len)
    if out.buffer is NULL:
        raise MemoryError('structure clone allocation failed')
    memcpy(out.buffer, src.buffer, src.buffer_len)
    _structure_blank(out)
    out.owns = True
    out.buffer_len = src.buffer_len
    out.total_len = src.buffer_len
    out.header = <StructureHeader *> out.buffer
    _structure_resolve_persistent(out)
    for seg in range(SEG_PERSISTENT_COUNT, SEG_COUNT):
        length = src._seg_len[seg]
        if length == 0:
            continue
        slot = seg - SEG_PERSISTENT_COUNT
        out._derived[slot] = PyMem_Malloc(length)
        if out._derived[slot] is NULL:
            raise MemoryError('derived segment clone allocation failed')
        memcpy(out._derived[slot], src._seg_base[seg], length)
        out._seg_base[seg] = out._derived[slot]
        out._seg_len[seg] = length
        out.total_len += length
    # A verbatim copy of the bonds has the same aromatic bonds. Carried rather than recounted only
    # because the buffer is copied rather than walked here; every path that WALKS the half-edges
    # recounts instead, so no path trusts a number it did not derive from the orders in front of it.
    out.aromatic_bond_count = src.aromatic_bond_count
    return out


cdef uint32_t _blob_validate(const char *data, const StructureHeader *src, int seg) except? 0:
    """Check one blob segment against its own header and return its record count.

    Returns the count so the S-group check can bound `strings_off` without re-reading it -- the same
    "derive it once and pass it" rule the rest of this validation follows.
    """
    cdef uint32_t off = src.segments[seg].offset
    cdef uint32_t seg_len = src.segments[seg].length
    if seg_len < 8:
        raise ValueError('blob segment %d is %d bytes; the header alone is 8'
                         % (seg, int(seg_len)))
    cdef const uint32_t *head = <const uint32_t *> (data + off)
    cdef uint32_t count = head[0]
    cdef uint32_t payload_off = head[1]
    # `count` is read from the buffer, so it is bounded BEFORE it is multiplied: 8 * count overflows a
    # u32 at 2^29 records and would then name a payload offset inside the header.
    if count > (seg_len - 8) // <uint32_t> sizeof(blob_rec_t):
        raise ValueError('blob segment %d declares %d records, which do not fit in %d bytes'
                         % (seg, int(count), int(seg_len)))
    if payload_off != 8 + align8(8 * <size_t> count) or payload_off > seg_len:
        raise ValueError('blob segment %d has a payload offset of %d, expected %d'
                         % (seg, int(payload_off), int(8 + align8(8 * <size_t> count))))
    cdef const blob_rec_t *rec = <const blob_rec_t *> (data + off + 8)
    cdef uint32_t i, room = seg_len - payload_off
    for i in range(count):
        if rec[i].off > room or rec[i].len > room - rec[i].off:
            raise ValueError('blob segment %d record %d runs past the end of the segment'
                             % (seg, int(i)))
        if rec[i].off & 7:
            raise ValueError('blob segment %d record %d is misaligned' % (seg, int(i)))
    return count


cdef int _sgroup_validate(const char *data, const StructureHeader *src,
                          uint16_t seg_count, uint32_t handles) except -1:
    """Walk the S-group records of an untrusted buffer before any accessor touches them.

    THE ORDER OF THE CHECKS IS THE POINT.  Each record's reference run and blob run are bounded first,
    because every later check indexes through them; then the atom indices inside the run are bounded
    against `atom_count`, because `structure_carry_sgroups` will use them to subscript `newidx`
    without a guard.  Reversing those two would read the very slots being validated.
    """
    cdef uint32_t rec_len = src.segments[SEG_SGROUP_RECORD].length
    if rec_len % sizeof(sgroup_t):
        raise ValueError('sgroup record segment is %d bytes, not a multiple of the %d-byte record'
                         % (int(rec_len), int(sizeof(sgroup_t))))
    cdef uint32_t n = rec_len // <uint32_t> sizeof(sgroup_t)
    cdef uint32_t idx_slots = 0
    if seg_count > SEG_SGROUP_INDEX:
        idx_slots = src.segments[SEG_SGROUP_INDEX].length // <uint32_t> sizeof(uint32_t)
    cdef const sgroup_t *rec = <const sgroup_t *> (data + src.segments[SEG_SGROUP_RECORD].offset)
    cdef const uint32_t *idx = <const uint32_t *> (data + src.segments[SEG_SGROUP_INDEX].offset)
    cdef uint32_t r, k, run, strings, pairs_from
    cdef set numbered = set()

    if n and not handles:
        raise ValueError('sgroup records need a blob for their type and payload; there is none')
    for r in range(n):
        if rec[r].spare:
            raise ValueError('sgroup record %d has a non-zero spare field' % int(r))
        if rec[r].flags & ~<uint16_t> SGROUP_FLAG_DEFINED:
            raise ValueError('sgroup record %d sets an undefined flag bit (0x%x)'
                             % (int(r), int(rec[r].flags)))
        if rec[r].bonds_len & 1 or rec[r].cstates_len & 1:
            raise ValueError('sgroup record %d has an odd endpoint-pair length' % int(r))
        if rec[r].fields_len & 1:
            raise ValueError('sgroup record %d has an odd keyword/value length' % int(r))
        # The four sub-lists are contiguous from `refs_off`, so one bound covers all four -- but the
        # sum is computed in a size_t: four u32 counts can each be legal and still overflow together.
        if (<size_t> rec[r].atoms_len + rec[r].patoms_len + rec[r].bonds_len
                + rec[r].cstates_len) > idx_slots:
            raise ValueError('sgroup record %d references more index slots than the segment holds'
                             % int(r))
        run = rec[r].atoms_len + rec[r].patoms_len + rec[r].bonds_len + rec[r].cstates_len
        if <size_t> rec[r].refs_off + run > idx_slots:
            raise ValueError('sgroup record %d\'s reference run ends past the index segment'
                             % int(r))
        strings = (4 + <uint32_t> rec[r].data_len + rec[r].fields_len + rec[r].log_len
                   + (rec[r].cstates_len >> 1))
        if <size_t> rec[r].strings_off + strings > handles:
            raise ValueError('sgroup record %d\'s string run ends past the blob' % int(r))
        # SGROUP_NO_REF IS LEGAL IN THE CSTATE SLOTS AND NOWHERE ELSE, and the boundary is checked
        # rather than assumed: an atoms list holding the sentinel would subscript `newidx` at
        # 0xFFFFFFFF in the carry, four gigabytes past the array.
        pairs_from = rec[r].refs_off + rec[r].atoms_len + rec[r].patoms_len + rec[r].bonds_len
        for k in range(rec[r].refs_off, pairs_from):
            if idx[k] >= src.atom_count:
                raise ValueError('sgroup record %d references atom index %d of %d'
                                 % (int(r), int(idx[k]), int(src.atom_count)))
        for k in range(pairs_from, rec[r].refs_off + run):
            if idx[k] != <uint32_t> SGROUP_NO_REF and idx[k] >= src.atom_count:
                raise ValueError('sgroup record %d references atom index %d of %d'
                                 % (int(r), int(idx[k]), int(src.atom_count)))
        if rec[r].flags & SGROUP_FLAG_ALIAS and rec[r].atoms_len != 1:
            raise ValueError('sgroup record %d is an alias with %d atoms; an alias labels exactly one'
                             % (int(r), int(rec[r].atoms_len)))
        if rec[r].index != <uint16_t> SGROUP_NO_INDEX:
            numbered.add(rec[r].index)

    # A DANGLING PARENT IS REFUSED AT THE BOUNDARY rather than found later in a writer loop, and an
    # UNNUMBERED record cannot BE a parent because the name is the number.  Checked in a second pass
    # because a parent may be declared before the record it names -- record order is file order and
    # file order is not hierarchy order.
    for r in range(n):
        if rec[r].parent != <uint16_t> SGROUP_NO_INDEX and rec[r].parent not in numbered:
            raise ValueError('sgroup record %d names parent %d, which no record in this buffer '
                             'carries as its index' % (int(r), int(rec[r].parent)))
    return 0


cdef Structure structure_from_bytes(const char *data, size_t length):
    if length < 24:
        raise ValueError('packed molecule is too short to hold a header')

    cdef const StructureHeader *src = <const StructureHeader *> data
    if src.magic != STRUCT_MAGIC:
        raise ValueError('bad structure magic')

    # Version dispatch. A v3 buffer needs no relocation -- v4 kept the segment table at offset 24
    # and persistent ids 0-4 unchanged, so a v3 payload is already where this reader looks. What it
    # needs instead is two corrections, and both are here rather than anywhere later:
    #
    #   * `seg_count` must be supplied, because bytes 20-23 of a v3 buffer hold its `total_len`;
    #   * table entries at or above SEG_PERSISTENT_COUNT must be ZEROED, because in v3 those name
    #     DERIVED segments while in v4 ids 5-7 are the persistent S-group segments. A v3 buffer
    #     serialised after any read carries live-looking offsets there -- that is the v3 defect this
    #     release fixes -- and read by this build they would claim S-groups that do not exist, at
    #     offsets past the end of the buffer. The compatibility suite runs on frozen WARM v3 bytes to
    #     keep this branch honest, because a cold v3 buffer would pass either way.
    cdef uint16_t seg_count
    # Initialised to 0 so Cython's flow analysis accepts the hoisted v4/v5 checks that set it; the
    # v3 arm overwrites it before either branch reaches `prev_end = <uint32_t> header_len` below.
    cdef size_t header_len = 0
    # UNSIGNED BECAUSE EVERY COMPARISON IT TAKES PART IN IS AGAINST AN UNSIGNED VALUE -- a `SEG_*` enum
    # member or `seg_count`.  As a signed `int` it made each of the three `persistent_limit > SEG_*`
    # guards below a `-Wsign-compare`, which is a warning the gate counts; there is no negative limit.
    cdef uint32_t persistent_limit
    if src.version == STRUCT_VERSION:
        persistent_limit = SEG_PERSISTENT_COUNT
        seg_count = src.seg_count
    elif src.version == STRUCT_VERSION_V5:
        # THE SAME SEGMENT SET AS THIS VERSION, differing only in the conformer record's width -- which
        # is checked at the buffer's own stride below and narrowed by a rebuild after the copy.
        persistent_limit = SEG_PERSISTENT_COUNT
        seg_count = src.seg_count
    elif src.version == STRUCT_VERSION_V4:
        # A version-4 buffer states a parity in `atom_t.flags` bits 1 and 7, and this layer leaves it
        # there: it reads the buffer verbatim and normalises only the header.  Nine persistent segments,
        # so entry 9 cannot exist and the unknown-segment check below has nothing to find.
        persistent_limit = V4_PERSISTENT_COUNT
        seg_count = src.seg_count
    elif src.version == STRUCT_VERSION_V3:
        # v3 had five persistent segments, 0-4. Entries 5-11 are its derived caches, and in a buffer
        # serialised after any read they hold offsets past `persistent_len` -- so they must be
        # excluded from validation as well as zeroed, or a warm v3 buffer is rejected as having an
        # out-of-bounds segment 5. Under v4 that id belongs to SEG_SGROUP_RECORD, which is why this
        # limit is a separate number and not SEG_PERSISTENT_COUNT.
        persistent_limit = V3_PERSISTENT_COUNT
        seg_count = V3_SEG_COUNT
        header_len = V3_HEADER_LEN
        if length < header_len:
            raise ValueError('packed molecule is too short to hold a v3 header')
    else:
        raise ValueError('unsupported structure version %d' % src.version)

    # THE FIRST THREE ENTRIES MUST BE IN THE TABLE, and this is a bounds check rather than a
    # formality.  Every reader below dereferences SEG_ATOMS, SEG_CSR_PTR and SEG_CSR_EDGE
    # unconditionally -- the size checks, the CSR walk, the atom-record pass -- and the validation
    # loop skips entries past `seg_count`, so a buffer declaring `seg_count = 1` would have its
    # entries 1 and 2 read out of the atom payload and used as offsets, unvalidated.  The writer
    # never emits fewer than three (see structure_alloc_full); this is what stops a hand-made
    # buffer from claiming it did.  Not run for v3: its three constants are known sound.
    if src.version >= STRUCT_VERSION_V4:
        if seg_count == 0:
            raise ValueError('packed molecule declares an empty segment table')
        if seg_count < SEG_CSR_EDGE + 1:
            raise ValueError('packed molecule declares a %d-entry segment table; the atom and CSR '
                             'segments need %d' % (int(seg_count), SEG_CSR_EDGE + 1))
        if src.reserved0:
            raise ValueError('packed molecule has a non-zero reserved header field')
        header_len = structure_header_len(seg_count)
        if length < header_len:
            raise ValueError('packed molecule is too short for its %d-entry segment table'
                             % int(seg_count))

    if src.flags & FLAG_TOPOLOGY_DIRTY:
        raise ValueError('packed molecule has a dirty topology flag')
    if src.persistent_len != length:
        raise ValueError('packed molecule length %d does not match its header (%d)'
                         % (int(length), int(src.persistent_len)))

    # `seg` is unsigned for the same reason `persistent_limit` above is: it only ever indexes the
    # segment table and only ever compares against unsigned bounds.
    cdef uint32_t seg
    cdef uint32_t i, offset, seg_len, prev_end
    # Validate the persistent block, which is the whole of the serialised table. Derived segments are
    # not in the table at all -- they live in the `Structure` object -- so a new derived segment has no
    # obligation in this loop and nothing to join at the bottom of this function.
    #
    # For each non-empty persistent segment: check it is within the buffer, that its OFFSET and its
    # LENGTH are both 8-aligned, and that it does not overlap its predecessor. The length is checked
    # here rather than argued from the writer: `structure_respan` memcpy's the DESTINATION's align8
    # length out of the source segment, so a source whose length is not a multiple of 8 is read past
    # its end. structure_alloc_full always lays segments out in ascending index order for every
    # seg_mask and var_len combination, so tracking prev_end is sound and catches both overlap and
    # re-use of the same region.
    prev_end = <uint32_t> header_len
    for seg in range(persistent_limit):
        if seg >= seg_count:
            break
        offset = src.segments[seg].offset
        seg_len = src.segments[seg].length
        if seg_len == 0:
            # structure_alloc_full writes a non-zero offset even for empty
            # segments (e.g. SEG_ATOMS with 0 atoms, SEG_CSR_EDGE with 0 bonds),
            # so we cannot require offset == 0 here — just skip absent segments.
            continue
        if offset < header_len or offset > length or seg_len > length - offset:
            raise ValueError('segment %d is out of bounds' % seg)
        if offset & 7:
            raise ValueError('segment %d is misaligned' % seg)
        if seg_len & 7:
            raise ValueError('segment %d has a length of %d, which is not 8-aligned' % (seg, seg_len))
        if offset < prev_end:
            raise ValueError('segment %d overlaps the previous segment' % seg)
        prev_end = offset + seg_len

    # Forward compatibility, as a rule rather than a hope: a reader accepts any `seg_count`, and an
    # entry it does not know must be EMPTY. A non-empty unknown segment means the buffer carries
    # information this build cannot model, and the honest answer is to refuse it -- writing the
    # molecule back would silently drop whatever it was. Only versioned buffers are held to this; a
    # v3 buffer's entries above the persistent block are its derived caches and are ignored by design.
    if src.version >= STRUCT_VERSION_V4:
        for seg in range(persistent_limit, seg_count):
            if src.segments[seg].length:
                raise ValueError('packed molecule carries unknown segment %d, which this build '
                                 'cannot model' % seg)

    if src.segments[SEG_ATOMS].length < src.atom_count * sizeof(atom_t):
        raise ValueError('atom segment is too small for %d atoms' % src.atom_count)
    if src.segments[SEG_CSR_PTR].length < (<size_t> src.atom_count + 1) * sizeof(uint32_t):
        raise ValueError('csr pointer segment is too small for %d atoms' % src.atom_count)
    if src.segments[SEG_CSR_EDGE].length < <size_t> 2 * src.bond_count * sizeof(halfedge_t):
        raise ValueError('csr edge segment is too small for %d bonds' % src.bond_count)
    # GUARDED ON `seg_count` and not only on the length: with a truncated table the bytes at entry
    # 3's position belong to the atom payload, so an unguarded read here would test a length that is
    # really part of an atom record -- and reject a perfectly good molecule whose first atoms happen
    # to look like a short coordinate segment.  The three checks above need no guard because the
    # table is known to hold at least three entries.
    #
    # AND ZERO READS AS ABSENT, WHICH IS A BLIND SPOT ON PURPOSE.  A declared entry of length 0 and no
    # entry at all describe the same molecule -- no coordinates -- so there is nothing here to reject.
    # A writer that emits the payload and then declares 0 loses it, and NO READER CHECK CAN SEE THAT:
    # the buffer is self-consistent, and the only trace is trailing bytes nothing addresses, which legal
    # padding also produces.  Under-declaring is a writer defect and is caught in the writer.  Note the
    # consequence for mutation testing: emptying a segment does NOT exercise this line -- only a payload
    # that is short WHILE DECLARED reaches it -- so the two rejections need separate cases.
    if seg_count > SEG_XY and src.segments[SEG_XY].length and \
            src.segments[SEG_XY].length < <size_t> src.atom_count * sizeof(xy_t):
        raise ValueError('coordinate segment is too small for %d atoms' % src.atom_count)
    if seg_count > SEG_STEREO_GROUPS and src.segments[SEG_STEREO_GROUPS].length and \
            src.segments[SEG_STEREO_GROUPS].length < <size_t> src.atom_count:
        raise ValueError('stereo group segment is too small for %d atoms' % src.atom_count)
    # `persistent_limit` AS WELL AS `seg_count`, for the reason the conformer check states: entry 9
    # is SEG_PARITY only in a version-5 buffer, and in a v3 one it is a derived cache.
    if persistent_limit > SEG_PARITY and seg_count > SEG_PARITY and src.segments[SEG_PARITY].length:
        if src.segments[SEG_PARITY].length < <size_t> src.atom_count:
            raise ValueError('parity segment is too small for %d atoms' % src.atom_count)
    # THE CONFORMER SEGMENT IS CHECKED AS AN EQUALITY, WHERE THE TWO ABOVE ARE CHECKED AS A FLOOR, and
    # that is the payoff for storing `count` in the payload rather than deriving it from the length.
    # A coordinate segment can only be asked "is it big enough", because its length is the only thing
    # that states its extent; this segment states its extent TWICE -- once as a length and once as a
    # count -- so a disagreement between them is detectable and is a corrupt buffer.  Under-declaring,
    # the blind spot the note above describes for SEG_XY, is therefore not a blind spot here.
    #
    # The count is read from the payload, which is safe at this point and only at this point: the
    # bounds loop above has already established that this segment lies inside the buffer, is 8-aligned
    # and does not overlap its predecessor.  Ordering, not luck -- a check moved above that loop would
    # be reading an offset nothing has validated.
    #
    # `persistent_limit` as well as `seg_count`, for the reason the S-group checks below spell out: a
    # v3 buffer reports seg_count 13 and its entry 8 is a DERIVED CACHE, not a conformer segment.
    #
    # THE STRIDE IS THE BUFFER'S, NOT THIS BUILD'S.  The check is an equality, so reading an older
    # buffer's table at this build's width would refuse it before the narrowing downstream can run.
    cdef uint32_t conf_count
    cdef size_t conf_record = (CONFORMER_RECORD_V5 if src.version < STRUCT_VERSION
                               else sizeof(conformer_t))
    if persistent_limit > SEG_CONFORMERS and seg_count > SEG_CONFORMERS \
            and src.segments[SEG_CONFORMERS].length:
        conf_count = (<const conformer_hdr_t *> (data + src.segments[SEG_CONFORMERS].offset)).count
        if conf_count == 0:
            raise ValueError('conformer segment is present but states zero models; an absent '
                             'segment is how a molecule with no conformers is written')
        if conf_count > <uint32_t> CONF_MAX_MODELS:
            raise ValueError('conformer segment states %d models, above the %d limit'
                             % (int(conf_count), int(CONF_MAX_MODELS)))
        if src.segments[SEG_CONFORMERS].length != conformer_seg_len_for(conf_count, src.atom_count,
                                                                       conf_record):
            raise ValueError('conformer segment is %d bytes, not the %d that %d models of %d atoms '
                             'come to'
                             % (int(src.segments[SEG_CONFORMERS].length),
                                int(conformer_seg_len_for(conf_count, src.atom_count, conf_record)),
                                int(conf_count), int(src.atom_count)))
        if (<const conformer_hdr_t *> (data + src.segments[SEG_CONFORMERS].offset)).reserved:
            raise ValueError('conformer segment sets its reserved word; this build cannot model '
                             'whatever it means')

    # THE S-GROUP SEGMENTS ARE THE FIRST PERSISTENT SEGMENTS WHOSE CONTENTS POINT AT EACH OTHER, so
    # the size checks above are not enough: a record's `refs_off + len` indexes the index segment, and
    # its `strings_off + len` indexes the blob's record table, and both are read inside `nogil` by
    # every accessor.  An untrusted buffer gets those walked here or not at all.  `_sgroup_validate`
    # is a separate function only because this one is already long; it is part of this validation and
    # not an optional extra.
    #
    # THE BLOB IS VALIDATED WHETHER OR NOT THERE ARE RECORDS, because handle 0 is the molecule title
    # and a titled molecule with no S-groups is the common case.  Sizing the check to its only caller
    # is how an untrusted segment goes unwalked.
    #
    # GATED ON `persistent_limit` AND NOT ON `seg_count`.  A v3 buffer reports
    # `seg_count = V3_SEG_COUNT = 13` and its entries 5-7 are DERIVED CACHES, not S-groups -- the very
    # re-use the version dispatch above zeroes them for.  Gated on `seg_count` alone these two checks
    # read a warm v3 buffer's ring bitmap as a blob header, find a payload offset of 0 where 8 is
    # required, and refuse every warm v3 fixture in the compatibility suite.  `persistent_limit` is 5
    # for v3, 9 for v4 and SEG_PERSISTENT_COUNT for this version, so it is the one number that already
    # means "which ids are S-groups here", and asking it is what keeps this from being a second place
    # that has to remember v3's layout.
    cdef uint32_t blob_handles = 0
    if persistent_limit > SEG_OPAQUE_BLOB and seg_count > SEG_OPAQUE_BLOB \
            and src.segments[SEG_OPAQUE_BLOB].length:
        blob_handles = _blob_validate(data, src, SEG_OPAQUE_BLOB)
    if persistent_limit > SEG_SGROUP_RECORD and seg_count > SEG_SGROUP_RECORD \
            and src.segments[SEG_SGROUP_RECORD].length:
        _sgroup_validate(data, src, seg_count, blob_handles)

    # The adjacency itself must be checked, not just the size of the segment holding it.
    # rebuild_derived walks csr_ptr and halfedge_t.to inside `nogil` with boundscheck=False
    # the moment this function returns, so a ptr entry past the edge array or a `to` past the
    # atom array is an out-of-bounds access — mark_bridges even writes through one. These loops
    # are the only thing between untrusted bytes and that walk.
    cdef const uint32_t *src_ptr = <const uint32_t *> (data + src.segments[SEG_CSR_PTR].offset)
    cdef const halfedge_t *src_edges = <const halfedge_t *> (data + src.segments[SEG_CSR_EDGE].offset)
    cdef const halfedge_t *e
    cdef size_t half_edges = <size_t> 2 * src.bond_count
    cdef uint32_t k
    # Counted here rather than in a pass of its own: this loop already visits every half-edge to
    # validate its order, so the aromatic bond count is free.  It is RECOMPUTED at every point where
    # the bonds can have changed -- here, and in `csr_build` -- rather than maintained incrementally,
    # so it cannot drift out of step with the orders it counts.
    cdef uint32_t aromatic_halves = 0
    if src.atom_count:
        if src_ptr[0] != 0:
            raise ValueError('csr pointer array does not start at zero')
        for i in range(src.atom_count):
            if src_ptr[i + 1] < src_ptr[i]:
                raise ValueError('csr pointer array is not monotonic at atom %d' % i)
        if src_ptr[src.atom_count] != half_edges:
            raise ValueError('csr pointer array ends at %d, not at 2 * bond_count (%d)'
                             % (int(src_ptr[src.atom_count]), int(half_edges)))
        for k in range(half_edges):
            e = &src_edges[k]
            if e.to >= src.atom_count:
                raise ValueError('csr half-edge %d points at atom %d, outside 0..%d'
                                 % (int(k), int(e.to), int(src.atom_count - 1)))
            # ALLOWED_ORDERS is now the stored domain as well as the stated one: order 4 is a
            # representation the caller chose and this layer keeps it.
            if e.order != 1 and e.order != 2 and e.order != 3 and e.order != 4 and e.order != 8:
                raise ValueError('half-edge %d carries bond order %d, not 1, 2, 3, 4 or 8'
                                 % (int(k), int(e.order)))
            if e.wedge > 3:
                raise ValueError('half-edge %d carries wedge %d, outside 0-3'
                                 % (int(k), int(e.wedge)))
            if e.flags & ~(<uint16_t> HE_FLAG_DEFINED):
                raise ValueError('half-edge %d sets reserved flag bits 0x%04x; only 0x%04x is '
                                 'defined and every other bit must be zero'
                                 % (int(k), int(e.flags & ~(<uint16_t> HE_FLAG_DEFINED)),
                                    int(HE_FLAG_DEFINED)))
            # THE FIELD IS WIDER THAN THE DOMAIN, so the width is not the check.  Three bits hold
            # 0..7 and four descriptors are defined, and a buffer stating 5 must be refused here
            # rather than surface later as a KeyError from a table lookup in Python.
            if ((e.flags & HE_CIP_MASK) >> HE_CIP_SHIFT) > BOND_CIP_MAX:
                raise ValueError('half-edge %d carries CIP code %d, outside 0..%d'
                                 % (int(k), int((e.flags & HE_CIP_MASK) >> HE_CIP_SHIFT),
                                    int(BOND_CIP_MAX)))
            # ORDER 4 AND HE_AROMATIC ARE ONE FACT WRITTEN TWICE, so they must agree.  Requiring it
            # here is what lets every consumer downstream test whichever of the two is convenient --
            # `e.order == 4` in an order switch, `e.flags & HE_AROMATIC` in a topology switch -- and
            # get the same answer.  A buffer that sets one without the other has been written by
            # something that understood half the format, and guessing which half it meant would be
            # the silent normalisation this release exists to refuse.
            if (e.order == 4) != <bint> (e.flags & HE_AROMATIC):
                raise ValueError('half-edge %d has order %d but %s the aromatic flag; order 4 and '
                                 'HE_AROMATIC must be set together'
                                 % (int(k), int(e.order),
                                    'sets' if e.flags & HE_AROMATIC else 'does not set'))
            if e.order == 4:
                aromatic_halves += 1
    elif src.bond_count:
        raise ValueError('csr edge segment claims %d bonds in an empty molecule' % src.bond_count)

    # CSR symmetry. Every half-edge must have its twin: mark_bridges searches for it and,
    # when the search fails, stores through the sentinel ptr[child + 1] -- one past the end
    # of the edge array for the last atom. That is a heap out-of-bounds write reachable by
    # editing one `to` field of a legitimate pack(). csr_build sorts each atom's half-edges
    # by `to` (csr_build, :677-684) and add_bond rejects self-loops and duplicates, so
    # requiring strictly increasing `to` costs nothing on real input and closes both.
    #
    # The twin search is a per-atom cursor rather than a scan, which keeps the whole pass
    # linear. Correctness of the cursor: the outer loop visits i in increasing order, so for
    # a fixed j the sequence of values looked up in adj(j) is increasing, and cursor[j] never
    # needs to move backwards. Consuming the matched half-edge means matched pairs consume
    # two distinct half-edges each, so `matched == bond_count` holds if and only if every
    # half-edge is paired -- which is what catches an unmatched half-edge pointing backwards.
    cdef uint32_t *cursor = NULL
    cdef size_t matched = 0
    cdef uint32_t j
    if src.atom_count:
        cursor = <uint32_t *> PyMem_Malloc(<size_t> src.atom_count * sizeof(uint32_t))
        if cursor is NULL:
            raise MemoryError('csr symmetry scratch allocation failed')
        try:
            for i in range(src.atom_count):
                cursor[i] = src_ptr[i]
            for i in range(src.atom_count):
                for k in range(src_ptr[i], src_ptr[i + 1]):
                    e = &src_edges[k]
                    j = e.to
                    if j == i:
                        raise ValueError('atom %d carries a self-loop half-edge' % i)
                    if k > src_ptr[i] and j <= src_edges[k - 1].to:
                        raise ValueError('atom %d half-edges are not strictly increasing' % i)
                    if j > i:
                        while cursor[j] < src_ptr[j + 1] and src_edges[cursor[j]].to < i:
                            cursor[j] += 1
                        if cursor[j] >= src_ptr[j + 1] or src_edges[cursor[j]].to != i:
                            raise ValueError('half-edge %d-%d has no twin' % (int(i), int(j)))
                        # A bond has one order, so both its half-edges must agree. csr_build
                        # writes the same value into both and no mutator can separate them, but
                        # an edited buffer can: the halves then disagree and order_of(a, b) and
                        # order_of(b, a) return different bonds, with hybridization, the feature
                        # words and the signature all derived from the inconsistent graph. Only
                        # `order` is checked. `wedge` is directional by design -- it names a
                        # narrow and a wide end -- and `flags` is recomputed by mark_bridges.
                        if src_edges[cursor[j]].order != e.order:
                            raise ValueError('half-edges of bond %d-%d disagree on order: %d vs %d'
                                             % (int(i), int(j), int(e.order),
                                                int(src_edges[cursor[j]].order)))
                        cursor[j] += 1
                        matched += 1
            if matched != <size_t> src.bond_count:
                raise ValueError('csr holds %d symmetric bonds, header claims %d'
                                 % (int(matched), int(src.bond_count)))
        finally:
            PyMem_Free(cursor)

    # Stable ids are the identity model: a duplicate would alias two atoms behind one key and a
    # zero is the never-issued sentinel. Both are cheap to reject here and impossible to detect
    # later, because the index dict unpack builds would simply be one entry short.
    # 0xFFFFFFFF is rejected for a third reason: unpack sets `_next_id = high + 1`, which wraps
    # to 0 on that value, so the very next add_atom would start reissuing ids from the bottom
    # and alias an atom the packed buffer already contains.
    # The element byte is validated in the same pass: `add_atom` guarantees 1-118 on the write
    # path; element 0 is reserved for R atoms (fragment/support markers) and is accepted here
    # but not admitted by `add_atom`, so it can only arrive via `from_bytes` or byte surgery.
    cdef const atom_t *src_atoms = <const atom_t *> (data + src.segments[SEG_ATOMS].offset)
    cdef const atom_t *a
    cdef uint32_t n
    cdef set seen_ids = set()
    for i in range(src.atom_count):
        a = &src_atoms[i]
        if a.element > 118:
            raise ValueError('atom %d carries element %d; 0 is R and 1-118 are the elements'
                             % (int(i), int(a.element)))
        if a.charge < CHARGE_MIN or a.charge > CHARGE_MAX:
            raise ValueError('atom %d carries charge %d, outside %d..%d'
                             % (int(i), int(a.charge), CHARGE_MIN, CHARGE_MAX))
        if a.map_number > MAP_NUMBER_MAX:
            raise ValueError('atom %d carries map number %d, outside 0..%d'
                             % (int(i), int(a.map_number), MAP_NUMBER_MAX))
        # A MASK CHECK, NOT A DELETED CHECK: everything above the CIP nibble is rejected, which is what
        # keeps the rest of `reserved` free -- the CGR / Query hook's bits stay as protected as the
        # whole word is.
        if a.reserved & ~(<uint32_t> ATOM_RESERVED_DEFINED):
            raise ValueError('atom %d sets reserved bits 0x%08x; only 0x%08x is defined and every '
                             'other bit must be zero'
                             % (int(i), int(a.reserved & ~(<uint32_t> ATOM_RESERVED_DEFINED)),
                                int(ATOM_RESERVED_DEFINED)))
        if (a.reserved & ATOM_CIP_MASK) > ATOM_CIP_MAX:
            raise ValueError('atom %d carries CIP code %d, outside 0..%d'
                             % (int(i), int(a.reserved & ATOM_CIP_MASK), int(ATOM_CIP_MAX)))
        if a.element != 0 and (a.reserved & <uint32_t> ATOM_R_INDEX_MASK):
            raise ValueError('atom %d carries an R index and element %d; the index is only meaningful '
                             'on an R (element 0)' % (int(i), int(a.element)))
        if ((a.reserved & <uint32_t> ATOM_R_INDEX_MASK) >> ATOM_R_INDEX_SHIFT) > R_INDEX_MAX:
            raise ValueError('atom %d carries R index %d, outside 0..%d'
                             % (int(i), int((a.reserved & <uint32_t> ATOM_R_INDEX_MASK) >> ATOM_R_INDEX_SHIFT),
                                int(R_INDEX_MAX)))
        # ONLY FOR A CURRENT-VERSION BUFFER.  In version 3 and version 4 these two bits ARE the
        # atom's parity, and `MoleculeContainer.from_bytes` adopts them into SEG_PARITY; from
        # version 5 they are reserved, so a buffer that sets one carries something this build does
        # not model.
        if src.version >= STRUCT_VERSION_V5 and (a.flags & ATOM_FLAGS_RESERVED):
            raise ValueError('atom %d sets reserved flag bits 0x%02x; ATOM_FLAGS_RESERVED must be zero'
                             % (int(i), int(a.flags & ATOM_FLAGS_RESERVED)))
        n = a.n
        if n == 0:
            raise ValueError('atom %d carries stable id 0, which is never issued' % i)
        if n == 0xFFFFFFFF:
            raise ValueError('atom %d carries the reserved stable id 0xFFFFFFFF' % i)
        if n in seen_ids:
            raise ValueError('stable id %d appears twice' % n)
        seen_ids.add(n)

    cdef const uint8_t *src_par
    if persistent_limit > SEG_PARITY and seg_count > SEG_PARITY and src.segments[SEG_PARITY].length:
        src_par = <const uint8_t *> (data + src.segments[SEG_PARITY].offset)
        for i in range(src.atom_count):
            # THE FIELD IS WIDER THAN THE DOMAIN, so the width is not the check -- the same rule the
            # half-edge CIP code above is held to.  A byte stating 3 has been written by something that
            # models a fourth parity, and reading it as odd would be the guess this layer refuses.
            if src_par[i] > 2:
                raise ValueError('parity byte %d states %d; the domain is 0 none, 1 even, 2 odd'
                                 % (int(i), int(src_par[i])))

    # READ BEFORE THE HEADER IS NORMALISED.  A version-4 or version-5 buffer's conformer records are
    # wider than this build's, and the normalisation below rewrites the version byte -- so afterwards
    # nothing in the buffer says which stride its table is at.  A v3 buffer never reaches the rebuild:
    # its entry 8 is a derived cache, which the clearing loop below empties.
    cdef bint wide_conformers = src.version < STRUCT_VERSION

    cdef Structure structure = Structure.__new__(Structure)
    structure.buffer = <char *> PyMem_Malloc(length)
    if structure.buffer is NULL:
        raise MemoryError('structure allocation failed')
    memcpy(structure.buffer, data, length)
    _structure_blank(structure)
    structure.owns = True
    structure.buffer_len = length
    structure.total_len = length
    structure.header = <StructureHeader *> structure.buffer

    # Normalise a v3 buffer into this version IN PLACE. No payload moves, because v4 chose its header
    # layout so that none would have to: only the version byte, the four bytes that were `total_len`
    # and the table entries above the persistent block are rewritten.
    if structure.header.version == STRUCT_VERSION_V3:
        structure.header.version = STRUCT_VERSION
        structure.header.seg_count = SEG_TABLE_MAX
        structure.header.reserved0 = 0
    elif (structure.header.version == STRUCT_VERSION_V4
            or structure.header.version == STRUCT_VERSION_V5):
        structure.header.version = STRUCT_VERSION
    # Bounded by the buffer's OWN `seg_count`, because a v4 or v5 buffer's table may stop short of the
    # persistent block and the bytes past it are payload -- clearing to SEG_TABLE_MAX unconditionally
    # would zero atom records.  A v3 buffer always has all thirteen entries (its header is 128 bytes
    # and the line above says so), which is the case this loop exists for; a versioned buffer's
    # entries above the persistent block were already checked to be empty by the unknown-segment rule,
    # so there the loop is a no-op either way.
    #
    # FROM `persistent_limit` AND NOT FROM SEG_PERSISTENT_COUNT.  v3's derived caches occupy ids 5-11,
    # and from v4 on ids 5, 6 and 7 are the three S-GROUP segments -- so clearing from the constant
    # leaves a v3 buffer's ring bitmap, relevant-ring table and feature words installed in the
    # normalised header as a live-looking S-group record segment, index segment and blob, and the
    # first edit of that molecule memcpy's a ring bitmap into a blob and walks S-group records out of
    # a ring table.  `persistent_limit` is 5 for v3, 9 for v4 and SEG_PERSISTENT_COUNT for this
    # version, and is already the answer to "which ids does this version own"; the constant is a
    # SECOND copy of that boundary, right for one version only.
    for seg in range(persistent_limit, structure.header.seg_count):
        structure.header.segments[seg].offset = 0
        structure.header.segments[seg].length = 0

    # Counted from the half-edges validated above, halved because every bond appears twice. The
    # divisor is exact rather than rounded: the twin-symmetry check earlier in this function has
    # already established that each half-edge has a matching twin with the same order.
    structure.aromatic_bond_count = aromatic_halves // 2

    _structure_resolve_persistent(structure)
    if wide_conformers and structure.header.segments[SEG_CONFORMERS].length:
        # NOT IN PLACE, where the header normalisation above is.  Narrowing the record table by 12
        # bytes per model moves the xyz block that follows it and SEG_PARITY after that, so the
        # payload has to be re-laid.  Compacting in place would be tractable only while
        # SEG_CONFORMERS is second to last in the allocation order, and would go wrong silently the
        # day a persistent segment 10 lands.
        return structure_migrate_conformers(structure)
    return structure


def _header_size():
    return sizeof(StructureHeader)


def _atom_record_size():
    return sizeof(atom_t)


def _halfedge_size():
    return sizeof(halfedge_t)


def _xy_size():
    return sizeof(xy_t)


def _xyz_size():
    return sizeof(xyz_t)


def _conformer_record_size():
    return sizeof(conformer_t)


def _conformer_header_size():
    return sizeof(conformer_hdr_t)


def _conformer_seg_len_probe(uint32_t models, uint32_t atom_count, size_t record):
    return conformer_seg_len_for(models, atom_count, record)


def _segment_count():
    return SEG_COUNT


def _persistent_segment_count():
    return SEG_PERSISTENT_COUNT


def _segment_table_max():
    return SEG_TABLE_MAX


def _alloc_probe(uint32_t atom_count, uint32_t bond_count, bint wide, uint32_t seg_mask = 0):
    cdef Structure structure = structure_alloc_full(atom_count, bond_count, wide, seg_mask)
    cdef int seg
    # An explicit loop rather than a comprehension: a comprehension gets its own scope in Cython, so
    # the `cdef int seg` above would not reach it and the index would be an implicitly declared Python
    # object -- two warnings out of `.pxi`, which the build treats as a gate.
    cdef list segments = []
    # `seg_count` entries and not SEG_TABLE_MAX: the table stops where the molecule's last used
    # segment does, and the bytes after it are the atom payload.  A probe that read thirteen entries
    # would report atom records as segment offsets.
    for seg in range(structure.header.seg_count):
        segments.append((structure.header.segments[seg].offset,
                         structure.header.segments[seg].length))
    # `offsets` and `lengths`: SEG_TABLE_MAX entries each, zero for segments absent from the table.
    # Indexed by segment id so a probe can look up a specific id without iterating `segments`.
    cdef list offsets = []
    cdef list lengths = []
    cdef int limit
    if structure.header.seg_count < <uint16_t> SEG_TABLE_MAX:
        limit = <int> structure.header.seg_count
    else:
        limit = SEG_TABLE_MAX
    for seg in range(SEG_TABLE_MAX):
        offsets.append(<uint32_t> 0)
        lengths.append(<uint32_t> 0)
    for seg in range(limit):
        offsets[seg] = structure.header.segments[seg].offset
        lengths[seg] = structure.header.segments[seg].length
    return {'magic': structure.header.magic, 'version': structure.header.version,
            'flags': structure.header.flags, 'atom_count': structure.header.atom_count,
            'bond_count': structure.header.bond_count,
            'persistent_len': structure.header.persistent_len,
            'seg_count': structure.header.seg_count,
            'buffer_len': <uint64_t> structure.buffer_len,
            'total_len': <uint64_t> structure.total_len,
            'segments': segments,
            'offsets': offsets,
            'lengths': lengths}


def _append_isolation_probe(uint32_t atom_count, uint32_t bond_count):
    """Does building a derived segment disturb the persistent buffer? (Ruling F60.)

    F60 says a raw pointer into the arena may not be held across anything that appends a derived
    segment, because the append reallocates the arena and the old block is freed. The reason that
    rule is a hazard rather than an inconvenience is that breaking it does not crash -- the stale
    pointer reads plausible garbage out of freed memory and the molecule simply answers wrong. It has
    been broken six times by four agents.

    This probe measures the property that makes the rule necessary, rather than the rule's symptoms:

      'moved'  -- for each derived segment appended in turn, did the persistent buffer's base
                  address change? Any True means a pointer taken before the append is now dangling.
      'shared' -- for each derived segment, does its payload lie inside the SAME ALLOCATION as the
                  persistent data, i.e. within [buffer, buffer + buffer_len)? Any True means one
                  realloc can move both, which is what makes F60 possible at all.

    'shared' is the deterministic half and is the one to assert on. 'moved' depends on whether the
    allocator happened to satisfy the request in place -- on small buffers it usually does, so
    'moved' can be all False while the design is entirely unsafe, and it is reported for information
    only. 'shared' cannot be False by luck: if no derived payload shares the persistent allocation
    then there is nothing a derived append could reallocate, so 'moved' is False by construction.

    Note that "outside [buffer, buffer + persistent_len)" would be the WRONG test and would pass
    vacuously on v3, because v3 appends derived segments PAST persistent_len while still inside the
    one block. The allocation, not the persistent region, is the unit that moves.
    """
    cdef Structure structure = structure_alloc(atom_count, bond_count, False)
    cdef uintptr_t base = <uintptr_t> structure.buffer
    cdef uintptr_t payload
    cdef list moved = [], shared = [], attached = []
    cdef tuple derived = (SEG_RING_BITS, SEG_RELEVANT_RINGS, SEG_FEATURES, SEG_ELEMENT_INDEX,
                          SEG_EDGE_WORD, SEG_COMPONENT_LABEL, SEG_STEREO_UNIT)
    cdef int seg
    # a size that is generous relative to the persistent block, so that an in-buffer append has
    # every opportunity to move it; the contents are never read
    cdef size_t length = 4096
    for seg in derived:
        structure_append(structure, seg, length)
        moved.append(<uintptr_t> structure.buffer != base)
        base = <uintptr_t> structure.buffer
        payload = <uintptr_t> structure.segment(seg)
        # `structure.buffer_len` is the number of bytes allocated at `structure.buffer`. In v3 it
        # grew with every derived append, because they landed in that same block. In v4 the whole
        # point is that it stays equal to persistent_len forever.
        shared.append(base <= payload < base + <uintptr_t> structure.buffer_len)
        attached.append(bool(structure_has(structure, seg)))
    return {'moved': moved, 'shared': shared, 'attached': attached,
            'persistent_len': structure.header.persistent_len,
            'buffer_len': <uint64_t> structure.buffer_len,
            'total_len': <uint64_t> structure.total_len}


cdef dict _segment_ids():
    return {'SEG_ATOMS': <int> SEG_ATOMS, 'SEG_CSR_PTR': <int> SEG_CSR_PTR,
            'SEG_CSR_EDGE': <int> SEG_CSR_EDGE, 'SEG_XY': <int> SEG_XY,
            'SEG_STEREO_GROUPS': <int> SEG_STEREO_GROUPS,
            'SEG_SGROUP_RECORD': <int> SEG_SGROUP_RECORD,
            'SEG_SGROUP_INDEX': <int> SEG_SGROUP_INDEX,
            'SEG_OPAQUE_BLOB': <int> SEG_OPAQUE_BLOB,
            'SEG_CONFORMERS': <int> SEG_CONFORMERS,
            'SEG_PARITY': <int> SEG_PARITY,
            'CONF_NO_INDEX': <uint32_t> CONF_NO_INDEX,
            'CONF_MAX_MODELS': <int> CONF_MAX_MODELS,
            'CONF_EXT_INDEX_MAX': <uint32_t> CONF_EXT_INDEX_MAX,
            'CONFORMER_RECORD_SIZE': <int> sizeof(conformer_t),
            'CONFORMER_RECORD_V5': <int> CONFORMER_RECORD_V5,
            'SEG_PERSISTENT_COUNT': <int> SEG_PERSISTENT_COUNT,
            'SEG_TABLE_MAX': <int> SEG_TABLE_MAX,
            'SGROUP_NO_INDEX': <int> SGROUP_NO_INDEX,
            'SGROUP_INDEX_MAX': <int> SGROUP_INDEX_MAX,
            'SGROUP_LIST_MAX': <int> SGROUP_LIST_MAX,
            'SGROUP_NO_REF': <uint32_t> SGROUP_NO_REF,
            'SGROUP_FLAG_DISP': <int> SGROUP_FLAG_DISP,
            'SGROUP_FLAG_ALIAS': <int> SGROUP_FLAG_ALIAS,
            'SGROUP_RECORD_SIZE': <int> sizeof(sgroup_t),
            'SEG_RING_BITS': <int> SEG_RING_BITS,
            'SEG_RELEVANT_RINGS': <int> SEG_RELEVANT_RINGS,
            'SEG_FEATURES': <int> SEG_FEATURES,
            'SEG_ELEMENT_INDEX': <int> SEG_ELEMENT_INDEX,
            'SEG_EDGE_WORD': <int> SEG_EDGE_WORD,
            'SEG_COMPONENT_LABEL': <int> SEG_COMPONENT_LABEL,
            'SEG_STEREO_UNIT': <int> SEG_STEREO_UNIT,
            'SEG_MASK_XY': <int> SEG_MASK_XY, 'SEG_MASK_STEREO': <int> SEG_MASK_STEREO,
            'SEG_MASK_PARITY': <int> SEG_MASK_PARITY,
            'STRUCT_VERSION': <int> STRUCT_VERSION,
            'STRUCT_VERSION_V5': <int> STRUCT_VERSION_V5,
            'STRUCT_VERSION_V4': <int> STRUCT_VERSION_V4, 'STRUCT_VERSION_V3': <int> STRUCT_VERSION_V3,
            'HE_IN_RING': <int> HE_IN_RING, 'HE_AROMATIC': <int> HE_AROMATIC}


globals().update(_segment_ids())

# `DEF SEG_DERIVED_COUNT` cannot be written as SEG_COUNT - SEG_PERSISTENT_COUNT (see the note there),
# so refuse to import a build where the literal has fallen out of step with the enum.  This is the
# only place the two can be compared, because one is a Cython compile-time constant and the other an
# enumerator, and getting it wrong sizes `_derived` and `_retired` one short.
if <int> SEG_COUNT - <int> SEG_PERSISTENT_COUNT != <int> SEG_DERIVED_COUNT:
    raise ImportError('arena built with SEG_DERIVED_COUNT=%d but %d derived segments declared'
                      % (<int> SEG_DERIVED_COUNT, <int> SEG_COUNT - <int> SEG_PERSISTENT_COUNT))


cdef inline bint structure_has(Structure structure, int seg) noexcept nogil:
    return structure._seg_len[seg] != 0


cdef inline uint32_t structure_seg_len(Structure structure, int seg) noexcept nogil:
    """A segment's byte length, for persistent and derived segments alike.

    Every caller outside this file goes through here rather than reading
    `header.segments[seg].length` directly. That is not a style preference: the header table no
    longer mentions derived segments at all, so a direct read of a derived entry would report zero
    for a cache that exists. Five call sites across `_features.pxi`, `_stereo.pxi` and
    `_isomorphism.pxi` read the table by hand under v3 and are the reason this exists.
    """
    return structure._seg_len[seg]


cdef inline void _emit_half(halfedge_t *edges, uint32_t *fill, uint32_t frm, uint32_t to,
                            uint8_t order) noexcept nogil:
    """Append one half-edge to atom `frm`'s row, advancing that row's fill cursor.

    A bond is two calls with frm/to swapped, which is the whole reason this is a function:
    the two halves must agree on order and both start with a cleared wedge and flags, and
    writing them out twice by hand is how they drift apart.

    HE_AROMATIC IS SET HERE, from the order, and nowhere else on the build path.  `structure_from_bytes`
    requires order 4 and HE_AROMATIC to agree, so if the flag were set anywhere but beside the order
    it came from, a round trip would reject a molecule this build had just created.  One function, one
    line, and the two facts cannot be written apart.
    """
    cdef uint32_t pos = fill[frm]
    cdef halfedge_t *e = &edges[pos]
    fill[frm] = pos + 1
    e.to = to
    e.order = order
    e.wedge = 0
    e.flags = HE_AROMATIC if order == 4 else 0


cdef int csr_build(Structure structure, const edge_edit_t *edits,
                   uint32_t bond_count) noexcept nogil:
    cdef uint32_t n = structure.header.atom_count
    cdef uint32_t *ptr = csr_ptr(structure)
    cdef halfedge_t *edges = csr_edges(structure)
    cdef uint32_t i, b, j, k, aromatic
    cdef uint32_t *fill = <uint32_t *> malloc((n if n else 1) * sizeof(uint32_t))
    cdef halfedge_t tmp
    cdef const edge_edit_t *ed
    if fill is NULL:
        return -1

    for i in range(n + 1):
        ptr[i] = 0
    for b in range(bond_count):
        ed = edits + b
        ptr[ed.src + 1] += 1
        ptr[ed.dst + 1] += 1
    for i in range(n):
        ptr[i + 1] += ptr[i]
    for i in range(n):
        fill[i] = ptr[i]

    aromatic = 0
    for b in range(bond_count):
        ed = edits + b
        _emit_half(edges, fill, ed.src, ed.dst, ed.order)
        _emit_half(edges, fill, ed.dst, ed.src, ed.order)
        if ed.order == 4:
            aromatic += 1
    # Recomputed from the orders this call just wrote, in the loop that wrote them. Every mutation
    # rebuilds the CSR through here, so there is no path by which the count and the orders diverge.
    structure.aromatic_bond_count = aromatic

    for i in range(n):
        for j in range(ptr[i] + 1, ptr[i + 1]):
            tmp = edges[j]
            k = j
            while k > ptr[i] and edges[k - 1].to > tmp.to:
                edges[k] = edges[k - 1]
                k -= 1
            edges[k] = tmp
    free(fill)
    return 0


cdef Structure _build_csr_from_list(uint32_t atom_count, list bonds):
    cdef uint32_t m = len(bonds)
    cdef Structure structure = structure_alloc(atom_count, m, False)
    cdef edge_edit_t *edits = <edge_edit_t *> malloc((m if m else 1) * sizeof(edge_edit_t))
    cdef edge_edit_t *edit
    cdef uint32_t b
    cdef int rc
    cdef object bond
    if edits is NULL:
        # `structure` is a cdef class: unwinding drops the last reference and __dealloc__ frees
        # its buffer, so the failure path owes it nothing.
        raise MemoryError('csr edit buffer allocation failed')
    try:
        for b in range(m):
            bond = bonds[b]
            edit = &edits[b]
            edit.src = bond[0]
            edit.dst = bond[1]
            edit.order = bond[2]
        with nogil:
            rc = csr_build(structure, edits, m)
        if rc:
            raise MemoryError('csr scratch allocation failed')
    finally:
        free(edits)
    return structure


def _csr_probe(uint32_t atom_count, list bonds):
    cdef Structure structure = _build_csr_from_list(atom_count, bonds)
    cdef uint32_t *ptr = csr_ptr(structure)
    cdef halfedge_t *edges = csr_edges(structure)
    cdef halfedge_t *e
    cdef uint32_t i, k
    cdef dict neighbors = {}
    cdef list canonical = []
    cdef list row
    for i in range(atom_count):
        row = []
        for k in range(ptr[i], ptr[i + 1]):
            e = &edges[k]
            row.append((e.to, e.order))
            if e.to > i:
                canonical.append((i, e.to, e.order))
        neighbors[i] = row
    cdef list ptr_out = []
    for i in range(atom_count + 1):
        ptr_out.append(ptr[i])
    return {'ptr': ptr_out, 'neighbors': neighbors, 'canonical': sorted(canonical)}


def _csr_find_probe(uint32_t atom_count, list bonds):
    cdef Structure structure = _build_csr_from_list(atom_count, bonds)
    cdef halfedge_t *found
    cdef uint32_t i, j
    cdef dict out = {}
    for i in range(atom_count):
        for j in range(atom_count):
            if i == j:
                continue
            found = csr_find(structure, i, j)
            out[(i, j)] = None if found is NULL else found.order
    return out


def _zero_page_probe(uint32_t atom_count, uint32_t bond_count):
    cdef Structure structure = structure_alloc(atom_count, bond_count, False)
    cdef int32_t *xy = <int32_t *> structure.segment(SEG_XY)
    cdef int i
    cdef list reads = []
    for i in range(2 * atom_count):
        reads.append(xy[i])
    return {'has_xy': structure_has(structure, SEG_XY), 'xy_reads': reads,
            'has_atoms': structure_has(structure, SEG_ATOMS)}


def _atom_field_probe():
    cdef Structure structure = structure_alloc(1, 0, False)
    cdef atom_t *a = structure.atoms()
    a.element = 6
    a.charge = -1
    a.isotope = 13
    a.map_number = 4095
    a.n = 7
    at_set_h(a, 3, 2)
    at_set_radical(a, True)
    at_set_h_pinned(a, True)
    at_set_in_ring(a, True)
    at_set_hybridization(a, 5)
    at_set_ring_counts(a, 2, 1)
    return {'element': a.element, 'charge': a.charge, 'isotope': a.isotope,
            'map_number': a.map_number, 'n': a.n,
            'implicit_h': at_implicit_h(a), 'explicit_h': at_explicit_h(a),
            'radical': at_radical(a), 'h_pinned': at_h_pinned(a), 'flags': a.flags,
            'in_ring': at_in_ring(a), 'hybridization': at_hybridization(a),
            'rings_count': at_ring_count(a),
            'aromatic_ring_count': at_aromatic_ring_count(a)}


def _atom_set_hybridization_probe(int z):
    if z < 1 or z > 6:
        raise ValueError('hybridization must be 1-6')
    cdef Structure structure = structure_alloc(1, 0, False)
    at_set_hybridization(structure.atoms(), <uint8_t> z)
    return at_hybridization(structure.atoms())


def _atom_h_nibble_probe(uint8_t implicit, uint8_t explicit):
    cdef Structure structure = structure_alloc(1, 0, False)
    cdef atom_t *a = structure.atoms()
    at_set_h(a, implicit, explicit)
    return {'implicit_h': at_implicit_h(a), 'explicit_h': at_explicit_h(a),
            'raw': a.hydrogens}


cdef int structure_append(Structure structure, int seg, size_t length) except -1:
    """Attach a DERIVED segment, in its own allocation.

    RULING F60 DOES NOT REACH THIS FUNCTION. It does not touch `structure.buffer`, it does not touch
    the header, and it cannot move anything: the new cache is a separate block, so every pointer any
    caller holds into the arena -- or into any other derived segment -- is still valid when this
    returns. `rebuild_derived`, `ensure_stereo_units` and `ensure_component_labels` may therefore be
    called in any order relative to any pointer fetch, and no hoist or warning is needed to protect
    that ordering.

    Persistent segments are refused. They are laid out once by `structure_alloc_full` and a growing
    persistent buffer is exactly the hazard this layout excludes, so making it unreachable here is
    what keeps it excluded -- a change cannot reintroduce it by calling this with a persistent id.
    """
    cdef size_t seg_len = align8(length)
    cdef int slot = seg - SEG_PERSISTENT_COUNT
    if not structure.owns:
        raise RuntimeError('cannot grow a borrowed structure')
    if seg < SEG_PERSISTENT_COUNT or seg >= SEG_COUNT:
        raise RuntimeError('segment %d is persistent; only derived segments are appended'
                           % int(seg))
    if structure._seg_len[seg]:
        raise RuntimeError('segment %d is already attached' % int(seg))
    if seg_len > 0xFFFFFFFF or structure.total_len + seg_len > 0xFFFFFFFF:
        raise OverflowError('structure exceeds the 4 GiB addressable limit')
    cdef void *block = PyMem_Malloc(seg_len if seg_len else 1)
    if block is NULL:
        raise MemoryError('structure segment allocation failed')
    memset(block, 0, seg_len if seg_len else 1)
    structure._derived[slot] = block
    structure._seg_base[seg] = block
    structure._seg_len[seg] = <uint32_t> seg_len
    structure.total_len += seg_len
    return 0


cdef int structure_retire(Structure structure, int seg) except -1:
    """Detach a derived segment so it can be rebuilt, keeping the old block readable.

    The block is RETIRED, not freed: it is held one deep and released at `__dealloc__`. So a pointer
    taken before this call still reads the previous, correct contents rather than freed memory --
    which is the same guarantee v3 gave by accident (it stranded the old table inside the buffer and
    documented that the bytes were not reclaimed) except that here the block is tracked instead of
    leaked, and the guarantee is the reason for the design rather than a side effect of it.

    One deep is enough because it makes a stale pointer read stale-but-valid data across a single
    replacement, and a caller that rebuilds the same cache twice while holding a pointer from before
    the first rebuild has a bug this layer cannot paper over.
    """
    cdef int slot = seg - SEG_PERSISTENT_COUNT
    if seg < SEG_PERSISTENT_COUNT or seg >= SEG_COUNT:
        raise RuntimeError('segment %d is persistent and cannot be retired' % int(seg))
    if structure._seg_len[seg] == 0:
        return 0
    if structure._retired[slot] is not NULL:
        PyMem_Free(structure._retired[slot])
    structure._retired[slot] = structure._derived[slot]
    structure._derived[slot] = NULL
    structure.total_len -= structure._seg_len[seg]
    structure._seg_base[seg] = <void *> _zero_page
    structure._seg_len[seg] = 0
    return 0


cdef void derive_scalars(Structure structure) noexcept nogil:
    """
    Fill the per-atom scalars that follow from the graph alone: explicit hydrogen count,
    heteroatom count and hybridization. Pure CSR derivations with no element table, which
    is why they live beside the layout they read rather than in a module of their own.

    Every field written here is written authoritatively on every atom. `_apply` memcpys
    each `atom_t` field forward from the old arena, so a field this pass only accumulates
    into would keep a previous edit's value.
    """
    cdef atom_t *atoms = structure.atoms()
    cdef uint32_t *ptr = csr_ptr(structure)
    cdef halfedge_t *edges = csr_edges(structure)
    cdef atom_t *a
    cdef halfedge_t *e
    cdef uint32_t n = structure.header.atom_count
    cdef uint32_t i, k
    cdef int doubles, triples, aromatics, hetero, explicit, order, z
    cdef uint8_t nb_element

    for i in range(n):
        a = &atoms[i]
        doubles = 0
        triples = 0
        aromatics = 0
        hetero = 0
        explicit = 0
        for k in range(ptr[i], ptr[i + 1]):
            e = &edges[k]
            order = e.order
            if order == 2:
                doubles += 1
            elif order == 3:
                triples += 1
            elif order == 4:
                aromatics += 1
            nb_element = atoms[e.to].element
            if nb_element == 1:
                explicit += 1
            elif element_is_heteroatom(nb_element):
                hetero += 1

        a.heteroatoms = <uint8_t> (hetero if hetero < 255 else 255)
        # Saturated at the EXPLICIT count's bound, which is the nibble's width because that nibble
        # has no sentinel to make room for.  Spelled H_EXPLICIT_MAX rather than H_NIBBLE_MAX so the
        # day the explicit nibble acquires a reserved value this line moves with it.
        if explicit > H_EXPLICIT_MAX:
            explicit = H_EXPLICIT_MAX

        # AROMATIC FIRST, and it is not a tie-break but the definition: 4 is the domain table's
        # reserved meaning for aromatic hybridization, and a stored order-4 bond is the only thing
        # that produces it.  An atom on a stored aromatic ring has no double and no triple bond, so
        # without this branch it falls through to `z = 1` and benzene answers sp3 -- matching
        # `[C;z1]` and refusing `[C;a]`, a plausible wrong answer rather than a visible refusal.
        #
        # An atom with an aromatic bond AND an exocyclic double bond (a quinone-like carbon written
        # part-aromatic, or a pyridine N-oxide) reports 4 as well.  Aromaticity is the stronger
        # statement about the atom's environment and the one a query asks about.
        if aromatics:
            z = 4
        elif triples == 0 and doubles == 0:
            z = 1
        elif triples == 0 and doubles == 1:
            z = 2
        elif triples == 1 and doubles == 0:
            z = 3
        elif triples == 0 and doubles == 2:
            z = 5
        else:
            z = 6
        at_set_hybridization(a, <uint8_t> z)

        # implicit_h is the caller's data -- read it back and write it unchanged, H_UNKNOWN included.
        # This is the one place the sentinel passes through a WRITE, and it must pass through: the
        # derivation owns the explicit nibble and nothing else, so re-deriving scalars on a record
        # with an unstated hydrogen count must not quietly resolve it to zero.
        at_set_h(a, at_implicit_h(a), <uint8_t> explicit)


cdef Py_ssize_t label_components(Structure structure, uint32_t *label) noexcept nogil:
    """Write each atom's component index into `label`, and return the component count.

    Components are numbered in order of their lowest-indexed member, so the numbering is a
    function of the arena alone. An isolated atom is a component of its own -- salt splitting
    and the circuit rank both need it counted. Returns -1 if the DFS stack cannot be allocated.

    Nothing caches this. The traversal is one pass over the CSR, cheaper than the invalidation
    bookkeeping a cached copy would need to stay honest across edits.
    """
    cdef uint32_t n = structure.header.atom_count
    if n == 0:
        return 0
    cdef uint32_t *ptr = csr_ptr(structure)
    cdef halfedge_t *edges = csr_edges(structure)
    cdef uint32_t NONE = 0xffffffff
    cdef uint32_t i, k, v, to, top
    cdef Py_ssize_t comps = 0
    for i in range(n):
        label[i] = NONE

    # each atom is pushed at most once, so n slots is the exact worst case
    cdef uint32_t *stack = <uint32_t *> malloc(<size_t> n * sizeof(uint32_t))
    if stack is NULL:
        return -1
    for i in range(n):
        if label[i] != NONE:
            continue
        label[i] = <uint32_t> comps
        stack[0] = i
        top = 1
        while top:
            top -= 1
            v = stack[top]
            for k in range(ptr[v], ptr[v + 1]):
                to = edges[k].to
                if label[to] == NONE:
                    label[to] = <uint32_t> comps
                    stack[top] = to
                    top += 1
        comps += 1
    free(stack)
    return comps


cdef int ensure_component_labels(Structure structure) except -1:
    """Fill SEG_COMPONENT_LABEL if it is not there yet; idempotent.

    Only queries that use component grouping need this, and those are rare, so the segment is
    filled on demand rather than by rebuild_derived.  Not nogil: it reallocates the arena through
    structure_append -> PyMem_Realloc, so every caller must invoke it *before* taking any pointer
    into the arena buffer.  In particular, the matcher must call this before caching
    structure_edge_words, the CSR edge array, or the atom array into a struct -- those pointers
    all become dangling after the reallocation.  A search that returns to Python between
    solutions cannot rely on that ordering, because the caller's loop body may append; those two
    generators call matcher_reseat on every resume instead.
    """
    if structure_has(structure, SEG_COMPONENT_LABEL):
        return 0
    cdef uint32_t n = structure.header.atom_count
    structure_append(structure, SEG_COMPONENT_LABEL,
                     sizeof(uint32_t) * (<size_t> n if n else 1))
    cdef uint32_t *label = structure_component_labels(structure)
    cdef Py_ssize_t comps
    with nogil:
        comps = label_components(structure, label)
    if comps < 0:
        raise MemoryError('component labelling scratch allocation failed')
    return 0
