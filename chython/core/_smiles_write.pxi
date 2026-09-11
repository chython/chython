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
# The SMILES writer: an arena to a SMILES string.
#
# Design: docs/superpowers/specs/2026-09-02-smiles-write-design.md.  Every module-level symbol
# here is prefixed `smw_`, agreed with the SMILES READER on 2026-09-02: the reader owns `smi_*` in
# `_smiles_read.pxi` and `arom_*` in `_kekule.pxi`, and this file touches neither.
#
# FOUR THINGS TO KNOW BEFORE EDITING.
#
# 1. THE WRITER PERCEIVES NOTHING AND MUTATES NOTHING.  It spells exactly what the molecule holds:
#    a stored order-4 bond becomes lowercase aromatic SMILES, a stored Kekule bond becomes an
#    explicit `=`.  It must never call `thiele()` or `kekule()` on the way -- a caller who wants the
#    aromatic spelling of a Kekule molecule calls `thiele()` first, deliberately.  A silent
#    representation change is the sin, and the two named operations are the only places a
#    representation may change.
#
#    ORDER 4 IS STORED as of arena v4, and `e.flags & HE_AROMATIC` is set iff `e.order == 4` by the
#    arena's own construction, so aromaticity here is a FACT TO READ and never something to infer.
#    Deciding it from an option standing in for a perception result writes an aromatic benzene as
#    `[CH]1[CH][CH][CH][CH][CH]1` -- cyclohexane, a silent representation change on the way OUT,
#    where the caller has no way to notice it.  `test_smiles_write_aromatic.py` pins it.
#
# 2. NO STORED PARITY BYTE REACHES THE OUTPUT (ruling F26).  A parity byte's reference frame is
#    the order the atoms were CREATED in, so a writer that emits `@`/`@@` from it is right for
#    some input orderings and silently wrong for others.  Every sign this file emits is computed
#    by `smw_sign_of`, which builds the direction list IN THE ORDER IT IS ABOUT TO WRITE
#    and translates the parity into that order.  If you add an emission path, route it there.
#
# 3. CANONICAL OUTPUT IS A FUNCTION OF THE CANONICAL POSITIONS AND OF NOTHING ELSE.  Component
#    order, start atom, neighbour order and closure numbering are all decided from `pos`, so
#    "same molecule, any creation order, same string" reduces to a property of `canonical_order`.
#    Any new decision in this file must be made from `pos` too -- a tie broken on a slot index
#    reintroduces the creation order through the side door, which is the shape of the ~8%
#    canonical-stereo oscillation this design exists to close.
#
# 4. NO TOKEN IS SUPPRESSED ON A PREDICTION, AND AN UNSPELLABLE FACT IS REPORTED.  Two halves of one
#    rule.  Note 1's silent loss is exactly a prediction: drop the `:` in the belief that the ATOM
#    will come out lowercase, let the mechanism that was to lowercase it not run, and the fact leaves
#    the string entirely.  So a decision to say less is made from WHAT WAS ACTUALLY EMITTED and never
#    from what another mechanism is expected to emit -- `smw_bond` asks `smw_both_lowercase`, not
#    `o.aromatic_bond`.  A redundant token is noise; a suppressed one is silent loss, and the
#    asymmetry is total.
#
#    When the notation cannot carry the fact AT ALL, the atom or unit goes into a REPORT rather than
#    quietly out of the string: `lost` for a configuration (an atropisomer, an axis under `!b`, a
#    contradictory cis/trans set, a frame whose hydrogen position is unknown) and `unknown_h` for the
#    arena's H_UNKNOWN sentinel, which a bracket cannot spell because an absent H term there means
#    zero.  Three instances of one principle, not three special cases.


# The writer's own domains.  H_IMPLICIT_MAX (14, one short of the nibble because H_UNKNOWN takes 15),
# H_EXPLICIT_MAX, CHARGE_MIN/MAX, ISOTOPE_MAX and MAP_NUMBER_MAX are declared in
# `_molecule_arena.pxi` and are not restated here (RULES.md §6).  H_NIBBLE_MAX is a LAYOUT width and
# not a bound on any count a caller may state -- a bound that admits the sentinel destroys it.
DEF SMW_MAX_CLOSURE = 99      # ring-closure numbers 1..99; `%NN` from 10 up
DEF SMW_BUF_MIN = 64          # initial output buffer, doubled on demand
# THE SMALLEST ATTACHMENT ID a detached fragment may use, and it is 10 BECAUSE OF A MEASUREMENT rather
# than for elbow room.  A detached attachment is always spelled `%NN` so that the token has a fixed
# width and can never be a bare digit, and the two-digit form with a LEADING ZERO -- `%05`, which
# OpenSMILES' `'%' DIGIT DIGIT` grammar plainly allows -- is refused by RDKit 2026.03.4 and by chython
# 2's own parser ("number starts with 0"), while Indigo and OpenBabel accept it.  So the fixed-width
# spelling and universal readability agree only above 9.  See `smw_apply_cuts` and spec §13.3.
DEF SMW_MIN_ATTACH = 10


# 0xFFFFFFFF as a typed global rather than a DEF, for the same reason as CANON_NO_SLOT: it is not
# representable in the `int` a C enumerator has to fit, and it is compared inside `nogil`.
cdef uint32_t SMW_NONE = 0xFFFFFFFF


# The element symbols again, as C string literals indexed by atomic number.  `SYMBOLS` in
# `_elements.pxi` is a Python tuple, so reading it in `smw_atom` would mean an `str` index, a
# `.encode()` and a temporary `bytes` per atom -- refcount traffic in the one function called once
# per atom, and a function that then cannot be `noexcept` honestly.  The duplication is checked:
# test_symbol_table_matches_elements compares this table to `SYMBOLS` entry by entry, so the two
# cannot drift.  The atomic number indexes directly.
cdef extern from *:
    """
    static const char SMW_SYMBOL[119][3] = {
    "R",   /* the fragment marker, element 0; its index, when nonzero, follows the symbol */
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S",
    "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga",
    "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd",
    "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm",
    "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os",
    "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa",
    "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg",
    "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"};
    """
    const char SMW_SYMBOL[119][3]


def smw_symbol_table():
    """Expose SMW_SYMBOL to the test suite as a tuple of 119 strings; index 0 is 'R', the marker."""
    cdef uint32_t i
    cdef list out = []
    for i in range(119):
        out.append((<bytes> SMW_SYMBOL[i]).decode('ascii'))
    return tuple(out)


cdef enum:
    # `he_kind` values.  A half-edge is classified exactly once, from whichever side reaches it
    # first, and the other side reads the mark and skips.
    SMW_HE_UNSEEN = 0
    SMW_HE_CHILD = 1          # tree edge, written from this side (this side is the parent)
    SMW_HE_PARENT = 2         # tree edge, this side is the child
    SMW_HE_CLOSURE = 3        # non-tree edge: a ring closure, both halves marked
    SMW_HE_ATTACH = 4         # a detached fragment's CUT bond: written as `%NN` on the retained side


cdef struct smw_opts_t:
    bint canonical          # canonical atom order; False = stored slot order
    bint random_order       # a fresh random atom order, refused together with `canonical = False`
    bint stereo             # emit @/@@ and /\
    bint aromatic_bond      # ':' bonds and UPPERCASE atoms; default is lowercase and no token
    bint mapping            # emit :N from the atom's map number
    bint hydrogens          # force brackets and an explicit H count on every atom
    bint bonds              # emit bond tokens at all
    bint charges            # emit charges
    bint cxsmiles           # append the CXSMILES tail
    bint asymmetric_closure  # write the bond token on the opening side of a closure only


# THE CUTS of a detached fragment (spec §13).  Fixed arrays and no allocation: an attachment id is
# 10..99, so 90 cuts is the ceiling the id space imposes and a struct that cannot fail to allocate has
# no failure path to get wrong.  A thousand bytes of stack, once per write.
#
# `keep`/`drop` are SLOTS, resolved from the caller's stable ids before anything else happens.  ORDERED,
# because the writer cannot infer which side of `C-C` the caller wants to keep.
cdef struct smw_cuts_t:
    uint32_t ncuts
    uint32_t keep[SMW_MAX_CLOSURE + 1]
    uint32_t drop[SMW_MAX_CLOSURE + 1]
    uint8_t ids[SMW_MAX_CLOSURE + 1]
    # Indexed BY ID and not by cut, because it holds more than this fragment's own attachments: a
    # caller joining three fragments reserves every id in the whole set on each of them, so that an
    # internal closure in one cannot collide with an attachment in another.
    uint8_t reserved[SMW_MAX_CLOSURE + 1]   # [id] 1 = withheld from the internal closure allocator


# THE STICKY ENDS (§14).  `sticky_smiles`'s contract is a string a caller GLUES rather
# than a fragment that re-joins by ring bonds -- so the two ends are the FIRST and LAST tokens of the
# string and one of them may be missing its atom.  Nothing here is a cut: both atoms stay in the
# traversal, and `remove_*` suppresses TOKENS.
#
# `left`/`right` are SLOTS, SMW_NONE when the caller named neither.  A sticky write is not canonical --
# the order is forced at both ends -- and is not a cache key; see `normalize_smiles_spec`.
cdef struct smw_sticky_t:
    uint32_t left           # written first, or SMW_NONE
    uint32_t right          # written last, or SMW_NONE
    bint remove_left        # suppress `left`'s atom token
    bint remove_right       # suppress `right`'s atom token
    bint keep_bond_left     # ... but keep its bond token, forced non-empty, for the caller to glue to
    bint keep_bond_right


cdef struct smw_buf_t:
    char *data
    size_t length
    size_t cap
    bint oom                # sticky: every put becomes a no-op and the caller raises once


cdef struct smw_scratch_t:
    uint32_t n
    uint32_t nhe            # half-edge count == csr_ptr[n] == 2 * bond count
    uint32_t nseq           # atoms placed in `seq` so far; ends equal to n
    uint32_t *pos           # per slot: canonical position, or the slot itself in stored mode
    uint32_t *bypos         # per position: the slot holding it
    uint32_t *seq           # emission order, as slots
    uint32_t *out_idx       # per slot: its index in `seq` (its position in the output)
    uint32_t *parent        # per slot: the slot it was reached from, or SMW_NONE
    uint32_t *nbr_off       # n + 1 offsets into `nbr` / `wnbr`; nbr_off[i+1]-nbr_off[i] = degree
    uint32_t *nbr           # half-edge indices, per atom ascending by `pos` of the far atom
    uint32_t *wnbr          # the same indices in WRITTEN order: parent, ring bonds, children
    uint32_t *stack         # DFS stack of slots
    uint32_t *cursor        # per slot: its frame's cursor into `nbr` / `wnbr`
    uint8_t *he_kind        # per half-edge: SMW_HE_*
    uint8_t *he_close       # per half-edge: ring-closure number, 0 when none
    uint8_t *he_dir         # per half-edge: 0 none, 1 '/', 2 '\' reading from this side
    uint8_t *visited        # per slot
    uint8_t *paren          # per slot: its subtree was opened with '(' and must close with ')'
    uint8_t *lost           # per slot: it anchors a unit whose configuration could not be written
    uint8_t *dropped        # per slot: cut away from this fragment; all zero for a whole-molecule write
    uint8_t *frame_cut      # per slot: a sticky removal took a neighbour away, so no sign may be written
    uint32_t *path_rank     # per slot: its index on the sticky left..right path, or SMW_NONE
    uint32_t *unit_at       # per slot: the anchor of the cis/trans unit it terminates, or SMW_NONE
    uint32_t *partner_at    # per slot: the OTHER terminal of that unit
    uint32_t *dq            # direction queue: (from slot, half-edge) pairs, 2 * nhe words
    uint32_t nlost          # how many units `smw_directions` could not write
    void *block             # the single allocation every pointer above points into


# ------------------------------------------------------------------------------------------------
# THE VALENCE MODEL.  One table, and it is a CONTRACT WITH THE READER, not a private helper.
#
# PREFIXED `smv_`, NOT `smw_`, DELIBERATELY.  Everything else in this fragment is `smw_` because the
# writer epic owns it; the reader owns `smi_*`.  These four functions are owned by neither: the
# reader calls `smv_default_h` to decide what a bare `N` means, and a reviewer who saw the reader
# reaching into `smw_*` would file a layering violation and would be reading the name correctly.
# The name would be claiming single ownership of something two directions depend on.
#
# It is also NOT the same table as `_valence.pxi`'s `val_*`, which is the reader epic's chemical
# valence model -- 1036 rules, consulted by MDL read, InChI and `check_valence`.  That one
# answers "is this structure chemically possible", takes the atom's ENVIRONMENT, and must never be
# called from here: strict output means spec-conformant SYNTAX, never validated content, or the
# library stops being able to show a user their own broken structure.  This one answers "will a
# reader infer the count I hold", takes a bond-order SUM, and cannot see an environment at all --
# a SMILES atom's syntax may not depend on what it is bonded to.  Different arity, so they are not
# one table with two entry points, and `test_the_two_valence_models_answer_differently` in this
# module's suite (mirrored by one of the same name in `test_valence.py`) fails if they are merged.
#
# An atom is written without brackets only when the hydrogen count the reader will INFER equals the
# count the arena holds.  THE CHECK IS THE RULE, so the cases a fixed list of bracket reasons has to
# name one at a time -- an aromatic B/N/P bearing hydrogen, an elemental B/C/P/S, a hypervalent P --
# fall out of it, and a carbon with two implicit hydrogens and no bonds cannot write a bare `C` that
# reads back as methane.
#
# TWO valence sets per element, and both must agree before a count may be omitted, because the
# models real readers use are not the same model.  Measured against RDKit 2026.03.4 by parsing
# `X(F)(F)...` and `X(=O)...` with sanitize off and reading GetTotalNumHs:
#
#            RDKit             OpenSMILES §3.1.5
#   B        3                 3
#   C        4                 4
#   N        3                 3, 5
#   O        2                 2
#   P        3, 5              3, 5
#   S        2, 4, 6           2, 4, 6
#   F Cl Br  1                 1
#   I        1, 3, 5, 7        1
#
# So a neutral four-bonded nitrogen bearing one hydrogen is `NH` to OpenSMILES and `N` to RDKit,
# and a two-bonded iodine bearing one hydrogen is the mirror image.  Writing either without
# brackets loses the hydrogen for half the world's readers.  Taking the NARROW set and the WIDE set
# and demanding they agree turns both cases into `[NH]` and `[IH]`, which every reader gets right.
# The disagreement is confined to those two elements; for the other eight the two sets are equal
# and the check costs nothing.
cdef inline void smv_valences(uint32_t element, uint32_t *narrow,
                              uint32_t *wide) noexcept nogil:
    """Bitmask of normal valences per element: bit v set means v is a normal valence.

    Zero for every element outside the SMILES organic subset, which is how the bracket predicate
    asks "is this element writable bare at all" -- one question, one table.
    """
    narrow[0] = 0
    wide[0] = 0
    if element == 5:                                    # B
        narrow[0] = 1u << 3
        wide[0] = 1u << 3
    elif element == 6:                                  # C
        narrow[0] = 1u << 4
        wide[0] = 1u << 4
    elif element == 7:                                  # N
        narrow[0] = 1u << 3
        wide[0] = (1u << 3) | (1u << 5)
    elif element == 8:                                  # O
        narrow[0] = 1u << 2
        wide[0] = 1u << 2
    elif element == 15:                                 # P
        narrow[0] = (1u << 3) | (1u << 5)
        wide[0] = narrow[0]
    elif element == 16:                                 # S
        narrow[0] = (1u << 2) | (1u << 4) | (1u << 6)
        wide[0] = narrow[0]
    elif element == 9 or element == 17 or element == 35:     # F, Cl, Br
        narrow[0] = 1u << 1
        wide[0] = 1u << 1
    elif element == 53:                                 # I
        narrow[0] = 1u << 1
        wide[0] = (1u << 1) | (1u << 3) | (1u << 5) | (1u << 7)


cdef inline uint32_t smv_h_from_mask(uint32_t mask, uint32_t order_sum) noexcept nogil:
    """Implicit hydrogens a reader with these normal valences infers for this bond-order sum.

    The smallest normal valence at or above the sum, minus the sum; zero when the sum is above
    every normal valence, which is how a hypervalent atom gets no hydrogens rather than a negative
    count.
    """
    cdef uint32_t v
    for v in range(order_sum, 16):
        if mask & (1u << v):
            return v - order_sum
    return 0


cdef inline bint smv_default_h(uint32_t element, uint32_t order_sum,
                               uint32_t *h_out) noexcept nogil:
    """Can this element's hydrogen count be left implicit, and what would a reader infer?

    Returns True with `h_out` set when the element is in the organic subset and both valence models
    above infer the same count; False when the element must be bracketed regardless -- either it is
    outside the subset, or the two models disagree at this bond-order sum and the count is not safe
    to omit.  The reader epic calls THIS FUNCTION rather than restating the numbers, so the two
    directions cannot drift (spec §9 item 2).
    """
    # Initialised although `smv_valences` writes both unconditionally: Cython cannot see
    # through a pointer out-parameter, so without this the maybe-uninitialized warning is a false
    # positive that would have to be waived, and the file's gate is zero warnings.
    cdef uint32_t narrow = 0, wide = 0, hn, hw
    smv_valences(element, &narrow, &wide)
    if narrow == 0:
        return False
    hn = smv_h_from_mask(narrow, order_sum)
    hw = smv_h_from_mask(wide, order_sum)
    if hn != hw:
        return False
    h_out[0] = hn
    return True


cdef inline bint smv_aromatic_h(uint32_t element, uint32_t order_sum,
                                uint32_t *h_out) noexcept nogil:
    """The same question for an atom written LOWERCASE, where the model above cannot be asked.

    `smv_default_h` takes a Kekule bond-order sum, and an aromatic bond has no order in that model:
    feeding it the stored 4 makes benzene's carbon look like a sum of 8, which the model reads as
    hypervalent and answers with zero hydrogens.  That is how a hand-built aromatic benzene came out
    as bare `C` while its Kekule twin came out as `[C]` -- the same atoms, two answers, because the
    question was nonsense in one of them.  So the aromatic case gets its own rule here rather than a
    new entry in `smv_valences`, which is the chemistry-facing table MDL and InChI consume.

    THE RULE, and `order_sum` must arrive with every aromatic bond counted as 1: add one for the
    atom's own share of the ring's pi system and subtract from the element's LOWEST normal valence.
    Lowest rather than smallest-at-or-above, which is what the Kekule rule uses, and thiophene is
    why: sulfur's normal valences are 2, 4 and 6, its aromatic sum is 3, and smallest-at-or-above
    would infer one hydrogen on an `s` that has none.  Taking the lowest and clamping at zero gives
    0 for `s` and for furan's `o`, 1 for benzene's `c`, 0 for pyridine's `n` and 0 for naphthalene's
    fusion carbons -- every one of them measured against RDKit 2026.03.4, and pyrrole's `n` falls
    out as a MISMATCH, which is exactly why pyrrole is spelled `[nH]` by everybody.

    Restricted to B, C, N, O, P and S: those are SMILES' aromatic organic subset, and an aromatic
    atom outside it (`se`, `as`, or a garbage aromatic chlorine the arena will happily store) gets
    False, so it is bracketed and its count is stated rather than guessed.
    """
    cdef uint32_t narrow = 0, wide = 0, v
    if not (element == 5 or element == 6 or element == 7
            or element == 8 or element == 15 or element == 16):
        return False
    smv_valences(element, &narrow, &wide)
    for v in range(16):
        if narrow & (1u << v):
            h_out[0] = (v - order_sum) if v > order_sum else 0
            return True
    return False


def smv_valence_model():
    """The valence model as {atomic number: (narrow valences, wide valences)}, for the test suite.

    Also the reference the reader epic reads: two tuples per element, ascending, and an element
    absent from this dict is one that is always bracketed.
    """
    cdef uint32_t element, v
    cdef uint32_t narrow = 0, wide = 0     # see smv_default_h
    cdef dict out = {}
    cdef list ln, lw
    for element in range(1, 119):
        smv_valences(element, &narrow, &wide)
        if narrow == 0:
            continue
        ln = []
        lw = []
        for v in range(16):
            if narrow & (1u << v):
                ln.append(v)
            if wide & (1u << v):
                lw.append(v)
        out[element] = (tuple(ln), tuple(lw))
    return out


# ------------------------------------------------------------------------------------------------
# THE OUTPUT BUFFER.  A growable char block with a STICKY out-of-memory flag rather than an error
# return on every put: the alternative is an `except -1` on a dozen one-line appenders and a check
# at each of forty call sites, all for a failure that is checked once at the end either way.  A put
# after the flag is set is a no-op, so a failed grow cannot make a later put write out of bounds.
cdef inline bint smw_reserve(smw_buf_t *b, size_t extra) noexcept:
    cdef size_t want
    cdef char *grown
    if b.oom:
        return False
    if b.length + extra <= b.cap:
        return True
    want = b.cap
    if want < SMW_BUF_MIN:
        want = SMW_BUF_MIN
    while want < b.length + extra:
        want *= 2
    grown = <char *> PyMem_Realloc(b.data, want)
    if grown is NULL:
        b.oom = True
        return False
    b.data = grown
    b.cap = want
    return True


cdef inline void smw_putc(smw_buf_t *b, char c) noexcept:
    if not smw_reserve(b, 1):
        return
    b.data[b.length] = c
    b.length += 1


cdef inline void smw_puts(smw_buf_t *b, const char *s, size_t k) noexcept:
    if not smw_reserve(b, k):
        return
    memcpy(b.data + b.length, s, k)
    b.length += k


cdef inline void smw_putu(smw_buf_t *b, uint32_t v) noexcept:
    """An unsigned decimal, with no snprintf: the values are element symbols' lengths away from
    huge and the buffer's growth is already handled."""
    cdef char tmp[12]
    cdef uint32_t k = 0
    if v == 0:
        smw_putc(b, c'0')
        return
    while v:
        tmp[k] = <char> (c'0' + (v % 10))
        v //= 10
        k += 1
    if not smw_reserve(b, k):
        return
    while k:
        k -= 1
        b.data[b.length] = tmp[k]
        b.length += 1


# ------------------------------------------------------------------------------------------------
# SCRATCH.  One struct, one malloc, one failure check, one free (RULES.md §5.2), sized with
# align8() so every uint32_t block is aligned however the uint8_t blocks fall.
cdef int smw_scratch_alloc(smw_scratch_t *s, uint32_t n,
                           uint32_t nhe) except -1:
    cdef size_t u32 = align8(<size_t> (11 * <size_t> n + 1) * sizeof(uint32_t)) \
                      + align8(<size_t> (4 * <size_t> nhe + 1) * sizeof(uint32_t))
    cdef size_t u8 = align8(<size_t> (3 * <size_t> nhe + 1) * sizeof(uint8_t)) \
                     + align8(<size_t> (5 * <size_t> n + 1) * sizeof(uint8_t))
    cdef char *p
    s.block = PyMem_Malloc(u32 + u8)
    if s.block is NULL:
        raise MemoryError()
    memset(s.block, 0, u32 + u8)
    s.n = n
    s.nhe = nhe
    s.nseq = 0
    p = <char *> s.block
    s.pos = <uint32_t *> p
    s.bypos = s.pos + n
    s.seq = s.bypos + n
    s.out_idx = s.seq + n
    s.parent = s.out_idx + n
    s.stack = s.parent + n
    s.cursor = s.stack + n
    s.unit_at = s.cursor + n
    s.partner_at = s.unit_at + n
    s.path_rank = s.partner_at + n
    s.nbr_off = s.path_rank + n             # n + 1 entries: the +1 is why the block carries 11n+1
    s.nbr = <uint32_t *> (p + align8(<size_t> (11 * <size_t> n + 1) * sizeof(uint32_t)))
    s.wnbr = s.nbr + nhe
    s.dq = s.wnbr + nhe                     # 2 * nhe: one (from slot, half-edge) pair per half-edge
    p = p + u32
    s.he_kind = <uint8_t *> p
    s.he_close = s.he_kind + nhe
    s.he_dir = s.he_close + nhe
    p = p + align8(<size_t> (3 * <size_t> nhe + 1) * sizeof(uint8_t))
    s.visited = <uint8_t *> p
    s.paren = s.visited + n
    s.lost = s.paren + n
    s.dropped = s.lost + n
    s.frame_cut = s.dropped + n
    s.nlost = 0
    return 0


cdef inline void smw_scratch_free(smw_scratch_t *s) noexcept:
    PyMem_Free(s.block)
    s.block = NULL


# ------------------------------------------------------------------------------------------------
# THE CUTS (spec §13).  Everything that turns "keep this side of that bond" into a `dropped` mask and a
# set of attachment half-edges, INCLUDING every refusal -- so a caller's mistake is a Python exception
# raised before a single character is written, naming stable ids they can act on.
#
# ONE RULE, TWO REFUSALS: THE CUT LIST MUST BE EXACTLY THE EDGE BOUNDARY BETWEEN RETAINED AND DROPPED.
# `dropped` is grown from the DROP atoms rather than from the keeps, which is not the same thing and is
# the reason a salt keeps its counter-ion: an unrelated component contains no drop atom, so nothing
# reaches it, so it stays.  Growing from the keeps would drop every component the caller did not name.
cdef int smw_apply_cuts(MoleculeContainer molecule, Structure structure, smw_scratch_t *s,
                        smw_cuts_t *c) except -1:
    """Mark the cut half-edges `SMW_HE_ATTACH`, their ids into `he_close`, and fill `dropped`.

    Runs before the canonical order, so a refusal costs nothing: the expensive part of a write is the
    order and there is no point computing one for a cut list that cannot be honoured.
    """
    cdef uint32_t *ptr = csr_ptr(structure)
    cdef halfedge_t *edges = csr_edges(structure)
    cdef list numbers = molecule._numbers
    cdef halfedge_t *e
    cdef uint32_t i, k, u, v, he, rev, head, tail, retained = 0
    cdef list path

    for i in range(c.ncuts):
        u = c.keep[i]
        v = c.drop[i]
        if c.ids[i] < SMW_MIN_ATTACH or c.ids[i] > SMW_MAX_CLOSURE:
            raise ValueError('attachment id %d is outside %d..%d; below %d it would be written as a '
                             'bare digit and be indistinguishable from an ordinary ring closure'
                             % (c.ids[i], SMW_MIN_ATTACH, SMW_MAX_CLOSURE, SMW_MIN_ATTACH))
        c.reserved[c.ids[i]] = 1
        if u == v:
            raise ValueError('cut %d names atom %r as both the retained and the dropped side'
                             % (c.ids[i], numbers[u]))
        e = csr_find_at(ptr, edges, u, v)
        if e is NULL:
            raise ValueError('cut %d names atoms %r and %r, which are not bonded'
                             % (c.ids[i], numbers[u], numbers[v]))
        he = <uint32_t> (e - edges)
        rev = <uint32_t> (csr_find_at(ptr, edges, v, u) - edges)
        if s.he_kind[he] == SMW_HE_ATTACH:
            raise ValueError('the bond %r-%r is named by two cuts (%d and %d)'
                             % (numbers[u], numbers[v], s.he_close[he] or s.he_close[rev], c.ids[i]))
        s.he_kind[he] = SMW_HE_ATTACH
        s.he_kind[rev] = SMW_HE_ATTACH
        s.he_close[he] = c.ids[i]           # the retained side carries the number
    # A drop atom named as a keep by ANOTHER cut is a contradiction, and it is checked before the walk
    # so that the walk's own answer cannot be blamed for it.
    for i in range(c.ncuts):
        for k in range(c.ncuts):
            if c.drop[i] == c.keep[k]:
                raise ValueError('atom %r is the dropped side of cut %d and the retained side of '
                                 'cut %d' % (numbers[c.drop[i]], c.ids[i], c.ids[k]))

    # DROPPED = everything reachable from a drop atom without crossing a cut.  Breadth-first with
    # `parent` recorded, because the interesting failure -- a cut that comes back on itself -- is best
    # explained by the path that came back, and `smw_traverse` reinitialises `parent` afterwards.
    head = 0
    tail = 0
    for i in range(s.n):
        s.parent[i] = SMW_NONE
    for i in range(c.ncuts):
        if not s.dropped[c.drop[i]]:
            s.dropped[c.drop[i]] = 1
            s.stack[tail] = c.drop[i]
            tail += 1
    while head < tail:
        u = s.stack[head]
        head += 1
        for k in range(ptr[u], ptr[u + 1]):
            if s.he_kind[k] == SMW_HE_ATTACH:
                continue
            v = edges[k].to
            if s.dropped[v]:
                continue
            s.dropped[v] = 1
            s.parent[v] = u
            s.stack[tail] = v
            tail += 1

    for i in range(c.ncuts):
        if s.dropped[c.keep[i]]:
            # A RING CUT, or a dropped fragment that hangs on by a second bond.  Dropping `drop` alone
            # would leave a valid fragment with TWO open valences, and one attachment id cannot carry
            # both ends -- an id written twice in one fragment is a ring closure, so it would silently
            # re-form the very ring the caller asked to break.  The path is in the message because the
            # fix is a second cut somewhere along it.
            path = []
            u = c.keep[i]
            while u != SMW_NONE:
                path.append(numbers[u])
                u = s.parent[u]
            raise ValueError('cut %d cannot be made: %r is still bonded to %r through %s, so the '
                             'bond is in a ring.  One attachment id cannot carry both ends of a ring '
                             'opening; cut a second bond along that path'
                             % (c.ids[i], numbers[c.keep[i]], numbers[c.drop[i]],
                                '-'.join(map(repr, path))))
    # THE INVARIANT, CHECKED RATHER THAN ASSUMED -- and both of its failures are UNREACHABLE given the
    # refusal above, by this function's own construction rather than by anything another mechanism
    # does.  The proof is two lines: `dropped` is closed under non-attachment edges, so a retained atom
    # with an unnamed edge into the dropped part would itself have been reached by the walk and be
    # dropped; and a `keep` atom is retained or the ring refusal already fired, so the retained set is
    # non-empty whenever there is one cut, and non-empty trivially when there is none.  Kept because
    # both proofs are about THE WALK: change the walk -- grow the set from the keeps instead, which is
    # the formulation that silently drops a mixture's counter-ion -- and these become live.
    for u in range(s.n):
        if s.dropped[u]:
            continue
        retained += 1
        for k in range(ptr[u], ptr[u + 1]):
            v = edges[k].to
            if s.dropped[v] and s.he_kind[k] != SMW_HE_ATTACH:
                raise ValueError('the bond %r-%r crosses into the dropped part and no cut names it; '
                                 'add (%r, %r) to the cuts' % (numbers[u], numbers[v], numbers[u], numbers[v]))
    if retained == 0:
        raise ValueError('every atom would be dropped; a fragment needs at least one atom')
    return 0


# ------------------------------------------------------------------------------------------------
# THE TRAVERSAL.
cdef void smw_sort_adjacency(Structure structure, smw_scratch_t *s) noexcept nogil:
    """Fill `nbr_off` and `nbr`: each atom's half-edge indices ascending by `pos` of the far atom.

    An insertion sort per atom, on `pos` and never on the slot index: this is decision 3 at the top
    of the file, and a slot tie-break here would be exactly the defect the design closes.  `pos` is
    a permutation, so no two neighbours tie and the order is total.
    """
    cdef uint32_t *ptr = csr_ptr(structure)
    cdef halfedge_t *edges = csr_edges(structure)
    cdef uint32_t i, k, j, he, key
    for i in range(s.n + 1):
        s.nbr_off[i] = ptr[i]
    for i in range(s.n):
        for k in range(ptr[i], ptr[i + 1]):
            he = k
            key = s.pos[edges[he].to]
            j = k
            while j > ptr[i] and s.pos[edges[s.nbr[j - 1]].to] > key:
                s.nbr[j] = s.nbr[j - 1]
                j -= 1
            s.nbr[j] = he


cdef void smw_traverse(Structure structure, smw_scratch_t *s) noexcept nogil:
    """Classify every half-edge and fill `seq`, `parent`, `out_idx` and `wnbr`.

    Components are entered in ascending minimum `pos`, each is entered at its own minimum-`pos`
    atom, and neighbours are taken in `nbr` order -- so the whole traversal is a function of `pos`.
    A depth-first walk with an explicit stack rather than recursion, because the recursion depth is
    the molecule's longest path and a 40,000-atom polymer would overflow the C stack.

    A DROPPED atom (spec §13) is never entered, and the walk cannot reach one either: the cut
    half-edges were classified `SMW_HE_ATTACH` before this ran, and an already-classified edge is
    skipped below.  So the mask is consulted at the component starts only, and `dropped` is all zero
    for a whole-molecule write -- one branch, no second traversal.
    """
    cdef uint32_t *ptr = csr_ptr(structure)
    cdef halfedge_t *edges = csr_edges(structure)
    cdef uint32_t p, i, u, v, he, rev, top
    for i in range(s.n):
        s.parent[i] = SMW_NONE
    for p in range(s.n):
        i = s.bypos[p]
        if s.visited[i] or s.dropped[i]:
            continue
        s.visited[i] = 1
        s.parent[i] = SMW_NONE
        s.cursor[i] = s.nbr_off[i]
        s.out_idx[i] = s.nseq
        s.seq[s.nseq] = i
        s.nseq += 1
        s.stack[0] = i
        top = 1
        while top:
            u = s.stack[top - 1]
            if s.cursor[u] < s.nbr_off[u + 1]:
                he = s.nbr[s.cursor[u]]
                s.cursor[u] += 1
                if s.he_kind[he] != SMW_HE_UNSEEN:
                    continue                # classified from the other side already
                v = edges[he].to
                rev = <uint32_t> (csr_find_at(ptr, edges, v, u) - edges)
                if s.visited[v]:
                    s.he_kind[he] = SMW_HE_CLOSURE
                    s.he_kind[rev] = SMW_HE_CLOSURE
                else:
                    s.he_kind[he] = SMW_HE_CHILD
                    s.he_kind[rev] = SMW_HE_PARENT
                    s.visited[v] = 1
                    s.parent[v] = u
                    s.cursor[v] = s.nbr_off[v]
                    s.out_idx[v] = s.nseq
                    s.seq[s.nseq] = v
                    s.nseq += 1
                    s.stack[top] = v
                    top += 1
            else:
                top -= 1

    smw_written_order(structure, s)


cdef inline void smw_sort_children(Structure structure, smw_scratch_t *s,
                                   uint32_t start, uint32_t stop) noexcept nogil:
    """Sort `wnbr[start:stop]` -- one atom's child edges -- ascending by the child's `out_idx`.

    An insertion sort, because the range is one atom's degree.  `out_idx` is a permutation of the
    written positions, so no two children tie and the order is total.
    """
    cdef halfedge_t *edges = csr_edges(structure)
    cdef uint32_t i, j, he, key
    for i in range(start + 1, stop):
        he = s.wnbr[i]
        key = s.out_idx[edges[he].to]
        j = i
        while j > start and s.out_idx[edges[s.wnbr[j - 1]].to] > key:
            s.wnbr[j] = s.wnbr[j - 1]
            j -= 1
        s.wnbr[j] = he


cdef void smw_written_order(Structure structure, smw_scratch_t *s) noexcept nogil:
    """Fill `wnbr` from the classified half-edges.  Every traversal ends here, and only here.

    WRITTEN order per atom: the parent first, then the ring bonds, then the children, each group
    keeping its `nbr` order.  The groups are separated because that is what the SMILES grammar
    says -- `atom ringbond* branch*` -- and the same list is the neighbour order a tetrahedral
    sign is read in, so the two cannot disagree.

    AN ATTACHMENT IS A RING BOND, in the grammar and here: `%12` after the atom, in the same group
    and taking its number from the same pool.  That identity is not a convenience -- it is what
    makes two fragments re-join by concatenation, since the reader has no other production that
    could accept a dangling number.
    """
    cdef uint32_t *ptr = csr_ptr(structure)
    cdef uint32_t i, j, k, w
    for i in range(s.n):
        w = s.nbr_off[i]
        for k in range(ptr[i], ptr[i + 1]):
            if s.he_kind[s.nbr[k]] == SMW_HE_PARENT:
                s.wnbr[w] = s.nbr[k]
                w += 1
        for k in range(ptr[i], ptr[i + 1]):
            if s.he_kind[s.nbr[k]] == SMW_HE_CLOSURE or s.he_kind[s.nbr[k]] == SMW_HE_ATTACH:
                s.wnbr[w] = s.nbr[k]
                w += 1
        # THE CHILDREN IN THE ORDER THEY WERE DISCOVERED, which is `out_idx` and NOT `nbr` order.  For
        # an ordinary traversal the two are the same thing -- children are entered in `nbr` order, so
        # the first in that order gets the smaller `out_idx` -- and this insertion sort does nothing.
        # A STICKY TRAVERSAL DEFERS ONE CHILD TO LAST, so there they differ, and `wnbr` has to follow
        # the walk: it is what `smw_emit` writes in, what `smw_closures` numbers in, and what a
        # tetrahedral sign is read in.  Sorting on the walk's own answer is what keeps all three equal to
        # each other; `nbr` order would silently make the string end at the wrong atom.
        j = w
        for k in range(ptr[i], ptr[i + 1]):
            if s.he_kind[s.nbr[k]] == SMW_HE_CHILD:
                s.wnbr[w] = s.nbr[k]
                w += 1
        smw_sort_children(structure, s, j, w)
        # A DROPPED atom's remaining edges are still `SMW_HE_UNSEEN`, and leaving the tail of its
        # `wnbr` window unwritten would leave stale indices there.  Nothing reads a dropped atom's
        # window, but "nothing reads it" is a property of other functions, and an uninitialised
        # window is the shape of bug this file has already been bitten by once.
        for k in range(ptr[i], ptr[i + 1]):
            if s.he_kind[s.nbr[k]] == SMW_HE_UNSEEN:
                s.wnbr[w] = s.nbr[k]
                w += 1


# ------------------------------------------------------------------------------------------------
# THE STICKY TRAVERSAL (§14).  A walk that STARTS at one named atom and ENDS at another, by
# construction rather than by retrying a randomised order until one lands.
cdef uint32_t smw_sticky_start(smw_scratch_t *s, smw_sticky_t *k) noexcept nogil:
    """The atom the string starts at: `left` when the caller named one, else the smallest-`pos` atom
    that is not `right`.

    Not simply `bypos[0]`, because that atom may BE `right` -- and then the walk would have to both
    begin and end at it.  Skipping it keeps the right-only call a function of
    `pos` without making it a special case anywhere else.  A one-atom molecule has no other atom, and
    then `start == right` is correct: one token is both the first and the last.
    """
    cdef uint32_t p
    if k.left != SMW_NONE:
        return k.left
    for p in range(s.n):
        if s.bypos[p] != k.right:
            return s.bypos[p]
    return s.bypos[0]


cdef void smw_sticky_path(Structure structure, smw_scratch_t *s, smw_sticky_t *k) noexcept nogil:
    """Fill `path_rank`: each atom's index on one shortest start..`right` path, SMW_NONE off it.

    A breadth-first pass from `right` for the distances, then a walk back from the start atom taking
    the SMALLEST-`pos` neighbour one step closer -- smallest-`pos` so that the path, and therefore
    the whole string, is a function of `pos` and the two named atoms and of nothing else.

    `cursor` holds the distances.  It is the DFS's per-frame cursor afterwards, and the DFS writes
    every frame's entry before reading it, so the two uses cannot collide; `stack` is the queue here
    for the same reason.  Both are why this needs no allocation of its own.
    """
    cdef uint32_t *ptr = csr_ptr(structure)
    cdef halfedge_t *edges = csr_edges(structure)
    cdef uint32_t i, j, u, v, head = 0, tail = 0, best, rank = 0
    for i in range(s.n):
        s.path_rank[i] = SMW_NONE
        s.cursor[i] = SMW_NONE
    if k.right == SMW_NONE:
        return
    s.cursor[k.right] = 0
    s.stack[0] = k.right
    tail = 1
    while head < tail:
        u = s.stack[head]
        head += 1
        for j in range(ptr[u], ptr[u + 1]):
            v = edges[j].to
            if s.cursor[v] == SMW_NONE:
                s.cursor[v] = s.cursor[u] + 1
                s.stack[tail] = v
                tail += 1
    u = smw_sticky_start(s, k)
    s.path_rank[u] = 0
    # The caller has already refused a `right` in another component (`smw_sticky_resolve`), so the
    # start atom has a finite distance and every step below finds a predecessor.  `SMW_NONE` on
    # `best` would be that refusal having failed, and the loop stops rather than walking off the
    # array -- a `noexcept nogil` function cannot report, and a silent stop leaves `right` merely
    # unforced, which the entry point's own check then catches.
    while u != k.right and s.cursor[u] != SMW_NONE:
        best = SMW_NONE
        for j in range(ptr[u], ptr[u + 1]):
            v = edges[j].to
            if s.cursor[v] == s.cursor[u] - 1 and (best == SMW_NONE or s.pos[v] < s.pos[best]):
                best = v
        if best == SMW_NONE:
            break
        u = best
        rank += 1
        s.path_rank[u] = rank


cdef void smw_sticky_traverse(Structure structure, smw_scratch_t *s, smw_sticky_t *k) noexcept nogil:
    """`smw_traverse`, with the first atom named and the last one forced.

    TWO RULES ON TOP OF THE DEPTH-FIRST WALK, and they are not the same rule:

    * an edge into an UNVISITED path atom is BLOCKED unless it is this atom's own successor.  The
      half-edge is left `SMW_HE_UNSEEN` and is classified later from the far side, which is what the
      classification's "from whichever side reaches it first" already allows for.
    * the successor edge is taken LAST, after everything else at this atom is exhausted.

    Blocking is the half that matters, and DEFERRING ALONE IS NOT ENOUGH.  `left-a-b-right` with a
    ring `a-c-b` and a pendant `c-e`: defer `b` at `a`, the walk enters `c`, `c` reaches `b` first,
    `b` defers `right` -- and `e` is written after `right`.  With `b` blocked from `c`, `c` finishes
    `e` first and `b` is entered only from `a`.

    `right` is then last, by induction on the path index: arriving at `p_i` every earlier path atom
    is visited and every later one is unreachable, so the side subtrees cannot discover one; the
    successor goes last; so `p_k == right` is discovered last of all.  And every atom is reached,
    because for any atom the last path atom on a walk to it from the start is some `p_m` and the rest
    of that walk crosses no path atom.  NOTHING IN THAT ARGUMENT NEEDS `right` TO BE TERMINAL -- that
    is a constraint of removing its TOKEN, not of the walk.
    """
    cdef uint32_t *ptr = csr_ptr(structure)
    cdef halfedge_t *edges = csr_edges(structure)
    cdef uint32_t p, i, u, v, he, rev, top, j, cand
    for i in range(s.n):
        s.parent[i] = SMW_NONE
    # The named atom's component FIRST and the rest in the ordinary order after it, so that `left` is
    # the first token of the whole string.  A `right` in another component is refused before this
    # runs, so the loop below only ever adds components that hold neither end.
    for p in range(s.n + 1):
        if p == 0:
            i = smw_sticky_start(s, k)
        else:
            i = s.bypos[p - 1]
            if s.visited[i]:
                continue
        s.visited[i] = 1
        s.parent[i] = SMW_NONE
        s.cursor[i] = s.nbr_off[i]
        s.out_idx[i] = s.nseq
        s.seq[s.nseq] = i
        s.nseq += 1
        s.stack[0] = i
        top = 1
        while top:
            u = s.stack[top - 1]
            he = SMW_NONE
            while s.cursor[u] < s.nbr_off[u + 1]:
                cand = s.nbr[s.cursor[u]]
                s.cursor[u] += 1
                if s.he_kind[cand] != SMW_HE_UNSEEN:
                    continue
                if s.path_rank[edges[cand].to] != SMW_NONE and not s.visited[edges[cand].to]:
                    continue            # blocked: a path atom is entered from its predecessor only
                he = cand
                break
            if he == SMW_NONE and s.path_rank[u] != SMW_NONE:
                # Everything else at this atom is written; now the successor.  Searched over `ptr`
                # rather than resumed from `cursor`, because the cursor has already walked past it.
                for j in range(ptr[u], ptr[u + 1]):
                    cand = s.nbr[j]
                    if s.he_kind[cand] != SMW_HE_UNSEEN:
                        continue
                    if s.path_rank[edges[cand].to] == s.path_rank[u] + 1:
                        he = cand
                        break
            if he == SMW_NONE:
                top -= 1
                continue
            v = edges[he].to
            rev = <uint32_t> (csr_find_at(ptr, edges, v, u) - edges)
            if s.visited[v]:
                s.he_kind[he] = SMW_HE_CLOSURE
                s.he_kind[rev] = SMW_HE_CLOSURE
            else:
                s.he_kind[he] = SMW_HE_CHILD
                s.he_kind[rev] = SMW_HE_PARENT
                s.visited[v] = 1
                s.parent[v] = u
                s.cursor[v] = s.nbr_off[v]
                s.out_idx[v] = s.nseq
                s.seq[s.nseq] = v
                s.nseq += 1
                s.stack[top] = v
                top += 1
    smw_written_order(structure, s)


cdef bint smw_sticky_severs(Structure structure, smw_scratch_t *s, uint32_t right) noexcept nogil:
    """Whether `right` is a CUT VERTEX, which is exactly when no walk can end there.

    A depth-first walk ends at `right` only if every other atom is discovered first, and an atom whose
    every route from the start passes through `right` cannot be.  So the question "is this end
    reachable last" is the question "does removing this atom disconnect what is left", and it is a
    property of the MOLECULE rather than of the walk -- which is why the entry point refuses on it and
    says so, instead of the walk failing and being retried.  Demanding a TERMINAL `right` is the
    sufficient condition; this is the necessary one, so an in-ring atom can be a sticky end.

    Breadth-first from one of `right`'s neighbours, over the whole molecule because the entry point has
    already refused a second component.  `visited` and `stack` are borrowed BEFORE the traversal owns
    them and `visited` is cleared by the caller; `smw_sticky_traverse` reads it as all-zero.
    """
    cdef uint32_t *ptr = csr_ptr(structure)
    cdef halfedge_t *edges = csr_edges(structure)
    cdef uint32_t i, u, v, head = 0, tail = 0, seen = 1
    if s.n < 3 or ptr[right + 1] == ptr[right]:
        # Two atoms cannot be separated by removing one of them, and a lone atom has nothing to sever.
        return False
    s.visited[right] = 1
    s.visited[edges[ptr[right]].to] = 1
    s.stack[0] = edges[ptr[right]].to
    tail = 1
    while head < tail:
        u = s.stack[head]
        head += 1
        for i in range(ptr[u], ptr[u + 1]):
            v = edges[i].to
            if s.visited[v]:
                continue
            s.visited[v] = 1
            seen += 1
            s.stack[tail] = v
            tail += 1
    return seen < s.n - 1


cdef void smw_sticky_frames(Structure structure, smw_scratch_t *s, smw_sticky_t *k) noexcept nogil:
    """Mark `frame_cut` at the neighbour of every end whose bond is removed OUTRIGHT.

    THE ONE STEREO DECISION A STICKY WRITE MAKES, and it is a refusal.  With `keep_bond_left=True` the
    removed atom's bond token is still written and a caller glues an atom onto it, so the neighbour's
    written order is the SAME sequence of positions before and after the glue -- the predecessor slot
    holds the removed atom now and the glued atom later -- and every sign keeps its meaning.  Nothing
    to do, and that is why this function only looks at the other case.

    With `keep_bond_left=False` the bond goes too, so the neighbour really has one fewer neighbour and
    the string shows it with one fewer position.  A stored parity describes FOUR directions; writing it
    onto an atom the string shows with three would state a configuration the record does not hold.  So
    the sign is refused (`smw_sign_of` returns 0 on this mark) and the atom is REPORTED -- `lost`, by
    `smw_directions` asking `smw_sign_of` the same question it asks about every other tetrahedral atom.
    The shape is `smw_h_frame_unknown`'s, which exists for the same reason.

    MUST RUN BEFORE `smw_directions`, so that its report sees the refusal rather than the sign the
    frame would have had.
    """
    cdef uint32_t *ptr = csr_ptr(structure)
    cdef halfedge_t *edges = csr_edges(structure)
    if k.left != SMW_NONE and k.remove_left and not k.keep_bond_left \
            and ptr[k.left + 1] > ptr[k.left]:
        s.frame_cut[edges[ptr[k.left]].to] = 1
    if k.right != SMW_NONE and k.remove_right and not k.keep_bond_right \
            and ptr[k.right + 1] > ptr[k.right]:
        s.frame_cut[edges[ptr[k.right]].to] = 1


cdef void smw_sticky_ends(Structure structure, smw_scratch_t *s, smw_sticky_t *k) noexcept nogil:
    """Clear the direction token on a bond that is removed outright, and report the unit that needed it.

    `/` and `\\` are the only tokens whose meaning is spread over two bonds, so the one case
    `smw_sticky_frames` cannot handle by refusing a sign is a cis/trans unit with a token on the bond
    that is going away.  Cleared rather than left set, because the token would otherwise be written by
    the OTHER end of that bond -- and reported at the unit's anchor, because a string with one `/` in it
    does not say which configuration was meant.

    A KEPT bond needs none of this: `smw_sticky_bond` writes the `/` itself, and a caller who glues an
    atom in front of it gets the token in exactly the position the unit put it.

    MUST RUN AFTER `smw_directions`, which is what fills `he_dir`.
    """
    cdef uint32_t *ptr = csr_ptr(structure)
    cdef halfedge_t *edges = csr_edges(structure)
    cdef uint32_t i, e, he, rev, far, anchor
    for i in range(2):
        e = k.left if i == 0 else k.right
        if e == SMW_NONE or ptr[e + 1] == ptr[e]:
            continue
        if i == 0:
            if not k.remove_left or k.keep_bond_left:
                continue
        elif not k.remove_right or k.keep_bond_right:
            continue
        he = ptr[e]
        far = edges[he].to
        rev = <uint32_t> (csr_find_at(ptr, edges, far, e) - edges)
        if not s.he_dir[he] and not s.he_dir[rev]:
            continue
        s.he_dir[he] = 0
        s.he_dir[rev] = 0
        # The direction sat on a single bond next to a double-bond terminal, and `e` has degree 1, so
        # the terminal is `far` and the unit it terminates is the one that has just lost a reference.
        anchor = s.unit_at[far]
        if anchor != SMW_NONE and not s.lost[anchor]:
            s.lost[anchor] = 1
            s.nlost += 1


cdef int smw_closures(Structure structure, smw_scratch_t *s, smw_cuts_t *cuts) except -1:
    """Assign ring-closure numbers into `he_close`, in written order.

    A number is taken when a closure OPENS and returned when it CLOSES, but the return is deferred
    to the end of the closing atom's list: an atom that closes 1 and immediately reopens 1 would
    write `C11`, which reads as closure 11.

    Smallest-free-first, so the numbers stay small and are a function of the written order -- which
    is a function of `pos`.

    ATTACHMENT AND RESERVED IDS ARE WITHHELD FOR THE WHOLE WRITE, never released: `%12` and `12` are
    one ring bond to every reader, so an internal closure that reused an attachment's number would
    bond to the wrong atom the moment two fragments were concatenated -- and it would do so silently,
    producing a valid molecule that is not the one anybody asked for.  Withholding for the whole
    string rather than for the attachment's lifetime is deliberate: a fragment's ids must mean the
    same thing at every position in it, because the join happens at the string level and knows
    nothing about where a number was in scope.
    """
    cdef uint32_t *ptr = csr_ptr(structure)
    cdef halfedge_t *edges = csr_edges(structure)
    cdef uint8_t inuse[SMW_MAX_CLOSURE + 1]
    cdef uint32_t released[SMW_MAX_CLOSURE + 1]
    cdef uint32_t idx, u, k, he, rev, c, nrel, i
    memset(inuse, 0, sizeof(inuse))
    if cuts is not NULL:
        for i in range(SMW_MIN_ATTACH, SMW_MAX_CLOSURE + 1):
            if cuts.reserved[i]:
                inuse[i] = 1
    for idx in range(s.nseq):
        u = s.seq[idx]
        nrel = 0
        for k in range(s.nbr_off[u], s.nbr_off[u + 1]):
            he = s.wnbr[k]
            if s.he_kind[he] != SMW_HE_CLOSURE:
                continue
            rev = <uint32_t> (csr_find_at(ptr, edges, edges[he].to, u) - edges)
            if s.he_close[rev]:                 # the far atom opened it; this is the closing end
                s.he_close[he] = s.he_close[rev]
                released[nrel] = s.he_close[he]
                nrel += 1
                continue
            c = 0
            for i in range(1, SMW_MAX_CLOSURE + 1):
                if not inuse[i]:
                    c = i
                    break
            if c == 0:
                raise ValueError('more than %d ring closures are open at once; the SMILES '
                                 'closure numbers are exhausted' % SMW_MAX_CLOSURE)
            inuse[c] = 1
            s.he_close[he] = <uint8_t> c
        for i in range(nrel):
            inuse[released[i]] = 0
    return 0


# ------------------------------------------------------------------------------------------------
# TOKENS.
cdef void smw_bond(smw_buf_t *b, halfedge_t *e, bint both_lower,
                   uint32_t direction, smw_opts_t *o) noexcept:
    """One bond token.  Empty for a plain single bond, which is the whole point of SMILES.

    `both_lower` is whether the two atoms this bond joins WILL BE WRITTEN LOWERCASE, and it decides
    both aromatic cases below.  NO TOKEN IS EVER SUPPRESSED ON A PREDICTION ABOUT ANOTHER MECHANISM --
    the rule is the file's, stated once at the top -- so the aromatic branch asks whether the
    aromaticity was actually carried elsewhere and writes `:` whenever it was not.  An extra `:` on a
    lowercase atom pair is noise; a missing one turns benzene into cyclohexane, and the asymmetry is
    total.
    """
    if not o.bonds:
        return
    if direction == 1:
        smw_putc(b, c'/')
        return
    elif direction == 2:
        smw_putc(b, c'\\')
        return
    if e.order == 2:
        smw_putc(b, c'=')
    elif e.order == 3:
        smw_putc(b, c'#')
    elif e.order == 8:
        # chython's own dialect and not standard SMILES, where `~` is SMARTS' any-bond.  Kept
        # because the arena can hold an order-8 bond and refusing to write the molecule is worse;
        # the reader epic accepts it (spec §9 item 3).
        smw_putc(b, c'~')
    elif e.flags & HE_AROMATIC:
        # `HE_AROMATIC` iff `order == 4`, guaranteed by the arena, so this reads the bond rather than
        # a perception result.  Silent by default because the LOWERCASE ATOMS carry the aromaticity and
        # `c1ccccc1` is what every reader expects -- but silent BECAUSE THEY DID, not because they were
        # expected to.  `A` makes them uppercase and the `:` appears here instead; so would any future
        # atom that could not be spelled lowercase.
        if not both_lower:
            smw_putc(b, c':')
    elif both_lower:
        # A single bond between two lowercase atoms MUST be written, or it reads back aromatic --
        # biphenyl's central bond, and `c1ccccc1c1ccccc1` is a different molecule from
        # `c1ccccc1-c1ccccc1`.  Under `A` the atoms are uppercase, so an unmarked single bond between
        # them cannot be mistaken for part of a ring system and `both_lower` is already False.
        smw_putc(b, c'-')


cdef void smw_sticky_bond(smw_buf_t *b, halfedge_t *e, uint32_t direction) noexcept:
    """A sticky end's KEPT bond token, and never the empty string.

    The whole point of `keep_bond_*`: the atom at this end of the bond is not in the string, so the
    caller glues one on and the token is what says how.  An empty token would be read as a single bond
    by the reader that eventually sees `X` + `-Y`, which is right for order 1 and WRONG FOR AN
    AROMATIC BOND.  Patching an empty token to `-` after the fact cannot tell the two apart; here the
    bond is in hand, so `:` is written.

    `o.bonds` is deliberately not consulted: under `!b` a caller asked for no bond tokens, but this
    token is not decoration -- it is the open valence, and a sticky end without it is a different
    request (`keep_bond_*=False`) that the caller can make directly.

    A DIRECTION WINS OVER THE ORDER TOKEN, and both are single-bond spellings: `/` says single as well
    as which way, so a cis/trans configuration whose token sits on this very bond SURVIVES the glue --
    `F/C=C/` with an `F` glued on is `F/C=C/F`.  The alternative would be writing `-` here and
    reporting the unit lost, which throws away a fact the notation can carry.
    """
    if direction == 1:
        smw_putc(b, c'/')
        return
    elif direction == 2:
        smw_putc(b, c'\\')
        return
    if e.order == 2:
        smw_putc(b, c'=')
    elif e.order == 3:
        smw_putc(b, c'#')
    elif e.order == 8:
        smw_putc(b, c'~')
    elif e.flags & HE_AROMATIC:
        smw_putc(b, c':')
    else:
        smw_putc(b, c'-')


cdef void smw_closure_number(smw_buf_t *b, uint32_t c) noexcept:
    if c >= 10:
        smw_putc(b, c'%')
    smw_putu(b, c)


cdef void smw_charge(smw_buf_t *b, int charge) noexcept:
    """`+`, `-`, `+2`, `-3`: the sign alone for a unit charge and sign-plus-digits above it."""
    if charge > 0:
        smw_putc(b, c'+')
        if charge > 1:
            smw_putu(b, <uint32_t> charge)
    else:
        smw_putc(b, c'-')
        if charge < -1:
            smw_putu(b, <uint32_t> (-charge))


cdef inline uint32_t smw_written_h(atom_t *a) noexcept nogil:
    """How many hydrogens this atom's token SHOWS: its implicit count, or none when that count is
    H_UNKNOWN.

    `at_implicit_h` is raw by the arena's design -- an untaught reader gets a visibly absurd 15 rather
    than a plausible 0 -- so every place in this file that turns the count into characters, or into the
    POSITIONS those characters occupy in the written order, goes through here instead.  Zero and not
    "unstated" because SMILES has no spelling for unstated inside a bracket; `smw_atom` argues that and
    reports the loss.
    """
    return 0 if at_implicit_h_unknown(a) else at_implicit_h(a)


cdef inline bint smw_h_frame_unknown(smw_scratch_t *s, atom_t *a, uint32_t slot,
                                     uint32_t fills) noexcept nogil:
    """True when an UNKNOWN implicit count makes this atom's written order unreadable.

    A four-direction frame can only be translated into a sign if the WRITTEN POSITIONS of its unnamed
    directions are known, and an unknown count says nothing about how many of them are hydrogens.
    `fills` is the degree at which the atom's own bonds leave no room for one -- 4 at a tetrahedral
    centre, 3 at an axial terminal, counting the chain bond it drops -- and at or above it the count
    must be zero whatever the nibble says, so the frame is readable after all.  Worth not refusing:
    a stated configuration on a fully substituted centre is exactly where an unknown count arrives
    from a query format.

    CALLED, AND ITS TRUE BRANCH IS UNREACHABLE.  Perception emits no unit where an unknown count could
    have mattered, which is this predicate spelled on the other side: a unit does arrive for the fully
    substituted centre and the tetrasubstituted axis, this is asked about them, and it answers False
    (measured 2026-09-02).  The refusal and the guard agreeing is exactly the state in which one looks
    redundant, and it stays anyway: deleting it would make the frame's correctness rest on another
    mechanism's current refusal, which is what note 4 at the top of the file forbids.  The failure it
    prevents is not a wrong character but a sign computed from a written order with fewer positions than
    the record has -- a different answer per creation order, which is what the sweeps look for.
    """
    return at_implicit_h_unknown(a) and s.nbr_off[slot + 1] - s.nbr_off[slot] < fills


cdef void smw_atom(smw_buf_t *b, atom_t *a, uint32_t order_sum, bint wild_bond,
                   uint32_t narom, uint32_t sign, smw_opts_t *o) noexcept:
    """One atom token, bracketed exactly when its hydrogen count would not survive being omitted.

    The eight reasons for a bracket are the spec's §5 predicate, in the same order.  Whenever the
    bracket goes on, the implicit hydrogen count is written in FULL -- that pair of rules together is
    the fidelity guarantee: the count is either stated or inferable, never guessed.

    `narom` is how many of this atom's bonds are aromatic and `order_sum` arrives with each of them
    counted as 1, because 4 is not a bond order any reader's valence model knows.  A nonzero `narom`
    means a lowercase symbol -- unless `A` asked for the aromaticity on the bonds instead, in which
    case the atom is written uppercase and ALWAYS bracketed: `C:C` is outside OpenSMILES (an aromatic
    bond between aliphatic atoms), so no rule says what hydrogen count it implies and a reader has to
    pick one.  Stating the count is the only unambiguous answer, so `A` costs brackets on every aromatic
    atom and that is the price of the dialect.

    `sign` is 0, 1 for `@` or 2 for `@@`, and it arrives already translated into the order this
    atom's neighbours are about to be written in (ruling F26; see `smw_sign_of`).
    """
    cdef uint32_t element = a.element
    cdef bint unknown = at_implicit_h_unknown(a)
    cdef uint32_t implicit = smw_written_h(a)
    cdef uint32_t inferred = 0
    cdef bint bare = False
    cdef bint aromatic = narom != 0 and not o.aromatic_bond and not smw_no_lowercase(a)
    cdef const char *csym
    cdef bint brackets

    if not wild_bond:
        if narom == 0:
            bare = smv_default_h(element, order_sum, &inferred)
        elif aromatic:
            bare = smv_aromatic_h(element, order_sum + 1, &inferred)
    # AN UNSTATED IMPLICIT COUNT REMOVES TWO OF THE EIGHT REASONS rather than changing what goes inside
    # the brackets, and the OpenSMILES rule behind that is worth stating: inside brackets an ABSENT H
    # term means ZERO hydrogens, not "unstated".  So there is no spelling of "unknown" in a bracket
    # atom at all -- `[13C]` states a hydrogen-free carbon exactly as `[13CH0]` would.  A bare `C`,
    # meanwhile, means precisely "the reader derives it", which is what the molecule says.  Hence:
    #
    #   * `implicit != inferred` cannot demand brackets: there is no stated count to preserve.
    #   * `h` cannot either.  Its contract is "state the count explicitly", and when there is nothing
    #     to state, bracketing the atom would state a zero -- the one thing it must not do.
    #
    # When some OTHER reason brackets the atom anyway (isotope, charge, radical, map number) the H term
    # is omitted and the string then says zero where the molecule said nothing.  That is a real loss,
    # it is unavoidable in this notation, and `smw_traversal` reports it in `unknown_h` rather than
    # letting it pass silently -- the same treatment as a configuration SMILES cannot spell.
    brackets = (not bare
                or a.isotope != 0
                or (a.charge != 0 and o.charges)
                or at_radical(a)
                or (o.mapping and a.map_number != 0)
                or sign != 0
                or (o.hydrogens and not unknown)
                or (implicit != inferred and not unknown))

    csym = SMW_SYMBOL[element]
    if brackets:
        smw_putc(b, c'[')
    if a.isotope:
        smw_putu(b, a.isotope)
    if aromatic:
        # A lowercase symbol is the aromatic spelling, and the first letter is the only one that
        # changes -- `[se]` and `[as]` are spelled that way too.  An element with no aromatic form at
        # all (a chlorine someone gave an order-4 bond) still comes out lowercase, because the arena
        # stored an aromatic bond on it and refusing to show the caller their own structure is the
        # one thing this writer may not do; `smv_aromatic_h` has already forced the brackets.
        smw_putc(b, <char> (csym[0] + 32))
        if csym[1]:
            smw_putc(b, csym[1])
    else:
        smw_putc(b, csym[0])
        if csym[1]:
            smw_putc(b, csym[1])
    if element == 0 and at_r_index(a):
        smw_putu(b, at_r_index(a))
    if sign == 1:
        smw_putc(b, c'@')
    elif sign == 2:
        smw_putc(b, c'@')
        smw_putc(b, c'@')
    if brackets and implicit:
        smw_putc(b, c'H')
        if implicit > 1:
            smw_putu(b, implicit)
    if a.charge and o.charges:
        smw_charge(b, a.charge)
    if o.mapping and a.map_number:
        smw_putc(b, c':')
        smw_putu(b, a.map_number)
    if brackets:
        smw_putc(b, c']')


# ------------------------------------------------------------------------------------------------
# EMISSION.
cdef inline void smw_atom_env(uint32_t *ptr, halfedge_t *edges, uint32_t u,
                              uint32_t *order_sum, bint *wild,
                              uint32_t *narom) noexcept nogil:
    """Everything the atom token needs from the atom's bonds, in one pass over its half-edges.

    AN AROMATIC BOND COUNTS AS 1 in `order_sum`, not as its stored 4.  The sum exists to be handed
    to a valence model, and no reader's model has a rule for 4 -- passing it through makes benzene's
    carbon look like a sum of 8 and every aromatic atom look hypervalent.  `narom` carries the
    aromaticity separately so `smw_atom` can ask `smv_aromatic_h` instead.

    One pass and three outputs rather than three passes, and a helper rather than the two inline
    copies this replaced: three outputs is where "written twice" stops being cheaper than a call.
    """
    cdef uint32_t k
    order_sum[0] = 0
    narom[0] = 0
    wild[0] = False
    for k in range(ptr[u], ptr[u + 1]):
        if edges[k].order == 4:
            order_sum[0] += 1
            narom[0] += 1
        else:
            order_sum[0] += edges[k].order
            if edges[k].order == 8:
                wild[0] = True


cdef inline bint smw_no_lowercase(atom_t *a) noexcept nogil:
    """Whether this atom's symbol has no lowercase spelling at all.

    Element 0's lowercase is `[r]`, the ring-count query primitive, so a marker keeps its case however
    its bonds are stored.  Read by the atom token's case AND by whether the bond token may be
    suppressed, which is the pair `smw_lowercase`'s docstring requires to agree.
    """
    return a.element == 0


cdef inline bint smw_lowercase(atom_t *atoms, uint32_t *ptr, halfedge_t *edges, uint32_t u,
                               smw_opts_t *o) noexcept nogil:
    """Whether atom `u`'s symbol will be WRITTEN lowercase.

    The same fact `smw_atom` decides its own case on -- an atom with at least one order-4 half-edge,
    unless `A` moved the aromaticity onto the bonds, and unless `smw_no_lowercase` rules it out --
    read from the bonds both times, so the two cannot disagree.  `at_hybridization(a) == 4` is the
    arena's own summary of the same thing and would be two byte reads instead of this scan, and it is
    NOT what either place uses: it is a maintained cache, the case of a letter is a fidelity decision,
    and a decision about what the string carries may not rest on something that could go stale.
    Degrees are small; this is a handful of comparisons per bond token.
    """
    cdef uint32_t k
    if o.aromatic_bond:
        return False
    if smw_no_lowercase(&atoms[u]):
        return False
    for k in range(ptr[u], ptr[u + 1]):
        if edges[k].order == 4:
            return True
    return False


cdef inline bint smw_both_lowercase(atom_t *atoms, uint32_t *ptr, halfedge_t *edges,
                                    uint32_t u, uint32_t v, smw_opts_t *o) noexcept nogil:
    return smw_lowercase(atoms, ptr, edges, u, o) and smw_lowercase(atoms, ptr, edges, v, o)


cdef int smw_emit(Structure structure, smw_scratch_t *s, smw_buf_t *b,
                  smw_opts_t *o, smw_sticky_t *sk) except -1:
    """Walk `wnbr` and write the string.

    The same depth-first shape as the traversal and for the same stack-depth reason.  A child gets
    parentheses when it is not the last child, so `C(N)O` and not `C(N)(O)`.

    `sk` is NULL unless this is a sticky write (§14), and all it does here is SUPPRESS TOKENS: the
    named atom's own token, and its bond token unless the caller kept it.  A suppression and not an
    edit -- the atom is still in the traversal, still holds its written position at its neighbour, and
    still contributes its bond's order to the token that replaces it.  That is what makes the string a
    thing a caller can glue an atom onto and get the molecule back, rather than characters cut off a
    finished string.
    """
    cdef atom_t *atoms = structure.atoms()
    cdef uint32_t *ptr = csr_ptr(structure)
    cdef halfedge_t *edges = csr_edges(structure)
    cdef uint32_t idx, i, u, v, he, top, k, sign
    cdef bint more
    # Initialised although `smw_atom_env` writes all three unconditionally: Cython cannot see through
    # a pointer out-parameter, so without this the maybe-uninitialized warning is a false positive
    # and this file's gate is zero warnings.  Same reason as `smv_default_h`'s locals.
    cdef uint32_t order_sum = 0, narom = 0
    cdef bint wild = False
    cdef uint8_t *emitted = s.visited     # re-used: the traversal's marks are spent by now

    memset(emitted, 0, s.n)
    for i in range(s.n):
        s.paren[i] = 0
    idx = 0
    while idx < s.nseq:
        u = s.seq[idx]
        if emitted[u]:
            idx += 1
            continue
        if idx:
            smw_putc(b, c'.')
        top = 0
        if sk is not NULL and u == sk.left and sk.remove_left:
            # No token, and no ring tokens either: the entry point refuses `remove_left` on an atom
            # whose degree is not 1, so there are none to write.
            pass
        else:
            smw_atom_env(ptr, edges, u, &order_sum, &wild, &narom)
            sign = smw_sign_of(structure, s, u, o)
            smw_atom(b, &atoms[u], order_sum, wild, narom, sign, o)
            smw_ring_tokens(structure, s, b, u, o)
        emitted[u] = 1
        s.cursor[u] = s.nbr_off[u]
        s.stack[0] = u
        top = 1
        while top:
            u = s.stack[top - 1]
            he = SMW_NONE
            while s.cursor[u] < s.nbr_off[u + 1]:
                if s.he_kind[s.wnbr[s.cursor[u]]] == SMW_HE_CHILD:
                    he = s.wnbr[s.cursor[u]]
                    s.cursor[u] += 1
                    break
                s.cursor[u] += 1
            if he == SMW_NONE:
                top -= 1
                if s.paren[u]:
                    smw_putc(b, c')')
                continue
            more = False
            for k in range(s.cursor[u], s.nbr_off[u + 1]):
                if s.he_kind[s.wnbr[k]] == SMW_HE_CHILD:
                    more = True
                    break
            v = edges[he].to
            if more:
                smw_putc(b, c'(')
            s.paren[v] = 1 if more else 0
            if sk is not NULL and u == sk.left and sk.remove_left:
                if sk.keep_bond_left:
                    smw_sticky_bond(b, &edges[he], s.he_dir[he])
            elif sk is not NULL and v == sk.right and sk.remove_right:
                if sk.keep_bond_right:
                    smw_sticky_bond(b, &edges[he], s.he_dir[he])
            else:
                smw_bond(b, &edges[he], smw_both_lowercase(atoms, ptr, edges, u, v, o), s.he_dir[he], o)
            if sk is not NULL and v == sk.right and sk.remove_right:
                pass
            else:
                smw_atom_env(ptr, edges, v, &order_sum, &wild, &narom)
                sign = smw_sign_of(structure, s, v, o)
                smw_atom(b, &atoms[v], order_sum, wild, narom, sign, o)
                smw_ring_tokens(structure, s, b, v, o)
            emitted[v] = 1
            s.cursor[v] = s.nbr_off[v]
            s.stack[top] = v
            top += 1
        idx += 1
    return 0


cdef void smw_ring_tokens(Structure structure, smw_scratch_t *s,
                          smw_buf_t *b, uint32_t u,
                          smw_opts_t *o) noexcept:
    """The ring bonds -- closures and attachments -- and the digits that follow atom `u`'s token."""
    cdef atom_t *atoms = structure.atoms()
    cdef uint32_t *ptr = csr_ptr(structure)
    cdef halfedge_t *edges = csr_edges(structure)
    cdef uint32_t k, he, rev
    for k in range(s.nbr_off[u], s.nbr_off[u + 1]):
        he = s.wnbr[k]
        if s.he_kind[he] == SMW_HE_ATTACH:
            # THE BOND TOKEN IS WRITTEN HERE, on a bond whose other half this string does not contain.
            # Both fragments write it, and that is measured safe: RDKit, Indigo, OpenBabel and chython 2
            # all accept `C=%10` + `C=%10` and all four reject a clash, so a token at one end only
            # would make the join order-dependent for no gain.  `o.asymmetric_closure` deliberately
            # does not reach here -- it is about which END of a closure carries the token, and a
            # fragment has only one end of this bond to speak for.
            #
            # DIRECTION 0: spec §13.4.  A cis/trans configuration across a cut is refused in
            # `smw_directions` and reported, so there is nothing to write and no chance of writing a
            # `/` whose partner is in another string.
            smw_bond(b, &edges[he], smw_both_lowercase(atoms, ptr, edges, u, edges[he].to, o), 0, o)
            smw_closure_number(b, s.he_close[he])
            continue
        if s.he_kind[he] != SMW_HE_CLOSURE:
            continue
        rev = <uint32_t> (csr_find_at(ptr, edges, edges[he].to, u) - edges)
        if o.asymmetric_closure and s.he_close[rev] and s.he_close[rev] == s.he_close[he]:
            # The opening side already carried the token; `he_close[rev]` is set only once the far
            # atom has been through here, so this is the closing side.
            if s.out_idx[edges[he].to] < s.out_idx[u]:
                smw_closure_number(b, s.he_close[he])
                continue
        # A directional token on a ring-closure bond is written at the OPENING side only.  Writing
        # it at both ends needs the two characters to be opposite, which readers disagree about;
        # one end is unambiguous everywhere.
        if s.he_dir[he] and s.out_idx[edges[he].to] < s.out_idx[u]:
            smw_bond(b, &edges[he], smw_both_lowercase(atoms, ptr, edges, u, edges[he].to, o), 0, o)
        else:
            smw_bond(b, &edges[he], smw_both_lowercase(atoms, ptr, edges, u, edges[he].to, o), s.he_dir[he], o)
        smw_closure_number(b, s.he_close[he])


# ------------------------------------------------------------------------------------------------
# STEREO.  THE ONLY PLACE IN THIS FILE THAT MAY READ A PARITY (ruling F26, note 2 at the top).
#
# WHAT `@` MEANS, and the three facts it rests on.  Two are external and measured, one is the arena's:
#
#   1. OpenSMILES: looking from the FIRST neighbour in the written order towards the centre, the
#      remaining three appear ANTICLOCKWISE for `@`.
#   2. Measured against RDKit 2026.03.4 on 2026-09-02: an implicit hydrogen occupies the position in
#      the written order WHERE IT IS WRITTEN -- immediately after the preceding atom, or first when
#      there is none.  `[C@H](F)(Cl)Br` and `F[C@@H](Cl)Br` are one molecule and
#      `[C@H](F)(Cl)Br` and `F[C@H](Cl)Br` are two, which is the transposition of the first two
#      positions and nothing else.  A writer that puts the hydrogen last unconditionally is right
#      for every atom with a parent and inverted for every component's first atom.
#   3. The core's parity, in the ruling-F26 refs frame: parity 2 (odd) is anticlockwise, hence `@`.
#      Measured by the MDL epic through V2, which is the only external anchor the value has -- the
#      core itself defines `even`/`odd` and nothing else, so this correspondence is a CONVENTION
#      SHARED WITH THE READER (`smi_*`) rather than something derivable here.  It is frame-relative:
#      "parity 2 means `@`" is only true of the refs IN THE ORDER F26 names them, which is exactly
#      why this function exists instead of a byte reaching the output.
#
# The 4-direction bookkeeping needs no bounds check beyond `n_refs == 4`: every neighbour is exactly
# one direction (a double bond included -- the second lobe is not a place a substituent can sit, and
# a triple bond refuses the atom outright), so `degree + implicit_h + lone_pair == 4` follows from
# `n_refs == 4` and the loops below cannot run off the end of a 4-entry list.
cdef void smw_direction_order(Structure structure, smw_scratch_t *s, uint32_t slot,
                              uint32_t *want) noexcept nogil:
    """`want[0:4]`: the atom's four directions IN THE ORDER THEY ARE ABOUT TO BE WRITTEN.

    The parent, then the implicit hydrogens, then the ring bonds and the branches in `wnbr`
    order, then the lone pair.  `wnbr` is already parent-ringbonds-children (see `smw_traverse`), so
    the only insertion is the hydrogen's, and the pad at the end is the lone pair's.

    AN ATTACHMENT (spec §13) IS ONE OF THESE DIRECTIONS and needs no special case, which is the whole
    reason a tetrahedral configuration survives a cut: `%12` occupies a written position exactly as a
    closure digit does, so a reader of the JOINED string counts the same four directions in the same
    order this loop did.  It holds because a join concatenates fragments as separate `.` components
    and never rewrites one, so a fragment's own text is read the way it was written.

    `SU_NO_REF` for a direction with no atom of its own.  Meaningful only for an atom that anchors a
    four-direction unit, which is the caller's business to establish; `smw_traversal`'s `directions`
    key exposes exactly this list for the atoms where it means something, so that a test can put the
    writer's own order through `MoleculeContainer.translate_stereo` and compare the answer to the
    sign in the string.  That comparison is the anti-drift test, and it only works because both
    sides read THIS function's output rather than each building an order.
    """
    cdef halfedge_t *edges = csr_edges(structure)
    cdef uint32_t i = s.nbr_off[slot]
    cdef uint32_t k = 0
    cdef uint32_t implicit
    if i < s.nbr_off[slot + 1] and s.he_kind[s.wnbr[i]] == SMW_HE_PARENT:
        want[0] = edges[s.wnbr[i]].to
        k = 1
        i += 1
    implicit = smw_written_h(structure.atoms() + slot)
    while implicit and k < 4:
        want[k] = SU_NO_REF
        k += 1
        implicit -= 1
    while i < s.nbr_off[slot + 1] and k < 4:
        want[k] = edges[s.wnbr[i]].to
        k += 1
        i += 1
    while k < 4:
        want[k] = SU_NO_REF
        k += 1


cdef uint32_t smw_sign_of(Structure structure, smw_scratch_t *s, uint32_t slot,
                          smw_opts_t *o) noexcept nogil:
    """`@` (1) or `@@` (2) for the atom in `slot`, or 0 for no sign.

    SU_CIS_TRANS and SU_ATROPISOMER return 0: the first one's configuration is `/` and `\\` on the
    two single bonds and lands in `he_dir`, the second one has no SMILES syntax at all.  SU_ALLENE
    is an atom sign like this one and goes to `smw_allene_sign_of`, which is where the axial case's
    two extra problems -- the four directions are not this atom's neighbours, and the string
    interleaves them -- are argued.

    The unnamed direction that is NOT a hydrogen -- the sulfur lone pair -- goes LAST and does not
    move when the atom leads its component.  Measured the same day: `[S@](=O)(C)CC` and
    `O=[S@](C)CC` are one molecule to RDKit, so unlike the hydrogen the lone pair has no positional
    rule, and V2 agrees (its frame is the three named substituents with the fourth direction fixed
    at the end).  Last rather than first is then the only remaining choice, and it is the one that
    makes the lone pair's position in the written order equal to its position in `refs`.
    """
    cdef atom_t *atoms
    cdef stereo_unit_t *u
    cdef uint32_t want[4]
    cdef uint32_t perm[4]
    cdef uint32_t norefs[4]     # the positions in `refs` that hold an unnamed direction
    cdef uint32_t i, j, nnoref = 0, nr = 0, used = 0
    cdef uint8_t parity

    if not o.stereo:
        return 0
    if s.frame_cut[slot]:
        # A sticky end took a neighbour away outright (`smw_sticky_frames`), so the string shows this
        # atom with one fewer direction than the parity describes.  Refused here rather than adjusted:
        # an atom that has lost a substituent is not the same stereocentre, and the writer does not get
        # to decide what the record would have said about the one that remains.
        return 0
    atoms = structure.atoms()
    if not structure_parity_at(structure, slot):
        return 0
    # A STATED configuration is written whether or not the unit is marked stereogenic, which is why
    # `smw_prepare` builds the table through `ensure_stereo_units_unmarked` (ruling F70).  Dropping
    # the sign of an unjustified parity would be a silent edit of the input, and the container
    # already has `stereo_rejections` for the caller who wants to know.
    u = stereo_unit_of(structure, slot)
    if u is NULL or u.n_refs != 4:
        return 0
    if u.kind == SU_ALLENE:
        # An axial sign is unreadable without the `=` tokens that show the axis: under `!b` the
        # string is `CCC`, and an `@` on it reads as a tetrahedral centre on a two-coordinate carbon.
        # Refusing here is what puts the unit in `lost` -- `smw_directions` asks this same function.
        if not o.bonds:
            return 0
        return smw_allene_sign_of(structure, s, slot, u)
    if u.kind != SU_TETRA:
        return 0
    if smw_h_frame_unknown(s, &atoms[slot], slot, 4):
        # The sign is over four POSITIONS in the written order and one of them would be a hydrogen
        # whose existence the arena does not claim.  Refused, and `smw_directions` reports it through
        # `lost` by asking this same predicate -- a bracket cannot say "unknown" either, so writing a
        # sign here would state a frame the string does not show.
        return 0
    perm[0] = 0; perm[1] = 0; perm[2] = 0; perm[3] = 0   # seatbelt, as in `translate_stereo`
    smw_direction_order(structure, s, slot, want)

    # `perm[i] = j` means want[i] is refs[j].  The unnamed positions are consumed in order, exactly
    # as `translate_stereo` consumes them, so the two cannot disagree about which `None` is which --
    # and the hydrogen is always before the lone pair in both lists, which is what makes that
    # in-order consumption the right correspondence rather than merely a determinate one.
    for j in range(4):
        if u.refs[j] == SU_NO_REF:
            norefs[nnoref] = j
            nnoref += 1
    for i in range(4):
        if want[i] == SU_NO_REF:
            perm[i] = norefs[nr]
            nr += 1
        else:
            for j in range(4):
                if u.refs[j] == want[i] and not (used & (1u << j)):
                    perm[i] = j
                    used |= 1u << j
                    break
    parity = structure_parity_at(structure, slot)
    return 1 if translate_parity(parity, perm) == 2 else 2


# ------------------------------------------------------------------------------------------------
# THE ALLENE'S SIGN.  `@` on the CENTRE atom, over the four directions of the two CHAIN TERMINALS.
#
# An axial unit is written like a tetrahedral one -- OpenSMILES calls it EXTENDED TETRAHEDRAL -- and
# the two differences are both about which four directions the sign is over:
#
#   * they are not the signed atom's neighbours.  The centre has two (the axis); the four are the
#     terminals' non-chain directions, and for a longer odd cumulene the terminals are several bonds
#     away.  So this function walks the chain instead of reading `wnbr` at `slot`.
#   * the string INTERLEAVES them, and it does not matter.  The far terminal's subtree is written
#     inside the near terminal's branch list, so the four direction tokens appear in one of exactly
#     three arrangements: `(n0 n1 f0 f1)`, `(n0 f0 f1 n1)` or `(f0 f1 n0 n1)`.  Every one of those is
#     an EVEN permutation of the others ([0,2,3,1] has two inversions; the pair exchange has four),
#     so the grouped tuple this function builds is the string's own order up to a sign-preserving
#     permutation.  That is also why the near/far CHOICE below is free: it is the pair exchange.
#     Only the order WITHIN each terminal can flip the sign, and that is the order `wnbr` gives.
#
# WHAT `@` MEANS HERE.  The same rule as `smw_sign_of`, applied to those four directions: looking from
# the first one towards the axis, the remaining three appear anticlockwise.  Measured on 2026-09-02
# against RDKit 2026.03.4 for the tetrahedral case, from a hand-built conformer rather than from a
# SMILES string so that the geometry is the input and not an inference: for CBrClFI with
# Br(1,1,1) F(1,-1,-1) Cl(-1,1,-1) I(-1,-1,1) the signed volume over `(Br, F, Cl, I)` is NEGATIVE and
# RDKit writes `F[C@](Cl)(Br)I`, which is `@` in the frame `(Br, F, Cl, I)` (its written order is an
# even permutation of it).  Mirroring z gives `+` and `F[C@@](Cl)(Br)I`.  So: NEGATIVE SIGNED VOLUME
# OVER THE WRITTEN ORDER == `@`, and the core's parity 2 (odd) is that same handedness.
#
# THE AXIAL SIGN HAS NO EXTERNAL ANCHOR IN THIS TREE, and `SMW_ALLENE_AT_FOR_ODD` below is the one
# place to correct it if a cross-format consumer ever says otherwise.  What was tried, 2026-09-02:
#
#   * RDKit 2026.03.4 DROPS every allene tag at sanitization -- `CC(F)=[C@]=C(F)C` parses to a
#     CHI_TETRAHEDRAL_CCW on atom 3 and `SanitizeMol` clears it, with or without
#     `SetUseLegacyStereoPerception(False)` or `SetAllowNontetrahedralChirality(True)`.  It also does
#     not perceive one from 3D coordinates.  So RDKit cannot state the axial case at all.
#   * OpenBabel 3.1.0 refuses it on read ("Ignoring stereochemistry.  Not enough connections").
#   * Indigo (InChI API 1.06) DOES read and write it, and its arithmetic agrees with the model above:
#     `CC(F)=[C@]=C(F)C` comes back as `CC(=[C@@]=C(C)F)F`, which is one within-pair swap and one
#     flipped tag.  But its InChI export drops the layer, so it cannot carry the sign to a second
#     authority, and its InChI READER drops it too.
#   * libinchi DOES perceive it from 3D: the two mirror conformers of BrC(F)=C=C(Cl)I give
#     `/t1-/m1/s1` and `/t1-/m0/s1`.  The InChI bridge cannot carry that into the arena, though --
#     `_inchi.pxi`'s SU_ALLENE and SU_CIS_TRANS branches pass the chain ATOMS to `translate_stereo`
#     in InChI's own `(X, A, B, Y)` neighbour order, and those are not the unit's refs, so both
#     branches raise `ValueError` before libinchi is reached.  Reported to that epic; when it is
#     fixed, `inchi_to_molecule` on those two strings is the measurement that pins this DEF.
#   * chython 2 reads and writes the tag, but its stored bool is in ITS frame, and nothing states
#     that frame's geometry -- asking V2 what its own convention is, is circular.
#
# Until then the sign is CHOSEN, not measured, and it is chosen to be the tetrahedral rule applied
# unchanged: one handedness convention for the whole file, so that a reader agreeing with `@` on a
# stereocentre agrees with `@` on an axis.  The round trip is what is actually pinned by tests -- two
# enantiomers are two strings, one configuration is one string over every creation order -- and those
# hold under either value of the DEF, which is exactly why the value needs saying out loud.
DEF SMW_ALLENE_AT_FOR_ODD = 1     # 1: parity 2 (odd) writes `@`.  0 writes `@@`.  THE FLIP POINT.


cdef inline void smw_terminal_pair_order(Structure structure, smw_scratch_t *s, uint32_t t,
                                         uint32_t chain, uint32_t *want) noexcept nogil:
    """Terminal `t`'s two NON-CHAIN directions, in the order their tokens are written.

    The same rule as `smw_direction_order` -- parent, implicit hydrogens, then closures and branches
    in `wnbr` order -- with the chain neighbour dropped wherever it sits.  It can sit anywhere: the
    axis is this atom's parent when the traversal came down the chain, a child when it came up one of
    the substituents, and a ring closure when the allene is in a macrocycle.

    `SU_NO_REF` pads, so a terminal carrying an implicit hydrogen reads `(heavy, None)` or
    `(None, heavy)` according to where the hydrogen is written -- which is inside the bracket, hence
    after the parent bond and before the closure digits.  That position is the measured one
    (`smw_sign_of`'s note 2) and it is the whole reason this cannot be "the heavy ones, ascending".
    """
    cdef halfedge_t *edges = csr_edges(structure)
    cdef uint32_t i = s.nbr_off[t]
    cdef uint32_t k = 0
    cdef uint32_t implicit
    cdef uint32_t to
    if i < s.nbr_off[t + 1] and s.he_kind[s.wnbr[i]] == SMW_HE_PARENT:
        to = edges[s.wnbr[i]].to
        if to != chain:
            want[k] = to
            k += 1
        i += 1
    implicit = smw_written_h(structure.atoms() + t)
    while implicit and k < 2:
        want[k] = SU_NO_REF
        k += 1
        implicit -= 1
    while i < s.nbr_off[t + 1] and k < 2:
        to = edges[s.wnbr[i]].to
        if to != chain:
            want[k] = to
            k += 1
        i += 1
    while k < 2:
        want[k] = SU_NO_REF
        k += 1


cdef bint smw_allene_order(Structure structure, smw_scratch_t *s, uint32_t slot,
                           stereo_unit_t *u, uint32_t *want, uint32_t *perm) noexcept nogil:
    """`want[0:4]` the axis's four directions grouped by terminal, `perm[i] = j` meaning want[i] is
    refs[j].  False when the record and the arena disagree, which is the caller's `lost`.

    `want[0:2]` is the terminal owning the STORED pair `refs[0:2]`, decided by membership rather than
    by re-deriving perception's "lower slot first" -- the answer has to be right for the record in
    hand, not for the record perception would build today.  A terminal with no named direction at all
    cannot vote and does not need to: with both of its slots unnamed the two orderings differ by the
    pair exchange, which is even.

    The unnamed slots are consumed WITHIN their pair, which is `translate_stereo`'s rule for a bond
    kind (a pinned slot is frozen at its offset within its pair, ruling F55) and not the global
    in-order consumption `smw_sign_of` uses for an atom kind.
    """
    cdef uint32_t *ptr = csr_ptr(structure)
    cdef halfedge_t *edges = csr_edges(structure)
    cdef uint32_t eo[4]             # end 0's two directions, then end 1's
    cdef uint32_t ends[2]
    cdef uint32_t chains[2]         # each end's own chain neighbour
    cdef uint32_t i, j, k, p, base, cur, prev, nxt, budget, v
    cdef uint32_t n_ends = 0
    cdef uint32_t first = 0
    cdef uint32_t used
    cdef bint decided = False
    cdef bint found

    # Both ways along the chain from the centre.  `budget` is the guard against a fully cumulated
    # ring, where the walk would never reach a terminal; perception cannot emit a unit for one, so
    # this is a forged-arena seatbelt and not a live case.
    for k in range(ptr[slot], ptr[slot + 1]):
        if not _is_chain_bond(&edges[k]):
            continue
        if n_ends == 2:
            return False            # three chain bonds at the centre: not an axis
        prev = slot
        cur = edges[k].to
        budget = s.n
        while budget:
            nxt = _chain_next(ptr, edges, cur, prev)
            if nxt == SU_NO_REF:
                break
            prev = cur
            cur = nxt
            budget -= 1
        if not budget:
            return False
        ends[n_ends] = cur
        chains[n_ends] = prev
        n_ends += 1
    if n_ends != 2 or ends[0] == ends[1]:
        return False
    # A THREE-ATOM AXIS ONLY, which is `chains[i] == slot` for both ends.  A longer odd cumulene is
    # axially chiral and the arena names it the same way, but SMILES HAS NO SYNTAX FOR IT: measured
    # 2026-09-02, Indigo (the one allene-capable reader in reach) refuses `FC=C=[C@@]=C=CF` outright
    # -- "chirality on atom 3 makes no sense" -- so writing the sign there would produce a string a
    # real reader rejects, which is worse than a string that says less.  chython 2 does write it, so a
    # longer cumulene's sign is dropped here, and the refusal reaches `lost`.
    if chains[0] != slot or chains[1] != slot:
        return False
    # The same refusal as `smw_sign_of`'s, one atom further out: a terminal whose implicit count is
    # unknown does not say whether its two directions are (heavy, H) or (heavy, nothing), and those are
    # different frames.  `fills` is 3 here -- the chain bond plus both directions -- so a terminal with
    # two heavy substituents is fine and only an under-substituted one refuses.
    for i in range(2):
        if smw_h_frame_unknown(s, structure.atoms() + ends[i], ends[i], 3):
            return False
        # A DROPPED TERMINAL means the axis itself was cut, and the sign is a claim about the relation
        # between the two ends' directions.  A dropped SUBSTITUENT at a retained terminal is NOT
        # refused: its attachment holds a written position at that terminal exactly as a closure digit
        # would, so the frame is still four positions in this string (spec §13.4, and the same argument
        # `smw_direction_order` makes for a tetrahedral centre).
        if s.dropped[ends[i]]:
            return False
    smw_terminal_pair_order(structure, s, ends[0], chains[0], &eo[0])
    smw_terminal_pair_order(structure, s, ends[1], chains[1], &eo[2])

    for i in range(2):
        for k in range(2):
            v = eo[i * 2 + k]
            if v == SU_NO_REF:
                continue
            if v == u.refs[0] or v == u.refs[1]:
                first = i
                decided = True
            elif v == u.refs[2] or v == u.refs[3]:
                first = 1 - i
                decided = True
            if decided:
                break
        if decided:
            break
    want[0] = eo[first * 2]
    want[1] = eo[first * 2 + 1]
    want[2] = eo[(1 - first) * 2]
    want[3] = eo[(1 - first) * 2 + 1]

    for p in range(2):
        base = 2 * p
        used = 0
        for k in range(2):
            found = False
            for j in range(2):
                if used & (1u << j):
                    continue
                if u.refs[base + j] == want[base + k]:
                    perm[base + k] = base + j
                    used |= 1u << j
                    found = True
                    break
            if not found:
                return False        # a direction the record does not have, or two for one slot
    return True


cdef uint32_t smw_allene_sign_of(Structure structure, smw_scratch_t *s, uint32_t slot,
                                 stereo_unit_t *u) noexcept nogil:
    """`@` (1) or `@@` (2) for the axis anchored at `slot`, or 0 when it cannot be written.

    Same shape as `smw_sign_of`'s tail, and deliberately the same constant: `translate_parity` is the
    arena's arithmetic for both kinds, and for a bond kind it reproduces `translate_stereo` exactly
    (a within-pair swap is one transposition, the pair exchange is two).
    """
    cdef uint32_t want[4]
    cdef uint32_t perm[4]
    cdef uint32_t odd_tag = 1 if SMW_ALLENE_AT_FOR_ODD else 2
    cdef uint8_t parity
    perm[0] = 0; perm[1] = 0; perm[2] = 0; perm[3] = 0   # seatbelt, as in `translate_stereo`
    if not smw_allene_order(structure, s, slot, u, want, perm):
        return 0
    parity = structure_parity_at(structure, slot)
    return odd_tag if translate_parity(parity, perm) == 2 else 3 - odd_tag


# ------------------------------------------------------------------------------------------------
# CIS/TRANS.  `/` and `\`, the one part of SMILES stereo that is not a property of a single atom.
#
# WHAT A DIRECTION MEANS.  Written `A/B` the bond rises left to right, so B is up relative to A and A
# is down relative to B.  `he_dir` therefore holds a token FOR ONE SIDE, and the two halves of a bond
# always hold opposite values.  That opposition is not bookkeeping: it IS the conjugation coupling,
# and it is why the shared single bond of a 1,3-diene cannot be given two independent tokens.
#
# THE CIS RULE, derived from the above rather than asserted: substituent `a` on terminal `n` and `b`
# on terminal `m` lie on the SAME side exactly when up(n->a) == up(m->b).  Check it against a known
# string.  `F/C=C/F` is trans.  Its first token says C is up from F, hence F is DOWN from C, so
# up(n->a) is down; its second says F is UP from C, so up(m->b) is up.  Different, and the molecule
# is trans.  The rule holds, and it is stated in terms of directions pointing AWAY from each
# terminal, which is the orientation the code uses everywhere below.
#
# WHICH PARITY IS CIS.  Parity 2 (odd) means refs[0] and refs[2] -- one named direction from each
# terminal, both always present by rulings F26 and F47 -- lie on the SAME side.  MEASURED 2026-09-02
# through chython 2, with RDKit 2026.03.4 confirming the geometry independently: `F/C=C\F` is Z, V2
# stores `bond.stereo` True for it, `_alkene_translate[(0, 1)]` is False so that bool is V2's answer
# for its own frame pair with no flip in between, and V2's True is core parity 2 -- the same mapping
# the tetrahedral sign already uses, taken uniformly so that task 12's bridge needs no per-kind flip.
#
# Like `@`, this correspondence HAS NO ANCHOR INSIDE THE CORE: the arena defines even and odd and
# nothing else.  It is a convention shared with the reader's `smi_cis_sign`, and inverting one
# inverts both.
#
# WHY A SOLVER AND NOT AN ASSIGNMENT.  A direction is a property of a BOND, but a configuration is a
# property of a double bond's two ends, so one bond can be constrained by two configurations at once
# and the constraints have to be solved together.  Three relations, all of them "same" or "opposite":
#
#   * a bond seen from its two ends: OPPOSITE, always;
#   * the two directions of one terminal: OPPOSITE, always (they are the two in-plane positions);
#   * refs[0] and refs[2] across the double bond: SAME for parity 2, OPPOSITE for parity 1.
#
# That is a 2-colouring, so it is solved by breadth-first search from one seeded bond per connected
# component -- and it can FAIL, on a cycle of constraints with odd total parity.  A cyclic polyene
# whose stated configurations cannot all hold at once is the real shape; small-ring alkenes never
# reach here because perception already refuses them (`_terminals_share_small_ring`).  On failure the
# whole component is unwound and every unit in it is reported through `lost`, because writing SOME of
# a contradictory set would produce a string that reads back as a molecule nobody stated.
cdef inline int smw_pair_slot(smw_scratch_t *s, Structure structure,
                              uint32_t t, uint32_t v) noexcept nogil:
    """Which of terminal `t`'s refs slots holds `v` -- 0/1 at the near end, 2/3 at the far one.

    -1 when `t` is not a live cis/trans terminal, or `v` is not one of its two directions.  Live
    means `smw_directions` phase 1 admitted it: parity known, partner found, every named direction on
    a single bond.
    """
    cdef stereo_unit_t *u
    cdef uint32_t base
    if s.unit_at[t] == SMW_NONE:
        return -1
    u = stereo_unit_of(structure, s.unit_at[t])
    base = 0 if s.unit_at[t] == t else 2
    if u.refs[base] == v:
        return <int> base
    if u.refs[base + 1] == v:
        return <int> (base + 1)
    return -1


cdef inline int smw_dir_offer(smw_scratch_t *s, uint32_t frm, uint32_t he, uint8_t val,
                              uint32_t *tail) noexcept nogil:
    """Assign `val` to half-edge `he`, enqueue it, and report whether that contradicts what is there.

    1 when `he` already holds the OTHER token: the stated configurations are unsatisfiable and the
    caller unwinds the component.  0 when the slot was empty (now assigned and queued) or already
    held this very token (nothing to do, and not a contradiction -- a 2-colouring reaches most nodes
    by several routes and they agree).
    """
    if s.he_dir[he]:
        return 1 if s.he_dir[he] != val else 0
    s.he_dir[he] = val
    s.dq[tail[0] * 2] = frm
    s.dq[tail[0] * 2 + 1] = he
    tail[0] += 1
    return 0


cdef void smw_directions(Structure structure, smw_scratch_t *s, smw_opts_t *o) noexcept nogil:
    """Fill `he_dir` with the `/` and `\\` tokens, and `lost` with the units that got none.

    Runs after the traversal, because the seed choice is made in emission order and because a
    ring-closure bond's token belongs to the opening side, which is a traversal fact.
    """
    cdef uint32_t *ptr = csr_ptr(structure)
    cdef halfedge_t *edges = csr_edges(structure)
    cdef stereo_unit_t *u
    cdef halfedge_t *e
    cdef uint32_t i, k, idx, he, frm, to, far, anchor, qfrm, qhe
    cdef uint32_t head = 0, tail = 0, cbase = 0
    cdef int bad

    if not o.stereo:
        return
    if not o.bonds:
        # NO BOND TOKENS AT ALL, so no configuration that needs one can be written: `/` has nowhere to
        # go and an axial `@` would read as a tetrahedral centre on a two-coordinate atom.  Reporting
        # them instead of returning early is what makes `!b` honest -- the TETRAHEDRAL signs still
        # come out, so the string carries some of the stereo, and a caller who cannot see which part
        # went would have to diff the molecule against the string to find out.
        for i in range(s.n):
            if s.dropped[i] or not structure_parity_at(structure, i):
                continue
            u = stereo_unit_of(structure, i)
            if u is not NULL and u.kind == SU_TETRA and smw_sign_of(structure, s, i, o) != 0:
                continue                   # an atom sign needs no bond token, so `!b` keeps it
            s.lost[i] = 1
            s.nlost += 1
        return
    for i in range(s.n):
        s.unit_at[i] = SMW_NONE
        s.partner_at[i] = SMW_NONE

    # PHASE 1 -- which units are writable at all, and where their terminals are.  A unit that fails
    # any test here is `lost` rather than silently skipped: a caller who asked for stereo and got a
    # string with no `/` in it deserves to be able to find out why.
    for i in range(s.n):
        if s.dropped[i]:
            # NOT REPORTED, unlike every other skip in this loop: a dropped atom's configuration is not
            # something this fragment failed to write, it is an atom the caller asked not to be here.
            # `lost` is about the string's coverage of what it claims to describe.
            continue
        if not structure_parity_at(structure, i):
            continue
        u = stereo_unit_of(structure, i)
        if u is NULL:
            # A STATED PARITY WITH NO UNIT UNDER IT.  Perception declines to name a frame it cannot
            # read, and an implicit count of H_UNKNOWN in the neighbourhood is one such refusal --
            # since arena 6835438, only where the count COULD have moved a written position (CHFClBr,
            # an allene with one unstated terminal), where before it was every such atom including
            # CFClBrI.  The narrowing changed which atoms arrive here and nothing else: the parity is
            # real, nothing will write it, so it is a loss and it is NAMED -- the writer does not get
            # to be quiet about a fact of the input because another mechanism was.
            s.lost[i] = 1
            s.nlost += 1
            continue
        if u.kind == SU_TETRA:
            # An atom sign, written during emission.  Asked here anyway, and for the same reason as
            # the axial branch below: `smw_sign_of` can refuse, and only its own answer can say
            # whether it did.  Re-deriving the conditions would be a second truth about what the
            # string carries, which is the shape of the defect note 4 exists for.
            if smw_sign_of(structure, s, i, o) == 0:
                s.lost[i] = 1
                s.nlost += 1
            continue
        if u.kind == SU_ALLENE:
            # An atom sign, like SU_TETRA -- but unlike SU_TETRA it can REFUSE (a chain that does not
            # end, a record whose refs the arena no longer has), and a refusal has to reach `lost`.
            # Asking the emitter's own function rather than re-deriving the conditions is the point:
            # the two answers cannot disagree about whether a sign was written.
            if smw_allene_sign_of(structure, s, i, u) == 0:
                s.lost[i] = 1
                s.nlost += 1
            continue
        if u.kind != SU_CIS_TRANS:
            # SU_ATROPISOMER, and it will never move out of here, because SMILES HAS NO SYNTAX FOR
            # IT.  Reporting it every time is the point -- a format that cannot carry a configuration
            # should say so rather than hand back a string that looks complete.
            s.lost[i] = 1
            s.nlost += 1
            continue
        far = stereo_unit_partner(structure, u)
        if far == SU_NO_REF or s.dropped[far]:
            # A CIS/TRANS UNIT DOES NOT SURVIVE A CUT (spec §13.4), and the dropped partner is one of
            # the two ways it can be split.  Unlike a tetrahedral sign, which is four positions at ONE
            # atom and so is entirely inside whichever fragment holds that atom, `/` and `\` are two
            # tokens on two different bonds whose meaning is their relation -- and after a join the
            # reader has no way to relate a token in one string to a token in another.  Refused and
            # REPORTED, because this half is a real loss of this fragment's own coverage.
            s.lost[i] = 1
            s.nlost += 1
            continue
        # Every named direction must sit on a SINGLE bond, because `/` replaces the single-bond
        # token and there is nowhere to write it otherwise.  A cumulene terminal's non-chain bonds
        # are single by construction (a second double bond would make it a chain atom, not a
        # terminal), so this is a guard against a forged arena rather than a live case -- and it is
        # cheaper than reasoning about one at the point where the token would be emitted.
        bad = 0
        for k in range(4):
            if u.refs[k] == SU_NO_REF:
                continue
            if s.dropped[u.refs[k]]:
                # The second way a cut splits this unit: the terminals are both here but a REFERENCE
                # atom is not, so one of the two tokens would have to go on a bond this string does not
                # contain.  Refused for the whole unit even when the other reference at that terminal
                # is retained: the frame is refs[0] against refs[2] and re-seating it onto a different
                # pair is a change of frame, not a smaller version of the same claim.
                bad = 1
                break
            e = csr_find_at(ptr, edges, i if k < 2 else far, u.refs[k])
            if e is NULL or e.order != 1:
                bad = 1
                break
        if bad:
            s.lost[i] = 1
            s.nlost += 1
            continue
        s.unit_at[i] = i
        s.unit_at[far] = i
        s.partner_at[i] = far
        s.partner_at[far] = i

    # PHASE 2 -- seed one bond per constraint component and propagate.  The seed is the first
    # directional half-edge in (emission index of the writing atom, written order at that atom),
    # which is a function of the canonical positions and therefore the same for every creation
    # order.  It is given `/`, so the sign of the whole component is decided by that one choice --
    # `F/C=C/F` and `F\C=C\F` are the same molecule and only one of them is the canonical string.
    for idx in range(s.nseq):
        frm = s.seq[idx]
        for k in range(s.nbr_off[frm], s.nbr_off[frm + 1]):
            he = s.wnbr[k]
            to = edges[he].to
            if s.he_kind[he] == SMW_HE_PARENT:
                continue                   # written from the other side; that side will seed it
            if s.he_kind[he] == SMW_HE_CLOSURE and s.out_idx[to] < s.out_idx[frm]:
                continue                   # the closing side of a closure carries no token
            if s.he_dir[he]:
                continue                   # an earlier component already reached this bond
            if smw_pair_slot(s, structure, frm, to) < 0 and \
                    smw_pair_slot(s, structure, to, frm) < 0:
                continue                   # not a direction of any live unit
            cbase = tail
            head = tail
            bad = smw_dir_offer(s, frm, he, 1, &tail)
            while head < tail and not bad:
                qfrm = s.dq[head * 2]
                qhe = s.dq[head * 2 + 1]
                head += 1
                bad = smw_dir_propagate(structure, s, ptr, edges, qfrm, qhe, &tail)
            if bad:
                # UNWIND THE WHOLE COMPONENT.  A contradiction cannot involve a half-edge assigned by
                # an EARLIER component: the constraint relation is symmetric and each component is
                # explored to closure, so any half-edge reachable from this seed was reachable from
                # that one too and would already be assigned -- and then this seed would have been
                # skipped above.  So everything to clear is in [cbase, tail).
                for i in range(cbase, tail):
                    s.he_dir[s.dq[i * 2 + 1]] = 0
                for i in range(cbase, tail):
                    anchor = s.unit_at[s.dq[i * 2]]
                    if anchor != SMW_NONE and not s.lost[anchor]:
                        s.lost[anchor] = 1
                        s.nlost += 1
                tail = cbase


cdef inline int smw_dir_propagate(Structure structure, smw_scratch_t *s,
                                  uint32_t *ptr, halfedge_t *edges,
                                  uint32_t frm, uint32_t he, uint32_t *tail) noexcept nogil:
    """Offer the three relations of one assigned half-edge to its neighbours.  1 on contradiction."""
    cdef stereo_unit_t *u
    cdef uint32_t to = edges[he].to
    cdef uint8_t val = s.he_dir[he]
    cdef uint32_t rev, other, oref, mate, base, anchor
    cdef int j

    # The same bond from the other end, always the opposite token.
    rev = <uint32_t> (csr_find_at(ptr, edges, to, frm) - edges)
    if smw_dir_offer(s, to, rev, 3 - val, tail):
        return 1
    j = smw_pair_slot(s, structure, frm, to)
    if j < 0:
        return 0                           # `frm` is not a terminal; `to` is handled through `rev`
    anchor = s.unit_at[frm]
    u = stereo_unit_of(structure, anchor)
    base = 0 if j < 2 else 2
    # The terminal's OTHER direction, opposite to this one.  Reached from either slot, which is what
    # lets the frame constraint below be stated for slot 0 alone.
    other = u.refs[base + 1] if <uint32_t> j == base else u.refs[base]
    if other != SU_NO_REF:
        if smw_dir_offer(s, frm,
                         <uint32_t> (csr_find_at(ptr, edges, frm, other) - edges),
                         3 - val, tail):
            return 1
    if <uint32_t> j != base:
        return 0                           # the frame is refs[0] against refs[2] and nothing else
    mate = s.partner_at[frm]
    oref = u.refs[base ^ 2]                # slot 0 of the other pair: 2 from 0, 0 from 2
    # parity 2 (odd, configured) is CIS, so the two frame directions get the SAME token.
    return smw_dir_offer(s, mate,
                         <uint32_t> (csr_find_at(ptr, edges, mate, oref) - edges),
                         val if structure_parity_at(structure, anchor) == 2 else 3 - val, tail)


# ------------------------------------------------------------------------------------------------
# CXSMILES.
cdef tuple smw_tail_parts(Structure structure, smw_scratch_t *s, dict groups, dict aliases):
    """The tail as STRUCTURE: `(radicals, abs_atoms, and_groups, or_groups, labels)`, indices into
    the emitted order.

    Separate from the text because a detached fragment's tail has to be RE-INDEXED when fragments are
    joined -- every index shifts by the atoms written before it, and two fragments' `&1` are not the
    same group.  Doing that by editing the string would mean parsing back out of it, and a formatter
    whose output is its own input is how the two sides drift apart.

    Indices are positions in the emitted order, which is what CXSMILES means by an atom index.
    `groups` is {n: (kind, number)} -- `write_smiles` INVERTS `canonical_stereo_groups()`,
    whose own shape is {(kind, canonical_id): [n, ...]}, so that this loop can ask about one
    atom at a time -- or None to read the stored bytes.  The canonical view is used for canonical
    output so the group NUMBERS are canonical too, and the stored bytes for stored-order output,
    whose point is to show what is stored.  An atom absent from `groups` carries no group.

    `labels` is UNLIKE the other four: a list POSITIONAL in the emitted order, one entry per atom,
    holding the alias bytes or None -- because that is what `$...$` is, and because a join then
    concatenates the lists instead of shifting indices.  Empty when no atom carries an alias, which
    is the common case and the one that must not pay for the field.  `aliases` is
    `MoleculeContainer.aliases`, keyed by stable id.
    """
    cdef atom_t *atoms = structure.atoms()
    cdef uint8_t *sg = structure_stereo_groups(structure)
    cdef uint32_t idx, u
    cdef list radicals = []
    cdef list abs_atoms = []
    cdef dict and_groups = {}
    cdef dict or_groups = {}
    cdef list labels = [None] * s.nseq if aliases else []
    cdef uint32_t kind, number
    cdef object entry

    for idx in range(s.nseq):
        u = s.seq[idx]
        if at_radical(&atoms[u]):
            radicals.append(idx)
        if aliases:
            labels[idx] = aliases.get(atoms[u].n)
        if groups is None:
            kind = sg_kind(sg[u])
            number = sg_group(sg[u])
        else:
            entry = groups.get(atoms[u].n)
            if entry is None:
                continue
            kind = <uint32_t> entry[0]
            number = <uint32_t> entry[1]
        if kind == 1:
            abs_atoms.append(idx)
        elif kind == 2:
            if number not in or_groups:
                or_groups[number] = []
            or_groups[number].append(idx)
        elif kind == 3:
            if number not in and_groups:
                and_groups[number] = []
            and_groups[number].append(idx)
    return (radicals, abs_atoms, and_groups, or_groups, labels)


cdef str smw_label_text(bytes raw):
    """One atom label as `$...$` can hold it: `&#NN;` for every character it cannot.

    `;` ends an entry, `$` ends the field, `|` ends the block and `&` starts a reference, so those
    four have no literal spelling; a space is escaped because a SMILES file's title column begins at
    one, and everything outside printable ASCII because a SMILES is ASCII.  Decimal references, and
    the code point rather than the byte: Marvin 25.1.3 writes `αβ` as `&#945;&#946;` and reads
    `&#32;` back as a space, so the wire form is measured and not chosen.

    Alias bytes that are not UTF-8 are read as Latin-1, which cannot fail -- the alternative is
    refusing to write a molecule over the encoding of a display label.
    """
    cdef str text
    cdef list out = []
    cdef Py_UCS4 c
    try:
        text = raw.decode('utf8')
    except UnicodeDecodeError:
        text = raw.decode('latin1')
    for c in text:
        if c in ';$|& ' or c < ' ' or c > '~':
            out.append('&#%d;' % ord(c))
        else:
            out.append(c)
    return ''.join(out)


cdef str smw_tail_text(tuple t, list fgroups=None):
    """`(radicals, abs_atoms, and_groups, or_groups, labels)` as ` |...|`, or `''` when there is
    nothing.

    The ONE formatter, used by `write_smiles`, by `DetachedSmiles.join` and by
    `write_reaction_smiles`: a join produces the same tail a direct write of the joined molecule
    would, and the only way to be sure of that is for the characters to come from one place.

    `fgroups` is the reaction writer's, and only the reaction writer's: a list of component-index
    lists for the `f:` field, which says "these components are one molecule".  A molecule has no use
    for it -- reading its string back gives one container whatever its components -- but a reaction
    does, and without it a salt reactant comes back as two reactants.  It goes LAST, after `^1:`.
    """
    cdef list radicals = <list> t[0]
    cdef list abs_atoms = <list> t[1]
    cdef dict and_groups = <dict> t[2]
    cdef dict or_groups = <dict> t[3]
    cdef list labels = <list> t[4]
    cdef list parts = []
    cdef list group
    cdef object key

    if labels and any(labels):
        # ONE ENTRY PER ATOM, trailing empties included, which is what Marvin writes -- a three-atom
        # molecule labelled on its first atom is `|$Me;;$|`.  The field goes FIRST, also as Marvin
        # writes it, so a reader that stops looking at an unknown field still sees the labels.
        group = []
        for key in labels:
            group.append('' if key is None else smw_label_text(<bytes> key))
        parts.append('$%s$' % ';'.join(group))
    if radicals:
        parts.append('^1:' + ','.join(map(str, sorted(radicals))))
    if abs_atoms and (and_groups or or_groups):
        # ONLY BESIDE AN AND OR OR COLLECTION, where it says which centres are NOT in one.  Alone it
        # says only "these centres are absolute", which is what a configured atom outside any
        # collection already means, so the field would separate two spellings of one structure -- a
        # molfile that named STEABS and a SMILES that had nowhere to name it.  The cost is that a
        # round trip through SMILES turns a lone explicit ABS kind into an unspecified one.
        parts.append('a:' + ','.join(map(str, sorted(abs_atoms))))
    for key in sorted(and_groups):
        parts.append('&%d:%s' % (key, ','.join(map(str, sorted(and_groups[key])))))
    for key in sorted(or_groups):
        parts.append('o%d:%s' % (key, ','.join(map(str, sorted(or_groups[key])))))
    if fgroups:
        group = []
        for key in fgroups:
            group.append('.'.join(map(str, <list> key)))
        parts.append('f:' + ','.join(group))
    if not parts:
        return ''
    return ' |%s|' % ','.join(parts)


cdef int smw_cxsmiles(Structure structure, smw_scratch_t *s, smw_buf_t *b,
                      dict groups, dict aliases) except -1:
    """The ` |...|` tail into the buffer: `$...$` labels, `^1:` radicals and the three
    enhanced-stereo group forms."""
    # BOUND TO A LOCAL FIRST.  `<char *> (...).encode()` casts a temporary whose last reference is
    # the cast itself, so the pointer is dangling before `smw_puts` reads it -- a use-after-free
    # that works in practice most of the time, which is the worst kind.
    cdef bytes tail = smw_tail_text(smw_tail_parts(structure, s, groups, aliases)).encode('ascii')
    if not tail:
        return 0
    smw_puts(b, <const char *> tail, <size_t> len(tail))
    return 0


# ------------------------------------------------------------------------------------------------
# THE CANONICAL ORDER, AND WHY IT HAS TO BE TOLD ABOUT STEREO.
#
# `canonical_order()` is a function of the CONSTITUTION.  On a molecule whose constitution is
# symmetric it therefore leaves the atoms an automorphism can exchange sharing a pair of positions,
# with nothing constitutional to say which takes which -- and the extremal search then settles it
# from the slot order, which is the creation order.  Note 3 at the top of this file says every
# decision here is a function of `pos`; that is true, and it is not enough on its own, because `pos`
# itself is only canonical up to that group.
#
# Where the exchanged atoms differ in STEREO, the string differs.  Measured: (2Z,4E)-hexa-2,4-diene
# writes as TWO strings over its 720 creation orders, one per end the order happens to start from.
# The constitution is symmetric end to end, the CONFIGURATION is not, and RDKit reads both strings
# as one molecule -- so it is a spelling defect and not a lost configuration, which is exactly what
# makes it dangerous: nothing downstream can see it except by comparing two strings that should have
# been equal.  It is no spelling nicety, because `__eq__`/`__hash__` rest on the canonical string:
# `MoleculeContainer.signature` carries ONE stereo bit for a whole molecule (cis- and trans-2-butene
# differ in it, but all sixteen anchor-by-parity combinations of this diene take just two values),
# so it can screen and cannot decide, and an oscillating canonical string would put one compound in
# a set twice.
#
# THE FIX IS A SEED, NOT A TIE-BREAK.  `compute_atoms_order` takes per-atom starting labels and
# refines them, and `mol_canonical_order` forwards them, so a colouring that already separates the
# two ends leaves the extremal search no tie to break.  What may go into that colouring is settled
# by the stereo epic's standing requirement (ruling F95, quoted in `_stereo.pxi` above
# `_frame_free_parity_code`): A SEED TERM MUST BE sigma-EQUIVARIANT, NOT MERELY ENCODING-INVARIANT.
# A stored parity byte is perfectly stable for one atom order and meaningless across orders -- ruling
# F26 again, from the other direction -- so it is `_frame_free_parity_code` that goes in: the parity
# RE-BASED onto the frame the current colouring names, which is a fact about the molecule.
#
# So this is not new machinery.  `canonical_stereo_group_ids` already runs this fixpoint to number
# stereo groups invariantly, and both call the same helper; what is dropped here is that function's
# group-membership label, which is about ids the writer does not spell.  Sharing the helper is the
# point: if the parity term were reimplemented here the two would drift, and the drift would look
# like a canonical string that disagrees with a canonical group id on one symmetric molecule.
cdef bint smw_stereo_seed(Structure structure, uint32_t *seed_out) except -1:
    """Per-slot seed labels for `mol_canonical_order`: the refinement class plus the parity in it.

    Returns False when the molecule carries no configured parity at all, and then writes nothing --
    the caller passes NULL instead, so a stereo-free molecule takes exactly the path it took before
    this function existed and cannot have its output moved by it.  `seed_out` is caller-owned and
    holds atom_count words.

    The unit table must already be built; `smw_prepare` does it under the same `o.stereo` that gates
    this call, and the UNMARKED build is the right one for the same reason it is there (ruling F70).

    THE FIXPOINT, and why one round is not enough: round 0 colours by the stereo-blind classes, and a
    parity can only be re-based onto a frame those classes can NAME -- two directions sharing a
    colour give `_frame_free_parity_code` nothing to order them by, and it answers "configured but
    unnamed" (code 1).  Feeding the round's own colouring back in names more frames each round, so
    the loop runs until the class count stops rising.  Termination is the class count: the class is
    the leading digit of the seed, so each round REFINES the last one's partition and can never merge
    two classes, `compute_atoms_order` refines its seed, and the count is bounded by n -- so round k
    is reachable only if the count rose k - 1 times, and the loop cannot reach round n.  The cap
    below makes that a branch rather than a promise, and the state it would leave is the state the
    top-of-loop break leaves.

    WHAT IS DELIBERATELY NOT IN THE SEED.  The stored group id and the stored parity byte, for
    ruling F95's reasons.  The group KIND, which `canonical_stereo_group_ids` does fold in: two atoms
    that differ only by ABS-versus-AND write the same SMILES, and the CXSMILES tail is written from
    the canonical group view rather than from `pos`, so seeding on a kind would split a tie no token
    depends on.  A configuration this writer cannot spell -- an atropisomer, or an axis it refuses --
    IS in the seed, because it is a fact about the molecule and dropping it would make
    the order of an atropisomer's two halves depend on the creation order again; the string says less
    than the seed knows, which is the safe direction.

    Arithmetic: `cls * 4 + code`, `code` in 0..3, so the largest label is 4n + 3 and this is exact in
    uint32 for any molecule under 2**30 atoms -- one quarter of the arena's own atom ceiling.
    """
    cdef uint32_t n = structure.header.atom_count
    cdef stereo_unit_t *units = structure_stereo_units(structure)
    cdef uint32_t nunits = structure_stereo_unit_count(structure)
    cdef uint32_t i
    cdef bint configured = False
    for i in range(nunits):
        if units[i].anchor < n and structure_parity_at(structure, units[i].anchor):
            configured = True
            break
    if not configured:
        return False

    # This round's classes, this round's parity code per atom, and each unit's partner terminal --
    # one allocation, and `nunits <= n` by the anchor invariant the unit table carries.
    cdef uint32_t *block = <uint32_t *> PyMem_Malloc(<size_t> 3 * <size_t> n * sizeof(uint32_t))
    if block is NULL:
        raise MemoryError()
    cdef uint32_t *cur = block
    cdef uint32_t *par = cur + n
    cdef uint32_t *partner = par + n
    cdef Py_ssize_t classes, prev
    cdef Py_ssize_t rounds = 0
    try:
        for i in range(nunits):
            partner[i] = stereo_unit_partner(structure, &units[i])
        with nogil:
            prev = compute_atoms_order(structure, cur, NULL)
        if prev < 0:
            raise MemoryError('atom order refinement failed to allocate')
        _frame_free_parity_seed(structure, units, nunits, partner, cur, par, n)
        for i in range(n):
            seed_out[i] = cur[i] * 4 + par[i]
        while True:
            with nogil:
                classes = compute_atoms_order(structure, cur, seed_out)
            if classes < 0:
                raise MemoryError('atom order refinement failed to allocate')
            if classes == prev:
                break                       # the seed already in `seed_out` is the fixpoint
            prev = classes
            # Re-borrowed per round rather than held across the call (ruling F60): nothing above
            # appends a segment today, and the rule is about what a reader may assume.
            units = structure_stereo_units(structure)
            _frame_free_parity_seed(structure, units, nunits, partner, cur, par, n)
            for i in range(n):
                seed_out[i] = cur[i] * 4 + par[i]
            rounds += 1
            if rounds >= <Py_ssize_t> n:
                break
    finally:
        PyMem_Free(block)
    return True


cdef int smw_canonical_positions(Structure structure, uint32_t *pos, bint stereo) except -1:
    """Fill `pos[slot]` with the atom's canonical position, seeded with stereo when `stereo`.

    Raises `AutomorphismBudgetExceeded` on a truncated search rather than returning SOME labelling:
    a canonical string built on a truncated order is not canonical and reads exactly like one that
    is.

    The seed is taken ONLY under `o.stereo`, so `!s` output stays a function of the constitution
    alone -- two molecules that differ only in configuration then write one string, which is what
    makes `format(mol, '!s')` usable as a constitution key.  Seeding it anyway would cost nothing in
    invariance and would quietly break that.

    `stereo` is forwarded to `mol_canonical_order` as well, and for the same reason rather than as a
    convenience: it also gates the search's own parity fold, which would reach `!s` even with no seed
    passed here.  Measured on the 393-record corpus of `test/`, folding it moved 63 `!s` strings and
    every one of the 63 then failed to survive `!s` -> read -> `!s`.  So the cost of the contract is
    that `!s` pays for the unfolded search; `_canon_order`'s docstring carries the division.
    """
    cdef uint32_t n = structure.header.atom_count
    cdef uint32_t *seed = NULL
    cdef uint32_t flags = 0
    if n == 0:
        return 0
    if stereo:
        seed = <uint32_t *> PyMem_Malloc(<size_t> n * sizeof(uint32_t))
        if seed is NULL:
            raise MemoryError()
    try:
        if seed is not NULL and not smw_stereo_seed(structure, seed):
            PyMem_Free(seed)
            seed = NULL                     # no configured parity: the unseeded path, exactly
        mol_canonical_order(structure, seed, pos, &flags, stereo)
    finally:
        PyMem_Free(seed)
    return 0


def smw_stereo_seed_labels(MoleculeContainer molecule not None):
    """The stereo seed as `{n: label}`, or None when the molecule carries no configured parity.

    Exposed for the test suite, and for one specific kind of test: the seed's whole job is to SPLIT a
    refinement class that the constitution leaves tied, and the only other way to see that happen is
    to compare two strings and infer it.  A test that reads the labels can assert both directions --
    that an asymmetrically configured diene's two ends get different labels, and that a symmetrically
    configured one's do NOT, because there the automorphism is real and breaking it would be the bug.

    Only the equality classes of the labels mean anything (`compute_atoms_order` documents that), so a
    test may compare two labels and must not read a label's value.
    """
    molecule._require_clean()
    cdef Structure structure = molecule._structure
    cdef uint32_t n = structure.header.atom_count
    cdef list numbers = molecule._numbers
    cdef uint32_t *seed
    cdef uint32_t i
    cdef dict out = None            # assigned at declaration: the `return None` inside the `try`
    if n == 0:                      # leaves Cython unable to see that `return out` is unreachable
        return None
    ensure_stereo_units_unmarked(structure)
    seed = <uint32_t *> PyMem_Malloc(<size_t> n * sizeof(uint32_t))
    if seed is NULL:
        raise MemoryError()
    try:
        if not smw_stereo_seed(structure, seed):
            return None
        out = {}
        for i in range(n):
            out[numbers[i]] = seed[i]
    finally:
        PyMem_Free(seed)
    return out


# ------------------------------------------------------------------------------------------------
# THE PYTHON SURFACE.
cdef int smw_parse_spec(str spec, smw_opts_t *o) except -1:
    """The format spec, one key per writer option.

    `''` canonical with everything on; `!s` no stereo; `A` aromatic bonds rather than lowercase
    atoms; `m` atom mapping; `h` every hydrogen count explicit; `!b` no bond tokens; `!x` no
    CXSMILES; `!z` no charges; `a` asymmetric ring-closure bonds; `i` stored slot order;
    `r` a random atom order.

    `i` AND `r` TOGETHER RAISE, rather than one winning: they are two answers to the one question of
    where the atom order comes from, and letting the later key win would make `ir` and `ri` different
    specs -- which `normalize_smiles_spec` promises they are not.
    """
    o.canonical = True
    o.random_order = False
    o.stereo = True
    o.aromatic_bond = False
    o.mapping = False
    o.hydrogens = False
    o.bonds = True
    o.charges = True
    o.cxsmiles = True
    o.asymmetric_closure = False
    cdef Py_ssize_t i = 0
    cdef Py_ssize_t k = len(spec)
    cdef bint negate
    cdef str c
    while i < k:
        negate = False
        if spec[i] == '!':
            negate = True
            i += 1
            if i == k:
                raise ValueError("format spec ends with a bare '!'")
        c = spec[i]
        i += 1
        if c == 's':
            o.stereo = not negate
        elif c == 'A':
            o.aromatic_bond = not negate
        elif c == 'm':
            o.mapping = not negate
        elif c == 'h':
            o.hydrogens = not negate
        elif c == 'b':
            o.bonds = not negate
        elif c == 'x':
            o.cxsmiles = not negate
        elif c == 'z':
            o.charges = not negate
        elif c == 'a':
            o.asymmetric_closure = not negate
        elif c == 'i':
            o.canonical = negate
        elif c == 'r':
            o.random_order = not negate
        else:
            raise ValueError('unknown format key %r' % c)
    # after the loop, so that the refusal does not depend on which key was written first
    if o.random_order and not o.canonical:
        raise ValueError("format keys 'i' and 'r' both name the atom order and cannot be combined")
    return 0


def normalize_smiles_spec(str spec=''):
    """One format spec, spelled canonically -- and usable directly as a cache key.

    `normalize_smiles_spec('sm') == normalize_smiles_spec('ms') == 'm'`: the keys are resolved into the
    option struct `write_smiles` actually uses and then re-spelled in a fixed order with the defaults
    left out, so two specs normalize equal EXACTLY WHEN they select the same writer behaviour.  That is
    the property a cache needs and the reason this function exists here rather than being reimplemented
    by the caller: a second copy of the grammar would be a second truth, and the failure mode is a
    cache that returns the string for `!s` when asked for `s`.

    An unknown or refused key raises, with the same message `format()` gives -- so a caller may
    normalize first and know the write will not fail on the spec.

    THE BYPASS LIST -- entry points whose output is NOT a function of molecule state and spec alone, so
    a normalized spec does not identify their result and none of them may be cached under one:

    * `write_smiles(mol, 'r')` -- a fresh random atom order per call, drawn from `random`.
    * `detached_smiles(mol, cuts, spec, reserve)` -- the cuts and the reserved ids are the caller's.
    * `sticky_smiles(mol, left, right, spec, ...)` -- the two named atoms force the atom order, the
      four `remove_*`/`keep_bond_*` flags decide which tokens appear, and it drops the CXSMILES tail.

    **`r` IS THE ONE OF THE THREE REACHABLE THROUGH A SPEC**, so a spec-keyed cache must refuse `'r'`
    rather than store the first string it happens to see; under every other spec the writer is a pure
    function of molecule state and a normalized spec is a complete cache key.  The other two are
    separate functions rather than spec keys for that reason, and `format()` deliberately has no letter
    for either.  Whoever adds a fourth writes it here in the same commit.
    """
    cdef smw_opts_t o
    smw_parse_spec(spec, &o)
    cdef list out = []
    # One entry per key, in the order `smw_parse_spec` documents them, each with the default the
    # struct is initialised to.  A key appears only when it differs from that default, so the empty
    # spec normalizes to the empty string and stays the cheapest key there is.
    if not o.canonical:
        out.append('i')
    if o.random_order:
        out.append('r')
    if not o.stereo:
        out.append('!s')
    if o.aromatic_bond:
        out.append('A')
    if o.mapping:
        out.append('m')
    if o.hydrogens:
        out.append('h')
    if not o.bonds:
        out.append('!b')
    if not o.charges:
        out.append('!z')
    if not o.cxsmiles:
        out.append('!x')
    if o.asymmetric_closure:
        out.append('a')
    return ''.join(out)


cdef int smw_random_positions(smw_scratch_t *s, uint32_t n) except -1:
    """A fresh uniformly random `pos`, drawn from the `random` module's shared generator.

    `random.shuffle` and not a Fisher-Yates written here: a caller who seeds `random` wants the same
    batch of strings back, which holds only while the draws come from that generator.  The import is
    local because the core has exactly one module-level Python import (`warn`, see `_core.pyx`) and a
    key nobody in a hot loop uses does not earn a second.
    """
    cdef uint32_t i
    cdef list order = list(range(n))
    from random import shuffle
    shuffle(order)
    for i in range(n):
        s.pos[i] = <uint32_t> <Py_ssize_t> order[i]
    return 0


cdef dict smw_prepare(MoleculeContainer molecule, smw_scratch_t *s,
                      smw_opts_t *o, smw_cuts_t *cuts, smw_sticky_t *sticky):
    """Fill `pos`, `bypos`, the adjacency, the traversal and the closure numbers.

    Everything between "here is a molecule" and "here is a fully decided traversal", so that
    `write_smiles` and the traversal probe cannot decide it differently.  Returns the stereo-group
    map `smw_cxsmiles` wants, or None when the molecule carries no groups.

    `cuts` is NULL for a whole-molecule write, which is the only difference between the two: a
    detached fragment is not a second pipeline but the same one with a `dropped` mask and four extra
    half-edges classified, so nothing can be canonical in one and not the other.

    `sticky` is NULL unless the caller named an end (§14), and it replaces THE TRAVERSAL ONLY -- the
    canonical positions are still computed and still decide every tie, so a sticky string is as much a
    function of `pos` as any other, with two atoms' places in it forced.  It is not canonical, because
    canonical means "the same for every way of building this molecule" and this one also depends on
    which atoms the caller named.
    """
    cdef Structure structure = molecule._structure
    cdef uint32_t n = s.n
    cdef uint32_t i
    cdef dict groups = None
    cdef object gkey, gsid
    cdef list gsids

    if cuts is not NULL:
        # FIRST, before the canonical order: the order is the expensive part of a write and a refused
        # cut list has no use for one.  Nothing here mutates the arena, so the pointers it takes
        # cannot go stale under it.
        smw_apply_cuts(molecule, structure, s, cuts)
    if o.stereo:
        # BEFORE any pointer into the arena is taken: this can append the stereo-unit segment and
        # therefore MOVE the buffer (ruling F60).  The UNMARKED build, because `smw_sign_of` reads
        # constitution -- kind, refs, the unnamed mask -- plus the anchor's own parity byte, and the
        # stereogenic mark is not one of its inputs (ruling F70).  Building it here rather than
        # lazily also keeps `smw_sign_of` honestly `noexcept nogil`.
        ensure_stereo_units_unmarked(structure)
    if o.random_order:
        # `r` REPLACES THE POSITIONS AND NOTHING ELSE.  Everything downstream -- the adjacency sort, the
        # traversal, the closure numbers, the stereo signs -- reads `pos` as a ranking and does not care
        # where the ranking came from, so a random one gives a valid SMILES of the same molecule by the
        # same code path.  `o.canonical` stays True under `r` (the two keys are mutually exclusive), so
        # the stereo groups below are still numbered canonically: a group is chemistry, not order.
        smw_random_positions(s, n)
    elif o.canonical:
        # `smw_canonical_positions` and not `molecule.canonical_order()`: the C entry point takes the
        # STEREO SEED argued at its own comment block, and it fills `pos` by slot rather than
        # building a {n: position} dict for this function to invert.
        #
        # CONTINGENT FOR AROMATIC MOLECULES, owned by the arena epic: with order 4 storable, an
        # aromatic molecule's labelling here is canonical only if the atom and bond invariants behind
        # this order treat order 4 as its own value -- and a wrong labelling would produce a
        # canonical-LOOKING string, which is indistinguishable from a correct one.  The
        # writer does not work around it (kekulising the input would make this a mutator, note 1 at
        # the top of the file); `format(mol, 'i')` is the stored-order escape that claims nothing.
        smw_canonical_positions(structure, s.pos, o.stereo)
    else:
        for i in range(n):
            s.pos[i] = i
    if o.canonical and structure.header.segments[SEG_STEREO_GROUPS].length:
        # {(kind, canonical_id): [n, ...]} inverted to {n: (kind, id)}: the
        # writer asks per atom, the view answers per group, and the inversion is total because
        # an atom belongs to at most one group.
        groups = {}
        for gkey, gsids in molecule.canonical_stereo_groups().items():
            for gsid in gsids:
                groups[gsid] = gkey
    for i in range(n):
        s.bypos[s.pos[i]] = i
    smw_sort_adjacency(structure, s)
    if sticky is not NULL:
        smw_sticky_path(structure, s, sticky)
        smw_sticky_traverse(structure, s, sticky)
        smw_sticky_frames(structure, s, sticky)
    else:
        smw_traverse(structure, s)
    smw_closures(structure, s, cuts)
    smw_directions(structure, s, o)
    if sticky is not NULL:
        smw_sticky_ends(structure, s, sticky)
    return groups


cdef int smw_build_cuts(MoleculeContainer molecule, object cuts, object reserve,
                        smw_cuts_t *out) except -1:
    """`{attachment_id: (keep_n, drop_n)}` and an id set as a filled `smw_cuts_t`.

    The stable-id-to-slot boundary, and nothing else: every question about whether the cuts make sense
    as a CUT is `smw_apply_cuts`', because those answers need the bonds.  What is checked here is what
    can be checked without them -- the id range, the atoms' existence, the pair's shape.

    A dict and not a list of triples, because an id used twice would be a ring closure rather than two
    attachments and a mapping cannot express one.
    """
    cdef dict index = molecule._index_of
    cdef object key, pair, keep, drop
    cdef uint32_t i = 0
    memset(out, 0, sizeof(smw_cuts_t))
    if not isinstance(cuts, dict):
        raise TypeError('cuts must be a {attachment_id: (keep_n, drop_n)} mapping')
    if len(<dict> cuts) > SMW_MAX_CLOSURE - SMW_MIN_ATTACH + 1:
        raise ValueError('at most %d cuts; the attachment ids are %d..%d and one id is one cut'
                         % (SMW_MAX_CLOSURE - SMW_MIN_ATTACH + 1, SMW_MIN_ATTACH, SMW_MAX_CLOSURE))
    for key in sorted(<dict> cuts):
        # SORTED so that a refusal names the same cut every time.  Nothing downstream depends on the
        # order -- the marks and the walk are order-independent -- but an error message that moved with
        # a dict's insertion order would be a test nobody can write.
        if not isinstance(key, int) or isinstance(key, bool):
            raise TypeError('attachment id %r is not an int' % (key,))
        if <int> key < SMW_MIN_ATTACH or <int> key > SMW_MAX_CLOSURE:
            raise ValueError('attachment id %d is outside %d..%d' % (<int> key, SMW_MIN_ATTACH,
                                                                     SMW_MAX_CLOSURE))
        pair = (<dict> cuts)[key]
        if not isinstance(pair, tuple) or len(<tuple> pair) != 2:
            raise TypeError('cut %d must be a (keep_n, drop_n) pair, not %r'
                            % (<int> key, pair))
        keep = (<tuple> pair)[0]
        drop = (<tuple> pair)[1]
        if keep not in index:
            raise KeyError('atom %r is not in this molecule' % (keep,))
        if drop not in index:
            raise KeyError('atom %r is not in this molecule' % (drop,))
        out.ids[i] = <uint8_t> <int> key
        out.keep[i] = <uint32_t> index[keep]
        out.drop[i] = <uint32_t> index[drop]
        i += 1
    out.ncuts = i
    if reserve is not None:
        # The ids of OTHER fragments in the same join.  Withheld from this fragment's internal closure
        # numbers, because a join is a concatenation and a number in scope anywhere in the result is in
        # scope everywhere in it.
        for key in reserve:
            if not isinstance(key, int) or isinstance(key, bool):
                raise TypeError('reserved id %r is not an int' % (key,))
            if <int> key < SMW_MIN_ATTACH or <int> key > SMW_MAX_CLOSURE:
                raise ValueError('reserved id %d is outside %d..%d' % (<int> key, SMW_MIN_ATTACH,
                                                                       SMW_MAX_CLOSURE))
            out.reserved[<int> key] = 1
    return 0


def smw_traversal(MoleculeContainer molecule not None, str spec='', cuts=None, reserve=None):
    """The traversal, before any token is written, for the test suite.

    `{'order': (n, ...), 'tree': ((parent, child), ...), 'closures': ((a, b, number), ...),
    'directions': {anchor_n: (n or None, ...)}, 'lost': (anchor_n, ...),
    'tokens': {(from_n, to_n): '/' or '\\'}, 'unknown_h': (n, ...)}`,
    with `order` the emission order and `tree`/`closures` each edge once.  `directions` is the
    four-direction list per anchor whose sign this writer computes, in stable ids -- from
    `smw_direction_order` for a tetrahedral centre and from `smw_allene_order` for an axis, which for
    the axis means the four are the TERMINALS' directions and not the anchor's own.  It is empty
    under `!s` because the unit table is not built then.  Exposed because the
    invariants of §3 -- every atom emitted exactly once, every bond classified exactly once as
    either a tree edge or a closure -- are properties of the traversal and not of the string, and a
    test that could only read the string would have to infer them.

    `tokens` is `{(from_n, to_n): '/' or '\\'}` for BOTH halves of every directional
    bond, so a test can read `up(terminal -> substituent)` straight out of it instead of parsing the
    string for a character whose meaning depends on which end came first.  That is what makes the
    cis/trans anti-drift test possible: it compares this against `translate_stereo`, and neither side
    reconstructs the other's answer.

    `lost` is the configurations the string does NOT carry, by anchor, sorted by emission order: an
    atropisomer always (SMILES has no syntax for one), an axis whose chain or refs the arena cannot
    confirm, an axis under `!b` (no `=` tokens, so no axis to read the sign against), and a cis/trans
    unit whose stated configuration contradicts another one it shares a bond with.  It is a
    report and not a warning -- nothing raises -- because a writer that refused would leave the
    caller unable to see the structure they actually hold.

    `unknown_h` is every atom whose implicit count is H_UNKNOWN, in emission order, because SMILES has
    no spelling for an unstated count ANYWHERE: a bracket's absent H term means zero and a bare symbol
    means "derive it from the valence model".  So both spellings degrade the fact, and they degrade it
    differently -- a bare atom to the valence-derived count, which is what a caller most likely wants,
    and a bracketed one to zero, which is a number the molecule never claimed.  Which of the two
    happened is visible in the string itself, so it is not reported twice here.  A stated configuration
    on such an atom is additionally in `lost` when the missing count moves a written position
    (`smw_h_frame_unknown`).

    `cuts` and `reserve` are `detached_smiles`' arguments, and they are here because the cut model's
    own invariants -- exactly the boundary edges are attachments, no dropped atom is emitted, the
    attachment numbers are disjoint from the closure numbers -- are properties of the TRAVERSAL and a
    test that could only read the string would have to infer them from it.  `attachments` is then
    `((keep_n, drop_n, id), ...)` in emission order of the retained atom.
    """
    molecule._require_clean()
    cdef Structure structure = molecule._structure
    cdef uint32_t n = structure.header.atom_count
    cdef smw_opts_t o
    smw_parse_spec(spec, &o)
    if n == 0:
        return {'order': (), 'tree': (), 'closures': (), 'directions': {}, 'lost': (),
                'tokens': {}, 'unknown_h': (), 'attachments': ()}

    cdef smw_scratch_t s
    cdef atom_t *atoms
    cdef halfedge_t *edges
    cdef stereo_unit_t *u
    cdef uint32_t want[4]
    cdef uint32_t perm[4]           # `smw_allene_order`'s out-parameter; the probe wants only `want`
    cdef uint32_t i, k, he
    cdef list order = []
    cdef list tree = []
    cdef list closures = []
    cdef list lost = []
    cdef list unknown_h = []
    cdef list attachments = []
    cdef dict directions = {}
    cdef dict tokens = {}
    cdef list one
    cdef smw_cuts_t cutbuf
    cdef smw_cuts_t *cp = NULL
    if cuts is not None:
        smw_build_cuts(molecule, cuts, reserve, &cutbuf)
        cp = &cutbuf
    smw_scratch_alloc(&s, n, csr_ptr(structure)[n])
    try:
        # Both pointers AFTER `smw_prepare`, which builds the stereo-unit table and can move the
        # arena (ruling F60).  `smw_emit` takes its own for the same reason.
        smw_prepare(molecule, &s, &o, cp, NULL)
        atoms = structure.atoms()
        edges = csr_edges(structure)
        for i in range(s.nseq):
            order.append(atoms[s.seq[i]].n)
        for i in range(n):
            for k in range(s.nbr_off[i], s.nbr_off[i + 1]):
                he = s.wnbr[k]
                if s.he_kind[he] == SMW_HE_CHILD:
                    tree.append((atoms[i].n, atoms[edges[he].to].n))
                elif s.he_kind[he] == SMW_HE_CLOSURE and \
                        s.out_idx[i] < s.out_idx[edges[he].to]:
                    closures.append((atoms[i].n, atoms[edges[he].to].n,
                                     s.he_close[he]))
        for i in range(s.nseq):          # emission order, and retained atoms only by construction
            for k in range(s.nbr_off[s.seq[i]], s.nbr_off[s.seq[i] + 1]):
                he = s.wnbr[k]
                if s.he_kind[he] == SMW_HE_ATTACH:
                    attachments.append((atoms[s.seq[i]].n, atoms[edges[he].to].n,
                                        s.he_close[he]))
        if o.stereo:
            for i in range(n):
                if s.dropped[i]:
                    continue
                u = stereo_unit_of(structure, i)
                if u is NULL or u.n_refs != 4:
                    continue
                if u.kind == SU_ALLENE:
                    # The axis's four, grouped by terminal -- `smw_allene_order`'s tuple, which is the
                    # one the sign was translated into.  A refusal is reported through `lost` and
                    # gets no entry here, so a test can tell "no order" from "an order with no sign".
                    if not smw_allene_order(structure, &s, i, u, want, perm):
                        continue
                elif u.kind == SU_TETRA:
                    smw_direction_order(structure, &s, i, want)
                else:
                    continue
                one = []
                for k in range(4):
                    one.append(None if want[k] == SU_NO_REF
                               else atoms[want[k]].n)
                directions[atoms[i].n] = tuple(one)
        for i in range(s.nseq):
            if s.lost[s.seq[i]]:
                lost.append(atoms[s.seq[i]].n)
            if at_implicit_h_unknown(&atoms[s.seq[i]]):
                unknown_h.append(atoms[s.seq[i]].n)
        for i in range(n):
            for k in range(s.nbr_off[i], s.nbr_off[i + 1]):
                he = s.wnbr[k]
                if s.he_dir[he]:
                    tokens[(atoms[i].n, atoms[edges[he].to].n)] = \
                        '/' if s.he_dir[he] == 1 else '\\'
    finally:
        smw_scratch_free(&s)
    return {'order': tuple(order), 'tree': tuple(tree), 'closures': tuple(closures),
            'directions': directions, 'lost': tuple(lost), 'tokens': tokens,
            'unknown_h': tuple(unknown_h), 'attachments': tuple(attachments)}


def write_smiles(MoleculeContainer molecule not None, str spec='', bint return_order=False):
    """The molecule as a SMILES string.

    Canonical by default: the same molecule written from any creation order gives the same string,
    because every choice the writer makes is a function of the canonical positions -- which under
    `s` are seeded with the configuration, so that a constitutional symmetry the CONFIGURATION breaks
    is broken in the order too (`smw_stereo_seed`).

    The representation is the molecule's own: a stored order-4 bond writes lowercase, a stored Kekule
    bond writes `=`, and nothing here converts between them.  The first note at the top of
    `_smiles_write.pxi` is why.

    `spec` is the `format()` spec described at `smw_parse_spec`.

    `return_order=True` answers `(string, (n, ...))` -- the atoms in the order the string
    writes them.  It exists because a caller that needs both must not run the traversal twice: the
    order is a function of the whole option set, so a second call with a different spec would return
    an order that does not describe the string in hand.  The consumers that cannot be served by the
    string alone are the CXSMILES tail of a REACTION, whose radical and stereo indices count atoms
    across every molecule in it, and `smiles_atoms_order`.
    """
    molecule._require_clean()
    cdef Structure structure = molecule._structure
    cdef uint32_t n = structure.header.atom_count
    cdef smw_opts_t o
    smw_parse_spec(spec, &o)
    if n == 0:
        return ('', ()) if return_order else ''

    cdef smw_scratch_t s
    cdef smw_buf_t b
    cdef dict groups
    cdef object out
    cdef list order
    cdef atom_t *atoms
    cdef uint32_t i

    smw_scratch_alloc(&s, n, csr_ptr(structure)[n])
    b.data = NULL
    b.length = 0
    b.cap = 0
    b.oom = False
    try:
        groups = smw_prepare(molecule, &s, &o, NULL, NULL)
        smw_emit(structure, &s, &b, &o, NULL)
        if o.cxsmiles:
            smw_cxsmiles(structure, &s, &b, groups, molecule.aliases)
        if b.oom:
            raise MemoryError()
        out = b.data[:b.length].decode('ascii')
        if return_order:
            # Re-borrowed here rather than held from before `smw_prepare`, which can append the
            # stereo-unit segment and move the buffer (ruling F60).
            atoms = structure.atoms()
            order = []
            for i in range(s.nseq):
                order.append(atoms[s.seq[i]].n)
            out = (out, tuple(order))
    finally:
        PyMem_Free(b.data)
        smw_scratch_free(&s)
    return out


# ------------------------------------------------------------------------------------------------
# THE REACTION.  Three sides, one string, ONE tail -- and the tail is why the reaction writer cannot
# be three calls to `write_smiles` with `>` between them.  A CXSMILES tail's atom indices count from
# the start of the whole string, so a per-molecule tail is stranded mid-string the moment anything is
# written after it, and two molecules' `&1` are two different AND groups that would silently merge.
#
# So each molecule is written WITHOUT a tail and hands back its tail's STRUCTURE, exactly as a
# detached fragment does for `DetachedSmiles.join`; the aggregation below is the same aggregation that
# method performs, for the same reason and with the same renumbering.

cdef tuple smw_reaction_part(MoleculeContainer molecule, smw_opts_t *o):
    """One molecule of a reaction: `(text, ids, tail parts, components, map numbers)`, ids and map
    numbers both in WRITTEN order.

    A REACTION-SHAPED `write_smiles`.  The text carries no tail of its own and the tail comes back
    unformatted, which is the one thing `write_smiles` cannot answer: by the time it returns, its tail
    is characters and its indices are local to it.  Everything between the two calls is `smw_prepare`
    and `smw_emit`, so no traversal decision is made twice and a molecule of a reaction is written by
    the same canonical writer as a molecule on its own.
    """
    molecule._require_clean()
    cdef Structure structure = molecule._structure
    cdef uint32_t n = structure.header.atom_count
    cdef uint32_t components = molecule.connected_components_count if n else 1
    if n == 0:
        return ('', (), ([], [], {}, {}, []), 1, ())

    cdef smw_scratch_t s
    cdef smw_buf_t b
    cdef dict groups
    cdef str out
    cdef list order, maps
    cdef tuple parts
    cdef atom_t *atoms
    cdef uint32_t i

    smw_scratch_alloc(&s, n, csr_ptr(structure)[n])
    b.data = NULL
    b.length = 0
    b.cap = 0
    b.oom = False
    try:
        groups = smw_prepare(molecule, &s, o, NULL, NULL)
        smw_emit(structure, &s, &b, o, NULL)
        if b.oom:
            raise MemoryError()
        out = b.data[:b.length].decode('ascii')
        parts = smw_tail_parts(structure, &s, groups, molecule.aliases) if o.cxsmiles \
            else ([], [], {}, {}, [])
        # re-borrowed after `smw_prepare`, which can append the stereo-unit segment and move the
        # buffer -- the same ruling `write_smiles` cites
        atoms = structure.atoms()
        order = []
        maps = []
        for i in range(s.nseq):
            order.append(atoms[s.seq[i]].n)
            maps.append(<int> atoms[s.seq[i]].map_number)
    finally:
        PyMem_Free(b.data)
        smw_scratch_free(&s)
    return (out, tuple(order), parts, components, tuple(maps))


cdef int smw_collect_groups(dict groups, tuple maps, int atom_base, list out) except -1:
    """One molecule's groups of one kind, appended as `(shifted indices, map numbers)`."""
    cdef object key, idx
    cdef list shifted
    cdef set numbers
    for key in sorted(groups):
        shifted = []
        numbers = set()
        for idx in <list> groups[key]:
            shifted.append(<int> idx + atom_base)
            if <int> maps[<int> idx]:
                numbers.add(<int> maps[<int> idx])
        out.append((shifted, numbers))
    return 0


cdef Py_ssize_t smw_group_root(list parent, Py_ssize_t i):
    """Union-find with path halving.  A rank is not worth carrying for a handful of groups."""
    while <Py_ssize_t> parent[i] != i:
        parent[i] = parent[<Py_ssize_t> parent[i]]
        i = <Py_ssize_t> parent[i]
    return i


cdef dict smw_merge_groups(list pending):
    """`[(indices, map numbers)]` collapsed to `{group id: indices}`, one id per correlated SET.

    Two molecules' `&1` are two different groups -- unless the atom-atom mapping says their members are
    the same atoms restated, which is exactly a shared map number.  Merging there is what makes one
    racemic centre carried across the arrow come out as one group; a group whose members are all
    unmapped shares with nothing and keeps an id of its own.  Ids are handed out in written order, so
    the string stays a function of the reaction.
    """
    cdef Py_ssize_t total = len(pending)
    cdef list parent = list(range(total))
    cdef dict owner = {}, ids = {}, out = {}
    cdef Py_ssize_t i, a, b
    cdef object number
    for i in range(total):
        for number in <set> (<tuple> pending[i])[1]:
            if number in owner:
                a = smw_group_root(parent, i)
                b = smw_group_root(parent, <Py_ssize_t> owner[number])
                if a != b:
                    parent[b] = a
            else:
                owner[number] = i
    for i in range(total):
        a = smw_group_root(parent, i)
        if a not in ids:
            ids[a] = len(ids) + 1
            out[<int> ids[a]] = []
        (<list> out[<int> ids[a]]).extend(<list> (<tuple> pending[i])[0])
    return out


def write_reaction_smiles(rxn not None, str spec=''):
    """The reaction as a reaction SMILES: `reactants>agents>products`, with one CXSMILES tail.

    **Each side\'s molecules are sorted by their own string**, so the same reaction assembled in any
    order gives one string -- the identifier property `write_smiles` already promises per molecule,
    lifted to the reaction.  `!c` keeps the container\'s order instead, and every other key is passed
    to `write_smiles` unchanged.

    The tail aggregates `^1:` radicals, the three enhanced-stereo group fields and `f:`, whose groups
    name the components of every molecule that has more than one -- so `[Na+].[Cl-]` on one side comes
    back as one reactant and not two.  A tail carrying `^1:` and `f:` alone drops enhanced stereo
    groups from every reaction it writes.  `!x` suppresses the whole block.
    """
    cdef smw_opts_t o
    cdef bint keep = '!c' in spec
    smw_parse_spec(spec.replace('!c', ''), &o)

    cdef list side_texts = []
    cdef list radicals = [], abs_atoms = [], fgroups = [], labels = []
    cdef list and_pending = [], or_pending = []
    cdef uint32_t atom_base = 0, component_base = 0, components, i, j
    cdef list rows, keys, texts, group, order
    cdef tuple row, tail, maps
    cdef object side, mol, key, other
    cdef str text

    for side in (rxn.reactants, rxn.agents, rxn.products):
        rows = []
        for mol in side:
            rows.append(smw_reaction_part(<MoleculeContainer> mol, &o))
        if not keep:
            # SORTED ON THE TEXT ALONE, with the position as the tiebreaker, and never on the rows
            # themselves: a row holds the tail\'s dicts, and two equal molecules in one side would
            # reach them and fail to compare.  The position keeps the sort total and, being the
            # container\'s order, keeps it stable for a repeated molecule.
            keys = []
            i = 0
            for row in rows:
                keys.append((<str> row[0], i))
                i += 1
            keys.sort()
            texts = rows
            rows = []
            for key in keys:
                rows.append(texts[<uint32_t> key[1]])

        texts = []
        for row in rows:
            text = <str> row[0]
            order = list(<tuple> row[1])
            tail = <tuple> row[2]
            components = <uint32_t> row[3]
            maps = <tuple> row[4]
            texts.append(text)
            if o.cxsmiles:
                if components > 1:
                    # the components of one molecule are written consecutively, so the group is a run
                    group = []
                    for j in range(components):
                        group.append(component_base + j)
                    fgroups.append(group)
                for key in <list> tail[0]:
                    radicals.append(<int> key + <int> atom_base)
                for key in <list> tail[1]:
                    abs_atoms.append(<int> key + <int> atom_base)
                # CONCATENATED, not shifted: the labels are positional, one entry per atom, so a
                # molecule that carries none still owes the field its own atoms' worth of blanks
                labels.extend(<list> tail[4] if <list> tail[4] else [None] * len(order))
                # MERGED WHERE THE MAPPING SAYS SO, RENUMBERED WHERE IT DOES NOT.  Two molecules\' `&1`
                # are two AND groups and keeping both numbers would claim their atoms invert together
                # -- but a centre carried across the arrow is ONE group, and a second id there asserts
                # two independently racemic centres where there is one.  Resolved after all three
                # sides, a reactant\'s partner being written later.
                smw_collect_groups(<dict> tail[2], maps, <int> atom_base, and_pending)
                smw_collect_groups(<dict> tail[3], maps, <int> atom_base, or_pending)
            component_base += components
            atom_base += <uint32_t> len(order)
        side_texts.append('.'.join(texts))

    cdef str body = '>'.join(side_texts)
    cdef dict and_groups, or_groups
    if o.cxsmiles:
        and_groups = smw_merge_groups(and_pending)
        or_groups = smw_merge_groups(or_pending)
        return body + smw_tail_text((radicals, abs_atoms, and_groups, or_groups, labels), fgroups)
    return body


# ------------------------------------------------------------------------------------------------
# DETACHED SMILES (spec §13).  A fragment whose cut bonds are RING BONDS, so that fragments re-join by
# CONCATENATION -- no string surgery, no re-parse, no search for a spelling that happens to work.
#
# Joining TEXT is what forces the alternative's constraints: the two attachment atoms have to land at
# the two ENDS of the string, which a randomised non-canonical order reaches only by chance, the cut is
# `smiles[2:]` and `smiles[:-2]`, exactly two attachments are expressible, a second component is not,
# and a terminal atom inside a ring has no spelling that survives the cut.  A ring bond is the
# notation's own mechanism for "this bond's other end is elsewhere", so none of that applies here: the
# attachment can be any atom, there can be up to 90 of them, and the fragment is written by the one
# canonical writer everything else uses.
cdef class DetachedSmiles:
    """A SMILES fragment with open attachment ids, and the pieces a join needs.

    `text` is the fragment WITHOUT its CXSMILES tail, because a join concatenates bodies and then
    writes ONE tail for the result -- an index in a tail counts atoms from the start of the whole
    string, so a tail is meaningless in the middle of one.  `tail` is therefore kept as STRUCTURE
    (`smw_tail_parts`' five collections) and formatted by `smw_tail_text`, the same function
    `write_smiles` uses.  Rebuilding a tail by editing its text would mean parsing the format this
    module also writes, and two readings of one syntax are how they come apart.

    `str(fragment)` IS NOT A SMILES while any id is open: `%12` with no partner is a dangling ring
    bond and every reader rejects it, which is the point -- a fragment cannot be mistaken for a
    molecule.  Once `open_ids` is empty the text is a valid SMILES for the joined molecule, though not
    a canonical one: canonical means "the same for every way of building this molecule", and a join
    keeps the fragments' own orders.  So a `DetachedSmiles` IS NOT A CACHE KEY and is not a spec key --
    `format(mol, ...)` has no letter for it, deliberately.  Ask for one and you get a function call.
    """
    cdef readonly str text
    cdef readonly tuple order          # stable ids, in the order the text writes them
    cdef readonly tuple open_ids       # attachment ids with no partner in this text, ascending
    cdef readonly tuple closure_ids    # ring numbers the text already uses for its OWN closures
    cdef readonly tuple lost           # anchors whose configuration the text does not carry
    cdef readonly tuple unknown_h      # atoms whose implicit count nobody stated
    cdef readonly tuple tail           # (radicals, abs, and_groups, or_groups, labels); do not mutate

    def __init__(self, str text not None, tuple order not None, tuple open_ids not None,
                 tuple closure_ids not None, tuple lost not None, tuple unknown_h not None,
                 tuple tail not None):
        self.text = text
        self.order = order
        self.open_ids = open_ids
        self.closure_ids = closure_ids
        self.lost = lost
        self.unknown_h = unknown_h
        self.tail = tail

    @property
    def atom_count(self):
        """The atoms this text writes -- the shift a following fragment's tail indices need."""
        return len(self.order)

    def __str__(self):
        return self.text + smw_tail_text(self.tail)

    def __repr__(self):
        return 'DetachedSmiles(%r, open=%r)' % (self.text, self.open_ids)

    @staticmethod
    def join(*fragments):
        """The fragments as one `DetachedSmiles`, bodies concatenated as `.` components.

        `.` and not nothing between them: each fragment keeps being its own component of the text, and
        the ring bonds are what make it one molecule.  That is also what makes the join safe atom by
        atom -- a fragment's first atom stays a component leader, so the implicit hydrogen in
        `[C@H](%12)F` sits in the same written position after the join as before it, and every
        tetrahedral sign keeps its meaning (`smw_direction_order` argues this at length).

        An id in exactly two fragments is CLOSED by the join and becomes an ordinary ring closure of
        the result; an id in one stays open, so a molecule can be built up in stages.  Three
        occurrences is refused: a third `%12` has nothing to bond to.
        """
        cdef DetachedSmiles f
        cdef list parts = [], order = [], lost = [], unknown_h = []
        cdef list radicals = [], abs_atoms = [], labels = []
        cdef dict and_groups = {}, or_groups = {}
        cdef dict counts = {}
        cdef set closures = set()
        cdef list still_open, clash, shifted
        cdef uint32_t shift = 0
        cdef object frag, key, other
        cdef int and_next = 0, or_next = 0

        if not fragments:
            raise ValueError('join needs at least one fragment')
        for frag in fragments:
            if not isinstance(frag, DetachedSmiles):
                raise TypeError('join takes DetachedSmiles, not %r' % (type(frag).__name__,))
            for key in (<DetachedSmiles> frag).open_ids:
                counts[key] = counts.get(key, 0) + 1
        # EXPLICIT LOOPS AND NOT COMPREHENSIONS throughout this method: a comprehension gets its own
        # scope in Cython 3, `warn.undeclared` is on for this module, and a target declared in the
        # enclosing function does not satisfy it.  The loops are the same length anyway.
        clash = []
        for key in sorted(counts):
            if <int> counts[key] > 2:
                clash.append(key)
        if clash:
            raise ValueError('attachment id(s) %s appear in more than two fragments; a ring bond '
                             'joins exactly two atoms' % ', '.join(map(str, clash)))
        for frag in fragments:
            f = <DetachedSmiles> frag
            # THE COLLISION `reserve` EXISTS FOR.  An internal closure numbered the same as an
            # attachment open anywhere in the join would be paired with it by the reader -- silently,
            # producing a valid molecule that is not this one.  Refused rather than renumbered: the
            # texts are already written, and renumbering would mean editing them.
            clash = sorted(set(f.closure_ids) & set(counts))
            if clash:
                raise ValueError('a fragment already uses ring number(s) %s for its own closures, and '
                                 'they are attachment ids in this join; pass reserve=%r when writing '
                                 'it' % (', '.join(map(str, clash)), tuple(sorted(counts))))
        for frag in fragments:
            f = <DetachedSmiles> frag
            parts.append(f.text)
            order.extend(f.order)
            lost.extend(f.lost)
            unknown_h.extend(f.unknown_h)
            closures.update(f.closure_ids)
            for key in <list> f.tail[0]:
                radicals.append(<int> key + <int> shift)
            for key in <list> f.tail[1]:
                abs_atoms.append(<int> key + <int> shift)
            # CONCATENATED, not shifted: `$...$` is positional, one entry per atom, so a fragment
            # with no label still owes the field a blank for each of its own atoms
            labels.extend(<list> f.tail[4] if <list> f.tail[4] else [None] * len(f.order))
            # RENUMBERED, not merged: two fragments' `&1` are two different AND groups, and keeping the
            # numbers would silently claim their atoms invert together.  Sequential in fragment order,
            # which is the only numbering available -- the result is not canonical anyway.
            for key in sorted(<dict> f.tail[2]):
                and_next += 1
                shifted = []
                for other in (<dict> f.tail[2])[key]:
                    shifted.append(<int> other + <int> shift)
                and_groups[and_next] = shifted
            for key in sorted(<dict> f.tail[3]):
                or_next += 1
                shifted = []
                for other in (<dict> f.tail[3])[key]:
                    shifted.append(<int> other + <int> shift)
                or_groups[or_next] = shifted
            shift += <uint32_t> len(f.order)
        still_open = []
        for key in sorted(counts):
            if <int> counts[key] == 1:
                still_open.append(key)
            else:
                closures.add(key)
        return DetachedSmiles('.'.join(parts), tuple(order), tuple(still_open),
                              tuple(sorted(closures)), tuple(lost), tuple(unknown_h),
                              (radicals, abs_atoms, and_groups, or_groups, labels))


def detached_smiles(MoleculeContainer molecule not None, cuts not None, str spec='',
                    reserve=None):
    """The molecule minus the dropped side of every cut, with the cut bonds as open ring bonds.

    `cuts` is `{attachment_id: (keep_n, drop_n)}`.  ORDERED pairs, because there is
    nothing in `C-C` to say which half the caller wants; `attachment_id` is 10..99 and is written
    `%NN`, never a bare digit, so that an attachment is machine-findable in the text and cannot be
    confused with an ordinary closure.  Ten as the floor is a MEASUREMENT: `%05` is refused by RDKit
    and by chython 2, so a fixed-width low spelling is not available.

    `reserve` is the other fragments' attachment ids when this one is written for a join, withheld from
    this fragment's own closure numbers.

    Refused, all naming atoms: a cut whose two atoms are not bonded, a bond in a RING (one id cannot
    carry both ends of a ring opening), any bond crossing into the dropped part that no cut names, and
    an empty retained set.  The last two are one rule -- the cut list must be exactly the boundary --
    and it is what makes the dropped side a decision of the caller's rather than a consequence of
    graph reachability, so a salt keeps its counter-ion.
    """
    molecule._require_clean()
    cdef Structure structure = molecule._structure
    cdef uint32_t n = structure.header.atom_count
    cdef smw_opts_t o
    cdef smw_cuts_t c
    smw_parse_spec(spec, &o)
    smw_build_cuts(molecule, cuts, reserve, &c)

    cdef smw_scratch_t s
    cdef smw_buf_t b
    cdef atom_t *atoms
    cdef uint32_t i, k, he
    cdef list order = [], lost = [], unknown_h = []
    cdef set open_ids = set(), closure_ids = set()
    cdef str text
    cdef dict groups
    cdef tuple tail

    if n == 0:
        return DetachedSmiles('', (), (), (), (), (), ([], [], {}, {}, []))
    smw_scratch_alloc(&s, n, csr_ptr(structure)[n])
    b.data = NULL
    b.length = 0
    b.cap = 0
    b.oom = False
    try:
        groups = smw_prepare(molecule, &s, &o, &c, NULL)
        smw_emit(structure, &s, &b, &o, NULL)
        if b.oom:
            raise MemoryError()
        text = b.data[:b.length].decode('ascii')
        # AFTER `smw_prepare`, which can move the arena (ruling F60).
        atoms = structure.atoms()
        tail = smw_tail_parts(structure, &s, groups, molecule.aliases) if o.cxsmiles \
            else ([], [], {}, {}, [])
        for i in range(s.nseq):
            order.append(atoms[s.seq[i]].n)
            if s.lost[s.seq[i]]:
                lost.append(atoms[s.seq[i]].n)
            if at_implicit_h_unknown(&atoms[s.seq[i]]):
                unknown_h.append(atoms[s.seq[i]].n)
            for k in range(s.nbr_off[s.seq[i]], s.nbr_off[s.seq[i] + 1]):
                he = s.wnbr[k]
                if s.he_close[he] == 0:
                    continue
                if s.he_kind[he] == SMW_HE_ATTACH:
                    open_ids.add(s.he_close[he])
                elif s.he_kind[he] == SMW_HE_CLOSURE:
                    closure_ids.add(s.he_close[he])
    finally:
        PyMem_Free(b.data)
        smw_scratch_free(&s)
    return DetachedSmiles(text, tuple(order), tuple(sorted(open_ids)), tuple(sorted(closure_ids)),
                          tuple(lost), tuple(unknown_h), tail)


# ------------------------------------------------------------------------------------------------
# STICKY SMILES (§14).  `sticky_smiles`' consumers live outside this repository, so the STRING SHAPE is
# a contract and not a design choice: the named atoms are the FIRST and LAST tokens and a caller GLUES
# strings together, `A.sticky_right + B.sticky_left`.  That is a weaker notation than
# `detached_smiles`' ring bonds -- it can only open the two ends of a chain, where a cut can open ninety
# bonds anywhere -- and it is the one this entry point owes its callers.
#
# The traversal is CONSTRAINED (`smw_sticky_traverse`, which proves its own termination) rather than
# randomised and retried until an atom lands where it is wanted, and the tokens are suppressed by the
# writer, which knows the bond's order and can therefore write `:` where text surgery writes `-`.
# `tries` is accepted in the Python signature and does nothing.
def sticky_smiles(MoleculeContainer molecule not None, left=None, right=None, str spec='', *,
                  bint remove_left=False, bint remove_right=False,
                  bint keep_bond_left=False, bint keep_bond_right=False, bint report=False):
    """The molecule as a SMILES whose first token is `left`'s and whose last is `right`'s.

    `left` and `right` are stable ids; either may be None, but not both.  `remove_*` suppresses that
    end's ATOM token, leaving the string open for a caller to glue an atom onto, and `keep_bond_*`
    keeps its BOND token -- forced non-empty, so `-`, `=`, `#`, `:`, `~`, or `/` and `\\` when the bond
    carries a direction.  `spec` is `write_smiles`' spec, minus the CXSMILES tail: an index in a tail
    counts atoms from the start of the whole string and a glued string has no such start.

    `keep_bond_*=False` drops the BOND as well as the atom, so what the caller glues on lands on an
    implicit SINGLE bond whatever the record held, and any configuration that referenced the removed
    atom is reported through `lost` rather than written.  An open fragment is also not a standalone
    molecule: an atom token is computed over the whole molecule, so the end's neighbour may re-parse
    with a different implicit-hydrogen count until the caller has glued the atom back.

    NOT CANONICAL and not a cache key: the order depends on the atoms the caller named.

    `report=True` answers `(text, order, lost)`: the written order as stable ids -- which is where
    "starts at `left`, ends at `right`" is CHECKABLE, since a removed end leaves no token to look at --
    and the anchors whose configuration the text does not carry, because a bare `str` has nowhere to put
    that and rule 4 at the top of this file says an unspellable fact is reported, not dropped.

    Refused, and each because the notation cannot say it rather than because the walk failed:

    * neither end named -- there is nothing to make sticky;
    * `left` and `right` the same atom -- one token cannot be both the first and the last;
    * `remove_*` without that end named;
    * `remove_*` on an atom whose degree is not 1 -- `L(A)B` minus `L` is `(A)B`, and an atom with a
      ring closure would leave the digits dangling;
    * `remove_*` at both ends of a two-atom molecule -- no atom would be left;
    * `right` named on a molecule with more than one component -- the string would have to end in the
      middle of it;
    * `right` a CUT VERTEX -- some atom is then reachable only through it, so no walk can leave it for
      last.  This, and not terminality, is the real precondition: a non-terminal `right` whose removal
      leaves the rest connected is fine, and `...C%12` is a good last token as long as nobody deletes it.
    """
    molecule._require_clean()
    cdef Structure structure = molecule._structure
    cdef uint32_t n = structure.header.atom_count
    cdef dict index = molecule._index_of
    cdef uint32_t *ptr
    cdef smw_opts_t o
    cdef smw_sticky_t k
    smw_parse_spec(spec, &o)
    o.cxsmiles = False

    if left is None and right is None:
        raise ValueError('either left or right atom should be specified')
    if left is not None and left not in index:
        raise KeyError('atom %r is not in this molecule' % (left,))
    if right is not None and right not in index:
        raise KeyError('atom %r is not in this molecule' % (right,))
    if left is not None and right is not None and left == right:
        raise ValueError('left and right name the same atom %r; one token cannot be both the first '
                         'and the last' % (left,))
    if remove_left and left is None:
        raise ValueError('remove_left needs a left atom to remove')
    if remove_right and right is None:
        raise ValueError('remove_right needs a right atom to remove')
    if right is not None and molecule.connected_components_count > 1:
        raise ValueError('right=%r was named on a molecule with %d components; the string would have '
                         'to end in the middle of it.  Write the components separately, or pass left '
                         'only' % (right, molecule.connected_components_count))
    k.left = <uint32_t> index[left] if left is not None else SMW_NONE
    k.right = <uint32_t> index[right] if right is not None else SMW_NONE
    k.remove_left = remove_left
    k.remove_right = remove_right
    k.keep_bond_left = keep_bond_left
    k.keep_bond_right = keep_bond_right
    ptr = csr_ptr(structure)
    if remove_left and ptr[k.left + 1] - ptr[k.left] != 1:
        raise ValueError('remove_left=True on atom %r, whose degree is %d: only a terminal atom can be '
                         'removed from the string, or its branches and ring closures would be left '
                         'dangling' % (left, ptr[k.left + 1] - ptr[k.left]))
    if remove_right and ptr[k.right + 1] - ptr[k.right] != 1:
        raise ValueError('remove_right=True on atom %r, whose degree is %d: only a terminal atom can '
                         'be removed from the string, or its branches and ring closures would be left '
                         'dangling' % (right, ptr[k.right + 1] - ptr[k.right]))
    if remove_left and remove_right and n == 2:
        raise ValueError('removing both ends of a two-atom molecule leaves no atom to write')

    cdef smw_scratch_t s
    cdef smw_buf_t b
    cdef object out
    cdef atom_t *atoms
    cdef list order, lost
    cdef uint32_t i
    smw_scratch_alloc(&s, n, ptr[n])
    b.data = NULL
    b.length = 0
    b.cap = 0
    b.oom = False
    try:
        if k.right != SMW_NONE and smw_sticky_severs(structure, &s, k.right):
            raise ValueError('right=%r is a cut vertex: removing it would break the molecule in two, '
                             'so some atom can only be reached THROUGH it and no string can end there.  '
                             'Name a terminal atom, or one whose removal leaves the rest connected'
                             % (right,))
        memset(s.visited, 0, n)     # borrowed by the check above; the traversal reads it as all-zero
        smw_prepare(molecule, &s, &o, NULL, &k)
        # THE TRAVERSAL'S OWN CLAIM, CHECKED.  `smw_sticky_traverse` proves that `right` is written
        # last, and `smw_sticky_path` stops rather than reporting when its input is not what it was
        # promised -- a `noexcept nogil` walk cannot raise.  So the claim is verified here, where it can
        # be: a string that does not end where the caller asked is worse than no string, because they
        # will glue something onto it.
        if k.right != SMW_NONE and s.seq[s.nseq - 1] != k.right:
            raise ValueError('the traversal did not end at right=%r; this is a writer bug, please '
                             'report the molecule' % (right,))
        smw_emit(structure, &s, &b, &o, &k)
        if b.oom:
            raise MemoryError()
        out = b.data[:b.length].decode('ascii')
        if report:
            # `(text, order, lost)`.  The ORDER because "the string starts at `left` and ends at
            # `right`" is not visible in the text -- a removed end leaves no token there to look at --
            # and `LOST` because rule 4 at the top of this file says an unspellable fact is reported
            # rather than dropped, and the plain string this entry point returns has nowhere to put it.
            # `detached_smiles` carries the same two in `DetachedSmiles`; here they are opt-in, because
            # the frozen `MoleculeContainer.sticky_smiles` signature returns a `str`.
            # Atoms re-borrowed after `smw_prepare`, which can move the arena (ruling F60).
            atoms = structure.atoms()
            order = []
            lost = []
            for i in range(s.nseq):
                order.append(atoms[s.seq[i]].n)
                if s.lost[s.seq[i]]:
                    lost.append(atoms[s.seq[i]].n)
            out = (out, tuple(order), tuple(lost))
    finally:
        PyMem_Free(b.data)
        smw_scratch_free(&s)
    return out
