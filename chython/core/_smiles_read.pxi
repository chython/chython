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
# The SMILES reader: a string becomes a MoleculeContainer, in one pass over the bytes and one
# pass over a flat parse graph, with a single malloc for the whole parse.
#
# THREE DECISIONS THIS FILE IS BUILT ON
#
# 1. Aromatic bonds are STORED AROMATIC.  A lowercase string produces order-4 bonds and this file
#    never calls `kekule()` to convert one.  `kekule` is a deliberate operation the caller runs; a
#    reader that ran it would make "what the string said" unrecoverable.  The single place this
#    file touches `kekule` at all is the promotion fallback below, where it asks a QUESTION of a
#    discarded copy -- see `smi_promote`.
#
# 2. Hydrogen counts come from the SMILES NOTATION model, `smv_default_h` in `_smiles_write.pxi`,
#    not from the chemistry rules in `_valence.pxi`.  The two models answer different questions
#    and must not be merged: `_valence.pxi` knows what an atom's environment permits, the notation
#    model knows what a reader is REQUIRED to infer when a bracket is absent.  Calling the writer's
#    function rather than restating its numbers is what keeps read and write from drifting.
#
# 3. Syntax errors raise; chemistry is stored and logged.  A malformed string has no molecule in it
#    and `IncorrectSmiles` carries the byte offset.  A well-formed string describing something
#    impossible -- a hypervalent atom, a five-membered all-lowercase ring with no Kekule form, an
#    element with no aromatic form -- produces a molecule plus log lines.  This is Ramil's standing
#    constraint: input by default is garbage, no exceptions; store it and say what was repaired.
#
# ATOM-CASE AROMATIC PROMOTION
#
# For each smallest ring whose every atom was written lowercase, every bond of that ring is
# aromatic, whatever order the string wrote.  `c1cccc-c-1` and `c1cccc-c1` and `c1ccccc-1` all
# read as benzene.  The aromatic set comes from atom case and never from bond order, so an
# explicit `-` inside an all-lowercase ring is a preference WITHIN the pi system, not a statement
# that leaves it.
#
# Promotion is a REPAIR and is logged as one, because the stored orders are then not the written
# ones and the caller has to be able to see that.  Biphenyl's inter-ring bond is in no smallest
# ring and is therefore never promoted -- that case is what makes the rule safe, and
# `test_biphenyl_inter_ring_bond_is_not_promoted` pins it.  If the promoted set has no Kekule form
# the promotion is reverted to the stated orders and the revert is logged; the stated set cannot do
# worse than itself, so the fallback is free insurance.
#
# THE CXSMILES TAIL
#
# `^N:` radicals are applied.  Every other field -- coordinates, atom labels, fragment grouping,
# enhanced stereo groups -- is named in the log and not applied, so a caller can see exactly what
# of its input this reader kept.  Nothing in the tail raises: the molecule in front of it is intact
# and refusing it would be the larger loss.  See `smi_cx`.
#
# WHAT IS NOT HERE YET
#
# Stereo (step 5).  Configuration marks in the string are parsed and held in the parse graph but not
# yet applied to the molecule, and this file logs one line saying so rather than dropping them in
# silence.  The tail's enhanced stereo groups (`a:`, `o1:`, `&1:`) land with them.


# Element atomic number by symbol, for the two-character lookup the tokeniser does per atom.
# Indexed (first_char - 'A') * 27 + (0 for a one-letter symbol, else second_char - 'a' + 1); 0
# means no element with that spelling.  A lowercase-first symbol -- `n`, `se`, `as` -- is looked up
# by uppercasing the first character into this same table, so aromatic spelling costs nothing.
# Generated from SYMBOLS in `_elements.pxi`; `test_smi_element_table_matches_symbols` compares the
# two entry by entry so they cannot drift.
cdef extern from *:
    """
    static const unsigned char SMI_ELEMENT[26 * 27] = {
      0,   0,   0,  89,   0,   0,   0,  47,   0,   0,   0,   0,  13,  95,   0,   0,   0,   0,  18,  33,  85,  79,   0,   0,   0,   0,   0,   /* A: Al Ar As Ag Au At Ac Am */
      5,  56,   0,   0,   0,   4,   0,   0, 107,  83,   0,  97,   0,   0,   0,   0,   0,   0,  35,   0,   0,   0,   0,   0,   0,   0,   0,   /* B: Be B Br Ba Bi Bk Bh */
      6,  20,   0,   0,  48,  58,  98,   0,   0,   0,   0,   0,  17,  96, 112,  27,   0,   0,  24,  55,   0,  29,   0,   0,   0,   0,   0,   /* C: C Cl Ca Cr Co Cu Cd Cs Ce Cm Cf Cn */
      0,   0, 105,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 110,   0,   0,   0,   0,   0,  66,   0,   /* D: Dy Db Ds */
      0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  68,  99,   0,  63,   0,   0,   0,   0,   0,   /* E: Eu Er Es */
      9,   0,   0,   0,   0,  26,   0,   0,   0,   0,   0,   0, 114, 100,   0,   0,   0,   0,  87,   0,   0,   0,   0,   0,   0,   0,   0,   /* F: F Fe Fr Fm Fl */
      0,  31,   0,   0,  64,  32,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   /* G: Ga Ge Gd */
      1,   0,   0,   0,   0,   2,  72,  80,   0,   0,   0,   0,   0,   0,   0,  67,   0,   0,   0, 108,   0,   0,   0,   0,   0,   0,   0,   /* H: H He Ho Hf Hg Hs */
     53,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  49,   0,   0,   0,  77,   0,   0,   0,   0,   0,   0,   0,   0,   /* I: In I Ir */
      0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   /* J: */
     19,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  36,   0,   0,   0,   0,   0,   0,   0,   0,   /* K: K Kr */
      0,  57,   0,   0,   0,   0,   0,   0,   0,   3,   0,   0,   0,   0,   0,   0,   0,   0, 103,   0,   0,  71, 116,   0,   0,   0,   0,   /* L: Li La Lu Lr Lv */
      0,   0,   0, 115, 101,   0,   0,  12,   0,   0,   0,   0,   0,   0,  25,  42,   0,   0,   0,   0, 109,   0,   0,   0,   0,   0,   0,   /* M: Mg Mn Mo Md Mt Mc */
      7,  11,  41,   0,  60,  10,   0,   0, 113,  28,   0,   0,   0,   0,   0, 102,  93,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   /* N: N Ne Na Ni Nb Nd Np No Nh */
      8,   0,   0,   0,   0,   0,   0, 118,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  76,   0,   0,   0,   0,   0,   0,   0,   /* O: O Os Og */
     15,  91,  82,   0,  46,   0,   0,   0,   0,   0,   0,   0,   0,  61,   0,  84,   0,   0,  59,   0,  78,  94,   0,   0,   0,   0,   0,   /* P: P Pd Pr Pm Pt Pb Po Pa Pu */
      0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   /* Q: */
      0,  88,  37,   0,   0,  75, 104, 111,  45,   0,   0,   0,   0,   0,  86,   0,   0,   0,   0,   0,   0,  44,   0,   0,   0,   0,   0,   /* R: Rb Ru Rh Re Rn Ra Rf Rg */
     16,   0,  51,  21,   0,  34,   0, 106,   0,  14,   0,   0,   0,  62,  50,   0,   0,   0,  38,   0,   0,   0,   0,   0,   0,   0,   0,   /* S: Si S Sc Se Sr Sn Sb Sm Sg */
      0,  73,  65,  43,   0,  52,   0,   0,  90,  22,   0,   0,  81,  69,   0,   0,   0,   0,   0, 117,   0,   0,   0,   0,   0,   0,   0,   /* T: Ti Tc Te Tb Tm Ta Tl Th Ts */
     92,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   /* U: U */
     23,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   /* V: V */
     74,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   /* W: W */
      0,   0,   0,   0,   0,  54,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   /* X: Xe */
     39,   0,  70,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   /* Y: Y Yb */
      0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,  30,   0,   0,   0,  40,   0,   0,   0,   0,   0,   0,   0,   0};  /* Z: Zn Zr */
    """
    const uint8_t SMI_ELEMENT[702]


cdef str SMI_UPPER = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
cdef str SMI_LOWER = 'abcdefghijklmnopqrstuvwxyz'


def smi_element_table():
    """Expose SMI_ELEMENT to the test suite as a dict {symbol: atomic number}."""
    cdef uint32_t a, b, z
    cdef dict out = {}
    for a in range(26):
        for b in range(27):
            z = SMI_ELEMENT[a * 27 + b]
            if z:
                if b:
                    out[SMI_UPPER[a] + SMI_LOWER[b - 1]] = z
                else:
                    out[SMI_UPPER[a]] = z
    return out


with cython.warn.undeclared(False):
    # bare so Python can import it, guarded so warn.undeclared stays quiet
    class IncorrectSmiles(ValueError):
        """The string is not a SMILES: the reader could not decide what molecule it names.

        Raised for SYNTAX only -- an unknown element spelling, an unbalanced parenthesis, a ring
        bond that never closes, a character with no meaning at that position, a field the arena
        cannot store.  The message ends in a byte offset into the input.

        NOT raised for chemistry.  A hypervalent atom, an aromatic ring with no Kekule form, an
        element with no aromatic form and a lowercase atom carrying no aromatic bond are all
        stored and reported through the log, because a reader that rejects them is a reader that
        cannot read what other tools emit.
        """


# "no atom" / "no bond" in the parse graph.  A typed global rather than a DEF for the reason
# `_kekule.pxi` gives for AROM_NO_ATOM: it is compared against inside `nogil`, where a literal
# 0xFFFFFFFF would be a Python int.
cdef uint32_t SMI_NONE = 0xFFFFFFFF

# ring-closure labels: `0`-`9` and `%00`-`%99` in slots 0..99, then `%(0)`-`%(99999)` -- the ChemAxon
# spelling for a label above 99 -- in the slots above, one per label LIVE AT ONCE and freed on close.
# Bare `0` is accepted, because OpenSMILES allows it and refusing a label another tool emits buys
# nothing.  A bracketed label is read and never written: the writer numbers its own closures, and 100
# concurrent rings in one string has never been the shape of a real record.
DEF SMI_PLAIN_CLOSURES = 100
DEF SMI_BIG_CLOSURES = 32
DEF SMI_CLOSURES = SMI_PLAIN_CLOSURES + SMI_BIG_CLOSURES

# the largest hydrogen count a bracket may STATE.  One below the arena's nibble maximum, because the
# top value is reserved for the "count unknown" sentinel `H_UNKNOWN`: an atom with no valence rule
# has no count, and a bracket that says 15 is not making that statement -- bounding by the nibble's
# WIDTH would let a file's literal `[CH15]` through as the sentinel, turning a molecule into an
# unanswered question.  No element carries fifteen hydrogens, so this costs no input.
#
# Named for the domain it bounds rather than derived from the storage width next door: the arena
# states the largest real count itself, and a reader that recomputed it from the nibble would go on
# agreeing with it right up to the day the nibble grew.
DEF SMI_H_MAX = H_IMPLICIT_MAX

# what a `@` token said.  UNKNOWN is OpenSMILES `@?`: a centre the writer declared unresolved,
# which is a different fact from a centre nobody mentioned.
cdef enum:
    SMI_CHIRAL_NONE = 0
    SMI_CHIRAL_AT = 1
    SMI_CHIRAL_ATAT = 2
    SMI_CHIRAL_UNKNOWN = 3

# what a `/` or `\` said about the bond, oriented from the atom written FIRST: UP means the second
# atom is up-right of the first, which is the `F/C` reading.
cdef enum:
    SMI_DIR_NONE = 0
    SMI_DIR_UP = 1
    SMI_DIR_DOWN = 2


cdef struct smi_atom_t:
    uint8_t element
    int8_t charge
    uint8_t chiral         # SMI_CHIRAL_*
    bint lower             # written lowercase: the atom-case aromatic statement
    bint radical           # from the CXSMILES tail's `^N:` field
    uint8_t sg_kind        # from the tail's `a:` / `o<n>:` / `&<n>:` fields; 0 is "not named"
    uint16_t sg_group      # as written, unbounded here: `set_stereo_group` owns the range
    uint8_t cip            # a CIP descriptor the tail stated, HELD AS THE ASCII LETTER ITSELF, and
                           # the case is the point: lowercase `r`/`s` are CIP's pseudo-asymmetric
                           # descriptors from the auxiliary rules, a different determination about a
                           # different kind of centre -- not a spelling of `R`/`S`.  Anything that
                           # upper-cases on the way in has lost information.  0 is "the tail named
                           # no descriptor for this atom", which is not the same as "no centre".
    uint16_t isotope
    uint16_t map_number
    int8_t stated_h        # the bracket's H count; AROM_H_UNSTATED for a bare atom
    int8_t implicit_h      # what `smi_read_h` decided
    uint32_t slots         # neighbour slots claimed, the implicit-H slot included
    uint32_t nbr_off       # offset of this atom's slots in `nbrs`
    uint32_t h_pos         # the neighbour position this atom's hydrogens take: 0 when it leads its
                           # component, 1 (just past the bond to the atom before it) otherwise.  A
                           # slot is CLAIMED there only when the bracket stated a count; for a bare
                           # atom nothing holds the position and `smi_written_pair` inserts it.
    bint is_r              # `[R]`: the marker, element 0.  A zeroed struct is not one, hence the flag
    uint8_t r_index        # the digits after `R`, 0 for a bare `[R]`
    uint32_t label_off     # the text a bracket gave where an element symbol belongs -- `[Pol]`,
                           # `[Resin]`, `[REG42]` -- as a byte offset, since the struct holds no
                           # Python object.  Resolved against the source in `smi_parse_build` and
                           # stored as the atom's ALIAS, which is where a label the arena cannot read
                           # as an element belongs.
    uint16_t label_len     # 0 when the atom carries no label
    bint label_tail        # the offset is into the CXSMILES TAIL's bytes and not the body's: the
                           # tail's `$...$` field names labels too, and it is a separate `bytes`.
                           # Which string to slice is the only difference; the alias is the same fact
    uint32_t sid           # stable id in the built molecule


cdef struct smi_bond_t:
    uint32_t u             # the atom written first
    uint32_t v
    uint32_t upos          # u's neighbour slot this bond occupies
    uint32_t vpos
    uint8_t order          # as the string stated it; promotion may raise it to 4
    uint8_t dir_           # SMI_DIR_*, oriented from u


cdef struct smi_parse_t:
    void *block
    smi_atom_t *atoms
    smi_bond_t *bonds
    uint32_t *nbrs         # per slot: the bond index, or SMI_NONE for the implicit-H slot
    uint32_t *stack        # branch stack: the atom each `(` returns to
    uint32_t n_atoms
    uint32_t n_bonds
    uint32_t n_slots
    uint32_t depth
    uint32_t promote_hint  # bonds between two lowercase atoms whose stated order is not 4
    uint32_t chiral_count  # `@` tokens seen, so step 3 can say what it is not applying yet
    # `/` and `\` tokens seen.  COUNTED AT LINK TIME rather than walked for later, because the count
    # decides whether the arena is laid out with a parity segment and that decision precedes the build.
    uint32_t dir_count
    uint32_t dative_arrows # `->` / `<-` tokens: order 8 is stored, the arrow's direction is not
    uint32_t sg_marks      # enhanced stereo group marks the tail placed, so the apply can be skipped
    uint32_t cip_marks     # CIP descriptors the tail stated, likewise
    uint32_t open_atom[SMI_CLOSURES]
    uint32_t open_pos[SMI_CLOSURES]
    uint32_t open_at[SMI_CLOSURES]     # byte offset of the opening label, for the error message
    uint8_t open_order[SMI_CLOSURES]   # 0 when the opening label stated no order
    uint8_t open_dir[SMI_CLOSURES]
    uint8_t open_live[SMI_CLOSURES]
    uint32_t big_label[SMI_BIG_CLOSURES]   # the bracketed label a slot above 99 stands for


cdef int smi_alloc(smi_parse_t *p, uint32_t n) except -1:
    """One struct, one malloc, one check, one free -- RULES.md 5.2, and no realloc at all.

    Every array is bounded by the input length, which is why there is no growth path: an atom
    needs at least one character, so does a bond (a chain bond consumes its atom's character, a
    ring bond consumes two label characters between them), a branch needs its `(`, and the slots
    are two per bond plus at most one implicit hydrogen per atom.
    """
    cdef size_t cap = <size_t> n + 1
    cdef size_t atoms_len = align8(cap * sizeof(smi_atom_t))
    cdef size_t bonds_len = align8(cap * sizeof(smi_bond_t))
    cdef size_t nbrs_len = align8(3 * cap * sizeof(uint32_t))
    cdef size_t stack_len = align8(cap * sizeof(uint32_t))
    cdef size_t total = atoms_len + bonds_len + nbrs_len + stack_len
    cdef char *block = <char *> PyMem_Malloc(total)
    if block is NULL:
        raise MemoryError('SMILES parse allocation failed')
    memset(block, 0, total)
    p.block = <void *> block
    cdef size_t off = 0
    p.atoms = <smi_atom_t *> (block + off); off += atoms_len
    p.bonds = <smi_bond_t *> (block + off); off += bonds_len
    p.nbrs = <uint32_t *> (block + off); off += nbrs_len
    p.stack = <uint32_t *> (block + off)
    p.n_atoms = 0
    p.n_bonds = 0
    p.n_slots = 0
    p.depth = 0
    p.promote_hint = 0
    p.chiral_count = 0
    p.dir_count = 0
    p.dative_arrows = 0
    p.sg_marks = 0
    p.cip_marks = 0
    cdef uint32_t i
    for i in range(SMI_CLOSURES):
        p.open_live[i] = 0
    return 0


cdef void smi_free(smi_parse_t *p) noexcept:
    PyMem_Free(p.block)
    p.block = NULL


# Character classes as byte comparisons rather than a lookup table: the tokeniser tests each
# character against at most a handful of these, and a 256-byte table would be a cache line spent
# to save an integer compare.
cdef inline bint smi_isdigit(char c) noexcept nogil:
    return 48 <= c <= 57


cdef inline bint smi_isupper(char c) noexcept nogil:
    return 65 <= c <= 90


cdef inline bint smi_islower(char c) noexcept nogil:
    return 97 <= c <= 122


cdef inline uint32_t smi_bracket_symbol(const char *s, uint32_t n, uint32_t *ip,
                                        bint *lower) noexcept nogil:
    """The element inside a bracket: greedy two letters, then one.  0 when it is not an element.

    Greedy is right here and only here, because inside a bracket the whole symbol is the element:
    `[Sc]` is scandium, never sulfur next to an aromatic carbon.  The character after the symbol is
    always one of `@ H + - : ]`, none of them lowercase, so the greedy match cannot run past the
    element into another field.
    """
    cdef uint32_t i = ip[0]
    cdef char c0 = s[i]
    cdef char c1
    cdef uint32_t a, z
    if smi_islower(c0):
        lower[0] = True
        a = <uint32_t> (c0 - 97)
    elif smi_isupper(c0):
        lower[0] = False
        a = <uint32_t> (c0 - 65)
    else:
        return 0
    if i + 1 < n:
        c1 = s[i + 1]
        if smi_islower(c1):
            z = SMI_ELEMENT[a * 27 + <uint32_t> (c1 - 97) + 1]
            if z:
                ip[0] = i + 2
                return z
    z = SMI_ELEMENT[a * 27]
    if z:
        ip[0] = i + 1
        return z
    return 0


cdef inline uint32_t smi_bare_symbol(const char *s, uint32_t n, uint32_t *ip,
                                     bint *lower) noexcept nogil:
    """The element of an unbracketed atom, or 0 when this character starts no atom.

    The organic subset only, and NOT the bracket table: outside brackets `SC` is sulfur bonded to
    carbon, so a greedy two-letter match would silently read scandium.  Only `Cl` and `Br` are two
    characters, which is the whole reason the subset exists.
    """
    cdef uint32_t i = ip[0]
    cdef char c = s[i]
    cdef uint32_t z
    lower[0] = False
    if c == 67:                                     # C, Cl
        if i + 1 < n and s[i + 1] == 108:
            ip[0] = i + 2
            return 17
        ip[0] = i + 1
        return 6
    if c == 66:                                     # B, Br
        if i + 1 < n and s[i + 1] == 114:
            ip[0] = i + 2
            return 35
        ip[0] = i + 1
        return 5
    if c == 78:
        z = 7                                       # N
    elif c == 79:
        z = 8                                       # O
    elif c == 80:
        z = 15                                      # P
    elif c == 83:
        z = 16                                      # S
    elif c == 70:
        z = 9                                       # F
    elif c == 73:
        z = 53                                      # I
    elif c == 98:
        lower[0] = True
        z = 5                                       # b
    elif c == 99:
        lower[0] = True
        z = 6                                       # c
    elif c == 110:
        lower[0] = True
        z = 7                                       # n
    elif c == 111:
        lower[0] = True
        z = 8                                       # o
    elif c == 112:
        lower[0] = True
        z = 15                                      # p
    elif c == 115:
        lower[0] = True
        z = 16                                      # s
    else:
        return 0
    ip[0] = i + 1
    return z


cdef inline uint32_t smi_digits(const char *s, uint32_t n, uint32_t *ip,
                                uint32_t limit) noexcept nogil:
    """Read a run of at most `limit` digits.  ip[0] does not move when there are none."""
    cdef uint32_t i = ip[0]
    cdef uint32_t value = 0
    cdef uint32_t count = 0
    while i < n and count < limit and smi_isdigit(s[i]):
        value = value * 10 + <uint32_t> (s[i] - 48)
        i += 1
        count += 1
    ip[0] = i
    return value


cdef int smi_bracket(smi_parse_t *p, const char *s, uint32_t n, uint32_t *ip,
                     smi_atom_t *a, object log) except -1:
    """Parse a bracket atom.  `s[ip[0]]` is the `[`; on return ip[0] is past the `]`.

    The fields after the element are read in a LOOP keyed on the leading character rather than in
    OpenSMILES' fixed order, so `[C-H3]` and `[CH3-]` both read.  Field order is a spelling
    preference and readers in the wild differ on it; a REPEATED field is a contradiction and raises.
    """
    cdef uint32_t start = ip[0]
    cdef uint32_t i = start + 1
    cdef uint32_t sym_at
    cdef uint32_t value, z
    cdef int sign
    cdef char c, cls0, cls1
    cdef str cname
    cdef bint lower = False
    cdef bint have_h = False
    cdef bint have_charge = False
    cdef bint have_map = False
    cdef bint have_chiral = False

    if i >= n:
        raise IncorrectSmiles('unterminated bracket atom at position %d' % start)
    if smi_isdigit(s[i]):
        value = smi_digits(s, n, &i, 5)
        if value < 1 or value > ISOTOPE_MAX:
            raise IncorrectSmiles('isotope %d is outside 1..%d at position %d'
                                  % (value, ISOTOPE_MAX, start + 1))
        a.isotope = <uint16_t> value
        if i >= n:
            raise IncorrectSmiles('unterminated bracket atom at position %d' % start)
    # R IS THE MARKER, and `R` + a lowercase letter is the only other reading: Rb, Ru, Rh, Rn, Re, Ra,
    # Rf, Rg.  The character after a bracket symbol is always one of `@ H + - : ]`, none of them
    # lowercase, so this test cannot mistake an element for a marker.  The digits AFTER `R` are the R
    # index -- an isotope precedes the symbol, so the two runs cannot collide.
    #
    # `*` IS THAT SAME MARKER, index 0, and not a wildcard query.  An atom matching nothing is what a
    # producer of `*` in a stored record means by it -- an attachment point, a polymer end, a resin
    # support -- and the rest of the bracket then reads as it does for `[R]`, so `[*:1]` keeps its map
    # number and `[*+]` its charge.  OpenSMILES gives `*` no index to carry; `[R<n>]` is the spelling
    # for one that matters.
    sym_at = i
    if s[i] == 42 or (s[i] == 82 and not (i + 1 < n and smi_islower(s[i + 1]))):   # '*', 'R'
        a.is_r = True
        z = 0
        if s[i] == 42:
            i += 1
        else:
            i += 1
            value = smi_digits(s, n, &i, 5)
            if value > R_INDEX_MAX:
                raise IncorrectSmiles('R index %d is outside 0..%d at position %d'
                                      % (value, R_INDEX_MAX, start + 1))
            a.r_index = <uint8_t> value
    else:
        z = smi_bracket_symbol(s, n, &i, &lower)
    # A LABEL, NOT AN ELEMENT, AND NOT A REFUSAL EITHER.  `[Pol]`, `[Resin]`, `[Rgp]`, `[REG42]`: a
    # bracket holding a display label, a polymer end or a registry identifier.  It becomes the marker
    # carrying that text as its alias -- element 0, matching nothing -- so the record is read and its
    # own words are kept.  An ABBREVIATION (`OMe`, `CF3`) lands here too and stays a marker until
    # `chython.chemistry` expands it: naming a fragment is chemical knowledge and this reader holds
    # none.
    #
    # Two ways in, and the second is the common one: no symbol matched at all, or a symbol matched and
    # what follows it CANNOT BEGIN A FIELD.  Only `H + - : @ ]` can, so `[Pol]` -- polonium, then an
    # `l` -- is a label and not a polonium with a syntax error after it.  The run is letters and digits
    # from where the symbol began, which leaves `[Pol+]` its charge and `[Pol:1]` its map number.
    if z == 0 and not a.is_r or (i < n and (smi_islower(s[i]) or smi_isdigit(s[i])
                                            or (smi_isupper(s[i]) and s[i] != 72))):
        i = sym_at
        while i < n and (smi_isupper(s[i]) or smi_islower(s[i]) or smi_isdigit(s[i])):
            i += 1
        if i == sym_at:
            raise IncorrectSmiles('unknown element symbol at position %d' % i)
        a.is_r = True
        a.r_index = 0
        a.label_off = sym_at
        a.label_len = <uint16_t> (i - sym_at)
        z = 0
        lower = False
        log.append(mc_record('smiles:label-as-marker', (),
                            'the bracket at position %d names `%s`, which is not an element; it is '
                            'stored as the marker carrying that text as its alias'
                            % (start, s[sym_at:i].decode('ascii', 'replace')), mc_info()))
    a.element = <uint8_t> z
    a.lower = lower
    # a bracket ALWAYS states the hydrogen count: no H token means zero, which is a different fact
    # from the silence of a bare atom.  `[n]` is pyridine and `n` is undecided, and `_kekule.pxi`
    # gets that difference through `stated_h` because it changes the Kekule form.
    a.stated_h = 0

    while i < n:
        c = s[i]
        if c == 93:                                 # ']'
            ip[0] = i + 1
            return 0
        elif c == 72:                               # 'H'
            if have_h:
                raise IncorrectSmiles('hydrogen count given twice at position %d' % i)
            have_h = True
            i += 1
            if i < n and smi_isdigit(s[i]):
                value = smi_digits(s, n, &i, 2)
            else:
                value = 1
            if value > SMI_H_MAX:
                raise IncorrectSmiles('hydrogen count %d is above the storable %d at position %d'
                                      % (value, SMI_H_MAX, i))
            a.stated_h = <int8_t> value
        elif c == 43 or c == 45:                    # '+', '-'
            if have_charge:
                raise IncorrectSmiles('charge given twice at position %d' % i)
            have_charge = True
            sign = 1 if c == 43 else -1
            i += 1
            if i < n and smi_isdigit(s[i]):
                value = smi_digits(s, n, &i, 2)
            else:
                # the repeated-sign spelling: `[Fe++]` is `[Fe+2]`
                value = 1
                while i < n and s[i] == c:
                    value += 1
                    i += 1
            # a charge outside the field is CLAMPED, not grounds for dropping the record.  `[Pt+10]`
            # and `[ZrH8+12]` come from a writer turning every dative contact into a formal charge
            # pair; the record's connectivity is still worth having.  LOST and not REPAIRED, because
            # +8 for a stated +12 is a different species and the line has to say so.
            if sign * <int> value < CHARGE_MIN:
                log.append(mc_record('smiles:charge-clamped', (),
                                     'atom at position %d states charge %d and the field holds %d..%d; '
                                     'it is stored as %d' % (start, sign * <int> value, CHARGE_MIN,
                                                             CHARGE_MAX, CHARGE_MIN), mc_lost()))
                a.charge = CHARGE_MIN
            elif sign * <int> value > CHARGE_MAX:
                log.append(mc_record('smiles:charge-clamped', (),
                                     'atom at position %d states charge %d and the field holds %d..%d; '
                                     'it is stored as %d' % (start, sign * <int> value, CHARGE_MIN,
                                                             CHARGE_MAX, CHARGE_MAX), mc_lost()))
                a.charge = CHARGE_MAX
            else:
                a.charge = <int8_t> (sign * <int> value)
        elif c == 64:                               # '@'
            if have_chiral:
                raise IncorrectSmiles('configuration given twice at position %d' % i)
            have_chiral = True
            i += 1
            if i < n and s[i] == 64:
                a.chiral = SMI_CHIRAL_ATAT
                i += 1
            elif i < n and s[i] == 63:              # '?', a writer's explicit "unresolved"
                a.chiral = SMI_CHIRAL_UNKNOWN
                i += 1
            else:
                a.chiral = SMI_CHIRAL_AT
            # OpenSMILES' extended classes.  `[C@H]` cannot be mistaken for one: a class needs two
            # UPPERCASE letters and the `H` field is followed by `]` or a digit.
            if i + 1 < n and smi_isupper(s[i]) and smi_isupper(s[i + 1]):
                cls0 = s[i]
                cls1 = s[i + 1]
                cname = s[i:i + 2].decode('ascii')
                i += 2
                value = smi_digits(s, n, &i, 2)
                if (cls0 == 84 and cls1 == 72) or (cls0 == 65 and cls1 == 76):   # TH, AL
                    a.chiral = SMI_CHIRAL_AT if value <= 1 else SMI_CHIRAL_ATAT
                else:
                    log.append(mc_record('smiles:chiral-class-unsupported', (),
                                        'atom at position %d: chirality class %s%d is not supported and '
                                        'the configuration is dropped' % (start, cname, value),
                                        mc_lost()))
                    a.chiral = SMI_CHIRAL_NONE
            if a.chiral != SMI_CHIRAL_NONE:
                p.chiral_count += 1
        elif c == 58:                               # ':'
            if have_map:
                raise IncorrectSmiles('atom map given twice at position %d' % i)
            have_map = True
            i += 1
            if i >= n or not smi_isdigit(s[i]):
                raise IncorrectSmiles('atom map `:` with no number at position %d' % i)
            value = smi_digits(s, n, &i, 5)
            if value > MAP_NUMBER_MAX:
                raise IncorrectSmiles('atom map %d is above the storable %d at position %d'
                                      % (value, MAP_NUMBER_MAX, i))
            a.map_number = <uint16_t> value
        else:
            raise IncorrectSmiles('unexpected %r inside a bracket atom at position %d'
                                  % (s[i:i + 1].decode('ascii'), i))
    raise IncorrectSmiles('unterminated bracket atom at position %d' % start)


cdef inline uint32_t smi_new_atom(smi_parse_t *p) noexcept nogil:
    """Take the next parse-graph atom.  The block is zeroed, so only the sentinels are set here."""
    cdef uint32_t idx = p.n_atoms
    cdef smi_atom_t *a = p.atoms + idx
    a.stated_h = AROM_H_UNSTATED
    a.implicit_h = -1
    p.n_atoms = idx + 1
    return idx


cdef inline uint32_t smi_claim(smi_parse_t *p, uint32_t a) noexcept nogil:
    """Reserve the next neighbour slot of atom `a`.

    Slot ORDER is the whole reason this is a separate step from making the bond: `@` is a parity
    over the neighbours in the order the string names them, and a ring-closure label names its
    partner where the LABEL stands, not where the ring closes.  Claiming at label time is what
    makes `[C@H]1...` come out the same as if the partner had been written inline.
    """
    cdef uint32_t pos = p.atoms[a].slots
    p.atoms[a].slots = pos + 1
    return pos


cdef int smi_big_label(smi_parse_t *p, uint32_t value, uint32_t at) except -1:
    """The slot a bracketed ring label `%(value)` occupies: its live slot, or the first free one.

    Slots are per label LIVE AT ONCE and not per label value, so a string may name any label in
    0..99999 and reuse it as often as it likes; only concurrency is bounded.
    """
    cdef uint32_t k, free_slot = SMI_BIG_CLOSURES
    for k in range(SMI_BIG_CLOSURES):
        if p.open_live[SMI_PLAIN_CLOSURES + k]:
            if p.big_label[k] == value:
                return <int> (SMI_PLAIN_CLOSURES + k)
        elif free_slot == SMI_BIG_CLOSURES:
            free_slot = k
    if free_slot == SMI_BIG_CLOSURES:
        raise IncorrectSmiles('more than %d bracketed ring labels are open at once, at position %d'
                              % (SMI_BIG_CLOSURES, at))
    p.big_label[free_slot] = value
    return <int> (SMI_PLAIN_CLOSURES + free_slot)


cdef inline uint32_t smi_label_shown(smi_parse_t *p, uint32_t label) noexcept nogil:
    """The label to name in a message: a slot above 99 stands for a bracketed label, not for itself."""
    if label >= SMI_PLAIN_CLOSURES:
        return p.big_label[label - SMI_PLAIN_CLOSURES]
    return label


cdef inline uint8_t smi_default_order(smi_parse_t *p, uint32_t u, uint32_t v) noexcept nogil:
    """An unwritten bond: aromatic between two lowercase atoms, single otherwise."""
    if p.atoms[u].lower and p.atoms[v].lower:
        return 4
    return 1


cdef inline void smi_link(smi_parse_t *p, uint32_t u, uint32_t v, uint32_t upos,
                          uint8_t order, uint8_t dir_) noexcept nogil:
    cdef smi_bond_t *b = p.bonds + p.n_bonds
    b.u = u
    b.v = v
    b.upos = upos
    b.vpos = smi_claim(p, v)
    b.order = order
    b.dir_ = dir_
    p.n_bonds += 1
    if dir_:
        p.dir_count += 1
    if order != 4 and p.atoms[u].lower and p.atoms[v].lower:
        # a candidate for atom-case promotion.  Counting here, where both ends are known, is what
        # lets the common string skip ring perception entirely: the hint is zero for every string
        # that writes its aromatic bonds aromatic, which is almost all of them.
        p.promote_hint += 1


cdef int smi_tokenize(smi_parse_t *p, const char *s, uint32_t n, object log) except -1:
    """One pass over the bytes, building the parse graph.  Raises `IncorrectSmiles` on syntax."""
    cdef uint32_t i = 0
    cdef uint32_t j
    cdef uint32_t prev = SMI_NONE
    cdef uint32_t a_idx = 0
    cdef uint32_t label, u, value
    cdef uint8_t order, dir_, cdir
    cdef uint8_t pend_order = 0
    cdef uint8_t pend_dir = SMI_DIR_NONE
    cdef uint32_t pend_at = 0                        # where the pending bond token was written
    cdef char c
    cdef bint lower = False
    cdef uint32_t z
    cdef smi_atom_t *a

    # every branch that is not an atom `continue`s; falling out of the chain means `a_idx` is a new
    # atom waiting to be joined to `prev`, which is the one place that joining happens
    while i < n:
        c = s[i]
        if c == 91:                                  # '['
            a_idx = smi_new_atom(p)
            smi_bracket(p, s, n, &i, p.atoms + a_idx, log)
        elif c == 40:                                # '('
            if prev == SMI_NONE:
                raise IncorrectSmiles('branch opens before any atom at position %d' % i)
            if pend_order or pend_dir:
                raise IncorrectSmiles('bond token immediately before `(` at position %d' % i)
            p.stack[p.depth] = prev
            p.depth += 1
            i += 1
            continue
        elif c == 41:                                # ')'
            if not p.depth:
                raise IncorrectSmiles('unbalanced `)` at position %d' % i)
            if pend_order or pend_dir:
                raise IncorrectSmiles('bond token immediately before `)` at position %d' % i)
            p.depth -= 1
            prev = p.stack[p.depth]
            i += 1
            continue
        elif c == 46:                                # '.', a component break
            if pend_order or pend_dir:
                raise IncorrectSmiles('bond token immediately before `.` at position %d' % i)
            prev = SMI_NONE
            i += 1
            continue
        elif c == 45 or c == 61 or c == 35 or c == 58:     # '-', '=', '#', ':'
            if pend_order or pend_dir:
                raise IncorrectSmiles('two bond tokens in a row at position %d' % i)
            if c == 45 and i + 1 < n and s[i + 1] == 62:   # '->', a dative bond
                pend_order = 8
                pend_at = i
                p.dative_arrows += 1
                i += 2
                continue
            if c == 45:
                pend_order = 1
            elif c == 61:
                pend_order = 2
            elif c == 35:
                pend_order = 3
            else:
                pend_order = 4
            pend_at = i
            i += 1
            continue
        elif c == 47 or c == 92:                     # '/', '\'
            if pend_order or pend_dir:
                raise IncorrectSmiles('two bond tokens in a row at position %d' % i)
            # a directional token states a DIRECTION, not an order: `C/C=C/C` is single bonds and
            # `c/c` is still aromatic
            pend_dir = SMI_DIR_UP if c == 47 else SMI_DIR_DOWN
            pend_at = i
            i += 1
            # `C/C=C\\C`: a doubled token is one token.  An unescape the producer lost, and it is the
            # backslash that doubles, so the pair is read as the single direction it stood for rather
            # than as the "two bond tokens in a row" it looks like.  Repaired and logged; a THIRD in a
            # row is still a refusal, since nothing plausible doubles a token twice.
            if i < n and s[i] == c:
                log.append(mc_record('smiles:doubled-bond-direction', (),
                                     'the directional token at position %d is written twice; the pair '
                                     'is read as one' % pend_at, mc_repaired()))
                i += 1
            continue
        elif c == 37 or smi_isdigit(c):              # '%NN' or a single-digit ring label
            if c == 37:
                j = i + 1
                if j < n and s[j] == 40:            # '%(NNNNN)', the bracketed label
                    j += 1
                    value = smi_digits(s, n, &j, 5)
                    if j == i + 2 or j >= n or s[j] != 41:
                        raise IncorrectSmiles('`%%(` needs one to five digits and a `)` at position %d'
                                              % i)
                    j += 1
                    label = <uint32_t> smi_big_label(p, value, i)
                else:
                    label = smi_digits(s, n, &j, 2)
                    if j != i + 3:
                        raise IncorrectSmiles('`%%` needs two digits at position %d' % i)
            else:
                label = <uint32_t> (c - 48)
                j = i + 1
            if prev == SMI_NONE:
                raise IncorrectSmiles('ring bond label before any atom at position %d' % i)
            if not p.open_live[label]:
                p.open_live[label] = 1
                p.open_atom[label] = prev
                p.open_pos[label] = smi_claim(p, prev)
                p.open_at[label] = i
                p.open_order[label] = pend_order
                p.open_dir[label] = pend_dir
            else:
                u = p.open_atom[label]
                if u == prev:
                    raise IncorrectSmiles('ring bond %d closes on its own atom at position %d'
                                          % (smi_label_shown(p, label), i))
                if p.open_order[label] and pend_order:
                    order = p.open_order[label]
                    if order != pend_order:
                        # both ends stated an order and they disagree.  The first statement wins
                        # and the second is reported, because a ring that closes is worth more than
                        # a refusal.
                        log.append(mc_record('smiles:ring-bond-order-conflict', (),
                                            'ring bond %d states order %d where it opens (position %d) and '
                                            'order %d where it closes (position %d); the opening order is '
                                            'kept' % (smi_label_shown(p, label), order,
                                                      p.open_at[label], pend_order, i),
                                            mc_repaired()))
                elif p.open_order[label]:
                    order = p.open_order[label]
                elif pend_order:
                    order = pend_order
                else:
                    order = smi_default_order(p, u, prev)
                # `/` and `\` read from the atom written first.  The closing label speaks from the
                # CLOSING atom, so its direction is inverted before it can be compared with, or
                # substituted for, the one the opening label gave.
                dir_ = p.open_dir[label]
                if pend_dir:
                    cdir = SMI_DIR_DOWN if pend_dir == SMI_DIR_UP else SMI_DIR_UP
                    if not dir_:
                        dir_ = cdir
                    elif dir_ != cdir:
                        log.append(mc_record('smiles:ring-bond-dir-conflict', (),
                                            'ring bond %d states opposite directions at its two ends '
                                            '(positions %d and %d); the opening direction is kept'
                                            % (smi_label_shown(p, label), p.open_at[label], i),
                                            mc_repaired()))
                smi_link(p, u, prev, p.open_pos[label], order, dir_)
                p.open_live[label] = 0
            pend_order = 0
            pend_dir = SMI_DIR_NONE
            i = j
            continue
        elif c == 42:                                # '*'
            # the marker, index 0, exactly as `[*]` -- see `smi_bracket`.  A bare marker states no
            # hydrogen count, so it keeps a bare atom's silence rather than the bracket's zero.
            a_idx = smi_new_atom(p)
            a = p.atoms + a_idx
            a.is_r = True
            i += 1
        elif c == 36:                                # '$'
            raise IncorrectSmiles('the quadruple bond `$` cannot be stored: in this core bond '
                                  'order 4 means aromatic, at position %d' % i)
        elif c == 126:                               # '~', the dative bond in chython's dialect
            if pend_order or pend_dir:
                raise IncorrectSmiles('two bond tokens in a row at position %d' % i)
            # NOT a SMARTS any-bond here.  `_smiles_write.pxi` writes an order-8 bond as `~`, so
            # refusing it means the core cannot read its own output: ammonia-borane goes out as
            # `[BH3]~[NH3]` and came back as an exception.  A reader that rejects its own writer is
            # not permissive-but-honest, it is just broken.
            pend_order = 8
            pend_at = i
            i += 1
            continue
        elif c == 60:                                # '<-', a dative bond written backwards
            if i + 1 >= n or s[i + 1] != 45:
                raise IncorrectSmiles('`<` is only meaningful as the dative bond `<-`, at '
                                      'position %d' % i)
            if pend_order or pend_dir:
                raise IncorrectSmiles('two bond tokens in a row at position %d' % i)
            pend_order = 8
            pend_at = i
            p.dative_arrows += 1
            i += 2
            continue
        elif c == 62:                                # '>'
            # UNREACHABLE FROM `read_smiles`, which now splits an arrow off before tokenizing anything
            # and returns a reaction.  Still reachable, and still a refusal, for a `>` INSIDE one side
            # -- a fourth arrow in a reaction, or one in a SMARTS or a SMIRKS side -- so the message
            # says which reading was tried rather than telling the caller to split the string.
            raise IncorrectSmiles('`>` is a reaction arrow and this is one side of a reaction '
                                  'already, at position %d' % i)
        else:
            j = i
            z = smi_bare_symbol(s, n, &j, &lower)
            if z == 0:
                if s[i] == 82:                       # 'R'
                    raise IncorrectSmiles('a fragment attachment point is `[R]`, or `[R<n>]` when the '
                                          'index matters; a bare `R` outside brackets is rubidium '
                                          'misspelt, at position %d' % i)
                raise IncorrectSmiles('unexpected %r at position %d'
                                      % (s[i:i + 1].decode('ascii'), i))
            a_idx = smi_new_atom(p)
            a = p.atoms + a_idx
            a.element = <uint8_t> z
            a.lower = lower
            i = j

        a = p.atoms + a_idx
        if prev == SMI_NONE:
            if pend_order or pend_dir:
                raise IncorrectSmiles('bond token starts a component at position %d' % pend_at)
            # an atom that leads its component puts its hydrogens in position 0
            a.h_pos = 0
            if a.stated_h >= 1:
                smi_claim(p, a_idx)
        else:
            order = pend_order if pend_order else smi_default_order(p, prev, a_idx)
            smi_link(p, prev, a_idx, smi_claim(p, prev), order, pend_dir)
            # ...otherwise they take the position immediately after the bond to the previous atom,
            # which is what makes `[C@H](F)(Cl)Br` and `F[C@H](Cl)Br` opposite configurations
            # rather than something the stereo pass has to patch up afterwards
            a.h_pos = 1
            if a.stated_h >= 1:
                smi_claim(p, a_idx)
        prev = a_idx
        pend_order = 0
        pend_dir = SMI_DIR_NONE

    if p.depth:
        raise IncorrectSmiles('unbalanced `(`: %d branch(es) never close' % p.depth)
    if pend_order or pend_dir:
        raise IncorrectSmiles('bond token at the end of the string, position %d' % pend_at)
    for label in range(SMI_CLOSURES):
        if p.open_live[label]:
            raise IncorrectSmiles('ring bond %d opens at position %d and never closes'
                                  % (smi_label_shown(p, label), p.open_at[label]))
    if not p.n_atoms:
        raise IncorrectSmiles('no atoms in the string')
    return 0


cdef int smi_flat(smi_parse_t *p) except -1:
    """Lay the neighbour slots out flat, one entry per slot holding its BOND index.

    Bond indices and not neighbour atom indices, because both readers of this array want the bond:
    the hydrogen pass needs its order and the stereo pass needs its direction, and the neighbour is
    one comparison away.  The implicit-hydrogen slot holds SMI_NONE.

    This is also where two bonds between the same pair of atoms are caught.  A per-atom quadratic
    scan is the right shape at these degrees, and it is why the tokeniser keeps no incremental
    adjacency: one positional array, built once, is a smaller thing to be wrong about.
    """
    cdef uint32_t i, k, k2, off = 0
    cdef uint32_t nb, nb2
    cdef smi_atom_t *a
    for i in range(p.n_atoms):
        a = p.atoms + i
        a.nbr_off = off
        off += a.slots
    p.n_slots = off
    for k in range(off):
        p.nbrs[k] = SMI_NONE
    for k in range(p.n_bonds):
        p.nbrs[p.atoms[p.bonds[k].u].nbr_off + p.bonds[k].upos] = k
        p.nbrs[p.atoms[p.bonds[k].v].nbr_off + p.bonds[k].vpos] = k
    for i in range(p.n_atoms):
        a = p.atoms + i
        for k in range(a.nbr_off, a.nbr_off + a.slots):
            if p.nbrs[k] == SMI_NONE:
                continue
            nb = p.bonds[p.nbrs[k]].v if p.bonds[p.nbrs[k]].u == i else p.bonds[p.nbrs[k]].u
            for k2 in range(k + 1, a.nbr_off + a.slots):
                if p.nbrs[k2] == SMI_NONE:
                    continue
                nb2 = (p.bonds[p.nbrs[k2]].v if p.bonds[p.nbrs[k2]].u == i
                       else p.bonds[p.nbrs[k2]].u)
                if nb == nb2:
                    raise IncorrectSmiles('two atoms are bonded twice; the ring bond labels '
                                          'around atom %d contradict each other' % (i + 1))
    return 0


cdef inline uint32_t smi_find_bond(smi_parse_t *p, uint32_t u, uint32_t v) noexcept nogil:
    """The bond between two parse-graph atoms, or SMI_NONE.  A slot scan, so O(degree)."""
    cdef uint32_t k, bi
    cdef smi_atom_t *a = p.atoms + u
    for k in range(a.nbr_off, a.nbr_off + a.slots):
        bi = p.nbrs[k]
        if bi == SMI_NONE:
            continue
        if p.bonds[bi].u == v or p.bonds[bi].v == v:
            return bi
    return SMI_NONE


cdef inline uint32_t smi_cx_field_end(const char *s, uint32_t n, uint32_t i) noexcept nogil:
    """One past the end of the CXSMILES field starting at `i`.

    A field's value may itself contain commas -- `^1:0,2`, `f:0.1,2.3`, `Sg:n:1,2` -- so a comma
    ends a field only when what follows it is not another number.  No CXSMILES field key begins with
    a digit, which is what makes that rule exact rather than a heuristic.

    `$...$` and `(...)` are delimited blocks and are consumed to their closing character instead:
    their contents are free text and coordinates, so a comma inside them means nothing here and a
    `^` inside them is not a radical field.
    """
    cdef uint32_t j = i
    if s[j] == 36:                                   # `$` atom labels and values
        j += 1
        while j < n and s[j] != 36:
            j += 1
        return j + 1 if j < n else n
    if s[j] == 40:                                   # `(` coordinates
        while j < n and s[j] != 41:
            j += 1
        return j + 1 if j < n else n
    while j < n:
        if s[j] == 124:                              # `|` closes the block
            return j
        if s[j] == 44 and (j + 1 >= n or not smi_isdigit(s[j + 1])):
            return j
        j += 1
    return n


cdef inline bint smi_cx_next_index(const char *s, uint32_t j, uint32_t *kp,
                                   uint32_t *out) noexcept nogil:
    """The next comma-separated atom index in `[*kp, j)`, advancing past it and its comma.

    False when what stands at `*kp` is not a digit, which is how both callers tell a finished list
    from a malformed one -- they differ only in the message, so the scanning is here and once.
    """
    cdef uint32_t k = kp[0]
    cdef uint32_t idx = 0
    if k >= j or not smi_isdigit(s[k]):
        return False
    while k < j and smi_isdigit(s[k]):
        idx = idx * 10 + <uint32_t> (s[k] - 48)
        k += 1
    if k < j and s[k] == 44:
        k += 1
    kp[0] = k
    out[0] = idx
    return True


# The three enhanced-stereo kinds `set_stereo_group` accepts, read at import from the names
# `_molecule_container.pxi` publishes rather than respelled as 1, 2 and 3 here (RULES.md 6): the
# domain is declared once, beside the method that validates it, and a C copy that drifted from it
# would be invisible.  Module init runs the fragments in include order and this one is last, so the
# names are bound by the time these lines execute.
cdef uint8_t SMI_SG_ABS = STEREO_ABS
cdef uint8_t SMI_SG_OR = STEREO_OR
cdef uint8_t SMI_SG_AND = STEREO_AND


# ------------------------------------------------------------------------------------------------
# ONE TAIL, MANY COMPONENTS
# ------------------------------------------------------------------------------------------------
#
# A molecule's tail indexes the atoms of one parse; a REACTION's tail indexes the atoms of every
# component of every side, in written order, as one space -- Daylight's reaction-CXSMILES rule.  So
# every function below takes an ARRAY of parse states and a count rather than one state, and resolves
# an index to the component that holds it.  `count` is 1 on the molecule path, where the resolution
# collapses to a single component.
#
# The alternative -- a second router beside `smi_cx`, the way `smk_cx` is a second router beside
# `sma_cx` -- would be a second copy of every message and every field's meaning.  `smk_cx`
# cannot avoid that, because a query side and a patch side do DIFFERENT things with the same
# field; three sides of a reaction all do the same thing, so one router serves them.

cdef inline uint32_t smi_cx_total(smi_parse_t *ps, uint32_t count) noexcept nogil:
    cdef uint32_t c, total = 0
    for c in range(count):
        total += ps[c].n_atoms
    return total


cdef inline uint32_t smi_cx_owner(smi_parse_t *ps, uint32_t count, uint32_t idx,
                                  uint32_t *local) noexcept nogil:
    """The component holding global atom index `idx`, or `count` when no component does."""
    cdef uint32_t c, base = 0
    for c in range(count):
        if idx - base < ps[c].n_atoms:
            local[0] = idx - base
            return c
        base += ps[c].n_atoms
    return count


cdef int smi_cx_group(smi_parse_t *ps, uint32_t count, const char *s, uint32_t k, uint32_t j,
                      uint8_t kind, uint32_t group, bytes block, uint32_t i, object log,
                      str dialect) except -1:
    """One `a:` / `o<n>:` / `&<n>:` field: mark every atom it names with that group.

    The mark is stored in the parse graph and applied after the molecule is built, because
    `set_stereo_group` needs a stable id and the tail is read before there are any.

    Two tails spell these three fields identically and only their separators differ, so this is
    shared and `dialect` is which tail is talking -- a log line about a brace block that called
    itself CXSMILES would send a reader to the wrong half of the string.
    """
    # `idx` is written through a pointer, which Cython cannot see, so its initialiser is here for the
    # same reason `smi_read_h`'s `hn` has one: this file's gate is zero warnings
    cdef uint32_t idx = 0
    cdef uint32_t taken = 0
    cdef uint32_t owner, local = 0
    while k < j:
        if not smi_cx_next_index(s, j, &k, &idx):
            log.append(mc_record('smiles:stereo-group-malformed', (),
                                'the %s stereo group field `%s` is malformed after %d index(es) and '
                                'the rest of it was dropped'
                                % (dialect, block[i:j].decode('ascii', 'replace'), taken),
                                mc_lost()))
            break
        taken += 1
        owner = smi_cx_owner(ps, count, idx, &local)
        if owner == count:
            log.append(mc_record('smiles:stereo-group-bad-index', (),
                                'the %s stereo group field names atom %d, but the string has %d '
                                'atom(s); the mark was dropped'
                                % (dialect, idx, smi_cx_total(ps, count)), mc_lost()))
        elif ps[owner].atoms[local].sg_kind:
            log.append(mc_record('smiles:stereo-group-duplicate', (),
                                'atom %d is put in two enhanced stereo groups by the %s tail; the '
                                'first one is kept' % (idx + 1, dialect), mc_lost()))
        else:
            ps[owner].atoms[local].sg_kind = kind
            ps[owner].atoms[local].sg_group = <uint16_t> group
            ps[owner].sg_marks += 1
    if not taken:
        log.append(mc_record('smiles:stereo-group-empty', (),
                            'the %s stereo group field `%s` names no atom'
                            % (dialect, block[i:j].decode('ascii', 'replace'))))
    return 0


cdef inline uint32_t smi_cx_ref_end(const char *s, uint32_t n, uint32_t i,
                                    uint32_t *value) noexcept nogil:
    """One past the `&#NN;` character reference at `i`, its code point through `value`.

    `i` itself when what stands there is not one: unterminated, empty, not a number, or past the
    Unicode range.  Used by the entry SPLITTER as well as by the decoder, and that is the point -- a
    reference's own `;` is not an entry separator, so `|$a&#59;b;;c$|` is three labels, not five.
    """
    cdef uint32_t k = i + 2, start, v = 0, digit
    cdef bint hexadecimal
    if i + 2 >= n or s[i] != 38 or s[i + 1] != 35:                  # `&#`
        return i
    hexadecimal = s[k] == 120 or s[k] == 88                          # `x`, `X`
    if hexadecimal:
        k += 1
    start = k
    while k < n and s[k] != 59:                                      # `;`
        if 48 <= s[k] <= 57:
            digit = <uint32_t> (s[k] - 48)
        elif hexadecimal and 97 <= s[k] <= 102:
            digit = <uint32_t> (s[k] - 87)
        elif hexadecimal and 65 <= s[k] <= 70:
            digit = <uint32_t> (s[k] - 55)
        else:
            return i
        v = v * (16 if hexadecimal else 10) + digit
        if v > 0x10FFFF:
            return i
        k += 1
    if k >= n or k == start:
        return i
    value[0] = v
    return k + 1


cdef bytes smi_cx_unescape(bytes raw):
    """A CXSMILES label with its `&#NN;` character references resolved.

    That is how ChemAxon spells what a label cannot hold literally -- `;` ends an entry, `$` ends the
    field, `|` ends the block, `&` starts a reference -- and how it spells every non-ASCII character,
    a SMILES being ASCII.  Measured against Marvin 25.1.3: `a;b|c$d` writes as
    `a&#59;b&#124;c&#36;d` and `αβ` as `&#945;&#946;`.  `,`, a space and `%` are written literally.

    Decimal and hexadecimal (`&#x3b1;`), because XML defines both.  A code point above 127 is stored
    as its UTF-8 bytes, an alias being bytes.  A reference that is unterminated, empty or not a number
    STAYS LITERAL: an alias is display text and a `&#` in it may be the text.
    """
    if b'&#' not in raw:
        return raw
    cdef const char *s = <const char *> raw
    cdef bytearray out = bytearray()
    cdef uint32_t n = <uint32_t> len(raw)
    cdef uint32_t i = 0, k
    cdef uint32_t value = 0
    while i < n:
        k = smi_cx_ref_end(s, n, i, &value)
        if k == i:
            out.append(raw[i])
            i += 1
        else:
            out += chr(value).encode('utf8')
            i = k
    return bytes(out)


cdef int smi_cx_labels(smi_parse_t *ps, uint32_t count, const char *s, uint32_t i, uint32_t j,
                       bytes block, object log) except -1:
    """The tail's `$...$` field: one label per atom, `;`-separated, in the tail's index space.

    A label is DISPLAY TEXT and becomes the atom's alias, whatever the atom is -- Marvin writes
    `CCC |$;;OMe$|`, a carbon carrying `OMe`, and `[OMe]C` reads as the marker carrying it, so the
    element is the file's statement and the alias is the label either way.  A tail label OUTRANKS a
    bracket label on the same atom: the tail is written second and by the same producer.

    Three entries are RESERVED and are not text.  `_R<n>` is ChemAxon's R-group spelling -- the atom
    is `*` in the body and its index is here -- and lands in `r_index`.  `_AP<n>` marks an attachment
    point, whose ordinal nothing in this arena holds; the marker itself is already the attachment.
    `star_e` says the star is a star, which the body already said.

    `$_AV:...$` is the atom-VALUES field and not this one: same delimiters, different meaning, so it
    is named in the log rather than read as a hundred labels.
    """
    cdef bytes content
    cdef uint32_t owner, local = 0, idx = 0
    cdef uint32_t off, end, ref, stop = j - 1
    cdef uint32_t value = 0
    cdef smi_atom_t *a
    if stop <= i or s[stop] != 36:
        log.append(mc_record('smiles:cx-labels-unterminated', (),
                             'the CXSMILES atom-label field is not terminated by `$` and was dropped',
                             mc_lost()))
        return 0
    if block[i + 1:i + 5] == b'_AV:':
        log.append(mc_record('smiles:cx-field-not-applied', (),
                             'the CXSMILES field `$_AV:...$` is not applied', mc_lost()))
        return 0
    off = i + 1
    while off <= stop:
        end = off
        while end < stop and s[end] != 59:                          # `;`
            ref = smi_cx_ref_end(s, stop, end, &value)              # never over a reference's own `;`
            end = ref if ref != end else end + 1
        if end > off:
            content = block[off:end]
            owner = smi_cx_owner(ps, count, idx, &local)
            if owner == count:
                log.append(mc_record('smiles:cx-label-bad-index', (),
                                     'the CXSMILES atom-label field gives a label for atom %d, but the '
                                     'string has %d atom(s); it was dropped'
                                     % (idx, smi_cx_total(ps, count)), mc_lost()))
            else:
                a = ps[owner].atoms + local
                if content == b'star_e':
                    log.append(mc_record('smiles:cx-label-star', (),
                                         'atom %d carries the CXSMILES label `star_e`, which says it is '
                                         'the star the string already writes; nothing is stored for it'
                                         % (idx + 1), mc_info()))
                elif content.startswith(b'_AP') and content[3:].isdigit():
                    log.append(mc_record('smiles:cx-label-attachment-point', (),
                                         'atom %d is attachment point %s; the marker is stored and the '
                                         'ordinal is not, nothing in this arena holding one'
                                         % (idx + 1, content[3:].decode('ascii')), mc_info()))
                elif content.startswith(b'_R') and content[2:].isdigit():
                    smi_cx_r_index(a, idx, int(content[2:]), log)
                elif end - off > 0xFFFF:
                    log.append(mc_record('smiles:cx-label-too-long', (),
                                         'the CXSMILES label of atom %d is %d bytes long and was dropped'
                                         % (idx + 1, end - off), mc_lost()))
                else:
                    if a.label_len:
                        log.append(mc_record('smiles:cx-label-outranks-bracket', (),
                                             'atom %d carries a label in its bracket and another in the '
                                             'CXSMILES tail; the tail\'s `%s` is stored'
                                             % (idx + 1, content.decode('ascii', 'replace')), mc_info()))
                    a.label_off = off
                    a.label_len = <uint16_t> (end - off)
                    a.label_tail = True
        idx += 1
        off = end + 1
    return 0


cdef int smi_cx_r_index(smi_atom_t *a, uint32_t idx, uint32_t value, object log) except -1:
    """One `_R<n>` label: the R index ChemAxon writes in the tail for a `*` in the body."""
    if not a.is_r:
        log.append(mc_record('smiles:cx-label-r-on-element', (),
                             'the CXSMILES tail gives atom %d the R index %d, but the string writes it as '
                             'an element; the index was dropped' % (idx + 1, value), mc_lost()))
    elif value > R_INDEX_MAX:
        log.append(mc_record('smiles:cx-label-r-index-too-wide', (),
                             'the CXSMILES tail gives atom %d the R index %d, past the %d this arena '
                             'holds; the marker is stored without an index'
                             % (idx + 1, value, R_INDEX_MAX), mc_lost()))
    elif a.r_index and a.r_index != <uint8_t> value:
        log.append(mc_record('smiles:cx-label-r-index-conflict', (),
                             'atom %d is written `[R%d]` and the CXSMILES tail calls it R%d; the tail\'s '
                             'index is stored' % (idx + 1, a.r_index, value), mc_repaired()))
        a.r_index = <uint8_t> value
    else:
        a.r_index = <uint8_t> value
    return 0


cdef str smi_cx_name(bytes block, uint32_t i, uint32_t j):
    """A short name for the field at `[i, j)`, for a log line a human is going to read.

    The key, which is what identifies the field, plus the value only when the value is short enough
    to be worth reading: a coordinate block for a fifty-atom molecule is six hundred bytes of digits
    and putting them in a log line hides the four other fields around it.
    """
    cdef const char *s = <const char *> block
    cdef uint32_t k = i
    if s[i] == 36:
        return '$...$'
    if s[i] == 40:
        return '(...)'
    if j - i <= 32:
        return block[i:j].decode('ascii', 'replace')
    while k < j and s[k] != 58:
        k += 1
    return block[i:k + 1].decode('ascii', 'replace') + '...'


cdef int smi_cx(smi_parse_t *ps, uint32_t count, bytes block, object log,
                bint f_handled) except -1:
    """Apply the CXSMILES tail.  Radicals, stereo groups and atom labels are applied; every other
    field is named in the log.

    Runs BEFORE `smi_read_h`, because a radical changes the hydrogen count of a bare atom and the
    count is an argument to `add_atom`.

    A radical costs its atom one unit of valence, exactly like a bond: `CC |^1:0|` is the ethyl
    radical and reads as CH2(.)-CH3, not as ethane wearing a flag.  That is measured and not
    assumed -- `valence_implicit_h(z, 0, True, k)` from the chemistry collection equals
    `valence_implicit_h(z, 0, False, k + 1)` for C, N, O and S, and RDKit reads these same strings
    the same way -- so charging one extra order in the NOTATION model reproduces the chemistry model
    without either file learning about the other.

    The arena stores one radical bit per atom, so `^2` through `^7` -- a carbene, a nitrene, a
    trivalent radical -- cannot be stored as what they are.  They become monoradicals and the
    narrowing is logged, because the alternative (charging two units of valence against a one-bit
    flag) returns a molecule whose own hydrogen count contradicts its own radical state, and nothing
    downstream could tell that from a real monoradical.  The writer emits only `^1:`, so a
    chython-written string round-trips exactly.

    Nothing in the tail RAISES.  A bad index or a malformed field is reported and dropped: the
    molecule in front of the tail is intact and refusing it would be the larger loss.  This is the
    one place the file's "syntax raises" rule does not reach, and the reason is that the tail is an
    extension whose failure does not make the SMILES unreadable.

    `f_handled` says the `f:` field has already been read and acted on, which is true only on the
    reaction path: `f:` names COMPONENTS rather than atoms, and a reaction has to know its component
    grouping before it tokenizes anything, so `smi_cx_fgroups` reads that one field first.  On the
    molecule path `f:` still lands in the log, because a molecule has no components to regroup.
    """
    cdef const char *s = <const char *> block
    cdef uint32_t n = <uint32_t> len(block)
    cdef uint32_t i = 1                              # past the opening `|`
    cdef uint32_t j, k, mult, taken
    cdef uint32_t idx = 0                # written through a pointer; see `smi_cx_group`
    cdef uint32_t owner, local = 0
    cdef uint32_t r_at = 0
    cdef bint have_r = False
    cdef bint have_group = False
    cdef char c
    while i + 1 < n:                                 # `n - 1` is the closing `|`
        c = s[i]
        if c == 44 or c == 32:                       # a separator, or space inside a coordinate list
            i += 1
            continue
        j = smi_cx_field_end(s, n, i)
        if j <= i:                                   # cannot happen; a field is at least one byte
            break
        if c == 94:                                  # `^N:` radicals
            k = i + 1
            if k + 1 < j and smi_isdigit(s[k]) and s[k + 1] == 58:
                mult = <uint32_t> (s[k] - 48)
                k += 2
                taken = 0
                while k < j:
                    if not smi_cx_next_index(s, j, &k, &idx):
                        log.append(mc_record('smiles:radical-field-malformed', (),
                                            'the CXSMILES radical field `%s` is malformed after %d index(es) '
                                            'and the rest of it was dropped'
                                            % (block[i:j].decode('ascii', 'replace'), taken),
                                            mc_lost()))
                        break
                    taken += 1
                    owner = smi_cx_owner(ps, count, idx, &local)
                    if owner == count:
                        log.append(mc_record('smiles:radical-bad-index', (),
                                            'the CXSMILES radical field names atom %d, but the string has '
                                            '%d atom(s); the mark was dropped'
                                            % (idx, smi_cx_total(ps, count)), mc_lost()))
                    elif ps[owner].atoms[local].radical:
                        log.append(mc_record('smiles:radical-duplicate', (),
                                            'atom %d is marked a radical twice in the CXSMILES tail; it is '
                                            'stored as one radical' % (idx + 1)))
                    else:
                        ps[owner].atoms[local].radical = True
                if not taken:
                    log.append(mc_record('smiles:radical-field-empty', (),
                                        'the CXSMILES radical field `%s` names no atom'
                                        % block[i:j].decode('ascii', 'replace')))
                elif mult == 0 or mult > 7:
                    # 1 to 7 is the whole defined range, so this is a spelling nobody meant; it
                    # still says "radical", which is the part this arena can store
                    log.append(mc_record('smiles:radical-class-out-of-range', (),
                                        'the CXSMILES radical class `^%d:` is not one of 1 to 7; the %d '
                                        'atom(s) it names are stored as radicals anyway' % (mult, taken),
                                        mc_repaired()))
                elif mult != 1:
                    # the class numbers name a carbene, a nitrene and the trivalent radicals; their
                    # electron counts are not the class number and this line does not claim they are
                    log.append(mc_record('smiles:radical-class-narrowed', (),
                                        'the CXSMILES radical class `^%d:` is not a monovalent radical; this '
                                        'arena stores one radical bit per atom, so the %d atom(s) it names '
                                        'are stored as monoradicals' % (mult, taken), mc_repaired()))
            else:
                log.append(mc_record('smiles:radical-field-malformed', (),
                                    'the CXSMILES radical field `%s` is malformed and was dropped'
                                    % block[i:j].decode('ascii', 'replace'), mc_lost()))
        elif c == 36:                                             # `$...$` atom labels
            smi_cx_labels(ps, count, s, i, j, block, log)
        elif c == 97 and i + 1 < j and s[i + 1] == 58:            # `a:` absolute
            smi_cx_group(ps, count, s, i + 2, j, SMI_SG_ABS, 0, block, i, log, 'CXSMILES')
        elif f_handled and c == 102 and i + 1 < j and s[i + 1] == 58:   # `f:`, read already
            pass
        elif j == i + 1 and c == 114:                             # `r`, the relative-stereo flag
            # a bare flag with no atom list, so what it means depends on the rest of the block: beside
            # an `&`/`o` group it restates what the group already says and the output is byte-identical
            # without it, alone it is the record's ONLY statement that the centres are relative and
            # they are stored absolute.  Decided after the loop, since the group may be written second.
            have_r = True
            r_at = i
        elif c == 111 or c == 38:                                 # `o<n>:` OR, `&<n>:` AND
            k = i + 1
            mult = smi_digits(s, j, &k, 4)
            if k == i + 1 or k >= j or s[k] != 58:
                log.append(mc_record('smiles:cx-field-not-applied', (),
                                     'the CXSMILES field `%s` is not applied' % smi_cx_name(block, i, j),
                                     mc_lost()))
            else:
                have_group = True
                # NOT bounded here.  Whether a group number is storable is `set_stereo_group`'s
                # question and it is answered there, once, in the class that owns the field -- so an
                # out-of-range number travels to the apply and is reported with that class's own
                # words rather than measured against a second copy of the range in this file.
                smi_cx_group(ps, count, s, k + 1, j, SMI_SG_OR if c == 111 else SMI_SG_AND, mult,
                             block, i, log, 'CXSMILES')
        else:
            log.append(mc_record('smiles:cx-field-not-applied', (),
                                 'the CXSMILES field `%s` is not applied' % smi_cx_name(block, i, j),
                                 mc_lost()))
        i = j
    if have_r:
        if have_group:
            log.append(mc_record('smiles:cx-relative-flag-redundant', (),
                                 'the CXSMILES `r` flag at position %d restates the enhanced stereo '
                                 'group in the same block and nothing is stored for it' % r_at,
                                 mc_info()))
        else:
            log.append(mc_record('smiles:cx-relative-flag-unbacked', (),
                                 'the CXSMILES `r` flag at position %d names no enhanced stereo group, '
                                 'so the configurations are stored absolute' % r_at, mc_lost()))
    return 0


cdef inline uint32_t smi_brace_field_end(const char *s, uint32_t n, uint32_t i) noexcept nogil:
    """One past the end of the brace-block field starting at `i`.

    A brace block separates its fields with `;` where CXSMILES uses `,`, but a field's own value
    still spells its atom list with commas -- `o1:19,22` -- so BOTH are terminators here and the
    comma is one only when what follows it is not another digit, for the same reason and by the same
    rule as `smi_cx_field_end`: no field key begins with a digit.  Taking the comma as well costs
    nothing and reads a block that mixes the two separators, which is what a converter between the
    dialects produces.

    There is no `$...$` or `(...)` case, unlike the CXSMILES scanner: those are that tail's delimited
    blocks and nothing here has been seen to spell them.  A block that does gets its field named in
    the log rather than guessed at.
    """
    cdef uint32_t j = i
    while j < n:
        if s[j] == 125 or s[j] == 59:                # `}` closes the block, `;` separates fields
            return j
        if s[j] == 44 and (j + 1 >= n or not smi_isdigit(s[j + 1])):
            return j
        j += 1
    return n


cdef int smi_brace_cip(smi_parse_t *p, const char *s, uint32_t i, uint32_t j, bytes block,
                       object log) except -1:
    """One `A<index>=<letter>` field: record `letter` as that atom's CIP descriptor.

    The letter travels into the parse graph VERBATIM, case included, and WHICH letters name a
    descriptor is not decided here: the arena owns that domain and refuses what it cannot store, so
    the apply reports a refusal in the words of the class that refused it.  Same reason the group
    numbers above travel unbounded -- the range is stated once, beside the field it bounds
    (RULES.md 6).  It is also why nothing in this function upper-cases: `r` and `R` are different
    determinations about different kinds of centre, so a reader that normalised the case would be
    destroying the distinction on the way in.

    What IS decided here is the field's SHAPE -- an index, an `=`, exactly one ASCII letter -- because
    that is this file's question.  `A19=rs` and `A19=` are malformed fields rather than unknown
    descriptors, and the log line that says so names the right problem.
    """
    cdef uint32_t k = i + 1
    cdef uint32_t idx = smi_digits(s, j, &k, 9)
    cdef char c
    if k == i + 1 or k + 2 != j or s[k] != 61:       # no index, no `=`, or not one letter after it
        log.append(mc_record('smiles:cip-field-malformed', (),
                            'the `{...}` CIP field `%s` is malformed and was dropped'
                            % block[i:j].decode('ascii', 'replace'), mc_lost()))
        return 0
    c = s[k + 1]
    # `A<i>=q` is not a determination.  ChemAxon writes it for a centre whose descriptor it did not
    # compute, so there is nothing to store and nothing lost by not storing it: INFO, not LOST, and
    # never `set_atom_cip('q')`, which would report ~35,000 refusals for input that stated no fact.
    if c == 113:                                     # 'q'
        log.append(mc_record('smiles:cip-undetermined', (),
                            'the `{...}` CIP field `%s` states no determination and nothing is stored'
                            % block[i:j].decode('ascii', 'replace'), mc_info()))
    elif not smi_isupper(c) and not smi_islower(c):
        log.append(mc_record('smiles:cip-field-bad-descriptor', (),
                            'the `{...}` CIP field `%s` does not name a descriptor and was dropped'
                            % block[i:j].decode('ascii', 'replace'), mc_lost()))
    elif idx >= p.n_atoms:
        log.append(mc_record('smiles:cip-bad-index', (),
                            'the `{...}` CIP field names atom %d, but the string has %d atom(s); the '
                            'descriptor was dropped' % (idx, p.n_atoms), mc_lost()))
    elif p.atoms[idx].cip:
        log.append(mc_record('smiles:cip-duplicate', (),
                            'atom %d is given two CIP descriptors by the `{...}` block; the first one is kept'
                            % (idx + 1), mc_lost()))
    else:
        p.atoms[idx].cip = <uint8_t> c
        p.cip_marks += 1
    return 0


cdef int smi_brace(smi_parse_t *p, bytes block, object log) except -1:
    """Apply a brace-delimited extension block: enhanced stereo groups and CIP descriptors.

    Some tools write the tail in braces instead of pipes.  It is the same idea as CXSMILES and its
    three enhanced-stereo fields are spelled character for character the same -- `a:`, `o<n>:`,
    `&<n>:`, atom indices ZERO-BASED -- so they route into `smi_cx_group` rather than into a second
    copy of the same scan.  Only the separator differs, `;` for `,`, and one field has no CXSMILES
    equivalent at all: `A<index>=<letter>`, an atom's CIP descriptor stated by whatever assigned it.

    That descriptor is worth reading rather than dropping because it is NOT derivable from the string
    around it.  Enhanced stereo groups say a centre's configuration is one of a set; a descriptor is
    somebody's completed determination about a specific centre, and for the `r`/`s` pseudo-asymmetric
    cases it is a determination this library cannot yet make for itself.  Reading it keeps the fact
    the input carried instead of asking the input to be re-derived from a form that lost it.

    Nothing here RAISES, for the reason the CXSMILES scanner does not either: the tail is an
    extension whose failure does not make the SMILES in front of it unreadable.  Every unrecognised
    or malformed field is named in the log and dropped.  There is no writer for this dialect -- the
    library emits CXSMILES -- so this is a read path only and no round-trip is claimed for it.
    """
    cdef const char *s = <const char *> block
    cdef uint32_t n = <uint32_t> len(block)
    cdef uint32_t i = 1                              # past the opening `{`
    cdef uint32_t j, k, mult
    cdef char c
    while i + 1 < n:                                 # `n - 1` is the closing `}`
        c = s[i]
        if c == 59 or c == 44 or c == 32:            # a separator, either dialect's, or a stray space
            i += 1
            continue
        j = smi_brace_field_end(s, n, i)
        if j <= i:                                   # cannot happen; a field is at least one byte
            break
        if c == 65:                                                # `A<index>=<letter>`
            smi_brace_cip(p, s, i, j, block, log)
        elif c == 97 and i + 1 < j and s[i + 1] == 58:             # `a:` absolute
            smi_cx_group(p, 1, s, i + 2, j, SMI_SG_ABS, 0, block, i, log, '`{...}`')
        elif c == 111 or c == 38:                                  # `o<n>:` OR, `&<n>:` AND
            k = i + 1
            mult = smi_digits(s, j, &k, 4)
            if k == i + 1 or k >= j or s[k] != 58:
                log.append(mc_record('smiles:brace-field-not-applied', (),
                                     'the `{...}` field `%s` is not applied' % smi_cx_name(block, i, j),
                                     mc_lost()))
            else:
                smi_cx_group(p, 1, s, k + 1, j, SMI_SG_OR if c == 111 else SMI_SG_AND, mult,
                             block, i, log, '`{...}`')
        else:
            log.append(mc_record('smiles:brace-field-not-applied', (),
                                 'the `{...}` field `%s` is not applied' % smi_cx_name(block, i, j),
                                 mc_lost()))
        i = j
    return 0


cdef int smi_read_h(smi_parse_t *p, object log) except -1:
    """Fill `implicit_h` for every atom.  Runs after the CXSMILES pass, which can set radicals.

    A bracket atom is already answered: the bracket states the count, zero included.

    A bare atom is the notation model's question, and for an AROMATIC bare atom it is one question
    and not two.  The bond-order sum a reader must charge the atom is its non-aromatic orders, plus
    one per aromatic bond, plus ONE MORE if this atom takes a ring double bond -- and whether it
    does is exactly what `arom_classify_atom` decides.  Thiophene's sulfur is why the sum cannot be
    guessed and then repaired: charge it an extra order and S{2,4,6} lands on 4 and gives it a
    hydrogen it does not have, and since that "found a valence" no second attempt would ever run.
    Asking the classifier first gives 0 for thiophene S, 0 for furan O, 0 for pyridine N, 1 for
    benzene C and 0 for a fusion carbon, from one table and one rule.

    The same table therefore decides the hydrogen count and, later, the Kekule form -- so the two
    cannot disagree about which atoms were spoken for.

    THE CLASSIFICATION HALF NOW LIVES IN `_hydrogens.pxi` as `hyd_arom_takes`, and this function asks
    it rather than restating it.  That half is what a Python caller could not reach -- which is how
    `chython.calc_implicit` came to be a weaker second answer to the same question, refusing on every
    aromatic atom because a `cdef` classifier is invisible from Python.

    THE VALENCE LOOKUP STAYS HERE, on `smv_default_h`, and that is deliberate rather than left over.
    The general derivation charges `takes` to the CHEMISTRY collection, because an MDL record may hold
    a charged aromatic or an element outside the organic subset.  A BARE SMILES atom can be neither:
    it is one of eleven elements by construction, and what its hydrogen count is was settled by
    OpenSMILES, not by chemistry.  Merging the two models is the thing decision #2 at the top of this
    file forbids, so the shared piece is the classification and the models stay apart.

    `count_stated=True` says the language always states the count -- a bare `n` means zero hydrogens
    by OpenSMILES' own rules -- so no SMILES atom is ever ambiguous in the pyrrole-versus-pyridine
    sense and the shared ambiguity gate stays out of the way.
    """
    cdef uint32_t i, k, bi, arom, osum, deg
    # initialised although `smv_default_h` writes it whenever it returns True: Cython cannot see
    # through a pointer out-parameter, so without this the maybe-uninitialized warning is a false
    # positive, and this file's gate is zero warnings.  Same reason as `_smiles_write.pxi:236`.
    cdef uint32_t hn = 0
    cdef bint exo
    cdef uint8_t o, takes = 0, flag = 0, rad
    cdef smi_atom_t *a
    for i in range(p.n_atoms):
        a = p.atoms + i
        if a.stated_h >= 0:
            a.implicit_h = a.stated_h
            continue
        arom = 0
        osum = 0
        deg = 0
        exo = False
        # an unpaired electron costs one unit of valence, exactly like a bond -- see `smi_cx`, which
        # is where that is measured against the chemistry collection rather than asserted
        rad = 1 if a.radical else 0
        for k in range(a.nbr_off, a.nbr_off + a.slots):
            bi = p.nbrs[k]
            if bi == SMI_NONE:                       # the implicit-hydrogen slot
                continue
            o = p.bonds[bi].order
            if o == 8:
                # a dative contact carries no electron pair, so it is not in the bond-order sum and
                # not a neighbour for the aromatic classifier -- `_valence.pxi` says the same and
                # `arom_classify_atom`'s `nbrs` is documented as "every bond of any order except 8".
                # It still holds its neighbour slot, because `@` counts positions, not electrons.
                continue
            deg += 1
            if o == 4:
                arom += 1
            else:
                osum += o
                if o >= 2:
                    exo = True
        if not arom:
            if a.lower:
                # every bond of a lowercase atom was written non-aromatic, so promotion did not
                # reach it either.  Read as written and say so: this is not a case where guessing
                # is better than reporting.
                log.append(mc_record('smiles:lowercase-no-aromatic-bond', (),
                                    'atom %d is written lowercase but carries no aromatic bond; its bonds '
                                    'are stored as the string wrote them' % (i + 1)))
            if smv_default_h(a.element, osum + rad, &hn):
                a.implicit_h = <int8_t> hn
            else:
                # H_UNKNOWN AND NOT 0.  An unbracketed atom does not state its count -- the notation
                # says "the valence model knows this one" -- so when no rule answers, the count is
                # genuinely unknown and 0 would be this reader inventing the answer "none".  A
                # consumer cannot tell an invented 0 from a real one; it can tell None.
                a.implicit_h = H_UNKNOWN
                log.append(mc_record('smiles:no-valence-rule', (),
                                    'atom %d: no valence rule for an unbracketed atom with bond order sum '
                                    '%d, so its hydrogen count is stored as unknown' % (i + 1, osum + rad),
                                    mc_lost()))
            continue
        # `count_stated=True`: see the docstring.  It withdraws the ambiguity gate, which is correct
        # here and only here -- every other format has to leave a bare aromatic nitrogen unknown.
        # `count_stated=True` with no number to go with it: SMILES fixes the count by its own rules
        # (a bare `n` has none, `[nH]` says otherwise) but that is the answer this function is on its
        # way to computing, so there is nothing to hand the classifier yet.
        hyd_arom_takes(a.element, a.charge, a.radical, deg, exo, True, AROM_H_UNSTATED,
                       &takes, &flag)
        if flag & HYD_NO_AROMATIC_FORM:
            log.append(mc_record('smiles:no-aromatic-form', (),
                                'atom %d: this element in this state has no aromatic form; it is stored as '
                                'written and takes no ring double bond' % (i + 1)))
        if smv_default_h(a.element, osum + arom + takes + rad, &hn):
            a.implicit_h = <int8_t> hn
        elif smv_default_h(a.element, osum + arom + (1 - takes) + rad, &hn):
            # the classification the ring implies has no valence, the other one does.  Report it:
            # the stored hydrogen count then disagrees with the class the kekuliser will pick.
            a.implicit_h = <int8_t> hn
            log.append(mc_record('smiles:aromatic-valence-fallback', (),
                                'atom %d: no valence rule for the aromatic reading with bond order sum %d, '
                                'so the other reading was used and gives %d hydrogen(s)'
                                % (i + 1, osum + arom + takes + rad, hn), mc_repaired()))
        else:
            # neither reading has a rule, so nothing here can answer the question either -- same
            # sentinel, same reason as the non-aromatic branch above
            a.implicit_h = H_UNKNOWN
            log.append(mc_record('smiles:no-valence-rule', (),
                                'atom %d: no valence rule for an unbracketed aromatic atom with bond order '
                                'sum %d, so its hydrogen count is stored as unknown'
                                % (i + 1, osum + arom + takes + rad), mc_lost()))
    return 0


cdef MoleculeContainer smi_build(smi_parse_t *p):
    """One edit scope, one arena build.  Atoms in string order, so a stable id names its token."""
    cdef MoleculeContainer mol = MoleculeContainer()
    cdef uint32_t i
    cdef smi_atom_t *a
    cdef smi_bond_t *b
    with mol.edit():
        if p.chiral_count or p.dir_count:
            # `smi_stereo` writes parities into the SEALED arena -- a parity is read against the frame
            # perception derives, so the frame must exist first -- and a persistent segment is laid out
            # once.  So the segment is asked for HERE, where the layout is still open, on the same
            # condition `smi_parse_build` uses to decide whether to call `smi_stereo` at all.
            mol.request_parity()
        for i in range(p.n_atoms):
            a = p.atoms + i
            a.sid = mol.add_atom(<int> a.element, charge=a.charge, isotope=a.isotope,
                                 radical=a.radical, map_number=a.map_number,
                                 implicit_h=a.implicit_h)
            if a.r_index:
                mol.set_r_index(a.sid, a.r_index)
        for i in range(p.n_bonds):
            b = p.bonds + i
            mol.add_bond(p.atoms[b.u].sid, p.atoms[b.v].sid, b.order)
    return mol


cdef int smi_promote(smi_parse_t *p, MoleculeContainer mol, object log) except -1:
    """Atom-case aromatic promotion.  Returns 1 when the molecule was changed, 0 when it was not.

    Reached only when the tokeniser saw a bond between two lowercase atoms whose stated order was
    not aromatic, which is why the ordinary string never pays for this at all.  Biphenyl DOES reach
    it -- `c1ccc(-c2ccccc2)cc1` is what every writer emits -- and pays ring perception and nothing
    else: the inter-ring bond lies in no smallest ring, so no ring here is all-lowercase-with-a-
    non-aromatic-bond and the function returns 0 having touched neither the graph nor the molecule.

    Ring perception needs a built molecule, and hydrogen counts are arguments to `add_atom`, so the
    order is forced: build once with the stated orders, then repair in a second edit scope.  The
    repair rewrites the parse graph first and re-derives EVERY hydrogen count from it, rather than
    patching the atoms around the promoted bonds, because a promoted bond changes the aromatic
    count of its two atoms and therefore their classification.
    """
    cdef list rings = mol.rings
    if not rings:
        return 0
    cdef dict index_of = {}
    cdef uint32_t i, u, v, bi, count
    for i in range(p.n_atoms):
        index_of[p.atoms[i].sid] = i
    cdef set targets = set()
    cdef list plog = []
    cdef list names, members
    cdef tuple ring
    cdef bint all_lower
    for ring in rings:
        all_lower = True
        members = []
        for i in range(<uint32_t> len(ring)):
            u = <uint32_t> index_of[ring[i]]
            members.append(u + 1)
            if not p.atoms[u].lower:
                all_lower = False
                break
        if not all_lower:
            continue
        names = []
        count = <uint32_t> len(ring)
        for i in range(count):
            u = <uint32_t> index_of[ring[i]]
            v = <uint32_t> index_of[ring[(i + 1) % count]]
            bi = smi_find_bond(p, u, v)
            if bi == SMI_NONE:
                raise AssertionError('ring bond %d-%d is not in the parse graph' % (u, v))
            if p.bonds[bi].order != 4:
                targets.add(bi)
                names.append('%d-%d' % (u + 1, v + 1))
        if names:
            plog.append(mc_record('smiles:aromatic-promoted', (),
                                  'every atom of ring %s is written lowercase, so its bond(s) %s are stored '
                                  'aromatic although the string wrote them otherwise'
                                  % (tuple(members), ', '.join(names)), mc_repaired()))
    if not targets:
        return 0

    cdef list old_orders = []
    cdef list old_h = []
    cdef tuple item
    for bi in sorted(targets):
        old_orders.append((bi, p.bonds[bi].order))
        p.bonds[bi].order = 4
    for i in range(p.n_atoms):
        old_h.append(p.atoms[i].implicit_h)
    cdef list hlog = []
    smi_read_h(p, hlog)
    with mol.edit():
        for item in old_orders:
            bi = <uint32_t> item[0]
            mol.set_order(p.atoms[p.bonds[bi].u].sid, p.atoms[p.bonds[bi].v].sid, 4)
        for i in range(p.n_atoms):
            if p.atoms[i].implicit_h != <int8_t> old_h[i]:
                mol.set_hydrogens(p.atoms[i].sid, p.atoms[i].implicit_h)

    # Does the promoted set have a Kekule form?  This is a QUESTION asked of a discarded copy, not
    # a conversion: the molecule this reader returns has never been kekulised and still holds
    # order-4 bonds.  `stated_h` is passed because it changes the answer -- `[nH]` cannot take a
    # ring double bond and `[n]` must.
    cdef dict stated = {}
    for i in range(p.n_atoms):
        if p.atoms[i].stated_h >= 0:
            stated[p.atoms[i].sid] = p.atoms[i].stated_h
    if kekule(mol.copy(), None, stated).unresolved:
        with mol.edit():
            for item in old_orders:
                bi = <uint32_t> item[0]
                p.bonds[bi].order = <uint8_t> item[1]
                mol.set_order(p.atoms[p.bonds[bi].u].sid, p.atoms[p.bonds[bi].v].sid,
                              <int> item[1])
            for i in range(p.n_atoms):
                if p.atoms[i].implicit_h != <int8_t> old_h[i]:
                    p.atoms[i].implicit_h = <int8_t> old_h[i]
                    mol.set_hydrogens(p.atoms[i].sid, p.atoms[i].implicit_h)
        log.append(mc_record('smiles:promotion-failed', (),
                            'promoting the all-lowercase rings gives an aromatic system with no Kekule '
                            'form, so the bond orders the string wrote are kept instead',
                            mc_refused()))
        return 0
    log.extend(plog)
    log.extend(hlog)
    return 1


cdef int smi_free_stereo_group(uint64_t taken) noexcept nogil:
    """The lowest group id `taken` does not claim, or 0 when all `STEREO_GROUP_MAX` are spent."""
    cdef int i
    for i in range(1, STEREO_GROUP_MAX + 1):
        if not (taken >> i) & 1:
            return i
    return 0


cdef int smi_marks(smi_parse_t *p, MoleculeContainer mol, object log) except -1:
    """Apply the per-atom marks the tail stated: enhanced stereo groups, and CIP descriptors.

    Both in ONE edit scope, and one pass over the atoms, for two reasons.  An atom carrying both a
    group and a descriptor is one atom, so its two refusals belong next to each other in the log --
    which is what a single loop in atom order gives without sorting anything.  And the arena replays
    a scope's descriptor statements after its structural edits, so a scope that states a descriptor
    keeps it wherever in the scope it was stated: two scopes would be two chances to get that
    ordering wrong for no gain.

    Before `smi_stereo`, because this edits and that writes into the arena the edit leaves behind.

    Whether a mark is storable is not asked here.  `set_stereo_group` owns the group field and its
    range, `set_atom_cip` owns the descriptor domain, so each call is made and its `ValueError`
    becomes the log line.  That keeps each domain declared once (RULES.md 6) and means the reader
    reports a refusal in the words of the class that refused it -- including which letters ARE
    descriptors, which is why `smi_brace_cip` accepts any letter and does not guess.  Collected and
    appended after the scope so the log reads in atom order whatever the apply does.

    ONE MARK IS REPAIRED RATHER THAN REFUSED: a group id above `STEREO_GROUP_MAX` is renumbered
    to a free id of its own kind.  The id is a label -- the partition is the statement, which is why
    the stored id is opaque (ruling F79) -- so `&1384:` costs nothing but its spelling, and the CTfile
    reader repairs `MDLV30/STERAC1384` the same way.  An id the tail itself uses is never stolen: the
    first pass takes what is in range before the second assigns.
    """
    cdef uint32_t i
    cdef smi_atom_t *a
    cdef list refused = []
    cdef object e                        # the except-as target; `warn.undeclared` counts it
    cdef uint64_t taken[4]               # per kind, bit <id> set when the tail already states that id
    cdef uint16_t remap_from[STEREO_GROUP_MAX]
    cdef uint8_t remap_kind[STEREO_GROUP_MAX], remap_to[STEREO_GROUP_MAX]
    cdef int n_remap = 0, j, free_id
    cdef int group

    for j in range(4):
        taken[j] = 0
    for i in range(p.n_atoms):
        a = p.atoms + i
        if a.sg_kind < 4 and 1 <= a.sg_group <= STEREO_GROUP_MAX:
            taken[a.sg_kind] |= (<uint64_t> 1) << a.sg_group

    with mol.edit():
        for i in range(p.n_atoms):
            a = p.atoms + i
            if a.sg_kind:
                group = <int> a.sg_group
                if group > STEREO_GROUP_MAX and a.sg_kind < 4:
                    group = 0
                    for j in range(n_remap):
                        if remap_kind[j] == a.sg_kind and remap_from[j] == a.sg_group:
                            group = remap_to[j]
                            break
                    if not group:
                        free_id = smi_free_stereo_group(taken[a.sg_kind])
                        if free_id and n_remap < STEREO_GROUP_MAX:
                            taken[a.sg_kind] |= (<uint64_t> 1) << free_id
                            remap_kind[n_remap] = a.sg_kind
                            remap_from[n_remap] = a.sg_group
                            remap_to[n_remap] = <uint8_t> free_id
                            n_remap += 1
                            group = free_id
                            refused.append(mc_record(
                                'smiles:stereo-group-renumbered', (),
                                'atom %d: the CXSMILES tail names group %d, outside 1..%d, renumbered to %d'
                                % (i + 1, a.sg_group, STEREO_GROUP_MAX, free_id), mc_repaired()))
                        else:
                            group = <int> a.sg_group    # nothing free; let the arena say so below
                try:
                    mol.set_stereo_group(a.sid, <int> a.sg_kind, group)
                except ValueError as e:
                    refused.append(mc_record('smiles:stereo-group-refused', (),
                                             'atom %d: the enhanced stereo group the CXSMILES tail '
                                             'names cannot be stored (%s)' % (i + 1, e), mc_lost()))
            if a.cip:
                # guarded on `a.cip` rather than passing `chr(a.cip)` unconditionally: 0 means the
                # tail named no descriptor for this atom, and `chr(0)` is `'\x00'`, which the arena
                # refuses as a descriptor rather than reading as a clear.  "No descriptor" is spelled
                # by not calling, and `None` if it ever needs to be spelled at all.
                try:
                    mol.set_atom_cip(a.sid, chr(a.cip))
                except ValueError as e:
                    refused.append(mc_record('smiles:cip-refused', (),
                                             'atom %d: the CIP descriptor the `{...}` block names '
                                             'cannot be stored (%s)' % (i + 1, e), mc_lost()))
    log.extend(refused)
    return 0


cdef str smi_kind_name(uint8_t kind):
    """What a stereo unit of this kind is, in the words a log line needs."""
    if kind == SU_CIS_TRANS:
        return 'a cis/trans terminal'
    if kind == SU_ALLENE:
        return 'an allene centre'
    if kind == SU_ATROPISOMER:
        return 'an atropisomer pivot'
    return 'a tetrahedral centre'


cdef inline uint8_t smi_dir_from(smi_parse_t *p, uint32_t bi, uint32_t t) noexcept nogil:
    """The bond's direction read from the terminal `t`: UP means the substituent is up-right of it.

    `dir_` is stored oriented from the atom written FIRST, so reading it from the other end is the
    same statement upside down.  This is the whole content of `/` and `\\`: a direction is a
    statement about a bond and a side, and which side you stand on is not part of the notation.
    """
    if p.bonds[bi].u == t:
        return p.bonds[bi].dir_
    return SMI_DIR_DOWN if p.bonds[bi].dir_ == SMI_DIR_UP else SMI_DIR_UP


cdef inline void smi_perm_of(stereo_unit_t *u, uint32_t *want, uint32_t *perm) noexcept nogil:
    """`perm[i] = j` meaning `want[i]` is `refs[j]`, over four directions.

    The unnamed positions -- an implicit hydrogen, a lone pair -- carry no atom to match on, so they
    are paired off IN ORDER against the unit's `SU_NO_REF` slots, which is what `smw_sign_of` and
    `translate_stereo` both do.  That is sound because `want` was built with each unnamed direction
    already at the position the string put it (ruling F26 and `_inchi.pxi`'s worked example): the
    order is the information, and matching by identity is only how the named ones find their slot.
    """
    cdef uint32_t j, k
    cdef uint32_t norefs[4]
    cdef uint32_t nnoref = 0
    cdef uint32_t nr = 0
    cdef uint32_t used = 0
    perm[0] = 0; perm[1] = 0; perm[2] = 0; perm[3] = 0
    for j in range(4):
        if u.refs[j] == SU_NO_REF:
            norefs[nnoref] = j
            nnoref += 1
    for k in range(4):
        if want[k] == SU_NO_REF:
            if nr < nnoref:
                perm[k] = norefs[nr]
                nr += 1
        else:
            for j in range(4):
                if u.refs[j] == want[k] and not (used & (1u << j)):
                    perm[k] = j
                    used |= 1u << j
                    break


cdef inline uint32_t smi_chain_terminal(smi_parse_t *p, uint32_t centre, uint32_t first,
                                        uint32_t *inward) noexcept nogil:
    """Walk one arm of a cumulene chain from `centre` through `first` and return its terminal.

    `inward` receives the chain atom one step inside that terminal -- `centre` itself for an allene,
    the last-but-one for a longer cumulene -- because the caller needs it to tell the terminal's chain
    bond from its substituents, and one walk is enough to answer both.

    The chain is the maximal run of consecutive order-2 bonds, which is the same chain
    `_stereo.pxi`'s `_cumulene_walk` found when it emitted the unit -- asked here of the parse graph
    because what this pass needs is the WRITTEN order at the far end, which the arena does not keep.
    """
    cdef uint32_t k, bi, nxt, other
    cdef uint32_t prev = centre
    cdef uint32_t cur = first
    cdef uint32_t left = p.n_atoms
    cdef smi_atom_t *a
    # bounded by the atom count rather than left to run until it finds an end, because a ring of
    # nothing but double bonds -- `C1=C=C=C=1`, which the notation permits and which perception's
    # own walk rejects -- has no end and this must not be the loop that discovers that by hanging
    while left:
        left -= 1
        a = p.atoms + cur
        nxt = SMI_NONE
        for k in range(a.nbr_off, a.nbr_off + a.slots):
            bi = p.nbrs[k]
            if bi == SMI_NONE or p.bonds[bi].order != 2:
                continue
            other = p.bonds[bi].v if p.bonds[bi].u == cur else p.bonds[bi].u
            if other != prev:
                nxt = other
                break
        if nxt == SMI_NONE:
            inward[0] = prev
            return cur
        prev = cur
        cur = nxt
    inward[0] = prev
    return cur


cdef int smi_written_pair(smi_parse_t *p, dict index_of, uint32_t t, uint32_t toward,
                          uint32_t *out) except -1:
    """A cumulene terminal's two substituents as arena slots, in the order the string wrote them.

    `toward` is the chain neighbour, whose slot is skipped.  Returns the number of directions found,
    so a terminal the notation has over-filled is refused by the caller rather than truncated here.

    THE UNNAMED DIRECTION HOLDS A WRITTEN POSITION AND A BARE TERMINAL CLAIMS NO SLOT FOR IT, so the
    hydrogen of a bare `C` is inserted at `h_pos` here.  Padding at the END instead is a within-pair
    swap for every terminal whose axis bond is its parent -- the far terminal of `CC=[C@]=CC` reads
    `(C, H)` that way and `(H, C)` this way -- and a within-pair swap inverts the axial parity, so
    the same molecule spelled `C[CH]=[C@]=[CH]C` would come back as the other enantiomer.

    An UNKNOWN count contributes no position, which is `smw_written_h`'s answer too: neither side can
    say how many of an unknown number of hydrogens the string shows, and both saying zero is what
    keeps the writer's order and this one the same order.
    """
    cdef uint32_t k, bi, other, pos
    cdef int n = 0
    cdef smi_atom_t *a = p.atoms + t
    # a stated count already claimed its slot below, so only a bare atom's hydrogens are inserted
    cdef uint32_t unnamed = (0 if a.stated_h >= 0 or a.implicit_h == <int8_t> H_UNKNOWN
                             else <uint32_t> a.implicit_h)
    out[0] = SU_NO_REF
    out[1] = SU_NO_REF
    # one turn past the last slot, so an `h_pos` at the end of a terminal's slots is still reached
    for k in range(a.nbr_off, a.nbr_off + a.slots + 1):
        pos = k - a.nbr_off
        if pos == a.h_pos:
            while unnamed:
                if n < 2:
                    out[n] = SU_NO_REF
                n += 1
                unnamed -= 1
        if pos == a.slots:
            break
        bi = p.nbrs[k]
        if bi == SMI_NONE:                       # the slot a stated hydrogen count claimed
            other = SU_NO_REF
        else:
            other = <uint32_t> index_of[p.atoms[p.bonds[bi].v if p.bonds[bi].u == t
                                                else p.bonds[bi].u].sid]
            if (p.bonds[bi].v if p.bonds[bi].u == t else p.bonds[bi].u) == toward:
                continue
        if n < 2:
            out[n] = other
        n += 1
    return n


cdef int smi_allene_want(smi_parse_t *p, dict index_of, uint32_t centre, uint32_t *want,
                         object log) except -1:
    """The four directions of an allene configuration in the order the string wrote them.

    Two pairs, the chain arm the centre names FIRST leading -- the same shape ruling F26 gives the
    stored `refs`, and `smi_perm_of` reconciles the two orders whichever way round they came out.
    The end exchange is the pair exchange and therefore EVEN (`_stereo.pxi`: "an allene's C2 axis
    performs it"), so a reader that got the two arms the other way round would still be right; the
    within-pair order is where the information is, and that is the one this takes from the string.

    Returns 0 having logged when the notation does not describe an axial frame this can order.
    """
    cdef uint32_t k, bi, arm
    cdef uint32_t inward = 0   # written through a pointer, which Cython cannot see; see `smi_read_h`
    cdef uint32_t arms[2]
    cdef uint32_t n = 0
    cdef int got
    cdef smi_atom_t *a = p.atoms + centre
    for k in range(a.nbr_off, a.nbr_off + a.slots):
        bi = p.nbrs[k]
        if bi == SMI_NONE or p.bonds[bi].order != 2:
            continue
        if n < 2:
            arms[n] = p.bonds[bi].v if p.bonds[bi].u == centre else p.bonds[bi].u
        n += 1
    if n != 2:
        # perception found a chain through this atom, so this cannot happen from a string; it can
        # from one whose promotion raised a chain bond to order 4, and then the arms are not ours
        log.append(mc_record('smiles:stereo-allene-bad-bonds', (),
                            'atom %d states a configuration and is an allene centre, but the string gives it '
                            '%d double bond(s) rather than two; it is not applied' % (centre + 1, n),
                            mc_lost()))
        return 0
    for k in range(2):
        arm = smi_chain_terminal(p, centre, arms[k], &inward)
        got = smi_written_pair(p, index_of, arm, inward, want + 2 * k)
        if got > 2:
            log.append(mc_record('smiles:stereo-allene-overfilled', (),
                                'atom %d states an allene configuration and its end at atom %d has %d '
                                'substituent(s); two is what an axial frame orders, so it is not applied'
                                % (centre + 1, arm + 1, got), mc_lost()))
            return 0
    return 1


cdef uint32_t smi_marked_ref(smi_parse_t *p, dict parse_of, stereo_unit_t *u, uint32_t base,
                             uint32_t terminal, uint8_t *dir_out, set consumed,
                             object log) except 0xFFFFFFFE:
    """Which of a pair's two refs carries a `/` or `\\`, as an index into `refs`.

    `base` is 0 or 2, the first index of the pair; `terminal` is that pair's atom, as a PARSE index.
    Returns SU_NO_REF when neither ref of the pair is marked, and writes the direction read from the
    terminal into `dir_out`.  Every directed bond it looks at goes into `consumed`, whether or not it
    is the one used, so the caller can tell a direction this reader spent from one it never reached.

    Both refs marked is legal and usual -- `C(/F)(\\Cl)=C/Br` says one thing twice -- so the two are
    compared and only a CONTRADICTION is reported.  The first marked ref is the answer either way:
    when they agree the second is redundant, and when they disagree the string is broken and taking
    the first is the only reading that does not invent a third.
    """
    cdef uint32_t k, bi, other
    cdef uint32_t found = SU_NO_REF
    cdef uint8_t d
    for k in range(base, base + 2):
        if u.refs[k] == SU_NO_REF:
            continue
        other = <uint32_t> parse_of[u.refs[k]]
        bi = smi_find_bond(p, terminal, other)
        if bi == SMI_NONE or not p.bonds[bi].dir_:
            continue
        consumed.add(bi)
        d = smi_dir_from(p, bi, terminal)
        if found == SU_NO_REF:
            found = k
            dir_out[0] = d
        elif d == dir_out[0]:
            # the two directions of one terminal must be opposite; equal means the string put both
            # substituents on the same side of the double bond, which no geometry has
            log.append(mc_record('smiles:stereo-same-side', (),
                                'atom %d puts both of its substituents on the same side of the double '
                                'bond; the first direction is used' % (terminal + 1), mc_repaired()))
    return found


cdef int smi_stereo(smi_parse_t *p, MoleculeContainer mol, object log) except -1:
    """Apply what the string said about configuration: `@`/`@@` on an atom, `/` and `\\` on a bond.

    THE PARITY IS WRITTEN INTO SEG_PARITY rather than through the journal, which is what
    `_replay_parities` does and for the same reason: the value is only meaningful against the frame
    the stereo table just derived, and an edit would rebuild the molecule underneath it.
    `refresh_parity_features` then re-derives feature word IV, the one word a parity reaches.

    THE SIGN CONVENTION IS `_smiles_write.pxi`'s, READ BACKWARDS.  `translate_parity` is XOR with
    the permutation's parity, so it is its own inverse for a fixed permutation: the writer's
    `parity -> sign` and this function's `sign -> parity` are one formula in two directions, and the
    anti-drift test is the round trip -- read a string this core wrote and the sign must come back.

    THE CIS/TRANS CONVENTION IS EXTERNAL and was pinned by the InChI epic: parity 1 (even) in a
    unit's own `refs` frame means `refs[0]` and `refs[2]` are TRANS, because a zero-coordinate
    even record of but-2-ene reproduces `InChI=1S/C4H8/c1-3-4-2/h3-4H,1-2H3/b4-3+`, the published
    standard InChI of the (E) isomer (`_inchi.pxi`'s `ICH_CIS_TRANS_FLIP`, and the absolute
    assertions in `test_inchi.py`).  Frame-relative, so the same sentence reads in any frame: parity
    1 means position 0 and position 2 are trans, and a within-pair swap flips it.

    What is NOT applied is named in the log: an allene or atropisomer configuration, `@?`, a
    direction on a bond no unit uses.  Nothing here raises.
    """
    cdef Structure structure = mol._structure
    cdef dict index_of = mol._index_of
    cdef dict parse_of = {}                      # arena slot -> parse index
    cdef set consumed = set()                    # parse bond indices whose direction was read
    cdef uint32_t i, k, j, slot, bi, nb, nd, count
    cdef uint32_t ix, iy, partner, tpar
    cdef uint32_t want[4]
    cdef uint32_t perm[4]
    cdef uint8_t dx = SMI_DIR_NONE
    cdef uint8_t dy = SMI_DIR_NONE
    cdef uint8_t parity
    cdef bint wrote = False
    cdef uint32_t total_dirs = 0
    cdef stereo_unit_t *units
    cdef stereo_unit_t *u
    cdef smi_atom_t *a

    for i in range(p.n_bonds):
        if p.bonds[i].dir_:
            total_dirs += 1
    for i in range(p.n_atoms):
        parse_of[<uint32_t> index_of[p.atoms[i].sid]] = i

    # UNMARKED (ruling F70): every question below is about constitution -- kind, refs, the unnamed
    # nibble -- and none of them is "is this stereogenic?".  A string that states a configuration on
    # an atom that cannot hold two of them is still storing what the string said, and reporting that
    # is `stereo_rejections`' job on a consumer's demand, not a parser's.
    ensure_stereo_units_unmarked(structure)

    for i in range(p.n_atoms):
        a = p.atoms + i
        if a.chiral == SMI_CHIRAL_NONE:
            continue
        if a.chiral == SMI_CHIRAL_UNKNOWN:
            # `@?` says "there is a centre here and I do not know which way it points".  The arena's
            # third state is "nobody said", which is a different sentence, and the stereogenicity
            # this one asserts is derived rather than stored -- so there is nowhere to put it.
            log.append(mc_record('smiles:stereo-unknown', (),
                                'atom %d states an unknown configuration (`@?`); this arena has no state for '
                                '"stated but unresolved" and none is stored' % (i + 1), mc_lost()))
            continue
        slot = <uint32_t> index_of[a.sid]
        u = stereo_unit_of(structure, slot)
        if u is NULL:
            log.append(mc_record('smiles:stereo-no-unit', (),
                                'atom %d states a configuration but nothing here can hold one: it anchors no '
                                'stereo unit in this molecule' % (i + 1), mc_lost()))
            continue
        if u.n_refs != 4:
            log.append(mc_record('smiles:stereo-unit-incomplete', (),
                                'atom %d states a configuration over %d direction(s); four is what %s '
                                'orders, so it is not applied'
                                % (i + 1, u.n_refs, smi_kind_name(u.kind)), mc_lost()))
            continue
        if u.kind == SU_ALLENE:
            # `NC(Br)=[C@]=C(O)C`.  OpenSMILES calls it allene-like and means it literally: the four
            # substituents of the two CHAIN ENDS stand in for the centre's own neighbours and the
            # ordinary tetrahedral sentence is then read over them.  So the sign below is the same
            # expression, and the only new work is finding the four in written order.
            if not smi_allene_want(p, index_of, i, want, log):
                continue
        elif u.kind == SU_TETRA:
            nd = a.slots
            if nd > 4:
                log.append(mc_record('smiles:stereo-too-many-dirs', (),
                                    'atom %d states a configuration and names %d directions; four is the '
                                    'most a tetrahedral frame orders, so it is not applied' % (i + 1, nd),
                                    mc_lost()))
                continue
            if a.stated_h > 1:
                # two directions with no atom of their own cannot be told apart, so an order over
                # them is not an order.  `[C@H2]` is the string that does this
                log.append(mc_record('smiles:stereo-ambiguous-h', (),
                                    'atom %d states %d hydrogens and a configuration; two unnamed directions '
                                    'cannot be ordered and it is not applied' % (i + 1, a.stated_h),
                                    mc_lost()))
                continue
            # the written order.  The parse graph already holds it: a slot per bond in the order the
            # string names them, with the implicit hydrogen's slot claimed where it was WRITTEN --
            # which is what makes `[C@H](F)(Cl)Br` and `F[C@H](Cl)Br` opposite configurations here
            # rather than something this pass has to correct.  The pad at the end is the lone pair's,
            # last for the same reason `_smiles_write.pxi` puts it last: it has no positional rule.
            k = 0
            for j in range(a.nbr_off, a.nbr_off + a.slots):
                bi = p.nbrs[j]
                if bi == SMI_NONE:
                    want[k] = SU_NO_REF
                else:
                    nb = p.bonds[bi].v if p.bonds[bi].u == i else p.bonds[bi].u
                    want[k] = <uint32_t> index_of[p.atoms[nb].sid]
                k += 1
            while k < 4:
                want[k] = SU_NO_REF
                k += 1
        else:
            log.append(mc_record('smiles:stereo-unsupported-unit', (),
                                'atom %d states a configuration and is %s, whose frame this reader does not '
                                'build yet; it is not applied' % (i + 1, smi_kind_name(u.kind)),
                                mc_lost()))
            continue

        smi_perm_of(u, want, perm)
        parity = translate_parity(2 if a.chiral == SMI_CHIRAL_AT else 1, perm)
        structure_set_parity(structure, slot, parity)
        wrote = True

    # the bond directions, driven by the perceived units rather than by a double-bond walk of this
    # reader's own: which bonds can carry a configuration is a question `_stereo.pxi` answers, and
    # asking it here means a cumulene, a ring too small to hold one and a plain alkene all arrive
    # through the same door.
    if total_dirs:
        count = structure_stereo_unit_count(structure)
        units = structure_stereo_units(structure)
        for k in range(count):
            u = &units[k]
            if u.kind != SU_CIS_TRANS or u.n_refs != 4:
                continue
            slot = u.anchor
            partner = stereo_unit_partner(structure, u)
            if partner == SU_NO_REF:
                continue
            tpar = <uint32_t> parse_of[slot]
            ix = smi_marked_ref(p, parse_of, u, 0, tpar, &dx, consumed, log)
            iy = smi_marked_ref(p, parse_of, u, 2, <uint32_t> parse_of[partner], &dy, consumed, log)
            if ix == SU_NO_REF and iy == SU_NO_REF:
                continue
            if ix == SU_NO_REF or iy == SU_NO_REF:
                log.append(mc_record('smiles:stereo-one-sided', (),
                                    'the double bond between atoms %d and %d has a direction on one side '
                                    'only, so no configuration is stored'
                                    % (tpar + 1, <uint32_t> parse_of[partner] + 1), mc_lost()))
                continue
            # In the frame `(refs[ix], refs[ix^1], refs[iy], refs[iy^1])` the two marked
            # substituents are at positions 0 and 2, so the convention reads directly: even is
            # trans.  Opposite directions read from their own terminals means opposite sides.
            perm[0] = ix; perm[1] = ix ^ 1; perm[2] = iy; perm[3] = iy ^ 1
            parity = translate_parity(1 if dx != dy else 2, perm)
            structure_set_parity(structure, slot, parity)
            wrote = True
        if <uint32_t> len(consumed) < total_dirs:
            log.append(mc_record('smiles:stereo-unused-dirs', (),
                                '%d bond direction(s) name no configuration this molecule can hold; a `/` or '
                                '`\\` states nothing on its own and they were dropped'
                                % (total_dirs - <uint32_t> len(consumed)), mc_lost()))
    if wrote:
        refresh_parity_features(structure)
    return 0


def read_smiles(text, log=None):
    """Read a SMILES string into a `MoleculeContainer`, or a reaction SMILES into a `ReactionContainer`.

    POLYMORPHIC ON THE ARROW.  A `>` that is not part of a dative `->` makes the string a reaction and
    the result a `ReactionContainer`; without one the result is a molecule.  `smiles` is the
    bidirectional door over this function, so a corpus that mixes the two shapes reads in one loop.
    `read_reaction_smiles` is the strict door beside this one, for a caller who would rather the wrong
    shape failed than came back as a molecule.  Everything below describes the molecule half; the
    reaction half is documented on `read_reaction_smiles`, applies this reader per component, and
    shares this `log`.

    Aromatic input stays aromatic: `c1ccccc1` gives six order-4 bonds and nothing here calls
    `kekule()` to convert them.  Kekule input stays Kekule.  Which one you get is which one the
    string said, and converting is the caller's call.

    EVERY LINE LANDS ON THE RETURNED CONTAINER'S `log`, always.  Every place this reader had to prefer
    one reading of a contradictory string over another -- a ring bond whose two labels disagree, an
    all-lowercase ring whose bonds were written single, an atom with no valence rule, an element with no
    aromatic form -- is one record on `mol.log`, under stage `read`.  Atoms are named by their 1-based
    position among the atoms of the string, which is also their stable id in the returned molecule.

    `log` is an optional list to append the same lines to, for a caller reading a whole file who wants
    one sequence for the run rather than one per record.  Omitting it costs nothing now: it is a second
    view, not the storage, and there is no arrangement in which a line is written nowhere.

    A ` |...|` CXSMILES tail is read for its `^N:` radicals, which change the hydrogen count of a
    bare atom, and for its `a:` / `o<n>:` / `&<n>:` enhanced stereo groups; its other fields are
    named in the log and not applied.  A ` {...}` block is the same idea in another dialect: the
    three stereo-group fields are spelled identically there, and `A<index>=<letter>` states an atom's
    CIP descriptor.

    Raises `IncorrectSmiles`, with a byte offset, for SYNTAX: an unknown element spelling, an
    unbalanced parenthesis, a ring label that never closes, a field the arena cannot store, a token
    with no meaning here (`*`, `$`).  `~` and the arrows `->` / `<-` DO have a meaning: they are
    the dative bond, order 8, which the arena stores and `write_smiles` emits.  Never for chemistry
    -- a hypervalent atom or an
    unkekulisable ring is stored and logged, because a reader that refuses what other tools emit is
    a reader nobody can put in front of a database -- and never for the tail, which is an extension
    whose failure does not make the SMILES unreadable.
    """
    cdef bytes raw
    if isinstance(text, str):
        try:
            raw = (<str> text).encode('ascii')
        except UnicodeEncodeError:
            raise IncorrectSmiles('the string contains a non-ASCII character') from None
    elif isinstance(text, bytes):
        raw = <bytes> text
    else:
        raise TypeError('read_smiles takes a str or bytes')
    raw = raw.strip()

    cdef const char *s = <const char *> raw
    cdef uint32_t n = <uint32_t> len(raw)
    cdef uint32_t k = 0
    while k < n and s[k] > 32:
        k += 1
    cdef bytes tail = raw[k:].strip()

    # POLYMORPHIC: an arrow makes this a reaction and the return a `ReactionContainer`.  `smiles` IS
    # this function, so a corpus of mixed records keeps reading through one call.  The dispatch is on
    # the SAME rule the split itself uses -- a `>` preceded by `-` is a dative bond and not the arrow
    # -- because two spellings of "is there an arrow" would eventually disagree, and the one that
    # disagreed would send `N->[Cu]` to the wrong reader.  `read_reaction_smiles` is the strict door
    # for a caller who wants the shape checked.
    cdef uint32_t i
    for i in range(k):
        if s[i] == c'>' and (i == 0 or s[i - 1] != c'-'):
            return read_reaction_smiles(raw, log)

    cdef object mylog = log if log is not None else []
    return smi_one(raw[:k], tail, mylog)


cdef MoleculeContainer smi_one(bytes src, bytes tail, object mylog):
    """One component, from bytes to a sealed molecule: tokenize, tail, hydrogens, build, stereo.

    EXTRACTED FROM `read_smiles` RATHER THAN COPIED, because `read_reaction_smiles` needs the same
    seven steps in the same order per component and a second spelling of that order is a second
    reader.  The order is load-bearing twice over -- the tail's `^N:` radicals change hydrogen counts
    so they precede `smi_read_h`, and the stereo groups edit the arena so they precede `smi_stereo` --
    and neither constraint is visible from a call site.

    `tail` is this component's own extension block, which for a reaction is EMPTY: a reaction has one
    tail for the whole string whose indices span every side, so it cannot be applied one component at
    a time.  `read_reaction_smiles` applies it instead, between the two halves below, which is the
    whole reason those halves have names.

    THE MOLECULE GETS EVERY LINE THIS READ PRODUCED, on `mol.log`, whether or not the caller passed a
    `log=`.  The list is still filled for a caller reading a whole file's records in one place; the
    container is where a molecule handed on alone carries its own.
    """
    cdef uint32_t start = <uint32_t> len(mylog)
    cdef MoleculeContainer built
    cdef smi_parse_t p
    smi_alloc(&p, <uint32_t> len(src))
    try:
        smi_parse_read(&p, src, mylog)
        # before the hydrogen pass, which the tail's radicals change
        if tail:
            if tail.startswith(b'|') and tail.endswith(b'|') and len(tail) > 1:
                smi_cx(&p, 1, tail, mylog, False)
            elif tail.startswith(b'{') and tail.endswith(b'}') and len(tail) > 1:
                smi_brace(&p, tail, mylog)
            elif tail.startswith(b'|') or tail.startswith(b'{'):
                mylog.append(mc_record('smiles:extension-unterminated', (),
                                      'the extension block after the SMILES is not terminated and was '
                                      'ignored: %s' % tail.decode('ascii', 'replace'), mc_lost()))
            else:
                mylog.append(mc_record('smiles:extension-ignored', (),
                                      'text after the SMILES is not part of it and was ignored: %s'
                                      % tail.decode('ascii', 'replace'), mc_lost()))
        built = smi_parse_build(&p, src, tail, mylog)
        # the `if` asks whether there is anything to say, never whether to say it: `mol.log` builds its
        # storage on first touch, 0.2 us against a 5 us parse, and a clean string should not pay for an
        # empty one.  `_ctab.py`'s fold is guarded the same way and for the same reason
        if <uint32_t> len(mylog) > start:
            built.log.absorb('read', mylog[start:])
        return built
    finally:
        smi_free(&p)


cdef int smi_parse_read(smi_parse_t *p, bytes src, object mylog) except -1:
    """The first half: tokenize and flatten, leaving the parse state ready for a tail."""
    smi_tokenize(p, <const char *> src, <uint32_t> len(src), mylog)
    smi_flat(p)
    return 0


cdef MoleculeContainer smi_parse_build(smi_parse_t *p, bytes src, bytes tail, object mylog):
    """The second half: hydrogens, build, promotion, marks, stereo, labels.  Runs after the tail.

    `src` is this component's own bytes and `tail` the whole string's extension block; the label
    offsets in the parse graph point into one or the other.
    """
    cdef list hlog
    cdef MoleculeContainer mol
    # when promotion may run, the first hydrogen pass describes orders that may be replaced, so its
    # lines are held back until it is known whether they survive
    hlog = mylog if not p.promote_hint else []
    smi_read_h(p, hlog)
    mol = smi_build(p)
    if p.promote_hint and not smi_promote(p, mol, mylog):
        mylog.extend(hlog)
    if p.dative_arrows:
        # the order is stored, the arrow is not: nothing in the arena distinguishes a donor from an
        # acceptor, and `write_smiles` spells every order-8 bond `~`.  One line rather than silence,
        # because which atom donated is information the string carried in
        mylog.append(mc_record('smiles:dative-arrow-direction-lost', (),
                               '%d dative bond(s) were written as an arrow; order 8 is stored but which '
                               'atom donates is not' % p.dative_arrows, mc_lost()))
    if p.sg_marks or p.cip_marks:
        smi_marks(p, mol, mylog)
    if p.chiral_count or p.dir_count:
        # after the groups, which edit, because this writes into the arena the edit leaves
        smi_stereo(p, mol, mylog)
    smi_labels(p, src, tail, mol)
    return mol


cdef int smi_labels(smi_parse_t *p, bytes src, bytes tail, MoleculeContainer mol) except -1:
    """Store every label as its atom's alias.  Last, because `set_aliases` wants a sealed molecule
    and because it REPLACES the set -- nothing else in this reader writes one.

    A label came either from a bracket, whose offset is into `src`, or from the tail's `$...$` field,
    whose offset is into `tail`; `label_tail` says which.  Only the tail's text can carry a `&#NN;`
    character reference, a bracket label being letters and digits.
    """
    cdef uint32_t i
    cdef smi_atom_t *a
    cdef dict labels = {}
    for i in range(p.n_atoms):
        a = p.atoms + i
        if a.label_len:
            if a.label_tail:
                labels[a.sid] = smi_cx_unescape(tail[a.label_off:a.label_off + a.label_len])
            else:
                labels[a.sid] = src[a.label_off:a.label_off + a.label_len]
    if labels:
        mol.set_aliases(labels)
    return 0


# The extension cannot import `reaction.py` -- `reaction.py` imports from `._core` -- so the container
# is injected the way `_ich_set_kekule_fn` injects the kekule pass, and for the same reason.
cdef object smi_reaction_factory = None


def _set_reaction_factory(factory):
    """Register `ReactionContainer` as what `read_reaction_smiles` builds.  Called by `core/reaction.py`."""
    global smi_reaction_factory
    smi_reaction_factory = factory


cdef int smi_rxn_own(object mylog, list owner_of, uint32_t base, int owner) except -1:
    """Attribute every log line appended since the last call: a component index, or -1 for the reaction.

    A reaction's components are all tokenized before any is built, so one component's lines are two runs
    with other components' in between.  Marking as it goes is what lets each molecule absorb its own and
    the reaction stamp the same records with the subject that makes their atom ids readable.
    """
    while base + <uint32_t> len(owner_of) < <uint32_t> len(mylog):
        owner_of.append(owner)
    return 0


cdef list smi_components(bytes side):
    """One side into its components, splitting on a `.` that is neither bracketed nor in a branch.

    A `.` inside `[...]` cannot occur in valid SMILES but can occur in the garbage this reader is
    required to accept, and one inside a branch belongs to the branch's own component -- so depth is
    tracked rather than assumed.  An empty side gives no components, which is how `CC>>` reads.
    """
    cdef const char *s = <const char *> side
    cdef uint32_t n = <uint32_t> len(side), i, start = 0
    cdef int depth = 0
    cdef bint bracket = False
    cdef list out = []
    if not n:
        return out
    for i in range(n):
        if bracket:
            if s[i] == c']':
                bracket = False
        elif s[i] == c'[':
            bracket = True
        elif s[i] == c'(':
            depth += 1
        elif s[i] == c')':
            depth -= 1
        elif s[i] == c'.' and depth == 0:
            out.append(side[start:i])
            start = i + 1
    out.append(side[start:])
    return out


cdef list smi_cx_fgroups(bytes block, object log):
    """The tail's `f:` fields, as lists of COMPONENT indices.  Read before anything is tokenized.

    `f:` is the only field that does not name atoms, and the only one a reaction has to know before it
    parses: a group says "these components are one molecule", and the way to make them one molecule is
    to tokenize their text together.  Doing it that way rather than by unioning built molecules is not
    a shortcut -- `union` RENUMBERS enhanced stereo group ids, so a tail that puts atoms of two grouped
    components in one `&1` would come back with them in two, which is a different chemical claim.

    Field scanning is `smi_cx_field_end`'s.  The index scanning is NOT `smi_cx_next_index`'s, because
    this is the one field with two separators: `.` continues a group and `,` starts the next.
    """
    cdef const char *s = <const char *> block
    cdef uint32_t n = <uint32_t> len(block)
    cdef uint32_t i = 1                              # past the opening `|`
    cdef uint32_t j, k, value
    cdef char c
    cdef list out = []
    cdef list group
    while i + 1 < n:                                 # `n - 1` is the closing `|`
        c = s[i]
        if c == 44 or c == 32:
            i += 1
            continue
        j = smi_cx_field_end(s, n, i)
        if j <= i:
            break
        if c == 102 and i + 1 < j and s[i + 1] == 58:             # `f:`
            k = i + 2
            group = []
            while k < j:
                if not smi_isdigit(s[k]):
                    log.append(mc_record('smiles:cx-fgroup-malformed', (),
                                        'the CXSMILES component group field `%s` is malformed and the rest of '
                                        'it was dropped' % block[i:j].decode('ascii', 'replace'), mc_lost()))
                    break
                value = 0
                while k < j and smi_isdigit(s[k]):
                    value = value * 10 + <uint32_t> (s[k] - 48)
                    k += 1
                group.append(value)
                if k >= j:
                    break
                elif s[k] == 46:                     # `.` -- the same group continues
                    k += 1
                elif s[k] == 44:                     # `,` -- the next group starts
                    k += 1
                    out.append(group)
                    group = []
                else:
                    log.append(mc_record('smiles:cx-fgroup-malformed', (),
                                        'the CXSMILES component group field `%s` is malformed and the rest of '
                                        'it was dropped' % block[i:j].decode('ascii', 'replace'), mc_lost()))
                    break
            if group:
                out.append(group)
        i = j
    return out


cdef list smi_rxn_merge(list comps, list owners, list fgroups, object log):
    """Apply the `f:` groups by joining each group's component text with `.`.

    Returns the merged component list; `owners` is edited in step with it.  A group is applied only
    when it names components that exist, that share a side, and that are CONSECUTIVE -- consecutive
    because the tail's atom indices are positions in the string, so joining text out of order would
    move atoms out from under them.  Every writer that emits `f:` emits its groups consecutively,
    chython's own included; a group that is not gets a line and is left ungrouped, which costs the
    reaction a molecule boundary and nothing else.
    """
    cdef uint32_t count = <uint32_t> len(comps)
    cdef list target = list(range(count))            # which component each one is folded into
    cdef list group
    cdef list spelled
    cdef set owned
    cdef bint taken
    cdef object member
    cdef uint32_t low, high
    cdef str spelling
    # EXPLICIT LOOPS AND NOT COMPREHENSIONS, for the reason `DetachedSmiles.join` states: a
    # comprehension gets its own scope in Cython 3 and `warn.undeclared` is an error here, so a target
    # declared in this function does not satisfy it.
    for group in fgroups:
        group = sorted(set(group))
        if not group:
            continue
        # `group[len(group) - 1]` and not `group[-1]`: this translation unit compiles with
        # `wraparound=False`, under which a negative index into a `cdef list` reads off the end
        low = <uint32_t> group[0]
        high = <uint32_t> group[len(group) - 1]
        spelled = []
        for member in group:
            spelled.append(str(member))
        spelling = '.'.join(spelled)
        if high >= count:
            log.append(mc_record('smiles:fgroup-bad-index', (),
                                'the CXSMILES component group `f:%s` names component %d, but the reaction has '
                                '%d; the group was not applied' % (spelling, high, count), mc_lost()))
            continue
        owned = set()
        taken = False
        for member in group:
            owned.add(owners[<uint32_t> member])
            if <uint32_t> target[<uint32_t> member] != <uint32_t> member:
                taken = True
        if len(owned) != 1:
            log.append(mc_record('smiles:fgroup-cross-side', (),
                                'the CXSMILES component group `f:%s` spans more than one side of the reaction; '
                                'a molecule has one role, so the group was not applied' % spelling,
                                mc_lost()))
        elif high - low != <uint32_t> len(group) - 1:
            log.append(mc_record('smiles:fgroup-non-consecutive', (),
                                'the CXSMILES component group `f:%s` names components that are not consecutive '
                                'in the string; the group was not applied' % spelling, mc_lost()))
        elif taken:
            log.append(mc_record('smiles:fgroup-already-taken', (),
                                'the CXSMILES component group `f:%s` names a component another group already '
                                'took; the second group was not applied' % spelling, mc_lost()))
        else:
            # `low` is itself a member, and `target[low] = low` is what it already held
            for member in group:
                target[<uint32_t> member] = low

    cdef list out = []
    cdef list out_owners = []
    cdef dict slot = {}
    cdef uint32_t i, k
    for i in range(count):
        if target[i] == i:
            slot[i] = len(out)
            out.append(comps[i])
            out_owners.append(owners[i])
        else:
            k = slot[target[i]]
            out[k] = out[k] + b'.' + <bytes> comps[i]
    owners[:] = out_owners
    return out


def read_reaction_smiles(text, log=None):
    """Read a reaction SMILES into a `ReactionContainer`.

    `reactants>agents>products`, each side split on `.` into one molecule per component and each
    component read by the same tokenizer `read_smiles` uses.  Any side may be empty, so `CC>>`, `>>CC`
    and `>>` all read.  A `>` that belongs to a dative `->` is not the arrow, so `N->[Cu]>>N` reads as
    one reactant and one product rather than a record with three separators.

    ONE ` |...|` CXSMILES tail for the whole string, at the very end, and its atom indices count every
    atom of every side in written order -- reactants, then agents, then products.  Its `f:` field is
    APPLIED here and not merely logged: `f:0.1` says two components are one molecule, and without it
    `[Na+].[Cl-]>>` comes back as two reactants instead of one salt.

    EVERY LINE LANDS ON `rxn.log`, and a line about one component lands on that component's own `log`
    as well -- the arrangement a reaction pass uses, so `rxn.log.by_subject('products[0]')` and
    `rxn.products[0].log` answer the same question.  The reaction-level copy carries the `subject`,
    since pooled atom ids name a different atom in each container; a line about the tail, an `f:` group
    or a dropped component belongs to no component and carries none.

    `log` is an optional list to append the same lines to, shared by every component; atoms in a line
    are named by their 1-based position in the WHOLE string, matching the tail's index space.  Refusals
    are for syntax only, exactly as in `read_smiles` -- a chemically impossible reactant is stored and
    logged, and nothing in the tail raises.
    """
    cdef bytes raw
    if isinstance(text, str):
        try:
            raw = (<str> text).encode('ascii')
        except UnicodeEncodeError:
            raise IncorrectSmiles('the string contains a non-ASCII character') from None
    elif isinstance(text, bytes):
        raw = <bytes> text
    else:
        raise TypeError('read_reaction_smiles takes a str or bytes')
    raw = raw.strip()

    if smi_reaction_factory is None:
        raise RuntimeError('no reaction container is registered; chython.core.reaction must be '
                           'imported so that it can call _set_reaction_factory')

    # the body is the first whitespace-free run and the rest is the tail, exactly as `read_smiles`
    # splits them -- so a `.smi` file's title column behaves the same way on both readers
    cdef const char *s = <const char *> raw
    cdef uint32_t n = <uint32_t> len(raw), i
    cdef uint32_t k = 0
    while k < n and s[k] > 32:
        k += 1
    cdef bytes tail = raw[k:].strip()
    cdef bytes body = raw[:k]
    n = k

    # A `>` PRECEDED BY `-` IS A DATIVE BOND, NOT A SEPARATOR.  `N->[Cu]>>N` has three `>` bytes and
    # one arrow, so a reader that counts bytes refuses the whole record.  The rule is safe in the other
    # direction too: a `-` immediately before the reaction arrow would be a dangling bond at the end of
    # a component, which is not something any writer emits and which the tokenizer refuses on its own
    # if it does appear.
    cdef list cuts = []
    for i in range(n):
        if s[i] == c'>' and (i == 0 or s[i - 1] != c'-'):
            cuts.append(i)
    if len(cuts) != 2:
        raise IncorrectSmiles('a reaction SMILES has two `>` separators, giving '
                              '`reactants>agents>products`; this string has %d' % len(cuts))

    cdef object mylog = log if log is not None else []
    # every line from here on is attributed as it is written: `owner_of[j]` is the component
    # `mylog[base + j]` is about, or -1 for a line about the reaction (a tail field, a dropped
    # component, an `f:` group).  See `smi_rxn_own`.
    cdef uint32_t base = <uint32_t> len(mylog)
    cdef list owner_of = []
    cdef bint is_cx = tail.startswith(b'|') and tail.endswith(b'|') and len(tail) > 1
    if tail and not is_cx:
        if tail.startswith(b'|') or tail.startswith(b'{'):
            mylog.append(mc_record('smiles:extension-unterminated', (),
                                   'the extension block after the reaction SMILES is not terminated and was '
                                   'ignored: %s' % tail.decode('ascii', 'replace'), mc_lost()))
        else:
            mylog.append(mc_record('smiles:extension-ignored', (),
                                   'text after the reaction SMILES is not part of it and was ignored: %s'
                                   % tail.decode('ascii', 'replace'), mc_lost()))

    # one flat component list in the tail's own order -- reactants, agents, products -- with the side
    # each component came from, so `f:` and the atom indices are read against the same numbering
    cdef list comps = []
    cdef list owners = []
    cdef list bodies = [body[:cuts[0]], body[cuts[0] + 1:cuts[1]], body[cuts[1] + 1:]]
    cdef uint32_t side
    cdef object part
    for side in range(3):
        for part in smi_components(<bytes> bodies[side]):
            comps.append(part)
            owners.append(side)
    if is_cx:
        comps = smi_rxn_merge(comps, owners, smi_cx_fgroups(tail, mylog), mylog)

    cdef uint32_t count = <uint32_t> len(comps)
    cdef smi_parse_t *ps = NULL
    cdef uint32_t ready = 0
    cdef list mols
    if count:
        ps = <smi_parse_t *> PyMem_Malloc(count * sizeof(smi_parse_t))
        if ps is NULL:
            raise MemoryError('reaction SMILES parse allocation failed')
    try:
        # EVERY COMPONENT IS TOKENIZED BEFORE ANY IS BUILT, because the tail is read between the two
        # and it names atoms of all of them.  That is also why each side's ring labels stay its own:
        # one parse state per component means `C1CC>>C1CC` cannot silently bond across the arrow, it
        # is refused twice for the unclosed label it actually has.
        for i in range(count):
            smi_alloc(&ps[i], <uint32_t> len(<bytes> comps[i]))
            ready = i + 1
        smi_rxn_own(mylog, owner_of, base, -1)     # the tail and `f:` lines above are the reaction's
        for i in range(count):
            smi_parse_read(&ps[i], <bytes> comps[i], mylog)
            smi_rxn_own(mylog, owner_of, base, <int> i)
        if is_cx:
            # the tail's indices span every side, so a line about it belongs to no one component
            smi_cx(ps, count, tail, mylog, True)
            smi_rxn_own(mylog, owner_of, base, -1)
        mols = []
        for i in range(count):
            mols.append(smi_parse_build(&ps[i], <bytes> comps[i], tail, mylog))
            smi_rxn_own(mylog, owner_of, base, <int> i)
    finally:
        for i in range(ready):
            smi_free(&ps[i])
        PyMem_Free(ps)

    cdef list sides = [[], [], []]
    cdef list side_names = ['reactants', 'agents', 'products']
    cdef list subject_of = []       # per component, the `reactants[0]` a record of it is stamped with
    for i in range(count):
        if (<MoleculeContainer> mols[i]).atom_count:
            subject_of.append('%s[%d]' % (<str> side_names[owners[i]], len(<list> sides[owners[i]])))
            (<list> sides[owners[i]]).append(mols[i])
        else:
            # a `..` in a side, or a side that is nothing but separators.  Stored nowhere and logged:
            # an atomless molecule is not a participant, and the empty component was never one
            subject_of.append('')
            mylog.append(mc_record('smiles:empty-component', (),
                                   'component %d of the reaction is empty and was dropped' % (i + 1),
                                   mc_lost()))
            smi_rxn_own(mylog, owner_of, base, -1)

    cdef object rxn = smi_reaction_factory(reactants=sides[0], products=sides[2], agents=sides[1])
    # THE COMPONENT KEEPS ITS OWN RECORDS AND THE REACTION GETS A COPY, the arrangement
    # `_reaction_passes._mirror` uses for a pass and for the same reason: `LogRecord.atoms` are stable
    # ids in ONE container, so the pooled copy needs the `subject` that says which.  In string order on
    # `rxn.log`, so a reader sees the record the way the record was read.
    cdef uint32_t j, lines = <uint32_t> len(owner_of)
    cdef int who
    cdef list bucket
    for i in range(count):
        bucket = []
        for j in range(lines):
            if <int> owner_of[j] == <int> i:
                bucket.append(mylog[base + j])
        if bucket:
            (<MoleculeContainer> mols[i]).log.absorb('read', bucket)
    for j in range(lines):
        who = <int> owner_of[j]
        rxn.log.absorb('read', [mylog[base + j]],
                       subject='' if who < 0 else <str> subject_of[who])
    return rxn
