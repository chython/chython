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
# The SMARTS reader: a string becomes a QueryContainer, in one pass over the bytes and with no
# parse graph at all.
#
# WHY NO PARSE GRAPH
#
# `_smiles_read.pxi` needs one because `@` is a parity over neighbour slots in the order the string
# names them, so nothing can be built until every slot is known.  A query has no such ordering
# problem: `_query_seal.pxi` gathers a stable id's QOP_ATOM_TOKENs in JOURNAL ORDER, filtered per
# atom, so tokens for different atoms may be interleaved freely.  This reader therefore calls
# `add_atom` / `add_bond` / `atom_primitive` as it reads and keeps only what it cannot yet place:
# the bond expression in front of an atom that has not been read, and the one in front of a ring
# closure label whose partner has not been reached.  Those live in a token pool.
#
# THE DIALECT
#
# `;` is AND between primitives, `,` is OR, `&` and bare juxtaposition are the high-precedence AND.
# Isotope, charge, atom map and configuration are ordinary primitives in that grammar, so `[13C@:7]`
# works by construction and each of them works in any position.
#
# TEN LETTERS ARE PRIMITIVES, NOT ELEMENTS
#
# `A D H M R` and `a h r x z` never begin a ONE-letter element symbol inside a bracket.  The rule is
# not a preference: SMI_ELEMENT spells `H` as hydrogen and `Xe` as xenon, so without it `[C;h2]`
# reads as a hydrogen atom followed by a stray `2` and `[C;x0]` cannot be told from a xenon that
# forgot its `e`.  Two-letter symbols are untouched for the uppercase five -- `Dy`, `Ho`, `Mg`, `Ag`,
# `Ru` all still read -- and for the lowercase five the letter wins outright, which costs the
# unreachable spellings `as`, `he`, `rn`, `xe`, `zn`.  `[As;a]` is how this codebase writes aromatic
# arsenic anyway.
#
# Three of the ten ARE the element when they are the first primitive in the bracket, which is what the
# corpus relies on: `[H]` is a hydrogen atom, `[A]` is any atom, `[M]` is any metal.  Later in the same
# bracket the same letters are `total_h`, nothing at all (a trailing `A` is accepted and ignored) and
# the masked flag.
#
# COMPONENT GROUPING, AND WHY IT LIVES IN THE SMARTS READER
#
# `.` has always parsed here, and it says only "not bonded": the two fragments may land in one
# molecule or in two, and the query does not care.  Daylight's component-level grouping is the way to
# care, and the query arena has had the field for it since the isomorphism kernel was written --
# `qcomp_t.group`, matched as "same group, same molecule component; different groups, different
# components; no group, unconstrained".  Nothing could reach it from a string: `set_group` was
# API-only.
#
# So `(` at a COMPONENT POSITION -- no atom to branch from, no branch open -- opens a group instead of
# refusing, and `A.B` / `(A).(B)` / `(A.B)` are three different questions.  The refusal it replaces was
# `branch opens before any atom`, which is still what a `(` gets anywhere a group cannot be meant.
# The reason this is here rather than in the SMIRKS reader above it is that there is ONE lexer: an
# intramolecular reaction is a template whose reactant side is grouped, and a grouped SMARTS is useful
# on its own the moment it can be written.
#
# SIX DECISIONS THIS FILE IS BUILT ON
#
# 1. `~` IS ANY BOND, AND IT MATCHES.  It compiles to the five-way disjunction `order 1 , order 2 ,
#    order 3 , aromatic , order 8`.  `!~` forbids every bond and is refused with a position rather
#    than compiled into something that matches nothing.  IN V2 `~` IS ORDER 8, the coordination bond
#    (`_query_boxes.pxi`'s `prim_apply`, which refuses to negate it), so a `~` in a ported template
#    asks a narrower question there than it does here.
#
# 2. A BARE LOWERCASE ATOM IS AROMATIC, BOND INCLUDED.  The atom gets `hybridization 4` and an
#    implicit bond between two atoms both WRITTEN lowercase gets `bond_aromatic`, exactly as in
#    SMILES.  V2 carries no aromatic constraint on either, so `smarts('c1ccccc1')` there is a
#    single-bonded carbocycle.
#
# 3. AN IMPLICIT BOND IS OTHERWISE SINGLE, with no exception for aromatic atoms: `[C;a]:[C;a]` needs
#    its colon, and 185 `a` primitives in this codebase are written that way.  Decision 2 does not
#    weaken this one -- it is about atoms written lowercase, which no template does.
#
# 4. SYNTAX RAISES, CHEMISTRY IS LOGGED, as in `_smiles_read.pxi`.  `IncorrectSmarts` subclasses
#    `IncorrectSmiles` so a caller that catches the SMILES error catches this one too, and every
#    message ends in a byte offset.  A primitive that says something the box layout cannot state --
#    `r1`, `!@@`, `!~` -- is syntax, because a located error beats a query that matches nothing.
#
# 5. WHAT THE STRING SAYS IS WHAT IT MEANS, primitive by primitive.  Three measured V2 readings are
#    why that needs stating, and each is pinned by an acceptance test in `test/test_smarts_read.py`:
#    `[N+,O]` there is `[N;+]`, the body split building one element list and hoisting the charge onto
#    the atom, so the oxygen alternative is not asked for; `[C;D1;D2]` there is `[C;D2]`, a second
#    clause about one field replacing the first; and `[0C]` there is plain carbon with the zero
#    dropped, where here it is the `no_isotope` primitive.
#
# 6. `^` IS THE DATIVE BOND, AND IT IS NOT `~`.  Decision 1 spends `~` on "any bond", which is what
#    Daylight means by it and what a query author reaching for a wildcard expects -- so the
#    coordination bond, which chython SMILES writes `~`, needed a character of its own.  `^` was the
#    only free one: `>` and `->` are unusable because `_smirks_read.pxi` counts `>` bytes and refuses
#    any string without exactly two, and `%` is a ring label, `/` and `\` are stereo, `$` would
#    foreclose recursive SMARTS.  So the two dialects spell this bond differently, on purpose, and
#    `[Fe]~N(C)(C)C` is matched by `[M]^[N;D3]`.
#
#    `!^` IS THE USEFUL HALF and it expands rather than refusing.  A box cannot state "not order 8"
#    (see the refusal in `_query_boxes.pxi`), so the reader writes the disjunction the box layer
#    tells callers to write by hand: `-,=,#,:`.  It shares one caveat with `~`, which expands the
#    same way: the alternatives are pushed into the token stream with `,` between them, so a
#    high-precedence AND written immediately after one binds to its LAST alternative only.  `!^;@`
#    is right and `!^&@` is not, exactly as for `~&@` -- write the parenthesis-free form with `;`.


with cython.warn.undeclared(False):
    # bare so Python can import it, guarded so warn.undeclared stays quiet
    class IncorrectSmarts(IncorrectSmiles):
        """The string is not a SMARTS: the reader could not decide what query it names.

        Raised for SYNTAX only -- a letter that names no primitive, an unbalanced bracket, a ring
        label that never closes, a value outside what the box layout can hold.  The message ends in
        a byte offset into the input.  Subclasses `IncorrectSmiles` because the two readers make the
        same promise and a pipeline should not have to catch both.
        """


# "no atom" in the chain.  A typed global rather than a DEF, for the reason SMI_NONE has one.
cdef uint32_t SMA_NONE = 0xFFFFFFFF

# ring-closure labels `0`-`9` and `%00`-`%99`, as in SMILES and for the same reason
DEF SMA_CLOSURES = 100

# Bond-expression tokens the pool must hold per input character.  `~` is the worst case: five
# primitives, the four ORs between them, and the one operator that may stand in front of the whole
# thing -- ten.  Every other character contributes at most one primitive and one operator.  Twelve is
# a proven bound rather than an estimate, which is what removes the growth path from `sma_alloc`;
# `sma_push` checks it anyway, because a bound nobody tests is a bound nobody maintains.
DEF SMA_TOKS_PER_CHAR = 12

# the operator waiting in front of a primitive that has not been read yet
cdef enum:
    SMA_OP_NONE = 0
    SMA_OP_AND_LOW = 1
    SMA_OP_OR = 2
    SMA_OP_AND_HIGH = 3


# The primitive and operator names `QueryContainer` takes.  Bound once at module level rather than
# spelled at each call site: `atom_primitive` looks the name up in PRIM_NAMES, and a literal in the
# loop would allocate a str per primitive read.
cdef str SMA_AND_LOW = 'and_low'
cdef str SMA_OR = 'or'
cdef str SMA_AND_HIGH = 'and_high'

cdef str SMA_P_ELEMENT = 'element'
cdef str SMA_P_ANY = 'any'
cdef str SMA_P_METAL = 'metal'
cdef str SMA_P_ISOTOPE = 'isotope'
# Daylight spells "this atom carries NO mass number" as a leading zero, and the box layout has a bit
# for it.
cdef str SMA_P_NO_ISOTOPE = 'no_isotope'
cdef str SMA_P_CHARGE = 'charge'
# `*`: withdraw the charge default instead of demanding a value.  One box; the alternative spelling is
# a thirteen-way OR over the whole span, which is what the chemistry layer builds today.
cdef str SMA_P_ANY_CHARGE = 'any_charge'
cdef str SMA_P_RADICAL = 'radical'
cdef str SMA_P_DEGREE = 'degree'
cdef str SMA_P_IMPLICIT_H = 'implicit_h'
cdef str SMA_P_TOTAL_H = 'total_h'
cdef str SMA_P_HETEROATOMS = 'heteroatoms'
cdef str SMA_P_HYBRIDIZATION = 'hybridization'
cdef str SMA_P_RING_SIZE = 'ring_size'
cdef str SMA_P_RING_COUNT = 'ring_count'
cdef str SMA_P_STEREO = 'stereo'
# `#0`: the R marker, element 0.  A patch spelling and not a test -- the seal refuses it, because an R
# matches nothing.
cdef str SMA_P_R_MARKER = 'r_marker'
# `@=`: the configuration the reactant had.  A patch spelling and not a test, for the same reason -- a
# query has no reactant to have had one.
cdef str SMA_P_STEREO_KEEP = 'stereo_keep'
cdef str SMA_P_STEREO_INVERT = 'stereo_invert'
cdef str SMA_P_BOND_ORDER = 'bond_order'
cdef str SMA_P_BOND_AROMATIC = 'bond_aromatic'
cdef str SMA_P_BOND_RING = 'bond_ring'


cdef struct sma_tok_t:
    uint8_t opcode         # SMA_OP_* for an operator, 0xFF for a primitive
    uint8_t kind           # BPRIM_ORDER / BPRIM_AROMATIC / BPRIM_RING
    int32_t value
    uint8_t negated


# a primitive rather than an operator, in sma_tok_t.opcode
DEF SMA_TOK_PRIM = 0xFF

# A `/` or `\` in sma_tok_t.kind, beside the BPRIM_* the box compiler knows.  Private to this file
# because it never reaches the compiler: `sma_emit_bond` turns it into the single bond it carries plus
# a `set_bond_direction` label, so what the box compiler sees is an ordinary order primitive.
DEF SMA_BPRIM_DIRECTION = 0xFE


cdef struct sma_parse_t:
    void *block
    sma_tok_t *btoks       # bond-expression pool; expressions are never reused, so never freed
    uint32_t *stack        # branch stack: the atom each `(` returns to
    uint8_t *low           # per STABLE ID: the atom was written lowercase
    uint32_t n_btoks
    uint32_t cap_btoks
    uint32_t depth
    uint32_t n_atoms
    int32_t group          # the component group being read, or -1 outside one
    int32_t next_group     # the next group number to hand out; groups are 0-based, as the arena's are
    uint32_t group_at      # byte offset of the `(` that opened `group`, for the error message
    uint32_t group_n0      # n_atoms when it opened, so an empty group can be told from a full one
    uint32_t open_atom[SMA_CLOSURES]
    uint32_t open_at[SMA_CLOSURES]     # byte offset of the opening label, for the error message
    uint32_t open_off[SMA_CLOSURES]    # its bond expression in the pool
    uint32_t open_count[SMA_CLOSURES]
    uint8_t open_live[SMA_CLOSURES]


cdef int sma_alloc(sma_parse_t *p, uint32_t n) except -1:
    """One struct, one malloc, one check, one free -- RULES.md 5.2, and no realloc at all.

    Every array is bounded by the input length: an atom needs at least one character, a branch needs
    its `(`, and a bond expression's token count is bounded per character by SMA_TOKS_PER_CHAR.
    """
    cdef size_t cap = <size_t> n + 2
    cdef size_t btoks_len = align8(SMA_TOKS_PER_CHAR * cap * sizeof(sma_tok_t))
    cdef size_t stack_len = align8(cap * sizeof(uint32_t))
    cdef size_t low_len = align8(cap * sizeof(uint8_t))
    cdef size_t total = btoks_len + stack_len + low_len
    cdef char *block = <char *> PyMem_Malloc(total)
    cdef size_t off = 0
    cdef uint32_t i
    if block is NULL:
        raise MemoryError('SMARTS parse allocation failed')
    memset(block, 0, total)
    p.block = <void *> block
    p.btoks = <sma_tok_t *> (block + off); off += btoks_len
    p.stack = <uint32_t *> (block + off); off += stack_len
    p.low = <uint8_t *> (block + off)
    p.n_btoks = 0
    p.cap_btoks = <uint32_t> (SMA_TOKS_PER_CHAR * cap)
    p.depth = 0
    p.n_atoms = 0
    p.group = -1
    p.next_group = 0
    p.group_at = 0
    p.group_n0 = 0
    for i in range(SMA_CLOSURES):
        p.open_live[i] = 0
    return 0


cdef void sma_free(sma_parse_t *p) noexcept:
    PyMem_Free(p.block)
    p.block = NULL


cdef inline uint32_t sma_bracket_element(const char *s, uint32_t n, uint32_t *ip,
                                         bint *lower) noexcept nogil:
    """The element symbol at `ip`, or 0 -- without moving `ip` -- when this is not one.

    Differs from `smi_bracket_symbol` in exactly one place, and the file header explains why: the ten
    letters chython SMARTS spends on primitives never begin a ONE-letter symbol.  The two-letter
    lookup runs first regardless, so `[Dy]`, `[Ho]`, `[Mg]`, `[Ag]` and `[Ru]` keep their elements.
    """
    cdef uint32_t i = ip[0]
    cdef char c0 = s[i]
    cdef char c1
    cdef uint32_t a, z
    if smi_islower(c0):
        # a h r x z: the primitive wins, and the two-letter lookup is not even tried -- otherwise
        # `[C;as]` would be aromatic arsenic instead of the nonsense it is, and, far worse, a future
        # `[C;ru]` would silently stop being a ring-size demand
        if c0 == 97 or c0 == 104 or c0 == 114 or c0 == 120 or c0 == 122:
            return 0
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
    if c0 == 65 or c0 == 68 or c0 == 72 or c0 == 77 or c0 == 82:      # A D H M R
        return 0
    z = SMI_ELEMENT[a * 27]
    if z:
        ip[0] = i + 1
        return z
    return 0


cdef struct sma_emit_t:
    uint32_t sid
    uint32_t isotope       # stated by the leading digits; 0 = none
    uint8_t pend           # SMA_OP_*, honoured only once a token has been emitted
    uint8_t any_tok
    uint8_t saw_element


cdef int sma_atom_tok(QueryContainer q, sma_emit_t *e, str name, int32_t value,
                      bint negated) except -1:
    """Emit one atom primitive, flushing the operator that was waiting in front of it.

    The operator is flushed HERE rather than where it was read, because several spellings emit no
    token at all -- `M` sets a flag, `:12` sets a number, a trailing `A` says nothing -- and a `;` in
    front of one of those would leave a dangling operator that `compile_term` rejects as a malformed
    stream.  Waiting until a token actually arrives is also what makes a leading `[;C]` harmless.
    """
    if e.any_tok:
        if e.pend == SMA_OP_AND_LOW:
            q.atom_operator(e.sid, SMA_AND_LOW)
        elif e.pend == SMA_OP_OR:
            q.atom_operator(e.sid, SMA_OR)
        elif e.pend == SMA_OP_AND_HIGH:
            q.atom_operator(e.sid, SMA_AND_HIGH)
    e.pend = SMA_OP_NONE
    e.any_tok = 1
    q.atom_primitive(e.sid, name, value, negated)
    return 0


cdef int sma_atom_element(QueryContainer q, sma_emit_t *e, uint32_t z, bint lower,
                          bint negated, uint32_t at) except -1:
    """An element primitive, plus the two things that must ride in the SAME box as it.

    The isotope, because `prim_apply` resolves a mass number into an offset from the element's common
    isotope and so refuses one in a box with no settled element -- which is why a leading `13` is held
    until here instead of being emitted where it was read.  The aromatic flag, because a lowercase
    symbol states hybridization about that element and not about the box's other alternatives.
    """
    cdef object exc                      # the except-as target; `warn.undeclared` counts it
    sma_atom_tok(q, e, SMA_P_ELEMENT, <int32_t> z, negated)
    e.saw_element = 1
    if e.isotope:
        q.atom_operator(e.sid, SMA_AND_HIGH)
        try:
            q.atom_primitive(e.sid, SMA_P_ISOTOPE, <int32_t> e.isotope, False)
        except ValueError as exc:
            raise IncorrectSmarts('%s, at position %d' % (exc, at))
    if lower:
        q.atom_operator(e.sid, SMA_AND_HIGH)
        q.atom_primitive(e.sid, SMA_P_HYBRIDIZATION, 4, False)
    return 0


cdef inline uint32_t sma_count(const char *s, uint32_t n, uint32_t *ip, uint32_t limit,
                               bint *seen) noexcept nogil:
    """A run of at most `limit` digits, with `seen` saying whether there were any.

    `smi_digits` returns 0 both for `D0` and for a `D` with nothing after it, and those are a demand
    and a syntax error; the flag is what separates them without the caller re-reading `ip`.
    """
    cdef uint32_t before = ip[0]
    cdef uint32_t value = smi_digits(s, n, ip, limit)
    seen[0] = ip[0] != before
    return value


cdef int sma_ring(QueryContainer q, sma_emit_t *e, int32_t value, bint negated,
                  uint32_t at) except -1:
    """`r<n>`: a ring-size demand, except that size zero is a ring-COUNT demand.

    The core keeps the two questions apart -- `ring_size` is a multi-hot span over sizes 3 and up,
    `ring_count` a one-hot span over 0..8 -- while "acyclic" is written as the ring size zero (V2's
    `!R` sets `ring_sizes = (0,)`, one set answering both questions), so that spelling is translated
    here rather than pushed onto every caller.  Sizes 1 and 2 have no bit and no meaning: they are
    refused.
    """
    if value == 0:
        sma_atom_tok(q, e, SMA_P_RING_COUNT, 0, negated)
    elif value < 3:
        raise IncorrectSmarts('a ring of size %d does not exist; sizes start at 3, at position %d'
                              % (value, at))
    else:
        sma_atom_tok(q, e, SMA_P_RING_SIZE, value, negated)
    return 0


cdef int sma_bracket(QueryContainer q, const char *s, uint32_t n, uint32_t *ip, uint32_t sid,
                     uint8_t *arom, object log) except -1:
    """Parse a bracket atom, emitting its tokens for stable id `sid`.

    `s[ip[0]]` is the `[`; on return `ip[0]` is past the `]`.  Fields are read in a LOOP keyed on the
    leading character and in any order, which is what makes `[13C@:7]`, `[C:7@13]` and `[N;+2;D3]`
    all read, with no field deleted from the body before the rest of it is split.

    `arom` reports whether the element symbol was WRITTEN lowercase, which is what the implicit-bond
    rule keys on.  `[N;a]` is an aromatic atom and sets it to 0 on purpose: a template that spells
    the atom uppercase spells its bonds `:` too, and 185 `a` primitives in this codebase depend on
    that reading.
    """
    cdef uint32_t start = ip[0]
    cdef uint32_t i = start + 1
    cdef uint32_t at = i
    cdef uint32_t value = 0
    cdef uint32_t z
    cdef int sign
    cdef char c
    cdef bint negated
    cdef bint seen = False
    cdef bint lower = False
    cdef bint first = True
    cdef bint sgroup
    cdef sma_emit_t e

    e.sid = sid
    e.isotope = 0
    e.pend = SMA_OP_NONE
    e.any_tok = 0
    e.saw_element = 0
    arom[0] = 0

    if i >= n:
        raise IncorrectSmarts('unterminated bracket atom at position %d' % start)
    if smi_isdigit(s[i]):
        value = smi_digits(s, n, &i, 5)
        if value == 0:
            # `[0C]` is Daylight's "no mass number stated", and unlike a real isotope it needs no
            # settled element in the box -- PRIM_NO_ISOTOPE forbids one span bit and nothing else --
            # so it is emitted here instead of being held until the element arrives.
            sma_atom_tok(q, &e, SMA_P_NO_ISOTOPE, 0, False)
        elif value > ISOTOPE_MAX:
            raise IncorrectSmarts('isotope %d is above the storable %d at position %d'
                                  % (value, ISOTOPE_MAX, at))
        else:
            e.isotope = value

    while i < n:
        c = s[i]
        if c == 93:                                  # ']'
            # the isotope check FIRST: `[13]` has no token either, and "states no primitive" would
            # send the writer looking for a missing primitive instead of the element the digits need
            if e.isotope and not e.saw_element:
                raise IncorrectSmarts('the isotope at position %d has no element symbol to attach '
                                      'to' % (start + 1))
            if not e.any_tok:
                raise IncorrectSmarts('the bracket at position %d states no primitive' % start)
            ip[0] = i + 1
            return 0
        elif c == 59:                                # ';'
            e.pend = SMA_OP_AND_LOW
            i += 1
            continue
        elif c == 44:                                # ','
            e.pend = SMA_OP_OR
            i += 1
            continue
        elif c == 38 and not (i + 1 < n and smi_isdigit(s[i + 1])):
            # '&', Daylight's high AND.  A DIGIT AFTER IT IS NOT THIS: `&1` is an enhanced-stereo AND
            # group, handled below with `o1`.  One character of lookahead settles it, because no
            # primitive name is a digit -- `[C&1]` has nothing else it could mean, and `[C&D1]` is
            # untouched.
            e.pend = SMA_OP_AND_HIGH
            i += 1
            continue

        # from here on this is a primitive.  Arriving with nothing pending is juxtaposition, which
        # binds tighter than `,`: `[N+,O]` is (N and +1) or O, not N and (+1 or O)
        if e.any_tok and e.pend == SMA_OP_NONE:
            e.pend = SMA_OP_AND_HIGH
        at = i
        negated = False
        if c == 33:                                  # '!'
            negated = True
            i += 1
            if i >= n:
                raise IncorrectSmarts('unterminated bracket atom at position %d' % start)
            c = s[i]

        # `&<n>` AND, `o<n>` OR: THE ENHANCED-STEREO GROUP, IN THE BRACKET.  CXSMILES spells the same
        # two kinds in a `|...|` tail, by zero-based index over the atoms as written -- fine for a
        # molecule serialized once, and miserable for a template, where the atoms ARE the pattern and
        # a reaction's tail indexes both sides end to end.  Written here it names one atom on one
        # side and needs no counting.  The lookahead is what keeps `[o]` aromatic oxygen and `[C&D1]`
        # a high AND; the cost is that an aromatic oxygen in a group has to be spelled `[o;o1]`.
        sgroup = c == 38 or (c == 111 and i + 1 < n and smi_isdigit(s[i + 1]))
        z = 0 if sgroup else sma_bracket_element(s, n, &i, &lower)
        if sgroup:
            if negated:
                raise IncorrectSmarts('an enhanced-stereo group is a label on the atom, not a test, '
                                      'and cannot be negated, at position %d' % at)
            i += 1
            value = sma_count(s, n, &i, 2, &seen)
            # `seen` cannot be false -- the lookahead above found the digit
            if value < 1 or value > 63:
                raise IncorrectSmarts('enhanced-stereo group %d is outside the storable 1..63 at '
                                      'position %d' % (value, at))
            q.set_stereo_group(sid, 3 if c == 38 else 2, <int> value)
        elif z:
            sma_atom_element(q, &e, z, lower, negated, at)
            if lower:
                arom[0] = 1
        elif first and c == 72:                      # 'H' leading the body is hydrogen ITSELF
            i += 1
            sma_atom_element(q, &e, 1, False, negated, at)
        elif first and c == 65:                      # 'A' leading the body is any atom
            i += 1
            sma_atom_tok(q, &e, SMA_P_ANY, 0, negated)
        elif first and c == 77:                      # 'M' leading the body is any metal
            i += 1
            sma_atom_tok(q, &e, SMA_P_METAL, 0, negated)
        elif c == 35:                                # '#<n>', the atomic number
            i += 1
            value = sma_count(s, n, &i, 3, &seen)
            if not seen:
                raise IncorrectSmarts('`#` needs an atomic number at position %d' % at)
            if value > 118:
                raise IncorrectSmarts('atomic number %d is outside 0..118 at position %d'
                                      % (value, at))
            if value:
                sma_atom_element(q, &e, value, False, negated, at)
            else:
                # `#0` IS THE R MARKER, and it is not an element: it takes no isotope and sets no
                # `saw_element`, so `[13#0]` is refused where `[13C]` is read.  A SMIRKS product side
                # builds one; every other reader of this token seals, and the seal refuses it.
                sma_atom_tok(q, &e, SMA_P_R_MARKER, 0, negated)
        elif c == 42:                                # '*', UNCONSTRAINED -- element and charge both
            i += 1
            if negated:
                raise IncorrectSmarts('`!*` at position %d is not a demand: `*` withdraws a default '
                                      'rather than stating a value.  Spell what you mean' % at)
            # The element half only bites when nothing else settles the element, because PRIM_ANY
            # forbids no element bit: `[*]` is any element, `[C;*]` is still carbon.  The charge half
            # always bites, and that is what makes `*` one token wherever it is written -- `[M;*]` is
            # a metal of any charge, and a lone `[*]` is the truly wild wildcard, which is what the
            # rule tables mean by it.  `[A]` is the NEUTRAL any-atom wildcard.
            #
            # The radical default is NOT withdrawn.  `^` states a radical and `!^` states none, so all
            # three readings stay spellable per atom; `*` would have to guess which was meant.
            sma_atom_tok(q, &e, SMA_P_ANY, 0, False)
            q.atom_operator(e.sid, SMA_AND_HIGH)
            q.atom_primitive(e.sid, SMA_P_ANY_CHARGE, 0, False)
        elif c == 68:                                # 'D', heavy-atom degree
            i += 1
            value = sma_count(s, n, &i, 2, &seen)
            if not seen:
                raise IncorrectSmarts('`D` needs a degree at position %d' % at)
            sma_atom_tok(q, &e, SMA_P_DEGREE, <int32_t> value, negated)
        elif c == 104:                               # 'h', implicit hydrogens
            i += 1
            value = sma_count(s, n, &i, 2, &seen)
            if not seen:
                raise IncorrectSmarts('`h` needs a hydrogen count at position %d' % at)
            sma_atom_tok(q, &e, SMA_P_IMPLICIT_H, <int32_t> value, negated)
        elif c == 72:                                # 'H' later in the body: TOTAL hydrogens
            i += 1
            value = sma_count(s, n, &i, 2, &seen)
            if not seen:
                value = 1                            # `[CH]` is one hydrogen, as everywhere else
            sma_atom_tok(q, &e, SMA_P_TOTAL_H, <int32_t> value, negated)
        elif c == 120:                               # 'x', heteroatom neighbours
            i += 1
            value = sma_count(s, n, &i, 2, &seen)
            if not seen:
                raise IncorrectSmarts('`x` needs a heteroatom count at position %d' % at)
            sma_atom_tok(q, &e, SMA_P_HETEROATOMS, <int32_t> value, negated)
        elif c == 122:                               # 'z', hybridization
            i += 1
            value = sma_count(s, n, &i, 1, &seen)
            if not seen:
                raise IncorrectSmarts('`z` needs a hybridization at position %d' % at)
            if value < 1 or value > 6:
                # 1..6 and NOT V2's 1..3: the core reports what it found instead of saturating, so 5
                # is two cumulated doubles and 6 is anything past that.  The number goes through
                # unchanged -- translating a template written against V2 is the template's business
                raise IncorrectSmarts('hybridization %d is outside 1..6 at position %d'
                                      % (value, at))
            sma_atom_tok(q, &e, SMA_P_HYBRIDIZATION, <int32_t> value, negated)
        elif c == 97:                                # 'a', aromatic
            i += 1
            sma_atom_tok(q, &e, SMA_P_HYBRIDIZATION, 4, negated)
        elif c == 114:                               # 'r', ring size
            i += 1
            value = sma_count(s, n, &i, 2, &seen)
            if not seen:
                raise IncorrectSmarts('`r` needs a ring size at position %d' % at)
            sma_ring(q, &e, <int32_t> value, negated, at)
        elif c == 82:                                # 'R', ring count; bare `R` is "in some ring"
            i += 1
            value = sma_count(s, n, &i, 2, &seen)
            if not seen:
                # `!R` is the only spelling of either question and means acyclic, so bare `R` has to
                # be its negation for the `!` above to land on the right side of it
                sma_atom_tok(q, &e, SMA_P_RING_COUNT, 0, not negated)
            else:
                sma_atom_tok(q, &e, SMA_P_RING_COUNT, <int32_t> value, negated)
        elif c == 65:                                # 'A' later in the body: accepted and ignored
            i += 1
        elif c == 77:                                # 'M' later in the body: the masked flag
            i += 1
            if negated:
                raise IncorrectSmarts('`M` is a flag on the atom, not a test, and cannot be '
                                      'negated, at position %d' % at)
            q.set_masked(sid)
        elif c == 43 or c == 45:                     # '+', '-'
            sign = 1 if c == 43 else -1
            i += 1
            value = sma_count(s, n, &i, 2, &seen)
            if not seen:
                value = 1
                while i < n and s[i] == c:           # the repeated-sign spelling: `[Fe++]`
                    value += 1
                    i += 1
            if sign * <int> value < CHARGE_MIN or sign * <int> value > CHARGE_MAX:
                raise IncorrectSmarts('charge %d is outside the storable %d..%d at position %d'
                                      % (sign * <int> value, CHARGE_MIN, CHARGE_MAX, at))
            sma_atom_tok(q, &e, SMA_P_CHARGE, sign * <int32_t> value, negated)
        elif c == 94:                                # '^', THIS ATOM IS A RADICAL
            # The same character as the `|^1:idx|` tail spells it with, and the same primitive -- the
            # tail addresses an atom by INDEX, which is fine for a molecule written out once and
            # miserable in a query where the atoms are the pattern.  `[N;D3;^]` says it in place.
            #
            # No collision with the dative bond, which is also `^`: a bond token is lexed BETWEEN
            # atoms and this loop only ever runs inside a bracket.  `!^` is the explicit spelling of
            # what an unstated radical already means, and it is here because negation is generic, not
            # because anything needs it -- though it is what makes `^,!^` say "radical or not".
            i += 1
            sma_atom_tok(q, &e, SMA_P_RADICAL, 0, negated)
        elif c == 64:                                # '@', a configuration
            i += 1
            if i < n and s[i] == 64:
                i += 1
                value = 2
            elif i < n and s[i] == 63:               # `@?`: declared unresolved
                i += 1
                value = 0
            elif i < n and s[i] == 61:               # `@=`: THE CONFIGURATION THE REACTANT HAD
                # A SMIRKS product-side token: `@=` says the configuration comes through the reaction
                # unchanged whatever it was, which is what a template needs where the reaction centre
                # is the very atom carrying it.  A patch drops the configuration at its reaction
                # centre otherwise, so this is the retention statement, and `@~` below is its twin.
                #
                # No collision with the bond order, which is also `=`: a bond token is lexed BETWEEN
                # atoms and this loop only ever runs inside a bracket, exactly as for `^`.
                i += 1
                sma_atom_tok(q, &e, SMA_P_STEREO_KEEP, 0, negated)
                first = False
                continue
            elif i < n and s[i] == 126:              # `@~`: THE OTHER CONFIGURATION
                # The inversion statement, and product-side like `@=`.  Every unit chython models has
                # exactly two states -- SU_TETRA, SU_CIS_TRANS, SU_ALLENE, SU_ATROPISOMER -- so "the
                # other one" is well defined for all four, and this one token covers a parity and a
                # geometry alike.  Where the reactant held no configuration it does nothing, which is
                # what makes it usable in a template that does not narrow its match to configured
                # atoms: a sign on the reactant side is selectivity and nothing here needs it.
                #
                # `~` outside a bracket stays "any bond" (decision 1); a bond token is never lexed
                # inside one, which is the same non-collision argument `@=` and `^` rest on.
                i += 1
                sma_atom_tok(q, &e, SMA_P_STEREO_INVERT, 0, negated)
                first = False
                continue
            else:
                value = 1
            if not value:
                log.append(mc_record('smarts:unresolved-stereo-dropped', (),
                                     'atom %d: `@?` declares a configuration nobody resolved, which is not '
                                     'something a query can test; it was dropped' % sid,
                                     mc_lost()))
            elif negated:
                # the negation of "even in the query's own frame" is "odd OR not configured at all",
                # and the second half is what leaving the primitive out already says
                raise IncorrectSmarts('a configuration cannot be negated: write the other one, or '
                                      'nothing at all, at position %d' % at)
            else:
                sma_atom_tok(q, &e, SMA_P_STEREO, <int32_t> value, False)
        elif c == 58:                                # ':<n>', the atom map
            i += 1
            value = sma_count(s, n, &i, 5, &seen)
            if not seen:
                raise IncorrectSmarts('atom map `:` with no number at position %d' % at)
            if value > MAP_NUMBER_MAX:
                raise IncorrectSmarts('atom map %d is above the storable %d at position %d'
                                      % (value, MAP_NUMBER_MAX, at))
            if negated:
                raise IncorrectSmarts('an atom map is a label on the atom, not a test, and cannot '
                                      'be negated, at position %d' % at)
            q.set_map_number(sid, <int> value)
        else:
            raise IncorrectSmarts('%r names no primitive inside a bracket atom, at position %d'
                                  % (s[i:i + 1].decode('ascii', 'replace'), i))
        first = False
    raise IncorrectSmarts('unterminated bracket atom at position %d' % start)


cdef inline int sma_push(sma_parse_t *p, uint8_t opcode, uint8_t kind, int32_t value,
                         bint negated) except -1:
    """Append one bond-expression token to the pool."""
    cdef sma_tok_t *t
    if p.n_btoks >= p.cap_btoks:
        raise ValueError('the SMARTS bond token pool overflowed; SMA_TOKS_PER_CHAR is too small')
    t = p.btoks + p.n_btoks
    t.opcode = opcode
    t.kind = kind
    t.value = value
    t.negated = 1 if negated else 0
    p.n_btoks += 1
    return 0


cdef inline bint sma_is_bond(char c) noexcept nogil:
    """Does this character begin a bond primitive?  `!` counts: it is a bond unit's first byte."""
    return (c == 45 or c == 61 or c == 35 or c == 58 or c == 126 or c == 64 or c == 33
            or c == 94 or c == 47 or c == 92)


cdef int sma_bond_expr(sma_parse_t *p, const char *s, uint32_t n, uint32_t *ip,
                       uint32_t *off, uint32_t *cnt) except -1:
    """Parse one bond expression into the pool, returning its `(offset, count)`.

    Held rather than emitted because neither end of the bond is necessarily known yet: a chain bond
    is read before its second atom and a ring closure's opening label before its partner.

    `/` and `\\` are the one token here that is not a test of the bond it sits on.  They state which
    SIDE of a double bond a substituent is on, so they carry the bond's order as well (single, as in
    SMILES) and stand ALONE: negating a side or or-ing two of them describes neither a geometry to
    build nor one to look for.  A pair of them becomes a geometry -- `_seal_geometries` on a query
    side, `smk_directions` on a product side -- and one alone is refused there.
    """
    cdef uint32_t i = ip[0]
    cdef uint32_t start = p.n_btoks
    cdef uint32_t at = i
    cdef uint32_t k
    cdef char c
    cdef bint negated
    cdef bint want = True

    while i < n:
        c = s[i]
        if not want:
            if c == 59:                              # ';'
                sma_push(p, SMA_OP_AND_LOW, 0, 0, False)
            elif c == 44:                            # ','
                sma_push(p, SMA_OP_OR, 0, 0, False)
            elif c == 38:                            # '&'
                sma_push(p, SMA_OP_AND_HIGH, 0, 0, False)
            elif sma_is_bond(c):
                # juxtaposition, which is the high AND here as it is inside a bracket.  The character
                # is NOT consumed: the loop comes back round and reads it as the next unit.
                sma_push(p, SMA_OP_AND_HIGH, 0, 0, False)
                want = True
                continue
            else:
                break                                # the expression ends where the next atom starts
            i += 1
            want = True
            continue

        at = i
        negated = False
        if c == 33:                                  # '!'
            negated = True
            i += 1
            if i >= n:
                raise IncorrectSmarts('a bond expression ends the string at position %d' % at)
            c = s[i]

        if c == 45:                                  # '-'
            sma_push(p, SMA_TOK_PRIM, BPRIM_ORDER, 1, negated)
        elif c == 61:                                # '='
            sma_push(p, SMA_TOK_PRIM, BPRIM_ORDER, 2, negated)
        elif c == 35:                                # '#'
            sma_push(p, SMA_TOK_PRIM, BPRIM_ORDER, 3, negated)
        elif c == 58:                                # ':'
            sma_push(p, SMA_TOK_PRIM, BPRIM_AROMATIC, 0, negated)
        elif c == 64:                                # '@', in a ring
            sma_push(p, SMA_TOK_PRIM, BPRIM_RING, 0, negated)
        elif c == 94:                                # '^', THE DATIVE BOND -- decision 6
            if negated:
                # "NOT a coordination bond" is the one demand a box cannot state: the exact set is
                # {1, 2, 3, aromatic}, which over the word-0 order bits is a DISJUNCTION, and a box
                # is a conjunction of forbidden bits (`_query_boxes.pxi` refuses it with the same
                # reasoning and tells the caller to spell the orders instead).  So spell them --
                # `!^` is `-,=,#,:`, built here so the caller never has to.
                sma_push(p, SMA_TOK_PRIM, BPRIM_ORDER, 1, False)
                sma_push(p, SMA_OP_OR, 0, 0, False)
                sma_push(p, SMA_TOK_PRIM, BPRIM_ORDER, 2, False)
                sma_push(p, SMA_OP_OR, 0, 0, False)
                sma_push(p, SMA_TOK_PRIM, BPRIM_ORDER, 3, False)
                sma_push(p, SMA_OP_OR, 0, 0, False)
                sma_push(p, SMA_TOK_PRIM, BPRIM_AROMATIC, 0, False)
            else:
                sma_push(p, SMA_TOK_PRIM, BPRIM_ORDER, 8, False)
        elif c == 126:                               # '~', ANY bond -- see decision 1 in the header
            if negated:
                raise IncorrectSmarts('`!~` forbids every bond there is and can never match, at '
                                      'position %d' % at)
            sma_push(p, SMA_TOK_PRIM, BPRIM_ORDER, 1, False)
            sma_push(p, SMA_OP_OR, 0, 0, False)
            sma_push(p, SMA_TOK_PRIM, BPRIM_ORDER, 2, False)
            sma_push(p, SMA_OP_OR, 0, 0, False)
            sma_push(p, SMA_TOK_PRIM, BPRIM_ORDER, 3, False)
            sma_push(p, SMA_OP_OR, 0, 0, False)
            sma_push(p, SMA_TOK_PRIM, BPRIM_AROMATIC, 0, False)
            sma_push(p, SMA_OP_OR, 0, 0, False)
            sma_push(p, SMA_TOK_PRIM, BPRIM_ORDER, 8, False)
        elif c == 47 or c == 92:                     # '/' and '\', a SIDE rather than a test
            if negated:
                raise IncorrectSmarts('`!/` and `!\\` name every side but one, of which there is one, '
                                      'at position %d' % at)
            sma_push(p, SMA_TOK_PRIM, SMA_BPRIM_DIRECTION,
                     SMI_DIR_UP if c == 47 else SMI_DIR_DOWN, False)
        else:
            raise IncorrectSmarts('%r names no bond at position %d'
                                  % (s[i:i + 1].decode('ascii', 'replace'), i))
        i += 1
        want = False

    if want:
        raise IncorrectSmarts('a bond expression ends with an operator at position %d' % at)
    cnt[0] = p.n_btoks - start
    if cnt[0] != 1:
        for k in range(start, p.n_btoks):
            if p.btoks[k].opcode == SMA_TOK_PRIM and p.btoks[k].kind == SMA_BPRIM_DIRECTION:
                raise IncorrectSmarts('a `/` or `\\` states a side and combines with nothing, at '
                                      'position %d; it carries the single bond itself' % at)
    ip[0] = i
    off[0] = start
    return 0


cdef int sma_emit_bond(QueryContainer q, sma_parse_t *p, uint32_t u, uint32_t v,
                       uint32_t off, uint32_t cnt) except -1:
    """Replay a held bond expression onto the bond `(u, v)`, which must already exist.

    `u` is the atom written FIRST, which is the only thing a direction needs beyond the two ids: `/`
    means "up, going from `u` to `v`", and the same statement read from `v` is `\\`.
    """
    cdef uint32_t k
    cdef sma_tok_t *t
    if not cnt:
        # An untokenised bond is single, EXCEPT between two atoms both written lowercase, where it is
        # aromatic for the reason it is in SMILES.  `query_seal` supplies the single itself when a
        # bond carries no token at all, so the common case costs nothing at all here.
        if p.low[u] and p.low[v]:
            q.bond_primitive(u, v, SMA_P_BOND_AROMATIC, 0, False)
        return 0
    for k in range(off, off + cnt):
        t = p.btoks + k
        if t.opcode == SMA_TOK_PRIM:
            if t.kind == SMA_BPRIM_DIRECTION:
                # the single bond the token carries, then the side as a label; a lone expression by
                # construction, so the operator stream gets exactly one operand out of this
                q.bond_primitive(u, v, SMA_P_BOND_ORDER, 1, False)
                q.set_bond_direction(u, v, t.value)
            elif t.kind == BPRIM_ORDER:
                q.bond_primitive(u, v, SMA_P_BOND_ORDER, t.value, t.negated)
            elif t.kind == BPRIM_AROMATIC:
                q.bond_primitive(u, v, SMA_P_BOND_AROMATIC, t.value, t.negated)
            else:
                q.bond_primitive(u, v, SMA_P_BOND_RING, t.value, t.negated)
        elif t.opcode == SMA_OP_AND_LOW:
            q.bond_operator(u, v, SMA_AND_LOW)
        elif t.opcode == SMA_OP_OR:
            q.bond_operator(u, v, SMA_OR)
        else:
            q.bond_operator(u, v, SMA_AND_HIGH)
    return 0


cdef inline bint sma_same_expr(sma_parse_t *p, uint32_t off1, uint32_t cnt1,
                               uint32_t off2, uint32_t cnt2) noexcept nogil:
    """Do two held bond expressions say the same thing, token for token?"""
    if cnt1 != cnt2:
        return False
    return memcmp(p.btoks + off1, p.btoks + off2, cnt1 * sizeof(sma_tok_t)) == 0


cdef inline uint32_t sma_new_atom(QueryContainer q, sma_parse_t *p, bint lower) except 0:
    """Take the next stable id.  `add_atom` hands them out in reading order from 1, which is what
    lets the CXSMARTS tail turn a 0-based atom index into a stable id by adding one.

    An atom read inside a parenthesised component position is assigned that group here, which is the
    only place it can be: the seal reduces a per-atom group to a per-component one and refuses a
    component that spans two, so every atom of a grouped component has to carry the number.
    """
    cdef uint32_t sid = q.add_atom()
    p.n_atoms += 1
    p.low[sid] = 1 if lower else 0
    if p.group >= 0:
        q.set_group(sid, p.group)
    return sid


cdef int sma_tokenize(QueryContainer q, sma_parse_t *p, const char *s, uint32_t n,
                      object log) except -1:
    """One pass over the bytes, emitting journal ops.  Raises `IncorrectSmarts` on syntax."""
    cdef uint32_t i = 0
    cdef uint32_t j = 0
    cdef uint32_t label, u, sid
    cdef uint32_t prev = SMA_NONE
    cdef uint32_t pend_off = 0
    cdef uint32_t pend_count = 0
    cdef uint32_t pend_at = 0
    cdef uint32_t off = 0
    cdef uint32_t cnt = 0
    cdef bint pend_bond = False
    cdef bint lower = False
    cdef char c
    cdef uint32_t z
    cdef uint8_t arom = 0

    # every branch that is not an atom `continue`s; falling out of the chain means `sid` is a new atom
    # waiting to be joined to `prev`, which is the one place that joining happens
    while i < n:
        c = s[i]
        if c == 91:                                  # '['
            sid = sma_new_atom(q, p, False)
            sma_bracket(q, s, n, &i, sid, &arom, log)
            p.low[sid] = arom
        elif c == 40:                                # '(' -- a branch, or a component group
            if pend_bond:
                raise IncorrectSmarts('bond expression immediately before `(` at position %d' % i)
            if prev == SMA_NONE and not p.depth:
                # A COMPONENT POSITION, so this parenthesis groups components rather than opening a
                # branch: `(A.B)>>...` demands A and B in ONE molecule, `(A).(B)` in two.  Legal only
                # here -- with no atom to branch from and no branch already open -- which is exactly
                # where a branch cannot be meant.
                if p.group >= 0:
                    raise IncorrectSmarts('component group opens inside the one at position %d, at '
                                          'position %d' % (p.group_at, i))
                p.group = p.next_group
                p.next_group += 1
                p.group_at = i
                p.group_n0 = p.n_atoms
                i += 1
                continue
            if prev == SMA_NONE:
                raise IncorrectSmarts('branch opens before any atom at position %d' % i)
            p.stack[p.depth] = prev
            p.depth += 1
            i += 1
            continue
        elif c == 41:                                # ')'
            if pend_bond:
                raise IncorrectSmarts('bond expression immediately before `)` at position %d' % i)
            if not p.depth:
                if p.group < 0:
                    raise IncorrectSmarts('unbalanced `)` at position %d' % i)
                if p.n_atoms == p.group_n0:
                    # keeps the group numbers dense as well as saying what is wrong: the matcher
                    # sizes its per-group anchor table by the highest number it sees
                    raise IncorrectSmarts('component group at position %d holds no atom'
                                          % p.group_at)
                # closing a component group also ends the component: what follows starts a new one,
                # grouped or not, and cannot bond back into this parenthesis
                p.group = -1
                prev = SMA_NONE
                i += 1
                continue
            p.depth -= 1
            prev = p.stack[p.depth]
            i += 1
            continue
        elif c == 46:                                # '.', a component break
            if pend_bond:
                raise IncorrectSmarts('bond expression immediately before `.` at position %d' % i)
            prev = SMA_NONE
            i += 1
            continue
        elif c == 37 or smi_isdigit(c):              # '%NN' or a single-digit ring label
            if c == 37:
                j = i + 1
                label = smi_digits(s, n, &j, 2)
                if j != i + 3:
                    raise IncorrectSmarts('`%%` needs two digits at position %d' % i)
            else:
                label = <uint32_t> (c - 48)
                j = i + 1
            if prev == SMA_NONE:
                raise IncorrectSmarts('ring bond label before any atom at position %d' % i)
            if not p.open_live[label]:
                p.open_live[label] = 1
                p.open_atom[label] = prev
                p.open_at[label] = i
                p.open_off[label] = pend_off
                p.open_count[label] = pend_count
            else:
                u = p.open_atom[label]
                if u == prev:
                    raise IncorrectSmarts('ring bond %d closes on its own atom at position %d'
                                          % (label, i))
                off = p.open_off[label]
                cnt = p.open_count[label]
                # A DIRECTION NEEDS THE ATOM IT WAS WRITTEN FROM, and a ring closure has two
                # candidates: the label's opening end and its closing one, whose statements are each
                # other upside down.  Refused at either end rather than picked, since a template
                # stating a geometry across a ring closure can restate the same bond in the chain.
                if ((cnt == 1 and p.btoks[off].kind == SMA_BPRIM_DIRECTION) or
                        (pend_count == 1 and p.btoks[pend_off].kind == SMA_BPRIM_DIRECTION)):
                    raise IncorrectSmarts('ring bond %d carries a `/` or `\\`, which states a side '
                                          'relative to the atom it is written from and so has two '
                                          'readings on a closure, at position %d' % (label, i))
                if cnt and pend_count:
                    if not sma_same_expr(p, off, cnt, pend_off, pend_count):
                        # Both ends state a bond and they disagree.  The opening statement wins and
                        # the second is reported, because a ring that closes is worth more than a
                        # refusal.
                        log.append(mc_record('smarts:ring-bond-conflict', (),
                                            'ring bond %d states one bond expression where it opens '
                                            '(position %d) and a different one where it closes (position '
                                            '%d); the opening one is kept' % (label, p.open_at[label], i)))
                elif pend_count:
                    off = pend_off
                    cnt = pend_count
                q.add_bond(u, prev)
                sma_emit_bond(q, p, u, prev, off, cnt)
                p.open_live[label] = 0
            pend_bond = False
            pend_count = 0
            i = j
            continue
        elif sma_is_bond(c):
            if pend_bond:
                raise IncorrectSmarts('two bond expressions in a row at position %d' % i)
            pend_at = i
            sma_bond_expr(p, s, n, &i, &pend_off, &pend_count)
            pend_bond = True
            continue
        elif c == 42:                                # '*', unconstrained: any element, any charge
            # Written out here rather than shared with `sma_bracket`'s branch because the two build
            # atoms differently (there is no operator stream to flush yet), but the PAIR of
            # primitives is the point: a bare `*` and a bracketed `[*]` must not drift apart.
            sid = sma_new_atom(q, p, False)
            q.atom_primitive(sid, SMA_P_ANY, 0, False)
            q.atom_operator(sid, SMA_AND_HIGH)
            q.atom_primitive(sid, SMA_P_ANY_CHARGE, 0, False)
            i += 1
        else:
            j = i
            z = smi_bare_symbol(s, n, &j, &lower)
            if z == 0:
                raise IncorrectSmarts('unexpected %r at position %d; a query primitive belongs '
                                      'inside a bracket'
                                      % (s[i:i + 1].decode('ascii', 'replace'), i))
            sid = sma_new_atom(q, p, lower)
            q.atom_primitive(sid, SMA_P_ELEMENT, <int32_t> z, False)
            if lower:
                # decision 2: lowercase outside a bracket is aromatic
                q.atom_operator(sid, SMA_AND_HIGH)
                q.atom_primitive(sid, SMA_P_HYBRIDIZATION, 4, False)
            i = j

        if prev == SMA_NONE:
            if pend_bond:
                raise IncorrectSmarts('a bond expression starts a component at position %d'
                                      % pend_at)
        else:
            q.add_bond(prev, sid)
            sma_emit_bond(q, p, prev, sid, pend_off, pend_count)
        prev = sid
        pend_bond = False
        pend_count = 0

    if p.depth:
        raise IncorrectSmarts('unbalanced `(`: %d branch(es) never close' % p.depth)
    if p.group >= 0:
        raise IncorrectSmarts('component group opens at position %d and never closes' % p.group_at)
    if pend_bond:
        raise IncorrectSmarts('a bond expression ends the string, at position %d' % pend_at)
    for label in range(SMA_CLOSURES):
        if p.open_live[label]:
            raise IncorrectSmarts('ring bond %d opens at position %d and never closes'
                                  % (label, p.open_at[label]))
    if not p.n_atoms:
        raise IncorrectSmarts('no atoms in the string')
    return 0


cdef int sma_cx(QueryContainer q, sma_parse_t *p, bytes block, object log) except -1:
    """The `|...|` tail.  `^N:` radicals are applied; every other field is named in the log.

    Field scanning is `_smiles_read.pxi`'s -- `smi_cx_field_end` and `smi_cx_next_index` -- because a
    tail is a tail and two implementations of "a comma ends a field only when a non-digit follows"
    would drift.  Scanning the whole tail rather than searching it for `^N:` is what lets every other
    field be named instead of passed over.
    """
    cdef const char *s = <const char *> block
    cdef uint32_t n = <uint32_t> len(block)
    cdef uint32_t i = 1
    cdef uint32_t j, k
    cdef uint32_t idx = 0
    cdef uint32_t taken
    if n and s[n - 1] == 124:                        # the closing `|` is not a field
        n -= 1
    while i < n:
        if s[i] == 44:
            i += 1
            continue
        j = smi_cx_field_end(s, n, i)
        if j <= i:
            break
        if s[i] == 94 and i + 2 < j and s[i + 2] == 58:      # `^<n>:`, a radical
            k = i + 3
            taken = 0
            while k < j:
                if not smi_cx_next_index(s, j, &k, &idx):
                    log.append(mc_record('smarts:radical-field-truncated', (),
                                        'the radical field `%s` is malformed after %d index(es) and the '
                                        'rest of it was dropped'
                                        % (block[i:j].decode('ascii', 'replace'), taken),
                                        mc_lost()))
                    break
                taken += 1
                if idx >= p.n_atoms:
                    log.append(mc_record('smarts:radical-out-of-range', (),
                                        'the radical field names atom %d, but the string has %d atom(s); '
                                        'the mark was dropped' % (idx, p.n_atoms),
                                        mc_lost()))
                else:
                    # appended at the end of the journal, which is correct: the seal gathers an
                    # atom's tokens in journal order, so an AND_LOW here crosses the demand into
                    # every box the bracket built
                    q.atom_operator(idx + 1, SMA_AND_LOW)
                    q.atom_primitive(idx + 1, SMA_P_RADICAL, 0, False)
            if not taken:
                log.append(mc_record('smarts:radical-field-empty', (),
                                    'the radical field `%s` names no atom'
                                    % (block[i:j].decode('ascii', 'replace'))))
        else:
            log.append(mc_record('smarts:inapplicable-field', (),
                                'the extension field %s says nothing a query can test and was not applied'
                                % smi_cx_name(block, i, j)))
        i = j
    return 0


def read_smarts(text, log=None):
    """Read a chython SMARTS string into a `QueryContainer`.

    The dialect is chython's and not Daylight's: `;` is the AND between primitives, `,` the OR within
    one, and there is no recursive `$(...)`.  An implicit bond is a SINGLE bond and matches nothing
    else -- an aromatic bond is written `:` -- with one exception, which is two atoms both written
    lowercase.

    `.` separates fragments without saying whether they share a molecule.  To say so, group them:
    `(A.B)` demands one molecule component, `(A).(B)` demands two, and a bare `A.B` demands neither.
    A `(` means this only at a component position; anywhere an atom precedes it, it opens a branch as
    before.

    `log` is a list to append lines to.  Every place this reader preferred one reading of a
    contradictory string over another, and every extension field it declined to apply, puts one
    human-readable line there.  Omitting the list discards them.

    A ` |...|` CXSMARTS tail is read for its `^N:` radicals; its other fields are named in the log
    and not applied.

    Raises `IncorrectSmarts`, with a byte offset, for SYNTAX: a letter naming no primitive, an
    unbalanced bracket or parenthesis, a ring label that never closes, a value the box layout cannot
    hold.  The query is sealed before it is returned, so a term that can never match anything -- and
    an operator with no primitive beside it -- is an error here rather than a silent failure to match
    at the first use.
    """
    cdef bytes raw
    cdef object exc                      # the except-as target; `warn.undeclared` counts it
    if isinstance(text, str):
        try:
            raw = (<str> text).encode('ascii')
        except UnicodeEncodeError:
            raise IncorrectSmarts('the string contains a non-ASCII character') from None
    elif isinstance(text, bytes):
        raw = <bytes> text
    else:
        raise TypeError('read_smarts takes a str or bytes')
    raw = raw.strip()

    cdef const char *s = <const char *> raw
    cdef uint32_t n = <uint32_t> len(raw)
    cdef uint32_t k = 0
    while k < n and s[k] > 32:
        k += 1
    cdef bytes tail = raw[k:].strip()
    n = k

    cdef sma_parse_t p
    cdef object mylog = log if log is not None else []
    cdef QueryContainer q = QueryContainer()
    sma_alloc(&p, n)
    try:
        sma_tokenize(q, &p, s, n, mylog)
        if tail:
            if tail.startswith(b'|') and tail.endswith(b'|') and len(tail) > 1:
                sma_cx(q, &p, tail, mylog)
            elif tail.startswith(b'|'):
                mylog.append(mc_record('smarts:unterminated-tail', (),
                                      'the extension block after the SMARTS is not terminated and was '
                                      'ignored: %s' % tail.decode('ascii', 'replace'),
                                      mc_lost()))
            else:
                mylog.append(mc_record('smarts:trailing-text', (),
                                      'text after the SMARTS is not part of it and was ignored: %s'
                                      % tail.decode('ascii', 'replace'),
                                      mc_lost()))
    finally:
        sma_free(&p)
    try:
        q.atom_count_sealed()
    except ValueError as exc:
        # the seal names the atom or bond it could not compile, which is as located as it gets: by
        # then the byte offsets are gone and the atom number is the thing a writer can act on
        raise IncorrectSmarts('%s' % exc) from None
    return q
