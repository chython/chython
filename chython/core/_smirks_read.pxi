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
# The SMIRKS reader: `reactants>>products` becomes a ReactionTemplate.
#
# WHAT THIS LAYER IS AND IS NOT
#
# It is a READER.  It splits the arrow, runs both sides through the ONE SMARTS lexer next door, seals
# the reactant side, and works out the facts that only a two-sided string can state: which atoms pair
# by map number, which are deleted, which are created, and what the extension tail says about each
# side.  It applies nothing to a molecule -- there is no patcher here yet and no `__call__` on the
# template, so a caller cannot mistake this for a working reactor.
#
# WHY THE TWO SIDES ARE NOT SYMMETRIC
#
# The reactant side is a QUERY: it is sealed, its primitives compile to boxes, and it is matched.  The
# product side is a PATCH: it is matched against nothing, so a primitive on it is either a thing to
# BUILD (element, charge, isotope, radical, stereo) or a thing to CHECK ON THE RESULT (`r`, `D`, `h`,
# `x`, `z`, `R`).  It is therefore read into a QueryContainer used purely as a PARSE BUFFER and never
# sealed -- one lexer, two meanings for what it produces.  Classifying every product primitive into
# exactly one of those two roles, and refusing any that fits neither, is `smk_classify_products`
# below -- see its own comment for why that refusal is the thing that keeps dead template surface
# from existing at all.
#
# That asymmetry is the point of the notation: a product side that is itself a query has no way to
# tell a thing to build from a thing to test, and code written for the first case never runs.
#
# THE INDEX SPACE OF THE EXTENSION TAIL
#
# One tail for the whole string, at the end, with indices running over the reactant atoms and then the
# product atoms -- Daylight's rule for reaction CXSMILES, and the only rule under which a single tail
# can address both sides.  So in
#
#     [C;@:1][Br;D1]>>[C;@:1][O;D1;h1] |&1:2|
#
# atom 2 is the product carbon: the reactant side spent 0 and 1.  A field naming a reactant atom keeps
# the reactant-side treatment exactly as `read_smarts` gives it -- radicals applied, a stereo group
# logged as something no query can test -- while the same field naming a product atom records a
# directive on the patch.  One field may name both; each index is routed on its own.
#
# WHITESPACE
#
# `read_smarts` can say "the SMARTS ends at the first space" because a SMARTS is one token.  A SMIRKS
# is not: `A >> B` is how a human writes it.  So the tail is split off from the END -- the last
# whitespace-delimited token, taken only when it begins with `|` -- and each side is stripped after
# the arrow split.  Whitespace left INSIDE a side is refused with its position, because at that point
# the string is not a spacing preference but a structure this reader cannot guess.
#
# Byte offsets in an error message are offsets into the SIDE the message names, not into the whole
# SMIRKS.  Two sides through one lexer is what buys the shared dialect; a shared offset space is not
# something either side can see.


with cython.warn.undeclared(False):
    # bare so Python can import it, guarded so warn.undeclared stays quiet
    class IncorrectSmirks(IncorrectSmarts):
        """The string is not a SMIRKS: the reader could not decide what template it names.

        Raised for SYNTAX -- a missing or malformed arrow, the three-part reaction form, whitespace
        inside a side, a map number naming two atoms of one side -- and for a side that is not a
        readable SMARTS, in which case the message names the side and carries the lexer's own offset
        into it.  Subclasses `IncorrectSmarts`, which subclasses `IncorrectSmiles`, so a pipeline that
        catches either of those catches this too.
        """


# The properties the mapped-pair lint compares, as bits of one mask per atom.  Only these three,
# because only these three are things the product side silently RESETS: an unstated charge is zero, an
# unstated isotope is none, an unstated radical is none.  An unstated `D` or `r` resets nothing -- it
# is a check the patched product either passes or does not.
cdef enum:
    SMK_STATED_CHARGE = 1
    SMK_STATED_ISOTOPE = 2
    SMK_STATED_RADICAL = 4


cdef dict smk_journal_maps(QueryContainer q):
    """{stable id: map number} straight off the journal, without sealing anything.

    `QueryContainer.map_numbers()` seals, and the product side is never sealed -- it is a patch, and a
    patch has no boxes to compile.  So this walks the ops instead.  A later op wins, which is what the
    seal does too.
    """
    cdef dict out = {}
    cdef uint32_t i
    for i in range(q._journal_len):
        if q._journal[i].op == QOP_SET_MAP:
            if q._journal[i].value:
                out[q._journal[i].a] = <int> q._journal[i].value
            else:
                out.pop(q._journal[i].a, None)
    return out


cdef dict smk_journal_stereo_groups(QueryContainer q):
    """{stable id: (kind, group)} for every `&<n>` / `o<n>` written in a bracket, off the journal.

    The same walk as `smk_journal_maps`, for the same reason: the product side never seals, and the
    reactant side's copy of this op is refused by the seal rather than compiled.  A later op wins.
    """
    cdef dict out = {}
    cdef uint32_t i
    for i in range(q._journal_len):
        if q._journal[i].op == QOP_SET_STEREO_GROUP:
            out[q._journal[i].a] = (<int> q._journal[i].kind, <int> q._journal[i].value)
    return out


cdef tuple smk_journal_bonds(QueryContainer q):
    """Every bond of one side, as a low-first `(a, b)` pair of stable ids, straight off the journal.

    `QueryContainer` HAS NO BOND ACCESSOR AT ALL, and that is not an omission: a sealed query answers
    about boxes and DFS positions, and its bonds have been compiled into an adjacency the matcher
    walks, not into a list a caller can ask "is 3 bonded to 4" of.  The journal is therefore the only
    place either side's bonds are written down in the ids the template pairs by, which is why
    `smk_classify_products` walks the same op for the product side.  Deletion by absence needs the
    reactant list too: a bond survives the patch when its two endpoints pair and the product side
    states the bond between their partners.
    """
    cdef list out = []
    cdef uint32_t i, u, v
    for i in range(q._journal_len):
        if q._journal[i].op == QOP_ADD_BOND:
            u = q._journal[i].a
            v = q._journal[i].b
            out.append((u, v) if u < v else (v, u))
    return tuple(out)


cdef dict smk_journal_directions(QueryContainer q):
    """`{(a, b): SMI_DIR_UP/SMI_DIR_DOWN}` for every `/` and `\\`, `a` the atom written FIRST.

    The one journal walk whose key is NOT normalised low-first: a direction states a side relative to
    the atom it is written from, so which end came first is the statement rather than a spelling of it.
    `smk_directions` reads it from whichever terminal it needs, exactly as `smi_dir_from` does.
    """
    cdef dict out = {}
    cdef uint32_t i
    for i in range(q._journal_len):
        if q._journal[i].op == QOP_SET_BOND_DIRECTION:
            out[(<object> q._journal[i].a, <object> q._journal[i].b)] = <int> q._journal[i].value
    return out


cdef dict smk_journal_stated(QueryContainer q):
    """{stable id: SMK_STATED_* mask} for the three properties the lint compares.

    NEGATED primitives do not count.  `[C;!+]` says something about charge, but what it says is
    satisfied by the neutral atom an unstated product charge produces, so there is no surprise in it
    to report.  Only a positive statement can be silently undone.
    """
    cdef dict out = {}
    cdef uint32_t i
    cdef int bit, held
    cdef object sid
    for i in range(q._journal_len):
        if q._journal[i].op != QOP_ATOM_TOKEN or q._journal[i].opcode != OPC_PRIM:
            continue
        if q._journal[i].negated:
            continue
        if q._journal[i].kind == PRIM_CHARGE:
            bit = SMK_STATED_CHARGE
        elif q._journal[i].kind == PRIM_ISOTOPE:
            bit = SMK_STATED_ISOTOPE
        elif q._journal[i].kind == PRIM_RADICAL:
            bit = SMK_STATED_RADICAL
        else:
            continue
        sid = q._journal[i].a
        held = out.get(sid, 0)
        out[sid] = held | bit
    return out


cdef int smk_side(QueryContainer q, bytes src, str side, object log) except -1:
    """Run one side through the SMARTS lexer, naming the side in anything it raises."""
    cdef sma_parse_t p
    cdef const char *s = <const char *> src
    cdef uint32_t n = <uint32_t> len(src)
    cdef list sublog = []
    cdef object exc                      # the except-as target; `warn.undeclared` counts it
    cdef object line
    sma_alloc(&p, n)
    try:
        try:
            sma_tokenize(q, &p, s, n, sublog)
        except IncorrectSmarts as exc:
            # the lexer's offsets are into `src`, and the message has to say so or a reader counts
            # characters from the wrong end of the arrow
            raise IncorrectSmirks('the %s side `%s` is not a readable SMARTS: %s'
                                  % (side, src.decode('ascii', 'replace'), exc)) from None
    finally:
        sma_free(&p)
    for line in sublog:
        log.append(mc_record('smirks:side-message', (), '%s side: %s' % (side, line)))
    return 0


cdef int smk_cx(QueryContainer reactants, bytes block, uint32_t n_react, uint32_t n_prod,
                set prod_radicals, dict prod_stereo_groups, object log) except -1:
    """The one `|...|` tail, routed per index between the two sides.

    Field scanning is `_smiles_read.pxi`'s -- `smi_cx_field_end`, `smi_cx_next_index`, `smi_cx_name`
    -- for the reason `sma_cx` gives: a tail is a tail, and two implementations of "a comma ends a
    field only when a non-digit follows" would drift.
    """
    cdef const char *s = <const char *> block
    cdef uint32_t n = <uint32_t> len(block)
    cdef uint32_t i = 1
    cdef uint32_t j, k, sid
    cdef uint32_t idx = 0
    cdef uint32_t taken, applied, declined
    cdef uint32_t group = 0
    cdef uint32_t n_total = n_react + n_prod
    cdef char c
    cdef uint8_t kind = 0
    cdef bint is_radical
    if n and s[n - 1] == 124:                        # the closing `|` is not a field
        n -= 1
    while i < n:
        if s[i] == 44:
            i += 1
            continue
        j = smi_cx_field_end(s, n, i)
        if j <= i:
            break
        c = s[i]
        is_radical = False
        if c == 94 and i + 2 < j and s[i + 2] == 58:              # `^<n>:`, a radical
            is_radical = True
            k = i + 3
        elif c == 97 and i + 1 < j and s[i + 1] == 58:            # `a:`, absolute
            kind = SMI_SG_ABS
            group = 0
            k = i + 2
        elif c == 111 or c == 38:                                 # `o<n>:` OR, `&<n>:` AND
            k = i + 1
            group = smi_digits(s, j, &k, 4)
            if k == i + 1 or k >= j or s[k] != 58:
                log.append(mc_record('smirks:malformed-field', (),
                                    'the extension field %s is malformed and was not applied'
                                    % smi_cx_name(block, i, j)))
                i = j
                continue
            kind = SMI_SG_OR if c == 111 else SMI_SG_AND
            k += 1
        else:
            log.append(mc_record('smirks:inapplicable-field', (),
                                'the extension field %s says nothing this reader can apply to either side'
                                % smi_cx_name(block, i, j)))
            i = j
            continue

        taken = 0
        applied = 0
        declined = 0
        while k < j:
            if not smi_cx_next_index(s, j, &k, &idx):
                log.append(mc_record('smirks:field-truncated', (),
                                    'the extension field `%s` is malformed after %d index(es) and the rest '
                                    'of it was dropped' % (block[i:j].decode('ascii', 'replace'), taken),
                                    mc_lost()))
                break
            taken += 1
            if idx >= n_total:
                log.append(mc_record('smirks:field-out-of-range', (),
                                    'the extension field names atom %d, but the reaction has %d atom(s) '
                                    '(%d reactant, %d product); the mark was dropped'
                                    % (idx, n_total, n_react, n_prod),
                                    mc_lost()))
            elif idx < n_react:
                sid = idx + 1
                if is_radical:
                    # appended at the end of the journal, which is correct: the seal gathers an
                    # atom's tokens in journal order, so an AND_LOW here crosses the demand into
                    # every box the bracket built
                    reactants.atom_operator(sid, SMA_AND_LOW)
                    reactants.atom_primitive(sid, SMA_P_RADICAL, 0, False)
                    applied += 1
                else:
                    declined += 1
            else:
                sid = idx - n_react + 1
                if is_radical:
                    prod_radicals.add(sid)
                elif sid in prod_stereo_groups:
                    raise IncorrectSmirks(
                        'the extension field %s names product atom %d, which already states an '
                        'enhanced-stereo group in its bracket; one atom, one group'
                        % (smi_cx_name(block, i, j), sid))
                else:
                    prod_stereo_groups[sid] = (<int> kind, <int> group)
                applied += 1
        if not taken:
            log.append(mc_record('smirks:field-empty', (),
                                'the extension field `%s` names no atom'
                                % block[i:j].decode('ascii', 'replace')))
        elif declined:
            # the reactant side is a query and a stereo group is not a testable thing, which is what
            # `read_smarts` says about the same field; only the product half is a directive
            log.append(mc_record('smirks:reactant-group-skipped', (),
                                'the extension field %s names %d reactant atom(s); a query cannot test a '
                                'stereo group, so that part was not applied'
                                % (smi_cx_name(block, i, j), declined)))
        i = j
    return 0


# ------------------------------------------------------------------------------------------------
# WHAT A PRODUCT PRIMITIVE MEANS, AND THE RULE THAT KEEPS THE DEAD ONES OUT
# ------------------------------------------------------------------------------------------------
#
# Every primitive the lexer can put on a product atom or bond has to be one of two things:
#
#   BUILD  -- it says what to make.  Element, isotope (including the leading `0` that says "no mass
#             number"), charge, radical, stereo sign, bond order.
#   CHECK  -- it says nothing about what to make and everything about what the result must turn out
#             to be: `D`, `h`, `H`, `x`, `z`, `r`, `R`, `M` (the metal test), a ring bond.  These are
#             the product-side post-filter, and `r` is how a cyclization template states its ring
#             size -- there is no keyword argument for it, because the notation is the SMARTS.
#
# `smk_prim_role` places every kind the lexer has.  A kind it cannot place is REFUSED, and that
# branch is the ratchet: adding a primitive to the SMARTS lexer without deciding what it means on a
# product side makes every template using it fail to read, instead of parsing and quietly doing
# nothing.
#
# Three more things are refused, all of them for the same reason -- a patch states, it does not
# wonder:
#
#   * a BUILD primitive inside a `,` alternative.  `[C,N:1]` asks the patcher to choose an element;
#     `~` on a product bond is the same refusal, because the lexer spells it as an OR of five orders.
#   * a NEGATED build primitive.  `[C;!+:1]` names no charge to set, and since an unstated product
#     charge is zero it is satisfied by construction anyway.
#   * the same field stated twice, however it is spelled: `[C;+;+2:1]`, `[C;+&+2:1]`.
#
# `M` in the masked sense is refused too: it protects a MATCHED atom from deletion, and the product
# side matches nothing.  (`[M]` as the first primitive in a bracket is the metal test, a different
# token that lands in CHECK.)
#
# And the element is the one field with no default to fall back on -- there is no neutral element the
# way there is a neutral charge.  So a product atom states an element, or it pairs with a reactant
# atom and INHERITS one; `[A:1]` is that inheritance written out loud.  An atom that does neither
# cannot be built and is refused.
#
# `#0` IS AN ELEMENT HERE AND NOWHERE ELSE: it builds the R marker, the attachment point a template
# leaves where a fragment was cut off, and it is the only spelling of element 0 in the dialect.  A
# query holding one is refused at the seal, because an R matches nothing -- so `#0` reads on both sides
# of the arrow and survives only on the product side, which is the one side that never seals.

cdef enum:
    SMK_ROLE_UNPLACED = 0
    SMK_ROLE_BUILD = 1
    SMK_ROLE_CHECK = 2

# The FIELD a build primitive writes, so that two spellings of one field collide: `13C` and the
# leading `0` are both the isotope, `C` and `A` are both the element, and `=` and `:` are both the
# bond's order.
cdef enum:
    SMK_FIELD_ELEMENT = 1
    SMK_FIELD_ISOTOPE = 2
    SMK_FIELD_CHARGE = 3
    SMK_FIELD_RADICAL = 4
    SMK_FIELD_STEREO = 5
    SMK_FIELD_ORDER = 6


cdef inline int smk_prim_role(uint32_t kind) noexcept nogil:
    if (kind == PRIM_ELEMENT or kind == PRIM_ANY or kind == PRIM_R_MARKER or kind == PRIM_ISOTOPE
            or kind == PRIM_NO_ISOTOPE or kind == PRIM_CHARGE or kind == PRIM_RADICAL
            or kind == PRIM_STEREO or kind == PRIM_STEREO_KEEP
            or kind == PRIM_STEREO_INVERT or kind == BPRIM_ORDER
            or kind == BPRIM_AROMATIC):
        return SMK_ROLE_BUILD
    if (kind == PRIM_METAL or kind == PRIM_DEGREE or kind == PRIM_IMPLICIT_H
            or kind == PRIM_TOTAL_H or kind == PRIM_HETEROATOMS or kind == PRIM_HYBRIDIZATION
            or kind == PRIM_RING_SIZE or kind == PRIM_RING_COUNT or kind == BPRIM_RING):
        return SMK_ROLE_CHECK
    return SMK_ROLE_UNPLACED


cdef inline int smk_prim_field(uint32_t kind) noexcept nogil:
    if kind == PRIM_ELEMENT or kind == PRIM_ANY or kind == PRIM_R_MARKER:
        return SMK_FIELD_ELEMENT
    if kind == PRIM_ISOTOPE or kind == PRIM_NO_ISOTOPE:
        return SMK_FIELD_ISOTOPE
    if kind == PRIM_CHARGE:
        return SMK_FIELD_CHARGE
    if kind == PRIM_RADICAL:
        return SMK_FIELD_RADICAL
    if kind == PRIM_STEREO or kind == PRIM_STEREO_KEEP or kind == PRIM_STEREO_INVERT:
        return SMK_FIELD_STEREO
    return SMK_FIELD_ORDER


# How a primitive is NAMED in a refusal.  A template author reads the message and has to find the
# token in their own string, so each entry says the token and not the internal kind.
cdef dict SMK_PRIM_NAMES = {
    <int> PRIM_ELEMENT: 'an element',
    <int> PRIM_ANY: '`A`, which inherits its element',
    <int> PRIM_R_MARKER: '`#0`, the R marker',
    <int> PRIM_METAL: '`M`, the metal test',
    <int> PRIM_ISOTOPE: 'an isotope',
    <int> PRIM_NO_ISOTOPE: 'a leading `0`, the no-isotope statement',
    <int> PRIM_ANY_CHARGE: '`*`, the any-charge wildcard',
    <int> PRIM_CHARGE: 'a charge',
    <int> PRIM_RADICAL: 'a radical',
    <int> PRIM_STEREO: 'a stereo sign',
    <int> PRIM_STEREO_KEEP: '`@=`, the configuration the reactant had',
    <int> PRIM_STEREO_INVERT: '`@~`, the other configuration',
    <int> PRIM_DEGREE: '`D`',
    <int> PRIM_IMPLICIT_H: '`h`',
    <int> PRIM_TOTAL_H: '`H`',
    <int> PRIM_HETEROATOMS: '`x`',
    <int> PRIM_HYBRIDIZATION: '`z`',
    <int> PRIM_RING_SIZE: '`r`',
    <int> PRIM_RING_COUNT: '`R`',
    <int> BPRIM_ORDER: 'a bond order',
    <int> BPRIM_AROMATIC: 'an aromatic bond',
    <int> BPRIM_RING: 'a ring bond',
}

cdef dict SMK_FIELD_NAMES = {
    <int> SMK_FIELD_ELEMENT: 'element',
    <int> SMK_FIELD_ISOTOPE: 'isotope',
    <int> SMK_FIELD_CHARGE: 'charge',
    <int> SMK_FIELD_RADICAL: 'radical',
    <int> SMK_FIELD_STEREO: 'stereo sign',
    <int> SMK_FIELD_ORDER: 'bond order',
}


cdef inline str smk_prim_name(uint32_t kind):
    cdef object name = SMK_PRIM_NAMES.get(<int> kind)
    if name is None:
        return 'primitive kind %d' % kind
    return <str> name


# How a CHECK primitive is spelled back at a template author -- as the token they wrote, value and
# all, because a post-filter rejection is only useful if they can find the primitive that did it.
# `SMK_PRIM_NAMES` above cannot serve: it names a kind without its value, and `r5` failing where
# `r6` would have held is the whole of what a chemist needs to see.
cdef dict SMK_CHECK_LETTERS = {
    <int> PRIM_DEGREE: 'D',
    <int> PRIM_IMPLICIT_H: 'h',
    <int> PRIM_TOTAL_H: 'H',
    <int> PRIM_HETEROATOMS: 'x',
    <int> PRIM_HYBRIDIZATION: 'z',
    <int> PRIM_RING_SIZE: 'r',
    <int> PRIM_RING_COUNT: 'R',
}


cdef int smk_refuse_unplaced(str what, uint32_t kind) except -1:
    """The refusal for a primitive with no product-side role, tailored where the fix is obvious.

    One function rather than two raise sites, because `smk_place` reaches this from both the
    single-alternative and the disjunction branch and the message is the same question either way.
    """
    if kind == PRIM_ANY_CHARGE:
        raise IncorrectSmirks('%s carries `*`, the any-charge wildcard.  On a reactant side that '
                              'withdraws the neutral default and so matches any charge; a patch has '
                              'no default to withdraw and no charge to pick, so `*` would build '
                              'nothing and test nothing.  State the charge outright, or leave it '
                              'unstated for zero' % what)
    raise IncorrectSmirks('%s carries %s, which this reader can neither build nor test on a '
                          'product; every product primitive has to be one or the other'
                          % (what, smk_prim_name(kind)))


cdef str smk_check_spelling(tuple clause):
    """One check clause as its own notation: `r5`, `!R0`, `D2;z1` re-joined, alternatives by `,`."""
    cdef list alts = [], names
    cdef object alt, prim
    cdef uint32_t kind
    cdef object letter
    for alt in clause:
        names = []
        for prim in <tuple> alt:
            kind = <uint32_t> <int> (<tuple> prim)[0]
            letter = SMK_CHECK_LETTERS.get(<int> kind)
            if letter is None:
                # the two valueless ones: `M` the metal test, and `@` the ring bond
                names.append(('!' if <int> (<tuple> prim)[2] else '')
                             + ('M' if kind == PRIM_METAL else '@'))
            else:
                names.append('%s%s%d' % ('!' if <int> (<tuple> prim)[2] else '',
                                         <str> letter, <int> (<tuple> prim)[1]))
        alts.append('&'.join(names))
    return ','.join(alts)


cdef inline str smk_atom_name(uint32_t sid, dict prod_maps):
    """Name a product atom the way the string names it: by position, and by map number when it has
    one, because a template author thinks in map numbers and counts atoms only when made to."""
    cdef object number = prod_maps.get(<object> sid)
    if number is None:
        return 'product atom %d' % sid
    return 'product atom %d (map number %d)' % (sid, <int> number)


cdef list smk_clauses(list toks):
    """One key's token run, split by precedence into `;`-clauses of `,`-alternatives of primitives.

    The dialect has three levels and no parentheses, so any expression in it is a conjunction (over
    `;`) of disjunctions (over `,`) of conjunctions (over `&`, or over juxtaposition).  That shape is
    the whole reason a product primitive can be classified at all: a primitive in a clause with no
    `,` in it is unconditionally true of the product, and one inside a `,` is true of only some
    alternative -- so the first can build and the second can only ever check.
    """
    cdef list clauses = []
    cdef list clause = [[]]
    cdef list alt
    cdef tuple tok
    cdef int opc
    for tok in toks:
        opc = <int> tok[0]
        if opc == OPC_PRIM:
            alt = <list> clause[len(clause) - 1]
            alt.append((tok[1], tok[2], tok[3]))
        elif opc == OPC_OR:
            clause.append([])
        elif opc == OPC_AND_LOW:
            clauses.append(clause)
            clause = [[]]
        # OPC_AND_HIGH joins two primitives inside one alternative, which is where they already are
    clauses.append(clause)
    return clauses


cdef tuple smk_freeze(list clause):
    """A clause as nested tuples, so what lands on the template cannot be edited from outside."""
    cdef list out = []
    cdef object alt
    for alt in clause:
        out.append(tuple(<list> alt))
    return tuple(out)


cdef int smk_place(list toks, str what, set fields, list build, list checks) except -1:
    """Put every primitive of one product atom or bond into exactly one of the two roles.

    `build` collects `(kind, value)` in the order written; `checks` collects whole clauses, each a
    tuple of alternatives, each alternative a tuple of `(kind, value, negated)`.  A clause with one
    alternative keeps only its check primitives -- its build primitives have moved to `build`, where
    they are unconditional, so what is left to test is the rest of the same conjunction.
    """
    cdef list clauses = smk_clauses(toks)
    cdef list clause, alt, keep
    cdef tuple prim
    cdef int role, field
    cdef uint32_t kind
    cdef object c, a
    for c in clauses:
        clause = <list> c
        if len(clause) > 1:
            # a `,` disjunction: no primitive in it is certain, so none of them may build
            for a in clause:
                for prim in <list> a:
                    kind = <uint32_t> <int> prim[0]
                    role = smk_prim_role(kind)
                    if role == SMK_ROLE_UNPLACED:
                        smk_refuse_unplaced(what, kind)
                    if role == SMK_ROLE_BUILD:
                        raise IncorrectSmirks('%s offers %s as one of several alternatives; a patch '
                                              'builds one thing and cannot choose between them, so '
                                              'the %s has to be stated outright'
                                              % (what, smk_prim_name(kind),
                                                 <str> SMK_FIELD_NAMES[smk_prim_field(kind)]))
            checks.append(smk_freeze(clause))
            continue
        alt = <list> clause[0]
        keep = []
        for prim in alt:
            kind = <uint32_t> <int> prim[0]
            role = smk_prim_role(kind)
            if role == SMK_ROLE_UNPLACED:
                smk_refuse_unplaced(what, kind)
            if role == SMK_ROLE_CHECK:
                keep.append(prim)
                continue
            if <int> prim[2]:
                raise IncorrectSmirks('%s negates %s; a patch states what to build, and a negation '
                                      'names nothing to build' % (what, smk_prim_name(kind)))
            field = smk_prim_field(kind)
            if field in fields:
                raise IncorrectSmirks('%s states the %s twice'
                                      % (what, <str> SMK_FIELD_NAMES[field]))
            fields.add(field)
            build.append((prim[0], prim[1]))
        if keep:
            checks.append((tuple(keep),))
    return 0


cdef int smk_classify_products(QueryContainer q, dict prod_maps, set paired,
                               dict atom_build, dict atom_check, list bonds,
                               dict bond_build, dict bond_check, set inherited) except -1:
    """Classify the whole product side, filling the six outputs, and refuse what will not classify.

    Walks the journal rather than a sealed arena, for the reason the whole layer does: the product
    side is a patch and never seals.  Bond keys are normalised low-first, so a caller does not have
    to know which way round the string wrote the bond.
    """
    cdef uint32_t i, u, v, kind
    cdef uint32_t op
    cdef dict runs_a = {}
    cdef dict runs_b = {}
    cdef object key
    cdef tuple bkey
    cdef list run
    cdef set fields
    cdef list build, checks
    cdef tuple prim
    cdef bint stated, inherits

    for i in range(q._journal_len):
        op = q._journal[i].op
        if op == QOP_SET_MASKED:
            raise IncorrectSmirks('%s carries `M`, which protects a MATCHED atom from deletion; the '
                                  'product side matches nothing, so there is nothing to protect'
                                  % smk_atom_name(q._journal[i].a, prod_maps))
        elif op == QOP_ADD_BOND:
            u = q._journal[i].a
            v = q._journal[i].b
            bonds.append((u, v) if u < v else (v, u))
        elif op == QOP_ATOM_TOKEN:
            key = <object> q._journal[i].a
            run = <list> runs_a.get(key)
            if run is None:
                run = []
                runs_a[key] = run
            run.append((<int> q._journal[i].opcode, <int> q._journal[i].kind,
                        <int> q._journal[i].value, <int> q._journal[i].negated))
        elif op == QOP_BOND_TOKEN:
            u = q._journal[i].a
            v = q._journal[i].b
            key = (u, v) if u < v else (v, u)
            run = <list> runs_b.get(key)
            if run is None:
                run = []
                runs_b[key] = run
            run.append((<int> q._journal[i].opcode, <int> q._journal[i].kind,
                        <int> q._journal[i].value, <int> q._journal[i].negated))

    for i in range(1, q.atom_count + 1):
        key = <object> i
        run = <list> runs_a.get(key)
        fields = set()
        build = []
        checks = []
        smk_place(run if run is not None else [], smk_atom_name(i, prod_maps),
                  fields, build, checks)
        # the element has no default to fall back on, so it is stated, inherited, or refused
        stated = False
        inherits = False
        for prim in build:
            kind = <uint32_t> <int> prim[0]
            if kind == PRIM_ELEMENT or kind == PRIM_R_MARKER:
                stated = True
            elif kind == PRIM_ANY:
                inherits = True
        if not stated:
            if i not in paired:
                if inherits:
                    raise IncorrectSmirks('%s writes `A`, which takes its element from the reactant '
                                          'atom it pairs with, and it pairs with none'
                                          % smk_atom_name(i, prod_maps))
                raise IncorrectSmirks('%s states no element and pairs with no reactant atom, so '
                                      'there is nothing to build it from: name an element, or map '
                                      'it onto one' % smk_atom_name(i, prod_maps))
            inherited.add(i)
        atom_build[key] = tuple(build)
        atom_check[key] = tuple(checks)

    for bkey in bonds:
        run = <list> runs_b.get(bkey)
        fields = set()
        build = []
        checks = []
        # An untokenised bond journals nothing at all and is a single -- the seal supplies that for a
        # query, and the patcher supplies it here, from an empty build list.
        smk_place(run if run is not None else [],
                  'the product bond between atoms %d and %d' % (bkey[0], bkey[1]),
                  fields, build, checks)
        bond_build[bkey] = tuple(build)
        bond_check[bkey] = tuple(checks)
    return 0


cdef int smk_bond_order(tuple build) noexcept:
    """The order one product bond builds.  An untokenised bond builds a single, exactly as the seal
    supplies one for a query -- `smk_place` hands this an empty list for it and nothing else has to
    know.

    Lives on the reading side rather than in the patcher because the frame check below has to ask the
    same question: whether a stereo sign's anchor could be a tetrahedral centre at all is decided from
    the orders its bonds BUILD, and that has to be answered while the string is still in hand.
    """
    cdef int order = 1
    cdef tuple prim
    cdef int kind
    for prim in build:
        kind = <int> prim[0]
        if kind == BPRIM_ORDER:
            order = <int> prim[1]
        elif kind == BPRIM_AROMATIC:
            order = 4
    return order


# ------------------------------------------------------------------------------------------------
# WHAT A PRODUCT SIDE CAN SAY ABOUT A CONFIGURATION
# ------------------------------------------------------------------------------------------------
#
# A patch DROPS the configuration of every stereo unit it wrote a part of (`smk_rc_stereo`): a bond
# the template changed at one of a unit's OWN atoms is a bond the reaction made or broke, and nothing
# carries a configuration across that unless a template says so.  Away from the reaction centre
# nothing is dropped and none of this is needed.
#
# FOUR STATEMENTS, every one of them on the product side:
#
#   `@=`              the configuration the reactant had, unchanged
#   `@~`              the other one -- every unit chython models has exactly two states
#   `&<n>` / `o<n>`   a mixture: the unit comes out configured AND grouped
#   a sign, in one group with another signed product atom
#                     a DRAWN configuration, relative to that other one, in the product's own frame
#
# `@=` AND `@~` ARE RELATIVE AND TAKE NO FRAME.  What the arena re-based is what is kept or flipped,
# and the arena took its frame from the MOLECULE, which is the only place a frame is ever fully known.
# So the SN2 inversion a template exists to express is
#
#     [C;z1:1][Br;D1] >> [C@~:1][O;D1:2]
#
# and the retention is the same string with `@=`.  Both refuse a created atom, which had no
# configuration to be relative to, and both are one field, so `smk_place` refuses a string stating
# both.  ONE TOKEN FOR EVERY KIND, and it has to be: a cis/trans, allene or atropisomer configuration
# has no sign spelling.  `@=` on either end of a double bond keeps that bond's geometry and `@~` flips
# it; on a tetrahedral centre they keep and flip the parity.
#
# A REACTANT-SIDE SIGN IS SELECTIVITY AND NOTHING ELSE -- it narrows the match to a configured atom,
# and only that, because a sign naming fewer than three directions is unenforceable as a value and
# `_isomorphism.pxi` widens it to "configured, either sign".  No product statement reads it.
#
# THERE IS NO ABSOLUTE SIGN, deliberately, and its absence is chemistry rather than economy.  A
# configuration cannot appear where no chiral influence acted, so a template that sets one out of
# nothing is describing a reaction that does not happen.  Worse than useless: in a
# substrate-controlled diastereoselection -- Felkin-Anh addition, directed epoxidation of an allylic
# alcohol, 1,2-trans glycosylation -- an absolute product sign is actively WRONG, because it would
# turn the enantiomeric substrate into the same absolute product.  What those reactions state is a
# RELATIVE configuration, and a correlated group is its notation.
#
# The signs of a correlated group fix its members' configurations relative to EACH OTHER and the group
# says the set is a mixture; where the substrate carries one member configured, the patcher mirrors the
# drawn set to agree with it and writes no group at all.  One template, and the answer follows the
# substrate -- which is the whole reason to have this over a pair of enantiomer-specific rules.
#
# A GEOMETRY IS NOT A CONFIGURATION IN THIS RESPECT.  `/` and `\` on the product side state an
# absolute cis/trans geometry, which needs no chiral influence: a Wittig makes its alkene E or Z by
# mechanism.  They are bond tokens and are read in `smk_directions`, not here.
#
# THE REFUSALS, all at the string:
#
#   * a sign with no correlated partner: it could only be an absolute setting;
#   * a correlated sign whose anchor names fewer than three or more than four directions.  A drawn
#     configuration is relative to an order of the centre's directions, and nothing else describes a
#     tetrahedron;
#   * any sign whose anchor carries a bond the product BUILDS as anything but single.  `@` on an atom
#     denotes a tetrahedral unit and only ever that; an allene centre, a cis/trans anchor and an
#     aromatic atom are all reachable this way and none of them has an atom-sign spelling;
#   * `@=` or `@~` on an atom the reaction creates.

cdef dict smk_correlated_members(dict signed, dict groups):
    """The product atoms whose sign is a DRAWN configuration: `{atom: (kind, number)}`.

    A group correlates when two or more of its members are signed.  One signed member is a
    configuration relative to nothing, so it is refused.  `abs` carries no number and never
    correlates: it is the absence of a mixture, so there is no set for a configuration to be one
    member of.
    """
    cdef dict counts = {}
    cdef dict out = {}
    cdef object key
    cdef tuple spec
    for key in signed:
        spec = <tuple> groups.get(key)
        if spec is None or <int> spec[0] == SMI_SG_ABS:
            continue
        counts[spec] = <int> counts.get(spec, 0) + 1
    for key in signed:
        spec = <tuple> groups.get(key)
        if spec is not None and <int> counts.get(spec, 0) > 1:
            out[key] = spec
    return out


cdef tuple smk_stereo_reading(dict prod_maps, dict pairs, dict groups, dict atom_build,
                              dict bond_build, list bonds):
    """`({atom: (sign, frame)}, {atom}, {atom})` -- correlated, `@=` and `@~`.

    `frame` is the atom's product-side directions in ascending id order -- ruling F26's order, the one
    the reactant side's `qstereo_t.refs` also uses -- `None`-padded to four where the fourth is an
    implicit hydrogen, which is the shape `MoleculeContainer.translate_stereo` takes.  A sign is 1 for
    `@` and 2 for `@@`, which is what `PRIM_STEREO` validates.

    `@=` and `@~` name no sign and so need neither a value nor a frame: two sets of atoms whose
    configuration the patch is told to carry through, or flip at, its own reaction centre.
    """
    cdef dict correlated = {}
    cdef set keep = set()
    cdef set invert = set()
    cdef dict signed = {}
    cdef dict neighbours = {}
    cdef dict prod_to_react = {}
    cdef tuple bkey, prim, pair
    cdef object key, number, sid, other
    cdef list nb
    cdef int order, sign, kind
    for number in pairs:
        pair = <tuple> pairs[number]
        prod_to_react[pair[1]] = pair[0]
    for bkey in bonds:
        for sid, other in ((bkey[0], bkey[1]), (bkey[1], bkey[0])):
            nb = <list> neighbours.get(sid)
            if nb is None:
                nb = []
                neighbours[sid] = nb
            nb.append(other)

    # every signed atom first, with the one refusal that is about the anchor's BONDS rather than about
    # which reading applies -- a sign on a non-single anchor means nothing under either
    for key in atom_build:
        sign = 0
        for prim in <tuple> atom_build[key]:
            kind = <int> prim[0]
            if kind == PRIM_STEREO:
                sign = <int> prim[1]
            elif kind == PRIM_STEREO_KEEP or kind == PRIM_STEREO_INVERT:
                # `@=`, `@~` and a sign are one field, and `smk_place` refuses the second statement of
                # a field, so these branches and the one above are exclusive.
                if key not in prod_to_react:
                    raise IncorrectSmirks(
                        '%s carries `%s` and pairs with no reactant atom, so there is no configuration '
                        'for it to be relative to.  A centre the reaction creates gets `&<n>` to '
                        'racemise it, or a sign in one group with another signed product atom to state '
                        'a configuration relative to that one'
                        % (smk_atom_name(<uint32_t> <int> key, prod_maps),
                           '@=' if kind == PRIM_STEREO_KEEP else '@~'))
                if kind == PRIM_STEREO_KEEP:
                    keep.add(key)
                else:
                    invert.add(key)
        if not sign:
            continue
        for bkey in bonds:
            if key != bkey[0] and key != bkey[1]:
                continue
            order = smk_bond_order(<tuple> bond_build[bkey])
            if order != 1:
                raise IncorrectSmirks(
                    '%s carries a stereo sign and a bond this patch builds as order %d; a sign on an '
                    'atom states a TETRAHEDRAL configuration and nothing else has that spelling. '
                    'Carrying a cis/trans, allene or atropisomer configuration through a patch is '
                    '`@=`, inverting one is `@~`, and stating a geometry outright is `/` and `\\`'
                    % (smk_atom_name(<uint32_t> <int> key, prod_maps), order))
        signed[key] = sign

    cdef dict members = smk_correlated_members(signed, groups)
    for key in signed:
        sign = <int> signed[key]
        if key not in members:
            raise IncorrectSmirks(
                '%s carries a stereo sign and shares no enhanced-stereo group with another signed '
                'product atom, so the only thing it could state is an absolute configuration -- and a '
                'configuration cannot appear where nothing chiral acted.  To carry the reactant\'s '
                'through write `@=`, to invert it `@~`, to racemise a centre `&<n>`, and to state a '
                'configuration relative to another centre put both in one group, signed'
                % smk_atom_name(<uint32_t> <int> key, prod_maps))
        nb = <list> neighbours.get(key)
        if nb is None or len(nb) < 3 or len(nb) > 4:
            raise IncorrectSmirks(
                '%s states one member of a correlated stereo group and names %d direction(s); a '
                'drawn configuration is relative to the order of its centre\'s directions, so it '
                'has to name three of them (the fourth being an implicit hydrogen) or all four'
                % (smk_atom_name(<uint32_t> <int> key, prod_maps), 0 if nb is None else len(nb)))
        nb = sorted(nb)
        correlated[key] = (sign, tuple(nb) if len(nb) == 4 else tuple(nb) + (None,))
    return correlated, keep, invert


# ------------------------------------------------------------------------------------------------
# A DRAWN GEOMETRY, WHICH IS THE ONE ABSOLUTE CONFIGURATION A PRODUCT SIDE MAY STATE
# ------------------------------------------------------------------------------------------------
#
# `/` and `\` are read here and nowhere else on this side.  What the header above `smk_stereo_reading`
# says about an absolute configuration having no spelling stops at cis/trans: a geometry needs no chiral
# influence to appear, so a stabilised Wittig stating E outright is describing what its mechanism does.
# The argument that makes an absolute SIGN wrong never arises here -- an alkene's two faces are not
# enantiomeric, so the enantiomeric substrate does not give the enantiomeric product.
#
# THE STATEMENT IS READ FROM THE CHAIN'S TERMINALS, not from the directed bonds.  One single bond
# between two chains states a side for both of them -- `C/C=C/C=C/C` is one `/` doing two jobs -- so a
# walk over the directions would see the second geometry as unstated.  Same reason `smi_stereo` drives
# its own pass off the perceived units.
#
# THE REFUSALS, all at the string:
#
#   * a direction on one terminal only: a geometry is a statement about both ends;
#   * both substituents of one terminal on the same side, which no geometry has;
#   * a direction on a chain whose terminal also carries `@=` or `@~`: two answers to one question;
#   * a direction that reaches no chain at all -- dead surface, refused as N4 refuses its own.
#
# The unit's KIND is not checked here.  A chain of an even number of double bonds is a cumulene and one
# of an odd number is an allene, and which the patched molecule holds is a fact about the molecule; the
# patcher asks `unit_of` and logs a skip, exactly as it does for a correlated sign.

cdef dict smk_directions(dict prod_maps, dict raw, list bonds, dict bond_build,
                         set keep, set invert):
    """`{(t1, t2): (sub1, sub2, trans)}` -- the geometries the product side draws.

    `t1 < t2` are the two terminals of one chain of double bonds, `sub1` and `sub2` the substituent of
    each that a `/` or `\\` named, and `trans` whether the two stand on opposite sides.
    """
    cdef dict out = {}
    cdef dict chain = {}                 # sid -> the atoms it is double-bonded to
    cdef dict nbrs = {}                  # sid -> every atom it is bonded to
    cdef set consumed = set()            # keys of `raw` this pass spent
    cdef set seen = set()                # terminals already answered, so a chain is read once
    cdef tuple bkey, pair
    cdef list arms, nb, marks
    cdef object sid, other, far, prev, t, rkey, sub, d, sub1, sub2
    cdef int order, hops, d1, d2
    cdef bint reached

    for bkey in bonds:
        order = smk_bond_order(<tuple> bond_build[bkey])
        for pair in ((bkey[0], bkey[1]), (bkey[1], bkey[0])):
            nb = <list> nbrs.get(pair[0])
            if nb is None:
                nb = []
                nbrs[pair[0]] = nb
            nb.append(pair[1])
            if order != 2:
                continue
            nb = <list> chain.get(pair[0])
            if nb is None:
                nb = []
                chain[pair[0]] = nb
            nb.append(pair[1])

    for sid in sorted(chain):
        if sid in seen or len(<list> chain[sid]) != 1:
            continue
        # walk to the far terminal.  An all-double ring has no atom with one arm and is never entered;
        # an atom with three is not a chain end this can order, and the walk gives up on it -- its
        # directions then fall to the "names no chain" refusal below.
        t = sid
        prev = None
        reached = False
        hops = 0
        while hops <= len(bonds):
            arms = []
            for other in <list> chain[t]:
                if other != prev:
                    arms.append(other)
            if not arms:
                reached = True
                break
            if len(arms) > 1:
                break
            prev = t
            t = <object> arms[0]
            hops += 1
        if not reached:
            continue
        far = t
        seen.add(sid)
        seen.add(far)

        marks = []
        for pair in ((sid, far), (far, sid)):
            t = pair[0]
            sub = None
            d = None
            for other in <list> nbrs[t]:
                if other in <list> chain[t]:
                    continue
                if (t, other) in raw:
                    rkey = (t, other)
                    d1 = <int> raw[rkey]
                elif (other, t) in raw:
                    rkey = (other, t)
                    d1 = 3 - <int> raw[rkey]     # the same statement read from this end, upside down
                else:
                    continue
                consumed.add(rkey)
                if sub is None:
                    sub = other
                    d = d1
                elif d1 == <int> d:
                    raise IncorrectSmirks(
                        '%s puts both of its substituents on the same side of the double bond it '
                        'terminates, and no geometry does that' % smk_atom_name(<uint32_t> <int> t,
                                                                               prod_maps))
            marks.append((sub, d))
        sub1 = (<tuple> marks[0])[0]
        sub2 = (<tuple> marks[1])[0]
        if sub1 is None and sub2 is None:
            continue
        if sub1 is None or sub2 is None:
            raise IncorrectSmirks(
                'the product double bond between %s and %s carries a direction on one end only; a '
                'geometry is a statement about both, so a `/` or `\\` is needed on a substituent of '
                'each' % (smk_atom_name(<uint32_t> <int> sid, prod_maps),
                          smk_atom_name(<uint32_t> <int> far, prod_maps)))
        if sid in keep or sid in invert or far in keep or far in invert:
            raise IncorrectSmirks(
                'the product double bond between %s and %s carries both a drawn geometry and `@=` or '
                '`@~`; one states the geometry outright and the other takes the reactant\'s, so they '
                'are two answers to one question' % (smk_atom_name(<uint32_t> <int> sid, prod_maps),
                                                     smk_atom_name(<uint32_t> <int> far, prod_maps)))
        d1 = <int> (<tuple> marks[0])[1]
        d2 = <int> (<tuple> marks[1])[1]
        # opposite directions, each read from its own terminal, means opposite sides
        if sid < far:
            out[(sid, far)] = (sub1, sub2, d1 != d2)
        else:
            out[(far, sid)] = (sub2, sub1, d1 != d2)

    if len(consumed) < len(raw):
        raise IncorrectSmirks(
            '%d product-side `/` or `\\` name no chain of double bonds, and a direction states nothing '
            'on its own: it says which side of a geometry a substituent is on, so there has to be a '
            'geometry for it to be part of' % (len(raw) - len(consumed)))
    return out


cdef class ReactionTemplate:
    """A SMIRKS template: a sealed query for the reactant side and a patch spec for the product side.

    Built only by `read_smirks`.  Every attribute is read-only, and the two number spaces the atoms
    live in are kept apart on purpose (see `_smirks_read.pxi`'s header and the design's N6): a STABLE
    ID identifies an atom of one side of this template, a MAP NUMBER pairs one side's atom with the
    other's.  Neither is the atom-atom mapping of a reaction the template produces.

    The reaction a template yields is mapped: contiguous from 1, one number per reactant-product pair,
    and 0 on an atom that exists on one side only.  The inputs' own map numbers are not carried -- two
    inputs each numbered from 1 cannot be made 1-1 by preserving them.
    """
    cdef readonly str smirks
    # what a log record calls this template.  A string read by hand has no name but its string; a row
    # of a table has `'reactions:13'`, and the row is what its author can go and edit.
    cdef readonly str rule_id
    cdef readonly QueryContainer reactants
    # the parse buffer, never sealed: see the header on why the product side is not a query
    cdef QueryContainer _products
    cdef readonly dict reactant_map_numbers
    cdef readonly dict product_map_numbers
    cdef readonly dict mapped_pairs
    cdef readonly frozenset deleted_atoms
    cdef readonly frozenset created_atoms
    # the reactant side's bonds, low-first, in the ids `reactant_map_numbers` keys by.  The product
    # side's are `product_bonds`; the two lists together are what deletion-by-absence compares.
    cdef readonly tuple reactant_bonds
    cdef readonly dict product_stereo_groups
    cdef readonly frozenset product_radicals
    cdef readonly uint32_t product_atom_count
    # Bonds get a tuple where atoms get a count, and that is not an inconsistency: product atom ids
    # run 1..product_atom_count, so the count IS the list, while a bond key is a pair no count
    # implies.  Keys are normalised low-first.
    cdef readonly tuple product_bonds
    # The classification: what each product atom and bond BUILDS and what it CHECKS.  `*_build` maps
    # a key to `(kind, value)` pairs in the order written, `*_check` to whole clauses -- a tuple of
    # alternatives, each a tuple of `(kind, value, negated)`, satisfied when any alternative is.
    # `kind` is the core's own PRIM_* / BPRIM_* constant.  Every product atom and bond has an entry
    # in both, empty when it states nothing of that role.
    cdef readonly dict product_atom_build
    cdef readonly dict product_atom_check
    cdef readonly dict product_bond_build
    cdef readonly dict product_bond_check
    # product atoms whose element comes from the reactant atom they pair with, either because they
    # state no element or because they state `A`
    cdef readonly frozenset product_inherited_elements
    # `{product atom: (sign, frame)}` -- the members of a correlated group, where a sign IS a drawn
    # configuration in the product side's own frame: `frame` is the atom's directions in ascending id
    # order, `None`-padded to four where the fourth is an implicit hydrogen.  The drawn set is a
    # RELATIVE configuration -- the patcher mirrors it as a whole where the substrate already decided
    # one member.  A sign with no correlated partner is refused; see the header above
    # `smk_stereo_reading` for why an absolute configuration has no spelling at all.
    cdef readonly dict product_stereo_correlated
    # The product atoms spelled `@=` and `@~`: the reactant's configuration comes through the patch
    # unchanged, or comes out as the unit's other state.  Whatever it was and whatever kind of unit
    # holds it, since every kind chython models is two-state.  Disjoint from each other and from
    # `product_stereo_correlated`, and the only way a configuration survives a reaction centre.
    cdef readonly frozenset product_stereo_keep
    cdef readonly frozenset product_stereo_invert
    # `{(t1, t2): (sub1, sub2, trans)}` -- the cis/trans geometries the product side DRAWS with `/` and
    # `\`, keyed by the two terminals of one chain of double bonds, `t1 < t2`.  Absolute, unlike every
    # other product statement about a configuration: a geometry needs no chiral influence, so a template
    # may name E or Z outright.  See the header above `smk_directions`.
    cdef readonly dict product_stereo_geometry

    @cython.warn.unused_arg(False)
    def __init__(self, *args, **kwargs):
        raise TypeError('a ReactionTemplate is read from a string: call read_smirks')

    def __repr__(self):
        return 'read_smirks(%r)' % self.smirks

    def __call__(self, *molecules, bint automorphism_filter=True, log=None, bint report=False):
        """Apply this template to one molecule or to several, yielding one `ReactionContainer` per
        distinct outcome.

        `template(a, b)`, and a collection is unpacked at the call: `template(*rxn.reactants)`.  The
        molecules are unioned into one working container, so an INTRAMOLECULAR template needs no
        special call: `A.B>>C` written with `.` matches two fragments of one input as readily as one
        fragment of each of two.

        Only the inputs the match TOUCHED appear as the reaction's reactants, and the products are the
        components those inputs became -- an input the template did not reach is on neither side.  A
        counter-ion sitting in a touched input comes out as its own product molecule; nothing drops a
        component for not being in the reaction centre.

        Hydrogen counts are recomputed for the atoms the patch actually WROTE and for the neighbours
        of what it deleted -- never for an atom that merely sat inside the match, whose stored count
        is still the count its input stated.  Where the valence collection has no answer, including
        anywhere an aromatic bond reaches the reaction centre, the count is stored as `H_UNKNOWN` and
        not as a guessed zero: `kekule()` and then `chython.chemistry.calc_implicit` are the repair,
        and they are the caller's to run.

        EVERY OUTCOME CARRIES WHAT ITS OWN PATCH REPORTED, on `rxn.log` under stage `react`.  Two things
        put a line there and neither is visible any other way: a configured parity the arena could not
        carry across the change (the design's N7 -- the group membership is cleared in the same
        operation, N8), and a candidate whose patch raised, which is logged and skipped rather than
        allowed to end the enumeration (N11).  The second has no outcome to be carried by, so it reaches
        the optional `log` list only -- which is otherwise a second view of the same records, for a
        caller enumerating a corpus who wants one sequence for the run.

        Its atoms carry the imposed mapping, not the inputs': contiguous from 1 over the atoms present
        on both sides, 0 for a leaving or an incoming one.  A template's own `:N` numbers pair its two
        sides and are a different space entirely (N6).

        With `report=True` each yield is `(reaction, {product-side map number: atom id})` instead of the
        reaction alone.  The ids are the yielded products' own, and they survive `copy()` and `split()`,
        so a caller that must edit "the atom the product side called :1" can reach it -- by map number,
        the template's `:N` space, which is still not the reaction's imposed mapping.
        """
        return smk_apply(self, molecules, automorphism_filter,
                         log if log is not None else [], report)


def read_smirks(text, log=None, *, rule_id=None):
    """Read a SMIRKS reaction template into a `ReactionTemplate`.

    `reactants>>products`, where both sides are chython SMARTS -- the dialect `read_smarts`
    documents, `;` for AND and `,` for OR, with no recursive `$(...)`.  Whitespace may surround the
    arrow and must precede a ` |...|` extension tail; anywhere else inside a side it is refused.

    Components and how they may be grouped are the SMARTS reader's: `A.B` says only "not bonded",
    `(A.B)` demands one molecule component and `(A).(B)` demands two.  An intramolecular template is
    written with the first spelling.

    Map numbers pair the sides. A number on both sides is one atom carried through; on the reactant
    side only, that atom is DELETED (`M` exempts it); on the product side only, the atom is created
    and the number pairs with nothing, which is logged.  V2's `:100` / `:200` leaving-group convention
    is not read here: absence says it, so a ported template drops those numbers.

    The product side is EXPLICIT-ONLY: an unstated property is the default, not the matched atom's
    value, so `[C:1]` builds a neutral carbon with derived hydrogens whatever it matched.  Because
    that is silent, every map number whose reactant side states a charge, isotope or radical the
    product side leaves unstated puts one line in `log`.  It is a line and not a refusal: neutralizing
    a cation is a legitimate thing for a template to mean.

    A product primitive is one of two things and never neither: it BUILDS (element, isotope, charge,
    radical, stereo sign, bond order) or it CHECKS the result (`D`, `h`, `H`, `x`, `z`, `r`, `R`,
    `M`, a ring bond) -- so a cyclization states its ring size as a product-side `r`.  Anything that
    cannot be placed is refused, as is a build primitive offered as one of several `,` alternatives,
    a negated one, and a field stated twice.  A product atom either names an element or pairs with a
    reactant atom to inherit one; `[A:1]` says that inheritance out loud.

    `log` is a list to append lines to; omitting it discards them.  Lines a side's own lexer produced
    are prefixed with that side.

    `rule_id` is what a log record this template writes will call it.  A string read by hand has no
    name but its string, which is the default (`'smirks:C>>C'`); a template a table row produced is
    named after the row (`'reactions:13'`) because the row, not the composed string, is what its
    author can go and edit.  A composer that assembles one string out of several table rows is the
    only caller with something better to say than the default.

    Raises `IncorrectSmirks` for syntax.  The three-part `reactants>agents>products` form is refused
    by design: an agent is matched and never patched, which is a third semantics this notation does
    not have.  Offsets in a message about one side are offsets into that side.
    """
    cdef bytes raw
    cdef bytes body
    cdef bytes tail
    cdef bytes react_src
    cdef bytes prod_src
    cdef object exc                      # the except-as target; `warn.undeclared` counts it
    if isinstance(text, str):
        try:
            raw = (<str> text).encode('ascii')
        except UnicodeEncodeError:
            raise IncorrectSmirks('the string contains a non-ASCII character') from None
    elif isinstance(text, bytes):
        raw = <bytes> text
    else:
        raise TypeError('read_smirks takes a str or bytes')
    raw = raw.strip()

    cdef list parts = raw.rsplit(None, 1)
    body = raw
    tail = b''
    if len(parts) == 2 and (<bytes> parts[1]).startswith(b'|'):
        body = (<bytes> parts[0]).rstrip()
        tail = <bytes> parts[1]

    cdef int arrows = body.count(b'>')
    cdef int first = body.find(b'>')
    cdef int last = body.rfind(b'>')
    if arrows == 0:
        raise IncorrectSmirks('no `>>`: this names one side only, and one side is a pattern -- read '
                              'it with read_smarts')
    elif arrows == 1:
        raise IncorrectSmirks('a single `>` at position %d is not the SMIRKS arrow; write `>>`'
                              % first)
    elif arrows > 2:
        raise IncorrectSmirks('%d `>` characters: a SMIRKS has exactly one `>>`' % arrows)
    elif last != first + 1:
        raise IncorrectSmirks('the three-part form `reactants>agents>products` is not read: an agent '
                              'is matched and never patched, which is a third semantics this '
                              'notation does not have.  Name the agents on both sides, or leave '
                              'them out')

    react_src = body[:first].strip()
    prod_src = body[first + 2:].strip()

    # Refused here rather than left to the lexer, whose answer would be `unexpected ' ' ... a query
    # primitive belongs inside a bracket` -- true of a space and no help at all about what to do.
    cdef const char *b
    cdef uint32_t i, k, m
    for i in range(2):
        body = react_src if i == 0 else prod_src
        b = <const char *> body
        m = <uint32_t> len(body)
        for k in range(m):
            if b[k] <= 32:
                raise IncorrectSmirks('the %s side holds whitespace at position %d; only the arrow '
                                      'and the extension tail may be spaced'
                                      % ('reactant' if i == 0 else 'product', k))

    cdef object mylog = log if log is not None else []
    cdef QueryContainer reactants = QueryContainer()
    cdef QueryContainer products = QueryContainer()
    smk_side(reactants, react_src, 'reactant', mylog)
    smk_side(products, prod_src, 'product', mylog)

    cdef uint32_t n_react = reactants.atom_count
    cdef uint32_t n_prod = products.atom_count
    cdef set prod_radicals = set()
    # The bracket groups FIRST, so a tail naming an atom that already carries one is a contradiction
    # the tail can see rather than one it silently wins.  The reactant side's brackets are not read
    # here at all: the seal below refuses that op, which is the one refusal both entry points share.
    cdef dict prod_stereo_groups = smk_journal_stereo_groups(products)
    if tail:
        if tail.endswith(b'|') and len(tail) > 1:
            smk_cx(reactants, tail, n_react, n_prod, prod_radicals, prod_stereo_groups, mylog)
        else:
            mylog.append(mc_record('smirks:unterminated-tail', (),
                                   'the extension block after the SMIRKS is not terminated and was ignored: %s'
                                   % tail.decode('ascii', 'replace'),
                                   mc_lost()))

    # sealed AFTER the tail, because a `^N:` field appends reactant primitives
    try:
        reactants.atom_count_sealed()
    except ValueError as exc:
        # by now the byte offsets are gone and the atom number is what a template author can act on
        raise IncorrectSmirks('the reactant side does not compile: %s' % exc) from None

    cdef dict react_maps = smk_journal_maps(reactants)
    cdef dict prod_maps = smk_journal_maps(products)
    cdef dict seen
    cdef object side, sid, number         # dict-iteration targets
    for i in range(2):
        seen = {}
        for sid, number in (react_maps if i == 0 else prod_maps).items():
            side = 'reactant' if i == 0 else 'product'
            if number in seen:
                raise IncorrectSmirks('map number %d is on two atoms of the %s side; a map number '
                                      'names one atom per side' % (number, side))
            seen[number] = sid

    # a comprehension would have its own scope, where `warn.undeclared` cannot see these two
    cdef dict by_number = {}
    for sid, number in react_maps.items():
        by_number[number] = sid
    cdef dict pairs = {}
    for sid, number in prod_maps.items():
        if number in by_number:
            pairs[number] = (<object> by_number[number], sid)
        else:
            mylog.append(mc_record('smirks:product-only-map', (),
                                   'map number %d is on the product side only, so it pairs with nothing; that '
                                   'atom is created' % number))

    # N4: every product primitive is a thing to build or a thing to check, and one that is neither
    # is refused HERE rather than ignored at patch time.  See the comment above
    # `smk_classify_products` -- this is what makes dead product-side surface impossible.
    cdef set paired = set()
    for number in pairs:
        paired.add((<tuple> pairs[number])[1])
    cdef dict atom_build = {}
    cdef dict atom_check = {}
    cdef list prod_bonds = []
    cdef dict bond_build = {}
    cdef dict bond_check = {}
    cdef set inherited = set()
    smk_classify_products(products, prod_maps, paired, atom_build, atom_check, prod_bonds,
                          bond_build, bond_check, inherited)
    # What the product side says about a configuration, and the refusals a statement that says nothing
    # sound gets.  Read time, so a template carrying one fails at the string.  The groups are already
    # collected, brackets and tail both, which is what lets the correlated set be found from the string
    # alone.
    cdef tuple prod_stereo = smk_stereo_reading(prod_maps, pairs, prod_stereo_groups, atom_build,
                                                bond_build, prod_bonds)
    cdef dict prod_correlated = <dict> prod_stereo[0]
    cdef set prod_keep = <set> prod_stereo[1]
    cdef set prod_invert = <set> prod_stereo[2]
    # the one ABSOLUTE configuration a product side may state, after `@=`/`@~` are known so that a
    # string stating both gets the refusal rather than one of them silently
    cdef dict prod_geometry = smk_directions(prod_maps, smk_journal_directions(products), prod_bonds,
                                             bond_build, prod_keep, prod_invert)

    # N1, the mapped-pair lint: the price of explicit-only product semantics, paid in log lines
    cdef dict react_stated = smk_journal_stated(reactants)
    cdef dict prod_stated = smk_journal_stated(products)
    cdef int lost, said, kept
    cdef tuple pair
    for number in sorted(pairs):
        pair = <tuple> pairs[number]
        said = react_stated.get(pair[0], 0)
        kept = prod_stated.get(pair[1], 0)
        lost = said & ~kept
        if lost & SMK_STATED_CHARGE:
            mylog.append(mc_record('smirks:implicit-neutral-charge', (),
                                   'map number %d states a charge on the reactant side and none on the '
                                   'product side, so the product atom is neutral' % number))
        if lost & SMK_STATED_ISOTOPE:
            mylog.append(mc_record('smirks:implicit-no-isotope', (),
                                   'map number %d states an isotope on the reactant side and none on the '
                                   'product side, so the product atom has no mass number' % number))
        if lost & SMK_STATED_RADICAL:
            mylog.append(mc_record('smirks:implicit-no-radical', (),
                                   'map number %d states a radical on the reactant side and none on the '
                                   'product side, so the product atom is not a radical' % number))

    # Deletion BY ABSENCE, and by nothing else: a reactant atom survives when its map number pairs,
    # and `M` exempts an atom named purely as context from being removed.  An unmapped reactant atom
    # pairs with nothing, so `0` is never a key of `pairs` and it falls to the same rule.
    cdef frozenset masked = reactants.masked_atoms()
    cdef set deleted = set()
    cdef object number_of
    for sid in range(1, n_react + 1):
        if sid in masked:
            continue
        number_of = react_maps.get(sid, 0)
        if number_of not in pairs:
            deleted.add(sid)
    cdef set created = set()
    for sid in range(1, n_prod + 1):
        number_of = prod_maps.get(sid, 0)
        if number_of not in pairs:
            created.add(sid)

    cdef ReactionTemplate t = ReactionTemplate.__new__(ReactionTemplate)
    t.smirks = raw.decode('ascii')
    t.rule_id = 'smirks:%s' % t.smirks if rule_id is None else <str> rule_id
    t.reactants = reactants
    t._products = products
    t.reactant_map_numbers = react_maps
    t.product_map_numbers = prod_maps
    t.mapped_pairs = pairs
    t.deleted_atoms = frozenset(deleted)
    t.created_atoms = frozenset(created)
    t.reactant_bonds = smk_journal_bonds(reactants)
    t.product_stereo_groups = prod_stereo_groups
    t.product_radicals = frozenset(prod_radicals)
    t.product_atom_count = n_prod
    t.product_bonds = tuple(prod_bonds)
    t.product_atom_build = atom_build
    t.product_atom_check = atom_check
    t.product_bond_build = bond_build
    t.product_bond_check = bond_check
    t.product_inherited_elements = frozenset(inherited)
    t.product_stereo_correlated = prod_correlated
    t.product_stereo_keep = frozenset(prod_keep)
    t.product_stereo_invert = frozenset(prod_invert)
    t.product_stereo_geometry = prod_geometry
    return t
