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
# Kekulisation: an aromatic edge set becomes bond orders 1 and 2.
#
# WHY THIS IS AN OPERATION AND NOT SOMETHING A READER DOES ON THE WAY IN
#
# Aromatic bonds are stored.  A string that says aromatic is stored aromatic, a file that says
# aromatic is stored aromatic, and `kekule()` is one of exactly two operations in the library
# allowed to change a molecule's representation (`thiele()` is the other).  No reader calls it --
# not to make a downstream computation work, not for convenience.  A molecule's representation is
# therefore a fact a caller can read and change deliberately, never something a parser decided on
# its behalf.
#
# WHY THIS IS A CORE COMPONENT AND NOT A PARSER'S PRIVATE HELPER
#
# Four callers need the same conversion for four unrelated reasons: InChI, which must hand libinchi
# Kekule orders; MDL and MRV, whose writers may need one; and any caller that wants a Kekule form
# of its own.  Written inside a parser it would be written four times and would disagree with
# itself on some charged heteroatom.  Its symbols are `arom_*`, not `smi_*`, and nothing in its
# signatures knows what a tokenizer is.
#
# WHY A MATCHING AND NOT A PATH SEARCH
#
# In any Kekule form of an aromatic system, an atom that contributes one electron to the pi
# system carries exactly one double bond among its aromatic edges and an atom that donates a
# lone pair carries none.  The set of double bonds is therefore a MATCHING that saturates
# exactly the pi contributors.  Three things follow that a DFS over pi paths does not give:
#
#   * "no Kekule form exists" is a real answer, distinguishable from "not found down the paths
#     tried".
#   * "prefer pyridine over pyrrole" is not a heuristic and needs no retry budget: a may-match atom
#     that CAN take a ring double bond does, because phase 1 saturates every must-match atom and an
#     odd count of must-match atoms in a ring forces the may-match one in.
#   * the code is a classification table (`arom_classify`) plus a search (`arom_match`), with no
#     backtrack sites open-coded against a path list.
#
# The search is complete backtracking with a minimum-remaining-values order, not Edmonds'
# blossom.  Blossom is asymptotically better and, on the 5-to-30-atom aromatic systems that
# actually occur, measurably slower to write correctly than it saves.  Completeness is what the
# argument above needs, and backtracking has it.
#
# WHY THE AROMATIC EDGE SET IS STORED OR STATED, NEVER PERCEIVED
#
# The ordinary call passes no edge set and gets the bonds the molecule is holding: with aromaticity
# stored, "this bond is aromatic" is a fact rather than an inference.  A caller may still state a
# set explicitly, and then this file believes it -- because MDL bond type 4 is a per-bond statement
# and a vendor writer can mark a subset of a ring aromatic, or mark a bond that is in no ring at
# all.  Either way nothing here perceives: it reports, through the log, every place where believing
# the input required a repair.
#
# FAILURE POSTURE
#
# Everything recoverable is a log line plus the best available assignment; the caller gets the
# unresolved systems back as data.  `AromaticKekulizeError` is raised only for a broken internal
# invariant, so a caller can tell it from a real bug in this file.  Nothing here raises a bare
# ValueError.


# atom's role in the matching.  A "must" atom takes exactly one ring double bond, a "must-not"
# atom takes none, and a "may" atom is the pyrrole/pyridine free choice.
cdef enum:
    AROM_MUST_NOT = 0
    AROM_MUST = 1
    AROM_MAY = 2

# candidate tries per aromatic system before the search gives up.  A budget is not a tuning
# knob: exhausting it is reported through the log and the system comes back unresolved, never
# silently half-assigned and called done.
#
# A typed global rather than a DEF, for the reason `_canonical.pxi:82` gives for CANON_NO_SLOT:
# the search compares against it inside `nogil`, where a DEF's value would be a Python int.
cdef uint64_t AROM_NODE_BUDGET = 2000000

# "no atom" / "no candidate count yet" sentinels.  Typed globals for the same reason as the
# budget above: the search compares against them inside `nogil`, where a literal 0xFFFFFFFF
# would be a Python int.
cdef uint32_t AROM_NO_ATOM = 0xFFFFFFFF

# the nine elements that have an aromatic form at all; anything else in an aromatic ring is an
# input this file repairs and logs
cdef frozenset AROM_ELEMENTS = frozenset((5, 6, 7, 8, 15, 16, 33, 34, 52))

cdef tuple AROM_CLASS_NAMES = ('must_not', 'must', 'may')

# atomic numbers this file reasons about by name
cdef enum:
    AROM_Z_B = 5
    AROM_Z_C = 6
    AROM_Z_N = 7
    AROM_Z_O = 8
    AROM_Z_P = 15
    AROM_Z_S = 16
    AROM_Z_AS = 33
    AROM_Z_SE = 34
    AROM_Z_TE = 52

# `stated_h` sentinel: the input said nothing about this atom's hydrogen count.  Distinct from 0,
# because `[n]` is a free choice between pyrrole and pyridine while `[nH0]` is pyridine and
# nothing else -- see the reader's design note.
#
# NOT PORTABLE BEYOND FORMATS THAT HAVE A HYDROGEN CONVENTION.  Treating "unstated" as a free
# choice is sound for SMILES because the SMILES language itself supplies the missing count: bare
# lowercase `n` means no hydrogen and `[nH]` means one, so an unstated count is genuinely a range
# the writer of the string chose not to narrow.  A format with no such convention -- MDL bond type
# 4 on a ring of bare N atoms -- has no witness at all, and the right answer there is to refuse
# rather than to pick.  That refusal was asked for and upheld.  Anyone reusing this classifier for
# a new format must supply `stated_h` from something the format actually said, or accept a
# "may-match" answer that is a guess wearing a free choice's clothes.
DEF AROM_H_UNSTATED = -1


with cython.warn.undeclared(False):
    # bare so Python can import it, guarded so warn.undeclared stays quiet
    class AromaticKekulizeError(RuntimeError):
        """An invariant inside the kekuliser broke.

        NOT raised for bad input.  A non-ring aromatic bond, a partially aromatic ring, a
        hypervalent aromatic heteroatom and an aromatic system with no Kekule form are all
        recoverable: they produce a log line and, where possible, a partial assignment, and the
        unresolved systems come back from `kekulize` as data.  This exception means the code is
        wrong, which is why it is a distinct type -- a caller that stubs against `kekulize`
        needs to tell a real bug from an input it should have expected.
        """


cdef struct arom_scratch_t:
    void *block
    uint32_t *e_u          # [m]     edge endpoint, atom index
    uint32_t *e_v          # [m]
    uint8_t  *e_alive      # [m]     0 once an edge is pruned out of the aromatic set
    uint32_t *aptr         # [n + 1] CSR over the aromatic subgraph
    uint32_t *aadj         # [2m]    neighbour atom index
    uint32_t *aeid         # [2m]    edge id of that adjacency slot
    uint8_t  *adeg         # [n]     live aromatic degree
    uint8_t  *cls          # [n]     AROM_MUST / AROM_MAY / AROM_MUST_NOT
    uint32_t *nbrs         # [n]     `arom_classify_atom`'s `nbrs`: every bond except order 8
    uint8_t  *exo          # [n]     1 when a double or triple bond leaves the aromatic set
    int32_t  *mate         # [n]     matched partner atom index, -1 for unmatched
    int32_t  *comp         # [n]     component label, -1 outside the aromatic subgraph
    uint32_t *clist        # [n]     one component's atoms
    uint32_t *sv           # [n]     search / dfs stack: atom
    uint32_t *ss           # [n]     search / dfs stack: next adjacency slot to try


cdef int arom_scratch_alloc(arom_scratch_t *sc, uint32_t n, uint32_t m) except -1:
    # one struct, one malloc, one check, one free -- RULES.md 5.2
    cdef size_t n_u32 = align8(<size_t> n * sizeof(uint32_t))
    cdef size_t n_i32 = align8(<size_t> n * sizeof(int32_t))
    cdef size_t n_u8 = align8(<size_t> n * sizeof(uint8_t))
    cdef size_t m_u32 = align8(<size_t> m * sizeof(uint32_t))
    cdef size_t m_u8 = align8(<size_t> m * sizeof(uint8_t))
    cdef size_t m2_u32 = align8(<size_t> 2 * m * sizeof(uint32_t))
    cdef size_t ptr_len = align8((<size_t> n + 1) * sizeof(uint32_t))

    cdef size_t total = (2 * m_u32 + m_u8 + ptr_len + 2 * m2_u32 + 3 * n_u8
                         + 2 * n_i32 + 4 * n_u32)
    cdef char *block = <char *> PyMem_Malloc(total)
    if block is NULL:
        raise MemoryError('kekulisation scratch allocation failed')
    memset(block, 0, total)
    sc.block = <void *> block

    cdef size_t off = 0
    sc.e_u = <uint32_t *> (block + off); off += m_u32
    sc.e_v = <uint32_t *> (block + off); off += m_u32
    sc.e_alive = <uint8_t *> (block + off); off += m_u8
    sc.aptr = <uint32_t *> (block + off); off += ptr_len
    sc.aadj = <uint32_t *> (block + off); off += m2_u32
    sc.aeid = <uint32_t *> (block + off); off += m2_u32
    sc.adeg = <uint8_t *> (block + off); off += n_u8
    sc.cls = <uint8_t *> (block + off); off += n_u8
    sc.nbrs = <uint32_t *> (block + off); off += n_u32
    sc.exo = <uint8_t *> (block + off); off += n_u8
    sc.mate = <int32_t *> (block + off); off += n_i32
    sc.comp = <int32_t *> (block + off); off += n_i32
    sc.clist = <uint32_t *> (block + off); off += n_u32
    sc.sv = <uint32_t *> (block + off); off += n_u32
    sc.ss = <uint32_t *> (block + off); off += n_u32
    return 0


cdef void arom_csr(arom_scratch_t *sc, uint32_t n, uint32_t m) noexcept nogil:
    """Build the aromatic subgraph's CSR from the edge list, and set the live degrees."""
    cdef uint32_t i, k
    cdef uint32_t *aptr = sc.aptr
    for i in range(n + 1):
        aptr[i] = 0
    for k in range(m):
        aptr[sc.e_u[k] + 1] += 1
        aptr[sc.e_v[k] + 1] += 1
    for i in range(n):
        aptr[i + 1] += aptr[i]
        sc.adeg[i] = 0
    # a second cursor is not needed: adeg doubles as the fill offset while it is being built
    for k in range(m):
        i = sc.e_u[k]
        sc.aadj[aptr[i] + sc.adeg[i]] = sc.e_v[k]
        sc.aeid[aptr[i] + sc.adeg[i]] = k
        sc.adeg[i] += 1
        i = sc.e_v[k]
        sc.aadj[aptr[i] + sc.adeg[i]] = sc.e_u[k]
        sc.aeid[aptr[i] + sc.adeg[i]] = k
        sc.adeg[i] += 1


cdef void arom_prune_acyclic(arom_scratch_t *sc, Structure structure, uint32_t m) noexcept nogil:
    """Kill every stated aromatic edge that lies on no cycle OF THE MOLECULE.

    An aromatic bond outside every ring cannot carry ring pi electrons, so it is a single bond
    however it was written.  That is the one rule behind three unrelated inputs: the biphenyl
    spelled `c1ccccc1c2ccccc2`, MDL bond type 4 on an acyclic bond, and an aromatic substituent
    bond.

    RING MEMBERSHIP IS A PROPERTY OF THE MOLECULE, NOT OF THE AROMATIC EDGE SUBSET.  Running Tarjan
    over the aromatic subgraph and calling every bridge of THAT graph acyclic asks a different
    question: `c1c-cccc1` states five aromatic bonds that form a path, every one of them a bridge of
    the path and every one of them inside the same six-ring, so that test throws the whole ring away
    and answers cyclohexane instead of benzene.  `c1ccccc1` formally has six single bonds too.  The
    aromatic edge set comes from atom case (or from a per-bond aromatic mark), never from bond order,
    so an explicit `-` inside an all-lowercase ring is a preference WITHIN the pi system and not a
    wall around it: it says "not this bond", and the matching is what decides which bonds carry the
    double bonds anyway.

    `HE_IN_RING` is exactly the datum wanted -- `mark_bridges` sets it on every half-edge that
    lies on a cycle -- and it costs one lookup per stated edge instead of a graph search.
    Biphenyl stays correct for the right reason rather than by accident: its inter-ring bond is on
    no cycle of the MOLECULE either, so it is still pruned and still logged.

    NOTHING IS PROMOTED INTO THE AROMATIC SKELETON HERE, meaning a non-aromatic bond of an SSSR ring
    all of whose atoms carry an aromatic bond -- `c1cccc-c-1` and `c1cc-c-cc1`.  Promotion cannot
    live here: this function is handed an EDGE SET, so it cannot know that an atom with no aromatic
    edge was nonetheless written lowercase.  Only the reader knows that, and a reader that wants the
    promoted reading can state the promoted set.  What arrives is believed; what is unbelievable is
    logged.
    """
    cdef uint32_t k
    cdef halfedge_t *e
    for k in range(m):
        e = csr_find(structure, sc.e_u[k], sc.e_v[k])
        # NULL cannot happen -- `arom_setup` rejects a stated pair that is not a bond, and the
        # stored set is read out of this same CSR -- but a deref here would be a segfault
        if e is NULL or not (e.flags & HE_IN_RING):
            sc.e_alive[k] = 0


cdef void arom_relive(arom_scratch_t *sc, uint32_t n) noexcept nogil:
    """Recount live aromatic degrees after pruning."""
    cdef uint32_t i, k
    for i in range(n):
        sc.adeg[i] = 0
        for k in range(sc.aptr[i], sc.aptr[i + 1]):
            if sc.e_alive[sc.aeid[k]]:
                sc.adeg[i] += 1


cdef uint8_t arom_classify_atom(uint32_t element, int charge, bint radical, uint32_t nbrs,
                                bint exo_double, int stated_h, uint8_t *invalid) noexcept nogil:
    """The one-atom classification: does this atom take a ring double bond?

    NO CASE RAISES.  A reader that must accept whatever another tool emitted cannot treat an
    impossible aromatic atom as an error, so each of the twenty-odd invalid cases sets `invalid[0]`
    -- the caller logs it by name rather than letting it pass silently -- and returns AROM_MUST_NOT.

    AROM_MUST_NOT is the right degradation and AROM_MUST is not, which is the one judgement call
    in this function.  A spurious must-match atom can make a system that does have a Kekule form
    come back unresolved, taking the rest of the ring down with it; a spurious must-not-match
    atom only over-hydrogenates the one atom that was already wrong.  Errors stay local.

    `nbrs` counts every bond of any order except 8, explicit hydrogens included, aromatic and not.
    `exo_double` is a double or triple bond outside the aromatic set -- a quinone carbonyl, an
    N-oxide, an aromatic-written pyridone -- which spends the atom's pi electron elsewhere and so
    saturates it.  `stated_h` is AROM_H_UNSTATED unless the input gave a count.
    """
    if exo_double:
        return AROM_MUST_NOT

    if element == AROM_Z_C:
        if charge == 0:
            if nbrs != 2 and nbrs != 3:
                invalid[0] = 1
                return AROM_MUST_NOT
            return AROM_MUST
        if charge == 1 or charge == -1:
            if radical:
                if nbrs == 2:
                    return AROM_MUST_NOT
                invalid[0] = 1
                return AROM_MUST_NOT
            if nbrs == 3:
                return AROM_MUST_NOT
            if nbrs == 2:
                return AROM_MAY               # benzene cation/anion, or a charged pyrrole
            invalid[0] = 1
            return AROM_MUST_NOT
        invalid[0] = 1
        return AROM_MUST_NOT
    if element == AROM_Z_N or element == AROM_Z_P or element == AROM_Z_AS:
        if charge == 0:
            if radical:
                if nbrs != 2:                 # only a pyrrole radical is meaningful
                    invalid[0] = 1
                return AROM_MUST_NOT
            if nbrs == 3:
                # N with three neighbours is pyrrole and nothing else; P and As can be P(III)
                # or P(V)H, so for them the choice is still open
                return AROM_MUST_NOT if element == AROM_Z_N else AROM_MAY
            if nbrs == 2:
                if stated_h == AROM_H_UNSTATED:
                    return AROM_MAY           # pyrrole or pyridine, the classic free choice
                if stated_h == 0:
                    return AROM_MUST          # pyridine, stated
                if stated_h == 1:
                    return AROM_MUST_NOT      # pyrrole NH, stated
                invalid[0] = 1                # too many hydrogens for an aromatic ring
                return AROM_MUST_NOT
            if nbrs == 4 and element != AROM_Z_N:
                return AROM_MUST              # P(V) in ring, [P;a](-R1)-R2
            invalid[0] = 1
            return AROM_MUST_NOT
        if charge == -1:
            if nbrs != 2 or radical:
                invalid[0] = 1
                return AROM_MUST_NOT
            return AROM_MUST_NOT              # pyrrolide
        if charge == 1:
            if radical:
                if nbrs != 2:                 # not a pyridine cation-radical
                    invalid[0] = 1
                    return AROM_MUST_NOT
                return AROM_MUST
            if nbrs == 2:
                return AROM_MAY               # pyrrole cation or protonated pyridine
            if nbrs == 3:
                return AROM_MUST              # pyridinium, pyridine N-oxide
            invalid[0] = 1
            return AROM_MUST_NOT
        invalid[0] = 1
        return AROM_MUST_NOT
    if element == AROM_Z_O:
        if nbrs != 2:
            invalid[0] = 1
            return AROM_MUST_NOT
        if charge == 0:
            if radical:
                invalid[0] = 1
            return AROM_MUST_NOT              # furan
        if charge == 1:
            return AROM_MUST_NOT if radical else AROM_MUST     # pyrylium
        invalid[0] = 1
        return AROM_MUST_NOT
    if element == AROM_Z_S or element == AROM_Z_SE or element == AROM_Z_TE:
        if nbrs == 2:
            if radical:
                if charge != 1:
                    invalid[0] = 1
                return AROM_MUST_NOT
            if charge == 0:
                return AROM_MUST_NOT          # thiophene
            if charge == 1:
                return AROM_MUST
            invalid[0] = 1
            return AROM_MUST_NOT
        if nbrs == 3:
            if radical:
                if charge:
                    invalid[0] = 1
                return AROM_MUST_NOT
            if charge == 1:
                return AROM_MUST_NOT
            if charge == 0:
                return AROM_MUST
            invalid[0] = 1
            return AROM_MUST_NOT
        invalid[0] = 1                        # hypervalent S, Se, Te in a ring
        return AROM_MUST_NOT
    if element == AROM_Z_B:
        if charge == 0:
            if nbrs == 2:
                if radical:
                    return AROM_MUST_NOT      # C=1O[B]OC=1
                if stated_h == AROM_H_UNSTATED:
                    return AROM_MAY           # b1ccccc1, C=1OBOC=1 or B1C=CC=N1
                if stated_h == 0:
                    return AROM_MAY
                if stated_h == 1:
                    return AROM_MUST_NOT      # C=1O[BH]OC=1 or [BH]1C=CC=N1
                invalid[0] = 1
                return AROM_MUST_NOT
            if radical:
                invalid[0] = 1
            return AROM_MUST_NOT
        if charge == 1:
            if nbrs == 2 and not radical:
                return AROM_MUST_NOT
            invalid[0] = 1
            return AROM_MUST_NOT
        if charge == -1:
            if nbrs == 2:
                if radical:
                    return AROM_MUST          # the anion-radical is benzene-like
                return AROM_MAY               # C=1O[B-]OC=1 or [bH-]1ccccc1
            if radical:
                return AROM_MUST_NOT          # C=1O[B-*](R)OC=1
            return AROM_MAY
        invalid[0] = 1
        return AROM_MUST_NOT
    invalid[0] = 1                            # not an element with an aromatic form at all
    return AROM_MUST_NOT


cdef bint arom_match(arom_scratch_t *sc, uint32_t *atoms, uint32_t count,
                     uint64_t *nodes) noexcept nogil:
    """Saturate every must-match atom of one aromatic system.  Complete, so False means no
    Kekule form exists rather than "not found".

    Minimum-remaining-values order: always branch on the unsaturated must-match atom with the
    fewest available partners, so a forced move is taken before a free one and a dead end is hit
    at the shallowest depth it exists at.
    """
    cdef uint32_t i, v, w, k, pick, avail, best
    cdef int32_t top = -1
    cdef bint found

    # the budget is the loop guard rather than a check at the bottom, so falling out of the loop
    # IS budget exhaustion and the caller's `nodes > AROM_NODE_BUDGET` test reads the same fact
    while nodes[0] <= AROM_NODE_BUDGET:
        # --- pick the next must-match atom, MRV
        pick = AROM_NO_ATOM
        best = AROM_NO_ATOM
        for i in range(count):
            v = atoms[i]
            if sc.cls[v] != AROM_MUST or sc.mate[v] >= 0:
                continue
            avail = 0
            for k in range(sc.aptr[v], sc.aptr[v + 1]):
                if not sc.e_alive[sc.aeid[k]]:
                    continue
                w = sc.aadj[k]
                if sc.mate[w] < 0 and sc.cls[w] != AROM_MUST_NOT:
                    avail += 1
            if avail < best:
                best = avail
                pick = v
                if avail == 0:
                    break
        if pick == AROM_NO_ATOM:
            return True                        # every must-match atom is saturated
        if best == 0:
            # dead end: unwind to the shallowest frame with an untried candidate
            found = False
            while top >= 0:
                v = sc.sv[top]
                w = <uint32_t> sc.mate[v]
                sc.mate[w] = -1
                sc.mate[v] = -1
                if sc.ss[top] < sc.aptr[v + 1]:
                    found = True
                    break
                top -= 1
            if not found:
                return False                   # the search space is exhausted: impossible
        else:
            top += 1
            sc.sv[top] = pick
            sc.ss[top] = sc.aptr[pick]

        # --- take the next candidate at the current frame
        while True:
            v = sc.sv[top]
            found = False
            for k in range(sc.ss[top], sc.aptr[v + 1]):
                if not sc.e_alive[sc.aeid[k]]:
                    continue
                w = sc.aadj[k]
                if sc.mate[w] < 0 and sc.cls[w] != AROM_MUST_NOT:
                    sc.ss[top] = k + 1
                    sc.mate[v] = <int32_t> w
                    sc.mate[w] = <int32_t> v
                    found = True
                    break
            if found:
                break
            top -= 1                           # this frame is spent
            if top < 0:
                return False
            v = sc.sv[top]
            w = <uint32_t> sc.mate[v]
            sc.mate[w] = -1
            sc.mate[v] = -1

        nodes[0] += 1
    return False


cdef void arom_extend(arom_scratch_t *sc, uint32_t *atoms, uint32_t count) noexcept nogil:
    """Greedily match the still-unsaturated may-match atoms to each other.

    Maximal, not maximum, and that is enough: the pyrrole/pyridine choice is already decided by
    `arom_match`, because an odd count of must-match atoms around a ring forces the may-match
    atom in.  What is left here is a may-may pair with no must-match atom to force it, where
    either answer is a valid Kekule form.
    """
    cdef uint32_t i, v, w, k
    for i in range(count):
        v = atoms[i]
        if sc.cls[v] != AROM_MAY or sc.mate[v] >= 0:
            continue
        for k in range(sc.aptr[v], sc.aptr[v + 1]):
            if not sc.e_alive[sc.aeid[k]]:
                continue
            w = sc.aadj[k]
            if sc.mate[w] < 0 and sc.cls[w] == AROM_MAY:
                sc.mate[v] = <int32_t> w
                sc.mate[w] = <int32_t> v
                break


cdef uint32_t arom_partial(arom_scratch_t *sc, uint32_t *atoms, uint32_t count) noexcept nogil:
    """The fallback for a system with no Kekule form: a maximal matching, must-match first.

    Returns the number of must-match atoms it could not saturate.  The point is that the caller
    still gets a molecule -- with the aromatic system as faithful as it can be made -- and a log
    line naming what did not work out, rather than an exception and no molecule at all.
    """
    cdef uint32_t i, v, w, k
    cdef uint32_t left = 0
    for i in range(count):
        sc.mate[atoms[i]] = -1
    for i in range(count):
        v = atoms[i]
        if sc.cls[v] != AROM_MUST or sc.mate[v] >= 0:
            continue
        for k in range(sc.aptr[v], sc.aptr[v + 1]):
            if not sc.e_alive[sc.aeid[k]]:
                continue
            w = sc.aadj[k]
            if sc.mate[w] < 0 and sc.cls[w] != AROM_MUST_NOT:
                sc.mate[v] = <int32_t> w
                sc.mate[w] = <int32_t> v
                break
    arom_extend(sc, atoms, count)
    for i in range(count):
        v = atoms[i]
        if sc.cls[v] == AROM_MUST and sc.mate[v] < 0:
            left += 1
    return left


cdef class _AromRun:
    """One kekulisation in flight: the scratch, the edge list and the log.

    A cdef class rather than a bare struct because the scratch must be freed on every exit path
    including an exception raised while classifying, and `__dealloc__` is the only place that is
    true without a `try/finally` around every caller.
    """
    cdef arom_scratch_t sc
    cdef uint32_t n
    cdef uint32_t m
    cdef list log

    def __cinit__(self):
        self.sc.block = NULL
        self.n = 0
        self.m = 0
        self.log = []

    def __dealloc__(self):
        if self.sc.block is not NULL:
            PyMem_Free(self.sc.block)
            self.sc.block = NULL


cdef class KekuleResult:
    """What `kekule()` did: `changed`, `log`, `unresolved`.

    Named attributes and not a tuple, because this return value has already grown once (it was
    `(log, unresolved)` for an hour) and a positional shape makes the next growth a breaking change
    at every call site.  Unpacking is deliberately not supported: a caller that writes
    `changed, log, unresolved = mol.kekule()` is the call site that breaks next time.
    """
    cdef readonly bint changed
    """False when no bond order moved: a molecule with no aromatic bonds, or a second call."""
    cdef readonly list log
    """One human-readable line per repair, empty when the input needed none.

    A COPY, not the only channel: the same records are on `mol.log`, which is where a composed pipeline
    reads them back.  This one is here because `.changed` and `.unresolved` need a result object anyway.
    """
    cdef readonly list unresolved
    """A tuple of stable ids per aromatic system with no Kekule form; empty when all were assigned."""

    def __repr__(self):
        return (f'KekuleResult(changed={bool(self.changed)}, log={self.log!r}, '
                f'unresolved={self.unresolved!r})')


cdef KekuleResult arom_result(MoleculeContainer mol, bint changed, list log, list unresolved):
    cdef KekuleResult r = KekuleResult.__new__(KekuleResult)
    r.changed = changed
    r.log = log
    r.unresolved = unresolved
    # UNCONDITIONAL.  The molecule stores what happened to it; `.changed` and `.unresolved` are why the
    # result object exists at all.  An empty `log` folds to nothing on its own -- it must not be a
    # BRANCH, because a branch is where "sometimes we record" comes back in.
    mol.log.absorb('kekule', log, rule='kekule')
    return r


cdef list arom_stored_bonds(MoleculeContainer mol):
    """The aromatic bonds the molecule is already holding, as `(i, j)` index pairs, i < j.

    `aromatic_bonds=None` means this set, which is what makes `kekule()` a self-contained
    operation: a caller that just wants a Kekule form does not have to tell the molecule what it
    is already storing.  An explicitly stated set stays available for a caller that wants a
    SUBSET kekulised -- a vendor file's per-bond aromatic mark, for instance.
    """
    cdef Structure structure = mol._structure
    cdef uint32_t n = structure.header.atom_count
    cdef uint32_t *ptr = csr_ptr(structure)
    cdef halfedge_t *edges = csr_edges(structure)
    cdef list pairs = []
    cdef uint32_t i, k
    for i in range(n):
        for k in range(ptr[i], ptr[i + 1]):
            # each undirected bond appears as two half-edges; the i < to test takes one of them
            if edges[k].order == 4 and i < edges[k].to:
                pairs.append((i, <uint32_t> edges[k].to))
    return pairs


cdef _AromRun arom_setup(MoleculeContainer mol, aromatic_bonds, stated_h):
    """Read the stated aromatic edge set, build the subgraph, prune it, classify it."""
    mol._require_clean()
    cdef Structure structure = mol._structure
    cdef uint32_t n = structure.header.atom_count
    cdef dict index_of = mol._index_of
    cdef list numbers = mol._numbers

    # dedupe first: the same bond stated twice is a caller's convenience, not an error
    cdef set seen = set()
    cdef list pairs = []
    cdef uint32_t ia, ib
    cdef object pair, sa, sb, key
    if aromatic_bonds is None:
        pairs = arom_stored_bonds(mol)
        aromatic_bonds = ()
    for pair in aromatic_bonds:
        sa, sb = pair
        if sa not in index_of:
            raise KeyError(sa)
        if sb not in index_of:
            raise KeyError(sb)
        ia = <uint32_t> index_of[sa]
        ib = <uint32_t> index_of[sb]
        if ia == ib:
            raise ValueError(f'aromatic bond {pair!r} is a self loop')
        if csr_find(structure, ia, ib) is NULL:
            raise KeyError((sa, sb))
        key = (ia, ib) if ia < ib else (ib, ia)
        if key in seen:
            continue
        seen.add(key)
        pairs.append(key)

    cdef _AromRun run = _AromRun.__new__(_AromRun)
    run.n = n
    run.m = <uint32_t> len(pairs)
    arom_scratch_alloc(&run.sc, n, run.m if run.m else 1)
    cdef arom_scratch_t *sc = &run.sc
    cdef uint32_t k
    for k in range(run.m):
        sc.e_u[k] = <uint32_t> pairs[k][0]
        sc.e_v[k] = <uint32_t> pairs[k][1]
        sc.e_alive[k] = 1
    if not run.m:
        return run

    arom_csr(sc, n, run.m)
    arom_prune_acyclic(sc, structure, run.m)
    for k in range(run.m):
        if not sc.e_alive[k]:
            run.log.append(mc_record('kekule:acyclic-bond',
                                     (numbers[sc.e_u[k]], numbers[sc.e_v[k]]),
                                     f'aromatic bond {numbers[sc.e_u[k]]}-{numbers[sc.e_v[k]]} is in no '
                                     f'ring; read as single',
                                     mc_repaired()))
    arom_relive(sc, n)

    # --- classify
    cdef atom_t *atoms = structure.atoms()
    cdef atom_t *a
    cdef uint32_t *ptr = csr_ptr(structure)
    cdef halfedge_t *edges = csr_edges(structure)
    cdef uint32_t i, j, nbrs
    cdef uint8_t order
    cdef uint8_t invalid
    cdef bint exo, triple, aromatic_slot
    cdef int h, ad, ch
    # `stated` MUST be declared.  An undeclared Cython local is answered with an inferred Python
    # object and a build warning, and nothing else notices -- see RULES.md §9.7.
    cdef object rad_note, h_note, stated
    cdef dict h_map = stated_h if stated_h is not None else {}
    for i in range(n):
        if not sc.adeg[i]:
            sc.cls[i] = AROM_MUST_NOT
            continue
        if sc.adeg[i] > 3:
            # widened to int for the message: `uint8_t` is `unsigned char`, and Cython formats a
            # char-typed value as a one-character string, so the degree would come out as a
            # control character rather than a number
            ad = sc.adeg[i]
            run.log.append(mc_record('kekule:hypercondensed', (numbers[i],),
                                     f'atom {numbers[i]} carries {ad} aromatic bonds; hypercondensed, read '
                                     f'as saturated',
                                     mc_repaired()))
            sc.cls[i] = AROM_MUST_NOT
            continue
        a = atoms + i
        nbrs = 0
        exo = False
        triple = False
        for j in range(ptr[i], ptr[i + 1]):
            order = edges[j].order
            if order == 8:
                continue                      # a dative bond is not a neighbour for this count
            nbrs += 1
            if order < 2 or order == 4:
                # A SINGLE BOND, OR AN AROMATIC ONE, AND NEITHER SPENDS THE PI ELECTRON ELSEWHERE.
                # `order == 4` is the load-bearing half.  A stored aromatic bond that is not a LIVE
                # edge of the aromatic set -- pruned for lying on no cycle, or simply outside the set
                # a caller stated -- is already announced "read as single" by the prune's own log
                # line, and reading it as a double here would contradict that line in the same
                # function.  Biphenyl written the way OpenSMILES requires it be read,
                # `c1ccccc1c1ccccc1`, is the measured case: the inter-ring bond is on no cycle, so it
                # is pruned and logged, and each ipso carbon would then look like a quinone carbonyl
                # carbon and answer must-not.  Five must-match atoms left in a six-ring is an odd
                # count with no perfect matching, so BOTH rings would come back unresolved.
                # Only a stated double or triple bond saturates an aromatic atom; order 4 is neither.
                continue
            aromatic_slot = False
            for k in range(sc.aptr[i], sc.aptr[i + 1]):
                if sc.aadj[k] == edges[j].to and sc.e_alive[sc.aeid[k]]:
                    aromatic_slot = True
                    break
            if not aromatic_slot:
                exo = True
                if order == 3:
                    triple = True
        if triple:
            run.log.append(mc_record('kekule:triple-bond-exo', (numbers[i],),
                                     f'atom {numbers[i]} has a triple bond into an aromatic ring; read as '
                                     f'saturated',
                                     mc_repaired()))
        # THE STORED NIBBLE IS THE DEFAULT SOURCE AND `stated_h` ONLY OVERRIDES IT.  Forcing an atom
        # absent from the dict to AROM_H_UNSTATED throws away a count the arena is already holding;
        # on a five-ring with two nitrogens that freedom picks the wrong one and imidazole comes back
        # with a four-valent neutral N.  Reading the store still answers UNSTATED where nothing was
        # said -- a builder mid-flight has H_UNKNOWN in every nibble it has not written, so the store
        # answers UNSTATED for exactly those atoms.
        stated = h_map.get(numbers[i])
        if stated is None:
            h = AROM_H_UNSTATED if at_implicit_h_unknown(a) else at_implicit_h(a)
        else:
            h = stated
        # Bounded by the IMPLICIT count's maximum, 14, and not by the nibble's width.  `stated_h` is
        # an implicit hydrogen count, so 15 is not a count it may carry: that value is H_UNKNOWN in
        # `atom_t.hydrogens`, and a caller who means "the input said nothing" already has
        # AROM_H_UNSTATED.  Bounding by H_NIBBLE_MAX here let a 15 through and classified it as a
        # count, which is the one thing a sentinel must never be mistaken for.
        if h != AROM_H_UNSTATED and (h < 0 or h > H_IMPLICIT_MAX):
            raise ValueError(f'stated_h[{numbers[i]}] = {h!r} is outside 0..{H_IMPLICIT_MAX}')
        # KEPT, not recomputed later.  `arom_h_candidates` asks the classification table the same
        # question about the same atom with one hydrogen instead of none, and the two arguments it
        # would otherwise have to derive a second time are these.  Deriving them twice is how the
        # two copies drift; `nbrs` in particular is not the arena's degree and not `adeg` either --
        # it skips order 8 and counts explicit hydrogens -- so a second spelling of it would be a
        # second place to get that wrong.
        sc.nbrs[i] = nbrs
        sc.exo[i] = 1 if exo else 0
        invalid = 0
        sc.cls[i] = arom_classify_atom(a.element, a.charge, at_radical(a), nbrs, exo, h, &invalid)
        if invalid:
            if a.element not in AROM_ELEMENTS:
                run.log.append(mc_record('kekule:non-aromatic-element', (numbers[i],),
                                         f'atom {numbers[i]} is {symbol_of(a)}, which has no '
                                         f'aromatic form; its aromatic bonds are read as single',
                                         mc_repaired()))
            else:
                rad_note = ', radical' if at_radical(a) else ''
                h_note = '' if h == AROM_H_UNSTATED else f', {h}H'
                ch = a.charge          # `int8_t` is `signed char`; widened as `ad` above is
                run.log.append(mc_record('kekule:invalid-aromatic-state', (numbers[i],),
                                         f'atom {numbers[i]} ({symbol_of(a)}, charge {ch}, '
                                         f'{nbrs} neighbours{rad_note}{h_note}) is not a valid aromatic '
                                         f'state; read as saturated',
                                         mc_repaired()))
    for i in range(n):
        sc.mate[i] = -1
    return run


cdef bint arom_is_ring_bond(arom_scratch_t *sc, uint32_t i, uint32_t j) noexcept nogil:
    """Is the bond `i`-`j` a live edge of the aromatic subgraph?

    NOT `order == 4`.  The aromatic edge set is what the caller stated and what pruning left of it,
    which is not the same as what the arena stores: a caller may hand in a subset of a ring, and a
    bond stated aromatic but lying outside every ring has already been dropped from the set while
    keeping its stored order.  Every question of the form "is this bond in the ring system" has to
    be asked of `sc`, and asking the arena instead is the bug this helper exists to make hard.
    """
    cdef uint32_t k
    for k in range(sc.aptr[i], sc.aptr[i + 1]):
        if sc.aadj[k] == j and sc.e_alive[sc.aeid[k]]:
            return True
    return False


cdef uint32_t arom_oxide_exo(arom_scratch_t *sc, atom_t *atoms, uint32_t *ptr, halfedge_t *edges,
                             uint32_t i, uint8_t order, int exo_charge) noexcept nogil:
    """The exocyclic O or N hanging off aromatic ring nitrogen `i`, or `AROM_NO_ATOM`.

    Recognises an `[N;a;D3]` with one substituent: two live aromatic bonds and exactly one
    bond out of the ring system, of `order`, to a terminal `[O;D1]` or `[N;D1,D2]` whose formal
    charge is `exo_charge`.  Both atoms must be non-radical, and the partner must be outside every
    aromatic ring -- an `n`-`n` bond between two rings is a biaryl and not an N-imide.

    Nitrogen and nothing else.  `p(=O)` and `[as](=O)` are left alone deliberately: a phosphinine
    oxide's P really is pentavalent, so the classifier's P/As arms already accept the double bond.
    Sulfur is left alone for the opposite reason -- see the note in `arom_separate_charges`.
    """
    cdef atom_t *a = atoms + i
    if a.element != AROM_Z_N or at_radical(a) or sc.adeg[i] != 2:
        return AROM_NO_ATOM
    cdef uint32_t k, j
    cdef uint32_t found = AROM_NO_ATOM
    for k in range(ptr[i], ptr[i + 1]):
        if edges[k].order == 8:
            continue                          # dative, as the classifier also skips
        j = edges[k].to
        if arom_is_ring_bond(sc, i, j):
            continue
        if found != AROM_NO_ATOM:
            return AROM_NO_ATOM               # two substituents: D4, not the shape
        if edges[k].order != order:
            return AROM_NO_ATOM
        found = j
    if found == AROM_NO_ATOM:
        return AROM_NO_ATOM                   # D2: a plain ring nitrogen
    cdef atom_t *b = atoms + found
    if b.charge != exo_charge or at_radical(b) or sc.adeg[found]:
        return AROM_NO_ATOM
    cdef uint32_t deg = 0
    for k in range(ptr[found], ptr[found + 1]):
        if edges[k].order != 8:
            deg += 1
    if b.element == AROM_Z_O:
        return found if deg == 1 else AROM_NO_ATOM
    if b.element == AROM_Z_N:
        return found if deg <= 2 else AROM_NO_ATOM
    return AROM_NO_ATOM


cdef list arom_separate_charges(MoleculeContainer mol, _AromRun run):
    """Rewrite `n(=O)` and `n(=N)` as the charge-separated N-oxide and N-imide, in place.

    THE ONE UNCONDITIONAL REPAIR IN THIS FILE, and it is unconditional because it needs no context
    to be right: a neutral aromatic nitrogen with two ring bonds has spent its lone pair on the
    ring, so there is nothing left to make a pi bond to a substituent with.  `O=n1ccccc1` is not an
    alternative spelling of pyridine N-oxide that this file happens to dislike -- it is a nitrogen
    with no valid electronic state, and every consumer downstream would have to decide what to do
    with it.  Deciding here, once, is what "repair belongs at the input boundary" means, and the
    shape is narrow enough to test directly: no pattern and no isomorphism.

    CHARGE IS CONSERVED, which is what separates this from the shifts in `arom_shift_candidates`:
    the ring nitrogen gains +1 and the substituent gains -1 in the same edit, so the molecule's
    total charge is untouched and the repair can never turn a neutral input into an ion.  The
    hydrogen counts are untouched for the same reason -- trading one bond order unit for one charge
    unit leaves the valence of both atoms exactly where it was, so `[O-]` needs no more hydrogens
    than `=O` did and `[N-]H` no more than `=NH`.  Nothing here has to know about hydrogens, which
    matters because this runs before `calc_implicit` on a molecule being built.

    THE SULFUR RULE IS DELIBERATELY ABSENT: `[S;a;D3;+]-[O;D1;-]` -> `S=O`.  It runs the
    separation BACKWARDS, and the spelling it destroys is the only one that could ever kekulise --
    a three-coordinate neutral S is the classifier's must-match sulfonium-ylide state, while the
    `S=O` it writes is a saturated must-not.  Applied to `[O-][s+]1ccccc1` it converts a system
    with no Kekule form into a different system with no Kekule form and loses the input's charges
    on the way.

    Returns one log line per repair; an empty list means nothing was touched, which is the common
    case and costs one pass over the atoms.
    """
    cdef Structure structure = mol._structure
    cdef arom_scratch_t *sc = &run.sc
    cdef atom_t *atoms = structure.atoms()
    cdef uint32_t *ptr = csr_ptr(structure)
    cdef halfedge_t *edges = csr_edges(structure)
    cdef list numbers = mol._numbers
    cdef list log = []
    cdef list edits = []
    cdef uint32_t i, exo
    cdef object pair
    for i in range(run.n):
        if not sc.adeg[i] or atoms[i].charge:
            continue
        exo = arom_oxide_exo(sc, atoms, ptr, edges, i, 2, 0)
        if exo == AROM_NO_ATOM:
            continue
        edits.append((numbers[i], numbers[exo]))
        log.append(mc_record('kekule:n-oxide-charge-sep',
                             (numbers[i], numbers[exo]),
                             f'aromatic nitrogen {numbers[i]} cannot carry a double bond to '
                             f'{symbol_of(&atoms[exo])} {numbers[exo]}; read charge-separated as '
                             f'{numbers[i]}(+)-{numbers[exo]}(-)',
                             mc_repaired()))
    if not edits:
        return log
    # REPRESENTATION, NOT STRUCTURE, for the reason `kekule`'s own emit gives: the charge-separated
    # form and the hypervalent one are the same molecule differently written, so a stored CIP
    # descriptor is not made false by the rewrite.  Only the order edit needs the flag -- `_apply`
    # does not treat a charge edit as invalidating -- but the flag covers the scope rather than the
    # op, because a scope that is exempt in part is not one anybody can reason about.
    mol._representation_change = True
    try:
        with mol.edit():
            for pair in edits:
                mol.set_charge(pair[0], 1)
                mol.set_charge(pair[1], -1)
                mol.set_order(pair[0], pair[1], 1)
    finally:
        mol._representation_change = False
    return log


cdef list arom_shift_candidates(arom_scratch_t *sc, atom_t *atoms, uint32_t *ptr,
                                halfedge_t *edges, uint32_t *comp_atoms, uint32_t count):
    """The charge shifts that could turn a must-not nitrogen of a failed system into a must-match one.

    Each candidate is `(ring index, substituent index, kind)`, kind 1 meaning "write -1 on the
    substituent and drop the bond to single" and kind 0 meaning "write +1 on the ring nitrogen".
    Two shapes, both of them an N-oxide written a bond order away from the charge-separated form
    that `arom_separate_charges` normalises to:

    * `[n+](=O)` -- the nitrogen is already cationic, so separating the double bond puts the -1 on
      the substituent with nothing to cancel it.  Furoxan arrives from real files spelled this way.
    * `n(-[O-])` -- the anion is already on the substituent, so the nitrogen needs the +1.

    THESE ARE NOT CHARGE-CONSERVING, which is the whole reason they are candidates offered to a
    failed system rather than repairs applied to every molecule.  Each moves the molecule's total
    charge by one, and there are inputs where that would be flatly wrong: `O=[n+]1cccc[c-]1` is a
    valid neutral ylide spelling of pyridine N-oxide whose ring HAS a Kekule form, and shifting its
    substituent to `[O-]` would hand back an anion.  Gating on "no Kekule form as written" is what
    tells the two apart without a single pattern, and it does so on the property that actually
    matters: applied to a five-ring, `[O-]n1cccc1` is a pyrrol-1-olate that kekulises as written
    and is left alone, while the same shape in a six-ring cannot kekulise at all and is repaired.
    Ring size never appears here; it does not have to.
    """
    cdef list out = []
    cdef uint32_t i, v, exo
    for i in range(count):
        v = comp_atoms[i]
        if sc.cls[v] != AROM_MUST_NOT:
            continue                          # a may-match atom is already free
        if atoms[v].charge == 1:
            exo = arom_oxide_exo(sc, atoms, ptr, edges, v, 2, 0)
            if exo != AROM_NO_ATOM:
                out.append((v, exo, 1))
        elif atoms[v].charge == 0:
            exo = arom_oxide_exo(sc, atoms, ptr, edges, v, 1, -1)
            if exo != AROM_NO_ATOM:
                out.append((v, exo, 0))
    return out


cdef list arom_h_candidates(arom_scratch_t *sc, atom_t *atoms, uint32_t *comp_atoms,
                            uint32_t count):
    """The must-match atoms of a failed system that one added hydrogen would set free.

    THE OPENSMILES HOLE.  A lowercase aromatic `n` states no hydrogen count, and the specification
    is explicit that the ring decides it: with one hydrogen the nitrogen donates its lone pair and
    takes no ring double bond, without one it contributes a single pi electron and must take one.
    Toolkits write `c1cncn1` for imidazole and `c1ccnc1` for pyrrole all the time -- the hydrogen
    that makes those rings aromatic at all is simply not in the string.  Neither ring has a Kekule
    form as written, and refusing them would mean refusing a large fraction of the aromatic SMILES
    in the world.

    So the count is DERIVED here rather than demanded from the input, and this is the file's third
    relaxation after the two charge shifts: offered only to a system that has already failed, tried
    in the scratch, and written back only if it produced a complete matching.  A ring that
    kekulises as written never reaches this function, so a correctly spelled pyridine cannot be
    handed a hydrogen it does not want.

    THE PREDICATE IS THE CLASSIFICATION TABLE AND NOT A PATTERN.  A candidate is any must-match
    atom that `arom_classify_atom` calls must-NOT-match, and valid, when asked again with
    `stated_h = 1`.  Today exactly one arm answers that way -- neutral N, P or As with two
    neighbours and a stated zero -- which is precisely the pyridine-versus-pyrrole choice and
    nothing else.  Spelling that arm out here instead would be a second copy of it, and the copy
    would be the one that goes stale; asking the table means a new arm is covered the day it is
    written.  It also rules out by construction every must-match state the table decides without
    consulting hydrogens: pyrylium O(+), thiophenium S(+), the sulfonium ylide, pyridinium and the
    N-oxide cation all keep their class under the question and so are never offered a hydrogen.

    Returns atom indices.  The caller relaxes them to AROM_MAY rather than to AROM_MUST_NOT, which
    is what keeps the repair minimal: a may-match atom is a partner the search uses when a ring
    forces it and leaves alone when it does not, so only the nitrogens that the matching could not
    reach end up with the hydrogen.  Pinning them must-not instead would need a search over subsets
    to find the fewest, and would still get pyrimidine wrong.
    """
    cdef list out = []
    cdef uint32_t i, v
    cdef uint8_t invalid
    cdef uint8_t relaxed
    for i in range(count):
        v = comp_atoms[i]
        if sc.cls[v] != AROM_MUST:
            continue
        invalid = 0
        relaxed = arom_classify_atom(atoms[v].element, atoms[v].charge, at_radical(atoms + v),
                                     sc.nbrs[v], sc.exo[v], 1, &invalid)
        if not invalid and relaxed == AROM_MUST_NOT:
            out.append(v)
    return out


cdef list arom_cation_candidates(arom_scratch_t *sc, atom_t *atoms, uint32_t *ptr,
                                 halfedge_t *edges, uint32_t *comp_atoms, uint32_t count):
    """The must-not atoms of a failed system that one added POSITIVE charge would set free.

    THE MISSING FORMAL CHARGE, and it is the single commonest defect in aromatic SMILES coming out
    of registration systems: `c1ccn(C)cc1` for N-methylpyridinium, `NC(=O)c1cccn(C)c1` for the
    N-methylnicotinamide half of NAD(+), `c1ccn(N)cc1` for 1-aminopyridinium.  A three-coordinate
    neutral aromatic nitrogen has no room for a ring double bond, so five carbons are left with an
    odd count and nothing pairs off.  The ring is not asking for a hydrogen -- there is nowhere to
    put one -- and not asking for a bond order.  It is asking for the cation the third substituent
    already implies, which the writer left out because lowercase `n` looked like enough.

    This shape and the surplus hydrogen below are the two the relaxations exist for; the net holding
    them is `CORPUS_SHAPES` in `chemistry/test/test_isomers.py`, public molecules one per shape.

    THE PREDICATE IS THE CLASSIFICATION TABLE AND NOT A PATTERN, as in `arom_h_candidates`: a
    candidate is any must-not atom that `arom_classify_atom` calls must-MATCH, and valid, when asked
    again at `charge + 1`.  Asking rather than spelling the arm out is what keeps this covered when
    the table grows one, and it rules out by construction every state the table answers `may` for --
    a two-coordinate nitrogen at +1 is a pyridinium that may take a double bond or not, so
    `c1cncn1` is never offered a charge and gets the hydrogen it actually wants instead.

    Each candidate is `(index, new charge, new hydrogen count, substituent index or AROM_NO_ATOM,
    kind)`, kind 0 meaning "the cation alone" and kind 1 "the cation and its substituent's proton",
    and the values are what to WRITE rather than deltas -- the emit runs inside an edit session
    where the arena may not be read.  The substituent is a terminal neutral O or N carrying exactly
    one hydrogen: with one there, `[n+][O-]` is charge conserving and is pyridine N-oxide, which is
    what `c1ccn(O)cc1` means and what the caller tries first.

    NOT CHARGE CONSERVING at kind 0, which is why this is a candidate offered to a system that has
    already failed rather than a repair applied to every molecule.  `C[n+]1ccccc1` spelled correctly
    kekulises as written and never reaches here.
    """
    cdef list out = []
    cdef uint32_t i, v, exo
    cdef int h
    cdef uint8_t invalid, relaxed
    for i in range(count):
        v = comp_atoms[i]
        if sc.cls[v] != AROM_MUST_NOT or at_radical(atoms + v):
            continue
        if at_implicit_h_unknown(atoms + v):
            continue                          # nothing to preserve, and nothing to state
        h = at_implicit_h(atoms + v)
        invalid = 0
        relaxed = arom_classify_atom(atoms[v].element, atoms[v].charge + 1, 0,
                                     sc.nbrs[v], sc.exo[v], h, &invalid)
        if invalid or relaxed != AROM_MUST:
            continue
        exo = arom_oxide_exo(sc, atoms, ptr, edges, v, 1, 0)
        if exo != AROM_NO_ATOM and (at_implicit_h_unknown(atoms + exo)
                                    or at_implicit_h(atoms + exo) != 1):
            exo = AROM_NO_ATOM                # no proton to move, so the pair would not conserve
        out.append((v, atoms[v].charge + 1, h, exo, 0 if exo == AROM_NO_ATOM else 1))
    return out


cdef list arom_surplus_h_candidates(arom_scratch_t *sc, atom_t *atoms, uint32_t *comp_atoms,
                                    uint32_t count):
    """The must-not atoms of a failed system that giving up ONE hydrogen would set free.

    THE EXACT MIRROR OF `arom_h_candidates`, and the reason it exists is the symmetry: that function
    derives a hydrogen the notation never carried, and this one withdraws a hydrogen the notation
    carried wrongly.  `Cn1cc[nH]c1` states two two-electron donors in a five-ring, which leaves
    three carbons and an odd count; `c1cc[nH]cc1` states one in a six-ring, which leaves five.
    Neither is a molecule.  The N-methyl nitrogen cannot give anything up, the `[nH]` can, and once
    it does the first ring is 1-methylimidazole and the second is pyridine.

    A HYDROGEN THE INPUT DID STATE IS NOT PRIVILEGED, which is the answer `arom_h_candidates` already
    gives for a stated ZERO: honouring the input buys a half-assigned ring and a hypovalent atom
    instead of the one molecule the input could have meant.  A stated hydrogen is garbage input
    exactly as often as a stated zero is, and privileging one over the other is an asymmetry with no
    chemistry behind it.  The gate is what keeps it safe: `c1cc[nH]c1` kekulises as written, so
    pyrrole never reaches here.

    THE PREDICATE IS THE CLASSIFICATION TABLE, once more: a candidate is a must-not atom holding at
    least one hydrogen that the table calls must-MATCH, and valid, at `stated_h - 1`.  Today that is
    the neutral two-coordinate N, P or As arm and nothing else, so an anionic `[n-]`, a furan O and
    a thiophene S are all declined without a word about them here.

    Candidates share the caller's tuple shape -- `(index, new charge, new hydrogen count,
    AROM_NO_ATOM, 2)` -- so one loop applies both pools.
    """
    cdef list out = []
    cdef uint32_t i, v
    cdef int h
    cdef uint8_t invalid, relaxed
    for i in range(count):
        v = comp_atoms[i]
        if sc.cls[v] != AROM_MUST_NOT or at_radical(atoms + v):
            continue
        if at_implicit_h_unknown(atoms + v):
            continue                          # an unstated count is the other function's business
        h = at_implicit_h(atoms + v)
        if h < 1:
            continue
        invalid = 0
        relaxed = arom_classify_atom(atoms[v].element, atoms[v].charge, 0,
                                     sc.nbrs[v], sc.exo[v], h - 1, &invalid)
        if invalid or relaxed != AROM_MUST:
            continue
        out.append((v, atoms[v].charge, h - 1, AROM_NO_ATOM, 2))
    return out


cdef _AromRun arom_prepare(MoleculeContainer mol, aromatic_bonds, stated_h):
    """`arom_setup`, plus the unconditional charge separation and the reclassification it forces."""
    cdef _AromRun run = arom_setup(mol, aromatic_bonds, stated_h)
    cdef list repairs = arom_separate_charges(mol, run)
    if not repairs:
        return run
    # A charge and a bond order moved, so every class that was computed from them is stale.
    # Classifying the repaired molecule from scratch rather than patching the table in place is
    # both cheaper to get right and the only version that stays right when the table grows an arm.
    # The first pass's log is DISCARDED and not carried: it describes the molecule as it arrived,
    # including the "is not a valid aromatic state" line that the repair has just made untrue.  The
    # lines that are still true -- an acyclic aromatic bond, a hypercondensed atom -- are
    # regenerated identically, because the repair cannot change ring membership.
    run = arom_setup(mol, aromatic_bonds, stated_h)
    run.log = repairs + run.log
    return run


def kekule(MoleculeContainer mol not None, aromatic_bonds=None, stated_h=None):
    """Turn a set of aromatic bonds into Kekule orders 1 and 2.  `MoleculeContainer.kekule`.

    This is a DELIBERATE operation and one of the two in the library allowed to change a
    molecule's representation (`thiele` is the other).  No reader calls it: a string that says
    aromatic is stored aromatic, a string that says Kekule is stored Kekule, and the caller
    decides when to convert.  Silent kekulisation on input is the thing this design refuses.

    `aromatic_bonds` is an iterable of `(n, m)` pairs, or `None` for the aromatic
    bonds the molecule is already storing -- which is the ordinary call.  An explicit set is the
    aromatic edge set as the INPUT stated it, never as anything perceived it: MDL bond type 4 is a
    per-bond fact and a vendor file may mark a subset of a ring, or a bond in no ring at all, so
    re-perceiving would silently disagree with the file.  Duplicate pairs are ignored; a pair that
    is not a bond of `mol` raises `KeyError`.

    `stated_h` is an optional `{n: count}` for the atoms whose hydrogen count the input
    gave, and it OVERRIDES the arena rather than supplying what the arena lacks.  An absent key
    falls back to the stored count, and to "the input said nothing" only where the store itself
    says so -- a nibble of H_UNKNOWN, which is what a molecule being built has in every atom it has
    not written yet.  Pass it when you are holding a count the arena has not been told about; the
    ordinary call on a parsed molecule does not need it, because the parser has already stored what
    the string said.

    THIS FUNCTION REPAIRS, and charges and hydrogen counts are what it may rewrite.  A ring with no
    Kekule form as written is offered five relaxations in order:

    1. and 2. the two charge shifts for an N-oxide spelled a bond order away from its
       charge-separated form (`arom_shift_candidates`);
    3. a DERIVED hydrogen on an aromatic nitrogen whose count the notation never carried -- `c1cncn1`
       is imidazole and `c1ccnc1` is pyrrole, and both come back kekulised with the hydrogen the
       string omitted (`arom_h_candidates`);
    4. a MISSING FORMAL CHARGE -- `c1ccn(C)cc1` is N-methylpyridinium and comes back as the cation
       its third substituent implies, and `c1ccn(O)cc1` comes back as pyridine N-oxide
       (`arom_cation_candidates`);
    5. a SURPLUS HYDROGEN -- `Cn1cc[nH]c1` is 1-methylimidazole and `c1cc[nH]cc1` is pyridine, both
       with the hydrogen the ring cannot afford dropped (`arom_surplus_h_candidates`).

    NONE OF THE FIVE TOUCHES A RING THAT KEKULISES AS WRITTEN, which is the whole safety argument and
    the only one there is: correct input costs nothing, and a defect is told from a correct spelling
    by the property that actually matters rather than by a pattern.  So `c1cc[nH]c1` keeps its
    hydrogen and `C[n+]1ccccc1` keeps its charge.

    A HYDROGEN THE INPUT STATED IS NOT PRIVILEGED, on the same reading relaxation 3 applies to a
    stated ZERO: honouring either buys a half-assigned ring and a hypovalent atom instead of the one
    molecule the input could have meant, so `c1cc[nH]cc1` comes back as pyridine rather than
    unresolved.  What cannot be repaired is still reported: five aromatic carbons is an odd count no
    charge and no hydrogen makes even, so `c1cccc1` comes back in `unresolved`.

    IT ALSO HEALS THE HYDROGEN COUNTS ITS OWN ORDERS MADE DERIVABLE, and this is not one of the five
    relaxations: it touches only atoms that claim NO count at all.  A format with no hydrogen channel
    leaves the pyrrole-versus-pyridine nitrogen as H_UNKNOWN, because a local look cannot tell those
    apart -- and once the ring holds definite orders it can, so the count goes in.  Fill-only, so a
    count the reader stored from something no valence row reproduces (ferrocene, diborane's bridges)
    is not overwritten; skipped for the atoms of an `unresolved` system, whose order sum is deficient
    and would be filled with hydrogens the input never had; and not logged, because it is the
    read-time derivation finishing rather than a repair of anything the input stated.  A caller
    running the stages by hand therefore gets exactly what `canonicalize()` gets, which is the point.

    Call this with the journal clean -- outside any edit scope, or inside one with nothing
    pending -- because classification reads the arena.  A caller building a molecule mid-flight
    adds every atom, adds each aromatic bond, then calls this: whatever it has not stated is
    H_UNKNOWN, so the free choice is available where it belongs.  This function opens its own edit
    scope, so on return the orders are applied unless the caller holds an outer scope, in which case
    they apply when that scope closes -- and in that one case the heal above does NOT run, since the
    orders it would read are still pending.  Such a caller calls `derive_hydrogens(fill_only=True)`
    itself once its scope has closed.

    Returns a `KekuleResult`.  `.changed` is False when nothing moved -- a molecule with no
    aromatic bonds, or a second call -- so a caller can see idempotence rather than take it on
    trust.  `.log` is a list of human-readable lines, one per repair: a bond in no ring, a
    hypercondensed atom, a non-aromatic element, a charge shift, a derived hydrogen, a system with
    no Kekule form.  `.unresolved` is a list of tuples of stable ids, one per aromatic system that
    could not be fully assigned.

    AN UNRESOLVED SYSTEM IS STILL REWRITTEN, and a caller who expects it back as drawn will be
    surprised: it gets the best matching found, so its aromatic bonds are written as orders like
    every other system's and the atoms the matching could not pair come back with an incomplete
    valence and NO radical flag -- `check_valence()` is how to find them, and reports `'violation'`.
    That is deliberate on both counts.  A partial answer keeps the record readable and the rest of
    the molecule usable, which is what an input-is-garbage posture requires; and a radical flag would
    be a claim about what the input meant, which one failed matching is no evidence for.
    Nothing recoverable raises.
    `AromaticKekulizeError` means a broken invariant in this file, never a bad input.
    """
    cdef _AromRun run = arom_prepare(mol, aromatic_bonds, stated_h)
    cdef arom_scratch_t *sc = &run.sc
    cdef uint32_t n = run.n
    cdef uint32_t m = run.m
    cdef list numbers = mol._numbers
    cdef list unresolved = []
    if not m:
        return arom_result(mol, False, run.log, unresolved)

    cdef Structure structure = mol._structure
    cdef atom_t *atoms = structure.atoms()
    cdef uint32_t *ptr = csr_ptr(structure)
    cdef halfedge_t *edges = csr_edges(structure)
    cdef uint32_t i, j, k, v, w, count, head, left, nlow
    cdef uint64_t nodes
    cdef list names, cand, trial, trials, saved, hcand, relaxed_h, pool
    cdef list shifts = []                 # (ring id, substituent id, kind), applied with the emit
    cdef list hydrogens = []              # ring ids that take one derived H, applied with the emit
    cdef list relaxed = []                # (ring id, charge, h, substituent id or -1), ditto
    cdef dict seen
    cdef dict system_changed = {}     # component root -> at least one of its bonds moved
    cdef dict system_names = {}
    cdef set unresolved_roots = set()
    cdef object picked, c, key, why, root
    cdef bint ok
    for i in range(n):
        sc.comp[i] = -1
    for i in range(n):
        if sc.comp[i] >= 0 or not sc.adeg[i]:
            continue
        # one aromatic system, breadth first over live aromatic edges
        head = 0
        sc.clist[0] = i
        sc.comp[i] = <int32_t> i
        count = 1
        while head < count:
            v = sc.clist[head]
            head += 1
            for k in range(sc.aptr[v], sc.aptr[v + 1]):
                if not sc.e_alive[sc.aeid[k]]:
                    continue
                w = sc.aadj[k]
                if sc.comp[w] < 0:
                    sc.comp[w] = <int32_t> i
                    sc.clist[count] = w
                    count += 1
        nodes = 0
        names = None
        relaxed_h = None
        ok = arom_match(sc, sc.clist, count, &nodes)
        if not ok:
            # sorted, not in traversal order: a system's identity is its atoms, and reporting
            # them in the order a search happened to reach them would make the caller's log
            # depend on the search
            names = []
            for j in range(count):
                names.append(numbers[sc.clist[j]])
            names.sort()
        if nodes <= AROM_NODE_BUDGET:
            # A charge shift that makes a must-not N-oxide nitrogen take a ring double bond.  Two
            # gates, and which one applies is decided by the arithmetic in `arom_shift_candidates`
            # rather than by chemistry re-litigated here:
            #
            # * a BALANCED candidate set -- as many nitrogens mis-spelled cationic as anionic -- is
            #   charge-conserving taken as a whole, so it is tried whether or not the system already
            #   kekulised.  Pyrazine N,N'-dioxide arrives from real files as `[O-]n1cc[n+](=O)cc1`,
            #   with each half wrong in the opposite direction; taken together they cancel.  Without
            #   this gate that ring kekulises as 1,4-dihydropyrazine -- every atom valence-legal,
            #   and the aromaticity the file asserted on all six bonds thrown away.
            # * an UNBALANCED set moves the molecule's total charge, so it is offered only to a
            #   system that has no Kekule form as written and would otherwise be handed back
            #   half-assigned.
            #
            # THE MOLECULE IS NOT TOUCHED TO FIND OUT.  A trial reclassifies in the scratch and
            # re-runs the match, so a shift that does not help costs one search and leaves no trace,
            # and only a trial that produced a complete matching is ever written back.  Mutating
            # first and re-running would be simpler to write and would corrupt the charges of every
            # molecule the retry failed on -- which is most of them, since a system reaching the
            # unbalanced gate is usually just broken.
            #
            # For the unbalanced gate: all candidates together first, then each alone.  Those three
            # subsets are all of them for the counts that occur; a system with three candidates gets
            # four of its seven subsets tried and then falls back to the partial assignment, which
            # is the answer it would have got with no retry at all.
            cand = arom_shift_candidates(sc, atoms, ptr, edges, sc.clist, count)
            trials = []
            if cand:
                nlow = 0
                for c in cand:
                    if c[2]:
                        nlow += 1
                if nlow and 2 * nlow == <uint32_t> len(cand):
                    trials.append(cand)
                if not ok:
                    if not trials:
                        trials.append(cand)
                    if len(cand) > 1:
                        for c in cand:
                            trials.append([c])
            # the matching a balanced trial is about to overwrite, kept only when there is one to
            # lose -- an unbalanced trial runs on a system that failed, so there is nothing.  A
            # snapshot and not a second search: restoring by re-running `arom_match` would work only
            # while it is deterministic, and would need an error path for the case where it somehow
            # was not, which is an error path no test could ever reach.
            saved = None
            if ok and trials:
                saved = []
                for j in range(count):
                    saved.append(sc.mate[sc.clist[j]])
            picked = None
            for trial in trials:
                for c in trial:
                    sc.cls[<uint32_t> c[0]] = AROM_MUST
                for j in range(count):
                    sc.mate[sc.clist[j]] = -1
                nodes = 0
                if arom_match(sc, sc.clist, count, &nodes):
                    picked = trial
                    break
                for c in trial:
                    sc.cls[<uint32_t> c[0]] = AROM_MUST_NOT
            if picked is None:
                if saved is not None:
                    for j in range(count):
                        sc.mate[sc.clist[j]] = <int32_t> saved[j]
            else:
                ok = True
                if names is None:
                    names = []
                    for j in range(count):
                        names.append(numbers[sc.clist[j]])
                    names.sort()
                # which gate accepted this trial is a property of the trial itself: a balanced set
                # is the one that conserves charge, and the log has to say which of the two things
                # it found rather than assert the stronger claim in both cases
                nlow = 0
                for c in picked:
                    if c[2]:
                        nlow += 1
                why = ('is mis-spelled in both directions at once'
                       if nlow and 2 * nlow == <uint32_t> len(picked)
                       else 'has no Kekule form as written')
                for c in picked:
                    v = <uint32_t> c[0]
                    w = <uint32_t> c[1]
                    shifts.append((numbers[v], numbers[w], c[2]))
                    if c[2]:
                        run.log.append(mc_record(
                            'kekule:oxide-anion',
                            (numbers[v], numbers[w]),
                            f'aromatic system {tuple(names)!r} {why}; '
                            f'{symbol_of(&atoms[w])} {numbers[w]} read as '
                            f'{numbers[w]}(-) single bonded to nitrogen {numbers[v]}',
                            mc_repaired()))
                    else:
                        run.log.append(mc_record(
                            'kekule:oxide-cation',
                            (numbers[v], numbers[w]),
                            f'aromatic system {tuple(names)!r} {why}; nitrogen {numbers[v]} read '
                            f'as {numbers[v]}(+), the cation its '
                            f'{symbol_of(&atoms[w])}(-) substituent {numbers[w]} needs',
                            mc_repaired()))
            # THE HYDROGEN RELAXATION, tried last and only on a system still without a Kekule form:
            # an aromatic nitrogen whose hydrogen count the input never stated.  `c1cncn1` and
            # `c1ccnc1` arrive that way constantly -- see `arom_h_candidates` for why the notation
            # loses the count and why deriving it here is the only reading that accepts them.
            #
            # ONE TRIAL AND NO SUBSET SEARCH, unlike the charge shifts above.  Every candidate goes
            # to AROM_MAY at once, which does not commit any of them: `arom_match` uses a may-match
            # atom as a partner where a ring forces it and leaves it alone where it does not, so a
            # single search finds both the matching and the smallest set of nitrogens that has to
            # carry a hydrogen.  Pyrimidine needs none of them and pyrazole needs one, and neither
            # answer is a subset the caller had to enumerate.
            #
            # WHICH nitrogens is not known yet -- `arom_extend` below still gets to pair leftover
            # may-match atoms with each other -- so this records the candidates and the count is
            # read off the matching afterwards.
            if not ok and nodes <= AROM_NODE_BUDGET:
                hcand = arom_h_candidates(sc, atoms, sc.clist, count)
                if hcand:
                    for c in hcand:
                        sc.cls[<uint32_t> c] = AROM_MAY
                    for j in range(count):
                        sc.mate[sc.clist[j]] = -1
                    nodes = 0
                    if arom_match(sc, sc.clist, count, &nodes):
                        ok = True
                        relaxed_h = hcand
                    else:
                        # no trace left: the classes go back and the caller gets the partial
                        # assignment it would have got with no retry at all
                        for c in hcand:
                            sc.cls[<uint32_t> c] = AROM_MUST
            # THE STATE RELAXATIONS, tried last of all and only on a system that survived every
            # other offer without a Kekule form: a must-not atom that one added positive charge, or
            # one withdrawn hydrogen, would turn into a must-match one.  A pyridinium written
            # neutral and an azole carrying a hydrogen it cannot afford -- the two shapes
            # `CORPUS_SHAPES` holds.  See the two generator docstrings.
            #
            # PINNED TO AROM_MUST AND NOT RELAXED TO AROM_MAY, unlike the hydrogen derivation above.
            # A may-match atom lets the search decide whether to use it, which is right when the
            # question is "does this nitrogen want its pi electron"; here the question is "is this
            # atom in a valence state at all", and the answer the table gave is must-match.  An atom
            # left may-match would be written a charge it then did not need.
            #
            # SINGLES FIRST, ORDERED BY WHAT THEY COST, then the whole pool.  A charge-conserving
            # candidate (a cation whose substituent has a proton to give up, or a surplus hydrogen)
            # is tried before one that moves the molecule's total charge, so `c1ccn(O)cc1` lands on
            # pyridine N-oxide rather than on an N-hydroxypyridinium cation.  A single candidate is
            # what frees every shape the net holds; the whole-pool trials are there for one where no
            # single candidate does, and the worst pool is five.  NO ATTEMPT
            # COUNTER and no subset search: a bounded list of trials is what keeps the answer for
            # one ring independent of how many other rings the record happens to carry.
            if not ok and nodes <= AROM_NODE_BUDGET:
                seen = {}
                for c in arom_surplus_h_candidates(sc, atoms, sc.clist, count):
                    seen[c[0]] = c
                for c in arom_cation_candidates(sc, atoms, ptr, edges, sc.clist, count):
                    if c[0] not in seen:
                        seen[c[0]] = c
                # EXPLICIT LOOPS, NOT COMPREHENSIONS, and it is not a style preference.  A
                # comprehension has its own scope in Cython 3, so its loop variable is a fresh
                # implicit local whatever the enclosing `cdef` block declares -- three build warnings
                # for the one line this replaced, and RULES.md 9.7 is exactly about not leaving an
                # inferred Python local where a declaration was meant.
                pool = []
                for key in sorted(seen):
                    pool.append(seen[key])
                # kind first, plain cation last: a relaxation that conserves charge (kind 1) or takes
                # away a hydrogen the ring could not afford (kind 2) is a smaller claim about the
                # input than inventing a cation (kind 0), so it is tried before one.
                trials = []
                for c in pool:
                    if c[4]:
                        trials.append([c])
                for c in pool:
                    if not c[4]:
                        trials.append([c])
                if len(pool) > 1:
                    trials.append(pool)
                picked = None
                for trial in trials:
                    for c in trial:
                        sc.cls[<uint32_t> c[0]] = AROM_MUST
                    for j in range(count):
                        sc.mate[sc.clist[j]] = -1
                    nodes = 0
                    if arom_match(sc, sc.clist, count, &nodes):
                        picked = trial
                        break
                    for c in trial:
                        sc.cls[<uint32_t> c[0]] = AROM_MUST_NOT
                if picked is not None:
                    ok = True
                    for c in picked:
                        v = <uint32_t> c[0]
                        w = <uint32_t> c[3]
                        relaxed.append((numbers[v], c[1], c[2],
                                        -1 if c[3] == AROM_NO_ATOM else numbers[w]))
                        if c[4] == 2:
                            run.log.append(mc_record(
                                'kekule:surplus-h',
                                (numbers[v],),
                                f'aromatic system {tuple(names)!r} has no Kekule form with the '
                                f'hydrogen counts as written; '
                                f'{symbol_of(&atoms[v])} {numbers[v]} read as the '
                                f'one-electron contributor it has to be, its surplus hydrogen '
                                f'dropped',
                                mc_repaired()))
                        elif c[4]:
                            run.log.append(mc_record(
                                'kekule:cation-oxide',
                                (numbers[v], numbers[w]),
                                f'aromatic system {tuple(names)!r} has no Kekule form as written; '
                                f'{symbol_of(&atoms[v])} {numbers[v]} read as '
                                f'{numbers[v]}(+) and its {symbol_of(&atoms[w])} '
                                f'substituent {numbers[w]} as {numbers[w]}(-), the charge-separated '
                                f'oxide the ring needs',
                                mc_repaired()))
                        else:
                            run.log.append(mc_record(
                                'kekule:missing-charge',
                                (numbers[v],),
                                f'aromatic system {tuple(names)!r} has no Kekule form as written; '
                                f'{symbol_of(&atoms[v])} {numbers[v]} read as '
                                f'{numbers[v]}(+), the cation its {sc.nbrs[v]} substituents imply',
                                mc_repaired()))
        if ok:
            arom_extend(sc, sc.clist, count)
            if relaxed_h is not None:
                for c in relaxed_h:
                    v = <uint32_t> c
                    if sc.mate[v] >= 0:
                        continue              # the ring wanted its pi electron after all
                    hydrogens.append(numbers[v])
                    run.log.append(mc_record(
                        'kekule:derived-h',
                        (numbers[v],),
                        f'aromatic system {tuple(names)!r} has no Kekule form with the hydrogen '
                        f'counts as written; {symbol_of(&atoms[v])} {numbers[v]} read as '
                        f'the two-electron donor it has to be, with one hydrogen',
                        mc_repaired()))
        else:
            left = arom_partial(sc, sc.clist, count)
            unresolved.append(tuple(names))
            unresolved_roots.add(i)          # `i` is the component root: `sc.comp[i] == i` by construction
            # THE MESSAGE SAYS THE SYSTEM WAS REWRITTEN, because it was.  "left unsaturated" on its
            # own read as though the system came back as drawn, and it does not: the fallback writes
            # the best matching it found, so the aromatic flags are gone either way and the atoms the
            # matching could not pair carry single bonds and an incomplete valence.  A caller reading
            # only the first half of that sentence goes looking for aromatic bonds that are no longer
            # there.  The deficient atoms are NOT given a radical flag -- an unresolved system is not
            # evidence that the input meant a radical -- so `check_valence()` reports them, and that
            # is the intended way to find them.
            if nodes > AROM_NODE_BUDGET:
                run.log.append(mc_record('kekule:budget-exhausted', tuple(names),
                                         f'aromatic system {tuple(names)!r} exhausted the kekulisation search '
                                         f'budget; written as the best matching found, which leaves {left} '
                                         f'atom(s) with an incomplete valence and no radical flag',
                                         mc_lost()))
            else:
                run.log.append(mc_record('kekule:no-kekule-form', tuple(names),
                                         f'aromatic system {tuple(names)!r} has no Kekule form; written as the '
                                         f'best matching there is, which leaves {left} atom(s) with an '
                                         f'incomplete valence and no radical flag',
                                         mc_lost()))

    # --- emit.  Every stated aromatic bond is written, the matched ones double and the rest
    # single, so calling this twice on the same molecule is the same as calling it once.
    cdef list doubles = []
    cdef list singles = []
    cdef object pair
    cdef halfedge_t *e
    cdef bint changed = bool(shifts) or bool(hydrogens) or bool(relaxed)
    cdef uint8_t target
    for k in range(m):
        v = sc.e_u[k]
        w = sc.e_v[k]
        if sc.e_alive[k] and sc.mate[v] == <int32_t> w:
            if sc.mate[w] != <int32_t> v:
                raise AromaticKekulizeError(
                    f'matching is not symmetric at {numbers[v]}-{numbers[w]}')
            doubles.append((numbers[v], numbers[w]))
            target = 2
        else:
            singles.append((numbers[v], numbers[w]))
            target = 1
        # `changed` is measured against what the arena holds rather than assumed from having
        # written something: writing order 1 over a bond that was already order 1 changes nothing,
        # which is what makes a second call observably a no-op
        e = csr_find(structure, v, w)
        if e is not NULL and e.order != target:
            changed = True
            if sc.comp[v] >= 0:
                system_changed[sc.comp[v]] = True
    # THE KEKULISATION ITSELF IS AN EVENT, and `INFO` is its severity -- "did what was asked", whose own
    # docstring names `'kekulized'`.  Every other record here is a repair of something the input stated,
    # so without this one a molecule whose every bond this function rewrote came back `changed=True`
    # with an empty log, and `mol.canonicalize()` left no trace of having rewritten a ring twice.
    #
    # ONE PER SYSTEM THAT MOVED, measured against the arena rather than assumed from having run: a
    # second call over the same edge set writes the orders that are already there, and a record then
    # would say work happened where none did.  An unresolved system is excluded because its own LOST
    # record already says what was written -- the best matching there is, which is not a Kekule form.
    for i in range(n):
        if sc.comp[i] >= 0 and sc.comp[i] in system_changed and sc.comp[i] not in unresolved_roots:
            if sc.comp[i] in system_names:
                system_names[sc.comp[i]].append(numbers[i])
            else:
                system_names[sc.comp[i]] = [numbers[i]]
    for root in sorted(system_names):
        names = system_names[root]
        names.sort()
        run.log.append(mc_record('kekule:kekulized', tuple(names),
                                 f'aromatic system {tuple(names)!r} written as a Kekule form'))

    # REPRESENTATION, NOT STRUCTURE, and `_apply` needs to be told so.  These orders reach the journal
    # as ordinary OP_SET_ORDER records, indistinguishable from a caller's `set_order`, and `_apply`
    # drops stored CIP descriptors on the latter.  It must not drop them here: the aromatic form and
    # the Kekule form are the same molecule differently written, and a stored descriptor is an
    # assertion the INPUT made about the molecule rather than a function of the orders being stored --
    # so nothing this function does can make it false.  See `MoleculeContainer._representation_change`.
    mol._representation_change = True
    try:
        with mol.edit():
            for pair in singles:
                mol.set_order(pair[0], pair[1], 1)
            for pair in doubles:
                mol.set_order(pair[0], pair[1], 2)
            # the accepted charge shifts, in the same scope as the orders they made possible: a
            # molecule that briefly held the new orders with the old charges would be one no reader
            # of the arena should ever be able to observe
            for pair in shifts:
                if pair[2]:
                    mol.set_charge(pair[1], -1)
                    mol.set_order(pair[0], pair[1], 1)
                else:
                    mol.set_charge(pair[0], 1)
            # and the derived hydrogen counts, in the same scope for the same reason: the count and
            # the single bonds that made it necessary are one reading of the ring, and a molecule
            # briefly holding one without the other is a hypovalent nitrogen no reader should see
            for c in hydrogens:
                mol.set_hydrogens(<uint32_t> c, 1)
            # and the accepted state relaxations, in the same scope and for the same reason.  These
            # are absolute values rather than deltas because the arena may not be read from inside a
            # session: what to write was decided out there, where reading was still allowed.
            for pair in relaxed:
                mol.set_charge(pair[0], pair[1])
                mol.set_hydrogens(pair[0], pair[2])
                if pair[3] >= 0:
                    mol.set_charge(pair[3], -1)
                    mol.set_hydrogens(pair[3], 0)
    finally:
        mol._representation_change = False

    # --- THE HYDROGEN HEAL.  `kekule()` heals implicit hydrogens itself where there are no valence
    # errors, and it does so HERE rather than in `canonicalize()` between step 1 and step 2: a pipeline
    # is never stronger than its parts, and healing from the pipeline would make `mol.kekule()` by hand
    # a WEAKER operation than the same call inside it -- an asymmetry a caller cannot see or guess.
    #
    # WHAT IT HEALS is exactly what this function just made answerable.  A format with no hydrogen
    # channel -- an MDL bond block of type 4, any drawn format -- cannot say whether an aromatic
    # pnictogen is a pyrrole or a pyridine, so the read-time derivation stores H_UNKNOWN there rather
    # than guessing (`HYD_AMBIGUOUS_AROMATIC`, and `arom_prepare` reads that nibble back as the free
    # choice it is).  The emit above has just replaced those aromatic bonds with definite orders, so
    # the ordinary valence rows answer now.  Nothing else changes: `fill_only` writes only atoms that
    # still claim nothing, which is what keeps it off the counts a reader derived from something no
    # valence row reproduces -- diborane's bridging hydrogens and ferrocene were both measured losing
    # a correct 0 to a blanket recompute.
    #
    # "IF NO VALENCE ERRORS" IS PER SYSTEM, not per molecule, and `unresolved` is the whole test.  An
    # atom the matching could not pair carries a deficient order sum, so a valence row asked about it
    # would answer with the hydrogens that fill the deficit -- inventing hydrogens the input never had
    # and hiding the very defect this function reports.  So those atoms are named and skipped, they
    # stay H_UNKNOWN, and `check_valence()` still finds them, which is what the docstring promises.  A
    # ring that failed does not cost a ring that succeeded its counts.
    #
    # NOT LOGGED, and that is a decision rather than an omission.  A log record here would fire for
    # every pyridine in every SDF -- and it would be describing the read-time derivation finishing, not
    # a repair of anything the file stated.  The counts this function REWRITES against a statement are
    # logged, above, one line each.
    #
    # ONLY WHEN THE ARENA IS READABLE.  A caller holding an outer edit scope has our orders pending in
    # its journal, so there is nothing to derive from yet; deriving is then that caller's business
    # after its scope closes, and this function says so in its docstring.
    cdef list deficient
    if not mol._journal_len:
        deficient = None
        if unresolved:
            deficient = []
            for pair in unresolved:
                deficient.extend(pair)
        derive_implicit_hydrogens(mol, deficient, True)
    return arom_result(mol, changed, run.log, unresolved)


def kekule_classify(MoleculeContainer mol not None, aromatic_bonds=None, stated_h=None):
    """`{n: 'must' | 'may' | 'must_not'}` for the atoms of the stated aromatic system.

    The classification table is the chemistry in this file, so it is testable on its own rather
    than only through the bond orders it produces -- a wrong class and a wrong search both come
    out as wrong orders, and telling them apart afterwards is guesswork.

    Reported for the molecule AS `kekule` WOULD CLASSIFY IT, which since the charge separation
    became part of that means on a copy: `O=n1ccccc1` answers `must` for its nitrogen, because that
    is the atom the search will see.  Classifying the argument as it stands would make this function
    disagree with the one it exists to explain.  What it still does not show is the last-resort
    charge shift, and neither is the derived hydrogen: those are not classes but relaxations the
    search reaches for on a system that has already failed.  `[O-]n1ccccc1` answers `must_not` here
    and kekulises anyway, and `c1cncn1` answers `must` for both nitrogens and kekulises anyway.
    """
    cdef MoleculeContainer work = mol.copy()
    cdef _AromRun run = arom_prepare(work, aromatic_bonds, stated_h)
    cdef list numbers = work._numbers
    cdef uint32_t i
    cdef dict out = {}
    for i in range(run.n):
        if run.sc.adeg[i]:
            out[numbers[i]] = AROM_CLASS_NAMES[run.sc.cls[i]]
    return out


def kekule_copy(MoleculeContainer mol not None, aromatic_bonds=None, stated_h=None):
    """A Kekule copy of `mol`, with `mol` itself untouched.

    The adapter the InChI bridge registers, and the shape every caller that must not mutate its
    input wants: `kekule()` proper works in place, because it is an operation the owner of a
    molecule performs on it deliberately, while a format bridge is holding somebody else's
    molecule and input fidelity is an invariant of `molecule_to_inchi`.

    Raises `AromaticKekulizeError` when some system has no Kekule form.  A partial answer is the
    right thing to hand a caller who can inspect the log; it is the wrong thing to hand libinchi,
    which would silently receive a molecule with the wrong bond orders and return a wrong InChI.

    The two optional arguments mean what they mean on `kekule`.  The InChI hook passes neither --
    it takes the stored aromatic bonds, which is the only thing it could know about.  They are in
    the signature because without them the raising path is unreachable while the arena's order-4
    gate is shut, and an untestable error path is one that is wrong when it finally runs.
    """
    cdef MoleculeContainer out = mol.copy()
    cdef KekuleResult result = kekule(out, aromatic_bonds, stated_h)
    if result.unresolved:
        raise AromaticKekulizeError(
            f'{len(result.unresolved)} aromatic system(s) have no Kekule form, so there are no '
            f'bond orders to hand a consumer that requires them: '
            f'{"; ".join(map(str, result.log))}')
    return out


# The InChI bridge needs Kekule orders and cannot perceive them itself.  `_ich_set_kekule_fn`
# exists so neither file imports the other: `_inchi.pxi` declares the hook, this file fills it at
# module init, and the direction of the dependency is the one that makes sense -- the kekuliser
# knows nothing about InChI.  This is the one registration; if a second consumer ever needs the
# same adapter it calls `kekule_copy` directly rather than growing another hook.
_ich_set_kekule_fn(kekule_copy)
