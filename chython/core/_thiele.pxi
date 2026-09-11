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
# Aromatisation: Kekule bond orders become order 4.  The inverse of `_kekule.pxi`.
#
# THE RULE, AND IT IS ONE SENTENCE
#
# An aromatic edge set is admissible when the molecule the caller already has is a Kekule form of
# it.  Not "looks aromatic": `kekule()` would have produced these very orders.  So this file calls
# `arom_classify_atom` -- the read epic's one-atom table, which says whether an atom takes exactly
# one ring double bond (MUST), none (MUST_NOT) or either (MAY) -- and checks the orders in hand
# against that assignment.  Two independent copies of the chemistry would disagree on some charged
# heteroatom and nobody would find out for weeks, which is the same argument `_kekule.pxi` makes for
# not writing the table inside a parser.
#
# WHY THE MATCHING CONDITION IS NOT ENOUGH
#
# A must/may/must-not assignment that saturates every must atom is NOT the sextet condition, and two
# molecules say so:
#
#   * p-benzoquinone.  Four MUST carbons, two ring double bonds, a perfect matching -- and 4 pi
#     electrons in the ring.
#   * cyclobutadiene.  Four MUST carbons, two doubles, matched, 4 pi electrons, antiaromatic.
#
# So the matching is necessary and not sufficient, and three further conditions are needed.  Each is
# carried by a fixture in test_thiele.py rather than asserted here:
#
#   1. RING SIZE 5..7.  Size 4 is what refuses cyclobutadiene and >= 8 refuses cyclooctatetraene,
#      both of which the matching admits.
#   2. NO EXOCYCLIC DOUBLE OR TRIPLE BOND on a ring atom.  An atom that spent its pi electron
#      outside the ring is not donating it to the ring; this is what refuses p-benzoquinone, fulvene
#      and both pyridones.
#   3. HUCKEL, BUT ONLY ON AN ISOLATED CYCLE.  pi = one per MUST + two per lone-pair donor + zero per
#      empty-orbital atom, and the count must be 2 mod 4.  PER RING IT IS WRONG -- azulene's
#      five-ring has five MUST atoms and no donor and azulene is aromatic -- because a fused system's
#      matching spans its rings.  On a component of the aromatic edge set that is a single cycle it is
#      exactly right, and it is the only thing that refuses 1H-azepine (six MUST + an N donor = 8 pi)
#      and the cyclopentadienyl cation (four MUST + an empty C+ = 4 pi).
#
# Condition 3 is stated in general and not as a special case: allowing a seven-ring whose one non-sp2
# atom is BORON is the same parity rule with the general case filed off, neutral boron being the
# empty-orbital atom.  Stating it generally is what makes borole come out NOT aromatic here, one of
# the four measured divergences from chython 2; spec section 14.2 has the table and the reasons.
#
# WHAT THIS FILE DOES NOT DO
#
# It does not shift a pyrrole hydrogen.  Moving a hydrogen is a tautomer decision wearing a
# representation change's clothes, and it belongs to the standardization pack.  `thiele()` changes bond
# orders and nothing else -- the atoms it cannot spell as aromatic come back in `.refused` instead.
#
# It carries no SMARTS list re-aromatising rings that a per-ring sp2 COUNT dropped, because there is no
# such count: `arom_thiele_ring_ok` below has no sp2 filter, so `N1C=CN2C=CC=C12` and the five other
# bicycles in `test/heterocycles_charges.smi` come out fully aromatic from the one pass.  Do not add
# the count.


# a candidate ring atom's pi contribution, for the Huckel count of an isolated cycle
cdef enum:
    THIELE_PI_EMPTY = 0        # an empty p orbital: C+, neutral B
    THIELE_PI_ONE = 1          # one electron, from a ring double bond
    THIELE_PI_PAIR = 2         # a lone pair: pyrrole N, furan O, thiophene S, C-, B-

# a class value outside AROM_MUST / AROM_MAY / AROM_MUST_NOT, for an atom `arom_classify_atom`
# reported as not a valid aromatic state.  A sentinel and not "read it as MUST", because a spurious
# MUST atom that happens to carry one double bond would PASS the matching check -- the invalid atom
# has to fail the component outright.
cdef enum:
    THIELE_INVALID = 3

# WHY THE PRE-FILTER DROPPED A RING, for the two of its five exits that owe the caller a reason.
# `arom_thiele_ring_ok` is `noexcept nogil` and cannot append to a Python list, so it names the state
# and the atom through out-parameters and `thiele` does the reporting.
#
# Only two of the exits are refusals in `.refused`'s sense -- the ring LOOKS like a Kekule aromatic
# and is declined for a state that is fixable.  The other three say "not a candidate at all": a
# non-aromatizable element, a bond that is already aromatic, and a spiro or over-coordinated atom.
# A line per non-candidate ring would put four of them in every steroid, so those stay silent.
cdef enum:
    THIELE_WHY_NONE = 0
    THIELE_WHY_RADICAL = 1
    THIELE_WHY_H_UNKNOWN = 2


cdef struct thiele_t:
    void *block
    uint8_t *ring_ok           # [nrings]  1 while the ring is still a candidate
    uint8_t *he_arom           # [nhalf]   1 when this half-edge is in the candidate edge set
    uint8_t *cand              # [n]       1 when the atom is in a candidate ring
    uint8_t *cls               # [n]       AROM_MUST / AROM_MAY / AROM_MUST_NOT / THIELE_INVALID
    int32_t *comp              # [n]       component label over the candidate edge set, -1 outside
    uint32_t *clist            # [n]       one component's atoms
    uint32_t *stack            # [n]       breadth-first queue


cdef int arom_thiele_alloc(thiele_t *t, uint32_t n, uint32_t nhalf, uint32_t nrings) except -1:
    # one struct, one malloc, one check, one free -- RULES.md 5.2
    cdef size_t n_u8 = align8(<size_t> n * sizeof(uint8_t))
    cdef size_t n_u32 = align8(<size_t> n * sizeof(uint32_t))
    cdef size_t n_i32 = align8(<size_t> n * sizeof(int32_t))
    # `or 1`: a molecule with no rings and a molecule with no bonds both reach here, and a zero-byte
    # PyMem_Malloc may return NULL, which this function would report as MemoryError
    cdef size_t half_u8 = align8(<size_t> (nhalf if nhalf else 1) * sizeof(uint8_t))
    cdef size_t ring_u8 = align8(<size_t> (nrings if nrings else 1) * sizeof(uint8_t))
    cdef size_t total = ring_u8 + half_u8 + 2 * n_u8 + n_i32 + 2 * n_u32
    cdef char *block = <char *> PyMem_Malloc(total if total else 1)
    if block is NULL:
        raise MemoryError('aromatisation scratch allocation failed')
    memset(block, 0, total)
    t.block = <void *> block
    cdef size_t off = 0
    t.ring_ok = <uint8_t *> (block + off); off += ring_u8
    t.he_arom = <uint8_t *> (block + off); off += half_u8
    t.cand = <uint8_t *> (block + off); off += n_u8
    t.cls = <uint8_t *> (block + off); off += n_u8
    t.comp = <int32_t *> (block + off); off += n_i32
    t.clist = <uint32_t *> (block + off); off += n_u32
    t.stack = <uint32_t *> (block + off); off += n_u32
    return 0


cdef class _ThieleRun:
    """One aromatisation in flight: the scratch and the log.

    A cdef class for the reason `_AromRun` is one -- `__dealloc__` is the only place the scratch is
    freed on every exit path, including an exception raised while a ring is being classified.
    """
    cdef thiele_t t
    cdef uint32_t n
    cdef uint32_t nrings
    cdef list log

    def __cinit__(self):
        self.t.block = NULL
        self.n = 0
        self.nrings = 0
        self.log = []

    def __dealloc__(self):
        if self.t.block is not NULL:
            PyMem_Free(self.t.block)
            self.t.block = NULL


cdef class ThieleResult:
    """What `thiele()` did: `changed`, `log`, `refused`.

    The same three questions `KekuleResult` answers and a distinct type, because the two operations
    fail differently: kekulisation is handed an edge set and may find no Kekule form for it, while
    aromatisation CHOOSES the edge set and its interesting answer is which candidate systems it
    declined.  A shared type would make one of the two docstrings a lie.
    """
    cdef readonly bint changed
    """False when no bond order moved: nothing aromatisable, or a second call."""
    cdef readonly list log
    """One human-readable line per candidate system declined, empty when none was.

    A COPY, exactly as `KekuleResult.log` is: the storage is `mol.log`."""
    cdef readonly list refused
    """A tuple of stable ids per candidate ring system that was declined; empty when none was."""

    def __repr__(self):
        return (f'ThieleResult(changed={bool(self.changed)}, log={self.log!r}, '
                f'refused={self.refused!r})')


cdef ThieleResult arom_thiele_result(MoleculeContainer mol, bint changed, list log, list refused):
    cdef ThieleResult r = ThieleResult.__new__(ThieleResult)
    r.changed = changed
    r.log = log
    r.refused = refused
    # Unconditional, for the reason `arom_result` states.
    mol.log.absorb('thiele', log, rule='thiele')
    return r


cdef inline uint8_t arom_thiele_pi(atom_t *a, bint carries_double) noexcept nogil:
    """The atom's pi contribution to an isolated cycle's Huckel count.

    Keyed on whether the atom ACTUALLY carries a candidate double bond, not on its class, because a
    MAY atom is exactly the one whose class does not say: a pyridinium `[nH+]` in a six-ring carries
    one and brings one electron, a pyrrolium `[nH2+]` in a five-ring carries none and brings zero,
    and `arom_classify_atom` answers MAY for both.

    Without a double bond the atom either holds a lone pair in the p orbital or holds nothing, and
    the sign of the charge decides every case that can occur: an anion donates, a cation cannot,
    neutral boron has an empty orbital, and every other neutral atom that classifies non-MUST is a
    pyrrole-type donor (N with three neighbours or an NH, furan O, thiophene S, Se, Te, phosphole P).
    A NEUTRAL CARBON CANNOT REACH THE ELSE BRANCH -- `arom_classify_atom` answers MUST or `invalid`
    for every neutral carbon -- which is why there is no carbon case here, and why a radical is
    refused by the pre-filter before this is ever asked.
    """
    if carries_double:
        return THIELE_PI_ONE
    if a.charge > 0:
        return THIELE_PI_EMPTY
    if a.charge < 0:
        return THIELE_PI_PAIR
    if a.element == AROM_Z_B:
        return THIELE_PI_EMPTY
    return THIELE_PI_PAIR


cdef inline void arom_thiele_bonds(Structure structure, uint32_t i, uint32_t *nbrs,
                                   uint32_t *doubles, uint32_t *triples,
                                   uint32_t *aromatics) noexcept nogil:
    """One atom's bond census: neighbours, doubles, triples, order-4 bonds.

    `nbrs` is `arom_classify_atom`'s argument: every bond of any order except 8, explicit hydrogens
    included.  Order 4 is counted separately and NOT as a double,
    because a molecule that already holds aromatic bonds is not a Kekule form of anything and the
    ring carrying them is skipped rather than read as unsaturated.
    """
    cdef uint32_t *ptr = csr_ptr(structure)
    cdef halfedge_t *edges = csr_edges(structure)
    cdef uint32_t k
    cdef uint8_t order
    nbrs[0] = 0
    doubles[0] = 0
    triples[0] = 0
    aromatics[0] = 0
    for k in range(ptr[i], ptr[i + 1]):
        order = edges[k].order
        if order == 8:
            continue                      # a dative bond is not a neighbour for this census
        nbrs[0] += 1
        if order == 2:
            doubles[0] += 1
        elif order == 3:
            triples[0] += 1
        elif order == 4:
            aromatics[0] += 1


cdef bint arom_thiele_ring_ok(Structure structure, uint32_t *ring, uint32_t size,
                              uint8_t *cls_out, uint8_t *why, uint32_t *culprit) noexcept nogil:
    """The per-ring pre-filter: is this ring worth putting into the candidate edge set at all?

    On a `False` return `why[0]` is a `THIELE_WHY_*` code, and `THIELE_WHY_NONE` for the three exits
    that mean "not a Kekule aromatic candidate at all" rather than a refusal; `culprit[0]` is the atom
    index behind a code.  See the enum for which is which and why the split is not cosmetic.

    THIS IS THE STEP THAT KEEPS TETRALIN AROMATIC.  Without it the saturated ring of a tetralin is a
    candidate too, its CH2 carbons classify MUST with no double bond to match, and the failure takes
    the FUSED benzene down with it -- the matching check runs over a whole component, so one bad ring
    poisons its neighbours.  The filter reads the BOND ORDERS directly rather than a count of
    hybridisation-2 atoms: hybridisation is a maintained cache, and the case of a bond may not rest on
    something that can go stale.

    Accepted: every atom carries exactly one double bond and no triple, or the ones that do not each
    classify MUST_NOT or MAY -- the pyrrole / furan / thiophene donor, or a charged carbon.  What
    refuses 1,4-dihydropyridine and tetralin is that gate and not a count: their sp3 CH2 is a neutral
    carbon with two hydrogens, which `arom_classify_atom` reports as no valid aromatic state at all.

    THERE IS DELIBERATELY NO CAP ON HOW MANY SUCH DONORS A RING MAY HOLD.  A cap of one -- justified by
    1,4-dihydropyridine, which the classify gate already refuses -- makes `N1C=CN2C=CC=C12` come out
    with five of its nine bonds aromatic.  That molecule is pyrrolo[1,2-a]imidazole: a 5-5 bicyclic
    whose EIGHT atoms share ten pi electrons -- three ring double bonds and a lone pair from each
    nitrogen -- so it is aromatic throughout, exactly as indolizine is.  Its imidazole ring holds two
    non-sp2 atoms because BOTH nitrogens are donors, and the bridging one spends its pair on the system
    rather than on either ring; per-ring bookkeeping cannot see that and reads it as
    1,4-dihydropyridine.

    Over-donation is still refused, one layer down: three donors in an isolated five-ring is eight pi
    electrons and `arom_thiele_check`'s Huckel test throws it out.  A FUSED component gets no Huckel
    test, and there the matching is the whole of the guarantee.

    WHICH LEAVES ONE THING THE CAP WAS DOING BY ACCIDENT, and it has to be said on purpose: a ring
    every one of whose atoms is a donor has NO double bond, and the matching check passes it
    vacuously -- MUST_NOT wants no double bond and finds none.  Borazine is that ring.  Its three
    nitrogens donate a pair each and its three borons contribute an empty orbital, which is six pi
    over six atoms and passes Huckel, so nothing downstream refuses it and `B1NBNBN1` came out as
    six aromatic bonds.  It is not a Kekule form of an aromatic six-ring -- an aromatic six-ring's
    Kekule form has three double bonds -- and the requirement below says so directly.  Note that
    this is the weaker per-RING statement and not per-component: a ring whose every double bond
    belongs to a fused neighbour is refused here even though the component has double bonds to
    spare.  No such ring appears in any fixture, and refusing an exotic one is the conservative
    direction; ring A of `N1C=CN2C=CC=C12` is not one of them, since its `C=C` is its own.
    """
    cdef uint32_t i, idx, k, nxt
    cdef uint32_t *ptr = csr_ptr(structure)
    cdef halfedge_t *edges = csr_edges(structure)
    cdef bint ring_double = False
    # initialised at the declaration although `arom_thiele_bonds` fills all four: it fills them
    # through pointers, which Cython cannot see, and the module compiles with warn.undeclared
    cdef uint32_t nbrs = 0, doubles = 0, triples = 0, aromatics = 0
    cdef atom_t *atoms = structure.atoms()
    cdef atom_t *a
    cdef uint8_t invalid
    why[0] = THIELE_WHY_NONE
    if size < 5 or size > 7:
        return False
    for i in range(size):
        idx = ring[i]
        a = atoms + idx
        if a.element != AROM_Z_B and a.element != AROM_Z_C and a.element != AROM_Z_N \
                and a.element != AROM_Z_O and a.element != AROM_Z_P and a.element != AROM_Z_S \
                and a.element != AROM_Z_AS and a.element != AROM_Z_SE and a.element != AROM_Z_TE:
            return False
        if at_radical(a):
            # a radical brings one electron and not two, so `arom_thiele_pi`'s charge rule does not
            # describe it; refused rather than guessed
            why[0] = THIELE_WHY_RADICAL
            culprit[0] = idx
            return False
        if at_implicit_h_unknown(a):
            # the pyrrole / pyridine choice IS the hydrogen count, so an unstated one cannot be
            # resolved into a stated aromatic form.  Note 4: no token on a prediction.
            why[0] = THIELE_WHY_H_UNKNOWN
            culprit[0] = idx
            return False
        arom_thiele_bonds(structure, idx, &nbrs, &doubles, &triples, &aromatics)
        if aromatics:
            return False                  # already aromatic here; not a Kekule ring
        if nbrs > 3:
            return False                  # over-coordinated for a ring atom, or a spiro atom
        if doubles == 1 and not triples:
            continue
        invalid = 0
        cls_out[idx] = arom_classify_atom(a.element, a.charge, False, nbrs, False,
                                          <int> at_implicit_h(a), &invalid)
        if invalid or cls_out[idx] == AROM_MUST:
            return False

    # at least one of the ring's OWN bonds is a double bond; see the docstring's borazine paragraph
    for i in range(size):
        idx = ring[i]
        nxt = ring[i + 1 if i + 1 < size else 0]
        for k in range(ptr[idx], ptr[idx + 1]):
            if edges[k].to == nxt and edges[k].order == 2:
                ring_double = True
                break
        if ring_double:
            break
    return ring_double


cdef void arom_thiele_mark(Structure structure, thiele_t *t, uint32_t *rings,
                           uint32_t nrings) noexcept nogil:
    """Rebuild the candidate edge set and its atom support from the surviving rings.

    From scratch every time rather than incrementally: dropping one ring can make another ring's
    bond exocyclic, so the set is a fixpoint and an incremental update would have to un-mark edges a
    dropped ring shared with a surviving one.
    """
    cdef uint32_t base = 2 + nrings
    cdef uint32_t i, k, u, v, size, slot
    cdef uint32_t *ptr = csr_ptr(structure)
    cdef halfedge_t *edges = csr_edges(structure)
    memset(t.cand, 0, structure.header.atom_count * sizeof(uint8_t))
    memset(t.he_arom, 0, ptr[structure.header.atom_count] * sizeof(uint8_t))
    for i in range(nrings):
        if not t.ring_ok[i]:
            continue
        size = rings[2 + i] - rings[1 + i]
        for k in range(size):
            u = rings[base + rings[1 + i] + k]
            v = rings[base + rings[1 + i] + (k + 1 if k + 1 < size else 0)]
            t.cand[u] = 1
            for slot in range(ptr[u], ptr[u + 1]):
                if edges[slot].to == v:
                    t.he_arom[slot] = 1
            for slot in range(ptr[v], ptr[v + 1]):
                if edges[slot].to == u:
                    t.he_arom[slot] = 1


cdef bint arom_thiele_exo(Structure structure, thiele_t *t, uint32_t i) noexcept nogil:
    """True when atom `i` carries a double or triple bond outside the candidate edge set."""
    cdef uint32_t *ptr = csr_ptr(structure)
    cdef halfedge_t *edges = csr_edges(structure)
    cdef uint32_t k
    cdef uint8_t order
    for k in range(ptr[i], ptr[i + 1]):
        order = edges[k].order
        if (order == 2 or order == 3) and not t.he_arom[k]:
            return True
    return False


cdef uint32_t arom_thiele_prune(Structure structure, thiele_t *t, uint32_t *rings,
                                uint32_t nrings) noexcept nogil:
    """Drop every candidate ring holding an atom whose pi electron went elsewhere, to a fixpoint.

    Returns the number of rings dropped.  THE RING IS DROPPED AND THE ATOM IS NOT: deleting the atom
    from the ring graph and re-running SSSR over what is left gives the same answer wherever the
    surviving aromatic ring is itself a member of the cycle basis -- naphthoquinone's benzene is, and
    so is every fused case in the fixtures.  Where the two could differ is a cycle that only appears
    once a quinone carbon is deleted, and such a cycle is longer than the basis rings it replaces,
    hence longer than 7, hence not a candidate here at all.  A second SSSR pass with no size bound
    can emit a ten-membered aromatic ring.
    """
    cdef uint32_t base = 2 + nrings
    cdef uint32_t i, k, size, dropped = 0
    cdef bint again = True
    while again:
        again = False
        for i in range(nrings):
            if not t.ring_ok[i]:
                continue
            size = rings[2 + i] - rings[1 + i]
            for k in range(size):
                if arom_thiele_exo(structure, t, rings[base + rings[1 + i] + k]):
                    t.ring_ok[i] = 0
                    dropped += 1
                    again = True
                    break
        if again:
            arom_thiele_mark(structure, t, rings, nrings)
    return dropped


cdef uint32_t arom_thiele_component(Structure structure, thiele_t *t, uint32_t seed,
                                    uint32_t *nedges) noexcept nogil:
    """One connected component of the candidate edge set, breadth first.  Atoms into `t.clist`."""
    cdef uint32_t *ptr = csr_ptr(structure)
    cdef halfedge_t *edges = csr_edges(structure)
    cdef uint32_t head = 0, count = 1, v, w, k
    t.clist[0] = seed
    t.comp[seed] = <int32_t> seed
    nedges[0] = 0
    while head < count:
        v = t.clist[head]
        head += 1
        for k in range(ptr[v], ptr[v + 1]):
            if not t.he_arom[k]:
                continue
            w = edges[k].to
            # each candidate bond is two half-edges; the v < w test counts it once
            if v < w:
                nedges[0] += 1
            if t.comp[w] < 0:
                t.comp[w] = <int32_t> seed
                t.clist[count] = w
                count += 1
    return count


cdef bint arom_thiele_check(Structure structure, thiele_t *t, uint32_t count,
                            uint32_t nedges, uint32_t *pi_out) noexcept nogil:
    """Is the molecule in hand a Kekule form of this component's edge set?

    Two questions.  The MATCHING: every MUST atom carries exactly one candidate double bond, every
    MUST_NOT carries none, a MAY atom either -- which is the statement that `kekule()` run on this
    edge set could have produced these orders.  Then HUCKEL, and only when the component is a single
    cycle (`nedges == count`, every atom of degree 2): pi must be 2 mod 4.  A fused component gets no
    Huckel test, because azulene's five-ring fails one and azulene is aromatic.
    """
    cdef uint32_t *ptr = csr_ptr(structure)
    cdef halfedge_t *edges = csr_edges(structure)
    cdef uint32_t j, i, k, deg, dbl
    cdef atom_t *atoms = structure.atoms()
    cdef bint cycle = nedges == count
    pi_out[0] = 0
    for j in range(count):
        i = t.clist[j]
        if t.cls[i] == THIELE_INVALID:
            return False
        deg = 0
        dbl = 0
        for k in range(ptr[i], ptr[i + 1]):
            if not t.he_arom[k]:
                continue
            deg += 1
            if edges[k].order == 2:
                dbl += 1
        if deg != 2:
            cycle = False
        if t.cls[i] == AROM_MUST:
            if dbl != 1:
                return False
        elif t.cls[i] == AROM_MUST_NOT:
            if dbl:
                # `kekule()` would never have put a ring double bond on a lone-pair donor, so these
                # orders are not a Kekule form of this edge set
                return False
        elif dbl > 1:
            return False                  # MAY takes one or none, never two
        pi_out[0] += arom_thiele_pi(atoms + i, dbl == 1)
    if cycle and pi_out[0] % 4 != 2:
        return False
    return True


def thiele(MoleculeContainer mol not None):
    """Turn Kekule bond orders into aromatic ones.  `MoleculeContainer.thiele`.

    A DELIBERATE operation and one of the two in the library allowed to change a molecule's
    representation (`kekule` is the other).  No writer calls it: a molecule that holds Kekule orders
    is written Kekule, and a caller who wants the aromatic spelling asks for it here first.

    The rule is that the molecule in hand must be a Kekule form of the aromatic edge set this
    returns, checked against `arom_classify_atom` -- the same table `kekule()` classifies with, so
    the two directions cannot drift apart.  Necessary but not sufficient, so three more conditions
    apply: ring size 5 to 7, no exocyclic double or triple bond on a ring atom, and Huckel's 4n+2 on
    a candidate system that is a single cycle.  The header of this file has the two molecules that
    prove the matching alone is not enough and the two that prove Huckel cannot be applied per ring.

    Tautomers are not touched.  Shifting a pyrrole hydrogen around a condensed ring while aromatising,
    or patching `N1C=Cn2cccc12` with a SMARTS list, moves an atom's hydrogen count: a tautomer decision,
    and the standardization pack's business.  A system this function cannot spell aromatic comes back in
    `.refused` instead.

    Call this with the journal clean.  It opens its own edit scope, so on return the orders are
    applied unless the caller holds an outer scope, in which case they apply when that closes.

    Returns a `ThieleResult`.  `.changed` is measured against what the arena holds, so a second call
    is observably a no-op rather than a promise of one.  `.refused` lists the candidate systems
    declined, one tuple of stable ids each, with a `.log` line naming why.

    A ring the PRE-FILTER declines is reported the same way, as its own tuple, because it never
    becomes a system for the component loop to name -- but only for the two states that are a refusal
    rather than a non-candidacy: a radical, and an unknown implicit hydrogen count.  Those two used
    to return an empty `ThieleResult`, which mattered most for the molecule this whole layer is
    built around: between a read and a `kekule()` an unknown count on the pyrrole-versus-pyridine
    atom is the ORDINARY state, and a caller got a Kekule molecule back with nothing saying why.  A
    non-aromatizable element, an already-aromatic bond and a spiro atom stay silent; see
    `THIELE_WHY_*`.
    """
    mol._require_clean()
    cdef Structure structure = mol._structure
    cdef uint32_t n = structure.header.atom_count
    cdef list numbers = mol._numbers
    cdef list log = []
    cdef list refused = []
    if not n or not structure_has(structure, SEG_RELEVANT_RINGS):
        return arom_thiele_result(mol, False, log, refused)

    cdef uint32_t *rings = structure_rings(structure)
    cdef uint32_t nrings = rings[0]
    if not nrings:
        return arom_thiele_result(mol, False, log, refused)

    cdef uint32_t *ptr = csr_ptr(structure)
    cdef halfedge_t *edges = csr_edges(structure)
    cdef _ThieleRun run = _ThieleRun.__new__(_ThieleRun)
    run.n = n
    run.nrings = nrings
    arom_thiele_alloc(&run.t, n, ptr[n], nrings)
    cdef thiele_t *t = &run.t

    # --- the per-ring pre-filter
    cdef uint32_t base = 2 + nrings
    cdef uint32_t i, j, k, size, slot, count
    cdef uint32_t nedges = 0, pi = 0
    cdef bint any_ring = False
    cdef uint8_t why = THIELE_WHY_NONE
    cdef uint32_t culprit = 0
    cdef list ring_ids
    for i in range(nrings):
        size = rings[2 + i] - rings[1 + i]
        if arom_thiele_ring_ok(structure, rings + base + rings[1 + i], size, t.cls,
                               &why, &culprit):
            t.ring_ok[i] = 1
            any_ring = True
        elif why:
            # A ring declined here never becomes a component, so the loop below will never name it
            # and this is the only place it can be reported.  Reported per RING and not per system
            # for the same reason: there is no system yet.
            ring_ids = []
            for j in range(size):
                ring_ids.append(numbers[(rings + base + rings[1 + i])[j]])
            ring_ids.sort()
            refused.append(tuple(ring_ids))
            if why == THIELE_WHY_RADICAL:
                log.append(mc_record('thiele:radical-ring', (numbers[culprit],),
                                     f'ring {tuple(ring_ids)!r} carries a radical at atom '
                                     f'{numbers[culprit]}; a radical brings one pi electron and not two, so '
                                     f'the ring is left as it is',
                                     mc_refused()))
            else:
                log.append(mc_record('thiele:unknown-h', (numbers[culprit],),
                                     f'ring {tuple(ring_ids)!r} has an unknown implicit hydrogen count at '
                                     f'atom {numbers[culprit]}; that count IS the atom\'s aromatic class, so '
                                     f'the ring is left as it is -- derive_hydrogens() first',
                                     mc_refused()))
    if not any_ring:
        return arom_thiele_result(mol, False, log, refused)

    arom_thiele_mark(structure, t, rings, nrings)
    arom_thiele_prune(structure, t, rings, nrings)

    # --- classify the surviving support.  `arom_thiele_ring_ok` already classified the non-sp2
    # atoms, but a ring it dropped may have left one behind, so every atom is asked here and the
    # pre-filter's answers are not carried over.
    cdef atom_t *atoms = structure.atoms()
    cdef atom_t *a
    cdef uint32_t nbrs = 0, doubles = 0, triples = 0, aromatics = 0
    cdef uint8_t invalid
    cdef bint support = False
    for i in range(n):
        t.comp[i] = -1
        if not t.cand[i]:
            continue
        support = True
        a = atoms + i
        arom_thiele_bonds(structure, i, &nbrs, &doubles, &triples, &aromatics)
        invalid = 0
        # `exo_double` is False by construction: `arom_thiele_prune` dropped every ring holding an
        # atom with a double bond outside the candidate set, so an atom still in the support has none
        t.cls[i] = arom_classify_atom(a.element, a.charge, False, nbrs, False,
                                      <int> at_implicit_h(a), &invalid)
        if invalid:
            t.cls[i] = THIELE_INVALID
    if not support:
        return arom_thiele_result(mol, False, log, refused)

    # --- one component at a time: accept it, or drop its edges and name it
    cdef list names
    for i in range(n):
        if not t.cand[i] or t.comp[i] >= 0:
            continue
        nedges = 0
        count = arom_thiele_component(structure, t, i, &nedges)
        if arom_thiele_check(structure, t, count, nedges, &pi):
            # ACCEPTED, AND SAYING SO IS THE POINT.  The refusal below was the only thing this pass
            # reported, so a molecule it aromatised came back `changed=True` with an empty log.  The
            # accepted system is INFO and the declined one REFUSED, which is the pair a caller reads.
            names = []
            for j in range(count):
                names.append(numbers[t.clist[j]])
            names.sort()
            log.append(mc_record('thiele:aromatized', tuple(names),
                                 f'ring system {tuple(names)!r} written aromatic, {pi} pi electrons'))
            continue
        names = []
        for j in range(count):
            k = t.clist[j]
            names.append(numbers[k])
            # un-mark every half-edge AT the atom, not only the candidate ones: the emit loop below
            # reads the set back, and one half left marked would write an order-4 bond for a system
            # this function has just declined.  Clearing all of them clears both halves of every
            # candidate bond, since the component is closed under them.
            for slot in range(ptr[k], ptr[k + 1]):
                t.he_arom[slot] = 0
            t.cand[k] = 0
        names.sort()
        refused.append(tuple(names))
        log.append(mc_record('thiele:not-kekule-form', tuple(names),
                             f'ring system {tuple(names)!r} is not a Kekule form of an aromatic system; '
                             f'left as it is',
                             mc_refused()))

    # --- emit
    cdef list pairs = []
    cdef bint changed = False
    cdef object pair
    for i in range(n):
        for k in range(ptr[i], ptr[i + 1]):
            if t.he_arom[k] and i < edges[k].to:
                pairs.append((numbers[i], numbers[edges[k].to]))
                if edges[k].order != 4:
                    changed = True
    if pairs:
        # The same exemption `kekule` takes, for the same reason and stated there: this changes the
        # representation and not the molecule, so stored CIP descriptors survive it.
        mol._representation_change = True
        try:
            with mol.edit():
                for pair in pairs:
                    mol.set_order(pair[0], pair[1], 4)
        finally:
            mol._representation_change = False
    return arom_thiele_result(mol, changed, log, refused)
