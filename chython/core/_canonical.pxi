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
# Canonical machinery over the molecule arena: the symmetry orbits, and the canonical atom order.
#
# This is the molecule-side twin of `_compute_automorphisms` in _query_seal.pxi. The two share
# their shape -- refine to a partition, branch over same-class candidates, verify every candidate
# exactly -- and they must stay recognisable as one algorithm. What differs is smaller than it
# looks: two comparators (`_wterm_equal` over the query arena's `wterm_t` against
# `_atom_colour_equal`/`_bond_colour_equal` over `atom_t`/`halfedge_t`) and one adjacency accessor
# (the query's csr_head/csr_nbr/csr_bnd triple against this side's csr_ptr/csr_edges half-edges).
# A shared kernel taking a context struct plus function pointers is possible; it is not worth it
# yet, because the two now differ in POLICY as well as in types -- the query side enumerates the
# group in slot order and keeps a prefix of the rows, this side asks one decisive question per
# unresolved pair (see below) -- and a kernel that parameterises policy too would be harder to
# read than either copy. If a third copy appears, revisit.
#
# THE VERIFICATION IS THE AUTHORITY; THE PARTITION IS ONLY A PRUNER. `compute_atoms_order`
# hashes, so its classes may be coarser than the true equitable partition on a collision; that
# costs candidates the verification then rejects and nothing else. It can never be too FINE,
# because equal inputs hash equal, so an automorphism-related pair is never split.
#
# What an automorphism must preserve here: each atom's intrinsic record -- element, isotope,
# charge, radical, implicit hydrogen count -- and each bond's order together with its aromatic bit.
# Deliberately NOT the map number (AAM is annotation, and stereogenicity is a property of the
# structure, not of how someone numbered it), not the stereo flag (the whole point of this
# fragment is to decide stereogenicity, so it cannot presuppose it), and not the derived fields
# -- degree, heteroatoms, ring sizes and counts, and the EXPLICIT hydrogen nibble, are all
# functions of the graph, so a graph automorphism preserves them for free.
#
# WHERE THIS DIVERGES FROM THE QUERY TWIN, AND WHY. The query side enumerates the whole group in
# slot order and keeps the first 1024 rows it finds. That is fine for its consumer, which only
# needs SOME subgroup to filter duplicate matches with, but it is wrong for the consumer here.
# Six waters in one record fill the row cap with permutations of the last few atoms and never
# reach the ones that move the first, so the union of the rows reports five orbits where there
# are two -- an under-report of symmetry, which for stereo means inventing stereocentres. So this
# side asks the question its consumer actually asks: for each pair of atoms not yet known to
# share an orbit, does SOME automorphism map one onto the other? Each such search is exhaustive
# and therefore decisive, each success merges orbits, and a pair already merged is skipped, so
# at most atom_count - 1 searches can succeed. THE ORBITS ARE THE ANSWER THOSE SEARCHES BUILD.
#
# THIS SIDE DELIBERATELY NEVER MATERIALISES THE GROUP. The only question it answers is "do these
# two atoms share an orbit", and the union-find that accumulates the answers is the whole
# computation. It hands back no permutations: the ones such a search finds are orbit-complete
# without being a generating set, so nothing could be trusted on top of them. A consumer that needs
# actual group elements must ADD a search that produces them, with its own termination argument --
# not read a by-product of this one.


# Two budgets, because they protect against two different things.
#
# CANON_MAX_NODES_SEARCH bounds ONE pair search. It is what makes the answer for a pair a
# function of the structure and that pair alone, and nothing else -- not of how many searches
# happened to run before it. That is the property orbits need: a shared per-call budget made the
# orbits depend on atom numbering (11 cyclopropanes + 11 cyclobutanes, n = 77: 45 orbits when the
# rings are built blocked, 2 when interleaved, truth 2).
#
# CANON_MAX_NODES_CALL bounds the whole call, because a call runs up to O(atom_count^2) searches
# and the per-search bound alone is not a wall-time bound. Hitting it stops the call early, which
# is order-dependent again -- so it is the honest ceiling of last resort, and like the per-search
# bound it sets CANON_BUDGET_EXCEEDED. Both are counted in candidates CONSIDERED (a candidate
# that survives the class prune and reaches the colour test).
DEF CANON_MAX_NODES_SEARCH = 100000
DEF CANON_MAX_NODES_CALL = 20000000

# A typed global rather than a DEF, like Q_NO_SLOT on the query side: the search compares against
# it inside `nogil`, and a DEF's value is a Python int there.
cdef uint32_t CANON_NO_SLOT = 0xFFFFFFFF


cdef enum:
    CANON_ASYMMETRIC = 1          # the group is trivial: only the identity
    CANON_BUDGET_EXCEEDED = 2     # a search was truncated, so this result is not exact


# THE STEREO SEAM, AND WHY THE SEARCH NEEDED ONE AT ALL.
#
# This fragment is included before `_stereo.pxi` and deliberately does not know what a parity is --
# `mol_certificate_words` already says so, and takes its stereo term as an opaque `extra` array from a
# caller who does.  That worked for as long as the certificate was the only thing the parity reached.
# It is not: the SEARCH has to see the parity too, or it picks between two labellings the constitution
# cannot separate by SLOT ORDER, and the whole extremal construction is a function of the caller's
# input order again.  The witness is cis-cyclobutane-1,3-diol.  Its two carbinol carbons are one
# refinement class, the constitutional automorphism that exchanges them INVERTS the parity at both, and
# the orbit prune below therefore drops one of the two labellings as "interchangeable" when it is not.
# Which one survives is the lower slot.  Measured before this seam existed, over four creation orders
# of the 294 usable records of `test/stereo.sdf`: 55 records returned two or more `canonical_bytes`,
# 43 two or more canonical SMILES, and 45 failed `write -> read -> ==` against themselves.
#
# So the parity enters through two function pointers, filled in by `_stereo.pxi` at import.  Pointers
# rather than parameters threaded through four call sites, because every caller wants the same answer:
# a canonical order that depended on whether the caller remembered to pass the hook would be two
# canonical orders, and `mol_identity_bytes` and the SMILES writer would disagree about one molecule.
# NULL is a working configuration and means "stereo-blind", which is what the fragment did before.
#
#   * `_canon_prepare_hook` is called ONCE per `_canon_order`, before any arena pointer is taken,
#     because building the unit table appends a segment and reallocates (ruling F60).  It holds the
#     GIL and may raise.
#   * `_canon_stereo_hook` is called inside the search, under `nogil`, and may not touch the arena's
#     shape.  Given a colouring it fills one parity digit per slot -- ruling F95's frame-free code,
#     read in the frame that colouring names -- and returns whether that colouring NAMES every
#     configured unit's frame.  The digits go in the certificate's tail; the boolean gates the orbit
#     prune, which is only entitled to call two candidates interchangeable when it can prove the
#     symmetry relating them preserves the configuration and not merely the graph.
#   * `unnamed_out`, when it is not NULL, is atom_count bytes the hook marks with WHICH ATOMS the
#     False answer is about: an unnamed unit's anchor, its partner and its directions.  That turns the
#     gate from per node into per symmetry -- `_canon_branch` prunes when the node's group fixes every
#     marked atom, since a symmetry that moves none of a frame's atoms cannot invert its parity.
#
# `scratch` is 2 * atom_count uint32_t owned by the search, so the hook allocates nothing on a path
# that runs once per tree node.
ctypedef bint (*canon_stereo_fn)(Structure structure, uint32_t *colour, uint32_t *digits_out,
                                 uint32_t *scratch, uint8_t *unnamed_out) noexcept nogil
ctypedef int (*canon_prepare_fn)(Structure structure) except -1

cdef canon_stereo_fn _canon_stereo_hook = NULL
cdef canon_prepare_fn _canon_prepare_hook = NULL


with cython.warn.undeclared(False):
    # bare so Python can import it, guarded so warn.undeclared stays quiet
    class AutomorphismBudgetExceeded(RuntimeError):
        """A symmetry search hit its node budget, so its answer is not exact.

        Orbits from a truncated search can be FINER than the truth -- atoms that a symmetry does
        relate reported as unrelated -- which for stereo perception means inventing stereocentres.
        There is no safe degraded answer to hand back, so `automorphism_orbits` raises this
        instead of returning one. `is_asymmetric` does not raise: False there already means "not
        known to be asymmetric".

        `canonical_order` raises it too, for a different reason with the same conclusion: a
        truncated extremal search returns SOME labelling in place of THE labelling, so every hash
        and equality derived from it would be silently wrong. One exception type, because from a
        caller's side the fact is the same one -- a symmetry search ran out of nodes and there is
        nothing to hand back.
        """


cdef inline uint32_t _uf_find(uint32_t *parent, uint32_t i) noexcept nogil:
    """Union-find root, with path halving. The orbits of a group are the connected components of
    any generating set, so mol_automorphisms never has to materialise a group closure."""
    while parent[i] != i:
        parent[i] = parent[parent[i]]
        i = parent[i]
    return i


cdef inline void _uf_emit(uint32_t n, uint32_t *parent, uint32_t *orbits) noexcept nogil:
    """Write the union-find components into `orbits` as dense 1-based ids."""
    cdef uint32_t i, r, next_id = 0
    for i in range(n):
        if _uf_find(parent, i) == i:
            next_id += 1
            orbits[i] = next_id
    for i in range(n):
        r = _uf_find(parent, i)
        if r != i:
            orbits[i] = orbits[r]


cdef inline bint _atom_colour_equal(atom_t *a, atom_t *b) noexcept nogil:
    """Do two atom records carry the same intrinsic colour? See the fragment comment for the
    field list and for why map number, the stereo flag and the derived fields are all absent.

    Only the IMPLICIT hydrogen nibble is compared: the explicit one counts H neighbours, and a
    bijection that preserves elements and bonds preserves that count already, so comparing it
    could never reject anything -- and comparing it anyway would contradict the rule this
    fragment states about derived fields.

    RAW on the sentinel, deliberately unbranched.  An equality test over the nibble gives H_UNKNOWN
    the semantics an automorphism needs for free: unknown maps onto unknown, and never onto a stated
    count.  Mapping the two together would claim that an atom whose hydrogens were recorded and one
    whose were not are the same atom, which is exactly the claim a symmetry may not make -- and it
    would make the canonical form of a half-recorded molecule depend on which half was recorded.
    """
    # The R index is part of the intrinsic colour too: an automorphism may not map an R1 onto an R2.
    return (a.element == b.element and a.isotope == b.isotope and a.charge == b.charge
            and at_implicit_h(a) == at_implicit_h(b) and at_radical(a) == at_radical(b)
            and at_r_index(a) == at_r_index(b))


cdef inline bint _bond_colour_equal(halfedge_t *e, halfedge_t *f) noexcept nogil:
    """Bond colour is the Kekule order plus the aromatic bit. HE_IN_RING is derived from the
    topology an automorphism already preserves, so it is not part of the colour."""
    return e.order == f.order and (e.flags & HE_AROMATIC) == (f.flags & HE_AROMATIC)


cdef void _canon_search_order(uint32_t n, uint32_t *ptr, halfedge_t *edges, uint32_t *pin,
                              uint32_t *order, uint32_t *anchor, uint8_t *seen) noexcept nogil:
    """The order the search assigns slots in: every PINNED slot first, then breadth-first from them.

    `pin` is the per-slot pin array the search is about to run with -- `pin[s]` is the required
    image of slot `s`, or CANON_NO_SLOT to leave it free -- and AT LEAST ONE SLOT MUST BE PINNED.
    That is a precondition and not a check: both callers pin at least one slot (the pair search pins
    its source, `mol_find_pinned_begin` pins an anchor), and a search with nothing pinned would find
    the identity and answer nothing.

    Pinned slots lead because that is what keeps the pin cheap: they are assigned at their own
    depths, from a one-element candidate list each, before any free slot spends a candidate. Their
    `anchor` is CANON_NO_SLOT even when an earlier pinned slot is adjacent -- a pinned depth draws
    from its pin and never from a neighbour's image, so the anchor would only be read for the early
    bond test, and the verification walk covers that bond anyway.

    Every slot after the pinned prefix is adjacent to one already assigned, and `anchor` names
    which: anchor[d] is the DEPTH of that already-assigned neighbour, so order[anchor[d]] is the
    slot itself. Component roots get CANON_NO_SLOT. The search reads it twice over -- to draw its
    candidates from the neighbour's image instead of from all n slots, and to reject a wrong
    candidate at the shallowest depth that can see it.

    Disconnected components follow, each rooted at the lowest slot left. `seen` is scratch and is
    left dirty.
    """
    cdef uint32_t head = 0, tail = 0, v, u, k, d
    memset(seen, 0, <size_t> n)
    for u in range(n):
        if pin[u] != CANON_NO_SLOT:
            seen[u] = 1
            order[tail] = u
            anchor[tail] = CANON_NO_SLOT
            tail += 1
    while True:
        while head < tail:
            d = head
            v = order[head]
            head += 1
            for k in range(ptr[v], ptr[v + 1]):
                u = edges[k].to
                if not seen[u]:
                    seen[u] = 1
                    order[tail] = u
                    anchor[tail] = d
                    tail += 1
        if tail == n:
            return
        for u in range(n):          # next component, rooted at the lowest slot not yet taken
            if not seen[u]:
                seen[u] = 1
                order[tail] = u
                anchor[tail] = CANON_NO_SLOT
                tail += 1
                break


# THE SEARCH'S RESUMABLE STATE. Every array is CALLER-OWNED scratch of `n` entries; the three
# fields after them are the cursor into the backtracking tree, which is what lets a caller ask for
# the NEXT automorphism consistent with the same pins instead of restarting the search.
#
# The pair search wants one answer and the stereogenicity predicate wants an enumeration -- a
# candidate automorphism that fails the stereo-consistency test is not a witness and its existence
# proves nothing, so only an EXHAUSTED enumeration is a decision there. One resumable kernel serves
# both, which is why there is no second copy of the walk below.
cdef struct pinned_search_t:
    uint32_t n
    uint32_t *cls       # refinement classes; a candidate must share its slot's class
    uint32_t *pin       # pin[s] is the required image of slot s, or CANON_NO_SLOT to leave it free
    uint32_t *order     # from _canon_search_order, built against this same `pin`
    uint32_t *anchor    # ditto
    uint32_t *sigma     # the assignment; complete and valid exactly when a step returned True
    uint32_t *cursor    # per-depth candidate cursor
    uint8_t *taken      # per-slot injectivity mark
    uint32_t depth
    bint live           # the tree still has unexplored nodes
    bint solved         # the last step returned a complete assignment, still standing in `sigma`
    bint truncated      # the budget ran out, so "no more" is not a decision


cdef void mol_find_pinned_begin(pinned_search_t *st) noexcept nogil:
    """Reset the state to the root of the tree. `st.order` and `st.anchor` must already be built
    (by `_canon_search_order`, against the same `st.pin`), and `st.pin` must not change afterwards:
    the order's pinned prefix is a function of it."""
    cdef uint32_t s
    memset(st.taken, 0, <size_t> st.n)
    for s in range(st.n):
        st.sigma[s] = CANON_NO_SLOT
    st.cursor[0] = 0
    st.depth = 0
    st.live = True
    st.solved = False
    st.truncated = False


cdef bint mol_find_pinned_next(uint32_t *ptr, halfedge_t *edges, atom_t *atoms,
                               pinned_search_t *st, uint32_t *budget) noexcept nogil:
    """The next automorphism respecting `st.pin`. Writes it into `st.sigma`, returns true.

    Backtracking over depths. Each accepted slot is verified against the edges it has to slots
    already assigned, so by the time every slot is assigned every edge of the graph has been
    checked exactly once -- and that is enough to make the assignment an automorphism, non-edges
    included, by the counting argument written out at the walk below. A complete assignment is
    therefore returned with no further check.

    Three candidate sources, one per depth:

      * a PINNED slot offers `pin[s]` and nothing else -- a one-element candidate list at the
        slot's OWN depth, which is how every constraint this search imposes gets imposed, and
        which is what keeps the counting argument below intact. See that argument for why
        pre-seeding `sigma` with the pins instead would silently break the search.
      * a depth with an anchor draws from the CSR adjacency of the anchor's image. Everything the
        BFS order reaches has an assigned neighbour, so nothing legal is outside that list, and
        the bond check to the anchor comes free with the half-edge the scan is already holding.
      * a component root has no anchor and scans the class over all n slots.

    A pinned depth ignores its anchor even when it has one, and takes the verification walk as its
    only bond test. That is not a weakening: the walk checks every edge from `s` to an assigned
    slot, the anchor bond among them.

    False means the enumeration is over: either no further automorphism exists, or the budget ran
    out and `st.truncated` says so. Exhausting the budget can only lose symmetry, never invent it.
    """
    cdef uint32_t n = st.n
    cdef uint32_t *cls = st.cls
    cdef uint32_t *pin = st.pin
    cdef uint32_t *order = st.order
    cdef uint32_t *anchor = st.anchor
    cdef uint32_t *sigma = st.sigma
    cdef uint32_t *cursor = st.cursor
    cdef uint8_t *taken = st.taken
    cdef uint32_t s, t, k, v, vimg, anc, cand = 0, img, base = 0, stop = 0
    cdef uint32_t depth = st.depth
    cdef halfedge_t *e2
    cdef halfedge_t *anchor_src = NULL     # the bond s -- order[anc], fixed for the whole depth
    cdef halfedge_t *anchor_dst = NULL     # its image candidate, the half-edge the scan holds
    cdef bint ok, found, pinned, exhausted = False, result = False

    if not st.live:
        return False
    if st.solved:
        # RESUME. The previous step left a complete assignment standing; undo just the deepest one
        # and carry on scanning that depth from the cursor it stopped at -- which is exactly what
        # the backtrack arm at the bottom of the loop does.
        st.solved = False
        depth = n - 1
        taken[sigma[order[depth]]] = 0
        sigma[order[depth]] = CANON_NO_SLOT

    while True:
        found = False
        s = order[depth]
        anc = anchor[depth]
        pinned = pin[s] != CANON_NO_SLOT
        if anc == CANON_NO_SLOT or pinned:
            anchor_src = NULL
        else:
            img = sigma[order[anc]]
            base = ptr[img]
            stop = ptr[img + 1]
            anchor_src = csr_find_at(ptr, edges, s, order[anc])  # exists: the BFS order says so
        t = cursor[depth]
        while True:
            if pinned:
                if t:                    # the pin was already tried, and it is the only candidate
                    break
                cand = pin[s]
            elif anc == CANON_NO_SLOT:
                if t >= n:
                    break
                cand = t
            else:
                if base + t >= stop:
                    break
                anchor_dst = &edges[base + t]
                cand = anchor_dst.to
            t += 1
            if taken[cand] or cls[cand] != cls[s]:
                continue
            if budget[0] == 0:
                exhausted = True
                break
            budget[0] -= 1
            if not _atom_colour_equal(&atoms[s], &atoms[cand]):
                continue
            # The edge walk below checks the anchor bond too; doing it here first costs nothing
            # (the half-edge is already in hand) and rejects most candidates before the walk. The
            # walk stays the authority -- verification is never narrowed to one bond.
            if anchor_src is not NULL and not _bond_colour_equal(anchor_src, anchor_dst):
                continue
            # VERIFY EACH EDGE ONCE. The walk is over s's OWN adjacency, and it checks only the
            # neighbours already assigned: O(degree) per candidate, where walking the assigned
            # prefix and looking both ways cost O(depth) lookups per candidate and made the whole
            # call O(n^3) on a symmetric record.
            #
            # WHY THAT IS NOT A WEAKENING, THOUGH IT LOOKS LIKE ONE. The case it appears to drop
            # is a NON-EDGE: s has no edge to some assigned slot v while cand has one to sigma[v].
            # The source-side walk never looks at that pair. It cannot survive to a complete
            # assignment, by counting. Over the full assignment the walk verifies every one of
            # the graph's m edges exactly once -- edge {a, b} at the depth of whichever endpoint
            # is assigned later, walked from that endpoint towards the earlier one, which is
            # assigned by then -- and it maps them into edges of the SAME graph, injectively,
            # because `taken` keeps sigma injective. An injective map from an m-set into an m-set
            # is onto, so the image edge set is exactly the edge set. A spurious edge
            # {sigma[a], sigma[b]} over a non-edge {a, b} would have to be the image of some real
            # edge {c, d}, and injectivity then forces {a, b} == {c, d}. So non-edges are verified
            # for free. This is why the walk must be over s's edges specifically, and why EVERY
            # edge to an already-assigned neighbour must be checked, not just the anchor's.
            #
            # TWO PREMISES THE COUNT RESTS ON, NEITHER ENFORCED HERE. The graph must be simple. A
            # self-loop {a, a} is verified by no depth at all -- it has no later endpoint -- so the
            # count would fall short of m. Parallel edges break edge-injectivity, since two edges
            # on the same pair of slots have the same image. Both are excluded upstream, not here:
            # `add_bond` rejects a self bond and a duplicate pair, the journal's CSR build emits
            # one half-edge per pair, and the import validator requires each atom's `to` list to be
            # strictly increasing. If any of that is ever relaxed, this argument is the casualty.
            #
            # WHY THE PIN IS A CANDIDATE LIST AND NOT A PRE-SEEDING OF SIGMA. The counting argument
            # needs every slot to be assigned at its own depth with this walk running for it, which
            # is why the pinned arm above is a ONE-ELEMENT CANDIDATE LIST at the pinned slot's own
            # depth. Pre-seeding `sigma` with the pins and starting the depth loop after them looks
            # equivalent and is not: the edges INTERNAL to the pinned set would be verified by no
            # depth at all, because the later-endpoint walk that covers them never runs. For a
            # stereo unit those internal edges are the anchor-to-direction bonds, so there are
            # always four of them and never zero. The count falls to m - p, an injection from an
            # (m - p)-set into an m-set need not be onto, spurious edges become possible, and the
            # search starts accepting non-automorphisms WITH NO SYMPTOM. Pinning at depth restores
            # the count. The pinned assignment is still CHECKED here -- atom colour and the walk --
            # and not merely installed.
            #
            # Pruning strength on EDGES is unchanged: the same edge is still checked at the same
            # depth, just once instead of once per later depth. The one thing that moves is WHEN a
            # spurious-edge candidate dies -- at the depth where the compensating missing edge shows
            # up rather than the one that first saw the non-edge, which is no later than the depth
            # that assigns s's last neighbour, since by then the
            # walk has checked all deg(s) of them into cand's own adjacency. Measured on the
            # symmetric-cycloalkane family the candidate count does not move at all: the per-call
            # budget starts firing at exactly the same atom count as before, 4200, while wall time
            # for the records under it drops 3-7x.
            ok = True
            for k in range(ptr[s], ptr[s + 1]):
                v = edges[k].to
                vimg = sigma[v]
                if vimg == CANON_NO_SLOT:
                    continue            # v's own depth will verify this edge, from the other end
                e2 = csr_find_at(ptr, edges, cand, vimg)
                if e2 is NULL or not _bond_colour_equal(&edges[k], e2):
                    ok = False
                    break
            if not ok:
                continue
            sigma[s] = cand
            taken[cand] = 1
            found = True
            break
        cursor[depth] = t
        if exhausted:
            break
        if found:
            depth += 1
            if depth == n:
                result = True
                break
            cursor[depth] = 0
        elif depth == 0:
            break
        else:
            depth -= 1
            taken[sigma[order[depth]]] = 0
            sigma[order[depth]] = CANON_NO_SLOT
    st.depth = depth
    st.solved = result
    if exhausted:
        st.truncated = True
        st.live = False     # a truncated enumeration cannot be resumed into a decision
    elif not result:
        st.live = False
    return result


cdef int mol_automorphisms(Structure structure, uint32_t *seed,
                           uint32_t *orbits_out, uint32_t *flags_out) except -1:
    """The molecule's symmetry orbits.

    `orbits_out` is the deliverable, and the only one. It is a CALLER-SUPPLIED array of atom_count
    uint32_t which receives the orbit partition as dense 1-based ids: two atoms share an id exactly
    when some automorphism maps one onto the other. Which orbit gets which number carries no
    meaning. Pass NULL to ask for flags_out only -- then no orbits are produced and the search
    stops at the first automorphism it finds, which is all the CANON_ASYMMETRIC bit needs.

    THE PERMUTATIONS THEMSELVES ARE NOT RETURNED, deliberately: see the fragment comment. Each
    successful search merges its whole permutation into the union-find and then drops it, so the
    partition is exact while no group element and no closure ever outlives the search that found
    it. A consumer needing group elements must add a search that produces them.

    `seed` is forwarded to compute_atoms_order: NULL starts the colouring from the atom records,
    otherwise it is per-atom starting labels, which is how stereo perception feeds its own
    distinctions back in.

    flags_out[0] receives CANON_ASYMMETRIC exactly when the group is known to be trivial, and
    CANON_BUDGET_EXCEEDED when any search was truncated -- see the two budget constants. A
    truncated result is not exact: its orbits may be FINER than the truth, never coarser, so
    CANON_ASYMMETRIC is withheld whenever the bit is set, even with no merge at all. Consumers
    must treat the bit as "do not trust this partition"; `automorphism_orbits` raises on it.
    """
    cdef uint32_t n = structure.header.atom_count
    cdef uint32_t *ptr
    cdef halfedge_t *edges
    cdef atom_t *atoms
    cdef uint32_t *cls = NULL
    cdef uint32_t *parent = NULL
    cdef uint32_t *order = NULL
    cdef uint32_t *anchor = NULL
    cdef uint32_t *sigma = NULL
    cdef uint32_t *cursor = NULL
    cdef uint32_t *pin = NULL
    cdef uint8_t *taken = NULL
    cdef uint32_t a, b, s, ra, rb, budget, successes = 0
    cdef uint64_t spent = 0
    cdef Py_ssize_t classes
    cdef bint ok = False, truncated = False, stop = False, ordered = False
    cdef pinned_search_t st

    flags_out[0] = 0
    if n < 2:
        flags_out[0] = CANON_ASYMMETRIC
        if orbits_out is not NULL and n:
            orbits_out[0] = 1
        return 0

    # Read-only from here on: nothing in this fragment appends to the arena, so caching the
    # segment pointers cannot outlive a reallocation.
    ptr = csr_ptr(structure)
    edges = csr_edges(structure)
    atoms = structure.atoms()

    try:
        cls = <uint32_t *> PyMem_Malloc(<size_t> n * sizeof(uint32_t))
        parent = <uint32_t *> PyMem_Malloc(<size_t> n * sizeof(uint32_t))
        order = <uint32_t *> PyMem_Malloc(<size_t> n * sizeof(uint32_t))
        anchor = <uint32_t *> PyMem_Malloc(<size_t> n * sizeof(uint32_t))
        sigma = <uint32_t *> PyMem_Malloc(<size_t> n * sizeof(uint32_t))
        cursor = <uint32_t *> PyMem_Malloc(<size_t> n * sizeof(uint32_t))
        pin = <uint32_t *> PyMem_Malloc(<size_t> n * sizeof(uint32_t))
        taken = <uint8_t *> PyMem_Malloc(<size_t> n * sizeof(uint8_t))
        if (cls is NULL or parent is NULL or order is NULL or anchor is NULL or sigma is NULL
                or cursor is NULL or pin is NULL or taken is NULL):
            raise MemoryError('automorphism scratch allocation failed')

        with nogil:
            classes = compute_atoms_order(structure, cls, seed)
        if classes < 0:
            raise MemoryError('atom order refinement failed to allocate')
        if classes == <Py_ssize_t> n:
            # Every class a singleton. An automorphism preserves the refinement, so only the
            # identity is left. The common case, and it costs one refinement and no search.
            flags_out[0] = CANON_ASYMMETRIC
            if orbits_out is not NULL:
                for a in range(n):
                    orbits_out[a] = a + 1
            return 0

        for a in range(n):
            parent[a] = a
            pin[a] = CANON_NO_SLOT
        st.n = n
        st.cls = cls
        st.pin = pin
        st.order = order
        st.anchor = anchor
        st.sigma = sigma
        st.cursor = cursor
        st.taken = taken
        for a in range(n):
            if stop:
                break
            ordered = False         # the walk order depends on the source slot only, so it is
            for b in range(a + 1, n):               # built once per source -- and only for a
                if cls[b] != cls[a]:                # source that some pair actually searches
                    continue
                ra = _uf_find(parent, a)
                rb = _uf_find(parent, b)
                if ra == rb:
                    # Already known to share an orbit, and an orbit is an equivalence class, so
                    # there is nothing an automorphism mapping a onto b could add to the
                    # partition. Skipping these is what keeps the successful searches to n - 1.
                    continue
                if spent >= CANON_MAX_NODES_CALL:
                    truncated = True
                    stop = True
                    break
                # ONE pinned slot, `a`, whose required image is `b`. That is Task 1's pair search
                # expressed in the general kernel: the order's pinned prefix is then `a` alone, so
                # it is a BFS rooted at `a` and depends on `a` only -- which is why `ordered` can
                # hoist it out of the `b` loop.
                pin[a] = b
                if not ordered:
                    with nogil:
                        _canon_search_order(n, ptr, edges, pin, order, anchor, taken)
                    ordered = True
                budget = CANON_MAX_NODES_SEARCH
                with nogil:
                    mol_find_pinned_begin(&st)
                    ok = mol_find_pinned_next(ptr, edges, atoms, &st, &budget)
                pin[a] = CANON_NO_SLOT
                spent += <uint64_t> CANON_MAX_NODES_SEARCH - budget
                if not ok:
                    if st.truncated:
                        # This pair is undecided. Every OTHER pair still gets its own budget, so
                        # one expensive pair no longer coarsens -- or refines -- the whole answer.
                        truncated = True
                    continue                # decisive: a and b are in different orbits
                successes += 1
                if orbits_out is NULL:
                    stop = True             # one automorphism settles the only bit asked for
                    break
                for s in range(n):          # merge every cycle, not just the pair asked about
                    ra = _uf_find(parent, s)
                    rb = _uf_find(parent, sigma[s])
                    if ra != rb:
                        parent[ra] = rb

        if truncated:
            flags_out[0] = CANON_BUDGET_EXCEEDED
        elif successes == 0:
            # The refinement left a class with more than one member, but every search was
            # decisive and found nothing: the group really is trivial. A regular graph whose
            # vertices are not equivalent lands here.
            flags_out[0] = CANON_ASYMMETRIC
        if orbits_out is not NULL:
            with nogil:
                _uf_emit(n, parent, orbits_out)
    finally:
        PyMem_Free(cls)
        PyMem_Free(parent)
        PyMem_Free(order)
        PyMem_Free(anchor)
        PyMem_Free(sigma)
        PyMem_Free(cursor)
        PyMem_Free(pin)
        PyMem_Free(taken)
    return 0


# THE CANONICAL ORDER: AN EXTREMAL LABELLING, NOT MERELY A DISCRETE ONE.
#
# What the rest of this file computes is a partition. What follows computes an ORDER, and the
# distinction is the whole reason it exists. Refinement gets to a discrete partition on most
# records, and a discrete partition already numbers every atom -- but on a symmetric record it
# stops short, and the obvious repair (individualise some atom in an unresolved cell, refine, and
# take the first candidate that lands discrete) yields A discrete labelling that depends on which
# atom the input happened to put first. Measured, that is sixty random relabelings of one cubane
# skeleton answering forty distinct canonical strings. Nothing can be hashed or compared on top of
# that.
#
# So the search below takes an EXTREMUM instead of the first member of the unresolved cell: the
# greatest certificate over the candidates that no invariant reason rules out, which is a subset of
# the cell and not all of it -- `_canon_branch` states exactly which subset and why narrowing it
# keeps the result a function of the structure. The tree is the standard one:
#
#   * refine to an equitable partition (`compute_atoms_order`);
#   * if it is discrete, that is the labelling -- one refinement, no search, and this is the
#     common case for real molecules;
#   * otherwise pick the target cell, and for every candidate in it that no INVARIANT REASON rules
#     out -- see the two prunes in `_canon_branch`, which are the node's own symmetry orbits and a
#     node invariant -- individualise the candidate, refine again, and recurse; keep the labelling
#     whose certificate is greatest.
#
# WHY THAT IS A FUNCTION OF THE STRUCTURE. Every choice the tree makes is computed from the
# coloured graph and nothing else, so relabelling the input relabels the tree and leaves its shape
# alone. Concretely: refinement class ids come from `_classify`, which sorts by key value, and the
# keys are hashes of the atom invariant and of neighbour class multisets -- no slot number reaches
# them. The target cell is the LOWEST CLASS ID with more than one member, which is a choice among
# invariant ids rather than among slots. The individualising seed is the parent colouring with one
# atom moved to a fresh label, so it too corresponds under any relabelling. And the certificate is
# read off positions, never off slots. The maximum of an invariant function over an invariant set
# is invariant, which is the entire argument.
#
# WHAT IS STILL FREE, AND WHY THAT IS FINE. When several leaves tie on the certificate, the first
# one found wins, and which one that is does depend on slot order. Two tying leaves differ by an
# automorphism of the molecule, so the RELABELLED GRAPH -- the canonical form, and everything
# hashed from it -- is identical either way; only which of two interchangeable atoms got which of
# their two interchangeable positions moves. That is inherent: a canonical form is unique, a
# canonical labelling is unique only up to the automorphism group. A caller that needs a labelling
# pinned down further must feed the distinction in through `seed`.
#
# ONE COMPONENT AT A TIME. A record with more than one connected component is canonicalised component
# by component, and the blocks are laid end to end in the order of the components' OWN certificates
# (`_canon_order_split`). Refinement colours are graph-local -- a methyl in component 1 and a methyl in
# component 200 share a class however long refinement runs -- so a whole-record search individualises a
# candidate in every component in turn and calls `mol_automorphisms` over the whole record at every one
# of those nodes. Measured on N copies of nitroethane, canonical order only:
#
#   | components | whole record | decomposed |
#   | 128        |     0.214 s  |  0.148 ms  |
#   | 512        |      45.4 s  |  0.669 ms  |
#   | 1000       |     639.7 s  |  1.350 ms  |
#   | 4000       |            - |  6.107 ms  |
#
# THE DECOMPOSITION IS NOT FREE ON A SMALL MIXTURE, and this is the whole of its cost: a component pays
# for its own arena, its `rebuild_derived` and its unit table, all of which the record had already built
# once. Sodium benzoate, canonical order only: 1.89 us whole record against 4.26 us decomposed, next to
# 1.62 us for benzoic acid alone. Paid on every multi-component record whose refinement does not land
# discrete, in exchange for the column above; a single-atom component -- the counter-ion, and water once
# its hydrogens are implicit -- is special-cased below and pays none of it.
#
# SOUND BECAUSE AN AUTOMORPHISM CARRIES COMPONENTS ONTO COMPONENTS. No automorphism relates two atoms
# whose components have different certificates, so the pairs a whole-record search spends its budget
# refuting are exactly the ones this never asks about. The block order is a function of the multiset of
# component certificates, hence of the record; ties are isomorphic components, which is a freedom the
# paragraph above already grants.
#
# THE LABELLING IS NOT THE ONE A WHOLE-RECORD SEARCH PRODUCES, and no decomposition could be: the
# whole-record extremum INTERLEAVES components -- two hydrogen molecules score (2 << 8) at position 0
# interleaved against (1 << 8) blocked -- while these blocks are contiguous. The canonical FORM is
# still unique and every guarantee above still holds; what moved is which bytes a MIXTURE canonicalises
# to. Identity bytes and canonical SMILES of a multi-component record are not comparable across this
# change.
#
# THE BUDGET POLICY IS THE OPPOSITE OF mol_automorphisms', DELIBERATELY. Up there, a truncated
# search only prunes less: it can lose symmetry, the flag says so, and a consumer may legitimately
# fall back to doing the work the symmetry would have saved. Here a truncated search would return
# SOME labelling where THE labelling was asked for -- indistinguishable from a correct answer at
# the call site, and silently wrong in every hash and every equality built on top of it. So the
# node budget is hard, the failure is loud, and NO PARTIAL ORDER IS EVER WRITTEN to order_out.
# Do not soften this into a best-effort result for very symmetric records.
cdef enum:
    CANON_NODE_BUDGET = 1000000   # refinement tree nodes per call; see the paragraph above


# Everything the tree walk carries that does not change from node to node.
cdef struct canon_ctx_t:
    uint32_t n
    uint32_t *ptr
    halfedge_t *edges
    atom_t *atoms
    size_t graph_len              # 2 * atom_count + bond_count: the certificate's graph part
    size_t cert_len               # graph_len, plus atom_count when `digits` is live
    uint64_t *best                # the greatest certificate seen, cert_len words
    uint64_t *cert                # scratch for the candidate certificate
    uint64_t *nb                  # scratch for one atom's neighbour row, maxdeg words
    uint32_t *best_pos            # the labelling that produced `best`, atom_count entries
    uint32_t *inv                 # scratch: position -> slot
    uint32_t *digits              # NULL for a stereo-blind certificate; else n parity digits
    uint32_t *sscratch            # 2n uint32_t lent to `_canon_stereo_hook`
    uint32_t *pcls                # n uint32_t: the colouring the orbit search is seeded with
    uint8_t *unnamed              # n bytes: the atoms of the frames this colouring cannot name
    uint64_t nodes                # nodes entered so far, against `budget`
    uint64_t budget               # CANON_NODE_BUDGET, except through the seam; see _canon_order
    bint have_best
    bint asymmetric               # the root's group is trivial, so the extremum is unique
    bint exceeded


cdef void _canon_certificate(Structure structure, canon_ctx_t *ctx, uint32_t *pos,
                             uint64_t *out) noexcept nogil:
    """The labelled molecule as one word string, comparable lexicographically.

    `pos` is a DISCRETE colouring, so pos[s] - 1 is slot s's canonical position. The string is
    written in position order and mentions no slot number, which is what makes two labellings of
    the same molecule comparable at all:

        for each position p: the atom's invariant word, then how many of its bonds go to a HIGHER
        position, then those bonds as (position << 8 | bond colour), ascending.

    Every bond appears exactly once, at its lower-positioned end, so the graph part is
    2 * atom_count + bond_count words long for every labelling of a given molecule -- equal
    lengths, hence a total order under word-by-word comparison.

    THEN THE PARITY TAIL, one word per position, when `ctx.digits` is live: the digit
    `_canon_stereo_hook` reads at the slot that took that position, in the frame `pos` names.  A
    DISCRETE colouring names every frame -- every atom has its own colour, so no unit's directions can
    tie -- which is why the tail is meaningful exactly where it is written, at a leaf.

    WHY IT IS A TAIL AND NOT INTERLEAVED WITH THE ATOM WORDS: the graph part has to keep deciding
    first.  Two labellings of one molecule always tie on it -- they are the same graph -- so nothing is
    lost, while two molecules that differ constitutionally are separated without the parity ever being
    consulted, and the tail is then the tie-break it is meant to be rather than a term that could
    reorder a constitutional comparison.  It is also the layout `mol_certificate_words` already writes
    for its `extra` array, so the string the search maximises and the string `mol_identity_bytes`
    publishes are the same string, which is a property worth having structurally instead of by
    coincidence: they must agree, or the labelling the search chose would not be the labelling the
    published bytes were read in.

    `ctx.digits` NULL leaves the graph part alone and writes no tail, which is what
    `mol_certificate_words` wants -- there the caller supplies `extra` itself.

    The atom word is `_atom_invariant`, the same packing the refinement starts from, so the
    certificate cannot separate less than the refinement does. It carries one field
    `_atom_colour_equal` deliberately omits -- the ring-membership bit -- which is harmless
    because it is a function of the graph: it shifts the certificate of every labelling of a
    molecule by the same amount and so cannot change which one is greatest.

    The bond colour is the Kekule order with the aromatic bit below it, matching
    `_bond_colour_equal`. HE_IN_RING is derived, so it is left out on both sides.
    """
    cdef uint32_t n = ctx.n
    cdef uint32_t *ptr = ctx.ptr
    cdef halfedge_t *edges = ctx.edges
    cdef uint32_t *inv = ctx.inv
    cdef uint64_t *nb = ctx.nb
    cdef halfedge_t *e
    cdef uint32_t s, p, q, k, d
    cdef size_t w = 0

    for s in range(n):
        inv[pos[s] - 1] = s
    for p in range(n):
        s = inv[p]
        out[w] = _atom_invariant(&ctx.atoms[s])
        w += 1
        d = 0
        for k in range(ptr[s], ptr[s + 1]):
            e = &edges[k]
            q = pos[e.to] - 1
            if q > p:
                nb[d] = ((<uint64_t> q << 8) | (<uint64_t> e.order << 1)
                         | (<uint64_t> 1 if e.flags & HE_AROMATIC else <uint64_t> 0))
                d += 1
        _sort_words(nb, d)          # the CSR's neighbour order must not reach the string
        out[w] = <uint64_t> d
        w += 1
        for k in range(d):
            out[w] = nb[k]
            w += 1
    if ctx.digits is not NULL:
        _canon_stereo_hook(structure, pos, ctx.digits, ctx.sscratch, NULL)
        for s in range(n):
            out[w + pos[s] - 1] = <uint64_t> ctx.digits[s]


cdef inline int _canon_cert_cmp(uint64_t *a, uint64_t *b, size_t count) noexcept nogil:
    """Lexicographic comparison of two equal-length certificates: 1 if a > b, -1 if a < b, 0 if
    equal. Which extremum is taken is arbitrary as long as it is fixed; this file takes the
    maximum, and `_canon_branch` is the only caller."""
    cdef size_t i
    for i in range(count):
        if a[i] != b[i]:
            return 1 if a[i] > b[i] else -1
    return 0


cdef uint64_t _canon_indicator(canon_ctx_t *ctx, uint32_t *cls, uint32_t classes,
                               uint32_t *counts, uint64_t *sizes) noexcept nogil:
    """A relabelling-invariant fingerprint of one coloured graph, used to compare sibling nodes.

    Two parts, both of them functions of the colouring and the graph and of nothing else:

      * the CELL SIZES read out in class-id order. Class ids come from `_classify`, which sorts by
        key value, so the order they are read in is invariant even though it is not meaningful.
      * the multiset over atoms of one further refinement round -- each atom's own class hashed
        together with its sorted (neighbour class, bond order) list. Summed, because a sum over
        slots cannot leak the slot order, and because the quotient graph is what distinguishes
        colourings that happen to share a cell-size profile.

    The result is a hash, not an order-preserving encoding, which is fine: the pruning below needs
    an invariant to select the extremal siblings by, and a hash of an invariant is an invariant. A
    collision costs siblings that are explored when they need not have been -- and it costs it
    identically for every relabelling of the record, so it cannot cost invariance.

    `counts` (n + 1 entries) and `sizes` (n entries) are caller-owned scratch.
    """
    cdef uint32_t n = ctx.n
    cdef uint32_t *ptr = ctx.ptr
    cdef halfedge_t *edges = ctx.edges
    cdef uint64_t *nb = ctx.nb
    cdef uint32_t s, c, k, d
    cdef uint64_t acc = 0
    cdef uint64_t pair[2]

    memset(counts, 0, <size_t> (n + 1) * sizeof(uint32_t))
    for s in range(n):
        counts[cls[s]] += 1
    for c in range(classes):
        sizes[c] = <uint64_t> counts[c + 1]
    for s in range(n):
        d = 0
        for k in range(ptr[s], ptr[s + 1]):
            nb[d] = (<uint64_t> cls[edges[k].to] << 8) | <uint64_t> edges[k].order
            d += 1
        _sort_words(nb, d)
        acc += _xxh64(nb, d, <uint64_t> cls[s])
    pair[0] = _xxh64(sizes, classes, <uint64_t> classes)
    pair[1] = acc
    return _xxh64(pair, 2, 0)


cdef inline uint32_t *_canon_pinned(uint32_t n, uint32_t *pcls, uint8_t *unnamed) noexcept nogil:
    """Give every atom of an unnamed frame a colour no other atom carries, and hand `pcls` back.

    An automorphism preserving a singleton colour class fixes that atom, so the group of the pinned
    colouring fixes every unnamed frame pointwise and inverts none of their parities. `4 * n + 8 + i`
    clears the largest `cls[s] * 4 + digit` a node can hold, which is `4 * n + 3`.
    """
    cdef uint32_t i
    for i in range(n):
        if unnamed[i]:
            pcls[i] = 4 * n + 8 + i
    return pcls


cdef int _canon_branch(Structure structure, canon_ctx_t *ctx, uint32_t *cls, uint32_t classes,
                       bint may_be_symmetric) except -1:
    """One node of the refinement tree. `cls` is an equitable colouring with `classes` classes.

    A discrete `cls` is a leaf: its certificate is built and kept if it beats the best so far.
    Otherwise the node branches over the target cell -- the lowest class id with more than one
    member, which exists by pigeonhole whenever classes < n -- individualising one candidate at a
    time and recursing on the refinement of the result.

    TERMINATION. Individualising an atom of a cell of size two or more, then refining, gives a
    partition strictly finer than `cls`, so `classes` rises by at least one per level and the
    depth is at most atom_count. Nodes are counted against CANON_NODE_BUDGET; the depth is
    bounded by the atom count rather than by the budget, and it is C recursion, so a record with
    tens of thousands of atoms AND symmetry deep enough to branch at every level would be a stack
    problem before it was a budget problem. Nothing chemical comes close; an explicit stack is the
    fix if something ever does.

    TWO PRUNES, AND WHY NEITHER OF THEM PICKS A CANDIDATE BY SLOT ORDER. Taking the extremum over
    the target cell is the point of this search -- taking the first candidate that refines is the
    defect it exists to fix -- so every candidate dropped here has to be dropped for a reason that
    survives relabelling the record.

    (1) BY THE ORBITS OF THIS NODE, not of the molecule. Two candidates related by a symmetry of
    the CURRENT colouring have subtrees that are images of each other, so their leaves carry the
    same certificates and either one stands for both. Those symmetries are the automorphisms that
    preserve `cls`, which is what `mol_automorphisms` returns when `cls` is handed to it as the
    seed; the orbits of the whole molecule would only be valid to prune with at the root. The
    orbit ARRAY is what is read -- mol_automorphisms hands back no group elements, since the rows
    such a search yields are orbit-complete without being a generating set and nothing could be
    derived from them. Which member of an orbit is kept IS decided by slot order, and that is
    the one place it may be: the members are interchangeable by construction.

    "INTERCHANGEABLE" IS A CLAIM ABOUT THE CERTIFICATE, AND THE CERTIFICATE NOW CARRIES A PARITY, so
    the seed handed to `mol_automorphisms` is `cls` refined by the parity digits and not `cls`, and
    where no colouring can name a configured unit's frame this prune is switched off entirely. The
    reasoning is at the call site; the short version is that a graph symmetry which inverts a
    configuration relates two candidates whose leaves DO NOT carry the same certificate, so dropping
    one of them by slot order is exactly the defect the extremal search exists to prevent, moved one
    level down. With no stereo the refined colouring is `cls` scaled by four and this prune is
    unchanged, which is the common path and the one that must not pay.

    (2) BY THE NODE INVARIANT, keeping only the children whose `_canon_indicator` is greatest. This
    is an INVARIANT-NARROWED EXTREMUM, and it is NOT nauty's Lambda pruning even though the
    indicator plays the same role there: nauty prunes branches it can prove cannot hold the
    extremum, so its canonical form is the whole-cell one, while narrowing to the indicator-maximal
    children DOES change which labelling comes out canonical. It is what makes the search finish.
    Without it, 40 cycloalkane rings in one record -- three ring sizes, so three candidates survive
    the orbit prune at each of ~30 levels -- is a tree of 3**30 nodes, and 120 atoms of ordinary
    chemistry blew the node budget in 18 s. With it the same record is a path.

    Narrowing is legitimate for the same reason (1) is: the indicator is a function of the coloured
    graph, so "the children whose indicator is maximal" is an invariant subset of the cell,
    relabelling maps it onto itself, and an extremum over an invariant subset is still invariant.
    So the canonical form here is the greatest certificate among the leaves the indicator leaves
    reachable -- a different, equally canonical choice of representative than the greatest over all
    of them, and NOT interchangeable with it: this file's canonical form must not be compared
    against one produced by a whole-cell search. Any invariant selection would do; this one is
    cheap and strong.

    The two compose without argument, because automorphic candidates have equal indicators: the
    maximum over orbit representatives is the maximum over the whole cell, so (1) never hides the
    winner of (2).

    Pruning never changes which leaf wins among the leaves it leaves reachable, so the two ways
    the orbit prune can be unavailable are both safe: a truncated orbit search
    (CANON_BUDGET_EXCEEDED, whose partition must not be trusted) prunes nothing at that node, and
    a node whose parent proved its group trivial (CANON_ASYMMETRIC) skips the orbit search
    entirely, since a subgroup of the trivial group could prune nothing anyway. What either can
    change is how many nodes the call spends, hence whether it reaches CANON_NODE_BUDGET -- so a
    record that FAILS is not guaranteed to fail identically under a different slot order. A record
    that succeeds succeeds with the same certificate either way, which is the property that
    matters.
    """
    cdef uint32_t n = ctx.n
    cdef uint64_t *hashes
    cdef uint64_t *sizes
    cdef uint32_t *seed
    cdef uint32_t *child
    cdef uint32_t *orbits
    cdef uint32_t *counts
    cdef uint8_t *seen
    cdef uint8_t *alive
    cdef char *block
    cdef uint64_t best_ind = 0
    cdef size_t u64_len = <size_t> 2 * n * sizeof(uint64_t)
    cdef size_t u32_len = (4 * <size_t> n + 1) * sizeof(uint32_t)
    cdef uint32_t v, s, c, target = 0, flags = 0
    cdef Py_ssize_t child_classes
    cdef bint use_orbits = False, child_symmetric = False, first = True, all_named = True

    ctx.nodes += 1
    if ctx.nodes > ctx.budget:
        ctx.exceeded = True
        raise AutomorphismBudgetExceeded(
            'canonical labelling exceeded its budget of %d refinement nodes; a truncated extremal '
            'search returns some labelling instead of the canonical one, so there is no partial '
            'answer to return' % ctx.budget)
    if not (ctx.nodes & 0xFFF):
        # cpython.exc declares this `except -1`; that clause is what turns the error indicator the
        # C call leaves behind into a raised KeyboardInterrupt. Without it a long search defers
        # SIGINT until it ends.
        PyErr_CheckSignals()

    if classes == n:
        with nogil:
            _canon_certificate(structure, ctx, cls, ctx.cert)
        if not ctx.have_best or _canon_cert_cmp(ctx.cert, ctx.best, ctx.cert_len) > 0:
            memcpy(ctx.best, ctx.cert, ctx.cert_len * sizeof(uint64_t))
            for s in range(n):
                ctx.best_pos[s] = cls[s] - 1
            ctx.have_best = True
        return 0

    # One block per node. u64 spans first so the u32 spans behind them stay aligned: the per
    # candidate indicator, and the cell-size buffer `_canon_indicator` writes through. Then the
    # individualising seed, the child colouring it refines to, this node's orbits, and the
    # per-class member counts that locate the target cell. Then two byte maps -- which orbits have
    # been used, indexed by orbit id, and which candidates survived the orbit prune, indexed by
    # slot. Class and orbit ids run 1 .. n, so those two spans need n + 1 entries.
    block = <char *> PyMem_Malloc(u64_len + u32_len + <size_t> 2 * n + 2)
    if block is NULL:
        raise MemoryError('canonical branch scratch allocation failed')
    hashes = <uint64_t *> block
    sizes = hashes + n
    seed = <uint32_t *> (block + u64_len)
    child = seed + n
    orbits = child + n
    counts = orbits + n
    seen = <uint8_t *> (block + u64_len + u32_len)
    alive = seen + n + 1
    try:
        memset(counts, 0, <size_t> (n + 1) * sizeof(uint32_t))
        for s in range(n):
            counts[cls[s]] += 1
        for c in range(1, classes + 1):
            if counts[c] > 1:
                target = c
                break
        if may_be_symmetric:
            if ctx.digits is NULL:
                mol_automorphisms(structure, cls, orbits, &flags)
                use_orbits = not (flags & CANON_BUDGET_EXCEEDED)
                child_symmetric = not (flags & CANON_ASYMMETRIC)
            else:
                # THE ORBIT PRUNE UNDER A CONFIGURATION, in two steps, and the second is the whole
                # fix.  Step one: the orbits are those of the colouring REFINED BY THE PARITY DIGITS
                # and not of `cls` -- so a symmetry the search may prune with has to carry the
                # configuration along, not just the graph.  That is sound where the digits are
                # meaningful and free where there is no stereo: every digit is then 0, `cls[s] * 4`
                # has exactly `cls`'s equality classes, and the group is the same group.
                with nogil:
                    all_named = _canon_stereo_hook(structure, cls, ctx.digits, ctx.sscratch,
                                                   ctx.unnamed)
                    for s in range(n):
                        ctx.pcls[s] = cls[s] * 4 + ctx.digits[s]
                # Step two: WHERE THIS COLOURING CANNOT NAME A CONFIGURED UNIT'S FRAME, DO NOT PRUNE
                # AT ALL.  An unnamed frame means two of the unit's directions share a colour, so the
                # digit is the "configured, value unreadable" code and the refined colouring is blind
                # to the very parity the automorphism might invert -- which is precisely the mirror
                # case, because the two ring branches leaving a carbinol carbon share a colour.  There
                # the prune would drop a labelling it cannot prove is stereo-equivalent, and it drops
                # it by slot order.  Both labellings reach a leaf instead, where the colouring is
                # discrete, every frame is named, and the certificate's parity tail separates them as
                # a function of the molecule.
                #
                # Correctness does not depend on this test being TIGHT, only on it never claiming
                # "named" when a frame is not: a false "unnamed" costs nodes and pruning, a false
                # "named" would cost the answer.  It is deliberately loose in one further way -- the
                # hook also reports unnamed when two units land a digit on one atom, where `max`
                # merges the two codes and neither is recoverable.
                #
                # AND "DO NOT PRUNE AT ALL" MEANS PRUNE WITH A SMALLER GROUP, NOT WITH NONE.  An
                # unnamed frame bars the symmetries that permute ITS OWN two shared-colour
                # directions, and nothing else: a symmetry fixing every atom the hook marks carries
                # each unnamed unit onto itself with its directions in place, inverts nothing, and is
                # entitled to prune.  Pinning those atoms to singleton colours before the search
                # selects exactly that subgroup, and pruning by ANY subgroup of the stabiliser is
                # sound -- the candidates it relates still have equal certificates, so which of two
                # tying leaves wins is free, the same freedom the target cell's own extremum has.
                #
                # Testing the group's support instead -- "are the marked atoms in singleton orbits" --
                # reads as the cheaper form of this and is very nearly vacuous: the marked atoms
                # INCLUDE the unit's two shared-colour directions, and the symmetry that exchanges
                # them is the reason the frame is unnamed, so the marked atoms are in a non-singleton
                # orbit whenever the mirror case is present, which is when the question is asked.
                #
                # Measured, a resin-bound protected peptide of 235 atoms carrying trityl groups spent
                # the whole 1,000,000-node budget with the prune standing down per node, because a
                # node's ~51 genuine 3-orbits of equivalent phenyls and methyls branch three ways each
                # and `_canon_indicator` scores automorphic candidates equally by construction. Those
                # orbits do not touch an unnamed frame, so pinning leaves them intact and the tree is a
                # path again. What pinning does give up is a symmetry carrying one unnamed unit onto
                # ANOTHER -- n copies of one symmetric component, where the copy swap is safe but is
                # not a subgroup element -- which stays a refusal, since proving that swap
                # parity-preserving needs the group elements the search does not hand back.
                mol_automorphisms(structure, ctx.pcls if all_named else
                                  _canon_pinned(n, ctx.pcls, ctx.unnamed), orbits, &flags)
                use_orbits = not (flags & CANON_BUDGET_EXCEEDED)
                # A subgroup being trivial says nothing about a child's, whose marked set may be
                # smaller, so only the unpinned group may answer for the subtree or for the root.
                child_symmetric = not all_named or not (flags & CANON_ASYMMETRIC)
            if ctx.nodes == 1 and not child_symmetric:
                ctx.asymmetric = True       # the root's own group, the only one worth reporting
        memset(seen, 0, <size_t> n + 1)
        memset(alive, 0, <size_t> n)

        # Pass one: refine every candidate the orbit prune keeps, and score it. `counts` is
        # scratch again from here -- the target cell has already been chosen from it.
        for v in range(n):
            if cls[v] != target:
                continue
            if use_orbits:
                if seen[orbits[v]]:
                    continue                # a symmetry of this node already tried this candidate
                seen[orbits[v]] = 1
            memcpy(seed, cls, <size_t> n * sizeof(uint32_t))
            # A label no other atom carries. Only equality classes of the seed reach the
            # partition, and the VALUE only reaches which class id the refinement gives which
            # cell -- and it is derived from the parent colouring, so it corresponds under any
            # relabelling just as the colouring does.
            seed[v] = classes + 1
            with nogil:
                child_classes = compute_atoms_order(structure, child, seed)
            if child_classes < 0:
                raise MemoryError('atom order refinement failed to allocate')
            alive[v] = 1
            with nogil:
                hashes[v] = _canon_indicator(ctx, child, <uint32_t> child_classes, counts, sizes)
            if first or hashes[v] > best_ind:
                best_ind = hashes[v]
                first = False

        # Pass two: recurse into the extremal children only. The refinement runs a second time
        # rather than being kept from pass one, which would cost a colouring per candidate held
        # live across the whole subtree walk below it.
        for v in range(n):
            if not alive[v] or hashes[v] != best_ind:
                continue
            memcpy(seed, cls, <size_t> n * sizeof(uint32_t))
            seed[v] = classes + 1
            with nogil:
                child_classes = compute_atoms_order(structure, child, seed)
            if child_classes < 0:
                raise MemoryError('atom order refinement failed to allocate')
            _canon_branch(structure, ctx, child, <uint32_t> child_classes, child_symmetric)
    finally:
        PyMem_Free(block)
    return 0


cdef inline int mol_canonical_order(Structure structure, uint32_t *seed, uint32_t *order_out,
                                    uint32_t *flags_out, bint stereo) except -1:
    """The molecule's canonical atom order, with the shipped node budget. See `_canon_order`, which
    is the same call with the budget spelled out; this is the entry point every consumer wants."""
    return _canon_order(structure, seed, order_out, flags_out, CANON_NODE_BUDGET, stereo)


cdef int _canon_order(Structure structure, uint32_t *seed, uint32_t *order_out,
                      uint32_t *flags_out, uint32_t budget, bint stereo) except -1:
    """The molecule's canonical atom order: the extremal labelling of the refinement tree.

    `order_out` is a CALLER-SUPPLIED array of atom_count uint32_t. Entry i receives atom slot i's
    0-based canonical position, so it is a permutation of 0 .. atom_count - 1. Two molecules that
    are isomorphic as coloured graphs get labellings that induce the SAME relabelled graph; see
    the fragment comment for what that does and does not pin down, and for why an isomorphic pair
    can still differ in which of two interchangeable atoms took which of their two positions.

    `seed` is forwarded to compute_atoms_order exactly as in mol_automorphisms: NULL starts from
    the atom records, otherwise it is per-atom starting labels, and any distinction it makes is
    respected by the order. That is how a caller pins down a labelling further than the structure
    alone can.

    flags_out[0] receives CANON_ASYMMETRIC when the molecule's group is known to be trivial --
    then the extremal labelling is the only one and no tie was broken anywhere.

    ON FAILURE THERE IS NO RESULT. If the tree exceeds `budget` nodes this sets
    CANON_BUDGET_EXCEEDED in flags_out[0], writes NOTHING to order_out, and returns -1 with
    AutomorphismBudgetExceeded raised. A labelling from a truncated extremal search is not an
    approximation of the canonical one, it is a different one, and it would corrupt every hash and
    every equality derived from it without any symptom. Contrast mol_automorphisms, where
    truncation only costs pruning and the caller is allowed to carry on with less symmetry.

    `stereo` asks for the parity fold below -- the stereo-aware START, as opposed to the parity tail
    on the leaf certificate, which every call gets and which is not optional. It costs nothing and
    changes nothing when `seed` is non-NULL, so the four call sites divide as follows and this is the
    only place the division is written down:

      * `MoleculeContainer.canonical_order`, seed NULL, stereo TRUE -- so it runs the same search
        `mol_identity_bytes` runs. The two must not disagree about WHAT the canonical order is.
      * `mol_identity_bytes` and `canonical_stereo_group_ids`, seed non-NULL, stereo TRUE -- the fold
        is already in their seed (this fold IS `mol_identity_bytes`' round zero), so the flag is a
        statement of intent and not a behaviour.
      * `smw_canonical_positions`, seed non-NULL under `o.stereo` and NULL without it, stereo passed
        straight through. FALSE IS LOAD-BEARING THERE: `format(mol, '!s')` is a constitution key, and
        a folded start makes the string a function of the configuration. Measured on the 393-record
        corpus of `test/`: folding it moved 63 strings, and all 63 stopped surviving a
        write-as-`!s`-then-read-then-write-as-`!s` round trip.

    `budget` is a parameter and not the constant inline BECAUSE OF THAT PARAGRAPH: the guarantee it
    states is only checkable by a test that can reach the failure, and no record small enough for a
    test comes anywhere near 1,000,000 nodes. Every real caller goes through
    `mol_canonical_order`, which passes CANON_NODE_BUDGET; the only other caller is the
    `_node_budget` seam on `MoleculeContainer.canonical_order`, which exists to test this path and
    nothing else. It is not a tuning knob -- a smaller budget does not buy a faster answer, it buys
    an exception.
    """
    cdef uint32_t n = structure.header.atom_count
    cdef canon_ctx_t ctx
    cdef char *block = NULL
    cdef uint32_t *cls = NULL
    cdef uint32_t *comp = NULL
    cdef uint32_t i, d, maxdeg = 0
    cdef Py_ssize_t classes, k
    cdef size_t u64_len

    flags_out[0] = 0
    if n == 0:
        return 0

    # BEFORE ANY POINTER IS TAKEN (ruling F60).  The stereo hook reads the unit table and the table is
    # built by appending a segment, which reallocates the arena -- so the one call that can move the
    # buffer happens here, at the top, and never from inside the search.  Idempotent and free after
    # the first call on a molecule; NULL when `_stereo.pxi` has not installed a hook, which is the
    # stereo-blind configuration.
    if _canon_prepare_hook is not NULL:
        _canon_prepare_hook(structure)

    # Read-only from here on: nothing in this fragment appends to the arena, so caching the
    # segment pointers cannot outlive a reallocation.
    ctx.n = n
    ctx.ptr = csr_ptr(structure)
    ctx.edges = csr_edges(structure)
    ctx.atoms = structure.atoms()
    ctx.nodes = 0
    ctx.budget = budget
    ctx.have_best = False
    ctx.asymmetric = False
    ctx.exceeded = False

    # Two allocations rather than one, against the usual rule, because they have different
    # lifetimes AND different odds: every call needs the colouring, while only a call that
    # actually searches needs the certificate machinery -- and on real molecules the refinement is
    # discrete, so the common path must not pay for a block it will not read.
    #
    # The stereo scratch moved into THIS block from the second one so that the parity fold below can
    # reach it before the certificate machinery exists -- the fold's whole purpose is to reach a
    # discrete colouring and return without allocating that machinery at all. 4n uint32_t on a call
    # that has already committed to a non-discrete refinement, and nothing on a stereo-blind build.
    try:
        # The byte map rides at the end of the word spans, so the words stay aligned by construction.
        if _canon_stereo_hook is not NULL:
            cls = <uint32_t *> PyMem_Malloc(<size_t> 5 * n * sizeof(uint32_t) + <size_t> n)
        else:
            cls = <uint32_t *> PyMem_Malloc(<size_t> n * sizeof(uint32_t))
        if cls is NULL:
            raise MemoryError('canonical order scratch allocation failed')
        if _canon_stereo_hook is not NULL:
            ctx.digits = cls + n
            ctx.sscratch = ctx.digits + n       # 2n
            ctx.pcls = ctx.sscratch + 2 * n
            ctx.unnamed = <uint8_t *> (ctx.pcls + n)
        else:
            ctx.digits = NULL
            ctx.sscratch = NULL
            ctx.pcls = NULL
            ctx.unnamed = NULL

        with nogil:
            classes = compute_atoms_order(structure, cls, seed)
        if classes < 0:
            raise MemoryError('atom order refinement failed to allocate')
        if classes == <Py_ssize_t> n:
            # Discrete already: an automorphism preserves the refinement, so the group is trivial
            # and this labelling is the only one the tree contains. One refinement, no search.
            flags_out[0] = CANON_ASYMMETRIC
            for i in range(n):
                order_out[i] = cls[i] - 1
            return 0

        # THE PARITY FOLD, AND WHY IT IS HERE RATHER THAN AT A CALL SITE.  The search refines by
        # parity at every node, so it reaches a stereo-distinguishing labelling from a bare
        # constitutional colouring -- but it reaches it by BRANCHING, and the branch it takes is
        # between candidates the constitution ties and the parity separates.  Where the orbit prune
        # stands down (it must, wherever a frame is unnamed) nothing else prunes them either: the
        # parity is a SUFFIX of the leaf certificate, so `_canon_indicator`, which is constitutional,
        # scores both candidates equally and both subtrees are walked in full.  Measured on a
        # 46-atom polychlorinated skeleton with twelve such pairs: 7413 tree nodes and 15592
        # refinements against 47 and 108 for the same molecule entered with the fold, a 100x gap.
        #
        # Folding the digits into the ROOT partition is exactly what `mol_identity_bytes` already did
        # for itself -- the same `_frame_free_parity_seed`, read in the frame the same colouring
        # names, folded the same way (`cls * 4 + digit`) -- so this is that optimisation moved down
        # one level to where every caller gets it, not a new mechanism.  It is sigma-equivariant by
        # ruling F95, so the order it produces is still a function of the molecule alone.
        #
        # ONLY WHEN THE CALLER PASSED NO SEED.  A caller who supplied one has already folded whatever
        # it wanted the refinement to start from (`smw_stereo_seed` runs the fixpoint, of which this
        # is round zero), and refining a supplied seed AGAIN would move that caller's answer.  And
        # only under `stereo`, for the reason the docstring's third bullet measures.
        if _canon_stereo_hook is not NULL and seed is NULL and stereo:
            with nogil:
                _canon_stereo_hook(structure, cls, ctx.digits, ctx.sscratch, NULL)
                for i in range(n):
                    ctx.pcls[i] = cls[i] * 4 + ctx.digits[i]
                classes = compute_atoms_order(structure, cls, ctx.pcls)
            if classes < 0:
                raise MemoryError('atom order refinement failed to allocate')
            if classes == <Py_ssize_t> n:
                flags_out[0] = CANON_ASYMMETRIC
                for i in range(n):
                    order_out[i] = cls[i] - 1
                return 0

        # THE DECOMPOSITION, AND WHY IT IS HERE AND NOT AT THE TOP. Everything above is O(n) and
        # answers most records without a search; a record whose refinement lands discrete needs no
        # component labelling at all, and a salt that reaches this line pays for one pass over its CSR.
        # See the fragment comment for what the split buys, what it costs and what it moves.
        comp = <uint32_t *> PyMem_Malloc(<size_t> n * sizeof(uint32_t))
        if comp is NULL:
            raise MemoryError('canonical order component labelling allocation failed')
        try:
            with nogil:
                k = label_components(structure, comp)
            if k < 0:
                raise MemoryError('component labelling failed to allocate')
            if k > 1:
                return _canon_order_split(structure, seed, order_out, flags_out, budget, stereo,
                                          comp, <uint32_t> k)
        finally:
            PyMem_Free(comp)

        for i in range(n):
            d = ctx.ptr[i + 1] - ctx.ptr[i]
            if d > maxdeg:
                maxdeg = d
        ctx.graph_len = <size_t> 2 * n + structure.header.bond_count
        # The parity tail costs n words of certificate here; its 4n of scratch is carved off `cls`
        # above, where the fold can reach it. Only a molecule that actually reaches the search pays
        # for this block -- the two discrete-refinement paths above return first. A NULL hook keeps
        # the old sizes exactly, so nothing here is a cost the stereo-blind build would not have had.
        if _canon_stereo_hook is not NULL:
            ctx.cert_len = ctx.graph_len + n
        else:
            ctx.cert_len = ctx.graph_len
        u64_len = (2 * ctx.cert_len + maxdeg) * sizeof(uint64_t)
        # u64 spans first so the u32 spans behind them stay aligned
        block = <char *> PyMem_Malloc(u64_len + <size_t> 2 * n * sizeof(uint32_t))
        if block is NULL:
            raise MemoryError('canonical order scratch allocation failed')
        ctx.best = <uint64_t *> block
        ctx.cert = ctx.best + ctx.cert_len
        ctx.nb = ctx.cert + ctx.cert_len
        ctx.best_pos = <uint32_t *> (block + u64_len)
        ctx.inv = ctx.best_pos + n

        _canon_branch(structure, &ctx, cls, <uint32_t> classes, True)
        if not ctx.have_best:
            # Unreachable: classes < n gives the root a target cell of at least two candidates,
            # each of which refines to a strictly finer partition, so the recursion reaches a
            # discrete colouring on every path. Checked anyway rather than handing back an
            # uninitialised order, which is the one failure this whole fragment exists to avoid.
            raise RuntimeError('canonical labelling produced no leaf')
        if ctx.asymmetric:
            flags_out[0] |= CANON_ASYMMETRIC
        memcpy(order_out, ctx.best_pos, <size_t> n * sizeof(uint32_t))
    finally:
        if ctx.exceeded:
            flags_out[0] |= CANON_BUDGET_EXCEEDED
        PyMem_Free(block)
        PyMem_Free(cls)
    return 0


cdef inline int _canon_key_cmp(const uint64_t *a, uint32_t alen,
                               const uint64_t *b, uint32_t blen) noexcept nogil:
    """Order two component keys: LENGTH FIRST, then the words. 1 if `a` sorts above `b`, -1 below.

    Length before content because two keys of different lengths are not prefixes to be compared word by
    word -- they describe components of different sizes, and the bigger component leading is the same
    convention `_canon_certificate` follows within a component.
    """
    cdef uint32_t i
    if alen != blen:
        return 1 if alen > blen else -1
    for i in range(alen):
        if a[i] != b[i]:
            return 1 if a[i] > b[i] else -1
    return 0


cdef void _canon_key_sort(uint32_t k, const uint64_t *keys, const uint32_t *koff,
                          uint32_t *idx, uint32_t *tmp) noexcept nogil:
    """Fill `idx[0:k]` with component ids ordered by key DESCENDING, stably.

    STABILITY IS WHAT MAKES THE TIE RULE STATABLE: equal keys mean isomorphic components, so which of
    them takes which block is a freedom the canonical labelling already has, and keeping arena order
    (components are numbered by lowest member) makes the choice deterministic rather than merely
    unspecified. Merge sort rather than an insertion sort because a plate of a thousand identical
    components is the case this whole decomposition exists for, and there every comparison runs the
    full key.
    """
    cdef uint32_t width = 1, lo, mid, hi, i, j, o
    for i in range(k):
        idx[i] = i
    while width < k:
        lo = 0
        while lo < k:
            mid = lo + width
            if mid > k:
                mid = k
            hi = mid + width
            if hi > k:
                hi = k
            i = lo
            j = mid
            o = lo
            while i < mid or j < hi:
                if j >= hi:
                    tmp[o] = idx[i]
                    i += 1
                elif i >= mid:
                    tmp[o] = idx[j]
                    j += 1
                elif _canon_key_cmp(keys + koff[idx[i]], koff[idx[i] + 1] - koff[idx[i]],
                                    keys + koff[idx[j]], koff[idx[j] + 1] - koff[idx[j]]) >= 0:
                    tmp[o] = idx[i]
                    i += 1
                else:
                    tmp[o] = idx[j]
                    j += 1
                o += 1
            lo = hi
        for i in range(k):
            idx[i] = tmp[i]
        width *= 2


cdef int _canon_order_split(Structure structure, uint32_t *seed, uint32_t *order_out,
                            uint32_t *flags_out, uint32_t budget, bint stereo,
                            const uint32_t *comp, uint32_t k) except -1:
    """`_canon_order` for a record of `k` > 1 connected components: each component's own order, blocked.

    `comp` is `label_components`' output. Contract identical to `_canon_order`'s -- `order_out` receives
    a permutation of 0 .. atom_count - 1, `flags_out[0]` gets CANON_ASYMMETRIC only when the whole
    record's group is trivial, and a truncated search writes nothing and raises. The fragment comment
    above holds the soundness argument and the one thing this moves.

    THE KEY A COMPONENT IS ORDERED BY is its own certificate (`mol_certificate_words`), with the parity
    tail included exactly when `stereo` -- so `format(mol, '!s')` stays a constitution key and keeps
    surviving its write-read-write round trip, which is the third bullet of `_canon_order`'s docstring
    applied to the block order. When the caller passed a `seed`, its values in canonical-position order
    are appended, because `canonical_order(seed=...)` promises that every distinction the seed makes is
    respected -- and a seed that separates two otherwise isomorphic components has to separate their
    blocks too, or the promise holds inside a component and breaks between them.

    EACH COMPONENT GETS THE FULL NODE BUDGET. A budget bounds ONE extremal search and each component is
    its own search; sharing one across components would make a component's labelling depend on how many
    components happened to be searched before it, which is the dependency `CANON_MAX_NODES_SEARCH`'s own
    comment refuses for the same reason.
    """
    cdef uint32_t n = structure.header.atom_count
    cdef uint32_t *sptr = csr_ptr(structure)
    cdef atom_t *atoms = structure.atoms()
    cdef uint32_t *block = NULL
    cdef uint64_t *keys = NULL
    cdef uint32_t *start
    cdef uint32_t *koff
    cdef uint32_t *idx
    cdef uint32_t *tmp
    cdef uint32_t *slots
    cdef uint32_t *lpos
    cdef uint32_t *local
    cdef uint32_t *sseed
    cdef uint32_t *pos1
    cdef uint32_t *digits
    cdef uint32_t *sscratch
    cdef uint32_t i, j, c, m, off, base, bonds, w
    cdef uint32_t sflags = 0
    cdef uint32_t asym = CANON_ASYMMETRIC
    cdef bint tail = _canon_stereo_hook is not NULL and stereo
    cdef Structure sub

    flags_out[0] = 0
    # One block for the component index (k + 1 offsets, k + 1 key offsets, k sort slots, k merge
    # scratch) and the per-component spans, all uint32_t: the slot lists and their local positions, the
    # inverse map `structure_component_graph` fills, the restricted seed, and the parity hook's 1-based
    # colouring, digits and 2m of scratch. Every span is sized at n rather than at the largest
    # component, which costs one array and removes a pass to find that size.
    block = <uint32_t *> PyMem_Malloc((<size_t> 4 * k + 2 + <size_t> 8 * n) * sizeof(uint32_t))
    if block is NULL:
        raise MemoryError('canonical order component scratch allocation failed')
    # 2n + bonds words of certificate per component, n of parity tail, n of seed -- so the whole record
    # fits in 4n + bond_count however it is divided.
    keys = <uint64_t *> PyMem_Malloc((<size_t> 4 * n + structure.header.bond_count)
                                    * sizeof(uint64_t))
    if keys is NULL:
        PyMem_Free(block)
        raise MemoryError('canonical order component key allocation failed')
    try:
        start = block
        koff = start + k + 1
        idx = koff + k + 1
        tmp = idx + k
        slots = tmp + k
        lpos = slots + n
        local = lpos + n
        sseed = local + n
        pos1 = sseed + n
        digits = pos1 + n
        sscratch = digits + n                       # 2n
        # Counting sort of the slots by component label. Ascending within a block, which is
        # `structure_component_graph`'s precondition: a monotone renumbering preserves every parity
        # frame, so no configuration has to be rewritten.
        for c in range(k + 1):
            start[c] = 0
        for i in range(n):
            start[comp[i] + 1] += 1
        for c in range(k):
            start[c + 1] += start[c]
        for c in range(k):
            koff[c] = start[c]                      # borrowed as the fill cursor
        for i in range(n):
            slots[koff[comp[i]]] = i
            koff[comp[i]] += 1

        koff[0] = 0
        for c in range(k):
            off = start[c]
            m = start[c + 1] - off
            if m == 1:
                # A LONE ION NEEDS NO SEARCH: its order is [0] and its key is the two words
                # `_canon_certificate` would write for it -- the atom's invariant, then a bond count of
                # zero -- with a zero parity digit under `tail`, a single atom having no stereo unit.
                # Worth the branch rather than the general path because it is the common component of a
                # salt: the counter-ion, and water once its hydrogens are implicit. It saves that
                # component its sub-structure, its `rebuild_derived` and its unit table.
                lpos[off] = 0
                keys[koff[c]] = _atom_invariant(&atoms[slots[off]])
                keys[koff[c] + 1] = 0
                w = 2
                if tail:
                    keys[koff[c] + 2] = 0
                    w = 3
                if seed is not NULL:
                    keys[koff[c] + w] = <uint64_t> seed[slots[off]]
                    w += 1
                koff[c + 1] = koff[c] + w
                continue
            bonds = 0
            for i in range(m):
                bonds += sptr[slots[off + i] + 1] - sptr[slots[off + i]]
            bonds //= 2
            sub = structure_component_graph(structure, slots + off, m, bonds, local)
            rebuild_derived(sub)
            if seed is not NULL:
                for i in range(m):
                    sseed[i] = seed[slots[off + i]]
            try:
                _canon_order(sub, sseed if seed is not NULL else NULL, lpos + off, &sflags,
                             budget, stereo)
            except:
                # The sub wrote its own flags before raising, and CANON_BUDGET_EXCEEDED is the one bit
                # the caller is promised to see beside the exception.
                flags_out[0] |= sflags & CANON_BUDGET_EXCEEDED
                raise
            asym &= sflags
            flags_out[0] |= sflags & CANON_BUDGET_EXCEEDED
            if tail:
                for i in range(m):
                    pos1[i] = lpos[off + i] + 1
                with nogil:
                    _canon_stereo_hook(sub, pos1, digits, sscratch, NULL)
                mol_certificate_words(sub, lpos + off, digits, keys + koff[c])
            else:
                mol_certificate_words(sub, lpos + off, NULL, keys + koff[c])
            w = <uint32_t> mol_certificate_len(sub, tail)
            if seed is not NULL:
                for i in range(m):
                    keys[koff[c] + w + lpos[off + i]] = <uint64_t> sseed[i]
                w += m
            koff[c + 1] = koff[c] + w

        with nogil:
            _canon_key_sort(k, keys, koff, idx, tmp)
            base = 0
            for j in range(k):
                c = idx[j]
                for i in range(start[c], start[c + 1]):
                    order_out[slots[i]] = base + lpos[i]
                base += start[c + 1] - start[c]
            # TWO EQUAL KEYS ARE A NON-TRIVIAL AUTOMORPHISM -- the one that swaps the two isomorphic
            # components -- so the record is asymmetric only when every component is and no two of them
            # are alike. Adjacent is enough: the sort put equal keys together.
            for j in range(1, k):
                if _canon_key_cmp(keys + koff[idx[j]], koff[idx[j] + 1] - koff[idx[j]],
                                  keys + koff[idx[j - 1]], koff[idx[j - 1] + 1] - koff[idx[j - 1]]) == 0:
                    asym = 0
                    break
        flags_out[0] |= asym
    finally:
        PyMem_Free(keys)
        PyMem_Free(block)
    return 0


cdef inline size_t mol_certificate_len(Structure structure, bint with_extra) noexcept nogil:
    """How many words `mol_certificate_words` writes: 2n + bonds, plus n when `extra` is supplied."""
    cdef size_t n = structure.header.atom_count
    return 2 * n + structure.header.bond_count + (n if with_extra else 0)


cdef int mol_certificate_words(Structure structure, uint32_t *order, uint32_t *extra,
                               uint64_t *out) except -1:
    """The canonical form of an ALREADY-LABELLED molecule, as one comparable word string.

    `order[s]` is slot s's 0-based canonical position, exactly as `mol_canonical_order` writes it;
    `out` is caller-owned and `mol_certificate_len` words long.  `extra`, when not NULL, is one
    per-SLOT word appended in POSITION order after the graph part -- the door through which a
    caller folds in state the graph does not carry.  Its only user is the stereo term in
    `mol_identity_bytes`; it is a parameter rather than a hardcoded read because `_canonical.pxi`
    is included before `_stereo.pxi` and must not know what a parity is.

    This is the piece `_canon_order` computes internally and throws away.  It is separated out
    because the certificate, not the labelling, is the thing an equality or a hash may rest on:
    the labelling is unique only up to the automorphism group (see the fragment comment), while
    the string below is read off positions and mentions no slot, so two labellings of one molecule
    produce the same words.  `_canon_order` cannot simply return its `best` either -- on a molecule
    whose refinement lands discrete it never builds one, which is the common case.

    NOT A SCREEN, AND NOT `signature`.  `MoleculeContainer.signature` is the OR of every atom's
    feature words: propane, butane and pentane share one, which is correct for a prefilter and
    catastrophic for an equality.  This string separates them, because it carries each position's
    own atom word and its own bond row rather than a disjunction over the molecule.
    """
    cdef uint32_t n = structure.header.atom_count
    cdef canon_ctx_t ctx
    cdef char *block = NULL
    cdef uint32_t *pos = NULL
    cdef uint32_t i, d, maxdeg = 0
    cdef size_t w
    if n == 0:
        return 0
    ctx.n = n
    ctx.ptr = csr_ptr(structure)
    ctx.edges = csr_edges(structure)
    ctx.atoms = structure.atoms()
    for i in range(n):
        d = ctx.ptr[i + 1] - ctx.ptr[i]
        if d > maxdeg:
            maxdeg = d
    # One block: the neighbour scratch `_canon_certificate` sorts into, the position -> slot
    # inverse it fills, and the 1-based colouring it expects in place of the 0-based order.
    block = <char *> PyMem_Malloc(<size_t> maxdeg * sizeof(uint64_t)
                                 + <size_t> 2 * n * sizeof(uint32_t))
    if block is NULL:
        raise MemoryError('canonical certificate scratch allocation failed')
    try:
        ctx.nb = <uint64_t *> block
        ctx.inv = <uint32_t *> (block + <size_t> maxdeg * sizeof(uint64_t))
        pos = ctx.inv + n
        # NULL, so `_canon_certificate` writes the graph part and stops: here the tail is the caller's
        # `extra`, filled below.  The two are the same words -- `mol_identity_bytes` passes the digits
        # the same hook produces -- but they arrive by different routes and this one must not double-
        # write them.
        ctx.digits = NULL
        ctx.sscratch = NULL
        for i in range(n):
            pos[i] = order[i] + 1
        _canon_certificate(structure, &ctx, pos, out)
        if extra is not NULL:
            w = mol_certificate_len(structure, False)
            for i in range(n):
                out[w + order[i]] = <uint64_t> extra[i]
    finally:
        PyMem_Free(block)
    return 0
