# -*- coding: utf-8 -*-
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
# GRAPH DESCRIPTORS: ARITHMETIC, AND NO PERCEPTION.  Every number in this file is computed from
# something the arena already holds -- the element index, the stored hybridization, the minimum cycle
# basis, the per-atom ring descriptors, and the one distance matrix in `_molecule_topology.pxi`.
# Nothing here walks for rings and nothing here runs a second BFS.  A descriptor that needed either
# would be a perception change wearing a descriptor's name, and the place to argue for it is the
# layer that owns the perception.
#
# THE GRAPH IS THE GRAPH AS STORED, ORDER 8 INCLUDED.  Degree is the CSR row length -- what
# `degree_of()` reports -- so a dative bond is an edge and an explicit hydrogen atom is a vertex.
# That is the same graph `distance_matrix`, `bond_count` and `connected_components` describe, and
# keeping one graph is what lets a distance-derived index and a degree-derived index be compared on
# one molecule.  It parts company with the SMARTS `D`, `x` and `z` primitives, which count
# SUBSTITUENTS and skip order 8: both counts are right and neither is a stale copy of the other.
#
# The consequence a caller has to know, and every docstring below repeats it: the classical indices
# are defined on the hydrogen-suppressed graph, which is exactly what chython stores when hydrogens
# are implicit.  A record carrying explicit hydrogen atoms is a different graph and gets different
# numbers.  Nothing here suppresses them on the caller's behalf.
#
# WHAT A DISCONNECTED MOLECULE GETS is settled per descriptor in the docstrings, never by accident:
# `distance_matrix` reports -1 for a pair in different components and nothing in this file sums a -1.
# Sums over vertices or bonds are additive over components and say nothing; the distance-derived ones
# say what they skip; `balaban_j` is the one refusal in the file, because Balaban's J is defined on a
# connected graph and its vertex distance sums are infinite otherwise.
#
# INCLUDED AFTER `_molecule_container.pxi`, `_molecule_topology.pxi` and `_fingerprints.pxi`: the
# wrappers in the class call forward into these `cdef` functions (RULES.md §7.1 says that direction is
# free), `mol_distance_matrix` and the `_NP_*` names are bound in the topology fragment, and two of
# these functions take a `MoleculeContainer` because they need `atoms_order` or have to raise.


# --- composition ----------------------------------------------------------------------------------

cdef inline uint32_t desc_carbon_count(Structure structure) noexcept nogil:
    """How many carbon atoms.  One subtraction on the element index; see `atoms_of_element`."""
    return element_bucket_end(structure, 6) - element_bucket_begin(structure, 6)


cdef uint32_t desc_carbon_sp3_count(Structure structure) noexcept nogil:
    """How many carbons store hybridization 1.

    z == 1 AND NOTHING ELSE.  The z scale is 1..6 and reports what it found rather than saturating,
    and 1 is sp3.  An aromatic carbon is 4, an allene's central carbon is 5, and neither is sp3.

    Walks the carbon bucket rather than every atom, so a metal complex with two carbons costs two
    iterations.
    """
    cdef uint32_t *idx = structure_element_index(structure) + 120
    cdef atom_t *atoms = structure.atoms()
    cdef uint32_t begin = element_bucket_begin(structure, 6)
    cdef uint32_t end = element_bucket_end(structure, 6)
    cdef uint32_t k
    cdef uint32_t count = 0
    for k in range(begin, end):
        if at_hybridization(&atoms[idx[k]]) == 1:
            count += 1
    return count


cdef inline uint32_t desc_heteroatoms_count(Structure structure) noexcept nogil:
    """How many atoms are neither carbon, hydrogen, nor R (element 0).

    A COUNT OF ATOMS, and deliberately not a sum of the per-atom `a.heteroatoms`, which counts an
    atom's heteroatom NEIGHBOURS -- summing that counts each heteroatom once per bond and answers a
    different question.  Hydrogen and R (element 0) are not heteroatoms, which is `derive_scalars`'
    rule for `a.heteroatoms` too, so the two descriptors agree.
    """
    return (structure.header.atom_count
            - (element_bucket_end(structure, 1) - element_bucket_begin(structure, 1))
            - (element_bucket_end(structure, 6) - element_bucket_begin(structure, 6))
            # Bucket 0 is the R atoms: a marker is not a heteroatom, on either side of the subtraction.
            - (element_bucket_end(structure, 0) - element_bucket_begin(structure, 0)))


# --- ring systems ---------------------------------------------------------------------------------
#
# ALL OF IT IS A FOLD OVER THE MINIMUM CYCLE BASIS the arena already perceived.  `structure_rings`
# gives the count, the per-ring member ranges and the members in cyclic walk order, which is what lets
# a ring's bonds be recovered with `csr_find` on consecutive pairs -- the same walk
# `aromatic_rings` does in Python, and the reason these counts cannot disagree with it.
#
# A COUNT OF BASIS RINGS IS BASIS-DEPENDENT and that is stated rather than hidden: the basis is a
# minimum cycle basis, so its ring SIZES are canonical, but which cycles were chosen among equal-size
# alternatives is not.  Every count here is invariant under that choice for the molecules chemistry
# cares about, and none of them is a claim about the exponential relevant-cycle set.


cdef void desc_ring_classes(Structure structure, uint32_t *out) noexcept nogil:
    """Five classifications of the ring basis into `out[0:5]`, which the caller owns.

    out[0] aromatic       -- every bond in the ring is stored order 4
    out[1] aliphatic      -- not aromatic, so out[0] + out[1] == rings_count always
    out[2] saturated      -- every bond in the ring is stored order 1
    out[3] heterocyclic   -- holds an atom that is neither carbon, hydrogen, nor R (element 0)
    out[4] aromatic heterocyclic -- both of the above

    THE ORDERS ARE READ, NOT PERCEIVED.  A kekulized benzene has no order-4 bond and answers 0
    aromatic rings, exactly as `aromatic_rings` does; `thiele()` is what changes the answer, and
    nothing here calls it.  Aromatic and saturated are not complements: tetralin's carbocycle is
    neither, because it shares one order-4 bond with the arene.
    """
    cdef uint32_t i
    for i in range(5):
        out[i] = 0
    if not structure_has(structure, SEG_RELEVANT_RINGS):
        return
    cdef uint32_t *r = structure_rings(structure)
    cdef uint32_t count = r[0]
    cdef uint32_t base = 2 + count
    cdef atom_t *atoms = structure.atoms()
    cdef halfedge_t *e
    cdef uint32_t k, cur, prev, z
    cdef bint aromatic, saturated, hetero
    for i in range(count):
        aromatic = True
        saturated = True
        hetero = False
        # the closing bond first, then the rest -- a carried `prev` and never `members[-1]`, since
        # this unit compiles with wraparound=False (RULES.md 7.6, and `aromatic_rings` says the same)
        prev = r[base + r[2 + i] - 1]
        for k in range(r[1 + i], r[2 + i]):
            cur = r[base + k]
            z = atoms[cur].element
            if element_is_heteroatom(z):
                hetero = True
            e = csr_find(structure, prev, cur)
            if e is NULL:      # unreachable with a consistent basis; a missing bond is not aromatic
                aromatic = False
                saturated = False
            else:
                if e.order != 4:
                    aromatic = False
                if e.order != 1:
                    saturated = False
            prev = cur
        if aromatic:
            out[0] += 1
        else:
            out[1] += 1
        if saturated:
            out[2] += 1
        if hetero:
            out[3] += 1
            if aromatic:
                out[4] += 1


cdef inline uint32_t _desc_uf_find(uint32_t *parent, uint32_t x) noexcept nogil:
    """Union-find root with halving path compression.  Used once, by `desc_ring_atoms`."""
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x


cdef int desc_ring_atoms(Structure structure, uint32_t *out) except -1:
    """Spiro atoms, bridgehead atoms and ring systems into `out[0:3]`, which the caller owns.

    SPIRO: an atom shared by two basis rings that share nothing else.  Sharing exactly one atom is
    the whole definition -- two rings cannot share one atom and any bond.

    BRIDGEHEAD: an atom shared by two basis rings that share AT LEAST TWO BONDS, and itself incident
    to at least three ring bonds.  Both clauses earn their place on norbornane, whose two five-rings
    share three atoms and two bonds: the middle atom of that shared path is the one-carbon bridge and
    carries only two ring bonds, so the second clause is what makes the answer 2 rather than 3.  The
    first clause is what makes naphthalene 0 -- its rings share one bond, which is fusion, not
    bridging, even though both fusion atoms do carry three ring bonds.

    RING SYSTEMS: connected components of the subgraph of ring bonds.  An isolated ring is one system;
    a spiro atom merges its two rings, because the bonds of both are incident to it.  Additive over
    components of the molecule, so a salt needs no special case.

    Ring pairs are compared by STAMPING rather than by building intersections: ring i writes `i + 1`
    into an atom slot and both half-edge slots of each of its bonds, and ring j then reads those slots.
    A stale stamp from an earlier ring can never match, so nothing is cleared between rings.  Costs
    O(rings^2 * ring size), which at chemical sizes is nothing and needs no index.
    """
    out[0] = 0
    out[1] = 0
    out[2] = 0
    cdef uint32_t n_atoms = structure.header.atom_count
    if n_atoms == 0:
        return 0
    cdef atom_t *atoms = structure.atoms()
    cdef uint32_t *ptr = csr_ptr(structure)
    cdef halfedge_t *edges = csr_edges(structure)
    cdef uint32_t n_edges = ptr[n_atoms]
    cdef uint32_t *r = NULL
    cdef uint32_t count = 0
    if structure_has(structure, SEG_RELEVANT_RINGS):
        r = structure_rings(structure)
        count = r[0]
    cdef uint32_t base = 2 + count

    # ONE BLOCK, FIVE REGIONS (RULES.md 5.2): the atom stamp, the union-find parent, the half-edge
    # stamp, and one byte per atom for each of the two verdicts.  The two uint8 regions go last so the
    # three uint32 regions stay naturally aligned.
    cdef void *block = PyMem_Malloc(<size_t> n_atoms * (2 * sizeof(uint32_t) + 2)
                                    + <size_t> n_edges * sizeof(uint32_t))
    if block is NULL:
        raise MemoryError('ring descriptor scratch allocation failed')
    cdef uint32_t *atom_stamp = <uint32_t *> block
    cdef uint32_t *parent = atom_stamp + n_atoms
    cdef uint32_t *edge_stamp = parent + n_atoms
    cdef uint8_t *spiro = <uint8_t *> (edge_stamp + n_edges)
    cdef uint8_t *bridge = spiro + n_atoms
    cdef halfedge_t *e
    cdef uint32_t i, j, k, m, cur, prev, ra, rb, shared_atoms, shared_bonds, ring_bonds, last_shared
    try:
        with nogil:
            for i in range(n_atoms):
                atom_stamp[i] = 0
                parent[i] = i
                spiro[i] = 0
                bridge[i] = 0
            for k in range(n_edges):
                edge_stamp[k] = 0

            # ring systems: union the endpoints of every ring bond, then count the roots that are
            # ring atoms.  A non-ring atom is never unioned, so it stays a singleton root and is
            # excluded by the `at_in_ring` test rather than by a second pass.
            for i in range(n_atoms):
                for k in range(ptr[i], ptr[i + 1]):
                    if edges[k].flags & HE_IN_RING:
                        ra = _desc_uf_find(parent, i)
                        rb = _desc_uf_find(parent, edges[k].to)
                        if ra != rb:
                            parent[ra] = rb
            for i in range(n_atoms):
                if at_in_ring(&atoms[i]) and _desc_uf_find(parent, i) == i:
                    out[2] += 1

            for i in range(count):
                prev = r[base + r[2 + i] - 1]
                for k in range(r[1 + i], r[2 + i]):
                    cur = r[base + k]
                    atom_stamp[cur] = i + 1
                    e = csr_find(structure, prev, cur)
                    if e is not NULL:
                        edge_stamp[<uint32_t> (e - edges)] = i + 1
                    e = csr_find(structure, cur, prev)
                    if e is not NULL:
                        edge_stamp[<uint32_t> (e - edges)] = i + 1
                    prev = cur
                for j in range(i + 1, count):
                    shared_atoms = 0
                    shared_bonds = 0
                    last_shared = 0
                    prev = r[base + r[2 + j] - 1]
                    for k in range(r[1 + j], r[2 + j]):
                        cur = r[base + k]
                        if atom_stamp[cur] == i + 1:
                            shared_atoms += 1
                            last_shared = cur
                        e = csr_find(structure, prev, cur)
                        if e is not NULL and edge_stamp[<uint32_t> (e - edges)] == i + 1:
                            shared_bonds += 1
                        prev = cur
                    if shared_atoms == 1:
                        spiro[last_shared] = 1
                    elif shared_bonds >= 2:
                        for k in range(r[1 + j], r[2 + j]):
                            cur = r[base + k]
                            if atom_stamp[cur] == i + 1:
                                ring_bonds = 0
                                for m in range(ptr[cur], ptr[cur + 1]):
                                    if edges[m].flags & HE_IN_RING:
                                        ring_bonds += 1
                                if ring_bonds >= 3:
                                    bridge[cur] = 1

            for i in range(n_atoms):
                out[0] += spiro[i]
                out[1] += bridge[i]
    finally:
        PyMem_Free(block)
    return 0


# --- distances ------------------------------------------------------------------------------------
#
# ONE MATRIX FOR FOUR DESCRIPTORS, and `balaban_j` makes it five.  `mol_distance_matrix` runs one BFS
# per atom and is the only shortest-path code in the core; every index below reads its output.  Adding
# a second traversal here would be two implementations of one question, which is how they drift.
#
# -1 IS NOT A DISTANCE.  The matrix reports -1 for a pair in different components (see its docstring
# for why that encoding and not a large sentinel), so a fold over it either skips -1 or is wrong.  The
# per-descriptor consequences are in the docstrings: W skips the pair, an eccentricity is the largest
# distance to an atom the vertex can REACH, and the radius of a molecule containing a lone counterion
# is therefore 0.


cdef object desc_eccentricities(Structure structure):
    """A `(n,)` int32 array: the largest distance from each atom to an atom it can reach.

    Indexed exactly like `distance_matrix`'s rows -- entry `i` belongs to `atom_numbers[i]`.  An atom
    that can reach nothing (a lone counterion, a one-atom molecule) gets 0, which is what makes the
    number within-component rather than infinite.
    """
    _numpy_load()
    cdef uint32_t n_atoms = structure.header.atom_count
    cdef object out = _NP_ZEROS(n_atoms, dtype='int32')
    if n_atoms == 0:
        return out
    cdef int32_t[:, ::1] dist = mol_distance_matrix(structure)
    cdef int32_t[::1] ecc = out
    cdef uint32_t i, j
    cdef int32_t best
    with nogil:
        for i in range(n_atoms):
            best = 0
            for j in range(n_atoms):
                if dist[i, j] > best:      # -1 loses to 0, so an unreachable pair is skipped
                    best = dist[i, j]
            ecc[i] = best
    return out


cdef int64_t desc_wiener(Structure structure) except? -1:
    """The Wiener index: the sum of topological distances over unordered pairs of atoms.

    Wiener, JACS 69 (1947) 17, where it is the "path number" w and is tabulated for the alkanes --
    n-butane 10, n-pentane 20, isopentane 18, neopentane 16.

    Defined on the hydrogen-suppressed graph, which is what chython stores when the hydrogens are
    implicit; a record carrying explicit hydrogen atoms is a different graph and gets a larger number
    (see the file header).  int64 because W grows as n^3 for a chain and a 65535-atom arena would
    overflow int32.

    A PAIR WITH NO PATH IS SKIPPED, which makes W additive over components: two butanes are 20.
    """
    cdef uint32_t n_atoms = structure.header.atom_count
    if n_atoms < 2:
        return 0
    cdef int32_t[:, ::1] dist = mol_distance_matrix(structure)
    cdef uint32_t i, j
    cdef int64_t acc = 0
    with nogil:
        for i in range(n_atoms):
            for j in range(i + 1, n_atoms):
                if dist[i, j] > 0:
                    acc += dist[i, j]
    return acc


cdef int desc_radius_diameter(Structure structure, int32_t *out) except -1:
    """The graph radius into `out[0]` and the diameter into `out[1]`; the caller owns the two slots.

    The minimum and the maximum eccentricity.  Both are 0 for an empty molecule, and the radius of any
    molecule with an isolated atom is 0 -- that atom's eccentricity is 0 and the minimum takes it.
    Computed off the matrix directly rather than off `desc_eccentricities`, so neither allocates for
    the other; the vector and these two agree because they fold the same rows the same way, which
    `test_radius_and_diameter_agree_with_the_eccentricity_vector` pins.
    """
    out[0] = 0
    out[1] = 0
    cdef uint32_t n_atoms = structure.header.atom_count
    if n_atoms == 0:
        return 0
    cdef int32_t[:, ::1] dist = mol_distance_matrix(structure)
    cdef uint32_t i, j
    cdef int32_t best
    cdef int32_t radius = -1
    cdef int32_t diameter = 0
    with nogil:
        for i in range(n_atoms):
            best = 0
            for j in range(n_atoms):
                if dist[i, j] > best:
                    best = dist[i, j]
            if radius < 0 or best < radius:
                radius = best
            if best > diameter:
                diameter = best
    out[0] = radius
    out[1] = diameter
    return 0


cdef int64_t desc_valence_electrons(MoleculeContainer mol) except? -1:
    """The molecule's valence electron count: sum of Zv - charge + implicit hydrogens.

    Zv is `el_valence_electrons`, the group number convention stated once in the header of
    `elements.tsv`.  The formal charge is subtracted because a cation has lost electrons, and each
    implicit hydrogen adds the one it brings.

    `explicit_h` IS NOT IN THE SUM, and that is what makes the arithmetic agree across the two
    spellings of one molecule: an explicit hydrogen is a vertex in this graph (see the file header) and
    contributes its own electron in its own iteration.  `read_smiles('C')` and
    `read_smiles('[H]C([H])([H])[H]')` both answer 8.

    REFUSES RATHER THAN GUESSING, twice.  An f-block atom has no stated Zv -- the 4f and 5f electrons
    are neither reliably core nor reliably valence -- and an atom with no implicit hydrogen count makes
    the sum undeterminable, exactly as `total_h_of` reports None per atom.  Both raise `ValueError`
    naming the atom, and the second names the exact `calc_implicit` call as the repair -- it is per-atom,
    not per-molecule, so a message that said "on the molecule" would earn a `TypeError`.  Additive over
    components, so
    a salt needs no special case.
    """
    cdef Structure structure = mol._structure
    cdef uint32_t n_atoms = structure.header.atom_count
    cdef atom_t *atoms = structure.atoms()
    cdef atom_t *a
    cdef uint32_t i
    cdef uint32_t zv
    cdef int64_t acc = 0
    for i in range(n_atoms):
        a = &atoms[i]
        zv = el_valence_electrons(a.element)
        if zv == VALENCE_ELECTRONS_UNKNOWN:
            raise ValueError(
                'atom %d is %s and chython states no valence electron count for the f block, so the '
                'molecule has none either; see the header of elements.tsv'
                % (mol._numbers[i], symbol_of(a)))
        if at_implicit_h_unknown(a):
            raise ValueError(
                'atom %d has no implicit hydrogen count, so the valence electron count is not '
                'derivable; run chython.chemistry.calc_implicit(molecule, %d) -- it takes one atom, '
                'and kekule() has to come first if that atom holds aromatic bonds'
                % (mol._numbers[i], mol._numbers[i]))
        acc += <int64_t> zv - <int64_t> a.charge + <int64_t> at_implicit_h(a)
    return acc


# --- degree indices -------------------------------------------------------------------------------
#
# DEGREE HERE IS THE CSR ROW LENGTH -- `ptr[i + 1] - ptr[i]`, which is what `degree_of()` reports and
# what the file header commits to.  An order-8 bond is a row entry and counts; the SMARTS `D`
# primitive skips it, and this is not that question.  Read off `ptr` rather than off `a.degree` for
# one reason: `ptr` is already the loop bound, so the length is free, and there is then no way for the
# two to disagree in this file.


cdef uint64_t desc_zagreb(Structure structure, bint second) noexcept nogil:
    """The first Zagreb index (`second` false) or the second (`second` true).

    Gutman and Trinajstic, Chem. Phys. Lett. 17 (1972) 535, where they arise as the two terms of a
    total pi-energy expansion:

        M1 = sum over atoms of deg(v)^2
        M2 = sum over bonds of deg(u) * deg(v)

    THE SELECTOR IS A `bint`, NOT THE PAPER'S ORDER NUMBER, and that is the point: a `nogil` function
    cannot raise, so an `order` parameter would have to treat some value as "everything that is not 1"
    and silently answer M2 for a caller that asked for order 3.  A two-valued parameter has no invalid
    value to mishandle.  `zagreb_index(order=...)` owns the 1-or-2 domain and refuses the rest, once,
    where a `ValueError` is reachable -- RULES.md 6.1, a field's domain declared exactly once.

    uint64 because M1 grows as the square of the maximum degree times n and there is no reason to make
    a caller think about a 32-bit edge.  Additive over components, so a salt needs no special case.
    """
    cdef uint32_t n_atoms = structure.header.atom_count
    if n_atoms == 0:
        return 0
    cdef uint32_t *ptr = csr_ptr(structure)
    cdef halfedge_t *edges = csr_edges(structure)
    cdef uint32_t i, k, to, di
    cdef uint64_t acc = 0
    if not second:
        for i in range(n_atoms):
            di = ptr[i + 1] - ptr[i]
            acc += (<uint64_t> di) * di
    else:
        for i in range(n_atoms):
            di = ptr[i + 1] - ptr[i]
            for k in range(ptr[i], ptr[i + 1]):
                to = edges[k].to
                if to > i:      # each bond once: the half-edge whose head has the larger index
                    acc += (<uint64_t> di) * (ptr[to + 1] - ptr[to])
    return acc


cdef double desc_randic(Structure structure) noexcept nogil:
    """The Randic branching index: sum over bonds of 1 / sqrt(deg(u) * deg(v)).

    Randic, JACS 97 (1975) 6609.  It is also the first-order connectivity index 1-chi of Kier and
    Hall, which is why `chi(1)` in this file answers the same number -- one definition reached by two
    names, and neither is a copy of the other's code.

    A degree-0 atom cannot reach the division: it is the endpoint of no bond.  0.0 for a molecule with
    no bonds, and additive over components.

    THE DEGREES ARE BOUND TO DOUBLES BEFORE THEY MEET, so the product is never a uint32 multiplication
    that could wrap.  Written as two locals rather than as one cast expression because a cast binds
    tighter than `*` and a reader who does not know that reads `<double> a * b` as covering the product;
    `chi()` copies this line, and the copy has to be unambiguous rather than merely correct.
    """
    cdef uint32_t n_atoms = structure.header.atom_count
    if n_atoms == 0:
        return 0.0
    cdef uint32_t *ptr = csr_ptr(structure)
    cdef halfedge_t *edges = csr_edges(structure)
    cdef uint32_t i, k, to
    cdef double di, dj
    cdef double acc = 0.0
    for i in range(n_atoms):
        di = ptr[i + 1] - ptr[i]
        for k in range(ptr[i], ptr[i + 1]):
            to = edges[k].to
            if to > i:
                dj = ptr[to + 1] - ptr[to]
                acc += 1.0 / sqrt(di * dj)
    return acc


# --- Balaban J ------------------------------------------------------------------------------------


cdef double desc_balaban_j(MoleculeContainer mol) except? -1.0:
    """Balaban's average distance sum connectivity index J.

        J = q / (mu + 1) * sum over bonds of 1 / sqrt(s(u) * s(v))

    q is the bond count, mu = q - n + 1 the cyclomatic number, and s(v) the sum of v's distances to
    every other atom -- a row sum of the distance matrix.  Balaban, Chem. Phys. Lett. 89 (1982) 399,
    "Highly discriminating distance-based topological index"; the short alkanes come out 1.000, 1.633,
    1.975, 2.191 with isobutane at 2.324, which is the discrimination the paper demonstrates.

    REFUSES A DISCONNECTED MOLECULE, and it is the only descriptor in this file that refuses anything.
    A row sum containing a -1 is not a distance sum; skipping the -1 would make s within-component
    while q and mu stayed global, which is a formula nobody defined.  So this raises ValueError naming
    `split()` -- a refusal at the answer boundary, where refusals belong -- and a caller who wants a J
    per part splits and asks each.  A ONE-ATOM MOLECULE IS NOT THAT CASE: it is connected, has no
    bonds, and its empty sum is honestly 0.0.  An empty molecule likewise.

    THERE IS NO DIVISION BY ZERO.  An acyclic molecule has mu = 0, so the denominator is mu + 1 = 1.
    A molecule with no bonds has q = 0, which makes the whole prefactor 0 before the (empty) sum is
    reached.

    TAKES THE CONTAINER so the message can quote `connected_components_count` -- the same number the
    caller would check -- rather than a count computed a second way here.
    """
    cdef Structure structure = mol._structure
    cdef uint32_t n_atoms = structure.header.atom_count
    if n_atoms < 2:                      # empty, or one atom with no bond to sum over
        return 0.0
    cdef uint32_t *ptr = csr_ptr(structure)
    cdef halfedge_t *edges = csr_edges(structure)
    cdef uint32_t n_bonds = ptr[n_atoms] // 2
    cdef int32_t[:, ::1] dist = mol_distance_matrix(structure)
    cdef uint32_t i, j, k, to, mu
    cdef int64_t row
    cdef bint disconnected = False
    cdef double acc = 0.0
    cdef double si, sj

    # ONE BLOCK, ONE REGION (RULES.md 5.2): the vertex distance sums, int64 because a row sum grows as
    # n^2 and the arena allows 65535 atoms
    cdef int64_t *s = <int64_t *> PyMem_Malloc(<size_t> n_atoms * sizeof(int64_t))
    if s is NULL:
        raise MemoryError('Balaban J scratch allocation failed')
    try:
        with nogil:
            for i in range(n_atoms):
                row = 0
                for j in range(n_atoms):
                    if dist[i, j] < 0:
                        disconnected = True
                        break
                    row += dist[i, j]
                if disconnected:
                    break
                s[i] = row
        if disconnected:
            raise ValueError(
                "Balaban's J is defined on a connected graph and this molecule has %d components; "
                'split() it and index each part' % mol.connected_components_count)
        # disconnected is False, so no -1 survives to the arithmetic below
        with nogil:
            for i in range(n_atoms):
                for k in range(ptr[i], ptr[i + 1]):
                    to = edges[k].to
                    if to > i:           # each bond once
                        si = <double> s[i]
                        sj = <double> s[to]
                        acc += 1.0 / sqrt(si * sj)
    finally:
        PyMem_Free(s)
    # connected, so q >= n - 1 >= 1 for n >= 2 and mu cannot underflow; mu + 1 >= 1 always
    mu = n_bonds - n_atoms + 1
    return <double> n_bonds / <double> (mu + 1) * acc


# --- Bertz CT -------------------------------------------------------------------------------------


cdef double desc_bertz_ct(Structure structure) except? -1.0:
    """Bertz's molecular complexity index CT.

        CT = 2*N*log2(N) - sum over connection classes of n*log2(n)
             + n_atoms*log2(n_atoms) - sum over elements of m*log2(m)

    Bertz, JACS 103 (1981) 3599, "The first general index of molecular complexity".  A CONNECTION is a
    pair of bonds sharing an atom -- a three-atom path -- so N is the sum of C(deg, 2) over atoms.

    THE PARTITION IS CHYTHON'S STATED READING, because the paper defines complexity in terms of
    equivalent connections without fixing how equivalence is decided, and every toolkit chose
    differently.  Here two connections are equivalent when their central atoms share a symmetry class
    and their outer atoms' classes agree as an UNORDERED PAIR.  The classes come from
    `compute_atoms_order`, which is `atoms_order` -- an equitable partition refined to a fixed point,
    which is a refinement of the orbit partition and equal to it for everything short of a
    strongly-regular pathology.  Keying on the centre alone would be wrong in an obvious way: phenol's
    ipso carbon carries three connections in two classes, not three in one.

    BOND ORDERS REACH THIS INDEX ONLY THROUGH THE RANKS, so benzene and cyclohexane get the same CT.
    That is a property of information-content indices -- two equally symmetric graphs of a size are
    equally complex -- and not something to patch by adding orders to the key: CT reads class SIZES,
    and a key that splits no class changes nothing.

    NO REFUSAL AND NO ADDITIVITY.  There is no distance in CT, so a disconnected molecule answers; two
    copies of one component enlarge N and its classes at once and the answer is not twice one copy's.

    ENUMERATES ONE CENTRE PER RANK, not all N connections.  Two atoms of the same rank have the same
    multiset of neighbour ranks -- that is what "equitable" means -- so one representative's pair
    histogram, multiplied by the rank's population, is the class size exactly.  The alternative, an
    array of N keys to sort, is bigger and slower for the same answer.  That shortcut is CHECKED, not
    trusted: the connections it accumulates are compared against the sum of C(deg, 2), which is the same
    number counted without reference to any partition.
    """
    cdef uint32_t n_atoms = structure.header.atom_count
    if n_atoms == 0:
        return 0.0
    cdef uint32_t *ptr = csr_ptr(structure)
    cdef halfedge_t *edges = csr_edges(structure)
    cdef uint32_t i, k, v, r, deg, run, run2, at, maxdeg = 0
    cdef uint32_t e
    cdef uint32_t population
    cdef uint32_t rank_i
    cdef uint32_t base
    cdef Py_ssize_t classes
    cdef uint64_t total = 0
    # SIZE AND PAIRS ARE 64-BIT BECAUSE A CLASS SIZE IS NOT BOUNDED BY THE ATOM COUNT.  A uint32 `size`
    # computes the product in 32 bits before widening into `total`, which is a silent wrong answer
    # rather than an overflow anyone sees: a star with 92,683 neighbours has
    # C(deg, 2) = 4,295,022,903 connections, just past 2**32, and CT came back 876,549.91 where the
    # formula gives 137,440,813,138.41.  The container accepts that graph without complaint, so the
    # descriptor has to survive it.  `pairs` is a separate local rather than a cast inside the expression
    # because a cast binds tighter than `*` and a reader should not have to know that to see the width.
    cdef uint64_t size
    cdef uint64_t pairs
    # the same C(deg, 2) sum, computed independently of the rank partition -- see the guard below
    cdef uint64_t expect = 0
    cdef double entropy = 0.0
    cdef double out = 0.0
    cdef size_t n_uint32
    for i in range(n_atoms):
        deg = ptr[i + 1] - ptr[i]
        if deg > maxdeg:
            maxdeg = deg
        if deg > 1:
            expect += <uint64_t> deg * (deg - 1) // 2

    # ONE BLOCK, FOUR REGIONS (RULES.md 5.2): the neighbour-rank buffer the sort works on (uint64
    # because `_sort_words` sorts uint64), then the per-atom rank, the population of each rank and one
    # representative atom per rank.  The last two are sized n_atoms + 1 rather than classes + 1
    # because the class count is only known after the refinement runs, and n_atoms is its bound.
    # uint64 region first so the uint32 regions stay aligned behind it.
    n_uint32 = <size_t> (3 * n_atoms + 2)         # hoisted so the cast covers the full product width
    cdef void *block = PyMem_Malloc(<size_t> maxdeg * sizeof(uint64_t)
                                    + n_uint32 * sizeof(uint32_t))
    if block is NULL:
        raise MemoryError('Bertz CT scratch allocation failed')
    cdef uint64_t *nbr = <uint64_t *> block
    cdef uint32_t *rank = <uint32_t *> (nbr + maxdeg)
    cdef uint32_t *population_of = rank + n_atoms
    cdef uint32_t *rep = population_of + n_atoms + 1
    try:
        with nogil:
            classes = compute_atoms_order(structure, rank, NULL)
        if classes < 0:
            raise MemoryError('atom order refinement failed to allocate')
        with nogil:
            for r in range(<uint32_t> classes + 1):
                population_of[r] = 0
                rep[r] = n_atoms          # sentinel: no representative seen yet
            for i in range(n_atoms):
                rank_i = rank[i]          # bound once -- RULES.md 2.1, three reads of one subscript
                population_of[rank_i] += 1
                if rep[rank_i] == n_atoms:
                    rep[rank_i] = i

            for r in range(1, <uint32_t> classes + 1):
                v = rep[r]
                base = ptr[v]
                deg = ptr[v + 1] - base
                if deg < 2:               # no pair of bonds to make a connection from
                    continue
                population = population_of[r]
                for k in range(deg):
                    nbr[k] = rank[edges[base + k].to]
                _sort_words(nbr, deg)
                # runs of equal neighbour rank: C(run, 2) connections pair a run with itself, and
                # run * run2 pair two different runs.  Each key's count at this centre, times the
                # rank's population, is the class size.
                i = 0
                while i < deg:
                    run = 1
                    while i + run < deg and nbr[i + run] == nbr[i]:
                        run += 1
                    if run > 1:
                        pairs = <uint64_t> run * (run - 1) // 2
                        size = pairs * population
                        total += size
                        entropy += <double> size * log2(<double> size)
                    at = i + run
                    while at < deg:
                        run2 = 1
                        while at + run2 < deg and nbr[at + run2] == nbr[at]:
                            run2 += 1
                        pairs = <uint64_t> run * run2
                        size = pairs * population
                        total += size
                        entropy += <double> size * log2(<double> size)
                        at += run2
                    i += run
    finally:
        PyMem_Free(block)

    # THE ONE REPRESENTATIVE PER RANK IS LOAD-BEARING, SO IT IS CHECKED RATHER THAN TRUSTED.  Summing one
    # representative's pair histogram times its rank's population reproduces every connection only while
    # the partition is equitable; `compute_atoms_order` also stops refining when a round produces no new
    # class, which an XXH64 collision can cause.  Total connections are the same number counted two ways,
    # so the disagreement is O(n) to detect -- and a wrong complexity nobody can spot is worse than a
    # refusal.  This also pins the assumption against any future refinement of `_morgan.pxi`.
    if total != expect:
        raise RuntimeError(
            'Bertz CT counted %d connections over the symmetry classes where the degrees give %d; the '
            'atom order partition is not equitable, which this index depends on' % (total, expect))

    if total > 0:
        out = 2.0 * <double> total * log2(<double> total) - entropy
    # the element diversity term.  Sizes are >= 1 so log2 never sees 0, and a molecule of one element
    # contributes exactly 0.
    out += <double> n_atoms * log2(<double> n_atoms)
    for e in range(1, 119):
        size = element_bucket_end(structure, e) - element_bucket_begin(structure, e)
        if size > 0:
            out -= <double> size * log2(<double> size)
    return out


# --- connectivity indices -------------------------------------------------------------------------
#
# THE ONE PATH ENUMERATOR IN THE FILE.  `desc_chi` walks simple paths of a stated edge count and
# weights each by the reciprocal square root of the product of its atoms' deltas.  Kappa (next section)
# needs the same paths UNWEIGHTED, and gets them by calling this with an all-ones delta -- every path
# then contributes exactly 1.0, and a sum of 1.0s is an exact integer well past any arena's size.  One
# enumerator, so a path counted for kappa is a path counted for chi.
#
# A PATH IS COUNTED ONCE, by requiring the far endpoint's index to exceed the start's.  For order 0
# there is no path and no endpoint: the sum is over atoms.
#
# A DELTA OF 0 CONTRIBUTES NOTHING, and this is chython's stated reading rather than the paper's: the
# term would be 1/sqrt(0).  It happens for the plain delta at an atom with no heavy neighbour, and for
# the valence delta at a fully hydrogenated one (methane's carbon, 4 - 4).  Any path through such an
# atom is skipped for the same reason, so the rule is one rule.


cdef void _desc_delta_plain(Structure structure, double *delta) noexcept nogil:
    """The simple delta: the atom's degree, which is its CSR row length."""
    cdef uint32_t n_atoms = structure.header.atom_count
    if n_atoms == 0:
        return
    cdef uint32_t *ptr = csr_ptr(structure)
    cdef uint32_t i
    for i in range(n_atoms):
        delta[i] = <double> (ptr[i + 1] - ptr[i])


cdef int _desc_delta_valence(MoleculeContainer mol, double *delta) except -1:
    """The valence delta: Zv - h, valence electrons less hydrogens (Kier and Hall).

    THE FORMAL CHARGE IS NOT SUBTRACTED.  Kier and Hall define delta-v over the element's valence
    electrons and the atom's hydrogens, and that is what this computes -- which is why it is not
    `desc_valence_electrons`' per-atom term, where the charge IS subtracted because that function counts
    electrons rather than free connections.  Two formulas, both stated, neither derived from the other.

    h IS THE IMPLICIT COUNT ONLY, and not the total count.  `desc_valence_electrons` states the reason:
    an explicit hydrogen is a vertex in this graph and contributes its own delta-v of 1 in the same sum,
    so adding it to its heavy neighbour's count subtracts it twice.  Kier and Hall define delta-v on the
    hydrogen-suppressed graph, where the explicit count is 0 and the question does not arise.

    Refuses exactly what it cannot state: an f-block element has no Zv, and an unknown implicit
    hydrogen count makes h unknown.  A fully hydrogenated atom is not a refusal -- Zv - h is 0 and the
    caller of this gets the "contributes nothing" rule.
    """
    cdef Structure structure = mol._structure
    cdef uint32_t n_atoms = structure.header.atom_count
    cdef atom_t *atoms = structure.atoms()
    cdef atom_t *a
    cdef uint32_t i, zv, h
    for i in range(n_atoms):
        a = &atoms[i]
        zv = el_valence_electrons(a.element)
        if zv == VALENCE_ELECTRONS_UNKNOWN:
            raise ValueError(
                'atom %d is %s and chython states no valence electron count for the f block, so its '
                'valence delta is not derivable; see the header of elements.tsv'
                % (mol._numbers[i], symbol_of(a)))
        if at_implicit_h_unknown(a):
            raise ValueError(
                'atom %d has no implicit hydrogen count, so its valence delta is not derivable; run '
                'chython.chemistry.calc_implicit on the molecule first' % mol._numbers[i])
        h = at_implicit_h(a)
        if h >= zv:
            delta[i] = 0.0
        else:
            delta[i] = <double> (zv - h)
    return 0


cdef double _desc_chi_walk(uint32_t *ptr, halfedge_t *edges, double *delta, uint8_t *visited,
                           uint32_t start, uint32_t v, uint32_t remaining, double prod) noexcept nogil:
    """Depth-first extension of one simple path; returns the summed weights of every completion.

    `remaining` is how many edges are still to be walked, `prod` the running product of the deltas of
    the atoms already on the path.  Recursion depth is `remaining`, which the caller caps at 4.
    """
    cdef double acc = 0.0
    cdef uint32_t k, to
    if remaining == 0:
        if v > start:      # each path once, from the end with the smaller index
            return 1.0 / sqrt(prod)
        return 0.0
    for k in range(ptr[v], ptr[v + 1]):
        to = edges[k].to
        if visited[to] or delta[to] <= 0.0:
            continue
        visited[to] = 1
        acc += _desc_chi_walk(ptr, edges, delta, visited, start, to, remaining - 1,
                             prod * delta[to])
        visited[to] = 0
    return acc


cdef double desc_chi(Structure structure, uint32_t order, double *delta,
                     uint8_t *visited) noexcept nogil:
    """The order-`order` connectivity index over the caller's delta vector.

    order 0 is the sum over atoms of 1/sqrt(delta); order m > 0 is the sum over simple paths of m edges
    of 1/sqrt(product of the deltas along the path).  Kier and Hall, Rev. Comput. Chem. 2 (1991)
    367-422; order 1 is Randic's branching index, which is why `randic_index` and `chi(1)` agree.

    `visited` is scratch of at least `n_atoms` bytes and its contents on entry are irrelevant.  THE
    CALLER VALIDATES `order`; nothing here caps the recursion.
    """
    cdef uint32_t n_atoms = structure.header.atom_count
    cdef uint32_t i
    cdef double acc = 0.0
    if n_atoms == 0:
        return 0.0
    if order == 0:
        for i in range(n_atoms):
            if delta[i] > 0.0:
                acc += 1.0 / sqrt(delta[i])
        return acc
    cdef uint32_t *ptr = csr_ptr(structure)
    cdef halfedge_t *edges = csr_edges(structure)
    for i in range(n_atoms):
        visited[i] = 0
    for i in range(n_atoms):
        if delta[i] <= 0.0:
            continue
        visited[i] = 1
        acc += _desc_chi_walk(ptr, edges, delta, visited, i, i, order, delta[i])
        visited[i] = 0
    return acc


cdef double desc_chi_of(MoleculeContainer mol, uint32_t order, bint valence) except? -1.0:
    """`chi()`'s body: own the scratch, build the requested delta, walk.

    Here rather than in the wrapper so the container keeps one line per descriptor, and so the delta
    choice sits beside the two builders.  ONE BLOCK, TWO REGIONS (RULES.md 5.2): the delta vector and
    the visited marks, doubles first for alignment.
    """
    cdef uint32_t n_atoms = mol._structure.header.atom_count
    if n_atoms == 0:
        return 0.0
    cdef void *block = PyMem_Malloc(<size_t> n_atoms * (sizeof(double) + 1))
    if block is NULL:
        raise MemoryError('connectivity index scratch allocation failed')
    cdef double *delta = <double *> block
    cdef uint8_t *visited = <uint8_t *> (delta + n_atoms)
    cdef double out
    try:
        if valence:
            _desc_delta_valence(mol, delta)
        else:
            _desc_delta_plain(mol._structure, delta)
        with nogil:
            out = desc_chi(mol._structure, order, delta, visited)
    finally:
        PyMem_Free(block)
    return out


# --- electrotopological state ---------------------------------------------------------------------
#
# Kier and Hall, Pharm. Res. 1990, 7, 801.  The intrinsic state is `I = ((2/N)**2 * dv + 1) / d` with N
# the period, dv the valence delta and d the plain delta; the electrotopological state adds every other
# atom's perturbation, `S_i = I_i + sum_j (I_i - I_j) / (dist_ij + 1)**2`.
#
# THE SAME TWO DELTAS AS CHI, AND THE SAME TWO REFUSALS.  `_desc_delta_valence` owns them, so an atom
# with no stated Zv or no stated hydrogen count has no intrinsic state either and the message it already
# writes is the one a caller sees.
#
# `d == 0` IS `nan`, NEVER 0.0.  An atom with no heavy neighbour puts a zero in the intrinsic state's
# denominator, and 0.0 is a perfectly ordinary EState value -- so a zero here would be a measurement
# nobody could tell from a gap.  `dv == 0` is NOT this case: methane's carbon has dv = 0 and I = 1/d,
# and `desc_chi`'s "a delta of 0 contributes nothing" rule does not reach here because dv is in a
# numerator rather than under a root.


cdef int _desc_estate_intrinsic(MoleculeContainer mol, double *out) except -1:
    """Kier and Hall's intrinsic state per atom, `nan` where the plain delta is 0.

    ONE BLOCK, TWO REGIONS (RULES.md 5.2): the valence delta and the plain delta, both doubles.
    """
    cdef Structure structure = mol._structure
    cdef uint32_t n_atoms = structure.header.atom_count
    if n_atoms == 0:
        return 0
    cdef void *block = PyMem_Malloc(<size_t> n_atoms * 2 * sizeof(double))
    if block is NULL:
        raise MemoryError('electrotopological state scratch allocation failed')
    cdef double *dv = <double *> block
    cdef double *dp = dv + n_atoms
    cdef atom_t *atoms = structure.atoms()
    cdef uint32_t i, period
    cdef double f
    try:
        _desc_delta_valence(mol, dv)
        _desc_delta_plain(structure, dp)
        with nogil:
            for i in range(n_atoms):
                if dp[i] == 0.0:
                    out[i] = NAN
                    continue
                period = el_period(atoms[i].element)
                f = 2.0 / <double> period
                out[i] = (f * f * dv[i] + 1.0) / dp[i]
    finally:
        PyMem_Free(block)
    return 0


cdef int _desc_estate(MoleculeContainer mol, double *out) except -1:
    """The electrotopological state: the intrinsic state plus every other atom's perturbation.

    A `nan` atom is skipped on BOTH sides -- it contributes to nobody's sum and gets none of its own --
    so a salt's organic component answers exactly what it answers on its own.  A pair in different
    components is skipped for the same reason with a different cause: `distance_matrix` says -1 and
    nothing here sums a -1.

    O(n^2) over the distance matrix, which is the definition and not an implementation choice.
    """
    cdef uint32_t n_atoms = mol._structure.header.atom_count
    if n_atoms == 0:
        return 0
    cdef void *block = PyMem_Malloc(<size_t> n_atoms * sizeof(double))
    if block is NULL:
        raise MemoryError('electrotopological state scratch allocation failed')
    cdef double *intrinsic = <double *> block
    cdef int32_t[:, ::1] dist
    cdef uint32_t i, j
    cdef double acc, step
    try:
        _desc_estate_intrinsic(mol, intrinsic)
        dist = mol_distance_matrix(mol._structure)
        with nogil:
            for i in range(n_atoms):
                if isnan(intrinsic[i]):
                    out[i] = NAN
                    continue
                acc = intrinsic[i]
                for j in range(n_atoms):
                    if j == i or isnan(intrinsic[j]) or dist[i, j] < 0:
                        continue
                    step = <double> dist[i, j] + 1.0
                    acc += (intrinsic[i] - intrinsic[j]) / (step * step)
                out[i] = acc
    finally:
        PyMem_Free(block)
    return 0


cdef object desc_estate_of(MoleculeContainer mol, bint intrinsic):
    """`estate_indices()`'s and `estate_intrinsic_states()`'s body: own the array, fill it, hand it over.

    Here rather than in the wrappers so the container keeps one line per descriptor and so the two
    spellings cannot come to disagree about the dtype or the order.
    """
    _numpy_load()
    cdef uint32_t n_atoms = mol._structure.header.atom_count
    cdef object out = _NP_EMPTY(n_atoms, dtype='float64')
    if n_atoms == 0:
        return out
    cdef double[::1] values = out
    if intrinsic:
        _desc_estate_intrinsic(mol, &values[0])
    else:
        _desc_estate(mol, &values[0])
    return out


# --- shape indices --------------------------------------------------------------------------------


cdef inline double desc_alpha_of(atom_t *a) noexcept nogil:
    """One atom's Hall-Kier alpha: its covalent radius relative to an sp3 carbon's, less one.

    Hall and Kier, Rev. Comput. Chem. 2 (1991) 367-422.  The values below are their published table to
    two decimals, and each one is r_cov / 0.77 - 1 for the radii of the period -- Csp2 0.67, Csp 0.60,
    Nsp3 0.74, Nsp2 0.62, Nsp 0.55, Osp3 0.74, Osp2 0.62, F 0.72, Psp3 1.10, Psp2 1.00, Ssp3 1.04,
    Ssp2 0.94, Cl 0.99, Br 1.14, I 1.33 -- so the table can be re-derived rather than trusted.

    RE-DERIVING IT REPRODUCES 12 OF THE 15 EXACTLY AND THE OTHER THREE ONE HUNDREDTH LOW, because BOTH
    columns are the paper's roundings and the radius is the one being fed to the arithmetic.  Nsp2 and
    Osp2 need 0.616 A to give -0.20 and F needs 0.716 to give -0.07; the paper prints 0.62 and 0.72, and
    0.62 / 0.77 - 1 is -0.19 while 0.72 / 0.77 - 1 is -0.06.  The three are named here so that a
    maintainer who takes the invitation above literally finds the discrepancy explained rather than
    "correcting" the code to values the paper does not publish -- which would shift alpha by 0.01 per
    sp2 nitrogen or oxygen, small enough to slip past a tolerance and large enough to be wrong.

    THE HYBRIDIZATION MAPPING IS STATED HERE because V3's z is 1..6 and the paper's rows are three.
    z1 is sp3; z4 (aromatic) reads as sp2; z3 and z6 read as sp.  z5 -- two cumulated double bonds --
    is sp for CARBON, where it is an allene centre and genuinely sp, and sp2 for NITROGEN, where it is a
    planar nitro group.

    SULFUR AND PHOSPHORUS ARE THE EXCEPTION, AND THE DISCRIMINATOR IS COORDINATION AND NOT z.  Alpha is a
    covalent-radius correction, and a radius tracks how many sigma bonds an atom holds, not how many
    formal double bonds were written on it.  A sulfone or sulfonamide S is z5 and a phosphate P is z2, but
    both are four-coordinate and tetrahedral and keep the sp3 radius -- 1.04 A and 1.10 A -- so reading
    either as sp2 understates alpha by 0.13 on every sulfonamide in a drug-like set.  A sulfoxide S is z2
    and three-coordinate and pyramidal, so it is sp3 too.

    The shortened sp2 radius belongs to a LOW-COORDINATE atom that is genuinely pi-bonded, so the rule is
    `hyb != 1 and degree <= 2`: thiophene S (z4, two bonds) and a thioketone S (z2, one bond) take 0.22,
    while a thioether S (z1) keeps 0.35 because it has no pi bond to shorten it.  Measured against
    chython's own perception rather than assumed -- CS(=O)(=O)C is z5/degree 4, CS(=O)C is z2/degree 3,
    OP(=O)(O)O is z2/degree 4, c1ccsc1 is z4/degree 2, CC(=S)C is z2/degree 1, CSC is z1/degree 2.

    The degree here is the structural one, so a dative contact counts towards coordination -- which is
    what a radius argument wants, unlike the SMARTS `D` primitive that deliberately excludes it.

    AN ELEMENT THE TABLE OMITS CONTRIBUTES 0.0, which is chython's reading and not a measurement: 0.0
    means "as far from an sp3 carbon as an sp3 carbon is", so a metal or an explicit hydrogen adds
    nothing.  It is what makes kappa answerable for an organometallic instead of a refusal, and the
    alternative -- extrapolating a radius the paper never published -- is invented chemistry.
    """
    cdef uint32_t z = a.element
    cdef uint32_t hyb = at_hybridization(a)
    # the sp2 row of the second period is a low-coordinate pi-bonded atom; see the docstring
    cdef bint shortened = hyb != 1 and a.degree <= 2
    if z == 6:
        if hyb == 1:
            return 0.0
        elif hyb == 2 or hyb == 4:
            return -0.13
        return -0.22                      # sp, and z5's allene centre is sp too
    elif z == 7:
        if hyb == 1:
            return -0.04
        elif hyb == 3 or hyb == 6:
            return -0.29
        return -0.20                      # sp2, aromatic, and nitro
    elif z == 8:
        if hyb == 1:
            return -0.04
        return -0.20                      # the paper has no sp oxygen
    elif z == 9:
        return -0.07
    elif z == 15:
        if shortened:
            return 0.30
        return 0.43                       # sp3, and a four-coordinate phosphate P is sp3 whatever z says
    elif z == 16:
        if shortened:
            return 0.22
        return 0.35                       # sp3, and that includes the sulfoxide, sulfone and sulfonamide
    elif z == 17:
        return 0.29
    elif z == 35:
        return 0.48
    elif z == 53:
        return 0.73
    return 0.0


cdef double desc_hall_kier_alpha(Structure structure) noexcept nogil:
    """The molecule's Hall-Kier alpha: the sum of its atoms' contributions.

    0.0 for a saturated hydrocarbon, negative for anything with sp2 or sp atoms or small heteroatoms,
    positive for the heavy halogens and the second-row heteroatoms.  Additive over components by
    construction.
    """
    cdef uint32_t n_atoms = structure.header.atom_count
    cdef atom_t *atoms = structure.atoms()
    cdef uint32_t i
    cdef double acc = 0.0
    for i in range(n_atoms):
        acc += desc_alpha_of(&atoms[i])
    return acc


cdef double desc_kappa(Structure structure, uint32_t order, bint alpha) except? -1.0:
    """Kier's kappa shape index of order 1, 2 or 3; `alpha` applies the Hall-Kier correction.

        kappa1 = n(n-1)^2 / P1^2
        kappa2 = (n-1)(n-2)^2 / P2^2
        kappa3 = (n-1)(n-3)^2 / P3^2   for n odd
                 (n-3)(n-2)^2 / P3^2   for n even

    Each compares the molecule's path count against the counts of the extremal graphs with the same
    atom count -- the star and the chain -- so it reads as "how linear is this".  Kier and Hall,
    Rev. Comput. Chem. 2 (1991) 367-422.  THE ORDER-3 PARITY SPLIT IS THE PUBLISHED ONE and is taken as
    published; it comes from the extremal graph for three-bond paths differing between odd and even n.

    With `alpha`, n becomes n + alpha and P becomes P + alpha, which is the published correction: it
    shrinks the effective atom count towards what a molecule of sp3 carbons of the same shape would
    have.  The parity still switches on the integer atom count.

    P_m COMES FROM `desc_chi` WITH AN ALL-ONES DELTA -- every path contributes 1/sqrt(1) = 1.0, so the
    sum is the count, exactly, and there is one path enumerator in this file rather than two.

    NO PATH OF THAT LENGTH MEANS 0.0, not a refusal and not a nan: isobutane has no three-bond path, so
    it has no three-bond shape.  A nan would poison every descriptor row this number lands in.
    """
    cdef uint32_t n_atoms = structure.header.atom_count
    cdef void *block
    cdef double *delta
    cdef uint8_t *visited
    cdef uint32_t i
    cdef double paths
    cdef double a
    cdef double n
    cdef double p
    if n_atoms == 0:
        return 0.0

    # one block, two regions (RULES.md 5.2): the unit delta and the visited marks
    block = PyMem_Malloc(<size_t> n_atoms * (sizeof(double) + 1))
    if block is NULL:
        raise MemoryError('kappa scratch allocation failed')
    delta = <double *> block
    visited = <uint8_t *> (delta + n_atoms)
    try:
        with nogil:
            for i in range(n_atoms):
                delta[i] = 1.0
            paths = desc_chi(structure, order, delta, visited)
    finally:
        PyMem_Free(block)
    if paths <= 0.0:
        return 0.0

    a = 0.0
    if alpha:
        a = desc_hall_kier_alpha(structure)
    n = <double> n_atoms + a
    p = paths + a
    if p == 0.0:
        return 0.0
    if order == 1:
        return n * (n - 1.0) * (n - 1.0) / (p * p)
    elif order == 2:
        return (n - 1.0) * (n - 2.0) * (n - 2.0) / (p * p)
    elif n_atoms % 2:
        return (n - 1.0) * (n - 3.0) * (n - 3.0) / (p * p)
    return (n - 3.0) * (n - 2.0) * (n - 2.0) / (p * p)
