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
# THE MOLECULE'S TOPOLOGY, as things a caller can hold: its connected components as separate
# molecules, the radius-bounded environment of a set of atoms, and the two matrix forms of the
# graph.  Four operations, one subject -- a walk over the CSR and nothing about chemistry.
#
# WHY A SEPARATE FRAGMENT.  `_molecule_container.pxi` is the class body and is edited by every epic
# at once; these four are free `cdef` functions over `Structure` (or, where they must build one, over
# `MoleculeContainer`) with three-line wrappers left behind in the class.  The split is not
# cosmetic: `mol_distance_matrix` is a nested loop with a queue, and a nested loop belongs where it
# can be read without scrolling past a journal.
#
# INCLUDED AFTER `_molecule_container.pxi`, because `mol_split` builds `MoleculeContainer`s and a
# `cdef class` -- unlike a `cdef` function or a `cdef struct` -- does not forward-declare across the
# translation unit.  The wrappers in the class call forward into this file, which is exactly the
# direction RULES.md §7.1 says is free.


cdef object _NP_EMPTY = None
cdef object _NP_ZEROS = None
cdef object _NP_ASCONTIGUOUS = None


cdef int _numpy_load() except -1:
    """Bind the three numpy names the core uses, on first use and only on first use.

    NOT A MODULE-LEVEL IMPORT, and now for two reasons rather than one.  The first is cost: the only
    things in the core that want numpy are the two matrix builders here and the fingerprint surface
    in `_fingerprints.pxi`, and paying numpy's import on every `import chython.core` to serve methods
    an ML data loader calls and a SMILES round trip never does is the wrong trade.  After the first
    call this costs one `is None` test.

    The second is that NUMPY IS AN OPTIONAL DEPENDENCY -- `chython[ml]` -- so on a minimal install
    this import is expected to fail, and it must fail at the call and not at `import chython`.  One
    module-level `import numpy` anywhere the façade imports eagerly is enough to break that, so the
    claim is not made in a comment: `chython/test/test_optional_numpy.py` asserts it on every run.

    ONE BINDER FOR THE WHOLE CORE.  A second binder in `_fingerprints.pxi`, the other numpy consumer,
    would be a second answer to "has numpy been imported yet",
    and the two would drift the moment one of them learned a fourth name.  It is also the one place
    that has to spell the install hint, which is the other half of why there is only one -- and
    `require_numpy` below extends that to the layers above rather than letting them spell it again.

    DECLARED BEFORE THEY ARE IMPORTED, the same way `_core.pyx` declares `warn`: this tree treats
    `implicit declaration of` as a build failure (RULES.md §9.7) and a bare `from numpy import ...`
    binds names Cython never saw declared.
    """
    cdef object empty
    cdef object zeros
    cdef object ascontiguousarray
    global _NP_EMPTY, _NP_ZEROS, _NP_ASCONTIGUOUS
    if _NP_EMPTY is None:
        try:
            from numpy import ascontiguousarray, empty, zeros
        except ImportError:
            # `from None`: the caller wants to know which feature needs numpy and how to get it, and
            # a chained "No module named 'numpy'" underneath that only buries the answer.
            #
            # THE LIST NAMES CATEGORIES WHERE IT CAN, and that is deliberate: a count of the methods
            # that reach here is wrong the first time someone adds one without rereading it, so the
            # message carries none.  Two of the categories surprise people, so those are spelled out
            # rather than left to "and related": every `morgan_*`/`linear_*` spelling needs numpy
            # including the ones that return a dict or a set, because they all go through the same
            # invariant vector; and the distance-derived descriptors need it because `distance_matrix`
            # is the only shortest-path code in the core and they all read its output.
            raise ImportError(
                'numpy is required and is not installed.  It is an optional dependency: install '
                '`chython[ml]`.  What needs it: every fingerprint method -- including the '
                '`*_hash_set`/`*_bit_set`/`*_hash_counts` spellings that answer a set or a dict, '
                'since all of them build the same invariant vector first -- `atom_invariants`, '
                '`adjacency_matrix`, `distance_matrix`, the descriptors derived from that matrix '
                '(`eccentricities`, `wiener_index`, `graph_radius`, `graph_diameter`, `balaban_j`), '
                '`pharmacophore_invariants`, and the four ML views -- `state_view` and '
                '`transition_view` on a molecule, and `transition_view` and `modeling_view` on a '
                'reaction.  `modeling_view` answers dicts rather than arrays, but it is a dict '
                'assembly over the same invariant arrays and needs numpy for exactly the same reason '
                'the `*_hash_set` spellings do.  Everything else works without it: readers and '
                'writers, standardize/kekule/thiele/canonicalize, stereo, depiction, the reactor, '
                'and the descriptors that are not distance-based (`tpsa`, `crippen_logp`, `bertz_ct`, '
                '`randic_index`, ring and atom counts).') from None
        _NP_EMPTY = empty
        _NP_ZEROS = zeros
        _NP_ASCONTIGUOUS = ascontiguousarray
    return 0


def require_numpy():
    """Bind numpy or raise the core's ImportError.  The layers above call this; the core does not.

    `chemistry/_pharmacophore.py` answers a numpy array and so has to fail the same way every
    numpy-backed method in the core does -- naming `chython[ml]` and what it is for.  It cannot call
    `_numpy_load` (a `cdef` function is not a module attribute), and a second copy of that message in
    `chemistry` would be a second thing to update, which is how the two would come to disagree about
    which extra to install.  So the message stays in exactly one place and this is the door to it.

    A `def` and not a `cpdef`: nothing in the core calls it, so there is no C signature worth having.
    """
    _numpy_load()


# --- connected components as molecules -----------------------------------------------------------

cdef list mol_split(MoleculeContainer mol):
    """One molecule per connected component, in first-seen atom order.  See `split`."""
    mol._require_clean()
    cdef list comps = mol.connected_components
    cdef Py_ssize_t n_comps = len(comps)
    if n_comps == 0:
        return []
    if n_comps == 1:
        # A LIST OF ONE, and a copy rather than `self`: a caller that edits `mol.split()[0]` must not
        # be editing the molecule it split.  `copy()` is O(1) here.
        return [mol.copy()]
    cdef set all_numbers = set(mol._numbers)
    cdef list out = []
    cdef MoleculeContainer part
    cdef object keep, n
    for keep in comps:
        part = mol.copy()
        with part.edit():
            for n in all_numbers - set(keep):
                part.delete_atom(n)
        out.append(part)
    return out


# --- the radius-bounded environment of a set of atoms --------------------------------------------

cdef list mol_augmented_levels(MoleculeContainer mol, object atoms, int deep):
    """The seed's atom numbers, then the seed plus each successive bond shell, out to `deep`.

    `deep` COUNTS BONDS, not atoms: level `d` is every atom within `d` bonds of the nearest seed
    atom, so level 0 is the seed itself and level 1 the seed plus its direct neighbours.  Levels are
    CUMULATIVE -- each contains the previous -- which is what makes `[-1]` the answer
    `augmented_substructure` wants and the whole list the answer `augmented_substructures` wants.

    THE LIST IS SHORTER THAN `deep + 1` WHEN THE SEED'S COMPONENTS RUN OUT.  The walk stops as soon
    as a shell adds nothing, so `deep=99` on ethanol gives three levels and not a hundred; BFS
    distances are contiguous, so an empty shell is the end and not a gap.  A caller must therefore
    read `len()` rather than trusting `deep`.

    One BFS from the whole seed at once, not one per seed atom: the distance that bounds a shell is
    the distance to the NEAREST seed atom, which is exactly what a multi-source BFS computes.
    """
    mol._require_clean()
    if deep < 0:
        raise ValueError('deep must be a non-negative number of bonds')
    cdef set seed = set()
    cdef object x
    for x in atoms:
        if x not in mol._index_of:
            raise KeyError(x)
        seed.add(x)
    if not seed:
        # the same refusal `substructure` makes, made here so that the message names this method
        raise ValueError('an augmented substructure of no atoms is not a molecule')

    cdef Structure structure = mol._structure
    cdef uint32_t n_atoms = structure.header.atom_count
    cdef uint32_t *ptr = csr_ptr(structure)
    cdef halfedge_t *edges = csr_edges(structure)
    cdef uint32_t NONE = 0xffffffff
    cdef uint32_t i, k, v, to, head = 0, fill = 0
    cdef uint32_t bound = <uint32_t> deep
    cdef list numbers = mol._numbers
    cdef list shells = []
    cdef list levels = []
    cdef list acc = []
    cdef uint32_t d

    # one block, two regions (RULES.md §5.2): the distances and the BFS queue.  Every atom is
    # queued at most once, so `n_atoms` slots is the exact worst case for the queue.
    cdef uint32_t *block = <uint32_t *> PyMem_Malloc(<size_t> 2 * n_atoms * sizeof(uint32_t))
    if block is NULL:
        raise MemoryError('augmented substructure scratch allocation failed')
    cdef uint32_t *dist = block
    cdef uint32_t *queue = block + n_atoms
    try:
        for i in range(n_atoms):
            dist[i] = NONE
        for x in seed:
            i = <uint32_t> mol._index_of[x]
            dist[i] = 0
            queue[fill] = i
            fill += 1
        with nogil:
            while head < fill:
                v = queue[head]
                head += 1
                if dist[v] == bound:
                    continue
                for k in range(ptr[v], ptr[v + 1]):
                    to = edges[k].to
                    if dist[to] == NONE:
                        dist[to] = dist[v] + 1
                        queue[fill] = to
                        fill += 1
        for d in range(bound + 1):
            shells.append([])
        for i in range(n_atoms):
            if dist[i] != NONE:
                (<list> shells[dist[i]]).append(numbers[i])
    finally:
        PyMem_Free(block)

    for d in range(bound + 1):
        if d and not <list> shells[d]:
            break
        acc = acc + <list> shells[d]
        levels.append(acc)
    return levels


# --- the matrix forms ----------------------------------------------------------------------------
#
# BOTH ARE INDEXED BY POSITION, not by atom number, and the position is the molecule's own atom
# order -- row `i` is `mol.atom_numbers[i]`, and a caller that needs to go back does it through
# `index_of`.  A matrix keyed by stable id is not expressible: ids are sparse after any deletion.
#
# ORDER 8 IS A BOND HERE.  Ring perception excludes the dative bond because a dative bond is not a
# ring closure; a WALK has no such argument, so a metal complex is one component for both of these
# and for `split`, exactly as `label_components` already has it.

cdef void csr_bfs_all(const uint32_t *ptr, const uint32_t *to, uint32_t n_atoms,
                      int32_t *out, uint32_t stride, uint32_t *queue) noexcept nogil:
    """One BFS per source over a plain CSR: 0 on the diagonal, -1 for a disconnected pair.

    A PLAIN `to` ARRAY AND NOT `halfedge_t`, so the reaction union can call it.  The union graph
    carries an order per side and has no half-edge to walk; a second BFS beside this one is the
    alternative, and two shortest-path implementations in one core is the thing worth avoiding.
    Materializing `to` costs 4 bytes per half-edge -- 224 B for a 27-atom molecule.

    `stride` IS THE ROW PITCH AND NOT THE ATOM COUNT.  `mol_distance_matrix` passes `n_atoms`, and the
    ML views pass their padded `width`, so the block can be the top-left corner of a bigger buffer.

    ALLOCATES NOTHING.  Both callers already own a scratch block, and `noexcept nogil` needs no
    failure path.  `out` needs `n_atoms * stride` int32; `queue` needs `n_atoms` uint32, which is the
    exact worst case because an atom is queued once per source.
    """
    cdef uint32_t i, k, v, w, head, fill, src
    cdef int32_t d
    cdef int32_t *row
    for src in range(n_atoms):
        row = out + <size_t> src * stride
        for i in range(n_atoms):
            row[i] = -1
        row[src] = 0
        queue[0] = src
        head = 0
        fill = 1
        while head < fill:
            v = queue[head]
            head += 1
            d = row[v] + 1
            for k in range(ptr[v], ptr[v + 1]):
                w = to[k]
                if row[w] < 0:
                    row[w] = d
                    queue[fill] = w
                    fill += 1


cdef object mol_adjacency_matrix(Structure structure, bint set_bonds):
    """`(n, n)` uint32: 1 where a bond exists, or its stored order when `set_bonds`.  Symmetric."""
    _numpy_load()
    cdef uint32_t n_atoms = structure.header.atom_count
    cdef object out = _NP_ZEROS((n_atoms, n_atoms), dtype='uint32')
    if n_atoms == 0:
        return out
    cdef uint32_t[:, ::1] adj = out
    cdef uint32_t *ptr = csr_ptr(structure)
    cdef halfedge_t *edges = csr_edges(structure)
    cdef uint32_t i, k
    with nogil:
        for i in range(n_atoms):
            for k in range(ptr[i], ptr[i + 1]):
                adj[i, edges[k].to] = edges[k].order if set_bonds else 1
    return out


cdef object mol_distance_matrix(Structure structure):
    """`(n, n)` int32 of topological distances: 0 on the diagonal, -1 for a disconnected pair.

    One BFS per atom, which is the right algorithm at this size and the reason the result is exact
    rather than a truncation: nothing here takes a cutoff, so a caller that wants one clamps.

    WHY -1 FOR "NO PATH", and it is the consumer that settles it rather than taste.  chytorch's
    `graph_distances` (`chytorch/utils/data/_utils.py`) documents the encoding it feeds a model:
    after `+ 2`, "1 marks a pair in different components, 2 an atom with itself, 3 neighbours", and
    0 is left free for padding.  -1 shifts to exactly its 1, 0 to its 2, and 1 to its 3 -- so the
    whole of that function becomes `mol.distance_matrix() + 2` under a `minimum`, and the
    `nan_to_num(posinf=1)` it needs today goes away with the float matrix that produced the inf.
    The two alternatives both lose: 0 collides with the diagonal, and a large sentinel is
    indistinguishable from "far but connected" after the consumer's own clamp.

    int32 for the same reason: the consumer's last line is `IntTensor(d)`, and `IntTensor` is int32.
    """
    _numpy_load()
    cdef uint32_t n_atoms = structure.header.atom_count
    cdef object out = _NP_EMPTY((n_atoms, n_atoms), dtype='int32')
    if n_atoms == 0:
        return out
    cdef int32_t[:, ::1] dist = out
    cdef uint32_t *ptr = csr_ptr(structure)
    cdef halfedge_t *edges = csr_edges(structure)
    cdef uint32_t n_half = ptr[n_atoms]
    cdef uint32_t k

    # one block, two regions (RULES.md §5.2): the half-edge targets flattened for `csr_bfs_all`, and
    # the BFS queue it reuses for every source.
    cdef uint32_t *block = <uint32_t *> PyMem_Malloc(<size_t> (n_half + n_atoms) * sizeof(uint32_t))
    if block is NULL:
        raise MemoryError('distance matrix scratch allocation failed')
    cdef uint32_t *to = block
    cdef uint32_t *queue = block + n_half
    try:
        with nogil:
            for k in range(n_half):
                to[k] = edges[k].to
            csr_bfs_all(ptr, to, n_atoms, &dist[0, 0], n_atoms, queue)
    finally:
        PyMem_Free(block)
    return out
