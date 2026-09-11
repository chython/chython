# cython: freethreading_compatible=True
# cython: undeclared_check_usage=error
# cython: warn.undeclared=True
# cython: warn.unused=True
# cython: warn.unused_arg=True
# cython: warn.maybe_uninitialized=True
# cython: boundscheck=False
# cython: wraparound=False
# cython: auto_pickle=False
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
"""
The core: one extension module, one translation unit.

Cython compiles one module per .pyx, so the twenty-six fragments below are .pxi files that this
file textually includes rather than twenty-six modules that cimport each other. The fragments are
still ordered -- the include order is the dependency order, and nothing in an earlier
fragment may reach into a later one except by forward reference -- but the C compiler sees
a single unit. That is the whole point: a cross-module `cdef` call goes through the
__pyx_capi__ function pointer table, so `perceive_rings`, `fill_features`, `csr_build` and
friends were indirect calls that no optimizer could see through. Merged, every one of them
is `static` and inlinable at the call site.

Fragments, in include order:

  _molecule_arena     the arena: header, segments, 24-byte atom record, half-edge CSR, accessors,
                      and the field domains every validator cites (bond orders, map number, etc.)
  _elements           the element symbol table and the MDL common-isotope table
  _valence            the query side of the valence rule collection, whose data lives in
                      valence_rules.tsv: how many implicit hydrogens an atom gets, and a
                      three-state verdict on a state it is asked about. The CHEMISTRY model --
                      the SMILES layer has its own model of the same shape answering a question
                      about the notation instead, and the two must not be merged
  _features           the derived layer: the 4x u64 feature words and the element index
  _rings              bridge detection and Vismara relevant-cycle perception
  _sssr               the dict-graph adapter onto that perception
  _morgan             refinement of the graph to an equitable partition: symmetry classes
  _canonical          the molecule's symmetry orbits, from an exhaustive search per unresolved
                      pair of atoms (mol_automorphisms), and the canonical atom order, from the
                      extremal labelling of the refinement tree (mol_canonical_order)
  _stereo             the derived stereo unit table (SEG_STEREO_UNIT) and the perception that
                      fills it: which atoms and bonds could carry a configuration
  _inchi              the libinchi bridge: the vendored structs, the loader, and the two directions
                      (molecule to InChI and InChIKey, InChI back to a molecule)
  _query_arena        the query arena: enums, structs, Query class, query_alloc, segment accessors
  _query_boxes        primitive constants, wbox_t, span table, and logic compiler
  _query_seal         journal-to-sealed-arena builder: _compute_automorphisms, query_seal
  _isomorphism        the subgraph isomorphism kernel: a resumable DFS over the query's plan
  _molecule_container MoleculeContainer -- the Python-facing molecule and its builder journal
  _molecule_topology  the topology surface behind four of the container's methods: connected
                      components as molecules, the radius-bounded environment of a set of atoms,
                      and the adjacency and distance matrices.  AFTER the container, because
                      `mol_split` builds one and a `cdef class` does not forward-declare
  _ml                 `TensorEncoding`, `mol_state_view`, `mol_transition_view`,
                      `reaction_transition_view`.  AFTER the topology fragment, because
                      `mol_state_view` calls `csr_bfs_all` from there
  _fingerprints       the fingerprint surface: per-atom labels, circular (Morgan/ECFP) and linear
                      (path) enumerators, and the three fold shapes (bit set, binary, counted).
                      AFTER the container and the topology fragment, whose numpy binder it reads
  _descriptors        the graph descriptors: composition counts, ring system counts and the classical
                      topological indices.  Arithmetic over what the arena already derived -- the
                      element index, the minimum cycle basis and the one distance matrix -- so it
                      perceives nothing.  AFTER the container and the topology fragment, whose
                      distance matrix and numpy binder it reads
  _molecule_views     Atom and Bond: borrowed, generation-checked handles into the arena
  _smiles_write       the SMILES writer: the canonical traversal, the atom/bond tokens, the
                      valence model that decides when brackets may be dropped, and the CXSMILES tail
  _query_container    QueryContainer: the journal surface for building a query
  _kekule             kekulisation: a stated aromatic edge set becomes bond orders 1 and 2
  _thiele             aromatisation: Kekule bond orders become order 4, the inverse direction, and
                      the caller of _kekule's one-atom classification table
  _smiles_read        the SMILES reader: the tokeniser, the parse graph, the notation-model
                      hydrogen counts, and atom-case aromatic promotion. AFTER _kekule, whose
                      atom classifier decides the hydrogen count of an aromatic atom and whose
                      kekuliser answers the one question promotion has to ask
  _smarts_read        the SMARTS reader: the tokeniser and the bracket-body grammar, emitting journal
                      ops into a QueryContainer. AFTER _smiles_read, whose element table, digit
                      scanners and CXSMILES field scanner it shares rather than restating them
  _smirks_read        the SMIRKS reader: the arrow, both sides through the SMARTS lexer next door, and
                      everything only a two-sided string can state -- which atoms pair, which are
                      deleted, which are created, and what each product primitive builds or checks
  _smirks_patch       the patcher: a ReactionTemplate applied to molecules on the container's edit
                      session, yielding ReactionContainers. AFTER _smirks_read, whose template it
                      consumes, and the only fragment that drives BOTH containers at once
  _pach               the legacy pach codec, at the end because it serialises the finished arena
                      These last six are last because they drive the containers rather than being
                      driven by them -- the only fragments here that are clients of
                      MoleculeContainer's and QueryContainer's public surface rather than layers
                      underneath them
"""
cimport cython
from cpython.exc cimport PyErr_CheckSignals
from cpython.mem cimport PyMem_Free, PyMem_Malloc, PyMem_Realloc
from libc.math cimport NAN, isnan, log2, round, sqrt
from libc.stdint cimport (int8_t, int32_t, int64_t, uint8_t, uint16_t, uint32_t, uint64_t,
                          uintptr_t)
from libc.stdlib cimport calloc, free, malloc, realloc
from libc.string cimport memcmp, memcpy, memset
from libc.time cimport time, time_t

# `round` ABOVE REBINDS THE PYTHON BUILTIN FOR THE WHOLE CORE, and it has to: `include` is textual, so
# one translation unit sees one `round`, and `_pach_f16_encode` is `noexcept nogil` -- the builtin
# cannot be called there at all.  The consequence to know before reading a call site: every `round()`
# in this core is C's, which breaks a tie away from zero where Python's breaks it to even.  The seven
# call sites all round a coordinate already scaled by XY_SCALE, where the two agree on every value an
# MDL file can state (four decimals, so no tie), and a display coordinate has no stake in the
# difference regardless.  A future cimport that shadows a builtin whose semantics DO matter gets an
# `as c_name` alias instead.
#
# THE CORE MAKES NO PYTHON-LEVEL IMPORT.  `warnings.warn` was the one, and nothing warns here now:
# an alternative spelling of a live name is not a deprecation.  A future import at this scope needs
# `cdef object <name>` above it, because `warn.undeclared` is on and this tree treats a Cython
# warning as a build failure.

include "_molecule_arena.pxi"
include "_elements.pxi"
include "_valence.pxi"
include "_features.pxi"
include "_rings.pxi"
include "_sssr.pxi"
include "_morgan.pxi"
include "_canonical.pxi"
include "_stereo.pxi"
include "_inchi.pxi"
include "_query_arena.pxi"
include "_query_boxes.pxi"
include "_query_seal.pxi"
include "_isomorphism.pxi"
include "_molecule_container.pxi"
include "_molecule_topology.pxi"
include "_ml.pxi"
include "_fingerprints.pxi"
include "_descriptors.pxi"
include "_molecule_views.pxi"
include "_smiles_write.pxi"
include "_query_container.pxi"
include "_kekule.pxi"
include "_thiele.pxi"
# AFTER `_kekule.pxi`, whose `arom_classify_atom` it calls, and `_valence.pxi`, whose
# `val_implicit_h` it calls; BEFORE `_smiles_read.pxi`, whose `smi_read_h` calls IT.  One translation
# unit, so a `cdef` function has to be declared above its caller.  `kekule()`'s own call back into
# this layer is not a cycle: it is a `def`, so the name is a module-globals lookup at run time.
include "_hydrogens.pxi"
include "_smiles_read.pxi"
include "_smarts_read.pxi"
include "_smirks_read.pxi"
include "_smirks_patch.pxi"
include "_pach.pxi"
include "_pach3.pxi"
