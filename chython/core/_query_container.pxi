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
# QueryContainer: the journal surface for building a query.
#
# What this is NOT: there is no property API (no atom.charge getters), no iteration protocol,
# no pretty-printing, no serialisation.  A query is built and matched; everything else is the
# molecule's job.  Task 13 appends inspection methods in the same style as closure_count()
# (Task 11) and automorphism_count() / automorphism_generation() (Task 12).


DEF QUERY_JOURNAL_MIN_CAP = 16


cdef class _QueryEditScope:
    """Context manager returned by QueryContainer.edit()."""
    cdef QueryContainer _q
    cdef uint32_t _saved_len
    cdef uint32_t _saved_next_id

    def __cinit__(self, QueryContainer q not None):
        self._q = q
        self._saved_len = 0
        self._saved_next_id = 1

    def __enter__(self):
        self._saved_len = self._q._journal_len
        self._saved_next_id = self._q._next_id
        self._q._scope_depth += 1
        return self._q

    @cython.warn.unused_arg(False)
    def __exit__(self, exc_type, exc_val, exc_tb):
        self._q._scope_depth -= 1
        if self._q._scope_depth == 0 and exc_type is not None:
            self._q._journal_len = self._saved_len
            self._q._next_id = self._saved_next_id
            self._q._invalidate()
        return False


cdef class QueryContainer:
    """A journal of construction ops that compiles to a sealed Query on first use.

    The sealed Query is cached; any mutation invalidates it so the next sealed() call
    recompiles from the updated journal.  Stable atom ids start at 1 and increment.
    """
    cdef Query _query
    cdef qop_t *_journal
    cdef uint32_t _journal_len
    cdef uint32_t _journal_cap
    cdef uint32_t _next_id
    cdef uint32_t *_position_to_n
    cdef int _scope_depth
    cdef uint32_t _seal_generation
    cdef uint32_t _automorphism_generation    # group computations, i.e. successful seals

    def __cinit__(self):
        self._query = None
        self._journal = NULL
        self._journal_len = 0
        self._journal_cap = 0
        self._next_id = 1
        self._position_to_n = NULL
        self._scope_depth = 0
        self._seal_generation = 0
        self._automorphism_generation = 0

    def __dealloc__(self):
        PyMem_Free(self._journal)
        self._journal = NULL
        PyMem_Free(self._position_to_n)
        self._position_to_n = NULL

    cdef int _invalidate(self) except -1:
        """Free the cached arena and its position map, forcing a reseal on next use.

        Resets _automorphism_generation, because the group lives IN the arena that just went away.
        It counts the group currently held, not the groups ever computed: a caller asking whether
        the filter it is about to use was computed for the query it is about to match wants the
        first, and _seal_generation next door already answers the second.
        """
        self._query = None
        PyMem_Free(self._position_to_n)
        self._position_to_n = NULL
        self._automorphism_generation = 0
        return 0

    cdef int _append(self, qop_t op) except -1:
        """Append one op to the journal, growing it geometrically from QUERY_JOURNAL_MIN_CAP.

        The doubling is unguarded: a uint32_t cap wraps to 0 at 2**31 ops, after which the next
        Realloc writes past a zero-sized block.  Same as _append (_molecule_container.pxi:175-176)
        -- no guard because 2**31 qop_t is 56 GiB of journal, unreachable in practice.
        """
        cdef uint32_t cap
        cdef qop_t *grown
        if self._journal_len == self._journal_cap:
            cap = QUERY_JOURNAL_MIN_CAP if self._journal_cap == 0 else self._journal_cap * 2
            grown = <qop_t *> PyMem_Realloc(self._journal, cap * sizeof(qop_t))
            if grown is NULL:
                raise MemoryError('journal reallocation failed')
            self._journal = grown
            self._journal_cap = cap
        self._journal[self._journal_len] = op
        self._journal_len += 1
        return 0

    cpdef uint32_t add_atom(self) except 0:
        """Allocate a new stable id, journal it as QOP_ADD_ATOM, and return the id."""
        cdef qop_t op
        cdef uint32_t n
        if self._next_id >= 0xFFFFFFFF:
            raise OverflowError('stable id space is exhausted; ids are never reused')
        n = self._next_id
        self._invalidate()
        memset(&op, 0, sizeof(qop_t))
        op.op = QOP_ADD_ATOM
        op.a = n
        self._append(op)
        self._next_id += 1
        return n

    cpdef int add_bond(self, uint32_t n, uint32_t m) except -1:
        """Journal a bond between two known stable ids.

        Rejects unknown ids (>= _next_id or 0) and self-loops immediately.
        """
        cdef qop_t op
        if n == m:
            raise ValueError('a self loop is not allowed')
        if n == 0 or n >= self._next_id:
            raise ValueError('unknown atom %d' % n)
        if m == 0 or m >= self._next_id:
            raise ValueError('unknown atom %d' % m)
        self._invalidate()
        memset(&op, 0, sizeof(qop_t))
        op.op = QOP_ADD_BOND
        op.a = n
        op.b = m
        self._append(op)
        return 0

    cpdef int atom_primitive(self, uint32_t n, str name, int32_t value=0,
                             bint negated=False) except -1:
        """Journal one SMARTS primitive for atom n.

        Raises KeyError for an unknown primitive name; raises ValueError for an unknown atom.
        """
        cdef qop_t op
        cdef uint32_t kind
        if n == 0 or n >= self._next_id:
            raise ValueError('unknown atom %d' % n)
        kind = <uint32_t> PRIM_NAMES[name]   # KeyError for an unknown name
        self._invalidate()
        memset(&op, 0, sizeof(qop_t))
        op.op = QOP_ATOM_TOKEN
        op.a = n
        op.opcode = OPC_PRIM
        op.kind = kind
        op.value = value
        op.negated = 1 if negated else 0
        self._append(op)
        return 0

    cpdef int atom_operator(self, uint32_t n, str name) except -1:
        """Journal a logic operator (or / and_low / and_high) for atom n.

        Raises KeyError for an unknown operator name; raises ValueError for an unknown atom.
        """
        cdef qop_t op
        cdef uint32_t opcode
        if n == 0 or n >= self._next_id:
            raise ValueError('unknown atom %d' % n)
        opcode = <uint32_t> OPC_NAMES[name]   # KeyError for an unknown name
        self._invalidate()
        memset(&op, 0, sizeof(qop_t))
        op.op = QOP_ATOM_TOKEN
        op.a = n
        op.opcode = opcode
        self._append(op)
        return 0

    cpdef int bond_primitive(self, uint32_t n, uint32_t m, str name, int32_t value=0,
                              bint negated=False) except -1:
        """Journal one SMARTS primitive for the bond between (n, m).

        Raises KeyError for an unknown primitive name.  Pair validation (that (n, m) is actually
        a bond) is deferred to seal — query_seal raises ValueError('unknown bond %d-%d') there.
        """
        cdef qop_t op
        cdef uint32_t kind
        kind = <uint32_t> PRIM_NAMES[name]   # KeyError for an unknown name
        self._invalidate()
        memset(&op, 0, sizeof(qop_t))
        op.op = QOP_BOND_TOKEN
        op.a = n
        op.b = m
        op.opcode = OPC_PRIM
        op.kind = kind
        op.value = value
        op.negated = 1 if negated else 0
        self._append(op)
        return 0

    cpdef int bond_operator(self, uint32_t n, uint32_t m, str name) except -1:
        """Journal a logic operator for the bond between (n, m).

        Raises KeyError for an unknown operator name.
        """
        cdef qop_t op
        cdef uint32_t opcode
        opcode = <uint32_t> OPC_NAMES[name]   # KeyError for an unknown name
        self._invalidate()
        memset(&op, 0, sizeof(qop_t))
        op.op = QOP_BOND_TOKEN
        op.a = n
        op.b = m
        op.opcode = opcode
        self._append(op)
        return 0

    cpdef int set_group(self, uint32_t n, int32_t group) except -1:
        """Journal a reaction-group assignment for atom n."""
        cdef qop_t op
        if n == 0 or n >= self._next_id:
            raise ValueError('unknown atom %d' % n)
        self._invalidate()
        memset(&op, 0, sizeof(qop_t))
        op.op = QOP_SET_GROUP
        op.a = n
        op.value = group
        self._append(op)
        return 0

    cpdef int set_masked(self, uint32_t n) except -1:
        """Journal a masked flag for atom n (protects it from deletion in Reactor)."""
        cdef qop_t op
        if n == 0 or n >= self._next_id:
            raise ValueError('unknown atom %d' % n)
        self._invalidate()
        memset(&op, 0, sizeof(qop_t))
        op.op = QOP_SET_MASKED
        op.a = n
        self._append(op)
        return 0

    cpdef int set_map_number(self, uint32_t n, int number) except -1:
        """Journal an atom-map number for atom n."""
        cdef qop_t op
        if n == 0 or n >= self._next_id:
            raise ValueError('unknown atom %d' % n)
        if number < 0 or number > MAP_NUMBER_MAX:
            raise ValueError('map number %d is out of range 0..%d' % (number, MAP_NUMBER_MAX))
        self._invalidate()
        memset(&op, 0, sizeof(qop_t))
        op.op = QOP_SET_MAP
        op.a = n
        op.value = <int32_t> number
        self._append(op)
        return 0

    cpdef int set_stereo_group(self, uint32_t n, int kind, int group) except -1:
        """Journal an enhanced-stereo group membership for atom n -- `&<group>` or `o<group>`.

        A LABEL, in the manner of `set_masked` and `set_map_number`, and not a primitive: it compiles
        to no box because there is nothing in a target it could be compared with.  `sealed()` refuses
        it outright for that reason, so this is only ever read back off the journal, by the SMIRKS
        reader for the one side that never seals.

        `kind` is 2 for OR and 3 for AND -- `MoleculeContainer.set_stereo_group`'s numbering, and its
        1..63 group range, because that is the field the number ends up in.  Absolute has no bracket
        spelling and none is accepted here: `a` inside a bracket is the aromatic flag.
        """
        cdef qop_t op
        if n == 0 or n >= self._next_id:
            raise ValueError('unknown atom %d' % n)
        if kind != 2 and kind != 3:
            raise ValueError('kind must be 2 or (`o<n>`) or 3 and (`&<n>`)')
        if group < 1 or group > 63:
            raise ValueError('OR and AND groups must have a group id in 1..63')
        self._invalidate()
        memset(&op, 0, sizeof(qop_t))
        op.op = QOP_SET_STEREO_GROUP
        op.a = n
        op.kind = <uint32_t> kind
        op.value = <int32_t> group
        self._append(op)
        return 0

    cpdef int set_bond_direction(self, uint32_t n, uint32_t m, int direction) except -1:
        """Journal a `/` or `\\` on the bond between (n, m), oriented from `n` to `m`.

        A LABEL like `set_stereo_group`, journalled and never compiled: `sealed()` refuses it, so it is
        only ever read back off the journal by the SMIRKS reader for the side that never seals.  Which
        end it is written from is not part of the notation -- the same statement read from the other
        atom is the other direction -- and this records the orientation rather than normalising it,
        because the atom written first is what the string's reader knows.

        `direction` is SMI_DIR_UP for `/` and SMI_DIR_DOWN for `\\`.  Pair validation is the seal's.
        """
        cdef qop_t op
        if n == 0 or n >= self._next_id or m == 0 or m >= self._next_id:
            raise ValueError('unknown bond %d-%d' % (n, m))
        if direction != SMI_DIR_UP and direction != SMI_DIR_DOWN:
            raise ValueError('direction must be 1 up (`/`) or 2 down (`\\`)')
        self._invalidate()
        memset(&op, 0, sizeof(qop_t))
        op.op = QOP_SET_BOND_DIRECTION
        op.a = n
        op.b = m
        op.value = <int32_t> direction
        self._append(op)
        return 0

    cdef Query sealed(self):
        """Seal the journal into a query arena on first use, then cache the result.

        Raises ValueError for an empty query or a contradictory one; raises
        NotImplementedError for a stereo primitive (not yet supported).  Does not catch
        either -- let them propagate to the caller.

        The returned Query's lifetime is independent of this container: _invalidate() drops
        the container's reference but a caller holding the old Query object keeps that arena
        alive.  Tasks 10-14 must re-fetch via sealed() after any mutation rather than caching
        the result across one.
        """
        if self._query is not None:
            return self._query
        if self._next_id == 1:
            raise ValueError('an empty query matches nothing')
        self._query = query_seal(self._journal, self._journal_len, self._next_id,
                                 &self._position_to_n)
        self._seal_generation += 1
        # query_seal computes the automorphism group inline, so one seal is one group computation.
        # Bumped after the call, not before: a seal that raises leaves no arena and no group.
        self._automorphism_generation += 1
        return self._query

    # -----------------------------------------------------------------------
    # Journal-level read-only properties — never seal
    # -----------------------------------------------------------------------

    @property
    def atom_count(self):
        """Number of atoms currently in the journal (mid-construction safe)."""
        cdef uint32_t i
        cdef uint32_t count
        count = 0
        for i in range(self._journal_len):
            if self._journal[i].op == QOP_ADD_ATOM:
                count += 1
        return count

    @property
    def bond_count(self):
        """Number of bonds currently in the journal (mid-construction safe).

        Counts journal ops, so a duplicate bond (which seal rejects) is counted twice; the
        sealed arena's bond_count may therefore disagree before the error is caught.
        """
        cdef uint32_t i
        cdef uint32_t count
        count = 0
        for i in range(self._journal_len):
            if self._journal[i].op == QOP_ADD_BOND:
                count += 1
        return count

    def __len__(self):
        return self.atom_count

    def edit(self):
        """Return a context manager that records the journal checkpoint and rolls back on error."""
        return _QueryEditScope(self)

    # -----------------------------------------------------------------------
    # Inspection surface — every method here forces a seal
    # -----------------------------------------------------------------------

    cpdef uint32_t atom_count_sealed(self) except 0:
        """Atom count from the sealed arena.  Raises ValueError for an empty query."""
        return self.sealed().header.atom_count

    cpdef list box_counts(self):
        """Box count per DFS position.  Forces a seal."""
        cdef Query q = self.sealed()
        cdef qatom_t *atoms = q.atoms()
        cdef list out = []
        cdef uint32_t i
        for i in range(q.header.atom_count):
            out.append(atoms[i].box_count)
        return out

    cpdef uint32_t closure_count(self) except? 0:
        """Number of closure bonds in the sealed query.  Forces a seal.

        A closure is a query bond both of whose endpoints are already placed when the DFS
        reaches it (Task 11).  Every acyclic query has zero closures, so `except? 0` is
        required -- `except 0` would misread a legitimate zero return as an exception.
        """
        cdef Query q = self.sealed()
        return <uint32_t> (q.header.segments[QSEG_CLOSURES].length // sizeof(qclosure_t))

    cpdef uint32_t automorphism_count(self) except? 0:
        """Stored automorphisms of the sealed query, the identity excluded.  Forces a seal.

        `except? 0` rather than `except 0`: an asymmetric query legitimately returns 0, and that is
        the common case, so a 0 must not be read as an exception.
        """
        cdef Query q = self.sealed()
        return q.header.automorphism_count

    cpdef uint32_t automorphism_generation(self):
        """How many times the automorphism group was computed for the group currently held.

        1 once the query is sealed, 0 before that and again after any mutation -- the group lives
        in the arena, so it goes away with it.  Deliberately does NOT seal: a caller asking
        whether the group has been computed must not cause it to be computed.
        """
        return self._automorphism_generation

    cpdef uint32_t seal_generation(self):
        """How many times sealed() actually built a new arena (0 before first seal)."""
        return self._seal_generation

    cpdef dict map_numbers(self):
        """{stable id: map number} for every atom whose map number is non-zero.  Forces a seal."""
        cdef Query q = self.sealed()
        cdef qatom_t *atoms = q.atoms()
        cdef dict out = {}
        cdef uint32_t i, n
        for i in range(q.header.atom_count):
            if atoms[i].map_number != 0:
                n = self._position_to_n[i]
                out[n] = <int> atoms[i].map_number
        return out

    cpdef tuple component_groups(self):
        """One entry per query component, in arena order: (frozenset of stable ids, group or None).

        Reports the group as the SEAL reduced it -- per component, not per atom -- which is the form
        the matcher constrains on: same group means one molecule component, different groups mean
        different ones, `None` means unconstrained.  Forces a seal, so a component that spans two
        groups raises here rather than reading back as one of them.
        """
        cdef Query q = self.sealed()
        cdef qcomp_t *comps = query_components(q)
        cdef list out = []
        cdef set ids
        cdef uint32_t i, s
        for i in range(q.header.component_count):
            ids = set()
            for s in range(comps[i].begin, comps[i].end):
                ids.add(self._position_to_n[s])
            out.append((frozenset(ids), None if comps[i].group < 0 else <int> comps[i].group))
        return tuple(out)

    cpdef dict wildcard_atoms(self):
        """{stable id: 'any' | 'metal'} for every atom that names no single element.  Forces a seal.

        `'any'` is an atom whose element span no box touched -- `[A]`, `[*]`, and a bracket that only
        counts, like `[D2]`.  `'metal'` is `[M]`: every box constrained to exactly the 93 metals.  An
        atom that names an element, a list of them, or a negation is absent from the dict.

        WHAT THIS IS FOR.  A rule table needs to know which atoms of a pattern are shared CONTEXT
        rather than the site being repaired, because two matches of one rule may legally overlap on
        context and must not overlap on the site.  The answer is read off the BOXES and not off a
        symbol, which is why `[A]` and `[*]` answer alike and a hand-written list of every metal does
        not answer `'metal'`.

        `'any'` deliberately does not distinguish `[A]` from `[*]`.  The two differ in the charge and
        radical spans, not the element one, and no caller has ever asked which token was typed --
        `*` is the WITHDRAWAL of the charge default, so asking the question in element terms would
        answer it wrong.
        """
        cdef Query q = self.sealed()
        cdef qatom_t *atoms = q.atoms()
        cdef dict out = {}
        cdef uint32_t i, n
        for i in range(q.header.atom_count):
            if atoms[i].flags & QATOM_ANY_ELEMENT:
                out[self._position_to_n[i]] = 'any'
            elif atoms[i].flags & QATOM_METAL_ELEMENT:
                out[self._position_to_n[i]] = 'metal'
        return out

    cpdef frozenset masked_atoms(self):
        """Stable ids of every masked atom.  Forces a seal."""
        cdef Query q = self.sealed()
        cdef qatom_t *atoms = q.atoms()
        cdef set out = set()
        cdef uint32_t i, n
        for i in range(q.header.atom_count):
            if atoms[i].flags & QATOM_MASKED:
                n = self._position_to_n[i]
                out.add(n)
        return frozenset(out)

    # -----------------------------------------------------------------------
    # Matching — every method here seals and then runs the kernel in _isomorphism.pxi
    # -----------------------------------------------------------------------

    cdef tuple _query_numbers(self, Query query):
        """Snapshot position -> query atom number as a tuple, given an already-sealed query.

        Holding a reference to this container does NOT protect _position_to_n: a mutation
        calls _invalidate(), which frees that array while leaving the container alive.  A
        generator must therefore copy it once up front rather than read it per match.
        """
        cdef uint32_t i
        cdef list out = []
        for i in range(query.header.atom_count):
            out.append(self._position_to_n[i])
        return tuple(out)

    cpdef tuple query_numbers(self):
        """Query atom numbers in DFS position order: slot i of a get_raw_mapping tuple is this atom.

        Forces a seal.  Without this the tuples get_raw_mapping returns are uninterpretable, which
        defeats that method's purpose: position order is query_seal's DFS order, NOT declaration
        order -- a rare element roots the DFS, so a query declared C-O-C hands back the OXYGEN in
        slot 0 -- and every other accessor here is keyed by atom number.
        """
        return self._query_numbers(self.sealed())

    def get_mapping(self, MoleculeContainer molecule not None, bint automorphism_filter=False):
        """Yield every embedding as {query stable id: molecule stable id}.

        With automorphism_filter, one embedding per automorphism orbit of the QUERY: a symmetric
        query stops reporting the same set of molecule atoms once per symmetry.
        """
        cdef Query query = self.sealed()
        cdef tuple query_ids = self._query_numbers(query)
        cdef Structure structure
        cdef list molecule_ids
        cdef matcher_t mt
        cdef uint32_t i, n_atoms
        cdef bint found
        cdef dict out
        molecule._require_clean()
        # Both of these are locals on purpose: the generator's frame keeps the arena and the
        # molecule's index -> stable id list alive for the whole search, so an edit to either
        # container mid-iteration cannot pull the buffers out from under the matcher.  Keeping the
        # object alive is not enough on its own -- see matcher_reseat after the yield.
        structure = molecule._structure
        molecule_ids = molecule._numbers
        n_atoms = query.header.atom_count
        matcher_init(&mt, query, structure, automorphism_filter)
        try:
            while True:
                with nogil:
                    found = matcher_next(&mt)
                if not found:
                    break
                out = {}
                for i in range(n_atoms):
                    out[query_ids[i]] = molecule_ids[mt.mapping[i]]
                yield out
                # the loop body just ran arbitrary Python, which may have appended a lazy derived
                # segment to this arena and moved the buffer the matcher points into
                matcher_reseat(&mt, structure)
        finally:
            matcher_free(&mt)

    def get_raw_mapping(self, MoleculeContainer molecule not None, bint automorphism_filter=False):
        """Yield each embedding as a tuple of molecule stable ids, in query DFS order.

        The query side is not translated and no dict is built, which is what a caller doing its
        own bookkeeping wants: the per-match dict dominates the cost of a fast search.  Stable
        ids rather than arena indices, because an index is invalidated by the next mutation while
        a stable id is the molecule's public identity, and one is an array lookup from the other.
        """
        cdef Query query = self.sealed()
        cdef Structure structure
        cdef list molecule_ids
        cdef matcher_t mt
        cdef uint32_t i, n_atoms
        cdef bint found
        cdef list row
        molecule._require_clean()
        structure = molecule._structure
        molecule_ids = molecule._numbers
        n_atoms = query.header.atom_count
        matcher_init(&mt, query, structure, automorphism_filter)
        try:
            while True:
                with nogil:
                    found = matcher_next(&mt)
                if not found:
                    break
                row = []
                for i in range(n_atoms):
                    row.append(molecule_ids[mt.mapping[i]])
                yield tuple(row)
                matcher_reseat(&mt, structure)   # see get_mapping: the yield can move the arena
        finally:
            matcher_free(&mt)

    def may_match(self, MoleculeContainer molecule not None):
        """Return False only when no embedding of this query into the molecule is possible.

        A sound lower bound: True means 'maybe', False means 'definitely not'.
        Forces a seal.  Returns a Python bool so that callers may use `is False` / `is True`.
        """
        cdef Query query = self.sealed()
        cdef Structure structure
        molecule._require_clean()
        structure = molecule._structure
        if query_may_match(query, structure):
            return True
        return False

    cpdef Py_ssize_t count(self, MoleculeContainer molecule,
                           bint automorphism_filter=False) except -1:
        """How many embeddings of this query the molecule admits.  Interruptible with Ctrl-C.

        The enumeration itself never touches Python, but it cannot be allowed to hold the GIL
        forever: eight unconstrained atoms against a 30-atom molecule is on the order of 1e11
        embeddings, and a `nogil` loop is not interruptible where Python bytecode is.  So the GIL comes
        back every 65536 embeddings purely to let pending signals raise -- one branch per embedding
        against the four AND tests per candidate.
        get_mapping and get_raw_mapping need no such thing: they yield to Python per embedding.
        """
        if molecule is None:
            raise TypeError('molecule must be a MoleculeContainer, not None')
        cdef Query query = self.sealed()
        cdef Structure structure
        cdef matcher_t mt
        cdef Py_ssize_t total = 0
        molecule._require_clean()
        structure = molecule._structure
        matcher_init(&mt, query, structure, automorphism_filter)
        try:
            with nogil:
                while matcher_next(&mt):
                    total += 1
                    if not (total & 0xFFFF):
                        with gil:
                            # cpython.exc declares this `except -1`, and that clause is what does the
                            # work: the C call leaves Python's error indicator set, and Cython turns
                            # that into a raised KeyboardInterrupt rather than a return value nobody
                            # inspects.  Without it an unbounded count() defers SIGINT until it ends.
                            PyErr_CheckSignals()
        finally:
            matcher_free(&mt)
        return total

    cpdef bint is_substructure(self, MoleculeContainer molecule,
                               bint automorphism_filter=False) except -1:
        """Does this query embed in the molecule at least once?  Stops at the first hit.

        automorphism_filter cannot change the answer -- an orbit is non-empty exactly when its
        members are -- and is accepted only so callers can pass the flag through uniformly.

        Unlike count this is NOT interruptible, and cannot be made so from here: it calls
        matcher_next exactly once, and a single search for a first embedding can itself be
        exponential.  A cancellation point would have to live inside matcher_next -- a node budget
        or a yield-control return -- which is a design change to a struct Tasks 11-13 are about to
        extend.  Parked deliberately until the kernel stops changing shape.
        """
        if molecule is None:
            raise TypeError('molecule must be a MoleculeContainer, not None')
        cdef Query query = self.sealed()
        cdef Structure structure
        cdef matcher_t mt
        cdef bint found
        molecule._require_clean()
        structure = molecule._structure
        matcher_init(&mt, query, structure, automorphism_filter)
        try:
            with nogil:
                found = matcher_next(&mt)
        finally:
            matcher_free(&mt)
        return found

    def __le__(self, other):
        """self matches inside other.  Returns NotImplemented for a non-molecule right side.

        Unfiltered, and it stays that way: an operator takes no keyword, and the filter cannot
        change a yes/no answer anyway (see is_substructure).
        """
        if not isinstance(other, MoleculeContainer):
            return NotImplemented
        return self.is_substructure(<MoleculeContainer> other)

    def __lt__(self, other):
        """self matches inside other and has strictly fewer atoms.

        The strict twin of `__le__`, with a `len()` guard that rejects before the kernel runs and
        makes `<` antisymmetric.  Only `<=` and `<` live on this class -- a query is the CONTAINED side
        of a containment test, and `q >= mol` would ask whether a molecule embeds in a pattern, which
        the kernel cannot answer in that direction.
        `mol >= q` and `mol > q` are the same two questions from the molecule's end.
        """
        if not isinstance(other, MoleculeContainer):
            return NotImplemented
        if len(self) >= len(<MoleculeContainer> other):
            return False
        return self.is_substructure(<MoleculeContainer> other)
