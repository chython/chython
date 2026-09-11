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
# The per-element views: `Atom`, `Bond` and `Conformer`.
#
# A view is a borrowed handle -- a stable id plus the molecule's generation counter -- so it
# raises rather than read an arena that someone else's edit has moved out from under it.  The
# query side has no counterpart to this file by decision; see RULES.md section 1.3.


cdef class Atom:
    cdef MoleculeContainer _molecule
    cdef uint32_t _n
    cdef uint32_t _gen

    cdef inline atom_t *_ptr(self) except NULL:
        if self._gen != self._molecule._gen:
            raise RuntimeError('stale Atom view: the molecule was mutated')
        return self._molecule._atom(self._n)

    cdef inline int _synced(self) except -1:
        # A setter goes through the molecule, so outside an edit scope it applies immediately and
        # bumps the generation. Re-arm this view against that: an atom handle must survive the
        # writes made through it, or `for a in mol.atoms(): a.charge = 0` would die on its
        # second read. Someone else's edit still invalidates the view, which is the point.
        self._gen = self._molecule._gen
        return 0

    @property
    def n(self):
        """This atom's number -- its stable id. Read-only: an id is issued, never assigned."""
        return self._n

    @property
    def element(self):
        return self._ptr().element

    @property
    def atomic_symbol(self):
        return symbol_of(self._ptr())

    @property
    def atomic_radius(self):
        """The calculated atomic radius in angstroms, and 0.0 for the R marker, which carries none.

        Element data on the view because its reader is a renderer already holding the atom.  See
        `el_atomic_radius` for which radius this is and where the published set stops.
        """
        return el_atomic_radius(self._ptr().element)

    @property
    def is_r(self):
        """True for the fragment marker, element 0.  A marker matches nothing and carries no mass."""
        return self._ptr().element == 0

    @property
    def r_index(self):
        """The R index, or 0 for a plain R and for every element."""
        return at_r_index(self._ptr())

    @property
    def charge(self):
        return self._ptr().charge

    @charge.setter
    def charge(self, int value):
        self._molecule.set_charge(self._n, value)
        self._synced()

    @property
    def isotope(self):
        return self._ptr().isotope

    @isotope.setter
    def isotope(self, int value):
        self._molecule.set_isotope(self._n, value)
        self._synced()

    @property
    def map_number(self):
        return self._ptr().map_number

    @map_number.setter
    def map_number(self, int value):
        self._molecule.set_map_number(self._n, value)
        self._synced()

    @property
    def degree(self):
        return self._ptr().degree

    @property
    def neighbors(self):
        # A COUNT, not a way to walk the graph.  The walk is mol.neighbors_of(id), which is where a
        # graph question belongs.
        return self._ptr().degree

    @property
    def radical(self):
        return at_radical(self._ptr())

    @radical.setter
    def radical(self, bint value):
        self._molecule.set_radical(self._n, value)
        self._synced()

    @property
    def is_radical(self):
        return at_radical(self._ptr())

    @is_radical.setter
    def is_radical(self, bint value):
        self._molecule.set_radical(self._n, value)
        self._synced()

    @property
    def stereo(self):
        self._ptr()                       # the generation guard; the read itself is by slot
        return self._molecule.stereo_of(self._n)

    @stereo.setter
    def stereo(self, bint value):
        self._molecule.set_stereo(self._n, value)
        self._synced()

    @property
    def parity(self):
        """0 no parity configured, 1 even, 2 odd -- `parity_of` for this atom, delegated to it."""
        self._ptr()
        return self._molecule.parity_of(self._n)

    @property
    def stereo_group(self):
        """`(kind, group)`, or `(STEREO_UNSPECIFIED, 0)` when this atom is in no collection.

        Delegated rather than re-read: the container's body is keyed by arena index, which is a lookup
        this view would have to do anyway, so a second copy would buy nothing and drift.  What section
        4.0 forbids is the CALLER doing this, not the view.
        """
        self._ptr()
        return self._molecule.stereo_group_of(self._n)

    @property
    def cip(self):
        """The stored CIP descriptor -- 'R' 'S' 'r' 's' 'M' 'P' 'm' 'p' -- or None.  Storage only."""
        return ATOM_CIP_CODES[at_cip(self._ptr())]

    @property
    def in_ring(self):
        return at_in_ring(self._ptr())

    @property
    def implicit_h(self):
        """This atom's implicit hydrogen count, or None when the record did not state one.

        Assign `H_UNKNOWN` to write that state; assigning None raises, because the setter's job is to
        write a value and None is not one.
        """
        cdef atom_t *a = self._ptr()
        if at_implicit_h_unknown(a):
            return None
        return at_implicit_h(a)

    @implicit_h.setter
    def implicit_h(self, int value):
        self._molecule.set_hydrogens(self._n, value)
        self._synced()

    @property
    def explicit_h(self):
        return at_explicit_h(self._ptr())

    @property
    def heteroatoms(self):
        return self._ptr().heteroatoms

    @property
    def hybridization(self):
        return at_hybridization(self._ptr())

    @property
    def total_h(self):
        """Implicit plus explicit, or None when the implicit count is unknown -- a sum with an
        unknown term is unknown.  `explicit_h` is always a number."""
        cdef atom_t *a = self._ptr()
        if at_implicit_h_unknown(a):
            return None
        return at_implicit_h(a) + at_explicit_h(a)

    @property
    def ring_count(self):
        return at_ring_count(self._ptr())

    @property
    def ring_sizes(self):
        cdef uint32_t w = self._ptr().ring_sizes
        cdef uint32_t size
        cdef list sizes = []
        for size in range(3, 25):
            if w >> size & 1:
                sizes.append(size)
        return frozenset(sizes)

    @property
    def macrocycle(self):
        """True when this atom lies on a ring larger than 24, whose size ring_sizes cannot hold."""
        return (self._ptr().ring_sizes & 7) != 0

    @property
    def x(self):
        # Guard: check generation then require_clean; resolve the arena index once and pass it
        # to xy_read_x so _index_of[n] is looked up exactly once.  Matches the bare-call
        # idiom of Bond.wedge except that the index is kept for the xy read, which is why _ptr()
        # (which would resolve it again inside _atom()) is not used here.
        if self._gen != self._molecule._gen:
            raise RuntimeError('stale Atom view: the molecule was mutated')
        cdef MoleculeContainer mol = self._molecule
        mol._require_clean()
        if not structure_has(mol._structure, SEG_XY):
            return None
        cdef uint32_t idx = <uint32_t> mol._index_of[self._n]
        return xy_read_x(structure_xy(mol._structure) + idx)

    @x.setter
    def x(self, double value):
        cdef object xy = self._molecule.xy_of(self._n)
        self._molecule.set_xy(self._n, value, 0.0 if xy is None else xy[1])
        self._synced()

    @property
    def y(self):
        # Guard: same single-resolution pattern as Atom.x above.
        if self._gen != self._molecule._gen:
            raise RuntimeError('stale Atom view: the molecule was mutated')
        cdef MoleculeContainer mol = self._molecule
        mol._require_clean()
        if not structure_has(mol._structure, SEG_XY):
            return None
        cdef uint32_t idx = <uint32_t> mol._index_of[self._n]
        return xy_read_y(structure_xy(mol._structure) + idx)

    @y.setter
    def y(self, double value):
        cdef object xy = self._molecule.xy_of(self._n)
        self._molecule.set_xy(self._n, 0.0 if xy is None else xy[0], value)
        self._synced()

    @property
    def xy(self):
        self._ptr()
        return self._molecule.xy_of(self._n)

    @xy.setter
    def xy(self, value):
        cdef double x, y
        x, y = value
        self._molecule.set_xy(self._n, x, y)
        self._synced()

    def __int__(self):
        """The atomic number, so `int(atom)` is `atom.element`.

        The twin of `int(bond)` being the bond order: in both cases the atom's or bond's ONE
        defining number.
        """
        return self._ptr().element

    def __eq__(self, other):
        """Compare against an element SYMBOL or an atomic NUMBER: `atom == 'C'`, `atom == 6`.

        The shape `if (atom := mol.atom(n)) == 'C'` is what makes `'C' in mol` read as English.
        Isotope, charge and radical state are NOT compared, deliberately: the question `atom == 'C'`
        asks is "is this a carbon", and an isotope-aware answer would make `atom == 'C'` False on a
        13-C, which no caller means.  Compare `atom.isotope` for that.

        `atom == 'R'` is True for any R and `atom == 'R7'` only for that index, which is the same two
        grains the isotope has: the index annotates the marker, it is not a different kind of atom.
        Both hold together with `atom == atom.atomic_symbol`, an unindexed R spelling `R`.

        Two `Atom` views compare by identity of the atom they name -- same molecule, same number --
        and NOT by their properties, because a view is a handle and two handles on one atom are the
        same atom.  Anything else returns `NotImplemented` so Python falls back to identity.
        """
        cdef uint32_t number
        if isinstance(other, str):
            if other == 'R' or (other.startswith('R') and other[1:].isdigit() and other != 'R0'):
                if self._ptr().element:
                    return False
                if other == 'R':
                    return True
                return int(other[1:]) <= R_INDEX_MAX and self.r_index == int(other[1:])
            number = SYMBOL_TO_NUMBER.get(other, NOT_AN_ELEMENT)
            return number != NOT_AN_ELEMENT and self._ptr().element == number
        if isinstance(other, int) and not isinstance(other, bool):
            return <int> self._ptr().element == <int> other
        if isinstance(other, Atom):
            return (self._molecule is (<Atom> other)._molecule
                    and self._n == (<Atom> other)._n)
        return NotImplemented

    def __hash__(self):
        """Hash the atom's IDENTITY -- its molecule and its number -- to match `__eq__`'s view of
        two `Atom` handles.  NOT hashed to agree with `atom == 'C'`: an atom is not a symbol, and
        making one hash equal to a string's would put atoms and strings in one bucket of every
        dict.  So `atom == 'C'` is True while `hash(atom) != hash('C')`, which breaks the hash
        invariant across TYPES on purpose, the alternative being worse."""
        return hash((id(self._molecule), self._n))

    def __repr__(self):
        return 'Atom(%s, n=%d)' % (symbol_of(self._ptr()), int(self._n))


cdef class Bond:
    cdef MoleculeContainer _molecule
    cdef uint32_t _n
    cdef uint32_t _m
    cdef uint32_t _gen

    cdef inline halfedge_t *_ptr(self) except NULL:
        if self._gen != self._molecule._gen:
            raise RuntimeError('stale Bond view: the molecule was mutated')
        cdef MoleculeContainer mol = self._molecule
        mol._require_clean()
        cdef halfedge_t *e = csr_find(mol._structure, <uint32_t> mol._index_of[self._n],
                                     <uint32_t> mol._index_of[self._m])
        if e is NULL:
            raise RuntimeError('stale Bond view: the bond no longer exists')
        return e

    @property
    def n(self):
        """Stable id of one endpoint. Ids, not atoms: fetch those with mol.atom(bond.n)."""
        return self._n

    @property
    def m(self):
        """Stable id of the other endpoint."""
        return self._m

    @property
    def order(self):
        return self._ptr().order

    @order.setter
    def order(self, int value):
        cdef MoleculeContainer mol = self._molecule
        mol.set_order(self._n, self._m, value)
        self._gen = mol._gen

    @property
    def in_ring(self):
        return (self._ptr().flags & HE_IN_RING) != 0

    @property
    def wedge(self):
        """`(narrow_n, code)` for the wedge drawn on this bond, or None when it carries none.

        Not a bare code: a wedge is DIRECTIONAL and a `Bond` is not (see `__eq__`, which compares
        endpoints unordered), so the answer to "which way does it point" is an atom.  Delegated to
        `wedge_between`, which is where the two half-edge reads live.
        """
        self._ptr()   # the generation guard, and the refusal on a bond that no longer exists
        return self._molecule.wedge_between(self._n, self._m)

    @property
    def cip(self):
        """The stored CIP descriptor of this bond -- 'E' 'Z' 'M' 'P' -- or None.  Storage only."""
        return BOND_CIP_CODES[he_cip(self._ptr())]

    def __int__(self):
        return self._ptr().order

    def __eq__(self, other):
        """Compare against a bond ORDER: `bond == 1`, `bond == 4`.

        An aromaticity test is written `if bond == 4`, and that spelling is the reason `__int__`
        exists next door.  `bond == 4` is therefore the aromatic
        test, and it is exact rather than a perception question: order 4 IS how the arena stores an
        aromatic bond.

        Two `Bond` views compare by the bond they name, unordered in their endpoints, so a view
        taken as (n, m) equals one taken as (m, n): a bond has no direction, and `bond(1, 2) ==
        bond(2, 1)` returning False would be a trap.  `set_wedge`'s narrow/wide pair is the one
        place a pair of atoms IS ordered, and it is not spelled with a Bond.
        """
        if isinstance(other, int) and not isinstance(other, bool):
            return <int> self._ptr().order == <int> other
        if isinstance(other, Bond):
            return (self._molecule is (<Bond> other)._molecule
                    and {self._n, self._m} == {(<Bond> other)._n, (<Bond> other)._m})
        return NotImplemented

    def __hash__(self):
        """Hash the bond's IDENTITY, endpoints unordered, to match `__eq__` between two views.  Not
        hashed to agree with `bond == 4`; see `Atom.__hash__` for why that asymmetry is chosen."""
        return hash((id(self._molecule), min(self._n, self._m), max(self._n, self._m)))

    def __repr__(self):
        return 'Bond(%d, %d-%d)' % (self._ptr().order, int(self._n), int(self._m))


cdef class Conformer:
    """One model of a molecule's geometry: its coordinates and the number the file gave it.

    A borrowed handle like `Atom` and `Bond`, holding a model INDEX rather than a stable id -- and the
    difference is real, because an index is a list position that `drop_conformer` moves.  The
    generation check covers it: any edit at all invalidates the view, so a compaction cannot be read
    through one.
    """
    cdef MoleculeContainer _molecule
    cdef uint32_t _index
    cdef uint32_t _gen

    cdef inline xyz_t *_ptr(self) except NULL:
        if self._gen != self._molecule._gen:
            raise RuntimeError('stale Conformer view: the molecule was mutated')
        cdef MoleculeContainer mol = self._molecule
        mol._require_clean()
        if self._index >= structure_conformer_count(mol._structure):
            raise RuntimeError('stale Conformer view: the model no longer exists')
        return structure_conformer_xyz(mol._structure, self._index)

    @property
    def index(self):
        """This model's position.  Read-only, and not a name: a drop renumbers the models above it."""
        return self._index

    @property
    def ext_index(self):
        """The number the file gave this model -- a PDB `MODEL`, an XYZ frame ordinal -- or None when
        nothing stated one.  Stored verbatim and never interpreted, so zero is a number and not a gap.
        """
        self._ptr()      # the generation check and the existence check, in that order
        cdef conformer_t *rec = structure_conformer_records(self._molecule._structure) + self._index
        return None if rec.ext_index == <uint32_t> CONF_NO_INDEX else int(rec.ext_index)

    def xyz_of(self, uint32_t n):
        """Atom `n`'s `(x, y, z)` in this model.

        Never None, where `MoleculeContainer.xyz_of` is None for a molecule with no conformer at all:
        this view exists only because a model does.  An atom the caller never placed reads the origin
        -- one dense column per model, and design D7 records that as the slice's known defect.
        """
        cdef xyz_t *base = self._ptr()
        cdef MoleculeContainer mol = self._molecule
        if n not in mol._index_of:
            raise KeyError(n)
        cdef xyz_t *p = base + <uint32_t> mol._index_of[n]
        return (xyz_read_x(p), xyz_read_y(p), xyz_read_z(p))

    @property
    def coordinates(self):
        """The whole model as a list of `(x, y, z)`, in `atom_numbers` order."""
        cdef xyz_t *base = self._ptr()
        cdef uint32_t i
        cdef list out = []
        for i in range(self._molecule._structure.header.atom_count):
            out.append((xyz_read_x(base + i), xyz_read_y(base + i), xyz_read_z(base + i)))
        return out

    def __eq__(self, other):
        """By IDENTITY -- molecule and index -- as two `Atom` handles are compared.  Never by
        coordinates: two models of one molecule that happen to agree are still two models."""
        if isinstance(other, Conformer):
            return (self._molecule is (<Conformer> other)._molecule
                    and self._index == (<Conformer> other)._index)
        return NotImplemented

    def __hash__(self):
        return hash((id(self._molecule), self._index))

    def __repr__(self):
        return 'Conformer(index=%d, of %d atom(s))' % (
            int(self._index), int(self._molecule._structure.header.atom_count))
