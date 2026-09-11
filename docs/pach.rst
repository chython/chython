Chython Binary Pack
===================

pach is chython's stored-structure format: a molecule or a reaction as a short byte string, small
enough to be a database column and cheap enough to decode in a loop over millions of rows.

It is not the arena.  ``to_bytes``/``from_bytes`` is the container's own memory written out verbatim,
lossless by construction and this release's own spelling; pach is the format a *different* chython
reads, and the one whose bytes may never change once a version byte ships.  Store with ``to_bytes``
what only this release will read back, and with ``pack`` what anything else will.

This page is the wire specification of every version byte there is, in enough detail to implement a
reader against.  The implementation-side statement of versions 0 and 2 is the header comment of
``chython/core/_pach.pxi``; versions 3 and 4 are in ``chython/core/_pach3.pxi``.  For the
surrounding I/O API see :doc:`io`.


Version space
-------------

======  ====================================================  ==============================
byte    record                                                status
======  ====================================================  ==============================
0       molecule, first generation                            read
1       reaction, first generation                            read
2       molecule, second generation                           read, and written on request
3       molecule, third generation, coordinates               read and written
4       molecule, third generation, coordinate-free           read and written
5       reaction, third generation                            read and written
======  ====================================================  ==============================

Molecules hold bytes 0, 2, 3 and 4; reactions hold bytes 1 and 5.  The two record kinds share one
version byte space and never overlap.

Versions 0 and 2 differ in one block, the bond orders, and are byte-identical everywhere else.
Versions 3 and 4 are one design and one reader; the version byte says whether the atom records carry
coordinates.  There is no version byte outside this table: 6 and up are unassigned, and a reader
handed one says so rather than guessing.


Framing
-------

A record starts at byte 0 with its version and carries no length field of its own.  Two consequences,
and both are the format's job rather than the caller's:

* **A record is self-delimiting.**  Its length is a function of its own header: arithmetic over the
  counts for versions 3 and 4, and a walk over the atom block's degree nibbles for 0 and 2, whose
  headers do not state the bond count.  ``pach_record_length`` answers it for the four *molecule*
  versions — 0, 2, 3 and 4 — so a caller can walk a stream of concatenated molecule records.  A
  reaction record is its four-byte header plus the molecule records it declares, each located by
  ``pach_record_length`` in turn; ``ReactionContainer.pack_len`` is the reaction-side reader that
  walks a body that way, and answers each molecule's atom count per side.
* **Trailing bytes are ignored.**  A decoder stops at the end of the record it was given and does not
  complain about what follows, which is what makes ``pach_load(stream)`` on a concatenation return the
  first record rather than an error.

**Compression is sniffed, not declared.**  ``pack()`` zlib-compresses by default.  A zlib stream's
first byte is its CMF, whose low nibble is the compression method and is always 8, so no version byte
above can be mistaken for one — ``0``, ``1``, ``2``, ``3``, ``4`` and ``5`` all have a low nibble that
is not 8, and so does the arena's own first byte ``0x33``.  Every door therefore accepts a compressed
record, an uncompressed one, or arena bytes without being told which.  ``compressed=True`` and
``compressed=False`` are for a caller who wants to be *told* rather than accommodated: either is a
``ValueError`` when the buffer disagrees.

**All multi-byte integers in versions 3, 4 and 5 are little-endian.**  Versions 0, 1 and 2 are
big-endian and bit-packed with no alignment anywhere.


The API
-------

.. testcode::

    from chython import smiles, MoleculeContainer, pach_load, pach_record_length

    benzene = smiles('c1ccccc1')
    raw = benzene.pack(compressed=False)

    print(raw[0], len(raw))                       # version 4: 12 header + 6*3 atom + 6*5 bond
    print(pach_record_length(raw, compressed=False))
    print(MoleculeContainer.unpack(raw) == benzene)

    back, problems = pach_load(raw, compressed=False)
    print(back == benzene, problems)

.. testoutput::

    4 60
    60
    True
    True []

``pach_load`` is the loop-safe door: it answers ``(molecule_or_None, problems)`` and never raises on
content.  A record it cannot use at all — including one whose graph the ring perception refuses — is
``None`` and a sentence saying so, never an exception.  ``pack``, ``MoleculeContainer.unpack``,
``pach_record_length`` and ``ReactionContainer.pack_len`` are the answer boundaries and raise
``ValueError`` — they return bytes, a molecule or a number and have no way to say "unknown".

``pach()`` is the short door over all of it, both directions on one name and both shapes through one
call; ``unpach()`` is its import half, ``unpack`` the same function under chython 2's name.  Byte 0
selects the decoder, so nothing about the buffer has to be declared — molecule pach, reaction pach and
``to_bytes`` output are told apart by it, and ``log=`` picks the error policy: absent it raises like
``unpack``, present it reports like ``pach_load``.

.. testcode::

    from chython import pach, unpach

    print(unpach(pach(benzene)) == benzene)
    print(unpach(pach(smiles('CCO>>CC=O'))))      # the reaction shape, same door

    log = []
    print(unpach(b'not a record', log=log), len(log))

.. testoutput::

    True
    C(C)O>>C(C)=O
    None 1

.. testcode::

    packed = benzene.pack()                       # zlib by default
    print(packed[0] & 0x0f)                       # a CMF's low nibble is always 8
    print(MoleculeContainer.unpack(packed) == benzene)

.. testoutput::

    8
    True

Walking a stream of records, each one located by the length of the one before it:

.. testcode::

    stream = smiles('CCO').pack(compressed=False) + smiles('CCN').pack(compressed=False)
    first = pach_record_length(stream, compressed=False)

    one, _ = pach_load(stream, compressed=False)
    two, _ = pach_load(stream[first:], compressed=False)
    print(one, two)

.. testoutput::

    C(C)O C(C)N

``version=`` chooses the layout: ``2``, ``3``, ``4``, or ``None`` for 3-or-4 by whether the molecule
has coordinates, which is the default.  ``2`` writes the legacy layout for a reader that predates this
release.  ``3`` asks for the coordinates the molecule *has*, so an undrawn molecule writes version 4
rather than a coordinate block of zeros stating a position nothing recorded; ``4`` on a drawn molecule
drops the drawing, which ``drop=['coordinates']`` says more plainly, and that waiver selects version 4
whatever ``version=`` asked for.  ``drop=`` waives the writer's refusals by name — ``map_number``,
``title``, ``meta``, ``sgroups``, ``cip``, ``wedges``, ``stereo_groups``, ``stereo``, ``coordinates``,
and ``conformers`` — or ``'*'`` for all of them, and an unrecognised name is refused rather than
ignored.

.. testcode::

    mol = smiles('CCO')
    print(mol.has_coordinates, mol.pack(compressed=False)[0])
    print(mol.pack(compressed=False, version=3)[0])        # no drawing, so version 4

    for i, n in enumerate(mol.atom_numbers):
        mol.set_xy(n, 1.5 * i, 0.0)

    raw = mol.pack(compressed=False)
    print(mol.has_coordinates, raw[0], len(raw))          # 12 + 3*9 + 2*5
    print(mol.pack(compressed=False, drop=['coordinates'])[0])
    print(mol.pack(compressed=False, version=3, drop=['coordinates'])[0])   # the waiver wins
    print(mol.pack(compressed=False, version=2)[0])

    mol.set_xyz(mol.atom_numbers[0], 0.0, 0.0, 1.0)       # a 3D conformer beside the drawing
    print(mol.pack(compressed=False, drop=['conformers'])[0])

.. testoutput::

    False 4
    4
    True 3 49
    4
    4
    2
    3


Versions 3 and 4
----------------

Header, 12 bytes
~~~~~~~~~~~~~~~~

.. code-block:: text

     0      version    u8    3 | 4
     1      flags      u8    bit 0 = map block present; bits 1-7 reserved, must be 0
     2-3    atoms      u16
     4-5    bonds      u16
     6-7    stereo     u16   configuration records
     8-9    sgroups    u16   enhanced-stereo entries
     10-11  reserved   u16   must be 0

The blocks follow in that order — atoms, bonds, stereo, enhanced stereo, map — each one a count from
the header times a constant stride, so a record's length is arithmetic:

.. code-block:: text

     length = 12
            + atoms   * (version == 3 ? 9 : 3)
            + bonds   * 5
            + stereo  * 9
            + sgroups * 3
            + (flags & 1 ? atoms * 2 : 0)

The header is the one place slack is deliberate: seven free flag bits and two reserved bytes.  A
future block needs a *count*, and the header is the only structure whose size cannot grow without
moving every offset.  A reader ignores a reserved bit it does not know, and reports that it did.

Atom block
~~~~~~~~~~

.. code-block:: text

     version 4, 3 bytes                  version 3, 9 bytes
       byte 0   el/R                       bytes 0-2  as version 4
       byte 1   iso 6 | rad 1 | pin 1      bytes 3-5  x, int24 LE two's complement
       byte 2   H 4 | chg 4                bytes 6-8  y, int24 LE two's complement

No field crosses a byte boundary.

``el/R`` (8 bits)
    Bit 7 clear: bits 6-0 are the atomic number, 1..118.  Bit 7 set: bits 6-0 are the R index,
    0..99 — ``0x80`` is a bare ``[R]`` and ``0x87`` is ``R7``.  37 codes are unused.

``iso`` (6 bits)
    ``MDL_ISOTOPE[z] - 32 + value`` for value 1..63, and 0 for unset: mass numbers from 31 below the
    element's MDL reference to 31 above.  The widest shift over chython's nuclides is 8.

``rad`` (1 bit)
    The radical flag.

``pin`` (1 bit)
    ``h_pinned`` — the hydrogen count was stated by the record that carried this atom, not derived.

``H`` (4 bits)
    The implicit hydrogen count, 0..14, and 15 for ``H_UNKNOWN``.  This is the container's own field
    with no translation, so a count nothing pins survives the round trip as unknown rather than as a
    zero.  The *explicit* hydrogen count needs no field: an explicit hydrogen is a real atom in the
    bond block.

``chg`` (4 bits)
    ``charge + 4``, admitting -4..+11 against the container's -4..+8.

``x``, ``y`` (24 bits each, version 3 only)
    The display coordinate at the container's own scale, x10000 fixed point — exact rather than
    rounded, reaching ±838.8607 Å.

**Atom order is load-bearing.**  Bonds, stereo slots, enhanced-stereo entries and map numbers all
address atoms by their *position* in this block, counting from 0, so a reader must keep the order
verbatim.  Stable ids are not stored: ``pack()`` renumbers, and an id is a container's private label
rather than chemistry.

Bond block, 5 bytes
~~~~~~~~~~~~~~~~~~~

.. code-block:: text

     bytes 0-1  a1     u16 LE, atom index
     bytes 2-3  a2     u16 LE, atom index
     byte 4     order 4 | wedge 4

``order``
    1 single, 2 double, 3 triple, 4 aromatic, 8 dative.  0 is invalid and 10 codes are reserved.
    ``order == 4`` is the whole of aromaticity — a record carries an unkekulised ring with nothing
    extra, and reads back as one.

``wedge``
    0 none, 1 up, 2 down, 3 either.  12 codes are reserved.

``a1 → a2`` is load-bearing: the wedge's narrow end is at ``a1``, the CTfile convention.  A bond
appears once — this is a bond list, not a bidirected table.

The block is identical in both versions, so a version 4 record has a wedge nibble and no drawing for
it to mean anything against.  The writer stores 0 there and logs the wedges it did not write; a reader
treats a nonzero wedge in a version 4 record as a problem and reads it as none.

Stereo block, 9 bytes
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: text

     bytes 0-1  slot0 u16     bytes 4-5  slot2 u16
     bytes 2-3  slot1 u16     bytes 6-7  slot3 u16
     byte 8     kind 3 | sign 1 | reserved 4

======  ==================  ==========  ========  ==========  ========  ================================
kind    name                slot 0      slot 1    slot 2      slot 3    implied by the graph
======  ==================  ==========  ========  ==========  ========  ================================
0       ``SU_TETRA``        centre      d0        d1          d2        the centre's fourth direction
1       ``SU_CIS_TRANS``    end A       a0        end B       b0        each end's other direction
2       ``SU_ALLENE``       end A       a0        end B       b0        each end's other direction
3       ``SU_ATROPISOMER``  pivot A     oA0       pivot B     oB0       each pivot's other neighbour
======  ==================  ==========  ========  ==========  ========  ================================

``kind`` uses four of eight codes.  ``sign`` is the parity, and every slot holds a real atom index:
there is no sentinel and no mask.

**A record states the owner of each direction list and all but one of its directions.**  The last is
implied by identity — "the owner's direction that is not already named" — not by an ordering
convention, which is what lets the direction ordering the container perceives change without
invalidating a stored record.  Four slots are complete because a direction list with *two* unnamed
directions is not stereogenic in the first place: a tetrahedral centre always has at least three named
directions, and each end of a cumulene or a biaryl axis has two directions of which at least one is
named.  So a record expresses a configuration exactly when that configuration means something.

Neither the allene chain's middle atoms nor an atropisomer's anchor is stored.  A reader finds the
unit from the slots: by the atom at slot 0 for a tetrahedral centre, by either pivot for an
atropisomer, and by the unique unit of matching kind holding both slot 1 and slot 3 among its
directions for the two bond kinds — slot 0 and slot 2 are then an adjacency cross-check.

.. testcode::

    from chython import STEREO_AND

    for text in ('N[C@@H](C)C(=O)O', 'C/C=C/C'):
        one = smiles(text)
        print(MoleculeContainer.unpack(one.pack()) == one)

    alanine = smiles('N[C@@H](C)C(=O)O')
    centre = alanine.atom_numbers[1]
    alanine.set_stereo_group(centre, STEREO_AND, 1)
    print(MoleculeContainer.unpack(alanine.pack()).stereo_group_of(centre))

    try:                                          # version 2 has no field for a stereo group
        alanine.pack(version=2)
    except ValueError as err:
        print('stereo_groups' in str(err))

.. testoutput::

    True
    True
    (3, 1)
    True

Enhanced stereo block, 3 bytes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: text

     bytes 0-1  atom  u16 LE, atom index
     byte 2     group 6 | kind 2

``kind`` is 0 unspecified, 1 abs, 2 or, 3 and; ``group`` is 1..63 for OR and AND and 0 otherwise.
The block is sparse and counted, so it costs nothing when unused.

It is its own block rather than a field in a configuration record because a group is a fact about an
*atom*: it can be set on an atom that owns no configuration, and for a bond configuration there is a
group at each end.

Map block, 2 bytes per atom
~~~~~~~~~~~~~~~~~~~~~~~~~~~

A ``u16`` per atom in atom order holding 0..9999, present when header flag bit 0 is set, which the
writer does when any atom has a nonzero map number.  The stored-structure corpus is unmapped and pays
nothing for the field; a **partially** mapped molecule is writable, and "mapped 1..N" is distinct from
"unmapped".

.. testcode::

    mapped = smiles('[CH3:1][OH:2]')
    raw = mapped.pack(compressed=False)
    print(raw[1] & 1, len(raw))                   # 12 + 2*3 + 1*5 + 2*2
    print(MoleculeContainer.unpack(raw) == mapped)

    partial = smiles('[CH3:1]O')                  # one atom mapped, one not
    back = MoleculeContainer.unpack(partial.pack())
    print([a.map_number for a in back.atoms()])

    try:                                          # version 2 has no field for a map number
        partial.pack(version=2)
    except ValueError as err:
        print('map_number' in str(err))

.. testoutput::

    1 27
    True
    [1, 0]
    True

What the writer refuses
~~~~~~~~~~~~~~~~~~~~~~~

A writer that drops a field quietly is discovered years later by whoever reads the record back, so
each of these raises ``ValueError`` naming what it refused, and the atom too where one atom is at
fault — the two count limits are a property of the whole molecule and name neither:

* more than 65535 atoms, or more than 65535 bonds
* a coordinate outside ±838.8607 Å, waivable with ``drop=['coordinates']``
* an isotope more than 31 mass units from its element's MDL reference
* an isotope on an ``R`` marker, which has no reference mass to shift from
* a formal charge outside -4..+11

**The format's other fields are at least as wide as the container's, so there is nothing to refuse.**
The hydrogen nibble is the arena's own 0..15 including ``H_UNKNOWN``; the R index field holds 0..127
against a maximum of 99; the map field holds 0..65535 against 0..9999; the enhanced-stereo group is the
arena's own six bits; and the configuration and enhanced-stereo counts are at most one per atom, so a
record of 65535 atoms or fewer cannot overfill either.  A *reader* must still accept each field's full
domain, because a record it did not write can state one.

What a version 3 or 4 record does not carry
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Each of these makes the writer refuse by name unless the matching ``drop=`` waiver is passed.  pach
has no text of any kind and no third coordinate; metadata ships in a form that does.

* ``title`` and ``meta``
* S-groups: superatoms, data labels, polymer brackets
* 3D conformers — a version 3 or 4 record carries 2D display coordinates only, and the refusal is on
  the models as a whole rather than on a first one: a molecule holding ten states them nowhere in this
  format.  ``drop=['conformers']`` waives it and keeps the drawing, so a drawn molecule still writes
  version 3; ``to_bytes()`` is the form that carries every model.
* CIP descriptors, which are derived and recomputed on demand
* stable ids, since ``pack()`` renumbers

Stereo is not on this list.  A version 3 or 4 record loses no configuration it is legal to hold, group
membership and atropisomers included.


Versions 0 and 2
----------------

The first two generations, written from chython 1.1 to 2.24 and read here.  Big-endian, bit-packed,
no alignment anywhere.  Version 0 is chython 1.1 through 1.44 and version 2 is 1.45 through 2.24; no
release wrote a version 1 molecule record.

.. code-block:: text

     header, 4 bytes    byte 0        format version, 0 | 2
                        12 bits       atom count             (data[1] << 4 | data[2] >> 4)
                        12 bits       cis/trans entry count  ((data[2] & 15) << 8 | data[3])

     atom block, 9 bytes per atom, in the writer's own atom order:
                        12 bits       atom number, 1..4095 -- the stable id, NOT an index
                        4 bits        neighbour count, 0..15
                        4 bits        stereo nibble (below)
                        5 bits        isotope shift, 0 = unset, else MDL reference - 16 + value
                        7 bits        atomic number, 1..118
                        16 + 16 bits  x, y as float16
                        3 bits        implicit hydrogen count, 7 = not stated
                        4 bits        formal charge + 4
                        1 bit         radical flag

     connection table   the atoms' neighbour lists concatenated in atom-block order, each entry a
                        12-bit atom NUMBER, two entries to three bytes.  Length is 3 * bond count,
                        and the bond count is half the sum of the neighbour counts: the table is
                        bidirected and names every bond twice.

     bond order block   one 3-bit value per bond, `order - 1`, in the order the connection table is
                        CONSUMED -- a bond's order appears where the table first names it, which is
                        at the lower-positioned of its two atoms.  THIS IS THE ONE BLOCK THE TWO
                        VERSIONS SPELL DIFFERENTLY: version 2 packs the values as a flat MSB-first
                        bitstream, ceil(3 * bonds / 8) bytes, and version 0 packs five values into
                        two bytes with one pad bit at the TOP, ceil(bonds / 5) * 2 bytes.

     cis/trans block    4 bytes per entry: two 12-bit atom numbers naming the two TERMINALS of a
                        cumulene chain, 7 pad bits, 1 sign bit.

The **stereo nibble** is two 2-bit fields, ``tetrahedron | allene``, and the writer chose between them
on the atom's neighbour count: an atom with exactly two neighbours got the allene field (``0b0010`` /
``0b0011``), anything else the tetrahedron field (``0b1000`` / ``0b1100``).  The second bit is the
sign.  The 2.24 reader collapses all four, so which field a record used is not recoverable from its
own answers, and this reader does not pretend otherwise: it reads the sign and lets the graph decide
which unit owns it.  Nothing is lost by that — an allene centre has two neighbours and a tetrahedral
centre has three or four.

Reading the header of a legacy record, since the counts are the same two a reader starts from:

.. testcode::

    legacy = smiles('C/C=C/C').pack(version=2, compressed=False)
    atoms = (legacy[1] << 4) | (legacy[2] >> 4)
    cis_trans = ((legacy[2] & 0x0f) << 8) | legacy[3]
    print(legacy[0], atoms, cis_trans)
    print(MoleculeContainer.unpack(legacy) == smiles('C/C=C/C'))

.. testoutput::

    2 4 1
    True

What versions 0 and 2 cannot carry
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Beyond the version 3 list above, all of which still applies: map numbers, wedges, enhanced stereo
groups, and every configuration that is not tetrahedral, allene or cis/trans — an atropisomer's parity
has nowhere to go.  Four hard ceilings have no ``drop=`` spelling because the record would be a lie
rather than a subset: an atom number above 4095, an atom with more than 15 neighbours, an implicit
hydrogen count above 6 (the field is 3 bits and 7 is the sentinel), and an isotope more than 15 mass
units from its element's MDL reference.

**The coordinate is degraded rather than refused**, since refusing it would refuse every drawn
molecule there is.  A float16 carries about three significant decimal digits, so ``1.2001953125`` comes
back as ``1.2002``; and the field is not a presence flag, so "this record has no drawing" and "every
atom sits at the origin" are the same four bytes.  A reader resolves that in favour of the first
reading when every coordinate in the record is zero.  Version 3's int24 is what removes both: it is
exact at the container's own scale, and version 4 says "no drawing" by being version 4.

A float16 spans a wider range than version 3's int24 field.  A coordinate whose absolute value exceeds
838.8607 Å is accepted into a version 0 or 2 record and refused when that record is re-encoded as
version 3.  ``drop=['coordinates']`` writes version 4 instead, storing the molecule without a drawing,
and it does so whether or not ``version=3`` was asked for beside it.


Reaction records
----------------

.. code-block:: text

     byte 0   version, 1 | 5
     byte 1   reactant count, u8
     byte 2   agent count,    u8
     byte 3   product count,  u8
     then     that many uncompressed molecule records, concatenated

The molecules are in ``molecules()`` order — reactants, then agents, then products — and that order is
load-bearing.  Each is a complete molecule record, located by ``pach_record_length`` and readable on
its own.

The two versions differ in the map numbers, and only there.  A version 1 body holds molecule records
of version 0 or 2, which have no map field, so the writer relabelled each molecule's atom *numbers*
onto its map numbers and the reader reverses that on load — which is why a molecule whose mapping was
partial could not be written at all.  A version 5 body holds version 3 or 4 records, where map numbers
are a field, so nothing is relabelled and nothing is refused.  ``pack(version=…)`` accepts ``None``,
``1`` or ``5``; any other value raises.

.. testcode::

    from chython import ReactionContainer

    rxn = smiles('[CH3:1][OH:2]>>[CH3:1][NH2:3]')
    raw = rxn.pack(compressed=False)

    print(raw[0], raw[1], raw[2], raw[3])
    print(ReactionContainer.unpack(raw) == rxn)
    print(rxn.pack(compressed=False, version=1)[0])

.. testoutput::

    5 1 0 1
    True
    1


Reading damaged records
-----------------------

Input is garbage by default, and the decode path never raises on content.  ``pach_load`` answers a
molecule together with a list of sentences saying what was wrong, because a loop over a million stored
records must not be stopped by one of them.  A record that is truncated, bit-flipped or holds a field
outside its domain yields the molecule the surviving bytes describe, and one whose graph cannot be built
at all yields ``None`` and the reason:

.. testcode::

    good = smiles('N[C@@H](C)C(=O)O').pack(compressed=False)
    mol, problems = pach_load(good[:-4], compressed=False)      # the last configuration record, cut
    print(mol is not None, len(problems) > 0)

.. testoutput::

    True True

**The header's counts define the record's layout.**  Every block starts where the declared counts put
it and stops at the earlier of its own declared count and the buffer's end.  An over-declared count
shifts the blocks behind it: a record declaring 5 bonds where 4 were written reads 5 bytes of what
follows the bond block as a fifth bond, and the stereo block starts 5 bytes late.  A truncated record
and an over-declaring one are byte-identical; there is nothing to distinguish them.  The header is the
layout, and a reader that tried to re-derive the offsets from what looked plausible would disagree with
``pach_record_length`` and with every other reader.

What a damaged record loses is bounded by the block the damage is in: a stereo record whose slots do
not resolve to a perceived configuration drops that configuration and says so, a bond naming an atom
that does not exist is dropped, and an enhanced-stereo entry naming one is dropped.  Two things are not
bounded that way and answer ``None`` instead: an atom block whose header count and bytes disagree, since
there are then no atoms to hang anything on, and a graph the ring perception refuses — its
relevant-cycle prototype limit and its deadline are reachable from a well-formed record, and a decode
reports them rather than raising.

``MoleculeContainer.unpack`` and ``ReactionContainer.unpack`` are the answer boundary: they raise
``ValueError`` naming what was wrong, with every problem the decoder found appended, including the ones
it recovered from.  A caller who cannot have a molecule is owed the whole story.
