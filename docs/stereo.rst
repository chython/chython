Enhanced stereo
===============

A molecule's standard stereochemistry says whether a centre is configured, and what
configuration it carries.  **Enhanced stereo** says something additional: whether the
configuration is absolute (one known enantiomer), part of a racemate (both enantiomers
present), or relative (one enantiomer, absolute configuration unknown).  CTfile calls
these collections; chython calls them stereo groups.


The three kinds
---------------

============  ====================  =====================  =========================================
Kind          Atom-member spelling  Bond-member spelling   Meaning
============  ====================  =====================  =========================================
absolute      ``STEABS``            ``STEBABS``            the configuration is as drawn
racemate      ``STERACn``           ``STEBRACn``           AND group *n*: both configurations present
relative      ``STERELn``           ``STEBRELn``           OR group *n*: one configuration, unknown which
============  ====================  =====================  =========================================

For an axis member, ``STERACn`` is used when the chain midpoint is an atom and ``STEBRACn`` when it
is a bond.

**STERAC is a racemate, which is AND; STEREL is relative, which is OR.**  The file spells
the chemistry and the API spells the logic, and the two vocabularies do not share a letter.


Anchor and owners
-----------------

A group is stored at its unit's **anchor** slot, one byte per slot, exactly as a parity is; the
converse holds equally: a byte at a slot is read as the collection of whatever unit anchors there.
pach's bond-group entry keys the collection on the unit's **owners** instead, and a decode resolves
owners back to the anchor once the unit table is built.  The two are not always the same atom.

=========================  ============================  =======================================
Unit                       Owners (how it is named)      Anchor (where the byte lives)
=========================  ============================  =======================================
tetrahedral centre         the centre                    the centre
cis/trans double bond      the two chain terminals       the lower terminal
extended cis/trans         the two chain terminals,      the lower terminal
                           which are not bonded
allene, cumulated even     the two chain terminals       the chain's midpoint ATOM, which is
                                                         neither owner
atropisomer axis           the two pivots                one pivot, not always the lower
=========================  ============================  =======================================

The anchor is a **slot choice**, not a name: a seal derives the unit table again, so an edit can move a
unit's anchor -- breaking the ring double bond at a biaryl pivot frees that pivot and the axis moves to
it -- and the group byte is migrated with it.  Owners are constitution and survive that, which is why a
wire keys on them.

``stereo_groups()`` spells a one-owner unit as an ``int`` and a two-owner unit as an ascending pair,
and ``bond_stereo_groups()`` is the pair-spelled subset of the same dict -- a view, not a namespace.
``stereo_group_anchor_of(element)`` converts either spelling to the anchor, which is what a format
whose member list is atom-only needs.

A collection stated on an axis is a collection, so both readers see it:

.. testcode::

    from chython import smiles

    mol = smiles('COC(=O)/C=C\\C(=O)OC')
    with mol.edit() as e:
        e.set_bond_stereo_group(5, 6, 2, 1)
    print(mol.stereo_groups())
    print(mol.bond_stereo_groups())
    print(mol.stereo_group_anchor_of((5, 6)))

.. testoutput::

    {(2, 1): [(5, 6)]}
    {(2, 1): [(5, 6)]}
    5


Group ids are labels
--------------------

An id identifies a collection within one record and carries no meaning across records: it is
not a rank, not a count, and two records' ``OR 1`` are unrelated.  Three consequences:

- An id above 63 read from a file is renumbered to the lowest free id of its kind rather than
  dropped, and dropped only when all 63 ids of that kind are already taken.
- There is **one** id space.  A centre and an axis in ``AND 1`` are one collection with two members.
  V3000 spells the two halves ``STERAC1`` and ``STEBRAC1``, so a file that used both for two different
  mixtures has them renumbered apart on read, and an AND or OR collection that needs both spellings is
  split on write with its bond half renumbered (ABS carries no id, so it never splits).  Both are
  logged.
- ``canonical_stereo_groups()`` and ``canonical_bond_stereo_groups()`` assign the ids a caller should
  compare, ranking by member count and then by members' refinement classes ascending.
  ``canonical_bond_stereo_groups()`` is the pair-spelled view of ``canonical_stereo_groups()`` and
  shares its ids.


What identity does not include
-------------------------------

``canonical_bytes``, ``hash(mol)`` and ``mol == other`` fold **no** group membership, atom
or bond.  Two molecules differing only in which AND group their centres belong to have the
same identity bytes.  A grouped molecule is equal to its ungrouped self, and the two group
kinds are indistinguishable to identity.  This holds across the whole namespace.

``canonical_bond_stereo_group_ambiguities()`` reports two groups the colouring cannot
separate, and it errs toward reporting: a caller may be told two ids are unreliable where
the molecule's symmetry could not in fact exchange them.  A caller loses a guarantee and
never gets a wrong answer.


Format support
--------------

==============  =========  =========  ==========================================================
Format          Read       Write      Note
==============  =========  =========  ==========================================================
V3000           both       both       ``MDLV30/STE*`` with ``ATOMS=``, ``MDLV30/STEB*`` with
                                      ``BONDS=``; an axis is named by its chain midpoint, an atom
                                      for an allene and a bond otherwise, and an AND or OR
                                      collection clears the record's chiral flag
V2000           neither    neither    the chiral flag is one bit for the whole record, with no way
                                      to state which collection an atom is in
pach v3/v4      both       both       an atom-group block and a bond-group block, split by owner
                                      count, one namespace
pach v0/v2      neither    neither    the format is frozen; ``pack(version=2)`` raises
                                      ``ValueError`` naming ``stereo_groups``
MRV             both       both       ``mrvStereoGroup`` is an ``atomArray`` column, so an axis is
                                      named by its anchor atom
CXSMILES        both       both       ``|&N:|`` and ``|oN:|`` take atom indices, so an axis is named
                                      by one owner: the owner of lowest canonical position that owns
                                      no second unit
==============  =========  =========  ==========================================================


An atom-only member list, and the one axis it cannot name
---------------------------------------------------------

CXSMILES ``|&N:|`` and ``|oN:|``, and MRV's ``mrvStereoGroup``, take atom indices and nothing else, so
an axis -- named by two owners -- has one atom to speak for it.  ``set_stereo_group(atom, kind,
group)`` resolves an index to the unit anchored there or, failing that, to a unit that owns the atom,
so either an owner or the anchor round trips and an allene, whose anchor is neither owner, needs no
second spelling.

Which of the two the CXSMILES writer names is decided by canonical position among the owners and never
by the anchor, because the anchor is the slot the arena happened to choose: muconic acid with both C=C
axes in ``AND 1`` writes ``|&1:2,5|`` from every spelling of it, and a tail read back is written the
same way again.

An atom that owns two units names neither of them.  Where both owners of an axis do, the lower one is
written anyway -- a token is never suppressed on a prediction -- and the axis is reported
``smiles:stereo-group-axis-ambiguous``: the index carries the collection, and which of that atom's two
elements the collection is on is what it cannot carry.  ``C12=CC=CC=CC=C1C.CC1=CC=CC=CC=C12`` is the
case, its pivots being atropisomer owners and ring cis/trans terminals at once.


The names a caller sees in the log
-----------------------------------

Two destinations, and which one a record reaches follows from who produced it: a reader, a writer and
``depict`` take a ``log=`` list and record into that, while a seal has no such argument and records into
``molecule.log``.

On write, in the writer's ``log=`` list, what a format cannot spell and what is written instead:

=============================================  ==========================================================
Rule                                           When
=============================================  ==========================================================
``v2000:enhanced-stereo-not-written``          V2000 write; the molecule carries any AND or OR group
``v3000:collection-split``                     V3000 write (AND and OR only -- ABS never splits);
                                               one collection holds both atom-spelled and bond-spelled
                                               members, so its bond half is written under a free id;
                                               when all 63 ids of the kind are taken both halves keep
                                               the id and read back as one collection
``v3000:collection-axis-unspelled``            V3000 write; no chain joins an axis member's owners,
                                               so it is named by the lower one as an atom
``smiles:stereo-group-axis-ambiguous``         CXSMILES write; both owners of an axis own a second
                                               stereo element, so the one index written for it names
                                               no single element; written on the lower owner anyway
=============================================  ==========================================================

Two records come from an edit session, in ``molecule.log``, and neither drops a collection:

=============================================  ==========================================================
Rule                                           When
=============================================  ==========================================================
``edit:stereo-group-not-an-axis``              a group is stated on a pair no unit owns; the byte is
                                               kept as a label on the lower atom
``edit:stereo-group-anchor-taken``             the seal re-anchored a unit that carries a group onto an
                                               atom already stating one; the byte keeps the atom it is
                                               on as a label for the original collection, while the
                                               axis is still spelled under the collection that was
                                               already at the new anchor
=============================================  ==========================================================

A byte otherwise follows its unit through the seal: an edit that re-anchors an axis -- breaking the
ring double bond at a biaryl pivot frees that pivot, and the axis moves to it (ruling F45) -- moves the
collection with it.  On this clean path the unit is identified by its owners, so an atom owning two
units retains the right collection for each.

On V3000 read, where the collection block is parsed:

=============================================  ==========================================================
Rule                                           When
=============================================  ==========================================================
``v3000:collection-unknown``                   a collection name that is not an ``MDLV30/STE*`` tag
``v3000:collection-no-group``                  a numbered tag with no number, read as group 1
``v3000:collection-wrong-object``              an atom collection carrying ``BONDS=``, or the reverse
``v3000:collection-atom-ref``                  an atom index the CTAB does not have
``v3000:collection-bond-ref``                  a bond index the CTAB does not have
``v3000:collection-group-renumbered``          an id above 63, renumbered to the lowest free id
``v3000:collection-group-dropped``             an id above 63 with all 63 ids of its kind taken
``v3000:collection-namespaces-merged``         V3000 read; ``STERACn`` and ``STEBRACn`` share an id,
                                               so the bond one is read under a free id; when all 63
                                               ids of the kind are taken the two are read as one
                                               collection
``ctab:stereo-out-of-range``                   a collection references an atom position past the block
``ctab:stereo-group-dropped``                  the setter refused the entry -- an invalid kind or
                                               group id, or a bond group whose two atom ids are equal;
                                               a pair no unit owns is **kept** as a label at the lower
                                               atom, not dropped, and logged
                                               ``edit:stereo-group-not-an-axis`` instead
=============================================  ==========================================================

On CXSMILES read, records the group field produces:

=============================================  ==========================================================
Rule                                           When
=============================================  ==========================================================
``smiles:stereo-group-malformed``              the group field cannot be parsed
``smiles:stereo-group-empty``                  the group names no atom
``smiles:stereo-group-bad-index``              the group names an out-of-range atom
``smiles:stereo-group-duplicate``              the group assigns one atom twice
``smiles:stereo-group-renumbered``             the id is outside 1..63, renumbered to the lowest free id
``smiles:stereo-group-refused``                the setter refused the group
=============================================  ==========================================================

The MRV sample carries an axis through a format whose member list is atom-only:

.. testcode::

    from chython import smiles, mrv

    mol = smiles('COC(=O)/C=C\\C(=O)OC')
    with mol.edit() as e:
        e.set_bond_stereo_group(5, 6, 2, 1)
    back, = mrv(mrv(mol))
    print(back.bond_stereo_groups())

.. testoutput::

    {(2, 1): [(5, 6)]}
