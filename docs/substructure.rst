Substructure Search, SMARTS & Fingerprints
==========================================

Isomorphism and the ``<``/``<=`` operators, the chython SMARTS dialect in full, query building, and
molecular fingerprints.

For the SMIRKS form -- where a reactant side is a query and a product side is a patch -- see
:doc:`reactions`.


Substructure Check
------------------

.. testcode::

    from chython import smiles

    benzene = smiles('c1ccccc1')
    toluene = smiles('Cc1ccccc1')

    # Operator-based: < is proper substructure, <= is substructure or equal
    assert benzene < toluene        # True: benzene is a substructure of toluene
    assert benzene <= toluene       # True
    assert not (benzene < benzene)  # False: not a strict substructure of itself
    assert benzene <= benzene       # True: equal to itself

    # Method-based
    assert benzene.is_substructure(toluene)   # True
    assert not (benzene == toluene)            # False: different molecules

``as_query()`` is the one door all three go through, and it is the object to hold when one fragment is
tested against many molecules -- ``<`` and ``is_substructure`` compile a query per call and cache
nothing:

.. testcode::

    query = benzene.as_query()                # a QueryContainer
    assert query <= toluene

What it demands is element, isotope where one is stated, formal charge and the radical flag per atom,
and the order per bond. What it deliberately does not demand is everything that cannot survive
embedding: degree, hydrogen counts, hybridization, ring membership and ring size, the map number, and
stereo. That is why ``smiles('CO') < smiles('COC')`` is True -- methanol's oxygen carries a hydrogen and
the ether's does not. A stored aromatic bond is demanded as aromatic and does not match its Kekulé twin,
so two records written differently are compared after ``kekule()`` on both. A molecule carrying an R
marker is refused as a query, the marker matching nothing; it remains a legitimate target.

``may_contain()`` is the screen to run before the search, and it reads the other way round --
``target.may_contain(fragment)``:

.. testcode::

    assert toluene.may_contain(benzene)          # worth searching
    assert not benzene.may_contain(toluene)      # definitive: no mapping exists

**False is the answer; True is permission to search.** It compares which elements, charges, isotopes,
radicals and bond orders occur, and whether ring bonds occur -- presence, not how many -- so a fragment
larger than its target passes the screen and fails the search:

.. testcode::

    assert smiles('CC').may_contain(smiles('CCCC'))   # True, and the search says False


Enumerating Matches
-------------------

``get_mapping()`` yields all substructure mappings as dicts ``{query_atom: target_atom}``.
It is a method of ``QueryContainer``, which ``smarts()`` returns:

.. testcode::

    from chython import smarts, smiles

    query = smarts('CC')
    target = smiles('CCC')

    # First match
    mapping = next(query.get_mapping(target))
    assert len(mapping) == 2                       # {query_atom: target_atom}

    # All matches (automorphism_filter=True by default skips symmetric duplicates)
    filtered = list(query.get_mapping(target))

    # All symmetry-equivalent matches
    every = list(query.get_mapping(target, automorphism_filter=False))

    assert filtered and len(every) >= len(filtered)


SMARTS Queries
--------------

``smarts()`` returns a ``QueryContainer`` with pattern-matching semantics:

.. testcode::

    from chython import smarts, smiles

    # Carbonyl
    query = smarts('[C]=[O]')
    acid = smiles('CC(=O)O')
    assert query <= acid                       # True
    assert list(query.get_mapping(acid))       # at least one match

    # Aromatic nitrogen — note the explicit aromatic bond :
    assert smarts('[N;a]') <= smiles('c1ccncc1')   # True

    # Element containment shortcut
    assert 'N' in smiles('c1ccncc1')       # True
    assert 'Br' not in smiles('c1ccncc1')  # False


The SMARTS Language
-------------------

Chython SMARTS is not Daylight SMARTS. Recursive SMARTS ``$(...)``, valence ``v`` and total
connectivity ``X`` are not supported; hybridization ``z``, heteroatom count ``x``, the metal wildcard
``M`` and component grouping are.

There is exactly one parser, in the core, and everything reads through it -- ``smarts()``,
``read_smirks()``, and the SMARTS columns of the rule tables in ``chython/chemistry/tables/`` and
``chython/reactions/tables/``. A primitive that works in one works in all of them.

.. testcode::

    from chython import smarts, read_smarts, smiles

    query = smarts('[N;D1;z1;x0][C;z1]')       # a QueryContainer
    assert query <= smiles('CCN')

    log = []
    query = read_smarts('[C;a]', log=log)      # the same reader, with a log you own


Logical Operators
~~~~~~~~~~~~~~~~~

Both AND levels exist, with Daylight's precedence: ``&`` (high AND) binds tighter than ``,`` (OR),
which binds tighter than ``;`` (low AND). Bare juxtaposition **is** the high AND, so ``[CD1]``,
``[C&D1]`` and ``[C;D1]`` are one query.

``,`` is therefore OR *between whole high-AND groups*, not OR within one primitive type. That is
what makes the following a single pattern saying "a four-coordinate sp3 phosphorus **or** a
three-coordinate sp3 sulfur, cationic either way":

.. testcode::

    q = smarts('[P&D4&z1,S&D3&z1;+]')

    assert q <= smiles('C[P+](C)(C)C')
    assert q <= smiles('C[S+](C)C')

Comma OR within one primitive type is the ordinary case: ``[C;D2,D3]`` is degree 2 or 3,
``[C;r5,r6]`` is a carbon in a five- or six-membered ring, ``[C,N,O]`` is an element list.

.. testcode::

    assert smarts('[C;r5,r6]') <= smiles('C1CCCC1')
    assert not smarts('[C;r5,r6]') <= smiles('C1CC1')
    assert smarts('[C,N,O]') <= smiles('CO')


Atom Primitives
~~~~~~~~~~~~~~~

===========  ==========================================================  ========
token        meaning                                                     range
===========  ==========================================================  ========
``D``        degree -- heavy-atom neighbours, excluding implicit H        0-14
             **and dative bonds**
``z``        hybridization                                               1-6
``x``        heteroatom neighbours (not C, not H), **excluding dative**   0-14
``h``        implicit hydrogen count                                     0-14
``r``        ring size membership -- ``r5`` is "in a 5-membered ring"     3-14
``a``        aromatic flag; the same thing as ``z4``                      --
``M``        metal wildcard, and in a template the deletion mask         --
``+``/``-``  formal charge -- ``+2``, ``-1``, and so on                  --
``*``        any charge: withdraws the charge default, states nothing     --
``^``        this atom is a radical; ``!^`` is "not a radical"            --
``R``        ring count -- bare ``R`` is "in at least one ring",         0-99
             ``R2`` is "in exactly two"; ``R0`` and ``!R`` are acyclic
``#N``       atomic number -- ``#6`` is carbon                           1-118
``:N``       atom map number, for reaction templates                     1-9999
===========  ==========================================================  ========

An isotope is a prefix number, ``[13C]`` or ``[2H]``. Tetrahedral stereo is ``[@]``/``[@@]``, with
limited support in queries; the other two spellings, ``[@=]`` and ``[@~]``, are SMIRKS product-side tokens
and a query holding either is refused -- see :doc:`reactions` for what a patch can say about a
configuration. Radicals
may also be set from a CXSMARTS tail appended to the whole string, ``|^1:idx,...|``, by zero-based atom
index -- which suits a molecule written out once rather than a query whose atoms *are* the pattern.

.. testcode::

    assert smarts('[13C]') <= smiles('[13CH4]')
    assert not smarts('[13C]') <= smiles('C')

    assert smarts('[C;h3]') <= smiles('CC')            # a methyl
    assert not smarts('[C;h3]') <= smiles('C1CC1')

    assert smarts('[C;!R]') <= smiles('Cc1ccccc1')      # the exocyclic methyl
    assert not smarts('[C;!R]') <= smiles('C1CC1')


Hybridization
^^^^^^^^^^^^^

``z`` is 1-6, and **chython 2's** ``z3`` **is not chython 3's** ``z3``. chython 2 saturated: it
started at sp3, promoted on each double bond and capped at 3, so a sulfone sulfur, a nitro nitrogen
and an allene's central carbon all reported ``z3`` alongside genuine sp carbons. The core reports
what it actually found:

======  ==============================================================================
``z``   meaning
======  ==============================================================================
1       sp3
2       sp2
3       sp -- **and nothing else**
4       aromatic
5       two cumulated double bonds, no triple (allene, sulfone, sulfonyl, sulfonamide,
        pentavalent nitro)
6       any other combination: three or more doubles, or a double plus a triple
======  ==============================================================================

Any chython 2 template being ported has to have every ``z3`` in it re-read atom by atom. ``z1`` and
``z2`` mean the same thing in both.

.. testcode::

    allene = smiles('C=C=C')
    assert allene.hybridization_of(2) == 5           # the central carbon, not sp

    sulfone = smiles('CS(=O)(=O)C')
    assert sulfone.hybridization_of(2) == 5

    assert smarts('[C;z3]') <= smiles('CC#C')        # sp, and nothing else
    assert not smarts('[C;z3]') <= allene


Charge, wildcards and radicals
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``[A]`` is any *neutral* atom, ``[M]`` is any *neutral* metal (all 93), and ``[*]`` is any atom of
any charge -- ``*`` is the element wildcard and the charge wildcard in one token.

**An unstated charge means charge zero, not "any charge"**, for ``[A]``, ``[M]`` and ``[C]`` alike.
``[M]`` matches ``[Na]`` and not ``[Na+]``; to catch a metal salt, write ``[M;*]`` for any charge or
spell the charge out as ``[M;+]``. In a product template this is the difference between carrying a
salt through and neutralising it.

.. testcode::

    assert smarts('[M]') <= smiles('[Na]')
    assert not smarts('[M]') <= smiles('[Na+]')
    assert smarts('[M;*]') <= smiles('[Na+]')

    assert smarts('[A]') <= smiles('C')
    assert not smarts('[A]') <= smiles('[CH3-]')
    assert smarts('[*]') <= smiles('[CH3-]')

``*`` frees the charge and nothing else. It is the *withdrawal* of a default rather than a demand:
the seal neutralises a charge span nobody touched, so ``*`` touches the span and forbids nothing. It
therefore composes like any other constraint, and a charge written beside it still wins -- ``[C;*;+2]``
is a dication. ``!*`` is refused, because there is nothing to negate.

A radical is ``^`` inside the bracket, and ``*`` does **not** free it. ``[C;^]`` is a radical carbon,
``[C;!^]`` spells the unstated default out loud, and ``[C;^,!^]`` is "either". So ``[M;*]`` is a
metal of any charge that is *not* a radical, ``[M;*;^]`` is the radical one, and ``[M;*;^,!^]`` is a
metal in any charge and radical state at all.

.. testcode::

    radical = smiles('C[CH2] |^1:1|')

    assert smarts('[C;^]') <= radical
    assert not smarts('[C;^]') <= smiles('CC')
    assert smarts('[C] |^1:0|') <= radical              # the same field, from the tail

    sodium, cation = smiles('[Na]'), smiles('[Na+]')
    sodium_radical = smiles('[Na] |^1:0|')

    assert not smarts('[M;*]') <= sodium_radical        # any charge, but not a radical
    assert smarts('[M;*;^]') <= sodium_radical
    assert not smarts('[M;*;^]') <= cation
    assert all(smarts('[M;*;^,!^]') <= m for m in (sodium, cation, sodium_radical))

    # the tail composes with everything else: an aromatic C in a 5- or 6-ring, joined by a
    # non-ring single bond to a radical sp2/sp3 carbon carrying 0-1 hydrogens
    smarts('[C;r5,r6;a]-;!@[C;h0,h1;z1,z2;x0] |^1:1|')

No pattern atom matches an R.  ``[A]``, ``[C]``, ``[#6]`` and every disjunction of them fail on an R
target, and ``[R]`` in a pattern goes on meaning ring count.  The marker's own spelling is ``[#0]``, and
it is a SMIRKS **product-side** one: it states an attachment point to build, and a query holding one is
refused at the seal because there is nothing for it to match.  See :doc:`reactions`.

There is no collision with the dative bond, also written ``^``: a bond token is lexed *between*
atoms and the atom token only ever inside a bracket.


Bond Primitives
~~~~~~~~~~~~~~~

===================  =========================================  ====================
token                meaning                                    internal order
===================  =========================================  ====================
(implicit / absent)  single bond **only**                       1
``-``                single bond                                1
``=``                double bond                                2
``#``                triple bond                                3
``:``                aromatic bond                              4
``^``                dative (coordination) bond                 8
``~``                any bond -- all five of the above          --
===================  =========================================  ====================

**An implicit (absent) bond matches single bonds only, NOT aromatic.** To match an aromatic bond you
must write ``:``. This is the single most common source of SMARTS bugs:

.. testcode::

    benzene = smiles('c1ccccc1')

    assert not smarts('[C;a][C;a]') <= benzene    # implicit bond: single only
    assert smarts('[C;a]:[C;a]') <= benzene       # explicit aromatic bond

``~`` matches, and it is not order 8: it is the disjunction of all five orders, so ``C~C`` matches a
single, double, triple, aromatic **and** dative bond. ``-,=,#,:`` remains the way to say "any
covalent bond except a dative one", and ``!^`` is exactly that -- "not order 8" is a disjunction, not
a box of forbidden bits, so the reader expands it for you.

.. testcode::

    for target in ('CC', 'C=C', 'C#C', 'c1ccccc1'):
        assert smarts('C~C') <= smiles(target)

    dative = smiles('[Fe]~N(C)(C)C')
    assert smarts('C~C') <= smiles('CC')
    assert smarts('C!^C') <= smiles('CC')
    assert not smarts('C!^C') <= dative

Bond OR is ``-,=`` (single or double) or ``-,:`` (single or aromatic). Ring membership is a bond
modifier: ``-;@`` is a single bond in a ring, ``-;!@`` a single bond not in one, and it combines with
any order (``=;!@``). One caveat shared by ``~`` and ``!^``: a high AND written straight after either
expansion binds to its *last* alternative only, so ring membership on such a bond is ``!^;@`` and
never ``!^&@``.

.. testcode::

    assert smarts('[C]-;@[C]') <= smiles('C1CCCCC1')
    assert not smarts('[C]-;@[C]') <= smiles('CC')
    assert smarts('[C]!^;@[C]') <= smiles('C1CCCCC1')

``/`` and ``\`` are the one token here that is not a test of the bond it sits on: they state which side of
a double bond a substituent is on, so they carry the bond's order as well (single, as in SMILES) and stand
**in pairs**. A pair flanking a chain of double bonds asks for a target holding that geometry -- E or Z,
and an unconfigured alkene answers neither. One alone is refused, as is a pair putting one terminal's two
substituents on the same side, and a direction beside no double bond at all. Each combines with no other
bond token (``!/``, ``/,\``, ``/;-`` all raise) and is refused on a ring closure, whose two ends would
state it from opposite atoms. The same pair on a SMIRKS product side states a geometry to *build*; see
:doc:`reactions`.

.. testcode::

    from chython.core import IncorrectSmarts

    assert smarts('C/C=C/C') <= smiles('C/C=C/C')
    assert not smarts('C/C=C/C') <= smiles('C/C=C\\C')
    assert not smarts('C/C=C/C') <= smiles('CC=CC')       # never said is not "either"

    assert smarts('C(/C)=C/C') <= smiles('C/C=C\\C')      # read from the atom written FIRST

    for pattern in ('C/C=CC', 'C/C(\\C)=CC', 'C/CC', '[C]!/[C]', '[C]/,\\[C]', 'C/1=C/C1'):
        try:
            smarts(pattern)
        except IncorrectSmarts:
            pass
        else:
            raise AssertionError(pattern)


The dative bond is ``^`` here and ``~`` in SMILES
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The two dialects spell coordination differently. A molecule is *written* ``[Fe]~N(C)(C)C`` and
*matched* by ``[M]^[N;D3]``:

.. testcode::

    dative = smiles('[Fe]~N(C)(C)C')
    assert smarts('[M]^[N;D3]') <= dative

``D``, ``x`` and ``z`` all ignore a dative bond, while ``degree_of()`` and ``heteroatoms_of()`` do
not. Both counts are correct and neither is a stale copy of the other: the primitives count
*substituents*, because a coordination contact is not one, while the stored scalars are *structural*.

.. testcode::

    assert dative.degree_of(2) == 4               # the nitrogen, structurally 4-coordinate
    assert smarts('[N;D3]') <= dative             # but 3 substituents
    assert not smarts('[N;D4]') <= dative

    assert dative.heteroatoms_of(2) == 1          # the iron, structurally
    assert smarts('[N;x0]') <= dative             # but no heteroatom substituent

That is what makes ``[N;D4;z1;+]`` a complete quaternary-ammonium test: it cannot be satisfied by a
three-coordinate amine donating to a metal or to a borane.

.. testcode::

    assert smarts('[N;D4;z1;+]') <= smiles('C[N+](C)(C)C')
    assert not smarts('[N;D4;z1;+]') <= dative


Component Grouping
~~~~~~~~~~~~~~~~~~

``.`` says only "not bonded" -- the two fragments may share a molecule or not. To say which, group
them with parentheses at a **component position**: no atom to branch from and no branch open.

=============  =============================================================================
notation       meaning
=============  =============================================================================
``A.B``        unconstrained: one molecule or two, the query does not care
``(A.B)``      **one** molecule component -- this is how an intramolecular pattern is written
``(A).(B)``    **different** components
=============  =============================================================================

.. testcode::

    together, apart = smiles('NCCO'), smiles('N.O')

    intramolecular = smarts('(N.O)')
    assert intramolecular <= together
    assert not intramolecular <= apart

    intermolecular = smarts('(N).(O)')
    assert not intermolecular <= together
    assert intermolecular <= apart

    either = smarts('N.O')
    assert either <= together and either <= apart

A ``(`` means this only at a component position; once an atom precedes it, it opens a branch exactly
as it always did, group or no group. So ``(CC(C)C.N)`` is an isobutane and an amine demanded in one
molecule. Closing a component group also ends the component, so ``(C)(N)`` needs no ``.``.
``QueryContainer.component_groups()`` reports the groups as the seal reduced them -- per component,
not per atom.

.. testcode::

    smarts('(CC(C)C.N)').component_groups()
    # ((frozenset({1, 2, 3, 4}), 0), (frozenset({5}), 0)) -- two components, one group


Examples
~~~~~~~~

.. testcode::

    patterns = {
        'primary amine (NH2 on sp3 carbon)': '[N;D1;z1;x0:1][C;z1:2]',
        'aldehyde (O=CH-)':                  '[O;z2;x0:2]=[C;D2;x1;z2:1]',
        'aromatic carbon':                   '[C;a:1]',
        'pyrrole NH (aromatic bonds!)':      '[N;D2;z4;x0:1](:[C;a:2]):[C;a:3]',
        'sulfonyl chloride':                 '[S;x3;D4:1](=[O:2])(=[O:3])[Cl;D1:4]',
        'aryl boronic acid':                 '[B;D3;z1;x2:4](-[O;D1])(-[O;D1])-;!@[C;a:1]',
    }

    targets = ['CCN', 'CC=O', 'c1ccccc1', 'c1cc[nH]c1', 'CS(=O)(=O)Cl', 'OB(O)c1ccccc1']

    for pattern, target in zip(patterns.values(), targets):
        mol = smiles(target)
        mol.canonicalize()                    # the patterns are written against the aromatic form
        assert smarts(pattern) <= mol

Two habits those examples encode. The aromatic bonds in the pyrrole pattern are written out because
an absent bond would match nothing aromatic; and ``canonicalize()`` runs before matching because a
pattern written against ``[C;a]`` cannot fire on a Kekule record.


In the rule tables
~~~~~~~~~~~~~~~~~~

The SMARTS columns of ``chython/chemistry/tables/*.tsv`` and ``chython/reactions/tables/*.tsv`` are
read by the same parser, so everything above is available in a table row. Two conventions there are
worth stating, because a row is read by a human as often as by the lexer:

* ``[M]`` is neutral and non-radical in a table too, exactly as in a hand-written query, so a row
  meant to catch a metal in any state spells it ``[M;*;^,!^]``;
* a ``~`` in a **SMARTS** column is "any of five orders" and not the dative bond -- a row that means
  coordination writes ``^``. The ``examples`` columns are SMILES, where ``~`` *is* dative, so the
  two spellings sit side by side in one file and mean different things.


Query Building API
------------------

The recommended way to build queries is via SMARTS strings, which cover every query use
case (including ``&``, ``,``, ``;``, ``~``, component grouping and all primitives):

.. testcode::

    from chython import smarts, smiles

    # Acyclic thia/oxa carbonyl with 3 neighbors
    q = smarts('[C;D3;z2;!R;h0]=[O,S]')
    assert q < smiles('CC(=O)O')    # True (acid)
    assert q < smiles('CC(=S)C')    # True (thioketone)
    assert not (q < smiles('CC=O')) # False (aldehyde — only 2 neighbors)

For low-level programmatic construction, ``QueryContainer`` exposes a journal API:
``add_atom()`` allocates a stable id; ``atom_primitive(id, name, value)`` adds a constraint;
``add_bond(n, m)`` connects two atoms; ``bond_primitive(n, m, name, value)`` constrains a bond.
Primitive names are ``'element'`` (atomic number), ``'degree'``, ``'hybridization'``,
``'implicit_h'``, ``'ring_size'``, ``'bond_order'`` and others; the full list is in
``PRIM_NAMES`` in ``chython.core._query_boxes``.  For most use cases SMARTS is shorter and clearer.


Fingerprints
------------

Two families, five spellings each. ``morgan_*`` enumerates circular fragments (comparable to
ECFP); ``linear_*`` enumerates simple paths. The bit values are chython's own — no other
toolkit's hashes are reproduced.

============================ ======= ======== =========================
method                       folded  counted  returns
============================ ======= ======== =========================
``*_hash_set``               no      no       ``set[int]``
``*_hash_counts``            no      yes      ``dict[int, int]``
``*_bit_set``                yes     no       ``set[int]``
``*_fingerprint``            yes     no       ``ndarray(length)`` uint8
``*_count_vector``           yes     yes      ``ndarray(length)`` uint32
============================ ======= ======== =========================

Every method in both families needs numpy, which is an optional dependency — ``pip install
chython[ml]``. That includes the ``*_hash_set``, ``*_hash_counts`` and ``*_bit_set`` spellings, which
answer a plain ``set`` or ``dict``: they build the same uint32 invariant vector on the way there.

.. testcode::
    :skipif: __import__('importlib').util.find_spec('numpy') is None

    from chython import smiles

    mol = smiles('c1ccccc1O')

    # binary fingerprint, shape (1024,), dtype uint8
    fp = mol.morgan_fingerprint(min_radius=1, max_radius=4, length=1024,
                                number_active_bits=2)

    # the same positions as a set -- cheaper when the answer is a similarity
    bits = mol.morgan_bit_set(min_radius=1, max_radius=4, length=1024)

    # fragment counts rather than presence, shape (1024,), dtype uint32
    counts = mol.morgan_count_vector(min_radius=1, max_radius=4, length=1024)

    # unfolded: the fragment hashes themselves, and how many times each occurs
    hashes = mol.morgan_hash_set(min_radius=1, max_radius=4)
    occurrences = mol.morgan_hash_counts(min_radius=1, max_radius=4)

    # linear paths take the same arguments, reading the radii as path lengths in ATOMS:
    # length 1 is a lone atom, length 2 is a bond
    paths = mol.linear_fingerprint(min_radius=1, max_radius=4, length=1024)

``bit_set``, ``fingerprint().nonzero()`` and ``count_vector().nonzero()`` are the same answer in
three shapes. ``length`` must be a power of two, and ``number_active_bits * log2(length)`` may not
exceed 64 — a fragment hash is 64 bits wide, and past that bound every further bit is a constant.


Fingerprinting a different atom typing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Every method takes ``invariants=``, a ``(atom_count,)`` uint32 vector in the molecule's own atom
order. ``atom_invariants()`` returns the default scheme — element, isotope, charge, radical, total
hydrogen count, degree and ring membership — and any vector of that shape may be substituted.

.. testcode::
    :skipif: __import__('importlib').util.find_spec('numpy') is None

    labels = mol.atom_invariants()          # shape (atom_count,), dtype uint32
    fp = mol.morgan_fingerprint(invariants=labels)


Tanimoto Similarity
-------------------

.. testcode::
    :skipif: __import__('importlib').util.find_spec('numpy') is None

    import numpy as np

    mol1 = smiles('c1ccccc1O')
    mol2 = smiles('c1ccccc1N')

    fp1 = mol1.morgan_fingerprint()
    fp2 = mol2.morgan_fingerprint()

    # Via numpy
    tanimoto = np.dot(fp1, fp2) / (fp1.sum() + fp2.sum() - np.dot(fp1, fp2))

    # Via bit sets (faster)
    bits1 = mol1.morgan_bit_set()
    bits2 = mol2.morgan_bit_set()
    tanimoto = len(bits1 & bits2) / len(bits1 | bits2)
