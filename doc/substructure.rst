Substructure Search & Fingerprints
===================================

Isomorphism, SMARTS queries, query building, and molecular fingerprints.


Substructure Check
------------------

.. code-block:: python

    from chython import smiles

    benzene = smiles('c1ccccc1')
    toluene = smiles('Cc1ccccc1')

    # Operator-based
    benzene < toluene    # True: benzene is substructure of toluene
    benzene <= toluene   # True: substructure or equal
    benzene < benzene    # False: not strict substructure of itself
    benzene <= benzene   # True: equal

    # Method-based
    benzene.is_substructure(toluene)  # True
    benzene.is_equal(toluene)         # False


Enumerating Matches
-------------------

``get_mapping()`` yields all substructure mappings as dicts ``{query_atom: target_atom}``:

.. code-block:: python

    query = smiles('CC')
    target = smiles('CCC')

    # First match
    mapping = next(query.get_mapping(target))

    # All matches (automorphism_filter=True by default skips symmetric duplicates)
    for mapping in query.get_mapping(target):
        print(mapping)

    # All symmetry-equivalent matches
    for mapping in query.get_mapping(target, automorphism_filter=False):
        print(mapping)

    # Restrict search to specific atoms
    for mapping in query.get_mapping(target, searching_scope=[1, 2]):
        print(mapping)


SMARTS Queries
--------------

``smarts()`` returns a ``QueryContainer`` with pattern-matching semantics:

.. code-block:: python

    from chython import smarts, smiles

    # Carbonyl
    query = smarts('[C]=[O]')
    query <= smiles('CC(=O)O')  # True

    for mapping in query.get_mapping(smiles('CC(=O)O')):
        print(mapping)

    # Aromatic nitrogen
    smarts('[N;a]') <= smiles('c1ccncc1')  # True

    # Element containment shortcut
    'N' in smiles('c1ccncc1')   # True
    'Br' in smiles('c1ccncc1')  # False


SMARTS Language
---------------

Chython's SMARTS differs from RDKit/OpenBabel in several ways.

Atom Primitives
~~~~~~~~~~~~~~~

Standard:

- ``#N`` - atomic number (``#6`` for carbon)
- ``D`` - degree / neighbor count (``D3``)
- ``h`` - implicit hydrogen count (``h1``)
- ``r`` - ring size membership (``r5``, ``r6``)
- ``!R`` - acyclic (not in any ring)
- ``a`` - aromatic
- ``A`` - any element (wildcard)
- Charge: ``+``, ``-``, ``+2``, ``-3``
- Isotope: ``[14C]``
- Stereo: ``@``, ``@@``

Chython extensions:

- ``z`` - hybridization: ``z1`` = sp3, ``z2`` = sp2, ``z3`` = sp, ``z4`` = aromatic
- ``x`` - heteroatom neighbor count: ``x0`` = none, ``x2`` = two
- ``M`` - any metal (d-element)

.. code-block:: python

    smarts('[C;z2;x0]')                     # sp2 carbon, no heteroatom neighbors
    smarts('[O;D1;z1;x0][C;D3;x2;z2]=O')   # carboxylic acid
    smarts('[M]')                            # any metal

NOT Supported
~~~~~~~~~~~~~

- Recursive SMARTS ``$(...)``
- Valence ``v``
- Total connectivity ``X``
- Ring count ``R`` without size (use ``r`` with explicit sizes or ``!R``)
- Implicit AND ``&`` (use ``;`` instead)

Logical Operators
~~~~~~~~~~~~~~~~~

- ``;`` = AND between primitives: ``[C;D3;r6]``
- ``,`` = OR within same primitive type: ``[r5,r6]``, ``[C,N]``

OR cannot mix different primitive types: ``[D1,h1]`` raises an error.

.. code-block:: python

    smarts('[C;r5,r6;a]')   # aromatic C in 5- or 6-membered ring
    smarts('[C,N]')          # carbon or nitrogen
    smarts('[C;!R]')         # acyclic carbon

Bond Queries
~~~~~~~~~~~~

- ``-`` single, ``=`` double, ``#`` triple, ``:`` aromatic, ``~`` any
- OR: ``-,=`` (single or double)
- Negation: ``!-`` (not single), ``!:`` (not aromatic)
- Ring bonds: ``-;@`` (single and in ring), ``-;!@`` (single and not in ring)

Ring bond modifiers (``@``, ``!@``) require an explicit bond order prefix.

.. code-block:: python

    smarts('[C]-;!@[C]')   # non-ring single bond between two carbons
    smarts('[C]-;@[C]')    # single bond in a ring
    smarts('[C]-,=[C]')    # single or double bond

CXSMARTS Extensions
~~~~~~~~~~~~~~~~~~~

Radicals via CXSMARTS notation. Hybridization and heteroatom count use
chython's ``z`` and ``x`` primitives instead of CXSMARTS ``atomProp``:

.. code-block:: python

    # Radical on atom at position 1
    q = smarts('[C] |^1:0|')

    # Aromatic C in 5/6-ring, non-ring single bond to radical C
    # with 0-1 hydrogens — use z (hybridization) and x (heteroatoms) primitives
    q = smarts('[C;r5,r6;a]-;!@[C;h0,h1;z1,z2;x0] |^1:1|')


Query Building API
------------------

The recommended way to build queries is via SMARTS strings:

.. code-block:: python

    from chython import smarts, smiles

    # Acyclic thia/oxa carbonyl with 3 neighbors
    q = smarts('[C;D3;z2;!R;h0]=[O,S]')
    q < smiles('CC(=O)O')   # True (acid)
    q < smiles('CC(=S)C')   # True (thioketone)
    q < smiles('CC=O')       # False (aldehyde - only 2 neighbors)

For queries that cannot be expressed in SMARTS, build programmatically.
Add atoms first, then set query constraints on the atom objects:

.. code-block:: python

    from chython import QueryContainer, smiles
    from chython.containers.bonds import QueryBond
    from chython.periodictable import ListElement

    q = QueryContainer('')
    q.add_atom('C')
    a = q.atom(1)
    a.neighbors = 3
    a.hybridization = 2
    a.ring_sizes = 0
    a.implicit_hydrogens = 0

    q.add_atom(ListElement(['O', 'S']), n=3)
    q.add_bond(1, 3, 2)

    q < smiles('CC(=O)O')   # True (acid)
    q < smiles('CC(=S)C')   # True (thioketone)

    # Ring-ring linker using QueryBond(order, in_ring)
    q = QueryContainer('')
    q.add_atom('C')
    q.add_atom('C')
    q.atom(1).ring_sizes = 6
    q.atom(1).hybridization = 4
    q.atom(2).ring_sizes = 6
    q.atom(2).hybridization = 4
    q.add_bond(1, 2, QueryBond(1, False))  # single bond, NOT in ring

    q < smiles('c1ccc(cc1)-c1ccccc1')    # True (biphenyl)
    q < smiles('C1CC(=O)CC1')            # False

Query atom constraint values can be single integers or tuples for OR logic:

.. code-block:: python

    a = q.atom(1)
    a.neighbors = (2, 3)         # 2 OR 3 neighbors
    a.hybridization = 4          # aromatic only
    a.ring_sizes = (5, 6)        # in 5- or 6-membered ring
    a.implicit_hydrogens = 0     # no implicit H
    a.heteroatoms = 0            # no heteroatom neighbors
    # ring_sizes = 0 means "not in any ring"


Automorphism
------------

.. code-block:: python

    mol = smiles('c1ccccc1')

    mol.is_automorphic()  # True for benzene

    for mapping in mol.get_automorphism_mapping():
        print(mapping)


Maximum Common Substructure (MCS)
----------------------------------

.. code-block:: python

    mol1 = smiles('c1ccccc1O')   # phenol
    mol2 = smiles('c1ccccc1N')   # aniline

    for mapping in mol1.get_mcs_mapping(mol2, limit=10000):
        print(mapping)  # {mol1_atom: mol2_atom}
        break           # first = largest MCS


Morgan Fingerprints
-------------------

Similar to ECFP / RDKit Morgan fingerprints:

.. code-block:: python

    from chython import smiles
    import numpy as np

    mol = smiles('c1ccccc1O')

    # Binary fingerprint as numpy array (shape: (1024,))
    fp = mol.morgan_fingerprint(
        min_radius=1,
        max_radius=4,
        length=1024,
        number_active_bits=2,
    )

    # Bit indices only (more memory efficient)
    bits = mol.morgan_bit_set(min_radius=1, max_radius=4, length=1024)

    # Raw hashes (no folding) - useful for exact fragment matching
    hashes = mol.morgan_hash_set(min_radius=1, max_radius=4)


Linear Fingerprints
-------------------

Based on linear path fragments (similar to RDKit RDKFingerprint):

.. code-block:: python

    fp = mol.linear_fingerprint(
        min_radius=1,
        max_radius=4,
        length=1024,
        number_active_bits=2,
        number_bit_pairs=4,    # count-sensitive bits
    )

    bits = mol.linear_bit_set(min_radius=1, max_radius=4, length=1024)
    hashes = mol.linear_hash_set(min_radius=1, max_radius=4)


Tanimoto Similarity
-------------------

.. code-block:: python

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
