Standardization
===============

Canonicalization, functional group normalization, aromaticity, tautomers, and deduplication.


Canonicalize
------------

``canonicalize()`` applies the full normalization pipeline: neutralize, standardize functional groups,
Kekule normalization, implicit hydrogen cleanup, aromatization, and charge standardization.

.. code-block:: python

    from chython import smiles

    mol = smiles('C(=O)(O)c1ccccc1')

    mol.canonicalize()
    print(str(mol))  # canonical SMILES

    # Options
    mol.canonicalize(
        fix_tautomers=True,    # canonical tautomer form (default)
        keep_kekule=False,     # return Kekule instead of aromatic
        logging=False,         # return list of changes made
        ignore=True,           # skip standardization bugs
    )

    # With logging: returns list of (atoms, rule_id, description) tuples
    log = mol.canonicalize(logging=True)
    for atoms, rule_id, description in log:
        print(f'{description} at atoms {atoms}')


Functional Group Standardization
---------------------------------

``standardize()`` normalizes functional groups (nitro, sulfoxide, etc.)
without changing aromaticity or tautomers. Over 80 rules applied:

.. code-block:: python

    mol = smiles('c1ccccc1N(=O)=O')  # nitro
    mol.standardize()                 # normalizes to [N+]([O-])=O form

    # With logging
    log = mol.standardize(logging=True)

    # Charge normalization (zwitterions)
    mol.standardize_charges()


Neutralize
----------

.. code-block:: python

    mol = smiles('[NH3+]CC(=O)[O-]')  # zwitterion
    mol.neutralize()                   # removes zwitterionic charges


Aromaticity
-----------

.. code-block:: python

    mol = smiles('c1ccccc1')

    # Convert aromatic (Thiele) to Kekule form
    mol.kekule()
    print(str(mol))  # C1=CC=CC=C1

    # Convert back to aromatic form
    mol.thiele()
    print(str(mol))  # c1ccccc1

    # Enumerate all Kekule structures
    for kekule_form in mol.enumerate_kekule():
        print(str(kekule_form))


Implicit / Explicit Hydrogens
------------------------------

.. code-block:: python

    mol = smiles('CCO')

    # Add explicit hydrogens
    added = mol.explicify_hydrogens()  # returns count of added H
    mol.clean2d()  # recalculate layout after adding atoms

    # Remove explicit hydrogens (make implicit)
    mol.implicify_hydrogens()

    # Fix implicit hydrogen counts
    mol.fix_structure()

``implicify_hydrogens`` works for aromatic rings only in Kekule form.
``explicify_hydrogens`` for aromatized forms requires ``kekule()`` then optionally ``thiele()`` afterward.


Tautomers
---------

.. code-block:: python

    mol = smiles('Oc1ccncc1')  # 4-pyridinol

    # Enumerate tautomers
    for tautomer in mol.enumerate_tautomers(limit=100):
        print(str(tautomer))

    # Include charge-shifted forms
    for tautomer in mol.enumerate_charged_tautomers(limit=100):
        print(str(tautomer))


Valence Checking
----------------

.. code-block:: python

    mol = smiles('C=N=Cc1ccccc1')

    # Check for valence problems (returns list of atom numbers with issues)
    errors = mol.check_valence()
    print('errors:', errors)

    # Aromatic rings must be kekulized first for accurate checking
    mol.canonicalize()
    errors = mol.check_valence()


Deduplication
-------------

Using Sets
~~~~~~~~~~

Molecules are hashable (based on canonical SMILES), so sets remove duplicates:

.. code-block:: python

    mols = [smiles('CCO'), smiles('OCC'), smiles('C(O)C')]

    for m in mols:
        m.canonicalize()

    unique = set(mols)  # 1 molecule (all three are ethanol)


Using Canonical SMILES
~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    seen = set()
    unique = []
    for mol in mols:
        mol.canonicalize()
        s = str(mol)
        if s not in seen:
            seen.add(s)
            unique.append(mol)


Graph Equality
~~~~~~~~~~~~~~

``is_equal()`` compares molecular graphs including all atom/bond properties:

.. code-block:: python

    mol1 = smiles('CCO')
    mol2 = smiles('OCC')

    mol1.is_equal(mol2)  # True
