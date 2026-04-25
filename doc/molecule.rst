Molecules
=========

Properties, atom/bond access, building, and stereochemistry.


Molecular Formula and Mass
--------------------------

.. code-block:: python

    from chython import smiles

    mol = smiles('CC(=O)Oc1ccccc1C(=O)O')  # aspirin

    mol.brutto                # {'C': 9, 'H': 8, 'O': 4}
    mol.brutto_formula        # 'C9H8O4'
    mol.brutto_formula_html   # 'C<sub>9</sub>H<sub>8</sub>O<sub>4</sub>'
    mol.molecular_mass        # 180.159... (average atomic masses)
    mol.molecular_charge      # 0 (total formal charge)


Drug-Likeness Descriptors
-------------------------

.. code-block:: python

    mol = smiles('CC(=O)Oc1ccccc1C(=O)O')

    mol.hydrogen_bond_donors_count     # N/O/S with H
    mol.hydrogen_bond_acceptors_count  # O/N/S with lone pairs
    mol.rotatable_bonds_count          # non-ring single bonds (excludes amide-like)
    mol.carbon_count
    mol.carbon_sp3_count
    mol.carbon_sp3_fraction            # sp3 carbons / total carbons


Ring Properties
---------------

.. code-block:: python

    mol = smiles('c1ccc2ccccc2c1')  # naphthalene

    mol.sssr                  # list of smallest rings as tuples of atom numbers
    mol.rings_count           # number of SSSR rings
    mol.atoms_rings           # dict: atom_number -> list of ring tuples
    mol.atoms_rings_sizes     # dict: atom_number -> set of ring sizes
    mol.aromatic_rings        # tuple of aromatic ring atom tuples


Other Properties
----------------

.. code-block:: python

    mol.atoms_count     # heavy atoms only
    mol.bonds_count
    mol.is_radical      # True if any atom is a radical


Iterating Atoms and Bonds
--------------------------

.. code-block:: python

    mol = smiles('CCO')

    # Iterate atom numbers
    for n in mol:
        print(n)

    # Iterate (atom_number, atom_object) pairs
    for n, atom in mol.atoms():
        print(n, atom.atomic_symbol, atom.atomic_number)

    # Iterate (n, m, bond) triples
    for n, m, bond in mol.bonds():
        print(n, m, int(bond))  # int(bond) = bond order: 1,2,3,4(aromatic)

    # Connected components
    components = mol.connected_components  # list of sets of atom numbers


Single Atom / Bond Access
--------------------------

.. code-block:: python

    atom = mol.atom(1)       # get atom by number
    bond = mol.bond(1, 2)    # get bond between atoms 1 and 2

    mol.has_atom(1)          # True/False
    mol.has_bond(1, 2)       # True/False


Atom Properties
---------------

.. code-block:: python

    atom = mol.atom(1)

    atom.atomic_symbol       # 'C', 'N', 'O', etc.
    atom.atomic_number       # 6, 7, 8, etc.
    atom.atomic_mass         # average atomic mass (float)
    atom.isotope             # isotope number or None
    atom.charge              # formal charge (int)
    atom.is_radical          # bool
    atom.implicit_hydrogens  # count of implicit H (int or None)
    atom.explicit_hydrogens  # count of explicit H neighbors
    atom.neighbors           # count of non-H neighbors
    atom.hybridization       # 1=sp3, 2=sp2, 3=sp, 4=aromatic
    atom.heteroatoms         # count of non-C, non-H neighbors
    atom.ring_sizes          # set of ring sizes containing this atom
    atom.x, atom.y           # 2D coordinates
    atom.xy                  # Vector(x, y) - supports tuple unpacking


Atom Neighbors / Environment
-----------------------------

.. code-block:: python

    # Full environment: (neighbor_num, bond, neighbor_atom)
    for n, bond, neighbor in mol.environment(atom_num):
        print(n, int(bond), neighbor.atomic_symbol)

    # Just neighbor numbers
    for n in mol.environment(atom_num, include_bond=False, include_atom=False):
        print(n)

    # (neighbor_num, bond) pairs
    for n, bond in mol.environment(atom_num, include_atom=False):
        print(n, int(bond))


Adjacency Matrix
-----------------

.. code-block:: python

    import numpy as np

    adj = mol.adjacency_matrix()       # 0/1 matrix
    adj = mol.adjacency_matrix(True)   # bond orders as values


Building Molecules
------------------

.. code-block:: python

    from chython import MoleculeContainer

    mol = MoleculeContainer()

    # Add atoms (returns atom number)
    n1 = mol.add_atom('C')          # from symbol
    n2 = mol.add_atom('C')
    n3 = mol.add_atom(8)            # from atomic number (oxygen)

    # Add bonds (bond order: 1=single, 2=double, 3=triple, 4=aromatic)
    mol.add_bond(n1, n2, 1)
    mol.add_bond(n2, n3, 2)

    print(str(mol))  # CC=O

    # Assign specific atom numbers
    mol = MoleculeContainer()
    mol.add_atom('C', n=10)     # atom number 10
    mol.add_atom('O', n=20)
    mol.add_bond(10, 20, 1)

    # Delete atom/bond
    mol.delete_bond(n2, n3)
    mol.delete_atom(n3)

    # Batch modifications (defer recalculation for performance)
    n4 = mol.add_atom('N', _skip_calculation=True)
    mol.add_bond(n2, n4, 1, _skip_calculation=True)
    mol.fix_structure()  # recalculate everything once


Merging and Splitting
---------------------

.. code-block:: python

    from chython import smiles

    # Split disconnected components
    anion, cation = smiles('[Cl-].[Na+]').split()
    print(anion, cation)

    # Merge molecules (union)
    salt = anion | cation
    salt = anion.union(cation, remap=True)  # fix atom number overlap

    # Extract substructure by atom numbers
    toluene = smiles('Cc1ccccc1')
    ring = toluene.substructure([2, 3, 4, 5, 6, 7])

    # Substructure with neighbors (1 bond deep)
    aug = toluene.augmented_substructure([2], deep=1)

    # Remap atom numbers (in-place; use copy() first to keep original)
    mol_copy = toluene.copy()
    mol_copy.remap({1: 10, 2: 20})

    # Copy
    mol_copy = mol.copy()


Stereochemistry
---------------

Inspecting
~~~~~~~~~~

.. code-block:: python

    mol = smiles('C/C=C/C')  # trans-2-butene

    mol.stereogenic_tetrahedrons   # dict: atom -> neighbors tuple
    mol.stereogenic_allenes        # dict: atom -> neighbors tuple
    mol.stereogenic_cis_trans      # dict: (n, m) bond -> substituents tuple
    mol.chiral_tetrahedrons        # set of atoms with assigned tetrahedral stereo
    mol.chiral_cis_trans           # set of bonds with assigned cis/trans stereo


Setting
~~~~~~~

.. code-block:: python

    mol = smiles('CC(O)F')

    # Add tetrahedral stereo
    # env = neighbor atom numbers defining chirality order
    # mark = True (counterclockwise / S) or False (clockwise / R)
    mol.add_atom_stereo(n=2, env=(1, 3, 4), mark=True)

    # Add cis/trans stereo to double bond
    # n, m = double bond atoms; n1, n2 = substituents
    # mark = True (cis) or False (trans)
    mol.add_cis_trans_stereo(n=2, m=3, n1=1, n2=4, mark=False)

    # Auto-detect cis/trans from 2D coordinates
    mol.calculate_cis_trans_from_2d()

    # Wedge/hash bond indicators
    mol.add_wedge(n=1, m=2, mark=1)   # 1 = wedge, -1 = hash

    # Clear all stereo
    mol.clean_stereo()

    # Recalculate stereo from current state
    mol.fix_stereo()


Hashing and Comparison
-----------------------

Molecules are hashable and comparable via canonical SMILES:

.. code-block:: python

    mol1 = smiles('CCO')
    mol2 = smiles('OCC')

    mol1 == mol2         # True (same canonical SMILES)
    hash(mol1) == hash(mol2)  # True

    # Use in sets and dicts
    unique = {smiles('CCO'), smiles('OCC'), smiles('c1ccccc1')}
    len(unique)  # 2

    # Cryptographic hash (SHA-512 based)
    sig = bytes(mol1)

**Warning**: Avoid modifying molecules (standardize, aromatize, add/remove atoms) after
placing them in sets or dicts. The hash will change and lookups will break.