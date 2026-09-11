Chython
=======

Library for processing molecules and reactions in Python.

.. code-block:: bash

   pip install chython          # no numpy; about 7.5 MB of runtime files
   pip install 'chython[ml]'    # adds numpy, and with it fingerprints and the array surface

numpy is optional. Everything below works without it except the entries marked *needs* ``chython[ml]``,
and those raise an :class:`ImportError` naming the extra when called rather than at import.

**Key capabilities:**

- Read and write SMILES, SMARTS and SMIRKS, InChI, MDL MOL/SDF/RXN/RDF in both V2000 and V3000,
  CML, Marvin MRV, and Tripos MOL2
- Read PDBx/mmCIF, legacy PDB and XYZ. Each returns a record of atoms and coordinates rather than a
  molecule; PDB and mmCIF state connectivity and no orders, and :func:`chython.saturate` is the
  separate, explicitly invoked pass that perceives the orders for a ligand. XYZ states no bond either,
  and :func:`chython.perceive_bonds` is the explicitly invoked pass that reads the connectivity out of
  a stored geometry
- Read IUPAC names through OPSIN
- Substructure search with chython SMARTS, including component grouping for intramolecular patterns
- Canonical form, standardization, kekulization and aromatization, resonance repair, salt stripping
- Morgan and linear fingerprints, as hash sets, folded bit vectors or count vectors
  -- *needs* ``chython[ml]``, the hash-set and bit-set spellings included: every one of them builds
  the same numpy invariant vector first, whatever it finally returns
- Descriptors: TPSA, Crippen logP and MR, hydrogen-bond donors and acceptors, rotatable bonds,
  ring counts, Bertz CT, Randić and Zagreb indices
- The 166 MACCS structural keys, one-based and read off the published key descriptions, and QED with
  its three published weight sets. Both state what they are transcriptions of and neither claims parity
  with another implementation's bits or score
- Graph descriptors over the topological distance matrix -- eccentricities, Wiener index, radius,
  diameter, Balaban J -- and the adjacency and distance matrices themselves. *Needs* ``chython[ml]``
- ML views: a molecule or a mapped reaction as ``int32`` numpy arrays -- element, hydrogen and degree
  columns per atom, a topological distance block, and a reaction's two sides in one union. *Needs*
  ``chython[ml]``
- Template-based reaction application from SMIRKS, functional-group detection, and protecting-group
  detection and removal
- Stereochemistry: tetrahedral, cis-trans, allene, atropisomer and helical units, with CIP labels
- 2D layout and depiction, with Jupyter support
- Interoperability with RDKit, CDK, Indigo, OpenBabel and CDPKit

Input is treated as unreliable by default: a reader stores and logs what a file says, and repair is a
pipeline you run afterwards. See :doc:`standardize`.


Cookbook
--------

.. toctree::
   :maxdepth: 2

   io
   pach
   molecule
   standardize
   substructure
   reactions
   ml
   glossary
   depiction
   config


Links
------

- `Source code <https://github.com/chython/chython>`_
- `PyPI <https://pypi.org/project/chython/>`_
- `Issues <https://github.com/chython/chython/issues>`_

Chython is a fork of `CGRtools <https://github.com/stsouko/CGRtools>`_.
