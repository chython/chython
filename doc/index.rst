Chython
=======

Library for processing molecules and reactions in Python.

.. code-block:: bash

   pip install chython

**Key capabilities:**

- Parse and write SMILES, InChI, IUPAC, MDL (SDF/RDF/MOL), MRV, XYZ, PDB
- Substructure and MCS search with SMARTS support
- Standardize, canonicalize, and enumerate tautomers
- Morgan and linear fingerprints with Tanimoto similarity
- Atom-to-atom mapping (neural + rule-based)
- Template-based reaction application and deprotection
- Stereo handling (tetrahedral, allene, cis-trans)
- 2D/3D depiction with Jupyter support
- RDKit interoperability


Cookbook
--------

.. toctree::
   :maxdepth: 2

   io
   molecule
   standardize
   substructure
   reactions
   depiction
   config


Links
------

- `Source code <https://github.com/chython/chython>`_
- `PyPI <https://pypi.org/project/chython/>`_
- `Issues <https://github.com/chython/chython/issues>`_

Chython is a fork of `CGRtools <https://github.com/stsouko/CGRtools>`_.
