Chython [ˈkʌɪθ(ə)n]
===================

Library for processing molecules and reactions in python way.

Features:
   - Read/write/convert formats: MDL .RDF (.RXN) and .SDF (.MOL), .MRV, SMILES, INCHI (Linux and Windows), .XYZ, .PDB
   - Standardize molecules and reactions and valid structures checker
   - Supported python-magic
   - Tetrahedron, Allene and CIS-TRANS stereo supporting
   - Perform subgraph search
   - Build/edit molecules and reactions
   - Produce template based reactions and molecules
   - Atom-to-atom mapping checker and rule-based fixer
   - Perform MCS search
   - 2d coordinates generation (based on `SmilesDrawer <https://github.com/reymond-group/smilesDrawer>`_)
   - 2d/3d depiction
   - Produce CGRs (Condensed Graph of Reaction)

Full documentation can be found `here <https://chython.readthedocs.io>`_.

Chython is fork of `CGRtools <https://github.com/stsouko/CGRtools>`_.


INSTALL
=======

Highly recommended to use python 3.8+. Python 3.6 and 3.7 deprecated.

Note: for using `clean2d` install NodeJS into system.

* **stable version available through PyPI**::

    pip install chython

* Install CGRtools library DEV version for features that are not well tested::

    pip install -U git+https://github.com/chython/chython.git@master#egg=chython

TESTS
=====

Run unit tests::

    git clone https://github.com/chython/chython.git && cd chython  # skip if already got sources
    pip install -e .
    pytest --pyargs chython

COPYRIGHT
=========

* 2014-2021 Ramil Nugmanov nougmanoff@protonmail.com main developer

CONTRIBUTORS
============

* Dinar Batyrshin batyrshin-dinar@mail.ru
* Timur Gimadiev timur.gimadiev@gmail.com
* Adelia Fatykhova adelik21979@gmail.com
* Tagir Akhmetshin tagirshin@gmail.com
* Ravil Mukhametgaleev sonic-mc@mail.ru
