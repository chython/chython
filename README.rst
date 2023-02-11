Chython [ˈkʌɪθ(ə)n]
===================

Library for processing molecules and reactions in python way.

Features:
   - Read/write/convert formats: MDL .RDF (.RXN) and .SDF (.MOL), .MRV, SMILES, INCHI (inchi-trust library), .XYZ, .PDB
   - Standardize molecules and reactions and valid structures checker
   - Supported python-magic
   - Tetrahedron, Allene and CIS-TRANS stereo supported
   - Perform subgraph search
   - Build/edit molecules and reactions with Python API
   - Produce template based reactions and molecules
   - Atom-to-atom mapping, checking and rule-based fixing
   - Perform MCS search
   - 2d coordinates generation (based on `SmilesDrawer <https://github.com/reymond-group/smilesDrawer>`_)
   - 2d/3d depiction with Jupyter support
   - SMARTS parser with restrictions
   - Protective groups remover
   - Common reaction templates collection

Full documentation can be found `here <https://chython.readthedocs.io>`_.

Chython is fork of `CGRtools <https://github.com/stsouko/CGRtools>`_.

Install
=======

Only python 3.8+.

Note: for using `clean2d` install NodeJS into system.

* **stable version available through PyPI**::

    pip install chython

* Install chython library DEV version for features that are not well tested::

    pip install -U git+https://github.com/chython/chython.git@master#egg=chython

Copyright
=========

* 2014-2023 Ramil Nugmanov nougmanoff@protonmail.com main developer

Contributors
============

* Adelia Fatykhova adelik21979@gmail.com
* Aleksandr Sizov murkyrussian@gmail.com
* Dinar Batyrshin batyrshin-dinar@mail.ru
* Dmitrij Zanadvornykh zandmitrij@gmail.com
* Ravil Mukhametgaleev sonic-mc@mail.ru
* Tagir Akhmetshin tagirshin@gmail.com
* Timur Gimadiev timur.gimadiev@gmail.com
* Zarina Ibragimova
