<p align="center">
  <img src="doc/logo256.png" width="256" alt="chython logo"/>
</p>

<h1 align="center">Chython [ˈkʌɪθ(ə)n]</h1>

<p align="center">
  <a href="https://pypi.org/project/chython/"><img src="https://img.shields.io/pypi/v/chython.svg" alt="PyPI version"/></a>
  <a href="https://pypi.org/project/chython/"><img src="https://img.shields.io/pypi/pyversions/chython.svg" alt="Python versions"/></a>
  <a href="https://github.com/chython/chython/blob/master/LICENSE"><img src="https://img.shields.io/pypi/l/chython.svg" alt="License: LGPLv3"/></a>
  <a href="https://chython.readthedocs.io"><img src="https://img.shields.io/readthedocs/chython.svg" alt="Documentation"/></a>
</p>

Library for processing molecules and reactions in a Python way.

## Features

**I/O and interoperability**

- Read/write/convert MDL RDF/RXN and SDF/MOL (V2000 & V3000), MRV, SMILES/SMARTS, InChI (inchi-trust library), IUPAC names (OPSIN), XYZ, PDB
- Compact binary (de)serialization and full pickle support
- RDKit interoperability (`to_rdkit` / `from_rdkit`)

**Molecules**

- Build and edit molecules and reactions with a pure-Python API
- Standardize, canonicalize, kekulize/aromatize, and check valences
- Tautomer enumeration
- Tetrahedral, allene, and cis/trans stereo
- Morgan and linear fingerprints with Tanimoto similarity

**Search**

- Subgraph isomorphism and maximum common substructure (MCS) search
- SMARTS parser with chython-specific query semantics

**Reactions**

- Condensed Graph of Reaction (CGR)
- Atom-to-atom mapping (neural attention + rule-based fixing) and template-based mapping reconstruction
- Template-based reaction application (`Reactor`, `Transformer`) and reaction enumeration
- Functional and protective group detection and deprotection
- Sticky fragment / linker enumeration for combinatorial reassembly

**Depiction**

- 2D coordinate generation (based on [SmilesDrawer](https://github.com/reymond-group/smilesDrawer)/RDKit/CDK/Open Babel/Indigo backends)
- 2D/3D SVG depiction with Jupyter support

Full documentation can be found [here](https://chython.readthedocs.io).

## Install

Only Python 3.10+.

```bash
pip install chython[racer-default]
```

## CGRtools

Chython is a fork of [CGRtools](https://github.com/stsouko/CGRtools).

## Copyright

- 2014-2026 Ramil Nugmanov <nougmanoff@protonmail.com> main developer

## Contributors

CGRtools contributors are included too.

- Adelia Fatykhova <adelik21979@gmail.com>
- Aigul Khakimova
- Aleksandr Sizov <murkyrussian@gmail.com>
- Alexandre Varnek <varnek@unistra.fr>
- Dinar Batyrshin <batyrshin-dinar@mail.ru>
- Dmitrij Zanadvornykh <zandmitrij@gmail.com>
- Philippe Gantzer
- Ravil Mukhametgaleev <sonic-mc@mail.ru>
- Tagir Akhmetshin <tagirshin@gmail.com>
- Timur Gimadiev <timur.gimadiev@gmail.com>
- Timur Madzhidov <tmadzhidov@gmail.com>
- Zarina Ibragimova
