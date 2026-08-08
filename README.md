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

**File formats**

- Read/write/convert MDL RDF/RXN and SDF/MOL (V2000 & V3000, including atom parity and enhanced stereo), MRV, SMILES/SMARTS, InChI (inchi-trust library), XYZ, PDB
- Compact binary (de)serialization and full pickle support

**Toolkit interoperability**

Conversions build the target structure directly from the graph, so atom order matches
`atoms()` and stereo is carried over without needing a 2D layout.

| Toolkit | API | Requires |
|---------|-----|----------|
| RDKit | `to_rdkit()` / `MoleculeContainer.from_rdkit()` | extra `rdkit` |
| Open Babel | `to_openbabel()` | extra `extra-clean2d` |
| Indigo | `to_indigo()` | extra `extra-clean2d` |
| CDK | `to_cdk()` | extra `extra-clean2d` + `cdk.jar` (`CDK_PATH`) |
| CDPKit | 3D conformers (`conformer_engine = 'cdpkit'`) | extra `extra-clean3d` |

Allene stereo is not portable through any of these toolkits. Indigo additionally omits
cis-trans, which it derives from 2D coordinates.

**IUPAC names, both directions**

```python
from chython import iupac

mol = iupac('ethanol')   # name -> structure, via OPSIN
mol.iupac                # 'ethanol' -- structure -> name, via openclatura
```

`iupac()` needs JPype and `opsin.jar` (`OPSIN_PATH`); the `.iupac` property needs the
`iupac` extra (Python >= 3.11) and returns `None` when the structure cannot be named.

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

- 2D coordinate generation, default [SmilesDrawer](https://github.com/reymond-group/smilesDrawer), switchable to RDKit/CDK/Open Babel/Indigo (`clean2d_engine`)
- 3D conformer generation with RDKit or CDPKit (`conformer_engine`)
- 2D/3D SVG depiction with Jupyter support

Full documentation can be found [here](https://chython.readthedocs.io).

## Install

Only Python 3.10+.

```bash
pip install chython[racer-default]
```

Optional extras, combinable (`chython[racer-default,rdkit,iupac]`):

| Extra | Enables |
|-------|---------|
| `racer-default` | JS engine for the default 2D layout backend |
| `rdkit` | `to_rdkit`/`from_rdkit`, RDKit 2D layout and 3D conformers |
| `mapping` | neural atom-to-atom mapping |
| `iupac` | `molecule.iupac` name generation (Python >= 3.11) |
| `extra-clean2d` | CDK, Open Babel and Indigo backends (CDK also needs `cdk.jar`) |
| `extra-clean3d` | CDPKit conformer engine |
| `png` | PNG output from depiction |

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
