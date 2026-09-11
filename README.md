<p align="center">
  <img src="docs/logo256.png" width="256" alt="chython logo"/>
</p>

<h1 align="center">Chython [ˈkʌɪθ(ə)n]</h1>

<p align="center">
  <a href="https://pypi.org/project/chython/"><img src="https://img.shields.io/pypi/v/chython.svg" alt="PyPI version"/></a>
  <a href="https://pypi.org/project/chython/"><img src="https://img.shields.io/pypi/pyversions/chython.svg" alt="Python versions"/></a>
  <a href="https://github.com/chython/chython/blob/master/LICENSE"><img src="https://img.shields.io/pypi/l/chython.svg" alt="License: LGPLv3"/></a>
  <a href="https://chython.readthedocs.io"><img src="https://img.shields.io/readthedocs/chython.svg" alt="Documentation"/></a>
  <a href="https://app.codecov.io/gh/chython/chython"><img src="https://img.shields.io/codecov/c/github/chython/chython/master?label=python%20coverage" alt="Coverage of the Python layers"/></a>
</p>

Library for processing molecules and reactions in a Python way.

## Features

**File formats**

- Read and write MDL RDF/RXN and SDF/MOL (V2000 and V3000, including atom parity and enhanced stereo), Marvin MRV, CML, SMILES, and InChI with InChIKey (InChI Trust library)
- Read SMARTS and SMIRKS, Tripos MOL2, PDBx/mmCIF, legacy PDB, XYZ, and IUPAC names through OPSIN
- Compact binary (de)serialization — `pach` for a wire record, `mol.to_bytes()` for the arena buffer — and full pickle support
- A coordinate format hands back a record of atoms and coordinates rather than a molecule, because it states no bond order: `build_molecule()` places the atoms, `perceive_bonds()` reads the connectivity out of a geometry, `saturate()` raises the orders, and all three are calls the caller makes

**Input is unreliable by default**

A reader stores and logs what a file says — an illegal valence, a nonsense charge, an underivable
hydrogen count — and never rejects a record for being chemically wrong. Repair is a pipeline you run
afterwards: `kekule()`, `standardize()`, `fix_resonance()`, `thiele()`.

**Toolkit interoperability**

Conversions build the target structure directly from the graph, so atom order matches
`atoms()` and stereo is carried over without needing a 2D layout. Each toolkit has one callable in
`chython.interop` that dispatches on its argument, and the export direction is also a container method.

| Toolkit | API | Requires |
|---------|-----|----------|
| RDKit | `mol.to_rdkit()`, `rxn.to_rdkit()`, `chython.interop.rdkit()` both ways | extra `rdkit` |
| Open Babel | `mol.to_openbabel()` | extra `extra-clean2d` |
| Indigo | `mol.to_indigo()` | extra `extra-clean2d` |
| CDK | `mol.to_cdk()` | extra `extra-clean2d` + `cdk.jar` (`CDK_PATH`) |
| CDPKit | `mol.to_cdpkit()`, and 3D conformers (`conformer_engine = 'cdpkit'`) | extra `extra-clean3d` |

RDKit is the only one of the five with a reaction form.

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

- Atoms and bonds in one contiguous buffer, addressed by an id that survives editing; edits go through
  an explicit session (`with mol.edit() as e:`) and derived data is recomputed on seal
- Standardize, canonicalize, kekulize/aromatize, repair resonance forms, neutralize, put a mobile
  hydrogen and charge where the canonical order says, check valences
- Split and decompose salts, expand contracted groups, derive implicit hydrogen counts
- Many 3D models per molecule in one conformer store
- Tetrahedral, cis-trans, allene, atropisomer and helical stereo, with CIP labels
- Descriptors: TPSA, Crippen logP and MR, hydrogen-bond donors and acceptors, rotatable bonds, ring
  counts, Bertz CT, Randić and Zagreb indices
- The 166 MACCS structural keys, one-based, and QED with its three published weight sets — both state
  what they transcribe and neither claims parity with another implementation's bits or score
- Morgan and linear fingerprints with Tanimoto similarity, as hash sets, folded bit vectors or count
  vectors, and the graph matrices and distance-derived descriptors — extra `ml`
- A molecule or a mapped reaction as `int32` arrays for a model: `mol.state_view()`,
  `rxn.transition_view()` — extra `ml`

**Search**

- Subgraph isomorphism
- SMARTS parser with chython-specific query semantics, including component grouping for intramolecular
  patterns

**Reactions**

- Template application from SMIRKS: `mol.react(template)`, or `mol @ other` for a two-component join
- Reaction enumeration over the shipped reaction corpus
- Functional and protective group detection and deprotection, with the whole corpus in the docs
- Sticky fragment / linker enumeration for combinatorial reassembly
- Atom-to-atom mapping reconstruction against a template corpus: `rxn.reconstruct_mapping()`
- Reaction-level standardization: each molecule pass once per molecule, plus the passes a loop cannot do

**Depiction**

- 2D coordinate generation, default [SmilesDrawer](https://github.com/reymond-group/smilesDrawer), switchable to RDKit/CDK/Open Babel/Indigo (`clean2d_engine`)
- SVG and SVGZ output with Jupyter support, and scalar data overlaid on the same scene
- 3D: a stored conformer as an X3DOM document (`mol.depict3d()`) or a notebook widget (`mol.view3d()`)
- 3D conformer generation with RDKit or CDPKit (`conformer_engine`)

Full documentation can be found [here](https://chython.readthedocs.io).

## Install

Only Python 3.10+.

```bash
pip install chython
```

The default 2D layout backend needs no extra: its JS engine (QuickJS) is a required dependency and
costs under 2.5 MB.

A plain install has **no numpy**, and that is deliberate — the base install is about 7.5 MB of runtime
files, small enough for a serverless bundle. Reading and writing every format, `standardize()`,
`kekule()`, `thiele()`, `canonicalize()`, stereo, substructure matching, template application and
depiction all work without it. What needs `chython[ml]` is the surface that answers a numpy array: the
fingerprints, `atom_invariants`, `adjacency_matrix`, `distance_matrix`, the distance-derived graph
descriptors, `pharmacophore_invariants`, `maccs_keys()`/`maccs_bit_set()` and the ML views. Each of
those raises an `ImportError` naming the extra when you call it, so nothing fails at import time.

Optional extras, combinable (`chython[rdkit,iupac]`):

| Extra | Enables |
|-------|---------|
| `ml` | numpy, and with it fingerprints, atom invariants, the graph matrices, the descriptors built on them and the ML views |
| `rdkit` | RDKit conversion both ways, RDKit 2D layout and 3D conformers |
| `iupac` | `molecule.iupac` name generation (Python >= 3.11) |
| `extra-clean2d` | CDK, Open Babel and Indigo backends (CDK also needs `cdk.jar`) |
| `extra-clean3d` | CDPKit conformer engine |

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
