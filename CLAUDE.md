# CLAUDE.md

Guidance for Claude Code (claude.ai/code) working in this repository.

## Project Overview

Chython is a Python library (LGPLv3) for processing molecules and reactions, historically a fork of
CGRtools. This branch is **chython 3**, a rewrite around a compiled arena-based core.

**chython 2 IS DELETED.** `chython/algorithms/`, `chython/containers/`, `chython/files/`,
`chython/reactor/`, `chython/periodictable/` and `chython/utils/` no longer exist, including the modules
that were never ported. Do not restore one to read a table or a mixin out of it — three ratchets fail if
any of those names becomes importable again (`chython/test/test_v2_boundary.py`,
`chython/formats/test/test_isolation.py`, `chython/interop/test/test_config.py`). Read them out of git:

```bash
git show 5e39eb5:chython/algorithms/tautomers/__init__.py   # any V2 file, at its last living commit
```

`docs/superpowers/research/2026-09-03-v2-ported-audit.md` is the per-file PORTED / PARTIAL / NOT PORTED
table with each V3 home. **Read it before assuming a V2 feature exists in V3**, and know it lags:
tautomers, MCS, CGR and *neural* atom-atom mapping are genuinely absent, while the reactor, mapping
reconstruction, isomer standardization, fingerprints, salts, PDB and MRV-as-an-XML-dialect have landed.

## Build & Development

Standard PEP 517 setuptools. Requires Cython >= 3.1 (for the `freethreading_compatible` directive).

```bash
pip install -e . --group dev        # editable install + dev group (PEP 735); or: uv sync --group dev
python setup.py build_ext --inplace  # rebuild in place — the developer loop
python -m build --wheel
```

`setup.py` holds only what `pyproject.toml` cannot compute: platform flags and the `Extension` list.
**That list has exactly one entry, `chython.core._core`**, ratcheted by `test_v2_boundary.py` — a second
extension needs a declared reason. `build_inchi.py` builds libinchi from the `INCHI` submodule into
`build/inchi/` and stages it into `chython/core/`, warning and skipping when cmake or the submodule is
absent. Runtime data files are named explicitly in `[tool.setuptools.package-data]` with
`include-package-data = false` — nothing is globbed and `test_packaging.py` fails if you forget one.

## Testing

```bash
pytest chython/                    # all of it — about 8700 tests in ~3 min
pytest chython/core/test/          # per package: core, chemistry, reactions, formats/*, depict, interop
pytest chython/test/               # tree-wide: packaging, libinchi staging, V2 boundary, doc samples
```

Tests live in `test/` subdirectories per package; no centralized conftest; data files in `/test/` at the
repo root. **Differential tests run chython 2.24 out of process** — `core/test/oracle.py` spawns a pinned
separate install under `-I`, the only sanctioned way to consult V2 since a subprocess is not an import;
they skip when that venv is absent. **Every Python sample in `docs/` is executed**: `test_doc_samples.py`
runs each `.. testcode::` block and forbids `.. code-block:: python`, which renders identically and runs
never.

## Architecture

```
core  <-  chemistry   reactions        formats   depict   interop
```

- **`core/`** — the one Cython extension: containers, arena storage, SMILES/SMARTS read and write,
  isomorphism, rings, kekule/thiele, canonical form, stereo, fingerprints, InChI, pach.
- **`chemistry/`** — chemical knowledge *and* the passes applying it: TSV tables in `tables/`, loaders in
  `_tables.py`, and `standardize`, `fix_resonance`, `calc_implicit`, `check_valence`, `saturate` beside
  them. `chython/transformation/` is **gone**; `chython.transformation.X` is now `chython.chemistry.X`.
- **`reactions/`** — SMIRKS templates, the reaction/functional/protective corpora, mapping reconstruction.
  A **sibling** of `chemistry`, not a layer below. Then **`formats/`**, **`depict/`**, **`interop/`**.

Two rules, ratcheted by `chython/chemistry/test/test_dependency_direction.py` and
`chython/formats/test/test_isolation.py`:

1. **Nothing imports the façade.** `chython/__init__.py` must never appear in an import inside `core`,
   `chemistry`, `reactions` or `formats` — production code or tests. The layers are provable in isolation,
   which is what keeps `lazy_object_proxy` out of the façade.
2. **The direction never reverses.** `chemistry` may see `core`; nothing below sees anything above. Inside
   `chemistry` a third rule holds per module: `_tables.py` and `_smarts.py` read the tables and may not
   import the passes, so the knowledge half stays reviewable alone.

Registration onto the sealed core container is by **injection**, not inheritance — `chemistry` calls
`_set_standardize_fn`, `reactions` calls `_set_reactions_fns`. `MoleculeContainer` is a `cdef class`: no
mixins, no `Graph` base, and a special method must be compiled in whatever package supplies its body.

### Core Container Model

Atoms and bonds live in a **single contiguous arena buffer**, not in dicts — no `_atoms: Dict[int, Atom]`
and no adjacency dict. Atoms are addressed by a **stable id** that survives edits; `map_number` is a
**separate** `uint16_t` field — so a V2 `remap()` is a V3 `set_map_number()` when it meant atom-atom
mapping, and nothing at all when it meant renumbering. Mutation goes through an explicit edit session
(`with molecule.edit() as e:`) and derived data is recomputed on seal.
**`chython/core/RULES.md` is the binding coding standard for anything in `core/` — read it before editing
a `.pxi`.** The core is **one translation unit**: `_core.pyx` `include`s each `.pxi`, so adding a layer
means adding an `include` line, not an `Extension`. `ReactionContainer.molecules()` yields
reactants→agents→products and that order is load-bearing; **`CGRContainer` does not exist in V3** and is
not planned — `ReactionModelingView` is the CGR-free view chytorch consumes.

Element data, isotopes and valence rules are **generated C tables** compiled into `_core` from
`core/elements.tsv`, `isotopes.tsv` and `valence_rules.tsv`; `core/test/` re-runs each generator's
`compile` verb against its TSV, so they cannot drift. There is no `periodictable/` and no
metaclass-generated element class. Edit the TSV and regenerate; **never hand-edit a generated table.**

### Rule Tables

`chython/chemistry/tables/` and `chython/reactions/tables/` hold knowledge as TSV. When adding one:

- put the TSV in `tables/` and name it in `[tool.setuptools.package-data]` **with the `tables/` prefix** —
  the pattern is matched against the path relative to the package, so a bare filename ships nothing;
- a NamedTuple for the compiled row, a `read_table()` call, a module-level `_*_CACHE` dict and an accessor
  function — tables load lazily on first use, never at import;
- `lazy_object_proxy` is **banned** and ratcheted against; rule ids are table-qualified (`'groups:13'`);
- a pass takes `(molecule, ...)` and **no `log=`**, returns `bool`, records through
  `with recording(molecule, stage='<name>') as log:` — unconditionally, since `molecule.log` is the one
  destination and a `if log is not None` branch is a defect — and is **all-or-nothing** per site. Only
  readers, writers and `depict` take a `log=` list.

## Input Posture

**Input is garbage by default, no exceptions.** A reader never rejects a record for being chemically wrong:
an illegal valence, a nonsense charge, an underivable hydrogen count are *stored* and *logged*. An
underivable fact becomes a registered reserved *unknown* (`H_UNKNOWN`), never a silent zero and never
grounds for dropping the record. Repair is an **explicit** pipeline run afterwards — `kekule()`,
`standardize()`, `fix_resonance()`, `thiele()`. `kekule()` repairs by design (it will separate charges on
an aromatic N-oxide to find a Kekulé form) and logs it; `thiele()` is single-purpose and **refuses** rather
than repairing. Refusals live at the answer boundary and nowhere else. `chython.chemistry.saturate()` is
the separate, explicitly invoked bond perception for a file that gave connectivity and no orders, and it
refuses where nothing pins an order sum; XYZ, PDB and mmCIF therefore return a record of atoms and
coordinates rather than a `MoleculeContainer`.

**A parsed molecule already has its implicit hydrogen counts** — at read time, by **one** algorithm every
reader shares, `core/_hydrogens.pxi`. Do not write a second copy; there were three and they had all
drifted. Exactly one class of atom is left unknown, the pnictogen whose class the ring decides (pyrrole
versus pyridine), tested by `arom_classify_atom` answering `AROM_MAY` and not by a pattern or an element
list. `kekule()` closes that class itself, in `fill_only` mode and only when no aromatic system is
unresolved — so `canonicalize()` is literally its own stages run by hand, a pipeline never being stronger
than its parts. `calc_implicit` and `check_valence` delegate to the same derivation.

## Domain References

Two reference bodies in the shipped documentation are the authority there. Designs, including what is
specified but unimplemented, are under `docs/superpowers/specs/`.

- **`docs/substructure.rst` — matching and the chython SMARTS dialect**: primitive tables, `&`/`,`/`;`
  precedence, `z` 1–6, `*`/`^` charge and radical semantics, `[A]`/`[M]` neutrality, component grouping,
  and the two rules most often got backwards — **an implicit bond matches single only, never aromatic**,
  and **V2's `z3` is not V3's `z3`**, so re-read every `z3` in a ported template. One parser, in the core;
  the tables use it too.
- **`docs/reactions.rst` — SMIRKS templates, the corpora, mapping reconstruction**: `read_smirks` and the
  callable template, `mol.react()`/`mol @ other`, `functional_groups()`, `protective_groups()`/
  `deprotect()`, one table and one method for the whole corpus, slots and the `i * 100` map-number
  stride, `ring_sizes` as what makes a row intramolecular, and `reconstruct_mapping()`'s five-rung ladder
  and its refusal of a multi-product record.

## Key Conventions

- Python 3.10+ (`tomllib` is 3.11, which is why `test_packaging.py` parses TOML by regex)
- `__slots__` on Python-level classes; `cdef class` with declared attributes in `core`
- Caching by `functools.cached_property` and lazy `_CACHE` dicts; `CachedMethods` is a V2 dependency, gone
- No enforced linter, but `chython/test/test_code_hygiene.py` ratchets the part of PEP 8 that is mechanical
  — final newline, no CRLF, no trailing whitespace, no blank line holding spaces, no run of three, and 120
  columns **including the `.pxi` files `pycodestyle` does not read**. Widen `WIDE_BY_CONSTRUCTION` there,
  with a reason, rather than re-running a sweep. What is left to `pycodestyle --max-line-length=120` and to
  judgement: aligned assignment blocks (E221) and `chemistry/test/test_tpsa_tsv.py`'s hand-aligned probe
  table (E131), `depict/style.py`'s trailing comments continued under their column (E114/E116), and the
  documented mid-file imports in `interop/__init__.py` and `formats/ctfile/test/conftest.py` (E402)
- `chython/__init__.py` is a thin façade; `__all__` is empty by design. Engine and depiction configuration
  is documented in `docs/config.rst`. **`torch_device` is gone** — when the neural mapper lands it runs on
  ONNX Runtime, not torch
- Only `chython/` and `docs/` are source. Everything else at the repo root — `scripts/`, `mapping/`,
  `java/`, `pach/`, benchmark scripts, `test/` data — is temporary or data. Scope sweeps accordingly
- Leave `docs/superpowers/` specs and plans untracked; `.gitignore` covers the directory
- A coverage or agreement claim ships with the harness that produced it, or it is not made

## Writing Style

Comments, docstrings and committed `.md`/`.tsv` prose are **crisp and factual**, never narrative. State the
rule, then at most one short concrete instance that makes it actionable. What does not belong: postmortem
storytelling, tallies of past mistakes, chronology, rhetorical wind-up, or the same fact restated in prose
after a table already gave it. Prefer a table to a paragraph. A stale path or symbol name in a comment is a
defect — comments and docstrings that look like code get swept with the code. A `# --- label ---…` block
rule is padded to one column per file; `test_code_hygiene.py` holds that column where a file's rules already
share one, and leaves ragged widths alone because they are labels rather than a block.

**A V3 rule states itself.** No test can enforce this — `_pach.pxi` says "v2" 52 times about a *wire format*
version — so it is read in the diff:

| Cut | Keep |
| --- | --- |
| "used to", "previously", "formerly", "historically" | a semantic incompatibility a reader would assume away: **V2's `z3` is not V3's `z3`** |
| before/after tallies and scores | a renamed concept whose meaning moved: a V2 `remap()` is `set_map_number()` for mapping, nothing for renumbering |
| commit archaeology, deleted-test history | `git show 5e39eb5:<path>` as the way to read a deleted V2 file |
| a verdict on V2 — "conflated", "naive", "a defect" | the differential harness, where chython 2.24 is the subject under test |
| a V3 rule stated by contrast when it stands alone | an assertion message that names what a ratchet forbids |

## What Never Ships

These rules apply to comments, docstrings, test fixtures, assertion messages and generated files alike.
No test in the tree enforces them — a scanner would have to spell out what it forbids — so read the diff:

- **No employer name and no internal host.** Not any spelling or transposition of it, not in a URL. Live
  case: `npm install` behind a mirror rewrote all 29 `resolved` URLs in `clean2d/package-lock.json` to an
  internal artifact repository. It must resolve to `registry.npmjs.org`, and regenerating it behind a
  mirror reintroduces this — check the diff.
- **Never disparage another toolkit.** State capability as fact and nothing more — *"Indigo cannot
  represent `H_UNKNOWN`"*, *"RDKit 2026.03.4 refuses `%05`"*, and *"chython is 3.7× slower on TPSA"* are
  all fine. What is not: calling another tool's behaviour a guess, a bug, broken, naive or wrong. Where a
  spec is silent and implementations differ, say the spec is silent and that readers differ. Benchmarks
  compare by measurement, never by verdict, and the reference toolkit's *fastest* spelling is the one to
  quote.
- **No internal corpus.** Test structures are public compounds. A real-looking scaffold with no citation
  is treated as internal until shown otherwise.

## Copyright Headers

When modifying files, update Ramil Nugmanov's copyright year to include 2026:

- **Keep the first year** of development (from git history or existing header)
- **Range for 3+ years**: `2019-2026`; **comma for exactly 2**: `2025, 2026`; **new file**: `2026`
- **Never drop other contributors** — all co-author copyright lines must remain unchanged
- **Only update years for Ramil** — don't modify other contributors' year ranges

```
Copyright 2019-2024 Ramil Nugmanov  →  Copyright 2019-2026 Ramil Nugmanov
Copyright 2025 Ramil Nugmanov       →  Copyright 2025, 2026 Ramil Nugmanov
Copyright 2023, 2024 Ramil Nugmanov →  Copyright 2023-2026 Ramil Nugmanov
```
