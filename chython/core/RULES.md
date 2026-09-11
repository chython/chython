# chython/core — coding rules

The binding coding standard for new and modified code in `chython/core/`.  Each rule carries one
short instance so it can be applied rather than interpreted.  A rule whose basis is a language
invariant rather than a measurement says so.

---

## 1. File naming

### 1.1 One file, one primary definition

**Rule:** A file is named for the one thing it primarily defines.  If naming a file takes a
conjunction (`_arena_and_features.pxi`), it holds two things and must be split.

### 1.2 Layout

The core's structural files:

```
_molecule_arena.pxi   _molecule_container.pxi   _molecule_views.pxi
_features.pxi         _elements.pxi             _query_arena.pxi
_query_boxes.pxi      _query_seal.pxi           _query_container.pxi
_rings.pxi            _sssr.pxi                 _morgan.pxi
_isomorphism.pxi
```

This matches the table in §1.3, which is the authority.

### 1.3 Naming table

| role | molecule side | query side |
|---|---|---|
| sealed arena | `_molecule_arena.pxi` (defines `Structure`) | `_query_arena.pxi` (defines `Query`) |
| primitive → box → DNF compiler | — | `_query_boxes.pxi` |
| builder journal → sealed arena | inside the container (`journal_t`, `_apply`) | `_query_seal.pxi` |
| Python-facing class | `_molecule_container.pxi` | `_query_container.pxi` |
| per-element views | `_molecule_views.pxi` (`Atom`, `Bond`) | — deliberately none |
| shared vocabulary | `_features.pxi` — the 4×u64 feature-word encoding, written by the molecule side and read by the query side | |
| standalone algorithms | `_elements.pxi`, `_rings.pxi`, `_sssr.pxi`, `_morgan.pxi`, `_isomorphism.pxi` | |

Four entries in that table are choices, not gaps:

- **The query side has no views file.**  A query needs fast matching, not a rich property API, so a
  `_query_views.pxi` would be wrong rather than merely empty.
- **The molecule side has no separate seal file.**  The builder journal lives inside `MoleculeContainer`
  and applies in place (`_apply`).  If it outgrows the container it becomes `_molecule_seal.pxi`.
- **Side-prefixed names wherever both sides have a counterpart**, because the parallel architecture must
  be visible in the file names: `_molecule_arena.pxi` / `_query_arena.pxi`, `_molecule_container.pxi` /
  `_query_container.pxi`.  Standalone algorithms keep single-subject names, naming a subject rather than
  a role.
- **`Structure` is not renamed `Molecule`.**  `Molecule` beside `MoleculeContainer` does not say which is
  the arena.  The file name carries the alignment, and each file's header comment states the primary type
  it defines.

### 1.4 The S-group segment, split by what an EDIT can invalidate

The split is **not** structured versus unstructured; it is *invalidatable* versus *opaque*:

| arena-native, because an edit can break it | opaque payload, stored without interpretation |
|---|---|
| `atoms`, `patoms` — stable ids | `type`, `subtype`, `name`, `disp` |
| `bonds` — ENDPOINT PAIRS | `fields` — unmodelled keyword → list of values, order significant |
| `cstates` — `((n, m), tail)` | `data` — list of raw **bytes** |
| | `log` |

`aliases` (V2000 `A  <n>` display labels) sits on the store as a stable-id-keyed map, not per record.

**Endpoint pairs, never bond indices.**  A bond index is a position, and a position does not survive
a reordered write — the same reason the arena addresses atoms by stable id.

Five invariants, each a test rather than a hope:

1. **A record with zero atom references SURVIVES.**  An empty DAT record still asserts that a field
   was attached to something; this is §6.3's empty-versus-absent distinction, and dropping such
   records breaks the reader's fidelity guarantee silently.
2. **On compaction, references are REMAPPED**, and a reference to a genuinely deleted atom is dropped
   **and reported** — never silently zeroed.  Id 0 must not become reachable by omission.
3. **Record order is file order**, and `fields` values keep their order within a key.
4. **`data` is bytes in, bytes out.**  Undecodable bytes are the fidelity requirement; a decode at
   store time loses them permanently.
5. **`NO_INDEX = 0xFFFF` is a sentinel for three fields** (`index`, `ext_index`, `parent`), so a
   `uint16_t` holding them has a real maximum of `0xFFFE`.  §6.3.

**A slot wider than the field it serves must not be range-checked as if it were the field.**  `xy_t` is
`int32_t` scaled by 10000 — exactly the fixed-point grid of a ten-column `F10.4` field, but not its
range: ten columns reach `99999.9999` upward and only `-9999.9999` downward, because the sign costs a
column.  `-214748.3647` fits the slot and cannot be written to the field, and a writer trusting the
slot's domain emits an eleven-column number, which shifts every later column so the *next* value
reparses as garbage instead of failing.  The arena keeps the wider slot — clamping on the way in loses a
coordinate the input stated — and the field's range is the writer's to enforce at emission.

### 1.5 The CIP descriptor fields, and what an EDIT can invalidate

Storage only.  Nothing in the core assigns a descriptor.  An atom's descriptor is the low nibble of
`atom_t.reserved` (9 values: none, `R`, `S`, `r`, `s`, `M`, `P`, `m`, `p`), a bond's is three bits of
`halfedge_t.flags` (5 values: none, `E`, `Z`, `M`, `P`).  Both words were already serialised, which
is why they were chosen: `to_bytes()` is a molecule identity in v4, so widening `atom_t` would
reprice every stored key.

**A stored descriptor is a claim the INPUT made about the molecule, not a function of the fields the
arena holds.**  Everything below follows from that.

1. **Case is meaning.**  `r`/`s` are the pseudo-asymmetric descriptors, a different determination about
   a different kind of centre, not spellings of `R`/`S`.  Nothing on the path calls `.upper()`, and
   `'e'` is refused rather than read as `'E'`, because the same case-folding turns `R` into `r`.
2. **Two domains, two tables.**  `M` and `P` are axial descriptors on both sides, so one merged table
   would give a bond descriptor sent to an atom silence instead of an error.
3. **A bond descriptor is one statement about one bond.**  Both half-edges are written from a single
   call site and either end reads the same answer.  Unlike a wedge, it is not directional.
4. **The drop rule keys on the OPERATION, not on the field.**  An edit changing which atoms exist,
   which are bonded, or a bond's order drops every descriptor the molecule held and logs a count.
   Charge, isotope, radical, map number, hydrogen count, coordinates, wedges, stereo flags and stereo
   groups keep them — isotope included, even though CIP Rule 2 ranks by mass, because this layer did
   not compute the stored descriptor.
5. **`kekule` and `thiele` are exempt, and the exemption list is closed.**  They are the only
   operations allowed to change a representation, and the aromatic and Kekulé spellings are one
   molecule.  They reach the journal as ordinary order changes, so they announce themselves with a flag
   nothing else sets.  A caller's own `set_order` is not exempt.
6. **The scope that states a descriptor wins, wherever in the scope it states it.**  The drop clears
   pre-scope descriptors and the replay applies the scope's own, so a parser emitting atoms, bonds and
   descriptors in file order need not care where in its scope a descriptor landed.
7. **A drop is recoverable only from the log.**  Storage cannot tell a never-labelled atom from a
   dropped one — both hold code 0 — so a labelled molecule and its unlabelled twin are byte-identical
   after an invalidating edit.  Same shape as `sgroup_log`.  `substructure` is the one stated exception:
   a cut is a caller asking for a smaller molecule, and its docstring says so.
8. **A descriptor is in the BYTES and out of the CANONICAL FORM.**  `to_bytes()` differs between a
   labelled molecule and its unlabelled twin; `==`, `hash` and `atoms_order` do not.  A dict keyed on
   molecules must not grow a second entry because someone annotated one of them.

**Ruling — the autolabeler must RECOMPUTE, must not read a stored descriptor as an input, and must
handle the aromatic form directly.**  The assignment algorithm is later work; these three are decided
now because each is cheaper to state than to retract.

* *Recompute.*  A surviving descriptor's only guarantee is that nothing has invalidated it since an
  input file stated it, and input is garbage by default.
* *Do not read one as an input* — not as a hint, a tie-break or a cache.  The moment a stored descriptor
  can influence a computed one, the drop rule becomes a correctness dependency and every exemption in
  item 4 has to be re-argued as a chemical claim rather than a storage one.
* *Handle the aromatic form directly, do not kekulise first.*  CIP treats a mancude ring system with
  **averaged** duplicate atomic numbers at the duplicated positions, so a descriptor computed from an
  arbitrary Kekulé form depends on which form was picked — and passes any test written against that same
  form.

---

## 2. Struct pointer binding

### 2.1 Rule

Bind a struct pointer once, then read members through it:

```cython
cdef journal_t *rec = self._journal + n
rec.op = ...
rec.a  = ...
```

Do not repeat `self._journal[n].` for each field of the same record.

**The threshold is two accesses to the same path in one basic block, reads included.**  That covers
writing four fields of one record, reading two fields of one record, testing one field twice in an
`if`/`elif` chain, and reading the same field at several subscripts (`box.neg[0..3]`).  A path touched
exactly once per basic block stays as written.

Multi-level walks matter most, because each level multiplies: `atom_terms[s].boxes[k].neg[0]` appeared
eleven times in twelve lines of one loop.  Bind the innermost record you use twice, not the outermost:

```cython
box = &term.boxes[k]          # not `wt = &atom_terms[s]` alone
for u in range(4):
    box.neg[u] = src.neg[u]
```

The rule is not about speed: `tokens[i].opcode = x` on a typed pointer already compiles to a plain
offset store, so the repeated form trips no performance alarm.  The cost is in §2.4.

### 2.2 Evidence — `_molecule_container.pxi`

`_append` (`:181–185`) once wrote four repeated index expressions for one 16-byte record
(`self._journal[self._journal_len].op = op`, and three more); it now binds
`cdef journal_t *rec = self._journal + self._journal_len` and writes four members.

`_apply` binds `cdef journal_t *jr = self._journal` at `:233` and `rec = jr + i` at the top of the
replay loop body (`:340`), replacing 23 `jr[i].` accesses.  Its pre-pass loop (`:235–245`) reads only
`jr[i].op` and is left alone.

`total_h_of` (`:556–558`) is the reference shape: `cdef atom_t *a = self._atom(sid)` bound once, because
`at_implicit_h(a) + at_explicit_h(a)` needs two fields.

### 2.3 Exception — never bind a pointer into a `packed` struct

`-Waddress-of-packed-member` fires on any address-of inside a `cdef packed struct`, whatever the
provable alignment.  This branch tracks the build's warning table as an invariant — currently
**five** `-Wunreachable-code` sites; the gate is "the table is unchanged", not "the build is silent",
and it covers every line matching `warn`, Cython's own included.  Never quote the count from memory;
**run §9.7's gate**, whose first clause is `rm -f chython/core/_core.c &&` — a plain `build_ext
--inplace` on a warm tree re-runs neither Cython nor clang and prints nothing, so it reads as a silent
build rather than as a gate that ran.  **Hoist the values instead**, as `_molecule_arena.pxi:from_bytes`
does:

```cython
# segment_t lives inside `cdef packed struct StructureHeader`
seg_off = structure.header.segments[j].offset      # read the values out
seg_len = structure.header.segments[j].length
```

The exception covers the address of a *member* inside a packed struct.  The address of a *whole*
packed struct is always valid — alignment 1 — so `&out_boxes[i]` on a `qbox_t *` is fine and
`_query_seal.pxi` relies on it.  In `chython/core` today the packed-member cases are `segment_t` inside
`StructureHeader` and `xy_t` inside `sgroup_t`, both in `_molecule_arena.pxi`, and `refs` inside
`stereo_unit_t` (`_stereo.pxi`), whose four entries `_pach3.pxi` hoists into a local array rather than
addressing; at the `xy_t` site
(`_molecule_container.pxi:_sgroup_dict`, cited by NAME because every edit above it moves the line — see
§10.3) the coordinate read divides in place rather than calling
`xy_read_x`/`xy_read_y`.  `atom_t`, `halfedge_t`, `journal_t`, `wedge_edit_t`, `wbox_t`, `qtoken_t`
and `qop_t` are unpacked and unrestricted; `qatom_t`, `qbox_t`, `qany_t`, `qbond_t`, `qclosure_t` and
`qcomp_t` are packed but only ever addressed whole.

### 2.4 Corollary — the repeated walk hides duplication, and that is the real cost

The verbose form is what lets two copies of one emitter look different enough to write twice.  Every
layer-3 extraction on this branch was found by hoisting a local and re-reading the neighbour, and
each pair had already survived review:

| the pair | what the walk hid |
|---|---|
| `_features.pxi:fill_features` / `fill_edge_words` | both built feature word 0's element span and topology/order triple; one documented the dormant `HE_AROMATIC` branch and the other carried it bare (now `w0_element_bits` + `w0_bond_bits`) |
| `_molecule_arena.pxi:csr_build` | the two halves of a bond, written out twice and required to agree on order, wedge and flags (now `_emit_half`) |
| `_query_boxes.pxi:boxes_merge` | the `diff_cnt == 0` and `diff_cnt == 1` arms shared a four-line tail (now written once) |
| `_query_seal.pxi:query_seal` | the atom-box and bond-box emitters were identical down to clearing `spare`; the sizing pass counted any-records twice (now `_emit_boxes` + `_term_any_total`) |

So the check is two-step: bind the pointer, then **re-read the neighbouring block**.  If it now looks
the same, it was the same — extract it.  A pair that must agree field for field belongs in one
function, where the agreement is structural rather than a promise in a docstring.

---

## 3. Do not rely on the optimiser to undo repetition

### 3.1 Rule

Never write repeated struct member accesses on the assumption that the compiler will collapse them.

### 3.2 Basis — a language rule, not a build measurement

Cython does not collapse the repeated `jr + i` arithmetic: before the pointer bind, `_apply` emitted 28
separate `(__pyx_v_jr[__pyx_v_i]).field` index expressions; after, the replay loop opens with one
`__pyx_v_rec = (__pyx_v_jr + __pyx_v_i)` and does 23 `__pyx_v_rec->field` member reads.  What was counted
is dereference *expressions* in generated C, not assembly loads; whether GCC or Clang removes them at
`-O2` was not measured.

The reason to write the explicit form is the language, not the build.  Cython emits no `restrict`, so a
store through any `uint8_t *` in scope may legally alias a struct field the compiler has already loaded,
and the standard permits reloading a field after any interleaved store.  Even where the optimiser does
collapse the loads, the repetition is noise the reader pays for on every review.

---

## 4. Struct layout: when AoS is right and when SoA is right

### 4.1 Rule

Fields that travel together (accessed in the same statement, or in immediately adjacent statements on
the same element) belong in one struct.  Fields scanned independently across many elements belong in
parallel arrays.

### 4.2 Evidence — `halfedge_t` is right as AoS

```cython
cdef packed struct halfedge_t:
    uint32_t to
    uint8_t  order
    uint8_t  wedge
    uint16_t flags
```

Eight bytes, packed.  Every CSR traversal reads `to` and at least one of `order`, `wedge`, `flags` in the
same pass; a bond is never visited without its order.

### 4.3 Evidence — feature words and query neg masks are right as flat arrays

The per-atom feature words (`structure_features`, returning `uint64_t *`) and the query box `neg` masks
(`qbox_t.neg[4]`) are 4×64-bit payloads scanned word by word: the isomorphism kernel tests
`f[1] & box.neg[1]` for three words in one `if` (`_isomorphism.pxi:_box_admits`), and `fill_features` writes
four consecutive words per atom.  Flat `uint64_t` arrays avoid the indirection a nested struct introduces
at those call sites and admit simple pointer arithmetic.

---

## 5. Scratch follows the loop; arena segments follow the format

### 5.1 Rule

A segment's layout is part of the serialised format and is round-tripped verbatim.  It must not be
reshaped to suit one loop's access pattern.  Scratch has no such constraint and should be reshaped to
whatever reduces allocation overhead.

### 5.2 Rule — one scratch struct, one malloc, one failure check, one free

When a function needs several scratch arrays, pack them into one struct, allocate the struct in one
`PyMem_Malloc`, check the one result, and free the one block in `finally`.  Size each region with
`align8()`.  `_molecule_arena.pxi`'s `structure_alloc_full` is the reference implementation.

### 5.3 Evidence — `_molecule_container.pxi:_apply`

`_apply` once allocated scratch in nine `PyMem_Malloc` calls, with a nine-clause NULL test and nine
matching `PyMem_Free` calls — a symptom of the nine allocations, not a style choice.  It now uses one
`apply_scratch_t` block (`:254–274`, aliases `:276-282`).  `esrc`/`edst`/`eord` became
`edge_edit_t *edits` (AoS); seven 8-aligned regions are carved from a single `PyMem_Malloc`; `wxy` and
`wsg` get zero bytes and a NULL pointer when their segment is not wanted:

```cython
cdef apply_sizes_t _sz = _scratch_sizes(n_max, e_max, w_max, want_xy, want_sg)
cdef apply_scratch_t scratch
scratch.block = PyMem_Malloc(_sz.work + _sz.live + _sz.newidx + _sz.edits
                             + _sz.wedge + _sz.wxy + _sz.wsg)
if scratch.block is NULL:
    self._discard()
    raise MemoryError('journal apply scratch allocation failed')
# ... carve seven typed pointers from scratch.block, then free the one block in `finally` ...
```

---

## 6. Every field's domain is declared exactly once

### 6.1 Rule

The valid range of every field is declared once, next to the field definition.  Every validator — in
parsers, setters, journal operations — cites that declaration rather than repeating the numbers.

### 6.2 Evidence — the map-number domain

`atom_t.map_number` is a `uint16_t`, but the atom-map domain in this library is 0–9999.  The
validators once disagreed: the molecule side validated 0–9999 while `QueryContainer.set_map_number`
and `query_seal` validated 0–65535, so a query could seal with a map number no molecule can carry.

The domain is now declared once as `DEF MAP_NUMBER_MAX = 9999` in `_molecule_arena.pxi`, above the
`atom_t` struct, and every validator interpolates from it — cited by symbol, since a line number drifts
out from under the claim (§10.3):

- `_molecule_container.pxi` — `MoleculeContainer.add_atom`, `MoleculeContainer.set_map_number`
- `_query_container.pxi` — `QueryContainer.set_map_number`
- `_query_seal.pxi` — `query_seal`
- `_molecule_arena.pxi` — `structure_from_bytes`
- the readers, added since: `_smiles_read.pxi:smi_bracket`, `_smarts_read.pxi:sma_bracket`,
  `_smirks_patch.pxi:smk_numbering`

Map number is the one domain resolved so far.  Six others — element, hybridization, stored bond order,
stable id, wedge and stereo group — still restate their numbers inline at each validator.  The rule
governs new and modified code; that list is the deferred work.

### 6.3 Corollary — a domain a caller needs must be EXPORTED, or it gets re-derived wrongly

Outside the core a caller cannot cite the declaration: `DEF` constants do not survive into the Python
surface, so an unexported domain is unavailable and gets re-derived from whatever the core does
publish — wrongly, in exactly the cases the domain existed to cover.

Measured: the arena published `H_UNKNOWN = 15` and nothing about the implicit count's *bound*, so the
CTfile reader derived `H_MAX = 15` from the nibble's width and admitted the sentinel as a count on
three write paths — the worst a valence-rule overflow clamp, putting a **computed** count onto the
value meaning "nobody could compute it".  No corpus reaches it: no real record exceeds four hydrogens.

So: **when a domain is needed at a package boundary, export the domain, not only the values it
excludes.**  `H_IMPLICIT_MAX` goes out beside `H_UNKNOWN`, and the test asserts the RELATION
(`H_UNKNOWN == H_IMPLICIT_MAX + 1`) rather than the two numbers, because two independent literals can
drift while a relation is the fact.

Consequences, each paid for:

- **A width is not a bound.**  `H_NIBBLE_MAX = 15` is a fact about the layout and bounds nothing a
  caller may state.  A bound that admits a sentinel destroys the sentinel.  Where a bound is named it
  must say WHICH range it bounds: `H_EXPLICIT_MAX` is spelled separately from `H_NIBBLE_MAX` even
  though the two are equal today, so the saturation moves with the nibble if it reserves a value.

  **Applied forward.**  The S-group segment carries `NO_INDEX = 0xFFFF` as the sentinel for `index`,
  `ext_index` and `parent`, so if those fields are `uint16_t` the greatest REAL value is **0xFFFE** and
  the validator must say `0xFFFE`, not "fits in the field".  Settle it before the struct exists;
  afterwards the wrong bound is already in a validator someone will cite.  Assert the relation once —
  `NO_INDEX == SG_INDEX_MAX + 1` — rather than three pairs of numbers.
- **A reserved value needs a test on the seam, not a comment.**  A correct `ZERO_VALENCE = 15`
  translation sat under a comment that misdescribed it, on a seam with no test — inviting a later
  reader to delete the translation as redundant and emit V3000 claiming a valence of fifteen.  A
  constant whose whole job is to be reserved has its translation pinned by a test in both directions.
- **For a translation table, the test that matters is on the member the corpus does NOT contain, and
  you find it by enumerating what the corpus does contain and subtracting.**  In V2000's charge code
  `ccc`, `4` means neutral doublet radical while `1..3` and `5..6` are `+3..+1` and `-1..-2`.  A corpus
  survey came out `{0: 8333, 1: 1, 2: 4, 3: 38, 5: 60, 6: 1}` — every member exercised except the
  reserved one, which is the general case, since a reserved value is reserved because it is rare.  The
  table also reads as arithmetic, so regenerating it from the pattern produces `4 → 0` and drops the
  radical; changing that row passed 247 tests before the new test existed.  Drive the search from the
  table's domain, not from the data.
- **An asymmetric translation must say it is asymmetric.**  A radical is written as `M  RAD` and never
  as `ccc=4`, because `ccc` holds one fact and an atom that is both charged and a radical has two.  Pin
  the write direction as well as the read direction, or the asymmetry gets "fixed".
- **A measurement that returns zero is not believable until it has been made to return non-zero.**
  Same family as the reserved value nothing reaches: the answer is absent for a reason unrelated to the
  question.  A probe counting stereo units on sentinel-carrying atoms reported 0 while calling a method
  that does not exist, inside a `try/except AttributeError` that turned every atom into "no
  configuration".  Before believing a zero, feed the probe a case that must come back positive; a clean
  number is the cue to check the instrument.  When the instrument is in another tree, **the retraction
  names the branch that holds it**, so the note reads as a handoff rather than a mystery.

  That handoff closed at 0 sentinel atoms over 5,536 molecules, with the probe `N(F)(F)(F)F` returning 1
  first, and two facts came out of it.  **The probe must be UNBRACKETED, or you are testing the wrong
  door**: a bracket STATES its hydrogen count, so `[N](F)(F)(F)F` never reaches the derivation that has
  no answer, and a control wired to the wrong door still comes back positive.  And **the sweep found a
  real defect**, which is the reason to run one: unbracketed atoms with no valence rule had been storing
  0, indistinguishable from a measured 0, and now store `H_UNKNOWN` with a log.  The same rule caught two
  more instruments: a stereo corpus of 4,990 strings containing no `@` at all, and a `/t`-only InChI
  comparison blind to the commonest case there is, because enantiomers of a one-centre molecule share
  `/t1-` and differ only in `/m0` vs `/m1`.
- **A passing suite in another epic's tree is evidence only for the invariants you can name a failing
  line for, and you find out which by breaking your own invariant in both directions.**  Worse-shaped
  than a zero, because a zero invites suspicion and a green does not.  Measured on the segment-table
  truncation:

  | mutation | what it breaks | tests failed | in the other tree |
  |---|---|---|---|
  | `seg_count` always `SEG_TABLE_MAX` | the truncation decision — waste, no data loss | 1 | **0** |
  | `SEG_XY` dropped when coordinates exist | reachability — payload silently lost | 30 | **21** |

  That package's only arena contact is named public constants and the container API, so it cannot witness
  a layout DECISION; it consumes what the table addresses, so it does witness layout REACHABILITY.  **A
  boundary package validates that nothing was lost; only the owning tree validates that nothing was
  wasted.**  A differently shaped mutation there (`segments[SEG_XY].length = 0`) fails the same 21 tests,
  diffed name by name, which makes the 21 a property of the boundary rather than of one edit.
- **Two failure counts are reconciled by diffing names, never by picking the likely one.**  Totals of 30
  against 29 were reconciled by naming the plausible missing test; that test fails under both mutations.
  The line that differs is `test_unpack_rejects_a_coordinate_segment_too_small_for_its_atoms`, failing
  only under truncation — because **a zero-length segment reads as ABSENT, not as TOO SMALL.**  An absent
  segment is defined to read as zero, so the length check `0 < atom_count * sizeof(xy_t)` sits behind a
  presence test and never runs.  `comm -13` on two sorted `grep '^FAILED'` lists answers *which*; a count
  only ever answers *how many*.

### 6.4 Measured — "unset" must not be spelled the same as "stated zero"

`add_atom(implicit_h=None)` stored a zero, defended in the docstring as a third statement.  It is not
one: a caller who omits an argument has said nothing.  It cost `kekule()` the ability to tell a
builder-made aromatic N nobody had counted from one stated to carry no hydrogen, and on a two-nitrogen
five-ring the free choice landed on the wrong nitrogen — imidazole kekulised to a four-valent neutral
N, pyrazole to a two-valent N with no hydrogen.  The default is now `H_UNKNOWN`.

**Review question, for every "unset": ask what a caller who means zero says instead, and if the answer
is "the same thing" the encoding is wrong.**  Same family as §6.3's `H_NIBBLE_MAX`, `NO_INDEX` vs
`0xFFFE`, and `H_UNKNOWN` vs `H_IMPLICIT_MAX`.

Four findings, each needing the measurement rather than the argument:

- **A ruling that says "change the default" is a claim about a consumer; check the consumer reads it.**
  Moving the default did nothing alone: `arom_setup` forced every atom absent from its `stated_h` dict
  to unstated and never read the nibble.  The store is now the default source and `stated_h` overrides
  it, and nothing was lost, because a builder mid-flight holds `H_UNKNOWN` in every unwritten nibble.
- **A default that keeps a sentinel out of a code path does not make that path correct, it makes it
  unmeasured.**  `_inchi.pxi` wrote the raw nibble as `num_iso_H[0] = <char> at_implicit_h(a)`, so an
  unknown count went to libinchi as a literal 15 and alanine came out `C3H37NO2` with a different
  InChIKey — live for MDL records with undeterminable counts, untested because the builder could not
  reach it.  It now writes -1, InChI's own "auto", which the import path already reads as `H_UNKNOWN`.
- **The value carries this and not a flag.**  `at_h_pinned` records "somebody wrote a count" but does
  not survive `copy()`: the copy path passes `implicit_h=at_implicit_h(src)` unconditionally, so every
  copied atom comes out pinned.  A sentinel in the nibble copies as itself.
- **A control that writes a property's empty value tests sensitivity; one that writes a different value
  tests that two values differ, which is a weaker claim in the same shape.**
  `test_stereo_bluebook.py`'s ruling-F102 control drops one property at a time and requires
  `dirty == carrying` exactly.  Its drop for `implicit_h` was `None`, which under the old default
  substituted zero instead of erasing, so seven all-zero records became sensitive to a sabotage that
  could not touch them.  The equality caught it; a `dirty >= carrying` would not have.

Blast radius, measured: 552 `add_atom` call sites, 406 omitting `implicit_h`, **48 failing tests in 13
files, all inside `chython/core/test/`** — the readers all state their counts.  The triage unit is the
atom, not the test: `_stereo.pxi` refuses a unit whose ANCHOR has an unknown count and asks nothing
about the substituents, and fixtures are shared, so the honest edit per fixture is two or three atoms
wide.  Ask "which of these atoms is asked", not "did this test want zero".

---

## 7. Include order: `DEF` binds textually, `cdef extern from *` does not

### 7.1 Rule

A `DEF` compile-time constant is substituted at parse time and is visible only from its point of
definition onward in the textual stream of the translation unit.  In a build made of `.pxi` includes:
**a `DEF` must be included before every file that uses it, and a fragment that moves takes its `DEF`s
with it.**

This is the only hard correctness constraint on include order this codebase has.  `cdef` functions and
`cdef struct` declarations are forward-declared module-wide and may appear in any order relative to
their users — measured against Cython 3.2.8, the same toolchain scope as §7.3.  `cdef extern from *`
blocks are hoisted (§7.3).

### 7.2 Evidence — `DEF W0_*` and `DEF JOURNAL_MIN_CAP`

The twelve `W0_*`/`W1_*` span constants are defined in `_features.pxi` (lines 49–60) and read by
`fill_features` and `fill_edge_words` in the same file and by the query primitive compiler in
`_query_boxes.pxi`, included after it.  When the split moved them out of `_structure.pxi`, every use had
to be checked to fall inside the fragment that became `_features.pxi`; had one fallen in the arena
fragment, the `DEF` block would have had to go to `_molecule_arena.pxi` instead.

`DEF JOURNAL_MIN_CAP` is the sharper case, being a boundary a plan got wrong: the plan started the
container's first fragment at the `OP_*` enum, three lines below `DEF JOURNAL_MIN_CAP = 64`, which
`MoleculeContainer._append` reads.  Taking the range literally would have left the `DEF` in the deleted
file and broken the parse.  The constant sits in `_molecule_container.pxi`'s own `DEF` block.

A `DEF` has no forward declaration to fall back on, so the failure surfaces in the file that *reads* it,
not at the boundary that dropped it.  Reason about this constraint first when moving a fragment.

### 7.3 Verified non-hazard — `cdef extern from *` blocks are hoisted, not emitted in place

Cython collects every `cdef extern from *` literal into the generated file's `/* Early includes */`
preamble, ahead of every function body:

```
1132:  /* Early includes */
1141:      static const unsigned short MDL_ISOTOPE[119] = {
1155:      static const unsigned long long SIG_MASK[4] = {
```

The first function body that reads either is near line 29,000.  Concretely: `sig_mask()` lives in
`_elements.pxi` at include position 2 and reads `SIG_MASK`, declared in `_features.pxi` at position 3.
It compiles.

Scope of the claim, as in §3.2: a **measurement of one toolchain** (Cython 3.2.8), not a language
guarantee; the declared floor is `cython>=3.1` (`pyproject.toml:73`).  If that floor moves, re-check in
one command — the first mention of the symbol in the generated C must be its definition:

```bash
grep -n MDL_ISOTOPE chython/core/_core.c | head -1
```

If that prints `static const unsigned short MDL_ISOTOPE[119] = {`, extern ordering is free.  If it
prints a use, §7.1's rule extends to extern blocks and each must precede its readers.

Until then, keeping an extern block ahead of its readers is a **convention**, followed by `_core.pyx`
because it costs nothing and is what a C reader expects.  It is not a correctness requirement and no
task fails for reordering two extern blocks.  The one real ordering among them: they are emitted in
include order, so an extern block referencing another must follow it.  None does.

A FUNCTION BODY may reference a later fragment's extern block, and one does deliberately:
`fill_features` in `_features.pxi` reads `SPAN_MASK[SPAN_IMPLICIT_H]` and `[SPAN_TOTAL_H]` from
`_query_boxes.pxi` to fill both hydrogen spans for an `H_UNKNOWN` atom, because those two masks are the
layout contract between the writer and the query side and a second spelling would drift.  A `cdef`
*function* defined later still does not resolve backwards.

```python
include "_molecule_arena.pxi"     # DEF CHARGE_MIN/CHARGE_MAX, ISOTOPE_MAX, H_NIBBLE_MAX,
                                  # H_EXPLICIT_MAX, H_IMPLICIT_MAX, H_UNKNOWN, MAP_NUMBER_MAX —
                                  # the §6 field domains, read downstream by _molecule_container,
                                  # _query_container and _query_seal.  H_NIBBLE_MAX is the nibble's
                                  # WIDTH and bounds nothing a caller states: an implicit count
                                  # stops at H_IMPLICIT_MAX because H_UNKNOWN takes the top value,
                                  # and a bound that admits a sentinel destroys it
include "_elements.pxi"           # MDL_ISOTOPE extern block
include "_features.pxi"           # SIG_MASK extern block, DEF W0_*/W1_*, _bit_of, fill_features
include "_rings.pxi"
include "_sssr.pxi"
include "_morgan.pxi"
include "_query_arena.pxi"        # DEF QUERY_MAGIC/QUERY_VERSION; enums, Query, query_alloc
include "_query_boxes.pxi"        # DEF Q_BOX_MAX_ANY, Q_ATOM_MAX_BOXES, METAL_*, SPAN_COUNT
                                  # reads DEF W0_*/W1_*, MDL_ISOTOPE (:225), _bit_of (7 sites)
include "_query_seal.pxi"         # DEF Q_AUTOMORPHISM_MAX_ROWS/NODES, Q_TERM_HASH_WORDS
                                  # reads DEF Q_BOX_MAX_ANY, Q_ATOM_MAX_BOXES, SPAN_COUNT
include "_isomorphism.pxi"
include "_molecule_container.pxi"
include "_molecule_views.pxi"
include "_query_container.pxi"
```

The annotations mark the `DEF`s and extern blocks carrying the ordering argument, not an index of the
tree's 37 `DEF`s.  Some annotated names (`QUERY_MAGIC`, `QUERY_VERSION`, `METAL_*`,
`Q_AUTOMORPHISM_MAX_ROWS`/`NODES`, `Q_TERM_HASH_WORDS`) are read only inside their own file, and eight
`DEF`s are absent entirely (`JOURNAL_MIN_CAP`, `QUERY_JOURNAL_MIN_CAP`, `ZERO_PAGE_SIZE`,
`XXH_P1`–`XXH_P5`) because they too are same-file-only.  For the current list,
`grep -n '^DEF ' chython/core/*.pxi`.

`MDL_ISOTOPE` is declared in `_elements.pxi`'s `cdef extern` literal and read by `_features.pxi`'s
`fill_features` and `_query_boxes.pxi`'s `prim_apply`.  `_bit_of` is a `cdef inline` in `_features.pxi`
with seven call sites in `_query_boxes.pxi`.  Neither
constrains include order; the `DEF` block in the same file does.

### 7.4 Rule — a comment inside a `cdef extern from *` literal reaches the generated C

**Rule:** the triple-quoted string in a `cdef extern from *` block is C source, not Cython, so a
`/* ... */` comment inside it is copied verbatim into the generated `_core.c`.  A `#` comment in Cython
code is not.  The boundary is the string literal, not the block — the `#` lines above
`cdef extern from *` in `_elements.pxi` (lines 65–67, recording where `MDL_ISOTOPE` came from) never
reach the C, while a comment three lines lower would.

Prose is therefore not free in a diff checked by comparing generated C.  Three sites broke that
premise by three mechanisms, and only the first is the rule above:

1. **The C comment heading the `SPAN_WORD`/`SPAN_MASK` extern block** named the old file and had to be
   retargeted to `_features.pxi`; that one word changed a line of generated C.  (The comment later
   moved verbatim into `_query_boxes.pxi`.)
2. **An ordinary `#` comment**, which does not reach the C as a comment — but Cython echoes source
   lines into its `/* "file":line */` position comments, so re-wrapping it changed the C anyway.
3. **`_core.pyx`'s module docstring** — a Python string constant, emitted as
   `static const char __pyx_k_The_core_one_extension_module_o[]`; its layer list had to be rewritten
   for the new file names.

Comments may be edited, but an edit to one is declared as an intentional change rather than found as a
diff.

### 7.5 Why the shared file must be named for its subject

Calling the shared vocabulary file `_molecule_features.pxi` would imply the molecule side owns the
feature encoding.  The seven `_bit_of` call sites in `_query_boxes.pxi` and the `MDL_ISOTOPE` read in
`prim_apply` make the query compiler an equally first-class reader, so a name that picks a side implies
an ownership that is not there.  `_features.pxi` names the subject.

### 7.6 Measured hazard — a misspelled field on a `cdef struct` local is NOT a compile error

Cython converts the struct to a Python dict and does a `setattr` on the temporary.  No error, no
warning, and the write goes nowhere.  Found by renaming `op.a`/`op.b` to `op.n`/`op.m` on a
`cdef qop_t op` whose fields are `op, a, b, opcode, kind, value, negated`:

```c
__pyx_t_4 = __pyx_convert__to_py_struct____pyx_t_..._qop_t(__pyx_v_op);   /* struct -> dict */
if (__Pyx_PyObject_SetAttrStr(__pyx_t_4, ..._n_u_n, __pyx_t_1) < 0) ...   /* dict.n = value */
```

Verified in both directions: with the defect in place a clean rebuild
(`rm -f chython/core/_core.c && python setup.py build_ext --inplace`) leaves the warning table unchanged
(§9.7) and emits the C above; reverting it restores a green suite.  The line right above,
`op.op = QOP_ADD_BOND`, compiles to a direct member access, so the failure is per-attribute, not
per-struct.

Three consequences:

- **A struct-field rename is not compiler-checked in this build.**  The opposite argument — "a missed site
  is a compile error, so renaming `matcher_t *m` across 279 sites is safe" — is withdrawn.  Any
  struct-field rename needs a test that EXECUTES each renamed line.
- **A rarely-taken branch would ship broken.**  This one was caught because 106 tests run through
  `add_bond`; a field write in an error path has no such luck, and the symptom is
  `AttributeError: 'dict' object has no attribute ...` from a function with no dict in it, which reads
  like a caller's bug rather than a typo three lines up.
- **Generic operand fields must keep generic names.**  `qop_t.a/.b` and `journal_t.a/.b` are NOT atom
  endpoints — `journal_t.b` holds a coordinate in one op and an S-group byte in another — so renaming them
  to `n`/`m` would be a lie the compiler cannot catch, in a file where it catches nothing about field
  names at all.

---

## 8. Summary reference

| rule | location in this document |
|---|---|
| One file, one primary definition | §1.1 |
| Target naming table (authority) | §1.3 |
| No query views file (decision) | §1.3 |
| No molecule seal file today (decision) | §1.3 |
| `Structure` is not renamed `Molecule` | §1.3 |
| Bind struct pointer once — threshold is two accesses in one basic block, reads included | §2.1 |
| Never bind a pointer into a `packed` struct; hoist the values | §2.3 |
| After hoisting, re-read the neighbouring block: if it now looks the same, extract it | §2.4 |
| Do not rely on optimiser | §3 |
| AoS vs SoA layout | §4 |
| Scratch: one malloc, one free, one check | §5 |
| Domain declared once | §6 |
| A domain a caller needs must be exported, or it gets re-derived wrongly | §6.3 |
| A width is not a bound; a bound that admits a sentinel destroys it | §6.3 |
| A reserved value's translation is pinned by a test, not by a comment | §6.3 |
| A translation table's test belongs on the member the corpus does NOT contain | §6.3 |
| An asymmetric translation must say so, and pin both directions | §6.3 |
| A zero is not believable until the measurement has returned non-zero | §6.3 |
| A negative control can be wired to the wrong door and still come back positive | §6.3 |
| An attribute name is not a type — count the class, not the spelling | §9.4 |
| A migration cost quoted to its payer gets two independent counts, both by grep | §9.4 |
| A line executing is not an assertion depending on its value — coverage ≠ correctness | §9.6 |
| An execution trace proves your list is live, never that it is complete | §9.8 |
| A completeness claim verified against your own enumeration is not verified | §9.8 |
| Scope a public-attribute sweep to the repository, not to the migrating package | §9.8 |
| A line that runs is not a line that is asserted — break it and name the failing test | §9.8 |
| An equality between two names for one value tests the aliasing, not the value | §9.6 |
| A `Bond`'s endpoint ORDER is load-bearing downstream, not an implementation detail | §9.6 |
| A hazard argued by composing two true statements has not been measured | §9.6 |
| The alarming version of a finding gets less scrutiny than the boring one | §9.6 |
| A shim blocked upstream of its owner must say so, or it reads as laziness | §9.4 |
| S-group segment: arena-native is what an EDIT can invalidate, the rest is opaque | §1.4 |
| S-group references are endpoint pairs, never bond indices | §1.4 |
| A zero-reference S-group record survives; compaction remaps and REPORTS drops | §1.4 |
| Another tree's green counts only for invariants you can name a failing line for | §6.3 |
| Two failure counts are reconciled by diffing names, never by picking the likely one | §6.3 |
| A zero-length segment reads as absent, so no size validator can catch it | §6.3 |
| A `DEF` precedes every reader; a moved fragment takes its `DEF`s | §7.1 |
| `cdef extern from *` blocks are hoisted — ordering them is convention, not correctness | §7.3 |
| A comment inside a `cdef extern from *` literal reaches the generated C | §7.4 |
| Shared file gets a subject name, not a side name | §7.5 |
| A misspelled `cdef struct` field compiles to a dict setattr — renames need executed tests | §7.6 |
| An identifier sweep never touches prose, and `a` is an English word before it is a name | §9 |
| One word for an atom: `n`, `n, m` for a pair; counts get named first | §9 |
| The warning gate covers Cython's warnings too; `implicit declaration of` is a failure | §9.7, §2.3 |
| A remembered warning count goes stale — grep the build, do not recall it | §9.7, §2.3 |

---

## 9. One word for an atom: `n`, and `n, m` for a pair

### 9.1 Rule

An atom's identity is `n`; the two endpoints of a bond are `n, m` — the spelling a caller types, and the
one `self._bonds[n][m]` reads as.  None of the alternatives may come back: `sid`, `sids`, `a`/`b`,
`sid_a`/`sid_b`, `qsid`, `anchor_sid`, `idx_to_sid`, or bare `number` and `index`.  Compound roles take
the same letter with the role in front: `anchor_n`, `parent_n`, `h_n`, `idx_to_n`.  A list of them is
`numbers`.

### 9.2 The letter was not free, and that is the expensive half

`n` and `m` were already in use 895 times in code, almost never for an atom:

| existing meaning | where | sites |
|---|---|---|
| `matcher_t *m` — the matcher receiver | `_isomorphism.pxi` | 279 |
| atom count / bond count / a length | 10 files | ~36 declarations, 73 parameters |
| a running fill cursor | `_query_seal.pxi` | 8 ranges |
| a `MoleculeContainer` | container, views | 9 |

So every file swept gets its counts named first — `n_atoms`, `n_bonds`, `n_persistent`, `n_refs`,
`n_nbrs`, `fill`, `want`, `total` — and the letter is handed to the atom afterwards.  A count called
`n` next to an atom called `n` is worse than either name alone.

### 9.3 Boundary: the letter follows the file's addressing mode

**`n`/`m` name atoms in every file where an atom has a name at all.**  Kernels addressing atoms by
DENSE INDEX only — `i`, `k`, `s`, `t`, `idx` — never name a stable id, and there `m` may keep naming a
struct receiver.  `_isomorphism.pxi` decides this: its 279 `m`s are the matcher struct, which is `self`
and not an atom, and the file contains no atom-valued identifier at all (every `sid` hit in it is the
word "inside" or "side").  `_query_container.pxi` is the opposite case: it names atoms in its public
methods AND ran a local `matcher_t m`, so there the matcher became `mt`.  The test is per FILE: no
single file may spell two things `n` or two things `m`.

### 9.4 Exemptions, each for a stated reason

- **`set_wedge(narrow, wide, ...)` and `wedge_of(narrow, wide)`.**  A wedge is DIRECTIONAL — the narrow
  end is at the stereocentre — so `n, m` would make a silently reversible mistake look symmetric.
  Ruled at dispatch; do not "fix" it.
- **`qop_t.a/.b` and `journal_t.a/.b`.**  Generic operand slots, not endpoints: `journal_t.b` holds a
  coordinate in one op and an S-group byte in another.  See §7.6 — renaming these is a lie the compiler
  cannot catch, and the attempt was caught by 106 failing tests.
- **`bond_a`/`bond_b` in `_query_seal.pxi`.**  Parallel arrays of SLOTS, not of numbers.
- **`Bond.a`/`Bond.b`.**  DISCHARGED — the properties are gone.  It was a read-only migration shim, not
  an exemption, because a public attribute rename cannot land atomically across worktrees.  The one rule
  it leaves behind: **an attribute name is not a type, so count the class, not the spelling.**  `.a`/`.b`
  in `formats/ctfile` is mostly `CtabBond`, whose `self.a`/`self.b` hold FILE INDICES, and `_ctab.py`'s
  `mol.set_wedge(sids[b.a], sids[b.b], ...)` puts the two one subscript apart.  `\.a\b|\.b\b` finds the
  SPELLING and over-includes; `\.bonds\(\)` finds the CLASS and decides, so a cost quoted to the epic
  that must pay it gets two greps rather than a reclassification by eye.  Two agreeing greps still missed
  a site, and the trace over the sites they found could not say so — §9.8.
- **A shim's deletion condition can have a prerequisite its owner does not control, and then the
  condition must say so.**  `bond.n` did not exist in the caller's tree, so the order was fixed: core
  lands, they merge, nine lines in one commit, then the properties go.  "Nobody got round to it" and
  "blocked upstream of the owner" are indistinguishable in a comment stating only the condition.
  DISCHARGED, with no assertion added by the migration commit: `assert bond.a == bond.n` tests the
  aliasing and not the value, since both properties returned the same two slots — the shim's own test made
  that mistake and survived transposing the endpoints at the source.  What replaced it asserts the old
  names no longer RESOLVE, because a name that quietly comes back is how a release ships two spellings
  for one endpoint.
- **`Atom.stable_id`, `MoleculeContainer.stable_ids`, `atom_t.stable_id`.**  Not swept by decision: the
  same concept under a third name, and unifying them breaks another package's imports, so it needed a
  ruling.  DISCHARGED as rename with no shim — `Atom.n`, `atom_t.n`, `mol.number_of(index)`,
  `mol.atom_numbers`, `QueryContainer.query_numbers()`, and `mol.stable_ids` DROPPED rather than renamed.
  Three things carry forward:

  - **The struct member was the cheap half and the public property the expensive one.**  Renaming
    `atom_t.stable_id` touched every `.stable_id` read in 6 files and Cython checked each, so a miss was a
    compile error; the ~700 Python-side lines had only the suite, which §9.6 measured as coverage without
    verification.  So: value tests on `n` and `number_of` first, then the sweep.
  - **Ask whether a name should be renamed or deleted.**  `stable_ids` and `atoms_numbers` had returned
    the same list since the container was written, so a rename would have shipped two spellings under new
    names.  Delete one, rename the other — and the survivor needed fixing too, since `atoms_numbers`
    carries a plural prefix where this core says `atom_count`, hence `atom_numbers`.  What replaced
    both asserts the old four do not RESOLVE; `atoms_numbers` was not added to the chython-2 alias block,
    and the comment there records why.
  - **A shape written in a docstring is code.**  The gate a sweep can state — no `\.stable_id`,
    `\b_stable_ids?\b`, `stable_id_of`, `atoms_numbers` anywhere — says nothing about the ~60 docstrings
    spelling a dict shape `{stable_id: (x, y)}`, which name the field without matching any pattern.
    §9.5 applies: substitute in the code half, read the rest by eye, and make the residual grep the bare
    word rather than the punctuation that makes it an access.

### 9.5 An identifier sweep never touches prose

Three corruptions from substituting a single letter into English:

- `for n in mol` in a docstring became `for n_atoms in mol` — damaging text that was already written in
  the new convention.
- `"""Journal a bond between two known stable ids"""` became `"""Journal n bond ..."""`, three times.
  **`a` is an English article before it is an identifier**, which makes `a` -> `n` the most dangerous
  rename in the set and the one that cannot be automated.
- A comment naming a renamed tuple field (`(parent_sid, isotope_kind, count)`) went stale silently,
  because it is prose that LOOKS like code.

So: substitute in the code half of the line only, then read every changed comment and docstring by eye,
then grep the diff for the new name appearing next to an English verb.  Neither the compiler nor the
test suite sees any of this.

### 9.6 Measured — the sweep's renamed lines were EXECUTED by hundreds of tests and verified by none

A `sys.settrace` run showed all nine `Bond.a`/`Bond.b` migration sites execute, and named the limit:
**a rename typo that swapped `n` and `m` would still execute all nine lines.**  Coverage and correctness
come from different tests, and a trace supplies only the first.  (That was the right limit to name and
the wrong one to worry about: there were **ten** sites, and the trace could not have said so — §9.8.)

Applied to the core sweep by mutation, it found a hole:

| mutation | tests failed | where |
|---|---|---|
| transpose `bd._n`/`bd._m` in `bond(n, m)` | **0 of 2467** | nowhere |
| transpose `bd._n`/`bd._m` in `bonds()` | 1 | `formats/ctfile`, **not the core** |

So a `Bond`'s endpoint orientation was pinned by exactly one test in another package — one about to
migrate off `Bond.a`/`Bond.b`, so the guard would have left with the migration.  Both epics ran the
mutation independently and got the same single failure.  It is now pinned in the core by
`test_bond_answers_the_endpoints_in_the_ORDER_ASKED` and
`test_bonds_yields_each_bond_once_with_the_lower_stable_id_FIRST`, both confirmed to fail under it.

**The explanation of WHY it mattered was wrong, which is the more instructive half.**  The claim was that
a transposition silently inverts a stereo parity in the CTfile wedge writers.  Measured instead: the
writer looks up `wedge_of[(n, m)]`, and **on a miss tries `(m, n)` and swaps** so the narrow end is
written first, which removes orientation from the answer before any parity is computed.  With the
mutation still applied, 264 corpus records and 1191 wedge bonds round-tripped with 0 parity mismatches,
against 243 of 264 under the positive control that swaps emitted stereo codes 1 and 6.  What the one
failing test guards is that the emitted bond block's endpoint order reflects what the molecule holds — a
real fidelity property, but not the silent-corruption one.

- **A hazard argued by composing two true statements is not measured.**  "A parity inversion is invisible
  to a round trip" is true; "this writer swaps endpoints on a wedge miss" is true.  They do not compose,
  because the swap happens *before* the parity is computed.
- **The alarming version of a finding gets less scrutiny than the boring one**, from its author most of
  all — and sharper, **a finding that makes your own work look important is the one to measure first.**
  Inventing a parity hazard made a correct fix rest on a false premise, which is how a later reader
  defends the wrong invariant.  It cuts both ways: the other epic accepted the mechanism unmeasured
  because it made one of its own tests load-bearing for a corruption class.  So the trigger is not "am I
  exaggerating" but "does this conclusion raise the stakes of something I own".

**An equality between two names for one value tests the aliasing, not the value.**  The shim's own test
asserted `(bd.a, bd.b) == (bd.n, bd.m)`, which says nothing about what either holds, so transposing the
source satisfied it.  A shim test needs a second assertion against something outside the shim.

### 9.7 Measured — an undeclared local is a WARNING, and nothing else in the pipeline sees it

Renaming uses without renaming declarations in `_query_seal.pxi` — the component loop variable and the
neighbour counter had both been spelled `n` — left them undeclared.  Cython's response is not an error:

```
warning: chython/core/_query_seal.pxi:958:12: implicit declaration of 'comp'
warning: chython/core/_query_seal.pxi:1270:16: implicit declaration of 'n_nbrs'
```

An implicitly declared local is a **Python object**, and every consumer still behaves: `comp_of[s] !=
comp` coerces, `comp_group_of[comp]` coerces back to a C index.  Which is why it got through everything:

| gate | verdict on the defect |
|---|---|
| Cython | warning only, exit 0 |
| C compiler | nothing — the generated C is valid |
| the 2467-test suite | **green**, before and after the fix |
| a 2000-iteration seal benchmark | **no measurable cost** — 6.4-7.0 us per build+seal either way |

Do not claim a speedup that was not measured: `comp_count` is 3 and `n_nbrs` tops out at 4, so a Python
object in those loops amortises to nothing.  The reason to fix it is **divergence on overflow** — a
`uint32_t` wraps and a Python int does not, so an undeclared counter silently stops being the type the
surrounding arithmetic assumes.  Not reachable at those two sites, and nothing downstream would tell you
when it becomes reachable.

**Consequence for the gate.**  §2.3's warning gate said "exactly N warnings and nothing else", and that
wording was wrong three times:

1. **C-compiler-only**.  Read literally it gated only the C compiler's output, which is how two Cython
   warnings sat unnoticed in a green tree.
2. **Wrong count**.  It said "three" while the truth was four, because the number was remembered rather
   than measured.
3. **`sort -u` destroys multiplicity**.  Normalizing line numbers and then deduplicating collapses every
   `-Wunreachable-code` site into one line, so a new site in a file that already has one would be
   invisible.  **Never use `sort -u` for the gate**; use `sort | uniq -c`, where a new site raises a
   count.

The current gate:

```bash
rm -f chython/core/_core.c && python setup.py build_ext --inplace 2>&1 | grep -i warn \
  | sed -E 's/:[0-9]+:[0-9]+:/:L:C:/g; s/:[0-9]+:/:L:/g' | sort | uniq -c | sort -rn \
  | tee /tmp/warn-table.txt
```

Measured on a macOS universal build (arm64 + x86_64) at `922e2d2`, 2026-09-09 — **15 distinct
normalized messages, 40 total occurrences** (the rows' sum, so re-derive it rather than recall it),
**10 of them `-Wunreachable-code` and 10 `-Wsign-compare`**:

```
10 chython/core/_core.c:L:C: warning: comparison of integers of different signs: 'int' and 'unsigned int' [-Wsign-compare]
10 chython/core/_core.c:L:C: warning: code will never be executed [-Wunreachable-code]
3  warning: chython/core/_pach.pxi:L:C: local variable 'far_prev' might be referenced before assignment
3  warning: chython/core/_pach.pxi:L:C: local variable 'far' might be referenced before assignment
2  warning: chython/core/_pach.pxi:L:C: local variable 'n_chain' might be referenced before assignment
2  warning: chython/core/_pach.pxi:L:C: implicit declaration of 'err'
2  10 warnings generated.
1  warning: chython/core/_smirks_patch.pxi:L:C: implicit declaration of 'exc'
1  warning: chython/core/_smirks_patch.pxi:L:C: Unused entry 'other'
1  warning: chython/core/_smiles_write.pxi:L:C: implicit declaration of 'shuffle'
1  warning: chython/core/_smiles_write.pxi:L:C: Unused entry 'other'
1  warning: chython/core/_query_container.pxi:L:C: Unused entry 'n'
1  warning: chython/core/_pach.pxi:L:C: Unused entry 'k'
1  warning: chython/core/_molecule_container.pxi:L:C: implicit declaration of 'err'
1  warning: chython/core/_molecule_container.pxi:L:C: implicit declaration of '_'
```

One `build_inchi.py` `UserWarning` appears in a checkout where the `INCHI` submodule is absent, and it
matches the gate on two lines — its own message, and the `warn(...)` source line Python echoes under it.
The call is `build_inchi.py:82`'s; `chython/core/*.pxi` and `_core.pyx` contain no `warn(` at all, so no
row of the table can be a `warnings.warn` from the core.  Both lines are environmental and are not rows
of the table.

Clang rows scale with the number of architectures: ten `-Wunreachable-code` occurrences on a universal
build are five sites × two arches, and a single-arch build sees five; the `2  10 warnings generated.`
row halves the same way.  Cython rows are arch-**independent**, and their counts are source-line
multiplicities: `2  implicit declaration of 'err'` counts two distinct source lines in `_pach.pxi`.

**A clang row is a count of occurrences, not of sites, so it must be attributed before it is frozen.**
Map each `_core.c` line back through the generated `/* "chython/core/….pxi":NNN */` comments — the gate's
own `sed` erases the line numbers, so tee the raw build log and read it for this. As attributed:

| row | occurrences | sites | where |
|---|---|---|---|
| `-Wsign-compare` | 10 | 5 comparisons at 3 source lines | `_molecule_arena.pxi`'s segment-table guards — `_structure_resolve_persistent` (`:660`), `structure_append` (`:2679`), `structure_retire` (`:2711`) — where an `int present`/`int seg` meets the `SEG_*` enum; a `seg < A or seg >= B` guard is two comparisons |
| `-Wunreachable-code` | 10 | 5 | `_rings.pxi:_seen_lookup` and `_isomorphism.pxi:matcher_next`, both `noexcept nogil` bodies that always return, so Cython's exit-code `__pyx_r = 0` is dead; three instantiations of Cython's `CIntToPyUnicode` utility (`uint8_t`, `uint16_t`, `unsigned char`), which is generated code with no `.pxi` origin |

**Which rows arrived when.**  The table above is the whole build's state, and a table that absorbs a
drift it did not measure cannot enforce "no new message and no count that rises".  `922e2d2` sits above a
merge, so both of its parents were measured with the same gate — `3d2fa85` for this branch and `3430963`
for master:

| row | `3d2fa85` (branch) | `3430963` (master) | merged |
|---|---|---|---|
| `-Wunreachable-code` | 10 (5 sites) | 8 (4 sites) | 10 (5 sites) |
| `N warnings generated.` | `10` | `9` | `10` |
| `_smirks_read.pxi:L:C: noexcept clause is ignored for function returning Python object` | 1 | — | — |
| `_smirks_patch.pxi:L:C: Unused entry 'other'` | — | 1 | 1 |
| every other row, `-Wsign-compare` included | identical | identical | identical |

So the only part of the table this branch owns is the fifth `-Wunreachable-code` site, the `unsigned
char` instantiation of `CIntToPyUnicode`.  `_smirks_patch.pxi`'s `Unused entry 'other'` is a **second**
row and not `_smiles_write.pxi`'s relocated — both files carry one, at `_smirks_patch.pxi:465` and
`_smiles_write.pxi:3124`.  Every Cython row in the merged build is present in master's own build, which
is what allows an `implicit declaration of` row to be attributed rather than fixed: in code this branch
writes it is a hard failure and gets a fix.

**To attribute a row, build the other side too.**  `git archive <rev> | tar -x -C /tmp/<dir>` and run the
gate there — a second checkout is not needed and a shared worktree must not be moved.  Nothing else
distinguishes a row you introduced from one you inherited.

The invariant is **"no new message and no count that rises"**, which is arch-independent.  The gate
covers **every** line matching `warn`, Cython's included, and `implicit declaration of` is a hard
failure, not a note.  Run the gate before the tests after any identifier sweep; the tests will not tell
you.

### 9.8 Measured — the trace proved the LIST was live and said nothing about completeness

§9.6 records the `sys.settrace` check that all nine `Bond.a`/`Bond.b` migration sites execute, and the
limit it named at the time.  **There were ten sites.**  The tenth was
`core/test/test_smiles_read.py:53`, a SMILES-epic helper inside the core's own test directory; both
epics' counts came from a grep scoped to `formats/ctfile`, the package doing the migration.  Two
independent greps agreeing is not either of them covering the tree.

**Scope a public-attribute sweep to the repository, not to the package doing the migration.**  An
attribute on an exported type has no owner; anything that can import the type is a call site.

> **An execution trace proves the sites in your list are live.  It never proves the list is complete —
> and it reads as more rigorous than the grep it rests on.**

The trace was pointed at nine hard-coded `(file, lineno)` pairs, so it could only confirm the enumeration
it was handed; "all nine execute, so the suite is a real gate on this diff" attached a measurement's
authority to an unmeasured premise.  **A completeness claim verified against your own enumeration is not
verified.**

Corollary: **a line that runs is not a line that is asserted.**  §9.6's mutation table is the evidence:
every renamed site executed, and transposing the endpoints under one of them failed nothing.  Closing
that takes the step a trace does not supply — break each line deliberately and name the test that fails,
in both failure modes: the reference vanishing (logged) and the reference resolving to the wrong bond
(silent).  Chasing the same probe into `_v3000.py`'s S-group bond lookups produced a coverage claim that
§10.3 measures and refutes; the lookups execute, and it is their `None` arms that do not.

---

## 10. A relayed claim about this tree is stale until re-measured

### 10.1 Rule

A coordinator's or peer epic's report names a line, a defect or a state in **their** view of the tree, and
between the observation and the instruction the tree moved.  Before acting on a relayed claim, locate it
by **content** and check whether what it describes is still true here.  Two instances, both cheap to check
and both wrong:

- **A line number.**  "`_molecule_container.pxi:1720-1723` documents the zero as deliberate."  That range
  is `set_hydrogens`' `ValueError`; the comment was at 1613.  Grepping the quoted phrase found it in one
  step, and reading 1720 instead would have edited the wrong guard.
- **A defect.**  "`_inchi.pxi`'s `translate_stereo` calls for `SU_CIS_TRANS`/`SU_ALLENE` raise
  `ValueError: order is not a permutation of the unit refs`."  True when the writer epic saw it, and
  **fixed at `eef2732` before the instruction arrived** — the parity goes through
  `_ich_bond_order_from_refs` in the unit's own refs frame, with 17 tests on both kinds in both
  directions, and the layer tracks the parity (`/b4-3+` ↔ `/b4-3-`, `/t1-/m0/s1` ↔ `/t1-/m1/s1`, no layer
  at parity 0), which cannot happen if the record is skipped.

**A report's shelf life is the interval in which nobody touched the file, and the reporter does not know
that interval.**  Confirm-then-act costs one grep.

### 10.2 Corollary — a stale report is still worth reading for what it was pointing at

Checking the discharged `translate_stereo` report surfaced an untested asymmetry in the same two branches.
Both InChI records need `neighbor[0]` (X) bonded to `neighbor[1]` (A).  The **allene** anchors on the
centre, so `refs[0:2]` is whichever terminal perception filled first and the export must orient A/B by
adjacency to X — pinned.  The **cis/trans** export writes `neighbor[0] = refs[0]`, `neighbor[1] = anchor`
with no adjacency check, and is right to, because `_stereo.pxi` pass 2 anchors the unit on the
lower-indexed terminal and fills `refs[0:2]` from that same terminal.  That guarantee lives in another
file, was asserted nowhere, and a green round trip does not imply it;
`test_a_cis_trans_units_refs_0_is_bonded_to_its_anchor_and_refs_2_is_not` now states it.

**Where two sibling branches handle the same requirement differently, the one that does nothing is the
one carrying an unstated assumption, and it is the one to assert.**

### 10.3 Measured — a drifted line number inverts the claim

The MDL merge was reported to land with a coverage hole: two lines in
`formats/ctfile/_v3000.py:_parse_sgroups` "never execute under 230 passing tests, because every V3000
S-group test uses a record type whose bond list is dropped on input."  The report cited line numbers that
were already several lines off, which is why the statements are named here instead.  Traced with
`sys.settrace` over the 235-test ctfile suite:

| statement | claim | when traced | now |
|---|---|---|---|
| `number = bond_position.get(pair)` (XBONDS) | never executes | **executes** | executes |
| `number = bond_position.get(pair)` (CSTATE) | never executes | **executes** | executes |
| its `None` arm, `log.append('... reference dropped')` | — | **never executes** | executes |
| its `None` arm, `log.append('... state dropped')` | — | **never executes** | executes |

`test_v3000.py` reads an `SRU` with `XBONDS=(1 1)` and a `CSTATE=(4 1 ...)` and asserts both come back, so
the lookups are covered and the premise about record types is wrong.  What was uncovered is the **`None`
branch of each lookup** — the two diagnostics for "a stored endpoint pair no longer names a bond".  A
handful of lines inverts the finding: *the feature is untested* versus *the feature is tested and its
failure diagnostic is not*.  **Cite a statement, not a line number**, or the next edit above it makes the
claim unfalsifiable.

The `now` column is the point.  Those branches are reachable only once a molecule carrying a
bond-referencing S-group can survive a `delete_bond`, which the accessor made possible, and
`test_v3000.py:test_a_bond_reference_whose_bond_was_deleted_is_dropped_and_reported` closes them.  **A
coverage hole in another epic's file can be a statement about a capability yours does not
have yet, and then closing it is your work and not theirs.**
