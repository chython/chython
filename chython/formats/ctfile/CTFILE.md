# CTFile Formats Specification Summary

Extracted from MDL CTFile Formats (October 2003) + BIOVIA 2020 additions. Covers molecular structure
and stereo information relevant to chython.

Sections 1–10 keep the reference's numbering, so an older citation still lands in the right place.
§11 is what this package stores for the constructs the reference describes loosely; §12 is the
dialect facts measured against other implementations.

---

## 1. MOL V2000

A molfile = Header Block + Connection Table (CTAB).

### 1.1 Header Block (3 lines)

```
Line 1: Molecule name (max 80 chars, unformatted)
Line 2: IIPPPPPPPPMMDDYYHHmmddSSssssssssssEEEEEEEEEEEERRRRRR
         (initials, program, date, dimensions, scaling, energy, registry)
Line 3: Comment (blank if none)
```

Lines 2-3 may be blank. Line 1 must NOT start with `$MDL`, `$$$$`, `$RXN`, or `$RDFILE`.

### 1.2 Counts Line

```
aaabbblllfffcccsssxxxrrrpppiiimmmvvvvvv
```

| Field | Width | Meaning                                                   |
|-------|-------|-----------------------------------------------------------|
| aaa   | 3     | Number of atoms                                           |
| bbb   | 3     | Number of bonds                                           |
| lll   | 3     | Number of atom lists (query)                              |
| fff   | 3     | Obsolete                                                  |
| ccc   | 3     | Chiral flag: 0=not chiral, 1=chiral                       |
| sss   | 3     | Number of stext entries                                   |
| xxx   | 3     | Obsolete                                                  |
| rrr   | 3     | Obsolete                                                  |
| ppp   | 3     | Obsolete                                                  |
| iii   | 3     | Obsolete                                                  |
| mmm   | 3     | Number of properties lines (always 999, read until M END) |
| vvvvvv| 6     | Version: ` V2000` or ` V3000`                             |

### 1.3 Atom Block

One line per atom, fixed-width columns:

```
xxxxx.xxxxyyyyy.yyyyzzzzz.zzzz aaaddcccssshhhbbbvvvHHHrrriiimmmnnneee
```

| Field | Cols | Meaning                                 | Values |
|-------|------|-----------------------------------------|--------|
| x,y,z | 10.4 each | Coordinates (Angstroms)                 | float |
| (space) | 1 | separator                               | |
| aaa | 3 | Atom symbol                             | Element, `L`, `A`, `Q`, `*`, `LP`, `R#` |
| dd | 2 | Mass difference (deprecated, use M ISO) | -3..+4 |
| ccc | 3 | Charge (deprecated, use M CHG)          | 0=none, 1=+3, 2=+2, 3=+1, 4=doublet, 5=-1, 6=-2, 7=-3 |
| sss | 3 | Atom stereo parity                      | 0=none, 1=odd, 2=even, 3=either |
| hhh | 3 | Hydrogen count+1 (query)                | |
| bbb | 3 | Stereo care box (query)                 | |
| vvv | 3 | Valence                                 | 0=default, 1-14, 15=zero valence |
| HHH | 3 | H0 designator                           | |
| rrr | 3 | Not used                                | |
| iii | 3 | Not used                                | |
| mmm | 3 | Atom-atom mapping number                | 0=no mapping, >0=mapped |
| nnn | 3 | Inversion/retention (reaction)          | 0=none, 1=inverts, 2=retained |
| eee | 3 | Exact change (reaction)                 | 0=none, 1=exact |

Fields `dd` and `ccc` are superseded by `M ISO`, `M CHG`, `M RAD` in the properties block.

### 1.4 Bond Block

One line per bond, fixed-width columns:

```
111222tttsssxxxrrrccc
```

| Field | Width | Meaning                         | Values |
|-------|-------|---------------------------------|--------|
| 111 | 3 | First atom number               | 1..natoms |
| 222 | 3 | Second atom number              | 1..natoms |
| ttt | 3 | Bond type                       | 1=single, 2=double, 3=triple, 4=aromatic, 5-8=query |
| sss | 3 | Bond stereo                     | **Single:** 0=none, 1=Up(wedge), 4=Either, 6=Down(hash). **Double:** 0=use coords, 3=either |
| xxx | 3 | Not used                        | |
| rrr | 3 | Bond topology (query)           | 0=either, 1=ring, 2=chain |
| ccc | 3 | Reacting center status          | 0=unmarked, 1=center, -1=not center, 2/4/8/12=changes |

**Stereo convention**: The wedge (pointed) end is at the FIRST atom (field 111).

### 1.5 Properties Block

Terminated by `M  END`. The per-atom properties:

```
M  CHGnn8 aaa vvv ...       Charge: vvv = -15..+15
M  RADnn8 aaa vvv ...       Radical: 0=none, 1=singlet, 2=doublet, 3=triplet
M  ISOnn8 aaa vvv ...       Isotope: absolute atomic mass (positive integer)
M  END                      End of CTAB
```

When `M CHG`/`M RAD` present, they supersede ALL atom block charge/radical values (forces 0 on unlisted atoms).

Two non-`M` lines belong to this block as well:

```
A  aaa                      Atom alias: the NEXT line is free label text for atom aaa
V  aaa text                 Atom value: free text attached to atom aaa
```

The line after `A  aaa` is data, not a property line, and may itself begin with `M  `, so it has to be
consumed as text before the property dispatch sees it. See §11.4 for how both are stored.

#### 1.5.1 The V2000 S-group facility

`M  STY` and its companions are the **general** V2000 S-group facility, not a stereo carrier: a group
is declared by `M  STY` and described by further lines that share only its S-group number. They may
arrive in any order, including data before declaration, so a reader creates the record on first
mention by whichever line mentions it.

```
M  STYnn8 sss ttt ...       Declare S-group sss of type ttt (SUP MUL SRU MON COP DAT GEN ...)
M  SSTnn8 sss ttt ...       Subtype
M  SLBnn8 sss vvv ...       External (display) number
M  SPLnn8 sss ppp ...       Parent S-group number
M  SAL sss nn8 aaa ...      Atoms of the group
M  SPA sss nn8 aaa ...      Parent-atom subset (MUL)
M  SBL sss nn8 bbb ...      Bonds of the group, as bond numbers
M  SBV sss bbb x y          One bond plus a display vector along it  (§11.2)
M  SMT sss text             Subscript / label
M  SDT sss name ...         DAT field definition: field name, then FIELDINFO/TYPE/QUERY columns
M  SDD sss <display>        DAT display position and styling  (§11.1)
M  SED sss data             DAT field data, terminating line
M  SCD sss data             DAT field data, continued
```

`SDS`, `SCN`, `SAP`, `SCL`, `SNC`, `SPS`, `CRS`, `MRV`, `LOG`, `APO` and the rest are read and
preserved without interpretation (§11).

#### 1.5.2 Enhanced stereo in V2000 (BIOVIA 2020 extension)

V2000 has no enhanced-stereo block. The 2020 convention states one as a **`DAT` S-group** whose
`SDT` field name is a reserved token, which is a *use* of §1.5.1 and not a facility of its own:

```
M  STY  1   1 DAT           Define S-group 1 as data type
M  SAL   1  n  a1 a2 ...    Atoms in the stereo group
M  SDT   1  MDLV30/STERAC1  Field name identifies the stereo type
M  SED   1                  (empty data)
```

Field names: `MDLV30/STEABS`, `MDLV30/STERAC{n}`, `MDLV30/STEREL{n}` — same semantics as the V3000
collection block (§2.6). Absent such a group, V2000 states relative stereochemistry only through
the chiral flag on the counts line, and readers differ in what they infer from it (§12).

---

## 2. MOL V3000

V3000 file = V2000 "no structure" header (with version stamp `V3000`) + extended CTAB blocks.

### 2.1 General Syntax

- Every line begins with `M  V30 ` (2 spaces after M, 1 after 30)
- Line continuation: `-` as last char, next line's `M  V30 ` prefix stripped and concatenated
- Max 80 chars per physical line
- Values: positional first, then `KEYWORD=value` optional
- List values: `KEYWORD=(N val1 val2 ... valN)` where N = count
- Strings with spaces/parens/quotes must be double-quoted; literal `"` doubled

### 2.2 Overall Structure

```
{V2000 header: name, program line, comment, counts "0  0  0     0  0            999 V3000"}
M  V30 BEGIN CTAB
M  V30 COUNTS na nb nsg n3d chiral
M  V30 BEGIN ATOM
...atoms...
M  V30 END ATOM
M  V30 BEGIN BOND
...bonds...
M  V30 END BOND
[M  V30 BEGIN SGROUP ... M  V30 END SGROUP]
[M  V30 BEGIN COLLECTION ... M  V30 END COLLECTION]
M  V30 END CTAB
M  END
```

### 2.3 Counts Line

```
M  V30 COUNTS na nb nsg n3d chiral [REGNO=regno]
```

| Field | Meaning |
|-------|---------|
| na | Number of atoms |
| nb | Number of bonds |
| nsg | Number of Sgroups |
| n3d | Number of 3D constraints |
| chiral | 1=chiral, 0=not |

### 2.4 Atom Block

```
M  V30 BEGIN ATOM
M  V30 index type x y z aamap [CHG=val] [RAD=val] [CFG=val] [MASS=val] [VAL=val] ...
M  V30 END ATOM
```

| Field | Meaning | Values |
|-------|---------|--------|
| index | Atom index (unique integer >0) | |
| type | Atom symbol | Element string, `R#`, `A`, `Q`, `*`, or `[NOT] [list]` |
| x, y, z | Coordinates | float (Angstroms) |
| aamap | Atom-atom mapping | 0=none, >0=mapped |
| CHG | Charge | integer (-15..+15) |
| RAD | Radical | 0=none, 1=singlet, 2=doublet, 3=triplet |
| CFG | Stereo configuration (parity) | 0=none, 1=odd, 2=even, 3=either |
| MASS | Isotope (absolute mass) | positive integer |
| VAL | Valence | >0 or -1=zero |

Query keywords: HCOUNT, STBOX, SUBST, UNSAT, RBCNT, ATTCHPT, RGROUPS, ATTCHORD.
Reaction keywords: INVRET (0/1/2), EXACHG (0/1).

### 2.5 Bond Block

```
M  V30 BEGIN BOND
M  V30 index type atom1 atom2 [CFG=val] [TOPO=val] [RXCTR=val] [STBOX=val]
M  V30 END BOND
```

| Field | Meaning | Values |
|-------|---------|--------|
| index | Bond index (unique integer >0) | |
| type | Bond type | 1=single, 2=double, 3=triple, 4=aromatic, 5-8=query |
| atom1 | First atom index | |
| atom2 | Second atom index | |
| CFG | Bond stereo | 0=none, **1=Up(wedge)**, 2=Either, **3=Down(hash)** |
| TOPO | Topology (query) | 0=default, 1=ring, 2=chain |
| RXCTR | Reacting center (reaction) | same as V2000 |
| STBOX | Stereo care box (query) | |

V2000 and V3000 disagree on the Down and Either values; see the mapping table in §7.4.

### 2.6 Collection Block (Enhanced Stereo)

```
M  V30 BEGIN COLLECTION
M  V30 MDLV30/STEABS ATOMS=(n a1 a2 ...)
M  V30 MDLV30/STERAC1 ATOMS=(n a1 a2 ...)
M  V30 MDLV30/STEREL1 ATOMS=(n a1 a2 ...)
M  V30 END COLLECTION
```

| Collection name | Meaning |
|-----------------|---------|
| `MDLV30/STEABS` | **Absolute** (ABS): configuration is exactly as drawn, a single known enantiomer |
| `MDLV30/STERACn` | **Racemic** (AND group n): relative config known, mixture of both enantiomers present |
| `MDLV30/STERELn` | **Relative** (OR group n): relative config known, one enantiomer present, which one unknown |

- `n` is an integer >= 1 identifying the group
- Multiple, independently flipping groups can coexist (e.g. STERAC1, STERAC2, STEREL1)
- Within one group all atoms flip together — their relative configuration is fixed
- Atoms not listed in any collection default to ABS
- ATOMS list format: `(count atom_index atom_index ...)`

### 2.7 Sgroup Block

```
M  V30 BEGIN SGROUP
M  V30 index type extindex [ATOMS=(...)] [FIELDNAME=name] [FIELDDATA=data] ...
M  V30 END SGROUP
```

Types: SUP(eratom), MUL(tiple), SRU, MON(omer), COP(olymer), DAT(a), GEN and the rest.

This package models `DAT` fully: `FIELDNAME`, `FIELDDATA` and `FIELDDISP` are parsed and re-emitted.
**Every other type is preserved rather than interpreted** — references (`ATOMS`, `PATOMS`, `CBONDS`,
`XBONDS`, `CSTATE`, `PARENT`, the external index) are translated into the owner's own alphabet, and
every keyword the release does not model rides in the record's `fields` as `keyword -> [values]`.
Nothing is dropped for being unrecognised.

`fields` is keyed by keyword, so keyword **order** across the line is not preserved: the writer emits
a canonical order. A keyword may legitimately repeat on one group (`CSTATE`, `BRKXYZ`, `FIELDDATA`),
which is why each value list is a list.

`FIELDNAME="MRV_IMPLICIT_H"` with `FIELDDATA=IMPL_Hn` is one ordinary `DAT` group among others.

---

## 3. SDF (Structure-Data File)

Multiple molecules + associated data. Format:

```
[Molfile]            <- V2000 or V3000 MOL block (header + CTAB + M END)
>  <FIELD_NAME>      <- Data header (field name in angle brackets)
data value           <- One or more lines of data
                     <- Blank line terminates data item
>  <ANOTHER_FIELD>   <- Repeat for each data field
more data

$$$$                 <- Record delimiter (separates molecules)
```

### Data header format:
```
> [registry_info] <field_name>
```
The `>` must be in column 1. Field name is in `<angle_brackets>`. A value may span several lines; the
blank line, not the line count, terminates it.

### V2000 vs V3000 detection:
The version stamp in the counts line (chars 34-39) is `V2000` or `V3000`. Pragmatically, the line
after the counts line starts with `M  V30 BEGIN CTAB` for V3000. An SD or RD file may mix the two
versions record by record, so the version is sniffed per CTAB and never per file.

---

## 4. RXN (Reaction File) V2000

```
$RXN
reaction name
      IIIIIIPPPPPPPPPMMDDYYYYHHmmRRRRRRR    (program info line)
comment line
rrrppp                                       (counts: 3-digit reactants, 3-digit products)
$MOL
[Molfile for reactant 1]
$MOL
[Molfile for reactant 2]
...
$MOL
[Molfile for product 1]
...
```

- Line 1: `$RXN` identifier
- Line 2: Reaction name (or blank)
- Line 3: Program info (or blank)
- Line 4: Comment (or blank)
- Line 5: `rrrppp` - number of reactants (3 chars) + number of products (3 chars)
  - A third 3-char field carries the agent count: `rrrpppaaa`. The 2003 reference states only two
    count fields; the third is what the ecosystem writes and reads (§12).
- Each `$MOL` delimiter followed by a complete Molfile (header + CTAB + M END)
- Order: all reactants first, then all products, then agents (if any)

---

## 5. RXN V3000 (Extended Reaction File)

```
$RXN V3000
reaction name
program info
comment
M  V30 COUNTS nreactants nproducts [nagents]
M  V30 BEGIN REACTANT
M  V30 BEGIN CTAB
...ctab for reactant 1...
M  V30 END CTAB
M  V30 BEGIN CTAB
...ctab for reactant 2...
M  V30 END CTAB
M  V30 END REACTANT
M  V30 BEGIN PRODUCT
M  V30 BEGIN CTAB
...ctab for product 1...
M  V30 END CTAB
M  V30 END PRODUCT
[M  V30 BEGIN AGENT
M  V30 BEGIN CTAB...END CTAB
M  V30 END AGENT]
M  END
```

- Line 1: `$RXN V3000` (the `V3000` token distinguishes from V2000)
- Lines 2-4: name, program info, comment (same as V2000 rxn header)
- Counts line: `M  V30 COUNTS nreactants nproducts [nagents]`
- Each molecule is a full CTAB block (same format as V3000 MOL, without the outer header)
- No `$MOL` delimiters — molecules are wrapped in `BEGIN/END CTAB` pairs
- No per-molecule headers (no name/program/comment per molecule)

---

## 6. RDF (Reaction-Data File)

Contains molecules OR reactions with associated data. More general than SDF.

```
$RDFILE 1                        <- Header (required, file start)
$DATM    MM/DD/YY HH:mm         <- Date stamp (treated as comment)
$RFMT [$RIREG regno]            <- Reaction record start (or $MFMT for molecule)
$RXN                             <- Embedded rxnfile
...rxnfile content...
M  END
$DTYPE field_name                <- Data field identifier
$DATUM data_value                <- Data value
$DTYPE another_field
$DATUM another_value
$RFMT                            <- Next record starts here
...
```

### Record identifiers:
- `$MFMT` — molecule record (followed by embedded molfile)
- `$RFMT` — reaction record (followed by embedded rxnfile starting with `$RXN`)
- `$MIREG` / `$RIREG` — internal registry reference
- `$MEREG` / `$REREG` — external registry reference

### Data format:
- `$DTYPE field_name` — field name (one per data item)
- `$DATUM value` — data value (can be multi-line for fields >80 chars)

### Key differences from SDF:
- No `$$$$` delimiter (records separated by next `$RFMT`/`$MFMT`)
- Can contain both molecules and reactions in same file
- Uses `$DTYPE`/`$DATUM` instead of `>  <field>`
- Has file-level header (`$RDFILE`, `$DATM`)
- No blank lines allowed except within embedded mol/rxn blocks

---

## 7. Stereo Conventions

### 7.1 Bond Stereo (Wedge Notation)

For **single bonds** at tetrahedral centers, the pointed end of the wedge is at the first atom in the
bond definition; Up (wedge) puts the second atom above the plane, Down (hash/dash) below it.

For **double bonds** (cis/trans), stereo comes from the 2D coordinates of the substituents: value
0 = use coordinates, value 3 = either (ignore stereo).

### 7.2 Atom Stereo Parity (V3000 CFG on atoms)

Calculated by viewing the center from behind the highest-numbered neighbor:
- **1 = odd parity**: atoms 1,2,3 in clockwise order
- **2 = even parity**: atoms 1,2,3 in counterclockwise order
- **3 = either**: unmarked or racemic at that center

In V2000 the atom-block `sss` field is redundant with the wedges: stereo is determined from bond
wedges plus coordinates. In V3000, atom `CFG` is informational; bond `CFG` (wedge/hash) is what
defines stereo.

### 7.3 Enhanced/Extended Stereo

Groups of stereocenters with a stated epistemic relationship: ABS, AND (`STERAC`) and OR (`STEREL`),
per the table in §2.6. Typical sources: a single known enantiomer for ABS, a racemate for AND, a
natural-product isolate of unassigned absolute configuration for OR.

### 7.4 V2000 vs V3000 Bond Stereo Value Mapping

| Meaning | V2000 bond `sss` | V3000 bond `CFG` |
|---------|-------------------|-------------------|
| None | 0 | 0 |
| Up (wedge) | 1 | 1 |
| Either | 4 | 2 |
| Down (hash) | 6 | 3 |

---

## 8. Special Atoms

| Symbol | Meaning | Handling |
|--------|---------|----------|
| `*` | Star atom (attachment point) | Track for ENDPTS bonds, skip as real atom |
| `D` | Deuterium | H with MASS=2 |
| `T` | Tritium | H with MASS=3 |
| `A` | Any atom (query) | Query construct: not a structure atom |
| `Q` | Any non-C non-H (query) | Query construct: not a structure atom |
| `L` | Atom list (query) | Query construct, see §12 |
| `R#` | R-group label | Query construct: not a structure atom |
| `LP` | Lone pair | Not a structure atom |

The reference lists exactly these tokens for the symbol field. It says nothing about free text there,
and implementations differ on what a free-text label in that column means (§12). chython reads such a
text as the display label it is: the atom is kept, the text becomes its alias (§11.4), the element is
carbon as a placeholder and the hydrogen count is `H_UNKNOWN`. Same in V3000's atom-type field, which
is the only channel a label has there, V3000 having no `A` line.

---

## 9. Bond Types

| Value | Meaning | Notes |
|-------|---------|-------|
| 1 | Single | |
| 2 | Double | |
| 3 | Triple | |
| 4 | Aromatic | Query in the V2000 reference, and written by structure writers in practice |
| 5 | Single or Double | Query only |
| 6 | Single or Aromatic | Query only |
| 7 | Double or Aromatic | Query only |
| 8 | Any | Query only / coordinate bond |
| 9 | Coordinate | BIOVIA extension |
| 10 | Hydrogen bond | BIOVIA extension |

---

## 10. Star Atom / ENDPTS Handling (V3000)

Star atoms (`*`) represent attachment points with multiple possible endpoints:

```
M  V30 5 1 1 7            <- bond from atom 1 to star atom 7
                              with ENDPTS on the bond or star atom
M  V30 5 1 1 7 ENDPTS=(3 2 3 4)   <- star atom 7 connects to atoms 2, 3, or 4
```

Format: `ENDPTS=(N atom1 atom2 ... atomN)` where N is count of endpoint atoms. In practice: create a
bond (type 8/special) from the real atom to each endpoint atom.

---

## 11. What this package stores

The reference describes these four constructs loosely or not at all. This is the layout the readers
in `_v2000.py` and `_v3000.py` parse and the writers re-emit, and the shape `_sgroup.py` holds it in.

### 11.1 `M  SDD` and V3000 `FIELDDISP`

One value in two spellings, byte-for-byte identical after the S-group number, so one model serves
both. The layout is fixed, not free text:

```
xxxxx.xxxxyyyyy.yyyy eeefgh i jjjkkk ll m noo
```

```
M  SDD   1     0.7661   -0.8250    DR    ALL  0       0
        num  x (F10.4)  y (F10.4)  <----- display styling ----->
```

`parse_fielddisp` splits it into `(x, y, rest)`; `format_fielddisp` renders it back. **The two F10.4
coordinates are followed by four blanks before the display flags**, which are the first four
characters of the round-tripped tail. They are display coordinates in the molecule's own frame, so
they are parsed out and re-emitted from the moved molecule; the styling tail rides verbatim.

F10.4 is asymmetric — the sign costs a column, so `99999.9999` and `-9999.9999` are ten characters
while `-10000.0000` is eleven. Hence `DISP_MAX = 99999.9999` and `DISP_MIN = -9999.9999`. A wider
value on read means the file was not written to the fixed layout, so the whole value is kept as
opaque text in `fields['FIELDDISP']`; on write a wider value is clamped and reported.

### 11.2 `M  SBV` and V3000 `CSTATE`

Also one value in two spellings: an S-group bond plus a display vector along it. Stored as
`[((a, b), tail), ...]`, structured and not verbatim: the first value in the file is a **bond
number**, a position in the bond block, so it becomes an endpoint pair on read and a fresh number on
write.

A reference no bond answers is kept as `(None, tail)` with the original text in the tail, and
reported — demoted rather than removed, the run length deriving from the number of `cstates`.

### 11.3 `PARENT` and the external index

Both are `uint16` arena slots, with `NO_INDEX = 0xFFFF` as the sentinel and `INDEX_MAX = 65534` as
the largest value distinguishable from it. `0xFFFF` and not `0`, because V2000 `M  STY` writes the
number in a 3-char field where `  0` is representable and nothing in CTfile forbids it.

A V3000 `Sgroup 65535` is legal — the keyword takes an unbounded integer — so a number outside
`0..65534` is **renumbered** on read and reported (`normalize_indices`), rather than demoted to the
sentinel, because two records must stay distinguishable.

A `PARENT` naming a number no record in the file carries is dropped and logged, never raised.

### 11.4 Atom aliases

V2000 `A  <n>` with the label on the following line, V2000 `V  <n> text`, MRV `mrvAlias`. Stored as
an **S-group record** carrying `SGROUP_FLAG_ALIAS`, exactly one atom, and the text in the name slot;
`type` is empty. `mol.aliases` materialises them as `{stable_id: bytes}` and `mol.set_aliases`
replaces them; `mol.sgroups` filters them out, so no caller sees the sharing.

An alias is not an S-group by any reading of the reference, but its storage requirement is an
S-group's exactly: a one-atom label that survives a remap, follows its atom into a substructure, and
is dropped and reported when its atom dies.

---

## 12. Measured dialect facts

Where the reference is silent, implementations differ; the difference is a fact about the file rather
than about any implementation. Measured by running both sides: chython 3 (this tree), RDKit 2026.03.4,
Indigo 1.45.0.0, CDK 2.12, ChemAxon Marvin 25.1.3 (`molconvert`). Harness:
`chython/formats/test/test_conformance.py` with `chython/formats/test/oracles.py`;
`python -m chython.formats.test.oracles matrix` prints the table and every oracle version present.

- **The third field of a V2000 `$RXN` counts line.** The reference states two counts. chython writes
  `  2  2  1` for a two-reactant, two-product, one-agent reaction, and all four read the third field
  as the agent count, reporting 2/1/2.
- **An SD or RD file mixing CTAB versions.** `MR.rdf` in the test corpus carries both V2000 and
  V3000 CTABs in one file, which is why the version is sniffed per CTAB.
- **The V2000 chiral flag and enhanced stereo.** V2000 has no enhanced-stereo block (§1.5.2 is a
  convention on top of `DAT`). For a molecule with two stereocentres, CDK 2.12 and Indigo 1.45 write
  chiral flag 0 with the parities in the atom block and read that back as one AND collection over
  those centres; chython reports a collection only where a file states one, which in CTfile means
  either the V3000 collection block or the `MDLV30/STERAC` `DAT` convention. The formats differ in
  what they can state and the readers differ in what they infer from the flag.
- **Bond type 4 is the de facto aromatic bond.** The reference lists 4 among the query bond types,
  and in practice it is how an aromatic bond is written and read.

  | Writer, on a molecule it perceives as aromatic (pyrrole) | Bond block |
  |---|---|
  | chython 3, Indigo 1.45, `molconvert mol` | `4` |
  | RDKit 2026.03.4 | Kekulé `1`/`2` by default, `4` on `kekulize=False` |
  | CDK 2.12 | Kekulé `1`/`2` |

  Readers, on an order-4 benzene: all five hand back an aromatic ring with one hydrogen per carbon —
  chython, RDKit, Indigo and `molconvert` directly, CDK after its own
  `percieveAtomTypesAndConfigureAtoms` plus `CDKHydrogenAdder` step, which is where CDK fills
  implicit counts. A Kekulé preference is a preference and not a disagreement about what 4 means.
- **`MRV_IMPLICIT_H` and the pyrrole nitrogen.** chython states an aromatic pnictogen's count
  in a `MRV_IMPLICIT_H` data S-group; RDKit 2026.03.4, Indigo 1.45 and `molconvert` consume it
  and report the N-H. CDK 2.12's `MDLV2000Reader` reads the group — `getSgroups()` returns
  `(MRV_IMPLICIT_H, IMPL_H1, 1)` — and states the count from its own configuration rather than from the
  group, the field being a ChemAxon convention and not a CTfile one; its `MDLV3000Reader` does not read
  the group at all (next bullet). A Kekulé form, written after `kekule()`, carries the count without it.
- **`DAT` in a V3000 CTAB.** CDK 2.12's `MDLV3000Reader` reports `Skipping unrecognized SGROUP type:
  DAT`; its `MDLV2000Reader` reads the same group. The same chython data S-group therefore survives
  into CDK through V2000 and not through V3000.
- **A free-text label in the atom-line symbol field.** CDK 2.12 states an atom label on an
  `IPseudoAtom` and writes it in that field — `    0.6495    1.1250    0.0000 Me ` in V2000,
  `M  V30 1 Me 0.64952 1.125 0 0` in V3000 — rather than as an `A  <n>` alias line. §8 lists the
  tokens for that field and says nothing about free text, and the readers differ: CDK reports
  `invalid symbol: Me` and returns a pseudo atom keeping the symbol it was built with; chython keeps
  the atom with the text as its alias, a placeholder element and no hydrogen count (§8). An `A  <n>`
  line for the same atom outranks the column, being a label stated as one.
- **Aliases in a V3000 CTAB.** For an aliased atom, RDKit 2026.03.4, Marvin 25.1.3 and chython all
  emit a V3000 CTAB stating no alias.
- **An atom list.** For one RDKit-written CTAB carrying `L` plus `M  ALS` (V2000) or `[C,N]`
  (V3000), the five readers hand back five kinds of answer. An atom list is a query construct and the
  reference does not state what object a reader builds from one.

  | Reader | Answer |
  |---|---|
  | RDKit 2026.03.4 | a molecule whose atom states both members |
  | Indigo 1.45 | a refusal from `loadMolecule` (atom lists being for queries); both members from `loadQueryMolecule` |
  | CDK 2.12 | a `QueryAtom` with no symbol from V2000; a `PseudoAtom` labelled `R` from V3000 |
  | Marvin 25.1.3 | both members through `-g smarts` |
  | chython 3 | an `UnsupportedCtfile` naming a query reader |
