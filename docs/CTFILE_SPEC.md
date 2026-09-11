# CTFile Formats Specification Summary

Extracted from MDL CTFile Formats (October 2003) + BIOVIA 2020 additions.
Focused on molecular structure + stereo information relevant to chython.

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
| lll   | 3     | Number of atom lists (query, ignore)                      |
| fff   | 3     | Obsolete                                                  |
| ccc   | 3     | Chiral flag: 0=not chiral, 1=chiral  (ignore)             |
| sss   | 3     | Number of stext entries (ignore)                          |
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
| sss | 3 | Atom stereo parity (ignored by readers) | 0=none, 1=odd, 2=even, 3=either |
| hhh | 3 | Hydrogen count+1 (query, ignore)        | |
| bbb | 3 | Stereo care box (query, ignore)         | |
| vvv | 3 | Valence                                 | 0=default, 1-14, 15=zero valence |
| HHH | 3 | H0 designator (ignore)                  | |
| rrr | 3 | Not used                                | |
| iii | 3 | Not used                                | |
| mmm | 3 | Atom-atom mapping number                | 0=no mapping, >0=mapped |
| nnn | 3 | Inversion/retention (reaction, ignore)  | 0=none, 1=inverts, 2=retained |
| eee | 3 | Exact change (reaction, ignore)         | 0=none, 1=exact |

**Important**: Fields `dd` and `ccc` are superseded by `M ISO`, `M CHG`, `M RAD` in properties block.

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
| rrr | 3 | Bond topology (query, ignore)   | 0=either, 1=ring, 2=chain |
| ccc | 3 | Reacting center status (ignore) | 0=unmarked, 1=center, -1=not center, 2/4/8/12=changes |

**Stereo convention**: The wedge (pointed) end is at the FIRST atom (field 111).

### 1.5 Properties Block

Terminated by `M  END`. Key properties for molecules:

```
M  CHGnn8 aaa vvv ...       Charge: vvv = -15..+15
M  RADnn8 aaa vvv ...       Radical: 0=none, 1=singlet, 2=doublet, 3=triplet
M  ISOnn8 aaa vvv ...       Isotope: absolute atomic mass (positive integer)
M  END                      End of CTAB
```

When `M CHG`/`M RAD` present, they supersede ALL atom block charge/radical values (forces 0 on unlisted atoms).

#### V2000 Enhanced Stereo (BIOVIA 2020 extension)

Not in the 2003 spec. Uses Sgroup-like property lines:

```
M  STY  1   1 DAT           Define Sgroup 1 as data type
M  SAL   1  n  a1 a2 ...    Atoms in stereo group
M  SDT   1  MDLV30/STERAC1  Field name identifies stereo type
M  SED   1                   (empty data)
```

Types: `MDLV30/STEABS`, `MDLV30/STERAC{n}`, `MDLV30/STEREL{n}` (same semantics as V3000).

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
| nsg | Number of Sgroups (ignore) |
| n3d | Number of 3D constraints (ignore) |
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

Query-only (ignore): HCOUNT, STBOX, SUBST, UNSAT, RBCNT, ATTCHPT, RGROUPS, ATTCHORD.
Reaction-only: INVRET (0/1/2), EXACHG (0/1).

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
| TOPO | Topology (query, ignore) | 0=default, 1=ring, 2=chain |
| RXCTR | Reacting center (reaction) | same as V2000 |
| STBOX | Stereo care box (query, ignore) | |

**Note V3000 vs V2000 bond stereo difference**: V2000 uses `1=Up, 6=Down`. V3000 uses `1=Up, 3=Down`.

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
| `MDLV30/STEABS` | **Absolute** stereo: these centers have known absolute configuration |
| `MDLV30/STERACn` | **Racemic** (AND group n): relative config is known, mixture of both enantiomers present. All atoms with same `n` flip together. |
| `MDLV30/STERELn` | **Relative** (OR group n): relative config is known, only one enantiomer present but which one is unknown. All atoms with same `n` flip together. |

- `n` is an integer >= 1 identifying the group
- Multiple groups can coexist (e.g., STERAC1, STERAC2, STEREL1)
- Atoms not listed in any collection default to ABS
- ATOMS list format: `(count atom_index atom_index ...)`

### 2.7 Sgroup Block (for reference, mostly ignored)

```
M  V30 BEGIN SGROUP
M  V30 index type extindex [ATOMS=(...)] [FIELDNAME=name] [FIELDDATA=data] ...
M  V30 END SGROUP
```

Types: SUP(eratom), MUL(tiple), SRU, MON(omer), COP(olymer), DAT(a), etc.

Only relevant for chython: `DAT` type SGROUPs with `FIELDNAME="MRV_IMPLICIT_H"` and `FIELDDATA` specifying implicit H count.

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
The `>` must be in column 1. Field name is in `<angle_brackets>`.

### V2000 vs V3000 detection:
Check line 4 (counts line): if it starts with `M  V30 BEGIN CTAB` then V3000, else V2000.
(Line indices: 0=name, 1=program, 2=comment, 3=counts/begin, 4=first CTAB line for V3000)

Actually, the standard detection: check if the version stamp in the counts line (chars 34-39) is `V3000`.
Or pragmatically: line[4] (0-indexed from start of MOL block) starting with `M  V30`.

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
  - Extended: `rrrpppaaa` with optional agent count
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

For **single bonds** at tetrahedral centers:
- **Up (wedge)**: bond goes from center atom (first atom) ABOVE the plane toward second atom
- **Down (hash/dash)**: bond goes from center atom BELOW the plane toward second atom
- The pointed end of the wedge is at the **first atom** in the bond definition

For **double bonds** (cis/trans):
- Stereo determined from 2D coordinates of substituents
- Value 0 = use coordinates; value 3 = either (ignore stereo)

### 7.2 Atom Stereo Parity (V3000 CFG on atoms)

Calculated by viewing the center from behind the highest-numbered neighbor:
- **1 = odd parity**: atoms 1,2,3 in clockwise order
- **2 = even parity**: atoms 1,2,3 in counterclockwise order
- **3 = either**: unmarked or racemic at that center

Note: In V2000, the `sss` field in atom block is IGNORED by readers. Stereo is determined from bond wedges + coordinates. In V3000, atom `CFG` is informational; bond `CFG` (wedge/hash) is what defines stereo.

### 7.3 Enhanced/Extended Stereo

Defines groups of stereocenters with specific epistemic relationships:

| Type | Meaning | Example |
|------|---------|---------|
| **ABS** (absolute) | Configuration is exactly as drawn | Single known enantiomer |
| **AND** (racemic, STERAC) | Relative config known, mixture of both enantiomers | Racemic drug |
| **OR** (relative, STEREL) | Relative config known, single enantiomer, but absolute config unknown | Natural product isolate |

Group numbering allows multiple independent groups:
- STERAC1, STERAC2 = two independent AND groups (flip independently)
- STEREL1, STEREL2 = two independent OR groups

Within a single group, all atoms flip together (their relative configuration is fixed).

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
| `D` | Deuterium | Convert to H with MASS=2 |
| `T` | Tritium | Convert to H with MASS=3 |
| `A` | Any atom (query) | Reject or skip |
| `Q` | Any non-C non-H (query) | Reject or skip |
| `L` | Atom list (query) | Reject or skip |
| `R#` | R-group label | Reject or skip |
| `LP` | Lone pair | Reject or skip |

---

## 9. Bond Types

| Value | Meaning | Notes |
|-------|---------|-------|
| 1 | Single | |
| 2 | Double | |
| 3 | Triple | |
| 4 | Aromatic | Query in V2000 spec, but used in practice |
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

Format: `ENDPTS=(N atom1 atom2 ... atomN)` where N is count of endpoint atoms.

In practice: create a bond (type 8/special) from the real atom to each endpoint atom.
