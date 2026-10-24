# -*- coding: utf-8 -*-
#
#  Copyright 2026 Ramil Nugmanov <nougmanoff@protonmail.com>
#  This file is part of chython.
#
#  chython is free software; you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with this program; if not, see <https://www.gnu.org/licenses/>.
#
"""The hand-derived CIP cases, and the two properties that need no reference implementation.

EVERY ROW OF `cip_cases.tsv` IS DERIVED RATHER THAN RECALLED, by one rule applied three times:

  * OpenSMILES `@` means the last three neighbours run anticlockwise viewed FROM the first, so the
    neighbour order `(lowest, r1, r2, r3)` with `@` is R -- CIP looks from the opposite side, which
    reverses the sense;
  * `(lowest, r1, r2, r3)` -> `(r1, r2, r3, lowest)` is a 4-cycle, odd, so in DESCENDING PRIORITY order
    `@` is S and `@@` is R;
  * therefore: write the SMILES neighbour order, write the descending-priority order, and take the parity
    of the permutation between them.  Even leaves `@` = S, odd flips it to `@` = R.

Worked, per compound, with the centre's neighbours in SMILES order (a bracket atom's implicit hydrogen
comes first when nothing precedes the atom, otherwise straight after the preceding neighbour):

| SMILES | neighbours, SMILES order | descending priority | as list indices | parity | `@` |
| --- | --- | --- | --- | --- | --- |
| `[C@H](F)(Cl)Br` | H, F, Cl, Br | Br, Cl, F, H | 4, 3, 2, 1 | even | S |
| `C[C@H](O)CC` | CH3, H, OH, C2H5 | OH, C2H5, CH3, H | 3, 4, 1, 2 | even | S |
| `OC[C@H](O)C=O` | CH2OH, H, OH, CHO | OH, CHO, CH2OH, H | 3, 4, 1, 2 | even | S |
| `C[C@H](N)C(=O)O` | CH3, H, NH2, COOH | NH2, COOH, CH3, H | 3, 4, 1, 2 | even | S |
| `C[C@H]([2H])O` | CH3, H, 2H, OH | OH, CH3, 2H, H | 4, 1, 3, 2 | even | S |
| `[C@H](C)(C=C)C(C)C` | H, CH3, vinyl, iPr | vinyl, iPr, CH3, H | 3, 4, 2, 1 | odd | R |
| `[C@H](C)(C#C)C(C)(C)C` | H, CH3, ethynyl, tBu | ethynyl, tBu, CH3, H | 3, 4, 2, 1 | odd | R |
| `[C@H](O)([13CH3])CC` | H, OH, 13CH3, C2H5 | OH, C2H5, 13CH3, H | 2, 4, 3, 1 | even | S |
| `[P@](=O)(F)(CC)OC` | =O, F, C2H5, OCH3 | F, OCH3, =O, C2H5 | 2, 4, 1, 3 | odd | R |
| `[S@@](=O)(=NC)(CC)OC` | =O, =NCH3, C2H5, OCH3 | OCH3, =O, =NCH3, C2H5 | 4, 1, 2, 3 | odd | `@@` S |

The rankings each row rests on: ethyl over methyl at sphere 2, (C,H,H) against (H,H,H); CHO over CH2OH at
sphere 2, (O, O-duplicate, H) against (O, H, H), which is what the digraph's duplicate atoms buy; COOH
over CH3 at sphere 2, (O, O, O-duplicate) against (H,H,H); and 2H over 1H by rule 2 alone, the one row
rules 1a and 1b are blind to.  The last three rows order the rules: vinyl over isopropyl and ethynyl over
tert-butyl by rule 1a at sphere 3, (C, C-duplicate, H) against (H,H,H), although at sphere 2 a duplicate
faces a real carbon (rule 1b); ethyl over 13CH3 by rule 1a at sphere 2 before rule 2 reads the 13C.  A
double bond at a four-coordinate P(V) or S(VI) counts as single (IUPAC 2013 P-93.2.4), so its =O carries no
duplicate: OCH3 over =O at sphere 2, (C) against nothing.
"""
from pathlib import Path

import pytest

from chython.core import read_smiles

CASES = Path(__file__).parent / 'cip_cases.tsv'
COLUMNS = ['smiles', 'atom', 'descriptor', 'rule', 'compound']


def _rows():
    with CASES.open() as f:
        header = next(f).rstrip('\n').split('\t')
        assert header == COLUMNS
        for line in f:
            if line.strip():
                yield dict(zip(header, line.rstrip('\n').split('\t')))


ROWS = list(_rows())


@pytest.mark.parametrize('row', ROWS, ids=lambda r: r['compound'])
def test_the_derived_descriptor_is_the_computed_one(row):
    mol = read_smiles(row['smiles'])
    mol.assign_cip()
    assert mol.atom_cips().get(int(row['atom'])) == row['descriptor']


@pytest.mark.parametrize('row', ROWS, ids=lambda r: r['compound'])
def test_reflecting_every_parity_flips_every_letter_within_its_case(row):
    """The mirror property, and it needs no oracle: R <-> S at every site the case labels."""
    mol = read_smiles(row['smiles'])
    mol.assign_cip()
    before = mol.atom_cips()
    mirror = read_smiles(row['smiles'])
    # READ EVERY PARITY BEFORE OPENING THE SESSION: a container refuses a read while a journal is
    # pending, so `parity_of` inside the `with` block would raise rather than answer.
    flipped = {n: (1 if mirror.parity_of(n) == 2 else 2)
               for n in mirror.chiral_atoms() if mirror.parity_of(n)}
    with mirror.edit() as e:
        for n, p in flipped.items():
            e.set_parity(n, p)
    mirror.assign_cip()
    after = mirror.atom_cips()
    assert set(before) == set(after)
    flip = {'R': 'S', 'S': 'R'}
    for n, code in before.items():
        assert after[n] == flip[code]


def test_a_computed_descriptor_is_in_the_bytes_and_out_of_the_canonical_form():
    """RULES.md 1.5 item 8, now for a COMPUTED descriptor: `assign_cip` writes the arena in place, and
    that write must not move the identity."""
    plain = read_smiles('[C@H](F)(Cl)Br')
    labelled = read_smiles('[C@H](F)(Cl)Br')
    labelled.assign_cip()
    assert labelled.to_bytes() != plain.to_bytes()
    assert labelled == plain
    assert hash(labelled) == hash(plain)
    assert labelled.atoms_order == plain.atoms_order
