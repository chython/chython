#!/usr/bin/env python3
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
"""Transpose the element and isotope tables into _elements.pxi, refusing a pair that disagree.

    python chython/core/test/gen_element_tables.py

THE INPUT IS THE AUTHORITY.  `chython/core/elements.tsv` and `chython/core/isotopes.tsv` are
maintained data files, edited by hand, and nothing regenerates them -- there is no upstream to
regenerate them from, because an element's mass is measured rather than computed.  This script
transposes them into the C arrays the core reads, and checks the relationships those arrays cannot
express.  It is the only direction that exists.

WHY THE DATA IS NOT SIMPLY WRITTEN AS C ARRAYS.  Six arrays have to agree: `ISOTOPE_OFFSETS[119]`
must be the exact prefix sum of `ISOTOPE_COUNTS[119]`, which must be the exact run lengths of three
parallel 436-entry arrays, and `MDL_ISOTOPE[119]` has to name a mass number that appears in them.
As C initialisers that is about 1500 magic numbers in which changing one abundance means counting
commas across a wrapped block, no relationship is stated anywhere, and nothing checks any of it.
`assert_invariants` states each of those relationships and refuses a pair that breaks one.
"""
from __future__ import annotations

import pathlib
import sys
from collections import Counter
from math import fsum


ROOT = pathlib.Path(__file__).resolve().parent.parent

# run as a script, sys.path[0] is this directory and the repo is not on it at all, so `chython`
# would resolve to whatever is installed -- a different checkout's tables, silently
if str(ROOT.parent.parent) not in sys.path:
    sys.path.insert(0, str(ROOT.parent.parent))

from chython.core.test.gen_valence_rules import emit_array   # noqa: E402  -- needs sys.path

ELEMENTS_TSV = ROOT / 'elements.tsv'
ISOTOPES_TSV = ROOT / 'isotopes.tsv'
PXI = ROOT / '_elements.pxi'
BEGIN = '# --- BEGIN GENERATED TABLES: python chython/core/test/gen_element_tables.py ---'
END = '# --- END GENERATED TABLES ---'

ELEMENTS_HEADER = ('z', 'symbol', 'mdl_isotope', 'valence_electrons', 'atomic_radius')
ISOTOPES_HEADER = ('z', 'symbol', 'mass_number', 'exact_mass', 'abundance')

# `exact_mass` for a row whose mass nobody has.  Written as a token rather than as 0.0 so that
# "nobody has measured this" and "this weighs nothing" are different strings in the file.
UNKNOWN = '?'


class Element:
    """One row of elements.tsv.  `valence_electrons` is None when the file said `?`.

    `atomic_radius` is never None: the column states a number for all 118 rows, per the header there.
    """
    __slots__ = ('z', 'symbol', 'mdl_isotope', 'valence_electrons', 'atomic_radius')

    def __init__(self, z, symbol, mdl_isotope, valence_electrons, atomic_radius):
        self.z = z
        self.symbol = symbol
        self.mdl_isotope = mdl_isotope
        self.valence_electrons = valence_electrons
        self.atomic_radius = atomic_radius

    def row(self):
        return (str(self.z), self.symbol, str(self.mdl_isotope),
                UNKNOWN if self.valence_electrons is None else str(self.valence_electrons),
                repr(self.atomic_radius))


class Isotope:
    """One row of isotopes.tsv.  `mass` is None when the file said `?`."""
    __slots__ = ('z', 'symbol', 'mass_number', 'mass', 'abundance')

    def __init__(self, z, symbol, mass_number, mass, abundance):
        self.z = z
        self.symbol = symbol
        self.mass_number = mass_number
        self.mass = mass
        self.abundance = abundance

    def row(self):
        return (str(self.z), self.symbol, str(self.mass_number),
                UNKNOWN if self.mass is None else repr(self.mass), repr(self.abundance))


def _rows(path, header):
    for lineno, line in enumerate(path.read_text(encoding='utf-8').splitlines(), 1):
        if not line.strip() or line.lstrip().startswith('#'):
            continue
        fields = line.split('\t')
        if tuple(fields) == header:
            continue
        if len(fields) != len(header):
            raise ValueError(f'{path}:{lineno}: {len(fields)} fields, expected {len(header)}')
        yield lineno, fields


def read_elements(path=ELEMENTS_TSV):
    out = []
    for lineno, (z, symbol, mdl, valence, radius) in _rows(path, ELEMENTS_HEADER):
        if radius == UNKNOWN:
            raise ValueError(f'{path}:{lineno}: {symbol} states no atomic radius.  The column has no '
                             f'unknown -- where the calculated set stops, the row carries the group '
                             f'analogue one period up, per the header of elements.tsv.')
        out.append(Element(int(z), symbol, int(mdl),
                           None if valence == UNKNOWN else int(valence), float(radius)))
    return out


def read_isotopes(path=ISOTOPES_TSV):
    out = []
    for lineno, (z, symbol, a, mass, abundance) in _rows(path, ISOTOPES_HEADER):
        out.append(Isotope(int(z), symbol, int(a), None if mass == UNKNOWN else float(mass),
                           float(abundance)))
    return out


def assert_invariants(elements, isotopes):
    """Every relationship between the six C arrays that the arrays themselves cannot state, as a
    refusal to compile.

    Two of them are the argument for generating the pair of files at all: nothing in a hand-written
    array stops `mdl_isotope` naming a mass number with no row, and nothing connects the offset
    array to the count array.
    """
    if [e.z for e in elements] != list(range(1, 119)):
        raise ValueError('elements.tsv must hold atomic numbers 1..118, once each, in order')
    symbols = {e.z: e.symbol for e in elements}
    if len(set(symbols.values())) != 118:
        dupes = [s for s, n in Counter(symbols.values()).items() if n > 1]
        raise ValueError(f'elements.tsv has a repeated symbol: {dupes}')
    mdl = {(e.z, e.mdl_isotope) for e in elements if e.mdl_isotope}

    # The 28 rows that state no valence electron count, as an extent rather than a scatter of `?`:
    # 4f (Ce..Lu) and 5f (Th..Lr).  La and Ac are group 3 in every layout and are not among them.
    f_block = frozenset(range(58, 72)) | frozenset(range(90, 104))
    for e in elements:
        if e.valence_electrons is None:
            if e.z not in f_block:
                raise ValueError(f'{e.symbol} (Z={e.z}) states `{UNKNOWN}` for its valence electron '
                                 f'count and is not in the f block, where the doubt lives; the '
                                 f'group number convention gives every other element a number')
        elif e.z in f_block:
            raise ValueError(f'{e.symbol} (Z={e.z}) is in the f block and states '
                             f'{e.valence_electrons} valence electrons; chython states no count '
                             f'there -- see the header of elements.tsv')
        elif not 1 <= e.valence_electrons <= 12:
            raise ValueError(f'{e.symbol} (Z={e.z}) states {e.valence_electrons} valence electrons; '
                             f'the group number convention bounds the column at 1..12')

    # The column is angstroms, and the one error it is exposed to is a value written in picometres --
    # a hundredfold, so any bound at all catches it.  0.3 is under helium's 0.31 and 3.0 is over
    # caesium's 2.98, the narrowest and the widest atom in the table.
    for e in elements:
        if not 0.3 <= e.atomic_radius <= 3.0:
            raise ValueError(f'{e.symbol} (Z={e.z}) states an atomic radius of {e.atomic_radius}; '
                             f'the column is angstroms and is bounded at 0.3..3.0, between helium '
                             f'and caesium -- {e.atomic_radius} is picometres or a typo')

    keys = [(i.z, i.mass_number) for i in isotopes]
    if keys != sorted(keys):
        raise ValueError('isotopes.tsv must be sorted by (z, mass_number) -- the sort IS the '
                         'compiled layout, so a reordered file is a different set of arrays')
    if len(set(keys)) != len(keys):
        dupes = [k for k, n in Counter(keys).items() if n > 1]
        raise ValueError(f'isotopes.tsv repeats a nuclide: {dupes}')

    by_z = {}
    for i in isotopes:
        # a row with no mass carries no mass and no weight, so the only thing it can be doing is
        # existing -- and the one caller that needs a row to merely exist is the MDL reference below
        if i.mass is None and (i.z, i.mass_number) not in mdl:
            raise ValueError(f'{i.symbol}-{i.mass_number}: `{UNKNOWN}` for the mass says nobody '
                             f'has measured this nuclide, and elements.tsv does not name it as an '
                             f'MDL reference mass number either, so nothing needs the row -- give '
                             f'it a mass or delete it')
        if i.symbol != symbols.get(i.z):
            raise ValueError(f'isotopes.tsv says Z={i.z} is {i.symbol!r}, elements.tsv says '
                             f'{symbols.get(i.z)!r}; there is one symbol table, not two')
        by_z.setdefault(i.z, []).append(i)

    for z, rows in by_z.items():
        # `fsum`: `sum` accumulates floats in extended precision from 3.12 and naively before it, and a
        # column that misses 1.0 by exactly the tolerance is then refused on one interpreter and
        # compiled on another.  A correctly rounded total gives every interpreter the same verdict.
        total = fsum(r.abundance for r in rows)
        if total != 0.0 and abs(total - 1.0) > 1e-6:
            raise ValueError(f'{symbols[z]} (Z={z}) abundances sum to {total!r}; a set of natural '
                             f'abundances sums to 1.0, and an element with none sums to 0.0')
        if len(rows) > 255:
            raise ValueError(f'{symbols[z]} has {len(rows)} rows; ISOTOPE_COUNTS is uint8_t')

    for e in elements:
        if e.mdl_isotope and not any(r.mass_number == e.mdl_isotope for r in by_z.get(e.z, ())):
            raise ValueError(
                f'elements.tsv gives {e.symbol} (Z={e.z}) the MDL reference mass number '
                f'{e.mdl_isotope}, and isotopes.tsv has no row for it.  A file is entitled to '
                f'state that isotope -- it is the one MDL itself hands out -- so `element_mass` '
                f'would answer 0.0 for it.  Add the row, with `{UNKNOWN}` for the mass if no '
                f'measured mass exists.')

    if len(isotopes) > 65535:
        raise ValueError(f'{len(isotopes)} rows; ISOTOPE_OFFSETS is uint16_t')


def compile_tables(elements, isotopes):
    """The generated block: the symbol tuple, the MDL reference table, and the flat isotope arrays.

    Layout is flat parallel arrays behind a 119-entry offset/count index because the access pattern
    is "give me all isotopes of element Z" -- one `range(offset, offset + count)` scan, no hashing.
    Index 0 of the per-element arrays is unused so that the index IS the atomic number.

    `ISOTOPE_OFFSETS` is the prefix sum of `ISOTOPE_COUNTS` and is emitted anyway rather than
    computed at load: it is read on the mass path, and one redundant 238-byte table generated from
    the same rows in the same pass cannot disagree with its source the way a hand-written one did.
    """
    offsets = [0] * 119
    counts = [0] * 119
    numbers = []
    masses = []
    abundances = []

    offset = 0
    for z in range(1, 119):
        rows = [i for i in isotopes if i.z == z]
        offsets[z] = offset
        counts[z] = len(rows)
        offset += len(rows)
        for r in rows:
            numbers.append(r.mass_number)
            masses.append(0.0 if r.mass is None else r.mass)
            abundances.append(r.abundance)

    assert offsets[1:] == [sum(counts[:z]) for z in range(1, 119)], 'offsets are the prefix sum'

    unknown = [f'{i.symbol}-{i.mass_number}' for i in isotopes if i.mass is None]

    lines = [
        BEGIN,
        '# Compiled from chython/core/elements.tsv and chython/core/isotopes.tsv, which are the',
        '# authority and are maintained by hand.  Do not edit here -- run the command above.',
        '# Both files are in compiled order, so the k-th entry here is the k-th row there.',
        '#',
        f'# {len(isotopes)} nuclides over 118 elements.  Nobody has a mass for '
        f'{len(unknown)} of them',
        f'# ({", ".join(unknown)}), and those compile to 0.0.',
        '#',
        '# Every mass number in MDL_ISOTOPE has a row among them, and the compile step refuses a',
        '# pair of tables where one does not: a file is entitled to state the isotope MDL itself',
        '# hands out, and a missing row makes `element_mass` answer 0.0 for it -- an atom of',
        '# bromine-80, the mass number every MDL bromine measures against, weighing nothing.',
        'cdef tuple SYMBOLS = (',
    ]
    for i in range(0, 118, 10):
        chunk = ', '.join(repr(e.symbol) for e in elements[i:i + 10])
        lines.append(f'    {chunk},' if i + 10 < 118 else f'    {chunk})')

    lines += ['', '', 'cdef extern from *:', '    """']
    lines += emit_array('MDL_ISOTOPE', 'unsigned short',
                        [0] + [e.mdl_isotope for e in elements], 12,
                        'mass number MDL measures its mass-difference field from; index 0 unused')
    lines += emit_array('VALENCE_ELECTRONS', 'unsigned char',
                        [0] + [0 if e.valence_electrons is None else e.valence_electrons
                               for e in elements], 20,
                        'group number convention; 0 is the f block, which states none; index 0 unused')
    lines += emit_array('ATOMIC_RADIUS', 'double',
                        [repr(0.0)] + [repr(e.atomic_radius) for e in elements], 8,
                        'calculated atomic radius in angstroms; index 0 unused')
    lines += emit_array('ISOTOPE_OFFSETS', 'unsigned short', offsets, 12,
                        'first row of element Z in the flat arrays; prefix sum of ISOTOPE_COUNTS')
    lines += emit_array('ISOTOPE_COUNTS', 'unsigned char', counts, 20,
                        'how many rows element Z has')
    lines += emit_array('ISOTOPE_NUMBERS', 'unsigned short', numbers, 15,
                        'mass number of the k-th row')
    lines += emit_array('ISOTOPE_MASSES', 'double', [repr(m) for m in masses], 6,
                        'exact mass in daltons; 0.0 where none is known')
    lines += emit_array('ISOTOPE_ABUNDANCES', 'double', [repr(a) for a in abundances], 6,
                        'natural terrestrial fraction; 0.0 for a nuclide with none')
    lines += [
        '    """',
        '    const uint16_t MDL_ISOTOPE[119]',
        '    const uint8_t  VALENCE_ELECTRONS[119]',
        '    const double   ATOMIC_RADIUS[119]',
        '    const uint16_t ISOTOPE_OFFSETS[119]',
        '    const uint8_t  ISOTOPE_COUNTS[119]',
        f'    const uint16_t ISOTOPE_NUMBERS[{len(numbers)}]',
        f'    const double   ISOTOPE_MASSES[{len(masses)}]',
        f'    const double   ISOTOPE_ABUNDANCES[{len(abundances)}]',
        '',
        f'DEF ISOTOPE_ROWS = {len(isotopes)}',
        END,
    ]
    return '\n'.join(lines)


def rewrite_pxi(block, path=PXI):
    text = path.read_text(encoding='utf-8')
    start = text.index(BEGIN)
    stop = text.index(END) + len(END)
    if text[start:stop] == block:
        return False
    path.write_text(text[:start] + block + text[stop:])
    return True


def main(argv):
    if argv:                                     # there is one direction, so there are no verbs
        print(__doc__)
        return 2
    elements = read_elements()
    isotopes = read_isotopes()
    assert_invariants(elements, isotopes)
    if rewrite_pxi(compile_tables(elements, isotopes)):
        print(f'{PXI.name} updated; rebuild the extension')
    else:
        print(f'{PXI.name} already matches the tables')
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
