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
"""Polymer types beyond peptides, end to end, with the shipped public library."""

from collections import Counter

from chython.core import smiles

from .. import helm


def formula(mol):
    counts = Counter()
    for n in mol.atom_numbers:
        counts[mol.atom(n).atomic_symbol] += 1
        counts['H'] += mol.implicit_h_of(n)
    return dict(counts)


def test_a_dinucleotide_reads_as_apu():
    # adenylyl-(3'->5')-uridine: adenosine C10H13N5O4 + uridine C9H12N2O6 + H3PO4, less two waters
    m = helm('RNA1{r(A)p.r(U)}$$$$')
    assert formula(m) == {'C': 19, 'H': 24, 'N': 7, 'O': 12, 'P': 1}
    assert m.connected_components_count == 1


def test_a_dinucleotide_round_trips():
    m = helm('RNA1{r(A)p.r(U)}$$$$')
    text = helm(m)
    assert text.startswith('RNA1{') and helm(helm(text)) == text
    assert formula(helm(text)) == formula(m)


def test_a_chem_monomer_joins_a_peptide():
    # an inline SMILES CHEM monomer, a propionyl, on the peptide's N-terminus
    m = helm('PEPTIDE1{A.G}|CHEM1{[CCC(=O)[*:1]]}$CHEM1,PEPTIDE1,1:R1-1:R1$$$')
    assert m == smiles('CCC(=O)N[C@@H](C)C(=O)NCC(=O)O')


def test_a_blob_is_a_marker():
    m = helm('PEPTIDE1{A.G}|BLOB1{Bead}$PEPTIDE1,BLOB1,2:R2-1:R1$$$')
    (marker,) = [a for a in m.atom_numbers if m.atom(a).is_r]
    assert m.aliases[marker] in ('Bead', b'Bead')
