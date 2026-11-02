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
"""`nucleotides()`: units, links and the 5' to 3' strand of public oligonucleotides."""
import pytest

from chython.core import smiles
from chython.core.chains import perceive
from chython.core.nucleotides import nucleotides

THYMIDINE = 'C[C@@H]1O[C@@H](n2cc(C)c(=O)[nH]c2=O)C[C@@H]1O'          # 5'-C, then the 3'-O of each unit
OMEU = 'C[C@@H]1O[C@@H](n2ccc(=O)[nH]c2=O)[C@@H](OC)[C@@H]1O'        # 2'-O-methyluridine
STRANDS = {
    'dT4': 'O' + 'P(=O)(O)O'.join([THYMIDINE] * 4),                 # tetrathymidylate, phosphodiester
    'mU4 PS': 'O' + 'P(=O)(S)O'.join([OMEU] * 4),                   # 2'-OMe uridines, phosphorothioate
}


@pytest.mark.parametrize('name', sorted(STRANDS))
def test_a_strand_is_its_units_from_5_to_3(name):
    mol = smiles(STRANDS[name])
    strand = nucleotides(mol)
    assert strand.kind == 'nucleic' and len(strand.residues) == len(strand.main) == 4
    assert set().union(*strand.members) == set(mol)
    first, *rest = strand.residues
    assert mol.atom(first.path[0]).atomic_symbol == 'O'               # the 5' end: no phosphorus before it
    assert all(mol.atom(r.path[0]).atomic_symbol == 'P' for r in rest)
    for a, b in zip(strand.residues, rest):                          # each link: O3' to the next unit's P
        assert b.path[0] in mol.neighbors_of(a.path[-1])


def test_a_strand_is_perceived_after_peptides():
    assert perceive(smiles(STRANDS['dT4'])).kind == 'nucleic'


def test_three_units_are_no_strand():
    assert nucleotides(smiles('O' + 'P(=O)(O)O'.join([THYMIDINE] * 3))) is None
