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
"""`helm()` and `biln()`: a string reads, a molecule writes."""

from pytest import raises

from chython.core import smiles

from ._fixtures import OXYTOCIN_HELM, OXYTOCIN_SMILES, rules, same
from .. import biln, helm


def test_helm_both_ways():
    m = helm(OXYTOCIN_HELM)
    assert same(m, smiles(OXYTOCIN_SMILES))
    assert helm(m) == OXYTOCIN_HELM


def test_biln_both_ways():
    text = biln(helm(OXYTOCIN_HELM))
    assert text == 'C(1,3)-Y-I-Q-N-C(1,3)-P-L-G-am'
    assert same(biln(text), smiles(OXYTOCIN_SMILES))


def test_a_read_logs_onto_the_molecule_and_the_list():
    log = []
    m = helm('PEPTIDE1{A.[NoSuch].G}$$$$', log=log)
    assert 'polymers:monomer-unknown' in rules(log)
    assert 'polymers:monomer-unknown' in {r.rule for r in m.log}


def test_a_write_logs_its_fallbacks():
    log = []
    text = helm(smiles('C[C@H](N)C(=O)NCC(=O)OCC1CCCCC1'), log=log)
    assert text.startswith('PEPTIDE1{A.[') and 'polymers:inline-monomer' in rules(log)


def test_neither_takes_anything_else():
    with raises(TypeError):
        helm(42)
    with raises(TypeError):
        biln(None)
