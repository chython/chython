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
"""Library monomers as lookup keys: the core with unlabelled attachment markers, and the bonds to cut."""

from time import perf_counter

from .._library import MonomerLibrary, get_monomer_library, load_library, DEFAULT_LIBRARY, piece_keys


def test_a_core_key_is_the_monomer_with_its_leaving_groups_as_markers():
    lib = get_monomer_library()
    key, constitution = piece_keys('*N[C@@H](C)C(*)=O')
    assert [(m.symbol, used) for m, used in lib.by_core(key)] == [('A', frozenset((1, 2)))]
    assert piece_keys('*C(=O)[C@H](C)N*')[0] == key                 # any atom order
    assert {m.symbol for m, _ in lib.by_constitution(constitution)} >= {'A', 'dA'}


def test_a_terminal_residue_keeps_its_cap_in_the_key():
    lib = get_monomer_library()
    assert ('A', frozenset((1,))) in [(m.symbol, u) for m, u in lib.by_core(piece_keys('*N[C@@H](C)C(=O)O')[0])]
    assert ('A', frozenset()) in [(m.symbol, u) for m, u in lib.by_core(piece_keys('N[C@@H](C)C(=O)O')[0])]


def test_keys_do_not_depend_on_the_aromatic_spelling():
    assert piece_keys('*c1ccccc1')[0] == piece_keys('*C1=CC=CC=C1')[0]


def test_the_labelled_core_names_each_marker():
    lib = get_monomer_library()
    core, labels = lib.labelled_core(lib.get('PEPTIDE', 'A'), frozenset((1, 2)))
    by_label = {label: core.neighbors_of(atom) for atom, label in labels.items()}
    assert sorted(by_label) == [1, 2]
    (n,) = by_label[1]
    assert core.atom(n).atomic_symbol == 'N'


def test_cut_signatures_hold_the_peptide_bond_ends():
    signatures = get_monomer_library().cut_signatures
    assert ('N', False, 0, 1, 1, 0) in signatures       # a secondary amide nitrogen's side
    assert ('C', False, 0, 0, 2, 1) in signatures       # the carbonyl carbon's side, its O the heteroatom


def test_the_default_index_builds_quickly():
    monomers, log = load_library(DEFAULT_LIBRARY)
    start = perf_counter()
    MonomerLibrary(monomers, log).cut_signatures
    assert perf_counter() - start < 2.0
