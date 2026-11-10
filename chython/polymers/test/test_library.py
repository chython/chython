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
"""The monomer library: the shipped public set, attaching another, and what loading does with damage."""

import json
from subprocess import run
from sys import executable

from .._library import get_monomer_library, load_library, set_monomer_library


def test_default_library_is_the_public_set():
    lib = get_monomer_library()
    a = lib.get('PEPTIDE', 'A')
    assert a.labels == (1, 2) and a.natural_analog == 'A'
    assert lib.get('RNA', 'r') is not None              # ribose; the public set spells sugars lower case
    assert lib.get('PEPTIDE', 'no-such') is None


def test_an_attached_library_overrides_by_symbol(tmp_path):
    p = tmp_path / 'mine.json'
    p.write_text(json.dumps([{'symbol': 'A', 'polymerType': 'PEPTIDE', 'monomerType': 'Backbone',
                              'smiles': 'N([H:1])CC([OH:2])=O', 'naturalAnalog': 'G', 'rgroups': []}]))
    set_monomer_library(p)
    try:
        assert get_monomer_library().get('PEPTIDE', 'A').natural_analog == 'G'
        assert get_monomer_library().get('PEPTIDE', 'L') is not None     # stacked on the default
    finally:
        set_monomer_library()
    assert get_monomer_library().get('PEPTIDE', 'A').natural_analog == 'A'


def test_damaged_and_duplicate_records_are_logged_not_raised(tmp_path):
    p = tmp_path / 'bad.json'
    p.write_text(json.dumps([
        {'symbol': 'X1', 'polymerType': 'PEPTIDE', 'smiles': 'not a smiles'},
        {'symbol': 'X2', 'polymerType': 'PEPTIDE', 'smiles': 'C[OH:1]'},
        {'symbol': 'X2', 'polymerType': 'PEPTIDE', 'smiles': 'CC[OH:1]'}]))
    monomers, log = load_library(p)
    assert [m.symbol for m in monomers] == ['X2'] and monomers[0].smiles == 'C[OH:1]'
    assert {r.rule for r in log} == {'polymers:monomer-unreadable', 'polymers:monomer-duplicate'}


def test_a_record_with_no_attachment_is_unreadable(tmp_path):
    p = tmp_path / 'flat.json'
    p.write_text(json.dumps([{'symbol': 'Z', 'polymerType': 'CHEM', 'smiles': 'CCO'}]))
    monomers, log = load_library(p)
    assert not monomers and [r.rule for r in log] == ['polymers:monomer-unreadable']


def test_nothing_loads_at_import():
    code = ('import chython.polymers, chython.polymers._library as l; '
            'print(len(l._LIBRARY_CACHE))')
    done = run([executable, '-c', code], capture_output=True, text=True)
    assert done.returncode == 0, done.stderr
    assert done.stdout.strip() == '0'
