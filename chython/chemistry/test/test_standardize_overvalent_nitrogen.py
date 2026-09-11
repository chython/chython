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
"""A neutral nitrogen drawn with four bonds' worth of valence, and which repair it gets.

Two repairs exist and the drawing does not say which is meant.  Where a *terminal oxygen* hangs off the
nitrogen by a single bond, the pair is one charge separation apart and the compound is neutral overall
-- an N-oxide, `C[N+](C)(C)[O-]`.  Otherwise the nitrogen is simply a cation whose charge the drawing
omitted -- N-methylpyridinium.  Choosing by the partner rather than by the spelling is the whole content
of this file, because the two repairs differ in the molecule's *net charge* and picking the wrong one
silently changes the compound.

Only oxygen takes the counter-charge, and that is measured rather than assumed: the aromatic reader
promotes `c1ccn(N)cc1` to the 1-aminopyridinium cation and does not separate it into an aminide, so a
terminal amino group is a substituent like any other.  A row written for the symmetry, separating onto
nitrogen as well, disagreed with the reader on that molecule and was dropped.

`kekule()` already makes this choice, in the reader, for a drawing that has an aromatic system to
resolve; `test_isomers.py::CORPUS_SHAPES` measures those.  The same compounds written Kekule reach no
aromatic system, so the choice has to exist in the rule table too, and both spellings have to land on
one key or `canonical_bytes` deduplicates by drawing again.
"""
import pytest

from .. import canonicalize, check_valence, standardize
from ...core import read_smiles


#: `(aromatic, kekule, label)` -- one compound, two drawings.  Public compounds only.
SPELLINGS = [
    ('c1ccn(O)cc1', 'ON1=CC=CC=C1', 'pyridine N-oxide'),
    ('c1ccn(C)cc1', 'CN1=CC=CC=C1', 'N-methylpyridinium'),
    ('c1ccn(N)cc1', 'NN1=CC=CC=C1', '1-aminopyridinium'),
    ('NC(=O)c1cccn(C)c1', 'NC(=O)C1=CN(C)=CC=C1', 'N-methylnicotinamide'),
]

#: `(smiles, net charge the repair must reach, label)`.  A terminal oxygen keeps the compound neutral;
#: without one the nitrogen was a cation all along and the charge is genuinely new.
NET_CHARGE = [
    ('ON1=CC=CC=C1', 0, 'pyridine N-oxide drawn Kekule'),
    ('[O-]N1=CC=CC=C1', 0, 'the same, with the minus already on the oxygen'),
    ('CN1=CC=CC=C1', 1, 'N-methylpyridinium -- nothing can take a minus'),
    ('NN1=CC=CC=C1', 1, '1-aminopyridinium -- an amino group is a substituent, not an acceptor'),
    ('CN(C)(C)C', 1, 'tetramethylammonium, four-coordinate'),
    ('C[N](C)(C)OC', 1, 'N-methoxy-trimethylammonium -- the oxygen is substituted, not terminal'),
    ('CN(=N)(C)C', 0, 'the four-coordinate double-bonded case groups:33 already separates'),
]

#: Left alone: a three-coordinate nitrogen whose own bonds are all single is `z1`, valence 3, and a
#: perfectly ordinary compound.  The repair keys on `z2`, so these must not move.
UNTOUCHED = [
    ('CC(=O)N(C)O', 'N-hydroxy-N-methylacetamide'),
    ('ON1CCCCC1', 'N-hydroxypiperidine'),
    ('CC(=O)N(C)N', '1-acetyl-1-methylhydrazine'),
    ('NN=C(C)C', 'acetone hydrazone'),
    ('C[N+](C)(C)O', 'N-hydroxy-trimethylammonium -- already a cation, and no violation to repair'),
]


def _violations(mol):
    """Valence violations on the Kekule form.  An aromatic atom answers `unknown` for want of a row."""
    probe = mol.copy()
    probe.kekule()
    return [n for n, kind in check_valence(probe) if kind == 'violation']


def _net_charge(mol):
    return sum(mol.atom(n).charge for n in mol.atom_numbers)


def test_a_kekule_drawn_overvalent_nitrogen_is_repaired():
    """`standardize()` alone, with no aromatic system to lean on, still resolves the valence."""
    for _, kekule, label in SPELLINGS:
        mol = read_smiles(kekule)
        assert _violations(mol), f'{label}: fixture no longer starts broken -- {kekule}'
        standardize(mol)
        assert not _violations(mol), f'{label}: still over-valent after standardize -- {kekule}'


def test_both_spellings_of_one_compound_reach_one_key():
    """The deduplication guarantee: how it was drawn must not survive `canonicalize()`."""
    for aromatic, kekule, label in SPELLINGS:
        a, k = read_smiles(aromatic), read_smiles(kekule)
        canonicalize(a)
        canonicalize(k)
        assert a.canonical_bytes == k.canonical_bytes, f'{label}: {aromatic} vs {kekule}'


def test_the_repair_reaches_the_right_net_charge():
    """A charge separation conserves the total; promoting a cation does not, and must not be applied
    where a separation was available."""
    for string, expected, label in NET_CHARGE:
        mol = read_smiles(string)
        standardize(mol)
        assert _net_charge(mol) == expected, f'{label}: {string} -> {mol} has {_net_charge(mol):+d}'


def test_a_separated_oxide_keeps_no_hydrogen_on_the_oxygen():
    """The proton leaves with the charge.  No rule writes a hydrogen count, so this is `calc_implicit`
    re-deriving the oxygen the patch touched -- and it is the observable half of the repair."""
    mol = read_smiles('ON1=CC=CC=C1')
    standardize(mol)
    oxygens = [n for n in mol.atom_numbers if mol.atom(n).atomic_symbol == 'O']
    assert len(oxygens) == 1
    assert mol.atom(oxygens[0]).charge == -1, str(mol)
    assert mol.implicit_h_of(oxygens[0]) == 0, str(mol)


def test_an_already_anionic_oxide_only_charges_the_nitrogen():
    """Half-separated input: the oxygen carries the minus and the nitrogen was left neutral."""
    mol = read_smiles('[O-]N1=CC=CC=C1')
    assert _net_charge(mol) == -1
    standardize(mol)
    assert not _violations(mol), str(mol)
    assert _net_charge(mol) == 0, str(mol)


def test_a_valid_nitrogen_is_left_alone():
    """The repair keys on `z2`, a nitrogen holding a double bond of its own.  An N-hydroxy amide and an
    N-amino amide are `z1`, valence 3, and ordinary."""
    for string, label in UNTOUCHED:
        mol = read_smiles(string)
        before = mol.canonical_bytes
        standardize(mol)
        assert mol.canonical_bytes == before, f'{label}: {string} -> {mol}'


@pytest.mark.xfail(reason='`groups:10` charges a four-coordinate nitrogen without consulting whether a '
                          'terminal oxygen could take the counter-charge, so trimethylamine N-oxide '
                          'drawn neutral comes out as the hydroxy-ammonium cation instead.  The '
                          'three-coordinate rows cannot help: they run later, and a row matching the '
                          'already-charged product would rewrite the legitimate cation '
                          '`C[N+](C)(C)O`, which has no violation to repair.  Fixing it means '
                          'splitting `groups:10`, a proven union of two chython 2 rules, so that its '
                          'nitrogen half sits after the separating rows -- a change to '
                          '`test_standardize_rules_merges.py` provenance and not to this table',
                   strict=True)
def test_a_four_coordinate_oxide_separates_rather_than_promoting():
    """Trimethylamine N-oxide, drawn without its charges.  The `D4` counterpart of `groups:82`."""
    mol = read_smiles('C[N](C)(C)O')
    standardize(mol)
    assert _net_charge(mol) == 0, str(mol)
