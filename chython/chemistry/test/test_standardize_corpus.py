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
"""The standardization pass over the drawn-wrong corpus: it repairs, it is idempotent, and the
dative-bond rules leave the donor neutral.
"""
import pytest

from chython.core import read_smiles
from ._corpus import CORPUS
import chython.chemistry                                                         # noqa: F401


def test_the_pass_repairs_most_of_the_corpus():
    """Every corpus molecule is drawn wrong, so the pass repairs most of them."""
    changed = 0
    for source in CORPUS:
        molecule = read_smiles(source)
        if molecule.standardize():
            changed += 1
    assert changed > len(CORPUS) // 2, (
        f'only {changed} of {len(CORPUS)} molecules were repaired; this corpus exists because every '
        'one of them is drawn wrong, so the pass has stopped working rather than the corpus having '
        'become clean')


def test_the_pass_is_idempotent():
    """Standardizing twice must be standardizing once.

    Catches the `D` translation: the core counts a dative bond in an atom's degree, so a rule meaning
    "four single bonds, therefore a formal charge" can fire on an atom whose fourth bond the pass
    itself just made dative.
    """
    drifting = []
    for source in CORPUS:
        molecule = read_smiles(source)
        molecule.standardize()
        once = molecule.smiles
        if molecule.standardize():
            drifting.append(f'{source}\n    after one pass:  {once}\n    after two: {molecule.smiles}')
    assert not drifting, (
        f'{len(drifting)} molecules changed on a second pass:\n' + '\n'.join(drifting))


@pytest.mark.parametrize('source,expected', [
    # A dative bond and a formal charge both say the nitrogen donated its lone pair.  Writing both
    # counts the same electron pair twice, and the hydrogen that follows from the charge is a phantom.
    ('C[N](C)(C)[Fe]', 'C[N](~[Fe+])(C)C'),
    # The same claim from the boron side: three covalent bonds and an incoming donation is neutral.
    ('C[N](C)(C)B(C)(C)C', 'CB(C)(C)~N(C)(C)C'),
    # And a phosphine ligand is not a phosphonium.
    ('C[P](C)(C)[Fe]', 'C[P](~[Fe])(C)C'),
])
def test_a_dative_bond_does_not_make_its_donor_charged(source, expected):
    """The `D` translation, stated as chemistry rather than as a primitive.

    chython 2's `D` describes the covalent neighbourhood (it skips order 8 before counting); the core
    counts every edge.  The three rules whose whole content is a covalent count therefore write their
    bonds out explicitly -- see `DEGREE_MAP` in `gen_standardize_rules.py`.
    """
    molecule = read_smiles(source)
    molecule.standardize()
    assert read_smiles(molecule.smiles) == read_smiles(expected), molecule.smiles


def test_a_metal_carbonyl_survives_the_hydrogen_recompute():
    """Three carbonyls on one iron: all three repaired, and the valence collection not asked to guess.

    Two constraints meet here.  The order-8 policy in `_implicit.py` must hold, or the valence
    collection refuses the dative environment the pass just created and raises.  And dedupe must not
    be keyed on the whole match, or only the first carbonyl is repaired -- all three share the iron.
    """
    molecule = read_smiles('C(=O)[Fe](C=O)C=O')
    assert molecule.standardize()
    # every carbon is now a carbonyl anion with a dative contact to the metal, three times over
    assert molecule.smiles.count('~') == 3, molecule.smiles
    assert read_smiles(molecule.smiles) == read_smiles('[O+]#[C-]~[Fe](~[C-]#[O+])~[C-]#[O+]'), \
        molecule.smiles


def test_two_cyclopentadienyl_rings_leave_the_iron_doubly_charged():
    """Ferrocene, where the shared-anchor rule earns its keep.

    Each ring takes one electron from the metal, so the metal's `+1` applies twice, once per ring.
    That needs `[M]` to stay a shared anchor even though the patch writes it.
    """
    molecule = read_smiles('[Fe]12345678(C9C1C6C4C39)C1C2C7C5C81')
    assert molecule.standardize()
    assert '[Fe+2]' in molecule.smiles, molecule.smiles
