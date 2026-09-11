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
"""Does each ported rule fire on the same molecules as a pinned out-of-tree chython 2.24, and produce
the same product?  A disagreement is not automatically V3 being wrong: triage it and name the
deliberate divergence below rather than relaxing the comparison.

Version pinning, isolation of the child and skipping when chython 2 is absent live in
`chython.chemistry.test._oracle`, along with the silent ways this comparison can be made worthless.
"""
import pytest

from chython.chemistry.test._oracle import ask, needs_oracle
from chython.core import read_smiles
from ._corpus import CORPUS
import chython.chemistry                                                         # noqa: F401


# Run inside the oracle interpreter.  Deliberately tiny: everything it prints is data, since analysis
# done there is analysis this repository cannot see or maintain.
_SCRIPT = """
from chython import smiles

for line in sys.stdin:
    line = line.strip()
    if not line:
        continue
    try:
        m = smiles(line)
        m.standardize()
        sys.stdout.write('%s\\t%s\\n' % (line, format(m, 's')))
    except Exception as e:
        sys.stdout.write('%s\\tERROR %s: %s\\n' % (line, type(e).__name__, e))
"""


@pytest.fixture(scope='module')
def oracle_answers():
    """`{input smiles: standardized smiles}` from the pinned chython 2, in one subprocess."""
    answers = dict(line.partition('\t')[::2] for line in ask(_SCRIPT, '\n'.join(CORPUS)))
    assert len(answers) == len(CORPUS)
    return answers


@needs_oracle
def test_the_oracle_is_the_version_the_rules_were_ported_against(oracle_answers):
    """Reaching the fixture at all is the assertion: `ask` checks the version and the import path.

    Named so a broken oracle fails once with a reason, not as 79 apparent chemistry regressions.
    """
    assert oracle_answers


@needs_oracle
def test_every_corpus_molecule_standardizes_to_what_chython_2_produces(oracle_answers):
    """All 79, compared AS MOLECULES: a string comparison would fail on atom order and on which
    resonance partner carries the charge, which are writer differences and not repair differences.
    """
    disagreements = []
    for source, expected in sorted(oracle_answers.items()):
        assert not expected.startswith('ERROR'), f'the oracle cannot read its own corpus: {source}'
        molecule = read_smiles(source)
        molecule.standardize()
        got = molecule.smiles
        if read_smiles(expected) != read_smiles(got):
            disagreements.append(f'{source}\n    chython 2: {expected}\n    chython 3: {got}')
    assert not disagreements, (
        f'{len(disagreements)} of {len(CORPUS)} molecules standardize differently.  Triage each one '
        'and decide which engine is right -- do not relax this test:\n' + '\n'.join(disagreements))


@needs_oracle
def test_the_comparison_can_fail():
    """Negative control: a pass that did nothing would agree with the oracle everywhere, so "no
    disagreements" is evidence only while something is actually being repaired.
    """
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
    """Standardizing twice must be standardizing once.  No oracle needed: a self-consistency.

    Catches the `D` translation: the core counts a dative bond in an atom's degree where chython 2
    does not, so a rule meaning "four single bonds, therefore a formal charge" can fire on an atom
    whose fourth bond the pass itself just made dative.
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
