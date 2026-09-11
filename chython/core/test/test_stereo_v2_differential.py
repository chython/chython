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
"""Differential stereo perception against a pinned chython 2.24, as an outside oracle.

chython 2 is a *behaviour oracle*, never a contract.  Where the two agree it is the cheapest available
check that this rewrite did not lose something; where V3 diverges it does so on purpose and the
divergence is listed below with its reason.  Nothing here treats agreement as obligatory.

**Why a subprocess and not an import.**  The comparison has to outlive the code it compares against:
an in-tree differential dies with the copy it imports, silently becoming a test of V3 against V3.
Shelling out to an isolated interpreter that holds its own pinned copy keeps the oracle meaningful and
pins WHICH V2 is being believed.

**What is compared.**  A count triple per molecule -- (tetrahedral, axial, planar) stereocentres --
and nothing finer.  Atom numbering, parity conventions and frame order all differ between the two
implementations, so any comparison of those would be measuring the translation and not the
chemistry.  Counts of each kind are the largest thing that means the same on both sides.

**Why the corpus carries no stereo descriptors.**  V2's `chiral_*` properties report centres that
are stereogenic AND still unlabelled, whereas V3's `stereogenic_units()` reports stereogenic
centres whether or not they carry a configuration.  On an unlabelled molecule the two coincide; on
a labelled one they cannot, and the difference would be bookkeeping rather than perception.
`C/C=C/CC` is the witness: it disagrees for exactly that reason and no other.
"""
from functools import lru_cache

import pytest

from chython.core import read_smiles, SU_TETRA, SU_CIS_TRANS, SU_ALLENE
from . import oracle
from .oracle import requires_oracle


# ------------------------------------------------------------------------------------------------
# The corpus.  Public, commonplace compounds only, and no stereo descriptors (see the module
# docstring).  Grouped by what each group is here to protect.

_ALLENES_AND_CUMULENES = [
    'CC=C=CC',            # 2,3-pentadiene -- the textbook chiral allene
    'CC=C=C=CC',          # 2,3,4-hexatriene -- even chain, so planar
    'CC=C=C=C=CC',        # heptatetraene -- odd again, so axial
    'CC=C=C=C=C=CC',      # octapentaene -- even
    'C=C=C',              # allene itself: symmetric, nothing to name
    'CC(C)=C=C',          # 3-methyl-1,2-butadiene: one terminal symmetric
    'FC(Cl)=C=C(Br)I',    # two heavy substituents at each terminal
    'CCC=C=CC',           # 2,3-hexadiene
    'OC=C=CO',            # 1,3-propadiene-1,3-diol
    'CC(F)=C=C(F)C',
    'C=CC=C',             # 1,3-butadiene: conjugated, NOT cumulated
    'CC=CC',              # 2-butene: the two-atom rung
]

_SPIRO = [
    'OC1CCCC12CCCC2O',    # spiro[4.4]nonane-1,6-diol
    'C1CCCC12CCCC2',      # spiro[4.4]nonane: symmetric
    'OC1CCC12CCC2O',      # spiro[3.3]heptane-diol
    'C1CC2(CC1)CCCCC2',   # spiro[4.5]decane
    'OC1CCCCC12CCCCC2O',  # spiro[5.5]undecane-diol
]

_AMIDES_AND_CARBONYLS = [
    'CC(=O)NC', 'CC(=O)N(C)C', 'CC(=O)N(C)CC', 'CC(=O)Nc1ccccc1',
    'CNC(=O)NC',          # N,N'-dimethylurea
    'CNC(=O)OC',          # methyl N-methylcarbamate
    'CC(=S)NC',           # N-methylthioacetamide
    'NC=O',               # formamide
    'CC(=O)CC', 'CC=O', 'CC(=O)O', 'CC(=O)OC',
]

_IMINES = ['CC=NO', 'CCC(C)=NO', 'c1ccccc1C=NO', 'CC=NC', 'CC=NN']

_RINGS = ['C1=CCCCC1', 'C1=CCCCCC1', 'C1=CCCCCCC1', 'C1=CCCCCCCC1', 'C1CCCCC1', 'C1=CC=CC1']

_TETRAHEDRAL = [
    'CC(N)C(=O)O',        # alanine
    'CC(O)CC',            # butan-2-ol
    'CC(N)Cc1ccccc1',     # amphetamine skeleton
    'OCC(O)C(O)C(O)C(O)CO',   # a hexitol
    'CC(C)C(N)C(=O)O',    # valine
    'C1CCC(O)CC1', 'CC1CCC(C)CC1', 'OC1CCCCC1O', 'CC(Cl)Br', 'CC(O)C(N)C',
]

_FUSED = ['C1CC2CCC1CC2', 'C1CC2CCCC(C1)C2', 'C12CCCCC1CCCC2']

PRESERVE = (_ALLENES_AND_CUMULENES + _SPIRO + _AMIDES_AND_CARBONYLS + _IMINES + _RINGS
            + _TETRAHEDRAL + _FUSED)


# Molecules where V3 is deliberately different.  Agreement here would be the failure.  Each entry
# carries the reason, because a divergence without a stated reason is indistinguishable from a
# regression.
DIVERGENT = [
    ('CCCC[N+](C)(CC)CCC',
     'quaternary ammonium: V2 perceives no centre at all, its candidate pass opening with '
     '`atom == C`.  A cationic nitrogen has no lone pair to invert through, so the centre is '
     'configurationally stable and genuinely resolvable.'),
    ('C[S](=O)c1ccccc1',
     'sulfoxide: V2 perceives no centre -- not carbon, and it has a double bond, so V2 refuses it '
     'twice over.  Sulfoxide chirality is ordinary published chemistry.'),
    ('CC[S](C)=NC',
     'sulfilimine: V2 perceives no tetrahedral centre and reads the S=N as a CIS/TRANS axis -- both '
     'S and N form double bonds and the order is 2, so its cumulene walk accepts them.  The sulfur '
     'is a tetrahedral centre bearing a lone pair, not a planar cumulene terminal, so V3 reports one '
     'tetrahedral unit and no planar one.'),
]


# ------------------------------------------------------------------------------------------------

@lru_cache(maxsize=1)
def _counts():
    """`{smiles: [tetrahedral, axial, planar]}` from the pinned chython 2.

    One subprocess for the whole corpus, cached, because the interpreter start-up dominates and
    there is nothing per-test about it.

    The isolation, the version pin and the not-this-tree check all live in `oracle.py`, and `ask` runs
    `verify` before it answers -- so a leaked or mis-pinned oracle fails here on the first call rather
    than being discovered by a reader wondering why every comparison agrees.
    """
    wanted = PRESERVE + [s for s, _ in DIVERGENT]
    return oracle.ask('''
from chython import smiles
out = {}
for s in _payload:
    m = smiles(s)
    out[s] = [len(m.chiral_tetrahedrons), len(m.chiral_allenes), len(m.chiral_cis_trans)]
_emit(out)
''', wanted)


def _v3_counts(smi):
    """The same triple from the arena: (tetrahedral, axial, planar) stereogenic units."""
    units = read_smiles(smi).stereogenic_units()
    return [sum(1 for u in units if u['kind'] == kind)
            for kind in (SU_TETRA, SU_ALLENE, SU_CIS_TRANS)]


@requires_oracle
def test_this_files_oracle_is_pinned_and_is_not_the_tree_under_test():
    """BEFORE ANY COMPARISON IS BELIEVED.  A differential against yourself passes and means nothing.

    Both guards live in `oracle.verify` and every `oracle.ask` calls it, so this test is not the only
    thing standing between a leaked oracle and a green file.  It stays because the guards are worth
    naming in the file whose whole content depends on them, and because a reader who sees this fail
    knows to go and read `oracle.py`.  Their own tests, including the negative control that drops `-I`
    and watches the leak happen, are in `test_oracle.py`.
    """
    oracle.verify()


@requires_oracle
@pytest.mark.parametrize('smi', PRESERVE)
def test_v3_perceives_the_same_stereocentres_as_chython_2(smi):
    """Counts of tetrahedral, axial and planar stereocentres must match.

    Where this fails, one of two things is true and the test cannot tell which: V3 lost something V2
    reports, or the divergence is the intended one and the entry belongs in DIVERGENT with a written
    reason.  Decide it by chemistry, not by whichever side is more convenient.
    """
    counts = _counts()
    assert _v3_counts(smi) == counts[smi], (
        f'{smi}: V3 (tetra, axial, planar) = {_v3_counts(smi)}, chython 2 = {counts[smi]}')


@requires_oracle
@pytest.mark.parametrize('smi,reason', DIVERGENT, ids=[s for s, _ in DIVERGENT])
def test_v3_diverges_from_chython_2_where_chython_2_is_wrong(smi, reason):
    """These must DISAGREE.  Agreement means the gap came back."""
    counts = _counts()
    assert _v3_counts(smi) != counts[smi], (
        f'{smi}: V3 now agrees with chython 2 ({counts[smi]}), but it should not.  {reason}')


@requires_oracle
def test_the_sulfone_negative_control_agrees():
    """The control for the sulfimide and sulfoxide divergences.

    Both of those rest on sulfur being granted exactly ONE lone-pair direction.  A sulfone has no
    lone pair left, so it is not a centre -- and here V2 and V3 must AGREE.  Without this, a bug
    that handed sulfur a spurious extra direction would show up only as the divergences above
    "still diverging", which is what the test for them checks.
    """
    # ISOLATION MATTERS MOST HERE: the control's whole job is to AGREE, so a child that imported the
    # tree under test would agree perfectly and look like the comparison's strongest result.  Routed
    # through `oracle.ask`, the isolation is not this file's to remember.
    v2 = oracle.ask('''
from chython import smiles
m = smiles('CC[S](=O)(=O)C')
_emit([len(m.chiral_tetrahedrons), len(m.chiral_allenes), len(m.chiral_cis_trans)])
''')
    assert _v3_counts('CC[S](=O)(=O)C') == v2 == [0, 0, 0]


@requires_oracle
def test_the_cumulene_ladder_agrees_rung_by_rung_on_which_kind_it_is():
    """The single most important row of the comparison, stated on its own.

    The odd/even split is the whole of the allene-versus-cis/trans decision, and it is the thing a
    reimplementation loses first by special-casing three-atom allenes.  V2 classifies an arbitrary
    chain length, so agreement rung by rung is real evidence.
    """
    counts = _counts()
    ladder = ['CC=CC', 'CC=C=CC', 'CC=C=C=CC', 'CC=C=C=C=CC', 'CC=C=C=C=C=CC']
    for smi in ladder:
        v3, v2 = _v3_counts(smi), counts[smi]
        assert v3 == v2, f'{smi}: V3 {v3} vs chython 2 {v2}'
        # and the classification really does alternate, so that agreeing on all-zeroes
        # (a writer that had stopped perceiving anything) cannot pass this
        assert v3[1] or v3[2], f'{smi}: neither axial nor planar -- the rung vanished'
