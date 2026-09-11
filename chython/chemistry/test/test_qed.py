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
"""QED: eight desirability functions, three weight sets and the weighted geometric mean."""
from pytest import approx, mark, raises

from chython.chemistry import qed
from chython.chemistry._maccs import aromatic_ring_count
from chython.chemistry._qed import (ADS_PARAMETERS, QED_PROPERTIES, QED_WEIGHTS, ads, alert_count,
                                    qed_properties)
from chython.core import read_smiles


def test_every_property_has_parameters_and_three_weight_sets():
    assert len(QED_PROPERTIES) == 8
    assert set(ADS_PARAMETERS) == set(QED_PROPERTIES)
    assert set(QED_WEIGHTS) == {'mean', 'max', 'unweighted'}
    for name, w in QED_WEIGHTS.items():
        assert set(w) == set(QED_PROPERTIES), name


@mark.parametrize('name', QED_PROPERTIES)
def test_every_ads_maximum_equals_its_published_dmax(name):
    """`dmax` IS the maximum of the un-normalised ADS, which makes it a free check on six coefficients.

    A mis-typed coefficient is otherwise invisible: the score stays in range and stays plausible.
    """
    a, b, c, d, e, f, dmax = ADS_PARAMETERS[name]
    lo, hi = {'MW': (0, 1000), 'ALOGP': (-10, 15), 'HBA': (0, 30), 'HBD': (0, 30),
              'PSA': (0, 400), 'ROTB': (0, 40), 'AROM': (0, 15), 'ALERTS': (0, 30)}[name]
    step = (hi - lo) / 200000
    best = max(ads(lo + n * step, a, b, c, d, e, f, 1.0) for n in range(200001))
    assert best == approx(dmax, rel=1e-6)


def test_the_normalised_ads_never_leaves_the_unit_interval():
    for name in QED_PROPERTIES:
        a, b, c, d, e, f, dmax = ADS_PARAMETERS[name]
        for x in (-50.0, 0.0, 1.0, 10.0, 100.0, 1000.0):
            v = ads(x, a, b, c, d, e, f, dmax)
            assert 0.0 <= v <= 1.0 + 1e-9, (name, x, v)


def test_the_eight_inputs_are_the_f3_quantities():
    """Every input is a quantity chython already answers, and QED derives none of its own."""
    m = read_smiles('CC(=O)Nc1ccc(O)cc1')                # paracetamol
    p = qed_properties(m)
    assert set(p) == set(QED_PROPERTIES)
    assert p['MW'] == approx(float(m))
    assert p['PSA'] == approx(m.tpsa)
    assert p['ALOGP'] == approx(m.crippen_logp)
    assert p['HBA'] == m.hydrogen_bond_acceptors_count == 2
    assert p['HBD'] == m.hydrogen_bond_donors_count == 2
    assert p['ROTB'] == m.rotatable_bonds_count == 1
    assert p['AROM'] == 1
    assert p['ALERTS'] == alert_count(m)


def test_the_score_is_in_the_unit_interval():
    for smi in ('CC(=O)Nc1ccc(O)cc1', 'CC(=O)Oc1ccccc1C(=O)O', 'CC(C)Cc1ccc(cc1)C(C)C(=O)O',
                'c1ccccc1', 'CCCCCCCCCCCCCCCCCC', 'C'):
        v = qed(read_smiles(smi))
        assert 0.0 <= v <= 1.0, (smi, v)


def test_a_drug_like_molecule_scores_above_a_long_alkane():
    assert qed(read_smiles('CC(=O)Nc1ccc(O)cc1')) > qed(read_smiles('CCCCCCCCCCCCCCCCCC'))


def test_the_three_weight_sets_give_three_scores_and_mean_is_the_default():
    m = read_smiles('CC(=O)Oc1ccccc1C(=O)O')
    assert qed(m) == approx(qed(m, weights='mean'))
    assert m.qed == approx(qed(m, weights='mean'))
    assert qed(m, weights='max') != approx(qed(m, weights='mean'))
    assert qed(m, weights='unweighted') != approx(qed(m, weights='mean'))


def test_an_unknown_weight_set_is_refused_by_name():
    with raises(ValueError, match='mean'):
        qed(read_smiles('CCO'), weights='lipinski')


def test_alert_count_counts_distinct_alerts_and_not_matches():
    """Two nitro groups are one alert kind, not two.

    The pair is the assertion: one nitro and two must give the same count and the same alert set, so
    the test fails if `alert_count` ever counts matches.  The absolute number is not asserted -- it is
    a property of `tables/qed_alerts.tsv`, which is ratcheted where it lives.
    """
    one = read_smiles('[O-][N+](=O)c1ccccc1')                    # nitrobenzene
    two = read_smiles('[O-][N+](=O)c1ccc([N+](=O)[O-])cc1')      # 1,4-dinitrobenzene
    assert alert_count(two) == alert_count(one) >= 1


def test_a_molecule_with_no_alert_scores_zero_alerts():
    assert alert_count(read_smiles('CC(C)Cc1ccc(cc1)C(C)C(=O)O')) == 0   # ibuprofen
    assert qed_properties(read_smiles('CC(C)Cc1ccc(cc1)C(C)C(=O)O'))['ALERTS'] == 0.0


def test_it_is_renumbering_invariant():
    assert qed(read_smiles('CC(=O)Nc1ccc(O)cc1')) == approx(qed(read_smiles('Oc1ccc(NC(C)=O)cc1')))


def test_the_arom_term_is_the_containers_own_aromatic_ring_count():
    """AROM is `MoleculeContainer.aromatic_rings_count`, reached through `_maccs.aromatic_ring_count`.

    Unconditional, not a `hasattr` guard: a conditional assertion is not an assertion.
    """
    for smi, want in (('c1ccc2ccccc2c1', 2), ('c1ccccc1', 1), ('C1CCCCC1', 0),
                      ('CC(=O)Nc1ccc(O)cc1', 1)):
        m = read_smiles(smi)
        assert aromatic_ring_count(m) == want == m.aromatic_rings_count, smi
        assert qed_properties(m)['AROM'] == float(want), smi
