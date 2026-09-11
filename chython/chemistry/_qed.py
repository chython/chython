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
"""Quantitative estimate of drug-likeness (QED), Bickerton, Paolini, Besnard, Muresan, Hopkins,
Nat. Chem. 2012, 4, 90.

A weighted geometric mean of eight asymmetric double sigmoidal desirability functions, one per
molecular property.  Seven of the eight inputs are chython's own published descriptors; the eighth is
the count of structural alerts, and that count is `tables/qed_alerts.tsv`'s -- 64 alerts, fewer than the
116 the alert set names, each proved to fire by its own probe.  So a chython score and another
implementation's score are two numbers with the same name, and neither reads the other.
"""
from math import exp, log as ln

from ._counts import (hydrogen_bond_acceptors_count, hydrogen_bond_donors_count,
                      rotatable_bonds_count)
from ._crippen import crippen_logp
from ._maccs import aromatic_ring_count
from ._tables import qed_alerts
from ._tpsa import tpsa


#: The eight properties QED weighs, in the paper's order.
QED_PROPERTIES = ('MW', 'ALOGP', 'HBA', 'HBD', 'PSA', 'ROTB', 'AROM', 'ALERTS')

#: Table 2 of the paper: the ADS coefficients `(a, b, c, d, e, f)` and the published `dmax`.
#: `dmax` is the maximum of the un-normalised function over the property's range, which is why
#: `test_qed.py` can check every coefficient by scanning: a typo moves the maximum.
ADS_PARAMETERS = {
    'MW':     (2.817065973, 392.5754953, 290.7489764, 2.419764353, 49.22325677, 65.37051707,
               104.9805805),
    'ALOGP':  (3.172690585, 137.8624751, 2.534937431, 4.581497897, 0.822739154, 0.576295591,
               131.3186604),
    'HBA':    (2.948620388, 160.4605972, 3.615294657, 4.435986202, 0.290141953, 1.300669958,
               148.7763046),
    'HBD':    (1.618662227, 1010.051101, 0.985094388, 0.000000001, 0.713820843, 0.920922555,
               258.1632616),
    'PSA':    (1.876861559, 125.2232657, 62.90773554, 87.83366614, 12.01999824, 28.51324732,
               104.5686167),
    'ROTB':   (0.010000000, 272.4121427, 2.558379970, 1.565547684, 1.271567166, 2.758063707,
               105.4420403),
    'AROM':   (3.217788970, 957.7374108, 2.274627939, 0.000000001, 1.317690384, 0.375760881,
               312.3372610),
    'ALERTS': (0.010000000, 1199.094025, -0.09002883, 0.000000001, 0.185904477, 0.875193782,
               417.7253140)}

#: The paper's three weight sets.  `'mean'` is QED_w,mo, its recommended default and `mol.qed`.
QED_WEIGHTS = {
    'mean': {'MW': 0.66, 'ALOGP': 0.46, 'HBA': 0.05, 'HBD': 0.61, 'PSA': 0.06, 'ROTB': 0.65,
             'AROM': 0.48, 'ALERTS': 0.95},
    'max': {'MW': 0.50, 'ALOGP': 0.25, 'HBA': 0.00, 'HBD': 0.50, 'PSA': 0.00, 'ROTB': 0.50,
            'AROM': 0.25, 'ALERTS': 1.00},
    'unweighted': {p: 1.0 for p in QED_PROPERTIES}}


def _sigmoid(t: float) -> float:
    """Logistic, clamped: `exp(-t)` overflows past |t| ~ 710 and the answer is 0 or 1 well before."""
    if t < -700.0:
        return 0.0
    elif t > 700.0:
        return 1.0
    return 1.0 / (1.0 + exp(-t))


def ads(x: float, a: float, b: float, c: float, d: float, e: float, f: float, dmax: float) -> float:
    """The asymmetric double sigmoidal desirability function, equation 1, normalised by `dmax`.

    Pass `dmax=1.0` for the un-normalised value, whose maximum over the property's range IS the
    published `dmax` -- which is how a mis-typed coefficient is caught, since the score alone stays
    plausible.
    """
    return (a + b * _sigmoid((x - c + d / 2.0) / e) * (1.0 - _sigmoid((x - c - d / 2.0) / f))) / dmax


def alert_count(molecule) -> int:
    """How many of `tables/qed_alerts.tsv`'s alerts this molecule contains.

    DISTINCT ALERTS, not matches: two nitro groups are one alert.  Alerts overlap by design -- a
    charge-separated nitro group is also an N-O single bond -- and both are counted, because the score
    weighs alert kinds and not sites.
    """
    return sum(1 for row in qed_alerts() if row.query.is_substructure(molecule))


def qed_properties(molecule) -> dict:
    """QED's eight raw inputs, before desirability.  Every one is a descriptor chython already answers.

    `MW` is `float(molecule)`, the average molecular mass; `AROM` is the container's own aromatic ring
    count reached through `_maccs.aromatic_ring_count`, so there is one answer to that question.
    """
    return {'MW': float(molecule),
            'ALOGP': crippen_logp(molecule),
            'HBA': float(hydrogen_bond_acceptors_count(molecule)),
            'HBD': float(hydrogen_bond_donors_count(molecule)),
            'PSA': tpsa(molecule),
            'ROTB': float(rotatable_bonds_count(molecule)),
            'AROM': float(aromatic_ring_count(molecule)),
            'ALERTS': float(alert_count(molecule))}


def qed(molecule, *, weights='mean') -> float:
    """Quantitative estimate of drug-likeness, in [0, 1].

    `weights='mean'` is the paper's QED_w,mo and the default, `'max'` is QED_w,max and `'unweighted'`
    is QED_w,u.  `mol.qed` is the default variant.  The ALERTS term counts `tables/qed_alerts.tsv`,
    which holds the alerts chython states rather than the whole published list, so this score is
    chython's and is not a reading of another implementation's.
    """
    if weights not in QED_WEIGHTS:
        raise ValueError(f'weights must be one of {tuple(QED_WEIGHTS)}, got {weights!r}')
    w = QED_WEIGHTS[weights]
    properties = qed_properties(molecule)
    weighted = total = 0.
    for name in QED_PROPERTIES:
        if not w[name]:
            continue                      # QED_w,max zeroes HBA and PSA, and ln() has no value there
        d = ads(properties[name], *ADS_PARAMETERS[name])
        if d <= 0.:
            return 0.                     # a geometric mean with a zero factor is zero
        weighted += w[name] * ln(d)
        total += w[name]
    return exp(weighted / total)
