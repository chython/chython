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
"""`tables/qed_alerts.tsv`: the structural alerts QED counts as its eighth input.

The shipped count is a RATCHET, not a claim of completeness -- see the table header and
`test_the_shipped_alert_count_is_a_ratchet`.
"""
from re import finditer

from chython.chemistry._tables import qed_alerts
from chython.core import read_smiles


#: The count this table ships today.  Raise it in the same commit that adds a row.
SHIPPED_ALERTS = 64


def implicit_junctions(pattern: str) -> list:
    """The `][` junctions in a pattern, i.e. every place a bond token was left out."""
    return [m.group('between') for m in finditer(r'\](?P<between>[()%\d]*)\[', pattern)]


def test_the_shipped_alert_count_is_a_ratchet():
    """The published list is longer than this table and the table says so.

    Brenk, Schipani, James, Krasowski, Gilbert, Frearson, Wyatt, ChemMedChem 2008, 3, 435 names 116
    alerts; this table transcribes the ones chython states, each with a probe that proves it fires.
    The assertion is an EQUALITY on the shipped number and not `>= 24` or `<= 116`: an inequality
    passes the day a row is dropped, which is the one change worth noticing.  `chython.chemistry.qed`
    documents that its ALERTS term is this table's count.
    """
    assert len(qed_alerts()) == SHIPPED_ALERTS


def test_every_alert_is_well_formed_and_uniquely_named():
    rows = qed_alerts()
    assert len({r.name for r in rows}) == len(rows), 'two alerts share a name'
    assert [r.id for r in rows] == [f'qed_alerts:{n}' for n in range(1, len(rows) + 1)]
    for r in rows:
        assert r.name and r.name != '-', r.id
        assert r.description and r.description != '-', r.id
        assert r.probe and r.probe != '-', r.id
        read_smiles(r.probe)                 # the probe is a SMILES and must parse


def test_every_alert_matches_its_own_probe():
    """An alert that matches nothing lowers no score and would sit here undetected."""
    for r in qed_alerts():
        assert r.query.is_substructure(read_smiles(r.probe)), f'{r.id} ({r.name}) misses its own probe'


def test_no_alert_matches_a_bare_alkane():
    """A runaway pattern is the other failure mode, and hexane must trip NOTHING.

    No escape hatch: the assertion is `== []`, never a disjunction over both possible answers.  Hexane
    has six chained sp3 carbons and the long-chain alert demands seven, so the answer is not in doubt.
    If a row is added that legitimately fires on hexane, change this test to name it.
    """
    tripped = [r.name for r in qed_alerts() if r.query.is_substructure(read_smiles('CCCCCC'))]
    assert tripped == [], tripped


def test_the_bare_alkane_gate_can_fail():
    """Positive control: the long-chain alert fires on heptane, and it is the only one that does.

    Without this, the sibling test passes the day the loader returns an empty tuple or every pattern
    stops matching -- the same silent green the probe test exists to prevent.  One carbon separates the
    two molecules, so the pair also pins the chain length the alert states.
    """
    tripped = [r.name for r in qed_alerts() if r.query.is_substructure(read_smiles('CCCCCCC'))]
    assert tripped == ['aliphatic long chain'], tripped


def test_no_two_alerts_share_a_pattern():
    patterns = [r.pattern for r in qed_alerts()]
    assert len(patterns) == len(set(patterns))


def test_no_alert_leaves_a_bond_implicit():
    """An implicit bond matches SINGLE ONLY, so `[C;*;z4][C;*;z4]` matches no aromatic ring at all.

    This is the dialect's most common transcription slip and the reason a row would miss its own probe
    for no visible reason, so it is refused outright rather than diagnosed twice.
    """
    for r in qed_alerts():
        assert not implicit_junctions(r.pattern), f'{r.id} ({r.name}) leaves a bond implicit'
