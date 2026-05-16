# -*- coding: utf-8 -*-
#
#  Copyright 2026 Tagir Akhmetshin <tagirshin@gmail.com>
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
"""Regression tests for K-component query matching against
<K-component targets.

The bug: ``Isomorphism._get_mapping`` enumerated target components via
``permutations(target.connected_components, len(query.connected_components))``
which is empty whenever ``target`` has fewer components than the query.
This silently dropped every mapping for the common case of "disconnected
query, single connected target".
"""
from chython import smiles, smarts
from chython.containers.query import QueryContainer
from chython.periodictable import QueryElement


def _two_disjoint_carbons() -> QueryContainer:
    q = QueryContainer("")
    q.add_atom(QueryElement.from_atomic_number(6)(), 1)
    q.add_atom(QueryElement.from_atomic_number(6)(), 2)  # separate component
    return q


def test_disjoint_query_into_single_component_target():
    # Two-carbon disconnected query vs. ethane (CC): the algorithm must
    # find at least one mapping. (Both ordered mappings have the same
    # target atom-set, so the default ``automorphism_filter=True`` keeps
    # only one — but the *bug* was zero mappings.)
    q = _two_disjoint_carbons()
    mappings = list(q.get_mapping(smiles("CC")))
    assert len(mappings) >= 1, (
        f"expected at least 1 mapping into ethane, got {len(mappings)}"
    )
    mappings = list(q.get_mapping(smiles("CC"), automorphism_filter=False))
    assert len(mappings) >= 2, (
        f"with automorphism_filter=False expected both ordered mappings, "
        f"got {len(mappings)}"
    )


def test_disjoint_query_into_single_component_target_atom_disjoint():
    # Returned mappings must never reuse a target atom for two query atoms.
    q = _two_disjoint_carbons()
    for m in q.get_mapping(smiles("CC")):
        target_atoms = list(m.values())
        assert len(target_atoms) == len(set(target_atoms)), (
            f"mapping {m} reuses a target atom — atom-disjointness violated"
        )


def test_disjoint_query_uses_all_query_components():
    # Each mapping must cover every atom of the query (i.e. no partial
    # mappings emitted from the multi-component branch).
    q = _two_disjoint_carbons()
    query_atoms = set(q._atoms)
    for m in q.get_mapping(smiles("CC")):
        assert set(m.keys()) == query_atoms


def test_kk_case_unchanged():
    # Two-component query against a two-component target must continue
    # to work — the new branch must subsume the old behaviour without
    # regression. Note: the new branch additionally yields *intra-
    # component* mappings (two query carbons inside the same CC), so
    # the total count is strictly larger than before; we just assert it
    # is at least the old count and that disjointness still holds.
    q = _two_disjoint_carbons()
    mappings = list(q.get_mapping(smiles("CC.CC")))
    # Old behaviour: 4 cross-component mappings (after automorphism filter).
    # New behaviour: 4 cross + 2 intra (same CC matched twice atom-disjoint)
    # = 6. We require the old set to still be present.
    assert len(mappings) >= 4
    for m in mappings:
        vs = list(m.values())
        assert len(vs) == len(set(vs))


def test_disjoint_query_into_larger_target_finds_all_pairs():
    # Two disconnected query carbons against propane (CCC): there are
    # 3*2 ordered choices of distinct target atoms — the matcher must
    # find each as a valid mapping when automorphism filtering is off.
    q = _two_disjoint_carbons()
    mappings = list(q.get_mapping(smiles("CCC"), automorphism_filter=False))
    distinct_ordered = {tuple(sorted(m.items())) for m in mappings}
    assert len(distinct_ordered) >= 6, (
        f"expected to find all 6 ordered distinct pairs, "
        f"got {len(distinct_ordered)}"
    )


def test_multicomponent_smarts_against_single_component_target():
    # A SMARTS query with two disconnected fragments (joined by '.') must
    # match a single-component target — this is the chython.smarts path
    # SynPlanner actually uses.
    q = smarts("[C].[O]")
    target = smiles("CO")
    mappings = list(q.get_mapping(target))
    assert len(mappings) >= 1, (
        f"disconnected SMARTS [C].[O] must match CO, got {len(mappings)} mappings"
    )
    for m in mappings:
        vs = list(m.values())
        assert len(vs) == len(set(vs))


def test_three_component_query_into_two_component_target():
    # K=3 query vs K=2 target — the algorithm must allow two of the
    # query components to share a target component as long as atoms remain
    # disjoint.
    q = QueryContainer("")
    for n in (1, 2, 3):
        q.add_atom(QueryElement.from_atomic_number(6)(), n)
    mappings = list(q.get_mapping(smiles("CC.CC")))
    # at least one mapping should exist (4 atoms available, 3 needed,
    # atom-disjoint).
    assert len(mappings) >= 1
    for m in mappings:
        vs = list(m.values())
        assert len(vs) == len(set(vs))


def test_d_constrained_bonded_query_into_single_component():
    # Mirrors the SP rule shape `[C;D1]-[C;D3].[C;D3]-[N;D1]`: two bonded
    # fragments with D-degree constraints. Target has two D3 carbons each
    # carrying their own D1 substituents so an atom-disjoint pairing exists.
    q = smarts('[C;D1]-[C;D3].[C;D3]-[N;D1]')
    mappings = list(q.get_mapping(smiles('CC(N)C(N)C')))
    assert len(mappings) >= 2, (
        f"expected ≥2 disjoint mappings for [C;D1]-[C;D3].[C;D3]-[N;D1] "
        f"into CC(N)C(N)C, got {len(mappings)}"
    )
    for m in mappings:
        vs = list(m.values())
        assert len(vs) == len(set(vs))


def test_d_constrained_query_rejects_overlap_only_target():
    # Same rule shape, but target has only one D3 carbon — both query
    # fragments would have to share it. Atom-disjoint enforcement must
    # return zero mappings.
    q = smarts('[C;D1]-[C;D3].[C;D3]-[N;D1]')
    assert list(q.get_mapping(smiles('CC(N)C'))) == []
