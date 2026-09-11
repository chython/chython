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
from chython.core import read_smiles
from chython.chemistry._tables import (CRIPPEN_CATCH_ALLS, CRIPPEN_ROLES, CRIPPEN_TYPES,
                                       crippen_rules, crippen_rules_by_role, first_match)


NO_PUBLISHED_MR = frozenset(('N10', 'N12', 'O12', 'Hal', 'Me2'))


def test_the_inventory_is_exactly_the_papers():
    # Wildman, Crippen, J. Chem. Inf. Comput. Sci. 1999, 39, 868, Table 1.
    rows = crippen_rules()
    # 72 published types in 73 rows: S2 needs two, its published alternatives being a bare charge test
    # and an S=X bond test, which chython cannot join without recursive SMARTS.  Both the length and the
    # set are asserted, since collapsing one into the other hides a duplicate or a missing type.
    assert len(rows) == len(CRIPPEN_TYPES) == 73
    assert [r.type for r in rows] == list(CRIPPEN_TYPES)      # file order IS inventory order
    assert len({r.type for r in rows}) == 72
    assert [r.type for r in rows].count('S2') == 2            # the one type spelled twice
    assert all(t == 'S2' for t in {t for t in CRIPPEN_TYPES
                                   if [r.type for r in rows].count(t) > 1})


def test_every_row_is_well_formed():
    for r in crippen_rules():
        assert r.id == f'crippen:{r.type}'
        assert r.role in CRIPPEN_ROLES
        assert isinstance(r.logp, float)
        assert isinstance(r.mr, float) and r.mr >= 0.0
        # asks whether the row wrote `:1`, the atom `first_match` types.  Not `r.numbers`, which
        # auto-numbers every atom and so cannot fail.
        assert 1 in set(r.query.map_numbers().values()), r.type
        assert r.description and r.description != '-'


def test_the_five_types_with_no_published_mr_say_so():
    # the MR cell is blank for exactly these five; the count is asserted too, so the test's name cannot
    # drift from its data
    assert len(NO_PUBLISHED_MR) == 5
    got = {r.type for r in crippen_rules() if not r.mr_published}
    assert got == NO_PUBLISHED_MR
    for r in crippen_rules():
        if not r.mr_published:
            assert r.mr == 0.0


def test_a_published_zero_mr_is_not_an_absent_one():
    # O7 and O9 carry a real published MR of exactly 0; without `mr_published` they would be
    # indistinguishable from the five blanks
    got = {r.type: (r.mr, r.mr_published) for r in crippen_rules()}
    assert got['O7'] == (0.0, True)
    assert got['O9'] == (0.0, True)
    assert got['O12'] == (0.0, False)


def test_every_type_wins_its_own_probe():
    """No row is shadowed by a row above it.

    A first-match table's order disambiguates deliberately overlapping types, so a row that is a
    special case of an earlier one never fires and its published contribution silently never reaches an
    answer.  Scoped by role: a hydrogen row's subject is the carrier atom, which a heavy row would
    claim first in an unscoped pass.
    """
    by_role = crippen_rules_by_role()
    for r in crippen_rules():
        if r.type in CRIPPEN_CATCH_ALLS:
            assert r.probe == '-', f'{r.type} is a catch-all and must not carry a probe'
            continue
        winners = first_match(by_role[r.role], read_smiles(r.probe))
        assert r in winners.values(), (
            f'{r.type} wins no atom of its own probe {r.probe!r}; it is shadowed by '
            f'{sorted({w.type for w in winners.values()})} and must move above them'
        )


def test_the_probe_gate_can_fail():
    """Negative control for `test_every_type_wins_its_own_probe`, whose false case is the paper's own
    numbering: with C1 ahead of C8, C8 wins nothing on its own probe, which is why the file order moves
    C8 up.  Both orders are resolved through the same `first_match` the gate uses.
    """
    rows = {r.type: r for r in crippen_rules()}
    c1, c8 = rows['C1'], rows['C8']
    probe = read_smiles(c8.probe)
    assert c8 in first_match([c8, c1], probe).values()        # the file's order: C8 is reachable
    # by type name, not by row: a CrippenRow carries an unhashable QueryContainer, so a set of rows
    # raises TypeError instead of failing the assertion.
    winners = {r.type for r in first_match([c1, c8], probe).values()}   # paper's order: C1 swallows it
    assert 'C8' not in winners and 'C1' in winners


def test_the_measured_shadowing_pairs_stay_fixed():
    """The seven measured orderings, pinned so a later re-sort by type name cannot undo them.

    Four are from our widening to one row per type: C1 first would swallow C8 (toluene's methyl) and C2
    (neopentane's quaternary carbon), H2's `[O;H1,H2;!z4]` subsumes H4's carboxyl OH, C6 subsumes C26.
    Two are the paper's own, marked "order flip here is intentional" in the transcription source: O7
    swallows O12's carboxylate and S1 swallows S2's sulfoxide sulfur.  The seventh runs the other way:
    broad O9 sits below O10 and O11, which reproduces all five of its published alternatives.
    """
    order = {r.type: n for n, r in enumerate(crippen_rules())}
    assert order['C8'] < order['C1']
    assert order['C2'] < order['C1']
    assert order['C26'] < order['C6']
    assert order['H4'] < order['H2']
    assert order['O12'] < order['O7']
    assert order['S2'] < order['S1']
    assert order['O10'] < order['O9'] and order['O11'] < order['O9']


def test_each_catch_all_is_last_in_its_block():
    order = {r.type: n for n, r in enumerate(crippen_rules())}
    assert order['CS'] > max(order[f'C{n}'] for n in range(1, 28))
    assert order['HS'] > max(order[f'H{n}'] for n in range(1, 5))
    assert order['NS'] > max(order[f'N{n}'] for n in range(1, 15))
    assert order['OS'] > max(order[f'O{n}'] for n in range(1, 13))


def test_the_hydrogen_rows_are_the_hydrogen_role_and_nothing_else():
    by_role = {}
    for r in crippen_rules():
        by_role.setdefault(r.role, []).append(r.type)
    assert sorted(by_role['hydrogen']) == ['H1', 'H2', 'H3', 'H4', 'HS']
    assert len(by_role['heavy']) == 68            # 67 heavy types, S2 spelled in two rows


def test_the_spot_check_values_are_the_published_ones():
    got = {r.type: (r.logp, r.mr) for r in crippen_rules()}
    assert got['C1'] == (0.1441, 2.503)
    assert got['C2'] == (0.0000, 2.433)
    assert got['C18'] == (0.1581, 3.350)
    assert got['H1'] == (0.1230, 1.057)
    assert got['H2'] == (-0.2677, 1.395)
    assert got['H3'] == (0.2142, 0.9627)
    assert got['H4'] == (0.2980, 1.805)
    assert got['HS'] == (0.1125, 1.112)
    assert got['N1'] == (-1.0190, 2.262)
    assert got['O1'] == (0.1552, 1.0800)     # O1 is the aromatic oxygen, not the alcohol
    assert got['O2'] == (-0.2893, 0.8238)    # the alcohol/water oxygen
    assert got['F'] == (0.4202, 1.108)
    assert got['Cl'] == (0.6895, 5.853)
    assert got['Br'] == (0.8456, 8.927)
    assert got['I'] == (0.8857, 14.02)
    assert got['P'] == (0.8612, 6.920)
    assert got['S1'] == (0.6482, 7.591)
    assert got['Me1'] == (-0.3808, 5.754)


def test_no_row_carries_a_transcription_placeholder():
    # a plan may not invent a published number; a table may not ship without one.
    for r in crippen_rules():
        assert (r.logp, r.mr) != (0.0, 0.0) or r.type in NO_PUBLISHED_MR
