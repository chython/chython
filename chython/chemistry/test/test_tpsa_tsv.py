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
from pytest import mark
from chython.chemistry._tables import TPSA_CLASSES, first_match, tpsa_rules
from chython.core import read_smarts, read_smiles


def test_the_table_has_ertls_forty_three_environments():
    # Ertl, Rohde, Selzer, J. Med. Chem. 2000, 43, 3714, Table 1: 26 nitrogen, 6 oxygen, 7 sulfur,
    # 4 phosphorus.  The class counts are the same table counted a second way (26 + 6 == 32 NO,
    # 7 + 4 == 11 SP), so a misfiled row fails the class check and a mistranscribed one the element.
    rows = tpsa_rules()
    assert len(rows) == 43
    counts = {}
    for r in rows:
        element = r.description.split()[0].rstrip(',')   # 'nitrogen cation, NH3' -> 'nitrogen'
        counts[element] = counts.get(element, 0) + 1
    assert counts == {'nitrogen': 26, 'oxygen': 6, 'sulfur': 7, 'phosphorus': 4}
    assert sum(1 for r in rows if r.element_class == 'NO') == 32
    assert sum(1 for r in rows if r.element_class == 'SP') == 11


def test_every_row_is_well_formed():
    for r in tpsa_rules():
        assert r.id.startswith('tpsa:')
        assert r.element_class in TPSA_CLASSES
        assert r.contribution >= 0.0
        # `map_numbers()` returns non-zero entries only, so this asks whether the row wrote `:1`.
        # Do not use `compile_smarts`'s `numbers` dict instead: it auto-numbers every unnumbered
        # atom from the lowest unclaimed integer, so `1 in numbers` can never be false.
        assert 1 in set(r.query.map_numbers().values()), r.id
        assert r.description and r.description != '-'


def test_the_subject_check_can_fail():
    """Negative control for the line above, which is unfailable if written the other way."""
    assert 1 not in set(read_smarts('[O;D2;z1;h0]').map_numbers().values())
    assert 1 in set(read_smarts('[O;D2;z1;h0:1]').map_numbers().values())


def test_the_ids_are_the_file_order():
    assert [r.id for r in tpsa_rules()] == [f'tpsa:{n}' for n in range(1, 44)]


def test_the_ring_rows_precede_their_generic_rows():
    # first-match-wins: an epoxide oxygen also matches the ether row, and an aziridine nitrogen the
    # generic sp3 amine row.  a table sorted the other way is a different descriptor.
    order = {r.pattern: n for n, r in enumerate(tpsa_rules())}
    assert order['[O;D2;z1;h0;r3:1]'] < order['[O;D2;z1;h0:1]']
    assert order['[N;D3;z1;h0;r3:1]'] < order['[N;D3;z1;h0:1]']
    assert order['[N;D2;z1;h1;r3:1]'] < order['[N;D2;z1;h1:1]']


def test_every_contribution_is_the_published_number():
    got = {r.pattern: r.contribution for r in tpsa_rules()}
    assert got['[N;D3;z1;h0:1]'] == 3.24
    assert got['[N;D1;z1;h2:1]'] == 26.02
    assert got['[N;D2;z4;h0:1]'] == 12.89
    assert got['[N;D2;z4;h1:1]'] == 15.79
    assert got['[N;D4;z1;h0;+:1]'] == 0.00
    assert got['[O;D2;z1;h0:1]'] == 9.23
    assert got['[O;D1;z2;h0:1]'] == 17.07
    assert got['[O;D1;z1;h1:1]'] == 20.23
    assert got['[O;D1;z1;h0;-:1]'] == 23.06
    assert got['[S;D4;z5;h0:1]'] == 8.38
    assert got['[P;D4;z2;h0:1]'] == 9.81


def test_no_two_rows_share_a_pattern():
    patterns = [r.pattern for r in tpsa_rules()]
    assert len(patterns) == len(set(patterns))


# Reachability gate (ruling F3-31).  In a first-match table a row can match its probe and still be
# dead, because an earlier row claimed the same atom; so the gate resolves each probe through the
# whole table in file order and asserts the row under test is the one that won.
PROBES = {
    1: 'C1CN1C',      2: 'C1CN1',        3: 'CN(C)C',        4: 'CC=NC',
    5: 'CC#N',        6: 'CN(=O)=O',     7: 'CC=N#N',        8: 'CNC',
    9: 'CC=N',       10: 'CN',          11: 'C[N+](C)(C)C', 12: 'CC=[N+](C)C',
   13: 'C[N+]#[C-]', 14: 'C[NH+](C)C',  15: 'CC=[NH+]C',    16: 'C[NH2+]C',
   17: 'CC=[NH2+]',  18: 'C[NH3+]',     19: 'c1ccncc1',     20: 'c1ccn2cccc2c1',
   21: 'Cn1cccc1',   22: 'O=n1ccccc1',  23: 'c1cc[nH]c1',   24: 'c1ccc2cccc[n+]2c1',
   25: 'C[n+]1ccccc1', 26: 'c1cc[nH+]cc1', 27: 'C1CO1',     28: 'COC',
   29: 'CC=O',       30: 'CCO',         31: 'CC(=O)[O-]',   32: 'c1ccoc1',
   33: 'CSC',        34: 'CC=S',        35: 'CS(=O)C',      36: 'CS(=O)(=O)C',
   37: 'CS',         38: 'c1ccsc1',     39: 'O=s1cccc1',    40: 'CP(C)C',
   41: 'CP=C',       42: 'COP(=O)(OC)OC', 43: 'COP(=O)OC',
}


def _resolve(molecule):
    """{stable id: rule id}, first match wins in file order -- the descriptor's own resolution.

    Deliberately calls `first_match` rather than re-implementing the loop, so the 43 probes below
    are also its coverage.  Re-implementing here leaves `first_match` untested.
    """
    return {atom: row.id for atom, row in first_match(tpsa_rules(), molecule).items()}


@mark.parametrize('n', list(PROBES))
def test_every_row_wins_its_own_probe(n):
    smi = PROBES[n]
    assert f'tpsa:{n}' in set(_resolve(read_smiles(smi)).values()), \
        f'row {n} never wins on {smi}: it is dead or shadowed by an earlier row'


def test_the_probe_table_covers_every_row():
    # without this, deleting a probe silently retires its reachability check.
    assert {f'tpsa:{n}' for n in PROBES} == {r.id for r in tpsa_rules()}


def test_the_reachability_gate_can_fail():
    # Negative control (ruling F3-32): a row appended after the generic ether row can never win on
    # an ether oxygen, which is the shadowing the gate above exists to detect.
    m = read_smiles('COC')
    assert set(_resolve(m).values()) == {'tpsa:28'}          # the ether row, and only it
    assert 'tpsa:27' not in set(_resolve(m).values())        # the epoxide row is not reached here
