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
"""CIP cost: rings in series are linear (`cip_across`), a fused system past CIP_WORK_MAX refuses.

| molecule | what it pins |
| --- | --- |
| `_chain(k)`: CH(OH) linkers between k 1,4-cyclohexylenes | two laps per ring: 2^k walks and frontier pairs |
| perhydro-circumcoronene, fusion CH signed | simple paths exponential in the cyclomatic number: it refuses |
"""
from time import perf_counter
from chython.core import read_smiles


def _chain(k):
    """CH3-CH(OH)-[C6H10-CH(OH)]k-CH3, the linkers' signs alternating."""
    return 'C[C@H](O)' + ''.join(f'C{1 + i % 2}CCC(CC{1 + i % 2})[C@{"@" * (i % 2)}H](O)' for i in range(k)) + 'C'


# perhydro-circumcoronene: 19 fused cyclohexanes, the 36 fusion carbons signed in alternation
FLAKE = ('C1C[C@H]2C[C@@H]3CC[C@@H]4C[C@H]5CC[C@H]6C[C@H]7CC[C@H]8C[C@H]9CC[C@H]%10C[C@@H]%11CC[C@@H]%12'
         'C[C@@H]1[C@@H]1[C@@H]2[C@@H]2[C@H]3[C@H]4[C@@H]3[C@H]5[C@H]6[C@@H]4[C@H]7[C@@H]8[C@@H]5[C@H]9'
         '[C@H]%10[C@@H]6[C@H]%11[C@@H]%12[C@@H]1[C@H]1[C@H]2[C@H]3[C@H]4[C@@H]5[C@@H]61')


def test_chain_letters_are_the_unshared_walk_letters():
    """The letters a walk of every lap gives: sharing a branch across a bridge changes the cost, never a letter."""
    mol = read_smiles(_chain(12))
    mol.assign_cip()
    assert mol.atom_cips() == {2: 'S', 10: 'S', 18: 'R', 26: 'S', 34: 'R', 42: 'S', 58: 'R', 66: 'S', 74: 'R',
                               82: 'S', 90: 'R', 98: 'S'}


def test_a_long_chain_of_rings_decides_every_linker_quickly():
    """40 rings in series: 2^40 laps without the bridge memo and the frontier dedup.  The chain is
    constitutionally symmetric, so its middle linker is a true rule-1/2 tie and the only refusal."""
    mol = read_smiles(_chain(40))
    start = perf_counter()
    mol.assign_cip()
    assert perf_counter() - start < 5
    assert len(mol.atom_cips()) == 40
    refused = [r for r in mol.log.by_stage('cip') if r.rule == 'cip:undecided']
    assert [r.message for r in refused] == ['atom 162: rules 1a/1b/2 do not settle the site']


def test_a_fused_system_past_the_budget_refuses_every_site_and_names_the_budget():
    mol = read_smiles(FLAKE)
    start = perf_counter()
    assert mol.assign_cip() is False
    assert perf_counter() - start < 5
    assert mol.atom_cips() == {}
    refused = [r for r in mol.log.by_stage('cip') if r.rule == 'cip:undecided']
    assert len(refused) == 36
    assert all('spent more than' in r.message for r in refused)


def test_an_inverted_ranked_letter_swaps():
    """`canonicalize_stereo_groups` inverts the AND group's phase, and rules 1a/1b/2 rank the site."""
    mol = read_smiles(_chain(12).replace('C[C@H](O)', 'C[C@@H](O)', 1) + ' |&1:1|')
    mol.assign_cip()
    assert mol.atom_cips()[2] == 'R'
    assert mol.canonicalize_stereo_groups() is True
    assert mol.atom_cips()[2] == 'S'


def test_an_inverted_letter_past_the_budget_is_dropped_not_kept():
    """`_carry_inversion` ranks before it swaps; a ranking refused by the budget drops the stated letter."""
    mol = read_smiles(FLAKE + ' |&1:2|')
    mol.set_atom_cip(3, 'R')
    start = perf_counter()
    assert mol.canonicalize_stereo_groups() is True
    assert perf_counter() - start < 5
    assert 3 not in mol.atom_cips()
    assert any(r.rule == 'stereo-groups:cip-dropped' for r in mol.log.by_stage('stereo-groups'))
