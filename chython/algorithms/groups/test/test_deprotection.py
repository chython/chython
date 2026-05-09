# -*- coding: utf-8 -*-
#
#  Copyright 2022-2026 Ramil Nugmanov <nougmanoff@protonmail.com>
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
import pytest
from chython import smiles
from chython.algorithms.groups._protective import rules


_test_cases = [(name, *rule) for name, rule in rules.items()]


@pytest.mark.parametrize('name,q,keep,add,test_smi,cleaved_smi,decoys', _test_cases, ids=[x[0] for x in _test_cases])
def test_deprotection(name, q, keep, add, test_smi, cleaved_smi, decoys):
    t = smiles(test_smi)
    t.canonicalize()
    # query matches protected form
    assert q < t, f'query does not match protected form'
    # deprotection produces expected product
    a = smiles(cleaved_smi)
    a.canonicalize()
    t.remove_protection(name)
    assert t == a, f'deprotection gave {t}, expected {a}'
    # decoys are not matched
    for d in decoys:
        d = smiles(d)
        d.canonicalize()
        assert not q < d, f'query matches decoy {d}'


def test_no_false_matches_from_ordering():
    """Ensure general patterns come after specific ones to prevent false atom claims.

    If pattern B matches within the protected molecule of pattern A, then A must come
    BEFORE B in the rules dict (so A claims atoms first). This test fails when B comes
    before A - indicating rules dict needs reordering.
    """
    names = list(rules.keys())
    for name_a, (_, _, _, test_smi, *_) in rules.items():
        mol = smiles(test_smi)
        mol.canonicalize()
        idx_a = names.index(name_a)
        for name_b in names[:idx_a]:  # patterns that come BEFORE name_a
            q_b = rules[name_b][0]
            if next(q_b.get_mapping(mol), None) is not None:
                pytest.fail(
                    f"Pattern '{name_b}' (idx {names.index(name_b)}) matches within "
                    f"protected molecule of '{name_a}' (idx {idx_a}), but comes before it. "
                    f"Move '{name_b}' after '{name_a}' in _protective.py rules dict."
                )


def test_protective_groups_no_overlap():
    """Ensure protective_groups doesn't report false positives due to sub-pattern overlap.

    Checks that deletable atoms (atoms not in keep list) of different detected PGs
    don't overlap - matching the logic in the protective_groups property.
    """
    for name, (q, keep, add, test_smi, *_) in rules.items():
        mol = smiles(test_smi)
        mol.canonicalize()
        pgs = mol.protective_groups
        assert name in pgs, f'{name} not detected in its own protected SMILES'
        # verify deletable atoms of detected PGs don't overlap
        for other in pgs:
            if other == name:
                continue
            q_other, keep_other = rules[other][0], rules[other][1]
            q_self = rules[name][0]
            self_atoms = set()
            for mp in q_self.get_mapping(mol, automorphism_filter=False):
                self_atoms.update(m for n, m in mp.items() if n not in keep)
            for mp in q_other.get_mapping(mol, automorphism_filter=False):
                other_atoms = {m for n, m in mp.items() if n not in keep_other}
                assert not self_atoms.intersection(other_atoms), (
                    f"'{other}' deletable atoms overlap with '{name}' on '{test_smi}'. "
                    f"Reorder rules in _protective.py so '{name}' comes before '{other}'."
                )


@pytest.mark.parametrize('smi,expected', [
    ('CC(C)(C)OC(=O)N(C(=O)OC(C)(C)C)c1ccccc1', {'amine_boc': 2}),
    ('CC(C)(C)OC(=O)Nc1ccc(OC(C)(C)C)cc1', {'amine_boc': 1, 'hydroxyl_tbu': 1}),
    ('c1ccc(COC(=O)Nc2ccc(OCc3ccccc3)cc2)cc1', {'amine_cbz': 1, 'hydroxyl_benzyl': 1}),
    ('CC(C)OCOCc1ccccc1', {'hydroxyl_bom': 1}),
    ('C(OCOC)(OCOC)C', {'hydroxyl_mom': 2}),
    ('OC(CO[Si](C)(C)C)O[Si](C)(C)C', {'hydroxyl_tms': 2}),
], ids=['double_boc', 'boc_and_tbu', 'cbz_and_benzyl', 'bom_no_benzyl', 'double_mom', 'double_tms'])
def test_multiple_pg_on_molecule(smi, expected):
    """Test correct detection of multiple PGs and no false positives from sub-pattern overlap."""
    mol = smiles(smi)
    mol.canonicalize()
    pgs = mol.protective_groups
    assert pgs == expected, f'got {pgs}'
