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
"""The MACCS engine: the one-based vector, the seven predicates and the permanently-zero key."""
from numpy import uint8

from chython.chemistry import maccs_bit_set, maccs_keys
from chython.chemistry._maccs import MACCS_PREDICATE_FNS, aromatic_ring_count, maccs_match_counts
from chython.chemistry._tables import MACCS_PREDICATES, MACCS_UNSET_KEYS
from chython.core import read_smiles


def test_shape_dtype_and_the_reserved_zero():
    v = maccs_keys(read_smiles('CC(=O)Nc1ccc(O)cc1'))
    assert v.dtype == uint8
    assert v.shape == (167,)
    assert v[0] == 0
    assert set(v.tolist()) <= {0, 1}


def test_the_vector_is_one_based():
    # keys[n] is published key n.  no caller writes n - 1.
    v = maccs_keys(read_smiles('c1ccccc1'))
    assert all(v[n] in (0, 1) for n in range(1, 167))


def test_the_bit_set_and_the_vector_are_the_same_information():
    m = read_smiles('CC(=O)Oc1ccccc1C(=O)O')
    v = maccs_keys(m)
    assert maccs_bit_set(m) == frozenset(n for n in range(1, 167) if v[n])
    assert 0 not in maccs_bit_set(m)


def test_the_reserved_zero_is_never_set_by_any_molecule():
    for smi in ('C', 'O', 'c1ccccc1', 'CC(=O)[O-].[Na+]', '[13CH4]', 'CC(=O)Nc1ccc(O)cc1'):
        assert maccs_keys(read_smiles(smi))[0] == 0


def test_the_predicate_registry_matches_the_tables_vocabulary():
    assert set(MACCS_PREDICATE_FNS) == set(MACCS_PREDICATES)
    assert all(callable(f) for f in MACCS_PREDICATE_FNS.values())


def test_a_key_with_no_published_definition_is_never_set():
    """`MACCS_UNSET_KEYS` is documented as permanently zero, so the engine must never set one.

    Key 44's published description is the placeholder `OTHER`; there is nothing to transcribe and a
    guessed pattern would be invented chemistry.  This is the test that turns "we chose not to
    implement it" into a checked property rather than a comment: the `unset` branch could be replaced
    by a catch-all that sets the bit, and this is what would catch that.
    """
    for smi in ('C', 'O', 'c1ccccc1', 'CC(=O)[O-].[Na+]', '[13CH4]', 'CC(=O)Nc1ccc(O)cc1',
                'NS(=O)(=O)N', 'c1ccc2ccccc2c1'):
        v = maccs_keys(read_smiles(smi))
        for key in MACCS_UNSET_KEYS:
            assert v[key] == 0, (smi, key)
        assert not set(MACCS_UNSET_KEYS) & maccs_bit_set(read_smiles(smi))
        assert not set(MACCS_UNSET_KEYS) & set(maccs_match_counts(read_smiles(smi)))


def test_each_predicate_answers_the_question_it_names():
    fns = MACCS_PREDICATE_FNS
    assert fns['isotope'](read_smiles('[13CH4]'))
    assert not fns['isotope'](read_smiles('C'))
    assert fns['charge'](read_smiles('CC(=O)[O-]'))
    assert not fns['charge'](read_smiles('CC(=O)O'))
    assert fns['fragments_gt_1'](read_smiles('CC(=O)[O-].[Na+]'))
    assert not fns['fragments_gt_1'](read_smiles('CC(=O)O'))
    assert fns['ring_present'](read_smiles('c1ccccc1'))
    assert not fns['ring_present'](read_smiles('CCCC'))
    assert fns['atomic_number_gt_103'](read_smiles('[Db]'))
    assert not fns['atomic_number_gt_103'](read_smiles('[U]'))

    # keys 125 and 145: RINGS, not atoms.  benzene is the near miss for both.
    assert fns['aromatic_rings_gt_1'](read_smiles('c1ccc2ccccc2c1'))          # naphthalene, 2
    assert fns['aromatic_rings_gt_1'](read_smiles('c1ccc(-c2ccccc2)cc1'))     # biphenyl, 2
    assert not fns['aromatic_rings_gt_1'](read_smiles('c1ccccc1'))            # benzene, 1
    assert not fns['aromatic_rings_gt_1'](read_smiles('C1CCCCC1'))            # cyclohexane, 0
    assert fns['six_rings_gt_1'](read_smiles('c1ccc2ccccc2c1'))               # naphthalene, 2
    assert not fns['six_rings_gt_1'](read_smiles('c1ccccc1'))                 # benzene, 1
    assert not fns['six_rings_gt_1'](read_smiles('C1CCCC1'))                  # cyclopentane, 0


def test_the_aromatic_ring_count_is_the_containers_own_answer():
    """A delegation, so it must never differ.  A second aromatic-ring answer in the tree would drift."""
    for smi, want in (('c1ccccc1', 1), ('c1ccc2ccccc2c1', 2), ('C1CCCCC1', 0),
                      ('c1ccc(-c2ccccc2)cc1', 2), ('c1ccncc1', 1)):
        m = read_smiles(smi)
        assert aromatic_ring_count(m) == want == m.aromatic_rings_count, smi


def test_it_is_renumbering_invariant():
    a = maccs_bit_set(read_smiles('CC(=O)Nc1ccc(O)cc1'))
    b = maccs_bit_set(read_smiles('Oc1ccc(NC(C)=O)cc1'))
    assert a == b


def test_the_vector_does_not_change_when_hydrogens_become_explicit():
    """A descriptor reports on the molecule, not on the drawing.

    Every wildcard in the table carries `!#1` for this reason: `[*]` matches an explicit hydrogen, so a
    bare-wildcard path key would answer differently on the same compound drawn two ways.
    """
    m = read_smiles('CC(=O)Nc1ccc(O)cc1')
    implicit = maccs_bit_set(m)
    with m.edit() as e:
        e.explicify_hydrogens()
    assert maccs_bit_set(m) == implicit


def test_the_match_counts_only_report_keys_that_matched():
    counts = maccs_match_counts(read_smiles('c1ccccc1'))
    assert all(v > 0 for v in counts.values())
    assert set(counts) <= set(range(1, 167))


def test_the_container_methods_agree_with_the_functions():
    m = read_smiles('c1ccccc1')
    assert (m.maccs_keys() == maccs_keys(m)).all()
    assert m.maccs_bit_set() == maccs_bit_set(m)
