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
from importlib.util import find_spec
from pytest import mark, raises
from chython.chemistry import pharmacophore_invariants
from chython.chemistry._pharmacophore import (PH_ACCEPTOR, PH_AROMATIC, PH_DONOR, PH_HYDROPHOBE,
                                              PH_NEGATIVE, PH_POSITIVE, PH_TYPES,
                                              pharmacophore_atoms)
from chython.core import read_smiles


# `pharmacophore_invariants` ANSWERS AN ARRAY AND `pharmacophore_atoms` ANSWERS STABLE IDS, and that is
# the whole of why this marker is per-test and not a module-level `pytestmark`: numpy is optional
# (`chython[ml]`), only the vector needs it, and the reachability gate below -- the nine rows each having
# to win a substrate -- is the part of this file most worth still running on a minimal install, since a
# dead row costs a feature type on real input and moves no other assertion.
#
# `numpy` IS NOT IMPORTED AT MODULE LEVEL EITHER, and a marker alone could not have fixed that: pytest
# imports a module to collect it and reads its marks afterwards, so `from numpy import uint32` up here
# was a collection ERROR that no skip can reach.  The one test that needs the dtype imports it itself.
# `find_spec` rather than `importorskip`, so collection does not import numpy at all -- the same
# reasoning as `interop/test/conftest.py` gives for the optional toolkits.
needs_numpy = mark.skipif(find_spec('numpy') is None,
                          reason='numpy is not installed; the invariant vector is an array')


@needs_numpy
def test_dtype_shape_and_order():
    from numpy import uint32
    m = read_smiles('CC(=O)Nc1ccc(O)cc1')  # paracetamol
    v = pharmacophore_invariants(m)
    assert v.dtype == uint32
    assert v.ndim == 1
    assert v.shape == (m.atom_count,)


@needs_numpy
def test_every_alkane_carbon_is_hydrophobic_and_nothing_else():
    m = read_smiles('CCCC')  # butane: hydrophobic carbons only, no donor/acceptor/charge/aromatic
    v = pharmacophore_invariants(m)
    assert set(v.tolist()) == {PH_HYDROPHOBE}


@needs_numpy
def test_an_atom_with_no_feature_is_zero():
    # The only test that observes the `0`, and butane cannot be it: all four of its carbons are `x0`,
    # so row 9 types every one.  This methyl carbon has a heteroatom neighbour, so row 9 rejects it.
    #
    # The nitrogen is `PH_DONOR | PH_POSITIVE`, not `PH_POSITIVE` alone: `[NH3+]` carries h=3, so
    # hbond.tsv's row 1 types it a donor, and the roles are ADDITIVE.  Do not "fix" a failure here by
    # narrowing the donor row -- `test_the_cations_donate_up_to_four_hydrogens_and_accept_none` pins
    # this same atom.
    m = read_smiles('C[NH3+]')                       # atom_numbers order: [C, N]
    assert pharmacophore_invariants(m).tolist() == [0, PH_DONOR | PH_POSITIVE]


@needs_numpy
def test_a_bare_methane_carbon_is_hydrophobic_and_nothing_else():
    assert pharmacophore_invariants(read_smiles('C')).tolist() == [PH_HYDROPHOBE]


@needs_numpy
def test_an_alcohol_oxygen_is_both_donor_and_acceptor():
    m = read_smiles('CCO')
    o = next(a.n for a in m.atoms() if a.element == 8)
    v = pharmacophore_invariants(m)
    assert v[m.atom_numbers.index(o)] == PH_DONOR | PH_ACCEPTOR


def test_a_carboxylate_oxygen_is_negative_and_an_ammonium_nitrogen_positive():
    a = pharmacophore_atoms(read_smiles('CC(=O)[O-]'))
    assert len(a['negative']) == 2          # both oxygens of the delocalised carboxylate
    b = pharmacophore_atoms(read_smiles('C[NH3+]'))
    assert len(b['positive']) == 1


def test_a_dative_oxide_oxygen_is_not_a_negative_feature():
    """A formal charge is not an ionisation state.

    `[O;D1;-]` alone made an N-oxide, a nitro group and a charge-separated sulfoxide anionic centres.
    None is an anion -- the charge is valence bookkeeping on a neutral molecule -- and a model that
    places a negative feature there is looking for a counter-ion that does not exist.  Row 4 now names
    the neutral carbon, phosphorus or sulfur a deprotonated oxygen hangs off, so a cationic centre
    excludes itself by carrying the charge.
    """
    for smi in ('C[N+](C)(C)[O-]', 'c1cc[n+]([O-])cc1', 'c1ccccc1[N+](=O)[O-]',
                'C[S+]([O-])C', 'C[P+](C)(C)[O-]'):
        assert pharmacophore_atoms(read_smiles(smi))['negative'] == frozenset(), smi


@mark.parametrize('smi,count', [('CC(=O)[O-]', 2),                  # carboxylate, both oxygens
                                ('c1ccccc1[O-]', 1),                # phenoxide
                                ('c1ccccc1S(=O)(=O)[O-]', 3),       # benzenesulfonate
                                ('COP(=O)([O-])[O-]', 3),           # methyl phosphate
                                ('C[S-]', 1),                       # thiolate, row 6
                                ('CC(=O)[N-]C', 1)])                # deprotonated amide, row 7
def test_the_real_anions_keep_their_negative_feature(smi, count):
    """The other half of the narrowing: what row 4 exists for must still reach it."""
    assert len(pharmacophore_atoms(read_smiles(smi))['negative']) == count


@needs_numpy
def test_a_benzene_carbon_is_aromatic_and_hydrophobic():
    v = pharmacophore_invariants(read_smiles('c1ccccc1'))
    assert set(v.tolist()) == {PH_AROMATIC | PH_HYDROPHOBE}


@needs_numpy
def test_a_pyridine_nitrogen_is_aromatic_and_an_acceptor_and_not_hydrophobic():
    m = read_smiles('c1ccncc1')
    n = next(a.n for a in m.atoms() if a.element == 7)
    got = pharmacophore_invariants(m)[m.atom_numbers.index(n)]
    assert got & PH_AROMATIC and got & PH_ACCEPTOR
    assert not got & PH_HYDROPHOBE


def test_every_type_name_has_a_bit_and_the_bits_are_distinct():
    bits = (PH_DONOR, PH_ACCEPTOR, PH_POSITIVE, PH_NEGATIVE, PH_AROMATIC, PH_HYDROPHOBE)
    assert len(PH_TYPES) == len(bits) == 6
    assert len(set(bits)) == 6
    assert all(b and not (b & (b - 1)) for b in bits)   # each is a single power of two


@needs_numpy
def test_the_vector_is_accepted_as_invariants_by_a_fingerprint():
    m = read_smiles('CC(=O)Nc1ccc(O)cc1')
    fp = m.morgan_fingerprint(invariants=pharmacophore_invariants(m))
    assert fp.any()


@needs_numpy
def test_it_is_renumbering_invariant_as_a_multiset():
    from collections import Counter
    a = pharmacophore_invariants(read_smiles('CC(=O)Nc1ccc(O)cc1'))
    b = pharmacophore_invariants(read_smiles('Oc1ccc(NC(C)=O)cc1'))
    assert Counter(a.tolist()) == Counter(b.tolist())


@needs_numpy
def test_the_container_method_agrees():
    m = read_smiles('CCO')
    assert (m.pharmacophore_invariants() == pharmacophore_invariants(m)).all()


def test_donor_and_acceptor_have_exactly_one_definition_in_the_tree():
    # pharmacophore.tsv must not re-spell what hbond.tsv already says, so the two tables are asserted
    # disjoint rather than merely believed to be.
    from chython.chemistry._counts import hbond_atoms
    from chython.chemistry._tables import (HBOND_ROLES, PHARMACOPHORE_ROLES,
                                           pharmacophore_rules)
    assert not set(PHARMACOPHORE_ROLES) & set(HBOND_ROLES)
    assert not {r.role for r in pharmacophore_rules()} & {'donor', 'acceptor'}
    # and the two keys really are hbond.tsv's answer, not a copy of it
    m = read_smiles('CC(=O)Nc1ccc(O)cc1')
    a = pharmacophore_atoms(m)
    assert a['donor'] == hbond_atoms(m, 'donor')
    assert a['acceptor'] == hbond_atoms(m, 'acceptor')


def test_all_six_keys_are_present_even_when_empty():
    a = pharmacophore_atoms(read_smiles('C'))
    assert set(a) == set(PH_TYPES)
    assert all(isinstance(v, frozenset) for v in a.values())


# The reachability gate.  Every other test above observes a molecule, so a row that matches nothing
# moves no assertion and is a silent defect -- it costs a feature type on real input and looks fine
# forever.  Each of the nine rows therefore names a substrate it must win.
PROBES = {
    'pharmacophore:1': ('C[N+](C)(C)C', 'positive'),      # quaternary ammonium
    'pharmacophore:2': ('C[NH3+]', 'positive'),           # protonated amine
    'pharmacophore:3': ('CC(N)=N', 'positive'),           # acetamidine -- D1 h1 imine N
    'pharmacophore:4': ('CC(=O)[O-]', 'negative'),        # carboxylate anion oxygen
    'pharmacophore:5': ('CC(=O)[O-]', 'negative'),        # its formally neutral partner oxygen
    'pharmacophore:6': ('C[S-]', 'negative'),             # thiolate
    'pharmacophore:7': ('CC(=O)[N-]C', 'negative'),       # deprotonated N-methylacetamide, D2
    'pharmacophore:8': ('c1ccccc1', 'aromatic'),          # benzene carbon
    'pharmacophore:9': ('CCCC', 'hydrophobe'),            # alkane carbon
}


@mark.parametrize('rule_id', list(PROBES))
def test_every_row_matches_its_own_probe(rule_id):
    from chython.chemistry._tables import pharmacophore_rules
    smi, role = PROBES[rule_id]
    row = next(r for r in pharmacophore_rules() if r.id == rule_id)
    assert row.role == role, 'the probe table and the TSV disagree about this row\'s role'
    # `is_substructure`, NOT `bool(get_mapping(...))`: `get_mapping` returns a generator, and a
    # generator object is truthy whether or not it will yield, which would make this gate unfailable.
    assert row.query.is_substructure(read_smiles(smi)), f'{rule_id} matches nothing: a dead row'


def test_the_reachability_gate_can_fail():
    # negative control for the line above: `is_substructure` returns a real `bool`, while
    # `bool(get_mapping(...))` is `True` even for a xenon query against butane.
    from chython.core import read_smarts
    dead, alkane = read_smarts('[Xe;*:1]'), read_smiles('CCCC')
    assert dead.is_substructure(alkane) is False
    assert bool(dead.get_mapping(alkane)) is True, 'the generator is truthy: this is what made the gate unfailable'


def test_the_probe_table_covers_every_row():
    # without this, deleting a row from PROBES silently retires its reachability check.
    from chython.chemistry._tables import pharmacophore_rules
    assert set(PROBES) == {r.id for r in pharmacophore_rules()}


def test_the_two_carboxylate_rows_type_different_oxygens():
    # rows 4 and 5 share a probe, so `test_every_row_matches_its_own_probe` cannot tell them apart;
    # this is what says they are two rows rather than one written twice.
    from chython.chemistry._tables import pharmacophore_rules
    m = read_smiles('CC(=O)[O-]')
    hit = {}
    for rid in ('pharmacophore:4', 'pharmacophore:5'):
        row = next(r for r in pharmacophore_rules() if r.id == rid)
        s = row.numbers[1]
        hit[rid] = {mapping[s] for mapping in row.query.get_mapping(m)}
    assert len(hit['pharmacophore:4']) == len(hit['pharmacophore:5']) == 1
    assert not hit['pharmacophore:4'] & hit['pharmacophore:5']
    assert pharmacophore_atoms(m)['negative'] == hit['pharmacophore:4'] | hit['pharmacophore:5']
