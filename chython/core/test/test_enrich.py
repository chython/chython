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
"""Merging two readings of one compound: `isomorphism`, `set_parity(order=)` and `enrich_from`.

The three are one feature and the middle one is why it is correct: a parity integer means "odd with
respect to THIS molecule's ref order", so copying it between two molecules inverts centres, and every
configuration here crosses through `order=` instead.

The hazard the fixtures exist for is cis-cyclobutane-1,3-diol -- two equivalent carbinol carbons, so
several isomorphisms and a stereo-blind labelling cannot say which one the source stated its parities
against.  A wrong choice there is a plausible diastereomer, which is why refusing the site is the
answer and a mapped count with no mismatches is not evidence of anything.
"""
from pytest import raises

from chython.core import (AutomorphismBudgetExceeded, INFO, LOST, MoleculeContainer, REFUSED,
                          REPAIRED, STEREO_ABS, STEREO_AND, STEREO_OR, STEREO_UNSPECIFIED,
                          read_smiles, read_smarts)


#: cis and trans cyclobutane-1,3-diol differ at ONE of two atoms that share an automorphism orbit.
CIS_DIOL = 'O[C@H]1C[C@@H](O)C1'
TRANS_DIOL = 'O[C@H]1C[C@H](O)C1'


def _rules(molecule):
    return [x.rule for x in molecule.log]


# --------------------------------------------------------------------------- isomorphism()

def test_the_same_constitution_maps_atom_for_atom_whatever_order_it_was_written_in():
    left = read_smiles('C[C@H](N)C(=O)O')
    right = read_smiles('OC(=O)[C@@H](C)N')
    mapping = left.isomorphism(right)
    assert sorted(mapping) == sorted(left)
    assert sorted(mapping.values()) == sorted(right)
    for n, m in mapping.items():
        assert left.element_of(n) == right.element_of(m)
    # and every bond travels with its order, which is what makes it an isomorphism rather than a
    # bag of atoms that happen to agree
    for bond in left.bonds():
        assert right.order_of(mapping[bond.n], mapping[bond.m]) == bond.order


def test_two_empty_molecules_map_to_the_empty_dict_and_not_to_none():
    """`None` is reserved for "not the same constitution", which is the one distinction the caller
    cannot afford to read as "took nothing"."""
    assert MoleculeContainer().isomorphism(MoleculeContainer()) == {}


def test_a_different_constitution_is_none_rather_than_a_wrong_map():
    assert read_smiles('CCO').isomorphism(read_smiles('CCN')) is None
    assert read_smiles('CCO').isomorphism(read_smiles('CCCO')) is None
    # same skeleton, one bond order apart
    assert read_smiles('CCC=O').isomorphism(read_smiles('CCCO')) is None


def test_the_three_normalizations_the_caller_owns_each_answer_none():
    """A kekule form, an explicit hydrogen and a protonation state are all differences in the RECORD
    rather than in the substance, and each makes a composition meaningless.  `kekule`/`thiele`,
    `implicify_hydrogens` and `neutralize` run before the call, not inside it."""
    assert read_smiles('c1ccccc1O').isomorphism(read_smiles('C1=CC=CC=C1O')) is None
    assert read_smiles('CCO').isomorphism(read_smiles('[H]C([H])([H])CO')) is None
    assert read_smiles('CC(=O)O').isomorphism(read_smiles('CC(=O)[O-]')) is None
    # and the first pair maps once the caller has done the normalization
    aromatic, kekule = read_smiles('c1ccccc1O'), read_smiles('C1=CC=CC=C1O')
    aromatic.kekule()
    assert aromatic.isomorphism(kekule) is not None


def test_an_isotope_a_charge_a_radical_and_an_r_index_each_break_the_map():
    assert read_smiles('CCO').isomorphism(read_smiles('[13C]CO')) is None
    assert read_smiles('CCO').isomorphism(read_smiles('CC[O-]')) is None
    assert read_smiles('CCO').isomorphism(read_smiles('[CH2]CO')) is None
    left, right = read_smiles('*CO'), read_smiles('*CO')
    with right.edit():
        right.set_r_index(1, 3)
    assert left.isomorphism(right) is None


def test_the_labelling_is_stereo_blind_so_two_diastereomers_still_map():
    """THE PREMISE OF THE FEATURE.  `canonical_order()` breaks its ties on the parities, so two
    readings that state different stereo get differently-tied labellings; this search runs with the
    stereo term off, which is what makes the composition a chosen isomorphism rather than any one."""
    cis, trans = read_smiles(CIS_DIOL), read_smiles(TRANS_DIOL)
    assert cis != trans
    assert cis.isomorphism(trans) is not None
    assert read_smiles('C/C=C/CO').isomorphism(read_smiles(r'C/C=C\CO')) is not None


def test_a_truncated_search_raises_rather_than_answering():
    """`canonical_order`'s posture, for `canonical_order`'s reason: a labelling from a truncated tree
    reads exactly like the right one."""
    cube = read_smiles('C12C3C4C1C5C4C3C25')
    with raises(AutomorphismBudgetExceeded):
        cube.isomorphism(read_smiles('C12C3C4C1C5C4C3C25'), _node_budget=2)
    assert cube.isomorphism(read_smiles('C12C3C4C1C5C4C3C25')) is not None


def test_a_query_is_not_a_source():
    with raises(TypeError):
        read_smiles('CCO').isomorphism(read_smarts('[C;a]'))


# -------------------------------------------------------------- set_parity(order=), primitive B

def test_translating_out_and_writing_back_in_the_same_order_is_a_fixed_point():
    """The pair `enrich_from` rests on.  The translation is an XOR by the presented permutation's
    parity, so the read and the write are one arithmetic run twice."""
    source = read_smiles('C[C@H](N)C(=O)O')
    unit = [x for x in source.stereo_units() if x['parity']][0]
    flat = read_smiles('CC(N)C(=O)O')
    assert flat.parity_of(unit['anchor']) == 0
    flat.set_parity(unit['anchor'], unit['parity'], order=unit['refs'])
    assert flat == source


def test_one_transposition_of_the_order_is_the_other_configuration():
    """WHY THE KEYWORD IS NOT COSMETIC.  A parity integer means "odd with respect to SOME order", and
    the same integer under two orders one transposition apart is two different compounds -- which is
    exactly what a parity copied bare between two molecules risks, since each molecule's refs follow
    its own neighbour ordering."""
    flat = read_smiles('CC(N)C(=O)O')
    refs = [x for x in read_smiles('C[C@H](N)C(=O)O').stereo_units() if x['parity']][0]['refs']
    assert refs == (1, 3, 4, None)
    straight, swapped = flat.copy(), flat.copy()
    straight.set_parity(2, 2, order=(1, 3, 4, None))
    swapped.set_parity(2, 2, order=(3, 1, 4, None))
    assert straight != swapped
    assert {straight.parity_of(2), swapped.parity_of(2)} == {1, 2}
    # the molecule's own order is what `order=None` means, so it agrees with the unpermuted one
    bare = flat.copy()
    bare.set_parity(2, 2)
    assert bare == straight


def test_order_none_is_the_meaning_the_method_already_had():
    left, right = read_smiles('CC(N)C(=O)O'), read_smiles('CC(N)C(=O)O')
    left.set_parity(2, 2)
    right.set_parity(2, 2, order=None)
    assert left == right


def test_clearing_a_parity_needs_no_frame_to_clear_it_against():
    """0 is 0 in every order, so `order=` on an unset does not go looking for a unit."""
    molecule = read_smiles('CCO')
    molecule.set_parity(2, 0, order=(1, 3, None, None))
    assert molecule.parity_of(2) == 0


def test_an_order_naming_a_direction_the_unit_has_not_got_is_refused():
    molecule = read_smiles('C[C@H](N)C(=O)O')
    with raises(ValueError, match='permutation'):
        molecule.set_parity(2, 1, order=(1, 3, 5, None))


def test_an_order_inside_a_structural_session_reads_a_frame_that_is_not_there_yet():
    molecule = read_smiles('C[C@H](N)C(=O)O')
    with raises(RuntimeError, match='changed the graph'):
        with molecule.edit():
            molecule.add_atom('Cl')
            molecule.set_parity(2, 1, order=(1, 3, 4, None))
    # a session of nothing but parities, wedges, collections and coordinates is fine, which is what
    # `enrich_from` is -- and is the reason the check is not `_require_clean`
    source = read_smiles('C[C@H](N)C(=O)O')
    unit = [x for x in source.stereo_units() if x['parity']][0]
    flat = read_smiles('CC(N)C(=O)O')
    with flat.edit():
        flat.set_stereo_group(2, STEREO_AND, 1)
        flat.set_parity(unit['anchor'], unit['parity'], order=unit['refs'])
    assert flat.parity_of(2) == source.parity_of(2)
    assert flat.stereo_groups() == {(STEREO_AND, 1): [2]}


# -------------------------------------------------------------------------------- enrich_from()

def test_a_configuration_travels_both_directions_and_arrives_unmirrored():
    for spelling in ('C[C@H](N)C(=O)O', 'C[C@@H](N)C(=O)O', 'C/C=C/CO', r'C/C=C\CO'):
        source = read_smiles(spelling)
        for flat_spelling in ('CC(N)C(=O)O', 'CC=CCO'):
            flat = read_smiles(flat_spelling)
            if flat.isomorphism(source) is None:
                continue
            assert flat.enrich_from(source)
            assert flat == source, spelling
            assert _rules(flat) == ['enrich:parity-borrowed']


def test_the_e_z_flip_is_the_test_that_a_pass_would_survive():
    """Both spellings of the same skeleton take, and they take DIFFERENT configurations -- a transfer
    that dropped the order would pass the first assertion and fail this one."""
    e, z = read_smiles('C/C=C/CO'), read_smiles(r'C/C=C\CO')
    from_e, from_z = read_smiles('CC=CCO'), read_smiles('CC=CCO')
    assert from_e.enrich_from(e) and from_z.enrich_from(z)
    assert from_e == e and from_z == z
    assert from_e != from_z


def test_nothing_is_taken_from_a_source_that_is_a_different_constitution():
    molecule = read_smiles('CC(N)C(=O)O')
    assert not molecule.enrich_from(read_smiles('C[C@H](N)C(=O)OC'))
    assert _rules(molecule) == ['enrich:not-the-same-constitution']
    assert molecule.log[0].severity == REFUSED
    assert 'nothing was taken' in molecule.log[0].message


def test_cis_cyclobutane_diol_is_refused_rather_than_quietly_mirrored():
    """THE REQUIRED FIXTURE.  Its two carbinol carbons share an automorphism orbit, so which of the
    source's two centres a given map calls "this one" is not settled by the constitution; the source
    states the opposite configuration at one of them and the site is refused."""
    cis, trans = read_smiles(CIS_DIOL), read_smiles(TRANS_DIOL)
    assert not cis.is_asymmetric(), 'the fixture stopped being the hazard it was chosen for'
    base = read_smiles(CIS_DIOL)
    assert not base.enrich_from(trans)
    assert _rules(base) == ['enrich:symmetry-unresolved']
    assert base.log[0].severity == REFUSED
    assert base == cis and base != trans
    # and with nothing of its own stated there is no ambiguity to resolve: any isomorphism carries
    # the source's whole configuration faithfully, because relabelling cannot change chirality
    flat = read_smiles('OC1CC(O)C1')
    assert flat.enrich_from(cis) and flat == cis
    flat = read_smiles('OC1CC(O)C1')
    assert flat.enrich_from(trans) and flat == trans


def test_where_the_orbits_leave_no_choice_the_disagreement_is_the_sources_and_the_policy_answers():
    """The other side of the same coin: an asymmetric constitution has exactly one isomorphism, so a
    disagreement there is a statement and not a tie.  `keep` logs it and writes nothing."""
    base, source = read_smiles('C[C@H](N)C(=O)O'), read_smiles('C[C@@H](N)C(=O)O')
    assert base.is_asymmetric()
    assert not base.enrich_from(source)
    assert _rules(base) == ['enrich:parity-conflict']
    assert base.log[0].severity == INFO
    assert base == read_smiles('C[C@H](N)C(=O)O')


def test_override_gives_the_disagreed_site_to_the_source_and_clear_gives_it_to_neither():
    source = read_smiles('C[C@@H](N)C(=O)O')
    won = read_smiles('C[C@H](N)C(=O)O')
    assert won.enrich_from(source, stereo='override')
    assert won == source
    assert _rules(won) == ['enrich:parity-overridden']

    wiped = read_smiles('C[C@H](N)C(=O)O')
    assert wiped.enrich_from(source, stereo='clear')
    assert wiped.parity_of(2) == 0
    assert wiped.stereo_group_of(2) == (STEREO_UNSPECIFIED, 0)
    assert _rules(wiped) == ['enrich:parity-cleared']
    assert wiped.log[0].severity == LOST


def test_clear_takes_the_collection_with_the_parity():
    """An AND membership holding no parity names a configuration that no longer exists, which is
    `clean_stereo()`'s own reason for wiping the two together."""
    base = read_smiles('C[C@H](N)C(=O)O')
    base.set_stereo_group(2, STEREO_AND, 1)
    assert base.enrich_from(read_smiles('C[C@@H](N)C(=O)O'), stereo='clear')
    assert base.parity_of(2) == 0
    assert not base.has_stereo_groups


def test_a_site_the_source_agrees_with_is_not_a_conflict_and_not_a_borrow():
    base, source = read_smiles('C[C@H](N)C(=O)O'), read_smiles('C[C@H](N)C(=O)O')
    assert not base.enrich_from(source)
    assert base.log == []


# ------------------------------------------------------------------------------- the collections

def test_a_collection_is_taken_whole_and_renumbered():
    """Two files numbering their `&1` differently must not have their racemates merged, so the id the
    source used is not the id that lands."""
    source = read_smiles('C[C@H](N)C(=O)O')
    source.set_stereo_group(2, STEREO_AND, 5)
    base = read_smiles('CC(N)C(=O)O')
    assert base.enrich_from(source)
    assert base.stereo_groups() == {(STEREO_AND, 1): [2]}
    assert 'enrich:group-borrowed' in _rules(base)


def test_a_collection_is_not_taken_onto_a_member_that_already_carries_one():
    source = read_smiles('C[C@H](N)C(=O)O')
    source.set_stereo_group(2, STEREO_OR, 2)
    base = read_smiles('C[C@H](N)C(=O)O')
    base.set_stereo_group(2, STEREO_AND, 1)
    assert not base.enrich_from(source)
    assert base.stereo_groups() == {(STEREO_AND, 1): [2]}
    assert _rules(base) == ['enrich:group-not-taken']
    assert base.log[0].severity == INFO


def test_an_abs_collection_keeps_its_id_because_abs_has_none():
    source = read_smiles('C[C@H](N)C(=O)O')
    source.set_stereo_group(2, STEREO_ABS)
    base = read_smiles('CC(N)C(=O)O')
    assert base.enrich_from(source)
    assert base.stereo_groups() == {(STEREO_ABS, 0): [2]}


# ------------------------------------------------------------------------------------ the plane

def test_the_plane_is_taken_when_this_record_has_none_and_kept_when_it_has_one():
    source = read_smiles('C[C@H](N)C(=O)O')
    source.clean2d()
    flat = read_smiles('C[C@H](N)C(=O)O')
    assert not flat.has_coordinates
    assert flat.enrich_from(source)
    mapping = flat.isomorphism(source)
    assert flat.coordinates() == {n: source.xy_of(mapping[n]) for n in flat}
    assert 'enrich:layout-borrowed' in _rules(flat)

    shifted = read_smiles('C[C@H](N)C(=O)O')
    shifted.clean2d()
    with shifted.edit():
        for n, (x, y) in list(shifted.coordinates().items()):
            shifted.set_xy(n, x + 10., y)
    drawn = read_smiles('C[C@H](N)C(=O)O')
    drawn.clean2d()
    before = drawn.coordinates()
    assert not drawn.enrich_from(shifted)
    assert drawn.coordinates() == before
    assert drawn.enrich_from(shifted, layout='override')
    assert drawn.coordinates() != before
    assert drawn.xy_box() == shifted.xy_box()


def test_an_unmapped_source_hands_over_no_plane_either():
    """A layout for a different graph is not a fallback, whatever `layout=` says."""
    drawn = read_smiles('CC(N)C(=O)OC')
    drawn.clean2d()
    flat = read_smiles('CC(N)C(=O)O')
    assert not flat.enrich_from(drawn, layout='override')
    assert not flat.has_coordinates


def test_a_wedge_drawn_against_the_old_plane_goes_with_the_plane():
    """A wedge is a statement about a drawing and asserts nothing about a different one.  Dropping it
    costs nothing: a writer derives a wedge from a parity when none is stored."""
    base = read_smiles('C[C@H](N)C(=O)O')
    base.clean2d()
    with base.edit():
        base.set_wedge(2, 3, 1)
    assert base.wedges()
    source = read_smiles('C[C@H](N)C(=O)O')
    source.clean2d()
    assert base.enrich_from(source, layout='override')
    assert base.wedges() == []
    assert 'wedge(s) dropped' in [x for x in base.log if x.rule == 'enrich:layout-borrowed'][0].message


def test_the_sources_wedges_come_only_where_its_configuration_won_as_well():
    """Carried in beside a kept parity they would make the CTfile writer emit the source's
    configuration, since it writes back the wedges a molecule holds rather than re-deriving them."""
    source = read_smiles('C[C@@H](N)C(=O)O')
    source.clean2d()
    with source.edit():
        source.set_wedge(2, 3, 1)
    kept = read_smiles('C[C@H](N)C(=O)O')
    assert kept.enrich_from(source, layout='override')
    assert kept.wedges() == []
    won = read_smiles('C[C@H](N)C(=O)O')
    assert won.enrich_from(source, layout='override', stereo='override')
    assert won.wedges() == [(2, 3, 1)] or won.wedges() == [(3, 2, 1)]


# -------------------------------------------------------------------------------- the signature

def test_sources_are_consumed_in_argument_order_and_the_base_still_wins():
    """A record carrying a molfile and a SMILES is one operation applied twice, and sequencing the
    calls is how a per-source policy is spelled."""
    drawn = read_smiles('CC(N)C(=O)O')
    drawn.clean2d()
    configured = read_smiles('C[C@H](N)C(=O)O')
    base = read_smiles('CC(N)C(=O)O')
    assert base.enrich_from(drawn, configured)
    assert base.has_coordinates and base.parity_of(2) == configured.parity_of(2)
    assert _rules(base) == ['enrich:layout-borrowed', 'enrich:parity-borrowed']


def test_a_policy_that_is_not_one_of_the_words_is_refused_rather_than_read_as_the_default():
    molecule = read_smiles('CCO')
    with raises(ValueError, match="'keep', 'override' or 'clear'"):
        molecule.enrich_from(read_smiles('CCO'), stereo='replace')
    with raises(ValueError, match="'keep' or 'override'"):
        molecule.enrich_from(read_smiles('CCO'), layout='clear')
    with raises(TypeError, match='MoleculeContainer'):
        molecule.enrich_from('C[C@H](N)C(=O)O')


def test_taking_nothing_is_false_and_says_so_in_the_return_value():
    molecule = read_smiles('CCO')
    assert molecule.enrich_from(read_smiles('CCO')) is False
    assert molecule.enrich_from() is False


def test_the_report_is_the_containers_own_log_and_nothing_is_conditional():
    """No `log=` parameter: `molecule.log` is the one destination, so a caller's rates -- borrowed
    parity, borrowed group, borrowed layout, conflicts -- are log analysis and not a return shape."""
    base = read_smiles('CC(N)C(=O)O')
    source = read_smiles('C[C@H](N)C(=O)O')
    source.clean2d()
    source.set_stereo_group(2, STEREO_AND, 7)
    assert base.enrich_from(source)
    assert sorted(_rules(base)) == ['enrich:group-borrowed', 'enrich:layout-borrowed',
                                    'enrich:parity-borrowed']
    assert {x.stage for x in base.log} == {'enrich'}
    assert {x.severity for x in base.log} == {REPAIRED}


def test_the_atom_count_does_not_move_which_is_what_separates_this_from_union():
    base = read_smiles('CC(N)C(=O)O')
    before = base.atom_count
    assert base.enrich_from(read_smiles('C[C@H](N)C(=O)O'))
    assert base.atom_count == before
    assert base.connected_components_count == 1
