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
"""The fingerprint substrate: the atom invariant, the folder, and the two enumerators."""
from importlib.util import find_spec
from pytest import mark, raises

from chython.core import MoleculeContainer, read_smiles


# numpy is an optional dependency (`chython[ml]`), and EVERY spelling in both families needs it -- the
# `*_hash_set`, `*_hash_counts` and `*_bit_set` ones included.  They answer a plain set or dict, so the
# return type says nothing; they all build the same uint32 invariant vector on the way there.  The
# marker therefore lands on all but the twelve tests that never reach the vector: the argument
# validators, which refuse before the walk begins, and the three surface tests that only read
# `dir(MoleculeContainer)`.  Which twelve those are was measured, not reasoned about.
#
# `find_spec` rather than `importorskip`, so collection does not import numpy at all -- the same
# reasoning as `interop/test/conftest.py` gives for the optional toolkits.
needs_numpy = mark.skipif(find_spec('numpy') is None,
                          reason='numpy is not installed; every fingerprint spelling needs it')


# --- atom_invariants -----------------------------------------------------------------------------

@needs_numpy
def test_atom_invariants_is_one_uint32_per_atom_in_atom_order():
    mol = read_smiles('CCO')
    inv = mol.atom_invariants()
    assert inv.shape == (3,)
    assert inv.dtype.name == 'uint32'


@needs_numpy
def test_atom_invariants_of_an_empty_molecule_is_empty():
    assert MoleculeContainer().atom_invariants().shape == (0,)


@needs_numpy
def test_atom_invariants_agree_across_equivalent_atoms():
    # every benzene carbon is the same atom: aromatic, one hydrogen, two heavy neighbours, in a ring
    inv = read_smiles('c1ccccc1').atom_invariants()
    assert len(set(inv.tolist())) == 1


@needs_numpy
def test_atom_invariants_separate_atoms_that_differ_only_in_hydrogen_count():
    # propane: the two methyls carry three hydrogens, the middle carbon two
    inv = read_smiles('CCC').atom_invariants()
    assert inv[0] == inv[2] != inv[1]


@needs_numpy
def test_atom_invariants_separate_atoms_that_differ_only_by_ring_membership():
    # a cyclohexane CH2 and a propane CH2 agree on element, charge, hydrogens and degree
    ring = read_smiles('C1CCCCC1').atom_invariants()
    chain = read_smiles('CCC').atom_invariants()
    assert ring[0] != chain[1]


@needs_numpy
def test_atom_invariants_read_an_unknown_hydrogen_count_as_zero():
    # the ONE place in the tree where the sentinel collapses, and it is deliberate: a fingerprint is
    # a screen, not a stored field, so an underivable count must not be a value of its own
    unknown = MoleculeContainer()
    with unknown.edit():
        unknown.add_atom('C')                    # the argument omitted says nothing
    stated = MoleculeContainer()
    with stated.edit():
        stated.add_atom('C', implicit_h=0)       # zero, stated out loud
    assert unknown.unknown_h_count == 1
    assert stated.unknown_h_count == 0
    assert unknown.atom_invariants()[0] == stated.atom_invariants()[0]


@needs_numpy
def test_a_record_with_underivable_hydrogens_labels_like_its_curated_form():
    # the case that decides the ruling above.  An MDL record of a Suzuki palladium catalyst whose
    # hydrogen counts nobody could derive has to screen against the curated form of the same
    # catalyst, where those counts are zero -- retrieval is what the screen is wanted for.
    raw = MoleculeContainer()
    with raw.edit():
        centre = raw.add_atom('Pd')                        # nothing said
        for _ in range(2):
            donor = raw.add_atom('P', implicit_h=0)
            raw.add_bond(centre, donor, 8)
    curated = MoleculeContainer()
    with curated.edit():
        centre = curated.add_atom('Pd', implicit_h=0)      # curated to zero
        for _ in range(2):
            donor = curated.add_atom('P', implicit_h=0)
            curated.add_bond(centre, donor, 8)
    assert raw.unknown_h_count == 1
    assert curated.unknown_h_count == 0
    assert raw.atom_invariants().tolist() == curated.atom_invariants().tolist()


@needs_numpy
def test_atom_invariants_sum_implicit_and_explicit_hydrogens():
    # `C` and its written-out spelling are one molecule, so methane's CARBON is one label.  The
    # reader keeps the four hydrogens as atoms -- measured -- so the carbon holds them as
    # `explicit_h=4, implicit_h=0` where bare methane holds `implicit_h=4, explicit_h=0`, and both
    # its degree of four and its degree of zero must read as the same heavy-atom degree of zero.
    #
    # THIS IS A CLAIM ABOUT THE LABEL AND NOT ABOUT A FINGERPRINT.  The enumerators walk the CSR,
    # where written-out methane genuinely has five nodes and bare methane has one; no atom label can
    # hide that and none should.  `implicify_hydrogens()` is the caller's normalisation step.
    implicit = read_smiles('C')
    explicit = read_smiles('[H]C([H])([H])[H]')
    assert explicit.element_of(2) == 6                    # atom 2 is the carbon; 1, 3, 4, 5 are H
    assert (implicit.implicit_h_of(1), implicit.explicit_h_of(1)) == (4, 0)
    assert (explicit.implicit_h_of(2), explicit.explicit_h_of(2)) == (0, 4)
    assert implicit.atom_invariants()[0] == explicit.atom_invariants()[1]


@needs_numpy
def test_atom_invariants_do_not_move_when_atom_numbers_do():
    mol = read_smiles('CCO')
    before = mol.atom_invariants().tolist()
    mol.remap({1: 7, 2: 8, 3: 9})
    assert mol.atom_invariants().tolist() == before


# --- morgan, unfolded ----------------------------------------------------------------------------

@needs_numpy
def test_morgan_one_shell_reports_every_atom_exactly_once():
    # a shell is a fragment centred on an atom, so the counts of one radius sum to the atom count
    mol = read_smiles('CC(=O)OCC')
    for radius in (1, 2, 3, 4):
        counts = mol.morgan_hash_counts(radius, radius)
        assert sum(counts.values()) == mol.atom_count, radius


@needs_numpy
def test_morgan_counts_sum_over_every_requested_shell():
    mol = read_smiles('CC(=O)OCC')
    counts = mol.morgan_hash_counts(1, 4)
    assert sum(counts.values()) == mol.atom_count * 4


@needs_numpy
def test_morgan_radius_one_is_the_atom_label_and_nothing_else():
    # benzene: one label, six atoms carrying it, so radius 1 is a single key of count six
    counts = read_smiles('c1ccccc1').morgan_hash_counts(1, 1)
    assert list(counts.values()) == [6]


@needs_numpy
def test_morgan_radius_one_counts_the_symmetry_of_neopentane():
    counts = read_smiles('CC(C)(C)C').morgan_hash_counts(1, 1)
    assert sorted(counts.values()) == [1, 4]


@needs_numpy
def test_morgan_growing_the_radius_only_adds_fragments():
    mol = read_smiles('CC(C)C(=O)Nc1ccccc1')
    assert mol.morgan_hash_set(1, 2) < mol.morgan_hash_set(1, 4)


@needs_numpy
def test_morgan_is_invariant_to_the_written_atom_order():
    a = read_smiles('CC(=O)OCC')
    b = read_smiles('CCOC(C)=O')
    assert a == b                                   # one molecule, two spellings
    assert a.morgan_hash_counts() == b.morgan_hash_counts()


@needs_numpy
def test_morgan_hash_set_is_the_keys_of_the_counts():
    mol = read_smiles('c1ccccc1O')
    assert mol.morgan_hash_set(1, 3) == set(mol.morgan_hash_counts(1, 3))


@needs_numpy
def test_morgan_of_an_empty_molecule_is_empty():
    assert MoleculeContainer().morgan_hash_counts() == {}
    assert MoleculeContainer().morgan_hash_set() == set()


@needs_numpy
def test_morgan_of_one_atom_has_one_fragment_per_shell():
    # nothing to expand into, so every shell is the same atom -- and the hashes still differ,
    # because each shell hashes the previous shell's identifier rather than the raw atom label
    counts = read_smiles('[Na+]').morgan_hash_counts(1, 4)
    assert sum(counts.values()) == 4


@needs_numpy
def test_morgan_separates_two_molecules_that_share_every_atom_label():
    # same brutto and the same multiset of atom labels; the difference is where the bonds go
    assert read_smiles('CCCCO').morgan_hash_set(1, 4) != read_smiles('CC(C)CO').morgan_hash_set(1, 4)


@needs_numpy
def test_morgan_reads_a_caller_supplied_label_vector():
    mol = read_smiles('CCO')
    flat = mol.atom_invariants()
    flat[:] = 1                                     # every atom the same label
    assert mol.morgan_hash_counts(1, 1, invariants=flat) != mol.morgan_hash_counts(1, 1)
    assert list(mol.morgan_hash_counts(1, 1, invariants=flat).values()) == [3]


@needs_numpy
def test_a_frozen_invariants_array_works_in_both_families():
    # a cached pharmacophore array is naturally read-only; the enumerators never write through the
    # view and must accept it rather than raising "buffer source array is read-only"
    from numpy import array
    mol = read_smiles('CCO')
    writable = mol.atom_invariants()
    frozen = array(writable)
    frozen.setflags(write=False)
    assert frozen.flags['WRITEABLE'] is False
    # both families must accept the frozen array and return the same result as the writable one
    assert mol.morgan_hash_counts(invariants=frozen) == mol.morgan_hash_counts(invariants=writable)
    assert mol.linear_hash_counts(invariants=frozen) == mol.linear_hash_counts(invariants=writable)


@needs_numpy
def test_morgan_refuses_a_label_vector_of_the_wrong_length():
    from numpy import zeros
    with raises(ValueError, match='one entry per atom'):
        read_smiles('CCO').morgan_hash_counts(invariants=zeros(4, dtype='uint32'))


@needs_numpy
def test_morgan_refuses_a_label_vector_of_the_wrong_dtype():
    from numpy import zeros
    with raises(ValueError, match='uint32'):
        read_smiles('CCO').morgan_hash_counts(invariants=zeros(3, dtype='int64'))


def test_morgan_refuses_a_radius_below_one():
    # unvalidated, min_radius=0 is silently 2 in the linear family, so the bound is checked
    with raises(ValueError, match='min_radius'):
        read_smiles('CCO').morgan_hash_counts(0, 4)


def test_morgan_refuses_a_max_radius_below_the_min():
    with raises(ValueError, match='max_radius'):
        read_smiles('CCO').morgan_hash_counts(4, 2)


# --- folding -------------------------------------------------------------------------------------

@needs_numpy
def test_morgan_fingerprint_is_a_binary_uint8_vector_of_the_requested_length():
    fp = read_smiles('c1ccccc1O').morgan_fingerprint(length=512)
    assert fp.shape == (512,)
    assert fp.dtype.name == 'uint8'
    assert set(fp.tolist()) <= {0, 1}
    assert fp.sum() > 0


@needs_numpy
def test_morgan_count_vector_is_a_uint32_vector_of_the_requested_length():
    cv = read_smiles('c1ccccc1O').morgan_count_vector(length=512)
    assert cv.shape == (512,)
    assert cv.dtype.name == 'uint32'


@needs_numpy
def test_the_three_folded_spellings_light_up_the_same_positions():
    mol = read_smiles('CC(C)C(=O)Nc1ccccc1')
    bits = mol.morgan_bit_set()
    assert set(mol.morgan_fingerprint().nonzero()[0].tolist()) == bits
    assert set(mol.morgan_count_vector().nonzero()[0].tolist()) == bits


@needs_numpy
def test_the_count_vector_carries_at_least_as_much_weight_as_the_binary_one():
    mol = read_smiles('CCCCCCCCCC')          # decane: one fragment repeated many times
    assert mol.morgan_count_vector().sum() > mol.morgan_fingerprint().sum()


# NOTE: the uint32 saturation branch in fp_fold_counted is deliberately untested.  Reaching it
# needs a position's accumulated count to exceed 4 294 967 295 -- four billion copies of one
# fragment in one molecule.  No public path produces that, and exposing a cdef function or a
# synthetic-dict entry point purely to exercise an unreachable branch would be a worse trade than
# the coverage gap.  The clamp is there because a wrap would turn the commonest fragment into the
# rarest; the comment in _fingerprints.pxi documents the same.


@needs_numpy
def test_the_count_vector_total_is_the_unfolded_total_times_the_active_bits():
    # nothing is dropped by folding and nothing is invented: a collision adds, it does not replace
    mol = read_smiles('CC(=O)OCC')
    unfolded = sum(mol.morgan_hash_counts().values())
    assert mol.morgan_count_vector(number_active_bits=2).sum() == unfolded * 2


@needs_numpy
def test_each_active_bit_reads_a_distinct_window_of_the_hash():
    # folds the unfolded dict by hand in Python and compares POSITIONS, not totals.
    # `test_the_count_vector_total_is_the_unfolded_total_times_the_active_bits` does not catch an
    # implementation that reads window 0 for every active bit -- the total still sums correctly.
    # length=1024, number_active_bits=3 gives three distinct 10-bit windows [0:10], [10:20], [20:30].
    # (confirmed: against a zeroed-shift defect the expected set is 98 positions vs 33 for the defect)
    mol = read_smiles('CC(=O)Nc1ccc(O)cc1')   # paracetamol
    length = 1024
    number_active_bits = 3
    width = length.bit_length() - 1           # 10
    mask = length - 1
    counts = mol.morgan_hash_counts(1, 4)
    expected = set()
    for h in counts:
        for i in range(number_active_bits):
            expected.add((h >> (i * width)) & mask)
    assert mol.morgan_bit_set(1, 4, length, number_active_bits) == expected
    cv = mol.morgan_count_vector(1, 4, length, number_active_bits)
    assert set(cv.nonzero()[0].tolist()) == expected


@needs_numpy
def test_folding_cannot_produce_more_positions_than_it_was_given_hashes():
    mol = read_smiles('CC(C)C(=O)Nc1ccccc1')
    assert len(mol.morgan_bit_set()) <= len(mol.morgan_hash_set()) * 2


@needs_numpy
def test_a_longer_fingerprint_collides_less():
    mol = read_smiles('CC(C)C(=O)Nc1ccccc1O')
    assert len(mol.morgan_bit_set(length=256)) <= len(mol.morgan_bit_set(length=4096))


@needs_numpy
def test_one_active_bit_lights_at_most_one_position_per_hash():
    mol = read_smiles('CC(C)C(=O)Nc1ccccc1')
    assert len(mol.morgan_bit_set(number_active_bits=1)) <= len(mol.morgan_hash_set())


@needs_numpy
def test_tanimoto_ranks_the_closer_pair_higher():
    phenol = read_smiles('c1ccccc1O')
    aniline = read_smiles('c1ccccc1N')
    hexane = read_smiles('CCCCCC')

    def tanimoto(a, b):
        x, y = a.morgan_bit_set(), b.morgan_bit_set()
        return len(x & y) / len(x | y)

    assert tanimoto(phenol, aniline) > tanimoto(phenol, hexane)


@needs_numpy
def test_the_folded_spellings_of_an_empty_molecule_are_empty_not_absent():
    mol = MoleculeContainer()
    assert mol.morgan_bit_set() == set()
    assert mol.morgan_fingerprint().shape == (1024,)
    assert mol.morgan_fingerprint().sum() == 0
    assert mol.morgan_count_vector().sum() == 0


@needs_numpy
def test_folding_is_invariant_to_the_written_atom_order():
    a = read_smiles('CC(=O)OCC')
    b = read_smiles('CCOC(C)=O')
    assert (a.morgan_fingerprint() == b.morgan_fingerprint()).all()
    assert (a.morgan_count_vector() == b.morgan_count_vector()).all()


# --- the folding validation ----------------------------------------------------------------------

def test_a_length_that_is_not_a_power_of_two_is_refused():
    with raises(ValueError, match='power of two'):
        read_smiles('CCO').morgan_fingerprint(length=1000)


def test_a_length_below_two_is_refused():
    with raises(ValueError, match='power of two'):
        read_smiles('CCO').morgan_fingerprint(length=1)


def test_zero_active_bits_is_refused():
    with raises(ValueError, match='number_active_bits'):
        read_smiles('CCO').morgan_fingerprint(number_active_bits=0)


def test_more_active_bits_than_the_hash_can_pay_for_is_refused():
    # seven ten-bit slices need seventy bits and a fragment hash is sixty-four wide, so an unvalidated
    # fold returns a fingerprint whose last bits are all zero
    with raises(ValueError, match='64 bits'):
        read_smiles('CCO').morgan_fingerprint(length=1024, number_active_bits=7)


@needs_numpy
def test_exactly_as_many_active_bits_as_the_hash_pays_for_is_allowed():
    fp = read_smiles('CC(C)C(=O)Nc1ccccc1').morgan_fingerprint(length=1024, number_active_bits=6)
    assert fp.sum() > 0


def test_the_folding_validation_runs_before_the_walk():
    # a bad length is cheap to notice, so it is noticed first -- not after a fused polycyclic walk
    with raises(ValueError, match='power of two'):
        read_smiles('CCO').morgan_bit_set(0, 4, 1000, 2)


# --- linear paths --------------------------------------------------------------------------------

@needs_numpy
def test_linear_radius_one_is_one_fragment_per_atom():
    mol = read_smiles('CC(=O)OCC')
    assert sum(mol.linear_hash_counts(1, 1).values()) == mol.atom_count


@needs_numpy
def test_linear_radius_two_counts_every_bond_exactly_once():
    # a path of two atoms IS a bond, and each is counted once rather than once per direction
    for smi in ('CCO', 'C1CCCCC1', 'CC(=O)OCC', 'c1ccccc1O'):
        mol = read_smiles(smi)
        assert sum(mol.linear_hash_counts(2, 2).values()) == mol.bond_count, smi


@needs_numpy
def test_linear_reading_a_path_from_either_end_gives_one_hash():
    a = read_smiles('CCO')
    b = read_smiles('OCC')
    assert a == b
    assert a.linear_hash_counts(1, 3) == b.linear_hash_counts(1, 3)
    # mixed bond orders: CC=O written head-to-tail as two components of one molecule.
    # Each is the same fragment, so path-hash must give ONE key with count TWO --
    # this fails if the reversal only swaps the label slots and leaves the order slots alone.
    mixed = read_smiles('CC=O.O=CC')
    counts = mixed.linear_hash_counts(3, 3)
    assert list(counts.values()) == [2], counts


@needs_numpy
def test_linear_pools_two_equivalent_bonds_into_one_key():
    # propane's two C-C bonds are the same fragment, so one key of count two
    counts = read_smiles('CCC').linear_hash_counts(2, 2)
    assert list(counts.values()) == [2]


@needs_numpy
def test_linear_separates_two_bonds_that_differ_only_in_order():
    # the C-C and the C=O of acetaldehyde are different fragments
    counts = read_smiles('CC=O').linear_hash_counts(2, 2)
    assert sorted(counts.values()) == [1, 1]


@needs_numpy
def test_linear_growing_the_length_only_adds_fragments():
    mol = read_smiles('CC(C)C(=O)Nc1ccccc1')
    assert mol.linear_hash_set(1, 2) < mol.linear_hash_set(1, 4)


@needs_numpy
def test_linear_is_invariant_to_the_written_atom_order():
    a = read_smiles('CC(=O)OCC')
    b = read_smiles('CCOC(C)=O')
    assert a == b
    assert a.linear_hash_counts() == b.linear_hash_counts()


@needs_numpy
def test_linear_walks_a_ring_without_revisiting_an_atom():
    # cyclopropane has three atoms, so a path of four atoms does not exist: a path is SIMPLE
    mol = read_smiles('C1CC1')
    assert mol.linear_hash_counts(4, 4) == {}
    assert sum(mol.linear_hash_counts(3, 3).values()) == 3


@needs_numpy
def test_linear_of_an_empty_molecule_is_empty():
    assert MoleculeContainer().linear_hash_counts() == {}


@needs_numpy
def test_linear_of_one_atom_has_only_the_one_atom_path():
    mol = read_smiles('[Na+]')
    assert sum(mol.linear_hash_counts(1, 4).values()) == 1


@needs_numpy
def test_linear_huge_max_radius_on_tiny_molecule_matches_diameter_result():
    # ethanol has three heavy atoms, so no path can be longer than 3 atoms.  A max_radius far
    # beyond the molecule's diameter must not allocate a giant scratch -- it must clamp to n_atoms
    # and return the same result as a max_radius that already covers the whole molecule.
    mol = read_smiles('CCO')   # ethanol: C-C-O, diameter = 2 bonds = 3 atoms
    normal = mol.linear_hash_counts(1, 3)     # covers everything; any larger max_radius adds nothing
    huge = mol.linear_hash_counts(1, 10 ** 6)
    assert huge == normal


@needs_numpy
def test_linear_and_morgan_disagree_because_they_enumerate_different_things():
    mol = read_smiles('CC(C)C(=O)Nc1ccccc1')
    assert mol.linear_hash_set() != mol.morgan_hash_set()


@needs_numpy
def test_linear_separates_two_isomers_that_share_every_bond_type():
    assert read_smiles('CCCCO').linear_hash_set(1, 4) != read_smiles('CC(C)CO').linear_hash_set(1, 4)


@needs_numpy
def test_linear_reads_a_caller_supplied_label_vector():
    mol = read_smiles('CCO')
    flat = mol.atom_invariants()
    flat[:] = 1
    assert list(mol.linear_hash_counts(1, 1, invariants=flat).values()) == [3]


@needs_numpy
def test_linear_folded_spellings_agree_with_each_other():
    mol = read_smiles('CC(C)C(=O)Nc1ccccc1')
    bits = mol.linear_bit_set()
    assert set(mol.linear_fingerprint().nonzero()[0].tolist()) == bits
    assert set(mol.linear_count_vector().nonzero()[0].tolist()) == bits
    assert mol.linear_count_vector().sum() == sum(mol.linear_hash_counts().values()) * 2


def test_linear_takes_the_same_validation_as_morgan():
    mol = read_smiles('CCO')
    with raises(ValueError, match='min_radius'):
        mol.linear_hash_counts(0, 4)
    with raises(ValueError, match='max_radius'):
        mol.linear_hash_counts(4, 2)
    with raises(ValueError, match='power of two'):
        mol.linear_fingerprint(length=1000)
    with raises(ValueError, match='64 bits'):
        mol.linear_fingerprint(length=1024, number_active_bits=7)


def test_linear_folding_validation_runs_before_the_walk():
    # a bad length is cheap to notice, so it is noticed first -- not after an exponential walk
    mol = read_smiles('c1ccccc1')
    with raises(ValueError, match='power of two'):
        mol.linear_bit_set(1, 20, 1000, 2)
    with raises(ValueError, match='power of two'):
        mol.linear_count_vector(1, 20, 1000, 2)


# --- the surface itself --------------------------------------------------------------------------

def test_both_families_offer_the_same_five_spellings():
    # derive the two families' suffixes from the class itself so an asymmetric addition fails by name
    _EXPECTED = {'hash_set', 'hash_counts', 'bit_set', 'fingerprint', 'count_vector'}
    morgan_suffixes = {n[len('morgan_'):] for n in dir(MoleculeContainer) if n.startswith('morgan_')}
    linear_suffixes = {n[len('linear_'):] for n in dir(MoleculeContainer) if n.startswith('linear_')}
    assert morgan_suffixes == linear_suffixes, (
        f'families are asymmetric: morgan-only={morgan_suffixes - linear_suffixes}, '
        f'linear-only={linear_suffixes - morgan_suffixes}'
    )
    assert _EXPECTED <= morgan_suffixes, f'missing spellings: {_EXPECTED - morgan_suffixes}'


def test_the_dropped_chython_two_spellings_stay_dropped():
    # the four *_smiles methods returned a human explanation of a bit and had no call site anywhere;
    # `number_bit_pairs` approximated counting and is replaced by the count_vector spellings
    for gone in ('morgan_hash_smiles', 'morgan_smiles_hash', 'linear_hash_smiles',
                 'linear_smiles_hash'):
        assert not hasattr(MoleculeContainer, gone), gone
    with raises(TypeError):
        read_smiles('CCO').linear_fingerprint(number_bit_pairs=4)


def test_features_is_not_offered_before_the_table_that_gives_it_meaning():
    # `features=True` is sugar for `invariants=pharmacophore_invariants()` and that table is F3's; a
    # keyword whose only behaviour is to raise is worse than one that arrives with its table.
    # Loop over all ten methods: adding features= to linear_fingerprint alone must fail this gate.
    _ALL_TEN = [
        'morgan_hash_set', 'morgan_hash_counts', 'morgan_bit_set',
        'morgan_fingerprint', 'morgan_count_vector',
        'linear_hash_set', 'linear_hash_counts', 'linear_bit_set',
        'linear_fingerprint', 'linear_count_vector',
    ]
    mol = read_smiles('CCO')
    for name in _ALL_TEN:
        with raises(TypeError, match='features'):
            getattr(mol, name)(features=True)


# --- the R marker --------------------------------------------------------------------------------

_R_SURFACES = ('linear_fingerprint', 'morgan_fingerprint', 'linear_bit_set', 'linear_hash_set',
               'morgan_bit_set', 'morgan_hash_set')


def _screen(mol, name):
    """One surface's answer, as something comparable: an array by its bytes, a set as itself."""
    out = getattr(mol, name)()
    return out.tobytes() if hasattr(out, 'tobytes') else out


@needs_numpy
@mark.parametrize('name', _R_SURFACES)
def test_a_marker_is_not_a_carbon_to_a_fingerprint(name):
    """Rule 3, the fingerprint half: an R is not carbon for identity, on every screening surface."""
    assert _screen(read_smiles('[R1]c1ccccc1'), name) != _screen(read_smiles('Cc1ccccc1'), name)


@needs_numpy
@mark.parametrize('name', _R_SURFACES)
def test_the_index_is_invisible_to_a_fingerprint(name):
    """The index does not reach the atom invariant, which is what a screen can carry.

    A fingerprint is a lossy prefilter over element, charge, isotope, radical and connectivity; two
    fragments differing only in which numbered handle they expose are the same substructure, and a
    screen that separated them would reject a genuine superstructure. The canonical form is where the
    index is a distinction.
    """
    assert _screen(read_smiles('[R1]c1ccccc1'), name) == _screen(read_smiles('[R2]c1ccccc1'), name)
    assert read_smiles('[R1]c1ccccc1') != read_smiles('[R2]c1ccccc1')
