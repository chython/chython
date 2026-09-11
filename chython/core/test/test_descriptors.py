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
"""The graph descriptors: composition counts, ring system counts, and the topological indices.

EVERY EXPECTED VALUE HERE IS EITHER QUOTED FROM THE PAPER THAT DEFINED THE INDEX OR DERIVED BY HAND
FROM THE FORMULA IN THE DOCSTRING, and each one says which in a comment.  No number in this file came
from another toolkit -- an oracle would pin chython to a foreign reading of Kier and Hall rather than
to Kier and Hall, and the epic forbids it.

The molecules are textbook: n-butane, isobutane, neopentane, benzene, cyclohexane, naphthalene,
toluene, phenol, norbornane.  That is not modesty, it is what makes the arithmetic checkable in the
comment beside the assertion.

A PINNED LITERAL CARRIES `rel=PIN` AND NOT `approx`'s DEFAULT, and the reason is a defect this file
shipped twice.  Many assertions come in pairs: a closed form, then the decimal it evaluates to.  The
closed form cannot catch a drift, because it recomputes with the same arithmetic the code uses -- that
is the literal's whole job.  But `approx`'s default relative tolerance is **1e-6**, and two literals
here were transcribed wrong by hand at the eighth and ninth significant digit and passed anyway,
pinning nothing at all: a Bertz CT of 28.364527700498157 against the true 28.364527976600278, and a
kappa of 5.482234032309491 against 5.482229779997874.  So every hand-transcribed decimal is asserted at
`rel=PIN` (1e-12), loose enough to survive a last-bit difference in `libm`'s `log2` across the two
architectures this wheel builds for, and tight enough that a wrong transcription fails immediately.
A value taken straight from a paper's own table keeps a `round(...)` at the paper's precision instead --
that is a different assertion and says so.

THE GRAPH IS THE GRAPH AS STORED (see `_descriptors.pxi`): a dative bond is an edge and an explicit
hydrogen is a vertex, so `read_smiles('C')` and `read_smiles('[H]C([H])([H])[H]')` are two different
graphs here and get different numbers.  Every molecule below is written with implicit hydrogens, which
is the hydrogen-suppressed graph the classical indices are defined on.
"""
from importlib.util import find_spec
from math import log2, sqrt
from pytest import approx, mark, raises

from chython.core import H_UNKNOWN, MoleculeContainer, read_smiles


#: the relative tolerance for a hand-transcribed decimal -- see the module docstring for why it is not
#: `approx`'s default 1e-6
PIN = 1e-12

# THE DISTANCE-DERIVED DESCRIPTORS ARE THE NUMPY HALF OF THIS FILE, and only they: numpy is an optional
# dependency (`chython[ml]`), and `distance_matrix` is the core's only shortest-path code, so
# `eccentricities`, `wiener_index`, `graph_radius`, `graph_diameter`, `balaban_j` and the two `estate_*`
# all read an array and everything else here -- the counts, the degree indices, chi, kappa,
# hall_kier_alpha, bertz_ct -- does not.  Which is which was measured against an install with numpy
# blocked, not reasoned about.  `estate_intrinsic_states` reads no distance and is marked anyway: its
# ANSWER is an array, which is the second way into this half and the one a reader has to be told about.
#
# `find_spec` rather than `importorskip`, so collection does not import numpy at all -- the same
# reasoning as `interop/test/conftest.py` gives for the optional toolkits.
needs_numpy = mark.skipif(find_spec('numpy') is None,
                          reason='numpy is not installed; the distance matrix these read is an array')


# --- composition counts ---------------------------------------------------------------------------


def test_carbon_count_is_the_element_bucket():
    assert read_smiles('Cc1ccccc1').carbon_count == 7            # toluene
    assert read_smiles('CC(=O)O').carbon_count == 2              # acetic acid
    assert read_smiles('O').carbon_count == 0                    # water
    assert MoleculeContainer().carbon_count == 0


def test_carbon_sp3_count_reads_the_stored_hybridization():
    """z == 1 and nothing else, so an aromatic carbon (z == 4) is not sp3 and neither is a nitrile's.

    chython 2 spelled this `a == C and a.hybridization == 1` and V3's z scale agrees with V2 on 1.
    """
    assert read_smiles('Cc1ccccc1').carbon_sp3_count == 1        # toluene: the methyl only
    assert read_smiles('C1CCCCC1').carbon_sp3_count == 6         # cyclohexane
    assert read_smiles('c1ccccc1').carbon_sp3_count == 0         # benzene
    assert read_smiles('CC#N').carbon_sp3_count == 1             # acetonitrile: methyl yes, nitrile no
    assert read_smiles('CC=C').carbon_sp3_count == 1             # propene: one sp3, two sp2


def test_carbon_sp3_fraction_is_the_ratio():
    assert read_smiles('Cc1ccccc1').carbon_sp3_fraction == approx(1 / 7)
    assert read_smiles('C1CCCCC1').carbon_sp3_fraction == 1.0
    assert read_smiles('c1ccccc1').carbon_sp3_fraction == 0.0


def test_carbon_sp3_fraction_of_a_molecule_with_no_carbon_is_zero():
    """0.0 and not nan and not a refusal -- chython 2 answered 0. and the number feeds a vector
    where a nan poisons the whole row."""
    assert read_smiles('O').carbon_sp3_fraction == 0.0
    assert read_smiles('[Na+].[Cl-]').carbon_sp3_fraction == 0.0
    assert MoleculeContainer().carbon_sp3_fraction == 0.0


def test_heteroatoms_count_is_every_atom_that_is_not_carbon_and_not_hydrogen():
    """A COUNT OF ATOMS, and deliberately not a sum of the per-atom `heteroatoms_of`.

    `heteroatoms_of(n)` counts an atom's heteroatom NEIGHBOURS, so summing it over the molecule counts
    each heteroatom once per bond it has -- a different question with a different answer.
    """
    assert read_smiles('CC(=O)O').heteroatoms_count == 2         # acetic acid: two oxygens
    assert read_smiles('c1ccncc1').heteroatoms_count == 1        # pyridine
    assert read_smiles('c1ccccc1').heteroatoms_count == 0
    assert read_smiles('[Na+].[Cl-]').heteroatoms_count == 2
    assert MoleculeContainer().heteroatoms_count == 0


def test_valence_electrons_count_sums_the_group_number_less_the_charge_plus_the_hydrogens():
    """Zv - charge + implicit hydrogens, over every stored atom.

    Derived by hand: benzene is six carbons at 4 + 1 = 5; toluene is 4 + 3 for the methyl, 4 + 0 for
    the ipso carbon and five ring CH at 4 + 1; nitrate is N at 5 - 1 = 4, the neutral O at 6, and two
    O- at 6 + 1 = 7, which is NO3-'s 24 the textbook way (5 + 18 + 1 for the charge).
    """
    assert read_smiles('c1ccccc1').valence_electrons_count == 30
    assert read_smiles('Cc1ccccc1').valence_electrons_count == 36
    assert read_smiles('C').valence_electrons_count == 8         # methane
    assert read_smiles('O').valence_electrons_count == 8         # water
    assert read_smiles('[O-][N+](=O)[O-]').valence_electrons_count == 24
    assert MoleculeContainer().valence_electrons_count == 0


def test_valence_electrons_count_does_not_double_count_an_explicit_hydrogen():
    """The one arithmetic that has to work in both spellings, and the reason `explicit_h` is NOT in
    the sum: an explicit hydrogen is a vertex contributing its own electron already."""
    implicit = read_smiles('C')
    explicit = read_smiles('[H]C([H])([H])[H]')
    assert implicit.valence_electrons_count == 8
    assert explicit.valence_electrons_count == 8


def test_valence_electrons_count_refuses_an_f_block_atom():
    """The reserved unknown reaches the surface as a refusal, never as a zero."""
    with raises(ValueError, match='no valence electron count'):
        read_smiles('[Ce]').valence_electrons_count


def test_valence_electrons_count_refuses_an_unknown_hydrogen_count():
    """A sum with an unknown term is unknown -- the same answer `total_h_of` gives per atom, and the
    message names the repair rather than guessing at zero."""
    m = MoleculeContainer()
    with m.edit():
        m.add_atom('C', implicit_h=H_UNKNOWN)
    with raises(ValueError, match='calc_implicit'):
        m.valence_electrons_count


# --- ring systems ---------------------------------------------------------------------------------


def test_ring_classes_of_the_textbook_rings():
    """Every count is a classification of the MINIMUM CYCLE BASIS -- `rings`, the same set `sssr`
    and `aromatic_rings` report.  Aromatic means every bond in the ring is stored order 4; saturated
    means every bond is stored order 1; heterocyclic means the ring holds an atom that is neither
    carbon nor hydrogen.
    """
    benzene = read_smiles('c1ccccc1')
    assert benzene.aromatic_rings_count == 1
    assert benzene.aliphatic_rings_count == 0
    assert benzene.saturated_rings_count == 0
    assert benzene.heterocycles_count == 0
    assert benzene.aromatic_heterocycles_count == 0

    cyclohexane = read_smiles('C1CCCCC1')
    assert cyclohexane.aromatic_rings_count == 0
    assert cyclohexane.aliphatic_rings_count == 1
    assert cyclohexane.saturated_rings_count == 1

    cyclohexene = read_smiles('C1CCCC=C1')
    assert cyclohexene.aliphatic_rings_count == 1
    assert cyclohexene.saturated_rings_count == 0     # one double bond is enough

    pyridine = read_smiles('c1ccncc1')
    assert pyridine.aromatic_rings_count == 1
    assert pyridine.heterocycles_count == 1
    assert pyridine.aromatic_heterocycles_count == 1

    morpholine = read_smiles('C1COCCN1')
    assert morpholine.aliphatic_rings_count == 1
    assert morpholine.saturated_rings_count == 1
    assert morpholine.heterocycles_count == 1
    assert morpholine.aromatic_heterocycles_count == 0


def test_ring_classes_of_the_fused_pairs():
    quinoline = read_smiles('c1ccc2ncccc2c1')
    assert quinoline.aromatic_rings_count == 2
    assert quinoline.heterocycles_count == 1          # the pyridine ring only
    assert quinoline.aromatic_heterocycles_count == 1

    tetralin = read_smiles('c1ccc2c(c1)CCCC2')
    assert tetralin.aromatic_rings_count == 1
    assert tetralin.aliphatic_rings_count == 1
    # the carbocycle shares an ORDER-4 bond with the arene, so it is not saturated either
    assert tetralin.saturated_rings_count == 0

    decalin = read_smiles('C1CCC2CCCCC2C1')
    assert decalin.aliphatic_rings_count == 2
    assert decalin.saturated_rings_count == 2


def test_aromatic_and_aliphatic_partition_the_ring_basis():
    """The internal consistency the two names owe each other: a ring is one or the other, never
    both and never neither, so they sum to `rings_count` for every molecule."""
    for smiles in ('c1ccccc1', 'C1CCCCC1', 'c1ccc2ccccc2c1', 'c1ccc2c(c1)CCCC2', 'CCO',
                   'C1CC2CCC1C2', 'c1ccc(-c2ccccc2)cc1'):
        m = read_smiles(smiles)
        assert m.aromatic_rings_count + m.aliphatic_rings_count == m.rings_count, smiles


def test_aromatic_rings_count_agrees_with_the_shipped_aromatic_rings():
    """`aromatic_rings` already ships and is a Python-level filter over `rings`; this count must be
    its length or one of the two is wrong."""
    for smiles in ('c1ccccc1', 'c1ccc2ccccc2c1', 'c1ccc2c(c1)CCCC2', 'C1CCCCC1', 'c1ccncc1'):
        m = read_smiles(smiles)
        assert m.aromatic_rings_count == len(m.aromatic_rings), smiles


def test_spiro_atoms_count_is_a_ring_pair_sharing_exactly_one_atom():
    """spiro[4.5]decane is a cyclohexane and a cyclopentane meeting at one carbon.

    Fused rings share two atoms and a bond, so naphthalene has no spiro atom, and an isolated ring
    has no pair to share with.
    """
    assert read_smiles('C1CCC2(CC1)CCCC2').spiro_atoms_count == 1   # spiro[4.5]decane
    assert read_smiles('c1ccc2ccccc2c1').spiro_atoms_count == 0     # naphthalene
    assert read_smiles('C1CCCCC1').spiro_atoms_count == 0
    assert read_smiles('CCO').spiro_atoms_count == 0


def test_bridgehead_atoms_count_needs_a_shared_path_and_three_ring_bonds():
    """Norbornane's two bridgeheads, and nothing else in the textbook set.

    The two five-rings of the basis share three atoms and TWO bonds -- a path, which is what makes
    the system bridged rather than merely fused.  The middle atom of that path is the one-carbon
    bridge and is NOT a bridgehead: it carries two ring bonds, where a bridgehead carries three.
    That second clause is what separates norbornane's 2 from an unqualified 3.

    Naphthalene and decalin share ONE bond, so they are fused and have no bridgehead at all, even
    though their two fusion atoms do carry three ring bonds each.
    """
    assert read_smiles('C1CC2CCC1C2').bridgehead_atoms_count == 2   # norbornane
    assert read_smiles('c1ccc2ccccc2c1').bridgehead_atoms_count == 0
    assert read_smiles('C1CCC2CCCCC2C1').bridgehead_atoms_count == 0
    assert read_smiles('C1CCC2(CC1)CCCC2').bridgehead_atoms_count == 0


def test_fused_ring_systems_count_counts_ring_SYSTEMS():
    """Connected components of the subgraph of ring bonds.

    AN ISOLATED RING IS ONE SYSTEM.  The name says "fused" because a system may be fused, not
    because fusion is required -- benzene answers 1, and a caller who wants "systems of more than one
    ring" subtracts the count of systems that are a single ring, which `rings_count` and this number
    do not give on their own.  A spiro atom joins its two rings into one system: the bonds of both
    meet at it.
    """
    assert read_smiles('c1ccccc1').fused_ring_systems_count == 1
    assert read_smiles('c1ccc2ccccc2c1').fused_ring_systems_count == 1     # naphthalene
    assert read_smiles('c1ccc(-c2ccccc2)cc1').fused_ring_systems_count == 2  # biphenyl
    assert read_smiles('C1CCC2CCCCC2C1').fused_ring_systems_count == 1     # decalin
    assert read_smiles('C1CCC2(CC1)CCCC2').fused_ring_systems_count == 1   # spiro[4.5]decane
    assert read_smiles('C1CC2CCC1C2').fused_ring_systems_count == 1        # norbornane
    assert read_smiles('CCO').fused_ring_systems_count == 0


def test_the_ring_counts_of_an_acyclic_and_an_empty_molecule_are_zero():
    for m in (read_smiles('CCCC'), read_smiles('[Na+].[Cl-]'), MoleculeContainer()):
        assert m.aromatic_rings_count == 0
        assert m.aliphatic_rings_count == 0
        assert m.saturated_rings_count == 0
        assert m.heterocycles_count == 0
        assert m.aromatic_heterocycles_count == 0
        assert m.spiro_atoms_count == 0
        assert m.bridgehead_atoms_count == 0
        assert m.fused_ring_systems_count == 0


def test_the_ring_counts_are_additive_over_components():
    """A salt of two ring systems is the sum of its parts, and `connected_components_count` is the
    number that says the molecule was disconnected in the first place."""
    salt = read_smiles('c1ccccc1.C1CCCCC1')
    assert salt.connected_components_count == 2
    assert salt.aromatic_rings_count == 1
    assert salt.aliphatic_rings_count == 1
    assert salt.saturated_rings_count == 1
    assert salt.fused_ring_systems_count == 2


# --- invariance -----------------------------------------------------------------------------------


DESCRIPTORS = (
    ('carbon_count', lambda m: m.carbon_count),
    ('carbon_sp3_count', lambda m: m.carbon_sp3_count),
    ('carbon_sp3_fraction', lambda m: m.carbon_sp3_fraction),
    ('heteroatoms_count', lambda m: m.heteroatoms_count),
    ('valence_electrons_count', lambda m: m.valence_electrons_count),
    ('aromatic_rings_count', lambda m: m.aromatic_rings_count),
    ('aliphatic_rings_count', lambda m: m.aliphatic_rings_count),
    ('saturated_rings_count', lambda m: m.saturated_rings_count),
    ('heterocycles_count', lambda m: m.heterocycles_count),
    ('aromatic_heterocycles_count', lambda m: m.aromatic_heterocycles_count),
    ('spiro_atoms_count', lambda m: m.spiro_atoms_count),
    ('bridgehead_atoms_count', lambda m: m.bridgehead_atoms_count),
    ('fused_ring_systems_count', lambda m: m.fused_ring_systems_count),
    ('randic_index', lambda m: m.randic_index),
    ('bertz_ct', lambda m: m.bertz_ct),
    ('hall_kier_alpha', lambda m: m.hall_kier_alpha),
    # the two that shipped before F2, swept with the rest because the sweep is about the surface and not
    # about which commit landed it
    ('is_radical', lambda m: m.is_radical),
    ('len(aromatic_rings)', lambda m: len(m.aromatic_rings)),
) + tuple(
    # order=order default is load-bearing: a closure over the loop variable would make every
    # entry read the last order and the sweep would pass while testing one case
    ('chi(%d)' % order, lambda m, order=order: m.chi(order)) for order in range(5)
) + tuple(
    # order=order default is load-bearing: same reason as above
    ('chi(%d, valence=True)' % order, lambda m, order=order: m.chi(order, valence=True))
    for order in range(5)
) + tuple(
    # order=order default is load-bearing: same reason as above
    ('kappa(%d)' % order, lambda m, order=order: m.kappa(order)) for order in (1, 2, 3)
) + tuple(
    # order=order default is load-bearing: same reason as above
    ('kappa(%d, alpha=True)' % order, lambda m, order=order: m.kappa(order, alpha=True))
    for order in (1, 2, 3)
) + tuple(
    # order=order default is load-bearing: same reason as above
    ('zagreb_index(%d)' % order, lambda m, order=order: m.zagreb_index(order)) for order in (1, 2)
)
# balaban_j stays out: it raises on a disconnected molecule and the corpus's last row is
# 'c1ccccc1.C1CCCCC1'; its own invariance test (test_balaban_j_is_invariant_under_renumbering)
# covers connected molecules.
# eccentricities stays out: the generic loop compares scalars and cannot compare sorted multisets;
# test_eccentricities_invariant_under_renumbering handles that.

#: The three that read `distance_matrix`, held in their own table so the sweep above still runs on an
#: install without numpy.  THE SPLIT IS BY WHAT THE DESCRIPTOR READS AND NOTHING ELSE -- they are swept
#: over the same corpus by the same body, one test below the other, and the only difference between the
#: two is a `needs_numpy`.  Folding them back in would cost the other twenty-nine their sweep on a
#: minimal install; skipping the whole sweep instead of splitting it would cost the same twenty-nine,
#: which is the reverse of the narrowest honest guard.
DISTANCE_DESCRIPTORS = (
    ('wiener_index', lambda m: m.wiener_index),
    ('graph_radius', lambda m: m.graph_radius),
    ('graph_diameter', lambda m: m.graph_diameter),
)

#: The corpus both renumbering sweeps read.  One list, so the two cannot drift apart.
RENUMBERING_CORPUS = [
    'c1ccc2[nH]ccc2c1',                              # indole -- a fused aromatic heterocycle
    'CC(=O)Oc1ccccc1C(=O)O',                         # aspirin
    'C1CCC2(CC1)CCCCC2',                             # spiro[5.5]undecane -- moves spiro_atoms_count
    'C1CC2CCC1C2',                                   # norbornane -- moves bridgehead_atoms_count
    'Cn1cnc2c1c(=O)n(C)c(=O)n2C',                    # caffeine
    'c1ccccc1.C1CCCCC1',                             # two components, so two ring systems
]


def _rebuilt_backwards(molecule):
    """The same graph with its atoms appended in reverse, so every stable id and arena slot moves."""
    out = MoleculeContainer()
    ids = {}
    with out.edit():
        for atom in reversed(list(molecule.atoms())):
            ids[atom.n] = out.add_atom(atom.atomic_symbol, charge=atom.charge,
                                       isotope=atom.isotope, radical=atom.is_radical,
                                       implicit_h=atom.implicit_h)
        for bond in molecule.bonds():
            out.add_bond(ids[bond.n], ids[bond.m], bond.order)
    return out


@mark.parametrize('smiles', RENUMBERING_CORPUS)
def test_every_descriptor_is_invariant_under_renumbering(smiles):
    """RENUMBERING INVARIANCE IS THE STRONGEST TEST THESE DESCRIPTORS HAVE, because a graph descriptor
    that reads an arena slot instead of the graph passes every value test and fails this one.

    `_rebuilt_backwards` reverses the append order, so no atom keeps its stable id or its index, and
    the element buckets, the CSR rows and the ring basis are all rebuilt from a different starting
    point.  A count that came out of `sssr`'s choice among equal-size cycles, or out of a union-find
    seeded by index order, is what this catches.
    """
    molecule = read_smiles(smiles)
    other = _rebuilt_backwards(molecule)
    assert other.connected_components_count == molecule.connected_components_count
    for name, getter in DESCRIPTORS:
        assert getter(other) == approx(getter(molecule)), name


@needs_numpy
@mark.parametrize('smiles', RENUMBERING_CORPUS)
def test_every_distance_descriptor_is_invariant_under_renumbering(smiles):
    """The same sweep for the three that read `distance_matrix`, and the same claim about them.

    A separate test only because numpy is optional -- see `DISTANCE_DESCRIPTORS`.  The body is the one
    above with the other table; the corpus is the same list, so a row added there is swept by both.
    """
    molecule = read_smiles(smiles)
    other = _rebuilt_backwards(molecule)
    assert other.connected_components_count == molecule.connected_components_count
    for name, getter in DISTANCE_DESCRIPTORS:
        assert getter(other) == approx(getter(molecule)), name


# --- distances ------------------------------------------------------------------------------------


@needs_numpy
def test_wiener_index_of_the_pentanes_and_butanes():
    """W = sum of d(i, j) over unordered pairs -- Wiener, JACS 69 (1947) 17, where he calls it the
    path number and tabulates exactly these.

    Derived by hand for the path graphs: n-butane's six pairs are 1 + 2 + 3 + 1 + 2 + 1 = 10, and
    n-pentane's ten are (1+2+3+4) + (1+2+3) + (1+2) + 1 = 20.  Isobutane is the star K(1,3): three
    pairs at distance 1 through the centre and three leaf pairs at 2, so 3 + 6 = 9.  Neopentane is
    K(1,4): 4 + 6*2 = 16.  Isopentane CC(C)CC: (1+2+2+3) + (1+1+2) + (2+3) + 1 = 18.
    """
    assert read_smiles('CCCC').wiener_index == 10
    assert read_smiles('CCCCC').wiener_index == 20
    assert read_smiles('CC(C)C').wiener_index == 9
    assert read_smiles('CC(C)(C)C').wiener_index == 16
    assert read_smiles('CC(C)CC').wiener_index == 18


@needs_numpy
def test_wiener_index_of_the_rings_and_the_fused_pair():
    """Benzene: every vertex sees 1, 1, 2, 2, 3 = 9, so the ordered sum is 54 and W is 27.  A
    cyclohexane is the same graph and the same 27 -- W reads distances, not bond orders.

    Toluene adds a methyl whose distances to the ring are 1 + 2 + 2 + 3 + 3 + 4 = 15, so 27 + 15 = 42;
    phenol is the same skeleton and the same 42.

    Naphthalene, derived by BFS from one atom of each of its three symmetry classes: an alpha carbon
    sums 2*1 + 3*2 + 3*3 + 4 = 21, a beta carbon 2 + 4 + 6 + 8 + 5 = 25, a fusion carbon
    3 + 8 + 6 = 17.  Four alphas, four betas, two fusions: (4*21 + 4*25 + 2*17) / 2 = 218 / 2 = 109.
    """
    assert read_smiles('c1ccccc1').wiener_index == 27
    assert read_smiles('C1CCCCC1').wiener_index == 27
    assert read_smiles('Cc1ccccc1').wiener_index == 42
    assert read_smiles('Oc1ccccc1').wiener_index == 42
    assert read_smiles('c1ccc2ccccc2c1').wiener_index == 109


@needs_numpy
def test_wiener_index_skips_a_pair_with_no_path():
    """A -1 in the distance matrix is not summed, which makes W additive over components -- two
    butanes are 10 + 10 and not "10 + 10 + something for the pairs that do not exist".

    That is the honest reading and it is stated rather than assumed: the alternative, an infinite or
    a sentinel-laden W, is a number no caller can use.
    """
    assert read_smiles('CCCC.CCCC').wiener_index == 20
    assert read_smiles('[Na+].[Cl-]').wiener_index == 0
    assert read_smiles('c1ccccc1.[Na+]').wiener_index == 27
    assert read_smiles('C').wiener_index == 0
    assert MoleculeContainer().wiener_index == 0


@needs_numpy
def test_eccentricities_are_within_the_component_and_indexed_like_distance_matrix():
    """Eccentricity is the largest distance from an atom TO AN ATOM IT CAN REACH.

    Row order is `distance_matrix`'s: entry i belongs to `atom_numbers[i]`.  n-butane's two ends see 3
    and its two middles see 2.
    """
    butane = read_smiles('CCCC')
    assert list(butane.eccentricities()) == [3, 2, 2, 3]
    assert butane.eccentricities().dtype == 'int32'
    assert butane.eccentricities().shape == (4,)

    benzene = read_smiles('c1ccccc1')
    assert list(benzene.eccentricities()) == [3, 3, 3, 3, 3, 3]


@needs_numpy
def test_eccentricity_of_a_lone_counterion_is_zero():
    """An atom that can reach nothing has eccentricity 0, because -1 is not a distance and is not
    folded into a maximum.  It is the same 0 a single-atom molecule gets, which is why
    `connected_components_count` and not the eccentricity is the way to learn a molecule is a salt.
    """
    salt = read_smiles('c1ccccc1.[Na+]')
    assert list(salt.eccentricities()) == [3, 3, 3, 3, 3, 3, 0]
    assert MoleculeContainer().eccentricities().shape == (0,)


@needs_numpy
def test_graph_radius_and_diameter_are_the_min_and_max_eccentricity():
    """Toluene's para carbon and its methyl both see 4, everything else 3 -- so radius 3, diameter 4.
    Naphthalene: a fusion carbon sees 3, a beta carbon 5.
    """
    butane = read_smiles('CCCC')
    assert butane.graph_radius == 2
    assert butane.graph_diameter == 3

    benzene = read_smiles('c1ccccc1')
    assert benzene.graph_radius == 3
    assert benzene.graph_diameter == 3

    toluene = read_smiles('Cc1ccccc1')
    assert toluene.graph_radius == 3
    assert toluene.graph_diameter == 4

    naphthalene = read_smiles('c1ccc2ccccc2c1')
    assert naphthalene.graph_radius == 3
    assert naphthalene.graph_diameter == 5


@needs_numpy
def test_graph_radius_of_a_disconnected_molecule_is_zero_and_the_diameter_is_the_widest_component():
    """The consequence of an eccentricity being within-component, spelled out: a lone counterion has
    eccentricity 0, so it is the minimum, so the radius of any salt containing one is 0.  The diameter
    is the widest component's diameter, which is a usable number.
    """
    salt = read_smiles('c1ccccc1.[Na+]')
    assert salt.graph_radius == 0
    assert salt.graph_diameter == 3
    assert read_smiles('CCCC.CCCCCCC').graph_diameter == 6
    assert MoleculeContainer().graph_radius == 0
    assert MoleculeContainer().graph_diameter == 0


@needs_numpy
def test_radius_and_diameter_agree_with_the_eccentricity_vector():
    """Internal consistency, over molecules whose numbers are asserted individually above and a few
    that are not."""
    for smiles in ('CCCC', 'c1ccccc1', 'Cc1ccccc1', 'c1ccc2ccccc2c1', 'C1CC2CCC1C2',
                   'c1ccc(-c2ccccc2)cc1', 'CC(=O)Oc1ccccc1C(=O)O'):
        m = read_smiles(smiles)
        ecc = list(m.eccentricities())
        assert m.graph_radius == min(ecc), smiles
        assert m.graph_diameter == max(ecc), smiles


@needs_numpy
def test_eccentricities_invariant_under_renumbering():
    """The multiset of eccentricities is graph-invariant: every atom's max distance depends only on
    the graph topology, not on the arena slot order.  Per-index values shift, so we compare sorted.
    """
    for smiles in ('CCCC', 'c1ccccc1', 'Cc1ccccc1', 'c1ccc2ccccc2c1', 'c1ccccc1.[Na+]'):
        molecule = read_smiles(smiles)
        other = _rebuilt_backwards(molecule)
        assert sorted(molecule.eccentricities()) == sorted(other.eccentricities()), smiles


# --- degree indices -------------------------------------------------------------------------------


def test_first_zagreb_index_is_the_sum_of_squared_degrees():
    """M1 = sum over atoms of deg(v)^2.  Gutman and Trinajstic, Chem. Phys. Lett. 17 (1972) 535.

    By hand: n-butane's degrees are 1, 2, 2, 1 so M1 = 1 + 4 + 4 + 1 = 10; isobutane is the star with
    degrees 3, 1, 1, 1 so 9 + 3 = 12; neopentane 16 + 4 = 20; benzene is six degree-2 atoms, 24.
    Toluene: one degree-1 methyl, one degree-3 ipso carbon and five degree-2 carbons, 1 + 9 + 20 = 30.
    """
    assert read_smiles('CCCC').zagreb_index() == 10
    assert read_smiles('CC(C)C').zagreb_index() == 12
    assert read_smiles('CC(C)(C)C').zagreb_index() == 20
    assert read_smiles('c1ccccc1').zagreb_index() == 24
    assert read_smiles('Cc1ccccc1').zagreb_index() == 30
    assert read_smiles('C').zagreb_index() == 0
    assert MoleculeContainer().zagreb_index() == 0


def test_second_zagreb_index_is_the_sum_over_bonds_of_the_degree_product():
    """M2 = sum over bonds of deg(u) * deg(v), same paper.

    By hand: n-butane's three bonds give 1*2 + 2*2 + 2*1 = 8; isobutane's three give 3*1 each = 9;
    neopentane 4*1 four times = 16; benzene 2*2 six times = 24.  Toluene's seven bonds:
    methyl-ipso 3, two ipso-ortho at 6, two ortho-meta at 4, two meta-para at 4 -- 3 + 12 + 8 + 8 = 31.
    """
    assert read_smiles('CCCC').zagreb_index(2) == 8
    assert read_smiles('CC(C)C').zagreb_index(2) == 9
    assert read_smiles('CC(C)(C)C').zagreb_index(2) == 16
    assert read_smiles('c1ccccc1').zagreb_index(2) == 24
    assert read_smiles('Cc1ccccc1').zagreb_index(2) == 31
    assert read_smiles('C').zagreb_index(2) == 0


def test_zagreb_index_defaults_to_the_first():
    m = read_smiles('Cc1ccccc1')
    assert m.zagreb_index() == m.zagreb_index(1) == 30


@mark.parametrize('order', [0, 3, 4, 100])
def test_zagreb_index_refuses_an_order_it_has_no_definition_for(order):
    """The paper defines two, so the method offers two and refuses the rest by name rather than
    returning something for an order nobody defined."""
    with raises(ValueError, match='order must be 1 or 2'):
        read_smiles('CCCC').zagreb_index(order)


def test_randic_index_reproduces_the_1975_branching_index():
    """chi = sum over bonds of 1 / sqrt(deg(u) * deg(v)).  Randic, JACS 97 (1975) 6609, where it is
    the branching index and the alkanes are tabulated to three decimals.

    n-butane: 1/sqrt(2) + 1/2 + 1/sqrt(2) = 1.914, which is his value.  n-hexane: two end bonds at
    1/sqrt(2) and three interior at 1/2 = 2.914, his value again.  2-methylpentane: 2/sqrt(3) for the
    two bonds off the branch carbon, 1/sqrt(6), 1/2 and 1/sqrt(2) = 2.770, his value.  Isobutane is
    3/sqrt(3) = sqrt(3) = 1.732 and neopentane 4/sqrt(4) = 2.
    """
    assert read_smiles('CCCC').randic_index == approx(2 / sqrt(2) + 0.5)
    assert round(read_smiles('CCCC').randic_index, 3) == 1.914

    assert read_smiles('CCCCCC').randic_index == approx(2 / sqrt(2) + 1.5)
    assert round(read_smiles('CCCCCC').randic_index, 3) == 2.914

    assert read_smiles('CC(C)CCC').randic_index == approx(
        2 / sqrt(3) + 1 / sqrt(6) + 0.5 + 1 / sqrt(2))
    assert round(read_smiles('CC(C)CCC').randic_index, 3) == 2.770

    assert read_smiles('CC(C)C').randic_index == approx(sqrt(3))
    assert read_smiles('CC(C)(C)C').randic_index == approx(2.0)
    assert read_smiles('c1ccccc1').randic_index == approx(3.0)   # six bonds at 1/2


def test_randic_index_of_a_molecule_with_no_bonds_is_zero():
    """A sum over bonds, so no bonds is 0.0 and a degree-0 atom never reaches the reciprocal square
    root -- there is no bond for it to be an endpoint of."""
    assert read_smiles('C').randic_index == 0.0
    assert read_smiles('[Na+].[Cl-]').randic_index == 0.0
    assert MoleculeContainer().randic_index == 0.0


def test_the_degree_indices_are_additive_over_components():
    """Both are sums over atoms or over bonds, so a salt is the sum of its parts -- no -1 to skip and
    nothing to refuse."""
    assert read_smiles('CCCC.CCCC').zagreb_index() == 20
    assert read_smiles('CCCC.CCCC').zagreb_index(2) == 16
    assert read_smiles('CCCC.CCCC').randic_index == approx(2 * (2 / sqrt(2) + 0.5))


# --- Balaban J ------------------------------------------------------------------------------------


@needs_numpy
def test_balaban_j_of_the_short_alkanes():
    """J = q / (mu + 1) * sum over bonds of 1 / sqrt(s(u) * s(v)), where q is the bond count, mu the
    cyclomatic number q - n + 1, and s(v) the sum of a vertex's distances to every other atom.
    Balaban, Chem. Phys. Lett. 89 (1982) 399.

    Derived by hand.  Ethane: s = 1, 1; q = 1, mu = 0; J = 1 * 1 = 1.  Propane: s = 3, 2, 3; two bonds
    at 1/sqrt(6); J = 2 * 2/sqrt(6) = 1.633.  n-Butane: s = 6, 4, 4, 6; bonds (6,4), (4,4), (4,6);
    J = 3 * (2/sqrt(24) + 1/4) = 1.975.  n-Pentane: s = 10, 7, 6, 7, 10; bonds (10,7), (7,6), (6,7),
    (7,10); J = 4 * (2/sqrt(70) + 2/sqrt(42)) = 2.191.  Isobutane: s = 3 for the centre and 5 for each
    leaf; three bonds; J = 3 * 3/sqrt(15) = 2.324.

    That ascending series -- 1.000, 1.633, 1.975, 2.191, with the branched isomer above its linear one
    at 2.324 -- is the discrimination the paper was written to demonstrate.
    """
    assert read_smiles('CC').balaban_j == approx(1.0)
    assert read_smiles('CCC').balaban_j == approx(2 * 2 / sqrt(6))
    assert round(read_smiles('CCC').balaban_j, 3) == 1.633
    assert read_smiles('CCCC').balaban_j == approx(3 * (2 / sqrt(24) + 0.25))
    assert round(read_smiles('CCCC').balaban_j, 3) == 1.975
    assert read_smiles('CCCCC').balaban_j == approx(4 * (2 / sqrt(70) + 2 / sqrt(42)))
    assert round(read_smiles('CCCCC').balaban_j, 3) == 2.191
    assert read_smiles('CC(C)C').balaban_j == approx(3 * 3 / sqrt(15))
    assert round(read_smiles('CC(C)C').balaban_j, 3) == 2.324


@needs_numpy
def test_balaban_j_of_a_ring_uses_the_cyclomatic_number():
    """Benzene: every vertex has s = 1 + 1 + 2 + 2 + 3 = 9, q = 6 and mu = 6 - 6 + 1 = 1, so the
    prefactor is 6/2 = 3 and J = 3 * 6 * (1/9) = 2.

    The mu in the denominator is why a ring does not simply out-score a chain of the same size: it is
    what makes J comparable across cyclic and acyclic molecules at all.  Cyclohexane is the same graph
    and the same 2.0 -- J reads distances, not bond orders.
    """
    assert read_smiles('c1ccccc1').balaban_j == approx(2.0)
    assert read_smiles('C1CCCCC1').balaban_j == approx(2.0)


@needs_numpy
def test_balaban_j_of_a_bicyclic_molecule_exercises_mu_above_one():
    """Naphthalene, and it is here because benzene cannot catch a denominator that goes wrong only above
    mu = 1.  Benzene's mu is 1, so `mu + 1` is 2 -- the same value a bare `2` would give, and the same
    value several plausible misreadings give.  Naphthalene has q = 11, n = 10, mu = 2 and a denominator
    of 3, which no off-by-one and no confusion of mu with the ring count reproduces.

    Distance sums by orbit, and they are checked rather than asserted: alpha (positions 1,4,5,8) sum to
    21, beta (2,3,6,7) to 25, the two fusion carbons to 17.  Those must reconcile with the Wiener index
    this file already pins, and they do -- 4*21 + 4*25 + 2*17 = 218 = 2 * 109, and 109 is naphthalene's
    Wiener value.  So the four bond classes are alpha-beta at 21*25 = 525 (four bonds), beta-beta at 625
    (two), alpha-fusion at 21*17 = 357 (four) and the fusion-fusion bond at 289.

    This phase does not have Balaban's own table entry in hand, so the hand derivation reconciling
    against the pinned Wiener index above is the whole assertion -- which it is, and it is stronger
    than a transcribed decimal.
    """
    naphthalene = read_smiles('c1ccc2ccccc2c1')
    assert naphthalene.balaban_j == approx(11 / 3 * (4 / sqrt(525) + 2 / 25 + 4 / sqrt(357) + 1 / 17))
    assert naphthalene.balaban_j == approx(1.9253677344386608, rel=PIN)
    # decalin is the same graph, so J is the same number -- the sibling of the benzene/cyclohexane pair
    # above, at the mu = 2 that pair cannot reach
    assert read_smiles('C1CCC2CCCCC2C1').balaban_j == approx(naphthalene.balaban_j)


@needs_numpy
def test_balaban_j_refuses_a_disconnected_molecule_and_names_split():
    """THE ONE REFUSAL IN F2, and it is at the answer boundary.

    A vertex distance sum over a disconnected graph is infinite -- the matrix says -1, which is not a
    distance and cannot be summed.  Every alternative is invented arithmetic: skipping the -1 makes s a
    within-component sum while q and mu stay global, which is a formula Balaban did not define and
    nobody has published.  So the answer is a refusal that names the repair, and the caller who wants
    a J per component runs `split()` and asks each part.
    """
    with raises(ValueError, match='split'):
        read_smiles('[Na+].[Cl-]').balaban_j
    with raises(ValueError, match='2 components'):
        read_smiles('CCCC.CCCC').balaban_j
    with raises(ValueError, match='split'):
        read_smiles('c1ccccc1.[Na+]').balaban_j


@needs_numpy
def test_balaban_j_of_each_half_of_a_salt_is_answerable():
    """The refusal is not a dead end: `split()` yields connected molecules and each one answers."""
    parts = read_smiles('CCCC.c1ccccc1').split()
    assert len(parts) == 2
    values = sorted(round(p.balaban_j, 3) for p in parts)
    assert values == [1.975, 2.0]


def test_balaban_j_of_a_single_atom_and_an_empty_molecule_is_zero():
    """A one-atom molecule is connected and has no bonds, so the sum is empty and J is 0.0 -- not a
    refusal, because nothing about it is disconnected.  An empty molecule has no components at all and
    answers 0.0 for the same reason."""
    assert read_smiles('C').balaban_j == 0.0
    assert read_smiles('O').balaban_j == 0.0
    assert MoleculeContainer().balaban_j == 0.0


@needs_numpy
@mark.parametrize('smiles', [
    'CC',                                             # ethane -- simplest non-trivial connected graph
    'CCCC',                                           # n-butane
    'CC(C)C',                                         # isobutane
    'c1ccccc1',                                       # benzene
    'C1CCCCC1',                                       # cyclohexane
    'c1ccc2[nH]ccc2c1',                               # indole
    'CC(=O)Oc1ccccc1C(=O)O',                          # aspirin
    'Cn1cnc2c1c(=O)n(C)c(=O)n2C',                     # caffeine
])
def test_balaban_j_is_invariant_under_renumbering(smiles):
    """balaban_j raises on disconnected molecules, so it is not in DESCRIPTORS (which is parametrized
    over a two-component salt).  Invariance is tested here, over connected molecules only.
    """
    molecule = read_smiles(smiles)
    other = _rebuilt_backwards(molecule)
    assert other.balaban_j == approx(molecule.balaban_j), smiles


# --- Bertz CT -------------------------------------------------------------------------------------


def test_bertz_ct_of_the_small_alkanes():
    """CT = [2N*log2(N) - sum over connection classes of n*log2(n)] + [n*log2(n) - sum over elements
    of m*log2(m)].  Bertz, JACS 103 (1981) 3599.

    A "connection" is a pair of bonds sharing an atom -- a path of three atoms -- and N is how many
    the molecule has, sum over atoms of C(deg, 2).  The classes are chython's stated reading: two
    connections are equivalent when their central atoms share a symmetry class AND their two outer
    atoms' classes match as an unordered pair (`atoms_order` supplies the classes).  The second
    bracket is the element diversity term.

    Derived by hand.  Ethane and propane: 0 and 1 connection, one class, all carbon -- both terms
    vanish and CT is 0.0, which is what an index of SYMMETRY-WEIGHTED SIZE says about a molecule with
    no diversity and nothing to distinguish.  n-Butane: two connections, both (middle | end, middle),
    one class of 2, so 2*2*1 - 2*1 = 2.  Isobutane: the centre has degree 3, so C(3,2) = 3 connections
    in one class: 2*3*log2(3) - 3*log2(3) = 3*log2(3) = 4.755, above n-butane -- branching is
    complexity, which is the paper's point.
    """
    assert read_smiles('C').bertz_ct == 0.0
    assert read_smiles('CC').bertz_ct == 0.0
    assert read_smiles('CCC').bertz_ct == 0.0
    assert read_smiles('CCCC').bertz_ct == approx(2.0)
    assert read_smiles('CC(C)C').bertz_ct == approx(3 * log2(3))
    assert read_smiles('CC(C)C').bertz_ct == approx(4.754887502163468, rel=PIN)
    assert MoleculeContainer().bertz_ct == 0.0


def test_bertz_ct_of_the_arenes():
    """Benzene: six connections, one class, so 2*6*log2(6) - 6*log2(6) = 6*log2(6) = 15.510.

    Toluene: eight connections in five classes of sizes 2, 1, 2, 2, 1 -- two ipso connections pairing
    the methyl with an ortho carbon, one pairing the two orthos, two at the orthos, two at the metas,
    one at the para.  2*8*3 - (2 + 0 + 2 + 2 + 0) = 48 - 6 = 42, and the element term is 0 because
    every atom is carbon.

    Naphthalene: fourteen connections in four classes of 4, 4, 2, 4 -- one at each alpha, one at each
    beta, one pairing the two alphas at each fusion carbon, and two pairing an alpha with the other
    fusion carbon.  2*14*log2(14) - (8 + 8 + 2 + 8) = 106.606 - 26 = 80.606.
    """
    assert read_smiles('c1ccccc1').bertz_ct == approx(6 * log2(6))
    assert read_smiles('c1ccccc1').bertz_ct == approx(15.509775004326936, rel=PIN)
    assert read_smiles('Cc1ccccc1').bertz_ct == approx(42.0)
    assert read_smiles('c1ccc2ccccc2c1').bertz_ct == approx(28 * log2(14) - 26)
    assert read_smiles('c1ccc2ccccc2c1').bertz_ct == approx(80.60593781761291, rel=PIN)


def test_bertz_ct_adds_the_element_diversity_term():
    """Phenol has toluene's skeleton, so the connection term is the same 42; its element term is
    7*log2(7) - 6*log2(6) - 1*log2(1) = 4.142, and CT is 46.142.

    Acetic acid: three connections at the carbonyl carbon, each its own class because the methyl
    carbon, the carbonyl oxygen and the hydroxyl oxygen are three distinct symmetry classes -- so
    2*3*log2(3) - 0 = 9.510 -- plus an element term of 4*2 - 2*1 - 2*1 = 4.  13.510 in all.
    """
    assert read_smiles('Oc1ccccc1').bertz_ct == approx(42 + 7 * log2(7) - 6 * log2(6))
    assert read_smiles('Oc1ccccc1').bertz_ct == approx(46.14170945007629, rel=PIN)
    assert read_smiles('CC(=O)O').bertz_ct == approx(6 * log2(3) + 4)
    assert read_smiles('CC(=O)O').bertz_ct == approx(13.509775004326936, rel=PIN)


def test_bertz_ct_reads_symmetry_and_not_bond_orders():
    """TWO PROPERTIES OF THE READING, both deliberate and both stated here so neither looks like a bug.

    Benzene and cyclohexane are the same graph with the same symmetry, so they get the same CT: bond
    orders reach this index only through `atoms_order`, and a vertex-transitive six-ring is
    vertex-transitive either way.  Neopentane also lands there -- C(4,2) = 6 connections in one class
    is arithmetically the same molecule as far as an information-content index is concerned.

    Phenol and chlorobenzene get the same CT too: the element term counts a partition, not which
    elements are in it.
    """
    assert read_smiles('C1CCCCC1').bertz_ct == approx(read_smiles('c1ccccc1').bertz_ct)
    assert read_smiles('CC(C)(C)C').bertz_ct == approx(6 * log2(6))
    assert read_smiles('Clc1ccccc1').bertz_ct == approx(read_smiles('Oc1ccccc1').bertz_ct)


def test_bertz_ct_of_a_salt_is_global_and_not_additive():
    """CT is defined for a disconnected molecule and needs no refusal -- there is no distance in it.

    It is NOT additive, and that is a property of every information-content index rather than a defect:
    two butanes have four connections in ONE class of 4, so 2*4*2 - 4*2 = 8, where one butane is 2.
    The duplicate raises N and enlarges the class at the same time, and the two do not cancel.
    """
    assert read_smiles('CCCC').bertz_ct == approx(2.0)
    assert read_smiles('CCCC.CCCC').bertz_ct == approx(8.0)
    assert read_smiles('[Na+].[Cl-]').bertz_ct == approx(2 * log2(2) - 2 * 0.0)
    assert read_smiles('[Na+].[Cl-]').bertz_ct == approx(2.0)


def test_bertz_ct_total_connections_matches_degree_formula():
    """Internal consistency: the N in the formula equals sum of C(deg, 2) over atoms.

    The two molecules exercise the two halves of the counting, which is why there are two.  Neopentane
    CC(C)(C)C puts all six connections at ONE centre -- a degree-4 atom, C(4,2) = 6 -- so it checks the
    pair count within one atom's neighbour list, on a rank whose population is 1.  Benzene spreads six
    connections over SIX centres of one rank, one each, so its N is recovered only if the population
    multiplier is applied: drop the multiplier and benzene's N falls to 1 while neopentane's stays 6.

    Both come out at CT = 6*log2(6), which is a coincidence of them having six connections in one class
    apiece, and the element term is zero for both -- all carbon.

    degree_of() is used here as the independent count and it saturates at 255, so this identity is
    stated for molecules well below that; the index itself reads the uncapped CSR row length.
    """
    for smiles, expect in (('CC(C)(C)C', 6), ('c1ccccc1', 6)):
        m = read_smiles(smiles)
        n_connections = sum(
            m.degree_of(a.n) * (m.degree_of(a.n) - 1) // 2
            for a in m.atoms()
        )
        assert n_connections == expect, smiles                      # derived by hand
        assert m.bertz_ct == approx(n_connections * log2(n_connections)), smiles


def test_bertz_ct_counts_a_dative_bond_as_a_connection():
    """THE GRAPH IS THE GRAPH AS STORED, ORDER 8 INCLUDED -- this file's global convention, tested here
    because a coordination contact changes N and there is no other index in F2 where it is this visible.

    Trimethylamine N(C)(C)C: the nitrogen's degree is 3, so N = C(3,2) = 3, and the three pairs are one
    class of 3 because the methyls are one rank.  2*3*log2(3) - 3*log2(3) = 3*log2(3) = 4.755, and the
    element term is 4*log2(4) - 3*log2(3) = 8 - 4.755 = 3.245.  CT is exactly 8.

    Give the nitrogen an iron to donate to and its degree becomes 4: N = 6 in two classes of 3, the
    methyl-methyl pairs and the iron-methyl pairs.  2*6*log2(6) - 2*3*log2(3) = 21.510, plus an element
    term of 5*log2(5) - 3*log2(3) = 6.855, so CT = 28.365.  The rise is the dative bond being counted.
    """
    assert read_smiles('N(C)(C)C').bertz_ct == approx(8.0)
    assert read_smiles('[Fe]~N(C)(C)C').bertz_ct == approx(12 * log2(6) - 6 * log2(3)
                                                           + 5 * log2(5) - 3 * log2(3))
    assert read_smiles('[Fe]~N(C)(C)C').bertz_ct == approx(28.364527976600278, rel=PIN)


# --- connectivity indices -------------------------------------------------------------------------


def test_chi_zero_is_the_sum_of_reciprocal_root_degrees():
    """0-chi = sum over atoms of 1/sqrt(delta), delta being the heavy-atom degree.  Kier and Hall,
    Rev. Comput. Chem. 2 (1991) 367-422, and the definition dates to their 1976 monograph.

    n-Butane's deltas are 1, 2, 2, 1: 1 + 2/sqrt(2) + 1 = 2 + sqrt(2) = 3.414, their tabulated value.
    Isobutane: three leaves and a degree-3 centre, 3 + 1/sqrt(3) = 3.577.  Benzene: six degree-2
    atoms, 6/sqrt(2) = 4.243.
    """
    assert read_smiles('CCCC').chi(0) == approx(2 + sqrt(2))
    assert round(read_smiles('CCCC').chi(0), 3) == 3.414
    assert read_smiles('CC(C)C').chi(0) == approx(3 + 1 / sqrt(3))
    assert read_smiles('c1ccccc1').chi(0) == approx(6 / sqrt(2))


def test_chi_one_is_the_randic_index():
    """1-chi is the sum over bonds of 1/sqrt(delta(u) * delta(v)) -- Randic's branching index by
    another name.  Two names, two code paths, one number: if these ever disagree, one is wrong."""
    for smiles in ('CCCC', 'CC(C)C', 'CC(C)(C)C', 'c1ccccc1', 'Cc1ccccc1', 'c1ccc2ccccc2c1',
                   'CC(=O)Oc1ccccc1C(=O)O'):
        m = read_smiles(smiles)
        assert m.chi(1) == approx(m.randic_index), smiles


def test_chi_two_and_three_walk_paths_of_three_and_four_atoms():
    """2-chi sums 1/sqrt(delta(u)*delta(v)*delta(w)) over three-atom paths, 3-chi over four-atom paths.
    Each path is counted ONCE.

    n-Butane has two three-atom paths, both with delta product 1*2*2 = 4, so 2-chi = 1/2 + 1/2 = 1.000
    -- Kier and Hall's tabulated value -- and one four-atom path with product 4, so 3-chi = 0.5.

    Isobutane has three three-atom paths (leaf-centre-leaf), each with product 1*3*1 = 3, so
    2-chi = 3/sqrt(3) = sqrt(3) = 1.732, and NO four-atom path at all: 3-chi = 0.0.

    Benzene has six three-atom paths (one centred at each atom, product 8) and six four-atom paths
    (one starting at each atom, product 16): 6/sqrt(8) = 2.121 and 6/4 = 1.5.
    """
    butane = read_smiles('CCCC')
    assert butane.chi(2) == approx(1.0)
    assert butane.chi(3) == approx(0.5)
    assert butane.chi(4) == approx(0.0)          # no five-atom path in four atoms

    isobutane = read_smiles('CC(C)C')
    assert isobutane.chi(2) == approx(sqrt(3))
    assert isobutane.chi(3) == approx(0.0)

    benzene = read_smiles('c1ccccc1')
    assert benzene.chi(2) == approx(6 / sqrt(8))
    assert benzene.chi(3) == approx(1.5)


def test_chi_of_a_hydrocarbon_is_the_same_valence_or_not():
    """delta-v = Zv - h, so a CH3 carbon is 4 - 3 = 1 and a CH2 is 4 - 2 = 2 -- exactly the heavy-atom
    degrees in a saturated hydrocarbon.  The two variants therefore agree on the alkanes, which is the
    check that the valence delta is built right before any heteroatom is involved.
    """
    for smiles in ('CCCC', 'CC(C)C', 'CC(C)(C)C', 'CCCCCC'):
        m = read_smiles(smiles)
        for order in (0, 1, 2, 3):
            assert m.chi(order, valence=True) == approx(m.chi(order)), (smiles, order)


def test_chi_valence_separates_a_heteroatom():
    """Ethanol: the deltas are 1, 2, 1 but the valence deltas are 1, 2 and 6 - 1 = 5 for the hydroxyl
    oxygen.

    0-chi = 1 + 1/sqrt(2) + 1 = 2.707 against 0-chi-v = 1 + 1/sqrt(2) + 1/sqrt(5) = 2.154.
    1-chi = 2/sqrt(2) = 1.414 against 1-chi-v = 1/sqrt(2) + 1/sqrt(10) = 1.023.

    THE FORMAL CHARGE IS NOT IN delta-v.  Kier and Hall define it as valence electrons less hydrogens,
    and that is what this is -- deliberately not `valence_electrons_count`'s per-atom term, which does
    subtract the charge because it is counting electrons rather than free connections.
    """
    ethanol = read_smiles('CCO')
    assert ethanol.chi(0) == approx(2 + 1 / sqrt(2))
    assert ethanol.chi(0, valence=True) == approx(1 + 1 / sqrt(2) + 1 / sqrt(5))
    assert ethanol.chi(1) == approx(2 / sqrt(2))
    assert ethanol.chi(1, valence=True) == approx(1 / sqrt(2) + 1 / sqrt(10))


def test_a_delta_of_zero_contributes_nothing():
    """CHYTHON'S STATED READING of the one case the formula cannot answer.

    1/sqrt(0) is not a number, so an atom whose delta is 0 contributes no term and no path through it
    contributes one.  It reaches the plain delta as an atom with no heavy neighbour -- water's oxygen,
    a lone counterion, methane's carbon -- and the valence delta as a fully hydrogenated atom, methane
    again (4 - 4 = 0).  The alternative is an infinite index, which no caller can use, and the
    alternative to stating the rule is two different silent answers for the same 1/sqrt(0).

    Water shows the two variants parting company: delta is 0 and delta-v is 6 - 2 = 4.
    """
    assert read_smiles('O').chi(0) == 0.0
    assert read_smiles('O').chi(0, valence=True) == approx(0.5)
    assert read_smiles('C').chi(0) == 0.0
    assert read_smiles('C').chi(0, valence=True) == 0.0
    assert read_smiles('[Na+].[Cl-]').chi(0) == 0.0
    assert read_smiles('CCCC.[Na+]').chi(0) == approx(2 + sqrt(2))   # the ion adds nothing
    assert MoleculeContainer().chi(0) == 0.0
    assert MoleculeContainer().chi(2) == 0.0


def test_the_valence_delta_does_not_double_count_an_explicit_hydrogen():
    """delta-v subtracts the IMPLICIT hydrogen count only, for the same reason
    `desc_valence_electrons` leaves `explicit_h` out of its sum: an explicit hydrogen is a vertex in
    this graph and carries its own delta-v of 1, so adding it to its neighbour's `h` as well subtracts
    it twice.  Under the double-counting spelling this molecule's carbon had delta-v 4 - 4 = 0, which
    zeroed every path through it and made a four-bonded carbon report no index at all.

    Kier and Hall define delta-v on the hydrogen-suppressed graph, where the explicit count is 0, so no
    published value in this file is computed from a changed expression -- only the explicit-H spelling
    moves, and it moves from a degenerate answer to a defined one.
    """
    explicit = read_smiles('[H]C([H])([H])[H]')
    # the carbon is Zv 4 with no implicit hydrogen left, so delta-v 4; each hydrogen is Zv 1, delta-v 1
    assert explicit.chi(1, valence=True) == approx(4 * (1 / sqrt(4.0)))
    assert explicit.chi(1, valence=True) == approx(explicit.chi(1))   # saturated C, so the two coincide
    assert explicit.chi(0, valence=True) == approx(1 / sqrt(4.0) + 4.0)
    # the suppressed spelling is unchanged and still the classical value: Zv - h = 4 - 4 = 0
    assert read_smiles('C').chi(0, valence=True) == 0.0


def test_chi_is_additive_over_components_when_every_delta_is_defined():
    assert read_smiles('CCCC.CCCC').chi(0) == approx(2 * (2 + sqrt(2)))
    assert read_smiles('CCCC.CCCC').chi(2) == approx(2.0)


@mark.parametrize('order', [5, 6, 20])
def test_chi_refuses_an_order_past_four(order):
    """Kier and Hall tabulate 0 through 4; past that the path enumeration grows exponentially and no
    published index uses it.  A refusal that says so beats an answer nobody can check."""
    with raises(ValueError, match='order must be 0-4'):
        read_smiles('CCCCCCCC').chi(order)


def test_chi_valence_refuses_what_the_valence_delta_cannot_state():
    """Same two refusals as `valence_electrons_count`, for the same two reasons -- and NOT for the
    plain variant, which needs neither the element's Zv nor a hydrogen count."""
    with raises(ValueError, match='no valence electron count'):
        read_smiles('[Ce]').chi(0, valence=True)
    assert read_smiles('[Ce]').chi(0) == 0.0

    m = MoleculeContainer()
    with m.edit():
        m.add_atom('C', implicit_h=H_UNKNOWN)
        m.add_atom('C')
        m.add_bond(1, 2, 1)
    with raises(ValueError, match='calc_implicit'):
        m.chi(1, valence=True)
    assert m.chi(1) == approx(1.0)


# --- electrotopological state ---------------------------------------------------------------------
#
# NUMPY FOR TWO REASONS RATHER THAN ONE, which is why every test in this section is marked and not just
# the ones that reach a distance: the perturbation sum reads `distance_matrix`, AND both answers are
# `(n,)` float64 arrays, so even the refusal path allocates one.  A marker per test in the section would
# be nine copies of the same reason.


@needs_numpy
def test_estate_ethane_is_two_by_symmetry():
    """DERIVED: each C has d=1, dv=4-3=1, N=2, so I = ((2/2)**2*1 + 1)/1 = 2.  Equal I, so the
    perturbation sum is 0 and S = I.

    `rel=PIN` AND NOT THE DEFAULT: a closed form compared against the same closed form needs only the
    slack of a differing arithmetic path, and 1e-6 would swallow a wrong period or a wrong delta.
    """
    assert read_smiles('CC').estate_indices() == approx([2.0, 2.0], rel=PIN)


@needs_numpy
def test_estate_propane_hand_derived():
    """DERIVED: terminal C -> I = ((1)*1 + 1)/1 = 2; middle C -> I = ((1)*2 + 1)/2 = 1.5.
    S(term) = 2 + (2-1.5)/2**2 + (2-2)/3**2 = 2.125
    S(mid)  = 1.5 + 2*(1.5-2)/2**2 = 1.25

    `rel=PIN`: both values are exact in binary, so the only difference a correct implementation can show
    is the last bit of a different summation order.
    """
    assert read_smiles('CCC').estate_indices() == approx([2.125, 1.25, 2.125], rel=PIN)


@needs_numpy
def test_estate_ethanol_hand_derived():
    """DERIVED: C1 I=2, C2 I=(2+1)/2=1.5, O I=((1)*(6-1)+1)/1=6; distances 1, 2, 1.
    S(C1) = 2 + (2-1.5)/4 + (2-6)/9   = 121/72
    S(C2) = 1.5 + (1.5-2)/4 + (1.5-6)/4 = 0.25
    S(O)  = 6 + (6-2)/9 + (6-1.5)/4   = 545/72

    `rel=PIN`, and this is the one where the tolerance earns its keep: 121/72 and 545/72 are NOT exact
    in binary, and `121 / 72` here versus `2 + 0.5/4 - 4/9` in the code are two paths to the same
    rational.  PIN passes their few-ULP difference and refuses anything larger.
    """
    assert read_smiles('CCO').estate_indices() == approx([121 / 72, 0.25, 545 / 72], rel=PIN)


@needs_numpy
def test_estate_sum_equals_intrinsic_sum():
    """STRUCTURAL: every perturbation term appears twice with opposite signs, so the sum of S over the
    molecule equals the sum of I.  Holds for any connected molecule and is the invariant a refactor
    would break first.

    NO `rel=PIN` HERE, DELIBERATELY, and the two reasons are why the rule in this file is per-assertion
    rather than per-file.  It compares two COMPUTED quantities, so what it pins is a relation and not a
    literal -- there is no transcription to protect.  And it sums a signed series whose terms cancel, so
    the cancellation earns it real slack: tightening this one to PIN would make it fail on a molecule
    large enough for the cancellation to lose digits, which is a false alarm about arithmetic and not a
    finding about EState.
    """
    for smiles in ('CCO', 'CCC', 'c1ccccc1O', 'CC(=O)Nc1ccccc1'):
        m = read_smiles(smiles)
        assert sum(m.estate_indices()) == approx(sum(m.estate_intrinsic_states())), smiles


@needs_numpy
def test_estate_zero_degree_atom_is_nan_not_zero():
    """0.0 is a legal EState value, so a lone atom cannot be reported as 0.0."""
    from math import isnan
    assert all(isnan(v) for v in read_smiles('O').estate_indices())
    values = read_smiles('CCO.[Na+]').estate_indices()
    assert isnan(values[3]) and not any(isnan(v) for v in values[:3])


@needs_numpy
def test_estate_ignores_a_pair_in_another_component():
    """A distance across components is -1 and nothing sums a -1, so a salt's organic part answers
    exactly what it answers alone.

    NO `rel=PIN` HERE EITHER, and for the first of the two reasons above: both sides are computed by the
    same code over the same three atoms, so this is a relation between two computations and the default
    is the honest tolerance for it.  (It would in fact pass at `PIN` today -- the two sides are
    bit-identical -- which is exactly why asserting at PIN would be asserting something this test does
    not mean to claim.)
    """
    assert read_smiles('CCO.[Na+]').estate_indices()[:3] == approx(
        read_smiles('CCO').estate_indices())


@needs_numpy
def test_estate_refuses_what_the_valence_delta_cannot_state():
    """The same two refusals `chi(valence=True)` has, and for the same reason: one derivation of Zv - h.

    An atom whose Zv is not stated has no intrinsic state, so this lets `_desc_delta_valence`'s two
    messages through unchanged rather than spelling a third.
    """
    with raises(ValueError, match='no valence electron count'):
        read_smiles('[Gd]C').estate_indices()

    m = MoleculeContainer()
    with m.edit():
        m.add_atom('C', implicit_h=H_UNKNOWN)
        m.add_atom('C')
        m.add_bond(1, 2, 1)
    with raises(ValueError, match='calc_implicit'):
        m.estate_intrinsic_states()


@needs_numpy
def test_estate_is_the_shape_and_dtype_the_other_per_atom_answers_are():
    """`(n,)` float64 in `atom_numbers` order, so it drops into the same slot `atom_invariants` fills."""
    m = read_smiles('CC(=O)Nc1ccccc1')
    for values in (m.estate_indices(), m.estate_intrinsic_states()):
        assert values.dtype.name == 'float64'
        assert values.shape == (m.atom_count,)
    assert MoleculeContainer().estate_indices().shape == (0,)


# --- shape indices --------------------------------------------------------------------------------


def test_hall_kier_alpha_sums_the_published_atom_contributions():
    """alpha(atom) = r_cov(atom) / r_cov(Csp3) - 1, with r_cov(Csp3) = 0.77 A; Hall and Kier,
    Rev. Comput. Chem. 2 (1991) 367-422, whose table this reproduces to its two decimals:

        Csp3  0.00   Csp2 -0.13   Csp  -0.22
        Nsp3 -0.04   Nsp2 -0.20   Nsp  -0.29
        Osp3 -0.04   Osp2 -0.20
        F    -0.07   Cl    0.29   Br    0.48   I  0.73
        Psp3  0.43   Psp2  0.30
        Ssp3  0.35   Ssp2  0.22

    Aromatic counts as sp2, so benzene is 6 * -0.13 = -0.78 and cyclohexane is exactly 0.0 -- the
    index measures how far the atoms are from an sp3 carbon, and cyclohexane is all sp3 carbon.
    """
    assert read_smiles('c1ccccc1').hall_kier_alpha == approx(-0.78)
    assert read_smiles('C1CCCCC1').hall_kier_alpha == 0.0
    assert read_smiles('Cc1ccccc1').hall_kier_alpha == approx(-0.78)   # the methyl adds 0
    assert read_smiles('CCO').hall_kier_alpha == approx(-0.04)
    assert read_smiles('CC(=O)O').hall_kier_alpha == approx(-0.37)     # 0 - 0.13 - 0.20 - 0.04
    assert read_smiles('Clc1ccccc1').hall_kier_alpha == approx(-0.78 + 0.29)
    assert read_smiles('c1ccncc1').hall_kier_alpha == approx(5 * -0.13 - 0.20)
    assert read_smiles('CC#N').hall_kier_alpha == approx(-0.22 - 0.29)  # sp carbon, sp nitrogen
    assert MoleculeContainer().hall_kier_alpha == 0.0


def test_hall_kier_alpha_of_an_element_the_table_omits_is_zero():
    """CHYTHON'S STATED READING: the table is the paper's and nothing is extrapolated for an element it
    does not list, so a metal contributes 0.0.

    That is a REFERENCE, not a measurement -- 0.0 means "treated as an sp3 carbon" and it is what makes
    kappa answerable for an organometallic rather than a refusal.  A caller who needs a metal's radius
    correction has to supply it; chython does not invent one.
    """
    assert read_smiles('[Na+].[Cl-]').hall_kier_alpha == approx(0.29)   # the chloride only
    assert read_smiles('[Fe]').hall_kier_alpha == 0.0


def test_hall_kier_alpha_reads_sulfur_and_phosphorus_by_coordination_not_by_z():
    """ALPHA IS A COVALENT-RADIUS CORRECTION, so what picks the paper's sp3 row over its sp2 row is how
    many sigma bonds the atom holds -- not how many formal double bonds someone wrote on it.

    S and P are the only elements where the two answers differ.  Chython's own perception makes a sulfone
    S `z5` and a phosphate P `z2`, and reading either as sp2 would take 0.94 A and 1.00 A where the atom
    is four-coordinate and tetrahedral and its radius is 1.04 A and 1.10 A.  The rule is instead
    `hybridization != sp3 and degree <= 2`: the shortened sp2 radius belongs to a LOW-COORDINATE atom
    that is genuinely pi-bonded.

    Each expectation below is written as its sum rather than as a decimal, so the arithmetic is the
    assertion:

        thioether CSC        S is z1              0.35
        thiophene            S is z4, 2 bonds     0.22, and four aromatic carbons at -0.13
        thioacetone          S is z2, 1 bond      0.22, and the thiocarbonyl carbon at -0.13
        DMSO                 S is z2, 3 bonds     0.35 -- pyramidal, so sp3
        dimethyl sulfone     S is z5, 4 bonds     0.35 -- tetrahedral, so sp3
        methanesulfonamide   the same S           0.35, with an sp3 N
        trimethylphosphine   P is z1              0.43
        phosphoric acid      P is z2, 4 bonds     0.43 -- tetrahedral, so sp3

    NITROGEN IS NOT AN EXCEPTION and the last row is why the rule is not "z5 is always sp3": a nitro N is
    z5 and three-coordinate, but it is planar and its radius is the sp2 one.  Reading coordination rather
    than z gets both cases right with one test.
    """
    assert read_smiles('CSC').hall_kier_alpha == approx(0.35)
    assert read_smiles('c1ccsc1').hall_kier_alpha == approx(4 * -0.13 + 0.22)
    assert read_smiles('CC(=S)C').hall_kier_alpha == approx(-0.13 + 0.22)

    assert read_smiles('CS(=O)C').hall_kier_alpha == approx(0.35 - 0.20)
    assert read_smiles('CS(=O)(=O)C').hall_kier_alpha == approx(0.35 - 2 * 0.20)
    assert read_smiles('CS(=O)(=O)N').hall_kier_alpha == approx(0.35 - 2 * 0.20 - 0.04)

    assert read_smiles('CP(C)C').hall_kier_alpha == approx(0.43)
    assert read_smiles('OP(=O)(O)O').hall_kier_alpha == approx(0.43 - 0.20 - 3 * 0.04)

    assert read_smiles('CN(=O)=O').hall_kier_alpha == approx(-0.20 - 2 * 0.20)


def test_kappa_one_and_two_of_the_butanes():
    """Kier's shape indices, from the extremal path counts of a graph with n atoms:

        kappa1 = n(n-1)^2 / P1^2
        kappa2 = (n-1)(n-2)^2 / P2^2
        kappa3 = (n-1)(n-3)^2 / P3^2   for n odd
                 (n-3)(n-2)^2 / P3^2   for n even

    P_m is the number of paths of m bonds -- the same paths `chi` walks, unweighted.  Kier and Hall,
    Rev. Comput. Chem. 2 (1991) 367-422.

    n-Butane: n = 4, P1 = 3, P2 = 2, P3 = 1, so kappa1 = 4*9/9 = 4, kappa2 = 3*4/4 = 3, and
    kappa3 = 1*4/1 = 4 by the even branch.  Isobutane has the same 4 bonds... the same P1 = 3, hence
    the same kappa1 = 4 -- kappa1 cannot see branching, which is exactly why kappa2 exists: its
    P2 = C(3, 2) = 3 gives 3*4/9 = 1.333 against n-butane's 3.
    """
    butane = read_smiles('CCCC')
    assert butane.kappa(1) == approx(4.0)
    assert butane.kappa(2) == approx(3.0)
    assert butane.kappa(3) == approx(4.0)

    isobutane = read_smiles('CC(C)C')
    assert isobutane.kappa(1) == approx(4.0)
    assert isobutane.kappa(2) == approx(4 / 3)


def test_kappa_of_a_ring():
    """Benzene: n = 6 and P1 = P2 = P3 = 6, so kappa1 = 6*25/36 = 4.167,
    kappa2 = 5*16/36 = 2.222 and kappa3 = 3*16/36 = 1.333 by the even branch.  Cyclohexane is the same
    graph and the same three numbers -- the kappas read paths, not bond orders, and the alpha variant is
    what makes an aromatic ring differ from a saturated one.
    """
    benzene = read_smiles('c1ccccc1')
    assert benzene.kappa(1) == approx(6 * 25 / 36)
    assert benzene.kappa(2) == approx(5 * 16 / 36)
    assert benzene.kappa(3) == approx(3 * 16 / 36)
    for order in (1, 2, 3):
        assert read_smiles('C1CCCCC1').kappa(order) == approx(benzene.kappa(order)), order


def test_kappa_alpha_shifts_every_term_by_the_hall_kier_alpha():
    """The alpha variant replaces n by n + alpha and P by P + alpha:

        kappa1_alpha = (n+a)(n+a-1)^2 / (P1+a)^2

    Benzene's alpha is -0.78, so kappa1_alpha = 5.22 * 4.22^2 / 5.22^2 = 3.412 against the plain 4.167.
    Cyclohexane's alpha is exactly 0.0 (all sp3 carbons, the reference point), so its alpha variant
    equals its plain one -- this is the NO-SHIFT CLAIM and carries no value check; the value assertions
    that confirm both alpha places are touched are in
    test_kappa_alpha_of_a_molecule_whose_atom_count_and_path_count_differ.
    """
    benzene = read_smiles('c1ccccc1')
    a = -0.78
    assert benzene.kappa(1, alpha=True) == approx((6 + a) * (6 + a - 1) ** 2 / (6 + a) ** 2)
    assert benzene.kappa(1, alpha=True) == approx(3.4115708812260537, rel=PIN)

    cyclohexane = read_smiles('C1CCCCC1')
    for order in (1, 2, 3):
        assert cyclohexane.kappa(order, alpha=True) == approx(cyclohexane.kappa(order)), order


def test_kappa_alpha_of_a_molecule_whose_atom_count_and_path_count_differ():
    """Benzene cannot catch a swap of n for P: it has n = P1 = P2 = P3 = 6, so reading the wrong one
    gives the same answer.  Naphthalene separates them -- n = 10, P1 = 11 (bonds), P2 = 14 and P3 = 18.

    alpha is 10 * -0.13 = -1.30 (ten aromatic carbons), so n + alpha = 8.70 and P1 + alpha = 9.70:

        kappa1_alpha = 8.70 * 7.70**2 / 9.70**2 = 515.823 / 94.09 = 5.4822298

    Reading n where P belongs gives 6.8126 and the reverse gives 5.9367, so all three are distinct.

    kappa2 alpha uses P2 = 14: (9+a)(8+a)^2 / (14+a)^2.  Plain value 9*64/196 = 2.939 confirms P2=14.
    kappa3 alpha uses P3 = 18 (even branch): (7+a)(8+a)^2 / (18+a)^2.  Plain 7*64/324 = 1.383 confirms P3=18.
    A substitution that misses either alpha place -- numerator or denominator -- gives a different number.

    THE LAST LITERAL IS TIGHTER THAN `approx`, deliberately.  Its predecessor was 5.482234032309491 --
    my own arithmetic slip, wrong from the sixth decimal -- and it passed, because `approx`'s default
    relative tolerance is 1e-6 and the error was 7.8e-7.  A pinned literal exists to catch a drift the
    closed form above cannot, so one that agrees only to the tolerance pins nothing; 94.09 * 5.4822298 is
    515.82300 and that is the check.
    """
    naphthalene = read_smiles('c1ccc2ccccc2c1')
    a = 10 * -0.13
    assert naphthalene.hall_kier_alpha == approx(a)
    assert naphthalene.kappa(1, alpha=True) == approx((10 + a) * (10 + a - 1) ** 2 / (11 + a) ** 2)
    assert naphthalene.kappa(1, alpha=True) == approx(5.482229779997874, rel=PIN)
    assert naphthalene.kappa(2, alpha=True) == approx((9 + a) * (8 + a) ** 2 / (14 + a) ** 2)
    assert naphthalene.kappa(2, alpha=True) == approx(2.143052886105772, rel=PIN)
    assert naphthalene.kappa(3, alpha=True) == approx((7 + a) * (8 + a) ** 2 / (18 + a) ** 2)
    assert naphthalene.kappa(3, alpha=True) == approx(0.9174692531105451, rel=PIN)


def test_kappa_of_a_molecule_with_no_path_of_that_length_is_zero():
    """P3 of isobutane is 0, so kappa3's denominator is 0.  0.0 rather than an exception or a nan: the
    molecule has no three-bond path, so there is no three-bond shape to report, and a nan poisons the
    descriptor row it goes into.  Same for a single atom and an empty molecule.
    """
    assert read_smiles('CC(C)C').kappa(3) == 0.0
    assert read_smiles('C').kappa(1) == 0.0
    assert read_smiles('CC').kappa(2) == 0.0
    assert MoleculeContainer().kappa(1) == 0.0


@mark.parametrize('order', [0, 4, 5])
def test_kappa_refuses_an_order_other_than_one_two_or_three(order):
    """Kier defines three, each with its own extremal graph, and kappa0 is not a thing."""
    with raises(ValueError, match='order must be 1, 2 or 3'):
        read_smiles('CCCC').kappa(order)


def test_kappa_values_invert_back_to_the_path_counts_counted_by_hand():
    """Recovers P_m from the kappa value and checks it against paths counted by hand.

    The internal consistency that keeps ONE path enumerator honest: `chi(m)` over a graph whose every
    delta is 1 is the path count, and that is exactly how kappa gets P_m.  Recover P1, P2, P3 from the
    kappa values and check them against the paths counted by hand.
    """
    butane = read_smiles('CCCC')
    # kappa1 = n(n-1)^2 / P1^2  =>  P1 = sqrt(n(n-1)^2 / kappa1)
    assert sqrt(4 * 9 / butane.kappa(1)) == approx(3.0)     # three bonds
    assert sqrt(3 * 4 / butane.kappa(2)) == approx(2.0)     # two three-atom paths
    assert sqrt(1 * 4 / butane.kappa(3)) == approx(1.0)     # one four-atom path


def test_kappa3_takes_the_odd_atom_count_branch():
    """Toluene has n = 7 atoms (odd) and P3 = 8 three-bond paths, so kappa3 uses the odd branch:

        kappa3 = (n-1)(n-3)^2 / P3^2 = 6 * 4^2 / 8^2 = 96 / 64 = 1.5

    The even branch would give (n-3)(n-2)^2 / P3^2 = 4 * 5^2 / 64 = 100 / 64 = 1.5625, so the
    assertion separates the two branches: 1.5 is the odd-branch answer and 1.5625 is the even one.

    Butane (n=4) and benzene (n=6) are both even, so neither exercises this branch.
    """
    assert read_smiles('Cc1ccccc1').kappa(3) == approx(1.5)


# --- the disconnected policy ----------------------------------------------------------------------


# Marked whole rather than split: the value of this test is that ONE body answers "what does chython do
# about salts" for every descriptor at once, and four of its rows are distance-derived.  A version that
# kept the additive rows running without numpy would be a second, shorter table of the same policy,
# which is precisely the duplication the docstring below argues against.
@needs_numpy
def test_the_disconnected_policy_in_one_table():
    """WHAT A SALT ANSWERS, per descriptor, in one place.  Benzene and a sodium ion:

      additive over components   every count, the degree indices, hall_kier_alpha
      skips the missing pairs    wiener_index
      within the component       eccentricities, and so graph_radius (0, the ion's) and graph_diameter
      global, not additive       bertz_ct, kappa
      refuses                    balaban_j -- the only one

    Each of those is asserted in its own descriptor's test too; the value of having them in one place is
    that a reader asking "what does chython do about salts" gets one answer instead of nine, and a future
    descriptor that picks a different policy has to change a table that says out loud what the others do.

    `connected_components_count` is what tells a caller the molecule was disconnected in the first place,
    which is why it is asserted first.
    """
    m = read_smiles('c1ccccc1.[Na+]')
    assert m.connected_components_count == 2

    assert m.carbon_count == 6
    assert m.heteroatoms_count == 1
    assert m.aromatic_rings_count == 1
    assert m.fused_ring_systems_count == 1
    assert m.zagreb_index() == 24                     # the ion has degree 0 and adds nothing
    assert m.randic_index == approx(3.0)
    assert m.hall_kier_alpha == approx(-0.78)

    assert m.wiener_index == 27                       # the six-by-six block only

    assert list(m.eccentricities()) == [3, 3, 3, 3, 3, 3, 0]
    assert m.graph_radius == 0                        # the ion's, and that is the point
    assert m.graph_diameter == 3

    # global rather than additive: neither is the value of the benzene alone, which is why they are
    # asserted against their own numbers rather than against a component's
    assert m.bertz_ct == approx(19.651484454403228, rel=PIN)
    assert m.kappa(1) == approx(7.0)

    with raises(ValueError, match='split'):
        m.balaban_j


# Marked whole because `wiener_index` is the half that carries the claim: a count recomputed after an
# edit is cheap either way, where a memoised distance matrix is the tempting cache this forbids.  Keeping
# the `carbon_count` half alive without numpy would leave the test passing while no longer testing the
# case it was written for.
@needs_numpy
def test_no_descriptor_is_cached():
    """An edit changes the answer, because nothing in this file memoises.  A cached derived number is a
    second truth an edit can contradict, which is why `functional_groups()` is a method and why none
    of these is a `cached_property`.
    """
    m = read_smiles('c1ccccc1')
    assert m.carbon_count == 6
    assert m.wiener_index == 27
    with m.edit():
        new = m.add_atom('C')
        m.add_bond(1, new, 1)
    assert m.carbon_count == 7
    assert m.wiener_index == 42                       # toluene's, and the arene test pins that
