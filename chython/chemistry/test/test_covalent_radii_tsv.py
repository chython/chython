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
"""`tables/covalent_radii.tsv`, CHECKED AGAINST ITS PURPOSE and not against a transcription.

A radius here exists to answer one question -- is this pair of atoms bonded -- so the table is probed
with measured bond lengths and measured nonbonded contacts of public compounds.  A mistyped digit
that matters shows up as a bond length the table rejects or a contact it accepts; one that does not
matter is not worth a test.  The structural tests above them catch a shifted or duplicated row.
"""
from inspect import signature

from .._perceive import perceive_bonds
from .._tables import covalent_radii
from ...core._core import element_symbols


#: The multiplier the pass applies, read from the pass rather than repeated here.
MULTIPLIER = signature(perceive_bonds).parameters['radius_multiplier'].default

#: Measured bond lengths, in Angstroms.  Every one must be at or under the threshold, or the pass
#: misses the bond.  The tightest of them, F-F, is what sets the multiplier's floor.
BONDED = {
    'H2 H-H': ('H', 'H', .741),
    'N2 N#N': ('N', 'N', 1.098),
    'O2 O=O': ('O', 'O', 1.208),
    'F2 F-F': ('F', 'F', 1.412),
    'Cl2 Cl-Cl': ('Cl', 'Cl', 1.988),
    'Br2 Br-Br': ('Br', 'Br', 2.281),
    'I2 I-I': ('I', 'I', 2.666),
    'HF H-F': ('H', 'F', .917),
    'HCl H-Cl': ('H', 'Cl', 1.275),
    'water O-H': ('O', 'H', .958),
    'ammonia N-H': ('N', 'H', 1.012),
    'methane C-H': ('C', 'H', 1.087),
    'ethane C-C': ('C', 'C', 1.535),
    'acetylene C#C': ('C', 'C', 1.203),
    'benzene C-C': ('C', 'C', 1.397),
    'carbon dioxide C=O': ('C', 'O', 1.163),
    'hydrogen peroxide O-O': ('O', 'O', 1.475),
    'hydrazine N-N': ('N', 'N', 1.447),
    'oxygen difluoride O-F': ('O', 'F', 1.405),
    'carbon disulfide C=S': ('C', 'S', 1.553),
    'hydrogen sulfide S-H': ('S', 'H', 1.336),
    'phosphine P-H': ('P', 'H', 1.420),
    'silane Si-H': ('Si', 'H', 1.480),
    'sulfur hexafluoride S-F': ('S', 'F', 1.564),
    'tetrafluoromethane C-F': ('C', 'F', 1.319),
    'tetrachloromethane C-Cl': ('C', 'Cl', 1.767),
    'tetrabromomethane C-Br': ('C', 'Br', 1.942),
    'iodomethane C-I': ('C', 'I', 2.132),
    'white phosphorus P-P': ('P', 'P', 2.210),
    'cyclooctasulfur S-S': ('S', 'S', 2.050),
    'disilane Si-Si': ('Si', 'Si', 2.330),
    'diborane B-H': ('B', 'H', 1.190),
    'ferrocene Fe-C': ('Fe', 'C', 2.064),
    'sodium chloride Na-Cl': ('Na', 'Cl', 2.361),
    'mercury(II) chloride Hg-Cl': ('Hg', 'Cl', 2.250),
    'tetramethyltin Sn-C': ('Sn', 'C', 2.144),
    'aluminium chloride dimer Al-Cl': ('Al', 'Cl', 2.060),
    'tetrachloroplatinate Pt-Cl': ('Pt', 'Cl', 2.320),
}

#: Measured nonbonded contacts.  Every one must be over the threshold, or the pass invents a bond.
#: The tightest, cyclobutadiene's transannular carbons, is what sets the multiplier's ceiling.
CONTACTS = {
    # From the accepted D2h rectangle, 1.344 and 1.441 A sides, rather than from a measured distance.
    'cyclobutadiene C...C (transannular)': ('C', 'C', 1.970),
    'water H...H (1,3)': ('H', 'H', 1.514),
    'methane H...H (1,3)': ('H', 'H', 1.775),
    'benzene H...H (ortho)': ('H', 'H', 2.481),
    'ethane C...H (1,3)': ('C', 'H', 2.160),
    'benzene C...C (meta)': ('C', 'C', 2.420),
    'water dimer H...O': ('H', 'O', 1.992),
    'water dimer O...O': ('O', 'O', 2.950),
    'graphite interlayer C...C': ('C', 'C', 3.354),
    'solid argon Ar...Ar': ('Ar', 'Ar', 3.760),
}


#: The nonbonded contact NO multiplier separates from a bond, and the bond it collides with.  A pair
#: this tight leaves the two corpora with no window at all, so it is stated here rather than in
#: `CONTACTS`: what the test below asserts is the collision, not a threshold that survives it.
UNSEPARABLE = ('bicyclo[1.1.1]pentane C1...C3', ('C', 'C', 1.845), 'F2 F-F')


def _threshold(a: str, b: str) -> float:
    radii = covalent_radii()
    numbers = {s: z for z, s in enumerate(element_symbols())}
    return (radii[numbers[a]] + radii[numbers[b]]) * MULTIPLIER


def test_the_table_states_a_radius_for_every_element_it_covers():
    """1..96 inclusive: the range the crystallographic survey covers, and no invented row past it."""
    radii = covalent_radii()
    assert set(radii) == set(range(1, 97))


def test_the_r_marker_has_no_row():
    """Element 0 is a marker rather than an element, so it has no radius and gets no bond."""
    assert 0 not in covalent_radii()


def test_every_radius_is_a_plausible_length():
    for z, radius in covalent_radii().items():
        assert .2 < radius < 2.7, f'element {z} states {radius} A'


def test_the_symbol_column_agrees_with_the_element_table():
    """A shifted row is a wrong radius for every element after it, and nothing else would notice."""
    from .._tables import read_table

    symbols = element_symbols()
    for row in read_table('covalent_radii.tsv'):
        assert symbols[int(row['z'])] == row['symbol'], row


def test_helium_is_the_smallest_and_the_alkali_metals_grow_down_the_group():
    """Helium below hydrogen and francium above caesium: two orderings a shifted column breaks."""
    radii = covalent_radii()
    assert min(radii, key=radii.get) == 2
    assert radii[2] < radii[1]
    assert radii[87] > radii[55] > radii[37] > radii[19] > radii[11] > radii[3]


def test_every_measured_bond_length_is_under_the_threshold():
    missed = {name: (length, round(_threshold(a, b), 3))
              for name, (a, b, length) in BONDED.items() if length > _threshold(a, b)}
    assert not missed, f'bond(s) the radii would miss, as length vs threshold: {missed}'


def test_every_measured_contact_is_over_the_threshold():
    invented = {name: (length, round(_threshold(a, b), 3))
                for name, (a, b, length) in CONTACTS.items() if length <= _threshold(a, b)}
    assert not invented, f'contact(s) the radii would bond, as length vs threshold: {invented}'


def test_the_multiplier_sits_inside_the_window_the_corpus_leaves():
    """The two corpora bracket the multiplier, and this states by how much.

    The floor is the longest bond over its radius sum, the ceiling the shortest contact over its
    own.  A default outside them cannot answer both corpora, whatever the radii are.
    """
    floor = max(length / (_threshold(a, b) / MULTIPLIER) for a, b, length in BONDED.values())
    ceiling = min(length / (_threshold(a, b) / MULTIPLIER) for a, b, length in CONTACTS.values())
    assert floor < ceiling, f'no multiplier answers both corpora: floor {floor:.3f}, ceiling {ceiling:.3f}'
    assert floor <= MULTIPLIER <= ceiling, f'{MULTIPLIER} is outside {floor:.3f}..{ceiling:.3f}'


def test_one_threshold_cannot_separate_a_strained_bridgehead_from_a_bond():
    """The limit of a single distance rule, stated as a measurement rather than left to be discovered.

    Rejecting bicyclo[1.1.1]pentane's 1.845 A bridgehead contact needs a multiplier under the one that
    reaches fluorine's 1.412 A bond, so the pass bonds that pair and the strain is what it reads as
    connectivity.  A caller for whom that matters passes a tighter ``radius_multiplier`` and loses F2.
    """
    name, (a, b, length), collides_with = UNSEPARABLE
    needed = length / (_threshold(a, b) / MULTIPLIER)
    ca, cb, bond = BONDED[collides_with]
    floor = bond / (_threshold(ca, cb) / MULTIPLIER)
    assert needed < floor, \
        f'{name} at {needed:.4f} no longer collides with {collides_with} at {floor:.4f}; the window ' \
        'has moved and the pass can now answer both'
    assert length <= _threshold(a, b), f'{name} is not bonded at {MULTIPLIER}, so the docstring is stale'
