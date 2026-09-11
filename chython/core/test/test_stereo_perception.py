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
"""Which candidate units are genuinely stereogenic.

`test_stereo_units.py` pins what perception EMITS -- every site that could carry a configuration.
This file pins which of those the automorphism group leaves standing: a unit is stereogenic unless
some automorphism of the constitution fixes its anchor, permutes its directions oddly, and is
consistent with every other unit's stored configuration.
"""
import re
from collections import Counter

import pytest

from chython.core import MoleculeContainer, _core

_SYMBOL = re.compile(r'[A-Z][a-z]?')


def _symbols(atoms):
    """'CCClCClCClC' -> 8 symbols. Never iterate the string: 'Cl' is two characters."""
    return _SYMBOL.findall(atoms)


def _mol(*, atoms, bonds, hydrogens=None, charges=None, parities=None):
    """Build a molecule from an element string, a bond list and an explicit hydrogen count.

    `hydrogens` is not optional in practice. The core never DERIVES an implicit hydrogen count
    (test_derive.py), so a carbon written with three heavy neighbours and no count has three
    directions, not four, and is refused before stereogenicity is ever asked about. Every record
    below that needs a hydrogen direction states it.
    """
    m = MoleculeContainer()
    with m.edit():
        sids = [m.add_atom(e,
                           charge=0 if charges is None else charges[k],
                           implicit_h=None if hydrogens is None else hydrogens[k])
                for k, e in enumerate(_symbols(atoms))]
        for i, j, o in bonds:
            m.add_bond(sids[i], sids[j], o)
        for k, p in (parities or {}).items():
            m.set_parity(sids[k], p)
    return m, sids


def _trichloropentane(arms_alike, middle=1):
    """CC(Cl)C(Cl)C(Cl)C -- 2,3,4-trichloropentane, the pseudo-asymmetry case.

    The middle carbon is stereogenic exactly when the outer two have OPPOSITE configurations: then
    its arms are (R)- and (S)-1-chloroethyl, four different ligands, CIP's lowercase r/s. When the
    outer two match, a C2 axis swaps two identical arms and the middle centre is not stereogenic
    at all. `arms_alike` selects which molecule this is, by configuration and not by label.

    THE RAW PARITY LABELS ARE NOT COMPARABLE BETWEEN THE TWO ARMS. Parity is stored against each
    anchor's own CSR-ascending `refs`, and here those two orders are mirror images of each other:
    atom 1's refs are (methyl, Cl, C3, H) and atom 5's are (C3, Cl, methyl, H), which differ by
    one transposition. So EQUAL raw labels mean OPPOSITE configurations, which is why the
    stereogenic case below stores 1 and 1 rather than 1 and 2.
    `test_the_two_pseudo_asymmetry_records_differ_in_configuration_not_in_label` measures exactly
    that, by translating both parities into one common direction order.

    `middle` is the raw label stated on the MIDDLE carbon, and it changes no chemistry: which
    molecule this is, and therefore whether that centre is stereogenic, is decided by the two arms
    alone. It exists because the parity's raw VALUE is what feature word IV screens (the SEG_PARITY
    byte is 2 for odd and 1 for even), so a test about the feature words needs the middle sign to
    be odd while `arms_alike` keeps it unjustified.
    """
    atoms = 'CCClCClCClC'
    #        0 1 2   3 4   5 6   7
    bonds = [(0, 1, 1), (1, 2, 1), (1, 3, 1), (3, 4, 1),
             (3, 5, 1), (5, 6, 1), (5, 7, 1)]
    outer_b = 2 if arms_alike else 1
    return _mol(atoms=atoms, bonds=bonds, hydrogens=[3, 1, 0, 1, 0, 1, 0, 3],
                parities={1: 1, 3: middle, 5: outer_b})


# ----------------------------------------------------------------------------------------------
# the predicate
# ----------------------------------------------------------------------------------------------

def test_pseudo_asymmetric_centre_is_stereogenic_when_outer_differ():
    m, sids = _trichloropentane(arms_alike=False)
    anchors = {u['anchor'] for u in m.stereogenic_units()}
    assert sids[3] in anchors, 'the middle centre is stereogenic; V2 and RDKit both drop it'
    assert anchors == {sids[1], sids[3], sids[5]}


def test_pseudo_asymmetric_centre_is_not_stereogenic_when_outer_match():
    m, sids = _trichloropentane(arms_alike=True)
    anchors = {u['anchor'] for u in m.stereogenic_units()}
    assert sids[3] not in anchors, 'the C2 axis makes the middle centre non-stereogenic'
    assert anchors == {sids[1], sids[5]}


def test_the_two_pseudo_asymmetry_records_differ_in_configuration_not_in_label():
    """The discriminating property of the pair above, measured rather than asserted by name.

    Both arms are translated into the same direction order -- (methyl, Cl, inner carbon, H) -- so the
    two numbers are finally comparable. Differing means the arms are enantiomeric, which is the
    molecule whose middle carbon is stereogenic. This is also the test that would catch the labels
    being read as if they were comparable: swap the two records and it fails.
    """
    for arms_alike in (True, False):
        m, sids = _trichloropentane(arms_alike=arms_alike)
        left = m.translate_stereo(sids[1], (sids[0], sids[2], sids[3], None))
        right = m.translate_stereo(sids[5], (sids[7], sids[6], sids[3], None))
        assert left and right, 'both outer centres are configured'
        assert (left == right) is arms_alike
        assert (sids[3] in m.chiral_atoms()) is not arms_alike
        # and the automorphism that decides this MOVES an anchor: the two outer centres share an
        # orbit, so the witness maps one onto the other and only the stored labels can refute it
        orbits = m.automorphism_orbits()
        assert orbits[sids[1]] == orbits[sids[5]]


def test_asymmetric_molecule_marks_every_candidate():
    m, sids = _mol(atoms='CFClBr', hydrogens=[1, 0, 0, 0],
                   bonds=[(0, 1, 1), (0, 2, 1), (0, 3, 1)])
    assert len(m.stereogenic_units()) == 1
    assert m.is_asymmetric(), 'this is the shortcut that answered it'


def test_tartaric_acid_centres_stay_stereogenic():
    # HOOC-CH(OH)-CH(OH)-COOH
    atoms = 'COOCOCOCOO'
    #        0 1 2 3 4 5 6 7 8 9   -> C0 and C7 carboxyl, C3 and C5 the stereocentres
    bonds = [(0, 1, 2), (0, 2, 1), (0, 3, 1), (3, 4, 1), (3, 5, 1),
             (5, 6, 1), (5, 7, 1), (7, 8, 2), (7, 9, 1)]
    m, sids = _mol(atoms=atoms, bonds=bonds, hydrogens=[0, 0, 1, 1, 1, 1, 1, 0, 0, 1],
                   parities={3: 1, 5: 2})
    anchors = {u['anchor'] for u in m.stereogenic_units()}
    assert sids[3] in anchors and sids[5] in anchors
    # meso or not, both centres are stereogenic -- an automorphism exchanging them is not a witness
    # for either one, because a witness has to FIX the anchor it is asked about
    m2, sids2 = _mol(atoms=atoms, bonds=bonds, hydrogens=[0, 0, 1, 1, 1, 1, 1, 0, 0, 1])
    assert {u['anchor'] for u in m2.stereogenic_units()} == {sids2[3], sids2[5]}


def test_symmetric_arms_suppress_a_centre():
    # CH(CH3)(CH3)Cl -- two identical methyls, so never stereogenic
    m, sids = _mol(atoms='CCCCl', hydrogens=[1, 3, 3, 0],
                   bonds=[(0, 1, 1), (0, 2, 1), (0, 3, 1)])
    assert m.stereogenic_units() == []
    assert m.stereo_units(), 'the candidates are still emitted; only the flag is off'


def test_stereogenic_flag_appears_on_every_unit_dict():
    m, sids = _mol(atoms='CFClBr', hydrogens=[1, 0, 0, 0],
                   bonds=[(0, 1, 1), (0, 2, 1), (0, 3, 1)])
    assert m.stereo_units()[0]['stereogenic'] is True


def test_marking_is_idempotent():
    m, sids = _trichloropentane(arms_alike=False)
    assert m.stereogenic_units() == m.stereogenic_units()


# ----------------------------------------------------------------------------------------------
# the surface
# ----------------------------------------------------------------------------------------------

def test_is_chiral_is_true_for_a_labelled_centre_too():
    # bromochlorofluoromethane: sign it, and it stays chiral -- a label does not end the site
    m, sids = _mol(atoms='CFClBr', hydrogens=[1, 0, 0, 0],
                   bonds=[(0, 1, 1), (0, 2, 1), (0, 3, 1)])
    assert m.is_chiral(sids[0])
    with m.edit():
        m.set_parity(sids[0], 1)
    assert m.is_chiral(sids[0]), 'is_chiral must not mean "still needs a sign"'
    assert sids[0] in m.chiral_atoms()


def test_the_unsigned_subset_is_a_comprehension():
    m, sids = _mol(atoms='CFClBr', hydrogens=[1, 0, 0, 0],
                   bonds=[(0, 1, 1), (0, 2, 1), (0, 3, 1)])
    assert [s for s in m.chiral_atoms() if m.parity_of(s) == 0] == [sids[0]]
    with m.edit():
        m.set_parity(sids[0], 1)
    assert [s for s in m.chiral_atoms() if m.parity_of(s) == 0] == []


def test_chiral_bonds_keys_on_the_bond_not_the_anchor():
    # 2-butene: the stereogenic unit is the C=C bond, and no atom anchors one
    m, sids = _mol(atoms='CCCC', hydrogens=[3, 1, 1, 3],
                   bonds=[(0, 1, 1), (1, 2, 2), (2, 3, 1)])
    assert list(m.chiral_bonds()) == [tuple(sorted((sids[1], sids[2])))]
    assert m.chiral_atoms() == {}
    # the atom-side question still answers at whichever terminal holds the bits, which is the
    # reason the bond-side key exists
    assert m.is_chiral(sids[1]) and not m.is_chiral(sids[2])


def test_chiral_atoms_values_are_the_unit_dicts():
    m, sids = _mol(atoms='CFClBr', hydrogens=[1, 0, 0, 0],
                   bonds=[(0, 1, 1), (0, 2, 1), (0, 3, 1)])
    assert m.chiral_atoms()[sids[0]] == m.unit_of(sids[0])


def test_is_chiral_on_a_direction_atom_and_on_a_stranger():
    m, sids = _mol(atoms='CFClBr', hydrogens=[1, 0, 0, 0],
                   bonds=[(0, 1, 1), (0, 2, 1), (0, 3, 1)])
    assert not m.is_chiral(sids[1]), 'a direction anchors nothing'
    with pytest.raises(KeyError):
        m.is_chiral(max(sids) + 100)


def test_signing_the_outer_centres_reveals_the_middle_one():
    """chiral_atoms must re-perceive after an edit, not replay a cached answer."""
    m, sids = _trichloropentane(arms_alike=False)
    assert set(m.chiral_atoms()) == {sids[1], sids[3], sids[5]}
    m2, sids2 = _trichloropentane(arms_alike=True)
    assert set(m2.chiral_atoms()) == {sids2[1], sids2[5]}
    # and the same molecule with the outer signs stripped loses the middle centre too:
    # with nothing configured, the C2 axis is unbroken
    m3, sids3 = _mol(atoms='CCClCClCClC', hydrogens=[3, 1, 0, 1, 0, 1, 0, 3],
                     bonds=[(0, 1, 1), (1, 2, 1), (1, 3, 1), (3, 4, 1),
                            (3, 5, 1), (5, 6, 1), (5, 7, 1)])
    assert sids3[3] not in m3.chiral_atoms()
    with m3.edit():
        m3.set_parity(sids3[1], 1)
        m3.set_parity(sids3[5], 1)
    assert sids3[3] in m3.chiral_atoms(), 'perception must re-run against the new parities'


# ----------------------------------------------------------------------------------------------
# spec gate 3: mutually dependent units
# ----------------------------------------------------------------------------------------------

def test_mutually_dependent_ring_centres_are_both_stereogenic():
    """cis/trans-1,4-dimethylcyclohexane, CC1CCC(C)CC1.

    Neither ring CH is stereogenic alone -- the ring's mirror swaps its two arms -- but that same
    mirror is odd at the OTHER CH, so it is a witness for neither.  Two real diastereomers.
    Chython 2 gets this; a stored-parities-only predicate loses it.
    """
    m, sids = _mol(atoms='CCCCCCCC',                       # 0 Me, 1 CH, 2 3 CH2, 4 CH, 5 Me, 6 7 CH2
                   bonds=[(0, 1, 1), (1, 2, 1), (2, 3, 1), (3, 4, 1),
                          (4, 5, 1), (4, 6, 1), (6, 7, 1), (7, 1, 1)],
                   hydrogens=[3, 1, 2, 2, 1, 3, 2, 2])
    assert set(m.chiral_atoms()) == {sids[1], sids[4]}
    assert m.chiral_bonds() == {}
    # the discriminating property: nothing is configured and the group is non-trivial, so this is
    # decided by the SELF-LOOP clause -- the candidate witness fixes both anchors (it swaps each
    # one's pair of ring arms) and is therefore odd at the other unit as well as at this one
    assert not m.is_asymmetric(), 'the asymmetry shortcut must not be what answered this'
    assert all(m.parity_of(s) == 0 for s in sids), 'no stored parity is available to refute it'
    orbits = m.automorphism_orbits()
    assert orbits[sids[2]] == orbits[sids[7]] and orbits[sids[3]] == orbits[sids[6]], \
        'each CH has its two ring arms exchanged by an automorphism'


def test_mutually_dependent_exocyclic_double_bonds_are_both_stereogenic():
    """1,4-bis(ethylidene)cyclohexane, CC=C1CCC(CC1)=CC -- the same argument one kind up."""
    m, sids = _mol(atoms='CCCCCCCCCC',                     # 0..5 ring, 6/8 =CH, 7/9 Me
                   bonds=[(0, 1, 1), (1, 2, 1), (2, 3, 1), (3, 4, 1), (4, 5, 1), (5, 0, 1),
                          (0, 6, 2), (6, 7, 1), (3, 8, 2), (8, 9, 1)],
                   hydrogens=[0, 2, 2, 0, 2, 2, 1, 3, 1, 3])
    assert set(m.chiral_bonds()) == {tuple(sorted((sids[0], sids[6]))),
                                     tuple(sorted((sids[3], sids[8])))}
    assert m.chiral_atoms() == {}


def test_a_spiro_system_resolves_all_three_of_its_units():
    """CC=C1CCC2(CCC(C)CC2)CC1 -- ethylidene on ring A, methyl on ring B of spiro[5.5]undecane.

    The spiro atom, ring B's methyl-bearing CH and the ethylidene bond are mutually dependent, and
    all three come out stereogenic.  The spiro atom is the case the epic is measured on.
    """
    m, sids = _mol(atoms='CCCCCCCCCCCCCC',                 # 0 A1, 1 2 A, 3 spiro, 4 5 A,
                                                           # 6 7 B, 8 B-CH, 9 10 B, 11 =CH, 12 13 Me
                   bonds=[(0, 1, 1), (1, 2, 1), (2, 3, 1), (3, 4, 1), (4, 5, 1), (5, 0, 1),
                          (3, 6, 1), (6, 7, 1), (7, 8, 1), (8, 9, 1), (9, 10, 1), (10, 3, 1),
                          (0, 11, 2), (11, 12, 1), (8, 13, 1)],
                   hydrogens=[0, 2, 2, 0, 2, 2, 2, 2, 1, 2, 2, 1, 3, 3])
    assert set(m.chiral_atoms()) == {sids[3], sids[8]}
    assert set(m.chiral_bonds()) == {tuple(sorted((sids[0], sids[11])))}
    # the dependent set here is LARGER THAN A PAIR and spans two kinds -- three units decided
    # together, which is the case a pairwise argument cannot reach
    assert len(m.chiral_atoms()) + len(m.chiral_bonds()) == 3


def test_one_ring_marker_alone_is_not_stereogenic():
    """methylcyclohexane, CC1CCCCC1: the only candidate is the CH itself, so the arm swap IS a
    witness for it and nothing refuses it.  The near-miss that keeps the rule honest."""
    m, sids = _mol(atoms='CCCCCCC',
                   bonds=[(0, 1, 1), (1, 2, 1), (2, 3, 1), (3, 4, 1),
                          (4, 5, 1), (5, 6, 1), (6, 1, 1)],
                   hydrogens=[3, 1, 2, 2, 2, 2, 2])
    assert m.chiral_atoms() == {}


def test_one_exocyclic_double_bond_alone_is_not_stereogenic():
    """ethylidenecyclohexane, CC=C1CCCCC1 -- the cis/trans twin of the test above."""
    m, sids = _mol(atoms='CCCCCCCC',
                   bonds=[(0, 1, 1), (1, 2, 1), (2, 3, 1), (3, 4, 1), (4, 5, 1), (5, 0, 1),
                          (0, 6, 2), (6, 7, 1)],
                   hydrogens=[0, 2, 2, 2, 2, 2, 1, 3])
    assert m.chiral_bonds() == {}


def test_trimethylcyclohexane_answers_per_diastereomer():
    """1,3,5-trimethylcyclohexane, CC1CC(C)CC(C)C1: all eight labellings of its three centres.

    The ring has exactly two diastereomers. In the all-cis one, inverting ANY of the three centres
    gives the other diastereomer, so all three are stereogenic. In cis,cis,trans only the odd centre
    out is: inverting either of the matching pair reproduces the same compound through the ring flip.
    Two of the eight labellings are the all-cis isomer and its mirror image; the other six are
    cis,cis,trans, in each of the three rotations and both hands.

    This is where the epic's convention shows. V2 reports all three centres for every record,
    because it asks "could some isomer distinguish this site" -- and answers that even for a record
    with no configuration at all. This predicate asks about THE RECORD IN FRONT OF IT, which is what
    makes 2,3,4-trichloropentane come out right, and here it is what makes six of the eight answers
    a single centre. Cross-checked against V2 on fifty public compounds, this molecule and the
    sulfoxide (which V2's carbon-centric tetrahedral perception does not see at all) are the only
    two disagreements out of fifty.

    The raw labels are not comparable between the three anchors -- their `refs` orders are not
    aligned -- so which centre survives is read off the table rather than predicted from the signs.
    """
    bonds = [(0, 1, 1), (1, 2, 1), (2, 3, 1), (3, 4, 1), (4, 5, 1), (5, 0, 1),
             (0, 6, 1), (2, 7, 1), (4, 8, 1)]
    hydrogens = [1, 2, 1, 2, 1, 2, 3, 3, 3]
    expected = {(1, 1, 1): (0,), (1, 1, 2): (2,), (1, 2, 1): (4,), (1, 2, 2): (0, 2, 4),
                (2, 1, 1): (0, 2, 4), (2, 1, 2): (4,), (2, 2, 1): (2,), (2, 2, 2): (0,)}
    for labels, centres in expected.items():
        m, sids = _mol(atoms='C' * 9, bonds=bonds, hydrogens=hydrogens,
                       parities={0: labels[0], 2: labels[1], 4: labels[2]})
        assert set(m.chiral_atoms()) == {sids[c] for c in centres}, labels
    assert sum(1 for c in expected.values() if len(c) == 3) == 2, 'one isomer and its mirror'

    # and with nothing configured there is nothing to be inconsistent with, so no centre is marked
    m, sids = _mol(atoms='C' * 9, bonds=bonds, hydrogens=hydrogens)
    assert m.chiral_atoms() == {}


# ----------------------------------------------------------------------------------------------
# the two chemistry gates, from the reachable side
# ----------------------------------------------------------------------------------------------

def test_a_methyl_carbon_is_refused_for_its_repeated_hydrogens():
    """propane. Gate 0a: two directions in one list that are both nameless cannot be told apart.

    A methyl is a tetrahedral candidate whose refs are one carbon and three unnamed hydrogens, so
    the mask has three bits set and no graph reasoning is needed to refuse it.
    """
    m, sids = _mol(atoms='CCC', hydrogens=[3, 2, 3], bonds=[(0, 1, 1), (1, 2, 1)])
    masks = {u['anchor']: u['unnamed_mask'] for u in m.stereo_units()}
    assert masks[sids[0]] == 0b1110 and masks[sids[1]] == 0b1100
    assert m.chiral_atoms() == {}


def test_a_protic_sulfonium_is_refused_by_its_hydrogen_and_not_by_symmetry():
    """S-methyl-S-ethylsulfonium with an explicit hydrogen, [S+](C)(CC)[H].

    Gate 0b: a group-15 or -16 anchor carrying a hydrogen inverts too fast to hold a configuration.
    Every other signal here says stereogenic -- FOUR distinguishable directions (methyl, ethyl, the
    hydrogen and the lone pair), exactly ONE nameless slot so gate 0a cannot fire, and the whole
    molecule is provably asymmetric so the asymmetry shortcut would have marked it. The hydrogen is
    the only thing that can be refusing it.
    """
    m, sids = _mol(atoms='SCCCH', charges=[1, 0, 0, 0, 0], hydrogens=[0, 3, 2, 3, 0],
                   bonds=[(0, 1, 1), (0, 2, 1), (2, 3, 1), (0, 4, 1)])
    assert m.unit_of(sids[0])['unnamed_mask'] == 0b1000, 'one nameless slot: not gate 0a'
    assert m.is_asymmetric(), 'the asymmetry shortcut would have marked this'
    assert not m.is_chiral(sids[0])

    # the same sulfonium with the hydrogen replaced by a propyl IS stereogenic, so the gate is the
    # hydrogen and not the element
    m2, sids2 = _mol(atoms='SCCCCCC', charges=[1] + [0] * 6, hydrogens=[0, 3, 2, 3, 2, 2, 3],
                     bonds=[(0, 1, 1), (0, 2, 1), (2, 3, 1), (0, 4, 1), (4, 5, 1), (5, 6, 1)])
    assert m2.is_chiral(sids2[0])


def test_a_silicon_hydride_centre_stays_stereogenic():
    """O[SiH](CCC)C, a silane. Gate 0b is a GROUP membership test, not `element != 6`.

    Silicon is neither group 15 nor 16 and its hydride does not invert, so an Si-H stereocentre is
    real. Written the lazy way -- refuse any hydrogen-bearing non-carbon anchor -- this molecule
    would be lost.
    """
    m, sids = _mol(atoms='OSiCCCC', hydrogens=[1, 1, 2, 2, 3, 3],
                   bonds=[(0, 1, 1), (1, 2, 1), (2, 3, 1), (3, 4, 1), (1, 5, 1)])
    assert m.unit_of(sids[1])['unnamed_mask'] == 0b1000, 'the hydrogen is the nameless slot'
    assert set(m.chiral_atoms()) == {sids[1]}


def test_a_sulfoxide_keeps_its_lone_pair_centre():
    """methyl ethyl sulfoxide, C[S](=O)CC: a sulfur with three ligands and a lone pair."""
    m, sids = _mol(atoms='CSOCC', hydrogens=[3, 0, 0, 2, 3],
                   bonds=[(0, 1, 1), (1, 2, 2), (1, 3, 1), (3, 4, 1)])
    assert set(m.chiral_atoms()) == {sids[1]}
    assert m.unit_of(sids[1])['unnamed_mask'] == 0b1000, 'the lone pair is the nameless direction'


def test_a_quaternary_ammonium_with_four_different_chains_is_stereogenic():
    """methyl(ethyl)(propyl)(butyl)ammonium. A nitrogen with no hydrogen cannot invert.

    Four DIFFERENT chains, because the obvious record -- trimethyl(ethyl)ammonium -- has three
    interchangeable methyls and is correctly not stereogenic, which would have tested nothing.
    """
    m, sids = _mol(atoms='NCCCCCCCCCC', charges=[1] + [0] * 10,
                   hydrogens=[0, 3, 2, 3, 2, 2, 3, 2, 2, 2, 3],
                   bonds=[(0, 1, 1), (0, 2, 1), (2, 3, 1), (0, 4, 1), (4, 5, 1), (5, 6, 1),
                          (0, 7, 1), (7, 8, 1), (8, 9, 1), (9, 10, 1)])
    assert set(m.chiral_atoms()) == {sids[0]}
    assert m.unit_of(sids[0])['unnamed_mask'] == 0, 'all four directions are named'

    m2, sids2 = _mol(atoms='NCCCCC', charges=[1] + [0] * 5, hydrogens=[0, 3, 3, 3, 2, 3],
                     bonds=[(0, 1, 1), (0, 2, 1), (0, 3, 1), (0, 4, 1), (4, 5, 1)])
    assert m2.chiral_atoms() == {}, 'three identical methyls: not a centre'


def test_an_allene_is_stereogenic_and_its_terminals_are_not():
    """penta-2,3-diene, CC=C=CC. The allene axis is real; the methyls are gate 0a."""
    m, sids = _mol(atoms='CCCCC', hydrogens=[3, 1, 0, 1, 3],
                   bonds=[(0, 1, 1), (1, 2, 2), (2, 3, 2), (3, 4, 1)])
    assert set(m.chiral_atoms()) == {sids[2]}
    assert m.unit_of(sids[2])['kind'] == 2
    assert m.unit_of(sids[2])['unnamed_mask'] == 0b1010, 'one nameless slot per pair, so not 0a'


# ----------------------------------------------------------------------------------------------
# the machinery: the odd half of S4, and the asymmetry shortcut
# ----------------------------------------------------------------------------------------------

def test_the_odd_permutation_table_is_exactly_the_odd_half_of_s4():
    """The enumeration searches twelve permutations; these twelve and no others."""
    table = _core._odd_permutation_table()
    assert len(table) == 12 and len(set(table)) == 12
    for p in table:
        assert sorted(p) == [0, 1, 2, 3]
        assert _core._permutation_parity_probe(p) == 1
    # and it is the complement of the even half within all 24
    from itertools import permutations
    assert set(table) == {p for p in permutations(range(4))
                          if _core._permutation_parity_probe(p) == 1}


def test_permutation_parity_agrees_with_the_pair_decomposition():
    """Ruling F56: one parity function serves every kind.

    A bond kind's parity translates by pair decomposition -- swap within pair 0, swap within pair 1
    -- and the composed permutation's parity must agree, because the wholesale pair exchange
    (0 2)(1 3) is EVEN and therefore contributes nothing.
    """
    for p0 in (False, True):
        for p1 in (False, True):
            perm = [1, 0] if p0 else [0, 1]
            perm += [3, 2] if p1 else [2, 3]
            assert _core._permutation_parity_probe(tuple(perm)) == (p0 != p1)
            # the same two swaps after exchanging the pairs wholesale: same parity
            exchanged = tuple(perm[2:]) + tuple(perm[:2])
            assert _core._permutation_parity_probe(exchanged) == (p0 != p1)


def test_the_asymmetry_shortcut_is_not_what_decides_the_hard_cases():
    """Marking every candidate when the graph is asymmetric is only sound one way round.

    So the two sides are asserted apart: where `is_asymmetric()` is True everything is marked, and
    the molecules this task exists for answer with it False -- their verdicts come from the search.
    """
    asymmetric, a = _mol(atoms='OSiCCCC', hydrogens=[1, 1, 2, 2, 3, 3],
                         bonds=[(0, 1, 1), (1, 2, 1), (2, 3, 1), (3, 4, 1), (1, 5, 1)])
    assert asymmetric.is_asymmetric()
    for u in asymmetric.stereo_units():
        # every candidate that survived the two chemistry gates is marked
        if bin(u['unnamed_mask']).count('1') < 2:
            assert u['stereogenic'], u

    for m, _ in (_trichloropentane(arms_alike=True), _trichloropentane(arms_alike=False)):
        assert not m.is_asymmetric()
        assert m.chiral_atoms(), 'and the search still finds centres'


def test_a_terminal_exchanging_witness_unmarks_a_macrocyclic_double_bond():
    """A 20-ring whose only witness against its C=C EXCHANGES the alkene's two terminals.

    Four identical arcs `-CH<, NH, CH2, CH2, CH2-` close a 20-membered macrocycle, and one `C=C`
    bridges it from ring position 0 to 10 on one terminal and 5 to 15 on the other. The arcs are not
    palindromic, so no reflection survives and the group is the Z4 rotation, `|Aut| = 4`. All four
    rings through both terminals are 8-membered, so the unit is emitted rather than refused for ring
    size.

    THE QUARTER ROTATION IS THE WITNESS AND IT DOES NOT FIX EITHER TERMINAL. It carries terminal 0
    onto terminal 1, so it induces `[2, 3, 1, 0]` on the four direction slots -- a 4-cycle, ODD --
    and it therefore carries one cis/trans pairing onto the other: the two configurations are the
    same molecule and the bond names nothing. The half rotation fixes both terminals and induces
    `[1, 0, 3, 2]`, which is even and says nothing either way, so a predicate that only searches for
    anchor-fixing automorphisms sees no witness at all and marks the bond. That is ruling F61 and it
    is why phase 4 runs a second pinned search with the terminals pinned ACROSS (`_stereo.pxi`).

    Constructed by hand, and not a compound out of anyone's catalogue: the point is the symmetry, and
    a 20-ring is the smallest place it fits -- the same skeleton with three-atom arcs puts both
    terminals in a 6-ring, where perception refuses the unit for ring size before any of this.
    """
    # atoms 0 and 1 are the alkene; 2..21 are the macrocycle in ring order, each arc being a
    # substituted CH followed by NH, CH2, CH2, CH2
    atoms = 'CC' + 'CNCCC' * 4
    hydrogens = [0, 0] + [1, 1, 2, 2, 2] * 4
    bonds = [(0, 1, 2)] + [(2 + i, 2 + (i + 1) % 20, 1) for i in range(20)]
    bonds += [(0, 2 + 0, 1), (0, 2 + 10, 1), (1, 2 + 5, 1), (1, 2 + 15, 1)]
    m, sids = _mol(atoms=atoms, hydrogens=hydrogens, bonds=bonds)

    orbits = m.automorphism_orbits()
    assert orbits[sids[0]] == orbits[sids[1]], 'the two terminals are exchangeable, which is the case'
    assert sorted(len(r) for r in m.rings if sids[0] in r and sids[1] in r) == [8, 8, 8, 8], \
        'and the unit is emitted: no ring through the axis is too small'
    assert [u['kind'] for u in m.stereo_units()].count(1) == 1, 'so the C=C is a candidate'
    assert m.chiral_bonds() == {}, 'and the terminal exchange refutes it'

    # the same macrocycle with ONE arc lengthened by a carbon: the rotation is gone, and with it the
    # only witness, so the very same axis is stereogenic again. Without this the test would pass on
    # a predicate that simply refused every macrocyclic double bond.
    bonds2 = list(bonds)
    bonds2[bonds2.index((2, 3, 1))] = (2, 22, 1)        # splice a CH2 into the first arc
    bonds2.append((22, 3, 1))
    m2, sids2 = _mol(atoms=atoms + 'C', hydrogens=hydrogens + [2], bonds=bonds2)
    assert m2.is_asymmetric(), 'the longer arc breaks the Z4 rotation'
    assert list(m2.chiral_bonds()) == [tuple(sorted((sids2[0], sids2[1])))]


def test_an_automorphism_that_acts_evenly_is_not_a_witness():
    """FC(Cl)(Br)C(CH3)3. The group is non-trivial and permutes nothing the centre can see.

    NOT A PHASE 4 FIXTURE, AND IT CANNOT BE ONE. Phase 4 enumerates the twelve ODD rows and nothing
    else, so it never asks whether an even action is a witness -- that dimension of the predicate is
    structurally uncoverable from the outside, and this test covers the SHORTCUT that makes it moot:
    F, Cl, Br and C fall in four distinct refinement classes, so `_directions_separated` decides the
    headline centre at shortcut 2 and the search never sees it. (Measured with a probe that makes
    phase 4 refuse everything: this test does not fail. It fails under the opposite probe through its
    OTHER unit, the tert-butyl carbon, whose three interchangeable methyls do need a witness.)

    Rotating the tert-butyl's three methyls fixes the anchor and fixes each of its four directions
    pointwise, so it is EVEN there and refutes nothing. The four directions lying in four distinct
    orbits is what makes that visible.

    The three dimensions phase 4 does decide are varied on both sides by the fixtures in this file:
    seven need the search to find NO witness (both trichloropentane records, 1,4-dimethylcyclohexane,
    the bis-ethylidene pair, spiro[5.5], 1,3,5-trimethylcyclohexane) and ten need it to FIND one
    (arms-alike trichloropentane, methylcyclohexane, the single ring marker, the single exocyclic
    double bond, trimethylethylammonium, the tert-butyl carbon here, and the macrocycle above, which
    is the only one whose witness exchanges a unit's two terminals).
    """
    m, sids = _mol(atoms='CFClBrCCCC', hydrogens=[0, 0, 0, 0, 0, 3, 3, 3],
                   bonds=[(0, 1, 1), (0, 2, 1), (0, 3, 1), (0, 4, 1),
                          (4, 5, 1), (4, 6, 1), (4, 7, 1)])
    assert not m.is_asymmetric(), 'the three methyls are interchangeable'
    orbits = m.automorphism_orbits()
    assert len({orbits[sids[1]], orbits[sids[2]], orbits[sids[3]], orbits[sids[4]]}) == 4
    assert set(m.chiral_atoms()) == {sids[0]}


def _tetramethylcyclooctane(copies):
    """`copies` disconnected 1,3,5,7-tetramethylcyclooctanes in ONE record. Returns (mol, first sids).

    A plain public compound, and the shape the whole component restriction exists for: four
    stereocentres per copy, each copy an 8-ring with a symmetry group of its own.  Unrestricted, k
    identical components MULTIPLY those groups instead of adding them and the search runs out of budget
    at 60 atoms.
    """
    m = MoleculeContainer()
    firsts = []
    with m.edit():
        for _ in range(copies):
            ring = [m.add_atom('C', implicit_h=1 if i % 2 == 0 else 2) for i in range(8)]
            for i in range(8):
                m.add_bond(ring[i], ring[(i + 1) % 8], 1)
            for i in range(0, 8, 2):
                m.add_bond(ring[i], m.add_atom('C', implicit_h=3), 1)
            firsts.append(ring[0])
    return m, firsts


def test_k_identical_components_mark_k_times_one_copy_and_do_not_truncate():
    """The k-copies identity, and it is the sharpest verdict-invariance check this predicate has.

    The witness search is restricted to the anchor's own connected component, because a witness
    stabilizes its unit setwise and components are blocks of the automorphism group, so it stabilizes
    the component too.  Repeating one fragment as k separate components therefore cannot change any
    single unit's answer.  Unrestricted, the same repetition multiplies the group by k! and burns the
    budget: at k = 5 (60 atoms) the record comes back FLAGGED rather than decided.

    `marked == marked(k=1) * k` is the identity that cannot go vacuous.  `marked(k=1) == 4` is asserted
    separately, so a build that marked nothing would fail the left side and the right side both; and the
    component count is asserted, so a build that silently fused the copies could not pass either.

    THE TRUNCATION ASSERTION IS THE DELICATE ONE.  k = 5 (60 atoms), k = 6 (72 atoms) and ring6 k = 4
    (72 atoms) all decide in under a millisecond here, and Ruling F62's contract is untouched by that:
    no reader raises, the flag is still exposed, and the second-read path is pinned by
    `test_a_forged_truncation_word_survives_every_later_read`.  The negative control is a CONNECTED
    record: `test_a_connected_record_can_still_exhaust_the_budget` truncates at 54 atoms, so the flag
    is reachable and this assertion is not a claim that nothing truncates.
    """
    one, _ = _tetramethylcyclooctane(1)
    assert one.atom_count == 12 and one.connected_components_count == 1
    base = len(one.chiral_atoms())
    assert base == 4, 'a single copy is decided decisively, and at four -- not at zero'

    for copies in (2, 4, 5, 6):
        m, firsts = _tetramethylcyclooctane(copies)
        assert m.atom_count == 12 * copies
        assert m.connected_components_count == copies, 'k separate components, not one fused record'
        units = m.stereo_units()
        assert len(units) == 12 * copies
        assert sum(1 for u in units if u['stereogenic']) == base * copies, \
            'every copy answers exactly as it does alone: the restriction changes no verdict'
        assert len(m.chiral_atoms()) == base * copies
        assert m.chiral_bonds() == {}
        assert m.stereo_truncated is False, \
            'and it is now DECIDED, not flagged -- k=5 and k=6 used to exhaust the budget'
        assert all(m.is_chiral(f) is True for f in firsts), 'one marked centre per copy, by name'


def _methylated_macrocycle_explicit_h(size):
    """cyclo[CH(CH3)CH2]_size with EVERY hydrogen an explicit atom.  4*size heavy + 6*size hydrogens.

    Named for the family, not for one member: `size=4` is 1,3,5,7-tetramethylcyclooctane (the same
    public compound `_tetramethylcyclooctane` builds with implicit hydrogens), `size=6` is
    1,3,5,7,9,11-hexamethylcyclododecane, `size=8` is the cyclohexadecane.  Written with explicit
    hydrogens on purpose: an explicit hydrogen is a real atom the automorphism search must place, so
    three per methyl and two per CH2 give a CONNECTED record a large group without a second component.
    """
    m = MoleculeContainer()
    with m.edit():
        ring = [m.add_atom('C', implicit_h=0) for _ in range(2 * size)]
        for i in range(2 * size):
            m.add_bond(ring[i], ring[(i + 1) % (2 * size)], 1)
        for i in range(0, 2 * size, 2):                    # CH, one methyl, one hydrogen
            methyl = m.add_atom('C', implicit_h=0)
            m.add_bond(ring[i], methyl, 1)
            for _ in range(3):
                m.add_bond(methyl, m.add_atom('H', implicit_h=0), 1)
            m.add_bond(ring[i], m.add_atom('H', implicit_h=0), 1)
        for i in range(1, 2 * size, 2):                    # CH2
            for _ in range(2):
                m.add_bond(ring[i], m.add_atom('H', implicit_h=0), 1)
    return m


def test_a_connected_record_can_still_exhaust_the_budget():
    """The flag is reachable, so the assertions above are not a claim that nothing truncates.

    cyclo[CH(CH3)CH2]6 -- 1,3,5,7,9,11-hexamethylcyclododecane -- written with every hydrogen as an
    EXPLICIT atom, which is what a V3000 or PDB record hands over.  54 atoms in ONE component: each
    methyl's three hydrogens permute freely and each CH2's two do, so the group is a product of small
    factorials over a connected graph, and the component restriction cannot touch it.  Measured
    identical before and after the restriction: 13.3 ms and 13.5 ms, `stereo_truncated` true in both.

    THERE IS NO SIZE THRESHOLD HERE -- IT IS RING-SIZE PARITY, and a later task looking for a cheaper
    truncating record needs to know that.  Measured, size / atoms / marked / truncated:
    3/27/0/False, 4/36/4/False, 5/45/0/False, 6/54/6/True, 7/63/0/False, 8/72/8/True.  Odd sizes mark
    NOTHING and decide in under a tenth of a millisecond -- a witness turns up at once on an odd ring,
    so no search runs long -- while even sizes mark every ring CH, which means every one of those
    searches has to exhaust itself to prove there is no witness.  So 5 and 7 are LARGER than the control
    below and still do not truncate; only 6 and 8 do.

    Size 4 is the negative control BECAUSE IT IS EVEN: it runs the same shape of exhaustive search as
    size 6, marks all four of its ring CH, and merely does not run out of budget doing it.  That is what
    makes it a control rather than a coincidence -- this test fails if the budget is raised out of reach
    (6 stops truncating) and it fails if the searches stop happening at all (4 stops marking 4).
    """
    small = _methylated_macrocycle_explicit_h(4)
    assert small.atom_count == 36 and small.connected_components_count == 1
    assert small.stereo_truncated is False, 'the same EVEN shape, two rings smaller, still finishes'
    assert len(small.chiral_atoms()) == 4

    m = _methylated_macrocycle_explicit_h(6)
    assert m.atom_count == 54 and m.connected_components_count == 1
    units = m.stereo_units()                       # does not raise -- ruling F62
    assert m.stereo_truncated is True, 'the budget is still reachable, on a CONNECTED record'
    assert len(units) == 18
    assert sum(1 for u in units if u['stereogenic']) == 6, \
        'six ring CH: conservative here, and exact because a single copy decides at six'
    assert len(m.chiral_atoms()) == 6
    m.automorphism_orbits()                        # succeeds on the very record the stereo search cut


SU_TETRA, SU_CIS_TRANS = 0, 1        # `stereo_unit_t.kind`, as `stereo_units()` reports it

# Public compounds `_disjoint` assembles into multi-component records, with the answer each one gives
# ALONE written next to it -- that answer is what the two tests below require the same fragment to give
# inside a record with other components.  The first four are two PAIRS, chosen so that the members of a
# pair differ in the one property under test: both reach the witness search, one is REFUSED by it and
# the other survives.  Both halves are needed, because a leaked complement pin can only over-restrict,
# and over-restricting turns a refusal into a mark -- so the refusals are the sensitive half and the
# marks are the control.  The rest carry a second kind, a marked ATOM outside a ring, or no unit at all.
_FRAGMENTS = {
    # methylcyclohexane: the ring's arm swap IS a witness for the one CH -- 0 marks (tetrahedral)
    'methylcyclohexane': dict(atoms='CCCCCCC', hydrogens=[3, 1, 2, 2, 2, 2, 2],
                              bonds=[(0, 1, 1), (1, 2, 1), (2, 3, 1), (3, 4, 1), (4, 5, 1),
                                     (5, 6, 1), (6, 1, 1)]),
    # 1,4-dimethylcyclohexane: the same swap is odd at the OTHER CH -- 2 marks (tetrahedral)
    'dimethylcyclohexane': dict(atoms='CCCCCCCC', hydrogens=[3, 1, 2, 2, 1, 3, 2, 2],
                                bonds=[(0, 1, 1), (1, 2, 1), (2, 3, 1), (3, 4, 1), (4, 5, 1),
                                       (4, 6, 1), (6, 7, 1), (7, 1, 1)]),
    # ethylidenecyclohexane: the cis/trans twin of methylcyclohexane -- 0 marks (bond kind, sets = 2)
    'ethylidenecyclohexane': dict(atoms='CCCCCCCC', hydrogens=[0, 2, 2, 2, 2, 2, 1, 3],
                                  bonds=[(0, 1, 1), (1, 2, 1), (2, 3, 1), (3, 4, 1), (4, 5, 1),
                                         (5, 0, 1), (0, 6, 2), (6, 7, 1)]),
    # 1,4-bis(ethylidene)cyclohexane: the twin of dimethylcyclohexane -- 2 marks (bond kind)
    'bisethylidenecyclohexane': dict(atoms='CCCCCCCCCC', hydrogens=[0, 2, 2, 0, 2, 2, 1, 3, 1, 3],
                                     bonds=[(0, 1, 1), (1, 2, 1), (2, 3, 1), (3, 4, 1), (4, 5, 1),
                                            (5, 0, 1), (0, 6, 2), (6, 7, 1), (3, 8, 2), (8, 9, 1)]),
    # butan-2-ol: an acyclic centre, marked, and nothing about it involves a ring -- 1 mark on atom 1
    'butan2ol': dict(atoms='CCOCC', hydrogens=[3, 1, 1, 2, 3],
                     bonds=[(0, 1, 1), (1, 2, 1), (1, 3, 1), (3, 4, 1)]),
    # propan-2-ol: the same site with two identical arms -- a candidate, and refused.  0 marks
    'propan2ol': dict(atoms='CCOC', hydrogens=[3, 1, 1, 3], bonds=[(0, 1, 1), (1, 2, 1), (1, 3, 1)]),
    # but-2-ene: a marked CIS/TRANS unit with no ring anywhere -- 1 mark, on the bond 1=2
    'but2ene': dict(atoms='CCCC', hydrogens=[3, 1, 1, 3], bonds=[(0, 1, 1), (1, 2, 2), (2, 3, 1)]),
    # water: a component with NO unit at all, which is what a solvate or a counter-ion usually is
    'water': dict(atoms='O', hydrogens=[2], bonds=[]),
}

# the two pairs, and only they, carry the refused/surviving property the both-orders test asserts
_PAIRED = ('methylcyclohexane', 'dimethylcyclohexane',
           'ethylidenecyclohexane', 'bisethylidenecyclohexane')


def _disjoint(*names):
    """One record holding each named fragment as a separate component, in the order given.

    Returns (mol, blocks), where `blocks[b]` is fragment `b`'s stable ids in its own atom order.
    """
    m = MoleculeContainer()
    blocks = []
    with m.edit():
        for name in names:
            spec = _FRAGMENTS[name]
            sids = [m.add_atom(e, implicit_h=spec['hydrogens'][k])
                    for k, e in enumerate(_symbols(spec['atoms']))]
            for i, j, o in spec['bonds']:
                m.add_bond(sids[i], sids[j], o)
            blocks.append(sids)
    return m, blocks


def _unit_keys(m, blocks, marked_only):
    """Every unit as (block, kind, its named refs as within-block indices, unnamed count).

    ANCHOR-FREE AND REFS-ORDER-FREE, deliberately (ruling F101).  Which end of a bond kind holds the
    stated parity, and in which order its four refs sit, are artifacts of the slot order the record
    happens to have -- so a key built on either would compare two records' creation orders rather than
    their verdicts.  The named refs as a SORTED SET of within-block indices, the kind, and how many
    slots are nameless are the parts that mean something.  Sorted ints, never a frozenset: `<` on a
    frozenset is subset containment, not an order, and sorting a list of them is not a comparison.
    """
    where = {s: (b, i) for b, sids in enumerate(blocks) for i, s in enumerate(sids)}
    out = []
    for u in m.stereo_units():
        if marked_only and not u['stereogenic']:
            continue
        named = sorted(where[r][1] for r in u['refs'] if r is not None)
        # THE ANCHOR IS IN THE SET TOO, not just the refs: "U's component" is defined by the anchor,
        # and the restriction pins the complement of `comp[anchor]`, so an anchor in a different
        # component from its own directions would be the one shape that breaks the argument outright.
        block = {where[r][0] for r in u['refs'] if r is not None} | {where[u['anchor']][0]}
        assert len(block) == 1, 'no unit may straddle two components, anchor included'
        out.append((block.pop(), u['kind'], tuple(named), 4 - len(named)))
    return Counter(out)


def test_two_units_in_different_components_are_decided_the_same_in_both_orders():
    """The complement pins must not leak from one unit's component into the next unit's search.

    Phase 4 pins every atom OUTSIDE the anchor's component to itself, and the pin array is reused
    across units.  The next unit's component is a different set, so a complement pin left standing pins
    that unit's OWN component too -- the whole record becomes the identity, the identity is even, no
    witness is found, and the unit is wrongly MARKED.  Silent: no crash, no truncation, just a
    stereocentre reported on a molecule that does not have one.

    Both orders, because the order decides whether the leak is visible at all.  Its first victim is the
    SECOND component processed, so `methylcyclohexane + dimethylcyclohexane` -- refusal first -- passes
    even with the reset deleted, while the reverse fails.  Measured, not reasoned: with the reset
    deleted five of the eight records below disagree and three do not.

    The property under test is the only thing separating the pairs.  Each pair is two cyclohexanes of
    the same kind that differ in exactly one thing -- whether the ring's arm swap is odd at a second
    unit and therefore refutes itself -- so one member is refused by the search and the other survives
    it.  The refused member is the sensitive one and the surviving member is the control; a build that
    marked everything would fail on the first and a build that marked nothing would fail on the second.
    """
    alone = {}
    every = {}
    for name in _PAIRED:
        m, blocks = _disjoint(name)
        assert m.connected_components_count == 1
        alone[name] = _unit_keys(m, blocks, True)
        every[name] = _unit_keys(m, blocks, False)

    # the discriminating property, stated rather than assumed: within each pair one is refused by the
    # witness search and the other is not, and the two pairs cover both kinds the search handles
    assert sum(alone['methylcyclohexane'].values()) == 0
    assert sum(alone['dimethylcyclohexane'].values()) == 2
    assert sum(alone['ethylidenecyclohexane'].values()) == 0
    assert sum(alone['bisethylidenecyclohexane'].values()) == 2
    assert {k[1] for k in alone['dimethylcyclohexane']} == {SU_TETRA}, 'tetrahedral: one pinned search'
    assert {k[1] for k in alone['bisethylidenecyclohexane']} == {SU_CIS_TRANS}, \
        'cis/trans: two pinned searches, ruling F61'

    cases = [
        ('methylcyclohexane', 'dimethylcyclohexane'),
        ('dimethylcyclohexane', 'methylcyclohexane'),
        ('ethylidenecyclohexane', 'bisethylidenecyclohexane'),
        ('bisethylidenecyclohexane', 'ethylidenecyclohexane'),
        ('methylcyclohexane',) * 3,                       # identical components, all refused
        ('dimethylcyclohexane',) * 3,                     # identical components, all marked
        ('methylcyclohexane', 'dimethylcyclohexane', 'ethylidenecyclohexane',
         'bisethylidenecyclohexane'),
        ('bisethylidenecyclohexane', 'ethylidenecyclohexane', 'dimethylcyclohexane',
         'methylcyclohexane'),
    ]
    seen_empty = seen_marked = False
    for names in cases:
        m, blocks = _disjoint(*names)
        assert m.connected_components_count == len(names), 'separate components, not one fused record'

        # RULING F100: perception's unit SET is itself order-dependent where two bond-kind units
        # contest one atom, so the set this record actually has is asserted rather than assumed --
        # otherwise a verdict comparison could be comparing two different unit tables.
        expected_units = Counter()
        expected_marked = Counter()
        for b, name in enumerate(names):
            for (_, kind, named, nameless), c in every[name].items():
                expected_units[(b, kind, named, nameless)] += c
            for (_, kind, named, nameless), c in alone[name].items():
                expected_marked[(b, kind, named, nameless)] += c
        assert _unit_keys(m, blocks, False) == expected_units, \
            'the same candidates as the fragments have alone'

        assert _unit_keys(m, blocks, True) == expected_marked, \
            'and the same verdicts: a component decides on its own, in either order'
        seen_empty |= not expected_marked
        seen_marked |= bool(expected_marked)
    assert seen_empty and seen_marked, 'both a fully refused record and a marked one were compared'


# Salt-shaped records: several components, at least one of which carries NO stereogenic unit, which is
# what a salt, a solvate, a reagent mixture or one side of a reaction actually looks like.  The expected
# answers are written out as (block, index) and cross-checked against the same fragment alone, so a
# typo in the table fails and a build that marked nothing fails too.
_SALTS = [
    # THE MINIMAL WITNESS.  1,4-dimethylcyclohexane plus methylcyclohexane: two marks, both on the
    # dimethyl ring, and a leaked complement pin reports THREE by inventing one on the methyl ring.
    (('dimethylcyclohexane', 'methylcyclohexane'), {(0, 1), (0, 4)}, set()),
    (('methylcyclohexane', 'dimethylcyclohexane'), {(1, 1), (1, 4)}, set()),
    # an acyclic centre, then two components that answer nothing -- one of them with no unit at all
    (('butan2ol', 'propan2ol', 'water'), {(0, 1)}, set()),
    # and the same three reversed, so the record STARTS with a component that has no unit: `comp[0]`
    (('water', 'propan2ol', 'butan2ol'), {(2, 1)}, set()),
    # NOTHING IS STEREOGENIC ANYWHERE, three components, three kinds of refusal.  The most sensitive
    # row there is: under a pin leak every component after the first marks, so 0 becomes 3.
    (('methylcyclohexane', 'ethylidenecyclohexane', 'propan2ol'), set(), set()),
    # the bond kind on both sides of the question: one marked cis/trans, one refused
    (('but2ene', 'ethylidenecyclohexane'), set(), {(0, 1, 2)}),
    (('bisethylidenecyclohexane', 'ethylidenecyclohexane', 'methylcyclohexane'), set(),
     {(0, 0, 6), (0, 3, 8)}),
    # both kinds marked at once, in different components, with two inert components between and after
    (('dimethylcyclohexane', 'but2ene', 'methylcyclohexane', 'water'), {(0, 1), (0, 4)}, {(1, 1, 2)}),
]


def test_a_multi_component_record_answers_component_by_component():
    """A salt's stereocentres are its components' stereocentres, and the readers must say so BY NAME.

    The invariance is asserted on `chiral_atoms()` and `chiral_bonds()` -- the answers a caller
    actually gets -- and by identity rather than by count, because the defect this guards against
    invents a centre on a DIFFERENT molecule than the one that has any.  Phase 4 restricts each unit's
    witness search to its anchor's component by pinning the complement, and the pin array is reused
    across units; a complement pin left standing pins the next component to the identity, the identity
    is even, no witness is found, and the unit is marked.  Measured with that reset deleted: a record
    holding 1,4-dimethylcyclohexane and methylcyclohexane reports THREE chiral atoms where two is
    correct, and `stereo_truncated` is False, so nothing tells the caller the answer is a guess.

    Order is part of the coverage.  The leak's first victim is the second component processed, so a
    record whose refusals all come before its marks survives it.  The CONVERSE DOES NOT HOLD, and this
    is measured rather than reasoned: with the reset deleted four of the eight cases below disagree, and
    two of the four survivors do carry a mark before a refusal.  Whether a component's unit is sensitive
    at all depends on the component and not only on its position -- propan-2-ol and
    ethylidenecyclohexane are unaffected wherever they sit.  So every case that has any mark carries a
    companion with the order reversed instead of relying on a rule about ordering, and two cases mark
    nothing at all -- those are the sharpest, since under the leak every component after the first turns
    into a mark.

    Also the shape a counter-ion or a solvate really has: `water` contributes no unit whatsoever, and
    it appears both after the marked component and as the FIRST component of a record, which is the one
    place `comp[0]` is read.
    """
    seen_atom = seen_bond = seen_empty = False
    for names, atoms, bonds in _SALTS:
        m, blocks = _disjoint(*names)
        assert m.connected_components_count == len(names), 'separate components, not one fused record'
        assert m.stereo_truncated is False, 'decided, so these marks are proven and not conservative'

        # the table cross-checked against the fragments alone: the same marks, component by component
        expect_atoms = set()
        expect_bonds = set()
        for b, name in enumerate(names):
            one, one_blocks = _disjoint(name)
            index = {s: i for i, s in enumerate(one_blocks[0])}
            expect_atoms |= {(b, index[s]) for s in one.chiral_atoms()}
            expect_bonds |= {(b,) + tuple(sorted(index[s] for s in pair))
                             for pair in one.chiral_bonds()}
        assert (expect_atoms, expect_bonds) == (atoms, bonds), \
            'the stated answer is the one each fragment gives on its own'

        where = {s: (b, i) for b, sids in enumerate(blocks) for i, s in enumerate(sids)}
        assert {where[s] for s in m.chiral_atoms()} == atoms, \
            'exactly these atoms, in exactly these components'
        assert {(where[pair[0]][0],) + tuple(sorted(where[s][1] for s in pair))
                for pair in m.chiral_bonds()} == bonds, 'and exactly these bonds'
        # THE ATOM-SIDE READER AGREES FOR EVERY ATOM, not only for the ones expected to be marked --
        # except at the two ends of a marked bond, where one of them anchors that unit and answers
        # True.  WHICH one is a slot-order artifact (ruling F101), so both are skipped rather than
        # predicted; every other atom in the record must answer False.
        ends = {(b, i) for b, *pair in bonds for i in pair}
        for s, key in where.items():
            if key in ends:
                continue
            assert m.is_chiral(s) is (key in atoms), f'is_chiral disagrees at {key}'

        seen_atom |= bool(atoms)
        seen_bond |= bool(bonds)
        seen_empty |= not atoms and not bonds
    assert seen_atom and seen_bond and seen_empty, \
        'marked atoms, marked bonds and a record with neither were all compared'


def test_a_forged_truncation_word_survives_every_later_read():
    """The second-read path: perception is not re-run, so the header word is all that is left.

    The record above tests the search reaching truncation. This one tests what a table already built
    and already flagged does on every read after that -- the word must keep being reported, and the
    readers must keep answering rather than acquiring a raise. Forged, deliberately: a real truncating
    record would test the search again instead of the state.
    """
    m, sids = _mol(atoms='CCClCClCClC', hydrogens=[3, 1, 0, 1, 0, 1, 0, 3],
                   bonds=[(0, 1, 1), (1, 2, 1), (1, 3, 1), (3, 4, 1), (3, 5, 1), (5, 6, 1),
                          (5, 7, 1)])
    assert len(m.chiral_atoms()) == 2, 'the table builds and decides, to start with'
    assert m.stereo_truncated is False

    _core._stereo_forge_truncation(m)
    assert m.stereo_truncated is True
    assert len(m.stereo_units()) == 5, 'the five carbons with four directions each'
    assert len(m.stereogenic_units()) == 2, 'the conservative table still reads back'
    assert len(m.chiral_atoms()) == 2
    assert m.is_chiral(sids[1]) is True
    assert m.unit_of(sids[1])['anchor'] == sids[1]
    assert m.stereo_truncated is True, 'and the flag did not evaporate on the way'

    with m.edit():                          # an edit that CHANGES the record drops the derived
        m.add_bond(sids[0], m.add_atom('C', implicit_h=3), 1)
    assert m.stereo_truncated is False, 'segment, so the next read perceives afresh and decides again'
    assert len(m.chiral_atoms()) == 3


def _tetramethylcyclooctane_with_a_spare_chlorine(copies):
    """`_tetramethylcyclooctane`, but with an unbonded chlorine at slot 0. Returns (mol, firsts, Cl).

    Slot 0 is the point: it is below every ring atom, so bonding it to a ring CH -- with that carbon's
    implicit hydrogen stated away in the SAME edit, since the core derives none -- puts a new direction
    at the FRONT of an existing row and permutes it. That is the only way to reach a permutation at all
    (`add_atom` appends), and it is what makes this record test the re-basing arithmetic rather than the
    identity.
    """
    m = MoleculeContainer()
    firsts = []
    with m.edit():
        chlorine = m.add_atom('Cl')
        for _ in range(copies):
            ring = [m.add_atom('C', implicit_h=1 if i % 2 == 0 else 2) for i in range(8)]
            for i in range(8):
                m.add_bond(ring[i], ring[(i + 1) % 8], 1)
            for i in range(0, 8, 2):
                m.add_bond(ring[i], m.add_atom('C', implicit_h=3), 1)
            firsts.append(ring[0])
    return m, firsts, chlorine


def test_an_apply_that_re_bases_a_parity_leaves_the_marking_to_the_next_reader():
    """Ruling F70: the apply builds the table UNMARKED, and every mark still comes out right.

    `_harvest_parities` and `_replay_parities` read `kind`, `refs` and the unnamed nibble and never the
    SU_STEREOGENIC marks, so they go through `ensure_stereo_units_unmarked` and the budgeted symmetry
    search does not run inside the edit at all -- which is what took twenty trivial edits on this very
    record from 3449 ms back to 0.3 ms once a single parity was stored on it.  The table the replay
    leaves in the arena carries a zero in its "marked" header word, so the FIRST READER that wants
    `stereogenic` runs `mark_stereogenic` over that table and writes the truncation word then.

    THE EDIT HERE IS A GENUINE ODD PERMUTATION, deliberately: the chlorine at slot 0 enters the row in
    front of all three named directions and the hydrogen is stated away in the same edit, so the row goes
    `(ring, ring, methyl, implicitH) -> (Cl, ring, ring, methyl)` -- the cycle (3 0 1 2), three
    transpositions, odd -- and the stored bit MUST come out flipped. An earlier version of this test
    added a lone water molecule instead, whose permutation is the identity, so it would have passed
    against a replay that did nothing at all; the parity assertion below is what makes the pairing of a
    real re-base with a marks read non-vacuous, and that pairing is what F70's split makes worth pinning.

    The marks are then read AFTER the apply, and they must be the ones the same record gives without
    one: 20, exactly as in `test_k_identical_components_mark_k_times_one_copy_and_do_not_truncate` --
    four copies untouched at four marks each, plus four on the chlorinated copy, whose own symmetry the
    chlorine can only lower.

    THE TRUNCATION EXPECTATION MOVED FROM TRUE TO FALSE, and the new value is the right one.  This
    record is five separate components plus a chlorine, which is the k-copies shape the component
    restriction removed; the search now decides it in under a millisecond instead of exhausting the
    budget at 25 ms.  Nothing about F70's split moved with it -- the reader still runs
    `mark_stereogenic` on its own account, which is exactly why the 20 marks below are readable at all,
    and the word it writes is still read back on the second look.  A connected record still reaches
    truncation (`test_a_connected_record_can_still_exhaust_the_budget`), so this is a fact about this
    fixture and not about the flag.
    """
    m, firsts, chlorine = _tetramethylcyclooctane_with_a_spare_chlorine(5)
    centre = firsts[0]
    with m.edit():
        m.set_parity(centre, 1)
    before = m.unit_of(centre)              # also leaves a MARKED table in the arena, which the
    assert before['refs'][3] is None        # rebuild below drops -- the third of the review's attacks
    assert before['unnamed_mask'] == 0b1000, "slot 3 is the ring carbon's implicit hydrogen"
    assert m.parity_of(centre) == 1

    with m.edit():                          # the odd permutation, in one batch
        m.add_bond(centre, chlorine, 1)
        m.set_hydrogens(centre, 0)
    unit = m.unit_of(centre)
    assert unit is not None, 'four directions again, so still a unit'
    assert unit['refs'] == (chlorine,) + before['refs'][:3], 'the row gained a front element'
    assert m.parity_of(centre) == 2, 'an odd permutation of the directions flips the stored bit'

    # `unit_of` above was the first reader after the apply, so IT is what ran `mark_stereogenic` over the
    # unmarked table the replay left -- and over the RE-BASED parities, which `_stereo_consistent` reads.
    # Everything below is what that pass decided.
    assert m.stereo_truncated is False, "the search runs on the reader's account, and it DECIDES"
    units = m.stereo_units()
    assert len(units) == 60, 'the chlorine anchors nothing'
    assert sum(1 for u in units if u['stereogenic']) == 20, 'the same marks as without the apply'
    assert unit['stereogenic'] is True, 'including the re-based centre itself'
    assert len(m.chiral_atoms()) == 20
    assert m.is_chiral(centre) is True
    assert sum(1 for u in m.stereo_units() if u['stereogenic']) == 20, \
        'and a second read of the same table gives the same answer'
    assert m.stereo_truncated is False, 'with the word it wrote still reading back'


# ----------------------------------------------------------------------------------------------
# what this predicate does NOT reach
# ----------------------------------------------------------------------------------------------

def test_a_symmetric_ortho_biaryl_is_still_marked_on_a_kekule_record():
    """2,6-dichloro-2'-chlorobiphenyl: not an atropisomer, and this branch cannot see that.

    Both ortho positions of the left ring carry a chlorine, so turning that ring over reproduces the
    molecule and the axis has no second configuration. The refutation needs the ring turn to BE an
    automorphism -- and on a Kekule record it is not one, because the pivot's two ring bonds have
    different orders. The two ortho chlorines come out in different orbits, so no witness exists and
    the axis is marked.

    Nothing in this task can fix that: the arena stores Kekule orders only (an order-4 apply raises),
    so the symmetry simply is not in the graph. When aromatic bonds or a resonance-invariant
    refinement arrive, this test flips -- which is exactly what it is here to announce.
    `test_stereo_units.test_a_symmetric_ortho_pair_is_still_a_candidate_here` predicted
    `mark_stereogenic` would remove this candidate; it does not.
    """
    bonds = [(0, 1, 2), (1, 2, 1), (2, 3, 2), (3, 4, 1), (4, 5, 2), (5, 0, 1),
             (6, 7, 2), (7, 8, 1), (8, 9, 2), (9, 10, 1), (10, 11, 2), (11, 6, 1), (0, 6, 1)]
    m, sids = _mol(atoms='C' * 12 + 'ClClCl',
                   hydrogens=[0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0],
                   bonds=bonds + [(1, 12, 1), (5, 13, 1), (7, 14, 1)])
    orbits = m.automorphism_orbits()
    assert orbits[sids[1]] != orbits[sids[5]], 'the Kekule record breaks the ring turn'
    assert list(m.chiral_bonds()) == [tuple(sorted((sids[0], sids[6])))]

    # the genuine atropisomer, for contrast: 2-chloro-2'-fluorobiphenyl has no such symmetry to lose
    real, s = _mol(atoms='C' * 12 + 'ClF', hydrogens=[0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0],
                   bonds=bonds + [(1, 12, 1), (7, 13, 1)])
    assert list(real.chiral_bonds()) == [tuple(sorted((s[0], s[6])))]


# ----------------------------------------------------------------------------------------------
# The small-ring cut for cis/trans units (ruling F63, spec §4.6)
# ----------------------------------------------------------------------------------------------
# SU_MIN_STEREO_RING = 8 is the single threshold in `_stereo.pxi`.  A cis/trans unit whose two
# terminals share a ring smaller than that is not emitted by perception.  These tests pin the
# observable boundary so that editing SU_MIN_STEREO_RING moves something measurable.
#
# WHY THE CUT IS IN PERCEPTION AND NOT IN A LATER MARKING STAGE.  §4.6 specified a separate
# SU_UNREALIZABLE bit set after stereogenicity is decided.  That stage is provably empty on a
# Kekulé arena: it would call `_terminals_share_small_ring` on exactly the units that already
# survived `_terminals_share_small_ring` in perception's pass 2.  The reason the candidate rule
# is the only place the test can live is also structural: on a Kekulé arena an aromatic ring bond
# is order 2 and indistinguishable from a small-ring alkene, so admitting the unit would flood
# `stereo_units()` with junk candidates on every aromatic molecule.  The full argument is in the
# `_terminals_share_small_ring` comment and at the refusal site in `_stereo.pxi`.

def _cyclic_alkene(ring_size):
    """A single carbocycle of `ring_size` atoms with one ring double bond.

    Each ring carbon carries one implicit hydrogen so that the double-bond terminals each have two
    distinguishable directions (the ring neighbour and the H).  Without it each terminal has only
    one non-chain direction, and the unit's emission depends on the ring-size test rather than on
    the terminal-pair test, so the assertions below would not probe the right boundary.
    """
    atoms = 'C' * ring_size
    bonds = [(i, (i + 1) % ring_size, 1) for i in range(ring_size)]
    bonds[0] = (0, 1, 2)
    return _mol(atoms=atoms, bonds=bonds, hydrogens=[1] * ring_size)


def test_cyclohexene_double_bond_is_excluded_by_small_ring_cut():
    """Cyclohexene's double bond is not admitted: ring size 6 < SU_MIN_STEREO_RING (8).

    The unit is not in `stereo_units()` at all.  On a Kekulé arena this is the realizability
    answer -- the unit is genuinely stereogenic in the graph but not realizable in 3D, and
    conflating those two is the defect §4.6 calls out.  When the arena can distinguish aromatic
    bonds from Kekulé order 2, the cut can be relaxed and an SU_UNREALIZABLE mark becomes
    meaningful (value 2 in the flag nibble is free for that purpose).
    """
    m, sids = _cyclic_alkene(6)
    assert not [u for u in m.stereo_units() if u['kind'] == 1]


def test_cyclodecene_double_bond_is_admitted():
    """Cyclodecene's double bond is admitted: ring size 10 >= SU_MIN_STEREO_RING (8).

    The cis/trans unit IS in `stereo_units()` and is stereogenic (the 10-ring has no
    automorphism that acts oddly on the double bond's four direction slots).
    """
    m, sids = _cyclic_alkene(10)
    units = [u for u in m.stereo_units() if u['kind'] == 1]
    assert units, 'the cis/trans unit must be emitted for a 10-ring'
    assert any(u['stereogenic'] for u in units)


def test_acyclic_alkene_is_admitted_and_stereogenic():
    """But-2-ene: a plain acyclic cis/trans bond is admitted and stereogenic.

    No ring constraint applies; this is the base case showing the cut does not affect acyclic
    alkenes.  The unit is in `stereo_units()` and is stereogenic.
    """
    m, sids = _mol(atoms='CCCC', bonds=[(0, 1, 1), (1, 2, 2), (2, 3, 1)],
                   hydrogens=[3, 1, 1, 3])
    units = [u for u in m.stereo_units() if u['kind'] == 1]
    assert units, 'but-2-ene must have a cis/trans stereo unit'
    assert any(u['stereogenic'] for u in units)


def test_the_small_ring_cut_admits_at_threshold_and_refuses_below():
    """SU_MIN_STEREO_RING is 8: a 7-ring is refused, an 8-ring is admitted.

    This is the test that makes editing `SU_MIN_STEREO_RING` observable: changing it to 7 would
    make the cycloheptene assertion fail, and changing it to 9 would make the cyclooctene
    assertion fail.  Both assertions are present so that a single threshold change causes exactly
    one failure rather than making both pass or both fail silently.
    """
    # cycloheptene: 7-ring, refused (7 < 8)
    m7, _ = _cyclic_alkene(7)
    assert not [u for u in m7.stereo_units() if u['kind'] == 1], \
        'ring size 7 must be refused (7 < SU_MIN_STEREO_RING=8)'
    # cyclooctene: 8-ring, admitted (8 >= 8)
    m8, _ = _cyclic_alkene(8)
    assert [u for u in m8.stereo_units() if u['kind'] == 1], \
        'ring size 8 must be admitted (8 >= SU_MIN_STEREO_RING=8)'


def test_tetrahedral_centres_are_not_suppressed_by_ring_size():
    """A tetrahedral centre in a 5-ring is not suppressed by the small-ring cut.

    The cut applies to SU_CIS_TRANS only.  A cyclopentane ring with Cl on C0 and F on C1 gives
    two tetrahedral stereocentres; neither is affected by the ring size.

    Atom 0: heavy neighbors 1, 4, 5(Cl) -> degree 3, needs 1 implicit H for four directions.
    Atom 1: heavy neighbors 0, 2, 6(F)  -> degree 3, needs 1 implicit H for four directions.
    Atoms 2-4: degree 2, need 2 implicit H each.
    Cl(5), F(6): degree 1, no H.
    """
    m, sids = _mol(atoms='CCCCCClF',
                   bonds=[(0, 1, 1), (1, 2, 1), (2, 3, 1), (3, 4, 1), (4, 0, 1),
                          (0, 5, 1), (1, 6, 1)],
                   hydrogens=[1, 1, 2, 2, 2, 0, 0])
    # both tetrahedral candidates are stereogenic; stereogenic_units() = chiral_atoms() for TETRA
    tetra_stereo = {u['anchor'] for u in m.stereogenic_units() if u['kind'] == 0}
    assert sids[0] in tetra_stereo or sids[1] in tetra_stereo, \
        'at least one cyclopentane stereocentre must be stereogenic despite the 5-ring'


# ----------------------------------------------------------------------------------------------
# deferred validation: every stated configuration is stored, and judged once afterwards
# ----------------------------------------------------------------------------------------------

def test_three_stated_configurations_are_all_stored():
    # three stated parities, all three justified: `arms_alike=False` is the molecule whose middle
    # centre IS stereogenic, and a perception that drops it keeps only two of the three.
    m, sids = _trichloropentane(arms_alike=False)
    stored = [s for s in m.atom_numbers if m.parity_of(s)]
    assert len(stored) == 3, 'every stated configuration must survive'
    assert m.validate_stereo() == []


def test_a_parity_on_a_non_stereogenic_atom_is_reported():
    m, sids = _trichloropentane(arms_alike=True)
    # the C2 axis makes the middle centre non-stereogenic, so its stated parity is unjustified
    assert m.validate_stereo() == [sids[3]]


def test_a_reported_parity_is_cleared():
    """The clear happens, it is idempotent, and the stale unit table does NOT travel in the clone.

    NOT `unit_of(sid)['parity'] == 0`: that value is read live from SEG_PARITY, so the
    assertion cannot fail.  What the clone can carry stale is the SU_STEREOGENIC marks -- and on this
    record, and on every record anyone has built so far, the marks do not MOVE across the clear, for
    a reason that looks structural: witnesses compose.  If clearing B's pin admits a witness sigma
    that is odd on unit A, and B has its own witness tau that is even on A, then sigma*tau is a
    witness for A that was already admissible WITH B pinned, so A was never marked in the first
    place.

    THE ARGUMENT IS SOUND BUT NARROWER THAN THAT SENTENCE.  Producing tau needs A to be CONFIGURED
    as well -- an unconfigured marked unit puts no constraint on tau, so there is nothing to say tau
    is even on A -- and `_stereo_consistent` can reach its contradiction through a CHAIN of pins
    rather than through B's own ground pin, which the composition step does not model.  It is a
    heuristic backed by measurement, not a proof.  The measurements, all at this behaviour: nine
    parity assignments on 3-chloro-2,4-dimethylpentane (two identical isopropyl arms on a candidate
    middle carbon, built for this) leave the middle unmarked every time; a random sweep of 2607
    records carrying a parity on an unmarked unit moved no mark; and three independent review sweeps
    -- 9409 random clearing cases, 91650 on deliberately C2-symmetric skeletons with exhaustive
    parity assignments, and 800 on truncated records -- moved zero marks and zero `stereo_truncated`
    values.  So the comparison below is a TRIPWIRE: it will not fail today, and it is what catches
    the first record that breaks the argument.

    THE INVALIDATION ITSELF IS PINNED BY THE ARENA, which is measurable today.  `structure_clone`
    copies every derived cache, so the clone arrives carrying the stale table and `total_len` includes
    it; the invalidate RETIRES that block, which drops `total_len` by exactly one table, and the next
    reader allocates a fresh one, which puts it back.  So validation must SHRINK the arena and the
    following read must restore it to the byte.  Delete `structure_invalidate_stereo_units` and the
    shrink assertion fails with the two numbers equal, and nothing else in the suite does.  (Under v3
    the same call read the other way round -- the stale bytes were stranded inside the persistent
    buffer rather than tracked, so `total_len` could not fall and the evidence was the GROWTH on the
    following read.  Same call, same guarantee, opposite sign.)
    """
    m, sids = _trichloropentane(arms_alike=True)
    m.stereo_units()                                    # the table exists, marked, before the clear
    grown = m.total_len

    assert m.validate_stereo() == [sids[3]]
    assert m.parity_of(sids[3]) == 0
    assert m.total_len < grown, \
        'the stale table was retired in the clone, so its bytes left the arena'

    m.stereo_units()
    assert m.total_len == grown, \
        'the next reader had to derive a new table, of the same size as the one that was dropped'
    assert m.validate_stereo() == [], 'validation is idempotent'

    fresh = MoleculeContainer.from_bytes(m.to_bytes())   # derived segments are not serialised
    assert ([u['stereogenic'] for u in m.stereo_units()]
            == [u['stereogenic'] for u in fresh.stereo_units()]), \
        'the clone must not keep marks computed against the cleared parity'
    assert len(m.stereo_units()) == 5, 'and the comparison above was over a non-empty table'


def test_a_shared_arena_keeps_its_parity_when_the_original_validates():
    m, sids = _trichloropentane(arms_alike=True)
    other = m.copy()
    assert m.validate_stereo() == [sids[3]]
    assert other.parity_of(sids[3]) == 1, 'the clear must go into a clone (ruling F65)'
    assert not m.shares_arena_with(other)


def test_a_cleared_parity_does_not_come_back_through_a_round_trip():
    """`validate_stereo` clears the byte, and `to_bytes` carries the segment, so the clear persists.

    The parity is ODD because word IV bit 6 is set for parity 2 and only for parity 2, which is what
    `test_a_cleared_parity_leaves_the_feature_words_matching_a_round_trip` measures.  An even parity
    never sets bit 6 and would pass with the clear removed.
    """
    m, sids = _mol(atoms='CClCl', bonds=[(0, 1, 1), (0, 2, 1)], parities={0: 2})
    assert m.parity_of(sids[0]) == 2, 'odd: the parity byte is 2 going in'
    assert m.validate_stereo() == [sids[0]]
    again = MoleculeContainer.from_bytes(m.to_bytes())
    assert again.parity_of(sids[0]) == 0, 'the segment carries the cleared byte, and a zero stays 0'
    assert again.stereo_of(sids[0]) is False, 'and `stereo_of` reads the same byte'

    # the same clear on a record that also carries justified signs: those must survive it
    m, sids = _trichloropentane(arms_alike=True)
    assert m.validate_stereo() == [sids[3]]
    again = MoleculeContainer.from_bytes(m.to_bytes())
    assert again.parity_of(sids[3]) == 0
    assert [again.parity_of(s) for s in (sids[1], sids[5])] == [1, 2], \
        'and the two justified signs round-tripped untouched'


def test_a_parity_on_an_atom_with_no_candidate_is_reported():
    m, sids = _mol(atoms='CClCl', bonds=[(0, 1, 1), (0, 2, 1)], parities={0: 1})
    assert m.validate_stereo() == [sids[0]]


def test_setting_a_parity_never_raises():
    # mid-edit a container may legitimately carry parities nothing justifies yet
    m = MoleculeContainer()
    with m.edit():
        c = m.add_atom('C')
        m.set_parity(c, 1)
    assert m.parity_of(c) == 1
    assert m.validate_stereo() == [c]


def test_validation_is_clean_for_a_plain_stereocentre():
    # the hydrogen count is not optional: the core never derives one, so without it this carbon has
    # three directions, anchors no unit, and its parity would be reported
    m, sids = _mol(atoms='CFClBr', hydrogens=[1, 0, 0, 0],
                   bonds=[(0, 1, 1), (0, 2, 1), (0, 3, 1)], parities={0: 1})
    other = m.copy()
    assert m.validate_stereo() == []
    assert m.shares_arena_with(other), 'an empty report is a pure read: no clone, no _gen bump'
    assert m.generation == other.generation


def test_the_report_is_ascending_by_stable_id():
    """Order is a promise, not an accident of the slot walk.  `remap` is what separates the two.

    Slots are walked ascending and the ids are sorted afterwards, so a molecule whose stable ids
    descend against its slots reports the two orders differently.  Without the sort this list would
    come back as [99, 7] -- the slot order -- and a caller comparing it against ids of its own would be
    reading a walk order as a report.
    """
    m, sids = _mol(atoms='CClClCClCl',
                   bonds=[(0, 1, 1), (0, 2, 1), (3, 4, 1), (3, 5, 1)],
                   parities={0: 1, 3: 2})
    m.remap({sids[0]: 99, sids[3]: 7})
    assert [m.parity_of(s) for s in (99, 7)] == [1, 2], 'both signs survived the relabelling'
    assert m.validate_stereo() == [7, 99], 'ascending by stable id, not by slot'


def test_a_truncated_record_still_validates_clean():
    """Decision 5 marks what it could not settle, so a parity on a truncated candidate survives.

    The natural worry is the opposite: that a search which ran out of budget makes this function
    delete real input.  It cannot, because the conservative direction of the approximation is to MARK
    -- and a marked unit justifies its parity.  Measured on a record that really does truncate, not
    on a forged word.

    THE RECORD IS CONNECTED, and it has to be: the witness search is restricted to the anchor's own
    component, so a fixture of k identical copies decides instead of truncating.  The
    explicit-hydrogen ring truncates on its own -- 54 atoms, 13 ms.
    """
    m = _methylated_macrocycle_explicit_h(6)
    marked = sorted(u['anchor'] for u in m.stereogenic_units())
    assert len(marked) == 6, 'six ring CH, and the parities below go on two of them'
    with m.edit():
        m.set_parity(marked[0], 1)
        m.set_parity(marked[1], 2)
    assert m.stereo_truncated is True, 'the record still exhausts the budget with two signs on it'
    again = {u['anchor'] for u in m.stereogenic_units()}
    assert {marked[0], marked[1]} <= again, 'both signs sit on marked units'
    other = m.copy()
    assert m.validate_stereo() == [], 'a conservative mark keeps input; it never discards it'
    assert m.shares_arena_with(other)


# ----------------------------------------------------------------------------------------------
# ruling F78: every parity writer outside `rebuild_derived` maintains feature word IV
# ----------------------------------------------------------------------------------------------
# Feature word IV screens the SEG_PARITY byte at its bit 6, and `features_of()` and
# `_union_feature_words` hand those words to Python verbatim.  Two writers touch a parity after
# `rebuild_derived` has already built the words -- `validate_stereo`'s clear in the clone and the
# apply's `_replay_parities` -- and both went in forgetting the words, so a cleared or dropped sign
# stayed visible in the union row and disagreed with `parity_of` on the same molecule.  A
# `to_bytes`/`from_bytes` round trip re-derives from scratch, so it is the oracle: whatever a writer
# leaves behind must equal what the round trip computes, word for word.  Both tests use an ODD parity;
# an even one never sets bit 6 and would pass with the maintenance removed.

def _features_against_a_round_trip(m):
    """Which atoms' feature words disagree with a `from_bytes` rebuild of the same molecule."""
    fresh = MoleculeContainer.from_bytes(m.to_bytes())
    return [s for s in m.atom_numbers if m.features_of(s) != fresh.features_of(s)], fresh._union_feature_words


def test_a_cleared_parity_leaves_the_feature_words_matching_a_round_trip():
    """`validate_stereo`'s clear re-bases word IV, so the words cannot outlive the sign they screen.

    Before this was maintained: `parity_of` 0 and `stereo_of` False, but `features_of()[3]` bit 6 and
    `_union_feature_words[3]` bit 6 both still 1, against a round trip's 0.

    TWO RECORDS, and the second one is the union row's only pin.  A record whose ONLY odd sign is the
    one being cleared cannot tell a rebuilt union row from an un-ORed one -- both give 0 -- so the
    second half keeps an odd sign on an atom that is NOT cleared, where the two answers differ: 1 for
    the rebuild, 0 for `feat[3] &= ~(1 << 6)`, which is the exact mistake the rebuild exists to avoid
    and which the first half passes under.
    """
    m, sids = _mol(atoms='CClCl', bonds=[(0, 1, 1), (0, 2, 1)], parities={0: 2})
    assert (m.features_of(sids[0])[3] >> 6) & 1 == 1, 'the odd sign is in word IV to begin with'
    assert m.validate_stereo() == [sids[0]], 'and it is the sign that gets cleared'
    assert m.parity_of(sids[0]) == 0
    assert (m.features_of(sids[0])[3] >> 6) & 1 == 0, 'so it is out of word IV afterwards'
    disagree, fresh_union_row = _features_against_a_round_trip(m)
    assert len(m.atom_numbers) == 3, 'the comparison below ran over three atoms, not over none'
    assert disagree == [], 'every atom\'s four words equal what a fresh derivation computes'
    assert m._union_feature_words == fresh_union_row, 'including the union row'

    # THE SECOND CONFIGURED ARM IS LOAD-BEARING: atom 5 keeps an odd sign through the clear, so bit 6
    # must still be set in the union row afterwards, and un-ORing it out of `feat[3]` -- instead of
    # rebuilding the row from the atom words, which is what an OR forces -- is then measurably wrong.
    m, sids = _trichloropentane(arms_alike=True, middle=2)
    assert [m.parity_of(s) for s in (sids[1], sids[3], sids[5])] == [1, 2, 2], \
        'the middle sign is odd and unjustified; arm 5 carries an odd sign of its own'
    assert (m._union_feature_words[3] >> 6) & 1 == 1
    assert m.validate_stereo() == [sids[3]], 'only the middle sign is cleared'
    assert m.parity_of(sids[3]) == 0 and m.parity_of(sids[5]) == 2
    assert (m.features_of(sids[3])[3] >> 6) & 1 == 0, 'the cleared atom loses its own bit 6'
    assert (m._union_feature_words[3] >> 6) & 1 == 1, \
        'but the union row keeps it: atom 5 still owns that bit, and an OR cannot be un-ORed'
    disagree, fresh_union_row = _features_against_a_round_trip(m)
    assert len(m.atom_numbers) == 8, 'the comparison below ran over eight atoms, not over none'
    assert disagree == [], 'every atom\'s four words still equal a fresh derivation'
    assert m._union_feature_words == fresh_union_row, 'and so does the union row, on all four words'


def test_a_dropped_parity_leaves_the_feature_words_matching_a_round_trip():
    """The apply's re-base maintains word IV through the same helper, on the same oracle.

    `_replay_parities` runs AFTER `rebuild_derived`, on every edit that re-bases or drops a sign --
    not only when a caller asks for validation.  Deleting one of the four directions of a stereocentre
    destroys the frame, `rebase_parity` answers RB_DROP, and the parity byte is cleared.
    Unmaintained, word IV bit 6 stays at 1 against a round trip's 0.
    """
    m, sids = _mol(atoms='CFClBr', hydrogens=[1, 0, 0, 0],
                   bonds=[(0, 1, 1), (0, 2, 1), (0, 3, 1)], parities={0: 2})
    assert m.validate_stereo() == [], 'a justified odd sign on a plain stereocentre'
    assert (m.features_of(sids[0])[3] >> 6) & 1 == 1
    with m.edit():
        m.delete_atom(sids[3])
    assert m.parity_of(sids[0]) == 0, 'the frame lost a direction, so the apply dropped the sign'
    assert (m.features_of(sids[0])[3] >> 6) & 1 == 0, 'and word IV lost it with the byte'
    disagree, fresh_union_row = _features_against_a_round_trip(m)
    assert len(m.atom_numbers) == 3, 'the comparison below ran over three atoms, not over none'
    assert disagree == [], 'every atom\'s four words equal what a fresh derivation computes'
    assert m._union_feature_words == fresh_union_row
