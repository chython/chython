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
"""Acceptance criteria for the stereo stack.

This file is the *contract*, not the unit tests.  `test_stereo_units.py`,
`test_stereo_perception.py` and `test_stereo_parity.py` pin the arena's internals; what is pinned
here is the set of behaviours the stereo epic was told to preserve or to deliver, each expressed
through the public surface only, so that a reimplementation of the interior cannot quietly lose one.

Five groups, in the order the epic states them:

1. **Amide** -- the deliberate handling.  An amide C-N carries restricted rotation but is NOT a
   cis/trans unit, because the discriminator is BOND ORDER and not sp2 character.  Anything that
   starts perceiving hybridisation instead of order will break exactly these tests.
2. **Quaternary ammonium and sulfimides** -- the two named coverage gaps.  Both are tetrahedral
   centres that no element test reaches; the discriminator is the count of DIRECTIONS.
3. **Allene / cumulene** -- the odd/even ladder, which is the regression surface.
4. **Spiro** -- the other regression surface.
5. **Retranslation and fidelity** -- that a parity is re-expressed against a caller's own direction
   order rather than copied as a raw bit, and that input survives being wrong.

Every molecule here is a published, commonplace compound.  Nothing proprietary, nothing invented.
"""
import pytest

from chython.core import (read_smiles, molecule_to_inchi, inchi_library_loaded,
                          SU_TETRA, SU_CIS_TRANS, SU_ALLENE, SU_ATROPISOMER)


# --------------------------------------------------------------------------------------------
# helpers.  Deliberately thin: a test that has to be read alongside a clever helper is a test
# that does not document anything.

def _units(smi, kind=None):
    """Every perceived unit of `smi`, optionally filtered to one kind."""
    units = read_smiles(smi).stereo_units()
    return [u for u in units if kind is None or u['kind'] == kind]


def _stereogenic(smi, kind=None):
    """Only the units a molecule can actually hold two configurations of."""
    units = read_smiles(smi).stereogenic_units()
    return [u for u in units if kind is None or u['kind'] == kind]


def _kinds(smi):
    return {u['kind'] for u in read_smiles(smi).stereogenic_units()}


def _roundtrip(smi):
    """(first written form, second written form).  Equal is the invariant; see group 5."""
    first = read_smiles(smi).smiles
    return first, read_smiles(first).smiles


# ============================================================================================
# 1. AMIDE.  The epic's words: "amide stereo handling must be preserved -- it is deliberate, not
#    an accident."  This group says what that handling IS, so that the deliberateness survives a
#    rewrite by someone who was not told.
# ============================================================================================

# Restricted rotation is a *conformational* fact and the perception pass is a *constitutional*
# one.  An amide C-N is a single bond, so it is not a cumulene link, so no unit is spelled across
# it -- and that is the answer we want even though the barrier is real, because a cis/trans label
# on it would claim a configuration that ordinary chemistry interconverts at room temperature.
@pytest.mark.parametrize('name,smi', [
    ('N-methylacetamide',      'CC(=O)NC'),
    ('N,N-dimethylacetamide',  'CC(=O)N(C)C'),
    ('N-ethyl-N-methyl amide', 'CC(=O)N(C)CC'),
    ('acetanilide',            'CC(=O)Nc1ccccc1'),
    ('urea, N,N-disubst.',     'CNC(=O)NC'),
    ('methyl carbamate',       'CNC(=O)OC'),
    ('thioacetamide, N-methyl', 'CC(=S)NC'),
    ('formamide',              'NC=O'),
])
def test_an_amide_bond_is_never_a_cis_trans_unit(name, smi):
    """No planar unit anywhere in an amide, thioamide, urea or carbamate.

    The C-N is order 1 so it is not a chain link, and the C=O is order 2 but its oxygen terminal
    has no second substituent to be arranged against -- so neither bond can carry one.  Asserting
    on the KIND rather than on a count keeps the test honest when a molecule also has a
    tetrahedral candidate (the N-methyl carbon of `CC(=O)N(C)CC`, for instance).
    """
    assert SU_CIS_TRANS not in _kinds(smi), f'{name} grew a cis/trans unit'


# The contrast that proves the rule is about order and not about sp2: swap the amide's C-N for a
# C=N and the very same nitrogen becomes a planar unit.  If these two groups ever agree, the
# discriminator has drifted to hybridisation.
@pytest.mark.parametrize('name,smi', [
    ('acetaldoxime',        'CC=NO'),
    ('butan-2-one oxime',   'CCC(C)=NO'),
    ('benzaldehyde oxime',  'c1ccccc1C=NO'),
    ('N-methyl imine',      'CC=NC'),
    ('acetaldehyde hydrazone', 'CC=NN'),
])
def test_a_carbon_nitrogen_double_bond_is_a_cis_trans_unit(name, smi):
    """C=N is planar-stereogenic; C-N is not.  Order, not sp2 character.

    Butan-2-one oxime and not ACETONE oxime: acetone's two methyls make its carbon terminal
    symmetric, so that molecule is a candidate and not stereogenic -- correctly, and it would
    make this test assert the opposite of what it means to.
    """
    assert SU_CIS_TRANS in _kinds(smi), f'{name} lost its C=N unit'


def test_the_amide_carbonyl_is_not_a_unit_for_want_of_a_second_substituent():
    """Stated separately from the parametrised sweep because the REASON differs.

    The C-N is excluded by bond order.  The C=O is a genuine order-2 chain link and is excluded
    by something else entirely: its oxygen end has one neighbour, so there is no pair of
    directions at that terminal to arrange.  A rewrite that fixed the order test but forgot the
    terminal-substituent test would pass the sweep above and fail here.
    """
    # a ketone, where the only order-2 bond in the molecule is the carbonyl
    assert _stereogenic('CC(=O)CC', SU_CIS_TRANS) == []
    assert _stereogenic('CC=O', SU_CIS_TRANS) == []
    # and the same carbonyl inside an amide
    assert _stereogenic('CC(=O)NC', SU_CIS_TRANS) == []


def test_an_amide_nitrogen_is_not_a_tetrahedral_centre():
    """A trisubstituted amide N has three named directions and a lone pair, which is the same
    direction count as a sulfoxide -- but nitrogen inverts, and the arena does not grant a
    neutral group-15 lone pair a slot the way it grants sulfur's.  N,N-dimethylacetamide has no
    tetrahedral unit on its nitrogen.

    Contrast with `test_a_quaternary_ammonium_with_four_different_groups_is_stereogenic`: the
    difference is the charge, which is what removes the lone pair and stops the inversion.
    """
    assert _stereogenic('CC(=O)N(C)CC', SU_TETRA) == []
    assert _stereogenic('CC(=O)N(C)C', SU_TETRA) == []


# ============================================================================================
# 2. QUATERNARY AMMONIUM AND SULFIMIDES -- the two coverage gaps the epic named.  A candidate pass
#    opening with `atom == C` reaches neither, as V2's does; the arena's counts DIRECTIONS, so both
#    arrive without a special case.  These tests are what stop an element test creeping back in.
# ============================================================================================

def test_a_quaternary_ammonium_with_four_different_groups_is_stereogenic():
    """N-butyl-N-ethyl-N-methylpropan-1-aminium -- four different chains on a cationic nitrogen.

    There is no lone pair to invert through, so unlike a neutral amine this centre is
    configurationally stable and genuinely resolvable.
    """
    units = _stereogenic('CCCC[N+](C)(CC)CCC', SU_TETRA)
    assert len(units) == 1, 'expected exactly the ammonium nitrogen'


def test_a_quaternary_ammonium_holds_and_reports_a_parity():
    """The centre is not merely perceived -- a stated configuration lands on it and comes back.

    Both tags, so that the test cannot pass by always answering the same way, and the two must
    disagree: they are enantiomers.
    """
    at = [u['parity'] for u in _units('CCCC[N@+](C)(CC)CCC', SU_TETRA) if u['stereogenic']]
    atat = [u['parity'] for u in _units('CCCC[N@@+](C)(CC)CCC', SU_TETRA) if u['stereogenic']]
    assert at and atat
    assert at[0] != 0 and atat[0] != 0, 'a stated ammonium configuration was not stored'
    assert at[0] != atat[0], 'the two ammonium enantiomers collapsed onto one parity'


def test_a_quaternary_ammonium_survives_a_smiles_round_trip():
    first, second = _roundtrip('CCCC[N@+](C)(CC)CCC')
    assert first == second
    assert '[N@' in first, 'the ammonium configuration was dropped on the way out'


def test_a_repeated_group_makes_an_ammonium_not_stereogenic():
    """Two ethyls on the nitrogen and the mirror is an automorphism, so there is nothing to name.

    Perception still emits the CANDIDATE -- four directions is four directions -- and it is the
    stereogenicity pass that removes it.  Keeping these two questions apart is the whole design.
    """
    assert _stereogenic('CC[N+](C)(CC)CCC', SU_TETRA) == []
    # the candidate is nonetheless emitted: exactly one atom in this molecule has four heavy
    # neighbours and it is the nitrogen, so a four-named-direction candidate must exist for it
    candidates = [u for u in _units('CC[N+](C)(CC)CCC', SU_TETRA)
                  if all(r is not None for r in u['refs'])]
    assert len(candidates) == 1, 'the ammonium candidate was refused at perception, not at ' \
                                 'stereogenicity -- the two passes have been conflated'


def test_a_neutral_tertiary_amine_is_not_stereogenic():
    """The counterpart that must NOT change when the ammonium starts working.

    N-ethyl-N-methylpropan-1-amine has three different substituents and a lone pair, which looks
    like a centre and is not one: pyramidal inversion at neutral nitrogen is fast, and the arena
    declines to grant the lone pair a direction.  A fix for R4N+ that reached this molecule too
    would be a regression.
    """
    assert _stereogenic('CCN(C)CCC', SU_TETRA) == []


def test_a_sulfilimine_is_a_stereogenic_tetrahedral_centre():
    """S-ethyl-S-methyl-N-methylsulfilimine.

    Sulfur carries two carbons, an imine nitrogen and a lone pair.  That is one tetrahedral unit, and
    the S=N is NOT a cis/trans unit -- the sulfur is a tetrahedral centre, not a cumulene terminal.
    """
    assert len(_stereogenic('CC[S](C)=NC', SU_TETRA)) == 1
    assert _stereogenic('CC[S](C)=NC', SU_CIS_TRANS) == []


def test_a_sulfilimine_holds_and_reports_a_parity():
    at = [u['parity'] for u in _units('CC[S@](C)=NC', SU_TETRA) if u['stereogenic']]
    atat = [u['parity'] for u in _units('CC[S@@](C)=NC', SU_TETRA) if u['stereogenic']]
    assert at and atat and at[0] and atat[0]
    assert at[0] != atat[0], 'the two sulfilimine enantiomers collapsed onto one parity'


def test_a_sulfilimine_survives_a_smiles_round_trip():
    first, second = _roundtrip('CC[S@](C)=NC')
    assert first == second
    assert '[S@' in first


def test_a_sulfoxide_is_stereogenic_and_uses_the_same_lone_pair_slot():
    """Methyl phenyl sulfoxide -- the compound the sulfimide machinery shares its mechanism with.

    Pinned alongside the sulfimide because both depend on sulfur being granted exactly ONE lone
    pair direction: two pairs would be the same direction twice and would make every sulfone a
    centre.  Dimethyl sulfone is the negative control.
    """
    assert len(_stereogenic('C[S](=O)c1ccccc1', SU_TETRA)) == 1
    first, second = _roundtrip('C[S@](=O)c1ccccc1')
    assert first == second and '[S@' in first
    # a sulfone has no lone pair left, so it is not a centre however its substituents differ
    assert _stereogenic('CC[S](=O)(=O)C', SU_TETRA) == []


# ============================================================================================
# 3. ALLENE / CUMULENE.  The odd/even ladder is the single rule, and it is the thing most easily
#    lost by a reimplementation that special-cases three-atom allenes.
# ============================================================================================

# Atom count along the chain decides the KIND, with no per-length special case: odd -> axial,
# anchored on the middle atom; even -> planar, anchored on a terminal.  Walking the ladder is the
# test, because a three-atom-only implementation passes the allene row and fails the rest.
@pytest.mark.parametrize('name,smi,kind', [
    ('2-butene,      2 atoms', 'CC=CC',        SU_CIS_TRANS),
    ('2,3-pentadiene, 3 atoms', 'CC=C=CC',      SU_ALLENE),
    ('hexatriene,    4 atoms', 'CC=C=C=CC',    SU_CIS_TRANS),
    ('heptatetraene, 5 atoms', 'CC=C=C=C=CC',  SU_ALLENE),
    ('octapentaene,  6 atoms', 'CC=C=C=C=C=CC', SU_CIS_TRANS),
])
def test_the_cumulene_ladder_alternates_axial_and_planar(name, smi, kind):
    kinds = _kinds(smi)
    assert kind in kinds, f'{name}: expected kind {kind}, got {kinds}'
    other = SU_ALLENE if kind is SU_CIS_TRANS else SU_CIS_TRANS
    assert other not in kinds, f'{name}: got the other kind too'


def test_an_axial_unit_anchors_on_the_chain_centre_and_a_planar_one_on_a_terminal():
    """The anchor is not decoration -- it is where the parity is stored, and for a long chain the
    centre atom is NOT adjacent to either substituted terminal.  An implementation that stored an
    allene's parity next to the centre instead of ON it gives the same answer for a three-atom
    chain and the wrong one for a five-atom chain, which is why the five-atom row is here.
    """
    # 2,3-pentadiene: C C = C = C C, slots 0..4 -> centre is the middle chain atom
    axial = _stereogenic('CC=C=CC', SU_ALLENE)[0]
    chain_centre = axial['anchor']
    # the anchor's two refs are the substituents of the two TERMINALS, not of the centre
    assert axial['refs'][0] != chain_centre and axial['refs'][2] != chain_centre
    # and for the five-atom chain the same holds, with the refs two bonds further out
    long_axial = _stereogenic('CC=C=C=C=CC', SU_ALLENE)[0]
    assert long_axial['refs'][0] != long_axial['anchor']
    assert long_axial['refs'][2] != long_axial['anchor']


def test_a_symmetric_allene_is_not_stereogenic():
    """Allene itself, and 1,1-dimethylallene.  Both have a mirror automorphism through the axis,
    so the two configurations are one molecule.  The candidate is still perceived.
    """
    assert _stereogenic('C=C=C', SU_ALLENE) == []
    assert _stereogenic('CC(C)=C=C', SU_ALLENE) == []
    # ...but 2,3-pentadiene, the textbook chiral allene, IS
    assert len(_stereogenic('CC=C=CC', SU_ALLENE)) == 1


def test_the_small_ring_cut_applies_to_planar_units_and_not_to_axial_ones():
    """Cyclohexene's double bond is not a cis/trans unit: the ring path holds the two terminals
    on one side and there is no second configuration to reach.  Cyclooctene's is, which is
    chemically right -- (Z)- and (E)-cyclooctene are both isolable.

    The cut is deliberately NOT applied to axial units, whose terminals are perpendicular rather
    than coplanar, so a ring cannot hold them together the same way.
    """
    assert _stereogenic('C1=CCCCC1', SU_CIS_TRANS) == []          # cyclohexene
    assert _stereogenic('C1=CCCCCC1', SU_CIS_TRANS) == []         # cycloheptene
    assert len(_stereogenic('C1=CCCCCCC1', SU_CIS_TRANS)) == 1    # cyclooctene, at threshold


def test_a_cumulene_terminal_needs_a_substituent_to_be_arranged():
    """1,3-butadiene's terminal CH2 has nothing but hydrogens on one side of nothing -- a terminal
    with a single non-chain direction cannot express two arrangements.  Propadiene's terminals
    likewise.  This is the same exclusion the amide carbonyl relies on, restated for chains.
    """
    assert _stereogenic('C=CC=C', SU_CIS_TRANS) == []
    assert _stereogenic('C=C=C', SU_ALLENE) == []


def test_an_axial_configuration_distinguishes_its_two_enantiomers():
    """The weakest thing that must be true of allene support, and the thing the round-trip defect
    below does NOT break: the two tags produce two different molecules.

    Stated through InChI as well as through the arena, because InChI is an outside witness and
    the arena's parity numbering is our own.
    """
    at = [u['parity'] for u in _units('CC=[C@]=CC', SU_ALLENE)]
    atat = [u['parity'] for u in _units('CC=[C@@]=CC', SU_ALLENE)]
    assert at[0] and atat[0] and at[0] != atat[0]
    if inchi_library_loaded():
        assert molecule_to_inchi(read_smiles('CC=[C@]=CC')) \
            != molecule_to_inchi(read_smiles('CC=[C@@]=CC'))


def test_an_axial_terminals_hydrogen_holds_the_position_the_bracket_would_give_it():
    """2,3-pentadiene spelled twice.  `[CH]` writes the hydrogen where OpenSMILES puts it -- right
    after the bond to the atom before it -- and a bare `C` leaves the same position implicit, so
    the two strings are ONE molecule and one axial configuration.

    The discriminator between reader and writer: this is the pair that says which of them orders a
    terminal's directions correctly, where a round trip only says that they disagree.
    """
    assert read_smiles('C[CH]=[C@]=[CH]C') == read_smiles('CC=[C@]=CC')
    assert read_smiles('C[CH]=[C@@]=[CH]C') == read_smiles('CC=[C@@]=CC')
    assert read_smiles('CC=[C@]=CC') != read_smiles('CC=[C@@]=CC')


@pytest.mark.parametrize('smi', ['CC=[C@]=CC', 'CC=[C@@]=CC', 'CCC=[C@]=CC',
                                 'CC=[C@]=CCl', 'OC=[C@]=CO'])
def test_an_axial_configuration_survives_a_smiles_round_trip(smi):
    """One configuration must be one string, whatever creation order produced it.

    Every string here has an allene terminal whose second direction is an implicit hydrogen, which
    is the case where the reader and the writer have to agree about a written position no token
    occupies; `smi_written_pair` is where that position is inserted.
    """
    first, second = _roundtrip(smi)
    assert first == second


# ---- the remaining axial defect, recorded as a failing acceptance criterion ------------------
#
# `strict=True` on purpose: when this is fixed it must start failing as xpass, so that nobody has
# to remember to come back and delete a marker.
#
# The tag sits on a chain atom that is not the chain's CENTRE -- atom 3 of a five-atom chain, whose
# unit the arena anchors on atom 4 -- so the reader finds no unit under it and reports
# `smiles:stereo-no-unit`.  The implicit hydrogen is not the discriminator:
# `CC(F)=[C@]=C=C=C(F)C` drops it too, and `CC=C=[C@]=C=CC` stores it.

@pytest.mark.xfail(strict=True, reason='V3 defect: an axial configuration stated on a chain atom '
                                       'other than the centre of an odd cumulene is discarded at '
                                       'read -- the unit is perceived on the centre and reported '
                                       'stereogenic, but its parity stays 0.  Input fidelity says '
                                       'a stated descriptor is stored.')
def test_a_stated_axial_configuration_on_a_long_chain_is_stored():
    unit = _units('CC=[C@]=C=C=CC', SU_ALLENE)[0]
    assert unit['stereogenic'], 'precondition: the unit is stereogenic'
    assert unit['parity'] != 0, 'the stated configuration was dropped'


# ============================================================================================
# 4. SPIRO.  A spiro atom has all four bonds in rings, so a candidate pass that looks for
#    out-of-ring substituents finds nothing and the atom is simply missed.  It is also the case
#    where "are these two branches different" has to be asked per ring rather than globally.
# ============================================================================================

def test_a_spiro_atom_is_stereogenic_when_neither_ring_is_symmetric():
    """Spiro[4.4]nonane-1,6-diol.  The quaternary spiro carbon carries no substituent of its own
    -- all four of its directions run into rings -- and it is nonetheless a stereocentre, because
    each ring is desymmetrised by its hydroxyl.
    """
    mol = read_smiles('OC1CCCC12CCCC2O')
    units = mol.stereogenic_units()
    anchors = {u['anchor'] for u in units if u['kind'] == SU_TETRA}
    # the spiro atom is the one whose every neighbour is in a ring; find it by degree 4 and no H
    spiro = [u['anchor'] for u in units
             if u['kind'] == SU_TETRA and u['unnamed_mask'] == 0
             and all(r is not None for r in u['refs'])]
    assert spiro, f'the spiro carbon is not among the stereogenic units {anchors}'


def test_a_symmetric_spiro_atom_is_not_stereogenic():
    """Spiro[4.4]nonane itself.  Both rings are symmetric about the spiro atom, so a mirror is an
    automorphism and there is nothing to name.  The per-ring question, not the global one: a
    global Morgan test alone gets this wrong in one direction or the other.
    """
    assert _stereogenic('C1CCCC12CCCC2', SU_TETRA) == []


def test_a_spiro_configuration_survives_a_smiles_round_trip():
    first, second = _roundtrip('O[C@@H]1CCC[C@]12CCC[C@@H]2O')
    assert first == second
    assert first.count('@') >= 3, 'spiro configurations were dropped on the way out'


# ============================================================================================
# 5. RETRANSLATION AND FIDELITY.
#
#    A patcher that asks whether an UNMATCHED atom was stereogenic in the OLD structure and then
#    stamps the raw parity value forward, with no retranslation, is the failure mode this group pins:
#    `C[C@@H]1CCC(=O)O1` raises `KeyError: 2` through V2's reactor while `CC1CCC(=O)O1` goes through.
#
#    The requirement it implies is the one tested here: a parity is meaningful only against the
#    frame of directions it was measured in, so the API must offer RETRANSLATION as an operation
#    a caller invokes -- never a bit for a caller to copy.
# ============================================================================================

# the fixture pair: one lactone, with and without a stated configuration
_LACTONE_ACHIRAL = 'CC1CCC(=O)O1'      # gamma-valerolactone, no configuration stated
_LACTONE_CHIRAL = 'C[C@@H]1CCC(=O)O1'  # (R)-gamma-valerolactone


def test_a_structural_edit_treats_a_stated_and_an_unstated_configuration_alike():
    """The acceptance criterion for the section above.

    Both members of the fixture pair go through the same edit -- opening the lactone, which is
    what a hydrolysis template does -- and the presence of a configuration must not change
    whether the edit SUCCEEDS.  Through V2's reactor the chiral member raises `KeyError`.
    """
    results = []
    for smi in (_LACTONE_ACHIRAL, _LACTONE_CHIRAL):
        mol = read_smiles(smi)
        ring_bond = _lactone_ring_bond(mol)
        with mol.edit():
            mol.delete_bond(*ring_bond)
        results.append(mol.smiles)
    assert len(results) == 2, 'an edit raised on one member of the pair'


def test_a_configuration_the_edit_destroyed_is_dropped_and_not_stamped_forward():
    """Opening the ring leaves the former stereocentre with two hydrogens, so it is no longer
    stereogenic and its parity names nothing.  The parity must be GONE, not carried over.

    This is the silent half of the failure mode: where stamping the bit forward raises nothing, it
    carries a sign whose frame does not exist.
    """
    mol = read_smiles(_LACTONE_CHIRAL)
    assert [u for u in mol.stereogenic_units() if u['parity']], 'precondition: a stated centre'
    with mol.edit():
        mol.delete_bond(*_lactone_ring_bond(mol))
    assert mol.stereogenic_units() == [], 'the opened ring left a stereogenic unit behind'
    assert '@' not in mol.smiles, 'a dead configuration was stamped forward'


def _lactone_ring_bond(mol):
    """The ring C(sp3)-O bond of a gamma-lactone.

    Found rather than hard-coded, because the reader assigns its own slots and a literal pair
    would silently start pointing at a different bond if canonical ordering changed.  The
    discriminators: both ends in a ring, one oxygen and one carbon, and the carbon carries a
    hydrogen -- which is what distinguishes the ring C-O from the carbonyl's ester C-O.
    """
    for bond in mol.bonds():
        a, b = mol.atom(bond.n), mol.atom(bond.m)
        if {a.element, b.element} != {6, 8}:      # element is an atomic number
            continue
        (oxygen, _), (carbon, catom) = ((bond.n, a), (bond.m, b)) if a.element == 8 \
            else ((bond.m, b), (bond.n, a))
        if bond.in_ring and catom.implicit_h:
            return oxygen, carbon
    raise AssertionError('no lactone ring C-O bond found')


def test_translate_stereo_answers_in_the_callers_own_direction_order():
    """`translate_stereo` is the first-class retranslation a reactor needs.

    Handed the unit's own frame it returns the stored parity; handed one transposition of it, it
    returns the other parity.  That is the whole contract, and it is what makes copying a raw bit
    unnecessary.
    """
    mol = read_smiles('C[C@H](N)CC')
    unit = [u for u in mol.stereo_units() if u['kind'] == SU_TETRA and u['parity']][0]
    refs = unit['refs']
    own = mol.translate_stereo(unit['anchor'], refs)
    assert own == unit['parity'], 'the unit\'s own frame did not reproduce its parity'
    swapped = (refs[1], refs[0]) + refs[2:]
    assert mol.translate_stereo(unit['anchor'], swapped) != own, \
        'one transposition did not flip the parity'
    # and swapping twice returns to the original: parity is a permutation sign, not a flag
    twice = (refs[1], refs[0], refs[3], refs[2])
    assert mol.translate_stereo(unit['anchor'], twice) == own


def test_translate_stereo_refuses_a_frame_that_is_not_the_units_own():
    """The operation is total on valid permutations and REFUSES otherwise, rather than answering
    something plausible.  A reactor that hands over the wrong frame gets an exception, which is
    the failure mode the raw-bit copy did not have.
    """
    mol = read_smiles('C[C@H](N)CC')
    unit = [u for u in mol.stereo_units() if u['kind'] == SU_TETRA and u['parity']][0]
    with pytest.raises(ValueError):
        mol.translate_stereo(unit['anchor'], unit['refs'][:2] + (unit['refs'][0], None))


def test_a_molecule_and_its_enantiomer_stay_distinguishable():
    """The fidelity invariant, across every kind that has two configurations."""
    pairs = [('C[C@H](N)CC', 'C[C@@H](N)CC'),                    # tetrahedral
             ('CCCC[N@+](C)(CC)CCC', 'CCCC[N@@+](C)(CC)CCC'),    # ammonium
             ('CC[S@](C)=NC', 'CC[S@@](C)=NC'),                  # sulfilimine
             ('C[S@](=O)c1ccccc1', 'C[S@@](=O)c1ccccc1'),        # sulfoxide
             ('CC=[C@]=CC', 'CC=[C@@]=CC'),                      # axial
             ('C/C=C/CC', 'C/C=C\\CC')]                          # planar
    for left, right in pairs:
        assert read_smiles(left).smiles != read_smiles(right).smiles, \
            f'{left} and {right} wrote the same string'


# Alanine, spelled six ways: three traversals of one enantiomer and three of the other.  Which
# group each spelling belongs to was checked against InChI's `/m` layer rather than by eye, because
# hand-deriving `@` from a written order is exactly the step that is easy to get backwards -- three
# of these six were, on the first attempt at this test.
_L_ALANINE = ['N[C@@H](C)C(O)=O', 'C([C@@H](N)C)(O)=O', 'C[C@H](N)C(=O)O']       # InChI /m0
_D_ALANINE = ['[C@@H](N)(C)C(O)=O', 'OC(=O)[C@H](N)C', 'OC(=O)[C@@H](C)N']       # InChI /m1


@pytest.mark.parametrize('name,spellings', [('L', _L_ALANINE), ('D', _D_ALANINE)])
def test_an_automorphic_relabelling_does_not_change_stereo_meaning(name, spellings):
    """The other half of the fidelity invariant: the SAME molecule presented in different atom
    orders is one configuration, so it must write one string.

    This is the test that would catch a canonical form whose tie-break depends on input order.
    """
    written = {read_smiles(s).smiles for s in spellings}
    assert len(written) == 1, \
        f'{name}-alanine: one configuration wrote {len(written)} strings: {written}'


def test_the_two_alanine_spelling_groups_are_the_two_enantiomers():
    """Guards the fixture above.  If both groups ever wrote the same string the test would pass
    for the wrong reason -- vacuously, on a writer that had stopped emitting stereo at all.
    """
    left = {read_smiles(s).smiles for s in _L_ALANINE}
    right = {read_smiles(s).smiles for s in _D_ALANINE}
    assert left != right
    if inchi_library_loaded():
        assert molecule_to_inchi(read_smiles(_L_ALANINE[0])) \
            != molecule_to_inchi(read_smiles(_D_ALANINE[0]))


def test_a_configuration_stated_on_a_centre_that_is_not_stereogenic_is_stored_not_rejected():
    """"Input by default is garbage" -- a descriptor on a centre with two identical substituents
    is contradictory, and the answer is to STORE it, never to refuse the record.

    Isopropanol's central carbon carries two methyls, so no configuration is nameable there.  The
    parse must succeed, the parity must be kept on the candidate, and the stereogenicity pass --
    not the parser -- is what declines to call it a centre.
    """
    for smi in ('C[C@H](C)C', 'C[C@@H](C)C', 'C[C@H](C)O'):
        mol = read_smiles(smi)                      # must not raise
        stated = [u for u in mol.stereo_units() if u['parity']]
        assert stated, f'{smi}: the stated configuration was discarded at parse'
        assert not any(u['stereogenic'] for u in stated), \
            f'{smi}: a non-stereogenic centre was reported stereogenic'


def test_a_contradictory_configuration_never_refuses_the_record():
    """The same rule for descriptors that are structurally impossible rather than merely
    redundant: a tag on an atom with too few directions to measure one.  Parsing succeeds and
    the answer boundary, not the parser, is where a caller may complain.
    """
    for smi in ('C[C@H]C', '[C@H](C)C', 'C[C@](C)(C)C', 'O=[C@H]C'):
        read_smiles(smi)     # the assertion IS that this does not raise
