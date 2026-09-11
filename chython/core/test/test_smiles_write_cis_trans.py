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
"""`/` and `\\`: the directional bonds, and the one convention behind them.

THE CONVENTION, and it is measured rather than chosen: core parity 2 (odd) means refs[0] and refs[2]
-- one named direction from each terminal, both always present by rulings F26 and F47 -- lie on the
SAME side.  The anchor is chython 2 with RDKit 2026.03.4 confirming the geometry, and it is pinned
here by `test_the_cis_trans_anchor_comes_back_out`, which asserts the character sequence and not a
derived quantity.  Like `@` it has no anchor inside the core, so it is shared -- with the reader's
`smi_cis_sign` and, since `eef2732`, with `_inchi.pxi`'s `ICH_CIS_TRANS_FLIP`.  THREE consumers of one
unanchored convention is exactly the shape in which two of them quietly drift apart, so the last
section of this file measures ours against libinchi's rather than against ourselves.

WHY THIS FILE IS SEPARATE from test_smiles_write_stereo.py: a direction is a property of a BOND, so
its failure modes are different in kind.  A tetrahedral sign can only be wrong about a permutation; a
direction can also be wrong about which of two ends it was read from, can contradict a direction
three bonds away, and can be unsatisfiable.  Those need their own fixtures.
"""
from itertools import permutations, product
from math import factorial
from random import Random

from pytest import importorskip, mark

from chython.core import MoleculeContainer
from chython.core._core import (inchi_to_molecule, inchi_library_loaded, molecule_to_inchi,
                                read_smiles, smw_stereo_seed_labels, smw_traversal, write_smiles)


# only the last section needs it, and it is a runtime fact rather than an import one: the module
# always exports the two functions and they raise when the binary is missing
needs_libinchi = mark.skipif(not inchi_library_loaded(), reason='libinchi is not loaded')


# ------------------------------------------------------------------------------------------------
# FIXTURE PLUMBING.  Same shape as the other two writer suites.
def build(atoms, bonds, order=None):
    m = MoleculeContainer()
    sids = {}
    for j in (range(len(atoms)) if order is None else order):
        element, hydrogens = atoms[j]
        sids[j] = m.add_atom(element, implicit_h=hydrogens)
    for a, b, o in bonds:
        m.add_bond(sids[a], sids[b], o)
    return m, sids


def configure(m, parities):
    """Set the given parity on each cis/trans unit, taken in ascending anchor order.

    ONLY SOUND WHERE EVERY TERMINAL HAS ONE HEAVY SUBSTITUENT, which is every fixture in this file
    except `CHLOROBUTENE`.  A stored parity means something only against the unit's refs, and refs[0]
    is the near terminal's heavy neighbours in ascending SLOT order -- so on a terminal with two of
    them the same number is two different molecules under two different creation orders, and a sweep
    using it would report a writer defect that is really a fixture defect.  Where the terminal has
    one, refs[0] is fixed by identity and the number means one molecule.  `cis_in` is the
    identity-stated form for the rest; this one stays because it keeps the symmetric fixtures short.

    Returns the anchors, so a caller can say which units it configured.
    """
    anchors = sorted(u['anchor'] for u in m.stereo_units() if u['kind'] == 1)
    assert len(anchors) == len(parities), (anchors, parities)
    for anchor, parity in zip(anchors, parities):
        m.set_parity(anchor, parity)
    return anchors


def cis_in(m, a, b, want_cis=True):
    """Store the parity that puts atoms `a` and `b` -- one per terminal -- on the given side.

    Searching over the two values and asking `translate_stereo` which one lands where we want, for
    the reason the tetrahedral suite's `configure` gives: computing the parity here would mean
    reimplementing the frame arithmetic the writer uses, so a sign error would cancel out and the
    test would pass on a broken writer.  `translate_stereo` is not the code under test.
    """
    for (anchor, partner), unit in m.chiral_bonds().items():
        refs = unit['refs']
        near, far = refs[:2], refs[2:]
        if a in near and b in far:
            frame = (a, _other(near, a), b, _other(far, b))
        elif b in near and a in far:
            frame = (b, _other(near, b), a, _other(far, a))
        else:
            continue
        for parity in (1, 2):
            m.set_parity(anchor, parity)
            if (m.translate_stereo(anchor, frame) == 2) == want_cis:
                return anchor, parity
        raise AssertionError('neither parity puts %r and %r on the wanted side' % (a, b))
    raise AssertionError('no cis/trans unit spans %r and %r' % (a, b))


def _other(pair, this):
    return pair[1] if pair[0] == this else pair[0]


def directions(smiles):
    """How many `/` and `\\` the string carries, as a count -- the tokens' positions are not the point."""
    return smiles.count('/') + smiles.count('\\')


# but-2-ene: the smallest cis/trans unit that has one, and the fixture the convention is anchored on.
BUTENE = ([(6, 3), (6, 1), (6, 1), (6, 3)], [(0, 1, 1), (1, 2, 2), (2, 3, 1)])
# 1,2-difluoroethene: the same shape with heteroatoms, which is the pair chython 2 was measured on.
DIFLUOROETHENE = ([(9, 0), (6, 1), (6, 1), (9, 0)], [(0, 1, 1), (1, 2, 2), (2, 3, 1)])
# hexa-2,4-diene: TWO units sharing one single bond, so the middle bond's single token has to satisfy
# both configurations at once.  This is the fixture the solver exists for.
DIENE = ([(6, 3), (6, 1), (6, 1), (6, 1), (6, 1), (6, 3)],
         [(0, 1, 1), (1, 2, 2), (2, 3, 1), (3, 4, 2), (4, 5, 1)])
# (Z)- and (E)-2-chlorobut-2-ene: one terminal DISUBSTITUTED, so a terminal's two directions are both
# named and the "opposite sides" relation within a terminal is under test rather than implied.
CHLOROBUTENE = ([(17, 0), (6, 0), (6, 1), (6, 3), (6, 3)],
                [(0, 1, 1), (1, 2, 2), (2, 3, 1), (1, 4, 1)])


def _ring(n, doubles, hydrogens=None):
    """A carbocycle of `n` atoms with double bonds at the given positions (i -> i+1)."""
    atoms = []
    for i in range(n):
        if hydrogens is not None:
            atoms.append((6, hydrogens))
        else:
            atoms.append((6, 1 if (i in doubles or (i - 1) % n in doubles) else 2))
    bonds = [(i, (i + 1) % n, 2 if i in doubles else 1) for i in range(n)]
    return atoms, bonds


# Cyclododecene: the double bond is IN a ring, so one of its directions lands on a ring-closure bond
# and the token has to be written at the opening.  Twelve-membered because perception refuses smaller
# ones (`_terminals_share_small_ring`), which is the boundary test_smiles_write_stereo.py's sibling
# suite pins from the other side.
CYCLODODECENE = _ring(12, {0})
# Cyclooctatetraene: four units around one ring, so the constraint graph has a CYCLE and an odd number
# of trans units makes the set unsatisfiable.  The only fixture that reaches the unwind path.
COT = _ring(8, {0, 2, 4, 6}, hydrogens=1)
# 1,3-difluoroallene: a bond kind whose configuration is NOT a direction, and is not written yet.
DIFLUOROALLENE = ([(6, 1), (6, 0), (6, 1), (9, 0), (9, 0)],
                  [(0, 1, 2), (1, 2, 2), (0, 3, 1), (2, 4, 1)])


# ------------------------------------------------------------------------------------------------
# THE ANCHOR.  Character sequences, no arithmetic in the way.
def test_the_cis_trans_anchor_comes_back_out():
    """Parity 2 is CIS, and the two strings it produces are asserted literally.

    but-2-ene's unit has refs (methyl, None, methyl', None), so the frame pair is the two methyls and
    "parity 2 means refs[0] and refs[2] on the same side" reads directly as "the methyls are cis".
    The measurement that fixes it is outside chython: `F/C=C\\F` is Z to RDKit 2026.03.4, chython 2
    stores True for that molecule, `_alkene_translate[(0, 1)]` is False so the stored bool is V2's own
    answer for its frame pair with no flip in between, and V2's True is core parity 2.

    Asserting the exact strings rather than a property is the point -- a test that recomputed "is this
    cis" from the tokens would use the writer's own rule and agree with an inverted one.
    """
    atoms, bonds = BUTENE
    m, sids = build(atoms, bonds)
    configure(m, (2,))
    assert write_smiles(m) == 'C(/C)=C/C'
    m, sids = build(atoms, bonds)
    configure(m, (1,))
    assert write_smiles(m) == 'C(/C)=C\\C'


def test_the_anchor_is_the_molecule_rdkit_says_it_is():
    """The other half of the anchor: the two strings above are Z and E to an oracle outside chython."""
    chem = importorskip('rdkit.Chem')
    for fixture in (BUTENE, DIFLUOROETHENE):
        atoms, bonds = fixture
        for parity, expected in ((2, 'STEREOZ'), (1, 'STEREOE')):
            m, sids = build(atoms, bonds)
            configure(m, (parity,))
            written = write_smiles(m)
            mol = chem.MolFromSmiles(written)
            assert mol is not None, written
            stereo = [str(b.GetStereo()) for b in mol.GetBonds()
                      if str(b.GetStereo()) != 'STEREONONE']
            assert stereo == [expected], (written, parity, stereo)


def test_the_two_parities_are_two_molecules():
    """Trivially necessary and worth pinning: an inverted convention would still pass the sweeps."""
    chem = importorskip('rdkit.Chem')
    seen = set()
    atoms, bonds = DIFLUOROETHENE
    for parity in (1, 2):
        m, sids = build(atoms, bonds)
        configure(m, (parity,))
        seen.add(chem.MolToSmiles(chem.MolFromSmiles(write_smiles(m))))
    assert seen == {'F/C=C/F', 'F/C=C\\F'}


# ------------------------------------------------------------------------------------------------
# ANTI-DRIFT.  The tokens against `translate_stereo`, with neither side rebuilding the other's answer.
def _frame_pair_agrees(m, unit, partner, tokens):
    """Whether the string puts refs[0] and refs[2] on the same side, read out of `tokens`.

    `tokens` holds BOTH halves of every directional bond, so `up(terminal -> substituent)` is a
    lookup and not an inference -- which matters, because the character in the STRING means opposite
    things depending on which end of the bond was written first.
    """
    near = tokens[(unit['anchor'], unit['refs'][0])]
    far = tokens[(partner, unit['refs'][2])]
    return near == far


@mark.parametrize('name,fixture,count', [('butene', BUTENE, 1), ('difluoroethene', DIFLUOROETHENE, 1),
                                         ('chlorobutene', CHLOROBUTENE, 1), ('diene', DIENE, 2)])
def test_the_side_in_the_string_is_translate_stereo_of_the_stored_parity(name, fixture, count):
    """For every creation order and every configuration: the string agrees with the translator.

    The two sides are genuinely independent.  `translate_stereo` reads the stored parity byte in the
    refs frame; `tokens` is what the solver decided, three relations deep in a breadth-first search
    that never looks at a parity except through `smw_dir_propagate`'s one call.  A sign error in
    either shows up here, and a sign error in BOTH would have to agree about the seed as well.
    """
    atoms, bonds = fixture
    for order in permutations(range(len(atoms))):
        for parities in product((1, 2), repeat=count):
            m, sids = build(atoms, bonds, order=list(order))
            configure(m, parities)
            probe = smw_traversal(m)
            assert not probe['lost'], (name, order, parities)
            bonds_map = m.chiral_bonds()
            assert len(bonds_map) == count, (name, order, bonds_map)
            for (anchor, partner), unit in bonds_map.items():
                same = _frame_pair_agrees(m, unit, partner, probe['tokens'])
                parity = m.translate_stereo(anchor, unit['refs'])
                assert same == (parity == 2), (name, order, parities, anchor, parity, same)


def test_the_shared_single_bond_of_a_diene_carries_exactly_one_token():
    """Two configurations, one bond, one character -- which is why the assignment has to be solved.

    Written naively, each unit would claim the middle bond and the second claim would overwrite the
    first, silently changing the configuration of the first double bond.  Three single bonds carry a
    token and `tokens` holds both halves of each, so six entries -- and the middle bond appearing
    ONCE is the whole point: two units, three characters, not four.
    """
    atoms, bonds = DIENE
    m, sids = build(atoms, bonds)
    configure(m, (2, 1))
    probe = smw_traversal(m)
    assert len(probe['tokens']) == 6, sorted(probe['tokens'])
    middle = (sids[2], sids[3])
    assert middle in probe['tokens'] and middle[::-1] in probe['tokens']
    assert probe['tokens'][middle] != probe['tokens'][middle[::-1]], \
        'a bond seen from its two ends is the opposite character, always'
    assert directions(write_smiles(m)) == 3, write_smiles(m)


def test_a_disubstituted_terminal_puts_its_two_directions_on_opposite_sides():
    """The second of the three relations, on the terminal that has two named directions.

    2-chlorobut-2-ene's C2 carries both a chlorine and a methyl, so the string must give those two
    bonds opposite characters when read from C2 -- they are the terminal's two in-plane positions and
    there is nowhere else for them to go.
    """
    atoms, bonds = CHLOROBUTENE
    m, sids = build(atoms, bonds)
    configure(m, (2,))
    tokens = smw_traversal(m)['tokens']
    assert tokens[(sids[1], sids[0])] != tokens[(sids[1], sids[4])], sorted(tokens.items())


# ------------------------------------------------------------------------------------------------
# RING CLOSURES.
def test_a_ring_closure_bond_carries_its_token_at_the_opening_only():
    """Cyclododecene: the token sits before the opening digit, and the closing digit is bare.

    Writing it at both ends needs the two characters to be OPPOSITE, which readers disagree about;
    one end is unambiguous everywhere.  Asserted through RDKit rather than by counting characters,
    because "does this string mean what we meant" is the only question that matters here.
    """
    chem = importorskip('rdkit.Chem')
    atoms, bonds = CYCLODODECENE
    for parity, expected in ((2, 'STEREOZ'), (1, 'STEREOE')):
        m, sids = build(atoms, bonds)
        configure(m, (parity,))
        written = write_smiles(m)
        assert directions(written) == 2, written
        mol = chem.MolFromSmiles(written)
        assert mol is not None, written
        stereo = [str(b.GetStereo()) for b in mol.GetBonds() if str(b.GetStereo()) != 'STEREONONE']
        assert stereo == [expected], (written, parity, stereo)


# ------------------------------------------------------------------------------------------------
# THE LOSS REPORT.
def test_contradictory_configurations_are_dropped_together_and_reported():
    """Cyclooctatetraene, where the constraint graph is a CYCLE and parity can make it odd.

    Around the ring each unit contributes one same/opposite relation and each of the four shared
    single bonds contributes one more (a bond is opposite to itself reversed), so the total is
    `xor(units) ^ 0` and the set is unsatisfiable exactly when an ODD number of the four units is
    trans.  That is a prediction of the model, not a description of the code, and both halves are
    asserted: three cis plus one trans loses ALL FOUR, and four cis writes all four.

    All four rather than the one that could not be satisfied: writing three of a contradictory set
    would hand back a string that reads as a molecule nobody stated, which is worse than a string
    that carries no configuration at all and says so.
    """
    atoms, bonds = COT
    m, sids = build(atoms, bonds)
    anchors = configure(m, (2, 2, 2, 1))
    written = write_smiles(m)
    assert directions(written) == 0, written
    assert sorted(smw_traversal(m)['lost']) == sorted(anchors)

    m, sids = build(atoms, bonds)
    configure(m, (2, 2, 2, 2))
    assert smw_traversal(m)['lost'] == ()
    assert directions(write_smiles(m)) == 4, write_smiles(m)


def test_an_even_number_of_trans_units_around_the_ring_is_satisfiable():
    """The other half of the prediction, so the test above cannot pass by refusing everything."""
    chem = importorskip('rdkit.Chem')
    atoms, bonds = COT
    for parities in ((1, 1, 1, 1), (1, 1, 2, 2), (2, 1, 1, 2)):
        m, sids = build(atoms, bonds)
        configure(m, parities)
        written = write_smiles(m)
        assert smw_traversal(m)['lost'] == (), (parities, written)
        mol = chem.MolFromSmiles(written)
        assert mol is not None, written
        assert sum(1 for b in mol.GetBonds() if str(b.GetStereo()) != 'STEREONONE') == 4, written


def test_an_atropisomer_is_reported_lost_every_time():
    """SMILES HAS NO SYNTAX FOR AN ATROPISOMER, so this one never stops being reported.

    Which makes it the entry in `lost` that documents what the report is FOR: not a bug to be fixed
    but a statement that this format cannot carry this configuration.  The molecule is still written,
    because refusing would leave the caller unable to see the structure they hold.
    """
    m = MoleculeContainer()
    rings = []
    for _ in range(2):
        a = [m.add_atom(6) for _ in range(6)]
        for i in range(6):
            m.add_bond(a[i], a[(i + 1) % 6], 2 if i % 2 == 0 else 1)
        rings.append(a)
    m.add_bond(rings[0][0], rings[1][0], 1)
    for a in rings:
        m.add_bond(a[1], m.add_atom(6, implicit_h=3), 1)
        m.add_bond(a[5], m.add_atom(6, implicit_h=3), 1)
        for i in (2, 3, 4):
            m.set_hydrogens(a[i], 1)
        # The pivot and the two ortho carbons carry no hydrogen -- three heavy neighbours each -- and
        # they have to SAY so.  `add_atom` stores H_UNKNOWN for an unsaid count and axis detection
        # refuses an anchor whose count is unknown, so leaving these three silent leaves the molecule
        # with no axis and this test with nothing to report as lost.
        for i in (0, 1, 5):
            m.set_hydrogens(a[i], 0)
    axes = [u['anchor'] for u in m.stereo_units() if u['kind'] == 3]
    assert len(axes) == 1, [u['kind'] for u in m.stereo_units()]
    m.set_parity(axes[0], 2)
    assert smw_traversal(m)['lost'] == (axes[0],)
    assert directions(write_smiles(m)) == 0


# ------------------------------------------------------------------------------------------------
# THE AXIS.  `@` on the centre of an allene, over the two TERMINALS' directions -- OpenSMILES calls it
# extended tetrahedral, the arena calls it SU_ALLENE, and it is a bond kind with an atom's syntax.
#
# THE ANCHOR IS NOT A MEASUREMENT OF THE AXIAL CASE, and saying so is the point of this comment.  No
# tool in reach can state an axial configuration in a frame the arena also states: RDKit 2026.03.4
# drops every allene tag at sanitization and does not perceive one from 3D, OpenBabel 3.1.0 refuses
# it on read, Indigo reads and writes it but its InChI export drops the layer, and the in-tree InChI
# bridge raises before libinchi is reached (`_inchi.pxi`'s SU_ALLENE branch hands `translate_stereo`
# the chain ATOMS, which are not the unit's refs).  What IS measured is the tetrahedral rule -- from a
# hand-built conformer, RDKit writes `@` exactly when the signed volume over the written order is
# negative -- and the writer applies that same rule to the axis.  `SMW_ALLENE_AT_FOR_ODD` in
# `_smiles_write.pxi` is the single flip point if a cross-format consumer ever disagrees.
#
# So these tests pin the two things that CAN regress: the round trip (two enantiomers are two strings,
# and Indigo reads ours back as the configuration we wrote) and the frame (one configuration is one
# string over every creation order).  Both hold under either value of the DEF.
def axial_in(m, a, b, odd=True):
    """Store the parity that is ODD -- or EVEN -- in the frame `(a, a', b, b')`, one atom per terminal.

    Identity-stated for the reason `test_the_axial_sweep_would_notice_a_frame_error` measures: the
    same stored NUMBER is two different molecules under two different creation orders, because refs
    are in slot order.  Searching the two values through `translate_stereo` rather than computing one
    keeps the fixture from reimplementing the writer's arithmetic, exactly as `cis_in` does.
    """
    for u in m.stereo_units():
        if u['kind'] != 2:
            continue
        refs = u['refs']
        near, far = refs[:2], refs[2:]
        if a in near and b in far:
            frame = (a, _other(near, a), b, _other(far, b))
        elif b in near and a in far:
            frame = (b, _other(near, b), a, _other(far, a))
        else:
            continue
        for parity in (1, 2):
            m.set_parity(u['anchor'], parity)
            if (m.translate_stereo(u['anchor'], frame) == 2) == odd:
                return u['anchor'], parity
        raise AssertionError('neither parity is %s for %r' % ('odd' if odd else 'even', frame))
    raise AssertionError('no axial unit spans %r and %r' % (a, b))


# 1-bromo-1-fluoro-3-chloro-3-iodoallene: four DISTINCT heavy substituents, so every direction is
# named, the four are told apart by element in any string, and no automorphism can hide a frame error.
TETRAHALOALLENE = ([(6, 0), (6, 0), (6, 0), (35, 0), (9, 0), (17, 0), (53, 0)],
                   [(0, 1, 2), (1, 2, 2), (0, 3, 1), (0, 4, 1), (2, 5, 1), (2, 6, 1)])
# A five-carbon cumulene: axially chiral, named by the arena, and UNWRITABLE (see the refusal test).
CUMULENE5 = ([(6, 1), (6, 0), (6, 0), (6, 0), (6, 1), (9, 0), (9, 0)],
             [(0, 1, 2), (1, 2, 2), (2, 3, 2), (3, 4, 2), (0, 5, 1), (4, 6, 1)])


def _frame_tag(written):
    """The tag `written` carries, re-expressed in the frame `(Br, F, Cl, I)`.

    Only for TETRAHALOALLENE, whose four directions are one halogen each, so the written order can be
    read off the string by scanning for element symbols -- no SMILES parser needed.  The point is to
    compare OUR string against ANOTHER writer's string for the same molecule without either one's atom
    order mattering: two strings agree iff their tags agree after this reduction.
    """
    ref = ['Br', 'F', 'Cl', 'I']
    seen = []
    i = 0
    while i < len(written):
        if written[i:i + 2] in ('Br', 'Cl'):
            seen.append(written[i:i + 2])
            i += 2
            continue
        if written[i] in ('F', 'I'):
            seen.append(written[i])
        i += 1
    assert sorted(seen) == sorted(ref), (written, seen)
    tag = '@@' if '@@' in written else ('@' if '@' in written else None)
    if tag is None:
        return None
    perm = [ref.index(x) for x in seen]
    swaps = sum(1 for a in range(4) for b in range(a + 1, 4) if perm[a] > perm[b])
    return ('@@' if tag == '@' else '@') if swaps % 2 else tag


def test_the_axial_anchor_comes_back_out():
    """The two enantiomers, as literal strings, with the arithmetic that makes them readable.

    Odd in the frame `(Br, F, Cl, I)` gives `C(I)(=[C@]=C(F)Br)Cl`.  Read the string's own order of
    the four directions -- `I`, `F`, `Br`, `Cl`, by appearance -- and against `(Br, F, Cl, I)` that is
    the permutation `[3, 1, 0, 2]`, four inversions, EVEN.  So the string's order carries the same
    parity as the stated frame, odd, and odd is `@` (`SMW_ALLENE_AT_FOR_ODD`).  The literal assert is
    the test; the arithmetic is here so that a reader can check the literal rather than trust it.
    """
    atoms, bonds = TETRAHALOALLENE
    m, sids = build(atoms, bonds)
    axial_in(m, sids[3], sids[5], odd=True)
    assert write_smiles(m) == 'C(I)(=[C@]=C(F)Br)Cl'
    m, sids = build(atoms, bonds)
    axial_in(m, sids[3], sids[5], odd=False)
    assert write_smiles(m) == 'C(I)(=[C@@]=C(F)Br)Cl'


def test_an_implicit_hydrogen_on_a_terminal_sits_where_it_is_written():
    """1,3-difluoroallene: each terminal is one F and one implicit H, and the H's POSITION matters.

    `FC=[C@@]=CF` for the configuration stated odd in `(F, H, F', H')`.  The near terminal's parent is
    the F, so its order is `(F, H)`; the far terminal's parent is the centre, so the H comes first and
    its order is `(H, F')`.  That is one within-pair swap away from the stated frame -- odd -- so the
    odd configuration writes `@@` here and `@` in `test_the_axial_anchor_comes_back_out`.  A writer
    that put the hydrogen last unconditionally would invert exactly this fixture and no other.
    """
    atoms, bonds = DIFLUOROALLENE
    m, sids = build(atoms, bonds)
    axial_in(m, sids[3], sids[4], odd=True)
    assert write_smiles(m) == 'FC=[C@@]=CF'
    m, sids = build(atoms, bonds)
    axial_in(m, sids[3], sids[4], odd=False)
    assert write_smiles(m) == 'FC=[C@]=CF'


@mark.parametrize('spec,fixture,ends,total', [('tetrahalo', TETRAHALOALLENE, (3, 5), 5040),
                                              ('difluoro', DIFLUOROALLENE, (3, 4), 120)])
def test_one_axial_configuration_is_one_string_over_every_creation_order(spec, fixture, ends, total):
    """The whole factorial, both enantiomers, identity-stated.  Two configurations, two strings."""
    atoms, bonds = fixture
    seen = {}
    count = 0
    for odd in (True, False):
        out = set()
        for order in permutations(range(len(atoms))):
            m, sids = build(atoms, bonds, order)
            axial_in(m, sids[ends[0]], sids[ends[1]], odd=odd)
            out.add(write_smiles(m))
            count += 1
        assert len(out) == 1, (spec, odd, sorted(out)[:4])
        seen[odd] = out.pop()
    assert count == 2 * total
    assert seen[True] != seen[False], seen
    assert seen[True].replace('@@', '@') == seen[False].replace('@@', '@'), seen


def test_the_axial_sweep_would_notice_a_frame_error():
    """Ruling F102: the sweep above is evidence only if it CAN fail, so here is it failing.

    The same stored NUMBER -- parity 2, not "odd in a named frame" -- over TETRAHALOALLENE's 5,040
    creation orders gives BOTH strings, 2,520 each, because refs are in slot order and the slot order
    is the creation order.  That is precisely the defect ruling F26 forbids reaching the output, so a
    writer that emitted the byte would fail the sweep by producing two strings, and the sweep's single
    string is a measurement rather than a tautology.
    """
    atoms, bonds = TETRAHALOALLENE
    counts = {}
    for order in permutations(range(7)):
        m, sids = build(atoms, bonds, order)
        anchor = [u['anchor'] for u in m.stereo_units() if u['kind'] == 2][0]
        m.set_parity(anchor, 2)
        written = write_smiles(m)
        counts[written] = counts.get(written, 0) + 1
    assert sorted(counts.values()) == [2520, 2520], counts
    assert set(counts) == {'C(I)(=[C@]=C(F)Br)Cl', 'C(I)(=[C@@]=C(F)Br)Cl'}


def test_the_axial_sign_agrees_with_translate_stereo_on_the_writers_own_order():
    """ANTI-DRIFT.  The tag in the string, against the parity in the order the writer says it used.

    `directions[anchor]` is `smw_allene_order`'s tuple, which for an axis is the two TERMINALS'
    directions grouped by terminal and not the anchor's own neighbours.  Putting it through
    `translate_stereo` and comparing to the character is the same trade the tetrahedral suite makes:
    neither side rebuilds the other's answer, so a sign error cannot cancel out.

    The parity-to-character mapping is LEARNED from the first case rather than written down, so that
    flipping `SMW_ALLENE_AT_FOR_ODD` breaks the two anchor tests above and nothing else.  What is under
    test here is that the mapping is the SAME for every creation order: a writer that reported an order
    it had not used would agree with itself on some orders and disagree on others.
    """
    atoms, bonds = TETRAHALOALLENE
    seen = {}
    for order in permutations(range(7)):
        if order[0] % 3:              # three of the seven first slots; the sweep above is the full one
            continue
        for odd in (True, False):
            m, sids = build(atoms, bonds, order)
            anchor, _ = axial_in(m, sids[3], sids[5], odd=odd)
            probe = smw_traversal(m)
            assert probe['lost'] == ()
            frame = probe['directions'][anchor]
            written = write_smiles(m)
            parity = m.translate_stereo(anchor, frame)
            tag = '@@' if '[C@@]' in written else ('@' if '[C@]' in written else None)
            assert tag is not None, (order, odd, written)
            assert seen.setdefault(parity, tag) == tag, (order, odd, parity, frame, written, seen)
    assert sorted(seen) == [1, 2] and len(set(seen.values())) == 2, seen


def test_indigo_reads_our_axial_string_as_the_configuration_we_wrote():
    """AN EXTERNAL ROUND TRIP, which is the strongest statement available about the axial sign.

    Indigo is the only reader in reach that keeps an allene tag at all (RDKit 2026.03.4 drops it at
    sanitization, OpenBabel 3.1.0 refuses it).  It re-writes our string in its OWN atom order, so the
    comparison is made after reducing both to the frame `(Br, F, Cl, I)` -- which also tests the
    writer's claim that the string's interleaving of the two terminals is an even permutation of the
    grouped tuple, because Indigo's interleaving is a different one and the tags still agree.

    Which tag means which enantiomer is not asserted here: our polarity lives in the anchor tests, and
    a shared convention would only be provable against a tool that states an axial configuration in a
    frame of its own -- there is none (see this section's header).  What IS proved is that a reader
    which does understand the syntax gets back what we put in, and tells the two apart.
    """
    indigo = importorskip('indigo')
    session = indigo.Indigo()
    atoms, bonds = TETRAHALOALLENE
    tags = {}
    for odd in (True, False):
        m, sids = build(atoms, bonds)
        axial_in(m, sids[3], sids[5], odd=odd)
        written = write_smiles(m)
        tags[odd] = _frame_tag(written)
        assert tags[odd] is not None, written
        assert _frame_tag(session.loadMolecule(written).smiles()) == tags[odd], written
    assert tags[True] != tags[False], tags


def test_a_longer_odd_cumulene_is_reported_lost():
    """SMILES HAS NO SYNTAX FOR A FIVE-CARBON AXIS, so this one is reported rather than written.

    The arena names it the same way it names an allene -- odd chain, anchor at the centre -- and the
    configuration is real.  But measured 2026-09-02, Indigo REFUSES `FC=C=[C@@]=C=CF` with "chirality
    on atom 3 makes no sense", and a string a reader rejects is worse than a string that says less.
    chython 2 writes the sign here; this is a deliberate divergence with a measurement behind it.
    """
    atoms, bonds = CUMULENE5
    m, sids = build(atoms, bonds)
    axial = [u['anchor'] for u in m.stereo_units() if u['kind'] == 2]
    assert len(axial) == 1
    m.set_parity(axial[0], 2)
    assert smw_traversal(m)['lost'] == (axial[0],)
    written = write_smiles(m)
    assert '@' not in written, written
    assert written == 'FC=C=C=C=CF'


def test_an_axial_configuration_is_lost_without_the_bond_tokens():
    """Under `!b` there are no `=` tokens, so there is no axis for a sign to be read against.

    `@` on a two-coordinate carbon is not a weaker statement, it is a WRONG one -- a reader takes it
    for a tetrahedral centre.  So the sign goes and the unit is reported, which is also what makes
    `!b` honest: tetrahedral signs still come out in that mode, so the caller has to be told which
    part of the stereo the string kept.
    """
    atoms, bonds = TETRAHALOALLENE
    m, sids = build(atoms, bonds)
    anchor, _ = axial_in(m, sids[3], sids[5], odd=True)
    written = write_smiles(m, '!b')
    assert '@' not in written, written
    assert smw_traversal(m, '!b')['lost'] == (anchor,)
    assert '@' in write_smiles(m)          # and it is the spec that dropped it, not the molecule


def test_an_axial_parity_that_was_never_stated_is_not_a_loss():
    """No parity, no sign, no report -- `lost` is about configurations the writer DROPPED.

    The ORDER is still reported, and that is the deliberate asymmetry: `directions` says which four
    directions an `@` on this atom would be read against, which is a fact about the axis and not about
    the parity.  So an absent entry means the writer REFUSED the axis (the cumulene above) while an
    entry with no `@` in the string means the molecule never said which enantiomer it was -- two
    different situations that would be indistinguishable if the order were withheld too.
    """
    atoms, bonds = TETRAHALOALLENE
    m, sids = build(atoms, bonds)
    assert [u['kind'] for u in m.stereo_units()] == [2]
    anchor = m.stereo_units()[0]['anchor']
    probe = smw_traversal(m)
    assert probe['lost'] == ()
    assert sorted(probe['directions'][anchor]) == sorted(m.stereo_units()[0]['refs'])
    assert '@' not in write_smiles(m)


def test_a_parity_that_was_never_stated_is_not_a_loss():
    """`lost` is about configurations the writer DROPPED, so an unset parity must not appear in it."""
    atoms, bonds = DIENE
    m, sids = build(atoms, bonds)
    probe = smw_traversal(m)
    assert probe['lost'] == () and probe['tokens'] == {}
    assert directions(write_smiles(m)) == 0


# ------------------------------------------------------------------------------------------------
# THE SWEEPS.  Ruling F26 again, for the bond kind this time.
def _stated(*parities):
    """A setup that states the configuration as raw parities -- symmetric terminals only."""
    return lambda m, sids: configure(m, parities)


def _chlorobutene_z(m, sids):
    """Chlorine cis to the far methyl, stated by IDENTITY because this fixture needs it to be."""
    cis_in(m, sids[0], sids[3], True)


CIS_TRANS_SWEEPS = [('butene', BUTENE, _stated(2), 3),
                    ('difluoroethene', DIFLUOROETHENE, _stated(1), 3),
                    ('chlorobutene', CHLOROBUTENE, _chlorobutene_z, 16),
                    ('diene-ZZ', DIENE, _stated(2, 2), 5),
                    # (2Z,4E): the entry the stereo seed exists for.  It was the strict xfail below
                    # until the seed landed, and it stays in the sweep because a regression in the
                    # seed shows up here as two strings and nowhere else.
                    ('diene-ZE', DIENE, _stated(2, 1), 6)]


def _sweep(fixture, setup, spec):
    atoms, bonds = fixture
    seen = set()
    orders = list(permutations(range(len(atoms))))
    for order in orders:
        m, sids = build(atoms, bonds, order=list(order))
        setup(m, sids)
        seen.add(write_smiles(m, spec))
    return seen, len(orders)


@mark.parametrize('name,fixture,setup,stored_count', CIS_TRANS_SWEEPS)
def test_canonical_direction_output_does_not_depend_on_the_creation_order(name, fixture, setup,
                                                                          stored_count):
    """One configuration, every creation order, ONE string -- including the seed choice.

    The seed is the extra thing this sweep tests over the tetrahedral one.  `F/C=C/F` and `F\\C=C\\F`
    are the same molecule, so the solver has a free bit per constraint component, and a canonical
    writer has to spend it the same way every time.  It does, because the seed is the first
    directional half-edge in emission order and emission order is a function of `canonical_order()`.
    """
    seen, count = _sweep(fixture, setup, '')
    assert len(seen) == 1, (name, count, sorted(seen)[:4])
    assert count == factorial(len(fixture[0]))


@mark.parametrize('name,fixture,setup,stored_count', CIS_TRANS_SWEEPS)
def test_stored_order_directions_are_not_creation_order_invariant(name, fixture, setup,
                                                                  stored_count):
    """The could-have-failed evidence (ruling F102), as a NUMBER rather than as `> 1`."""
    seen, _ = _sweep(fixture, setup, 'i')
    assert len(seen) == stored_count, (name, sorted(seen))


@mark.parametrize('name,fixture,setup,stored_count', CIS_TRANS_SWEEPS)
def test_every_stored_order_spelling_is_the_same_molecule_to_rdkit(name, fixture, setup,
                                                                   stored_count):
    """The oracle half: the many spellings are ONE molecule, so the frame is right and not just stable.

    For a direction this is stronger than it is for a tetrahedral sign, because a spelling can put
    the two ends of a double bond in either order and can reach a bond from either side.  Internal
    agreement cannot tell a consistent reversal from a correct one; this can.
    """
    chem = importorskip('rdkit.Chem')
    seen, _ = _sweep(fixture, setup, 'i')
    assert len(seen) == stored_count, name
    canonical = {chem.MolToSmiles(chem.MolFromSmiles(s)) for s in seen}
    assert len(canonical) == 1, (name, sorted(canonical))


def test_a_symmetric_diene_with_an_asymmetric_configuration_gives_one_string():
    """WAS A DEFECT, and the reason the stereo seed exists: two strings over 720 orders, 660 to 60.

    Found by the sweep above, which is what a sweep is for, and it stood as a strict xfail against
    task 6 until the seed landed.  hexa-2,4-diene's constitution is symmetric end to end, so
    `canonical_order()` had an automorphism to break and broke it from the graph alone -- while the
    two ends are constitutionally identical and stereochemically not, one Z and one E, so whichever
    end the order happened to put first decided the string.  Never a lost configuration: RDKit read
    both spellings as `C/C=C\\C=C\\C`, which is what made it dangerous, because nothing downstream
    could see it except by comparing two strings that should have been equal.

    Not the direction solver's defect either: give both units the same configuration and every one of
    the 720 orders already agreed, which is the neighbouring sweep entry.  `smw_stereo_seed` feeds the
    parities into the canonical order, so the automorphism is now broken by the thing that actually
    distinguishes the two halves.
    """
    seen, count = _sweep(DIENE, _stated(2, 1), '')
    assert count == 720
    assert len(seen) == 1, sorted(seen)


def test_the_asymmetric_diene_spelling_is_the_configured_molecule_to_rdkit():
    """The oracle half: the one string means the molecule that was configured.

    Invariance without this would be satisfied by a writer that dropped both signs, or emitted them
    consistently reversed.  Kept as its own test after the seed landed because it is a different
    claim: the sweep says "one string", this says "the right one".
    """
    chem = importorskip('rdkit.Chem')
    seen, _ = _sweep(DIENE, _stated(2, 1), '')
    assert len(seen) == 1
    assert chem.MolToSmiles(chem.MolFromSmiles(seen.pop())) == 'C/C=C\\C=C\\C'


def test_the_four_diene_configurations_are_three_compounds():
    """Four parity combinations, THREE strings: (2Z,4E) and (2E,4Z) are one compound, numbered
    from the other end.

    The collapse is as much a requirement as the separation, and it is the test a seed that
    over-separated would fail: a seed keyed on anything that distinguishes "the Z end came first"
    from "the E end came first" would hand these two four strings, and each one would be invariant
    over creation orders, so the sweeps above would all pass.  chython 2 agrees on the count -- its
    own spellings differ from ours, which is what an oracle is for -- and RDKit is asserted below to
    partition the four the same way.
    """
    chem = importorskip('rdkit.Chem')
    ours = {}
    for combo in ((2, 2), (1, 1), (2, 1), (1, 2)):
        seen, _ = _sweep(DIENE, _stated(*combo), '')
        assert len(seen) == 1, (combo, sorted(seen))
        ours[combo] = seen.pop()
    assert len(set(ours.values())) == 3, ours
    assert ours[(2, 1)] == ours[(1, 2)], ours
    # And the same partition to RDKit, which is the check that our three are THEIR three and not
    # three of ours that happen to be two of theirs plus a mistake.
    theirs = {combo: chem.MolToSmiles(chem.MolFromSmiles(s)) for combo, s in ours.items()}
    assert len(set(theirs.values())) == 3, theirs
    assert theirs[(2, 1)] == theirs[(1, 2)], theirs


# ------------------------------------------------------------------------------------------------
# THE SEED ITSELF, read directly rather than inferred from two strings.
def test_the_seed_splits_a_class_the_constitution_leaves_tied():
    """(2Z,4E): the stereo-blind refinement ties the two ends, the seed does not.

    The whole mechanism in one assertion pair.  `atoms_order` is the constitutional refinement and it
    puts C1 with C6, C2 with C5 and C3 with C4 -- three classes for six atoms, which is the
    automorphism the extremal search then had to break on slot order.  The seed gives all six atoms
    different labels, so there is no tie left to break.
    """
    m, sids = build(*DIENE)
    _stated(2, 1)(m, sids)
    classes = [m.atoms_order[sids[j]] for j in range(6)]
    assert classes[:3] == classes[5:2:-1], classes      # the palindrome IS the automorphism
    assert len(set(classes)) == 3, classes
    labels = smw_stereo_seed_labels(m)
    assert len({labels[sids[j]] for j in range(6)}) == 6, labels


def test_the_seed_leaves_a_real_symmetry_alone():
    """(2Z,4Z): the two ends ARE interchangeable, and the seed's labels stay palindromic.

    The other direction, and the one a seed built out of anything slot-shaped would fail: this
    molecule's mirror is an automorphism of the CONFIGURED molecule too, so a seed that separated its
    ends would be reporting an asymmetry the molecule does not have -- and would then pin the
    canonical order to the creation order it read the ends in, which is the defect wearing a
    different hat (ruling F95's sigma-equivariance requirement, stated in `_stereo.pxi`).
    """
    m, sids = build(*DIENE)
    _stated(2, 2)(m, sids)
    labels = [smw_stereo_seed_labels(m)[sids[j]] for j in range(6)]
    assert labels[:3] == labels[5:2:-1], labels
    assert len(set(labels)) == 3, labels


def test_there_is_no_seed_without_a_configured_parity():
    """No configured parity, no seed -- and then the order is the one it was before the seed existed.

    The early-out is a behavioural promise and not an optimisation: every stereo-free molecule in the
    suite would otherwise have its canonical order recomputed from a seed, and if that seed ever
    disagreed with the unseeded refinement the change would land on molecules that have nothing to do
    with stereo at all.
    """
    m, _ = build(*DIENE)
    assert smw_stereo_seed_labels(m) is None


def test_the_seed_is_not_taken_without_the_stereo_key():
    """`!s` stays a function of the constitution: the two configurations write ONE string.

    That is what makes `format(mol, '!s')` usable as a constitution key, and it is not automatic --
    seeding unconditionally would cost nothing in invariance and would quietly give two molecules
    with one constitution two different `!s` strings.
    """
    seen = set()
    for combo in ((2, 2), (2, 1), (1, 1)):
        got, _ = _sweep(DIENE, _stated(*combo), '!s')
        seen |= got
    assert len(seen) == 1, sorted(seen)


def test_a_random_creation_order_sample_of_a_ring_fixture_gives_one_string():
    """Cyclododecene has 12! creation orders, so this one is sampled and says so.

    Sampled with a FIXED seed: a flaky invariance test is worse than none, because the failure is
    reported against whichever order the clock happened to pick.
    """
    atoms, bonds = CYCLODODECENE
    rng = Random(20260902)
    seen = set()
    for _ in range(200):
        order = list(range(len(atoms)))
        rng.shuffle(order)
        m, sids = build(atoms, bonds, order=order)
        configure(m, (2,))
        seen.add(write_smiles(m))
    assert len(seen) == 1, sorted(seen)[:4]


# ------------------------------------------------------------------------------------------------
# WHAT SUPPRESSES A DIRECTION.
def test_no_direction_under_the_no_stereo_key():
    atoms, bonds = BUTENE
    m, sids = build(atoms, bonds)
    configure(m, (2,))
    assert write_smiles(m, '!s') == 'C(C)=CC'
    probe = smw_traversal(m, '!s')
    assert probe['tokens'] == {} and probe['lost'] == ()


def test_no_direction_when_bond_tokens_are_suppressed():
    """`!b` drops the `=` too, so there is nothing for a direction to be written on -- AND IT IS LOST.

    `!b` is a caller saying they do not want bond tokens, so the dropped direction is not a writer
    defect.  It is still a dropped CONFIGURATION, and `lost` is the writer's list of those, so it is
    reported: under `!b` a tetrahedral sign still comes out, so "the string carries stereo" stays true
    while ceasing to be the whole truth, and a caller comparing two strings for identity has to know
    which part went.  Reporting is also the only way the axial refusal under `!b` is visible at all.

    A unit with no stated parity is NOT reported here -- nothing was dropped.  So the list tracks
    configurations, not units, in this mode exactly as in every other.
    """
    atoms, bonds = BUTENE
    m, sids = build(atoms, bonds)
    anchors = configure(m, (2,))
    assert directions(write_smiles(m, '!b')) == 0
    probe = smw_traversal(m, '!b')
    assert probe['tokens'] == {} and probe['lost'] == tuple(anchors)

    m, sids = build(atoms, bonds)          # same molecule, configuration never stated
    probe = smw_traversal(m, '!b')
    assert probe['tokens'] == {} and probe['lost'] == ()


# ------------------------------------------------------------------------------------------------
# THE THIRD CONSUMER.  Our polarity against libinchi's, not against ourselves.
#
# The arena epic's `eef2732` gave the core's cis/trans parity a stated geometric meaning through
# `ICH_CIS_TRANS_FLIP`, adopting InChI's rule because the signed-volume argument that pins the
# tetrahedral and allene conventions is degenerate for a planar bond.  That constant and this file's
# convention are now two independent statements about the same parity byte, arrived at from opposite
# directions -- ours measured through chython 2 with RDKit confirming the geometry, theirs adopted
# from libinchi -- and the ONLY way to find out whether they agree is to run the parity out through
# one and read it back through the other.  They agree.  This is the test that says so, and the test
# that fails if either side is ever flipped alone.
@needs_libinchi
@mark.parametrize('text,layer,label', [('F/C=C\\F',   '-', 'Z-1,2-difluoroethene'),
                                       ('F/C=C/F',    '+', 'E-1,2-difluoroethene'),
                                       ('C/C=C\\C',   '-', 'Z-2-butene'),
                                       ('C/C=C/C',    '+', 'E-2-butene'),
                                       ('Cl/C=C\\Br', '-', 'Z-1-bromo-2-chloroethene'),
                                       ('Cl/C=C/Br',  '+', 'E-1-bromo-2-chloroethene')])
def test_our_cis_trans_polarity_is_the_one_libinchi_means(text, layer, label):
    """`/b...-` is Z and `/b...+` is E, and our string round-trips to the same InChI.

    ABSOLUTE, which is the point: `-` versus `+` in InChI's `/b` layer is documented and external, so
    this is not a self-consistency check.  `F/C=C\\F` is Z as a fact about the world and not as a fact
    about our conventions, which is why the parametrisation names the isomer.

    BE PRECISE ABOUT WHAT THIS PINS, because two chains are easy to confuse: this one is
    reader-then-export, so flipping the reader's `smi_cis_sign` or `ICH_CIS_TRANS_FLIP` alone breaks
    it and flipping BOTH would not.  The WRITER's polarity is not in this chain at all -- it is pinned
    by `test_the_cis_trans_anchor_comes_back_out` and by the literal string in the next test.  The
    three-way agreement is the conjunction of those, not a claim either one makes alone.
    """
    m = read_smiles(text)
    inchi = molecule_to_inchi(m)
    b = [p for p in inchi.split('/') if p.startswith('b') and p[1:2].isdigit()]
    assert len(b) == 1, (label, inchi)
    assert b[0].endswith(layer), '%s: expected /b...%s, got %r' % (label, layer, b[0])

    # and the geometry survives the trip out and back, so the export and the import agree with each
    # other as well as with us -- a pair of matching flips would show up here and nowhere else
    assert molecule_to_inchi(inchi_to_molecule(inchi)) == inchi, (label, inchi)


@needs_libinchi
def test_the_writers_own_respelling_does_not_move_the_geometry():
    """`F/C=C\\F` and `C(/F)=C/F` are the same molecule, and libinchi is asked rather than told.

    The writer chooses its own start atom, so its output for a cis double bond is frequently spelled
    from the other end than the input was -- which is the case where "the character in the string
    means opposite things depending on which end came first" turns into a real sign error.  Two
    spellings, one InChI.
    """
    a, b = read_smiles('F/C=C\\F'), read_smiles('C(/F)=C/F')
    assert molecule_to_inchi(a) == molecule_to_inchi(b)
    assert write_smiles(a) == write_smiles(b) == 'C(/F)=C/F'
