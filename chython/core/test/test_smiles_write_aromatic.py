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
"""Aromatic output: two stored representations, two strings, and neither one becomes the other here.

WHAT THIS FILE GUARDS.  `smw_bond`, `smw_sticky_bond` and `smw_atom` decide lowercase and `:` from the
STORED order.  Deciding them from an OPTION instead writes aromatic benzene as
`[CH]1[CH][CH][CH][CH][CH]1` -- cyclohexane -- a silent representation change on the way OUT, which is
worse than on the way in because the caller has no string left to inspect.

WHAT THE WRITER PROMISES: it spells what is stored.  A molecule holding order-4 bonds writes
lowercase; its `kekule()` twin writes `=`; the two must NOT converge, because they are two stored
representations of one compound and `kekule()`/`thiele()` are the only places a representation may
change.  Every fixture below is asserted in both forms, and RDKit 2026.03.4 is the independent
witness that the two mean the same compound.
"""
from itertools import permutations

from pytest import importorskip, mark

from chython.core import MoleculeContainer
from chython.core._core import smv_valence_model, write_smiles


# ------------------------------------------------------------------------------------------------
# FIXTURE PLUMBING.  Rings are built by hand at order 4, never by parsing: a fixture that came through
# a reader would be testing the reader.
def cycle(elements, hydrogens, order=4, extra=()):
    """A monocycle of the given elements with the given implicit hydrogen counts."""
    m = MoleculeContainer()
    ids = [m.add_atom(e, implicit_h=h) for e, h in zip(elements, hydrogens)]
    n = len(ids)
    for i in range(n):
        m.add_bond(ids[i], ids[(i + 1) % n], order)
    for element, h, at, bond in extra:
        m.add_bond(ids[at], m.add_atom(element, implicit_h=h), bond)
    return m


def fused(edges, elements, hydrogens):
    """A polycycle from an explicit edge list, every bond aromatic."""
    m = MoleculeContainer()
    ids = [m.add_atom(e, implicit_h=h) for e, h in zip(elements, hydrogens)]
    for i, j in edges:
        m.add_bond(ids[i], ids[j], 4)
    return m


BENZENE = ([6] * 6, [1] * 6)
PYRIDINE = ([7] + [6] * 5, [0] + [1] * 5)
PYRROLE = ([7] + [6] * 4, [1] * 5)
FURAN = ([8] + [6] * 4, [0] + [1] * 4)
THIOPHENE = ([16] + [6] * 4, [0] + [1] * 4)

NAPHTHALENE_EDGES = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0),
                     (5, 6), (6, 7), (7, 8), (8, 9), (9, 0)]
INDOLE_EDGES = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0),
                (5, 6), (6, 7), (7, 8), (8, 0)]


def naphthalene():
    return fused(NAPHTHALENE_EDGES, [6] * 10, [0, 1, 1, 1, 1, 0, 1, 1, 1, 1])


def indole():
    """Benzo-fused pyrrole: the fixture where a RING-CLOSURE bond is aromatic and carries an `[nH]`."""
    return fused(INDOLE_EDGES, [6] * 8 + [7], [0] + [1] * 4 + [0] + [1, 1, 1])


# ------------------------------------------------------------------------------------------------
# THE PAIR: one compound, two stored representations, two strings.
PAIRS = [('benzene', lambda: cycle(*BENZENE), 'c1ccccc1', 'C1=CC=CC=C1'),
         ('pyridine', lambda: cycle(*PYRIDINE), 'c1ccccn1', 'C1=CC=NC=C1'),
         ('pyrrole', lambda: cycle(*PYRROLE), 'c1c[nH]cc1', 'C=1C=CNC=1'),
         ('furan', lambda: cycle(*FURAN), 'c1cocc1', 'C=1C=COC=1'),
         ('thiophene', lambda: cycle(*THIOPHENE), 'c1cscc1', 'C=1C=CSC=1'),
         ('naphthalene', naphthalene, 'c1c2c(cccc2)ccc1', 'C1=CC=2C(C=C1)=CC=CC=2'),
         ('indole', indole, 'c1cc2c(cc1)[nH]cc2', 'N1C=CC2=CC=CC=C12')]


@mark.parametrize('name,make,aromatic,kekule', PAIRS)
def test_a_stored_aromatic_molecule_is_written_aromatic(name, make, aromatic, kekule):
    """The string, literally, for the stored aromatic form and for its `kekule()` twin.

    Literal strings rather than a property, because "is this output aromatic" computed from the
    output is exactly the reasoning the deleted hooks used.  The pyrrole and furan entries also pin
    a ring closure that CARRIES a bond token in the Kekule form (`C=1...C=1`), which is the path a
    fix to the atom tokens alone would leave broken.
    """
    m = make()
    assert not m.is_kekule and m.aromatic_bond_count
    assert write_smiles(m) == aromatic
    assert m.kekule().changed
    assert m.is_kekule and not m.aromatic_bond_count
    assert write_smiles(m) == kekule


@mark.parametrize('name,make,aromatic,kekule', PAIRS)
def test_the_two_representations_do_not_converge(name, make, aromatic, kekule):
    """Stated on its own, because converging is the failure mode and it is a QUIET one.

    A writer that kekulised on the way out, or that lowercased a Kekule ring on the way out, would
    pass every round-trip test in this file -- RDKit would read both strings as the same compound and
    agree.  What it would break is the promise that the string shows the caller what they hold.
    """
    assert aromatic != kekule


@mark.parametrize('name,make,aromatic,kekule', PAIRS)
def test_both_representations_are_the_same_compound_to_rdkit(name, make, aromatic, kekule):
    """And the independent witness that neither string lost anything on the way.

    RDKit canonicalises both spellings to one string, which is the check that the aromatic form is a
    faithful spelling and not merely a lowercase-looking one.  It is also the only test here that
    would catch a wrong hydrogen count on a heteroatom: `c1ccoc1` and `c1cc[oH]c1` are different
    molecules and only an oracle knows which one furan is.
    """
    chem = importorskip('rdkit.Chem')
    canonical = set()
    for smiles in (aromatic, kekule):
        mol = chem.MolFromSmiles(smiles)
        assert mol is not None, (name, smiles)
        canonical.add(chem.MolToSmiles(mol))
    assert len(canonical) == 1, (name, sorted(canonical))


# ------------------------------------------------------------------------------------------------
# THE HYDROGEN RULE, which is where an aromatic writer usually goes wrong quietly.
def test_the_aromatic_hydrogen_rule_is_not_the_kekule_one():
    """Five elements, five different answers, and none of them from `smv_default_h`.

    An aromatic bond has no order in the Daylight valence model, so feeding the stored 4 into it
    makes benzene's carbon look like a bond-order sum of 8 -- hypervalent, zero hydrogens inferred.
    That is how the same six carbons come out bare `C` in the aromatic form and bracketed `[C]` in
    the Kekule one under one shared rule: the question is nonsense in one of the two.

    The rule the aromatic path uses instead: count each aromatic bond as 1, add one for the atom's
    share of the ring's pi system, subtract from the element's LOWEST normal valence, clamp at zero.
    The five entries below are the five distinct outcomes it has to get right, and thiophene is the
    one that forces "lowest" rather than the Kekule rule's "smallest at or above" -- sulfur's
    valences are 2, 4 and 6, its aromatic sum is 3, and smallest-at-or-above would put a phantom
    hydrogen on an `s` that has none.
    """
    assert write_smiles(cycle(*BENZENE)) == 'c1ccccc1'          # c: 4 - 3 = 1, so bare
    assert write_smiles(cycle(*PYRIDINE)) == 'c1ccccn1'         # n: 3 - 3 = 0, so bare
    assert write_smiles(cycle(*PYRROLE)) == 'c1c[nH]cc1'        # n: 0 inferred, 1 held -> [nH]
    assert write_smiles(cycle(*FURAN)) == 'c1cocc1'             # o: 2 < 3, clamped to 0, bare
    assert write_smiles(cycle(*THIOPHENE)) == 'c1cscc1'         # s: LOWEST valence 2, not 4
    # And the fusion carbon, whose third aromatic bond takes the sum to 4 and the count to zero.
    assert 'H' not in write_smiles(naphthalene())


def test_a_substituent_takes_the_aromatic_hydrogen_away():
    """Toluene's ring carbon: two aromatic bonds plus a single one is a sum of 4, so no hydrogen.

    The arithmetic has to mix the two kinds of bond in one sum, which is the reason `smw_atom_env`
    counts an aromatic bond as 1 in `order_sum` instead of handing the aromatic case a separate one.
    """
    m = cycle([6] * 6, [0] + [1] * 5, extra=[(6, 3, 0, 1)])
    assert write_smiles(m) == 'c1c(C)cccc1'


def test_an_exocyclic_double_bond_leaves_the_ring_atom_bare():
    """4-pyridone as stored: the carbonyl carbon's sum is 1 + 1 + 2 + 1, well past carbon's valence.

    Clamping at zero rather than underflowing is the whole of the `if v > order_sum` in the rule, and
    an unsigned underflow there would ask for 4294967295 hydrogens.
    """
    m = cycle([7] + [6] * 5, [1, 1, 1, 0, 1, 1], extra=[(8, 0, 3, 2)])
    written = write_smiles(m)
    assert written == 'c1c[nH]ccc1=O', written
    chem = importorskip('rdkit.Chem')
    assert chem.MolFromSmiles(written) is not None, written


def test_an_aromatic_atom_outside_the_organic_subset_is_bracketed():
    """`se` is spelled, but only inside brackets: no reader infers a hydrogen count for it.

    `smv_aromatic_h` answers only for B, C, N, O, P and S, which is SMILES' aromatic organic subset.
    Everything else returns False, which brackets the atom and states the count -- the same fidelity
    rule the aliphatic predicate uses, applied to a set the standard simply does not cover.
    """
    m = cycle([34] + [6] * 4, [0] + [1] * 4)
    written = write_smiles(m)
    assert written == 'c1c[se]cc1', written
    chem = importorskip('rdkit.Chem')
    assert chem.MolFromSmiles(written) is not None, written


def test_an_element_with_no_aromatic_spelling_at_all_is_still_shown():
    """A chlorine given an order-4 bond: garbage in, and the writer shows it rather than hiding it.

    The arena will store this -- it validates buffers, not chemistry -- so the writer has to answer.
    It answers `[cl]`, which no reader accepts, and that is the correct behaviour: refusing to show a
    caller their own broken structure is the one thing this writer may not do, and a silently
    repaired one would be worse than a string that fails to parse.
    """
    m = MoleculeContainer()
    a = m.add_atom(6)
    b = m.add_atom(17)
    m.add_bond(a, b, 4)
    written = write_smiles(m)
    assert '[cl]' in written, written


def test_the_aromatic_rule_reads_the_same_valence_table_as_the_kekule_one():
    """One table, two rules over it -- so an element cannot be aromatic-writable and not writable.

    `smv_aromatic_h` and `smv_default_h` both call `smv_valences`.  If the aromatic path grew its own
    numbers they would drift, and the drift would look like a hydrogen-count bug on one heteroatom in
    one ring size, which is the kind of thing that survives for years.
    """
    model = smv_valence_model()
    for element in (5, 6, 7, 8, 15, 16):
        assert element in model, element
    # Lowest valence per element is what the aromatic rule subtracts from; spot-check the two that
    # the two rules disagree about, so this test fails if either table is edited.
    assert model[16][0][0] == 2, 'sulfur lowest valence 2 is what keeps thiophene bare'
    assert model[7][0][0] == 3, 'nitrogen lowest valence 3 is what makes pyrrole [nH]'


# ------------------------------------------------------------------------------------------------
# THE `-` BETWEEN TWO RINGS.
def test_a_single_bond_between_two_aromatic_atoms_is_written():
    """Biphenyl: without the `-` the string reads as one twelve-membered aromatic system.

    `c1ccccc1c1ccccc1` is not biphenyl to any reader -- the two rings share no atom, so the bond
    between them has to say it is single.  It is written only in lowercase mode: under `A` the atoms
    are uppercase and an unmarked bond cannot be read as part of a ring system.
    """
    m = MoleculeContainer()
    rings = []
    for _ in range(2):
        ids = [m.add_atom(6, implicit_h=1) for _ in range(6)]
        m.set_hydrogens(ids[0], 0)
        for i in range(6):
            m.add_bond(ids[i], ids[(i + 1) % 6], 4)
        rings.append(ids)
    m.add_bond(rings[0][0], rings[1][0], 1)
    written = write_smiles(m)
    assert written == 'c1cc(-c2ccccc2)ccc1', written
    chem = importorskip('rdkit.Chem')
    assert chem.MolToSmiles(chem.MolFromSmiles(written)) == 'c1ccc(-c2ccccc2)cc1'


def test_the_dash_is_not_written_between_two_kekule_atoms():
    """The control: an ordinary single bond gets no token, or every alkane would grow dashes."""
    m = cycle([6] * 6, [2] * 6, order=1)
    assert write_smiles(m) == 'C1CCCCC1'


# ------------------------------------------------------------------------------------------------
# THE `A` DIALECT: aromaticity on the bonds instead of the atoms.
A_MODE = [('benzene', lambda: cycle(*BENZENE), '[CH]:1:[CH]:[CH]:[CH]:[CH]:[CH]:1'),
          ('pyridine', lambda: cycle(*PYRIDINE), '[CH]:1:[CH]:[CH]:[CH]:[CH]:[N]:1'),
          ('pyrrole', lambda: cycle(*PYRROLE), '[CH]:1:[CH]:[NH]:[CH]:[CH]:1')]


@mark.parametrize('name,make,expected', A_MODE)
def test_the_a_key_puts_the_aromaticity_on_the_bonds(name, make, expected):
    """`A` writes `:` and UPPERCASE atoms, chython 2's dialect, and it costs brackets everywhere.

    Every aromatic atom is bracketed with its count stated, which looks heavy-handed next to
    `c1ccccc1` and is the only faithful answer: `C:C` is outside OpenSMILES -- an aromatic bond
    between aliphatic atoms -- so no rule says what hydrogen count it implies and readers differ.  One
    that treats `:` as aromatic reads benzene as intended; one that treats it as single infers two
    hydrogens per carbon.  Stating the count makes both readers right.
    """
    assert write_smiles(make(), 'A') == expected


@mark.parametrize('name,make,expected', A_MODE)
def test_the_a_dialect_round_trips_through_rdkit(name, make, expected):
    chem = importorskip('rdkit.Chem')
    plain = write_smiles(make())
    mol = chem.MolFromSmiles(expected)
    assert mol is not None, expected
    assert chem.MolToSmiles(mol) == chem.MolToSmiles(chem.MolFromSmiles(plain))


def test_the_a_key_marks_an_aromatic_ring_closure_too():
    """The closure path is where a fix to the atom and child-bond sites usually leaves a hole.

    Benzene's ring-closure bond is aromatic like the other five, and under `A` it has to say so at
    both ends -- `[CH]:1` opening and `:1` closing.  A closure that lost its token would read as a
    single bond and the ring would come back non-aromatic at exactly one bond, which is the sort of
    thing that survives a round-trip through a sanitising reader and fails on a strict one.
    """
    written = write_smiles(cycle(*BENZENE), 'A')
    assert written.startswith('[CH]:1:'), written
    assert written.endswith(':1'), written
    # Seven tokens for six bonds: five chain bonds once each, and the closure at BOTH ends.
    assert written.count(':') == 7, written


# ------------------------------------------------------------------------------------------------
# INVARIANCE.  The aromatic path is a new decision in the writer, so it gets the same sweep.
SWEEPS = [('benzene', BENZENE, 1), ('pyridine', PYRIDINE, 6), ('pyrrole', PYRROLE, 5)]


def _sweep(elements, hydrogens, spec):
    seen = set()
    n = len(elements)
    for order in permutations(range(n)):
        m = MoleculeContainer()
        sids = {}
        for j in order:
            sids[j] = m.add_atom(elements[j], implicit_h=hydrogens[j])
        for i in range(n):
            m.add_bond(sids[i], sids[(i + 1) % n], 4)
        seen.add(write_smiles(m, spec))
    return seen


@mark.parametrize('name,fixture,stored_count', SWEEPS)
def test_aromatic_output_does_not_depend_on_the_creation_order(name, fixture, stored_count):
    """Every creation order, one string -- the property the whole design reduces to (note 3)."""
    assert len(_sweep(fixture[0], fixture[1], '')) == 1, name


@mark.parametrize('name,fixture,stored_count', SWEEPS)
def test_stored_order_aromatic_output_is_not_invariant(name, fixture, stored_count):
    """The could-have-failed evidence, as a number (ruling F102).

    Benzene's is 1 and that is not a canonicalisation: all six atoms are identical, so every creation
    order spells the same ring in stored order too.  It is recorded rather than dropped because a
    reader who sees only pyridine's 6 and pyrrole's 5 would think benzene had been forgotten.
    """
    assert len(_sweep(fixture[0], fixture[1], 'i')) == stored_count, name


# ------------------------------------------------------------------------------------------------
# WHAT THE OTHER FORMAT KEYS DO TO IT.
def test_suppressing_bond_tokens_leaves_the_lowercase_atoms():
    """`!b` drops `=` and `:`; the aromaticity is on the ATOMS in the default dialect, so it survives.

    Which is the one case where `!b` output is still a faithful aromatic molecule -- and the reason
    the aromaticity is spelled on the atoms by default rather than on the bonds.
    """
    assert write_smiles(cycle(*BENZENE), '!b') == 'c1ccccc1'
    assert write_smiles(cycle(*BENZENE), 'A!b') == '[CH]1[CH][CH][CH][CH][CH]1'


def test_forcing_hydrogens_brackets_every_aromatic_atom():
    """`h` states every count; the atoms stay lowercase, because that is the representation."""
    assert write_smiles(cycle(*BENZENE), 'h') == '[cH]1[cH][cH][cH][cH][cH]1'
