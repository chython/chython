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
"""`_hydrogens.pxi`: one implicit-hydrogen derivation, and the same answer whatever read the molecule.

WHAT THIS FILE IS GUARDING.  One derivation against one valence table, reached by every reader and by
`chython.calc_implicit` alike, and no aromatic atom is refused for being aromatic -- a benzene CH
answers 1.  The tests below pin the two halves that would otherwise drift apart:

* the derivation reaches the CHEMISTRY collection, not the SMILES notation model, so a charged
  aromatic and an element outside the organic subset both answer;
* the only atoms it refuses are the ones whose class the RING settles, and `kekule()` plus a
  fill-only second pass closes even those.

THE AMBIGUOUS CLASS IS NAMED BY THE CLASSIFIER AND NOT BY THIS FILE.  No test here lists elements: a
case is expected to answer `HYD_AMBIGUOUS_AROMATIC` because `arom_classify_atom` returns AROM_MAY for
it, and that is the whole criterion.  If a new arm of that table becomes undecidable, these tests are
what should start reporting it.
"""
from pytest import mark
from chython.core import H_UNKNOWN, MoleculeContainer, read_smiles
from chython.core._core import (HYD_AMBIGUOUS_AROMATIC, HYD_DERIVED, HYD_NO_AROMATIC_FORM,
                                HYD_NO_VALENCE_RULE, HYD_REASON_MASK, derive_implicit_hydrogen,
                                derive_implicit_hydrogens)


def _wipe(molecule):
    """Every count unknown, which is what a format with no hydrogen channel hands us."""
    for atom in molecule.atoms():
        molecule.set_hydrogens(atom.n, H_UNKNOWN)
    return molecule


def _derive_all(source):
    """`[(element, count, reason)]` in arena order, deriving each atom from scratch."""
    molecule = read_smiles(source)
    out = []
    for atom in molecule.atoms():
        count, reason = derive_implicit_hydrogen(molecule, atom.n)
        out.append((atom.element, count, reason & HYD_REASON_MASK))
    return out


# THE FIVE ATOMS THE READER'S OWN DOCSTRING NAMES.  `smi_read_h` promises "0 for thiophene S, 0 for
# furan O, 0 for pyridine N, 1 for benzene C and 0 for a fusion carbon", and the general derivation
# has to agree with it wherever a format states as much as SMILES does -- which is what
# `count_stated=True` means.
@mark.parametrize('source,index,expected', [
    ('c1ccccc1', 0, 1),             # benzene CH
    ('c1ccc2ccccc2c1', 3, 0),       # naphthalene fusion carbon, no hydrogen whatever the Kekule form
    ('c1ccncc1', 3, 0),             # pyridine N
    ('c1ccsc1', 3, 0),              # thiophene S
    ('c1ccoc1', 3, 0),              # furan O
])
def test_the_reader_s_five_atoms(source, index, expected):
    molecule = read_smiles(source)
    n = list(molecule.atoms())[index].n
    count, reason = derive_implicit_hydrogen(molecule, n, True)
    assert count == expected
    assert reason & HYD_REASON_MASK == HYD_DERIVED


def test_the_smiles_reader_still_answers_them_itself():
    """The reader's stored counts, not a re-derivation -- the refactor had to leave these alone.

    `smi_read_h` shares the CLASSIFICATION with `_hydrogens.pxi` and keeps its own valence lookup on
    `smv_default_h`, because a bare SMILES atom's count is fixed by OpenSMILES rather than by
    chemistry.  So this asserts the reader's output directly: if the shared half ever starts deciding
    the valence half too, this is the test that notices.
    """
    assert [a.implicit_h for a in read_smiles('c1ccccc1').atoms()] == [1] * 6
    assert [a.implicit_h for a in read_smiles('c1ccncc1').atoms()] == [1, 1, 1, 0, 1, 1]
    assert [a.implicit_h for a in read_smiles('c1cc[nH]c1').atoms()] == [1, 1, 1, 1, 1]
    assert [a.implicit_h for a in read_smiles('c1ccsc1').atoms()] == [1, 1, 1, 0, 1]
    assert [a.implicit_h for a in read_smiles('c1ccc2ccccc2c1').atoms()][3] == 0


# --- the two things the SMILES NOTATION model cannot do, and the reason the collection is consulted

def test_a_charged_aromatic_derives():
    """N-methylpyridinium: `smv_default_h` ignores charge entirely, so it could not answer this.

    Build the derivation on the notation model and every MDL record with a charged aromatic loses its
    counts -- which is most of the azolium and pyridinium chemistry in any corpus.
    """
    result = _derive_all('[n+]1(C)ccccc1')
    assert result[0] == (7, 0, HYD_DERIVED)
    assert [count for _, count, _ in result] == [0, 3, 1, 1, 1, 1, 1]


def test_an_element_outside_the_organic_subset_derives():
    """Tetramethylstannane and selenophene: eleven elements have a notation model, ~1036 have rows."""
    assert _derive_all('[Sn](C)(C)(C)C')[0] == (50, 0, HYD_DERIVED)
    assert _derive_all('c1cc[se]c1')[3] == (34, 0, HYD_DERIVED)


def test_a_radical_is_not_charged_twice():
    """A nitroxide oxygen: radical, one single bond, 0 hydrogens.  The regression this file caught.

    `val_implicit_h` takes `radical` as its own argument, so adding a unit of bond order for the
    unpaired electron as well counts it twice -- and the collection answers 0 at sum 1 and NO RULE at
    sum 2, so the double count turned a derivable atom into an undecidable one.  `smv_default_h` has
    no radical parameter, which is why the SMILES reader legitimately does add it to its sum.
    """
    result = _derive_all('CN([O])C |^1:2|')
    assert result[2] == (8, 0, HYD_DERIVED)


def test_a_dative_bond_is_outside_the_count():
    """`[Fe]~N(C)(C)C`: the nitrogen answers as the three-coordinate amine it is.

    Order 8 is in neither the bond-order sum nor the classifier's neighbour count, so a donated lone
    pair cannot make a metal complex look like a valence problem.  Same policy as
    `chemistry/_implicit.py`'s `environment_of`, now applied by the C walk.
    """
    result = _derive_all('[Fe]~N(C)(C)C')
    assert result[1] == (7, 0, HYD_DERIVED)
    assert result[0] == (26, 0, HYD_DERIVED)


# --- the ambiguous class: what it is, and that it is nothing wider

@mark.parametrize('source,index', [
    ('c1ccncc1', 3),                # pyridine's N, and pyrrole's presents identical arguments
    ('c1cc[nH]c1', 3),              # so does this one, once the bracket's statement is set aside
    ('c1ccc2[nH]cnc2c1', 4),        # benzimidazole, both nitrogens
    ('c1ccc2[nH]cnc2c1', 6),
])
def test_the_ring_decides_these_and_a_table_cannot(source, index):
    """`None` with `HYD_AMBIGUOUS_AROMATIC`, because the two readings are both valences.

    With one hydrogen the nitrogen donates its lone pair and takes no ring double bond; without one it
    contributes a single pi electron and must take one.  A local look sees the same atom either way,
    so the honest answer is that the count is not derivable here -- not a guessed zero, which would
    read as a fact and would take `kekule()`'s freedom away.
    """
    molecule = read_smiles(source)
    n = list(molecule.atoms())[index].n
    count, reason = derive_implicit_hydrogen(molecule, n)
    assert count is None
    assert reason & HYD_REASON_MASK == HYD_AMBIGUOUS_AROMATIC


@mark.parametrize('source', ['c1ccccc1', 'c1ccsc1', 'c1ccoc1', 'c1ccc2ccccc2c1', '[n+]1(C)ccccc1',
                             'c1cc[se]c1', 'CC(=O)O', 'OC(=O)c1ccccc1', 'Clc1ccc(F)cc1'])
def test_nothing_else_is_ambiguous(source):
    """The refusal is narrow by construction: exactly the atoms `arom_classify_atom` calls AROM_MAY.

    This is the test that would fail if the gate were ever widened into "any aromatic heteroatom",
    which is the shape the regression takes.
    """
    for element, count, reason in _derive_all(source):
        assert reason != HYD_AMBIGUOUS_AROMATIC, f'element {element} was refused'
        assert count is not None


def test_a_format_that_states_its_counts_may_withdraw_the_gate():
    """`count_stated=True` collapses AROM_MAY to the no-hydrogen reading, which is what SMILES means.

    A bare aromatic `n` is a nitrogen with no hydrogen by OpenSMILES' own rules, so the language has
    already answered and the gate would be refusing a question that is not open.  No other format may
    say this, which is why it is an argument and not the default.
    """
    molecule = read_smiles('c1ccncc1')
    n = list(molecule.atoms())[3].n
    assert derive_implicit_hydrogen(molecule, n) == (None, HYD_AMBIGUOUS_AROMATIC)
    assert derive_implicit_hydrogen(molecule, n, True) == (0, HYD_DERIVED)


# --- the whole-molecule sweep

def test_the_sweep_writes_the_counts_and_reports_the_rest():
    molecule = _wipe(read_smiles('c1ccncc1'))
    assert [a.implicit_h for a in molecule.atoms()] == [None] * 6
    undecided = molecule.derive_hydrogens()
    counts = [a.implicit_h for a in molecule.atoms()]
    assert counts == [1, 1, 1, None, 1, 1]
    nitrogen = list(molecule.atoms())[3].n
    assert undecided == {nitrogen: HYD_AMBIGUOUS_AROMATIC}


def test_the_sweep_leaves_a_stated_count_alone():
    """`stated` is the reader saying "the record gave this one"; nothing here second-guesses it.

    Ferrocene's iron is the case that matters: a reader that has derived 0 from something this pass
    cannot see must be able to keep it.
    """
    molecule = read_smiles('c1ccncc1')
    carbon = list(molecule.atoms())[0].n
    molecule.set_hydrogens(carbon, 7)
    molecule.derive_hydrogens([carbon])
    assert list(molecule.atoms())[0].implicit_h == 7


def test_fill_only_touches_nothing_that_claims_a_count():
    molecule = read_smiles('c1ccncc1')
    carbon = list(molecule.atoms())[0].n
    molecule.set_hydrogens(carbon, 7)
    molecule.derive_hydrogens(fill_only=True)
    assert list(molecule.atoms())[0].implicit_h == 7, 'fill-only overwrote a stated count'
    molecule.derive_hydrogens()
    assert list(molecule.atoms())[0].implicit_h == 1, 'the overwriting mode did not overwrite'


def test_kekule_settles_the_ambiguous_atom_BY_ITSELF():
    """THE POINT OF THE WHOLE ARRANGEMENT: read leaves it open, and `kekule()` ALONE closes it.

    After kekulisation the ring holds definite orders, there is no aromatic bond left to be ambiguous
    about, and the ordinary valence rows answer.  So the atoms the reader was right to refuse are
    exactly the atoms `kekule()` can settle -- and it settles them itself, with no second call from
    anybody.  A fill run by `canonicalize()` as a line of its own would make a hand-run
    `mol.kekule()` weaker than the same call inside the pipeline; this test is what forbids that
    asymmetry.
    """
    molecule = _wipe(read_smiles('c1ccncc1'))
    molecule.derive_hydrogens()
    assert [a.implicit_h for a in molecule.atoms()] == [1, 1, 1, None, 1, 1]
    molecule.kekule()
    assert [a.implicit_h for a in molecule.atoms()] == [1, 1, 1, 0, 1, 1]


def test_kekule_s_heal_does_not_touch_a_count_that_is_stated():
    """Fill-only, so the reader's own answer wins wherever it has one.

    Ferrocene is the case that matters and this is its shape: a count no valence row reproduces, held
    by an atom in an aromatic system, which a blanket recompute inside `kekule()` would destroy.
    """
    molecule = read_smiles('c1ccncc1')
    carbon = list(molecule.atoms())[0].n
    molecule.set_hydrogens(carbon, 7)
    molecule.kekule()
    assert list(molecule.atoms())[0].implicit_h == 7


def test_kekule_leaves_an_unresolved_system_alone():
    """"If no valence errors" is per system: five aromatic carbons get no invented hydrogens.

    `c1cccc1` has an odd atom count, so no charge and no hydrogen makes it even and the matching comes
    back deficient.  A valence row asked about a deficient atom answers with the hydrogens that fill
    the deficit -- which would invent them and hide the very defect `kekule()` reports.  So those
    atoms stay unknown and `check_valence()` still finds them.
    """
    molecule = _wipe(read_smiles('c1cccc1'))
    result = molecule.kekule()
    assert result.unresolved
    assert [a.implicit_h for a in molecule.atoms()] == [None] * 5


def test_kekule_heals_one_ring_while_another_fails():
    """A ring that failed does not cost a ring that succeeded its counts.  Two systems, one molecule."""
    molecule = _wipe(read_smiles('c1cccc1.c1ccncc1'))
    result = molecule.kekule()
    assert result.unresolved
    assert [a.implicit_h for a in molecule.atoms()] == [None] * 5 + [1, 1, 1, 0, 1, 1]


def test_kekule_s_heal_waits_for_a_caller_s_own_edit_scope():
    """Inside an outer scope the orders are still pending, so there is nothing to derive from yet.

    The heal cannot run there and must not raise there either -- a caller building a molecule
    mid-flight is a supported caller.  It derives once its own scope has closed, which is what
    `kekule()`'s docstring tells it to do.
    """
    molecule = _wipe(read_smiles('c1ccncc1'))
    with molecule.edit():
        molecule.kekule()
    assert [a.implicit_h for a in molecule.atoms()] == [None] * 6
    molecule.derive_hydrogens(fill_only=True)
    assert [a.implicit_h for a in molecule.atoms()] == [1, 1, 1, 0, 1, 1]


def test_the_sweep_survives_an_empty_molecule():
    """An empty container, not an empty string -- the SMILES reader refuses that one.

    It is not a corner nobody reaches: a reader that has just built a container and hit a record with
    no atom block calls the sweep before it knows.
    """
    assert derive_implicit_hydrogens(MoleculeContainer()) == {}


def test_an_element_with_no_aromatic_form_is_flagged_and_still_answered():
    """The flag is a bit beside the outcome, not a fifth outcome, because both can be true.

    A carbon written aromatic but carrying an exocyclic double bond -- a quinone carbonyl -- is read
    as saturated and has a perfectly good count.  Reporting "no aromatic form" instead of the count
    would throw the count away; reporting the count without the observation would hide the input's
    problem.  So the caller gets both and decides which to log.
    """
    molecule = read_smiles('[se]1cccc1')
    for atom in molecule.atoms():
        count, reason = derive_implicit_hydrogen(molecule, atom.n)
        assert count is not None
        assert reason & HYD_REASON_MASK in (HYD_DERIVED, HYD_NO_VALENCE_RULE)


def test_no_valence_rule_is_not_the_same_refusal_as_ambiguity():
    """Two reasons for `None`, and conflating them is how a coverage hole reads as bad input.

    `HYD_NO_VALENCE_RULE` says the collection describes nothing here -- our gap, and no claim about
    the molecule.  `HYD_AMBIGUOUS_AROMATIC` says the collection describes it fine and the RING has to
    pick.  The first is closed by writing a row; the second by kekulising.  A caller that cannot tell
    them apart cannot tell which of those to do.
    """
    assert HYD_NO_VALENCE_RULE != HYD_AMBIGUOUS_AROMATIC
    assert HYD_NO_AROMATIC_FORM & HYD_REASON_MASK == 0, 'the flag bit must not collide with a reason'
