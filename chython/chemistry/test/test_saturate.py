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
"""`saturate`: bond-order recovery on skeletons flattened to single bonds.

Every case is built from SMILES, kekulised, every double and triple bond flattened to single, and
handed back for the orders to be derived again.  Charges and hydrogen counts stay as written -- the
count is the only signal that can force a multiple bond on a skeleton with no coordinates."""
from ast import (AnnAssign, Assign, Attribute, Call, For, Import, ImportFrom, Name, parse, Set,
                 SetComp, walk)
from inspect import signature
from pathlib import Path

import pytest

from .. import _saturate
from .._saturate import saturate
from ...core import H_UNKNOWN, INFO, LOST, MoleculeContainer, read_smiles, REFUSED


#: The module under test as SOURCE: two claims below are about the code, not its output -- no RNG,
#: and no loop over an unordered container.  Neither is provable from the outside, since a seeded RNG
#: is deterministic run to run and a set of small integer tuples does not move under a hash seed.
SOURCE = Path(__file__).resolve().parent.parent / '_saturate.py'


#: The corpus, by class.  Public compounds only.
CORPUS = {
    'ketone': 'CC(C)=O',
    'aldehyde': 'CCC=O',
    'carboxylic acid': 'CC(=O)O',
    'carboxylate': 'CC(=O)[O-]',
    'ester': 'CCOC(C)=O',
    'amide': 'CC(=O)NC',
    'urea': 'NC(=O)N',
    'thiourea': 'NC(=S)N',
    'nitrile': 'CC#N',
    'alkyne': 'CC#CC',
    'alkene': 'CC=CC',
    'allene': 'C=C=C',
    'carbon dioxide': 'O=C=O',
    'nitro': 'C[N+](=O)[O-]',
    'nitroso': 'CN=O',
    'oxime': 'CC(C)=NO',
    'imine': 'CC=NC',
    'azide': 'CN=[N+]=[N-]',
    'diazonium': 'c1ccccc1[N+]#N',
    'isocyanate': 'CN=C=O',
    'sulfone': 'CS(C)(=O)=O',
    'sulfonamide': 'CS(=O)(=O)N',
    'sulfoxide': 'CS(C)=O',
    'sulfonate': 'CS(=O)(=O)[O-]',
    'sulfonyl chloride': 'CS(=O)(=O)Cl',
    'thioamide': 'CC(=S)N',
    'phosphate': 'OP(=O)(O)O',
    'phosphate ester': 'CCOP(=O)(O)O',
    'phosphine oxide': 'CP(C)(C)=O',
    'ammonium': 'C[NH3+]',
    'guanidinium': 'NC(N)=[NH2+]',
    'benzene': 'c1ccccc1',
    'naphthalene': 'c1ccc2ccccc2c1',
    'pyridine': 'c1ccncc1',
    'pyrrole': 'c1cc[nH]c1',
    'imidazole': 'c1c[nH]cn1',
    'thiophene': 'c1ccsc1',
    'furan': 'c1ccoc1',
    'pyrimidine': 'c1cncnc1',
    'indole': 'c1ccc2[nH]ccc2c1',
    'tetrazole': 'c1nnn[nH]1',
    'phenol': 'Oc1ccccc1',
    'aniline': 'Nc1ccccc1',
    'benzoic acid': 'OC(=O)c1ccccc1',
    'acetophenone': 'CC(=O)c1ccccc1',
    'quinone': 'C1=CC(=O)C=CC1=O',
    'nicotinamide': 'NC(=O)c1cccnc1',
    'caffeine': 'CN1C=NC2=C1C(=O)N(C)C(=O)N2C',
    'aspirin': 'CC(=O)Oc1ccccc1C(=O)O',
    'paracetamol': 'CC(=O)Nc1ccc(O)cc1',
    'ibuprofen': 'CC(C)Cc1ccc(cc1)C(C)C(=O)O',
    'warfarin': 'CC(=O)CC(c1ccccc1)c1c(O)c2ccccc2oc1=O',
    'benzamidine': 'NC(=[NH2+])c1ccccc1',
    'glucose': 'OCC1OC(O)C(O)C(O)C1O',
    'biotin': 'OC(=O)CCCCC1SCC2NC(=O)NC12',
    'penicillin skeleton': 'CC1(C)SC2C(NC(=O)Cc3ccccc3)C(=O)N2C1C(=O)O',
    'adenosine phosphate': 'Nc1ncnc2n(cnc12)C1OC(COP(=O)(O)OP(=O)(O)O)C(O)C1O',
    'sildenafil': 'CCCc1nn(C)c2c1nc(nc2=O)c1cc(ccc1OCC)S(=O)(=O)N1CCN(C)CC1',
    'imatinib': 'Cc1ccc(cc1Nc1nccc(n1)c1cccnc1)C(=O)Nc1ccc(CN2CCN(C)CC2)cc1',
    'atropine skeleton': 'CN1C2CCC1CC(C2)OC(=O)C(CO)c1ccccc1',
    'porphine': 'C1=CC2=NC1=CC1=CC=C(N1)C=C1C=CC(=N1)C=C1C=CC(=C2)N1',
    'acetate salt': 'CC(=O)[O-].[Na+]',
    'iron amine complex': 'CN(C)C~[Fe]',
    'ferrous ion': '[Fe+2]',
}

#: The 22 of the corpus whose orders come back as a different Kekule form of the same molecule.  Named
#: rather than counted, so a case moving in or out of the set is a visible diff.  Every one contains an
#: aromatic ring and nothing else in the corpus does.
KEKULE_ALTERNATIVE = frozenset((
    'diazonium', 'benzene', 'naphthalene', 'pyridine', 'pyrimidine', 'indole', 'phenol', 'aniline',
    'benzoic acid', 'acetophenone', 'nicotinamide', 'aspirin', 'paracetamol', 'ibuprofen', 'warfarin',
    'benzamidine', 'penicillin skeleton', 'adenosine phosphate', 'sildenafil', 'imatinib',
    'atropine skeleton', 'porphine'))

#: The two entries an orders-comparison cannot score, because there was nothing in them to derive:
#: `[Fe+2]` has no bonds and the iron complex's only non-single bond is the dative one, which `flatten`
#: deliberately does not touch.  They pass because the pass declined to act, which is a different
#: outcome from recovery; `test_a_dative_bond_and_its_metal_are_left_alone` asserts it directly.
NOTHING_TO_DERIVE = frozenset(('iron amine complex', 'ferrous ion'))


def flatten(molecule: MoleculeContainer) -> MoleculeContainer:
    """Every double and triple bond down to single, in place: the input contract, manufactured.

    Aromatic bonds cannot appear (the caller kekulises first) and a dative bond is left alone: order 8
    is outside valence bookkeeping here, so flattening one would invent a covalent bond.
    """
    bonds = [(bond.n, bond.m) for bond in molecule.bonds() if bond.order in (2, 3)]
    with molecule.edit():
        for n, m in bonds:
            molecule.set_order(n, m, 1)
    return molecule


def prepared(smiles: str):
    """`(flattened molecule, the molecule it should come back as)`, both kekulised."""
    answer = read_smiles(smiles)
    answer.kekule()
    molecule = read_smiles(smiles)
    molecule.kekule()
    return flatten(molecule), answer


def orders(molecule: MoleculeContainer):
    return {(min(b.n, b.m), max(b.n, b.m)): b.order for b in molecule.bonds()}


def aromatic_form(molecule: MoleculeContainer) -> str:
    """The canonical SMILES of the aromatised molecule -- the answer as a molecule, not a spelling."""
    molecule.thiele()
    return molecule.smiles


# --- the gate: recovery on molecules whose answer is already known -------------------------------- #

@pytest.mark.parametrize('name', sorted(CORPUS))
def test_recovers_the_molecule(name):
    """Recovery as a molecule: 64 of 64, compared after aromatisation, every one satisfied."""
    molecule, answer = prepared(CORPUS[name])
    assert saturate(molecule) is True
    assert aromatic_form(molecule) == aromatic_form(answer)


@pytest.mark.parametrize('name', sorted(set(CORPUS) - KEKULE_ALTERNATIVE))
def test_recovers_the_exact_bond_orders(name):
    """The stricter score: the 42 whose raw orders come back identical, bond for bond.

    Two of the 42 are in `NOTHING_TO_DERIVE` and score nothing -- they are here for the non-mutation
    half of the claim, not the recovery half.
    """
    molecule, answer = prepared(CORPUS[name])
    saturate(molecule)
    assert orders(molecule) == orders(answer)


def test_a_dative_bond_and_its_metal_are_left_alone():
    """The distinction `NOTHING_TO_DERIVE` names: correctly untouched, not recovered.

    Order 8 exempts an atom and triggers nothing, contributing no order to either end's valence, so
    every atom settles at zero headroom and no fragment is opened.  What must not happen is the iron
    acquiring a bond order because a search noticed spare capacity on a loosely described metal.
    """
    molecule, answer = prepared('CN(C)C~[Fe]')
    log = molecule.log
    assert saturate(molecule) is True
    assert not log, [record.rule for record in log]
    assert orders(molecule) == orders(answer) == {(1, 2): 1, (2, 3): 1, (2, 4): 1, (4, 5): 8}
    assert molecule.charge_of(5) == 0 and molecule.element_of(5) == 26


@pytest.mark.parametrize('name', sorted(KEKULE_ALTERNATIVE))
def test_a_kekule_alternative_is_reported_and_not_hidden(name):
    """The other 22.  A different Kekule form is allowed -- claiming it is the only one is not."""
    molecule, answer = prepared(CORPUS[name])
    log = molecule.log
    saturate(molecule)
    assert orders(molecule) != orders(answer), (
        f'{name} now recovers its exact orders, so it belongs outside KEKULE_ALTERNATIVE and the '
        'two recovery numbers in this file need restating')
    assert any(record.rule == 'saturate:ambiguous' for record in log), (
        f'{name} came back as a different Kekule form and nothing in the log said the choice was '
        'not unique.  That is the failure mode this pass exists to avoid')


# --- determinism ---------------------------------------------------------------------------------- #

def test_the_answer_does_not_change_between_runs():
    """The same input saturated repeatedly is byte-identical, log included.

    The whole corpus is walked rather than one molecule sampled: a single molecule with one forced
    answer would pass while a shuffled search over an ambiguous ring still varied.
    """
    def run():
        out = []
        for name in sorted(CORPUS):
            molecule, _ = prepared(CORPUS[name])
            log = molecule.log
            result = saturate(molecule)
            out.append((name, result, tuple(sorted(orders(molecule).items())),
                        tuple((r.rule, r.atoms, r.message) for r in log)))
        return tuple(out)

    first = run()
    for attempt in range(8):
        assert run() == first, f'run {attempt + 2} disagreed with run 1'


def test_the_module_contains_no_randomness():
    """No RNG: not seeded, not optional, absent.  Asserted over the source, because that is the claim.

    A behavioural test cannot prove it: a seeded RNG passes the test above and still makes the answer
    depend on where the seed came from.  It parses rather than greps because it must catch a use and
    not a mention -- names come off the syntax tree, so `random` in prose stays free.
    """
    tree = parse(SOURCE.read_text(encoding='utf-8'))
    imported = set()
    for node in walk(tree):
        if isinstance(node, Import):
            imported.update(alias.name.split('.')[0] for alias in node.names)
        elif isinstance(node, ImportFrom):
            imported.add((node.module or '').split('.')[0])
            imported.update(alias.name for alias in node.names)
        elif isinstance(node, Name):
            imported.add(node.id)
        elif isinstance(node, Attribute):
            imported.add(node.attr)
    for banned in ('random', 'shuffle', 'sample', 'choice', 'randint', 'seed', 'getrandbits'):
        assert banned not in imported, (
            f'{banned!r} is used in saturate.py.  chython 2 shuffled twice and one molecule came '
            'back differently on different runs; there is no seeded version of that which is '
            'acceptable in a perception primitive')


def test_no_loop_reads_an_unordered_container():
    """The other half of determinism, asserted over the source because a hash seed cannot test it.

    Every loop in the module must run over atoms sorted by stable id or bonds sorted by their
    `(low, high)` pair; walking the open-bond set as laid out changes the answer on thirteen corpus
    molecules.  But `PYTHONHASHSEED` does not move a set of small integer tuples, so the behavioural
    test above passes on a version that iterates it raw.  Here a set-valued name may not be the subject
    of a `for`, a comprehension, or a `list()`/`tuple()` that freezes its layout; `sorted()` is fine.
    """
    tree = parse(SOURCE.read_text(encoding='utf-8'))

    def is_set(node) -> bool:
        """A set literal, a set comprehension, or a `set()`/`frozenset()` call."""
        return (isinstance(node, (Set, SetComp))
                or isinstance(node, Call) and isinstance(node.func, Name)
                and node.func.id in ('set', 'frozenset'))

    # every name this module binds to a set, by annotation (`x: set[...] = ...`) or by value
    unordered = set()
    for node in walk(tree):
        if isinstance(node, AnnAssign) and isinstance(node.target, Name):
            annotation = node.annotation
            if (isinstance(annotation, Name) and annotation.id in ('Set', 'set')
                    or getattr(getattr(annotation, 'value', None), 'id', None) in ('Set', 'set')
                    or node.value is not None and is_set(node.value)):
                unordered.add(node.target.id)
        elif isinstance(node, Assign) and is_set(node.value):
            unordered.update(target.id for target in node.targets if isinstance(target, Name))

    def named(node) -> str:
        return node.id if isinstance(node, Name) else ''

    offences = []
    for node in walk(tree):
        if isinstance(node, For) and named(node.iter) in unordered:
            offences.append(f'{node.lineno}: for ... in {named(node.iter)}')
        elif isinstance(node, Call) and named(node.func) in ('list', 'tuple') and node.args \
                and named(node.args[0]) in unordered:
            offences.append(f'{node.lineno}: {named(node.func)}({named(node.args[0])})')
        for generator in getattr(node, 'generators', ()):
            if named(generator.iter) in unordered:
                offences.append(f'{node.lineno}: comprehension over {named(generator.iter)}')

    assert unordered, 'the scanner found no set at all, so it is proving nothing'
    assert not offences, (
        'a loop in saturate.py reads an unordered container, so the answer is a function of set '
        'layout rather than of the stable ids:\n  ' + '\n  '.join(sorted(offences))
        + '\n\nSort it once and read the sorted sequence: chython 2 produced a different molecule on '
        'different runs and this is the other way to get there.')


def test_the_search_order_is_the_stable_id_order():
    """The documented tie-break -- open bonds walked in sorted `(low id, high id)` order.

    These two fixtures are measured, not obvious: every other corpus entry comes back the same under
    both mutations.  Toluene catches dropping the `sorted()` and walking the open-bond set as laid out;
    pyrene is the only molecule whose answer is not symmetric under reversing the sort.
    """
    toluene, _ = prepared('Cc1ccccc1')
    assert saturate(toluene)
    assert orders(toluene) == {(1, 2): 1, (2, 3): 1, (2, 7): 2, (3, 4): 2, (4, 5): 1, (5, 6): 2,
                               (6, 7): 1}

    pyrene, _ = prepared('c1cc2ccc3cccc4ccc(c1)c2c34')
    assert saturate(pyrene)
    assert orders(pyrene) == {
        (1, 2): 1, (1, 14): 2, (2, 3): 2, (3, 4): 1, (3, 15): 1, (4, 5): 2, (5, 6): 1, (6, 7): 1,
        (6, 16): 2, (7, 8): 2, (8, 9): 1, (9, 10): 2, (10, 11): 1, (10, 16): 1, (11, 12): 2,
        (12, 13): 1, (13, 14): 1, (13, 15): 2, (15, 16): 1}


# --- what it must never do to the molecule it was given ------------------------------------------- #

@pytest.mark.parametrize('name', sorted(CORPUS))
def test_charges_and_radicals_are_never_touched(name):
    """Charges and radicals are input, never rewritten to make the search succeed."""
    molecule, _ = prepared(CORPUS[name])
    charges = {n: molecule.charge_of(n) for n in molecule.atom_numbers}
    radicals = {n: molecule.radical_of(n) for n in molecule.atom_numbers}
    saturate(molecule)
    assert {n: molecule.charge_of(n) for n in molecule.atom_numbers} == charges
    assert {n: molecule.radical_of(n) for n in molecule.atom_numbers} == radicals


@pytest.mark.parametrize('name', sorted(CORPUS))
def test_no_bond_is_created_or_deleted(name):
    """The bond set is untouched: no bond is added and none removed to help the balance."""
    molecule, _ = prepared(CORPUS[name])
    before = set(orders(molecule))
    saturate(molecule)
    assert set(orders(molecule)) == before


@pytest.mark.parametrize('name', sorted(CORPUS))
def test_hydrogen_counts_are_read_and_never_written(name):
    """The counts are the input.  Deriving one is `calc_implicit`'s job and the caller's call."""
    molecule, _ = prepared(CORPUS[name])
    before = {n: molecule.implicit_h_of(n) for n in molecule.atom_numbers}
    saturate(molecule)
    assert {n: molecule.implicit_h_of(n) for n in molecule.atom_numbers} == before


def test_no_order_is_ever_lowered():
    """A multiple bond already in the molecule is a stated fact, so the pass only raises.

    Two halves: a fully ordered molecule has no open bond, so the first alone cannot see the write
    path.  The second puts an existing double bond inside an open fragment (`CC=N` with both stated
    counts zero), so the answer is only reachable by raising 2 to 3 -- a write that replaced the order
    instead of adding to it would leave the double bond alone and the molecule unsatisfied.
    """
    molecule = read_smiles('CC(=O)C=CC#N')
    before = orders(molecule)
    assert saturate(molecule)
    assert orders(molecule) == before

    molecule = read_smiles('CC=N')
    with molecule.edit():
        molecule.set_hydrogens(2, 0)
        molecule.set_hydrogens(3, 0)
    assert saturate(molecule)
    assert orders(molecule) == {(1, 2): 1, (2, 3): 3}


def test_a_partly_ordered_skeleton_is_not_refused():
    """A molecule that already carries multiple bonds is answered, not refused.  Nothing raises."""
    molecule = read_smiles('CC=CC=O')
    log = molecule.log
    assert saturate(molecule)
    assert molecule.smiles == 'O=CC=CC'
    assert not log, [str(record) for record in log]


# --- honest failure ------------------------------------------------------------------------------- #

def test_an_unsatisfiable_fragment_returns_false():
    """An unbalanced fragment is `False` plus a log, never a silently unbalanced answer.

    Five carbons in a ring, each stated to carry one hydrogen: every one needs exactly one unit of
    extra bond order and five is odd, so no assignment exists.
    """
    molecule = read_smiles('C1CCCC1')
    with molecule.edit():
        for n in molecule.atom_numbers:
            molecule.set_hydrogens(n, 1)
    before = orders(molecule)
    log = molecule.log
    assert saturate(molecule) is False
    assert orders(molecule) == before, 'a refused fragment must keep every one of its bonds'

    refusals = [r for r in log if r.rule == 'saturate:no-valence-state']
    assert len(refusals) == 1 and refusals[0].atoms == (1, 2, 3, 4, 5)
    assert refusals[0].severity == REFUSED
    assert {r.atoms[0] for r in log if r.rule == 'saturate:unsatisfied'} == {1, 2, 3, 4, 5}


def test_an_atom_no_valence_row_admits_is_reported_and_frozen():
    """The other branch of the same rule id: an atom refused before any neighbour is looked at.

    Trimethyl oxonium drawn neutral -- oxygen at charge 0 has no row at bond order sum 3 whatever its
    neighbours do.  Reported exactly once is what pins the freeze: an atom left in the search would be
    named again under `saturate:unsatisfied` and its bonds would stay open.
    """
    molecule = read_smiles('O(C)(C)C')
    before = orders(molecule)
    log = molecule.log
    assert saturate(molecule) is False
    assert orders(molecule) == before
    assert [(r.rule, r.atoms, r.severity) for r in log] == [
        ('saturate:no-valence-state', (1,), REFUSED)]
    assert 'no valence row accepts O' in log[0].message


def test_a_fragment_the_environment_column_refuses_keeps_every_bond():
    """The exact test is asked inside the search, so this fragment is refused rather than written.

    A chlorine drawn with two oxygens and no charge: order-sum arithmetic offers `Cl(=O)=O`, but that
    row's environment column demands three `=O` neighbours.  Accepting on the arithmetic and writing
    would give half an answer -- written orders plus an unsatisfied verdict -- so the orders are
    asserted bond for bond and not merely that the return is `False`.
    """
    molecule, _ = prepared('C[Cl](=O)=O')
    log = molecule.log
    assert saturate(molecule) is False
    assert orders(molecule) == {(1, 2): 1, (2, 3): 1, (2, 4): 1}
    refusals = [r for r in log if r.rule == 'saturate:no-valence-state']
    assert len(refusals) == 1 and refusals[0].atoms == (2, 3, 4)
    assert refusals[0].severity == REFUSED
    assert not [r for r in log if r.rule == 'saturate:orders'], 'a refused fragment was written'


def test_asking_the_collection_exactly_lets_the_search_keep_looking():
    """The same test that refuses the fragment above makes this one succeed.

    Perchloric acid's skeleton: rejecting an assignment the collection does not accept is not a veto on
    the fragment, the search backtracks and goes on.  A rollback after the write could only undo.
    """
    molecule, answer = prepared('OCl(=O)(=O)=O')
    log = molecule.log
    assert saturate(molecule) is True
    assert orders(molecule) == orders(answer)
    assert [r.rule for r in log] == ['saturate:orders']


def test_propagation_runs_to_a_fixpoint_and_not_a_single_sweep():
    """Closing one atom's bonds can starve a neighbour, whose bonds then close in turn.

    Hydroxylamine's skeleton with the nitrogen stating no hydrogens.  The first sweep closes the
    oxygens; only the second sees that nothing incident to the nitrogen is left to raise.  The rule id
    is asserted rather than the return, because a single sweep would report `unsatisfied` instead.
    """
    molecule = read_smiles('ONO')
    with molecule.edit():
        molecule.set_hydrogens(2, 0)
    before = orders(molecule)
    log = molecule.log
    assert saturate(molecule) is False
    assert orders(molecule) == before
    assert [(r.rule, r.atoms, r.severity) for r in log] == [
        ('saturate:no-valence-state', (2,), REFUSED)]
    assert 'more bond order than its neighbours can accept' in log[0].message


def test_the_search_budget_refuses_when_it_found_nothing(monkeypatch):
    """Cut the budget below the first solution and the fragment is refused whole, like any other.

    `_NODES_MAX` is sized for input nobody has yet seen, so lowering it is the only way to execute the
    branch.  A truncated search is a refusal and not a truncated answer: every bond stays as it was.
    """
    monkeypatch.setattr(_saturate, '_NODES_MAX', 1)
    molecule, _ = prepared('c1ccccc1')
    before = orders(molecule)
    log = molecule.log
    assert saturate(molecule) is False
    assert orders(molecule) == before
    budget = [r for r in log if r.rule == 'saturate:budget']
    assert len(budget) == 1 and budget[0].severity == REFUSED
    assert 'left as it was' in budget[0].message
    assert not [r for r in log if r.rule == 'saturate:orders']


def test_the_search_budget_keeps_an_answer_it_could_not_prove_unique(monkeypatch):
    """The sibling branch, which is why the two severities differ.

    With enough budget to reach a solution but not to look for a second, the answer is written and the
    log says only that uniqueness went unchecked -- `LOST`, not `REFUSED`, since a caller filtering for
    refusals wants fragments it has to deal with itself.
    """
    monkeypatch.setattr(_saturate, '_NODES_MAX', 8)
    molecule, answer = prepared('c1ccccc1')
    log = molecule.log
    assert saturate(molecule) is True
    assert aromatic_form(molecule) == aromatic_form(answer)
    budget = [r for r in log if r.rule == 'saturate:budget']
    assert len(budget) == 1 and budget[0].severity == LOST
    assert 'not as the only answer' in budget[0].message
    assert [r.severity for r in log if r.rule == 'saturate:orders'] == [INFO]


def test_an_impossible_charge_state_is_reported():
    """A mis-drawn pentavalent nitro: the stated charges make saturation impossible.

    It is a finding, not a rewrite -- the two oxygens are named, the nitrogen keeps its charge, and
    `standardize()` is the pass that turns this drawing into `C[N+]([O-])=O`.  Both records come from
    propagation (starved), not from the per-atom read, which is the other branch of the same rule id in
    `test_an_atom_no_valence_row_admits_is_reported_and_frozen`.
    """
    molecule, _ = prepared('CN(=O)=O')
    log = molecule.log
    assert saturate(molecule) is False
    assert {n: molecule.charge_of(n) for n in molecule.atom_numbers} == {1: 0, 2: 0, 3: 0, 4: 0}
    starved = [r for r in log if r.rule == 'saturate:no-valence-state']
    assert {r.atoms for r in starved} == {(3,), (4,)}
    assert all('its neighbours can accept' in r.message for r in starved)
    # and the charged spelling of the same group is recovered exactly
    charged, answer = prepared('C[N+](=O)[O-]')
    assert saturate(charged)
    assert orders(charged) == orders(answer)


def test_a_failed_fragment_does_not_block_another():
    """All-or-nothing per fragment: one fragment's refusal does not cancel another's answer."""
    molecule, _ = prepared('CN(=O)=O.CC(C)=O')
    log = molecule.log
    assert saturate(molecule) is False
    assert molecule.smiles == 'CN([O])[O].C(C)(=O)C', molecule.smiles
    assert [r.rule for r in log].count('saturate:orders') == 1


def test_an_aromatic_bond_is_reported_rather_than_valued():
    """Order 4 has no valence row, so an atom carrying one is left alone and named.  Kekulise first."""
    molecule = read_smiles('c1ccccc1')
    before = orders(molecule)
    log = molecule.log
    assert saturate(molecule) is False
    assert orders(molecule) == before
    assert {r.atoms[0] for r in log if r.rule == 'saturate:aromatic-bond'} == {1, 2, 3, 4, 5, 6}
    assert all(r.severity == LOST for r in log if r.rule == 'saturate:aromatic-bond')


def test_a_state_the_collection_does_not_describe_is_a_gap_and_not_a_violation():
    """`'unknown'` is a hole in the collection and no claim about the molecule.  Iron at charge +1."""
    molecule = read_smiles('C[Fe+]')
    log = molecule.log
    assert saturate(molecule) is False
    gaps = [r for r in log if r.rule == 'saturate:collection-gap']
    assert len(gaps) == 1 and gaps[0].atoms == (2,) and gaps[0].severity == LOST
    assert 'gap in the collection' in gaps[0].message


def test_an_oversized_fragment_is_refused_and_says_why():
    """The ligand-scale boundary, executable.  A long polyene is not what this pass is for.

    Every carbon states one hydrogen, so every bond stays open and the fragment is one problem 300
    bonds wide.  What must not happen is an unbounded search or a recursion limit.
    """
    molecule, _ = prepared('C' + '=CC' * 150)
    before = orders(molecule)
    log = molecule.log
    assert saturate(molecule) is False
    assert orders(molecule) == before
    oversized = [r for r in log if r.rule == 'saturate:oversized']
    assert len(oversized) == 1 and oversized[0].severity == REFUSED
    assert 'ligand' in oversized[0].message


def test_a_long_chain_whose_counts_are_stated_never_reaches_the_size_cap():
    """The cap bounds the search, not the input.

    660 carbons each stating a count: every one is already at the only order sum its count admits, so
    propagation closes every bond and no fragment is left to search.
    """
    molecule = read_smiles('C' * 660)
    log = molecule.log
    assert saturate(molecule) is True
    assert set(orders(molecule).values()) == {1}
    assert not log, [str(record) for record in log]


def test_a_long_chain_with_no_stated_counts_is_answered_unforced_at_any_size():
    """The other reading of a chain, the one a PDB delivers: nothing in it demands a raise.

    The unforced test deliberately runs BEFORE the size cap -- otherwise a correct cheap answer would
    be traded for a refusal on a 659-bond fragment.  This pins that ordering.
    """
    molecule = read_smiles('C' * 660)
    with molecule.edit():
        for n in molecule.atom_numbers:
            molecule.set_hydrogens(n, H_UNKNOWN)
    log = molecule.log
    assert saturate(molecule) is True
    assert set(orders(molecule).values()) == {1}
    assert [(r.rule, len(r.atoms), r.severity) for r in log] == [('saturate:unforced', 660, INFO)]


# --- ambiguity, and the hydrogen count as the only signal ----------------------------------------- #

def test_ambiguity_is_reported_with_the_atoms_that_differ():
    """Benzene has two Kekule forms.  One is returned and the log says the choice was not unique."""
    molecule, _ = prepared('c1ccccc1')
    log = molecule.log
    assert saturate(molecule)
    ambiguous = [r for r in log if r.rule == 'saturate:ambiguous']
    assert len(ambiguous) == 1
    assert ambiguous[0].atoms == (1, 2, 3, 4, 5, 6)
    assert ambiguous[0].severity == LOST
    assert 'thiele()' in ambiguous[0].message


def test_a_unique_answer_is_not_reported_as_ambiguous():
    """Without this the test above passes for a pass that cries ambiguity on everything."""
    for smiles in ('CC(C)=O', 'CC#N', 'CS(=O)(=O)N', 'OP(=O)(O)O', 'c1cc[nH]c1', 'c1c[nH]cn1'):
        molecule, _ = prepared(smiles)
        log = molecule.log
        assert saturate(molecule)
        assert not [r for r in log if r.rule == 'saturate:ambiguous'], smiles


def test_a_skeleton_with_no_stated_hydrogens_is_reported_and_left_alone():
    """The CONECT-only PDB ligand: with no counts and no coordinates nothing demands a raise.

    Benzene's skeleton is also cyclohexane's, and this says so instead of guessing.  The log must be
    exactly one line: an orders comparison cannot tell "reported and skipped" from "searched and found
    nothing to do", and a search would also add a `saturate:ambiguous` line, turning a fragment nobody
    can derive into one reported as merely not unique.
    """
    molecule, _ = prepared('c1ccccc1')
    with molecule.edit():
        for n in molecule.atom_numbers:
            molecule.set_hydrogens(n, H_UNKNOWN)
    assert saturate(molecule) is True               # every bond single is a legal state
    assert orders(molecule) == {(1, 2): 1, (2, 3): 1, (3, 4): 1, (4, 5): 1, (5, 6): 1, (1, 6): 1}
    # the container's log accumulates, and `prepared()` kekulised: this pass's own lines are its stage
    log = [r for r in molecule.log if r.stage == 'saturate']
    assert [r.rule for r in log] == ['saturate:unforced']
    assert log[0].atoms == (1, 2, 3, 4, 5, 6)
    # INFO: nothing was lost and nothing refused; the line tells a caller their ligand needs
    # hydrogens, a residue template or human eyes
    assert log[0].severity == INFO


def test_an_unforced_fragment_can_still_hold_an_atom_nothing_satisfies():
    """`unforced` says no atom demands a raise.  It does not say every assignment would be legal.

    The neutral chlorine with two oxygens, now with no stated counts: the fragment is answered unforced
    and left single, and the chlorine is still in a state no valence row accepts.  The verdict's word
    cannot come from `valence_check` -- with no count there is no complete question -- so it is
    `violation`, which is the caller's difference between bad input and a coverage hole.
    """
    molecule, _ = prepared('C[Cl](=O)=O')
    with molecule.edit():
        for n in molecule.atom_numbers:
            molecule.set_hydrogens(n, H_UNKNOWN)
    log = molecule.log
    assert saturate(molecule) is False
    assert orders(molecule) == {(1, 2): 1, (2, 3): 1, (2, 4): 1}
    assert [(r.rule, r.atoms) for r in log] == [('saturate:unforced', (1, 2, 3, 4)),
                                                ('saturate:unsatisfied', (2,))]
    assert log[1].severity == LOST
    assert 'an unstated hydrogen count' in log[1].message
    assert 'calls a violation' in log[1].message


def test_one_stated_hydrogen_count_is_enough_to_force_a_bond():
    """The signal is per atom, so a partly-annotated skeleton is partly derivable.

    Acetone's skeleton with only the carbonyl oxygen's count stated: it must reach order sum two, which
    no hydrogen of its own can supply, so the C=O is forced while the three carbons stay slack.
    """
    molecule, answer = prepared('CC(C)=O')
    with molecule.edit():
        for n in (1, 2, 3):
            molecule.set_hydrogens(n, H_UNKNOWN)
    log = molecule.log
    assert saturate(molecule)
    assert orders(molecule) == orders(answer)
    # INFO: what was assigned is not damage, so a caller filtering for what it must look at itself
    # does not have to read every success line
    assert [(r.rule, r.severity) for r in log] == [('saturate:orders', INFO)]


# --- the shape of the pass ------------------------------------------------------------------------ #

def test_the_signature_asks_for_nothing_but_a_molecule():
    """One parameter only: no expected charge, no radical count, no flag that hides a defect -- and no
    `log=`, the molecule being where the records go."""
    assert list(signature(saturate).parameters) == ['molecule']


def test_implicit_hydrogens_are_enough():
    """No explicit hydrogen atom is required anywhere, and the pass does not quietly add one.

    A PDB ligand never has them, and the whole corpus above is implicit-hydrogen input.
    """
    molecule, answer = prepared('CC(=O)Nc1ccc(O)cc1')
    assert saturate(molecule)
    assert not any(molecule.explicit_h_of(n) for n in molecule.atom_numbers)
    assert aromatic_form(molecule) == aromatic_form(answer)


def test_explicit_hydrogens_work_too():
    """The XYZ shape works too -- hydrogens as atoms, heavy atoms stating none of their own."""
    molecule = read_smiles('CC(=O)Nc1ccc(O)cc1')
    molecule.kekule()
    molecule.explicify_hydrogens()
    answer = read_smiles('CC(=O)Nc1ccc(O)cc1')
    answer.kekule()
    answer.explicify_hydrogens()
    flatten(molecule)
    assert saturate(molecule)
    assert aromatic_form(molecule) == aromatic_form(answer)


def test_the_log_is_optional_and_the_return_is_the_verdict():
    """`log=None` is the default and the return says whether every atom ended up satisfied."""
    molecule, _ = prepared('CC(C)=O')
    assert saturate(molecule) is True
    assert saturate(molecule) is True                        # idempotent: nothing left to raise

    broken = read_smiles('c1ccccc1')
    assert saturate(broken) is False                         # and False without a log to say why
