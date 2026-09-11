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
"""The twelve reaction-level passes: a loop over molecules, plus the four things a loop cannot do.

`ReactionContainer`'s methods are one line each and forward here.  They live in their own module and
not in `reaction.py` because that file is an argument about a container's IDENTITY -- what equality
compares, what pach can store, what the ML view reports -- and a standardization chapter appended to
it would bury that.  Nothing here is imported by `reaction.py` at its own import time except this
module, and this module imports only the extension, so the layering is unchanged.

EIGHT OF THE TWELVE ARE A LOOP.  `standardize`, `kekule`, `thiele`, `canonicalize`, `neutralize`,
`implicify_hydrogens`, `clean_isotopes` and `clean_stereo` do to a reaction exactly what the molecule
method does to a molecule, once per molecule, and their whole content is the aggregation of the answer
and the log.  Seven fold their answer into one `bool` or one count; `clean_stereo` keys the molecule's
report by location instead, because five readers do not fold into a bool without the caller losing what
was wiped.  The other four are where a reaction is not a bag of molecules:

  * `explicify_hydrogens` has to give a hydrogen added on the left and the matching hydrogen added on
    the right ONE map number, or the mapping it hands back says a C-H bond broke and an identical one
    formed;
  * `remove_reagents` is a statement about the relation between the sides and has no molecule-level
    counterpart at all;
  * `contract_ions` pairs a cation on one side with an anion on the same side;
  * `reset_mapping` allocates from one counter across every molecule, which is the only reason it
    cannot be `mol.reset_mapping()` three times.

THE LOG.  NO PASS HERE TAKES A `log=`, and neither does the molecule pass it calls: a pass writes to
the container it was given, `rxn.log` and `mol.log` are how a caller reads it, and there is nothing to
switch off.  A component keeps its own records and `rxn.log` gets a copy -- `_mirror` below is the only
mechanism -- so `rxn.products[0].log` and `rxn.log.by_subject('products[0]')` answer the same question
from the two ends.

The copy is what needs a subject.  A `LogRecord`'s `atoms` are stable ids *in one container*, so a
reaction log pooling three sides' records unstamped would hand back numbers that name a different atom
depending on which molecule you read them against.  `LogRecord.subject` is that field; `_mirror` sets it
to the molecule's location -- `'reactants[0]'`, `'agents[1]'`.  Nothing is encoded into `rule`, which
would be a second spelling of a fact the substrate already has a field for.

WHAT THIS MODULE ASKS OF `MoleculeContainer`, in one place so that the pass still waiting on it has one
contract to satisfy rather than four guesses:

    standardize(*, fix_hydrogens=True, fix_tautomers=True) -> bool               landed
    kekule() -> KekuleResult(changed, log, unresolved)                           landed
    thiele() -> ThieleResult(changed, log, refused)                              landed
    canonicalize(*, fix_tautomers=True, keep_kekule=False) -> bool               landed
    implicify_hydrogens() -> int   (hydrogen ATOMS removed)                      landed
    explicify_hydrogens() -> int   (hydrogen ATOMS added)                        landed
    neutralize(*, keep_charge=True) -> bool                                      landed

`kekule`, `thiele`, `clean_isotopes` and `clean_stereo` are the core's own; the rest are registered onto
the container by `chython.chemistry`, so with only the core imported a reaction pass raises `ImportError`
naming that package rather than `AttributeError` on a method that looks absent.

The two hydrogen ones are asked for a COUNT and nothing else, and there is no `_return_map=True` /
`start_map=n` protocol to hand back allocated map numbers.  An atom's `n` and `map_number` are separate
fields, a new hydrogen arrives unmapped, and the reaction layer finds the ones it has to number by
comparing the molecule's stable ids before and after.  So the molecule pass needs no private keyword and
no knowledge that a reaction exists -- and a test in the standardization pack pins its signature at
exactly `{molecule}` so none can creep back.
"""
from collections.abc import Sequence
from contextlib import contextmanager
from itertools import chain, count

from ._core import MoleculeContainer
from ._log import LogRecord


__all__ = ['reaction_canonicalize', 'reaction_clean_isotopes', 'reaction_clean_stereo',
           'reaction_contract_ions', 'reaction_explicify_hydrogens', 'reaction_implicify_hydrogens',
           'reaction_kekule', 'reaction_neutralize', 'reaction_number_new_hydrogens',
           'reaction_remove_reagents', 'reaction_reset_mapping', 'reaction_standardize',
           'reaction_thiele']


_RULE_HYDROGENS = 'reaction:explicify_hydrogens'


def _located(rxn):
    """`(location, molecule)` in `molecules()` order, which is the order every log is written in.

    The location string is the `subject` every record from that molecule is stamped with, and it is
    also how a caller addresses the molecule again: `'products[0]'` is `rxn.products[0]`.
    """
    out = []
    for side, molecules in (('reactants', rxn.reactants), ('agents', rxn.agents),
                            ('products', rxn.products)):
        for index, molecule in enumerate(molecules):
            out.append(('%s[%d]' % (side, index), molecule))
    return out


@contextmanager
def _mirror(rxn, molecule, subject):
    """Copy onto `rxn.log`, stamped with `subject`, whatever the block wrote to `molecule.log`.

    THE COMPONENT KEEPS ITS OWN RECORDS AND THE REACTION GETS A COPY, which is the only arrangement in
    which both `rxn.log.by_subject('products[0]')` and `rxn.products[0].log` answer.  The copy is what
    needs a subject: `LogRecord.atoms` are stable ids in ONE container, so a reaction log pooling three
    sides' records unstamped hands back numbers that name a different atom depending on which molecule
    you read them against.  `stage` is not touched -- the molecule pass already named it, and it is
    more precise than anything this layer knows.
    """
    mine = molecule.log
    start = len(mine)
    try:
        yield mine
    finally:
        with rxn.log.stage('', subject=subject) as log:
            log.extend(mine[start:])


# --------------------------------------------------------------------------------------------------
# the eight that are a loop

def reaction_standardize(rxn, *, fix_hydrogens: bool = True, fix_tautomers: bool = True) -> bool:
    """`ReactionContainer.standardize`."""
    changed = False
    for where, molecule in _located(rxn):
        with _mirror(rxn, molecule, where):
            if molecule.standardize(fix_hydrogens=fix_hydrogens, fix_tautomers=fix_tautomers):
                changed = True
    return changed


def reaction_canonicalize(rxn, *, fix_tautomers: bool = True, keep_kekule: bool = False) -> bool:
    """`ReactionContainer.canonicalize`."""
    changed = False
    for where, molecule in _located(rxn):
        with _mirror(rxn, molecule, where):
            if molecule.canonicalize(fix_tautomers=fix_tautomers, keep_kekule=keep_kekule):
                changed = True
    return changed


def reaction_neutralize(rxn, *, keep_charge: bool = True) -> bool:
    """`ReactionContainer.neutralize`."""
    changed = False
    for where, molecule in _located(rxn):
        with _mirror(rxn, molecule, where):
            if molecule.neutralize(keep_charge=keep_charge):
                changed = True
    return changed


def reaction_kekule(rxn) -> bool:
    """`ReactionContainer.kekule`."""
    changed = False
    for where, molecule in _located(rxn):
        with _mirror(rxn, molecule, where):
            if molecule.kekule().changed:
                changed = True
    return changed


def reaction_thiele(rxn) -> bool:
    """`ReactionContainer.thiele`."""
    changed = False
    for where, molecule in _located(rxn):
        with _mirror(rxn, molecule, where):
            if molecule.thiele().changed:
                changed = True
    return changed


def reaction_implicify_hydrogens(rxn) -> int:
    """`ReactionContainer.implicify_hydrogens`."""
    total = 0
    for where, molecule in _located(rxn):
        with _mirror(rxn, molecule, where):
            total += molecule.implicify_hydrogens()
    return total


def reaction_clean_isotopes(rxn) -> bool:
    """`ReactionContainer.clean_isotopes`."""
    changed = False
    for where, molecule in _located(rxn):
        with _mirror(rxn, molecule, where):
            if molecule.clean_isotopes():
                changed = True
    return changed


def reaction_clean_stereo(rxn) -> dict:
    """`ReactionContainer.clean_stereo`.

    THE ONE PASS HERE WHOSE ANSWER IS NOT A BOOL OR A COUNT.  The molecule's report names five readers'
    worth of state and its ids are ids in ONE container, so the aggregation is the same dict the log
    copy uses -- keyed by location, and a molecule that carried no stereo is ABSENT rather than mapped
    to `{}`, exactly as an empty reader is absent from the molecule's own report.
    """
    report = {}
    for where, molecule in _located(rxn):
        with _mirror(rxn, molecule, where):
            wiped = molecule.clean_stereo()
        if wiped:
            report[where] = wiped
    return report


# --------------------------------------------------------------------------------------------------
# hydrogens: the count is the molecule's, the numbering is the reaction's

def reaction_explicify_hydrogens(rxn) -> int:
    """`ReactionContainer.explicify_hydrogens`."""
    located = _located(rxn)
    before = [set(molecule.atom_numbers) for _, molecule in located]
    total = 0
    for (where, molecule), old in zip(located, before):
        with _mirror(rxn, molecule, where):
            total += molecule.explicify_hydrogens()
    if total:
        reaction_number_new_hydrogens(rxn, [set(molecule.atom_numbers) - old
                                            for (_, molecule), old in zip(located, before)])
    return total


def reaction_number_new_hydrogens(rxn, new_ids: Sequence[set[int]]) -> int:
    """Give each newly added hydrogen a map number, pairing across the arrow.  How many were numbered?

    `new_ids` is one set of stable ids per molecule, in `molecules()` order -- the atoms that were just
    added.  Split out from :func:`reaction_explicify_hydrogens` because it is the part that belongs to
    the reaction rather than to any molecule, and splitting it is what makes it testable before the
    molecule pass that feeds it exists.

    A HYDROGEN IS NUMBERED ONLY IF ITS MOLECULE IS ALREADY FULLY MAPPED, and that single condition
    replaces the guesswork.  Numbering a hydrogen in an unmapped molecule would invent a mapping the
    record never claimed; leaving one unmapped in a fully mapped molecule would make that molecule
    PARTIALLY mapped, which `reaction_pach_dump` refuses by name at version 1, whose format has one
    number field per atom.  So the rule is per molecule, and an unmapped reaction comes out of here untouched.

    THE PAIRING IS BY THE HEAVY ATOM'S MAP NUMBER and nothing else.  Hydrogens on one heavy atom are
    interchangeable -- they differ by an automorphism of that atom's environment -- so which product
    hydrogen inherits which reactant's number is arbitrary and taking them in order is as good as any
    other choice.  A hydrogen on a heavy atom that the other side does not have gets a fresh number,
    which is what "this hydrogen is not the same hydrogen" is spelled as.
    """
    located = _located(rxn)
    if len(new_ids) != len(located):
        raise ValueError('new_ids holds %d entry/entries and the reaction has %d molecule(s); it is '
                         'one set of stable ids per molecule in molecules() order'
                         % (len(new_ids), len(located)))

    highest = 0
    for _, molecule in located:
        for atom in molecule.atoms():
            if atom.map_number > highest:
                highest = atom.map_number
    fresh = count(highest + 1)

    left = len(rxn.reactants)
    middle = left + len(rxn.agents)
    # heavy map number -> the numbers given to hydrogens hanging off it on the REACTANT side
    pool: dict[int, list[int]] = {}
    numbered = 0
    for position, ((where, molecule), added) in enumerate(zip(located, new_ids)):
        if not added or not _fully_mapped(molecule, added):
            continue
        products = position >= middle
        assignments = []
        for n in sorted(added):
            if molecule.element_of(n) != 1:
                continue
            heavy = _lone_heavy_neighbour(molecule, n)
            if heavy is None:
                continue
            key = molecule.map_number_of(heavy)
            if not key:
                continue
            if products and pool.get(key):
                assignments.append((n, pool[key].pop(0), key, True))
            else:
                number = next(fresh)
                assignments.append((n, number, key, False))
                if position < left:
                    pool.setdefault(key, []).append(number)
        if not assignments:
            continue
        with molecule.edit():
            for n, number, _, _ in assignments:
                molecule.set_map_number(n, number)
        numbered += len(assignments)
        with _mirror(rxn, molecule, where) as mine:
            with mine.stage('number_new_hydrogens'):
                for n, number, key, paired in assignments:
                    if paired:
                        message = ('new hydrogen on the atom mapped %d takes map number %d from the '
                                   'matching new hydrogen on the reactant side' % (key, number))
                    else:
                        message = ('new hydrogen on the atom mapped %d takes the fresh map number %d; '
                                   'the other side has no hydrogen to pair it with' % (key, number))
                    mine.append(LogRecord(_RULE_HYDROGENS, (n,), message))
    return numbered


def _fully_mapped(molecule, added) -> bool:
    """Does every atom of this molecule except the ones just added carry a map number?"""
    for atom in molecule.atoms():
        if atom.n in added:
            continue
        if not atom.map_number:
            return False
    return True


def _lone_heavy_neighbour(molecule, n) -> int | None:
    """The one atom `n` hangs off, or None when it has no neighbour or more than one."""
    neighbours = list(molecule.neighbors_of(n))
    if len(neighbours) != 1:
        return None
    return neighbours[0]


# --------------------------------------------------------------------------------------------------
# reset_mapping

def reaction_reset_mapping(rxn) -> bool:
    """`ReactionContainer.reset_mapping`."""
    molecules = list(rxn.molecules())
    numbers = [a.map_number for m in molecules for a in m.atoms()]
    if len(set(numbers)) == len(numbers) and 0 not in numbers:
        return False
    fresh = count(1)
    for molecule in molecules:
        ids = list(molecule.atom_numbers)
        if not ids:
            continue
        with molecule.edit():
            for n in ids:
                molecule.set_map_number(n, next(fresh))
    return bool(numbers)


# --------------------------------------------------------------------------------------------------
# the reaction centre, and remove_reagents on top of it

def _side_states(molecules):
    """`({map number: state}, {colliding map numbers})` for one side.

    A state is `(element, charge, radical, implicit h, {neighbour map number: order})` -- everything a
    CGR's dynamic atom and dynamic bond carried between them, which is what the reaction centre is
    defined against.  Bonds to an unmapped atom are left out: they cannot be compared across the arrow
    because there is nothing to compare them to.
    """
    states = {}
    collisions = set()
    for molecule in molecules:
        local = {}
        for atom in molecule.atoms():
            # `mn`/`mm` are MAP numbers here; `n`/`m` are atom numbers, and this function holds both.
            mn = atom.map_number
            if not mn:
                continue
            if mn in states:
                collisions.add(mn)
            local[atom.n] = mn
            states[mn] = (atom.element, atom.charge, atom.is_radical, atom.implicit_h, {})
        for bond in molecule.bonds():
            mn, mm = local.get(bond.n), local.get(bond.m)
            if mn is None or mm is None:
                continue
            states[mn][4][mm] = bond.order
            states[mm][4][mn] = bond.order
    return states, collisions


def reaction_center(rxn) -> set[int]:
    """The map numbers of the atoms this reaction changes.  Empty for a reaction with no mapping.

    An atom is in the centre when it appears on BOTH sides and something about it differs: its element,
    charge, radical state, implicit hydrogen count, or the map numbers and orders of its bonds.  It is
    computed without building a CGR: there is no CGR container on this release, and the question does
    not need one.

    AN ATOM PRESENT ON ONLY ONE SIDE IS NOT IN THE CENTRE, and that is a decision.  Calling it dynamic
    instead -- a bond that exists on one side and not the other -- keeps sodium hydroxide a REACTANT in
    `[Na+:1].[OH-:2].MeOAc >> AcOH`, because its atoms go missing from the products.  The documented
    purpose of `remove_reagents` is that NaOH becomes an agent there, so this reading is the one that
    makes the documented example work: a molecule the mapping does not follow through the arrow is
    present but untracked, which is what an agent is.  A genuine leaving group is unaffected -- it is
    part of a molecule that also has retained atoms, and the retained atoms put that molecule in the
    centre.

    A COLLIDING MAP NUMBER IS TREATED AS ACTIVE.  Two atoms on one side sharing a number make the
    comparison meaningless for that number, and the conservative answer keeps the molecule where the
    record put it rather than demoting it on evidence that does not exist.
    """
    reactants, r_collisions = _side_states(rxn.reactants)
    products, p_collisions = _side_states(rxn.products)
    active = r_collisions | p_collisions
    for n, before in reactants.items():
        after = products.get(n)
        if after is not None and before != after:
            active.add(n)
    return active


def _touches(molecule, center) -> bool:
    for atom in molecule.atoms():
        if atom.map_number in center:
            return True
    return False


def reaction_remove_reagents(rxn, *, keep_reagents: bool = False, mapping: bool = True,
                             common: Sequence[MoleculeContainer] | None = None) -> bool:
    """`ReactionContainer.remove_reagents`."""
    if mapping:
        return _remove_reagents_mapping(rxn, keep_reagents)
    return _remove_reagents_rules(rxn, keep_reagents, common)


def _apply_sides(rxn, reactants, products, demoted, keep_reagents) -> bool:
    """Write the three sides back, or refuse to when it would empty one of them.

    `demoted` goes AFTER the agents the record already had, and it is NOT DEDUPLICATED HERE: collecting
    reagents into a `set` drops the second equivalent of a solvent that appeared twice, and a reaction's
    stoichiometry is data this pass is not the place to lose.  The callers decide how many copies to
    demote and this function only writes down what they decided.
    """
    if not demoted:
        return False
    if not reactants or not products:
        # a reaction whose every reactant or every product looks like a reagent is a record this pass
        # cannot improve, and emptying a side would turn a bad reaction into a broken one.
        return False
    rxn._reactants = tuple(reactants)
    rxn._products = tuple(products)
    rxn._agents = tuple(chain(rxn._agents, demoted)) if keep_reagents else ()
    return True


def _remove_reagents_mapping(rxn, keep_reagents) -> bool:
    center = reaction_center(rxn)
    if not center:
        raise ValueError('this reaction has no reaction centre according to its atom-to-atom mapping, '
                         'so there is nothing to tell a reagent from a reactant; pass mapping=False '
                         'to use the rule-based door, which needs no mapping')
    reactants, products, demoted = [], [], []
    for molecule in rxn.reactants:
        (reactants if _touches(molecule, center) else demoted).append(molecule)
    for molecule in rxn.products:
        (products if _touches(molecule, center) else demoted).append(molecule)
    return _apply_sides(rxn, reactants, products, demoted, keep_reagents)


def _remove_reagents_rules(rxn, keep_reagents, common) -> bool:
    """The door for an unmapped record: a molecule on both sides, or one the caller calls common.

    `common` IS AN ARGUMENT AND NOT A TABLE HERE.  A set of solvent SMILES inside this pass -- water,
    the halogen acids, benzene, toluene, hexane, the low alcohols, formic and acetic acid, ethyl
    acetate, the ethers -- is chemistry knowledge, and chemistry knowledge in this release is a row in a
    table owned by the package that holds the tables, not a literal in `core`.  So the algorithm is here
    and the data is the caller's until the reactions package lands one; `None` runs the stage that needs
    no data.

    Matching is by MOLECULE EQUALITY, so a caller's `common` has to be in the same representation as
    the record: `c1ccccc1` and `C1=CC=CC=C1` are two molecules to `==` and one to a chemist.  Kekulise
    or aromatise both sides first.
    """
    if not rxn.reactants or not rxn.products:
        return False

    # stage 1: the same molecule on both sides is not part of the transformation.
    #
    # COUNTED, not set-membership, and the counting is the whole subtlety.  A molecule present twice on
    # the left and once on the right passed through ONCE -- the second equivalent was consumed -- so one
    # copy leaves each side and ONE agent is produced, not two.  Demoting per occurrence would report
    # two equivalents of a solvent where the record shows one, and demoting by set membership would
    # take the consumed equivalent with it.
    left_counts, right_counts = {}, {}
    for molecule in rxn.reactants:
        left_counts[molecule] = left_counts.get(molecule, 0) + 1
    for molecule in rxn.products:
        right_counts[molecule] = right_counts.get(molecule, 0) + 1
    #: how many copies of each molecule passed through untouched
    shared = {m: min(n, right_counts.get(m, 0)) for m, n in left_counts.items()}

    stage1_r, stage1_p, demoted = [], [], []
    budget = dict(shared)
    for molecule in rxn.reactants:
        if budget.get(molecule, 0):
            budget[molecule] -= 1
            demoted.append(molecule)
        else:
            stage1_r.append(molecule)
    budget = dict(shared)
    for molecule in rxn.products:
        if budget.get(molecule, 0):
            budget[molecule] -= 1  # its one agent was already demoted off the reactant side
        else:
            stage1_p.append(molecule)
    if not stage1_r or not stage1_p:
        return False  # every molecule appears on both sides; keep the bad record as it is

    if common is None:
        return _apply_sides(rxn, stage1_r, stage1_p, demoted, keep_reagents)

    # stage 2: and the ones the caller calls common, rolled back when it would empty a side
    common = list(common)
    stage2_r = [m for m in stage1_r if m not in common]
    stage2_p = [m for m in stage1_p if m not in common]
    if not stage2_r or not stage2_p:
        return _apply_sides(rxn, stage1_r, stage1_p, demoted, keep_reagents)
    demoted.extend(m for m in stage1_r if m in common)
    demoted.extend(m for m in stage1_p if m in common)
    return _apply_sides(rxn, stage2_r, stage2_p, demoted, keep_reagents)


# --------------------------------------------------------------------------------------------------
# contract_ions

def _sift_ions(molecules) -> tuple[list, list, list, int]:
    """`(neutral, cations, anions, the side's total charge)`."""
    neutral, cations, anions = [], [], []
    total = 0
    for molecule in molecules:
        charge = int(molecule)
        total += charge
        if charge > 0:
            cations.append(molecule)
        elif charge < 0:
            anions.append(molecule)
        else:
            neutral.append(molecule)
    return neutral, cations, anions, total


def _contract(anions, cations, total) -> list[MoleculeContainer] | None:
    """The salts one side's ions make, or None when which pairs with which is not determined.

    THREE CASES ARE THE WHOLE OF IT.  A side with a charge surplus can only be contracted when the
    minority is a single molecule -- otherwise there is no way to say which counter-ion belongs to which
    -- and a balanced side needs the anions all alike or the cations all alike for the same reason.
    Everything else is left as separate molecules, which is a refusal to guess and not a failure.
    """
    if not anions or not cations:
        return None
    if total > 0:
        if len(cations) > 1:
            return None
        salt = cations[0]
        for other in anions:
            salt = salt | other
        return [salt]
    elif total < 0:
        if len(anions) > 1:
            return None
        salt = anions[0]
        for other in cations:
            salt = salt | other
        return [salt]
    elif len(set(anions)) > 1 and len(set(cations)) > 1:
        return None

    salts = []
    anions = list(anions)
    cations = list(cations)
    while anions:
        salt = cations.pop() | anions.pop()
        while True:
            charge = int(salt)
            if charge > 0:
                salt = salt | anions.pop()
            elif charge < 0:
                salt = salt | cations.pop()
            else:
                break
        salts.append(salt)
    return salts


def reaction_contract_ions(rxn) -> bool:
    """`ReactionContainer.contract_ions`.

    EACH SIDE IS CONTRACTED ON ITS OWN, and no attempt is made to pair the reactant side's salts with
    the product side's.  The three sides are independent statements about what was in the flask, a salt
    that survives the reaction unchanged is the same molecule on both sides and therefore contracts the
    same way anyway, and the only cases where the sides could disagree are exactly the cases
    :func:`_contract` refuses on both.
    """
    changed = False
    for name in ('_reactants', '_agents', '_products'):
        neutral, cations, anions, total = _sift_ions(getattr(rxn, name))
        salts = _contract(anions, cations, total)
        if salts:
            setattr(rxn, name, tuple(chain(neutral, salts)))
            changed = True
    return changed
