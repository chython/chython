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
"""Template-driven mapping reconstruction: the recorded product, rebuilt from the recorded inputs.

A LADDER ORDERED BY STRENGTH OF EVIDENCE, not by cost: the first rung whose claim reproduces the
recorded product wins, so a record both a direct coupling and a protection could explain is read as the
coupling.  Every rung yields `_Explanation`s and writes nothing; the orchestrator applies exactly one,
so a half-match cannot leave half a mapping behind for the next rung.
"""
from collections.abc import Callable, Iterator, Sequence
from itertools import combinations
from typing import NamedTuple

from ._enumerate import _run, deprotect, EnumeratedReaction
from ._numbering import fast_mapping
from ._tables import reaction_rules
from ..core import INFO, LOST, LogRecord, MoleculeContainer, REFUSED, ReactionContainer


class _Options(NamedTuple):
    max_size_ratio: float
    min_filter_size: int


class _Explanation(NamedTuple):
    label: str
    write: Callable[[MoleculeContainer], bool]
    rule: str = ''          # the table-qualified row id, when a row produced this; logged, not returned


def reconstruct_mapping(reaction: ReactionContainer, *, max_size_ratio: float = 5.,
                        min_filter_size: int = 42) -> tuple[str, ...]:
    """See `ReactionContainer.reconstruct_mapping`, whose body this is."""
    log = reaction.log
    inputs = _inputs(reaction)
    if not inputs or not reaction.products:
        _refuse(log, 'reconstruct:empty',
                'a record with no inputs or no products has nothing to reconstruct from')
        return ()
    if len(reaction.products) > 1:
        _refuse(log, 'reconstruct:multiproduct',
                'a record with %d products is refused: which product a given input atom went to is a '
                'choice this makes silently and wrongly, and a wrong mapping is worse than none'
                % len(reaction.products))
        return ()

    # Canonicalize in place, then map: a template written against aromatic bonds cannot fire on a Kekule
    # record, and a mapping over a structure the caller is about to normalize is a mapping of something
    # else.
    reaction.canonicalize()
    inputs = _inputs(reaction)
    recorded = reaction.products[0]
    _clear(recorded)
    _number_inputs(inputs)

    options = _Options(max_size_ratio, min_filter_size)
    unbalanced = _grossly_unbalanced(recorded, inputs, options)
    if unbalanced:
        _refuse(log, 'reconstruct:unbalanced',
                'the recorded product has %d atoms against %d in every input together: the rungs that '
                'search the corpus are skipped, since the inputs do not contain the product'
                % (len(recorded.atom_numbers),
                   sum(len(molecule.atom_numbers) for molecule in inputs)))
    for phase in _PHASES:
        if unbalanced and phase not in _FILTER_EXEMPT:
            continue
        found = list(phase(recorded, inputs, options))
        if not found:
            continue
        if not found[0].write(recorded):
            _record(log, 'reconstruct:partial', INFO,
                    'some components of the recorded product were not reproduced and are left '
                    'unnumbered')
        if found[0].rule:
            # The label is the chemistry, the id is the row; only the id says which spelling earned the
            # hit.  The APPLIED rule only -- logging one merely considered is a false attribution.
            _record(log, 'reconstruct:rules', INFO, 'explained by %s' % found[0].rule)
        if len(found) > 1:
            alternatives = sorted({e.label for e in found[1:]})
            _record(log, 'reconstruct:alternatives', INFO,
                    'alternative explanations not applied: %s' % ', '.join(alternatives))
        _compact(reaction)
        return (found[0].label,)
    # Incoming numbers are deliberately not restored: a number the record arrived with is the claim
    # under test, and canonicalizing in place may already have invalidated it atom-for-atom.  And the
    # inputs' working numbers go with them -- a number that pairs with nothing is not a mapping.
    _compact(reaction)
    _record(log, 'reconstruct:unexplained', LOST,
            'no rung explained the record; the product\'s map numbers were cleared and not replaced')
    return ()


# --- the rungs ------------------------------------------------------------------------------------

def _purification(recorded, inputs, options) -> Iterator[_Explanation]:
    """The product went in and came out.  The mapping it implies is the identity.

    The one explanation whose mapping is certain, and the only label with no namespace -- no rule and no
    table produced it.
    """
    if _number_product(recorded.copy(), inputs):
        yield _Explanation('purification', lambda target: _number_product(target, inputs))


def _translate(reaction, sources: Sequence[MoleculeContainer]) -> list[MoleculeContainer]:
    """The reaction's products carrying `sources`' map numbers in place of the reactor's own.

    The two schemes compose through ATOM numbers: a reactor reactant is a copy of its source at the
    same atom numbers, so reactant map -> atom number -> source map.  `sources` is positional, which
    `_run` guarantees -- an outcome names the inputs its match touched, in input order.  An atom the
    reactor left at 0, and one no source claims, stays 0.
    """
    table = {}
    for reactant, source in zip(reaction.reactants, sources):
        for n in reactant.atom_numbers:
            number = reactant.map_number_of(n)
            if number:
                table[number] = source.map_number_of(n)
    out = []
    for product in reaction.products:
        product = product.copy()
        writes = [(n, table.get(product.map_number_of(n), 0)) for n in product.atom_numbers]
        with product.edit():
            for n, number in writes:
                product.set_map_number(n, number)
        out.append(product)
    return out


def _react(recorded, inputs, options) -> Iterator[_Explanation]:
    """A corpus row, applied to the inputs as they arrived.  The strongest evidence there is."""
    for pool, outcome in _applications(inputs):
        if not _reproduces(recorded, outcome.reaction.products):
            continue
        sources = [*_translate(outcome.reaction, pool), *inputs]
        yield _Explanation('react:%s' % outcome.name,
                           lambda target, s=sources: _number_product(target, s),
                           outcome.rule_id)


def _deprotect(recorded, inputs, options) -> Iterator[_Explanation]:
    """An input, unmasked.  The recorded product IS the recorded input minus a protecting group.

    `partial=True` because incomplete cleavage is ordinary and `R-N(Boc)2 -> R-NHBoc` is a record this
    rung must read.  The generator is largest-first, so nothing beyond what is taken is computed.
    """
    for molecule in inputs:
        for outcome in deprotect(molecule, partial=True):
            if not _reproduces(recorded, outcome.reaction.products):
                continue
            sources = [*_translate(outcome.reaction, [molecule]), *inputs]
            yield _Explanation('deprotect:%s' % '+'.join(outcome.names),
                               lambda target, s=sources: _number_product(target, s),
                               '+'.join(outcome.rule_ids))


def _deprotect_then_react(recorded, inputs, options) -> Iterator[_Explanation]:
    """Strip what can be stripped, then let the corpus fire on what is left.

    ONE all-stripped pass and not every raw/stripped combination, which would be exponential in the
    number of protected inputs: a deliberate lower bound on what composition buys.  A stripped form
    keeps the stable ids it came in with, so numbering survives the composition unaided.
    """
    pool = []
    changed = False
    for molecule in inputs:
        outcome = next(deprotect(molecule), None)
        if outcome is None:
            pool.append(molecule)
        else:
            pool.extend(_translate(outcome.reaction, [molecule]))
            changed = True
    if not changed:
        return
    for subset, outcome in _applications(pool):
        if not _reproduces(recorded, outcome.reaction.products):
            continue
        sources = [*_translate(outcome.reaction, subset), *pool, *inputs]
        yield _Explanation('deprotect+react:%s' % outcome.name,
                           lambda target, s=sources: _number_product(target, s),
                           outcome.rule_id)


def _protect(recorded, inputs, options) -> Iterator[_Explanation]:
    """The recorded product is a recorded input, masked.

    THE WEAKEST RUNG AND SO THE LAST: an amide, an ester and a carbamate are all protecting groups as
    well as products, so offered first this reads every acylation as a protection.
    """
    for outcome in deprotect(recorded, partial=True):
        pairs = _pair_with_inputs(outcome.reaction.products, inputs)
        if not pairs:
            continue
        yield _Explanation('protect:%s' % '+'.join(outcome.names),
                           lambda target, p=pairs: _number_from_pairs(target, p),
                           '+'.join(outcome.rule_ids))


_PHASES: tuple[Callable[..., Iterator[_Explanation]], ...] = (_purification, _react, _deprotect,
                                                              _deprotect_then_react, _protect)
_FILTER_EXEMPT = frozenset((_purification, _protect))


# --- writing the numbers --------------------------------------------------------------------------

def _number_product(recorded: MoleculeContainer, sources: Sequence[MoleculeContainer]) -> bool:
    """Number the components of `recorded` from the components of the numbered `sources`.

    Per connected component, which is what keeps a salt alive: a component no template touched is
    matched against the input it arrived in and numbered from there.  `split()` preserves stable ids, so
    the write goes straight onto `recorded`.  Each source component is consumed at most once.

    Returns True when every component was numbered; a partial answer is still written and the caller
    logs the shortfall.

    Writes onto the RECORDED containers, and `mapping_agrees` depends on it: that comparison excuses an
    automorphic swap by an orbit lookup keyed on the input's own stable ids, so rebuilt inputs would
    make every legitimate swap count as a disagreement.
    """
    parts = list(recorded.split())
    pool: list[MoleculeContainer] = [part for source in sources for part in source.split()]
    used = set()
    numbered = 0
    for part in parts:
        for i, candidate in enumerate(pool):
            if i in used:
                continue
            pairs = fast_mapping(candidate, part)
            if pairs is None:
                continue
            for source_id, target_id in pairs.items():
                recorded.set_map_number(target_id, candidate.map_number_of(source_id))
            used.add(i)
            numbered += 1
            break
    return numbered == len(parts)


def _pair_with_inputs(stripped, inputs):
    """`[(stripped component, the input component it equals)]`, each input component used once.

    Empty when nothing lines up, which is this rung's whole test.
    """
    pool = [part for molecule in inputs for part in molecule.split()]
    used = set()
    pairs = []
    for product in stripped:
        for part in product.split():
            for i, candidate in enumerate(pool):
                if i in used or candidate != part:
                    continue
                used.add(i)
                pairs.append((part, candidate))
                break
    return pairs


def _number_from_pairs(target, pairs) -> bool:
    """Copy map numbers from each input component onto the stripped component it matched.

    The other direction from `_number_product`: here the recorded product has been stripped, and the
    stripped fragment's stable ids are mostly a subset of the recorded product's own, so writing at those
    ids writes onto the right atoms.  Mostly -- `deprotect()` also BUILDS, e.g. the `=O` an acetal row
    restores, and a built atom's id belongs to nothing in `target`.  Hence the guard on `ids`; such an
    atom has no source and must stay at zero.

    Always True.  A protecting group's atoms are genuinely new, so zero is the right answer and not a
    shortfall; a cleaved fragment matching no input never enters `pairs` at all.
    """
    ids = set(target.atom_numbers)
    for part, candidate in pairs:
        mapped = fast_mapping(part, candidate)
        if mapped is None:                              # equal components, so this cannot happen
            continue
        for target_id, source_id in mapped.items():
            if target_id in ids:                        # an id the deprotection ADDED, see above
                target.set_map_number(target_id, candidate.map_number_of(source_id))
    return True


# --- enumeration ----------------------------------------------------------------------------------

def _applications(inputs: Sequence[MoleculeContainer],
                  rules=None) -> Iterator[tuple[list[MoleculeContainer], 'EnumeratedReaction']]:
    """Every way a corpus row applies to a SUBSET of `inputs`, with the subset it applied to.

    Subsets and not the whole list, because `_run` requires every molecule it is handed to be touched
    while a recorded record files its base, solvent and catalyst among the inputs.  Bounded by the widest
    row's slot count, so this is a small fixed number of combinations and not a power set.

    The subset comes back because the reactor's products carry the reactor's own numbering: `_translate`
    needs the inputs, positionally, to cross back to theirs.

    TODO: the bare `except` is defensive, not load-bearing -- `_run` has not been observed to raise on
    the current corpus, but a row that consistently fails is invisible here.  Narrow it to the observed
    type once there is one.
    """
    if rules is None:
        rules = reaction_rules()
    widest = max((len(rule.groups) for family in rules.values() for rule in family), default=0)
    order = range(len(inputs))
    for size in range(1, min(widest, len(inputs)) + 1):
        for subset in combinations(order, size):
            pool = [inputs[i] for i in subset]
            try:
                for outcome in _run(pool, rules):
                    yield pool, outcome
            except Exception:
                continue


def _reproduces(recorded: MoleculeContainer, products: Sequence[MoleculeContainer]) -> bool:
    """True when a connected component of `recorded` is a component of one of `products`.

    Per component, because a template answers the reaction centre while the record carries the salt too;
    and by container equality, never by SMILES -- `__eq__` is the canonical form.
    """
    want = set(recorded.split())
    return any(part in want for product in products for part in product.split())


# --- housekeeping ---------------------------------------------------------------------------------

def _inputs(reaction) -> list[MoleculeContainer]:
    """Reactants and agents both: an agent atom that ends up in the product is one the mapping owes an
    answer for."""
    return [*reaction.reactants, *reaction.agents]


def _clear(molecule) -> None:
    """Zero every map number on the recorded product; a leftover one would read as an answer."""
    for n in molecule.atom_numbers:
        molecule.set_map_number(n, 0)


def _compact(reaction: ReactionContainer) -> None:
    """Reduce the record's numbering to the 1-1 mapping, in place.

    The same rule the reactor imposes, at the other producer of a mapping: contiguous from 1 over the
    numbers present on BOTH sides, 0 for everything else.  A spectator input and a leaving group are on
    one side only, so a number there pairs with nothing and is dropped rather than left to read as a
    pairing.  Walked in `molecules()` order, so 1 is on the first reactant.
    """
    left = set()
    for molecule in (*reaction.reactants, *reaction.agents):
        left |= {molecule.map_number_of(n) for n in molecule.atom_numbers}
    right = set()
    for molecule in reaction.products:
        right |= {molecule.map_number_of(n) for n in molecule.atom_numbers}
    keep = (left & right) - {0}
    table = {}
    for molecule in reaction.molecules():
        for n in molecule.atom_numbers:
            number = molecule.map_number_of(n)
            if number in keep and number not in table:
                table[number] = len(table) + 1
    for molecule in reaction.molecules():
        writes = [(n, table.get(molecule.map_number_of(n), 0)) for n in molecule.atom_numbers]
        with molecule.edit():
            for n, number in writes:
                molecule.set_map_number(n, number)


def _number_inputs(inputs) -> None:
    """Number every input atom 1..N from one counter, in input order.

    Unconditional, unlike `ReactionContainer.reset_mapping()`: a number the record arrived with is the
    claim under test, not an input to reconstruction.
    """
    number = 0
    for molecule in inputs:
        for n in molecule.atom_numbers:
            number += 1
            molecule.set_map_number(n, number)


def _grossly_unbalanced(recorded, inputs, options) -> bool:
    """True when the recorded product is too much bigger than everything that went in to be explained.

    Off entirely at `max_size_ratio <= 0`, never consulted below `min_filter_size` heavy atoms, otherwise
    a ratio against the input total -- every input, agents included.  A bound on the SEARCH and not a
    judgement about the chemistry: the record still gets its ordinary unexplained line.
    """
    if options.max_size_ratio <= 0:
        return False
    product_size = len(recorded.atom_numbers)
    if product_size < options.min_filter_size:
        return False
    total = sum(len(molecule.atom_numbers) for molecule in inputs)
    return product_size >= options.max_size_ratio * total


def _refuse(log, rule, message) -> None:
    _record(log, rule, REFUSED, message)


def _record(log, rule, severity, message) -> None:
    """`log` is always `reaction.log`; there is no `log=` to pass and nothing to switch off."""
    log.append(LogRecord(rule, (), message, severity, 'reconstruct'))
