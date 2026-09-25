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

THE CONSTITUTION AND THE CONFIGURATION ARE TWO QUESTIONS.  A rung reproduces the constitution or it does
not; a configuration that disagrees is a fact about the record, reported and optionally repaired, and
never a reason to hand back no mapping.  `_Options.loose` is that separation, and `_Link` is what carries
each rung's paired components out to it.
"""
from collections import Counter
from collections.abc import Callable, Container, Iterator, Mapping, Sequence
from itertools import combinations
from typing import NamedTuple

from ._enumerate import (_deprotect_toward, _fitting, _outcomes, deprotect, EnumeratedReaction,
                         functional_groups)
from ._numbering import fast_mapping
from ._tables import reaction_rules
from ..core import (INFO, LOST, LogRecord, MoleculeContainer, REFUSED, ReactionContainer, STEREO_ABS,
                    STEREO_UNSPECIFIED)


class _Options(NamedTuple):
    max_size_ratio: float
    min_filter_size: int
    loose: bool = False     # is the configuration part of the question this pass asks?


class _Link(NamedTuple):
    """One recorded component and the component a rung paired it with, plus the atom correspondence.

    `source` carries the PREDICTION -- the rebuilt product for a rung that built one, the input itself
    for the identity rungs -- and `recorded` carries the OBSERVATION.  `pairs` is `{source id: recorded
    id}`.  Both the stereo-mismatch line and the heal read these and nothing else, so the comparison is
    written once and every rung feeds it by returning links from its write.
    """
    source: MoleculeContainer
    recorded: MoleculeContainer
    pairs: dict


class _Explanation(NamedTuple):
    label: str
    write: Callable[[MoleculeContainer], tuple[bool, list[_Link]]]
    rule: str = ''          # the table-qualified row id, when a row produced this; logged, not returned


def reconstruct_mapping(reaction: ReactionContainer, *, max_size_ratio: float = 5.,
                        min_filter_size: int = 42, stereo: str = 'loose',
                        heal_stereo=False) -> tuple[str, ...]:
    """See `ReactionContainer.reconstruct_mapping`, whose body this is."""
    if stereo != 'loose' and stereo != 'strict':
        raise ValueError("stereo is 'loose' or 'strict', not %r" % (stereo,))
    if heal_stereo is not False and heal_stereo != 'product' and heal_stereo != 'reactant':
        raise ValueError("heal_stereo is False, 'product' or 'reactant', not %r" % (heal_stereo,))
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

    unbalanced = _grossly_unbalanced(recorded, inputs, _Options(max_size_ratio, min_filter_size))
    if unbalanced:
        _refuse(log, 'reconstruct:unbalanced',
                'the recorded product has %d atoms against %d in every input together: the rungs that '
                'search the corpus are skipped, since the inputs do not contain the product'
                % (len(recorded.atom_numbers),
                   sum(len(molecule.atom_numbers) for molecule in inputs)))
    # THE STRICT WALK FIRST, ALWAYS.  The corpus discriminates rows by configuration -- one row carries
    # a centre through and another turns it over -- so a record that states its product configuration is
    # entitled to the row that actually explains it.  The loose walk is what a record whose
    # configuration explains nothing falls through to, and it reaches only records that would otherwise
    # have come back unexplained.
    #
    # ONE MEMO FOR THE CALL AND BOTH WALKS SHARE IT: a rung enumerates what the inputs and the corpus
    # allow, which is not a question about strictness, so the loose walk reads the strict walk's outcomes
    # back rather than putting the same question again.
    memo = {}
    for loose in ((False, True) if stereo == 'loose' else (False,)):
        options = _Options(max_size_ratio, min_filter_size, loose)
        for phase in _PHASES:
            if unbalanced and phase not in _FILTER_EXEMPT:
                continue
            found = list(phase(recorded, inputs, options, memo))
            if not found:
                continue
            complete, links = found[0].write(recorded)
            if not complete:
                _record(log, 'reconstruct:partial', INFO,
                        'some components of the recorded product were not reproduced and are left '
                        'unnumbered')
            if found[0].rule:
                # The label is the chemistry, the id is the row; only the id says which spelling earned
                # the hit.  The APPLIED rule only -- logging one merely considered is a false
                # attribution.
                _record(log, 'reconstruct:rules', INFO, 'explained by %s' % found[0].rule)
            if len(found) > 1:
                alternatives = sorted({e.label for e in found[1:]})
                _record(log, 'reconstruct:alternatives', INFO,
                        'alternative explanations not applied: %s' % ', '.join(alternatives))
            _report_stereo(log, links)
            if heal_stereo:
                _heal(log, links, heal_stereo, reaction)
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

def _once(memo: dict, key: str, compute: Callable[[], object]):
    """`compute()`'s value, computed on the first walk of a call and read back on the second.

    WHAT A RUNG ENUMERATES IS NOT A QUESTION ABOUT STRICTNESS.  `_applications` and `deprotect` take no
    `_Options`: both walks put the same question to the corpus and get the same outcomes back, and what
    strictness decides is which of those outcomes `_reproduces` accepts.  So the walk is paid once.

    `compute` materializes -- the driver exhausts every rung with `list()`, so nothing is lost -- and the
    outcomes are held, not the generator, because a generator replays nothing.  Sound only because every
    consumer of an outcome reads it: `_reproduces` compares, `_translate` copies before it edits.
    """
    try:
        return memo[key]
    except KeyError:
        got = memo[key] = compute()
        return got


def _purification(recorded, inputs, options, memo) -> Iterator[_Explanation]:
    """The product went in and came out.  The mapping it implies is the identity.

    The one explanation whose mapping is certain, and the only label with no namespace -- no rule and no
    table produced it.  Nothing to memo: the identity test IS the strictness question.
    """
    if _number_product(recorded.copy(), inputs, options.loose)[0]:
        yield _Explanation('purification',
                           lambda target, o=options: _number_product(target, inputs, o.loose))


def _heavy(molecule: MoleculeContainer) -> dict[int, int]:
    """Element counts without hydrogen: what canonicalization cannot change."""
    return {element: n for element, n in molecule.element_counts.items() if element != 1}


def _settle(products: Sequence[MoleculeContainer], targets: list[dict[int, int]]
            ) -> list[MoleculeContainer]:
    """`products`, each one with a component the heavy-atom formula of a `targets` entry canonicalized.

    The recorded side is canonicalized; a reactor or a strip returns its patch raw -- mobile hydrogens
    where the patch left them, an unmasked azole nitrogen at `H_UNKNOWN`.  So the comparison is
    canonical against canonical.  The formula screen is there because `canonicalize()` is the cost: an
    outcome that cannot be a recorded component is never compared, so it is never settled either.
    Copies; an outcome is read by both walks.
    """
    out = []
    for product in products:
        if any(_heavy(part) in targets for part in product.split()):
            product = product.copy()
            product.canonicalize()
        out.append(product)
    return out


def _settled(recorded: MoleculeContainer, found: Iterator) -> list:
    """`(source, outcome, settled products)` for every `(source, outcome)` a rung enumerated."""
    targets = [_heavy(part) for part in recorded.split()]
    return [(source, outcome, _settle(outcome.reaction.products, targets)) for source, outcome in found]


def _translate(reaction, sources: Sequence[MoleculeContainer],
               products: Sequence[MoleculeContainer] | None = None) -> list[MoleculeContainer]:
    """The reaction's products carrying `sources`' map numbers in place of the reactor's own.

    The two schemes compose through ATOM numbers: a reactor reactant is a copy of its source at the
    same atom numbers, so reactant map -> atom number -> source map.  `sources` is positional, which
    `_run` guarantees -- an outcome names the inputs its match touched, in input order.  An atom the
    reactor left at 0, and one no source claims, stays 0.  `products` replaces the reaction's own, for a
    settled copy still wearing the reactor's numbers.
    """
    table = {}
    for reactant, source in zip(reaction.reactants, sources):
        for n in reactant.atom_numbers:
            number = reactant.map_number_of(n)
            if number:
                table[number] = source.map_number_of(n)
    out = []
    for product in (reaction.products if products is None else products):
        product = product.copy()
        writes = [(n, table.get(product.map_number_of(n), 0)) for n in product.atom_numbers]
        with product.edit():
            for n, number in writes:
                product.set_map_number(n, number)
        out.append(product)
    return out


def _react(recorded, inputs, options, memo) -> Iterator[_Explanation]:
    """A corpus row, applied to the inputs as they arrived.  The strongest evidence there is."""
    for pool, outcome, products in _once(memo, 'react',
                                         lambda: _settled(recorded, _applications(recorded, inputs))):
        if not _reproduces(recorded, products, options.loose):
            continue
        sources = [*_translate(outcome.reaction, pool, products), *inputs]
        yield _Explanation('react:%s' % outcome.name,
                           lambda target, s=sources, o=options: _number_product(target, s, o.loose),
                           outcome.rule_id)


def _deprotect(recorded, inputs, options, memo) -> Iterator[_Explanation]:
    """An input, unmasked.  The recorded product IS the recorded input minus a protecting group.

    `partial=True` because incomplete cleavage is ordinary and `R-N(Boc)2 -> R-NHBoc` is a record this
    rung must read.  Largest-first, which is the order the explanations come out in.  Only the subsets
    that can give a recorded component are stripped: `_deprotect_toward`.
    """
    parts = list(recorded.split())
    for molecule, outcome, products in _once(memo, 'deprotect',
                                             lambda: _settled(recorded, ((m, o) for m in inputs
                                                                         for o in _deprotect_toward(m, parts)))):
        if not _reproduces(recorded, products, options.loose):
            continue
        sources = [*_translate(outcome.reaction, [molecule], products), *inputs]
        yield _Explanation('deprotect:%s' % '+'.join(outcome.names),
                           lambda target, s=sources, o=options: _number_product(target, s, o.loose),
                           '+'.join(outcome.rule_ids))


def _strip(inputs) -> tuple[list[MoleculeContainer], set[int]]:
    """`inputs` with every protecting group taken off, and which members of the result a strip produced.

    ONE all-stripped pass and not every raw/stripped combination, which would be exponential in the
    number of protected inputs: a deliberate lower bound on what composition buys.  A stripped form
    keeps the stable ids it came in with, so numbering survives the composition unaided.  A stripped
    form is canonicalized, as the inputs were: the corpus is written against that form, and a strip
    leaves an unmasked azole nitrogen at `H_UNKNOWN`.
    """
    pool = []
    stripped = set()
    for molecule in inputs:
        outcome = next(deprotect(molecule), None)
        if outcome is None:
            pool.append(molecule)
        else:
            for product in _translate(outcome.reaction, [molecule]):
                product.canonicalize()
                stripped.add(len(pool))
                pool.append(product)
    return pool, stripped


def _deprotect_then_react(recorded, inputs, options, memo) -> Iterator[_Explanation]:
    """Strip what can be stripped, then let the corpus fire on what is left.

    Only the pools a stripping reached are walked -- `stripped`, handed on as `_applications`' `required`.
    A pool of unstripped members is one `_react` already put to the corpus, and this rung runs after it.
    """
    pool, stripped = _once(memo, 'strip', lambda: _strip(inputs))
    if not stripped:
        return
    for subset, outcome, products in _once(memo, 'deprotect+react',
                                           lambda: _settled(recorded,
                                                            _applications(recorded, pool, required=stripped))):
        if not _reproduces(recorded, products, options.loose):
            continue
        sources = [*_translate(outcome.reaction, subset, products), *pool, *inputs]
        yield _Explanation('deprotect+react:%s' % outcome.name,
                           lambda target, s=sources, o=options: _number_product(target, s, o.loose),
                           outcome.rule_id)


def _protect(recorded, inputs, options, memo) -> Iterator[_Explanation]:
    """The recorded product is a recorded input, masked.

    THE WEAKEST RUNG AND SO THE LAST: an amide, an ester and a carbamate are all protecting groups as
    well as products, so offered first this reads every acylation as a protection.  Only the subsets
    that can give an input component are stripped.
    """
    parts = [part for molecule in inputs for part in molecule.split()]
    for outcome in _once(memo, 'protect', lambda: list(_deprotect_toward(recorded, parts))):
        pairs = _pair_with_inputs(outcome.reaction.products, inputs, options.loose)
        if not pairs:
            continue
        yield _Explanation('protect:%s' % '+'.join(outcome.names),
                           lambda target, p=pairs, o=options: _number_from_pairs(target, p, o.loose),
                           '+'.join(outcome.rule_ids))


_PHASES: tuple[Callable[..., Iterator[_Explanation]], ...] = (_purification, _react, _deprotect,
                                                              _deprotect_then_react, _protect)
_FILTER_EXEMPT = frozenset((_purification, _protect))


# --- writing the numbers --------------------------------------------------------------------------

def _same(a: MoleculeContainer, b: MoleculeContainer, loose: bool) -> bool:
    """Are these two components the same structure?  `loose` drops the configuration from the question.

    `__eq__` is the canonical form and carries every parity.  `isomorphism()` asks the same question with
    the invariant's stereo term off and still requires the whole constitution -- element, charge,
    isotope, radical, R index, implicit hydrogen count and every bond order -- so loose relaxes the
    configuration and NOTHING else.
    """
    return a.isomorphism(b) is not None if loose else a == b


def _pairs_of(a: MoleculeContainer, b: MoleculeContainer, loose: bool) -> dict | None:
    """`{a id: b id}` when the two are the same structure under this pass's question, else None."""
    return a.isomorphism(b) if loose else fast_mapping(a, b)


def _number_product(recorded: MoleculeContainer, sources: Sequence[MoleculeContainer],
                    loose: bool) -> tuple[bool, list[_Link]]:
    """Number the components of `recorded` from the components of the numbered `sources`.

    Per connected component, which is what keeps a salt alive: a component no template touched is
    matched against the input it arrived in and numbered from there.  `split()` preserves stable ids, so
    the write goes straight onto `recorded`.  Each source component is consumed at most once.

    Returns whether every component was numbered, and the links it consumed; a partial answer is still
    written and the caller logs the shortfall.  THE ORDER OF `sources` IS THE PREDICTION ORDER: a rung
    that rebuilt the product puts the rebuilt components first, so a recorded component pairs with the
    row's own answer where there is one and with the input it arrived in where there is not.

    Writes onto the RECORDED containers, and `mapping_agrees` depends on it: that comparison excuses an
    automorphic swap by an orbit lookup keyed on the input's own stable ids, so rebuilt inputs would
    make every legitimate swap count as a disagreement.
    """
    parts = list(recorded.split())
    pool: list[MoleculeContainer] = [part for source in sources for part in source.split()]
    used = set()
    links = []
    for part in parts:
        for i, candidate in enumerate(pool):
            if i in used:
                continue
            pairs = _pairs_of(candidate, part, loose)
            if pairs is None:
                continue
            for source_id, target_id in pairs.items():
                recorded.set_map_number(target_id, candidate.map_number_of(source_id))
            used.add(i)
            links.append(_Link(candidate, part, pairs))
            break
    return len(links) == len(parts), links


def _pair_with_inputs(stripped, inputs, loose: bool):
    """`[(stripped component, the input component it equals)]`, each input component used once.

    Empty when nothing lines up, which is this rung's whole test.
    """
    pool = [part for molecule in inputs for part in molecule.split()]
    used = set()
    pairs = []
    for product in stripped:
        for part in product.split():
            for i, candidate in enumerate(pool):
                if i in used or not _same(candidate, part, loose):
                    continue
                used.add(i)
                pairs.append((part, candidate))
                break
    return pairs


def _number_from_pairs(target, pairs, loose: bool) -> tuple[bool, list[_Link]]:
    """Copy map numbers from each input component onto the stripped component it matched.

    The other direction from `_number_product`: here the recorded product has been stripped, and the
    stripped fragment's stable ids are mostly a subset of the recorded product's own, so writing at those
    ids writes onto the right atoms.  Mostly -- `deprotect()` also BUILDS, e.g. the `=O` an acetal row
    restores, and a built atom's id belongs to nothing in `target`.  Hence the guard on `ids`; such an
    atom has no source and must stay at zero.

    Complete is always True.  A protecting group's atoms are genuinely new, so zero is the right answer
    and not a shortfall; a cleaved fragment matching no input never enters `pairs` at all.

    THE LINKS RUN THE OTHER WAY ROUND TOO.  This rung reconstructs backwards, so the stripped fragment
    carries the record's own configuration and the input carries what went in -- which is the same pair
    of roles `_number_product` hands over, reached from the other end.
    """
    ids = set(target.atom_numbers)
    links = []
    for part, candidate in pairs:
        mapped = _pairs_of(part, candidate, loose)
        if mapped is None:                              # equal components, so this cannot happen
            continue
        for target_id, source_id in mapped.items():
            if target_id in ids:                        # an id the deprotection ADDED, see above
                target.set_map_number(target_id, candidate.map_number_of(source_id))
        links.append(_Link(candidate, part, {s: t for t, s in mapped.items()}))
    return True, links


# --- the configuration ----------------------------------------------------------------------------

class _Site(NamedTuple):
    """One stereo unit both sides of a link reach, with each side's reading of it in ONE order.

    `order` is the recorded unit's own direction order carried into the source's ids, and both parities
    are stated in it, which is the only way two molecules' configurations compare at all.

    `numbers` is the owners as MAP NUMBERS.  A link holds COPIES -- `split()` hands out its own
    containers and a rebuilt product never had a home -- so a heal writing onto one writes into nothing;
    the numbering this pass just established is the route from a site back to the reaction's own
    molecules, and the only name the two sides share.
    """
    link: _Link
    unit: dict                  # the recorded side's unit
    source_unit: dict
    order: tuple
    recorded_parity: int
    source_parity: int
    numbers: tuple | int | None


def _units_by_owners(molecule) -> dict:
    """`{owners: unit}` over every stereo unit.

    OWNERS AND NOT THE ANCHOR, because the key has to survive the crossing between two molecules: an
    allene's anchor is its chain midpoint and either terminal can anchor a cis/trans axis, while the
    owners are the atoms the configuration is named on whichever side names them.
    """
    return {unit['owners']: unit for unit in molecule.stereo_units()}


def _carry(owners, table):
    """`owners` in the other molecule's ids, or None where the correspondence does not reach them.

    An ascending pair stays ascending: `stereo_units()` spells an axis's owners sorted, so the mapped
    unit's key is the sorted image of the pair and not the image of the sorted pair.
    """
    if isinstance(owners, tuple):
        if owners[0] not in table or owners[1] not in table:
            return None
        return tuple(sorted((table[owners[0]], table[owners[1]])))
    return table.get(owners)


def _parity_across(molecule, anchor, order) -> int:
    """`molecule`'s parity at the unit anchored at `anchor`, read in the caller's direction `order`.

    NEVER THE TWO STORED INTEGERS SIDE BY SIDE: each molecule reads its parity against its own refs, so
    assigning one molecule's parity to another's atom inverts the centre (`set_parity`).  One order both
    molecules can state a parity in is the whole of the comparison.
    """
    try:
        return molecule.translate_stereo(anchor, order)
    except ValueError:
        # A BOND KIND THE TWO SIDES ANCHOR FROM OPPOSITE TERMINALS.  The wholesale pair exchange is a
        # legal spelling of the same order and swapping within a pair is what is forbidden (Ruling F55),
        # so the axis read from its other end is this and not a second configuration.
        return molecule.translate_stereo(anchor, order[2:] + order[:2])


def _sites(link: _Link) -> Iterator[_Site]:
    """Every stereo unit of the recorded component that the component it paired with also holds."""
    table = {target: source for source, target in link.pairs.items()}
    source_units = _units_by_owners(link.source)
    for unit in link.recorded.stereo_units():
        owners = _carry(unit['owners'], table)
        if owners is None:
            continue
        source_unit = source_units.get(owners)
        if source_unit is None:
            continue
        order = []
        for ref in unit['refs']:
            if ref is None:
                order.append(None)
            elif ref in table:
                order.append(table[ref])
            else:                                       # a direction outside the correspondence
                order = None
                break
        if order is None:
            continue
        order = tuple(order)
        yield _Site(link, unit, source_unit, order, unit['parity'],
                    _parity_across(link.source, source_unit['anchor'], order),
                    _numbers_of(link.source, owners))


def _numbers_of(molecule, owners):
    """`owners` as MAP NUMBERS, or None when one of them carries none.

    Read off the SOURCE side, which carries input numbers on every rung: an identity rung pairs with the
    input itself and a building rung with `_translate`'s copy of the row's product, which is the reactor's
    answer wearing the inputs' numbers.  Unsorted -- `_carry` sorts again in whatever frame it lands in.
    """
    if isinstance(owners, tuple):
        a = molecule.map_number_of(owners[0])
        b = molecule.map_number_of(owners[1])
        return (a, b) if a and b else None
    return molecule.map_number_of(owners) or None


def _locate(molecules, numbers) -> tuple:
    """The molecule carrying every map number in `numbers` and `{map number: its stable id}`.

    `(None, None)` when none does.  At most one can: `_number_inputs` numbers the inputs 1..N from one
    counter, and a product atom carries the number of the input atom it came from.
    """
    wanted = set(numbers) if isinstance(numbers, tuple) else {numbers}
    for molecule in molecules:
        table = {}
        for n in molecule.atom_numbers:
            number = molecule.map_number_of(n)
            if number in wanted:
                table[number] = n
        if len(table) == len(wanted):
            return molecule, table
    return None, None


def _target(site: _Site, molecules) -> tuple:
    """The container a heal writes into and ITS OWN unit at the site, or `(None, None)`.

    Its own, because a container's unit is the only thing its `set_parity` and `set_stereo_group` accept:
    the anchor a link's copy names can be the other terminal of the same axis, and an id the copy holds
    can be one a deprotection built and no molecule of the record has.
    """
    if site.numbers is None:
        return None, None
    molecule, table = _locate(molecules, site.numbers)
    if molecule is None:
        return None, None
    return molecule, _units_by_owners(molecule).get(_carry(site.numbers, table))


def _order_in(site: _Site, molecule) -> tuple | None:
    """`site.order` -- the recorded unit's direction order, in the SOURCE's ids -- in `molecule`'s ids.

    Through the map numbers again, and None where one of them does not reach `molecule`: a ref the row
    replaced belongs to another input, or to nothing that went in.
    """
    numbers = {}
    for n in molecule.atom_numbers:
        number = molecule.map_number_of(n)
        if number:
            numbers[number] = n
    out = []
    for ref in site.order:
        if ref is None:
            out.append(None)
            continue
        number = site.link.source.map_number_of(ref)
        if number not in numbers:
            return None
        out.append(numbers[number])
    return tuple(out)


def _state_parity(molecule, anchor, parity, order) -> bool:
    """State at `anchor` a `parity` read in `order`.  False when `order` is not that unit's own frame.

    False is the whole reaction-centre test for `_heal_input`: an order that is not a permutation of the
    unit's refs names a neighbour the unit does not have, which is a bond the row replaced.
    """
    try:
        molecule.set_parity(anchor, parity, order=order)
    except ValueError:
        try:
            # The axis read from its other end -- a wholesale pair exchange, Ruling F55, see
            # `_parity_across`.
            molecule.set_parity(anchor, parity, order=order[2:] + order[:2])
        except ValueError:
            return False
    return True


def _collection(molecule, owners) -> int:
    """The collection kind stated at `owners`, ABS and unstated read as the one statement they are.

    A configured atom in no collection already means absolutely configured -- which is why the writer
    emits `a:` bare -- so only OR and AND are a second claim.
    """
    kind = molecule.stereo_group_of(owners)[0]
    return STEREO_UNSPECIFIED if kind == STEREO_ABS else kind


def _report_stereo(log, links) -> None:
    """One line naming every site where the record's configuration is not the one the rung predicts.

    INFO, not LOST: the mapping is complete and correct, and the disagreement is a fact about the record.
    THREE WAYS TO DISAGREE and all three are reported, because each one is a record the strict walk
    missed: the two sides state different configurations, the record states none where the explanation
    predicts one, or the record states one the explanation does not predict -- that last being a row
    silent at the centre it wrote.  A collection can differ on a strict hit too, a collection not being
    part of the canonical form.
    """
    said = {}
    for link in links:
        for site in _sites(link):
            owners = site.unit['owners']
            if site.recorded_parity and site.source_parity:
                if site.recorded_parity != site.source_parity:
                    said.setdefault('is not the one predicted at', []).append(owners)
            elif site.source_parity:
                said.setdefault('is unstated where one is predicted at', []).append(owners)
            elif site.recorded_parity:
                said.setdefault('is stated where none is predicted at', []).append(owners)
            if _collection(link.recorded, owners) != _collection(link.source,
                                                                 site.source_unit['owners']):
                said.setdefault('is in a different collection at', []).append(owners)
    if not said:
        return
    _record(log, 'reconstruct:stereo-mismatch', INFO,
            'the record reproduces constitutionally and its configuration %s'
            % '; '.join('%s %s' % (verb, _owners_text(owners)) for verb, owners in said.items()),
            _anchors([o for owners in said.values() for o in owners]))


def _heal(log, links, source: str, reaction) -> None:
    """Write one side's configuration onto the other, through what the explanation predicts.

    `source='product'` takes the record's own product as the truth and repairs the inputs; `'reactant'`
    takes the inputs plus the row and repairs the product.  Either way the prediction mediates, which is
    what makes the reaction type count: the rung's paired component was built from the inputs THROUGH the
    row, so it already states whatever course the row states.

    Every write lands on a molecule of `reaction`, reached by map number -- see `_Site.numbers`.  Each
    line is stage `'heal'`, so the repair report is one filter away from the mapping's own log.
    """
    inputs = _inputs(reaction)
    product = reaction.products[0]
    for link in links:
        for site in _sites(link):
            if source == 'product':
                _heal_input(log, site, inputs)
            else:
                _heal_product(log, site, product)


def _heal_input(log, site: _Site, inputs) -> None:
    """The record's product decides; the input's own unit takes what follows from the record.

    TWO MECHANISMS, and which one applies is what "the reaction type on the reaction centre" comes to:

    * The input states a configuration and the prediction contradicts the record -- TURN THE INPUT OVER.
      The prediction's integer cannot simply be copied: the input's unit is read against the input's own
      neighbours, which the row changed.  What does cross is the DISAGREEMENT.  Every unit chython models
      is two-state, so the row is a bijection on those two states, and a prediction the record
      contradicts means the input holds the state's other member.  A row that carries a centre through
      and one that turns it over are repaired by this one line.
    * The input is flat and so is the prediction -- STATE THE RECORD'S OWN CONFIGURATION, in the record's
      direction order carried onto the input's atoms.  Sound exactly when that carrying succeeds: a
      stereocentre only turns over when a bond to it is broken, so an order that is still a permutation
      of the input unit's refs is a centre the row left alone, and the record's configuration is the
      input's.  A row that DID replace a bond there fails the carry and is refused.
    """
    if not site.recorded_parity:
        return                                          # the record states nothing to take
    molecule, unit = _target(site, inputs)
    if unit is None:
        return
    if site.recorded_parity == site.source_parity:
        _heal_groups(log, site, molecule, unit, from_source=False)
        return
    if not unit['stereogenic']:
        _record(log, 'heal:not-stereogenic', INFO,
                'the input unit at %s holds no second configuration, so the record\'s disagreement '
                'cannot be about it' % _owners_text([unit['owners']]),
                _anchors([unit['owners']]), 'heal')
        return
    if site.source_parity and unit['parity']:
        molecule.set_parity(unit['anchor'], 3 - unit['parity'])
        _record(log, 'heal:parity', INFO,
                'the record configures %s against what the explanation predicts from the input, so the '
                'input unit was turned over -- every unit is two-state, so the other state is the one '
                'the record implies' % _owners_text([unit['owners']]),
                _anchors([unit['owners']]), 'heal')
        _heal_groups(log, site, molecule, unit, from_source=False)
        return
    order = None if site.source_parity or unit['parity'] else _order_in(site, molecule)
    if order is not None and _state_parity(molecule, unit['anchor'], site.recorded_parity, order):
        _record(log, 'heal:parity', INFO,
                'the record configures %s where nothing went in configured, and the row leaves that '
                'centre\'s bonds alone, so the record\'s own configuration is the input\'s'
                % _owners_text([unit['owners']]), _anchors([unit['owners']]), 'heal')
        _heal_groups(log, site, molecule, unit, from_source=False)
        return
    # NO ROUTE FROM THE RECORD BACK TO AN INPUT UNIT: the row states the course outright, or is silent at
    # a centre the input configures -- the common one, and this line is a survey of which rows do not
    # state their stereochemical course -- or nothing went in configured at a centre the row rebuilt, and
    # the course is the fact that is missing.
    _record(log, 'heal:course-unstated', INFO,
            'the record configures %s and %s, so no input unit follows from the record'
            % (_owners_text([unit['owners']]),
               'the explanation predicts a configuration there that no input unit decides'
               if site.source_parity else
               'the explanation predicts none there' if unit['parity'] else
               'the row rebuilt that centre without stating a course'),
            _anchors([unit['owners']]), 'heal')


def _heal_product(log, site: _Site, product) -> None:
    """The inputs and the row decide; the recorded product takes what the prediction states.

    The direction that STATES a configuration rather than turning one over, because the record that
    dropped its product stereo is the ordinary one and a turn-over has nothing to turn.  Stated with
    `order=`, which is what makes it safe: `source_parity` is read in the recorded unit's own direction
    order, and re-basing it onto the target's refs is the difference between the configuration and its
    mirror.
    """
    if not site.source_parity:
        return                                          # the prediction states nothing to take
    molecule, unit = _target(site, [product])
    if unit is None:
        return
    if site.source_parity != site.recorded_parity:
        if not unit['stereogenic']:
            _record(log, 'heal:not-stereogenic', INFO,
                    'the recorded product holds no second configuration at %s'
                    % _owners_text([unit['owners']]), _anchors([unit['owners']]), 'heal')
            return
        if not _state_parity(molecule, unit['anchor'], site.source_parity, site.unit['refs']):
            # THE TWO SIDES READ THE SITE AGAINST DIFFERENT NEIGHBOURS, which the protect rung can
            # produce: its recorded side is the record STRIPPED, so a ref there can be an atom the
            # deprotection built and the record never had.  A parity that cannot be re-based is not one.
            _record(log, 'heal:frame-mismatch', INFO,
                    'the configuration predicted at %s is read against neighbours the recorded product '
                    'does not have, so it does not cross' % _owners_text([unit['owners']]),
                    _anchors([unit['owners']]), 'heal')
            return
        _record(log, 'heal:parity', INFO,
                'the recorded product at %s took the configuration the explanation predicts from the '
                'inputs' % _owners_text([unit['owners']]), _anchors([unit['owners']]), 'heal')
    _heal_groups(log, site, molecule, unit, from_source=True)


def _heal_groups(log, site, target, unit, from_source) -> None:
    """Carry the collection that goes with a configuration, or say why it stayed behind.

    A FRESH ID, ALWAYS.  Two records numbering their `&1` differently must not have their racemates
    merged, so the kind crosses and the number does not.  And the collection goes with the parity: a
    group is a statement about a configured unit, so one is never written where no parity is established.
    """
    if from_source:
        giver, give_owners = site.link.source, site.source_unit['owners']
    else:
        giver, give_owners = site.link.recorded, site.unit['owners']
    owners = unit['owners']
    kind = _collection(giver, give_owners)
    if kind == _collection(target, owners):
        return
    if not target.parity_of(unit['anchor']):
        _record(log, 'heal:group-not-anchored', INFO,
                'a collection was stated at %s where no configuration is established, and a promotion '
                'never invents a parity' % _owners_text([owners]), _anchors([owners]), 'heal')
        return
    if kind == STEREO_UNSPECIFIED:
        target.set_stereo_group(owners, STEREO_UNSPECIFIED)
        message = ('the collection at %s was dropped: the side the caller named states none there, and a '
                   'configured unit in no collection is already an absolute one'
                   % _owners_text([owners]))
    else:
        target.set_stereo_group(owners, kind, _fresh_group(target))
        message = ('the collection at %s crossed with the configuration, under a fresh id'
                   % _owners_text([owners]))
    _record(log, 'heal:group', INFO, message, _anchors([owners]), 'heal')


def _fresh_group(molecule) -> int:
    """The lowest collection id this molecule does not use."""
    used = {group for _, group in molecule.stereo_groups()}
    number = 1
    while number in used:
        number += 1
    return number


def _owners_text(owners) -> str:
    """`owners` as prose: an atom is its number and an axis is its two, the list ascending."""
    return ', '.join('%d-%d' % o if isinstance(o, tuple) else '%d' % o
                     for o in sorted(owners, key=lambda o: o if isinstance(o, tuple) else (o, o)))


def _anchors(owners) -> tuple:
    """The atoms a log line's `atoms` field names: every owner, an axis contributing both."""
    out = []
    for o in owners:
        out.extend(o) if isinstance(o, tuple) else out.append(o)
    return tuple(sorted(set(out)))


# --- enumeration ----------------------------------------------------------------------------------

def _supplies(supply: Counter, targets: Sequence[Mapping[str, int]], rule) -> bool:
    """Could this pool, plus whatever the row creates, hold the atoms of SOME recorded component?

    A NECESSARY CONDITION AND NOT A MATCH.  `_reproduces` asks whether a component of the record IS a
    component of what the row built, and a built component's atoms are the pool's, less what the row
    deleted, plus what it created -- so a record component wanting more of an element than the pool
    carries is reachable only through a created atom, and a row creates a fixed number of those.  One
    component sufficing is enough, the record being matched per component; deletions are not subtracted,
    which only makes the bound weaker and never wrong.

    THE SCREEN THE RECORDED PRODUCT AFFORDS AND THE GROUP PREFILTER CANNOT.  That prefilter reads the
    inputs alone, so it passes every row whose groups are in the pot however little the pot could build
    -- and the rows it passes are the ones whose cost is the isomorphism search.
    """
    created = max(len(template.created_atoms) for template in rule.templates)
    for counts in targets:
        deficit = 0
        for element, n in counts.items():
            short = n - supply.get(element, 0)
            if short > 0:
                deficit += short
                if deficit > created:
                    break
        else:
            return True
    return False


def _applications(recorded: MoleculeContainer, inputs: Sequence[MoleculeContainer],
                  rules=None, required: Container[int] | None = None
                  ) -> Iterator[tuple[list[MoleculeContainer], 'EnumeratedReaction']]:
    """Every way a corpus row applies to a SUBSET of `inputs`, with the subset it applied to.

    Subsets and not the whole list, because `_run` requires every molecule it is handed to be touched
    while a recorded record files its base, solvent and catalyst among the inputs.  Bounded by the widest
    row's slot count, so the count is a polynomial in the number of inputs and not a power set -- but
    `C(n, 1..4)` is still `n**4/24`, which is why the walk is over `reachable` and not over every input.

    THE INPUTS A ROW COULD REACH, AND NOT ALL OF THEM.  An input carrying none of the groups any row
    names a slot for cannot be one of an outcome's reactants, and `_run` keeps an outcome only when every
    molecule it was handed was touched -- so every pool containing such an input is empty by
    construction.  Enumerating them anyway lets a record's own solvents dominate its cost: a record files
    its base and solvent among its reagents, and they are typically 2 of the 4 a median record carries.

    Each input is scanned for its groups ONCE and the result handed to `_fitting`, rather than rescanned
    per pool: a molecule's groups do not depend on what it is enumerated beside.  `recorded` is here for
    `_supplies`, which is the one filter in the walk that reads the answer rather than the inputs.

    `required` is a set of indices a pool must draw at least one member from, and the caller's way of
    saying that the pools it omits were already walked: `_deprotect_then_react` runs after `_react` in the
    same walk, so a pool of its stripped list holding no stripped molecule IS a pool `_react` offered the
    corpus and `_react` not having returned is that pool's verdict.

    The subset comes back because the reactor's products carry the reactor's own numbering: `_translate`
    needs the inputs, positionally, to cross back to theirs.

    TODO: the bare `except` is defensive, not load-bearing -- applying a row has not been observed to
    raise on the current corpus, but a row that consistently fails is invisible here.  Narrow it to the
    observed type once there is one.
    """
    if rules is None:
        rules = reaction_rules()
    widest = 0
    slots = set()
    for family in rules.values():
        for rule in family:
            if len(rule.groups) > widest:
                widest = len(rule.groups)
            slots.update(rule.groups)
    carried = [functional_groups(molecule) for molecule in inputs]
    reachable = [i for i, groups in enumerate(carried) if not slots.isdisjoint(groups)]
    targets = [part.element_counts for part in recorded.split()]
    for size in range(1, min(widest, len(reachable)) + 1):
        for subset in combinations(reachable, size):
            if required is not None and not any(i in required for i in subset):
                continue
            pool = [inputs[i] for i in subset]
            supply = Counter()
            for molecule in pool:
                supply.update(molecule.element_counts)
            try:
                for rule in _fitting(pool, rules, carried=[carried[i] for i in subset]):
                    if not _supplies(supply, targets, rule):
                        continue
                    for outcome in _outcomes(rule, pool):
                        yield pool, outcome
            except Exception:
                continue


def _reproduces(recorded: MoleculeContainer, products: Sequence[MoleculeContainer],
                loose: bool) -> bool:
    """True when a connected component of `recorded` is a component of one of `products`.

    Per component, because a template answers the reaction centre while the record carries the salt too;
    and by container equality, never by SMILES -- `__eq__` is the canonical form.  The loose pass asks
    `_same` instead, which costs the set: components per record are two or three, so the pairwise scan
    is the same work.
    """
    parts = list(recorded.split())
    if not loose:
        want = set(parts)
        return any(part in want for product in products for part in product.split())
    return any(_same(candidate, part, True)
               for product in products for candidate in product.split() for part in parts)


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


def _record(log, rule, severity, message, atoms=(), stage='reconstruct') -> None:
    """`log` is always `reaction.log`; there is no `log=` to pass and nothing to switch off.

    `atoms` are stable ids in the molecule the line is about, which for a heal is the molecule it wrote.
    """
    log.append(LogRecord(rule, atoms, message, severity, stage))
