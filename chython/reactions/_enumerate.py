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
"""`mol.react()`, `mol @ mol`, `mol.functional_groups()`, `mol.protective_groups()`, `mol.deprotect()`.

This file decides only WHICH templates to try; applying one is the core's `template(*molecules)`.  Every
input goes to the matcher at once, so argument order cannot matter and a mixture handed in as one
container works; an outcome is then kept only when it touched every input.  Deprotection is the second
half, and the primitive it adds is the CLAIM -- see `_claims`.
"""
from collections import Counter
from collections.abc import Iterable, Iterator, Mapping, Sequence
from itertools import chain, combinations
from typing import NamedTuple
from ._tables import (PROTECTS_ELEMENTS, ProtectiveGroup, ReactionRule,
                      functional_rules as _group_table, protective_rules, reaction_rules)
from ..core import MoleculeContainer, ReactionContainer, ReactionTemplate


__all__ = ['EnumeratedDeprotection', 'EnumeratedReaction', 'GroupHit', 'deprotect',
           'functional_group_hits', 'functional_groups', 'protective_group_hits', 'protective_groups',
           'react']


class EnumeratedReaction(NamedTuple):
    """One enumerated outcome: the rule's name, the reaction it produced, and the row it came from.

    `name` is the chemistry (`'suzuki'`), shared by every row spelling a variant of it, and is what
    `reaction=` selects on.  `rule_id` is the row (`'reactions:8'`): only the id says which spelling
    earned the hit, which is what a coverage measurement counts on.
    """
    name: str
    reaction: ReactionContainer
    rule_id: str = ''


class GroupHit(NamedTuple):
    """One group a molecule carries: the row it came from, the row's name, and how many times.

    The id is here rather than reachable through a second `functional_rules()` lookup because an id is
    what a consumer persists.  A name is the API name and an id is the row: both are stable, and only the
    id is unique across the four corpora.
    """
    id: str
    name: str
    count: int


def functional_group_hits(molecule: MoleculeContainer) -> tuple[GroupHit, ...]:
    """Every group `functional.tsv` names and this molecule carries, in the table's order.

    The count is distinct SITES -- how many sets of atoms the pattern covers -- so a symmetric diester
    reports two acid groups and not four.  Counting the mappings instead over-reports by the pattern's
    OWN symmetry, which is a property of how the row is spelled rather than of the molecule:
    `trifluoromethyl` writes its three fluorines out, so one CF3 admits 3! = 6 mappings of the same four
    atoms.  Two vicinal diols in glycerol and two acids in terephthalic acid are genuinely two sites and
    still report 2.

    `functional_groups()` is this folded to `{name: count}`.
    """
    found = []
    for name, group in _group_table().items():
        sites = {frozenset(mapping.values()) for mapping in group.query.get_mapping(molecule)}
        if sites:
            found.append(GroupHit(group.id, name, len(sites)))
    return tuple(found)


def functional_groups(molecule: MoleculeContainer) -> dict[str, int]:
    """`{name: count}` for the groups `functional.tsv` names and this molecule carries.

    An absent group is absent from the dict rather than present with a zero, which makes
    `name in mol.functional_groups()` the presence test.  `functional_group_hits()` is the same answer
    with each row's id beside its count.
    """
    return {hit.name: hit.count for hit in functional_group_hits(molecule)}


def _selected(rules: Mapping[str, tuple[ReactionRule, ...]],
              reaction: str | None) -> Iterator[ReactionRule]:
    """The rows a call may use: the family with the right name, or all of them.

    A LOOKUP AND NOT A SCAN, which is what keying the corpus on the name buys: `reaction=` names a row
    FAMILY, so selection is one `dict` hit and the known-name set the failure message needs is the keys.

    An unknown `reaction=` raises rather than yielding nothing -- the one refusal in this file, since
    "no such reaction" and "that reaction does not apply here" are the same empty generator to a caller
    and only one of them is their bug.  The message lists every name there is.
    """
    if reaction is not None:
        if reaction not in rules:
            raise ValueError('unknown reaction %r; the corpus names %d reactions: %s'
                             % (reaction, len(rules), ', '.join(sorted(rules))))
        return iter(rules[reaction])
    return chain.from_iterable(rules.values())


def _run(molecules: Sequence[MoleculeContainer], rules: Mapping[str, tuple[ReactionRule, ...]],
         reaction: str | None = None) -> Iterator[EnumeratedReaction]:
    """The enumeration itself.  One presence scan per input, then every row that fits.

    Separate from `react()` only so a test can hand it a rule fixture the corpus has no row for.

    Two filters.  The functional-group multiset over the union of the inputs is a PREFILTER and nothing
    more -- the groups are separate subgraph matches, so all of them present is not a promise that one
    match contains them all.  On the way out, every input must be touched
    (`len(rxn.reactants) == inputs`), so a three-molecule question is never answered by a row that
    ignores one; an untouched COMPONENT of a touched input is different and survives, which is the salt
    rule.  Nothing between them knows how many molecules a row "expects".
    """
    available = Counter()
    for molecule in molecules:
        available.update(functional_groups(molecule))
    inputs = len(molecules)

    for rule in _selected(rules, reaction):
        if Counter(rule.groups) - available:
            continue
        for template in rule.templates:
            for rxn in template(*molecules):
                if len(rxn.reactants) == inputs:
                    yield EnumeratedReaction(rule.name, rxn, rule.id)


def react(molecule: MoleculeContainer, others=(), reaction: str | None = None
          ) -> Iterator[EnumeratedReaction]:
    """Enumerate reactions of `molecule` with `others`.  The one enumeration entry point.

    The argument order is not the slot order: a row's slots are chemical roles, so `amine.react(acid)`
    and `acid.react(amine)` are the same question.

    No partner is a complete question -- `mol.react()` is every single-molecule row the corpus has.  A
    call WITH partners never yields a one-slot row, since such a row can touch only one input and every
    input must be touched.

    `reaction` selects rows by chemical name; an unknown one raises.
    """
    return _run((molecule, *others), reaction_rules(), reaction)


# --- deprotection ---------------------------------------------------------------------------------

class EnumeratedDeprotection(NamedTuple):
    """One enumerated deprotection: which groups came off, the reaction, and the rows behind it.

    `names` and `rule_ids` are tuples because one outcome applies a whole set of rules, in the order
    they fired (most-specific-first).  `reaction` has the untouched molecule as its single reactant and
    the stripped one as its products, mapped 1-1 from 1 with the cleaved group at 0; the caller's
    molecule is never mutated.
    """
    names: tuple[str, ...]
    reaction: ReactionContainer
    rule_ids: tuple[str, ...]


class _Claim(NamedTuple):
    """One rule matched at one site that nothing more specific had already taken.

    THE CLAIM IS OVER THE ATOMS THE RULE DELETES, not over its whole match: the atoms it merely reads --
    above all the one it reveals -- are shared context.  That is what makes an N,N-di-Boc amine two Boc
    groups, since the two matches share that nitrogen and nothing else.
    """
    rule: ProtectiveGroup
    atoms: frozenset[int]
    deleted: frozenset[int]


def _claims(molecule: MoleculeContainer, rules: Iterable[ProtectiveGroup]) -> list[_Claim]:
    """Every protecting group in `molecule`, one claim per site, most specific first.

    `rules` arrives sorted by reactant atom count descending (`protective_rules().values()`), and the
    ORDER IS LOAD-BEARING: a site is claimed by the first rule to reach it, and every later match whose
    DELETED atoms meet a standing claim is refused.  That refusal is also the automorphism filter -- a
    tert-butyl's six mappings all delete the same four atoms, so one claim comes out -- and it is what
    keeps `hydroxyl_tbu` out of a Boc.

    A revealed atom must keep a substituent.  Each row guards its own site with a degree primitive
    (`amine_boc` demands `[N;D2,D3]`), but two rows together can consume every neighbour of the atom
    they reveal, and then the product is a bare heteroatom rather than a deprotection.  Only one reading
    of a doubly substituted O can be true, so the second claim is refused; a genuinely nested group is
    seen by the next pass, once the outer one is gone.
    """
    claimed = set()
    out = []
    for rule in rules:
        deleted_query_atoms = rule.template.deleted_atoms
        for mapping in rule.template.reactants.get_mapping(molecule):
            deleted = frozenset(mapping[n] for n in deleted_query_atoms)
            if not claimed.isdisjoint(deleted):
                continue
            after = claimed | deleted
            if any(all(n in after for n in molecule.neighbors_of(atom))
                   for atom in mapping.values() if atom not in deleted):
                continue                    # nothing of the substrate would be left on a revealed atom
            claimed = after
            out.append(_Claim(rule, frozenset(mapping.values()), deleted))
    return out


def protective_group_hits(molecule: MoleculeContainer) -> tuple[GroupHit, ...]:
    """Every protecting group this molecule carries, MOST SPECIFIC FIRST -- the accessor's own order.

    The count is CLAIMS and not matches, so a Boc-protected alcohol reports one `hydroxyl_boc` and no
    `hydroxyl_tbu`, and both a bis-Boc diamine and an N,N-di-Boc amine report two.

    `protective_groups()` is this folded to `{name: count}`.
    """
    counted = Counter()
    rows = {}
    for claim in _claims(molecule, protective_rules().values()):
        counted[claim.rule.name] += 1
        rows[claim.rule.name] = claim.rule.id
    return tuple(GroupHit(rows[name], name, count) for name, count in counted.items())


def protective_groups(molecule: MoleculeContainer) -> dict[str, int]:
    """`{name: count}` for the protecting groups this molecule carries.

    An absent group is absent from the dict, so `'amine_boc' in mol.protective_groups()` is the presence
    test.  A method and not a cached property, because `deprotect()` reads it and a cache goes stale on
    the first edit.  `protective_group_hits()` is the same answer with each row's id beside its count.
    """
    return {hit.name: hit.count for hit in protective_group_hits(molecule)}


def _selected_protective(rules: tuple[ProtectiveGroup, ...], names: Sequence[str],
                         protects: Iterable[str] | None) -> tuple[ProtectiveGroup, ...]:
    """The rows a `deprotect()` call may ACT on, filtered by name and by what they reveal.

    Act on, not claim with: claims are always computed against the whole table, so asking for tert-butyl
    ethers on a molecule whose only tert-butyl is half of a Boc answers nothing rather than cleaving the
    Boc down to a carbonate.

    Both filters refuse an unknown value rather than yielding nothing.  The table's order is preserved,
    because it is the specificity order.
    """
    out = rules
    if names:
        known = {rule.name for rule in out}
        unknown = sorted(set(names) - known)
        if unknown:
            raise ValueError('unknown protecting group%s %s; protective.tsv names %d: %s'
                             % ('s' if len(unknown) > 1 else '', ', '.join(repr(n) for n in unknown),
                                len(known), ', '.join(sorted(known))))
        wanted = set(names)
        out = tuple(rule for rule in out if rule.name in wanted)
    if protects is not None:
        if isinstance(protects, str):
            protects = (protects,)
        wanted = set(protects)
        unknown = sorted(wanted - set(PROTECTS_ELEMENTS))
        if unknown:
            raise ValueError('unknown protects %s; the column takes %s'
                             % (', '.join(repr(p) for p in unknown),
                                ', '.join(sorted(PROTECTS_ELEMENTS))))
        out = tuple(rule for rule in out if wanted & set(rule.protects))
    return out


def _patch_within(template: ReactionTemplate, molecule: MoleculeContainer,
                  allowed: frozenset[int]) -> MoleculeContainer | None:
    """Apply `template` once, at a site inside `allowed`, and return the whole patched container.

    How a rule is held to its claim without the core needing a "match here" argument: stable ids survive
    a patch, so the atoms an outcome removed are a set difference, and an outcome reaching outside the
    claim was applied at a site this rule does not own.  `allowed` is the claim's DELETED set.

    The products are re-`union`ed rather than left split, because the working molecule has to stay one
    container for the next pass -- otherwise a counter-ion falls out between two deprotections.  The
    products are disjoint components of one container, so the union cannot overlap and `remap=False`
    keeps every atom number -- which is what the reaction's mapping is paired on.
    """
    before = set(molecule.atom_numbers)
    for rxn in template(molecule):
        after = set()
        for product in rxn.products:
            after |= set(product.atom_numbers)
        if before - after <= allowed:
            working = rxn.products[0]
            for product in rxn.products[1:]:
                working = working.union(product, remap=False)
            return working
    return None


def _products_key(outcome: EnumeratedDeprotection) -> frozenset:
    """What makes two outcomes the same outcome: the products, as a multiset.

    Not the names -- two site subsets of one rule are the same answer whenever a symmetry relates them.
    `__hash__` is the canonical form, so this compares structures and never SMILES strings.
    """
    return frozenset(Counter(outcome.reaction.products).items())


def _numbered(reactant: MoleculeContainer,
              products: Sequence[MoleculeContainer]) -> ReactionContainer:
    """One deprotection with the imposed 1-1 mapping: contiguous from 1 over the atoms on both sides.

    Paired by ATOM number, which survives a patch and a non-remapping `union`, so no structural search
    is needed.  An atom the strip removed is absent from the products and stays 0; deprotection creates
    none.  Ascending atom number, so the numbering is a function of the reactant and not of the order
    the claims fired.
    """
    kept = set()
    for product in products:
        kept |= set(product.atom_numbers)
    numbers = {}
    for n in reactant.atom_numbers:
        if n in kept:
            numbers[n] = len(numbers) + 1
    return ReactionContainer((_write_numbers(reactant, numbers),),
                             tuple(_write_numbers(p, numbers) for p in products))


def _write_numbers(molecule: MoleculeContainer, numbers: dict[int, int]) -> MoleculeContainer:
    """A copy of `molecule` whose map numbers are `numbers`, and 0 wherever `numbers` says nothing.

    Read before the scope opens: a container with a pending journal refuses a read.
    """
    molecule = molecule.copy()
    writes = [(n, numbers.get(n, 0)) for n in molecule.atom_numbers]
    with molecule.edit():
        for n, number in writes:
            molecule.set_map_number(n, number)
    return molecule


def _strip_sites(molecule: MoleculeContainer, chosen: Sequence[int],
                 rules: tuple[ProtectiveGroup, ...]) -> EnumeratedDeprotection | None:
    """Remove exactly the claims at the indices in `chosen`, and nothing else.

    A site is addressed by its index in the claim walk, which needs no id bookkeeping: the walk is
    deterministic and removing one site cannot renumber a lower one, so the indices are walked HIGH TO
    LOW and each is still itself when its turn comes.  The claims are recomputed from the working
    molecule each time rather than carried, because `_patch_within` may `union` and renumber.
    """
    claims = _claims(molecule, rules)
    if not all(0 <= index < len(claims) for index in chosen):
        return None
    working = molecule.copy()
    acted: list[_Claim] = []
    for index in sorted(chosen, reverse=True):
        claim = _claims(working, rules)[index]
        patched = _patch_within(claim.rule.template, working, claim.deleted)
        if patched is None:                 # the claim stood but no outcome stayed inside it
            continue
        working = patched
        acted.append(claim)
    if not acted:
        return None
    acted.reverse()                         # report in claim order, which is most-specific-first
    return EnumeratedDeprotection(tuple(claim.rule.name for claim in acted),
                                  _numbered(molecule, tuple(working.split())),
                                  tuple(claim.rule.id for claim in acted))


def _strip(molecule: MoleculeContainer, chosen: frozenset[str],
           rules: tuple[ProtectiveGroup, ...]) -> EnumeratedDeprotection | None:
    """Apply the rules named in `chosen` at every site they claim, until none claims anything.

    `rules` is the WHOLE TABLE and `chosen` the subset, and keeping them separate is the correctness
    argument for partial deprotection: claiming with the subset re-opens every shadow the subset
    excludes, and `hydroxyl_tbu` alone would claim a Boc's tert-butyl half and hand back a carbonate.

    One pass per claim, with the claims recomputed each time -- which makes sequential application its
    own overlap filter, since a site the previous pass deleted cannot be claimed again and one the
    previous cleavage exposed can be.  The atom-count bound is unreachable (every row deletes at least
    one atom, so the molecule strictly shrinks); it guards a future row that regenerates its own site.

    Returns `None` when nothing was claimed at all.
    """
    working = molecule.copy()
    names: list[str] = []
    ids: list[str] = []
    for _ in range(len(molecule)):
        claim = next((c for c in _claims(working, rules) if c.rule.name in chosen), None)
        if claim is None:
            break
        patched = _patch_within(claim.rule.template, working, claim.deleted)
        if patched is None:                 # the claim stood but no outcome stayed inside it
            break
        working = patched
        names.append(claim.rule.name)
        ids.append(claim.rule.id)
    if not names:
        return None
    return EnumeratedDeprotection(tuple(names), _numbered(molecule, tuple(working.split())),
                                  tuple(ids))


def deprotect(molecule: MoleculeContainer, names: Sequence[str] = (), *,
              protects: Iterable[str] | None = None, partial: bool = False
              ) -> Iterator[EnumeratedDeprotection]:
    """Enumerate deprotections of `molecule`.  Always an iterator; `partial` only widens it.

    `partial=False` (the default) yields at most one outcome, the full strip: every protecting group the
    molecule carries, removed.  `partial=True` yields that first and then every non-empty subset of the
    SITES, largest first, deduped by product multiset -- so a bis-Boc diamine's two symmetry-equivalent
    sites answer once.  The unit is the site and not the rule because chemistry does not go to
    completion on request: `R-N(Boc)2 -> R-NHBoc` is a real outcome.  `2^sites - 1` needs no cap, the
    subsets being generated lazily.

    `names` selects rows by name and `protects` by what they reveal (`'amine'`, `('hydroxyl', 'thiol')`).
    Neither narrows what CLAIMS a site -- claims come from the whole table on every pass -- so
    `deprotect(names=['hydroxyl_tbu'])` on a Boc-protected alcohol yields nothing rather than a
    carbonate.  An unknown name or `protects` value raises.
    """
    rules = tuple(protective_rules().values())
    actionable = {rule.name for rule in _selected_protective(rules, names, protects)}
    sites = [index for index, claim in enumerate(_claims(molecule, rules))
             if claim.rule.name in actionable]
    if not sites:
        return
    # `actionable` and not the names claimed on this pass: a group masked by the one above it is still
    # one the caller asked for, and the alternative would read as a policy about nesting.
    full = _strip(molecule, frozenset(actionable), rules)
    if full is None:
        return
    yield full
    if not partial:
        return
    seen = {_products_key(full)}
    for size in range(len(sites) - 1, 0, -1):
        for subset in combinations(sites, size):
            outcome = _strip_sites(molecule, subset, rules)
            if outcome is None:
                continue
            key = _products_key(outcome)
            if key in seen:                 # a symmetry-equivalent site, or a shadow that fell anyway
                continue
            seen.add(key)
            yield outcome
