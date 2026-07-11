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
from itertools import combinations, permutations, product as iproduct
from ..groups._oxidations import rules as oxidation_rules
from ..groups._reactions import rules as reaction_rules
from ..groups._reductions import rules as reduction_rules
from ..groups._transformations import rules as transformation_rules


class Reconstruct:
    __slots__ = ()

    def reconstruct_mapping(self, *, max_size_ratio: float = 5., min_filter_size: int = 42) -> list[str]:
        """
        Annotate reaction by trying to reconstruct the product from reactants
        using predefined reaction templates.

        Expects all structures (reactants and products) to be already
        standardized and canonicalized before calling this method.

        Tries in order:
        1. Standalone protection (kept separate; protection is normally its own reaction)
        2. Protective-group balancing: strip each reactant's reactant-only PGs once, then
        3. Standalone deprotection (a deprotected reactant equals the product)
        4. Single-molecule transforms (oxidize/reduce/transform) on deprotected reactants
        5. Multi-component reactions on deprotected reactant subsets

        Steps 3-5 run on the deprotected reactants, so deprotection composes with both transforms
        and couplings; their labels are prefixed with the ``deprotect:`` steps that produced them.

        :param max_size_ratio: reject grossly unbalanced records where the product has at least
            this many times the total reactant atom count (purifications, missing reactants, bad
            data). Lower it to reject more aggressively, raise it to be more permissive. The
            protection path is always attempted regardless (protections may graft large groups from
            an unlisted reagent). Set to a non-positive value or ``float('inf')`` to disable.
        :param min_filter_size: only apply ``max_size_ratio`` when the product has at least this
            many atoms. Small products are cheap to enumerate, so we skip the filter and go full
            compute below this threshold; the filter exists to save time on large, slow records.
            Set to 0 to always apply the ratio filter.

        Returns list of matched reaction/deprotection names.
        Empty list if no match found.
        If found, updates atom-to-atom mapping.
        """
        assert self.reactants, 'No reactants in reaction'
        assert len(self.products) == 1, 'Only single product reactions supported'
        # reset_mapping makes every atom number globally unique, so the reactant and product number
        # spaces are disjoint. The reconstruction remaps surviving reactant atoms onto their product
        # numbers, which is only safe (and correct) when the two spaces don't overlap.
        self.reset_mapping()

        reactants = self.reactants
        product = self.products[0]
        product_size = len(product)

        # Precheck: purification (product is literally one of the reactants, i.e. recrystallization
        # or workup logged as a reaction). Cheap structural equality; bail before the expensive
        # functional_groups / template enumeration below.
        if product in reactants:
            return []

        # 1. Standalone protection (separate reaction; never composed with transforms/couplings).
        # Tried before the size filter: a protection legitimately grafts a large group from an
        # unlisted reagent, so the product can dwarf the listed reactants yet still be a real match.
        for r in reactants:
            if (result := _try_protect(r, product)) is not None:
                mapping, labels = result
                r.remap(mapping)
                return labels

        # Reject grossly unbalanced records (missing reactants, bad data): a product with
        # >=max_size_ratio x the total reactant atoms cannot be reconstructed by any remaining
        # template. Skipped for small products (< min_filter_size atoms): enumeration there is cheap,
        # so full compute is preferred over the risk of rejecting a real reaction.
        total_reactant_atoms = sum(len(r) for r in reactants)
        if max_size_ratio > 0 and product_size >= min_filter_size \
                and product_size >= max_size_ratio * total_reactant_atoms:
            return []

        # 2+3. Deprotection: strip each reactant's reactant-only PGs and, the moment a deprotected
        # reactant equals the product, return without deprotecting the rest. `x` maps the surviving
        # atoms (reactant numbers) onto product numbers; the original reactant shares those numbers,
        # so we remap it directly (removed-PG atoms aren't in `x` and keep their own numbers).
        #
        # Each reactant then contributes one or two candidate forms to the chemistry phases: always
        # the raw reactant, plus its deprotected form when a deprotection happened. We run chemistry
        # on both because the listed PG is not always cleaved as a clean deprotection — it may be
        # transformed in place by the reaction, in which case the raw form is the right substrate.
        # `forms[i]` is a list of (mol, labels) for reactant i, raw first so it is preferred.
        forms = []
        for r in reactants:
            candidates = [(r, [])]
            if (result := _deprotect_excess(r, product)) is not None:
                m, labels = result
                if (mapping := m.get_fast_mapping(product)) is not None:
                    r.remap(mapping)
                    return labels
                candidates.append((m, labels))
            forms.append(candidates)

        # 4. Single-molecule transforms (largest reactant first). FG-screened, so tiny reagents cost
        # almost nothing; this also covers reagent-adding transforms the size-based subset gate skips.
        for i in sorted(range(len(forms)), key=lambda j: len(reactants[j]), reverse=True):
            for m, dep in forms[i]:
                if (result := _try_transforms(m, product)) is not None:
                    mapping, labels = result
                    reactants[i].remap(mapping)
                    return dep + labels

        # 5. Multi-component reactions over reactant subsets (size-ordered). For each subset we try
        # every combination of its members' candidate forms (raw / deprotected).
        for subset_indices in _candidate_subsets(reactants, product_size):
            if len(subset_indices) == 1:
                continue  # singles already covered by phases 3-4
            for combo in iproduct(*(forms[i] for i in subset_indices)):
                subset = [m for m, _ in combo]
                if (result := _try_multi(subset, product)) is not None:
                    mapping, labels = result
                    for i in subset_indices:
                        reactants[i].remap(mapping)
                    return [l for _, dep in combo for l in dep] + labels

        return []


def _candidate_subsets(reactants, product_size):
    """
    Yield reactant index tuples ordered by likelihood of match.

    1. Singles (sorted by atom count descending, threshold: size >= product_size * 0.5)
    2. Pairs (sorted by combined size descending, threshold: combined >= product_size * 0.7)
    3. Triples+ (sorted by combined size descending)
    """
    sizes = [len(r) for r in reactants]
    n = len(reactants)

    # Singles: largest first, threshold 50% of product
    threshold_single = product_size * 0.5
    singles = [(i,) for i in range(n) if sizes[i] >= threshold_single]
    singles.sort(key=lambda t: sizes[t[0]], reverse=True)
    yield from singles

    # Pairs: largest combined first, threshold 70% of product
    if n >= 2:
        threshold_pair = product_size * 0.7
        pairs = []
        for combo in combinations(range(n), 2):
            combined = sum(sizes[i] for i in combo)
            if combined >= threshold_pair:
                pairs.append(combo)
        pairs.sort(key=lambda t: sum(sizes[i] for i in t), reverse=True)
        yield from pairs

    # Triples
    if n >= 3:
        triples = list(combinations(range(n), 3))
        triples.sort(key=lambda t: sum(sizes[i] for i in t), reverse=True)
        yield from triples

    # Quads
    if n >= 4:
        quads = list(combinations(range(n), 4))
        quads.sort(key=lambda t: sum(sizes[i] for i in t), reverse=True)
        yield from quads


def _strip(mol, groups, start=None):
    """
    Return a copy of ``mol`` with the named protective groups removed, or None on failure.

    ``remove_protection`` canonicalizes the result (kekule/thiele/...), which can raise on a
    structure that has no valid kekule form. Reconstruction must never crash, so a failed strip is
    reported as "no match" (None) rather than propagated.

    ``start``, when given, is passed to ``remove_protection`` as the starting atom number for any
    newly added atoms (e.g. the carbonyl O restored from an acetal), avoiding collisions with
    other molecules in a multi-component reaction context.
    """
    m = mol.copy(keep_sssr=True, keep_components=True)
    try:
        cursor = start
        for p in groups:
            before = set(m._atoms)
            m.remove_protection(p, start=cursor)
            if cursor is not None:
                cursor += len(m) - len(before)
    except Exception:  # noqa: BLE001 - an unkekulizable deprotection product is simply not a match
        return None
    return m


def _deprotect_excess(r, product):
    """
    Strip the protective groups that ``r`` carries but the product does not (a deprotection step).

    Returns ``(molecule, labels)`` when a deprotection happened — ``molecule`` is the deprotected
    copy of ``r`` (keeping ``r``'s atom numbers minus the removed PG atoms, so it stays directly
    comparable to the disjoint-numbered product) and ``labels`` are the ``deprotect:<name>`` steps.
    Returns ``None`` when there is nothing to strip or the strip fails.

    PGs shared with the product, or present only in the product, are left untouched: extra
    product-side PGs mean a protection happened and are handled by ``_try_protect`` separately.
    """
    rpg = r.protective_groups
    ppg = product.protective_groups
    # only strip reactant-only PGs whose shared instances already balance; if a shared PG has
    # mismatched counts, the bookkeeping is ambiguous, so deprotect nothing and let matching fail.
    for p in rpg:
        if p in ppg and rpg[p] != ppg[p]:
            return None
    to_remove = [p for p in rpg if p not in ppg]
    if not to_remove:
        return None
    if (m := _strip(r, to_remove, max(product) + 1)) is None:
        return None
    return m, [f'deprotect:{name}' for name in to_remove]


def _try_protect(r, product):
    """
    Try standalone protection (reverse direction): the product is the reactant with PGs added, i.e.
    the product carries PGs the reactant does not. Verify by deprotecting the product back to the
    reactant.

    Returns ``(mapping, labels)`` on a match — ``mapping`` is ``{reactant_atom: product_atom}`` and
    ``labels`` are the ``protect:<name>`` steps — else None.
    """
    rpg = r.protective_groups
    ppg = product.protective_groups
    # reactant must have NO PGs that product doesn't have
    for p in rpg:
        if p not in ppg:
            return None
    # shared PGs must have same counts
    for p in ppg:
        if p in rpg and ppg[p] != rpg[p]:
            return None
    # PGs to add (present in product, absent from reactant)
    to_add = [p for p in ppg if p not in rpg]
    if not to_add:
        return None
    # deprotect a copy of product and compare against reactant. `m` keeps the product's numbering
    # minus the grafted PG atoms, so reactant -> m is {reactant_atom: product_atom}; the grafted PG
    # atoms live only on the product side and stay unmapped (keeping their own unique numbers).
    if (m := _strip(product, to_add)) is None:
        return None
    if x := r.get_fast_mapping(m):
        return x, [f'protect:{name}' for name in to_add]
    return None


def _reactor_products(reactor, *mols):
    """
    Iterate a reactor's product molecules, tolerating chemistry failures.

    Template application ends with canonicalization (kekule/thiele/...), which can raise on a
    generated structure that has no valid kekule form or otherwise breaks aromaticity. That
    exception surfaces from the reactor generator's ``__next__`` and, once raised, kills the
    generator. Reconstruction must never crash on such a template, so we stop iterating this reactor
    on any exception and move on to the next rule.
    """
    gen = iter(reactor(*mols))
    while True:
        try:
            rxn = next(gen)
        except StopIteration:
            return
        except Exception:  # noqa: BLE001 - a template yielding invalid chemistry is not a match
            return
        yield rxn.products[0]


def _try_transforms(r, product):
    """
    Unified loop over oxidation, reduction, and transformation rules.
    Reactor products are fully canonicalized via reactor pipeline.

    Returns (mapping, labels) or None. The mapping is restricted to ``r``'s atoms, so reactor-created
    atoms (numbered outside ``r``, e.g. the new Br of a bromination) are dropped and the result is
    directly remappable onto the reactant.
    """
    rfgs = r.functional_groups
    pfgs = product.functional_groups  # cached on product

    # Oxidation rules
    for name, fg_name, output_fg, reactor in oxidation_rules:
        if fg_name not in rfgs:
            continue
        if output_fg is not None and pfgs.get(output_fg, 0) - rfgs.get(output_fg, 0) != 1:
            continue
        if output_fg is None and rfgs[fg_name] - pfgs.get(fg_name, 0) != 1:
            continue
        for p in _reactor_products(reactor, r):
            if x := p.get_fast_mapping(product):
                return {n: x[n] for n in r if n in x}, [f'oxidize:{name}']

    # Reduction rules
    for name, fg_name, output_fg, reactor in reduction_rules:
        if fg_name not in rfgs:
            continue
        if output_fg is not None and pfgs.get(output_fg, 0) - rfgs.get(output_fg, 0) != 1:
            continue
        if output_fg is None and rfgs[fg_name] - pfgs.get(fg_name, 0) != 1:
            continue
        for p in _reactor_products(reactor, r):
            if x := p.get_fast_mapping(product):
                return {n: x[n] for n in r if n in x}, [f'reduce:{name}']

    # Transformation rules
    for name, fg_name, output_fg, reactor in transformation_rules:
        if fg_name not in rfgs:
            continue
        if output_fg is not None and pfgs.get(output_fg, 0) - rfgs.get(output_fg, 0) != 1:
            continue
        if output_fg is None and rfgs[fg_name] - pfgs.get(fg_name, 0) != 1:
            continue
        for p in _reactor_products(reactor, r):
            if x := p.get_fast_mapping(product):
                return {n: x[n] for n in r if n in x}, [f'transform:{name}']

    return None


def _try_multi(subset, product):
    """
    Multi-component reactions: try all reaction templates that match
    the subset length and FG requirements.

    Returns (mapping, labels) or None. The mapping is restricted to the fed reactants' atoms, so
    reactor-created atoms (numbered outside the subset) are dropped and the result is directly
    remappable onto the reactants.
    """
    fgs = [r.functional_groups for r in subset]
    n = len(subset)
    source = set().union(*subset)

    for name, fg_names, reactor in reaction_rules:
        if len(fg_names) != n:
            continue
        for perm in permutations(range(n)):
            if all(fg_names[i] in fgs[j] for i, j in enumerate(perm)):
                for p in _reactor_products(reactor, *(subset[j] for j in perm)):
                    if x := p.get_fast_mapping(product):
                        return {a: x[a] for a in source if a in x}, [f'react:{name}']
                break

    return None


__all__ = ['Reconstruct']
