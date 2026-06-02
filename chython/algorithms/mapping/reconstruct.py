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
from itertools import combinations, permutations
from ..groups._oxidations import rules as oxidation_rules
from ..groups._reactions import rules as reaction_rules
from ..groups._reductions import rules as reduction_rules
from ..groups._transformations import rules as transformation_rules


class Reconstruct:
    __slots__ = ()

    def reconstruct_mapping(self) -> list[str]:
        """
        Annotate reaction by trying to reconstruct the product from reactants
        using predefined reaction templates.

        Tries in order:
        1. Standalone deprotection
        2. Standalone protection (reverse)
        3. Single-molecule transforms (oxidize/reduce/transform)
        4. Deprotection + transform composition
        5. Multi-component reactions (subset-based)

        Returns list of matched reaction/deprotection names.
        Empty list if no match found.
        If found, updates atom-to-atom mapping.
        """
        assert self.reactants, 'No reactants in reaction'
        assert len(self.products) == 1, 'Only single product reactions supported'
        self.reset_mapping()

        # Prepare clean copies for FG detection and reactor input
        reactants = [_prepare(m) for m in self.reactants]

        # Prepare product for comparison (canonicalize ensures canonical form)
        product = _prepare(self.products[0])
        pfgs = product.functional_groups
        product_size = len(product)

        # Iterate over candidate subsets by atom count
        for subset_indices in _candidate_subsets(reactants, product_size):
            if len(subset_indices) == 1:
                r = reactants[subset_indices[0]]
                result = _try_single(r, product, pfgs)
            else:
                subset = [reactants[i] for i in subset_indices]
                result = _try_multi(subset, product)
            if result:
                mapping, labels = result
                _safe_remap(self.products[0], mapping)
                return labels

        return []


def _safe_remap(mol, mapping):
    """Remap mol atoms using mapping, handling potential overlaps via temp numbers."""
    try:
        mol.remap(mapping)
    except ValueError:
        # Overlap: remap all atoms to safe temp space first
        existing = list(mol._atoms)
        targets = set(mapping.values())
        base = max(max(existing), max(targets), max(mapping)) + 1
        # Step 1: all atoms to temp
        temp = {n: base + i for i, n in enumerate(existing)}
        mol.remap(temp)
        # Step 2: temp to final (mapped atoms get target, unmapped get unique new numbers)
        used = set(mapping.values())
        counter = base + len(existing)
        final = {}
        for i, n in enumerate(existing):
            if n in mapping:
                final[base + i] = mapping[n]
            else:
                # Pick a number not in use
                while counter in used:
                    counter += 1
                final[base + i] = counter
                used.add(counter)
                counter += 1
        mol.remap(final)


def _prepare(mol):
    """Clean copy for FG detection and reactor input. Applies full canonicalization."""
    c = mol.copy(keep_sssr=True, keep_components=True)
    c.clean_stereo()
    c.clean_isotopes()
    c.canonicalize()
    return c


def _match(generated, product):
    """Compare generated product to expected product using fast mapping."""
    return product.get_fast_mapping(generated)


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


def _try_single(r, product, pfgs):
    """
    Try all single-reactant paths in order:
    1. Standalone deprotection
    2. Standalone protection (reverse)
    3. Single-molecule transforms (unified oxidize/reduce/transform)
    4. Deprotection + transform composition

    Returns (mapping, labels) or None.
    """
    # 1. Standalone deprotection
    result = _try_deprotect(r, product)
    if result:
        return result

    # 2. Standalone protection (reverse: product is deprotected form of reactant)
    result = _try_protect(r, product)
    if result:
        return result

    # 3. Single-molecule transforms with FG screening
    result = _try_transforms(r, product, pfgs)
    if result:
        return result

    # 4. Deprotection + transform composition
    result = _try_deprotect_then_transform(r, product, pfgs)
    if result:
        return result

    return None


def _try_deprotect(r, product):
    """
    Try standalone deprotection: reactant has PGs, product has none of those.
    Returns (mapping, labels) or None.
    """
    rpg = r.protective_groups
    ppg = product.protective_groups
    # product must have NO PGs that reactant doesn't have
    for p in ppg:
        if p not in rpg:
            return None
    # shared PGs must have same counts
    for p in rpg:
        if p in ppg and rpg[p] != ppg[p]:
            return None
    # PGs to remove: present in reactant, absent from product
    to_remove = [p for p in rpg if p not in ppg]
    if not to_remove:
        return None
    # remove all instances of each PG type
    m = r.copy(keep_sssr=True, keep_components=True)
    for p in to_remove:
        m.remove_protection(p)
    if x := product.get_fast_mapping(m):
        return x, [f'deprotect:{name}' for name in to_remove]
    return None


def _try_protect(r, product):
    """
    Try standalone protection (reverse direction):
    product is reactant with PGs added.
    Equivalent to: product has PGs that reactant doesn't.
    We check if deprotecting product yields reactant.

    Returns (mapping, labels) or None.
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
    # deprotect a copy of product and compare against reactant
    m = product.copy(keep_sssr=True, keep_components=True)
    for p in to_add:
        m.remove_protection(p)
    if x := r.get_fast_mapping(m):
        return x, [f'protect:{name}' for name in to_add]
    return None


def _try_transforms(r, product, pfgs):
    """
    Unified loop over oxidation, reduction, and transformation rules.
    Reactor products are fully canonicalized via reactor pipeline.

    Returns (mapping, labels) or None.
    """
    rfgs = r.functional_groups

    # Oxidation rules
    for name, fg_name, output_fg, reactor in oxidation_rules:
        if fg_name not in rfgs:
            continue
        if output_fg is not None and pfgs.get(output_fg, 0) - rfgs.get(output_fg, 0) != 1:
            continue
        if output_fg is None and rfgs[fg_name] - pfgs.get(fg_name, 0) != 1:
            continue
        for rxn in reactor(r):
            p = rxn.products[0]
            if x := _match(p, product):
                return x, [f'oxidize:{name}']

    # Reduction rules
    for name, fg_name, output_fg, reactor in reduction_rules:
        if fg_name not in rfgs:
            continue
        if output_fg is not None and pfgs.get(output_fg, 0) - rfgs.get(output_fg, 0) != 1:
            continue
        if output_fg is None and rfgs[fg_name] - pfgs.get(fg_name, 0) != 1:
            continue
        for rxn in reactor(r):
            p = rxn.products[0]
            if x := _match(p, product):
                return x, [f'reduce:{name}']

    # Transformation rules
    for name, fg_name, output_fg, reactor in transformation_rules:
        if fg_name not in rfgs:
            continue
        if output_fg is not None and pfgs.get(output_fg, 0) - rfgs.get(output_fg, 0) != 1:
            continue
        if output_fg is None and rfgs[fg_name] - pfgs.get(fg_name, 0) != 1:
            continue
        for rxn in reactor(r):
            p = rxn.products[0]
            if x := _match(p, product):
                return x, [f'transform:{name}']

    return None


def _try_deprotect_then_transform(r, product, pfgs):
    """
    Composition: deprotect first, then apply single-molecule transforms.
    Returns (mapping, labels) or None.
    """
    rpg = r.protective_groups
    ppg = product.protective_groups
    # Validate PG delta (same logic as standalone deprotect)
    for p in ppg:
        if p not in rpg:
            return None
    for p in rpg:
        if p in ppg and rpg[p] != ppg[p]:
            return None
    to_remove = [p for p in rpg if p not in ppg]
    if not to_remove:
        return None
    # Deprotect on a raw copy (preserves atom numbers for transform)
    deprot = r.copy(keep_sssr=True, keep_components=True)
    for p in to_remove:
        deprot.remove_protection(p)
    # Use deprotected molecule's FGs for screening (raw, unpatched)
    result = _try_transforms(deprot, product, pfgs)
    if result:
        mapping, transform_labels = result
        labels = [f'deprotect:{n}' for n in to_remove] + transform_labels
        return mapping, labels
    return None


def _try_multi(subset, product):
    """
    Multi-component reactions: try all reaction templates that match
    the subset length and FG requirements.

    Returns (mapping, labels) or None.
    """
    fgs = [r.functional_groups for r in subset]
    n = len(subset)

    for name, fg_names, reactor in reaction_rules:
        if len(fg_names) != n:
            continue
        for perm in permutations(range(n)):
            if all(fg_names[i] in fgs[j] for i, j in enumerate(perm)):
                for rxn in reactor(*(subset[j] for j in perm)):
                    p = rxn.products[0]
                    if x := _match(p, product):
                        return x, [f'react:{name}']
                break

    return None


__all__ = ['Reconstruct']
