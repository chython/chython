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
from itertools import combinations


class Reconstruct:
    __slots__ = ()

    def reconstruct_mapping(self) -> list[str]:
        """
        Annotate reaction by trying to reconstruct the product from reactants
        using predefined reaction templates.

        Tries in order:
        1. Standalone deprotection (reactant -> product via PG removal)
        2. Standalone protection (product -> reactant via PG removal, i.e. reverse)
        3. Single-molecule reactions (oxidize/reduce/transform)
           with optional one-pot deprotection of the product
        4. Multi-component reactions with optional deprotection

        Returns list of matched reaction/deprotection names.
        Empty list if no match found.
        If found, updates atom-to-atom mapping.
        """
        assert self.reactants, 'No reactants in reaction'
        assert len(self.products) == 1, 'Only single product reactions supported'
        self.reset_mapping()

        reactants = []
        for m in self.reactants:
            c = m.copy()
            c.clean_stereo()
            c.clean_isotopes()
            c.standardize_tautomers(prepare_molecule=False)
            reactants.append(c)

        product = self.products[0].copy()
        product.clean_stereo()
        product.clean_isotopes()
        product.standardize_tautomers(prepare_molecule=False)

        # 1. standalone deprotection: reactant -> product by removing PGs
        for r in reactants:
            if result := _deprotect(r, product):
                mapping, removed = result
                self.products[0].remap(mapping)
                return [f'deprotect:{name}' for name in removed]

        # 2. standalone protection: product -> reactant by removing PGs (reverse direction)
        for i, r in enumerate(reactants):
            if result := _deprotect(product, r):
                mapping, removed = result
                self.reactants[i].remap(mapping)
                return [f'protect:{name}' for name in removed]

        # 3. single-molecule transforms
        for i, rg in enumerate(reactants):
            # oxidations
            for name, rxn in rg.oxidize():
                if x := _match(rxn.products[0], product):
                    self.products[0].remap(x)
                    return [f'oxidize:{name}']
                if result := _deprotect(rxn.products[0], product):
                    mapping, removed = result
                    self.products[0].remap(mapping)
                    return [f'oxidize:{name}'] + [f'deprotect:{p}' for p in removed]

            # reductions
            for name, rxn in rg.reduce():
                if x := _match(rxn.products[0], product):
                    self.products[0].remap(x)
                    return [f'reduce:{name}']
                if result := _deprotect(rxn.products[0], product):
                    mapping, removed = result
                    self.products[0].remap(mapping)
                    return [f'reduce:{name}'] + [f'deprotect:{p}' for p in removed]

            # transformations
            for name, rxn in rg.transform():
                if x := _match(rxn.products[0], product):
                    self.products[0].remap(x)
                    return [f'transform:{name}']
                if result := _deprotect(rxn.products[0], product):
                    mapping, removed = result
                    self.products[0].remap(mapping)
                    return [f'transform:{name}'] + [f'deprotect:{p}' for p in removed]

        # 4. multi-component reactions
        for size in range(min(len(reactants), 4), 1, -1):
            for mols in combinations(reactants, size):
                for name, rxn in mols[0].react(*mols[1:]):
                    if len(rxn.products) != 1:
                        continue
                    if x := _match(rxn.products[0], product):
                        self.products[0].remap(x)
                        return [f'react:{name}']
                    if result := _deprotect(rxn.products[0], product):
                        mapping, removed = result
                        self.products[0].remap(mapping)
                        return [f'react:{name}'] + [f'deprotect:{p}' for p in removed]
        return []


def _match(generated, reference):
    """Match generated product to reference, with tautomer standardization fallback."""
    if x := reference.get_fast_mapping(generated):
        return x
    # tautomer fallback: standardize generated in-place (throwaway object)
    generated.standardize_tautomers(prepare_molecule=False)
    return reference.get_fast_mapping(generated)


def _deprotect(reactant, product) -> tuple[dict[int, int], list[str]] | None:
    """Try to match reactant to product by removing PGs.

    Returns (mapping, list_of_removed_pg_names) or None.
    """
    rpg = reactant.protective_groups
    ppg = product.protective_groups
    # check if all PGs in product are present in reactant
    for p, c in ppg.items():
        if p not in rpg or c > rpg[p]:
            return
    # find PGs in reactant that are not in product
    missing = {}
    for p, c in rpg.items():
        if p not in ppg or c > ppg[p]:
            missing[p] = c - ppg.get(p, 0)
    if not missing: return

    # remove PGs from reactant
    m = reactant.copy()
    for p in missing:
        m.remove_protection(p)
    m.kekule()
    m.thiele()
    if x := _match(m, product):
        return x, list(missing)


__all__ = ['Reconstruct']
