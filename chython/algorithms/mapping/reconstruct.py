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

    def reconstruct_mapping(self) -> bool:
        """
        Annotate reaction by trying to reconstruct the product from reactants
        using predefined reaction templates.

        Tries in order:
        1. Standalone deprotection (reactant -> product via PG removal)
        2. Standalone protection (product -> reactant via PG removal, i.e. reverse)
        3. Multi/single-component reactions (react/oxidize/reduce/transform)
           with optional one-pot deprotection of the product

        Returns True if a reaction was found, False otherwise.
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
            reactants.append(c)

        product = self.products[0].copy()
        product.clean_stereo()
        product.clean_isotopes()

        # 1. standalone deprotection: reactant -> product by removing PGs
        for r in reactants:
            if x := _deprotect(r, product):
                # remap product atom to reactant atom
                self.products[0].remap(x)
                return True

        # 2. standalone protection: product -> reactant by removing PGs (reverse direction)
        for i, r in enumerate(reactants):
            if x := _deprotect(product, r):
                self.reactants[i].remap(x)
                return True

        # 3. single-molecule transforms
        for i, rg in enumerate(reactants):
            # oxidations
            for name, rxn in rg.oxidize():
                if x := product.get_fast_mapping(rxn.products[0]):
                    self.products[0].remap(x)
                    return True
                if x := _deprotect(rxn.products[0], product):
                    self.products[0].remap(x)
                    return True

            # reductions
            for name, rxn in rg.reduce():
                if x := product.get_fast_mapping(rxn.products[0]):
                    self.products[0].remap(x)
                    return True
                if x := _deprotect(rxn.products[0], product):
                    self.products[0].remap(x)
                    return True

            # transformations
            for name, rxn in rg.transform():
                if x := product.get_fast_mapping(rxn.products[0]):
                    self.products[0].remap(x)
                    return True
                if x := _deprotect(rxn.products[0], product):
                    self.products[0].remap(x)
                    return True

        # 4. multi-component reactions
        for size in range(min(len(reactants), 4), 1, -1):
            for mols in combinations(reactants, size):
                for name, rxn in mols[0].react(*mols[1:]):
                    if len(rxn.products) != 1:
                        continue
                    if x := product.get_fast_mapping(rxn.products[0]):
                        self.products[0].remap(x)
                        return True
                    if x := _deprotect(rxn.products[0], product):
                        self.products[0].remap(x)
                        return True
        return False


def _deprotect(reactant, product) -> dict[int, int] | None:
    rpg = reactant.protective_groups
    ppg = product.protective_groups
    # check if all PGs in product are present in reactant
    for p, c in ppg.items():
        if p not in rpg or c >= rpg[p]:
            return
    # find PGs in reactant that are not in product
    missing = {}
    for p, c in rpg.items():
        if p not in ppg or c > ppg[p]:
            missing[p] = c - ppg.get(p, 0)
    if not missing: return

    # remove PGs from reactant
    m = reactant.copy()
    for p, c in missing.items():
        m.remove_protection(p)
    if m == product: return product.get_fast_mapping(m)


__all__ = ['Reconstruct']
