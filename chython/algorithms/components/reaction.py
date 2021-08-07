# -*- coding: utf-8 -*-
#
#  Copyright 2020, 2021 Ramil Nugmanov <nougmanoff@protonmail.com>
#  Copyright 2020 Ravil Mukhametgaleev <sonic-mc@mail.ru>
#  Copyright 2021 Timur Gimadiev <timur.gimadiev@gmail.com>
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
from CachedMethods import cached_property
from collections import ChainMap
from itertools import chain, product
from typing import Tuple, Iterator, Union, TYPE_CHECKING
from ...containers import molecule  # cyclic imports resolve
from ...exceptions import MappingError


if TYPE_CHECKING:
    from chython import ReactionContainer


class ReactionComponents:
    __slots__ = ()

    @cached_property
    def centers_list(self: 'ReactionContainer') -> Tuple[Tuple[int, ...], ...]:
        """
        Reaction centers with leaving and coming groups.
        """
        if not self:
            return ()  # no rc

        cgr = ~self
        bonds = cgr._bonds
        all_groups = {x for x in self.reactants for x in x} ^ {x for x in self.products for x in x}
        all_groups = cgr._connected_components({n: bonds[n].keys() & all_groups for n in all_groups})
        centers_list = list(cgr.centers_list)

        for x in all_groups:
            intersection = []
            for i, y in enumerate(centers_list):
                if not x.isdisjoint(y):
                    intersection.append(i)
            if intersection:
                for i in reversed(intersection):
                    x.update(centers_list.pop(i))
                centers_list.append(x)
        return tuple(tuple(x) for x in centers_list)

    @cached_property
    def extended_centers_list(self: Union['ReactionContainer', 'ReactionComponents']) -> Tuple[Tuple[int, ...], ...]:
        """
        Additionally to `centers_list` include:
        * First environment of dynamic atoms.
        * Whole formed cycles. For condensed cycles smallest is taken.
        * Whole aromatic cycle with at least one dynamic atom.
        * Whole small (3, 4) cycle with at least one dynamic atom.
        * Double or triple bonds connected to previous atoms.

        Note for multiple RCs intersection possible. Use `enumerate_centers` to prevent unobvious RCs.
        """
        cgr = ~self
        bonds = cgr._bonds
        center_atoms = set(cgr.center_atoms)

        formed_rings = {}
        small_aromatic_rings = set()
        for r in cgr.sssr:
            if len(r) < 5 and not center_atoms.isdisjoint(r):
                small_aromatic_rings.add(frozenset(r))

            n, m = r[0], r[-1]
            if bonds[n][m].order is None:
                nm = frozenset((n, m))
                if nm not in formed_rings or len(formed_rings[nm]) > len(r):
                    formed_rings[nm] = r
            for n, m in zip(r, r[1:]):
                if bonds[n][m].order is None:
                    nm = frozenset((n, m))
                    if nm not in formed_rings or len(formed_rings[nm]) > len(r):
                        formed_rings[nm] = r

        for r in cgr.aromatic_rings:
            if not center_atoms.isdisjoint(r):
                small_aromatic_rings.add(frozenset(r))

        out = []
        for rc in self.centers_list:
            c = center_atoms.intersection(rc)
            fe = {m for n in c for m in bonds[n]}  # add first environment
            # add double or triple bond to first env
            fe |= {m for n in fe for m, b in bonds[n].items() if m not in fe and b.order in (2, 3)}
            fe.update(rc)  # add leaving|coming groups and alone center atoms

            for rk, r in formed_rings.items():  # add formed rings to RC
                if c.issuperset(rk):
                    fe.update(r)
            for r in small_aromatic_rings:  # add small or aromatic rings with dyn atoms
                if not c.isdisjoint(r):
                    fe.update(r)
            out.append(tuple(fe))
        return tuple(out)

    def enumerate_centers(self: 'ReactionContainer') -> Iterator['ReactionContainer']:
        """
        Get all possible single stage reactions from multistage.
        Note multicomponent molecules (salts etc) can be treated incorrectly.
        """
        if len(self.centers_list) > 1:
            centers_list = self.centers_list

            charges = ChainMap(*(x._charges for x in self.reactants))
            radicals = ChainMap(*(x._radicals for x in self.reactants))
            bonds = ChainMap(*(x._bonds for x in self.reactants))
            atoms = ChainMap(*(x._atoms for x in self.reactants))
            p_charges = ChainMap(*(x._charges for x in self.products))
            p_radicals = ChainMap(*(x._radicals for x in self.products))
            p_bonds = ChainMap(*(x._bonds for x in self.products))
            p_atoms = ChainMap(*(x._atoms for x in self.products))

            centers = {x for x in centers_list for x in x}
            common = {x for x in chain(self.reactants, self.products) for x in x if x not in centers}
            reactants = {x for x in self.reactants for x in x}
            products = {x for x in self.products for x in x}

            common_molecule = molecule.MoleculeContainer()
            for n in common:
                common_molecule.add_atom(atoms[n].copy(), n, charge=charges[n], is_radical=radicals[n])
            seen = set()
            for n in common:
                seen.add(n)
                for m, b in bonds[n].items():
                    if m not in seen and m in common:
                        common_molecule.add_bond(n, m, b.copy())

            products_bonds = {}
            reactants_bonds = {}
            common_bonds = []
            seen = set()
            p_seen = set()
            for c in centers_list:
                not_rc = centers.difference(c)
                reactants_bonds[c] = (c_bonds, c_atoms) = [], reactants.intersection(c)
                for n in c_atoms:
                    seen.add(n)
                    for m, b in bonds[n].items():
                        if m not in seen and m in reactants:
                            if m in not_rc:
                                common_bonds.append((n, m, b))
                            else:
                                c_bonds.append((n, m, b))
                products_bonds[c] = (c_bonds, c_atoms) = [], products.intersection(c)
                for n in c_atoms:
                    p_seen.add(n)
                    for m, b in p_bonds[n].items():
                        if m not in p_seen and m in products and m not in not_rc:
                            c_bonds.append((n, m, b))

            for rc in range(len(centers_list)):
                not_rc = centers_list[:rc] + centers_list[rc + 1:]
                rc = centers_list[rc]
                for combo in list(product((False, True), repeat=len(not_rc))):
                    r = common_molecule.copy()
                    p = common_molecule.copy()

                    for n in reactants_bonds[rc][1]:
                        r.add_atom(atoms[n].copy(), n, charge=charges[n], is_radical=radicals[n])
                    for n in products_bonds[rc][1]:
                        p.add_atom(p_atoms[n].copy(), n, charge=p_charges[n], is_radical=p_radicals[n])

                    for is_p, center in zip(combo, not_rc):
                        if is_p:
                            c_bonds, c_atoms = products_bonds[center]
                            for n in c_atoms:
                                r.add_atom(p_atoms[n].copy(), n, charge=p_charges[n], is_radical=p_radicals[n])
                                p.add_atom(p_atoms[n].copy(), n, charge=p_charges[n], is_radical=p_radicals[n])
                        else:
                            c_bonds, c_atoms = reactants_bonds[center]
                            for n in c_atoms:
                                r.add_atom(atoms[n].copy(), n, charge=charges[n], is_radical=radicals[n])
                                p.add_atom(atoms[n].copy(), n, charge=charges[n], is_radical=radicals[n])
                        for n, m, b in c_bonds:
                            r.add_bond(n, m, b.copy())
                            p.add_bond(n, m, b.copy())

                    for n, m, b in products_bonds[rc][0]:
                        p.add_bond(n, m, b.copy())
                    for n, m, b in reactants_bonds[rc][0]:
                        r.add_bond(n, m, b.copy())
                    for n, m, b in common_bonds:
                        r.add_bond(n, m, b.copy())
                        p.add_bond(n, m, b.copy())
                    yield self.__class__(r.split(), p.split(), [x.copy() for x in self.reagents])
        else:
            cp = self.copy()
            cp.meta.clear()
            yield cp

    def remove_reagents(self: 'ReactionContainer', *, keep_reagents: bool = False) -> bool:
        """
        Preprocess reaction according to mapping, using the following idea: molecules(each separated graph) will be
        placed to reagents if it is not changed in the reaction (no bonds, charges reorders)

        Return True if any reagent found.
        """
        cgr = ~self
        if cgr.center_atoms:
            active = set(cgr.center_atoms)
            reactants = []
            products = []
            reagents = set(self.reagents)
            for i in self.reactants:
                if not active.isdisjoint(i):
                    reactants.append(i)
                else:
                    reagents.add(i)
            for i in self.products:
                if not active.isdisjoint(i):
                    products.append(i)
                else:
                    reagents.add(i)
            if keep_reagents:
                tmp = []
                for m in self.reagents:
                    if m in reagents:
                        tmp.append(m)
                        reagents.discard(m)
                tmp.extend(reagents)
                reagents = tuple(tmp)
            else:
                reagents = ()

            if len(reactants) != len(self.reactants) or len(products) != len(self.products) or \
                    len(reagents) != len(self.reagents):
                self._ReactionContainer__reactants = tuple(reactants)
                self._ReactionContainer__products = tuple(products)
                self._ReactionContainer__reagents = reagents
                self.flush_cache()
                self.fix_positions()
                return True
            return False
        raise MappingError("Reaction center is absent according to mapping")

    def contract_ions(self: 'ReactionContainer') -> bool:
        """
        Contract ions into salts (Molecules with disconnected components).
        Note: works only for unambiguous cases. e.g. equal anions/cations and different or equal cations/anions.

        Return True if any ions contracted.
        """
        neutral, cations, anions, total = _sift_ions(self.reagents)
        salts = _contract_ions(anions, cations, total)
        if salts:
            neutral.extend(salts)
            self._ReactionContainer__reagents = tuple(neutral)
            changed = True
        else:
            changed = False

        neutral, cations, anions, total = _sift_ions(self.reactants)
        salts = _contract_ions(anions, cations, total)
        if salts:
            anions_order = {frozenset(m): n for n, m in enumerate(anions)}
            cations_order = {frozenset(m): n for n, m in enumerate(cations)}
            neutral.extend(salts)
            self._ReactionContainer__reactants = tuple(neutral)
            changed = True
        else:
            anions_order = cations_order = {}

        neutral, cations, anions, total = _sift_ions(self.products)
        if cations and anions:
            anions.sort(key=lambda x: anions_order.get(frozenset(x), -1))
            cations.sort(key=lambda x: cations_order.get(frozenset(x), -1))
        salts = _contract_ions(anions, cations, total)
        if salts:
            neutral.extend(salts)
            self._ReactionContainer__products = tuple(neutral)
            changed = True

        if changed:
            self.flush_cache()
            self.fix_positions()
            return True
        return False


def _sift_ions(mols):
    anions = []
    cations = []
    neutral = []
    total = 0
    for m in mols:
        c = int(m)
        total += c
        if c > 0:
            cations.append(m)
        elif c < 0:
            anions.append(m)
        else:
            neutral.append(m)
    return neutral, cations, anions, total


def _contract_ions(anions, cations, total):
    salts = []
    if not anions or not cations:  # nothing to contract
        return
    # check ambiguous cases
    if total > 0:
        if len(cations) > 1:  # deficit of anions
            return
    elif total < 0:
        if len(anions) > 1:  # deficit of cations
            return
    elif len(set(anions)) > 1 and len(set(cations)) > 1:  # different anions and cations
        return

    anions = anions.copy()
    cations = cations.copy()
    while anions:
        ct = cations.pop()
        an = anions.pop()
        shift_x = ct._fix_plane_mean(0) + 1
        shift_x = an._fix_plane_mean(shift_x)
        salt = ct | an
        while True:
            c = int(salt)
            if c > 0:
                an = anions.pop()
                shift_x = an._fix_plane_mean(shift_x) + 1
                salt = salt | an
            elif c < 0:
                ct = cations.pop()
                shift_x = ct._fix_plane_mean(shift_x) + 1
                salt = salt | ct
            else:
                break
        salts.append(salt)
    return salts


__all__ = ['ReactionComponents']
