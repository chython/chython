# -*- coding: utf-8 -*-
#
#  Copyright 2018-2025 Ramil Nugmanov <nougmanoff@protonmail.com>
#  Copyright 2021 Timur Gimadiev <timur.gimadiev@gmail.com>
#  Copyright 2024 Philippe Gantzer <p.gantzer@icredd.hokudai.ac.jp>
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
from collections import defaultdict
from typing import List, Tuple, TYPE_CHECKING, Union
from ._reagents import *
from ...exceptions import MappingError


if TYPE_CHECKING:
    from chython import ReactionContainer


class StandardizeReaction:
    __slots__ = ()

    def canonicalize(self: 'ReactionContainer', *, fix_mapping: bool = True, logging=False, fix_tautomers=True) -> \
            Union[bool, List[Tuple[int, Tuple[int, ...], int, str]]]:
        """
        Convert molecules to canonical forms of functional groups and aromatic rings without explicit hydrogens.
        Return True if in any molecule found not canonical group.

        :param fix_mapping: Search AAM errors of functional groups.
        :param logging: return log from molecules with index of molecule.
            Otherwise, return True if these groups found in any molecule.
        :param fix_tautomers: convert tautomers to canonical forms.
        """
        total = []
        for n, m in enumerate(self.molecules()):
            total.extend((n, *x) for x in m.canonicalize(logging=True, fix_tautomers=fix_tautomers))

        if fix_mapping:
            total.extend((-1, x, -1, m) for m, x in self.fix_groups_mapping(logging=True))

        if total:
            self.flush_cache(keep_molecule_cache=True)
        if logging:
            return total
        return bool(total)

    def standardize(self: 'ReactionContainer', *, fix_mapping: bool = True, logging=False, fix_tautomers=True) -> \
            Union[bool, List[Tuple[int, Tuple[int, ...], int, str]]]:
        """
        Fix functional groups representation.
        Return True if in any molecule fixed group.

        Deprecated method. Use `canonicalize` directly.

        :param fix_mapping: Search AAM errors of functional groups.
        :param logging: return log from molecules with index of molecule.
            Otherwise, return True if these groups found in any molecule.
        :param fix_tautomers: convert tautomers to canonical forms.
        """
        total = []
        for n, m in enumerate(self.molecules()):
            total.extend((n, *x) for x in m.standardize(logging=True, fix_tautomers=fix_tautomers))

        if fix_mapping:
            total.extend((-1, x, -1, m) for m, x in self.fix_groups_mapping(logging=True))

        if total:
            self.flush_cache(keep_molecule_cache=True)
        if logging:
            return total
        return bool(total)

    def thiele(self: 'ReactionContainer', *, fix_tautomers=True) -> bool:
        """
        Convert structures to aromatic form.
        Return True if in any molecule found kekule ring

        :param fix_tautomers: convert tautomers to canonical forms.
        """
        total = False
        for m in self.molecules():
            if m.thiele(fix_tautomers=fix_tautomers):
                total = True
        if total:
            self.flush_cache(keep_molecule_cache=True)
        return total

    def kekule(self: 'ReactionContainer', *, buffer_size=7, ignore_pyrrole_hydrogen=False) -> bool:
        """
        Convert structures to a kekule form.
        Return True if in any molecule found aromatic ring

        :param buffer_size: number of attempts of pyridine form searching.
        :param ignore_pyrrole_hydrogen: ignore hydrogen on pyrrole to fix invalid rings like Cn1cc[nH]c1.
        """
        total = False
        for m in self.molecules():
            if m.kekule(buffer_size=buffer_size, ignore_pyrrole_hydrogen=ignore_pyrrole_hydrogen):
                total = True
        if total:
            self.flush_cache(keep_molecule_cache=True)
        return total

    def clean_isotopes(self: 'ReactionContainer') -> bool:
        """
        Clean isotope marks for all molecules in reaction.
        Returns True if in any molecule found isotope.
        """
        flag = False
        for m in self.molecules():
            if m.clean_isotopes():
                flag = True
        if flag:
            self.flush_cache(keep_molecule_cache=True)
        return flag

    def clean_stereo(self: 'ReactionContainer'):
        """
        Remove stereo data
        """
        for m in self.molecules():
            m.clean_stereo()
        self.flush_cache(keep_molecule_cache=True)

    def check_valence(self: 'ReactionContainer') -> List[Tuple[int, Tuple[int, ...]]]:
        """
        Check valences of all atoms of all molecules.

        Works only on molecules with aromatic rings in Kekule form.
        :return: list of invalid molecules with invalid atoms lists
        """
        out = []
        for n, m in enumerate(self.molecules()):
            if c := m.check_valence():
                out.append((n, tuple(c)))
        return out

    def implicify_hydrogens(self: 'ReactionContainer') -> int:
        """
        Remove explicit hydrogens if possible.

        :return: number of removed hydrogens.
        """
        total = 0
        for m in self.molecules():
            total += m.implicify_hydrogens()
        if total:
            self.flush_cache(keep_molecule_cache=True)
        return total

    def explicify_hydrogens(self: 'ReactionContainer') -> int:
        """
        Add explicit hydrogens to atoms

        :return: number of added atoms
        """
        total = 0
        start_map = 0
        for m in self.molecules():
            map_ = max(m, default=0)
            if map_ > start_map:
                start_map = map_

        mapping = defaultdict(list)
        for m in self.reactants:
            maps = m.explicify_hydrogens(_return_map=True, start_map=start_map + 1)
            if maps:
                for n, h in maps:
                    mapping[n].append(h)
                start_map = maps[-1][1]
                total += len(maps)

        for m in self.reagents:
            maps = m.explicify_hydrogens(_return_map=True, start_map=start_map + 1)
            if maps:
                start_map = maps[-1][1]
                total += len(maps)

        for m in self.products:
            maps = m.explicify_hydrogens(_return_map=True, start_map=start_map + 1)
            if maps:
                total += len(maps)
                remap = {}
                free = []
                for n, h in maps:
                    if n in mapping and mapping[n]:
                        remap[h] = mapping[n].pop()
                        free.append(h)
                    elif free:
                        remap[h] = start_map = free.pop(0)
                    else:
                        start_map = h
                m.remap(remap)

        if total:
            self.flush_cache(keep_molecule_cache=True)
        return total

    def remove_reagents(self, *, keep_reagents: bool = False, mapping: bool = True) -> bool:
        """
        Place molecules, except reactants, to reagents list. Reagents - molecules which atoms not presented in products.
        Mapping based approach remove molecules without reaction center.
        Rule based approach remove equal molecules in reactants and products, and predefined reactants.

        :param mapping: use atom-to-atom mapping to detect reagents, otherwise use predefined list of common reagents.
        :param keep_reagents: delete reagents if False

        Return True if any reagent found.
        """
        if mapping:
            return self.__remove_reagents_mapping(keep_reagents)
        return self.__remove_reagents_rules(keep_reagents)

    def __remove_reagents_rules(self: 'ReactionContainer', keep_reagents):
        if not self.reactants or not self.products:  # there is no reaction
            return False

        reactants_st1 = []
        products_st1 = []
        reagents_st1 = set(self.reagents)

        # find the same molecules in reactants and products
        for m in self.reactants:
            if m in self.products:
                reagents_st1.add(m)
            else:
                reactants_st1.append(m)
        for m in self.products:
            if m in self.reactants:
                reagents_st1.add(m)
            else:
                products_st1.append(m)
        if not reactants_st1 or not products_st1:
            return False  # keep bad reaction as is

        reactants_st2 = []
        products_st2 = []
        reagents_st2 = reagents_st1.copy()

        # filter out predefined reagents
        for m in reactants_st1:
            if m in reagents_set:
                reagents_st2.add(m)
            else:
                reactants_st2.append(m)
        for m in products_st1:
            if m in reagents_set:
                reagents_st2.add(m)
            else:
                products_st2.append(m)
        if not reactants_st2 or not products_st2:  # reaction contains only simple molecules. roll-back to step 1
            reactants_st2 = reactants_st1
            products_st2 = products_st1
            reagents_st2 = reagents_st1

        # remove reagents from reactants
        tmp = []
        for m in self.reagents:
            tmp.append(m)
            if m in reagents_st2:
                reagents_st2.discard(m)
        tmp.extend(reagents_st2)
        reagents = tuple(tmp) if keep_reagents else ()

        self._reactants = tuple(reactants_st2)
        self._products = tuple(products_st2)
        self._reagents = reagents
        self.fix_positions()
        return True

    def __remove_reagents_mapping(self: 'ReactionContainer', keep_reagents):
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

            # remove reagents from reactants
            tmp = []
            for m in self.reagents:
                tmp.append(m)
                if m in reagents:
                    reagents.discard(m)
            tmp.extend(reagents)
            reagents = tuple(tmp) if keep_reagents else ()

            if len(reactants) != len(self.reactants) or len(products) != len(self.products) or len(reagents) != len(self.reagents):
                self._reactants = tuple(reactants)
                self._products = tuple(products)
                self._reagents = reagents
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
            self._reagents = tuple(neutral)
            changed = True
        else:
            changed = False

        neutral, cations, anions, total = _sift_ions(self.reactants)
        salts = _contract_ions(anions, cations, total)
        if salts:
            anions_order = {frozenset(m): n for n, m in enumerate(anions)}
            cations_order = {frozenset(m): n for n, m in enumerate(cations)}
            neutral.extend(salts)
            self._reactants = tuple(neutral)
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
            self._products = tuple(neutral)
            changed = True

        if changed:
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
    if not anions or not cations:  # nothing to contract
        return
    # check ambiguous cases
    if total > 0:
        if len(cations) > 1:  # deficit of anions
            # we have an excess of cations. we can't assign anions univocally
            return  # unite is ambiguous
        salt = cations[0]
        shift_x = salt._fix_plane_mean(0) + 1
        for x in anions:
            shift_x = x._fix_plane_mean(shift_x) + 1
            salt = salt | x
        return [salt]
    elif total < 0:
        if len(anions) > 1:  # deficit of cations
            # we have an excess of anions. we can't assign cations univocally
            return  # unite is ambiguous
        salt = anions[0]
        shift_x = salt._fix_plane_mean(0) + 1
        for x in cations:
            shift_x = x._fix_plane_mean(shift_x) + 1
            salt = salt | x
        return [salt]
    elif len(set(anions)) > 1 and len(set(cations)) > 1:  # different anions and cations
        return

    salts = []
    anions = anions.copy()
    cations = cations.copy()
    while anions:
        ct = cations.pop()
        an = anions.pop()
        shift_x = ct._fix_plane_mean(0) + 1
        shift_x = an._fix_plane_mean(shift_x) + 1
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


__all__ = ['StandardizeReaction']
