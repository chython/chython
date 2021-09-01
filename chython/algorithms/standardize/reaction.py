# -*- coding: utf-8 -*-
#
#  Copyright 2018-2021 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from CachedMethods import class_cached_property
from collections import defaultdict
from itertools import count
from typing import List, Tuple, Type, TYPE_CHECKING, Union
from ._mapping import rules as mapping_rules
from ...exceptions import MappingError


if TYPE_CHECKING:
    from chython import ReactionContainer


class StandardizeReaction:
    __slots__ = ()

    def canonicalize(self: 'ReactionContainer', fix_mapping: bool = True, *, logging=False) -> \
            Union[bool, List[Tuple[int, Tuple[int, ...], int, str]]]:
        """
        Convert molecules to canonical forms of functional groups and aromatic rings without explicit hydrogens.
        Return True if in any molecule found not canonical group.

        :param fix_mapping: Search AAM errors of functional groups.
        :param logging: return log from molecules with index of molecule.
            Otherwise return True if these groups found in any molecule.
        """
        if logging:
            total = []
        else:
            total = False
        for n, m in enumerate(self.molecules()):
            out = m.canonicalize(logging=logging)
            if out:
                if logging:
                    total.extend((n, *x) for x in out)
                else:
                    total = True

        if fix_mapping and self.fix_mapping():
            if logging:
                total.append((-1, (), -1, 'mapping fixed'))
                return total
            return True

        if total:
            self.flush_cache()
        return total

    def standardize(self: 'ReactionContainer', fix_mapping: bool = True, *, logging=False) -> \
            Union[bool, Tuple[int, Tuple[int, ...], int, str]]:
        """
        Standardize functional groups.

        :param fix_mapping: Search AAM errors of functional groups.
        :param logging: return log from molecules with index of molecule at first position.
            Otherwise return True if these groups found in any molecule.
        """
        if logging:
            total = []
        else:
            total = False
        for n, m in enumerate(self.molecules()):
            out = m.standardize(logging=logging)
            if out:
                if logging:
                    total.extend((n, *x) for x in out)
                else:
                    total = True

        if fix_mapping and self.fix_mapping():
            if logging:
                total.append((-1, (), -1, 'mapping fixed'))
                return total
            return True

        if total:
            self.flush_cache()
        return total

    def standardize_charges(self: 'ReactionContainer', *, logging=False) -> \
            Union[bool, List[Tuple[int, Tuple[int, ...]]]]:
        """
        Set canonical positions of charges in heterocycles and some nitrogen compounds.

        :param logging: return list of changed molecules with tuple of changed atoms.
        """
        if logging:
            total = []
        else:
            total = False
        for n, m in enumerate(self.molecules()):
            out = m.standardize_charges(logging=logging)
            if out:
                if logging:
                    total.append((n, tuple(out)))
                else:
                    total = True
        if total:
            self.flush_cache()
        return total

    def thiele(self: 'ReactionContainer') -> bool:
        """
        Convert structures to aromatic form.
        Return True if in any molecule found kekule ring
        """
        total = False
        for m in self.molecules():
            if m.thiele() and not total:
                total = True
        if total:
            self.flush_cache()
        return total

    def kekule(self: 'ReactionContainer') -> bool:
        """
        Convert structures to kekule form. Works only for Molecules.
        Return True if in any molecule found aromatic ring
        """
        total = False
        for m in self.molecules():
            if m.kekule() and not total:
                total = True
        if total:
            self.flush_cache()
        return total

    def clean_isotopes(self: 'ReactionContainer') -> bool:
        """
        Clean isotope marks for all molecules in reaction.
        Returns True if in any molecule found isotope.
        """
        flag = False
        for m in self.molecules():
            if m.clean_isotopes() and not flag:
                flag = True
        if flag:
            self.flush_cache()
        return flag

    def clean_stereo(self: 'ReactionContainer'):
        """
        Remove stereo data
        """
        for m in self.molecules():
            m.clean_stereo()
        self.flush_cache()

    def check_valence(self: 'ReactionContainer') -> List[Tuple[int, Tuple[int, ...]]]:
        """
        Check valences of all atoms of all molecules.

        Works only on molecules with aromatic rings in Kekule form.
        :return: list of invalid molecules with invalid atoms lists
        """
        out = []
        for n, m in enumerate(self.molecules()):
            c = m.check_valence()
            if c:
                out.append((n, tuple(c)))
        return out

    def fix_mapping(self: Union['ReactionContainer', 'StandardizeReaction'], *, logging: bool = False) -> \
            Union[bool, List[Union[Tuple[str, List[int]], Tuple[int, str, str, Tuple[int, ...]]]]]:
        """
        Fix atom-to-atom mapping of some functional groups. Return True if found AAM errors.
        """
        if logging:
            log = []
        seen = set()
        if not self:
            return False

        for r_pattern, p_pattern, fix in mapping_rules:
            found = []
            for m in self.reactants:
                for mapping in r_pattern.get_mapping(m, optimize=False, automorphism_filter=False):
                    if mapping[1] not in seen:
                        found.append(({fix.get(k, k): v for k, v in mapping.items()},
                                      {mapping[k]: mapping[v] for k, v in fix.items()}))

            if not found:
                continue
            for m in self.products:
                for mapping in p_pattern.get_mapping(m, optimize=False, automorphism_filter=False):
                    atom = mapping[1]
                    if atom in seen:
                        continue
                    for n, (k, v) in enumerate(found):
                        if k == mapping:
                            break
                    else:
                        continue

                    del found[n]
                    m.remap(v)
                    seen.add(atom)
                    if logging:
                        log.append(('group remap', list(v)))
        if seen:
            self.flush_cache()
            flag = True
            seen = set()
        else:
            flag = False

        for rule_num, (bad_query, good_query, fix, valid) in enumerate(self.__remapping_compiled_rules):
            cgr = ~self
            first_free = max(cgr) + 1
            free_number = count(first_free)
            cgr_c = set(cgr.center_atoms)
            del self.__dict__['__cached_method_compose']

            flag_m = False
            for mapping in bad_query.get_mapping(cgr, optimize=False, automorphism_filter=False):
                if not seen.isdisjoint(mapping.values()):  # prevent matching same RC
                    continue
                mapping = {mapping[n]: next(free_number) if m is None else mapping[m] for n, m in fix.items()}
                flag_m = True
                reverse = {m: n for n, m in mapping.items()}
                for m in self.products:
                    m.remap(mapping)

                check = ~self
                check_c = set(check.center_atoms)
                delta = check_c - cgr_c

                for m in good_query.get_mapping(check, automorphism_filter=False):
                    if valid.issubset(m) and delta.issubset(m.values()):
                        seen.update(mapping)
                        if logging:
                            log.append((rule_num, str(bad_query), str(good_query), tuple(mapping.values())))
                        flag = True
                        break
                else:
                    # restore old mapping
                    for m in self.products:
                        m.remap(reverse)
                    del self.__dict__['__cached_method_compose']
                    free_number = count(first_free)
                    continue
                break
            else:
                if logging and flag_m:
                    log.append((rule_num, str(bad_query), str(good_query), ()))
        if seen:
            self.flush_cache()
        if logging:
            return log
        return flag

    @classmethod
    def load_remapping_rules(cls: Type['ReactionContainer'], reactions):
        """
        Load AAM fixing rules. Required pairs of bad mapped and good mapped reactions.
        Reactants in pairs should be fully equal (equal molecules and equal atom orders).
        Products should be equal but with different atom numbers.
        """
        rules = []
        for bad, good in reactions:
            if str(bad) != str(good):
                raise ValueError('bad and good reaction should be equal')

            cgr_good, cgr_bad = ~good, ~bad
            gc = cgr_good.substructure({y for x in cgr_good.center_atoms for y in cgr_good._bonds[x]} |
                                       set(cgr_good.center_atoms))
            bc = cgr_bad.substructure({y for x in cgr_bad.center_atoms for y in cgr_bad._bonds[x]} |
                                      set(cgr_bad.center_atoms))
            atoms = set(bc) | set(gc)

            pr_g, pr_b, re_g, re_b = set(), set(), set(), set()
            for pr in good.products:
                pr_g.update(pr)
            for pr in bad.products:
                pr_b.update(pr)
            for pr in good.reactants:
                re_g.update(pr)
            for pr in bad.reactants:
                re_b.update(pr)
            atoms.update((re_b.difference(pr_b)).intersection(pr_g))

            strange_atoms = pr_b.difference(pr_g)
            atoms.update(strange_atoms)

            bad_query = cgr_bad.substructure(atoms.intersection(cgr_bad))
            good_query = cgr_good.substructure(atoms.intersection(cgr_good))

            fix = {}
            for mb, mg in zip(bad.products, good.products):
                fix.update({k: v for k, v in zip(mb, mg) if k != v and k in atoms})

            valid = set(fix).difference(strange_atoms)
            rules.append((bad_query, good_query, fix, valid))

        cls.__class_cache__[cls] = {'_StandardizeReaction__remapping_compiled_rules': tuple(rules)}

    @class_cached_property
    def __remapping_compiled_rules(self):
        return ()

    def implicify_hydrogens(self: 'ReactionContainer') -> int:
        """
        Remove explicit hydrogens if possible

        :return: number of removed hydrogens
        """
        total = 0
        for m in self.molecules():
            total += m.implicify_hydrogens()
        if total:
            self.flush_cache()
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
            maps = m.explicify_hydrogens(return_maps=True, start_map=start_map + 1)
            if maps:
                for n, h in maps:
                    mapping[n].append(h)
                start_map = maps[-1][1]
                total += len(maps)

        for m in self.reagents:
            maps = m.explicify_hydrogens(return_maps=True, start_map=start_map + 1)
            if maps:
                start_map = maps[-1][1]
                total += len(maps)

        for m in self.products:
            maps = m.explicify_hydrogens(return_maps=True, start_map=start_map + 1)
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
            self.flush_cache()
        return total

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


__all__ = ['StandardizeReaction']
