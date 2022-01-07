# -*- coding: utf-8 -*-
#
#  Copyright 2022 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from itertools import count
from typing import List, Tuple, Type, TYPE_CHECKING, Union
from ._mapping import rules


if TYPE_CHECKING:
    from chython import ReactionContainer


class Mapping:
    __slots__ = ()

    def reset_mapping(self) -> bool:
        """
        Do atom-to-atom mapping. Return True if mapping changed.
        """

    def fix_mapping(self: Union['ReactionContainer', 'Mapping'], *, logging: bool = False) -> \
            Union[bool, List[Union[Tuple[str, List[int]], Tuple[int, str, str, Tuple[int, ...]]]]]:
        """
        Fix atom-to-atom mapping of some functional groups. Return True if found AAM errors.
        """
        log = []
        seen = set()
        if not self:
            return False

        for r_pattern, p_pattern, fix in rules:
            found = []
            for m in self.reactants:
                for mapping in r_pattern.get_mapping(m, automorphism_filter=False):
                    if mapping[1] not in seen:
                        found.append(({fix.get(k, k): v for k, v in mapping.items()},
                                      {mapping[k]: mapping[v] for k, v in fix.items()}))

            if not found:
                continue
            for m in self.products:
                for mapping in p_pattern.get_mapping(m, automorphism_filter=False):
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
            for mapping in bad_query.get_mapping(cgr, automorphism_filter=False):
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
                if flag_m:
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
        fix_rules = []
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
            fix_rules.append((bad_query, good_query, fix, valid))

        cls.__class_cache__[cls] = {'_Mapping__remapping_compiled_rules': tuple(fix_rules)}

    @class_cached_property
    def __remapping_compiled_rules(self):
        return ()


__all__ = ['Mapping']
