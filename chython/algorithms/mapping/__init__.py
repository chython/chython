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
from itertools import count, repeat
from numpy import ix_, ones_like, errstate, unravel_index, nanargmax, zeros, array
from pkg_resources import resource_stream
from typing import List, Tuple, Type, TYPE_CHECKING, Union
from ._mapping import rules


if TYPE_CHECKING:
    from chython import ReactionContainer


class Mapping:
    __slots__ = ()

    def reset_mapping(self: Union['ReactionContainer', 'Mapping'], *, multiplier=90.,
                      keep_reactants_numbering=False) -> bool:
        """
        Do atom-to-atom mapping. Return True if mapping changed.
        """
        r_map, p_map, r_atoms, p_atoms, r_adj, p_adj, ram, pam, ra, pa = self.__prepare_remapping()
        am = self.__get_attention()

        am = am[ix_(pam, ram)] + am[ix_(ram, pam)].T  # sum of reactants to products attention and vice-versa
        am *= p_atoms[:, None] == r_atoms  # equal atom type attention

        mm = ones_like(am)
        mapping = {}
        for _ in range(pa):  # iteratively map each product atom to reactant
            # todo: optimize
            nam = am * mm
            with errstate(invalid='ignore'):
                nam /= nam[:, None, :].sum(2)  # normalized attention

            # select highest attention
            try:
                i, j = unravel_index(nanargmax(nam), nam.shape)
            except ValueError:  # no more products atoms in reactants
                # mark as unmapped
                for n in set(p_map).difference(mapping):
                    mapping[n] = 0
                break

            mapping[p_map[i]] = r_map[j]
            # update mm
            mm[ix_(p_adj[i], r_adj[j])] *= multiplier
            mm[i] = 0  # mask already mapped product atom
            mm[:, j] = 0  # mask already mapped reactant. todo: work with side products

        if any(n != m for n, m in mapping.items()):  # old mapping changed
            if keep_reactants_numbering:
                r_mapping = {n: n for n in r_map}
            else:
                r_mapping = {m: n for n, m in enumerate(r_map, 1)}  # remap reactants to contiguous range
                for m in self.reactants:
                    m.remap(r_mapping)

            p_mapping = {}
            nm = ra
            for n, m in mapping.items():
                if m := r_mapping.get(m):
                    p_mapping[n] = m
                else:  # not found in reactants atoms. set unique numbers.
                    nm += 1
                    p_mapping[n] = nm

            for m in self.products:
                m.remap(p_mapping)
            self.flush_cache()
            return True
        return False

    def fix_mapping(self: Union['ReactionContainer', 'Mapping'], *, logging: bool = False) -> \
            Union[bool, List[Union[Tuple[str, List[int]], Tuple[int, str, str, Tuple[int, ...]]]]]:
        """
        Fix atom-to-atom mapping of some functional groups. Return True if found AAM errors.
        """
        log = []
        seen = set()
        if not self:
            return False

        for r_pattern, p_pattern, _map, fix in rules:
            found = []
            for m in self.reactants:
                for mapping in r_pattern.get_mapping(m, automorphism_filter=False):
                    if mapping[1] not in seen:
                        found.append((m.atom(mapping[1]).atomic_number, mapping[1],  # get matched first atom
                                      # prepare expected in product mapping
                                      {(v, mapping[k]) for k, v in _map.items()}))

            if not found:
                continue
            for m in self.products:
                for mapping in p_pattern.get_mapping(m, automorphism_filter=False):
                    atom = mapping[1]
                    if atom in seen:
                        continue
                    an = m.atom(atom).atomic_number
                    for i, (a, n, k) in enumerate(found):
                        if atom == n and a == an and k.issubset(mapping.items()):
                            break
                    else:
                        continue

                    del found[i]
                    v = {mapping[k]: mapping[v] for k, v in fix.items()}
                    m.remap(v)
                    seen.add(atom)
                    log.append(('group remap', list(v.values())))
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

    @class_cached_property
    def __attention_model(self):
        from chython import torch_device
        from chytorch.nn import ReactionEncoder
        from torch import load

        model = ReactionEncoder(d_model=1024, dim_feedforward=3072)
        model.load_state_dict(load(resource_stream(__package__, 'mapping.pt'), map_location=torch_device))
        model.eval()
        return model

    def __prepare_remapping(self: 'ReactionContainer'):
        r_map = [n for m in self.reactants for n in m]
        p_map = [n for m in self.products for n in m]
        ra = len(r_map)  # number of reactants atoms
        pa = len(p_map)  # number of products atoms

        ram = [False]  # reactants atoms mask
        r_atoms = []
        r_adj = zeros((ra, ra), dtype=bool)
        i = 0
        for m in self.reactants:
            ram.append(False)
            ram.extend(repeat(True, len(m)))
            a = m.adjacency_matrix()
            j = i + len(m)
            r_adj[i:j, i:j] = a
            i = j
            r_atoms.extend(a.atomic_number for _, a in m.atoms())
        r_atoms = array(r_atoms, dtype=int)

        pam = [False] * len(ram)  # products atoms mask
        p_atoms = []
        p_adj = zeros((pa, pa), dtype=bool)
        i = 0
        for m in self.products:
            pam.append(False)
            pam.extend(repeat(True, len(m)))
            a = m.adjacency_matrix()
            j = i + len(m)
            p_adj[i:j, i:j] = a
            i = j
            p_atoms.extend(a.atomic_number for _, a in m.atoms())
        p_atoms = array(p_atoms, dtype=int)

        ram.extend(repeat(False, len(pam) - len(ram)))
        ram = array(ram, dtype=bool)
        pam = array(pam, dtype=bool)
        return r_map, p_map, r_atoms, p_atoms, r_adj, p_adj, ram, pam, ra, pa

    def __get_attention(self):
        from chython import torch_device
        from chytorch.utils.data import ReactionDataset
        from torch import no_grad

        converter = ReactionDataset((self,), (None,), property_type=bool)
        a, n, i, e, r, _ = converter.collate((converter[0],))
        a.to(torch_device)
        n.to(torch_device)
        i.to(torch_device)
        e.to(torch_device)
        r.to(torch_device)

        with no_grad():
            am = self.__attention_model(a, n, i, e, r, True)[1].squeeze(0).cpu().numpy()  # raw attention matrix
        return am


__all__ = ['Mapping']
