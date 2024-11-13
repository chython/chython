# -*- coding: utf-8 -*-
#
#  Copyright 2022-2024 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from itertools import combinations, product
from typing import TYPE_CHECKING, Union, List
from ._acid import rules as acid_rules, stripped_rules as stripped_acid_rules
from ._base import rules as base_rules, stripped_rules as stripped_base_rules


if TYPE_CHECKING:
    from chython import MoleculeContainer


class AcidBase:
    __slots__ = ()

    def neutralize(self: 'MoleculeContainer', *, keep_charge=True, logging=False,
                   _fix_stereo=True) -> Union[bool, List[int]]:
        """
        Convert organic salts to neutral form if possible. Only one possible form used for charge unbalanced structures.

        :param keep_charge: do partial neutralization to keep total charge of molecule.
        :param logging: return changed atoms list.
        """
        try:
            mol, changed = next(self._neutralize(keep_charge))
        except StopIteration:
            if logging:
                return []
            return False

        self._atoms = mol._atoms
        self.flush_cache(keep_sssr=True, keep_components=True)
        if _fix_stereo:
            self.fix_stereo()
        if logging:
            return list(changed)
        return True

    def enumerate_charged_forms(self: 'MoleculeContainer', *, deep: int = 4, limit: int = 1000):
        """
        Enumerate protonated and deprotonated ions. Use on neutralized molecules.

        :param deep: Maximum amount of added or removed protons.
        :param limit: Maximum amount of generated structures.
        """
        if limit < 1:
            raise ValueError('limit should be greater or equal 1')

        donors = set()
        acceptors = set()
        for q in acid_rules:
            for mapping in q.get_mapping(self, automorphism_filter=False):
                donors.add(mapping[1])
        for q in base_rules:
            for mapping in q.get_mapping(self, automorphism_filter=False):
                acceptors.add(mapping[1])

        h_source = [None] * deep + list(donors)
        h_drain = list(acceptors) + [None] * deep

        seen = set()
        seen_combo = set()
        for r in range(min(len(acceptors), deep), 0, -1):  # number of acceptors to protonate
            for ac in combinations(acceptors, r):  # acceptors atoms to protonate
                uniq = set()
                for dc in combinations(h_source, r):  # H donors
                    if dc in uniq:  # filter None-containing duplicates
                        continue
                    uniq.add(dc)
                    seen_combo.add((dc, ac))
                    mol = self.copy(keep_sssr=True, keep_components=True)
                    for n in ac:
                        a = mol._atoms[n]
                        a._implicit_hydrogens += 1
                        a._charge += 1
                    for n in dc:
                        if n is not None:
                            a = mol._atoms[n]
                            a._implicit_hydrogens -= 1
                            a._charge -= 1
                    if mol not in seen:
                        seen.add(mol)
                        yield mol
                    limit -= 1
                    if not limit:
                        return

        for r in range(1, min(len(donors), deep) + 1):
            for dc in combinations(donors, r):
                uniq = set()
                for ac in combinations(h_drain, r):  # H acceptors
                    if ac in uniq:
                        continue
                    uniq.add(ac)
                    if (dc, ac) in seen_combo:
                        continue
                    mol = self.copy(keep_sssr=True, keep_components=True)
                    for n in ac:
                        if n is not None:
                            a = mol._atoms[n]
                            a._implicit_hydrogens += 1
                            a._charge += 1
                    for n in dc:
                        if n is not None:
                            a = mol._atoms[n]
                            a._implicit_hydrogens -= 1
                            a._charge -= 1
                    if mol not in seen:
                        seen.add(mol)
                        yield mol
                    limit -= 1
                    if not limit:
                        return

    def _neutralize(self: 'MoleculeContainer', keep_charge=True):
        donors = set()
        acceptors = set()
        for q in stripped_acid_rules:
            for mapping in q.get_mapping(self, automorphism_filter=False):
                donors.add(mapping[1])
        for q in stripped_base_rules:
            for mapping in q.get_mapping(self, automorphism_filter=False):
                acceptors.add(mapping[1])

        if keep_charge:
            if not donors or not acceptors:
                return  # neutralization impossible
            elif len(donors) > len(acceptors):
                copy = self.copy(keep_sssr=True, keep_components=True)
                for n in acceptors:
                    a = copy._atoms[n]
                    a._implicit_hydrogens += 1
                    a._charge += 1
                for c in combinations(donors, len(acceptors)):
                    mol = copy.copy(keep_sssr=True, keep_components=True)
                    for n in c:
                        a = mol._atoms[n]
                        a._implicit_hydrogens -= 1
                        a._charge -= 1
                    yield mol, acceptors.union(c)
            elif len(donors) < len(acceptors):
                copy = self.copy(keep_sssr=True, keep_components=True)
                for n in donors:
                    a = copy._atoms[n]
                    a._implicit_hydrogens -= 1
                    a._charge -= 1
                for c in combinations(acceptors, len(donors)):
                    mol = copy.copy(keep_sssr=True, keep_components=True)
                    for n in c:
                        a = mol._atoms[n]
                        a._implicit_hydrogens += 1
                        a._charge += 1
                    yield mol, donors.union(c)
            else:  # balanced!
                mol = self.copy(keep_sssr=True, keep_components=True)
                for n in donors:
                    a = mol._atoms[n]
                    a._implicit_hydrogens -= 1
                    a._charge -= 1
                for n in acceptors:
                    a = mol._atoms[n]
                    a._implicit_hydrogens += 1
                    a._charge += 1
                yield mol, donors | acceptors
        elif donors or acceptors:
            mol = self.copy(keep_sssr=True, keep_components=True)
            for n in donors:
                a = mol._atoms[n]
                a._implicit_hydrogens -= 1
                a._charge -= 1
            for n in acceptors:
                a = mol._atoms[n]
                a._implicit_hydrogens += 1
                a._charge += 1
            yield mol, donors | acceptors

    def _enumerate_zwitter_tautomers(self: 'MoleculeContainer'):
        donors = set()
        acceptors = set()
        for q in acid_rules:
            for mapping in q.get_mapping(self, automorphism_filter=False):
                donors.add(mapping[1])
        for q in base_rules:
            for mapping in q.get_mapping(self, automorphism_filter=False):
                acceptors.add(mapping[1])

        for d, a in product(donors, acceptors):
            mol = self.copy(keep_sssr=True, keep_components=True)
            d = mol._atoms[d]
            a = mol._atoms[a]
            d._implicit_hydrogens -= 1
            a._implicit_hydrogens += 1
            d._charge -= 1
            a._charge += 1
            yield mol


__all__ = ['AcidBase']
