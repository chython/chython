# -*- coding: utf-8 -*-
#
#  Copyright 2020-2022 Ramil Nugmanov <nougmanoff@protonmail.com>
#  Copyright 2020 Nail Samikaev <samikaevn@yandex.ru>
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
from collections import deque, defaultdict
from functools import cached_property
from itertools import product, chain, repeat, combinations
from typing import TYPE_CHECKING, Iterator, Union, List
from ._acid import rules as acid_rules, stripped_rules as stripped_acid_rules
from ._base import rules as base_rules, stripped_rules as stripped_base_rules
from ._keto_enol import *
from ..aromatics.kekule import _kekule_component
from ...exceptions import InvalidAromaticRing


if TYPE_CHECKING:
    from chython import MoleculeContainer


class Tautomers:
    """
    Oxides and sulphides ignored.
    """
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

        self._charges.update(mol._charges)
        self._hydrogens.update(mol._hydrogens)
        self.flush_cache()
        if _fix_stereo:
            self.fix_stereo()
        if logging:
            return list(changed)
        return True

    def enumerate_tautomers(self: Union['MoleculeContainer', 'Tautomers'], *, prepare_molecules=True, zwitter=True,
                            partial=False, increase_aromaticity=True, keep_sugars=True,
                            limit: int = 1000) -> Iterator['MoleculeContainer']:
        """
        Enumerate all possible tautomeric forms of molecule.

        :param prepare_molecules: Standardize structures for correct processing.
        :param zwitter: Do zwitter-ions enumeration.
        :param partial: Allow OC=CC=C>>O=CCC=C or O=CC=CC>>OC=C=CC
        :param increase_aromaticity: prevent aromatic ring destruction
        :param keep_sugars: prevent carbonyl moving in sugars
        :param limit: Maximum amount of generated structures.
        """
        if limit < 3:
            raise ValueError('limit should be greater or equal 3')

        has_stereo = bool(self._atoms_stereo or self._allenes_stereo or self._cis_trans_stereo)
        counter = 1

        copy = self.copy()
        copy.clean_stereo()
        # sssr, neighbors and heteroatoms are same for all tautomers.
        # prevent recalculation by sharing cache.
        self.__set_cache(copy)
        if prepare_molecules:  # transform to kekule form without hydrogens
            k = copy.kekule()
            i = copy.implicify_hydrogens(_fix_stereo=False)
            if k or i:  # reset cache after flush
                self.__set_cache(copy)

        thiele = copy.copy()  # transform to thiele to prevent duplicates and dearomatization
        self.__set_cache(thiele)
        if thiele.thiele(fix_tautomers=False):
            self.__set_cache(thiele)

        # return origin structure as first tautomer
        if has_stereo:
            yield self.__set_stereo(thiele.copy())
        else:
            yield thiele

        seen = {thiele: None}  # value is parent molecule - required for preventing migrations in sugars.

        # first try to neutralize
        if copy.neutralize(_fix_stereo=False):  # found neutral form
            thiele = copy.copy()
            self.__set_cache(copy)  # restore cache
            self.__set_cache(thiele)
            if thiele.thiele(fix_tautomers=False):
                self.__set_cache(thiele)

            if has_stereo:
                yield self.__set_stereo(thiele.copy())
            else:
                yield thiele
            counter += 1
            seen[thiele] = None

        # lets iteratively do keto-enol transformations.
        rings_count = len(thiele.aromatic_rings)  # increase rings strategy.
        queue = deque([(copy, thiele)])
        new_queue = [thiele]  # new_queue - molecules suitable for hetero-arenes enumeration.
        # store aromatic form to seen. kekule forms not suitable for duplicate checking.

        while queue:
            current, thiele_current = queue.popleft()
            for mol, ket in current._enumerate_keto_enol_tautomers(partial):
                thiele = mol.copy()
                self.__set_cache(mol)
                self.__set_cache(thiele)
                if thiele.thiele(fix_tautomers=False):  # reset cache after flush_cache.
                    self.__set_cache(thiele)

                if thiele not in seen:
                    seen[thiele] = current
                    rc = len(thiele.aromatic_rings)
                    if increase_aromaticity:
                        if rc < rings_count:  # skip aromatic rings destruction
                            continue
                        elif rc > rings_count:  # higher aromaticity found. flush old queues.
                            rings_count = rc
                            queue = deque([(mol, thiele)])
                            new_queue = [thiele]
                            copy = mol  # new entry point.
                            if has_stereo:
                                yield self.__set_stereo(thiele.copy())
                            else:
                                yield thiele
                            counter += 1
                            if counter == limit:
                                return
                            break
                    if keep_sugars and current is not copy and ket:
                        # prevent carbonyl migration in sugars. skip entry point.
                        # search alpha hydroxy ketone inversion
                        before = seen[thiele_current]._sugar_groups
                        if any((k, e) in before for e, k in mol._sugar_groups):
                            continue

                    queue.append((mol, thiele))
                    new_queue.append(thiele)
                    if has_stereo:
                        yield self.__set_stereo(thiele.copy())
                    else:
                        yield thiele
                    counter += 1
                    if counter == limit:
                        return

        queue = deque(new_queue)
        while queue:
            current = queue.popleft()
            for mol in current._enumerate_hetero_arene_tautomers():
                self.__set_cache(mol)
                if mol not in seen:
                    seen[mol] = None
                    queue.append(mol)
                    new_queue.append(mol)  # new hetero-arenes also should be included to this list.
                    if has_stereo:
                        yield self.__set_stereo(mol.copy())
                    else:
                        yield mol
                    counter += 1
                    if counter == limit:
                        return

        if not zwitter:
            return
        queue = deque(new_queue)
        while queue:
            current = queue.popleft()
            for mol in current._enumerate_zwitter_tautomers():
                self.__set_cache(mol)
                if mol not in seen:
                    seen[mol] = None
                    queue.append(mol)
                    if has_stereo:
                        yield self.__set_stereo(mol.copy())
                    else:
                        yield mol
                    counter += 1
                    if counter == limit:
                        return

    def __set_cache(self: 'MoleculeContainer', mol):
        try:
            neighbors = self.__dict__['__cached_args_method_neighbors']
        except KeyError:
            neighbors = self.__dict__['__cached_args_method_neighbors'] = {}
        try:
            heteroatoms = self.__dict__['__cached_args_method_heteroatoms']
        except KeyError:
            heteroatoms = self.__dict__['__cached_args_method_heteroatoms'] = {}
        try:
            is_ring_bond = self.__dict__['__cached_args_method_is_ring_bond']
        except KeyError:
            is_ring_bond = self.__dict__['__cached_args_method_is_ring_bond'] = {}

        mol.__dict__['sssr'] = self.sssr  # thiele/kekule
        mol.__dict__['ring_atoms'] = self.ring_atoms  # morgan
        mol.__dict__['_connected_components'] = self._connected_components  # isomorphism
        mol.__dict__['atoms_rings_sizes'] = self.atoms_rings_sizes  # isomorphism
        mol.__dict__['__cached_args_method_neighbors'] = neighbors  # isomorphism
        mol.__dict__['__cached_args_method_heteroatoms'] = heteroatoms  # isomorphism
        mol.__dict__['__cached_args_method_is_ring_bond'] = is_ring_bond  # isomorphism

    def __set_stereo(self: 'MoleculeContainer', mol):
        mol._atoms_stereo.update(self._atoms_stereo)
        mol._allenes_stereo.update(self._allenes_stereo)
        mol._cis_trans_stereo.update(self._cis_trans_stereo)
        mol.fix_stereo()
        return mol

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
                copy = self.copy()
                for a in acceptors:
                    copy._hydrogens[a] += 1
                    copy._charges[a] += 1
                for c in combinations(donors, len(acceptors)):
                    mol = copy.copy()
                    for d in c:
                        mol._hydrogens[d] -= 1
                        mol._charges[d] -= 1
                    yield mol, acceptors.union(c)
            elif len(donors) < len(acceptors):
                copy = self.copy()
                for d in donors:
                    copy._hydrogens[d] -= 1
                    copy._charges[d] -= 1
                for c in combinations(acceptors, len(donors)):
                    mol = copy.copy()
                    for a in c:
                        mol._hydrogens[a] += 1
                        mol._charges[a] += 1
                    yield mol, donors.union(c)
            else:  # balanced!
                mol = self.copy()
                for d in donors:
                    mol._hydrogens[d] -= 1
                    mol._charges[d] -= 1
                for a in acceptors:
                    mol._hydrogens[a] += 1
                    mol._charges[a] += 1
                yield mol, donors | acceptors
        elif donors or acceptors:
            mol = self.copy()
            for d in donors:
                mol._hydrogens[d] -= 1
                mol._charges[d] -= 1
            for a in acceptors:
                mol._hydrogens[a] += 1
                mol._charges[a] += 1
            yield mol, donors | acceptors

    def _enumerate_keto_enol_tautomers(self: Union['MoleculeContainer', 'Tautomers'], partial=False):
        for fix, ket in self.__enumerate_bonds(partial):
            if ket:
                a = fix[-1][1]
                d = fix[0][0]
            else:
                a = fix[0][0]
                d = fix[-1][1]

            mol = self.copy()
            m_bonds = mol._bonds
            for n, m, b in fix:
                m_bonds[n][m]._Bond__order = b

            mol._hydrogens[a] += 1
            mol._hydrogens[d] -= 1
            yield mol, ket

    def _enumerate_zwitter_tautomers(self: 'MoleculeContainer'):
        donors = []
        acceptors = []
        for q in acid_rules:
            for mapping in q.get_mapping(self, automorphism_filter=False):
                donors.append(mapping[1])
        for q in base_rules:
            for mapping in q.get_mapping(self, automorphism_filter=False):
                acceptors.append(mapping[1])

        for d, a in product(donors, acceptors):
            mol = self.copy()
            mol._hydrogens[d] -= 1
            mol._hydrogens[a] += 1
            mol._charges[d] -= 1
            mol._charges[a] += 1
            yield mol

    def _enumerate_hetero_arene_tautomers(self: 'MoleculeContainer'):
        atoms = self._atoms
        bonds = self._bonds
        hydrogens = self._hydrogens
        charges = self._charges
        radicals = self._radicals

        rings = defaultdict(list)  # aromatic skeleton
        for n, m_bond in bonds.items():
            for m, bond in m_bond.items():
                if bond.order == 4:
                    rings[n].append(m)
        if not rings:
            return

        acceptors = set()
        donors = set()
        single_bonded = set()
        for n, ms in rings.items():
            if len(ms) == 2:
                if atoms[n].atomic_number in (5, 7, 15) and not charges[n] and not radicals[n]:
                    # only neutral B, N, P
                    if hydrogens[n]:  # pyrrole
                        donors.add(n)
                    elif len(bonds[n]) == 2:  # pyridine
                        acceptors.add(n)
                    else:
                        single_bonded.add(n)
                elif charges[n] == -1 and atoms[n].atomic_number == 6:  # ferrocene
                    single_bonded.add(n)
            elif len(ms) == 3 and atoms[n].atomic_number in (5, 7, 15) and not charges[n] and not radicals[n]:
                single_bonded.add(n)
        if not donors or not acceptors:
            return

        atoms = set(rings)
        components = []
        while atoms:
            start = atoms.pop()
            component = {start: rings[start]}
            queue = deque([start])
            while queue:
                current = queue.popleft()
                for n in rings[current]:
                    if n not in component:
                        queue.append(n)
                        component[n] = rings[n]

            atoms.difference_update(component)
            if donors.isdisjoint(component) or acceptors.isdisjoint(component):
                continue
            components.append(component)

        if not components:
            return
        for component in components:
            for d, a in product(component.keys() & donors, component.keys() & acceptors):
                sb = component.keys() & single_bonded
                sb.add(a)  # now pyrrole
                try:
                    next(_kekule_component(component, sb, (), 0))
                except InvalidAromaticRing:
                    continue
                mol = self.copy()
                mol._hydrogens[d] = 0
                mol._hydrogens[a] = 1
                yield mol

    @cached_property
    def _sugar_groups(self):
        ek = []
        for mapping in sugar_group.get_mapping(self, automorphism_filter=False):
            e, k = mapping[1], mapping[2]
            ek.append((e, k))
        return ek

    def __enumerate_bonds(self: 'MoleculeContainer', partial):
        atoms = self._atoms
        bonds = self._bonds
        hydrogens = self._hydrogens
        hybridization = self.hybridization
        rings = self.atoms_rings_sizes

        # search neutral oxygen and nitrogen
        donors = defaultdict(set)
        acceptors = defaultdict(set)
        for q in enol_rules:
            for mapping in q.get_mapping(self, automorphism_filter=False):
                donors[mapping[1]].add(mapping[2])
        for q in keto_rules:
            for mapping in q.get_mapping(self, automorphism_filter=False):
                acceptors[mapping[1]].add(mapping[2])

        for (atom, dirs), hydrogen, anti in chain(zip(donors.items(), repeat(True), repeat(acceptors)),
                                                  zip(acceptors.items(), repeat(False), repeat(donors))):
            path = []
            seen = {atom}
            stack = [(atom, n, 2 if hydrogen else 1, 0) for n in dirs]
            while stack:
                last, current, bond, depth = stack.pop()

                if partial and path and not len(path) % 2 and \
                        (hydrogen or  # enol > ketone
                         hydrogens[(x := path[-1][1])] and (x not in rings or all(x > 7 for x in rings[x]))):  # ketone>
                    # return partial hops. ignore allenes in small rings.
                    yield path, hydrogen
                if len(path) > depth:  # fork found
                    if not partial and not len(path) % 2 and (hydrogen or hydrogens[path[-1][1]]):
                        # end of path found. return it and start new one.
                        yield path, hydrogen
                    seen.difference_update(x for _, x, _ in path[depth:])
                    path = path[:depth]

                path.append((last, current, bond))

                # adding neighbors
                depth += 1
                seen.add(current)
                if bond == 2:
                    next_bond = 1
                else:
                    next_bond = 2

                for n, b in bonds[current].items():
                    if n == last:
                        continue
                    elif n in seen:  # aromatic ring destruction. pyridine double bonds shift
                        continue
                    elif n in anti:  # enol-ketone switch
                        if current in anti[n]:
                            if hydrogens:
                                if b.order == 2:
                                    cp = path.copy()
                                    cp.append((current, n, 1))
                                    yield cp, True
                            elif b.order == 1:
                                cp = path.copy()
                                cp.append((current, n, 2))
                                yield cp, False
                    elif b.order == bond and atoms[n].atomic_number == 6:  # classic keto-enol route
                        hb = hybridization(n)
                        if hb == 2:  # grow up
                            stack.append((current, n, next_bond, depth))
                        elif hydrogen:
                            if hb == 3:  # OC=CC=C=C case
                                cp = path.copy()
                                cp.append((current, n, 1))
                                yield cp, True  # ketone found
                        elif hb == 1 and hydrogens[n]:  # ketone >> enol
                            cp = path.copy()
                            cp.append((current, n, 2))
                            yield cp, False

            if path and not len(path) % 2 and \
                    (hydrogen or  # enol > ketone
                     hydrogens[(x := path[-1][1])] and (x not in rings or all(x > 7 for x in rings[x]))):
                yield path, hydrogen


__all__ = ['Tautomers']
