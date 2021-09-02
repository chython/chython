# -*- coding: utf-8 -*-
#
#  Copyright 2020, 2021 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from lazy_object_proxy import Proxy
from typing import TYPE_CHECKING, Iterator, Union, List
from ._acid import rules as acid_rules, stripped_rules as stripped_acid_rules
from ._base import rules as base_rules, stripped_rules as stripped_base_rules
from ...exceptions import InvalidAromaticRing
from ...periodictable import ListElement


if TYPE_CHECKING:
    from chython import MoleculeContainer


def _sugar_group():
    from ...containers import QueryContainer
    q = QueryContainer()
    q.add_atom(ListElement(['O', 'N']), hydrogens=(1, 2))  # enol
    q.add_atom(ListElement(['O', 'N']))  # ketone
    q.add_atom('C', hybridization=1)
    q.add_atom('C', hybridization=2)
    q.add_bond(1, 3, 1)
    q.add_bond(3, 4, 1)
    q.add_bond(2, 4, 2)
    return q


sugar_group = Proxy(_sugar_group)


class Tautomers:
    __slots__ = ()

    def neutralize(self: 'MoleculeContainer', *, fix_stereo=True, logging=False) -> Union[bool, List[int]]:
        """
        Convert organic salts to neutral form if possible. Only one possible form used for charge unbalanced structures.

        :param logging: return changed atoms list
        """
        try:
            mol, changed = next(self._neutralize())
        except StopIteration:
            if logging:
                return []
            return False

        self._charges.update(mol._charges)
        self._hydrogens.update(mol._hydrogens)
        self.flush_cache()
        if fix_stereo:
            self.fix_stereo()
        if logging:
            return changed
        return True

    def tautomerize(self: 'MoleculeContainer', *, prepare_molecules=True, limit: int = 1000) -> bool:
        """
        Convert structure to canonical tautomeric form. Return True if structure changed.
        """
        def key(m):
            # more aromatic rings is better
            # less charged atoms is better
            # lower Huckel energy is better
            # smiles alphabetically sorted (descent). bad rule!
            return len(m.aromatic_rings), -sum(x != 0 for x in m._charges.values()), -m.huckel_energy, str(m)

        # do zwitter-ions enumeration for charged structure - search canonic ion.
        canon = max(self.enumerate_tautomers(prepare_molecules=prepare_molecules, full=bool(int(self)), limit=limit),
                    key=key)
        if canon != self:  # attach state of canonic tautomer to self
            # atoms, radicals state, parsed_mapping and plane are unchanged
            self._bonds = canon._bonds
            self._charges = canon._charges  # for zwitter-ionic tautomers
            self._hydrogens = canon._hydrogens
            self._atoms_stereo = canon._atoms_stereo
            self._allenes_stereo = canon._allenes_stereo
            self._cis_trans_stereo = canon._cis_trans_stereo
            self._conformers.clear()  # flush 3d
            self.flush_cache()
            return True
        return False

    def enumerate_tautomers(self: 'MoleculeContainer', *, prepare_molecules=True, full=True, limit: int = 1000) -> \
            Iterator['MoleculeContainer']:
        """
        Enumerate all possible tautomeric forms of molecule.

        :param prepare_molecules: Standardize structures for correct processing.
        :param full: Do zwitter-ions enumeration.
        :param limit: Maximum amount of generated structures.
        """
        if limit < 3:
            raise ValueError('limit should be greater or equal 3')

        atoms_stereo = self._atoms_stereo
        allenes_stereo = self._allenes_stereo
        cis_trans_stereo = self._cis_trans_stereo

        # sssr, neighbors and heteroatoms are same for all tautomers.
        # prevent recalculation by sharing cache.
        if '__cached_args_method_neighbors' in self.__dict__:
            neighbors = self.__dict__['__cached_args_method_neighbors']
        else:
            neighbors = self.__dict__['__cached_args_method_neighbors'] = {}
        if '__cached_args_method_heteroatoms' in self.__dict__:
            heteroatoms = self.__dict__['__cached_args_method_heteroatoms']
        else:
            heteroatoms = self.__dict__['__cached_args_method_heteroatoms'] = {}

        has_stereo = bool(atoms_stereo or allenes_stereo or cis_trans_stereo)
        counter = 1

        copy = self.copy()
        copy.clean_stereo()
        # share cache
        copy.__dict__['sssr'] = self.sssr  # thiele
        copy.__dict__['ring_atoms'] = self.ring_atoms  # morgan
        copy.__dict__['connected_components'] = self.connected_components  # isomorphism
        copy.__dict__['atoms_rings_sizes'] = self.atoms_rings_sizes  # isomorphism
        copy.__dict__['__cached_args_method_neighbors'] = neighbors  # isomorphism
        copy.__dict__['__cached_args_method_heteroatoms'] = heteroatoms  # isomorphism
        if prepare_molecules:
            k = copy.kekule()
            i = copy.implicify_hydrogens(fix_stereo=False)
            if k or i:  # reset after flush
                copy.__dict__['sssr'] = self.sssr
                copy.__dict__['ring_atoms'] = self.ring_atoms
                copy.__dict__['connected_components'] = self.connected_components
                copy.__dict__['atoms_rings_sizes'] = self.atoms_rings_sizes
                copy.__dict__['__cached_args_method_neighbors'] = neighbors
                copy.__dict__['__cached_args_method_heteroatoms'] = heteroatoms

        thiele = copy.copy()
        thiele.__dict__['sssr'] = self.sssr
        thiele.__dict__['__cached_args_method_neighbors'] = neighbors
        if thiele.thiele(fix_tautomers=False):
            thiele.__dict__['sssr'] = self.sssr
            thiele.__dict__['__cached_args_method_neighbors'] = neighbors  # isomorphism
        thiele.__dict__['__cached_args_method_heteroatoms'] = heteroatoms  # isomorphism
        thiele.__dict__['ring_atoms'] = self.ring_atoms  # morgan
        thiele.__dict__['connected_components'] = self.connected_components  # isomorphism
        thiele.__dict__['atoms_rings_sizes'] = self.atoms_rings_sizes  # isomorphism

        if has_stereo:
            out = thiele.copy()
            out._atoms_stereo.update(atoms_stereo)
            out._allenes_stereo.update(allenes_stereo)
            out._cis_trans_stereo.update(cis_trans_stereo)
            out.fix_stereo()
            yield out  # no cache
        else:
            yield thiele

        seen = {thiele: None}  # value is parent molecule - required for preventing migrations in sugars.

        # first of all try to neutralize
        if copy.neutralize(fix_stereo=False):
            thiele = copy.copy()

            copy.__dict__['sssr'] = self.sssr
            copy.__dict__['ring_atoms'] = self.ring_atoms
            copy.__dict__['connected_components'] = self.connected_components
            copy.__dict__['atoms_rings_sizes'] = self.atoms_rings_sizes
            copy.__dict__['__cached_args_method_neighbors'] = neighbors
            copy.__dict__['__cached_args_method_heteroatoms'] = heteroatoms

            thiele.__dict__['sssr'] = self.sssr
            thiele.__dict__['__cached_args_method_neighbors'] = neighbors
            if thiele.thiele(fix_tautomers=False):
                thiele.__dict__['sssr'] = self.sssr
                thiele.__dict__['__cached_args_method_neighbors'] = neighbors
            thiele.__dict__['__cached_args_method_heteroatoms'] = heteroatoms
            thiele.__dict__['ring_atoms'] = self.ring_atoms
            thiele.__dict__['connected_components'] = self.connected_components
            thiele.__dict__['atoms_rings_sizes'] = self.atoms_rings_sizes

            if has_stereo:
                out = thiele.copy()
                out._atoms_stereo.update(atoms_stereo)
                out._allenes_stereo.update(allenes_stereo)
                out._cis_trans_stereo.update(cis_trans_stereo)
                out.fix_stereo()
                yield out
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
            for mol, ket in current._enumerate_keto_enol_tautomers():
                thiele = mol.copy()
                # set cache to child structures.
                mol.__dict__['sssr'] = self.sssr
                mol.__dict__['ring_atoms'] = self.ring_atoms
                mol.__dict__['connected_components'] = self.connected_components
                mol.__dict__['atoms_rings_sizes'] = self.atoms_rings_sizes
                mol.__dict__['__cached_args_method_neighbors'] = neighbors
                mol.__dict__['__cached_args_method_heteroatoms'] = heteroatoms

                thiele.__dict__['sssr'] = self.sssr
                thiele.__dict__['__cached_args_method_neighbors'] = neighbors
                if thiele.thiele(fix_tautomers=False):  # reset cache after flush_cache.
                    thiele.__dict__['sssr'] = self.sssr
                    thiele.__dict__['__cached_args_method_neighbors'] = neighbors
                thiele.__dict__['__cached_args_method_heteroatoms'] = heteroatoms
                thiele.__dict__['ring_atoms'] = self.ring_atoms
                thiele.__dict__['connected_components'] = self.connected_components
                thiele.__dict__['atoms_rings_sizes'] = self.atoms_rings_sizes

                if thiele not in seen:
                    seen[thiele] = current
                    rc = len(thiele.aromatic_rings)
                    if rc < rings_count:  # skip aromatic rings destruction
                        continue
                    elif rc > rings_count:  # higher aromaticity found. flush old queues.
                        rings_count = rc
                        queue = deque([(mol, thiele)])
                        new_queue = [thiele]
                        copy = mol  # new entry point.
                        if has_stereo:
                            thiele = thiele.copy()
                            thiele._atoms_stereo.update(atoms_stereo)
                            thiele._allenes_stereo.update(allenes_stereo)
                            thiele._cis_trans_stereo.update(cis_trans_stereo)
                            thiele.fix_stereo()
                        yield thiele
                        break
                    elif current is not copy and ket:  # prevent carbonyl migration in sugars. skip entry point.
                        # search alpha hydroxy ketone inversion
                        before = seen[thiele_current]._sugar_groups
                        if any((k, e) in before for e, k in mol._sugar_groups):
                            continue

                    queue.append((mol, thiele))
                    new_queue.append(thiele)
                    if has_stereo:
                        thiele = thiele.copy()
                        thiele._atoms_stereo.update(atoms_stereo)
                        thiele._allenes_stereo.update(allenes_stereo)
                        thiele._cis_trans_stereo.update(cis_trans_stereo)
                        thiele.fix_stereo()
                    yield thiele
                    counter += 1
                    if counter == limit:
                        return

        queue = deque(new_queue)
        while queue:
            current = queue.popleft()
            for mol in current._enumerate_hetero_arene_tautomers():
                mol.__dict__['sssr'] = self.sssr
                mol.__dict__['__cached_args_method_neighbors'] = neighbors
                mol.__dict__['__cached_args_method_heteroatoms'] = heteroatoms
                mol.__dict__['ring_atoms'] = self.ring_atoms
                mol.__dict__['connected_components'] = self.connected_components
                mol.__dict__['atoms_rings_sizes'] = self.atoms_rings_sizes
                if mol not in seen:
                    seen[mol] = None
                    queue.append(mol)
                    new_queue.append(mol)  # new hetero-arenes also should be included to this list.
                    if has_stereo:
                        mol = mol.copy()
                        mol._atoms_stereo.update(atoms_stereo)
                        mol._allenes_stereo.update(allenes_stereo)
                        mol._cis_trans_stereo.update(cis_trans_stereo)
                        mol.fix_stereo()
                    yield mol
                    counter += 1
                    if counter == limit:
                        return

        if not full:
            return

        queue = deque(new_queue)
        while queue:
            current = queue.popleft()
            for mol in current._enumerate_zwitter_tautomers():
                mol.__dict__['sssr'] = self.sssr
                mol.__dict__['__cached_args_method_neighbors'] = neighbors
                mol.__dict__['__cached_args_method_heteroatoms'] = heteroatoms
                mol.__dict__['ring_atoms'] = self.ring_atoms
                mol.__dict__['connected_components'] = self.connected_components
                mol.__dict__['atoms_rings_sizes'] = self.atoms_rings_sizes
                if mol not in seen:
                    seen[mol] = None
                    queue.append(mol)
                    if has_stereo:
                        mol = mol.copy()  # prevent seen hashtable breaking
                        mol._atoms_stereo.update(atoms_stereo)
                        mol._allenes_stereo.update(allenes_stereo)
                        mol._cis_trans_stereo.update(cis_trans_stereo)
                        mol.fix_stereo()
                    yield mol
                    counter += 1
                    if counter == limit:
                        return

    def _neutralize(self: 'MoleculeContainer'):
        donors = []
        acceptors = []
        for q in stripped_acid_rules:
            for mapping in q.get_mapping(self, automorphism_filter=False):
                donors.append(mapping[1])
        for q in stripped_base_rules:
            for mapping in q.get_mapping(self, automorphism_filter=False):
                acceptors.append(mapping[1])

        if not donors or not acceptors:
            return
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
                yield mol, acceptors + list(c)
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
                yield mol, donors + list(c)
        else:
            mol = self.copy()
            for d in donors:
                mol._hydrogens[d] -= 1
                mol._charges[d] -= 1
            for a in acceptors:
                mol._hydrogens[a] += 1
                mol._charges[a] += 1
            yield mol, donors + acceptors

    def _enumerate_keto_enol_tautomers(self: Union['MoleculeContainer', 'Tautomers']):
        for fix, ket in self.__enumerate_bonds():
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
                    next(self._kekule_component(component, sb, ()))
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

    def __enumerate_bonds(self: 'MoleculeContainer'):
        atoms = self._atoms
        bonds = self._bonds
        charge = self._charges
        hydrogens = self._hydrogens
        hybridization = self.hybridization
        heteroatoms = self.heteroatoms
        neighbors = self.neighbors

        # search neutral oxygen and nitrogen
        donors = []
        acceptors = []
        for n, a in atoms.items():
            if charge[n] or not neighbors(n) or heteroatoms(n):
                continue
            if a.atomic_number in (7, 8):
                if hybridization(n) == 2:  # imine!
                    acceptors.append(n)
                elif hydrogens[n]:  # amine!
                    donors.append(n)

        for atom, hydrogen in chain(zip(donors, repeat(True)), zip(acceptors, repeat(False))):
            path = []
            seen = {atom}
            # stack is neighbors
            if hydrogen:  # enol
                stack = [(atom, n, 2, 0) for n in bonds[atom] if hybridization(n) == 2]
            else:  # ketone
                stack = [(atom, n, 1, 0) for n, b in bonds[atom].items() if hybridization(n) == 2 and b.order == 2]

            while stack:
                last, current, bond, depth = stack.pop()

                if len(path) > depth:
                    if len(path) % 2 == 0:
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
                        if not stack:
                            path = []
                    elif b.order == bond and atoms[n].atomic_number == 6:
                        hb = hybridization(n)
                        if hb == 2:  # grow up
                            stack.append((current, n, next_bond, depth))
                        elif hydrogen:
                            if hb == 3:  # OC=CC=C=A case
                                cp = path.copy()
                                cp.append((current, n, 1))
                                yield cp, True  # ketone found
                        elif hb == 1 and hydrogens[n]:  # ketone >> enole
                            cp = path.copy()
                            cp.append((current, n, 2))
                            yield cp, False

            if path and len(path) % 2 == 0:  # time to yield
                yield path, hydrogen


__all__ = ['Tautomers']
