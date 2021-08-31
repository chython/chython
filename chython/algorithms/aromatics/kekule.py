# -*- coding: utf-8 -*-
#
#  Copyright 2021 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from collections import defaultdict, deque
from typing import List, Optional, Tuple, TYPE_CHECKING, Union
from ._rules import rules
from ..._functions import lazy_product
from ...exceptions import InvalidAromaticRing


if TYPE_CHECKING:
    from chython import MoleculeContainer


class Kekule:
    __slots__ = ()

    def kekule(self: Union['Kekule', 'MoleculeContainer']) -> bool:
        """
        Convert structure to kekule form. Return True if found any aromatic ring. Set implicit hydrogen count and
        hybridization marks on atoms.

        Only one of possible double/single bonds positions will be set.
        For enumerate bonds positions use `enumerate_kekule`.
        """
        kekule = next(self.__kekule_full(), None)
        if kekule:
            bonds = self._bonds
            atoms = set()
            for n, m, b in kekule:
                bonds[n][m]._Bond__order = b
                atoms.add(n)
                atoms.add(m)
            for n in atoms:
                self._calc_implicit(n)
            self.flush_cache()
            return True
        return False

    def enumerate_kekule(self: Union['Kekule', 'MoleculeContainer']):
        """
        Enumerate all possible kekule forms of molecule.
        """
        for form in self.__kekule_full():
            copy = self.copy()
            bonds = copy._bonds
            atoms = set()
            for n, m, b in form:
                bonds[n][m]._Bond__order = b
                atoms.add(n)
                atoms.add(m)
            for n in atoms:
                copy._calc_implicit(n)
            yield copy

    def check_thiele(self, fast=True) -> bool:
        """
        Check basic aromaticity errors of molecule.

        :param fast: don't try to solve kekule form
        """
        try:
            if fast:
                self.__prepare_rings()
            else:
                next(self.__kekule_full(), None)
        except InvalidAromaticRing:
            return False
        return True

    def __fix_rings(self: 'MoleculeContainer'):
        bonds = self._bonds
        seen = set()
        for q, af, bf in rules:
            for mapping in q.get_mapping(self, automorphism_filter=False):
                match = set(mapping.values())
                if not match.isdisjoint(seen):  # prevent double patching of atoms
                    continue
                seen.update(match)

                for n, fix in af.items():
                    n = mapping[n]
                    for key, value in fix.items():
                        getattr(self, key)[n] = value
                for n, m, b in bf:
                    n = mapping[n]
                    m = mapping[m]
                    bonds[n][m]._Bond__order = b
        if seen:
            self.flush_cache()

    def __prepare_rings(self: 'MoleculeContainer'):
        atoms = self._atoms
        charges = self._charges
        radicals = self._radicals
        bonds = self._bonds
        hydrogens = self._hydrogens
        neighbors = self.neighbors

        rings = defaultdict(list)  # aromatic skeleton
        pyroles = set()

        double_bonded = defaultdict(list)
        triple_bonded = set()
        for n, m_bond in bonds.items():
            for m, bond in m_bond.items():
                bo = bond.order
                if bo == 4:
                    rings[n].append(m)
                elif bo == 2:
                    double_bonded[n].append(m)
                elif bo == 3:
                    triple_bonded.add(n)

        if not rings:
            return rings, pyroles, set()
        elif not triple_bonded.isdisjoint(rings):
            raise InvalidAromaticRing('triple bonds connected to rings')

        copy_rings = {n: ms.copy() for n, ms in rings.items()}
        for r in self.sssr:
            if set(r).issubset(rings):
                n, *_, m = r
                if n not in rings[m]:  # fix invalid structures: c1ccc-cc1
                    # remove inner ring double bonds: c1ccc=cc1
                    if n in double_bonded and m in double_bonded and m in double_bonded[n]:
                        double_bonded[n].remove(m)
                        double_bonded[m].remove(n)
                    rings[m].append(n)
                    rings[n].append(m)
                elif m in copy_rings[n]:
                    copy_rings[n].remove(m)
                    copy_rings[m].remove(n)
                for n, m in zip(r, r[1:]):
                    if n not in rings[m]:
                        if n in double_bonded and m in double_bonded and m in double_bonded[n]:
                            double_bonded[n].remove(m)
                            double_bonded[m].remove(n)
                        rings[m].append(n)
                        rings[n].append(m)
                    elif m in copy_rings[n]:
                        copy_rings[n].remove(m)
                        copy_rings[m].remove(n)

        if any(len(ms) not in (2, 3) for ms in rings.values()):
            raise InvalidAromaticRing('not in ring aromatic bond or hypercondensed rings: '
                                      f'{{{", ".join(str(n) for n, ms in rings.items() if len(ms) not in (2, 3))}}}')

        # fix invalid smiles: c1ccccc1c2ccccc2 instead of c1ccccc1-c2ccccc2
        seen = set()
        for n, ms in copy_rings.items():
            if ms:
                seen.add(n)
                for m in ms:
                    if m not in seen:
                        rings[n].remove(m)
                        rings[m].remove(n)
                        bonds[n][m]._Bond__order = 1

        # get double bonded ring atoms
        double_bonded = {n for n, ms in double_bonded.items() if ms and n in rings}
        if any(len(rings[n]) != 2 for n in double_bonded):  # double bonded never condensed
            raise InvalidAromaticRing('quinone valence error')
        for n in double_bonded:
            if atoms[n].atomic_number == 7:
                if charges[n] != 1:
                    raise InvalidAromaticRing('quinone should be charged N atom')
            elif atoms[n].atomic_number not in (6, 15, 16, 33, 34, 52) or charges[n]:
                raise InvalidAromaticRing('quinone should be neutral S, Se, Te, C, P, As atom')

        for n in rings:
            an = atoms[n].atomic_number
            ac = charges[n]
            ab = neighbors(n)
            if an == 6:  # carbon
                if ac == 0:
                    if ab not in (2, 3):
                        raise InvalidAromaticRing
                elif ac in (-1, 1):
                    if radicals[n]:
                        if ab == 2:
                            double_bonded.add(n)
                        else:
                            raise InvalidAromaticRing
                    elif ab == 3:
                        double_bonded.add(n)
                    elif ab == 2:  # benzene an[cat]ion or pyrole
                        pyroles.add(n)
                    else:
                        raise InvalidAromaticRing
                else:
                    raise InvalidAromaticRing
            elif an in (7, 15, 33):
                if ac == 0:  # pyrole or pyridine. include radical pyrole
                    if radicals[n]:
                        if ab != 2:
                            raise InvalidAromaticRing
                        double_bonded.add(n)
                    elif ab == 3:
                        if an == 7:  # pyrole only possible
                            double_bonded.add(n)
                        else:  # P(III) or P(V)H
                            pyroles.add(n)
                    elif ab == 2:
                        ah = hydrogens[n]
                        if ah is None:
                            pyroles.add(n)
                        elif ah == 1:  # only pyrole
                            double_bonded.add(n)
                        elif ah:
                            raise InvalidAromaticRing
                    elif ab != 4 or an not in (15, 33):  # P(V) in ring
                        raise InvalidAromaticRing
                elif ac == -1:  # pyrole only
                    if ab != 2 or radicals[n]:
                        raise InvalidAromaticRing
                    double_bonded.add(n)
                elif ac != 1:
                    raise InvalidAromaticRing
                elif radicals[n]:
                    if ab != 2:  # not cation-radical pyridine
                        raise InvalidAromaticRing
                elif ab == 2:  # pyrole cation or protonated pyridine
                    pyroles.add(n)
                elif ab != 3:  # not pyridine oxyde
                    raise InvalidAromaticRing
            elif an == 8:  # furan
                if ab == 2:
                    if ac == 0:
                        if radicals[n]:
                            raise InvalidAromaticRing('radical oxygen')
                        double_bonded.add(n)
                    elif ac == 1:
                        if radicals[n]:  # furan cation-radical
                            double_bonded.add(n)
                        # pyrylium
                    else:
                        raise InvalidAromaticRing('invalid oxygen charge')
                else:
                    raise InvalidAromaticRing('Triple-bonded oxygen')
            elif an in (16, 34, 52):  # thiophene
                if n not in double_bonded:  # not sulphoxyde or sulphone
                    if ab == 2:
                        if radicals[n]:
                            if ac == 1:
                                double_bonded.add(n)
                            else:
                                raise InvalidAromaticRing('S, Se, Te cation-radical expected')
                        if ac == 0:
                            double_bonded.add(n)
                        elif ac != 1:
                            raise InvalidAromaticRing('S, Se, Te cation in benzene like ring expected')
                    elif ab == 3:
                        if radicals[n]:
                            if ac:
                                raise InvalidAromaticRing('S, Se, Te ion-radical ring')
                            double_bonded.add(n)
                        elif ac == 1:
                            double_bonded.add(n)
                        elif ac:
                            raise InvalidAromaticRing('S, Se, Te invalid charge ring')
                    else:
                        raise InvalidAromaticRing('S, Se, Te hypervalent ring')
            elif an == 5:  # boron
                if ac == 0:
                    if ab == 2:
                        if radicals[n]:  # C=1O[B]OC=1
                            double_bonded.add(n)
                        else:
                            ah = hydrogens[n]
                            if ah is None:  # b1ccccc1, C=1OBOC=1 or B1C=CC=N1
                                pyroles.add(n)
                            elif ah == 1:  # C=1O[BH]OC=1 or [BH]1C=CC=N1
                                double_bonded.add(n)
                            elif ah:
                                raise InvalidAromaticRing
                    elif not radicals[n]:
                        double_bonded.add(n)
                    else:
                        raise InvalidAromaticRing
                elif ac == 1:
                    if ab == 2 and not radicals[n]:
                        double_bonded.add(n)
                    else:
                        raise InvalidAromaticRing
                elif ac == -1:
                    if ab == 2:
                        if not radicals[n]:  # C=1O[B-]OC=1 or [bH-]1ccccc1
                            pyroles.add(n)
                        # anion-radical is benzene like
                    elif radicals[n]:  # C=1O[B-*](R)OC=1
                        double_bonded.add(n)
                    else:
                        pyroles.add(n)
                else:
                    raise InvalidAromaticRing
            else:
                raise InvalidAromaticRing(f'only B, C, N, P, O, S, Se, Te possible, not: {atoms[n].atomic_symbol}')
        return rings, pyroles, double_bonded

    def __kekule_full(self):
        self.__fix_rings()  # fix bad aromatic rings
        rings, pyroles, double_bonded = self.__prepare_rings()
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

            components.append(component)
            atoms.difference_update(component)

        for keks in lazy_product(*(_kekule_component(c, double_bonded & c.keys(), pyroles & c.keys())
                                   for c in components)):
            yield [x for x in keks for x in x]


def _kekule_component(rings, double_bonded, pyroles):
    # (current atom, previous atom, bond between cp atoms, path deep for cutting [None if cut impossible])
    stack: List[List[Tuple[int, int, int, Optional[int]]]]
    if double_bonded:  # start from double bonded if exists
        start = next(iter(double_bonded))
        stack = [[(next(iter(rings[start])), start, 1, 0)]]
    else:  # select not pyrole not condensed atom
        try:
            start = next(n for n, ms in rings.items() if len(ms) == 2 and n not in pyroles)
        except StopIteration:  # all pyroles. select not condensed atom.
            try:
                start = next(n for n, ms in rings.items() if len(ms) == 2)
            except StopIteration:  # fullerene?
                start = next(iter(rings))
                double_bonded.add(start)
                stack = [[(next_atom, start, 2, 0)] for next_atom in rings[start]]
            else:
                stack = [[(next_atom, start, 1, 0)] for next_atom in rings[start]]
        else:
            stack = [[(next_atom, start, 1, 0)] for next_atom in rings[start]]

    size = sum(len(x) for x in rings.values()) // 2
    path = []
    hashed_path = set()
    nether_yielded = True

    while stack:
        atom, prev_atom, bond, _ = stack[-1].pop()
        path.append((atom, prev_atom, bond))
        hashed_path.add(atom)

        if len(path) == size:
            yield path
            if nether_yielded:
                nether_yielded = False
            del stack[-1]
            if stack:
                path = path[:stack[-1][-1][-1]]
                hashed_path = {x for x, *_ in path}
        elif atom != start:
            for_stack = []
            closures = []
            loop = 0
            for next_atom in rings[atom]:
                if next_atom == prev_atom:  # only forward. behind us is the homeland
                    continue
                elif next_atom == start:
                    loop = next_atom
                elif next_atom in hashed_path:  # closure found
                    closures.append(next_atom)
                else:
                    for_stack.append(next_atom)

            if loop:  # we found starting point.
                if bond == 2:  # finish should be single bonded
                    if double_bonded:  # ok
                        stack[-1].insert(0, (loop, atom, 1, None))
                    else:
                        del stack[-1]
                        if stack:
                            path = path[:stack[-1][-1][-1]]
                            hashed_path = {x for x, *_ in path}
                        continue
                elif double_bonded:  # we in quinone ring. finish should be single bonded
                    # side-path for storing double bond or atom is quinone or pyrole
                    if for_stack or atom in double_bonded or atom in pyroles:
                        stack[-1].insert(0, (loop, atom, 1, None))
                    else:
                        del stack[-1]
                        if stack:
                            path = path[:stack[-1][-1][-1]]
                            hashed_path = {x for x, *_ in path}
                        continue
                else:  # finish should be double bonded
                    stack[-1].insert(0, (loop, atom, 2, None))
                    bond = 2  # grow should be single bonded

            if bond == 2 or atom in double_bonded:  # double in - single out. quinone has two single bonds
                for next_atom in closures:
                    path.append((next_atom, atom, 1))  # closures always single-bonded
                    stack[-1].remove((atom, next_atom, 1, None))  # remove fork from stack
                for next_atom in for_stack:
                    stack[-1].append((next_atom, atom, 1, None))
            elif len(for_stack) == 1:  # easy path grow. next bond double or include single for pyroles
                next_atom = for_stack[0]
                if next_atom in double_bonded:  # need double bond, but next atom quinone
                    if atom in pyroles:
                        stack[-1].append((next_atom, atom, 1, None))
                    else:
                        del stack[-1]
                        if stack:
                            path = path[:stack[-1][-1][-1]]
                            hashed_path = {x for x, *_ in path}
                else:
                    if atom in pyroles:  # try pyrole and pyridine
                        opposite = stack[-1].copy()
                        opposite.append((next_atom, atom, 2, None))
                        stack[-1].append((next_atom, atom, 1, len(path)))
                        stack.append(opposite)
                    else:
                        stack[-1].append((next_atom, atom, 2, None))
                        if closures:
                            next_atom = closures[0]
                            path.append((next_atom, atom, 1))  # closures always single-bonded
                            stack[-1].remove((atom, next_atom, 1, None))  # remove fork from stack
            elif for_stack:  # fork
                next_atom1, next_atom2 = for_stack
                if next_atom1 in double_bonded:  # quinone next from fork
                    if next_atom2 in double_bonded:  # bad path
                        del stack[-1]
                        if stack:
                            path = path[:stack[-1][-1][-1]]
                            hashed_path = {x for x, *_ in path}
                    else:
                        stack[-1].append((next_atom1, atom, 1, None))
                        stack[-1].append((next_atom2, atom, 2, None))
                elif next_atom2 in double_bonded:  # quinone next from fork
                    stack[-1].append((next_atom2, atom, 1, None))
                    stack[-1].append((next_atom1, atom, 2, None))
                else:  # new path
                    opposite = stack[-1].copy()
                    stack[-1].append((next_atom1, atom, 1, None))
                    stack[-1].append((next_atom2, atom, 2, len(path)))
                    opposite.append((next_atom2, atom, 1, None))
                    opposite.append((next_atom1, atom, 2, None))
                    stack.append(opposite)
            elif closures and atom not in pyroles:  # need double bond, but closure should be single bonded
                del stack[-1]
                if stack:
                    path = path[:stack[-1][-1][-1]]
                    hashed_path = {x for x, *_ in path}

    if nether_yielded:
        raise InvalidAromaticRing(f'kekule form not found for: {list(rings)}')


__all__ = ['Kekule']
