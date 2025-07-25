# -*- coding: utf-8 -*-
#
#  Copyright 2021-2025 Ramil Nugmanov <nougmanoff@protonmail.com>
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


# atomic number constants
B = 5
C = 6
N = 7
O = 8
P = 15
S = 16
As = 33
Se = 34
Te = 52


class Kekule:
    __slots__ = ()

    def kekule(self: Union['Kekule', 'MoleculeContainer'], *, buffer_size=7, ignore_pyrrole_hydrogen=False) -> bool:
        """
        Convert structure to kekule form. Return True if found any aromatic ring. Set implicit hydrogen count and
        hybridization marks on atoms.

        Only one of possible double/single bonds positions will be set.
        For enumerate bonds positions use `enumerate_kekule`.

        :param buffer_size: number of attempts of pyridine form searching.
        :param ignore_pyrrole_hydrogen: ignore hydrogen on pyrrole to fix invalid rings like Cn1cc[nH]c1.
        """
        fixed = self.__fix_rings()  # fix bad aromatic rings
        kekule = next(self.__kekule_full(buffer_size, ignore_pyrrole_hydrogen), None)
        if kekule:
            bonds = self._bonds
            atoms = set()
            for n, m, b in kekule:
                bonds[n][m]._order = b
                atoms.add(n)
                atoms.add(m)
            for n in atoms:
                self.calc_implicit(n)
            self.flush_cache(keep_sssr=True, keep_components=True)
            self.calc_labels()
            return True
        return fixed

    def enumerate_kekule(self: Union['Kekule', 'MoleculeContainer'], ignore_pyrrole_hydrogen=False):
        """
        Enumerate all possible kekule forms of molecule.
        """
        self.__fix_rings()  # fix bad aromatic rings
        for form in self.__kekule_full(0, ignore_pyrrole_hydrogen):
            copy = self.copy(keep_sssr=True, keep_components=True)
            bonds = copy._bonds
            atoms = set()
            for n, m, b in form:
                bonds[n][m]._order = b
                atoms.add(n)
                atoms.add(m)
            for n in atoms:
                copy.calc_implicit(n)
            copy.calc_labels()
            yield copy

    def __fix_rings(self: 'MoleculeContainer'):
        atoms = self._atoms
        bonds = self._bonds
        seen = set()
        keep = True
        for q, af, bf, mm in rules:
            for mapping in q.get_mapping(self, automorphism_filter=False):
                match = set(mapping.values())
                if not mm and not match.isdisjoint(seen):  # prevent double patching of atoms
                    continue
                seen.update(match)

                for n, c in af.items():
                    n = mapping[n]
                    atoms[n]._charge = c
                for n, m, b in bf:
                    n = mapping[n]
                    m = mapping[m]
                    bond = bonds[n][m]
                    bond._order = b
                    bond._stereo = None  # prevent potential errors with cis-trans-labelled double-bonds
                    if b == 8:
                        # flush sssr and components cache
                        keep = False
        if seen:
            self.flush_cache(keep_sssr=keep, keep_components=keep)
            self.calc_labels()
            return True
        return False

    def __prepare_rings(self: 'MoleculeContainer', ignore_pyrrole_hydrogen):
        atoms = self._atoms
        bonds = self._bonds

        rings = defaultdict(list)  # aromatic skeleton
        pyrroles = set()

        double_bonded = defaultdict(list)
        triple_bonded = set()
        for n, m_bond in bonds.items():
            for m, bond in m_bond.items():
                if bond == 4:
                    rings[n].append(m)
                elif bond == 2:
                    double_bonded[n].append(m)
                elif bond == 3:
                    triple_bonded.add(n)

        if not rings:
            return rings, pyrroles, set()
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

        # fix invalid smiles: c1ccccc1c2ccccc2 instead of c1ccccc1-c2ccccc2
        seen = set()
        for n, ms in copy_rings.items():
            if ms:
                seen.add(n)
                for m in ms:
                    if m not in seen:
                        rings[n].remove(m)
                        rings[m].remove(n)
                        bonds[n][m]._order = 1

        if any(len(ms) not in (2, 3) for ms in rings.values()):
            raise InvalidAromaticRing('not in ring aromatic bond or hypercondensed rings: '
                                      f'{{{", ".join(str(n) for n, ms in rings.items() if len(ms) not in (2, 3))}}}')

        # get double bonded ring atoms
        double_bonded = {n for n, ms in double_bonded.items() if ms and n in rings}
        if any(len(rings[n]) != 2 for n in double_bonded):  # double bonded never condensed
            raise InvalidAromaticRing('quinone valence error')
        for n in double_bonded:
            if (atom := atoms[n]) == N:
                if atom.charge != 1:
                    raise InvalidAromaticRing('quinone should be charged N atom')
            elif atom not in (C, P, S, As, Se, Te) or atom.charge:
                raise InvalidAromaticRing('quinone should be neutral S, Se, Te, C, P, As atom')

        for n in rings:
            if (atom := atoms[n]) == C:  # carbon
                if atom.charge == 0:
                    if atom.neighbors not in (2, 3):
                        raise InvalidAromaticRing
                elif atom.charge in (-1, 1):
                    if atom.is_radical:
                        if atom.neighbors == 2:
                            double_bonded.add(n)
                        else:
                            raise InvalidAromaticRing
                    elif atom.neighbors == 3:
                        double_bonded.add(n)
                    elif atom.neighbors == 2:  # benzene (an|cat)ion or pyrrole
                        pyrroles.add(n)
                    else:
                        raise InvalidAromaticRing
                else:
                    raise InvalidAromaticRing
            elif atom in (N, P, As):
                if atom.charge == 0:  # pyrrole or pyridine. include radical pyrrole
                    if atom.is_radical:
                        if atom.neighbors != 2:  # only pyrrole radical
                            raise InvalidAromaticRing
                        double_bonded.add(n)
                    elif atom.neighbors == 3:
                        if atom == N:  # pyrrole only possible
                            double_bonded.add(n)
                        else:  # P(III) or P(V)H
                            pyrroles.add(n)
                    elif atom.neighbors == 2:
                        if atom.implicit_hydrogens is None:  # pyrrole or pyridine
                            pyrroles.add(n)
                        elif atom.implicit_hydrogens == 1:  # only pyrrole
                            if ignore_pyrrole_hydrogen:
                                pyrroles.add(n)
                            else:
                                double_bonded.add(n)
                        elif atom.implicit_hydrogens:  # too many hydrogens for aromatic rings
                            raise InvalidAromaticRing
                    elif atom.neighbors != 4 or atom not in (P, As):  # P(V) in ring [P;a](-R1)-R2
                        raise InvalidAromaticRing
                elif atom.charge == -1:  # pyrrole only
                    if atom.neighbors != 2 or atom.is_radical:
                        raise InvalidAromaticRing
                    double_bonded.add(n)
                elif atom.charge != 1:
                    raise InvalidAromaticRing
                elif atom.is_radical:
                    if atom.neighbors != 2:  # not cation-radical pyridine
                        raise InvalidAromaticRing
                elif atom.neighbors == 2:  # pyrrole cation or protonated pyridine
                    pyrroles.add(n)
                elif atom.neighbors != 3:  # not pyridine oxyde
                    raise InvalidAromaticRing
            elif atom == O:  # furan
                if atom.neighbors == 2:
                    if atom.charge == 0:
                        if atom.is_radical:
                            raise InvalidAromaticRing('radical oxygen')
                        double_bonded.add(n)
                    elif atom.charge == 1:
                        if atom.is_radical:  # furan cation-radical
                            double_bonded.add(n)
                        # pyrylium
                    else:
                        raise InvalidAromaticRing('invalid oxygen charge')
                else:
                    raise InvalidAromaticRing('Triple-bonded oxygen')
            elif atom in (S, Se, Te):  # thiophene
                if n not in double_bonded:  # not sulphoxyde nor sulphone
                    if atom.neighbors == 2:
                        if atom.is_radical:
                            if atom.charge == 1:
                                double_bonded.add(n)
                            else:
                                raise InvalidAromaticRing('S, Se, Te cation-radical expected')
                        if atom.charge == 0:
                            double_bonded.add(n)
                        elif atom.charge != 1:
                            raise InvalidAromaticRing('S, Se, Te cation in benzene like ring expected')
                    elif atom.neighbors == 3:
                        if atom.is_radical:
                            if atom.charge:
                                raise InvalidAromaticRing('S, Se, Te ion-radical ring')
                            double_bonded.add(n)
                        elif atom.charge == 1:
                            double_bonded.add(n)
                        elif atom.charge:
                            raise InvalidAromaticRing('S, Se, Te invalid charge ring')
                    else:
                        raise InvalidAromaticRing('S, Se, Te hypervalent ring')
            elif atom == B:
                if atom.charge == 0:
                    if atom.neighbors == 2:
                        if atom.is_radical:  # C=1O[B]OC=1
                            double_bonded.add(n)
                        else:
                            if atom.implicit_hydrogens is None:  # b1ccccc1, C=1OBOC=1 or B1C=CC=N1
                                pyrroles.add(n)
                            elif atom.implicit_hydrogens == 1:  # C=1O[BH]OC=1 or [BH]1C=CC=N1
                                double_bonded.add(n)
                            elif atom.implicit_hydrogens:
                                raise InvalidAromaticRing
                    elif not atom.is_radical:
                        double_bonded.add(n)
                    else:
                        raise InvalidAromaticRing
                elif atom.charge == 1:
                    if atom.neighbors == 2 and not atom.is_radical:
                        double_bonded.add(n)
                    else:
                        raise InvalidAromaticRing
                elif atom.charge == -1:
                    if atom.neighbors == 2:
                        if not atom.is_radical:  # C=1O[B-]OC=1 or [bH-]1ccccc1
                            pyrroles.add(n)
                        # anion-radical is benzene like
                    elif atom.is_radical:  # C=1O[B-*](R)OC=1
                        double_bonded.add(n)
                    else:
                        pyrroles.add(n)
                else:
                    raise InvalidAromaticRing
            else:
                raise InvalidAromaticRing(f'only B, C, N, P, O, S, Se, Te possible, not: {atoms[n].atomic_symbol}')
        return rings, pyrroles, double_bonded

    def __kekule_full(self, buffer_size, ignore_pyrrole_hydrogen):
        rings, pyrroles, double_bonded = self.__prepare_rings(ignore_pyrrole_hydrogen)
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

        for keks in lazy_product(*(_kekule_component(c, double_bonded & c.keys(), pyrroles & c.keys(), buffer_size)
                                   for c in components)):
            yield [x for x in keks for x in x]


def _kekule_component(rings, double_bonded, pyrroles, buffer_size):
    # (current atom, previous atom, bond between cp atoms, path deep for cutting [None if cut impossible])
    stack: List[List[Tuple[int, int, int, Optional[int]]]]
    if double_bonded:  # start from double bonded if exists
        start = next(iter(double_bonded))
        stack = [[(next(iter(rings[start])), start, 1, 0)]]
    else:  # select not pyrrole not condensed atom
        try:
            start = next(n for n, ms in rings.items() if len(ms) == 2 and n not in pyrroles)
        except StopIteration:  # all pyrroles. select not condensed atom.
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
    buffer = []

    while stack:
        atom, prev_atom, bond, _ = stack[-1].pop()
        path.append((atom, prev_atom, bond))
        hashed_path.add(atom)

        if len(path) == size:
            if nether_yielded:
                nether_yielded = False
            if pyrroles and buffer_size:  # prioritize pyridine over pyrrole
                g = defaultdict(int)
                for n, m, b in path:
                    g[n] += b
                    g[m] += b
                # should be pairs of pyrrole atoms
                if sum(b == 2 and n in pyrroles for n, b in g.items()) >= 2:
                    if len(buffer) == buffer_size:  # optimization. try only few times to prevent freezes.
                        buffer_size = 0  # disable bufferization
                        yield from buffer
                        yield path
                        buffer = []
                    else:
                        buffer.append(path)
                else:
                    yield path
                    buffer_size = 0  # disable bufferization
                    if buffer:  # empty buffer
                        yield from buffer
                        buffer = []
            else:
                yield path

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
                    # side-path for storing double bond or atom is quinone or pyrrole
                    if for_stack or atom in double_bonded or atom in pyrroles:
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
            elif len(for_stack) == 1:  # easy path grow. next bond double or include single for pyrroles
                next_atom = for_stack[0]
                if next_atom in double_bonded:  # need double bond, but next atom quinone
                    if atom in pyrroles:
                        stack[-1].append((next_atom, atom, 1, None))
                    else:
                        del stack[-1]
                        if stack:
                            path = path[:stack[-1][-1][-1]]
                            hashed_path = {x for x, *_ in path}
                elif atom in pyrroles:  # try pyrrole and pyridine
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
                    if next_atom2 in double_bonded:
                        if atom in pyrroles:  # shit like O=C1C=CC2=CC=CC3=C2P1C(=O)C=C3
                            stack[-1].append((next_atom1, atom, 1, None))
                            stack[-1].append((next_atom2, atom, 1, None))
                        else:  # bad path
                            del stack[-1]
                            if stack:
                                path = path[:stack[-1][-1][-1]]
                                hashed_path = {x for x, *_ in path}
                    elif atom in pyrroles:  # O=C1C=CC2=CC=CC3=C2P1C=C3 or O=C1C=CC2=CC=CC3=C2P1=CC=C3
                        opposite = stack[-1].copy()
                        opposite.append((next_atom1, atom, 1, None))
                        opposite.append((next_atom2, atom, 2, None))
                        stack[-1].append((next_atom1, atom, 1, None))
                        stack[-1].append((next_atom2, atom, 1, len(path)))
                        stack.append(opposite)  # pyridine first
                    else:  # normal condensed ring
                        stack[-1].append((next_atom1, atom, 1, None))
                        stack[-1].append((next_atom2, atom, 2, None))
                elif next_atom2 in double_bonded:  # quinone next from fork
                    if atom in pyrroles:
                        opposite = stack[-1].copy()
                        opposite.append((next_atom2, atom, 1, None))
                        opposite.append((next_atom1, atom, 2, None))
                        stack[-1].append((next_atom1, atom, 1, None))
                        stack[-1].append((next_atom2, atom, 1, len(path)))
                        stack.append(opposite)
                    else:
                        stack[-1].append((next_atom2, atom, 1, None))
                        stack[-1].append((next_atom1, atom, 2, None))
                elif atom in pyrroles:  # C1=CC2=CC=CC3=C2P1C=C3 or C1=CP2=CC=CC3=C2C1=CC=C3
                    opposite1 = stack[-1].copy()
                    opposite1.append((next_atom2, atom, 1, None))
                    opposite1.append((next_atom1, atom, 2, len(path)))
                    opposite2 = stack[-1].copy()
                    opposite2.append((next_atom1, atom, 1, None))
                    opposite2.append((next_atom2, atom, 2, None))

                    stack[-1].append((next_atom1, atom, 1, None))
                    stack[-1].append((next_atom2, atom, 1, len(path)))
                    stack.append(opposite1)
                    stack.append(opposite2)
                else:  # new path
                    opposite = stack[-1].copy()
                    stack[-1].append((next_atom1, atom, 1, None))
                    stack[-1].append((next_atom2, atom, 2, len(path)))  # double bond on top of stack
                    opposite.append((next_atom2, atom, 1, None))
                    opposite.append((next_atom1, atom, 2, None))
                    stack.append(opposite)
            elif closures and atom not in pyrroles:  # need double bond, but closure should be single bonded
                del stack[-1]
                if stack:
                    path = path[:stack[-1][-1][-1]]
                    hashed_path = {x for x, *_ in path}

    if nether_yielded:
        raise InvalidAromaticRing(f'kekule form not found for: {list(rings)}')
    elif buffer:  # optimal solution not found. return available.
        yield from buffer


__all__ = ['Kekule']
