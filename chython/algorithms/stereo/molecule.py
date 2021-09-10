# -*- coding: utf-8 -*-
#
#  Copyright 2019-2021 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from functools import cached_property
from itertools import combinations
from logging import info
from typing import Dict, Set, Tuple, Union, TYPE_CHECKING
from .graph import Stereo
from ..morgan import _morgan
from ...exceptions import AtomNotFound, IsChiral, NotChiral


if TYPE_CHECKING:
    from chython import MoleculeContainer


def _pyramid_sign(n, u, v, w):
    #
    #  |   n /
    #  |   |\
    #  |   | \
    #  |  /|  \
    #  | / u---v
    #  |/___\_/___
    #        w
    #
    nx, ny, nz = n
    ux, uy, uz = u
    vx, vy, vz = v
    wx, wy, wz = w

    q1x = ux - nx
    q1y = uy - ny
    q1z = uz - nz
    q2x = vx - nx
    q2y = vy - ny
    q2z = vz - nz
    q3x = wx - nx
    q3y = wy - ny
    q3z = wz - nz

    vol = q1x * (q2y * q3z - q2z * q3y) + q1y * (q2z * q3x - q2x * q3z) + q1z * (q2x * q3y - q2y * q3x)
    if vol > 0:
        return 1
    elif vol < 0:
        return -1
    return 0


def _cis_trans_sign(n, u, v, w):
    # n      w
    #  \    /
    #   u--v
    #  /    \
    # x      x
    nx, ny = n
    ux, uy = u
    vx, vy = v
    wx, wy = w

    q1x = ux - nx
    q1y = uy - ny
    q2x = vx - ux
    q2y = vy - uy
    q3x = wx - vx
    q3y = wy - vy

    # cross vectors
    q1q2z = q1x * q2y - q1y * q2x
    q2q3z = q2x * q3y - q2y * q3x

    dot = q1q2z * q2q3z
    if dot > 0:
        return 1
    elif dot < 0:
        return -1
    return 0


def _allene_sign(n, u, v, w):
    # n    w
    # |   /
    # u--v
    nx, ny, nz = n
    ux, uy = u
    vx, vy = v
    wx, wy, wz = w

    q1x = ux - nx
    q1y = uy - ny
    q1z = -nz
    q2x = vx - ux
    q2y = vy - uy
    q3x = wx - vx
    q3y = wy - vy
    q3z = wz

    # cross vectors
    q1q2x = -q1z * q2y
    q1q2y = q1z * q2x
    q1q2z = q1x * q2y - q1y * q2x
    q2q3x = q2y * q3z
    q2q3y = -q2x * q3z
    q2q3z = q2x * q3y - q2y * q3x

    q1q2q3x = q1q2y * q2q3z - q1q2z * q2q3y
    q1q2q3y = q1q2z * q2q3x - q1q2x * q2q3z

    dot = q1q2q3x * q2x + q1q2q3y * q2y
    if dot > 0:
        return 1
    elif dot < 0:
        return -1
    return 0


class MoleculeStereo(Stereo):
    __slots__ = ()

    def add_wedge(self: 'MoleculeContainer', n: int, m: int, mark: bool, *, clean_cache=True):
        """
        Add stereo data by wedge notation of bonds. Use it for tetrahedrons of allenes.

        :param n: number of atom from which wedge bond started
        :param m: number of atom to which wedge bond coming
        :param mark: up bond is True, down is False
        """
        if n not in self._atoms:
            raise AtomNotFound
        if n in self._atoms_stereo:
            raise IsChiral

        plane = self._plane
        if n in self._chiral_tetrahedrons:
            if m not in self._bonds[n]:
                raise AtomNotFound

            if self._atoms[m].atomic_number == 1:
                s = _pyramid_sign((*plane[m], mark), *((*plane[x], 0) for x in self._stereo_tetrahedrons[n]))
            else:
                order = [(*plane[x], mark if x == m else 0) for x in self._stereo_tetrahedrons[n]]
                if len(order) == 3:
                    s = _pyramid_sign((*plane[n], 0), *order)
                else:
                    s = _pyramid_sign(order[-1], *order[:3])
            if s:
                self._atoms_stereo[n] = s > 0
                if clean_cache:
                    self.flush_cache()
        else:
            c = self._stereo_allenes_centers.get(n)
            if c:
                if c in self._allenes_stereo:
                    raise IsChiral
                elif c not in self._chiral_allenes:
                    raise NotChiral

                order = self._stereo_allenes[c]
                t1, t2 = self._stereo_allenes_terminals[c]
                w = order.index(m)
                if w == 0:
                    m1 = order[1]
                    r = False
                elif w == 1:
                    m1 = order[0]
                    t1, t2 = t2, t1
                    r = False
                elif w == 2:
                    m1 = order[1]
                    r = True
                else:
                    m1 = order[0]
                    t1, t2 = t2, t1
                    r = True
                s = _allene_sign((*plane[m], mark), plane[t1], plane[t2], (*plane[m1], 0))
                if s:
                    self._allenes_stereo[c] = s < 0 if r else s > 0
                    if clean_cache:
                        self.flush_cache()
            else:
                # only tetrahedrons and allenes supported
                raise NotChiral

    def calculate_cis_trans_from_2d(self: 'MoleculeContainer', *, clean_cache=True):
        """
        Calculate cis-trans stereo bonds from given 2d coordinates. Unusable for SMILES and INCHI.
        """
        cis_trans_stereo = self._cis_trans_stereo
        plane = self._plane
        flag = False
        while self._chiral_cis_trans:
            stereo = {}
            for nm in self._chiral_cis_trans:
                n, m = nm
                n1, m1, *_ = self._stereo_cis_trans[nm]
                s = _cis_trans_sign(plane[n1], plane[n], plane[m], plane[m1])
                if s:
                    stereo[nm] = s > 0
            if stereo:
                cis_trans_stereo.update(stereo)
                flag = True
                self.flush_stereo_cache()
            else:
                break
        if flag and clean_cache:
            self.flush_cache()

    def add_atom_stereo(self: 'MoleculeContainer', n: int, env: Tuple[int, ...], mark: bool, *, clean_cache=True):
        """
        Add stereo data for specified neighbors bypass. Use it for tetrahedrons of allenes.

        :param n: number of tetrahedron atom or central atom of allene.
        :param env: numbers of atoms with specified bypass
        :param mark: clockwise or anti bypass.

        See <https://www.daylight.com/dayhtml/doc/theory/theory.smiles.html> and <http://opensmiles.org/opensmiles.html>
        """
        if n not in self._atoms:
            raise AtomNotFound
        if n in self._atoms_stereo or n in self._allenes_stereo:
            raise IsChiral
        if not isinstance(mark, bool):
            raise TypeError('stereo mark should be bool')

        if n in self._chiral_tetrahedrons:
            self._atoms_stereo[n] = self._translate_tetrahedron_sign(n, env, mark)
            if clean_cache:
                self.flush_cache()
        elif n in self._chiral_allenes:
            self._allenes_stereo[n] = self._translate_allene_sign(n, *env, mark)
            if clean_cache:
                self.flush_cache()
        else:  # only tetrahedrons supported
            raise NotChiral

    def add_cis_trans_stereo(self: 'MoleculeContainer', n: int, m: int, n1: int, n2: int, mark: bool, *,
                             clean_cache=True):
        """
        Add stereo data to cis-trans double bonds (not allenes).

        n1/n=m/n2

        :param n: number of starting atom of double bonds chain (alkenes of cumulenes)
        :param m: number of ending atom of double bonds chain (alkenes of cumulenes)
        :param n1: number of neighboring atom of starting atom
        :param n2: number of neighboring atom of ending atom
        :param mark: cis or trans

        See <https://www.daylight.com/dayhtml/doc/theory/theory.smiles.html> and <http://opensmiles.org/opensmiles.html
        """
        atoms = self._atoms
        if n not in atoms or m not in atoms or n1 not in atoms or n2 not in atoms:
            raise AtomNotFound
        if not isinstance(mark, bool):
            raise TypeError('stereo mark should be bool')
        if (n, m) in self._cis_trans_stereo or (m, n) in self._cis_trans_stereo:
            raise IsChiral

        if (n, m) in self._chiral_cis_trans:
            self._cis_trans_stereo[(n, m)] = self._translate_cis_trans_sign(n, m, n1, n2, mark)
            if clean_cache:
                self.flush_cache()
        elif (m, n) in self._chiral_cis_trans:
            self._cis_trans_stereo[(m, n)] = self._translate_cis_trans_sign(m, n, n2, n1, mark)
            if clean_cache:
                self.flush_cache()
        else:
            raise NotChiral

    def flush_stereo_cache(self):
        """
        Flush chiral morgan and chiral centers cache.
        """
        self.__dict__.pop('_chiral_morgan', None)
        self.__dict__.pop('_MoleculeStereo__chiral_centers', None)

    def fix_stereo(self: 'MoleculeContainer'):
        """
        Reset stereo marks.
        """
        if self._atoms_stereo:  # filter tetrahedrons
            stereo_tetrahedrons = self._stereo_tetrahedrons
            atoms_stereo = {k: v for k, v in self._atoms_stereo.items() if k in stereo_tetrahedrons}
            self._atoms_stereo = self_atoms_stereo = {}
        else:
            atoms_stereo = {}

        if self._allenes_stereo:  # filter allenes
            stereo_allenes = self._stereo_allenes
            allenes_stereo = {k: v for k, v in self._allenes_stereo.items() if k in stereo_allenes}
            self._allenes_stereo = self_allenes_stereo = {}
        else:
            allenes_stereo = {}

        if self._cis_trans_stereo:  # filter cis-trans
            stereo_cis_trans = self._stereo_cis_trans
            cis_trans_stereo = {k: v for k, v in self._cis_trans_stereo.items() if k in stereo_cis_trans}
            self._cis_trans_stereo = self_stereo_cis_trans = {}
        else:
            cis_trans_stereo = {}

        old_stereo = len(atoms_stereo) + len(allenes_stereo) + len(cis_trans_stereo)
        while old_stereo:
            chiral_tetrahedrons = self._chiral_tetrahedrons
            chiral_allenes = self._chiral_allenes
            chiral_cis_trans = self._chiral_cis_trans

            tmp = {}
            for n, s in atoms_stereo.items():
                if n in chiral_tetrahedrons:
                    self_atoms_stereo[n] = s
                else:
                    tmp[n] = s
            atoms_stereo = tmp

            tmp = {}
            for n, s in allenes_stereo.items():
                if n in chiral_allenes:
                    self_allenes_stereo[n] = s
                else:
                    tmp[n] = s
            allenes_stereo = tmp

            tmp = {}
            for n, s in cis_trans_stereo.items():
                if n in chiral_cis_trans:
                    self_stereo_cis_trans[n] = s
                else:
                    tmp[n] = s
            cis_trans_stereo = tmp

            fail_stereo = len(atoms_stereo) + len(allenes_stereo) + len(cis_trans_stereo)
            if fail_stereo == old_stereo:
                break
            old_stereo = fail_stereo
            self.flush_stereo_cache()

    @cached_property
    def _wedge_map(self: 'MoleculeContainer'):
        plane = self._plane
        atoms_stereo = self._atoms_stereo
        allenes_centers = self._stereo_allenes_centers
        atoms = self._atoms
        used = set()
        wedge = []
        for n, s in self._allenes_stereo.items():
            env = self._stereo_allenes[n]
            term = self._stereo_allenes_terminals[n]
            order = [(*env[:2], *term), (*env[1::-1], *term[::-1])]
            if env[2]:
                order.append((env[2], env[1], *term))
                order.append((env[1], env[2], *term[::-1]))
            if env[3]:
                order.append((env[3], env[0], *term[::-1]))
                order.append((env[0], env[3], *term))
            order = sorted(order, key=lambda x: (x[0] in atoms_stereo, x[0] in allenes_centers,
                                                 -atoms[x[0]].atomic_number))
            while (order[0][0], order[0][2]) in used:
                order.append(order.pop(0))
            order = order[0]
            used.add((order[2], order[0]))
            s = self._translate_allene_sign(n, *order[:2])
            v = _allene_sign((*plane[order[0]], 1), plane[order[2]], plane[order[3]], (*plane[order[1]], 0))
            if not v:
                info(f'need 2d clean. wedge stereo ambiguous for atom {{{n}}}')
            if s:
                wedge.append((order[2], order[0], v))
            else:
                wedge.append((order[2], order[0], -v))

        for n, s in atoms_stereo.items():
            order = sorted(self._stereo_tetrahedrons[n], key=lambda x: (x in atoms_stereo, x in allenes_centers,
                                                                        -atoms[x].atomic_number, atoms[x].in_ring))
            while (order[0], n) in used:
                order.append(order.pop(0))
            used.add((n, order[0]))

            s = self._translate_tetrahedron_sign(n, order)
            # need recalculation if XY changed
            if len(order) == 3:
                v = _pyramid_sign((*plane[n], 0),
                                  (*plane[order[0]], 1), (*plane[order[1]], 0), (*plane[order[2]], 0))
            else:
                v = _pyramid_sign((*plane[order[3]], 0),
                                  (*plane[order[0]], 1), (*plane[order[1]], 0), (*plane[order[2]], 0))
            if not v:
                info(f'need 2d clean. wedge stereo ambiguous for atom {{{n}}}')
            if s:
                wedge.append((n, order[0], v))
            else:
                wedge.append((n, order[0], -v))
        return tuple(wedge)

    @property
    def _chiral_tetrahedrons(self) -> Set[int]:
        return self.__chiral_centers[0]

    @property
    def _chiral_cis_trans(self) -> Set[Tuple[int, int]]:
        return self.__chiral_centers[1]

    @property
    def _chiral_allenes(self) -> Set[int]:
        return self.__chiral_centers[2]

    @cached_property
    def _chiral_morgan(self: Union['MoleculeContainer', 'MoleculeStereo']) -> Dict[int, int]:
        if not self._atoms_stereo and not self._allenes_stereo and not self._cis_trans_stereo:
            return self.atoms_order
        morgan = self.atoms_order.copy()
        atoms_stereo = set(self._atoms_stereo)
        cis_trans_stereo = set(self._cis_trans_stereo)
        allenes_stereo = set(self._allenes_stereo)
        while True:
            # try iteratively differentiate stereo atoms.
            morgan, atoms_stereo, cis_trans_stereo, allenes_stereo, atoms_groups, cis_trans_groups, allenes_groups = \
                self.__differentiation(morgan, atoms_stereo, cis_trans_stereo, allenes_stereo)
            if not atoms_groups and not cis_trans_groups and not allenes_groups:
                break
            # for some rings differentiation by morgan impossible. try randomly set new weights.
            # sometimes this will lead to pseudo chiral centers and non-unique morgan.
            for group in atoms_groups:
                for n in group[:len(group) // 2]:  # set new weight in half of group randomly.
                    morgan[n] = -morgan[n]
            for group in cis_trans_groups:
                for n, _ in group[:len(group) // 2]:  # set new weight in half of group randomly.
                    morgan[n] = -morgan[n]
            for group in allenes_groups:
                for n in group[:len(group) // 2]:  # set new weight in half of group randomly.
                    morgan[n] = -morgan[n]
            morgan = _morgan(morgan, self.int_adjacency)
        return morgan

    @cached_property
    def _rings_tetrahedrons_linkers(self: 'MoleculeContainer') -> Dict[int, Tuple[int, int, int, int]]:
        """
        Ring-linkers tetrahedrons.

        Values are neighbors in first and second rings.
        """
        out = {}
        tetrahedrons = self._stereo_tetrahedrons
        for n, r in self.atoms_rings.items():
            if n in tetrahedrons:
                for nr, mr in combinations(r, 2):
                    if len(set(nr).intersection(mr)) == 1:
                        ni = nr.index(n)
                        mi = mr.index(n)
                        out[n] = (nr[ni - 1], nr[ni - len(nr) + 1], mr[mi - 1], mr[mi - len(mr) + 1])
                        break
        return out

    @cached_property
    def _rings_tetrahedrons(self: 'MoleculeContainer') -> Dict[int, Union[Tuple[int, int], Tuple[int], Tuple]]:
        """
        Tetrahedrons in rings, except ring-linkers.

        Values are out of ring atoms.
        """
        out = {}
        atoms_rings = self.atoms_rings
        tetrahedrons = self._stereo_tetrahedrons
        points = self._rings_tetrahedrons_linkers
        environment = self.not_special_connectivity
        for n, r in atoms_rings.items():
            if n in tetrahedrons and n not in points:
                out[n] = tuple(environment[n].difference(atoms_rings))
        return out

    @cached_property
    def _rings_cumulenes_linkers(self: 'MoleculeContainer') -> Dict[Tuple[int, int], Tuple[int, int, int, int]]:
        """
        Ring-linkers cumulenes except chords.

        Values are neighbors in first and second rings.
        """
        out = {}
        ar = self.atoms_rings
        chord = self._rings_cumulenes
        for (n, *_, m), (n1, m1, n2, m2) in self._stereo_cumulenes.items():
            if n in ar and m in ar and (n, m) not in chord:
                out[(n, m)] = (n1, n2, m1, m2)
        return out

    @cached_property
    def _rings_cumulenes(self: 'MoleculeContainer') -> Set[Tuple[int, int]]:
        """
        Cumulenes in rings always chiral.
        """
        out = set()
        ar = self.atoms_rings
        for n, *_, m in self._stereo_cumulenes:
            if n in ar and m in ar and not set(ar[n]).isdisjoint(ar[m]):
                out.add((n, m))
        return out

    @cached_property
    def _rings_cumulenes_attached(self: 'MoleculeContainer') -> Dict[Tuple[int, int],
                                                                     Union[Tuple[int, int], Tuple[int]]]:
        """
        Cumulenes attached to rings.

        Values are out of ring atoms.
        """
        ar = self.atoms_rings
        out = {}
        for (n, *_, m), (n1, m1, n2, m2) in self._stereo_cumulenes.items():
            if n in ar:
                if m in ar:
                    continue
                if m2:
                    out[(n, m)] = (m1, m2)
                else:
                    out[(n, m)] = (m1,)
            elif m in ar:
                if n2:
                    out[(n, m)] = (n1, n2)
                else:
                    out[(n, m)] = (n1,)
        return out

    @cached_property
    def __chiral_centers(self: Union['MoleculeStereo', 'MoleculeContainer']):
        atoms_rings = self.atoms_rings
        tetrahedrons = self._stereo_tetrahedrons
        cis_trans = self._stereo_cis_trans
        allenes_centers = self._stereo_allenes_centers
        cis_trans_terminals = self._stereo_cis_trans_terminals
        morgan = self._chiral_morgan

        # find new chiral atoms and bonds.
        # tetrahedron is chiral if all its neighbors are unique.
        chiral_t = {n for n, env in tetrahedrons.items() if len({morgan[x] for x in env}) == len(env)}
        # tetrahedrons-linkers is chiral if in each rings neighbors are unique.
        chiral_t.update(n for n, (n1, n2, m1, m2) in self._rings_tetrahedrons_linkers.items()
                        if morgan[n1] != morgan[n2] and morgan[m1] != morgan[m2])

        # required for axes detection.
        graph = {}
        stereogenic = set()
        pseudo = {}

        # double bond is chiral if neighbors of each terminal atom is unique.
        # ring-linkers and rings-attached also takes into account.
        chiral_c = set()
        chiral_a = set()
        for path, (n1, m1, n2, m2) in self._stereo_cumulenes.items():
            if morgan[n1] != morgan.get(n2, 0) and morgan[m1] != morgan.get(m2, 0):
                n, m = path[0], path[-1]
                if len(path) % 2:
                    chiral_a.add(path[len(path) // 2])
                else:
                    chiral_c.add((n, m))
                stereogenic.add(n)
                stereogenic.add(m)
        # ring cumulenes always chiral. can be already added.
        for nm in self._rings_cumulenes:
            n, m = nm
            if any(len(x) < 8 for x in atoms_rings[n]):  # skip small rings.
                if nm in chiral_c:  # remove already added small rings cumulenes.
                    chiral_c.discard(nm)
                elif (c := allenes_centers[n]) in chiral_a:
                    chiral_a.discard(c)
                continue
            elif nm in cis_trans:
                chiral_c.add(nm)
            else:
                chiral_a.add(allenes_centers[n])
            pseudo[m] = n
            graph[n] = set()
            stereogenic.add(n)

        # find chiral axes. build graph of stereogenic atoms in rings.
        # atoms connected then located in same ring or cumulene.
        for n, env in self._rings_tetrahedrons.items():
            if len(env) == 2:  # one or zero non-ring neighbors stereogenic.
                n1, n2 = env
                if morgan[n1] == morgan[n2]:  # only unique non-ring members required.
                    continue
            graph[n] = set()
            stereogenic.add(n)  # non-linker tetrahedrons in rings - stereogenic.
        for n, (n1, n2, m1, m2) in self._rings_tetrahedrons_linkers.items():
            graph[n] = set()
            if morgan[n1] != morgan[n2] or morgan[m1] != morgan[m2]:
                stereogenic.add(n)  # linkers with at least one unsymmetric ring.
        for n, m in self._rings_cumulenes_linkers:
            graph[n] = {m}
            graph[m] = {n}
            # stereogenic atoms already found.
        for (n, m), env in self._rings_cumulenes_attached.items():
            if len(env) == 2:
                n1, n2 = env
                if morgan[n1] == morgan[n2]:  # only unique non-ring members required.
                    continue
            if n in atoms_rings:
                graph[n] = set()  # non ring endpoints not required.
                stereogenic.add(n)  # mark as stereogenic
            else:
                graph[m] = set()
                stereogenic.add(m)

        if len(graph) > 1:  # add bonds to graph. bonds connects atoms in same rings and terminal atoms of cumulenes.
            for n, ms in graph.items():
                for r in atoms_rings[n]:
                    for m in r:
                        if n == m:
                            continue
                        elif m in graph:
                            ms.add(m)
                        elif m in pseudo and (m := pseudo[m]) != n:
                            ms.add(m)
            # remove not stereogenic terminals.
            while True:
                try:
                    n = next(n for n, ms in graph.items() if not ms or len(ms) == 1 and n not in stereogenic)
                except StopIteration:
                    break
                for m in graph.pop(n):
                    graph[m].discard(n)
            # update chiral atoms.
            for n in graph:
                if n in tetrahedrons:
                    chiral_t.add(n)
                elif n in allenes_centers:
                    chiral_a.add(allenes_centers[n])
                else:
                    chiral_c.add(cis_trans_terminals[n])

        # skip already marked.
        chiral_t.difference_update(self._atoms_stereo)
        chiral_a.difference_update(self._allenes_stereo)
        chiral_c.difference_update(self._cis_trans_stereo)
        return chiral_t, chiral_c, chiral_a

    def __differentiation(self: Union['MoleculeStereo', 'MoleculeContainer'], morgan,
                          atoms_stereo, cis_trans_stereo, allenes_stereo):
        bonds = self.int_adjacency

        tetrahedrons = self._stereo_tetrahedrons
        cis_trans = self._stereo_cis_trans
        allenes = self._stereo_allenes

        translate_tetrahedron = self._translate_tetrahedron_sign
        translate_cis_trans = self._translate_cis_trans_sign
        translate_allene = self._translate_allene_sign

        while True:
            morgan_update = {}
            atoms_groups = []
            cis_trans_groups = []
            allenes_groups = []
            # recalculate morgan weights with taking into account existing stereo marks.
            if atoms_stereo:
                grouped_stereo = defaultdict(list)
                for n in atoms_stereo:
                    grouped_stereo[morgan[n]].append(n)  # collect equal stereo atoms.
                for group in grouped_stereo.values():
                    if not len(group) % 2:  # only even number of equal stereo atoms give new stereo center.
                        # process only truly stereogenic.
                        if len(env := tetrahedrons[group[0]]) == len({morgan[x] for x in env}):
                            s = [n for n in group if translate_tetrahedron(n, sorted(tetrahedrons[n], key=morgan.get))]
                            if 0 < len(s) < len(group):  # RS pair required.
                                for m in s:
                                    morgan_update[m] = -morgan[m]
                            for n in group:  # prevent checks repeating.
                                atoms_stereo.discard(n)
                        else:  # stereo group in rings. unambiguous environment order impossible.
                            atoms_groups.append(group)

            if cis_trans_stereo:
                grouped_stereo = defaultdict(list)
                for nm in cis_trans_stereo:
                    n, m = nm
                    if (mn := morgan[n]) <= (mm := morgan[m]):
                        grouped_stereo[mn].append((n, nm))
                    else:
                        grouped_stereo[mm].append((m, nm))
                for group in grouped_stereo.values():
                    if not len(group) % 2:  # only even number of equal stereo bonds give new stereo center.
                        n1, m1, n2, m2 = cis_trans[group[0][1]]
                        if morgan[n1] != morgan.get(n2, 0) and morgan[m1] != morgan.get(m2, 0):
                            s = []
                            for x, nm in group:
                                n, m = nm
                                n1, m1, n2, m2 = cis_trans[nm]
                                if n2 is None:
                                    a = n1
                                else:
                                    a = min(n1, n2, key=morgan.get)
                                if m2 is None:
                                    b = m1
                                else:
                                    b = min(m1, m2, key=morgan.get)
                                if translate_cis_trans(n, m, a, b):
                                    s.append(x)
                            if 0 < len(s) < len(group):  # RS pair required.
                                for n in s:
                                    morgan_update[n] = -morgan[n]
                                for _, nm in group:
                                    cis_trans_stereo.discard(nm)
                        else:
                            cis_trans_groups.append(group)

            if allenes_stereo:
                grouped_stereo = defaultdict(list)
                for c in allenes_stereo:
                    grouped_stereo[morgan[c]].append(c)
                for group in grouped_stereo.values():
                    if not len(group) % 2:  # only even number of equal stereo bonds give new stereo center.
                        n1, m1, n2, m2 = allenes[group[0]]
                        if morgan[n1] != morgan.get(n2, 0) and morgan[m1] != morgan.get(m2, 0):
                            s = []
                            for c in group:
                                n1, m1, n2, m2 = allenes[c]
                                if n2 is None:
                                    a = n1
                                else:
                                    a = min(n1, n2, key=morgan.get)
                                if m2 is None:
                                    b = m1
                                else:
                                    b = min(m1, m2, key=morgan.get)
                                if translate_allene(c, a, b):
                                    s.append(c)
                            if 0 < len(s) < len(group):  # RS pair required.
                                for c in s:
                                    morgan_update[c] = -morgan[c]
                                for c in group:
                                    allenes_stereo.discard(c)
                        else:
                            allenes_groups.append(group)
            if not morgan_update:
                break
            morgan = _morgan({**morgan, **morgan_update}, bonds)
        return morgan, atoms_stereo, cis_trans_stereo, allenes_stereo, atoms_groups, cis_trans_groups, allenes_groups


__all__ = ['MoleculeStereo']
