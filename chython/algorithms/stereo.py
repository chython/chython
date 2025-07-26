# -*- coding: utf-8 -*-
#
#  Copyright 2019-2024 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from itertools import combinations, product
from logging import getLogger, INFO
from typing import Dict, Set, Tuple, Union, List, Optional, TYPE_CHECKING
from .morgan import _morgan
from ..exceptions import AtomNotFound, IsChiral, NotChiral


logger = getLogger('chython.stereo')
logger.setLevel(INFO)


if TYPE_CHECKING:
    from chython import MoleculeContainer


# atomic number constants
H = 1
C = 6

# 1  2
#  \ |
#   \|
#    n---3
#   /
#  /
# 0
_tetrahedron_translate = {(0, 1, 2): False, (1, 2, 0): False, (2, 0, 1): False,
                          (0, 2, 1): True, (1, 0, 2): True, (2, 1, 0): True,
                          (0, 3, 1): False, (3, 1, 0): False, (1, 0, 3): False,
                          (0, 1, 3): True, (1, 3, 0): True, (3, 0, 1): True,
                          (0, 2, 3): False, (2, 3, 0): False, (3, 0, 2): False,
                          (0, 3, 2): True, (3, 2, 0): True, (2, 0, 3): True,
                          (1, 3, 2): False, (3, 2, 1): False, (2, 1, 3): False,
                          (1, 2, 3): True, (2, 3, 1): True, (3, 1, 2): True}
# 2       1
#  \     /
#   n---m
#  /     \
# 0       3
_alkene_translate = {(0, 1): False, (1, 0): False, (0, 3): True, (3, 0): True,
                     (2, 3): False, (3, 2): False, (2, 1): True, (1, 2): True}


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


def _allene_sign(mark, u, v, w):
    # n    w
    # |   /
    # u--v
    ux, uy = u
    vx, vy = v
    wx, wy = w

    q2x = vx - ux
    q2y = vy - uy
    q3x = wx - vx
    q3y = wy - vy

    # cross vectors
    q2q3z = q2x * q3y - q2y * q3x

    dot = -mark * q2q3z
    if dot > 0:
        return 1
    elif dot < 0:
        return -1
    return 0


class MoleculeStereo:
    __slots__ = ()

    def clean_stereo(self: 'MoleculeContainer'):
        """
        Remove stereo data.
        """
        for _, a in self.atoms():
            a._stereo = None
        for *_, b in self.bonds():
            b._stereo = None
        self.flush_cache(keep_sssr=True, keep_components=True)

    @cached_property
    def tetrahedrons(self: 'MoleculeContainer') -> Tuple[int, ...]:
        """
        Carbon sp3 atom numbers.
        """
        tetra = []
        for n, atom in self.atoms():
            if atom == C and not atom.charge and not atom.is_radical:
                env = self._bonds[n]
                if all(b == 1 for b in env.values()):
                    if sum(int(b) for b in env.values()) > 4:
                        continue
                    tetra.append(n)
        return tuple(tetra)

    @cached_property
    def cumulenes(self: 'MoleculeContainer') -> List[Tuple[int, ...]]:
        """
        All double-bonds chains (e.g. alkenes, allenes, cumulenes).
        """
        atoms = self._atoms
        bonds = self._bonds

        adj = defaultdict(set)  # double bonds adjacency matrix
        for n, atom in atoms.items():
            if atom.is_forming_double_bonds:
                adj_n = adj[n].add
                for m, bond in bonds[n].items():
                    if bond == 2 and atoms[m].is_forming_double_bonds:
                        adj_n(m)
        if not adj:
            return []

        terminals = [x for x, y in adj.items() if len(y) == 1]  # list to keep atoms order!
        cumulenes = []
        while terminals:
            n = terminals.pop(0)
            m = adj[n].pop()
            path = [n, m]
            while m not in terminals:
                if len(bonds[m]) > 2:  # not cumulene. SO3, SO4- etc.
                    cumulenes.extend(zip(path, path[1:]))  # keep single double bonds instead of cumulene chain.
                    break
                adj_m = adj[m]
                adj_m.discard(n)
                n, m = m, adj_m.pop()
                path.append(m)
            else:
                terminals.remove(m)
                adj[m].pop()
                cumulenes.append(tuple(path))
        return cumulenes

    @cached_property
    def stereogenic_tetrahedrons(self: 'MoleculeContainer') -> Dict[int, Union[Tuple[int, int, int], Tuple[int, int, int, int]]]:
        """
        Tetrahedrons which contains at least 3 non-hydrogen neighbors and corresponding neighbors order.
        """
        #    2
        #    |
        # 1--K--3
        #    |
        #    4?
        atoms = self._atoms
        bonds = self._bonds
        tetrahedrons = {}
        for n in self.tetrahedrons:
            if any(not atoms[x].is_forming_single_bonds for x in bonds[n]):
                continue  # skip metal-carbon complexes
            env = tuple(x for x in bonds[n] if atoms[x] != H)
            if len(env) in (3, 4):
                tetrahedrons[n] = env
        return tetrahedrons

    @cached_property
    def stereogenic_cumulenes(self: 'MoleculeContainer') -> Dict[Tuple[int, ...], Tuple[int, int, Optional[int], Optional[int]]]:
        """
        Cumulenes which contains at least one non-hydrogen neighbor on both ends and corresponding neighbors order.
        """
        # 5       4
        #  \     /
        #   2---3
        #  /     \
        # 1       6
        bonds = self._bonds
        atoms = self._atoms
        cumulenes = {}
        for path in self.cumulenes:
            nf = bonds[path[0]]
            nl = bonds[path[-1]]
            n1, m1 = path[1], path[-2]
            if any(b == 3 or not atoms[m].is_forming_single_bonds and b != 8
                   for m, b in nf.items() if m != n1):
                continue  # skip X=C=C structures and metal-carbon complexes
            if any(b == 3 or not atoms[m].is_forming_single_bonds and b != 8
                   for m, b in nl.items() if m != m1):
                continue  # skip X=C=C structures and metal-carbon complexes
            nn = [x for x, b in nf.items() if x != n1 and atoms[x] != H and b != 8]
            mn = [x for x, b in nl.items() if x != m1 and atoms[x] != H and b != 8]
            if nn and mn:
                sn = nn[1] if len(nn) == 2 else None
                sm = mn[1] if len(mn) == 2 else None
                cumulenes[path] = (nn[0], mn[0], sn, sm)
        return cumulenes

    @cached_property
    def stereogenic_allenes(self) -> Dict[int, Tuple[int, int, Optional[int], Optional[int]]]:
        """
        Allenes which contains at least one non-hydrogen neighbor on both ends and corresponding neighbors order.
        """
        return {path[len(path) // 2]: env for path, env in self.stereogenic_cumulenes.items() if len(path) % 2}

    @cached_property
    def stereogenic_cis_trans(self) -> Dict[Tuple[int, int], Tuple[int, int, Optional[int], Optional[int]]]:
        """
        Cis-trans bonds which contains at least one non-hydrogen neighbor on both ends and corresponding neighbors order.
        """
        stereo = {}
        for path, env in self.stereogenic_cumulenes.items():
            if len(path) % 2:
                continue
            stereo[(path[0], path[-1])] = env
        return stereo

    @cached_property
    def ring_tetrahedrons(self: 'MoleculeContainer') -> Dict[int, Union[Tuple[int, int], Tuple[int], Tuple]]:
        """
        Tetrahedrons in rings, except ring-linkers. Values are non-ring atoms.
        """
        out = {}
        atoms_rings = self.atoms_rings
        tetrahedrons = self.stereogenic_tetrahedrons
        points = self.rings_linker_tetrahedrons
        environment = self.not_special_connectivity
        for n, r in atoms_rings.items():
            if n in tetrahedrons and n not in points:
                out[n] = tuple(environment[n].difference(atoms_rings))
        return out

    @cached_property
    def rings_linker_tetrahedrons(self: 'MoleculeContainer') -> Dict[int, Tuple[int, int, int, int]]:
        """
        A dictionary where the keys are tetrahedron atoms shared between two rings (not condensed rings) and the values
        are tuples representing their neighbors in the first and second rings respectively.
        """
        out = {}
        tetrahedrons = self.stereogenic_tetrahedrons
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
    def ring_cumulenes_terminals(self: 'MoleculeContainer') -> Set[Tuple[int, int]]:
        """
        Terminal atoms of inside ring cumulenes.
        """
        out = set()
        ar = self.atoms_rings
        for n, *_, m in self.stereogenic_cumulenes:
            if n in ar and m in ar and not set(ar[n]).isdisjoint(ar[m]):
                out.add((n, m))
        return out

    @cached_property
    def rings_linker_cumulenes_terminals(self: 'MoleculeContainer') -> Dict[Tuple[int, int], Tuple[int, int, int, int]]:
        """
        Terminal atoms of cumulenes connecting two rings. Values are neighbors in first and second rings.
        """
        out = {}
        ar = self.atoms_rings
        chord = self.ring_cumulenes_terminals
        for (n, *_, m), (n1, m1, n2, m2) in self.stereogenic_cumulenes.items():
            if n in ar and m in ar and (n, m) not in chord:
                out[(n, m)] = (n1, n2, m1, m2)
        return out

    @cached_property
    def ring_attached_cumulenes(self: 'MoleculeContainer') -> Dict[Tuple[int, int], Union[Tuple[int, int], Tuple[int]]]:
        """
        Cumulenes attached to rings from one side. Values are out of ring neighbor atoms.
        """
        ar = self.atoms_rings
        out = {}
        for (n, *_, m), (n1, m1, n2, m2) in self.stereogenic_cumulenes.items():
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

    @property
    def chiral_tetrahedrons(self) -> Set[int]:
        """
        Chiral tetrahedrons except already labeled ones.
        """
        return self.__chiral_centers[0]

    @property
    def chiral_cis_trans(self) -> Set[Tuple[int, int]]:
        """
        Chiral cis-trans bonds except already labeled ones.
        """
        return self.__chiral_centers[1]

    @property
    def chiral_allenes(self) -> Set[int]:
        """
        Chiral allenes except already labeled ones.
        """
        return self.__chiral_centers[2]

    def add_wedge(self: 'MoleculeContainer', n: int, m: int, mark: int, *, clean_cache=True):
        """
        Add stereo data by wedge notation of bonds. Use it for tetrahedrons of allenes.

        :param n: number of atom from which wedge bond started
        :param m: number of atom to which wedge bond coming
        :param mark: up bond is 1, down is -1
        """
        atoms = self._atoms
        if n not in atoms or m not in atoms or m not in self._bonds[n]:
            raise AtomNotFound
        elif atoms[n].stereo is not None:
            raise IsChiral
        elif c := self._stereo_allenes_centers.get(n):
            # allenes
            if atoms[c].stereo is not None:
                raise IsChiral
            elif c not in self.chiral_allenes:
                raise NotChiral

            t1, t2 = self._stereo_allenes_terminals[c]
            order = self.stereogenic_allenes[c]
            if atoms[m] == H:
                if t1 == n:
                    m1 = order[1]
                else:
                    t1, t2 = t2, t1
                    m1 = order[0]
                r = True
            else:
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
            if s := _allene_sign(mark, atoms[t1].xy, atoms[t2].xy, atoms[m1].xy):
                atoms[c]._stereo = s < 0 if r else s > 0
                if clean_cache:
                    self.flush_cache(keep_sssr=True, keep_components=True)
        # tetrahedrons
        elif n in self.chiral_tetrahedrons:
            th = self.stereogenic_tetrahedrons[n]
            am = atoms[m]
            if am == H:
                order = []
                for x in th:
                    ax = atoms[x]
                    order.append((ax.x, ax.y, 0))
                s = _pyramid_sign((am.x, am.y, mark), *order)
            else:
                order = []
                for x in th:
                    ax = atoms[x]
                    order.append((ax.x, ax.y, mark if x == m else 0))
                if len(order) == 3:
                    if len(self._bonds[n]) == 4:  # explicit hydrogen
                        x = next(x for x in self._bonds[n] if x not in th)
                        ax = atoms[x]
                        s = _pyramid_sign((ax.x, ax.y, 0), *order)
                    else:
                        an = atoms[n]
                        s = _pyramid_sign((an.x, an.y, 0), *order)
                else:
                    s = _pyramid_sign(order[-1], *order[:3])
            if s:
                atoms[n]._stereo = s > 0
                if clean_cache:
                    self.flush_cache(keep_components=True, keep_sssr=True)
        else:
            raise NotChiral

    def calculate_cis_trans_from_2d(self: 'MoleculeContainer', *, clean_cache=True):
        """
        Calculate cis-trans stereo bonds from given 2d coordinates. Unusable for SMILES and INCHI.
        """
        atoms = self._atoms
        flag = False
        while self.chiral_cis_trans:
            stereo = False
            for nm in self.chiral_cis_trans:
                n, m = nm
                n1, m1, *_ = self.stereogenic_cis_trans[nm]
                s = _cis_trans_sign(atoms[n1].xy, atoms[n].xy, atoms[m].xy, atoms[m1].xy)
                if s:
                    stereo = True
                    i, j = self._stereo_cis_trans_centers[n]
                    self._bonds[i][j]._stereo = s > 0
            if stereo:
                flag = True
                self.flush_stereo_cache()
            else:
                break
        if flag and clean_cache:
            self.flush_cache(keep_components=True, keep_sssr=True)

    def add_atom_stereo(self: 'MoleculeContainer', n: int, env: Tuple[int, ...], mark: bool, *, clean_cache=True):
        """
        Add stereo data for specified neighbors bypass. Use it for tetrahedrons or allenes.

        :param n: number of tetrahedron atom or central atom of allene.
        :param env: numbers of atoms with specified bypass
        :param mark: clockwise or anti bypass.

        See <https://www.daylight.com/dayhtml/doc/theory/theory.smiles.html> and <http://opensmiles.org/opensmiles.html>
        """
        try:
            atom = self._atoms[n]
        except KeyError:
            raise AtomNotFound
        if atom.stereo is not None:
            raise IsChiral
        if not isinstance(mark, bool):
            raise TypeError('stereo mark should be bool')

        if n in self.chiral_tetrahedrons:
            atom._stereo = self._translate_tetrahedron_sign(n, env, mark)
            if clean_cache:
                self.flush_cache(keep_components=True, keep_sssr=True)
        elif n in self.chiral_allenes:
            atom._stereo = self._translate_allene_sign(n, *env, mark)
            if clean_cache:
                self.flush_cache(keep_components=True, keep_sssr=True)
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

        if n not in self._stereo_cis_trans_counterpart or self._stereo_cis_trans_counterpart[n] != m:
            raise NotChiral
        i, j = self._stereo_cis_trans_centers[n]
        if self._bonds[i][j].stereo is not None:
            raise IsChiral

        if (n, m) in self.chiral_cis_trans:
            self._bonds[i][j]._stereo = self._translate_cis_trans_sign(n, m, n1, n2, mark)
            if clean_cache:
                self.flush_cache(keep_components=True, keep_sssr=True)
        elif (m, n) in self.chiral_cis_trans:
            self._bonds[i][j]._stereo = self._translate_cis_trans_sign(m, n, n2, n1, mark)
            if clean_cache:
                self.flush_cache(keep_components=True, keep_sssr=True)
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
        stereo_tetrahedrons = self.stereogenic_tetrahedrons
        stereo_allenes = self.stereogenic_allenes
        stereo_cis_trans = self._stereo_cis_trans_terminals
        atoms_stereo = []
        allenes_stereo = []
        cis_trans_stereo = []
        for n, a in self.atoms():
            if a.stereo is None:
                continue
            elif n in stereo_tetrahedrons:
                atoms_stereo.append((n, a, a.stereo))
            elif n in stereo_allenes:
                allenes_stereo.append((n, a, a.stereo))
            a._stereo = None  # flush stereo label

        for n, m, b in self.bonds():
            if b.stereo is None:
                continue
            elif (ta := stereo_cis_trans.get(n)) and ta == stereo_cis_trans.get(m):
                cis_trans_stereo.append((ta, b, b.stereo))
            b._stereo = None  # flush stereo label
        self.flush_stereo_cache()

        old_stereo = len(atoms_stereo) + len(allenes_stereo) + len(cis_trans_stereo)
        while old_stereo:
            chiral_tetrahedrons = self.chiral_tetrahedrons
            chiral_allenes = self.chiral_allenes
            chiral_cis_trans = self.chiral_cis_trans

            # filter out resolved
            tmp = []
            for n, a, s in atoms_stereo:
                if n in chiral_tetrahedrons:
                    a._stereo = s  # restore stereo
                else:
                    tmp.append((n, a, s))
            atoms_stereo = tmp

            tmp = []
            for n, a, s in allenes_stereo:
                if n in chiral_allenes:
                    a._stereo = s  # restore stereo
                else:
                    tmp.append((n, a, s))
            allenes_stereo = tmp

            tmp = []
            for ta, b, s in cis_trans_stereo:
                if ta in chiral_cis_trans:
                    b._stereo = s
                else:
                    tmp.append((ta, b, s))
            cis_trans_stereo = tmp

            fail_stereo = len(atoms_stereo) + len(allenes_stereo) + len(cis_trans_stereo)
            if fail_stereo == old_stereo:
                break
            old_stereo = fail_stereo
            self.flush_stereo_cache()

    @cached_property
    def _cis_trans_count(self) -> int:
        return sum(b.stereo is not None for *_, b in self.bonds())

    @cached_property
    def _stereo_cis_trans_centers(self) -> Dict[int, Tuple[int, int]]:
        """
        Cis-Trans terminal atoms to cis-trans key mapping. Key is central double bond in a cumulene chain.
        """
        terminals = {}
        for path in self.stereogenic_cumulenes:
            if len(path) % 2:
                continue
            n, m = path[0], path[-1]
            i = len(path) // 2
            terminals[n] = terminals[m] = (path[i - 1], path[i])
        return terminals

    @cached_property
    def _stereo_cis_trans_terminals(self) -> Dict[int, Tuple[int, int]]:
        """
        Cis-Trans terminal and central atoms to terminal pair mapping.
        """
        terminals = {}
        for path in self.stereogenic_cumulenes:
            if len(path) % 2:
                continue
            n, m = path[0], path[-1]
            i = len(path) // 2
            terminals[n] = terminals[m] = terminals[path[i]] = terminals[path[i - 1]] = (n, m)
        return terminals

    @cached_property
    def _stereo_cis_trans_counterpart(self) -> Dict[int, int]:
        """
        Cis-Trans terminal atoms counterparts
        """
        counterpart = {}
        for path in self.stereogenic_cumulenes:
            if len(path) % 2:
                continue
            n, m = path[0], path[-1]
            counterpart[n] = m
            counterpart[m] = n
        return counterpart

    @cached_property
    def _stereo_allenes_centers(self) -> Dict[int, int]:
        """
        Allene terminal atom to center mapping
        """
        terminals = {}
        for c, (n, m) in self._stereo_allenes_terminals.items():
            terminals[n] = terminals[m] = c
        return terminals

    @cached_property
    def _stereo_allenes_terminals(self) -> Dict[int, Tuple[int, int]]:
        """
        Allene center atom to terminals mapping
        """
        return {path[len(path) // 2]: (path[0], path[-1]) for path in self.stereogenic_cumulenes if len(path) % 2}

    def _translate_tetrahedron_sign(self: 'MoleculeContainer', n, env, s=None):
        """
        Get sign of chiral tetrahedron atom for specified neighbors order

        :param n: stereo atom
        :param env: neighbors order
        :param s: if None, use existing sign else translate given to molecule
        """
        if s is None:
            s = self._atoms[n].stereo
            if s is None:
                raise KeyError

        order = self.stereogenic_tetrahedrons[n]
        if len(order) == 3:
            if len(env) == 4:  # hydrogen atom passed to env
                # hydrogen always last in order
                try:
                    order = (*order, next(x for x in env if self._atoms[x] == H))  # see translate scheme
                except StopIteration:
                    raise KeyError
            elif len(env) != 3:  # pyramid or tetrahedron expected
                raise ValueError('invalid atoms list')
        elif len(env) not in (3, 4):  # pyramid or tetrahedron expected
            raise ValueError('invalid atoms list')

        translate = tuple(order.index(x) for x in env[:3])
        if _tetrahedron_translate[translate]:
            return not s
        return s

    def _translate_cis_trans_sign(self: 'MoleculeContainer', n, m, nn, nm, s=None):
        """
        Get sign for specified opposite neighbors

        :param n: first double bonded atom
        :param m: last double bonded atom
        :param nn: neighbor of first atom
        :param nm: neighbor of last atom
        :param s: if None, use existing sign else translate given to molecule
        """
        try:
            n0, n1, n2, n3 = self.stereogenic_cis_trans[(n, m)]
        except KeyError:
            n0, n1, n2, n3 = self.stereogenic_cis_trans[(m, n)]
            n, m = m, n  # in alkenes sign not order depended
            nn, nm = nm, nn

        if s is None:
            i, j = self._stereo_cis_trans_centers[n]
            s = self._bonds[i][j].stereo
            if s is None:
                raise KeyError

        if nn == n0:  # same start
            t0 = 0
            if nm == n1:
                t1 = 1
            elif nm == n3 or n3 is None and self._atoms[nm] == H:
                t1 = 3
            else:
                raise KeyError
        elif nn == n1:
            t0 = 1
            if nm == n0:
                t1 = 0
            elif nm == n2 or n2 is None and self._atoms[nm] == H:
                t1 = 2
            else:
                raise KeyError
        elif nn == n2 or n2 is None and self._atoms[nn] == H:
            t0 = 2
            if nm == n1:
                t1 = 1
            elif nm == n3 or n3 is None and self._atoms[nm] == H:
                t1 = 3
            else:
                raise KeyError
        elif nn == n3 or n3 is None and self._atoms[nn] == H:
            t0 = 3
            if nm == n0:
                t1 = 0
            elif nm == n2 or n2 is None and self._atoms[nm] == H:
                t1 = 2
            else:
                raise KeyError
        else:
            raise KeyError

        if _alkene_translate[(t0, t1)]:
            return not s
        return s

    def _translate_allene_sign(self: 'MoleculeContainer', c, nn, nm, s=None):
        """
        get sign for specified opposite neighbors

        :param c: central double bonded atom
        :param nn: neighbor of first double bonded atom
        :param nm: neighbor of last double bonded atom
        :param s: if None, use existing sign else translate given to molecule
        """
        if s is None:
            s = self._atoms[c].stereo
            if s is None:
                raise KeyError

        n0, n1, n2, n3 = self.stereogenic_allenes[c]
        if nn == n0:  # same start
            t0 = 0
            if nm == n1:
                t1 = 1
            elif nm == n3 or n3 is None and self._atoms[nm] == H:
                t1 = 3
            else:
                raise KeyError
        elif nn == n1:
            t0 = 1
            if nm == n0:
                t1 = 0
            elif nm == n2 or n2 is None and self._atoms[nm] == H:
                t1 = 2
            else:
                raise KeyError
        elif nn == n2 or n2 is None and self._atoms[nn] == H:
            t0 = 2
            if nm == n1:
                t1 = 1
            elif nm == n3 or n3 is None and self._atoms[nm] == H:
                t1 = 3
            else:
                raise KeyError
        elif nn == n3 or n3 is None and self._atoms[nn] == H:
            t0 = 3
            if nm == n0:
                t1 = 0
            elif nm == n2 or n2 is None and self._atoms[nm] == H:
                t1 = 2
            else:
                raise KeyError
        else:
            raise KeyError

        if _alkene_translate[(t0, t1)]:
            return not s
        return s

    @cached_property
    def _wedge_map(self: Union['MoleculeContainer', 'MoleculeStereo']):
        atoms = self._atoms

        overlap = set()
        space = []
        solved = []
        seen = set()
        for n, env in self.stereogenic_allenes.items():
            if atoms[n].stereo is None:
                continue
            term = self._stereo_allenes_terminals[n]
            overlap.update(term)  # don't allow incoming wedge to allenes terminals
            orders = [(*env[:2], *term, n, True), (*env[1::-1], *term[::-1], n, True)]
            if env[2]:
                orders.append((env[2], env[1], *term, n, True))
            if env[3]:
                orders.append((env[3], env[0], *term[::-1], n, True))
            space.append(orders)
        for n, env in self.stereogenic_tetrahedrons.items():
            if atoms[n].stereo is None:
                continue
            overlap.add(n)  # don't allow incoming wedge to stereo tetrahedrons
            order = list(env)
            orders = [(*order, n, False)]
            for _ in range(1, len(order)):
                order = order.copy()
                order.append(order.pop(0))
                orders.append((*order, n, False))
            space.append(orders)

        while space:
            ls = len(space)
            unsolved = []
            for orders in space:
                good = []
                if orders[0][-1]:
                    for x in orders:
                        x0 = x[0]
                        if x0 in seen or x0 not in overlap:
                            good.append(x)
                            seen.add(x[2])
                    if good:
                        solved.append(max(good, key=lambda x: (not atoms[x[0]].in_ring, atoms[x[0]].atomic_number)))
                    else:
                        unsolved.append(orders)
                else:
                    for x in orders:
                        x0 = x[0]
                        if x0 in seen or x0 not in overlap:
                            good.append(x)
                    if good:
                        seen.add(x[-2])
                        solved.append(max(good, key=lambda x: (not atoms[x[0]].in_ring, atoms[x[0]].atomic_number)))
                    else:
                        unsolved.append(orders)
            space = unsolved
            if len(unsolved) == ls:
                break

        solved = [y for x in solved if (y := self.__wedge_sign(x))]
        if not space:
            return solved

        for orders in product(*space):
            used = set()
            wedge = []
            for order in orders:
                if order[-1]:  # allene
                    if (order[0], order[2]) in used:
                        break
                    used.add((order[2], order[0]))
                    wedge.append(self.__wedge_sign(order))
                else:  # TH
                    n = order[-2]
                    if (order[0], n) in used:
                        break
                    used.add((n, order[0]))
                    wedge.append(self.__wedge_sign(order))
            else:  # found
                solved.extend(wedge)
                return solved
        logger.info('wedge stereo mapping failed')
        return solved

    def __wedge_sign(self: 'MoleculeContainer', order):
        if order[-1]:  # allene
            s = self._translate_allene_sign(order[-2], order[0], order[1])
            v = _allene_sign(1, self._atoms[order[2]].xy, self._atoms[order[3]].xy, self._atoms[order[1]].xy)
            if not v:
                logger.info(f'need 2d clean. allenes wedge stereo ambiguous for atom {order[-2]}')
            if s:
                return order[2], order[0], v
            else:
                return order[2], order[0], -v
        else:  # TH
            n = order[-2]
            s = self._translate_tetrahedron_sign(n, order[:-2])
            # need recalculation if XY changed
            ao0 = self._atoms[order[0]]
            ao1 = self._atoms[order[1]]
            ao2 = self._atoms[order[2]]
            if len(order) == 5:
                an = self._atoms[n]
                v = _pyramid_sign((an.x, an.y, 0),
                                  (ao0.x, ao0.y, 1),
                                  (ao1.x, ao1.y, 0),
                                  (ao2.x, ao2.y, 0))
            else:
                ao3 = self._atoms[order[3]]
                v = _pyramid_sign((ao3.x, ao3.y, 0),
                                  (ao0.x, ao0.y, 1),
                                  (ao1.x, ao1.y, 0),
                                  (ao2.x, ao2.y, 0))
            if not v:
                logger.info(f'need 2d clean. tetrahedron wedge stereo ambiguous for atom {n}')
            if s:
                return n, order[0], v
            else:
                return n, order[0], -v

    @cached_property
    def _chiral_morgan(self: Union['MoleculeContainer', 'MoleculeStereo']) -> Dict[int, int]:
        stereo_atoms = {n for n, a in self.atoms() if a.stereo is not None}
        stereo_bonds = {n for n, mb in self._bonds.items() if any(b.stereo is not None for m, b in mb.items())}
        if not stereo_atoms and not stereo_bonds:
            return self.atoms_order

        morgan = self.atoms_order.copy()
        atoms_stereo = stereo_atoms.intersection(self.tetrahedrons)
        allenes_stereo = stereo_atoms - atoms_stereo

        cis_trans_terminals = self._stereo_cis_trans_terminals
        cis_trans_stereo = {cis_trans_terminals[n] for n in stereo_bonds}

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
    def __chiral_centers(self: Union['MoleculeStereo', 'MoleculeContainer']):
        atoms_rings = self.atoms_rings
        tetrahedrons = self.stereogenic_tetrahedrons
        cis_trans = self.stereogenic_cis_trans
        allenes_centers = self._stereo_allenes_centers
        cis_trans_terminals = self._stereo_cis_trans_terminals
        cis_trans_centers = self._stereo_cis_trans_centers
        morgan = self._chiral_morgan

        # find new chiral atoms and bonds.
        # tetrahedron is chiral if all its neighbors are unique.
        chiral_t = {n for n, env in tetrahedrons.items() if len({morgan[x] for x in env}) == len(env)}
        # tetrahedrons-linkers is chiral if in each rings neighbors are unique.
        chiral_t.update(n for n, (n1, n2, m1, m2) in self.rings_linker_tetrahedrons.items()
                        if morgan[n1] != morgan[n2] and morgan[m1] != morgan[m2])

        # required for axes detection.
        graph = {}
        stereogenic = set()
        pseudo = {}

        # double bond is chiral if neighbors of each terminal atom is unique.
        # ring-linkers and rings-attached also takes into account.
        chiral_c = set()
        chiral_a = set()
        for path, (n1, m1, n2, m2) in self.stereogenic_cumulenes.items():
            if morgan[n1] != morgan.get(n2, 0) and morgan[m1] != morgan.get(m2, 0):
                n, m = path[0], path[-1]
                if len(path) % 2:
                    chiral_a.add(path[len(path) // 2])
                else:
                    chiral_c.add(n)
                stereogenic.add(n)
                stereogenic.add(m)
        # ring cumulenes always chiral. can be already added.
        for nm in self.ring_cumulenes_terminals:
            n, m = nm
            if any(len(x) < 8 for x in atoms_rings[n]):  # skip small rings.
                if n in chiral_c:  # remove already added small rings cumulenes.
                    chiral_c.discard(n)
                if m in chiral_c:
                    chiral_c.discard(m)
                elif n in allenes_centers and (c := allenes_centers[n]) in chiral_a:
                    chiral_a.discard(c)
                continue
            elif nm in cis_trans:
                chiral_c.add(n)
            else:
                chiral_a.add(allenes_centers[n])
            pseudo[m] = n
            graph[n] = set()
            stereogenic.add(n)

        # find chiral axes. build graph of stereogenic atoms in rings.
        # atoms connected then located in same ring or cumulene.
        for n, env in self.ring_tetrahedrons.items():
            if len(env) == 2:  # one or zero non-ring neighbors stereogenic.
                n1, n2 = env
                if morgan[n1] == morgan[n2]:  # only unique non-ring members required.
                    continue
            graph[n] = set()
            stereogenic.add(n)  # non-linker tetrahedrons in rings - stereogenic.
        for n, (n1, n2, m1, m2) in self.rings_linker_tetrahedrons.items():
            graph[n] = set()
            if morgan[n1] != morgan[n2] or morgan[m1] != morgan[m2]:
                stereogenic.add(n)  # linkers with at least one unsymmetric ring.
        for n, m in self.rings_linker_cumulenes_terminals:
            graph[n] = {m}
            graph[m] = {n}
            # stereogenic atoms already found.
        for (n, m), env in self.ring_attached_cumulenes.items():
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
                    chiral_c.add(n)

        # skip already marked.
        stereo_atoms = {n for n, a in self.atoms() if a.stereo is not None}
        chiral_t.difference_update(stereo_atoms)
        chiral_a.difference_update(stereo_atoms)
        diff = set()
        for n in chiral_c:
            i, j = cis_trans_centers[n]
            if self._bonds[i][j].stereo is None:
                diff.add(cis_trans_terminals[n])
        return chiral_t, diff, chiral_a

    def __differentiation(self: Union['MoleculeStereo', 'MoleculeContainer'], morgan,
                          atoms_stereo, cis_trans_stereo, allenes_stereo):
        bonds = self.int_adjacency

        tetrahedrons = self.stereogenic_tetrahedrons
        cis_trans = self.stereogenic_cis_trans
        allenes = self.stereogenic_allenes

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
