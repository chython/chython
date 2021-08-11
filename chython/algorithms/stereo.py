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
from CachedMethods import cached_property
from collections import defaultdict, deque
from logging import info
from typing import Dict, Optional, Set, Tuple, Union, TYPE_CHECKING
from ..exceptions import AtomNotFound, IsChiral, NotChiral


if TYPE_CHECKING:
    from chython import MoleculeContainer, QueryContainer
    Container = Union[MoleculeContainer, QueryContainer]


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


class Stereo:
    __slots__ = ()

    def clean_stereo(self: 'Container'):
        """
        Remove stereo data
        """
        self._atoms_stereo.clear()
        self._allenes_stereo.clear()
        self._cis_trans_stereo.clear()
        self.flush_cache()

    def get_mapping(self: 'Container', other: 'Container', **kwargs):
        atoms_stereo = self._atoms_stereo
        allenes_stereo = self._allenes_stereo
        cis_trans_stereo = self._cis_trans_stereo
        if atoms_stereo or allenes_stereo or cis_trans_stereo:
            other_atoms_stereo = other._atoms_stereo
            other_allenes_stereo = other._allenes_stereo
            other_cis_trans_stereo = other._cis_trans_stereo
            other_translate_tetrahedron_sign = other._translate_tetrahedron_sign
            other_translate_allene_sign = other._translate_allene_sign
            other_translate_cis_trans_sign = other._translate_cis_trans_sign

            tetrahedrons = self._stereo_tetrahedrons
            cis_trans = self._stereo_cis_trans
            allenes = self._stereo_allenes

            for mapping in super().get_mapping(other, **kwargs):
                for n, s in atoms_stereo.items():
                    m = mapping[n]
                    if m not in other_atoms_stereo:  # self stereo atom not stereo in other
                        break
                    # translate stereo mark in other in order of self tetrahedron
                    if other_translate_tetrahedron_sign(m, [mapping[x] for x in tetrahedrons[n]]) != s:
                        break
                else:
                    for n, s in allenes_stereo.items():
                        m = mapping[n]
                        if m not in other_allenes_stereo:  # self stereo allene not stereo in other
                            break
                        # translate stereo mark in other in order of self allene
                        nn, nm, *_ = allenes[n]
                        if other_translate_allene_sign(m, mapping[nn], mapping[nm]) != s:
                            break
                    else:
                        for nm, s in cis_trans_stereo.items():
                            n, m = nm
                            on, om = mapping[n], mapping[m]
                            if (on, om) not in other_cis_trans_stereo:
                                if (om, on) not in other_cis_trans_stereo:
                                    break  # self stereo cis_trans not stereo in other
                                else:
                                    nn, nm, *_ = cis_trans[nm]
                                    if other_translate_cis_trans_sign(om, on, mapping[nm], mapping[nn]) != s:
                                        break
                            else:
                                nn, nm, *_ = cis_trans[nm]
                                if other_translate_cis_trans_sign(on, om, mapping[nn], mapping[nm]) != s:
                                    break
                        else:
                            yield mapping
        else:
            yield from super().get_mapping(other, **kwargs)

    def _translate_tetrahedron_sign(self: 'Container', n, env):
        """
        Get sign of chiral tetrahedron atom for specified neighbors order

        :param n: stereo atom
        :param env: neighbors order
        """
        s = self._atoms_stereo[n]
        return self._translate_tetrahedron_sign_reversed(n, env, s)

    def _translate_tetrahedron_sign_reversed(self: 'Container', n, env, s):
        order = self._stereo_tetrahedrons[n]
        if len(order) == 3:
            if len(env) == 4:  # hydrogen atom passed to env
                atoms = self._atoms
                # hydrogen always last in order
                try:
                    order = (*order, next(x for x in env if atoms[x].atomic_number == 1))  # see translate scheme
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

    def _translate_cis_trans_sign(self: 'Container', n, m, nn, nm):
        """
        Get sign for specified opposite neighbors

        :param n: first double bonded atom
        :param m: last double bonded atom
        :param nn: neighbor of first atom
        :param nm: neighbor of last atom
        """
        try:
            s = self._cis_trans_stereo[(n, m)]
        except KeyError:
            s = self._cis_trans_stereo[(m, n)]
            n, m = m, n  # in alkenes sign not order depended
            nn, nm = nm, nn
        return self._translate_cis_trans_sign_reversed(n, m, nn, nm, s)

    def _translate_cis_trans_sign_reversed(self: 'Container', n, m, nn, nm, s):
        atoms = self._atoms
        n0, n1, n2, n3 = self._stereo_cis_trans[(n, m)]
        if nn == n0:  # same start
            t0 = 0
            if nm == n1:
                t1 = 1
            elif nm == n3 or n3 is None and atoms[nm].atomic_number == 1:
                t1 = 3
            else:
                raise KeyError
        elif nn == n1:
            t0 = 1
            if nm == n0:
                t1 = 0
            elif nm == n2 or n2 is None and atoms[nm].atomic_number == 1:
                t1 = 2
            else:
                raise KeyError
        elif nn == n2 or n2 is None and atoms[nn].atomic_number == 1:
            t0 = 2
            if nm == n1:
                t1 = 1
            elif nm == n3 or n3 is None and atoms[nm].atomic_number == 1:
                t1 = 3
            else:
                raise KeyError
        elif nn == n3 or n3 is None and atoms[nn].atomic_number == 1:
            t0 = 3
            if nm == n0:
                t1 = 0
            elif nm == n2 or n2 is None and atoms[nm].atomic_number == 1:
                t1 = 2
            else:
                raise KeyError
        else:
            raise KeyError

        if _alkene_translate[(t0, t1)]:
            return not s
        return s

    def _translate_allene_sign(self: 'Container', c, nn, nm):
        """
        get sign for specified opposite neighbors

        :param c: central double bonded atom
        :param nn: neighbor of first double bonded atom
        :param nm: neighbor of last double bonded atom
        """
        s = self._allenes_stereo[c]
        return self._translate_allene_sign_reversed(c, nn, nm, s)

    def _translate_allene_sign_reversed(self: 'Container', c, nn, nm, s):
        atoms = self._atoms

        n0, n1, n2, n3 = self._stereo_allenes[c]
        if nn == n0:  # same start
            t0 = 0
            if nm == n1:
                t1 = 1
            elif nm == n3 or n3 is None and atoms[nm].atomic_number == 1:
                t1 = 3
            else:
                raise KeyError
        elif nn == n1:
            t0 = 1
            if nm == n0:
                t1 = 0
            elif nm == n2 or n2 is None and atoms[nm].atomic_number == 1:
                t1 = 2
            else:
                raise KeyError
        elif nn == n2 or n2 is None and atoms[nn].atomic_number == 1:
            t0 = 2
            if nm == n1:
                t1 = 1
            elif nm == n3 or n3 is None and atoms[nm].atomic_number == 1:
                t1 = 3
            else:
                raise KeyError
        elif nn == n3 or n3 is None and atoms[nn].atomic_number == 1:
            t0 = 3
            if nm == n0:
                t1 = 0
            elif nm == n2 or n2 is None and atoms[nm].atomic_number == 1:
                t1 = 2
            else:
                raise KeyError
        else:
            raise KeyError

        if _alkene_translate[(t0, t1)]:
            return not s
        return s

    @cached_property
    def _stereo_cumulenes(self: 'Container') -> Dict[Tuple[int, ...], Tuple[int, int, Optional[int], Optional[int]]]:
        """
        Cumulenes which contains at least one non-hydrogen neighbor on both ends
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
            if any(b not in (1, 4) for m, b in nf.items() if m != n1) or \
                    any(b not in (1, 4) for m, b in nl.items() if m != m1):
                continue
            nn = [x for x in nf if x != n1 and atoms[x].atomic_number != 1]
            mn = [x for x in nl if x != m1 and atoms[x].atomic_number != 1]
            if nn and mn:
                sn = nn[1] if len(nn) == 2 else None
                sm = mn[1] if len(mn) == 2 else None
                cumulenes[path] = (nn[0], mn[0], sn, sm)
        return cumulenes

    @cached_property
    def _stereo_tetrahedrons(self: 'Container') -> Dict[int, Union[Tuple[int, int, int], Tuple[int, int, int, int]]]:
        """
        Tetrahedrons which contains at least 3 non-hydrogen neighbors
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
            env = tuple(x for x in bonds[n] if atoms[x].atomic_number != 1)
            if len(env) in (3, 4):
                tetrahedrons[n] = env
        return tetrahedrons

    @cached_property
    def _stereo_cis_trans(self) -> Dict[Tuple[int, int], Tuple[int, int, Optional[int], Optional[int]]]:
        """
        Cis-trans bonds which contains at least one non-hydrogen neighbor on both ends
        """
        return {(n, m): env for (n, *mid, m), env in self._stereo_cumulenes.items() if not len(mid) % 2}

    @cached_property
    def _stereo_cis_trans_paths(self) -> Dict[Tuple[int, int], Tuple[int, ...]]:
        return {(path[0], path[-1]): path for path in self._stereo_cumulenes if not len(path) % 2}

    @cached_property
    def _stereo_cis_trans_terminals(self) -> Dict[int, Tuple[int, int]]:
        """
        Cis-Trans terminal atoms to cis-trans key mapping
        """
        terminals = {}
        for nm in self._stereo_cis_trans_paths:
            n, m = nm
            terminals[n] = terminals[m] = nm
        return terminals

    @cached_property
    def _stereo_allenes(self) -> Dict[int, Tuple[int, int, Optional[int], Optional[int]]]:
        """
        Allenes which contains at least one non-hydrogen neighbor on both ends
        """
        return {path[len(path) // 2]: env for path, env in self._stereo_cumulenes.items() if len(path) % 2}

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
        return {c: (path[0], path[-1]) for c, path in self._stereo_allenes_paths.items()}

    @cached_property
    def _stereo_allenes_paths(self) -> Dict[int, Tuple[int, ...]]:
        return {path[len(path) // 2]: path for path in self._stereo_cumulenes if len(path) % 2}


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
                del self.__dict__['_MoleculeStereo__chiral_centers']
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
            self._atoms_stereo[n] = self._translate_tetrahedron_sign_reversed(n, env, mark)
            if clean_cache:
                self.flush_cache()
        elif n in self._chiral_allenes:
            self._allenes_stereo[n] = self._translate_allene_sign_reversed(n, *env, mark)
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
            self._cis_trans_stereo[(n, m)] = self._translate_cis_trans_sign_reversed(n, m, n1, n2, mark)
            if clean_cache:
                self.flush_cache()
        elif (m, n) in self._chiral_cis_trans:
            self._cis_trans_stereo[(m, n)] = self._translate_cis_trans_sign_reversed(m, n, n2, n1, mark)
            if clean_cache:
                self.flush_cache()
        else:
            raise NotChiral

    def _fix_stereo(self: 'MoleculeContainer'):
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
            # flush cache
            del self.__dict__['_MoleculeStereo__chiral_centers']

    @cached_property
    def _wedge_map(self: 'Container'):
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

    @property
    def _chiral_morgan(self: Union['MoleculeContainer', 'MoleculeStereo']) -> Dict[int, int]:
        if self._atoms_stereo or self._allenes_stereo or self._cis_trans_stereo:
            return self.__chiral_centers[3]
        return self.atoms_order

    @cached_property  # todo: speedup
    def _stereo_axises(self: 'MoleculeContainer') -> Tuple[Tuple[Tuple[int, ...], ...], Tuple[Tuple[int, ...], ...]]:
        """
        Get all stereogenic axises in rings with attached cumulenes. Stereogenic axises has only stereogenic atoms.
        """
        morgan = self.atoms_order
        bonds = self._bonds
        stereo_tetrahedrons = self._stereo_tetrahedrons
        cumulenes_terminals = {}
        for n, *_, m in self._stereo_cumulenes:
            cumulenes_terminals[n] = m
            cumulenes_terminals[m] = n

        out = []
        env = []
        axises = set()
        for c in self.connected_rings_cumulenes:
            out_c = []
            env_c = []
            adj = {n: {m: 1 for m, b in bonds[n].items() if m in c} for n in c}
            w_atoms = {n: morgan[n] for n in c}
            for mapping in self._get_automorphism_mapping(w_atoms, adj, w_atoms):
                sym = {k for k, v in mapping.items() if k != v}
                if not sym:
                    continue
                # get self-matched atoms connected with automorphic atoms
                ax = {k for k in mapping.keys() - sym if not sym.isdisjoint(bonds[k])}
                # add terminal stereogenic cumulenes
                c = 0
                tmp = set()
                for n in ax:
                    if n in stereo_tetrahedrons:
                        c += 1
                        tmp.add(n)
                    elif n in cumulenes_terminals:
                        tmp.add(n)
                        m = cumulenes_terminals[n]
                        if m not in ax:
                            tmp.add(m)
                            c += 1
                if c < 2:
                    continue
                ax = frozenset(tmp)
                if ax in axises:
                    continue
                axises.add(ax)
                if len(ax) > 2:  # get order of atoms
                    start = next(iter(ax))
                    # prepare distances
                    dist = {start: 0}
                    queue = deque([(start, 1)])
                    while queue:
                        n, d = queue.popleft()
                        d1 = d + 1
                        for m in adj[n].keys() - dist.keys():
                            queue.append((m, d1))
                            dist[m] = d
                    # get paths
                    seen = {start}
                    path = None
                    cut = {start}
                    stack = deque()
                    for n in adj[start]:
                        if mapping[n] not in cut:
                            stack.append((n, 0))
                            cut.add(n)
                    full_path = None
                    while stack:
                        n, d = stack.pop()
                        if not d:  # lets start
                            if path:
                                full_path = path[::-1]
                            path = [start]
                            d = 1
                        elif dist[n] != d:
                            continue
                        elif len(path) > d:
                            path = path[:d]
                        path.append(n)
                        if n in ax:
                            seen.add(n)
                            if seen == ax:
                                break
                        d += 1
                        for m in adj[n]:
                            if m in cut:
                                if m not in path:
                                    stack.append((m, d))
                            elif mapping[m] not in cut:
                                cut.add(m)
                                if m not in path:
                                    stack.append((m, d))

                    if full_path:
                        path = full_path + path[1:]
                    ax = tuple(n for n in path if n in ax)
                else:
                    ax = tuple(ax)
                out_c.append(ax)
                env_c.append(tuple(n for n in sym if n in bonds[ax[0]] or n in bonds[ax[-1]]))
            for ax, e in sorted(zip(out_c, env_c), key=lambda x: len(x[0]), reverse=True):
                out.append(ax)
                env.append(e)
        return tuple(out), tuple(env)

    @cached_property
    def __stereo_axises(self: 'MoleculeContainer'):
        bonds = self._bonds
        tetrahedrons = self._stereo_tetrahedrons
        cumulenes = self._stereo_cumulenes

        cumulenes_terminals = {}
        for path in cumulenes:
            cumulenes_terminals[path[0]] = cumulenes_terminals[path[-1]] = path
        axises = []
        for ax, env in zip(*self._stereo_axises):
            ax_t, ax_a, ax_c = set(), set(), set()
            checks = []
            axises.append((ax_t, ax_a, ax_c, checks))
            for n in ax:
                if n in tetrahedrons:
                    ax_t.add(n)
                else:
                    path = cumulenes_terminals[n]
                    if len(path) % 2:
                        ax_a.add(path)
                    else:
                        ax_c.add(path)
            for n in (ax[0], ax[-1]):
                if n in tetrahedrons:
                    ngb = tuple(m for m in tetrahedrons[n] if m not in env)
                else:
                    path = cumulenes_terminals[n]
                    ngb = tuple(m for m in bonds[n] if m not in path)
                if len(ngb) == 2:  # only these atoms should be checked
                    checks.append(ngb)
        return axises

    @cached_property
    def __chiral_centers(self: Union['MoleculeContainer', 'MoleculeStereo']):
        atoms_stereo = self._atoms_stereo
        cis_trans_stereo = self._cis_trans_stereo
        allenes_stereo = self._allenes_stereo

        bonds = {n: {m: int(b) for m, b in mb.items()} for n, mb in self._bonds.items()}
        morgan = self.atoms_order
        tetrahedrons = self._stereo_tetrahedrons.copy()
        cumulenes = self._stereo_cumulenes.copy()
        axises = self.__stereo_axises

        morgan_update = {}
        while True:
            # tetrahedron is chiral if all its neighbors are unique.
            chiral_t = {n for n, env in tetrahedrons.items() if len({morgan[x] for x in env}) == len(env)}
            # double bond is chiral if neighbors of each terminal atom is unique.
            chiral_c = set()
            chiral_a = set()
            for path, (n1, m1, n2, m2) in cumulenes.items():
                if morgan[n1] != morgan.get(n2, 0) and morgan[m1] != morgan.get(m2, 0):
                    if len(path) % 2:
                        chiral_a.add(path)
                    else:
                        chiral_c.add(path)
            # axis with 2 terminal chiral atoms is chiral
            for ax in axises:
                ax_t, ax_a, ax_c, check = ax
                if not chiral_t.isdisjoint(ax_t) or not ax_a.isdisjoint(chiral_a) or not ax_c.isdisjoint(chiral_c):
                    continue  # self chiral centers can't be in axises
                elif check and any(morgan[n] == morgan[m] for n, m in check):  # need additional check
                    continue  # not chiral
                chiral_t.update(ax_t)
                chiral_a.update(ax_a)
                chiral_c.update(ax_c)

            # separate equal constitutionally but unique by stereo type chiral centers
            # need for searching depended chiral centers
            if atoms_stereo:
                grouped_stereo = defaultdict(list)
                for n in chiral_t:
                    if n in atoms_stereo:
                        grouped_stereo[morgan[n]].append(n)  # collect equal stereo atoms
                for group in grouped_stereo.values():
                    if not len(group) % 2:  # only even number of equal stereo atoms give new stereo center
                        s = [n for n in group
                             if self._translate_tetrahedron_sign(n, sorted(tetrahedrons[n], key=morgan.get))]
                        if 0 < len(s) < len(group):  # RS pair required
                            for n in s:
                                morgan_update[n] = -morgan[n]
                    for n in group:  # remove seen stereo atoms
                        del tetrahedrons[n]
                        chiral_t.discard(n)

            if cis_trans_stereo:
                grouped_stereo = defaultdict(list)
                for path in chiral_c:
                    n, *_, m = path
                    mn, mm = morgan[n], morgan[m]
                    if (n, m) in cis_trans_stereo or (m, n) in cis_trans_stereo:
                        if mn <= mm:
                            grouped_stereo[mn].append((n, m, cumulenes[path], path))
                        else:
                            grouped_stereo[mm].append((m, n, cumulenes[path], path))
                for group in grouped_stereo.values():
                    if not len(group) % 2:  # only even number of equal stereo bonds give new stereo center
                        s = []
                        for n, m, (n1, m1, n2, m2), _ in group:
                            if n2 is None:
                                a = n1
                            else:
                                a = min(n1, n2, key=morgan.get)
                            if m2 is None:
                                b = m1
                            else:
                                b = min(m1, m2, key=morgan.get)
                            if self._translate_cis_trans_sign(n, m, a, b):
                                s.append(n)
                        if 0 < len(s) < len(group):  # RS pair required
                            for n in s:
                                morgan_update[n] = -morgan[n]
                    for *_, path in group:  # remove seen stereo atoms
                        del cumulenes[path]
                        chiral_c.discard(path)

            if allenes_stereo:
                grouped_stereo = defaultdict(list)
                for path in chiral_a:
                    c = path[len(path) // 2]
                    grouped_stereo[morgan[c]].append((c, cumulenes[path], path))
                for group in grouped_stereo.values():
                    if not len(group) % 2:  # only even number of equal stereo bonds give new stereo center
                        s = []
                        for c, (n1, m1, n2, m2), _ in group:
                            if n2 is None:
                                a = n1
                            else:
                                a = min(n1, n2, key=morgan.get)
                            if m2 is None:
                                b = m1
                            else:
                                b = min(m1, m2, key=morgan.get)
                            if self._translate_allene_sign(c, a, b):
                                s.append(c)
                        if 0 < len(s) < len(group):  # RS pair required
                            for c in s:
                                morgan_update[c] = -morgan[c]
                    for *_, path in group:  # remove seen stereo atoms
                        del cumulenes[path]
                        chiral_a.discard(path)

            tmp = []
            for ax in axises:
                ax_t, ax_a, ax_c, _ = ax
                # remove fully labeled axises
                if ax_t.issubset(tetrahedrons) and ax_a.issubset(cumulenes) and ax_c.issubset(cumulenes):
                    tmp.append(ax)
            axises = tmp

            if morgan_update:
                morgan = self._morgan({**morgan, **morgan_update}, bonds)
                morgan_update = {}
            else:
                return chiral_t, {(n, m) for n, *_, m in chiral_c}, {path[len(path) // 2] for path in chiral_a}, morgan


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


__all__ = ['MoleculeStereo', 'Stereo']
