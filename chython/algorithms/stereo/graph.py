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
from typing import Dict, Optional, Tuple, TYPE_CHECKING, Union


if TYPE_CHECKING:
    from chython import MoleculeContainer, QueryContainer
    Container = Union[MoleculeContainer, QueryContainer]


_heteroatoms = {5, 6, 7, 8, 14, 15, 16, 17, 33, 34, 35, 52, 53}

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

# allowed atoms. these atoms have stable covalent bonds.
_organic_subset = {1, 5, 6, 7, 8, 9, 14, 15, 16, 17, 33, 34, 35, 52, 53}


class Stereo:
    __slots__ = ()

    @cached_property
    def cumulenes(self) -> Tuple[Tuple[int, ...], ...]:
        """
        Alkenes, allenes and cumulenes atoms numbers.
        """
        return tuple(self._cumulenes())

    @cached_property
    def tetrahedrons(self: 'Container') -> Tuple[int, ...]:
        """
        Carbon sp3 atoms numbers.
        """
        atoms = self._atoms
        bonds = self._bonds
        charges = self._charges
        radicals = self._radicals

        tetra = []
        for n, atom in atoms.items():
            if atom.atomic_number == 6 and not charges[n] and not radicals[n]:
                env = bonds[n]
                if all(x == 1 for x in env.values()):
                    if sum(int(x) for x in env.values()) > 4:
                        continue
                    tetra.append(n)
        return tuple(tetra)

    def clean_stereo(self: 'Container'):
        """
        Remove stereo data.
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

    def _translate_tetrahedron_sign(self: 'Container', n, env, s=None):
        """
        Get sign of chiral tetrahedron atom for specified neighbors order

        :param n: stereo atom
        :param env: neighbors order
        :param s: if None, use existing sign else translate given to molecule
        """
        if s is None:
            s = self._atoms_stereo[n]

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

    def _translate_cis_trans_sign(self: 'Container', n, m, nn, nm, s=None):
        """
        Get sign for specified opposite neighbors

        :param n: first double bonded atom
        :param m: last double bonded atom
        :param nn: neighbor of first atom
        :param nm: neighbor of last atom
        :param s: if None, use existing sign else translate given to molecule
        """
        if s is None:
            try:
                s = self._cis_trans_stereo[(n, m)]
            except KeyError:
                s = self._cis_trans_stereo[(m, n)]
                n, m = m, n  # in alkenes sign not order depended
                nn, nm = nm, nn

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

    def _translate_allene_sign(self: 'Container', c, nn, nm, s=None):
        """
        get sign for specified opposite neighbors

        :param c: central double bonded atom
        :param nn: neighbor of first double bonded atom
        :param nm: neighbor of last double bonded atom
        :param s: if None, use existing sign else translate given to molecule
        """
        if s is None:
            s = self._allenes_stereo[c]

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

    def _cumulenes(self: 'Container', heteroatoms=False):
        atoms = self._atoms
        bonds = self._bonds

        adj = defaultdict(set)  # double bonds adjacency matrix
        if heteroatoms:
            for n, atom in atoms.items():
                if atom.atomic_number in _heteroatoms:
                    adj_n = adj[n].add
                    for m, bond in bonds[n].items():
                        if bond == 2 and atoms[m].atomic_number in _heteroatoms:
                            adj_n(m)
        else:
            for n, atom in atoms.items():
                if atom.atomic_number == 6:
                    adj_n = adj[n].add
                    for m, bond in bonds[n].items():
                        if bond == 2 and atoms[m].atomic_number == 6:
                            adj_n(m)
        if not adj:
            return ()

        terminals = [x for x, y in adj.items() if len(y) == 1]
        cumulenes = []
        while terminals:
            n = terminals.pop(0)
            m = adj[n].pop()
            path = [n, m]
            while m not in terminals:
                adj_m = adj[m]
                if len(adj_m) > 2:  # not cumulene. SO3 etc.
                    cumulenes.extend(zip(path, path[1:]))  # keep single double bonds.
                    break
                adj_m.discard(n)
                n, m = m, adj_m.pop()
                path.append(m)
            else:
                terminals.remove(m)
                adj[m].pop()
                cumulenes.append(tuple(path))
        return cumulenes

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
            if any(b not in (1, 4) and atoms[m].atomic_number not in _organic_subset for m, b in nf.items() if m != n1)\
                    or any(b not in (1, 4) and atoms[m].atomic_number not in _organic_subset
                           for m, b in nl.items() if m != m1):
                continue  # skip X=C=C structures and metal-carbon complexes
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
            if any(atoms[x].atomic_number not in _organic_subset for x in bonds[n]):
                continue  # skip metal-carbon complexes
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


__all__ = ['Stereo']
