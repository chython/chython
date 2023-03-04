# -*- coding: utf-8 -*-
#
#  Copyright 2018-2023 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from itertools import chain, product
from typing import Dict, List, Tuple, Union
from .bonds import Bond, QueryBond
from .graph import Graph
from ..algorithms.isomorphism import QueryIsomorphism
from ..algorithms.smiles import QuerySmiles
from ..algorithms.stereo import Stereo
from ..periodictable import Element, ListElement, QueryElement
from ..periodictable.element import Query


def _validate_neighbors(neighbors):
    if neighbors is None:
        neighbors = ()
    elif isinstance(neighbors, int):
        if neighbors < 0 or neighbors > 14:
            raise ValueError('neighbors should be in range [0, 14]')
        neighbors = (neighbors,)
    elif isinstance(neighbors, (tuple, list)):
        if not all(isinstance(n, int) for n in neighbors):
            raise TypeError('neighbors should be list or tuple of ints')
        if any(n < 0 or n > 14 for n in neighbors):
            raise ValueError('neighbors should be in range [0, 14]')
        if len(set(neighbors)) != len(neighbors):
            raise ValueError('neighbors should be unique')
        neighbors = tuple(sorted(neighbors))
    else:
        raise TypeError('neighbors should be int or list or tuple of ints')
    return neighbors


class QueryContainer(Stereo, Graph[Query, QueryBond], QueryIsomorphism, QuerySmiles):
    __slots__ = ('_neighbors', '_hybridizations', '_hydrogens', '_rings_sizes', '_heteroatoms', '_masked')

    _neighbors: Dict[int, Tuple[int, ...]]
    _hybridizations: Dict[int, Tuple[int, ...]]
    _hydrogens: Dict[int, Tuple[int, ...]]
    _rings_sizes: Dict[int, Tuple[int, ...]]
    _heteroatoms: Dict[int, Tuple[int, ...]]
    _masked: Dict[int, bool]

    def __init__(self):
        super().__init__()
        self._neighbors = {}
        self._hybridizations = {}
        self._hydrogens = {}
        self._rings_sizes = {}
        self._heteroatoms = {}
        self._masked = {}

    def add_atom(self, atom: Union[Query, Element, int, str], *args,
                 neighbors: Union[int, List[int], Tuple[int, ...], None] = None,
                 hybridization: Union[int, List[int], Tuple[int, ...], None] = None,
                 hydrogens: Union[int, List[int], Tuple[int, ...], None] = None,
                 rings_sizes: Union[int, List[int], Tuple[int, ...], None] = None,
                 heteroatoms: Union[int, List[int], Tuple[int, ...], None] = None,
                 masked: bool = False, **kwargs):
        if hybridization is None:
            hybridization = ()
        elif isinstance(hybridization, int):
            if hybridization < 1 or hybridization > 4:
                raise ValueError('hybridization should be in range [1, 4]')
            hybridization = (hybridization,)
        elif isinstance(hybridization, (tuple, list)):
            if not all(isinstance(h, int) for h in hybridization):
                raise TypeError('hybridizations should be list or tuple of ints')
            if any(h < 1 or h > 4 for h in hybridization):
                raise ValueError('hybridizations should be in range [1, 4]')
            if len(set(hybridization)) != len(hybridization):
                raise ValueError('hybridizations should be unique')
            hybridization = tuple(sorted(hybridization))
        else:
            raise TypeError('hybridization should be int or list or tuple of ints')

        if rings_sizes is None:
            rings_sizes = ()
        elif isinstance(rings_sizes, int):
            if rings_sizes < 3 and rings_sizes != 0:
                raise ValueError('rings should be greater or equal 3. ring equal to zero is no ring atom mark')
            rings_sizes = (rings_sizes,)
        elif isinstance(rings_sizes, (tuple, list)):
            if not all(isinstance(n, int) for n in rings_sizes):
                raise TypeError('rings should be list or tuple of ints')
            if any(n < 3 for n in rings_sizes):
                raise ValueError('rings should be greater or equal 3')
            if len(set(rings_sizes)) != len(rings_sizes):
                raise ValueError('rings should be unique')
            rings_sizes = tuple(sorted(rings_sizes))
        else:
            raise TypeError('rings should be int or list or tuple of ints')

        neighbors = _validate_neighbors(neighbors)
        hydrogens = _validate_neighbors(hydrogens)
        heteroatoms = _validate_neighbors(heteroatoms)

        if not isinstance(atom, Query):
            if isinstance(atom, Element):
                atom = QueryElement.from_atomic_number(atom.atomic_number)(atom.isotope)
            elif isinstance(atom, str):
                atom = QueryElement.from_symbol(atom)()
            elif isinstance(atom, int):
                atom = QueryElement.from_atomic_number(atom)()
            else:
                raise TypeError('QueryElement object expected')

        n = super().add_atom(atom, *args, **kwargs)
        self._neighbors[n] = neighbors
        self._hybridizations[n] = hybridization
        self._hydrogens[n] = hydrogens
        self._rings_sizes[n] = rings_sizes
        self._heteroatoms[n] = heteroatoms
        self._masked[n] = masked
        return n

    def add_bond(self, n, m, bond: Union[QueryBond, Bond, int, Tuple[int, ...]]):
        if isinstance(bond, Bond):
            bond = QueryBond.from_bond(bond)
        elif not isinstance(bond, QueryBond):
            bond = QueryBond(bond)

        sct = self._stereo_cis_trans_paths  # save
        sa = self._stereo_allenes_paths

        super().add_bond(n, m, bond)
        # remove stereo marks on bonded atoms and all its bonds
        if n in self._atoms_stereo:
            del self._atoms_stereo[n]
        if m in self._atoms_stereo:
            del self._atoms_stereo[m]
        if self._cis_trans_stereo:
            for nm, path in sct.items():
                if (n in path or m in path) and nm in self._cis_trans_stereo:
                    del self._cis_trans_stereo[nm]
        if self._allenes_stereo:
            for c, path in sa.items():
                if (n in path or m in path) and c in self._allenes_stereo:
                    del self._allenes_stereo[c]

    def copy(self) -> 'QueryContainer':
        copy = super().copy()

        copy._bonds = cb = {}
        for n, m_bond in self._bonds.items():
            cb[n] = cbn = {}
            for m, bond in m_bond.items():
                if m in cb:  # bond partially exists. need back-connection.
                    cbn[m] = cb[m][n]
                else:
                    cbn[m] = bond.copy()

        copy._neighbors = self._neighbors.copy()
        copy._hybridizations = self._hybridizations.copy()
        copy._hydrogens = self._hydrogens.copy()
        copy._heteroatoms = self._heteroatoms.copy()
        copy._rings_sizes = self._rings_sizes.copy()
        copy._atoms_stereo = self._atoms_stereo.copy()
        copy._allenes_stereo = self._allenes_stereo.copy()
        copy._cis_trans_stereo = self._cis_trans_stereo.copy()
        copy._masked = self._masked.copy()
        return copy

    def union(self, other: 'QueryContainer', *, remap: bool = False, copy: bool = True) -> 'QueryContainer':
        if not isinstance(other, QueryContainer):
            raise TypeError('QueryContainer expected')
        u, o = super().union(other, remap=remap, copy=copy)

        ub = u._bonds
        for n, m_bond in o._bonds.items():
            ub[n] = ubn = {}
            for m, bond in m_bond.items():
                if m in ub:  # bond partially exists. need back-connection.
                    ubn[m] = ub[m][n]
                else:
                    ubn[m] = bond.copy()

        u._neighbors.update(o._neighbors)
        u._hybridizations.update(o._hybridizations)
        u._hydrogens.update(o._hydrogens)
        u._rings_sizes.update(o._rings_sizes)
        u._heteroatoms.update(o._heteroatoms)
        u._masked.update(o._masked)
        return u

    def remap(self, mapping: Dict[int, int], *, copy=False) -> 'QueryContainer':
        atoms = self._atoms  # keep original atoms dict
        h = super().remap(mapping, copy=copy)

        mg = mapping.get
        hydrogens = self._hydrogens
        neighbors = self._neighbors
        hybridizations = self._hybridizations
        heteroatoms = self._heteroatoms
        rings_sizes = self._rings_sizes
        masked = self._masked

        if copy:
            hb = h._bonds
            hhg = h._hydrogens
            hn = h._neighbors
            hh = h._hybridizations
            hx = h._heteroatoms
            hrs = h._rings_sizes
            hm = h._masked

            # deep copy of bonds
            for n, m_bond in self._bonds.items():
                n = mg(n, n)
                hb[n] = hbn = {}
                for m, bond in m_bond.items():
                    m = mg(m, m)
                    if m in hb:  # bond partially exists. need back-connection.
                        hbn[m] = hb[m][n]
                    else:
                        hbn[m] = bond.copy()
        else:
            hb = {}
            hhg = {}
            hn = {}
            hh = {}
            hx = {}
            hrs = {}
            hm = {}

            for n, m_bond in self._bonds.items():
                n = mg(n, n)
                hb[n] = hbn = {}
                for m, bond in m_bond.items():
                    m = mg(m, m)
                    if m in hb:  # bond partially exists. need back-connection.
                        hbn[m] = hb[m][n]
                    else:
                        hbn[m] = bond

        for n in atoms:
            m = mg(n, n)
            hhg[m] = hydrogens[n]
            hn[m] = neighbors[n]
            hh[m] = hybridizations[n]
            hx[m] = heteroatoms[n]
            hrs[m] = rings_sizes[n]
            hm[m] = masked[n]

        if copy:
            return h  # noqa

        self._bonds = hb
        self._hydrogens = hhg
        self._neighbors = hn
        self._hybridizations = hh
        self._heteroatoms = hx
        self._rings_sizes = hrs
        self._masked = hm
        return self

    def enumerate_queries(self, *, enumerate_marks: bool = False):
        """
        Enumerate complex queries into multiple simple ones. For example `[N,O]-C` into `NC` and `OC`.

        :param enumerate_marks: enumerate multiple marks to separate queries
        """
        atoms = [(n, a._numbers) for n, a in self._atoms.items() if isinstance(a, ListElement)]
        bonds = [(n, m, b.order) for n, m, b in self.bonds() if len(b.order) > 1]
        for combo in product(*(x for *_, x in chain(atoms, bonds))):
            copy = self.copy()
            for (n, _), a in zip(atoms, combo):
                copy._atoms[n] = a = QueryElement.from_atomic_number(a)()
                a._attach_graph(copy, n)
            for (n, m, _), b in zip(bonds, combo[len(atoms):]):
                copy._bonds[n][m]._QueryBond__order = (b,)  # noqa

            if enumerate_marks:
                c = 0
                slices = []
                data = []
                for attr in ('_neighbors', '_hybridizations', '_hydrogens', '_heteroatoms', '_rings_sizes'):
                    tmp = [(n, v) for n, v in getattr(self, attr).items() if len(v) > 1]
                    if tmp:
                        data.extend(tmp)
                        slices.append((attr, c, c + len(tmp)))
                        c += len(tmp)

                for combo2 in product(*(x for _, x in data)):
                    copy2 = copy.copy()
                    for attr, i, j in slices:
                        attr = getattr(copy2, attr)
                        for (n, _), v in zip(data[i: j], combo2[i: j]):
                            attr[n] = (v,)
                    yield copy2
            else:
                yield copy

    def __getstate__(self):
        return {'atoms_stereo': self._atoms_stereo, 'allenes_stereo': self._allenes_stereo,
                'cis_trans_stereo': self._cis_trans_stereo, 'neighbors': self._neighbors,
                'hybridizations': self._hybridizations, 'hydrogens': self._hydrogens, 'masked': self._masked,
                'rings_sizes': self._rings_sizes, 'heteroatoms': self._heteroatoms, **super().__getstate__()}

    def __setstate__(self, state):
        super().__setstate__(state)
        self._atoms_stereo = state['atoms_stereo']
        self._allenes_stereo = state['allenes_stereo']
        self._cis_trans_stereo = state['cis_trans_stereo']
        self._neighbors = state['neighbors']
        self._hybridizations = state['hybridizations']
        self._hydrogens = state['hydrogens']
        self._rings_sizes = state['rings_sizes']
        self._heteroatoms = state['heteroatoms']
        self._masked = state['masked']


__all__ = ['QueryContainer']
