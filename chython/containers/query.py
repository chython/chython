# -*- coding: utf-8 -*-
#
#  Copyright 2018-2022 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from array import array
from functools import cached_property
from itertools import chain, product
from typing import Dict, List, Tuple, Union
from . import molecule  # cyclic imports resolve
from .bonds import Bond, QueryBond
from .graph import Graph
from ..algorithms.smiles import QuerySmiles
from ..algorithms.stereo import Stereo
from ..periodictable import AnyElement, AnyMetal, Element, ListElement, QueryElement
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


class QueryContainer(Stereo, Graph[Query, QueryBond], QuerySmiles):
    __slots__ = ('_neighbors', '_hybridizations', '_atoms_stereo', '_cis_trans_stereo', '_allenes_stereo',
                 '_hydrogens', '_rings_sizes', '_heteroatoms', '_masked')

    _neighbors: Dict[int, Tuple[int, ...]]
    _hybridizations: Dict[int, Tuple[int, ...]]
    _atoms_stereo: Dict[int, bool]
    _allenes_stereo: Dict[int, bool]
    _cis_trans_stereo: Dict[Tuple[int, int], bool]
    _hydrogens: Dict[int, Tuple[int, ...]]
    _rings_sizes: Dict[int, Tuple[int, ...]]
    _heteroatoms: Dict[int, Tuple[int, ...]]
    _masked: Dict[int, bool]

    def __init__(self):
        super().__init__()
        self._neighbors = {}
        self._hybridizations = {}
        self._atoms_stereo = {}
        self._allenes_stereo = {}
        self._cis_trans_stereo = {}
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

    def union(self, other: 'QueryContainer') -> 'QueryContainer':
        if not isinstance(other, QueryContainer):
            raise TypeError('QueryContainer expected')
        u = super().union(other)

        ub = u._bonds
        for n, m_bond in other._bonds.items():
            ub[n] = ubn = {}
            for m, bond in m_bond.items():
                if m in ub:  # bond partially exists. need back-connection.
                    ubn[m] = ub[m][n]
                else:
                    ubn[m] = bond.copy()

        u._neighbors.update(other._neighbors)
        u._hybridizations.update(other._hybridizations)
        u._hydrogens.update(other._hydrogens)
        u._rings_sizes.update(other._rings_sizes)
        u._atoms_stereo.update(other._atoms_stereo)
        u._allenes_stereo.update(other._allenes_stereo)
        u._cis_trans_stereo.update(other._cis_trans_stereo)
        u._heteroatoms.update(other._heteroatoms)
        u._masked.update(other._masked)
        return u

    def remap(self, mapping: Dict[int, int], *, copy=False) -> 'QueryContainer':
        """
        Change atom numbers.

        :param copy: keep original graph.
        """
        if len(mapping) != len(set(mapping.values())) or \
                not (self._atoms.keys() - mapping.keys()).isdisjoint(mapping.values()):
            raise ValueError('mapping overlap')

        mg = mapping.get
        charges = self._charges
        radicals = self._radicals
        hydrogens = self._hydrogens
        neighbors = self._neighbors
        hybridizations = self._hybridizations
        heteroatoms = self._heteroatoms
        rings_sizes = self._rings_sizes
        masked = self._masked

        if copy:
            h = self.__class__()
            ha = h._atoms
            hb = h._bonds
            hc = h._charges
            hr = h._radicals
            hhg = h._hydrogens
            hn = h._neighbors
            hh = h._hybridizations
            hx = h._heteroatoms
            hrs = h._rings_sizes
            hm = h._masked
            has = h._atoms_stereo
            hal = h._allenes_stereo
            hcs = h._cis_trans_stereo

            for n, atom in self._atoms.items():
                m = mg(n, n)
                atom = atom.copy()
                ha[m] = atom
                atom._attach_graph(h, m)

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
            ha = {}
            hb = {}
            hc = {}
            hr = {}
            hhg = {}
            hn = {}
            hh = {}
            hx = {}
            hrs = {}
            hm = {}
            has = {}
            hal = {}
            hcs = {}

            for n, atom in self._atoms.items():
                m = mg(n, n)
                ha[m] = atom
                atom._change_map(m)  # change mapping number

            for n, m_bond in self._bonds.items():
                n = mg(n, n)
                hb[n] = hbn = {}
                for m, bond in m_bond.items():
                    m = mg(m, m)
                    if m in hb:  # bond partially exists. need back-connection.
                        hbn[m] = hb[m][n]
                    else:
                        hbn[m] = bond

        for n in self._atoms:
            m = mg(n, n)
            hc[m] = charges[n]
            hr[m] = radicals[n]
            hhg[m] = hydrogens[n]
            hn[m] = neighbors[n]
            hh[m] = hybridizations[n]
            hx[m] = heteroatoms[n]
            hrs[m] = rings_sizes[n]
            hm[m] = masked[n]

        for n, stereo in self._atoms_stereo.items():
            has[mg(n, n)] = stereo
        for n, stereo in self._allenes_stereo.items():
            hal[mg(n, n)] = stereo
        for (n, m), stereo in self._cis_trans_stereo.items():
            hcs[(mg(n, n), mg(m, m))] = stereo

        if copy:
            return h

        self._atoms = ha
        self._bonds = hb
        self._charges = hc
        self._radicals = hr
        self._hydrogens = hhg
        self._neighbors = hn
        self._hybridizations = hh
        self._heteroatoms = hx
        self._rings_sizes = hrs
        self._masked = hm
        self._atoms_stereo = has
        self._allenes_stereo = hal
        self._cis_trans_stereo = hcs
        self.flush_cache()
        return self

    def get_mapping(self, other: Union['QueryContainer', 'molecule.MoleculeContainer'], /, *, _cython=True, **kwargs):
        # _cython - by default cython implementation enabled.
        # disable it by overriding method if Query Atoms or Containers logic changed.
        # Lv, Ts and Og in cython optimized mode treated as equal.
        if isinstance(other, QueryContainer):
            return super().get_mapping(other, **kwargs)
        elif isinstance(other, molecule.MoleculeContainer):
            return super().get_mapping(other, _cython=_cython, **kwargs)
        raise TypeError('MoleculeContainer or QueryContainer expected')

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
                copy._bonds[n][m]._QueryBond__order = (b,)

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

    @cached_property
    def _cython_compiled_query(self):
        # long I:
        # bond: single, double, triple, aromatic, special = 5 bit
        # bond in ring: 2 bit
        # atom: H-Ba: 56 bit
        # transfer bit

        # long II:
        # atom La-Mc: 59 bit
        # Lv-Ts-Og: 3 elements packed into 1 bit.
        # hybridizations: 1-4 = 4 bit

        # long III:
        # isotope: not specified, isotope - common_isotope = -8 - +8 = 18 bit
        # is_radical: 2 bit
        # charge: -4 - +4: 9 bit
        # implicit_hydrogens: 0-4 = 5 bit
        # neighbors: 0-14 = 15 bit
        # heteroatoms: 0-14 = 15 bit

        # long IV:
        # ring_sizes: not-in-ring bit, 3-atom ring, 4-...., 65-atom ring

        # int V: bonds closures
        # padding: 1 bit
        # bond: single, double, triple, aromatic, special = 5 bit
        # bond in ring: 2 bit
        from ..files._mdl.mol import common_isotopes

        _components, _closures = self._compiled_query
        components = []
        for c in _components:
            mapping = {n: i for i, (n, *_) in enumerate(c)}
            masks1 = []
            masks2 = []
            masks3 = []
            masks4 = []
            for *_, a, b in c:
                if isinstance(a, AnyMetal):  # isotope, radical, charge, hydrogens and heteroatoms states ignored
                    # elements except 1, 2, 5, 6, 7, 8, 9, 10, 14, 15, 16, 17, 18, 32, 33, 34, 35, 36, 51, 52, 53, 54
                    v1 = 0x0060707ffc1fff87
                    v2 = 0xfffffffffffffff0
                    v3 = 0xffffffffc0007fff
                    v4 = 0xffffffffffffffff
                else:
                    if isinstance(a, AnyElement):
                        v1 = 0x01ffffffffffffff
                        v2 = 0xfffffffffffffff0
                    else:
                        if isinstance(a, ListElement):
                            v1 = v2 = 0
                            for n in a._numbers:
                                if n > 56:
                                    if n > 116:  # Ts, Og
                                        n = 116
                                    v1 |= 1  # set transfer bit
                                    v2 |= 1 << (120 - n)
                                else:
                                    v1 |= 1 << (57 - n)
                        elif (n := a.atomic_number) > 56:
                            if n > 116:  # Ts, Og
                                n = 116
                            v1 = 1  # transfer bit
                            v2 = 1 << (120 - n)
                        else:
                            v1 = 1 << (57 - n)
                            v2 = 0
                    if a.isotope:
                        v3 = 1 << (a.isotope - common_isotopes[a.atomic_symbol] + 54)
                        if a.is_radical:
                            v3 |= 0x200000000000
                        else:
                            v3 |= 0x100000000000
                    elif a.is_radical:  # any isotope
                        v3 = 0xffffe00000000000
                    else:
                        v3 = 0xffffd00000000000

                    v3 |= 1 << (a.charge + 39)

                    if not a.implicit_hydrogens:
                        v3 |= 0x7c0000000
                    else:
                        for h in a.implicit_hydrogens:
                            v3 |= 1 << (h + 30)

                    if not a.heteroatoms:
                        v3 |= 0x7fff
                    else:
                        for n in a.heteroatoms:
                            v3 |= 1 << n

                    if a.ring_sizes:
                        if a.ring_sizes[0]:
                            v4 = 0
                            for r in a.ring_sizes:
                                if r > 65:  # big rings not supported
                                    continue
                                v4 |= 1 << (65 - r)
                            if not v4:  # only 65+ rings. set as rings-free.
                                v4 = 0x8000000000000000
                        else:  # not in rings
                            v4 = 0x8000000000000000
                    else:  # any rings
                        v4 = 0xffffffffffffffff

                if not a.neighbors:
                    v3 |= 0x3fff8000
                else:
                    for n in a.neighbors:
                        v3 |= 1 << (n + 15)

                if not a.hybridization:
                    v2 |= 0xf
                else:
                    for n in a.hybridization:
                        v2 |= 1 << (n - 1)

                if b is not None:
                    for o in b.order:
                        if o == 1:
                            v1 |= 0x0800000000000000
                        elif o == 4:
                            v1 |= 0x4000000000000000
                        elif o == 2:
                            v1 |= 0x1000000000000000
                        elif o == 3:
                            v1 |= 0x2000000000000000
                        else:
                            v1 |= 0x8000000000000000
                    if b.in_ring is None:
                        v1 |= 0x0600000000000000
                    elif b.in_ring:
                        v1 |= 0x0400000000000000
                    else:
                        v1 |= 0x0200000000000000

                masks1.append(v1)
                masks2.append(v2)
                masks3.append(v3)
                masks4.append(v4)

            closures = [0] * len(c)  # closures amount
            q_from = [0] * len(c)
            q_to = [0] * len(c)
            indices = [0] * sum(len(ms) for n, ms in _closures.items() if n in mapping)
            bonds = indices.copy()

            start = 0
            for n, ms in _closures.items():
                if (i := mapping.get(n)) is not None:
                    closures[i] = len(ms)
                    q_from[i] = start
                    for j, (m, b) in enumerate(ms, start):
                        v = 0x01ffffffffffffff  # atom doesn't matter.
                        for o in b.order:
                            if o == 1:
                                v |= 0x0800000000000000
                            elif o == 4:
                                v |= 0x4000000000000000
                            elif o == 2:
                                v |= 0x1000000000000000
                            elif o == 3:
                                v |= 0x2000000000000000
                            else:
                                v |= 0x8000000000000000
                        if b.in_ring is None:
                            v |= 0x0600000000000000
                        elif b.in_ring:
                            v |= 0x0400000000000000
                        else:
                            v |= 0x0200000000000000
                        bonds[j] = v
                        indices[j] = mapping[m]
                    start += len(ms)
                    q_to[i] = start
            components.append((array('L', [n for n, *_ in c]), array('I', [0] + [mapping[x] for _, x, *_ in c[1:]]),
                               array('Q', masks1), array('Q', masks2), array('Q', masks3), array('Q', masks4),
                               array('I', closures), array('I', q_from), array('I', q_to),
                               array('I', indices), array('Q', bonds)))
        return components

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
