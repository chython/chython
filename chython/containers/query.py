# -*- coding: utf-8 -*-
#
#  Copyright 2018-2021 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from CachedMethods import cached_args_method
from collections import defaultdict
from functools import cached_property
from itertools import chain, product
from typing import Dict, List, Set, Tuple, Union
from . import molecule  # cyclic imports resolve
from .bonds import Bond, QueryBond
from .graph import Graph
from ..algorithms.calculate2d import Calculate2DQuery
from ..algorithms.depict import DepictQuery
from ..algorithms.smiles import QuerySmiles
from ..algorithms.stereo import Stereo
from ..periodictable import AnyElement, AnyMetal, Element, ListElement, QueryElement
from ..periodictable.element import Query


class QueryContainer(Stereo, Graph[Query, QueryBond], QuerySmiles, DepictQuery, Calculate2DQuery):
    __slots__ = ('_neighbors', '_hybridizations', '_atoms_stereo', '_cis_trans_stereo', '_allenes_stereo',
                 '_hydrogens', '_rings_sizes', '_heteroatoms')

    _neighbors: Dict[int, Tuple[int, ...]]
    _hybridizations: Dict[int, Tuple[int, ...]]
    _atoms_stereo: Dict[int, bool]
    _allenes_stereo: Dict[int, bool]
    _cis_trans_stereo: Dict[Tuple[int, int], bool]
    _hydrogens: Dict[int, Tuple[int, ...]]
    _rings_sizes: Dict[int, Tuple[int, ...]]
    _heteroatoms: Dict[int, Tuple[int, ...]]

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

    def add_atom(self, atom: Union[Query, Element, int, str], *args,
                 neighbors: Union[int, List[int], Tuple[int, ...], None] = None,
                 hybridization: Union[int, List[int], Tuple[int, ...], None] = None,
                 hydrogens: Union[int, List[int], Tuple[int, ...], None] = None,
                 rings_sizes: Union[int, List[int], Tuple[int, ...], None] = None,
                 heteroatoms: Union[int, List[int], Tuple[int, ...], None] = None,
                 **kwargs):
        neighbors = self._validate_neighbors(neighbors)
        hybridization = self._validate_hybridization(hybridization)
        hydrogens = self._validate_neighbors(hydrogens)
        rings_sizes = self._validate_rings(rings_sizes)
        heteroatoms = self._validate_neighbors(heteroatoms)

        if not isinstance(atom, Query):
            if isinstance(atom, Element):
                atom = QueryElement.from_atomic_number(atom.atomic_number)(atom.isotope)
            elif isinstance(atom, str):
                atom = QueryElement.from_symbol(atom)()
            elif isinstance(atom, int):
                atom = QueryElement.from_atomic_number(atom)()
            else:
                raise TypeError('QueryElement object expected')

        _map = super().add_atom(atom, *args, **kwargs)
        self._neighbors[_map] = neighbors
        self._hybridizations[_map] = hybridization
        self._hydrogens[_map] = hydrogens
        self._rings_sizes[_map] = rings_sizes
        self._heteroatoms[_map] = heteroatoms
        return _map

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

    def delete_atom(self, n):
        bonds = set(self._bonds[n])  # save
        sct = self._stereo_cis_trans_paths
        sa = self._stereo_allenes_paths

        super().delete_atom(n)

        del self._neighbors[n]
        del self._hybridizations[n]
        del self._hydrogens[n]
        del self._rings_sizes[n]
        del self._heteroatoms[n]

        sas = self._atoms_stereo
        if n in sas:
            del sas[n]
        for m in bonds:
            if m in sas:
                del sas[m]
        if self._cis_trans_stereo:
            for nm, path in sct.items():
                if not bonds.isdisjoint(path) and nm in self._cis_trans_stereo:
                    del self._cis_trans_stereo[nm]
        if self._allenes_stereo:
            for c, path in sa.items():
                if not bonds.isdisjoint(path) and c in self._allenes_stereo:
                    del self._allenes_stereo[c]

    def delete_bond(self, n, m):
        sct = self._stereo_cis_trans_paths  # save
        sa = self._stereo_allenes_paths

        super().delete_bond(n, m)

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
        copy._neighbors = self._neighbors.copy()
        copy._hybridizations = self._hybridizations.copy()
        copy._hydrogens = self._hydrogens.copy()
        copy._heteroatoms = self._heteroatoms.copy()
        copy._rings_sizes = self._rings_sizes.copy()
        copy._atoms_stereo = self._atoms_stereo.copy()
        copy._allenes_stereo = self._allenes_stereo.copy()
        copy._cis_trans_stereo = self._cis_trans_stereo.copy()
        return copy

    def union(self, other: 'QueryContainer') -> 'QueryContainer':
        if not isinstance(other, QueryContainer):
            raise TypeError('QueryContainer expected')
        u = super().union(other)
        u._neighbors.update(other._neighbors)
        u._hybridizations.update(other._hybridizations)
        u._hydrogens.update(other._hydrogens)
        u._rings_sizes.update(other._rings_sizes)
        u._atoms_stereo.update(other._atoms_stereo)
        u._allenes_stereo.update(other._allenes_stereo)
        u._cis_trans_stereo.update(other._cis_trans_stereo)
        u._heteroatoms.update(other._heteroatoms)
        return u

    def get_mapping(self, other: Union['QueryContainer', 'molecule.MoleculeContainer'], **kwargs):
        if isinstance(other, (QueryContainer, molecule.MoleculeContainer)):
            return super().get_mapping(other, **kwargs)
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
                a._attach_to_graph(copy, n)
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
    def _screen_fingerprints(self) -> Tuple[Dict[int, Set[int]], ...]:
        """
        Fingerprints of all possible simple queries. Used for isomorphism tests filtering.
        """
        params = molecule.MoleculeContainer._fingerprint_config
        if not params or any(isinstance(a, (AnyElement, AnyMetal)) for _, a in self.atoms()):
            return ()  # skip queries with Any element

        fps = []
        for q in self.enumerate_queries():
            mol = molecule.MoleculeContainer()
            for n, a in q.atoms():
                mol.add_atom(Element.from_atomic_number(a.atomic_number)(a.isotope), n,
                             charge=a.charge, is_radical=a.is_radical)
            for n, m, b in q.bonds():
                mol.add_bond(n, m, Bond(b.order[0]))

            # extract only longest chains
            chains = []
            for c in sorted(mol._chains(**params), reverse=True, key=len):
                sc = set(c)
                if any(sc.issubset(x) for x in chains):
                    continue
                chains.append(c)
            if not chains:
                continue

            atoms = {idx: hash((atom.isotope or 0, atom.atomic_number, atom.charge, atom.is_radical))
                     for idx, atom in mol.atoms()}
            bonds = mol._bonds
            out = defaultdict(list)

            for frag in chains:
                var = [atoms[frag[0]]]
                for x, y in zip(frag, frag[1:]):
                    var.append(int(bonds[x][y]))
                    var.append(atoms[y])
                var = tuple(var)
                rev_var = var[::-1]
                out[hash(var if var > rev_var else rev_var)].extend(frag)

            fps.append({k: set(v) for k, v in out.items()})
        return tuple(fps)

    @cached_args_method
    def _component_fingerprints(self, component):
        """
        Fingerprints of specific component.
        """
        scope = set(self.connected_components[component])
        return tuple({k for k, v in fp.items() if not v.isdisjoint(scope)} for fp in self._screen_fingerprints)

    def _isomorphism_candidates(self, other, self_component, other_component):
        if self._screen_fingerprints and isinstance(other, molecule.MoleculeContainer):
            other_fingerprint = other._component_fingerprint(other_component)
            scope = set()
            for fp in self._component_fingerprints(self_component):
                if fp <= other_fingerprint.keys():
                    scope.update(x for k in fp for x in other_fingerprint[k])
            return scope
        return super()._isomorphism_candidates(other, self_component, other_component)

    @staticmethod
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

    @staticmethod
    def _validate_rings(rings):
        if rings is None:
            rings = ()
        elif isinstance(rings, int):
            if rings < 3 and rings != 0:
                raise ValueError('rings should be greater or equal 3. ring equal to zero is no ring atom mark')
            rings = (rings,)
        elif isinstance(rings, (tuple, list)):
            if not all(isinstance(n, int) for n in rings):
                raise TypeError('rings should be list or tuple of ints')
            if any(n < 3 for n in rings):
                raise ValueError('rings should be greater or equal 3')
            if len(set(rings)) != len(rings):
                raise ValueError('rings should be unique')
            rings = tuple(sorted(rings))
        else:
            raise TypeError('rings should be int or list or tuple of ints')
        return rings

    @staticmethod
    def _validate_hybridization(hybridization):
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
        return hybridization

    def __getstate__(self):
        return {'atoms_stereo': self._atoms_stereo, 'allenes_stereo': self._allenes_stereo,
                'cis_trans_stereo': self._cis_trans_stereo, 'neighbors': self._neighbors,
                'hybridizations': self._hybridizations, 'hydrogens': self._hydrogens,
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


__all__ = ['QueryContainer']
