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
from abc import ABC, abstractmethod
from CachedMethods import cached_args_method
from functools import cached_property
from typing import Dict, Optional, Tuple, Iterable, Iterator, Union, List, Type
from .bonds import Bond, DynamicBond, QueryBond
from ..algorithms.components import GraphComponents
from ..algorithms.isomorphism import Isomorphism
from ..algorithms.morgan import Morgan
from ..algorithms.sssr import SSSR
from ..exceptions import AtomNotFound
from ..periodictable import AnyAtom


class Graph(GraphComponents, Morgan, SSSR, Isomorphism, ABC):
    __slots__ = ('_atoms', '_bonds', '_plane', '_charges', '_radicals', '__meta', '__name', '__dict__', '__weakref__')
    __class_cache__ = {}

    def __init__(self):
        self._charges: Dict[int, int] = {}
        self._radicals: Dict[int, bool] = {}
        self._plane: Dict[int, Tuple[float, float]] = {}
        self.__meta = {}
        self.__name = ''

    def __getstate__(self):
        state = {'atoms': self._atoms, 'bonds': self._bonds, 'meta': self.__meta, 'plane': self._plane,
                 'charges': self._charges, 'radicals': self._radicals, 'name': self.__name}
        if self.__class_cache__.get('save_cache', False):
            state['cache'] = {k: v for k, v in self.__dict__.items() if k != '__cached_method___hash__'}
        return state

    def __setstate__(self, state):
        self._atoms = state['atoms']
        for n, a in state['atoms'].items():
            a._attach_to_graph(self, n)
        self._charges = state['charges']
        self._radicals = state['radicals']
        self._plane = state['plane']
        self._bonds = state['bonds']
        self.__meta = state['meta']
        self.__name = state['name']
        if 'cache' in state:  # >= 4.1.15
            self.__dict__.update(state['cache'])

    @classmethod
    def pickle_save_cache(cls, arg: bool):
        """
        Store cache of Graph into pickle for speedup loading
        """
        cls.__class_cache__['save_cache'] = arg

    def __len__(self):
        return len(self._atoms)

    def __iter__(self):
        return iter(self._atoms)

    def __contains__(self, n: int):
        return n in self._atoms

    def __bool__(self):
        return bool(self._atoms)

    def atom(self, n: int) -> AnyAtom:
        return self._atoms[n]

    def has_atom(self, n: int) -> bool:
        return n in self._atoms

    def atoms(self) -> Iterator[Tuple[int, AnyAtom]]:
        """
        iterate over all atoms
        """
        return iter(self._atoms.items())

    @cached_property
    def atoms_count(self) -> int:
        return len(self._atoms)

    @cached_property
    def atoms_numbers(self) -> Tuple[int, ...]:
        return tuple(self._atoms)

    @cached_args_method
    def environment(self, atom: int, include_bond: bool = True, include_atom: bool = True) -> \
            Tuple[Union[Tuple[int, Union[Bond, DynamicBond], AnyAtom],
                        Tuple[int, AnyAtom],
                        Tuple[int, Union[Bond, DynamicBond]],
                        int], ...]:
        """
        groups of (atom_number, bond, atom) connected to atom or
        groups of (atom_number, bond) connected to atom or
        groups of (atom_number, atom) connected to atom or
        neighbors atoms connected to atom

        :param atom: number
        :param include_atom: include atom object
        :param include_bond: include bond object
        """
        if include_atom:
            atoms = self._atoms
            if include_bond:
                return tuple((n, bond, atoms[n]) for n, bond in self._bonds[atom].items())
            return tuple((n, atoms[n]) for n in self._bonds[atom])
        elif include_bond:
            return tuple(self._bonds[atom].items())
        return tuple(self._bonds[atom])

    def bond(self, n: int, m: int) -> Union[Bond, DynamicBond]:
        return self._bonds[n][m]

    def has_bond(self, n: int, m: int) -> bool:
        try:
            self._bonds[n]  # check if atom exists
            return n in self._bonds[m]
        except KeyError:
            raise AtomNotFound

    def bonds(self) -> Iterator[Tuple[int, int, Union[Bond, QueryBond, DynamicBond]]]:
        """
        iterate other all bonds
        """
        seen = set()
        for n, m_bond in self._bonds.items():
            seen.add(n)
            for m, bond in m_bond.items():
                if m not in seen:
                    yield n, m, bond

    @cached_property
    def bonds_count(self) -> int:
        return sum(len(x) for x in self._bonds.values()) // 2

    @property
    def meta(self) -> Dict:
        return self.__meta

    @property
    def name(self) -> str:
        return self.__name

    @name.setter
    def name(self, name):
        if not isinstance(name, str):
            raise TypeError('name should be string up to 80 symbols')
        self.__name = name

    @abstractmethod
    def add_atom(self, atom, _map: Optional[int] = None, *, charge: int = 0,
                 is_radical: bool = False, xy: Tuple[float, float] = (0., 0.)) -> int:
        """
        new atom addition
        """
        if _map is None:
            _map = max(self._atoms, default=0) + 1
        elif not isinstance(_map, int):
            raise TypeError('mapping should be integer')
        elif _map in self._atoms:
            raise ValueError('atom with same number exists')

        if not isinstance(xy, tuple) or len(xy) != 2 or not isinstance(xy[0], float) or not isinstance(xy[1], float):
            raise TypeError('XY should be tuple with 2 float')

        charge = self._validate_charge(charge)
        is_radical = self._validate_radical(is_radical)

        self._atoms[_map] = atom
        self._charges[_map] = charge
        self._radicals[_map] = is_radical
        self._plane[_map] = xy
        self._bonds[_map] = {}
        atom._attach_to_graph(self, _map)
        self.__dict__.clear()
        return _map

    @abstractmethod
    def add_bond(self, n: int, m: int, bond):
        """
        new bond addition
        """
        if n == m:
            raise ValueError('atom loops impossible')
        if n not in self._bonds or m not in self._bonds:
            raise AtomNotFound('atoms not found')
        if n in self._bonds[m]:
            raise ValueError('atoms already bonded')

        self._bonds[n][m] = self._bonds[m][n] = bond
        self.__dict__.clear()

    @abstractmethod
    def delete_atom(self, n: int):
        """
        implementation of atom removing
        """
        del self._atoms[n]
        del self._charges[n]
        del self._radicals[n]
        del self._plane[n]
        sb = self._bonds
        for m in sb.pop(n):
            del sb[m][n]
        self.__dict__.clear()

    def delete_bond(self, n: int, m: int):
        """
        implementation of bond removing
        """
        del self._bonds[n][m]
        del self._bonds[m][n]
        self.__dict__.clear()

    @abstractmethod
    def remap(self, mapping: Dict[int, int], *, copy: bool = False):
        if len(mapping) != len(set(mapping.values())) or \
                not (self._atoms.keys() - mapping.keys()).isdisjoint(mapping.values()):
            raise ValueError('mapping overlap')

        mg = mapping.get
        sp = self._plane
        sc = self._charges
        sr = self._radicals

        if copy:
            h = self.__class__()
            h.meta.update(self.__meta)
            h._Graph__name = self.__name
            hb = h._bonds
            ha = h._atoms
            hc = h._charges
            hr = h._radicals
            hp = h._plane
            for n, atom in self._atoms.items():
                m = mg(n, n)
                hc[m] = sc[n]
                hr[m] = sr[n]
                hp[m] = sp[n]
                atom = atom.copy()
                ha[m] = atom
                atom._attach_to_graph(h, m)

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
            return h

        hb = {}
        ha = {}
        hc = {}
        hr = {}
        hp = {}
        for n, atom in self._atoms.items():
            m = mg(n, n)
            hc[m] = sc[n]
            hr[m] = sr[n]
            hp[m] = sp[n]
            ha[m] = atom
            atom._change_map(m)  # change mapping number
        self._atoms = ha
        self._charges = hc
        self._radicals = hr
        self._plane = hp

        for n, m_bond in self._bonds.items():
            hb[mg(n, n)] = {mg(m, m): b for m, b in m_bond.items()}
        self._bonds = hb
        self.__dict__.clear()
        return self

    @abstractmethod
    def copy(self, *, meta: bool = True):
        """
        copy of graph

        :param meta: include metadata
        """
        copy = object.__new__(self.__class__)
        if meta:
            copy._Graph__meta = self.__meta.copy()
            copy._Graph__name = self.__name
        else:
            copy._Graph__meta = {}
            copy._Graph__name = ''

        copy._charges = self._charges.copy()
        copy._radicals = self._radicals.copy()
        copy._plane = self._plane.copy()

        copy._bonds = cb = {}
        for n, m_bond in self._bonds.items():
            cb[n] = cbn = {}
            for m, bond in m_bond.items():
                if m in cb:  # bond partially exists. need back-connection.
                    cbn[m] = cb[m][n]
                else:
                    cbn[m] = bond.copy()

        copy._atoms = ca = {}
        for n, atom in self._atoms.items():
            atom = atom.copy()
            ca[n] = atom
            atom._attach_to_graph(copy, n)
        return copy

    @abstractmethod
    def substructure(self, atoms: Iterable[int], *, meta: bool = False, graph_type: Type['Graph'] = None,
                     atom_type: Type[AnyAtom] = None,
                     bond_type: Type[Union[Bond, DynamicBond, QueryBond]] = None):
        if not atoms:
            raise ValueError('empty atoms list not allowed')
        if set(atoms) - self._atoms.keys():
            raise ValueError('invalid atom numbers')
        atoms = tuple(n for n in self._atoms if n in atoms)  # save original order
        sub = object.__new__(graph_type)

        sa = self._atoms
        sc = self._charges
        sr = self._radicals
        sp = self._plane
        sb = self._bonds

        if meta:
            sub._Graph__meta = self.__meta.copy()
            sub._Graph__name = self.__name
        else:
            sub._Graph__meta = {}
            sub._Graph__name = ''

        sub._charges = {n: sc[n] for n in atoms}
        sub._radicals = {n: sr[n] for n in atoms}
        sub._plane = {n: sp[n] for n in atoms}

        sub._atoms = ca = {}
        for n in atoms:
            atom = sa[n]
            atom = atom_type.from_atom(atom)
            ca[n] = atom
            atom._attach_to_graph(sub, n)

        sub._bonds = cb = {}
        for n in atoms:
            cb[n] = cbn = {}
            for m, bond in sb[n].items():
                if m in cb:  # bond partially exists. need back-connection.
                    cbn[m] = cb[m][n]
                elif m in atoms:
                    cbn[m] = bond_type.from_bond(bond)
        return sub, atoms

    def __and__(self, other):
        """
        Substructure of graph with given nodes.
        """
        return self.substructure(other)

    def __sub__(self, other):
        """
        Given nodes excluded substructure of graph.
        """
        atoms = set(other)
        if atoms - self._atoms.keys():
            raise ValueError('invalid atom numbers')
        atoms = self._atoms.keys() - atoms
        if atoms:
            return self.substructure(atoms)
        raise ValueError('full substitution not allowed')

    def augmented_substructure(self, atoms: Iterable[int], deep: int = 1, **kwargs) -> 'Graph':
        """
        create substructure containing atoms and their neighbors

        :param atoms: list of core atoms in graph
        :param deep: number of bonds between atoms and neighbors
        :param meta: copy metadata to each substructure
        :param as_query: return Query object based on graph substructure. for Molecule and CGR only
        """
        return self.substructure(self._augmented_substructure(atoms, deep)[-1], **kwargs)

    def augmented_substructures(self, atoms: Iterable[int], deep: int = 1, **kwargs) -> List['Graph']:
        """
        create list of substructures containing atoms and their neighbors

        :param atoms: list of core atoms in graph
        :param deep: number of bonds between atoms and neighbors
        :param meta: copy metadata to each substructure
        :param as_query: return Query object based on graph substructure. for Molecule and CGR only
        :return: list of graphs containing atoms, atoms + first circle, atoms + 1st + 2nd,
            etc up to deep or while new nodes available
        """
        return [self.substructure(a, **kwargs) for a in self._augmented_substructure(atoms, deep)]

    def union(self, other: 'Graph', *, remap=False,
              atom_type: Type[AnyAtom] = None, bond_type: Type[Union[Bond, DynamicBond, QueryBond]] = None):
        """
        Merge Graphs into one.

        :param remap: if atoms has collisions then remap other graph atoms else raise exception.
        """
        if self._atoms.keys() & other._atoms.keys():
            if remap:
                other = other.remap({n: i for i, n in enumerate(other, start=max(self._atoms) + 1)}, copy=True)
            else:
                raise ValueError('mapping of graphs is not disjoint')

        u = self.copy(meta=False)
        u._charges.update(other._charges)
        u._radicals.update(other._radicals)
        u._plane.update(other._plane)

        ua = u._atoms
        for n, atom in other._atoms.items():
            atom = atom_type.from_atom(atom)
            ua[n] = atom
            atom._attach_to_graph(u, n)

        ub = u._bonds
        for n, m_bond in other._bonds.items():
            ub[n] = ubn = {}
            for m, bond in m_bond.items():
                if m in ub:  # bond partially exists. need back-connection.
                    ubn[m] = ub[m][n]
                else:
                    ubn[m] = bond_type.from_bond(bond)
        return u, other

    def __or__(self, other):
        """
        G | H is union of graphs
        """
        return self.union(other)

    def flush_cache(self):
        self.__dict__.clear()

    @staticmethod
    def _validate_charge(charge):
        if not isinstance(charge, int):
            raise TypeError('formal charge should be int in range [-4, 4]')
        if charge > 4 or charge < -4:
            raise ValueError('formal charge should be in range [-4, 4]')
        return charge

    @staticmethod
    def _validate_radical(is_radical):
        if not isinstance(is_radical, bool):
            raise TypeError('radical state should be bool')
        return is_radical


__all__ = ['Graph']
