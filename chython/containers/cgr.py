# -*- coding: utf-8 -*-
#
#  Copyright 2017-2021 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from functools import cached_property
from typing import Dict, Iterator, List, Tuple
from . import molecule  # cyclic imports resolve
from .bonds import DynamicBond
from ..algorithms.calculate2d import Calculate2DCGR
from ..algorithms.depict import DepictCGR
from ..algorithms.isomorphism import Isomorphism
from ..algorithms.morgan import Morgan
from ..algorithms.rings import Rings
from ..algorithms.smiles import CGRSmiles
from ..algorithms.x3dom import X3domCGR
from ..periodictable import DynamicElement, Element


class CGRContainer(CGRSmiles, Morgan, Rings, Isomorphism, DepictCGR, Calculate2DCGR, X3domCGR):
    __slots__ = ('_atoms', '_bonds', '_charges', '_radicals', '_p_charges', '_p_radicals', '_plane', '_conformers',
                 '__dict__', '__weakref__')
    _atoms: Dict[int, DynamicElement]
    _bonds: Dict[int, Dict[int, DynamicBond]]
    _charges: Dict[int, int]
    _radicals: Dict[int, bool]
    _p_charges: Dict[int, int]
    _p_radicals: Dict[int, bool]
    _plane: Dict[int, Tuple[float, float]]
    _conformers: List[Dict[int, Tuple[float, float, float]]]

    def __init__(self):
        super().__init__()
        self._atoms = {}
        self._bonds = {}
        self._charges = {}
        self._radicals = {}
        self._p_charges = {}
        self._p_radicals = {}
        self._plane = {}
        self._conformers = []

    def bonds(self) -> Iterator[Tuple[int, int, DynamicBond]]:
        """
        Iterate other all bonds
        """
        seen = set()
        for n, m_bond in self._bonds.items():
            seen.add(n)
            for m, bond in m_bond.items():
                if m not in seen:
                    yield n, m, bond

    @cached_property
    def center_atoms(self) -> Tuple[int, ...]:
        """ Get list of atoms of reaction center (atoms with dynamic: bonds, charges, radicals).
        """
        radicals = self._radicals
        p_charges = self._p_charges
        p_radicals = self._p_radicals

        center = set()
        for n, c in self._charges.items():
            if c != p_charges[n] or radicals[n] != p_radicals[n]:
                center.add(n)

        for n, m_bond in self._bonds.items():
            if any(bond.order != bond.p_order for bond in m_bond.values()):
                center.add(n)

        return tuple(center)

    @cached_property
    def aromatic_rings(self) -> Tuple[Tuple[int, ...], ...]:
        """
        existed or formed aromatic rings atoms numbers
        """
        adj = self._bonds
        return tuple(ring for ring in self.sssr if
                     adj[ring[0]][ring[-1]].order == 4 and all(adj[n][m].order == 4 for n, m in zip(ring, ring[1:])) or
                     adj[ring[0]][ring[-1]].p_order == 4 and all(adj[n][m].p_order == 4 for n, m in zip(ring, ring[1:]))
                     )

    def substructure(self, atoms) -> 'CGRContainer':
        """
        Create substructure containing atoms from atoms list

        :param atoms: list of atoms numbers of substructure
        """
        atoms = set(atoms)
        sa = self._atoms
        sc = self._charges
        sr = self._radicals
        sp = self._plane
        sb = self._bonds
        spc = self._p_charges
        spr = self._p_radicals

        sub = object.__new__(self.__class__)
        sub._charges = {n: sc[n] for n in atoms}
        sub._radicals = {n: sr[n] for n in atoms}
        sub._p_charges = {n: spc[n] for n in atoms}
        sub._p_radicals = {n: spr[n] for n in atoms}
        sub._plane = {n: sp[n] for n in atoms}
        sub._conformers = [{n: c[n] for n in atoms} for c in self._conformers]

        sub._atoms = ca = {}
        for n in atoms:
            ca[n] = atom = sa[n].copy()
            atom._attach_to_graph(sub, n)

        sub._bonds = cb = {}
        for n in atoms:
            cb[n] = cbn = {}
            for m, bond in sb[n].items():
                if m in cb:  # bond partially exists. need back-connection.
                    cbn[m] = cb[m][n]
                elif m in atoms:
                    cbn[m] = bond.copy()
        return sub

    def decompose(self) -> Tuple['molecule.MoleculeContainer', 'molecule.MoleculeContainer']:
        """
        decompose CGR to pair of Molecules, which represents reactants and products state of reaction

        :return: tuple of two molecules
        """
        charges = self._charges
        p_charges = self._p_charges
        radicals = self._radicals
        p_radicals = self._p_radicals
        plane = self._plane

        reactants = molecule.MoleculeContainer()
        products = molecule.MoleculeContainer()

        for n, atom in self._atoms.items():
            atom = Element.from_atomic_number(atom.atomic_number)(atom.isotope)
            reactants.add_atom(atom, n, charge=charges[n], is_radical=radicals[n], xy=plane[n])
            products.add_atom(atom.copy(), n, charge=p_charges[n], is_radical=p_radicals[n], xy=plane[n])

        for n, m, bond in self.bonds():
            if bond.order:
                reactants.add_bond(n, m, bond.order)
            if bond.p_order:
                products.add_bond(n, m, bond.p_order)
        return reactants, products

    def get_mapping(self, other: 'CGRContainer', **kwargs):
        if isinstance(other, CGRContainer):
            return super().get_mapping(other, **kwargs)
        raise TypeError('CGRContainer expected')

    def __invert__(self):
        """
        decompose CGR
        """
        return self.decompose()

    def __getstate__(self):
        return {'atoms': self._atoms, 'bonds': self._bonds, 'plane': self._plane, 'conformers': self._conformers,
                'charges': self._charges, 'radicals': self._radicals,
                'p_charges': self._p_charges, 'p_radicals': self._p_radicals}

    def __setstate__(self, state):
        self._atoms = state['atoms']
        for n, a in state['atoms'].items():
            a._attach_to_graph(self, n)
        self._charges = state['charges']
        self._radicals = state['radicals']
        self._plane = state['plane']
        self._bonds = state['bonds']
        self._p_charges = state['p_charges']
        self._p_radicals = state['p_radicals']
        self._conformers = state['conformers']


__all__ = ['CGRContainer']
