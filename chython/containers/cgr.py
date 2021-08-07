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
from CachedMethods import cached_args_method
from collections import defaultdict
from typing import List, Union, Tuple, Dict, Optional
from . import molecule  # cyclic imports resolve
from .bonds import Bond, DynamicBond
from .graph import Graph
from ..algorithms.calculate2d import Calculate2DCGR
from ..algorithms.components import CGRComponents
from ..algorithms.depict import DepictCGR
from ..algorithms.smiles import CGRSmiles
from ..algorithms.x3dom import X3domCGR
from ..exceptions import MappingError
from ..periodictable import DynamicElement, Element


class CGRContainer(Graph, CGRSmiles, CGRComponents, DepictCGR, Calculate2DCGR, X3domCGR):
    __slots__ = ('_conformers', '_p_charges', '_p_radicals', '_hybridizations', '_p_hybridizations')

    def __init__(self):
        super().__init__()
        self._atoms: Dict[int, DynamicElement] = {}
        self._bonds: Dict[int, Dict[int, DynamicBond]] = {}
        self._conformers: List[Dict[int, Tuple[float, float, float]]] = []
        self._p_charges: Dict[int, int] = {}
        self._p_radicals: Dict[int, bool] = {}
        self._hybridizations: Dict[int, int] = {}
        self._p_hybridizations: Dict[int, int] = {}

    def add_atom(self, atom: Union[DynamicElement, Element, int, str], *args, p_charge: int = 0,
                 p_is_radical: bool = False, **kwargs):
        p_charge = self._validate_charge(p_charge)
        p_is_radical = self._validate_radical(p_is_radical)

        if not isinstance(atom, DynamicElement):
            if isinstance(atom, Element):
                atom = DynamicElement.from_atomic_number(atom.atomic_number)(atom.isotope)
            elif isinstance(atom, str):
                atom = DynamicElement.from_symbol(atom)()
            elif isinstance(atom, int):
                atom = DynamicElement.from_atomic_number(atom)()
            else:
                raise TypeError('DynamicElement object expected')

        _map = super().add_atom(atom, *args, **kwargs)
        self._p_charges[_map] = p_charge
        self._p_radicals[_map] = p_is_radical
        self._hybridizations[_map] = 1
        self._p_hybridizations[_map] = 1
        self._conformers.clear()  # clean conformers. need full recalculation for new system
        return _map

    def add_bond(self, n, m, bond: Union[DynamicBond, Bond, int]):
        if isinstance(bond, DynamicBond):
            order = bond.order
            p_order = bond.p_order
        elif isinstance(bond, Bond):
            order = p_order = bond.order
            bond = DynamicBond.from_bond(bond)
        else:
            order = p_order = bond
            bond = DynamicBond(order, order)

        super().add_bond(n, m, bond)
        self._conformers.clear()

        if order != 1 or p_order != 1:  # 1 is neutral.
            self._calc_hybridization(n)
            self._calc_hybridization(m)

    def delete_atom(self, n):
        old_bonds = self._bonds[n]  # save bonds
        super().delete_atom(n)

        del self._p_charges[n]
        del self._p_radicals[n]
        del self._hybridizations[n]
        del self._p_hybridizations[n]

        for m in old_bonds:
            self._calc_hybridization(m)
        self._conformers.clear()

    def delete_bond(self, n, m):
        super().delete_bond(n, m)
        self._conformers.clear()
        self._calc_hybridization(n)
        self._calc_hybridization(m)

    @cached_args_method
    def neighbors(self, n: int) -> Tuple[int, int]:
        """number of neighbors atoms excluding any-bonded"""
        s = p = 0
        for b in self._bonds[n].values():
            if b.order is not None:
                if b.order != 8:
                    s += 1
                if b.p_order is not None and b.p_order != 8:
                    p += 1
            elif b.p_order != 8:
                p += 1
        return s, p

    def remap(self, mapping, *, copy=False) -> 'CGRContainer':
        h = super().remap(mapping, copy=copy)
        mg = mapping.get
        spr = self._p_radicals
        sh = self._hybridizations
        sph = self._p_hybridizations

        if copy:
            hpc = h._p_charges
            hpr = h._p_radicals
            hh = h._hybridizations
            hc = h._conformers
            hph = h._p_hybridizations
        else:
            hpc = {}
            hpr = {}
            hh = {}
            hc = []
            hph = {}

        for n, c in self._p_charges.items():
            m = mg(n, n)
            hpc[m] = c
            hpr[m] = spr[n]
            hh[m] = sh[n]
            hph[m] = sph[n]

        hc.extend({mg(n, n): x for n, x in c.items()} for c in self._conformers)

        if copy:
            return h

        self._p_charges = hpc
        self._p_radicals = hpr
        self._hybridizations = hh
        self._conformers = hc
        self._p_hybridizations = hph
        return self

    def copy(self, **kwargs) -> 'CGRContainer':
        copy = super().copy(**kwargs)
        copy._hybridizations = self._hybridizations.copy()
        copy._conformers = [c.copy() for c in self._conformers]
        copy._p_hybridizations = self._p_hybridizations.copy()
        copy._p_radicals = self._p_radicals.copy()
        copy._p_charges = self._p_charges.copy()
        return copy

    def substructure(self, atoms, **kwargs) -> 'CGRContainer':
        """
        create substructure containing atoms from atoms list

        :param atoms: list of atoms numbers of substructure
        :param meta: if True metadata will be copied to substructure
        """
        sub, atoms = super().substructure(atoms, graph_type=self.__class__, atom_type=DynamicElement,
                                          bond_type=DynamicBond, **kwargs)
        spc = self._p_charges
        spr = self._p_radicals
        sub._p_charges = {n: spc[n] for n in atoms}
        sub._p_radicals = {n: spr[n] for n in atoms}
        sub._conformers = [{n: c[n] for n in atoms} for c in self._conformers]
        # recalculate query marks
        sub._hybridizations = {}
        sub._p_hybridizations = {}
        for n in sub._atoms:
            sub._calc_hybridization(n)
        return sub

    def union(self, other, **kwargs) -> 'CGRContainer':
        if isinstance(other, CGRContainer):
            u, other = super().union(other, atom_type=DynamicElement, bond_type=DynamicBond, **kwargs)
            u._conformers.clear()
            u._p_charges.update(other._p_charges)
            u._p_radicals.update(other._p_radicals)
            u._hybridizations.update(other._hybridizations)
            u._p_hybridizations.update(other._p_hybridizations)
            return u
        elif isinstance(other, molecule.MoleculeContainer):
            u, other = super().union(other, atom_type=DynamicElement, bond_type=DynamicBond, **kwargs)
            u._conformers.clear()
            u._p_charges.update(other._charges)
            u._p_radicals.update(other._radicals)
            hyb = {n: other.hybridization(n) for n in other}
            u._hybridizations.update(hyb)
            u._p_hybridizations.update(hyb)
            return u
        else:
            raise TypeError('CGRContainer or MoleculeContainer expected')

    def compose(self, other: Union['molecule.MoleculeContainer', 'CGRContainer']) -> 'CGRContainer':
        """
        compose 2 graphs to CGR

        :param other: Molecule or CGR Container
        :return: CGRContainer
        """
        sa = self._atoms
        sc = self._charges
        sr = self._radicals
        spc = self._p_charges
        spr = self._p_radicals
        sp = self._plane
        sb = self._bonds

        bonds = []
        adj: Dict[int, Dict[int, List[Optional[int]]]] = defaultdict(lambda: defaultdict(lambda: [None, None]))
        h = self.__class__()  # subclasses support
        atoms = h._atoms

        if isinstance(other, molecule.MoleculeContainer):
            oa = other._atoms
            oc = other._charges
            or_ = other._radicals
            op = other._plane
            ob = other._bonds
            common = sa.keys() & other

            for n in sa.keys() - common:  # cleavage atoms
                h.add_atom(sa[n].copy(), n, charge=sc[n], is_radical=sr[n], xy=sp[n], p_charge=spc[n],
                           p_is_radical=spr[n])
                for m, bond in sb[n].items():
                    if m not in atoms:
                        if m in common:  # bond to common atoms is broken bond
                            order = bond.order
                            if order:  # skip formed bond. None>X => None>None
                                bond = object.__new__(DynamicBond)
                                bond._DynamicBond__order, bond._DynamicBond__p_order = order, None
                                bonds.append((n, m, bond))
                        else:
                            bonds.append((n, m, bond))
            for n in other._atoms.keys() - common:  # coupling atoms
                h.add_atom(oa[n], n, charge=oc[n], is_radical=or_[n], xy=op[n], p_charge=oc[n], p_is_radical=or_[n])
                for m, bond in ob[n].items():
                    if m not in atoms:
                        if m in common:  # bond to common atoms is formed bond
                            order = bond.order
                            bond = object.__new__(DynamicBond)
                            bond._DynamicBond__order, bond._DynamicBond__p_order = None, order
                        bonds.append((n, m, bond))
            for n in common:
                an = adj[n]
                for m, bond in ob[n].items():
                    if m in common:
                        an[m][1] = bond.order
                for m, bond in sb[n].items():
                    if m in an or m in common and bond.order:
                        an[m][0] = bond.order
            for n in common:
                san = sa[n]
                if san.atomic_number != oa[n].atomic_number or san.isotope != oa[n].isotope:
                    raise MappingError(f'atoms with number {{{n}}} not equal')
                h.add_atom(san.copy(), n, charge=sc[n], is_radical=sr[n], xy=sp[n], p_charge=oc[n], p_is_radical=or_[n])
                for m, (o1, o2) in adj[n].items():
                    if m not in atoms:
                        bond = object.__new__(DynamicBond)
                        bond._DynamicBond__order, bond._DynamicBond__p_order = o1, o2
                        bonds.append((n, m, bond))
        elif isinstance(other, CGRContainer):
            oa = other._atoms
            oc = other._charges
            or_ = other._radicals
            opc = other._p_charges
            opr = other._p_radicals
            op = other._plane
            ob = other._bonds
            common = sa.keys() & other

            for n in sa.keys() - common:  # cleavage atoms
                h.add_atom(sa[n].copy(), n, charge=sc[n], is_radical=sr[n], xy=sp[n], p_charge=spc[n],
                           p_is_radical=spr[n])
                for m, bond in sb[n].items():
                    if m not in atoms:
                        if m in common:  # bond to common atoms is broken bond
                            order = bond.order
                            if order:  # skip formed bond. None>X => None>None
                                bond = object.__new__(DynamicBond)
                                bond._DynamicBond__order, bond._DynamicBond__p_order = order, None
                                bonds.append((n, m, bond))
                        else:
                            bonds.append((n, m, bond))
            for n in other._atoms.keys() - common:  # coupling atoms
                h.add_atom(oa[n].copy(), n, charge=oc[n], is_radical=or_[n], xy=op[n], p_charge=opc[n],
                           p_is_radical=opr[n])
                for m, bond in ob[n].items():
                    if m not in atoms:
                        if m in common:  # bond to common atoms is formed bond
                            order = bond.p_order
                            if order:  # skip broken bond. X>None => None>None
                                bond = object.__new__(DynamicBond)
                                bond._DynamicBond__order, bond._DynamicBond__p_order = None, order
                                bonds.append((n, m, bond))
                        else:
                            bonds.append((n, m, bond))
            for n in common:
                an = adj[n]
                for m, bond in sb[n].items():
                    if m in common and bond.order:  # skip formed bonds
                        an[m][0] = bond.order
                for m, bond in ob[n].items():
                    if m in an or m in common and bond.p_order:
                        # self has nm bond or other bond not broken
                        an[m][1] = bond.p_order
            for n in common:
                san = sa[n]
                if san.atomic_number != oa[n].atomic_number or san.isotope != oa[n].isotope:
                    raise MappingError(f'atoms with number {{{n}}} not equal')
                h.add_atom(san.copy(), n, charge=sc[n], is_radical=sr[n], xy=sp[n], p_charge=opc[n],
                           p_is_radical=opr[n])
                for m, (o1, o2) in adj[n].items():
                    if m not in atoms:
                        bond = object.__new__(DynamicBond)
                        bond._DynamicBond__order, bond._DynamicBond__p_order = o1, o2
                        bonds.append((n, m, bond))
        else:
            raise TypeError('MoleculeContainer or CGRContainer expected')

        for n, m, bond in bonds:
            h.add_bond(n, m, bond)
        return h

    def __xor__(self, other):
        """
        G ^ H is CGR generation
        """
        return self.compose(other)

    def get_mapping(self, other: 'CGRContainer', **kwargs):
        if isinstance(other, CGRContainer):
            return super().get_mapping(other, **kwargs)
        raise TypeError('CGRContainer expected')

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

    def __invert__(self):
        """
        decompose CGR
        """
        return self.decompose()

    def _calc_hybridization(self, n: int):
        hybridization = p_hybridization = 1
        for bond in self._bonds[n].values():
            order = bond.order
            p_order = bond.p_order
            if order and hybridization != 4:
                if order == 4:
                    hybridization = 4
                elif order == 3:
                    if hybridization != 3:
                        hybridization = 3
                elif order == 2:
                    if hybridization == 2:
                        hybridization = 3
                    elif hybridization == 1:
                        hybridization = 2
            if p_order and p_hybridization != 4:
                if p_order == 4:
                    p_hybridization = 4
                elif p_order == 3:
                    if p_hybridization != 3:
                        p_hybridization = 3
                elif p_order == 2:
                    if p_hybridization == 2:
                        p_hybridization = 3
                    elif p_hybridization == 1:
                        p_hybridization = 2
        self._hybridizations[n] = hybridization
        self._p_hybridizations[n] = p_hybridization

    def __getstate__(self):
        return {'conformers': self._conformers, 'p_charges': self._p_charges, 'p_radicals': self._p_radicals,
                **super().__getstate__()}

    def __setstate__(self, state):
        super().__setstate__(state)
        self._p_charges = state['p_charges']
        self._p_radicals = state['p_radicals']
        self._conformers = state['conformers']
        # restore query marks
        self._hybridizations = {}
        self._p_hybridizations = {}
        for n in state['bonds']:
            self._calc_hybridization(n)


__all__ = ['CGRContainer']
