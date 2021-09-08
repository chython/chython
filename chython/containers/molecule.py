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
from collections import Counter, defaultdict
from functools import cached_property
from itertools import zip_longest
from math import ceil
from numpy import uint, zeros
from struct import pack_into, unpack_from
from typing import Dict, Iterable, List, Optional, Set, Tuple, Union
from weakref import ref
from zlib import compress, decompress
from . import cgr, query  # cyclic imports resolve
from .bonds import Bond, DynamicBond, QueryBond
from .graph import Graph
from ..algorithms.aromatics import Aromatize
from ..algorithms.calculate2d import Calculate2DMolecule
from ..algorithms.depict import DepictMolecule
from ..algorithms.fingerprints import Fingerprints
from ..algorithms.huckel import Huckel
from ..algorithms.mcs import MCS
from ..algorithms.smiles import MoleculeSmiles
from ..algorithms.standardize import Saturation, StandardizeMolecule
from ..algorithms.stereo import MoleculeStereo
from ..algorithms.tautomers import Tautomers
from ..algorithms.x3dom import X3domMolecule
from ..exceptions import MappingError, ValenceError
from ..periodictable import DynamicElement, Element, QueryElement


class MoleculeContainer(MoleculeStereo, Graph[Element, Bond], Aromatize, StandardizeMolecule, MoleculeSmiles,
                        DepictMolecule, Calculate2DMolecule, Fingerprints, Tautomers, MCS, Huckel,
                        Saturation, X3domMolecule):
    __slots__ = ('_conformers', '_atoms_stereo', '_hydrogens', '_cis_trans_stereo', '_allenes_stereo',
                 '_parsed_mapping', '_backup', '__meta', '__name')

    def __init__(self):
        super().__init__()
        self._conformers: List[Dict[int, Tuple[float, float, float]]] = []
        self._hydrogens: Dict[int, Optional[int]] = {}
        self._atoms_stereo: Dict[int, bool] = {}
        self._allenes_stereo: Dict[int, bool] = {}
        self._cis_trans_stereo: Dict[Tuple[int, int], bool] = {}
        self._parsed_mapping: Dict[int, int] = {}
        self.__meta = None
        self.__name = None

    @property
    def meta(self) -> Dict:
        if self.__meta is None:
            self.__meta = {}  # lazy
        return self.__meta

    @property
    def name(self) -> str:
        return self.__name or ''

    @name.setter
    def name(self, name):
        if not isinstance(name, str):
            raise TypeError('name should be string up to 80 symbols')
        self.__name = name

    @cached_args_method
    def environment(self, atom: int, include_bond: bool = True, include_atom: bool = True) -> \
            Tuple[Union[Tuple[int, Bond, Element],
                        Tuple[int, Element],
                        Tuple[int, Bond],
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

    @cached_args_method
    def neighbors(self, n: int) -> int:
        """number of neighbors atoms excluding any-bonded"""
        return sum(b.order != 8 for b in self._bonds[n].values())

    @cached_args_method
    def hybridization(self, n: int) -> int:
        """
        Atom hybridization.

        1 - if atom has zero or only single bonded neighbors, 2 - if has only one double bonded neighbor and any amount
        of single bonded, 3 - if has one triple bonded and any amount of double and single bonded neighbors or
        two and more double bonded and any amount of single bonded neighbors, 4 - if atom in aromatic ring.
        """
        hybridization = 1
        for bond in self._bonds[n].values():
            order = bond.order
            if order == 4:
                return 4
            elif order == 3:
                if hybridization != 3:
                    hybridization = 3
            elif order == 2:
                if hybridization == 1:
                    hybridization = 2
                elif hybridization == 2:
                    hybridization = 3
        return hybridization

    @cached_args_method
    def heteroatoms(self, n: int) -> int:
        """
        Number of neighbored heteroatoms (not carbon or hydrogen)
        """
        atoms = self._atoms
        return sum(atoms[m].atomic_number not in (1, 6) for m in self._bonds[n])

    def implicit_hydrogens(self, n: int) -> Optional[int]:
        """
        Number of implicit hydrogen atoms connected to atom.

        Returns None if count are ambiguous.
        """
        return self._hydrogens[n]

    @cached_args_method
    def explicit_hydrogens(self, n: int) -> int:
        """
        Number of explicit hydrogen atoms connected to atom.

        Take into account any type of bonds with hydrogen atoms.
        """
        atoms = self._atoms
        return sum(atoms[m].atomic_number == 1 for m in self._bonds[n])

    @cached_args_method
    def total_hydrogens(self, n: int) -> int:
        """
        Number of hydrogen atoms connected to atom.

        Take into account any type of bonds with hydrogen atoms.
        """
        return self._hydrogens[n] + self.explicit_hydrogens(n)

    def adjacency_matrix(self, set_bonds=False):
        """
        Adjacency matrix of Graph.

        :param set_bonds: if True set bond orders instead of 1.
        """
        adj = zeros((len(self), len(self)), dtype=uint)
        mapping = {n: x for x, n in enumerate(self._atoms)}
        if set_bonds:
            for n, ms in self._bonds.items():
                n = mapping[n]
                for m, b in ms.items():
                    adj[n, mapping[m]] = int(b)
        else:
            for n, ms in self._bonds.items():
                n = mapping[n]
                for m, b in ms.items():
                    adj[n, mapping[m]] = 1
        return adj

    @cached_property
    def molecular_charge(self) -> int:
        """
        Total charge of molecule
        """
        return sum(self._charges.values())

    @cached_property
    def is_radical(self) -> bool:
        """
        True if at least one atom is radical
        """
        return any(self._radicals.values())

    @cached_property
    def molecular_mass(self) -> float:
        return sum(x.atomic_mass for x in self._atoms.values())

    @cached_property
    def brutto(self) -> Dict[str, int]:
        """Counted atoms dict"""
        return Counter(x.atomic_symbol for x in self._atoms.values())

    def add_atom(self, atom: Union[Element, int, str], *args, charge=0, is_radical=False, **kwargs):
        """
        Add new atom.
        """
        if not isinstance(atom, Element):
            if isinstance(atom, str):
                atom = Element.from_symbol(atom)()
            elif isinstance(atom, int):
                atom = Element.from_atomic_number(atom)()
            else:
                raise TypeError('Element object expected')

        _map = super().add_atom(atom, *args, charge=charge, is_radical=is_radical, **kwargs)
        self._conformers.clear()  # clean conformers. need full recalculation for new system

        if atom.atomic_number != 1:
            try:
                rules = atom.valence_rules(charge, is_radical, 0)
            except ValenceError:
                self._hydrogens[_map] = None
            else:
                self._hydrogens[_map] = rules[0][2]  # first rule without neighbors
        else:
            self._hydrogens[_map] = 0
        return _map

    def add_bond(self, n, m, bond: Union[Bond, int]):
        """
        Connect atoms with bonds.

        For Thiele forms of molecule causes invalidation of internal state.
        Implicit hydrogens marks will not be set if atoms in aromatic rings.
        Call `kekule()` and `thiele()` in sequence to fix marks.
        """
        if not isinstance(bond, Bond):
            bond = Bond(bond)

        super().add_bond(n, m, bond)
        self._conformers.clear()  # clean conformers. need full recalculation for new system

        self._calc_implicit(n)
        self._calc_implicit(m)

        if self._atoms[n].atomic_number != 1 and self._atoms[m].atomic_number != 1:  # not hydrogen
            # fix stereo if formed not to hydrogen bond
            self.fix_stereo()

    def delete_atom(self, n):
        """
        Remove atom.

        For Thiele forms of molecule causes invalidation of internal state.
        Implicit hydrogens marks will not be set if atoms in aromatic rings.
        Call `kekule()` and `thiele()` in sequence to fix marks.
        """
        old_bonds = self._bonds[n]  # save bonds
        isnt_hydrogen = self._atoms[n].atomic_number != 1
        super().delete_atom(n)

        del self._hydrogens[n]
        self._conformers.clear()  # clean conformers. need full recalculation for new system
        try:
            del self._parsed_mapping[n]
        except KeyError:
            pass

        for m in old_bonds:
            self._calc_implicit(m)

        if isnt_hydrogen:  # hydrogen atom not used for stereo coding
            self.fix_stereo()

    def delete_bond(self, n, m):
        """
        Disconnect atoms.

        For Thiele forms of molecule causes invalidation of internal state.
        Implicit hydrogens marks will not be set if atoms in aromatic rings.
        Call `kekule()` and `thiele()` in sequence to fix marks.
        """
        super().delete_bond(n, m)
        self._conformers.clear()  # clean conformers. need full recalculation for new system

        self._calc_implicit(n)
        self._calc_implicit(m)

        if self._atoms[n].atomic_number != 1 and self._atoms[m].atomic_number != 1:
            self.fix_stereo()

    def remap(self, mapping: Dict[int, int], *, copy: bool = False) -> 'MoleculeContainer':
        if len(mapping) != len(set(mapping.values())) or \
                not (self._atoms.keys() - mapping.keys()).isdisjoint(mapping.values()):
            raise ValueError('mapping overlap')

        mg = mapping.get
        sp = self._plane
        sc = self._charges
        sr = self._radicals
        shg = self._hydrogens

        if copy:
            h = self.__class__()
            h._MoleculeContainer__name = self.__name
            if self.__meta is not None:
                h._MoleculeContainer__meta = self.__meta.copy()
            hb = h._bonds
            ha = h._atoms
            hc = h._charges
            hr = h._radicals
            hp = h._plane
            hhg = h._hydrogens
            hcf = h._conformers
            has = h._atoms_stereo
            hal = h._allenes_stereo
            hcs = h._cis_trans_stereo
            hm = h._parsed_mapping

            for n, atom in self._atoms.items():
                m = mg(n, n)
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
        else:
            hb = {}
            ha = {}
            hc = {}
            hr = {}
            hp = {}
            hhg = {}
            hcf = []
            has = {}
            hal = {}
            hcs = {}
            hm = {}
            for n, atom in self._atoms.items():
                m = mg(n, n)
                ha[m] = atom
                atom._change_map(m)  # change mapping number

            for n, m_bond in self._bonds.items():
                hb[mg(n, n)] = {mg(m, m): b for m, b in m_bond.items()}

        for n in self._atoms:
            m = mg(n, n)
            hc[m] = sc[n]
            hr[m] = sr[n]
            hp[m] = sp[n]
            hhg[m] = shg[n]

        hcf.extend({mg(n, n): x for n, x in c.items()} for c in self._conformers)
        for n, m in self._parsed_mapping.items():
            hm[mg(n, n)] = m
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
        self._plane = hp
        self._hydrogens = hhg
        self._conformers = hcf
        self._atoms_stereo = has
        self._allenes_stereo = hal
        self._cis_trans_stereo = hcs
        self._parsed_mapping = hm
        self.flush_cache()
        return self

    def copy(self) -> 'MoleculeContainer':
        copy = super().copy()
        copy._MoleculeContainer__name = self.__name
        if self.__meta is None:
            copy._MoleculeContainer__meta = None
        else:
            copy._MoleculeContainer__meta = self.__meta.copy()
        copy._hydrogens = self._hydrogens.copy()
        copy._parsed_mapping = self._parsed_mapping.copy()
        copy._conformers = [c.copy() for c in self._conformers]
        copy._atoms_stereo = self._atoms_stereo.copy()
        copy._allenes_stereo = self._allenes_stereo.copy()
        copy._cis_trans_stereo = self._cis_trans_stereo.copy()
        return copy

    def union(self, other: 'MoleculeContainer', *, remap=False) -> 'MoleculeContainer':
        """
        :param remap: if atoms has collisions then remap other graph atoms else raise exception.
        """
        if not isinstance(other, MoleculeContainer):
            raise TypeError('MoleculeContainer expected')
        elif self._atoms.keys() & other._atoms.keys():
            if remap:
                other = other.remap({n: i for i, n in enumerate(other, start=max(self._atoms) + 1)}, copy=True)
            else:
                raise MappingError('mapping of graphs is not disjoint')
        u = super().union(other)
        u._MoleculeContainer__name = u._MoleculeContainer__meta = None
        u._conformers.clear()
        u._hydrogens.update(other._hydrogens)
        u._parsed_mapping.update(other._parsed_mapping)
        u._atoms_stereo.update(other._atoms_stereo)
        u._allenes_stereo.update(other._allenes_stereo)
        u._cis_trans_stereo.update(other._cis_trans_stereo)
        return u

    def substructure(self, atoms: Iterable[int], *, as_query: bool = False, recalculate_hydrogens=True,
                     skip_neighbors_marks=False, skip_hybridizations_marks=False, skip_hydrogens_marks=False,
                     skip_rings_sizes_marks=False,) -> Union['MoleculeContainer', 'query.QueryContainer']:
        """
        Create substructure containing atoms from atoms list.

        For Thiele forms of molecule In Molecule substructure causes invalidation of internal state.
        Implicit hydrogens marks will not be set if atoms in aromatic rings.
        Call `kekule()` and `thiele()` in sequence to fix marks.

        :param atoms: list of atoms numbers of substructure
        :param as_query: return Query object based on graph substructure
        :param recalculate_hydrogens: calculate implicit H count in substructure
        :param skip_neighbors_marks: Don't set neighbors count marks on substructured queries
        :param skip_hybridizations_marks: Don't set hybridizations marks on substructured queries
        :param skip_hydrogens_marks: Don't set hydrogens count marks on substructured queries
        :param skip_rings_sizes_marks: Don't set rings_sizes marks on substructured queries
        """
        if not atoms:
            raise ValueError('empty atoms list not allowed')
        if set(atoms) - self._atoms.keys():
            raise ValueError('invalid atom numbers')
        atoms = tuple(n for n in self._atoms if n in atoms)  # save original order
        if as_query:
            atom_type = QueryElement
            bond_type = QueryBond
            sub = object.__new__(query.QueryContainer)
        else:
            atom_type = Element
            bond_type = Bond
            sub = object.__new__(self.__class__)
            sub._MoleculeContainer__name = sub._MoleculeContainer__meta = None

        sa = self._atoms
        sb = self._bonds
        sc = self._charges
        sr = self._radicals
        sp = self._plane

        sub._charges = {n: sc[n] for n in atoms}
        sub._radicals = {n: sr[n] for n in atoms}
        sub._plane = {n: sp[n] for n in atoms}

        sub._atoms = ca = {}
        for n in atoms:
            ca[n] = atom = atom_type.from_atom(sa[n])
            atom._attach_to_graph(sub, n)

        sub._bonds = cb = {}
        for n in atoms:
            cb[n] = cbn = {}
            for m, bond in sb[n].items():
                if m in cb:  # bond partially exists. need back-connection.
                    cbn[m] = cb[m][n]
                elif m in atoms:
                    cbn[m] = bond_type.from_bond(bond)

        if as_query:
            lost = {n for n, a in sa.items() if a.atomic_number != 1} - set(atoms)  # atoms not in substructure
            not_skin = {n for n in atoms if lost.isdisjoint(sb[n])}
            sub._atoms_stereo = {n: s for n, s in self._atoms_stereo.items() if n in not_skin}
            sub._allenes_stereo = {n: s for n, s in self._allenes_stereo.items()
                                   if not_skin.issuperset(self._stereo_allenes_paths[n]) and
                                      not_skin.issuperset(x for x in self._stereo_allenes[n] if x)}
            sub._cis_trans_stereo = {nm: s for nm, s in self._cis_trans_stereo.items()
                                     if not_skin.issuperset(self._stereo_cis_trans_paths[nm]) and
                                        not_skin.issuperset(x for x in self._stereo_cis_trans[nm] if x)}

            sub._heteroatoms = {n: () for n in atoms}

            if skip_hybridizations_marks:
                sub._hybridizations = {n: () for n in atoms}
            else:
                sh = self.hybridization
                sub._hybridizations = {n: (sh(n),) for n in atoms}
            if skip_neighbors_marks:
                sub._neighbors = {n: () for n in atoms}
            else:
                sn = self.neighbors
                sub._neighbors = {n: (sn(n),) for n in atoms}
            if skip_hydrogens_marks:
                sub._hydrogens = {n: () for n in atoms}
            else:
                shg = self._hydrogens
                sub._hydrogens = {n: () if shg[n] is None else (shg[n],) for n in atoms}
            if skip_rings_sizes_marks:
                sub._rings_sizes = {n: () for n in atoms}
            else:
                rs = self.atoms_rings_sizes
                sub._rings_sizes = {n: rs.get(n, ()) for n in atoms}
        else:
            sub._conformers = [{n: c[n] for n in atoms} for c in self._conformers]

            if recalculate_hydrogens:
                sub._hydrogens = {}
                for n in atoms:
                    sub._calc_implicit(n)
            else:
                hg = self._hydrogens
                sub._hydrogens = {n: hg[n] for n in atoms}

            sub._parsed_mapping = {n: m for n, m in self._parsed_mapping.items() if n in atoms}

            # fix_stereo will repair data
            sub._atoms_stereo = self._atoms_stereo.copy()
            sub._allenes_stereo = self._allenes_stereo.copy()
            sub._cis_trans_stereo = self._cis_trans_stereo.copy()
            sub.fix_stereo()
        return sub

    def augmented_substructure(self, atoms: Iterable[int], deep: int = 1, **kwargs) -> 'MoleculeContainer':
        """
        Create substructure containing atoms and their neighbors

        :param atoms: list of core atoms in graph
        :param deep: number of bonds between atoms and neighbors
        """
        return self.substructure(self._augmented_substructure(atoms, deep)[-1], **kwargs)

    def augmented_substructures(self, atoms: Iterable[int], deep: int = 1, **kwargs) -> List['MoleculeContainer']:
        """
        Create list of substructures containing atoms and their neighbors

        :param atoms: list of core atoms in graph
        :param deep: number of bonds between atoms and neighbors
        :return: list of graphs containing atoms, atoms + first circle, atoms + 1st + 2nd,
            etc up to deep or while new nodes available
        """
        return [self.substructure(a, **kwargs) for a in self._augmented_substructure(atoms, deep)]

    def split(self) -> List['MoleculeContainer']:
        """
        Split disconnected structure to connected substructures
        """
        return [self.substructure(c, recalculate_hydrogens=False) for c in self.connected_components]

    def compose(self, other: 'MoleculeContainer') -> 'cgr.CGRContainer':
        """
        Compose 2 graphs to CGR.
        """
        if not isinstance(other, MoleculeContainer):
            raise TypeError('MoleculeContainer expected')
        sa = self._atoms
        sc = self._charges
        sr = self._radicals
        sp = self._plane
        sb = self._bonds

        bonds = []
        adj = defaultdict(lambda: defaultdict(lambda: [None, None]))

        oa = other._atoms
        oc = other._charges
        or_ = other._radicals
        op = other._plane
        ob = other._bonds

        common = sa.keys() & oa.keys()

        h = cgr.CGRContainer()
        ha = h._atoms
        hb = h._bonds
        hc = h._charges
        hpc = h._p_charges
        hr = h._radicals
        hpr = h._p_radicals
        hp = h._plane

        for n in sa.keys() - common:  # cleavage atoms
            hc[n] = hpc[n] = sc[n]
            hr[n] = hpr[n] = sr[n]
            hp[n] = sp[n]
            hb[n] = {}
            ha[n] = a = DynamicElement.from_atom(sa[n])
            a._attach_to_graph(h, n)

            for m, bond in sb[n].items():
                if m not in ha:
                    if m in common:  # bond to common atoms is broken bond
                        bond = DynamicBond(bond.order, None)
                    else:
                        bond = DynamicBond(bond.order, bond.order)
                    bonds.append((n, m, bond))
        for n in oa.keys() - common:  # coupling atoms
            hc[n] = hpc[n] = oc[n]
            hr[n] = hpr[n] = or_[n]
            hp[n] = op[n]
            hb[n] = {}
            ha[n] = a = DynamicElement.from_atom(oa[n])
            a._attach_to_graph(h, n)

            for m, bond in ob[n].items():
                if m not in ha:
                    if m in common:  # bond to common atoms is formed bond
                        bond = DynamicBond(None, bond.order)
                    else:
                        bond = DynamicBond(bond.order, bond.order)
                    bonds.append((n, m, bond))
        for n in common:
            an = adj[n]
            for m, bond in sb[n].items():
                if m in common:
                    an[m][0] = bond.order
            for m, bond in ob[n].items():
                if m in common:
                    an[m][1] = bond.order
        for n in common:
            san = sa[n]
            if san.atomic_number != oa[n].atomic_number or san.isotope != oa[n].isotope:
                raise MappingError(f'atoms with number {{{n}}} not equal')

            hc[n] = sc[n]
            hpc[n] = oc[n]
            hr[n] = sr[n]
            hpr[n] = or_[n]
            hp[n] = sp[n]
            hb[n] = {}
            ha[n] = a = DynamicElement.from_atom(san)
            a._attach_to_graph(h, n)

            for m, (o1, o2) in adj[n].items():
                if m not in ha:
                    bonds.append((n, m, DynamicBond(o1, o2)))

        for n, m, bond in bonds:
            hb[n][m] = hb[m][n] = bond
        return h

    def get_fast_mapping(self, other: 'MoleculeContainer') -> Optional[Dict[int, int]]:
        """
        Get self to other fast (suboptimal) structure mapping.
        Only one possible atoms mapping returned.
        Effective only for big molecules.
        """
        if isinstance(other, MoleculeContainer):
            if len(self) != len(other):
                return
            so = self.smiles_atoms_order
            oo = other.smiles_atoms_order
            if self != other:
                return
            return dict(zip(so, oo))
        raise TypeError('MoleculeContainer expected')

    def get_mapping(self, other: 'MoleculeContainer', **kwargs):
        if isinstance(other, MoleculeContainer):
            return super().get_mapping(other, **kwargs)
        raise TypeError('MoleculeContainer expected')

    def pack(self) -> bytes:
        """
        Pack into compressed bytes.
        Note:
            * Less than 4096 atoms supported. Atoms mapping should be in range 1-4095.
            * Implicit hydrogens count should be in range 0-7
            * Isotope shift should be in range -15 - 15 relatively mdl.common_isotopes
            * Atoms neighbors should be in range 0-15

        Format specification:
        Big endian bytes order
        8 bit - empty byte for future extending
        12 bit - number of atoms
        12 bit - cis/trans stereo block size
        Atom block 9 bytes (repeated):
        12 bit - atom number
        4 bit - number of neighbors
        2 bit tetrahedron sign (00 - not stereo, 10 or 11 - has stereo)
        2 bit - allene sign
        5 bit - isotope (00000 - not specified, over = isotope - common_isotope + 16)
        7 bit - atomic number (<=118)
        32 bit - XY float16 coordinates
        3 bit - hydrogens (0-7)
        4 bit - charge (charge + 4. possible range -4 - 4)
        1 bit - radical state
        Connection table: flatten list of neighbors. neighbors count stored in atom block.
        For example CC(=O)O - {1: [2], 2: [1, 3, 4], 3: [2], 4: [2]} >> [2, 1, 3, 4, 2, 2].
        Repeated block (equal to bonds count).
        24 bit - paired 12 bit numbers.
        Bonds order block (repeated):
        16 bit - 5 bonds grouped (3 bit each). 1 bit unused. Zero padding used than bonds count not proportional to 5.
        Cis/trans data block (repeated):
        24 bit - atoms pair
        7 bit - zero padding. in future can be used for extra bond-level stereo, like atropoisomers.
        1 bit - sign
        """
        bonds = self._bonds
        if max(bonds) > 4095:
            raise ValueError('Big molecules not supported')
        if any(len(x) > 15 for x in bonds.values()):
            raise ValueError('To many neighbors not supported')
        from ..files._mdl.mol import common_isotopes

        plane = self._plane
        charges = self._charges
        radicals = self._radicals
        hydrogens = self._hydrogens
        atoms_stereo = self._atoms_stereo
        allenes_stereo = self._allenes_stereo
        cis_trans_stereo = self._cis_trans_stereo

        data = bytearray(4 +  # extension byte + atoms count + cis/trans bit
                         9 * self.atoms_count +  # atoms data
                         3 * self.bonds_count +  # connection table
                         2 * ceil(self.bonds_count / 5) +  # bonds order
                         4 * len(cis_trans_stereo))
        pack_into('>HB', data, 1, (self.atoms_count << 4) | (len(cis_trans_stereo) >> 8), len(cis_trans_stereo) & 0xff)
        shift = 4

        neighbors = []
        bonds_pack = []
        seen = set()
        hold = []
        for o, (n, a) in enumerate(self._atoms.items()):
            bs = bonds[n]
            neighbors.extend(bs)
            seen.add(n)
            for m, b in bs.items():
                if m not in seen:
                    bonds_pack.append(b.order - 1)  # 8 - 4 bit, but 7 - 3 bit

            # 3 bit - hydrogens (0-7) | 4 bit - charge | 1 bit - radical
            hcr = (charges[n] + 4) << 1
            if radicals[n]:
                hcr |= 1
            hcr |= hydrogens[n] << 5

            # 2 bit tetrahedron sign | 2 bit - allene sign | 5 bit - isotope | 7 bit - atomic number (<=118)
            sia = a.atomic_number
            if a.isotope:
                sia |= (a.isotope - common_isotopes[a.atomic_symbol] + 16) << 7

            if n in atoms_stereo:
                if atoms_stereo[n]:
                    sia |= 0xc000
                else:
                    sia |= 0x8000
            if n in allenes_stereo:
                if allenes_stereo[n]:
                    sia |= 0x3000
                else:
                    sia |= 0x2000

            hold.append((shift + 9 * o, (n << 4) | len(bs), sia, *plane[n], hcr))

        shift += 9 * self.atoms_count + 3 * self.bonds_count - 4
        ngb = iter(reversed(neighbors))
        for o, (n2, n1) in enumerate(zip_longest(ngb, ngb)):
            # 12 bit + 12 bit
            pack_into('>I', data, shift - 3 * o, (n1 << 12) | n2)

        # pack after connection table for preventing override!
        for x in hold:
            pack_into('>2H2eB', data, *x)

        # 16 bit - 5 bonds packing. 1 bit empty.
        shift += 4
        bp = iter(bonds_pack)
        for o, (b1, b2, b3, b4, b5) in enumerate(zip_longest(bp, bp, bp, bp, bp, fillvalue=0)):
            pack_into('>H', data, shift + 2 * o, (b1 << 12) | (b2 << 9) | (b3 << 6) | (b4 << 3) | b5)

        shift += 2 * ceil(self.bonds_count / 5)
        for o, ((n, m), s) in enumerate(cis_trans_stereo.items()):
            pack_into('>I', data, shift + 4 * o, (n << 20) | (m << 8) | s)

        return compress(bytes(data), 9)

    @classmethod
    def unpack(cls, data: bytes) -> 'MoleculeContainer':
        """
        Unpack from compressed bytes.
        """
        try:  # windows? ;)
            from ._unpack import unpack
        except ImportError:
            return cls.pure_unpack(data)
        (mapping, atom_numbers, isotopes, charges, radicals, hydrogens, plane, bonds,
         atoms_stereo, allenes_stereo, cis_trans_stereo) = unpack(decompress(data))

        mol = object.__new__(cls)
        mol._bonds = bonds
        mol._plane = plane
        mol._charges = charges
        mol._radicals = radicals
        mol._hydrogens = hydrogens
        mol._atoms_stereo = atoms_stereo
        mol._allenes_stereo = allenes_stereo
        mol._cis_trans_stereo = cis_trans_stereo

        mol._conformers = []
        mol._parsed_mapping = {}
        mol._MoleculeContainer__meta = None
        mol._MoleculeContainer__name = None
        mol._atoms = atoms = {}

        for n, a, i in zip(mapping, atom_numbers, isotopes):
            atoms[n] = a = object.__new__(Element.from_atomic_number(a))
            a._Core__isotope = i
            a._graph = ref(mol)
            a._map = n
        return mol

    @classmethod
    def pure_unpack(cls, data: bytes) -> 'MoleculeContainer':
        """
        Unpack from compressed bytes. Python implementation.
        """
        from ..files._mdl.mol import common_isotopes

        data = memoryview(decompress(data))
        mol = cls()
        atoms = mol._atoms
        bonds = mol._bonds
        plane = mol._plane
        charges = mol._charges
        radicals = mol._radicals
        hydrogens = mol._hydrogens
        atoms_stereo = mol._atoms_stereo
        allenes_stereo = mol._allenes_stereo
        cis_trans_stereo = mol._cis_trans_stereo

        neighbors = {}
        acs = int.from_bytes(data[1:4], 'big')
        shift = 4
        for o in range(acs >> 12):
            nn, sia, x, y, hcr = unpack_from('>2H2eB', data, shift + 9 * o)
            n = nn >> 4
            neighbors[n] = nn & 0x0f
            # stereo
            s = sia >> 14
            if s:
                atoms_stereo[n] = s == 3
            s = (sia >> 12) & 3
            if s:
                allenes_stereo[n] = s == 3

            # atoms
            a = Element.from_atomic_number(sia & 0x7f)
            ai = (sia >> 7) & 0x1f
            if ai:
                ai += common_isotopes[a.__name__] - 16
            else:
                ai = None
            atoms[n] = a = a(ai)
            a._attach_to_graph(mol, n)

            charges[n] = ((hcr >> 1) & 0x0f) - 4
            radicals[n] = bool(hcr & 0x01)
            hydrogens[n] = hcr >> 5
            plane[n] = (x, y)

        bc = sum(neighbors.values()) // 2
        shift += 9 * len(neighbors) - 1
        connections = []
        for o in range(bc):
            nn = unpack_from('>I', data, shift + 3 * o)[0]
            connections.append((nn >> 12) & 0x0fff)
            connections.append(nn & 0x0fff)

        shift += 1 + 3 * bc
        orders = []
        for o in range(ceil(bc / 5)):
            bb = unpack_from('>H', data, shift + 2 * o)[0]
            orders.append(((bb >> 12) & 0x07) + 1)
            orders.append(((bb >> 9) & 0x07) + 1)
            orders.append(((bb >> 6) & 0x07) + 1)
            orders.append(((bb >> 3) & 0x07) + 1)
            orders.append((bb & 0x07) + 1)
        orders = orders[:bc]  # skip padding

        con = iter(connections)
        ords = iter(orders)
        for n, ms in neighbors.items():
            bonds[n] = cbn = {}
            for _ in range(ms):
                m = next(con)
                if m in bonds:  # bond partially exists. need back-connection.
                    cbn[m] = bonds[m][n]
                else:
                    cbn[m] = Bond(next(ords))

        shift += 2 * ceil(bc / 5)
        for o in range(acs & 0x0fff):  # cis/trans
            ct = unpack_from('>I', data, shift + 4 * o)[0]
            cis_trans_stereo[(ct >> 20, (ct >> 8) & 0x0fff)] = ct & 0x01
        return mol

    def _augmented_substructure(self, atoms: Iterable[int], deep: int):
        atoms = set(atoms)
        bonds = self._bonds
        if atoms - self._atoms.keys():
            raise ValueError('invalid atom numbers')
        nodes = [atoms]
        for _ in range(deep):
            n = {y for x in nodes[-1] for y in bonds[x]} | nodes[-1]
            if n in nodes:
                break
            nodes.append(n)
        return nodes

    @cached_property
    def _screen_fingerprint(self) -> Dict[int, Set[int]]:
        """
        Fingerprint of available linear fragments with set of mapped atoms.
        Required for isomorphism tests filtering speedup.
        Parameters can be modified globally in `MoleculeContainer._fingerprint_config`.
        """
        if self._fingerprint_config:
            return {hash(k): {x for x in v for x in x} for k, v in self._fragments(**self._fingerprint_config).items()}
        return {}

    @cached_args_method
    def _component_fingerprint(self, component):
        """
        Fingerprint of specific component.
        """
        scope = set(self.connected_components[component])
        return {k: v & scope for k, v in self._screen_fingerprint.items() if not v.isdisjoint(scope)}

    def _calc_implicit(self, n: int):
        atoms = self._atoms
        atom = atoms[n]
        if atom.atomic_number != 1:
            charge: int = self._charges[n]
            is_radical = self._radicals[n]
            explicit_sum = 0
            explicit_dict = defaultdict(int)
            for m, bond in self._bonds[n].items():
                order = bond.order
                if order == 4:  # aromatic rings not supported
                    self._hydrogens[n] = None
                    return
                elif order != 8:  # any bond used for complexes
                    explicit_sum += order
                    explicit_dict[(order, atoms[m].atomic_number)] += 1
            try:
                rules = atom.valence_rules(charge, is_radical, explicit_sum)
            except ValenceError:
                self._hydrogens[n] = None
                return
            for s, d, h in rules:
                if s.issubset(explicit_dict) and all(explicit_dict[k] >= c for k, c in d.items()):
                    self._hydrogens[n] = h
                    return
        self._hydrogens[n] = 0

    def __int__(self):
        """
        Total charge of molecule
        """
        return self.molecular_charge

    def __float__(self):
        return self.molecular_mass

    def __xor__(self, other):
        """
        G ^ H is CGR generation
        """
        return self.compose(other)

    def __and__(self, other: Iterable[int]):
        """
        Substructure of graph with given nodes.
        """
        return self.substructure(other)

    def __sub__(self, other: Iterable[int]):
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

    def __enter__(self):
        """
        Transaction of changes. Keep current state for restoring on errors.
        """
        atoms = {}
        for n, atom in self._atoms.items():
            atom = atom.copy()
            atoms[n] = atom
            atom._attach_to_graph(self, n)

        bonds = {}
        for n, m_bond in self._bonds.items():
            bonds[n] = cbn = {}
            for m, bond in m_bond.items():
                if m in bonds:  # bond partially exists. need back-connection.
                    cbn[m] = bonds[m][n]
                else:
                    cbn[m] = bond.copy()
        self._backup = {'atoms': atoms, 'bonds': bonds, 'parsed_mapping': self._parsed_mapping.copy(),
                        'plane': self._plane.copy(), 'charges': self._charges.copy(), 'radicals': self._radicals.copy(),
                        'hydrogens': self._hydrogens.copy(), 'conformers': [x.copy() for x in self._conformers],
                        'atoms_stereo': self._atoms_stereo.copy(), 'allenes_stereo': self._allenes_stereo.copy(),
                        'cis_trans_stereo': self._cis_trans_stereo.copy()}
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if exc_type:  # restore state
            backup = self._backup
            self._atoms = backup['atoms']
            self._bonds = backup['bonds']
            self._parsed_mapping = backup['parsed_mapping']
            self._plane = backup['plane']
            self._charges = backup['charges']
            self._radicals = backup['radicals']
            self._hydrogens = backup['hydrogens']
            self._conformers = backup['conformers']
            self._atoms_stereo = backup['atoms_stereo']
            self._allenes_stereo = backup['allenes_stereo']
            self._cis_trans_stereo = backup['cis_trans_stereo']
            self.flush_cache()
        del self._backup

    def __getstate__(self):
        return {'conformers': self._conformers, 'hydrogens': self._hydrogens, 'atoms_stereo': self._atoms_stereo,
                'allenes_stereo': self._allenes_stereo, 'cis_trans_stereo': self._cis_trans_stereo,
                'parsed_mapping': self._parsed_mapping, 'meta': self.__meta, 'name': self.__name,
                **super().__getstate__()}

    def __setstate__(self, state):
        super().__setstate__(state)
        self._conformers = state['conformers']
        self._atoms_stereo = state['atoms_stereo']
        self._allenes_stereo = state['allenes_stereo']
        self._cis_trans_stereo = state['cis_trans_stereo']
        self._hydrogens = state['hydrogens']
        self._parsed_mapping = state['parsed_mapping']
        self.__meta = state['meta']
        self.__name = state['name']

    _fingerprint_config = {'min_radius': 2, 'max_radius': 4}  # set empty for disable screening


__all__ = ['MoleculeContainer']
