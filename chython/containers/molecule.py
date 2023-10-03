# -*- coding: utf-8 -*-
#
#  Copyright 2017-2023 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from numpy import uint, zeros
from typing import Dict, Iterable, List, Optional, Tuple, Union
from weakref import ref
from zlib import compress, decompress
from .bonds import Bond, DynamicBond, QueryBond
from .cgr import CGRContainer
from .graph import Graph
from .query import QueryContainer
from ..algorithms.aromatics import Aromatize
from ..algorithms.calculate2d import Calculate2DMolecule
from ..algorithms.depict import DepictMolecule
from ..algorithms.isomorphism import MoleculeIsomorphism
from ..algorithms.fingerprints import Fingerprints
from ..algorithms.mcs import MCS
from ..algorithms.smiles import MoleculeSmiles
from ..algorithms.standardize import StandardizeMolecule
from ..algorithms.stereo import MoleculeStereo
from ..algorithms.tautomers import Tautomers
from ..algorithms.x3dom import X3domMolecule
from ..exceptions import MappingError, ValenceError
from ..periodictable import DynamicElement, Element, QueryElement, H


class MoleculeContainer(MoleculeStereo, Graph[Element, Bond], MoleculeIsomorphism, Aromatize, StandardizeMolecule,
                        MoleculeSmiles, DepictMolecule, Calculate2DMolecule, Fingerprints, Tautomers, MCS,
                        X3domMolecule):
    __slots__ = ('_plane', '_conformers', '_hydrogens', '_parsed_mapping', '_backup', '__meta', '__name')

    _conformers: List[Dict[int, Tuple[float, float, float]]]
    _hydrogens: Dict[int, Optional[int]]
    _parsed_mapping: Dict[int, int]
    _plane: Dict[int, Tuple[float, float]]

    def __init__(self):
        super().__init__()
        self._conformers = []
        self._hydrogens = {}
        self._parsed_mapping = {}
        self._plane = {}
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
        Number of neighbored heteroatoms (not carbon or hydrogen) except any-bond connected.
        """
        atoms = self._atoms
        return sum(atoms[m].atomic_number not in (1, 6) for m, b in self._bonds[n].items() if b.order != 8)

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

    @cached_args_method
    def adjacency_matrix(self, set_bonds=False, /):
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
        return sum(x.atomic_mass for x in self._atoms.values()) + sum(self._hydrogens.values()) * H().atomic_mass

    @cached_property
    def brutto(self) -> Dict[str, int]:
        """Counted atoms dict"""
        c = Counter(x.atomic_symbol for x in self._atoms.values())
        c['H'] += sum(self._hydrogens.values())
        return dict(c)

    @cached_property
    def aromatic_rings(self) -> Tuple[Tuple[int, ...], ...]:
        """
        Aromatic rings atoms numbers
        """
        bonds = self._bonds
        return tuple(ring for ring in self.sssr if bonds[ring[0]][ring[-1]] == 4
                     and all(bonds[n][m] == 4 for n, m in zip(ring, ring[1:])))

    def add_atom(self, atom: Union[Element, int, str], *args, charge=0, is_radical=False,
                 xy: Tuple[float, float] = (0., 0.), _skip_hydrogen_calculation=False, **kwargs):
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
        if not isinstance(xy, tuple) or len(xy) != 2 or not isinstance(xy[0], float) or not isinstance(xy[1], float):
            raise TypeError('XY should be tuple with 2 float')

        n = super().add_atom(atom, *args, charge=charge, is_radical=is_radical, **kwargs)
        self._plane[n] = xy
        self._conformers.clear()  # clean conformers. need full recalculation for new system

        if _skip_hydrogen_calculation:
            self._hydrogens[n] = None
        elif atom.atomic_number != 1:
            try:
                rules = atom.valence_rules(charge, is_radical, 0)
            except ValenceError:
                self._hydrogens[n] = None
            else:
                self._hydrogens[n] = rules[0][2]  # first rule without neighbors
        else:
            self._hydrogens[n] = 0
        return n

    def add_bond(self, n, m, bond: Union[Bond, int], *, _skip_hydrogen_calculation=False):
        """
        Connect atoms with bonds.

        For Thiele forms of molecule causes invalidation of internal state.
        Implicit hydrogens marks will not be set if atoms in aromatic rings.
        Call `kekule()` and `thiele()` in sequence to fix marks.
        """
        if not isinstance(bond, Bond):
            bond = Bond(bond)

        bond._attach_graph(self, n, m)
        super().add_bond(n, m, bond)
        self._conformers.clear()  # clean conformers. need full recalculation for new system

        if _skip_hydrogen_calculation:  # skip stereo fixing too
            return

        self._calc_implicit(n)
        self._calc_implicit(m)

        if self._atoms[n].atomic_number != 1 and self._atoms[m].atomic_number != 1:  # not hydrogen
            # fix stereo if formed not to hydrogen bond
            self.fix_stereo()

    def delete_atom(self, n: int, *, _skip_hydrogen_calculation=False):
        """
        Remove atom.

        For Thiele forms of molecule causes invalidation of internal state.
        Implicit hydrogens marks will not be set if atoms in aromatic rings.
        Call `kekule()` and `thiele()` in sequence to fix marks.
        """
        ngb = self._bonds.pop(n)
        fix = self._atoms.pop(n).atomic_number != 1 and ngb and not _skip_hydrogen_calculation

        del self._charges[n]
        del self._radicals[n]
        del self._hydrogens[n]
        del self._plane[n]

        for m in ngb:
            del self._bonds[m][n]
            if not _skip_hydrogen_calculation:
                self._calc_implicit(m)

        self._conformers.clear()  # clean conformers. need full recalculation for new system
        try:
            del self._parsed_mapping[n]
        except KeyError:
            pass

        if fix:  # hydrogen atom not used for stereo coding
            self.fix_stereo()
        self.flush_cache()

    def delete_bond(self, n: int, m: int, *, _skip_hydrogen_calculation=False):
        """
        Disconnect atoms.

        For Thiele forms of molecule causes invalidation of internal state.
        Implicit hydrogens marks will not be set if atoms in aromatic rings.
        Call `kekule()` and `thiele()` in sequence to fix marks.
        """
        del self._bonds[n][m]
        del self._bonds[m][n]
        self._conformers.clear()  # clean conformers. need full recalculation for new system

        if not _skip_hydrogen_calculation:
            self._calc_implicit(n)
            self._calc_implicit(m)

        if self._atoms[n].atomic_number != 1 and self._atoms[m].atomic_number != 1 and not _skip_hydrogen_calculation:
            self.fix_stereo()
        self.flush_cache()

    def remap(self, mapping: Dict[int, int], *, copy: bool = False) -> 'MoleculeContainer':
        atoms = self._atoms  # keep original atoms dict
        h = super().remap(mapping, copy=copy)

        mg = mapping.get
        sp = self._plane
        shg = self._hydrogens

        if copy:
            h._MoleculeContainer__name = self.__name
            if self.__meta is not None:
                h._MoleculeContainer__meta = self.__meta.copy()
            hb = h._bonds
            hp = h._plane
            hhg = h._hydrogens
            hcf = h._conformers
            hm = h._parsed_mapping

            # deep copy of bonds
            for n, m_bond in self._bonds.items():
                n = mg(n, n)
                hb[n] = hbn = {}
                for m, bond in m_bond.items():
                    m = mg(m, m)
                    if m in hb:  # bond partially exists. need back-connection.
                        hbn[m] = hb[m][n]
                    else:
                        hbn[m] = bond = bond.copy()
                        bond._attach_graph(h, n, m)
        else:
            hb = {}
            hp = {}
            hhg = {}
            hcf = []
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
                        bond._change_map(n, m)

        for n in atoms:
            m = mg(n, n)
            hp[m] = sp[n]
            hhg[m] = shg[n]

        hcf.extend({mg(n, n): x for n, x in c.items()} for c in self._conformers)
        for n, m in self._parsed_mapping.items():
            hm[mg(n, n)] = m

        if copy:
            return h

        self._bonds = hb
        self._plane = hp
        self._hydrogens = hhg
        self._conformers = hcf
        self._parsed_mapping = hm
        return self

    def copy(self) -> 'MoleculeContainer':
        copy = super().copy()

        copy._bonds = cb = {}
        for n, m_bond in self._bonds.items():
            cb[n] = cbn = {}
            for m, bond in m_bond.items():
                if m in cb:  # bond partially exists. need back-connection.
                    cbn[m] = cb[m][n]
                else:
                    cbn[m] = bond = bond.copy()
                    bond._attach_graph(copy, n, m)

        copy._MoleculeContainer__name = self.__name
        if self.__meta is None:
            copy._MoleculeContainer__meta = None
        else:
            copy._MoleculeContainer__meta = self.__meta.copy()
        copy._plane = self._plane.copy()
        copy._hydrogens = self._hydrogens.copy()
        copy._parsed_mapping = self._parsed_mapping.copy()
        copy._conformers = [c.copy() for c in self._conformers]
        copy._atoms_stereo = self._atoms_stereo.copy()
        copy._allenes_stereo = self._allenes_stereo.copy()
        copy._cis_trans_stereo = self._cis_trans_stereo.copy()
        return copy

    def union(self, other: 'MoleculeContainer', *, remap: bool = False, copy: bool = True) -> 'MoleculeContainer':
        if not isinstance(other, MoleculeContainer):
            raise TypeError('MoleculeContainer expected')
        u, o = super().union(other, remap=remap, copy=copy)

        ub = u._bonds
        for n, m_bond in o._bonds.items():
            ub[n] = ubn = {}
            for m, bond in m_bond.items():
                if m in ub:  # bond partially exists. need back-connection.
                    ubn[m] = ub[m][n]
                else:
                    ubn[m] = bond = bond.copy()
                    bond._attach_graph(u, n, m)

        u._MoleculeContainer__name = u._MoleculeContainer__meta = None
        u._conformers.clear()
        u._plane.update(o._plane)
        u._hydrogens.update(o._hydrogens)
        u._parsed_mapping.update(o._parsed_mapping)
        return u

    def substructure(self, atoms: Iterable[int], *, as_query: bool = False, recalculate_hydrogens=True,
                     skip_neighbors_marks=False, skip_hybridizations_marks=False, skip_hydrogens_marks=False,
                     skip_rings_sizes_marks=False, skip_heteroatoms_marks=False) -> \
            Union['MoleculeContainer', 'QueryContainer']:
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
        :param skip_heteroatoms_marks: Don't set heteroatoms count marks
        """
        if not atoms:
            raise ValueError('empty atoms list not allowed')
        if set(atoms) - self._atoms.keys():
            raise ValueError('invalid atom numbers')
        atoms = tuple(n for n in self._atoms if n in atoms)  # save original order
        if as_query:
            atom_type = QueryElement
            bond_type = QueryBond
            sub = object.__new__(QueryContainer)
        else:
            atom_type = Element
            bond_type = Bond
            sub = object.__new__(self.__class__)
            sub._MoleculeContainer__name = sub._MoleculeContainer__meta = None

        sa = self._atoms
        sb = self._bonds
        sc = self._charges
        sr = self._radicals

        sub._charges = {n: sc[n] for n in atoms}
        sub._radicals = {n: sr[n] for n in atoms}

        sub._atoms = ca = {}
        for n in atoms:
            ca[n] = atom = atom_type.from_atom(sa[n])
            atom._attach_graph(sub, n)

        sub._bonds = cb = {}
        for n in atoms:
            cb[n] = cbn = {}
            for m, bond in sb[n].items():
                if m in cb:  # bond partially exists. need back-connection.
                    cbn[m] = cb[m][n]
                elif m in atoms:
                    cbn[m] = bond = bond_type.from_bond(bond)
                    if not as_query:
                        bond._attach_graph(sub, n, m)

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

            sub._masked = {n: False for n in atoms}
            if skip_heteroatoms_marks:
                sub._heteroatoms = {n: () for n in atoms}
            else:
                sha = self.heteroatoms
                sub._heteroatoms = {n: (sha(n),) for n in atoms}

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

            sp = self._plane
            sub._plane = {n: sp[n] for n in atoms}
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

    def compose(self, other: 'MoleculeContainer') -> 'CGRContainer':
        """
        Compose 2 graphs to CGR.
        """
        if not isinstance(other, MoleculeContainer):
            raise TypeError('MoleculeContainer expected')
        sa = self._atoms
        sc = self._charges
        sr = self._radicals
        sb = self._bonds

        bonds = []
        adj = defaultdict(lambda: defaultdict(lambda: [None, None]))

        oa = other._atoms
        oc = other._charges
        or_ = other._radicals
        ob = other._bonds

        common = sa.keys() & oa.keys()

        h = CGRContainer()
        ha = h._atoms
        hb = h._bonds
        hc = h._charges
        hpc = h._p_charges
        hr = h._radicals
        hpr = h._p_radicals

        for n in sa.keys() - common:  # cleavage atoms
            hc[n] = hpc[n] = sc[n]
            hr[n] = hpr[n] = sr[n]
            hb[n] = {}
            ha[n] = a = DynamicElement.from_atom(sa[n])
            a._attach_graph(h, n)

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
            hb[n] = {}
            ha[n] = a = DynamicElement.from_atom(oa[n])
            a._attach_graph(h, n)

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
                raise MappingError(f'atoms with number {n} not equal')

            hc[n] = sc[n]
            hpc[n] = oc[n]
            hr[n] = sr[n]
            hpr[n] = or_[n]
            hb[n] = {}
            ha[n] = a = DynamicElement.from_atom(san)
            a._attach_graph(h, n)

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

    def pack(self, *, compressed=True, check=True, version=2, order: List[int] = None) -> bytes:
        """
        Pack into compressed bytes.

        Note:

        * Less than 4096 atoms supported. Atoms mapping should be in range 1-4095.
        * Implicit hydrogens count should be in range 0-6 or unspecified.
        * Isotope shift should be in range -15 - 15 relatively chython.files._mdl.mol.common_isotopes
        * Atoms neighbors should be in range 0-15

        Format V2 specification::

            Big endian bytes order
            8 bit - 0x02 (format specification version)
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
            3 bit - hydrogens (0-7). Note: 7 == None
            4 bit - charge (charge + 4. possible range -4 - 4)
            1 bit - radical state
            Connection table: flatten list of neighbors. neighbors count stored in atom block.
            For example CC(=O)O - {1: [2], 2: [1, 3, 4], 3: [2], 4: [2]} >> [2, 1, 3, 4, 2, 2].
            Repeated block (equal to bonds count).
            24 bit - paired 12 bit numbers.
            Bonds order block 3 bit per bond zero-padded to full byte at the end.
            Cis/trans data block (repeated):
            24 bit - atoms pair
            7 bit - zero padding. in future can be used for extra bond-level stereo, like atropoisomers.
            1 bit - sign

        Format V3 specification::

            Big endian bytes order
            8 bit - 0x03 (format specification version)
            Atom block 3 bytes (repeated):
            1 bit - atom entrance flag (always 1)
            7 bit - atomic number (<=118)
            3 bit - hydrogens (0-7). Note: 7 == None
            4 bit - charge (charge + 4. possible range -4 - 4)
            1 bit - radical state
            1 bit padding
            3 bit tetrahedron/allene sign
                (000 - not stereo or unknown, 001 - pure-unknown-enantiomer, 010 or 011 - has stereo)
            4 bit - number of following bonds and CT blocks (0-15)

            Bond block 2 bytes (repeated 0-15 times)
            12 bit - negative shift from current atom to connected (e.g. 0x001 = -1 - connected to previous atom)
            4 bit - bond order: 0000 - single, 0001 - double, 0010 - triple, 0011 - aromatic, 0111 - special

            Cis-Trans 2 bytes
            12 bit - negative shift from current atom to connected (e.g. 0x001 = -1 - connected to previous atom)
            4 bit - CT sign: 1000 or 1001 - to avoid overlap with bond

        V2 format is faster than V3. V3 format doesn't include isotopes, atom numbers and XY coordinates.

        :param compressed: return zlib-compressed pack.
        :param check: check molecule for format restrictions.
        :param version: format version
        :param order: atom order in V3
        """
        from ._pack import pack

        if check:
            bonds = self._bonds
            if not bonds:
                raise ValueError('Empty molecules not supported')
            if max(bonds) > 4095:
                raise ValueError('Big molecules not supported')
            if any(len(x) > 15 for x in bonds.values()):
                raise ValueError('To many neighbors not supported')

        if version == 2:
            data = pack(self)
        elif version == 3:
            data = self._cpack(order, check)
        else:
            raise ValueError('invalid specification version')
        if compressed:
            return compress(data, 9)
        return data

    @classmethod
    def pack_len(cls, data: bytes, /, *, compressed=True) -> int:
        """
        Returns atoms count in molecule pack.
        """
        if compressed:
            data = decompress(data)
        if data[0] not in (0, 2):
            raise ValueError('invalid pack header')
        return int.from_bytes(data[1:3], 'big') >> 4

    @classmethod
    def unpack(cls, data: Union[bytes, memoryview], /, *, compressed=True,
               _return_pack_length=False) -> 'MoleculeContainer':
        """
        Unpack from compressed bytes.

        :param compressed: decompress data before processing.
        """
        from ._unpack import unpack
        from ._cpack import unpack as cpack

        if compressed:
            data = decompress(data)
        if data[0] in (0, 2):
            (mapping, atom_numbers, isotopes, charges, radicals, hydrogens, plane, bonds,
             atoms_stereo, allenes_stereo, cis_trans_stereo, pack_length, bonds_flat) = unpack(data)
        elif data[0] == 3:
            (mapping, atom_numbers, isotopes, charges, radicals, hydrogens, plane, bonds,
             atoms_stereo, allenes_stereo, cis_trans_stereo, pack_length, bonds_flat) = cpack(data)
        else:
            raise ValueError('invalid pack header')

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
            a._n = n
        for b in bonds_flat:
            b._Bond__graph = ref(mol)

        if _return_pack_length:
            return mol, pack_length
        return mol

    def _cpack(self, order=None, check=True):
        if order is None:
            order = list(self._atoms)
        elif check:
            if not isinstance(order, (list, tuple)):
                raise TypeError('invalid atoms order')
            elif len(so := set(order)) != len(order) or not so.issubset(self._atoms):
                raise ValueError('invalid atoms order')

        atoms = self._atoms
        bonds = self._bonds
        charges = self._charges
        radicals = self._radicals
        hydrogens = self._hydrogens
        atoms_stereo = self._atoms_stereo
        allenes_stereo = self._allenes_stereo
        allenes_terminals = self._stereo_allenes_terminals

        cumulenes = {}
        ct_map = {}
        for n, m in self._cis_trans_stereo:
            ct_map[n] = m
            ct_map[m] = n
            cumulenes[n] = [x for x, b in bonds[n].items() if b.order in (1, 4)]
            cumulenes[m] = [x for x, b in bonds[m].items() if b.order in (1, 4)]

        for c in self._allenes_stereo:
            n, m = allenes_terminals[c]
            cumulenes[n] = [x for x, b in bonds[n].items() if b.order in (1, 4)]
            cumulenes[m] = [x for x, b in bonds[m].items() if b.order in (1, 4)]

        seen = {}
        data = [b'\x03']
        for i, n in enumerate(order):
            seen[n] = i
            env = bonds[n]

            data.append((0x80 | atoms[n].atomic_number).to_bytes(1, 'big'))

            # 3 bit - hydrogens (0-6, None) | 4 bit - charge | 1 bit - radical
            hcr = (charges[n] + 4) << 1 | radicals[n]
            if (h := hydrogens[n]) is None:
                hcr |= 0b11100000
            else:
                hcr |= h << 5
            data.append(hcr.to_bytes(1, 'big'))

            if n in atoms_stereo:
                if self._translate_tetrahedron_sign(n, [x for x in order if x in env]):
                    s = 0b0011_0000
                else:
                    s = 0b0010_0000
            elif n in allenes_stereo:
                t1, t2 = allenes_terminals[n]
                nn = None
                for x in order:
                    if nn is None:
                        if x in cumulenes[t1]:
                            nn = x
                            flag = True
                        elif x in cumulenes[t2]:
                            flag = False
                            nn = x
                    elif flag:  # noqa
                        if x in cumulenes[t2]:
                            nm = x
                            break
                    elif x in cumulenes[t1]:
                        nm = x
                        break
                if self._translate_allene_sign(n, nn, nm):  # noqa
                    s = 0b0011_0000
                else:
                    s = 0b0010_0000
            else:
                s = 0

            tmp = []
            for m in order[:i]:
                if (b := env.get(m)) is not None:
                    tmp.append(((i - seen[m]) << 4 | b.order - 1).to_bytes(2, 'big'))
            if n in ct_map and (m := ct_map[n]) in seen:  # only right atom codes stereo sign
                nm = None
                for x in order:
                    if nm is None:
                        if x in cumulenes[n]:
                            nm = x
                            flag = True
                        elif x in cumulenes[m]:
                            nm = x
                            flag = False
                    elif flag:  # noqa
                        if x in cumulenes[m]:
                            nn = x
                            break
                    elif x in cumulenes[n]:
                        nn = x
                        break
                if self._translate_cis_trans_sign(m, n, nm, nn):  # noqa
                    cs = 0b1001
                else:
                    cs = 0b1000
                tmp.append(((i - seen[m]) << 4 | cs).to_bytes(2, 'big'))

            data.append((s | len(tmp)).to_bytes(1, 'big'))
            data.extend(tmp)
        return b''.join(data)

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

    def _calc_implicit(self, n: int):
        """
        Set firs possible hydrogens count based on rules
        """
        atoms = self._atoms
        atom = atoms[n]
        if (an := atom.atomic_number) == 1:  # hydrogen nether has implicit H
            self._hydrogens[n] = 0
            return

        charge: int = self._charges[n]
        is_radical = self._radicals[n]
        explicit_sum = 0
        explicit_dict = defaultdict(int)
        aroma = 0
        for m, bond in self._bonds[n].items():
            order = bond.order
            if order == 4:  # only neutral carbon aromatic rings supported
                if not charge and not is_radical and an == 6:
                    aroma += 1
                else:  # use `kekule()` to calculate proper implicit hydrogens count
                    self._hydrogens[n] = None
                    return
            elif order != 8:  # any bond used for complexes
                explicit_sum += order
                explicit_dict[(order, atoms[m].atomic_number)] += 1

        if aroma == 2:
            if explicit_sum == 0:  # H-Ar
                self._hydrogens[n] = 1
            elif explicit_sum == 1:  # R-Ar
                self._hydrogens[n] = 0
            else:  # invalid aromaticity
                self._hydrogens[n] = None
            return
        elif aroma == 3:  # condensed rings
            if explicit_sum:  # invalid aromaticity
                self._hydrogens[n] = None
            else:
                self._hydrogens[n] = 0
            return
        elif aroma:
            self._hydrogens[n] = None
            return

        try:
            rules = atom.valence_rules(charge, is_radical, explicit_sum)
        except ValenceError:
            self._hydrogens[n] = None
            return
        for s, d, h in rules:
            if s.issubset(explicit_dict) and all(explicit_dict[k] >= c for k, c in d.items()):
                self._hydrogens[n] = h
                return
        self._hydrogens[n] = None  # rule not found

    def _check_implicit(self, n: int, h: int) -> bool:
        atoms = self._atoms
        atom = atoms[n]
        if atom.atomic_number == 1:  # hydrogen nether has implicit H
            return h == 0

        explicit_sum = 0
        explicit_dict = defaultdict(int)

        for m, bond in self._bonds[n].items():
            order = bond.order
            if order == 4:  # can't check aromatic rings
                return False
            elif order != 8:  # any bond used for complexes
                explicit_sum += order
                explicit_dict[(order, atoms[m].atomic_number)] += 1

        try:
            rules = atom.valence_rules(self._charges[n], self._radicals[n], explicit_sum)
        except ValenceError:
            return False
        for s, d, _h in rules:
            if h == _h and s.issubset(explicit_dict) and all(explicit_dict[k] >= c for k, c in d.items()):
                return True
        return False

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
            atom._attach_graph(self, n)

        bonds = {}
        for n, m_bond in self._bonds.items():
            bonds[n] = cbn = {}
            for m, bond in m_bond.items():
                if m in bonds:  # bond partially exists. need back-connection.
                    cbn[m] = bonds[m][n]
                else:
                    cbn[m] = bond = bond.copy()
                    bond._attach_graph(self, n, m)

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
                'plane': self._plane, **super().__getstate__()}

    def __setstate__(self, state):
        super().__setstate__(state)
        self._conformers = state['conformers']
        self._atoms_stereo = state['atoms_stereo']
        self._allenes_stereo = state['allenes_stereo']
        self._cis_trans_stereo = state['cis_trans_stereo']
        self._hydrogens = state['hydrogens']
        self._parsed_mapping = state['parsed_mapping']
        self._plane = state['plane']
        self.__meta = state['meta']
        self.__name = state['name']

        # attach bonds to graph
        for n, m, b in self.bonds():
            b._attach_graph(self, n, m)


__all__ = ['MoleculeContainer']
