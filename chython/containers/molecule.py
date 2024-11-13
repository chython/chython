# -*- coding: utf-8 -*-
#
#  Copyright 2017-2024 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from ..algorithms.morgan import Morgan
from ..algorithms.rings import Rings
from ..algorithms.smiles import MoleculeSmiles
from ..algorithms.standardize import StandardizeMolecule
from ..algorithms.stereo import MoleculeStereo
from ..algorithms.tautomers import Tautomers
from ..algorithms.x3dom import X3domMolecule
from ..exceptions import ValenceError
from ..periodictable import DynamicElement, Element, QueryElement, H


class MoleculeContainer(MoleculeStereo, Graph[Element, Bond], Morgan, Rings, MoleculeIsomorphism,
                        Aromatize, StandardizeMolecule, MoleculeSmiles, DepictMolecule, Calculate2DMolecule,
                        Fingerprints, Tautomers, MCS, X3domMolecule):
    __slots__ = ('_meta', '_name', '_conformers', '_changed', '_backup')

    def __init__(self):
        super().__init__()
        self._meta = None
        self._name = None
        self._changed = None
        self._backup = None

    @property
    def meta(self) -> Dict:
        if self._meta is None:
            self._meta = {}  # lazy
        return self._meta

    @property
    def name(self) -> str:
        return self._name or ''

    @name.setter
    def name(self, name):
        if not isinstance(name, str):
            raise TypeError('name should be a string preferably up to 80 symbols')
        self._name = name

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
    def adjacency_matrix(self, set_bonds=False, /):
        """
        Adjacency matrix of Graph.

        :param set_bonds: if True set bond orders instead of 1.
        """
        from numpy import uint, zeros

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
        return sum(a.charge for a in self._atoms.values())

    @cached_property
    def is_radical(self) -> bool:
        """
        True if at least one atom is radical
        """
        return any(a.is_radical for a in self._atoms.values())

    @cached_property
    def molecular_mass(self) -> float:
        h = H().atomic_mass
        return sum(a.atomic_mass + a.implicit_hydrogens * h for a in self._atoms.values())

    @cached_property
    def brutto(self) -> Dict[str, int]:
        """Counted atoms dict"""
        c = Counter(a.atomic_symbol for a in self._atoms.values())
        c['H'] += sum(a.implicit_hydrogens for a in self._atoms.values())
        return dict(c)

    @cached_property
    def aromatic_rings(self) -> Tuple[Tuple[int, ...], ...]:
        """
        Aromatic rings atoms numbers
        """
        bonds = self._bonds
        return tuple(ring for ring in self.sssr if bonds[ring[0]][ring[-1]] == 4
                     and all(bonds[n][m] == 4 for n, m in zip(ring, ring[1:])))

    def add_atom(self, atom: Union[Element, int, str], *args, _skip_calculation=False, **kwargs):
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

        n = super().add_atom(atom, *args, **kwargs)
        if self._changed is None:
            self._changed = {n}
        else:
            self._changed.add(n)
        if not _skip_calculation and self._backup is None:
            self.fix_structure()
        return n

    def add_bond(self, n, m, bond: Union[Bond, int], *, _skip_calculation=False):
        """
        Connect atoms with bonds.

        For Thiele forms of molecule causes invalidation of internal state.
        Implicit hydrogens marks will not be set if atoms in aromatic rings.
        Call `kekule()` and `thiele()` in sequence to fix marks.
        """
        if not isinstance(bond, Bond):
            bond = Bond(bond)

        super().add_bond(n, m, bond)
        if bond.order == 8:
            return  # any bond doesn't change anything
        if self._changed is None:
            self._changed = {n, m}
        else:
            self._changed.add(n)
            self._changed.add(m)
        if not _skip_calculation and self._backup is None:
            self.fix_structure()
            self.fix_stereo()

    def delete_atom(self, n: int, *, _skip_calculation=False):
        """
        Remove atom.

        For Thiele forms of molecule causes invalidation of internal state.
        Implicit hydrogens marks will not be set if atoms in aromatic rings.
        Call `kekule()` and `thiele()` in sequence to fix marks.
        """
        del self._atoms[n]
        for m, bond in self._bonds.pop(n).items():
            del self._bonds[m][n]
            if bond.order == 8:
                continue
            if self._changed is None:
                self._changed = {m}
            else:
                self._changed.add(m)
        if not _skip_calculation and self._backup is None:
            self.fix_structure()
            self.fix_stereo()

    def delete_bond(self, n: int, m: int, *, _skip_calculation=False):
        """
        Disconnect atoms.

        For Thiele forms of molecule causes invalidation of internal state.
        Implicit hydrogens marks will not be set if atoms in aromatic rings.
        Call `kekule()` and `thiele()` in sequence to fix marks.
        """
        del self._bonds[n][m]
        if self._bonds[m].pop(n).order != 8:
            if self._changed is None:
                self._changed = {n, m}
            else:
                self._changed.add(n)
                self._changed.add(m)
        if not _skip_calculation and self._backup is None:
            self.fix_structure()
            self.fix_stereo()

    def copy(self) -> 'MoleculeContainer':
        copy = super().copy()
        copy._name = self._name
        if self._meta is None:
            copy._meta = None
        else:
            copy._meta = self._meta.copy()
        return copy

    def union(self, other: 'MoleculeContainer', *, remap: bool = False, copy: bool = True) -> 'MoleculeContainer':
        if not isinstance(other, MoleculeContainer):
            raise TypeError('MoleculeContainer expected')
        return super().union(other, remap=remap, copy=copy)

    def substructure(self, atoms: Iterable[int], *, as_query: bool = False, recalculate_hydrogens=True,
                     skip_neighbors_marks=False, skip_hybridizations_marks=False, skip_hydrogens_marks=False,
                     skip_rings_sizes_marks=False, skip_heteroatoms_marks=False, skip_in_ring_bond_marks=False,
                     skip_stereo_marks=False) -> \
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
        :param skip_in_ring_bond_marks: Don't set in_ring bond marks
        :param skip_stereo_marks: Don't set stereo marks on substructured queries
        """
        if not atoms:
            raise ValueError('empty atoms list not allowed')
        if set(atoms) - self._atoms.keys():
            raise ValueError('invalid atom numbers')
        atoms = tuple(n for n in self._atoms if n in atoms)  # save original order
        if as_query:
            sub = object.__new__(QueryContainer)

            lost = {n for n, a in self._atoms.items() if a.atomic_number != 1} - set(atoms)  # atoms not in substructure
            # atoms with fully present neighbors
            not_skin = {n for n in atoms if lost.isdisjoint(self._bonds[n])}

            # check for full presence of cumulene chains and terminal attachments
            for p in self.stereogenic_cumulenes.values():
                if not not_skin.issuperset(p):
                    not_skin.difference_update(p)

            sub._atoms = {n: QueryElement.from_atom(self._atoms[n],
                                                    neighbors=not skip_neighbors_marks,
                                                    hybridization=not skip_hybridizations_marks,
                                                    hydrogens=not skip_hydrogens_marks,
                                                    ring_sizes=not skip_rings_sizes_marks,
                                                    heteroatoms=not skip_heteroatoms_marks,
                                                    stereo=not skip_stereo_marks and n in not_skin)
                          for n in atoms}
            sub._bonds = sb = {}
            for n in atoms:
                sb[n] = sbn = {}
                for m, bond in self._bonds[n].items():
                    if m in sb:  # bond partially exists. need back-connection.
                        sbn[m] = sb[m][n]
                    elif m in atoms:
                        sbn[m] = QueryBond.from_bond(bond,
                                                     in_ring=not skip_in_ring_bond_marks,
                                                     stereo=not skip_stereo_marks and n in not_skin and m in not_skin)
            return sub

        # molecule substructure
        sub = object.__new__(self.__class__)
        sub._name = sub._meta = sub._changed = None
        sub._atoms = {n: self._atoms[n].copy(hydrogens=not recalculate_hydrogens, stereo=True) for n in atoms}
        sub._bonds = sb = {}
        for n in atoms:
            sb[n] = sbn = {}
            for m, bond in self._bonds[n].items():
                if m in sb:  # bond partially exists. need back-connection.
                    sbn[m] = sb[m][n]
                elif m in atoms:
                    sbn[m] = bond.copy(stereo=True)
        sub.fix_structure(recalculate_hydrogens=recalculate_hydrogens)
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
        bonds = []
        adj = defaultdict(lambda: defaultdict(lambda: [None, None]))
        common = self._atoms.keys() & other._atoms.keys()

        h = CGRContainer()
        ha = h._atoms
        hb = h._bonds

        for n in self._atoms.keys() - common:  # cleavage atoms
            ha[n] = DynamicElement.from_atom(self._atoms[n])
            hb[n] = {}
            for m, bond in self._bonds[n].items():
                if m not in ha:
                    if m in common:  # bond to common atoms is broken bond
                        bond = DynamicBond(bond.order, None)
                    else:
                        bond = DynamicBond.from_bond(bond)
                    bonds.append((n, m, bond))
        for n in other._atoms.keys() - common:  # coupling atoms
            ha[n] = DynamicElement.from_atom(other._atoms[n])
            hb[n] = {}

            for m, bond in other._bonds[n].items():
                if m not in ha:
                    if m in common:  # bond to common atoms is formed bond
                        bond = DynamicBond(None, bond.order)
                    else:
                        bond = DynamicBond.from_bond(bond)
                    bonds.append((n, m, bond))
        for n in common:
            an = adj[n]
            for m, bond in self._bonds[n].items():
                if m in common:
                    an[m][0] = bond.order
            for m, bond in other._bonds[n].items():
                if m in common:
                    an[m][1] = bond.order
        for n in common:
            ha[n] = DynamicElement.from_atoms(self._atoms[n], other._atoms[n])
            hb[n] = {}

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
        if atoms - bonds.keys():
            raise ValueError('invalid atom numbers')
        nodes = [atoms]
        for _ in range(deep):
            n = {y for x in nodes[-1] for y in bonds[x]} | nodes[-1]
            if n in nodes:
                break
            nodes.append(n)
        return nodes

    def fix_structure(self, recalculate_hydrogens=True):
        """
        Fix molecule internal representation
        """
        self.calc_labels()  # refresh all labels

        if recalculate_hydrogens:
            for n in (self._changed or self._atoms):
                self.calc_implicit(n)  # fix Hs count
        self._changed = None

    def calc_labels(self):
        atoms = self._atoms
        atoms_rings_sizes = self.atoms_rings_sizes  # expensive: sssr based
        atoms_rings = {n: set(r) for n, r in self.atoms_rings.items()}

        for n, m_bond in self._bonds.items():
            neighbors = 0
            heteroatoms = 0
            hybridization = 1
            explicit_hydrogens = 0
            anr = atoms_rings.get(n) or False
            for m, bond in m_bond.items():
                bond._in_ring = anr and (amr := atoms_rings.get(m) or False) and not anr.isdisjoint(amr)  # have common rings

                order = bond.order
                if order == 8:
                    continue
                elif order == 4:
                    hybridization = 4
                elif hybridization != 4:
                    if order == 3:
                        hybridization = 3
                    elif order == 2:
                        if hybridization == 1:
                            hybridization = 2
                        elif hybridization == 2:
                            hybridization = 3

                neighbors += 1
                an = atoms[m].atomic_number
                if an == 1:
                    explicit_hydrogens += 1
                elif an != 6:
                    heteroatoms += 1
            atom = atoms[n]
            atom._neighbors = neighbors
            atom._heteroatoms = heteroatoms
            atom._hybridization = hybridization
            atom._explicit_hydrogens = explicit_hydrogens

            atom._in_ring = n in atoms_rings_sizes
            atom._ring_sizes = atoms_rings_sizes.get(n) or set()

    def calc_implicit(self, n: int):
        """
        Set firs possible hydrogens count based on rules
        """
        atom = self._atoms[n]
        if atom.atomic_number == 1:  # hydrogen nether has implicit H
            atom._implicit_hydrogens = 0
            return

        explicit_sum = 0
        explicit_dict = defaultdict(int)
        aroma = 0
        for m, bond in self._bonds[n].items():
            order = bond.order
            if order == 4:  # only neutral carbon aromatic rings supported
                if not atom.charge and not atom.is_radical and atom.atomic_number == 6:
                    aroma += 1
                else:  # use `kekule()` to calculate proper implicit hydrogens count
                    atom._implicit_hydrogens = None
                    return
            elif order != 8:  # any bond used for complexes
                explicit_sum += order
                explicit_dict[(order, self._atoms[m].atomic_number)] += 1

        if aroma == 2:
            if explicit_sum == 0:  # H-Ar
                atom._implicit_hydrogens = 1
            elif explicit_sum == 1:  # R-Ar
                atom._implicit_hydrogens = 0
            else:  # invalid aromaticity
                atom._implicit_hydrogens = None
            return
        elif aroma == 3:  # condensed rings
            if explicit_sum:  # invalid aromaticity
                atom._implicit_hydrogens = None
            else:
                atom._implicit_hydrogens = 0
            return
        elif aroma:
            atom._implicit_hydrogens = None
            return

        try:
            rules = atom.valence_rules(explicit_sum)
        except ValenceError:
            atom._implicit_hydrogens = None
            return
        for s, d, h in rules:
            if s.issubset(explicit_dict) and all(explicit_dict[k] >= c for k, c in d.items()):
                atom._implicit_hydrogens = h
                return
        atom._implicit_hydrogens = None  # rule not found

    def check_implicit(self, n: int, h: int) -> bool:
        atom = self._atoms[n]
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
                explicit_dict[(order, self._atoms[m].atomic_number)] += 1

        try:
            rules = atom.valence_rules(explicit_sum)
        except ValenceError:
            return False
        for s, d, _h in rules:
            if h == _h and s.issubset(explicit_dict) and all(explicit_dict[k] >= c for k, c in d.items()):
                return True
        return False

    def flush_cache(self, *, keep_sssr=False, keep_components=False):
        backup = {}
        if keep_sssr:
            # good to keep if no new bonds or bonds deletions or bonds to/from any change
            if 'sssr' in self.__dict__:
                backup['sssr'] = self.sssr
            if 'atoms_rings' in self.__dict__:
                backup['atoms_rings'] = self.atoms_rings
            if 'atoms_rings_sizes' in self.__dict__:
                backup['atoms_rings_sizes'] = self.atoms_rings_sizes
            if 'ring_atoms' in self.__dict__:
                backup['ring_atoms'] = self.ring_atoms
            if 'not_special_connectivity' in self.__dict__:
                backup['not_special_connectivity'] = self.not_special_connectivity
            if 'rings_count' in self.__dict__:
                backup['rings_count'] = self.rings_count
        if keep_components:
            # good to keep if no new bonds or bonds deletions
            if 'connected_components' in self.__dict__:
                backup['connected_components'] = self.connected_components
        self.__dict__ = backup

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
        self._backup = self.copy()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if exc_type:  # restore state
            backup = self._backup
            self._atoms = backup._atoms
            self._bonds = backup._bonds
            self._meta = backup._meta
            self._name = backup._name
            self.flush_cache()
        else:  # update internal state
            self.fix_structure()
            self.fix_stereo()
        self._backup = None  # drop backup


__all__ = ['MoleculeContainer']
