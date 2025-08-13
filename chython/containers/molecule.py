# -*- coding: utf-8 -*-
#
#  Copyright 2017-2025 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from typing import Dict, Iterable, List, Tuple, Union
from zlib import compress, decompress
from .bonds import Bond, DynamicBond
from .cgr import CGRContainer
from .graph import Graph
from .rdkit import RDkit
from ..algorithms.aromatics import Aromatize
from ..algorithms.calculate2d import Calculate2DMolecule
from ..algorithms.conformers import Conformers
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
from ..periodictable import DynamicElement, Element, H as _H


# atomic number constants
H = 1
C = 6


class MoleculeContainer(MoleculeStereo, Graph[Element, Bond], Morgan, Rings, MoleculeIsomorphism,
                        Aromatize, StandardizeMolecule, MoleculeSmiles, DepictMolecule, Calculate2DMolecule,
                        Conformers, Fingerprints, Tautomers, RDkit, MCS, X3domMolecule):
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
        return sum(a.charge for _, a in self.atoms())

    @cached_property
    def is_radical(self) -> bool:
        """
        True if at least one atom is radical
        """
        return any(a.is_radical for _, a in self.atoms())

    @cached_property
    def molecular_mass(self) -> float:
        h = _H().atomic_mass
        return sum(a.atomic_mass + (a.implicit_hydrogens or 0) * h for _, a in self.atoms())

    @cached_property
    def brutto(self) -> Dict[str, int]:
        """Counted atoms dict"""
        c = Counter()
        # make an order
        c['C'] = 0
        c['H'] = 0
        c['O'] = 0
        c['N'] = 0
        c['B'] = 0
        c.update(a.atomic_symbol for a in sorted((a for _, a in self.atoms()), key=lambda a: a.atomic_number))
        c['H'] += sum(a.implicit_hydrogens or 0 for _, a in self.atoms())
        return {k: v for k, v in c.items() if v}

    @cached_property
    def brutto_formula(self) -> str:
        return ''.join(f'{a}{c}' if c > 1 else a for a, c in self.brutto.items())

    @cached_property
    def brutto_formula_html(self) -> str:
        return ''.join(f'{a}<sub>{c}</sub>' if c > 1 else a for a, c in self.brutto.items())

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
        if bond == 8:
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
            if bond == 8:
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
        if self._bonds[m].pop(n) != 8:
            if self._changed is None:
                self._changed = {n, m}
            else:
                self._changed.add(n)
                self._changed.add(m)
        if not _skip_calculation and self._backup is None:
            self.fix_structure()
            self.fix_stereo()

    def copy(self, *, keep_sssr=False, keep_components=False) -> 'MoleculeContainer':
        copy = super().copy()
        copy._name = self._name
        if self._meta is None:
            copy._meta = None
        else:
            copy._meta = self._meta.copy()

        if keep_sssr:
            for k,  v in self.__dict__.items():
                if k in ('sssr', 'atoms_rings', 'atoms_rings_sizes', 'not_special_connectivity', 'rings_count'):
                    copy.__dict__[k] = v
        if keep_components:
            if 'connected_components' in self.__dict__:
                copy.__dict__['connected_components'] = self.connected_components
        return copy

    def union(self, other: 'MoleculeContainer', *, remap: bool = False, copy: bool = True) -> 'MoleculeContainer':
        if not isinstance(other, MoleculeContainer):
            raise TypeError('MoleculeContainer expected')
        return super().union(other, remap=remap, copy=copy)

    def substructure(self, atoms: Iterable[int], *, recalculate_hydrogens=True) -> 'MoleculeContainer':
        """
        Create substructure containing atoms from atoms list.

        For Thiele forms of molecule In Molecule substructure causes invalidation of internal state.
        Implicit hydrogens marks will not be set if atoms in aromatic rings.
        Call `kekule()` and `thiele()` in sequence to fix marks.

        :param atoms: list of atoms numbers of substructure
        :param recalculate_hydrogens: calculate implicit H count in substructure
        """
        if not atoms:
            raise ValueError('empty atoms list not allowed')
        if set(atoms) - self._atoms.keys():
            raise ValueError('invalid atom numbers')
        atoms = tuple(n for n in self if n in atoms)  # save original order
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

    def compose(self, other: 'MoleculeContainer', dynamic=True) -> Union['CGRContainer', 'MoleculeContainer']:
        """
        Compose 2 graphs to CGR.

        :param dynamic: produce CGR with dynamic bonds and atoms,
            overwise keep reactants' electronic and implicit hydrogens state and label dynamic bonds as "any".
            This representation can't catch atom-only changes like (de)protonation, red-ox, etc;
            or ambiguous bond changes like triple to double bond reduction.
        """
        if not isinstance(other, MoleculeContainer):
            raise TypeError('MoleculeContainer expected')
        bonds = []
        adj = defaultdict(lambda: defaultdict(lambda: [None, None]))
        common = self._atoms.keys() & other._atoms.keys()

        if dynamic:
            h = CGRContainer()
            from_atom = DynamicElement.from_atom
            from_atoms = DynamicElement.from_atoms
            from_bond = DynamicBond.from_bond
            dynamic_bond = DynamicBond
        else:
            h = self.__class__()
            from_atom = from_atoms = lambda x, *_: x.copy(hydrogens=True)
            from_bond = lambda x: x.copy()
            dynamic_bond = lambda x, y: Bond(8 if x != y else x)

        ha = h._atoms
        hb = h._bonds

        for n in self._atoms.keys() - common:  # cleavage atoms
            ha[n] = from_atom(self._atoms[n])
            hb[n] = {}
            for m, bond in self._bonds[n].items():
                if m not in ha:
                    if m in common:  # bond to common atoms is broken bond
                        bond = dynamic_bond(bond.order, None)
                    else:
                        bond = from_bond(bond)
                    bonds.append((n, m, bond))
        for n in other._atoms.keys() - common:  # coupling atoms
            ha[n] = from_atom(other._atoms[n])
            hb[n] = {}

            for m, bond in other._bonds[n].items():
                if m not in ha:
                    if m in common:  # bond to common atoms is formed bond
                        bond = dynamic_bond(None, bond.order)
                    else:
                        bond = from_bond(bond)
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
            ha[n] = from_atoms(self._atoms[n], other._atoms[n])
            hb[n] = {}

            for m, (o1, o2) in adj[n].items():
                if m not in ha:
                    bonds.append((n, m, dynamic_bond(o1, o2)))

        for n, m, bond in bonds:
            hb[n][m] = hb[m][n] = bond

        if not dynamic:
            h.calc_labels()
        return h

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

        :param compressed: return zlib-compressed pack.
        :param check: check molecule for format restrictions.
        :param version: format version. Only V2 is supported.
        :param order: atom order in V3
        """
        from ._pack_v2 import pack as pack_v2

        if check:
            bonds = self._bonds
            if not bonds:
                raise ValueError('Empty molecules not supported')
            if max(bonds) > 4095:
                raise ValueError('Big molecules not supported')
            if any(len(x) > 15 for x in bonds.values()):
                raise ValueError('To many neighbors not supported')

        if version == 2:
            data = pack_v2(self)
        else:
            raise ValueError('invalid specification version')
        if compressed:
            return compress(data, 9)
        return data

    def pach(self, *, compressed=True, check=True, version=2, order: List[int] = None) -> bytes:
        return self.pack(compressed=compressed, check=check, version=version, order=order)

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
    def unpack(cls, data: Union[bytes, memoryview], /, *, compressed=True, skip_labels_calculation=False,
               _return_pack_length=False) -> 'MoleculeContainer':
        """
        Unpack from compressed bytes.

        :param compressed: decompress data before processing.
        """
        from ._unpack_v0v2 import unpack as unpack_v0v2

        if compressed:
            data = decompress(data)
        if data[0] in (0, 2):
            mol, cis_trans, pack_length = unpack_v0v2(data)
            for n, m, s in cis_trans:
                if n in mol._stereo_cis_trans_centers:  # check for invalid CT data
                    mol.bond(*mol._stereo_cis_trans_centers[n])._stereo = s
        else:
            raise ValueError('invalid pack header')

        if not skip_labels_calculation:
            mol.calc_labels()

        if _return_pack_length:
            return mol, pack_length
        return mol

    @classmethod
    def unpach(cls, data: Union[bytes, memoryview], /, *, compressed=True) -> 'MoleculeContainer':
        """
        Unpack from compressed bytes.
        """
        return cls.unpack(data, compressed=compressed)

    def __bytes__(self):
        return self.pack()

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

                if bond == 8:
                    continue
                elif bond == 4:
                    hybridization = 4
                elif hybridization != 4:
                    if bond == 3:
                        hybridization = 3
                    elif bond == 2:
                        if hybridization == 1:
                            hybridization = 2
                        elif hybridization == 2:
                            hybridization = 3

                neighbors += 1
                if (a := atoms[m]) == H:
                    explicit_hydrogens += 1
                elif a != C:
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
        if (atom := self._atoms[n]) == H:  # hydrogen nether has implicit H
            atom._implicit_hydrogens = 0
            return

        explicit_sum = 0
        explicit_dict = defaultdict(int)
        aroma = 0
        for m, bond in self._bonds[n].items():
            if bond == 4:  # only neutral carbon aromatic rings supported
                if not atom.charge and not atom.is_radical and atom == C:
                    aroma += 1
                else:  # use `kekule()` to calculate proper implicit hydrogens count
                    atom._implicit_hydrogens = None
                    return
            elif bond != 8:  # any bond used for complexes
                explicit_sum += bond.order
                explicit_dict[(bond.order, self._atoms[m].atomic_number)] += 1

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
        if (atom := self._atoms[n]) == H:  # hydrogen nether has implicit H
            return h == 0

        explicit_sum = 0
        explicit_dict = defaultdict(int)

        for m, bond in self._bonds[n].items():
            if bond == 4:  # can't check aromatic rings
                return False
            elif bond != 8:  # any bond used for complexes
                explicit_sum += bond.order
                explicit_dict[(bond.order, self._atoms[m].atomic_number)] += 1

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
            for k,  v in self.__dict__.items():
                if k in ('sssr', 'atoms_rings', 'atoms_rings_sizes', 'not_special_connectivity', 'rings_count'):
                    backup[k] = v
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
        self._backup = self.copy(keep_sssr=True, keep_components=True)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if exc_type:  # restore state
            backup = self._backup
            self._atoms = backup._atoms
            self._bonds = backup._bonds
            self._meta = backup._meta
            self._name = backup._name
            self.__dict__ = backup.__dict__
        else:  # update internal state
            self.fix_structure()
            self.fix_stereo()
        self._backup = None  # drop backup


__all__ = ['MoleculeContainer']
