# -*- coding: utf-8 -*-
#
#  Copyright 2024 Timur Gimadiev <timur.gimadiev@gmail.com>
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
from collections import defaultdict
from itertools import product
from typing import Dict, Iterable, Optional

from .molecule import Bond, MoleculeContainer

var_atoms = ["X", "Y", "Z"]
var_groups = ["R"]


class MarkushContainer:
    __slots__ = (
        "__meta",
        "__tree",
        "__r_groups",
        "__r_groups_map",
        "__substituents",
        "__initial_mol",
    )

    def __init__(self, substituents: Optional[list["MoleculeContainer"]] = None):
        super().__init__()
        self.__tree = {}
        self.__meta = None
        self.__r_groups = {}
        self.__substituents = substituents
        self.__r_groups_map = {}
        self.__initial_mol = None

    @staticmethod
    def connect_w_bond(first: MoleculeContainer, other: MoleculeContainer, variables=var_groups):
        new = first.union(other, remap=True, copy=True)
        print(new)
        new._changed = None  # dirty fix
        new._backup = None  # dirty fix
        first, other = new.split()  # instead of remap
        self_atoms = new.connected_components[: first.connected_components_count][0]
        first_r_groups = MarkushContainer.r_groups_search(first)
        substituent_r_groups = MarkushContainer.r_groups_search(new, exclude=self_atoms)
        for group, self_num in first_r_groups.items():
            if group[0] in variables:
                if other_num := substituent_r_groups.get(group):
                    self_neigbour = [x for x in new.int_adjacency[self_num]][0]
                    other_neigbour = [x for x in new.int_adjacency[other_num]][0]
                    new.delete_atom(self_num)
                    new.delete_atom(other_num)
                    new.add_bond(self_neigbour, other_neigbour, 1)
                    return new
        return first

    @staticmethod
    def connect_no_bond(first: MoleculeContainer, other: MoleculeContainer, variables=var_atoms):
        new = first.union(other, remap=True, copy=True)
        new._changed = None  # dirty fix
        new._backup = None  # dirty fix
        first, other = new.split()  # instead of remap
        first_r_groups = MarkushContainer.r_groups_search(first)
        self_atoms = new.connected_components[: first.connected_components_count][0]
        substituent_r_groups = MarkushContainer.r_groups_search(new, exclude=self_atoms)
        for group, self_num in first_r_groups.items():
            if group[0] in variables:
                if other_num := substituent_r_groups.get(group):
                    other_neigbours = [x for x in new.int_adjacency[other_num]]
                    if len(other_neigbours) != 1:
                        raise ValueError("X groups should have exactly one neighbour")
                    other_X_atom_num = other_neigbours[0]
                    self_bonds = [
                        (x[0], new.bond(self_num, x[0]))
                        for x in new.int_adjacency[self_num].items()
                    ]

                    new.delete_atom(self_num)
                    new.delete_atom(other_num)
                    for atom, bond in self_bonds:
                        new.add_bond(other_X_atom_num, atom, Bond(bond.order))
                return new
        return first

    @property
    def meta(self) -> Dict:
        if self.__meta is None:
            self.__meta = {}  # lazy
        return self.__meta

    @classmethod
    def from_molecule(
        cls, molecule: MoleculeContainer, substituents: list[MoleculeContainer] = None
    ) -> "MarkushContainer":
        obj = cls()
        obj.__initial_mol = molecule
        obj.r_groups = cls.r_groups_search(molecule)
        obj.substituents = substituents
        return obj

    @staticmethod
    def r_groups_search(molecule, exclude: Iterable[int] = []):
        r_groups = {}
        for num, atom in molecule.atoms():
            if atom.atomic_symbol in var_atoms + var_groups and num not in exclude:
                isotope = 0 if atom.isotope is None else atom.isotope
                if 0 <= isotope <= 99:
                    r_groups[(atom.atomic_symbol, isotope)] = num
        return r_groups

    @property
    def r_groups(self):
        return self.__r_groups

    @r_groups.setter
    def r_groups(self, rgroups: dict):
        self.__r_groups = rgroups

    @property
    def subsituents(self):
        return self.__subsituents

    @property
    def initial_mol(self):
        return self.__initial_mol

    @initial_mol.setter
    def initial_mol(self, mol):
        self.__initial_mol = mol

    def generate_molecules(
        self, substituents: list["MoleculeContainer"], exclude_groups: list[str] = []
    ):
        mapping = defaultdict(list)
        for group in self.r_groups:
            if group not in exclude_groups:
                for n_sub, mol in enumerate(substituents):
                    sub_groups = MarkushContainer.r_groups_search(mol)
                    if n_atom := sub_groups.get(group):
                        mapping[group].append(
                            {"group": group, "mol_position": n_sub, "atom_position": n_atom}
                        )
        # forming new molecule without R groups (if possible)
        tmp = []
        for group in self.r_groups:
            if groups := mapping[group]:
                tmp.append(groups)
        for combintaion in product(*tmp):
            new = self.initial_mol
            for modification in combintaion:
                if modification and modification["group"][0] in var_groups:
                    new = MarkushContainer.connect_w_bond(
                        new,
                        self.substituents[modification["mol_position"]],
                        variables=modification["group"],
                    )
                elif modification and modification["group"][0] in var_atoms:
                    new = MarkushContainer.connect_no_bond(
                        new,
                        self.substituents[modification["mol_position"]],
                        variables=modification["group"],
                    )
                else:
                    continue
            yield new

    def copy(self) -> "MarkushContainer":
        copy = super(MoleculeContainer, self).copy()

        copy._bonds = cb = {}
        for n, m_bond in self._bonds.items():
            cb[n] = cbn = {}
            for m, bond in m_bond.items():
                if m in cb:  # bond partially exists. need back-connection.
                    cbn[m] = cb[m][n]
                else:
                    cbn[m] = bond = bond.copy()
                    bond._attach_graph(copy, n, m)

        copy._MarkushiContainer__name = self.__name
        if self.__meta is None:
            copy._MarkushiContainer__meta = None
        else:
            copy._MarkushiContainer__meta = self.__meta.copy()
        copy._plane = self._plane.copy()
        copy._hydrogens = self._hydrogens.copy()
        copy._parsed_mapping = self._parsed_mapping.copy()
        copy._conformers = [c.copy() for c in self._conformers]
        copy._atoms_stereo = self._atoms_stereo.copy()
        copy._allenes_stereo = self._allenes_stereo.copy()
        copy._cis_trans_stereo = self._cis_trans_stereo.copy()
        return copy

    @property
    def substituents(self):
        return self.__substituents

    @substituents.setter
    def substituents(self, substituents: list["MoleculeContainer"]):
        self.__substituents = substituents

    def __str__(self):
        return ".".join([str(self.initial_mol), *[str(x) for x in self.substituents]])


__all__ = ["MarkushContainer"]
