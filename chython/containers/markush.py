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
from typing import Iterable, Dict
from .molecule import MoleculeContainer


class MarkushContainer:
    __slots__ = ('__meta', '__tree')

    def __init__(self):
        super().__init__()
        self.__tree = {}
        self.__meta = None

    @property
    def meta(self) -> Dict:
        if self.__meta is None:
            self.__meta = {}  # lazy
        return self.__meta

    @classmethod
    def from_molecule(cls, molecule: MoleculeContainer) -> 'MarkushContainer':
        obj = cls()
        ...
        return obj

    def find_matchs(self, mol):
        R_mapping = {}
        for atom_num, atom in mol.atoms():
            if matched_group := self.R_groups.get(atom.atomic_symbol):
                R_mapping[matched_group] = atom_num
        return R_mapping

    def generate_mapping(self, fragments: list['MarkushiContainer']):
        R_mapping = defaultdict(list)
        for num, fragment in enumerate(fragments):
            for core_atom, frag_atom in self.find_matchs(fragment).items():
                R_mapping[core_atom].append((num, frag_atom))
        return R_mapping

    def sequence(self, fragments: list['MarkushiContainer']):
        return product(*[x for x in self.generate_mapping(fragments).values()])

    def generate_molecules(self, fragments: list['MarkushiContainer']):
        for groups in self.sequence(fragments):
            # forming new molecule without R groups (if possible)
            new = self.copy()
            for _, (num_frag, _) in zip(self.R_groups.items(), groups):
                new = new+fragments[num_frag]
            yield new

    def copy(self) -> 'MarkushContainer':
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
        return self._substituents

    @substituents.setter
    def substituents(self, substituents: list['MarkushiContainer']):
        self._substituents = [str(x) for x in substituents]

    def __str__(self):
        return ".".join([self_str, *self.substituents])


__all__ = ['MarkushContainer']
