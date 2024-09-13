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
from ..periodictable import MarkushiElement
from .bonds import Bond
from .molecule import MoleculeContainer
from collections import defaultdict
from itertools import product
from typing import Iterable, Dict


class MarkushiContainer(MoleculeContainer):
    #__slots__ = ('__meta')
    __meta = None
    __name = None
    _substituents: list['str'] = []

    @property
    def meta(self) -> Dict:
        if self.__meta is None:
            self.__meta = {}  # lazy
        return self.__meta

    @property
    def R_groups(self):
        groups = {}
        for num, atom in self.atoms():
            if isinstance(atom, MarkushiElement):  # searching  other_R.update({atom.atomic_symbol: num})
                groups.update({atom.atomic_symbol: num})
        return groups

    @property
    def R_groups_count(self):
        return len(self.R_groups)

    def __add__(self, other):
        new = self | other
        self_R = {}
        other_R = {}
        self_atoms = new.connected_components[:self.connected_components_count]
        for num in [x for x in self_atoms for x in x]:
            atom = new.atom(num)
            if isinstance(atom, MarkushiElement):  # searching for R element
                self_R.update({atom.atomic_symbol: num})
        other_atoms = new.connected_components[self.connected_components_count:]
        for num in [x for x in other_atoms for x in x]:
            atom = new.atom(num)
            if isinstance(atom, MarkushiElement):  # searching for R element
                other_R.update({atom.atomic_symbol: num})

        if matching_groups := set(x for x in self_R).intersection(x for x in other_R):
            for R in matching_groups:
                n_neigbour = [x for x in new.int_adjacency[self_R[R]]][0]
                m_neigbour = [x for x in new.int_adjacency[other_R[R]]][0]
                new.delete_atom(self_R[R])
                new.delete_atom(other_R[R])
                new.add_bond(n_neigbour, m_neigbour, 1)
            return new
        else:
            return self

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

    def copy(self) -> 'MarkushiContainer':
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
        self_str = super(MarkushiContainer, self).__str__()
        return ".".join([self_str, *self.substituents])

    def substructure(self, atoms: Iterable[int], *, as_query: bool = False, recalculate_hydrogens=True,
                     skip_neighbors_marks=False, skip_hybridizations_marks=False, skip_hydrogens_marks=False,
                     skip_rings_sizes_marks=False, skip_heteroatoms_marks=False) -> 'MarkushiContainer':
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
        atoms_type = {x: self.atom(x) for x in atoms}
        if as_query:
            raise ValueError("No query for MarkushiContainer")
        bond_type = Bond
        sub = object.__new__(self.__class__)
        sub._MarkushiContainer__name = sub._MarkushiContainer__meta = None

        sa = self._atoms
        sb = self._bonds
        sc = self._charges
        sr = self._radicals

        sub._charges = {n: sc[n] for n in atoms}
        sub._radicals = {n: sr[n] for n in atoms}

        sub._atoms = ca = {}
        for n in atoms:
            ca[n] = atom = atoms_type[n].from_atom(sa[n])
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


__all__ = ['MarkushiContainer']