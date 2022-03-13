# -*- coding: utf-8 -*-
#
#  Copyright 2021 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from typing import TYPE_CHECKING, Union, List, Dict
from ..tautomers._base import stripped_rules

if TYPE_CHECKING:
    from chython import MoleculeContainer


class Salts:
    __slots__ = ()

    def neutralize_metal_salts(self: 'MoleculeContainer', *, logging=False) -> Union[bool, Dict]:
        """
        Convert metal salts to mixture of metal base and neutral anion form. Works only for stripped salts.

        Example: [K+].CC(=O)[O-] >> [K+].[OH-].CC(=O)O
        Note: do '.neutralize()' procedure before for preventing ambiguous results.

        :param logging: return changed atoms list
        """
        charges = self._charges
        bonds = self._bonds
        hydrogens = self._hydrogens

        metals = []
        for n, a in self._atoms.items():
            if a.is_metal and charges[n] > 0 and not bonds[n]:
                metals.append(n)

        if metals:
            acceptors = set()
            for q in stripped_rules[:-1]:  # except halogenides and hydroxy group.
                for mapping in q.get_mapping(self, automorphism_filter=False):
                    acceptors.add(mapping[1])

            # for imbalanced structures neutralize only part.
            acceptors = list(acceptors)[:sum(charges[n] for n in metals)]
            new_atoms = []
            for n in acceptors:
                hydrogens[n] += 1
                charges[n] += 1
                new_atoms.append(self.add_atom('O', charge=-1))

            if logging:
                return {"charge_changed": acceptors, "new_atoms": new_atoms}
            return bool(acceptors)
        elif logging:
            return []
        return False

    def remove_metals(self: 'MoleculeContainer', *, logging=False) -> Union[bool, List]:
        deleted_atoms = []
        for n, atom in self.atoms():
            if atom.is_metal and not atom.neighbors:
                deleted_atoms.append((n, atom.atomic_symbol))
        if deleted_atoms:
            for n, _ in deleted_atoms:
                self.delete_atom(n)
        if logging:
            return deleted_atoms
        else:
            return bool(deleted_atoms)

    def split_metals(self: 'MoleculeContainer', *, logging=False) -> Union[bool, List]:
        changed = []
        for n, atom in self.atoms():
            if atom.is_metal:
                for neighbour_n, neighbour_atom in list(self.augmented_substructure([n]).atoms()):
                    if not neighbour_atom.is_metal:
                        # Grignard reactivs and Li organics
                        if atom.atomic_symbol in ["Mg", "Zn", "Li"] and neighbour_atom.atomic_symbol == "C":
                            continue
                        # keep M-O-O- type compounds unchanged (peroxides, sulfides, selenides)
                        if neighbour_atom.atomic_symbol in ["S", "Se", "O"] and \
                                self.atom(neighbour_n).neighbors > 1 and \
                                all([neighbour_atom.atomic_symbol == atom.atomic_symbol for i, atom in
                                     list(self.augmented_substructure([neighbour_n]).atoms()) if i != n]):
                            continue
                        atom.charge += 1
                        changed.append(n)
                        self.atom(neighbour_n).charge -= 1
                        changed.append(neighbour_n)
                        self.delete_bond(n, neighbour_n)
        if logging:
            return changed
        else:
            return bool(changed)


__all__ = ['Salts']
