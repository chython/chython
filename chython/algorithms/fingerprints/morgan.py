# -*- coding: utf-8 -*-
#
#  Copyright 2021 Aleksandr Sizov <murkyrussian@gmail.com>
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
from math import log2
from numpy import zeros, uint8
from typing import Set, TYPE_CHECKING


if TYPE_CHECKING:
    from chython import MoleculeContainer


class MorganFingerprint:
    __slots__ = ()

    def morgan_fingerprint(self, min_radius: int = 1, max_radius: int = 4,
                           length: int = 1024, number_active_bits: int = 2,
                           include_hydrogens: bool = True, with_pharmacophores: bool = False):
        """
        Transform structures into array of binary features.
        Morgan fingerprints. Similar to RDkit implementation.

        :param min_radius: minimal radius of EC
        :param max_radius: maximum radius of EC
        :param length: bit string's length. Should be power of 2
        :param number_active_bits: number of active bits for each hashed tuple
        :param include_hydrogens: take into account hydrogen atoms
        :param with_pharmacophores: use pharmacophoric features to identify and distinguish each atom in an molecule

        :return: array(n_features)
        """
        bits = self.morgan_bit_set(min_radius, max_radius, length, number_active_bits, include_hydrogens,
                                   with_pharmacophores)
        fingerprints = zeros(length, dtype=uint8)
        fingerprints[list(bits)] = 1
        return fingerprints

    def morgan_bit_set(self, min_radius: int = 1, max_radius: int = 4,
                       length: int = 1024, number_active_bits: int = 2,
                       include_hydrogens: bool = True, with_pharmacophores: bool = False) -> Set[int]:
        """
        Transform structures into set of indexes of True-valued features.

        :param min_radius: minimal radius of EC
        :param max_radius: maximum radius of EC
        :param length: bit string's length. Should be power of 2
        :param number_active_bits: number of active bits for each hashed tuple
        :param include_hydrogens: take into account hydrogen atoms
        :param with_pharmacophores: use pharmacophoric features to identify and distinguish each atom in an molecule
        """
        mask = length - 1
        log = int(log2(length))

        active_bits = set()
        for tpl in self.morgan_hash_set(min_radius, max_radius, include_hydrogens, with_pharmacophores):
            active_bits.add(tpl & mask)
            if number_active_bits == 2:
                active_bits.add(tpl >> log & mask)
            elif number_active_bits > 2:
                for _ in range(1, number_active_bits):
                    tpl >>= log
                    active_bits.add(tpl & mask)
        return active_bits

    def morgan_hash_set(self: 'MoleculeContainer', min_radius: int = 1, max_radius: int = 4,
                        include_hydrogens: bool = True, with_pharmacophores: bool = False) -> Set[int]:
        """
        Transform structures into integer hashes of atoms with EC.

        :param min_radius: minimal radius of EC
        :param max_radius: maximum radius of EC
        :param include_hydrogens: take into account hydrogen atoms
        :param with_pharmacophores: use pharmacophoric features to identify and distinguish each atom in an molecule
        """
        if with_pharmacophores:
            identifiers = self.pharmacophores
        elif include_hydrogens:
            identifiers = {idx: int(atom) for idx, atom in self.atoms()}
        else:
            identifiers = {idx: hash((atom.isotope or 0, atom.atomic_number, atom.charge, atom.is_radical))
                           for idx, atom in self.atoms()}

        bonds = self._bonds
        arr = set()
        for step in range(1, max_radius + 1):
            if step >= min_radius:
                arr.update(identifiers.values())
            identifiers = {idx: hash((tpl, *(x for x in
                                             sorted((int(b), identifiers[ngb]) for ngb, b in bonds[idx].items())
                                             for x in x)))
                           for idx, tpl in identifiers.items()}

        if max_radius > 1:  # add last ring
            arr.update(identifiers.values())
        return arr


__all__ = ['MorganFingerprint']
