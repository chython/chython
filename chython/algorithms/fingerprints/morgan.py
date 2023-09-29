# -*- coding: utf-8 -*-
#
#  Copyright 2021-2023 Ramil Nugmanov <nougmanoff@protonmail.com>
#  Copyright 2023 Timur Gimadiev <timur.gimadiev@gmail.com>
#  Copyright 2021 Aleksandr Sizov <murkyrussian@gmail.com>
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
from math import log2
from numpy import uint8, zeros
from typing import Dict, List, Set, TYPE_CHECKING


if TYPE_CHECKING:
    from chython import MoleculeContainer


class MorganFingerprint:
    __slots__ = ()

    def morgan_fingerprint(self, min_radius: int = 1, max_radius: int = 4,
                           length: int = 1024, number_active_bits: int = 2):
        """
        Transform structures into array of binary features.
        Morgan fingerprints. Similar to RDkit implementation.

        :param min_radius: minimal radius of EC
        :param max_radius: maximum radius of EC
        :param length: bit string's length. Should be power of 2
        :param number_active_bits: number of active bits for each hashed tuple

        :return: array(n_features)
        """
        bits = self.morgan_bit_set(min_radius, max_radius, length, number_active_bits)
        fingerprints = zeros(length, dtype=uint8)
        fingerprints[list(bits)] = 1
        return fingerprints

    def morgan_bit_set(self, min_radius: int = 1, max_radius: int = 4,
                       length: int = 1024, number_active_bits: int = 2) -> Set[int]:
        """
        Transform structures into set of indexes of True-valued features.

        :param min_radius: minimal radius of EC
        :param max_radius: maximum radius of EC
        :param length: bit string's length. Should be power of 2
        :param number_active_bits: number of active bits for each hashed tuple
        """
        mask = length - 1
        log = int(log2(length))

        active_bits = set()
        for tpl in self.morgan_hash_set(min_radius, max_radius):
            active_bits.add(tpl & mask)
            if number_active_bits == 2:
                active_bits.add(tpl >> log & mask)
            elif number_active_bits > 2:
                for _ in range(1, number_active_bits):
                    tpl >>= log
                    active_bits.add(tpl & mask)
        return active_bits

    def morgan_hash_set(self: 'MoleculeContainer', min_radius: int = 1, max_radius: int = 4) -> Set[int]:
        """
        Transform structures into integer hashes of atoms with EC.

        :param min_radius: minimal radius of EC
        :param max_radius: maximum radius of EC
        """
        return {x for x in self._morgan_hash_dict(min_radius, max_radius) for x in x.values()}

    def morgan_hash_smiles(self: 'MoleculeContainer', min_radius: int = 1, max_radius: int = 4) -> Dict[int, List[str]]:
        """
        Transform structures into dictionary of hashes of atoms with EC and corresponding SMILES.

        :param min_radius: minimal radius of EC
        :param max_radius: maximum radius of EC
        """
        smiles_dict = defaultdict(set)
        for radius, hash_dict in enumerate(self._morgan_hash_dict(min_radius, max_radius), min_radius - 1):
            for atom, morgan_hash in hash_dict.items():
                smiles_dict[morgan_hash].add(format(self.augmented_substructure((atom,), deep=radius), 'A'))
        return {k: list(v) for k, v in smiles_dict.items()}

    def morgan_smiles_hash(self: 'MoleculeContainer', min_radius: int = 1, max_radius: int = 4) -> Dict[str, List[int]]:
        """
        Transform structures into dictionary of smiles and corresponding hashes of atoms with EC.

        :param min_radius: minimal radius of EC
        :param max_radius: maximum radius of EC
        """
        out = defaultdict(list)
        for k, sl in self.morgan_hash_smiles(min_radius, max_radius).items():
            for s in sl:
                out[s].append(k)
        return dict(out)

    def _morgan_hash_dict(self: 'MoleculeContainer', min_radius: int = 1, max_radius: int = 4) -> List[Dict[int, int]]:
        """
        Transform structures into integer hashes of atoms with EC.
        Returns list of atom-hash pairs for different radii.

        :param min_radius: minimal radius of EC
        :param max_radius: maximum radius of EC
        """
        assert min_radius >= 1, 'min_radius should be positive'
        assert max_radius >= min_radius, 'max_radius should be greater or equal to min_radius'
        identifiers = self._atom_identifiers

        bonds = self._bonds
        out = [identifiers]
        for _ in range(1, max_radius):
            identifiers = {idx: hash((tpl, *(x for x in sorted((int(b), identifiers[ngb])
                                                               for ngb, b in bonds[idx].items()) for x in x)))
                           for idx, tpl in identifiers.items()}
            out.append(identifiers)
        return out[-(max_radius - min_radius + 1):]  # slice [min, max] radii range

    @property
    def _atom_identifiers(self) -> Dict[int, int]:
        raise NotImplementedError


__all__ = ['MorganFingerprint']
