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
from typing import TYPE_CHECKING, Any, Dict, Set, Union

from numpy import uint8, zeros

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

    def morgan_smiles_bit(self, min_radius: int = 1, max_radius: int = 4,
                       length: int = 1024, number_active_bits: int = 2) -> defaultdict[str, set]:
        """
        Transform structures into set of indexes of True-valued features.

        :param min_radius: minimal radius of EC
        :param max_radius: maximum radius of EC
        :param length: bit string's length. Should be power of 2
        :param number_active_bits: number of active bits for each hashed tuple
        """
        mask = length - 1
        log = int(log2(length))

        active_bits = defaultdict(set)
        for smi, tpls in self.morgan_smiles_hash(min_radius, max_radius).items():
            for tpl in tpls:
                active_bits[smi].add(tpl & mask)
                if number_active_bits == 2:
                    active_bits[smi].add(tpl >> log & mask)
                elif number_active_bits > 2:
                    for _ in range(1, number_active_bits):
                        tpl >>= log
                        active_bits[smi].add(tpl & mask)
        return active_bits

    def morgan_hash_set(self: Union["MoleculeContainer"],
                        min_radius: int = 1,  max_radius: int = 4) -> Set[int]:
        """
        Transform structures into integer hashes of atoms with EC.

        :param min_radius: minimal radius of EC
        :param max_radius: maximum radius of EC
        """
        return self.morgan_smiles_hash(min_radius=min_radius, max_radius=max_radius, only_hashes=True)

    def morgan_smiles_hash(self: "MoleculeContainer", min_radius: int = 1, max_radius: int = 4,
                           only_hashes: bool = False) -> Union[Set[int], defaultdict[Any, list[int]]]:
        """
        Transform structures into integer hashes of atoms with EC.

        :param min_radius: minimal radius of EC
        :param max_radius: maximum radius of EC
        :param only_hashes: return dict of substructures and hashes or return only set of hashes
        """
        identifiers = self._atom_identifiers
        bonds = self._bonds
        hash2smi = {}
        arr = set()
        for step in range(1, max_radius + 1):
            tmp_identifiers = {}
            if step >= min_radius:
                if not only_hashes:
                    hash2smi.update(self._hash2smi(identifiers, arr, step))
                arr.update(identifiers.values())
            for atom, tpl in identifiers.items():
                hashes = tuple(x for x in sorted((int(b), identifiers[ngb])
                                                 for ngb, b in bonds[atom].items()) for x in x)
                tmp_identifiers.update({atom: hash((tpl, *hashes))})
            identifiers = tmp_identifiers
        if max_radius > 1:  # add last ring
            if not only_hashes:
                hash2smi.update(self._hash2smi(identifiers, arr, max_radius))
            arr.update(identifiers.values())
        if only_hashes:
            return arr
        else:
            smiles2hash = defaultdict(list)
            for h, smi in hash2smi.items():
                smiles2hash[smi].append(h)
            return smiles2hash

    def _hash2smi(self: "MoleculeContainer", identifiers, arr, radius):
        hash2smi = {}
        for atom, h in identifiers.items():
            if h not in arr:
                smi = format(self.augmented_substructure((atom,), deep=radius - 1), "A")
                hash2smi[h] = smi
        return hash2smi

    @property
    def _atom_identifiers(self) -> Dict[int, int]:
        raise NotImplementedError


__all__ = ["MorganFingerprint"]
