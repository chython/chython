# -*- coding: utf-8 -*-
#
#  Copyright 2023 Ramil Nugmanov <nougmanoff@protonmail.com>
#  Copyright 2023 Timur Gimadiev <timur.gimadiev@gmail.com>
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
from abc import ABC
from collections import Counter, defaultdict
from math import log2
from numpy import uint8, zeros
from typing import Dict, List, Set, TYPE_CHECKING
from .morgan import MorganFingerprint


if TYPE_CHECKING:
    from chython import MoleculeContainer


class CircusFingerprint(MorganFingerprint):

    __slots__ = ()

    def circus_hash_bit(self, min_radius: int = 1, max_radius: int = 4, length: int = 1024,
                          number_active_bits: int = 2, number_bit_pairs: int = 4) -> Dict[str, set]:
        """
        Transform structures into integer hashes of atoms with EC.

        :param min_radius: minimal radius of EC
        :param max_radius: maximum radius of EC
        """
        assert min_radius >= 1, 'min_radius should be positive'

        smiles_dict = defaultdict(list)
        for radius, hash_dict in enumerate(self._morgan_hash_dict(min_radius, max_radius), min_radius - 1):
            for atom, morgan_hash in hash_dict.items():
                smiles_dict[format(self.augmented_substructure((atom,), deep=radius), 'A')].append(morgan_hash)
        return dict(smiles_dict)

    def circus_smiles_bit(self, min_radius: int = 1, max_radius: int = 4, length: int = 1024,
                          number_active_bits: int = 2, number_bit_pairs: int = 4) -> Dict[str, set]:
        """
        Transform structures into dictionary of smiles and corresponding set of indexes of True-valued features.

        :param min_radius: minimal radius of EC
        :param max_radius: maximum radius of EC
        :param length: bit string's length. Should be power of 2
        :param number_active_bits: number of active bits for each hashed tuple
        :param number_bit_pairs: describe how much repeating hashes we can count in hashable fingerprint (if
            number of fragment in molecule greater or equal this number, we will activate only this number of
            fragments). To take into account all repeating fragments put 0 as a value.
        """
        if not number_bit_pairs:
            number_bit_pairs = 999_999_999  # unreachable count

        assert number_bit_pairs >= 1
        mask = length - 1
        log = int(log2(length))

        active_bits = defaultdict(set)
        for smi, tpls in self.circus_smiles_hash(min_radius, max_radius).items():
            for tpl, cnt in Counter(tpls).items():
                for item in range(min(cnt, number_bit_pairs)):
                    tpl = hash((tpl, item))
                    active_bits[smi].add(tpl & mask)
                    if number_active_bits == 2:
                        active_bits[smi].add(tpl >> log & mask)
                    elif number_active_bits > 2:
                        for _ in range(1, number_active_bits):
                            tpl >>= log
                            active_bits[smi].add(tpl & mask)
        return active_bits

    def circus_smiles_hash(self: 'MoleculeContainer', min_radius: int = 1, max_radius: int = 4) -> \
            Dict[str, list[int]]:
        """
        Transform structures into dictionary of smiles and corresponding hashes of atoms with EC.

        :param min_radius: minimal radius of EC
        :param max_radius: maximum radius of EC
        """
        assert min_radius >= 1, 'min_radius should be positive'

        smiles_dict = defaultdict(list)
        for radius, hash_dict in enumerate(self._morgan_hash_dict(min_radius, max_radius), min_radius-1):
            for atom, morgan_hash in hash_dict.items():
                smiles_dict[format(self.augmented_substructure((atom,), deep=radius), 'A')].append(morgan_hash)
        return dict(smiles_dict)

    def circus_smiles_count(self: 'MoleculeContainer', min_radius: int = 1, max_radius: int = 4) -> \
            Dict[str, list[int]]:
        """
        Transform structures into dictionary of smiles and count of corresponding fragments.

        :param min_radius: minimal radius of EC
        :param max_radius: maximum radius of EC
        """
        return {smi: max([x for x in Counter(hashes).values()]) for smi, hashes in
                self.circus_smiles_hash(min_radius, max_radius).items()}

    def circus_fingerprint(self, min_radius: int = 1, max_radius: int = 4,
                           length: int = 1024, number_active_bits: int = 2):
        """
        Transform structures into array of binary features.
        Morgan fingerprints. Similar to RDkit implementation.

        :param min_radius: minimal radius of EC
        :param max_radius: maximum radius of EC
        :param length: bit string's length. Should be power of 2
        :param number_active_bits: number of active bits for each hashed tuple

        :return: arra
        y(n_features)
        """
        bits = {x for x in self.circus_smiles_bit(min_radius=min_radius, max_radius=max_radius,
                                                  number_active_bits=number_active_bits).values() for x in x}
        fingerprints = zeros(length, dtype=uint8)
        fingerprints[list(bits)] = 1
        return fingerprints