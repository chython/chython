# -*- coding: utf-8 -*-
#
#  Copyright 2021-2023 Ramil Nugmanov <nougmanoff@protonmail.com>
#  Copyright 2021 Aleksandr Sizov <murkyrussian@gmail.com>
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
from collections import defaultdict, deque
from math import log2
from numpy import uint8, zeros
from typing import Deque, Dict, List, Set, Tuple, TYPE_CHECKING


if TYPE_CHECKING:
    from chython import MoleculeContainer


class LinearFingerprint:
    __slots__ = ()
    """
    Linear fragments fingerprints.
    Transform structures into fingerprints based on linear fragments descriptors.
    Also count of fragments takes into account by activating multiple bits,
    but less or equal to `number_bit_pairs`.To take into account
    all repeating fragments put 0 as a value of `number_bit_pairs` parameter.

    For example `CC` fragment found 4 times and `number_bit_pairs` is set to 3.
    In this case will be activated 3 bits: for count 1, for count 2 and for count 3.
    This gives intersection in bits with another structure with only 2 `CC` fragments.
    """

    def linear_fingerprint(self, min_radius: int = 1, max_radius: int = 4,
                           length: int = 1024, number_active_bits: int = 2,
                           number_bit_pairs: int = 4):
        """
        Transform structures into array of binary features.

        :param min_radius: minimal length of fragments
        :param max_radius: maximum length of fragments
        :param length: bit string's length. Should be power of 2
        :param number_active_bits: number of active bits for each hashed tuple
        :param number_bit_pairs: describe how much repeating fragments we can count in hashable fingerprint (if
            number of fragment in molecule greater or equal this number, we will activate only this number of
            fragments). To take into account all repeating fragments put 0 as a value.

        :return: array(n_features)
        """
        bits = self.linear_bit_set(min_radius, max_radius, length, number_active_bits,
                                   number_bit_pairs)
        fingerprints = zeros(length, dtype=uint8)
        fingerprints[list(bits)] = 1
        return fingerprints

    def linear_bit_set(self, min_radius: int = 1, max_radius: int = 4, length: int = 1024, number_active_bits: int = 2,
                       number_bit_pairs: int = 4) -> Set[int]:
        """
        Transform structure into set of indexes of True-valued features.

        :param min_radius: minimal length of fragments
        :param max_radius: maximum length of fragments
        :param length: bit string's length. Should be power of 2
        :param number_active_bits: number of active bits for each hashed tuple
        :param number_bit_pairs: describe how much repeating fragments we can count in hashable fingerprint (if
            number of fragment in molecule greater or equal this number, we will activate only this number of
            fragments). To take into account all repeating fragments put 0 as a value.
        """
        mask = length - 1
        log = int(log2(length))

        hashes = self.linear_hash_set(min_radius, max_radius, number_bit_pairs)
        active_bits = set()
        for tpl in hashes:
            active_bits.add(tpl & mask)
            if number_active_bits == 2:
                active_bits.add((tpl >> log) & mask)
            elif number_active_bits > 2:
                for _ in range(1, number_active_bits):
                    tpl >>= log  # shift
                    active_bits.add(tpl & mask)
        return active_bits

    def linear_hash_set(self, min_radius: int = 1, max_radius: int = 4, number_bit_pairs: int = 4) -> Set[int]:
        """
        Transform structure into set of integer hashes of fragments with count information.

        :param min_radius: minimal length of fragments
        :param max_radius: maximum length of fragments
        :param number_bit_pairs: describe how much repeating fragments we can count in hashable fingerprint (if
            number of fragment in molecule greater or equal this number, we will activate only this number of
            fragments). To take into account all repeating fragments put 0 as a value.
        """
        if not number_bit_pairs:
            number_bit_pairs = 999_999_999  # unreachable count

        return {hash((*tpl, cnt)) for tpl, count in
                self._fragments(min_radius, max_radius).items()
                for cnt in range(min(len(count), number_bit_pairs))}

    def linear_hash_smiles(self: 'MoleculeContainer', min_radius: int = 1, max_radius: int = 4,
                           number_bit_pairs: int = 4) -> Dict[int, List[str]]:
        """
        Transform structure into dict of integer hashes of fragments with count information and
            corresponding fragment SMILES.

        :param min_radius: minimal length of fragments
        :param max_radius: maximum length of fragments
        :param number_bit_pairs: describe how much repeating fragments we can count in hashable fingerprint (if
            number of fragment in molecule greater or equal this number, we will activate only this number of
            fragments). To take into account all repeating fragments put 0 as a value.
        """
        if not number_bit_pairs:
            number_bit_pairs = 999_999_999  # unreachable count

        out = defaultdict(set)
        for frg, chains in self._fragments(min_radius, max_radius).items():
            chain = chains[0]
            smiles = [self._format_atom(chain[0], None, stereo=False)]
            for x, y in zip(chain, chain[1:]):
                smiles.append(self._format_bond(x, y, None, stereo=False, aromatic=False))
                smiles.append(self._format_atom(y, None, stereo=False))
            smiles = ''.join(smiles)

            for cnt in range(min(len(chains), number_bit_pairs)):
                out[hash((*frg, cnt))].add(smiles)  # collisions possible
        return {k: list(v) for k, v in out.items()}

    def linear_smiles_hash(self, min_radius: int = 1, max_radius: int = 4,
                           number_bit_pairs: int = 4) -> Dict[str, List[int]]:
        """
        Transform structure into dict of fragment SMILES and list of corresponding integer hashes of fragments.

        :param min_radius: minimal length of fragments
        :param max_radius: maximum length of fragments
        :param number_bit_pairs: describe how much repeating fragments we can count in hashable fingerprint (if
            number of fragment in molecule greater or equal this number, we will activate only this number of
            fragments). To take into account all repeating fragments put 0 as a value.
        """
        out = defaultdict(list)
        for k, sl in self.linear_hash_smiles(min_radius, max_radius, number_bit_pairs).items():
            for s in sl:
                out[s].append(k)
        return dict(out)

    def _chains(self: 'MoleculeContainer', min_radius: int = 1, max_radius: int = 4) -> Set[Tuple[int, ...]]:
        queue: Deque[Tuple[int, ...]]  # typing
        atoms = self._atoms
        bonds = self._bonds

        if min_radius == 1:
            arr = {(x,) for x in atoms}
            if max_radius == 1:  # special case
                return arr
            else:
                queue = deque(arr)
        else:
            arr = set()
            queue = deque((x,) for x in atoms)

        while queue:
            now = queue.popleft()
            var = [now + (x,) for x in bonds[now[-1]] if x not in now]
            if var:
                if len(var[0]) < max_radius:
                    queue.extend(var)
                if len(var[0]) >= min_radius:
                    for frag in var:
                        rev = frag[::-1]
                        arr.add(frag if frag > rev else rev)
        return arr

    def _fragments(self: 'MoleculeContainer', min_radius: int = 1,
                   max_radius: int = 4) -> Dict[Tuple[int, ...], List[Tuple[int, ...]]]:
        atoms = self._atom_identifiers
        bonds = self._bonds
        out = defaultdict(list)

        for frag in self._chains(min_radius, max_radius):
            var = [atoms[frag[0]]]
            for x, y in zip(frag, frag[1:]):
                var.append(int(bonds[x][y]))
                var.append(atoms[y])
            var = tuple(var)
            rev_var = var[::-1]
            if var > rev_var:
                out[var].append(frag)
            else:
                out[rev_var].append(frag[::-1])
        return dict(out)

    @property
    def _atom_identifiers(self) -> Dict[int, int]:
        raise NotImplementedError


__all__ = ['LinearFingerprint']
