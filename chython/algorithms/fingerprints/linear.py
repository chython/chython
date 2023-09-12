# -*- coding: utf-8 -*-
#
#  Copyright 2021, 2022 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from collections import defaultdict, deque
from joblib import Parallel, delayed
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
    Also count of fragments takes into account by activating multiple bits, but less or equal to 
    `number_bit_pairs`.

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
        :param number_bit_pairs: describe how much repeating fragments we can count in hashable
               fingerprint (if number of fragment in molecule greater or equal this number,
               we will be
               activate only this number of fragments

        :return: array(n_features)
        """
        bits = self.linear_bit_set(min_radius, max_radius, length, number_active_bits,
                                   number_bit_pairs)
        fingerprints = zeros(length, dtype=uint8)
        fingerprints[list(bits)] = 1
        return fingerprints

    def linear_bit_set(self, min_radius: int = 1, max_radius: int = 4,
                       length: int = 1024, number_active_bits: int = 2,
                       number_bit_pairs: int = 4) -> Set[int]:
        """
        Transform structure into set of indexes of True-valued features.

        :param min_radius: minimal length of fragments
        :param max_radius: maximum length of fragments
        :param length: bit string's length. Should be power of 2
        :param number_active_bits: number of active bits for each hashed tuple
        :param number_bit_pairs: describe how much repeating fragments we can count in hashable
               fingerprint (if number of fragment in molecule greater or equal this number,
               we will be activate only this number of fragments
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

    def linear_hash_set(self, min_radius: int = 1, max_radius: int = 4,
                        number_bit_pairs: int = 4) -> Set[int]:
        """
        Transform structure into set of integer hashes of fragments with count information.

        :param min_radius: minimal length of fragments
        :param max_radius: maximum length of fragments
        :param number_bit_pairs: describe how much repeating fragments we can count in hashable
               fingerprint (if number of fragment in molecule greater or equal this number,
               we will be activate only this number of fragments
        """
        return {hash((*tpl, cnt)) for tpl, count in
                self._fragments(min_radius, max_radius).items()
                for cnt in range(min(len(count), number_bit_pairs))}

    def _chains(self: 'MoleculeContainer', min_radius: int = 1, max_radius: int = 4) -> Set[
        Tuple[int, ...]]:
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

    @staticmethod
    def _fragment(atoms: dict[int, int], bonds: dict[int, dict[int, 'Bond']], fragment):
        var = [atoms[fragment[0]]]
        for x, y in zip(fragment, fragment[1:]):
            var.append(int(bonds[x][y]))
            var.append(atoms[y])
        var = tuple(var)
        rev_var = var[::-1]
        if var <= rev_var:
            var = rev_var
        return var

    def _fragments(self: 'MoleculeContainer', min_radius: int = 1, max_radius: int = 4) -> \
            Dict[Tuple[int, ...],
                 List[Tuple[int, ...]]]:
        atoms = {idx: hash((atom.isotope or 0, atom.atomic_number, atom.charge, atom.is_radical))
                 for idx, atom in self.atoms()}

        bonds = self._bonds
        out = defaultdict(list)

        for frag in self._chains(min_radius, max_radius):
            var = self._fragment(atoms, bonds, frag)
            out[var].append(frag)
        return dict(out)

    def _fragments_smiles(self: 'MoleculeContainer', min_radius: int = 1, max_radius: int = 4) -> \
            Dict[str, List[Tuple[int, ...]]]:

        atoms = {idx: hash((atom.isotope or 0, atom.atomic_number, atom.charge, atom.is_radical))
                 for idx, atom in self.atoms()}

        bonds = self._bonds
        out = defaultdict(list)
        smiles_map = {}

        for frag in self._chains(min_radius, max_radius):
            var = self._fragment(atoms, bonds, frag)
            out[var].append(frag)
            if var in out:
                smi_direct = [self._format_atom(frag[0], 0, stereo=False)]
                for x, y in zip(frag, frag[1:]):
                    smi_direct.append(self._format_bond(x, y, 0, stereo=False, aromatic=False))
                    smi_direct.append(self._format_atom(y, 0, stereo=False))
                smi_back = smi_direct[::-1]
                smi_direct = "".join(smi_direct).upper()
                smi_back = "".join(smi_back).upper()
                smi_frag = sorted([smi_direct, smi_back])[0]
                smiles_map[var] = smi_frag
        out = dict(out)
        return {smiles_map[k]: v for k, v in out.items()}


class Fragmentor:
    """
    Fragmentor object can store all seen fragments by .fit() procedure and then produce descriptors
    using dictionary of seen fragments
    """
    def __init__(self, min_radius: int = 1, max_radius: int = 4):
        self.fragments = set()
        self.fragments_map = {}
        self.smi2num = {}
        self.min_radius = min_radius
        self.max_radius = max_radius

    def get_keys(self, molecule):
        desc = molecule._fragments_smiles(self.min_radius, self.max_radius)
        return tuple(desc)

    def get_desciptor_vector(self, molecule, max_count=0):
        desc = molecule._fragments_smiles(self.min_radius, self.max_radius)
        fp = zeros(len(self.fragments), dtype=uint8)
        for key, val in desc.items():
            if key in self.fragments:
                if max_count > 0:
                    fp[self.smi2num[key]] = len(val) if len(val) < max_count else max_count
                else:
                    fp[self.smi2num[key]] = len(val)
        return fp

    def fit(self, molecules: list['MoleculeContainer'], n_jobs=1):
        self.clear()
        self.partial_fit(molecules, n_jobs=n_jobs)

    def partial_fit(self, molecules: list['MoleculeContainer'], n_jobs=1):
        if n_jobs < 2:
            for mol in molecules:
                self.fragments.update(self.get_keys(mol))
        else:
            for calc in Parallel(n_jobs=n_jobs)(delayed(self.get_keys)(molecule) for molecule in
                                                molecules):
                self.fragments.update(calc)
        self.fragments_map = {num: val for num, val in enumerate(self.fragments)}
        self.smi2num = {val: num for num, val in enumerate(self.fragments)}

    def transform(self, molecules, n_jobs=1, max_count=0):
        results = Parallel(n_jobs=n_jobs)(
            delayed(self.get_desciptor_vector)(molecule, max_count) for
            molecule in molecules)
        return results

    def clear(self):
        self.fragments = set()
        self.fragments_map = {}
        self.smi2num = {}

    @property
    def fragments_map_smi_ordered(self):
        return {x: self.fragments_map[x] for x in sorted(self.fragments_map, key=lambda x:
                self.fragments_map[x])}


__all__ = ['LinearFingerprint', 'Fragmentor']
