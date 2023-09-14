# -*- coding: utf-8 -*-
#
#  Copyright 2023 Timur Gimadiev <timur.gimadiev@gmail.com>
#  Copyright 2023 Pavel Sidorov <pavel.o.sidorov@gmail.com>
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
from typing import TYPE_CHECKING
from collections import defaultdict

if TYPE_CHECKING:
    from chython import MoleculeContainer


class CircusFingerprint:
    __slots__ = ()

    def circus_hash_set(self: 'MoleculeContainer', min_radius: int = 1, max_radius: int = 4,
                        history: bool = False):
        identifiers = {
            idx: hash((atom.isotope or 0, atom.atomic_number, atom.charge, atom.is_radical))
            for idx, atom in self.atoms()}
        bonds = self._bonds
        smiles2hash = defaultdict(list)
        hash2smi = {}
        arr = set()
        history_log = defaultdict(list)
        for step in range(1, max_radius + 1):
            tmp_identifiers = {}
            if step >= min_radius:
                hash2smi.update(self._hash2smi(identifiers, arr, step))
                arr.update(identifiers.values())
            for atom, tpl in identifiers.items():
                hashes = tuple(x for x in sorted((int(b), identifiers[ngb]) for ngb, b in
                                                 bonds[atom].items()) for x in x)
                if history:
                    history_log[step].append({atom: hashes})
                h = hash((tpl, *hashes))
                tmp_identifiers.update({atom: h})
            identifiers = tmp_identifiers
        if max_radius > 1:  # add last ring
            hash2smi.update(self._hash2smi(identifiers, arr, max_radius))
            arr.update(identifiers.values())
            if history:
                history_log[-1].append(identifiers)
        for h, smi in hash2smi.items():
            smiles2hash[smi].append(h)
        if history:
            return smiles2hash, history_log
        else:
            return smiles2hash

    def _hash2smi(self: 'MoleculeContainer', identifiers, arr, radius):
        hash2smi = {}
        for atom, h in identifiers.items():
            if h not in arr:
                smi = format(self.augmented_substructure((atom,), deep=radius-1), "A")
                hash2smi[h] = smi
        return hash2smi


__all__ = ['CircusFingerprint']
