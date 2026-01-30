# -*- coding: utf-8 -*-
#
#  Copyright 2021-2023 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from typing import TYPE_CHECKING
from .linear import *
from .morgan import *


if TYPE_CHECKING:
    from chython import MoleculeContainer, CGRContainer


class Fingerprints(LinearFingerprint, MorganFingerprint):
    __slots__ = ()

    @property
    def _atom_identifiers(self: 'MoleculeContainer'):
        return {idx: hash((atom.isotope or 0, atom.atomic_number, atom.charge, atom.is_radical))
                for idx, atom in self.atoms()}


class FingerprintsCGR(LinearFingerprint, MorganFingerprint):
    __slots__ = ()

    @property
    def _atom_identifiers(self: 'CGRContainer'):
        return {idx: hash((atom.isotope or 0, atom.atomic_number, atom.charge, atom.p_charge,
                           atom.is_radical, atom.p_is_radical))
                for idx, atom in self._atoms.items()}


__all__ = ['Fingerprints', 'FingerprintsCGR']
