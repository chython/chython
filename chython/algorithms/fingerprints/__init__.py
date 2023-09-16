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
from .circus import *
from .linear import *
from .morgan import *


class FingerprintsMol(LinearFingerprint, MorganFingerprint, CircusFingerprint):
    __slots__ = ()

    def _atom2identifiers(self, atom):
        return atom.isotope or 0, atom.atomic_number, atom.charge, atom.is_radical


class FingerprintsCGR(LinearFingerprint, MorganFingerprint, CircusFingerprint):
    __slots__ = ()

    def _atom2identifiers(self, atom):
        return atom.isotope or 0, atom.atomic_number, atom.charge, atom.is_radical, \
               atom.p_is_radical, atom.p_charge


__all__ = ['FingerprintsMol', 'FingerprintsCGR']
