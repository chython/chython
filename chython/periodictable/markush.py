# -*- coding: utf-8 -*-
#
#  Copyright 2024 Timur Gimadiev <timur.gimadiev@gmail.com>
#  Copyright 2024 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from .base import Element


class R(Element):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 6  # same as carbon for proper stereo treatment

    @property
    def isotopes_distribution(self):
        return {n: 0. for n in range(1, 100)}

    @property
    def isotopes_masses(self):
        return {n: 0. for n in range(1, 100)}

    @property
    def _common_valences(self):
        return 0,

    @property
    def _valences_exceptions(self):
        return ()

    @property
    def atomic_radius(self):
        return 1.

    @property
    def mdl_isotope(self):
        return 0

    @property
    def group_number(self):
        return self.__group_number


class X(Element):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 6  # same as carbon for proper stereo treatment

    @property
    def isotopes_distribution(self):
        return {n: 0. for n in range(1, 100)}

    @property
    def isotopes_masses(self):
        return {n: 0. for n in range(1, 100)}

    @property
    def _common_valences(self):
        return 0,

    @property
    def _valences_exceptions(self):
        return ()

    @property
    def atomic_radius(self):
        return 1.

    @property
    def mdl_isotope(self):
        return 0

class Z(Element):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 6  # same as carbon for proper stereo treatment

    @property
    def isotopes_distribution(self):
        return {n: 0. for n in range(1, 100)}

    @property
    def isotopes_masses(self):
        return {n: 0. for n in range(1, 100)}

    @property
    def _common_valences(self):
        return 0,

    @property
    def _valences_exceptions(self):
        return ()

    @property
    def atomic_radius(self):
        return 1.

    @property
    def mdl_isotope(self):
        return 0

__all__ = ['R', 'X', "Z"]
