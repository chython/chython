# -*- coding: utf-8 -*-
#
#  Copyright 2019-2024 Ramil Nugmanov <nougmanoff@protonmail.com>
#  Copyright 2019 Alexander Nikanshin <17071996sasha@gmail.com>
#  Copyright 2019 Tagir Akhmetshin <tagirshin@gmail.com>
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
from .base.groups import GroupV
from .base.periods import PeriodIV, PeriodV, PeriodVI, PeriodVII


class V(Element, PeriodIV, GroupV):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 23

    @property
    def isotopes_distribution(self):
        return {50: 0.0025, 51: 0.9975}

    @property
    def isotopes_masses(self):
        return {50: 49.947163, 51: 50.943964}

    @property
    def _common_valences(self):
        return 0, 2

    @property
    def _valences_exceptions(self):
        return ((2, False, 0, ()), (3, False, 0, ()), (2, False, 0, ((2, 'O'),)),
                (0, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'))),
                (0, False, 0, ((1, 'Cl'), (1, 'Cl'), (1, 'Cl'))),
                (0, False, 0, ((1, 'Br'), (1, 'Br'), (1, 'Br'))),
                (0, False, 0, ((2, 'O'), (1, 'O'))),
                (0, False, 0, ((2, 'O'), (1, 'Cl'))),

                (0, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'))),
                (0, False, 0, ((1, 'Cl'), (1, 'Cl'), (1, 'Cl'), (1, 'Cl'))),
                (0, False, 0, ((1, 'Br'), (1, 'Br'), (1, 'Br'), (1, 'Br'))),
                (0, False, 0, ((2, 'O'), (2, 'O'))),
                (0, False, 0, ((2, 'O'), (1, 'Cl'), (1, 'Cl'))),

                (0, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'))),
                (0, False, 0, ((2, 'O'), (2, 'O'), (1, 'O'))),
                (0, False, 0, ((2, 'S'), (2, 'S'), (1, 'S'))),
                (0, False, 0, ((2, 'O'), (1, 'O'), (1, 'O'), (1, 'O'))),
                (0, False, 0, ((2, 'O'), (1, 'F'), (1, 'F'), (1, 'F'))),
                (0, False, 0, ((2, 'O'), (1, 'Cl'), (1, 'Cl'), (1, 'Cl'))))

    @property
    def atomic_radius(self):
        return 1.71

    @property
    def mdl_isotope(self):
        return 51


class Nb(Element, PeriodV, GroupV):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 41

    @property
    def isotopes_distribution(self):
        return {93: 1.0}

    @property
    def isotopes_masses(self):
        return {93: 92.906378}

    @property
    def _common_valences(self):
        return 0,

    @property
    def _valences_exceptions(self):
        return ((0, False, 0, ((2, 'O'),)),
                (0, False, 0, ((2, 'O'), (1, 'O'))),
                (0, False, 0, ((3, 'N'),)),

                (0, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'))),
                (0, False, 0, ((1, 'Cl'), (1, 'Cl'), (1, 'Cl'), (1, 'Cl'))),
                (0, False, 0, ((1, 'Br'), (1, 'Br'), (1, 'Br'), (1, 'Br'))),
                (0, False, 0, ((2, 'O'), (2, 'O'))),
                (0, False, 0, ((2, 'S'), (2, 'S'))),

                (0, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'))),
                (0, False, 0, ((1, 'Cl'), (1, 'Cl'), (1, 'Cl'), (1, 'Cl'), (1, 'Cl'))),
                (0, False, 0, ((2, 'O'), (1, 'F'), (1, 'F'), (1, 'F'))),
                (0, False, 0, ((2, 'O'), (1, 'Cl'), (1, 'Cl'), (1, 'Cl'))),
                (0, False, 0, ((2, 'O'), (2, 'O'), (1, 'O'))),
                (0, False, 0, ((2, 'O'), (1, 'O'), (1, 'O'), (1, 'O'))))

    @property
    def atomic_radius(self):
        return 1.98

    @property
    def mdl_isotope(self):
        return 93


class Ta(Element, PeriodVI, GroupV):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 73

    @property
    def isotopes_distribution(self):
        return {180: 0.00012, 181: 0.99988}

    @property
    def isotopes_masses(self):
        return {180: 179.947466, 181: 180.947996}

    @property
    def _common_valences(self):
        return 0,

    @property
    def _valences_exceptions(self):
        return ((0, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'))),
                (0, False, 0, ((1, 'Cl'), (1, 'Cl'), (1, 'Cl'), (1, 'Cl'), (1, 'Cl'))),
                (0, False, 0, ((2, 'O'), (1, 'F'), (1, 'F'), (1, 'F'))),
                (0, False, 0, ((2, 'O'), (1, 'Cl'), (1, 'Cl'), (1, 'Cl'))),
                (0, False, 0, ((2, 'O'), (2, 'O'), (1, 'O'))),
                (0, False, 0, ((2, 'O'), (1, 'O'), (1, 'O'), (1, 'O'))))

    @property
    def atomic_radius(self):
        return 2.0

    @property
    def mdl_isotope(self):
        return 181


class Db(Element, PeriodVII, GroupV):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 105

    @property
    def isotopes_distribution(self):
        return {268: 1.0}

    @property
    def isotopes_masses(self):
        return {268: 268.125676}

    @property
    def _common_valences(self):
        return 0,

    @property
    def _valences_exceptions(self):
        return ()

    @property
    def atomic_radius(self):
        return 2.0  # unknown, taken radius of previous element in group

    @property
    def mdl_isotope(self):
        return 270


__all__ = ['V', 'Nb', 'Ta', 'Db']
