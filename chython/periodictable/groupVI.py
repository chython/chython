# -*- coding: utf-8 -*-
#
#  Copyright 2019-2024 Ramil Nugmanov <nougmanoff@protonmail.com>
#  Copyright 2019 Tagir Akhmetshin <tagirshin@gmail.com>
#  Copyright 2019 Dayana Bashirova <dayana.bashirova@yandex.ru>
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
from .base.groups import GroupVI
from .base.periods import PeriodIV, PeriodV, PeriodVI, PeriodVII


class Cr(Element, PeriodIV, GroupVI):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 24

    @property
    def isotopes_distribution(self):
        return {50: 0.04345, 51: 0., 52: 0.83789, 53: 0.09501, 54: 0.02365}

    @property
    def isotopes_masses(self):
        return {50: 49.946050, 51: 50.944767, 52: 51.940512, 53: 52.940654, 54: 53.938885}

    @property
    def _common_valences(self):
        return 0, 2, 3

    @property
    def _valences_exceptions(self):
        return ((2, False, 0, ()), (3, False, 0, ()),
                (0, False, 0, ((2, 'O'), (2, 'O'))),  # CrO2
                (0, False, 0, ((2, 'O'), (1, 'O'), (1, 'O'))),
                (0, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'))),  # CrF4
                (0, False, 0, ((1, 'C'), (1, 'C'), (1, 'C'), (1, 'C'))),  # CrC4

                (0, False, 0, ((2, 'O'), (2, 'O'), (2, 'O'))),
                (0, False, 0, ((2, 'O'), (2, 'O'), (1, 'O'), (1, 'O'))),
                (0, False, 0, ((2, 'O'), (2, 'O'), (1, 'O'), (1, 'Cl'))),
                (0, False, 0, ((2, 'O'), (2, 'O'), (1, 'Cl'), (1, 'Cl'))))

    @property
    def atomic_radius(self):
        return 1.66

    @property
    def mdl_isotope(self):
        return 52


class Mo(Element, PeriodV, GroupVI):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 42

    @property
    def isotopes_distribution(self):
        return {92: 0.1484, 94: 0.0925, 95: 0.1592, 96: 0.1668, 97: 0.0955, 98: 0.2413, 99: 0., 100: 0.0963}

    @property
    def isotopes_masses(self):
        return {92: 91.906810, 94: 93.905088, 95: 94.905841, 96: 95.904679, 97: 96.906021, 98: 97.905408, 99: 98.907712,
                100: 99.907477}

    @property
    def _common_valences(self):
        return 0,

    @property
    def _valences_exceptions(self):
        return ((0, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'))),
                (0, False, 0, ((1, 'Cl'), (1, 'Cl'), (1, 'Cl'), (1, 'Cl'))),
                (0, False, 0, ((1, 'Br'), (1, 'Br'), (1, 'Br'), (1, 'Br'))),
                (0, False, 0, ((2, 'O'), (2, 'O'))),
                (0, False, 0, ((2, 'S'), (2, 'S'))),

                (0, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'))),
                (0, False, 0, ((1, 'Cl'), (1, 'Cl'), (1, 'Cl'), (1, 'Cl'), (1, 'Cl'))),
                (0, False, 0, ((1, 'Br'), (1, 'Br'), (1, 'Br'), (1, 'Br'), (1, 'Br'))),

                (0, False, 0, ((2, 'O'), (2, 'O'), (2, 'O'))),
                (0, False, 0, ((2, 'O'), (2, 'O'), (1, 'O'), (1, 'O'))),

                (0, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'))),
                (0, False, 0, ((2, 'O'), (1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'))))

    @property
    def atomic_radius(self):
        return 1.90

    @property
    def mdl_isotope(self):
        return 96


class W(Element, PeriodVI, GroupVI):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 74

    @property
    def isotopes_distribution(self):
        return {180: 0.0012, 182: 0.265, 183: 0.1431, 184: 0.3064, 186: 0.2843}

    @property
    def isotopes_masses(self):
        return {180: 179.946706, 182: 181.948206, 183: 182.950224, 184: 183.950933, 186: 185.954362}

    @property
    def _common_valences(self):
        return 0,

    @property
    def _valences_exceptions(self):
        return ((0, False, 0, ((2, 'O'), (2, 'O'), (2, 'O'))),
                (0, False, 0, ((2, 'O'), (2, 'O'), (1, 'O'), (1, 'O'))),

                (0, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'))),
                (0, False, 0, ((1, 'Cl'), (1, 'Cl'), (1, 'Cl'), (1, 'Cl'), (1, 'Cl'), (1, 'Cl'))),
                (0, False, 0, ((2, 'O'), (1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'))))

    @property
    def atomic_radius(self):
        return 1.93

    @property
    def mdl_isotope(self):
        return 184


class Sg(Element, PeriodVII, GroupVI):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 106

    @property
    def isotopes_distribution(self):
        return {269: 1.0}

    @property
    def isotopes_masses(self):
        return {269: 269.128634}

    @property
    def _common_valences(self):
        return 0,

    @property
    def _valences_exceptions(self):
        return ()

    @property
    def atomic_radius(self):
        return 1.93  # unknown, taken radius of previous element in group

    @property
    def mdl_isotope(self):
        return 269


__all__ = ['Cr', 'Mo', 'W', 'Sg']
