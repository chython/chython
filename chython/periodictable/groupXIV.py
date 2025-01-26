# -*- coding: utf-8 -*-
#
#  Copyright 2019-2024 Ramil Nugmanov <nougmanoff@protonmail.com>
#  Copyright 2019 Dayana Bashirova <dayana.bashirova@yandex.ru>
#  Copyright 2019 Tansu Nasyrova <tansu.nasyrova@gmail.com>
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
from .base.groups import GroupXIV
from .base.periods import PeriodII, PeriodIII, PeriodIV, PeriodV, PeriodVI, PeriodVII


class C(Element, PeriodII, GroupXIV):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 6

    @property
    def isotopes_distribution(self):
        return {11: 0., 12: 0.9893, 13: 0.0107, 14: 0.}

    @property
    def isotopes_masses(self):
        return {11: 11.011432, 12: 12.0, 13: 13.003355, 14: 14.003242}

    @property
    def _common_valences(self):
        return 4,

    @property
    def _valences_exceptions(self):
        return (0, True, 3, ()), (1, False, 3, ()), (-1, False, 3, ()), (0, False, 0, ())

    @property
    def atomic_radius(self):
        return .67

    @property
    def mdl_isotope(self):
        return 12

    @property
    def is_forming_single_bonds(self):
        return True

    @property
    def is_forming_double_bonds(self):
        return True


class Si(Element, PeriodIII, GroupXIV):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 14

    @property
    def isotopes_distribution(self):
        return {28: 0.922297, 29: 0.046832, 30: 0.030872}

    @property
    def isotopes_masses(self):
        return {28: 27.976927, 29: 28.976495, 30: 29.97377}

    @property
    def _common_valences(self):
        return 4,

    @property
    def _valences_exceptions(self):
        return (-2, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'))), (0, False, 0, ())

    @property
    def atomic_radius(self):
        return 1.11

    @property
    def mdl_isotope(self):
        return 28

    @property
    def is_forming_single_bonds(self):
        return True

    @property
    def is_forming_double_bonds(self):
        return True


class Ge(Element, PeriodIV, GroupXIV):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 32

    @property
    def isotopes_distribution(self):
        return {70: 0.2084, 72: 0.2754, 73: 0.0773, 74: 0.3628, 76: 0.0761}

    @property
    def isotopes_masses(self):
        return {70: 69.92425, 72: 71.922076, 73: 72.923459, 74: 73.921178, 76: 75.921403}

    @property
    def _common_valences(self):
        return 4,

    @property
    def _valences_exceptions(self):
        return (-2, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'))), (0, False, 0, ())

    @property
    def atomic_radius(self):
        return 1.25

    @property
    def mdl_isotope(self):
        return 73

    @property
    def is_forming_single_bonds(self):
        return True

    @property
    def is_forming_double_bonds(self):
        return True


class Sn(Element, PeriodV, GroupXIV):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 50

    @property
    def isotopes_distribution(self):
        return {112: 0.0097, 113: 0., 114: 0.0066, 115: 0.0034, 116: 0.1454, 117: 0.0768, 118: 0.2422, 119: 0.0859,
                120: 0.3258, 122: 0.0463, 124: 0.0579}

    @property
    def isotopes_masses(self):
        return {112: 111.904821, 113: 112.905171, 114: 113.902782, 115: 114.903346, 116: 115.901744, 117: 116.902954,
                118: 117.901606, 119: 118.903309, 120: 119.902197, 122: 121.903440, 124: 123.905275}

    @property
    def _common_valences(self):
        return 0, 4

    @property
    def _valences_exceptions(self):
        return ((2, False, 0, ()),
                (0, False, 0, ((2, 'O'),)), (0, False, 0, ((1, 'O'), (1, 'O'))),
                (0, False, 0, ((2, 'S'),)),
                (0, False, 0, ((1, 'Cl'), (1, 'Cl'))),

                (1, False, 0, ((1, 'C'), (1, 'C'), (1, 'C'))),
                (0, False, 1, ((1, 'C'), (1, 'C'), (1, 'C'))),

                (-2, False, 0, ((1, 'O'), (1, 'O'), (1, 'O'), (1, 'O'), (1, 'O'), (1, 'O'))))

    @property
    def atomic_radius(self):
        return 1.45

    @property
    def mdl_isotope(self):
        return 119


class Pb(Element, PeriodVI, GroupXIV):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 82

    @property
    def isotopes_distribution(self):
        return {204: 0.014, 206: 0.241, 207: 0.221, 208: 0.524, 210: 0.}

    @property
    def isotopes_masses(self):
        return {204: 203.973029, 206: 205.974449, 207: 206.975881, 208: 207.976636, 210: 209.984189}

    @property
    def _common_valences(self):
        return 0, 2

    @property
    def _valences_exceptions(self):
        return ((2, False, 0, ()),
                (-2, False, 0, ((1, 'O'), (1, 'O'), (1, 'O'), (1, 'O'))),
                (-2, False, 0, ((2, 'O'), (1, 'O'), (1, 'O'))),

                (0, False, 0, ((2, 'O'), (2, 'O'))),
                (0, False, 0, ((1, 'O'), (1, 'O'), (1, 'O'), (1, 'O'))),
                (0, False, 0, ((2, 'O'), (1, 'O'), (1, 'O'))),

                (0, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'))),
                (0, False, 0, ((1, 'Cl'), (1, 'Cl'), (1, 'Cl'), (1, 'Cl'))),
                (0, False, 0, ((1, 'C'), (1, 'C'), (1, 'C'), (1, 'C'))))

    @property
    def atomic_radius(self):
        return 1.54

    @property
    def mdl_isotope(self):
        return 207


class Fl(Element, PeriodVII, GroupXIV):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 114

    @property
    def isotopes_distribution(self):
        return {289: 1.0}

    @property
    def isotopes_masses(self):
        return {289: 289.190444}

    @property
    def _common_valences(self):
        return 0,

    @property
    def _valences_exceptions(self):
        return ()

    @property
    def atomic_radius(self):
        return 1.54  # unknown, taken radius of previous element in group

    @property
    def mdl_isotope(self):
        return 289


__all__ = ['C', 'Si', 'Ge', 'Sn', 'Pb', 'Fl']
