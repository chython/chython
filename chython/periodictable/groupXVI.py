# -*- coding: utf-8 -*-
#
#  Copyright 2019-2024 Ramil Nugmanov <nougmanoff@protonmail.com>
#  Copyright 2019 Dayana Bashirova <dayana.bashirova@yandex.ru>
#  Copyright 2019 Tagir Akhmetshin <tagirshin@gmail.com>
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
from .base.groups import GroupXVI
from .base.periods import PeriodII, PeriodIII, PeriodIV, PeriodV, PeriodVI, PeriodVII


class O(Element, PeriodII, GroupXVI):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 8

    @property
    def isotopes_distribution(self):
        return {15: 0., 16: 0.99757, 17: 0.00038, 18: 0.00205}

    @property
    def isotopes_masses(self):
        return {15: 15.003065, 16: 15.994915, 17: 16.999132, 18: 17.99916}

    @property
    def _common_valences(self):
        return 2,

    @property
    def _valences_exceptions(self):
        return (-1, False, 1, ()), (-2, False, 0, ()), (0, True, 1, ()), (1, False, 3, ())

    @property
    def atomic_radius(self):
        return .48

    @property
    def mdl_isotope(self):
        return 16

    @property
    def is_forming_single_bonds(self):
        return True

    @property
    def is_forming_double_bonds(self):
        return True


class S(Element, PeriodIII, GroupXVI):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 16

    @property
    def isotopes_distribution(self):
        return {32: 0.9493, 33: 0.0076, 34: 0.0429, 35: 0., 36: 0.0002}

    @property
    def isotopes_masses(self):
        return {32: 31.972071, 33: 32.971458, 34: 33.967867, 35: 34.969032, 36: 35.967081}

    @property
    def _common_valences(self):
        return 2,

    @property
    def _valences_exceptions(self):
        return ((-1, False, 1, ()), (-2, False, 0, ()),  # anions
                (0, True, 1, ()), (0, True, 3, ()),  # radical OGB-dataset
                (0, False, 0, ()),  # elemental

                (1, False, 0, ((2, 'C'), (1, 'C'))),
                (1, False, 0, ((2, 'C'), (1, 'S'))),
                (1, False, 0, ((2, 'N'), (1, 'C'))),

                (1, False, 0, ((1, 'C'), (1, 'C'), (1, 'C'))),
                (1, False, 0, ((1, 'C'), (1, 'C'), (1, 'B'))),
                (1, False, 0, ((1, 'C'), (1, 'C'), (1, 'O'))),
                (1, False, 0, ((1, 'C'), (1, 'C'), (1, 'N'))),

                (1, False, 0, ((2, 'O'), (1, 'C'), (1, 'C'), (1, 'C'))),
                (1, False, 0, ((2, 'O'), (1, 'C'), (1, 'C'), (1, 'N'))),

                (0, False, 0, ((2, 'O'), (2, 'O'))),
                (0, False, 0, ((2, 'O'), (2, 'N'))),
                (0, False, 0, ((2, 'N'), (2, 'N'))),
                (0, False, 0, ((2, 'O'), (2, 'C'))),
                (0, False, 0, ((2, 'C'), (2, 'C'))),
                (0, False, 0, ((2, 'C'), (2, 'N'))),

                (0, False, 0, ((2, 'O'), (1, 'O'), (1, 'O'))),
                (0, False, 0, ((2, 'O'), (1, 'O'), (1, 'S'))),
                (0, False, 0, ((2, 'O'), (1, 'O'), (1, 'N'))),
                (0, False, 0, ((2, 'O'), (1, 'O'), (1, 'F'))),
                (0, False, 0, ((2, 'O'), (1, 'O'), (1, 'Cl'))),
                (0, False, 0, ((2, 'O'), (1, 'O'), (1, 'Br'))),
                (0, False, 0, ((2, 'O'), (1, 'O'), (1, 'C'))),

                (0, False, 0, ((2, 'O'), (1, 'N'), (1, 'N'))),
                (0, False, 0, ((2, 'O'), (1, 'N'), (1, 'Cl'))),
                (0, False, 0, ((2, 'O'), (1, 'N'), (1, 'C'))),
                (0, False, 0, ((2, 'O'), (1, 'N'), (1, 'S'))),

                (0, False, 0, ((2, 'O'), (1, 'Cl'), (1, 'Cl'))),
                (0, False, 0, ((2, 'O'), (1, 'Br'), (1, 'Br'))),

                (0, False, 0, ((2, 'O'), (1, 'S'), (1, 'S'))),

                (0, False, 0, ((2, 'O'), (1, 'C'), (1, 'C'))),
                (0, False, 0, ((2, 'O'), (1, 'C'), (1, 'F'))),
                (0, False, 0, ((2, 'O'), (1, 'C'), (1, 'Cl'))),
                (0, False, 0, ((2, 'O'), (1, 'C'), (1, 'Br'))),
                (0, False, 0, ((2, 'O'), (1, 'C'), (1, 'S'))),
                (0, False, 0, ((2, 'O'), (1, 'C'), (1, 'P'))),

                (0, False, 0, ((2, 'N'), (1, 'F'), (1, 'F'))),
                (0, False, 0, ((2, 'N'), (1, 'C'), (1, 'C'))),
                (0, False, 0, ((2, 'N'), (1, 'C'), (1, 'O'))),
                (0, False, 0, ((2, 'N'), (1, 'C'), (1, 'Cl'))),
                (0, False, 0, ((2, 'N'), (1, 'C'), (1, 'S'))),
                (0, False, 0, ((2, 'N'), (1, 'N'), (1, 'C'))),
                (0, False, 0, ((2, 'N'), (1, 'N'), (1, 'N'))),
                (0, False, 0, ((2, 'N'), (1, 'N'), (1, 'O'))),
                (0, False, 0, ((2, 'N'), (1, 'O'), (1, 'O'))),

                (0, False, 0, ((2, 'C'), (1, 'C'), (1, 'C'))),
                (0, False, 0, ((2, 'C'), (1, 'C'), (1, 'F'))),
                (0, False, 0, ((2, 'C'), (1, 'C'), (1, 'S'))),
                (0, False, 0, ((2, 'C'), (1, 'C'), (1, 'N'))),
                (0, False, 0, ((2, 'C'), (1, 'S'), (1, 'S'))),
                (0, False, 0, ((2, 'C'), (1, 'S'), (1, 'N'))),
                (0, False, 0, ((2, 'C'), (1, 'N'), (1, 'N'))),
                (0, False, 0, ((2, 'C'), (1, 'O'), (1, 'O'))),

                (0, False, 0, ((2, 'S'), (1, 'C'), (1, 'C'))),
                (0, False, 0, ((2, 'S'), (1, 'C'), (1, 'O'))),
                (0, False, 0, ((2, 'S'), (1, 'O'), (1, 'O'))),
                (0, False, 0, ((2, 'S'), (1, 'C'), (1, 'N'))),
                (0, False, 0, ((2, 'S'), (1, 'C'), (1, 'S'))),

                (0, False, 0, ((1, 'N'), (1, 'F'), (1, 'F'), (1, 'F'))),
                (0, False, 0, ((1, 'F'), (1, 'C'), (1, 'C'), (1, 'C'))),
                (0, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'), (1, 'C'))),
                (0, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'))),

                (0, False, 0, ((1, 'O'), (1, 'C'), (1, 'C'), (1, 'C'))),
                (0, False, 0, ((1, 'O'), (1, 'O'), (1, 'C'), (1, 'C'))),
                (0, False, 0, ((1, 'O'), (1, 'O'), (1, 'O'), (1, 'C'))),
                (0, False, 0, ((1, 'O'), (1, 'O'), (1, 'O'), (1, 'N'))),
                (0, False, 0, ((1, 'O'), (1, 'O'), (1, 'N'), (1, 'C'))),
                (0, False, 0, ((1, 'O'), (1, 'N'), (1, 'C'), (1, 'C'))),

                (0, False, 0, ((1, 'C'), (1, 'C'), (1, 'Cl'), (1, 'Cl'))),
                (0, False, 0, ((1, 'C'), (1, 'C'), (1, 'N'), (1, 'N'))),

                (0, False, 0, ((1, 'S'), (1, 'C'), (1, 'C'), (1, 'C'))),
                (0, False, 0, ((1, 'C'), (1, 'C'), (1, 'C'), (1, 'C'))),
                (0, False, 0, ((1, 'C'), (1, 'C'), (1, 'C'), (1, 'I'))),

                (0, False, 0, ((2, 'O'), (2, 'O'), (2, 'O'))),
                (0, False, 0, ((2, 'O'), (2, 'O'), (2, 'C'))),

                # sulfat derivatives
                (0, False, 0, ((2, 'O'), (2, 'O'), (1, 'O'), (1, 'O'))),
                (0, False, 0, ((2, 'O'), (2, 'O'), (1, 'O'), (1, 'N'))),
                (0, False, 0, ((2, 'O'), (2, 'O'), (1, 'O'), (1, 'C'))),
                (0, False, 0, ((2, 'O'), (2, 'O'), (1, 'O'), (1, 'S'))),
                (0, False, 0, ((2, 'O'), (2, 'O'), (1, 'O'), (1, 'F'))),
                (0, False, 0, ((2, 'O'), (2, 'O'), (1, 'O'), (1, 'Cl'))),
                (0, False, 0, ((2, 'O'), (2, 'O'), (1, 'O'), (1, 'Br'))),
                (0, False, 0, ((2, 'O'), (2, 'O'), (1, 'O'), (1, 'I'))),
                (0, False, 0, ((2, 'O'), (2, 'O'), (1, 'N'), (1, 'N'))),
                (0, False, 0, ((2, 'O'), (2, 'O'), (1, 'N'), (1, 'C'))),
                (0, False, 0, ((2, 'O'), (2, 'O'), (1, 'N'), (1, 'S'))),
                (0, False, 0, ((2, 'O'), (2, 'O'), (1, 'N'), (1, 'F'))),
                (0, False, 0, ((2, 'O'), (2, 'O'), (1, 'N'), (1, 'Cl'))),
                (0, False, 0, ((2, 'O'), (2, 'O'), (1, 'C'), (1, 'C'))),
                (0, False, 0, ((2, 'O'), (2, 'O'), (1, 'C'), (1, 'S'))),
                (0, False, 0, ((2, 'O'), (2, 'O'), (1, 'C'), (1, 'F'))),
                (0, False, 0, ((2, 'O'), (2, 'O'), (1, 'C'), (1, 'Cl'))),
                (0, False, 0, ((2, 'O'), (2, 'O'), (1, 'C'), (1, 'Br'))),
                (0, False, 0, ((2, 'O'), (2, 'O'), (1, 'C'), (1, 'I'))),
                (0, False, 0, ((2, 'O'), (2, 'O'), (1, 'F'), (1, 'F'))),
                (0, False, 0, ((2, 'O'), (2, 'O'), (1, 'Cl'), (1, 'Cl'))),
                (0, False, 0, ((2, 'O'), (2, 'O'), (1, 'Cl'), (1, 'F'))),

                (0, False, 0, ((2, 'O'), (2, 'N'), (1, 'N'), (1, 'O'))),
                (0, False, 0, ((2, 'O'), (2, 'N'), (1, 'N'), (1, 'N'))),
                (0, False, 0, ((2, 'O'), (2, 'N'), (1, 'N'), (1, 'C'))),
                (0, False, 0, ((2, 'O'), (2, 'N'), (1, 'O'), (1, 'O'))),
                (0, False, 0, ((2, 'O'), (2, 'N'), (1, 'O'), (1, 'C'))),
                (0, False, 0, ((2, 'O'), (2, 'N'), (1, 'C'), (1, 'C'))),
                (0, False, 0, ((2, 'O'), (2, 'N'), (1, 'C'), (1, 'Cl'))),
                (0, False, 0, ((2, 'O'), (2, 'N'), (1, 'C'), (1, 'F'))),

                (0, False, 0, ((2, 'N'), (2, 'N'), (1, 'C'), (1, 'C'))),

                # aci forms of tautomers
                (0, False, 0, ((2, 'O'), (2, 'C'), (1, 'C'), (1, 'C'))),
                (0, False, 0, ((2, 'O'), (2, 'C'), (1, 'O'), (1, 'C'))),
                (0, False, 0, ((2, 'O'), (2, 'C'), (1, 'O'), (1, 'N'))),
                (0, False, 0, ((2, 'O'), (2, 'C'), (1, 'O'), (1, 'O'))),
                (0, False, 0, ((2, 'O'), (2, 'C'), (1, 'N'), (1, 'C'))),
                (0, False, 0, ((2, 'O'), (2, 'C'), (1, 'N'), (1, 'N'))),

                (0, False, 0, ((2, 'N'), (2, 'C'), (1, 'C'), (1, 'O'))),

                (0, False, 0, ((2, 'O'), (2, 'S'), (1, 'O'), (1, 'O'))),  # [S2O3]2-
                (0, False, 0, ((2, 'O'), (2, 'S'), (1, 'O'), (1, 'C'))),
                (0, False, 0, ((2, 'O'), (2, 'S'), (1, 'C'), (1, 'C'))),

                (0, False, 0, ((2, 'S'), (2, 'S'), (1, 'O'), (1, 'O'))),

                (0, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'))),
                (0, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'), (1, 'O'))),
                (0, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'), (1, 'C'))))

    @property
    def atomic_radius(self):
        return .87

    @property
    def mdl_isotope(self):
        return 32

    @property
    def is_forming_single_bonds(self):
        return True

    @property
    def is_forming_double_bonds(self):
        return True


class Se(Element, PeriodIV, GroupXVI):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 34

    @property
    def isotopes_distribution(self):
        return {73: 0., 74: 0.0089, 75: 0., 76: 0.0937, 77: 0.0763, 78: 0.2377, 80: 0.4961, 82: 0.0873}

    @property
    def isotopes_masses(self):
        return {73: 72.926765, 74: 73.922477, 75: 74.922523, 76: 75.919214, 77: 76.919915, 78: 77.917310, 80: 79.916522,
                82: 81.916700}

    @property
    def _common_valences(self):
        return 2,

    @property
    def _valences_exceptions(self):
        return ((-1, False, 1, ()), (-2, False, 0, ()),
                (0, False, 0, ()),  # elemental
                (1, False, 0, ((1, 'C'), (1, 'C'), (1, 'C'))),
                (1, False, 0, ((1, 'C'), (2, 'C'))),

                (0, False, 0, ((2, 'O'), (2, 'O'))),
                (0, False, 0, ((2, 'S'), (2, 'S'))),
                (0, False, 0, ((2, 'N'), (2, 'N'))),

                (0, False, 0, ((2, 'O'), (1, 'O'), (1, 'O'))),
                (0, False, 0, ((2, 'O'), (1, 'Cl'), (1, 'Cl'))),
                (0, False, 0, ((2, 'O'), (1, 'F'), (1, 'F'))),
                (0, False, 0, ((2, 'O'), (1, 'C'), (1, 'C'))),
                (0, False, 0, ((2, 'O'), (1, 'N'), (1, 'C'))),
                (0, False, 0, ((2, 'O'), (1, 'O'), (1, 'C'))),

                (0, False, 0, ((2, 'C'), (1, 'C'), (1, 'C'))),
                (0, False, 0, ((2, 'C'), (1, 'O'), (1, 'O'))),

                (0, False, 0, ((1, 'C'), (1, 'C'), (1, 'Cl'), (1, 'Cl'))),
                (0, False, 0, ((1, 'C'), (1, 'C'), (1, 'Cl'), (1, 'O'))),
                (0, False, 0, ((1, 'C'), (1, 'C'), (1, 'Br'), (1, 'Br'))),
                (0, False, 0, ((1, 'C'), (1, 'C'), (1, 'Br'), (1, 'O'))),

                (0, False, 0, ((2, 'O'), (2, 'O'), (1, 'O'), (1, 'O'))),
                (0, False, 0, ((2, 'O'), (2, 'O'), (1, 'Cl'), (1, 'Cl'))),
                (0, False, 0, ((2, 'O'), (2, 'O'), (1, 'C'), (1, 'C'))),
                (0, False, 0, ((2, 'O'), (2, 'O'), (1, 'C'), (1, 'O'))),
                (0, False, 0, ((2, 'O'), (2, 'N'), (1, 'C'), (1, 'C'))),
                (0, False, 0, ((2, 'O'), (2, 'N'), (1, 'C'), (1, 'O'))),
                (0, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'))))

    @property
    def atomic_radius(self):
        return 1.03

    @property
    def mdl_isotope(self):
        return 79

    @property
    def is_forming_single_bonds(self):
        return True

    @property
    def is_forming_double_bonds(self):
        return True


class Te(Element, PeriodV, GroupXVI):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 52

    @property
    def isotopes_distribution(self):
        return {120: 0.0009, 122: 0.0255, 123: 0.0089, 124: 0.0474, 125: 0.0707, 126: 0.1884, 128: 0.3174, 130: 0.3408}

    @property
    def isotopes_masses(self):
        return {120: 119.90402, 122: 121.903047, 123: 122.904273, 124: 123.90282, 125: 124.904425, 126: 125.903306,
                128: 127.904461, 130: 129.906223}

    @property
    def _common_valences(self):
        return 2,

    @property
    def _valences_exceptions(self):
        return ((0, False, 0, ()),  # elemental,
                (1, False, 0, ((1, 'C'), (1, 'C'), (1, 'C'))),
                (-1, False, 0, ((1, 'C'), (1, 'O'), (1, 'Cl'), (1, 'Cl'), (1, 'Cl'))),
                (1, False, 0, ((1, 'C'), (2, 'C'))),
                (0, False, 0, ((2, 'O'), (2, 'O'))),
                (0, False, 0, ((2, 'O'), (1, 'O'), (1, 'O'))),
                (0, False, 0, ((2, 'O'), (1, 'O'), (1, 'C'))),
                (0, False, 0, ((2, 'O'), (1, 'C'), (1, 'C'))),
                (0, False, 0, ((1, 'C'), (1, 'Cl'), (1, 'Cl'), (1, 'Cl'))),
                (0, False, 0, ((1, 'C'), (1, 'C'), (1, 'Cl'), (1, 'Cl'))),
                (0, False, 0, ((1, 'C'), (1, 'C'), (1, 'Br'), (1, 'Br'))),
                (0, False, 0, ((1, 'C'), (1, 'C'), (1, 'O'), (1, 'Cl'))),
                (0, False, 0, ((1, 'C'), (1, 'C'), (1, 'O'), (1, 'O'))),
                (0, False, 0, ((1, 'C'), (1, 'C'), (1, 'C'), (1, 'C'))),
                (0, False, 0, ((1, 'C'), (1, 'O'), (1, 'Cl'), (1, 'Cl'))),

                (0, False, 0, ((1, 'O'), (1, 'O'), (1, 'Cl'), (1, 'Cl'))),

                (0, False, 0, ((2, 'O'), (2, 'O'), (1, 'O'), (1, 'O'))),
                (0, False, 0, ((1, 'O'), (1, 'O'), (1, 'O'), (1, 'O'), (1, 'O'), (1, 'O'))),
                (0, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'))),
                (0, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'), (1, 'O'))))

    @property
    def atomic_radius(self):
        return 1.23

    @property
    def mdl_isotope(self):
        return 128

    @property
    def is_forming_single_bonds(self):
        return True

    @property
    def is_forming_double_bonds(self):
        return True


class Po(Element, PeriodVI, GroupXVI):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 84

    @property
    def isotopes_distribution(self):
        return {210: 1.0}

    @property
    def isotopes_masses(self):
        return {210: 209.982874}

    @property
    def _common_valences(self):
        return 0, 2

    @property
    def _valences_exceptions(self):
        return ((0, False, 0, ((2, 'O'), (2, 'O'), (2, 'O'))),
                (0, False, 0, ((2, 'O'), (2, 'O'))),
                (0, False, 0, ((1, 'Cl'), (1, 'Cl'), (1, 'Cl'), (1, 'Cl'))),
                (0, False, 0, ((1, 'Br'), (1, 'Br'), (1, 'Br'), (1, 'Br'))),
                (0, False, 0, ((1, 'I'), (1, 'I'), (1, 'I'), (1, 'I'))),
                (0, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'))))

    @property
    def atomic_radius(self):
        return 1.35

    @property
    def mdl_isotope(self):
        return 209


class Lv(Element, PeriodVII, GroupXVI):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 116

    @property
    def isotopes_distribution(self):
        return {293: 1.0}

    @property
    def isotopes_masses(self):
        return {293: 293.204555}

    @property
    def _common_valences(self):
        return 0,

    @property
    def _valences_exceptions(self):
        return ()

    @property
    def atomic_radius(self):
        return 1.35  # unknown, taken radius of previous element in group

    @property
    def mdl_isotope(self):
        return 293


__all__ = ['O', 'S', 'Se', 'Te', 'Po', 'Lv']
