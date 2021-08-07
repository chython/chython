# -*- coding: utf-8 -*-
#
#  Copyright 2019-2021 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from .element import Element
from .groups import GroupXIII
from .periods import PeriodII, PeriodIII, PeriodIV, PeriodV, PeriodVI, PeriodVII


class B(Element, PeriodII, GroupXIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 5

    @property
    def isotopes_distribution(self):
        return {10: 0.199, 11: 0.801}

    @property
    def isotopes_masses(self):
        return {10: 10.012937, 11: 11.009305}

    @property
    def _common_valences(self):
        return 3,

    @property
    def _valences_exceptions(self):
        return (-1, False, 4, ()),

    @property
    def atomic_radius(self):
        return .87


class Al(Element, PeriodIII, GroupXIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 13

    @property
    def isotopes_distribution(self):
        return {27: 1.0}

    @property
    def isotopes_masses(self):
        return {27: 26.981538}

    @property
    def _common_valences(self):
        return 0,  # 3 omitted. added to exceptions

    @property
    def _valences_exceptions(self):
        return ((3, False, 0, ()), (2, False, 1, ()), (1, False, 2, ()), (0, False, 3, ()), (-1, False, 4, ()),
                (-3, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'))))

    @property
    def atomic_radius(self):
        return 1.18


class Ga(Element, PeriodIV, GroupXIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 31

    @property
    def isotopes_distribution(self):
        return {67: 0., 69: 0.60108, 71: 0.39892}

    @property
    def isotopes_masses(self):
        return {67: 66.928202, 69: 68.925581, 71: 70.924705}

    @property
    def _common_valences(self):
        return 0, 3

    @property
    def _valences_exceptions(self):
        return ((3, False, 0, ()),
                (0, False, 0, ((1, 'Cl'),)), (0, False, 0, ((1, 'Br'),)), (0, False, 0, ((1, 'I'),)),
                (-1, False, 0, ((1, 'H'), (1, 'H'), (1, 'H'), (1, 'H'))),
                (-1, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'))),
                (-1, False, 0, ((1, 'Cl'), (1, 'Cl'), (1, 'Cl'), (1, 'Cl'))),
                (-1, False, 0, ((1, 'Br'), (1, 'Br'), (1, 'Br'), (1, 'Br'))),
                (-1, False, 0, ((1, 'I'), (1, 'I'), (1, 'I'), (1, 'I'))))

    @property
    def atomic_radius(self):
        return 1.36


class In(Element, PeriodV, GroupXIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 49

    @property
    def isotopes_distribution(self):
        return {111: 0., 113: 0.0429, 115: 0.9571}

    @property
    def isotopes_masses(self):
        return {111: 110.905103, 113: 112.904061, 115: 114.903878}

    @property
    def _common_valences(self):
        return 0, 3

    @property
    def _valences_exceptions(self):
        return ((3, False, 0, ()),
                (0, False, 0, ((1, 'Cl'),)), (0, False, 0, ((1, 'Br'),)), (0, False, 0, ((1, 'I'),)),
                (0, False, 0, ((1, 'O'),)))

    @property
    def atomic_radius(self):
        return 1.56


class Tl(Element, PeriodVI, GroupXIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 81

    @property
    def isotopes_distribution(self):
        return {203: 0.29524, 205: 0.70476}

    @property
    def isotopes_masses(self):
        return {203: 202.972329, 205: 204.974412}

    @property
    def _common_valences(self):
        return 0, 1

    @property
    def _valences_exceptions(self):
        return ((1, False, 0, ()), (3, False, 0, ()),
                (-3, False, 0, ((1, 'Cl'), (1, 'Cl'), (1, 'Cl'), (1, 'Cl'), (1, 'Cl'), (1, 'Cl'))),
                (1, False, 0, ((1, 'C'), (1, 'C'))))

    @property
    def atomic_radius(self):
        return 1.56


class Nh(Element, PeriodVII, GroupXIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 113

    @property
    def isotopes_distribution(self):
        return {286: 1.0}

    @property
    def isotopes_masses(self):
        return {286: 286.182555}

    @property
    def _common_valences(self):
        return 0,

    @property
    def _valences_exceptions(self):
        return ()

    @property
    def atomic_radius(self):
        return 1.56  # unknown, taken radius of previous element in group


__all__ = ['B', 'Al', 'Ga', 'In', 'Tl', 'Nh']
