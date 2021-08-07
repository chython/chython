# -*- coding: utf-8 -*-
#
#  Copyright 2019-2021 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from .element import Element
from .groups import GroupXV
from .periods import PeriodII, PeriodIII, PeriodIV, PeriodV, PeriodVI, PeriodVII


class N(Element, PeriodII, GroupXV):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 7

    @property
    def isotopes_distribution(self):
        return {14: 0.99632, 15: 0.00368}

    @property
    def isotopes_masses(self):
        return {14: 14.003074, 15: 15.000109}

    @property
    def _common_valences(self):
        return 3,

    @property
    def _valences_exceptions(self):
        return (-1, False, 2, ()), (1, False, 4, ()), (0, True, 0, ((2, 'O'),))  # *NO

    @property
    def atomic_radius(self):
        return .56


class P(Element, PeriodIII, GroupXV):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 15

    @property
    def isotopes_distribution(self):
        return {31: 1.0, 32: 0., 33: 0.}

    @property
    def isotopes_masses(self):
        return {31: 30.973762, 32: 31.973908, 33: 32.971726}

    @property
    def _common_valences(self):
        return 3, 5

    @property
    def _valences_exceptions(self):
        return ((-1, False, 2, ()), (1, False, 4, (),),
                (-1, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'))))

    @property
    def atomic_radius(self):
        return .98


class As(Element, PeriodIV, GroupXV):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 33

    @property
    def isotopes_distribution(self):
        return {75: 1.0, 76: 0., 77: 0.}

    @property
    def isotopes_masses(self):
        return {75: 74.921596, 76: 75.922394, 77: 76.920647}

    @property
    def _common_valences(self):
        return 0, 3, 5

    @property
    def _valences_exceptions(self):
        return (1, False, 4, ()),

    @property
    def atomic_radius(self):
        return 1.14


class Sb(Element, PeriodV, GroupXV):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 51

    @property
    def isotopes_distribution(self):
        return {121: 0.5721, 123: 0.4279}

    @property
    def isotopes_masses(self):
        return {121: 120.903818, 123: 122.904216}

    @property
    def _common_valences(self):
        return 0, 3, 5

    @property
    def _valences_exceptions(self):
        return (1, False, 4, ()),

    @property
    def atomic_radius(self):
        return 1.33


class Bi(Element, PeriodVI, GroupXV):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 83

    @property
    def isotopes_distribution(self):
        return {207: 0., 209: 1.0, 210: 0.}

    @property
    def isotopes_masses(self):
        return {207: 206.978471, 209: 208.980383, 210: 209.984120}

    @property
    def _common_valences(self):
        return 0, 3

    @property
    def _valences_exceptions(self):
        return ((0, False, 0, ((1, 'Cl'),)),
                (0, False, 0, ((1, 'Br'),)),

                (0, False, 0, ((1, 'Cl'), (1, 'Cl'))),
                (0, False, 0, ((1, 'Br'), (1, 'Br'))),
                (0, False, 0, ((1, 'I'), (1, 'I'))),
                (0, False, 0, ((1, 'S'), (1, 'S'))),
                (0, False, 0, ((2, 'S'),)),
                (0, False, 0, ((1, 'Se'), (1, 'Se'))),
                (0, False, 0, ((2, 'Se'),)),

                (0, False, 0, ((1, 'Cl'), (1, 'Cl'), (1, 'Cl'), (1, 'Cl'))),
                (0, False, 0, ((2, 'O'), (2, 'O'))),

                (0, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'))),
                (0, False, 0, ((1, 'O'), (2, 'O'), (2, 'O'))),
                (0, False, 0, ((1, 'O'), (1, 'O'), (1, 'O'), (2, 'O'))))

    @property
    def atomic_radius(self):
        return 1.43


class Mc(Element, PeriodVII, GroupXV):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 115

    @property
    def isotopes_distribution(self):
        return {289: 1.0}

    @property
    def isotopes_masses(self):
        return {289: 289.0}

    @property
    def _common_valences(self):
        return 0,

    @property
    def _valences_exceptions(self):
        return ()

    @property
    def atomic_radius(self):
        return 1.43  # unknown, taken radius of previous element in group


__all__ = ['N', 'P', 'As', 'Sb', 'Bi', 'Mc']
