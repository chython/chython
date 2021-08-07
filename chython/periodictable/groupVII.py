# -*- coding: utf-8 -*-
#
#  Copyright 2019-2021 Ramil Nugmanov <nougmanoff@protonmail.com>
#  Copyright 2019 Tagir Akhmetshin <tagirshin@gmail.com>
#  Copyright 2019 Alexander Nikanshin <17071996sasha@gmail.com>
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
from .groups import GroupVII
from .periods import PeriodIV, PeriodV, PeriodVI, PeriodVII


class Mn(Element, PeriodIV, GroupVII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 25

    @property
    def isotopes_distribution(self):
        return {52: 0., 55: 1.0}

    @property
    def isotopes_masses(self):
        return {52: 51.945566, 55: 54.938050}

    @property
    def _common_valences(self):
        return 0, 2

    @property
    def _valences_exceptions(self):
        return ((2, False, 0, ()), (3, False, 0, ()),
                (0, False, 0, ((2, 'O'), (2, 'O'))),
                (0, False, 0, ((2, 'O'),)),  # MnO
                (0, False, 0, ((2, 'O'), (1, 'O'))),  # Mn2O3
                (0, False, 0, ((1, 'O'), (1, 'O'), (1, 'O'))),
                (0, False, 0, ((2, 'O'), (2, 'O'))),  # MnO2
                (0, False, 0, ((2, 'O'), (2, 'O'), (1, 'O'), (1, 'O'))),  # [MnO4]2-
                (0, False, 0, ((2, 'O'), (2, 'O'), (2, 'O'), (1, 'O'))))  # [MnO4]-

    @property
    def atomic_radius(self):
        return 1.61


class Tc(Element, PeriodV, GroupVII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 43

    @property
    def isotopes_distribution(self):
        return {98: 0., 99: 1.0}

    @property
    def isotopes_masses(self):
        return {98: 97.907216, 99: 98.906255}

    @property
    def _common_valences(self):
        return 0,

    @property
    def _valences_exceptions(self):
        return ((0, False, 0, ((2, 'O'), (2, 'O'))),  # TcO2
                (0, False, 0, ((2, 'S'), (2, 'S'))))  # TcS2

    @property
    def atomic_radius(self):
        return 1.83


class Re(Element, PeriodVI, GroupVII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 75

    @property
    def isotopes_distribution(self):
        return {185: 0.374, 186: 0., 187: 0.626, 188: 0.}

    @property
    def isotopes_masses(self):
        return {185: 184.952956, 186: 185.954986, 187: 186.955751, 188: 187.958114}

    @property
    def _common_valences(self):
        return 0,

    @property
    def _valences_exceptions(self):
        return ()

    @property
    def atomic_radius(self):
        return 1.88


class Bh(Element, PeriodVII, GroupVII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 107

    @property
    def isotopes_distribution(self):
        return {270: 1.0}

    @property
    def isotopes_masses(self):
        return {270: 270.133363}

    @property
    def _common_valences(self):
        return 0,

    @property
    def _valences_exceptions(self):
        return ()

    @property
    def atomic_radius(self):
        return 1.88  # unknown, taken radius of previous element in group


__all__ = ['Mn', 'Tc', 'Re', 'Bh']
