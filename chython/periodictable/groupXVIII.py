# -*- coding: utf-8 -*-
#
#  Copyright 2019-2024 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from .base.groups import GroupXVIII
from .base.periods import *


class He(Element, PeriodI, GroupXVIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 2

    @property
    def isotopes_distribution(self):
        return {3: 1e-06, 4: 0.999999}

    @property
    def isotopes_masses(self):
        return {3: 3.016029, 4: 4.002603}

    @property
    def _common_valences(self):
        return 0,

    @property
    def _valences_exceptions(self):
        return ()

    @property
    def atomic_radius(self):
        return .31

    @property
    def mdl_isotope(self):
        return 4


class Ne(Element, PeriodII, GroupXVIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 10

    @property
    def isotopes_distribution(self):
        return {20: 0.9048, 21: 0.0027, 22: 0.0925}

    @property
    def isotopes_masses(self):
        return {20: 19.99244, 21: 20.993847, 22: 21.991386}

    @property
    def _common_valences(self):
        return 0,

    @property
    def _valences_exceptions(self):
        return ()

    @property
    def atomic_radius(self):
        return .38

    @property
    def mdl_isotope(self):
        return 20


class Ar(Element, PeriodIII, GroupXVIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 18

    @property
    def isotopes_distribution(self):
        return {36: 0.003365, 38: 0.000632, 40: 0.996003}

    @property
    def isotopes_masses(self):
        return {36: 35.967546, 38: 37.962732, 40: 39.962383}

    @property
    def _common_valences(self):
        return 0,

    @property
    def _valences_exceptions(self):
        return ()

    @property
    def atomic_radius(self):
        return .71

    @property
    def mdl_isotope(self):
        return 40


class Kr(Element, PeriodIV, GroupXVIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 36

    @property
    def isotopes_distribution(self):
        return {78: 0.0035, 80: 0.0228, 81: 0., 82: 0.1158, 83: 0.1149, 84: 0.57, 86: 0.173}

    @property
    def isotopes_masses(self):
        return {78: 77.920386, 80: 79.916378, 81: 80.916592, 82: 81.913485, 83: 82.914136, 84: 83.911507, 86: 85.91061}

    @property
    def _common_valences(self):
        return 0,

    @property
    def _valences_exceptions(self):
        return ()

    @property
    def atomic_radius(self):
        return .87

    @property
    def mdl_isotope(self):
        return 84


class Xe(Element, PeriodV, GroupXVIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 54

    @property
    def isotopes_distribution(self):
        return {124: 0.0009, 126: 0.0009, 127: 0., 128: 0.0192, 129: 0.2644, 130: 0.0408, 131: 0.2118, 132: 0.2689,
                133: 0., 134: 0.1044, 136: 0.0887}

    @property
    def isotopes_masses(self):
        return {124: 123.905896, 126: 125.904269, 127: 126.905184, 128: 127.90353, 129: 128.904779, 130: 129.903508,
                131: 130.905082, 132: 131.904155, 133: 132.905911, 134: 133.905394, 136: 135.90722}

    @property
    def _common_valences(self):
        return 0,

    @property
    def _valences_exceptions(self):
        # XeF2, XeF4, XeF6, XeO3, XeO4, XeO2F2, XeOF4, XeO3F2, [XeO6]4-
        return ((0, False, 0, ((1, 'F'), (1, 'F'))),
                (0, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'))),
                (0, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'))),
                (0, False, 0, ((2, 'O'), (2, 'O'), (2, 'O'))),
                (0, False, 0, ((2, 'O'), (2, 'O'), (2, 'O'), (2, 'O'))),
                (0, False, 0, ((2, 'O'), (2, 'O'), (1, 'F'), (1, 'F'))),
                (0, False, 0, ((2, 'O'), (1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'))),
                (0, False, 0, ((2, 'O'), (2, 'O'), (2, 'O'), (1, 'F'), (1, 'F'))),
                (0, False, 0, ((1, 'O'), (1, 'O'), (1, 'O'), (1, 'O'), (2, 'O'), (2, 'O'))))

    @property
    def atomic_radius(self):
        return 1.08

    @property
    def mdl_isotope(self):
        return 131


class Rn(Element, PeriodVI, GroupXVIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 86

    @property
    def isotopes_distribution(self):
        return {222: 1.0}

    @property
    def isotopes_masses(self):
        return {222: 222.017578}

    @property
    def _common_valences(self):
        return 0,

    @property
    def _valences_exceptions(self):
        return (0, False, 0, ((1, 'F'), (1, 'F'))), (1, False, 0, ((1, 'F'),))

    @property
    def atomic_radius(self):
        return 1.2

    @property
    def mdl_isotope(self):
        return 222


class Og(Element, PeriodVII, GroupXVIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 118

    @property
    def isotopes_distribution(self):
        return {294: 1.0}

    @property
    def isotopes_masses(self):
        return {294: 294.0}

    @property
    def _common_valences(self):
        return 0,

    @property
    def _valences_exceptions(self):
        return ()

    @property
    def atomic_radius(self):
        return 1.2  # unknown, taken radius of previous element in group

    @property
    def mdl_isotope(self):
        return 294


__all__ = ['He', 'Ne', 'Ar', 'Kr', 'Xe', 'Rn', 'Og']
