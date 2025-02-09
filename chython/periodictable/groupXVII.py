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
from .base.groups import GroupXVII
from .base.periods import PeriodII, PeriodIII, PeriodIV, PeriodV, PeriodVI, PeriodVII


class F(Element, PeriodII, GroupXVII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 9

    @property
    def isotopes_distribution(self):
        return {17: 0., 18: 0., 19: 1.0}

    @property
    def isotopes_masses(self):
        return {17: 17.002095, 18: 18.000938, 19: 18.998403}

    @property
    def _common_valences(self):
        return 1,

    @property
    def _valences_exceptions(self):
        return (-1, False, 0, ()),

    @property
    def atomic_radius(self):
        return .42

    @property
    def mdl_isotope(self):
        return 19

    @property
    def is_forming_single_bonds(self):
        return True


class Cl(Element, PeriodIII, GroupXVII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 17

    @property
    def isotopes_distribution(self):
        return {35: 0.7578, 36: 0., 37: 0.2422}

    @property
    def isotopes_masses(self):
        return {35: 34.968853, 36: 35.968307, 37: 36.965903}

    @property
    def _common_valences(self):
        return 1,

    @property
    def _valences_exceptions(self):
        return ((-1, False, 0, ()),
                (-1, False, 0, ((1, 'Cl'), (1, 'I'))),  # [I-Cl-Cl]-

                (0, False, 0, ((1, 'O'), (2, 'O'))),  # HClO2
                (0, False, 0, ((1, 'O'), (2, 'O'), (2, 'O'))),  # HClO3
                (0, False, 0, ((1, 'O'), (2, 'O'), (2, 'O'), (2, 'O'))),  # HClO4

                (0, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'))),  # ClF3

                (0, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'))),  # ClF5
                (0, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'), (2, 'O'))),  # ClOF3
                (0, False, 0, ((1, 'F'), (2, 'O'), (2, 'O'))))  # ClO2F

    @property
    def atomic_radius(self):
        return .79

    @property
    def mdl_isotope(self):
        return 35

    @property
    def is_forming_single_bonds(self):
        return True

    @property
    def is_forming_double_bonds(self):
        return True


class Br(Element, PeriodIV, GroupXVII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 35

    @property
    def isotopes_distribution(self):
        return {76: 0., 77: 0., 79: 0.5069, 81: 0.4931, 82: 0.}

    @property
    def isotopes_masses(self):
        return {76: 75.924541, 77: 76.921379, 79: 78.918338, 81: 80.916291, 82: 81.916804}

    @property
    def _common_valences(self):
        return 1,

    @property
    def _valences_exceptions(self):
        return ((-1, False, 0, ()),
                (-1, False, 0, ((1, 'Br'), (1, 'I'))),  # [I-Br-Br]-
                (-1, False, 0, ((1, 'Br'), (1, 'Br'))),  # [Br-Br-Br]-
                (-1, False, 0, ((1, 'Br'), (1, 'Cl'))),  # [Br-Br-Cl]-
                (-1, False, 0, ((1, 'Cl'), (1, 'Cl'))),  # [Cl-Br-Cl]-
                (-1, False, 0, ((1, 'I'), (1, 'I'))),  # [I-Br-I]-

                (0, False, 0, ((1, 'O'), (2, 'O'))),  # HBrO2
                (0, False, 0, ((1, 'O'), (2, 'O'), (2, 'O'))),  # HBrO3
                (0, False, 0, ((1, 'O'), (2, 'O'), (2, 'O'), (2, 'O'))),  # HBrO4

                (0, False, 0, ((1, 'O'), (1, 'O'), (1, 'O'))),  # Br(OX)3
                (0, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'))),  # BrF3

                (0, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'))),  # BrF5
                (0, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'), (2, 'O'))),  # BrOF3
                (0, False, 0, ((1, 'F'), (2, 'O'), (2, 'O'))),  # BrO2F

                (0, False, 0, ((1, 'F'), (2, 'O'), (2, 'O'), (2, 'O'))))  # BrO3F

    @property
    def atomic_radius(self):
        return 0.94

    @property
    def mdl_isotope(self):
        return 80

    @property
    def is_forming_single_bonds(self):
        return True

    @property
    def is_forming_double_bonds(self):
        return True


class I(Element, PeriodV, GroupXVII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 53

    @property
    def isotopes_distribution(self):
        return {123: 0., 124: 0., 125: 0., 127: 1.0, 129: 0., 131: 0., 135: 0.}

    @property
    def isotopes_masses(self):
        return {123: 122.905589, 124: 123.906210, 125: 124.904630, 127: 126.904468, 129: 128.904988, 131: 130.906125,
                135: 134.910048}

    @property
    def _common_valences(self):
        return 1,

    @property
    def _valences_exceptions(self):
        return ((-1, False, 0, ()),
                (-1, False, 0, ((1, 'I'), (1, 'I'))),  # [I-I-I]-
                (-1, False, 0, ((1, 'I'), (1, 'Br'))),  # [I-I-Br]-
                (-1, False, 0, ((1, 'Cl'), (1, 'Cl'))),  # [Cl-I-Cl]-
                (1, False, 0, ((1, 'C'), (1, 'C'))),

                (0, False, 0, ((1, 'O'), (2, 'O'))),  # HIO2
                (0, False, 0, ((1, 'C'), (2, 'O'))),
                (0, False, 0, ((1, 'C'), (2, 'C'))),
                (0, False, 0, ((1, 'C'), (2, 'N'))),
                (0, False, 0, ((1, 'O'), (2, 'O'), (2, 'O'))),  # HIO3
                (0, False, 0, ((1, 'O'), (2, 'O'), (2, 'O'), (2, 'O'))),  # HIO4
                (0, False, 0, ((1, 'O'), (1, 'O'), (1, 'O'), (2, 'O'), (2, 'O'))),  # H3IO5
                (0, False, 0, ((1, 'O'), (1, 'O'), (1, 'O'), (1, 'O'), (1, 'O'), (2, 'O'))),  # H5IO6

                (0, False, 0, ((1, 'O'), (1, 'O'), (1, 'O'))),  # I(OX)3
                (0, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'))),  # IHal3
                (0, False, 0, ((1, 'Cl'), (1, 'Cl'), (1, 'Cl'))),
                (0, False, 0, ((1, 'Br'), (1, 'Br'), (1, 'Br'))),
                (0, False, 0, ((1, 'C'), (1, 'O'), (1, 'Cl'))),
                (0, False, 0, ((1, 'C'), (1, 'O'), (1, 'C'))),
                (0, False, 0, ((1, 'C'), (1, 'O'), (1, 'O'))),
                (0, False, 0, ((1, 'C'), (1, 'O'), (1, 'N'))),
                (0, False, 0, ((1, 'C'), (1, 'F'), (1, 'F'))),
                (0, False, 0, ((1, 'C'), (1, 'Cl'), (1, 'Cl'))),
                (0, False, 0, ((1, 'C'), (1, 'C'), (1, 'Cl'))),
                (0, False, 0, ((1, 'C'), (1, 'C'), (1, 'C'))),
                (0, False, 0, ((1, 'C'), (1, 'C'), (1, 'N'))),

                (0, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'))),  # IF5
                (0, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'), (2, 'O'))),  # IOF3
                (0, False, 0, ((1, 'F'), (2, 'O'), (2, 'O'))),  # IO2F
                (0, False, 0, ((1, 'C'), (2, 'O'), (2, 'O'))),
                (0, False, 0, ((1, 'C'), (1, 'O'), (1, 'O'), (2, 'O'))),
                (0, False, 0, ((1, 'C'), (1, 'O'), (1, 'O'), (1, 'O'), (1, 'O'))),

                (0, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'))),  # IF7
                (0, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'), (2, 'O'))),  # IOF5
                (0, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'), (2, 'O'), (2, 'O'))),  # IO2F3
                (0, False, 0, ((1, 'F'), (2, 'O'), (2, 'O'), (2, 'O'))))  # IO3F

    @property
    def atomic_radius(self):
        return 1.15

    @property
    def mdl_isotope(self):
        return 127

    @property
    def is_forming_single_bonds(self):
        return True

    @property
    def is_forming_double_bonds(self):
        return True


class At(Element, PeriodVI, GroupXVII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 85

    @property
    def isotopes_distribution(self):
        return {210: 1.0, 211: 0.}

    @property
    def isotopes_masses(self):
        return {210: 209.987155, 211: 210.987496}

    @property
    def _common_valences(self):
        return 0, 1

    @property
    def _valences_exceptions(self):
        return ((1, False, 0, ()), (-1, False, 0, ()),
                (0, False, 0, ((1, 'O'), (2, 'O'), (2, 'O'))))

    @property
    def atomic_radius(self):
        return 1.27

    @property
    def mdl_isotope(self):
        return 210

    @property
    def is_forming_single_bonds(self):
        return True

    @property
    def is_forming_double_bonds(self):
        return True


class Ts(Element, PeriodVII, GroupXVII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 117

    @property
    def isotopes_distribution(self):
        return {293: 1.0}

    @property
    def isotopes_masses(self):
        return {293: 293.0}

    @property
    def _common_valences(self):
        return 0,

    @property
    def _valences_exceptions(self):
        return ()

    @property
    def atomic_radius(self):
        return 1.27  # unknown, taken radius of previous element in group

    @property
    def mdl_isotope(self):
        return 297


__all__ = ['F', 'Cl', 'Br', 'I', 'At', 'Ts']
