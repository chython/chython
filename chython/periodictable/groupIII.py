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
from .base.groups import GroupIII
from .base.periods import PeriodIV, PeriodV, PeriodVI, PeriodVII


class Sc(Element, PeriodIV, GroupIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 21

    @property
    def isotopes_distribution(self):
        return {44: 0., 45: 1.0}

    @property
    def isotopes_masses(self):
        return {44: 43.959403, 45: 44.955910}

    @property
    def _common_valences(self):
        return 0, 3

    @property
    def _valences_exceptions(self):
        return (3, False, 0, ()), (-3, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'), (1, 'F')))

    @property
    def atomic_radius(self):
        return 1.84

    @property
    def mdl_isotope(self):
        return 45


class Y(Element, PeriodV, GroupIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 39

    @property
    def isotopes_distribution(self):
        return {86: 0., 89: 1.0, 90: 0.}

    @property
    def isotopes_masses(self):
        return {86: 85.914886, 89: 88.905848, 90: 89.907152}

    @property
    def _common_valences(self):
        return 0, 3

    @property
    def _valences_exceptions(self):
        return (3, False, 0, ()),

    @property
    def atomic_radius(self):
        return 2.12

    @property
    def mdl_isotope(self):
        return 89


class La(Element, PeriodVI, GroupIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 57

    @property
    def isotopes_distribution(self):
        return {138: 0.0009, 139: 0.9991}

    @property
    def isotopes_masses(self):
        return {138: 137.907107, 139: 138.906348}

    @property
    def _common_valences(self):
        return 0, 3

    @property
    def _valences_exceptions(self):
        return (3, False, 0, ()),

    @property
    def atomic_radius(self):
        return 2.12  # unknown, taken radius of previous element in group

    @property
    def mdl_isotope(self):
        return 139


class Ce(Element, PeriodVI, GroupIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 58

    @property
    def isotopes_distribution(self):
        return {136: 0.00185, 138: 0.00251, 140: 0.8845, 142: 0.11114}

    @property
    def isotopes_masses(self):
        return {136: 135.90714, 138: 137.905986, 140: 139.905434, 142: 141.90924}

    @property
    def _common_valences(self):
        return 0, 3

    @property
    def _valences_exceptions(self):
        return ((3, False, 0, ()),
                (0, False, 0, ((1, 'O'), (1, 'O'), (1, 'O'), (1, 'O'))),
                (0, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'))),
                (0, False, 0, ((1, 'Cl'), (1, 'Cl'), (1, 'Cl'), (1, 'Cl'))),
                (0, False, 0, ((2, 'O'), (2, 'O'))))

    @property
    def atomic_radius(self):
        return 2.12  # unknown, taken radius of previous element in group

    @property
    def mdl_isotope(self):
        return 140


class Pr(Element, PeriodVI, GroupIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 59

    @property
    def isotopes_distribution(self):
        return {141: 1.0}

    @property
    def isotopes_masses(self):
        return {141: 140.907648}

    @property
    def _common_valences(self):
        return 0, 3

    @property
    def _valences_exceptions(self):
        return ((3, False, 0, ()),
                (0, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'))),
                (0, False, 0, ((2, 'O'), (2, 'O'))))

    @property
    def atomic_radius(self):
        return 2.47

    @property
    def mdl_isotope(self):
        return 141


class Nd(Element, PeriodVI, GroupIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 60

    @property
    def isotopes_distribution(self):
        return {142: 0.272, 143: 0.122, 144: 0.238, 145: 0.083, 146: 0.172, 148: 0.057, 150: 0.056}

    @property
    def isotopes_masses(self):
        return {142: 141.907719, 143: 142.90981, 144: 143.910083, 145: 144.912569, 146: 145.913112, 148: 147.916889,
                150: 149.920887}

    @property
    def _common_valences(self):
        return 0, 3

    @property
    def _valences_exceptions(self):
        return ((3, False, 0, ()),
                (0, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'))),
                (0, False, 0, ((2, 'O'),)),
                (0, False, 0, ((1, 'O'), (1, 'O'))),
                (0, False, 0, ((1, 'F'), (1, 'F'))),
                (0, False, 0, ((1, 'Cl'), (1, 'Cl'))),
                (0, False, 0, ((1, 'Br'), (1, 'Br'))),
                (0, False, 0, ((1, 'I'), (1, 'I'))),
                (0, False, 0, ((1, 'F'), (1, 'Cl'))),
                (0, False, 0, ((1, 'F'), (1, 'Br'))),
                (0, False, 0, ((1, 'F'), (1, 'I'))),
                (0, False, 0, ((1, 'C'), (1, 'C'))),
                (0, False, 0, ((1, 'H'), (1, 'H'))))

    @property
    def atomic_radius(self):
        return 2.06

    @property
    def mdl_isotope(self):
        return 144


class Pm(Element, PeriodVI, GroupIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 61

    @property
    def isotopes_distribution(self):
        return {145: 1.0}

    @property
    def isotopes_masses(self):
        return {145: 144.912749}

    @property
    def _common_valences(self):
        return 0, 3

    @property
    def _valences_exceptions(self):
        return (3, False, 0, ()),

    @property
    def atomic_radius(self):
        return 2.05

    @property
    def mdl_isotope(self):
        return 145


class Sm(Element, PeriodVI, GroupIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 62

    @property
    def isotopes_distribution(self):
        return {144: 0.0307, 145: 0., 147: 0.1499, 148: 0.1124, 149: 0.1382, 150: 0.0738, 152: 0.2675, 153: 0.,
                154: 0.2275}

    @property
    def isotopes_masses(self):
        return {144: 143.911995, 145: 144.913410, 147: 146.914893, 148: 147.914818, 149: 148.917180, 150: 149.917271,
                152: 151.919728, 153: 152.922097, 154: 153.922205}

    @property
    def _common_valences(self):
        return 0, 3

    @property
    def _valences_exceptions(self):
        return ((3, False, 0, ()),
                (0, False, 0, ((1, 'O'), (1, 'O'))),
                (0, False, 0, ((1, 'F'), (1, 'F'))),
                (0, False, 0, ((1, 'Cl'), (1, 'Cl'))),
                (0, False, 0, ((1, 'Br'), (1, 'Br'))),
                (0, False, 0, ((1, 'I'), (1, 'I'))),
                (0, False, 0, ((1, 'F'), (1, 'Cl'))),
                (0, False, 0, ((1, 'F'), (1, 'Br'))),
                (0, False, 0, ((1, 'F'), (1, 'I'))),
                (0, False, 0, ((1, 'C'), (1, 'C'))),
                (0, False, 0, ((1, 'H'), (1, 'H'))),
                (0, False, 0, ((2, 'O'),)))

    @property
    def atomic_radius(self):
        return 2.38

    @property
    def mdl_isotope(self):
        return 150


class Eu(Element, PeriodVI, GroupIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 63

    @property
    def isotopes_distribution(self):
        return {151: 0.4781, 152: 0., 153: 0.5219}

    @property
    def isotopes_masses(self):
        return {151: 150.919846, 152: 151.921744, 153: 152.921226}

    @property
    def _common_valences(self):
        return 0, 3

    @property
    def _valences_exceptions(self):
        return ((3, False, 0, ()),
                (0, False, 0, ((1, 'O'), (1, 'O'))),
                (0, False, 0, ((1, 'F'), (1, 'F'))),
                (0, False, 0, ((1, 'Cl'), (1, 'Cl'))),
                (0, False, 0, ((1, 'Br'), (1, 'Br'))),
                (0, False, 0, ((1, 'I'), (1, 'I'))),
                (0, False, 0, ((1, 'F'), (1, 'Cl'))),
                (0, False, 0, ((1, 'F'), (1, 'Br'))),
                (0, False, 0, ((1, 'F'), (1, 'I'))),
                (0, False, 0, ((1, 'C'), (1, 'C'))),
                (0, False, 0, ((1, 'H'), (1, 'H'))),
                (0, False, 0, ((2, 'O'),)))

    @property
    def atomic_radius(self):
        return 2.31

    @property
    def mdl_isotope(self):
        return 152


class Gd(Element, PeriodVI, GroupIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 64

    @property
    def isotopes_distribution(self):
        return {152: 0.002, 153: 0., 154: 0.0218, 155: 0.148, 156: 0.2047, 157: 0.1565, 158: 0.2484, 160: 0.2186}

    @property
    def isotopes_masses(self):
        return {152: 151.919788, 153: 152.921750, 154: 153.920862, 155: 154.922619, 156: 155.922120, 157: 156.923957,
                158: 157.924101, 160: 159.927051}

    @property
    def _common_valences(self):
        return 0, 3

    @property
    def _valences_exceptions(self):
        return (3, False, 0, ()),

    @property
    def atomic_radius(self):
        return 2.33

    @property
    def mdl_isotope(self):
        return 157


class Tb(Element, PeriodVI, GroupIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 65

    @property
    def isotopes_distribution(self):
        return {159: 1.0, 160: 0.}

    @property
    def isotopes_masses(self):
        return {159: 158.925343, 160: 159.927168}

    @property
    def _common_valences(self):
        return 0, 3

    @property
    def _valences_exceptions(self):
        return ((3, False, 0, ()),
                (0, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'))),
                (0, False, 0, ((2, 'O'), (2, 'O'))))

    @property
    def atomic_radius(self):
        return 2.25

    @property
    def mdl_isotope(self):
        return 159


class Dy(Element, PeriodVI, GroupIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 66

    @property
    def isotopes_distribution(self):
        return {156: 0.0006, 158: 0.001, 160: 0.0234, 161: 0.1891, 162: 0.2551, 163: 0.249, 164: 0.2818}

    @property
    def isotopes_masses(self):
        return {156: 155.924278, 158: 157.924405, 160: 159.925194, 161: 160.92693, 162: 161.926795, 163: 162.928728,
                164: 163.929171}

    @property
    def _common_valences(self):
        return 0, 3

    @property
    def _valences_exceptions(self):
        return ((3, False, 0, ()),
                (0, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'))),
                (0, False, 0, ((2, 'O'), (2, 'O'))))

    @property
    def atomic_radius(self):
        return 2.28

    @property
    def mdl_isotope(self):
        return 163


class Ho(Element, PeriodVI, GroupIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 67

    @property
    def isotopes_distribution(self):
        return {165: 1.0, 166: 0.}

    @property
    def isotopes_masses(self):
        return {165: 164.930319, 166: 165.932284}

    @property
    def _common_valences(self):
        return 0, 3

    @property
    def _valences_exceptions(self):
        return ((3, False, 0, ()),
                (0, False, 0, ((1, 'O'), (1, 'O'))),
                (0, False, 0, ((1, 'F'), (1, 'F'))),
                (0, False, 0, ((1, 'Cl'), (1, 'Cl'))),
                (0, False, 0, ((1, 'Br'), (1, 'Br'))),
                (0, False, 0, ((1, 'I'), (1, 'I'))),
                (0, False, 0, ((1, 'F'), (1, 'Cl'))),
                (0, False, 0, ((1, 'F'), (1, 'Br'))),
                (0, False, 0, ((1, 'F'), (1, 'I'))),
                (0, False, 0, ((1, 'C'), (1, 'C'))),
                (0, False, 0, ((1, 'H'), (1, 'H'))),
                (0, False, 0, ((2, 'O'),)))

    @property
    def atomic_radius(self):
        return 2.26

    @property
    def mdl_isotope(self):
        return 165


class Er(Element, PeriodVI, GroupIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 68

    @property
    def isotopes_distribution(self):
        return {162: 0.0014, 164: 0.0161, 166: 0.3361, 167: 0.2293, 168: 0.2678, 170: 0.1493}

    @property
    def isotopes_masses(self):
        return {162: 161.928775, 164: 163.929197, 166: 165.93029, 167: 166.932045, 168: 167.932368, 170: 169.93546}

    @property
    def _common_valences(self):
        return 0, 3

    @property
    def _valences_exceptions(self):
        return (3, False, 0, ()),

    @property
    def atomic_radius(self):
        return 2.26

    @property
    def mdl_isotope(self):
        return 167


class Tm(Element, PeriodVI, GroupIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 69

    @property
    def isotopes_distribution(self):
        return {169: 1.0, 170: 0.}

    @property
    def isotopes_masses(self):
        return {169: 168.934211, 170: 169.935801}

    @property
    def _common_valences(self):
        return 0, 3

    @property
    def _valences_exceptions(self):
        return ((3, False, 0, ()),
                (0, False, 0, ((1, 'O'), (1, 'O'))),
                (0, False, 0, ((1, 'F'), (1, 'F'))),
                (0, False, 0, ((1, 'Cl'), (1, 'Cl'))),
                (0, False, 0, ((1, 'Br'), (1, 'Br'))),
                (0, False, 0, ((1, 'I'), (1, 'I'))),
                (0, False, 0, ((1, 'F'), (1, 'Cl'))),
                (0, False, 0, ((1, 'F'), (1, 'Br'))),
                (0, False, 0, ((1, 'F'), (1, 'I'))),
                (0, False, 0, ((1, 'C'), (1, 'C'))),
                (0, False, 0, ((1, 'H'), (1, 'H'))),
                (0, False, 0, ((2, 'O'),)))

    @property
    def atomic_radius(self):
        return 2.22

    @property
    def mdl_isotope(self):
        return 169


class Yb(Element, PeriodVI, GroupIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 70

    @property
    def isotopes_distribution(self):
        return {168: 0.0013, 169: 0., 170: 0.0304, 171: 0.1428, 172: 0.2183, 173: 0.1613, 174: 0.3183, 176: 0.1276}

    @property
    def isotopes_masses(self):
        return {168: 167.933894, 169: 168.935190, 170: 169.934759, 171: 170.936322, 172: 171.936378, 173: 172.938207,
                174: 173.938858, 176: 175.942568}

    @property
    def _common_valences(self):
        return 0, 3

    @property
    def _valences_exceptions(self):
        return ((3, False, 0, ()),
                (0, False, 0, ((1, 'O'), (1, 'O'))),
                (0, False, 0, ((1, 'F'), (1, 'F'))),
                (0, False, 0, ((1, 'Cl'), (1, 'Cl'))),
                (0, False, 0, ((1, 'Br'), (1, 'Br'))),
                (0, False, 0, ((1, 'I'), (1, 'I'))),
                (0, False, 0, ((1, 'F'), (1, 'Cl'))),
                (0, False, 0, ((1, 'F'), (1, 'Br'))),
                (0, False, 0, ((1, 'F'), (1, 'I'))),
                (0, False, 0, ((1, 'C'), (1, 'C'))),
                (0, False, 0, ((1, 'H'), (1, 'H'))),
                (0, False, 0, ((2, 'O'),)))

    @property
    def atomic_radius(self):
        return 2.22

    @property
    def mdl_isotope(self):
        return 173


class Lu(Element, PeriodVI, GroupIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 71

    @property
    def isotopes_distribution(self):
        return {175: 0.9741, 176: 0.0259, 177: 0.}

    @property
    def isotopes_masses(self):
        return {175: 174.940768, 176: 175.942682, 177: 176.943758}

    @property
    def _common_valences(self):
        return 0, 3

    @property
    def _valences_exceptions(self):
        return (3, False, 0, ()),

    @property
    def atomic_radius(self):
        return 2.17

    @property
    def mdl_isotope(self):
        return 175


class Ac(Element, PeriodVII, GroupIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 89

    @property
    def isotopes_distribution(self):
        return {225: 0., 227: 1.0}

    @property
    def isotopes_masses(self):
        return {225: 225.023230, 227: 227.027752}

    @property
    def _common_valences(self):
        return 0, 3

    @property
    def _valences_exceptions(self):
        return (3, False, 0, ()),

    @property
    def atomic_radius(self):
        return 2.17  # unknown, taken radius of previous element in group

    @property
    def mdl_isotope(self):
        return 227


class Th(Element, PeriodVII, GroupIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 90

    @property
    def isotopes_distribution(self):
        return {227: 0., 232: 1.0}

    @property
    def isotopes_masses(self):
        return {227: 227.027704, 232: 232.038050}

    @property
    def _common_valences(self):
        return 0, 4

    @property
    def _valences_exceptions(self):
        return ((4, False, 0, ()),
                (0, False, 0, ((1, 'Br'), (1, 'Br'), (1, 'Br'))),
                (0, False, 0, ((1, 'I'), (1, 'I'), (1, 'I'))),
                (0, False, 0, ((1, 'Br'), (1, 'Br'))),
                (0, False, 0, ((1, 'I'), (1, 'I'))),
                (0, False, 0, ((1, 'H'), (1, 'H'))))

    @property
    def atomic_radius(self):
        return 2.17  # unknown, taken radius of previous element in group

    @property
    def mdl_isotope(self):
        return 232


class Pa(Element, PeriodVII, GroupIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 91

    @property
    def isotopes_distribution(self):
        return {231: 1.0, 233: 0.}

    @property
    def isotopes_masses(self):
        return {231: 231.035879, 233: 233.040247}

    @property
    def _common_valences(self):
        return 0, 4, 5

    @property
    def _valences_exceptions(self):
        return ((4, False, 0, ()),
                (0, False, 0, ((1, 'H'), (1, 'H'), (1, 'H'))),
                (0, False, 0, ((2, 'O'),)))

    @property
    def atomic_radius(self):
        return 2.17  # unknown, taken radius of previous element in group

    @property
    def mdl_isotope(self):
        return 231


class U(Element, PeriodVII, GroupIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 92

    @property
    def isotopes_distribution(self):
        return {234: 5.5e-05, 235: 0.0072, 238: 0.992745}

    @property
    def isotopes_masses(self):
        return {234: 234.040946, 235: 235.043923, 238: 238.050783}

    @property
    def _common_valences(self):
        return 0, 3, 4, 5, 6

    @property
    def _valences_exceptions(self):
        return ((3, False, 0, ()), (4, False, 0, ()),
                (2, False, 0, ((2, 'O'), (2, 'O'))))

    @property
    def atomic_radius(self):
        return 2.17  # unknown, taken radius of previous element in group

    @property
    def mdl_isotope(self):
        return 238


class Np(Element, PeriodVII, GroupIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 93

    @property
    def isotopes_distribution(self):
        return {237: 1.0}

    @property
    def isotopes_masses(self):
        return {237: 237.048173}

    @property
    def _common_valences(self):
        return 0, 2, 3, 4, 5, 6, 7

    @property
    def _valences_exceptions(self):
        return ((3, False, 0, ()), (4, False, 0, ()),
                (1, False, 0, ((2, 'O'), (2, 'O'))),
                (2, False, 0, ((2, 'O'), (2, 'O'))))

    @property
    def atomic_radius(self):
        return 2.17  # unknown, taken radius of previous element in group

    @property
    def mdl_isotope(self):
        return 237


class Pu(Element, PeriodVII, GroupIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 94

    @property
    def isotopes_distribution(self):
        return {239: 1.0, 242: 0.}

    @property
    def isotopes_masses(self):
        return {239: 239.052163, 242: 242.058743}

    @property
    def _common_valences(self):
        return 0, 3, 4, 5, 6

    @property
    def _valences_exceptions(self):
        return ((3, False, 0, ()), (4, False, 0, ()),
                (1, False, 0, ((2, 'O'), (2, 'O'))),
                (2, False, 0, ((2, 'O'), (2, 'O'))),
                (0, False, 0, ((2, 'Se'), )),
                (0, False, 0, ((2, 'S'),)),
                (0, False, 0, ((2, 'Te'),)),
                (0, False, 0, ((2, 'O'),)),
                (0, False, 0, ((1, 'Cl'), (1, 'Cl'))),
                (0, False, 0, ((1, 'Br'), (1, 'Br'))),
                (0, False, 0, ((1, 'I'), (1, 'I'))),
                (0, False, 0, ((1, 'H'), (1, 'H'))))

    @property
    def atomic_radius(self):
        return 2.17  # unknown, taken radius of previous element in group

    @property
    def mdl_isotope(self):
        return 244


class Am(Element, PeriodVII, GroupIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 95

    @property
    def isotopes_distribution(self):
        return {241: 1.0, 243: 0.}

    @property
    def isotopes_masses(self):
        return {241: 241.056829, 243: 243.061380}

    @property
    def _common_valences(self):
        return 0, 2, 3, 4

    @property
    def _valences_exceptions(self):
        return (3, False, 0, ()),

    @property
    def atomic_radius(self):
        return 2.17  # unknown, taken radius of previous element in group

    @property
    def mdl_isotope(self):
        return 243


class Cm(Element, PeriodVII, GroupIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 96

    @property
    def isotopes_distribution(self):
        return {243: 0., 244: 1.0, 248: 0.}

    @property
    def isotopes_masses(self):
        return {243: 243.061389, 244: 244.062753, 248: 248.072349}

    @property
    def _common_valences(self):
        return 0, 3, 4

    @property
    def _valences_exceptions(self):
        return (0, False, 0, ((2, 'O'),)), (0, False, 0, ((1, 'H'), (1, 'H')))

    @property
    def atomic_radius(self):
        return 2.17  # unknown, taken radius of previous element in group

    @property
    def mdl_isotope(self):
        return 247


class Bk(Element, PeriodVII, GroupIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 97

    @property
    def isotopes_distribution(self):
        return {249: 1.0}

    @property
    def isotopes_masses(self):
        return {249: 249.074987}

    @property
    def _common_valences(self):
        return 0, 2, 3, 4

    @property
    def _valences_exceptions(self):
        return (3, False, 0, ()), (4, False, 0, ())

    @property
    def atomic_radius(self):
        return 2.17  # unknown, taken radius of previous element in group

    @property
    def mdl_isotope(self):
        return 247


class Cf(Element, PeriodVII, GroupIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 98

    @property
    def isotopes_distribution(self):
        return {249: 1.0}

    @property
    def isotopes_masses(self):
        return {249: 249.074854}

    @property
    def _common_valences(self):
        return 0, 2, 3, 4

    @property
    def _valences_exceptions(self):
        return (3, False, 0, ()),

    @property
    def atomic_radius(self):
        return 2.17  # unknown, taken radius of previous element in group

    @property
    def mdl_isotope(self):
        return 251


class Es(Element, PeriodVII, GroupIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 99

    @property
    def isotopes_distribution(self):
        return {252: 1.0}

    @property
    def isotopes_masses(self):
        return {252: 252.08298}

    @property
    def _common_valences(self):
        return 0, 2, 3

    @property
    def _valences_exceptions(self):
        return (3, False, 0, ()),

    @property
    def atomic_radius(self):
        return 2.17  # unknown, taken radius of previous element in group

    @property
    def mdl_isotope(self):
        return 252


class Fm(Element, PeriodVII, GroupIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 100

    @property
    def isotopes_distribution(self):
        return {257: 1.0}

    @property
    def isotopes_masses(self):
        return {257: 257.095106}

    @property
    def _common_valences(self):
        return 0, 2, 3

    @property
    def _valences_exceptions(self):
        return (3, False, 0, ()),

    @property
    def atomic_radius(self):
        return 2.17  # unknown, taken radius of previous element in group

    @property
    def mdl_isotope(self):
        return 257


class Md(Element, PeriodVII, GroupIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 101

    @property
    def isotopes_distribution(self):
        return {258: 1.0}

    @property
    def isotopes_masses(self):
        return {258: 258.098431}

    @property
    def _common_valences(self):
        return 0, 2, 3

    @property
    def _valences_exceptions(self):
        return (3, False, 0, ()),

    @property
    def atomic_radius(self):
        return 2.17  # unknown, taken radius of previous element in group

    @property
    def mdl_isotope(self):
        return 258


class No(Element, PeriodVII, GroupIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 102

    @property
    def isotopes_distribution(self):
        return {259: 1.0}

    @property
    def isotopes_masses(self):
        return {259: 259.10103}

    @property
    def _common_valences(self):
        return 0, 2, 3

    @property
    def _valences_exceptions(self):
        return (2, False, 0, ()),

    @property
    def atomic_radius(self):
        return 2.17  # unknown, taken radius of previous element in group

    @property
    def mdl_isotope(self):
        return 259


class Lr(Element, PeriodVII, GroupIII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 103

    @property
    def isotopes_distribution(self):
        return {266: 1.0}

    @property
    def isotopes_masses(self):
        return {266: 266.11983}

    @property
    def _common_valences(self):
        return 0, 3

    @property
    def _valences_exceptions(self):
        return (3, False, 0, ()),

    @property
    def atomic_radius(self):
        return 2.17  # unknown, taken radius of previous element in group

    @property
    def mdl_isotope(self):
        return 260


__all__ = ['Sc', 'Y',
           'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu',
           'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr']
