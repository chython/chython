# -*- coding: utf-8 -*-
#
#  Copyright 2024 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from ..files import smarts


_groups = {
    'Alkene': '[C;z2;x0]=[C;x0;z2]',
    'Alkene hetero': '[C;z2;x1,x2]=[C;z2]',
    'Alkene terminal': '[C;z2;x0;D1]=[C;x0;z2]',
    'Alkene hetero terminal': '[C;z2;x0;D1]=[C;x1,x2;z2]',
    'Alkyne': '[C;z3;x0]#[C;x0]',
    'Alkyne hetero': '[C;z3;x1]#[C]',
    'Alkyne terminal': '[C;z3;x0;D1]#[C;x0]',
    'Alkyne hetero terminal': '[C;z3;x0;D1]#[C;x1]',

    'Alkyl Halide': '[F,Cl,Br,I;D1][C;x1;z1]',
    'Cyclopropyl Halide': '[F,Cl,Br,I;D1][C;x1;z1;r3]',

    'Aryl Halide': '[F,Cl,Br,I;D1]-[C;a]',
    'Aryl Fluoride': '[F;D1]-[C;a]',
    'Aryl Chloride': '[Cl;D1]-[C;a]',
    'Aryl Bromide': '[Br;D1]-[C;a]',
    'Aryl Iodide': '[I;D1]-[C;a]',

    'Aryl Halide SNAr alpha': '[F,Cl,Br,I;D1][C;a]:N',
    'Aryl Halide SNAr gamma': '[F,Cl,Br,I;D1][C;a]:C:C:N',

    'Aryl Fluoride SNAr alpha': '[F;D1][C;a]:N',
    'Aryl Fluoride SNAr gamma': '[F;D1][C;a]:C:C:N',
    'Aryl Chloride SNAr alpha': '[Cl;D1][C;a]:N',
    'Aryl Chloride SNAr gamma': '[Cl;D1][C;a]:C:C:N',
    'Aryl Bromide SNAr alpha': '[Br;D1][C;a]:N',
    'Aryl Bromide SNAr gamma': '[Br;D1][C;a]:C:C:N',
    'Aryl Iodide SNAr alpha': '[I;D1][C;a]:N',
    'Aryl Iodide SNAr gamma': '[I;D1][C;a]:C:C:N',

    'Alcohol aliphatic': '[O;D1;x0;z1][C;x1;z1]',
    'Alcohol primary or secondary aliphatic': '[O;D1;x0;z1][C;D1,D2,D3;x1;z1]',
    'Alcohol tertiary aliphatic': '[O;D1;x0;z1][C;D4;x1;z1]',

    'Alcohol aromatic': '[O;D1;x0;z1][C;a]',

    'Aldehyde': '[O;z2;x0]=[C;D1,D2;x1;z2]',
    'Aldehyde aliphatic': '[O;z2;x0]=[C;D2;x1;z2][C;z1]',
    'Aldehyde aromatic': '[O;z2;x0]=[C;D2;x1;z2][C;a]',
    'Ketone': '[O;z2;x0]=[C;D3;x1;z2]',

    'Carboxylic Acid': '[O;D1;z1;x0][C;D3;x2;z2]=O',
    'Carboxylic Acid aliphatic': '[O;D1;z1;x0][C;D3;x2;z2](=O)[C;z1]',
    'Carboxylic Acid aromatic': '[O;D1;z1;x0][C;D3;x2;z2](=O)[C;a]',

    'Carboxylic Acid Ester': '[O;z2;x0]=[C;D3;x2;z2][O;D2;x0]',
    'Carboxylic Acid Halide': '[F,Cl,Br,I;D1][C;D3;x2;z2]=O',

    'Amine primary': '[N;D1;x0;z1][C;z1,z4;x1]',
    'Amine primary aliphatic': '[N;D1;x0;z1][C;z1;x1]',
    'Amine primary aromatic': '[N;D1;x0;z1][C;a]',
    'Amine secondary': '[N;D2;x0;z1]([C;z1,z4;x1])[C;z1,z4;x1]',
    'Amine secondary aliphatic': '[N;D2;x0;z1]([C;z1;x1])[C;z1;x1]',
    'Amine secondary aromatic': '[N;D2;x0;z1]([C;a])[C;z1,z4;x1]',
    'Amine tertiary': '[N;D3;x0;z1]([C;z1,z4;x1])([C;z1,z4;x1])[C;z1,z4;x1]',
    'Amine cyclic': '[N;D2;x0;z1;r4,r5,r6,r7,r8]([C;z1;x1])[C;z1;x1]',

    'Alpha-aminoacid': '[N;D1;x0;z1][C;z1;x1][C;D3;z2;x2](=O)[O;D1]',
    'Beta-aminoacid': '[N;D1;x0;z1][C;z1;x1][C;z1][C;D3;z2;x2](=O)[O;D1]',
    'Gamma-aminoacid': '[N;D1;x0;z1][C;z1;x1][C;z1][C;z1][C;D3;z2;x2](=O)[O;D1]',

    'Alpha-aminoacid N-protected': '[N;D2;x0;z1]([C;z1;x1][C;D3;z2;x2](=O)[O;D1])[C;D3;z2;x3](=O)[O;D2;x0]C',
    'Beta-aminoacid N-protected': '[N;D2;x0;z1]([C;z1;x1][C;z1][C;D3;z2;x2](=O)[O;D1])[C;D3;z2;x3](=O)[O;D2;x0]C',
    'Gamma-aminoacid N-protected': '[N;D2;x0;z1]([C;z1;x1][C;z1][C;z1][C;D3;z2;x2](=O)[O;D1])[C;D3;z2;x3](=O)[O;D2;x0]C',

    'Alpha-aminoacid O-protected': '[N;D1;x0;z1][C;z1;x1][C;D3;z2;x2](=O)[O;D2;x0][C;z1;x1]',
    'Beta-aminoacid O-protected': '[N;D1;x0;z1][C;z1;x1][C;z1][C;D3;z2;x2](=O)[O;D2;x0][C;z1;x1]',
    'Gamma-aminoacid O-protected': '[N;D1;x0;z1][C;z1;x1][C;z1][C;z1][C;D3;z2;x2](=O)[O;D2;x0][C;z1;x1]',

    'Aryl Sulfone': '[S;D4;z3;x2](=O)(=O)(-[C;a])-[C;x1]',
    'Azide': '[N-;D1;x1;z2]=[N+;D2;x2;z3]=[N;D2;x1;z2]',

    'Boronic Acid': '[B;D3;x2;z1]([O;D1])[O;D1]',
    'Boronic Acid Ester': '[B;D3;x2;z1]([O;D2;x1])([O;D2;x1])-;!@C',
    'Boronic Acid Ester aliphatic': '[B;D3;x2;z1]([O;D2;x1])([O;D2;x1])-;!@[C;x1;z1]',
    'Boronic Acid Ester aromatic': '[B;D3;x2;z1]([O;D2;x1])([O;D2;x1])-;!@[C;a]',
    'Trifluoroborate': '[B;D4;x3;z1;-](F)(F)(F)',

    'Thiol aliphatic': '[S;D1;x0;z1][C;x1;z1]',
    'Thiol aromatic': '[S;D1;x0;z1][C;a]',

    'Hydrazine aromatic': '[N;D1;x1;z1][N;D2;x1;z1][C;a]',
    'Hydrazine aliphatic': '[N;D1;x1;z1][N;D2;x1;z1][C;x1;z1]',

    'Isocyanate': '[O;D1;z2;x0]=[C;D2;x2;z3]=[N;D2;x0;z2]',
    'Isothiocyanate': '[S;D1;z2;x0]=[C;D2;x2;z3]=[N;D2;x0;z2]',
    'Nitrile': '[N;D1;z3;x0]#[C;D2;x1]',
    'Isonitrile': '[C-;D1;x1;z3]#[N+;D2;x0]',
    'Sulfonyl Halide': '[S;D4;z3;x3]([F,Cl,Br,I;D1])(=O)=O',
    'Lactam': '[O;D1;x0;z2]=[C;D3;x2](-;@[N;x0;z1]-;@[C;x1;z1])[C;x0;z1]',
    'Cyclic Anhydride': '[O;D1;x0;z2]=[C;D3;x2;z2;r4,r5,r6,r7,r8][O;D2;x0][C;D3;x2;z2]=O'
}

_smarts = {n.lower().replace(' ', '_'): smarts(s) for n, s in _groups.items()}
__all__ = list(_smarts)


def __getattr__(name):
    try:
        return _smarts[name]
    except KeyError:
        raise AttributeError


def __dir__():
    return __all__
