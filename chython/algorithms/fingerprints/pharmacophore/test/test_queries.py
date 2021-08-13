# -*- coding: utf-8 -*-
#
#  Copyright 2020, 2021 Aleksandr Sizov <murkyrussian@gmail.com>
#  Copyright 2020, 2021 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from CGRtools import smiles
from .._acceptor import queries as acceptor_queries, banned as acceptor_banned
from .._acidic import queries as acidic_queries
from .._basic import queries as basic_queries, banned as basic_banned
from .._donor import queries as donor_queries


def main_atom_ids(molecule, query):
    return {x[1] for x in query.get_mapping(molecule)}


mol = smiles('[cH:20]1[cH:21][c:22]([cH:24][cH:25][c:19]1[CH:12]([c:13]2[cH:18][cH:17][cH:16][cH:15][cH:14]2)'
             '[N:11]3[CH2:10][CH2:9][N:8]([CH2:27][CH2:26]3)[CH2:7][CH2:6][O:5][CH2:4][C:2](=[O:1])[OH:3])[Cl:23]')
frankensteins_fiend = smiles('[c:5]1[o:6][c:7][c:8][o:1][n:2]1[CH:3]=[S:4]')
pyrrole = smiles('[cH:3]1[cH:4][cH:5][nH:1][cH:2]1')
more_nitro = smiles('[C:10](=[O:11])([NH2:12])[c:3]1[cH:4][c:5]([N:7]=[N+:8]=[N-:9])[cH:6][c:1]([NH2:13])[cH:2]1')


def test_banned_acceptor():
    f = frankensteins_fiend
    carboxyl_q, amide_like, os_n_arom, os_c_n_arom = acceptor_banned
    assert main_atom_ids(mol, carboxyl_q) == {3, }
    assert main_atom_ids(f, os_n_arom) == {1, }
    assert main_atom_ids(f, os_c_n_arom) == {6, }


def test_allowed_acceptor():
    aromatic_n, aromatic_ox_sulf, *others = acceptor_queries
    assert main_atom_ids(frankensteins_fiend, aromatic_n) == {2, }
    assert main_atom_ids(frankensteins_fiend, aromatic_ox_sulf) == {1, 6}
    assert main_atom_ids(mol, acceptor_queries[2]) == {3, 5}
    assert main_atom_ids(mol, acceptor_queries[3]) == {8, 11}


def test_acidic():
    q = acidic_queries[0]
    assert main_atom_ids(mol, q) == {2, }


def test_basic():
    plus_n, aniline_n, sp3_n = basic_queries
    amide = basic_banned[0]
    assert main_atom_ids(more_nitro, plus_n) == {8, }
    assert main_atom_ids(more_nitro, aniline_n) == {7, 13}
    assert main_atom_ids(more_nitro, sp3_n) == {12, }
    assert main_atom_ids(more_nitro, amide) == {12, }
    assert main_atom_ids(mol, sp3_n) == {8, 11}


def test_donor():
    aromatic_n, n_withH, oh_sh = donor_queries
    assert main_atom_ids(pyrrole, aromatic_n) == {1, }
    assert main_atom_ids(more_nitro, n_withH) == {12, 13}
    assert main_atom_ids(mol, oh_sh) == {3, }
