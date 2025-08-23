# -*- coding: utf-8 -*-
# cython: language_level=3
#
#  Copyright 2021-2025 Ramil Nugmanov <nougmanoff@protonmail.com>
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
cimport cython
from cpython.mem cimport PyMem_Malloc, PyMem_Free
from libc.math cimport ldexp

from chython.containers import MoleculeContainer
from chython.containers.bonds import Bond
from chython.periodictable import (H, He, Li, Be, B, C, N, O, F, Ne, Na, Mg, Al, Si, P, S, Cl, Ar, K, Ca, Sc, Ti, V, Cr,
                                   Mn, Fe, Co, Ni, Cu, Zn, Ga, Ge, As, Se, Br, Kr, Rb, Sr, Y, Zr, Nb, Mo, Tc, Ru, Rh,
                                   Pd, Ag, Cd, In, Sn, Sb, Te, I, Xe, Cs, Ba, La, Ce, Pr, Nd, Pm, Sm, Eu, Gd, Tb, Dy,
                                   Ho, Er, Tm, Yb, Lu, Hf, Ta, W, Re, Os, Ir, Pt, Au, Hg, Tl, Pb, Bi, Po, At, Rn, Fr,
                                   Ra, Ac, Th, Pa, U, Np, Pu, Am, Cm, Bk, Cf, Es, Fm, Md, No, Lr, Rf, Db, Sg, Bh, Hs,
                                   Mt, Ds, Rg, Cn, Nh, Fl, Mc, Lv, Ts, Og)
from chython.periodictable.base.vector import Vector


# Format specification::
#
# Big endian bytes order
# 8 bit - 0x02 (current format specification)
# 12 bit - number of atoms
# 12 bit - cis/trans stereo block size
#
# Atom block 9 bytes (repeated):
# 12 bit - atom number
# 4 bit - number of neighbors
# 2 bit tetrahedron sign (00 - not stereo, 10 or 11 - has stereo)
# 2 bit - allene sign
# 5 bit - isotope (00000 - not specified, over = isotope - common_isotope + 16)
# 7 bit - atomic number (<=118)
# 32 bit - XY float16 coordinates
# 3 bit - hydrogens (0-7). Note: 7 == None
# 4 bit - charge (charge + 4. possible range -4 - 4)
# 1 bit - radical state
# Connection table: flatten list of neighbors. neighbors count stored in atom block.
# For example CC(=O)O - {1: [2], 2: [1, 3, 4], 3: [2], 4: [2]} >> [2, 1, 3, 4, 2, 2].
# Repeated block (equal to bonds count).
# 24 bit - paired 12 bit numbers.
# Bonds order block 3 bit per bond zero-padded to full byte at the end.
# Cis/trans data block (repeated):
# 24 bit - atoms pair
# 7 bit - zero padding. in future can be used for extra bond-level stereo, like atropoisomers.
# 1 bit - sign

@cython.nonecheck(False)
@cython.boundscheck(False)
@cython.cdivision(True)
@cython.wraparound(False)
def unpack(const unsigned char[::1] data not None):
    cdef unsigned char a, b, c, d, isotope, atomic_number, neighbors_count, s = 0, version, stereo, hydrogens
    cdef unsigned char *neighbors, *orders
    cdef unsigned short atoms_count, bonds_count = 0, cis_trans_count, order_count
    cdef unsigned short i, j, k = 0, n, m, buffer_b, shift = 0
    cdef unsigned short *mapping, *connections
    cdef unsigned int size, atoms_shift = 4, bonds_shift, order_shift, cis_trans_shift
    cdef unsigned char[4096] seen

    cdef object py_mol, py_bond, py_n, py_m, py_atom, py_nan_bool, py_vector
    cdef dict py_atoms, py_bonds, py_ngb
    cdef list py_cis_trans

    # read header
    version = data[0]
    a, b, c = data[1], data[2], data[3]
    atoms_count = a << 4| b >> 4
    cis_trans_count = (b & 0x0f) << 8 | c

    # allocate memory
    neighbors = <unsigned char*> PyMem_Malloc(atoms_count * sizeof(unsigned char))
    mapping = <unsigned short*> PyMem_Malloc(atoms_count * sizeof(unsigned short))

    if not neighbors or not mapping:
        raise MemoryError()

    py_mol = MoleculeContainer()
    py_atoms = py_mol._atoms
    py_bonds = py_mol._bonds
    py_cis_trans = []

    # unpack atom block to separate attributes arrays
    for i in range(atoms_count):
        a, b = data[atoms_shift], data[atoms_shift + 1]
        mapping[i] = n = a << 4 | b >> 4
        seen[n] = 0  # erase random value
        neighbors[i] = neighbors_count = b & 0x0f
        bonds_count += neighbors_count

        a, b = data[atoms_shift + 2], data[atoms_shift + 3]
        stereo = a >> 4
        if stereo == 0:
            py_nan_bool = None
        elif stereo == 0b0010:
            py_nan_bool = False
        elif stereo == 0b0011:
            py_nan_bool = True
        elif stereo == 0b1000:
            py_nan_bool = False
        else:  # if stereo == 0b1100:
            py_nan_bool = True

        atomic_number = b & 0x7f
        py_atom = object.__new__(elements[atomic_number])
        py_atom._extended_stereo = None
        py_n = n
        py_atoms[py_n] = py_atom
        py_bonds[py_n] = {}

        py_atom._stereo = py_nan_bool

        isotope = (a & 0x0f) << 1 | b >> 7
        if isotope:
            py_atom._isotope = common_isotopes[atomic_number] + isotope
        else:
            py_atom._isotope = None

        py_vector = object.__new__(Vector)
        a, b = data[atoms_shift + 4], data[atoms_shift + 5]
        py_vector.x = double_from_bytes(a, b)
        a, b = data[atoms_shift + 6], data[atoms_shift + 7]
        py_vector.y = double_from_bytes(a, b)
        py_atom._xy = py_vector

        a = data[atoms_shift + 8]
        hydrogens = a >> 5
        if hydrogens == 7:
            py_atom._implicit_hydrogens = None
        else:
            py_atom._implicit_hydrogens = hydrogens

        py_atom._charge = ((a >> 1) & 0x0f) - 4
        if a & 0x01:
            py_atom._is_radical = True
        else:
            py_atom._is_radical = False
        atoms_shift += 9

    # calculate bonds count and pack sections
    bonds_count /= 2

    if version == 2:
        order_count = bonds_count * 3
        if order_count % 8:
            order_count = order_count / 8 + 1
        else:
            order_count /= 8
    else:  # if version == 0:
        order_count = bonds_count / 5
        if bonds_count % 5:
            order_count += 1
        order_count *= 2

    bonds_shift = atoms_shift
    order_shift = bonds_shift + 3 * bonds_count
    cis_trans_shift = order_count + order_shift
    size = cis_trans_shift + 4 * cis_trans_count

    if bonds_count:
        # keep 4 extra cells for padding
        orders = <unsigned char*> PyMem_Malloc((bonds_count + 4) * sizeof(unsigned char))
        connections = <unsigned short*> PyMem_Malloc(2 * bonds_count * sizeof(unsigned short))

        # connection table is bidirected
        for i in range(0, 2 * bonds_count, 2):
            a, b, c = data[bonds_shift], data[bonds_shift + 1], data[bonds_shift + 2]
            connections[i] = a << 4| b >> 4
            connections[i + 1] = (b & 0x0f) << 8 | c
            bonds_shift += 3

        # collect flat bond order list
        i = 0
        if version == 2:
            for j in range(order_shift, cis_trans_shift):
                # 3 3 2 | 1 3 3 1 | 2 3 3
                a = data[j]
                if s == 1:
                    orders[i] = buffer_b | a >> 7
                    orders[i + 1] = (a >> 4) & 0x7
                    orders[i + 2] = (a >> 1) & 0x7
                    buffer_b = (a & 1) << 2
                    s = 2
                    i += 3
                elif s == 2:
                    orders[i] = buffer_b | a >> 6
                    orders[i + 1] = (a >> 3) & 0x7
                    orders[i + 2] = a & 0x7
                    i += 3
                    s = 0
                else:
                    orders[i] = a >> 5
                    orders[i + 1] = (a >> 2) & 0x7
                    buffer_b = (a & 0x3) << 1
                    s = 1
                    i += 2
        else:  # if version == 0:
            for j in range(order_shift, cis_trans_shift, 2):
                # 0 3 3 1 | 2 3 3
                a, b = data[j], data[j + 1]
                orders[i] = a >> 4
                orders[i + 1] = (a >> 1) & 0x7
                orders[i + 2] = (a & 0x1) << 2 | b >> 6
                orders[i + 3] = (b >> 3) & 0x7
                orders[i + 4] = b & 0x7
                i += 5

        for i in range(atoms_count):
            n = mapping[i]
            py_n = n  # shared py int obj

            py_ngb = py_bonds[py_n]
            seen[n] = 1

            neighbors_count = neighbors[i]
            for j in range(shift, shift + neighbors_count):
                m = connections[j]
                py_m = m
                if seen[m]:  # bond partially exists. need back-connection.
                    py_ngb[py_m] = py_bonds[py_m][py_n]
                else:
                    py_bond = object.__new__(Bond)
                    py_bond._order = orders[k] + 1
                    py_bond._stereo = None
                    py_ngb[py_m] = py_bond
                    k += 1
            shift += neighbors_count

        PyMem_Free(orders)
        PyMem_Free(connections)

    if cis_trans_count:
        for i in range(cis_trans_count):
            a, b = data[cis_trans_shift], data[cis_trans_shift + 1]
            c, d = data[cis_trans_shift + 2], data[cis_trans_shift + 3]
            py_n = a << 4 | b >> 4
            py_m = (b & 0x0f) << 8 | c
            if d:
                py_cis_trans.append((py_n, py_m, True))
            else:
                py_cis_trans.append((py_n, py_m, False))
            cis_trans_shift += 4

    PyMem_Free(neighbors)
    PyMem_Free(mapping)
    return py_mol, py_cis_trans, size


cdef short[119] common_isotopes
common_isotopes[:] = [0, -15, -12, -9, -7, -5, -4, -2, 0, 3, 4, 7, 8, 11, 12, 15, 16, 19, 24, 23, 24, 29,
                      32, 35, 36, 39, 40, 43, 43, 48, 49, 54, 57, 59, 63, 64, 68, 69, 72, 73, 75, 77,
                      80, 82, 85, 87, 90, 92, 96, 99, 103, 106, 112, 111, 115, 117, 121, 123, 124, 125,
                      128, 129, 134, 136, 141, 143, 147, 149, 151, 153, 157, 159, 162, 165, 168, 170,
                      174, 176, 179, 181, 185, 188, 191, 193, 193, 194, 206, 207, 210, 211, 216, 215,
                      222, 221, 228, 227, 231, 231, 235, 236, 241, 242, 243, 244, 245, 254, 253, 254,
                      254, 262, 265, 265, 269, 262, 273, 273, 277, 281, 278]

cdef list elements
elements = [None, H, He, Li, Be, B, C, N, O, F, Ne, Na, Mg, Al, Si, P, S, Cl, Ar, K, Ca, Sc, Ti, V, Cr, Mn, Fe, Co,
            Ni, Cu, Zn, Ga, Ge, As, Se, Br, Kr, Rb, Sr, Y, Zr, Nb, Mo, Tc, Ru, Rh, Pd, Ag, Cd, In, Sn, Sb, Te, I, Xe,
            Cs, Ba, La, Ce, Pr, Nd, Pm, Sm, Eu, Gd, Tb, Dy, Ho, Er, Tm, Yb, Lu, Hf, Ta, W, Re, Os, Ir, Pt, Au, Hg,
            Tl, Pb, Bi, Po, At, Rn, Fr, Ra, Ac, Th, Pa, U, Np, Pu, Am, Cm, Bk, Cf, Es, Fm, Md, No, Lr, Rf, Db, Sg,
            Bh, Hs, Mt, Ds, Rg, Cn, Nh, Fl, Mc, Lv, Ts, Og]


cdef double double_from_bytes(unsigned char a, unsigned char b):
    cdef bint sign
    cdef int e
    cdef unsigned int f
    cdef double x

    sign = a >> 7
    e = (a >> 2) & 0x1f
    f = ((a & 0x03) << 8) | b

    x = f / 1024.
    if e:
        x += 1.
        e -= 15
    else:
        e = -14

    x = ldexp(x, e)
    if sign:
        return -x
    return x
