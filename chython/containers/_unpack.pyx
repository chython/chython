# -*- coding: utf-8 -*-
#
#  Copyright 2021 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from libc.math cimport ldexp
from chython.containers.bonds import Bond


@cython.boundscheck(False)
def unpack(const unsigned char[::1] data not None):
    cdef int isotope_shift
    cdef unsigned char a, b, c, d
    cdef unsigned int na, nct, i, j, n, m, pack_length, shift = 4, order_shift = 0, nb = 0

    cdef unsigned int[4095] mapping, atom, isotopes, hydrogens, neighbors, orders, cis_trans_1, cis_trans_2
    cdef unsigned int[8190] connections
    cdef int[4095] charges
    cdef bint[4095] radicals, is_tet, is_all, tet_sign, all_sign, ct_sign
    cdef double[4095] x, y
    cdef bint[4096] seen

    cdef object bond, py_n, py_m
    cdef dict py_charges, py_radicals, py_hydrogens, py_plane, py_bonds, tmp
    cdef dict py_atoms_stereo, py_allenes_stereo, py_cis_trans_stereo
    cdef list py_mapping, py_atoms, py_isotopes, py_bonds_flat

    # lets extract data
    a, b, c = data[1], data[2], data[3]
    na = a << 4| (b & 0xf0) >> 4
    nct = (b & 0x0f) << 8 | c

    for i in range(na):
        a, b = data[shift], data[shift + 1]
        mapping[i] = n = a << 4 | (b & 0xf0) >> 4
        seen[n] = False
        neighbors[i] = b & 0x0f
        nb += b & 0x0f

        a, b = data[shift + 2], data[shift + 3]
        if a & 0x80:
            is_tet[i] = True
            tet_sign[i] = a & 0x40
        else:
            is_tet[i] = False
        if a & 0x20:
            is_all[i] = True
            all_sign[i] = a & 0x10
        else:
            is_all[i] = False

        atom[i] = b & 0x7f
        isotope_shift = (a & 0x0f) << 1 | b >> 7
        if isotope_shift:
            isotopes[i] = common_isotopes[b & 0x7f] + isotope_shift
        else:
            isotopes[i] = 0

        a, b = data[shift + 4], data[shift + 5]
        x[i] = double_from2bytes(a, b)
        a, b = data[shift + 6], data[shift + 7]
        y[i] = double_from2bytes(a, b)

        a = data[shift + 8]
        hydrogens[i] = a >> 5
        charges[i] = ((a >> 1) & 0x0f) - 4
        radicals[i] = a & 0x01

        shift += 9

    nb //= 2
    for i in range(nb):
        a, b, c = data[shift], data[shift + 1], data[shift + 2]
        connections[i * 2] = a << 4| (b & 0xf0) >> 4
        connections[i * 2 + 1] = (b & 0x0f) << 8 | c
        shift += 3

    for i in range((nb // 5 + 1) if nb % 5 else (nb // 5)):
        a, b = data[shift], data[shift + 1]
        orders[i * 5] = (a >> 4) + 1
        orders[i * 5 + 1] = ((a >> 1) & 0x07) + 1
        orders[i * 5 + 2] = ((a & 0x01) << 2 | b >> 6) + 1
        orders[i * 5 + 3] = ((b >> 3) & 0x07) + 1
        orders[i * 5 + 4] = (b & 0x07) + 1
        shift += 2

    for i in range(nct):
        a, b, c, d = data[shift], data[shift + 1], data[shift + 2], data[shift + 3]
        cis_trans_1[i] = a << 4 | (b & 0xf0) >> 4
        cis_trans_2[i] = (b & 0x0f) << 8 | c
        ct_sign[i] = d & 0x01
        shift += 4
    pack_length = shift

    # define returned data
    py_mapping = []
    py_atoms = []
    py_isotopes = []
    py_bonds = {}
    py_charges = {}
    py_radicals = {}
    py_hydrogens = {}
    py_plane = {}
    py_atoms_stereo = {}
    py_allenes_stereo = {}
    py_cis_trans_stereo = {}
    py_bonds_flat = []

    shift = 0
    for i in range(na):
        n = mapping[i]
        py_n = n

        # fill intermediate data
        py_mapping.append(py_n)
        py_atoms.append(atom[i])
        py_isotopes.append(isotopes[i] or None)

        py_charges[py_n] = charges[i]
        py_radicals[py_n] = radicals[i]
        if hydrogens[i] == 7:
            py_hydrogens[py_n] = None
        else:
            py_hydrogens[py_n] = hydrogens[i]
        py_plane[py_n] = (x[i], y[i])

        if is_tet[i]:
            py_atoms_stereo[py_n] = tet_sign[i]
        if is_all[i]:
            py_allenes_stereo[py_n] = all_sign[i]

        tmp = {}
        py_bonds[py_n] = tmp
        seen[n] = True
        for j in range(shift, shift + neighbors[i]):
            m = connections[j]
            py_m = m
            if seen[m]:  # bond partially exists. need back-connection.
                tmp[py_m] = py_bonds[py_m][py_n]
            else:
                bond = object.__new__(Bond)
                bond._Bond__order = orders[order_shift]
                bond._Bond__n = py_n
                bond._Bond__m = py_m
                tmp[py_m] = bond
                py_bonds_flat.append(bond)
                order_shift += 1

        shift += neighbors[i]

    for i in range(nct):
        py_cis_trans_stereo[(cis_trans_1[i], cis_trans_2[i])] = ct_sign[i]

    return (py_mapping, py_atoms, py_isotopes,
            py_charges, py_radicals, py_hydrogens, py_plane, py_bonds,
            py_atoms_stereo, py_allenes_stereo, py_cis_trans_stereo, pack_length, py_bonds_flat)


cdef int[119] common_isotopes
common_isotopes[:] = [0, -15, -12, -9, -7, -5, -4, -2, 0, 3, 4, 7, 8, 11, 12, 15, 16, 19, 24, 23, 24, 29,
                      32, 35, 36, 39, 40, 43, 43, 48, 49, 54, 57, 59, 63, 64, 68, 69, 72, 73, 75, 77,
                      80, 82, 85, 87, 90, 92, 96, 99, 103, 106, 112, 111, 115, 117, 121, 123, 124, 125,
                      128, 129, 134, 136, 141, 143, 147, 149, 151, 153, 157, 159, 162, 165, 168, 170,
                      174, 176, 179, 181, 185, 188, 191, 193, 193, 194, 206, 207, 210, 211, 216, 215,
                      222, 221, 228, 227, 231, 231, 235, 236, 241, 242, 243, 244, 245, 254, 253, 254,
                      254, 262, 265, 265, 269, 262, 273, 273, 277, 281, 278]


cdef double double_from2bytes(unsigned char a, unsigned char b):
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
