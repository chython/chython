# -*- coding: utf-8 -*-
#
#  Copyright 2022 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from libc.math cimport ldexp, frexp


@cython.nonecheck(False)
@cython.boundscheck(False)
@cython.cdivision(True)
@cython.wraparound(False)
def pack(object molecule):
    cdef bint b  # binary flag
    cdef char charge
    cdef unsigned char atomic_number, isotope, bond, s = 0, buffer_b, buffer_o
    cdef unsigned char *p, *data
    cdef unsigned short atoms_count, bonds_count = 0, cis_trans_count, n, m
    cdef unsigned int size, atoms_shift = 4, bonds_shift, order_shift, cis_trans_shift  # can be > 2^16
    cdef unsigned char[4096] stereo, hcr, seen
    cdef unsigned int[4096] xy  # 2 * 16bit

    cdef bytes py_pack
    cdef dict py_ngb, py_atoms, py_bonds, py_charges, py_radicals, py_hydrogens, py_plane
    cdef dict py_cis_trans_stereo, py_atoms_stereo, py_allenes_stereo
    cdef tuple py_tuple
    cdef object py_atom, py_bond, py_nan_int, py_obj

    # map molecule to vars
    py_atoms = molecule._atoms
    py_bonds = molecule._bonds
    py_charges = molecule._charges
    py_radicals = molecule._radicals
    py_hydrogens = molecule._hydrogens
    py_cis_trans_stereo = molecule._cis_trans_stereo
    py_atoms_stereo = molecule._atoms_stereo
    py_allenes_stereo = molecule._allenes_stereo
    py_plane = molecule._plane

    # calculate elements count
    atoms_count = len(py_atoms)
    cis_trans_count = len(py_cis_trans_stereo)

    for py_ngb in py_bonds.values():
        bonds_count += len(py_ngb)
    bonds_count /= 2  # graph is bidirected

    # calculate pack blocks entries
    size = bonds_count * 3  # bonds bits
    if size % 8:  # partial byte fill
        size = size / 8 + 1
    else:
        size /= 8
    bonds_shift = 4 + 9 * atoms_count  # connection table starting byte
    order_shift = bonds_shift + 3 * bonds_count  # bond orders block starting byte
    cis_trans_shift = size + order_shift  # cis-trans block starting byte
    size = cis_trans_shift + 4 * cis_trans_count  # total pack size

    # allocate pack in memory
    data = <unsigned char *> PyMem_Malloc(size * sizeof(unsigned char))
    if not data:
        raise MemoryError()

    # precalculate atom attrs
    # should be done independently, due to possible randomness in dicts order.
    # 3 bit - hydrogens (0-7) | 4 bit - charge | 1 bit - radical
    for n, py_nan_int in py_hydrogens.items():
        if py_nan_int is None:
            hcr[n] = 0xe0  # 0b11100000
        else:
            hcr[n] = <unsigned char> py_nan_int << 5
    for n, charge in py_charges.items():
        hcr[n] |= (charge + 4) << 1
    for n, b in py_radicals.items():
        if b:  # lazy memory access
            hcr[n] |= 1

    # 2 float16 big endian
    for n, py_tuple in py_plane.items():
        p = <unsigned char *> &xy[n]
        double_to_float16(py_tuple[0], &p[0])
        double_to_float16(py_tuple[1], &p[2])

        # erase random data
        seen[n] = 0
        stereo[n] = 0

    # 2 bit tetrahedron | 2 bit allene | 0000
    for n, b in py_atoms_stereo.items():
        stereo[n] = 0xc0 if b else 0x80
    for n, b in py_allenes_stereo.items():
        stereo[n] = 0x30 if b else 0x20

    # start pack collection
    data[0] = 2  # header. specification version 2
    data[1] = atoms_count >> 4  # 5-12b of atom count value
    data[2] = atoms_count << 4 | cis_trans_count >> 8  # 1-4b of atom count value, 9-12b of cis-trans count value
    data[3] = cis_trans_count  # 1-8b of cis-trans count value

    b = True  # init connection table flag
    for py_obj, py_atom in py_atoms.items():
        py_ngb = py_bonds[py_obj]
        n = py_obj  # cast to C
        seen[n] = 1
        p = <unsigned char *> &xy[n]  # XY
        atomic_number = py_atom.atomic_number
        py_nan_int = py_atom._Core__isotope  # direct access
        if py_nan_int is None:
            isotope = 0
        else:
            isotope = <short> py_nan_int - common_isotopes[atomic_number]

        data[atoms_shift] = n >> 4  # 5-12b AN
        data[atoms_shift + 1] = n << 4 | len(py_ngb)  # 1-4b AN, 4b NC
        data[atoms_shift + 2] = stereo[n] | isotope >> 1  # TS , AS , 4b I
        data[atoms_shift + 3] = isotope << 7 | atomic_number  # 1bI , A
        data[atoms_shift + 4] = p[0]
        data[atoms_shift + 5] = p[1]
        data[atoms_shift + 6] = p[2]
        data[atoms_shift + 7] = p[3]
        data[atoms_shift + 8] = hcr[n]
        atoms_shift += 9

        # collect connection table
        for m, py_bond in py_ngb.items():
            if b:  # 8 + 4
                data[bonds_shift] = m >> 4
                bonds_shift += 1
                buffer_b = m << 4
                b = False  # switch
            else:  # 4 + 8
                data[bonds_shift] = buffer_b | m >> 8
                bonds_shift += 1
                data[bonds_shift] = m
                bonds_shift += 1  # next free 3 bytes block
                b = True

            if not seen[m]:
                bond = <unsigned char> py_bond._Bond__order - 1
                # 3 3 2 | 1 3 3 1 | 2 3 3
                if s == 0:
                    buffer_o = bond << 5
                    s = 1
                elif s == 1:
                    buffer_o |= bond << 2
                    s = 2
                elif s == 2:
                    data[order_shift] = buffer_o | bond >> 1
                    order_shift += 1
                    buffer_o = bond << 7
                    s = 3
                elif s == 3:
                    buffer_o |= bond << 4
                    s = 4
                elif s == 4:
                    buffer_o |= bond << 1
                    s = 5
                elif s == 5:
                    data[order_shift] = buffer_o | bond >> 2
                    order_shift += 1
                    buffer_o = bond << 6
                    s = 6
                elif s == 6:
                    buffer_o |= bond << 3
                    s = 7
                else:  # 7
                    data[order_shift] = buffer_o | bond
                    order_shift += 1
                    s = 0

    if s:  # flush buffer
        data[order_shift] = buffer_o

    for py_tuple, b in py_cis_trans_stereo.items():
        n, m = py_tuple
        data[cis_trans_shift] = n >> 4
        data[cis_trans_shift + 1] = n << 4 | m >> 8
        data[cis_trans_shift + 2] = m
        data[cis_trans_shift + 3] = b
        cis_trans_shift += 4

    try:
        py_pack = data[:size]
    finally:
        PyMem_Free(data)
    return py_pack


cdef short[119] common_isotopes
common_isotopes[:] = [0, -15, -12, -9, -7, -5, -4, -2, 0, 3, 4, 7, 8, 11, 12, 15, 16, 19, 24, 23, 24, 29,
                      32, 35, 36, 39, 40, 43, 43, 48, 49, 54, 57, 59, 63, 64, 68, 69, 72, 73, 75, 77,
                      80, 82, 85, 87, 90, 92, 96, 99, 103, 106, 112, 111, 115, 117, 121, 123, 124, 125,
                      128, 129, 134, 136, 141, 143, 147, 149, 151, 153, 157, 159, 162, 165, 168, 170,
                      174, 176, 179, 181, 185, 188, 191, 193, 193, 194, 206, 207, 210, 211, 216, 215,
                      222, 221, 228, 227, 231, 231, 235, 236, 241, 242, 243, 244, 245, 254, 253, 254,
                      254, 262, 265, 265, 269, 262, 273, 273, 277, 281, 278]


cdef void double_to_float16(double x, unsigned char* p):
    # adopted from cpython source code
    cdef unsigned char sign
    cdef int e
    cdef double f
    cdef unsigned short bits

    if x == 0.:
        p[0] = p[1] = 0
        return

    sign = x < 0.
    if sign:
        x = -x
    f = frexp(x, &e)
    e -= 1
    if f < .5 or f >= 1. or e >= 16 or e < -25:
        p[0] = p[1] = 0
        return  # ignore big values

    f *= 2.0
    if e < -14:
        f = ldexp(f, 14 + e)
        e = 0
    else:
        e += 15
        f -= 1.

    f *= 1024.
    bits = <unsigned short> f | (e << 10) | (sign << 15)
    p[0] = bits >> 8
    p[1] = bits
