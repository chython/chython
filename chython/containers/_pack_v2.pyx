# -*- coding: utf-8 -*-
#
#  Copyright 2022-2025 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from libc.string cimport memset

# Format V2 specification::
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
def pack(object molecule):
    cdef bint b = True # binary flag
    cdef char charge
    cdef unsigned char atomic_number, ngb_count, isotope, bond, s = 0, buffer_b, buffer_o, stereo, hcr
    cdef unsigned char *data
    cdef unsigned short atoms_count, bonds_count = 0, cis_trans_count, n, m, tn, tm
    cdef unsigned int size, atoms_shift = 4, bonds_shift, order_shift, cis_trans_shift  # can be > 2^16
    cdef unsigned char[4096] seen  # atom number is 12 bit, thus, can be any value up to 4095. numbers are not continuous

    cdef bytes py_pack
    cdef dict py_ngb, py_atoms, py_bonds, py_stereo
    cdef tuple py_tuple
    cdef object py_atom, py_bond, py_nan_int, py_obj

    # map molecule to vars
    py_atoms = molecule._atoms
    py_bonds = molecule._bonds
    py_stereo = molecule._stereo_cis_trans_terminals

    # calculate elements count
    atoms_count = len(py_atoms)
    cis_trans_count = molecule._cis_trans_count

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

    memset(seen, 0, 4096 * sizeof(unsigned char))  # erase random data

    # start pack collection
    data[0] = 2  # header. specification version 2
    data[1] = atoms_count >> 4  # 5-12b of atom count value
    data[2] = atoms_count << 4 | cis_trans_count >> 8  # 1-4b of atom count value, 9-12b of cis-trans count value
    data[3] = cis_trans_count  # 1-8b of cis-trans count value

    for py_obj, py_atom in py_atoms.items():
        py_ngb = py_bonds[py_obj]
        ngb_count = len(py_ngb)
        n = py_obj  # cast to C
        seen[n] = 1
        atomic_number = py_atom.atomic_number

        py_nan_int = py_atom._isotope  # direct access
        if py_nan_int is None:
            isotope = 0
        else:
            isotope = <short> py_nan_int - common_isotopes[atomic_number]

        py_nan_int = py_atom._stereo
        if py_nan_int is None:
            stereo = 0
        # V2 specification
        # 2 bit tetrahedron | 2 bit allene | 0000
        elif py_nan_int:
            if ngb_count == 2:  # allene
                stereo = 0x30
            else:
                stereo = 0xc0
        else:
            if ngb_count == 2:  # allene
                stereo = 0x20
            else:
                stereo = 0x80

        # precalculate atom attrs
        # should be done independently, due to possible randomness in dicts order.
        # 3 bit - hydrogens (0-7) | 4 bit - charge | 1 bit - radical
        py_nan_int = py_atom._implicit_hydrogens
        if py_nan_int is None:
            hcr = 0xe0  # 0b11100000
        else:
            hcr = <unsigned char> py_nan_int << 5

        charge = py_atom._charge
        hcr |= (charge + 4) << 1
        if py_atom._is_radical:
            hcr |= 1

        data[atoms_shift] = n >> 4  # 5-12b AN
        data[atoms_shift + 1] = n << 4 | ngb_count  # 1-4b AN, 4b NC
        data[atoms_shift + 2] = stereo | isotope >> 1  # TS , AS , 4b I
        data[atoms_shift + 3] = isotope << 7 | atomic_number  # 1bI , A

        # 2 float16 big endian
        double_to_float16(py_atom.x, &data[atoms_shift + 4])
        double_to_float16(py_atom.y, &data[atoms_shift + 6])

        data[atoms_shift + 8] = hcr
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
                bond = <unsigned char> py_bond._order - 1
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

                py_nan_int = py_bond._stereo
                if py_nan_int is not None:
                    py_tuple = py_stereo[py_obj]
                    tn, tm = py_tuple
                    data[cis_trans_shift] = tn >> 4
                    data[cis_trans_shift + 1] = tn << 4 | tm >> 8
                    data[cis_trans_shift + 2] = tm
                    data[cis_trans_shift + 3] = py_nan_int
                    cis_trans_shift += 4

    if s:  # flush buffer
        data[order_shift] = buffer_o

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
