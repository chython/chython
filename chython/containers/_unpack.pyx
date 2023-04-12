# -*- coding: utf-8 -*-
# cython: language_level=3
#
#  Copyright 2021-2023 Ramil Nugmanov <nougmanoff@protonmail.com>
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

from chython.containers.bonds import Bond

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
    cdef char *charges
    cdef unsigned char a, b, c, d, isotope, atomic_number, neighbors_count, s = 0, nc, version
    cdef unsigned char *atoms, *hydrogens, *neighbors, *orders, *is_tet, *is_all
    cdef bint *stereo_sign, *ct_sign, *radicals
    cdef unsigned short atoms_count, bonds_count = 0, cis_trans_count, order_count
    cdef unsigned short i, j, k = 0, n, m, buffer_b, shift = 0
    cdef unsigned short *mapping, *isotopes, *cis_trans_1, *cis_trans_2, *connections
    cdef unsigned int size, atoms_shift = 4, bonds_shift, order_shift, cis_trans_shift
    cdef double *x_coord, *y_coord
    cdef unsigned char[4096] seen

    cdef object bond, py_n, py_m
    cdef dict py_charges, py_radicals, py_hydrogens, py_plane, py_bonds, py_ngb
    cdef dict py_atoms_stereo, py_allenes_stereo, py_cis_trans_stereo
    cdef list py_mapping, py_atoms, py_isotopes, py_bonds_flat

    # read header
    version = data[0]
    a, b, c = data[1], data[2], data[3]
    atoms_count = a << 4| b >> 4
    cis_trans_count = (b & 0x0f) << 8 | c

    # allocate memory
    charges = <char*> PyMem_Malloc(atoms_count * sizeof(char))
    radicals = <bint*> PyMem_Malloc(atoms_count * sizeof(bint))
    atoms = <unsigned char*> PyMem_Malloc(atoms_count * sizeof(unsigned char))
    hydrogens = <unsigned char*> PyMem_Malloc(atoms_count * sizeof(unsigned char))
    neighbors = <unsigned char*> PyMem_Malloc(atoms_count * sizeof(unsigned char))
    is_tet = <unsigned char*> PyMem_Malloc(atoms_count * sizeof(unsigned char))
    is_all = <unsigned char*> PyMem_Malloc(atoms_count * sizeof(unsigned char))
    stereo_sign = <bint*> PyMem_Malloc(atoms_count * sizeof(bint))
    mapping = <unsigned short*> PyMem_Malloc(atoms_count * sizeof(unsigned short))
    isotopes = <unsigned short*> PyMem_Malloc(atoms_count * sizeof(unsigned short))
    x_coord = <double*> PyMem_Malloc(atoms_count * sizeof(double))
    y_coord = <double*> PyMem_Malloc(atoms_count * sizeof(double))

    if not charges or not radicals or not atoms or not hydrogens or not neighbors or not is_tet or not is_all:
        raise MemoryError()
    if not stereo_sign or not mapping or not isotopes or not x_coord or not y_coord:
        raise MemoryError()

    # unpack atom block to separate attributes arrays
    for i in range(atoms_count):
        a, b = data[atoms_shift], data[atoms_shift + 1]
        mapping[i] = n = a << 4 | b >> 4
        seen[n] = 0  # erase random value
        neighbors[i] = neighbors_count = b & 0x0f
        bonds_count += neighbors_count

        a, b = data[atoms_shift + 2], data[atoms_shift + 3]
        if a >> 7:  # tetrahedron bit set
            is_tet[i] = 1
            is_all[i] = 0
            stereo_sign[i] = a & 0x40  # mask th bit
        else:
            is_tet[i] = 0
            if a >> 5:  # allene bit set
                is_all[i] = 1
                stereo_sign[i] = a & 0x10  # mask al bit
            else:
                is_all[i] = 0

        atoms[i] = atomic_number = b & 0x7f
        isotope = (a & 0x0f) << 1 | b >> 7
        if isotope:
            isotopes[i] = common_isotopes[atomic_number] + isotope
        else:
            isotopes[i] = 0

        a, b = data[atoms_shift + 4], data[atoms_shift + 5]
        x_coord[i] = double_from_bytes(a, b)
        a, b = data[atoms_shift + 6], data[atoms_shift + 7]
        y_coord[i] = double_from_bytes(a, b)

        a = data[atoms_shift + 8]
        hydrogens[i] = a >> 5
        charges[i] = ((a >> 1) & 0x0f) - 4
        radicals[i] = a & 0x01
        atoms_shift += 9

    # calculate bonds count and pack sections
    bonds_count /= 2

    if version == 2:
        order_count = bonds_count * 3
        if order_count % 8:
            order_count = order_count / 8 + 1
        else:
            order_count /= 8
    elif version == 0:
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
        elif version == 0:
            for j in range(order_shift, cis_trans_shift, 2):
                # 0 3 3 1 | 2 3 3
                a, b = data[j], data[j + 1]
                orders[i] = a >> 4
                orders[i + 1] = (a >> 1) & 0x7
                orders[i + 2] = (a & 0x1) << 2 | b >> 6
                orders[i + 3] = (b >> 3) & 0x7
                orders[i + 4] = b & 0x7
                i += 5

        if cis_trans_count:
            cis_trans_1 = <unsigned short *> PyMem_Malloc(cis_trans_count * sizeof(unsigned short))
            cis_trans_2 = <unsigned short *> PyMem_Malloc(cis_trans_count * sizeof(unsigned short))
            ct_sign = <bint *> PyMem_Malloc(cis_trans_count * sizeof(bint))
            if not cis_trans_1 or not cis_trans_2 or not ct_sign:
                raise MemoryError()

            for i in range(cis_trans_count):
                a, b = data[cis_trans_shift], data[cis_trans_shift + 1]
                c, d = data[cis_trans_shift + 2], data[cis_trans_shift + 3]
                cis_trans_1[i] = a << 4 | b >> 4
                cis_trans_2[i] = (b & 0x0f) << 8 | c
                ct_sign[i] = d  # d = 0x01 or 0x00
                cis_trans_shift += 4

    # define returned data
    py_mapping = []
    py_atoms = []
    py_isotopes = []
    py_charges = {}
    py_radicals = {}
    py_hydrogens = {}
    py_plane = {}
    py_atoms_stereo = {}
    py_allenes_stereo = {}
    py_cis_trans_stereo = {}
    py_bonds = {}
    py_bonds_flat = []

    for i in range(atoms_count):
        n = mapping[i]
        py_n = n  # shared py int obj

        # fill intermediate data
        py_mapping.append(py_n)
        py_atoms.append(atoms[i])
        py_isotopes.append(isotopes[i] or None)

        py_charges[py_n] = charges[i]
        py_radicals[py_n] = radicals[i]
        if hydrogens[i] == 7:
            py_hydrogens[py_n] = None
        else:
            py_hydrogens[py_n] = hydrogens[i]

        py_plane[py_n] = (x_coord[i], y_coord[i])

        if is_tet[i]:
            py_atoms_stereo[py_n] = stereo_sign[i]
        elif is_all[i]:
            py_allenes_stereo[py_n] = stereo_sign[i]

        py_bonds[py_n] = py_ngb = {}
        seen[n] = 1

        nc = neighbors[i]
        for j in range(shift, shift + nc):
            m = connections[j]
            py_m = m
            if seen[m]:  # bond partially exists. need back-connection.
                py_ngb[py_m] = py_bonds[py_m][py_n]
            else:
                bond = object.__new__(Bond)
                bond._Bond__order = orders[k] + 1
                bond._Bond__n = py_n
                bond._Bond__m = py_m
                py_ngb[py_m] = bond
                py_bonds_flat.append(bond)
                k += 1
        shift += nc

    for i in range(cis_trans_count):
        py_cis_trans_stereo[(cis_trans_1[i], cis_trans_2[i])] = ct_sign[i]

    PyMem_Free(charges)
    PyMem_Free(radicals)
    PyMem_Free(atoms)
    PyMem_Free(hydrogens)
    PyMem_Free(neighbors)
    PyMem_Free(is_tet)
    PyMem_Free(is_all)
    PyMem_Free(stereo_sign)
    PyMem_Free(mapping)
    PyMem_Free(isotopes)
    PyMem_Free(x_coord)
    PyMem_Free(y_coord)
    if bonds_count:
        PyMem_Free(connections)
        PyMem_Free(orders)
        if cis_trans_count:
            PyMem_Free(cis_trans_1)
            PyMem_Free(cis_trans_2)
            PyMem_Free(ct_sign)
    return (py_mapping, py_atoms, py_isotopes,
            py_charges, py_radicals, py_hydrogens, py_plane, py_bonds,
            py_atoms_stereo, py_allenes_stereo, py_cis_trans_stereo, size, py_bonds_flat)


cdef short[119] common_isotopes
common_isotopes[:] = [0, -15, -12, -9, -7, -5, -4, -2, 0, 3, 4, 7, 8, 11, 12, 15, 16, 19, 24, 23, 24, 29,
                      32, 35, 36, 39, 40, 43, 43, 48, 49, 54, 57, 59, 63, 64, 68, 69, 72, 73, 75, 77,
                      80, 82, 85, 87, 90, 92, 96, 99, 103, 106, 112, 111, 115, 117, 121, 123, 124, 125,
                      128, 129, 134, 136, 141, 143, 147, 149, 151, 153, 157, 159, 162, 165, 168, 170,
                      174, 176, 179, 181, 185, 188, 191, 193, 193, 194, 206, 207, 210, 211, 216, 215,
                      222, 221, 228, 227, 231, 231, 235, 236, 241, 242, 243, 244, 245, 254, 253, 254,
                      254, 262, 265, 265, 269, 262, 273, 273, 277, 281, 278]


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
