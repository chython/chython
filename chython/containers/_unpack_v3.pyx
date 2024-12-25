# -*- coding: utf-8 -*-
# cython: language_level=3
#
#  Copyright 2023 Ramil Nugmanov <nougmanoff@protonmail.com>
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

from chython.containers.bonds import Bond

# Format specification::
#
# Big endian bytes order
# 8 bit - 0x03 (format specification version)
# Atom block 3 bytes (repeated):
# 1 bit - atom entrance flag (always 1)
# 7 bit - atomic number (<=118)
# 3 bit - hydrogens (0-7). Note: 7 == None
# 4 bit - charge (charge + 4. possible range -4 - 4)
# 1 bit - radical state
# 1 bit padding
# 3 bit tetrahedron/allene sign
#     (000 - not stereo or unknown, 001 - pure-unknown-enantiomer, 010 or 011 - has stereo)
# 4 bit - number of following bonds and CT blocks (0-15)
#
# Bond block 2 bytes (repeated 0-15 times)
# 12 bit - negative shift from current atom to connected (e.g. 0x001 = -1 - connected to previous atom)
# 4 bit - bond order: 0000 - single, 0001 - double, 0010 - triple, 0011 - aromatic, 0111 - special
#
# Cis-Trans 2 bytes
# 12 bit - negative shift from current atom to connected (e.g. 0x001 = -1 - connected to previous atom)
# 4 bit - CT sign: 1000 or 1001 - to avoid overlap with bond

@cython.nonecheck(False)
@cython.boundscheck(False)
@cython.cdivision(True)
@cython.wraparound(False)
def unpack(const unsigned char[::1] data not None):
    cdef char *charges
    cdef unsigned char *atoms, *hydrogens, *radicals, *is_chiral, *neighbors, **orders, *seen
    cdef unsigned short **connections, *ct_stereo
    cdef bint *stereo_sign, *ct_sign

    cdef unsigned char a, b, i
    cdef unsigned short size, shift = 1, n, m, bond_shift, atoms_count, ct_count = 0, ct_shift = 0

    cdef tuple py_xy
    cdef object bond, py_n, py_m
    cdef list py_mapping, py_atoms, py_isotopes, py_bonds_flat
    cdef dict py_charges, py_radicals, py_hydrogens, py_plane, py_bonds, py_ngb
    cdef dict py_atoms_stereo, py_allenes_stereo, py_cis_trans_stereo

    # allocate memory
    size = len(data)
    atoms = <unsigned char*> PyMem_Malloc(size / 3 * sizeof(unsigned char))
    charges = <char*> PyMem_Malloc(size / 3 * sizeof(char))
    radicals = <unsigned char*> PyMem_Malloc(size / 3 * sizeof(unsigned char))
    hydrogens = <unsigned char*> PyMem_Malloc(size / 3 * sizeof(unsigned char))
    is_chiral = <unsigned char*> PyMem_Malloc(size / 3 * sizeof(unsigned char))
    stereo_sign = <bint*> PyMem_Malloc(size / 3 * sizeof(bint))
    ct_stereo = <unsigned short*> PyMem_Malloc(size / 3 * sizeof(unsigned short))
    ct_sign = <bint*> PyMem_Malloc(size / 6 * sizeof(bint))
    seen = <unsigned char*> PyMem_Malloc(size / 3 * sizeof(unsigned char))
    neighbors = <unsigned char *> PyMem_Malloc(size / 3 * sizeof(unsigned char))
    connections = <unsigned short**> PyMem_Malloc(size / 3 * sizeof(unsigned short*))
    orders = <unsigned char**> PyMem_Malloc(size / 3 * sizeof(unsigned char *))
    for n in range(size / 3):
        connections[n] = <unsigned short*> PyMem_Malloc(15 * sizeof(unsigned short))
        orders[n] = <unsigned char*> PyMem_Malloc(15 * sizeof(unsigned char))

    # unpack atom block to separate attributes arrays
    n = 0
    while shift < size:
        seen[n] = 0  # erase randomness
        a = data[shift]
        if a & 0x80 == 0:  # end of pack
            break
        atoms[n] = a & 0x7f

        a = data[shift + 1]
        hydrogens[n] = a >> 5
        charges[n] = ((a >> 1) & 0x0f) - 4
        radicals[n] = a & 0x01

        a = data[shift + 2]
        bond_shift = a & 0x0f
        b = a >> 4
        if b == 0b0011:
            is_chiral[n] = 1
            stereo_sign[n] = True
        elif b == 0b0010:
            is_chiral[n] = 1
            stereo_sign[n] = False
        else:
            is_chiral[n] = 0

        shift += 3
        neighbors[n] = 0
        for i in range(bond_shift):
            a, b = data[shift], data[shift + 1]
            shift += 2

            m = n - (a << 4 | b >> 4)  # second atom index
            b &= 0x0f

            if b < 8:
                connections[n][neighbors[n]] = m
                connections[m][neighbors[m]] = n
                orders[m][neighbors[m]] = b + 1  # only single direction
                neighbors[n] += 1
                neighbors[m] += 1
            else:  # CT stereo
                ct_stereo[ct_shift] = m + 1
                ct_stereo[ct_shift + 1] = n + 1
                ct_sign[ct_count] = b & 0x01
                ct_count += 1
                ct_shift += 2
        n += 1
    atoms_count = n

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
    py_xy = (0., 0.)

    for n in range(atoms_count):
        seen[n] = 1
        py_n = n + 1  # shared py int obj

        # fill intermediate data
        py_mapping.append(py_n)
        py_atoms.append(atoms[n])
        py_isotopes.append(None)

        py_charges[py_n] = charges[n]
        py_radicals[py_n] = bool(radicals[n])
        if hydrogens[n] == 7:
            py_hydrogens[py_n] = None
        else:
            py_hydrogens[py_n] = hydrogens[n]

        py_plane[py_n] = py_xy

        if is_chiral[n]:
            if neighbors[n] == 2:  # allene
                py_allenes_stereo[py_n] = stereo_sign[n]
            else:
                py_atoms_stereo[py_n] = stereo_sign[n]

        py_bonds[py_n] = py_ngb = {}
        for i in range(neighbors[n]):
            m = connections[n][i]
            py_m = m + 1
            if seen[m]:  # bond partially exists. need back-connection.
                py_ngb[py_m] = py_bonds[py_m][py_n]
            else:
                bond = object.__new__(Bond)
                bond._Bond__order = orders[n][i]
                bond._Bond__n = py_n
                bond._Bond__m = py_m
                py_ngb[py_m] = bond
                py_bonds_flat.append(bond)

    ct_shift = 0
    for n in range(ct_count):
        py_cis_trans_stereo[(ct_stereo[ct_shift], ct_stereo[ct_shift + 1])] = ct_sign[n]
        ct_shift += 2

    PyMem_Free(atoms)
    PyMem_Free(charges)
    PyMem_Free(radicals)
    PyMem_Free(hydrogens)
    PyMem_Free(is_chiral)
    PyMem_Free(stereo_sign)
    PyMem_Free(ct_stereo)
    PyMem_Free(ct_sign)
    PyMem_Free(neighbors)
    PyMem_Free(seen)
    for n in range(size / 3):
        PyMem_Free(connections[n])
        PyMem_Free(orders[n])
    PyMem_Free(connections)
    PyMem_Free(orders)

    return (py_mapping, py_atoms, py_isotopes,
            py_charges, py_radicals, py_hydrogens, py_plane, py_bonds,
            py_atoms_stereo, py_allenes_stereo, py_cis_trans_stereo, shift, py_bonds_flat)
