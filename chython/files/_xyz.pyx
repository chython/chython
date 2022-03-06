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
from libc.math cimport sqrt


cdef extern from "Python.h":
    dict _PyDict_NewPresized(Py_ssize_t minused)


@cython.boundscheck(False)
@cython.wraparound(False)
def possible_bonds(double[:, ::1] xyz not None, double[::1] radii not None, double multiplier):
    cdef unsigned int size, max_bonds, c = 0, b = 0, n, m, k
    cdef double d, nx, ny, nz, rn, mx, my, mz
    cdef dict py_bonds, tmp

    size = xyz.shape[0]
    max_bonds = size * 10  # each atom has less then 10 neighbors approximately

    cdef unsigned int *ns = <unsigned int *> PyMem_Malloc(size * sizeof(unsigned int))
    cdef unsigned int *mi = <unsigned int *> PyMem_Malloc(max_bonds * sizeof(unsigned int))
    cdef double *ds = <double *> PyMem_Malloc(max_bonds * sizeof(double))

    if not ns or not mi or not ds:
        raise MemoryError()

    for n in range(size - 1):
        nx, ny, nz = xyz[n, 0], xyz[n, 1], xyz[n, 2]
        rn = radii[n]
        for m in range(n + 1, size):
            mx, my, mz = nx - xyz[m, 0], ny - xyz[m, 1], nz - xyz[m, 2]
            d = sqrt(mx * mx + my * my + mz * mz)
            if d <= (rn + radii[m]) * multiplier:
                mi[c] = m
                ds[c] = d
                c += 1
        ns[n] = c

    # prepare dict of dicts
    py_bonds = _PyDict_NewPresized(size)
    for n in range(size):
        py_bonds[n + 1] = {}

    for n in range(size - 1):
        c = ns[n]
        n += 1
        tmp = py_bonds[n]
        for m in range(b, c):
            k = mi[m] + 1
            tmp[k] = py_bonds[k][n] = ds[m]
        b = c

    PyMem_Free(ns)
    PyMem_Free(mi)
    PyMem_Free(ds)
    return py_bonds
