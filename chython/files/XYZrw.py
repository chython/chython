# -*- coding: utf-8 -*-
#
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
from fileinput import FileInput
from io import StringIO, TextIOWrapper
from math import sqrt
from numpy import array, uint16, empty
from numba import njit, f8, u2, u4
from numba.core.types import Tuple as nTuple
from pathlib import Path
from traceback import format_exc
from typing import List, Iterable, Tuple, Optional
from ._mdl import parse_error
from ..containers import MoleculeContainer
from ..periodictable import Element


@njit(nTuple((u2[:, :], f8[:]))(f8[:, :], f8[:], f8),
      {'size': u2, 'max_bonds': u4, 'c': u4, 'n': u2, 'm': u2, 'rn': f8, 'r': f8, 'd': f8,
       'nx': f8, 'ny': f8, 'nz': f8, 'mx': f8, 'my': f8, 'mz': f8}, cache=True)
def _get_possible_bonds(xyz, radii, multiplier):
    size = len(xyz)
    max_bonds = size * 10  # each atom has less then 10 neighbors approximately
    nm = empty((max_bonds, 2), dtype=uint16)
    ds = empty(max_bonds)
    c = 0
    for n in range(size - 1):
        nx, ny, nz = xyz[n]
        rn = radii[n]
        for m in range(n + 1, size):
            mx, my, mz = xyz[m]
            d = sqrt((nx - mx) ** 2 + (ny - my) ** 2 + (nz - mz) ** 2)
            r = (rn + radii[m]) * multiplier
            if d <= r:
                nm[c] = n + 1, m + 1
                ds[c] = d
                c += 1
    return nm[:c], ds[:c]


def get_possible_bonds(atoms, conformer, multiplier):
    possible_bonds = {n: {} for n in atoms}  # distance matrix
    radii = array([a.atomic_radius for a in atoms.values()])
    xyz = array(list(conformer.values()))
    nm, ds = _get_possible_bonds(xyz, radii, multiplier)
    for (n, m), d in zip(nm.tolist(), ds.tolist()):
        possible_bonds[n][m] = possible_bonds[m][n] = d
    return possible_bonds


class XYZ:
    """
    Override class below then inheritance used.
    """
    MoleculeContainer = MoleculeContainer

    def __init__(self, radius_multiplier=1.25, store_log=False):
        """
        :param radius_multiplier: Multiplier of sum of covalent radii of atoms which has bonds
        :param store_log: Store parser log if exists messages to `.meta` by key `ParserLog`.
        """
        self.__radius = radius_multiplier
        self._store_log = store_log
        self._log_buffer = []

    def _info(self, msg):
        self._log_buffer.append(msg)

    def _flush_log(self):
        self._log_buffer.clear()

    def _format_log(self):
        log = '\n'.join(self._log_buffer)
        self._log_buffer.clear()
        return log

    def close(self, force=False):
        """
        Close opened file

        :param force: Force closing of externally opened file or buffer
        """
        if not self._is_buffer or force:
            self._file.close()

    def __enter__(self):
        return self

    def __exit__(self, _type, value, traceback):
        self.close()

    def read(self) -> List[MoleculeContainer]:
        """
        Parse whole file

        :return: list of parsed molecules
        """
        return list(iter(self))

    def __iter__(self) -> Iterable[MoleculeContainer]:
        return (x for x in self._data if not isinstance(x, parse_error))

    def __next__(self) -> MoleculeContainer:
        return next(iter(self))

    def _convert_molecule(self, matrix: Iterable[Tuple[str, Optional[int], float, float, float]], charge=0, radical=0):
        mol = self.MoleculeContainer()
        conformer = {}

        atoms = mol._atoms
        plane = mol._plane
        charges = {}
        for n, (a, c, x, y, z) in enumerate(matrix, 1):
            atoms[n] = atom = Element.from_symbol(a)()
            atom._attach_to_graph(mol, n)
            plane[n] = (x, y)
            conformer[n] = (x, y, z)
            charges[n] = c
        if all(x is not None for x in charges.values()):
            charge = sum(charges.values())
        bonds = get_possible_bonds(atoms, conformer, self.__radius)

        self._log_buffer.extend(mol.saturate(bonds, expected_charge=charge, expected_radicals_count=radical,
                                             logging=True))
        mol._conformers.append(conformer)
        if self._store_log:
            log = self._format_log()
            if log:
                mol.meta['ParserLog'] = log
        return mol


class XYZRead(XYZ):
    """XYZ files reader. Works similar to opened file object. Support `with` context manager.
    On initialization accept opened in text mode file, string path to file,
    pathlib.Path object or another buffered reader object.

    Supported multiple structures in same file. In second line possible to store total charge of system. Example::

        2
        charge=-1
        O 0.0 0.0 0.0
        H 1.0 0.0 0.0

    """
    def __init__(self, file, **kwargs):
        if isinstance(file, str):
            self._file = open(file)
            self._is_buffer = False
        elif isinstance(file, Path):
            self._file = file.open()
            self._is_buffer = False
        elif isinstance(file, (TextIOWrapper, StringIO, FileInput)):
            self._file = file
            self._is_buffer = True
        else:
            raise TypeError('invalid file. TextIOWrapper, StringIO subclasses possible')
        super().__init__(**kwargs)
        self.__file = iter(self._file.readline, '')
        self._data = self.__reader()

    def __reader(self):
        failkey = True
        meta = False
        xyz = charge = size = radical = None
        file = self._file
        try:
            seekable = file.seekable()
        except AttributeError:
            seekable = False
        pos = None
        count = -1
        for n, line in enumerate(self.__file):
            if failkey:
                try:
                    size = int(line)
                except ValueError:
                    pass
                else:
                    if seekable:
                        pos = file.tell() - len(line)
                    count += 1
                    meta = True
                    failkey = False
                    xyz = []
            elif meta:  # second line
                charge = 0
                radical = 0
                for x in line.split():
                    if x.startswith('charge='):
                        try:
                            charge = int(line[7:])
                        except ValueError:
                            failkey = True
                            self._info(f'Line [{n}]: consist errors in charge definition')
                            yield parse_error(count, pos, self._format_log(), line)
                            break
                    elif x.startswith('radical='):
                        try:
                            radical = int(line[8:])
                        except ValueError:
                            failkey = True
                            self._info(f'Line [{n}]: consist errors in radical atoms count definition')
                            yield parse_error(count, pos, self._format_log(), line)
                            break
                else:
                    meta = False
            elif len(xyz) < size:  # XYZ block
                try:
                    symbol, x, y, z = line.split()
                    xyz.append((symbol, None, float(x), float(y), float(z)))
                except ValueError:
                    failkey = True
                    self._info(f'Line [{n}]: consist errors in xyz atom coordinates definition')
                    yield parse_error(count, pos, self._format_log(), line)
                else:
                    if len(xyz) == size:
                        try:
                            container = self._convert_molecule(xyz, charge, radical)
                        except ValueError:
                            self._info(f'record consist errors:\n{format_exc()}')
                            yield parse_error(count, pos, self._format_log(), None)
                        else:
                            yield container
                        failkey = True  # trigger end of XYZ
        if not failkey:  # cut XYZ
            self._info('Last structure not finished')
            yield parse_error(count, pos, self._format_log(), {})

    def parse(self, matrix: Iterable[Tuple[str, float, float, float]], charge: int = 0, radical: int = 0) -> \
            Optional[MoleculeContainer]:
        try:
            container = self._convert_molecule([(e, None, x, y, z) for e, x, y, z in matrix], charge, radical)
        except ValueError:
            self._flush_log()
        else:
            return container

    @classmethod
    def create_parser(cls, *args, **kwargs):
        """
        Create XYZ parser function configured same as XYZRead object.
        """
        obj = object.__new__(cls)
        super(XYZRead, obj).__init__(*args, **kwargs)
        return obj.parse


__all__ = ['XYZRead']
