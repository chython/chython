# -*- coding: utf-8 -*-
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
from collections import defaultdict
from io import StringIO, TextIOWrapper
from pathlib import Path
from ...containers import MoleculeContainer


charge_map = {-4: '  0', -3: '  7', -2: '  6', -1: '  5', 0: '  0', 1: '  3', 2: '  2', 3: '  1', 4: '  0'}


class IO:
    def __init__(self, file, *, mapping: bool = True, append: bool = False):
        """
        :param mapping: write atom mapping.
        :param append: open file path in append mode.
        """
        self._mapping = mapping

        if isinstance(file, str):
            self._file = open(file, 'a' if append else 'w')
            self._is_buffer = False
        elif isinstance(file, Path):
            self._file = file.open('a' if append else 'w')
            self._is_buffer = False
        elif isinstance(file, (TextIOWrapper, StringIO)):
            self._file = file
            self._is_buffer = True
        else:
            raise TypeError('invalid file. TextIOWrapper, StringIO subclasses or path to file expected')

    def close(self, force=False):
        """
        Close opened file

        :param force: force closing of externally opened file or buffer
        """
        self.write = self.__write_closed

        if not self._is_buffer or force:
            self._file.close()

    def __enter__(self):
        return self

    def __exit__(self, _type, value, traceback):
        self.close()

    @staticmethod
    def __write_closed(_):
        raise ValueError('I/O operation on closed writer')


class EMOLWrite(IO):
    def _write_molecule(self, g, write3d=None):
        if not isinstance(g, MoleculeContainer):
            raise TypeError('MoleculeContainer expected')

        if write3d is not None:
            xyz = g._conformers[write3d]
        else:
            z = 0

        gc = g._charges
        gr = g._radicals
        gp = g._plane
        gb = g._bonds

        file = self._file
        file.write(f'M  V30 BEGIN CTAB\nM  V30 COUNTS {g.atoms_count} {g.bonds_count} 0 0 0\nM  V30 BEGIN ATOM\n')

        for n, (m, a) in enumerate(g._atoms.items(), start=1):
            if write3d is not None:
                x, y, z = xyz[m]
                z = f'{z:.4f}'
            else:
                x, y = gp[m]

            c = gc[m]
            c = f' CHG={c}' if c else ''
            r = ' RAD=2' if gr[m] else ''
            i = f' MASS={a.isotope}' if a.isotope else ''

            if not self._mapping:
                m = 0
            file.write(f'M  V30 {n} {a.atomic_symbol} {x:.4f} {y:.4f} {z} {m}{c}{r}{i}\n')

        file.write('M  V30 END ATOM\nM  V30 BEGIN BOND\n')

        mapping = {m: n for n, m in enumerate(g, start=1)}
        wedge = defaultdict(set)
        i = 0  # trick for empty wedge_map
        for i, (n, m, s) in enumerate(g._wedge_map, start=1):
            file.write(f'M  V30 {i} {gb[n][m].order} {mapping[n]} {mapping[m]} CFG={s == 1 and "1" or "3"}\n')
            wedge[n].add(m)
            wedge[m].add(n)

        for n, m, b in g.bonds():
            if m not in wedge[n]:
                i += 1
                file.write(f'M  V30 {i} {b.order} {mapping[n]} {mapping[m]}\n')
        file.write('M  V30 END BOND\nM  V30 END CTAB\n')


class MOLWrite(IO):
    def _write_molecule(self, g, write3d=None):
        if not isinstance(g, MoleculeContainer):
            raise TypeError('MoleculeContainer expected')
        elif max(g) > 999:
            raise ValueError('MOL file support only small molecules')

        if write3d is not None:
            xyz = g._conformers[write3d]
        else:
            z = 0.

        gc = g._charges
        gr = g._radicals
        gp = g._plane
        gb = g._bonds

        file = self._file
        file.write(f'{g.name}\n\n\n{g.atoms_count:3d}{g.bonds_count:3d}  0  0  0  0            999 V2000\n')

        for n, (m, a) in enumerate(g._atoms.items(), start=1):
            if write3d is not None:
                x, y, z = xyz[m]
            else:
                x, y = gp[m]

            c = charge_map[gc[m]]
            if not self._mapping:
                m = 0
            file.write(f'{x:10.4f}{y:10.4f}{z:10.4f} {a.atomic_symbol:3s} 0{c}  0  0  0  0  0  0  0{m:3d}  0  0\n')

        atoms = {m: n for n, m in enumerate(g._atoms, start=1)}
        wedge = defaultdict(set)
        for n, m, s in g._wedge_map:
            file.write(f'{atoms[n]:3d}{atoms[m]:3d}  {gb[n][m].order}  {s == 1 and "1" or "6"}  0  0  0\n')
            wedge[n].add(m)
            wedge[m].add(n)
        for n, m, b in g.bonds():
            if m not in wedge[n]:
                file.write(f'{atoms[n]:3d}{atoms[m]:3d}  {b.order}  0  0  0  0\n')

        for n, (m, a) in enumerate(g._atoms.items(), start=1):
            if a.isotope:
                file.write(f'M  ISO  1 {n:3d} {a.isotope:3d}\n')
            if gr[m]:
                file.write(f'M  RAD  1 {n:3d}   2\n')  # invalid for carbenes
            c = gc[m]
            if c in (-4, 4):
                file.write(f'M  CHG  1 {n:3d} {c:3d}\n')
        file.write('M  END\n')


__all__ = ['MOLWrite', 'EMOLWrite']
