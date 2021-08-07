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
from collections import defaultdict
from ._write import _MDLWrite
from ...containers import MoleculeContainer


charge_map = {-4: '  0', -3: '  7', -2: '  6', -1: '  5', 0: '  0', 1: '  3', 2: '  2', 3: '  1', 4: '  0'}


class MDLWrite(_MDLWrite):
    def _convert_molecule(self, g, write3d=None):
        if not isinstance(g, MoleculeContainer):
            raise TypeError('MoleculeContainer expected')
        elif max(g) > 999:
            raise ValueError('MOL file support only small molecules')

        out = [f'{g.name}\n\n\n{g.atoms_count:3d}{g.bonds_count:3d}  0  0  0  0            999 V2000\n']
        if write3d is not None:
            out.append(self.__convert_atoms3d(g, g._conformers[write3d]))
        else:
            out.append(self.__convert_atoms2d(g))
        out.append(self.__convert_bonds(g))
        gc = g._charges
        gr = g._radicals
        for n, (m, a) in enumerate(g._atoms.items(), start=1):
            if a.isotope:
                out.append(f'M  ISO  1 {n:3d} {a.isotope:3d}\n')
            if gr[m]:
                out.append(f'M  RAD  1 {n:3d}   2\n')  # invalid for carbenes
            c = gc[m]
            if c in (-4, 4):
                out.append(f'M  CHG  1 {n:3d} {c:3d}\n')
        out.append('M  END\n')
        return ''.join(out)

    def __convert_atoms2d(self, g):
        gc = g._charges
        gp = g._plane
        out = []
        for n, (m, a) in enumerate(g._atoms.items(), start=1):
            x, y = gp[m]
            c = gc[m]
            if c in (-4, 4):
                if self._mapping:
                    out.append(f'{x:10.4f}{y:10.4f}    0.0000 {a.atomic_symbol:3s} '
                               f'0  0  0  0  0  0  0  0  0{m:3d}  0  0\n')
                else:
                    out.append(f'{x:10.4f}{y:10.4f}    0.0000 {a.atomic_symbol:3s} '
                               f'0  0  0  0  0  0  0  0  0  0  0  0\n')
            else:
                if self._mapping:
                    out.append(f'{x:10.4f}{y:10.4f}    0.0000 {a.atomic_symbol:3s} 0{charge_map[c]}  0  0  0  0'
                               f'  0  0  0{m:3d}  0  0\n')
                else:
                    out.append(f'{x:10.4f}{y:10.4f}    0.0000 {a.atomic_symbol:3s} 0{charge_map[c]}  0  0  0  0'
                               f'  0  0  0  0  0  0\n')
        return ''.join(out)

    def __convert_atoms3d(self, g, xyz):
        gc = g._charges
        out = []
        for n, (m, a) in enumerate(g._atoms.items(), start=1):
            x, y, z = xyz[m]
            c = gc[m]
            if c in (-4, 4):
                if self._mapping:
                    out.append(f'{x:10.4f}{y:10.4f}{z:10.4f} {a.atomic_symbol:3s} '
                               f'0  0  0  0  0  0  0  0  0{m:3d}  0  0\n')
                else:
                    out.append(f'{x:10.4f}{y:10.4f}{z:10.4f} {a.atomic_symbol:3s} 0  0  0  0  0  0  0  0  0  0  0  0\n')
            else:
                if self._mapping:
                    out.append(f'{x:10.4f}{y:10.4f}{z:10.4f} {a.atomic_symbol:3s} 0{charge_map[c]}  0  0  0  0'
                               f'  0  0  0{m:3d}  0  0\n')
                else:
                    out.append(f'{x:10.4f}{y:10.4f}{z:10.4f} {a.atomic_symbol:3s} 0{charge_map[c]}  0  0  0  0'
                               f'  0  0  0  0  0  0\n')
        return ''.join(out)

    @staticmethod
    def __convert_bonds(g):
        bonds = g._bonds
        atoms = {m: n for n, m in enumerate(g._atoms, start=1)}
        wedge = defaultdict(set)
        out = []
        for n, m, s in g._wedge_map:
            out.append(f'{atoms[n]:3d}{atoms[m]:3d}  {bonds[n][m].order}  {s == 1 and "1" or "6"}  0  0  0\n')
            wedge[n].add(m)
            wedge[m].add(n)
        for n, m, b in g.bonds():
            if m not in wedge[n]:
                out.append(f'{atoms[n]:3d}{atoms[m]:3d}  {b.order}  0  0  0  0\n')
        return ''.join(out)


__all__ = ['MDLWrite']
