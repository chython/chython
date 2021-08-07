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


class EMDLWrite(_MDLWrite):
    def _convert_molecule(self, g, write3d=None):
        if not isinstance(g, MoleculeContainer):
            raise TypeError('MoleculeContainer expected')

        head = f'M  V30 BEGIN CTAB\nM  V30 COUNTS {g.atoms_count} {g.bonds_count} 0 0 0'

        bond = ['M  V30 BEGIN BOND']
        mapping = {m: n for n, m in enumerate(g, start=1)}
        wedge = defaultdict(set)
        bonds = g._bonds
        i = 0
        for i, (n, m, s) in enumerate(g._wedge_map, start=1):
            bond.append(f'M  V30 {i} {bonds[n][m].order} {mapping[n]} {mapping[m]} CFG={s == 1 and "1" or "3"}')
            wedge[n].add(m)
            wedge[m].add(n)

        for i, (n, m, b) in enumerate(g.bonds(), start=i + 1):
            if m not in wedge[n]:
                bond.append(f'M  V30 {i} {b.order} {mapping[n]} {mapping[m]}')
        bond.append('M  V30 END BOND\nM  V30 END CTAB\n')
        bond = '\n'.join(bond)

        if write3d is not None:
            return '\n'.join((head, self.__convert_atoms3d(g, g._conformers[write3d]), bond))
        else:
            return '\n'.join((head, self.__convert_atoms2d(g), bond))

    def __convert_atoms2d(self, g):
        gc = g._charges
        gr = g._radicals
        gp = g._plane

        out = ['M  V30 BEGIN ATOM']
        for n, (m, a) in enumerate(g._atoms.items(), start=1):
            x, y = gp[m]
            c = gc[m]
            c = f' CHG={c}' if c else ''
            r = ' RAD=2' if gr[m] else ''
            i = f' MASS={a.isotope}' if a.isotope else ''

            if self._mapping:
                out.append(f'M  V30 {n} {a.atomic_symbol} {x:.4f} {y:.4f} 0 {m}{c}{r}{i}')
            else:
                out.append(f'M  V30 {n} {a.atomic_symbol} {x:.4f} {y:.4f} 0 0{c}{r}{i}')
        out.append('M  V30 END ATOM')
        return '\n'.join(out)

    def __convert_atoms3d(self, g, xyz):
        gc = g._charges
        gr = g._radicals

        out = ['M  V30 BEGIN ATOM']
        for n, (m, a) in enumerate(g._atoms.items(), start=1):
            x, y, z = xyz[m]
            c = gc[m]
            c = f' CHG={c}' if c else ''
            r = ' RAD=2' if gr[m] else ''
            i = f' MASS={a.isotope}' if a.isotope else ''

            if self._mapping:
                out.append(f'M  V30 {n} {a.atomic_symbol} {x:.4f} {y:.4f} {z:.4f} {m}{c}{r}{i}')
            else:
                out.append(f'M  V30 {n} {a.atomic_symbol} {x:.4f} {y:.4f} {z:.4f} 0{c}{r}{i}')
        out.append('M  V30 END ATOM')
        return '\n'.join(out)


__all__ = ['EMDLWrite']
