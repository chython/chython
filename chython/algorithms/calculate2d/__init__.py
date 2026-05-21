# -*- coding: utf-8 -*-
#
#  Copyright 2019-2025 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from .molecule import *
from .reaction import *
from math import sqrt
from importlib.resources import files
from random import random
from ...exceptions import ImplementationError
from ...periodictable import Element


class Calculate2D:
    __slots__ = ()

    def _fix_plane_component(self, component):
        plane = self._plane
        atoms = self._atoms
        if component is None:
            component = atoms

        component = tuple(n for n in component if n in atoms)
        for n in component:
            if n not in plane:
                atom = atoms[n]
                try:
                    plane[n] = tuple(atom.xy)
                except (AttributeError, TypeError, KeyError):
                    plane[n] = (0., 0.)
        return component

    def clean2d(self):
        """
        Calculate 2d layout of graph. https://pubs.acs.org/doi/10.1021/acs.jcim.7b00425 JS implementation used.
        """
        if ctx is None:
            raise ImportError('mini-racer is not installed or broken')
        plane = {}
        for _ in range(5):
            smiles, order = self._clean2d_prepare()
            try:
                xy = ctx.call('$.clean2d', smiles)
            except JSEvalException:
                continue
            break
        else:
            raise ImplementationError

        shift_x, shift_y = xy[0]
        for n, (x, y) in zip(order, xy):
            plane[n] = (x - shift_x, shift_y - y)

        bonds = []
        for n, m, _ in self.bonds():
            xn, yn = plane[n]
            xm, ym = plane[m]
            bonds.append(sqrt((xm - xn) ** 2 + (ym - yn) ** 2))
        if bonds:
            bond_reduce = sum(bonds) / len(bonds) / .825
        else:
            bond_reduce = 1.

        self_plane = self._plane
        atoms = self._atoms
        for n, (x, y) in plane.items():
            xy = (x / bond_reduce, y / bond_reduce)
            self_plane[n] = xy
            if n in atoms and hasattr(atoms[n], 'xy'):
                atoms[n].xy = xy

        self.__dict__.pop('__cached_method__repr_svg_', None)

    def _fix_plane_mean(self, shift_x, shift_y=0, component=None):
        plane = self._plane
        atoms = self._atoms
        component = self._fix_plane_component(component)
        if not component:
            return shift_x

        left_atom = min(component, key=lambda x: plane[x][0])
        right_atom = max(component, key=lambda x: plane[x][0])

        min_x = plane[left_atom][0] - shift_x
        if len(atoms[left_atom].atomic_symbol) == 2:
            min_x -= .2

        max_x = plane[right_atom][0] - min_x
        min_y = min(plane[x][1] for x in component)
        max_y = max(plane[x][1] for x in component)
        mean_y = (max_y + min_y) / 2 - shift_y
        for n in component:
            x, y = plane[n]
            xy = (x - min_x, y - mean_y)
            plane[n] = xy
            if n in atoms and hasattr(atoms[n], 'xy'):
                atoms[n].xy = xy

        return max_x

    def _fix_plane_min(self, shift_x, shift_y=0, component=None):
        plane = self._plane
        atoms = self._atoms
        component = self._fix_plane_component(component)
        if not component:
            return shift_x

        right_atom = max(component, key=lambda x: plane[x][0])
        min_x = min(plane[x][0] for x in component) - shift_x
        max_x = plane[right_atom][0] - min_x
        min_y = min(plane[x][1] for x in component) - shift_y

        for n in component:
            x, y = plane[n]
            xy = (x - min_x, y - min_y)
            plane[n] = xy
            if n in atoms and hasattr(atoms[n], 'xy'):
                atoms[n].xy = xy

        return max_x


class Calculate2DQuery(Calculate2D):
    __slots__ = ()

    def _clean2d_prepare(self):
        from ...containers.molecule import MoleculeContainer
        mol = MoleculeContainer()
        for n, atom in self._atoms.items():
            atom = Element.from_atomic_number(atom.atomic_number or 6)()
            mol.add_atom(atom, n)
        for n, m, bond in self.bonds():
            mol.add_bond(n, m, bond.order[0])
        for n in mol._atoms:
            mol._atoms[n]._implicit_hydrogens = 0
        smiles, order = mol._smiles(lambda x: random(), _return_order=True)
        return ''.join(smiles).replace('~', '-'), order


class Calculate2DCGR(Calculate2D):
    __slots__ = ()

    def _clean2d_prepare(self):
        from ...containers.molecule import MoleculeContainer
        mol = MoleculeContainer()
        for n, atom in self._atoms.items():
            atom = Element.from_atomic_number(atom.atomic_number or 6)()
            mol.add_atom(atom, n)
        for n, m, bond in self.bonds():
            mol.add_bond(n, m, bond.order or 1)
        for n in mol._atoms:
            mol._atoms[n]._implicit_hydrogens = 0
        smiles, order = mol._smiles(lambda x: random(), _return_order=True)
        return ''.join(smiles).replace('~', '-'), order


ctx = None

try:
    from py_mini_racer import MiniRacer, JSEvalException
except (ImportError, RuntimeError):
    try:
        from mini_racer import MiniRacer, JSEvalException
    except (ImportError, RuntimeError):
        MiniRacer = None
        JSEvalException = Exception

if MiniRacer is not None:
    try:
        ctx = MiniRacer()
        ctx.eval('const self = this')
        ctx.eval(files(__package__).joinpath('clean2d.js').read_text())
    except RuntimeError:
        ctx = None


__all__ = ['Calculate2DMolecule', 'Calculate2DReaction', 'Calculate2DCGR', 'Calculate2DQuery']
