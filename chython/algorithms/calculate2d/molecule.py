# -*- coding: utf-8 -*-
#
#  Copyright 2019-2025 Ramil Nugmanov <nougmanoff@protonmail.com>
#  Copyright 2019, 2020 Dinar Batyrshin <batyrshin-dinar@mail.ru>
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
from math import sqrt
from random import random
from typing import TYPE_CHECKING, Union, Dict, Literal
from ...exceptions import ImplementationError
from ...periodictable.base.vector import Vector


if TYPE_CHECKING:
    from chython import MoleculeContainer

try:
    from py_mini_racer import MiniRacer, JSEvalException
    try:
        from importlib.resources import files
    except ImportError:  # python3.8
        from importlib_resources import files

    ctx = MiniRacer()
    ctx.eval('const self = this')
    ctx.eval(files(__package__).joinpath('clean2d.js').read_text())
except RuntimeError:
    ctx = None


class Calculate2DMolecule:
    __slots__ = ()
    _atoms: Dict[int, 'Element']
    _bonds: Dict[int, Dict[int, 'Bond']]

    def clean2d(self: Union['MoleculeContainer', 'Calculate2DMolecule'],
                *, engine: Literal['rdkit', 'smilesdrawer'] = None):
        """
        Calculate 2d layout of graph.

        By default, https://pubs.acs.org/doi/10.1021/acs.jcim.7b00425 JS implementation is used.
        Can be changed globally with the `chython.clean2d_engine` parameter.

        :param engine: override globally set engine
        """
        if engine is None:
            from chython import clean2d_engine as engine

        plane = {}
        if engine == 'rdkit':
            from rdkit.Chem.AllChem import Compute2DCoords

            mol = self.to_rdkit(keep_mapping=False)
            Compute2DCoords(mol)
            # set coordinates from the first rdkit conformer. usually it's 2d layout
            for n, (x, y, _) in zip(self, mol.GetConformers()[0].GetPositions()):
                plane[n] = (x, y)
        elif engine == 'smilesdrawer':
            if ctx is None:
                raise ImportError('mini_racer is not installed or broken')
            entry = iter(sorted(self, key=lambda n: len(self._bonds[n])))
            for _ in range(min(5, len(self))):
                smiles, order = self.__clean2d_prepare(next(entry))
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
        else: raise ValueError(f'Invalid clean2d engine: {engine}')

        bonds = []
        for n, m, _ in self.bonds():
            xn, yn = plane[n]
            xm, ym = plane[m]
            bonds.append(sqrt((xm - xn) ** 2 + (ym - yn) ** 2))
        if bonds:
            bond_reduce = sum(bonds) / len(bonds) / .825
        else:
            bond_reduce = 1.

        atoms = self._atoms
        for n, (x, y) in plane.items():
            atoms[n].xy = (x / bond_reduce,  y / bond_reduce)

        if self.connected_components_count > 1:
            shift_x = 0.
            for c in self.connected_components:
                shift_x = self._fix_plane_mean(shift_x, component=c) + .9
        self.__dict__.pop('__cached_method__repr_svg_', None)

    def _fix_plane_mean(self, shift_x: float, shift_y=0., component=None) -> float:
        atoms = self._atoms
        if component is None:
            component = atoms

        left_atom = atoms[min(component, key=lambda x: atoms[x].x)]
        right_atom = atoms[max(component, key=lambda x: atoms[x].x)]

        min_x = left_atom.x - shift_x
        if len(left_atom.atomic_symbol) == 2:
            min_x -= .2

        max_x = right_atom.x - min_x
        min_y = min(atoms[x].y for x in component)
        max_y = max(atoms[x].y for x in component)
        mean_y = (max_y + min_y) / 2 - shift_y
        delta = Vector(min_x, mean_y)
        for n in component:
            atoms[n].xy -= delta

        if -.18 <= right_atom.y <= .18:
            factor = right_atom.implicit_hydrogens
            if factor == 1:
                max_x += .15
            elif factor:
                max_x += .25
        return max_x

    def _fix_plane_min(self, shift_x: float, shift_y=0., component=None) -> float:
        atoms = self._atoms
        if component is None:
            component = atoms

        right_atom = atoms[max(component, key=lambda x: atoms[x].x)]
        min_x = min(atoms[x].x for x in component) - shift_x
        max_x = right_atom.x - min_x
        min_y = min(atoms[x].y for x in component) - shift_y
        delta = Vector(min_x, min_y)
        for n in component:
            atoms[n].xy -= delta

        if shift_y - .18 <= right_atom.y <= shift_y + .18:
            factor = right_atom.implicit_hydrogens
            if factor == 1:
                max_x += .15
            elif factor:
                max_x += .25
        return max_x

    def __clean2d_prepare(self: 'MoleculeContainer', entry):
        w = {n: random() for n in self._atoms}
        w[entry] = -1
        smiles, order = self._smiles(w.__getitem__, random=True, charges=False, stereo=False, _return_order=True)
        return ''.join(smiles).replace('~', '-'), order


__all__ = ['Calculate2DMolecule']
