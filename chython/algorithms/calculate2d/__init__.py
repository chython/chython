# -*- coding: utf-8 -*-
#
#  Copyright 2019-2024 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from typing import TYPE_CHECKING, Union
from ._templates import rules


if TYPE_CHECKING:
    from chython import ReactionContainer, MoleculeContainer


BL = .825


class Calculate2DMolecule:
    __slots__ = ()

    def clean2d(self: Union['MoleculeContainer', 'Calculate2DMolecule']):
        """
        Calculate 2d layout of graph.
        https://pubs.acs.org/doi/10.1021/acs.jcim.7b00425 JS implementation used as a reference.
        """
        shift_x = 0
        groups = None
        for component in self.connected_components:
            if len(component) == 2:  # 2-atom mols always stored horizontally
                n, m = component
                a = self._atoms[n]
                a._x = a._y = 0
                a = self._atoms[m]
                a._x, a._y = BL, 0
            elif len(component) > 2:
                if groups is None:
                    # apply templates with predefined layout: rings and hard cases.
                    groups, seen = self._apply_2d_templates()
                    super_rings = [r for r in self.sssr if any(n not in seen for n in r)]



            # else len == 1: just a dot. no need for layout calculation
            shift_x = self._fix_plane_mean(shift_x, component=component) + .9
        self.__dict__.pop('__cached_method__repr_svg_', None)

    def _apply_2d_templates(self):
        atoms = self._atoms
        seen = set()
        groups = []
        for q, layout in rules:
            for m in q.get_mapping(self, automorphism_filter=False):
                if not seen.isdisjoint(m.values()):  # avoid any overlap. rules preordered from complex to simple
                    continue
                seen.update(m.values())
                groups.append(list(m.values()))
                for n, xy in zip(m.values(), layout):
                    atom = atoms[n]
                    atom._x, atom._y = xy
        return groups, seen

    def _apply_kamada_kawai(self, group):
        pass

    def _fix_plane_mean(self: 'MoleculeContainer', shift_x: float, shift_y=0., component=None) -> float:
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
        for n in component:
            a = atoms[n]
            a._x -= min_x
            a._y -= mean_y

        if -.18 <= right_atom.y <= .18:
            factor = right_atom.implicit_hydrogens
            if factor == 1:
                max_x += .15
            elif factor:
                max_x += .25
        return max_x

    def _fix_plane_min(self: 'MoleculeContainer', shift_x: float, shift_y=0., component=None) -> float:
        atoms = self._atoms
        if component is None:
            component = atoms

        right_atom = atoms[max(component, key=lambda x: atoms[x].x)]
        min_x = min(atoms[x].x for x in component) - shift_x
        max_x = right_atom.x - min_x
        min_y = min(atoms[x].y for x in component) - shift_y

        for n in component:
            a = atoms[n]
            a._x -= min_x
            a._y -= min_y

        if shift_y - .18 <= right_atom.y <= shift_y + .18:
            factor = right_atom.implicit_hydrogens
            if factor == 1:
                max_x += .15
            elif factor:
                max_x += .25
        return max_x


class Calculate2DReaction:
    __slots__ = ()

    def clean2d(self: 'ReactionContainer'):
        """
        Recalculate 2d coordinates
        """
        for m in self.molecules():
            m.clean2d()
        self.fix_positions()

    def fix_positions(self: 'ReactionContainer'):
        """
        Fix coordinates of molecules in reaction
        """
        shift_x = 0
        reactants = self.reactants
        amount = len(reactants) - 1
        signs = []
        for m in reactants:
            max_x = m._fix_plane_mean(shift_x)
            if amount:
                max_x += .2
                signs.append(max_x)
                amount -= 1
            shift_x = max_x + 1
        arrow_min = shift_x

        if self.reagents:
            shift_x += .4
            for m in self.reagents:
                max_x = m._fix_plane_min(shift_x, .5)
                shift_x = max_x + 1
            shift_x += .4
            if shift_x - arrow_min < 3:
                shift_x = arrow_min + 3
        else:
            shift_x += 3
        arrow_max = shift_x - 1

        products = self.products
        amount = len(products) - 1
        for m in products:
            max_x = m._fix_plane_mean(shift_x)
            if amount:
                max_x += .2
                signs.append(max_x)
                amount -= 1
            shift_x = max_x + 1
        self._arrow = (arrow_min, arrow_max)
        self._signs = tuple(signs)
        self.flush_cache()


__all__ = ['Calculate2DMolecule', 'Calculate2DReaction']
