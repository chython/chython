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
from math import isnan, nan, radians
from typing import TYPE_CHECKING, Union
from ._templates import rules
from ...exceptions import ImplementationError
from ...periodictable.base.vector import Vector


if TYPE_CHECKING:
    from chython import ReactionContainer, MoleculeContainer


SINGLE = 1
DOUBLE = 2  # double bond
BL = .825
D0 = 0
D30 = radians(30)
D60 = radians(60)
D90 = radians(90)
D120 = radians(120)
D180 = radians(180)
D360 = radians(360)


class Calculate2DMolecule:
    __slots__ = ()

    def clean2d(self: Union['MoleculeContainer', 'Calculate2DMolecule']):
        """
        Calculate 2d layout of graph.
        https://pubs.acs.org/doi/10.1021/acs.jcim.7b00425 JS implementation used as a reference.
        """
        atoms = self._atoms
        bonds = self._bonds
        components = []
        tail = []  # small components pushed to the right for better visuality
        for component in self.connected_components:
            if len(component) == 2:  # 2-atom mols always stored horizontally
                n, m = component
                atoms[n].xy = (0., 0.)
                atoms[m].xy = (BL, 0.)
                tail.insert(0, component)
            elif len(component) > 2:
                components.append(component)
                # mark atoms as non-positioned
                for n in component:
                    atoms[n].xy = (nan, nan)
            else:  # len == 1: just a dot. no need for layout calculation
                tail.append(component)

        if components:
            # preset coordinates with templates
            groups = self._apply_2d_templates()
            # apply KK to fix environment or align groups or process unmatched rings
            for kk in self._kamada_kawai_candidates(groups):
                self._apply_kamada_kawai(kk, [g for g in groups if not g.isdisjoint(kk)])

            for component in components:
                if any(isnan(atoms[n].x) for n in component):
                    self._position_atoms(component)

        components.extend(tail)
        shift_x = 0
        for component in components:
            shift_x = self._fix_plane_mean(shift_x, component=component) + .9
        self.__dict__.pop('__cached_method__repr_svg_', None)

    def _position_atoms(self: 'MoleculeContainer', component):
        atoms = self._atoms
        bonds = self._bonds
        ctc = self._stereo_cis_trans_centers
        ctt = self._stereo_cis_trans_terminals
        cte = self.stereogenic_cis_trans
        # prepare starting points. pick previous atom for angle calculation
        stack = []
        for n in component:
            an = atoms[n]
            if an.in_ring or isnan(an.x):
                continue
            # KK/template layouted non-ring atom
            # we always have at least 1 layouted neighbor
            m = next(m for m in bonds[n] if not isnan(atoms[m].x))
            am = atoms[m]
            angle = am.xy.angle(an.xy)
            # high priority
            stack.append((n, m, angle, -1))
        if not stack:
            # pick any terminal atom. we have tree-like molecule without rings.
            # we for sure have at least 3 atoms in a row, thus, we have to layout at least 1 extra atom.
            for n in component:
                ms = bonds[n]
                if len(ms) == 1:
                    m = next(iter(ms))
                    atoms[n].xy = (0., 0.)
                    atoms[m].xy = Vector(BL, 0).rotate(D30)  # place second atom always top-right
                    stack.append((m, n, D30, -1))
                    break  # 1 is enough

        while stack:
            current, previous, angle, sign = stack.pop()

            env = bonds[current]
            if len(env) == 1:
                # layouting of the branch/molecule is finished
                continue

            ac = atoms[current]

            # chiral cis-trans case.
            if env[previous] == DOUBLE and (b := ctc.get(current)) and (s := self.bond(*b).stereo) is not None:
                # cis-trans case. we came from a cumulene chain, thus, we have one layouted end
                t1, t2 = ts = ctt[current]
                n11, n21, n12, n22 = cte[ts]

                if len(env) == 3:
                    n1, n2 = (n for n in env if n != previous)
                else:  # env == 2
                    n1 = next(n for n in env if n != previous)
                    n2 = None

                if n1 == n11:  # picked 1st atom. no need to switch stereo sigh
                    if not isnan(atoms[n21].x):
                        m = n21  # picked 1st atom. no need to switch stereo sign
                    elif not isnan(atoms[n22].x):  # stereo sign switch
                        m = n22
                        s = not s
                    else:
                        raise ImplementationError
                    counter = t2
                elif n1 == n12:  # picked 2nd atom. stereo sign switch
                    if not isnan(atoms[n21].x):
                        m = n21
                        s = not s
                    elif not isnan(atoms[n22].x):  # picked 2nd atom. double stereo-switch. keep as is.
                        m = n22
                    else:
                        raise ImplementationError
                    counter = t2
                elif n1 == n21:
                    if not isnan(atoms[n11].x):
                        m = n11
                    elif not isnan(atoms[n12].x):
                        m = n12
                        s = not s
                    else:
                        raise ImplementationError
                    counter = t1
                else:
                    if not isnan(atoms[n11].x):
                        m = n11
                        s = not s
                    elif not isnan(atoms[n12].x):
                        m = n12
                    else:
                        raise ImplementationError
                    counter = t1

                vt = atoms[counter].xy
                if (atoms[m].xy - vt) @ (ac.xy - vt) > 0:
                    sign = 1 if s else -1
                else:
                    sign = -1 if s else 1

                angle -= sign * D60
                stack.append((n1, current, angle, sign))
                xy = ac.xy + Vector(BL, 0).rotate(angle)
                an = atoms[n1]
                if not isnan(an.x):
                    # reached layouted fragment. rotate the whole fragment and drop chain.
                    stack.pop()
                    raise NotImplementedError
                else:
                    an.xy = xy
                if n2:
                    angle += sign * D120
                    stack.append((n2, current, angle, -sign))
                    xy = ac.xy + Vector(BL, 0).rotate(angle)
                    an = atoms[n2]
                    if not isnan(an.x):
                        # reached layouted fragment. rotate the whole fragment and drop chain.
                        stack.pop()
                        raise NotImplementedError
                    else:
                        an.xy = xy

            # simple non-chiral cases
            elif len(env) == 2:
                n = next(n for n in env if n != previous)
                if ac.hybridization == 3:
                    # keep the same direction
                    stack.append((n, current, angle, sign))
                else:
                    angle += sign * D60
                    stack.append((n, current, angle, -sign))
                xy = ac.xy + Vector(BL, 0).rotate(angle)

                an = atoms[n]
                if not isnan(an.x):
                    # reached layouted fragment. rotate the whole fragment and drop chain.
                    stack.pop()
                    raise NotImplementedError
                else:
                    an.xy = xy
            elif len(env) == 3:
                n, m = (n for n in env if n != previous)
                # continue to grow to the same direction
                angle += sign * D60
                stack.append((n, current, angle, -sign))
                xy = ac.xy + Vector(BL, 0).rotate(angle)
                an = atoms[n]
                if not isnan(an.x):
                    # reached layouted fragment. rotate the whole fragment and drop chain.
                    stack.pop()
                    raise NotImplementedError
                else:
                    an.xy = xy

                # make a side branch
                angle -= sign * D120
                stack.append((m, current, angle, sign))
                xy = ac.xy + Vector(BL, 0).rotate(angle)
                am = atoms[m]
                if not isnan(am.x):
                    # reached layouted fragment. rotate the whole fragment and drop chain.
                    stack.pop()
                    raise NotImplementedError
                else:
                    am.xy = xy
            else:  # 4+ neighbors. position on circle
                delta = D360 / len(env)
                angle += D180
                for n in env:
                    if n != previous:
                        angle += delta
                        xy = ac.xy + Vector(BL, 0).rotate(angle)
                        stack.append((n, current, angle, sign))  # keep sign to minimize overlaps
                        an = atoms[n]
                        if not isnan(an.x):
                            # reached layouted fragment. rotate the whole fragment and drop chain.
                            stack.pop()
                            raise NotImplementedError
                        else:
                            an.xy = xy

    def _apply_2d_templates(self):
        atoms = self._atoms
        seen = set()
        groups = []
        for q, layout in rules:
            for m in q.get_mapping(self, automorphism_filter=False):
                if not seen.isdisjoint(m.values()):  # avoid any overlap
                    continue
                seen.update(m.values())
                groups.append(set(m.values()))
                for n, xy in zip(m.values(), layout):
                    atoms[n].xy = xy
        return groups

    def _kamada_kawai_candidates(self, groups):
        atoms = self._atoms
        bonds = self._bonds
        clusters = [{n} | bonds[n].keys() for n, a in atoms.items() if a.in_ring]
        clusters.extend(g | {m for n in g for m in bonds[n]} for g in groups)  # add layouted groups
        solved = []
        while clusters:
            c1 = clusters.pop()
            for c2 in clusters:
                if not c1.isdisjoint(c2):
                    c2.update(c1)
                    break
            else:
                if c1 not in groups:
                    solved.append(c1)
        return solved

    def _apply_kamada_kawai(self, system, groups):
        atoms = self._atoms
        bonds = self._bonds

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
            a = atoms[n].xy
            a.x -= min_x
            a.y -= mean_y

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
            a = atoms[n].xy
            a.x -= min_x
            a.y -= min_y

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
