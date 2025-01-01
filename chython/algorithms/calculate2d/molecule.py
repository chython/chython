# -*- coding: utf-8 -*-
#
#  Copyright 2019-2024 Ramil Nugmanov <nougmanoff@protonmail.com>
#  Copyright 2024 Denis Lipatov <denis.lipatov163@gmail.com>
#  Copyright 2024 Vyacheslav Grigorev <slavick2000@yandex.ru>
#  Copyright 2024 Timur Gimadiev <timur.gimadiev@gmail.com>
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
from itertools import combinations
from math import isnan, nan, radians
from numpy import zeros, linspace, column_stack, sin, cos, sqrt, nan_to_num, argmax, errstate
from scipy.sparse.csgraph import shortest_path
from typing import TYPE_CHECKING, Union, Dict
from ._templates import rules
from ...exceptions import ImplementationError
from ...periodictable.base.vector import Vector


if TYPE_CHECKING:
    from chython import MoleculeContainer
    from chython.containers import Bond
    from chython.periodictable import Element


SINGLE = 1
DOUBLE = 2  # double bond
BL = .825
RADIUS = 500
D0 = 0
D30 = radians(30)
D60 = radians(60)
D90 = radians(90)
D120 = radians(120)
D180 = radians(180)
D360 = radians(360)


class Calculate2DMolecule:
    __slots__ = ()
    _atoms: Dict[int, 'Element']
    _bonds: Dict[int, Dict[int, 'Bond']]

    def clean2d(self: Union['MoleculeContainer', 'Calculate2DMolecule'], *,
                kk_outer_iterations: int = 1000, kk_outer_threshold: float =.1,
                kk_inner_iterations: int = 50, kk_inner_threshold: float =.1):
        """
        Calculate 2d layout of graph.
        https://pubs.acs.org/doi/10.1021/acs.jcim.7b00425 JS implementation used as a reference.
        """
        atoms = self._atoms
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
            # apply KK to process bridged rings
            # groups = self._apply_kamada_kawai(groups, kk_outer_iterations, kk_inner_iterations,
            #                                   kk_outer_threshold, kk_inner_threshold)

            for component in components:
                if component not in groups:
                    self._position_atoms(component, [x for x in groups if not component.isdisjoint(x)])

        components.extend(tail)
        shift_x = 0
        for component in components:
            shift_x = self._fix_plane_mean(shift_x, component=component) + .9
        self.__dict__.pop('__cached_method__repr_svg_', None)

    def _initialize_positioning(self, component, groups):
        """
        Prepare starting point and previous atom for angle calculation
        """
        atoms = self._atoms
        bonds = self._bonds

        if not groups:
            # nothing pre-layouted. pick any ring atom if exists or any terminal atom.
            for n in component:
                if atoms[n].in_ring:
                    m = next(m for m in bonds[n] if atoms[m].in_ring)
                    atoms[n].xy = (0., 0.)
                    atoms[m].xy = (0., BL)  # place 1st and second ring atoms always vertical
                    return [(m, n, D90, -1)]

            # pick any terminal atom. we have tree-like molecule without rings.
            # we for sure have at least 3 atoms in a row, thus, we have to layout at least 1 extra atom.
            for n in component:
                ms = bonds[n]
                if len(ms) == 1:
                    m = next(iter(ms))
                    atoms[n].xy = (0., 0.)
                    atoms[m].xy = Vector(BL, 0).rotate(D30)  # place second atom always top-right
                    return [(m, n, D30, -1)]

        # we have pre-layouted groups. let's prepare them for extention.
        stack = []
        layouted = set()
        seen = set()
        for group in groups:
            seen.update(group)
            for n in group:
                env = bonds[n]
                if env.keys() <= seen:  # neighbors already pre-layouted
                    continue
                an = atoms[n]

                # treat atom in group as star where layouted atoms are close to each other and non-layouted on opposite side
                #
                #   L   N
                #    \ /
                # L - L
                #    / \
                #   L   N
                #
                v, c = Vector(0, 0), -1
                for m in env:
                    if m in seen:
                        v += (atoms[m].xy - an.xy).normalise()
                        c += 1
                delta = D360 / len(env)
                angle = v.angle() + delta * c / 2  # ideal position of frontal layouted atom

                for m in env:
                    if m in seen:  # already layouted.
                        continue

                    angle += delta
                    xy = an.xy + Vector(BL, 0).rotate(angle)

                    for other in groups:
                        if not other.isdisjoint(seen):
                            continue
                        elif m in other:
                            # Ring-bond-Ring case
                            # atom is part of another pre-layouted group.
                            # let's reposition the whole group.
                            am, emv = atoms[m], bonds[m]

                            self._rotate_group(other, an.xy, angle - an.xy.angle(am.xy))  # fix angle g-n-m
                            self._shift_group(other, xy - am.xy)  # fix position of m ang the whole group

                            # calculate opposite angle
                            v, c = Vector(0, 0), 1
                            for k in emv:
                                if k in other:
                                    v += (atoms[k].xy - am.xy).normalise()
                                    c += 1
                            # fix angle o-m-n
                            self._rotate_group(other, am.xy, v.angle() + D360 / len(emv) * c / 2)
                            break
                        elif m in layouted:
                            # Ring - Linker Atom - Ring case
                            ...
                            break
                    else:  # non-layouted atom.
                        layouted.add(m)
                        atoms[m].xy = xy
                        stack.append((m, n, angle, -1))
                        continue
        return stack

    def _rotate_group(self, group, point: Vector, angle):
        """
        Rotate the whole group around given point
        """
        atoms = self._atoms

        for n in group:
            an = atoms[n]
            an.xy = an.xy.rotate(angle, point)

    def _shift_group(self, group, shift: Vector):
        """
        Shift the whole group by given vector
        """
        atoms = self._atoms

        for n in group:
            atoms[n].xy += shift

    def _position_atoms(self, component, groups):
        atoms = self._atoms
        bonds = self._bonds
        ctc = self._stereo_cis_trans_centers

        seen = set()
        stack = self._initialize_positioning(component, groups)
        while stack:
            current, previous, angle, sign = stack.pop()
            if current in seen:
                continue
            seen.add(current)

            env = bonds[current]
            if len(env) == 1:
                # layouting of the branch/molecule is finished
                continue

            # chiral cis-trans case.
            if env[previous] == DOUBLE:
                if b := ctc.get(current):
                    if (s := self.bond(*b).stereo) is not None:
                        for n, current, angle, sign, xy in self._position_cis_trans(current, previous, angle, s):
                            an = atoms[n]
                            if not isnan(an.x):
                                # reached layouted fragment. rotate the whole fragment and drop chain.
                                ...
                                # prevent double processing in cases of 2 KK layouted fragments linked by chain
                                seen.add(n)
                                raise NotImplementedError
                            else:
                                an.xy = xy
                                stack.append((n, current, angle, sign))
                        continue

            ac = atoms[current]
            # simple non-chiral cases
            if len(env) == 2:
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
                    seen.add(n)  # prevent double processing in cases of 2 KK layouted fragments linked by chain
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
                    seen.add(n)
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
                    seen.add(n)
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
                            seen.add(n)
                            raise NotImplementedError
                        else:
                            an.xy = xy

    def _position_cis_trans(self, current, previous, angle, s):
        """
        cis-trans case. we came from a cumulene chain, thus, we have one layouted end
        """
        atoms = self._atoms
        ctt = self._stereo_cis_trans_terminals
        cte = self.stereogenic_cis_trans

        t1, t2 = ts = ctt[current]
        n11, n21, n12, n22 = cte[ts]

        ac = atoms[current]
        env = self._bonds[current]
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
        xy = ac.xy + Vector(BL, 0).rotate(angle)
        yield n1, current, angle, sign, xy
        if n2:
            angle += sign * D120
            xy = ac.xy + Vector(BL, 0).rotate(angle)
            yield n2, current, angle, -sign, xy

    def _apply_2d_templates(self):
        """
        Use predefined templates to layout atoms.
        """
        atoms = self._atoms
        seen = set()
        groups = []
        for q, layout in rules:
            for m in q.get_mapping(self, automorphism_filter=False):
                if not seen.isdisjoint(m.values()):  # avoid any overlap
                    continue
                seen.update(m.values())
                groups.append(set(m.values()))
                for i, n in m.items():
                    atoms[n].xy = layout[i - 1]
        return groups

    def _apply_kamada_kawai(self, groups, outer_iterations, inner_iterations, outer_threshold, inner_threshold):
        atoms = self._atoms

        solved = []
        for cluster, length, strength, coordinates, mapping in self._initialize_kamada_kawai(groups):
            pi = -1
            for _ in range(outer_iterations):
                diff = coordinates[:, None, :] - coordinates[None, :, :]  # NxNx2
                sdiff = diff * diff
                energy = diff * (strength * (1 - length / (sqrt(sdiff.sum(-1)) + 1e-5)))[:, :, None]  # NxNx2
                forces = energy.sum(1)  # Nx2
                total = (forces ** 2).sum(-1)  # N

                # pick an atom with the highest force/energy
                i = argmax(total)
                if i == pi:
                    total[i] = 0
                    i = argmax(total)
                pi = i
                if total[i] <= outer_threshold:
                    # if it less than threshold, we have solved system. finish.
                    break

                li = length[i]  # N
                si = strength[i]  # N
                diff_i = diff[i]  # Nx2
                sdiff_i = sdiff[i]  # Nx2
                for _ in range(inner_iterations):
                    norm = li / (sdiff_i.sum(-1) ** 1.5 + 1e-5)
                    dxx, dyy = (si[:, None] * (1 - norm[:, None] * sdiff_i)).sum(0).tolist()
                    dxy = float((si * norm * diff_i.prod(-1)).sum())
                    if abs(dxy) < 0.1:
                        dxy = 0.1 if dxy > 0 else -0.1
                    if abs(dxx) < 0.1:
                        dxx = 0.1 if dxx > 0 else -0.1

                    d_ex, d_ey = forces[i].tolist()
                    dy = (d_ex / dxx + d_ey / dxy) / (dxy / dxx - dyy / dxy)
                    dx = -(dxy * dy + d_ex) / dxx
                    coordinates[i] += (dx, dy)

                    # update forces
                    diff_i = coordinates[i] - coordinates  # Nx2
                    sdiff_i = diff_i * diff_i  # Nx2
                    energy_i = diff_i * (si * (1 - li / (sqrt(sdiff_i.sum(-1)) + 1e-5)))[:, None]  # Nx2
                    forces[i] = energy_i.sum(0)  # 2
                    total[i] = (forces[i] ** 2).sum()  # 1

                    if total[i] <= inner_threshold:
                        # local minima for i-th atom found.
                        break

            for n in cluster:
                atoms[n].xy = coordinates[mapping[n]].tolist()
            solved.append(cluster)
        return solved

    def _initialize_kamada_kawai(self, groups):
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
                if c1 in groups:
                    continue
                solved.append(c1)

        for cluster in solved:
            mapping = {n: i for i, n in enumerate(cluster)}

            adj = zeros((len(cluster), len(cluster)))
            angles = linspace(0, D360, len(cluster) + 1)[:-1]
            coordinates = column_stack([cos(angles), sin(angles)]) * RADIUS

            # create adjacency matrix
            for n in cluster:
                i = mapping[n]
                for m in bonds[n].keys() & cluster:
                    j = mapping[m]
                    adj[i, j] = BL

            layouted = []
            for g in groups:
                if g.isdisjoint(cluster):
                    continue
                # for pre-layouted groups calc pairwise distances
                for n, m in combinations(g, 2):
                    d = atoms[n].xy | atoms[m].xy
                    i, j = mapping[n], mapping[m]
                    adj[i, j] = adj[j, i] = d
                layouted.extend(g)
            length = shortest_path(adj, method='FW', directed=False)
            # originally used BL / (topological distance)**2
            # here distance is already BL scaled: BL**3 / (BL*TD)**2 = BL**3 / BL**2 / TD **2 = BL / TD ** 2
            # but we have prelayouted atoms with bonds != BL. Let's just assume they are close enough.
            # adj magic here to reset strength of layouted groups to actual distances.
            with errstate(divide='ignore'):
                strength = nan_to_num(BL**3 / (length ** 2), posinf=0) * (adj == 0) + adj

            if layouted:
                center = sum((atoms[n].xy for n in layouted), Vector(0, 0)) / len(layouted)
                for n in layouted:
                    coordinates[mapping[n], :] = tuple(atoms[n].xy - center)
            yield cluster, length, strength, coordinates, mapping

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

    def _fix_plane_min(self, shift_x: float, shift_y=0., component=None) -> float:
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


__all__ = ['Calculate2DMolecule']
