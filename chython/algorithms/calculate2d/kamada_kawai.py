# -*- coding: utf-8 -*-
#
#  Copyright 2025 Ramil Nugmanov <nougmanoff@protonmail.com>
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
import numpy as np
from itertools import combinations
from math import isnan, nan, radians, sin, cos, sqrt, pi
from numpy.linalg import norm
from scipy.optimize import minimize
from scipy.spatial.distance import pdist, squareform
from typing import Set, Dict, TYPE_CHECKING


if TYPE_CHECKING:
    from chython.containers import Bond
    from chython.periodictable import Element


SINGLE = 1
DOUBLE = 2  # double bond
BL = .825
BLC = 1.32
BL25 = BL * 2.5
BL120 = BL * sqrt(3)


class KamadaKawai:
    """
    Based on https://pubs.acs.org/doi/10.1021/acs.jcim.6b00391
    """
    __slots__ = ()
    _atoms: Dict[int, 'Element']
    _bonds: Dict[int, Dict[int, 'Bond']]

    def _apply_kamada_kawai(self, component: Set[int]):
        atoms = self._atoms
        pos0, adj, bond_force, angle_force, rep_force, mapping = self._initialize_kamada_kawai(component)

        defined = ~np.isnan(adj)
        upper = np.triu(np.ones_like(adj, dtype=bool), k=1)
        pair_mask = defined & upper
        rep_mask = upper

        # Avoid diagonal artifacts in repulsion weights
        np.fill_diagonal(rep_force, 0.0)

        def _dist(pos_flat):
            pos = pos_flat.reshape(-1, 3)
            return pos, squareform(pdist(pos))

        def _bond_energy(dist, k_mat):
            # E = 1/2 * k * (r - r0)^2 over selected pairs
            diff = dist[pair_mask] - adj[pair_mask]
            k = k_mat[pair_mask]
            return 0.5 * np.sum(k * diff * diff)

        def _rep_energy(dist, c_rep):
            # If F = c_rep / r, then U = -c_rep * ln(r) (+const)
            r = np.clip(dist, 0.1, None)
            w = rep_force[rep_mask]
            return -c_rep * np.sum(w * np.log(r[rep_mask]))

        def stage1(pos_flat):
            pos, dist = _dist(pos_flat)
            return _bond_energy(dist, bond_force) + _rep_energy(dist, 0.04)

        def stage2(pos_flat):
            pos, dist = _dist(pos_flat)
            return _bond_energy(dist, bond_force) + _rep_energy(dist, 1.0)

        def stage3(pos_flat):
            pos, dist = _dist(pos_flat)
            # flatten: penalize z^2 smoothly (constant force can't be represented as a simple potential)
            z_pen = 0.2 * np.sum(pos[:, 2] * pos[:, 2])
            return _bond_energy(dist, bond_force) + _rep_energy(dist, 1.0) + z_pen

        def stage4(pos_flat):
            pos, dist = _dist(pos_flat)
            z_pen = 0.2 * np.sum(pos[:, 2] * pos[:, 2])
            return _bond_energy(dist, angle_force) + _rep_energy(dist, 1.0) + z_pen

        res = minimize(stage1, pos0, method='L-BFGS-B', options={'maxiter': 1000})
        res = minimize(stage2, res.x, method='L-BFGS-B', options={'maxiter': 1000})
        res = minimize(stage3, res.x, method='L-BFGS-B', options={'maxiter': 1000})
        res = minimize(stage4, res.x, method='L-BFGS-B', options={'maxiter': 1000})

        final_pos = res.x.reshape(-1, 3).tolist()
        for n, i in mapping.items():
            x, y, z = final_pos[i]
            atoms[n].xy = x, y

    def _initialize_kamada_kawai(self, component):
        atoms = self._atoms
        bonds = self._bonds

        mapping = {n: i for i, n in enumerate(component)}
        # add "lone pair" to each atom with two neighbours, except allenes/alkynes
        lone = {n: bs for n in component if len(bs := bonds[n]) == 2 and atoms[n].hybridization != 3}

        size = len(component) + len(lone)
        adj = np.full((size, size), nan)
        bond_force = np.zeros((size, size))
        angle_force = np.zeros((size, size))
        rep_force = np.ones((size, size))

        for i, (n, bs) in enumerate(lone.items(), len(component)):
            # set bond length to the lone pair
            j = mapping[n]
            adj[i, j] = adj[j, i] = BL
            bond_force[i, j] = bond_force[j, i] = .1
            angle_force[i, j] = angle_force[j, i] = .2  # doubled force at stage 4
            m1, m2 = bs
            j1, j2 = mapping[m1], mapping[m2]
            adj[i, j1] = adj[j1, i] = BL120  # set optimal length from LP to neighbours
            angle_force[i, j1] = angle_force[j1, i] = .3
            adj[i, j2] = adj[j2, i] = BL120
            angle_force[i, j2] = angle_force[j2, i] = .3
            # set angle spring between neighbours
            adj[j1, j2] = adj[j2, j1] = BL120
            angle_force[j1, j2] = angle_force[j2, j1] = .3

        for n in component:
            bs = bonds[n]
            if len(bs) == 2 and atoms[n].hybridization == 3:
                # set angle force for alkynes/allenes
                m1, m2 = bs
                j1, j2 = mapping[m1], mapping[m2]
                adj[j1, j2] = adj[j2, j1] = BL25
                angle_force[j1, j2] = angle_force[j2, j1] = .3
            elif len(bs) == 3:
                # set angle force for triangles
                for m1, m2 in combinations(bs, 2):
                    j1, j2 = mapping[m1], mapping[m2]
                    adj[j1, j2] = adj[j2, j1] = BL120
                    angle_force[j1, j2] = angle_force[j2, j1] = .3
            elif len(bs) == 4:
                # set doubled repulsion force for 90C
                for m1, m2 in combinations(bs, 2):
                    j1, j2 = mapping[m1], mapping[m2]
                    rep_force[j1, j2] = rep_force[j2, j1] = 2

        # set bonds. order dependent, to be sure angle constants overridden
        for n in component:
            i = mapping[n]
            bs = bonds[n]
            for m in bs:
                j = mapping[m]
                adj[i, j] = BL
                bond_force[i, j] = .1
                angle_force[i, j] = .2

        # extend hypervalent atoms' bond length
        for n in component:
            bs = bonds[n]
            if len(bs) > 4:
                i = mapping[n]
                for m in bs:
                    j = mapping[m]
                    adj[i, j] = adj[j, i] = BLC

        pos0 = np.random.normal(0., max(3., sqrt(size)) * BL, size=(size, 3)).flatten()
        return pos0, adj, bond_force, angle_force, rep_force, mapping


__all__ = ['KamadaKawai']
