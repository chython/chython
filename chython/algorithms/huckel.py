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
from functools import cached_property
from numpy import zeros
from numpy.linalg import eig
from typing import TYPE_CHECKING


if TYPE_CHECKING:
    from chython import MoleculeContainer


# atom, charge, radical, non sp3 : h, k, ne
basis = {(5, 0, False, False): (-1., .75, 0),  # X3B
         (6, 0, False, True): (0., 1., 1),
         (7, 0, False, True): (.5, 1., 1),  # pyridine
         (7, 0, False, False): (1.49, .8, 2),  # pyrole
         (7, 1, False, True): (2.3, .7, 1),  # nitro
         (8, 0, False, True): (1., 1., 1),  # X=O
         (8, 0, False, False): (2., .8, 2),  # X-O-Y
         (8, -1, False, False): (2., .8, 2),  # X-[O-]
         (16, 0, False, True): (.4, 1., 1),  # X=S
         (16, 0, False, False): (1.3, .6, 2),  # X-S-Y
         (16, -1, False, False): (1.3, .6, 2),  # X-[S-]
         }


class Huckel:
    __slots__ = ()

    @cached_property
    def huckel_energy(self: 'MoleculeContainer') -> float:
        """
        Huckel method based Pi electrons energy calculator.
        Parametrized for B C N O S.
        """
        hyb = self.hybridization
        charge = self._charges
        radical = self._radicals

        adj = {}
        alpha = {}
        beta = {}
        electrons = {}
        for n, a in self._atoms.items():  # collect Pi atoms
            an = a.atomic_number
            ac = charge[n]
            ar = radical[n]
            ah = hyb(n) > 1
            if an == 6:
                if ah or ac or ar:  # unsaturated carbon
                    adj[n] = {}
                    alpha[n], beta[n], electrons[n] = basis[(6, ac, ar, ah)]
            elif an in (5, 7, 8, 16):
                try:
                    alpha[n], beta[n], electrons[n] = basis[(an, ac, ar, ah)]
                except KeyError:  # not parametrized or don't have Pi orbitals
                    continue
                adj[n] = {}
        if not adj:
            return 0.
        for n, m, _ in self.bonds():
            if n in adj and m in adj:
                adj[n][m] = adj[m][n] = min(beta[n], beta[m])

        energy = 0.
        for comp in self._connected_components(adj):
            h_matrix = zeros((len(comp), len(comp)))
            mapping = {}
            e = 0
            for i, n in enumerate(comp):
                h_matrix[i, i] = alpha[n]
                mapping[n] = i
                e += electrons[n]
            for n in comp:
                i = mapping[n]
                for m, b in adj[n].items():
                    j = mapping[m]
                    h_matrix[i, j] = h_matrix[j, i] = b

            orbs = sorted(list(eig(h_matrix)[0]))
            paired = e // 2
            energy += sum(x * 2 for x in orbs[:paired])
            if e % 2:  # unpaired
                energy += orbs[paired]
        return energy


__all__ = ['Huckel']
