# -*- coding: utf-8 -*-
#
#  Copyright 2022-2024 Ramil Nugmanov <nougmanoff@protonmail.com>
#  Copyright 2024 Philippe Gantzer <p.gantzer@icredd.hokudai.ac.jp>
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
from CachedMethods import class_cached_property
from itertools import chain, count, repeat
from logging import getLogger, INFO
from numpy import ix_, unravel_index, argmax, zeros, array, isclose, nonzero, ones, mean
from typing import TYPE_CHECKING, Union


if TYPE_CHECKING:
    from chython import ReactionContainer

logger = getLogger('chython.attention')
logger.setLevel(INFO)


class Attention:
    __slots__ = ()
    __class_cache__ = {}

    def reset_mapping(self: Union['ReactionContainer', 'Attention'], *, return_score: bool = False, multiplier=1.75,
                      keep_reactants_numbering=False) -> Union[bool, float]:
        """
        Do atom-to-atom mapping. Return True if mapping changed.
        """
        if any(len(bs) > 14 for m in self.molecules() for bs in m._bonds.values()):
            logger.info('atom-to-atom mapping not supported for hypervalent compounds')
            return False
        fixed = self.__fix_collisions()
        equal_atoms, p2r, r2p, r_adj, p_adj, r_map, p_map, pa, rg_map = self.__prepare_remapping()

        # rxnmapper-inspired algorithm
        am = self.__get_attention()
        # sum of reactants to products attention and vice-versa for equal atom types only
        am = (am[p2r] + am[r2p].T) * equal_atoms
        amc = am.copy()

        mapping = {}
        scope = zeros(pa, dtype=bool)
        seen = ones(pa, dtype=bool)
        score = []
        for x in range(pa):  # iteratively map each product atom to reactant
            # select highest attention
            # todo: optimize
            if not x:
                i, j = unravel_index(argmax(am), am.shape)
            else:
                ams = am[scope]
                if ams.size:
                    i, j = unravel_index(argmax(ams), ams.shape)
                    i = nonzero(scope)[0][i]
                else:
                    i, j = unravel_index(argmax(am), am.shape)
            if isclose(am[i, j], 0.):  # no more products atoms in reactants
                # mark as unmapped
                for n in set(p_map).difference(mapping):
                    mapping[n] = 0
                break
            else:
                score.append(amc[i, j])
            mapping[p_map[i]] = r_map[j]
            am[ix_(p_adj[i], r_adj[j])] *= multiplier  # highlight neighbors
            am[i] = am[:, j] = 0  # mask mapped product and reactant atoms
            seen[i] = False
            scope[i] = False
            scope[p_adj[i] & seen] = True

        score = float(mean(score)) if score else 0.
        # mapping done.
        if any(n != m for n, m in mapping.items()):  # old mapping changed
            if keep_reactants_numbering or fixed:
                r_mapping = {n: n for n in r_map}
            else:
                r_mapping = {m: n for n, m in enumerate(r_map, 1)}  # remap reactants to contiguous range
                for m in self.reactants:
                    m.remap(r_mapping)

            p_mapping = {}
            nm = max(r_mapping.values())
            for n, m in mapping.items():
                if m := r_mapping.get(m):
                    p_mapping[n] = m
                else:  # not found in reactants atoms. set unique numbers.
                    nm += 1
                    p_mapping[n] = nm

            for m in self.products:
                m.remap(p_mapping)

            if not keep_reactants_numbering and not fixed:
                rg_mapping = {m: n for n, m in enumerate(rg_map, nm+1)}  # remap reagents to contiguous range without overlapping
                for m in self.reagents:
                    m.remap(rg_mapping)

            self.flush_cache()
            fixed = True

        if self.fix_groups_mapping():  # fix carboxy etc
            fixed = True
        if self.fix_mapping():  # fix common mistakes in mechanisms
            fixed = True
        if return_score:
            return score
        return fixed

    def __fix_collisions(self: 'ReactionContainer'):
        r = [n for m in chain(self.reactants, self.reagents) for n in m._atoms]
        p = [n for m in self.products for n in m._atoms]
        c = count(1)
        if len(r) != len(set(r)):
            for m in chain(self.reactants, self.reagents):
                m.remap({n: next(c) for n in m._atoms})
        if len(p) != len(set(p)):
            for m in self.products:
                m.remap({n: next(c) for n in m._atoms})
        if next(c) != 1:
            self.flush_cache()
            return True
        return False

    def __prepare_remapping(self: 'ReactionContainer'):
        r_map = [n for m in self.reactants for n in m]
        p_map = [n for m in self.products for n in m]
        rg_map = [n for m in self.reagents for n in m]
        ra = len(r_map)  # number of reactants atoms
        pa = len(p_map)  # number of products atoms

        ram = [False]  # reactants atoms mask
        r_atoms = []
        r_adj = zeros((ra, ra), dtype=bool)
        i = 0
        for m in self.reactants:
            ram.append(False)
            ram.extend(repeat(True, len(m)))
            a = m.adjacency_matrix()
            j = i + len(m)
            r_adj[i:j, i:j] = a
            i = j
            r_atoms.extend(a.atomic_number for _, a in m.atoms())
        r_atoms = array(r_atoms, dtype=int)

        pam = [False] * len(ram)  # products atoms mask
        p_atoms = []
        p_adj = zeros((pa, pa), dtype=bool)
        i = 0
        for m in self.products:
            pam.append(False)
            pam.extend(repeat(True, len(m)))
            a = m.adjacency_matrix()
            j = i + len(m)
            p_adj[i:j, i:j] = a
            i = j
            p_atoms.extend(a.atomic_number for _, a in m.atoms())
        p_atoms = array(p_atoms, dtype=int)

        ram.extend(repeat(False, len(pam) - len(ram)))
        ram = array(ram, dtype=bool)
        pam = array(pam, dtype=bool)
        return p_atoms[:, None] == r_atoms, ix_(pam, ram), ix_(ram, pam), r_adj, p_adj, r_map, p_map, pa, rg_map

    @class_cached_property
    def __attention_model(self):
        from chython import torch_device
        from chytorch.zoo.rxnmap import Model

        return Model().to(torch_device)

    @class_cached_property
    def __autocast(self):
        from chython import torch_device

        if torch_device.startswith('cuda'):
            try:
                from torch import autocast
            except ImportError:  # torch 1.8 ad-hoc
                from torch.cuda.amp import autocast

                return autocast()
            else:
                return autocast('cuda')
        return autocast_filler()

    def __get_attention(self):
        from torch import no_grad

        with no_grad(), self.__autocast:
            am = self.__attention_model(self).float().cpu().numpy()
        return am


class autocast_filler:
    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        ...


__all__ = ['Attention']
