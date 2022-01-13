# -*- coding: utf-8 -*-
#
#  Copyright 2022 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from itertools import repeat
from numpy import ix_, ones_like, errstate, unravel_index, nanargmax, zeros, array
from pkg_resources import resource_stream
from typing import TYPE_CHECKING, Union


if TYPE_CHECKING:
    from chython import ReactionContainer


class Attention:
    __slots__ = ()

    def reset_mapping(self: Union['ReactionContainer', 'Attention'], *, multiplier=90.,
                      keep_reactants_numbering=False) -> bool:
        """
        Do atom-to-atom mapping. Return True if mapping changed.
        """
        r_map, p_map, r_atoms, p_atoms, r_adj, p_adj, ram, pam, ra, pa = self.__prepare_remapping()
        am = self.__get_attention()

        am = am[ix_(pam, ram)] + am[ix_(ram, pam)].T  # sum of reactants to products attention and vice-versa
        am *= p_atoms[:, None] == r_atoms  # equal atom type attention

        mm = ones_like(am)
        mapping = {}
        for _ in range(pa):  # iteratively map each product atom to reactant
            # todo: optimize
            nam = am * mm
            with errstate(invalid='ignore'):
                nam /= nam[:, None, :].sum(2)  # normalized attention

            # select highest attention
            try:
                i, j = unravel_index(nanargmax(nam), nam.shape)
            except ValueError:  # no more products atoms in reactants
                # mark as unmapped
                for n in set(p_map).difference(mapping):
                    mapping[n] = 0
                break

            mapping[p_map[i]] = r_map[j]
            # update mm
            mm[ix_(p_adj[i], r_adj[j])] *= multiplier
            mm[i] = 0  # mask already mapped product atom
            mm[:, j] = 0  # mask already mapped reactant. todo: work with side products

        if any(n != m for n, m in mapping.items()):  # old mapping changed
            if keep_reactants_numbering:
                r_mapping = {n: n for n in r_map}
            else:
                r_mapping = {m: n for n, m in enumerate(r_map, 1)}  # remap reactants to contiguous range
                for m in self.reactants:
                    m.remap(r_mapping)

            p_mapping = {}
            nm = ra
            for n, m in mapping.items():
                if m := r_mapping.get(m):
                    p_mapping[n] = m
                else:  # not found in reactants atoms. set unique numbers.
                    nm += 1
                    p_mapping[n] = nm

            for m in self.products:
                m.remap(p_mapping)
            self.flush_cache()
            return True
        return False

    @class_cached_property
    def __attention_model(self):
        from chython import torch_device
        from chytorch.nn import ReactionEncoder
        from torch import load

        model = ReactionEncoder(d_model=1024, dim_feedforward=3072)
        model.load_state_dict(load(resource_stream(__package__, 'mapping.pt'), map_location=torch_device))
        model.eval()
        return model

    def __prepare_remapping(self: 'ReactionContainer'):
        r_map = [n for m in self.reactants for n in m]
        p_map = [n for m in self.products for n in m]
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
        return r_map, p_map, r_atoms, p_atoms, r_adj, p_adj, ram, pam, ra, pa

    def __get_attention(self):
        from chython import torch_device
        from chytorch.utils.data import ReactionDataset
        from torch import no_grad

        converter = ReactionDataset((self,), (None,), property_type=bool)
        a, n, i, e, r, _ = converter.collate((converter[0],))
        a.to(torch_device)
        n.to(torch_device)
        i.to(torch_device)
        e.to(torch_device)
        r.to(torch_device)

        with no_grad():
            am = self.__attention_model(a, n, i, e, r, True)[1].squeeze(0).cpu().numpy()  # raw attention matrix
        return am


__all__ = ['Attention']
