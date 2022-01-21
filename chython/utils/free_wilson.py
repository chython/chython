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
from collections import defaultdict
from itertools import product
from typing import Union, List, Dict
from ..containers import MoleculeContainer, QueryContainer


def fw_prepare_groups(core: Union[MoleculeContainer, QueryContainer], molecule: MoleculeContainer):
    """
    Prepare dict of core atoms in keys and list of connected groups in values.
    Hydrogens added to groups for marking connection point.
    Groups connected multiple times (rings) - contains multiple hydrogens and presented in each key.

    :param core: core structure for searching
    :param molecule: target structure
    """
    try:
        core_map = next(core.get_mapping(molecule))
    except StopIteration:
        return {}

    reverse = {v: k for k, v in core_map.items()}
    cs = set(core_map.values())
    groups = molecule - cs
    gs = set(groups)

    mapping = {}
    for n, m, b in molecule.bonds():
        if n in cs:
            if m in gs:
                x = groups.add_atom('H')
                mapping[x] = n
                groups.add_bond(x, m, b.copy())
        elif m in cs and n in gs:
            x = groups.add_atom('H')
            mapping[x] = m
            groups.add_bond(x, n, b.copy())

    r_groups = defaultdict(list)
    for g in groups.split():
        for n in set(g).intersection(mapping):
            r_groups[reverse[mapping[n]]].append(g)
    return dict(r_groups)


def fw_onehot_groups(groups: List[Dict[int, List[MoleculeContainer]]]):
    """
    Prepare one-hot matrix with found groups separated by connection points.

    :param groups: List of groups returned from `fw_prepare_groups`
    """
    unique = set()
    for kgs in groups:
        unique.update((k, g) for k, gs in kgs.items() for g in gs)

    mapping = {x: {i} for i, x in enumerate(sorted(unique, key=lambda x: x[0]))}
    q_mapping = [(k, g.substructure(set(g), as_query=True, skip_neighbors_marks=True, skip_hybridizations_marks=False,
                                    skip_hydrogens_marks=True, skip_rings_sizes_marks=True),
                  v)
                 for (k, g), v in mapping.items()]
    for (a, xa, ia), (b, xb, ib) in product(q_mapping, repeat=2):
        if a == b and xa < xb:
            ib.update(ia)

    out = [[k for k, _ in mapping], [g for _, g in mapping]]  # legend
    for kgs in groups:
        onehot = [0] * len(mapping)
        for k, gs in kgs.items():
            for g in gs:
                for i in mapping[(k, g)]:
                    onehot[i] = 1
        out.append(onehot)
    return out


__all__ = ['fw_prepare_groups', 'fw_onehot_groups']
