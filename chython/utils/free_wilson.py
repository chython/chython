# -*- coding: utf-8 -*-
#
#  Copyright 2022-2024 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from itertools import product
from typing import Union, List
from .retro import Tree
from ..containers import MoleculeContainer, QueryContainer
from ..periodictable import H


def fw_prepare_groups(core: Union[MoleculeContainer, QueryContainer], molecule: MoleculeContainer) -> \
       List[MoleculeContainer]:
    """
    Prepare list of core with connected groups. Hydrogens added to groups for marking connection point.
    Hydrogens have isotope marks equal to mapping of core atoms.
    Groups connected multiple times (rings) - contains multiple hydrogens.

    :param core: core structure for searching
    :param molecule: target structure
    """
    try:
        core_map = next(core.get_mapping(molecule))
    except StopIteration:
        return []

    reverse = {v: k for k, v in core_map.items()}
    cs = set(core_map.values())
    groups = molecule.substructure(molecule._atoms.keys() - cs, recalculate_hydrogens=False)
    gs = set(groups)

    cf = molecule.substructure(cs, recalculate_hydrogens=False)

    for n, m, b in molecule.bonds():
        if n in cs:
            if m in gs:
                a = molecule.atom(n)
                h = H(x=a.x, y=a.y)
                h._isotope = reverse[n]  # mark mapping to isotope
                groups.add_bond(groups.add_atom(h, _skip_calculation=True), m, b.copy(), _skip_calculation=True)

                a = molecule.atom(m)
                h = H(x=a.x, y=a.y)
                h._isotope = reverse[n]  # mark mapping to isotope
                cf.add_bond(cf.add_atom(h, _skip_calculation=True), n, b.copy(), _skip_calculation=True)
        elif m in cs and n in gs:
            a = molecule.atom(m)
            h = H(x=a.x, y=a.y)
            h._isotope = reverse[m]
            groups.add_bond(groups.add_atom(h, _skip_calculation=True), n, b.copy(), _skip_calculation=True)

            a = molecule.atom(n)
            h = H(x=a.x, y=a.y)
            h._isotope = reverse[m]  # mark mapping to isotope
            cf.add_bond(cf.add_atom(h.copy(), _skip_calculation=True), n, b.copy(), _skip_calculation=True)
    groups = groups.split()
    groups.insert(0, cf)
    return groups


def fw_decomposition_tree(groups: List[MoleculeContainer]) -> Tree:
    assert len(groups) == len(set(groups))

    pred = {}  # directed graph from substructures to superstructures
    succ = {}
    for m in groups:
        pred[m] = set()
        succ[m] = set()

    for m in groups:
        for n in groups:
            if m < n:
                pred[n].add(m)
                succ[m].add(n)

    # break triangles
    scope = {m for m, ns in succ.items() if len(ns) > 1}
    while scope:
        m = sorted(scope, key=lambda x: len(pred[x]))[0]
        s = succ[m]
        scope.discard(m)
        while True:
            for x, y in product((x for x in s if succ[x]), (x for x in s if len(pred[x]) > 1)):
                if y in succ[x]:
                    s.discard(y)
                    pred[y].discard(m)
                    break
            else:
                break

    def _rec_tree(x):
        return x, [_rec_tree(y) for y in succ[x]]

    m = MoleculeContainer()
    m.add_atom('H')

    return m, [_rec_tree(x) for x, p in pred.items() if not p]


__all__ = ['fw_prepare_groups', 'fw_decomposition_tree']
