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
from typing import Union, List
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
    hs = molecule._hydrogens
    hgs = groups._hydrogens
    plane = molecule._plane

    for n, m, b in molecule.bonds():
        if n in cs:
            if m in gs:
                h = H()
                h._Core__isotope = reverse[n]  # mark mapping to isotope
                groups.add_bond(groups.add_atom(h, xy=plane[n]), m, b.copy())
                hgs[m] = hs[m]  # restore H count
        elif m in cs and n in gs:
            h = H()
            h._Core__isotope = reverse[m]
            groups.add_bond(groups.add_atom(h, xy=plane[m]), n, b.copy())
            hgs[n] = hs[n]
    groups = groups.split()
    groups.insert(0, molecule.substructure(cs, recalculate_hydrogens=False))
    return groups


__all__ = ['fw_prepare_groups']
