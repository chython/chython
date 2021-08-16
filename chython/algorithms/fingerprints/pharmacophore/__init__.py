# -*- coding: utf-8 -*-
#
#  Copyright 2021 Aleksandr Sizov <murkyrussian@gmail.com>
#  Copyright 2021 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from typing import Dict, TYPE_CHECKING
from ._acceptor import queries as acceptor_queries, banned as acceptor_banned
from ._acidic import queries as acidic_queries
from ._basic import queries as basic_queries, banned as basic_banned
from ._donor import queries as donor_queries

if TYPE_CHECKING:
    from ....containers import MoleculeContainer

queries = {
    'acc_allow': acceptor_queries,
    'acc_ban': acceptor_banned,
    'acid_allow': acidic_queries,
    'base_allow': basic_queries,
    'base_ban': basic_banned,
    'donor_allow': donor_queries,
}


class Pharmacophore:
    __slots__ = ()
    """
    Matcher of atoms in an molecule to some integer number.
    Each integer represents 6 bit, each bit is enable/disable some pharmacophoric feature:
    H acceptor | acidic prop. | aromatic | basic prop. | H donor | some halogen atom
    After the constructing, binary number will be converted to 10-based integer 
    """

    @cached_property
    def pharmacophores(self: MoleculeContainer) -> Dict[int, int]:
        """
        Match id of each atom of some molecule with integer number which define some features using in FCFP
        """
        id_sets = {}
        for name, qrs in queries.items():
            first, others = qrs[0], qrs[1:]
            first = {dct[1] for dct in first.get_mapping(self)}
            first = first.union(*({dct[1] for dct in q.get_mapping(self)} for q in others))
            id_sets[name] = first

        # rules for looking for aromatic and halogen atoms used here with direct getting ids atoms from molecule's info
        bins = [
            id_sets['acc_allow'] - id_sets['acc_ban'],
            id_sets['acid_allow'],
            {idx for idx in self._atoms if self.hybridization(idx) == 4},
            id_sets['base_allow'] - id_sets['base_ban'],
            id_sets['donor_allow'],
            {idx for idx, atom in self.atoms() if atom.atomic_number in {9, 17, 35, 53} and self.neighbors(idx) < 2},
        ]

        out = {idx: 0 for idx in self._atoms}
        for pos, bin_ in zip((1, 2, 4, 8, 16, 32), bins):
            for idx in bin_:
                out[idx] |= pos

        return out


__all__ = ['Pharmacophore']
