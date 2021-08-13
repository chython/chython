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
from typing import Iterator, Tuple, TYPE_CHECKING
from .acceptor import queries as acceptor_queries, banned as acceptor_banned
from .acidic import queries as acidic_queries
from .aromatic import queries as aromatic_queries
from .basic import queries as basic_queries, banned as basic_banned
from .donor import queries as donor_queries
from .halogen import queries as halogen_queries

if TYPE_CHECKING:
    from ....containers import MoleculeContainer, QueryContainer

queries = {
    'acc_allow': acceptor_queries,
    'acc_ban': acceptor_banned,
    'acid_allow': acidic_queries,
    'ar_allow': aromatic_queries,
    'base_allow': basic_queries,
    'base_ban': basic_banned,
    'donor_allow': donor_queries,
    'hal_allow': halogen_queries,
}


class Pharmacophore:
    __slots__ = ()
    """
    Matcher of atoms in an molecule to some integer number.
    Each integer represents 6 bit, each bit is enable/disable some pharmacophoric feature:
    H acceptor | acidic prop. | aromatic | basic prop. | H donor | some halogen atom
    After the constructing, binary number will be converted to 10-based integer 
    """

    def features(self: MoleculeContainer) -> Iterator[Tuple[int, int]]:
        """
        Match id of each atom of some molecule with integer number which define some features using in FCFP
        """
        id_sets = {}
        for name, qrs in queries.items():
            if len(qrs) == 1:
                id_sets[name] = {dct[1] for dct in qrs[0].get_mapping(self)}
                continue
            first, others = qrs[0], qrs[1:]
            first = {dct[1] for dct in first.get_mapping(self)}
            first = first.union(*({dct[1] for dct in q.get_mapping(self)} for q in others))
            id_sets[name] = first

        bins = [
            id_sets['acc_allow'] - id_sets['acc_ban'],
            id_sets['acid_allow'],
            id_sets['ar_allow'],
            id_sets['base_allow'] - id_sets['base_ban'],
            id_sets['donor_allow'],
            id_sets['hal_allow'],
        ]

        for idx, _ in self.atoms():
            yield idx, int(''.join('1' if idx in set_ else '0' for set_ in bins), base=2)


__all__ = ['Pharmacophore']
