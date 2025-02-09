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
from collections import ChainMap
from itertools import count
from typing import List, Tuple, TYPE_CHECKING, Union
from ._reactions import rules


if TYPE_CHECKING:
    from chython import ReactionContainer


class FixMapper:
    __slots__ = ()

    def fix_mapping(self: 'ReactionContainer', *, logging: bool = False) -> \
            Union[bool, List[Tuple[int, str, Tuple[int, ...]]]]:
        """
        Fix mapping by using loaded rules.
        """
        if not self:
            if logging:
                return []
            return False

        cgr = ~self
        if not cgr.center_atoms:
            if logging:
                return []
            return False
        del self.__dict__['__cached_method_compose']

        log = []
        free_number = count(max(cgr) + 1)
        components = [(cgr.substructure(c),
                       cgr.augmented_substructure(c, 2),  # deep DEPENDS on rules!
                       c)
                      for c in cgr.substructure(cgr.center_atoms).connected_components]

        r_atoms = ChainMap(*(x._atoms for x in self.reactants))
        for c, ac, cs in components:
            for rule_num, (query, signature, restrict, fix) in enumerate(rules):
                if str(c) == signature:
                    for mapping in query.get_mapping(ac, automorphism_filter=False):
                        if not cs.issubset(mapping.values()):
                            continue
                        if restrict is not None and any(a != r_atoms.get(mapping[n]) for n, a in restrict.atoms()):
                            continue
                        mapping = {mapping[n]: next(free_number) if m is None else mapping[m] for n, m in fix.items()}
                        for m in self.products:
                            m.remap(mapping)
                        log.append((rule_num, signature, tuple(mapping.values())))
                        break
                    else:
                        continue
                    break  # component remapped!

        if log:
            self.flush_cache()
            if logging:
                return log
            return True
        elif logging:
            return log
        return False


__all__ = ['FixMapper']
