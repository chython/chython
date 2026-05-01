# -*- coding: utf-8 -*-
#
#  Copyright 2026 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from ._protective import rules as protective_rules


class FunctionalGroups:
    __slots__ = ()

    @cached_property
    def protective_groups(self) -> list[str]:
        """
        List of protective group found in the molecule.
        """

        return sorted({pg for pg, (q, *_) in protective_rules.items() if q < self})

    def remove_protection(self, name=None) -> bool:
        """
        Remove protective groups from the given molecule if applicable.
        """
        seen = set()
        to_delete = set()
        to_add = []
        if name is None:
            rules = protective_rules.values()
        elif name in protective_rules:
            rules = [protective_rules[name]]
        else:
            raise ValueError(f'Unknown protective group: {name}')

        for q, keep, add, *_ in rules:
            for mp in q.get_mapping(self, automorphism_filter=False):
                if not seen.isdisjoint(mp.values()):
                    continue
                seen.update(mp.values())
                for n, m in mp.items():
                    if n not in keep:
                        to_delete.add(m)
                for n, a, b in add:
                    to_add.append((mp[n], a, b))

        for n, a, b in to_add:
            m = self.add_atom(a, _skip_calculation=True)
            self.add_bond(m, n, b, _skip_calculation=True)
        for n in to_delete:
            self.delete_atom(n, _skip_calculation=True)
        if to_delete or to_add:
            self._changed.difference_update(to_delete)
            self.fix_structure()
            self.fix_stereo()
            return True
        return False


__all__ = ['FunctionalGroups']
