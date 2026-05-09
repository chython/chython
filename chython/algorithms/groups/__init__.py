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
from collections.abc import Iterator
from functools import cached_property
from itertools import permutations
from ._functional import rules as functional_rules
from ._oxidations import rules as oxidation_rules
from ._protective import rules as protective_rules
from ._reactions import rules as reaction_rules
from ._reductions import rules as reduction_rules
from ._transformations import rules as transformation_rules


class FunctionalGroups:
    __slots__ = ()

    @cached_property
    def functional_groups(self) -> dict[str, int]:
        """
        Dict of functional group names to their count in the molecule.
        """
        found = {}
        for name, q in functional_rules.items():
            c = sum(1 for _ in q.get_mapping(self))
            if c:
                found[name] = c
        return found

    @cached_property
    def protective_groups(self) -> dict[str, int]:
        """
        Dict of protective group names to their count in the molecule.
        """
        found = {}
        for name, (q, *_) in protective_rules.items():
            c = sum(1 for _ in q.get_mapping(self))
            if c:
                found[name] = c
        return found

    def remove_protection(self, name=None) -> bool:
        """
        Remove protective groups from the given molecule if applicable.
        """
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
                delete = {m for n, m in mp.items() if n not in keep}
                if not to_delete.isdisjoint(delete):
                    continue
                to_delete.update(delete)
                for n, a, b in add:
                    to_add.append((mp[n], a, b))

        for n, a, b in to_add:
            m = self.add_atom(a, _skip_calculation=True)
            self.add_bond(m, n, b, _skip_calculation=True)
        for n in to_delete:
            self.delete_atom(n, _skip_calculation=True)
        if to_delete or to_add:
            self.fix_structure()
            self.fix_stereo()
            return True
        return False

    def react(self, *others, reaction=None) -> Iterator[tuple[str, 'ReactionContainer']]:
        """
        Enumerate possible reaction products between molecules.

        mol1.react(mol2) -> [(reaction_name, ReactionContainer), ...]
        mol1.react(mol2, mol3) -> [(reaction_name, ReactionContainer), ...]  # multi-component
        mol1.react(mol2, reaction='suzuki') -> only suzuki coupling

        :param reaction: optional reaction name to apply selectively.
        """
        mols = [self, *others]

        for name, fg_names, reactor in reaction_rules:
            if reaction is not None and name != reaction:
                continue
            if len(fg_names) != len(mols):
                continue
            for perm in permutations(mols):
                if all(fg in mol.functional_groups for mol, fg in zip(perm, fg_names)):
                    for rxn in reactor(*perm):
                        yield name, rxn
                    break

    def oxidize(self, reaction=None) -> Iterator[tuple[str, 'ReactionContainer']]:
        """
        Enumerate possible single-step oxidation products.

        mol.oxidize() -> [(reaction_name, ReactionContainer), ...]
        mol.oxidize(reaction='alcohol_to_aldehyde') -> only this oxidation

        :param reaction: optional reaction name to apply selectively.
        """
        fgs = self.functional_groups
        for name, fg_name, reactor in oxidation_rules:
            if reaction is not None and name != reaction:
                continue
            if fg_name in fgs:
                for rxn in reactor(self):
                    yield name, rxn

    def reduce(self, reaction=None) -> Iterator[tuple[str, 'ReactionContainer']]:
        """
        Enumerate possible single-step reduction products.

        mol.reduce() -> [(reaction_name, ReactionContainer), ...]
        mol.reduce(reaction='ketone_to_alcohol') -> only this reduction

        :param reaction: optional reaction name to apply selectively.
        """
        fgs = self.functional_groups
        for name, fg_name, reactor in reduction_rules:
            if reaction is not None and name != reaction:
                continue
            if fg_name in fgs:
                for rxn in reactor(self):
                    yield name, rxn

    def transform(self, reaction=None) -> Iterator[tuple[str, 'ReactionContainer']]:
        """
        Enumerate possible single-molecule functional group interconversions
        (ring formations from open-chain precursors with implicit reagents).

        mol.transform() -> [(reaction_name, ReactionContainer), ...]
        mol.transform(reaction='appel') -> only Appel reaction

        :param reaction: optional reaction name to apply selectively.
        """
        fgs = self.functional_groups
        for name, fg_name, reactor in transformation_rules:
            if reaction is not None and name != reaction:
                continue
            if fg_name in fgs:
                for rxn in reactor(self):
                    yield name, rxn

    def __invert__(self):
        """
        Enumerate all possible single-step molecular transformations
        (oxidations, reductions, and functional group interconversions).

        ~mol -> [(reaction_name, ReactionContainer), ...]
        """
        yield from self.oxidize()
        yield from self.reduce()
        yield from self.transform()

    def __matmul__(self, other):
        """
        Enumerate possible reaction products between molecules.

        mol1 @ mol2 -> [(reaction_name, ReactionContainer), ...]
        mol1 @ [mol2, mol3] -> [(reaction_name, ReactionContainer), ...]  # multi-component
        """
        if isinstance(other, (list, tuple)):
            return self.react(*other)
        return self.react(other)


__all__ = ['FunctionalGroups']
