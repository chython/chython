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

# Descriptor bucket boundaries (MW in Da)
_MW_BOUNDS = (150, 250, 350, 450, 550)


def fingerprint_schema() -> dict[str, int]:
    """
    Return the stable mapping from feature name to bit index used by
    ``functional_fingerprint()``.  Import this once in chemder to build
    the selector's ``mapping`` array.

    Bit layout
    ----------
    0 – 5   rings_count   buckets 0,1,2,3,4,5+
    6 – 11  hba_count     buckets 0,1,2,3,4,5+
    12 – 17 hbd_count     buckets 0,1,2,3,4,5+
    18 – 22 sp3_fraction  quintiles 0-20%, 20-40%, 40-60%, 60-80%, 80-100%
    23 – 28 molecular_mass buckets <150,150-250,250-350,350-450,450-550,550+
    29 – 34 rotatable_bonds buckets 0,1,2,3,4-6,7+
    35 +    functional group presence (one bit per pattern in functional_rules)
    """
    schema: dict[str, int] = {}
    for i, name in enumerate(('rings_0', 'rings_1', 'rings_2', 'rings_3', 'rings_4', 'rings_5p')):
        schema[name] = i
    for i, name in enumerate(('hba_0', 'hba_1', 'hba_2', 'hba_3', 'hba_4', 'hba_5p')):
        schema[name] = 6 + i
    for i, name in enumerate(('hbd_0', 'hbd_1', 'hbd_2', 'hbd_3', 'hbd_4', 'hbd_5p')):
        schema[name] = 12 + i
    for i, name in enumerate(('sp3_q0', 'sp3_q1', 'sp3_q2', 'sp3_q3', 'sp3_q4')):
        schema[name] = 18 + i
    for i, name in enumerate(('mw_lt150', 'mw_150', 'mw_250', 'mw_350', 'mw_450', 'mw_550p')):
        schema[name] = 23 + i
    for i, name in enumerate(('rbc_0', 'rbc_1', 'rbc_2', 'rbc_3', 'rbc_4_6', 'rbc_7p')):
        schema[name] = 29 + i
    for offset, fg_name in enumerate(functional_rules):
        schema[fg_name] = 35 + offset
    return schema


class FunctionalGroups:
    __slots__ = ()

    @cached_property
    def functional_fingerprint(self) -> bytes:
        """
        256-bit (32-byte) binary feature vector for fast AVX2 bitmap filtering.

        Bit layout matches ``fingerprint_schema()``:
          0– 5  rings_count   buckets  0 / 1 / 2 / 3 / 4 / 5+
          6–11  hba           buckets  0 / 1 / 2 / 3 / 4 / 5+
         12–17  hbd           buckets  0 / 1 / 2 / 3 / 4 / 5+
         18–22  sp3 fraction  quintiles 0–20% / 20–40% / 40–60% / 60–80% / 80–100%
         23–28  molecular_mass  <150 / 150–250 / 250–350 / 350–450 / 450–550 / 550+
         29–34  rotatable_bonds  0 / 1 / 2 / 3 / 4–6 / 7+
         35+    one bit per pattern in functional_rules (presence only)
        """
        bits = bytearray(32)

        def _set(idx):
            bits[idx >> 3] |= 1 << (idx & 7)

        # rings_count bucket (0–5+)
        rc = min(self.rings_count, 5)
        _set(rc)

        # hba bucket (0–5+)
        hba = min(self.hydrogen_bond_acceptors_count, 5)
        _set(6 + hba)

        # hbd bucket (0–5+)
        hbd = min(self.hydrogen_bond_donors_count, 5)
        _set(12 + hbd)

        # sp3 fraction quintile (0–4)
        sp3q = min(int(self.carbon_sp3_fraction * 5), 4)
        _set(18 + sp3q)

        # molecular_mass bucket
        mw = self.molecular_mass
        if mw < 150:
            _set(23)
        elif mw < 250:
            _set(24)
        elif mw < 350:
            _set(25)
        elif mw < 450:
            _set(26)
        elif mw < 550:
            _set(27)
        else:
            _set(28)

        # rotatable_bonds bucket: 0/1/2/3/4-6/7+
        rbc = self.rotatable_bonds_count
        if rbc == 0:
            _set(29)
        elif rbc == 1:
            _set(30)
        elif rbc == 2:
            _set(31)
        elif rbc == 3:
            _set(32)
        elif rbc <= 6:
            _set(33)
        else:
            _set(34)

        # functional group bits (one bit per named pattern, presence only)
        fgs = self.functional_groups
        for offset, name in enumerate(functional_rules):
            if name in fgs:
                _set(35 + offset)

        return bytes(bits)

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
        seen = set()
        for name, (q, keep, *_) in protective_rules.items():
            c = 0
            for mp in q.get_mapping(self, automorphism_filter=False):
                atoms = {m for n, m in mp.items() if n not in keep}
                if seen.isdisjoint(atoms):
                    seen.update(atoms)
                    c += 1
            if c:
                found[name] = c
        return found

    def remove_protection(self, name=None, canonicalize=True,
                          fix_tautomers=True, ignore_pyrrole_hydrogen=False) -> bool:
        """
        Remove protective groups from the given molecule if applicable.

        :param name: Specific protective group name to remove (None = all).
        :param canonicalize: Run full canonicalization after removal.
        :param fix_tautomers: Canonicalize tautomer forms. Passed to canonicalize().
        :param ignore_pyrrole_hydrogen: Fix invalid rings like Cn1cc[nH]c1. Passed to canonicalize().
        """
        to_delete = set()
        to_add = []
        if name is None:
            rules = protective_rules.values()
        elif name in protective_rules:
            rules = [protective_rules[name]]
        else:
            raise ValueError(f'Unknown protective group: {name}')

        kept_atoms = set()
        for q, keep, add, *_ in rules:
            for mp in q.get_mapping(self, automorphism_filter=False):
                delete = {m for n, m in mp.items() if n not in keep}
                if not to_delete.isdisjoint(delete):
                    continue
                to_delete.update(delete)
                for n in keep:
                    kept_atoms.add(mp[n])
                for n, a, b in add:
                    to_add.append((mp[n], a, b))

        for n, a, b in to_add:
            m = self.add_atom(a, _skip_calculation=True)
            self.add_bond(m, n, b, _skip_calculation=True)
        for n in to_delete:
            self.delete_atom(n, _skip_calculation=True)
        if to_delete or to_add:
            self.fix_structure()
            if canonicalize:
                self.canonicalize(fix_tautomers=fix_tautomers,
                                 ignore_pyrrole_hydrogen=ignore_pyrrole_hydrogen)
            else:
                # fix implicit H on aromatic N freed from PG
                for n in kept_atoms:
                    a = self.atom(n)
                    if a.atomic_symbol == 'N' and a.hybridization == 4 and a.implicit_hydrogens is None:
                        a._implicit_hydrogens = 1
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
        for name, fg_name, _, reactor in oxidation_rules:
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
        for name, fg_name, _, reactor in reduction_rules:
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
        for name, fg_name, _, reactor in transformation_rules:
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


__all__ = ['FunctionalGroups', 'fingerprint_schema']
