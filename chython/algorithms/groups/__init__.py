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
from typing import NamedTuple, Optional
from ._functional import rules as functional_rules
from ._oxidations import rules as oxidation_rules
from ._protective import rules as protective_rules
from ._reactions import rules as reaction_rules
from ._reductions import rules as reduction_rules
from ._roles import roles as role_rules, transformers as role_transformers
from ._transformations import rules as transformation_rules


class EnumeratedReaction(NamedTuple):
    name: str                    # reaction rule name (e.g. 'suzuki', 'appel')
    reaction: 'ReactionContainer'  # the enumerated reaction


class StickyFragment(NamedTuple):
    role: str              # role name of the single attachment point
    sticky_left: str       # open bond at the start ('-...'), glues onto the left
    sticky_right: str      # open bond at the end ('...-'), glues onto the right
    canonical_smiles: str  # canonical SMILES of the [At]-capped fragment (dedup key)


class StickyLinker(NamedTuple):
    role_left: str         # role at the [210At] end
    role_right: str        # role at the [211At] end
    sticky_left: str       # '-A...B-' role_left end open first (glues that end left)
    sticky_right: str      # '-B...A-' role_right end open first (linker flipped)
    canonical_smiles: str  # canonical SMILES of the capped linker (dedup key)


def _present_transformers(role, fgs):
    """
    Resolve a role selector to the (role_name, transformer) pairs whose handle
    functional group is actually present in ``fgs``.

    role=None -> every role; a role name -> just that role. Raises ValueError on
    an unknown role name. Filtering against ``fgs`` (the molecule's cached
    functional_groups) up front means we never run a transformer for a handle the
    molecule does not carry -- capping a cut never creates a new coupling handle,
    so a handle absent from the source is absent from every intermediate too.
    """
    if role is None:
        items = role_transformers.items()
    elif role in role_rules:
        items = [(role, role_transformers[role])]
    else:
        raise ValueError(f'Unknown role: {role}')
    return [(role_name, transformer) for role_name, entries in items
            for fg_name, transformer in entries if fg_name in fgs]


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

    def remove_protection(self, names=None, canonicalize=True,
                          fix_tautomers=True, ignore_pyrrole_hydrogen=False,
                          *, start=None, logging=False) -> 'bool | list':
        """
        Remove protective groups from the given molecule if applicable.

        :param names: protective group names to remove, as a list. ``None``
            removes all; a list removes exactly that set (useful to strip
            incidental protecting groups while leaving an intrinsic one in
            place). Unknown names raise ``ValueError``.
        :param canonicalize: Run full canonicalization after removal.
        :param fix_tautomers: Canonicalize tautomer forms. Passed to canonicalize().
        :param ignore_pyrrole_hydrogen: Fix invalid rings like Cn1cc[nH]c1. Passed to canonicalize().
        :param start: Starting atom number for newly added atoms (e.g. restored carbonyl O).
            When None, defaults to max(existing) + 1.
        :param logging: return the list of freed (kept) atom numbers instead of a
            bool. These are the heteroatoms retained from each protecting-group
            match (e.g. the revealed amine nitrogen) — the atoms that become newly
            reactive, ready to feed a sticky_fragments/sticky_linkers ``masked`` set.
        """
        to_delete = set()
        to_add = []
        if names is None:
            rules = protective_rules.values()
        else:
            unknown = set(names) - protective_rules.keys()
            if unknown:
                raise ValueError(f'Unknown protective group(s): {sorted(unknown)}')
            rules = [protective_rules[n] for n in names]

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

        m = start
        for n, a, b in to_add:
            m = self.add_atom(a, n=m, _skip_calculation=True)
            self.add_bond(m, n, b, _skip_calculation=True)
            m += 1
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
            return list(kept_atoms) if logging else True
        return [] if logging else False

    def react(self, *others, reaction=None) -> Iterator['EnumeratedReaction']:
        """
        Enumerate possible reaction products between molecules.

        mol1.react(mol2) -> [EnumeratedReaction(name, ReactionContainer), ...]
        mol1.react(mol2, mol3) -> [EnumeratedReaction(...), ...]  # multi-component
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
                        yield EnumeratedReaction(name, rxn)
                    break

    def oxidize(self, reaction=None) -> Iterator['EnumeratedReaction']:
        """
        Enumerate possible single-step oxidation products.

        mol.oxidize() -> [EnumeratedReaction(name, ReactionContainer), ...]
        mol.oxidize(reaction='alcohol_to_aldehyde') -> only this oxidation

        :param reaction: optional reaction name to apply selectively.
        """
        fgs = self.functional_groups
        for name, fg_name, _, reactor in oxidation_rules:
            if reaction is not None and name != reaction:
                continue
            if fg_name in fgs:
                for rxn in reactor(self):
                    yield EnumeratedReaction(name, rxn)

    def reduce(self, reaction=None) -> Iterator['EnumeratedReaction']:
        """
        Enumerate possible single-step reduction products.

        mol.reduce() -> [EnumeratedReaction(name, ReactionContainer), ...]
        mol.reduce(reaction='ketone_to_alcohol') -> only this reduction

        :param reaction: optional reaction name to apply selectively.
        """
        fgs = self.functional_groups
        for name, fg_name, _, reactor in reduction_rules:
            if reaction is not None and name != reaction:
                continue
            if fg_name in fgs:
                for rxn in reactor(self):
                    yield EnumeratedReaction(name, rxn)

    def transform(self, reaction=None) -> Iterator['EnumeratedReaction']:
        """
        Enumerate possible single-molecule functional group interconversions
        (ring formations from open-chain precursors with implicit reagents).

        mol.transform() -> [EnumeratedReaction(name, ReactionContainer), ...]
        mol.transform(reaction='appel') -> only Appel reaction

        :param reaction: optional reaction name to apply selectively.
        """
        fgs = self.functional_groups
        for name, fg_name, _, reactor in transformation_rules:
            if reaction is not None and name != reaction:
                continue
            if fg_name in fgs:
                for rxn in reactor(self):
                    yield EnumeratedReaction(name, rxn)

    def sticky_fragments(self, role: Optional[str] = None, *,
                         masked=None, tries: int = 10,
                         hydrogens: bool = False) -> Iterator['StickyFragment']:
        """
        Enumerate mono-attachment sticky fragments.

        mol.sticky_fragments() -> yields StickyFragment for all known roles.
        mol.sticky_fragments('aryl_halide') -> only that role.

        Each result carries the attachment point in three forms: ``sticky_left``
        (open bond first, ``-...`` — glues onto the left of another fragment),
        ``sticky_right`` (open bond last, ``...-`` — glues onto the right) and
        ``canonical_smiles`` (the [At]-capped fragment, the dedup key).

        One result is yielded per match; dedup by canonical_smiles downstream.

        :param role: optional role name to restrict enumeration.
        :param masked: optional atom numbers barred from becoming the attachment
            point. A fragment whose fresh [At] cap hangs off a masked atom is
            skipped. Use it after ``remove_protection(logging=True)`` so a group
            revealed by an incidental deprotection ("cleaved FG") never becomes a
            standalone coupling handle. Numbers refer to this molecule's atoms and
            survive into the transformer products unchanged.
        :param tries: attempts per sticky-SMILES generation (forwarded to
            ``sticky_smiles``). The right-terminal traversal is a bounded random
            search; raise this for topologies where 10 tries occasionally fail.
        :param hydrogens: show implicit hydrogens in the generated sticky SMILES.
        """
        # no cap rule disconnects a connected molecule; only a salt/mixture input
        # would leave a detached counter-ion, which sticky_smiles cannot serialize.
        if self.connected_components_count != 1:
            return
        masked = frozenset(masked or ())
        for role_name, transformer in _present_transformers(role, self.functional_groups):
            for product in transformer(self):
                # the transformer guarantees exactly one freshly created [At]
                # cap, and a new atom is always the highest atom number.
                n_at = max(product)
                # a masked atom (a group revealed by incidental deprotection) must
                # never participate in the coupling this fragment represents. It
                # participates in two ways, both barred:
                #   - as the attachment site: the fresh [At] cap hangs off it
                #     (addition/substitution roles, e.g. the freed amine as an
                #     alkyl_amine handle);
                #   - as a leaving group: the transform consumes (deletes) it while
                #     capping a neighbour (deaminative/deoxy/decarboxy roles, e.g.
                #     the freed amine as an alkyl_deamino leaving group).
                # Atom numbers survive into products unchanged, so a masked atom is
                # "consumed" iff it is absent from the product's atom set.
                if next(iter(product._bonds[n_at])) in masked or not masked <= set(product):
                    continue
                # the bond that joins two pieces belongs to the LEFT partner: the
                # left form keeps its leading bond, the right form drops its
                # trailing bond, so ``A.sticky_right + B.sticky_left`` emits exactly
                # one bond token (two would be "2 bonds in a row" -> invalid SMILES).
                left = product.sticky_smiles(left=n_at, remove_left=True, keep_bond_left=True,
                                             tries=tries, hydrogens=hydrogens)
                right = product.sticky_smiles(right=n_at, remove_right=True,
                                              tries=tries, hydrogens=hydrogens)
                yield StickyFragment(role_name, left, right, str(product))

    def sticky_linkers(self, role_left: Optional[str] = None,
                       role_right: Optional[str] = None, *,
                       masked=None, tries: int = 10,
                       hydrogens: bool = False) -> Iterator['StickyLinker']:
        """
        Enumerate bi-attachment sticky linkers. The first (left) role is capped
        with [210At], the second (right) role with [211At].

        mol.sticky_linkers() -> all role pairs.
        mol.sticky_linkers('aryl_halide', 'aryl_acyl') -> restrict both ends.

        Each result carries the whole linker in two open-bond traversals ready to
        concatenate: ``sticky_left`` (role_left end open first, ``-A...B-``) and
        ``sticky_right`` (role_right end open first, the flipped ``-B...A-``), plus
        ``canonical_smiles`` (the [210At]/[211At]-capped structure, the dedup key,
        always 210=left / 211=right).

        One result per (left match, right match); dedup downstream.

        :param role_left: optional role for the [210At] end (None = all roles).
        :param role_right: optional role for the [211At] end (None = all roles).
        :param masked: optional atom numbers barred from the LEFT (step-1, [210At])
            end only; the RIGHT (step-2, [211At]) end is exempt. A masked atom is
            one revealed by an intrinsic deprotection whose sole coupling role is
            the deferred second step (e.g. a Boc-masked amine): it may sit on the
            right but never the left. This also removes the (amine_left, X_right) /
            (X_left, amine_right) ordering duplication for such reserved ends.
        :param tries: attempts per sticky-SMILES generation (forwarded to
            ``sticky_smiles``). The terminal-atom traversal is a bounded random
            search; raise this for topologies where 10 tries occasionally fail.
        :param hydrogens: show implicit hydrogens in the generated sticky SMILES.
        """
        # as in sticky_fragments: no cap rule disconnects a connected molecule,
        # so only a salt/mixture input could yield a detached component.
        if self.connected_components_count != 1:
            return
        masked = frozenset(masked or ())
        fgs = self.functional_groups
        # capping a cut never creates a coupling handle, so every right handle
        # must already be present in the source. Resolve both ends up front and
        # bail before any transform work if either end has no present handle --
        # this avoids running left transforms only to find no right FG exists.
        left_transformers = _present_transformers(role_left, fgs)
        right_transformers = _present_transformers(role_right, fgs)
        if not left_transformers or not right_transformers:
            return

        for left_name, left_t in left_transformers:
            for inter in left_t(self):
                # the transformer guarantees one fresh [At] cap, always the
                # highest atom number; label it 210 before the second stage.
                n210 = max(inter)
                # source atom the 210 cap hangs off (caps are terminal -> one neighbor)
                core210 = next(iter(inter._bonds[n210]))
                # a masked atom must not drive the step-1 (left) coupling, whether
                # as the attachment site (cap hangs off it) or as a leaving group
                # consumed by the left transform (absent from the intermediate).
                # The step-2 (right) end stays exempt from `masked` by design.
                if core210 in masked or not masked <= set(inter):
                    continue
                inter.atom(n210).isotope = 210

                for right_name, right_t in right_transformers:
                    for product in right_t(inter):
                        n211 = max(product)
                        # both cuts landing on the same atom would collapse the
                        # linker to a single core (e.g. bromoacetic acid halide +
                        # decarboxy -> [At]C[At]); require >=1 atom between caps.
                        # the right (step-2) end is exempt from `masked` by design.
                        if next(iter(product._bonds[n211])) == core210:
                            continue
                        product.atom(n211).isotope = 211
                        # two traversal orientations of the same linker: 210-first
                        # (-A...B-) and 211-first (-B...A-), both open-bond both-ends.
                        # same bond-ownership rule as sticky_fragments: each open
                        # end keeps its LEADING bond and drops its TRAILING one, so
                        # a frag-linker-frag chain (``A.sticky_right + linker +
                        # B.sticky_left``) emits exactly one bond per junction.
                        sticky_left = product.sticky_smiles(left=n210, right=n211,
                                                            remove_left=True, keep_bond_left=True,
                                                            remove_right=True, keep_bond_right=False,
                                                            tries=tries, hydrogens=hydrogens)
                        sticky_right = product.sticky_smiles(left=n211, right=n210,
                                                             remove_left=True, keep_bond_left=True,
                                                             remove_right=True, keep_bond_right=False,
                                                             tries=tries, hydrogens=hydrogens)
                        yield StickyLinker(left_name, right_name,
                                           sticky_left, sticky_right, str(product))

    def __invert__(self):
        """
        Enumerate all possible single-step molecular transformations
        (oxidations, reductions, and functional group interconversions).

        ~mol -> [EnumeratedReaction(name, ReactionContainer), ...]
        """
        yield from self.oxidize()
        yield from self.reduce()
        yield from self.transform()

    def __matmul__(self, other):
        """
        Enumerate possible reaction products between molecules.

        mol1 @ mol2 -> [EnumeratedReaction(name, ReactionContainer), ...]
        mol1 @ [mol2, mol3] -> [EnumeratedReaction(...), ...]  # multi-component
        """
        if isinstance(other, (list, tuple)):
            return self.react(*other)
        return self.react(other)


__all__ = ['FunctionalGroups', 'fingerprint_schema', 'EnumeratedReaction',
           'StickyFragment', 'StickyLinker']
