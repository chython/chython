# -*- coding: utf-8 -*-
#
#  Copyright 2018-2026 Ramil Nugmanov <nougmanoff@protonmail.com>
#  Copyright 2021 Dmitrij Zanadvornykh <zandmitrij@gmail.com>
#  Copyright 2018 Tagir Akhmetshin <tagirshin@gmail.com>
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
from ._groups import rules as fragment_rules
from ._isomers import rules as isomer_rules
from ._metal_organics import rules as metal_rules
from ...containers.bonds import Bond
from ...exceptions import ValenceError, ImplementationError
from ...periodictable import H as _H


H = 1


class Standardize:
    __slots__ = ()

    def canonicalize(self, *, fix_tautomers=True, keep_kekule=False,
                     ignore_pyrrole_hydrogen=False, buffer_size=7,
                     logging=False, ignore=True) -> bool | list:
        """
        Convert molecule to canonical forms of functional groups and aromatic rings without explicit hydrogens.

        :param logging: return log.
        :param ignore: ignore standardization bugs.
        :param fix_tautomers: convert tautomers to canonical forms.
        :param keep_kekule: return kekule form.
        :param buffer_size: number of attempts of pyridine form searching.
        :param ignore_pyrrole_hydrogen: ignore hydrogen on pyrrole to fix invalid rings like Cn1cc[nH]c1.
        """
        k = self.kekule(buffer_size=buffer_size, ignore_pyrrole_hydrogen=ignore_pyrrole_hydrogen)
        s = self.standardize(logging=True, ignore=ignore, fix_tautomers=fix_tautomers)
        h, changed = self.implicify_hydrogens(logging=True)

        if logging or keep_kekule:  # thiele can change tautomeric form
            hgs = {n: a.implicit_hydrogens for n, a in self.atoms()}
        if keep_kekule:  # save bond orders
            bonds = [(b, b.order) for _, _, b in self.bonds()]

        t = self.thiele(fix_tautomers=fix_tautomers)
        a = self.standardize_isomers(logging=True)

        if keep_kekule and t:  # restore
            # check ring charge/hydrogen moving
            if a or any(hgs[n] != aa.implicit_hydrogens for n, aa in self.atoms()):
                self.kekule()  # we need to do full kekule again
            else:
                for b, o in bonds:  # noqa
                    b._order = o
                self.flush_cache()
                self.calc_labels()

        if logging:
            if k:
                s.insert(0, ((), -1, 'kekulized'))
            if h:
                s.append((tuple(changed), -1, 'implicified'))
            if t:
                s.append(((), -1, 'aromatized'))
                if fix_tautomers and (x := tuple(n for n, aa in self.atoms() if hgs[n] != aa.implicit_hydrogens)):
                    s.append((x, -1, 'aromatic tautomer found'))
            if a:
                s.append((tuple(a), -1, 'isomers standardized'))
            if keep_kekule and t:
                if a or any(hgs[n] != aa.implicit_hydrogens for n, aa in self.atoms()):
                    s.append(((), -1, 'kekulized again'))
                else:
                    s.append(((), -1, 'kekule form restored'))
            return s
        return bool(k or s or h or t or a)

    def standardize(self, *, logging=False, ignore=True, fix_tautomers=True) -> bool | list:
        """
        Standardize functional groups. Return True if any non-canonical group found.

        :param fix_tautomers: convert tautomers to canonical forms.
        :param logging: return list of fixed atoms with matched rules.
        :param ignore: ignore standardization bugs.
        """
        log, fixed = [], set()

        l, f = self.__standardize(fragment_rules, fix_tautomers)
        log.extend(l)
        fixed.update(f)
        l, f = self.__standardize(metal_rules, fix_tautomers)  # metal-organics fix
        log.extend(l)
        fixed.update(f)

        if b := fixed.intersection(n for n, a in self.atoms() if a.implicit_hydrogens is None):
            if ignore:
                log.append((tuple(b), -1, 'standardization failed'))
            else:
                raise ImplementationError(f'standardization leads to invalid valences: {b}')

        if fixed:
            self.fix_stereo()

        if logging:
            if fixed:
                log.append((tuple(fixed), -1, 'standardized atoms'))
            return log
        return bool(fixed)

    def standardize_isomers(self, *, logging=False) -> bool | list:
        """
        Set canonical positions of charges and hydrogens.

        :param logging: return list of changed atoms.
        """
        changed: list[int] = []
        bonds = self._bonds
        atoms = self._atoms

        seen = set()
        charge_morgan_pairs = []
        fcr = []
        tautomer_morgan_pairs = []
        nitrogens = set()
        doubles = set()
        protonated = set()
        carbons = set()

        for q, rule_type in isomer_rules:
            for mapping in q.get_mapping(self, automorphism_filter=False):
                match = set(mapping.values())
                if len(match.intersection(seen)) > 2:
                    continue
                seen.update(match)

                if rule_type == 'charge_fixed':
                    atom_1, atom_2 = mapping[1], mapping[2]
                    atoms[atom_1]._charge = 0
                    atoms[atom_2]._charge = 1
                    changed.append(atom_1)
                    changed.append(atom_2)

                elif rule_type == 'charge_fixed_bridge':
                    atom_1, atom_2, atom_3 = mapping[1], mapping[2], mapping[3]
                    atoms[atom_3]._charge = 0
                    atoms[atom_2]._charge = 1
                    changed.append(atom_3)
                    changed.append(atom_2)

                elif rule_type == 'tautomer_fixed':
                    atom_1, atom_2 = mapping[1], mapping[2]
                    atoms[atom_1]._implicit_hydrogens = 1
                    atoms[atom_2]._implicit_hydrogens = 0
                    changed.append(atom_1)
                    changed.append(atom_2)

                elif rule_type == 'charge_morgan':
                    atom_1, atom_2 = mapping[1], mapping[2]
                    atoms[atom_1]._charge = 0
                    charge_morgan_pairs.append((atom_1, atom_2, atom_1))

                elif rule_type == 'charge_morgan_bridge':
                    atom_1, atom_2, atom_3 = mapping[1], mapping[2], mapping[3]
                    atoms[atom_3]._charge = 0
                    charge_morgan_pairs.append((atom_1, atom_2, atom_3))

                elif rule_type == 'ferrocene':
                    atoms[mapping[1]]._charge = 0
                    fcr.append([mapping[i] for i in range(1, 6)])
                    changed.append(mapping[1])

                elif rule_type == 'tautomer_morgan':
                    atom_1, atom_2 = mapping[1], mapping[2]
                    atoms[atom_1]._implicit_hydrogens = 0
                    tautomer_morgan_pairs.append((atom_1, atom_2, atom_1))

                elif rule_type == 'tautomer_morgan_donor':
                    atom_1, atom_2, atom_3 = mapping[1], mapping[2], mapping[3]
                    atoms[atom_3]._implicit_hydrogens = 0
                    tautomer_morgan_pairs.append((atom_1, atom_2, atom_3))

                elif rule_type == 'amidine':
                    atom_1, atom_2 = mapping[1], mapping[2]
                    nitrogens.add(atom_2)
                    if atom_1 not in nitrogens:
                        atom_3 = mapping[3]
                        carbons.add(atom_3)
                        nitrogens.add(atom_1)
                        bonds[atom_1][atom_3]._order = 1
                        atoms[atom_1]._implicit_hydrogens += 1
                        doubles.add((atom_1, atom_3))
                        protonated.add(atom_1)

        if charge_morgan_pairs or fcr or tautomer_morgan_pairs or carbons:
            self.__dict__.pop('atoms_order', None)
            self.__dict__.pop('int_adjacency', None)

            for atom_1, atom_2, source in charge_morgan_pairs:
                if self.atoms_order[atom_1] > self.atoms_order[atom_2]:
                    atoms[atom_2]._charge = 1
                    if source != atom_2:
                        changed.append(source)
                        changed.append(atom_2)
                else:
                    atoms[atom_1]._charge = 1
                    if source != atom_1:
                        changed.append(source)
                        changed.append(atom_1)

            for ca in fcr:
                n = min(ca, key=self.atoms_order.get)
                atoms[n]._charge = -1
                changed.append(n)

            for atom_1, atom_2, source in tautomer_morgan_pairs:
                if self.atoms_order[atom_1] >= self.atoms_order[atom_2]:
                    atoms[atom_1]._implicit_hydrogens = 1
                    if source != atom_1:
                        changed.append(source)
                        changed.append(atom_1)
                else:
                    atoms[atom_2]._implicit_hydrogens = 1
                    if source != atom_2:
                        changed.append(source)
                        changed.append(atom_2)

            if carbons:
                nitrogens_sorted = sorted(nitrogens, key=self.atoms_order.get)
                while nitrogens_sorted:
                    atom_1 = nitrogens_sorted.pop()
                    for atom_3 in sorted(bonds[atom_1], key=self.atoms_order.get):
                        if atom_3 in carbons:
                            carbons.discard(atom_3)
                            atoms[atom_1]._implicit_hydrogens -= 1
                            bonds[atom_1][atom_3]._order = 2
                            if (atom_1, atom_3) not in doubles:
                                changed.append(atom_1)
                                changed.append(atom_3)
                            break
                    else:
                        if atom_1 in protonated:
                            changed.append(atom_1)

            del self.__dict__['atoms_order']
            del self.__dict__['int_adjacency']

        if changed:
            self.flush_cache(keep_sssr=True, keep_components=True, keep_special_connectivity=True)
            self.fix_stereo()
            if logging:
                return changed
            return True
        if logging:
            return []
        return False

    def remove_coordinate_bonds(self, *, keep_to_terminal=True) -> int:
        """Remove coordinate (or hydrogen) bonds marked with 8 (any) bond

        :param keep_to_terminal: Keep any bonds to terminal hydrogens
        :return: removed bonds count
        """
        ab = [(n, m) for n, m, b in self.bonds() if b == 8]

        if keep_to_terminal:
            skeleton = self.not_special_connectivity
            hs = {n for n, a in self.atoms() if a == H and not skeleton[n]}
            ab = [(n, m) for n, m in ab if n not in hs and m not in hs]

        for n, m in ab:
            self.delete_bond(n, m, _skip_calculation=True)

        if ab:
            self.flush_cache(keep_sssr=True, keep_special_connectivity=True)
            self.calc_labels()
            self.fix_stereo()
        return len(ab)

    def implicify_hydrogens(self, *, logging=False) -> int | tuple:
        """
        Remove explicit hydrogen if possible. Return number of removed hydrogens.
        Works only with Kekule forms of aromatic structures.
        Keeps isotopes of hydrogen.

        :param logging: return list of changed atoms.
        """
        atoms = self._atoms
        bonds = self._bonds

        explicit = defaultdict(list)
        for n, atom in atoms.items():
            if atom == H and (atom.isotope is None or atom.isotope == 1):
                if len(bonds[n]) > 1:
                    raise ValenceError(f'Hydrogen atom {n} has invalid valence. Try to use remove_coordinate_bonds()')
                for m, b in bonds[n].items():
                    if b == 1:
                        if atoms[m] != H:  # not H-H
                            explicit[m].append(n)
                    elif b != 8:
                        raise ValenceError(f'Hydrogen atom {n} has invalid valence {b.order}.')

        to_remove = set()
        fixed = {}
        for n, hs in explicit.items():
            atom = atoms[n]
            len_h = len(hs)
            for i in range(len_h, 0, -1):
                hi = hs[:i]
                explicit_sum = 0
                explicit_dict = defaultdict(int)
                for m, bond in bonds[n].items():
                    if m not in hi and bond != 8:
                        explicit_sum += bond.order
                        explicit_dict[(bond.order, atoms[m].atomic_number)] += 1
                try:
                    # aromatic rings don't match any rule
                    rules = atom.valence_rules(explicit_sum)
                except ValenceError:
                    break
                for s, d, h in rules:
                    if s.issubset(explicit_dict) and all(explicit_dict[k] >= c for k, c in d.items()) and h >= i:
                        to_remove.update(hi)
                        fixed[n] = h
                        break
                else:
                    continue
                break

        for n in to_remove:
            del atoms[n]
            for m in bonds.pop(n):
                del bonds[m][n]

        for n, h in fixed.items():
            atoms[n]._implicit_hydrogens = h

        if to_remove:
            self.flush_cache(keep_sssr=True)
            self.calc_labels()
            self.fix_stereo()

        if logging:
            return len(to_remove), list(fixed)
        return len(to_remove)

    def explicify_hydrogens(self, *, start_map=None, _return_map=False) -> int | list:
        """
        Add explicit hydrogens to atoms.

        :return: number of added atoms
        """
        atoms = self._atoms
        to_add = []
        for n, a in atoms.items():
            try:
                to_add.extend([n] * a.implicit_hydrogens)
            except TypeError:
                raise ValenceError(f'atom {n} has valence error')

        if to_add:
            log = []
            bonds = self._bonds
            m = start_map if start_map is not None else max(atoms) + 1
            for n in to_add:
                atoms[m] = _H(implicit_hydrogens=0)
                bonds[n][m] = b = Bond(1)
                bonds[m] = {n: b}
                atoms[n]._implicit_hydrogens = 0
                log.append((n, m))
                m += 1

            self.flush_cache(keep_sssr=True)
            self.calc_labels()
            self.fix_stereo()
            if _return_map:
                return log
            return len(to_add)
        elif _return_map:
            return []
        return 0

    def check_valence(self) -> list[int]:
        """
        Check valences of all atoms.

        :return: list of invalid atoms
        """
        # only invalid atoms have None hydrogens.
        return [n for n, a in self.atoms() if a.implicit_hydrogens is None]

    def clean_isotopes(self) -> bool:
        """
        Clean isotope marks from molecule.
        Return True if any isotope found.
        """
        isotopes = [a for _, a in self.atoms() if a.isotope]
        if isotopes:
            for i in isotopes:
                i._isotope = None
            self.flush_cache(keep_sssr=True, keep_components=True, keep_special_connectivity=True)
            self.fix_stereo()
            return True
        return False

    def __standardize(self, rules, fix_tautomers):
        atoms = self._atoms
        bonds = self._bonds

        log = []
        fixed = set()
        for r, (pattern, atom_fix, bonds_fix, any_atoms, is_tautomer) in enumerate(rules):
            if not fix_tautomers and is_tautomer:
                continue
            keep_sssr = keep_components = True
            hs = set()
            seen = set()
            for mapping in pattern.get_mapping(self, automorphism_filter=False):
                match = set(mapping.values())
                if not match.isdisjoint(seen):  # skip intersected groups
                    continue
                if any_atoms:  # accept overlapping of Any-atoms
                    seen.update(match - {mapping[n] for n in any_atoms})
                else:
                    seen.update(match)
                for n, (ch, ir) in atom_fix.items():
                    n = mapping[n]
                    hs.add(n)
                    a = atoms[n]
                    a._charge += ch
                    if a.charge > 4:
                        a._charge -= ch
                        log.append((tuple(match), r, f'bad charge formed. changes omitted: {pattern}'))
                        break  # skip changes
                    if ir is not None:
                        a._is_radical = ir
                else:
                    for n, m, bo in bonds_fix:
                        n = mapping[n]
                        m = mapping[m]
                        hs.add(n)
                        hs.add(m)
                        if m in bonds[n]:
                            b = bonds[n][m]
                            if b == 8 or bo == 8:
                                keep_sssr = False
                            b._order = bo
                        else:  # new bond
                            keep_sssr = keep_components = False
                            bonds[n][m] = bonds[m][n] = Bond(bo)
                    log.append((tuple(match), r, str(pattern)))

            if not hs:  # not matched
                continue
            self.flush_cache(keep_sssr=keep_sssr, keep_components=keep_components, keep_special_connectivity=keep_sssr)
            # recalculate isomorphism labels
            self.calc_labels()
            for n in hs:  # hydrogens count recalculation
                self.calc_implicit(n)
            fixed.update(hs)
        return log, fixed


__all__ = ['Standardize']
