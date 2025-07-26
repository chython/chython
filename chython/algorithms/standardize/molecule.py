# -*- coding: utf-8 -*-
#
#  Copyright 2018-2025 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from typing import List, TYPE_CHECKING, Union, Tuple
from ._charged import fixed_rules, morgan_rules
from ._groups import *
from ._tautomers import rules as tautomers_rules
from ._metal_organics import rules as metal_rules
from ...containers.bonds import Bond
from ...exceptions import ValenceError, ImplementationError
from ...periodictable import H as _H


if TYPE_CHECKING:
    from chython import MoleculeContainer


# atomic number constants
H = 1
C = 6


class Standardize:
    __slots__ = ()

    def canonicalize(self: 'MoleculeContainer', *, fix_tautomers=True, keep_kekule=False,
                     ignore_pyrrole_hydrogen=False, buffer_size=7,
                     logging=False, ignore=True) -> Union[bool, List[Tuple[Tuple[int, ...], int, str]]]:
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
        s = self.standardize(_fix_stereo=False, logging=True, ignore=ignore, fix_tautomers=fix_tautomers)
        h, changed = self.implicify_hydrogens(_fix_stereo=False, logging=True)

        if fix_tautomers and (logging or keep_kekule):  # thiele can change tautomeric form
            hgs = {n: a.implicit_hydrogens for n, a in self.atoms()}
        if keep_kekule:  # save bond orders
            bonds = [(b, b.order) for _, _, b in self.bonds()]

        t = self.thiele(fix_tautomers=fix_tautomers)
        if not t and fix_tautomers:  # after standardizations keto-enols stereo should be fixed
            self.fix_stereo()
        c = self.standardize_charges(prepare_molecule=False, logging=True)
        a = self.standardize_tautomers(prepare_molecule=False, logging=True)

        if keep_kekule and t:  # restore
            # check ring charge/hydrogen moving
            if c or fix_tautomers and any(hgs[n] != a.implicit_hydrogens for n, a in self.atoms()):
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
                if fix_tautomers and (x := tuple(n for n, a in self.atoms() if hgs[n] != a.implicit_hydrogens)):
                    s.append((x, -1, 'aromatic tautomer found'))
            if c:
                s.append((tuple(c), -1, 'recharged'))
            if a:
                s.append((tuple(a), -1, 'tautomers standardized'))
            if keep_kekule and t:
                if c or fix_tautomers and any(hgs[n] != a.implicit_hydrogens for n, a in self.atoms()):
                    s.append(((), -1, 'kekulized again'))
                else:
                    s.append(((), -1, 'kekule form restored'))
            return s
        return bool(k or s or h or t or c or a)

    def standardize(self: Union['MoleculeContainer', 'Standardize'], *, logging=False, ignore=True, fix_tautomers=True,
                    _fix_stereo=True) -> Union[bool, List[Tuple[Tuple[int, ...], int, str]]]:
        """
        Standardize functional groups. Return True if any non-canonical group found.

        :param fix_tautomers: convert tautomers to canonical forms.
        :param logging: return list of fixed atoms with matched rules.
        :param ignore: ignore standardization bugs.
        """
        r = self.fix_resonance(logging=True, _fix_stereo=False)
        if r:
            log = [(tuple(r), -1, 'resonance fixed')]
            fixed = set(r)
        else:
            log, fixed = [], set()

        l, f = self.__standardize(double_rules, fix_tautomers)
        log.extend(l)
        fixed.update(f)
        if f:
            l, f = self.__standardize(double_rules, fix_tautomers)  # double shot rules for overlapped groups
            log.extend(l)
            fixed.update(f)
        l, f = self.__standardize(single_rules, fix_tautomers)
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

        if fixed and _fix_stereo:
            self.fix_stereo()

        if logging:
            if fixed:
                log.append((tuple(fixed), -1, 'standardized atoms'))
            return log
        return bool(fixed)

    def standardize_charges(self: 'MoleculeContainer', *, logging=False, prepare_molecule=True,
                            _fix_stereo=True) -> Union[bool, List[int]]:
        """
        Set canonical positions of charges in heterocycles and ferrocenes.

        :param logging: return list of changed atoms.
        :param prepare_molecule: do thiele procedure.
        """
        changed: List[int] = []
        bonds = self._bonds
        nsc = self.not_special_connectivity
        atoms = self._atoms

        if prepare_molecule:
            self.thiele()

        seen = set()
        # not morgan
        for q, fix in fixed_rules:
            for mapping in q.get_mapping(self, automorphism_filter=False):
                match = set(mapping.values())
                if len(match.intersection(seen)) > 2:  # if matched more than 2 atoms
                    continue
                seen.update(match)
                # if not 2 neighbors and 1 hydrogen or 3 neighbors within 1st and second atoms - break
                atom_1, atom_2 = mapping[1], mapping[2]
                if len(bonds[atom_1]) == 2:
                    if not atoms[atom_1].implicit_hydrogens:
                        continue
                elif all(x == 4 for x in bonds[atom_1].values()):
                    continue

                if len(bonds[atom_2]) == 2:
                    if not atoms[atom_2].implicit_hydrogens:
                        continue
                elif all(x == 4 for x in bonds[atom_2].values()):
                    continue

                if fix:
                    atom_3 = mapping[3]
                    atoms[atom_3]._charge = 0
                    changed.append(atom_3)
                else:
                    atoms[atom_1]._charge = 0
                    changed.append(atom_1)
                atoms[atom_2]._charge = 1
                changed.append(atom_2)  # add atoms to changed

        # morgan
        pairs = []
        for q, fix in morgan_rules:
            for mapping in q.get_mapping(self, automorphism_filter=False):
                match = set(mapping.values())
                if len(match.intersection(seen)) > 2:  # if matched more than 2 atoms
                    continue
                seen.update(match)
                atom_1, atom_2 = mapping[1], mapping[2]
                if len(bonds[atom_1]) == 2:
                    if not atoms[atom_1].implicit_hydrogens:
                        continue
                elif all(x == 4 for x in bonds[atom_1].values()):
                    continue

                if len(bonds[atom_2]) == 2:
                    if not atoms[atom_2].implicit_hydrogens:
                        continue
                elif all(x == 4 for x in bonds[atom_2].values()):
                    continue

                if fix:
                    atom_3 = mapping[3]
                    atoms[atom_3]._charge = 0
                    changed.append(atom_3)
                else:
                    # remove charge from 1st N atom
                    atoms[atom_1]._charge = 0
                pairs.append((atom_1, atom_2, fix))

        # ferrocene
        fcr = []
        for r in self.sssr:
            if len(r) != 5 or not all(atoms[n].hybridization == 4 for n in r):
                continue
            ch = [(n, x) for n in r if (x := atoms[n].charge)]
            if len(ch) != 1 or ch[0][1] != -1:
                continue
            ch = ch[0][0]
            ca = [n for n in r if atoms[n] == C and
                  (len(bs := nsc[n]) == 2 or len(bs) == 3 and any(b == 1 for b in bonds[n].values()))]
            if len(ca) < 2 or ch not in ca:
                continue
            atoms[ch]._charge = 0  # reset charge for morgan recalculation
            fcr.append(ca)
            changed.append(ch)

        if pairs or fcr:
            self.__dict__.pop('atoms_order', None)  # remove cached morgan
            for atom_1, atom_2, fix in pairs:
                if self.atoms_order[atom_1] > self.atoms_order[atom_2]:
                    atoms[atom_2]._charge = 1
                    changed.append(atom_2)
                    if not fix:
                        changed.append(atom_1)
                else:
                    atoms[atom_1]._charge = 1
                    if fix:
                        changed.append(atom_1)

            for ca in fcr:
                n = min(ca, key=self.atoms_order.get)
                atoms[n]._charge = -1
                changed.append(n)
            del self.__dict__['atoms_order']  # remove invalid morgan

        if changed:
            self.flush_cache(keep_sssr=True, keep_components=True)  # clear cache
            if _fix_stereo:
                self.fix_stereo()
            if logging:
                return changed
            return True
        if logging:
            return []
        return False

    def standardize_tautomers(self: 'MoleculeContainer', *, logging=False, prepare_molecule=True,
                              _fix_stereo=True) -> Union[bool, List[int]]:
        """
        Set canonical positions of hydrogens in azoles, guanidines, etc.

        :param logging: return a list of changed atoms.
        :param prepare_molecule: apply thiele procedure.
        """
        changed: List[int] = []
        atoms = self._atoms
        bonds = self._bonds

        if prepare_molecule:
            self.thiele()

        seen = set()
        pairs = []
        nitrogens = set()
        doubles = set()
        protonated = set()
        carbons = set()
        for q, r_type in tautomers_rules:
            for mapping in q.get_mapping(self, automorphism_filter=False):
                match = set(mapping.values())
                if len(match.intersection(seen)) > 2:  # if matched more than 2 atoms
                    continue
                seen.update(match)
                atom_1, atom_2 = mapping[1], mapping[2]
                if r_type == 0:  # move hydrogen
                    atoms[atom_1]._implicit_hydrogens = 1
                    atoms[atom_2]._implicit_hydrogens = 0
                    changed.append(atom_1)
                    changed.append(atom_2)
                elif r_type == 1:
                    atoms[atom_1]._implicit_hydrogens = 0  # remove hydrogen for morgan calculation
                    pairs.append((atom_1, atom_2, False))
                elif r_type == 2:
                    atom_3 = mapping[3]
                    atoms[atom_3]._implicit_hydrogens = 0
                    pairs.append((atom_1, atom_2, True))
                    changed.append(atom_3)
                else:  # r_type == 3
                    nitrogens.add(atom_2)
                    if atom_1 not in nitrogens:
                        atom_3 = mapping[3]
                        carbons.add(atom_3)
                        nitrogens.add(atom_1)
                        # make symmetric form for morgan
                        bonds[atom_1][atom_3]._order = 1  # reset double bond
                        atoms[atom_1]._implicit_hydrogens += 1  # add hydrogen
                        doubles.add((atom_1, atom_3))
                        protonated.add(atom_1)

        if pairs or carbons:
            self.__dict__.pop('atoms_order', None)  # remove cached morgan
            for atom_1, atom_2, fix in pairs:
                if self.atoms_order[atom_1] >= self.atoms_order[atom_2]:
                    # keep as is
                    atoms[atom_1]._implicit_hydrogens = 1
                    if fix:  # type 2 is always changing
                        changed.append(atom_1)
                else:  # move hydrogen
                    atoms[atom_2]._implicit_hydrogens = 1
                    changed.append(atom_2)
                    if not fix:
                        changed.append(atom_1)
            if carbons:
                nitrogens = sorted(nitrogens, key=self.atoms_order.get)
                while nitrogens:
                    atom_1 = nitrogens.pop()
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
            del self.__dict__['atoms_order']  # remove invalid morgan

        if changed:
            self.flush_cache(keep_sssr=True, keep_components=True)  # clear cache
            if _fix_stereo:
                self.fix_stereo()
            if logging:
                return changed
            return True
        if logging:
            return []
        return False

    def remove_coordinate_bonds(self: 'MoleculeContainer', *, keep_to_terminal=True, _fix_stereo=True) -> int:
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
            self.flush_cache(keep_sssr=True)
            self.calc_labels()
            if _fix_stereo:
                self.fix_stereo()
        return len(ab)

    def implicify_hydrogens(self: 'MoleculeContainer', *, logging=False, _fix_stereo=True) -> \
            Union[int, Tuple[int, List[int]]]:
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
            if _fix_stereo:
                self.fix_stereo()

        if logging:
            return len(to_remove), list(fixed)
        return len(to_remove)

    def explicify_hydrogens(self: 'MoleculeContainer', *, start_map=None, _return_map=False, _fix_stereo=True) -> \
            Union[int, List[Tuple[int, int]]]:
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
            if _fix_stereo:
                self.fix_stereo()
            if _return_map:
                return log
            return len(to_add)
        elif _return_map:
            return []
        return 0

    def check_valence(self: 'MoleculeContainer') -> List[int]:
        """
        Check valences of all atoms.

        :return: list of invalid atoms
        """
        # only invalid atoms have None hydrogens.
        return [n for n, a in self.atoms() if a.implicit_hydrogens is None]

    def clean_isotopes(self: 'MoleculeContainer') -> bool:
        """
        Clean isotope marks from molecule.
        Return True if any isotope found.
        """
        isotopes = [a for _, a in self.atoms() if a.isotope]
        if isotopes:
            for i in isotopes:
                i._isotope = None
            self.flush_cache(keep_sssr=True, keep_components=True)
            self.fix_stereo()
            return True
        return False

    def __standardize(self: 'MoleculeContainer', rules, fix_tautomers):
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
            self.flush_cache(keep_sssr=keep_sssr, keep_components=keep_components)
            # recalculate isomorphism labels
            self.calc_labels()
            for n in hs:  # hydrogens count recalculation
                self.calc_implicit(n)
            fixed.update(hs)
        return log, fixed


__all__ = ['Standardize']
