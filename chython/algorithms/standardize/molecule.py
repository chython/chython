# -*- coding: utf-8 -*-
#
#  Copyright 2018-2023 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from ._metal_organics import rules as metal_rules
from ...containers.bonds import Bond
from ...exceptions import ValenceError, ImplementationError
from ...periodictable import H


if TYPE_CHECKING:
    from chython import MoleculeContainer


class Standardize:
    __slots__ = ()

    def canonicalize(self: 'MoleculeContainer', *, fix_tautomers=True, keep_kekule=False,
                     logging=False, ignore=True) -> Union[bool, List[Tuple[Tuple[int, ...], int, str]]]:
        """
        Convert molecule to canonical forms of functional groups and aromatic rings without explicit hydrogens.

        :param logging: return log.
        :param ignore: ignore standardization bugs.
        :param fix_tautomers: convert tautomers to canonical forms.
        :param keep_kekule: return kekule form.
        """
        k = self.kekule()
        s = self.standardize(_fix_stereo=False, logging=True, ignore=ignore, fix_tautomers=fix_tautomers)
        h, changed = self.implicify_hydrogens(_fix_stereo=False, logging=True)

        if fix_tautomers and (logging or keep_kekule):  # thiele can change tautomeric form
            hgs = self._hydrogens.copy()
        if keep_kekule:  # save bond orders
            bonds = [(b, b.order) for _, _, b in self.bonds()]

        t = self.thiele(fix_tautomers=fix_tautomers)
        c = self.standardize_charges(prepare_molecule=False, logging=True)

        if keep_kekule and t:  # restore
            # check ring charge/hydrogen moving
            if c or fix_tautomers and hgs != self._hydrogens:  # noqa
                self.kekule()  # we need to do full kekule again
            else:
                for b, o in bonds:  # noqa
                    b._Bond__order = o  # noqa
            self.flush_cache()

        if logging:
            if k:
                s.insert(0, ((), -1, 'kekulized'))
            if h:
                s.append((tuple(changed), -1, 'implicified'))
            if t:
                s.append(((), -1, 'aromatized'))
                if fix_tautomers and hgs != self._hydrogens:
                    s.append((tuple(x for x, y in self._hydrogens.items() if hgs[x] != y),
                              -1, 'aromatic tautomer found'))
            if c:
                s.append((tuple(c), -1, 'recharged'))
            if keep_kekule and t:
                if c or fix_tautomers and hgs != self._hydrogens:
                    s.append(((), -1, 'kekulized again'))
                else:
                    s.append(((), -1, 'kekule form restored'))
            return s
        return bool(k or s or h or t or c)

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

        if b := fixed.intersection(n for n, h in self._hydrogens.items() if h is None):
            if ignore:
                log.append((tuple(b), -1, 'standardization failed'))
            else:
                raise ImplementationError(f'standardization leads to invalid valences: {b}')

        if fixed:
            self.flush_cache()
            if _fix_stereo:
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
        hydrogens = self._hydrogens
        charges = self._charges
        atoms = self._atoms
        hybridization = self.hybridization

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
                # if not 2 neighbors and 1 hydrogen  or 3 neighbors within 1st and second atoms - break
                atom_1, atom_2 = mapping[1], mapping[2]
                if len(bonds[atom_1]) == 2:
                    if not hydrogens[atom_1]:
                        continue
                elif all(x == 4 for x in bonds[atom_1].values()):
                    continue

                if len(bonds[atom_2]) == 2:
                    if not hydrogens[atom_2]:
                        continue
                elif all(x == 4 for x in bonds[atom_2].values()):
                    continue

                if fix:
                    atom_3 = mapping[3]
                    charges[atom_3] = 0
                    changed.append(atom_3)
                else:
                    charges[atom_1] = 0
                    changed.append(atom_1)
                charges[atom_2] = 1
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
                    if not hydrogens[atom_1]:
                        continue
                elif all(x == 4 for x in bonds[atom_1].values()):
                    continue

                if len(bonds[atom_2]) == 2:
                    if not hydrogens[atom_2]:
                        continue
                elif all(x == 4 for x in bonds[atom_2].values()):
                    continue

                if fix:
                    atom_3 = mapping[3]
                    charges[atom_3] = 0
                    changed.append(atom_3)
                else:
                    # remove charge from 1st N atom
                    charges[atom_1] = 0
                pairs.append((atom_1, atom_2, fix))

        if pairs:
            self.__dict__.pop('atoms_order', None)  # remove cached morgan
            for atom_1, atom_2, fix in pairs:
                if self.atoms_order[atom_1] > self.atoms_order[atom_2]:
                    charges[atom_2] = 1
                    changed.append(atom_2)
                    if not fix:
                        changed.append(atom_1)
                else:
                    charges[atom_1] = 1
                    if fix:
                        changed.append(atom_1)
            del self.__dict__['atoms_order']  # remove invalid morgan

        # ferrocene
        fcr = []
        for r in self.sssr:
            if len(r) != 5 or not all(hybridization(n) == 4 for n in r):
                continue
            ch = [(n, x) for n in r if (x := charges[n])]
            if len(ch) != 1 or ch[0][1] != -1:
                continue
            ch = ch[0][0]
            ca = [n for n in r if atoms[n].atomic_number == 6 and
                  (len(bs := nsc[n]) == 2 or len(bs) == 3 and any(b.order == 1 for b in bonds[n].values()))]
            if len(ca) < 2 or ch not in ca:
                continue
            charges[ch] = 0  # reset charge for morgan recalculation
            fcr.append(ca)
            changed.append(ch)
        if fcr:
            self.__dict__.pop('atoms_order', None)  # remove cached morgan
            for ca in fcr:
                n = min(ca, key=self.atoms_order.get)
                charges[n] = -1
                changed.append(n)
            del self.__dict__['atoms_order']  # remove invalid morgan

        if changed:
            self.flush_cache()  # clear cache
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
        bonds = self._bonds

        ab = [(n, m) for n, m, b in self.bonds() if b.order == 8]

        if keep_to_terminal:
            skeleton = self.not_special_connectivity
            hs = {n for n, a in self._atoms.items() if a.atomic_number == 1 and not skeleton[n]}
            ab = [(n, m) for n, m in ab if n not in hs and m not in hs]

        for n, m in ab:
            del bonds[n][m], bonds[m][n]

        if ab:
            self.flush_cache()
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
        charges = self._charges
        radicals = self._radicals
        bonds = self._bonds
        plane = self._plane
        hydrogens = self._hydrogens
        parsed_mapping = self._parsed_mapping

        explicit = defaultdict(list)
        for n, atom in atoms.items():
            if atom.atomic_number == 1 and (atom.isotope is None or atom.isotope == 1):
                if len(bonds[n]) > 1:
                    raise ValenceError(f'Hydrogen atom {n} has invalid valence. Try to use remove_coordinate_bonds()')
                for m, b in bonds[n].items():
                    if b.order == 1:
                        if atoms[m].atomic_number != 1:  # not H-H
                            explicit[m].append(n)
                    elif b.order != 8:
                        raise ValenceError(f'Hydrogen atom {n} has invalid valence {b.order}.')

        to_remove = set()
        fixed = {}
        for n, hs in explicit.items():
            atom = atoms[n]
            charge = charges[n]
            is_radical = radicals[n]
            len_h = len(hs)
            for i in range(len_h, 0, -1):
                hi = hs[:i]
                explicit_sum = 0
                explicit_dict = defaultdict(int)
                for m, bond in bonds[n].items():
                    if m not in hi and bond.order != 8:
                        explicit_sum += bond.order
                        explicit_dict[(bond.order, atoms[m].atomic_number)] += 1
                try:
                    # aromatic rings don't match any rule
                    rules = atom.valence_rules(charge, is_radical, explicit_sum)
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
            del charges[n]
            del radicals[n]
            del plane[n]
            del hydrogens[n]
            for m in bonds.pop(n):
                del bonds[m][n]
            try:
                del parsed_mapping[n]
            except KeyError:
                pass

        for n, h in fixed.items():
            hydrogens[n] = h

        if to_remove:
            self.flush_cache()
            self._conformers = [{x: y for x, y in c.items() if x not in to_remove} for c in self._conformers]  # noqa
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
        hydrogens = self._hydrogens
        to_add = []
        for n, h in hydrogens.items():
            try:
                to_add.extend([n] * h)
            except TypeError:
                raise ValenceError(f'atom {n} has valence error')

        if to_add:
            log = []
            bonds = self._bonds
            m = start_map
            for n in to_add:
                m = self.add_atom(H(), m)
                bonds[n][m] = bonds[m][n] = b = Bond(1)
                b._attach_graph(self, n, m)
                hydrogens[n] = 0
                log.append((n, m))
                m += 1

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
        return [n for n, h in self._hydrogens.items() if h is None]  # only invalid atoms have None hydrogens.

    def clean_isotopes(self: 'MoleculeContainer') -> bool:
        """
        Clean isotope marks from molecule.
        Return True if any isotope found.
        """
        atoms = self._atoms
        isotopes = [x for x in atoms.values() if x.isotope]
        if isotopes:
            for i in isotopes:
                i._Core__isotope = None
            self.flush_cache()
            self.fix_stereo()
            return True
        return False

    def __standardize(self: 'MoleculeContainer', rules, fix_tautomers):
        bonds = self._bonds
        charges = self._charges
        radicals = self._radicals
        calc_implicit = self._calc_implicit

        log = []
        fixed = set()
        flush = False
        for r, (pattern, atom_fix, bonds_fix, any_atoms, is_tautomer) in enumerate(rules):
            if not fix_tautomers and is_tautomer:
                continue
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
                    charges[n] += ch
                    if charges[n] > 4:
                        charges[n] -= ch
                        log.append((tuple(match), r, f'bad charge formed. changes omitted: {pattern}'))
                        break  # skip changes
                    if ir is not None:
                        radicals[n] = ir
                else:
                    for n, m, b in bonds_fix:
                        n = mapping[n]
                        m = mapping[m]
                        hs.add(n)
                        hs.add(m)
                        if m in bonds[n]:
                            bonds[n][m]._Bond__order = b  # noqa
                            if b == 8:
                                # expected original molecule don't contain `any` bonds or these bonds not changed
                                flush = True
                        else:
                            if b != 8:
                                flush = True
                            bonds[n][m] = bonds[m][n] = b = Bond(b)
                            b._attach_graph(self, n, m)
                    log.append((tuple(match), r, str(pattern)))

            if not hs:  # not matched
                continue
            # flush cache only for changed atoms.
            if flush:  # neighbors count changed
                ngb = self.__dict__['__cached_args_method_neighbors']
                for n in hs:
                    try:
                        del ngb[(n,)]
                    except KeyError:
                        pass
                del self.__dict__['bonds_count']
                flush = False
            # need hybridization recalculation
            hyb = self.__dict__['__cached_args_method_hybridization']
            for n in hs:
                try:
                    del hyb[(n,)]
                except KeyError:  # already flushed before
                    pass
            for n in hs:  # hydrogens count recalculation
                calc_implicit(n)
            del self.__dict__['_cython_compiled_structure']
            fixed.update(hs)
        return log, fixed


__all__ = ['Standardize']
