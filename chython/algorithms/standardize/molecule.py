# -*- coding: utf-8 -*-
#
#  Copyright 2018-2021 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from ...exceptions import ValenceError
from ...periodictable import H


if TYPE_CHECKING:
    from chython import MoleculeContainer


class StandardizeMolecule:
    __slots__ = ()

    def canonicalize(self: 'MoleculeContainer', *, logging=False) -> \
            Union[bool, List[Tuple[Tuple[int, ...], int, str]]]:
        """
        Convert molecule to canonical forms of functional groups and aromatic rings without explicit hydrogens.

        :param logging: return log.
        """
        k = self.kekule()
        s = self.standardize(fix_stereo=False, logging=logging)
        h = self.implicify_hydrogens(fix_stereo=False)
        t = self.thiele()
        c = self.standardize_charges(prepare_molecule=False, logging=logging)
        if logging:
            if k:
                s.insert(0, ((), -1, 'kekulized'))
            if h:
                s.append(((), -1, 'implicified'))
            if t:
                s.append(((), -1, 'aromatized'))
            if c:
                s.append((tuple(c), -1, 'recharged'))
            return s
        return k or s or h or t or c

    def standardize(self: Union['MoleculeContainer', 'StandardizeMolecule'], *, fix_stereo=True, logging=False) -> \
            Union[bool, List[Tuple[Tuple[int, ...], int, str]]]:
        """
        Standardize functional groups. Return True if any non-canonical group found.

        :param logging: return list of fixed atoms with matched rules.
        """
        neutralized = self.__fix_resonance(logging=logging)
        if neutralized:
            self.flush_cache()
        log = self.__standardize(double_rules)
        log.extend(self.__standardize(double_rules))  # double shot rules for overlapped groups
        log.extend(self.__standardize(single_rules))
        log.extend(self.__standardize(metal_rules))  # metal-organics fix
        if log:
            self.flush_cache()
        if fix_stereo:
            self.fix_stereo()
        if logging:
            if neutralized:
                log.append((tuple(neutralized), -1, 'resonance'))
            return log
        return neutralized or bool(log)

    def standardize_charges(self: 'MoleculeContainer', *, fix_stereo=True,
                            logging=False, prepare_molecule=True) -> Union[bool, List[int]]:
        """
        Set canonical positions of charges in heterocycles and some nitrogen compounds.

        :param logging: return list of changed atoms.
        :param prepare_molecule: do thiele procedure.
        """
        changed: List[int] = []
        bonds = self._bonds
        hydrogens = self._hydrogens
        charges = self._charges

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
        if changed:
            self.flush_cache()  # clear cache
            if fix_stereo:
                self.fix_stereo()
            if logging:
                return changed
            return True
        if logging:
            return []
        return False

    def remove_hydrogen_bonds(self: 'MoleculeContainer', *, keep_to_terminal=True, fix_stereo=True) -> int:
        """Remove hydrogen bonds marked with 8 (any) bond

        :param keep_to_terminal: Keep any bonds to terminal hydrogens
        :return: removed bonds count
        """
        bonds = self._bonds
        hg = defaultdict(set)
        for n, atom in self._atoms.items():
            if atom.atomic_number == 1:
                for m, b in bonds[n].items():
                    if b.order == 8:
                        hg[n].add(m)

        if keep_to_terminal:
            for n, ms in hg.items():
                if len(bonds[n]) == len(ms):  # H~A or A~H~A etc case
                    m = ms.pop()
                    if m in hg:  # H~H case
                        hg[m].discard(n)

        seen = set()
        c = 0
        for n, ms in hg.items():
            seen.add(n)
            for m in ms:
                if m in seen:
                    continue
                c += 1
                del bonds[n][m], bonds[m][n]
        if c:
            self.flush_cache()
            if fix_stereo:
                self.fix_stereo()
        return c

    def implicify_hydrogens(self: 'MoleculeContainer', *, fix_stereo=True) -> int:
        """
        Remove explicit hydrogen if possible.
        Works only with Kekule forms of aromatic structures.
        Keeps isotopes of hydrogen.

        :return: number of removed hydrogens
        """
        atoms = self._atoms
        charges = self._charges
        radicals = self._radicals
        bonds = self._bonds
        explicit = defaultdict(list)
        for n, atom in atoms.items():
            if atom.atomic_number == 1 and (atom.isotope is None or atom.isotope == 1):
                if len(bonds[n]) > 1:
                    raise ValenceError(f'Hydrogen atom {{{n}}} has invalid valence. Try to use remove_hydrogen_bonds()')
                for m, b in bonds[n].items():
                    if b.order == 1:
                        if atoms[m].atomic_number != 1:  # not H-H
                            explicit[m].append(n)
                    elif b.order != 8:
                        raise ValenceError(f'Hydrogen atom {{{n}}} has invalid valence {{{b.order}}}.')

        to_remove = set()
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
                    rules = atom.valence_rules(charge, is_radical, explicit_sum)
                except ValenceError:
                    break
                if any(s.issubset(explicit_dict) and all(explicit_dict[k] >= c for k, c in d.items()) and h >= i
                       for s, d, h in rules):
                    to_remove.update(hi)
                    break
        for n in to_remove:
            self.delete_atom(n)
        if to_remove and fix_stereo:
            self.fix_stereo()
        return len(to_remove)

    def explicify_hydrogens(self: 'MoleculeContainer', *, fix_stereo=True, start_map=None, return_maps=False) -> \
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
                raise ValenceError(f'atom {{{n}}} has valence error')

        if return_maps:
            log = []
        if to_add:
            bonds = self._bonds
            m = start_map
            for n in to_add:
                m = self.add_atom(H(), _map=m)
                bonds[n][m] = bonds[m][n] = Bond(1)
                hydrogens[n] = 0
                if return_maps:
                    log.append((n, m))
                m += 1

            if fix_stereo:
                self.fix_stereo()
            if return_maps:
                return log
            return len(to_add)
        if return_maps:
            return log
        return 0

    def check_valence(self: 'MoleculeContainer') -> List[int]:
        """
        Check valences of all atoms.

        Works only on molecules with aromatic rings in Kekule form.
        :return: list of invalid atoms
        """
        atoms = self._atoms
        charges = self._charges
        radicals = self._radicals
        bonds = self._bonds
        hydrogens = self._hydrogens
        errors = set(atoms)
        for n, atom in atoms.items():
            charge = charges[n]
            is_radical = radicals[n]
            explicit_sum = 0
            explicit_dict = defaultdict(int)
            for m, bond in bonds[n].items():
                order = bond.order
                if order == 4:  # aromatic rings not supported
                    break
                elif order != 8:  # any bond used for complexes
                    explicit_sum += order
                    explicit_dict[(order, atoms[m].atomic_number)] += 1
            else:
                try:
                    rules = atom.valence_rules(charge, is_radical, explicit_sum)
                except ValenceError:
                    pass
                else:
                    hs = hydrogens[n]
                    for s, d, h in rules:
                        if h == hs and s.issubset(explicit_dict) and all(explicit_dict[k] >= c for k, c in d.items()):
                            errors.discard(n)
                            break
        return list(errors)

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

    def __standardize(self: 'MoleculeContainer', rules):
        bonds = self._bonds
        charges = self._charges
        radicals = self._radicals

        log = []
        flush = False
        for r, (pattern, atom_fix, bonds_fix) in enumerate(rules):
            hs = set()
            seen = set()
            # collect constrain Any-atoms
            any_atoms = [n for n, a in pattern.atoms() if a.atomic_symbol == 'A' and n not in atom_fix]
            # AnyMetal can match multiple times
            any_atoms.extend(n for n, a in pattern.atoms() if a.atomic_symbol == 'M')
            for mapping in pattern.get_mapping(self, optimize=False, automorphism_filter=False):
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
                            bonds[n][m]._Bond__order = b
                            if b == 8:
                                # expected original molecule don't contain `any` bonds or these bonds not changed
                                flush = True
                        else:
                            bonds[n][m] = bonds[m][n] = Bond(b)
                            if b != 8:
                                flush = True
                    log.append((tuple(match), r, str(pattern)))
            # flush cache only for changed atoms.
            if flush:  # neighbors count changed
                if '__cached_args_method_neighbors' in self.__dict__:
                    ngb = self.__dict__['__cached_args_method_neighbors']
                    for n in hs:
                        try:
                            del ngb[(n,)]
                        except KeyError:
                            pass
                flush = False
            # need hybridization recalculation
            if '__cached_args_method_hybridization' in self.__dict__:
                hyb = self.__dict__['__cached_args_method_hybridization']
                for n in hs:
                    try:
                        del hyb[(n,)]
                    except KeyError:  # already flushed before
                        pass
            for n in hs:  # hydrogens count recalculation
                self._calc_implicit(n)
        return log

    def __fix_resonance(self: Union['MoleculeContainer', 'StandardizeMolecule'], *, logging=False) -> \
            Union[bool, List[int]]:
        """
        Transform biradical or dipole resonance structures into neutral form. Return True if structure form changed.

        :param logging: return list of changed atoms.
        """
        atoms = self._atoms
        charges = self._charges
        radicals = self._radicals
        bonds = self._bonds
        entries, exits, rads, constrains = self.__entries()
        hs = set()
        while len(rads) > 1:
            n = rads.pop()
            for path in self.__find_delocalize_path(n, rads, constrains):
                l, m, b = path[-1]
                if b == 1:  # required pi-bond
                    continue
                try:
                    atoms[m].valence_rules(charges[m], False, sum(int(y) for x, y in bonds[m].items() if x != l) + b)
                except ValenceError:
                    continue
                radicals[n] = radicals[m] = False
                rads.discard(m)
                hs.add(n)
                hs.update(x for _, x, _ in path)
                for n, m, b in path:
                    bonds[n][m]._Bond__order = b
                break  # path found
            # path not found. atom n keep as is
        while entries and exits:
            n = entries.pop()
            for path in self.__find_delocalize_path(n, exits, constrains):
                l, m, b = path[-1]
                try:
                    atoms[m].valence_rules(charges[m] - 1, radicals[m],
                                           sum(int(y) for x, y in bonds[m].items() if x != l) + b)
                except ValenceError:
                    continue
                charges[n] = charges[m] = 0
                exits.discard(m)
                hs.add(n)
                hs.update(x for _, x, _ in path)
                for n, m, b in path:
                    bonds[n][m]._Bond__order = b
                break  # path from negative atom to positive atom found.
            # path not found. keep negative atom n as is
        if hs:
            for n in hs:
                self._calc_implicit(n)
            if logging:
                return list(hs)
            return True
        if logging:
            return []
        return False

    def __find_delocalize_path(self: 'MoleculeContainer', start, finish, constrains):
        bonds = self._bonds
        stack = [(start, n, 0, b.order + 1) for n, b in bonds[start].items() if n in constrains and b.order < 3]
        path = []
        seen = {start}
        while stack:
            last, current, depth, order = stack.pop()
            if len(path) > depth:
                seen.difference_update(x for _, x, _ in path[depth:])
                path = path[:depth]

            path.append((last, current, order))

            if current in finish:
                if depth:  # one bonded ignored. we search double bond transfer! A=A-A >> A-A=A.
                    yield path
                continue  # stop grow

            depth += 1
            seen.add(current)
            diff = -1 if depth % 2 else 1
            stack.extend((current, n, depth, bo) for n, b in bonds[current].items()
                         if n not in seen and n in constrains and 1 <= (bo := b.order + diff) <= 3)

    def __entries(self: 'MoleculeContainer'):
        charges = self._charges
        radicals = self._radicals
        atoms = self._atoms
        bonds = self._bonds

        transfer = set()
        entries = set()
        exits = set()
        rads = set()
        for n, a in atoms.items():
            if a.atomic_number not in {5, 6, 7, 8, 14, 15, 16, 33, 34, 52}:
                # filter non-organic set, halogens and aromatics
                continue
            if charges[n] == -1:
                if len(bonds[n]) == 4 and a.atomic_number == 5:  # skip boron
                    continue
                entries.add(n)
            elif charges[n] == 1:
                if len(bonds[n]) == 4 and a.atomic_number == 7:  # skip ammonia
                    continue
                exits.add(n)
            elif radicals[n]:
                rads.add(n)
            transfer.add(n)
        return entries, exits, rads, transfer


__all__ = ['StandardizeMolecule']
