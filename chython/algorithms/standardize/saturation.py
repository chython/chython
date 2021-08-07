# -*- coding: utf-8 -*-
#
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
from collections import defaultdict
from itertools import product
from operator import itemgetter
from random import shuffle
from typing import TYPE_CHECKING, Dict, Optional, Union, List
from ...containers.bonds import Bond
from ...exceptions import ValenceError


if TYPE_CHECKING:
    from chython import MoleculeContainer

charge_priority = {(7, 0): -3, (7, 1): -2, (7, -1): -1,
                   (8, 0): -3, (8, -1): -2, (8, 1): -1}
common_charge_priority = {0: 0, -1: 1, 1: 2, 2: 3, 3: 4, -2: 5, 4: 6, -3: 7, -4: 8}


class Saturation:
    __slots__ = ()

    def saturate(self: 'MoleculeContainer', neighbors_distances: Optional[Dict[int, Dict[int, float]]] = None,
                 reset_electrons: bool = True, expected_charge: int = 0, expected_radicals_count: int = 0,
                 allow_errors: bool = True, logging: bool = False) -> Union[bool, List[str]]:
        """
        Saturate molecules with double and triple bonds and charges and radical states to correct valences of atoms.
        Note: works only with fully explicit hydrogens!

        :param neighbors_distances: If given longest bonds can be removed if need.
        :param reset_electrons: Can change charges and radicals if need.
        :param expected_charge: Reset charge to given. Works only with reset_electrons=True.
        :param expected_radicals_count: Reset radical atoms count to given. Works only with reset_electrons=True.
        :param allow_errors: allow unbalanced result.
        :param logging: return log.
        """
        if any(b != 1 for _, _, b in self.bonds()):
            raise ValenceError('only single bonded skeleton can be saturated')
        atoms = self._atoms
        if not reset_electrons:
            expected_radicals_count = any(self._radicals.values())
            expected_charge = int(self)

        sat, adjacency = _find_possible_valences(atoms, neighbors_distances or self._bonds,
                                                 {x: None for x in self._atoms} if reset_electrons else self._charges,
                                                 {x: None for x in self._atoms} if reset_electrons else self._radicals,
                                                 neighbors_distances is not None)
        charges = {}  # new charge states
        radicals = {}  # new radical states
        bonds = defaultdict(dict)  # new bonds

        seen = set()
        unsaturated = {}
        for n, env in adjacency.items():  # set single bonds in molecule. collect unsaturated atoms
            s = sat[n]
            if len(s) == 1:
                c, r, h = s.pop()
                if not h:
                    seen.add(n)
                    charges[n] = c
                    radicals[n] = r
                    for m in env:
                        if m not in seen:
                            bonds[n][m] = bonds[m][n] = Bond(1)
                else:
                    unsaturated[n] = [(c, r, h)]
            else:
                # radicals have the lowest priority
                # multiple bonds have the highest priority
                # charges have special priority for some atoms and common otherwise: 0>-1>1>2>3>-2>4>-3>-4
                unsaturated[n] = sorted(s, key=lambda x: (x[1], -x[2],
                                                          charge_priority.get((atoms[n].atomic_number, x[0])) or
                                                          common_charge_priority[x[0]]))

        # create graph of unsaturated atoms
        bonds_graph = {n: {m for m in env if m in unsaturated} for n, env in adjacency.items() if n in unsaturated}
        ua, sb, sa = _saturate(bonds_graph, unsaturated)
        for n, m, b in sb:
            bonds[n][m] = bonds[m][n] = Bond(b)
        for n, c, r in sa:
            charges[n] = c
            radicals[n] = r

        combo_ua = []  # possible single atoms electron states
        for n, s in ua.items():
            if len(s) == 1:
                c, r = s[0]
                charges[n] = c
                radicals[n] = r
            elif s:
                combo_ua.append([(n, c, r) for c, r in s])

        log = []
        # try randomly set charges and radicals.
        # first pick required radical states.
        # second try to minimize charge delta.
        if combo_ua:
            need_radical = expected_radicals_count - sum(radicals.values())
            for attempt in range(1, len(combo_ua) + 1):
                shuffle(combo_ua)
                rad = []
                chg = []
                for atom in combo_ua:
                    if len(rad) < need_radical:  # pick radicals
                        r = next((x for x in atom if x[2]), None)
                        if r:  # pick random radical states
                            rad.append(r)
                        else:  # not radical
                            chg.append(atom)
                    else:  # pick not radical states
                        c = [x for x in atom if not x[2]]
                        if len(c) > 1:
                            chg.append(c)
                        elif c:
                            n, c, r = c[0]
                            charges[n] = c
                            radicals[n] = r
                        elif attempt == len(combo_ua):  # all states has radical. balancing impossible
                            chg.append(atom)  # fuck it horse. we in last attempt
                            log.append('Radical state not balanced.')
                        else:  # do next attempt
                            break
                else:
                    for n, c, r in rad:
                        charges[n] = c
                        radicals[n] = r

                    current_charge = sum(charges.values())
                    for x in chg:
                        n, c, r = min(x, key=lambda x: abs(current_charge + x[1]))
                        charges[n] = c
                        radicals[n] = r
                        current_charge += c
                    if sum(charges.values()) == expected_charge:
                        break
            else:
                log.append('Charge state not balanced.')
        elif sum(radicals.values()) != expected_radicals_count or sum(charges.values()) != expected_charge:
            log.append('Charge or radical state not balanced.')

        if not allow_errors and log:
            if logging:
                return log
            return False

        # reset molecule
        self._bonds = dict(bonds)
        self._radicals = radicals
        self._charges = charges
        self._hydrogens = {x: 0 for x in atoms}
        self.flush_cache()
        if logging:
            log.append('saturated successfully')
            return log
        return True


def _find_possible_valences(atoms, neighbors_distances, charges, radicals, allow_deleting=True):
    if allow_deleting:
        possible_bonds = {n: md.copy() for n, md in neighbors_distances.items()}
    else:
        possible_bonds = {n: list(md) for n, md in neighbors_distances.items()}
    while True:
        saturation = defaultdict(set)
        for n, env in possible_bonds.items():
            el = len(env)
            dc = charges[n]
            dr = radicals[n]
            for (charge, is_radical, valence), rules in atoms[n]._compiled_valence_rules.items():
                if valence < el or dc is not None and dc != charge or dr is not None and dr != is_radical:
                    continue  # skip impossible rules
                for _, explicit_dict, h in rules:
                    if explicit_dict:
                        env_atoms = defaultdict(int)
                        for m in env:
                            env_atoms[atoms[m].atomic_number] += 1
                        bonds = 0
                        for (b, a), c in explicit_dict.items():  # stage 1. find explicit valence
                            if env_atoms[a] < c:  # `c` always > 0
                                break  # rule not matched
                            env_atoms[a] -= c
                            bonds += b * c
                        else:  # stage 2. find possible valence
                            unmatched = sum(env_atoms.values())  # atoms outside rule
                            implicit = valence - bonds + h  # implicit H in rule
                            if unmatched:
                                if implicit >= unmatched:
                                    # number of implicit H should be greater or equal to number of neighbors
                                    # excess of implicit H saved as unsaturated atom
                                    saturation[n].add((charge, is_radical, implicit - unmatched))
                            else:  # pattern fully matched. save implicit H count as unsaturated atom.
                                saturation[n].add((charge, is_radical, implicit))
                    elif el == valence:   # unspecific rule. found possible valence
                        saturation[n].add((charge, is_radical, h))
            if n not in saturation:  # valence not found
                break
        else:  # all atoms passed
            break
        if allow_deleting:
            out = max(env.items(), key=itemgetter(1))[0]
            del possible_bonds[out][n]
            del possible_bonds[n][out]
        else:
            raise ValenceError('Structure has invalid atoms neighbors count and electron states')
    return saturation, possible_bonds


def _saturate(bonds, atoms):
    dots = {}
    saturation = []
    electrons = []
    while True:
        # get isolated atoms. atoms should be charged or radical
        new_dots = {n: [(c, r) for c, r, h in atoms[n] if not h] for n, env in bonds.items() if not env}
        for n in new_dots:
            del bonds[n]
        dots.update(new_dots)
        if not bonds:
            break

        try:  # get terminal atom
            n = next(n for n, ms in bonds.items() if len(ms) == 1)
        except StopIteration:
            # get ring or linker atom
            n, _ = min(bonds.items(), key=lambda x: len(x[1]))
            m = bonds[n].pop()
            bonds[m].discard(n)

            for (nc, nr, nh), (i, (mc, mr, mh)) in product(atoms[n], enumerate(atoms[m])):
                if nh == mh:
                    saturation.append((n, m, nh + 1))
                    electrons.append((n, nc, nr))
                    electrons.append((m, mc, mr))

                    for x in bonds.pop(n):
                        saturation.append((n, x, 1))
                        bonds[x].discard(n)
                    for x in bonds.pop(m):
                        saturation.append((m, x, 1))
                        bonds[x].discard(m)
                    break
                elif nh < mh:
                    electrons.append((n, nc, nr))
                    saturation.append((n, m, nh + 1))
                    atoms[m].pop(i)
                    atoms[m].insert(i, (mc, mr, mh - nh))

                    for x in bonds.pop(n):
                        saturation.append((n, x, 1))
                        bonds[x].discard(n)
                    break
                elif nh > mh:
                    electrons.append((m, mc, mr))
                    saturation.append((n, m, mh + 1))
                    atoms[n].pop(i)
                    atoms[n].insert(i, (nc, nr, nh - mh))

                    for x in bonds.pop(m):
                        saturation.append((m, x, 1))
                        bonds[x].discard(m)
                    break
        else:
            m = bonds.pop(n).pop()
            bonds[m].discard(n)

            for (nc, nr, nh), (i, (mc, mr, mh)) in product(atoms[n], enumerate(atoms[m])):
                if nh == mh:
                    saturation.append((n, m, nh + 1))
                    electrons.append((n, nc, nr))
                    electrons.append((m, mc, mr))
                    for x in bonds.pop(m):
                        saturation.append((m, x, 1))
                        bonds[x].discard(m)
                    break
                elif nh < mh and bonds[m]:
                    electrons.append((n, nc, nr))
                    saturation.append((n, m, nh + 1))
                    atoms[m].pop(i)
                    atoms[m].insert(i, (mc, mr, mh - nh))
                    break
            else:
                saturation.append((n, m, 1))
                if not bonds[m]:
                    del bonds[m]
    return dots, saturation, electrons


__all__ = ['Saturation']
