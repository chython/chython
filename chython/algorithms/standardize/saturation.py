# -*- coding: utf-8 -*-
#
#  Copyright 2021-2024 Ramil Nugmanov <nougmanoff@protonmail.com>
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

# atom, charge, unsaturation
tuned_priority = {(7, 0, 0): -3,  # amine
                  (7, 0, 1): -3,  # X=N-X
                  (7, 0, 2): -3,  # X#N
                  (7, 1, 1): -2,  # X=[N+](-X)-X
                  (7, 1, 2): -2,  # X=[N+]=X
                  (7, -1, 0): -1,  # X-[N-]-X
                  (7, -1, 1): -1,  # X=[N-]
                  (8, 0, 0): -2,  # X-O-X
                  (8, 0, 1): -2,  # X=O
                  (8, -1, 0): -1,  # X-[O-]
                  (8, 1, 1): -1,  # X=[O+]-X
                  (16, 0, 0): -3,  # X-S-X
                  (16, 0, 2): -2,  # X-S(-X)(=X)=X
                  (16, 0, 1): -2,  # X-S(-X)(=X)
                  (16, 0, 2): -2,  # X=S=X
                  (16, 0, 0): -1,  # X-S(-X)(-X)-X
                  (16, 0, 1): -1,  # X-[S+](-X)(-X)
                  (16, 1, 1): -1,  # X=[S+]-X
                  }
charge_priority = {0: 0, -1: 1, 1: 2, 2: 3, 3: 4, -2: 5, 4: 6, -3: 7, -4: 8}


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
            expected_radicals_count = sum(a.is_radical for a in atoms.values())
            expected_charge = int(self)

        if reset_electrons:
            charges = {x: None for x in self}
            radicals = {x: None for x in self}
        else:
            charges = {n: a.charge for n, a in self.atoms()}
            radicals =  {n: a.is_radical for n, a in self.atoms()}
        sat, adjacency = _find_possible_valences(atoms, neighbors_distances or self._bonds,
                                                 charges, radicals, neighbors_distances is not None)
        charges = {}  # new charge states
        radicals = {}  # new radical states
        bonds = {n: {} for n in atoms}  # new bonds

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
                # tuned priority
                # multiple bonds have higher priority
                # charges priority: 0>-1>1>2>3>-2>4>-3>-4
                unsaturated[n] = sorted(s, key=lambda x: (x[1],
                                                          tuned_priority.get((atoms[n].atomic_number, x[0], x[2]), 0),
                                                          -x[2], charge_priority[x[0]]))

        log = []
        if (need_radicals := expected_radicals_count - sum(radicals.values())) < 0:
            log.append('Radical state not balanced')
            if not allow_errors:
                if logging:
                    return log
                return False
            need_radicals = 0  # reset to zero

        if unsaturated:
            # create graph of unsaturated atoms
            bonds_graph = {n: {m for m in adjacency[n] if m in unsaturated} for n in unsaturated}
            order = list(unsaturated)
            # try to saturate with different random states
            for _ in range(len(unsaturated)):
                shuffle(order)
                sb, sa, log_ = _saturate({n: bonds_graph[n].copy() for n in order}, unsaturated, need_radicals,
                                         expected_charge - sum(charges.values()))
                if not log_:  # success
                    break
            else:  # failed
                if log_ == 1:
                    log.append('Charge state not balanced')
                elif log_ == 2:
                    log.append('Radical state not balanced')
                else:
                    log.append('Charge state not balanced')
                    log.append('Radical state not balanced')
                if not allow_errors:  # all attempts failed
                    if logging:
                        return log
                    return False

            for n, m, b in sb:
                bonds[n][m] = bonds[m][n] = Bond(b)
            for n, c, r in sa:
                charges[n] = c
                radicals[n] = r
        elif expected_charge != sum(charges.values()):  # check charge for saturated case
            log.append('Charge state not balanced')
            if not allow_errors:
                if logging:
                    return log
                return False
        # reset molecule
        self._bonds = bonds
        for n, r in radicals.items():
            atoms[n]._is_radical = r
        for n, c in charges.items():
            atoms[n]._charge = c
        for a in atoms.values():
            a._implicit_hydrogens = 0  # reset invalid hydrogens counts.
        self.flush_cache()
        self.calc_labels()
        if logging:
            if not log:  # check for errors
                log.append('Saturated successfully')
            else:
                log.append('Saturated with errors')
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
            env_atoms = None
            el = len(env)
            dc = charges[n]
            dr = radicals[n]
            for charge, is_radical, valence, implicit, explicit_dict in atoms[n]._compiled_saturation_rules:
                if valence < el or dc is not None and dc != charge or dr is not None and dr != is_radical:
                    continue  # skip impossible rules
                if explicit_dict:
                    if env_atoms is None:  # lazy caching
                        env_atoms = defaultdict(int)
                        for m in env:
                            env_atoms[atoms[m].atomic_number] += 1
                    env_atoms_copy = env_atoms.copy()
                    for (b, a), c in explicit_dict.items():  # stage 1. find explicit valence
                        if env_atoms_copy[a] < c:  # `c` always > 0
                            break  # rule not matched
                        env_atoms_copy[a] -= c
                    else:  # stage 2. find possible valence
                        if unmatched := sum(env_atoms_copy.values()):  # number of atoms outside rule
                            if implicit >= unmatched:
                                # number of implicit H should be greater or equal to number of neighbors
                                saturation[n].add((charge, is_radical, valence - el))
                        else:  # pattern fully matched. difference bw valence and connectivity is unsaturation.
                            saturation[n].add((charge, is_radical, valence - el))
                else:   # unspecific rule. found possible valence
                    saturation[n].add((charge, is_radical, valence - el))
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


def _saturate(bonds, atoms, expected_radicals_count, expected_charge):
    atoms = {k: v.copy() for k, v in atoms.items()}
    dots = []
    saturation = []
    electrons = []
    while True:
        # get isolated atoms. atoms should be charged or radical
        to_del = []
        for n, env in bonds.items():
            if not env:
                es = [(n, c, r) for c, r, h in atoms[n] if not h]
                if not es:
                    raise ValenceError('Saturation impossible. '
                                       f"Isolated atom ({n}) doesn't have appropriate charge-radical state")
                to_del.append(n)
                dots.append(es)
        for n in to_del:
            del bonds[n]
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

    combo_ua = []  # possible single atoms electron states
    for s in dots:
        if len(s) == 1:
            electrons.extend(s)
        elif s:
            combo_ua.append(s)

    # if < 0 - we already in bad situation
    # if > 0 - we need more radicals
    need_radical = expected_radicals_count - sum(x for _, _, x in electrons)
    need_charge = expected_charge - sum(x for _, x, _ in electrons)
    if combo_ua:
        # try randomly set charges and radicals.
        # first pick required radical states.
        # second try to minimize charge delta.
        for attempt in range(1, len(combo_ua) + 1):
            shuffle(combo_ua)
            charges_radicals = []
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
                        charges_radicals.extend(c)
                    elif attempt == len(combo_ua):  # all states has radical. balancing impossible
                        chg.append(atom)  # fuck it horse. we in last attempt
                    else:  # do next attempt
                        break
            else:
                charges_radicals.extend(rad)
                current_charge = need_charge - sum(x for _, x, _ in charges_radicals)
                current_radical = need_radical - len(rad)
                for x in chg:
                    n, c, r = min(x, key=lambda x: abs(current_charge - x[1]))
                    charges_radicals.append((n, c, r))
                    current_charge -= c
                    current_radical -= r

                if current_radical:  # radical unbalanced
                    if current_charge:
                        log = 3
                    else:
                        log = 2
                elif current_charge:
                    log = 1
                else:  # balanced!
                    log = 0
                    break

        electrons.extend(charges_radicals)
    elif need_radical:
        log = 3 if need_charge else 2
    elif need_charge:
        log = 1
    else:
        log = 0
    return saturation, electrons, log


__all__ = ['Saturation']
