# -*- coding: utf-8 -*-
#
#  Copyright 2021 Aleksandr Sizov <murkyrussian@gmail.com>
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
from functools import cached_property
from typing import Dict, TYPE_CHECKING
from ...tautomers._acid import rules as acidic_rules
from ...tautomers._base import rules as basic_rules


if TYPE_CHECKING:
    from ....containers import MoleculeContainer


class Pharmacophore:
    __slots__ = ()
    """
    Matcher of atoms in an molecule to some integer number.
    Each integer represents 10 bit, each bit is enable/disable some pharmacophoric feature:
    | H acceptor | H donor | acidic | basic | positive charged | negative charged | ->
    -> | fluorine | any halogen atom | hydrophobic | aromatic |
    After the constructing, binary number will be converted to 10-based integer
    """

    @cached_property
    def pharmacophores(self: 'MoleculeContainer') -> Dict[int, int]:
        """
        Match id of each atom of some molecule with integer number which define some features using in FCFP
        """
        acid = acidic_rules[:-1]  # discard a rule with halogen acids
        base = basic_rules[:3] + basic_rules[4:]  # there is no need a rule with halogen ions

        atoms = self._atoms
        bonds = self._bonds
        hydrogens = self._hydrogens
        neighbors = self.neighbors
        hybridization = self.hybridization

        acceptors_donors = {n for n, a in atoms.items()
                            if (num := a.atomic_number) in (7, 8) or  # any N,O
                            num == 16 and hybridization(n) == 1 and neighbors(n) <= 2}  # or S-sp3
        halogens = {n for n, atom in atoms.items() if atom.atomic_number in (17, 35, 53) and neighbors(n) < 2}
        charged = {n: c for n, c in self._charges.items() if c}

        bins = (
            acceptors_donors,  # H acceptors
            {n for n in acceptors_donors if hydrogens[n]},  # H donors

            {dct[1] for q in acid for dct in q.get_mapping(self, automorphism_filter=False)},  # acidic
            {dct[1] for q in base for dct in q.get_mapping(self, automorphism_filter=False)},  # basic

            (positive := {n for n, c in charged.items() if c > 0}),  # positive charged
            charged.keys() - positive,  # negative charged

            {n for n, a in atoms.items() if a.atomic_number == 9},  # fluorine
            halogens,  # halogens except fluorine

            {n for n, a in atoms.items() if a.atomic_number == 6 and
             all(atoms[x].atomic_number in (1, 5, 6, 9, 14) for x in bonds[n])},  # hydrophobic

            {n for n in atoms if hybridization(n) == 4}  # aromatic
        )

        out = {idx: 0 for idx in atoms}
        for pos, bin_ in zip((512, 256, 128, 64, 32, 16, 8, 4, 2, 1), bins):
            for idx in bin_:
                out[idx] |= pos
        return out


__all__ = ['Pharmacophore']
