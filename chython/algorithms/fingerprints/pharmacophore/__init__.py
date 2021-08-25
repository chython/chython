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

        hydrogens = self._hydrogens
        atoms = self._atoms
        bonds = self._bonds
        charges = self._charges

        acceptors_donors = {n for n, a in atoms.items() if a.atomic_number in {7, 8, 16} and a.neighbors <= 3}
        halogens = {idx: num for idx, atom in atoms.items()
                    if (num := atom.atomic_number) in {17, 35, 53} and self.neighbors(idx) < 2}
        charged = {n: c for n, c in atoms.items() if c}

        bins = [
            (acceptors := {n for n in acceptors_donors if not hydrogens[n]}),  # H acceptors
            acceptors_donors - acceptors,  # H donors

            {n for q in acid for dct in q.get_mapping(self, automorphism_filter=False) for n in dct.values()
             if hydrogens[n]},  # acidic

            {n for q in base for dct in q.get_mapping(self, automorphism_filter=False) for n in dct.values()
             if not hydrogens[n]},  # basic

            (positive := {n for n, c in charged.items() if c > 0}),  # positive charged
            charged.keys() - positive,  # negative charged

            {n for n, a in atoms.items() if a.atomic_number == 9},  # fluorine
            halogens,  # halogens except fluorine

            {n for n, a in atoms.items() if a.atomic_number == 6 and
             all(atoms[x].atomic_number in {1, 5, 6, 9, 14} for x in bonds[n])},  # hydrophobic

            {n for n in atoms if self.hybridization(n) == 4},  # aromatic
        ]

        out = {idx: 0 for idx in self._atoms}
        for pos, bin_ in zip((512, 256, 128, 64, 32, 16, 8, 4, 2, 1), bins):
            for idx in bin_:
                out[idx] |= pos

        return out


__all__ = ['Pharmacophore']
