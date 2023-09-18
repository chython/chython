# -*- coding: utf-8 -*-
#
#  Copyright 2023 Timur Gimadiev <timur.gimadiev@gmail.com>
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
from joblib import Parallel, delayed
from numpy import zeros, uint8
from typing import TYPE_CHECKING, Union

if TYPE_CHECKING:
    from chython import MoleculeContainer, CGRContainer


class Fragmentor:
    """
    Fragmentor object can store all seen fragments by .fit() procedure and then produce descriptors
    using dictionary of seen fragments.
    """
    def __init__(self, min_radius: int = 1, max_radius: int = 4, circus: bool = False):
        self.fragments = set()
        self.fragments_map = {}
        self.smi2num = {}
        self.min_radius = min_radius
        self.max_radius = max_radius
        self.circus = circus

    def get_keys_linear(self, molecule: Union['MoleculeContainer', 'CGRContainer']):
        """
        returns text strings of fragments within min_radius and max_radius.

        :param molecule MolecularContainer of molecule
        """
        desc = molecule.linear_fragments_smiles(self.min_radius, self.max_radius)
        return desc

    def get_keys_circus(self, molecule: Union['MoleculeContainer', 'CGRContainer']):
        desc = molecule.circus_hash_dict(self.min_radius, self.max_radius)
        return desc

    def get_desciptor_vector(self, molecule, max_count=0):
        if self.circus:
            keys_getter = self.get_keys_circus
        else:
            keys_getter = self.get_keys_linear
        desc = keys_getter(molecule)
        fp = zeros(len(self.fragments), dtype=uint8)
        for key, val in desc.items():
            if key in self.fragments:
                if max_count > 0:
                    fp[self.smi2num[key]] = len(val) if len(val) < max_count else max_count
                else:
                    fp[self.smi2num[key]] = len(val)
        return fp

    def fit(self, molecules: list[Union['MoleculeContainer', 'CGRContainer']], n_jobs: int = 1):
        """
        Method collects all unique fragments form the set of molecules, but before clears up
        remembered fragments.

        :param molecules List[MolecularContainer] of molecules for collecting fragments
        :param n_jobs number of parallel jobs
        """
        self.clear()
        self.partial_fit(molecules, n_jobs=n_jobs)

    def partial_fit(self, molecules: list[Union['MoleculeContainer', 'CGRContainer']],
                    n_jobs: int = 1):
        """
        Method collects all unique fragments form the set of molecules and append them to the known
        fragments. Designed for feeding data by portions.

        :param molecules List[MolecularContainer] of molecules for collecting fragments
        :param n_jobs Number of jobs in parallel
        """
        if self.circus:
            keys_getter = self.get_keys_circus
        else:
            keys_getter = self.get_keys_linear
        if n_jobs < 2:
            for mol in molecules:
                self.fragments.update(keys_getter(mol))
        else:
            for calc in Parallel(n_jobs=n_jobs)(delayed(keys_getter)(molecule) for molecule in
                                                molecules):
                self.fragments.update(calc)
        self.fragments_map = {num: val for num, val in enumerate(self.fragments)}
        self.smi2num = {val: num for num, val in enumerate(self.fragments)}

    def transform(self, molecules: list[Union['MoleculeContainer', 'CGRContainer']],
                  n_jobs: int = 1, max_count: int = 0):
        """
         Method produce list of descriptor vectors according to dictionary of seen fragments

          :param molecules List[MolecularContainer] of molecules for collecting fragments
          :param n_jobs Number of jobs in parallel
          :param max_count max value of each position in descriptor vector
        """
        results = Parallel(n_jobs=n_jobs)(
            delayed(self.get_desciptor_vector)(molecule, max_count) for
            molecule in molecules)
        return results

    def clear(self):
        self.fragments = set()
        self.fragments_map = {}
        self.smi2num = {}

    @property
    def fragments_map_smi_ordered(self):
        return {x: self.fragments_map[x] for x in sorted(self.fragments_map, key=lambda x:
                self.fragments_map[x])}
