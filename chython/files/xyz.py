# -*- coding: utf-8 -*-
#
#  Copyright 2020-2024 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from numpy import array
from typing import Sequence, Tuple, Optional
from ..containers import MoleculeContainer
from ..periodictable import Element


def xyz(matrix: Sequence[Tuple[str, float, float, float]], charge=0, radical=0, radius_multiplier=1.25,
        atom_charge: Optional[Sequence[int]] = None, _cls=MoleculeContainer) -> MoleculeContainer:
    from ._xyz import possible_bonds  # windows?

    if atom_charge and len(atom_charge) != len(matrix):
        raise ValueError('atom_charge should be None or the same size as matrix')

    mol = _cls()
    conformer = {}
    mol._conformers = [conformer]

    atoms = mol._atoms
    bonds = mol._bonds
    for n, (a, x, y, z) in enumerate(matrix, 1):
        atoms[n] = Element.from_symbol(a)(x=x, y=y, implicit_hydrogens=0)
        bonds[n] = {}
        conformer[n] = (x, y, z)

    if atom_charge is not None and None not in atom_charge:
        for n, c in enumerate(atom_charge, 1):
            atoms[n]._charge = c
        charge = sum(atom_charge)

    mol.calc_labels()
    pb = possible_bonds(array(list(conformer.values())), array([a.atomic_radius for a in atoms.values()]),
                        radius_multiplier)

    log = mol.saturate(pb, expected_charge=charge, expected_radicals_count=radical, logging=True)
    mol.meta['saturation_log'] = log
    return mol


def xyz_file(data) -> MoleculeContainer:
    data = data.splitlines()

    size = int(data[0])
    charge = 0
    radical = 0

    for x in data[1].split():
        if x.startswith('charge='):
            charge = int(x[7:])
        elif x.startswith('radical='):
            radical = int(x[8:])

    _xyz = []
    for n, line in enumerate(data[2: size + 2]):
        symbol, x, y, z = line.split()
        _xyz.append((symbol, float(x), float(y), float(z)))

    return xyz(_xyz, charge, radical)


__all__ = ['xyz', 'xyz_file']
