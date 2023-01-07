# -*- coding: utf-8 -*-
#
#  Copyright 2020-2023 Ramil Nugmanov <nougmanoff@protonmail.com>
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
    mol._conformers.append(conformer)

    atoms = mol._atoms
    bonds = mol._bonds
    plane = mol._plane
    hydrogens = mol._hydrogens
    radicals = mol._radicals
    for n, (a, x, y, z) in enumerate(matrix, 1):
        atoms[n] = atom = Element.from_symbol(a)()
        atom._attach_graph(mol, n)
        bonds[n] = {}
        plane[n] = (x, y)
        conformer[n] = (x, y, z)
        hydrogens[n] = 0  # implicit hydrogens not supported.
        radicals[n] = False  # set default value

    if atom_charge is None or None in atom_charge:
        mol._charges = {n: 0 for n in atoms}  # reset charges
    else:
        mol._charges = dict(enumerate(atom_charge, 1))
        charge = sum(atom_charge)

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
