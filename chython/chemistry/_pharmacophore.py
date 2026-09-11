# -*- coding: utf-8 -*-
#
#  Copyright 2026 Ramil Nugmanov <nougmanoff@protonmail.com>
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
"""Per-atom pharmacophore invariants over `tables/pharmacophore.tsv`.

Six 2D types after Kutlushina, Khakimova, Madzhidov, Polishchuk, Molecules 2018, 23, 3094.
`donor` and `acceptor` are read through `hbond_atoms()` rather than duplicated here; the four
remaining types (`positive`, `negative`, `aromatic`, `hydrophobe`) come from `pharmacophore.tsv`.
"""
from ..core._core import require_numpy
from ._counts import hbond_atoms
from ._tables import pharmacophore_rules_by_role


PH_DONOR = 1
PH_ACCEPTOR = 2
PH_POSITIVE = 4
PH_NEGATIVE = 8
PH_AROMATIC = 16
PH_HYDROPHOBE = 32

PH_TYPES = ('donor', 'acceptor', 'positive', 'negative', 'aromatic', 'hydrophobe')
_PH_BITS = {'donor': PH_DONOR, 'acceptor': PH_ACCEPTOR, 'positive': PH_POSITIVE,
            'negative': PH_NEGATIVE, 'aromatic': PH_AROMATIC, 'hydrophobe': PH_HYDROPHOBE}


def pharmacophore_atoms(molecule) -> dict:
    """Stable ids per 2D pharmacophore feature type.  Six keys, always all six, possibly empty."""
    out = {'donor': hbond_atoms(molecule, 'donor'),
           'acceptor': hbond_atoms(molecule, 'acceptor')}
    for role, rows in pharmacophore_rules_by_role().items():
        found = set()
        for row in rows:
            subject = row.numbers[1]      # stable id of :1, from compile_smarts at load time
            for mapping in row.query.get_mapping(molecule):
                found.add(mapping[subject])
        out[role] = frozenset(found)
    return out


def pharmacophore_invariants(molecule):
    """Per-atom pharmacophore feature bitmask, `uint32[atom_count]`, in `atom_numbers` order.

    Six 2D types after Kutlushina, Khakimova, Madzhidov, Polishchuk, Molecules 2018, 23, 3094.  An
    atom with no feature is 0, and that is deliberate: the point of a pharmacophore fingerprint is to
    lose element identity.  Pass it straight to any `morgan_*` or `linear_*` method as `invariants=`.

    NUMPY IS IMPORTED HERE AND NOT AT MODULE LEVEL, and this line is the reason the whole dependency
    was not optional.  It is `chython[ml]`, and this module is reached eagerly from
    `chython.chemistry.__init__`, so a module-level `from numpy import uint32, zeros` made `import
    chython` fail outright on a minimal install -- one featurizer most callers never touch, deciding
    for the entire façade.  Before numpy was optional at all the same line quietly made every `import
    chython` pay numpy's import cost, and falsified the measured claim in the core's binder docstring.
    `pharmacophore_atoms` above answers stable ids and needs no array, so only this function pays.

    `require_numpy()` first, so the failure is the core's one message naming the extra rather than a
    bare "No module named 'numpy'" from two lines down.
    """
    require_numpy()
    from numpy import uint32, zeros

    index = {i: n for n, i in enumerate(molecule.atom_numbers)}
    out = zeros(molecule.atom_count, dtype=uint32)
    for role, ids in pharmacophore_atoms(molecule).items():
        bit = _PH_BITS[role]
        for i in ids:
            out[index[i]] |= bit
    return out
