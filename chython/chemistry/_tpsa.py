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
"""Topological polar surface area (TPSA) over `tables/tpsa.tsv`.

Ertl, Rohde, Selzer, J. Med. Chem. 2000, 43, 3714.
"""
from ._standardize import LogRecord
from ._tables import TPSA_CLASSES, first_match, tpsa_rules
from ..core import recording


_POLAR = {'NO': frozenset((7, 8)), 'SP': frozenset((15, 16))}  # the elements each class types


def tpsa_contributions(molecule, *, sulfur_phosphorus=False) -> dict:
    """Per-atom TPSA contribution in A^2, keyed by stable id.

    Ertl, Rohde, Selzer, J. Med. Chem. 2000, 43, 3714.  An atom with no contribution is absent from the
    dict, not present as zero.  A polar atom that matches no published environment contributes nothing
    and puts one `tpsa:unmatched` record on `molecule.log`.
    """
    # `TPSA_CLASSES` rather than a second copy of the class names, so a class added to `tpsa.tsv` needs
    # no edit here.  `_POLAR` cannot be derived and is pinned against the same constant by
    # `test_tpsa.py`, so a class arriving without its elements fails a test rather than logging noise.
    wanted = TPSA_CLASSES if sulfur_phosphorus else ('NO',)
    rows = [r for r in tpsa_rules() if r.element_class in wanted]
    matched = first_match(rows, molecule)
    out = {}
    for i, row in matched.items():
        if row.contribution:
            out[i] = row.contribution
    # the element set comes from `wanted`, the same tuple that selected the rows: under
    # `sulfur_phosphorus=False` no SP row is loaded, so reporting every S and P would state the
    # caller's own choice back as a defect.
    polar = frozenset().union(*(_POLAR[c] for c in wanted))
    with recording(molecule, stage='tpsa') as log:
        for atom in molecule.atoms():
            # `atom.element` is the atomic number; `atom.atomic_symbol` is the string.
            if atom.element in polar and atom.n not in matched:
                log.append(LogRecord('tpsa:unmatched', (atom.n,),
                                     'no published TPSA environment matches this atom; it '
                                     'contributes zero'))
    return out


def tpsa(molecule, *, sulfur_phosphorus=False) -> float:
    """Topological polar surface area in A^2.

    Ertl, Rohde, Selzer, J. Med. Chem. 2000, 43, 3714.  Nitrogen and oxygen only, as published;
    `sulfur_phosphorus=True` adds the paper's optional S and P contributions, which is a different
    published quantity rather than a refinement of this one.
    """
    return sum(tpsa_contributions(molecule, sulfur_phosphorus=sulfur_phosphorus).values())
