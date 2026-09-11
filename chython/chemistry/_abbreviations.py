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
"""Expansion of a contracted group drawn as one labelled atom.

A file that draws one atom and writes `OMe` on it has stated a methoxy group.  The reader stores the
label as the atom's alias and the atom as an R -- neither invents a structure -- and this pass turns the
ones `tables/abbreviations.tsv` knows into atoms.  A label the table does not know is left alone with its
alias, which is what makes the table safe to grow.

THE LABELLED ATOM IS TRANSMUTED, NOT REPLACED.  `set_element` keeps its stable id, its bonds and its
neighbours' parities; deleting it and adding the fragment's attachment atom instead would reorder a
neighbouring stereocentre's references and invalidate a parity nothing here restated.
"""
from ._tables import abbreviation_row
from ..core import H_UNKNOWN, LogRecord, MoleculeContainer, REFUSED, REPAIRED, recording


__all__ = ['expand_abbreviations']


#: Neither exactly one neighbour nor a single bond to it: the marker in the table states the one single
#: bond a contracted group hangs by, and a site that is not that is not this group.
_RULE_ATTACHMENT = 'abbreviations:attachment'

#: The file stated a charge, a radical or an isotope on the labelled atom, and the label states one too.
#: Nothing here ranks the two, so the site keeps the label.
_RULE_STATED = 'abbreviations:stated-atom'


def expand_abbreviations(molecule: MoleculeContainer) -> bool:
    """Replace every atom whose alias names a row of `tables/abbreviations.tsv` with that fragment.

    All-or-nothing per site.  A site that survives every check is expanded and its alias dropped; a site
    that fails one is left exactly as the file drew it, alias included, and the reason is a REFUSED
    record.  Returns True when at least one site was expanded.

    | Outcome | Rule | Severity |
    | --- | --- | --- |
    | expanded | the row's own id, `abbreviations:OMe` | REPAIRED |
    | not one single bond to one neighbour | `abbreviations:attachment` | REFUSED |
    | charge, radical or isotope stated on the labelled atom | `abbreviations:stated-atom` | REFUSED |

    The grafted atoms take the labelled atom's coordinates, so a record with a depiction needs
    `clean2d()` afterwards; the message says so where there was one to disturb.
    """
    with recording(molecule, stage='abbreviations') as log:
        aliases = molecule.aliases
        if not aliases:
            return False

        sites = []
        for n, text in aliases.items():
            try:
                label = text.decode('utf-8')
            except UnicodeDecodeError:                 # not a spelling any table holds
                continue
            row = abbreviation_row(label)
            if row is None:
                continue

            neighbors = list(molecule.neighbors_of(n))
            if len(neighbors) != 1 or molecule.order_of(n, neighbors[0]) != 1:
                log.append(LogRecord(_RULE_ATTACHMENT, (n,),
                                     f'atom {n}: {label} hangs by one single bond and this atom has '
                                     f'{len(neighbors)} neighbours; the label is kept', REFUSED))
                continue
            atom = molecule.atom(n)
            if atom.charge or atom.is_radical or atom.isotope:
                log.append(LogRecord(_RULE_STATED, (n,),
                                     f'atom {n}: the record states charge {atom.charge}, radical '
                                     f'{atom.is_radical} and isotope {atom.isotope} here, and {label} '
                                     f'states its own; the label is kept', REFUSED))
                continue
            sites.append((n, label, row, atom.r_index, molecule.xy_of(n)))

        if not sites:
            return False

        with molecule.edit():
            for n, _, row, r_index, xy in sites:
                fragment = row.fragment
                anchor = fragment.atom(row.attachment)
                if r_index:                            # an R index is only settable on element 0
                    molecule.set_r_index(n, 0)
                molecule.set_element(n, anchor.element)
                molecule.set_charge(n, anchor.charge)
                molecule.set_radical(n, anchor.is_radical)
                molecule.set_isotope(n, anchor.isotope)
                molecule.set_hydrogens(n, H_UNKNOWN if anchor.implicit_h is None else anchor.implicit_h)

                grafted = {row.attachment: n}
                for a in fragment.atoms():
                    if a.n == row.marker or a.n == row.attachment:
                        continue
                    grafted[a.n] = molecule.add_atom(a.element, charge=a.charge, isotope=a.isotope,
                                                     radical=a.is_radical, implicit_h=a.implicit_h)
                    if xy is not None:
                        molecule.set_xy(grafted[a.n], xy[0], xy[1])
                for bond in fragment.bonds():
                    if bond.n == row.marker or bond.m == row.marker:
                        continue
                    molecule.add_bond(grafted[bond.n], grafted[bond.m], bond.order)

        molecule.set_aliases({n: text for n, text in aliases.items()
                              if n not in {n for n, _, _, _, _ in sites}})
        for n, label, row, _, xy in sites:
            drawn = ', and the grafted atoms share its coordinates, so the record needs a 2D clean' \
                if xy is not None else ''
            log.append(LogRecord(row.id, (n,), f'atom {n}: the label {label} named a contracted group '
                                              f'and was expanded to {row.smiles}{drawn}', REPAIRED))
        return True
