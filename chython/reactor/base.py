# -*- coding: utf-8 -*-
#
#  Copyright 2014-2025 Ramil Nugmanov <nougmanoff@protonmail.com>
#  Copyright 2019 Adelia Fatykhova <adelik21979@gmail.com>
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
from typing import Union
from ..containers import MoleculeContainer, QueryContainer
from ..containers.bonds import Bond
from ..periodictable import Element, AnyElement, QueryElement


class BaseReactor:
    def __init__(self, pattern, replacement, delete_atoms, fix_rings, fix_tautomers, fix_broken_pyrroles):
        if isinstance(replacement, QueryContainer):
            for n, a in replacement.atoms():
                if not isinstance(a, (AnyElement, QueryElement)):
                    raise TypeError('Unsupported query atom type')
                elif len(a.implicit_hydrogens) > 1:
                    raise ValueError('Query element in patch has more than one implicit hydrogen clause')
            for *_, b in replacement.bonds():
                if len(b.order) > 1:
                    raise ValueError('Variable bond in replacement')

        self._to_delete = {n for n, a in pattern.atoms() if not a.masked} - set(replacement) if delete_atoms else ()
        self._replacement = replacement
        self._fix_rings = fix_rings
        self._fix_tautomers = fix_tautomers
        self._fix_broken_pyrroles = fix_broken_pyrroles

    def _get_deleted(self, structure, mapping):
        if not self._to_delete:
            return set()

        bonds = structure._bonds
        to_delete = {mapping[x] for x in self._to_delete}
        # if deleted atoms have another path to remain fragment, the path is preserved
        remain = set(mapping.values()).difference(to_delete)
        delete, global_seen = set(), set()
        for x in to_delete:
            for n in bonds[x]:
                if n in global_seen or n in remain:
                    continue
                seen = {n}
                global_seen.add(n)
                stack = [x for x in bonds[n] if x not in global_seen]
                while stack:
                    current = stack.pop()
                    if current in remain:
                        break
                    if current in to_delete:
                        continue
                    seen.add(current)
                    global_seen.add(current)
                    stack.extend([x for x in bonds[current] if x not in global_seen])
                else:
                    delete.update(seen)

        to_delete.update(delete)
        return to_delete

    def _patcher(self, structure: MoleculeContainer, mapping):
        satoms = structure._atoms
        sbonds = structure._bonds

        to_delete = self._get_deleted(structure, mapping)
        new = structure.__class__()
        natoms = new._atoms
        nbonds = new._bonds
        max_atom = max(satoms)
        stereo_atoms = []
        stereo_bonds = []

        # let's preserve connectivity order from replacement to keep stereo signs as is.
        # stereo labels from original structure will be recalculated after full molecule construction.
        for n, ra in self._replacement.atoms():
            if isinstance(ra, AnyElement):
                if m := mapping.get(n):
                    # keep matched atom type and isotope
                    sa = satoms[m]
                    a = sa.copy()
                    a.charge = ra.charge
                    a.is_radical = ra.is_radical
                    if ra.stereo is not None:  # override stereo
                        a._stereo = ra.stereo
                    elif sa.stereo is not None:  # keep original stereo
                        stereo_atoms.append(m)  # mark for stereo fix
                else:
                    raise ValueError("AnyElement doesn't match to pattern")
            else:  # QueryElement or Element
                ra: Union[QueryElement, Element]  # typehint
                e = Element.from_atomic_number(ra.atomic_number)
                a = e(ra.isotope, charge=ra.charge, is_radical=ra.is_radical)
                if not (m := mapping.get(n)):  # new atom
                    m = max_atom + 1
                    mapping[n] = max_atom = m
                    a._stereo = ra.stereo  # keep stereo from patch for new atoms
                    if isinstance(ra, Element):
                        a._implicit_hydrogens = ra.implicit_hydrogens  # keep H count from patch
                        a.xy = ra.xy  # keep coordinates from patch
                    elif ra.implicit_hydrogens:  # keep H count from patch
                        a._implicit_hydrogens = ra.implicit_hydrogens[0]
                else:  # existing atoms
                    sa = satoms[m]
                    a.xy = sa.xy  # preserve existing coordinates
                    if ra.stereo is not None:
                        a._stereo = ra.stereo
                    elif sa.stereo is not None:  # keep original stereo
                        stereo_atoms.append(m)
            natoms[m] = a
            nbonds[m] = {}

        # preserve connectivity order
        for n, bs in self._replacement._bonds.items():
            n = mapping[n]
            for m, rb in bs.items():
                m = mapping[m]
                if n in nbonds[m]:  # back-link
                    nbonds[n][m] = nbonds[m][n]
                else:
                    nbonds[n][m] = b = Bond(int(rb))
                    if rb.stereo is not None:  # override stereo
                        b._stereo = rb.stereo
                    # check bond exists in source and has stereo label and the same order
                    elif (sbn := sbonds.get(n)) is None or (sb := sbn.get(m)) is None or sb.stereo is None or sb != b:
                        continue
                    else:  # original structure has stereo bond
                        stereo_bonds.append((n, m))

        patched_atoms = set(new)
        for n, sa in satoms.items():  # add unmatched or masked atoms
            if n not in patched_atoms and n not in to_delete:
                natoms[n] = a = sa.copy(hydrogens=True)
                nbonds[n] = {}
                if sa.stereo is not None:
                    # in case of allenes label can disappear/change, thus, requires recalculation
                    # for tetrahedrons label can be stored as is
                    if n in structure.stereogenic_tetrahedrons:
                        a._stereo = sa.stereo
                    else:
                        stereo_atoms.append(n)

        for n, bs in sbonds.items():  # preserve connectivity order for keeping stereo labels as is
            if n in to_delete:  # atoms for removing
                continue
            for m, b in bs.items():
                # ignore deleted atoms and patch atoms
                if m in to_delete or n in patched_atoms and m in patched_atoms:
                    continue
                elif n in nbonds[m]:  # back-link
                    nbonds[n][m] = nbonds[m][n]
                else:
                    nbonds[n][m] = b.copy()
                    if b.stereo is not None:
                        # stereo label should be recalculated
                        stereo_bonds.append((n, m))

        for n, a in new.atoms():
            if a.implicit_hydrogens is None:
                new.calc_implicit(n)
        new.calc_labels()

        # translate stereo sign from old order to new order
        for n in stereo_atoms:
            if n in new.stereogenic_tetrahedrons:
                if sbonds[n].keys() == nbonds[n].keys():
                    # flush stereo from reaction center. should be explicitly set in replacement.
                    s = new._translate_tetrahedron_sign(n, structure.stereogenic_tetrahedrons[n], satoms[n].stereo)
                    natoms[n]._stereo = s
            elif n in new.stereogenic_allenes:
                if set(new.stereogenic_allenes[n]) == set(structure.stereogenic_allenes[n]):
                    # flush stereo for changed allene substituents
                    s = new._translate_allene_sign(n, *structure.stereogenic_allenes[n][:2], satoms[n].stereo)
                    natoms[n]._stereo = s
            # else: ignore label

        for n, m in stereo_bonds:
            # check if bond is center of cumulene
            if (n12 := new._stereo_cis_trans_terminals.get(n, True)) != new._stereo_cis_trans_terminals.get(m, False):
                continue
            s12 = structure._stereo_cis_trans_terminals[n]
            # check if cumulene terminals are the same
            if set(n12) != set(s12):
                continue
            if set(new.stereogenic_cis_trans[n12]) == set(env := structure.stereogenic_cis_trans[s12]):
                # connected to cumulenes atoms should be the same
                s = new._translate_cis_trans_sign(*n12, *env[:2], sbonds[n][m].stereo)
                nbonds[n][m]._stereo = s
            # else: ignore label

        if self._fix_rings:
            new.kekule(ignore_pyrrole_hydrogen=self._fix_broken_pyrroles)  # keeps stereo as is
            if not new.thiele(fix_tautomers=self._fix_tautomers):  # fixes stereo if any ring aromatized
                new.fix_stereo()
        else:
            new.fix_stereo()
        return new


__all__ = ['BaseReactor']
