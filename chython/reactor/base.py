# -*- coding: utf-8 -*-
#
#  Copyright 2014-2024 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from collections import defaultdict
from itertools import product
from typing import Union
from ..containers import MoleculeContainer, QueryContainer
from ..containers.bonds import Bond
from ..periodictable import Element, ListElement, AnyElement, QueryElement, AnyMetal


class BaseReactor:
    def __init__(self, pattern, replacement, delete_atoms, fix_rings, fix_tautomers):
        if isinstance(replacement, QueryContainer):
            for n, a in replacement.atoms():
                if not isinstance(a, (AnyElement, QueryElement)):
                    raise TypeError('Unsupported query atom type')
            for *_, b in replacement.bonds():
                if len(b.order) > 1:
                    raise ValueError('Variable bond in replacement')

        self._to_delete = {n for n, a in pattern.atoms() if not a.masked} - set(replacement) if delete_atoms else ()
        self._replacement = replacement
        self._fix_rings = fix_rings
        self._fix_tautomers = fix_tautomers

    def _patcher(self, structure: MoleculeContainer, mapping):
        new = self._prepare_skeleton(structure, mapping)
        self._fix_stereo(new, structure, mapping)

        if self._fix_rings:
            new.kekule()  # keeps stereo as is
            if not new.thiele(fix_tautomers=self._fix_tautomers):  # fixes stereo if any ring aromatized
                new.fix_stereo()
        else:
            new.fix_stereo()
        yield new

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

    def _prepare_skeleton(self, structure, mapping):
        atoms = structure._atoms
        bonds = structure._bonds

        to_delete = self._get_deleted(structure, mapping)
        new = structure.__class__()
        natoms = new._atoms
        nbonds = new._bonds
        max_atom = max(atoms)
        stereo_atoms = []
        stereo_bonds = []

        for n, a in self._replacement.atoms():
            if isinstance(a, AnyElement):
                if n := mapping.get(n):
                    # keep matched atom type and isotope
                    e = atoms[n].copy(stereo=True)
                    e.charge = a.charge
                    e.is_radical = a.is_radical
                    if a.stereo is not None:  # override stereo
                        e._stereo = a.stereo
                    elif e.stereo is not None:  # keep original stereo
                        stereo_atoms.append(n)  # mark for stereo fix
                    natoms[n] = e
                    nbonds[n] = {}
                else:
                    raise ValueError("AnyElement doesn't match to pattern")
            else:  # QueryElement or Element
                a: Union[QueryElement, Element]  # typehint
                e = Element.from_atomic_number(a.atomic_number)
                e = e(a.isotope, charge=a.charge, is_radical=a.is_radical, stereo=a.stereo)
                if not (m := mapping.get(n)):  # new atom
                    m = max_atom + 1
                    max_atom += 1
                    mapping[n] = m
                    if isinstance(a, Element):
                        e._implicit_hydrogens = a.implicit_hydrogens  # keep H count from patch
                        e.x = a.x  # keep coordinates from patch
                        e.y = a.y
                    elif len(a.implicit_hydrogens) == 1:
                        e._implicit_hydrogens = a.implicit_hydrogens[0]
                    elif a.implicit_hydrogens:
                        raise ValueError('Query element in patch has more than one implicit hydrogen')
                else:  # existing atoms
                    b = atoms[m]
                    e.x = b.x  # preserve existing coordinates
                    e.y = b.y
                    if a.stereo is None and b.stereo is not None:  # keep original stereo
                        e._stereo = b.stereo
                        stereo_atoms.append(m)
                natoms[m] = e
                nbonds[m] = {}

        # preserve connectivity order
        for n, bs in self._replacement._bonds.items():
            n = mapping[n]
            for m, b in bs.items():
                m = mapping[m]
                if n in nbonds[m]:
                    nbonds[n][m] = nbonds[m][n]
                else:
                    nbonds[n][m] = b = Bond(int(b), stereo=b.stereo)
                    if b.stereo is None:
                        if not (nb := bonds.get(n)):
                            continue
                        if not (mb := nb.get(m)):
                            continue
                        if mb.stereo is None:
                            continue
                        # original structure has stereo bond
                        b._stereo = mb.stereo
                        stereo_bonds.append((n, m))

        patch_atoms = set(new)  # don't move!
        for n, a in atoms.items():  # add unmatched or masked atoms
            if n not in patch_atoms and n not in to_delete:
                natoms[n] = a.copy(hydrogens=True, stereo=True)
                nbonds[n] = {}

        for n, bs in bonds.items():
            if n in to_delete:  # atoms for removing
                continue
            for m, b in bs.items():
                # ignore deleted atoms and patch atoms
                if m in to_delete or n in patch_atoms and m in patch_atoms:
                    continue
                elif n in nbonds[m]:
                    nbonds[n][m] = nbonds[m][n]
                else:
                    nbonds[n][m] = b.copy(stereo=True)
                    if b.stereo is not None and (n in patch_atoms or m in patch_atoms):
                        stereo_bonds.append((n, m))

        for n, a in new.atoms():
            if a.implicit_hydrogens is None:
                new.calc_implicit(n)
        new.calc_labels()
        return new

    def _fix_stereo(self, new, structure, mapping):
        products = self.__products
        stereo_override = set()
        r_mapping = {m: n for n, m in mapping.items()}

        # set patch atoms stereo
        for n, s in products._atoms_stereo.items():
            m = mapping[n]
            new._atoms_stereo[m] = products._translate_tetrahedron_sign(n, [r_mapping[x] for x in
                                                                            new.stereogenic_tetrahedrons[m]], s)
            stereo_override.add(m)

        for n, s in products._allenes_stereo.items():
            m = mapping[n]
            t1, t2, *_ = new.stereogenic_allenes[m]
            new._allenes_stereo[m] = products._translate_allene_sign(n, r_mapping[t1], r_mapping[t2], s)
            stereo_override.add(m)

        for (n, m), s in products._cis_trans_stereo.items():
            nm = (mapping[n], mapping[m])
            try:
                t1, t2, *_ = new.stereogenic_cis_trans[nm]
            except KeyError:
                nm = nm[::-1]
                t2, t1, *_ = new.stereogenic_cis_trans[nm]
            new._cis_trans_stereo[nm] = products._translate_cis_trans_sign(n, m, r_mapping[t1], r_mapping[t2], s)
            stereo_override.update(nm)

        # set unmatched part stereo and not overridden by patch.
        for n, s in structure._atoms_stereo.items():
            if n in stereo_override or n not in new.stereogenic_tetrahedrons or \
                    new._bonds[n].keys() != structure._bonds[n].keys():
                # skip atoms with changed neighbors
                continue
            new._atoms_stereo[n] = structure._translate_tetrahedron_sign(n, new.stereogenic_tetrahedrons[n], s)

        for n, s in structure._allenes_stereo.items():
            if n in stereo_override or n not in new.stereogenic_allenes or \
                    set(new.stereogenic_allenes[n]) != set(structure.stereogenic_allenes[n]):
                # skip changed allenes
                continue
            t1, t2, *_ = new.stereogenic_allenes[n]
            new._allenes_stereo[n] = structure._translate_allene_sign(n, t1, t2, s)

        for nm, s in structure._cis_trans_stereo.items():
            n, m = nm
            if n in stereo_override or m in stereo_override:
                continue
            env = structure.stereogenic_cis_trans[nm]
            try:
                new_env = new.stereogenic_cis_trans[nm]
            except KeyError:
                nm = nm[::-1]
                try:
                    new_env = new.stereogenic_cis_trans[nm]
                except KeyError:
                    continue
                t2, t1, *_ = new_env
            else:
                t1, t2, *_ = new_env
            if set(env) != set(new_env):
                continue
            new._cis_trans_stereo[nm] = structure._translate_cis_trans_sign(n, m, t1, t2, s)


__all__ = ['BaseReactor']
