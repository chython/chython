# -*- coding: utf-8 -*-
#
#  Copyright 2014-2022 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from ..containers import MoleculeContainer, QueryContainer
from ..containers.bonds import Bond
from ..periodictable import Element


class BaseReactor:
    def __init__(self, reactants, products, delete_atoms, fix_rings, fix_tautomers):
        self.__to_delete = reactants.difference(products) if delete_atoms else ()

        # prepare atoms patch
        self.__elements = elements = {}
        self.__hydrogens = hydrogens = {}

        atoms = defaultdict(dict)
        for n, atom in products.atoms():
            atoms[n].update(charge=atom.charge, is_radical=atom.is_radical)
            if atom.atomic_number:  # replace atom
                elements[n] = Element.from_atomic_number(atom.atomic_number)(atom.isotope)
                if n not in reactants and isinstance(products, MoleculeContainer):
                    atoms[n]['xy'] = atom.xy
                    if atom.hybridization == 4 and atom.implicit_hydrogens is not None:
                        hydrogens[n] = atom.implicit_hydrogens  # keep it for new aromatic atoms
            elif n not in reactants:
                raise ValueError('New atom should be defined')
            else:  # use atom from reactant
                elements[n] = None

        if isinstance(products, QueryContainer):
            bonds = []
            for n, m, b in products.bonds():
                if len(b.order) > 1:
                    raise ValueError('bond list in patch not supported')
                else:
                    bonds.append((n, m, Bond(b.order[0])))
        else:
            bonds = [(n, m, b.copy()) for n, m, b in products.bonds()]

        self.__bonds = bonds
        self.__atom_attrs = dict(atoms)
        self.__products = products
        self.__fix_rings = fix_rings
        self.__fix_tautomers = fix_tautomers

    def _patcher(self, structure: MoleculeContainer, mapping):
        elements = self.__elements
        products = self.__products
        patch_hydrogens = self.__hydrogens
        patch_bonds = self.__bonds

        atoms = structure._atoms
        plane = structure._plane
        bonds = structure._bonds
        charges = structure._charges
        radicals = structure._radicals
        hydrogens = structure._hydrogens

        to_delete = {mapping[x] for x in self.__to_delete}
        if to_delete:
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

        new = structure.__class__()
        keep_hydrogens = {}
        max_atom = max(atoms)
        for n, atom in self.__atom_attrs.items():
            if n in mapping:  # add matched atoms
                m = mapping[n]
                e = elements[n]
                if e is None:
                    e = atoms[m]
                new.add_atom(e.copy(), m, xy=plane[m], **atom)
            else:  # new atoms
                max_atom += 1
                mapping[n] = new.add_atom(elements[n].copy(), max_atom, **atom)
                if n in patch_hydrogens:  # keep patch aromatic atoms hydrogens count
                    keep_hydrogens[max_atom] = patch_hydrogens[n]

        patch_atoms = set(new)  # don't move!
        for n, atom in structure.atoms():  # add unmatched atoms
            if n not in patch_atoms and n not in to_delete:
                new.add_atom(atom.copy(), n, charge=charges[n], is_radical=radicals[n], xy=plane[n])
                if structure.hybridization(n) == 4:  # keep hydrogens on unmatched aromatic atoms as is.
                    keep_hydrogens[n] = hydrogens[n]

        for n, m, bond in patch_bonds:  # add patch bonds
            new.add_bond(mapping[n], mapping[m], bond.copy())

        for n, m_bond in bonds.items():
            if n in to_delete:  # atoms for removing
                continue
            to_delete.add(n)  # reuse to_delete set for seen atoms
            for m, bond in m_bond.items():
                # ignore deleted atoms and patch atoms
                if m in to_delete or n in patch_atoms and m in patch_atoms:
                    continue
                new.add_bond(n, m, bond.copy())

        # fix hydrogens count.
        new._hydrogens.update(keep_hydrogens)
        if self.__fix_rings:
            new.kekule()
            new.thiele(fix_tautomers=self.__fix_tautomers)

        # check needs of stereo calculations
        if structure._atoms_stereo or structure._allenes_stereo or structure._cis_trans_stereo or \
                products._atoms_stereo or products._allenes_stereo or products._cis_trans_stereo:
            stereo_override = set()
            r_mapping = {m: n for n, m in mapping.items()}

            # set patch atoms stereo
            for n, s in products._atoms_stereo.items():
                m = mapping[n]
                new._atoms_stereo[m] = products._translate_tetrahedron_sign(n, [r_mapping[x] for x in
                                                                                new._stereo_tetrahedrons[m]], s)
                stereo_override.add(m)

            for n, s in products._allenes_stereo.items():
                m = mapping[n]
                t1, t2, *_ = new._stereo_allenes[m]
                new._allenes_stereo[m] = products._translate_allene_sign(n, r_mapping[t1], r_mapping[t2], s)
                stereo_override.add(m)

            for (n, m), s in products._cis_trans_stereo.items():
                nm = (mapping[n], mapping[m])
                try:
                    t1, t2, *_ = new._stereo_cis_trans[nm]
                except KeyError:
                    nm = nm[::-1]
                    t2, t1, *_ = new._stereo_cis_trans[nm]
                new._cis_trans_stereo[nm] = products._translate_cis_trans_sign(n, m, r_mapping[t1], r_mapping[t2], s)
                stereo_override.update(nm)

            # set unmatched part stereo and not overridden by patch.
            for n, s in structure._atoms_stereo.items():
                if n in stereo_override or n not in new._stereo_tetrahedrons or \
                        new._bonds[n].keys() != structure._bonds[n].keys():
                    # skip atoms with changed neighbors
                    continue
                new._atoms_stereo[n] = structure._translate_tetrahedron_sign(n, new._stereo_tetrahedrons[n], s)

            for n, s in structure._allenes_stereo.items():
                if n in stereo_override or n not in new._stereo_allenes or \
                        set(new._stereo_allenes[n]) != set(structure._stereo_allenes[n]):
                    # skip changed allenes
                    continue
                t1, t2, *_ = new._stereo_allenes[n]
                new._allenes_stereo[n] = structure._translate_allene_sign(n, t1, t2, s)

            for nm, s in structure._cis_trans_stereo.items():
                n, m = nm
                if n in stereo_override or m in stereo_override:
                    continue
                env = structure._stereo_cis_trans[nm]
                try:
                    new_env = new._stereo_cis_trans[nm]
                except KeyError:
                    nm = nm[::-1]
                    try:
                        new_env = new._stereo_cis_trans[nm]
                    except KeyError:
                        continue
                    t2, t1, *_ = new_env
                else:
                    t1, t2, *_ = new_env
                if set(env) != set(new_env):
                    continue
                new._cis_trans_stereo[nm] = structure._translate_cis_trans_sign(n, m, t1, t2, s)

            new.fix_stereo()
        return new


__all__ = ['BaseReactor']
