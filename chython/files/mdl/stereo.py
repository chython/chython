# -*- coding: utf-8 -*-
#
#  Copyright 2020-2026 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from ...exceptions import NotChiral, IsChiral, ValenceError


def postprocess_molecule(molecule, data, *, ignore_stereo=False, calc_cis_trans=False):
    if ignore_stereo:
        return
    mapping = data['mapping']
    log = []

    if calc_cis_trans:
        molecule.calculate_cis_trans_from_2d(clean_cache=False)

    # a stereo centre only becomes chiral once its neighbours are resolved, so both signals need
    # iterative resolution (e.g. CC(O)C(O)C(C)O needs several passes). merge them into one loop as
    # (method, *args) work items - the pattern the SMILES/InChI parsers use. coordinate-free atom
    # parity is listed first so it wins on shared centres: it sets the sign, then add_wedge raises
    # IsChiral and is dropped.
    stereo = []
    if data.get('stereo_atoms'):
        # neighbours are numbered by ascending atom-block position, not chython atom number.
        file_of = {cn: fi for fi, cn in enumerate(mapping)}
        bonds = molecule._bonds
        for i, mark in data['stereo_atoms']:
            n = mapping[i]
            env = sorted(bonds[n], key=file_of.__getitem__)
            stereo.append((molecule.add_atom_stereo, n, env, mark))
    stereo.extend((molecule.add_wedge, mapping[n], mapping[m], s) for n, m, s in data['stereo'])

    while stereo:
        fail_stereo = []
        old_stereo = len(stereo)
        for f, *args in stereo:
            try:
                f(*args, clean_cache=False)
            except NotChiral:
                fail_stereo.append((f, *args))
            except IsChiral:
                pass
            except ValenceError:
                log.append('structure has errors, stereo data skipped')
                molecule.flush_cache()
                break
        else:
            stereo = fail_stereo
            if len(stereo) == old_stereo:
                break
            molecule.flush_stereo_cache()
            if calc_cis_trans:
                molecule.calculate_cis_trans_from_2d(clean_cache=False)
            continue
        break

    if log:
        if 'chython_parsing_log' not in molecule.meta:
            molecule.meta['chython_parsing_log'] = log
        else:
            molecule.meta['chython_parsing_log'].extend(log)


__all__ = ['postprocess_molecule']
