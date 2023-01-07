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
from ...exceptions import NotChiral, IsChiral, ValenceError


def postprocess_molecule(molecule, data, *, ignore=True, ignore_stereo=False, calc_cis_trans=False):
    mapping = data['mapping']
    hydrogens = molecule._hydrogens
    for n, h in data['hydrogens'].items():
        n = mapping[n]
        hc = hydrogens[n]
        if hc is None:  # aromatic rings or valence errors. just store given H count.
            hydrogens[n] = h
        elif hc != h:
            if ignore:
                hydrogens[n] = h  # set parsed hydrogens count
                if 'chython_parsing_log' not in molecule.meta:
                    molecule.meta['chython_parsing_log'] = []
                molecule.meta['chython_parsing_log'].append(f'implicit hydrogen count ({h}) mismatch with calculated '
                                                            f'({hc}) on atom {n}. calculated count replaced.')
            else:
                raise ValueError(f'implicit hydrogen count ({h}) mismatch with calculated ({hc}) on atom {n}.')

    if ignore_stereo:
        return

    if calc_cis_trans:
        molecule.calculate_cis_trans_from_2d()

    stereo = [(mapping[n], mapping[m], s) for n, m, s in data['stereo']]
    while stereo:
        fail_stereo = []
        old_stereo = len(stereo)
        for n, m, s in stereo:
            try:
                molecule.add_wedge(n, m, s, clean_cache=False)
            except NotChiral:
                fail_stereo.append((n, m, s))
            except IsChiral:
                pass
            except ValenceError:
                if 'chython_parsing_log' not in molecule.meta:
                    molecule.meta['chython_parsing_log'] = []
                molecule.meta['chython_parsing_log'].append('structure has errors, stereo data skipped')
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


__all__ = ['postprocess_molecule']
