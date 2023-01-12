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


def postprocess_molecule(molecule, data, *, ignore=True, ignore_stereo=False, calc_cis_trans=False,
                         keep_implicit=False):
    mapping = data['mapping']
    hydrogens = molecule._hydrogens
    hyb = molecule.hybridization

    implicit_mismatch = {}
    if 'chython_parsing_log' in molecule.meta:
        log = molecule.meta['chython_parsing_log']
    else:
        log = []

    for n, h in data['hydrogens'].items():
        n = mapping[n]
        if keep_implicit:  # override any calculated hydrogens count.
            hydrogens[n] = h
        if (hc := hydrogens[n]) is None:  # aromatic rings or valence errors
            if hyb(n) == 4:  # this is aromatic rings. just store given H count.
                hydrogens[n] = h
        elif hc != h:
            if hyb(n) == 4:
                if ignore:
                    implicit_mismatch[n] = h
                    log.append(f'implicit hydrogen count ({h}) mismatch with calculated on atom {n}')
                else:
                    raise ValueError(f'implicit hydrogen count ({h}) mismatch with calculated on atom {n}')
            elif molecule._check_implicit(n, h):  # set another possible implicit state. probably Al, P
                hydrogens[n] = h
            elif ignore:  # just ignore it
                implicit_mismatch[n] = h
                log.append(f'implicit hydrogen count ({h}) mismatch with calculated on atom {n}')
            else:
                raise ValueError(f'implicit hydrogen count ({h}) mismatch with calculated on atom {n}')

    if implicit_mismatch:
        molecule.meta['chython_implicit_mismatch'] = implicit_mismatch
    if log and 'chython_parsing_log' not in molecule.meta:
        molecule.meta['chython_parsing_log'] = log
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

    if log and 'chython_parsing_log' not in molecule.meta:
        molecule.meta['chython_parsing_log'] = log


__all__ = ['postprocess_molecule']
