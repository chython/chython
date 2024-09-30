# -*- coding: utf-8 -*-
#
#  Copyright 2023 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from ..containers import MoleculeContainer, ReactionContainer
from ..containers.bonds import Bond
from ..exceptions import AtomNotFound
from ..periodictable import Element


def create_molecule(data, *, skip_calc_implicit=False, ignore_bad_isotopes=False, _cls=MoleculeContainer):
    g = object.__new__(_cls)
    pm = {}
    atoms = {}
    plane = {}
    charges = {}
    radicals = {}
    bonds = {}
    mapping = data['mapping']
    for n, atom in enumerate(data['atoms']):
        n = mapping[n]
        e = Element.from_symbol(atom['element'])
        try:
            atoms[n] = e(atom['isotope'])
        except ValueError:
            if not ignore_bad_isotopes:
                raise
            atoms[n] = e()  # reset isotope mark on errors.
        bonds[n] = {}
        if (charge := atom['charge']) > 4 or charge < -4:
            raise ValueError('formal charge should be in range [-4, 4]')
        charges[n] = charge
        radicals[n] = atom['is_radical']
        plane[n] = (atom['x'], atom['y'])
        pm[n] = atom['mapping']
    for n, m, b in data['bonds']:
        n, m = mapping[n], mapping[m]
        if n == m:
            raise ValueError('atom loops impossible')
        if n not in bonds or m not in bonds:
            raise AtomNotFound('atoms not found')
        if n in bonds[m]:
            raise ValueError('atoms already bonded')
        bonds[n][m] = bonds[m][n] = Bond(b)
    if any(a['z'] for a in data['atoms']):
        conformers = [{mapping[n]: (a['x'], a['y'], a['z']) for n, a in enumerate(data['atoms'])}]
    else:
        conformers = []

    if data['log']:  # store log to the meta
        if data['meta'] is None:
            data['meta'] = {}
        data['meta']['chython_parsing_log'] = data['log']

    g.__setstate__({'atoms': atoms, 'bonds': bonds, 'meta': data['meta'], 'plane': plane, 'parsed_mapping': pm,
                    'charges': charges, 'radicals': radicals, 'name': data['title'], 'conformers': conformers,
                    'atoms_stereo': {}, 'allenes_stereo': {}, 'cis_trans_stereo': {}, 'hydrogens': {}})
    if not skip_calc_implicit:
        for n in atoms:
            g._calc_implicit(n)
    return g


def create_reaction(data, *, ignore=True, skip_calc_implicit=False, ignore_bad_isotopes=False,
                    _r_cls=ReactionContainer, _m_cls=MoleculeContainer):
    rc, pr, rg = [], [], []
    for ms, pms, gr in ((rc, data['reactants'], 'reactant'),
                        (pr, data['products'], 'products'),
                        (rg, data['reagents'], 'reagent')):
        tdl = []
        for n, m in enumerate(pms):
            try:
                ms.append(create_molecule(m, skip_calc_implicit=skip_calc_implicit,
                                          ignore_bad_isotopes=ignore_bad_isotopes, _cls=_m_cls))
            except ValueError as e:
                if not ignore:
                    raise
                data['log'].append(f'ignored {gr} molecule {n} with {e}')
                tdl.append(n)
        if tdl:  # ad-hoc for later postprocessing
            for n in reversed(tdl):
                del pms[n]

    if data['log']:  # store log to the meta
        if data['meta'] is None:
            data['meta'] = {}
        data['meta']['chython_parsing_log'] = data['log']
    return _r_cls(rc, pr, rg, meta=data['meta'], name=data['title'])


__all__ = ['create_molecule']
