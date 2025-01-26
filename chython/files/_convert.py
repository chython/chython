# -*- coding: utf-8 -*-
#
#  Copyright 2023, 2024 Ramil Nugmanov <nougmanoff@protonmail.com>
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


# atomic number constants
B = 5
C = 6
N = 7
P = 15


def create_molecule(data, *, ignore_bad_isotopes=False, skip_calc_implicit=False,
                    keep_implicit=False, keep_radicals=True, ignore_aromatic_radicals=True, ignore=True,
                    ignore_carbon_radicals=False, _cls=MoleculeContainer):
    g = _cls()
    g._name = data.get('title')
    atoms = g._atoms
    bonds = g._bonds
    mapping = data['mapping']

    if any(a.get('z') for a in data['atoms']):
        # store conformer
        g._conformers = [{n: (a['x'], a['y'], a['z']) for n, a in zip(mapping, data['atoms'])}]

    for n, atom in zip(mapping, data['atoms']):
        e = Element.from_symbol(atom.pop('element'))
        atom.pop('z', None)  # clean up MDL
        try:
            atoms[n] = e(**atom)
        except (ValueError, TypeError):
            if not ignore_bad_isotopes:
                raise
            del atom['isotope']  # reset isotope mark on errors and try again.
            atoms[n] = e(**atom)
        bonds[n] = {}

    for n, m, b in data['bonds']:
        n, m = mapping[n], mapping[m]
        if n == m:
            raise ValueError('atom loops impossible')
        if n not in bonds or m not in bonds:
            raise AtomNotFound('atoms not found')
        if n in bonds[m]:
            raise ValueError('atoms already bonded')
        bonds[n][m] = bonds[m][n] = Bond(b)

    g.calc_labels()  # set all labels except rings

    if data.get('log'):  # store log to the meta
        if data.get('meta') is None:
            data['meta'] = {}
        data['meta']['chython_parsing_log'] = data['log']
    g._meta = data.get('meta') or None

    if skip_calc_implicit:  # don't calc Hs. e.g. INCHI
        return g

    implicit_mismatch = {}
    radicalized = []
    # precalculate Hs
    for n, a in atoms.items():
        if a.implicit_hydrogens is None:
            # let's try to calculate. in case of errors just keep as is. radicals in smiles should be in [brackets],
            # thus has implicit Hs value
            g.calc_implicit(n)
        elif keep_implicit:
            # keep given Hs count as is
            continue
        else:  # recheck given Hs count
            h = a.implicit_hydrogens  # parsed Hs
            g.calc_implicit(n)  #  recalculate
            if a.implicit_hydrogens is None:  # atom has invalid valence or aromatic ring.
                if a.hybridization == 4:
                    # this is aromatic ring. just restore given H count.
                    a._implicit_hydrogens = h
                    # rare H0 case
                    if (not keep_radicals and not ignore_aromatic_radicals
                        and not h and not a.charge and not a.is_radical and a in (B, C, N, P)
                        and sum(b != 8 for b in bonds[n].values()) == 2):
                        # c[c]c - aromatic B,C,N,P radical
                        a._is_radical = True
                        radicalized.append(n)
                elif not keep_radicals and not a.is_radical:  # CXSMILES radical not set.
                    # SMILES doesn't code radicals. so, let's try to guess.
                    a._is_radical = True
                    if g.check_implicit(n, h):  # radical form is valid
                        radicalized.append(n)
                        a._implicit_hydrogens = h
                    elif ignore:  # radical state also has errors.
                        a._is_radical = False  # reset radical state
                        implicit_mismatch[n] = h
                        if data.get('log') is None:
                            data['log'] = []
                        data['log'].append(f'implicit hydrogen count ({h}) mismatch with calculated on atom {n}')
                    else:
                        raise ValueError(f'implicit hydrogen count ({h}) mismatch with calculated on atom {n}')
            elif h != a.implicit_hydrogens:  # H count mismatch.
                if a.hybridization == 4:
                    if (not keep_radicals
                        and not h and not a.charge and not a.is_radical and a in (B, C, N, P)
                        and sum(b != 8 for b in bonds[n].values()) == 2):
                        # c[c]c - aromatic B,C,N,P radical
                        a._implicit_hydrogens = 0
                        a._is_radical = True
                        radicalized.append(n)
                    elif ignore:
                        implicit_mismatch[n] = h
                        if data.get('log') is None:
                            data['log'] = []
                        data['log'].append(f'implicit hydrogen count ({h}) mismatch with calculated on atom {n}')
                    else:
                        raise ValueError(f'implicit hydrogen count ({h}) mismatch with calculated on atom {n}')
                elif g.check_implicit(n, h):  # set another possible implicit state. probably Al, P
                    a._implicit_hydrogens = h
                elif not keep_radicals and not a.is_radical:  # CXSMILES radical is not set. try radical form
                    a._is_radical = True
                    if g.check_implicit(n, h):
                        a._implicit_hydrogens = h
                        radicalized.append(n)
                    # radical state also has errors.
                    elif ignore:
                        a._is_radical = False  # reset radical state
                        implicit_mismatch[n] = h
                        if data.get('log') is None:
                            data['log'] = []
                        data['log'].append(f'implicit hydrogen count ({h}) mismatch with calculated on atom {n}')
                    else:
                        raise ValueError(f'implicit hydrogen count ({h}) mismatch with calculated on atom {n}')
                elif ignore:  # just ignore it
                    implicit_mismatch[n] = h
                    if data.get('log') is None:
                        data['log'] = []
                    data['log'].append(f'implicit hydrogen count ({h}) mismatch with calculated on atom {n}')
                else:
                    raise ValueError(f'implicit hydrogen count ({h}) mismatch with calculated on atom {n}')

    if ignore_carbon_radicals:
        for n in radicalized:
            if (a := atoms[n]) == C:
                a._is_radical = False
                a._implicit_hydrogens += 1
                if data.get('log') is None:
                    data['log'] = []
                data['log'].append(f'carbon radical {n} replaced with implicit hydrogen')
    elif radicalized:
        g.meta['chython_radicalized_atoms'] = radicalized
    if data.get('log') and 'chython_parsing_log' not in g.meta:
        g.meta['chython_parsing_log'] = data['log']
    if implicit_mismatch:
        g.meta['chython_implicit_mismatch'] = implicit_mismatch
    return g


def create_reaction(data, *, ignore=True, skip_calc_implicit=False, ignore_bad_isotopes=False,
                    keep_implicit=False, keep_radicals=True, ignore_aromatic_radicals=True,
                    ignore_carbon_radicals=False, _r_cls=ReactionContainer, _m_cls=MoleculeContainer):
    rc, pr, rg = [], [], []
    for ms, pms, gr in ((rc, data['reactants'], 'reactant'),
                        (pr, data['products'], 'products'),
                        (rg, data['reagents'], 'reagent')):
        tdl = []
        for n, m in enumerate(pms):
            try:
                ms.append(create_molecule(m, skip_calc_implicit=skip_calc_implicit,
                                          ignore_bad_isotopes=ignore_bad_isotopes, keep_implicit=keep_implicit,
                                          keep_radicals=keep_radicals,
                                          ignore_aromatic_radicals=ignore_aromatic_radicals, ignore=ignore,
                                          ignore_carbon_radicals=ignore_carbon_radicals, _cls=_m_cls))
            except ValueError as e:
                if not ignore:
                    raise
                if data.get('log') is None:
                    data['log'] = []
                data['log'].append(f'ignored {gr} molecule {n} with {e}')
                tdl.append(n)
        if tdl:  # ad-hoc for later postprocessing
            for n in reversed(tdl):
                del pms[n]

    if data.get('log'):  # store log to the meta
        if data.get('meta') is None:
            data['meta'] = {}
        data['meta']['chython_parsing_log'] = data['log']
    return _r_cls(rc, pr, rg, meta=data.get('meta') or None, name=data.get('title'))


__all__ = ['create_molecule']
