# -*- coding: utf-8 -*-
#
#  Copyright 2022-2025 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from itertools import chain
from re import compile
from typing import List, Union, Optional
from .parser import parser
from .tokenize import smiles_tokenize
from .._convert import create_molecule, create_reaction
from .._mapping import postprocess_parsed_molecule, postprocess_parsed_reaction
from ...containers import MoleculeContainer, ReactionContainer
from ...exceptions import IsChiral, NotChiral, ValenceError


cx_fragments = compile(r'f:(?:[0-9]+(?:\.[0-9]+)+)(?:,(?:[0-9]+(?:\.[0-9]+)+))*')
cx_radicals = compile(r'\^[1-7]:(\d+(?:,\d+)*)')
cx_stereo_rel = compile(r'([o&])(\d+):(\d+(?:,\d+)*)')
cx_stereo_abs = compile(r'a:(\d+(?:,\d+)*)')


def smiles(data, /, *, ignore: bool = True, remap: bool = False, ignore_stereo: bool = False,
           ignore_bad_isotopes: bool = False, keep_implicit: bool = False, ignore_carbon_radicals: bool = False,
           ignore_aromatic_radicals: bool = True,
           _r_cls=ReactionContainer, _m_cls=MoleculeContainer) -> Union[MoleculeContainer, ReactionContainer]:
    """
    SMILES string parser

    :param ignore: Skip some checks of data or try to fix some errors.
    :param remap: Remap atom numbers started from one.
    :param ignore_stereo: Ignore stereo data.
    :param keep_implicit: keep given in smiles implicit hydrogen count, otherwise ignore on valence error.
    :param ignore_bad_isotopes: reset invalid isotope mark to non-isotopic.
    :param ignore_carbon_radicals: fill carbon radicals with hydrogen (X[C](X)X case).
    :param ignore_aromatic_radicals: don't treat aromatic tokens like c[c]c as radicals.
    """
    if not isinstance(data, str):
        raise TypeError('Must be a SMILES string')
    elif not data:
        raise ValueError('Empty string')

    contract: Optional[List[List[int]]] = None
    radicals = []
    extended_stereo = {}
    log = []
    smi, *data = data.split()
    if data and (cxs := data[0]).startswith(('|', '{')) and cxs.endswith(('|', '}')):
        if cxs.startswith('|'):  # cx smiles
            fr = cx_fragments.search(cxs)
            if fr is not None:
                contract = [sorted(int(x) for x in x.split('.')) for x in fr.group()[2:].split(',')]
                if len({x for x in contract for x in x}) < len([x for x in contract for x in x]):
                    log.append(f'collisions in cxsmiles fragments description: {cxs}')
                    contract = None

            radicals = [int(x) for x in cx_radicals.findall(cxs) for x in x.split(',')]
            if radicals and len(set(radicals)) != len(radicals):
                log.append(f'collisions in cxsmiles radicals description: {cxs}')
                radicals = []

        # cx smiles and jnj smiles
        for st, sg, sa in cx_stereo_rel.findall(cxs):
            sg = int(sg)
            if st == 'o':
                sg = -sg
            for x in sa.split(','):
                x = int(x)
                if x in extended_stereo:
                    log.append(f'enhanced stereo description collision: {cxs}')
                    extended_stereo = {}  # ignore bad notation
                    break
                extended_stereo[x] = sg
            else:
                continue
            break

    if '>' in smi:
        record = {'reactants': [], 'reagents': [], 'products': [], 'log': log}
        try:
            reactants, reagents, products = smi.split('>')
        except ValueError as e:
            raise ValueError('invalid reaction smiles') from e

        mol_count = 0
        for k, d in zip(('reactants', 'products', 'reagents'), (reactants, products, reagents)):
            if not d:
                continue
            for x in d.split('.'):
                if not x:
                    if ignore:
                        log.append('two dots in line ignored')
                    else:
                        raise ValueError('invalid reaction smiles. two dots in line')
                else:
                    record[k].append(x)
                    mol_count += 1

        if contract:
            if max(x for x in contract for x in x) >= mol_count:
                log.append(f'skipped invalid contract data: {contract}')
            lr = len(record['reactants'])
            lp = len(record['products'])
            reactants = set(range(lr))
            reagents = set(range(lr, mol_count - lp))
            products = set(range(mol_count - lp, mol_count))
            new_molecules: List[Optional[str]] = [None] * mol_count
            for c in contract:
                if reactants.issuperset(c):
                    new_molecules[c[0]] = '.'.join(record['reactants'][x] for x in c)
                    reactants.difference_update(c)
                elif products.issuperset(c):
                    new_molecules[c[0]] = '.'.join(record['products'][x - mol_count] for x in c)
                    products.difference_update(c)
                elif reagents.issuperset(c):
                    new_molecules[c[0]] = '.'.join(record['reagents'][x - lr] for x in c)
                    reagents.difference_update(c)
                else:
                    log.append(f'impossible to contract different parts of reaction: {contract}')
            for x in reactants:
                new_molecules[x] = record['reactants'][x]
            for x in products:
                new_molecules[x] = record['products'][x - mol_count]
            for x in reagents:
                new_molecules[x] = record['reagents'][x - lr]

            record['reactants'] = [x for x in new_molecules[:lr] if x is not None]
            record['products'] = [x for x in new_molecules[-lp:] if x is not None]
            record['reagents'] = [x for x in new_molecules[lr: -lp] if x is not None]

        for k in ('reactants', 'products', 'reagents'):
            tmp = []
            for x in record[k]:
                tmp.append(parser(smiles_tokenize(x), not ignore))
            record[k] = tmp

        if radicals or extended_stereo:
            atom_map = dict(enumerate(a for m in chain(record['reactants'], record['reagents'], record['products'])
                                      for a in m['atoms']))
            for x in radicals:
                atom_map[x]['is_radical'] = True
            for x, y in extended_stereo.items():
                atom_map[x]['extended_stereo'] = y

        postprocess_parsed_reaction(record, remap=remap, ignore=ignore)
        rxn = create_reaction(record, ignore_bad_isotopes=ignore_bad_isotopes, keep_radicals=False,
                              ignore_carbon_radicals=ignore_carbon_radicals, keep_implicit=keep_implicit,
                              ignore_aromatic_radicals=ignore_aromatic_radicals, ignore=ignore,
                              _r_cls=_r_cls, _m_cls=_m_cls)
        for mol, tmp in zip(rxn.molecules(), chain(record['reactants'], record['reagents'], record['products'])):
            postprocess_molecule(mol, tmp, ignore_stereo=ignore_stereo)
        return rxn
    else:
        record = parser(smiles_tokenize(smi), not ignore)
        for x in radicals:
            record['atoms'][x]['is_radical'] = True
        for x, y in extended_stereo.items():
            record['atoms'][x]['extended_stereo'] = y
        record['log'].extend(log)

        postprocess_parsed_molecule(record, remap=remap, ignore=ignore)
        mol = create_molecule(record, ignore_bad_isotopes=ignore_bad_isotopes, keep_radicals=False,
                              ignore_carbon_radicals=ignore_carbon_radicals, keep_implicit=keep_implicit,
                              ignore_aromatic_radicals=ignore_aromatic_radicals, ignore=ignore,
                              _cls=_m_cls)
        postprocess_molecule(mol, record, ignore_stereo=ignore_stereo)
        return mol


def postprocess_molecule(molecule, data, *, ignore_stereo=False):
    mapping = data['mapping']

    if ignore_stereo:
        return
    elif not data['stereo_atoms'] and not data['stereo_bonds']:
        return

    atoms = molecule._atoms
    st = molecule.stereogenic_tetrahedrons
    sa = molecule.stereogenic_allenes
    sat = molecule._stereo_allenes_terminals
    ctc = molecule._stereo_cis_trans_counterpart

    order = {mapping[n]: [mapping[m] for m in ms] for n, ms in data['order'].items()}

    log = []
    stereo = []
    for i, s in data['stereo_atoms'].items():
        n = mapping[i]
        if not i and atoms[n].implicit_hydrogens:  # first atom in smiles has reversed chiral mark
            s = not s

        if n in st:
            stereo.append((molecule.add_atom_stereo, n, order[n], s))
        elif n in sa:
            t1, t2 = sat[n]
            env = sa[n]
            n1 = next(x for x in order[t1] if x in env)
            n2 = next(x for x in order[t2] if x in env)
            stereo.append((molecule.add_atom_stereo, n, (n1, n2), s))
        else:
            log.append(f'non chiral atom {n} has stereo label in smiles')

    stereo_bonds = {mapping[n]: {mapping[m]: s for m, s in ms.items()}
                    for n, ms in data['stereo_bonds'].items()}
    seen = set()
    for n, ns in stereo_bonds.items():
        if n in seen:
            continue
        if n in ctc:
            m = ctc[n]
            if m in stereo_bonds:
                seen.add(m)
                n2, s2 = stereo_bonds[m].popitem()
                n1, s1 = ns.popitem()
                stereo.append((molecule.add_cis_trans_stereo, n, m, n1, n2, s1 == s2))

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
            continue
        break

    if log:
        if 'chython_parsing_log' not in molecule.meta:
            molecule.meta['chython_parsing_log'] = log
        else:
            molecule.meta['chython_parsing_log'].extend(log)


__all__ = ['smiles']
