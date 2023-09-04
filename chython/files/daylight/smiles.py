# -*- coding: utf-8 -*-
#
#  Copyright 2022, 2023 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from re import compile, findall, search
from typing import List, Union, Optional
from .parser import parser
from .tokenize import smiles_tokenize
from .._convert import create_molecule, create_reaction
from .._mapping import postprocess_parsed_molecule, postprocess_parsed_reaction
from ...containers import MoleculeContainer, ReactionContainer
from ...exceptions import IsChiral, NotChiral, ValenceError


cx_fragments = compile(r'f:(?:[0-9]+(?:\.[0-9]+)+)(?:,(?:[0-9]+(?:\.[0-9]+)+))*')
cx_radicals = compile(r'\^[1-7]:[0-9]+(?:,[0-9]+)*')


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
    if not data:
        raise ValueError('Empty string')

    contract: Optional[List[List[int]]]  # typing
    log = []
    smi, *data = data.split()
    if data and (cxs := data[0]).startswith('|') and cxs.endswith('|'):
        fr = search(cx_fragments, cxs)
        if fr is not None:
            contract = [sorted(int(x) for x in x.split('.')) for x in fr.group()[2:].split(',')]
            if len({x for x in contract for x in x}) < len([x for x in contract for x in x]):
                log.append(f'collisions in cxsmiles fragments description: {cxs}')
                contract = None

            radicals = [int(x) for x in findall(cx_radicals, cxs) for x in x[3:].split(',')]
            if radicals and len(set(radicals)) != len(radicals):
                log.append(f'collisions in cxsmiles radicals description: {cxs}')
                radicals = []
        else:
            radicals = [int(x) for x in findall(cx_radicals, cxs) for x in x[3:].split(',')]
            if radicals:
                if len(set(radicals)) != len(radicals):
                    log.append(f'collisions in cxsmiles radicals description: {cxs}')
                    radicals = []
            contract = None
    else:
        radicals = []
        contract = None

    if '>' in smi:
        record = {'reactants': [], 'reagents': [], 'products': [], 'log': log, 'meta': None, 'title': None}
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

        if radicals:
            atom_map = dict(enumerate(a for m in chain(record['reactants'], record['reagents'], record['products'])
                                      for a in m['atoms']))
            for x in radicals:
                atom_map[x]['is_radical'] = True

        postprocess_parsed_reaction(record, remap=remap, ignore=ignore)
        rxn = create_reaction(record, ignore_bad_isotopes=ignore_bad_isotopes, _r_cls=_r_cls, _m_cls=_m_cls)
        for mol, tmp in zip(rxn.molecules(), chain(record['reactants'], record['reagents'], record['products'])):
            postprocess_molecule(mol, tmp, ignore=ignore, ignore_stereo=ignore_stereo,
                                 ignore_carbon_radicals=ignore_carbon_radicals, keep_implicit=keep_implicit,
                                 ignore_aromatic_radicals=ignore_aromatic_radicals)
        return rxn
    else:
        record = parser(smiles_tokenize(smi), not ignore)
        for x in radicals:
            record['atoms'][x]['is_radical'] = True
        record['log'].extend(log)

        postprocess_parsed_molecule(record, remap=remap, ignore=ignore)
        mol = create_molecule(record, ignore_bad_isotopes=ignore_bad_isotopes, _cls=_m_cls)
        postprocess_molecule(mol, record, ignore=ignore, ignore_stereo=ignore_stereo,
                             ignore_carbon_radicals=ignore_carbon_radicals, keep_implicit=keep_implicit,
                             ignore_aromatic_radicals=ignore_aromatic_radicals)
        return mol


def postprocess_molecule(molecule, data, *, ignore=True, ignore_stereo=False, ignore_carbon_radicals=False,
                         keep_implicit=False, ignore_aromatic_radicals=True):
    mapping = data['mapping']

    atoms = molecule._atoms
    bonds = molecule._bonds
    charges = molecule._charges
    hydrogens = molecule._hydrogens
    radicals = molecule._radicals
    hyb = molecule.hybridization
    radicalized = []

    implicit_mismatch = {}
    if 'chython_parsing_log' in molecule.meta:
        log = molecule.meta['chython_parsing_log']
    else:
        log = []

    for n, a in enumerate(data['atoms']):
        h = a['hydrogen']
        if h is None:  # simple atom token
            continue
        # bracket token should always contain implicit hydrogens count.
        n = mapping[n]
        if keep_implicit:  # override any calculated hydrogens count.
            hydrogens[n] = h
        elif (hc := hydrogens[n]) is None:  # atom has invalid valence or aromatic ring.
            if hyb(n) == 4:  # this is aromatic rings. just store given H count.
                hydrogens[n] = h
                # rare H0 case
                if (not ignore_aromatic_radicals and not h and not charges[n] and not radicals[n] and
                    atoms[n].atomic_number in (5, 6, 7, 15) and sum(b.order != 8 for b in bonds[n].values()) == 2):
                    # c[c]c - aromatic B,C,N,P radical
                    radicals[n] = True
                    radicalized.append(n)
            elif not radicals[n]:  # CXSMILES radical not set.
                # SMILES doesn't code radicals. so, let's try to guess.
                radicals[n] = True
                if molecule._check_implicit(n, h):  # radical form is valid
                    radicalized.append(n)
                    hydrogens[n] = h
                elif ignore:  # radical state also has errors.
                    radicals[n] = False  # reset radical state
                    implicit_mismatch[n] = h
                    log.append(f'implicit hydrogen count ({h}) mismatch with calculated on atom {n}')
                else:
                    raise ValueError(f'implicit hydrogen count ({h}) mismatch with calculated on atom {n}')
        elif hc != h:  # H count mismatch.
            if hyb(n) == 4:
                if not h and not charges[n] and not radicals[n] and atoms[n].atomic_number in (5, 6, 7, 15) and \
                        sum(b.order != 8 for b in bonds[n].values()) == 2:
                    # c[c]c - aromatic B,C,N,P radical
                    hydrogens[n] = 0
                    radicals[n] = True
                    radicalized.append(n)
                elif ignore:
                    implicit_mismatch[n] = h
                    log.append(f'implicit hydrogen count ({h}) mismatch with calculated on atom {n}')
                else:
                    raise ValueError(f'implicit hydrogen count ({h}) mismatch with calculated on atom {n}')
            elif molecule._check_implicit(n, h):  # set another possible implicit state. probably Al, P
                hydrogens[n] = h
            elif not radicals[n]:  # CXSMILES radical is not set. try radical form
                radicals[n] = True
                if molecule._check_implicit(n, h):
                    hydrogens[n] = h
                    radicalized.append(n)
                # radical state also has errors.
                elif ignore:
                    radicals[n] = False  # reset radical state
                    implicit_mismatch[n] = h
                    log.append(f'implicit hydrogen count ({h}) mismatch with calculated on atom {n}')
                else:
                    raise ValueError(f'implicit hydrogen count ({h}) mismatch with calculated on atom {n}')
            elif ignore:  # just ignore it
                implicit_mismatch[n] = h
                log.append(f'implicit hydrogen count ({h}) mismatch with calculated on atom {n}')
            else:
                raise ValueError(f'implicit hydrogen count ({h}) mismatch with calculated on atom {n}')

    if ignore_carbon_radicals:
        for n in radicalized:
            if atoms[n].atomic_number == 6:
                radicals[n] = False
                hydrogens[n] += 1
                log.append(f'carbon radical {n} replaced with implicit hydrogen')

    if implicit_mismatch:
        molecule.meta['chython_implicit_mismatch'] = implicit_mismatch
    if log and 'chython_parsing_log' not in molecule.meta:
        molecule.meta['chython_parsing_log'] = log
    if ignore_stereo:
        return

    stereo_atoms = [(n, s) for n, a in enumerate(data['atoms']) if (s := a['stereo']) is not None]
    if not stereo_atoms and not data['stereo_bonds']:
        return

    st = molecule._stereo_tetrahedrons
    sa = molecule._stereo_allenes
    sat = molecule._stereo_allenes_terminals
    ctt = molecule._stereo_cis_trans_terminals

    order = {mapping[n]: [mapping[m] for m in ms] for n, ms in data['order'].items()}

    stereo = []
    for i, s in stereo_atoms:
        n = mapping[i]
        if not i and hydrogens[n]:  # first atom in smiles has reversed chiral mark
            s = not s

        if n in st:
            stereo.append((molecule.add_atom_stereo, n, order[n], s))
        elif n in sa:
            t1, t2 = sat[n]
            env = sa[n]
            n1 = next(x for x in order[t1] if x in env)
            n2 = next(x for x in order[t2] if x in env)
            stereo.append((molecule.add_atom_stereo, n, (n1, n2), s))

    stereo_bonds = {mapping[n]: {mapping[m]: s for m, s in ms.items()}
                    for n, ms in data['stereo_bonds'].items()}
    seen = set()
    for n, ns in stereo_bonds.items():
        if n in seen:
            continue
        if n in ctt:
            nm = ctt[n]
            m = nm[1] if nm[0] == n else nm[0]
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

    if log and 'chython_parsing_log' not in molecule.meta:
        molecule.meta['chython_parsing_log'] = log


__all__ = ['smiles']
