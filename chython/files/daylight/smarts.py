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
from functools import partial
from itertools import count
from re import compile, findall
from typing import Union
from .parser import parser
from .tokenize import smarts_tokenize
from ...containers import QueryContainer, QueryBond, ReactionContainer
from ...periodictable import ListElement, QueryElement


cx_radicals = compile(r'\^[1-7]:[0-9]+(?:,[0-9]+)*')

# AD-HOC for masked atoms in SMARTS
# not thread safe
global_free_masked = count(10 ** 9 + 1)


def _parse_mol_smarts(data: str) -> QueryContainer:
    """Parse a single molecule SMARTS string into a QueryContainer."""
    smr, *cx = data.split()

    parsed = parser(smarts_tokenize(smr), False)

    if cx and cx[0].startswith('|') and cx[0].endswith('|'):
        for x in findall(cx_radicals, cx[0]):
            for i in x[3:].split(','):
                parsed['atoms'][int(i)]['is_radical'] = True

    for i, s in parsed['stereo_atoms'].items():
        parsed['atoms'][i]['stereo'] = s
    stereo_bonds = parsed['stereo_bonds']

    g = QueryContainer(data)

    mapping = {}
    default_aliphatic = {}  # internal-atom-number -> True if hybridization defaulted from uppercase symbol
    free = count(max(a.get('parsed_mapping', 0) for a in parsed['atoms']) + 1)
    for i, a in enumerate(parsed['atoms']):
        mapping[i] = n = a.pop('parsed_mapping', 0) or next(global_free_masked if a.get('masked') else free)
        e = a.pop('element')
        charge_not = a.pop('charge_not', None)
        recursive_smarts_raw = a.pop('recursive_smarts', None)
        excluded_elements = a.pop('excluded_elements', None)
        if a.pop('_default_aliphatic', False):
            default_aliphatic[n] = True
        if isinstance(e, int):
            e = QueryElement.from_atomic_number(e)
        elif isinstance(e, str):
            e = QueryElement.from_symbol(e)
        else:
            e = partial(ListElement, e)
        g.add_atom(e(**a), n)
        if charge_not:
            g.atom(n)._charge_not = charge_not
        if excluded_elements:
            g.atom(n)._excluded_elements = excluded_elements
        if recursive_smarts_raw:
            compiled = []
            for positive, inner_str in recursive_smarts_raw:
                sub_query = smarts(inner_str)
                root = next(iter(sub_query))
                compiled.append((positive, sub_query, root))
            g.atom(n)._recursive_smarts = compiled

    for n, m, b in parsed['bonds']:
        if n in stereo_bonds and m in stereo_bonds:
            if m not in stereo_bonds[n]:  # only simple cis-trans supported, not cumulenes
                _, s1 = stereo_bonds[n].popitem()
                _, s2 = stereo_bonds[m].popitem()
                if isinstance(b, int):
                    b = QueryBond(b, stereo=s1 == s2)
                else:
                    b.stereo = s1 == s2
        g.add_bond(mapping[n], mapping[m], b)

    # Relax aliphatic-default if the atom has an aromatic bond — preserves
    # legacy ``[C]:[C]`` idioms while keeping strict ``[C]`` semantics.
    if default_aliphatic:
        for atom_n in default_aliphatic:
            for neighbour, bond in g._bonds[atom_n].items():
                bond_orders = bond.order
                if bond_orders == (4,):
                    g._atoms[atom_n]._hybridization = ()
                    break
    return g


def smarts(data: str) -> Union[QueryContainer, ReactionContainer]:
    """
    Parse SMARTS or reaction SMARTS string.

    For molecule SMARTS returns QueryContainer.
    For reaction SMARTS (containing `>>`) returns ReactionContainer composed of QueryContainers.

    Reaction SMARTS format: ``reactants>reagents>products`` or ``reactants>>products``

    * only D, a, h, r and !R atom primitives supported.
    * bond order list (max 2) and not bond supported.
    * [not]ring bond supported only in combination with explicit bonds, not bonds and bonds orders lists.
    * mapping, charge and isotopes supported.
    * list of elements supported.
    * A - treats as any element. <A> primitive (aliphatic) ignored.
    * M - treats as any metal
    * <&> logic operator unsupported.
    * <;> logic operator is mandatory except (however preferable) for charge, isotope, stereo marks.
    * CXSMARTS radicals supported.
    * masked atom - `chython.Reactor` specific mark for masking reactant atoms from deletion.

    For example::

        [C;r5,r6;a]-;!@[C;h1,h2;z2,z4] |^1:1| - aromatic C member of 5 or 6 atoms ring
        connected with non-ring single bond to aromatic or SP2 radical C with 1 or 2 hydrogens.

    * primitive <xN> - heteroatoms (e.g. x2 - two heteroatoms)
    * primitive <zN> - hybridization (N = 1 - sp3, 2 - sp2, 3 - sp, 4 - aromatic)
    * primitive <M> - masked atom

    Note: atom numbers greater than 10 ** 9 forbidden for usage and reserved for masked atoms numbering.
    In multiprocess mode has potential bugs in reaction enumeration task then used templates prepared from components
    from different processes. For avoiding, prepare templates on single process and then share it.
    """
    if not isinstance(data, str):
        raise TypeError('Must be a SMARTS string')

    if '>' in data:
        # Reaction-level trailing CX block: indices are absolute across
        # reactants → reagents → products.
        body, rxn_rad_atoms = _split_trailing_cx(data)
        try:
            reactants_str, reagents_str, products_str = body.split('>')
        except ValueError as e:
            raise ValueError('invalid reaction SMARTS') from e

        reactants = _parse_smarts_side(reactants_str)
        reagents = _parse_smarts_side(reagents_str)
        products = _parse_smarts_side(products_str)

        if rxn_rad_atoms:
            cum = 0
            for side_mols in (reactants, reagents, products):
                for mol in side_mols:
                    atom_ids = list(mol._atoms.keys())
                    for pos, aid in enumerate(atom_ids):
                        if (cum + pos) in rxn_rad_atoms:
                            mol._atoms[aid]._is_radical = True
                    cum += len(atom_ids)
        return ReactionContainer(reactants, products, reagents)

    return _parse_mol_smarts(data)


def _split_trailing_cx(s: str):
    """Strip trailing CX block; return (body, radical_indices_set)."""
    s = s.strip()
    if not s.endswith('|'):
        return s, set()
    open_pipe = s.rfind('|', 0, -1)
    if open_pipe <= 0 or s[open_pipe - 1] != ' ':
        return s, set()
    cx_block = s[open_pipe:]
    body = s[:open_pipe - 1]
    if '|' in body:  # legacy per-fragment CX present; leave body alone
        return s, set()
    rad = set()
    for x in findall(cx_radicals, cx_block):
        for i in x[3:].split(','):
            rad.add(int(i))
    return body, rad


def _parse_smarts_side(side: str) -> list[QueryContainer]:
    """Parse one reaction side. Supports legacy per-fragment CX
    (``[A] |^1:0|.[B]``) and side-local CX (``[A].[B] |^1:N|``)."""
    side, rad_atoms = _split_trailing_cx(side)
    if not side:
        return []
    frags = [f.strip() for f in side.split('.') if f.strip()]
    mols = [_parse_mol_smarts(f) for f in frags]
    if rad_atoms:
        cum = 0
        for mol in mols:
            atom_ids = list(mol._atoms.keys())
            for pos, aid in enumerate(atom_ids):
                if (cum + pos) in rad_atoms:
                    mol._atoms[aid]._is_radical = True
            cum += len(atom_ids)
    return mols


__all__ = ['smarts']
