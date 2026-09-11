# -*- coding: utf-8 -*-
#
#  Copyright 2026 Ramil Nugmanov <nougmanoff@protonmail.com>
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
"""Read `tables/residues.tsv` and hand back the hardcoded connectivity of one residue.

A PDB-family file states which residue a group of atoms belongs to and almost never states a bond
between them; a distance cutoff being prohibited here, the residues where the answer is a fact are
written down instead.  Knowledge only, applying nothing: the pass that consumes a `PDBRecord` decides
what to do when a file and a template disagree.  Loading is lazy and cached.
"""
from collections.abc import Mapping
from types import MappingProxyType
from typing import NamedTuple
from ._tables import read_table


__all__ = ['RESIDUE_KINDS', 'ResidueTemplate', 'normalize_atom_name', 'residue_template',
           'residue_templates']


#: The closed vocabulary of the `kind` column; anything else is a load-time error.  A caller branches
#: on `kind` to decide whether a residue links into a chain, so an unknown spelling falls through.
RESIDUE_KINDS = ('amino_acid', 'nucleotide', 'water', 'ion')

#: The kinds that are polymer residues, and therefore the ones obliged to name both link atoms.
_POLYMER_KINDS = frozenset({'amino_acid', 'nucleotide'})

#: Legacy spellings of nucleotide phosphate oxygen names that no character rule reaches.  Applied
#: only for nucleotide rows: `O1P`, `O2P` and `O3P` are live atom names in phosphorylated residues
#: such as SEP, TPO and AMP, so a table-wide alias would corrupt any such row the table grows.
_ALIASES = {'O1P': 'OP1', 'O2P': 'OP2', 'O3P': 'OP3'}


class ResidueTemplate(NamedTuple):
    """One row of `tables/residues.tsv`, compiled.

    `atoms` maps the atom name to `(element symbol, formal charge)`, heavy atoms only -- the implicit
    count is derived from the valence rules afterwards, which is what makes a mid-chain and a terminal
    residue come out right off one row.

    `bonds` are `(name_a, name_b, order)` and are Kekule: the aromatic rings of HIS, PHE, TRP, TYR and
    the nucleobases carry alternating 1/2, so an aromatic form is `thiele()`'s to produce.

    `link_in` and `link_out` name the atoms bonding to the preceding and following residue of a chain,
    two names and no leaving group.  Both are `None` for water and for an ion.
    """
    name: str
    kind: str
    atoms: Mapping[str, tuple[str, int]]
    bonds: tuple[tuple[str, str, int], ...]
    link_in: str | None
    link_out: str | None


_RESIDUES_CACHE: dict[str, Mapping[str, ResidueTemplate]] = {}


def normalize_atom_name(name: str, kind: str | None = None) -> str:
    """The one spelling of an atom name this table is keyed by.

    Upper-cased, stripped, `*` folded to `'` (legacy files spell the ribose oxygens `O3*`, `O5*`), and
    then `_ALIASES` for nucleotide rows only.  The `kind` gate is load-bearing: `O1P`, `O2P` and `O3P`
    are current CCD names in phosphorylated residues such as SEP and TPO, so a table-wide alias would
    rename a real atom to one its own row does not have.

    `OW` (GROMACS) and `OT1`/`OT2` (CHARMM) are deliberately not aliased -- neither can be renamed
    safely across every row, so the consuming pass handles those vocabularies itself.

    The table's own keys are normalised through this, so there is one spelling and not two that agree.
    """
    name = name.strip().upper().replace('*', "'")
    if kind == 'nucleotide':
        return _ALIASES.get(name, name)
    return name


def _compile() -> Mapping[str, ResidueTemplate]:
    out: dict[str, ResidueTemplate] = {}
    for row in read_table('residues.tsv'):
        name = row['name'].strip().upper()
        if name in out:
            raise ValueError(f'residues.tsv: {name} appears twice; a component id is the only handle '
                             'a caller has on a row and must name one')
        kind = row['kind']
        if kind not in RESIDUE_KINDS:
            raise ValueError(f'{name}: kind {kind!r} is not one of {", ".join(RESIDUE_KINDS)}')

        atoms: dict[str, tuple[str, int]] = {}
        for entry in row['atoms'].split(','):
            parts = entry.split(':')
            if len(parts) == 2:
                atom_name, element = parts
                charge = 0
            elif len(parts) == 3:
                atom_name, element, charge_text = parts
                charge = int(charge_text)
            else:
                raise ValueError(f'{name}: atom entry {entry!r} is not NAME:ELEMENT or '
                                 'NAME:ELEMENT:CHARGE')
            atom_name = normalize_atom_name(atom_name, kind)
            if atom_name in atoms:
                raise ValueError(f'{name}: atom {atom_name} appears twice.  Names are how a file\'s '
                                 'atoms are matched to this row, so a duplicate makes one of the two '
                                 'unreachable')
            if kind in _POLYMER_KINDS and charge:
                raise ValueError(f'{name}: {atom_name} carries charge {charge}, but a polymer residue '
                                 'is written as its neutral free component -- a PDB file states no '
                                 'protonation state, so a charge here would be invented')
            atoms[atom_name] = (element, charge)

        bonds = []
        seen_pairs = set()
        for entry in row['bonds'].split(',') if row['bonds'] else ():
            parts = entry.split('-')
            if len(parts) != 3:
                raise ValueError(f'{name}: bond entry {entry!r} is not NAME_A-NAME_B-ORDER')
            a, b, order_text = parts
            a = normalize_atom_name(a, kind)
            b = normalize_atom_name(b, kind)
            for atom_name in (a, b):
                if atom_name not in atoms:
                    raise ValueError(f'{name}: bond {entry!r} names atom {atom_name}, which the '
                                     'atoms column does not declare')
            if a == b:
                raise ValueError(f'{name}: bond {entry!r} joins an atom to itself')
            pair = (a, b) if a < b else (b, a)
            if pair in seen_pairs:
                raise ValueError(f'{name}: atoms {pair[0]} and {pair[1]} are bonded twice; a second '
                                 'order for one pair is two claims about it and this row must make '
                                 'one')
            seen_pairs.add(pair)
            order = int(order_text)
            # order 4 is refused, not merely unused: aromatizing at read time is `thiele()`'s decision
            if order not in (1, 2, 3):
                raise ValueError(f'{name}: bond order {order} is not one of 1 2 3.  Rings are Kekule '
                                 'in this table; an aromatic form is what thiele() produces')
            bonds.append((a, b, order))

        links = []
        for column in ('link_in', 'link_out'):
            cell = row[column].strip()
            if not cell:
                if kind in _POLYMER_KINDS:
                    raise ValueError(f'{name}: a {kind} row must name {column}.  A residue that links '
                                     'into no chain cannot be joined to its neighbours, and the join '
                                     'is the whole reason a polymer needs this table')
                links.append(None)
                continue
            link = normalize_atom_name(cell, kind)
            if link not in atoms:
                raise ValueError(f'{name}: {column} names {link}, which the atoms column does not '
                                 'declare')
            if kind not in _POLYMER_KINDS:
                raise ValueError(f'{name}: {column} names {link}, but a {kind} row has no neighbour '
                                 'to link to -- only a polymer residue does')
            links.append(link)

        out[name] = ResidueTemplate(name, kind, MappingProxyType(atoms), tuple(bonds), links[0], links[1])
    return out


def residue_templates() -> Mapping[str, ResidueTemplate]:
    """Every row of `tables/residues.tsv`, keyed by component id, loaded lazily on first use."""
    if 'rows' not in _RESIDUES_CACHE:
        _RESIDUES_CACHE['rows'] = MappingProxyType(_compile())
    return _RESIDUES_CACHE['rows']


def residue_template(name: str) -> ResidueTemplate | None:
    """The template for one component id, or `None` when the table does not have one.

    A miss and not a raise: an unrecognised residue is the common case in any real structure, and one
    ligand must not abort a file that read perfectly well.

    The id is accepted in any case and stripped, a legacy PDB residue-name field arriving padded.
    """
    return residue_templates().get(name.strip().upper())
