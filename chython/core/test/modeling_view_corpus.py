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
"""The reactions the frozen `modeling_view()` answer is recorded from, and the loader for it.

THE FIXTURE IS A PIN AND NOT A REGENERABLE ARTEFACT.  `gen_modeling_view_corpus.py` was run against
the dict-of-dicts implementation; `test_modeling_view_frozen.py` asserts the array kernel behind
`modeling_view()` still gives those answers.  Regenerating it makes the kernel agree with itself.

PUBLIC, TEXTBOOK REACTIONS, each named for the branch of the union it reaches: a fully mapped
transformation, a leaving fragment (reactant-only atoms), an arriving one (product-only), an atom with
no map number, two atoms sharing one, agents that must be excluded, a multi-molecule side, an aromatic
ring, a bond-order change and a ring closure.

ONE RECORD IS BUILT AND NOT PARSED.  An atom whose implicit hydrogen count is unknown cannot carry a
map number in SMILES: a map number needs brackets, and a bracket states the count.  So the H_UNKNOWN
arm of the state derivation is unreachable from a string and is fed by `_unknown_h_record` instead.
"""
import gzip
import json
from pathlib import Path

from chython.core import H_UNKNOWN, MoleculeContainer, ReactionContainer, read_reaction_smiles


__all__ = ['PATH', 'RECORDS', 'SMILES', 'load']

PATH = Path(__file__).parent / 'modeling_view_corpus.json.gz'

# name -> reaction SMILES.  Each comment names what the record pins.
SMILES = {
    # every atom mapped on both sides, one molecule per side
    'sn2_bromide_to_alcohol': '[CH3:1][Br:2].[OH-:3]>>[CH3:1][OH:3].[Br-:2]',
    # two reactants, two products, all mapped: the union spans four containers
    'esterification': '[CH3:1][C:2](=[O:3])[OH:4].[CH3:5][OH:6]'
                      '>>[CH3:1][C:2](=[O:3])[O:4][CH3:5].[OH2:6]',
    # an agent, which must contribute no atom and no bond
    'amidation_with_ethanol_agent': '[CH3:1][C:2](=[O:3])[OH:4].[NH2:5][CH3:6]'
                                    '>[CH3:10][CH2:11][OH:12]'
                                    '>[CH3:1][C:2](=[O:3])[NH:5][CH3:6].[OH2:4]',
    # a leaving fragment: 5..8 are mapped on the left and absent on the right, so they keep their
    # hydrogen counts and count only their bonds to each other
    'tert_butyl_ester_hydrolysis': '[CH3:1][C:2](=[O:3])[O:4][C:5]([CH3:6])([CH3:7])[CH3:8].[OH2:9]'
                                   '>>[CH3:1][C:2](=[O:3])[OH:9]',
    # the mirror: 5..8 exist only on the product side
    'tert_butyl_esterification': '[CH3:1][C:2](=[O:3])[OH:9]'
                                 '>>[CH3:1][C:2](=[O:3])[O:4][C:5]([CH3:6])([CH3:7])[CH3:8].[OH2:9]',
    # one unmapped atom on each side, counted and left out of the union
    'esterification_partly_mapped': '[CH3:1][C:2](=[O:3])[OH:4].CO'
                                    '>>[CH3:1][C:2](=[O:3])[O:4]C.O',
    # two atoms sharing map number 1 on the reactant side: the union merges them and says so
    'colliding_map_numbers': '[CH3:1][CH3:1]>>[CH3:2][CH3:3]',
    # an aromatic ring through the union, and a bond that changes order
    'benzene_nitration': '[cH:1]1[cH:2][cH:3][cH:4][cH:5][cH:6]1.[N+:7](=[O:8])([O-:9])[OH:10]'
                         '>>[c:1]1([N+:7](=[O:8])[O-:9])[cH:2][cH:3][cH:4][cH:5][cH:6]1.[OH2:10]',
    # two ring closures and four order changes in one record
    'diels_alder': '[CH2:1]=[CH:2][CH:3]=[CH2:4].[CH2:5]=[CH2:6]'
                   '>>[CH2:1]1[CH:2]=[CH:3][CH2:4][CH2:5][CH2:6]1',
    # a salt on the product side: two molecules, one of them a single mapped ion
    'saponification': '[CH3:1][C:2](=[O:3])[O:4][CH3:5].[OH-:6].[Na+:7]'
                      '>>[CH3:1][C:2](=[O:3])[O-:4].[Na+:7].[CH3:5][OH:6]',
    # no mapping at all: every atom is unmapped, the union is empty, and that is a record
    'unmapped_hydrogenation': 'C=C.[H][H]>>CC',
    # an empty product side, which is a record and not an error
    'empty_products': '[CH3:1][CH3:2]>>',
}


def _unknown_h_record():
    """A mapped reaction carrying an atom with no derivable hydrogen count.

    Built rather than parsed for the reason the module docstring gives.  Tetrafluoroammonium has no
    valence rule, so the reader leaves the nitrogen's count unknown; here the same state is written
    directly, and the map number makes the atom reach the union.
    """
    left = MoleculeContainer()
    with left.edit():
        n = left.add_atom('N', implicit_h=H_UNKNOWN)
        for _ in range(4):
            left.add_bond(n, left.add_atom('F', implicit_h=0), 1)
    left.set_map_number(n, 1)
    right = left.copy()
    return ReactionContainer([left], [right])


RECORDS = {name: read_reaction_smiles(line) for name, line in SMILES.items()}
RECORDS['tetrafluoroammonium_unknown_h'] = _unknown_h_record()


def load():
    """The recorded answers, keyed the same way as `RECORDS`."""
    with gzip.open(PATH, 'rt', encoding='utf8') as f:
        return json.load(f)
