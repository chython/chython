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
"""Atom-atom mapping from a transformer's attention over the two sides of a reaction.

THE MODEL AND NOTHING ELSE.  No rule-based repair runs after it, so the mapper's own accuracy is a
number a caller can obtain; composing it with a fixer is the caller's next line.  chython 2 ended this
function with `if self.fix_mapping(): fixed = True` and offered no flag to stop it, which left the
neural result and two rule fixers' results indistinguishable in one bool.

THE ONLY WRITE IS `set_map_number`.  Atom ids, bonds, charges and the structures come back untouched.
chython 2 mapped by RENUMBERING atoms, so it juggled three `remap()` calls and a flag to keep the
reactant side's numbering; here a map number is its own field and an atom id is stable, so numbering and
mapping are not the same operation and only the second one happens.

NUMPY AND ONNX RUNTIME ARE IMPORTED INSIDE THE FUNCTION.  `chython.reactions` imports this module at its
own init to register the body, and an 80 MiB model behind a module-level import would be loaded by
`import chython`.  `chython/reactions/test/test_attention_isolation.py` is the gate.
"""
from itertools import count

from ._session import attention_available, default_threads
from ...core import LOST, MappingResult, REFUSED, recording


def attention_mapping(reaction, *, multiplier: float = 1.75, keep_reactant_mapping: bool = False,
                      threads: int | None = None) -> MappingResult:
    """See `ReactionContainer.attention_mapping`, whose body this is."""
    from ._assign import greedy_mapping, side_adjacency
    from ._encode import MAX_NEIGHBORS, encode_reaction, run_model

    reactants, products, agents = reaction.reactants, reaction.products, reaction.agents
    with recording(reaction, stage='attention_mapping') as log:
        if not reactants or not products:
            log.record('a record with no reactants or no products has no correspondence to find',
                       severity=REFUSED, rule='attention:empty')
            return MappingResult(False, 0., (), 'empty')

        # HEAVY DEGREE ALONE, and 14 here is not the 14 the `neighbors` column clamps at: that clamp is
        # on degree plus hydrogens and merely saturates, while an atom past this bound has no token the
        # weights ever saw.  Agents are included because chython 2 included them.
        hypervalent = sum(1 for molecule in reaction.molecules() for n in molecule.atom_numbers
                          if molecule.degree_of(n) > MAX_NEIGHBORS)
        if hypervalent:
            log.record('%d atom(s) carry more than %d heavy neighbours, which is outside the domain the '
                       'weights were trained on; no map number is written'
                       % (hypervalent, MAX_NEIGHBORS), severity=REFUSED, rule='attention:hypervalent')
            return MappingResult(False, 0., (), 'hypervalent')

        r_ids = [molecule.atom_numbers for molecule in reactants]
        p_ids = [molecule.atom_numbers for molecule in products]
        before = _snapshot(reaction)

        encoded = encode_reaction(reactants, products)
        attention = run_model(encoded, default_threads() if threads is None else threads)
        assignment, score = greedy_mapping(attention, side_adjacency(reactants),
                                           side_adjacency(products), multiplier)

        if keep_reactant_mapping:
            r_numbers = [molecule.map_number_of(n) for molecule, ids in zip(reactants, r_ids)
                         for n in ids]
            fresh = count(max(r_numbers, default=0) + 1)
            # A REACTANT ATOM CARRYING NO NUMBER STILL GETS ONE.  Leaving the hole and giving the product
            # atom matched to it a fresh number writes a correspondence to an atom that has no such
            # number -- a mapping to nowhere, which is worse than the hole it preserves.
            r_numbers = [number or next(fresh) for number in r_numbers]
        else:
            r_numbers = list(range(1, sum(len(ids) for ids in r_ids) + 1))
            fresh = count(len(r_numbers) + 1)

        p_numbers = []
        unplaced = []
        position = 0
        for index, ids in enumerate(p_ids):
            for n in ids:
                matched = assignment[position]
                position += 1
                if matched < 0:
                    p_numbers.append(0)
                    unplaced.append((index, n))
                else:
                    p_numbers.append(r_numbers[matched])

        _write(reactants, r_ids, r_numbers)
        _write(products, p_ids, p_numbers)
        if not keep_reactant_mapping:
            # Numbered although never modelled: a record may carry a mapped catalyst, and an agent with
            # no number at all in an otherwise mapped record is a hole a consumer has to interpret.
            a_ids = [molecule.atom_numbers for molecule in agents]
            _write(agents, a_ids, [next(fresh) for ids in a_ids for _ in ids])

        log.record('%d of %d product atoms placed, mean attention %.3f'
                   % (len(p_numbers) - len(unplaced), len(p_numbers), score),
                   rule='attention:score')
        if unplaced:
            log.record('%d product atom(s) had no correspondence left and keep map number 0: %s'
                       % (len(unplaced),
                          ', '.join('products[%d] atom %d' % pair for pair in unplaced)),
                       severity=LOST, rule='attention:unplaced')
        return MappingResult(_snapshot(reaction) != before, score, tuple(unplaced), None)


def _snapshot(reaction):
    """Every map number in the record, in `molecules()` order -- what `changed` is measured against."""
    return [[molecule.map_number_of(n) for n in molecule.atom_numbers]
            for molecule in reaction.molecules()]


def _write(molecules, ids, numbers):
    """`numbers` is one flat list in the order `ids` walks the molecules.

    One edit session per molecule, so the arena reseals once rather than once per atom.
    """
    position = 0
    for molecule, own in zip(molecules, ids):
        with molecule.edit():
            for n in own:
                molecule.set_map_number(n, numbers[position])
                position += 1


__all__ = ['attention_available', 'attention_mapping']
