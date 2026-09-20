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
"""An organozinc or Grignard drawn apart from its halide: `CC[Zn+].[Cl-]` is one `CC[Zn]Cl`.

A ONE-COORDINATE ZINC OR MAGNESIUM HOLDING A CARBON IS AN INCOMPLETE DRAWING, and this stage completes
it two ways.  A free halide in the record is the halide that belongs on the metal, so it is bonded there.
A metal left without one is charged `+`: the halide is missing from the drawing rather than from the
compound, and a neutral one-coordinate metal is not a species anybody meant.  Both readings apply to
zinc and magnesium alike.

A `standardize()` stage rather than a `standardize_metals.tsv` row, and for one reason: a rule table
applies whichever match the isomorphism search returned first, so a drawing offering more than one
candidate would be settled by its atom order.  Pairing is a question about all the candidates at once.

Two orderings answer it, and neither can see the input's spelling:

* **halides** by `Cl > Br > I > F`, a fact about the reagents rather than about the graph;
* **metals** by canonical rank, which is `_isomers._ranks` without the stripping step -- the bond being
  placed is the one that is absent, so the molecule as it arrived already is the placement-free frame.

Ties are left tied.  Two candidates of equal rank are automorphic, so which one is taken is not
observable in the result, and the sort being stable makes the arbitrary half of the choice cheap.
Charging needs no order at all, being one atom's charge and nothing else's.

NET CHARGE MOVES HERE, and every spelling that moves it is accepted: a lone `[Zn+]` beside a neutral
halogen sits at `+1`, a neutral metal beside `[Cl-]` at `-1`, and a charged metal with no halide
leaves `0` for `+1`.  A dropped sign is the premise of the stage, so the record states which way the
total went rather than the pass declining to move it.
"""
from collections.abc import MutableSequence
from ..core import LogRecord, MoleculeContainer


__all__ = ['unite_organometallics']

#: Table-qualified, as every rule id must be -- never a bare index.  Two and not one: bonding a halide
#: the record HAS and charging a metal whose halide the record LACKS are different claims about the
#: drawing, and a consumer filtering the log can decline the second while keeping the first.
_RULE = 'organometallics:unite'
_RULE_CHARGE = 'organometallics:charge'

#: `Cl > Br > I > F`, low wins.  Chloride and bromide are the reagent halides; fluoride is last because
#: a free `[F-]` beside a magnesium is more often a separate salt than a bond somebody forgot to draw.
_PREFERENCE = {17: 0, 35: 1, 53: 2, 9: 3}

#: Zinc and magnesium only.  Widening this is a chemical claim per element, not a constant to grow.
_METALS = frozenset({12, 30})


def _candidates(molecule: MoleculeContainer) -> tuple[list[int], list[int]]:
    """The one-bonded metals holding a carbon, and the free halides, in arena order.

    A metal is a candidate at charge `0` or `+1` and at exactly one single bond to carbon: two bonds
    means it is already satisfied, and an alkoxide or amide is a different compound rather than a
    reagent drawn apart.  A halide is a candidate at charge `0` or `-1` with no bonds at all -- bonded,
    it belongs to whatever it is bonded to.  Radicals are nobody's dropped sign.
    """
    metals: list[int] = []
    halides: list[int] = []
    for atom in molecule.atoms():
        n = atom.n
        if atom.is_radical:
            continue
        if atom.element in _METALS:
            if atom.degree == 1 and atom.charge in (0, 1):
                m = next(iter(molecule.neighbors_of(n)))
                if molecule.atom(m).element == 6 and molecule.order_of(n, m) == 1:
                    metals.append(n)
        elif atom.element in _PREFERENCE and not atom.degree and atom.charge in (0, -1):
            halides.append(n)
    return metals, halides


def unite_organometallics(molecule: MoleculeContainer, log: MutableSequence) -> set[int]:
    """Bond each candidate metal to one candidate halide, and charge whichever went without.  Written ids.

    Spare halides are left as drawn -- a second chloride beside a satisfied metal is a counterion -- so
    only the metal side is completed either way.
    """
    metals, halides = _candidates(molecule)
    if not metals:
        return set()

    pairs: list[tuple[int, int, int]] = []
    if halides:
        # before the edit scope, which answers reads from the pre-scope arena.  `atoms_order` is the
        # expensive word here, so it is asked for only when there is a pairing to decide.
        ranks = molecule.atoms_order
        metals.sort(key=lambda n: ranks[n])
        halides.sort(key=lambda n: (_PREFERENCE[molecule.atom(n).element], molecule.atom(n).isotope,
                                    ranks[n]))
        pairs = [(n, m, molecule.charge_of(n) + molecule.charge_of(m)) for n, m in zip(metals, halides)]

    # a metal that took a halide is neutral by the join; one that did not is charged, unless the drawing
    # already said `+`.  Ranking cannot matter to this half -- every leftover is treated alike.
    bonded = {n for n, _, _ in pairs}
    stranded = [n for n in metals if n not in bonded and not molecule.charge_of(n)]
    if not pairs and not stranded:
        return set()

    with molecule.edit():
        for n, m, _ in pairs:
            molecule.add_bond(n, m, 1)
            molecule.set_charge(n, 0)
            molecule.set_charge(m, 0)
        for n in stranded:
            molecule.set_charge(n, 1)

    for n, m, before in pairs:
        moved = f', and the net charge of the pair moved {before:+d} -> 0' if before else ''
        log.append(LogRecord(_RULE, (n, m), f'atoms {n} and {m} are one organometallic reagent drawn '
                                            f'apart; joined by a single bond{moved}'))
    for n in stranded:
        log.append(LogRecord(_RULE_CHARGE, (n,), f'atom {n} holds one carbon and no halide; charged +1, '
                                                 f'the halide being absent from the drawing rather than '
                                                 f'from the compound, so the net charge moves 0 -> +1'))
    return bonded | {m for _, m, _ in pairs} | set(stranded)
