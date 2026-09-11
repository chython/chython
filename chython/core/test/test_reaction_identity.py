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
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with this program; if not, see <https://www.gnu.org/licenses/>.
#
"""
`ReactionContainer.__eq__` / `__hash__`.

Every assertion here is a SEMANTIC decision, not an implementation detail, which is why they are
tested one per decision with the decision named: the module docstring of `core/reaction.py` argues
that a nearly-right reaction equality is worse than none, and a test file that only checked "equal
things are equal" would let any of the three answers be reversed silently.

The decisions, and the test that pins each:

  * agents PARTICIPATE                        `test_agents_participate`
  * atom-to-atom mapping does NOT             `test_mapping_does_not_participate`
  * `title` and `meta` do NOT                 `test_title_and_meta_do_not_participate`
  * sides are compared SEPARATELY             `test_sides_are_not_pooled`
  * multiplicity counts (multiset, not set)   `test_multiplicity_participates`
  * order within a side does not              `test_order_within_a_side_does_not_participate`

Public compounds only, and small ones: the canonical search runs per molecule per comparison.
"""
from chython.core._core import read_smiles
from chython.core.reaction import ReactionContainer


def _esterification(*, catalyst=None, mapped=False, title=b'', meta=None):
    """Acetic acid + ethanol -> ethyl acetate + water, optionally with a catalyst and a mapping.

    A textbook Fischer esterification: public, small, and it has a genuine agent to put on the
    middle side, which is what the agents decision needs.
    """
    if mapped:
        reactants = [read_smiles('[CH3:1][C:2](=[O:3])[OH:4]'), read_smiles('[CH3:5][CH2:6][OH:7]')]
        products = [read_smiles('[CH3:1][C:2](=[O:3])[O:7][CH2:6][CH3:5]'), read_smiles('[OH2:4]')]
    else:
        reactants = [read_smiles('CC(=O)O'), read_smiles('CCO')]
        products = [read_smiles('CC(=O)OCC'), read_smiles('O')]
    agents = [read_smiles(catalyst)] if catalyst else []
    return ReactionContainer(reactants, products, agents, title=title, meta=meta)


# --- the decisions ---------------------------------------------------------------------------------

def test_agents_participate():
    """Ruled 2026-09-03: a reaction is a record of what was done, so a catalyst makes it a different
    reaction.  Sulfuric acid is the classic esterification catalyst."""
    plain = _esterification()
    catalysed = _esterification(catalyst='OS(=O)(=O)O')

    assert plain != catalysed
    assert hash(plain) != hash(catalysed)

    # and two runs with the SAME catalyst are one reaction
    assert catalysed == _esterification(catalyst='OS(=O)(=O)O')


def test_the_transformation_alone_is_still_spellable():
    """The reading agents-are-circumstance is not lost by the ruling above -- the sides are compared
    separately, so dropping the middle one is one expression.  This is the escape hatch the docstring
    promises a caller who wants transformation identity rather than record identity."""
    plain = _esterification()
    catalysed = _esterification(catalyst='OS(=O)(=O)O')

    assert plain != catalysed                                        # as records
    assert plain._identity()[0::2] == catalysed._identity()[0::2]    # as transformations


def test_mapping_does_not_participate():
    """Ruled 2026-09-03: "mapping is not needed for reaction comparison."

    This holds for free rather than by arrangement -- a map number is not in the atom invariant word
    `mol_identity_bytes` reads -- and that is exactly why it needs a test: nothing in `__eq__`
    mentions mapping, so nothing in `__eq__` would break if the molecule layer started counting it.
    """
    unmapped = _esterification()
    mapped = _esterification(mapped=True)

    # the mapping really is there, or this test proves nothing
    assert any(a.map_number for m in mapped.molecules() for a in m.atoms())
    assert not any(a.map_number for m in unmapped.molecules() for a in m.atoms())

    assert unmapped == mapped
    assert hash(unmapped) == hash(mapped)


def test_mapping_metric_cannot_route_through_equality():
    """The cost of the ruling above, stated as a test so nobody rediscovers it as a bug.

    Two DIFFERENT mappings of one reaction compare equal, so a mapping-quality harness written as
    `produced == reference` would report a perfect score while measuring nothing at all.  The mapping
    epic must compare `map_number` explicitly.
    """
    reference = _esterification(mapped=True)

    # the same reaction, mapped differently: every map number shifted by 10
    other = _esterification(mapped=True)
    for mol in other.molecules():
        for atom in mol.atoms():
            if atom.map_number:
                mol.atom(atom.n).map_number = atom.map_number + 10

    numbers_ref = sorted(a.map_number for m in reference.molecules() for a in m.atoms())
    numbers_other = sorted(a.map_number for m in other.molecules() for a in m.atoms())
    assert numbers_ref != numbers_other      # the annotations differ
    assert reference == other                # the reactions do not


def test_title_and_meta_do_not_participate():
    """Not chemistry.  One record read from two files under two names is one reaction."""
    a = _esterification(title=b'ester-001', meta={'source': 'file-a'})
    b = _esterification(title=b'a completely different name', meta={'source': 'file-b'})

    assert a == b
    assert hash(a) == hash(b)


def test_sides_are_not_pooled():
    """Three tuples and not one pooled multiset: which side a molecule is on is chemistry.

    Run forwards and backwards, the same two molecules are two different reactions -- a pooled
    comparison would call them equal.
    """
    forward = ReactionContainer([read_smiles('CC=O')], [read_smiles('CCO')])
    reverse = ReactionContainer([read_smiles('CCO')], [read_smiles('CC=O')])

    assert forward != reverse
    assert hash(forward) != hash(reverse)


def test_multiplicity_participates():
    """A multiset, not a set: `2 A -> B` is not `A -> B`.  This is what rules out `frozenset`."""
    once = ReactionContainer([read_smiles('CCO')], [read_smiles('CCOCC')])
    twice = ReactionContainer([read_smiles('CCO'), read_smiles('CCO')], [read_smiles('CCOCC')])

    assert once != twice
    assert hash(once) != hash(twice)


def test_order_within_a_side_does_not_participate():
    """`A + B` and `B + A` are one reaction: a side is a multiset, so it is sorted before comparison.

    Two molecules whose canonical bytes sort in a known-unequal order, so this cannot pass by the
    two happening to be identical.
    """
    a, b = read_smiles('CC(=O)O'), read_smiles('CCO')
    assert a.canonical_bytes != b.canonical_bytes

    ab = ReactionContainer([a, b], [read_smiles('CC(=O)OCC'), read_smiles('O')])
    ba = ReactionContainer([b, a], [read_smiles('CC(=O)OCC'), read_smiles('O')])

    assert ab == ba
    assert hash(ab) == hash(ba)


# --- protocol --------------------------------------------------------------------------------------

def test_usable_as_a_dict_key_and_set_member():
    """The point of the whole decision: a reaction can be deduplicated.

    Four records, of which two are the same reaction under different titles and one differs only by
    its catalyst -- so a correct set holds three.
    """
    records = [_esterification(title=b'first'),
               _esterification(title=b'second'),
               _esterification(catalyst='OS(=O)(=O)O'),
               ReactionContainer([read_smiles('CC=O')], [read_smiles('CCO')])]

    assert len(set(records)) == 3
    counts = {}
    for rxn in records:
        counts[rxn] = counts.get(rxn, 0) + 1
    assert sorted(counts.values()) == [1, 1, 2]


def test_identity_is_stable_across_copy():
    """`copy()` rebuilds every molecule, so equality must survive it or nothing above is reliable."""
    original = _esterification(catalyst='OS(=O)(=O)O', title=b'x')
    duplicate = original.copy()

    assert original == duplicate
    assert hash(original) == hash(duplicate)


def test_comparison_with_a_non_reaction():
    """`NotImplemented`, not `False`, so Python can try the other operand and `!=` stays consistent."""
    rxn = _esterification()

    assert rxn.__eq__('not a reaction') is NotImplemented
    assert rxn != 'not a reaction'
    assert not rxn == 'not a reaction'
    assert rxn != None  # noqa: E711 -- the operator is the subject of the test


def test_identical_object_is_equal_without_a_canonical_search():
    """`self is other` short-circuits.  An empty reaction is the cheapest witness that the fast path
    is taken at all: it is equal to itself either way, so this only documents the intent -- the
    measurable part is that it does not raise on a record with no sides."""
    rxn = ReactionContainer()

    assert rxn == rxn
    assert hash(rxn) == hash(ReactionContainer())
