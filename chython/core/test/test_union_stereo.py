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
"""`union` carries the SECOND molecule's stereo, which for a while it did not.

`union` is `copy()` of the left side plus a rebuild of the right one through `add_atom`/`add_bond`,
so everything the left side had came along for free and everything the right side had was silently
dropped: `reactant1.union(reactant2)` racemised reactant 2.  A reaction driver unions its inputs
before it patches them, so that loss was one call away from every stereospecific reaction.

The assertions compare CANONICAL BYTES of the split components against the originals rather than
SMILES strings, per the design's N9: a formatted string is not the identity of a stereo-bearing
molecule.
"""
import pytest

from chython.core import MoleculeContainer
from chython.core._core import read_smiles as smiles


def _round_trip(left, right):
    """Union then split, and the two pieces have to be the two molecules that went in."""
    a, b = smiles(left), smiles(right)
    parts = a.union(b).split()
    assert len(parts) == 2
    return {p.canonical_bytes for p in parts} == {a.canonical_bytes, b.canonical_bytes}


@pytest.mark.parametrize('left,right', [
    ('C[C@H](N)O', 'C[C@@H](N)Cl'),          # a tetrahedral centre on each side
    ('C[C@@H](N)O', 'C[C@H](N)Cl'),          # and both signs the other way round
    ('CCO', 'C[C@H](N)Cl'),                  # only the RIGHT side is configured: the regression
    ('C[C@H](N)Cl', 'CCO'),                  # only the left, which always worked
    ('F/C=C/F', 'C[C@H](N)Cl'),              # a cis/trans axis meeting a centre
    ('C[C@H](N)O', 'F/C=C\\F'),
    ('F/C=C/F', 'F/C=C\\F'),                 # two axes, opposite configurations
])
def test_both_sides_come_through_a_union_and_split_unchanged(left, right):
    assert _round_trip(left, right)


def test_the_right_sides_parity_is_carried_not_recomputed():
    """The loudest form of the regression: before the fix this atom read 0."""
    b = smiles('C[C@H](N)Cl')
    u = smiles('CCO').union(b)
    # the right side's atoms are appended, so its centre is atom 3 + 2 = 5
    assert u.parity_of(5) == b.parity_of(2)
    assert u.parity_of(5) != 0


def test_or_and_and_group_ids_are_renumbered_so_two_mixtures_stay_two():
    """THE ONE PIECE OF STEREO WHOSE MEANING IS NOT LOCAL TO AN ATOM.  Both sides number from 1, so
    carrying the ids verbatim would declare one mixture where the input declared two."""
    a = smiles('C[C@H](O)CC |&1:1|')
    b = smiles('C[C@H](O)CC |&1:1|')
    groups = a.union(b).stereo_groups()
    assert len(groups) == 2
    assert {k for k, _ in groups} == {3}                      # both AND
    assert sorted(g for _, g in groups) == [1, 2]
    assert sorted(len(v) for v in groups.values()) == [1, 1]


def test_an_or_group_and_an_and_group_are_numbered_independently():
    a = smiles('C[C@H](O)CC |&1:1|')
    b = smiles('C[C@H](O)CC |o1:1|')
    groups = a.union(b).stereo_groups()
    assert sorted(groups) == [(2, 1), (3, 1)]                 # OR 1 and AND 1 do not collide


def test_an_abs_group_needs_no_renumbering_and_keeps_its_zero():
    a = smiles('C[C@H](O)CC')
    b = smiles('C[C@H](O)CC |a:1|')
    assert a.union(b).stereo_groups() == {(1, 0): [7]}


def test_sixty_four_groups_of_one_kind_is_refused_rather_than_merged():
    def chain_of_groups(count):
        m = MoleculeContainer()
        with m.edit():
            for _ in range(count):
                m.add_atom(6, implicit_h=0)
        with m.edit():
            for i in range(1, count + 1):
                m.set_stereo_group(i, 3, i)
        return m

    a, b = chain_of_groups(32), chain_of_groups(32)
    with pytest.raises(ValueError, match='63 AND stereo groups'):
        a.union(b)
    # and 63 all told fits, which is what makes the refusal a boundary rather than a policy
    assert len(chain_of_groups(32).union(chain_of_groups(31)).stereo_groups()) == 63


def test_coordinates_come_from_either_side():
    a = MoleculeContainer()
    with a.edit():
        a.add_atom(6, implicit_h=4)
        a.set_xy(1, 1.5, -2.5)
    b = MoleculeContainer()
    with b.edit():
        b.add_atom(8, implicit_h=2)
        b.set_xy(1, 3.25, 4.0)
    u = a.union(b)
    assert u.has_coordinates
    assert u.xy_of(1) == (1.5, -2.5)
    assert u.xy_of(2) == (3.25, 4.0)
    # and the right side alone is enough: a flat left side must not suppress the coordinates it has
    plain = MoleculeContainer()
    with plain.edit():
        plain.add_atom(6, implicit_h=4)
    v = plain.union(b)
    assert v.has_coordinates
    assert v.xy_of(2) == (3.25, 4.0)


def test_a_union_of_two_flat_molecules_has_no_coordinates_and_no_groups():
    u = smiles('CCO').union(smiles('CCN'))
    assert not u.has_coordinates
    assert not u.has_stereo_groups
    assert u.wedges() == []


def test_wedges_come_from_both_sides_with_their_direction_intact():
    a = smiles('CC(N)O')
    b = smiles('CC(N)Cl')
    with a.edit():
        a.set_wedge(2, 3, 1)
    with b.edit():
        b.set_wedge(2, 4, 2)
    u = a.union(b)
    assert u.wedge_of(2, 3) == 1
    assert u.wedge_of(6, 8) == 2                              # b's atom 2 -> 6, its atom 4 -> 8
    assert u.wedge_of(8, 6) == 0                              # a wedge is narrow -> wide, one way


def test_remap_false_keeps_the_stereo_it_puts_the_numbers_back_on():
    a = smiles('CCO')
    b = smiles('C[C@H](N)Cl')
    with b.edit():                                            # move b clear of a's 1..3
        b.remap({1: 11, 2: 12, 3: 13, 4: 14})
    u = a.union(b, remap=False)
    assert u.parity_of(12) == b.parity_of(12) != 0
    assert {p.canonical_bytes for p in u.split()} == {a.canonical_bytes, b.canonical_bytes}


def test_a_molecule_unioned_with_itself_keeps_both_copies_configured():
    a = smiles('C[C@H](N)O')
    u = a.union(a)
    parts = u.split()
    assert len(parts) == 2
    assert {p.canonical_bytes for p in parts} == {a.canonical_bytes}
    assert u.parity_of(2) == u.parity_of(6) == a.parity_of(2)
