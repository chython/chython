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
"""`union` carries the SECOND molecule's stereo, which is the half that has to be carried.

`union` is `copy()` of the left side plus a rebuild of the right one through `add_atom`/`add_bond`, so
the left side's stereo arrives with its buffer and the right side's has to be written across
explicitly.  A gap there racemises reactant 2 of every `reactant1.union(reactant2)`, and a reaction
driver unions its inputs before it patches them, which puts that loss one call from every
stereospecific reaction.

The assertions compare CANONICAL BYTES of the split components against the originals rather than
SMILES strings, per the design's N9: a formatted string is not the identity of a stereo-bearing
molecule.
"""
import pytest

from chython.core import MoleculeContainer
from chython.core._core import read_smiles as smiles, read_smarts as smarts


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


def test_a_centre_and_an_axis_of_the_same_number_are_one_collection():
    """OR 1 is OR 1: one id namespace, so a centre and an axis given the same number are one group.

    The two members are spelled differently -- a centre is its own atom, an axis is its owner pair --
    and both live in one dict under one key.  A query that states no geometry is unaffected by the
    collection either way, which is why the merge is a storage change and not a matching change.
    """
    # `C[C@H](O)/C=C/C`: atom 2 is the tetrahedral centre, atoms 4-5 are the double bond.
    query = smarts('[C;z2]=[C;z2]')
    mol = smiles('C[C@H](O)/C=C/C')
    with mol.edit() as e:
        e.set_stereo_group(2, 2, 1)          # OR 1 on the tetrahedral centre
        e.set_bond_stereo_group(4, 5, 2, 1)  # OR 1 on the double bond
    assert mol.stereo_groups() == {(2, 1): [2, (4, 5)]}
    assert mol.bond_stereo_groups() == {(2, 1): [(4, 5)]}
    assert list(query.get_mapping(mol))


def test_allene_group_is_the_byte_at_the_units_anchor():
    """SU_ALLENE has no partner bond -- `chiral_atoms` already yields it beside SU_TETRA.

    Its anchor is the chain midpoint, so a group named on that atom is the group of the axis, and
    reads back as the terminal pair the axis is named on.
    """
    # `CC=C=CC`: the chain centre (atom 3) is the allene anchor; `chiral_atoms()` keys on it.
    mol = smiles('CC=C=CC')
    centre = next(n for n in mol.chiral_atoms())
    with mol.edit() as e:
        e.set_stereo_group(centre, 2, 1)
    assert mol.stereo_group_of(centre) == (2, 1)
    assert mol.bond_stereo_groups() == {(2, 1): [(2, 4)]}


def test_cis_trans_unit_with_an_or_group_on_its_anchor_is_a_mixture():
    """A cis/trans unit's group is the byte at its anchor, whatever spelling stated it.

    CXSMILES `|o1:|` and V3000 `STEREL ATOMS=` state the OR collection on atoms, and an atom that
    anchors the unit is exactly where the byte belongs, so `_sg_byte_raw` reads it there.
    Both query geometries must match because the OR group makes the molecule a mixture.
    """
    mol = smiles('C/C=C/C')
    with mol.edit() as e:
        e.set_stereo_group(2, 2, 1)   # atom OR 1 on the double-bond anchor
    assert smarts('C/C=C/C').is_substructure(mol)    # same geometry matches
    assert smarts('C/C=C\\C').is_substructure(mol)   # opposite geometry also matches


def test_bond_only_or_group_gates_stereo_groups_admit():
    """A molecule whose only OR groups are stated on AXES must still run `stereo_groups_admit`.

    `has_or_group` scans the one group segment, and an axis's byte sits at its unit's anchor like any
    other, so the scan finds it: that is the gate deciding whether `stereo_groups_admit` is consulted
    at all.  A gate left False on this molecule wrongly accepts a query that can only satisfy one
    member of a multi-member OR group.

    The molecule is `F/C=C/C=C/F` (E,E diene) with both double bonds in OR group 1, stated through
    the pair spelling; the mixture is {E,E} and {Z,Z}.  A query demanding E on one bond and Z on the
    other has no consistent assignment and must be refused.
    """
    mol = smiles('F/C=C/C=C/F')    # E,E diene; double bonds 2-3 and 4-5
    with mol.edit() as e:
        e.set_bond_stereo_group(2, 3, 2, 1)   # OR 1 on the first axis
        e.set_bond_stereo_group(4, 5, 2, 1)   # OR 1 on the second axis (same group)
    # A consistent assignment (both E) must succeed.
    assert smarts('F/C=C/C=C/F').is_substructure(mol)
    # No consistent assignment for E on bond 1 and Z on bond 2: must be refused.
    assert not smarts('F/C=C/C=C\\F').is_substructure(mol)


def test_every_spelling_of_one_cumulene_axis_reads_the_same_group():
    """ONE AXIS, ONE BYTE, at the unit's anchor -- reached by every honest way of naming the axis.

    `F/C(=C=C=C/C)Cl` has a stereogenic SU_CIS_TRANS unit whose owners are the chain terminals C2 and
    C5, three bonds apart: atoms F=1, C=2, C=3, C=4, C=5, C=6, Cl=7 with chain double bonds 2-3, 3-4
    and 4-5.  The owner pair, any chain bond of it, and the anchor atom alone all name that unit, so
    all four spellings store one group at one slot and the molecule becomes a mixture either way.

    Control first: without a group the unit is configured and the opposite geometry is refused.
    """
    assert not smarts('F/C(=C=C=C\\C)Cl').is_substructure(smiles('F/C(=C=C=C/C)Cl'))
    for spelling in ((2, 5), (2, 3), (3, 4), (4, 5)):
        mol = smiles('F/C(=C=C=C/C)Cl')
        with mol.edit() as e:
            e.set_bond_stereo_group(spelling[0], spelling[1], 2, 1)
        assert mol.stereo_groups() == {(2, 1): [(2, 5)]}, spelling
        assert smarts('F/C(=C=C=C\\C)Cl').is_substructure(mol), spelling
    mol = smiles('F/C(=C=C=C/C)Cl')
    with mol.edit() as e:
        e.set_stereo_group(2, 2, 1)               # the anchor atom, a bare int
    assert mol.stereo_groups() == {(2, 1): [(2, 5)]}
    assert smarts('F/C(=C=C=C\\C)Cl').is_substructure(mol)


def test_canonical_bond_ids_correct_for_disconnected_molecule_n_gt_n_edges():
    """Canonical group ids are permutation-invariant when n > n_edges (disconnected molecule).

    `F/C=C/Cl.C/C=C/C` plus 15 isolated [NH4+] gives n=23 atoms and 12 half-edges.  A member is one
    group byte at one anchor slot, so `canonical_stereo_group_ids` writes at most n members into a
    region of 2n class slots; a shape with more grouped atoms than half-edges is the one where a
    region sized off the edge count would run into `par`, the parity-code array laid out after it.

    The assertion is PERMUTATION INVARIANCE -- the two axes' canonical ids must not depend on which
    stored id each carries.  It is the guarantee `test_canonical_bond_ids_are_permutation_invariant`
    states for a connected molecule, on the n > n_edges shape.  The axes' stored ids stay clear of the
    ammonium ions' OR 1: with one namespace an axis given id 1 would JOIN that collection, which is a
    different molecule rather than a permutation of this one.
    """
    base = 'F/C=C/Cl.C/C=C/C.' + '.'.join(['[NH4+]'] * 15)
    answers = set()
    for a, b in ((2, 3), (3, 2), (7, 19), (63, 4)):
        mol = smiles(base)
        with mol.edit() as e:
            for i in range(9, 24):            # atoms 9-23 are the 15 [NH4+] ions
                e.set_stereo_group(i, 2, 1)   # OR 1 on each
            e.set_bond_stereo_group(2, 3, 2, a)   # OR on the F/C=C/Cl axis
            e.set_bond_stereo_group(6, 7, 2, b)   # OR on the C/C=C/C axis
        answers.add((tuple(sorted((k, tuple(sorted(v))) for k, v in
                                  mol.canonical_bond_stereo_groups().items())),
                     tuple(sorted((k, tuple(sorted(map(str, v)))) for k, v in
                                  mol.canonical_stereo_groups().items()))))
    assert len(answers) == 1, answers
