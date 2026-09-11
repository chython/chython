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
"""A view answers for itself.

The rule these tests pin: a caller holding an `Atom` or a `Bond` never goes back to the container to
ask about that atom or that bond.  Before them a renderer wrote `element_symbols()[a.element]` and
`mol.parity_of(a.n)` while already holding `a`, which is data in one object and every question
about it answered somewhere else.  The bulk readers below exist for the same reason one step up: a
caller that wants the whole plane, or the aromatic rings, should not assemble either.
"""
from pytest import approx, raises
# `chython.smiles` is the bidirectional door over this function, and nothing under `chython/core/` may
# import the facade (`test_no_chython_two_imports.py` ratchets it).  Aliased rather than spelled out at
# each call so the bodies below read the way a caller writes them.
from chython.core import STEREO_UNSPECIFIED, WEDGE_UP, read_smiles as smiles


def test_atomic_symbol():
    mol = smiles('CCl')
    assert [a.atomic_symbol for a in mol.atoms()] == ['C', 'Cl']


def test_atomic_symbol_agrees_with_the_table_for_every_element():
    """one symbol table in the library, and the 0-vs-1 index base is handled in exactly one place"""
    from chython.core._core import element_symbols

    table = element_symbols()
    mol = smiles('C')
    sid = next(iter(mol))
    for z in range(1, 119):
        with mol.edit():
            mol.set_element(sid, z)
        assert mol.atom(sid).atomic_symbol == table[z], z


def test_atomic_radius_agrees_with_the_table_for_every_element():
    """Element data on the view, so a renderer holding an `Atom` does not go looking for a table."""
    from chython.core._core import atomic_radius_table

    table = atomic_radius_table()
    mol = smiles('C')
    sid = next(iter(mol))
    for z in range(1, 119):
        with mol.edit():
            mol.set_element(sid, z)
        assert mol.atom(sid).atomic_radius == table[z], z


def test_atomic_radius_matches_the_container():
    mol = smiles('CCl')
    for a in mol.atoms():
        assert a.atomic_radius == mol.radius_of(a.n)
    assert [a.atomic_radius for a in mol.atoms()] == [0.67, 0.79]


def test_the_r_marker_has_no_radius():
    """0.0, the answer it gives for mass as well: a marker is not an element and carries neither.

    Not None, because the reader of this is a renderer sizing a sphere -- a marker gets no sphere, and
    an arithmetic zero says so where None would raise inside the drawing loop.
    """
    mol = smiles('[R]C')
    assert [a.atomic_radius for a in mol.atoms()] == [0.0, 0.67]


def test_repr_reads_as_chemistry():
    mol = smiles('c1ccccc1')
    sid = next(iter(mol))
    assert repr(mol.atom(sid)) == f'Atom(C, n={sid})'


def test_parity_matches_the_container():
    mol = smiles('C[C@H](N)O')
    for a in mol.atoms():
        assert a.parity == mol.parity_of(a.n)
    assert any(a.parity for a in mol.atoms()), 'no parity to compare'


def test_stereo_group_matches_the_container():
    mol = smiles('C[C@H](N)O')
    sid = next(iter(mol))
    assert mol.atom(sid).stereo_group == (STEREO_UNSPECIFIED, 0)
    with mol.edit():
        mol.set_stereo_group(sid, 2, 1)
    assert mol.atom(sid).stereo_group == (2, 1)


def test_atom_cip_matches_the_container():
    mol = smiles('C[C@H](N)O')
    sid = [a.n for a in mol.atoms() if a.atomic_symbol == 'C'][1]
    assert mol.atom(sid).cip is None
    with mol.edit():
        mol.set_atom_cip(sid, 'R')
    assert mol.atom(sid).cip == 'R'


def test_xy_still_answers_none_without_a_plane():
    mol = smiles('CC')
    a = next(iter(mol.atoms()))
    assert a.x is None and a.y is None


def test_xy_reads_what_set_xy_wrote():
    mol = smiles('CC')
    sid = next(iter(mol))
    with mol.edit():
        mol.set_xy(sid, 1.25, -2.5)
    a = mol.atom(sid)
    assert (a.x, a.y) == approx((1.25, -2.5))
    assert a.xy == approx((1.25, -2.5))


def test_a_stale_view_still_refuses():
    """the speed-ups must not lose the generation guard; a delegating property must not lose it either"""
    mol = smiles('C[C@H](N)O')
    a = next(iter(mol.atoms()))
    with mol.edit():
        mol.delete_atom([sid for sid in mol if sid != a.n][0])
    with raises(RuntimeError):
        a.x
    with raises(RuntimeError):
        a.atomic_symbol
    with raises(RuntimeError):
        a.parity
    with raises(RuntimeError):
        a.stereo
    with raises(RuntimeError):
        a.stereo_group
    with raises(RuntimeError):
        a.xy


def test_a_stale_bond_view_still_refuses():
    """Bond.wedge and Bond.cip both delegate to _ptr(), which guards the generation"""
    mol = smiles('CCO')
    b = next(iter(mol.bonds()))
    with mol.edit():
        mol.delete_atom(b.n)           # deletes one endpoint, invalidating the bond
    with raises(RuntimeError):
        b.wedge
    with raises(RuntimeError):
        b.cip


def test_bond_wedge_is_none_when_nothing_was_drawn():
    mol = smiles('C[C@H](N)O')
    assert all(b.wedge is None for b in mol.bonds())


def test_bond_wedge_names_the_narrow_end_from_either_side():
    """a wedge is DIRECTIONAL and a Bond is not, so the answer has to be an atom

    This is the try-then-swap dance -- look the pair up, look it up reversed -- moved into the core
    once.  A caller that gets a bare code back cannot tell which way the wedge points and has to ask
    twice, which is exactly what three call sites did.
    """
    mol = smiles('C[C@H](N)O')
    a, b = [sid for sid in mol][:2]
    with mol.edit():
        mol.set_wedge(b, a, WEDGE_UP)         # narrow end at b

    assert mol.wedge_between(a, b) == (b, WEDGE_UP)
    assert mol.wedge_between(b, a) == (b, WEDGE_UP), 'the answer must not depend on the argument order'
    assert mol.bond(a, b).wedge == (b, WEDGE_UP)
    assert mol.bond(b, a).wedge == (b, WEDGE_UP)


def test_wedge_between_refuses_a_pair_that_is_not_a_bond():
    mol = smiles('CCO')
    a, _, c = list(mol)
    with raises(KeyError):
        mol.wedge_between(a, c)


def test_bond_cip_matches_the_container():
    mol = smiles('C/C=C/C')
    n, m = [(b.n, b.m) for b in mol.bonds() if b.order == 2][0]
    assert mol.bond(n, m).cip is None
    with mol.edit():
        mol.set_bond_cip(n, m, 'E')
    assert mol.bond(n, m).cip == 'E'
    assert mol.bond(m, n).cip == 'E', 'a bond CIP is not directional'


def test_aromatic_rings_is_a_filter_over_rings():
    mol = smiles('c1ccc2ccccc2c1')          # naphthalene: two aromatic rings
    assert len(mol.aromatic_rings) == 2
    assert all(ring in mol.rings for ring in mol.aromatic_rings)

    mol = smiles('c1ccccc1C1CCCCC1')        # one aromatic, one saturated
    assert len(mol.rings) == 2
    assert len(mol.aromatic_rings) == 1
    assert all(mol.order_of(ring[i - 1], ring[i]) == 4
               for ring in mol.aromatic_rings for i in range(len(ring)))


def test_aromatic_rings_is_empty_after_kekule():
    mol = smiles('c1ccccc1')
    assert len(mol.aromatic_rings) == 1
    mol.kekule()
    assert mol.aromatic_rings == [], 'a kekulized molecule has no order-4 bonds; say so'


def test_has_layout_separates_a_plane_from_a_segment():
    mol = smiles('CCO')
    assert not mol.has_coordinates
    assert not mol.has_layout

    with mol.edit():                       # a segment, every atom at the origin
        for sid in mol:
            mol.set_xy(sid, 0., 0.)
    assert mol.has_coordinates, 'expected a molecule with an XY segment'
    assert not mol.has_layout, 'a degenerate plane is not a layout'

    with mol.edit():
        for i, sid in enumerate(mol):
            mol.set_xy(sid, .825 * i, 0.)
    assert mol.has_layout, 'span in one axis is a layout; a linear molecule is not degenerate'


def test_one_atom_needs_no_layout():
    mol = smiles('C')
    with mol.edit():
        mol.set_xy(next(iter(mol)), 0., 0.)
    assert mol.has_layout


def test_has_layout_empty_molecule_is_false():
    """no atoms means no layout even when the segment housekeeping would allow one"""
    mol = smiles('C')
    sid = next(iter(mol))
    with mol.edit():
        mol.set_xy(sid, 1., 1.)
    assert mol.has_layout, 'pre-condition: one atom with a segment is True'
    with mol.edit():
        mol.delete_atom(sid)           # deletes the only atom; the XY segment disappears with it
    assert not mol.has_layout, 'empty molecule has no layout'


def test_has_layout_straddles_the_threshold():
    """span = 99 fixed-point units is below the threshold; span = 100 is at or above"""
    mol = smiles('CC')
    a, b = list(mol)
    # span 99: 0.0099 molecule units apart -- below 0.01 (= 100 in fixed-point)
    with mol.edit():
        mol.set_xy(a, 0., 0.)
        mol.set_xy(b, 0.0099, 0.)
    assert not mol.has_layout, 'span 99 fixed-point is below the 100-unit threshold'
    # span 100: exactly 0.01 molecule units -- at the threshold
    with mol.edit():
        mol.set_xy(b, 0.01, 0.)
    assert mol.has_layout, 'span 100 fixed-point meets the threshold'


def test_has_layout_span_does_not_overflow_int32():
    """two atoms at opposite ends of the int32 range differ by > INT32_MAX -- the span check must
    use int64_t arithmetic or the subtraction overflows to a negative value and has_layout lies.

    Verified to FAIL against the int32_t arithmetic by checking the condition max_x - min_x >= 100
    with int32_t arithmetic: round(214748 * 10000) = 2147480000, and
    2147480000 - (-2147480000) overflows int32_t, producing a negative result < 100.
    """
    mol = smiles('CC')
    a, b = list(mol)
    with mol.edit():
        mol.set_xy(a,  214748., 0.)    # near-max positive coordinate
        mol.set_xy(b, -214748., 0.)    # near-max negative coordinate
    assert mol.has_layout, 'span far exceeds the threshold; wrong only if int32 overflow occurred'


def test_coordinates_reads_the_whole_plane_once():
    mol = smiles('CCO')
    assert mol.coordinates() == {}, 'no plane stated, so no plane reported -- not a plane of zeros'

    sids = list(mol)
    with mol.edit():
        for i, sid in enumerate(sids):
            mol.set_xy(sid, float(i), -float(i))
    plane = mol.coordinates()
    # Assert the literal (i, -i) values the loop wrote -- not against Atom.x/y (same code path
    # after the consolidation) so the test has an independent witness.
    for i, sid in enumerate(sids):
        assert plane[sid] == approx((float(i), -float(i))), f'atom {i}: expected ({i}, {-i})'
    # Agreement assertion kept as a cross-check, but the literals above are what makes it independent.
    assert plane == {sid: approx((mol.atom(sid).x, mol.atom(sid).y)) for sid in mol}
    assert list(plane) == list(mol), 'arena order, so a caller can zip it against iteration'


def test_aromatic_rings_checks_the_bond_that_CLOSES_the_ring():
    """the walk wraps, and the closing bond is the one an off-by-one silently drops

    `rings` yields a cycle as a tuple, so one of its bonds -- `(ring[-1], ring[0])` -- is not between
    two adjacent entries.  A filter that walked only the adjacent pairs would call this ring aromatic.
    Every other ring in this file has all of its bonds aromatic or none of them and would pass either
    way, so benzene with just its closing bond set to order 1 is what tells the two implementations
    apart.  Measured rather than argued: under `prev = ring[0]` with `range(1, size)` -- the walk that
    visits the adjacent pairs and nothing else -- this is the ONLY test of 2751 that fails.
    """
    mol = smiles('c1ccccc1')
    ring = mol.rings[0]
    with mol.edit():
        mol.set_order(ring[-1], ring[0], 1)
    assert mol.rings == [ring], 'the ring must survive the edit, or this tests nothing'
    assert mol.aromatic_rings == [], 'one non-aromatic bond is enough, wherever in the ring it sits'
