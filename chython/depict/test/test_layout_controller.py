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
"""The layout controller: score, tree embedding, edits, search."""
from math import atan2, cos, degrees, hypot, isfinite, pi, radians, sin

import pytest

from chython.core import MoleculeContainer, smiles
from chython.depict.layout import search as search_mod
from chython.depict.layout.edits import TARGETS, VARIANTS, Move, apply, candidates, embed, masked, targets, variants
from chython.depict.layout.search import controller_layout
from chython.depict.layout.lattice import L
from chython.depict.layout.quality import defects, perm, rescaled, shape
from chython.depict.layout.rings import ring_blocks
from chython.depict.layout.score import Scope, add, better, delta, frontier, full, key, open_closures
from chython.depict.layout.search import controller_layout
from chython.depict.layout.tree import tree_embedding


def _chain(n, turn=120):
    """`n` atoms of an alkane zigzag at `L`, as a plane keyed 1..n."""
    out, x, y, a = {}, 0., 0., radians(30)
    for i in range(1, n + 1):
        out[i] = (x, y)
        a = a - radians(180 - turn) if i % 2 else a + radians(180 - turn)
        x, y = x + L * cos(a), y + L * sin(a)
    return out


def test_a_zigzag_scores_zero():
    mol = smiles('CCCCCC')
    s = Scope(mol, set(mol), mol.canonical_order(), [])
    assert full(s, _chain(6)) == (0, 0., 0.)


def test_a_closure_off_L_costs_its_error_and_is_open():
    mol = smiles('C1CCCCC1')
    s = Scope(mol, set(mol), mol.canonical_order(), [(1, 6)])
    plane = _chain(6)
    err = abs(hypot(plane[1][0] - plane[6][0], plane[1][1] - plane[6][1]) - L)
    assert full(s, plane)[1] == pytest.approx(err)
    assert open_closures(s, plane) == [(1, 6)]


def test_delta_equals_the_difference_of_full_scores():
    mol = smiles('CCCCCC')
    s = Scope(mol, set(mol), mol.canonical_order(), [])
    plane = _chain(6)
    moved = {6: (plane[4][0] + .1, plane[4][1])}               # atom 6 dropped next to atom 4
    after = full(s, {**plane, **moved})
    assert add(full(s, plane), delta(s, plane, moved)) == pytest.approx(after)


def test_defects_outrank_closure_and_closure_outranks_conventions():
    assert better((0, 9., 9.), (1, 0., 0.))
    assert better((0, 0., 9.), (0, .5, 0.))
    assert not better(key((0, 0., 1.)), key((0, 0., 1. + 1e-12)))


def test_frontier_is_the_trouble_and_two_bonds_around_it():
    mol = smiles('CCCCCCCCCC')
    s = Scope(mol, set(mol), mol.canonical_order(), [])
    plane = _chain(10)
    assert frontier(s, plane) == set()
    plane[10] = plane[8]                                    # 10 folds back onto 8, and crowds 7
    assert frontier(s, plane) == {5, 6, 7, 8, 9, 10}        # 4 is 3 bonds from 7, the nearest trouble


# --- tree embedding --------------------------------------------------------------------------------------

def _tree(text):
    mol = smiles(text)
    return mol, tree_embedding(Scope(mol, set(mol), mol.canonical_order(), []))


def _angle(plane, x, a, b):
    u = atan2(plane[a][1] - plane[x][1], plane[a][0] - plane[x][0])
    v = atan2(plane[b][1] - plane[x][1], plane[b][0] - plane[x][0])
    d = abs(degrees(u - v)) % 360
    return min(d, 360 - d)


def test_every_tree_bond_is_L_and_a_ring_leaves_one_closure():
    mol, t = _tree('c1ccccc1CCO')
    for x, p in t.parent.items():
        if p is not None:
            assert hypot(t.plane[x][0] - t.plane[p][0], t.plane[x][1] - t.plane[p][1]) == pytest.approx(L)
    assert len(t.closures) == 1


def test_a_chain_zigzags_and_scores_zero():
    mol, t = _tree('CCCCCCCC')
    assert full(Scope(mol, set(mol), mol.canonical_order(), t.closures), t.plane) == (0, 0., 0.)


def test_a_triple_bond_is_straight_and_a_branch_is_120():
    mol, t = _tree('CC#CC(C)C')
    assert _angle(t.plane, 2, 1, 3) == pytest.approx(180)
    assert _angle(t.plane, 4, 5, 6) == pytest.approx(120)


def test_a_quaternary_carbon_is_a_cross():
    mol, t = _tree('CC(C)(C)C')
    assert sorted(round(_angle(t.plane, 2, 1, y)) for y in (3, 4, 5)) == [90, 90, 180]


def test_a_gem_dimethyl_ring_atom_keeps_its_ring_angle():
    mol, t = _tree('CC1(C)CCCCC1')
    assert _angle(t.plane, 2, 4, 8) == pytest.approx(120)
    assert _angle(t.plane, 2, 1, 3) == pytest.approx(80)


@pytest.mark.parametrize('text', ['C/C=C/C', 'C/C=C\\C', 'CC/C=C(/C)CC', 'OC/C=C\\CCC'])
def test_a_stated_cis_trans_is_kept(text):
    mol, t = _tree(text)
    assert shape(mol, t.plane)['ct'] == 0


@pytest.mark.parametrize('text', ['c1ccc2ccccc2c1', 'c1ccc2cc3ccccc3cc2c1', 'c1ccc2c(c1)ccc1ccccc12',
                                  'c1cc2ccc3cccc4ccc(c1)c2c34'])
def test_fused_hexagons_tile_the_honeycomb_closed(text):
    mol, t = _tree(text)
    scope = Scope(mol, set(mol), mol.canonical_order(), t.closures)
    defects, closure, _ = full(scope, t.plane)
    assert defects == 0 and closure == pytest.approx(0, abs=1e-6)


@pytest.mark.parametrize('text', ['c1ccc2[nH]ccc2c1', 'C1CCC2C(C1)CCC1CCCC12', pytest.param(
    'CC12CCC3C(CCC4=CC(=O)CCC34C)C1CCC2O', marks=pytest.mark.xfail(reason='one-level lookahead: a steroid starts with '
                                                                       'defects the policy must repair')), pytest.param(
    'CN1CCC23C4C1CC5=C2C(=C(C=C5)O)OC3C(C=C4)O', marks=pytest.mark.xfail(reason='one-level lookahead'))])
def test_a_fused_system_starts_without_a_defect(text):
    mol, t = _tree(text)
    scope = Scope(mol, set(mol), mol.canonical_order(), t.closures)
    assert full(scope, t.plane)[0] == 0


def test_a_long_chain_is_no_recursion():
    mol, t = _tree('C' * 400)
    assert len(t.plane) == 400


def test_a_metal_of_six_neighbours_is_60_apart():
    mol, t = _tree('N#C[Fe](C#N)(C#N)(C#N)(C#N)C#N')
    fe = next(x for x in mol if mol.atom(x).atomic_symbol == 'Fe')
    nb = sorted(mol.neighbors_of(fe))
    assert min(_angle(t.plane, fe, a, b) for a in nb for b in nb if a < b) == pytest.approx(60)


# --- edits -----------------------------------------------------------------------------------------------

def _drawn(text):
    mol = smiles(text)
    return mol, Scope(mol, set(mol), mol.canonical_order(), []), controller_layout(mol)


def _dir(plane, x, y):
    return atan2(plane[y][1] - plane[x][1], plane[y][0] - plane[x][0]) % (2 * pi)


def _near(a, b, tol=radians(1)):
    d = (a - b) % (2 * pi)
    return min(d, 2 * pi - d) < tol


def test_an_angular_methyl_targets_the_line_of_the_fused_bond():
    mol, s, plane = _drawn('CC12CCCCC1CCCC2')           # 1 the methyl on fusion atom 2; 2-7 the fused bond
    opposite = (_dir(plane, 2, 7) + pi) % (2 * pi)
    off = opposite + radians(30)                        # turned off it, halfway to a ring bond
    plane[1] = (plane[2][0] + L * cos(off), plane[2][1] + L * sin(off))
    assert any(_near(t, opposite) for t in targets(s, plane, 2, 1))


def test_a_gem_dimethyl_targets_the_thirds_of_its_outside_gap():
    mol, s, plane = _drawn('CC1(C)CCCCC1')
    a, b = _dir(plane, 2, 4), _dir(plane, 2, 8)
    start, width = (a, (b - a) % (2 * pi)) if (b - a) % (2 * pi) > pi else (b, (a - b) % (2 * pi))
    thirds = [start + width / 3, start + 2 * width / 3]
    found = targets(s, plane, 2, 1) + [_dir(plane, 2, 1), _dir(plane, 2, 3)]
    assert all(any(_near(t, h) for t in found) for h in thirds)


def test_at_most_four_targets_none_on_a_bond():
    mol, s, plane = _drawn('CC(C)(C)CC(C)(C)C')
    for x in mol:
        for y in mol.neighbors_of(x):
            ts = targets(s, plane, x, y)
            assert len(ts) <= TARGETS
            for t in ts:
                assert all(not _near(t, _dir(plane, x, z), radians(29)) for z in mol.neighbors_of(x) if z != y)


def test_a_clean_plane_has_no_candidates_and_a_folded_one_has_some():
    mol, s, plane = _drawn('CCCCCCCC')
    assert candidates(s, plane, frontier(s, plane)) == []
    plane[8] = plane[6]
    moves = candidates(s, plane, frontier(s, plane))
    assert {m.kind for m in moves} == {'flip', 'rotate', 'stretch'}
    assert all(m.part and m.pivot not in m.part and m.other in m.part for m in moves)


def test_a_stretch_reaches_its_length():
    mol, s, plane = _drawn('CCCC')
    m = Move('stretch', 2, 1, 1.5 * L, frozenset({1}))
    p = {**plane, **apply(s, plane, m)}
    assert hypot(p[1][0] - p[2][0], p[1][1] - p[2][1]) == pytest.approx(1.5 * L)


def test_bending_a_linear_atom_is_masked():
    mol, s, plane = _drawn('CC#CCC')
    m = Move('rotate', 3, 2, pi / 3, frozenset({1, 2}))
    assert masked(s, plane, apply(s, plane, m))


def test_breaking_a_stated_trans_is_masked():
    mol, s, plane = _drawn('C/C=C/C')
    m = Move('rotate', 3, 4, 2 * pi / 3, frozenset({4}))
    moved = apply(s, plane, m)
    assert masked(s, plane, moved)


def test_bends_only_when_violating():
    mol, s, plane = _drawn('CCCCCCCC')
    plane[8] = plane[6]
    front = frontier(s, plane)
    assert all(m.kind != 'bend' for m in candidates(s, plane, front))
    assert {m.kind for m in candidates(s, plane, front, violate=True)} == {'bend'}


# --- ring systems ----------------------------------------------------------------------------------------

def _embedded(text):
    """`mol`, its scope and the tree plane with every ring system drawn as its first variant."""
    mol = smiles(text)
    rank = mol.canonical_order()
    tree = tree_embedding(Scope(mol, set(mol), rank, []))
    s = Scope(mol, set(mol), rank, tree.closures)
    plane = dict(tree.plane)
    for block in ring_blocks(mol):
        plane.update(embed(s, plane, block[0], variants(s, block)[0], tree.root))
    return mol, s, plane


def test_variants_are_few_distinct_and_at_L():
    mol = smiles('c1ccc2ccccc2c1')
    s = Scope(mol, set(mol), mol.canonical_order(), [])
    vs = variants(s, ring_blocks(mol)[0])
    assert 1 <= len(vs) <= VARIANTS
    for v in vs:
        for b in mol.bonds():
            assert hypot(v[b.n][0] - v[b.m][0], v[b.n][1] - v[b.m][1]) == pytest.approx(L, rel=.02)


def test_an_embedded_ring_closes():
    mol, s, plane = _embedded('c1ccccc1')
    assert open_closures(s, plane) == []
    assert full(s, plane)[0] == 0


def test_a_substituent_leaves_on_the_ring_bisector_rigidly():
    mol, s, plane = _embedded('CCCCc1ccccc1')           # butyl 1-4 on ring atom 5, ring neighbours 6 and 10
    tree = tree_embedding(Scope(mol, set(mol), mol.canonical_order(), []))
    assert _angle(plane, 5, 4, 6) == pytest.approx(120, abs=1)
    assert _angle(plane, 5, 4, 10) == pytest.approx(120, abs=1)
    for a, b in ((1, 3), (1, 4), (2, 4)):               # the butyl's own shape is untouched
        assert hypot(plane[a][0] - plane[b][0], plane[a][1] - plane[b][1]) == \
            pytest.approx(hypot(tree.plane[a][0] - tree.plane[b][0], tree.plane[a][1] - tree.plane[b][1]))


def test_an_angular_methyl_leaves_along_the_fused_bond():
    mol, s, plane = _embedded('CC12CCCCC1CCCC2')        # cholesterol's C10 methyl, on a decalin
    assert _angle(plane, 2, 1, 7) == pytest.approx(180, abs=1)


def test_cholesterol_embeds_closed():
    log = []                                            # its pentagon's lattice chord is regular after the post-steps
    controller_layout(smiles('CC(C)CCCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C'), log)
    assert not [r for r in log if r.rule == 'depict:open-ring']


def test_two_ring_systems_on_a_chain_both_close():
    mol, s, plane = _embedded('c1ccccc1CCCC1CCCCC1')
    assert open_closures(s, plane) == []


# --- search ----------------------------------------------------------------------------------------------

def _clean(mol, plane):
    return all(v == 0 for k, v in defects(mol, rescaled(mol, plane)).items() if k != 'ovl61')


def _whole(mol, plane):
    return set(plane) == set(mol) and all(isfinite(v) for p in plane.values() for v in p)


@pytest.mark.parametrize('text', ['c1ccccc1', 'c1ccc2ccccc2c1', 'CCCCCCCC', 'CC(C)(C)C', 'CC12CCCCC1CCCC2',
                                  'CC(C)CCCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C', 'c1ccccc1-c1ccccc1',
                                  'C1CCC2(C1)CCCCC2', 'CC(=O)Oc1ccccc1C(=O)O'])
def test_public_compounds_draw_clean_closed_and_order_blind(text):
    mol = smiles(text)
    log = []
    plane = controller_layout(mol, log)
    assert _whole(mol, plane)
    assert _clean(mol, plane)
    assert not [r for r in log if r.rule == 'depict:open-ring']
    assert perm(mol, controller_layout, k=3) == 0


@pytest.mark.parametrize('text', ['C/C=C/C=C/CCC(/C)=C/C', 'OC/C=C\\CCC/C=C/C'])
def test_stated_cis_trans_survive_the_search(text):
    mol = smiles(text)
    assert shape(mol, controller_layout(mol))['ct'] == 0


def test_empty_single_atom_and_salt():
    assert controller_layout(MoleculeContainer()) == {}
    assert len(controller_layout(smiles('[Na+]'))) == 1
    mol = smiles('[Na+].[O-]C(=O)c1ccccc1')
    plane = controller_layout(mol)
    assert _whole(mol, plane)
    xs = sorted({round(plane[x][0], 3) for x in mol})
    assert len(xs) > 1


def test_a_metal_of_six_cyanides_is_finite():
    mol = smiles('N#C[Fe-4](C#N)(C#N)(C#N)(C#N)C#N')
    assert _whole(mol, controller_layout(mol))


def test_a_long_chain_ends_within_its_budget():
    mol = smiles('C' * 40)
    stats = []
    plane = controller_layout(mol, stats=stats)
    assert _whole(mol, plane) and stats[0]['moves'] <= 16


def test_explicit_hydrogens_are_terminal_atoms():
    mol = smiles('[H]C([H])([H])C([H])([H])O[H]')
    plane = controller_layout(mol)
    assert _whole(mol, plane) and _clean(mol, plane)


def test_a_ring_no_block_function_draws_stays_a_tree_and_is_logged(monkeypatch):
    monkeypatch.setattr(search_mod, 'variants', lambda scope, block: [])
    mol = smiles('C1CC2CCC1CC2')                         # bicyclo[2.2.2]octane: no tree closes a bridge
    log = []
    plane = controller_layout(mol, log)
    assert _whole(mol, plane)
    assert [r.rule for r in log] == ['depict:open-ring']


def test_the_search_is_deterministic():
    mol = smiles('CC(C)CCCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C')
    assert controller_layout(mol) == controller_layout(mol)


# --- review fixes ----------------------------------------------------------------------------------------

def _scoped_tree(text):
    mol = smiles(text)
    rank = mol.canonical_order()
    tree = tree_embedding(Scope(mol, set(mol), rank, []))
    return mol, Scope(mol, set(mol), rank, tree.closures), tree


@pytest.mark.parametrize('text', ['CC(C)(C)CC(C)(C)OP(=O)(OC)OC(C)(C)CC(C)(C)C', 'CCCCCCCCC(CCCC)(CCCC)CCCCCCCC'])
def test_every_successor_delta_matches_the_full_score(text):
    mol, s, tree = _scoped_tree(text)
    plane = tree.plane
    base = full(s, plane)
    moves = candidates(s, plane, frontier(s, plane)) + candidates(s, plane, frontier(s, plane), violate=True)
    assert moves
    for m in moves:
        moved = apply(s, plane, m)
        assert add(base, delta(s, plane, moved)) == pytest.approx(full(s, {**plane, **moved})), m


def test_every_embed_delta_matches_the_full_score():
    mol, s, tree = _scoped_tree('CC(C)CCCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C')
    plane = tree.plane
    base = full(s, plane)
    for block in ring_blocks(mol):
        for v in variants(s, block):
            moved = embed(s, plane, block[0], v, tree.root)
            assert add(base, delta(s, plane, moved)) == pytest.approx(full(s, {**plane, **moved}))


@pytest.mark.parametrize('text', ['CC/C=C1\\CCC[C@H]1C', 'OC/C=C1/CCOC1', 'CC/C=C1/CCCCC1C'])
def test_an_exocyclic_cis_trans_is_kept_by_the_ring_embedding(text):
    mol = smiles(text)                                  # the tree's ring path cannot hold it: the variant must
    assert shape(mol, controller_layout(mol))['ct'] == 0


@pytest.mark.parametrize('text', ['FS(F)(F)(F)(F)F', 'F[P-](F)(F)(F)(F)F', 'N#C[Fe-4](C#N)(C#N)(C#N)(C#N)C#N'])
def test_the_ligands_of_one_centre_do_not_crowd_each_other(text):
    mol = smiles(text)
    plane = controller_layout(mol)
    lengths = [hypot(plane[b.n][0] - plane[b.m][0], plane[b.n][1] - plane[b.m][1]) for b in mol.bonds()]
    assert max(lengths) == pytest.approx(min(lengths), rel=1e-6)


def test_a_macrocycle_exit_points_out_of_the_ring_at_a_reflex_corner():
    from chython.depict.layout.edits import outside_angles
    from chython.depict.layout.quality import _inside
    mol = smiles('CC1C(C)C(C)C(C)C(C)C(C)CCCCCCCCCCCCC1')
    s = Scope(mol, set(mol), mol.canonical_order(), [])
    (atoms, rings), = ring_blocks(mol)
    ring = max(rings, key=len)
    for v in variants(s, (atoms, rings)):
        poly = [v[x] for x in ring]
        for x in ring:
            if len(mol.neighbors_of(x)) == 3:
                a, = outside_angles(s, v, x, 1, atoms)
                assert not _inside((v[x][0] + .5 * L * cos(a), v[x][1] + .5 * L * sin(a)), poly), x


def test_artemisinin_draws_the_same_from_any_atom_order():
    mol = smiles('C[C@@H]1CC[C@H]2[C@H](C(=O)O[C@H]3[C@@]24[C@H]1CC[C@](O3)(OO4)C)C')
    assert perm(mol, controller_layout, k=10) == 0


def test_a_cis_trans_between_two_ring_systems_is_kept():
    # PubChem 5614546: an exocyclic C=C joining a pyrrolopyrazole to a cyclohexadienone
    mol = smiles('CCOC1=C(C=C(C=C1)[C@H]2C\\3=C(C(=O)N2CCO)NN/C3=C\\4/C=C(C(=CC4=O)C)Cl)OCC')
    log = []
    plane = controller_layout(mol, log)
    assert shape(mol, plane)['ct'] == 0
    assert not [r for r in log if r.rule == 'depict:open-ring']


def test_the_beam_prefers_the_drawing_that_keeps_every_cis_trans():
    # PubChem 9987107: two exocyclic C=C, one between two ring systems drawn one after the other
    mol = smiles('CN\\1C2=C(C=CC3=CC=CC=C32)S/C1=C\\4/C(=O)N(/C(=C\\C5=CC=CC=[N+]5CC=C)/S4)CC6CCCO6')
    assert shape(mol, controller_layout(mol))['ct'] == 0


# --- moves of the removed passes -------------------------------------------------------------------------

def _turned(plane, part, o, a):
    c, s = cos(a), sin(a)
    for x in part:
        dx, dy = plane[x][0] - plane[o][0], plane[x][1] - plane[o][1]
        plane[x] = plane[o][0] + c * dx - s * dy, plane[o][1] + s * dx + c * dy


def test_a_cross_targets_its_right_angles_off_every_held_bond():
    mol, s, tree = _scoped_tree('CC(C)(C)CC')
    plane = dict(tree.plane)
    slot = _dir(plane, 2, 1)
    _turned(plane, {1}, 2, radians(40))                 # the methyl off its slot
    assert any(_near(t, slot) for t in targets(s, plane, 2, 1))


def test_an_even_move_spaces_terminal_ligands_at_their_slots():
    mol = smiles('F[P-](F)(F)(F)(F)F')
    p = next(x for x in mol if mol.atom(x).atomic_symbol == 'P')
    s = Scope(mol, set(mol), mol.canonical_order(), [])
    plane = dict(controller_layout(mol))
    f = min((x for x in mol if x != p), key=s.rank.__getitem__)
    _turned(plane, {f}, p, radians(25))
    moves = [m for m in candidates(s, plane, frontier(s, plane) | {p}) if m.kind == 'even']
    assert len(moves) == 1
    new = {**plane, **apply(s, plane, moves[0])}
    angles = sorted(round(degrees(_dir(new, p, x))) % 360 for x in mol if x != p)
    gaps = [(angles[(i + 1) % 6] - angles[i]) % 360 for i in range(6)]
    assert gaps == [60] * 6
    assert add(full(s, plane), delta(s, plane, apply(s, plane, moves[0]))) == pytest.approx(full(s, new))


def test_a_pivot_move_turns_a_spiro_half_opposite_the_other():
    mol = smiles('C1CCC2(C1)CCCCC2')                    # spiro atom 4: rings 1-5 and 4, 6-10
    s = Scope(mol, set(mol), mol.canonical_order(), [])
    plane = dict(controller_layout(mol))
    _turned(plane, {6, 7, 8, 9, 10}, 4, radians(50))

    def bisector(ys):
        u, v = (_dir(plane, 4, y) for y in ys)
        return (u + ((v - u) % (2 * pi)) / 2) % (2 * pi) if (v - u) % (2 * pi) < pi else \
            (v + ((u - v) % (2 * pi)) / 2) % (2 * pi)
    moves = [m for m in candidates(s, plane, set(mol)) if m.kind == 'pivot']
    assert moves
    moved = apply(s, plane, moves[0])
    assert add(full(s, plane), delta(s, plane, moved)) == pytest.approx(full(s, {**plane, **moved}))
    plane.update(moved)
    a, b = bisector((3, 5)), bisector((6, 10))
    assert _near(a, (b + pi) % (2 * pi))


def test_a_fused_system_offers_several_polyhexes_as_variants():
    mol = smiles('C1CCC2C(C1)CCC1C2CCC2CCCC12')          # gonane: its pentagon sits on several polyhexes
    s = Scope(mol, set(mol), mol.canonical_order(), [])
    assert len(variants(s, ring_blocks(mol)[0])) == VARIANTS


def test_a_ring_seated_into_a_macrocycle_draws_the_same_from_any_atom_order():
    # PubChem 44574314: a small ring fused into a macrocycle, seated by the post-steps
    mol = smiles(
        'C[C@H]1[C@H]([C@@](C[C@@H](O1)O[C@H]2[C@@H]3[C@H]4C5=NC(=CS5)C(=O)N[C@@H](COC(=O)C6=C(CO3)'
        'C7=C(COC2=O)C=CC=C7N6O)C8=NC(=CS8)C9=N/C(=C/1\\NC(=CS1)C(=O)NCCN(C)C)/C(=O)C=C9C1=NC(=CS1)C'
        '(=O)N[C@H](C(=O)N/C(=C(\\C)/OC)/C1=NC(=CS1)C(=O)N4)[C@@H](C)O)(C)O)N(C)C')
    assert perm(mol, controller_layout, k=3) == 0


@pytest.mark.xfail(strict=False, reason='regularize, wheels and fair sum over neighbours in storage order; the noise '
                   'flips a seat placement in some atom orders')
def test_a_seat_after_noisy_ring_post_steps_draws_the_same_from_any_atom_order():
    # PubChem 25086875: a ring seated into a macrocycle after the ring post-steps
    mol = smiles(
        'CCC1=C(C2=CC3=NC(=CC4=NC(=C5C6=C(C7=C(C(=CC1=N2)N=C75)C)NC8=CC9=C(C=C8N6)NC1=C2C(=C5C=C6C('
        '=C(C(=N6)C=C6C(=C(C(=CC7=NC(=C(C2=N5)C(=O)N91)C(C7C)CCC(=O)OC)N6)C)C=C)C)C)C)C(C4C)CCC(=O)'
        'OC)C(=C3C=C)C)C')
    assert perm(mol, controller_layout, k=6) == 0


def test_a_tree_with_conjugated_cis_trans_ignores_atom_order():
    # PubChem 13907103: a cyanine whose stated double bonds are checked as the tree places them
    mol = smiles(
        'CCN\\1C2=CC(=C(C=C2N(/C1=C\\C=C\\C3=[N+](C4=C(O3)C=CC(=C4)OC)CCCCS(=O)(=O)O)CCCCS(=O)(=O)O)Cl'
        ')Cl')

    def placed(m):
        r = m.canonical_order()
        inv = {v: k for k, v in r.items()}
        plane = tree_embedding(Scope(m, set(m), r, [])).plane
        return [(round(plane[inv[k]][0], 6), round(plane[inv[k]][1], 6)) for k in range(len(m))]
    ref = placed(mol)
    assert all(placed(smiles(format(mol, 'r'))) == ref for _ in range(12))
