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
"""The primitive environment: honeycomb tree start, flips, rotations, stretches and whole-drawing turns."""
from math import cos, degrees, hypot, pi, sin

import pytest

from chython.core import smiles
from chython.depict.layout.features import ATOM, GLOBAL, STEP, atom_features, global_features, step_features
from chython.depict.layout.greedy import best_step, drive
from chython.depict.layout.lattice import L
from chython.depict.layout.primitive import Initial, KINDS
from chython.depict.layout.quality import shape
from chython.depict.layout.score import TERMS, terms


def _initial(text):
    mol = smiles(text)
    return mol, Initial(mol, set(mol))


def test_only_primitive_steps_are_offered():
    mol, init = _initial('CC(C)CCCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C')
    env = init.env()
    assert {s.move.kind for s in env.legal()} <= set(KINDS) and set(KINDS) == {'flip', 'rotate', 'stretch', 'turn',
                                                                                'stop'}


def test_the_initial_layout_is_built_once_and_every_rollout_starts_from_it():
    mol, init = _initial('c1ccccc1CCO')
    a, b = init.env(), init.env()
    s = next(s for s in a.legal() if s.move.kind == 'rotate')
    a.step(s)
    assert b.plane == init.plane and a.plane != init.plane and a.scope is b.scope


def test_a_ring_closes_by_rotations_alone():
    mol, init = _initial('C1CCCCC1')
    env = drive(init.env(), best_step)
    closure = terms(env.scope, env.plane)[TERMS.index('closure')]
    assert closure < .1 * L


def test_a_vertical_drawing_turns_horizontal():
    mol, init = _initial('CCCCCCCCCC')
    env = init.env()
    c, s = cos(pi / 2), sin(pi / 2)                     # stand the zigzag on end
    env.plane = {x: (c * p[0] - s * p[1], s * p[0] + c * p[1]) for x, p in env.plane.items()}
    env.terms = terms(env.scope, env.plane)
    env._legal = None
    tilt = TERMS.index('tilt')
    assert env.terms[tilt] > .5
    drive(env, best_step)
    assert terms(env.scope, env.plane)[tilt] < .2


@pytest.mark.parametrize('text', ['C/C=C/C=C/CCC(/C)=C/C', 'OC/C=C\\CCC/C=C/C', 'CC#CCC'])
def test_no_primitive_breaks_a_cis_trans_or_bends_a_linear_atom(text):
    mol, init = _initial(text)
    env = drive(init.env(), best_step)
    assert shape(mol, env.plane)['ct'] == 0 and shape(mol, env.plane)['lin'] == 0


def test_features_read_a_primitive_environment():
    mol, init = _initial('CC(=O)Oc1ccccc1C(=O)O')
    env = init.env()
    assert {len(r) for r in atom_features(env)} == {ATOM} and len(global_features(env)) == GLOBAL
    assert all(len(step_features(env, s)[2]) == STEP for s in env.legal())


def test_a_primitive_fork_is_independent():
    mol, init = _initial('CC(=O)Oc1ccccc1C(=O)O')
    env = init.env()
    f = env.fork()
    f.step(next(s for s in f.legal() if s.move.kind != 'stop'))
    assert env.plane == init.plane and f.steps == 1


def test_every_step_kind_has_its_own_token():
    mol, init = _initial('CCCCCC')
    kinds = {}
    for s in init.env().legal():
        kinds.setdefault(s.move.kind, set()).add(step_features(init.env(), s)[0])
    ids = [i for v in kinds.values() for i in v]
    assert len(ids) == len(set(ids)) == len(kinds)


def test_a_primitive_chooser_lays_out_from_the_tree_with_no_post_step():
    from chython.depict.layout.search import controller_layout
    seen = []

    def choose(env, legal):
        seen.append(type(env).__name__)
        return best_step(env, legal)
    choose.primitive = True
    mol = smiles('C1CCCCC1CCO')
    plane = controller_layout(mol, choose=choose)
    assert set(plane) == set(mol) and seen and set(seen) == {'PrimitiveEnv'}


def test_stop_waits_for_a_drawing_without_a_defect_until_the_budget_ends():
    mol, init = _initial('C1CCCCCCCCCCC1')              # a macrocycle: no template, the tree start overlaps
    env = init.env()
    assert env.terms[0] > 0 and all(s.move.kind != 'stop' for s in env.legal())
    env.steps = env.cap
    env._legal = None
    assert [s.move.kind for s in env.legal()] == ['stop']
    clean = drive(init.env(), best_step)
    assert clean.terms[0] == 0 and clean.moves[-1].kind == 'stop'


def test_the_step_budget_is_bounded_whatever_the_size():
    from chython.depict.layout.primitive import LIMIT
    mol, init = _initial('C' * 120)
    assert init.env().cap == LIMIT


@pytest.mark.parametrize('size', range(3, 10))
def test_a_rotation_targets_the_angle_of_every_ring_from_3_to_9_and_straight(size):
    from math import degrees
    from chython.depict.layout.primitive import ANGLES
    assert any(abs(degrees(t) - (180 - 360 / size)) < 1e-6 for t in ANGLES) and any(
        abs(degrees(t) - 180) < 1e-9 for t in ANGLES)


def _at(plane, a, p, x):
    from math import atan2, degrees
    u = atan2(plane[a][1] - plane[p][1], plane[a][0] - plane[p][0])
    v = atan2(plane[x][1] - plane[p][1], plane[x][0] - plane[p][0])
    d = abs(degrees(u - v)) % 360
    return min(d, 360 - d)


def test_a_rotation_sets_its_angle_against_the_anchor_whatever_came_before():
    mol, init = _initial('CCCCCC')
    env = init.env()
    first = next(s for s in env.legal() if s.move.kind == 'rotate' and abs(s.move.arg) < 1.6)
    env.step(first)                                     # any earlier edit
    for s in env.legal():
        if s.move.kind == 'rotate':
            m = s.move
            a = init.anchor[(m.pivot, m.other)]
            f = env.fork()
            f.step(s)
            assert _at(f.plane, a, m.pivot, m.other) == pytest.approx(abs(degrees(m.arg)), abs=1e-6)


def test_a_rotation_already_met_is_not_offered():
    mol, init = _initial('CCCCCC')
    env = init.env()
    for s in env.legal():
        if s.move.kind == 'rotate':
            m = s.move
            assert _at(env.plane, init.anchor[(m.pivot, m.other)], m.pivot, m.other) != \
                pytest.approx(abs(degrees(m.arg)), abs=1e-6) or not _same_side(env, m)


def _same_side(env, m):
    a = env.initial.anchor[(m.pivot, m.other)]
    (ax, ay), (px, py), (xx, xy) = env.plane[a], env.plane[m.pivot], env.plane[m.other]
    return ((ax - px) * (xy - py) - (ay - py) * (xx - px) > 0) == (m.arg > 0)


@pytest.mark.parametrize('text', ['C1CCCC1', 'C1CCCCCC1', 'C1CCC1'])
def test_a_ring_of_another_size_closes_by_rotations(text):
    mol, init = _initial(text)
    env = drive(init.env(), best_step)
    assert terms(env.scope, env.plane)[TERMS.index('closure')] < .1 * L


def test_the_whole_drawing_turns_both_ways_by_30_and_90():
    from math import degrees
    from chython.depict.layout.primitive import SCENE
    assert sorted(round(degrees(t)) for t in SCENE) == [-90, -30, 30, 90]


def test_a_stretch_reaches_6_L_and_no_further():
    from chython.depict.layout.primitive import LONGEST
    assert LONGEST == pytest.approx(6 * L)
    mol, init = _initial('CCCC')
    env = init.env()
    p, x, part = next(e for e in init.edges if len(e[2]) == 1)
    ux, uy = env.plane[x][0] - env.plane[p][0], env.plane[x][1] - env.plane[p][1]
    env.plane[x] = (env.plane[p][0] + ux * 5.5, env.plane[p][1] + uy * 5.5)
    env.terms = terms(env.scope, env.plane)
    env._legal = None
    out = [s for s in env.legal() if s.move.kind == 'stretch' and s.move.other == x]
    assert sorted(s.move.arg > 0 for s in out) == [False, True]
    env.plane[x] = (env.plane[p][0] + ux * 6, env.plane[p][1] + uy * 6)
    env._legal = None
    assert [s.move.arg < 0 for s in env.legal() if s.move.kind == 'stretch' and s.move.other == x] == [True]


def test_a_step_through_a_defect_is_offered():
    mol, init = _initial('CC(C)CCCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C')
    assert any(s.terms[0] > 0 for s in init.env().legal())


FUSED = ['c1ccc2[nH]ccc2c1', 'c1ccc2c(c1)Cc1ccccc12', 'CC12CCC3C(CCC4=CC(=O)CCC34C)C1CCC2O']
BRIDGED = ['C1CC2CCC1CC2', 'C1C2CC3CC1CC(C2)C3', 'CN1CCC23C4C1CC5=C2C(=C(C=C5)O)OC3C(C=C4)O']


@pytest.mark.parametrize('text', FUSED + BRIDGED)
def test_a_small_ring_system_starts_as_its_template_without_a_defect_of_its_own(text):
    mol, init = _initial(text)
    assert init.env().defects == 0


@pytest.mark.parametrize('text', FUSED)
def test_a_fused_template_stretches_no_bond(text):
    mol, init = _initial(text)
    assert init.terms[TERMS.index('stretched')] == 0


def test_no_step_reshapes_a_template():
    mol, init = _initial('CC12CCC3C(CCC4=CC(=O)CCC34C)C1CCC2O')
    inside = [b for b in init.rigid]
    assert inside
    for s in init.env().legal():
        m = s.move
        if m.kind in ('flip', 'rotate', 'stretch'):
            assert not any(m.pivot in b and m.other in b for b in init.rigid)


def test_a_macrocycle_stays_on_the_tree_with_every_edit():
    mol, init = _initial('C1CCCCCCCCCCC1')
    assert not init.rigid and {s.move.kind for s in init.env().legal()} >= {'flip', 'rotate', 'stretch'}


def _distances(plane, atoms):
    pts = [plane[x] for x in atoms]
    return sorted(round(hypot(p[0] - q[0], p[1] - q[1]) / L, 3) for i, p in enumerate(pts) for q in pts[:i])


@pytest.mark.parametrize('text, name', [('C1CC2CCC1C2', 'norbornane'), ('CC1(C)C2CCC1(C)C(=O)C2', 'norbornane'),
                                        ('C1CN2CCC1CC2', 'bicyclo[2.2.2]octane'),
                                        ('CN1C2CCC1CC(O)C2', 'bicyclo[3.2.1]octane'),
                                        ('C1C2CC3CC1CC(C2)C3', 'adamantane')])
def test_a_bridged_system_starts_as_its_predefined_shape(text, name):
    from chython.depict.layout.shapes import _shapes, skeleton
    mol, init = _initial(text)
    (block,) = init.rigid
    skel, _ = skeleton(mol, block)
    wants = [_distances({i: (x * L, y * L) for i, (x, y) in pts.items()}, list(pts))
             for _, _, pts, _, _, _, _ in _shapes()[skel.canonical_bytes]]
    assert _distances(init.plane, block) in wants


@pytest.mark.parametrize('wedged', [0, 1])
def test_a_cage_view_keeps_a_wedged_atom_off_its_bond_crossing(monkeypatch, wedged):
    """6-hydroxytropinone's ethano bridge: whichever of its two carbons carries the wedge, the tropane view chosen
    keeps that one away from where the bridges cross, at the bridge's outer end."""
    from chython.depict.layout import shapes
    mol = smiles('CN1[C@@H]2C[C@H](O)[C@H]1CC(=O)C2')
    rings = {x for x in mol if mol.atom(x).in_ring}
    heads = {x for x in rings if sum(y in rings for y in mol.neighbors_of(x)) == 3}

    def by_head(x):
        return x not in heads and any(y in heads for y in mol.neighbors_of(x))
    ethano = sorted(x for x in rings if by_head(x) and any(by_head(y) for y in mol.neighbors_of(x) if y in rings))
    x = ethano[wedged]
    monkeypatch.setattr(shapes, '_wedged', lambda m: {x})
    plane, _ = shapes.shape(mol, rings)
    skel, _ = shapes.skeleton(mol, rings)
    (_, _, pts, _, _, crossed, _), = shapes._shapes()[skel.canonical_bytes]
    spots = [(pts[i][0] * L, pts[i][1] * L) for i in crossed]
    assert min(hypot(plane[x][0] - a, plane[x][1] - b) for a, b in spots) > .5 * L


def test_porphine_starts_as_its_template_and_a_crown_ether_as_its_zigzag_circle():
    mol, init = _initial('c1cc2cc3ccc(cc4ccc(cc5ccc(cc1n2)[nH]5)n4)[nH]3')
    assert init.rigid and init.env().defects == 0
    mol, init = _initial('C1COCCOCCOCCOCCOCCO1')
    assert init.rigid and init.env().defects == 0


@pytest.mark.parametrize('text', ['CN1C2CCC1CC(O)C2', 'CC12CC3CC(C)(C1)CC(N)(C3)C2'])
def test_a_shape_s_substituents_leave_away_from_it(text):
    mol, init = _initial(text)
    assert init.env().defects == 0


@pytest.mark.parametrize('text', ['C1CC2CCC1C2', 'C1C2CC3CC1CC(C2)C3', 'CN1C2CCC1CC(O)C2',
                                  'CC12CC3CC(C)(C1)CC(N)(C3)C2', 'CC1(C)C2CCC1(C)C(=O)C2'])
def test_a_shaped_system_stands_upright_as_its_reference_drawing(text):
    from chython.depict.layout.shapes import shape
    mol, init = _initial(text)
    (block,) = init.rigid
    ref = shape(mol, block)[0]
    up = max(ref, key=lambda x: ref[x][1])                  # the reference's top atom
    assert max(block, key=lambda x: init.plane[x][1]) == up


def test_camphor_is_drawn_as_camphor_and_norbornane_as_norbornane():
    from chython.depict.layout.shapes import DRAWN
    for text, name in (('CC1(C)C2CCC1(C)C(=O)C2', 'camphor'), ('C1CC2CCC1C2', 'norbornane')):
        mol, init = _initial(text)
        (block,) = init.rigid
        want = _distances({i: (x * L, y * L) for i, (x, y) in enumerate(DRAWN[name][1])}, range(len(DRAWN[name][1])))
        assert _distances(init.plane, block) == want, name


def test_a_spiro_hinge_swaps_the_far_rings_two_bonds():
    """Spirapril: the spiro atom's hinge mirrors the dithiolane on its bisector there, the two ring neighbours trading
    places and every bond keeping its length."""
    from chython.depict.layout.primitive import PrimitiveEnv
    mol, init = _initial('CCOC(=O)[C@H](CCc1ccccc1)N[C@@H](C)C(=O)N1CC2(C[C@H]1C(O)=O)SCCS2')
    (a, x, y, part), = init.hinges
    assert {mol.atom(z).atomic_symbol for z in (x, y)} == {'S'}
    env = PrimitiveEnv(init)
    (move,) = [m for m in env._moves() if m.kind == 'flip' and m.pivot == a and m.part == part]
    moved = {**env.plane, **env._moved(move)}
    for p, q in ((x, y), (y, x)):
        assert hypot(moved[p][0] - env.plane[q][0], moved[p][1] - env.plane[q][1]) < 1e-9
    for b in mol.bonds():
        before = hypot(env.plane[b.n][0] - env.plane[b.m][0], env.plane[b.n][1] - env.plane[b.m][1])
        assert abs(hypot(moved[b.n][0] - moved[b.m][0], moved[b.n][1] - moved[b.m][1]) - before) < 1e-9


def test_every_bridgehead_of_a_drawn_cage_has_its_models_handedness():
    """A drawn shape's depth turns each bridgehead its model's way, all alike: a bridgehead drawn nearly flat, as
    pseudopelletierine's, would otherwise read against its partner."""
    from itertools import islice
    from chython.depict.layout import shapes

    def volume(p, a, skel):
        (a1, a2, a3), (b1, b2, b3), (c1, c2, c3) = [[p[n][k] - p[a][k] for k in range(3)]
                                                    for n in sorted(skel.neighbors_of(a))]
        return a1 * (b2 * c3 - b3 * c2) - a2 * (b1 * c3 - b3 * c1) + a3 * (b1 * c2 - b2 * c1)
    models = shapes._models()
    for key, found in shapes._shapes().items():
        if key not in models:
            continue
        model, xyz = models[key]
        for name, skel, pts, _, _, _, depth in found:
            if name not in shapes.DRAWN:
                continue
            heads = [x for x in pts if len(skel.neighbors_of(x)) == 3]
            drawn = {x: (*pts[x], depth[x]) for x in pts}
            signs = set()
            for mapping in islice(model.as_query().get_mapping(skel), 1):
                solid = {j: xyz[i] for i, j in mapping.items()}
                signs = {volume(drawn, x, skel) * volume(solid, x, skel) > 0 for x in heads}
            assert len(signs) == 1, name
