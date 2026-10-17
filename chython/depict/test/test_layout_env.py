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
"""The controller's environment: score terms, legal steps, arc moves, peptide pins, the greedy driver, features."""
from math import cos, hypot, sin

import pytest

from chython.core import smiles
from chython.core.monomers import monomers
from chython.core.test.peptides import molecule
from chython.depict.layout.edits import _arc_carry, apply, arc_moves, candidates, masked
from chython.depict.layout.env import Env, PeptideEnv
from chython.depict.layout.features import (ATOM, GLOBAL, PAIR, SPARSE, STEP, atom_features, global_features,
                                            pair_features, step_features)
from chython.depict.layout.greedy import best_step, drive, peptide_greedy
from chython.depict.layout.peptide import _Context, peptide_layout
from chython.depict.layout.rings import ring_blocks
from chython.depict.layout.score import Scope, combine, delta, delta_terms, frontier, full, terms
from chython.depict.layout.search import controller_layout
from chython.depict.layout.tree import tree_embedding
from chython.depict.peptide import routed_links
from chython.depict.test.test_peptide_layout import LAYOUTS, _Tile, _audit


def test_the_five_terms_combine_into_the_score():
    mol = smiles('CC(C)(C)CC(C)(C)OP(=O)(OC)OC(C)(C)CC(C)(C)C')
    rank = mol.canonical_order()
    tree = tree_embedding(Scope(mol, set(mol), rank, []))
    s = Scope(mol, set(mol), rank, tree.closures)
    assert combine(terms(s, tree.plane)) == full(s, tree.plane)
    for m in candidates(s, tree.plane, frontier(s, tree.plane)):
        moved = apply(s, tree.plane, m)
        d = delta_terms(s, tree.plane, moved)
        assert combine(d) == pytest.approx(delta(s, tree.plane, moved))
        assert tuple(a + b for a, b in zip(terms(s, tree.plane), d)) == pytest.approx(
            terms(s, {**tree.plane, **moved}))


def _env(text):
    mol = smiles(text)
    env = Env(mol, set(mol), mol.canonical_order(), ring_blocks(mol))
    while any(d == -1 for d in env.drawn):                  # draw every ring system first
        env.step(min((s for s in env.legal() if s.move.kind == 'embed'), key=lambda s: s.rank))
    return mol, env


def _lengths(mol, plane):
    return sorted(round(hypot(plane[b.n][0] - plane[b.m][0], plane[b.n][1] - plane[b.m][1]), 9) for b in mol.bonds())


def test_a_crankshaft_flip_keeps_every_bond_length():
    mol, env = _env('C1CCCCCCCCCCC1')
    flips = [m for m in arc_moves(env.scope, env.plane, set(mol)) if m.kind == 'arc_flip']
    assert flips
    for m in flips:
        assert _lengths(mol, {**env.plane, **apply(env.scope, env.plane, m)}) == _lengths(mol, env.plane)


def test_an_arc_whose_bridge_reaches_its_ring_cannot_move():
    mol, env = _env('C12CCC(CC1)CC2')                   # bicyclo[2.2.2]octane: bridgeheads 1 and 4
    ring = next(r for r in env.scope.rings if {1, 2, 4, 8} <= set(r))
    k = ring.index(1)
    i, j = ring[k - 1], ring[(k + 1) % len(ring)]
    assert _arc_carry(env.scope, ring, (1,), i, j) is None   # bridgehead 1 carries 6-5 to bridgehead 4


def test_an_arc_holding_one_fused_atom_carries_the_other_ring_on_its_end():
    mol, env = _env('c1ccc2ccccc2c1')                   # naphthalene: fused atoms 4 and 9
    ring = next(r for r in env.scope.rings if 4 in r and 3 in r)
    k = ring.index(4)
    assert _arc_carry(env.scope, ring, (4,), ring[k - 1], ring[(k + 1) % len(ring)]) == {4, 5, 6, 7, 8}


def test_every_arc_delta_matches_the_full_terms():
    mol, env = _env('CC1CCCCC(C)CCCCCC1')
    base = terms(env.scope, env.plane)
    moves = arc_moves(env.scope, env.plane, set(mol))
    assert {m.kind for m in moves} == {'arc_flip', 'arc_rotate', 'arc_stretch'}
    for m in moves:
        moved = apply(env.scope, env.plane, m)
        assert tuple(a + d for a, d in zip(base, delta_terms(env.scope, env.plane, moved))) == \
            pytest.approx(terms(env.scope, {**env.plane, **moved}))


def test_an_arc_turn_that_would_shorten_its_absorber_below_L_is_masked():
    mol, env = _env('C1CCCCCCCCCCC1')
    turns = [m for m in arc_moves(env.scope, env.plane, set(mol)) if m.kind == 'arc_rotate']
    bad = [m for m in turns if masked(env.scope, env.plane, apply(env.scope, env.plane, m))]
    assert bad and len(bad) < len(turns)
    for m in turns:
        new = {**env.plane, **apply(env.scope, env.plane, m)}
        short = any(hypot(new[b.n][0] - new[b.m][0], new[b.n][1] - new[b.m][1]) < .825 - 1e-9 for b in mol.bonds())
        assert masked(env.scope, env.plane, apply(env.scope, env.plane, m)) or not short


def test_the_greedy_run_stops_within_its_budget():
    mol = smiles('CC(C)CCCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C')
    env = drive(Env(mol, set(mol), mol.canonical_order(), ring_blocks(mol)), best_step)
    assert env.done and env.steps <= env.cap and env.moves[-1].kind == 'stop'


def test_stop_is_never_offered_before_every_ring_system_is_drawn():
    mol = smiles('c1ccccc1CCc1ccccc1')
    env = Env(mol, set(mol), mol.canonical_order(), ring_blocks(mol))
    assert all(s.move.kind == 'embed' for s in env.legal())


@pytest.mark.parametrize('name', LAYOUTS)
def test_a_greedy_peptide_is_clean_and_keeps_its_chords(name):
    mol = molecule(name)
    pep = monomers(mol)
    plane = peptide_layout(mol, pep, _Tile(), choose=peptide_greedy())
    beam = peptide_layout(mol, pep, _Tile())
    assert sum(v for k, v in _audit(mol, pep, plane).items() if k != 'angle') == 0
    kinds = sorted(k.kind for k in routed_links(mol, plane))
    assert kinds.count('teleport') <= sorted(k.kind for k in routed_links(mol, beam)).count('teleport')


def test_a_peptide_no_slot_takes_is_none_under_greedy():
    mol = molecule('4-aminobenzoyl')                    # its ring hangs on two non-adjacent vertices
    assert peptide_layout(mol, monomers(mol), _Tile(), choose=peptide_greedy()) is None


def test_a_peptide_environment_stops_within_its_budget():
    mol = molecule('linaclotide')
    env = drive(PeptideEnv(_Context(mol, monomers(mol), _Tile(), None)), best_step)
    assert env.done and env.steps <= env.cap and env.moves[-1].kind == 'stop'


def _stepped(mol):
    env = Env(mol, set(mol), mol.canonical_order(), ring_blocks(mol))
    for _ in range(3):                                  # a few steps in: ring systems drawn, a frontier perhaps
        if env.done:
            break
        env.step(best_step(env, env.legal()))
    return env


def _turned(env, a=.7, dx=3., dy=-2.):
    c, s = cos(a), sin(a)
    env.plane = {x: (c * p[0] - s * p[1] + dx, s * p[0] + c * p[1] + dy) for x, p in env.plane.items()}
    return env


def _rounded(rows):
    return [[round(v, 6) for v in r] for r in rows]


@pytest.mark.parametrize('text', ['CC(C)CCCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C', 'CC(=O)Oc1ccccc1C(=O)O'])
def test_features_ignore_a_turn_and_a_shift(text):
    env = _stepped(smiles(text))
    a, (pairs, p), g = atom_features(env), pair_features(env), global_features(env)
    _turned(env)
    assert _rounded(atom_features(env)) == _rounded(a)
    assert pairs == pair_features(env)[0] and _rounded(pair_features(env)[1]) == _rounded(p)
    assert _rounded([global_features(env)]) == _rounded([g])


@pytest.mark.parametrize('text', ['CC(C)CCCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C', 'OC/C=C\\CCCN'])
def test_features_ignore_atom_order(text):
    mol = smiles(text)
    a = _stepped(mol)
    b = _stepped(smiles(format(mol, 'r')))
    assert _rounded(atom_features(a)) == _rounded(atom_features(b))
    assert _rounded(pair_features(a)[1]) == _rounded(pair_features(b)[1])
    sa = sorted((k, t, tuple(round(v, 6) for v in vs)) for k, t, vs, _ in (step_features(a, s) for s in a.legal()))
    sb = sorted((k, t, tuple(round(v, 6) for v in vs)) for k, t, vs, _ in (step_features(b, s) for s in b.legal()))
    assert sa == sb


def test_every_vector_has_its_width():
    env = _stepped(smiles('CC(C)CCCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C'))
    assert {len(r) for r in atom_features(env)} == {ATOM}
    assert {len(r) for r in pair_features(env)[1]} == {PAIR}
    assert len(global_features(env)) == GLOBAL
    assert all(len(step_features(env, s)[2]) == STEP for s in env.legal())


def test_a_large_record_keeps_its_pairs_sparse():
    mol = smiles('C' * (SPARSE + 20))
    env = Env(mol, set(mol), mol.canonical_order(), [])
    pairs, _ = pair_features(env)
    assert len(pairs) < 40 * len(mol)


def test_layout2d_passes_a_chooser_to_the_controller():
    from chython.depict.layout.molecule import layout2d
    mol = smiles('CC(=O)Oc1ccccc1C(=O)O')
    a = layout2d(mol, engine='chython', force=True, choose=best_step)
    b = controller_layout(mol, choose=best_step)

    def centred(p):
        cx, cy = sum(x for x, _ in p.values()) / len(p), sum(y for _, y in p.values()) / len(p)
        return {k: (round(x - cx, 6), round(y - cy, 6)) for k, (x, y) in p.items()}
    assert centred(a) == centred(b)
    with pytest.raises(ValueError):
        layout2d(mol, engine='smilesdrawer', force=True, choose=best_step)


def test_layout2d_drives_a_peptide_backbone_with_the_chooser():
    from chython.depict.layout.molecule import layout2d
    mol = molecule('oxytocin')
    plane = layout2d(mol, engine='chython', force=True, peptide=True, choose=best_step)
    assert set(plane) == set(mol) and routed_links(mol, plane)


# --- review fixes ----------------------------------------------------------------------------------------

@pytest.mark.parametrize('text', ['C1CCCCCCCCCCC1', 'OC1C2CC3CC1CC(C2)C3'])
def test_an_arc_that_stretches_a_ring_bond_pays_for_it(text):
    from chython.depict.layout.score import STRETCHED
    mol, env = _env(text)

    def long(plane):
        return sum(hypot(plane[n][0] - plane[m][0], plane[n][1] - plane[m][1]) > STRETCHED
                   for n, m in env.scope.drawn_bonds)
    paid = 0
    for s in env._steps(arc_moves(env.scope, env.plane, set(mol))):
        gained = long({**env.plane, **s.moved}) - long(env.plane)
        assert s.terms[3] == gained                     # `stretched` counts ring bonds too
        paid += gained > 0
    assert paid


def test_a_macrocycle_is_offered_long_crankshafts_not_only_single_atoms():
    mol, env = _env('CCC1OC(=O)C(C)C(O)C(C)C(O)C(C)(O)CC(C)C(=O)C(C)C(O)C1(C)O')    # erythronolide B
    ring = max(env.scope.rings, key=len)
    flips = [m for m in arc_moves(env.scope, env.plane, set(mol)) if m.kind == 'arc_flip' and m.part & set(ring)]
    sizes = [len(m.part & set(ring)) for m in flips]
    assert len(sizes) == 8 and sum(k > 1 for k in sizes) >= 4


def _peptide_env(name):
    mol = molecule(name)
    return mol, PeptideEnv(_Context(mol, monomers(mol), _Tile(), None))


def test_a_peptide_environment_has_features_like_any_other():
    mol, env = _peptide_env('oxytocin')
    assert {len(r) for r in atom_features(env)} == {ATOM} and len(global_features(env)) == GLOBAL
    assert pair_features(env)[1]
    for s in env.legal():
        kind, target, values, atoms = step_features(env, s)
        assert len(values) == STEP and all(0 <= a < len(mol) for a in atoms if a >= 0)
        assert all(isinstance(p, tuple) for p in s.moved.values())       # atom -> point, nothing else


def test_a_pin_never_breaks_a_cis_trans_and_reports_what_is_broken():
    from chython.depict.layout.relax import _ct_broken, _frames
    mol, env = _peptide_env('alkene staple')
    frames = _frames(mol, {x: None for x in mol})
    for s in env.legal():
        if s.move.kind == 'pin':
            plane = s.state['plane']
            assert s.ct == _ct_broken(mol, plane, frames, set(mol)) and s.ct <= env.ct
    assert env.ct == _ct_broken(mol, env.plane, frames, set(mol))


def test_every_legal_step_has_its_own_token():
    for text in ('CC(C)CCCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C', 'CC(=O)Oc1ccccc1C(=O)O'):
        mol = smiles(text)
        env = Env(mol, set(mol), mol.canonical_order(), ring_blocks(mol))
        for _ in range(3):
            if env.done:
                break
            tokens = [step_features(env, s) for s in env.legal()]
            keys = [(k, t, tuple(round(v, 9) for v in vs), tuple(a)) for k, t, vs, a in tokens]
            assert len(set(keys)) == len(keys)
            env.step(best_step(env, env.legal()))


def test_clean2d_takes_the_chooser_too():
    from chython.depict.layout.molecule import clean2d
    mol = smiles('CC(=O)Oc1ccccc1C(=O)O')
    clean2d(mol, engine='chython', force=True, choose=best_step)
    assert mol.has_layout


def test_stop_orders_after_every_move():
    from chython.depict.layout.env import STOP
    mol = smiles('CCCCCC')
    env = Env(mol, set(mol), mol.canonical_order(), [])
    assert all(env.tie(s.move) < env.tie(STOP) for s in env.legal() if s.move.kind != 'stop')
    assert min(env.legal(), key=lambda s: (s.rank, env.tie(s.move))) is not None


def test_a_policy_may_edit_a_drawing_that_has_no_trouble():
    mol = smiles('CCCCCCCC')
    quiet = Env(mol, set(mol), mol.canonical_order(), [])
    assert [s.move.kind for s in quiet.legal()] == ['stop']        # the model-free chooser edits only trouble
    loud = Env(mol, set(mol), mol.canonical_order(), [], everywhere=True)
    kinds = {s.move.kind for s in loud.legal()}
    assert {'flip', 'rotate', 'stop'} <= kinds


def test_greedy_keeps_a_cis_trans_between_two_ring_systems():
    # PubChem 9987107: the unit joining two ring systems is settled only once both are drawn
    from chython.depict.layout.quality import shape
    mol = smiles('CN\\1C2=C(C=CC3=CC=CC=C32)S/C1=C\\4/C(=O)N(/C(=C\\C5=CC=CC=[N+]5CC=C)/S4)CC6CCCO6')
    assert shape(mol, controller_layout(mol, choose=best_step))['ct'] == 0


# --- feedback round 1 ------------------------------------------------------------------------------------

@pytest.mark.parametrize('text', ['CC(C)CCCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C', 'CC(=O)Oc1ccccc1C(=O)O',
                                  'CCCCCCCCCCCCCCCC(=O)OCC(COC(=O)CCCCCCC)OC(=O)CCCCCCCC'])
def test_no_legal_edit_adds_a_defect(text):
    mol = smiles(text)
    env = Env(mol, set(mol), mol.canonical_order(), ring_blocks(mol), everywhere=True)
    for _ in range(4):
        if env.done:
            break
        legal = env.legal()
        assert all(s.terms[0] <= 0 for s in legal if s.move.kind not in ('embed', 'stop'))
        env.step(best_step(env, legal))


def test_stop_waits_while_a_defect_can_be_repaired():
    mol = smiles('CCCCCCCCCC')
    env = Env(mol, set(mol), mol.canonical_order(), [], everywhere=True)
    env.plane[10] = env.plane[8]                        # fold the chain's end onto atom 8: a defect to repair
    from chython.depict.layout.score import terms
    env.terms = terms(env.scope, env.plane)
    env._legal = None
    legal = env.legal()
    assert any(s.terms[0] < 0 for s in legal)
    assert all(s.move.kind != 'stop' for s in legal)


def test_a_distorted_hexagon_costs_its_irregularity():
    from chython.depict.layout.score import TERMS
    mol, env = _env('c1ccccc1')
    ring = TERMS.index('ring')
    assert env.terms[ring] == pytest.approx(0, abs=1e-9)
    x = env.scope.atoms[0]
    moved = {x: (env.plane[x][0] + .3, env.plane[x][1])}
    assert delta_terms(env.scope, env.plane, moved)[ring] > .05


def test_three_collinear_ring_atoms_cost_a_bend():
    from chython.depict.layout.score import TERMS
    mol = smiles('C1CC2CCC1C2')                         # norbornane: bridge 7 between bridgeheads 3 and 6
    rank = mol.canonical_order()
    s = Scope(mol, set(mol), rank, [])
    plane = {x: (float(i), 0.) for i, x in enumerate(mol)}
    plane[7] = ((plane[3][0] + plane[6][0]) / 2, 0.)       # the bridge on the bridgeheads' line
    straight = TERMS.index('straight')
    assert terms(s, plane)[straight] >= 1


def test_a_macrocycle_is_offered_lattice_drawings_before_any_circle():
    from chython.depict.layout.edits import variants
    from chython.depict.layout.quality import _polygon_rmsd
    mol = smiles('CCC1OC(=O)C(C)C(O)C(C)C(O)C(C)(O)CC(C)C(=O)C(C)C(O)C1(C)O')    # erythronolide B
    s = Scope(mol, set(mol), mol.canonical_order(), [])
    block = ring_blocks(mol)[0]
    ring = max(block[1], key=len)
    for v in variants(s, block):
        assert _polygon_rmsd([v[x] for x in ring]) > .5          # a meander, far from the regular 14-gon


# --- the compiled score kernel ---------------------------------------------------------------------------

@pytest.mark.parametrize('text', ['CC(C)CCCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C', 'CC(=O)Oc1ccccc1C(=O)O',
                                  'C1CC2CCC1C2', 'CCC1OC(=O)C(C)C(O)C(C)C(O)C(C)(O)CC(C)C(=O)C(C)C(O)C1(C)O',
                                  'N#C[Fe-4](C#N)(C#N)(C#N)(C#N)C#N', 'C1CCCCCCCCCCC1'])
def test_the_kernel_scores_as_the_python_oracle(text):
    from chython.depict.layout.score import _oracle_terms, _oracle_delta, kernel_delta, kernel_terms
    mol = smiles(text)
    env = Env(mol, set(mol), mol.canonical_order(), ring_blocks(mol), everywhere=True)
    for _ in range(3):
        assert kernel_terms(env.scope, env.plane) == pytest.approx(_oracle_terms(env.scope, env.plane), abs=1e-9)
        legal = env.legal()
        for s in legal:
            if s.moved:
                assert kernel_delta(env.scope, env.plane, s.moved) == \
                    pytest.approx(_oracle_delta(env.scope, env.plane, s.moved), abs=1e-9), s.move
        if env.done:
            break
        env.step(best_step(env, legal))


def test_a_fork_is_independent_and_cheap():
    from time import perf_counter
    mol = smiles('C' * 60 + 'c1ccccc1' + 'C' * 60)
    env = Env(mol, set(mol), mol.canonical_order(), ring_blocks(mol), everywhere=True)
    env.step(next(s for s in env.legal() if s.move.kind == 'embed'))
    env.legal()
    t = perf_counter()
    copy = env.fork()
    assert perf_counter() - t < .01
    s = next(s for s in copy.legal() if s.move.kind != 'stop')
    before = dict(env.plane)
    copy.step(s)
    assert env.plane == before and env.steps == copy.steps - 1 and env.moves != copy.moves
