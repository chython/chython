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
"""The cascade: deterministic layout of the covered molecules, and its freeze."""
import pytest

from chython.core import smiles
from chython.depict.layout.cascade import GOLDEN, cascade_layout, covered, drawing_hash, fit, golden
from chython.depict.layout.primitive import Initial
from chython.depict.layout.quality import shape
from chython.depict.layout.score import terms


@pytest.mark.parametrize('text', ['CC(=O)Oc1ccccc1C(=O)O', 'CC1(C)C2CCC1(C)C(=O)C2', 'CN1C2CCC1CC(O)C2',
                                  'CC12CC3CC(C)(C1)CC(N)(C3)C2', 'CC(C)CCCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C',
                                  'OCC(O)C(O)C(O)C(O)CO'])
def test_a_drawing_ignores_atom_order(text):
    mol = smiles(text)
    ref = drawing_hash(mol, cascade_layout(mol))
    for _ in range(5):
        other = smiles(format(mol, 'r'))
        assert drawing_hash(other, cascade_layout(other)) == ref


def test_the_fit_clears_a_defect_of_the_start():
    mol = smiles('CC(C)COP(=O)(C(C)(C)N)OCC(C)C')        # pubchem 22213697
    init = Initial(mol, set(mol))
    assert init.env().defects > 0
    assert terms(init.scope, fit(init))[0] - init.inherent == 0


def test_the_fit_repairs_a_stated_cis_trans_the_tree_breaks():
    mol = smiles('CN1C(=CC=C2CCC(=CC=C3N(C)c4ccc5ccccc5c4C3(C)C)C2Cl)C(C)(C)c2c1ccc1ccccc21')
    plane = cascade_layout(mol)
    assert shape(mol, plane)['ct'] == 0


@pytest.mark.parametrize('text, ok', [('c1ccc2ccccc2c1', True), ('C1CC2CCC1C2', True), ('C1CCCCCCCCCCC1', False),
                                      ('C1CC2CC1C1CCC21', False)])
def test_covered_is_shapes_small_rings_and_edge_fused_systems(text, ok):
    assert covered(smiles(text)) is ok


def test_the_covered_corpus_drawings_are_frozen():
    want = dict(line.split('\t') for line in GOLDEN.read_text().splitlines()[1:])
    got = golden()
    assert set(got) == set(want), 'the covered set moved: bless deliberately'
    changed = sorted(i for i in want if got[i] != want[i])
    assert not changed, f'{len(changed)} frozen drawings changed, first {changed[:5]}: bless deliberately'


@pytest.mark.parametrize('text', ['C1CCC2(CC1)CCCC2', 'C1CCC2(C1)CCC(C)(C)C2', 'O=C1CCC2(CCCCC2)N1'])
def test_a_spiro_atom_keeps_both_rings_whole(text):
    from math import hypot
    from chython.depict.layout.lattice import L
    from chython.depict.layout.molecule import layout2d
    mol = smiles(text)
    plane = layout2d(mol, engine='chython', force=True)
    lengths = [hypot(plane[b.n][0] - plane[b.m][0], plane[b.n][1] - plane[b.m][1]) for b in mol.bonds()]
    assert max(lengths) < 1.1 * min(lengths)


def test_every_approved_drawing_is_kept():
    """The user's approvals are drawn as approved, whatever the algorithm: a change that alters one is wrong."""
    from chython.depict.layout.approved import read_library, same_picture
    from chython.depict.layout.molecule import layout2d
    changed = []
    for i, (text, points) in read_library().items():
        mol = smiles(text)
        plane = layout2d(mol, engine='chython', force=True, peptide=True)          # the review page's drawing
        if not same_picture(mol, dict(zip(mol, map(tuple, points))), plane):
            changed.append(i)
    assert not changed, f'{len(changed)} approved drawings changed: {changed[:8]}'


def test_a_full_analog_is_held_to_its_approved_twin(tmp_path):
    from chython.depict.layout.approved import approve, held
    path = tmp_path / 'approved.tsv'
    ref = smiles('O=C(O)c1ccccc1')
    plane = cascade_layout(ref)
    approve('x', 'O=C(O)c1ccccc1', [list(plane[x]) for x in ref], path)
    analog = smiles('S=C(N)c1ccccc1')                    # thiobenzamide: the same skeleton
    drawn = cascade_layout(analog)
    assert held(analog, drawn, path) is True
    bent = {x: (p[0], p[1] + (.5 if x == max(analog) else 0)) for x, p in drawn.items()}
    assert held(analog, bent, path) is False
    assert held(smiles('CCCCO'), cascade_layout(smiles('CCCCO')), path) is None


@pytest.mark.parametrize('text', ['CC12CCC(=O)C=C1CCC1C2C(O)CC2(C)C1CCC2(O)C(=O)CO',            # cortisol
                                  'C1CCC(CC1)(CNS(=O)(=O)C2=C(C=CC=C2F)F)[NH+]3CCOCC3'])     # pubchem 8351277
def test_the_fit_leaves_no_atom_inside_a_ring(text):
    from chython.depict.layout.relax import _inr
    from chython.depict.layout.rings import ring_blocks
    mol = smiles(text)
    plane = cascade_layout(mol)
    rings = [r for _, rs in ring_blocks(mol) for r in rs]
    assert _inr(plane, rings) == 0


@pytest.mark.parametrize('text', ['CC(C)(c1ccc(O)cc1)c1ccc(O)cc1', 'CC(C)(C1CCCCC1)C1CCCCC1'])
def test_a_cross_puts_its_two_ring_branches_on_one_line(text):
    from math import atan2, degrees, pi
    mol = smiles(text)
    plane = cascade_layout(mol)
    centre = next(x for x in mol if len(mol.neighbors_of(x)) == 4)
    ns = mol.neighbors_of(centre)
    methyls = [y for y in ns if len(mol.neighbors_of(y)) == 1]
    heavy = [y for y in ns if y not in methyls]
    def angle(a, b):
        u = atan2(plane[a][1] - plane[centre][1], plane[a][0] - plane[centre][0])
        v = atan2(plane[b][1] - plane[centre][1], plane[b][0] - plane[centre][0])
        return abs(degrees((u - v + pi) % (2 * pi) - pi))
    assert angle(*heavy) == pytest.approx(180, abs=1) and angle(*methyls) == pytest.approx(180, abs=1)


def test_a_mixture_lays_out_every_component():
    mol = smiles('CC1(C)C2CCC1(C)C(=O)C2.c1ccc2ccccc2c1.CCO')
    assert set(cascade_layout(mol)) == set(mol)


@pytest.mark.parametrize('text', ['CCOCC', 'CCCCCC', 'CCOCC.C1=CC=[C-]C=C1.N=[W](Cl)(Cl)(Cl)Cl'])   # pubchem 174819204
def test_an_unbranched_chain_is_one_zigzag(text):
    mol = smiles(text)
    plane = cascade_layout(mol)
    chain = list(mol)[:5]
    for a, b, c, d in zip(chain, chain[1:], chain[2:], chain[3:]):
        (bx, by), (cx, cy) = plane[b], plane[c]
        side = [(cx - bx) * (plane[u][1] - by) - (cy - by) * (plane[u][0] - bx) for u in (a, d)]
        assert side[0] * side[1] < 0, f'{a}-{b}-{c}-{d} folds'


@pytest.mark.parametrize('text', ['O[W](=O)(=O)O[W](=O)(=O)O[W](=O)([O-])([O-])[O-]',                # pubchem 6369825
                                  'O[Mo](=O)(O)(O)O[Mo](=O)(O)(O[Mo](=O)(O)(O)O)O[Mo](=O)(O)(O)O'])   # pubchem 6396108
def test_a_metal_oxo_chain_runs_straight_through_each_metal(text):
    from math import atan2, degrees
    mol = smiles(text)
    plane = cascade_layout(mol)

    def bridging(x, y):
        return len(mol.neighbors_of(y)) == 2 and all(mol.atom(z).is_metal for z in mol.neighbors_of(y))
    for x in mol:
        if not mol.atom(x).is_metal:
            continue
        chain = [y for y in mol.neighbors_of(x) if bridging(x, y)]
        if len(chain) != 2:
            continue
        (a, b), (px, py) = chain, plane[x]
        u = atan2(plane[a][1] - py, plane[a][0] - px)
        v = atan2(plane[b][1] - py, plane[b][0] - px)
        assert abs(degrees(u - v)) % 360 == pytest.approx(180, abs=1)
        rest = [y for y in mol.neighbors_of(x) if y not in chain]
        cx, cy = plane[b][0] - plane[a][0], plane[b][1] - plane[a][1]
        sides = [(cx * (plane[y][1] - plane[a][1]) - cy * (plane[y][0] - plane[a][0])) > 0 for y in rest]
        assert abs(sum(sides) - (len(rest) - sum(sides))) <= 1, 'terminal ligands split over both sides'


@pytest.mark.parametrize('text, x, angle', [('OP(=O)(O)OP(=O)(O)O', 5, 180),               # pyrophosphate
                                            ('O[W](=O)(=O)O[W](=O)(=O)O', 5, 180),
                                            ('O=P1(O)OP(=O)(O)OP(=O)(O)O1', 4, 120)])   # a ring oxygen keeps its ring
def test_an_acyclic_oxygen_between_two_phosphorus_or_metals_is_straight(text, x, angle):
    from math import acos, degrees, hypot
    mol = smiles(text)
    plane = cascade_layout(mol)
    a, b = mol.neighbors_of(x)
    u = (plane[a][0] - plane[x][0], plane[a][1] - plane[x][1])
    v = (plane[b][0] - plane[x][0], plane[b][1] - plane[x][1])
    assert degrees(acos(max(-1., min(1., (u[0] * v[0] + u[1] * v[1]) / hypot(*u) / hypot(*v))))) == \
        pytest.approx(angle, abs=1)


@pytest.mark.parametrize('text, x', [('CCCCCCCCCCCC(=O)C(C(=O)OCC)(F)F', 14),        # pubchem 9857711
                                     ('NC(Cl)(Cl)C(Cl)(Cl)C(Cl)(Cl)Cl', 5)])
def test_a_dihalo_carbon_of_a_chain_is_a_cross(text, x):
    from math import atan2, degrees
    mol = smiles(text)
    plane = cascade_layout(mol)
    a, b = (y for y in mol.neighbors_of(x) if len(mol.neighbors_of(y)) > 1)
    t = degrees(atan2(plane[a][1] - plane[x][1], plane[a][0] - plane[x][0]) -
                atan2(plane[b][1] - plane[x][1], plane[b][0] - plane[x][0])) % 360
    assert t == pytest.approx(180, abs=1)


def test_a_substituent_in_the_way_of_an_exit_flips_away():
    from math import hypot
    mol = smiles('CC1=CC(=O)C2=C(C=C3C(=C2O1)C(=O)C4=C(C3=O)C=C(C=C4OC(=O)C)OC(=O)C)OC(=O)C')   # pubchem 10742773
    plane = cascade_layout(mol)
    # the pyranone C=O on its bisector, the peri acetoxy turned off it
    assert hypot(plane[5][0] - plane[6][0], plane[5][1] - plane[6][1]) == \
        pytest.approx(hypot(plane[5][0] - plane[3][0], plane[5][1] - plane[3][1]), abs=1e-6)


@pytest.mark.parametrize('text, x', [('CCS(=O)(=O)Oc1ccccc1', 3),
                                     ('CCCCCC[Se](=O)(=O)CC(=O)C(F)(F)F', 7),                    # pubchem 44367920
                                     ('CCCC1=NN(C)C2=C1N=C(NC2=O)c1cc(ccc1OCC)S(=O)(=O)N1CCN(C)CC1', 24)])  # sildenafil
def test_a_sulfonyl_is_a_cross(text, x):
    from math import atan2, degrees
    mol = smiles(text)
    plane = cascade_layout(mol)
    a, b = (y for y in mol.neighbors_of(x) if mol.atom(y).atomic_symbol != 'O' or len(mol.neighbors_of(y)) > 1)
    t = degrees(atan2(plane[a][1] - plane[x][1], plane[a][0] - plane[x][0]) -
                atan2(plane[b][1] - plane[x][1], plane[b][0] - plane[x][0])) % 360
    assert t == pytest.approx(180, abs=1)


def test_peri_substituents_keep_their_bisectors():
    from math import atan2, degrees
    mol = smiles('C1=CC2=C(C3=C(C=C(C=C31)[N+](=O)[O-])[N+](=O)[O-])C(=O)N=C(N2)N')     # pubchem 23507831
    plane = cascade_layout(mol)

    def angle(a, x, b):
        return degrees(atan2(plane[a][1] - plane[x][1], plane[a][0] - plane[x][0]) -
                       atan2(plane[b][1] - plane[x][1], plane[b][0] - plane[x][0])) % 360
    for x, (a, b), y in ((6, (5, 7), 14), (17, (4, 19), 18)):       # the peri nitro and the C=O
        assert min(angle(a, x, y), 360 - angle(a, x, y)) == pytest.approx(120, abs=1)
        assert min(angle(b, x, y), 360 - angle(b, x, y)) == pytest.approx(120, abs=1)


def test_the_fit_never_pinches_two_bonds_of_an_acyclic_atom():
    from math import atan2, degrees
    mol = smiles('CC(C)(C)[Si](C)(C(C)(C)C)[Si-][Si](C)(C(C)(C)C)C(C)(C)C')            # pubchem 12173831
    plane = cascade_layout(mol)
    for x in mol:
        ns = mol.neighbors_of(x)
        if len(ns) < 3:
            continue
        t = sorted(degrees(atan2(plane[y][1] - plane[x][1], plane[y][0] - plane[x][0])) % 360 for y in ns)
        assert min([v - u for u, v in zip(t, t[1:])] + [360 - t[-1] + t[0]]) >= 60 - 1e-6


def test_a_triphosphate_runs_on_one_line():
    from math import atan2, degrees
    mol = smiles('C1=NC(=C2C(=N1)N(C=N2)[C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O)N')   # ATP
    plane = cascade_layout(mol)
    for x in mol:
        if mol.atom(x).atomic_symbol != 'P':
            continue
        chain = [y for y in mol.neighbors_of(x) if len(mol.neighbors_of(y)) > 1]
        if len(chain) == 2:
            a, b = chain
            t = degrees(atan2(plane[a][1] - plane[x][1], plane[a][0] - plane[x][0]) -
                        atan2(plane[b][1] - plane[x][1], plane[b][0] - plane[x][0])) % 360
            assert t == pytest.approx(180, abs=1)


@pytest.mark.parametrize('text', ['CC(C)/C(=C(\\C)/C(F)(F)F)/C(C(F)(F)F)(F)F',                      # pubchem 22665005
                                  'C(=C(/N(C(C(F)(F)F)(F)F)C(C(F)(F)F)(F)F)\\F)(\\C(C(C(F)(F)F)(F)F)(F)F)/C(F)(F)F'])
def test_the_fit_stretches_rather_than_bends_a_crowded_acyclic_atom(text):
    from math import atan2, degrees
    mol = smiles(text)
    plane = cascade_layout(mol)
    for x in mol:
        ns = mol.neighbors_of(x)
        if len(ns) < 3:
            continue
        t = sorted(degrees(atan2(plane[y][1] - plane[x][1], plane[y][0] - plane[x][0])) % 360 for y in ns)
        gaps = [v - u for u, v in zip(t, t[1:])] + [360 - t[-1] + t[0]]
        assert min(gaps) >= 360 / len(ns) - 5, x


def test_a_short_chain_through_a_root_is_one_zigzag():
    mol = smiles('CCCNC(=O)CNC(=O)C[C@@H]([C@H](CC(C)C)NC(=O)OC(C)(C)C)O')                 # pubchem 13992573
    plane = cascade_layout(mol)
    path = [1, 2, 3, 4, 5, 7, 8, 9, 11, 12, 13, 18, 19, 21, 22]
    for a, b, c, d in zip(path, path[1:], path[2:], path[3:]):
        (bx, by), (cx, cy) = plane[b], plane[c]
        side = [(cx - bx) * (plane[u][1] - by) - (cy - by) * (plane[u][0] - bx) for u in (a, d)]
        assert side[0] * side[1] < 0, f'{a}-{b}-{c}-{d} folds'


def test_a_long_chain_folds_to_pack():
    from chython.depict.layout.molecule import layout2d
    mol = smiles('CCCCCCCC/C=C\\CCCCCCCC(=O)OCC(OC(=O)CCCCCCC/C=C\\CCCCCCCC)COC(=O)CCCCCCC/C=C\\CCCCCCCC')  # triolein
    plane = layout2d(mol, engine='chython', force=True)
    xs, ys = [p[0] for p in plane.values()], [p[1] for p in plane.values()]
    assert max(xs) - min(xs) < 3 * (max(ys) - min(ys))      # a straight zigzag of three chains runs about 4:1


def test_a_crowded_cage_atom_is_nudged_clear():
    from chython.depict.layout.quality import defects
    mol = smiles('C1C2C3C4CC5C3C1C67C28C49C56C\\1C2C7C(=O)CC8C2C9C(=O)/C1=C/O')            # pubchem 13806692
    assert defects(mol, cascade_layout(mol))['ovl'] == 0


def test_a_cached_system_drawing_ignores_which_molecule_computed_it():
    import chython.depict.layout.cascade as cascade
    from chython.depict.layout.approved import same_picture
    indole = smiles('c1ccc2c(c1)cc(C)[nH]2')
    cascade._DRAWINGS.clear()
    alone = cascade_layout(indole)
    cascade._DRAWINGS.clear()
    cascade_layout(smiles('CCOC(=O)C(=O)NC1=CC=C(C=C1)/C=C/C2=NC3=CC=CC=C3O2'))   # pubchem 13555408, the same key
    assert same_picture(indole, alone, cascade_layout(indole))


def test_a_fused_five_ring_is_a_regular_pentagon():
    from math import hypot
    from chython.depict.layout.lattice import L
    mol = smiles('CCN1C2=CC(=C(C=C2N(/C1=C\\C=C\\C3=[N+](C4=C(O3)C=CC(=C4)OC)CCCCS(=O)(=O)O)CCCCS(=O)(=O)O)Cl)Cl')
    plane = cascade_layout(mol)                                                             # pubchem 13907103
    for ring in ([15, 19, 18, 17, 16], [9, 10, 11, 3, 4]):
        for a, b in zip(ring, ring[1:] + ring[:1]):
            assert hypot(plane[a][0] - plane[b][0], plane[a][1] - plane[b][1]) == pytest.approx(L, rel=.02)


def test_a_steroid_has_ring_d_to_the_upper_right():
    from chython.depict.layout.rings import ring_blocks
    mol = smiles('CC12CCC(=O)C=C1CCC1C2C(O)CC2(C)C1CCC2(O)C(=O)CO')                       # cortisol
    plane = cascade_layout(mol)
    rings = [r for _, rs in ring_blocks(mol) for r in rs]
    centres = [(sum(plane[x][0] for x in r) / len(r), sum(plane[x][1] for x in r) / len(r)) for r in rings]
    mx, my = sum(c[0] for c in centres) / 4, sum(c[1] for c in centres) / 4
    (dx, dy), = [c for c, r in zip(centres, rings) if len(r) == 5]
    assert dx > mx and dy > my


def test_a_polysilane_runs_straight():
    from math import atan2, degrees
    mol = smiles('CC(C)(C)[Si](C)(C(C)(C)C)[Si-][Si](C)(C(C)(C)C)C(C)(C)C')            # pubchem 12173831
    plane = cascade_layout(mol)
    a, b = mol.neighbors_of(11)
    t = degrees(atan2(plane[a][1] - plane[11][1], plane[a][0] - plane[11][0]) -
                atan2(plane[b][1] - plane[11][1], plane[b][0] - plane[11][0])) % 360
    assert t == pytest.approx(180, abs=1)


def test_crowded_atoms_are_spaced_by_a_flip_or_stretch():
    from math import hypot
    from chython.depict.layout.cascade import SPACE
    mol = smiles('C[C@@H]([C@H](CC[C@@H](CC(=O)N1[C@@H](CCC1=O)COC(C2=CC=CC=C2)(C3=CC=CC=C3)C4=CC=CC=C4)'
                 '[Si](C)(C)C5=CC=CC=C5)[Si](C)(C)C6=CC=CC=C6)C(=O)OCC7=CC=CC=C7')                 # pubchem 10653263
    plane = cascade_layout(mol)
    within2 = {x: {z for y in mol.neighbors_of(x) for z in (y, *mol.neighbors_of(y))} for x in mol}
    assert min(hypot(plane[x][0] - plane[y][0], plane[x][1] - plane[y][1]) for x in mol for y in mol
               if x < y and y not in within2[x]) >= SPACE - 1e-6


def test_a_fused_backbone_five_ring_is_a_pentagon():
    from math import hypot
    from chython.depict.layout.lattice import L
    from chython.depict.layout.molecule import layout2d
    from chython.depict.layout.scorecard import CORPUS, read_corpus
    text = next(r['smiles'] for r in read_corpus(CORPUS) if r['id'] == 'pubchem:11629453')
    mol = smiles(text)
    plane = layout2d(mol, engine='chython', force=True)
    for ring in ([39, 40, 47, 48, 38], [25, 35, 34, 27, 26]):
        for a, b in zip(ring, ring[1:] + ring[:1]):
            assert hypot(plane[a][0] - plane[b][0], plane[a][1] - plane[b][1]) == pytest.approx(L, rel=.02)


@pytest.mark.parametrize('text', ['C1=CC2=C(C3=C(C=C(C=C31)[N+](=O)[O-])[N+](=O)[O-])C(=O)N=C(N2)N',  # pubchem 23507831
                                  'CC1=CC(=C(C=C1COC2=CC=C(C=C2)OC)C3C(=C(N(C4=C3C(=O)CCC4)C5=C(C=C(C=C5)OC)'
                                  '[N+](=O)[O-])N)C#N)C'])                                        # pubchem 4104121
def test_a_nitro_stretches_its_bond_to_the_ring_not_n_o(text):
    from math import hypot
    from chython.depict.layout.lattice import L
    mol = smiles(text)
    plane = cascade_layout(mol)
    for x in mol:
        if mol.atom(x).atomic_symbol == 'N' and mol.atom(x).charge == 1:
            for y in mol.neighbors_of(x):
                if mol.atom(y).atomic_symbol == 'O':
                    assert hypot(plane[x][0] - plane[y][0], plane[x][1] - plane[y][1]) == pytest.approx(L, rel=.01)


def test_the_fit_stretches_rather_than_bends_an_sp2_chain_atom():
    from math import atan2, degrees
    mol = smiles('CCCCCCCC/C(=C/S(=O)(=O)O)/S(=O)(=O)O')                                    # pubchem 20387168
    plane = cascade_layout(mol)
    a, b = mol.neighbors_of(10)
    t = degrees(atan2(plane[a][1] - plane[10][1], plane[a][0] - plane[10][0]) -
                atan2(plane[b][1] - plane[10][1], plane[b][0] - plane[10][0])) % 360
    assert min(t, 360 - t) == pytest.approx(120, abs=5)


@pytest.mark.parametrize('text, x, y', [
    ('C[C@H]1[C@@H]([C@H]([C@H]([C@@H](O1)OCC=C)O[C@H]2[C@@H]([C@@H]([C@H]([C@@H](O2)C)OCC3=CC=CC=C3)'
     'OCC4=CC=CC=C4)OC(=O)C)OCC5=CC=CC=C5)OCC6=CC=CC=C6', 17, 19),                               # pubchem 21633268
    ('CCC1=C(C2=NC1=CC3=C(C4=C([N-]3)C(=C5[C@H]([C@@H](C(=N5)C=C6C(=C(C(=C2)[N-]6)C=C)C)C)CCC(=O)OC/C=C(\\C)/CCC'
     '[C@H](C)CCC[C@H](C)CCCC(C)C)[C@H](C4=O)C(=O)OC)C)C.[Mg+2]', 17, 30)])                     # chlorophyll a
def test_a_ring_substituent_in_the_way_settles_on_its_bisector(text, x, y):
    from math import atan2, degrees
    mol = smiles(text)
    plane = cascade_layout(mol)
    a, b = (z for z in mol.neighbors_of(x) if z != y)
    def angle(z):
        t = degrees(atan2(plane[z][1] - plane[x][1], plane[z][0] - plane[x][0]) -
                    atan2(plane[y][1] - plane[x][1], plane[y][0] - plane[x][0])) % 360
        return min(t, 360 - t)
    assert angle(a) == pytest.approx(angle(b), abs=10) and angle(a) > 90


@pytest.mark.parametrize('text', ['CCOP(=O)(COC1=C(C2=C(C=C1)C3=C(P2(=O)OC)C(=C(C=C3)OCP(=O)(OCC)OCC)C4=C('
                                  'C(=CC(=C4)C(C)(C)C)C)CP(=O)(OC)OC)C5=C(C(=CC(=C5)C(C)(C)C)C)CP(=O)(OC)OC)OCC',
                                  'C[C@@]12CC[C@@H](C([C@@H]1CC[C@@]3([C@@H]2CC=C4[C@]3(CC[C@@]5([C@H]4CC(CC5)(C)C)'
                                  'C(=O)OC6[C@@H]([C@H]([C@@H](C(O6)CO)O)O)O)C)C)(C)C)OC7[C@@H]([C@H]([C@@H]('
                                  '[C@H](O7)C(=O)O)O)O)O'], ids=['pubchem 10820014', 'pubchem 25245326'])
def test_no_substituent_lies_inside_its_atom_ring_angle(text):
    from math import atan2
    from chython.depict.layout.cascade import in_ring_angle
    from chython.depict.layout.score import Scope
    mol = smiles(text)
    plane = cascade_layout(mol)
    scope = Scope(mol, set(mol), mol.canonical_order(), [])
    for x in mol:
        for y in mol.neighbors_of(x):
            if not any(x in r and y in r for r in scope.rings):
                assert not in_ring_angle(scope, plane, x, atan2(plane[y][1] - plane[x][1],
                                                                plane[y][0] - plane[x][0])), (x, y)


def test_a_dimethylsilyl_of_a_carbon_chain_is_a_cross():
    from math import atan2, degrees
    mol = smiles('C[C@@H]([C@H](CC[C@@H](CC(=O)N1[C@@H](CCC1=O)COC(C2=CC=CC=C2)(C3=CC=CC=C3)C4=CC=CC=C4)'
                 '[Si](C)(C)C5=CC=CC=C5)[Si](C)(C)C6=CC=CC=C6)C(=O)OCC7=CC=CC=C7')                 # pubchem 10653263
    plane = cascade_layout(mol)
    for x in mol:
        if mol.atom(x).atomic_symbol == 'Si':
            a, b = (y for y in mol.neighbors_of(x) if len(mol.neighbors_of(y)) > 1)
            t = degrees(atan2(plane[a][1] - plane[x][1], plane[a][0] - plane[x][0]) -
                        atan2(plane[b][1] - plane[x][1], plane[b][0] - plane[x][0])) % 360
            assert t == pytest.approx(180, abs=1)


def test_a_calixarene_hydroxyl_takes_its_bisector():
    from math import atan2, degrees
    from chython.depict.layout.molecule import layout2d
    mol = smiles('CS(=O)(=O)C1=CC=C(C=C1)C23C[C@H]4C[C@@H](C2)CC(C4)(C3)C5=CC6=C(C(=C5)CC7=CC(=CC(=C7O)CC8=CC=CC(=C8O)'
                 'CC9=C(C(=CC=C9)C6)O)C12C[C@@H]3C[C@H](C1)CC(C3)(C2)C1=CC=C(C=C1)S(=O)(=O)C)O')   # pubchem 10865993
    plane = layout2d(mol, engine='chython', force=True)                       # the chain route
    for o in mol:
        if mol.atom(o).atomic_symbol != 'O' or mol.bond(o, next(iter(mol.neighbors_of(o)))).order != 1:
            continue
        x, = mol.neighbors_of(o)
        a, b = (z for z in mol.neighbors_of(x) if z != o)
        def angle(z):
            t = degrees(atan2(plane[z][1] - plane[x][1], plane[z][0] - plane[x][0]) -
                        atan2(plane[o][1] - plane[x][1], plane[o][0] - plane[x][0])) % 360
            return min(t, 360 - t)
        assert angle(a) == pytest.approx(angle(b), abs=10), o


def test_a_crowded_polysilane_stretches_si_si_not_si_c():
    from math import hypot
    from chython.depict.layout.lattice import L
    mol = smiles('CC(C)(C)[Si](C)(C(C)(C)C)[Si-][Si](C)(C(C)(C)C)C(C)(C)C')            # pubchem 12173831
    plane = cascade_layout(mol)
    for x in (5, 12):
        for y in mol.neighbors_of(x):
            if mol.atom(y).atomic_symbol == 'C':
                assert hypot(plane[x][0] - plane[y][0], plane[x][1] - plane[y][1]) == pytest.approx(L, rel=.01)
    a, b = mol.neighbors_of(11)
    assert hypot(plane[a][0] - plane[11][0], plane[a][1] - plane[11][1]) == \
        pytest.approx(hypot(plane[b][0] - plane[11][0], plane[b][1] - plane[11][1]), rel=1e-6)


def test_an_adamantane_top_exit_is_vertical():
    from math import atan2, degrees
    mol = smiles('CC12CC3CC(C)(C1)CC(C)(C3)C2')
    plane = cascade_layout(mol)
    top = max((x for x in mol if len(mol.neighbors_of(x)) == 4), key=lambda x: plane[x][1])
    me = next(y for y in mol.neighbors_of(top) if len(mol.neighbors_of(y)) == 1)
    assert degrees(atan2(plane[me][1] - plane[top][1], plane[me][0] - plane[top][0])) == pytest.approx(90, abs=1)


@pytest.mark.parametrize('text', ['OC12CC3CC(CC(C3)C1)C2',                                         # 1-adamantanol
                                  'CN(C)S(=O)(=O)C1=CC=C(S1)CNC(=O)CN(C)CC23CC4CC(C2)CC(C4)C3'])   # pubchem 35286896
def test_an_adamantyl_bond_bisects_its_outline_bonds(text):
    from math import atan2, degrees, hypot
    from chython.depict.layout.rings import ring_blocks
    from chython.depict.layout.shapes import _inside
    mol = smiles(text)
    plane = cascade_layout(mol)
    cage = next(b for b, rs in ring_blocks(mol) if len(b) == 10)
    pts = sorted(cage)
    hidden = {pts[i] for i in _inside([plane[x] for x in pts])}
    x = next(x for x in cage if any(y not in cage for y in mol.neighbors_of(x)))
    out, = (y for y in mol.neighbors_of(x) if y not in cage)
    ux = uy = 0.
    for y in mol.neighbors_of(x):
        if y in cage and y not in hidden:
            d = hypot(plane[x][0] - plane[y][0], plane[x][1] - plane[y][1])
            ux, uy = ux + (plane[x][0] - plane[y][0]) / d, uy + (plane[x][1] - plane[y][1]) / d
    t = degrees(atan2(plane[out][1] - plane[x][1], plane[out][0] - plane[x][0]) - atan2(uy, ux))
    assert (t + 180) % 360 - 180 == pytest.approx(0, abs=1)


@pytest.mark.parametrize('text', ['c1ccc2c(c1)cc(S(=O)(=O)C)[nH]2', 'CS(=O)(=O)OCC1CN(CC1OS(=O)(=O)C)CC2=CC=CC=C2',
                                  'CC1=C(C(=C(C=C1)S(=O)(=O)N)C)S(=O)(=O)N', 'C1=CC(=C(C=C1O)SOS(=O)(=O)[O-])O'])
def test_o_x_o_of_a_dioxo_atom_is_on_one_line(text):
    from math import atan2, degrees
    mol = smiles(text)
    plane = cascade_layout(mol)
    for x in mol:
        oxo = [y for y in mol.neighbors_of(x) if len(mol.neighbors_of(y)) == 1 and mol.bond(x, y).order == 2 and
               mol.atom(y).atomic_symbol == 'O']
        if len(oxo) == 2 and len(mol.neighbors_of(x)) == 4:
            a, b = oxo
            t = degrees(atan2(plane[a][1] - plane[x][1], plane[a][0] - plane[x][0]) -
                        atan2(plane[b][1] - plane[x][1], plane[b][0] - plane[x][0])) % 360
            assert t == pytest.approx(180, abs=1)


def test_a_stuck_fit_bends_a_bridge_and_repairs_round_it():
    mol = smiles('CC(C1CCC(C(O1)OC2C(CC(C(C2O)N(C)C(=O)CNC(=O)C(CCCCNC(=O)OCC3=CC=CC=C3)NC(=O)OCC4=CC=C(C=C4)OC)O)'
                 'NC(=O)OCC5=CC=CC=C5)NC(=O)OCC6=CC=CC=C6)NC(=O)OCC7=CC=CC=C7')                     # pubchem 13475459
    from math import hypot
    from chython.depict.layout.cascade import _fit, _off
    from chython.depict.layout.lattice import L
    initial = Initial(mol, set(mol), mol.canonical_order())
    plane, left = _fit(initial, True)
    scope = initial.scope
    assert left == (0, 0)
    assert min(hypot(plane[x][0] - plane[y][0], plane[x][1] - plane[y][1]) for x in mol for y in mol
               if x < y and y not in scope.within2[x]) > .7 * L


@pytest.mark.parametrize('text', ['CC(C)(C1C=CN(C=C1)[N+](=O)[O-])[N+](=O)[O-]',  # pubchem 12228655
                                  'CC(C)[Si](C)(C)OCC=C', 'CCC(C)(C)CC'])
def test_two_end_atoms_of_a_cross_are_on_one_line(text):
    from math import atan2, degrees
    mol = smiles(text)
    plane = cascade_layout(mol)
    for x in mol:
        ends = [y for y in mol.neighbors_of(x) if len(mol.neighbors_of(y)) == 1]
        if len(mol.neighbors_of(x)) == 4 and len(ends) == 2:
            a, b = ends
            t = degrees(atan2(plane[a][1] - plane[x][1], plane[a][0] - plane[x][0]) -
                        atan2(plane[b][1] - plane[x][1], plane[b][0] - plane[x][0])) % 360
            assert t == pytest.approx(180, abs=1)


def test_the_fit_keeps_every_chain_angle_even_by_flips_and_stretches():
    from chython.depict.layout.cascade import _off
    from chython.depict.layout.score import Scope
    mol = smiles('CCOC1=CC=C(C=C1)C(C(=O)NCC2=CC=CC=C2)N(CC3=CC=CO3)C(=O)C4=C(C(=NS4)C(=O)N)N')  # pubchem 3185560
    plane = cascade_layout(mol)
    scope = Scope(mol, set(mol), mol.canonical_order(), [])
    assert not [x for x in mol if _off(scope, plane, x)]


def test_a_chain_atom_keeps_120():
    from math import atan2, degrees
    mol = smiles('CC[C@@H](C1=NN=NN1C2CCCC2)N(CC3=CC=CC=C3)CC4=CC5=CC=CC(=C5NC4=O)C')       # pubchem 25776000
    plane = cascade_layout(mol)
    t = degrees(atan2(plane[1][1] - plane[2][1], plane[1][0] - plane[2][0]) -
                atan2(plane[3][1] - plane[2][1], plane[3][0] - plane[2][0])) % 360
    assert min(t, 360 - t) == pytest.approx(120, abs=5)


def test_a_flip_opens_an_overloaded_corner():
    from math import hypot
    from chython.depict.layout.lattice import L
    mol = smiles('CCOC(=O)C1=C(C(=C(N1)C)C(=O)OCC(=O)NC2=C(C=C(C=C2)[N+](=O)[O-])C(F)(F)F)C')  # pubchem 4285435
    plane = cascade_layout(mol)
    within2 = {x: {z for y in mol.neighbors_of(x) for z in (y, *mol.neighbors_of(y))} for x in mol}
    assert min(hypot(plane[x][0] - plane[y][0], plane[x][1] - plane[y][1]) for x in mol for y in mol
               if x < y and y not in within2[x]) > .99 * L


def test_core_macrocycle_templates_follow_the_drawn_templates():
    from chython.core.chains import MACRO, TEMPLATES
    from chython.depict.layout.templates import TEMPLATES as DRAWN
    macro = {text for _, text, *_ in DRAWN if any(len(r) >= MACRO for r in smiles(text).sssr)}
    assert macro == set(TEMPLATES)


@pytest.mark.parametrize('text, angle', [('C1COCCOCCOCCOCCOCCO1', 120), ('C1COCCOCCOCCOCCO1', 132),
                                         ('C1COCCOc2ccccc2OCCOCCOc2ccccc2O1', 120)],
                         ids=['18-crown-6', '15-crown-5', 'dibenzo-18-crown-6'])
def test_a_crown_is_a_circle_of_hexagon_arcs(text, angle):
    from math import atan2, degrees, hypot
    from chython.core.chains import cycles
    from chython.depict.layout.lattice import L
    mol = smiles(text)
    ring = cycles(mol).ring
    plane = cascade_layout(mol)
    assert covered(mol)
    for x in mol:
        for y in mol.neighbors_of(x):
            assert hypot(plane[x][0] - plane[y][0], plane[x][1] - plane[y][1]) == pytest.approx(L, rel=.02)
    n = len(ring)
    for i in range(n):
        a, x, b = ring[i - 1], ring[i], ring[(i + 1) % n]
        t = degrees(atan2(plane[a][1] - plane[x][1], plane[a][0] - plane[x][0]) -
                    atan2(plane[b][1] - plane[x][1], plane[b][0] - plane[x][0])) % 360
        want = 120 if mol.atom(x).atomic_symbol == 'C' else angle
        assert min(t, 360 - t) == pytest.approx(want, abs=1)


def test_a_crown_turns_its_oxygens_inward():
    mol = smiles('C1COCCOCCOCCOCCOCCO1')
    plane = cascade_layout(mol)
    cx = sum(p[0] for p in plane.values()) / len(plane)
    cy = sum(p[1] for p in plane.values()) / len(plane)
    r = {x: ((plane[x][0] - cx) ** 2 + (plane[x][1] - cy) ** 2) ** .5 for x in mol}
    assert max(r[x] for x in mol if mol.atom(x).atomic_symbol == 'O') < \
        min(r[x] for x in mol if mol.atom(x).atomic_symbol == 'C')
