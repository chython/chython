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
"""Layout-quality metrics: each is shown to fail on a plane built to fail it before it is trusted."""
import time
from math import cos, pi, sin

import pytest

from chython.core import ReactionContainer, smiles
from chython.depict.layout.quality import L, defects, perm, pose, rescaled, rxn, ser, shape, shift_rmsd
from chython.depict.layout.scorecard import CORPUS, layout_fn, read_corpus, score_row


def _hexagon(mol, r=L):
    return {x: (r * cos(k * pi / 3), r * sin(k * pi / 3)) for k, x in enumerate(mol)}


def test_regular_hexagon_has_no_defect():
    mol = smiles('C1CCCCC1')
    assert defects(mol, _hexagon(mol)) == {'ovl': 0, 'ovl61': 0, 'crs': 0, 'aob': 0, 'lbl': 0}


def test_collapsed_plane_fires_overlap_and_atom_on_bond():
    mol = smiles('C1CCCCC1')
    d = defects(mol, {x: (0., 0.) for x in mol})
    assert d['ovl'] == 15 and d['aob'] > 0


def test_scrambled_hexagon_crosses():
    mol = smiles('C1CCCCC1')
    xs = list(mol)
    pts = _hexagon(mol)
    plane = {xs[k]: pts[xs[j]] for k, j in enumerate((0, 2, 4, 1, 3, 5))}
    assert defects(mol, plane)['crs'] > 0


def test_collinear_disjoint_bonds_do_not_cross():
    """Two bonds on one lattice line, as a rescaled plane places them: rounding is not a crossing."""
    mol = smiles('CC.CC')
    xs = list(mol)
    pts = ((-1.5086406448240892, 1.2375), (-0.7941696867019274, 0.8250000000000001),
           (2.7781851039088816, -1.2374999999999998), (3.4926560620310436, -1.65))
    assert defects(mol, dict(zip(xs, pts)))['crs'] == 0


def test_bond_through_a_label():
    """Ethane drawn through a chlorine's box: `lbl` fires, the geometric counts stay clean."""
    mol = smiles('CC.Cl')
    a, b, c = mol
    d = defects(mol, {a: (-L / 2, 0.), b: (L / 2, 0.), c: (0., .15)})
    assert d['lbl'] == 1 and d['crs'] == 0


@pytest.mark.parametrize('smi', ['[Na+]', '[Na+].[Cl-]', 'C'])
def test_no_bond_scores_zero(smi):
    mol = smiles(smi)
    plane = {x: (2. * k, 0.) for k, x in enumerate(mol)}
    assert defects(mol, plane) == {'ovl': 0, 'ovl61': 0, 'crs': 0, 'aob': 0, 'lbl': 0}
    assert rescaled(mol, plane) == plane


def test_regular_hexagon_is_regular():
    mol = smiles('C1CCCCC1')
    s = shape(mol, _hexagon(mol))
    assert s['cv'] == pytest.approx(0, abs=1e-12) and s['rng'] == pytest.approx(0, abs=1e-12)
    assert s['sqz'] == 0 and s['ct'] == 0


def test_collapsed_ring_is_off_its_polygon():
    mol = smiles('C1CCCCC1')
    assert shape(mol, {x: (0., 0.) for x in mol})['rng'] == pytest.approx(1.)   # the circumradius of side L is L


def test_squashed_hexagon_moves_uniformity_only():
    mol = smiles('C1CCCCC1')
    plane = {x: (1.6 * p[0], .5 * p[1]) for x, p in _hexagon(mol).items()}
    s, d = shape(mol, plane), defects(mol, plane)
    assert s['cv'] > .2 and s['rng'] > .2 and not any(d.values())


def test_a_substituent_inside_its_ring_is_counted():
    mol = smiles('CC1CCCCC1')
    ring = [x for x in mol if mol.in_ring_of(x)]
    plane = _hexagon(mol.substructure(ring))
    x = next(x for x in mol if x not in plane)
    p = plane[next(iter(mol.neighbors_of(x)))]
    assert shape(mol, {**plane, x: (2 * p[0], 2 * p[1])})['inr'] == 0
    assert shape(mol, {**plane, x: (0., 0.)})['inr'] == 1


def test_zigzag_is_on_the_lattice():
    """An odd zigzag: an even one tilts its principal axis by 2.9 degrees."""
    mol = smiles('CCCCCCC')
    plane = {x: (k * L * 3 ** .5 / 2, (k % 2) * L * .5) for k, x in enumerate(mol)}
    s = shape(mol, plane)
    assert s['ang'] == pytest.approx(0, abs=1e-6) and s['axis'] == pytest.approx(0, abs=1e-6)


def test_straight_chain_is_off_the_lattice():
    mol = smiles('CCCC')
    assert shape(mol, {x: (k * L, 0.) for k, x in enumerate(mol)})['ang'] == pytest.approx(60)


def test_a_straight_azide_is_on_the_lattice():
    mol = smiles('CN=[N+]=[N-]')                            # C1 N2 N3 N4: N3 drawn straight, the rest at 120
    plane = {1: (0., 0.), 2: (L * 3 ** .5 / 2, L * .5), 3: (L * 3 ** .5, 0.), 4: (L * 3 ** .5 * 1.5, -L * .5)}
    s = shape(mol, plane)
    assert s['ang'] == pytest.approx(0, abs=1e-6) and s['lin'] == 0


def test_drawn_cis_of_a_trans_bond_is_counted():
    mol = smiles('C/C=C/C')
    a, b, c, d = mol
    trans = {a: (-1.2, .7), b: (-.4, 0.), c: (.4, 0.), d: (1.2, -.7)}
    cis = {**trans, d: (1.2, .7)}
    assert shape(mol, trans)['ct'] == 0 and shape(mol, cis)['ct'] == 1


def _canonical_lay(mol):
    """A layout depending on constitution only: atoms on a circle in canonical order."""
    order = sorted(mol, key=mol.canonical_order().__getitem__)
    n = len(order)
    return {x: (cos(2 * pi * k / n), sin(2 * pi * k / n)) for k, x in enumerate(order)}


def test_order_blind_layout_scores_zero_perm_and_arena_order_does_not():
    mol = smiles('CC(O)CNc1ccccc1')
    assert perm(mol, _canonical_lay) == 0
    assert perm(mol, lambda m: {x: (float(k), 0.) for k, x in enumerate(m)}) > 0


def test_symmetric_layout_is_perm_invariant():
    mol = smiles('c1ccccc1')
    assert perm(mol, _hexagon) == 0


def test_shift_rmsd_ignores_translation_not_rotation():
    a = [(0., 0.), (L, 0.)]
    assert shift_rmsd([(p, (p[0] + 5, p[1] - 3)) for p in a]) == pytest.approx(0, abs=1e-12)
    assert shift_rmsd([(p, (-p[0], p[1])) for p in a]) == pytest.approx(1.)


def test_pose_is_automorphism_invariant_and_orientation_aware():
    mol = smiles('C1CCCCC1')
    ref = _hexagon(mol)
    xs = list(mol)
    rolled = {xs[k]: ref[xs[(k + 1) % 6]] for k in range(6)}   # the same drawing, atoms renamed
    turned = {x: (-p[1], p[0]) for x, p in ref.items()}        # the drawing turned 90 degrees
    assert pose(mol, rolled, ref) == pytest.approx([0], abs=1e-9)
    assert pose(mol, turned, ref)[0] > .15


def test_a_mirror_misses_only_when_it_is_not_a_symmetry():
    for smi, miss in (('CCC', False), ('CC(O)CC', True)):
        mol = smiles(smi)
        ref = {x: (k * L * .866, (k % 2) * L * .5) for k, x in enumerate(mol)}
        mirrored = {x: (-p[0], p[1]) for x, p in ref.items()}
        assert (pose(mol, mirrored, ref)[0] > .15) is miss, smi


def _reaction(smi):
    r = smiles(smi)
    assert isinstance(r, ReactionContainer)
    return r


def test_rxn_counts_unchanged_atoms_only():
    r = _reaction('[CH3:1][CH2:2][OH:3]>>[CH3:1][CH2:2][O:3]C')
    (a,), (b,) = r.reactants, r.products
    pa = {x: (float(k), 0.) for k, x in enumerate(a)}
    pb = {x: (float(k) + 7, 2.) for k, x in enumerate(b)}          # a translated copy, plus the methyl
    assert rxn(r, [pa, pb]) == pytest.approx(0, abs=1e-12)
    pb[next(x for x in b if b.map_number_of(x) == 1)] = (7., 3.)
    assert rxn(r, [pa, pb]) > 0


@pytest.mark.parametrize('smi', ['CCO>>CCOC', '[CH3:1][CH3:1]>>[CH3:1]O'])
def test_rxn_without_a_usable_mapping_is_none(smi):
    r = _reaction(smi)
    assert rxn(r, [{x: (float(k), 0.) for k, x in enumerate(m)} for m in r.molecules()]) is None


def test_ser_scores_the_scaffold_against_the_first_member():
    core = smiles('c1ccccc1')
    a, b = smiles('Cc1ccccc1'), smiles('Oc1ccccc1')
    pa = {x: (float(k), 0.) for k, x in enumerate(a)}
    pb = {x: (float(k) + 3, 1.) for k, x in enumerate(b)}
    assert ser(core, [(a, pa), (b, pb)]) == pytest.approx(0, abs=1e-12)
    assert ser(core, [(a, pa)]) is None


SLICE_STRATA = ('hard', 'drug', 'metal', 'rxn')         # peptide and large exceed the budget on smilesdrawer


def _slice():
    rows = read_corpus(CORPUS)
    by = {s: [r for r in rows if r['stratum'] == s] for s in SLICE_STRATA}
    return [r for group in zip(*by.values()) for r in group][:50]


def test_corpus_slice_scores_on_smilesdrawer():
    """Fifty corpus rows scored end to end on the default engine: every metric computes, inside 10 s."""
    lay = layout_fn('smilesdrawer')
    try:
        lay(smiles('CC'))
    except ImportError:
        pytest.skip('quickjs is not installed')
    rows = _slice()
    assert len(rows) == 50
    t = time.perf_counter()
    recs = [score_row(r, lay) for r in rows]
    assert time.perf_counter() - t < 10
    ok = [r for r in recs if r['ok']]
    assert len(ok) >= 45, [(r['id'], r['error']) for r in recs if not r['ok']]
    keys = {'ovl', 'ovl61', 'crs', 'aob', 'lbl', 'cv', 'sqz', 'ang', 'rng', 'axis', 'ct', 'lin', 'inr', 'perm'}
    for r in ok:
        assert r['metrics'].keys() == ({'rxn'} if r['stratum'] == 'rxn' else keys), r['id']


def _components_in_arena_order(mol):
    """Each component order-blind, the components placed left to right in arena order."""
    plane = {}
    for k, c in enumerate(mol.connected_components):
        plane.update({x: (p[0] + 10 * k, p[1]) for x, p in _canonical_lay(mol.substructure(c)).items()})
    return plane


def test_perm_does_not_score_component_placement():
    assert perm(smiles('c1ccccc1.CCO'), _components_in_arena_order) == 0


def test_a_small_unit_plane_is_rescaled():
    mol = smiles('CCCC')
    plane = rescaled(mol, {x: (.3 * k, 0.) for k, x in enumerate(mol)})
    assert [plane[x][0] for x in mol] == pytest.approx([0., L, 2 * L, 3 * L])


def test_pose_scores_each_ring_system():
    """A ring drawn as reviewed with its side chain flipped misses nothing; an acyclic molecule is one unit."""
    mol = smiles('C1CCCCC1CCO')
    xs = list(mol)
    ref = _hexagon(xs[:6])
    ref.update({xs[6]: (2 * L, 0.), xs[7]: (3 * L, .5 * L), xs[8]: (4 * L, 0.)})
    flipped = {**ref, xs[7]: (3 * L, -.5 * L), xs[8]: (4 * L, 0.)}
    assert pose(mol, flipped, ref) == pytest.approx([0.], abs=1e-9)
    line = smiles('CCC')
    plane = {x: (k * L, 0.) for k, x in enumerate(line)}
    assert pose(line, plane, plane) == pytest.approx([0.], abs=1e-9)
