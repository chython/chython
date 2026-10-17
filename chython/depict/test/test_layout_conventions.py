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
"""Conventions of the drawn plane on `engine='chython'`: an xfail is one the controller has not learnt yet; the hard
constraints carry none."""
from math import atan2, degrees, hypot, isfinite, pi, radians, tan

import pytest

from chython.core import smiles
from chython.depict.layout.lattice import DIRS, L, project
from chython.depict.layout.molecule import layout2d
from chython.depict.layout.quality import _polygon_rmsd, defects, perm, rescaled, shape
from chython.depict.layout.relax import count_defects
from chython.depict.layout.rings import (cage_block, embed_block, macro_block, planar_block, ring_blocks, roundness,
                                        sphere_block)
from chython.depict.layout.cascade import cascade_layout as chython_layout
from chython.depict.layout.templates import template_block


CRAMMED = .85 * L      # two labelled atoms this near crowd each other


B12 = ('CC1=CC2=C(C=C1C)N(C=N2)[C@@H]3C([C@@H]([C@H](O3)CO)OP(=O)([O-])O[C@H](C)CNC(=O)CC[C@@]\\4([C@H]([C@@H'
       ']5[C@]6([C@@]([C@@H](C(=N6)/C(=C\\7/[C@@]([C@@H](C(=N7)/C=C\\8/C([C@@H](C(=N8)/C(=C4\\[N-]5)/C)CCC(=O)N'
       ')(C)C)CCC(=O)N)(C)CC(=O)N)/C)CCC(=O)N)(C)CC(=O)N)C)CC(=O)N)C)O.[C-]#N.[Co+3]')     # vitamin B12, cyanocobalamin


CLEAN = ['CCCCCC', 'C/C=C\\C', 'C/C=C/C', 'C/C=C/C=C/C=C\\CCO', 'c1ccccc1-c1ccccc1', 'CC(C)(C)C',
         'FS(F)(F)(F)(F)F', 'F[Pt](F)(F)(F)(F)F', 'CC(=O)Oc1ccccc1C(=O)O', 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',
         'COc1cccc(OC)c1-c1ccccc1P(C1CCCCC1)C1CCCCC1', 'C1CC2CCC1C2', 'C1CN2CCC1CC2',
         '[H]C([H])([H])C([H])([H])[H]', 'C1CC2CC3CCC4CCC(C1)C2C34', 'C1CC12CC2', 'C1CCC2(C1)CCCCC2',
         'COC1=C(C=C(C=C1)Cl)C(=O)OCC(=O)NC23CC4CC(C2)CC(C4)C3', 'CCC/C=C(\\C(=O)OCC)/F',
         'CC1(C)C2CCC1(C)C(=O)C2',
         'C(=C(/N(C(C(F)(F)F)(F)F)C(C(F)(F)F)(F)F)\\F)(\\C(C(C(F)(F)F)(F)F)(F)F)/C(F)(F)F']     # PubChem 13196961


def _ring_angles(plane, ring):
    out = []
    for i, b in enumerate(ring):
        a, c = plane[ring[i - 1]], plane[ring[(i + 1) % len(ring)]]
        t = abs(degrees(atan2(a[1] - plane[b][1], a[0] - plane[b][0]) - atan2(c[1] - plane[b][1],
                                                                               c[0] - plane[b][0]))) % 360
        out.append(min(t, 360 - t))
    return out


def _gaps(mol, plane, x):
    angs = sorted(degrees(atan2(plane[y][1] - plane[x][1], plane[y][0] - plane[x][0])) % 360
                  for y in mol.neighbors_of(x))
    return [(angs[(q + 1) % len(angs)] - angs[q]) % 360 for q in range(len(angs))]


CEPHEM_ESTER = ('CC(C)(C)OC(=O)C[N+]1=CC=C(C=C1)SCC2=C(N3[C@@H]([C@@H](C3=O)NC(=O)/C(=N\\OCC4=CC(=O)'
                'C(=CN4OC(C5=CC=CC=C5)C6=CC=CC=C6)OC(C7=CC=CC=C7)C8=CC=CC=C8)/C9=CSC(=N9)NC(C1=CC=CC=C1)'
                '(C1=CC=CC=C1)C1=CC=CC=C1)SC2)C(=O)OCC1=CC=C(C=C1)OC')     # PubChem 16172150


def _convex(plane, ring):
    n = len(ring)
    turns = []
    for i in range(n):
        a, b, c = plane[ring[i - 1]], plane[ring[i]], plane[ring[(i + 1) % n]]
        turns.append((b[0] - a[0]) * (c[1] - b[1]) - (b[1] - a[1]) * (c[0] - b[0]))
    return all(t > 1e-6 for t in turns) or all(t < -1e-6 for t in turns)


PORPHINE = 'C1=CC2=CC3=CC=C(N3)C=C4C=CC(=N4)C=C5C=CC(=N5)C=C1N2'                  # PubChem 66868


CALIXARENE = 'CC(C)(C)c1cc2Cc3cc(C(C)(C)C)cc(c3O)Cc3cc(C(C)(C)C)cc(c3O)Cc3cc(C(C)(C)C)cc(Cc(c1)c2O)c3O'  # tBu-calix[4]


def _inside(p, poly):
    c = False
    for (x1, y1), (x2, y2) in zip(poly[-1:] + poly[:-1], poly):
        if (y1 > p[1]) != (y2 > p[1]) and p[0] < (x2 - x1) * (p[1] - y1) / (y2 - y1) + x1:
            c = not c
    return c


RIFAMYCIN = ('CCOC1=CC=C(C=C1)NC=C2C(=C3C4=C(C2=NC(=O)C(=CC=C[C@@H]([C@@H]([C@H]([C@H]([C@H]([C@@H]([C@@H]([C@H](C=CO'
             '[C@@]5(C(=O)C3=C(O5)C(=C4O)C)C)OC)C)OC(=O)C)C)O)C)O)C)C)O)O')


RAPAMYCIN = ('C[C@@H]1CC[C@H]2C[C@@H](/C(=C/C=C/C=C/[C@H](C[C@H](C(=O)[C@@H]([C@@H](/C(=C/[C@H](C(=O)C[C@H](OC(=O)'
             '[C@@H]3CCCCN3C(=O)C(=O)[C@@]1(O2)O)[C@H](C)C[C@@H]4CC[C@H]([C@@H](C4)OC)O)C)/C)O)OC)C)C)/C)OC')


C60 = ('C12=C3C4=C5C6=C1C7=C8C9=C1C%10=C%11C(=C29)C3=C2C3=C4C4=C5C5=C9C6=C7C6=C7C8=C1C1=C8C%10=C%10C%11=C2C2=C3C3=C4'
       'C4=C5C5=C%11C%12=C(C6=C95)C7=C1C1=C%12C5=C%11C4=C3C3=C5C(=C81)C%10=C23')        # PubChem 123591


_TAXOL = ('CC1=C2[C@H](C(=O)[C@@]3([C@H](C[C@@H]4[C@]([C@H]3[C@@H]([C@@](C2(C)C)(C[C@@H]1OC(=O)[C@@H]([C@H]'
          '(C5=CC=CC=C5)NC(=O)C6=CC=CC=C6)O)O)OC(=O)C7=CC=CC=C7)(CO4)OC(=O)C)O)C)OC(=O)C')


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
def test_a_macrocycle_turns_its_crowded_atoms_to_outward_corners():
    # PubChem 23573650: a spiro piperidine on a 23-membered lactam ring
    mol = smiles('CCCCC(=O)NC1CC(=O)NCCCCC(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C2(CCN(CC2)C3=CC=C(C=C3)C)NC1=O)CC4=CC=CC=C4)'
                 'CCCN=C(N)N)CC5=CNC6=CC=CC=C65)C(=O)N')
    rank = mol.canonical_order()
    atoms, rings = max(ring_blocks(mol), key=lambda b: len(b[0]))
    pts, conflicts = embed_block(mol, atoms, rings, rank)
    ring, = rings
    poly = [project(pts[x]) for x in ring]

    def outside(x, d):
        px, py = project((pts[x][0] + d[0] / 2, pts[x][1] + d[1] / 2))
        inside = False
        for (x1, y1), (x2, y2) in zip(poly, poly[1:] + poly[:1]):
            if (y1 > py) != (y2 > py) and px < (x2 - x1) * (py - y1) / (y2 - y1) + x1:
                inside = not inside
        return not inside
    spiro, = (x for x in ring if len(mol.neighbors_of(x)) == 4)
    assert sum(outside(spiro, d) for d in DIRS) >= 2
    d = defects(mol, chython_layout(mol))
    assert d['ovl'] == 0


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
def test_an_odd_macrocycle_keeps_its_exits_off_the_bulge():
    # PubChem 12715322: a 37-membered ring bulged at one edge; the tosyl N sits on a 120-degree corner
    mol = smiles('CC1=CC=C(C=C1)S(=O)(=O)N2CCCCCCCCCCC(=O)N(CCCCCCCCCCCCN(C(=O)CCCCCCCCCC2)C)C')
    pts = chython_layout(mol)
    for x in (11, 12, 49):
        a = sorted(atan2(pts[y][1] - pts[x][1], pts[y][0] - pts[x][0]) % (2 * pi) for y in mol.neighbors_of(x))
        gaps = [(a[(i + 1) % len(a)] - a[i]) % (2 * pi) for i in range(len(a))]
        assert all(abs(g - 2 * pi / 3 * round(g / (2 * pi / 3))) < .03 for g in gaps)     # 120 or 240 degrees


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
@pytest.mark.parametrize('s', CLEAN)
def test_a_plane_with_no_defect_and_honoured_cis_trans(s):
    mol = smiles(s)
    plane = chython_layout(mol)
    assert set(plane) == set(mol)
    d = defects(mol, plane)
    assert (d['ovl'], d['aob']) == (0, 0) and _own_crossings(mol, plane) == 0
    assert shape(mol, plane)['ct'] == 0


def _own_crossings(mol, plane):
    """Crossing bond pairs, those inside one pseudo-3D shape left out: such a drawing crosses by design."""
    from chython.depict.layout.relax import _cross
    from chython.depict.layout.shapes import shape as predefined
    shaped = [set(a) for a, rs in ring_blocks(mol) if len(rs) > 1 and predefined(mol, set(a)) is not None]
    bonds = [(b.n, b.m) for b in mol.bonds()]
    return sum(1 for k, (a, b) in enumerate(bonds) for c, e in bonds[:k] if len({a, b, c, e}) == 4 and
               _cross(plane[a], plane[b], plane[c], plane[e]) and not any({a, b, c, e} <= x for x in shaped))


@pytest.mark.parametrize('s', ['C1CCC/C=C\\CCCC1', 'C1CCCCC/C=C/CCCCC1', 'C1CC/C=C\\CCC/C=C/CCC1',
                               'COC(=O)C12C/C=C\\CC3=CC=CC=C3C1CNC2=O'])
def test_a_ring_double_bond_keeps_its_cis_trans(s):
    mol = smiles(s)
    plane = chython_layout(mol)
    assert shape(mol, plane)['ct'] == 0
    d = defects(mol, plane)
    assert (d['ovl'], d['crs'], d['aob']) == (0, 0, 0)


@pytest.mark.parametrize('s', ['CC(=O)Oc1ccccc1C(=O)O', 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C', 'C12C3C4C1C5C2C3C45',
                               'COc1cccc(OC)c1-c1ccccc1P(C1CCCCC1)C1CCCCC1', 'O.O.[Na+].[Cl-]',
                               'CCC/C=C(\\C(=O)OCC)/F', 'C1CCC2(C1)CCCCC2',
                               'CC(=O)N1C2(N(C(=N1)C3=CC=C(C=C3)OC)C)N(N=C(S2)C4=CC=CC=C4)C5=CC=CC=C5', B12])
def test_the_plane_ignores_atom_order(s):
    assert perm(smiles(s), chython_layout, k=10) == 0


def test_a_symmetric_cage_draws_the_same_from_either_atom_order():
    # twistane: its candidate outer cycles tie but for float noise, which must not decide between them
    one, two = smiles('C1CC2CC3C1CC2CC3'), smiles('C1C2C3CC(C(C2)CC3)C1')

    def by_rank(mol):
        rank = mol.canonical_order()
        plane = chython_layout(mol)
        return [c for _, x, y in sorted((rank[x], p[0], p[1]) for x, p in plane.items()) for c in (x, y)]
    assert by_rank(one) == pytest.approx(by_rank(two), abs=1e-4)


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
def test_hexane_is_a_horizontal_zigzag():
    mol = smiles('CCCCCC')
    plane = chython_layout(mol)
    ys = sorted({round(y, 6) for _, y in plane.values()})
    assert len(ys) == 2
    assert shape(mol, plane)['ang'] == pytest.approx(0, abs=1e-6)


@pytest.mark.parametrize('s', ['', 'C', '[Na+]', '[H][H]', 'O.O.O', 'C.CC.CCC', '[H]', 'C1CCCCC1.[Cl-]'])
def test_every_atom_gets_a_finite_point(s):
    mol = smiles(s) if s else smiles('C')
    if not s:
        with mol.edit() as e:
            e.delete_atom(next(iter(mol)))
    plane = chython_layout(mol)
    assert set(plane) == set(mol)
    assert all(isfinite(v) for p in plane.values() for v in p)
    assert defects(mol, plane)['ovl'] == 0


def test_components_do_not_overlap():
    mol = smiles('c1ccccc1.c1ccccc1.CCCC')
    xs = []
    for c in mol.connected_components:
        plane = chython_layout(mol)
        xs.append((min(plane[x][0] for x in c), max(plane[x][0] for x in c)))
    xs.sort()
    assert all(a[1] < b[0] for a, b in zip(xs, xs[1:]))


def test_layout2d_dispatches_to_the_engine():
    mol = smiles('c1ccccc1O')
    plane = layout2d(mol, engine='chython', force=True, peptide=False)
    assert set(plane) == set(mol)


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
@pytest.mark.parametrize('s', [
    'CC#CC', 'N#Cc1ccccc1', 'CC(C)(C#N)C#N', 'CC#CC#CC#CC', 'C=C=C', 'CN=[N+]=[N-]', 'O=C=NC1CCCCC1',
    'N#CC(C#N)(C#N)C#N',                                    # PubChem 12716: tetracyanomethane
    'C#CC1(O)CCC2C3CCc4cc(O)ccc4C3CCC21C',                  # PubChem 5991: ethinylestradiol
    'N#CCC(c1ccccc1)(c1ccccc1)C(=O)N1CCCC1C#N',
])
def test_a_triple_bond_or_a_cumulated_carbon_is_straight(s):
    mol = smiles(s)
    assert shape(mol, chython_layout(mol))['lin'] == 0


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
@pytest.mark.parametrize('s', [
    'CC(C)(C)C',                                            # PubChem 10041: neopentane
    'CS(N)(=O)=O',                                          # PubChem 72879: methanesulfonamide
    'CCOP(=O)(OCC)OCC',                                     # PubChem 6535: triethyl phosphate
    'CCCC(C)(C)CCC',
    'CC(C)(C)c1ccccc1',                                     # PubChem 7366: tert-butylbenzene
    'OC(=O)CS(=O)(=O)c1ccccc1',
])
def test_an_acyclic_atom_of_four_neighbours_is_a_right_angled_cross(s):
    mol = smiles(s)
    plane = chython_layout(mol)
    four = [x for x in mol if len(mol.neighbors_of(x)) == 4 and not mol.in_ring_of(x)]
    assert four
    for x in four:
        assert all(abs(g - 90) < 5 for g in _gaps(mol, plane, x)), (x, _gaps(mol, plane, x))
    d = defects(mol, plane)
    assert (d['ovl61'], d['crs'], d['aob']) == (0, 0, 0)


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
@pytest.mark.parametrize('s, centres', [
    ('CC(C)(C)[Si](C)(C)OCC(=C)O[Si](Cl)(Cl)Cl', (2, 5, 13)),                      # PubChem 10496032
    ('CC(C)(C)CC(C)(C)OP(=O)(OC)OC(C)(C)CC(C)(C)C', (2, 6, 10, 15, 19)),            # PubChem 14324920
])
def test_neighbouring_crosses_are_squared_together(s, centres):
    mol = smiles(s)
    plane = chython_layout(mol)
    for x in centres:
        assert all(abs(g - 90) < 1e-6 for g in _gaps(mol, plane, x)), (x, _gaps(mol, plane, x))


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
@pytest.mark.parametrize('s, centres', [
    ('CCCCCCCCCCCCS(=O)(=O)N[C@@H](CC(=O)[O-])C[N+](C)(C)C', (13,)),                  # PubChem 11731938
    ('CCCC[Sn](CCCC)(CCCC)C(C)N(C(=O)OC)C(=O)OC(C)(C)C', (5,)),                      # PubChem 11123834
])
def test_a_cross_turns_a_crowding_branch_aside(s, centres):
    mol = smiles(s)
    plane = chython_layout(mol)
    for x in centres:
        assert all(abs(g - 90) < 1e-6 for g in _gaps(mol, plane, x)), (x, _gaps(mol, plane, x))
    d = defects(mol, plane)
    assert (d['ovl61'], d['crs'], d['aob']) == (0, 0, 0)


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
@pytest.mark.parametrize('s, pairs', [
    ('CC(C)O[Al](OC(=O)CS)OC(=O)CS', ((4, 5), (5, 6), (5, 11))),                  # PubChem 16697539
    ('CCC[C@H](C=C)O[Si](C)(C)C(C)(C)C', ((7, 8),)),                              # PubChem 10932826
])
def test_a_bond_between_two_labels_one_of_two_letters_is_stretched(s, pairs):
    mol = smiles(s)
    plane = chython_layout(mol)

    def length(n, m):
        return hypot(plane[n][0] - plane[m][0], plane[n][1] - plane[m][1])
    for n, m in pairs:
        assert length(n, m) > 1.2 * L, (n, m)
    stretched = {frozenset(p) for p in pairs}
    for b in mol.bonds():
        if frozenset((b.n, b.m)) not in stretched:
            assert length(b.n, b.m) == pytest.approx(L, abs=1e-6), (b.n, b.m)
    if (4, 5) in pairs:                                     # the isopropoxide keeps its 120 degrees at oxygen
        assert sorted(_gaps(mol, plane, 4))[0] == pytest.approx(120, abs=1e-6)
    d = defects(mol, plane)
    assert (d['ovl61'], d['crs'], d['aob']) == (0, 0, 0)


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
@pytest.mark.parametrize('s, hubs', [
    ('CCCCCCCC/C(=C/S(=O)(=O)O)/S(=O)(=O)O', ((9, 15),)),                          # PubChem 20387168
    ('CCCCCCCCCCCC(=O)OCCC(C)(O[Si](C)(C)C)P(=O)(OC)OC', ((17, 24),)),            # PubChem 15666550
])
def test_a_crowded_labelled_hub_is_pushed_out(s, hubs):
    mol = smiles(s)
    plane = chython_layout(mol)
    for n, h in hubs:
        assert hypot(plane[n][0] - plane[h][0], plane[n][1] - plane[h][1]) > 1.2 * L, (n, h)
        side, todo = {h}, [h]
        while todo:
            for y in mol.neighbors_of(todo.pop()):
                if y != n and y not in side:
                    side.add(y)
                    todo.append(y)
        rest = [y for y in mol if y not in side]
        room = min(hypot(plane[x][0] - plane[y][0], plane[x][1] - plane[y][1])
                   for x in side for y in rest if y not in mol.neighbors_of(x))
        assert room > 1.1 * L, (h, room / L)
    d = defects(mol, plane)
    assert (d['ovl61'], d['crs'], d['aob']) == (0, 0, 0)


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
@pytest.mark.parametrize('s', ['[N-]=[N+]=[N-]', '[N-]=O', 'CC(C)(C)C#N'])     # azide, nitric oxide, pivalonitrile
def test_a_short_linear_molecule_lies_horizontal(s):
    mol = smiles(s)
    plane = chython_layout(mol)
    a, b = next((b.n, b.m) for b in mol.bonds() if mol.hybridization_of(b.m) in (3, 5) or len(mol) == 2)
    assert abs(plane[a][1] - plane[b][1]) < 1e-6


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
@pytest.mark.parametrize('s, a, b', [('C1CCC2(CC1)CCCC2', 1, 4),        # spiro[4.5]decane, across its cyclohexane
                                      ('CCN(CC)P(CC)CC', 3, 6)])            # PubChem 2785275, its N-P bond
def test_a_mirror_through_a_spiro_atom_or_an_acyclic_bond_lies_horizontal(s, a, b):
    plane = chython_layout(smiles(s))
    assert abs(plane[a][1] - plane[b][1]) < 1e-6


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
def test_a_benzene_keeps_its_vertical_bonds():
    mol = smiles('Cc1ccccc1')
    plane = chython_layout(mol)
    assert any(abs(plane[b.n][0] - plane[b.m][0]) < 1e-6 for b in mol.bonds())


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
def test_a_lone_four_arm_centre_is_an_upright_cross():
    mol = smiles('Br[Mo](Br)(Br)Br')                                     # PubChem 160388133, its MoBr4
    plane = chython_layout(mol)
    x = next(x for x in mol if len(mol.neighbors_of(x)) == 4)
    assert sorted(round(degrees(atan2(plane[y][1] - plane[x][1], plane[y][0] - plane[x][0]))) % 360
                  for y in mol.neighbors_of(x)) == [0, 90, 180, 270]


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
@pytest.mark.parametrize('s', ['C1=CC=C(C=C1)N=[Mo](Cl)(Cl)(Cl)Cl',             # PubChem 141304503
                               'O[W](=O)(=O)O[W](=O)(=O)O[W](=O)([O-])([O-])[O-]'])   # PubChem 6369825
def test_a_five_arm_centre_spaces_its_arms_equally(s):
    mol = smiles(s)
    plane = chython_layout(mol)
    for x in mol:
        if len(mol.neighbors_of(x)) == 5:
            assert all(abs(g - 72) < 1 for g in _gaps(mol, plane, x)), x


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
@pytest.mark.parametrize('s', ['C[C@](F)(Cl)CC',
                               'C[C@](CC(=C)F)(C(=O)NC)N'])        # PubChem 11182770
def test_a_stereocentre_is_a_cross(s):
    mol = smiles(s)
    plane = chython_layout(mol)
    x = next(x for x in mol if len(mol.neighbors_of(x)) == 4)
    assert all(abs(g - 90) < 5 for g in _gaps(mol, plane, x))


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
def test_two_bonded_crosses_both_square_on_a_stretched_bond():
    # PubChem 19087371: C7 and C11 square only with their bond stretched and a crowded arm mirrored or pushed alone
    mol = smiles('C1CCC(CC1)C(C(=O)O)(C(C2CCCCC2)(C(=O)O)OC3CCCCC3)OC4CCCCC4')
    plane = chython_layout(mol)
    assert all(abs(g - 90) < 5 for x in (7, 11) for g in _gaps(mol, plane, x))


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
def test_a_splayed_atom_spreads_on_a_stretched_bond_behind_it():
    # PubChem 16172150: the benzhydryl C39 of the N-O ether spreads onto thirds, the N-O bond stretched, O38 at 120
    mol = smiles(CEPHEM_ESTER)
    plane = chython_layout(mol)
    assert all(abs(g - 120) < 5 for g in _gaps(mol, plane, 39))
    assert min(_gaps(mol, plane, 38)) == pytest.approx(120, abs=5)


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
def test_a_splayed_atom_on_a_ring_stem_spreads_on_the_stretched_bond_into_the_ring():
    # PubChem 13038664: the ylide C15 keeps 120 with its pyridine bond stretched; the ring side, the smaller end, moves
    mol = smiles('COC(=O)C1=C(C(=NC(=C1C(=O)OC)C(=P(C2=CC=CC=C2)(C3=CC=CC=C3)C4=CC=CC=C4)C5=CC=CC=C5)C(=C(Cl)Cl)Cl)'
                 'C(=C(Cl)Cl)Cl')
    plane = chython_layout(mol)
    assert all(abs(g - 120) < 5 for g in _gaps(mol, plane, 15))


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
def test_a_crowded_ring_branch_turns_outward_on_a_stretched_bridge():
    # PubChem 25242259: the OH of C21 sits on its ring's bisector; the next sugar, crowding it, moves out on O25-C26
    mol = smiles('C[C@@H]1[C@H]([C@@H]([C@H](O[C@H]1O[C@H]2[C@@H]([C@H](O[C@H]([C@@H]2O)O[C@@H]3[C@H](O[C@@H]'
                 '([C@@H]([C@H]3O)O)O[C@@H]4[C@H]([C@@H](CO[C@H]4OC5CC[C@]6(C(C5(C)C)CCC7C6=C[C@@H](C89[C@]7(CC'
                 '[C@@]8([C@](OC9=O)(C)CCCC(C)(C)O)O)C)O)C)OS(=O)(=O)O)O)C)C)O)CO)O)OC')
    plane = chython_layout(mol)
    assert all(abs(g - 120) < 5 for g in _gaps(mol, plane, 21))


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
def test_a_crowded_bulky_ring_branch_keeps_the_bisector_on_a_stretched_bond():
    # PubChem 13038664: the trichlorovinyl on C6, crowded by the one on C7, sits at 120 degrees one longer bond out
    mol = smiles('COC(=O)C1=C(C(=NC(=C1C(=O)OC)C(=P(C2=CC=CC=C2)(C3=CC=CC=C3)C4=CC=CC=C4)C5=CC=CC=C5)'
                 'C(=C(Cl)Cl)Cl)C(=C(Cl)Cl)Cl')
    plane = chython_layout(mol)
    assert all(abs(g - 120) < 5 for g in _gaps(mol, plane, 6))
    assert hypot(plane[6][0] - plane[46][0], plane[6][1] - plane[46][1]) > 1.2 * L


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
def test_a_cross_with_a_long_arm_moves_out_along_its_held_bond():
    # PubChem 13038664: the ylide P16 leaves C15 on a long P=C so that its three phenyls sit at L, square
    mol = smiles('COC(=O)C1=C(C(=NC(=C1C(=O)OC)C(=P(C2=CC=CC=C2)(C3=CC=CC=C3)C4=CC=CC=C4)C5=CC=CC=C5)'
                 'C(=C(Cl)Cl)Cl)C(=C(Cl)Cl)Cl')
    plane = chython_layout(mol)
    assert all(abs(hypot(plane[16][0] - plane[y][0], plane[16][1] - plane[y][1]) - L) < .05 * L for y in (17, 23, 29))
    assert all(abs(g - 90) < 5 for g in _gaps(mol, plane, 16))


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
def test_a_bend_stretches_its_stem_bond_to_keep_120():
    # PubChem 11731938: C18 bends 120 between C17 and the carboxylate, the C17-C18 bond pushed out past the NMe3
    mol = smiles('CCCCCCCCCCCCS(=O)(=O)N[C@@H](CC(=O)[O-])C[N+](C)(C)C')
    plane = chython_layout(mol)
    assert hypot(plane[17][0] - plane[18][0], plane[17][1] - plane[18][1]) > 1.2 * L
    assert min(_gaps(mol, plane, 18)) > 115


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
def test_a_cis_sulfonyl_keeps_120_on_a_stretched_double_bond():
    # PubChem 20387168: C10 stays at 120 degrees, the C9=C10 bond pushed out instead
    mol = smiles('CCCCCCCC/C(=C/S(=O)(=O)O)/S(=O)(=O)O')
    plane = chython_layout(mol)
    assert hypot(plane[9][0] - plane[10][0], plane[9][1] - plane[10][1]) > 1.2 * L
    assert min(_gaps(mol, plane, 10)) > 115


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
def test_a_repair_flip_leaves_an_alkene_carbon_at_120():
    # PubChem 13196961: pulling back a stretched C-F arm flips the chain, C2 of the C=C staying at 120 degrees
    mol = smiles('C(=C(/N(C(C(F)(F)F)(F)F)C(C(F)(F)F)(F)F)\\F)(\\C(C(C(F)(F)F)(F)F)(F)F)/C(F)(F)F')
    plane = chython_layout(mol)
    assert min(_gaps(mol, plane, 2)) > 115


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
@pytest.mark.parametrize('text, x', [('CC(=O)NC1CC1', 5), ('c1ccccc1C1CC1', 7), ('CC(=O)NC1CCC1', 5)])
def test_a_small_ring_exit_bisects_the_angle_outside_the_ring(text, x):
    # PubChem 47220652: the bond off a cyclopropyl or cyclobutyl atom leaves at equal angles to both ring bonds
    mol = smiles(text)
    plane = chython_layout(mol)
    gaps = sorted(_gaps(mol, plane, x))
    assert abs(gaps[1] - gaps[2]) < 2


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
def test_crowded_heteroatoms_spread_by_a_mirror():
    # PubChem 42560872: the amide and the ester mirror across their bonds, so no two heteroatoms sit a bond apart;
    # the anilide N10 keeps 120, the ester stem C18-C19 stretching instead
    mol = smiles('CC(=O)C1=CC=C(C=C1)NC(=O)[C@H]2CC=CC[C@@H]2C(=O)OC[C@H]3CCCO3')
    plane = chython_layout(mol)
    het = [n for n in mol if mol.atom(n).atomic_symbol != 'C']
    near = min(hypot(plane[a][0] - plane[b][0], plane[a][1] - plane[b][1])
               for i, a in enumerate(het) for b in het[i + 1:] if b not in mol.neighbors_of(a))
    assert near > 1.3 * L, near / L
    assert abs(min(_gaps(mol, plane, 10)) - 120) < 5
    assert all(hypot(plane[a][0] - plane[b][0], plane[a][1] - plane[b][1]) < 1.55 * L
               for a in mol for b in mol.neighbors_of(a))


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
def test_a_ring_metal_fans_its_exits_evenly_outside_the_ring():
    # PubChem 6368903: Mo14 and Mo17 of the middle Mo2O2 square space their terminal oxygens evenly about the bridge
    mol = smiles('O[Mo](=O)(=O)O[Mo]1(=O)(O[Mo](=O)(O1)(O)O[Mo]2(=O)(O[Mo](=O)(O2)(O)O[Mo]3(=O)(O[Mo](=O)(=O)O3)O)O)O')
    plane = chython_layout(mol)
    for x in (14, 17):
        gaps = sorted(_gaps(mol, plane, x))
        assert gaps[-2] - gaps[0] < 5, (x, gaps)


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
def test_a_pentagon_with_two_heavy_atoms_turns_regular_on_their_bond():
    # PubChem 27591152: the pyrrolidine and the oxazolidine are regular, each held on another atom or bond
    mol = smiles('CC1(COCN1C(=O)C2=CC3=CC=CC=C3N2C[C@@H]4CCCN4CC5=CC(=CC=C5)C(F)(F)F)C')
    plane = chython_layout(mol)
    for ring in ((19, 20, 21, 22, 23), (2, 3, 4, 5, 6)):
        assert all(abs(hypot(plane[a][0] - plane[b][0], plane[a][1] - plane[b][1]) - L) < .05 * L
                   for a, b in zip(ring, ring[1:] + ring[:1]))


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
def test_a_fused_pentagon_seeds_its_long_lattice_edge_off_the_fused_bond():
    # PubChem 11179070: the furanose and the dioxolane fused on C3-C4 both draw regular
    mol = smiles('C[C@H]1[C@@H]2[C@H](O[C@H]1[C@@H](CCC(C)OCC3=CC=CC=C3)CO[Si](C)(C)C(C)(C)C)OC(O2)(C)C')
    plane = chython_layout(mol)
    for ring in mol.sssr:
        assert all(abs(hypot(plane[a][0] - plane[b][0], plane[a][1] - plane[b][1]) - L) < .05 * L
                   for a, b in zip(ring, ring[1:] + ring[:1]))


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
def test_a_spiro_hub_of_five_stretches_its_ring_bonds_for_room():
    # PubChem 14953932: the spirophosphorane's four ring bonds at P reach 1.75 L, the gaps between its halves
    # over 60 degrees
    mol = smiles('CN1C(=S)N(P2(OC(C(O2)(C(F)(F)F)C(F)(F)F)(C(F)(F)F)C(F)(F)F)(OC1(C(F)(F)F)C(F)(F)F)N(C3CCCCC3)'
                 'C4CCCCC4)C')
    plane = chython_layout(mol)
    p = next(x for x in mol if len(mol.neighbors_of(x)) == 5)
    px, py = plane[p]
    ring = [y for y in mol.neighbors_of(p) if mol.in_ring_of(y)]
    assert all(abs(hypot(plane[y][0] - px, plane[y][1] - py) - 1.75 * L) < .05 * L for y in ring)
    turns = sorted((atan2(plane[y][1] - py, plane[y][0] - px), y) for y in mol.neighbors_of(p))
    rings = {frozenset(r) & {y for _, y in turns} for r in mol.sssr if p in r}
    assert all((b - a) % (2 * pi) > radians(60) for (a, y), (b, z) in zip(turns, turns[1:] + turns[:1])
               if frozenset((y, z)) not in rings)
    assert not sum(defects(mol, rescaled(mol, plane)).values())


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
def test_two_bonded_crosses_crowding_each_other_stretch_their_bond():
    # PubChem 19087371: C7 and C11 each square, the bond between them pushed out to 2 L
    mol = smiles('C1CCC(CC1)C(C(=O)O)(C(C2CCCCC2)(C(=O)O)OC3CCCCC3)OC4CCCCC4')
    plane = chython_layout(mol)
    assert hypot(plane[7][0] - plane[11][0], plane[7][1] - plane[11][1]) > 1.9 * L
    assert all(abs(g - 90) < 5 for x in (7, 11) for g in _gaps(mol, plane, x))


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
def test_two_pentagons_on_a_hexagon_atom_share_its_angle():
    # PubChem 25242259, a holostane: C47 holds ring C, the D ring and the lactone; neither pentagon squashed
    mol = smiles('C[C@@H]1[C@H]([C@@H]([C@H](O[C@H]1O[C@H]2[C@@H]([C@H](O[C@H]([C@@H]2O)O[C@@H]3[C@H](O[C@@H]'
                 '([C@@H]([C@H]3O)O)O[C@@H]4[C@H]([C@@H](CO[C@H]4OC5CC[C@]6(C(C5(C)C)CCC7C6=C[C@@H](C89[C@]7(CC'
                 '[C@@]8([C@](OC9=O)(C)CCCC(C)(C)O)O)C)O)C)OS(=O)(=O)O)O)C)C)O)CO)O)OC')
    plane = chython_layout(mol)
    for ring in ((47, 48, 49, 50, 51), (47, 51, 52, 53, 54)):
        assert all(80 < t < 130 for t in _ring_angles(plane, ring)), ring


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
def test_a_bend_stretches_the_chain_bond_not_the_ring_attachment():
    # PubChem 44348845: benzyl C18 bends 60 degrees with C17-C18 pushed out, C18-C19 into the aryl left at L
    mol = smiles('CN(C1=CC=CC=C1N2CCN(CC2)C(=O)[C@@H](CC3=CC=C(C=C3)Cl)NC(=O)[C@H]4CC5=CC=CC=C5CN4)S(=O)(=O)C')
    plane = chython_layout(mol)
    assert abs(hypot(plane[18][0] - plane[19][0], plane[18][1] - plane[19][1]) - L) < .05 * L
    assert hypot(plane[17][0] - plane[18][0], plane[17][1] - plane[18][1]) > 1.2 * L


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
def test_ortho_branches_one_bond_apart_stretch_their_smaller_stem():
    # PubChem 44348845: C1 on N2 and C14 on N9 close a hexagon with no bond; C3-N2 stretches, N2 on its thirds
    mol = smiles('CN(C1=CC=CC=C1N2CCN(CC2)C(=O)[C@@H](CC3=CC=C(C=C3)Cl)NC(=O)[C@H]4CC5=CC=CC=C5CN4)S(=O)(=O)C')
    plane = chython_layout(mol)
    assert hypot(plane[1][0] - plane[14][0], plane[1][1] - plane[14][1]) > 1.2 * L
    assert hypot(plane[2][0] - plane[3][0], plane[2][1] - plane[3][1]) > 1.2 * L
    assert all(abs(g - 120) < 5 for g in _gaps(mol, plane, 2))


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
def test_a_traded_cross_takes_its_shortest_clear_square():
    # the silicon of PubChem 13529450 squares on its ring bond; of the clear squares, none reaching 3 L is kept
    mol = smiles('CC(C)[Si]1(N[Si](N1[Si](N(C)[Si](C)(C)C)(N(C)[Si](C)(C)C)F)(C(C)C)C(C)C)C(C)C')
    plane = chython_layout(mol)
    x = 8
    assert all(abs(g - 90) < 5 for g in _gaps(mol, plane, x))
    assert max(hypot(plane[y][0] - plane[x][0], plane[y][1] - plane[x][1]) for y in mol.neighbors_of(x)) < 2.1 * L


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
@pytest.mark.parametrize('s, x', [(r'CC(C)/C(=C(\C)/C(F)(F)F)/C(C(F)(F)F)(F)F', 7),              # PubChem 22665005
                                  ('CCCCCCCC/C(=C/S(=O)(=O)O)/S(=O)(=O)O', 11)])                # PubChem 20387168
def test_a_crowded_cx3_is_pushed_out_square(s, x):
    # the CF3 and the SO3H find no square room at a unit bond; they square pushed out along it
    mol = smiles(s)
    plane = chython_layout(mol)
    assert all(abs(g - 90) < 5 for g in _gaps(mol, plane, x))


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
@pytest.mark.parametrize('s', ['COC1=C(C=CC(=C1)C=CC2C3=CCC4C(C3CC5C2(C(=O)N(C5=O)NC6=CC=C(C=C6)F)C7=CC=C(C=C7)Cl)'
                               'C(=O)N(C4=O)C8=CC=CC(=C8)C(=O)O)O',                             # PubChem 3423984
                               'CCOP(=O)(COC1=C(C2=C(C=C1)C3=C(P2(=O)OC)C(=C(C=C3)OCP(=O)(OCC)OCC)C4=C(C(=CC(=C4)C(C'
                               ')(C)C)C)CP(=O)(OC)OC)C5=C(C(=CC(=C5)C(C)(C)C)C)CP(=O)(OC)OC)OCC'])  # PubChem 10820014
def test_no_block_is_placed_through_a_ring(s):
    # a crowded ring block leaves by its outside slot, stretched, never across the ring it hangs on
    mol = smiles(s)
    assert defects(mol, chython_layout(mol))['aob'] == 0


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
def test_crowded_crosses_square_with_their_arms_stretched():
    # PubChem 13196961: perfluoro enamine; the CF2 crosses square only once their turned arms stretch
    mol = smiles(r'C(=C(/N(C(C(F)(F)F)(F)F)C(C(F)(F)F)(F)F)\F)(\C(C(C(F)(F)F)(F)F)(F)F)/C(F)(F)F')
    plane = chython_layout(mol)
    for x in (11, 12, 19, 20, 21):
        assert all(abs(g - 90) < 5 for g in _gaps(mol, plane, x)), x


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
@pytest.mark.parametrize('s', [
    'CCOP(=O)(C1=CC=C(C=C1)C[C@@](C)(C(=O)N[C@@H](CC(=O)N)C(=O)NCCCC2=CC=CC3=CC=CC=C32)NC(=O)[C@@H]'
    '(CC4=CC=C(C=C4)CP(=O)(OC(C)(C)C)OC(C)(C)C)CC(=O)OC(C)(C)C)OCC',                         # PubChem 15974362
    'CCOC(=O)N1C=C(C2=C(C=CC=C21)C3=C(C(=CC=C3)I)OC(=O)/C(=C/C)/C4=CC=CC=C4)C5=CN=C(O5)CO[Si]'
    '(C6=CC=CC=C6)(C7=CC=CC=C7)C(C)(C)C'])                                                   # PubChem 11018374
def test_every_tert_butyl_is_square(s):
    # a stretched phosphonate leaves its tert-butyls room; a tert-butyl between two phenyls on Si is pushed out
    mol = smiles(s)
    plane = chython_layout(mol)
    for x in mol:
        if len(mol.neighbors_of(x)) == 4 and not mol.in_ring_of(x) and mol.atom(x).atomic_symbol == 'C' and \
                sum(len(mol.neighbors_of(y)) == 1 for y in mol.neighbors_of(x)) == 3:
            assert all(abs(g - 90) < 5 for g in _gaps(mol, plane, x)), x


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
def test_a_phosphate_chain_is_a_straight_run_of_crosses():
    # PubChem 5957: ATP, stereo left out
    mol = smiles('Nc1ncnc2c1ncn2C1OC(COP(=O)(O)OP(=O)(O)OP(=O)(O)O)C(O)C1O')
    plane = chython_layout(mol)
    ps = [x for x in mol if mol.atom(x).atomic_symbol == 'P']
    bridges = [x for x in mol if mol.atom(x).atomic_symbol == 'O' and
               sum(mol.atom(y).atomic_symbol == 'P' for y in mol.neighbors_of(x)) == 2]
    assert len(ps) == 3 and len(bridges) == 2
    for x in ps:
        assert all(abs(g - 90) < 5 for g in _gaps(mol, plane, x)), (x, _gaps(mol, plane, x))
    for x in bridges:
        assert abs(max(_gaps(mol, plane, x)) - 180) < 5
    d = defects(mol, plane)
    assert (d['ovl61'], d['crs'], d['aob']) == (0, 0, 0)


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
@pytest.mark.parametrize('s', [
    'C1C2=CC(=CC=C2)CC3=CC(=CC=C3)CC4=CC=CC(=C4)CC5=CC=CC1=C5',      # the corpus's hard:calix[4]arene core
])
def test_a_small_ring_of_a_macrocycle_is_convex(s):
    mol = smiles(s)
    plane = chython_layout(mol)
    assert all(_convex(plane, r) for r in mol.sssr if len(r) <= 8)


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
@pytest.mark.parametrize('s', [PORPHINE, CALIXARENE])
def test_a_macrocycle_of_rigid_units_draws_as_one_closed_polygon(s):
    mol = smiles(s)
    rank = mol.canonical_order()
    block = max(ring_blocks(mol), key=lambda b: len(b[0]))
    assert macro_block(mol, *block, rank) is not None
    plane = chython_layout(mol)
    assert all(_polygon_rmsd([plane[x] for x in r]) < 1e-6 for r in mol.sssr if len(r) <= 8)
    d = defects(mol, plane)
    assert (d['ovl61'], d['crs'], d['aob']) == (0, 0, 0)


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
def test_porphine_is_fourfold():
    mol = smiles(PORPHINE)
    plane = chython_layout(mol)
    ns = [x for x in mol if mol.atom(x).atomic_symbol == 'N']
    cx, cy = sum(plane[x][0] for x in ns) / 4, sum(plane[x][1] for x in ns) / 4
    r = [hypot(plane[x][0] - cx, plane[x][1] - cy) for x in ns]
    assert max(r) - min(r) < 1e-6


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
def test_a_ring_branch_bisects_the_outside_angle_of_its_atom():
    mol = smiles(CALIXARENE)
    plane = chython_layout(mol)
    rings = [set(r) for r in mol.sssr]
    seen = 0
    for x in mol:
        nb = mol.neighbors_of(x)
        if len(nb) != 3:
            continue
        branch = [y for y in nb if not any(x in r and y in r for r in rings)]
        if len(branch) != 1:
            continue
        a, b = (y for y in nb if y != branch[0])
        t = degrees(atan2(plane[branch[0]][1] - plane[x][1], plane[branch[0]][0] - plane[x][0]))
        ga, gb = (abs((t - degrees(atan2(plane[y][1] - plane[x][1], plane[y][0] - plane[x][0])) + 180) % 360 - 180)
                  for y in (a, b))
        assert ga == pytest.approx(gb, abs=1)
        seen += 1
    assert seen == 8


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
@pytest.mark.parametrize('s, x', [('CC1=CC=C(C=C1)C(C2C(=O)O[Sn](OC3=CC=CC=C3C=N2)(C4=CC=CC=C4)C5=CC=CC=C5)O',
                                   29)])                                                        # PubChem 16713265
def test_a_pendant_hexagon_turns_onto_its_bond(s, x):
    # the phenyl on the tin turns about its ipso, the larger side held, until the tin bond bisects it
    mol = smiles(s)
    plane = chython_layout(mol)
    assert all(abs(g - 120) < 2 for g in _gaps(mol, plane, x))


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
def test_a_saturated_bridge_of_two_units_holds_120_degrees():
    mol = smiles('C1C2=CC(=CC=C2)CC3=CC(=CC=C3)CC4=CC=CC(=C4)CC5=CC=CC1=C5')        # the calix[4]arene core
    plane = chython_layout(mol)
    bridges = [x for x in mol if not any(x in r for r in mol.sssr if len(r) <= 8)]
    assert len(bridges) == 4
    for x in bridges:
        assert abs(max(_gaps(mol, plane, x)) - 240) < 2, _gaps(mol, plane, x)
    d = defects(mol, plane)
    assert (d['ovl61'], d['crs'], d['aob']) == (0, 0, 0)


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
@pytest.mark.parametrize('s', [
    'C1COCCOCCOCCOCCO1',                                    # PubChem 36336: 15-crown-5
    'C1CCCCCCCCCCCCCC1',                                    # PubChem 9200: cyclopentadecane
    'C1CCCCCCCCCCCC(C)CC1',                                 # methylcyclopentadecane: the spike turns inward
])
def test_an_odd_macrocycle_has_no_spike(s):
    mol = smiles(s)
    plane = chython_layout(mol)
    ring = max(mol.sssr, key=len)
    n = len(ring)
    for i in range(n):
        a, b, c = plane[ring[i - 1]], plane[ring[i]], plane[ring[(i + 1) % n]]
        t = degrees(atan2(a[1] - b[1], a[0] - b[0]) - atan2(c[1] - b[1], c[0] - b[0])) % 360
        assert 100 < t < 260, (ring[i], t)
    lengths = [hypot(plane[ring[i]][0] - plane[ring[i - 1]][0], plane[ring[i]][1] - plane[ring[i - 1]][1])
               for i in range(n)]
    assert max(lengths) - min(lengths) < .01 * L
    d = defects(mol, plane)
    assert (d['ovl61'], d['crs'], d['aob']) == (0, 0, 0)


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
@pytest.mark.parametrize('s', [
    'C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2O)CCC4=CC(=O)CC[C@]34C',      # PubChem 6013: testosterone
    'C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C',    # PubChem 5997
])
def test_an_angular_methyl_stays_outside_its_rings(s):
    mol = smiles(s)
    plane = chython_layout(mol)
    for r in mol.sssr:
        poly = [plane[x] for x in r]
        assert not [x for x in mol if x not in r and _inside(plane[x], poly)], r


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
def test_a_five_ring_of_a_macrocycle_unit_is_regular():
    mol = smiles(RIFAMYCIN)
    plane = chython_layout(mol)
    ring = next(r for r in mol.sssr if len(r) == 5)
    assert _polygon_rmsd([plane[x] for x in ring]) < .01


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
def test_two_fused_five_rings_are_both_regular():
    # PubChem 21597303: a dioxolane on the furanose of a 6-5-5 block
    mol = smiles('C[C@@H]([C@@H]1[C@@H]2[C@@H]([C@@H]3[C@H](O2)OC(O3)(C)C)OC(O1)(C)C)O[Si](C)(C)C(C)(C)C')
    plane = chython_layout(mol)
    fives = [r for r in mol.sssr if len(r) == 5]
    assert len(fives) == 2
    for ring in fives:
        assert _polygon_rmsd([plane[x] for x in ring]) < .01, ring


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
def test_a_ring_fused_into_a_macrocycle_is_seated_on_unit_bonds():
    # PubChem 10424223: a dihydrofuran on a cembranoid ring; both macrocycle bonds into it L, their angles even
    mol = smiles('C/C/1=C/CC[C@@H]2[C@@H](O2)CC/C(=C/[C@H]3C(=C(CO3)C)CC1)/C')
    plane = chython_layout(mol)
    for b in mol.bonds():
        assert hypot(plane[b.n][0] - plane[b.m][0], plane[b.n][1] - plane[b.m][1]) == pytest.approx(L, abs=1e-6)
    for x in (12, 19):
        assert min(_gaps(mol, plane, x)) == pytest.approx(120, abs=1)
    for x in (13, 14):
        assert min(_gaps(mol, plane, x)) > 107
    d = defects(mol, plane)
    assert (d['ovl61'], d['crs'], d['aob']) == (0, 0, 0)


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
@pytest.mark.parametrize('s', [
    # PubChem 36314, paclitaxel: the benzoate and the acetate on two neighbours of the oxetane
    'CC1=C2[C@H](C(=O)[C@@]3([C@H](C[C@@H]4[C@]([C@H]3[C@@H]([C@@](C2(C)C)(C[C@@H]1OC(=O)[C@@H]([C@H](C5=CC=CC=C5)'
    'NC(=O)C6=CC=CC=C6)O)O)OC(=O)C7=CC=CC=C7)(CO4)OC(=O)C)O)C)OC(=O)C',
    # PubChem 12635218: a sulfonate between two ortho nitro groups
    'CCS(C)(CC(O)O)OS(=O)(=O)C1=C(C=C(C=C1[N+](=O)[O-])[N+](=O)[O-])[N+](=O)[O-]'])
def test_labelled_atoms_breathe_apart(s):
    mol = smiles(s)
    plane = chython_layout(mol)
    labelled = [x for x in mol if mol.atom(x).atomic_symbol != 'C']
    near = min(hypot(plane[a][0] - plane[b][0], plane[a][1] - plane[b][1])
               for i, a in enumerate(labelled) for b in labelled[i + 1:] if b not in mol.neighbors_of(a))
    assert near >= CRAMMED
    d = defects(mol, plane)
    assert (d['ovl61'], d['crs'], d['aob']) == (0, 0, 0)


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
def test_a_ribbon_macrocycle_is_redrawn_round():
    # rapamycin's 29-ring walked from its pyran is one hexagon high; the redrawn ring is round and clean
    mol = smiles(RAPAMYCIN)
    plane = chython_layout(mol)
    macro = max(mol.sssr, key=len)
    assert roundness([plane[x] for x in macro]) > .5
    d = defects(mol, plane)
    assert not d['ovl61'] and not d['crs'] and not d['aob']
    assert all(hypot(plane[b.n][0] - plane[b.m][0], plane[b.n][1] - plane[b.m][1]) < 1.1 * L for b in mol.bonds())


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
def test_a_fullerene_is_a_sphere_seen_down_a_pentagon():
    mol = smiles(C60)
    rank = mol.canonical_order()
    (atoms, rings), = ring_blocks(mol)
    pts = sphere_block(mol, atoms, rings, rank)
    assert pts is not None and set(pts) == atoms
    plane = chython_layout(mol)
    assert defects(mol, plane)['ovl'] == 0
    assert perm(mol, chython_layout, k=3) == 0


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
def test_a_small_cage_stays_planar():
    mol = smiles('C12C3C4C5C1C6C7C2C8C3C9C4C1C5C6C2C7C8C9C12')                 # dodecahedrane
    (atoms, rings), = ring_blocks(mol)
    assert sphere_block(mol, atoms, rings, mol.canonical_order()) is None
    d = defects(mol, chython_layout(mol))
    assert (d['ovl'], d['crs'], d['aob']) == (0, 0, 0)


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
def test_a_dense_cage_is_a_projected_model():
    mol = smiles('C1C2C3C4CC5C3C1C67C28C49C56C\\1C2C7C(=O)CC8C2C9C(=O)/C1=C/O')          # PubChem 13806692
    atoms, rings = max(ring_blocks(mol), key=lambda b: len(b[0]))
    assert cage_block(mol, atoms, rings, mol.canonical_order()) is not None
    plane = chython_layout(mol)
    ls = [hypot(plane[x][0] - plane[y][0], plane[x][1] - plane[y][1])
          for x in atoms for y in mol.neighbors_of(x) if y in atoms and y > x]
    assert max(ls) < 2 * min(ls)                    # its planar drawing's bonds span 0.55 to 4.1


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
def test_a_strained_cage_is_a_projected_model():
    mol = smiles('CC12CC3C4(C56C1C(=O)C(O5)(C7CC(C8(CC=CC(=O)C8(C7CCC6(C(=O)O4)O)C)O)O)OCC2C(=O)O3)C')  # PubChem 431071
    atoms, _ = max(ring_blocks(mol), key=lambda b: len(b[0]))
    plane = chython_layout(mol)
    ls = [hypot(plane[x][0] - plane[y][0], plane[x][1] - plane[y][1])
          for x in atoms for y in mol.neighbors_of(x) if y in atoms and y > x]
    assert max(ls) < 2 * min(ls)                    # its lattice drawing stretches an ether bridge to 3 bonds


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
def test_a_ring_fused_to_a_nonagon_takes_an_edge_with_room_outside():
    mol = smiles('C1CCCc2ccccc2CCC1')                                                 # PubChem 57720623
    plane = layout2d(mol, engine='chython', force=True)
    ls = {(b.n, b.m): hypot(plane[b.n][0] - plane[b.m][0], plane[b.n][1] - plane[b.m][1]) for b in mol.bonds()}
    median = sorted(ls.values())[len(ls) // 2]
    fused = {x for x in mol if len(list(mol.neighbors_of(x))) == 3}
    assert not any(v > 1.2 * median and fused.intersection(k) for k, v in ls.items())  # the chord sits clear of them


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
@pytest.mark.parametrize('s', [
    'CN1CC[C@]23[C@@H]4[C@H]1CC5=C2C(=C(C=C5)O)O[C@H]3[C@H](C=C4)O',                 # PubChem 5288826: morphine
    'CN1CCC23CCCCC2C1CC4=C3C=C(C=C4)OC',                                            # dextromethorphan, no stereo
])
def test_a_morphinan_takes_its_template(s):
    mol = smiles(s)
    rank = mol.canonical_order()
    block = max(ring_blocks(mol), key=lambda b: len(b[0]))
    assert template_block(mol, *block, rank) is not None
    plane = chython_layout(mol)
    assert sum(_polygon_rmsd([plane[x] for x in r]) < 1e-6 for r in mol.sssr if len(r) == 6) == 3
    d = defects(mol, plane)
    assert (d['ovl61'], d['crs'], d['aob']) == (0, 1, 0)                          # the bridge crosses C8-C14
    assert perm(mol, chython_layout, k=5) == 0


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
@pytest.mark.parametrize('s', [
    'C1=CC2=CC3=CC=C(N3)C=C4C=CC(=N4)C=C5C=CC(=N5)C=C1N2',                          # PubChem 66868: porphine
    'CC1=C(C2=CC3=NC(=CC4=C(C(=C([N-]4)C=C5C(=C(C(=N5)C=C1[N-]2)C)C=C)C)C=C)C(=C3CCC(=O)O)C)CCC(=O)O.[Fe+2]',  # heme b
    'CCC1=C(C2=NC1=CC3=C(C4=C([N-]3)C(=C5[C@H]([C@@H](C(=N5)C=C6C(=C(C(=C2)[N-]6)C=C)C)C)CCC(=O)OC/C=C(\\C)/CCC'
    '[C@H](C)CCC[C@H](C)CCCC(C)C)[C@H](C4=O)C(=O)OC)C)C.[Mg+2]',                  # chlorophyll a
])
def test_a_porphyrin_is_the_square_of_its_template(s):
    mol = smiles(s)
    rank = mol.canonical_order()
    block = max(ring_blocks(mol), key=lambda b: len(b[0]))
    assert template_block(mol, *block, rank) is not None
    plane = chython_layout(mol)
    fives = [r for r in mol.sssr if len(r) == 5]
    assert all(_polygon_rmsd([plane[x] for x in r]) < 1e-6 for r in fives if any(mol.atom(x).atomic_symbol == 'N'
                                                                                   for x in r))
    ns = [plane[x] for x in block[0] if mol.atom(x).atomic_symbol == 'N']
    cx, cy = sum(p[0] for p in ns) / 4, sum(p[1] for p in ns) / 4
    d = sorted(hypot(p[0] - cx, p[1] - cy) for p in ns)
    assert d[-1] - d[0] < 1e-6
    d = defects(mol, plane)
    assert (d['ovl61'], d['crs'], d['aob']) == (0, 0, 0)


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
@pytest.mark.parametrize('s', [
    _TAXOL,                                                                             # PubChem 36314: paclitaxel
    _TAXOL.replace('NC(=O)C6=CC=CC=C6', 'NC(=O)OC(C)(C)C')[:-7] + 'O',                  # PubChem 148124: docetaxel
    _TAXOL.replace('OC(=O)[C@@H]([C@H](C5=CC=CC=C5)NC(=O)C6=CC=CC=C6)O', 'O'),          # PubChem 65366: baccatin III
])
def test_a_taxane_keeps_ring_a_regular_and_its_gem_dimethyl_in_ring_b(s):
    mol = smiles(s)
    rank = mol.canonical_order()
    block = max(ring_blocks(mol), key=lambda b: len(b[0]))
    assert template_block(mol, *block, rank) is not None
    plane = chython_layout(mol)
    c15 = next(x for x in block[0] if sum(len(mol.neighbors_of(y)) == 1 for y in mol.neighbors_of(x)) == 2)
    ring_a = next(r for r in mol.sssr if len(r) == 6 and c15 in r)
    poly = [plane[x] for x in ring_a]
    assert _polygon_rmsd(poly) < 1e-6
    assert not any(_inside(plane[x], poly) for x in mol if x not in ring_a)
    d = defects(mol, plane)
    assert (d['ovl61'], d['crs'], d['aob']) == (0, 0, 0)


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
@pytest.mark.parametrize('s', [
    'C[C@@H]1CC[C@H]2[C@H](C(=O)O[C@H]3[C@@]24[C@H]1CC[C@](O3)(OO4)C)C',               # PubChem 68827: artemisinin
    'C[C@@H]1CC[C@H]2[C@H]([C@@H](O[C@H]3[C@@]24[C@H]1CC[C@](O3)(OO4)C)OC)C',          # PubChem 68911: artemether
])
def test_an_artemisinin_draws_its_peroxide_inside_the_cycloheptane(s):
    mol = smiles(s)
    rank = mol.canonical_order()
    block = max(ring_blocks(mol), key=lambda b: len(b[0]))
    assert template_block(mol, *block, rank) is not None
    plane = chython_layout(mol)
    carbocycles = [r for r in mol.sssr if len(r) == 6 and sum(mol.atom(x).atomic_symbol == 'O' for x in r) < 2]
    assert len(carbocycles) == 2 and all(_polygon_rmsd([plane[x] for x in r]) < 1e-6 for r in carbocycles)
    d = defects(mol, plane)
    assert (d['ovl61'], d['crs'], d['aob']) == (0, 0, 0)
    assert shape(mol, plane)['inr'] == 0
    assert perm(mol, chython_layout, k=5) == 0


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
def test_strychnine_is_three_hexagons_and_a_pentagon_on_the_honeycomb():
    # PubChem 441071: strychnine
    mol = smiles('C1CN2CC3=CCO[C@H]4CC(=O)N5[C@H]6[C@H]4[C@H]3C[C@H]2[C@@]61C7=CC=CC=C75')
    rank = mol.canonical_order()
    block = max(ring_blocks(mol), key=lambda b: len(b[0]))
    assert template_block(mol, *block, rank) is not None
    plane = chython_layout(mol)
    six = [r for r in mol.sssr if len(r) == 6 and 3 not in r]                     # benzene, lactam, cyclohexane
    assert len(six) == 3 and all(_polygon_rmsd([plane[x] for x in r]) < 1e-6 for r in six)
    assert plane[20][1] == pytest.approx(plane[19][1])                             # the indoline's chord is level
    d = defects(mol, plane)
    assert (d['ovl61'], d['crs'], d['aob']) == (0, 0, 0)
    assert perm(mol, chython_layout, k=5) == 0


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
def test_two_aryls_on_vicinal_atoms_stagger_clear_of_each_other():
    # PubChem 3598851: phenyl and dimethoxyphenol hang side by side from a succinimide's ring fusion
    mol = smiles('COC1=CC(=CC(=C1O)OC)C2C3=CCC4C(C3CC5C2(C(=O)N(C5=O)C6=CC(=CC=C6)Cl)C7=CC=CC=C7)C(=O)N(C4=O)'
                 'CCCCCC(=O)O')
    plane = chython_layout(mol)
    d = defects(mol, plane)
    assert (d['ovl61'], d['crs'], d['aob'], d['lbl']) == (0, 0, 0, 0)
    imide = next(r for r in mol.sssr if len(r) == 5 and 24 in r)                  # 20 21 22 N24 25
    sides = [hypot(plane[a][0] - plane[b][0], plane[a][1] - plane[b][1]) for a, b in zip(imide, imide[1:] + imide[:1])]
    assert max(sides) - min(sides) < 1e-6
    assert perm(mol, chython_layout, k=5) == 0


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
def test_a_branch_leaves_its_siblings_their_slots():
    # PubChem 14240641: the isobutyl branch, grown first, leaves the carbamate oxygen its slot
    mol = smiles('CC(C)CC(CC(C=C)O)OC(=O)N(C)C')
    plane = chython_layout(mol)
    for x in mol:
        nb = mol.neighbors_of(x)
        if len(nb) == 3:
            px, py = plane[x]
            a = sorted(atan2(plane[y][1] - py, plane[y][0] - px) for y in nb)
            assert min((v - u) % (2 * pi) for u, v in zip(a, a[1:] + a[:1])) > 119 * pi / 180, x


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
def test_a_ring_s_branches_leave_each_other_their_slots():
    # PubChem 14441915: the imine's ring system grows first; the fluorine keeps its bisector on the aryl
    mol = smiles('C=CCOC1=C(C=C(C(=C1)N=C2N3CCCCC3CCS2)F)Cl')
    plane = chython_layout(mol)
    f = next(x for x in mol if mol.atom(x).atomic_symbol == 'F')
    c, = mol.neighbors_of(f)
    px, py = plane[c]
    a = [atan2(plane[y][1] - py, plane[y][0] - px) for y in mol.neighbors_of(c)]
    gaps = sorted((v - u) % (2 * pi) for u, v in zip(sorted(a), sorted(a)[1:] + sorted(a)[:1]))
    assert degrees(gaps[0]) > 110


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
@pytest.mark.parametrize('s', ['C1=CC2CC1C1C2C(=O)OC1=O', 'C1CC2C3CCC(C3)C2C1', 'C1CC2CC1C1CCCCC21'])
def test_a_one_atom_bridge_spans_a_regular_hexagon(s):
    # nadic anhydride, tricyclo[5.2.1.0^2,6]decane, a benzo-fused norbornane: no shape draws these systems
    mol = smiles(s)
    plane = chython_layout(mol)
    ls = [hypot(plane[b.n][0] - plane[b.m][0], plane[b.n][1] - plane[b.m][1]) for b in mol.bonds()
          if mol.in_ring_of(b.n) and mol.in_ring_of(b.m)]
    assert max(ls) / min(ls) < 1.05
    d = defects(mol, plane)
    assert (d['ovl61'], d['crs'], d['aob']) == (0, 0, 0)


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
@pytest.mark.parametrize('s', ['C1C2=CC=CC=C2C3=CC=CC=C31', 'N1C2=CC=CC=C2C3=CC=CC=C31', 'O1C2=CC=CC=C2C3=CC=CC=C31'])
def test_a_mirror_symmetric_ring_system_stands_on_its_axis(s):
    # fluorene, carbazole, dibenzofuran: atom 1 is on the mirror line, which is drawn vertical
    mol = smiles(s)
    plane = chython_layout(mol)
    cx = sum(x for x, _ in plane.values()) / len(plane)
    assert plane[1][0] == pytest.approx(cx, abs=.01 * L)
    orbit = mol.automorphism_orbits()
    for x, (px, py) in plane.items():
        assert any(orbit[y] == orbit[x] and hypot(2 * cx - px - qx, py - qy) < .02 * L
                   for y, (qx, qy) in plane.items())


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
@pytest.mark.parametrize('s', ['C1CCCCCCC1', 'C1=CC=CC=CC=C1'])
def test_an_octagon_stands_on_an_edge(s):
    # cyclooctane, cyclooctatetraene: two bonds vertical, two horizontal
    mol = smiles(s)
    plane = chython_layout(mol)
    dx = [abs(plane[b.n][0] - plane[b.m][0]) for b in mol.bonds()]
    dy = [abs(plane[b.n][1] - plane[b.m][1]) for b in mol.bonds()]
    assert sum(d < .01 * L for d in dx) == sum(d < .01 * L for d in dy) == 2


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
def test_the_rings_of_a_spiro_atom_face_opposite_ways():
    # spiro[4.5]decane: the bisectors of the spiro atom's two ring bonds on each side are one line
    mol = smiles('C1CCC2(C1)CCCCC2')
    plane = chython_layout(mol)
    x = next(a for a in plane if len(mol.neighbors_of(a)) == 4)
    sides = []
    for r in mol.sssr:
        u, v = (sum(plane[y][k] - plane[x][k] for y in mol.neighbors_of(x) if y in r) for k in (0, 1))
        sides.append(atan2(v, u))
    assert abs(abs(sides[0] - sides[1]) - pi) < 1e-3


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
@pytest.mark.parametrize('s', ['Cl[Re](Cl)(Cl)(Cl)Cl', 'ClP(Cl)(Cl)(Cl)Cl'])
def test_five_terminal_ligands_are_a_regular_star(s):
    # rhenium(V) chloride, phosphorus pentachloride: 72 degrees apart, one straight up
    mol = smiles(s)
    plane = chython_layout(mol)
    x = next(a for a in plane if len(mol.neighbors_of(a)) == 5)
    angs = sorted(degrees(atan2(plane[y][1] - plane[x][1], plane[y][0] - plane[x][0])) % 360
                  for y in mol.neighbors_of(x))
    assert [b - a for a, b in zip(angs, angs[1:])] == pytest.approx([72] * 4, abs=.1)
    assert any(abs(a - 90) < .1 for a in angs)
    assert all(hypot(plane[y][0] - plane[x][0], plane[y][1] - plane[x][1]) == pytest.approx(L)
               for y in mol.neighbors_of(x))


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
@pytest.mark.parametrize('s', ['C1=CC2=C3C(=C1)C=CC3=CC=C2', 'C1CC2=CC=CC3=C2C1=CC=C3'])
def test_a_peri_fused_pentagon_is_an_arc_of_unit_bonds_on_top(s):
    # acenaphthylene, acenaphthene: the naphthalene holds, the pentagon's free bond is horizontal above it
    mol = smiles(s)
    plane = chython_layout(mol)
    ring = next(r for r in mol.sssr if len(r) == 5)
    for a, b in zip(ring, ring[1:] + ring[:1]):
        assert hypot(plane[a][0] - plane[b][0], plane[a][1] - plane[b][1]) == pytest.approx(L, rel=.02)
    a, b = [x for x in ring if len([r for r in mol.sssr if x in r]) == 1]
    cy = sum(y for _, y in plane.values()) / len(plane)
    assert plane[a][1] == pytest.approx(plane[b][1], abs=.01 * L) and plane[a][1] > cy


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
@pytest.mark.parametrize('s', ['C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C',
                               'C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2O)CCC4=CC(=O)CC[C@]34C',
                               'C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2O)CCC4=C3C=CC(=C4)O'])
def test_a_steroid_lies_horizontal_ring_a_on_the_left(s):
    # cholesterol, testosterone, estradiol: ring A left of ring D, the A-D axis nearer horizontal than 30 degrees
    mol = smiles(s)
    plane = chython_layout(mol)
    rings = sorted(mol.sssr, key=len)
    d = rings[0]
    a = next(r for r in rings if not set(r) & set(d) and
             any(mol.atom(y).atomic_symbol == 'O' for x in r for y in mol.neighbors_of(x) if not mol.in_ring_of(y)))
    centre = lambda r: (sum(plane[x][0] for x in r) / len(r), sum(plane[x][1] for x in r) / len(r))
    (ax, ay), (dx, dy) = centre(a), centre(d)
    assert ax < dx and abs(dy - ay) < tan(radians(30)) * (dx - ax)


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
@pytest.mark.parametrize('s', [
    'CC(C)(C)CC(C)(C)OP(=O)(OC)OC(C)(C)CC(C)(C)C',          # PubChem 14324920
    'CC1=CN(C(=O)NC1=O)[C@H]2C[C@@H]([C@H](O2)CO)O'])      # thymidine, PubChem 5789
def test_a_long_drawing_lies_horizontal(s):
    mol = smiles(s)
    plane = chython_layout(mol)
    n = len(plane)
    cx = sum(x for x, _ in plane.values()) / n
    cy = sum(y for _, y in plane.values()) / n
    xx = sum((x - cx) ** 2 for x, _ in plane.values())
    yy = sum((y - cy) ** 2 for _, y in plane.values())
    xy = sum((x - cx) * (y - cy) for x, y in plane.values())
    assert abs(degrees(atan2(2 * xy, xx - yy) / 2)) < 15


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
def test_a_corrin_is_mirror_symmetric_about_its_direct_bond():
    # cyanocobalamin, PubChem 5311498: the corrin's axis runs through its A-D bond and the opposite meso carbon
    mol = smiles(
           'CC1=CC2=C(C=C1C)N(C=N2)[C@@H]3[C@@H]([C@@H]([C@H](O3)CO)OP(=O)([O-])O[C@H](C)CNC(=O)CC[C@@]4([C@H](['
           'C@@H]5[C@]6([C@@]([C@@H](/C(=C(/C7=N/C(=C\C8=N/C(=C(\C4=N5)/C)/[C@H](C8(C)C)CCC(=O)N)/[C@H]([C@]7(C)'
           'CC(=O)N)CCC(=O)N)\C)/[N-]6)CCC(=O)N)(C)CC(=O)N)C)CC(=O)N)C)O.[C-]#N.[Co+3]')
    plane = chython_layout(mol)
    core = next(a for a, rs in ring_blocks(mol) if len(a) == 23)
    ring = {x for r in mol.sssr if len(r) == 5 and set(r) <= core for x in r}
    n, m = next((b.n, b.m) for b in mol.bonds() if b.n in ring and b.m in ring and
                not any(b.n in r and b.m in r for r in mol.sssr if len(r) == 5))
    ox, oy = (plane[n][0] + plane[m][0]) / 2, (plane[n][1] + plane[m][1]) / 2
    far = max(core - ring, key=lambda x: hypot(plane[x][0] - ox, plane[x][1] - oy))
    ux, uy = plane[far][0] - ox, plane[far][1] - oy
    d = hypot(ux, uy)
    ux, uy = ux / d, uy / d
    for x in core:
        px, py = plane[x]
        t = 2 * ((px - ox) * ux + (py - oy) * uy)
        rx, ry = ox + t * ux - (px - ox), oy + t * uy - (py - oy)
        assert min(hypot(rx - plane[y][0], ry - plane[y][1]) for y in core) < .05 * L


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
def test_a_relaxed_drawing_taller_than_wide_turns_level():
    # NAD+, PubChem 5892: regularizing the furanoses swings the diphosphate chain off the lattice pose
    plane = chython_layout(smiles('C1=CC(=C[N+](=C1)[C@H]2[C@@H]([C@@H]([C@H](O2)COP(=O)([O-])OP(=O)(O)OC[C@@H]3'
                                  '[C@H]([C@H]([C@@H](O3)N4C=NC5=C(N=CN=C54)N)O)O)O)O)C(=O)N'))
    xs = [p[0] for p in plane.values()]
    ys = [p[1] for p in plane.values()]
    assert max(xs) - min(xs) > 1.5 * (max(ys) - min(ys))


# --- relax outcomes --------------------------------------------------------------------------------------

def _rmsd(plane, ring):
    return _polygon_rmsd([plane[x] for x in ring])


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
def test_a_substituent_turns_with_its_ring_atom():
    mol = smiles('Cn1ccc2ccccc21')
    plane = chython_layout(mol)
    methyl = next(x for x in mol if not mol.neighbors_of(x)[1:] and mol.atom(x).atomic_symbol == 'C')
    n = mol.neighbors_of(methyl)[0]
    assert hypot(plane[n][0] - plane[methyl][0], plane[n][1] - plane[methyl][1]) == pytest.approx(L)
    assert defects(mol, plane)['ovl61'] == 0


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
def test_a_fused_pentagon_turns_its_squeezed_branch_onto_the_exit():
    # PubChem 27591152: an N-alkyl indole-2-carboxamide; the pyrrole ring held off a pentagon by its N branch
    mol = smiles('CC1(COCN1C(=O)C2=CC3=CC=CC=C3N2C[C@@H]4CCCN4CC5=CC(=CC=C5)C(F)(F)F)C')
    plane = chython_layout(mol)
    ring = (9, 10, 11, 16, 17)
    ls = [hypot(plane[a][0] - plane[b][0], plane[a][1] - plane[b][1]) for a, b in zip(ring, ring[1:] + ring[:1])]
    assert ls == pytest.approx([L] * 5, rel=.03)


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
def test_a_clean_layout_logs_nothing():
    log = []
    chython_layout(smiles('c1ccc2[nH]ccc2c1'), log=log)
    assert log == []


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
def test_a_fused_pentagon_carrying_a_spiro_ring_is_regular():
    mol = smiles('C1CCC2(C1)CCC1=CC=CC=C12')                # spiro[cyclopentane-1,1'-indane]
    plane = chython_layout(mol)
    five = next(r for r in mol.sssr if len(r) == 5 and any(len(mol.neighbors_of(x)) == 3 for x in r))
    assert _rmsd(plane, five) < 1e-6


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
@pytest.mark.parametrize('s', ['CC1(C)CCCCC1', 'CC1(C)CCCC1', 'OC1(CCN(C)CC1)C1=CC=C(Cl)C=C1'])
def test_two_branches_of_a_ring_atom_split_its_outside_gap_in_thirds(s):
    mol = smiles(s)
    plane = chython_layout(mol)
    x = next(x for x in mol if len(mol.neighbors_of(x)) == 4)
    px, py = plane[x]
    a = sorted(atan2(plane[y][1] - py, plane[y][0] - px) for y in mol.neighbors_of(x))
    gaps = sorted((v - u) % (2 * pi) for u, v in zip(a, a[1:] + a[:1]))
    assert gaps[2] - gaps[0] < 1e-6 and gaps[3] < pi


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
def test_a_ring_straight_on_the_lattice_turns_its_branch_outward():
    # PubChem 49276653: the diazepane's N is straight on the lattice, its chain outside; regular, the chain stays out
    mol = smiles('CC(C)CCOCCN1CCCN2CCCC2C1')
    plane = chython_layout(mol)
    seven = next(r for r in mol.sssr if len(r) == 7)
    assert _rmsd(plane, seven) < 1e-6
    assert defects(mol, plane)['ovl61'] == 0


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
@pytest.mark.parametrize('s', [
    'COC1=CC(=CC(=C1OC)OC)C(=O)N(CC2=CC=CC=C2)C3=NC4=C(S3)C=CC=C4SC',                 # PubChem 41679762
    'CC(C)(C)OC(=O)N[C@H](CC1=CC=C(C=C1)[N+](=O)[O-])CC(=O)[O-]',                     # PubChem 7021580
    # PubChem 44348845
    'CN(C1=CC=CC=C1N2CCN(CC2)C(=O)[C@@H](CC3=CC=C(C=C3)Cl)NC(=O)[C@H]4CC5=CC=CC=C5CN4)S(=O)(=O)C',
])
def test_a_straight_chain_atom_bends_where_it_has_room(s):
    # 60 degrees either way where it fits, 30 where only that does: no chain atom left straight
    mol = smiles(s)
    plane = chython_layout(mol)
    assert shape(mol, plane)['ang'] < 45
    d = defects(mol, plane)
    assert (d['ovl61'], d['crs'], d['aob']) == (0, 0, 0)


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
def test_a_straight_atom_stretches_a_bond_to_bend_120_degrees():
    # PubChem 14240641: the isobutyl CH2 has no room for 120 degrees at one bond, only at a stretched one
    mol = smiles('CC(C)CC(CC(C=C)O)OC(=O)N(C)C')
    plane = chython_layout(mol)
    (ax, ay), (bx, by), (cx, cy) = plane[2], plane[4], plane[5]
    turn = abs(atan2((ax - bx) * (cy - by) - (ay - by) * (cx - bx), (ax - bx) * (cx - bx) + (ay - by) * (cy - by)))
    assert turn == pytest.approx(2 * pi / 3, abs=.02)
    d = defects(mol, plane)
    assert (d['ovl61'], d['crs'], d['aob']) == (0, 0, 0)


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
@pytest.mark.parametrize('s', ['C1CCC(C1)C1=CC=CC=C1','C1CCCC(CC1)C1=CC=CC=C1', 'C1=NN=NN1C1CCCC1'])
def test_a_ring_fused_to_nothing_is_regular_at_one_bond(s):
    mol = smiles(s)
    plane = chython_layout(mol)
    ls = [hypot(plane[b.n][0] - plane[b.m][0], plane[b.n][1] - plane[b.m][1]) for b in mol.bonds()]
    assert max(ls) - min(ls) < 1e-6


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
def test_an_ester_carbon_crowded_to_one_side_spreads_its_bonds():
    # PubChem 10033120: a peracetylated lactoside; two acetyl carbons carry all three bonds within 240 degrees
    mol = smiles('CCN(CC)C(=S)S[C@H]1[C@@H]([C@H]([C@@H]([C@H](O1)COC(=O)C)O[C@H]2[C@@H]([C@H]([C@@H]([C@H](O2)'
                 'COC(=O)C)OC(=O)C)OC(=O)C)OC(=O)C)OC(=O)C)OC(=O)C')
    plane = chython_layout(mol)
    assert shape(mol, plane)['ang'] < 90
    d = defects(mol, plane)
    assert (d['ovl61'], d['crs'], d['aob']) == (0, 0, 0)


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
@pytest.mark.parametrize('s', ['c1ccc2c(c1)[nH]c1ccccc12', 'C1c2ccccc2-c2ccccc21', 'c1ccc2c(c1)oc1ccccc12',
                               'O=C1c2ccccc2-c2ccccc21', 'CC1(C)c2ccccc2-c2ccc(N)cc21', 'CC(=O)OCC1c2ccccc2-c2ccccc21'])
def test_a_ring_fused_on_two_bonds_is_regular(s):
    # carbazole, fluorene, dibenzofuran, fluorenone: the pentagon regular, its benzo rings with it
    mol = smiles(s)
    plane = chython_layout(mol)
    assert max(_rmsd(plane, r) for r in mol.sssr) < 1e-6
    d = defects(mol, plane)
    assert (d['ovl61'], d['crs'], d['aob']) == (0, 0, 0)


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
@pytest.mark.parametrize('s', [
    'C[C@]12CCC(=O)N(C1CCC3C2CC[C@]4(C3CC[C@]4(C#CCCCI)O)C)C',                         # PubChem 44266403
    'C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@]2(C#C)O)CCC4=C3C=CC(=C4)O',                     # ethinylestradiol, CID 5991
])
def test_a_pentagon_atom_carrying_two_branches_is_regular_with_them_on_the_thirds(s):
    mol = smiles(s)
    plane = chython_layout(mol)
    assert max(_rmsd(plane, r) for r in mol.sssr) < 1e-6
    d = defects(mol, plane)
    assert (d['ovl61'], d['crs'], d['aob']) == (0, 0, 0)


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
@pytest.mark.parametrize('s', [
    'CCOC(=O)N(CC1=C(C=C(C=C1)OC)OC)C[C@@]2(CC[C@@H]3[C@@]2(CCC=C(CC[C@@H](CC4=CC(=C3C=C4)C(=O)C5=CC=CC=C5C6=CC=CC=C6)'
    'O)C)C)O',                                                                          # PubChem 6728223
    'CCN[C@@H](CC1=CC=C(C=C1)O)CN[C@@H](CC2=CC=C(C=C2)O)CN[C@@H](CC3=CC=CC=C3)CN4CCC[C@H]4CN[C@@H](CN5CCC[C@H]5CN)'
    '[C@@H](C)O',                                                                       # PubChem 44415029
    # PubChem 44373163
    'CC(C)C(C(=O)C(F)(F)F)NC(=O)CN(C1CCCCCC1)C(=O)C(C(C)C)NC(=O)C2=CC=C(C=C2)C(=O)NS(=O)(=O)C3=CC=C(C=C3)Cl',
])
def test_a_splayed_atom_whose_thirds_are_taken_turns_with_its_branches_about_its_stem(s):
    mol = smiles(s)
    plane = chython_layout(mol)
    assert shape(mol, plane)['ang'] < 31
    d = defects(mol, plane)
    assert (d['ovl61'], d['crs'], d['aob']) == (0, 0, 0)


@pytest.mark.xfail(strict=False, reason='a convention the cascade does not yet keep')
@pytest.mark.parametrize('s', [
    'CC1(C)CC[C@]2(CC[C@]3(C)C(=CC[C@@H]4[C@@]5(C)CC[C@H](O)C(C)(C)[C@@H]5CC[C@]43C)[C@@H]2C1)C(=O)O',  # CID 10494
    # CID 10114
    'CC1(C)[C@@H](O)CC[C@]2(C)[C@H]1CC[C@@]1(C)[C@@H]2C(=O)C=C2[C@@H]3C[C@@](C)(C(=O)O)CC[C@]3(C)CC[C@]21C',
])
def test_two_branches_sharing_a_bay_split_it(s):
    # oleanolic and glycyrrhetinic acids: the 8- and 14-methyls face one bay; unsplit, one lies inside a ring
    mol = smiles(s)
    plane = chython_layout(mol)
    assert not shape(mol, plane)['inr']
    d = defects(mol, plane)
    assert (d['ovl61'], d['crs'], d['aob']) == (0, 0, 0)
