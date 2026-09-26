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
"""Peptide cross-links drawn as brackets, head-to-tail lines and teleports."""
from math import hypot

import pytest

from chython.core import smiles
from chython.core.monomers import monomers
from chython.core.test.peptides import CASES, PEPTIDES, molecule
from chython.depict import figure
from chython.depict.bonds import segment_hits_box
from chython.depict.figure import molecule_scene
from chython.depict.label import labels
from chython.depict.layout import molecule as layout
from chython.depict.layout.molecule import clean2d, layout2d
from chython.depict.peptide import BADGE, _free_directions, _on_lattice, link_paths, rounded, routed_links
from chython.depict.scene import Path
from chython.depict.style import get_depict_style


def _kinds(name):
    mol = molecule(name)
    return [k.kind for k in routed_links(mol, layout2d(mol, force=True))]


@pytest.mark.parametrize('name, brackets, teleports', [('oxytocin', 1, 0), ('linaclotide', 2, 1),
                                                       ('lactam', 1, 0)])
def test_counts(name, brackets, teleports):
    kinds = _kinds(name)
    assert (kinds.count('bracket'), kinds.count('teleport')) == (brackets, teleports)


@pytest.mark.parametrize('name', ['gramicidin S', 'cyclosporin A'])
def test_head_to_tail(name):
    assert _kinds(name) == ['head_to_tail']


@pytest.mark.parametrize('name', ['leu-enkephalin', 'semaglutide', 'branched'])
def test_no_links_nothing_routed(name):
    assert _kinds(name) == []


@pytest.mark.parametrize('name', ['oxytocin', 'gramicidin S'])        # the engine takes 18 s on linaclotide
def test_engine_plane_draws_plain(name):
    mol = molecule(name)
    assert routed_links(mol, layout2d(mol, force=True, peptide=False)) == ()


def _points(path):
    return [seg[-2:] for sub in path.subpaths for seg in sub if seg[0] != 'Z']


@pytest.mark.parametrize('name', PEPTIDES + ('lactam',))
def test_symbol_clearance(name):
    """No bond, link line or badge touches a label box it does not belong to."""
    mol = molecule(name)
    plane = layout2d(mol, force=True)
    style = get_depict_style()
    boxes = labels(mol, plane, style)
    ink = {x: lb.box for x, lb in boxes.items() if lb.text is not None}
    links = routed_links(mol, plane)
    routed = {frozenset((k.a, k.b)) for k in links}
    hits = []
    for b in mol.bonds():
        if frozenset((b.n, b.m)) not in routed:
            hits.extend(('bond', b.n, b.m, x) for x, bx in ink.items()
                        if x not in (b.n, b.m) and segment_hits_box(plane[b.n], plane[b.m], bx))
    own = {x for k in links for x in (k.a, k.b)}
    for node in link_paths(mol, plane, links, boxes, style):
        if not isinstance(node, Path):
            continue
        if node.fill is not None:                             # a badge: its centre clear of every symbol
            pts = _points(node)
            cx, cy = sum(p[0] for p in pts) / len(pts), sum(p[1] for p in pts) / len(pts)
            hits.extend(('badge', x) for x, bx in ink.items()
                        if bx.min_x - BADGE < cx < bx.max_x + BADGE and bx.min_y - BADGE < cy < bx.max_y + BADGE)
            continue
        pts = _points(node)
        for p, q in zip(pts, pts[1:]):
            hits.extend(('link', x) for x, bx in ink.items() if x not in own and segment_hits_box(p, q, bx))
    assert hits == []


def test_small_molecule_untouched(monkeypatch):
    mol = smiles('CC(=O)Oc1ccccc1C(=O)O')
    plane = layout2d(mol, force=True)
    drawn = molecule_scene(mol, plane=plane).to_svg()
    monkeypatch.setattr(figure, 'routed_links', lambda mol, plane: ())
    assert molecule_scene(mol, plane=plane).to_svg() == drawn


def test_routed_bond_leaves_the_bond_layer(monkeypatch):
    mol = molecule('oxytocin')
    plane = layout2d(mol, force=True)
    routed = molecule_scene(mol, plane=plane).to_svg()
    monkeypatch.setattr(figure, 'routed_links', lambda mol, plane: ())
    assert molecule_scene(mol, plane=plane).to_svg() != routed


def test_stored_plane_draws_brackets(monkeypatch):
    mol = molecule('oxytocin')
    clean2d(mol)
    monkeypatch.setattr(layout, 'peptide_layout', None)       # a recomputation would raise
    assert [k.kind for k in routed_links(mol, mol.coordinates())] == ['bracket']
    molecule_scene(mol).to_svg()


def test_rounded_corner_is_one_cubic():
    assert [s[0] for s in rounded([(0., 0.), (1., 0.), (1., 1.)])] == ['M', 'L', 'C', 'L']
    assert [s[0] for s in rounded([(0., 0.), (1., 0.), (2., 0.)])] == ['M', 'L', 'L']


def test_rounded_arc_radius():
    """A right-angle corner of radius .25 starts .25 before the vertex and ends .25 after it."""
    segs = rounded([(0., 0.), (1., 0.), (1., 1.)])
    assert segs[1][1:] == pytest.approx((.75, 0.))
    assert segs[2][-2:] == pytest.approx((1., .25))


def _bracket_ends(name):
    """`(mol, plane, link, leaving)`: the bracket path's first direction at each end atom."""
    mol = molecule(name)
    plane = layout2d(mol, force=True)
    style = get_depict_style()
    boxes = labels(mol, plane, style)
    (link,) = routed_links(mol, plane)
    (path,) = [n for n in link_paths(mol, plane, (link,), boxes, style) if isinstance(n, Path)]
    pts = _points(path)
    a, b = (link.a, link.b) if plane[link.a][0] <= plane[link.b][0] else (link.b, link.a)
    leaving = {}
    for x, (p, q) in ((a, pts[:2]), (b, pts[:-3:-1])):
        d = hypot(q[0] - p[0], q[1] - p[1])
        leaving[x] = ((q[0] - p[0]) / d, (q[1] - p[1]) / d)
    return mol, plane, link, leaving


@pytest.mark.parametrize('name', ['triazole staple', 'alkene staple', 'lactam'])
def test_bracket_leg_leaves_on_a_free_direction(name):
    """A ring end leaves on its exterior bisector, a chain end on a 120-degree slot: never a forced vertical."""
    mol, plane, link, leaving = _bracket_ends(name)
    keys = {frozenset((link.a, link.b))}
    for x, (dx, dy) in leaving.items():
        want = _free_directions(mol, plane, x, keys)
        assert max(dx * u + dy * v for u, v in want) == pytest.approx(1, abs=1e-4)


def test_lactam_carbonyl_gives_its_slot_to_the_link():
    """The acyl end's oxygen takes the slot facing the backbone, so its leg drops straight into the bar."""
    mol, plane, link, leaving = _bracket_ends('lactam')
    (c,) = [x for x in (link.a, link.b) if mol.element_of(x) == 6]
    (o,) = [y for y in mol.neighbors_of(c) if mol.element_of(y) == 8]
    assert link.side == -1 and plane[o][1] > plane[c][1]
    assert leaving[c] == pytest.approx((0., -1.), abs=1e-6)


def test_a_vertical_bond_leg_carries_straight_on():
    """The lactam's amine end hangs on a vertical bond: its leg continues that line into the bar."""
    mol, plane, link, leaving = _bracket_ends('lactam')
    (n,) = [x for x in (link.a, link.b) if mol.element_of(x) == 7]
    (y,) = [y for y in mol.neighbors_of(n) if y not in (link.a, link.b)]
    assert plane[n][0] == pytest.approx(plane[y][0], abs=1e-6) and plane[n][1] < plane[y][1]
    assert leaving[n] == pytest.approx((0., -1.), abs=1e-6)


@pytest.mark.parametrize('name', ['oxytocin', 'linaclotide', 'lactam'])
def test_a_leaf_link_end_leans_toward_its_partner(name):
    """A bracket's leaf ends never lean out, so its legs mirror each other rather than run parallel."""
    mol = molecule(name)
    plane = layout2d(mol, force=True)
    for k in routed_links(mol, plane):
        if k.kind != 'bracket':
            continue
        for x, o in ((k.a, k.b), (k.b, k.a)):
            ns = [y for y in mol.neighbors_of(x) if y != o]
            if len(ns) == 1:
                dx, to = plane[x][0] - plane[ns[0]][0], plane[o][0] - plane[ns[0]][0]
                assert dx * to > -1e-6, (name, x)


@pytest.mark.parametrize('other, kinds', [('leu-enkephalin', ['bracket']),
                                          ('linaclotide', ['bracket'] * 3 + ['teleport'])])
def test_each_peptide_component_on_its_own_lattice(other, kinds):
    """Two peptides in one record: each backbone on the lattice, each link routed within its own component."""
    mol = smiles(CASES['oxytocin'] + '.' + CASES[other])
    plane = layout2d(mol, force=True)
    for c in mol.connected_components:
        assert _on_lattice(plane, monomers(mol.substructure(c)).backbone)
    assert sorted(k.kind for k in routed_links(mol, plane)) == kinds
    molecule_scene(mol, plane=plane).to_svg()


def test_a_counter_ion_goes_to_the_engine(monkeypatch):
    """A peptide salt: the peptide on its lattice, the acid laid out by the engine alone."""
    mol = smiles(CASES['oxytocin'] + '.CC(=O)O')
    seen = []
    engine = layout._engine_layout
    monkeypatch.setattr(layout, '_engine_layout', lambda m, e: seen.append(len(m)) or engine(m, e))
    plane = layout2d(mol, force=True)
    assert 4 in seen and max(seen) < 10
    assert [k.kind for k in routed_links(mol, plane)] == ['bracket']
