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
from chython.depict.layout import molecule as layout, peptide as peptide_layout
from chython.depict.layout.lattice import _seg_dist
from chython.depict.layout.molecule import clean2d, layout2d
from chython.depict.peptide import (BAR, BADGE, CLOSE, _cross, _free_directions, _on_lattice, link_paths, rounded,
                                    routed_links)
from chython.depict.scene import Path
from chython.depict.style import get_depict_style
from chython.depict.test.test_peptide_layout import VANCOMYCIN


def _kinds(name):
    mol = molecule(name)
    return [k.kind for k in routed_links(mol, layout2d(mol, force=True))]


@pytest.fixture(autouse=True)
def _chains():
    """Chain layout is off by default; every drawing here is of one."""
    from chython.depict import get_chain_layout, set_chain_layout
    was = get_chain_layout()
    set_chain_layout(True)
    yield
    set_chain_layout(was)


@pytest.mark.parametrize('name, brackets, teleports', [('oxytocin', 1, 0), ('linaclotide', 2, 1),
                                                       ('lactam', 1, 0)])
def test_counts(name, brackets, teleports):
    kinds = _kinds(name)
    assert (kinds.count('bracket'), kinds.count('teleport')) == (brackets, teleports)


@pytest.mark.parametrize('name', ['gramicidin S', 'cyclosporin A'])
def test_head_to_tail(name):
    assert _kinds(name) == ['head_to_tail']


def test_a_tail_drops_clear_of_every_atom():
    """Gramicidin S's tail leaves its carbonyl at 120 degrees, lengthened until its drop passes the leucine
    methyl beside it."""
    mol = molecule('gramicidin S')
    plane = layout2d(mol, force=True)
    style = get_depict_style()
    (link,) = routed_links(mol, plane)
    (path,) = [n for n in link_paths(mol, plane, (link,), labels(mol, plane, style), style) if isinstance(n, Path)]
    pts = _points(path)
    for x, (cx, cy) in plane.items():
        if x not in (link.a, link.b):
            assert min(_seg_dist((cx, cy), p, q) for p, q in zip(pts, pts[1:])) > CLOSE, x


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


@pytest.mark.parametrize('name', ['triazole staple', 'alkene staple', 'lactam'])
def test_a_bracket_bar_stands_bar_past_what_it_spans(name):
    """The bar runs `BAR` past the outermost atom between its legs, the link ends among them."""
    mol = molecule(name)
    plane = layout2d(mol, force=True)
    style = get_depict_style()
    (link,) = routed_links(mol, plane)
    (path,) = [n for n in link_paths(mol, plane, (link,), labels(mol, plane, style), style) if isinstance(n, Path)]
    pts = _points(path)
    s = link.side
    x0, x1 = min(x for x, _ in pts), max(x for x, _ in pts)
    outer = max([s * y for x, y in plane.values() if x0 - .1 <= x <= x1 + .1] +
                [s * plane[link.a][1], s * plane[link.b][1]])
    assert max(s * y for _, y in pts) - outer == pytest.approx(BAR, abs=.01)


@pytest.mark.parametrize('name', ['oxytocin', 'linaclotide', 'lactam'])
def test_a_leaf_link_end_leans_toward_its_partner(name):
    """A bracket's leaf ends never lean out, so its legs mirror each other rather than run parallel; a partner on
    the backbone leaves its end as the bond drew it."""
    mol = molecule(name)
    plane = layout2d(mol, force=True)
    backbone = set(monomers(mol).backbone)
    for k in routed_links(mol, plane):
        if k.kind != 'bracket':
            continue
        for x, o in ((k.a, k.b), (k.b, k.a)):
            ns = [y for y in mol.neighbors_of(x) if y != o]
            if len(ns) == 1 and o not in backbone:
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
    monkeypatch.setattr(layout, '_engine_layout', lambda m, e, c=None: seen.append(len(m)) or engine(m, e, c))
    plane = layout2d(mol, force=True)
    assert 4 in seen and max(seen) < 10
    assert [k.kind for k in routed_links(mol, plane)] == ['bracket']


# PubChem 16185392, a protected glycopeptide aglycone: four aryl link ends on rings hosted across the backbone
AGLYCONE = ('CC(C)(C)C(=O)OCOC(=O)C1C2=C(C(=CC(=C2)O)O)C3=C(C=CC(=C3)C4C(=O)NC(C(C5=CC(=C(C=C5)OC6=CC7=CC(=C6O)OC8=C('
            'C=C(CC9C(=O)NC(C2=CC(=CC(=C2)OC2=C(C=CC(=C2)C(C(=O)N9)NC(=O)OCC2=CC=CC=C2)O)O)C(=O)NC7C(=O)N4)C=C8)Cl)Cl'
            ')O)C(=O)N1)O')


@pytest.mark.parametrize('smi', [VANCOMYCIN, AGLYCONE], ids=['vancomycin', 'aglycone'])
def test_a_bracket_leg_crosses_no_bond(smi, monkeypatch):
    """A leg blocked on its way to the bar steps sideways past the obstruction, through no sugar and no own ring;
    a link no leg can clear teleports.  Links left open, as closing them leaves no leg here."""
    monkeypatch.setattr(peptide_layout, '_close_links', lambda ctx, best: best[0])
    mol = smiles(smi)
    plane = layout2d(mol, engine='chython', force=True)
    style = get_depict_style()
    links = [k for k in routed_links(mol, plane) if k.kind == 'bracket']
    brackets = [n for n in link_paths(mol, plane, links, labels(mol, plane, style), style)
                if isinstance(n, Path) and n.fill is None and n.dashes is None]
    assert brackets
    keys = {frozenset((k.a, k.b)) for k in links}
    crossed = []
    for path in brackets:
        pts = _points(path)
        for p, q in zip(pts, pts[1:]):
            crossed.extend((b.n, b.m) for b in mol.bonds()
                           if frozenset((b.n, b.m)) not in keys and _cross(p, q, plane[b.n], plane[b.m]))
    assert crossed == []


@pytest.mark.parametrize('smi', [VANCOMYCIN, AGLYCONE], ids=['vancomycin', 'aglycone'])
def test_links_between_side_chain_rings_close_into_bonds(smi, monkeypatch):
    """Every biaryl and diaryl ether link is drawn as a bond: its two side-chain groups turned and stretched
    toward each other, no atom brought within `CLOSE` of another."""
    mol = smiles(smi)
    atoms = list(mol)

    def near(plane):
        return {(a, b) for i, a in enumerate(atoms) for b in atoms[i + 1:] if b not in mol.neighbors_of(a) and
                hypot(plane[a][0] - plane[b][0], plane[a][1] - plane[b][1]) < peptide_layout.CLOSE}
    plane = layout2d(mol, engine='chython', force=True)
    assert routed_links(mol, plane) == ()
    monkeypatch.setattr(peptide_layout, '_close_links', lambda ctx, best: best[0])
    assert near(plane) <= near(layout2d(mol, engine='chython', force=True))


# PubChem 11309298, a side-chain lactam linking a lysine amine to a glutamyl carbonyl
SIDE_LACTAM = ('C[C@H]([C@H]1C(=O)N[C@H](C(=O)N[C@@H](CCCCNC(=O)CC[C@@H](C(=O)N[C@@H](C(=O)N[C@H](C(=O)N1)'
               'CC2=CC=C(C=C2)CNC(C)C)CC3=CNC4=CC=CC=C43)NC(=O)[C@H](CC5=CC=CC=C5)N)C(=O)N[C@@H](CO)C(=O)O)'
               'CC6=CC=CC=C6)O')
# PubChem 13056982, a keto-chain closing on the N-terminal amine
KETO_CHAIN = ('CCC(C)[C@H]1C(=O)N[C@H](C(=O)N[C@H](C(=O)NC(CCCCCC(=O)N[C@H](C(=O)N1)CC2=CC=C(C=C2)O)C(=O)NN)'
              'CC(=O)N)CCC(=O)N')


def test_an_end_level_with_its_backbone_takes_its_partners_side():
    """A link whose ends both sit level with their backbone atoms brackets rather than teleports."""
    mol = smiles(KETO_CHAIN)
    plane = layout2d(mol, engine='chython', force=True)
    assert [k.kind for k in routed_links(mol, plane)] == ['bracket']


def test_a_cross_link_carbon_end_bearing_hydrogens_is_written():
    mol = molecule('alkene staple')
    plane = layout2d(mol)
    svg = molecule_scene(mol, plane=plane).to_svg()
    (a, b), = monomers(mol).crosslinks
    end = a if mol.atom(a).implicit_h else b
    assert mol.atom(end).atomic_symbol == 'C' and svg.count('>C</tspan>') == 1


def test_a_leg_never_leaves_into_the_angle_of_its_two_bonds():
    from math import atan2, pi
    mol = smiles(SIDE_LACTAM)
    plane = layout2d(mol, engine='chython', force=True)
    for k in routed_links(mol, plane):
        links = {frozenset((k.a, k.b))}
        for x in (k.a, k.b):
            ns = [y for y in mol.neighbors_of(x) if frozenset((x, y)) not in links]
            if len(ns) != 2:
                continue
            a, b = (atan2(plane[y][1] - plane[x][1], plane[y][0] - plane[x][0]) for y in ns)
            if (b - a) % (2 * pi) > pi:
                a, b = b, a                             # the inner arc runs from a to b anticlockwise
            for dx, dy in _free_directions(mol, plane, x, links):
                assert (atan2(dy, dx) - a) % (2 * pi) > (b - a) % (2 * pi) - 1e-6


@pytest.mark.parametrize('name', ['lactam', 'alkene staple', 'triazole staple', 'oxytocin', 'linaclotide'])
def test_every_routed_link_starts_on_the_backbone(name):
    mol = molecule(name)
    plane = layout2d(mol, force=True)
    backbone = set(monomers(mol).backbone)
    for k in routed_links(mol, plane):
        if k.kind != 'head_to_tail':
            assert k.a in backbone or k.b in backbone, (name, k)


def test_on_lattice_takes_a_bead_step_at_any_length():
    from math import cos, radians, sin
    from chython.depict.layout.lattice import L
    plane = {1: (0., 0.), 2: (L * cos(radians(30)), L * sin(radians(30)))}
    plane[3] = (plane[2][0] + 2.5 * L * cos(radians(-30)), plane[2][1] + 2.5 * L * sin(radians(-30)))
    assert _on_lattice(plane, (1, 2, 3), (1,)) and not _on_lattice(plane, (1, 2, 3))


def test_chain_layout_is_the_peptide_switch():
    import chython
    from chython.depict import get_chain_layout, get_peptide_layout, set_chain_layout, set_peptide_layout
    assert get_chain_layout is get_peptide_layout and set_chain_layout is set_peptide_layout
    chython.chain_layout = False
    try:
        assert chython.peptide_layout is False
    finally:
        chython.chain_layout = True


@pytest.mark.parametrize('name', ['alpha-cd', 'calix4'])
def test_a_cycle_closure_leaves_along_free_slots_and_crosses_nothing(name):
    from math import atan2, pi
    from chython.core.test.test_chains import ALPHA_CD, CALIX4
    from chython.depict.peptide import _hits
    mol = smiles({'alpha-cd': ALPHA_CD, 'calix4': CALIX4}[name])
    plane = layout2d(mol, engine='chython', force=True)
    style = get_depict_style()
    boxes = labels(mol, plane, style)
    (link,) = routed_links(mol, plane)
    assert link.kind == 'head_to_tail'
    keys = {frozenset((link.a, link.b))}
    (path,) = [n for n in link_paths(mol, plane, (link,), boxes, style) if isinstance(n, Path)]
    pts = _points(path)
    ink = {x: lb.box for x, lb in boxes.items() if lb.text is not None}
    ends = {link.a, link.b}
    for p, q in zip(pts, pts[1:]):
        assert _hits(mol, plane, None, p, q, ink, keys) <= sum(_hits(mol, plane, None, p, q, ink, keys) >
                                                                _hits(mol, plane, x, p, q, ink, keys) for x in ends)
    for x, (p, q) in ((min(ends, key=lambda x: (plane[x][0] - pts[0][0]) ** 2 + (plane[x][1] - pts[0][1]) ** 2),
                       (pts[0], pts[1])),
                      (min(ends, key=lambda x: (plane[x][0] - pts[-1][0]) ** 2 + (plane[x][1] - pts[-1][1]) ** 2),
                       (pts[-1], pts[-2]))):
        u = atan2(q[1] - p[1], q[0] - p[0])
        assert min(abs((u - atan2(dy, dx) + pi) % (2 * pi) - pi) for dx, dy in _free_directions(mol, plane, x, keys)) \
            < .1
