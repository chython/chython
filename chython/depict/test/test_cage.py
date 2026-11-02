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
"""A perspective cage's depth: its stereocentres read from their faces, the rear bond of a crossing broken."""
import pytest

from chython.core import smiles
from chython.core.wedge import wedges_for_write
from chython.depict.bonds import bond_paths
from chython.depict.figure import molecule_scene
from chython.depict.layout.cage import _units, cages
from chython.depict.layout.molecule import layout2d

CAGED = {                                       # public compounds, a stereocentre on each kind of cage atom
    'cocaine': 'CN1[C@H]2CC[C@@H]1[C@H]([C@H](C2)OC(=O)c1ccccc1)C(=O)OC',
    'teloidine': 'CN1[C@@H]2C[C@@H](O)C[C@H]1[C@@H](O)[C@H]2O',
    'valtropine': 'CN1[C@@H]2C[C@@H](C[C@H]1[C@@H](C2)OC(=O)CC(C)C)OC(=O)CC(C)C',
    '6-hydroxytropinone': 'CN1[C@@H]2C[C@H](O)[C@H]1CC(=O)C2',
    'borneol': 'C[C@@]12CC[C@@H](C1(C)C)C[C@@H]2O',
}


def _drawn(text):
    mol = smiles(text)
    return mol, layout2d(mol, engine='chython', force=True)


UNREAD = {'teloidine': {11},                    # a hydroxy carbon on the crossing of two cage bonds
          'cocaine': {8}}                       # the benzoate between its two faces: wedged


@pytest.mark.parametrize('name', sorted(CAGED))
def test_every_stereocentre_of_a_cage_reads_from_its_face(name):
    mol, plane = _drawn(CAGED[name])
    (cage,) = cages(mol, plane)
    assert cage.shows == set(_units(mol)) & cage.atoms - UNREAD.get(name, set())


def test_a_substituted_stereocentre_on_a_crossing_takes_a_wedge():
    """Teloidine's C6 hydroxy carbon is drawn on the crossing of two cage bonds: its bond is lost there, so it is
    not read from its face but wedged."""
    from math import hypot
    from chython.depict.layout.cage import CROSS
    from chython.depict.layout.lattice import L
    mol, plane = _drawn(CAGED['teloidine'])
    (cage,) = cages(mol, plane)
    near = {a for a in cage.atoms for _, (x, y) in cage.gaps if hypot(plane[a][0] - x, plane[a][1] - y) < CROSS * L
            and any(y not in cage.atoms for y in mol.neighbors_of(a))}
    assert near and not near & cage.shows
    chosen, _ = wedges_for_write(mol, plane=plane)
    assert near <= {n for n, _, _ in chosen}


def test_a_shown_stereocentre_takes_no_wedge_and_the_rear_bond_breaks():
    """Valtropine's four stereocentres all read from their faces: no wedge is drawn, the scene logs none missing,
    and the rear bond of the crossing is drawn in two pieces."""
    from chython.depict.label import labels
    from chython.depict.style import get_depict_style
    from chython.depict.wedge import wedge_paths
    mol, plane = _drawn(CAGED['valtropine'])
    (cage,) = cages(mol, plane)
    style = get_depict_style()
    boxes = labels(mol, plane, style)
    assert wedge_paths(mol, plane, boxes, style)[0], 'without the depth, the cage is wedged'
    paths, claimed = wedge_paths(mol, plane, boxes, style, shown=cage.shows)
    assert not paths and not claimed
    log = []
    molecule_scene(mol, plane=plane, log=log)
    assert not [r for r in log if r.code == 'depict:unwedged']
    (rear, _), = cage.gaps
    (ax, ay), (bx, by) = plane[rear[0]], plane[rear[1]]

    def on_rear(sub):                           # every vertex of the subpath on the rear bond's line, inside it
        pts = [cmd[1:3] for cmd in sub]
        return all(abs((bx - ax) * (y - ay) - (by - ay) * (x - ax)) < 1e-6 and
                   min(ax, bx) - 1e-6 <= x <= max(ax, bx) + 1e-6 for x, y in pts)
    pieces = [sub for p in bond_paths(mol, plane, boxes, style, gaps={k: [q] for k, q in cage.gaps})
              for sub in p.subpaths if on_rear(sub)]
    assert len(pieces) == 2, 'the rear bond is drawn in two pieces'


def test_a_mirrored_drawing_reads_with_its_depth_turned_round():
    """Reflected, the drawing is the same compound only with the other bond in front: the reading follows."""
    mol, plane = _drawn(CAGED['teloidine'])
    (cage,) = cages(mol, plane)
    mirror = {x: (-p[0], p[1]) for x, p in plane.items()}
    (back,) = cages(mol, mirror)
    assert back.shows == cage.shows
    assert {g for g, _ in back.gaps} != {g for g, _ in cage.gaps}


def test_a_written_cage_keeps_its_wedges():
    """Only the picture reads depth: a file states a cage's configuration with flat wedges, as the format does."""
    mol, plane = _drawn(CAGED['cocaine'])
    chosen, _ = wedges_for_write(mol, plane=plane)
    assert set(_units(mol)) <= {n for n, _, _ in chosen}


def test_a_benzo_fused_cage_is_drawn_as_its_cage_with_the_ring_on_one_edge():
    """Pentazocine's benzomorphan: the bicyclo[3.3.1]nonane and its benzene in one perspective, the benzene a
    regular hexagon seen at an angle: opposite atoms about one centre."""
    from math import hypot
    from chython.depict.layout.lattice import L
    mol, plane = _drawn('CC(C)=CCN1CC[C@]2(C)[C@@H](C)[C@@H]1Cc1ccc(O)cc12')
    (cage,) = cages(mol, plane)
    assert len(cage.atoms) == 9 and cage.shows
    ring = next(r for r in mol.sssr if len(r) == 6 and all(mol.atom(x).hybridization in (2, 4) for x in r))
    mids = [((plane[ring[k]][0] + plane[ring[k + 3]][0]) / 2, (plane[ring[k]][1] + plane[ring[k + 3]][1]) / 2)
            for k in range(3)]
    assert max(hypot(m[0] - mids[0][0], m[1] - mids[0][1]) for m in mids) < .02 * L
    sides = [hypot(plane[a][0] - plane[b][0], plane[a][1] - plane[b][1]) for a, b in zip(ring, ring[1:] + ring[:1])]
    assert min(sides) > .6 * L, 'seen at an angle, not edge on'


def test_a_fused_cages_one_atom_bridge_stands_out_of_its_outline():
    """Pentazocine draws from the pseudopelletierine drawing: its methano bridge on the outline, as camphor's is."""
    from chython.depict.layout.shapes import _inside
    mol, plane = _drawn('CC(C)=CCN1CC[C@]2(C)[C@@H](C)[C@@H]1Cc1ccc(O)cc12')
    (cage,) = cages(mol, plane)
    heads = {x for x in cage.atoms if sum(y in cage.atoms for y in mol.neighbors_of(x)) == 3}
    (bridge,) = [x for x in cage.atoms if x not in heads and sum(y in heads for y in mol.neighbors_of(x)) == 2]
    ring = {x for r in mol.sssr for x in r}
    atoms = sorted(ring)
    assert bridge not in {atoms[k] for k in _inside([plane[x] for x in atoms])}


@pytest.mark.parametrize('name', sorted(CAGED) + ['pentazocine'])
def test_coupled_bridgeheads_are_shown_together(name):
    """A cage states its two coupled bridgeheads at once: one shown means its partner takes no wedge either."""
    text = CAGED.get(name, 'CC(C)=CCN1CC[C@]2(C)[C@@H](C)[C@@H]1Cc1ccc(O)cc12')
    mol, plane = _drawn(text)
    (cage,) = cages(mol, plane)
    units = _units(mol)
    for a in cage.shows:
        assert set(units[a]['coupled']) & cage.atoms <= cage.shows


@pytest.mark.parametrize('text', ['C1CC2CCC1O2', 'C1CC2CCC1C2'])
def test_a_one_atom_bridge_passes_in_front(text):
    """7-Oxanorbornane and norbornane: read with no stereocentre to decide, the bridge is toward the viewer, so the
    bond broken at the crossing is never one of the bridge's."""
    mol, plane = _drawn(text)
    (cage,) = cages(mol, plane)
    heads = {x for x in cage.atoms if sum(y in cage.atoms for y in mol.neighbors_of(x)) == 3}
    (bridge,) = [x for x in cage.atoms if x not in heads and set(mol.neighbors_of(x)) <= heads]
    assert cage.gaps and all(bridge not in g for g, _ in cage.gaps)


def test_a_cage_fused_to_a_saturated_ring_draws_in_perspective():
    """PubChem 5055954's norbornene imide: the succinimide fused on the cage's edge is drawn in its perspective too,
    a regular pentagon seen as the cage is."""
    mol, plane = _drawn('C1CN(CCN1C2=C(C=C(C=C2)N3C(=O)C4C5CC(C4C3=O)C=C5)[N+](=O)[O-])C6=CC=CC=C6Cl')
    (cage,) = cages(mol, plane)
    assert len(cage.atoms) == 7
    imide = next(r for r in mol.sssr if len(r) == 5 and any(mol.atom(x).atomic_symbol == 'N' for x in r))
    mids = [plane[x] for x in imide]
    assert len(set(imide) & cage.atoms) == 2 and all(isinstance(q[0], float) for q in mids)


def test_a_bridgeless_cage_brings_its_substituted_bridge_forward():
    """Bicyclo[2.2.2]oct-5-ene-2,3-dicarboxylic acid, no stereocentre to decide and no one-atom bridge: the bridge
    bearing the acids is read in front, so the broken bond at the crossing is the alkene's side."""
    mol, plane = _drawn('OC(=O)C1C2CCC(C=C2)C1C(=O)O')
    (cage,) = cages(mol, plane)
    substituted = {x for x in cage.atoms if any(y not in cage.atoms for y in mol.neighbors_of(x))}
    assert cage.gaps and all(not set(g) & substituted for g, _ in cage.gaps)


def test_a_projected_rings_substituent_points_out_from_its_centre():
    """Pentazocine's phenol: the benzene seen at an angle with the cage, its hydroxy points from the projected ring's
    centre, parallel to its ring bonds as in 3D, not on the bisector of the foreshortened angle."""
    from math import atan2, degrees
    mol, plane = _drawn('CC(C)=CCN1CC[C@]2(C)[C@@H](C)[C@@H]1Cc1ccc(O)cc12')
    ring = next(r for r in mol.sssr if len(r) == 6 and all(mol.atom(x).hybridization in (2, 4) for x in r))
    cx = sum(plane[x][0] for x in ring) / 6
    cy = sum(plane[x][1] for x in ring) / 6
    (c, o), = [(x, y) for x in ring for y in mol.neighbors_of(x) if mol.atom(y).atomic_symbol == 'O']
    radial = degrees(atan2(plane[c][1] - cy, plane[c][0] - cx))
    bond = degrees(atan2(plane[o][1] - plane[c][1], plane[o][0] - plane[c][0]))
    assert abs((bond - radial + 180) % 360 - 180) < 1


def test_a_substituent_between_its_faces_is_not_read():
    """Cocaine's benzoate leaves C3 nearly as near one face as the other: it is not read from the drawing but
    wedged, and every substituent that is read lies `READ` nearer its own face."""
    from chython.depict.layout.cage import READ, faces
    from math import atan2, pi
    mol, plane = _drawn(CAGED['cocaine'])
    (cage,) = cages(mol, plane)
    units = _units(mol)
    chosen, _ = wedges_for_write(mol, plane=plane)
    assert 8 not in cage.shows and 8 in {n for n, _, _ in chosen}
    for a in cage.shows:
        exo = [y for y in mol.neighbors_of(a) if y not in cage.atoms]
        fs = faces(mol, plane, {x: cage.z[x] for x in cage.atoms}, a)
        if a not in units or len(exo) != 1 or fs is None or len(fs) != 2:
            continue
        t = atan2(plane[exo[0]][1] - plane[a][1], plane[exo[0]][0] - plane[a][0])
        d = sorted(abs((t - atan2(f[1], f[0]) + pi) % (2 * pi) - pi) for f in fs)
        assert d[1] - d[0] >= READ
