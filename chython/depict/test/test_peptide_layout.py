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
"""The peptide layout: the lattice solver, the placement invariants, the collision index, the tile cache,
and the `layout2d(peptide=)` switch."""
from math import atan2, degrees, hypot

import pytest

from chython.core import ReactionContainer, smiles
from chython.core.monomers import monomers
from chython.core.test.peptides import CASES, PEPTIDES, explicit, molecule, peptide
from chython.core.wedge import cis_trans_parity
from chython.depict import get_clean2d_engine, get_peptide_layout, set_peptide_layout
from chython.depict.layout.molecule import _engine_layout, _rescale_plane, layout2d
from chython.depict.layout.peptide import L, _backbone_dirs, _cross, _ring_systems, _seg_dist, _solve, \
    peptide_layout


LAYOUTS = PEPTIDES + ('tetrapeptide', 'Ac3c', 'branched', 'lactam', 'triazole staple', 'alkene staple',
                      'side-chain branch', 'Pro-Pro', 'collagen', 'C-terminal Pro-Pro', 'Pro3', 'Pro4')


class _Tile:
    """The engine tile, counting its calls."""
    def __init__(self):
        self.calls = 0

    def __call__(self, sub):
        self.calls += 1
        plane = _engine_layout(sub, get_clean2d_engine())
        _rescale_plane(sub, plane)
        return plane


def _layout(name, **kw):
    mol = molecule(name)
    pep = monomers(mol)
    tile = _Tile()
    return mol, pep, peptide_layout(mol, pep, tile, **kw), tile


def _audit(mol, pep, plane):
    """Invariant violations, cross-links excluded: they are routed at draw time."""
    links = {frozenset(x) for x in pep.crosslinks}
    rings = _ring_systems(mol, pep.members)
    bonds = [(b.n, b.m) for b in mol.bonds() if frozenset((b.n, b.m)) not in links]
    backbone = {frozenset(x) for x in zip(pep.backbone, pep.backbone[1:])}
    atoms = list(mol)

    def d(a, b):
        return hypot(plane[a][0] - plane[b][0], plane[a][1] - plane[b][1])
    out = {'close': 0, 'cross': 0, 'on_bond': 0, 'long': 0, 'angle': 0.}
    for i, a in enumerate(atoms):
        for b in atoms[i + 1:]:
            if b not in mol.neighbors_of(a) and d(a, b) < .5 * L:
                out['close'] += 1
    for i, (a, b) in enumerate(bonds):
        if d(a, b) > 1.15 * L and not (a in rings and b in rings and abs(d(a, b) - 3 ** .5 * L) < .01) \
                and not (frozenset((a, b)) in backbone and abs(d(a, b) - 2 * L) < .01):
            out['long'] += 1                                # nor are a five-ring's chord and a bridge
        for c, e in bonds[i + 1:]:
            if len({a, b, c, e}) == 4 and _cross(plane[a], plane[b], plane[c], plane[e]):
                out['cross'] += 1
        for x in atoms:
            if x != a and x != b and _seg_dist(plane[x], plane[a], plane[b]) < .4 * L:
                out['on_bond'] += 1
    for x in atoms:
        nb = [y for y in mol.neighbors_of(x) if frozenset((x, y)) not in links]
        if x in rings or mol.hybridization_of(x) == 3 or len(nb) not in (2, 3):
            continue
        angs = sorted(degrees(atan2(plane[y][1] - plane[x][1], plane[y][0] - plane[x][0])) % 360 for y in nb)
        gaps = [(angs[(q + 1) % len(angs)] - angs[q]) % 360 for q in range(len(angs))]
        dev = abs(min(gaps) - 120) if len(nb) == 2 else max(abs(g - 120) for g in gaps)
        out['angle'] = max(out['angle'], dev)
    return out


# --- the side solver --------------------------------------------------------------------------------------
@pytest.mark.parametrize('n', [2, 3, 8])
def test_solve_default_is_empty(n):
    assert _solve(n, {}, set()) == set()


def _runs(inv):
    runs, cur = [], []
    for k in sorted(inv):
        if cur and k != cur[-1] + 1:
            runs.append(cur)
            cur = []
        cur.append(k)
    return runs + [cur] if cur else runs


@pytest.mark.parametrize('pins', [{3: True}, {2: True, 6: True}, {4: True, 5: False}, {1: True, 8: True}])
def test_solve_inverts_in_even_runs(pins):
    inv = _solve(10, pins, set())
    assert inv and all(len(r) % 2 == 0 or r[-1] == 8 for r in _runs(inv))
    assert all(1 <= k <= 8 for k in inv)


def test_even_runs_keep_the_axis():
    """A step shifts the chain; the last bond keeps the default direction."""
    for inv in ({2, 3}, {2, 3, 4, 5}, {1, 2, 6, 7}):
        assert _backbone_dirs(10, inv)[-1] == _backbone_dirs(10, set())[-1]


def test_equal_side_pair():
    inv = _solve(10, {}, {(4, 5)})
    assert (4 in inv) != (5 in inv)


def test_same_side_pairs():
    """A proline: N and C-alpha differ, the carbonyls either side follow them."""
    inv = _solve(12, {}, {(4, 5)}, {(3, 4), (5, 6)})
    assert (4 in inv) != (5 in inv) and (3 in inv) == (4 in inv) and (5 in inv) == (6 in inv)


def test_an_odd_run_only_ends_the_chain():
    inv = _solve(10, {8: True}, set())
    assert inv == {8}
    assert _backbone_dirs(10, inv)[:-1] == _backbone_dirs(10, set())[:-1]


# --- placement --------------------------------------------------------------------------------------------
@pytest.mark.parametrize('name', LAYOUTS)
def test_layout_invariants(name):
    mol, pep, plane, _ = _layout(name)
    assert set(plane) == set(mol)
    assert _audit(mol, pep, plane) == {'close': 0, 'cross': 0, 'on_bond': 0, 'long': 0,
                                       'angle': pytest.approx(0, abs=1)}


@pytest.mark.parametrize('name', LAYOUTS)
def test_grid_agrees_with_brute_force(name):
    """Every candidate placement gets one collision count from the grid and from the full scan."""
    seen = []
    _layout(name, _check=lambda grid, brute: seen.append((grid, brute)))
    assert seen and all(g == b for g, b in seen)


def test_tile_cache_one_call_per_ring_system():
    mol, pep, _, tile = _layout('semaglutide')
    rings = _ring_systems(mol, pep.members)
    kinds = {mol.substructure(sorted(s | {y for x in s for y in mol.neighbors_of(x)})).canonical_bytes
             for s in set(rings.values())}
    assert tile.calls == len(kinds) == 4                  # His, Phe, Tyr, Trp


def test_no_side_chain_ring_calls_no_engine():
    assert _layout('cyclosporin A')[3].calls == 0


def test_disulfide_ring_is_no_tile():
    """A ring closed through a cross-link spans residues, so it is routed rather than tiled."""
    mol = molecule('oxytocin')
    rings = _ring_systems(mol, monomers(mol).members)
    assert all(len(s) <= 9 for s in rings.values())


def test_spiro_anchor():
    mol, pep, plane, _ = _layout('Ac3c')
    ring = next(s for s in _ring_systems(mol, pep.members).values() if len(s) == 3)
    a, b, c = (plane[x] for x in ring)
    sides = [hypot(p[0] - q[0], p[1] - q[1]) for p, q in ((a, b), (b, c), (c, a))]
    assert max(sides) - min(sides) < .05 * L


@pytest.mark.parametrize('name', ['Pro-Pro', 'Pro3', 'collagen'])
def test_proline_is_an_open_hexagon(name):
    """Four ring bonds of `L`, one horizontal chord of sqrt(3) L, 120 degrees at N and C-alpha."""
    mol, pep, plane, _ = _layout(name)
    idx = set(pep.backbone)
    for ring in mol.sssr:
        if len(ring) != 5 or len(idx.intersection(ring)) != 2:
            continue
        sides = sorted(hypot(plane[a][0] - plane[b][0], plane[a][1] - plane[b][1])
                       for a, b in zip(ring, ring[1:] + ring[:1]))
        assert sides[:4] == pytest.approx([L] * 4, abs=1e-6) and sides[4] == pytest.approx(3 ** .5 * L)
        a, b = max(zip(ring, ring[1:] + ring[:1]),
                   key=lambda e: hypot(plane[e[0]][0] - plane[e[1]][0], plane[e[0]][1] - plane[e[1]][1]))
        assert plane[a][1] == pytest.approx(plane[b][1], abs=1e-6)
        for x in idx.intersection(ring):
            angs = sorted(degrees(atan2(plane[y][1] - plane[x][1], plane[y][0] - plane[x][0])) % 360
                          for y in mol.neighbors_of(x))
            gaps = sorted((angs[(q + 1) % len(angs)] - angs[q]) % 360 for q in range(len(angs)))
            want = [120., 240.] if len(gaps) == 2 else [120.] * 3
            assert gaps == pytest.approx(want, abs=1e-6)


@pytest.mark.parametrize('name, count', [('Pro-Pro', 1), ('Pro3', 2), ('Pro4', 3), ('collagen', 3)])
def test_consecutive_prolines_bridge(name, count):
    """The amide between two backbone rings is 2 L on its lattice direction; the chain stays level."""
    mol, pep, plane, _ = _layout(name)
    bb = pep.backbone
    long = [(a, b) for a, b in zip(bb, bb[1:])
            if hypot(plane[a][0] - plane[b][0], plane[a][1] - plane[b][1]) > 1.15 * L]
    assert len(long) == count
    assert all(any(mol.element_of(x) == 7 and mol.in_ring_of(x) for x in ab) for ab in long)
    ys = [plane[x][1] for x in bb]
    assert max(ys) - min(ys) < 4 * L + 1e-6


@pytest.mark.parametrize('seq, height', [('GAPPPAG', 2.5), ('AGPPPPGA', 2.5), ('GPPGPPGPPG', 3.5)])
def test_backbone_is_flattest_found(seq, height):
    """A level backbone costs steps: the solver buys a lower one, and a `same` pair gives way when that stays
    clean."""
    mol = smiles(peptide(seq))
    pep = monomers(mol)
    plane = peptide_layout(mol, pep, _Tile())
    ys = [plane[x][1] for x in pep.backbone]
    assert max(ys) - min(ys) == pytest.approx(height * L, abs=1e-6)
    assert _audit(mol, pep, plane) == {'close': 0, 'cross': 0, 'on_bond': 0, 'long': 0,
                                       'angle': pytest.approx(0, abs=1)}


def test_branch_hangs_below_or_above():
    mol, pep, plane, _ = _layout('branched')
    br = pep.branches[0]
    bb_ys = [plane[x][1] for x in pep.backbone]
    far = [plane[x][1] for x in br.path[-3:]]
    assert all(y < min(bb_ys) for y in far) or all(y > max(bb_ys) for y in far)


@pytest.mark.parametrize('bond', ['/C=C/', '/C=C\\'])
def test_cis_trans_side_chain_kept(bond):
    mol = smiles(peptide('AGFL', special={1: f'NC(CC{bond}C)C(=O)'}))
    units = [u for u in mol.stereo_units() if u['kind'] == 1 and u['parity']]
    assert units
    plane = peptide_layout(mol, monomers(mol), _Tile())
    assert all(cis_trans_parity(mol, u, plane=plane) == u['parity'] for u in units)


# --- layout2d and the switch ------------------------------------------------------------------------------
def _is_lattice(mol, plane):
    """Every backbone bond at 30 deg + a multiple of 60, to .05 deg: a stored plane is rounded to 1e-4."""
    bb = monomers(mol).backbone
    for a, b in zip(bb, bb[1:]):
        d = degrees(atan2(plane[b][1] - plane[a][1], plane[b][0] - plane[a][0])) - 30
        if abs(d - 60 * round(d / 60)) > .05:
            return False
    return True


def test_switch():
    mol = molecule('oxytocin')
    assert _is_lattice(mol, layout2d(mol, force=True))
    assert not _is_lattice(mol, layout2d(mol, force=True, peptide=False))
    set_peptide_layout(False)
    try:
        assert get_peptide_layout() is False
        assert not _is_lattice(mol, layout2d(mol, force=True))
        assert _is_lattice(mol, layout2d(mol, force=True, peptide=True))
    finally:
        set_peptide_layout(True)


def test_threshold():
    """Three residues are no peptide: the engine plane, exactly."""
    mol = molecule('tripeptide')
    assert layout2d(mol, force=True) == layout2d(mol, force=True, peptide=False)
    assert _is_lattice(molecule('tetrapeptide'), layout2d(molecule('tetrapeptide'), force=True))


@pytest.mark.parametrize('flag', ['yes', 1, None])
def test_set_rejects_non_bool(flag):
    with pytest.raises(TypeError):
        set_peptide_layout(flag)


def test_salt():
    """A counter-ion is its own component, placed right of the peptide."""
    mol = smiles(CASES['oxytocin'] + '.OC(=O)C(F)(F)F')
    plane = layout2d(mol, force=True)
    assert set(plane) == set(mol)
    pep, acid = sorted(mol.connected_components, key=len, reverse=True)
    assert min(plane[x][0] for x in acid) > max(plane[x][0] for x in pep) + .8
    assert _is_lattice(mol, plane)


def test_facade_forwards():
    import chython
    chython.peptide_layout = False
    try:
        assert get_peptide_layout() is False and chython.peptide_layout is False
    finally:
        chython.peptide_layout = True
    with pytest.raises(TypeError):
        chython.peptide_layout = 'off'


def test_explicit_hydrogens():
    """Withheld like an engine's and placed at one bond, clear of every other atom."""
    mol = explicit(molecule('oxytocin'))
    heavy = [x for x in mol if mol.element_of(x) != 1]
    plane = layout2d(mol, force=True)
    assert set(plane) == set(mol) and _is_lattice(mol, plane)
    for x in mol:
        if mol.element_of(x) == 1:
            p, (y,) = plane[x], mol.neighbors_of(x)
            assert hypot(p[0] - plane[y][0], p[1] - plane[y][1]) == pytest.approx(L)
            assert min(hypot(p[0] - plane[z][0], p[1] - plane[z][1]) for z in heavy if z != y) > .3


def test_a_ring_in_the_backbone_path_goes_to_the_engine():
    """4-aminobenzoyl's ring hangs on two non-adjacent vertices: no slot takes it, the engine does."""
    mol = molecule('4-aminobenzoyl')
    assert peptide_layout(mol, monomers(mol), _Tile()) is None
    plane = layout2d(mol, force=True)
    assert set(plane) == set(mol)
    assert max(hypot(plane[b.n][0] - plane[b.m][0], plane[b.n][1] - plane[b.m][1]) for b in mol.bonds()) < 1.15 * L


def test_per_call_switch_rejects_non_bool():
    with pytest.raises(TypeError):
        layout2d(molecule('oxytocin'), force=True, peptide='no')


def test_container_method_follows_switch():
    mol = molecule('leu-enkephalin')
    mol.clean2d()
    assert _is_lattice(mol, mol.coordinates())
    assert not _is_lattice(mol, mol.layout2d(force=True, peptide=False))
    mol.clean2d(force=True, peptide=False)
    assert not _is_lattice(mol, mol.coordinates())


def test_reaction_forwards_switch():
    """A reaction's members follow `peptide=` like a molecule does."""
    rxn = ReactionContainer([molecule('leu-enkephalin')], [molecule('oxytocin')])
    planes, _, _ = rxn.layout2d(force=True)
    assert all(_is_lattice(m, p) for m, p in zip(rxn.molecules(), planes))
    planes, _, _ = rxn.layout2d(force=True, peptide=False)
    assert not any(_is_lattice(m, p) for m, p in zip(rxn.molecules(), planes))
    rxn.clean2d(force=True, peptide=False)
    assert not any(_is_lattice(m, m.coordinates()) for m in rxn.molecules())


@pytest.mark.parametrize('engine', ['smilesdrawer', 'rdkit', 'cdk', 'obabel', 'indigo'])
def test_engine_independence(engine):
    """Only ring tiles come from the engine, so the invariants hold whichever lays them out."""
    try:
        _engine_layout(smiles('c1ccccc1C'), engine)
    except Exception as e:                                  # the toolkit is absent
        pytest.skip(f'{engine}: {e}')
    for name in ('semaglutide', 'linaclotide'):
        mol = molecule(name)
        pep = monomers(mol)
        plane = layout2d(mol, engine=engine, force=True)
        assert _audit(mol, pep, plane) == {'close': 0, 'cross': 0, 'on_bond': 0, 'long': 0,
                                           'angle': pytest.approx(0, abs=1)}, name
