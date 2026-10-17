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
"""The peptide layout: the backbone's moves, the placement invariants, the collision index, the tile cache,
and the `layout2d(peptide=)` switch."""
from math import atan2, degrees, hypot

import pytest

from chython.core import ReactionContainer, smiles
from chython.core.monomers import monomers
from chython.core.test.peptides import CASES, PEPTIDES, explicit, molecule, peptide
from chython.core.wedge import cis_trans_parity
from chython.depict import get_clean2d_engine, get_peptide_layout, set_peptide_layout
from chython.depict.layout.molecule import _engine_layout, _chain_plane, _rescale_plane, layout2d
from chython.depict.layout.peptide import (L, ROUTE, _backbone_dirs, _cross, _ring_systems, _seg_dist, peptide_layout,
                                          toggles, valid)
from chython.depict.layout.quality import perm
from chython.depict.peptide import routed_links


@pytest.fixture(autouse=True)
def _chains():
    """Chain layout is off by default; the layout under test here is the chain one."""
    was = get_peptide_layout()
    set_peptide_layout(True)
    yield
    set_peptide_layout(was)


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


# --- the backbone's moves ---------------------------------------------------------------------------------
class _Ctx:
    """A stand-in context holding only a backbone of `n` atoms."""
    def __init__(self, n):
        self.peptide = type('P', (), {'backbone': list(range(n))})()


def _runs(inv):
    runs, cur = [], []
    for k in sorted(inv):
        if cur and k != cur[-1] + 1:
            runs.append(cur)
            cur = []
        cur.append(k)
    return runs + [cur] if cur else runs


@pytest.mark.parametrize('inv, ok', [((), True), ((3, 4), True), ((1, 2, 6, 7), True), ((8,), True),
                                     ((3,), False), ((0, 1), False), ((8, 9), False), ((2, 3, 4), False)])
def test_inverted_vertices_are_even_runs_inside_the_chain(inv, ok):
    assert valid(_Ctx(10), frozenset(inv)) is ok


def test_every_toggle_keeps_runs_even():
    ctx = _Ctx(10)
    seen, todo = {frozenset()}, [frozenset()]
    while todo:
        for new in toggles(ctx, todo.pop()):
            assert all(len(r) % 2 == 0 or r[-1] == 8 for r in _runs(new))
            if new not in seen:
                seen.add(new)
                todo.append(new)
    assert len(seen) > 20                               # the moves reach every even-run set


def test_even_runs_keep_the_axis():
    """A step shifts the chain; the last bond keeps the default direction."""
    for inv in ({2, 3}, {2, 3, 4, 5}, {1, 2, 6, 7}):
        assert _backbone_dirs(10, inv)[-1] == _backbone_dirs(10, set())[-1]


def test_an_odd_run_only_ends_the_chain():
    assert _backbone_dirs(10, {8})[:-1] == _backbone_dirs(10, set())[:-1]


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
def test_proline_is_a_regular_pentagon(name):
    """Five ring bonds of `L`, every ring angle 108 degrees."""
    mol, pep, plane, _ = _layout(name)
    idx = set(pep.backbone)
    for ring in mol.sssr:
        if len(ring) != 5 or len(idx.intersection(ring)) != 2:
            continue
        for a, x, b in zip(ring[-1:] + ring[:-1], ring, ring[1:] + ring[:1]):
            assert hypot(plane[a][0] - plane[x][0], plane[a][1] - plane[x][1]) == pytest.approx(L, abs=1e-6)
            t = degrees(atan2(plane[a][1] - plane[x][1], plane[a][0] - plane[x][0]) -
                        atan2(plane[b][1] - plane[x][1], plane[b][0] - plane[x][0])) % 360
            assert min(t, 360 - t) == pytest.approx(108, abs=1e-6)


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
    """Off: the engine plane; on for one call or for the process."""
    mol = molecule('oxytocin')
    set_peptide_layout(False)
    assert not _is_lattice(mol, layout2d(mol, force=True))
    assert _is_lattice(mol, layout2d(mol, force=True, peptide=True))
    set_peptide_layout(True)
    try:
        assert _is_lattice(mol, layout2d(mol, force=True))
        assert not _is_lattice(mol, layout2d(mol, force=True, peptide=False))
    finally:
        set_peptide_layout(False)


def test_off_by_default():
    import subprocess
    import sys
    out = subprocess.run([sys.executable, '-c', 'import chython.depict as d; print(d.get_peptide_layout())'],
                         capture_output=True, text=True, check=True).stdout
    assert out.strip() == 'False'


def test_threshold():
    """Three residues are no peptide: the engine plane, exactly."""
    mol = molecule('tripeptide')
    assert layout2d(mol, force=True, peptide=True) == layout2d(mol, force=True, peptide=False)
    assert _is_lattice(molecule('tetrapeptide'), layout2d(molecule('tetrapeptide'), force=True, peptide=True))


@pytest.mark.parametrize('flag', ['yes', 1, None])
def test_set_rejects_non_bool(flag):
    with pytest.raises(TypeError):
        set_peptide_layout(flag)


def test_salt():
    """A counter-ion is its own component, placed right of the peptide."""
    mol = smiles(CASES['oxytocin'] + '.OC(=O)C(F)(F)F')
    plane = layout2d(mol, force=True, peptide=True)
    assert set(plane) == set(mol)
    pep, acid = sorted(mol.connected_components, key=len, reverse=True)
    assert min(plane[x][0] for x in acid) > max(plane[x][0] for x in pep) + .8
    assert _is_lattice(mol, plane)


def test_facade_forwards():
    import chython
    chython.peptide_layout = False
    try:
        assert get_peptide_layout() is False and chython.chain_layout is False
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


# constitutions of Ser-Leu-Ser-Asn (PubChem CID 18741409) and gramicidin S (CID 71777)
BLIND = {'linear': 'CC(C)CC(C(=O)NC(CO)C(=O)NC(CC(=O)N)C(=O)O)NC(=O)C(CO)N',
         'cyclic': 'N%99C(C(C)C)C(=O)NC(CCCN)C(=O)NC(CC(C)C)C(=O)NC(Cc1ccccc1)C(=O)N%10CCCC%10C(=O)NC(C(C)C)C(=O)'
                   'NC(CCCN)C(=O)NC(CC(C)C)C(=O)NC(Cc1ccccc1)C(=O)N%11CCCC%11C(=O)%99'}


@pytest.mark.parametrize('kind', BLIND)
def test_layout_is_blind_to_atom_order(kind):
    """Random-order re-reads of one peptide lay out as one drawing: every choice goes by canonical rank."""
    mol = smiles(BLIND[kind])
    assert monomers(mol) is not None
    assert perm(mol, lambda m: layout2d(m, engine='chython', force=True, peptide=True)) == 0


# PubChem 14969
VANCOMYCIN = ('C[C@H]1[C@H]([C@@](C[C@@H](O1)O[C@@H]2[C@H]([C@@H]([C@H](O[C@H]2OC3=C4C=C5C=C3OC6=C(C=C(C=C6)[C@H]([C@H]'
              '(C(=O)N[C@H](C(=O)N[C@H]5C(=O)N[C@@H]7C8=CC(=C(C=C8)O)C9=C(C=C(C=C9O)O)[C@H](NC(=O)[C@H]([C@@H](C1=CC(=C'
              '(O4)C=C1)Cl)O)NC7=O)C(=O)O)CC(=O)N)NC(=O)[C@@H](CC(C)C)NC)O)Cl)CO)O)O)(C)N)O')


def test_a_ring_holding_two_link_ends_faces_each_partner():
    """Vancomycin's central aryl ether is reflected toward its two partner rings: one ether closes at a bond."""
    mol = smiles(VANCOMYCIN)
    plane = _chain_plane(mol, 'chython')
    assert plane is not None
    assert len(routed_links(mol, plane)) < len(monomers(mol).crosslinks)


def test_a_closed_biaryl_stretches_the_bulkier_group():
    """Vancomycin's biaryl closes as a bond with the 4-hydroxyphenyl ring (7 atoms) upright on its slot at one bond
    length; the 3,5-dihydroxyphenyl ring (8 atoms) takes the stretch.  Atoms by SMILES order."""
    mol = smiles(VANCOMYCIN)
    plane = _chain_plane(mol, 'chython')
    (ax, ay), (rx, ry) = plane[42], plane[43]
    assert hypot(rx - ax, ry - ay) == pytest.approx(L) and rx == pytest.approx(ax)
    assert hypot(plane[45][0] - plane[50][0], plane[45][1] - plane[50][1]) <= ROUTE
    assert hypot(plane[51][0] - plane[58][0], plane[51][1] - plane[58][1]) > 1.5 * L


def test_a_closed_alkene_staple_keeps_its_cis_trans():
    from chython.depict.layout.quality import shape
    mol, pep, plane, _ = _layout('alkene staple')
    assert shape(mol, plane)['ct'] == 0


def test_a_link_end_whose_partner_is_not_drawn_yet_keeps_its_leaves():
    from chython.depict.layout.molecule import layout2d
    mol = smiles('CC1=C(C(=C(C=C1)C(=O)NC(C(C)OC(=O)C(C(C)C)N(C)C(=O)OC(C)(C)C)C(=O)N[C@H](C(C)C)C(=O)N2CCCC2C(=O)'
                 'N(CCC(=O)OC)CC(=O)OC(C)(C)C)[N+](=O)[O-])OCC3=CC=CC=C3')       # pubchem 14845251
    plane = layout2d(mol, engine='chython', force=True, peptide=True)
    assert set(plane) == set(mol)


def test_a_cyclodextrin_lays_out_as_a_chain_of_regular_beads():
    from chython.core.chains import perceive
    from chython.core.test.test_chains import ALPHA_CD
    from chython.depict.peptide import _on_lattice
    mol = smiles(ALPHA_CD)
    chain = perceive(mol)
    plane = layout2d(mol, engine='chython', force=True)
    assert _on_lattice(plane, chain.backbone, chain.beads)
    for r in mol.sssr:
        if len(r) == 6:
            for a, b in zip(r, r[1:] + r[:1]):
                assert hypot(plane[a][0] - plane[b][0], plane[a][1] - plane[b][1]) == pytest.approx(L, rel=.05)


def test_a_chain_the_lattice_refuses_falls_back(monkeypatch):
    from chython.core.test.test_chains import ALPHA_CD
    from chython.depict.layout import molecule as entry
    monkeypatch.setattr(entry, 'peptide_layout', lambda *a, **k: None)
    mol = smiles(ALPHA_CD)
    assert set(layout2d(mol, engine='chython', force=True)) == set(mol)


def test_a_calixarene_zigzags_through_three_atoms_of_each_ring():
    from chython.core.chains import perceive
    from chython.core.test.test_chains import CALIX4
    from chython.depict.peptide import _on_lattice
    mol = smiles(CALIX4)
    chain = perceive(mol)
    plane = layout2d(mol, engine='chython', force=True)
    assert _on_lattice(plane, chain.backbone, chain.beads)
    for r in chain.residues:
        if r.kind == 'arc':
            ring = next(q for q in mol.sssr if set(q) == set(r.atoms))
            for a, b in zip(ring, ring[1:] + ring[:1]):
                assert hypot(plane[a][0] - plane[b][0], plane[a][1] - plane[b][1]) == pytest.approx(L, rel=.02)


@pytest.mark.parametrize('name', ['alpha-cd', 'calix4'])
def test_a_chain_of_beads_runs_level(name):
    from chython.core.chains import perceive
    from chython.core.test.test_chains import ALPHA_CD, CALIX4
    mol = smiles({'alpha-cd': ALPHA_CD, 'calix4': CALIX4}[name])
    bb = perceive(mol).backbone
    plane = layout2d(mol, engine='chython', force=True)
    assert abs(_angle(plane, bb[0], bb[-1])) < 10


@pytest.mark.parametrize('name', ['alpha-cd'])
def test_the_backbone_enters_and_leaves_a_bead_along_its_bisectors(name):
    from chython.core.chains import perceive
    from chython.core.test.test_chains import ALPHA_CD, CALIX4
    mol = smiles({'alpha-cd': ALPHA_CD, 'calix4': CALIX4}[name])
    chain = perceive(mol)
    plane = layout2d(mol, engine='chython', force=True)
    links = {frozenset(k) for k in chain.crosslinks}
    for k in chain.beads:
        system = next(r.atoms for r in chain.residues if r.kind == 'bead' and r.entry == chain.backbone[k])
        for x in chain.backbone[k:k + 2]:
            ring = [y for y in mol.neighbors_of(x) if y in system]
            outs = [y for y in mol.neighbors_of(x) if y not in system and frozenset((x, y)) not in links]
            if not outs:
                continue                                # its outgoing bond is the closure, routed apart
            a, b = (abs((_angle(plane, x, outs[0]) - _angle(plane, x, y) + 180) % 360 - 180) for y in ring)
            assert a == pytest.approx(b, abs=2), (k, x)


def _angle(plane, x, y):
    from math import atan2, degrees
    return degrees(atan2(plane[y][1] - plane[x][1], plane[y][0] - plane[x][0]))
