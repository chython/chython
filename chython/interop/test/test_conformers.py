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
"""
Conformer generation, an external-tool call.

Two properties carry the weight: it stores models on the molecule and returns only a count, and
chirality survives into the geometry.
"""
from pytest import mark, raises

from .conftest import requires_cdpkit, requires_rdkit


def _signed_volume(mol, coords, center):
    """
    The signed volume of the tetrahedron on `center`'s first four neighbours.

    Geometric rather than re-perceived CIP: a perception step can hide a mirrored embedding by agreeing
    with itself.  Adjacency order is stable between two alike-numbered isomorphic inputs.
    """
    env = mol.neighbors_of(center)[:4]
    assert len(env) == 4, 'need four neighbours to measure handedness'
    p0, p1, p2, p3 = (coords[x] for x in env)
    u = [p1[i] - p0[i] for i in range(3)]
    v = [p2[i] - p0[i] for i in range(3)]
    w = [p3[i] - p0[i] for i in range(3)]
    cross = (u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0])
    return sum(cross[i] * w[i] for i in range(3))


@requires_rdkit
def test_stores_a_model_per_conformer_over_the_input_atoms():
    """
    Each generated conformer is a model of the caller's molecule, over the caller's atoms.

    Over the heavy atoms and by the caller's numbers, not a toolkit index: the engine completes the
    hydrogens itself and they are not stored.
    """
    from chython import smiles
    from chython.interop.conformers import generate_conformers

    mol = smiles('CCO')
    stored = generate_conformers(mol, limit=3, engine='rdkit')

    assert stored, 'no conformers generated for ethanol'
    assert len(mol.conformers) == stored
    for c in mol.conformers:
        assert c.ext_index is None, 'nothing generated came from a file'
        coords = c.coordinates
        assert len(coords) == len(mol.atom_numbers)
        assert all(len(xyz) == 3 and all(isinstance(v, float) for v in xyz) for xyz in coords)


@requires_rdkit
def test_generation_changes_no_chemistry():
    """
    Generating conformers writes geometry and nothing else.

    Measured as the canonical form plus the implicit counts: geometry is out of the identity, so a
    molecule that gained two models is still the same compound.
    """
    from chython import smiles
    from chython.interop.conformers import generate_conformers

    mol = smiles('CCO')
    before = bytes(mol.canonical_bytes)
    before_h = {n: mol.implicit_h_of(n) for n in mol.atom_numbers}

    generate_conformers(mol, limit=2, engine='rdkit')

    assert bytes(mol.canonical_bytes) == before, 'the generator changed the caller\'s molecule'
    assert {n: mol.implicit_h_of(n) for n in mol.atom_numbers} == before_h, 'implicit hydrogens changed'
    assert len(mol.conformers) == 2


@requires_rdkit
def test_a_second_run_replaces_rather_than_appends():
    """A generation run is the whole answer, and two runs of one limit would leave models no field on a
    conformer could tell apart."""
    from chython import smiles
    from chython.interop.conformers import generate_conformers

    mol = smiles('CCO')
    generate_conformers(mol, limit=3, engine='rdkit')
    assert generate_conformers(mol, limit=2, engine='rdkit') == 2
    assert len(mol.conformers) == 2


@requires_rdkit
def test_a_layout_is_left_alone():
    """A generated geometry is not a layout; `set_xy` and `set_xyz` are independent by design."""
    from chython import smiles
    from chython.interop.conformers import generate_conformers

    mol = smiles('CCO')
    mol.clean2d()
    before = [mol.xy_of(n) for n in mol.atom_numbers]
    generate_conformers(mol, limit=1, engine='rdkit')
    assert mol.has_3d
    assert [mol.xy_of(n) for n in mol.atom_numbers] == before


def test_an_engine_that_produced_nothing_leaves_the_count_where_it_was(monkeypatch):
    """Nothing generated is nothing dropped: the models the molecule came with are still there."""
    from chython import smiles
    from chython.interop import conformers

    mol = smiles('CCO')
    numbers = mol.atom_numbers          # a clean read, and the session below refuses one
    with mol.edit():
        model = mol.add_conformer()
        for n in numbers:
            mol.set_xyz(n, 1., 1., 1., model=model)
    monkeypatch.setattr(conformers, '_rdkit_conformers', lambda *a, **k: [])
    assert conformers.generate_conformers(mol, engine='rdkit') == 0
    assert len(mol.conformers) == 1


def test_the_coordinate_map_alias_is_gone():
    """With no map returned, nothing is typed by it, and the name belongs to the core's view."""
    from chython.interop import conformers

    assert conformers.__all__ == ['generate_conformers']
    assert not hasattr(conformers, 'Conformer')


def test_unknown_engine_is_refused_by_name():
    """An engine nobody implements is a refusal naming the value, not a silent zero conformers."""
    from chython import smiles
    from chython.interop.conformers import generate_conformers

    with raises(ValueError, match='no-such-engine'):
        generate_conformers(smiles('CCO'), engine='no-such-engine')


def test_the_engine_default_is_read_from_interop_config(monkeypatch):
    """With no `engine=`, the default comes from `interop.config`, read at call time not at import."""
    from chython.interop import config, conformers

    seen = []
    monkeypatch.setattr(conformers, '_cdpkit_conformers',
                        lambda *args, **kwargs: seen.append('cdpkit') or [])
    monkeypatch.setattr(conformers, '_rdkit_conformers',
                        lambda *args, **kwargs: seen.append('rdkit') or [])

    monkeypatch.setattr(config, 'conformer_engine', 'cdpkit')
    conformers.generate_conformers(_ethanol())
    monkeypatch.setattr(config, 'conformer_engine', 'rdkit')
    conformers.generate_conformers(_ethanol())

    assert seen == ['cdpkit', 'rdkit'], f'engine not taken from config at call time: {seen}'


def _ethanol():
    from chython import smiles

    return smiles('CCO')


@mark.parametrize('engine', ['rdkit', 'cdpkit'])
def test_geometry_is_three_dimensional_and_chemically_sane(engine):
    """Bond lengths land in a plausible range, so a flat or collapsed embedding fails."""
    from importlib.util import find_spec
    from pytest import skip

    if find_spec('rdkit' if engine == 'rdkit' else 'CDPL') is None:
        skip(f'{engine} is not installed')

    from chython import smiles
    from chython.interop.conformers import generate_conformers

    mol = smiles('c1ccccc1CO')  # benzyl alcohol
    assert generate_conformers(mol, limit=1, engine=engine), f'{engine} generated nothing'

    c = {n: mol.conformer(0).xyz_of(n) for n in mol.atom_numbers}
    spread = max(max(abs(a[i] - b[i]) for i in range(3)) for a in c.values() for b in c.values())
    assert spread > 1., f'{engine} produced a collapsed embedding'

    for bond in mol.bonds():
        n, m = bond.n, bond.m
        d = sum((c[n][i] - c[m][i]) ** 2 for i in range(3)) ** .5
        assert 1.1 < d < 1.7, f'{engine}: bond {n}-{m} is {d:.2f} A, not a plausible bond length'


@mark.parametrize('engine', ['rdkit', 'cdpkit'])
def test_chirality_survives_into_the_geometry(engine):
    """Two enantiomers embed with opposite handedness, measured on the coordinates directly."""
    from importlib.util import find_spec
    from pytest import skip

    if find_spec('rdkit' if engine == 'rdkit' else 'CDPL') is None:
        skip(f'{engine} is not installed')

    from chython import smiles
    from chython.interop.conformers import generate_conformers

    # Isovaline and not alanine: the maps cover heavy atoms only, so only a quaternary centre puts all
    # four neighbours in the answer.
    volumes = []
    for smi in ('CC[C@](C)(N)C(=O)O', 'CC[C@@](C)(N)C(=O)O'):  # (R)- and (S)-isovaline
        mol = smiles(smi)
        assert generate_conformers(mol, limit=1, engine=engine), f'{engine} generated nothing for {smi}'
        c = {n: mol.conformer(0).xyz_of(n) for n in mol.atom_numbers}
        # `parity_of` is the three-state read (0 unconfigured, 1 even, 2 odd); `stereo_of` is one bit
        # and answers False for the even enantiomer as well as for no parity at all.
        center = next(n for n in mol.atom_numbers if mol.parity_of(n))
        volumes.append(_signed_volume(mol, c, center))

    assert volumes[0] * volumes[1] < 0, (
        f'{engine} embedded both enantiomers with the same handedness ({volumes}); chirality was lost'
    )


@requires_cdpkit
def test_cdpkit_path_goes_through_the_shared_converter(monkeypatch):
    """
    The CDPKit molecule is built by `interop.cdpkit`, not by a second builder living here.

    Reintroducing a local builder makes this fail; a duplicate is what lets the two drift on stereo.
    """
    from chython import smiles
    from chython.interop import _cdpkit, conformers

    calls = []
    original = _cdpkit.to_cdpkit
    monkeypatch.setattr(_cdpkit, 'to_cdpkit',
                        lambda mol, **kw: (calls.append(mol) or original(mol, **kw)))

    mol = smiles('CCO')
    conformers.generate_conformers(mol, limit=1, engine='cdpkit')

    assert calls, 'the cdpkit engine did not go through interop.cdpkit'
