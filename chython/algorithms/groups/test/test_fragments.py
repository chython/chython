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
from chython import smiles
from chython.algorithms.groups import FragmentResult


def _by_role(mol, role=None):
    return list(mol.sticky_fragments(role))


def test_aryl_bromide_yields_capped_fragment():
    mol = smiles('Brc1ccc(cc1)C')
    results = _by_role(mol, 'aryl_halide')
    assert len(results) == 1
    r = results[0]
    assert isinstance(r, FragmentResult)
    assert r.role == 'aryl_halide'
    # left-aligned form: open bond first ('-...'), glues onto the left
    assert r.sticky_left.startswith('-')
    assert '[At]' not in r.sticky_left
    # right-aligned form: open bond last ('...-'), glues onto the right
    assert r.sticky_right.rstrip().endswith('-')
    assert '[At]' not in r.sticky_right
    # canonical smiles is a stable string containing the [At] cap and toluene core
    assert '[At]' in r.canonical_smiles
    assert 'C' in r.canonical_smiles


def test_left_and_right_forms_are_glueable():
    # two copies of the same fragment glue right-open onto left-open into a
    # single connected biaryl (no leftover cap, one ring-fused product).
    r = _by_role(smiles('Brc1ccc(cc1)C'), 'aryl_halide')[0]
    glued = r.sticky_right + r.sticky_left[1:]  # drop the leading '-'
    mol = smiles(glued)
    assert mol.connected_components_count == 1
    assert '[At]' not in str(mol)


def test_canonical_key_stable_across_sources():
    # same fragment reached from two different source molecules -> same key
    a = smiles('Brc1ccc(cc1)C')      # 4-bromotoluene
    b = smiles('Ic1ccc(cc1)C')       # 4-iodotoluene
    ka = _by_role(a, 'aryl_halide')[0].canonical_smiles
    kb = _by_role(b, 'aryl_halide')[0].canonical_smiles
    assert ka == kb


def test_role_none_runs_all_roles():
    # bromobenzoic acid has both an aryl_halide and an aryl_acyl handle
    mol = smiles('Brc1ccc(cc1)C(=O)O')
    roles_found = {r.role for r in _by_role(mol)}
    assert 'aryl_halide' in roles_found
    assert 'aryl_acyl' in roles_found


def test_no_match_yields_nothing():
    mol = smiles('CCCC')
    assert _by_role(mol, 'aryl_halide') == []


def test_salt_input_does_not_abort_enumeration():
    # sodium 4-bromobenzoate: the Na+ counter-ion is a separate component, so
    # every cut yields a disconnected product. These must be skipped rather than
    # raising out of the generator (a raise would drop the whole enumeration).
    mol = smiles('[Na+].[O-]C(=O)c1ccc(Br)cc1')
    assert _by_role(mol) == []  # skipped cleanly, no exception


def test_unknown_role_raises():
    import pytest
    mol = smiles('Brc1ccccc1')
    with pytest.raises(ValueError):
        list(mol.sticky_fragments('not_a_role'))
