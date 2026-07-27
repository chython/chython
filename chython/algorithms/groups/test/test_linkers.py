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
from chython.algorithms.groups import LinkerResult


def test_bifunctional_yields_dual_capped_linker():
    # 4-bromobenzoic acid: aryl_halide on the ring, aryl_acyl on the COOH
    mol = smiles('Brc1ccc(cc1)C(=O)O')
    results = list(mol.sticky_linkers('aryl_halide', 'aryl_acyl'))
    assert results, 'expected at least one linker'
    r = results[0]
    assert isinstance(r, LinkerResult)
    assert r.role_left == 'aryl_halide'
    assert r.role_right == 'aryl_acyl'
    # sticky_smiles is open-bond at both ends ('-...-'): left end (210) first,
    # right end (211) last, ready to concatenate on either side.
    assert r.sticky_smiles.startswith('-')
    assert r.sticky_smiles.rstrip().endswith('-')
    assert '[At]' not in r.sticky_smiles
    # isotope caps live in the canonical dedup key, ordered 210 (left) / 211 (right)
    assert '[210At]' in r.canonical_smiles and '[211At]' in r.canonical_smiles


def test_ordering_controls_isotope_assignment():
    mol = smiles('Brc1ccc(cc1)C(=O)O')
    forward = list(mol.sticky_linkers('aryl_halide', 'aryl_acyl'))[0]
    reverse = list(mol.sticky_linkers('aryl_acyl', 'aryl_halide'))[0]
    assert forward.role_left == 'aryl_halide' and forward.role_right == 'aryl_acyl'
    assert reverse.role_left == 'aryl_acyl' and reverse.role_right == 'aryl_halide'
    # different ordering -> different canonical key (accepted per spec)
    assert forward.canonical_smiles != reverse.canonical_smiles


def test_both_none_runs_all_role_pairs():
    mol = smiles('Brc1ccc(cc1)C(=O)O')
    pairs = {(r.role_left, r.role_right) for r in mol.sticky_linkers()}
    assert ('aryl_halide', 'aryl_acyl') in pairs


def test_monofunctional_yields_no_linker():
    mol = smiles('Brc1ccccc1')  # only one handle
    assert list(mol.sticky_linkers('aryl_halide', 'aryl_acyl')) == []


def test_no_single_atom_linker_from_overlapping_cuts():
    # bromoacetic acid: alkyl_halide caps C2 (drops Br) and alkyl_decarboxy also
    # anchors on C2 (drops the whole COOH). Both cuts on one atom would collapse
    # to [At]C[At] -- a degenerate zero-span linker that must be suppressed.
    for smi in ('BrCC(=O)O', 'OCC(=O)O', 'OC(=O)CC(=O)O'):
        for r in smiles(smi).sticky_linkers():
            core = r.canonical_smiles.replace('[210At]', '').replace('[211At]', '')
            heavy = ''.join(c for c in core if c.isalpha())
            assert heavy != 'C', f'{smi} produced a single-atom linker: {r.canonical_smiles}'


def test_unknown_role_raises():
    import pytest
    mol = smiles('Brc1ccc(cc1)C(=O)O')
    with pytest.raises(ValueError):
        list(mol.sticky_linkers('aryl_halide', 'nope'))
    with pytest.raises(ValueError):
        list(mol.sticky_linkers('nope', 'aryl_acyl'))
