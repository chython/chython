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
from chython.algorithms.groups import StickyLinker


def test_bifunctional_yields_dual_capped_linker():
    # 4-bromobenzoic acid: aryl_halide on the ring, aryl_acyl on the COOH
    mol = smiles('Brc1ccc(cc1)C(=O)O')
    results = list(mol.sticky_linkers('aryl_halide', 'aryl_acyl'))
    assert results, 'expected at least one linker'
    r = results[0]
    assert isinstance(r, StickyLinker)
    assert r.role_left == 'aryl_halide'
    assert r.role_right == 'aryl_acyl'
    # two open-bond traversals ('-...-'): sticky_left leads with the role_left end,
    # sticky_right leads with the role_right end (the flipped orientation).
    for s in (r.sticky_left, r.sticky_right):
        assert s.startswith('-') and s.rstrip().endswith('-')
        assert '[At]' not in s
    assert r.sticky_left != r.sticky_right          # genuinely different orientations
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


def test_deamino_acyl_linker_from_aminobenzoic_acid():
    # 4-aminobenzoic acid: aryl_deamino (NH2 via diazonium) + aryl_acyl (COOH).
    mol = smiles('Nc1ccc(cc1)C(=O)O')
    results = list(mol.sticky_linkers('aryl_deamino', 'aryl_acyl'))
    assert results, 'expected at least one linker'
    r = results[0]
    assert isinstance(r, StickyLinker)
    assert r.role_left == 'aryl_deamino' and r.role_right == 'aryl_acyl'
    for s in (r.sticky_left, r.sticky_right):
        assert s.startswith('-') and s.rstrip().endswith('-')
    assert '[210At]' in r.canonical_smiles and '[211At]' in r.canonical_smiles


def test_alkyl_deamino_halide_linker():
    # 4-bromobenzylamine: alkyl_deamino (CH2NH2, drops N) + aryl_halide (Br).
    mol = smiles('NCc1ccc(Br)cc1')
    results = list(mol.sticky_linkers('alkyl_deamino', 'aryl_halide'))
    assert results, 'expected at least one linker'
    r = results[0]
    assert r.role_left == 'alkyl_deamino' and r.role_right == 'aryl_halide'
    # nitrogen dropped on the deaminative end
    assert 'N' not in r.canonical_smiles
    assert '[210At]' in r.canonical_smiles and '[211At]' in r.canonical_smiles


def test_unknown_role_raises():
    import pytest
    mol = smiles('Brc1ccc(cc1)C(=O)O')
    with pytest.raises(ValueError):
        list(mol.sticky_linkers('aryl_halide', 'nope'))
    with pytest.raises(ValueError):
        list(mol.sticky_linkers('nope', 'aryl_acyl'))


def test_sticky_left_and_right_are_reversed_traversals():
    # both traversals open the bond at BOTH ends ('-...-'); sticky_left leads with
    # the role_left (210) end, sticky_right leads with the role_right (211) end.
    mol = smiles('Brc1ccc(cc1)C(=O)O')
    r = list(mol.sticky_linkers('aryl_halide', 'aryl_acyl'))[0]
    # halide end is the aromatic ring carbon; acyl end is the carbonyl carbon.
    # sticky_left starts on the ring, sticky_right starts on the C=O.
    assert r.sticky_left.startswith('-c') and r.sticky_left.rstrip().endswith('-')
    assert r.sticky_right.startswith('-C(=O)') and r.sticky_right.rstrip().endswith('-')


def test_masked_atom_barred_from_left_but_allowed_on_right():
    # Boc-benzylamine: reveal the amine, then mask it. As a "cleaved FG" it may
    # only be the deferred step-2 (right, 211) end of a linker, never step-1
    # (left, 210) -- which also kills the (amine_left,X)/(X,amine_right) ordering
    # duplicate.
    mol = smiles('O=C(OC(C)(C)C)NCc1ccc(Br)cc1')
    mol.canonicalize()
    freed = mol.remove_protection(logging=True)
    unmasked = {(l.role_left, l.role_right) for l in mol.sticky_linkers()}
    masked = {(l.role_left, l.role_right) for l in mol.sticky_linkers(masked=freed)}
    # unmasked surfaces the amine at BOTH orientations
    assert ('alkyl_amine', 'aryl_halide') in unmasked
    assert ('aryl_halide', 'alkyl_amine') in unmasked
    # masked: amine barred as LEFT, kept as RIGHT
    assert ('alkyl_amine', 'aryl_halide') not in masked
    assert ('aryl_halide', 'alkyl_amine') in masked


def test_masked_none_and_empty_are_equivalent():
    mol = smiles('Brc1ccc(cc1)C(=O)O')
    default = {l.canonical_smiles for l in mol.sticky_linkers()}
    explicit = {l.canonical_smiles for l in mol.sticky_linkers(masked=None)}
    empty = {l.canonical_smiles for l in mol.sticky_linkers(masked=())}
    assert default == explicit == empty
