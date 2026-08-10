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
from chython.algorithms.groups import StickyFragment


def _by_role(mol, role=None):
    return list(mol.sticky_fragments(role))


def test_aryl_bromide_yields_capped_fragment():
    mol = smiles('Brc1ccc(cc1)C')
    results = _by_role(mol, 'aryl_halide')
    assert len(results) == 1
    r = results[0]
    assert isinstance(r, StickyFragment)
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


def test_alkyl_deamino_drops_nitrogen_caps_carbon():
    # deaminative coupling handle: propylamine -> propyl cap (N gone, alpha C capped)
    results = _by_role(smiles('NCCC'), 'alkyl_deamino')
    assert len(results) == 1
    r = results[0]
    assert isinstance(r, StickyFragment)
    assert r.role == 'alkyl_deamino'
    # no nitrogen left in the capped fragment
    assert 'N' not in r.canonical_smiles
    assert '[At]' in r.canonical_smiles
    assert r.sticky_left.startswith('-') and r.sticky_right.rstrip().endswith('-')


def test_aryl_deamino_from_aniline():
    # diazonium-route deaminative handle: aniline -> aryl cap (N gone)
    results = _by_role(smiles('Nc1ccc(cc1)C'), 'aryl_deamino')
    assert len(results) == 1
    r = results[0]
    assert r.role == 'aryl_deamino'
    assert 'N' not in r.canonical_smiles
    assert '[At]' in r.canonical_smiles


def test_tertiary_alcohol_deoxy_fragment():
    # tert-butanol under the deoxygenative handle: O dropped, quaternary C capped
    results = _by_role(smiles('CC(C)(O)C'), 'alkyl_deoxy')
    assert len(results) == 1
    r = results[0]
    assert r.role == 'alkyl_deoxy'
    assert 'O' not in r.canonical_smiles
    assert '[At]' in r.canonical_smiles


def test_amine_carries_both_nucleophile_and_deamino_roles():
    # a primary alkyl amine exposes an N-nucleophile handle (alkyl_amine, keeps N)
    # AND a deaminative C-handle (alkyl_deamino, drops N). role=None yields both.
    roles_found = {r.role for r in _by_role(smiles('NCCC'))}
    assert 'alkyl_amine' in roles_found
    assert 'alkyl_deamino' in roles_found


def test_carbamoyl_halide_keeps_nitrogen_and_carbonyl():
    # carbamoyl chloride (R2N-C(=O)-Cl): the halide leaves, the carbonyl carbon is
    # capped, and both the amine N and the carbonyl O are retained.
    results = _by_role(smiles('ClC(=O)N1CCCCC1'), 'carbamoyl_halide')
    assert len(results) == 1
    r = results[0]
    assert isinstance(r, StickyFragment)
    assert r.role == 'carbamoyl_halide'
    assert 'N' in r.canonical_smiles and 'O' in r.canonical_smiles
    assert 'Cl' not in r.canonical_smiles
    assert '[At]' in r.canonical_smiles


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


def test_masked_atom_barred_as_handle_and_as_leaving_group():
    # Boc-benzylamine: strip the Boc, then bar the just-revealed amine N (a
    # "cleaved FG") from driving any coupling this fragment represents. A masked
    # atom participates in two ways, both barred:
    #   - as the attachment site: alkyl_amine caps the N itself -> vanishes;
    #   - as a leaving group: alkyl_deamino consumes the N while capping the
    #     alpha-carbon -> also vanishes (the deaminative coupling would still be
    #     spending the freshly-revealed amine).
    # aryl_halide neither caps nor consumes the amine, so it survives.
    mol = smiles('O=C(OC(C)(C)C)NCc1ccc(Br)cc1')
    mol.canonicalize()
    freed = mol.remove_protection(logging=True)
    assert freed, 'expected the Boc removal to free the amine nitrogen'
    unmasked = {f.role for f in mol.sticky_fragments()}
    masked = {f.role for f in mol.sticky_fragments(masked=freed)}
    assert {'alkyl_amine', 'alkyl_deamino', 'aryl_halide'} <= unmasked
    assert 'alkyl_amine' not in masked          # freed amine barred as a handle
    assert 'alkyl_deamino' not in masked        # ...and as a leaving group
    assert 'aryl_halide' in masked              # untouched role survives


def test_masked_none_and_empty_are_equivalent():
    mol = smiles('Nc1ccc(Br)cc1')
    default = {f.canonical_smiles for f in mol.sticky_fragments()}
    explicit = {f.canonical_smiles for f in mol.sticky_fragments(masked=None)}
    empty = {f.canonical_smiles for f in mol.sticky_fragments(masked=())}
    assert default == explicit == empty
