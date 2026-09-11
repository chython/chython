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
"""Sticky fragments: the R-capped cut, its two open-bond spellings, and the glue rule."""
from pytest import mark, raises
from ...core import read_smiles as smiles


def test_an_aryl_bromide_yields_one_fragment():
    mol = smiles('Brc1ccccc1')
    mol.canonicalize()
    out = list(mol.sticky_fragments('aryl_halide'))
    assert len(out) == 1
    assert out[0].role == 'aryl_halide'


def test_the_dedup_key_carries_the_marker():
    mol = smiles('Brc1ccccc1')
    mol.canonicalize()
    fragment = next(iter(mol.sticky_fragments('aryl_halide')))
    assert '[R]' in fragment.canonical_smiles
    assert 'Br' not in fragment.canonical_smiles


def test_the_two_spellings_carry_the_bond_on_the_left():
    mol = smiles('Brc1ccccc1')
    mol.canonicalize()
    fragment = next(iter(mol.sticky_fragments('aryl_halide')))
    assert fragment.sticky_left.startswith('-')
    assert not fragment.sticky_right.endswith('-')


def test_glueing_two_fragments_gives_the_coupled_product():
    left = smiles('Brc1ccccc1')
    left.canonicalize()
    right = smiles('OB(O)c1ccncc1')
    right.canonicalize()
    a = next(iter(left.sticky_fragments('aryl_halide')))
    b = next(iter(right.sticky_fragments('aryl_boron')))
    glued = smiles(a.sticky_right + b.sticky_left)
    glued.canonicalize()
    expected = smiles('c1ccc(-c2ccncc2)cc1')
    expected.canonicalize()
    assert glued == expected


def test_every_role_is_enumerated_when_none_is_named():
    mol = smiles('Brc1ccc(C(=O)O)cc1')
    mol.canonicalize()
    found = {f.role for f in mol.sticky_fragments()}
    assert 'aryl_halide' in found
    assert 'aryl_acyl' in found


def test_an_unknown_role_is_refused():
    mol = smiles('Brc1ccccc1')
    mol.canonicalize()
    with raises(ValueError, match='role'):
        list(mol.sticky_fragments('not_a_role'))


def test_a_mixture_yields_nothing():
    mol = smiles('Brc1ccccc1.O')
    mol.canonicalize()
    assert not list(mol.sticky_fragments())


def test_masked_bars_the_attachment_site():
    mol = smiles('Brc1ccccc1')
    mol.canonicalize()
    site = next(n for n in mol if mol.atom(n).atomic_symbol == 'C'
                and any(mol.atom(x).atomic_symbol == 'Br' for x in mol.neighbors_of(n)))
    assert not list(mol.sticky_fragments('aryl_halide', masked=[site]))


def test_masked_bars_a_consumed_leaving_group():
    # `alkyl_deamino` consumes the nitrogen; masking it bars the role.
    mol = smiles('NCc1ccccc1')
    mol.canonicalize()
    nitrogen = next(n for n in mol if mol.atom(n).atomic_symbol == 'N')
    assert not list(mol.sticky_fragments('alkyl_deamino', masked=[nitrogen]))
    assert list(mol.sticky_fragments('alkyl_deamino'))


def test_the_capped_neighbour_reads_the_r_as_carbon():
    mol = smiles('Brc1ccccc1')
    mol.canonicalize()
    fragment = next(iter(mol.sticky_fragments('aryl_halide')))
    capped = smiles(fragment.canonical_smiles)
    r = next(a for a in capped.atoms() if a.is_r)
    carrier = next(capped.atom(n) for n in capped.neighbors_of(
        next(n for n in capped if capped.atom(n).is_r)))
    assert carrier.implicit_h == 0
    assert carrier.heteroatoms == 0
    assert r.is_r


def test_a_bifunctional_molecule_yields_a_linker():
    mol = smiles('Brc1ccc(C(=O)O)cc1')
    mol.canonicalize()
    out = [x for x in mol.sticky_linkers('aryl_halide', 'aryl_acyl')]
    assert out
    assert out[0].role_left == 'aryl_halide'
    assert out[0].role_right == 'aryl_acyl'


def test_the_linker_key_indexes_left_as_one_and_right_as_two():
    mol = smiles('Brc1ccc(C(=O)O)cc1')
    mol.canonicalize()
    linker = next(iter(mol.sticky_linkers('aryl_halide', 'aryl_acyl')))
    assert '[R1]' in linker.canonical_smiles
    assert '[R2]' in linker.canonical_smiles


def test_the_two_indices_are_not_interchangeable():
    mol = smiles('Brc1ccc(C(=O)O)cc1')
    mol.canonicalize()
    forward = next(iter(mol.sticky_linkers('aryl_halide', 'aryl_acyl')))
    reverse = next(iter(mol.sticky_linkers('aryl_acyl', 'aryl_halide')))
    assert forward.canonical_smiles != reverse.canonical_smiles


def test_both_spellings_are_open_at_both_ends():
    mol = smiles('Brc1ccc(C(=O)O)cc1')
    mol.canonicalize()
    linker = next(iter(mol.sticky_linkers('aryl_halide', 'aryl_acyl')))
    for spelling in (linker.sticky_left, linker.sticky_right):
        assert spelling.startswith('-')
        assert not spelling.endswith('-')


def test_a_fragment_linker_fragment_chain_re_reads():
    left = smiles('Brc1ccccc1')
    left.canonicalize()
    middle = smiles('Brc1ccc(C(=O)O)cc1')
    middle.canonicalize()
    right = smiles('NCc1ccccc1')
    right.canonicalize()
    a = next(iter(left.sticky_fragments('aryl_halide')))
    linker = next(iter(middle.sticky_linkers('aryl_halide', 'aryl_acyl')))
    b = next(iter(right.sticky_fragments('alkyl_amine')))
    chain = smiles(a.sticky_right + linker.sticky_left + b.sticky_left)
    assert chain.atom_count > left.atom_count


def test_two_caps_on_the_same_atom_are_skipped():
    # Bromoacetic acid: the halide handle and the decarboxylative handle share the methylene.
    mol = smiles('OC(=O)CBr')
    mol.canonicalize()
    assert not list(mol.sticky_linkers('alkyl_halide', 'alkyl_decarboxy'))


def test_masked_applies_to_the_left_end_only():
    mol = smiles('Brc1ccc(N)cc1')
    mol.canonicalize()
    nitrogen = next(n for n in mol if mol.atom(n).atomic_symbol == 'N')
    assert not list(mol.sticky_linkers('aryl_amine', 'aryl_halide', masked=[nitrogen]))
    assert list(mol.sticky_linkers('aryl_halide', 'aryl_amine', masked=[nitrogen]))


def test_a_mixture_yields_no_linker():
    mol = smiles('Brc1ccc(C(=O)O)cc1.O')
    mol.canonicalize()
    assert not list(mol.sticky_linkers())


def test_a_fragment_cannot_be_used_as_a_query():
    """The door a stickers caller walks into: an R matches nothing, so `as_query` refuses.

    Lives here rather than beside the other refusals in `chython/core/test/test_r_query.py`, because
    `sticky_fragments` is injected onto the container by this package and that directory may import
    nothing from this distribution but `chython.core`.
    """
    mol = smiles('Brc1ccccc1')
    mol.canonicalize()
    fragment = next(iter(mol.sticky_fragments('aryl_halide')))
    with raises(ValueError, match='matches nothing'):
        smiles(fragment.canonical_smiles) <= smiles('Cc1ccccc1')


@mark.parametrize('salt, role, left, right, fragment', [
    ('[K+].[B-](F)(F)(F)c1ccccc1', 'aryl_boron',    '-c1ccccc1', 'c(cccc1)c1', 'c1c([R])cccc1'),
    ('[K+].[B-](F)(F)(F)CCCC',     'alkyl_boron',   '-CCCC',     'C(C)CC',     'C(C)CC[R]'),
    ('[K+].[B-](F)(F)(F)C=C',      'alkenyl_boron', '-C=C',      'C(=C)',      '[R]C=C'),
    ('[K+].[B-](F)(F)(F)C#CC',     'alkynyl_boron', '-C#CC',     'C(C)#C',     'C(C)#C[R]'),
])
def test_a_molander_salt_enumerates_once_its_counter_ion_is_gone(salt, role, left, right, fragment):
    """The four `*_molander_salt` rows are reachable only through a caller who splits the salt.

    They are also the corpus's only multi-product patches, so this is the enumerator's live exercise of
    picking the product that HOLDS the site rather than the first one.  `decompose_salts` is the name
    for reducing a salt to its compound; a component walk stands in for it here, this suite being barred
    from a sibling package.
    """
    mol = smiles(salt)
    mol.canonicalize()
    assert not list(mol.sticky_fragments(role))

    anion = next(c for c in mol.split() if c.connected_components_count == 1 and len(c) > 1)
    anion.canonicalize()
    out = list(anion.sticky_fragments(role))
    assert len(out) == 1
    assert out[0].role == role
    assert out[0].sticky_left == left
    assert out[0].sticky_right == right
    assert out[0].canonical_smiles == fragment
