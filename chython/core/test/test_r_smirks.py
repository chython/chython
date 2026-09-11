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
"""`#0` is the R marker, and it is a SMIRKS product-side BUILD spelling.

One lexer serves both sides of the arrow, so the token READS on either side and is refused where it
would have to compile into a box: an R matches nothing, and the product side is the one side that
never seals.  What it buys is a template that states its own attachment point -- the centre carrying
the cap exists inside the patch, so `@=` and `@~` on it are statements the patcher can apply and the
arena's re-base carries every other configuration across.
"""
from pytest import mark, raises
from chython.core import IncorrectSmarts, read_smarts, read_smiles, read_smirks


@mark.parametrize('pattern', ['[#0]', '[C;#0]', '[!#0]', '[#0]C', '[#0,N]'])
def test_a_query_holding_the_marker_is_refused(pattern):
    """The refusal is the seal's, so it names the atom rather than a byte offset."""
    with raises(IncorrectSmarts, match='matches nothing'):
        read_smarts(pattern)


def test_the_reactant_side_is_refused_and_names_the_side():
    with raises(IncorrectSmarts, match='the reactant side does not compile'):
        read_smirks('[C:1][#0:2]>>[C:1]')


def test_the_marker_is_not_an_element_and_takes_no_isotope():
    """`[13C]` reads and `[13#0]` does not: the digits have nothing to attach to."""
    with raises(IncorrectSmarts, match='no element symbol'):
        read_smarts('[13#0]')


@mark.parametrize('smirks,message', [
    ('[C:1][Br;D1]>>[C:1][#0,N]', 'as one of several alternatives'),
    ('[C:1][Br;D1]>>[C:1][!#0]', 'names nothing to build'),
    ('[C:1][Br;D1]>>[C:1][C;#0]', 'states the element twice'),
])
def test_the_product_side_refusals_are_the_element_field_s(smirks, message):
    """`#0` IS the element field, so it collides with `C` and with `A` exactly as two elements do."""
    with raises(IncorrectSmarts, match=message):
        read_smirks(smirks)


def _isopropyl():
    mol = read_smiles('CC(C)Br')
    mol.canonicalize()
    return mol


def test_a_created_marker_is_the_attachment_point():
    """The cap is built by the template, not bolted on afterwards, and its neighbour's hydrogen count
    comes out of the ordinary recompute -- an R reads as carbon there."""
    product = next(iter(read_smirks('[C:1][Br;D1]>>[C:1][#0]')(_isopropyl()))).products[0]
    marker = next(n for n in product if product.atom(n).is_r)
    assert product.element_of(marker) == 0
    assert product.atom(marker).r_index == 0
    assert product.implicit_h_of(marker) == 0
    site = next(iter(product.neighbors_of(marker)))
    assert product.implicit_h_of(site) == 1
    assert str(product) == '[R]C(C)C'


def test_a_paired_marker_keeps_the_leaving_group_s_id():
    """`[Br:2]>>[#0:2]` turns the matched atom INTO the marker, so the id a caller was holding
    survives -- which is the whole reason `set_element` is on the mutation surface."""
    template = read_smirks('[C:1][Br;D1:2]>>[C:1][#0:2]')
    reaction, where = next(iter(template(_isopropyl(), report=True)))
    product = reaction.products[0]
    assert product.atom(where[2]).is_r
    assert product.implicit_h_of(where[2]) == 0
    assert product.implicit_h_of(where[1]) == 1


@mark.parametrize('smirks,expected', [
    ('[C:1][Br;D1]>>[C@=:1][#0]', 'CC[C@@H]([R])C'),       # `@=`: a cut keeps what it cut
    ('[C;@:1][Br;D1]>>[C@=:1][#0]', 'CC[C@@H]([R])C'),     # and a reactant sign only narrows the match
    ('[C:1][Br;D1]>>[C@~:1][#0]', 'CC[C@H]([R])C'),        # `@~`: the other configuration
    ('[C:1][Br;D1]>>[C:1][#0]', 'CCC([R])C'),              # unstated: the reaction centre's drop
])
def test_the_capped_centre_keeps_its_configuration(smirks, expected):
    """The centre exists in the patch, so its configuration is the patcher's business: a template
    states retention or inversion at the atom the cap hangs off, and states `@=` where a cut keeps
    whatever it found -- which is what every `roles.tsv` row says."""
    mol = read_smiles('C[C@H](Br)CC')
    mol.canonicalize()
    product = next(iter(read_smirks(smirks)(mol))).products[0]
    assert str(product) == expected


def test_a_cis_trans_unit_survives_a_cap():
    """The cut is allylic here, so the double bond is no part of the reaction centre and the arena's
    re-base stands with nothing stated.  A cut AT the alkene needs `@=`, which is the kind's only
    spelling."""
    mol = read_smiles('C/C=C/CBr')
    mol.canonicalize()
    product = next(iter(read_smirks('[C:1][Br;D1]>>[C:1][#0]')(mol))).products[0]
    assert str(product) == 'C(=C\\C)/C[R]'

    mol = read_smiles('C/C=C/Br')
    mol.canonicalize()
    assert str(next(iter(read_smirks('[C:1][Br;D1]>>[C:1][#0]')(mol))).products[0]) == 'C(C)=C[R]'
    assert str(next(iter(read_smirks('[C:1][Br;D1]>>[C@=:1][#0]')(mol))).products[0]) == 'C(/C)=C\\[R]'
