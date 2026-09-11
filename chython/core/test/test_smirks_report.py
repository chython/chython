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
"""`report=True` says which atom id each product-side map number landed on.

Needed because a template's map numbers and the reaction's imposed mapping are unrelated: an enumerator
that must edit "the atom the product side called :1" cannot find it by map number afterwards.
"""
from chython.core import read_smirks, read_smiles as smiles


def test_report_names_every_product_map_number():
    # Hydrolysis of an acyl chloride: the product side states the carbon and its oxygen.
    template = read_smirks('[Cl;D1][C:1]=[O:2]>>[C:1]=[O:2]')
    reaction, where = next(template(smiles('CC(=O)Cl'), report=True))
    assert sorted(where) == [1, 2]
    product = next(iter(reaction.products))
    assert product.atom(where[1]).atomic_symbol == 'C'
    assert product.atom(where[2]).atomic_symbol == 'O'


def test_without_report_the_yield_is_unchanged():
    template = read_smirks('[Cl;D1][C:1]=[O:2]>>[C:1]=[O:2]')
    reaction = next(template(smiles('CC(=O)Cl')))
    assert len(reaction.products) == 1


def test_report_survives_a_created_atom():
    template = read_smirks('[Cl;D1][C:1]=[O:2]>>[O:3][C:1]=[O:2]')
    reaction, where = next(template(smiles('CC(=O)Cl'), report=True))
    product = next(iter(reaction.products))
    assert product.atom(where[3]).atomic_symbol == 'O'


def test_the_reported_id_is_usable_in_an_edit_session():
    template = read_smirks('[Cl;D1][C:1]=[O:2]>>[C:1]=[O:2]')
    reaction, where = next(template(smiles('CC(=O)Cl'), report=True))
    product = next(iter(reaction.products))
    with product.edit() as e:
        e.set_map_number(where[1], 77)
    assert product.atom(where[1]).map_number == 77


def test_a_reported_id_lives_in_exactly_one_product_component():
    template = read_smirks('[Cl;D1][C:1]=[O:2]>>[C:1]=[O:2]')
    reaction, where = next(template(smiles('CC(=O)Cl.[Na+].[Cl-]'), report=True))
    holders = [p for p in reaction.products if where[1] in p.atom_numbers]
    assert len(holders) == 1
    assert holders[0].atom(where[1]).atomic_symbol == 'C'


def test_an_unmapped_product_atom_is_not_in_the_report():
    # The report is keyed by map number, so an atom the product side never numbered has no key.  It is
    # still built -- the product gains it -- and a caller that must reach it numbers it in the template.
    template = read_smirks('[C:1]=[O:2]>>[C:1]([O:2])O')
    reaction, where = next(template(smiles('CC=O'), report=True))
    assert sorted(where) == [1, 2]
    assert next(iter(reaction.products)).atom_count == 4
