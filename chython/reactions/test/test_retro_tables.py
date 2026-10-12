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
"""`functional.tsv`'s `precursor` column and the reverse templates composed from it.

A spelling is the concrete group `retro()` rebuilds; the gates here hold it to the group it spells.
"""
import re
import pytest
from .._tables import compose_retro, functional_rules, reaction_rules, retro_templates
from ...core import read_smiles, read_smirks


REFERENCED = sorted({name for family in reaction_rules().values()
                     for rule in family for name in rule.groups})


def _normalized(molecule):
    molecule = molecule.copy()
    molecule.kekule()
    molecule.thiele()
    return molecule


def test_every_referenced_group_has_a_spelling():
    """`retro_templates()` refuses a row whose group has none; this keeps the refusal out of reach."""
    missing = [name for name in REFERENCED if not functional_rules()[name].precursors]
    assert not missing, f'functional.tsv groups a reactions.tsv row uses with no `precursor`: {missing}'


@pytest.mark.parametrize('name', REFERENCED)
def test_a_spelling_rebuilds_its_own_group(name):
    """`<group smarts>>><spelling>` on the group's `example` is still the group.

    The products are unioned first: a salt's counter-ion is a separate component of the answer.
    """
    group = functional_rules()[name]
    example = _normalized(read_smiles(group.example))
    for spelling in group.precursors:
        patch = re.sub(r';h\d+', '', spelling)    # a patch `h` is a post-check, see `RetroTemplate`
        template = read_smirks('%s>>%s' % (group.smarts, patch))
        rebuilt = []
        for rxn in template(example):
            product = rxn.products[0]
            for part in rxn.products[1:]:
                product = product.union(part)
            rebuilt.append(_normalized(product))
        assert rebuilt, f'{group.id} ({name}): {spelling!r} does not apply to its own example {group.example!r}'
        assert any(group.query.is_substructure(m) for m in rebuilt), (
            f'{group.id} ({name}): {spelling!r} rebuilds {[format(m) for m in rebuilt]}, which is no longer '
            f'{name}.  A spelling must be an instance of the group it spells.')


@pytest.mark.parametrize('name', REFERENCED)
def test_a_spelling_names_an_element_unless_the_group_does_not(name):
    """`[A:n]` takes its element from the target, which a re-created atom does not have."""
    group = functional_rules()[name]
    alternatives = {int(m.group(2)) for m in re.finditer(r'\[([^\]:]*):(\d+)\]', group.smarts)
                    if m.group(1).split(';')[0] in ('A', '*') or ',' in m.group(1).split(';')[0]}
    for spelling in group.precursors:
        for m in re.finditer(r'\[A[;\]:][^\]]*?:(\d+)\]', spelling):
            assert int(m.group(1)) in alternatives, (
                f'{group.id} ({name}): {spelling!r} writes `A` for :{m.group(1)}, whose group atom names '
                'one element.  State it.')


RULES = {rule.id: rule for family in reaction_rules().values() for rule in family}


def _composed(rule_id):
    rule = RULES[rule_id]
    return compose_retro(rule.groups, rule.product)


def test_the_product_column_is_the_query_and_the_spellings_are_the_patch():
    """`[A:n]` takes slot n's element; a deleted atom the row re-creates loses its number."""
    smirks, hydrogens = _composed('reactions:1')[0]         # amidation: carboxylic_acid + primary_amine
    assert smirks == '[C:1](=[O:2])-[N:101]-[C:102]>>[O]-[C:1]=[O:2].[N:101]-[C:102]'
    assert hydrogens == ()


def test_one_template_per_spelling_combination_in_slot_order():
    composed = [smirks for smirks, _ in _composed('reactions:8')]    # suzuki: aryl_halide + aryl_boronic_acid
    assert [s.partition('>>')[2].partition('-')[0] for s in composed] == ['[Cl]', '[Br]', '[I]']


def test_a_relative_stereo_token_moves_from_the_query_to_the_patch():
    """A query refuses `@~`; it is its own inverse, so the reverse patch states it."""
    smirks, _ = _composed('reactions:17')[0]                 # williamson: alkyl_halide + phenol
    query, _, patch = smirks.partition('>>')
    assert '@' not in query
    assert patch.startswith('[Cl]-[C@~:1].')


def test_a_bond_between_two_ring_pinned_atoms_also_matches_aromatic():
    """A `product` column spells a new ring Kekule; `thiele()` aromatizes it in the target."""
    smirks, _ = _composed('reactions:393')[0]                # mitsunobu: methanol + azinone
    assert smirks.partition('>>')[0] == '[C:2]-[O:101]-[C:102](-,:[C,N:103])=,:[N:104]'


def test_a_bond_to_an_unpinned_atom_stays_literal():
    """Widening every bond lets a coupling row match every aromatic C-C bond of the target."""
    assert _composed('reactions:8')[0][0].partition('>>')[0] == '[C:1]-[C:101]'


def test_a_spelled_hydrogen_count_leaves_the_patch_and_is_reported():
    smirks, hydrogens = _composed('reactions:181')[0]        # carbamoylation: imidazolyl_carbonate + imidazole
    assert ';h' not in smirks
    assert hydrogens == ((101, 1), (103, 0))


def test_a_group_without_a_spelling_is_refused_naming_it():
    bare = next(name for name, group in functional_rules().items()
                if not group.precursors and ':1]' in group.smarts)
    with pytest.raises(ValueError, match=bare):
        compose_retro((bare,), '[A:1]')


def test_every_row_composes_and_reads():
    templates = retro_templates()
    assert set(templates) == set(RULES)
    assert all(templates[rule_id] for rule_id in RULES)
    assert all(t.template.rule_id == rule_id for rule_id, ts in templates.items() for t in ts)


def test_the_reverse_templates_are_cached():
    assert retro_templates() is retro_templates()
