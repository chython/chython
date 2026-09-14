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
"""Completing a one-coordinate organozinc or Grignard: `[Zn+].[Cl-]` joins, and a halide-less one charges.

The invariant under test is that THE RESULT DOES NOT DEPEND ON INPUT ATOM ORDER.  Two orderings settle
it -- halides by `Cl > Br > I > F`, metals by canonical rank -- and every test that could see an order
effect permutes its input and demands one answer.  Charging needs no order, so what is tested there is
that net charge moves, that it moves once, and that it is a record of its own.
"""
# `__all__` and not the package object: `from ... import chemistry` would execute the facade, which
# `test_dependency_direction.py` ratchets against.
from itertools import permutations
from .. import canonicalize, standardize
from ...core import read_smiles as smiles


def _std(s):
    """`standardize()` then a canonical string, so two spellings of one answer compare equal."""
    m = smiles(s)
    standardize(m)
    m.canonicalize()
    return format(m, '')


def _one(spellings):
    """The single canonical answer every spelling gives, or an assertion naming the ones that differ."""
    got = {s: _std(s) for s in spellings}
    assert len(set(got.values())) == 1, got
    return next(iter(got.values()))


# what gets joined


def test_every_charge_spelling_of_a_drawn_apart_reagent_reaches_the_covalent_form():
    """All four spellings are the same reagent: ethylzinc chloride, and phenylmagnesium bromide.

    Two of the four do not conserve net charge -- a lone `[Zn+]` beside a neutral halogen is at `+1`, a
    neutral metal beside `[Cl-]` at `-1` -- and both are joined anyway, the premise being that the
    drawing dropped a sign rather than that a cation and a radical were meant.
    """
    assert _one(['CC[Zn+].[Cl-]', 'CC[Zn+].[Cl]', 'CC[Zn].[Cl-]', 'CC[Zn].[Cl]']) == _std('CC[Zn]Cl')
    assert _one(['c1ccccc1[Mg+].[Br-]', 'c1ccccc1[Mg].[Br-]']) == _std('c1ccccc1[Mg]Br')


def test_all_four_halogens_are_joined():
    for x in ('F', 'Cl', 'Br', 'I'):
        assert _std(f'CC[Zn+].[{x}-]') == _std(f'CC[Zn]{x}'), x


def test_the_joined_form_is_a_fixed_point_and_a_plain_molecule_is_untouched():
    assert _std('CC[Zn]Cl') == _std(_std('CC[Zn]Cl'))
    m = smiles('c1ccccc1O')
    assert not standardize(m)


# what is left alone


def test_only_a_one_bonded_metal_holding_a_carbon_takes_a_halide():
    for s in ('C[Zn](C)C.[Cl-]',        # already two-coordinate: no room
              'CO[Zn].[Cl-]',           # bonded to oxygen, not to carbon: an alkoxide, not a reagent
              'CC[Zn]Cl.[Cl-]',         # the metal is satisfied; the second chloride is a counterion
              'CC[Zn+].C[Cl]',          # the halogen is bonded, so it is not a free halide
              'CC[Cu+].[Cl-]'):         # not Zn or Mg
        m = smiles(s)
        before = format(m, '')
        standardize(m)
        assert format(m, '') == before, s


def test_a_metal_with_no_halide_anywhere_is_charged_instead():
    """A neutral one-coordinate metal is not a species: the halide is missing from the DRAWING.

    Zinc and magnesium alike, and it is the one thing `standardize()` does that moves net charge.
    """
    for s, expect in [('CC[Zn]', 'CC[Zn+]'), ('CC[Mg]', 'CC[Mg+]'),
                      ('c1ccccc1[Mg]', 'c1ccccc1[Mg+]')]:
        m = smiles(s)
        assert standardize(m), s
        assert _std(s) == _std(expect), s
        assert sum(m.charge_of(n) for n in m) == 1, s
        assert not m.check_valence(), s


def test_a_metal_the_drawing_already_charged_is_not_charged_twice():
    m = smiles('CC[Zn+]')
    assert not standardize(m)
    assert sum(m.charge_of(n) for n in m) == 1


def test_a_metal_left_over_when_the_halides_run_short_is_charged():
    """One chloride between two reagents: rank decides who takes it, and the other is charged."""
    assert _one(['C[Zn+].CC[Zn].[Cl-]', '[Cl-].CC[Zn].C[Zn+]', 'CC[Zn].[Cl-].C[Zn+]'])


# the two orderings, each proved by permuting the input


def test_a_mixed_halide_set_is_decided_by_element_and_not_by_input_order():
    """`Cl > Br > I > F`: the metal takes the preferred halogen and the rest stay as counterions."""
    assert _one(['CC[Zn+].[Cl-].[Br-]', 'CC[Zn+].[Br-].[Cl-]',
                 '[Br-].[Cl-].CC[Zn+]']) == _std('CC[Zn]Cl.[Br-]')
    assert _one(['CC[Zn+].[Br-].[I-]', '[I-].CC[Zn+].[Br-]']) == _std('CC[Zn]Br.[I-]')
    assert _one(['CC[Zn+].[I-].[F-]', '[F-].[I-].CC[Zn+]']) == _std('CC[Zn]I.[F-]')


def test_interchangeable_halides_give_one_answer_whichever_is_picked():
    """Two identical chlorides are the same choice twice, so the tie needs no breaking."""
    assert _one(['CC[Zn+].[Cl-].[Cl-]', '[Cl-].CC[Zn+].[Cl-]',
                 '[Cl-].[Cl-].CC[Zn+]']) == _std('CC[Zn]Cl.[Cl-]')


def test_which_metal_takes_the_single_halide_is_decided_by_canonical_rank():
    """Two different reagents drawn with one chloride between them: rank picks, and picks the same way.

    This is the case that cannot be a rule-table row -- the table applies whichever match the
    isomorphism search returned first, which is what the input order decides.
    """
    assert _one(['C[Zn+].CC[Zn+].[Cl-]', 'CC[Zn+].C[Zn+].[Cl-]', '[Cl-].CC[Zn+].C[Zn+]',
                 'C[Zn+].[Cl-].CC[Zn+]'])


def test_equivalent_metals_are_automorphic_so_the_tie_is_not_observable():
    assert _one(['C[Zn+].C[Zn+].[Cl-]', '[Cl-].C[Zn+].C[Zn+]', 'C[Zn+].[Cl-].C[Zn+]'])


def test_several_metals_and_several_halides_are_paired_one_each():
    """Two reagents, two chlorides: each metal takes one, and no permutation changes the pairing."""
    assert _one([f'{a}.{b}.{c}.{d}' for a, b, c, d in
                 permutations(['C[Zn+]', 'CC[Zn+]', '[Cl-]', '[Cl-]'])]) == _std('C[Zn]Cl.CC[Zn]Cl')


# what the join is worth downstream


def test_the_joined_form_is_the_one_the_reaction_corpus_names():
    """The ion pair answers with the generic carbanion; joined, each reagent names itself."""
    for s, group in [('CC[Zn+].[Cl-]', 'alkyl_zinc'), ('c1ccccc1[Mg+].[Br-]', 'aryl_grignard')]:
        m = smiles(s)
        assert 'metalate_carbanion' in m.functional_groups()
        standardize(m)
        assert group in m.functional_groups(), s


def test_the_join_leaves_a_clean_valence_and_recomputes_the_hydrogen_count():
    for s in ('CC[Zn+].[Cl-]', 'CC[Zn].[Cl]', 'c1ccccc1[Mg+].[Br-]'):
        m = smiles(s)
        standardize(m)
        assert not m.check_valence(), s


# the record


def test_the_join_is_recorded_against_a_table_qualified_id():
    m = smiles('CC[Zn+].[Cl-]')
    standardize(m)
    records = [r for r in m.log if r.rule.startswith('organometallics:')]
    assert len(records) == 1, [r.rule for r in m.log]
    assert records[0].rule == 'organometallics:unite'
    assert m.log.by_stage('standardize')


def test_charging_a_stranded_metal_is_a_separate_record_from_joining_one():
    """A consumer can decline the weaker claim -- the halide nobody drew -- and keep the join."""
    m = smiles('CC[Zn]')
    standardize(m)
    records = [r for r in m.log if r.rule.startswith('organometallics:')]
    assert [r.rule for r in records] == ['organometallics:charge']
    assert 'no halide' in records[0].message

    both = smiles('C[Zn].CC[Zn].[Cl-]')
    standardize(both)
    assert {r.rule for r in both.log if r.rule.startswith('organometallics:')} == \
        {'organometallics:unite', 'organometallics:charge'}

    # and the two halves compose to one answer: whether the leftover metal arrived neutral or already
    # `+`, the halide goes to the same reagent and the other ends up charged either way.
    assert _std('C[Zn].CC[Zn].[Cl-]') == _std('C[Zn+].CC[Zn].[Cl-]') == _std('CC[Zn]Cl.C[Zn+]')


def test_a_join_that_moves_net_charge_says_so_and_one_that_does_not_stays_quiet_about_it():
    def charge(m):
        return sum(m.charge_of(n) for n in m)

    conserving, moving = smiles('CC[Zn+].[Cl-]'), smiles('CC[Zn+].[Cl]')
    assert charge(conserving) == 0 and charge(moving) == 1
    for m in (conserving, moving):
        standardize(m)
    assert charge(conserving) == 0 and charge(moving) == 0
    said = [r for r in moving.log if r.rule == 'organometallics:unite'][0].message
    assert 'charge' in said
    assert 'charge' not in [r for r in conserving.log if r.rule == 'organometallics:unite'][0].message


# the pipeline


def test_canonicalize_joins_it_too():
    m = smiles('CC[Zn+].[Cl-]')
    canonicalize(m)
    assert format(m, '') == _std('CC[Zn]Cl')
