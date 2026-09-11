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
from pytest import raises

from ...core import ReactionContainer, read_smiles as smiles


def test_the_accessor_names_the_package_when_unregistered():
    from ...core._core import _reaction_reconstruct_fn, _set_reconstruct_fn
    try:
        kept = _reaction_reconstruct_fn()       # another test module may have registered it already
    except ImportError:
        kept = None
    _set_reconstruct_fn(None)
    try:
        with raises(ImportError, match='chython.reactions'):
            _reaction_reconstruct_fn()
    finally:
        _set_reconstruct_fn(kept)


def test_purification_is_the_identity_mapping():
    rxn = ReactionContainer([smiles('CCO'), smiles('O')], [smiles('CCO')])
    assert rxn.reconstruct_mapping() == ('purification',)
    product = rxn.products[0]
    numbers = {product.map_number_of(n) for n in product.atom_numbers}
    assert 0 not in numbers                                     # every product atom is numbered
    ethanol = rxn.reactants[0]
    assert numbers == {ethanol.map_number_of(n) for n in ethanol.atom_numbers}

    # PER ATOM, and not only as a set: the set assertion passes for any consistent permutation inside
    # ethanol.  `(element, degree, implicit_h)` separates all three of its atoms.
    source = {ethanol.map_number_of(n): n for n in ethanol.atom_numbers}
    for n in product.atom_numbers:
        m = source[product.map_number_of(n)]
        assert (product.atom(n).element, product.degree_of(n), product.atom(n).implicit_h) == \
               (ethanol.atom(m).element, ethanol.degree_of(m), ethanol.atom(m).implicit_h)


def test_a_multiproduct_record_is_refused():
    from ...core import Log                     # severities are `INFO`/`REFUSED` constants, not an enum
    rxn = ReactionContainer([smiles('CCOC(C)=O'), smiles('O')], [smiles('CCO'), smiles('CC(=O)O')])
    assert rxn.reconstruct_mapping() == ()
    log = rxn.log
    refused = log.refused()
    assert len(refused) == 1
    assert refused[0].rule == 'reconstruct:multiproduct'
    assert refused[0].stage == 'reconstruct'        # the pass names the stage; nothing above it does


def test_an_empty_side_is_refused():
    from ...core import Log
    rxn = ReactionContainer([smiles('CCO')], [])
    assert rxn.reconstruct_mapping() == ()
    log = rxn.log
    assert [r.rule for r in log.refused()] == ['reconstruct:empty']


def test_an_empty_reactant_side_is_refused():
    from ...core import Log
    rxn = ReactionContainer([], [smiles('CCO')])
    assert rxn.reconstruct_mapping() == ()
    log = rxn.log
    assert [r.rule for r in log.refused()] == ['reconstruct:empty']


def test_incoming_map_numbers_are_replaced_unconditionally():
    # 10/20/30 rather than 1/2/3 keeps the test non-vacuous: `_number_inputs` assigns 1..N, so the
    # incoming numbers cannot be mistaken for the derived ones even if `_clear` is removed.
    ethanol_in = smiles('CCO')
    for atom_n, mn in zip(list(ethanol_in.atom_numbers), [10, 20, 30]):
        ethanol_in.set_map_number(atom_n, mn)
    ethanol_out = smiles('CCO')
    for atom_n, mn in zip(list(ethanol_out.atom_numbers), [10, 20, 30]):
        ethanol_out.set_map_number(atom_n, mn)
    rxn = ReactionContainer([ethanol_in], [ethanol_out])
    assert rxn.reconstruct_mapping() == ('purification',)
    product = rxn.products[0]
    product_numbers = {product.map_number_of(n) for n in product.atom_numbers}
    assert 10 not in product_numbers
    assert 20 not in product_numbers
    assert 30 not in product_numbers
    assert 0 not in product_numbers   # every product atom has been numbered
    reactant = rxn.reactants[0]
    reactant_numbers = {reactant.map_number_of(n) for n in reactant.atom_numbers}
    assert product_numbers == reactant_numbers


def test_unexplained_record_emits_lost_and_clears_product():
    # Ethanol -> propane is not a reaction: an underivable mapping comes back with a LOST line rather
    # than silently disappearing, and the product's numbers are cleared and not replaced.
    from ...core import Log, LOST
    rxn = ReactionContainer([smiles('[CH3:1][CH2:2][OH:3]')], [smiles('[CH3:4][CH2:5][CH3:6]')])
    result = rxn.reconstruct_mapping()
    log = rxn.log
    assert result == ()
    lost = [r for r in log if r.severity == LOST]
    assert len(lost) == 1
    assert lost[0].rule == 'reconstruct:unexplained'
    product = rxn.products[0]
    product_numbers = {product.map_number_of(n) for n in product.atom_numbers}
    assert product_numbers == {0}  # all cleared, none replaced


def test_the_recorded_reaction_comes_back_canonical():
    # a Kekule aromatic input is rewritten in place: the answer is a normalized, mapped record
    rxn = ReactionContainer([smiles('C1=CC=CC=C1O')], [smiles('C1=CC=CC=C1O')])
    assert rxn.reconstruct_mapping() == ('purification',)
    assert all(rxn.reactants[0].atom(n).hybridization == 4
               for n in rxn.reactants[0].atom_numbers if rxn.reactants[0].atom(n).element == 6)


def test_translate_carries_the_inputs_numbering_onto_a_reactor_product():
    # The reactor numbers its own output; every rung of the ladder reads the INPUTS' numbering.
    # Asserted as set membership, not as a literal, so the test survives whichever scheme the reactor
    # imposes: the product's carbonyl carbon has to come back with a number the ACID input carries.
    from .._enumerate import _run
    from .._reconstruct import _number_inputs, _translate
    from .._tables import reaction_rules

    acid, amine = smiles('CC(=O)O'), smiles('CCN')
    acid.canonicalize()
    amine.canonicalize()
    _number_inputs([acid, amine])
    outcome = next(o for o in _run([acid, amine], reaction_rules()) if o.name == 'amidation')
    product = _translate(outcome.reaction, [acid, amine])[0]

    acid_numbers = {acid.map_number_of(n) for n in acid.atom_numbers}
    amine_numbers = {amine.map_number_of(n) for n in amine.atom_numbers}
    carbonyl = next(n for n in product.atom_numbers
                    if product.element_of(n) == 6
                    and any(product.order_of(n, m) == 2 and product.element_of(m) == 8
                            for m in product.neighbors_of(n)))
    nitrogen = next(n for n in product.atom_numbers if product.element_of(n) == 7)
    assert product.map_number_of(carbonyl) in acid_numbers
    assert product.map_number_of(nitrogen) in amine_numbers


def test_a_direct_reaction_is_found_and_mapped():
    rxn = ReactionContainer([smiles('CC(=O)O'), smiles('CCN')], [smiles('CC(=O)NCC')])
    assert rxn.reconstruct_mapping() == ('react:amidation',)
    product = rxn.products[0]
    assert all(product.map_number_of(n) for n in product.atom_numbers
               if product.atom(n).element != 1)
    acid, amine = rxn.reactants
    acid_numbers = {acid.map_number_of(n) for n in acid.atom_numbers}
    amine_numbers = {amine.map_number_of(n) for n in amine.atom_numbers}
    assert acid_numbers.isdisjoint(amine_numbers)
    carried = {product.map_number_of(n) for n in product.atom_numbers}
    assert carried & acid_numbers and carried & amine_numbers
    # the carbonyl carbon carries the acid's number and not the amine's -- the literal claim the
    # assertion above only implies
    carbonyl_c = next(n for n in product.atom_numbers
                      if product.atom(n).element == 6
                      and any(product.bond(n, nb).order == 2 and product.atom(nb).element == 8
                              for nb in product.neighbors_of(n)))
    assert product.map_number_of(carbonyl_c) in acid_numbers
    assert product.map_number_of(carbonyl_c) not in amine_numbers


def test_a_spectator_input_does_not_block_the_match():
    # toluene is a solvent no row touches.  This is the subset enumeration: without it `_run` refuses
    # the whole list, every input having to be touched, and nothing is found.
    rxn = ReactionContainer([smiles('CC(=O)O'), smiles('CCN'), smiles('Cc1ccccc1')],
                            [smiles('CC(=O)NCC')])
    assert rxn.reconstruct_mapping() == ('react:amidation',)


def test_an_untouched_salt_component_survives_and_is_numbered():
    # the recorded product carries a chloride the template never sees; it must come back numbered
    # from the input it arrived in, not dropped and not left at zero
    rxn = ReactionContainer([smiles('CC(=O)O'), smiles('CCN.Cl')], [smiles('CC(=O)NCC.Cl')])
    assert rxn.reconstruct_mapping() == ('react:amidation',)
    amine_input = rxn.reactants[1]
    product = rxn.products[0]
    chloride = [n for n in product.atom_numbers if product.atom(n).element == 17]
    assert len(chloride) == 1
    input_cl = next(n for n in amine_input.atom_numbers if amine_input.atom(n).element == 17)
    assert product.map_number_of(chloride[0]) == amine_input.map_number_of(input_cl)


def test_the_winning_rule_id_is_logged():
    from ...core import Log
    rxn = ReactionContainer([smiles('CC(=O)O'), smiles('CCN')], [smiles('CC(=O)NCC')])
    assert rxn.reconstruct_mapping()
    log = rxn.log
    ids = [r.message for r in log if r.rule == 'reconstruct:rules']
    # pin the value and not the count, which is trivially 1 when only one row matches
    assert ids == ['explained by reactions:1']


def test_a_partial_product_emits_the_partial_log_line():
    from ...core import Log, INFO
    # the water in the recorded product came from nowhere, so `_number_product` cannot account for that
    # component and the orchestrator emits `reconstruct:partial`
    rxn = ReactionContainer([smiles('CC(=O)O'), smiles('CCN')], [smiles('CC(=O)NCC.O')])
    result = rxn.reconstruct_mapping()
    log = rxn.log
    assert result == ('react:amidation',)                     # rung fires despite the unaccounted component
    partial = [r for r in log if r.rule == 'reconstruct:partial']
    assert len(partial) == 1 and partial[0].severity == INFO
    product = rxn.products[0]
    components = list(product.split())
    water = next(c for c in components if len(list(c.atom_numbers)) == 1)
    amide = next(c for c in components if len(list(c.atom_numbers)) > 1)
    assert all(product.map_number_of(n) == 0 for n in water.atom_numbers)
    assert all(product.map_number_of(n) != 0 for n in amide.atom_numbers
               if product.atom(n).element != 1)


def test_a_standalone_deprotection_is_explained():
    # Boc removal: the recorded product is the recorded input, minus the protecting group
    rxn = ReactionContainer([smiles('CC(C)(C)OC(=O)NCc1ccccc1'), smiles('Cl')],
                            [smiles('NCc1ccccc1')])
    labels = rxn.reconstruct_mapping()
    assert labels and all(label.startswith('deprotect:') for label in labels)
    product = rxn.products[0]
    assert all(product.map_number_of(n) for n in product.atom_numbers)
    revealed = rxn.reactants[0]
    numbers = {revealed.map_number_of(n) for n in revealed.atom_numbers}
    assert {product.map_number_of(n) for n in product.atom_numbers} <= numbers

    # PER ATOM, and not only as a set: the set assertion passes for a swap of two ring CH atoms or of
    # the benzylic CH2 for a ring CH.  The revealed N is the one exception -- Boc removal changes its
    # degree and implicit_h -- so it is checked on element alone.
    source = {revealed.map_number_of(n): n for n in revealed.atom_numbers}
    for n in product.atom_numbers:
        m = source[product.map_number_of(n)]
        if product.atom(n).element == 7:      # the revealed atom; update if the substrate changes
            assert product.atom(n).element == revealed.atom(m).element
        else:
            assert (product.atom(n).element, product.degree_of(n), product.atom(n).implicit_h) == \
                   (revealed.atom(m).element, revealed.degree_of(m), revealed.atom(m).implicit_h)


def test_a_partial_deprotection_is_explained():
    # one of two Boc groups comes off: reachable only because the rung enumerates partial strips
    rxn = ReactionContainer([smiles('CC(C)(C)OC(=O)N(C(=O)OC(C)(C)C)Cc1ccccc1')],
                            [smiles('CC(C)(C)OC(=O)NCc1ccccc1')])
    labels = rxn.reconstruct_mapping()
    assert labels and all(label.startswith('deprotect:') for label in labels)


def test_deprotect_then_react_composes():
    # the amine arrives Boc-protected; strip it, then amidate.  Neither rung alone explains this.
    rxn = ReactionContainer([smiles('CC(=O)O'), smiles('CC(C)(C)OC(=O)NCC')],
                            [smiles('CC(=O)NCC')])
    labels = rxn.reconstruct_mapping()
    assert labels and all(label.startswith('deprotect+react:') for label in labels)


def test_a_direct_reaction_outranks_the_composed_rung():
    # NOT an ordering pin: CCN carries no protecting group, so `_deprotect_then_react` returns at once
    # whatever its position in `_PHASES`.  `test_phases_ladder_order_is_pinned` is the ordering pin.
    rxn = ReactionContainer([smiles('CC(=O)O'), smiles('CCN')], [smiles('CC(=O)NCC')])
    assert rxn.reconstruct_mapping() == ('react:amidation',)


def test_phases_ladder_order_is_pinned():
    # The ladder is ordered by STRENGTH OF EVIDENCE, and protection is last because an amide, an ester
    # and a carbamate are products as well as protecting groups.  Behavioural counterpart:
    # `test_an_acylation_is_a_reaction_and_not_a_protection`.
    from .._reconstruct import _PHASES
    assert [f.__name__ for f in _PHASES] == ['_purification', '_react', '_deprotect',
                                             '_deprotect_then_react', '_protect']


def test_two_surviving_candidates_yield_one_label():
    # Two inputs both deprotect to benzylamine, so `_deprotect` yields two candidates that reproduce the
    # recorded product.  Exactly one is applied: applying both writes two conflicting numberings.
    from ...core import Log
    rxn = ReactionContainer(
        [smiles('CC(C)(C)OC(=O)NCc1ccccc1'), smiles('O=C(OCc1ccccc1)NCc1ccccc1')],
        [smiles('NCc1ccccc1')],
    )
    labels = rxn.reconstruct_mapping()
    log = rxn.log
    assert len(labels) == 1 and labels[0].startswith('deprotect:')
    rules_logged = [r for r in log if r.rule == 'reconstruct:rules']
    assert len(rules_logged) == 1
    # pin the content and not only the count: the message must name a row and not a label
    assert rules_logged[0].message.startswith('explained by protective:')


def test_a_protection_is_explained():
    # Boc protection of benzylamine with Boc anhydride.  The rung deprotects the recorded product,
    # pairs the revealed fragment against the input amine, and numbers back.
    rxn = ReactionContainer([smiles('NCc1ccccc1'), smiles('CC(C)(C)OC(=O)OC(=O)OC(C)(C)C')],
                            [smiles('CC(C)(C)OC(=O)NCc1ccccc1')])
    labels = rxn.reconstruct_mapping()
    assert labels and all(label.startswith('protect:') for label in labels)
    product = rxn.products[0]
    amine = rxn.reactants[0]
    numbers = {amine.map_number_of(n) for n in amine.atom_numbers}
    carried = {product.map_number_of(n) for n in product.atom_numbers if product.map_number_of(n)}
    assert carried == numbers                       # the whole revealed fragment is numbered
    # the protecting group's own atoms are NEW and correctly carry no number
    assert any(product.map_number_of(n) == 0 for n in product.atom_numbers)

    # PER ATOM, and not only as a set.  `(element, degree, implicit_h)` separates the distinguishable
    # atoms; the ring CH positions share a triple and cannot be separated, being one symmetry orbit.
    source = {amine.map_number_of(n): n for n in amine.atom_numbers}
    for n in product.atom_numbers:
        mn = product.map_number_of(n)
        if mn == 0:        # protecting group atom, genuinely new -- no source to check against
            continue
        m = source[mn]
        if product.atom(n).element == 7:      # the reaction site; update if the substrate changes
            assert product.atom(n).element == amine.atom(m).element
        else:
            assert (product.atom(n).element, product.degree_of(n), product.atom(n).implicit_h) == \
                   (amine.atom(m).element, amine.degree_of(m), amine.atom(m).implicit_h)


def test_an_acylation_is_a_reaction_and_not_a_protection():
    # Why protection is the LAST rung: offered first, this reads as `protect:amine_benzoate`.
    rxn = ReactionContainer([smiles('OC(=O)c1ccccc1'), smiles('CCN')],
                            [smiles('CCNC(=O)c1ccccc1')])
    labels = rxn.reconstruct_mapping()
    assert labels
    assert not any(label.startswith('protect:') for label in labels)


def test_canonicalize_does_not_reconstruct():
    # NEVER IMPLICITLY COMPOSED: a caller asking for a canonical representation has not asked for a
    # mapping to be invented.  This would fail if `canonicalize()` called `reconstruct_mapping`, the
    # amidation row firing and writing nonzero numbers onto the product.
    rxn = ReactionContainer([smiles('CC(=O)O'), smiles('CCN')], [smiles('CC(=O)NCC')])
    rxn.canonicalize()
    product = rxn.products[0]
    assert not any(product.map_number_of(n) for n in product.atom_numbers)


def test_a_grossly_larger_product_is_refused():
    from ...core import Log
    # a 60-atom product from one 2-atom input: the corpus has nothing honest to say about this
    rxn = ReactionContainer([smiles('CO')], [smiles('C' * 60)])
    assert rxn.reconstruct_mapping() == ()
    log = rxn.log
    assert [r.rule for r in log.refused()] == ['reconstruct:unbalanced']
    # a bound on the SEARCH and not a rejection: the record still gets its ordinary unexplained line
    from ...core import LOST
    assert any(r.rule == 'reconstruct:unexplained' and r.severity == LOST for r in log)


def test_the_filter_is_off_below_the_floor():
    from ...core import Log
    # under `min_filter_size` the ratio is not consulted at all, however lopsided it looks
    rxn = ReactionContainer([smiles('CO')], [smiles('CCCCCCCCCC')])
    rxn.reconstruct_mapping()
    log = rxn.log
    assert not [r for r in log.refused() if r.rule == 'reconstruct:unbalanced']


def test_the_filter_can_be_disabled():
    from ...core import Log
    rxn = ReactionContainer([smiles('CO')], [smiles('C' * 60)])
    rxn.reconstruct_mapping(max_size_ratio=0.)
    log = rxn.log
    assert not [r for r in log.refused() if r.rule == 'reconstruct:unbalanced']


def test_a_protection_survives_the_filter():
    # Trityl protection of decanol.  THE THRESHOLDS ARE PASSED EXPLICITLY: at the defaults the filter
    # never engages and the test would pass with `_FILTER_EXEMPT` deleted.  Forced on, only the
    # exemption can let a `protect:` label out.
    rxn = ReactionContainer([smiles('OCCCCCCCCCC')],
                            [smiles('C(c1ccccc1)(c1ccccc1)(c1ccccc1)OCCCCCCCCCC')])
    assert len(rxn.products[0].atom_numbers) == 30       # 30 >= 10, and 30 >= 1.5 * 11
    labels = rxn.reconstruct_mapping(max_size_ratio=1.5, min_filter_size=10)
    assert labels and all(label.startswith('protect:') for label in labels)


def test_a_purification_survives_the_filter():
    # A purification's product IS one of its inputs, but the arithmetic does not know that: at a low
    # enough ratio it clears the bound.  Dropping `_purification` from `_FILTER_EXEMPT` turns this
    # answer into `()` with a `reconstruct:unexplained` line.
    from ...core import Log
    rxn = ReactionContainer([smiles('CCO'), smiles('O')], [smiles('CCO')])
    assert rxn.reconstruct_mapping(max_size_ratio=.5, min_filter_size=2) == ('purification',)
    log = rxn.log
    # confirm the filter genuinely engaged, or the test passes with `_purification` unreachable
    assert [r.rule for r in log.refused()] == ['reconstruct:unbalanced']


def test_a_reference_mapping_is_reproduced_exactly():
    from .._numbering import mapping_agrees
    for reference in _reference_records():
        reference.canonicalize()
        probe = reference.copy()
        assert probe.reconstruct_mapping(), 'nothing explained a record the corpus should explain'
        agreed, disagreed, missing = mapping_agrees(probe, reference)
        assert agreed, 'no product atom was traced back to an input'
        assert (disagreed, missing) == (0, 0)


def _reference_records():
    """Three public reactions, written with the mapping a chemist would draw.

    Hand-written rather than lifted from `mapping/golden.rdf`, a test not depending on a data file outside
    the installed package.  One per rung that can reproduce a whole product.

    THE COUPLING IS DELIBERATELY ASYMMETRIC -- 4-bromotoluene, not bromobenzene.  Biphenyl's two rings are
    one automorphism orbit of the PRODUCT while arriving from two different INPUTS, and `mapping_agrees`
    excuses a swap only within one input's orbits; for the symmetric spelling no mapping is the answer, so
    a reference naming one would state a convention rather than a fact.
    """
    return [
        ReactionContainer([smiles('[CH3:1][C:2](=[O:3])[OH:4]'), smiles('[CH3:5][CH2:6][NH2:7]')],
                          [smiles('[CH3:1][C:2](=[O:3])[NH:7][CH2:6][CH3:5]')]),
        ReactionContainer([smiles('[CH3:1][C:2]([CH3:3])([CH3:4])[O:5][C:6](=[O:7])'
                                  '[NH:8][CH2:9][c:10]1[cH:11][cH:12][cH:13][cH:14][cH:15]1')],
                          [smiles('[NH2:8][CH2:9][c:10]1[cH:11][cH:12][cH:13][cH:14][cH:15]1')]),
        ReactionContainer([smiles('[Br:1][c:2]1[cH:3][cH:4][c:5]([CH3:6])[cH:7][cH:8]1'),
                           smiles('[OH:9][B:10]([OH:11])[c:12]1[cH:13][cH:14][cH:15][cH:16][cH:17]1')],
                          [smiles('[CH3:6][c:5]1[cH:4][cH:3][c:2]([c:12]2[cH:13][cH:14][cH:15][cH:16]'
                                  '[cH:17]2)[cH:8][cH:7]1')]),
    ]


def test_the_reconstructed_mapping_is_1_1_and_starts_at_one():
    rxn = ReactionContainer([smiles('CC(=O)O'), smiles('CCN')], [smiles('CC(=O)NCC')])
    assert rxn.reconstruct_mapping() == ('react:amidation',)
    product = rxn.products[0]
    numbers = sorted(product.map_number_of(n) for n in product.atom_numbers)
    assert numbers == list(range(1, len(numbers) + 1))
    left = [m.map_number_of(n) for m in rxn.reactants for n in m.atom_numbers]
    assert sorted(x for x in left if x) == numbers        # 1-1, and the same set on both sides
    assert left.count(0) == 1                             # the acid's leaving OH, and only it
    assert rxn.modeling_view().collisions == {'reactants': (), 'products': ()}


def test_a_spectator_input_comes_back_unmapped():
    # A number on an atom the product never received is not part of a 1-1 mapping, and leaving it there
    # invites a reader to treat the toluene as a reagent that contributed atoms.
    rxn = ReactionContainer([smiles('CC(=O)O'), smiles('CCN'), smiles('Cc1ccccc1')],
                            [smiles('CC(=O)NCC')])
    assert rxn.reconstruct_mapping() == ('react:amidation',)
    toluene = rxn.reactants[2]
    assert all(toluene.map_number_of(n) == 0 for n in toluene.atom_numbers)
