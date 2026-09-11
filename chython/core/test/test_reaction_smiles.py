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
"""Reaction SMILES, both directions.

The invariants these tests exist to protect, in order of how expensive they are to lose:

  * a `>` that belongs to a dative `->` is NOT the reaction arrow, so a record using one is read
    rather than refused
  * a side splits on `.` into one molecule per component, and `f:` is the ONLY thing that puts them
    back together: without it a salt reactant silently becomes two reactants
  * one CXSMILES tail for the whole string, indices counting atoms across reactants, then agents,
    then products
  * the three-part form READS here and is refused by `read_smirks`, and both messages say why
  * syntax raises; chemistry never does
"""
from pytest import raises

from chython.core import MoleculeContainer, ReactionContainer
from chython.core._core import (IncorrectSmiles, read_reaction_smiles, read_smiles,
                                write_reaction_smiles)
from . import oracle


def test_a_two_part_reaction_smiles_reads_into_a_reaction_container():
    r = read_reaction_smiles('CC>>CO')
    assert isinstance(r, ReactionContainer)
    assert [format(m, 'A') for m in r.reactants] == ['CC']
    assert [format(m, 'A') for m in r.products] == ['CO']
    assert r.agents == ()


def test_a_dative_arrow_is_not_the_reaction_arrow():
    """`N->[Cu]>>N` has three `>` bytes and one arrow."""
    r = read_reaction_smiles('N->[Cu]>>N')
    assert len(r.reactants) == 1
    assert len(r.products) == 1
    reactant = r.reactants[0]
    assert [b.order for b in reactant.bonds()] == [8]


def test_a_side_splits_into_one_molecule_per_component():
    r = read_reaction_smiles('[Na+].[Cl-]>>CC')
    # the string's order, not a sorted one -- sorting is what the WRITER does to make an identifier,
    # and a reader that reordered would lose the record as it was written
    assert [format(m, 'A') for m in r.reactants] == ['[Na+]', '[Cl-]']


def test_agents_read_from_the_middle_side():
    r = read_reaction_smiles('CC>O>CO')
    assert [format(m, 'A') for m in r.agents] == ['O']


def test_the_tail_indexes_atoms_across_every_side_in_order():
    """Daylight's reaction-CXSMILES rule: reactants, then agents, then products, one index space."""
    r = read_reaction_smiles('CC>N>CO |^1:0,4|')
    assert list(r.reactants[0].atoms())[0].is_radical
    assert not list(r.agents[0].atoms())[0].is_radical
    assert list(r.products[0].atoms())[1].is_radical


def test_f_regroups_components_into_one_molecule():
    """Without `f:` a salt reactant silently becomes two reactants, and the round trip cannot survive."""
    r = read_reaction_smiles('[Na+].[Cl-]>>CC |f:0.1|')
    assert len(r.reactants) == 1
    assert r.reactants[0].connected_components_count == 2


def test_any_side_may_be_empty():
    assert read_reaction_smiles('CC>>').products == ()
    assert read_reaction_smiles('>>CC').reactants == ()
    empty = read_reaction_smiles('>>')
    assert empty.reactants == empty.agents == empty.products == ()


def test_read_smiles_returns_a_reaction_when_the_string_has_an_arrow():
    """The polymorphic door: `smiles` IS `read_smiles`, so one function has to answer both shapes."""
    assert isinstance(read_smiles('CC>>CO'), ReactionContainer)
    assert isinstance(read_smiles('CC'), MoleculeContainer)


def test_read_smiles_does_not_mistake_a_dative_bond_for_an_arrow():
    """`N->[Cu]` is a molecule, and the dispatch has to know that before it splits anything."""
    m = read_smiles('N->[Cu]')
    assert isinstance(m, MoleculeContainer)
    assert [b.order for b in m.bonds()] == [8]


def test_read_smiles_shares_its_log_with_the_reaction_reader():
    log = []
    r = read_smiles('[Na+].[Cl-]>>CC |f:0.9|', log=log)
    assert isinstance(r, ReactionContainer)
    assert any('f:0.9' in line for line in log)


def test_read_reaction_smiles_refuses_a_string_with_no_arrow():
    """The strict door's whole reason: the wrong shape fails rather than coming back as a molecule."""
    with raises(IncorrectSmiles) as e:
        read_reaction_smiles('CC')
    assert 'two `>`' in str(e.value)


def test_read_reaction_smiles_refuses_one_arrow_and_three_arrows():
    for text in ('CC>CO', 'CC>>CO>N'):
        with raises(IncorrectSmiles) as e:
            read_reaction_smiles(text)
        assert 'reactants>agents>products' in str(e.value), text


def test_a_reaction_writes_as_three_sides_separated_by_arrows():
    r = read_reaction_smiles('CC>N>CO')
    assert write_reaction_smiles(r) == 'CC>N>CO'


def test_str_and_format_and_the_property_are_the_one_writer():
    r = read_reaction_smiles('CC>>CO')
    assert str(r) == 'CC>>CO'
    assert r.smiles == 'CC>>CO'
    assert format(r, '') == 'CC>>CO'


def test_each_side_is_sorted_so_the_string_is_an_identifier():
    """The same reaction built in either order gives one string, as `write_smiles` promises per molecule."""
    one = read_reaction_smiles('CC.O>>CO')
    other = read_reaction_smiles('O.CC>>CO')
    assert str(one) == str(other)


def test_c_keeps_the_records_own_order():
    """V2's key, reused: `!c` is the documented off-switch for the sort above."""
    assert format(read_reaction_smiles('O.CC>>CO'), '!c') == 'O.CC>>CO'
    assert format(read_reaction_smiles('CC.O>>CO'), '!c') == 'CC.O>>CO'


def test_a_multi_component_molecule_writes_an_f_group():
    """Without `f:` the salt comes back as two reactants; this is the round trip's load-bearing field."""
    r = read_reaction_smiles('[Na+].[Cl-]>>CC |f:0.1|')
    written = str(r)
    assert '|f:0.1|' in written
    assert len(read_reaction_smiles(written).reactants) == 1


def test_the_tail_aggregates_radicals_across_every_side():
    """ONE field for the whole string, not one per side -- a reactant's radical and a product's share
    an index space.

    The indices are the WRITER's atom order and need not be the reader's: `CC>N>CO |^1:0,4|` comes
    back as `C[CH2]>N>C[O] |^1:1,4|`, because each molecule is written in its own canonical order.
    So what is asserted here is aggregation, not the literal string: one `^1:` field, two indices in
    it, and the radicals landing on the same sides when it is read back.
    """
    written = str(read_reaction_smiles('CC>N>CO |^1:0,4|'))
    assert written.count('|') == 2 and written.count('^1:') == 1
    back = read_reaction_smiles(written)
    assert sum(a.is_radical for a in back.reactants[0].atoms()) == 1
    assert sum(a.is_radical for a in back.agents[0].atoms()) == 0
    assert sum(a.is_radical for a in back.products[0].atoms()) == 1


def test_the_labels_span_every_side_in_one_field():
    """`$...$` counts atoms across the whole string exactly as `^1:` does, and it is POSITIONAL -- so
    a side with no label still holds its own atoms' worth of blanks."""
    rxn = read_reaction_smiles('CCO>>CC |$Me;;OH;;$|')
    assert [m.aliases for m in rxn.molecules()][1] == {}
    assert sorted(rxn.reactants[0].aliases.values()) == [b'Me', b'OH']
    written = str(rxn)
    assert written.count('$') == 2 and written.split('|')[1].count(';') == 4
    back = read_reaction_smiles(written)
    assert sorted(back.reactants[0].aliases.values()) == [b'Me', b'OH']
    assert not back.products[0].aliases


def test_every_documented_spec_key_is_accepted_and_nothing_else_is():
    """The key set is V2's, and this is the test that keeps `docs/reactions.rst` honest about it.

    A key named in `docs/reactions.rst` prose is executed by nothing, and a documented key that raises
    looks exactly like a documented key that works until somebody calls it -- `format(rxn, '!C')`, for
    instance, which this writer refuses.  So the set is pinned here, where the suite runs it.

    `r` reaches every molecule's writer like any other key, so a reaction is augmented the same way one
    molecule is; `ir` is in the second list because the two atom-order keys are mutually exclusive.
    """
    rxn = read_reaction_smiles('[CH3:1][OH:2].[Na+:4].[Cl-:5]>[H+:9]>[CH3:1][NH2:3]')
    for key in ('', '!c', 'a', '!s', 'A', 'm', 'h', '!b', '!x', '!z', 'r'):
        assert format(rxn, key), key

    for key in ('!C', 'C', 'Q', 'ir'):
        with raises(ValueError):
            format(rxn, key)


def test_the_tail_is_suppressed_by_the_writers_own_key():
    r = read_reaction_smiles('CC>>CO |^1:0|')
    assert '|' not in format(r, '!x')


def test_a_reaction_round_trips_under_a_random_atom_order():
    """`r` is forwarded to each molecule, and the CXSMILES tail indexes atoms by their position across
    the whole reaction -- so an order the tail did not follow would put a radical or an `f:` component
    group on another atom, and the string would read back as a different reaction."""
    for text in ('[Na+].[Cl-]>>CC |f:0.1|', 'CC>>CO |^1:0|', 'C[C@H](N)O>>C[C@@H](N)O',
                 'C/C=C/C>>C/C=C\\C', '[CH3:1][CH3:2]>>[CH3:1][OH:2]'):
        r = read_reaction_smiles(text)
        for _ in range(20):
            assert read_reaction_smiles(format(r, 'rm')) == r, text


def test_a_reaction_round_trips_through_its_own_string():
    for text in ('CC>>CO', 'CC>N>CO', 'CC>>', '>>CC',
                 '[Na+].[Cl-]>>CC |f:0.1|', 'C[C@H](N)O>>C[C@@H](N)O',
                 'C/C=C/C>>C/C=C\\C', 'N->[Cu]>>N', '[CH3:1][CH3:2]>>[CH3:1][OH:2]'):
        r = read_reaction_smiles(text)
        assert read_reaction_smiles(str(r)) == r, text


def test_enhanced_stereo_groups_survive_a_reaction_round_trip():
    """The case chython 2 CANNOT pass: its reaction tail carried `^1:` and `f:` and nothing else, so
    every reaction SMILES it wrote silently dropped `&n:` / `on:` / `a:`.

    Asserted on the groups' shape rather than on their ids: an id is opaque, and the writer renumbers
    them per side exactly as `DetachedSmiles.join` does, so `&1` on the product side of the input is
    not promised to still be spelled `&1`.
    """
    r = read_reaction_smiles('C[C@H](N)O.C[C@@H](N)O>>CC=O |&1:1,&2:6|')
    written = str(r)
    assert '&1:' in written and '&2:' in written
    back = read_reaction_smiles(written)
    for before, after in zip(r.reactants, back.reactants):
        assert (sorted(len(v) for v in before.stereo_groups().values())
                == sorted(len(v) for v in after.stereo_groups().values()))
    assert sum(m.has_stereo_groups for m in back.reactants) == 2


def test_a_carried_stereo_group_is_one_group_across_the_arrow():
    # The mapping says the two `&1`s are one centre restated, so the tail says so too.  Before this it
    # came out `|&1:0,&2:4|` -- two independently racemic centres where the record states one.
    r = read_reaction_smiles('[CH3:1][C@H:2]([NH2:3])[OH:4]'
                             '>>[CH3:1][C@H:2]([NH2:3])[Cl:5] |&1:1,&2:5|')
    assert str(r) == '[C@@H](C)(N)O>>[C@@H](C)(N)Cl |&1:0,4|'


def test_two_unmapped_stereo_groups_stay_two_groups():
    # The negative control: nothing links them, so nothing may merge them.
    r = read_reaction_smiles('C[C@H](N)O.C[C@@H](N)O>>CC=O |&1:1,&2:6|')
    assert str(r) == '[C@@H](C)(N)O.[C@H](C)(N)O>>C(C)=O |&1:0,&2:6|'


def test_an_and_group_and_an_or_group_never_merge():
    # Two kinds, two fields; a shared map number cannot make an AND set into an OR set.
    r = read_reaction_smiles('[CH3:1][C@H:2]([NH2:3])[OH:4]'
                             '>>[CH3:1][C@H:2]([NH2:3])[Cl:5] |&1:1,o1:5|')
    assert str(r) == '[C@@H](C)(N)O>>[C@@H](C)(N)Cl |&1:0,o1:4|'


def test_read_smirks_still_refuses_the_three_part_form_the_reader_accepts():
    """The two refusals sit next to each other on purpose, and each names the other reading: a
    reaction RECORD has agents, and a TEMPLATE cannot patch one."""
    from chython.core._core import IncorrectSmirks, read_smirks

    assert len(read_reaction_smiles('CC>N>CO').agents) == 1
    with raises(IncorrectSmirks) as e:
        read_smirks('[C:1]C>N>[C:1]O')
    assert 'agent' in str(e.value)


# --------------------------------------------------------------------------------------------------
# DIFFERENTIAL against chython 2.24, run out of process.  Both directions in one comparison: V2 reads
# the original string, V2 reads OUR rewrite of it, and the two readings must agree.  That is a stronger
# claim than comparing our string to V2's -- it says our writer emits something V2 accepts AND that
# means the same reaction -- and it never compares one library's canonical text to the other's, which
# this repository has ruled unsound.
#
# NO STEREOCENTRES IN THIS CORPUS.  V2's canonical writer is order-dependent on symmetric centres, so
# a disagreement there would be V2's and not ours; the round-trip test above is where stereo is stated.

ORACLE_SPLIT = r'''
from chython import smiles

out = []
for s in _payload:
    try:
        r = smiles(s)
    except Exception as e:
        out.append({'error': '%s: %s' % (type(e).__name__, e)})
        continue
    out.append({'reactants': sorted(str(m) for m in r.reactants),
                'agents': sorted(str(m) for m in r.reagents),
                'products': sorted(str(m) for m in r.products)})
_emit(out)
'''

DIFFERENTIAL = ['CCO.CC(=O)O>[H+]>CCOC(C)=O.O',
                '[CH3:1][C:2](=[O:3])[OH:4].[NH2:5][CH3:6]>>'
                '[CH3:1][C:2](=[O:3])[NH:5][CH3:6].[OH2:4]',
                'Brc1ccccc1.OB(O)c1ccccc1>>c1ccc(-c2ccccc2)cc1',
                'CCO>>CC=O',
                'CC(=O)Cl.CCN>CCN(CC)CC>CC(=O)NCC']


@oracle.requires_oracle
def test_what_we_write_reads_back_in_chython_two_as_the_same_reaction():
    ours = [str(read_reaction_smiles(text)) for text in DIFFERENTIAL]
    answers = oracle.ask(ORACLE_SPLIT, DIFFERENTIAL + ours)
    for i, text in enumerate(DIFFERENTIAL):
        original, rewritten = answers[i], answers[i + len(DIFFERENTIAL)]
        assert 'error' not in rewritten, (ours[i], rewritten.get('error'))
        assert 'error' not in original, (text, original.get('error'))
        assert original == rewritten, (text, ours[i])


@oracle.requires_oracle
def test_chython_two_refuses_what_this_reader_accepts_and_that_is_the_improvement():
    """The two records V2 cannot read at all, kept as a differential in the negative direction: this
    is the compatibility gap the reader closes, and a future V2 that could read them would mean the
    corpus assumption changed."""
    answers = oracle.ask(ORACLE_SPLIT, ['N->[Cu]>>N', 'CC>>CO |f:0.1|'])
    assert 'error' in answers[0], answers[0]
    assert len(read_reaction_smiles('N->[Cu]>>N').reactants) == 1
