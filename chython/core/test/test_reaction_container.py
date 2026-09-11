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
"""The chython 3 reaction container: what it holds, what it refuses to refuse, and what it does not carry.

Public reactions only -- esterification, amide formation, a Suzuki coupling, all textbook.

The chython 2 differential at the foot runs in a SEPARATE INTERPRETER with `-I`.  Not decoration:
without `-I` the child puts this working tree on `sys.path` first and imports the very code under
test, so the comparison passes by comparing the tree to itself.
"""
from warnings import catch_warnings, simplefilter

from pytest import mark, raises, warns

from chython.core import H_UNKNOWN, read_reaction_smiles, read_smiles
from chython.core.reaction import ReactionContainer, ReactionModelingView
from . import oracle


def _esterification():
    """Ethanol + acetic acid, acid-catalysed, mapped.  Water leaves; the ester oxygen is the alcohol's."""
    ethanol = read_smiles('[CH3:1][CH2:2][OH:3]')
    acid = read_smiles('[CH3:4][C:5](=[O:6])[OH:7]')
    ester = read_smiles('[CH3:1][CH2:2][O:3][C:5](=[O:6])[CH3:4]')
    water = read_smiles('[OH2:7]')
    return ReactionContainer([ethanol, acid], [ester, water], [read_smiles('[H+]')])


def _suzuki():
    """Bromobenzene + phenylboronic acid -> biphenyl.  Textbook, public, and mapped.

    The biaryl bond is written `-` ON PURPOSE: between two aromatic atoms an ABSENT bond in SMILES
    reads as aromatic, so omitting it would make the new bond order 4 and quietly turn this into a
    test of the reader rather than of the view.
    """
    ar_br = read_smiles('[cH:1]1[cH:2][cH:3][cH:4][cH:5][c:6]1[Br:7]')
    ar_b = read_smiles('[cH:11]1[cH:12][cH:13][cH:14][cH:15][c:16]1[B:17]([OH:18])[OH:19]')
    biaryl = read_smiles('[cH:1]1[cH:2][cH:3][cH:4][cH:5][c:6]1-[c:16]1[cH:15][cH:14][cH:13][cH:12][cH:11]1')
    return ReactionContainer([ar_br, ar_b], [biaryl])


def _quiet(obj, name):
    with catch_warnings():
        simplefilter('ignore', DeprecationWarning)
        return getattr(obj, name)


# --- the three sides -----------------------------------------------------------------------------

def test_the_three_sides_are_the_tuples_they_were_given():
    rxn = _esterification()
    assert len(rxn.reactants) == 2 and len(rxn.products) == 2 and len(rxn.agents) == 1
    assert len(rxn) == 5


def test_molecules_walks_reactants_then_agents_then_products():
    rxn = _esterification()
    assert list(rxn.molecules()) == [*rxn.reactants, *rxn.agents, *rxn.products]


def test_reagents_is_the_same_tuple_agents_is():
    rxn = _esterification()
    assert _quiet(rxn, 'reagents') is rxn.agents, 'a pure rename returns the object, not a copy'


def test_the_reagents_keyword_still_constructs_and_says_what_to_spell():
    a = read_smiles('[H+]')
    with warns(DeprecationWarning, match='agents'):
        rxn = ReactionContainer([read_smiles('CCO')], [read_smiles('CC=O')], reagents=[a])
    assert rxn.agents == (a,)


def test_giving_both_spellings_is_an_error_rather_than_a_silent_winner():
    """A caller who passes both has two different intentions in one call and no way to know which
    one took effect."""
    with raises(TypeError, match='not both'):
        ReactionContainer([read_smiles('CCO')], [read_smiles('CC=O')],
                          [read_smiles('[H+]')], reagents=[read_smiles('O')])


def test_agents_are_empty_not_absent():
    rxn = ReactionContainer([read_smiles('CCO')], [read_smiles('CC=O')])
    assert rxn.agents == () and _quiet(rxn, 'reagents') == ()


# --- garbage in, stored anyway -------------------------------------------------------------------

def test_an_empty_reaction_is_stored_rather_than_refused():
    """chython 2 RAISES HERE AND THIS CONTAINER STORES.  `ValueError('At least one graph object
    required')` leaves a reader meeting an empty `$RFMT` record in a real RDF with a choice between
    crashing and dropping the record -- and dropping a record is the one thing a reader may never do.
    Input is garbage by default; a container holds it."""
    rxn = ReactionContainer()
    assert len(rxn) == 0 and not rxn
    assert rxn.reactants == () and rxn.products == () and rxn.agents == ()
    assert list(rxn.molecules()) == []


def test_a_one_sided_reaction_is_stored_and_is_falsey():
    """`len` says how many molecules are held; `bool` says whether a transformation is described.
    They are different questions and a half-record answers them differently."""
    rxn = ReactionContainer([read_smiles('CCO')])
    assert len(rxn) == 1 and not rxn
    assert bool(_esterification())


def test_a_non_molecule_is_a_programming_error_and_is_refused():
    """The line between this and the test above: an empty side is a RECORD, a string in the
    reactants list is a CALLER'S MISTAKE, and deferring it hides the line that caused it."""
    with raises(TypeError, match='reactants'):
        ReactionContainer(['CCO'], [read_smiles('CC=O')])
    with raises(TypeError, match='agents'):
        ReactionContainer([read_smiles('CCO')], [read_smiles('CC=O')], [None])


# --- title ---------------------------------------------------------------------------------------

def test_title_is_str_and_empty_means_absent():
    rxn = _esterification()
    assert rxn.title == ''
    rxn.set_title(b'esterification, run 3')
    assert rxn.title == 'esterification, run 3'


def test_a_title_that_is_not_utf8_survives_verbatim():
    """THE REASON THE TYPE USED TO BE BYTES, kept as a promise instead.  A name line is a fixed-width
    field in a text file whose encoding nobody recorded; a byte that does not decode is still what the
    file said, and `surrogateescape` is what carries it through a `str`."""
    ugly = b'run \xff\xfe 3 \x00 <-- ugly on purpose'
    rxn = ReactionContainer(title=ugly)
    assert rxn.title.encode('utf8', 'surrogateescape') == ugly
    assert rxn.copy().title == rxn.title


def test_a_bytes_title_is_decoded_as_a_convenience():
    assert ReactionContainer(title=b'plain ascii').title == 'plain ascii'
    assert ReactionContainer(title='mixed éè').title == 'mixed éè'


def test_a_title_that_is_neither_bytes_nor_str_is_refused():
    with raises(TypeError, match='title'):
        ReactionContainer(title=42)


def test_the_container_has_no_name():
    """The rename, as an absence.  A bytes-valued `name` would make `rxn.name == 'x'` silently False;
    `chython/core/test/test_alternative_spellings.py` carries the full reasoning."""
    assert not hasattr(_esterification(), 'name')


# --- metadata and copies -------------------------------------------------------------------------

def test_meta_is_lazy_and_independent_per_copy():
    rxn = ReactionContainer([read_smiles('CCO')], [read_smiles('CC=O')], meta={'temperature': '80'})
    copy = rxn.copy()
    copy.meta['temperature'] = '120'
    assert rxn.meta['temperature'] == '80'


def test_a_copy_copies_the_molecules_too():
    rxn = _esterification()
    copy = rxn.copy()
    assert all(a is not b for a, b in zip(rxn.molecules(), copy.molecules()))
    assert [m.atom_count for m in copy.molecules()] == [m.atom_count for m in rxn.molecules()]


# --- what was deliberately dropped ---------------------------------------------------------------

def test_the_cgr_surface_is_gone():
    """A CGR existed for reaction ML and the ML consumer never built one -- it read a reaction and
    wanted per-atom, per-side numbers, which `modeling_view` gives directly.  So the overlay
    container and everything that only served it are not ported.  Asserted rather than merely
    omitted, so that re-adding one is a decision with a failing test attached."""
    rxn = _esterification()
    for gone in ('compose', 'decompose', 'centers_list', '__invert__'):
        assert not hasattr(rxn, gone), f'{gone} came back; CGR is not a general-purpose container'


def test_equality_is_a_per_side_multiset_and_the_DEPENDENCY_is_discharged():
    """A reaction is exactly as comparable as a molecule and no more.

    The three semantic questions are answered -- agents participate, mapping does not, `title` does
    not -- so what is here is the multiset comparison and the dependency measurement; the case per
    answer is `test_reaction_identity.py`, one test per decision.

    Why the dependency measurement lives here: it is the reason equality was withheld in the first
    place.  Hashing the canonical reaction SMILES, as V2 does, ties identity to a string, and a
    canonical writer that oscillates on a symmetric stereocentre lets a container compare unequal to
    itself over a round trip.  Reaction equality is sound here only because the mirror automorphism
    behind that oscillation is closed in the canonical search.  The last block MEASURES that on the
    molecules this very reaction is made of, so a regression in molecule identity surfaces as a failure
    next to the thing it would break.
    """
    a, b = _esterification(), _esterification()
    assert a is not b and a == b and hash(a) == hash(b)
    assert len({a, b, a}) == 1

    # The dependency, measured on the sides this reaction is made of.
    left, right = _esterification(), _esterification()
    for x, y in zip(left.reactants + left.products, right.reactants + right.products):
        assert x is not y and x == y and hash(x) == hash(y)
        assert len({x, y}) == 1


def test_repr_is_not_a_chemical_identifier():
    """There is no reaction SMILES here yet, and a nearly-right one would be worse than none: a
    per-molecule CXSMILES tail is meaningless in the middle of a longer string, so concatenating
    three sides' output produces something that looks like a reaction SMILES and is not one."""
    text = repr(_esterification())
    assert '2 reactants' in text and '1 agents' in text and '2 products' in text
    assert '>' not in text


# --- the ML view ---------------------------------------------------------------------------------

def test_the_view_is_keyed_by_map_number_and_carries_per_side_hydrogens():
    """THE WHOLE POINT OF THE VIEW.  The alcohol oxygen (map 3) loses its hydrogen becoming an ester
    oxygen: 1 before, 0 after, on the atom that is the same atom on both sides."""
    view = _esterification().modeling_view()
    element, h_before, n_before, h_after, n_after = view.states[3]
    assert element == 8
    assert (h_before, h_after) == (1, 0), 'the per-side hydrogen count, which is what ML consumes'
    assert (n_before, n_after) == (1, 2), 'and the per-side heavy-atom degree'


def test_an_atom_that_moves_between_molecules_is_still_one_atom():
    """Map 7 is the acid's hydroxyl oxygen and it ends up in the water.  It is MAPPED on both sides,
    so it is not a leaving atom -- it is one atom whose molecule changed, and its per-side numbers say
    so: it loses its bond to the carbonyl carbon and gains a hydrogen.  Pinned because keying on map
    number rather than on the molecule is what makes this work at all."""
    view = _esterification().modeling_view()
    element, h_before, n_before, h_after, n_after = view.states[7]
    assert element == 8
    assert (n_before, n_after) == (1, 0), 'the C-O bond is gone'
    assert (h_before, h_after) == (1, 2), 'and a hydrogen arrived'


def test_a_leaving_fragment_keeps_its_internal_bonds_and_its_hydrogens():
    """A Suzuki coupling's boronic acid leaves as a fragment.  The convention -- a modelling choice,
    not chemistry -- is that the fragment departs intact: boron keeps its two B-O bonds, loses only
    the bond to the ring it left, and its hydrogen count is not ours to invent because the record
    never says what the fragment became."""
    view = _suzuki().modeling_view()
    element, h_before, n_before, h_after, n_after = view.states[17]
    assert element == 5
    assert (n_before, n_after) == (3, 2), 'three bonds before, the two B-O after'
    assert h_after == h_before
    assert view.states[7][2:] == (1, 0, 0), 'and a lone leaving bromine keeps nothing'


def test_union_bonds_shows_a_bond_broken_and_a_bond_formed():
    view = _esterification().modeling_view()
    assert view.union_bonds[(5, 7)] == (1, 0), 'C-OH broken'
    assert view.union_bonds[(3, 5)] == (0, 1), 'C-O formed'
    assert view.union_bonds[(2, 3)] == (1, 1), 'and one that does not change'


def test_agents_contribute_nothing_to_the_view():
    """An agent is by definition not consumed, so it carries no signal and would make the atom count
    depend on how the record's author split the left side."""
    with_agent = _esterification()
    without = ReactionContainer(with_agent.reactants, with_agent.products)
    assert with_agent.modeling_view().states == without.modeling_view().states
    assert with_agent.modeling_view().union_bonds == without.modeling_view().union_bonds


def test_an_unstated_hydrogen_count_reports_the_sentinel_not_zero():
    """`h or 0` would train a molecule with an unstated hydrogen count as though it had none -- a
    PLAUSIBLE wrong number.  H_UNKNOWN is outside the 0..14 domain, so ignoring it produces an
    obviously broken value instead."""
    mol = read_smiles('[CH3:1][CH:2]=[O:3]')
    mol.set_hydrogens(2, H_UNKNOWN)  # the core's own way of writing "the record did not say"
    view = ReactionContainer([mol], [read_smiles('[CH3:1][C:2](=[O:3])[OH:4]')]).modeling_view()
    assert view.states[2][1] == H_UNKNOWN
    assert view.states[1][1] == 3, 'and a stated count is still the stated count'


def test_unmapped_atoms_are_counted_rather_than_silently_merged():
    """Leaving unmapped atoms undefined means, in a training set, silent corruption at whatever rate
    the corpus happens to contain.  Counted here, so a pipeline can refuse at ITS boundary."""
    rxn = ReactionContainer([read_smiles('[CH3:1][CH2:2]O')], [read_smiles('[CH3:1][CH:2]=O')])
    view = rxn.modeling_view()
    assert view.unmapped == {'reactants': 1, 'products': 1}
    assert 0 not in view.states, 'map number 0 is not an atom identity'
    assert view.collisions == {'reactants': (), 'products': ()}


def test_a_bond_to_an_unmapped_atom_is_left_out_of_the_union():
    rxn = ReactionContainer([read_smiles('[CH3:1][CH2:2]O')], [read_smiles('[CH3:1][CH:2]=O')])
    view = rxn.modeling_view()
    assert set(view.union_bonds) == {(1, 2)}


def test_colliding_map_numbers_are_reported():
    """Two atoms on one side claiming one map number means the union merged them.  A record is
    allowed to say that; a container is not allowed to hide it."""
    rxn = ReactionContainer([read_smiles('[CH3:1][CH3:1]')], [read_smiles('[CH3:1][CH2:2][OH:3]')])
    view = rxn.modeling_view()
    assert view.collisions['reactants'] == (1,)
    assert view.collisions['products'] == ()


def test_the_view_of_an_empty_reaction_is_empty_and_does_not_raise():
    view = ReactionContainer().modeling_view()
    assert view.states == {} and view.union_bonds == {}
    assert isinstance(view, ReactionModelingView) and 'atoms' in repr(view)


def test_an_unchanged_atom_far_from_the_reaction_centre_reads_as_unchanged():
    """A Suzuki coupling: the boron leaves, the biaryl bond forms, and the ring carbons that took no
    part must not acquire a spurious change."""
    view = _suzuki().modeling_view()
    element, h_before, n_before, h_after, n_after = view.states[3]
    assert element == 6 and (h_before, h_after) == (1, 1) and (n_before, n_after) == (2, 2)
    assert view.union_bonds[(6, 7)] == (1, 0), 'C-Br broken'
    assert view.union_bonds[(6, 16)] == (0, 1), 'biaryl bond formed'
    assert view.union_bonds[(1, 2)] == (4, 4), 'and an aromatic bond stays aromatic on both sides'


# --- the announcement ----------------------------------------------------------------------------

def test_reagents_warns_and_names_agents():
    rxn = _esterification()
    with warns(DeprecationWarning, match='agents') as record:
        rxn.reagents
    assert len(record) == 1 and 'reagents' in str(record[0].message)


def test_agents_and_title_are_silent():
    rxn = _esterification()
    with catch_warnings():
        simplefilter('error', DeprecationWarning)
        rxn.agents, rxn.reactants, rxn.products, rxn.title, rxn.meta
        rxn.set_title(b'ported')
        rxn.modeling_view()


def test_the_warning_is_charged_to_the_caller():
    """MEASURED, NOT ASSUMED -- `stacklevel` counts Python frames, and the same intent needs 1 inside
    the compiled core and 3 here."""
    from inspect import currentframe

    rxn = _esterification()
    with warns(DeprecationWarning) as record:
        here = currentframe().f_lineno + 1
        rxn.reagents
    assert record[0].filename == __file__ and record[0].lineno == here, \
        f'blamed {record[0].filename}:{record[0].lineno}, wanted this file:{here}'


# --- chython 2 as an independent witness ---------------------------------------------------------
#
# A SEPARATE INTERPRETER, AND `-I` IS LOAD-BEARING.  Without it the child prepends this working tree
# to `sys.path` and imports the code under test, so the differential compares the tree to itself and
# always agrees.
#
# This is a witness, not an oracle: it catches the container disagreeing with chython 2 on the things
# that are meant to mean the same on both sides.  It goes when the version pin stops resolving.

ORACLE_SCRIPT = r'''
from chython import smiles

out = []
for s in _payload:
    r = smiles(s)
    out.append({'reactants': [str(m) for m in r.reactants],
                'agents': [str(m) for m in r.reagents],
                'products': [str(m) for m in r.products],
                'title': r.name,
                'maps': [sorted(n for n in m) for m in r.molecules()]})
_emit(out)
'''

REACTIONS = ['[CH3:1][CH2:2][OH:3].[CH3:4][C:5](=[O:6])[OH:7]>[H+]>'
             '[CH3:1][CH2:2][O:3][C:5](=[O:6])[CH3:4].[OH2:7]',
             '[CH3:1][C:2](=[O:3])[OH:4].[NH2:5][CH3:6]>>[CH3:1][C:2](=[O:3])[NH:5][CH3:6].[OH2:4]',
             'CCO>>CC=O']


def _oracle(*reactions):
    """chython 2's split of each reaction string.  The reactions travel as the payload, so there is
    no temporary script file on disk and no argv quoting to get wrong.

    The isolation flag, the version pin and the not-this-tree check are `oracle.py`'s and not this
    file's.  Spelled out per test file, any one of them missing makes the differential pass vacuously
    instead of fail.
    """
    return oracle.ask(ORACLE_SCRIPT, list(reactions))


def test_the_oracle_is_really_a_different_library():
    """FIRST, A NEGATIVE CONTROL.  The version pin and the absence of this tree, checked before any
    comparison is believed -- a differential against yourself passes and means nothing.

    `oracle.verify` is what every `_oracle` call above already runs, so this is a statement of intent
    in the file that depends on it rather than the only thing enforcing it.  The guards' own tests --
    including the control that drops `-I` and watches an unisolated child import this tree -- are in
    `test_oracle.py`.
    """
    oracle.require()
    oracle.verify()


@mark.parametrize('index', range(len(REACTIONS)))
def test_the_three_sides_agree_with_chython_two(index):
    """Same reaction string, same split into three sides.  What is compared is the PARTITION, not the
    SMILES text -- the two libraries' canonical writers are different code and string identity is
    not the claim (and this repository has already ruled that it is not sound anyway)."""
    reference = _oracle(REACTIONS[index])[0]
    rxn = _rebuild(REACTIONS[index])
    assert len(rxn.reactants) == len(reference['reactants'])
    assert len(rxn.agents) == len(reference['agents'])
    assert len(rxn.products) == len(reference['products'])


def test_a_fully_stated_atom_map_agrees_with_chython_two():
    """The amide formation, where every atom in the string carries a map.  Where the two libraries can
    both be believed, they agree; the test below is about where only one of them can."""
    reference = _oracle(REACTIONS[1])[0]
    rxn = _rebuild(REACTIONS[1])
    ours = [sorted(a.map_number for a in m.atoms()) for m in rxn.molecules()]
    assert ours == reference['maps']


@mark.parametrize('index,invented', [(0, [[8]]), (2, [[1, 2, 3], [4, 5, 6]])])
def test_chython_two_invents_map_numbers_where_this_container_reports_none(index, invented):
    """A DIVERGENCE THE DIFFERENTIAL PINS, AND THE PORTING NOTE THAT GOES WITH IT.

    chython 2 has one integer per atom doing two jobs: the atom's identity within its container and
    its atom-atom map number.  So an atom the string left unmapped still comes back carrying a number,
    continuing whatever sequence the mapped atoms established, and a consumer cannot tell that number
    from a stated one -- which for a reaction is the difference between "these two are the same atom"
    and "nobody said".  The core keeps the two apart: `n` is identity, `map_number` is what the record
    stated, 0 when it stated nothing.

    Both shapes are here.  Reaction 0 is PARTIALLY mapped -- only its `[H+]` agent is bare, and
    chython 2 hands it 8, the next number after the seven that were stated, indistinguishable from part
    of the map.  Reaction 2 is mapped nowhere and comes back mapped throughout.

    This is why `modeling_view` counts unmapped atoms rather than folding them into the union.  If the
    numbers were believable there would be nothing to count.
    """
    reference = _oracle(REACTIONS[index])[0]
    rxn = _rebuild(REACTIONS[index])
    ours = [sorted(a.map_number for a in m.atoms()) for m in rxn.molecules()]
    theirs = reference['maps']

    bare = [(o, t) for o, t in zip(ours, theirs) if set(o) == {0}]
    assert [t for _, t in bare] == invented, 'chython 2 numbered what the record left bare'

    # and every bare atom is ACCOUNTED FOR: either counted by the view, or on the agents side, which
    # the view excludes by design.  Nothing goes missing without a number attached to it.
    agent_atoms = sum(m.atom_count for m in rxn.agents)
    assert sum(len(o) for o, _ in bare) == sum(rxn.modeling_view().unmapped.values()) + agent_atoms


def test_the_title_arrives_identically_including_when_it_is_empty():
    """The rename, checked against the library it renames.  chython 2's `name` normalises absent to
    `''`, and so does this `title` -- one type and one spelling for the same statement."""
    reference = _oracle(*REACTIONS)
    for record in reference:
        assert record['title'] == ''
    assert _rebuild(REACTIONS[0]).title == ''


def _rebuild(reaction_smiles):
    """Build the chython 3 container from a reaction SMILES.

    One line, calling the core's reaction reader, so the differential tests above compare V2's whole
    reader against V3's rather than against a splitter written in a test file.  The name is kept
    because these tests read better with it than with the reader's.
    """
    return read_reaction_smiles(reaction_smiles)
