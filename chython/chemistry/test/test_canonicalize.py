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
"""`canonicalize()`, `implicify_hydrogens()` and `explicify_hydrogens()`.

`canonical_bytes` is computed from what a molecule stores, so two drawings of one compound do not hash
equal until this pass has run.  Deduplication is therefore tested first, and the properties below it --
idempotence, never raising, a legible log -- are what that purpose needs."""
from pytest import raises

from .. import canonicalize, explicify_hydrogens, implicify_hydrogens   # the import injects the method
from ...core import LOST, MoleculeContainer, REFUSED, REPAIRED, Log, read_smiles as smiles


#: Pairs measured to store different bytes and mean the same compound.  Each names the stage that
#: closes it.
DEDUP_PAIRS = [
    ('c1ccccc1O', 'C1=CC=CC=C1O'),        # aromatic vs Kekule -- thiele()
    ('c1ccncc1', 'C1=CC=NC=C1'),          # ... with a heteroatom
    ('c1ccc2ccccc2c1', 'C1=CC=CC2=CC=CC=C12'),   # ... fused
    ('[CH4]', '[H]C([H])([H])[H]'),       # hydrogen count vs hydrogen atoms -- implicify_hydrogens()
    ('CCO', 'CCO[H]'),                    # ... one hydroxyl hydrogen, as MDL records carry it
    ('C[N+](=O)[O-]', 'CN(=O)=O'),        # nitro drawn two ways -- standardize()
]

#: The lines a plain rewrite writes: every aromatic system that comes out of `kekule()` or `thiele()`
#: says so, once per system per call.  The pipeline runs both twice on a molecule that needs a second
#: round, so the assertions below count repairs and not log entries.
ROUTINE = frozenset(('kekule:kekulized', 'thiele:aromatized'))


def repairs(log):
    """What the PIPELINE repaired: the routine rewrite notices out, and the reader's own lines with them.

    `mol.log` is one storage and the reader writes there too, so a record's stage is what says which
    call produced it.  Every assertion below is about `canonicalize()`.
    """
    return [r for r in log if r.stage != 'read' and r.rule not in ROUTINE]


def test_two_drawings_of_one_compound_hash_equal_only_after_canonicalizing():
    """The point of the pass: without it, `canonical_bytes` answers "same drawing"."""
    for a, b in DEDUP_PAIRS:
        ma, mb = smiles(a), smiles(b)
        assert ma.canonical_bytes != mb.canonical_bytes, f'{a} / {b} already store the same bytes'
        ma.canonicalize()
        mb.canonicalize()
        assert ma.canonical_bytes == mb.canonical_bytes, f'{a} / {b} still differ after canonicalize'


def test_local_tautomers_do_unify_because_the_group_rules_reach_them():
    """`standardize()`'s `fix_tautomers` rules unify pairs whose forms differ by a hydrogen moving
    between two heavy atoms one bond apart, so these dedup.  The gap is narrower than "tautomers".
    """
    for a, b in [('Oc1ccccn1', 'O=c1cccc[nH]1'),        # 2-hydroxypyridine / 2-pyridone
                 ('CC(=O)C', 'CC(O)=C'),                # keto / enol
                 ('CC(=O)CC(=O)C', 'CC(O)=CC(=O)C'),    # ... of a 1,3-diketone
                 ('OC=CC', 'O=CCC')]:
        ma, mb = smiles(a), smiles(b)
        ma.canonicalize()
        mb.canonicalize()
        assert ma.canonical_bytes == mb.canonical_bytes, f'{a} / {b} no longer dedup'


def test_a_prototropic_shift_around_a_ring_now_agrees_and_the_table_is_closed():
    """The two N-H forms of 4-methylimidazole hash equal, and it takes two mechanisms to get there.

    No local tautomer rule sees this pair -- the hydrogen moves three bonds and the double bonds move
    with it.  What closes it is `thiele()` giving back the aromatic form and `standardize_isomers()`
    choosing the placement in a frame that does not depend on which form arrived.
    """
    ma, mb = smiles('Cc1cnc[nH]1'), smiles('Cc1c[nH]cn1')
    assert ma.canonical_bytes != mb.canonical_bytes
    ma.canonicalize()
    mb.canonicalize()
    assert ma.canonical_bytes == mb.canonical_bytes, (str(ma), str(mb))
    assert ma == mb


def test_the_ring_shift_agrees_for_a_fused_system_and_for_an_anion_too():
    for x, y in [('c1cc2[nH]ncc2cn1', 'c1cc2n[nH]cc2cn1'),      # a pyrazolo-fused pyridine
                 ('c1ccc2[nH]ncc2c1', 'c1ccc2n[nH]cc2c1'),      # indazole
                 ('Cc1cnc[n-]1', 'Cc1c[n-]cn1')]:               # a mobile CHARGE, not a hydrogen
        ma, mb = smiles(x), smiles(y)
        ma.canonicalize()
        mb.canonicalize()
        assert ma.canonical_bytes == mb.canonical_bytes, (x, y, str(ma), str(mb))


def test_keep_kekule_changes_the_spelling_and_never_the_compound():
    """The flag must not decide which tautomer you get, only how the bonds are written.

    Which is why step 6 undoes step 4's aromatic form with a second `kekule()` instead of the pipeline
    skipping the placement stage.
    """
    ma, mb = smiles('Cc1c[nH]cn1'), smiles('Cc1cnc[nH]1')
    ma.canonicalize(keep_kekule=True)
    mb.canonicalize(keep_kekule=True)
    assert ma.canonical_bytes == mb.canonical_bytes, (str(ma), str(mb))
    plain = smiles('Cc1c[nH]cn1')
    plain.canonicalize()
    aromatised = ma.copy()
    aromatised.thiele()
    assert aromatised.canonical_bytes == plain.canonical_bytes


def test_fix_tautomers_false_does_not_reach_the_placement_stage():
    """The flag withholds the local repair rules; placement is not a repair, so it is not withheld.

    The local rules decide whether a drawing was wrong, the placement stage decides which of two right
    drawings to keep, and a caller distrusting the former has said nothing about the latter.
    """
    ma, mb = smiles('Cc1c[nH]cn1'), smiles('Cc1cnc[nH]1')
    ma.canonicalize(fix_tautomers=False)
    mb.canonicalize(fix_tautomers=False)
    assert ma.canonical_bytes == mb.canonical_bytes


def test_the_placement_stage_names_itself_in_the_log_and_runs_after_implicify():
    """Placement runs after implicify, which is why the fixture spells the hydrogen out.

    An explicit hydrogen atom makes the ring nitrogen three-coordinate, and a three-coordinate atom is
    not a placement site.  So placement sees one site and declines unless `implicify_hydrogens()` has
    already folded that atom into a count.
    """
    mol = smiles('Cc1cnc[n]1[H]')
    log = mol.log
    mol.canonicalize()
    stages = [r.stage for r in log]
    assert 'isomers' in stages, stages
    assert stages.index('isomers') > stages.index('implicify'), stages

    other = smiles('Cc1c[nH]cn1')                # and the compound the two stages agree on is the same
    other.canonicalize()
    assert mol.canonical_bytes == other.canonical_bytes


def test_turning_fix_tautomers_off_gives_up_part_of_the_dedup_guarantee():
    """Which is why it defaults to on.  The flag is not free and the caller should know the price."""
    ma, mb = smiles('Oc1ccccn1'), smiles('O=c1cccc[nH]1')
    ma.canonicalize(fix_tautomers=False)
    mb.canonicalize(fix_tautomers=False)
    assert ma.canonical_bytes != mb.canonical_bytes


def test_the_bool_says_changed_and_a_second_pass_says_no():
    """The output is a hash input, so a second pass must be a no-op that says so.

    Steps 1 and 5 are a round trip, so the bool cannot just forward `thiele()`'s -- that one is truthy
    whenever an aromatic ring exists.
    """
    unchanged = ['c1ccccc1O', 'c1ccncc1', 'CC(=O)[O-].[Na+]', 'C[N+](=O)[O-]', 'B1[H]B1']
    changed = ['C1=CC=CC=C1O', '[H]C([H])([H])[H]', 'CCO[H]', 'CN(=O)=O', '[H]c1cc[n](=O)cc1']
    for s in unchanged:
        m = smiles(s)
        assert m.canonicalize() is False, f'{s} is already canonical but reported a change'
    for s in changed:
        m = smiles(s)
        assert m.canonicalize() is True, f'{s} is not canonical but reported no change'
        assert m.canonicalize() is False, f'{s} reported a change on an idempotent second pass'


#: One compound drawn twice, differing only in which ring nitrogen holds the mobile hydrogen, where one
#: drawing puts it on the very nitrogen the hydroxy-azine rows need free.  These need the loop: repair
#: alone declines on the second member, and placement alone leaves the enol standing.
UNBLOCKED_PAIRS = [
    ('Sc1ncnc2[nH]cnc12', 'Sc1[nH]cnc2ncnc1-2', '6-mercaptopurine'),
    ('Oc1nc(O)c2nn[nH]c2n1', 'Oc1nc(O)c2nnnc-2[nH]1', '8-azaxanthine'),
]


def test_a_placement_that_unblocks_a_repair_is_repaired():
    """Steps 2 to 5 are a loop because step 5 can hand step 2 work it declined.

    A `tautomer` row turning a hydroxy-azine into the oxo form needs the ring nitrogen beside the C-OH
    free, and on the second member of each pair the mobile hydrogen is sitting on it.  Run once, that
    drawing keeps its hydroxy form while its twin gets the oxo form, and one compound gets two keys.
    """
    for a, b, label in UNBLOCKED_PAIRS:
        ma, mb = smiles(a), smiles(b)
        ma.canonicalize()
        mb.canonicalize()
        assert ma.canonical_bytes == mb.canonical_bytes, f'{label}: {ma} vs {mb}'


def test_the_second_round_kekulizes_before_it_repairs():
    """Not a spelling detail -- without it the loop is a guaranteed no-op.

    The `tautomer` rows are written against definite bond orders, so on the aromatic form step 4 leaves
    behind they match nothing.  Reaching the oxo form here is the only observable proof that the second
    round re-kekulized rather than re-running the rules against an aromatic ring.
    """
    m = smiles('Sc1[nH]cnc2ncnc1-2')
    m.canonicalize()
    thione = [n for n in m.atom_numbers
              if m.element_of(n) == 16 and m.order_of(n, next(iter(m.neighbors_of(n)))) == 2]
    assert thione, f'the thiol was never repaired to the thione: {m}'


def test_the_pipeline_is_a_fixed_point_and_not_merely_ordered():
    """The property the loop exists for, and the one a hash key actually needs.

    Idempotence was already pinned for shapes that converge in one pass; these converge in two, so before
    the loop a second `canonicalize()` moved them again -- which means the first answer was not the
    canonical form and which of two callers got it depended on how many times they had asked.
    """
    for a, b, label in UNBLOCKED_PAIRS:
        for string in (a, b):
            m = smiles(string)
            m.canonicalize()
            first = m.canonical_bytes
            assert m.canonicalize() is False, f'{label}: {string} still moved on a second pass'
            assert m.canonical_bytes == first, f'{label}: {string} is not a fixed point'


def test_the_extra_round_costs_no_duplicate_repair_record():
    """A repair is reported once however many rounds saw the molecule.

    Step 1 already reported every repair the kekuliser can find, so a caller counting `repaired()`
    records must not see a stage run twice rather than a molecule repaired twice.  The routine rewrite
    notices are per call by design -- each round really did rewrite the representation.
    """
    m = smiles('Sc1[nH]cnc2ncnc1-2')
    log = m.log
    m.canonicalize()
    rules = [r.rule for r in repairs(log)]
    assert len(rules) == len(set(rules)), f'a rule fired into the log twice: {rules}'
    assert not [r for r in log if r.rule == 'canonicalize:rounds'], 'the round cap was reached'


def test_a_bridging_hydride_is_a_record_and_not_an_exception():
    """A repair pass is not an answer boundary, so it has nothing to refuse the caller with.

    The molecule comes back untouched, the reason is in the log, and the other stages still ran.
    """
    m = smiles('B1[H]B1')
    log = m.log
    assert m.canonicalize() is False
    assert len(m) == 3, 'the bridging hydrogen was consumed'
    refused = log.refused()
    assert len(refused) == 1
    assert refused[0].rule == 'hydrogens:implicify-bridging'
    assert refused[0].atoms == (2,)
    assert 'bridges 2 atoms' in refused[0]


def test_more_explicit_hydrogens_than_a_count_holds_is_a_record_and_not_an_exception():
    """The count field reaches 14, and a record may draw more hydrogen ATOMS than that on one atom.

    Fourteen fold in, the rest stay atoms with a record each, and the decision is taken before the edit
    session opens -- a write refused from inside it would leave the molecule half-implicified.
    """
    m = MoleculeContainer()
    with m.edit():
        boron = m.add_atom('B', implicit_h=0)
        for _ in range(17):
            m.add_bond(boron, m.add_atom(1, implicit_h=0), 1)
    log = m.log
    assert m.canonicalize() is True
    assert m.implicit_h_of(boron) == 14, 'the count is filled, not overflowed'
    assert len(m) == 4, 'the three hydrogens the count cannot hold are still atoms'
    refused = log.refused()
    assert len(refused) == 3
    assert {r.rule for r in refused} == {'hydrogens:implicify-count-full'}
    assert 'records at most 14' in refused[0]


def test_the_pass_states_the_type_it_takes():
    """`smiles()` returns whichever container its string describes, so a `>>` in a structure column
    reaches this pass as a reaction.  An `AttributeError` naming a private attribute is not a refusal.
    """
    reaction = smiles('C1=CC=CC=C1O.CO>>CC(=O)OC')
    with raises(TypeError, match='takes a MoleculeContainer'):
        canonicalize(reaction)
    assert reaction.canonicalize() is True, 'the reaction canonicalizes through its own method'


def test_the_stages_run_in_the_designed_order_and_the_log_proves_it():
    """`standardize()` must see a Kekule structure, so `kekule()` precedes it.

    This is an order pin, not an outcome pin: no candidate molecule is unkekulizable before standardize
    and kekulizable after, because the rules that repair a ring atom's drawing are not in the table
    yet.  TODO: upgrade to an outcome pin when a rule lands whose repair unblocks kekulization.
    """
    # Two components in one container so all three logging stages have something to say: an N-oxide
    # pyridine that kekulizes only charge-separated, a pentavalent nitro for the group table, and one
    # explicit hydrogen to fold.
    m = smiles('[H]c1cc[n](=O)cc1.CN(=O)=O')
    log = m.log
    assert m.canonicalize() is True
    assert [r.stage for r in repairs(log)] == ['kekule', 'standardize', 'implicify']


def test_an_unkekulizable_ring_is_logged_and_the_rest_of_the_pipeline_still_runs():
    """Garbage in, logged, not refused, and the degradation is visible rather than silent.

    `C[n+]1cccc1` has no Kekule form as drawn: it comes back with a `LOST` record naming the ring
    system, and step 1's failure did not cancel steps 2 and 3.  The fixture cannot be `Cn1cc[nH]c1`,
    which the kekuliser now repairs by dropping the surplus hydrogen; a cationic three-coordinate
    nitrogen is must-match, so five must-match atoms is an odd count nothing makes even.
    """
    m = smiles('C[n+]1cccc1')
    log = m.log
    assert m.canonicalize() is True
    lost = log.lost()
    assert len(lost) == 1
    assert lost[0].stage == 'kekule'
    assert lost[0].atoms == (2, 3, 4, 5, 6)
    assert 'no Kekule form' in lost[0]
    assert len(m) == 6, 'the molecule was dropped rather than reported'


def test_keep_kekule_returns_the_kekule_form_and_costs_the_log_nothing_but_the_rewrite():
    """It is a correctness mechanism and nothing else: it repairs nothing and refuses nothing.

    The one line it does cost is the notice for the closing rewrite back to a Kekule form, which is
    work that happened and so is reported like any other.
    """
    m = smiles('c1ccccc1O')
    log = m.log
    m.canonicalize(keep_kekule=True)
    assert format(m) == 'C1(=CC=CC=C1)O'
    plain = smiles('c1ccccc1O')
    plain_log = plain.log
    plain.canonicalize()
    assert [str(r) for r in repairs(log)] == [str(r) for r in repairs(plain_log)], \
        'a repair depended on keep_kekule'
    assert [r.rule for r in log] == [r.rule for r in plain_log] + ['kekule:kekulized']


def test_a_plain_list_gets_what_it_gets_today():
    """A caller who never heard of `Log` is not broken by it.

    The substring idiom must keep working whether the entry is the `str` a reader appends or the
    `LogRecord` a pass builds.
    """
    m = smiles('CCO[H]')
    log = m.log
    assert m.canonicalize() is True
    assert len(log) == 1
    assert 'folded into its count' in log[0]


def test_a_log_gives_one_filterable_record_type_out_of_five_stages():
    """Three log mechanics meet here -- bare strings, `LogRecord`s, and result objects that return
    their own log -- and the caller still sees one type, stage-tagged and filterable by severity.
    """
    m = smiles('[H]c1cc[n](=O)cc1.CN(=O)=O')
    log = m.log
    m.canonicalize()

    assert len({type(r) for r in log}) == 1, 'the caller needs isinstance to read their own log'
    # kekule's record came off a KekuleResult, standardize's is a LogRecord the pass built; both
    # arrive stamped.  The rule is the kekuliser's own and NOT `canonicalize:kekule`: `absorb` fills
    # only blank provenance, so a pass that named its rule keeps it and the stage is what says which
    # pipeline it ran in.
    assert log.by_stage('kekule')[0].rule.startswith('kekule:')
    # the table only, not the row: `standardize_groups.tsv` ids are positional and get renumbered
    # whenever rows merge, so a pinned `groups:13` would fail on somebody else's table edit
    assert log.by_stage('standardize')[0].rule.startswith('groups:')
    assert log.by_stage('implicify')[0].rule == 'hydrogens:implicify'
    assert len(log.repaired()) == 2
    assert log.atoms_touched('implicify') == {2}
    assert log.atoms_touched('standardize') == {9, 10, 11, 12}


def test_nothing_is_logged_when_no_log_is_asked_for():
    """`log=None` is forwarded as `None`, so no stage builds a record it cannot deliver.

    There is no assertion available for "did not allocate", so this is the observable half: the bool is
    identical either way.
    """
    for s in [s for pair in DEDUP_PAIRS for s in pair]:
        a, b = smiles(s), smiles(s)
        assert a.canonicalize() == b.canonicalize()


# --- `implicify_hydrogens` on its own: the stage that touches the graph, so its refusals and its
#     stereo handling are pinned here rather than through the façade.


def test_an_ordinary_hydrogen_is_folded():
    for s, expect, went in [('CCO[H]', 'C(C)O', 1), ('[H]OCC', 'C(C)O', 1),
                            ('[H]C([H])([H])[H]', 'C', 4)]:
        m = smiles(s)
        assert implicify_hydrogens(m) == went
        assert format(m) == expect


def test_the_five_kinds_that_are_not_a_hydrogen_count():
    """Each would lose information a count cannot carry, so each is left as an atom."""
    for s, why in [('[2H]C', 'an isotope is a label a count has no room for'),
                   ('[H-].[Na+]', 'a hydride is a compound'),
                   ('[H+].[Cl-]', 'a proton is a compound'),
                   ('[H]C |^1:0|', 'a hydrogen radical is a compound'),
                   ('[H][H]', 'neither H2 atom has a heavy neighbour to fold into'),
                   ('B1[H]B1', 'a bridging hydride is not any one atom\'s count')]:
        m = smiles(s)
        n = len(m)
        assert implicify_hydrogens(m) == 0, why
        assert len(m) == n, why


def test_a_hydrogen_held_by_something_other_than_a_single_bond_is_refused():
    """Garbage input can spell one, and a count cannot record the order."""
    m = smiles('C(=[H])C')
    log = m.log
    assert implicify_hydrogens(m) == 0
    assert len(log.refused()) == 1
    assert 'a hydrogen count cannot record that' in log.refused()[0]


def test_folding_into_an_unknown_count_is_refused_rather_than_laundered():
    """`[H]I(C)(C)C`'s tetravalent iodine derives no count, so folding into it would report a total
    that is not known.  A `LOST` record, not a number.
    """
    m = smiles('[H]I(C)(C)C')
    assert m.implicit_h_of(2) is None, 'the premise of this test moved: iodine now derives a count'
    log = m.log
    assert implicify_hydrogens(m) == 0
    assert len(m) == 5
    lost = [r for r in log.lost() if r.stage == 'implicify']   # the reader lost the valence rule first
    assert len(lost) == 1
    assert lost[0].rule == 'hydrogens:implicify-unknown-count'
    assert lost[0].severity == LOST


def test_a_tetrahedral_centre_is_not_racemised_by_losing_its_hydrogen_atom():
    """`delete_atom` alone clears the parity, so the pass captures it first and writes it back.

    Both enantiomers, because a bug that maps both to one answer passes a one-sided test.
    """
    for s, ref in [('F[C@]([H])(Cl)Br', 'F[C@H](Cl)Br'), ('F[C@@]([H])(Cl)Br', 'F[C@@H](Cl)Br')]:
        m = smiles(s)
        assert implicify_hydrogens(m) == 1
        assert m.canonical_bytes == smiles(ref).canonical_bytes, f'{s} lost or inverted its parity'
    # and the two answers are still different from each other
    a, b = smiles('F[C@]([H])(Cl)Br'), smiles('F[C@@]([H])(Cl)Br')
    implicify_hydrogens(a)
    implicify_hydrogens(b)
    assert a.canonical_bytes != b.canonical_bytes


def test_double_bond_stereo_survives_without_being_restored():
    """Unlike a tetrahedral centre: the core re-derives bond stereo on seal and gets it right.

    Which is why there is no bond-parity restore in the pass and no bond-parity setter to reach for.
    """
    for s, ref in [('F/C=C([H])\\F', 'F/C=C\\F'), ('F/C=C([H])/F', 'F/C=C/F')]:
        m = smiles(s)
        assert implicify_hydrogens(m) == 1
        assert m.canonical_bytes == smiles(ref).canonical_bytes, f'{s} changed diastereomer'


def test_the_return_is_a_count_of_atoms_whether_or_not_a_log_was_given():
    """The return type does not depend on `log=`.

    It counts atoms REMOVED, not anchors touched, which is why it is returned at all rather than
    recovered from the log -- methane's four hydrogens are one record naming one carbon.
    """
    m = smiles('[H]C([H])([H])[H]')
    log = m.log
    assert implicify_hydrogens(m) == 4
    assert log.repaired()[0].severity == REPAIRED
    assert len(log) == 1, 'one record per anchor, not one per hydrogen and not a summary'
    assert log.atoms_touched() == {2}
    assert implicify_hydrogens(m) == 0
    # and the count is not a bool wearing an int's clothes
    assert implicify_hydrogens(smiles('CCO[H]')) is not True


def test_a_refusal_is_never_an_exception_and_takes_no_flag_to_become_one():
    """There is no `ignore=` parameter, by design: the event decides what it is, not the caller."""
    from inspect import signature
    assert 'ignore' not in signature(implicify_hydrogens).parameters
    m = smiles('B1[H]B1')
    log = m.log
    implicify_hydrogens(m)      # must not raise
    assert log.refused()[0].severity == REFUSED


# --- `explicify_hydrogens`, the other direction.  Not a canonicalize stage -- nothing canonical wants
#     five atoms where one will do -- so it is tested only here.


def test_a_count_becomes_atoms():
    for s, expect, arrived in [('C', '[H]C([H])([H])[H]', 4),
                               ('CCO', 'O(C([H])(C([H])([H])[H])[H])[H]', 6)]:
        m = smiles(s)
        assert explicify_hydrogens(m) == arrived
        assert format(m) == expect


def test_the_two_passes_are_inverse_over_a_corpus():
    """Explicify then implicify returns the molecule it started from, bytes for bytes.

    Stereocentres included, which is where a hydrogen round trip goes wrong if it goes wrong at all.
    """
    corpus = ['C', 'CCO', 'c1ccccc1O', 'c1ccncc1', 'CC(=O)[O-].[Na+]', 'C[N+](=O)[O-]',
              'C[C@H](N)C(=O)O', 'C[C@@H](N)C(=O)O', 'N[C@@H](Cc1ccccc1)C(O)=O',
              'F[C@H](Cl)Br', 'F[C@@H](Cl)Br', 'F/C=C/F', 'F/C=C\\F', 'OC=CC',
              'c1ccc2ccccc2c1', 'CC(=O)CC(=O)C', 'B1[H]B1', '[2H]C', '[H-].[Na+]']
    for s in corpus:
        m = smiles(s)
        before = m.canonical_bytes
        added = explicify_hydrogens(m)
        removed = implicify_hydrogens(m)
        assert m.canonical_bytes == before, f'{s} did not survive the round trip'
        assert added == removed, f'{s} added {added} hydrogens and gave back {removed}'


def test_a_stereocentre_survives_explicification_without_a_parity_restore():
    """The measured asymmetry: `delete_atom` clears a parity, `add_atom` does not.

    So implicify captures and restores and explicify does not, and the missing restore is deliberate.
    Against a hand-written explicit reference and against its epimer, because a pass that racemised
    both enantiomers to one answer would pass a one-sided test.
    """
    for s, ref, epimer in [('C[C@H](N)C(=O)O', '[H]C([H])([H])[C@]([H])(N([H])[H])C(=O)O[H]',
                            '[H]C([H])([H])[C@@]([H])(N([H])[H])C(=O)O[H]'),
                           ('C[C@@H](N)C(=O)O', '[H]C([H])([H])[C@@]([H])(N([H])[H])C(=O)O[H]',
                            '[H]C([H])([H])[C@]([H])(N([H])[H])C(=O)O[H]')]:
        m = smiles(s)
        explicify_hydrogens(m)
        assert m.canonical_bytes == smiles(ref).canonical_bytes, f'{s} lost or inverted its parity'
        assert m.canonical_bytes != smiles(epimer).canonical_bytes


def test_an_unknown_count_yields_no_atoms_and_a_lost_record():
    """`[H]I(C)(C)C`'s iodine derives no count, so there is no number to write out.

    Inventing zero would answer a question the record never answered.
    """
    m = smiles('[H]I(C)(C)C')
    assert m.implicit_h_of(2) is None, 'the premise of this test moved: iodine now derives a count'
    log = m.log
    added = explicify_hydrogens(m)
    lost = [r for r in log.lost() if r.atoms == (2,)]
    assert len(lost) == 1
    assert lost[0].rule == 'hydrogens:explicify-unknown-count'
    assert 'unknown implicit hydrogen count' in lost[0]
    # the other atoms were still served: the three methyls got their nine hydrogens
    assert added == 9


def test_the_new_hydrogens_are_unmapped_and_there_is_no_keyword_to_number_them():
    """A hydrogen this pass invented has no counterpart on the other side of anything, so a map number
    would assert a correspondence that does not exist.  The caller who needs one assigns it.
    """
    from inspect import signature
    params = signature(explicify_hydrogens).parameters
    assert set(params) == {'molecule'}, 'a numbering keyword came back'

    m = smiles('[CH3:1][OH:2]')
    assert explicify_hydrogens(m) == 4
    mapped = {n: m.map_number_of(n) for n in m.atom_numbers}
    assert mapped == {1: 1, 2: 2, 3: 0, 4: 0, 5: 0, 6: 0}


def test_a_new_hydrogen_states_zero_hydrogens_of_its_own():
    """And not `H_UNKNOWN`, which is what `add_atom`'s default stores.

    `[H]C([H])([H])[H]` where each H answers None is a worse record than `C` was.
    """
    m = smiles('C')
    explicify_hydrogens(m)
    assert [m.implicit_h_of(n) for n in m.atom_numbers] == [0, 0, 0, 0, 0]


def test_explicifying_is_information_not_a_repair():
    """`INFO`, not `REPAIRED`: nothing was wrong, both spellings are true of the same compound.

    The one severity difference between the two passes, and deliberate -- implicify runs inside
    `canonicalize()`, where folding an MDL record's spelled-out hydrogens is the repair.
    """
    m = smiles('CCO')
    log = m.log
    explicify_hydrogens(m)
    assert len(log) == 3, 'one record per anchor'
    assert not log.repaired()
    assert {r.rule for r in log} == {'hydrogens:explicify'}


def test_nothing_to_do_is_zero_and_touches_nothing():
    for s in ['[H]C([H])([H])[H]', '[Na+].[Cl-]', 'ClC(Cl)(Cl)Cl', 'O=C=O']:
        m = smiles(s)
        n = len(m)
        log = m.log
        assert explicify_hydrogens(m) == 0, s
        assert len(m) == n
        assert not [r for r in log if r.rule == 'hydrogens:explicify']
