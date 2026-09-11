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
"""The valence rule collection: the data file, the compiled tables, and the three questions.

Four kinds of test here and the distinction matters.

  * The data file is what chython 2 states -- `derive()` again and compare.
  * The compiled tables are what the data file says -- the generator's proof, and the reason the
    tables are committed rather than built.
  * The answers reproduce chython 2's, swept exhaustively against an oracle lifted out of the
    molecule methods the rules are embedded in.
  * The answers stated directly, which need no oracle at all. Those are the cases somebody would
    otherwise have to rediscover.

The first and third kinds ask chython 2 through `oracle`: an INSTALLED chython 2 in another
interpreter, one subprocess per sweep, never an import.  Unprovisioned, those tests skip and the other
two kinds still run.

A DISAGREEMENT IS A QUESTION, NOT A VERDICT.  V2's collection contradicts itself in at least one
place (`test_the_phosphorous_acid_row_is_unreachable_in_chython_two`) and V2's reader answers `[O]`
differently (see `test_alternative_spellings.py`).  Neither is smoothed over and neither may be: a
comparison made to pass by changing the core moves the port to match the divergence.
"""
from collections import Counter, defaultdict

from pytest import raises

from chython.core._core import (valence_check, valence_has_rules, valence_implicit_h,
                                valence_rules)

from .gen_valence_rules import (PREAMBLE, PXI, TSV, canonical_order, compile_tables, coverage,
                                derive, hydrogen_tables, read_tsv)


def rows(rules):
    """A rule list as comparable tuples, provenance dropped -- the shape `valence_rules()` gives."""
    out = []
    for rule in rules:
        out.append((rule.z, rule.charge, rule.radical, rule.bonds, rule.h, rule.env))
    return out


# chython 2's three answers, lifted out of the molecule methods the rules are embedded in and run
# IN THE ORACLE INTERPRETER -- see `oracle` for why the second opinion is a subprocess and not an
# import.  One call answers a whole sweep: the round trip costs more than 23,364 lookups do.
#
#   implicit_h   V2's `MoleculeContainer.calc_implicit` -- `valence_rules` for the bond-order sum,
#                then the FIRST rule whose environment is contained in the atom's, and None on either
#                miss.  The rest of that method is the aromatic special case and the element-H short
#                circuit, neither of which is the collection's business: V2 answers 0 hydrogens for
#                any H atom before it ever reaches a rule.
#   accepts      V2's `check_implicit`, the same lift: ANY matching rule wins, not the first.
#   has_rules    where V2 raises ValenceError, which is the boundary between "no rule" and "this
#                element in this state was never described".  The environment is not consulted.
ORACLE = """
from collections import defaultdict

from chython.exceptions import ValenceError
from chython.periodictable.base.element import _elements_map

MAP = _elements_map()


def rules_for(z, charge, radical, order_sum):
    atom = MAP[z](charge=charge, is_radical=radical)
    try:
        return atom.valence_rules(order_sum)
    except ValenceError:
        return None


def counted(environment):
    counts = defaultdict(int)
    for order, number in environment:
        counts[(order, number)] += 1
    return counts


def implicit_h(z, charge, radical, order_sum, environment):
    rules = rules_for(z, charge, radical, order_sum)
    if rules is None:
        return None
    counts = counted(environment)
    for _set, needed, h in rules:
        if all(counts[k] >= c for k, c in needed.items()):
            return h
    return None


def accepts(z, charge, radical, order_sum, environment, want):
    rules = rules_for(z, charge, radical, order_sum)
    if rules is None:
        return False
    counts = counted(environment)
    for _set, needed, h in rules:
        if h == want and all(counts[k] >= c for k, c in needed.items()):
            return True
    return False


out = []
for question in _payload:
    verb = question[0]
    if verb == 'implicit_h':
        z, charge, radical, bonds, env = question[1:]
        out.append(implicit_h(z, charge, radical, bonds, [tuple(t) for t in env]))
    elif verb == 'accepts':
        z, charge, radical, bonds, env, want = question[1:]
        out.append(accepts(z, charge, radical, bonds, [tuple(t) for t in env], want))
    elif verb == 'has_rules':
        z, charge, radical, bonds = question[1:]
        out.append(rules_for(z, charge, radical, bonds) is not None)
    else:
        raise ValueError(verb)
_emit(out)
"""


def oracle(questions):
    """Ask chython 2 a batch of questions. `questions` is a list of `(verb, *args)` tuples."""
    from .oracle import ask

    answers = ask(ORACLE, [list(q) for q in questions])
    assert len(answers) == len(questions), (len(answers), len(questions))
    return answers


# --- the data file

def test_the_shipped_file_is_what_chython_two_states():
    """Re-derive every `common` and `curated` row through the oracle and compare."""
    shipped = read_tsv()
    assert rows(shipped) == rows(canonical_order(derive())), \
        'chython/core/valence_rules.tsv has drifted from chython 2; run gen_valence_rules.py derive'
    assert {r.provenance for r in shipped} == {'common', 'curated'}, \
        'a row with another provenance is fine, but then this test must stop asserting the set'


def test_the_shape_of_the_collection_is_what_was_measured():
    # a size regression is the cheapest way to notice the oracle changed under us, and these are
    # the numbers the generated header states
    shipped = read_tsv()
    provenance = Counter(r.provenance for r in shipped)
    assert len(shipped) == 1036
    assert provenance == {'curated': 780, 'common': 256}
    assert len({r.z for r in shipped}) == 118
    assert len({(r.z, r.charge, r.radical, r.bonds) for r in shipped}) == 576
    assert len({r.env for r in shipped}) == 219
    assert max(len(r.env) for r in shipped) == 7
    assert {r.charge for r in shipped} == set(range(-4, 5))
    assert {r.h for r in shipped} == set(range(5))
    assert {order for r in shipped for order, _ in r.env} == {1, 2, 3}
    # the hot table's extents are these numbers, so a row outside them is a table change
    assert max(r.bonds for r in shipped) == 8
    index, flat, patterns = hydrogen_tables(canonical_order(shipped))
    assert (patterns, len(index), len(flat)) == (65, 2142, 585)


def test_the_shipped_file_is_already_in_scan_order():
    """The file's row order is the order rows are consulted in, so it must be canonical on disk.

    `compile` reorders a file that is not, which means an out-of-order file is a commit that
    skipped the generator -- and the hydrogen count of fifteen keys depends on the order.
    """
    shipped = read_tsv()
    assert rows(shipped) == rows(canonical_order(shipped))


def test_the_shipped_files_comment_block_is_the_generators_preamble():
    """The TSV's header comment is `PREAMBLE`'s output, so it must still BE `PREAMBLE`.

    `write_tsv` emits the constant and `read_tsv` skips every `#` line, which means nothing else in
    this file would notice the two drifting apart -- and the comment block is the only place a
    row-ordering convention is explained to whoever edits the rows.  A reason that can go stale
    silently is worse than no reason, because it will be trusted.

    Drift is possible in one direction only: hand-editing the TSV's comment, or editing the constant
    without re-running the generator.  Both are exactly what someone recording a new convention does.
    """
    assert TSV.read_text(encoding='utf-8').startswith(PREAMBLE), \
        'valence_rules.tsv\'s comment block has drifted from PREAMBLE in gen_valence_rules.py; ' \
        'edit the constant and re-run gen_valence_rules.py derive'


def test_the_committed_tables_are_the_shipped_file():
    """The generated block in _valence.pxi, byte for byte, from the TSV.

    This is what "generated and committed" buys: a wheel build does not depend on chython 2
    importing, and a table produced during a build is invisible in review.
    """
    assert compile_tables(canonical_order(read_tsv())) in PXI.read_text(encoding='utf-8'), \
        'the tables are stale; run gen_valence_rules.py compile and rebuild'


def test_the_compiled_module_is_the_shipped_file():
    """And the extension actually loaded carries those rows, in that order."""
    assert valence_rules() == rows(read_tsv())


# --- the sweeps against chython 2

def test_every_row_reproduces_the_oracle_on_its_own_environment():
    """Each row, asked about exactly the state it describes, in both directions.

    The sweep that catches a wrong offset, a truncated environment or a lost row -- every row is
    reached through its own key, so no row is merely present but unreachable.
    """
    shipped = valence_rules()
    questions = []
    for z, charge, radical, bonds, h, env in shipped:
        questions.append(('implicit_h', z, charge, radical, bonds, env))
        questions.append(('accepts', z, charge, radical, bonds, env, h))
    answers = oracle(questions)
    for i, (z, charge, radical, bonds, h, env) in enumerate(shipped):
        assert valence_implicit_h(z, charge, radical, bonds, env) == answers[2 * i], \
            (z, charge, radical, bonds, env)
        assert (valence_check(z, charge, radical, bonds, h, env) == 'valid') == answers[2 * i + 1], \
            (z, charge, radical, bonds, env)


def test_the_environment_free_sweep_reproduces_the_oracle():
    """Every element, every charge and radical state in the collection, every valence 0..10.

    Environment left empty on purpose: it is the case the `env=*` rows answer and the one every
    caller hits most, and it is where a mis-sorted key array shows up as a wrong element's answer.
    """
    states = [(z, charge, radical, order_sum)
              for z in range(1, 119) for charge in range(-4, 5)
              for radical in (False, True) for order_sum in range(11)]
    answers = oracle([('implicit_h', z, charge, radical, order_sum, ())
                      for z, charge, radical, order_sum in states])
    for state, expected in zip(states, answers):
        assert valence_implicit_h(*state) == expected, state


# --- two tables, one collection

def first_matching_row(z, charge, radical, order_sum, environment):
    """The hydrogen answer read straight off the shipped rows: the FIRST match wins.

    Deliberately not the compiled tables and not chython 2 -- this is the definition the dense
    hot table has to be a projection of, written in six lines so that "one indexed load" has
    something to be checked against that is obviously right.
    """
    have = Counter(environment)
    for row_z, row_charge, row_radical, bonds, h, env in valence_rules():
        if (row_z, row_charge, row_radical, bonds) != (z, charge, radical, order_sum):
            continue
        if all(have[token] >= count for token, count in Counter(env).items()):
            return h
    return None


def test_the_dense_table_is_the_full_scan():
    """The hot artifact is a projection of the collection, over its whole domain and past its edges.

    `val_implicit_h` is one indexed load into a generated table, so the thing that can go wrong is
    not a wrong rule but a wrong index: an off-by-one in the charge bias, a pattern shared between
    two states that only agree on the first few bond counts, a `-2` slot that should have been a
    `-1`.  This sweep is the only reason precomputing is allowed at all.  Note the ranges run one
    step outside the table's extents in every direction -- the short circuit that skips the load has
    to answer "no rule" there and not read someone else's row.
    """
    for z in (1, 5, 6, 7, 8, 15, 16, 17, 26, 33, 34, 50, 53, 78, 92, 118):
        for charge in range(-5, 6):
            for radical in (False, True):
                for order_sum in range(10):
                    assert valence_implicit_h(z, charge, radical, order_sum) == \
                        first_matching_row(z, charge, radical, order_sum, ()), \
                        (z, charge, radical, order_sum)


def test_the_dense_table_is_the_full_scan_with_an_environment_too():
    """And every row's own environment reaches the same answer through the `-2` fall-through."""
    for z, charge, radical, bonds, _, env in valence_rules():
        assert valence_implicit_h(z, charge, radical, bonds, env) == \
            first_matching_row(z, charge, radical, bonds, env), (z, charge, radical, bonds, env)


def test_an_environment_can_still_decide_a_hydrogen_count():
    """The `-2` state, and why the hot table could not be "the 43 rules that mention hydrogens".

    A sulfone is not an exotic input.  Its sulfur has no `env=*` row at six bonds, so the only
    thing that answers 0 rather than "no rule" is a row demanding two double-bonded oxygens -- and
    an MDL reader that got "no rule" here would store an unknown hydrogen count for every sulfone,
    sulfonamide, nitro group and perchlorate in the corpus.  Phosphorus at four bonds is the one
    place an environment produces a NON-zero count.
    """
    assert valence_implicit_h('S', 0, False, 6, ((2, 'O'), (2, 'O'), (1, 'C'), (1, 'C'))) == 0
    assert valence_implicit_h('S', 0, False, 6) is None
    assert valence_implicit_h('P', 0, False, 4, ((1, 'O'), (1, 'O'), (2, 'O'))) == 1
    assert valence_implicit_h('P', 0, False, 4) is None


def test_the_two_tables_are_not_one_table_queried_twice():
    """The hot table answers the hydrogen question ONLY, and merging the two is a bug.

    Sulfur at three bonds is described by the collection -- neutral sulfur is thoroughly described
    -- and no row covers three bonds, so the verdict is a violation while the hydrogen count is
    "no rule".  A single table serving both would have to answer one of those two questions wrong:
    a check that consulted the hot table would call a sulfone unknown, and a hydrogen count that
    consulted the cold table would pay 576 keys of binary search on every carbon in every molecule.
    """
    assert valence_implicit_h('S', 0, False, 3) is None
    assert valence_check('S', 0, False, 3, 0) == 'violation', \
        'the check question must not be answered out of the dense hydrogen table'
    assert valence_check('S', 0, False, 6, 0,
                         ((2, 'O'), (2, 'O'), (1, 'C'), (1, 'C'))) == 'valid'


def test_has_rules_is_the_valence_error_boundary():
    """`valence_has_rules` must be exactly where chython 2 raises, environment ignored.

    Swept over every element rather than a named sample: the oracle answers a whole sweep in one call,
    so there is no reason to sample.
    """
    states = [(z, charge, radical, valence)
              for z in range(1, 119) for charge in range(-2, 3)
              for radical in (False, True) for valence in range(9)]
    answers = oracle([('has_rules', *state) for state in states])
    for state, expected in zip(states, answers):
        assert valence_has_rules(*state) is expected, state


# --- row order

def test_row_order_within_a_key_is_observable():
    """Find every key whose answer depends on the order of its rows, and pin the count.

    Two rows at one key can both match one neighbourhood exactly when the multiset union of what
    they demand still fits inside the key's bond-order sum. Where that happens and their hydrogen
    counts differ, the first row wins and the order is semantics. This searches for those keys
    rather than asserting a belief about them, and the count is pinned so that a change to the
    collection has to look at this test.
    """
    by_key = defaultdict(list)
    for z, charge, radical, bonds, h, env in valence_rules():
        by_key[(z, charge, radical, bonds)].append((h, env))

    observable = []
    for (z, charge, radical, bonds), entries in by_key.items():
        for i in range(len(entries)):
            for j in range(i + 1, len(entries)):
                if entries[i][0] == entries[j][0]:
                    continue
                union = Counter(entries[i][1])
                other = Counter(entries[j][1])
                for token, count in other.items():
                    union[token] = max(union[token], count)
                if sum(order * count for (order, _), count in union.items()) <= bonds:
                    observable.append((z, charge, radical, bonds))
    assert len(observable) == 15, \
        'the number of keys whose hydrogen answer depends on row order changed; if you edited the ' \
        f'collection, decide whether that was intended: {sorted(set(observable))}'
    # eleven of them are bare atoms, where a `common` row and a `curated` row disagree; the other
    # four are bonded states of phosphorus and sulfur where two curated rows overlap
    assert sum(1 for *_, bonds in observable if bonds == 0) == 11


def test_a_bare_atom_takes_its_common_valence_and_not_the_atomic_row():
    # `[C]` in chython 2 is methane's four hydrogens: the common-valence row sits at the same key
    # as an atomic-carbon row with none, and it is scanned first.  This is the largest class of
    # order-observable key and the one a reordering would silently change
    assert valence_implicit_h('C', 0, False, 0) == 4
    assert valence_check('C', 0, False, 0, 0) == 'valid', \
        'zero hydrogens on a bare carbon is still a described state -- the count is not chosen, ' \
        'but it is accepted, which is the whole difference between the two questions'
    assert valence_implicit_h('B', 0, False, 0) == 3
    assert valence_implicit_h('S', 0, False, 0) == 2


def test_a_LONE_metal_is_the_metal_and_only_a_SUBSTITUTED_one_takes_hydrides():
    """Aluminium reverses the row order of the test above, and that is the chemistry, not a slip.

    Nine elements carry both a hydride ladder and a bare-atom row at the same key.  Eight of them --
    B C Ge P S Se Si Te -- put the hydride first, so `[C]` is methane.  Aluminium is the ninth and
    puts the bare atom first, which makes its ladder look self-contradictory: nothing on it is zero
    hydrogens, one carbon on it is two.  ALUMINIUM IS A METAL AND THE OTHER EIGHT ARE NOT.  A lone
    Al in a connectivity file is the metal or its ion -- a counterion, a coordination centre, a
    charge somebody forgot to draw -- while a substituted Al is an organoaluminium, where the
    hydride is ordinary and load-bearing: DIBAL-H is a two-coordinate aluminium with one.

    Ramil's ruling of 2026-09-04, in his words: *"C-Al is hydride.  Al alone is metal."*  And the
    corollary, which is why this costs nothing: *"[AlH3] is valid smiles.  ctfile and other formats
    without hydrogens should treat it as just metal."*  A notation that means alane can say so, and
    is believed; a file that states nothing is stating a metal.  So the eight are not a rule Al
    breaks, they are the non-metal half of one convention.

    Pinned because a reader who finds the inversion by grep -- as I did -- will read it as the one
    ordering defect in 1036 rows and move the row.  The whole ladder is asserted, since moving that
    row is exactly the edit that keeps every other rung passing.
    """
    assert valence_implicit_h('Al', 0, False, 0) == 0, 'a lone aluminium is aluminium'
    assert valence_implicit_h('Al', 0, False, 1, [(1, 'C')]) == 2
    assert valence_implicit_h('Al', 0, False, 2, [(1, 'C')] * 2) == 1, 'DIBAL-H'
    assert valence_implicit_h('Al', 0, False, 3, [(1, 'C')] * 3) == 0

    # ...and zero on the lone atom is a CHOICE among described states, not the absence of a row:
    # three is legal, it is simply not what an unstated count derives to
    assert valence_check('Al', 0, False, 0, 3) == 'valid', 'alane is a described state'
    assert valence_check('Al', 0, False, 0, 0) == 'valid'

    # the eight non-metals, for the contrast the docstring rests on -- all hydride-first
    for element, hydrogens in [('B', 3), ('C', 4), ('Ge', 4), ('P', 3),
                               ('S', 2), ('Se', 2), ('Si', 4), ('Te', 2)]:
        assert valence_implicit_h(element, 0, False, 0) == hydrogens, element


def test_the_phosphorous_acid_row_is_unreachable_in_chython_two():
    """chython 2's own collection contradicts itself here, and the port reproduces it exactly.

    `P` at three bonds with one single and one double oxygen is phosphorous acid, and the curated
    row for it grants two hydrogens. The common valence 3 row sits at the same key with zero and
    is scanned first, so the curated row's two-hydrogen state cannot be reached. Recorded rather
    than smoothed: fixing it means changing the collection on a chemistry argument, in the TSV,
    not quietly reordering a scan.
    """
    acid = [(1, 'O'), (2, 'O')]
    assert valence_implicit_h('P', 0, False, 3, acid) == 0
    # the row is in the collection, and the state it describes is still ACCEPTED -- only the
    # hydrogen count is shadowed, because that question stops at the first match and this one does
    # not
    assert valence_check('P', 0, False, 3, 2, acid) == 'valid'
    assert (15, 0, False, 3, 2, ((1, 8), (2, 8))) in valence_rules()


# --- the three states

def test_the_three_states_are_what_they_say():
    # valid: a row accepts it
    assert valence_check('C', 0, False, 4, 0, [(1, 'C')] * 4) == 'valid'
    # violation: neutral carbon is thoroughly described and there is no row for five bonds
    assert valence_check('C', 0, False, 5, 0, [(1, 'C')] * 5) == 'violation'
    # unknown: nobody wrote anything about carbon at charge -4, so the collection makes no claim
    assert not valence_has_rules('C', -4, False, 0)
    assert valence_check('C', -4, False, 0, 0) == 'unknown'


def test_unknown_is_a_gap_in_the_collection_and_violation_is_a_claim():
    """The boundary is (element, charge, radical), and this is the pair that shows why.

    A described element in a described charge and radical state has a positive list behind it, so
    a state absent from that list is a claim about the molecule. An element the collection has
    never spoken about has no list, so absence says nothing -- and that is the number the coverage
    report counts and a `mined:` row fixes.
    """
    assert valence_check('N', 0, False, 5, 0, [(1, 'C')] * 3 + [(2, 'C')]) == 'violation', \
        'neutral nitrogen is described, so five bonds on it is the collection disagreeing with ' \
        'the molecule, not the collection being silent'
    described = {(z, charge, radical) for z, charge, radical, *_ in valence_rules()}
    assert (7, 0, False) in described
    # a radical lanthanide: nothing in the collection, so no verdict may be invented
    assert (60, 0, True) not in described
    assert valence_check('Nd', 0, True, 3, 0, [(1, 'Cl')] * 3) == 'unknown'


def test_a_verdict_is_never_a_rejection():
    # every state, however broken, comes back as one of three strings -- no exception, nothing to
    # catch, and therefore nothing a reader can be tempted to reject on
    for state in ((6, 0, False, 9, 0), (7, 4, True, 8, 4), (92, -4, True, 0, 0),
                  (1, 0, False, 3, 3), (8, 3, False, 7, 2)):
        assert valence_check(*state) in ('valid', 'violation', 'unknown')


def test_an_impossible_hydrogen_count_is_a_verdict_and_not_an_error():
    # a caller checking a count it read out of a file must get an answer
    assert valence_check('C', 0, False, 0, 5) == 'violation'
    assert valence_check('C', 0, False, 0, 4) == 'valid'
    assert valence_check('C', 0, False, 0, 16) == 'violation'      # outside the nibble entirely


# --- the answers that must survive the oracle's removal

def test_the_answers_everyone_knows():
    assert valence_implicit_h('C', 0, False, 4, [(1, 'C')] * 4) == 0
    assert valence_implicit_h('C', 0, False, 3, [(1, 'C')] * 3) == 1
    assert valence_implicit_h('N', 0, False, 0) == 3          # ammonia
    assert valence_implicit_h('O', 0, False, 1, [(1, 'C')]) == 1      # an alcohol
    assert valence_implicit_h('O', 0, False, 2, [(1, 'C')] * 2) == 0  # an ether
    assert valence_implicit_h('H', 0, False, 1, [(1, 'C')]) == 0
    # and hydrogen never grows a partner: chython 2 excludes it from the implicit-H path by atomic
    # number, so the collection has no row for an unbonded neutral non-radical H at all -- the
    # states it does have are H+, H- and the radical.  A caller wanting V2's molecule-level answer
    # for an H atom must apply V2's short circuit itself
    assert valence_implicit_h('H', 0, False, 0) is None
    assert valence_implicit_h('H', 0, True, 0) == 0
    assert valence_implicit_h('H', 1, False, 0) == 0


def test_an_explicit_hydrogen_is_an_ordinary_neighbour():
    # chython 2 counts explicit H in both the order sum and the environment, and a row written
    # against a hydrogen neighbour is unreachable otherwise
    assert valence_implicit_h('C', 0, False, 4, [(1, 'H')] * 4) == 0
    assert valence_implicit_h('C', 0, False, 2, [(1, 'H')] * 2) == 2


def test_no_rule_is_not_zero():
    # the distinction the return type exists for: a pentavalent carbon is not a carbon with no
    # hydrogens, and the arena cannot store the difference
    assert valence_implicit_h('C', 0, False, 5, [(1, 'C')] * 5) is None
    assert valence_implicit_h('C', 0, False, 4, [(1, 'C')] * 4) == 0
    assert not valence_has_rules('C', 0, False, 5)
    assert valence_has_rules('C', 0, False, 4)


def test_the_environment_is_a_lower_bound_and_not_a_match():
    # one single-bonded carbon is the alcohol row, and it must also cover the ether, the ester and
    # anything else that adds neighbours on top of it
    assert valence_implicit_h('O', 0, False, 1, [(1, 'C')]) == 1
    assert valence_implicit_h('O', 0, False, 1, [(1, 'Si')]) == 1     # a wider row still applies


def test_a_multiplicity_is_not_a_set():
    """A sulfone row demands two double-bonded oxygens and must not fire on a sulfoxide.

    The half of the environment test that chython 2 writes as a count dict on top of a set. With
    the set alone, one oxygen would satisfy `=O =O`.
    """
    assert valence_implicit_h('S', 0, False, 6, [(2, 'O'), (2, 'O'), (1, 'C'), (1, 'C')]) == 0
    assert valence_implicit_h('S', 0, False, 6, [(2, 'O'), (1, 'C'), (1, 'C'), (1, 'C'),
                                                 (1, 'C')]) is None


def test_the_two_valence_models_answer_differently():
    """The one input where the chemistry model and the SMILES notation model MUST disagree.

    The mirror of this test lives in the SMILES writer's suite against `smv_default_h`, which for
    the same atom answers 0 -- a bare `S` with six bonds' worth of neighbours is legal SMILES
    needing no brackets, because sulfur's wide valence set in the Daylight specification includes
    6. Here the answer is "no rule", because no compound with six carbons on one neutral sulfur is
    in the collection.

    Same element, same charge, same valence, and the two cannot be one table: this model's answer
    changes with the environment (dimethyl sulfone is fine, above) and the notation model has
    nowhere to put an environment, because a string's syntax cannot depend on what its atoms are
    bonded to. Merging them would either make the reader reject strings RDKit emits or make MDL
    accept structures chython 2 rejects.
    """
    six_carbons = [(1, 'C')] * 6
    assert valence_implicit_h('S', 0, False, 6, six_carbons) is None, \
        'the chemistry model must answer "no rule" for six carbons on a neutral sulfur; the ' \
        'SMILES notation model answers 0 for the same atom, and that is why there are two tables'
    assert valence_check('S', 0, False, 6, 0, six_carbons) == 'violation'
    # and the sharper half: the KEY exists -- sulfur does reach valence 6 -- so this is not a
    # missing valence, it is a missing environment.  A merged table has no way to say both
    assert valence_has_rules('S', 0, False, 6)

    # the second witness, which needs no environment at all: neutral five-valent nitrogen is
    # spellable bare in SMILES and has no chemistry row in any environment
    assert valence_implicit_h('N', 0, False, 5, [(1, 'C')] * 3 + [(2, 'C')]) is None
    assert not valence_has_rules('N', 0, False, 5)


def test_the_two_questions_are_not_each_others_inverse():
    """A legal hydrogen count that `valence_implicit_h` would not have chosen.

    Why there are two entry points rather than one plus a comparison. Several rows can sit at one
    key; the hydrogen count is the first match's, the verdict accepts any match. The state is
    found from the collection itself rather than asserted from belief.
    """
    assert valence_implicit_h('S', 0, False, 1, [(1, 'C')]) == 1
    assert valence_check('S', 0, False, 1, 1, [(1, 'C')]) == 'valid'
    found = None
    for z, charge, radical, bonds, h, env in valence_rules():
        pick = valence_implicit_h(z, charge, radical, bonds, env)
        if pick is not None and pick != h and \
                valence_check(z, charge, radical, bonds, h, env) == 'valid':
            found = (z, charge, radical, bonds, h, pick)
            break
    assert found is not None, 'no state is legal-but-not-chosen, so the two questions would be one'


def test_a_radical_state_is_a_different_question():
    assert valence_implicit_h('N', 0, True, 0) == 2       # aminyl radical
    assert valence_implicit_h('N', 0, False, 0) == 3
    assert valence_implicit_h('C', 0, True, 3, [(1, 'C')] * 3) == 0


def test_a_charged_state_is_a_different_question():
    assert valence_implicit_h('N', 1, False, 0) == 4      # ammonium
    assert valence_implicit_h('N', -1, False, 0) == 2     # amide anion
    assert valence_implicit_h('O', -1, False, 0) == 1     # hydroxide
    assert valence_implicit_h('O', -1, False, 1, [(1, 'C')]) == 0     # an alkoxide


def test_a_metal_has_legal_valences_and_no_hydrogens():
    # every alkali metal states common valence 0 first, which is chython 2's way of saying "these
    # valences are legal, none of them grants a hydrogen"
    assert valence_implicit_h('Na', 0, False, 0) == 0
    assert valence_implicit_h('Na', 0, False, 1, [(1, 'Cl')]) == 0
    assert valence_implicit_h('Na', 0, False, 2, [(1, 'Cl')] * 2) is None


# --- contract

def test_aromatic_and_dative_orders_are_refused_not_counted():
    # there are no aromatic rows, so a caller must state its policy rather than get a plausible
    # number back
    with raises(ValueError, match='order 4'):
        valence_implicit_h('C', 0, False, 3, [(4, 'C'), (4, 'C'), (1, 'C')])
    with raises(ValueError, match='order 8'):
        valence_implicit_h('N', 0, False, 3, [(8, 'Pd')])
    with raises(ValueError, match='outside 1..3'):
        valence_implicit_h('C', 0, False, 3, [(5, 'C')])


def test_an_unknown_element_is_rejected():
    # an element outside the periodic table is a caller error, unlike an undescribed valence state
    for element in ('Xx', 119):
        with raises(ValueError):
            valence_implicit_h(element, 0, False, 1)


def test_r_has_no_valence_rules():
    """R (element 0) is a defined element with no valence rules; the answer is None, not a raise."""
    assert valence_implicit_h(0, 0, False, 0) is None
    assert valence_implicit_h('R', 0, False, 0) is None
    assert not valence_has_rules(0, 0, False, 0)
    assert valence_check(0, 0, False, 0, 0) == 'unknown'


def test_a_charge_that_cannot_be_packed_answers_rather_than_corrupts():
    """A charge outside the key's field must miss, not carry into the atomic number's bits.

    The one failure mode of a packed key that no test on realistic input would ever see: with the
    guard removed, `charge=+9` on carbon would find nitrogen's rows.
    """
    assert valence_implicit_h('C', 9, False, 0) is None
    assert valence_implicit_h('C', -9, False, 0) is None
    assert not valence_has_rules('C', 9, False, 0)
    assert valence_check('C', 9, False, 0, 0) == 'unknown'
    assert valence_implicit_h('N', 1, False, 0) == 4        # what it would have found instead


def test_the_charge_domain_is_the_rules_span_not_the_storage_range():
    """The arena stores -4..+8; the collection describes -4..+4 and nothing above it.

    The two ranges sit close enough to mislead whoever keys a new row, so this states the difference
    rather than leaving it to be inferred from a DEF.  A charge in +5..+8 is STORABLE AND
    UNDESCRIBED: it must get the no-rule answer, never a confident 0.  Answering 0 for a state
    nobody looked at invents chemistry and erases the gap the coverage report exists to count.
    """
    charges = [charge for _, charge, _, _, _, _ in valence_rules()]
    assert min(charges) == -4 and max(charges) == 4
    assert Counter(charges) == {-4: 1, -3: 10, -2: 28, -1: 60, 0: 791, 1: 69, 2: 30, 3: 39, 4: 8}
    # storable, undescribed, and therefore unanswered -- across the whole band, not just one element
    for z in (6, 7, 8, 16, 26):
        for charge in (5, 6, 7, 8):
            assert valence_implicit_h(z, charge, False, 0) is None, (z, charge)
            assert not valence_has_rules(z, charge, False, 0), (z, charge)
            assert valence_check(z, charge, False, 0, 0) == 'unknown', (z, charge)
    # and the top of the domain IS described, so the guard is not simply refusing everything up
    # there: all eight +4 rows are bare cations, Ti(IV) and the early actinides
    assert valence_has_rules('Ti', 4, False, 0)
    assert [z for z, charge, _, _, _, _ in valence_rules() if charge == 4] == [22, 90, 91, 92, 93,
                                                                               94, 97, 104]


def test_a_symbol_and_a_number_are_the_same_question():
    assert valence_implicit_h('C', 0, False, 2) == valence_implicit_h(6, 0, False, 2)
    assert valence_has_rules('S', 0, False, 6) == valence_has_rules(16, 0, False, 6)


# --- the coverage tool

def test_the_coverage_report_counts_gaps_and_not_violations(tmp_path, capsys):
    """The data-driven-development hook, on molecules whose verdicts are known by hand."""
    corpus = tmp_path / 'corpus.smi'
    corpus.write_text('CCO\nc1ccccc1\n[Fe+]\n[H][N]([H])([H])[H]\n')
    unknown, violation = coverage([str(corpus)])
    out = capsys.readouterr().out
    assert '4 molecules' in out
    assert unknown == {'Fe': 1}, \
        'a monocation of iron is a gap in the collection -- nobody wrote a row for it -- and the ' \
        f'report exists to say so, by element, so a mined row can fill it: {unknown}'
    assert violation == {'N': 1}, \
        'four hydrogens on a neutral nitrogen is the collection disagreeing with the molecule, ' \
        f'which is not a gap and must not be counted as one: {violation}'
    # benzene contributes nothing to either column: an aromatic atom is not undescribed, it is a
    # question this collection does not take
    assert '6 aromatic atoms skipped' in out
    assert 'C' not in unknown and 'C' not in violation


def test_both_readers_measure_the_same_coverage(tmp_path, capsys):
    """The pin on switching the coverage tool's parser from chython 2's reader to the core's.

    A parser swap under a number people quote is worse than a broken comparison, because nothing
    announces it: the corpus changes and the report still prints.  So both readers stay wired up and
    this requires them to agree, atom for atom and element for element, over strings they both
    accept.  The V2 half runs in the oracle interpreter, so this stays live after chython 2 leaves
    the tree.
    """
    corpus = tmp_path / 'corpus.smi'
    # public compounds spanning what the collection is asked about: charged, radical, hypervalent,
    # aromatic, metal, explicit hydrogens, a dative contact
    corpus.write_text('\n'.join((
        'CCO', 'CC(=O)O', 'CC(=O)[O-]', 'c1ccccc1', 'c1ccncc1', 'c1cc[nH]c1', 'Cc1ccccc1',
        'CS(=O)(=O)C', 'C[N+](C)(C)C', 'O=[N+]([O-])c1ccccc1', 'OP(=O)(O)O', 'FC(F)(F)S(=O)(=O)O',
        'CC[Si](C)(C)C', 'B(O)(O)c1ccccc1', 'ClC(Cl)(Cl)Cl', 'BrCCBr', 'ICI',
        '[Fe+]', '[Fe+2]', '[Na+].[Cl-]', '[Cu+2].[O-]S(=O)(=O)[O-]',
        '[H][N]([H])([H])[H]', '[H]O[H]', 'N', 'NO', 'ON=O', 'C#N', '[C-]#[O+]',
        'O=C=O', 'S=C=S', 'CN=[N+]=[N-]', 'C1CC1', 'C1CCCCC1', 'O1CCOCC1',
        'c1ccc2ccccc2c1', 'c1ccc(-c2ccccc2)cc1', 'CC(C)(C)OC(=O)N', 'CSC', 'CS(C)=O',
    )) + '\n')

    core_unknown, core_violation = coverage([str(corpus)], 'core')
    core_out = capsys.readouterr().out
    v2_unknown, v2_violation = coverage([str(corpus)], 'v2')
    v2_out = capsys.readouterr().out

    assert core_unknown == v2_unknown, \
        f'the two readers disagree about which states have no rule: {core_unknown} != {v2_unknown}'
    assert core_violation == v2_violation, \
        f'the two readers disagree about which states violate one: {core_violation} != {v2_violation}'
    # and the same corpus, not merely the same verdicts on a smaller one: a reader that accepted
    # fewer strings could agree on every atom it saw and still be measuring something else
    assert core_out.splitlines()[:2] == v2_out.splitlines()[:2], (core_out, v2_out)
    # the comparison is evidence only if it could have failed: something must be in the columns
    assert core_unknown or core_violation


# --- the harness

def test_the_harness_can_fail():
    # the sweeps above are evidence only if the comparison could have disagreed: make the oracle
    # and the port answer different questions on purpose and check the assertion would fire
    asked, = oracle([('implicit_h', 6, 0, False, 4, ())])
    assert valence_implicit_h('C', 0, False, 3) != asked
    assert valence_rules() != rows(read_tsv())[:-1]
    # and the bridge itself carries a real answer rather than a default: a broken subprocess that
    # returned None for everything would make every sweep above pass vacuously
    assert asked == 0
    assert oracle([('has_rules', 6, 0, False, 4), ('has_rules', 6, 0, False, 5)]) == [True, False]


def test_the_derivation_harness_can_fail():
    # separate from the rest because it needs the oracle to hand over 1036 rows rather than answer
    # a question, and a skip here must not hide the cheaper checks above
    assert rows(read_tsv()) != rows(canonical_order(derive()))[:-1]
