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
"""`Log` and `LogRecord`: the substrate that gives a composed pipeline one record type.

The design constraint these tests exist to hold is that the substrate is **free**.  There are ~268
sites in `chython/` that append a bare sentence and ~917 log assertions in the suite, and none of them
moved to land this.  So the two things most worth pinning are not features: they are that a plain `[]`
-- which is what a reader's `log=` still is -- is unaffected, and that the substring idiom those 917
assertions use keeps working when the entry is a `LogRecord` instead of a `str`.

A pass writes to `container.log` and no pass takes a `log=`, so a pass names its own stage; `stage()`
nesting is what makes that compose.  `test_a_stage_restores_the_previous_origin_and_nests` is the
ratchet, since an inner name must win and the outer must come back.
"""
from pytest import raises
from .. import INFO, LOST, REFUSED, REPAIRED, Log, LogRecord, read_smiles


def test_the_three_argument_construction_that_already_exists_still_works():
    """The five existing `LogRecord(...)` sites pass three positional arguments and must keep working.

    The three new fields are appended and defaulted for exactly this reason, and nothing in the tree
    unpacks a record positionally into three names -- which is what makes appending safe.
    """
    r = LogRecord('groups:13', (1, 2), 'nitro group re-drawn')
    assert r.rule == 'groups:13'
    assert r.atoms == (1, 2)
    assert r.message == 'nitro group re-drawn'
    assert r.severity == INFO
    assert r.stage == ''
    assert r.subject == ''


def test_a_record_answers_substring_containment_from_its_message():
    """`'phrase' in record` is the idiom ~917 assertions use, and it must not see the field tuple.

    Without this the whole substrate is a 917-assertion migration.  With it, an assertion written
    against a bare string passes unchanged when a `Log` wraps that string into a record.
    """
    r = LogRecord('groups:13', (), 'atom 7 carries 3 aromatic bonds; hypercondensed')
    assert 'hypercondensed' in r
    assert 'not in the message' not in r
    # and specifically NOT field membership, which is what a NamedTuple would answer: the rule id is
    # a field VALUE and is nowhere in the sentence, so a plain NamedTuple would say yes here.
    assert 'groups:13' not in r


def test_a_record_prints_as_the_sentence_it_replaced():
    r = LogRecord('r', (), 'the sentence')
    assert str(r) == 'the sentence'
    assert f'{r}' == 'the sentence'
    # repr still shows the whole record -- debugging did not get worse
    assert 'LogRecord' in repr(r) and 'the sentence' in repr(r)


def test_a_bare_string_append_is_wrapped_and_a_record_append_is_not():
    log = Log()
    log.append('a sentence')
    log.append(LogRecord('groups:1', (3,), 'a record', REPAIRED))
    assert len(log) == 2
    assert isinstance(log[0], LogRecord) and log[0].message == 'a sentence'
    assert log[0].severity == INFO
    assert log[1].rule == 'groups:1' and log[1].severity == REPAIRED


def test_a_plain_list_caller_is_untouched():
    """A reader's `log=[]` behaves exactly as a list does.

    Nothing wraps, nothing is stamped, and a reader appending a sentence gets a sentence back.  This
    is why zero of the 268 emit sites had to change.
    """
    log = []
    log.append('a sentence')
    assert log == ['a sentence']
    assert isinstance(log[0], str)


def test_a_stage_stamps_every_record_appended_inside_it():
    log = Log()
    with log.stage('kekule'):
        log.append('charge-separated')
    with log.stage('standardize', rule='groups:13'):
        log.append('nitro re-drawn')
    assert [r.stage for r in log] == ['kekule', 'standardize']
    assert [r.rule for r in log] == ['', 'groups:13']


def test_a_stage_restores_the_previous_origin_and_nests():
    log = Log()
    with log.stage('outer', rule='a', subject='reactant[0]'):
        log.append('one')
        with log.stage('inner'):
            log.append('two')
        log.append('three')
    log.append('four')
    assert [(r.stage, r.rule, r.subject) for r in log] == [
        ('outer', 'a', 'reactant[0]'),
        ('inner', 'a', 'reactant[0]'),
        ('outer', 'a', 'reactant[0]'),
        ('', '', ''),
    ]


def test_a_record_that_names_its_own_rule_keeps_it():
    """Only BLANK provenance is filled in.  A pass that knows its rule id is not overridden."""
    log = Log()
    with log.stage('standardize', rule='stage-level'):
        log.append(LogRecord('groups:77', (1,), 'msg'))
    assert log[0].rule == 'groups:77'
    assert log[0].stage == 'standardize'


def test_severity_is_never_guessed_from_an_incoming_record():
    """`INFO` is a real value, so "unset" and "deliberately INFO" are indistinguishable.

    Filling it in would silently relabel events, so `append` leaves it exactly as given.
    """
    log = Log()
    with log.stage('s'):
        log.append(LogRecord('r', (), 'm'))
    assert log[0].severity == INFO


def test_the_four_severities_filter():
    log = Log()
    log.append(LogRecord('r', (1,), 'did it', INFO))
    log.append(LogRecord('r', (2,), 'differs from the input', REPAIRED))
    log.append(LogRecord('r', (3,), 'could not derive', LOST))
    log.append(LogRecord('r', (4,), 'declined', REFUSED))
    assert [r.atoms for r in log.repaired()] == [(2,)]
    assert [r.atoms for r in log.lost()] == [(3,)]
    assert [r.atoms for r in log.refused()] == [(4,)]
    assert len(log.by_severity(INFO)) == 1


def test_provenance_survives_concatenation():
    """The reason the fields are on the RECORD and not on the `Log` around it.

    Merge two rule tables' logs and a bare index cannot say which table it came from -- `13` names two
    different rules.  Anything held by the container is lost at the `+`, and merging logs is what an
    orchestrator does.
    """
    a, b = Log(), Log()
    with a.stage('kekule'):
        a.append('from kekule')
    with b.stage('standardize'):
        b.append('from standardize')
    merged = list(a) + list(b)
    assert [r.stage for r in merged] == ['kekule', 'standardize']


def test_by_stage_and_by_subject_answer_the_questions_chython_two_could_not():
    log = Log()
    with log.stage('implicify', subject='reactant[0]'):
        log.append('one')
    with log.stage('implicify', subject='product[0]'):
        log.append('two')
    assert len(log.by_stage('implicify')) == 2
    assert [str(r) for r in log.by_subject('product[0]')] == ['two']


def test_a_total_is_computed_from_the_events_and_not_appended_as_a_summary():
    """A record holding the union of the others' atoms makes iterating the log double-count.

    The union is a question the reader asks instead.
    """
    log = Log()
    with log.stage('standardize'):
        log.append(LogRecord('groups:1', (1, 2), 'first'))
        log.append(LogRecord('groups:2', (2, 3), 'second'))
    assert len(log) == 2, 'a summary record would make this 3'
    assert log.atoms_touched() == {1, 2, 3}
    assert log.atoms_touched('standardize') == {1, 2, 3}
    assert log.atoms_touched('kekule') == set()


def test_a_sink_streams_instead_of_accumulating():
    """A log that is only a return value cannot stream."""
    seen = []
    log = Log(sink=seen.append)
    with log.stage('read'):
        log.append('one')
        log.append('two')
    assert len(log) == 0, 'a sink must not also accumulate'
    assert [str(r) for r in seen] == ['one', 'two']
    assert [r.stage for r in seen] == ['read', 'read']


def test_absorb_folds_in_a_log_that_came_back_on_a_result_object():
    """`kekule()`/`thiele()` return their log on a result object as well as writing it here.

    Both carry `.changed` and one more field beside `.log`, which is why the answer is an object; the
    lines land on the container's log through `absorb`, which is the one place that names the stage.
    """
    log = Log()
    log.absorb('kekule', ['charge-separated', 'one hydrogen'], rule='canonicalize:kekule',
               severity=REPAIRED)
    assert len(log) == 2
    assert all(r.stage == 'kekule' and r.rule == 'canonicalize:kekule' for r in log)
    assert all(r.severity == REPAIRED for r in log)


def test_extend_wraps_every_element():
    log = Log()
    with log.stage('read'):
        log.extend(['one', 'two', LogRecord('r', (), 'three')])
    assert len(log) == 3
    assert all(isinstance(r, LogRecord) and r.stage == 'read' for r in log)


def test_record_is_the_shape_a_new_pass_reaches_for():
    log = Log()
    with log.stage('implicify', rule='hydrogens:implicify'):
        log.record('folded', [3, 1, 2], severity=REPAIRED)
    assert log[0] == LogRecord('hydrogens:implicify', (3, 1, 2), 'folded', REPAIRED, 'implicify', '')


def test_a_log_is_a_list():
    """Everything that works on the caller's plain list works here -- that is the whole point."""
    log = Log(['one', 'two'])
    assert len(log) == 2 and bool(log)
    assert isinstance(log, list)
    assert [str(r) for r in log[:1]] == ['one']
    assert 'one' in log[0]


def test_the_constructor_wraps_what_it_is_seeded_with():
    log = Log(['a sentence', LogRecord('r', (), 'a record')])
    assert all(isinstance(r, LogRecord) for r in log)


def test_severity_values_are_distinct_strings():
    """Strings and not an Enum: `_log.py` is pure Python, reached from the extension by a lazy hook,
    and a `LogRecord` gets packed, compared and repr'd -- a string keeps all three trivial."""
    assert len({INFO, REPAIRED, LOST, REFUSED}) == 4
    assert all(isinstance(s, str) for s in (INFO, REPAIRED, LOST, REFUSED))


def test_a_record_is_hashable_and_comparable_as_a_tuple():
    """It is still a NamedTuple: `__contains__` is overridden, equality and hashing are not."""
    a = LogRecord('r', (1,), 'm', REPAIRED, 'stage', 'subject')
    b = LogRecord('r', (1,), 'm', REPAIRED, 'stage', 'subject')
    assert a == b and hash(a) == hash(b)
    assert len({a, b}) == 1
    with raises(AttributeError):
        a.rule = 'x'  # noqa -- immutable, as a record of what happened should be


def test_a_smiles_reader_record_is_findable_by_rule_and_severity():
    """The point of the conversion, on the shortest real emitter: a reader's line is now queryable.

    `c1ccc-c1` promotes to an aromatic system with no Kekule form, so the reader keeps the written
    orders and refuses -- the severity a reader may state without repairing anything.
    """
    log = Log()
    read_smiles('c1ccc-c1', log)
    assert len(log) == 1
    assert log[0].rule.startswith('smiles:')
    assert log.refused() == list(log)


def test_a_reader_states_no_stage_of_its_own():
    """A reader writes into the caller's list, and only the caller knows what to call that block.

    Unlike a pass, which owns `container.log` and names its own stage, a reader is handed a sequence it
    knows nothing about -- so it states its rule and leaves `stage` blank.
    """
    log = Log()
    read_smiles('c-c', log)
    assert len(log) == 2
    assert {x.stage for x in log} == {''}
    with log.stage('read'):
        read_smiles('c-c', log)
    assert [x.stage for x in log] == ['', '', 'read', 'read']


def test_a_plain_list_still_comes_back_from_the_reader():
    """A reader takes either shape and cannot tell: `cdef list` would have refused a subclass, so the
    parameter is `object` and both reach the same `append`."""
    log = []
    read_smiles('c1ccc-c1', log)
    assert len(log) == 1
    assert 'no Kekule form' in log[0]
    assert 'no Kekule form' in str(log[0])
