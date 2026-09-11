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
"""`LogRecord` -- one thing a pass did or declined to do, and the ONE definition of it.  Plus `Log`,
the container behind `molecule.log` and `reaction.log`, which gives a composed pipeline one record type.

WHY `Log` EXISTS.  The pipeline the input posture requires -- `read_smiles(log=log)`, `kekule()`,
`standardize()`, `thiele()` -- has ONE destination for everything it reports.  With more than one record
mechanic the caller ends up holding a list of bare `str` mixed with records plus further lists hanging
off result objects: it cannot filter by rule, cannot ask which atoms an event touched (the atom numbers
are prose inside the sentence), and needs `isinstance` to read its own log.  `canonicalize()` is exactly
the function that has to merge all of them.  Every emit site in `chython/` answers with a `LogRecord`,
and `chython/test/test_log_records.py` keeps them that way.

THE SENTENCES ARE AN ASSET, not debt: each names the atom, what was found and what was done instead.
The structure sits around them and none of them is rewritten -- the rule, the atoms and the severity
live in fields, and the prose is left alone.

This lives in `core` and not beside `standardize()`, the first pass to report in this shape, because it
is not that pass's alone: the SMIRKS patcher lives in `core` and reports the parities it dropped in the
same shape (the design's N7), and `core` cannot import `chemistry`.  The alternative was a second
NamedTuple with the same three fields, which is a second TYPE -- `isinstance` would answer no, a caller
merging two logs would get two record classes out of one list, and the two would drift.

`chemistry` re-exports this name, so `from chython import LogRecord` and
`from chython.chemistry import LogRecord` are the same class.

This module is pure Python and imports nothing from the extension on purpose: `_core.pyx` reaches it
through a lazy import, and a dependency the other way would be a cycle.
"""
from contextlib import contextmanager
from typing import NamedTuple


__all__ = ['LogRecord', 'Log', 'recording', 'INFO', 'REPAIRED', 'LOST', 'REFUSED']


# Severity.  FOUR VALUES, READ OFF THE EXISTING MESSAGES rather than invented: every one of them
# says one of these four things, and each answers a different question a caller actually has.
INFO = 'info'
"""Did what was asked.  `'kekulized'`, a rule that applied cleanly."""
REPAIRED = 'repaired'
"""The input was wrong and the answer DIFFERS FROM WHAT THE INPUT SAID.

The load-bearing one, and the reason severity is worth a field.  The whole input posture is store the
garbage, log it, repair only when asked -- and until this field existed a caller had no way to ask
*did anything get repaired?* except by reading message strings by hand.  This is what makes the
posture auditable.  `'read as ...'`, `'clamped to ...'`, `'recomputed'`.
"""
LOST = 'lost'
"""A fact could not be derived and is now unknown.  Never a silent zero, never a guess, never a
raise: `'no Kekule form to derive a hydrogen count'`, `'unknown implicit hydrogen count'`."""
REFUSED = 'refused'
"""Declined to act; the molecule is UNCHANGED where a rule matched.

Always a record and never an exception.  An `ignore=` flag would let the CALLER decide whether a
refused patch is a string in a list or a raised error, so the flag rather than the event would
determine what a defect is and no field would say which had happened.  There is no `ignore=` parameter
anywhere in this design: refusals raise at the answer boundary, and a repair pass is not one.
"""


class LogRecord(NamedTuple):
    """One thing a pass did or declined to do.

    `rule` is the QUALIFIED identity of whatever produced the record, never a bare index: two rule
    tables' logs get concatenated, and an untagged `13` then names two different rules.  A table-driven
    pass writes `'groups:13'`; the SMIRKS patcher, which has no table until the corpus lands, writes
    `'smirks:<the template string>'` -- the template is its own identity.

    `atoms` are stable ids in the container the record is about, and `message` is one sentence a
    human reads.

    THE LAST THREE FIELDS ARE PROVENANCE, AND THEY ARE ON THE RECORD RATHER THAN ON THE `Log` AROUND
    IT FOR ONE REASON: a log gets concatenated, sliced and merged, and anything held by the container
    is lost the moment it is, and a merged list that can no longer say which table an index came from
    cannot be repaired afterwards.  Six fields is wide for something emitted a million times; losing
    the answer on `log_a + log_b` is worse.

    `atoms` defaults to `()` because a reader reporting a file-level observation genuinely has no
    atoms, and `-1`-as-"no rule" must not come back wearing a different hat.

    `subject` is a STRING and never a container reference: a log outlives the molecules it describes,
    and a record holding one strong reference turns a million-record log into a leak.
    """
    rule: str
    atoms: tuple[int, ...] = ()
    message: str = ''
    severity: str = INFO
    stage: str = ''
    subject: str = ''

    def __contains__(self, item):
        """`'some phrase' in record` asks the MESSAGE, not the field tuple.

        The dominant test idiom over the emit sites is substring containment -- `assert 'cannot carry
        a double bond' in log[0]` -- and with this, that assertion passes whether `log[0]` is a `str`
        or a `LogRecord`.  The log-reading assertions in `chython/*/test/` carry no `str()` for that
        reason.

        It is a deliberate surprise: `in` on a NamedTuple normally means field membership.  The
        mitigation is that field membership over `(rule, atoms, message, severity, stage, subject)` --
        "is this string one of those six values" -- is a question nobody has ever wanted to ask, while
        "does the sentence mention this" is asked in the hundreds.
        """
        return item in self.message

    def __str__(self):
        """The sentence, so `f'{record}'` and `print(record)` read like the string it replaced."""
        return self.message


class Log(list):
    """A list of `LogRecord`s that a bare-string `append` still works on.

    ONE PER CONTAINER, AND IT IS THE DESTINATION.  `molecule.log` and `reaction.log` are `Log`s; no
    pass takes a `log=`, so there is no second sequence for a record to go to and no arrangement in
    which a repair happens and nothing records it.  A reader still takes a `log=` list -- a parse has
    no container to write to until it produces one -- and its records are folded in on seal.

    A PASS NAMES ITS OWN STAGE, with `recording(container, stage=...)`, because there is no
    orchestrator handing it a log to be stamped by.  `stage()` nests and restores, so an inner pass's
    name wins inside its block and it is the more precise of the two: `canonicalize()` opens no stage
    of its own for the passes it runs, and `rxn.standardize()` adds `subject` while leaving the
    molecule pass's `stage` alone.  Only BLANK provenance is filled, so nothing overwrites it.

    Subclasses `list` on purpose: `append`, `extend`, `__len__`, indexing, slicing, iteration and
    truthiness are all inherited, so anything that works on the caller's plain list works on this.
    """
    __slots__ = ('_rule', '_stage', '_subject', '_sink')

    def __init__(self, records=(), *, sink=None):
        """`sink`, when given, is called with each record INSTEAD of storing it.

        A log that is only ever a return value has to be held in full until the work finishes; a
        million-record SDF pass streams through here instead.
        """
        super().__init__()
        self._rule = ''
        self._stage = ''
        self._subject = ''
        self._sink = sink
        for r in records:
            self.append(r)

    def append(self, record):
        """Store one record, wrapping a bare sentence and filling in blank provenance."""
        if isinstance(record, LogRecord):
            # Only BLANK fields are filled: a pass that named its own rule keeps it.  `severity` is
            # not touched at all -- its default is a real value, so "unset" is indistinguishable from
            # "deliberately INFO", and guessing between them would silently relabel events.
            if not record.rule and self._rule:
                record = record._replace(rule=self._rule)
            if not record.stage:
                record = record._replace(stage=self._stage)
            if not record.subject and self._subject:
                record = record._replace(subject=self._subject)
        else:
            record = LogRecord(self._rule, (), str(record), INFO, self._stage, self._subject)
        if self._sink is not None:
            self._sink(record)
        else:
            super().append(record)

    def extend(self, records):
        for r in records:
            self.append(r)

    def record(self, message, atoms=(), *, severity=INFO, rule=''):
        """Emit one record in the current stage.  The shape a new pass should reach for."""
        self.append(LogRecord(rule or self._rule, tuple(atoms), message, severity,
                              self._stage, self._subject))

    @contextmanager
    def stage(self, name, *, rule=None, subject=None):
        """Everything appended inside the block is stamped with this origin.  Nests and restores."""
        saved = (self._rule, self._stage, self._subject)
        self._stage = name
        if rule is not None:
            self._rule = rule
        if subject is not None:
            self._subject = subject
        try:
            yield self
        finally:
            self._rule, self._stage, self._subject = saved

    def absorb(self, stage, lines, *, rule='', subject='', severity=INFO):
        """Fold in a log that came back ON A RESULT OBJECT, stamping the stage it belongs to.

        `kekule()` and `thiele()` return `KekuleResult` / `ThieleResult`, whose `.log` is a list of
        sentences, and the result object stays: each carries `.changed` and (`.unresolved` /
        `.refused`) beside `.log`, which is why the answer is an object and not a list.  The pass folds
        the same lines onto the container's log through here, unconditionally -- the result is a
        convenience, `mol.log` is the storage, and a branch on whether to absorb is where "sometimes we
        record" comes back in.
        """
        with self.stage(stage, rule=rule, subject=subject):
            for line in lines:
                if isinstance(line, LogRecord):
                    self.append(line)
                else:
                    self.record(str(line), severity=severity)

    # -- reading it back.  These are the deliverable; the fields exist to make them possible.

    def by_severity(self, severity):
        return [r for r in self if r.severity == severity]

    def repaired(self):
        """Every record where the answer differs from what the input said.  See `REPAIRED`."""
        return self.by_severity(REPAIRED)

    def lost(self):
        return self.by_severity(LOST)

    def refused(self):
        return self.by_severity(REFUSED)

    def by_stage(self, name):
        return [r for r in self if r.stage == name]

    def by_subject(self, name):
        """Which molecule of a reaction a record is about."""
        return [r for r in self if r.subject == name]

    def atoms_touched(self, stage=None):
        """The union of the atoms the EVENTS name.

        THERE ARE NO SUMMARY RECORDS: a record carrying the union of the atoms other records already
        name double-counts on iteration, so a total is computed by the reader from the events.
        """
        return {a for r in self if stage is None or r.stage == stage for a in r.atoms}

    def __repr__(self):
        return f'Log({list(self)!r})'


@contextmanager
def recording(container, *, stage='', subject=None):
    """The container's own log is where a pass writes.  Yields it, stage stamped.

    THE HELPER HAS NO OFF SWITCH, AND THAT IS THE ENTIRE POINT.  A pass taking `log=` and appending to a
    caller-supplied sequence skips the work when nobody supplies one, so `mol.canonicalize()` with no
    arguments repairs the molecule and leaves no trace of it anywhere.  There is no `log=` parameter on a
    pass: the container is the destination, `mol.log` and `rxn.log` are how a caller reads it, and
    `if log is not None` on a recording path is a defect.

    `subject` names which molecule of a reaction the block is about; see `Log.stage`.  A pass never
    passes either field -- the orchestrator that knows the provenance does.
    """
    with container.log.stage(stage, subject=subject) as log:
        yield log
