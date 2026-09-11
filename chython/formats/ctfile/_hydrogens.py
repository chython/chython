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
"""Implicit hydrogens for a molecule built from a CTfile.

The derivation itself is ``chython.core._core.derive_implicit_hydrogen``, shared by every reader in
the tree; this module only ranks it against what a CTfile states, in decreasing authority:

1. a stated count -- the ``MRV_IMPLICIT_H`` data S-group, or MRV's ``hydrogenCount``;
2. the core derivation, the ordinary case;
3. a stated total valence -- V2000 ``vvv``, V3000 ``VAL=``, MRV's ``mrvValence`` -- read only on an
   atom holding no aromatic bond, and there only where it equals the drawn orders (so: no implicit
   hydrogens), or on an atom with no bonds drawn, where it outranks even a derived count.

Rank 3 is last because ``vvv - drawn`` is a hydrogen count only where the drawing is complete; an
atom with no bonds has no drawing to be incomplete, which is the exception.  ``HCOUNT=`` is a query
minimum ("n or more") and is not read.  An undeterminable count is ``H_UNKNOWN``, registered in
``Ctab.unknown_hydrogens`` and logged per atom.  ``HYD_NO_VALENCE_RULE`` is the one core reason code
carrying ``unsupported: `` -- that gap is chython's; every other code describes the file.
"""

from typing import NamedTuple

from ._sgroup import NO_INDEX, SGroup, UNSUPPORTED
from ...core import H_IMPLICIT_MAX, H_UNKNOWN, LogRecord, LOST, REPAIRED
# From `_core` by name: `chython.core` re-exports the container surface and not the derivation or its
# reason codes, which are a format reader's interface.
from ...core._core import (HYD_AMBIGUOUS_AROMATIC, HYD_DERIVED_OTHER_READING, HYD_NO_AROMATIC_FORM,
                           HYD_NO_VALENCE_RULE, HYD_REASON_MASK, derive_implicit_hydrogen)


__all__ = ['apply_mrv_implicit_h', 'calc_implicit', 'HydrogenResult', 'implicit_for_atom',
           'implicit_from_valence', 'implicit_h_records', 'valence_for_write', 'H_MAX',
           'ZERO_VALENCE', 'MRV_IMPLICIT_H', 'StatedChannels', 'MOLFILE_CHANNELS']


# The data S-group that states a hydrogen count, its datum spelled `IMPL_H<n>`.  A spec extension,
# and the only channel in either CTAB version stating a count rather than a query minimum, so it is
# what this library reads and writes when the valence rules would not reproduce the number.
MRV_IMPLICIT_H = 'MRV_IMPLICIT_H'
_IMPL_H_PREFIX = 'IMPL_H'


class StatedChannels(NamedTuple):
    """How the dialect being read spells the two things that answer an unknown hydrogen count.

    Advice has to name a channel the file in front of the reader has: a molfile states a count with an
    ``MRV_IMPLICIT_H`` data S-group and an ``.mrv`` with ``hydrogenCount``.  `count` states a count
    outright; `valence` states a total valence and is ``None`` for a dialect with no such channel
    (CML), which tells the reader not to offer it.
    """
    count: str
    valence: str | None


#: The molfile family's channels, and the default: `Ctab` is the molfile intermediate before it is
#: anything else.
MOLFILE_CHANNELS = StatedChannels(count=f'an {MRV_IMPLICIT_H} data S-group', valence='a `vvv`/`VAL=` '
                                  'total valence')


# The largest implicit hydrogen count this module will store: the core's own bound, imported rather
# than restated.  The nibble holds 15 values and the bound is 14 because the top value is
# `H_UNKNOWN`; admitting 15 would make a file stating fifteen hydrogens and an atom nobody could
# compute the same atom.  15 is reserved on both sides of this reader -- `H_UNKNOWN` in the nibble
# written to, `ZERO_VALENCE` in the `vvv` field read from.
H_MAX = H_IMPLICIT_MAX

# Bond orders that do not contribute to a valence sum.  Order 8 is chython's "special" bond, used
# for donor-acceptor and ionic contacts in metal complexes: it carries no electron pair to count,
# and summing it would make every ferrocene iron a valence error.
_UNCOUNTED_ORDERS = frozenset((8,))

# V2000's spelling of "stated zero valence" in `vvv`, which the spec documents as
# `0=default, 1-14, 15=zero valence`: the field's own 0 means "not stated", so a valence that really
# is nothing needs its own token.
#
# V3000 does NOT spell it this way -- `VAL` is `>0, or -1 for zero` -- so this constant is the
# internal value and each version translates at its own boundary (`_v2000.parse_atom` and `_v3000`
# read their spelling onto `atom.valence = 0`; `emit_v3000` writes `VAL=-1` back).  That translation
# is the only thing standing between a stated zero valence and a V3000 file claiming fifteen.
ZERO_VALENCE = 15


def _drawn_sum(mol, sid):
    """The sum of the bond orders actually drawn at `sid`, with uncounted orders left out.

    An aromatic bond is counted at face value, 4.  A total valence needs a class-dependent subtrahend
    there, so callers that must not do that arithmetic gate on :func:`_holds_aromatic_bond` rather
    than making this function class-aware.
    """
    total = 0
    for m in mol.neighbors_of(sid):
        order = mol.order_of(sid, m)
        if order not in _UNCOUNTED_ORDERS:
            total += order
    return total


def _needs_statement(mol, sid):
    """Does `sid`'s stored hydrogen count have to be stated for a reader to get it back?

    The predicate both write channels share, so a record cannot state one and not the other and then
    have its two answers disagree.  ``False`` for an atom holding ``H_UNKNOWN``: there is no count to
    state and stating one would invent it.
    """
    stored = mol.implicit_h_of(sid)
    if stored is None:
        return False
    # The reason code is dropped on purpose: every outcome that produces a count round-trips, because
    # the reader on the other side runs this same derivation.
    rule, _ = implicit_for_atom(mol, sid)
    return rule != stored


def _holds_aromatic_bond(mol, sid):
    """Is any bond drawn at `sid` an aromatic one?

    The atoms this module may not do valence arithmetic about.  :func:`_drawn_sum` counts an aromatic
    bond at face value, 4, while a total valence needs a class-dependent number -- a pyrrole
    nitrogen's two aromatic bonds contribute 2 and a benzene carbon's 3 -- which takes the valence
    tables and the aromatic classifier, both in ``core``.  So a format module may neither *write* a
    total valence for such an atom (:func:`valence_for_write`) nor *compare against* one on the way
    in, and both directions ask this one predicate rather than two spellings of it.
    """
    return any(mol.order_of(sid, m) == 4 for m in mol.neighbors_of(sid))


def implicit_from_valence(mol, sid, valence):
    """Hydrogen count implied by a stated total `valence`, or ``None`` when it implies none.

    ``None`` for a statement that cannot be one: a valence below what the file itself drew, or one
    implying more hydrogens than an atom can hold.  Both occur in real files.

    Do not ask this about an atom holding an aromatic bond.  It subtracts ``_drawn_sum``, and a total
    valence cannot be compared with a face-value sum containing a 4; :func:`calc_implicit` declines
    the whole channel there.
    """
    h = valence - _drawn_sum(mol, sid)
    if h < 0 or h > H_MAX:
        return None
    return h


def implicit_for_atom(mol, sid):
    """The implicit hydrogen count for one atom, as CTfile asks the question.  ``(h, reason)``.

    `h` is ``None`` when the count cannot be determined, which is a real answer and not a failure, and
    `reason` is one of the core's ``HYD_*`` codes -- the low three bits an outcome, ``HYD_REASON_MASK``
    wide, plus ``HYD_NO_AROMATIC_FORM`` as a flag that may be OR-ed onto any of them.  A caller that
    wants the outcome must mask; :func:`calc_implicit` is the worked example.

    All this adds to the core derivation is ``count_stated=False``, which is a fact about the format: a
    bond block entry of type 4 has no channel for pyrrole-versus-pyridine, so an atom whose class only
    the ring can settle comes back ``None`` and ``kekule()`` settles it.
    """
    return derive_implicit_hydrogen(mol, sid, count_stated=False)


def valence_for_write(mol, sid):
    """The total valence a writer should state for `sid`, or ``None`` when the derivation gives it.

    True for every reader, and the ecosystem's readers honour it -- but on an atom with a bond drawn it
    is not what carries the count back into *this* library, which ranks a stated valence below the
    derivation; :func:`implicit_h_records` is the channel that always survives a round trip, firing on
    the same :func:`_needs_statement` predicate.

    ``None`` for an atom holding an aromatic bond because the number is not computable here, not
    because it does not exist -- reaching it needs the core's valence tables and aromatic classifier,
    see :func:`_holds_aromatic_bond` -- and ``vvv`` 0 is the spec's "no marking", so saying nothing is
    legal and true.  A bond-free atom is the exception: the reader takes a stated valence over the
    derivation there, so this ``vvv`` reads back as the count it was computed from.  ``ZERO_VALENCE``
    is returned for a genuine zero, the field's own 0 meaning "not stated".
    """
    if not _needs_statement(mol, sid):
        return None
    if _holds_aromatic_bond(mol, sid):
        return None
    total = _drawn_sum(mol, sid) + mol.implicit_h_of(sid)
    return total if total else ZERO_VALENCE


def implicit_h_records(mol, store=None):
    """``MRV_IMPLICIT_H`` data S-groups for every atom whose count the rules would not reproduce.

    Returns ``(records, log)``.  The write half of the top-authority read channel, and what makes a
    round trip exact for an atom whose stored count the valence tables disagree with -- an aromatic
    pyrrole nitrogen, a hypervalent sulfur.

    An atom whose count is not known holds ``H_UNKNOWN``, so :func:`_needs_statement` skips it in both
    channels; the loss is still logged, because a downstream reader will re-derive the count and may
    land elsewhere.  An atom holding an aromatic bond is written here and *not* in ``vvv``, that being
    a total valence this module cannot compute (:func:`valence_for_write`).

    Any ``MRV_IMPLICIT_H`` record already in `store` is *replaced* rather than added to, so read-write
    round trips stay idempotent and cannot contradict the stored counts.
    """
    log = []
    records = []
    used = {r.index for r in store.records if r.index != NO_INDEX} if store is not None else set()
    number = 1
    # `SGroup.atoms` holds stable ids on this write path and *position* indices on the read path
    # (`apply_mrv_implicit_h`); one production reader of the field in each sense is a rename away from
    # a silent misread.
    aromatic = 0
    for sid in mol.atom_numbers:
        if not _needs_statement(mol, sid):
            continue
        if _holds_aromatic_bond(mol, sid):
            aromatic += 1
        while number in used:
            number += 1
        used.add(number)
        sg = SGroup('DAT', index=number)
        sg.atoms.append(sid)
        sg.name = MRV_IMPLICIT_H
        sg.fields['FIELDTYPE'] = ['N']
        sg.data.append(f'{_IMPL_H_PREFIX}{mol.implicit_h_of(sid)}'.encode('latin-1'))
        records.append(sg)
    if records:
        # The caveat is earned per record: an atom holding an aromatic bond has no valence field for
        # this channel to agree with, so a molecule holding none must not be told about one.
        if aromatic:
            log.append(LogRecord('mdl-h:h-count-written', (),
                                 f'{len(records)} atom(s) carry a hydrogen count the valence rules do not '
                                 f'reproduce; written as {MRV_IMPLICIT_H} data S-group(s), and in the valence '
                                 f'field for the {len(records) - aromatic} of them whose bonds all have an '
                                 f'integral order -- the total valence of an atom holding an aromatic bond is not '
                                 f'derivable here'))
        else:
            log.append(LogRecord('mdl-h:h-count-written', (),
                                 f'{len(records)} atom(s) carry a hydrogen count the valence rules do not '
                                 f'reproduce; written as {MRV_IMPLICIT_H} data S-group(s) and in the valence '
                                 f'field'))
    if mol.unknown_h_count:
        log.append(LogRecord('mdl-h:unknown-h-count-not-written', (),
                             f'{mol.unknown_h_count} atom(s) have no known implicit hydrogen count; neither a '
                             f'valence nor an {MRV_IMPLICIT_H} group is written for them, leaving the file as '
                             f'silent on the point as the one they came from', LOST))
    return records, log


def apply_mrv_implicit_h(ctab, log):
    """Turn the ``MRV_IMPLICIT_H`` data S-groups in `ctab` into stated hydrogen counts.

    A spec extension, honoured because it is the only channel in either version stating a hydrogen
    count as a *count* rather than as a query minimum.  It is written for exactly the atoms whose
    count cannot be re-derived -- an aromatic pyrrole nitrogen -- so a reader that ignores it gets
    those atoms wrong.

    Lives here rather than in either parser because *both* versions carry these groups.  The groups
    stay in the record: they are what the file says, and the writer re-derives its own set from the
    stored counts (:func:`implicit_h_records`), so nothing accumulates.
    """
    for sg in ctab.sgroups:
        if not sg.is_data() or sg.name.upper() != MRV_IMPLICIT_H:
            continue
        elif sg.name != MRV_IMPLICIT_H:
            # Folded rather than skipped: a `FIELDNAME` is free text and case is meaning in general,
            # but this name is a keyword of an extension, so another spelling states the same thing.
            log.append(LogRecord('mdl-h:mrv-implicit-h-case-folded', (),
                                 f'{MRV_IMPLICIT_H} spelled {sg.name!r}; read as the same field', REPAIRED))
        if len(sg.atoms) != 1 or not sg.data:
            # Malformed S-group (wrong atom/data count): a file defect, not an unmodelled construct.
            log.append(LogRecord('mdl-h:mrv-implicit-h-bad-sgroup', (),
                                 f'{MRV_IMPLICIT_H} on {len(sg.atoms)} atoms with {len(sg.data)} data lines; '
                                 f'expected one of each, ignored', LOST))
            continue
        text = sg.data[0].decode('latin-1').strip()
        # The datum is spelled `IMPL_H<n>`, e.g. `IMPL_H1`.
        digits = text[len(_IMPL_H_PREFIX):].strip() \
            if text.upper().startswith(_IMPL_H_PREFIX) else ''
        if not digits.isdigit():
            log.append(LogRecord('mdl-h:mrv-implicit-h-bad-value', (),
                                 f'unsupported: {MRV_IMPLICIT_H} value {text!r} is not {_IMPL_H_PREFIX}<n>, ignored',
                                 LOST))
            continue
        position = sg.atoms[0]
        if 0 <= position < len(ctab.atoms):
            ctab.atoms[position].stated_h = int(digits)


class HydrogenResult:
    """What :func:`calc_implicit` did: ``.log`` and ``.unknown``.

    An object rather than a tuple, and deliberately not unpackable: the answer has already grown from
    one field to two, and a positional read has to fail now rather than misread later.
    """
    __slots__ = ('log', 'unknown')

    def __init__(self, log, unknown=()):
        self.log = log
        #: stable ids whose implicit hydrogen count could not be determined.  They hold ``H_UNKNOWN``
        #: in the arena, so this is the repair pipeline's input rather than the only record of them.
        self.unknown = unknown

    def __repr__(self):
        return f'HydrogenResult({len(self.log)} log lines, {len(self.unknown)} unknown)'


def calc_implicit(mol, stated=None, valences=None, channels=MOLFILE_CHANNELS):
    """Set the implicit hydrogen count of every atom in `mol`.  Returns a :class:`HydrogenResult`.

    `stated` maps stable id to a count the file gave outright; `valences` maps stable id to a stated
    total valence.  Both are sparse -- an absent atom simply had nothing stated about it.

    An atom whose element the file never named is the marker (:data:`LABEL_ELEMENT`) by the time it
    gets here, so it needs no exemption: element 0 holds no hydrogens and the derivation says zero.

    `channels` is a :class:`StatedChannels` naming how *the dialect being read* spells those two
    inputs, since every "count not known" line ends in advice; defaults to :data:`MOLFILE_CHANNELS`.

    This function does not raise: an atom whose count no source determines is reported in
    ``result.unknown`` and gets ``H_UNKNOWN``.  Computed in full, then written in one edit scope --
    the container's readers need a clean arena, and ``set_hydrogens`` outside a scope would cost one
    buffer copy per atom.
    """
    log = []
    stated = stated or {}
    valences = valences or {}
    counts = {}
    unknown = []
    for sid in mol.atom_numbers:
        # Asked for every atom, including one the record already answered for: the derivation carries
        # an observation about the input (the no-aromatic-form flag) that a stated count does not
        # displace.  The count derived for such an atom is then discarded.
        h, reason = implicit_for_atom(mol, sid)
        # The low three bits are the outcome and `HYD_NO_AROMATIC_FORM` is a flag OR-ed on top of any
        # of them, so switching on the raw value would stop recognising every flagged outcome.
        outcome = reason & HYD_REASON_MASK
        valence = valences.get(sid)
        if valence is not None and _holds_aromatic_bond(mol, sid):
            # The stated-valence channel is scoped to atoms whose bonds all have an integral order: a
            # total valence cannot be compared with `_drawn_sum`'s face-value 4, the correct
            # subtrahend being class-dependent (see `_holds_aromatic_bond`).  Gated at the one point
            # the field is read.  `unsupported: `, because the file made a legal statement this reader
            # does not read.
            log.append(LogRecord('mdl-h:aromatic-valence-not-read', (sid,),
                                 f'{UNSUPPORTED}atom {sid}: stated valence {valence} not read on an atom '
                                 f'holding an aromatic bond, where a total valence cannot be compared with the '
                                 f'bond orders drawn; the count comes from the derivation alone. State it with '
                                 f'{channels.count}, or kekulise before writing', LOST))
            valence = None

        if reason & HYD_NO_AROMATIC_FORM:
            # An observation about the input, reported separately because it can be true beside a
            # perfectly good count -- hence above the stated-count short-circuit below, which must not
            # silence it.  The bond count is part of the claim: "element 6 has no aromatic form" alone
            # would be false and would send a reader looking for a missing table row.
            drawn_aromatic = sum(1 for o in mol.neighbors_of(sid) if mol.order_of(sid, o) == 4)
            log.append(LogRecord('mdl-h:no-aromatic-form', (sid,),
                                 f'atom {sid}: element {mol.element_of(sid)} charge {mol.charge_of(sid)} has no '
                                 f'aromatic form with {drawn_aromatic} aromatic bond(s) drawn on it, so they '
                                 f'were read as those of a saturated atom'))

        if sid in stated:
            # Rank 1, the one source above the derivation.  A count out of range is not a number, so
            # it falls through to be recomputed; anything else is stored as given.
            given = stated[sid]
            if given is None or given < 0 or given > H_MAX:
                log.append(LogRecord('mdl-h:h-count-out-of-range', (sid,),
                                     f'atom {sid}: stated hydrogen count {given} out of range, recomputed',
                                     REPAIRED))
            else:
                counts[sid] = given
                continue

        if outcome == HYD_DERIVED_OTHER_READING:
            # Answered, but not by the reading the aromatic classifier picked: the class the ring
            # implies has no valence row and the other one does.  Worth a line, without
            # `unsupported: `, because `kekule()` will pick the classifier's class and the stored
            # count will then disagree with it.
            log.append(LogRecord('mdl-h:other-reading', (sid,),
                                 f'atom {sid}: the aromatic class implied for element {mol.element_of(sid)} '
                                 f'charge {mol.charge_of(sid)} has no valence; {h} implicit hydrogen(s) read '
                                 f'from the other reading, which kekule() will not pick'))

        if h is not None:
            # A stated valence that disagrees is overridden, but never silently.
            if valence is not None:
                drawn = _drawn_sum(mol, sid)
                from_valence = implicit_from_valence(mol, sid, valence)
                if from_valence is None:
                    # Impossible stated valence (below drawn sum): a file defect, not an unmodelled construct.
                    log.append(LogRecord('mdl-h:valence-below-drawn', (sid,),
                                         f'atom {sid}: stated valence {valence} is below the {drawn} drawn, '
                                         f'which cannot be a total valence; ignored, {h} implicit hydrogens '
                                         f'from the derivation', LOST))
                elif from_valence != h and mol.degree_of(sid):
                    log.append(LogRecord('mdl-h:valence-overridden', (sid,),
                                         f'atom {sid}: stated valence {valence} implies {from_valence} '
                                         f'implicit hydrogens against the {drawn} bond order(s) drawn, but the '
                                         f'derivation gives {h}; using {h}. A file whose drawn bond orders are '
                                         f'incomplete states a valence the drawing does not account for'))
                elif from_valence != h:
                    # An atom with no bonds is the one place a stated valence outranks the
                    # derivation: `vvv` is disqualified elsewhere because it counts orders the file
                    # does not always draw, and there are no undrawn orders here, so the statement
                    # can only be counting hydrogen.  Only this layer can see it: `valence_rules.tsv`
                    # gives every metal and metalloid a free-atom row, which is correct and is also
                    # the first match for "how many hydrogens", with 0.
                    log.append(LogRecord('mdl-h:bond-free-valence', (sid,),
                                         f'atom {sid}: stated valence {valence} with no bonds drawn, so it can '
                                         f'only be {from_valence} implicit hydrogen(s); read as that in '
                                         f'preference to the {h} the derivation gives, there being no undrawn '
                                         f'bond orders for the valence to be counting instead', REPAIRED))
                    h = from_valence
            if h > H_MAX:  # unreachable from the shipped tables; a guard, set_hydrogens would raise
                log.append(LogRecord('mdl-h:h-count-clamped', (sid,),
                                     f'atom {sid}: the derivation implies {h} hydrogens, clamped to {H_MAX}',
                                     REPAIRED))
                h = H_MAX
            counts[sid] = h
            continue

        # The derivation is silent, so the stated valence is the only information there is, and it is
        # taken in one direction only: equal to the drawn sum says "nothing is left over for
        # hydrogen", which reads the same under any valence model, while exceeding it is a claim the
        # difference is hydrogen and depends on the writer having meant the spec's total valence.
        #
        # One outcome reaches the message below, so its prefix is constant: `HYD_AMBIGUOUS_AROMATIC`
        # requires an aromatic bond, and such an atom no longer arrives with a `valence` at all.
        if valence is not None:
            drawn = _drawn_sum(mol, sid)
            if valence == drawn:
                counts[sid] = 0
                continue
            log.append(LogRecord('mdl-h:no-derivable-count-with-valence', (sid,),
                                 f'{UNSUPPORTED}atom {sid}: no derivable hydrogen count for element '
                                 f'{mol.element_of(sid)} charge {mol.charge_of(sid)}, and the stated valence '
                                 f'{valence} exceeds the {drawn} drawn by {valence - drawn}. Not read as '
                                 f'{valence - drawn} hydrogen(s): a valence stated on an atom the tables do '
                                 f'not model is commonly a coordination or oxidation-state marking rather than '
                                 f'a total valence. Count not known; state it with {channels.count}', LOST))
            unknown.append(sid)
            counts[sid] = H_UNKNOWN
            continue

        # Nothing determines the count, so it is UNKNOWN -- a third state, neither a number nor a
        # reason to reject the record.  Two ways to arrive here, and the core's reason code says
        # which: `HYD_NO_VALENCE_RULE` (no row for this element in this state) or
        # `HYD_AMBIGUOUS_AROMATIC` (drawn aromatic, only the ring can settle its class).  Both are
        # stored with a receipt -- the sentinel in the arena, the id in `unknown`, and a log line
        # saying which, since only that decides what repair a caller runs.
        unknown.append(sid)
        if outcome == HYD_AMBIGUOUS_AROMATIC:
            # No prefix: chython models an aromatic pyrrole nitrogen perfectly well, and what happened
            # is that a CTfile has no channel for saying which class the ring means.
            log.append(LogRecord('mdl-h:ambiguous-aromatic', (sid,),
                                 f'atom {sid}: aromatic bond(s) and no Kekule form to derive a hydrogen count '
                                 f'from -- element {mol.element_of(sid)} charge {mol.charge_of(sid)} may either '
                                 f'carry a hydrogen and donate its lone pair or take a ring double bond, and '
                                 f'only the ring decides; count not known. kekule() then read again, or state '
                                 f'the count with {channels.count}', LOST))
        elif outcome == HYD_NO_VALENCE_RULE:
            # `unsupported: `, the one message here that earns it: the valence collection has no row
            # for this element in this state, so the limitation is ours and the file may be fine.
            # The advice half is composed because this is the only line naming both channels, and CML
            # states a hydrogen count and has no total-valence field at all.
            if channels.valence is None:
                lacks = 'the file states no hydrogen count'
                advice = f'State it with {channels.count}'
            else:
                lacks = 'the file states neither a hydrogen count nor a usable valence'
                advice = f'State it with {channels.count} or {channels.valence} to fix the file'
            log.append(LogRecord('mdl-h:no-valence-rule', (sid,),
                                 f'{UNSUPPORTED}atom {sid}: no valence rule for element {mol.element_of(sid)} '
                                 f'charge {mol.charge_of(sid)} radical {bool(mol.radical_of(sid))} with '
                                 f'{_drawn_sum(mol, sid)} drawn bond order(s), and {lacks}; count not known. '
                                 f'{advice}', LOST))
        else:
            # An unrecognised code is not a statement about whose gap it is, so no prefix -- the
            # sibling decision fails open too.  Unreachable with today's four codes.
            log.append(LogRecord('mdl-h:unknown-reason', (sid,),
                                 f'atom {sid}: the derivation was silent about element {mol.element_of(sid)} '
                                 f'charge {mol.charge_of(sid)} and reported reason {reason}, which this reader '
                                 f'has no ruling for; count not known', LOST))
        counts[sid] = H_UNKNOWN

    with mol.edit():
        for sid, h in counts.items():
            mol.set_hydrogens(sid, h)
    if unknown:
        # One summary line as well as the per-atom lines, so a single grep of the log answers "does
        # this record have unknown hydrogens".
        log.append(LogRecord('mdl-h:unknown-h-summary', (),
                             f'{len(unknown)} atom(s) with an unknown implicit hydrogen count: '
                             f'{", ".join(str(x) for x in unknown)}', LOST))
    return HydrogenResult(log, tuple(unknown))
