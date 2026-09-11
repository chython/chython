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
"""`mol()` and `rxn()`: one callable per format, both directions.

The direction is the argument's type: a container in means export, a string in means import.  Import
takes no version argument -- the version is sniffed per CTAB, because one RDfile can mix V2000 and
V3000 CTABs (`test/MR.rdf` does).
"""

from ._errors import MalformedCtfile
from ._rxn import emit_rxn, parse_rxn
from ._sdf import RECORD_SEPARATOR, V2000_STAMP, V3000_STAMP, emit_record, parse_record
from ...core import LogRecord, MoleculeContainer
from ...core.reaction import ReactionContainer


__all__ = ['needs_v3000', 'mol', 'rxn']


#: `version=` accepts the two ints people say out loud, and the two stamps the emitters take.
_VERSIONS = {2000: V2000_STAMP, 3000: V3000_STAMP,
             V2000_STAMP: V2000_STAMP, V3000_STAMP: V3000_STAMP}


def _stamp(version):
    try:
        return _VERSIONS[version]
    except (KeyError, TypeError):
        raise ValueError(f'version must be 2000, 3000 or None (auto), got {version!r}')


def needs_v3000(molecule):
    """The reason `molecule` cannot be written as V2000, or `None` if it can.

    A predicate rather than a try/except: an enhanced stereo group does not make `emit_v2000` raise,
    it logs and writes the record anyway, so a dispatcher built on catching the failure would drop
    the AND/OR collections.  A charge is deliberately not on the list -- one past the reach of the
    `ccc` column writes 0 there and the truth in `M  CHG`, which is correct V2000 losing nothing.
    """
    count = molecule.atom_count
    if count > 999:
        return f'{count} atoms do not fit the V2000 3-character count field'
    if molecule.bond_count > 999:
        return f'{molecule.bond_count} bonds do not fit the V2000 3-character count field'
    if molecule.has_coordinates:
        for sid in molecule.atom_numbers:
            atom = molecule.atom(sid)
            # The V2000 coordinate column is 10 characters wide (F10.4), so the negative bound
            # overflows at -10000.0 against 100000.0 positive: measure the formatted width, not the
            # magnitude.
            for value in (atom.x, atom.y):
                if value is not None and len(f'{value:10.4f}') > 10:
                    return f'the coordinate {value} does not fit the V2000 10-character column'
    if molecule.has_stereo_groups and any(kind != 1 for kind, _ in
                                          molecule.canonical_stereo_groups()):
        return 'an AND/OR stereo group has no V2000 spelling'
    return None


def mol(data, *, version=None, meta=None, ignore_stereo=False, log=None):
    """A molecule from a molfile string, or a molfile string from a molecule.

    :param version: export only.  `2000`, `3000`, or `None` (default) to write V2000 unless the
        molecule states something V2000 cannot spell, in which case V3000 with a log line naming the
        reason.  Given explicitly it is obeyed: `version=2000` on a structure that does not fit
        raises rather than switching.
    :param meta: export only.  `None` (default) writes the molecule's own data fields, `{}` writes
        none, a mapping or an iterable of `(name, value)` pairs writes those.  A molfile has no
        data-field section, so they go after `M  END` in SD framing -- what `SDFWrite` writes minus
        the `$$$$` terminator, plus a log line saying so.
    :param ignore_stereo: import only; skip the stereo step, read the constitution either way.
    :param log: a list to append damage reports to.
    """
    log = [] if log is None else log
    if isinstance(data, MoleculeContainer):
        if version is None:
            reason = needs_v3000(data)
            if reason is None:
                stamp = V2000_STAMP
            else:
                stamp = V3000_STAMP
                log.append(LogRecord('ctfile:written-as-v3000', (), f'written as V3000: {reason}'))
        else:
            stamp = _stamp(version)
        lines, _ = emit_record(data, None, meta, version=stamp, separator=False, log=log)
        # Counted off the emitted lines, not off `meta`: the default is the molecule's own, and
        # `emit_record` drops the unparsed bucket on the way out.
        written = sum(1 for x in lines if x.startswith('>  <'))
        if written:
            log.append(LogRecord('ctfile:data-fields-after-end', (),
                                 f'{written} data field(s) written after M  END; mol() writes no '
                                 f'$$$$, so SDFWrite is the call for a record another reader will '
                                 f'walk'))
        return '\n'.join(lines)
    if isinstance(data, ReactionContainer):
        raise TypeError('mol() writes a molecule; use rxn() for a reaction')
    lines = [x.rstrip('\r') for x in data.split('\n') if not x.startswith(RECORD_SEPARATOR)]
    return parse_record(lines, log, ignore_stereo=ignore_stereo)


def rxn(data, *, version=None, ignore_stereo=False, log=None):
    """A reaction from an RXN string, or an RXN string from a reaction.

    :param version: export only; `2000`, `3000` or `None` (auto), as for :func:`mol`.  Auto asks
        :func:`needs_v3000` of every component and escalates if any one answers.

        An agent is not a reason to escalate: V2000's counts line officially carries two fields, but
        an agent count in the third is what the ecosystem writes and reads, so it is expressible and
        `emit_rxn` logs the convention.
    :param ignore_stereo: import only; skip the stereo step, read the constitution either way.
    :param log: a list to append damage reports to.
    """
    log = [] if log is None else log
    if isinstance(data, ReactionContainer):
        if version is None:
            stamp = V2000_STAMP
            for molecule in data.molecules():
                reason = needs_v3000(molecule)
                if reason is not None:
                    stamp = V3000_STAMP
                    log.append(LogRecord('ctfile:written-as-v3000', (), f'written as V3000: {reason}'))
                    break
        else:
            stamp = _stamp(version)
        lines, _ = emit_rxn(data, version=stamp, log=log)
        return '\n'.join(lines)
    if isinstance(data, MoleculeContainer):
        raise TypeError('rxn() writes a reaction; use mol() for a molecule')
    lines = [x.rstrip('\r') for x in data.split('\n')]
    # Locate the first $RXN line.  Input may come from an RDfile paste where $RFMT and $RXN
    # precede the reaction block; slicing from $RXN is the repair, with the skip count logged.
    # Collected in `own` rather than appended straight to `log`, because the reaction that will hold
    # it does not exist yet: it is folded onto `reaction.log` below.
    own = []
    for i, line in enumerate(lines):
        if line.startswith('$RXN'):
            if i > 0:
                own.append(LogRecord('ctfile:rfmt-skipped', (),
                                     f'skipped {i} leading line(s) before $RXN '
                                     f'(e.g. $RFMT header from an RDfile)'))
            lines = lines[i:]
            break
    else:
        raise MalformedCtfile('expected $RXN at the start of the reaction record')
    log.extend(own)
    reaction = parse_rxn(lines, log, ignore_stereo=ignore_stereo)
    reaction.log.absorb('read', own)
    return reaction
