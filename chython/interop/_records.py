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
"""Where an import bridge's records go: onto the container the import produced.

THE DIRECTION DECIDES THE DESTINATION.  An import returns a chython container, so `container.log` owns
its records -- unconditionally, whether or not the caller passed a `log=` list.  An export returns a
foreign object; there is no chython container coming out of it, so its records have nowhere to go but
the caller's list, and the export half of every converter here keeps its `log=` for that reason.

The stage is `'interop'` and does not name the toolkit: every rule id already does
(`'rdkit:note'`, `'cdk:coordinates-not-imported'`), so `log.by_stage('interop')` asks "what did a
toolkit conversion say" and the rule says which toolkit said it.
"""

__all__ = ['STAGE', 'deliver', 'mirror']


STAGE = 'interop'


def deliver(container, records, log=None):
    """Put `records` on `container.log` and copy them to the caller's `log` list when there is one.

    PER CONTAINER, NEVER POOLED: a converter reading three molecules calls this three times, so the
    third molecule's records are on the third molecule and its atom numbers are read against it.

    The caller's copy is taken back off the container's log rather than from `records`, so both hold
    the same stamped records instead of two spellings of one event.  `if log is not None` guards only
    that copy; the write to the container is not conditional on anything.
    """
    start = len(container.log)
    container.log.absorb(STAGE, records)
    if log is not None:
        log.extend(container.log[start:])


def mirror(reaction, subject, records):
    """Copy one component's records onto the reaction's log, stamped with which component they name.

    The component keeps its own -- `rxn.log.by_subject('products[0]')` and `rxn.products[0].log` answer
    the same question from the two ends, which they can only do if both hold the records.  `subject` is
    what makes the copy readable: `LogRecord.atoms` are stable ids in ONE container.  `stage` is left
    alone; the import already named it.
    """
    with reaction.log.stage('', subject=subject) as log:
        log.extend(records)
