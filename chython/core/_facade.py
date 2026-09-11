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
"""`smiles()` and `pach()`: one callable per format, both directions.

The direction is the argument's type -- a container in means export, a string or a buffer in means
import -- which is `ctfile`'s `mol()`/`rxn()` and `inchi()` written for the two formats the core owns.
A keyword serving one direction is documented as serving it and ignored by the other, as `mol()`'s
`version=` is: `spec=`, `drop=` and `version=` are export, `log=` is import.

`unpach` is the import half under its own name, `unpack` its chython 2 spelling.  Both dispatch on the
record's version byte rather than trying the molecule door and catching the failure, so a damaged
molecule record is reported as one instead of being retried as a reaction.
"""

from zlib import decompress
from ._core import (MoleculeContainer, QueryContainer, pach_dump, pach_load, read_smiles,
                    write_reaction_smiles, write_smiles)
from .reaction import ReactionContainer, reaction_pach_dump, reaction_pach_load


__all__ = ['smiles', 'pach', 'unpach', 'unpack']


#: Byte 0 of a raw record, per era.  The arena's is `0x33`, the low byte of its little-endian magic.
_MOLECULE_VERSIONS = frozenset((0, 2, 3, 4))
_REACTION_VERSIONS = frozenset((1, 5))
_ARENA_MAGIC = 0x33
_RAW_FIRST_BYTES = _MOLECULE_VERSIONS | _REACTION_VERSIONS | {_ARENA_MAGIC}


def smiles(data, log=None, *, spec=''):
    """A container from a SMILES string, or a SMILES string from a container.

    Import is :func:`~chython.core.read_smiles` and is polymorphic on the arrow: a string with one
    gives a `ReactionContainer`, a string without gives a `MoleculeContainer`.  Export is
    :func:`~chython.core.write_smiles` for a molecule and :func:`~chython.core.write_reaction_smiles`
    for a reaction, which are what `str(x)`, `format(x, spec)` and `x.smiles` call.

    :param log: import only; a list to append the reader's damage reports to, positional because
        `smiles(text, log)` is how a file loop spells it.  Every line lands on the returned
        container's own `log` either way.
    :param spec: export only; the `format()` spec, documented at `docs/io.rst`.

    A query has no export direction -- `read_smarts` is the door in and a `QueryContainer` has no
    SMARTS round trip -- and says so rather than writing a molecule's SMILES for a pattern.
    """
    if isinstance(data, MoleculeContainer):
        return write_smiles(data, spec)
    if isinstance(data, ReactionContainer):
        return write_reaction_smiles(data, spec)
    if isinstance(data, QueryContainer):
        raise TypeError('a query has no SMILES form; print what you want to know about it instead')
    return read_smiles(data, log)


def pach(data, *, log=None, compressed=None, drop=None, version=None):
    """A container from a pach record, or a pach record from a container.

    Import is :func:`unpach`, which reads a molecule record, a reaction record and `to_bytes` output
    alike.  Export is `x.pack()`: :func:`~chython.core.pach_dump` for a molecule and
    :func:`~chython.core.reaction_pach_dump` for a reaction.

    :param compressed: both directions.  On import `None` sniffs, `True` and `False` state it and are
        a `ValueError` when the buffer disagrees.  On export `None` and `True` both compress, which is
        what `pack()` does, and `False` writes the raw record.
    :param drop: export only; the fields whose loss is waived, or `'*'`.  Without it a field the
        format cannot carry is a `ValueError` naming it.
    :param version: export only; `None` for the current layout -- 3 with coordinates and 4 without for
        a molecule, 5 for a reaction -- or a version outright.
    :param log: import only; see :func:`unpach`.

    pach is small and lossy.  `mol.to_bytes()` is the arena verbatim and lossless, and this reads one
    back, so a store may hold both.
    """
    if isinstance(data, MoleculeContainer):
        return pach_dump(data, compressed=compressed is not False, drop=drop, version=version)
    if isinstance(data, ReactionContainer):
        return reaction_pach_dump(data, compressed=compressed is not False, drop=drop,
                                  version=version)
    if isinstance(data, QueryContainer):
        raise TypeError('pach has no record for a query; a pattern is not a stored structure')
    return unpach(data, compressed=compressed, log=log)


def unpach(data, /, *, compressed=None, log=None):
    """A container from a pach record of any version, or from `to_bytes` output.

    NOTHING HAS TO BE DECLARED ABOUT THE BUFFER.  Byte 0 says which era wrote it -- 0, 2, 3 and 4 are
    molecule pach, 1 and 5 reaction pach, `0x33` the arena magic -- and a zlib header's low nibble is
    its compression method, always 8, so none of the six can be one.  `compressed` states it instead
    for a caller who would rather hear that its store is not what it thought.

    `log` chooses between the two error policies.  Without it this is an answer boundary and raises
    `ValueError`, naming what was wrong with the record including the damage that WAS recovered, since
    a caller who cannot have a structure is owed the whole story.  With a list it is the loop-safe
    door: the complaints are appended and the recovered structure returned, or `None` when nothing
    could be built at all -- what a caller walking forty thousand stored records needs, since one bad
    record must not end the loop.

    `data` is positional-only, as chython 2's was.  `unpack` is the same function under that name.
    """
    if isinstance(data, (MoleculeContainer, ReactionContainer, QueryContainer)):
        raise TypeError('unpach reads a record; pach() writes one')
    obj, problems = _load(data, compressed)
    if log is not None:
        log.extend(problems)
        return obj
    if obj is None:
        raise ValueError('this is not a readable pach record: %s' % '; '.join(problems))
    if problems:
        raise ValueError('this pach record is damaged: %s. pass log=[] to take the structure that '
                         'could be recovered from it along with these problems' % '; '.join(problems))
    return obj


#: chython 2's spelling of :func:`unpach`, and the same object -- so `unpack is unpach` and neither is
#: a wrapper that could drift from the other.
unpack = unpach


def _load(data, compressed):
    """`(container or None, problems)` for any pach or arena buffer.

    The version byte is read once, here, and dispatched on.  chython 2 tried the molecule door and
    fell through to the reaction one on `ValueError`, which reads a damaged molecule record as a
    reaction record and reports it as a bad reaction.
    """
    problems = []
    raw = bytes(data)
    if not raw:
        problems.append('the buffer is empty; a pach record is at least a 4 byte header')
        return None, problems
    looks_raw = raw[0] in _RAW_FIRST_BYTES
    if compressed is True and looks_raw:
        problems.append('compressed=True was stated and the buffer begins with %d, which is a raw '
                        'record and not a zlib header' % raw[0])
        return None, problems
    if compressed is False and not looks_raw:
        problems.append('compressed=False was stated and the buffer begins with %d, which is '
                        'neither a pach version nor the arena magic' % raw[0])
        return None, problems
    if not looks_raw:
        try:
            raw = decompress(raw)
        except Exception as err:
            problems.append('the buffer begins with %d, so it is neither a raw record nor a readable '
                            'zlib stream: %s' % (raw[0], err))
            return None, problems
        if not raw:
            problems.append('the buffer decompressed to nothing')
            return None, problems
    if raw[0] in _REACTION_VERSIONS:
        return reaction_pach_load(raw, compressed=False)
    if raw[0] == _ARENA_MAGIC:
        try:
            return MoleculeContainer.from_bytes(raw), problems
        except Exception as err:
            # `from_bytes` is an answer boundary of its own and raises; caught because THIS door has
            # a caller who asked for the complaints, and an arena too short to hold its header is one.
            problems.append('the buffer begins with the arena magic and is not a readable arena: %s'
                            % err)
            return None, problems
    if raw[0] not in _MOLECULE_VERSIONS:
        problems.append('byte 0 is %d, which is neither a pach version -- 0, 1, 2, 3, 4, 5 -- nor '
                        'the arena magic 0x33' % raw[0])
        return None, problems
    return pach_load(raw, compressed=False)
