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
"""Reader for the committed REACTION pach fixture corpus.

`reaction_pach_v2_corpus.bin.gz` is the SPECIFICATION of the reaction-level wire format on this
branch.  It is the sibling of `pach_corpus.py`, one layer up: those corpora pin a MOLECULE record and
this one pins the four bytes chython 2 put in front of a concatenation of them.

PROVENANCE, exactly.  Every record was produced by `gen_reaction_pach_corpus.py`, which runs an
INSTALLED chython 2.24 in a separate interpreter (`oracle.py`) and captures two things per reaction:
the bytes `ReactionContainer.pack(compressed=False)` wrote, and the answers
`ReactionContainer.unpack` gave when handed those same bytes back.  So the answers are chython 2's
reading of chython 2's writing, and nothing in this tree contributed to either.  Nothing was
hand-assembled; nothing was transcribed from a docstring.

DO NOT REGENERATE THIS AGAINST A LATER CHYTHON.  A fixture whose answers came from the code under
test asserts only that the code agrees with itself.  Regenerating needs a chython 2 interpreter and
the diff has to be read record by record.

ONE BEHAVIOUR IS DELIBERATELY ABSENT FROM THE CORPUS, because chython 2 cannot express it.
`ReactionContainer.unpack` slices the product side as `molecules[-products:]`, so a record with ZERO
products comes back with the whole molecule list in the products: the writer emits a `(1, 1, 0, 0)`
header and the reader then reports `CCO>>CCO`.  There is therefore no chython 2 ANSWER for an empty
product side to freeze, so every record here has a non-empty product side and the empty-side
behaviour is asserted directly in `test_reaction_pach.py::test_empty_sides_round_trip` instead.

Container format, little-endian throughout -- the same shape as `pach_corpus.py`, so that a reader of
one recognises the other:

    <I  record count
    then per record:
    <III  name length, data length, answers length
    name (utf8) | reaction pach record bytes, UNCOMPRESSED | answers (JSON)

The answers of one record:

    'smiles'     the reaction SMILES chython 2 was handed
    'counts'     [reactant count, AGENT count, product count] as chython 2's unpacker reported them
    'molecules'  one entry per molecule, in `molecules()` order -- reactants, then agents, then
                 products -- each `{'atoms': [...], 'bonds': [...]}` where an atom is
                 `[number, atomic_number, isotope_or_null, charge, radical, implicit_h, degree]`
                 sorted ascending, and a bond is `[low, high, order]` sorted ascending.  `number` is
                 chython 2's atom number, which for a reaction record IS the atom-to-atom mapping:
                 chython 2 had no separate field for one.
    'atom_counts' [[per reactant], [per agent], [per product]] -- the answer `pack_len` owes.
"""
import gzip
import json
from pathlib import Path
from struct import unpack_from


__all__ = ['load_corpus', 'V2_PATH']


V2_PATH = Path(__file__).parent / 'reaction_pach_v2_corpus.bin.gz'


def load_corpus(path=V2_PATH):
    """`[{'name': str, 'data': bytes, 'answers': dict}, ...]` in the order the generator emitted them.

    Dicts and not tuples, unlike `pach_corpus.load_corpus`: a reaction record's answers have four
    keys and a positional triple at the call site would have to be unpacked into names anyway.
    """
    raw = gzip.open(path, 'rb').read()
    count, = unpack_from('<I', raw, 0)
    at = 4
    out = []
    for _ in range(count):
        name_len, data_len, ans_len = unpack_from('<III', raw, at)
        at += 12
        name = raw[at:at + name_len].decode()
        at += name_len
        data = raw[at:at + data_len]
        at += data_len
        answers = json.loads(raw[at:at + ans_len])
        at += ans_len
        out.append({'name': name, 'data': data, 'answers': answers})
    if at != len(raw):
        raise ValueError('%s has %d trailing bytes; the corpus is truncated or over-long'
                         % (path.name, len(raw) - at))
    return out
