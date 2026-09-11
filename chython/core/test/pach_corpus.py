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
"""Reader for the committed pach fixture corpora.

`pach_v2_corpus.bin.gz` and `pach_v0_corpus.bin.gz` are the SPECIFICATION of the legacy format on
this branch, not a convenience.  Every record in them was written by chython 2's own packer and
every answer was produced by chython 2's own unpacker, so where chython 2's writer lost or wrapped
something the loss is in the answers too -- a V3 decoder is held to reproducing what the format
really carried, not to a fidelity it never had.

DO NOT REGENERATE THEM AGAINST A LATER BUILD OF ANYTHING.  A fixture whose answers came from the
code under test asserts that the code agrees with itself.  The generator lives outside the package
and needs a chython 2 interpreter to run; if the corpora are ever rebuilt, they must be rebuilt from
chython 2 and the diff has to be read record by record.

Container format, little-endian throughout:

    <I  record count
    then per record:
    <III  name length, data length, answers length
    name (utf8) | pach record bytes, UNCOMPRESSED | answers (JSON)

The answers of one record:

    'atoms'   [[number, atomic_number, isotope_or_null, charge, radical, implicit_h_or_null,
               stereo (-1 unset / 0 / 1), repr(x), repr(y), degree], ...]
    'bonds'   sorted unique [[low, high, order], ...]
    'ct'      the cis/trans table the record carries: [[n, m, sign], ...]
    'size'    the record's own byte length as its header states it
    'smiles'  present only when the record was built by parsing a string -- the string.  Stereo
              parity is stated in pach against chython 2's neighbour order and in the arena against
              a different one, so a bit alone cannot say which arena parity is right; re-reading the
              same string with the V3 SMILES reader supplies the answer, and both readers number
              atoms in token order, which makes the comparison atom-by-atom sound.

THREE CORPORA, AND THE THIRD ONE IS THE ONE THAT SETTLES VERSION 0.

`pach_v2_corpus.bin.gz`   2492 records written by chython 2.24's packer.
`pach_v0_corpus.bin.gz`   2471 of those SYNTHESISED as version 0 -- the bond order block re-laid out,
                          five 3-bit orders to two bytes instead of eight to three, everything else a
                          byte-for-byte splice -- and then fed back through `_unpack_v0v2.pyx` and
                          kept only where its answers matched the v2 record's.
`pach_v0_native_corpus.bin.gz`
                          236 records written by chython 1.42's packer, which is a REAL v0 WRITER.

Version 0 is chython's own first format: every release from 1.1 through 1.44 wrote a 0 in byte 0 and
the five-orders-to-two-bytes block.  Records in stored data therefore came from a writer that can
still be installed, and the third corpus is that writer's output, decoded by 2.24's decoder.

The two v0 corpora agree byte for byte on 161 of the 236 records they share.  74 of the 75 disagreements
are ONE FIELD: 1.42's SMILES reader leaves an aromatic ring atom's implicit hydrogen count UNKNOWN
where 2.24 computes 1.  Bonds, cis/trans and record length are identical in every one of those, which
is what makes the splice's layout validated rather than merely self-consistent -- the difference is
two readers disagreeing about a molecule, not two writers disagreeing about a format.  The seventy-
fifth is the same kind of difference one field further out: for `c1ccc(/N=N/c2ccccc2)cc1` 1.42 does not
perceive the azo N=N as a cis/trans centre and writes no cis/trans entry at all, where 2.24 writes one,
so those two records differ in LENGTH by four bytes.  Again a perception disagreeing, not a layout.  It also gives
the corpus a real population of the implicit-H sentinel on aromatic atoms, which no 2.24 record has.
"""
import gzip
import json
from pathlib import Path
from struct import unpack_from


__all__ = ['load_corpus', 'V0_PATH', 'V0_NATIVE_PATH', 'V2_PATH']


V2_PATH = Path(__file__).parent / 'pach_v2_corpus.bin.gz'
V0_PATH = Path(__file__).parent / 'pach_v0_corpus.bin.gz'
V0_NATIVE_PATH = Path(__file__).parent / 'pach_v0_native_corpus.bin.gz'


def load_corpus(path):
    """`[(name, record_bytes, answers_dict), ...]` in the order the generator emitted them."""
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
        out.append((name, data, answers))
    if at != len(raw):
        raise ValueError('%s has %d trailing bytes; the corpus is truncated or over-long'
                         % (path.name, len(raw) - at))
    return out
