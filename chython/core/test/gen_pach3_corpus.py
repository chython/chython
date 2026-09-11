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
"""Freeze the version 3 and version 4 records this tree's writer emits for `pach3_corpus.BUILDERS`.

Run from the repository root:  python -m chython.core.test.gen_pach3_corpus

The output is committed.  Regenerate ONLY for a deliberate layout change, and say in the commit which
byte moved -- a fixture regenerated to make a test green pins nothing at all.
"""
import gzip
import json
from struct import pack

from .pach3_corpus import BUILDERS, V3_PATH, V4_PATH, answers, drawn


def _write(path, version):
    out = bytearray(pack('<I', len(BUILDERS)))
    for name, build in BUILDERS:
        mol = build()
        if version == 3:
            mol = drawn(mol)
        record = mol.pack(compressed=False, version=version)
        blob = json.dumps(answers(mol, version == 3), sort_keys=True).encode()
        payload = name.encode()
        out += pack('<III', len(payload), len(record), len(blob)) + payload + record + blob
    gzip.open(path, 'wb').write(bytes(out))
    print('%s: %d records, %d bytes' % (path.name, len(BUILDERS), len(out)))


if __name__ == '__main__':
    _write(V3_PATH, 3)
    _write(V4_PATH, 4)
