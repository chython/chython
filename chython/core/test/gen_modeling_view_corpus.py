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
"""Write `modeling_view_corpus.json.gz` from the CURRENT `modeling_view()`.

    python -m chython.core.test.gen_modeling_view_corpus

RUN ONCE, BEFORE THE UNION KERNEL REPLACES THE DICT IMPLEMENTATION.  The committed file is the pin;
re-running it after the rewrite replaces the pin with the thing being pinned.  A byte that differs
afterwards is either a deliberate change to the modelling convention, which is announced in the commit
message, or the defect this fixture exists to catch.
"""
import gzip
import json

from .modeling_view_corpus import PATH, RECORDS


def main():
    out = {}
    for name, rxn in RECORDS.items():
        view = rxn.modeling_view()
        out[name] = {
            'states': [[n, *state] for n, state in view.states.items()],
            'union_bonds': [[n, m, before, after]
                            for (n, m), (before, after) in view.union_bonds.items()],
            'unmapped': dict(view.unmapped),
            'collisions': {side: list(numbers) for side, numbers in view.collisions.items()},
        }
    # mtime=0: regenerating produces identical bytes; text wrapper via GzipFile (Python 3.10
    # gzip.open does not forward mtime= in text mode).
    content = json.dumps(out, indent=1, sort_keys=False, separators=(',', ': ')).encode('utf8')
    with gzip.GzipFile(PATH, 'wb', mtime=0) as gz:
        gz.write(content)
    print(f'wrote {PATH} -- {len(out)} records, '
          f'{sum(len(r["states"]) for r in out.values())} union atoms')


if __name__ == '__main__':
    main()
