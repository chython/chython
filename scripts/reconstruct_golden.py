#!/usr/bin/env python3
"""How much of a reference-mapped corpus does `reconstruct_mapping()` explain, and does it AGREE?

NOT PART OF THE LIBRARY AND NOT COMMITTED: `scripts/` is temporary by the tree's own rule.  Prints
COUNTS ONLY -- never a structure, a name or an id.

Two numbers, and the second is the one that matters.  COVERAGE is how many records some rung explained;
AGREEMENT is how many of those it mapped the way the corpus says.  A harness that reported only the
first would be measuring the corpus; a harness that compared containers would be measuring nothing at
all, since `__eq__` excludes map numbers and `produced == reference` holds for every record.

    python scripts/reconstruct_golden.py mapping/golden.rdf [records]
"""
import sys as _sys
import os as _os
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))

from collections import Counter
from sys import argv
from time import monotonic

from chython import RDFRead
from chython.reactions import mapping_agrees


def main(path, limit):
    labels = Counter()
    read = skipped = explained = exact = 0
    atoms = Counter()
    started = monotonic()

    with RDFRead(path) as f:
        for i, record in enumerate(f):
            if i >= limit:
                break
            if len(record.products) != 1:
                skipped += 1
                continue
            read += 1
            reference = record.copy()
            try:
                reference.canonicalize()
            except Exception:
                skipped += 1
                read -= 1
                continue
            probe = reference.copy()
            found = probe.reconstruct_mapping()
            if not found:
                continue
            explained += 1
            for label in found:
                labels[label.split(':')[0]] += 1
            agreed, disagreed, missing = mapping_agrees(probe, reference)
            atoms['agreed'] += agreed
            atoms['disagreed'] += disagreed
            atoms['missing'] += missing
            if not disagreed and not missing:
                exact += 1

    elapsed = monotonic() - started
    print('records read              %d  (skipped %d, %.1f rec/s)'
          % (read, skipped, read / max(elapsed, 1e-9)))
    print('explained by some rung    %d  (%.1f%%)' % (explained, 100. * explained / max(read, 1)))
    print('mapping reproduced exactly %d  (%.1f%% of all, %.1f%% of explained)'
          % (exact, 100. * exact / max(read, 1), 100. * exact / max(explained, 1)))
    print()
    print('product atoms: agreed %d, disagreed %d, missing %d'
          % (atoms['agreed'], atoms['disagreed'], atoms['missing']))
    print('  (a symmetric product whose halves came from different inputs disagrees on every atom while')
    print('   being the same answer: automorphic swaps are excused only WITHIN one input.  Inspect, do')
    print('   not read `disagreed` as a defect rate.)')
    print()
    print('rungs that fired:')
    for rung, n in labels.most_common():
        print('  %-16s %6d' % (rung, n))


if __name__ == '__main__':
    main(argv[1] if len(argv) > 1 else 'mapping/golden.rdf',
         int(argv[2]) if len(argv) > 2 else 200)
