#!/usr/bin/env python3
"""How well does `rxn.attention_mapping()` reproduce a reference mapping, and how fast?

NOT PART OF THE LIBRARY AND NOT COMMITTED: a benchmark script at the repo root is temporary by the
tree's own rule.  `chython/reactions/test/test_attention.py` holds the 25-record subset with a floor;
this runs the full set.  Prints COUNTS AND TIMINGS ONLY -- never a structure or a record.

    python bench_attention_mapping.py [mapping/golden.smiles] [records] [--kekule]

FOUR OUTCOMES PER RECORD, and the difference between the middle two is the reason this script exists:

    exact      every scored product atom took the reference's answer
    partial    some did not
    declined   the mapper refused -- an empty side, or an atom past 14 heavy neighbours
    unreadable the record did not parse, so nothing was measured

`thiele()` FIRST, ON BOTH SIDES, unless `--kekule` is passed.  `mapping_agrees` excuses a disagreement
when the two candidate atoms lie in one automorphism orbit, and a Kekulé ring has no mirror automorphism
-- alternating bond orders distinguish two ortho carbons the aromatic form makes equivalent.  Pass
`--kekule` to see the gap: on the committed 25-record subset it is 12 exact against 23, on identical
mappings.

The score column is the MODEL'S CONFIDENCE AND NOT AN ACCURACY: the mean raw attention at the accepted
cells.  It is printed beside the agreement so the two can be compared, which is the only way to learn
whether a confidence threshold would filter anything useful.
"""
from statistics import mean
from sys import argv
from time import monotonic

from chython import smiles
from chython.reactions import attention_available, mapping_agrees


def main(path, limit, kekule):
    if not attention_available():
        raise SystemExit('needs `chython[mapping]`: pip install onnxruntime chython-rxnmap')

    exact = partial = declined = unreadable = 0
    agreed = disagreed = missing = 0
    scores, exact_scores, partial_scores = [], [], []
    elapsed = 0.

    with open(path) as f:
        for i, line in enumerate(f):
            if i >= limit:
                break
            line = line.strip()
            if not line:
                continue
            try:
                reference = smiles(line.split('\t')[-1])
                if not kekule:
                    reference.thiele()
            except Exception:
                unreadable += 1
                continue

            probe = reference.copy()
            started = monotonic()
            result = probe.attention_mapping()
            elapsed += monotonic() - started

            if result.skipped:
                declined += 1
                continue
            scores.append(result.score)
            a, d, m = mapping_agrees(probe, reference)
            agreed += a
            disagreed += d
            missing += m
            if d or m:
                partial += 1
                partial_scores.append(result.score)
            else:
                exact += 1
                exact_scores.append(result.score)

    read = exact + partial + declined
    atoms = agreed + disagreed + missing
    print('form                      %s' % ('Kekule (as read)' if kekule else 'aromatic (thiele first)'))
    print('records measured          %d  (unreadable %d)' % (read, unreadable))
    print('  exact                   %d  (%.1f%%)' % (exact, 100. * exact / max(read, 1)))
    print('  partial                 %d  (%.1f%%)' % (partial, 100. * partial / max(read, 1)))
    print('  declined                %d  (%.1f%%)' % (declined, 100. * declined / max(read, 1)))
    print()
    print('product atoms scored      %d' % atoms)
    print('  agreed                  %d  (%.2f%%)' % (agreed, 100. * agreed / max(atoms, 1)))
    print('  disagreed               %d  (%.2f%%)' % (disagreed, 100. * disagreed / max(atoms, 1)))
    print('  missing                 %d  (%.2f%%)' % (missing, 100. * missing / max(atoms, 1)))
    print('  (an atom the reference numbers from an input the reference itself does not carry is scored')
    print('   by neither side, so `agreed + disagreed + missing` can fall short of the product atom')
    print('   count.  An incomplete reference is a property of the corpus.)')
    print()
    if scores:
        print('mean model score          %.3f  (exact %.3f, partial %.3f)'
              % (mean(scores), mean(exact_scores) if exact_scores else 0.,
                 mean(partial_scores) if partial_scores else 0.))
    print('time in the mapper        %.1f s  (%.1f ms/record, %.1f rec/s)'
          % (elapsed, 1000. * elapsed / max(read, 1), read / max(elapsed, 1e-9)))


if __name__ == '__main__':
    args = [a for a in argv[1:] if not a.startswith('--')]
    main(args[0] if args else 'mapping/golden.smiles',
         int(args[1]) if len(args) > 1 else 1 << 30,
         '--kekule' in argv)
