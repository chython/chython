"""Throwaway: the reaction half. `modeling_view()` is pure Python dicts -- how much does it cost?

Corpus is a local file of mapped reaction SMILES; only counts and timings are reported.
"""
from time import perf_counter

from numpy import empty, int32
from chython import smiles


def timeit(fn, arg, repeats=3):
    best = []
    for _ in range(repeats):
        t = perf_counter()
        fn(arg)
        best.append(perf_counter() - t)
    return min(best)


def corpus(limit=500):
    out = []
    with open('mapping/golden.smiles') as f:
        for line in f:
            if len(out) >= limit:
                break
            try:
                r = smiles(line.split()[0])
            except Exception:
                continue
            out.append(r)
    return out


def path_view(rxns):
    for r in rxns:
        r.modeling_view()


def path_view_to_arrays(rxns):
    for r in rxns:
        v = r.modeling_view()
        states = v.states
        n = len(states)
        order = {m: i for i, m in enumerate(states)}
        atoms = empty(n, dtype=int32)
        for i, s in enumerate(states.values()):
            atoms[i] = s[0]
        adj = empty((n, n), dtype=int32)
        adj[:] = 0
        for (a, b) in v.union_bonds:
            i, j = order[a], order[b]
            adj[i, j] = adj[j, i] = 1


def path_parse(lines):
    for line in lines:
        smiles(line)


def main():
    rxns = corpus()
    n_atoms = sum(sum(len(m) for m in r.molecules()) for r in rxns)
    print(f'{len(rxns)} mapped reactions, {n_atoms} atoms total, '
          f'{n_atoms / len(rxns):.0f} atoms/reaction\n')

    lines = [line.split()[0] for line in open('mapping/golden.smiles')][:len(rxns)]
    rows = [
        ('smirks/smiles parse (reference)', timeit(path_parse, lines)),
        ('modeling_view() alone', timeit(path_view, rxns)),
        ('modeling_view() + arrays', timeit(path_view_to_arrays, rxns)),
    ]
    w = max(len(r[0]) for r in rows)
    for label, t in rows:
        print(f'{label:<{w}}  {t * 1e6 / len(rxns):8.1f} us/rxn   {len(rxns) / t / 1000:7.1f} k rxn/s')


if __name__ == '__main__':
    main()
