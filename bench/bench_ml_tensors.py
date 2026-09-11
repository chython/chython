"""Throwaway: where the time goes on pach bytes -> (atoms, neighbors, distances).

Corpus: pach/lipophilicity.csv (MoleculeNet lipophilicity, ChEMBL ids).
Baseline: chytorch's compiled `_unpack.unpack`, which reads a pach v2 record straight to arrays.
"""
import sys
from os import environ
from statistics import median
from time import perf_counter

# chytorch is not a dependency and not in this tree: point CHYTORCH at a checkout, or install it.
if (_chytorch := environ.get('CHYTORCH')):
    sys.path.insert(0, _chytorch)

from numpy import empty, int32, minimum
from chython import smiles
from chython.core import pach_load
from chytorch.utils.data.molecule._unpack import unpack as ct_unpack

MAX_D = 10
MAX_N = 14


def corpus(limit=2000):
    out = []
    with open('pach/lipophilicity.csv') as f:
        next(f)
        for line in f:
            if len(out) >= limit:
                break
            try:
                m = smiles(line.rstrip().rsplit(',', 1)[1])
            except Exception:
                continue
            out.append(m)
    return out


def timeit(fn, arg, repeats=5):
    best = []
    for _ in range(repeats):
        t = perf_counter()
        fn(arg)
        best.append(perf_counter() - t)
    return min(best)


# --- the candidate paths ---------------------------------------------------------------------------

def path_chytorch_v2(packs):
    for p in packs:
        ct_unpack(p, 0, 0, 1, MAX_N, MAX_D)


def path_chython_container(packs):
    """pach -> container -> per-atom python loop + C distance matrix."""
    for p in packs:
        mol, _ = pach_load(p, compressed=False)
        n = len(mol)
        atoms = empty(n, dtype=int32)
        neighbors = empty(n, dtype=int32)
        for i, a in enumerate(mol.atoms()):
            atoms[i] = a.element + 2
            h = a.implicit_h
            nb = a.degree + (0 if h is None else h)
            neighbors[i] = (nb if nb < MAX_N else MAX_N) + 2
        d = mol.distance_matrix() + 2
        minimum(d, MAX_D + 2, out=d)


def path_decode_only(packs):
    for p in packs:
        pach_load(p, compressed=False)


def path_distance_only(mols):
    for m in mols:
        m.distance_matrix()


def path_atomloop_only(mols):
    for m in mols:
        n = len(m)
        atoms = empty(n, dtype=int32)
        neighbors = empty(n, dtype=int32)
        for i, a in enumerate(m.atoms()):
            atoms[i] = a.element + 2
            h = a.implicit_h
            nb = a.degree + (0 if h is None else h)
            neighbors[i] = (nb if nb < MAX_N else MAX_N) + 2


def main():
    mols = corpus()
    n_atoms = sum(len(m) for m in mols)
    print(f'{len(mols)} molecules, {n_atoms} atoms, '
          f'median {median([len(m) for m in mols]):.0f} atoms/mol\n')

    v2 = [m.pack(compressed=False, version=2, drop=['cip', 'wedges', 'stereo_groups']) for m in mols]
    v4 = [m.pack(compressed=False, version=4, drop=['cip', 'wedges', 'stereo_groups']) for m in mols]
    print(f'pach v2 {sum(len(x) for x in v2) / len(v2):.1f} B/rec, '
          f'v4 {sum(len(x) for x in v4) / len(v4):.1f} B/rec\n')

    rows = [
        ('chytorch _unpack (v2 -> arrays, C)', timeit(path_chytorch_v2, v2)),
        ('chython pach_load v2 + arrays', timeit(path_chython_container, v2)),
        ('chython pach_load v4 + arrays', timeit(path_chython_container, v4)),
        ('  of which: pach_load v2 alone', timeit(path_decode_only, v2)),
        ('  of which: pach_load v4 alone', timeit(path_decode_only, v4)),
        ('  of which: distance_matrix alone', timeit(path_distance_only, mols)),
        ('  of which: python atom loop alone', timeit(path_atomloop_only, mols)),
    ]
    w = max(len(r[0]) for r in rows)
    for label, t in rows:
        print(f'{label:<{w}}  {t * 1e6 / len(mols):8.1f} us/mol   {len(mols) / t / 1000:7.1f} k mol/s')


if __name__ == '__main__':
    main()
