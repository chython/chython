"""Throwaway: size scaling of the distance matrix, and what the batch buffer costs.

Two questions the per-molecule benchmark cannot answer:
  1. chytorch's Floyd-Warshall is O(V^3) and chython's BFS is O(V*E) -- where do they cross?
  2. a [n, n] matrix per molecule, then padded into [B, S, S], is a malloc and a copy per molecule.
     How much of a batch's wall clock is that?
"""
import sys
from os import environ
from time import perf_counter

# chytorch is not a dependency and not in this tree: point CHYTORCH at a checkout, or install it.
if (_chytorch := environ.get('CHYTORCH')):
    sys.path.insert(0, _chytorch)

from numpy import eye, int32, zeros
from chython import smiles
from chytorch.utils.data.molecule._unpack import unpack as ct_unpack

MAX_D = 10
MAX_N = 14


def timeit(fn, arg, repeats=5):
    return min(_t(fn, arg) for _ in range(repeats))


def _t(fn, arg):
    t = perf_counter()
    fn(arg)
    return perf_counter() - t


# --- size scaling ----------------------------------------------------------------------------------

def chain(n):
    """A linear alkane: n atoms, n-1 bonds -- the sparsest connected graph of its size."""
    return smiles('C' * n)


def sweep():
    print('distance matrix, one linear alkane, per molecule:\n')
    print(f'{"atoms":>6}  {"chytorch FW":>12}  {"chython BFS":>12}  {"ratio":>6}')
    for n in (10, 27, 50, 100, 200, 500):
        m = chain(n)
        p2 = m.pack(compressed=False, version=2)
        reps = max(3, 2000 // n)
        t_fw = min(_t(lambda _: [ct_unpack(p2, 0, 0, 1, MAX_N, MAX_D) for _ in range(reps)], None)
                   for _ in range(3)) / reps
        t_bfs = min(_t(lambda _: [m.distance_matrix() for _ in range(reps)], None)
                    for _ in range(3)) / reps
        print(f'{n:>6}  {t_fw * 1e6:9.1f} us  {t_bfs * 1e6:9.1f} us  {t_fw / t_bfs:5.1f}x')


# --- batch assembly --------------------------------------------------------------------------------

def corpus(limit=1024):
    out = []
    with open('pach/lipophilicity.csv') as f:
        next(f)
        for line in f:
            if len(out) >= limit:
                break
            try:
                out.append(smiles(line.rstrip().rsplit(',', 1)[1]))
            except Exception:
                pass
    return out


def batch_from_matrices(mats):
    """What collate does today: per-molecule [n, n] already built, copied into [B, S, S]."""
    s = max(m.shape[0] for m in mats)
    out = eye(s, dtype=int32)[None].repeat(len(mats), 0)
    for i, d in enumerate(mats):
        n = d.shape[0]
        out[i, :n, :n] = d
    return out


def batch_alloc_only(shape):
    b, s = shape
    out = zeros((b, s, s), dtype=int32)
    out[:, range(s), range(s)] = 1
    return out


def batches():
    mols = corpus()
    print('\n\nbatch of 1024 lipophilicity molecules:\n')
    mats = [m.distance_matrix() for m in mols]
    s = max(m.shape[0] for m in mats)
    print(f'padded to S={s}, buffer {1024 * s * s * 4 / 1e6:.1f} MB int32')

    t_d = timeit(lambda ms: [m.distance_matrix() for m in ms], mols)
    t_c = timeit(batch_from_matrices, mats)
    t_a = timeit(batch_alloc_only, (1024, s))
    print(f'  distance_matrix x1024        {t_d * 1e3:7.2f} ms')
    print(f'  pad+copy into [B, S, S]      {t_c * 1e3:7.2f} ms')
    print(f'  of which zeros+diag alloc    {t_a * 1e3:7.2f} ms')
    print(f'  total                        {(t_d + t_c) * 1e3:7.2f} ms'
          f'   -> {1024 / (t_d + t_c) / 1000:.0f} k mol/s, one thread')


if __name__ == '__main__':
    sweep()
    batches()
