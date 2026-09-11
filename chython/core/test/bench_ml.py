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
"""Timings for the ML views.  Hand-run; not collected by pytest.

    python -m chython.core.test.bench_ml

THE CORPUS IS GENERATED AND NOT READ FROM A PATH.  A benchmark whose input is a local file measures
nothing on another machine, and the numbers in `docs/ml.rst` are quoted with this module beside them.
Public compounds, repeated to a workable count: the point is the per-structure cost of each path, and
the size distribution is stated rather than sampled.
"""
from statistics import median
from time import perf_counter

from chython.core import TensorEncoding, read_reaction_smiles, read_smiles, unpach


MOLECULES = [
    'CC(=O)Oc1ccccc1C(=O)O',                              # aspirin, 13 atoms
    'CN1C=NC2=C1C(=O)N(C)C(=O)N2C',                       # caffeine, 14
    'CC(C)Cc1ccc(cc1)C(C)C(=O)O',                         # ibuprofen, 15
    'OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O',             # glucose, 12
    'CC(N)C(=O)NC(C)C(=O)NC(C)C(=O)NC(C)C(=O)O',          # a tetraalanine peptide, 21
    'Clc1ccc(cc1)C(c1ccccc1)n1ccnc1',                     # clotrimazole, 22
    'CCOC(=O)c1ccc(N)cc1',                                # benzocaine, 12
    'C1CC2CCC1C2',                                        # norbornane, 7
]

REACTIONS = [
    '[CH3:1][Br:2].[OH-:3]>>[CH3:1][OH:3].[Br-:2]',
    '[CH3:1][C:2](=[O:3])[OH:4].[CH3:5][OH:6]>>[CH3:1][C:2](=[O:3])[O:4][CH3:5].[OH2:6]',
    '[CH2:1]=[CH:2][CH:3]=[CH2:4].[CH2:5]=[CH2:6]>>[CH2:1]1[CH:2]=[CH:3][CH2:4][CH2:5][CH2:6]1',
    '[cH:1]1[cH:2][cH:3][cH:4][cH:5][cH:6]1.[N+:7](=[O:8])([O-:9])[OH:10]'
    '>>[c:1]1([N+:7](=[O:8])[O-:9])[cH:2][cH:3][cH:4][cH:5][cH:6]1.[OH2:10]',
]

REPEATS = 5


def timeit(fn, arg):
    best = None
    for _ in range(REPEATS):
        start = perf_counter()
        fn(arg)
        elapsed = perf_counter() - start
        if best is None or elapsed < best:
            best = elapsed
    return best


def report(rows, unit, count):
    width = max(len(label) for label, _ in rows)
    for label, seconds in rows:
        per = seconds * 1e6 / count
        print(f'  {label:<{width}}  {per:8.2f} us/{unit}   {count / seconds / 1000:8.1f} k {unit}/s')


def main():
    mols = [read_smiles(line) for line in MOLECULES] * 250
    rxns = [read_reaction_smiles(line) for line in REACTIONS] * 125
    sizes = sorted(len(m) for m in mols)
    print(f'{len(mols)} molecules, median {median(sizes):.0f} atoms, max {sizes[-1]} atoms')

    plain = TensorEncoding()
    chytorch_like = TensorEncoding(element_shift=2, neighbor_shift=2, distance_shift=2,
                                   disconnected=1, max_distance=10)
    padded = TensorEncoding(width=64, pad=0, pad_diagonal=1)
    vocabulary = {}
    for mol in mols[:len(MOLECULES)]:
        view = mol.state_view()
        for z, h, n in zip(view.elements, view.hydrogens, view.neighbors):
            vocabulary.setdefault((int(z), int(h), int(n), int(h), int(n)), len(vocabulary) + 1)
    tokenizing = TensorEncoding(vocabulary=vocabulary, unknown=999)

    report([
        ('distance_matrix() alone', timeit(lambda ms: [m.distance_matrix() for m in ms], mols)),
        ('state_view(), physical', timeit(lambda ms: [m.state_view(plain) for m in ms], mols)),
        ('state_view(), shifted and clamped',
         timeit(lambda ms: [m.state_view(chytorch_like) for m in ms], mols)),
        ('state_view(), padded to 64',
         timeit(lambda ms: [m.state_view(padded) for m in ms], mols)),
        ('state_view(), with a vocabulary',
         timeit(lambda ms: [m.state_view(tokenizing) for m in ms], mols)),
        ('transition_view(), molecule',
         timeit(lambda ms: [m.transition_view(plain) for m in ms], mols)),
    ], 'mol', len(mols))

    n_atoms = sum(sum(len(m) for m in r.molecules()) for r in rxns)
    print(f'\n{len(rxns)} mapped reactions, {n_atoms / len(rxns):.0f} atoms/reaction')
    report([
        ('transition_view(), reaction', timeit(lambda rs: [r.transition_view() for r in rs], rxns)),
        ('modeling_view(), dicts over it',
         timeit(lambda rs: [r.modeling_view() for r in rs], rxns)),
    ], 'rxn', len(rxns))

    print('\npach record to arrays, by wire version:')
    for version in (2, 3, 4):
        packed = [m.pack(compressed=False, version=version) for m in mols]
        median_bytes = median(sorted(len(p) for p in packed))
        report([
            (f'unpach(v{version}) + state_view(), median {median_bytes:.0f} B',
             timeit(lambda ps: [unpach(p, compressed=False).state_view(plain) for p in ps], packed)),
        ], 'mol', len(mols))


if __name__ == '__main__':
    main()
