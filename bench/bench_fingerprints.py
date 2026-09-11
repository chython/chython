"""Per-molecule fingerprint cost against parse cost.  Not source; see CLAUDE.md."""
from statistics import median
from time import perf_counter

from chython.core import read_smiles

# public compounds only, and a spread of sizes: no internal scaffolds anywhere in this tree
COMPOUNDS = [
    'CCO',
    'CC(=O)Nc1ccc(O)cc1',                                     # paracetamol
    'CN1C=NC2=C1C(=O)N(C)C(=O)N2C',                           # caffeine
    'CC(C)Cc1ccc(cc1)C(C)C(=O)O',                             # ibuprofen
    'CC1=C(C(=O)Nc2ccccc2)S(=O)(=O)c2ccccc21',
    'OCC1OC(O)C(O)C(O)C1O',                                   # glucose
    'CC(=O)Oc1ccccc1C(=O)O',                                  # aspirin
    'CN1CCC[C@H]1c1cccnc1',                                   # nicotine
    'C[C@]12CC[C@H]3[C@@H](CC[C@@H]4CC(=O)CC[C@]34C)[C@@H]1CCC2=O',
    'CC(C)(C)NC[C@H](O)c1ccc(O)c(CO)c1',                      # salbutamol
]

REPEATS = 200


def timed(fn, payload):
    best = []
    for _ in range(REPEATS):
        start = perf_counter()
        fn(payload)
        best.append(perf_counter() - start)
    return median(best) * 1e6


def main():
    print(f'{"compound":<12} {"atoms":>5} {"parse":>9} {"morgan":>9} {"linear":>9} '
          f'{"m/parse":>8} {"l/parse":>8}')
    for smi in COMPOUNDS:
        mol = read_smiles(smi)
        parse = timed(read_smiles, smi)
        morgan = timed(lambda m: m.morgan_fingerprint(), mol)
        linear = timed(lambda m: m.linear_fingerprint(), mol)
        print(f'{smi[:12]:<12} {mol.atom_count:>5} {parse:>8.1f}u {morgan:>8.1f}u '
              f'{linear:>8.1f}u {morgan / parse:>8.2f} {linear / parse:>8.2f}')


if __name__ == '__main__':
    main()
