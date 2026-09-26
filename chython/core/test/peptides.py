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
"""Public peptides built from sequence, shared by the segmentation, layout and drawing tests.

`peptide()` writes a SMILES: one-letter residues from `SIDE`, `Me` prefixes N-methylate, `G`/`Sar` carry no
side chain, `P` closes its pyrrolidine, `C` pairs by the `ss` index pairs, `cyclic` closes head to tail.
"""
from chython.core import smiles


__all__ = ['CASES', 'PEPTIDES', 'explicit', 'molecule', 'peptide']

SIDE = {'A': 'C', 'V': 'C(C)C', 'L': 'CC(C)C', 'I': 'C(C)CC', 'F': 'Cc1ccccc1', 'Y': 'Cc1ccc(O)cc1',
        'W': 'Cc1c[nH]c2ccccc12', 'S': 'CO', 'T': 'C(C)O', 'M': 'CCSC', 'N': 'CC(N)=O', 'Q': 'CCC(N)=O',
        'D': 'CC(=O)O', 'E': 'CCC(=O)O', 'K': 'CCCCN', 'R': 'CCCNC(=N)N', 'H': 'Cc1c[nH]cn1', 'O': 'CCCN',
        'Abu': 'CC', 'Bmt': 'C(O)C(C)CC=CC'}


def peptide(seq, cyclic=False, ss=(), nterm='', cterm='O', special=None):
    special = special or {}
    closures = iter(range(10, 99))
    tags = {}
    for i, j in ss:
        t = next(closures)
        tags.setdefault(i, []).append(t)
        tags.setdefault(j, []).append(t)
    out = [nterm]
    for i, r in enumerate(seq):
        me = r.startswith('Me')
        r = r[2:] if me else r
        n = 'N(C)' if me else 'N'
        if cyclic and i == 0:
            n = n[0] + '%99' + n[1:]
        if i in special:
            out.append(special[i])
        elif r in ('G', 'Sar'):
            out.append(n + 'CC(=O)')
        elif r == 'P':
            t = next(closures)
            out.append(f'N%{t}CCCC%{t}C(=O)')
        elif r == 'Aib':
            out.append(n + 'C(C)(C)C(=O)')
        elif r == 'Ac3c':
            t = next(closures)
            out.append(n + f'C%{t}(CC%{t})C(=O)')
        elif r == 'C':
            out.append(n + 'C(CS' + ''.join(f'%{t}' for t in tags[i]) + ')C(=O)')
        else:
            out.append(n + 'C(' + SIDE[r] + ')C(=O)')
    return ''.join(out) + ('%99' if cyclic else cterm)


# semaglutide's Lys20 side chain: two OEG spacers, gamma-Glu, the C18 diacid
LINKER = 'CCCCNC(=O)COCCOCCNC(=O)COCCOCCNC(=O)CCC(NC(=O)CCCCCCCCCCCCCCCCC(=O)O)C(=O)O'
SEMA = ['H', 'Aib', *'EGTFTSDVSSYLEGQAA', 'K', *'EFIAWLVRGRG']
CASES = {
    'leu-enkephalin': peptide('YGGFL'),
    'oxytocin': peptide('CYIQNCPLG', ss=[(0, 5)], cterm='N'),
    'linaclotide': peptide('CCEYCCNPACTGCY', ss=[(0, 5), (1, 9), (4, 12)]),
    'gramicidin S': peptide(['V', 'O', 'L', 'F', 'P'] * 2, cyclic=True),
    'cyclosporin A': peptide(['MeBmt', 'Abu', 'MeSar', 'MeL', 'V', 'MeL', 'A', 'A', 'MeL', 'MeL', 'MeV'],
                             cyclic=True),
    'semaglutide': peptide(SEMA, special={19: 'NC(' + LINKER + ')C(=O)'}),
    'beta/gamma mix': 'NCC(C)C(=O)NC(Cc1ccccc1)CC(=O)NCCCC(=O)NC(C)C(=O)NC(CC(C)C)CC(=O)O',
}
PEPTIDES = tuple(CASES)
# not in PEPTIDES: below the threshold, or edge topology
CASES['tripeptide'] = peptide('GFL')
CASES['tetrapeptide'] = peptide('GFLA')
CASES['Ac3c'] = peptide(['A', 'Ac3c', 'F', 'Ac3c', 'L'])
# a 3-residue branch on Lys3's epsilon-N
CASES['branched'] = peptide('AGKFL', special={2: 'NC(CCCCNC(=O)C(C)NC(=O)CNC(=O)C(N)C)C(=O)'})
# a Lys-Asp side-chain lactam, i to i+4
CASES['lactam'] = peptide('AKGFLDA', special={1: 'NC(CCCCN%50)C(=O)', 5: 'NC(CC%50=O)C(=O)'})
# staples: a triazole (Lys(N3)-Pra click) and an i to i+4 all-hydrocarbon E-alkene
CASES['triazole staple'] = peptide('AGAFLAGA', special={2: 'NC(CCCn8cc%51nn8)C(=O)', 6: 'NC(CC%51)C(=O)'})
CASES['alkene staple'] = peptide('AGAFLAGA', special={2: 'NC(C)(CCC/C=C/%50)C(=O)', 6: 'NC(C)(CCC%50)C(=O)'})
# Glu3's side chain acylates the epsilon-N of Lys, the N-terminus of a Lys-Gly-Ala-Val branch
CASES['side-chain branch'] = peptide('AGEFL', special={
    2: 'NC(CCC(=O)NCCCCC(N)C(=O)NCC(=O)NC(C)C(=O)NC(C(C)C)C(=O)O)C(=O)'})
CASES['Pro-Pro'] = peptide('APPA')
CASES['collagen'] = peptide('GPPGPPGPPG')
CASES['C-terminal Pro-Pro'] = peptide('PAPP')
CASES['Pro3'] = peptide('GAPPPAG')
CASES['Pro4'] = peptide('AGPPPPGA')
# 4-aminobenzoic acid inside the chain: the backbone path runs through the ring
CASES['4-aminobenzoyl'] = peptide('AGAFLA', special={2: 'Nc1ccc(cc1)C(=O)'})


def molecule(name):
    return smiles(CASES[name])


def explicit(mol):
    """`mol` with every implicit hydrogen an atom."""
    counts = {n: mol.implicit_h_of(n) for n in mol}
    with mol.edit():
        for n, k in counts.items():
            for _ in range(k):
                mol.add_bond(n, mol.add_atom('H'), 1)
    return mol
