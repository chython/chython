# -*- coding: utf-8 -*-
"""Agreement and timing of `salt_estimator_proto` against the RDKit code it replaces.

The RDKit half is Ivan's code verbatim, so the comparison is against what runs today and not against a
paraphrase of it.  Probes are public drugs and amino acids.
"""
from time import perf_counter
from rdkit import Chem
from rdkit import RDLogger
from chython import smiles as chython_smiles
from salt_estimator_proto import salt_equivalents, ionizable_sites

RDLogger.DisableLog('rdApp.*')

# --- Ivan's code, verbatim ------------------------------------------------------------------------ #
FG_FOR_BASIC = Chem.MolFromSmarts("[#6]-[#6](=[#8])-[#8H]")
FG_FOR_FA_LIST = [
    "[CX4]-[#7H2]", "[CX4]-[#7H]-[CX4]", "[CX4]-[#7](-[CX4])-[CX4]",
    "[#6]1=,:[#6][#6]=,:[#7][#6]=,:[#6]1", "[#6]1=,:[#6][#7][#6]=,:[#7]1",
    "[#7]=[#6](-[#7])-[#7]", "[#7H3]", "C-[#7]=[#6](-[#6])-[#6]",
    "[#6]-C(=[#7])-[#7](-[#6X4])-[#6X4]",
]
FG_FOR_TFA_LIST = FG_FOR_FA_LIST[:6] + [
    "[#6]=,:1[#6]=,:[#7][#6]=,:[#7][#6]=,:1", "[#6]=,:1[#6]=,:[#7][#7][#6]=,:1",
    "[#6]=,:1[#6]=,:[#7][#6]=,:[#6][#7]=,:1",
] + FG_FOR_FA_LIST[6:]
_FA = [Chem.MolFromSmarts(s) for s in FG_FOR_FA_LIST]
_TFA = [Chem.MolFromSmarts(s) for s in FG_FOR_TFA_LIST]


def rdkit_estimate(s):
    mol = Chem.MolFromMolBlock(s) or Chem.MolFromSmiles(s)
    if not mol:
        raise ValueError(f'Invalid molecule string: {s}')
    return (sum(len(mol.GetSubstructMatches(fg, uniquify=True)) for fg in _TFA),
            sum(len(mol.GetSubstructMatches(fg, uniquify=True)) for fg in _FA),
            len(mol.GetSubstructMatches(FG_FOR_BASIC, uniquify=True)))


def chython_estimate(s):
    m = chython_smiles(s)
    m.thiele()                        # the patterns are aromatic; a Kekule record needs this and no more
    return tuple(salt_equivalents(m, t) for t in ('TFA', 'HCOOH', 'NH3'))


PROBES = {
    'aspirin':          'CC(=O)Oc1ccccc1C(=O)O',
    'caffeine':         'Cn1cnc2c1c(=O)n(C)c(=O)n2C',
    'nicotine':         'CN1CCC[C@H]1c1cccnc1',
    'lysine':           'NCCCC[C@H](N)C(=O)O',
    'histidine':        'N[C@@H](Cc1c[nH]cn1)C(=O)O',
    'glycine':          'NCC(=O)O',
    'gabapentin':       'NCC1(CC(=O)O)CCCCC1',
    'metformin':        'CN(C)C(=N)NC(N)=N',
    'diphenhydramine':  'CN(C)CCOC(c1ccccc1)c1ccccc1',
    'ciprofloxacin':    'O=C(O)c1cn(C2CC2)c2cc(N3CCNCC3)c(F)cc2c1=O',
    'losartan':         'CCCCc1nc(Cl)c(CO)n1Cc1ccc(-c2ccccc2-c2nn[nH]n2)cc1',
    'celecoxib':        'Cc1ccc(-c2cc(C(F)(F)F)nn2-c2ccc(S(N)(=O)=O)cc2)cc1',
    '4-aminopyridine':  'Nc1ccncc1',
    'pyrimidine':       'c1cncnc1',
    'pyrazole':         'c1cc[nH]n1',
    'aniline':          'Nc1ccccc1',
    'imatinib':         'Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1',
    'acetazolamide':    'CC(=O)Nc1nnc(S(N)(=O)=O)s1',
    'benzenesulfonic':  'OS(=O)(=O)c1ccccc1',
    'phenol':           'Oc1ccccc1',
    'piperazine':       'C1CNCCN1',
    'ethylenediamine':  'NCCN',
    'adenine':          'Nc1ncnc2[nH]cnc12',
}

print(f'{"probe":18} {"TFA":>10} {"FA":>10} {"NH3":>10}   {"chython sites"}')
print(f'{"":18} {"rd/chy":>10} {"rd/chy":>10} {"rd/chy":>10}')
diff = 0
for name, s in PROBES.items():
    r, c = rdkit_estimate(s), chython_estimate(s)
    m = chython_smiles(s); m.thiele()
    sites = ','.join(f'{rid}' for _, _, _, rid in ionizable_sites(m)) or '-'
    flag = '' if r == c else '  <-- DIFF'
    diff += r != c
    print(f'{name:18} {r[0]:4}/{c[0]:<5} {r[1]:4}/{c[1]:<5} {r[2]:4}/{c[2]:<5}   {sites}{flag}')
print(f'\n{len(PROBES) - diff}/{len(PROBES)} probes agree on all three numbers')

# --- timing -------------------------------------------------------------------------------------- #
strings = list(PROBES.values()) * 20
for label, fn in (('rdkit  ', rdkit_estimate), ('chython', chython_estimate)):
    fn(strings[0])
    t = perf_counter()
    for s in strings:
        fn(s)
    dt = perf_counter() - t
    print(f'{label}  {dt / len(strings) * 1e6:8.1f} us/molecule  (parse + all three counts)')
