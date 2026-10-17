"""Build test/layout_corpus.tsv. Throwaway: the TSV is the artefact. Needs the network."""
import json
import random
import subprocess
import sys
import time
import urllib.parse
from pathlib import Path

from chython.core import ReactionContainer, smiles
from chython.core.test.peptides import CASES

OUT = Path(__file__).resolve().parents[1] / 'test' / 'layout_corpus.tsv'
PUG = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug'
UA = {'User-Agent': 'chython-layout-corpus/1.0'}

HARD = ['benzene', 'naphthalene', 'anthracene', 'phenanthrene', 'pyrene', 'triphenylene', 'coronene',
        'azulene', 'indole', 'purine', 'pteridine', 'biphenylene', 'fluorene', 'acenaphthylene',
        'cyclopropane', 'cyclobutane', 'cyclooctane', 'cyclooctatetraene', 'adamantane', 'cubane',
        'norbornane', 'bicyclo[2.2.2]octane', 'twistane', 'spiro[4.5]decane', 'spiropentane', 'decalin',
        'cholesterol', 'testosterone', 'estradiol', 'cortisol', 'morphine', 'codeine', 'quinine',
        'strychnine', 'reserpine', 'taxol', 'artemisinin', 'camphor', 'menthol', 'glucose', 'sucrose',
        'alpha-cyclodextrin', 'adenosine', 'ATP', 'NAD+', 'thymidine', 'porphine', 'chlorophyll a', 'heme b',
        'vitamin B12', '18-crown-6', '15-crown-5', '[2.2.2]cryptand', '[2.2]paracyclophane',
        'calix[4]arene', 'ferrocene', 'cisplatin', 'tris(bipyridine)ruthenium(II) chloride',
        'ethylenediaminetetraacetic acid', 'hexamethylenetetramine', 'buckminsterfullerene', 'caffeine',
        'aspirin', 'ibuprofen', 'penicillin G', 'erythromycin', 'vancomycin', 'rapamycin', 'cyclosporin A',
        'tetracycline', 'doxorubicin', 'warfarin', 'diazepam', 'sildenafil', 'imatinib', 'atorvastatin',
        'retinol', 'beta-carotene', 'squalene', 'triolein']

SCAFFOLDS = {                                           # `({R})` is the attachment; `core` drops it
    'benzene': 'c1ccc({R})cc1', 'indole': 'c1ccc2c(c1)cc({R})[nH]2', 'quinoline': 'c1ccc2ncc({R})cc2c1',
    'piperazine': 'C1CN(c2ccccc2)CCN1({R})', 'biphenyl': 'c1ccc(-c2ccc({R})cc2)cc1',
    'pyrimidine': 'Nc1nc({R})ccn1', 'benzimidazole': 'c1ccc2[nH]c({R})nc2c1',
    'amide': 'O=C(N({R}))c1ccc2ccccc2c1', 'purine': 'Nc1ncnc2c1ncn2({R})',
    'androstane': 'CC12CCC3C(CCC4CC({R})CCC34C)C1CCC2'}
RGROUPS = ['C', 'CC', 'CCC', 'C(C)C', 'C(C)(C)C', 'C9CC9', 'C9CCCC9', 'C9CCCCC9', 'c9ccccc9', 'c9ccncc9',
           'OC', 'N(C)C', 'C(=O)O', 'C(=O)N', 'C#N', 'C(F)(F)F', 'S(=O)(=O)C', 'CCO', 'C9CCOCC9',
           'Cc9ccccc9']


def get(url, tries=6):
    """GET through curl: from this host PubChem answers urllib 503 ServerBusy and curl 200."""
    for k in range(tries):
        r = subprocess.run(['curl', '-s', '-f', '-m', '120', '-A', UA['User-Agent'], url], capture_output=True)
        if r.returncode == 0:
            return r.stdout
        if k == tries - 1 or r.returncode == 22 and b'404' in r.stderr:
            raise OSError(f'curl {r.returncode} on {url[:80]}')
        print('retry', k, url[:80], r.returncode, file=sys.stderr)
        time.sleep(2 ** k)


def props(cids):
    url = f'{PUG}/compound/cid/{",".join(map(str, cids))}/property/SMILES,HeavyAtomCount/JSON'
    try:
        return json.loads(get(url))['PropertyTable']['Properties']
    except Exception:
        return []


def ring_systems(mol):
    systems = []
    for ring in mol.sssr:
        ring = set(ring)
        for s in [s for s in systems if s & ring]:
            systems.remove(s)
            ring |= s
        systems.append(ring)
    return len(systems)


def parses(smi):
    try:
        return smiles(smi)
    except Exception:
        return None


def hard():
    for name in HARD:
        try:
            p, *_ = json.loads(get(f'{PUG}/compound/name/{urllib.parse.quote(name)}/property/SMILES/JSON'))[
                'PropertyTable']['Properties']
        except Exception as e:
            print('hard: skipped', name, e, file=sys.stderr)
            continue
        yield f'hard:{name}', 'hard', p['SMILES'], ''
        time.sleep(.25)
    yield 'hard:40-fragment salt', 'hard', '.'.join(['[Na+].[Cl-]'] * 20), ''


def random_pubchem(rng):
    buckets = {k: [] for k in range(5)}
    large = []
    while any(len(b) < 200 for b in buckets.values()) or len(large) < 200:
        for p in props(rng.sample(range(1, 50_000_000), 100)):
            n, smi = p.get('HeavyAtomCount', 0), p.get('SMILES')
            mol = smi and '.' not in smi and parses(smi)
            if not mol:
                continue
            if 10 <= n < 50:
                b = buckets[min(ring_systems(mol), 4)]
                if len(b) < 200:
                    b.append((f'pubchem:{p["CID"]}', 'drug', smi, ''))
            elif 50 <= n <= 150 and len(large) < 200:
                large.append((f'pubchem:{p["CID"]}', 'large', smi, ''))
        time.sleep(.25)
    for b in buckets.values():
        yield from b
    yield from large


METALS = {21, 22, 23, 24, 25, 26, 27, 28, 29, 40, 41, 42, 43, 44, 45, 46, 72, 73, 74, 75, 76, 77, 78}


def metal(rng):
    cids = set()                                            # one query per element: each answers ~400 hits at most
    for z in sorted(METALS):
        q = urllib.parse.quote(f'[#{z};D5,D6]')
        cids.update(json.loads(get(f'{PUG}/compound/fastsubstructure/smarts/{q}/cids/JSON?MaxRecords=5000'))[
            'IdentifierList'].get('CID', []))
        time.sleep(.25)
    cids = sorted(cids)
    rng.shuffle(cids)
    out = []
    for k in range(0, len(cids), 100):
        for p in props(cids[k:k + 100]):
            mol = parses(p.get('SMILES', ''))                # PubChem's SMILES may disconnect the centre it matched
            if mol and p.get('HeavyAtomCount', 0) <= 80 and any(
                    mol.element_of(x) in METALS and len(mol.neighbors_of(x)) >= 5 for x in mol):
                out.append((f'pubchem:{p["CID"]}', 'metal', p['SMILES'], ''))
        if len(out) >= 100:
            return out[:100]
        time.sleep(.25)
    return out


RETROSIM = 'https://raw.githubusercontent.com/connorcoley/retrosim/master/retrosim/data/data_processed.csv'


def uspto(rng):
    """USPTO-50k (Schneider et al. 2016, drawn from Lowe's CC0 USPTO extraction), as shipped by retrosim (MIT)."""
    lines = get(RETROSIM).decode().splitlines()[1:]
    rng.shuffle(lines)
    out = []
    for line in lines:
        k, _, _, _, smi = line.split(',')[:5]
        r = parses(smi)
        if isinstance(r, ReactionContainer) and len(r.products) == 1 and 8 <= len(r.products[0]) <= 50:
            out.append((f'rxn:{k}', 'rxn', smi, ''))
            if len(out) == 300:
                return out
    return out


def series():
    for name, tpl in SCAFFOLDS.items():
        core = tpl.replace('({R})', '')
        for k, r in enumerate(RGROUPS):
            yield f'ser:{name}:{k}', 'ser', tpl.replace('{R}', r), core


def main():
    rng = random.Random(0)
    rows = [*hard(), *random_pubchem(rng), *metal(rng),
            *((f'peptide:{n}', 'peptide', s, '') for n, s in CASES.items()),
            *uspto(rng), *series()]
    bad = [r[0] for r in rows if parses(r[2]) is None or (r[3] and parses(r[3]) is None)]
    assert not bad, bad
    OUT.write_text('id\tstratum\tsmiles\tcore\n' + ''.join('\t'.join(r) + '\n' for r in rows))
    print(OUT, len(rows))


if __name__ == '__main__':
    main()
