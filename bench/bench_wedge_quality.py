# Scratch metric harness for the wedge write path.  Prints the table the report quotes.
from collections import Counter

from chython.core import WEDGE_NONE
from chython.formats.ctfile._sdf import split_records
from chython.formats.ctfile._v2000 import parse_v2000
from chython.formats.ctfile._v3000 import V3000_STAMP, parse_v3000
from chython.formats.ctfile._sdf import sniff_version
from chython.formats.ctfile._wedge import SU_TETRA, tetrahedral_parity, wedges_for_write


def build(record):
    parse = parse_v3000 if sniff_version(record, []) == V3000_STAMP else parse_v2000
    return parse(record, []).build()


def load(path):
    with open(path, encoding='utf8', errors='replace') as f:
        return list(split_records(f))


def strip_wedges(mol):
    existing = list(mol.wedges())
    if existing:
        with mol.edit():
            for narrow, wide, _ in existing:
                mol.set_wedge(narrow, wide, WEDGE_NONE)
    return mol


def measure(records, plane=None):
    m = Counter()
    detail = []
    for record in records:
        try:
            mol, _, _ = build(record)
        except Exception as e:
            m['build_failed'] += 1
            continue
        title = (record[0] or '').strip()
        configured = [u for u in mol.stereo_units()
                      if u['kind'] == SU_TETRA and mol.parity_of(u['anchor'])]
        if not configured:
            continue
        m['molecules'] += 1
        m['centres'] += len(configured)
        target = {u['anchor']: mol.parity_of(u['anchor']) for u in configured}

        strip_wedges(mol)
        wedges, log = wedges_for_write(mol)

        ring = sum(1 for a, b, _ in wedges if mol.bond_in_ring(a, b))
        m['ring_wedges'] += ring
        # adjacent pair: two wedges sharing any atom
        adj = 0
        for i in range(len(wedges)):
            for j in range(i + 1, len(wedges)):
                if set(wedges[i][:2]) & set(wedges[j][:2]):
                    adj += 1
        m['adjacent_pairs'] += adj
        m['wedges'] += len(wedges)

        encoded = {a for a, b, w in wedges}
        unenc = sum(1 for anchor in target if anchor not in encoded)
        m['unencoded'] += unenc

        # round trip: apply the wedges and read the parities back
        with mol.edit():
            for a, b, w in wedges:
                mol.set_wedge(a, b, w)
        bad = 0
        for unit in mol.stereo_units():
            if unit['kind'] != SU_TETRA or unit['anchor'] not in target:
                continue
            if tetrahedral_parity(mol, unit) != target[unit['anchor']]:
                bad += 1
        if bad:
            m['roundtrip_failed'] += 1
            detail.append((title, 'roundtrip', bad))
        if unenc:
            detail.append((title, 'unencoded', unenc))
    return m, detail


if __name__ == '__main__':
    import sys
    recs = load(sys.argv[1] if len(sys.argv) > 1 else 'test/wedge_stereo.sdf')
    m, detail = measure(recs)
    print(f'{"metric":24} value')
    for k in ('molecules', 'centres', 'wedges', 'ring_wedges', 'adjacent_pairs',
              'unencoded', 'roundtrip_failed', 'build_failed'):
        print(f'{k:24} {m[k]}')
    if detail:
        print('\n-- detail --')
        for t in detail:
            print(' ', t)
