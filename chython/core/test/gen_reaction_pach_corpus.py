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
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with this program; if not, see <https://www.gnu.org/licenses/>.
#
"""Write `reaction_pach_v2_corpus.bin.gz` from an INSTALLED chython 2.24.

Run from the repository root, AS A MODULE -- it reaches `oracle.py` through a relative import:

    python -m chython.core.test.gen_reaction_pach_corpus

THE chython 2 CODE IS A SUBPROCESS PAYLOAD AND NOT AN IMPORT, and that is not a style choice: nothing
under `chython/core/` may import the chython 2 facade -- see `test_no_chython_two_imports.py` -- so the
oracle is reached through `oracle.run`, which puts `-I` on the child and checks that the child did not
import the tree under test.  A generator that said `from chython import smiles` would be the very
blocker that test exists to keep removed.

THE OUTPUT IS COMMITTED AND MUST NOT BE REGENERATED FROM A LATER CHYTHON.  See
`reaction_pach_corpus.py` for what the records mean and why.

THE REACTIONS ARE TEXTBOOK AND PUBLIC, every one of them, and they were chosen to cover the layer this
corpus pins rather than to be interesting chemistry: all three sides populated and only two of them;
one, two and five molecules on a side; a mapping that relates the sides and a reaction with no mapping
at all; charge, isotope, aromatic rings and a tetrahedral centre, so the molecule layer is exercised
through the reaction wrapper.  NO RECORD HAS AN EMPTY PRODUCT SIDE -- chython 2's reader cannot read
one back, so there is no chython 2 answer for it to freeze.
"""
import base64
import gzip
import json
from struct import pack as struct_pack

from . import oracle


# name -> reaction SMILES.  Public, textbook, and each one named for what it pins.
REACTIONS = {
    # a mapped amide coupling with an agent: three populated sides, 2 / 1 / 2 molecules
    'amidation_mapped_with_agent':
        '[CH3:1][C:2](=[O:3])[OH:4].[NH2:5][CH3:6]'
        '>[CH3:10][CH2:11][OH:12]'
        '>[CH3:1][C:2](=[O:3])[NH:5][CH3:6].[OH2:4]',
    # Fischer esterification, mapped, no agents: the middle count is 0 and the header must say so
    'esterification_mapped_no_agent':
        '[CH3:1][C:2](=[O:3])[OH:4].[CH3:5][OH:6]>>[CH3:1][C:2](=[O:3])[O:4][CH3:5].[OH2:6]',
    # the same reaction with no mapping at all: chython 2 numbers the atoms 1..N across the whole
    # record, which is what makes "the number field is the mapping" ambiguous, and the ambiguity is
    # pinned here rather than argued about
    'esterification_unmapped':
        'CC(=O)O.CO>>CC(=O)OC.O',
    # Suzuki coupling: two aromatic rings through the wrapper, mapped across the arrow
    'suzuki_mapped':
        '[cH:1]1[cH:2][cH:3][c:4]([Br:20])[cH:5][cH:6]1.[cH:7]1[cH:8][cH:9][c:10]([B:21]([OH:22])'
        '[OH:23])[cH:11][cH:12]1>[Pd]>[cH:1]1[cH:2][cH:3][c:4]([c:10]2[cH:9][cH:8][cH:7][cH:12]'
        '[cH:11]2)[cH:5][cH:6]1',
    # nitration of toluene: an aromatic ring, a charge-separated nitro group, five molecules on the
    # left so the reactant count is not 1 or 2
    'toluene_nitration':
        '[cH:1]1[cH:2][cH:3][c:4]([CH3:5])[cH:6][cH:7]1.[N+:8](=[O:9])([O-:10])[OH:11].[OH2:12]'
        '.[OH2:13].[OH2:14]'
        '>[S:30](=[O:31])(=[O:32])([OH:33])[OH:34]'
        '>[c:1]1([N+:8](=[O:9])[O-:10])[cH:2][cH:3][c:4]([CH3:5])[cH:6][cH:7]1.[OH2:11]',
    # a tetrahedral centre either side of the arrow: (S)-lactic acid esterified
    'lactic_acid_esterification_stereo':
        '[CH3:1][C@H:2]([OH:3])[C:4](=[O:5])[OH:6].[CH3:7][OH:8]'
        '>>[CH3:1][C@H:2]([OH:3])[C:4](=[O:5])[O:6][CH3:7].[OH2:8]',
    # a cis/trans centre: maleic to fumaric acid, the isomerisation textbooks use
    'maleic_to_fumaric_cis_trans':
        '[OH:1][C:2](=[O:3])/[CH:4]=[CH:5]\\[C:6](=[O:7])[OH:8]'
        '>>[OH:1][C:2](=[O:3])/[CH:4]=[CH:5]/[C:6](=[O:7])[OH:8]',
    # charge and isotope: sodium acetate from 1-13C acetic acid
    'isotope_and_charge':
        '[13CH3:1][C:2](=[O:3])[OH:4].[Na+:5].[OH-:6]'
        '>>[13CH3:1][C:2](=[O:3])[O-:4].[Na+:5].[OH2:6]',
    # a lone metal cation, which chython 2's own `pack(check=True)` refused for having no bonds and
    # which the format has always been able to hold: one atom, zero bonds, on the agent side
    'lone_cation_agent':
        '[CH3:1][CH2:2][Br:3].[OH-:4]>[K+:40]>[CH3:1][CH2:2][OH:4].[Br-:3]',
    # the radical bit, on a one-atom molecule with no bonds: homolysis of ethane, the textbook
    # illustration of it.  The CXSMILES tail indexes atoms across the WHOLE reaction string.
    'ethane_homolysis_radical':
        '[CH3:1][CH3:2]>>[CH3:1].[CH3:2] |^1:2,3|',
}


# Runs INSIDE the chython 2 interpreter. `check=False` on `pack` so that the lone-cation record can be
# written at all: `pack(check=True)` there refuses a molecule with no bonds, which is a restriction the
# format itself does not carry.
PAYLOAD = r'''
import base64, json, sys
from chython import smiles, ReactionContainer

out = []
for name, spec in json.loads(sys.stdin.read()):
    rxn = smiles(spec)
    data = rxn.pack(compressed=False, check=False)
    back = ReactionContainer.unpack(data, compressed=False)
    molecules = []
    for mol in back.molecules():
        atoms = sorted([n, a.atomic_number, a.isotope, a.charge, int(a.is_radical),
                        a.implicit_hydrogens, len(mol._bonds[n])] for n, a in mol.atoms())
        bonds = sorted([min(n, m), max(n, m), b.order] for n, m, b in mol.bonds())
        molecules.append({'atoms': atoms, 'bonds': bonds})
    counts = [len(back.reactants), len(back.reagents), len(back.products)]
    lens = ReactionContainer.pack_len(data, compressed=False)
    out.append({'name': name, 'smiles': spec, 'data': base64.b64encode(data).decode(),
                'counts': counts, 'molecules': molecules,
                'atom_counts': [list(lens[0]), list(lens[1]), list(lens[2])]})
sys.stdout.write('\x1e' + json.dumps(out))
'''


def main():
    if oracle.interpreter() is None:
        raise SystemExit('the chython 2 oracle is not provisioned; see chython.core.test.oracle')
    oracle.verify()
    result = oracle.run('-c', PAYLOAD, input=json.dumps(list(REACTIONS.items())))
    if result.returncode:
        raise SystemExit('the chython 2 oracle failed:\n%s' % result.stderr)
    payload = json.loads(result.stdout.rsplit('\x1e', 1)[1])

    blob = bytearray(struct_pack('<I', len(payload)))
    for record in payload:
        data = base64.b64decode(record.pop('data'))
        name = record.pop('name').encode()
        assert data[0] == 1, name
        assert record['counts'][2], 'no record may have an empty product side; see the module docstring'
        answers = json.dumps(record, sort_keys=True).encode()
        blob += struct_pack('<III', len(name), len(data), len(answers))
        blob += name + data + answers
        print('%-42s %4d bytes  %s' % (name.decode(), len(data), record['counts']))

    path = oracle.ROOT / 'chython' / 'core' / 'test' / 'reaction_pach_v2_corpus.bin.gz'
    with gzip.GzipFile(path, 'wb', mtime=0) as f:      # mtime=0: the artefact is reproducible
        f.write(bytes(blob))
    print('\nwrote %s: %d records, %d bytes gzipped' % (path.name, len(payload), path.stat().st_size))


if __name__ == '__main__':
    main()
