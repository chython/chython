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
"""An atropisomer axis in the CXSMILES tail: `wU:a.b` / `wD:a.b`.

`a` is the wedge's narrow end and `b` the bond, both 0-based; bonds count in string order, a ring
closure at its closing digit.  Measured against RDKit 2026.03.4 on PubChem 3-D records of BINOL
(CID 11762) and BINAP (CID 634876):

| tail                   | RDKit reads the axis from                                                      |
| ---------------------- | ------------------------------------------------------------------------------ |
| with a `(...)` field   | the coordinates; the wedge only marks the bond as a candidate axis             |
| without coordinates    | the wedge alone: at each pivot the ring bond to its lower-index neighbour      |
|                        | wedged, or the other one hashed, is WEDGE; the converse is DASH.  WEDGE at the |
|                        | pivot written first, or DASH at the other, is CCW; the converse is CW          |

Against chython's parity, with `refs[0]` and `refs[2]` each pivot's lower-index ring neighbour: parity 1
is CCW and parity 2 is CW.  Each pivot whose `refs` pair is the other way round flips it once.

chython writes no coordinates, so the wedge is always read by the second row, and a tail carrying a
coordinates field has its wedges logged and not applied.  Compounds: BINOL (CID 11762), BINAP
(CID 634876), 6,6'-dinitrodiphenic acid (CID 95234); bromochlorofluoromethane and biphenyl for the
wedges that name no axis.
"""

from pytest import importorskip, mark

from .. import read_smiles, write_smiles
from ..wedge import atropisomer_parity


SU_ATROPISOMER = 3
COMPOUNDS = {
    'BINOL': 'Oc1ccc2ccccc2c1-c1c(O)ccc2ccccc12',
    'BINAP': 'c1ccc(cc1)P(c1ccccc1)c1ccc2ccccc2c1-c1c(ccc2ccccc12)P(c1ccccc1)c1ccccc1',
    'dinitrodiphenic acid': 'OC(=O)c1cccc(c1-c1c(cccc1[N+](=O)[O-])C(=O)O)[N+](=O)[O-]',
}


def _axis(mol):
    (unit,) = [u for u in mol.stereo_units() if u['kind'] == SU_ATROPISOMER and u['stereogenic']]
    return unit


def _configured(name, parity):
    mol = read_smiles(COMPOUNDS[name])
    mol.set_parity(_axis(mol)['anchor'], parity)
    return mol


def _rules(log):
    return [r.rule for r in log]


@mark.parametrize('name', sorted(COMPOUNDS))
def test_an_axis_survives_a_chython_cx_round_trip(name):
    texts = []
    for parity in (1, 2):
        mol = _configured(name, parity)
        log = []
        text = write_smiles(mol, log=log)
        assert ' |w' in text, text
        assert log == [], log
        back = read_smiles(text)
        assert list(back.log) == [], list(back.log)
        assert back.parity_of(_axis(back)['anchor']), text
        assert write_smiles(back) == text, 'the written string is a fixed point'
        assert back == mol
        texts.append(text)
    assert texts[0] != texts[1]
    assert read_smiles(texts[0]) != read_smiles(texts[1]), 'the two atropisomers are two compounds'


def test_plain_smiles_logs_the_axis_it_cannot_write():
    log = []
    text = write_smiles(_configured('BINOL', 1), '!x', log=log)
    assert '|' not in text
    assert _rules(log) == ['smiles:stereo-not-written'], log


def _parity_of_conformer(text, rd):
    """chython's parity of the axis in RDKit's conformer `rd`, whose atoms are `text`'s in order."""
    mol = read_smiles(text.split(' ')[0])
    unit = _axis(mol)
    conformer = rd.GetConformer()
    z = {}
    with mol.edit():
        for i, n in enumerate(mol.atom_numbers):
            p = conformer.GetAtomPosition(i)
            mol.set_xy(n, p.x, p.y)
            z[n] = p.z
    return atropisomer_parity(mol, unit, z=z)


@mark.parametrize('name', sorted(COMPOUNDS))
@mark.parametrize('parity', (1, 2))
def test_rdkit_reads_the_axis_chython_writes(name, parity):
    """RDKit reads chython's string, embeds it in 3-D, and chython reads the conformer back to the same
    parity; RDKit's own coordinate-free CXSMILES of that molecule reads back to it as well."""
    Chem = importorskip('rdkit.Chem')
    AllChem = importorskip('rdkit.Chem.AllChem')
    legacy = Chem.GetUseLegacyStereoPerception()
    Chem.SetUseLegacyStereoPerception(False)             # the perception that reads atropisomers
    try:
        _rdkit_round_trip(Chem, AllChem, name, parity)
    finally:
        Chem.SetUseLegacyStereoPerception(legacy)


def _rdkit_round_trip(Chem, AllChem, name, parity):
    mol = _configured(name, parity)
    text = write_smiles(mol)
    rd = Chem.MolFromSmiles(text)
    assert rd is not None, text
    assert [b.GetIdx() for b in rd.GetBonds() if 'ATROP' in str(b.GetStereo())], text
    rh = Chem.AddHs(rd)
    assert AllChem.EmbedMolecule(rh, randomSeed=42) == 0
    heavy = Chem.RemoveHs(rh)
    first = read_smiles(text)
    assert _parity_of_conformer(text, heavy) == first.parity_of(_axis(first)['anchor'])

    params = Chem.SmilesWriteParams()
    theirs = Chem.MolToCXSmiles(rd, params, Chem.CXSmilesFields.CX_ALL_BUT_COORDS)
    back = read_smiles(theirs)
    assert list(back.log) == [], (theirs, list(back.log))
    assert back == mol, theirs


def test_an_unmodelled_wedge_is_logged_not_applied():
    """A wedge on a tetrahedral centre, on a biaryl with no axis, and on an axis whose tail also
    carries coordinates: each stays in the log and stores nothing."""
    mol = read_smiles('F[C@H](Cl)Br |wU:1.0|')
    assert mol.parity_of(2), 'the atom sign stands'
    assert 'smiles:cx-wedge-not-applied' in _rules(mol.log), list(mol.log)

    mol = read_smiles('c1ccccc1-c1ccccc1 |wU:5.6|')
    assert all(not mol.parity_of(n) for n in mol.atom_numbers)
    assert 'smiles:cx-wedge-not-applied' in _rules(mol.log), list(mol.log)

    text = write_smiles(_configured('BINOL', 1))
    body, tail = text.split(' ')
    coords = '(' + ';'.join(['0,0,'] * 22) + ')'
    mol = read_smiles('%s |%s,%s' % (body, coords, tail[1:]))
    assert not mol.parity_of(_axis(mol)['anchor'])
    assert 'smiles:cx-wedge-not-applied' in _rules(mol.log), list(mol.log)
