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
"""CDXML both ways against RDKit's ``rdChemDraw``, which wraps Revvity's open-source ChemDraw library.

RDKit sits on both sides of every comparison, so only the CDXML layer is under test:

| Direction | Compared |
|---|---|
| out | RDKit reading our CDXML, against RDKit reading our molfile of the same molecule |
| in | chython reading RDKit's CDXML, against chython reading RDKit's molfile of the same molecule |

Binary CDX is not compared: ``rdChemDraw`` in RDKit 2026.03.4 returns no molecule for a CDX file,
Revvity-written ones included.  Absent RDKit, or an RDKit without ``rdChemDraw``, skips.
"""

from pytest import fixture, importorskip, mark

from chython.core import smiles

from .. import cdxml
from ...ctfile import mol
from ...test.oracles import require


#: Public compounds; each exercises one thing the CDXML layer has to carry.
COMPOUNDS = {
    'L-alanine': 'C[C@H](N)C(=O)O',
    '(R)-bromochlorofluoromethane': 'F[C@@H](Cl)Br',
    'caffeine (aromatic)': 'Cn1cnc2c1c(=O)n(C)c(=O)n2C',
    '(E)-but-2-ene': 'C/C=C/C',
    'sodium acetate': 'CC(=O)[O-].[Na+]',
    'deuterochloroform': '[2H]C(Cl)(Cl)Cl',
    'pyridinium': 'c1cc[nH+]cc1',
    'cisplatin (dative)': 'N->[Pt](<-N)(Cl)Cl',
    'ibuprofen': 'CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O',
    'threonine, AND group': 'C[C@@H](O)[C@H](N)C(=O)O |&1:1,3|',
}


@fixture(scope='module')
def chemdraw():
    require('rdkit')
    module = importorskip('rdkit.Chem.rdChemDraw')
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem
    RDLogger.DisableLog('rdApp.*')
    return Chem, AllChem, module


def _drawn(Chem, AllChem, text):
    """The compound with RDKit's 2D layout, as chython reads RDKit's V3000 of it."""
    rd = Chem.MolFromSmiles(text)
    AllChem.Compute2DCoords(rd)
    return rd, mol(Chem.MolToV3KMolBlock(rd))


def _rd_smiles(Chem, molecules):
    return '.'.join(sorted(Chem.MolToSmiles(m) for m in molecules))


@mark.parametrize('name', sorted(COMPOUNDS))
def test_rdkit_reads_our_cdxml_as_it_reads_our_molfile(chemdraw, name):
    Chem, AllChem, cd = chemdraw
    _, ours = _drawn(Chem, AllChem, COMPOUNDS[name])
    via_cdxml = _rd_smiles(Chem, cd.MolsFromChemDrawBlock(cdxml(ours)))
    via_molfile = _rd_smiles(Chem, [Chem.MolFromMolBlock(mol(ours, version=3000))])
    if 'dative' in name:
        # Our molfile states order 8 as type 8, which RDKit reads as the query "any"; CDXML says dative.
        via_molfile = via_molfile.replace('~', '')
        via_cdxml = via_cdxml.replace('->', '').replace('<-', '')
    assert via_cdxml == via_molfile


@mark.parametrize('name', sorted(COMPOUNDS))
def test_we_read_rdkit_cdxml_as_we_read_its_molfile(chemdraw, name):
    Chem, AllChem, cd = chemdraw
    rd, _ = _drawn(Chem, AllChem, COMPOUNDS[name])
    via_cdxml = cdxml(cd.MolToChemDrawBlock(rd, cd.CDXFormat.CDXML))
    via_molfile = mol(Chem.MolToV3KMolBlock(rd))
    for m in (via_cdxml, via_molfile):
        m.kekule()
        m.canonicalize_stereo_groups()
    assert str(via_cdxml) == str(via_molfile)


def test_rdkit_reads_our_reaction(chemdraw):
    Chem, AllChem, cd = chemdraw
    reaction = smiles('[CH3:1][C:2](=O)O.C[OH:3]>>[CH3:1][C:2](=O)[O:3]C')
    for m in reaction.molecules():
        rd = Chem.MolFromMolBlock(mol(m, version=3000))
        AllChem.Compute2DCoords(rd)
        layout = mol(Chem.MolToV3KMolBlock(rd))
        with m.edit():
            for sid, other in zip(m.atom_numbers, layout.atom_numbers):
                m.set_xy(sid, *layout.xy_of(other))
    (theirs,) = cd.ReactionsFromChemDrawBlock(cdxml(reaction))
    assert theirs.GetNumReactantTemplates() == 2 and theirs.GetNumProductTemplates() == 1
    product = theirs.GetProductTemplate(0)
    assert sorted(a.GetAtomMapNum() for a in product.GetAtoms() if a.GetAtomMapNum()) == [1, 2, 3]
