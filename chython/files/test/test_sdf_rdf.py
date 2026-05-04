# -*- coding: utf-8 -*-
from io import StringIO
from chython import smiles
from chython.files import SDFRead, SDFWrite, RDFRead, RDFWrite


def test_sdf_read():
    with SDFRead('test/implicit.sdf') as f:
        mols = f.read()
    assert len(mols) > 0
    for mol in mols:
        assert len(mol) > 0


def test_sdf_read_stereo():
    with SDFRead('test/stereo.sdf') as f:
        mols = f.read(amount=10)
    assert len(mols) == 10


def test_sdf_write_read_roundtrip():
    mol = smiles('CCO')
    mol.clean2d()
    mol.meta['name'] = 'ethanol'

    buf = StringIO()
    with SDFWrite(buf) as w:
        w.write(mol)

    buf.seek(0)
    with SDFRead(buf) as r:
        mol2 = next(r)

    assert str(mol) == str(mol2)
    assert mol2.meta.get('name') == 'ethanol'


def test_rdf_read():
    from chython import ReactionContainer
    with RDFRead('test/MR.rdf') as f:
        records = f.read()
    assert len(records) > 0
    rxns = [r for r in records if isinstance(r, ReactionContainer)]
    assert len(rxns) > 0
    for rxn in rxns:
        assert len(rxn.reactants) > 0
        assert len(rxn.products) > 0


def test_rdf_write_read_roundtrip():
    rxn = smiles('[CH3:1][OH:2]>>[CH3:1][NH2:3]')
    for m in rxn.molecules():
        m.clean2d()

    buf = StringIO()
    with RDFWrite(buf) as w:
        w.write(rxn)

    buf.seek(0)
    with RDFRead(buf) as r:
        rxn2 = next(r)

    assert len(rxn2.reactants) == 1
    assert len(rxn2.products) == 1
