# -*- coding: utf-8 -*-
from io import BytesIO, StringIO
from chython import smiles
from chython.files import MRVRead, MRVWrite


def test_mrv_read_molecule():
    with MRVRead('test/implicit.mrv') as f:
        mol = next(f)
    assert mol is not None
    assert len(mol) > 0
    assert 'N' in {a.atomic_symbol for _, a in mol.atoms()}


def test_mrv_read_all():
    with MRVRead('test/implicit.mrv') as f:
        mols = f.read()
    assert len(mols) == 1


def test_mrv_write_read_roundtrip():
    mol = smiles('c1ccccc1')
    mol.clean2d()

    buf = StringIO()
    with MRVWrite(buf) as w:
        w.write(mol)

    xml_bytes = buf.getvalue().encode()
    with MRVRead(BytesIO(xml_bytes)) as r:
        mol2 = next(r)

    assert str(mol) == str(mol2)


def test_mrv_write_reaction():
    rxn = smiles('[CH3:1][OH:2]>>[CH3:1][NH2:3]')
    for m in rxn.molecules():
        m.clean2d()

    buf = StringIO()
    with MRVWrite(buf) as w:
        w.write(rxn)

    assert '<reaction>' in buf.getvalue()
    assert 'reactantList' in buf.getvalue()
    assert 'productList' in buf.getvalue()
