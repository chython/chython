# -*- coding: utf-8 -*-
from io import StringIO
from chython import smiles
from chython.files import SDFRead, SDFWrite, RDFRead, RDFWrite, mdl_mol


# V3000 propane with a DAT SGROUP whose FIELDDISP wraps across two lines via a
# trailing `-` continuation marker. mdl_mol feeds parse_mol_v3000 via
# str.splitlines() (no trailing newline), while SDFRead keeps newlines; the
# continuation detection must work for both. Regression for a crash when the
# continued line was not rejoined (ValueError in the SGROUP kv split).
_V3000_CONTINUATION = '''
  test              2D

  0  0  0  0  0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 3 2 1 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0 0 0 0
M  V30 2 C 1 0 0 0
M  V30 3 C 2 0 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 DAT 1 ATOMS=(1 1) FIELDDISP="    0.0000    0.0000    DA    ALL  1    -
M  V30    5" FIELDDATA=R FIELDNAME=STEREOLABEL
M  V30 END SGROUP
M  V30 END CTAB
M  END
'''


def test_mdl_mol_v3000_line_continuation():
    # mdl_mol splits on '\n'-stripped lines; the continued SGROUP must still rejoin
    mol = mdl_mol(_V3000_CONTINUATION)
    assert str(mol) == 'CCC'


def test_sdf_read_v3000_line_continuation():
    # SDFRead feeds lines with newlines kept; same record must parse identically
    with SDFRead(StringIO(_V3000_CONTINUATION + '$$$$\n')) as f:
        mol = next(f)
    assert str(mol) == 'CCC'


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
