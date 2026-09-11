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
"""Every import bridge puts its records on the container it returned, with nothing passed in.

`mol.log` and `rxn.log` are the destination, so `interop.rdkit(rd_mol)` with no keyword must leave a
readable trace of what the conversion dropped.  Each test below calls the bridge bare and asks the
result; the `log=` list is checked only where it is checked as a COPY of the same records.

The export direction is not tested here and has nothing to deliver: it returns a foreign object, so
the caller's list is the only place its records can go.
"""
from pathlib import Path

from pytest import mark

from chython.core import LogRecord, read_smiles
from chython.interop._records import STAGE

from .conftest import requires_indigo, requires_jpype, requires_openbabel, requires_rdkit


_INTEROP = Path(__file__).resolve().parent.parent


def _check(container):
    """Records landed, every one is a `LogRecord`, and every one names this stage.  Returns them."""
    records = list(container.log)
    assert records, 'the import left no record on the container it returned'
    assert all(isinstance(r, LogRecord) for r in records)
    assert {r.stage for r in records} == {STAGE}
    return records


def _rules(records):
    return {r.rule for r in records}


# -- one per bridge, called with nothing but the foreign object ------------------------------------


@requires_rdkit
def test_rdkit_import_records_on_the_molecule():
    """An `RWMol` with no valence cache and an unstorable charge: two records, both on the molecule."""
    from rdkit.Chem import Atom, RWMol

    from chython.interop import rdkit

    rd = RWMol()
    a = Atom(6)
    a.SetFormalCharge(9)
    rd.AddAtom(a)

    mol = rdkit(rd)
    records = _check(mol)
    assert _rules(records) == {'rdkit:note'}
    assert any('valence cache' in r for r in records)
    assert any('clamped' in r for r in records)


@requires_rdkit
def test_rdkit_reaction_import_records_per_component():
    """A component keeps its own records and the reaction gets them stamped with which component.

    A reaction log pooling three sides would hand back atom numbers that name a different atom
    depending on which molecule they are read against; `subject` is what makes the copy readable.
    """
    from rdkit.Chem import Atom, RWMol
    from rdkit.Chem.rdChemReactions import ChemicalReaction

    from chython.interop import rdkit

    def broken():
        rd = RWMol()
        a = Atom(6)
        a.SetFormalCharge(9)
        rd.AddAtom(a)
        return rd

    rr = ChemicalReaction()
    rr.AddReactantTemplate(broken())
    rr.AddProductTemplate(broken())

    rxn = rdkit(rr)
    assert [r.subject for r in rxn.log] == ['reactants[0]'] * 2 + ['products[0]'] * 2
    assert {r.stage for r in rxn.log} == {STAGE}
    for where, molecule in (('reactants[0]', rxn.reactants[0]), ('products[0]', rxn.products[0])):
        mine = _check(molecule)
        # the same events, read from the two ends
        assert [r.message for r in rxn.log.by_subject(where)] == [r.message for r in mine]


@requires_indigo
def test_indigo_import_records_on_the_molecule():
    from chython.interop import indigo

    mol = indigo(indigo(read_smiles('CCO')))
    assert 'indigo:coordinates-not-imported' in _rules(_check(mol))


@requires_openbabel
def test_openbabel_import_records_on_the_molecule():
    from chython.interop import openbabel

    mol = openbabel(openbabel(read_smiles('CCO')))
    assert 'openbabel:coordinates-not-imported' in _rules(_check(mol))


@requires_jpype
def test_cdk_import_records_on_the_molecule(cdk):
    from chython.interop import cdk as interop_cdk

    mol = interop_cdk(interop_cdk(read_smiles('CCO')))
    assert 'cdk:coordinates-not-imported' in _rules(_check(mol))


def test_iupac_import_records_on_the_molecule(opsin):
    """The SMILES OPSIN produced is a record: it is what chython parsed, and the result does not
    otherwise say what it was."""
    from chython.interop import iupac

    mol = iupac('ethanol')
    records = _check(mol)
    assert 'iupac:parsed-by-opsin' in _rules(records)
    assert any('OPSIN' in r for r in records)


# -- the properties that hold across bridges ------------------------------------------------------


@requires_indigo
def test_the_caller_list_gets_a_copy_of_what_the_molecule_got():
    """`log=` is a copy for the caller, not an alternative destination: both hold the same records."""
    from chython.interop._indigo import from_indigo, to_indigo

    log = []
    mol = from_indigo(to_indigo(read_smiles('CCO')), log=log)
    assert log and log == list(mol.log)


@requires_indigo
def test_records_are_not_pooled_across_two_imports():
    """Molecule two's records go on molecule two.  A shared destination would be unreadable: the
    atom numbers in a record are stable ids in ONE container."""
    from chython.interop._indigo import from_indigo, to_indigo

    first = from_indigo(to_indigo(read_smiles('CCO')))
    before = len(first.log)
    second = from_indigo(to_indigo(read_smiles('c1ccccc1')))
    assert len(first.log) == before
    assert second.log and second is not first


# -- the ratchet, which does not need a toolkit ----------------------------------------------------


@mark.parametrize('filename,entry', [('_rdkit.py', 'def from_rdkit'),
                                     ('_indigo.py', 'def from_indigo'),
                                     ('_openbabel.py', 'def from_openbabel'),
                                     ('_cdk.py', 'def _from_cdk'),
                                     ('_iupac.py', 'def from_iupac')])
def test_every_importer_delivers_to_the_container(filename, entry):
    """Source-level, so it holds for a bridge whose toolkit is not installed on this machine.

    `_cdpkit.py` is absent by decision: `from_cdpkit` raises `DirectionNotImplemented` and returns no
    container.  Whoever builds that direction adds both the `deliver` call and a row here.
    """
    source = (_INTEROP / filename).read_text(encoding='utf8')
    start = source.index(entry)
    body = source[start:]
    end = body.find('\ndef ', len(entry))
    if end != -1:
        body = body[:end]
    assert 'deliver(' in body, f'{filename}: {entry} records nothing on the container it returns'
