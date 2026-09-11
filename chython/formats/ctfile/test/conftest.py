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
"""Fixtures for the CTfile tests.

The corpus is the repository's own ``test/*.sdf`` -- files other tools wrote, carrying the
malformations the recovery code exists for.  Nothing here builds a molecule from SMILES to test a
file reader: a hand-built structure cannot reproduce a vendor's idea of a valence column.
"""

from pathlib import Path

from pytest import fixture

from chython.core.test import oracle

from .._sdf import split_records


#: Files with more than a handful of records, used for the sweeps.  Named rather than globbed so a
#: new file appearing in `test/` cannot silently change what a passing suite means.
CORPUS = ('arenes.sdf', 'cycle.sdf', 'isomorphism.sdf', 'mcs.sdf', 'peptide.sdf',
          'standardize.sdf', 'stereo.sdf', 'implicit.sdf', 'depict.sdf', 'hbonds.sdf',
          'morgan_ruiner.sdf')


#: Is any bond drawn at `sid` an aromatic one?  An alias for the production predicate, deliberately
#: not a second copy: two suites use the answer as the licence for committing a hydrogen count, and a
#: local re-implementation could widen that escape hatch in one file while the other stayed narrow.
from .._hydrogens import _holds_aromatic_bond as holds_an_aromatic_bond


def _root():
    """The repository root, found by walking up from this file to the directory holding ``test/``."""
    for parent in Path(__file__).resolve().parents:
        if (parent / 'test').is_dir() and (parent / 'chython').is_dir():
            return parent
    raise RuntimeError('cannot locate the repository root from ' + __file__)


@fixture(scope='session')
def root():
    return _root()


@fixture(scope='session')
def corpus(root):
    """``{filename: [record lines, ...]}`` for every file in :data:`CORPUS` that exists."""
    out = {}
    for name in CORPUS:
        path = root / 'test' / name
        if not path.exists():
            continue
        with path.open(encoding='utf8', errors='replace') as f:
            out[name] = list(split_records(f))
    return out


@fixture(scope='session')
def oracle_session():
    """One chython 2 interpreter for the whole CTfile suite.  Skips when it is not provisioned.

    A subprocess and not an import: inside chython 2's own interpreter this tree's package does not
    exist, so an oracle cannot silently become a mirror of the reader it is checking.  See
    ``chython/core/test/oracle.py``.
    """
    live = oracle.session()
    yield live
    live.close()


@fixture(scope='session')
def v2_reader(oracle_session):
    """``f(record_text) -> chython 2 record`` for one hand-written MDL record.

    What comes back is an ``oracle.Record``, answering under chython 2's own attribute names
    (``atom(n).implicit_hydrogens``, ``bonds()``, iteration over atom numbers).  Not a molecule.
    """
    def read(text):
        record, = oracle_session.read_mdl([text])
        assert record is not None, 'chython 2 could not read this fixture record'
        return record
    return read


@fixture(scope='session')
def v2_molecules(root, oracle_session):
    """``{filename: [chython 2 record, ...]}`` -- the oracle side, read in one call.

    A ``None`` hydrogen count is a real answer ("unknown") to be matched, not a gap to be filled.  A
    file chython 2 cannot read at all is dropped rather than raised.
    """
    paths = {name: root / 'test' / name for name in CORPUS
             if (root / 'test' / name).exists()}
    return {name: records
            for name, records in oracle_session.read_sdf(paths, tolerant=True).items()
            if records is not None}
