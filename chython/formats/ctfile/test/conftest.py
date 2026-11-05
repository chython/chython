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

from .._sdf import split_records


#: Files with more than a handful of records, used for the sweeps.  Named rather than globbed so a
#: new file appearing in `test/` cannot silently change what a passing suite means.
CORPUS = ('arenes.sdf', 'cycle.sdf', 'isomorphism.sdf', 'mcs.sdf', 'peptide.sdf',
          'standardize.sdf', 'stereo.sdf', 'implicit.sdf', 'depict.sdf', 'hbonds.sdf',
          'morgan_ruiner.sdf')


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
