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
"""The read-only facades share one signature: `f(text, *, log=None) -> list`.

XYZ, PDB, mmCIF and MOL2 have no writer -- these formats state coordinates, or state bonds without
stating enough to check them -- so one direction each.  Parametrized because the point is that they agree.
"""

from pytest import mark, raises

from ..mol2 import mol2
from ..pdb import mmcif, pdb
from ..xyz import xyz


_XYZ = '1\nwater\nO 0.0 0.0 0.0\n'
_PDB = ('ATOM      1  O   HOH A   1       0.000   0.000   0.000  1.00 10.00           O  \n'
        'END\n')
_MMCIF = '\n'.join(['data_T', 'loop_',
                    '_atom_site.group_PDB', '_atom_site.id', '_atom_site.type_symbol',
                    '_atom_site.label_atom_id', '_atom_site.label_comp_id',
                    '_atom_site.label_asym_id', '_atom_site.Cartn_x', '_atom_site.Cartn_y',
                    '_atom_site.Cartn_z',
                    'HETATM 1 O O HOH A 0.000 0.000 0.000', ''])
_MOL2 = '\n'.join(['@<TRIPOS>MOLECULE', 'water', ' 1 0 1 0 0', 'SMALL', 'NO_CHARGES', '', '',
                   '@<TRIPOS>ATOM', '      1 O          0.0000    0.0000    0.0000 O.3   1 HOH  0.0000',
                   ''])

_FACADES = [(xyz, _XYZ), (pdb, _PDB), (mmcif, _MMCIF), (mol2, _MOL2)]
_IDS = ['xyz', 'pdb', 'mmcif', 'mol2']


@mark.parametrize('facade,text', _FACADES, ids=_IDS)
def test_a_facade_answers_a_list_of_one_for_a_one_record_document(facade, text):
    assert len(facade(text)) == 1


@mark.parametrize('facade,text', _FACADES, ids=_IDS)
def test_a_facade_takes_a_log_list(facade, text):
    log = []
    facade(text, log=log)
    assert isinstance(log, list)


@mark.parametrize('facade,_', _FACADES, ids=_IDS)
def test_bytes_are_refused_by_name(facade, _):
    """Named, because "'int' object has no attribute 'rstrip'" tells a caller nothing."""
    with raises(TypeError, match='str'):
        facade(b'anything')


@mark.parametrize('facade,_', _FACADES, ids=_IDS)
def test_a_path_shaped_string_is_read_as_text_not_opened(facade, _):
    """These take text.  `read_mol2`/`read_pdb` take files; a facade that guessed would be neither."""
    log = []
    assert facade('molecules.mol2', log=log) == [] or log


@mark.parametrize('facade,_', _FACADES, ids=_IDS)
def test_an_empty_document_is_an_empty_list_and_not_an_error(facade, _):
    log = []
    assert facade('', log=log) == []


def test_mol2_returns_every_record_not_just_the_first():
    assert len(mol2(_MOL2 + _MOL2)) == 2


def test_mol2_substitutes_a_failed_record_rather_than_raising():
    """An ATOM line short of its six required fields is the one thing MOL2 refuses to guess at.

    A record with no ATOM block at all is NOT that: it is an empty molecule, and the reader stores it.
    """
    from ..mol2 import FailedRecord

    broken = '\n'.join(['@<TRIPOS>MOLECULE', 'broken', ' 1 0 1 0 0', 'SMALL', 'NO_CHARGES', '', '',
                        '@<TRIPOS>ATOM', '      1 O   0.0000', ''])
    records = mol2(_MOL2 + broken)
    assert len(records) == 2
    assert isinstance(records[1], FailedRecord)
