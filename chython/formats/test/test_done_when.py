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
"""The five sentences one `Meta` and one `Log` are finished by, each as one assertion.

Here and not in `core/test/` because four of the five are about a reader: the claim is that a container
comes back from a file carrying what the file said, which is a statement neither layer can make alone.
Imports reach `chython.core` and `chython.formats.*` directly -- `test_isolation.py` is the ratchet.
"""
from chython.core import LogRecord, ReactionContainer
from chython.formats.ctfile import RDFRead, RDFWrite, SDFRead, SDFWrite, mol, parse_record
from chython.formats.xml import read_cml, write_cml


_BUTANE = ['butane', '', '',
           '  4  3  0  0  0  0            999 V2000',
           '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
           '    0.8660    0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
           '    1.7320    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
           '    2.5980    0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
           '  1  2  1  0  0  0  0',
           '  2  3  1  0  0  0  0',
           '  3  4  1  0  0  0  0',
           'M  END']


def test_sdf_data_fields_land_in_meta_and_come_back_out(tmp_path):
    text = '\n'.join(_BUTANE + ['>  <NAME>', 'butane', '', '$$$$']) + '\n'
    (tmp_path / 'in.sdf').write_text(text)
    with SDFRead(tmp_path / 'in.sdf') as f:
        molecule = next(iter(f))
    assert molecule.meta == {'NAME': 'butane'}
    with SDFWrite(tmp_path / 'out.sdf') as w:
        w.write(molecule)
    with SDFRead(tmp_path / 'out.sdf') as f:
        assert next(iter(f)).meta == {'NAME': 'butane'}


def test_rdf_fields_land_in_meta_on_both_kinds_of_record(tmp_path):
    molecule = parse_record(_BUTANE)
    molecule.meta['K'] = 'v'
    reaction = ReactionContainer([molecule.copy()], [molecule.copy()])
    reaction.meta['K'] = 'v'
    with RDFWrite(tmp_path / 'a.rdf') as w:
        w.write(molecule)
        w.write(reaction)
    with RDFRead(tmp_path / 'a.rdf') as f:
        assert [x.meta for x in f] == [{'K': 'v'}, {'K': 'v'}]


def test_a_cml_property_list_lands_in_meta():
    document = ('<molecule title="x"><propertyList><property dictRef="k">'
                '<scalar>v</scalar></property></propertyList>'
                '<atomArray><atom id="a1" elementType="C"/></atomArray></molecule>')
    assert read_cml(document)[0].meta == {'k': 'v'}
    assert 'dictRef="k"' in write_cml(read_cml(document))


def test_no_wrapper_type_exposes_a_second_meta():
    import chython.formats.ctfile as ctfile
    for name in ('CtfileRecord', 'ReactionRecord', 'FieldsView', 'DataField'):
        assert not hasattr(ctfile, name), f'{name} is back'


def test_the_log_shows_what_the_reader_repaired_with_no_log_passed():
    """No `log=` anywhere in the call.  This is the sentence the whole Log half is for."""
    molecule = mol('\n'.join(_BUTANE[:3] + ['  4  3  0  0  0  0            999     '] + _BUTANE[4:]))
    assert molecule.log and all(isinstance(x, LogRecord) for x in molecule.log)
    assert molecule.log.repaired(), 'an unstamped counts line is read as V2000, which is a repair'
