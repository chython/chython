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
"""`mrv()` and `cml()`: one callable per XML dialect, both directions.

Direction is decided by the argument's type, as in `mol()`.  A read always answers a LIST -- an MRV
document states any number of molecules, and one is not a special case.
"""

from pytest import raises

from .._facade import cml, mrv
from ...ctfile import mol
from chython.core import ReactionContainer, read_smiles


def test_mrv_writes_a_molecule_and_reads_it_back():
    m = read_smiles('c1ccccc1')
    text = mrv(m)
    assert text.startswith('<?xml') or text.lstrip().startswith('<')
    back, = mrv(text)
    assert str(back) == str(m)


def test_cml_writes_a_molecule_and_reads_it_back():
    m = read_smiles('CCO')
    back, = cml(cml(m))
    assert str(back) == str(m)


def test_a_list_of_molecules_is_a_write_and_answers_a_list_of_the_same_length():
    ms = [read_smiles('CCO'), read_smiles('c1ccccc1')]
    assert len(mrv(mrv(ms))) == 2
    assert len(cml(cml(ms))) == 2


def test_a_title_is_stated_on_the_write_side():
    m = read_smiles('CCO')
    assert 'ethanol' in mrv(m, title='ethanol')          # `title` is `str` after Stream 1


def test_data_fields_survive_a_cml_round_trip():
    """CML `<propertyList>` is `mol.meta` from Stream 1 on, so this is the round trip `mol()` gets."""
    m = mol('\n'.join(['e', '', '', '  1  0  0  0  0  0            999 V2000',
                       '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0',
                       'M  END']))
    m.meta['ACTIVITY'] = '5.0'
    back, = cml(cml(m))
    assert back.meta['ACTIVITY'] == '5.0'


def test_a_reaction_is_refused_by_name():
    """Neither dialect is modelled for reactions here -- an MRV `<reaction>` read is logged
    `unsupported:`, so a write that silently produced a molecule document would be worse than a refusal.
    """
    rxn = ReactionContainer(reactants=(read_smiles('CCO'),), products=(read_smiles('CC=O'),))
    with raises(TypeError, match='reaction'):
        mrv(rxn)
    with raises(TypeError, match='reaction'):
        cml(rxn)


def test_a_generator_is_read_as_text_and_refused_rather_than_written():
    """A generator cannot be told from an unread stream without consuming it, so it is not a write."""
    with raises(TypeError):
        mrv(read_smiles(x) for x in ('CCO',))


def test_an_empty_document_is_an_empty_list_and_not_an_error():
    """Nothing to read is not damage: `<cml/>` says no molecules and the reader agrees, silently."""
    log = []
    assert mrv('<cml></cml>', log=log) == []
    assert not log


def test_damage_is_logged_not_raised():
    """A `<reaction>` is the construct this dialect does not model, and it costs a log line, not a raise."""
    log = []
    assert mrv('<cml><MDocument><MChemicalStruct><reaction/>'
               '</MChemicalStruct></MDocument></cml>', log=log) == []
    assert any('unsupported: ' in str(x) for x in log), log
