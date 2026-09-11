# -*- coding: utf-8 -*-
#
#  Copyright 2026 Ramil Nugmanov <nougmanoff@protonmail.com>
#
#  This file is part of chython.
#
#  chython is free software; you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with this program; if not, see <https://www.gnu.org/licenses/>.
#
"""`title` is `str`, and the byte that is not UTF-8 still comes back.

The arena is unchanged -- the name line is `bytes` at handle 0 of the blob.  What changed is the type at
the Python surface, and `errors='surrogateescape'` is what lets ONE name carry both promises: text for
every caller, and the exact byte for the file it came from.
"""
from pytest import raises
from chython.core import MoleculeContainer, ReactionContainer


#: Not valid UTF-8, and with a NUL and a newline in it: `test_sgroups.RAW`, which the blob stores exactly.
RAW = b'\xff\xfe\x80 caf\xe9 latin-1, a NUL:\x00 and a newline:\n'


def _mol():
    mol = MoleculeContainer()
    with mol.edit() as e:
        e.add_atom('C')
    return mol


def test_a_title_is_str():
    mol = _mol()
    mol.set_title('ethanol')
    assert mol.title == 'ethanol' and isinstance(mol.title, str)
    assert MoleculeContainer().title == '', 'absent and empty are the same answer for a title'


def test_an_undecodable_byte_survives_as_a_surrogate_and_re_encodes_exactly():
    """THE FIDELITY PROMISE, now a promise and not a type."""
    mol = _mol()
    mol.set_title(RAW)
    assert isinstance(mol.title, str)
    assert mol.title.encode('utf8', 'surrogateescape') == RAW
    assert MoleculeContainer.from_bytes(mol.to_bytes()).title == mol.title


def test_a_str_title_comes_back_as_itself():
    mol = _mol()
    mol.set_title('café')
    assert mol.title == 'café'


def test_bytes_and_str_titles_do_not_collide():
    """Two different name lines stay two different name lines, which is what the handler buys."""
    a, b = _mol(), _mol()
    a.set_title(b'caf\xe9')
    b.set_title('café')
    assert a.title != b.title
    assert a.title.encode('utf8', 'surrogateescape') == b'caf\xe9'
    assert b.title.encode('utf8', 'surrogateescape') == b'caf\xc3\xa9'


def test_set_title_takes_str_bytes_and_a_buffer_and_refuses_the_rest():
    mol = _mol()
    mol.set_title(bytearray(b'buf'))
    assert mol.title == 'buf'
    mol.set_title(memoryview(b'view'))
    assert mol.title == 'view'
    with raises(TypeError):
        mol.set_title(42)


def test_set_sgroups_does_not_disturb_the_title():
    """It re-serialises the blob, so it must move the title as bytes and never through the property."""
    mol = _mol()
    mol.set_title(RAW)
    mol.set_sgroups([{'type': b'DAT', 'atoms': (mol.atom_numbers[0],)}])
    assert mol.title.encode('utf8', 'surrogateescape') == RAW and len(mol.sgroups) == 1


def test_the_reaction_makes_the_same_promise_with_the_same_spelling():
    rxn = ReactionContainer([_mol()], [_mol()], title=RAW)
    assert isinstance(rxn.title, str) and rxn.title.encode('utf8', 'surrogateescape') == RAW
    rxn.set_title('named')
    assert rxn.title == 'named'
    assert ReactionContainer([_mol()], [_mol()]).title == ''
    with raises(TypeError):
        ReactionContainer([_mol()], [_mol()], title=42)


def test_the_cost_is_real_and_narrow():
    """A surrogate-escaped title is not UTF-8-encodable WITHOUT the handler.

    Asserted rather than hidden: it is the documented trade, and a test is where a trade stops being a
    surprise.  `json.dumps` on such a title raises, and that needs an undecodable byte in the file.
    """
    mol = _mol()
    mol.set_title(b'caf\xe9')
    with raises(UnicodeEncodeError):
        mol.title.encode('utf8')
