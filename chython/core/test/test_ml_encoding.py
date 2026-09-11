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
"""`TensorEncoding`: the knobs, their defaults, and the domains they refuse."""
from pickle import dumps, loads

from pytest import mark, raises

from chython.core import TensorEncoding, read_smiles


def test_the_default_encoding_is_identity_everywhere_but_the_unknown_hydrogen_count():
    """A bare encoding must not move a number, so a physical readout needs no argument."""
    enc = TensorEncoding()
    assert enc.element_shift == 0
    assert enc.hydrogen_shift == 0
    assert enc.neighbor_shift == 0
    assert enc.distance_shift == 0
    assert enc.max_distance == 0, 'zero is off, not a clamp to zero'
    assert enc.max_neighbors == 0
    assert enc.width == 0
    assert enc.pad == 0
    assert enc.pad_diagonal == 0
    assert enc.disconnected == -1, 'the physical answer for a pair with no path'
    assert enc.unknown_h == 0, 'the one non-identity default; unknown_h=15 keeps H_UNKNOWN'
    assert enc.vocabulary is None
    assert enc.unknown == -1


def test_every_knob_reads_back_as_given():
    enc = TensorEncoding(element_shift=2, hydrogen_shift=1, neighbor_shift=2, distance_shift=2,
                         disconnected=1, unknown_h=15, max_distance=10, max_neighbors=14,
                         width=65, pad=0, pad_diagonal=1, vocabulary={(6, 3, 1, 3, 1): 42},
                         unknown=999)
    assert enc.element_shift == 2
    assert enc.hydrogen_shift == 1
    assert enc.neighbor_shift == 2
    assert enc.distance_shift == 2
    assert enc.disconnected == 1
    assert enc.unknown_h == 15
    assert enc.max_distance == 10
    assert enc.max_neighbors == 14
    assert enc.width == 65
    assert enc.pad == 0
    assert enc.pad_diagonal == 1
    assert enc.vocabulary == {(6, 3, 1, 3, 1): 42}
    assert enc.unknown == 999


def test_the_fields_are_readonly_so_a_compiled_table_cannot_go_stale():
    """The vocabulary is compiled at construction; a reassigned knob would not reach the table."""
    enc = TensorEncoding(element_shift=2)
    with raises(AttributeError):
        enc.element_shift = 3


def test_an_encoding_crosses_a_dataloader_worker():
    """A DataLoader pickles its dataset, and the encoding travels with it.

    All thirteen fields carry a DISTINCT non-default value so a transposition of any two in
    `__reduce__` or `_rebuild_tensor_encoding` fails rather than silently passing.
    """
    enc = loads(dumps(TensorEncoding(element_shift=2, hydrogen_shift=3, neighbor_shift=4,
                                     distance_shift=5, disconnected=6, unknown_h=7,
                                     max_distance=8, max_neighbors=9, width=10, pad=11,
                                     pad_diagonal=12, vocabulary={(6, 3, 1, 3, 1): 42},
                                     unknown=13)))
    assert enc.element_shift == 2
    assert enc.hydrogen_shift == 3
    assert enc.neighbor_shift == 4
    assert enc.distance_shift == 5
    assert enc.disconnected == 6
    assert enc.unknown_h == 7
    assert enc.max_distance == 8
    assert enc.max_neighbors == 9
    assert enc.width == 10
    assert enc.pad == 11
    assert enc.pad_diagonal == 12
    assert enc.vocabulary == {(6, 3, 1, 3, 1): 42}
    assert enc.unknown == 13


@mark.parametrize('kwargs, message', [
    ({'unknown_h': 16}, 'unknown_h'),
    ({'unknown_h': -1}, 'unknown_h'),
    ({'width': -1}, 'width'),
    ({'max_distance': -1}, 'max_distance'),
    ({'max_neighbors': -1}, 'max_neighbors'),
    ({'max_neighbors': 256}, 'max_neighbors'),
])
def test_a_value_outside_its_declared_domain_is_refused_by_name(kwargs, message):
    """A width is not a bound (RULES.md 6.1): the domain is stated here and nowhere else."""
    with raises(ValueError, match=message):
        TensorEncoding(**kwargs)


def test_a_vocabulary_key_is_five_numbers_and_a_wrong_shape_says_so():
    with raises(ValueError, match='five'):
        TensorEncoding(vocabulary={(6, 3, 1): 42})


def test_the_repr_shows_only_what_was_set():
    """Twelve zeros in a repr hide the one knob that is not zero."""
    assert repr(TensorEncoding()) == 'TensorEncoding()'
    assert repr(TensorEncoding(element_shift=2, width=65)) == 'TensorEncoding(element_shift=2, width=65)'


def test_the_all_zero_key_is_refused_because_it_equals_the_empty_slot_marker():
    """ML_VOCAB_KEY(0,0,0,0,0) == 0 == ML_KEY_EMPTY; the probe uses 0 as its terminator."""
    with raises(ValueError, match='empty-slot marker'):
        TensorEncoding(vocabulary={(0, 0, 0, 0, 0): 1})


def test_a_vocabulary_value_outside_int32_is_refused():
    with raises(ValueError, match='token'):
        TensorEncoding(vocabulary={(6, 3, 1, 3, 1): 2 ** 31})


def test_a_vocabulary_key_component_outside_its_column_is_refused_by_name():
    """The key packs into 31 bits; a value that does not fit would collide with another key."""
    with raises(ValueError, match='element'):
        TensorEncoding(vocabulary={(200, 3, 1, 3, 1): 1})
    with raises(ValueError, match='hydrogen'):
        TensorEncoding(vocabulary={(6, 16, 1, 3, 1): 1})
    with raises(ValueError, match='neighbor'):
        TensorEncoding(vocabulary={(6, 3, 300, 3, 1): 1})


def test_a_compiled_vocabulary_survives_a_pickle_round_trip():
    enc = loads(dumps(TensorEncoding(vocabulary={(6, 3, 1, 3, 1): 42, (6, 2, 2, 2, 2): 43,
                                                 (8, 1, 1, 1, 1): 44}, unknown=999)))
    assert enc.vocabulary == {(6, 3, 1, 3, 1): 42, (6, 2, 2, 2, 2): 43, (8, 1, 1, 1, 1): 44}
    assert enc.unknown == 999
    # a __reduce__ that carried the dict but never recompiled would also pass the two asserts above;
    # producing tokens confirms the compiled table was rebuilt on the far side.
    assert read_smiles('CCO').state_view(enc).tokens.tolist() == [42, 43, 44]
