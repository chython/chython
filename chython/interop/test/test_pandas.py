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
"""
`patch_pandas`.

The patch is global -- it rebinds a name inside pandas -- so the `restored` fixture is not optional:
every test here must put the original back and reset the idempotence flag.
"""
from pytest import fixture, importorskip


importorskip('pandas', reason='patch_pandas patches pandas, which is an optional dependency')


@fixture
def restored():
    """Undo the patch, whatever it did, and let the next test start from an unpatched pandas."""
    from pandas.io.formats import printing

    from .. import _pandas

    original, was = printing.is_sequence, _pandas._patched
    _pandas._patched = False
    try:
        yield printing
    finally:
        printing.is_sequence = original
        _pandas._patched = was


def _molecules():
    """The molecules the patch has to cover; the label names which one a failure classified."""
    from ...core import read_smiles

    return [('core', read_smiles('CCO'))]


def test_the_defect_reproduces_without_the_patch(restored):
    """A molecule iterates and has a length, so pandas calls it a sequence and prints its atoms."""
    for name, mol in _molecules():
        assert restored.is_sequence(mol), \
            f'{name}: pandas no longer calls a molecule a sequence, so re-derive this patch first'


def test_a_patched_pandas_treats_both_generations_as_one_value(restored):
    from .. import patch_pandas

    patch_pandas()
    for name, mol in _molecules():
        assert not restored.is_sequence(mol), name


def test_a_real_sequence_is_still_a_sequence(restored):
    """The patch narrows one answer; it must not answer False for everything."""
    from .. import patch_pandas

    patch_pandas()
    assert restored.is_sequence([1, 2, 3])
    assert restored.is_sequence((1, 2))
    assert not restored.is_sequence('CCO')      # pandas excludes strings, and still must


def test_patching_twice_does_not_wrap_the_wrapper(restored):
    """It is documented as a notebook's first line, and a notebook's first cell gets re-run."""
    from .. import patch_pandas

    patch_pandas()
    once = restored.is_sequence
    patch_pandas()
    assert restored.is_sequence is once
    assert restored.is_sequence([1, 2, 3])      # and the behaviour survives the second call


def test_a_dataframe_of_molecules_renders_the_molecules(restored):
    """The end the patch exists for, asserted end to end rather than through the predicate alone."""
    from pandas import DataFrame

    from .. import patch_pandas

    mols = [mol for _, mol in _molecules()]
    # `(1, 2, 3)` and not `[1, 2, 3]`: pandas formats a cell it calls a sequence as a tuple, so what a
    # reader sees in the cell is the atom numbers.
    assert '(1, 2, 3)' in DataFrame({'mol': mols}).to_string(), \
        'expected the unpatched frame to print atom numbers; the patch is being tested against nothing'

    patch_pandas()
    rendered = DataFrame({'mol': mols}).to_string()
    assert '(1, 2, 3)' not in rendered
    assert 'C(C)O' in rendered or 'CCO' in rendered, rendered
