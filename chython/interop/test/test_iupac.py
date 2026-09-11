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
Tests for `interop._iupac`.

`from_iupac` uses OPSIN (Java/JPype); its tests skip when the jar or jpype is missing.  `to_iupac`
needs openclatura, which requires Python >= 3.11, so those tests skip on 3.10.
"""
import pytest

from .conftest import absent_openclatura, requires_openclatura


def test_from_iupac_ethanol(opsin):
    """A well-known name round-trips through OPSIN and core's reader."""
    from chython.interop._iupac import from_iupac

    mol = from_iupac('ethanol')
    assert str(mol) is not None
    # CCO and C(C)O are both valid canonical forms, so only the atoms are checked.
    assert len(mol) == 3  # C, C, O


def test_from_iupac_benzene(opsin):
    """Benzene parses and returns a molecule with the right atom count."""
    from chython.interop._iupac import from_iupac

    mol = from_iupac('benzene')
    assert len(mol) == 6  # six carbons


def test_from_iupac_failure_raises_toolkit_error(opsin):
    """A name OPSIN cannot parse raises `ToolkitError`, not a bare `ValueError`."""
    from chython.interop._iupac import from_iupac
    from chython.exceptions import ToolkitError

    with pytest.raises(ToolkitError):
        from_iupac('xyzzy not a chemical name abc123')


def test_from_iupac_error_message_contains_opsin_output(opsin):
    """
    OPSIN's own message is preserved inside the `ToolkitError`.

    It names what OPSIN objected to, which chython does not know independently.
    """
    from chython.interop._iupac import from_iupac
    from chython.exceptions import ToolkitError

    bad_name = 'not a chemical xyz987'
    with pytest.raises(ToolkitError) as exc_info:
        from_iupac(bad_name)
    # OPSIN's message echoes the name it could not parse.
    assert bad_name in str(exc_info.value)


def test_from_iupac_returns_v3_molecule(opsin):
    """
    `from_iupac` returns the core's `MoleculeContainer`, parsed by the core's own SMILES reader.

    OPSIN hands over a SMILES string, so a regression here returns the string.
    """
    from chython.interop._iupac import from_iupac
    from chython.core import MoleculeContainer as V3Molecule

    mol = from_iupac('acetic acid')
    assert isinstance(mol, V3Molecule), f'expected a core MoleculeContainer, got {type(mol)}'


def test_from_iupac_stereo_preserved(opsin):
    """
    Tetrahedral stereo in OPSIN's SMILES output survives into the V3 molecule.

    OPSIN returns `[C@H]`/`[C@@H]` tokens; core's `read_smiles` preserves them.
    """
    from chython.interop._iupac import from_iupac

    mol = from_iupac('(R)-lactic acid')
    smiles = str(mol)
    assert '@' in smiles, (
        f'stereo center expected in SMILES for (R)-lactic acid, got {smiles!r}'
    )


def test_from_iupac_non_string_raises_unconvertible_type():
    """
    Passing a non-string to `from_iupac` raises `UnconvertibleType`, never `DirectionNotImplemented`.

    The latter is a `NotImplementedError` and would say "right argument, missing half" instead.  No
    `opsin` fixture: the check fires before any Java call.
    """
    from chython.interop._iupac import from_iupac
    from chython.exceptions import DirectionNotImplemented, UnconvertibleType

    for bad in (object(), 42, None, b'ethanol', [1]):
        with pytest.raises(UnconvertibleType):
            from_iupac(bad)

    with pytest.raises(UnconvertibleType) as e:
        from_iupac(object())
    assert not isinstance(e.value, DirectionNotImplemented)


@absent_openclatura
def test_to_iupac_raises_import_error_without_openclatura(opsin):
    """
    When openclatura is absent, `to_iupac` raises `ImportError` naming the install command.

    `ImportError` and not `DirectionNotImplemented`: the direction exists, the dependency does not.
    Exact because openclatura is checked before the rdkit conversion.
    """
    from chython.interop._iupac import from_iupac, to_iupac

    mol = from_iupac('ethanol')

    with pytest.raises(ImportError, match='openclatura is not installed'):
        to_iupac(mol)


@requires_openclatura
def test_to_iupac_ethanol():
    """`to_iupac` names ethanol correctly."""
    from chython.interop._iupac import from_iupac, to_iupac

    mol = from_iupac('ethanol')
    name = to_iupac(mol)
    assert name == 'ethanol'


# `to_iupac` returning None for a structure openclatura declines is deliberately not tested: which
# structures it declines changes between its releases, so the test would pin that library's coverage.
