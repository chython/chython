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
IUPAC names, both directions, over two vendors: OPSIN (Java, via `_java.py:get_opsin`) parses a name
and openclatura writes one.

`from_iupac` parses OPSIN's SMILES with core's `read_smiles`; `to_iupac` exports through
`interop.rdkit` and forwards `log` to it.  openclatura requires Python >= 3.11.
"""
from ..core import LogRecord
from ..exceptions import ToolkitError, UnconvertibleType
from ._records import deliver


def to_iupac(mol, /, *, log=None):
    """
    Name a chython container, returning `None` when openclatura cannot name it.

    Runs through `interop.rdkit`, since openclatura names an RDKit molecule, so a loss in that export
    is a loss here and `log` is forwarded to it.  openclatura is optional and requires Python >= 3.11;
    its absence raises `ImportError`, not `DirectionNotImplemented` -- the direction exists.

    :param log: optional list receiving one human-readable line per reportable loss.
    """
    # openclatura is checked before the rdkit conversion: it is what makes this direction possible, so
    # its absence is the actionable answer even when the rdkit half is also missing.
    try:
        from openclatura import name_rdkit_mol
    except ImportError:
        raise ImportError('openclatura is not installed. '
                          '`pip install openclatura` (requires Python >= 3.11)')
    # forward log so any loss in the rdkit export surfaces as a loss of this conversion too.
    # `keep_mapping=False`: a name has no atom-atom mapping, so a mapped molecule is named as the
    # molecule it is rather than handed to openclatura with labels it has no use for.
    from . import rdkit as _rdkit
    return name_rdkit_mol(_rdkit(mol, log=log, keep_mapping=False)) or None


def from_iupac(name, /, *, log=None):
    """
    Parse an IUPAC name into a chython 3 molecule, using OPSIN.

    Returns a `chython.core.MoleculeContainer`: it calls core's `read_smiles` directly, so the result
    is a V3 molecule regardless of which SMILES facade `import chython` exposes.  Failures raise
    `ToolkitError` with OPSIN's own message preserved, as does an unreadable SMILES from OPSIN.
    A non-string argument raises `UnconvertibleType`.

    EVERY RECORD LANDS ON THE RETURNED MOLECULE'S `.log`, in stage `'interop'`, with nothing passed in.
    THE SMILES OPSIN PRODUCED IS ONE OF THEM: the molecule is OPSIN's reading of the name, the string
    chython actually parsed is not recoverable from the result, and every other record of this import
    is about that string rather than about the name.  `log`, when given, receives a copy of the same
    records.

    :param log: optional list receiving a copy of the records put on the molecule.
    """
    if not isinstance(name, str):
        raise UnconvertibleType(
            'iupac reads string IUPAC names only; '
            f'got {type(name).__name__!r} -- to export a molecule to a name, pass a container'
        )
    from ._java import get_opsin
    from ..core._core import IncorrectSmiles, read_smiles

    result = get_opsin().parseChemicalName(name)
    if str(result.getStatus()) == 'FAILURE':
        raise ToolkitError(f'OPSIN failed to parse {name!r}: {result.getMessage()}')
    smiles_str = str(result.getSmiles())
    records = [LogRecord('iupac:parsed-by-opsin', (),
                         f'name {name!r} was read as OPSIN\'s SMILES {smiles_str!r}')]
    try:
        # the reader's own lines are about OPSIN's SMILES, so they belong to this conversion.
        mol = read_smiles(smiles_str, records)
    except IncorrectSmiles as e:
        raise ToolkitError(
            f'OPSIN returned unreadable SMILES {smiles_str!r} for name {name!r}: {e}'
        ) from e
    deliver(mol, records, log)
    return mol
