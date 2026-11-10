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
"""HELM and BILN: biopolymer notations against an attachable monomer library.

One intermediate, the monomer graph (:mod:`._graph`), sits between the notations and the molecule:
:mod:`._helm` and :mod:`._biln` parse into it, :mod:`._assemble` builds the molecule from it,
:mod:`._decompose` cuts a molecule back into it and :mod:`._write` spells it canonically.  The monomers
come from :mod:`._library`: the shipped public set plus whatever :func:`set_monomer_library` attached.
"""

from ._assemble import assemble
from ._biln import parse_biln
from ._decompose import decompose
from ._errors import BilnError, HelmError, StereoMismatch
from ._helm import parse_helm
from ._library import get_monomer_library, set_monomer_library
from ._write import write_biln, write_helm
from ..core import MoleculeContainer


__all__ = ['helm', 'biln', 'HelmError', 'BilnError', 'StereoMismatch', 'get_monomer_library',
           'set_monomer_library']


def helm(data, *, log=None):
    """A molecule from a HELM string, or the canonical HELM string of a molecule.

    :param log: a list to append records to; a read also lands them on the molecule's own log.

    Unknown monomers read as markers carrying their symbol; a piece matching no monomer writes as an
    inline SMILES monomer; a piece matching a monomer by constitution only raises
    :class:`StereoMismatch`.  Each is logged where it does not raise.
    """
    return _convert(data, log, parse_helm, write_helm)


def biln(data, *, log=None):
    """A molecule from a BILN string, or the canonical BILN string of a molecule.  As :func:`helm`; a
    molecule needing an inline SMILES monomer or a non-peptide chain raises :class:`BilnError`."""
    return _convert(data, log, parse_biln, write_biln)


def _convert(data, log, parse, write):
    log = [] if log is None else log
    library = get_monomer_library()
    if isinstance(data, MoleculeContainer):
        return write(decompose(data, library, log))
    if not isinstance(data, str):
        raise TypeError(f'expected a notation string or a MoleculeContainer, got {type(data).__name__}')
    start = len(log)
    mol = assemble(parse(data), library, log)
    if len(log) > start:
        mol.log.absorb('read', log[start:])
    return mol
