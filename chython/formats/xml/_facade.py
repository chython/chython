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
"""One callable per XML dialect: `mrv()` and `cml()`.

Both directions in one name, as `mol()` and `rxn()` do.  `read_mrv`/`write_mrv` and their CML twins stay
exported -- a caller who wants to state the direction can, and `parse_xml_document`'s sniffing door is
unaffected.
"""
from ._cml import read_cml, write_cml
from ._mrv import read_mrv, write_mrv
from ._dialect import Record
from ...core import MoleculeContainer
from ...core.reaction import ReactionContainer


__all__ = ['cml', 'mrv']

#: What counts as "the caller handed me structures" rather than "the caller handed me a document".  A
#: `list`/`tuple` only: a generator would have to be consumed to find out, and consuming it to decide is
#: not a decision a facade gets to make.
_WRITABLE = (MoleculeContainer, Record)


def _is_write(data):
    if isinstance(data, ReactionContainer):
        raise TypeError('this dialect is modelled for molecules; a reaction has no writer here')
    if isinstance(data, _WRITABLE):
        return True
    if isinstance(data, (list, tuple)):
        if not data:
            raise ValueError('nothing to write and nothing to read')
        if all(isinstance(x, _WRITABLE) for x in data):
            return True
        raise TypeError('a sequence must hold molecules only')
    return False


def mrv(data, *, log=None, title=None, indent='  '):
    """Molecules from an MRV document, or an MRV document from molecules.

    :param data: MRV XML text, or a molecule, or a list of molecules.
    :param log: a list to append damage reports to.
    :param title: export only; the document title.
    :param indent: export only; one level of indentation, `''` for one line.

    A read answers a list, however many molecules the document holds.  A record chython cannot build is
    logged, not raised -- input is garbage by default.
    """
    if _is_write(data):
        return write_mrv(data, log=log, title=title, indent=indent)
    return read_mrv(data, log=log)


def cml(data, *, log=None, title=None, indent='  '):
    """Molecules from a CML document, or a CML document from molecules.

    Arguments are :func:`mrv`'s.  Data fields ride in `<propertyList>`.
    """
    if _is_write(data):
        return write_cml(data, log=log, title=title, indent=indent)
    return read_cml(data, log=log)
