# -*- coding: utf-8 -*-
#
#  Copyright 2017-2022 Ramil Nugmanov <nougmanoff@protonmail.com>
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


class EmptyMolecule(ValueError):
    """
    Molecule without atoms
    """


class MappingError(ValueError):
    """
    Atom-to-Atom mapping invalid
    """


class AtomNotFound(KeyError):
    """
    Bad atom number
    """


class BondNotFound(KeyError):
    """
    Bad atoms numbers
    """


class NotChiral(KeyError):
    """
    Atom not chiral
    """


class IsChiral(KeyError):
    """
    Atom already chiral
    """


class InvalidAromaticRing(ValueError):
    """
    Aromatic ring has impossible Kekule structure
    """


class IsConnectedAtom(Exception):
    """
    Atom is already attached to graph
    """


class IsNotConnectedAtom(Exception):
    """
    Atom is not attached to graph
    """


class IsConnectedBond(Exception):
    """
    Bond is already attached to graph
    """


class IsNotConnectedBond(Exception):
    """
    Bond is not attached to graph
    """


class ValenceError(Exception):
    """
    Atom has valence error
    """


class IncorrectSmiles(ValueError):
    """
    SMILES string invalid
    """


class IncorrectSmarts(IncorrectSmiles):
    """
    SMARTS string invalid or unsupported
    """


class ImplementationError(Exception):
    """
    Algorithm has errors. Please send example of structure to author for analyze.
    """


class ParseError(ValueError):
    """
    File parsing error
    """
    __slots__ = ('number', 'position', 'log', 'meta', 'structure', 'errors')

    def __init__(self, number, position, log, meta, structure=None, errors=None):
        super().__init__()
        self.number = number
        self.position = position
        self.log = log
        self.meta = meta
        self.structure = structure
        self.errors = errors


class EmptyV2000(ValueError):
    """Empty V2000 in EMol parser"""


class ParseReactionError(ValueError):
    __slots__ = ('structure', 'errors')

    def __init__(self, structure, errors):
        self.structure = structure
        self.errors = errors
