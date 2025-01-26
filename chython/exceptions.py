# -*- coding: utf-8 -*-
#
#  Copyright 2017-2024 Ramil Nugmanov <nougmanoff@protonmail.com>
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


class EmptyReaction(ValueError):
    """
    Reaction without molecules
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


class BufferOverflow(BufferError):
    """
    Parser buffer overflow
    """


class InvalidV2000(ValueError):
    """
    Invalid V2000
    """


class InvalidCharge(ValueError):
    """
    Invalid MDL charge
    """


class InvalidMolBlock(ValueError):
    """
    Invalid MDL MOL
    """
