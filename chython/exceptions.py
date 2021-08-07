# -*- coding: utf-8 -*-
#
#  Copyright 2017-2021 Ramil Nugmanov <nougmanoff@protonmail.com>
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
    bad files parsing
    """


class MappingError(ValueError):
    """
    bad files parsing
    """


class AtomNotFound(KeyError):
    """
    bad atom number
    """


class NotChiral(KeyError):
    """
    bad atom number
    """


class IsChiral(KeyError):
    """
    bad atom number
    """


class InvalidAromaticRing(ValueError):
    """
    aromatic ring has impossible kekule structure
    """


class IsConnectedAtom(Exception):
    """
    atom already attached to graph
    """


class IsNotConnectedAtom(Exception):
    """
    atom already attached to graph
    """


class ValenceError(Exception):
    """
    atom has error in valence
    """


class IncorrectSmiles(ValueError):
    """
    SMILES string invalid
    """


class ImplementationError(Exception):
    """
    Algorithm has errors. Please send example of structure to author for analyze.
    """
