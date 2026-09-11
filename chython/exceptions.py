# -*- coding: utf-8 -*-
#
#  Copyright 2017-2026 Ramil Nugmanov <nougmanoff@protonmail.com>
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


class DirectionNotImplemented(NotImplementedError):
    """
    The format is known; this direction of it is not built.

    Every format callable is bidirectional -- `smiles(str)` parses and `smiles(mol)` writes -- and
    the direction is chosen by the argument.  Some halves do not exist: nothing writes an xyz matrix,
    and chython reads no CDPKit molecule.  Those raise this rather than `TypeError`, because the two
    say different things to a caller.  `TypeError` says "you passed the wrong thing"; this says "you
    passed the right thing and that half is not built yet", which is the difference between a bug in
    the caller and a gap in the library.

    A `NotImplementedError` so that a caller probing for capability can catch the standard exception
    without importing chython's own.
    """


class UnconvertibleType(TypeError):
    """
    Object is neither a chython container nor a type this converter reads.

    The direction of a format callable is decided by testing for a chython container first, because
    that is the only test available without importing a third-party toolkit.  Anything else is
    offered to the import direction, and this is what that direction raises when the object is not
    something it knows how to read.  The message names both types it accepts, so the answer to "what
    was I supposed to pass?" is in the traceback rather than in the documentation.
    """


class ToolkitError(RuntimeError):
    """
    A third-party toolkit refused a conversion chython asked it to make.

    Distinct from a reportable loss, which is the normal case: if a molecule can be converted at all
    it is converted and what was dropped goes to the caller's `log`.  This is the other case -- the
    toolkit itself rejected the structure or failed to sanitize it, so there is no molecule to hand
    back and no partial answer worth inventing.  The toolkit's own message is preserved, because it
    knows what it objected to and chython does not.
    """
