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
Interoperation with other cheminformatics toolkits: one callable per tool, named after the tool.

The argument decides the direction -- `rdkit(mol)` exports, `rdkit(rd_mol)` imports -- by testing for
a chython container first; anything else goes to the import side, which raises `UnconvertibleType` if
it does not recognise it.  Every toolkit is an optional dependency and every import of one is lazy,
inside the function that needs it: importing `chython` must not start a JVM or load RDKit.

THE EXPORT DIRECTION IS ALSO A METHOD -- `mol.to_rdkit()`, `mol.to_indigo()`, `mol.to_openbabel()`,
`mol.to_cdk()`, `mol.to_cdpkit()`, the `mol.iupac` property, and `rxn.to_rdkit()`, which is the only
side of a reaction any of the five has a form for.  A `cdef class` cannot be extended from outside, so
the methods reach the core through its registration hooks: the bodies registered through
`_set_interop_fns` at the foot of this file ARE the dispatchers below, and there is one implementation
per tool, reached two ways.  The import direction has no method: there is no `self` to hang "read this
foreign object" on.

THE DIRECTION ALSO DECIDES WHERE THE RECORDS GO.  An import returns a chython container, so what it
clamped or dropped is on that container's `.log` in stage `'interop'`, unconditionally; an export
returns a foreign object and keeps its `log=` list, there being no container to write to.  See
`_records.py`.
"""
from ..core import MoleculeContainer, QueryContainer, ReactionContainer


_CONTAINERS = (MoleculeContainer, QueryContainer, ReactionContainer)


def is_container(x, /) -> bool:
    """
    Is this one of chython's own containers, i.e. does a converter export it rather than read it?

    The dispatch predicate for every callable in this package, in one place so that adding a container
    type cannot fix four converters and miss the fifth.  `ReactionContainer` counts even where no
    converter exports a reaction yet, so the export side reports the gap instead of the import side
    complaining that a reaction is not an RDKit object.
    """
    return isinstance(x, _CONTAINERS)


def rdkit(x, /, **kwargs):
    """
    Convert to or from an RDKit molecule, or a reaction and a `ChemicalReaction`.

    The one converter here with a reaction form on both sides: a `ReactionContainer` exports as a
    `rdChemReactions.ChemicalReaction` and one imports back, carrying the atom-atom mapping each way.

    :param x: a chython container to export, or an RDKit `Mol`/`RWMol`/`ChemicalReaction` to import.
    :param kwargs: forwarded to the chosen direction.  The two take different keywords -- the exporter
        takes `keep_mapping`, `keep_numbers`, `keep_hydrogens`, `keep_coordinates` and `absolute`, the
        importer none -- so a keyword offered to the wrong direction raises `TypeError` naming the
        direction it reached rather than being silently ignored.
    """
    if is_container(x):
        from ._rdkit import to_rdkit

        return to_rdkit(x, **kwargs)
    from ._rdkit import from_rdkit

    return from_rdkit(x, **kwargs)


def indigo(x, /, **kwargs):
    """
    Convert to or from an Indigo molecule.

    :param x: a chython container to export, or an Indigo object to import.
    :param kwargs: forwarded to the chosen direction; see `rdkit` on keywords for the other direction.
    """
    if is_container(x):
        from ._indigo import to_indigo

        return to_indigo(x, **kwargs)
    from ._indigo import from_indigo

    return from_indigo(x, **kwargs)


def openbabel(x, /, **kwargs):
    """
    Convert to or from an OpenBabel `OBMol`.

    :param x: a chython container to export, or an `OBMol` to import.
    :param kwargs: forwarded to the chosen direction; see `rdkit` on keywords for the other direction.
    """
    if is_container(x):
        from ._openbabel import to_openbabel

        return to_openbabel(x, **kwargs)
    from ._openbabel import from_openbabel

    return from_openbabel(x, **kwargs)


def cdk(x, /, **kwargs):
    """
    Convert to or from a CDK `IAtomContainer`.  Starts a JVM through JPype on first use.

    :param x: a chython container to export, or an `IAtomContainer` to import.
    :param kwargs: forwarded to the chosen direction; see `rdkit` on keywords for the other direction.
    """
    if is_container(x):
        from ._cdk import to_cdk

        return to_cdk(x, **kwargs)
    from ._cdk import from_cdk

    return from_cdk(x, **kwargs)


def cdpkit(x, /, **kwargs):
    """
    Convert a chython container to a CDPKit molecule.

    Export only.  The import direction raises `DirectionNotImplemented` rather than `TypeError`, so a
    caller probing for capability can tell "that half is not built" from "wrong argument type".

    :param x: a chython container to export.
    :param kwargs: forwarded to the exporter.
    """
    if is_container(x):
        from ._cdpkit import to_cdpkit

        return to_cdpkit(x, **kwargs)
    from ._cdpkit import from_cdpkit

    return from_cdpkit(x, **kwargs)


def iupac(x, /, **kwargs):
    """
    Convert to or from an IUPAC name.

    Reading a name uses OPSIN (Java, through JPype); writing one uses openclatura, which needs
    Python >= 3.11 and names an RDKit molecule, so the export direction goes through `rdkit`.

    :param x: a chython container to name, or a `str` name to parse.
    :param kwargs: forwarded to the chosen direction; see `rdkit` on keywords for the other direction.
    """
    if is_container(x):
        from ._iupac import to_iupac

        return to_iupac(x, **kwargs)
    from ._iupac import from_iupac

    return from_iupac(x, **kwargs)


# At the bottom, after `is_container`, which `_pandas` reaches back for (inside a function, so this is
# reading order and not a cycle).
from ._pandas import patch_pandas
from ..core._core import _set_interop_fns


# The container methods, whose bodies are the six dispatchers above -- so a method is the export half of
# the published callable and cannot answer differently from it.  Registered here rather than compiled
# into the core for the reason every hook in `_molecule_container.pxi` gives: the direction is
# `core <- ... <- interop`, and a converter loads a toolkit the core may not name.
_set_interop_fns(rdkit=rdkit, indigo=indigo, openbabel=openbabel, cdk=cdk, cdpkit=cdpkit, iupac=iupac)


__all__ = ['rdkit', 'indigo', 'openbabel', 'cdk', 'cdpkit', 'iupac', 'is_container', 'patch_pandas']
