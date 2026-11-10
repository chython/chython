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
"""ChemDraw CDX and CDXML: one object tree in two spellings, read the way a molfile is.

:mod:`._binary` and :mod:`._cdxml` tokenize into :class:`~._tree.CdxObject`, every value in its CDXML
spelling; :mod:`._build` turns the tree into the :class:`~chython.formats.ctfile._ctab.Ctab` every MDL
path goes through, so hydrogens, wedges, enhanced stereo and S-groups are the V3000 code.  The
vocabulary is generated from Revvity's open-source ChemDraw library (BSD-3-Clause).
"""

from ._build import has_reaction, read_document, read_reaction
from ._binary import parse_cdx
from ._cdxml import parse_cdxml
from ._errors import CdxError, MalformedCdx, UnsupportedCdx
from ._write import write_cdx, write_cdxml
from ...core import MoleculeContainer
from ...core.reaction import ReactionContainer


__all__ = ['cdx', 'cdxml', 'CdxError', 'MalformedCdx', 'UnsupportedCdx']


def cdx(data, *, step=None, ignore_stereo=False, log=None):
    """A molecule or reaction from CDX bytes, or CDX bytes from either.

    :param step: import only; which reaction step to read, in document order across schemes, negative
        from the end.  ``None`` (default) reads the first step of a document that has one and a molecule
        otherwise; an index the document does not hold raises ``IndexError``.
    :param ignore_stereo: import only; skip the stereo step, read the constitution either way.
    :param log: a list to append damage reports to.

    A document holding a reaction step reads as a `ReactionContainer`: the step's reactants, products
    and the structures above and below its arrow as agents, its atom map as map numbers.  Any other
    document reads as one molecule, every structure on every page as components.
    """
    log = [] if log is None else log
    if isinstance(data, (MoleculeContainer, ReactionContainer)):
        return write_cdx(data, log)
    root = parse_cdx(data, log)
    return _build(root, log, ignore_stereo, step)


def cdxml(data, *, step=None, ignore_stereo=False, log=None):
    """A molecule or reaction from CDXML text, or CDXML text from either.  Arguments are :func:`cdx`'s."""
    log = [] if log is None else log
    if isinstance(data, (MoleculeContainer, ReactionContainer)):
        return write_cdxml(data, log)
    root = parse_cdxml(data, log)
    return _build(root, log, ignore_stereo, step)


def _build(root, log, ignore_stereo, step):
    if step is not None and not has_reaction(root):
        raise IndexError(f'step {step} requested; the document holds no reaction step')
    if not has_reaction(root):
        ctab = read_document(root, log)
        molecule, _, built = ctab.build(ignore_stereo=ignore_stereo)
        log[:] = built
        return molecule
    document = list(log)        # the tokenizer's lines: the document's, so the reaction's
    own = []
    reactants, products, agents, meta = read_reaction(root, own, 0 if step is None else step)
    sides = []
    for side in (reactants, products, agents):
        molecules = []
        for ctab in side:
            molecule, _, built = ctab.build(ignore_stereo=ignore_stereo)
            log.extend(built)
            molecules.append(molecule)
        sides.append(molecules)
    reaction = ReactionContainer(*sides, meta=meta)
    if document or own:
        reaction.log.absorb('read', document + own)
    log.extend(own)
    return reaction
