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
"""Public peptides and the helpers the polymer tests share."""

from .._assemble import assemble
from .._helm import parse_helm
from .._library import get_monomer_library


#: Oxytocin, CYIQNCPLG-NH2 with the Cys1-Cys6 disulfide; the isomeric SMILES is PubChem CID 439302's.
OXYTOCIN_HELM = 'PEPTIDE1{C.Y.I.Q.N.C.P.L.G.[am]}$PEPTIDE1,PEPTIDE1,1:R3-6:R3$$$'
OXYTOCIN_SMILES = ('CC[C@H](C)[C@H]1C(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@@H](CSSC[C@@H](C(=O)N[C@H](C(=O)N1)'
                   'CC2=CC=C(C=C2)O)N)C(=O)N3CCC[C@H]3C(=O)N[C@@H](CC(C)C)C(=O)NCC(=O)N)CC(=O)N)CCC(=O)N')


def rules(log):
    return {r.rule for r in log}


def same(a, b):
    """One structure, whichever of aromatic or Kekule each is written in."""
    a, b = a.copy(), b.copy()
    a.kekule()
    b.kekule()
    return a == b


def assemble_text(text, log=None):
    """The molecule of a HELM string, with the default library."""
    return assemble(parse_helm(text), get_monomer_library(), [] if log is None else log)
