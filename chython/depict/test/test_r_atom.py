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
"""An R depicts as its own label, in its own colour."""
from chython import smiles
from chython.depict._config import R_COLOUR, cpk
from chython.depict.label import element_colour, is_labelled, labels, to_hex
from chython.depict.style import DepictStyle


def _laid_out(smi):
    """The molecule, its plane, and the id of its one R."""
    mol = smiles(smi)
    mol.clean2d()
    r = next(sid for sid in mol if mol.atom(sid).is_r)
    return mol, mol.coordinates(), r


def test_an_r_is_labelled():
    mol, plane, r = _laid_out('[R]c1ccccc1')
    assert is_labelled(mol.atom(r), DepictStyle())


def test_the_label_reads_r():
    mol, plane, r = _laid_out('[R]c1ccccc1')
    runs = labels(mol, plane, DepictStyle())[r].text.runs
    assert [run.text for run in runs] == ['R']


def test_an_indexed_r_reads_its_index_in_one_upright_full_size_run():
    # `R3` is a label, not a descriptor: the index is not a subscript.  `label.py`'s own rule for the
    # stereo group id, and `atomic_symbol` already answers the two characters together.
    mol, plane, r = _laid_out('[R3]c1ccccc1')
    style = DepictStyle()
    runs = labels(mol, plane, style)[r].text.runs
    assert [run.text for run in runs] == ['R3']
    assert runs[0].size == style.label.size
    assert runs[0].dy == 0


def test_an_r_draws_no_unknown_hydrogen_mark():
    # Cross-checks Task 7: an R's implicit count is 0, not unknown, so `_compose` writes no `?`.
    mol, plane, r = _laid_out('[R]c1ccccc1')
    style = DepictStyle().tuned(**{'atom.hydrogens': True, 'atom.unknown_h_marks': True})
    assert '?' not in ''.join(run.text for run in labels(mol, plane, style)[r].text.runs)


def test_the_colour_is_r_colour_and_not_the_end_of_the_palette():
    # Both sides go through `to_hex`: it lowercases, and `cpk` is spelt in upper case, so comparing a
    # returned colour against a raw table entry would pass whatever the palette holds.
    mol, plane, r = _laid_out('[R]C')
    assert element_colour(mol.atom(r), DepictStyle()) == to_hex(R_COLOUR)
    assert element_colour(mol.atom(r), DepictStyle()) != to_hex(cpk[-1])


def test_an_r_colour_is_not_an_element_colour():
    # `R_COLOUR` has to be distinguishable from every CPK entry, or the marker reads as an element.
    assert to_hex(R_COLOUR) not in [to_hex(c) for c in cpk]


def test_a_molecule_with_an_r_renders():
    mol, plane, r = _laid_out('[R]c1ccccc1')
    svg = mol.depict()
    assert R_COLOUR.lower() in svg.lower()
    assert '>R<' in svg
