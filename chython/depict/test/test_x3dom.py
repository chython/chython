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
"""`mol.depict3d()` draws a stored conformer, and `mol.view3d()` puts that drawing in a notebook.

A conformer is the subject: the drawing reads the model the caller asks for and stores nothing, so a
molecule with no geometry is refused rather than drawn against an invented one -- the 2D side's
temporary-layout fallback has no counterpart here, because there is no such thing as a plausible
guess at a conformer.
"""
from pytest import raises

from chython import smiles
from chython.depict._config import R_COLOUR, cpk


def _placed(smi, *xyz):
    """The molecule with one model, its atoms placed in `atom_numbers` order."""
    mol = smiles(smi)
    for n, (x, y, z) in zip(mol.atom_numbers, xyz):
        mol.set_xyz(n, x, y, z)
    return mol


def _ethanol():
    return _placed('CCO', (0., 0., 0.), (1.5, 0., 0.), (2.2, 1.2, 0.))


# --- what it draws -------------------------------------------------------------------------- #

def test_one_sphere_per_atom():
    xml = _ethanol().depict3d()
    assert xml.count('<sphere') == 3


def test_the_document_is_one_x3d_scene():
    xml = _ethanol().depict3d()
    assert xml.startswith('<x3d ')
    assert xml.count('<scene>') == 1
    assert xml.rstrip().endswith('</x3d>')


def test_the_sphere_radius_is_the_atomic_radius():
    """The module's `atom_radius` is negative, which means "a multiplier, not a fixed size" -- so a
    carbon and an oxygen get different spheres, and that is the whole reason this file waited on
    `Atom.atomic_radius`."""
    xml = _ethanol().depict3d()
    assert "radius='0.13'" in xml       # carbon, 0.67 * .2
    assert "radius='0.10'" in xml       # oxygen, 0.48 * .2


def test_an_atom_is_drawn_in_its_cpk_colour():
    xml = _ethanol().depict3d()
    assert cpk[5] in xml                # carbon
    assert cpk[7] in xml                # oxygen


def test_each_bond_becomes_one_cylinder():
    xml = _ethanol().depict3d()
    assert xml.count('<cylinder') == 2


def test_a_double_bond_becomes_two_cylinders():
    xml = _placed('C=O', (0., 0., 0.), (1.2, 0., 0.)).depict3d()
    assert xml.count('<cylinder') == 2


def test_a_triple_bond_becomes_three_cylinders():
    xml = _placed('C#N', (0., 0., 0.), (1.2, 0., 0.)).depict3d()
    assert xml.count('<cylinder') == 3


def test_an_aromatic_ring_is_drawn_with_its_inner_dashes():
    """Benzene's six aromatic bonds are six cylinders plus a dashed inner circle, so the count is well
    above six; a kekulized benzene has no order-4 bond and draws no dashes at all."""
    ring = [(1.4, 0., 0.), (.7, 1.21, 0.), (-.7, 1.21, 0.), (-1.4, 0., 0.), (-.7, -1.21, 0.),
            (.7, -1.21, 0.)]
    aromatic = _placed('c1ccccc1', *ring).depict3d()
    kekule = _placed('C1=CC=CC=C1', *ring).depict3d()
    assert aromatic.count('<cylinder') > kekule.count('<cylinder')


def test_the_r_marker_is_labelled_and_not_a_sphere():
    """An R has no radius, so a sphere would be a point.  It is drawn as its own text in its own
    colour, which is what the 2D side does with it."""
    xml = _placed('[R]C', (0., 0., 0.), (1.5, 0., 0.)).depict3d()
    assert xml.count('<sphere') == 1
    assert "string='R'" in xml
    assert R_COLOUR in xml


# --- which model, and the refusals ---------------------------------------------------------- #

def test_the_index_selects_the_model():
    mol = _ethanol()
    second = mol.add_conformer()
    for n in mol.atom_numbers:
        mol.set_xyz(n, 5., 5., 5., second)
    assert mol.depict3d(0) != mol.depict3d(second)


def test_a_molecule_with_no_conformer_is_refused():
    with raises(ValueError, match='no conformer'):
        smiles('CCO').depict3d()


def test_a_model_that_does_not_exist_is_an_index_error():
    with raises(IndexError):
        _ethanol().depict3d(3)


def test_drawing_stores_nothing():
    mol = _ethanol()
    before = bytes(mol)
    mol.depict3d()
    assert bytes(mol) == before


def test_the_scene_is_centred_on_the_model():
    """Translating every atom by the same vector draws the same picture: the model is centred on its
    own centroid before rendering, so a conformer read out of a crystal file is not off-screen."""
    here = _ethanol().depict3d()
    there = _placed('CCO', (10., 10., 10.), (11.5, 10., 10.), (12.2, 11.2, 10.)).depict3d()
    assert here == there


# --- the notebook wrapper ------------------------------------------------------------------- #

def test_view3d_wraps_the_drawing_at_the_requested_size():
    widget = _ethanol().view3d(width='300px', height='200px')
    html = widget._repr_html_()
    assert 'width: 300px' in html
    assert 'height: 200px' in html
    assert '<x3d ' in html


def test_view3d_takes_the_same_model_index():
    mol = _ethanol()
    mol.add_conformer()
    assert mol.view3d(1)._repr_html_() != mol.view3d(0)._repr_html_()


def test_the_widget_renders_the_same_html_both_ways():
    """`__html__` is what a template engine calls and `_repr_html_` is what Jupyter calls."""
    widget = _ethanol().view3d()
    assert widget.__html__() == widget._repr_html_()


# --- the registration ----------------------------------------------------------------------- #

def test_the_core_owns_the_two_names_and_depict_supplies_the_bodies():
    """Same injection every other depiction method uses: the container is a `cdef class`, so a mixin
    cannot reach it and `chython.depict` registers the bodies at import."""
    from chython.core._core import MoleculeContainer

    assert MoleculeContainer.depict3d.__doc__
    assert MoleculeContainer.view3d.__doc__


def test_the_pair_is_one_registration():
    from chython.core._core import _set_depict_fns

    with raises(ValueError, match='depict3d'):
        _set_depict_fns(view3d=lambda *a, **kw: None)
