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
"""The repair a SMIRKS product's unknown hydrogen count is documented to have.

An atom the patch wrote inside an aromatic ring may come out with `H_UNKNOWN`, never a guessed zero.
The repair is `kekule()`, which closes the count itself.  Only the pyrrole-versus-pyridine nitrogen
still needs it; an atom in a state no valence row describes stays unknown, which is the honest answer.
"""
from chython.chemistry import calc_implicit
from chython.core import read_smiles, kekule
from chython.core._core import read_smirks


def test_kekule_resolves_an_aromatic_centre_by_itself():
    """The documented repair, in one call.

    N-demethylating N-methylpyrrole leaves a two-coordinate aromatic nitrogen with nothing stating its
    hydrogens -- pyrrole versus pyridine, the one class no local look can settle, because both
    readings are valences and the ring chooses.  `kekule()` chooses and writes the count that follows,
    so a caller running the pipeline's stages by hand gets what `canonicalize()` gives.
    """
    t = read_smirks('[N;a;D3:1]-[C;D1;z1]>>[N;a;D2:1]')
    product = next(iter(t(read_smiles('Cn1cccc1')))).products[0]
    assert product.implicit_h_of(2) is None
    assert product.unknown_h_count == 1

    kekule(product)
    assert product.implicit_h_of(2) == 1, 'kekule() heals the counts its own orders make derivable'
    assert product.unknown_h_count == 0

    # `calc_implicit` finds nothing left: one shared derivation means the two passes cannot disagree,
    # so the order they run in does not matter.
    assert calc_implicit(product, 2) == 1
    assert product.implicit_h_of(2) == 1


def test_the_patcher_needs_no_repair_where_every_kekule_form_agrees():
    """An aromatic atom whose class every Kekule form agrees on never goes unknown.

    Compared against the same molecule read straight from SMILES, because a number instead of None is
    only an improvement if it is the right number.
    """
    t = read_smirks('[C;a:1][Br;D1]>>[C;a:1][O;D1:2]')
    product = next(iter(t(read_smiles('c1ccccc1Br')))).products[0]
    assert product.unknown_h_count == 0

    reference = read_smiles('c1ccccc1O')
    assert ([product.implicit_h_of(n) for n in product.atom_numbers]
            == [reference.implicit_h_of(n) for n in reference.atom_numbers])


def test_an_underivable_state_stays_unknown_after_the_repair():
    """A carbon with six single bonds matches no valence rule, so it stays `H_UNKNOWN` after both
    passes -- "unknown" being the only honest answer left to a pass that repairs in place.
    """
    t = read_smirks('[C;D0:1]>>[C:1](-[F;D1:2])(-[F;D1:3])(-[F;D1:4])(-[F;D1:5])(-[F;D1:6])-[F;D1:7]')
    product = next(iter(t(read_smiles('C')))).products[0]
    assert product.implicit_h_of(1) is None

    kekule(product)
    calc_implicit(product, 1)
    assert product.implicit_h_of(1) is None
