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
from chython.core import MoleculeContainer


def test_unstated_hydrogens_are_unknown_and_NOT_zero():
    """A bare `add_atom` states no hydrogen count, so the count reads as unknown.

    A ZERO HERE WOULD BE THE ARENA'S MEMSET SHOWING THROUGH, and a statement the builder has no
    standing to make: it tells `kekule()` that a two-coordinate aromatic nitrogen nobody has counted
    carries no hydrogen, which is pyridine, and pyrrole then has no kekule form at all.  The core
    derives nothing here either way -- that is what the rest of this file is about; what "nobody said"
    answers is unknown, not zero.
    """
    m = MoleculeContainer()
    c = m.add_atom(6)
    assert m.implicit_h_of(c) is None
    assert m.explicit_h_of(c) == 0, 'explicit hydrogens are counted from the graph, not stated'
    assert m.total_h_of(c) is None, 'a total over an unknown addend is unknown, not the other addend'
    assert m.hybridization_of(c) == 1


def test_explicit_hydrogens_are_counted_from_the_graph():
    m = MoleculeContainer()
    c = m.add_atom(6)
    h1, h2 = m.add_atom(1), m.add_atom(1)
    m.add_bond(c, h1, 1)
    m.add_bond(c, h2, 1)
    assert m.explicit_h_of(c) == 2
    assert m.implicit_h_of(c) is None    # nothing was stated
    # The two bonded hydrogens are COUNTED and the implicit count is UNKNOWN, and the total may not
    # quietly report the half it happens to have.  2 here would read as "this carbon has two
    # hydrogens" when the truthful answer is "it has two plus an unknown number".
    assert m.total_h_of(c) is None
    assert m.heteroatoms_of(c) == 0      # hydrogen is not a heteroatom
    assert m.degree_of(c) == 2


def test_stated_and_explicit_hydrogens_add_up():
    m = MoleculeContainer()
    c = m.add_atom(6, implicit_h=3)
    h = m.add_atom(1)
    m.add_bond(c, h, 1)
    assert m.implicit_h_of(c) == 3
    assert m.explicit_h_of(c) == 1
    assert m.total_h_of(c) == 4


def test_stated_hydrogen_count_survives_the_pass_untouched():
    m = MoleculeContainer()
    c1 = m.add_atom(6, implicit_h=0)
    c2 = m.add_atom(6, implicit_h=9)     # chemically absurd, deliberately
    assert m.implicit_h_of(c1) == 0
    assert m.implicit_h_of(c2) == 9


def test_heteroatoms_exclude_carbon_and_hydrogen():
    m = MoleculeContainer()
    c1, c2, o = m.add_atom(6), m.add_atom(6), m.add_atom(8)
    h = m.add_atom(1)
    m.add_bond(c1, c2, 1)
    m.add_bond(c2, o, 1)
    m.add_bond(o, h, 1)
    assert m.heteroatoms_of(c1) == 0
    assert m.heteroatoms_of(c2) == 1     # the O
    assert m.heteroatoms_of(o) == 0      # a C and an H, neither counts
    assert m.explicit_h_of(o) == 1


def test_hybridization_cases():
    m = MoleculeContainer()
    sp3 = m.add_atom(6)
    sp2 = m.add_atom(6)
    o1 = m.add_atom(8)
    sp = m.add_atom(6)
    n1 = m.add_atom(7)
    cumul = m.add_atom(6)
    o2, o3 = m.add_atom(8), m.add_atom(8)
    weird = m.add_atom(6)
    n2 = m.add_atom(7)
    o4 = m.add_atom(8)
    m.add_bond(sp3, sp2, 1)
    m.add_bond(sp2, o1, 2)
    m.add_bond(sp, n1, 3)
    m.add_bond(cumul, o2, 2)
    m.add_bond(cumul, o3, 2)
    m.add_bond(weird, n2, 3)
    m.add_bond(weird, o4, 2)
    assert m.hybridization_of(sp3) == 1
    assert m.hybridization_of(sp2) == 2
    assert m.hybridization_of(sp) == 3
    assert m.hybridization_of(cumul) == 5    # CO2, cumulated double bonds
    assert m.hybridization_of(weird) == 6    # unmatched combination


def test_hybridization_never_reports_aromatic():
    m = MoleculeContainer()
    ring = [m.add_atom(6) for _ in range(6)]
    for k in range(6):
        m.add_bond(ring[k], ring[(k + 1) % 6], 2 if k % 2 == 0 else 1)
    # a Kekule benzene: every carbon has exactly one double bond
    assert [m.hybridization_of(x) for x in ring] == [2] * 6


def test_dative_bond_contributes_no_pi_bond():
    m = MoleculeContainer()
    n = m.add_atom(7)
    bo = m.add_atom(5)
    m.add_bond(n, bo, 8)
    assert m.hybridization_of(n) == 1
    assert m.hybridization_of(bo) == 1
    assert m.heteroatoms_of(n) == 1      # boron is a heteroatom
    assert m.heteroatoms_of(bo) == 1     # so is nitrogen


def test_degree_and_heteroatoms_on_a_crowded_atom():
    m = MoleculeContainer()
    c = m.add_atom(6)
    ns = [m.add_atom(7) for _ in range(5)]
    for x in ns:
        m.add_bond(c, x, 1)
    assert m.heteroatoms_of(c) == 5
    assert m.degree_of(c) == 5
    assert m.hybridization_of(c) == 1


def test_isolated_atom_of_any_element_derives_cleanly():
    m = MoleculeContainer()
    u = m.add_atom(92)
    fe = m.add_atom(26)
    for x in (u, fe):
        assert m.heteroatoms_of(x) == 0
        assert m.explicit_h_of(x) == 0
        assert m.hybridization_of(x) == 1


def test_explicit_h_saturates_rather_than_wrapping():
    m = MoleculeContainer()
    c = m.add_atom(6)
    with m.edit():
        for _ in range(16):
            m.add_bond(c, m.add_atom(1), 1)
    # 16 must clamp to 15, not wrap to 0 — the nibble is four bits wide
    assert m.explicit_h_of(c) == 15


def test_hydrogen_atom_derives_its_own_scalars():
    m = MoleculeContainer()
    c = m.add_atom(6)
    h = m.add_atom(1)
    m.add_bond(c, h, 1)
    assert m.heteroatoms_of(h) == 0       # its only neighbour is carbon
    assert m.explicit_h_of(h) == 0        # and carbon is not a hydrogen


def test_derive_handles_an_empty_molecule():
    m = MoleculeContainer()
    c = m.add_atom(6)
    m.delete_atom(c)                     # the fold runs with zero survivors
    assert m.atom_count == 0
    assert m.atom_numbers == []
