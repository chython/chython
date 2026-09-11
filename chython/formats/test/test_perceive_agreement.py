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
"""`perceive_bonds` against connectivity a file stated, on the tree's own 3D records.

`chemistry/test/test_perceive.py` probes the threshold with built geometries of small molecules and
`test_covalent_radii_tsv.py` with tabulated bond lengths and contacts.  What neither can show is a whole
structure: 60-odd atoms, fused rings, two coordinated metals, and every bond written down by whoever
stored the record.  The comparison lives HERE rather than beside the pass because reading an SDF needs
`formats`, which `chemistry` may not see.

`test/cycle.sdf` is the ring-perception corpus and read-only here.
"""
from pathlib import Path

from ..ctfile import SDFRead


_DATA = Path(__file__).resolve().parent.parent.parent.parent / 'test'


def _records():
    """Every record of `cycle.sdf` whose first model states a z coordinate."""
    with SDFRead(str(_DATA / 'cycle.sdf')) as f:
        return [m for m in f
                if m.has_3d and any(m.conformer(0).xyz_of(n)[2] for n in m.atom_numbers)]


def _geometries():
    """Those of them that are one connected molecule -- a geometry rather than a placed lattice.

    One record of the file holds five free atoms -- three Ru, one Ni, one hydride -- laid out a flat
    1 A apart along the x axis.  A distance rule reads a 1 A pair as bonded, correctly, so that record
    says nothing about agreement and is not in this corpus.  It is still in the one below it, which only
    asks that nothing STATED is missed.
    """
    return [m for m in _records() if len(m.connected_components) == 1]


def _perceived(molecule):
    """The bond set `perceive_bonds` finds from `molecule`'s first model alone, keyed on its own ids.

    A bondless copy, so the answer is the geometry's and not the record's: the pass leaves a bond that
    is already there alone, and comparing against a molecule that has them all would assert nothing.
    """
    from ...chemistry import perceive_bonds
    from ...core import MoleculeContainer

    conformer = molecule.conformer(0)
    bare = MoleculeContainer()
    ids = [(n, bare.add_atom(molecule.atom(n).atomic_symbol, implicit_h=0))
           for n in molecule.atom_numbers]
    for n, m in ids:
        bare.set_xyz(m, *conformer.xyz_of(n))
    perceive_bonds(bare)

    back = {m: n for n, m in ids}
    return {frozenset((back[b.n], back[b.m])) for b in bare.bonds()}


def _stated(molecule):
    return {frozenset((b.n, b.m)) for b in molecule.bonds()}


def test_the_corpus_holds_what_this_file_claims_it_does():
    """Three 3D records, two of them one molecule each -- asserted, so a shrunk corpus is not silence."""
    assert len(_records()) == 3
    assert [len(m) for m in _geometries()] == [62, 57]


def test_no_bond_a_record_states_is_missed():
    """Over every 3D record, the lattice included: the threshold reaches every bond that was written.

    A missed bond is the failure with no recovery downstream -- ``saturate()`` never adds one -- so this
    is the half of agreement that matters most.
    """
    missed = {}
    for molecule in _records():
        gap = _stated(molecule) - _perceived(molecule)
        if gap:
            missed[str(molecule)] = sorted(sorted(pair) for pair in gap)
    assert not missed, f'bond(s) stated by a record and not perceived from its geometry: {missed}'


def test_a_connected_geometry_is_reproduced_bond_for_bond():
    """Exact agreement, not a tolerance: 81 and 60 bonds, nothing missed and nothing invented.

    Both records are metal-organic -- porphyrin-like macrocycles around two Ru centres -- so this covers
    long coordination bonds, which a covalent-radius threshold reaches only because the metal radii are
    large, and fused aromatic rings, whose 1,3 contacts are what a loose threshold bonds.

    [mutant: `radius_multiplier` at 1.3 -- this test fails with two invented bonds, the diagonals of the
    62-atom record's four-membered carbocycle, whose 1.39 A sides put them 1.966 A apart.]
    """
    disagreement = {}
    for molecule in _geometries():
        stated, perceived = _stated(molecule), _perceived(molecule)
        if stated != perceived:
            disagreement[str(molecule)] = {'missed': sorted(sorted(p) for p in stated - perceived),
                                           'invented': sorted(sorted(p) for p in perceived - stated)}
    assert not disagreement, f'geometry and record disagree: {disagreement}'
