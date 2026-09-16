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
"""WHICH ATOM A STEREO UNIT IS KEYED AT, per kind, as a ratchet.

A unit's stereo facts are keyed at its ANCHOR slot -- `SEG_PARITY` stores the parity there -- so a
reader that reports a unit as an OWNER PAIR derives that pair from the anchor, never from arithmetic
over the pair.  These are the rules that derivation inverts:

| Kind | Owners | Anchor |
| --- | --- | --- |
| 0 `SU_TETRA` | the centre | the centre |
| 1 `SU_CIS_TRANS` | the two chain terminals | the lower-index terminal |
| 2 `SU_ALLENE` | the two chain terminals | the chain's midpoint ATOM |
| 3 `SU_ATROPISOMER` | the two pivots | a pivot, and NOT always the lower-index one |

The atropisomer row is the reason a group's write path resolves an owner pair through the unit table
rather than through `min()`: `_perceive_stereo_units`'s pass 3 relocates an axis to the higher pivot when
the lower one is already anchored, which `test_an_atropisomer_may_anchor_at_either_pivot` reaches.

Stable ids equal slot order for a molecule built by the SMILES reader, which is what lets "the
lower-index terminal" be asserted as `min(pair)` over stable ids here.

`stereo_units()` reports CANDIDATES, so a methyl carbon is a kind-0 record in it; every assertion
about a kind either filters on `stereogenic` or reads `chiral_atoms()` / `chiral_bonds()`, which
filter already.
"""
import pytest

from chython.core._core import read_smiles as smiles


def _chain_atoms(mol, terminal):
    """The double-bond chain from `terminal`, as stable ids in walk order."""
    chain = [terminal]
    prev = None
    cur = terminal
    while True:
        nxt = next((n for n in mol.neighbors_of(cur) if mol.order_of(cur, n) == 2 and n != prev), None)
        if nxt is None:
            return chain
        prev, cur = cur, nxt
        chain.append(cur)


def _units(mol, kind):
    return [u for u in mol.stereo_units() if u['kind'] == kind]


def test_a_tetrahedral_centre_anchors_at_itself():
    mol = smiles('C[C@H](N)C(=O)O')                     # alanine
    assert list(mol.chiral_atoms()) == [2]
    assert [u['anchor'] for u in _units(mol, 0) if u['stereogenic']] == [2]


@pytest.mark.parametrize('text, axes', [
    ('F/C=C/F', 1),                                     # 1,2-difluoroethene
    ('OC(=O)/C=C\\C(=O)O', 1),                          # maleic acid
    ('C/C=C/C=C/C', 2),                                 # (2E,4E)-hexa-2,4-diene: two axes
])
def test_a_cis_trans_axis_anchors_at_its_lower_owner(text, axes):
    mol = smiles(text)
    bonds = mol.chiral_bonds()
    checked = [p for p, u in bonds.items() if u['kind'] == 1]
    assert len(checked) == axes, bonds
    for pair in checked:
        assert bonds[pair]['anchor'] == min(pair)


@pytest.mark.parametrize('text', [
    'CC=C=CC',                                          # penta-2,3-diene: 3 chain atoms
    'CC=C=C=C=CC',                                      # hepta-2,3,4,5-tetraene: 5 chain atoms
])
def test_an_allene_anchors_at_its_chain_midpoint_atom(text):
    mol = smiles(text)
    units = _units(mol, 2)
    assert len(units) == 1
    terminals = [n for n in mol.atom_numbers
                 if sum(1 for m in mol.neighbors_of(n) if mol.order_of(n, m) == 2) == 1]
    chain = _chain_atoms(mol, min(terminals))
    assert len(chain) & 1, 'an allene chain has an odd atom count'
    assert units[0]['anchor'] == chain[len(chain) // 2]
    # THE MIDPOINT IS NOT AN OWNER, which is why no arithmetic over the owner pair reaches the anchor
    assert units[0]['anchor'] not in (chain[0], chain[-1])
    # `chiral_atoms` keys on the midpoint; the axis is not in `chiral_bonds` at all
    assert units[0]['anchor'] in mol.chiral_atoms()
    assert mol.chiral_bonds() == {}


@pytest.mark.parametrize('text', [
    'CC=C=C=CC',                                        # hexa-2,3,4-triene: 4 chain atoms
    'CC=C=C=C=C=CC',                                    # octa-2,3,4,5,6-pentaene: 6 chain atoms
])
def test_an_even_cumulene_anchors_at_its_lower_owner_and_owns_a_non_bond_pair(text):
    mol = smiles(text)
    units = _units(mol, 1)
    assert len(units) == 1
    bonds = mol.chiral_bonds()
    assert len(bonds) == 1
    (lo, hi), unit = next(iter(bonds.items()))
    assert unit['anchor'] == lo
    # THE DEFECT THE HARMONIZED KEY CLOSES: the owner pair is not a bond
    assert hi not in mol.neighbors_of(lo)


def test_an_atropisomer_anchors_at_a_pivot():
    mol = smiles('Oc1ccc2ccccc2c1-c1c(O)ccc2ccccc12')   # 1,1'-binaphthyl-2,2'-diol (BINOL)
    axes = [(p, u) for p, u in mol.chiral_bonds().items() if u['kind'] == 3]
    assert len(axes) == 1
    pair, unit = axes[0]
    assert unit['anchor'] in pair


def test_an_atropisomer_may_anchor_at_either_pivot():
    """RULING 1's WITNESS.  A pivot claimed by another unit sends the axis to the OTHER pivot, so the
    anchor is not a function of the owner pair and `min(owners)` is not the resolution.

    2,2'-dimethyl-1,1'-bi(cycloocta-1,3,5,7-tetraenyl), written so that the lower pivot opens its
    ring and is therefore the lower terminal of a chain bond: the cis/trans cut admits a ring of
    eight, so that unit claims the lower pivot and the axis lands on the higher one.  The `.` is a
    separator, not a disconnection — ring closure 2 spans it and is the axis itself.
    """
    mol = smiles('C12=CC=CC=CC=C1C.CC1=CC=CC=CC=C12')
    axes = [(p, u) for p, u in mol.chiral_bonds().items() if u['kind'] == 3]
    assert len(axes) == 1
    pair, unit = axes[0]
    assert pair == (1, 18)
    # the lower pivot is a cis/trans anchor, so the axis is keyed on the HIGHER one
    assert min(pair) in {u['anchor'] for u in _units(mol, 1)}
    assert unit['anchor'] == max(pair)


def test_every_unit_anchor_is_unique():
    """The invariant the storage layout depends on, over one molecule carrying three kinds at once."""
    mol = smiles('C[C@H](O)/C=C/C=C=CC1=CC=CC=C1C1=CC=CC=C1')
    anchors = [u['anchor'] for u in mol.stereo_units()]
    assert len(anchors) == len(set(anchors))
