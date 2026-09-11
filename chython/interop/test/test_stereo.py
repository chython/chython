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
The shared stereo translation the converters read through, `interop._stereo`.

Tests marked LATENT pin the spec's reading of an invariant that today's core happens to satisfy
anyway; the rest pin live behaviour.
"""
from chython.core import SU_ALLENE, SU_ATROPISOMER, SU_CIS_TRANS, SU_HELICAL, SU_TETRA, read_smiles
from chython.interop import _stereo


def _unit(mol, kind):
    for u in mol.stereo_units():
        if u['kind'] == kind:
            return u
    return None


def test_kind_names_cover_every_kind_the_core_can_emit():
    """Every kind the core can emit is named, so no unit can fall off a converter's if/elif chain."""
    for kind in (SU_TETRA, SU_CIS_TRANS, SU_ALLENE, SU_ATROPISOMER, SU_HELICAL):
        assert kind in _stereo.KIND_NAMES
        assert _stereo.KIND_NAMES[kind]


def test_kind_names_are_read_from_the_core():
    """Values come from `chython.core`, not from a hand-copied table that can drift out of step."""
    assert set(_stereo.KIND_NAMES) == {SU_TETRA, SU_CIS_TRANS, SU_ALLENE, SU_ATROPISOMER, SU_HELICAL}


def test_first_of_each_pair_falls_back_past_an_unnamed_slot():
    """A half whose first slot is unnamed still names an atom in its second.

    LATENT: `_stereo.pxi` permits a bond kind's unnamed slot to sit inside a half, though today's
    perception always names the first slot.
    """
    assert _stereo.first_of_each_pair((None, 7, None, 9)) == (7, 9)
    assert _stereo.first_of_each_pair((7, None, 9, None)) == (7, 9)
    assert _stereo.first_of_each_pair((None, None, None, None)) == (None, None)


def test_cis_trans_partner_is_found_through_the_graph():
    """The far end of the bond is the neighbour the far half hangs off."""
    mol = read_smiles('F/C=C/F')
    u = _unit(mol, SU_CIS_TRANS)
    assert u is not None
    partner = _stereo.cis_trans_partner(mol, u['anchor'], u['refs'])
    assert partner in mol.neighbors_of(u['anchor'])
    far = _stereo.first_of_each_pair(u['refs'])[1]
    assert far in mol.neighbors_of(partner)


def test_cis_trans_partner_ignores_a_second_double_bond_on_the_anchor():
    """A cumulated system must not resolve to whichever double bond comes first in adjacency order.

    LATENT: no molecule is known where "first neighbour with order 2" disagrees, but the partner is
    defined by the refs, and an order-based rule would land the descriptor on the wrong bond silently.
    """
    mol = read_smiles('F/C=C/C=C/F')
    for u in mol.stereo_units():
        if u['kind'] != SU_CIS_TRANS:
            continue
        anchor, refs = u['anchor'], u['refs']
        partner = _stereo.cis_trans_partner(mol, anchor, refs)
        assert partner is not None
        far = _stereo.first_of_each_pair(refs)[1]
        # the partner is on the same bond as the far reference, which is the whole definition
        assert far in mol.neighbors_of(partner)
        naive = next((nb for nb in mol.neighbors_of(anchor) if mol.order_of(anchor, nb) == 2), None)
        assert partner == naive or far not in mol.neighbors_of(naive)


def test_cis_trans_frame_returns_a_translatable_order():
    """The framed order is one `translate_stereo` accepts, and the answer is 1 or 2."""
    mol = read_smiles('F/C=C/F')
    u = _unit(mol, SU_CIS_TRANS)
    frame = _stereo.cis_trans_frame(mol, u['anchor'], u['refs'])
    assert frame is not None
    assert mol.translate_stereo(u['anchor'], frame.order) in (1, 2)


def test_cis_trans_frame_agrees_with_the_stored_configuration():
    """E and Z of the same skeleton frame to opposite parities -- the calibration, not just the shape."""
    e = read_smiles('F/C=C/F')
    z = read_smiles('F/C=C\\F')
    parities = []
    for mol in (e, z):
        u = _unit(mol, SU_CIS_TRANS)
        frame = _stereo.cis_trans_frame(mol, u['anchor'], u['refs'])
        parities.append(mol.translate_stereo(u['anchor'], frame.order))
    assert parities[0] != parities[1]


def test_set_parity_by_probe_round_trips():
    """Writing a foreign parity by probing reads back as the parity that was asked for.

    Probing and not an arithmetic inverse: `set_parity` writes in the unit's refs order and
    `translate_stereo` reads in the caller's, and a sign error there mirrors a molecule silently.
    """
    mol = read_smiles('N[C@@H](C)C(=O)O')
    u = _unit(mol, SU_TETRA)
    anchor = u['anchor']
    order = tuple(u['refs'])
    for want in (1, 2):
        assert _stereo.set_parity_by_probe(mol, anchor, order, want)
        assert mol.translate_stereo(anchor, order) == want


def test_set_parity_by_probe_clears_and_reports_on_a_rejected_order():
    """A frame the core will not accept leaves the centre unset and says so, rather than half-set."""
    mol = read_smiles('N[C@@H](C)C(=O)O')
    u = _unit(mol, SU_TETRA)
    anchor = u['anchor']
    assert not _stereo.set_parity_by_probe(mol, anchor, (1, 2, 3, 999), 1)
    assert next(x['parity'] for x in mol.stereo_units() if x['anchor'] == anchor) == 0


def test_set_parity_by_probe_does_not_touch_an_unrelated_centre():
    """The probe is local: a failure at one anchor leaves every other configuration alone."""
    mol = read_smiles('N[C@@H](C)[C@H](O)C(=O)O')
    units = [u for u in mol.stereo_units() if u['kind'] == SU_TETRA and u['parity']]
    assert len(units) == 2
    before = {u['anchor']: u['parity'] for u in units}
    _stereo.set_parity_by_probe(mol, units[0]['anchor'], (1, 2, 3, 999), 1)
    after = {u['anchor']: u['parity'] for u in mol.stereo_units() if u['kind'] == SU_TETRA}
    assert after[units[1]['anchor']] == before[units[1]['anchor']]
