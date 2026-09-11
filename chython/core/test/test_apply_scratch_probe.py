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
from chython.core import _core as _structure


def _req(count, esz):
    """Required bytes for a region: same (count if count else 1) guard that _apply uses."""
    return (count if count else 1) * esz


def _present_offsets(p):
    """Ordered list of (off, key) pairs for regions that are present (non-None)."""
    keys = ['work_off', 'live_off', 'newidx_off', 'edits_off', 'wedge_off',
            'wxy_off', 'wsg_off', 'wpar_off', 'cip_off']
    return [(p[k], k) for k in keys if p[k] is not None]


# ---------------------------------------------------------------------------
# geometry tests (alignment, ordering, presence/absence)
# ---------------------------------------------------------------------------

def test_probe_all_offsets_are_8_aligned():
    p = _structure._apply_scratch_probe(4, 3, 2, True, True)
    for key in ('work_off', 'live_off', 'newidx_off', 'edits_off', 'wedge_off',
                'wxy_off', 'wsg_off', 'cip_off'):
        assert p[key] % 8 == 0, f'{key} = {p[key]} is not 8-aligned'
    assert p['total'] % 8 == 0


def test_probe_offsets_are_monotonically_increasing():
    p = _structure._apply_scratch_probe(4, 3, 2, True, True)
    offs = _present_offsets(p)
    assert len(offs) == 8, f'expected 8 present regions, got {len(offs)}'
    for i in range(len(offs) - 1):
        assert offs[i][0] < offs[i + 1][0], (
            f'{offs[i][1]}={offs[i][0]} not < {offs[i+1][1]}={offs[i+1][0]}'
        )


def test_probe_wxy_absent_when_not_wanted():
    p = _structure._apply_scratch_probe(4, 3, 2, False, True)
    assert p['wxy_off'] is None
    assert p['wxy_esz'] is None
    # 7 present: work, live, newidx, edits, wedge, wsg, cip
    offs = _present_offsets(p)
    assert len(offs) == 7, f'expected 7 present regions, got {len(offs)}'
    assert p['wsg_off'] is not None


def test_probe_wsg_absent_when_not_wanted():
    p = _structure._apply_scratch_probe(4, 3, 2, True, False)
    assert p['wsg_off'] is None
    assert p['wsg_esz'] is None
    offs = _present_offsets(p)
    assert len(offs) == 7, f'expected 7 present regions, got {len(offs)}'
    assert p['wxy_off'] is not None


def test_probe_both_optional_segments_absent():
    p = _structure._apply_scratch_probe(4, 3, 2, False, False)
    assert p['wxy_off'] is None
    assert p['wsg_off'] is None
    assert p['wxy_esz'] is None
    assert p['wsg_esz'] is None
    # cip has no gate, so it is present even here: work, live, newidx, edits, wedge, cip
    offs = _present_offsets(p)
    assert len(offs) == 6, f'expected 6 present regions, got {len(offs)}'


def test_probe_zero_counts_keep_monotonic_offsets():
    """With all zeros, else-1 guards keep every offset distinct and monotonic."""
    p = _structure._apply_scratch_probe(0, 0, 0, False, False)
    assert p['work_off'] == 0
    offs = _present_offsets(p)
    assert len(offs) == 6
    for i in range(len(offs) - 1):
        assert offs[i][0] < offs[i + 1][0]


def test_probe_total_grows_with_counts():
    small = _structure._apply_scratch_probe(1, 1, 1, False, False)
    large = _structure._apply_scratch_probe(100, 100, 100, False, False)
    assert large['total'] > small['total']


# ---------------------------------------------------------------------------
# sufficiency tests — the test computes expected bytes from n/e/w and the
# element sizes returned by the probe; it does NOT call _scratch_sizes.
# Two independent derivations that agree is the invariant.
# ---------------------------------------------------------------------------

def test_probe_sufficiency_each_region_fits_its_data():
    """Each region is wide enough for the data _apply writes into it.

    n=4, e=3, w=2 ensures the three counts differ, so a region indexed by the
    wrong count would produce a wrong size and be caught here.
    """
    n, e, w = 4, 3, 2
    p = _structure._apply_scratch_probe(n, e, w, True, True)

    r_work   = _req(n, p['work_esz'])    # atom_t per atom
    r_live   = _req(n, p['live_esz'])    # uint8_t per atom
    r_newidx = _req(n, p['newidx_esz'])  # int32_t per atom
    r_edits  = _req(e, p['edits_esz'])   # edge_edit_t per bond
    r_wedge  = _req(w, p['wedge_esz'])   # wedge_edit_t per wedge slot
    r_wxy    = _req(n, p['wxy_esz'])     # xy_t per atom
    r_wsg    = _req(n, p['wsg_esz'])     # uint8_t per atom

    assert p['work_off']   + r_work   <= p['live_off'],   \
        f"work region too small: off={p['work_off']}, req={r_work}, next={p['live_off']}"
    assert p['live_off']   + r_live   <= p['newidx_off'], \
        f"live region too small: off={p['live_off']}, req={r_live}, next={p['newidx_off']}"
    assert p['newidx_off'] + r_newidx <= p['edits_off'],  \
        f"newidx region too small: off={p['newidx_off']}, req={r_newidx}, next={p['edits_off']}"
    assert p['edits_off']  + r_edits  <= p['wedge_off'],  \
        f"edits region too small: off={p['edits_off']}, req={r_edits}, next={p['wedge_off']}"
    assert p['wedge_off']  + r_wedge  <= p['wxy_off'],    \
        f"wedge region too small: off={p['wedge_off']}, req={r_wedge}, next={p['wxy_off']}"
    assert p['wxy_off']    + r_wxy    <= p['wsg_off'],    \
        f"wxy region too small: off={p['wxy_off']}, req={r_wxy}, next={p['wsg_off']}"
    assert p['wsg_off']    + r_wsg    <= p['cip_off'],    \
        f"wsg region too small: off={p['wsg_off']}, req={r_wsg}, next={p['cip_off']}"
    # The CIP region is sized by BONDS and is unconditional -- no `want_` gate, because a region that
    # is sometimes absent is how an index that is right for one version goes wrong for another.
    assert p['cip_off']    + _req(e, p['cip_esz']) <= p['total'], \
        f"cip region too small: off={p['cip_off']}, total={p['total']}"


def test_probe_sufficiency_without_optional_segments():
    """Sufficiency with no xy or sg — the last present region is then cip, which has no gate."""
    n, e, w = 4, 3, 2
    p = _structure._apply_scratch_probe(n, e, w, False, False)

    r_work   = _req(n, p['work_esz'])
    r_live   = _req(n, p['live_esz'])
    r_newidx = _req(n, p['newidx_esz'])
    r_edits  = _req(e, p['edits_esz'])
    r_wedge  = _req(w, p['wedge_esz'])

    assert p['work_off']   + r_work   <= p['live_off']
    assert p['live_off']   + r_live   <= p['newidx_off']
    assert p['newidx_off'] + r_newidx <= p['edits_off']
    assert p['edits_off']  + r_edits  <= p['wedge_off']
    assert p['wedge_off']  + r_wedge  <= p['cip_off']
    assert p['cip_off']    + _req(e, p['cip_esz']) <= p['total']


def test_probe_region_ends_are_within_total():
    """Every region end (off + required) is inside total, not just the start."""
    n, e, w = 10, 8, 6
    p = _structure._apply_scratch_probe(n, e, w, True, True)

    ends = [
        p['work_off']   + _req(n, p['work_esz']),
        p['live_off']   + _req(n, p['live_esz']),
        p['newidx_off'] + _req(n, p['newidx_esz']),
        p['edits_off']  + _req(e, p['edits_esz']),
        p['wedge_off']  + _req(w, p['wedge_esz']),
        p['wxy_off']    + _req(n, p['wxy_esz']),
        p['wsg_off']    + _req(n, p['wsg_esz']),
        p['cip_off']    + _req(e, p['cip_esz']),
    ]
    for end in ends:
        assert end <= p['total'], f'region end {end} exceeds total {p["total"]}'


def test_probe_zero_counts_sufficiency():
    """Else-1 guard must produce at least one element's worth of space per region."""
    p = _structure._apply_scratch_probe(0, 0, 0, True, True)
    # With count=0 the else-1 guard gives room for 1 element
    assert p['work_off']   + p['work_esz']   <= p['live_off']
    assert p['live_off']   + p['live_esz']   <= p['newidx_off']
    assert p['newidx_off'] + p['newidx_esz'] <= p['edits_off']
    assert p['edits_off']  + p['edits_esz']  <= p['wedge_off']
    assert p['wedge_off']  + p['wedge_esz']  <= p['wxy_off']
    assert p['wxy_off']    + p['wxy_esz']    <= p['wsg_off']
    assert p['wsg_off']    + p['wsg_esz']    <= p['cip_off']
    assert p['cip_off']    + p['cip_esz']    <= p['total']


# ---------------------------------------------------------------------------
# the parity region — gated like wxy and wsg, sized one byte per work slot
# ---------------------------------------------------------------------------

def test_probe_wpar_absent_when_not_wanted():
    p = _structure._apply_scratch_probe(4, 3, 2, True, True)
    assert p['wpar_off'] is None
    assert p['wpar_esz'] is None


def test_probe_wpar_sits_between_wsg_and_cip():
    """Its position in the carve is what the `_bp` chain in _apply must agree with."""
    p = _structure._apply_scratch_probe(4, 3, 2, True, True, False, True)
    assert p['wpar_off'] is not None
    assert p['wpar_esz'] == 1
    assert p['wsg_off'] < p['wpar_off'] < p['cip_off']
    assert p['wpar_off'] % 8 == 0
    offs = _present_offsets(p)
    assert len(offs) == 9, f'expected 9 present regions, got {len(offs)}'
    for i in range(len(offs) - 1):
        assert offs[i][0] < offs[i + 1][0], (
            f'{offs[i][1]}={offs[i][0]} not < {offs[i+1][1]}={offs[i+1][0]}'
        )


def test_probe_wpar_sufficiency_and_total():
    n, e, w = 4, 3, 2
    p = _structure._apply_scratch_probe(n, e, w, True, True, False, True)
    assert p['wsg_off']  + _req(n, p['wsg_esz'])  <= p['wpar_off']
    assert p['wpar_off'] + _req(n, p['wpar_esz']) <= p['cip_off']
    assert p['cip_off']  + _req(e, p['cip_esz'])  <= p['total']
    assert p['total'] % 8 == 0


def test_probe_wpar_zero_counts_keep_room_for_one():
    p = _structure._apply_scratch_probe(0, 0, 0, False, False, False, True)
    assert p['wsg_off'] is None
    assert p['wpar_off'] is not None
    assert p['wpar_off'] + p['wpar_esz'] <= p['cip_off']
