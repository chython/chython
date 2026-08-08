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
import re

from chython.algorithms.depict import _dasharray, _render_aromatic_bond, _render_config, depict_settings


_line = re.compile(r'x1="(-?[\d.]+)" y1="(-?[\d.]+)" x2="(-?[\d.]+)" y2="(-?[\d.]+)"')
_dash_re = re.compile(r'stroke-dasharray="([\d.]+) ([\d.]+)"')


def _endpoints(svg):
    x1, y1, x2, y2 = (float(v) for v in _line.search(svg).groups())
    # svg y axis is inverted on write
    return (x1, -y1), (x2, -y2)


def _assert_inside_bond(svg, n, m):
    """inner aromatic dash must stay within the bond footprint"""
    (ax, ay), (bx, by) = _endpoints(svg)
    lo, hi = min(n[0], m[0]), max(n[0], m[0])
    assert lo - 1e-6 <= ax <= hi + 1e-6, f'dash start {ax} escapes bond span [{lo}, {hi}]'
    assert lo - 1e-6 <= bx <= hi + 1e-6, f'dash end {bx} escapes bond span [{lo}, {hi}]'


def test_flat_ring_both_orientations():
    # centroid nearly collinear with the bond (overstretched/flattened ring).
    # must be rejected regardless of which side the centroid falls on.
    n, m = (0., 0.), (1., 0.)
    assert _render_aromatic_bond(*n, *m, .5, .05) is None
    assert _render_aromatic_bond(*n, *m, .5, -.05) is None


def test_dash_inside_bond_both_orientations():
    n, m = (0., 0.), (1., 0.)
    for cy in (.3, -.3, .87, -.87):
        svg = _render_aromatic_bond(*n, *m, .5, cy)
        assert svg is not None
        _assert_inside_bond(svg, n, m)


def test_dash_inside_bond_skewed_centroid():
    # centroid projected far off the bond segment must not slide the dash out
    n, m = (0., 0.), (1., 0.)
    for cx in (2., 5., -4.):
        for cy in (.3, -.3):
            svg = _render_aromatic_bond(*n, *m, cx, cy)
            if svg is not None:
                _assert_inside_bond(svg, n, m)


def test_distorted_rings_keep_dashes_inside():
    """end-to-end: no aromatic dash may escape the ring bounding box, however bad the layout"""
    from random import Random
    from chython import smiles

    rnd = Random(7)
    pad = .25
    for _ in range(2000):
        mol = smiles('c1ccccc1')
        mol.clean2d()
        for a in mol._atoms.values():
            a.x += rnd.uniform(-1.6, 1.6)
            a.y += rnd.uniform(-1.6, 1.6)

        xs = [a.x for a in mol._atoms.values()]
        ys = [a.y for a in mol._atoms.values()]
        lo_x, hi_x = min(xs) - pad, max(xs) + pad
        lo_y, hi_y = min(ys) - pad, max(ys) + pad

        for line in mol.depict().splitlines():
            if 'dasharray' not in line:
                continue
            (ax, ay), (bx, by) = _endpoints(line)
            assert lo_x <= ax <= hi_x and lo_x <= bx <= hi_x, f'dash escapes ring x range: {line.strip()}'
            assert lo_y <= ay <= hi_y and lo_y <= by <= hi_y, f'dash escapes ring y range: {line.strip()}'


def test_normal_aromatics_unchanged():
    from chython import smiles

    # inner dash count for well-formed layouts must not regress
    for smi, expected in (('c1ccccc1', 6), ('c1ccc2ccccc2c1', 12), ('c1cc[nH]c1', 5), ('c1ccncc1', 6)):
        mol = smiles(smi)
        mol.clean2d()
        assert mol.depict().count('stroke-dasharray') == expected, smi


def test_bonds_have_round_linecaps():
    from chython import smiles

    # silent-carbon joints are smoothed by round caps on the bond group
    for smi in ('CCCC', 'c1ccccc1', 'CC(C)C'):
        mol = smiles(smi)
        mol.clean2d()
        for line in mol.depict().splitlines():
            if '<g fill="none" stroke=' in line and 'stroke-width' in line:
                assert 'stroke-linecap="round"' in line, smi
                break
        else:
            raise AssertionError(f'no bond group found for {smi}')


def test_dash_pattern_compensates_round_caps():
    # round caps add bond_width/2 per end. emitted pattern must render at the configured lengths.
    cap = _render_config['bond_width'] / 2
    for dash, gap in ((.15, .05), (.2, .1)):
        emitted = _dash_re.search(_dasharray(dash, gap)).groups()
        e_dash, e_gap = (float(v) for v in emitted)
        assert abs(e_dash + 2 * cap - dash) < .005, f'painted length off for {dash}/{gap}'
        assert abs(e_gap - 2 * cap - gap) < .005, f'gap length off for {dash}/{gap}'


def test_dash_pattern_stays_renderable_on_thick_bonds():
    # caps wider than the dash must not produce a negative/zero dasharray
    try:
        depict_settings(bond_width=.2)
        e_dash = float(_dash_re.search(_dasharray(.15, .05)).group(1))
        assert e_dash > 0., 'dasharray collapsed to non-positive length'
    finally:
        depict_settings()


def test_dash_on_centroid_side():
    # the inner dash belongs between the bond and the ring centroid
    n, m = (0., 0.), (1., 0.)
    for cy in (.5, -.5):
        (_, ay), (_, by) = _endpoints(_render_aromatic_bond(*n, *m, .5, cy))
        assert ay * cy > 0, 'dash offset to the wrong side of the bond'
        assert by * cy > 0, 'dash offset to the wrong side of the bond'
