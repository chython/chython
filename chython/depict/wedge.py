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
"""Stereo-bond geometry: solid wedge, hashed wedge, either-bond, and the function that draws
whichever ones the arena stores.

``bond.wedges='stored'`` (the default) honours a stored wedge, but only in the plane it describes;
a centre with no stored wedge is chosen for by ``core.wedge.wedges_for_write``, and
``'recompute'`` discards the stored wedges and chooses for every centre.
"""
from math import ceil, floor, hypot

from ..core import (WEDGE_DOWN, WEDGE_EITHER, WEDGE_NONE, WEDGE_UP, LogRecord,
                    SU_ALLENE, SU_CIS_TRANS, SU_TETRA)
# `_Planar` is `core.wedge`'s proxy for asking a parity function about a layout the molecule does not
# store, which is exactly what a renderer needs; deriving the sign here would be a second copy of the
# convention `core/wedge.py` owns.
from ..core.wedge import _Planar, allene_parity, cis_trans_frame, cis_trans_parity, \
    tetrahedral_parity, wedges_for_write
from .bonds import has_ink, trim_ink, trim_per_end
from .scene import Box, Path, close, curve, line, move, polyline
from .style import DepictStyle


__all__ = ['either_bond', 'hashed_wedge', 'solid_wedge', 'wedge_paths']


Point = tuple[float, float]
Segment = tuple

# Below this a bond has no derivable direction, in molecule units.  Same threshold as bonds.py.
_MIN_LENGTH = 1e-9


def solid_wedge(narrow: Point, wide: Point, width: float) -> tuple[Segment, ...]:
    """A filled triangle: point at ``narrow``, base of ``width`` centred at ``wide``.

    Returns four segments: M (the point), L, L (the base corners), Z.  Filled, not stroked.
    """
    dx, dy = wide[0] - narrow[0], wide[1] - narrow[1]
    length = hypot(dx, dy)
    ux, uy = dx / length, dy / length
    # perpendicular unit vector (90° counter-clockwise)
    px, py = -uy, ux
    hw = width / 2.
    return (
        move(*narrow),
        line(wide[0] + px * hw, wide[1] + py * hw),
        line(wide[0] - px * hw, wide[1] - py * hw),
        close(),
    )


def hashed_wedge(narrow: Point, wide: Point, width: float,
                 step: float) -> list[tuple[Segment, ...]]:
    """Rungs that widen from ``narrow`` toward ``wide``, spaced ``step`` apart.

    Rungs sit at ``i * step`` for i = 1 .. floor(length/step), so no rung falls at the point.
    A rung of zero width at the narrow end is an invisible line and leaves an ugly gap there.
    Each rung is a two-segment tuple: (M point_a, L point_b).
    """
    dx, dy = wide[0] - narrow[0], wide[1] - narrow[1]
    length = hypot(dx, dy)
    ux, uy = dx / length, dy / length
    px, py = -uy, ux
    n = int(floor(length / step))
    rungs = []
    for i in range(1, n + 1):
        t = i * step
        cx = narrow[0] + t * ux
        cy = narrow[1] + t * uy
        hw = (t / length) * (width / 2.)
        rungs.append((
            move(cx + px * hw, cy + py * hw),
            line(cx - px * hw, cy - py * hw),
        ))
    return rungs


def either_bond(p: Point, q: Point, amplitude: float,
                period: float) -> tuple[Segment, ...]:
    """A smooth wavy bond: ``ceil(length/period)`` cubic Bezier segments of alternating sign.

    Each cubic covers ``length/n`` along the bond axis, with handles at ``period/6`` from each
    end -- the handle ratio that makes a sine-like wave.  Starts and ends on the bond axis.
    """
    dx, dy = q[0] - p[0], q[1] - p[1]
    length = hypot(dx, dy)
    ux, uy = dx / length, dy / length
    px, py = -uy, ux
    n = ceil(length / period)
    seg_len = length / n
    handle = period / 6.
    segs: list[Segment] = [move(*p)]
    for i in range(n):
        sign = 1. if i % 2 == 0 else -1.
        t0 = i * seg_len
        t1 = (i + 1) * seg_len
        # start and end on the bond axis; handles offset perpendicular
        x0, y0 = p[0] + t0 * ux, p[1] + t0 * uy
        x1, y1 = p[0] + t1 * ux, p[1] + t1 * uy
        # control points: offset by amplitude in the perpendicular direction
        h1x = x0 + handle * ux + sign * amplitude * px
        h1y = y0 + handle * uy + sign * amplitude * py
        h2x = x1 - handle * ux + sign * amplitude * px
        h2y = y1 - handle * uy + sign * amplitude * py
        segs.append(curve(h1x, h1y, h2x, h2y, x1, y1))
    # ensure the last point is exactly q (floating-point safety)
    last = list(segs[-1])
    last[-2], last[-1] = q[0], q[1]
    segs[-1] = tuple(last)
    return tuple(segs)


def _wedge_trim(p_narrow: Point, p_wide: Point, box_narrow: Box, box_wide: Box,
                trim_clearance: float) -> tuple[Point, Point] | None:
    """Trim a wedge's segment: narrow end to the label edge (clearance 0), wide end with clearance.

    The narrow end is the atom the configuration belongs to, so its point must touch that atom -- or its
    label box edge -- with no extra clearance; a floating tip leaves the reader guessing which atom the
    wedge belongs to.  The wide end follows the ink-conditional rule.
    """
    return trim_per_end(p_narrow, p_wide, box_narrow, box_wide,
                        0., trim_clearance if has_ink(box_wide) else 0.)


def _bond_key(n: int, m: int) -> tuple[int, int]:
    """Low-first bond key -- the same convention as ``bonds.py`` so ``skip=`` works."""
    return (n, m) if n < m else (m, n)


def _covers(unit, drawn) -> bool:
    """Does any wedge in `drawn` state this unit's configuration?  `drawn` is ``[(narrow, wide), ...]``.

    Coverage is per kind, and an allene is not covered through its anchor.

    * **Tetrahedral.**  ``narrow == anchor``, and nothing looser: a wedge whose narrow end is elsewhere
      says nothing about this centre, however close it lies.
    * **Allene.**  The anchor is the centre of the cumulene chain, every bond on it is double, and no
      notation wedges a double bond -- so the chooser puts the wedge on a bond from a chain terminal to
      one of the unit's own ``refs``, where the end the unit names is the *wide* one.
    """
    if unit['kind'] == SU_TETRA:
        return any(narrow == unit['anchor'] for narrow, _ in drawn)
    refs = {r for r in unit['refs'] if r is not None}
    return any(wide in refs and narrow not in refs for narrow, wide in drawn)


def _stored_for_plane(mol, plane, stored, units, log):
    """`stored` minus the wedges of every centre whose stored codes state the wrong configuration here.

    A wedge code renders a configuration in one plane, so a stored wedge drawn into a plane it does not
    describe is a picture of the enantiomer.  This only drops; `new_for_unchosen` then covers the freed
    centre with the chooser's answer, correct for this plane by construction.

    Parity 0 and ``WEDGE_EITHER`` are left alone -- neither states a configuration a plane could
    contradict.  The drop is keyed by ``(narrow, wide)`` pairs and not by narrow alone: an allene's
    narrow end may legitimately carry another centre's wedge on a different bond.
    """
    if not stored:
        return stored
    probe_table = {(n, w): c for n, w, c in stored}

    def probe(a, b):
        return probe_table.get((a, b), WEDGE_NONE)

    planar = _Planar(mol, plane)
    dropped: set[tuple[int, int]] = set()
    for unit in sorted(units, key=lambda u: u['anchor']):
        anchor = unit['anchor']
        if unit['kind'] == SU_TETRA:
            codes = [c for (n, _), c in probe_table.items() if n == anchor]
            if not codes or WEDGE_EITHER in codes:
                continue
            target = mol.parity_of(anchor)
            if not target:
                continue
            read = tetrahedral_parity(planar, unit, None, probe)
            if read == target:
                continue
            dropped.update((n, w) for (n, w) in probe_table if n == anchor)
            if log is not None:
                log.append(LogRecord(
                    'depict:rewedged', (anchor,),
                    f'atom {anchor}: the stored wedge states parity {read} in this layout where the '
                    f'arena holds {target}; redrawn from the chosen wedge rather than asserting the '
                    f'opposite configuration. The wedge was drawn for another layout, not the molecule'))
        elif unit['kind'] == SU_ALLENE:
            refs = {r for r in unit['refs'] if r is not None}
            own = [(n, w, c) for n, w, c in stored if w in refs and n not in refs]
            if not own or WEDGE_EITHER in (c for _, _, c in own):
                continue
            target = mol.parity_of(anchor)
            if not target:
                continue
            read = allene_parity(planar, unit, None, probe)
            if read == target:
                continue
            dropped.update((n, w) for n, w, _ in own)
            if log is not None:
                log.append(LogRecord(
                    'depict:rewedged', (anchor,),
                    f'atom {anchor}: the stored wedge states parity {read} in this layout where the '
                    f'arena holds {target}; redrawn from the chosen wedge rather than asserting the '
                    f'opposite configuration. The wedge was drawn for another layout, not the molecule'))
    if not dropped:
        return stored
    return [(n, w, c) for n, w, c in stored if (n, w) not in dropped]


def _stroke_path(subpaths, width: float, colour: str, bond_style) -> Path:
    """A stroked path with cap/join from the style."""
    return Path(subpaths, stroke=colour, width=width, cap=bond_style.cap,
                join=bond_style.join, miter_limit=bond_style.miter_limit)


def wedge_paths(mol, plane, boxes, style: DepictStyle, *,
                log=None) -> tuple[list[Path], set[tuple[int, int]]]:
    """All stereo-bond paths for ``mol`` at ``plane``, and the bond keys they claimed.

    Returns ``(paths, claimed)``, ``claimed`` being ``(low_id, high_id)`` pairs.  The caller passes it as
    ``bond_paths(skip=...)`` so a wedge bond is not also drawn as a plain line underneath.

    Under ``'stored'``, a stored wedge is kept only where its centre reads back as the parity the arena
    holds in this plane (``depict:rewedged`` otherwise), and ``wedges_for_write`` fills the rest; under
    ``'recompute'`` the stored wedges are discarded.  A configured unit left uncovered logs
    ``depict:unwedged``; a cis/trans unit the plane disagrees with is drawn crossed and logged
    ``depict:crossed``, or ``depict:crossed-contradiction`` where the layout states the opposite, which
    is the actionable one.  ``plane`` reaches ``wedges_for_write`` explicitly, so the chooser answers for
    this layout and not for whatever coordinates the molecule stores.
    """
    bond_style = style.bond
    claimed: set[tuple[int, int]] = set()
    paths: list[Path] = []

    # the whole units are kept, not just their anchors: which atom a wedge has to touch to cover a unit
    # depends on the unit's kind -- see `_covers`
    all_units = mol.stereo_units()
    wedged_units = [u for u in all_units
                    if u['kind'] in (SU_TETRA, SU_ALLENE) and mol.parity_of(u['anchor']) != 0]
    cistrans_units = [u for u in all_units
                      if u['kind'] == SU_CIS_TRANS and mol.parity_of(u['anchor']) != 0]

    if bond_style.wedges == 'recompute':
        chosen, _ = wedges_for_write(mol, plane=plane)
        wedge_list = chosen
    else:  # 'stored'
        stored = _stored_for_plane(mol, plane, list(mol.wedges()), wedged_units, log)
        stored_centres = {narrow for narrow, wide, code in stored}
        # "a bond carries at most one wedge" holds only within one `wedges_for_write` call -- the chooser
        # is never told what is already stored -- so a chosen wedge whose narrow atom is unstored can
        # still land on a bond a stored wedge occupies from the other end.  Exclude by bond as well as by
        # centre; the centre then reports `depict:unwedged`, which is true of the resulting picture.
        stored_bonds = {_bond_key(n, w) for n, w, _ in stored}
        # wedges_for_write always chooses when plane= is given: the stored shortcut is bypassed
        chosen, _ = wedges_for_write(mol, plane=plane)
        new_for_unchosen = [(n, w, c) for n, w, c in chosen
                            if n not in stored_centres and _bond_key(n, w) not in stored_bonds]
        wedge_list = stored + new_for_unchosen

    drawn: list[tuple[int, int]] = []
    for narrow, wide, code in wedge_list:
        p_narrow = plane[narrow]
        p_wide = plane[wide]
        box_narrow = boxes[narrow].box
        box_wide = boxes[wide].box

        segment = _wedge_trim(p_narrow, p_wide, box_narrow, box_wide, bond_style.trim)
        if segment is None:
            continue
        p_n, p_w = segment
        key = _bond_key(narrow, wide)

        if code == WEDGE_UP:
            segs = solid_wedge(p_n, p_w, bond_style.wedge_width)
            path = Path([segs], fill=bond_style.colour)
        elif code == WEDGE_DOWN:
            rungs = hashed_wedge(p_n, p_w, bond_style.wedge_width, bond_style.hash_step)
            if not rungs:
                continue
            path = _stroke_path(rungs, bond_style.width, bond_style.colour, bond_style)
        elif code == WEDGE_EITHER:
            wave = either_bond(p_n, p_w, bond_style.either_amplitude, bond_style.either_period)
            path = _stroke_path([wave], bond_style.width, bond_style.colour, bond_style)
        else:
            continue

        paths.append(path)
        claimed.add(key)
        drawn.append((narrow, wide))

    # configured units no drawn wedge covers
    for unit in sorted(wedged_units, key=lambda u: u['anchor']):
        if not _covers(unit, drawn):
            if log is not None:
                log.append(LogRecord('depict:unwedged', (unit['anchor'],),
                                     f'atom {unit["anchor"]}: no wedge could be drawn for this '
                                     f'stereocentre'))

    # crossed double bonds for cis/trans units the plane does not agree with
    if cistrans_units:
        for unit in cistrans_units:
            anchor = unit['anchor']
            # the test is agreement, not readability: a plain double bond asserts whatever geometry the
            # plane draws, so a layout reading back as cis for a trans molecule draws the wrong compound.
            # Parity 0 never equals a configured parity, so this subsumes the unrepresentable case.  The
            # crossed bond commits to neither, which is the least a picture can assert.
            read = cis_trans_parity(mol, unit, plane=plane)
            stored = mol.parity_of(anchor)
            if read == stored:
                continue
            # the partner is `cis_trans_frame`'s third slot -- the same expression the parity was measured
            # over, so the bond drawn crossed cannot differ from the bond the parity was read on
            frame = cis_trans_frame(mol, unit)
            if frame is None:
                continue   # no frame is also what `cis_trans_parity` returned 0 for
            partner = frame[2]
            if log is not None:
                # two rule ids: an unrepresentable layout is a property of the molecule and the caller can
                # do nothing about it, while a contradicted one means somebody's layout is wrong
                if read == 0:
                    log.append(LogRecord(
                        'depict:crossed', (anchor, partner),
                        f'bond {anchor}-{partner}: this layout cannot represent the stored cis/trans '
                        f'configuration (collinear or degenerate); drawn crossed'))
                else:
                    log.append(LogRecord(
                        'depict:crossed-contradiction', (anchor, partner),
                        f'bond {anchor}-{partner}: this layout draws parity {read} where the arena '
                        f'holds {stored}; drawn crossed rather than asserting the opposite geometry. '
                        f'The layout is wrong, not the molecule'))
            p0, p1 = plane[anchor], plane[partner]
            box_p, box_q = boxes[anchor].box, boxes[partner].box
            # trimmed the same way a plain double bond would be
            axis = trim_ink(p0, p1, box_p, box_q, bond_style.trim, log=log,
                            atoms=(anchor, partner) if anchor < partner else (partner, anchor))
            if axis is None:
                continue
            (ax, ay), (bx, by) = axis
            ldx, ldy = bx - ax, by - ay
            axis_len = hypot(ldx, ldy)
            if axis_len < _MIN_LENGTH:
                continue
            ux, uy = ldx / axis_len, ldy / axis_len
            px_v, py_v = -uy, ux
            half_space = bond_style.spacing / 2.
            # same geometry as an acyclic double bond
            straddle = [
                polyline([(ax + px_v * half_space, ay + py_v * half_space),
                          (bx + px_v * half_space, by + py_v * half_space)]),
                polyline([(ax - px_v * half_space, ay - py_v * half_space),
                          (bx - px_v * half_space, by - py_v * half_space)]),
            ]
            paths.append(_stroke_path(straddle, bond_style.width, bond_style.colour, bond_style))
            # X at the midpoint, arms at 45° to the bond axis, each tip wedge_width/2 from the centre in
            # both the unit and perpendicular directions
            mx, my = (ax + bx) / 2., (ay + by) / 2.
            hw = bond_style.wedge_width / 2.
            x_arms = [
                polyline([(mx - ux * hw - px_v * hw, my - uy * hw - py_v * hw),
                          (mx + ux * hw + px_v * hw, my + uy * hw + py_v * hw)]),
                polyline([(mx + ux * hw - px_v * hw, my + uy * hw - py_v * hw),
                          (mx - ux * hw + px_v * hw, my - uy * hw + py_v * hw)]),
            ]
            paths.append(_stroke_path(x_arms, bond_style.width, bond_style.colour, bond_style))
            claimed.add(_bond_key(anchor, partner))

    return paths, claimed
