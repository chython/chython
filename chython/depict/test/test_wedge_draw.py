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
"""Stereo bonds: the only bonds whose shape depends on which end is which.

Named apart from ``formats/ctfile/test/test_wedge.py``, which pytest would collide with.  A stored wedge
is honoured as drawn (``bond.wedges='stored'``, the default); a centre with none is chosen for by
``core.wedge.wedges_for_write``, and ``'recompute'`` discards the stored wedges and chooses for every
centre.  A centre the chooser cannot serve is drawn bare and logged -- not refused, not invented for.
"""
from math import ceil, hypot, sqrt
from pathlib import Path
from pytest import approx, mark
from chython import smiles
from chython.core import SU_ALLENE, SU_CIS_TRANS, SU_TETRA, WEDGE_DOWN, WEDGE_EITHER, WEDGE_NONE, WEDGE_UP
from chython.core.wedge import _Planar, allene_parity, cis_trans_parity, tetrahedral_parity, wedges_for_write
from chython.depict.bonds import has_ink
from chython.depict.label import labels
from chython.depict.style import DepictStyle
from chython.depict.wedge import _stored_for_plane, either_bond, hashed_wedge, solid_wedge, wedge_paths

_STEREO_SDF = Path(__file__).resolve().parents[3] / 'test' / 'stereo.sdf'


def _read_stereo_record(idx):
    """One molecule from ``test/stereo.sdf`` at 0-based ``idx``, or pytest-skip if absent."""
    from pytest import skip
    from chython.formats import SDFRead
    if not _STEREO_SDF.is_file():
        skip('test/stereo.sdf not found')
    with SDFRead(str(_STEREO_SDF)) as f:
        for i, mol in enumerate(f):
            if i == idx:
                return mol
    skip(f'record {idx} not found in stereo.sdf')


def _setup(text, style=None):
    style = style or DepictStyle()
    mol = smiles(text)
    mol.clean2d()
    plane = mol.coordinates()
    return mol, plane, labels(mol, plane, style), style


def _key(n, m):
    """The low-first bond key ``wedge_paths`` returns in `claimed`."""
    return (n, m) if n < m else (m, n)


def _tetra_parity_in(mol, unit, plane, wedges):
    """What `wedges` read back as at `unit`, in `plane`, through the CORE's reader and nothing else.

    ``[(narrow, wide, code), ...]`` in, parity out.  Every claim in this file about a stored wedge
    stating a configuration is measured through here rather than assumed.
    """
    def probe(a, b):
        return next((c for n, w, c in wedges if (n, w) == (a, b)), 0)

    return tetrahedral_parity(_Planar(mol, plane), unit, None, probe)


def _wedged(code, style=None):
    """:func:`_wedged_all` without the pieces most callers here do not look at."""
    paths, _, _, plane, centre, wide, boxes, _ = _wedged_all(code, style)
    return paths, plane, centre, wide, boxes


def _wedged_all(code, style=None):
    """``C[C@H](N)O`` with one stored wedge of `code` on the centre's first bond, ready to draw.

    The wide end is the methyl carbon, so NEITHER end carries a glyph, which makes an unspent trim
    visible.  Returns ``(paths, claimed, log, plane, centre, wide, boxes, mol)``.

    The plane is chosen so `code` STATES the molecule's own configuration: the other code asks for the
    enantiomer, which ``wedge_paths`` declines (``depict:rewedged``).  Reflecting x negates the reader's
    determinant, so one of the mirror pair states the arena's parity; ``WEDGE_EITHER`` states none.
    """
    style = style or DepictStyle()
    mol = smiles('C[C@H](N)O')
    mol.clean2d()
    unit = next(u for u in mol.stereo_units()
                if u['kind'] == SU_TETRA and mol.parity_of(u['anchor']))
    centre = unit['anchor']
    wide = next(iter(mol.neighbors_of(centre)))
    with mol.edit():
        mol.set_wedge(centre, wide, code)
    plane = mol.coordinates()
    if code != WEDGE_EITHER:
        stored = [(centre, wide, code)]
        if _tetra_parity_in(mol, unit, plane, stored) != mol.parity_of(centre):
            plane = {sid: (-x, y) for sid, (x, y) in plane.items()}
        assert _tetra_parity_in(mol, unit, plane, stored) == mol.parity_of(centre), \
            'the fixture must hand over a plane the stored code is TRUE in'
    boxes = labels(mol, plane, style)
    log = []
    paths, claimed = wedge_paths(mol, plane, boxes, style, log=log)
    return paths, claimed, log, plane, centre, wide, boxes, mol


def _collinear_butene(style=None):
    """``C/C=C/C`` on a hand-written COLLINEAR plane, which no cis/trans geometry can be drawn in.

    Every atom on the x-axis, so ``cis_trans_parity`` has no sign to read and the unit is drawn crossed.
    No atom is labelled, so the trimmed axis is the whole bond and the geometry below can be asserted
    against ``plane`` directly.  Returns ``(paths, claimed, log, plane, anchor, partner, mol, unit)``.
    """
    style = style or DepictStyle()
    mol = smiles('C/C=C/C')
    plane = {sid: (float(i), 0.) for i, sid in enumerate(sorted(mol))}
    boxes = labels(mol, plane, style)
    unit = next(u for u in mol.stereo_units()
                if u['kind'] == SU_CIS_TRANS and mol.parity_of(u['anchor']))
    anchor = unit['anchor']
    partner = next(m for m in mol.neighbors_of(anchor) if unit['refs'][2] in mol.neighbors_of(m))
    log = []
    paths, claimed = wedge_paths(mol, plane, boxes, style, log=log)
    return paths, claimed, log, plane, anchor, partner, mol, unit


def _read_the_picture_back(paths, plane, claimed):
    """``[(narrow, wide, code), ...]`` recovered from the DRAWING -- shape for the code, geometry for
    the ends.

    Asking the chooser instead would assert that two functions agree, not that the page is right.  A solid
    wedge is the one FILLED path, an either bond holds cubics, anything else stroked is a fan of rungs.
    The bond is the claimed pair the extremes align with best, not the nearest atom to each, a wedge into
    a labelled atom being trimmed by that label's box.  Wedge-only pictures: a crossed bond's two straddle
    lines would read as a two-rung hashed wedge.
    """
    def distance(point, sid):
        return hypot(plane[sid][0] - point[0], plane[sid][1] - point[1])

    out = []
    for path in paths:
        subpaths = path.subpaths
        if path.fill is not None and path.stroke is None:
            triangle = subpaths[0]
            narrow_point = (triangle[0][1], triangle[0][2])
            base = [(triangle[i][1], triangle[i][2]) for i in (1, 2)]
            wide_point = ((base[0][0] + base[1][0]) / 2., (base[0][1] + base[1][1]) / 2.)
            code = WEDGE_UP
        elif any(seg[0] == 'C' for sub in subpaths for seg in sub):
            first, last = subpaths[0][0], subpaths[0][-1]
            narrow_point = (first[1], first[2])
            wide_point = (last[-2], last[-1])
            code = WEDGE_EITHER
        else:
            rungs = [((s[0][1], s[0][2]), (s[1][1], s[1][2])) for s in subpaths]
            widths = [hypot(a[0] - b[0], a[1] - b[1]) for a, b in rungs]
            middles = [((a[0] + b[0]) / 2., (a[1] + b[1]) / 2.) for a, b in rungs]
            narrow_point = middles[min(range(len(widths)), key=lambda i: widths[i])]
            wide_point = middles[max(range(len(widths)), key=lambda i: widths[i])]
            code = WEDGE_DOWN

        def cost(pair):
            return distance(narrow_point, pair[0]) + distance(wide_point, pair[1])

        bond = min(claimed, key=lambda key: min(cost(key), cost(key[::-1])))
        narrow, wide = bond if cost(bond) <= cost(bond[::-1]) else bond[::-1]
        out.append((narrow, wide, code))
    return out


def _straddle_and_x(paths):
    """The crossed bond's two paths, as point lists: ``([line, line], [arm, arm])``."""
    def points(subpath):
        return [(seg[1], seg[2]) for seg in subpath]

    straddle, x_arms = paths
    return [points(s) for s in straddle.subpaths], [points(s) for s in x_arms.subpaths]


def test_a_solid_wedge_is_a_filled_triangle_pointing_at_the_narrow_end():
    triangle = solid_wedge((0., 0.), (1., 0.), .16)
    assert [s[0] for s in triangle] == ['M', 'L', 'L', 'Z']
    assert (triangle[0][1], triangle[0][2]) == approx((0., 0.)), 'the point is the narrow end'
    wide = [(s[1], s[2]) for s in triangle[1:3]]
    assert hypot(wide[0][0] - wide[1][0], wide[0][1] - wide[1][1]) == approx(.16)


def test_a_hashed_wedge_is_rungs_that_widen_toward_the_wide_end():
    rungs = hashed_wedge((0., 0.), (1., 0.), .16, .1)
    assert len(rungs) > 3
    first = rungs[0]
    last = rungs[-1]
    assert hypot(first[0][1] - first[1][1], first[0][2] - first[1][2]) < \
           hypot(last[0][1] - last[1][1], last[0][2] - last[1][2]), 'the rungs widen'


def test_the_rung_pitch_is_the_style_step():
    rungs = hashed_wedge((0., 0.), (1., 0.), .16, .1)
    assert len(rungs) == approx(10, abs=1), 'a unit bond at a .1 pitch is about ten rungs'


def test_a_hashed_wedge_never_puts_a_rung_at_the_point():
    """a rung of zero width at the narrow end is an invisible line and an ugly gap"""
    rungs = hashed_wedge((0., 0.), (1., 0.), .16, .1)
    assert all(hypot(r[0][1] - r[1][1], r[0][2] - r[1][2]) > 1e-6 for r in rungs)


def test_an_either_bond_is_a_wave_of_cubics():
    wave = either_bond((0., 0.), (1., 0.), .07, .2)
    assert wave[0][0] == 'M'
    assert all(s[0] == 'C' for s in wave[1:]), 'cubics, so it is smooth at every print size'
    assert len(wave) > 3


def test_the_wave_starts_and_ends_on_the_bond_axis():
    wave = either_bond((0., 0.), (1., 0.), .07, .2)
    assert (wave[0][1], wave[0][2]) == approx((0., 0.))
    assert (wave[-1][-2], wave[-1][-1]) == approx((1., 0.), abs=1e-12)


def test_the_wave_ends_EXACTLY_on_the_far_atom():
    """"starts and ends on the axis" is a guarantee, so the assertion is `==` and not a tolerance

    The endpoint snap removes an error of one ULP, which any tolerance would already satisfy.  One
    endpoint cannot show it either: `n * (length / n)` is often exactly `length` and the snap a no-op, so
    a spread is swept and a third of the pairs below are inexact without it.
    """
    for i in range(1, 6):
        for j in range(1, 6):
            q = (i / 7., j / 11.)
            wave = either_bond((0., 0.), q, .07, .2)
            assert (wave[-1][-2], wave[-1][-1]) == q, f'the wave did not land on {q}'


def test_the_wave_alternates_from_side_to_side_of_the_bond_axis():
    """a wave, not a bulge -- and only the SIGN tells the two apart

    Forcing every half-wave to one sign leaves the segment count, the smoothness and the endpoints
    untouched; an "either" bond that does not cross the axis reads as a curved bond.
    """
    wave = either_bond((0., 0.), (1., 0.), .07, .2)
    # p->q is the x-axis, so the perpendicular offset of a control point IS its y.
    handles = [(seg[2], seg[4]) for seg in wave[1:]]
    assert len(handles) > 3
    for first, second in handles:
        assert first == approx(second), 'both handles of one half-wave sit on the same side'
    signs = [1 if first > 0. else -1 for first, _ in handles]
    assert signs == [(-1) ** i * signs[0] for i in range(len(signs))], \
        'successive half-waves must fall on OPPOSITE sides of the axis'


def test_a_centre_with_no_stored_wedge_is_drawn_from_the_chosen_one():
    """'stored' with no stored wedges falls through to assignment, rather than drawing a plain vertex"""
    mol, plane, boxes, style = _setup('C[C@H](N)O')
    assert not list(mol.wedges())
    paths, claimed = wedge_paths(mol, plane, boxes, style)
    assert len(paths) == 1
    assert len(claimed) == 1


def test_the_chosen_wedge_is_the_one_chemistry_chose_and_not_a_second_opinion():
    mol, plane, boxes, style = _setup('C[C@H](N)O')
    _, claimed = wedge_paths(mol, plane, boxes, style)
    wedges, _ = wedges_for_write(mol, plane=plane)
    expected = {(min(n, m), max(n, m)) for n, m, _ in wedges}
    assert claimed == expected


def test_a_molecule_with_no_stereocentre_gets_no_wedge_at_all():
    mol, plane, boxes, style = _setup('CCO')
    assert wedge_paths(mol, plane, boxes, style) == ([], set())


def test_recompute_discards_a_stored_wedge_and_chooses_afresh():
    """for a file whose wedges were drawn against a layout that is not the one being rendered

    The stored wedge goes on a bond the chooser did NOT pick, so the two modes cannot agree by luck.
    """
    mol, plane, boxes, style = _setup('C[C@H](N)O')
    centre = [a.n for a in mol.atoms() if a.parity][0]
    chosen, _ = wedges_for_write(mol, plane=plane)
    chosen_keys = {(min(n, m), max(n, m)) for n, m, _ in chosen}

    candidate_bonds = {(min(centre, n), max(centre, n)) for n in mol.neighbors_of(centre)}
    stored_bond = next(b for b in sorted(candidate_bonds) if b not in chosen_keys)
    # the narrow end is always the stereocentre
    stored_wide = stored_bond[0] if stored_bond[1] == centre else stored_bond[1]
    with mol.edit():
        mol.set_wedge(centre, stored_wide, WEDGE_UP)
    plane = mol.coordinates()
    boxes = labels(mol, plane, style)

    honoured_paths, honoured_claimed = wedge_paths(mol, plane, boxes, style)
    recomputed_paths, recomputed_claimed = wedge_paths(
        mol, plane, boxes, style.tuned(**{'bond.wedges': 'recompute'}))

    assert len(honoured_paths) == len(recomputed_paths) == 1
    assert recomputed_claimed == chosen_keys
    assert honoured_claimed != recomputed_claimed


def test_a_centre_the_chooser_cannot_serve_is_logged_and_drawn_bare():
    """a record, not a refusal, and not an invented wedge

    The chooser is mocked to return ``([], [])`` -- its real shape -- standing in for a centre whose
    geometry it cannot satisfy.
    """
    from unittest.mock import patch

    mol, plane, boxes, style = _setup('C[C@H](N)O')
    log = []
    with patch('chython.depict.wedge.wedges_for_write', return_value=([], [])):
        paths, claimed = wedge_paths(mol, plane, boxes, style, log=log)
    assert paths == [] and claimed == set()
    assert any(record.rule == 'depict:unwedged' for record in log)


def test_a_configured_allene_drawn_correctly_is_not_reported_as_unwedged():
    """penta-2,3-diene, and the wedge is nowhere near the anchor -- which is not a failure

    An allene's anchor is the CENTRE of the cumulene chain and every bond on it is double, so the chooser
    puts the wedge on a bond from a chain TERMINAL: an unwedged check written as `narrow == anchor` can
    never hold for an allene and reports every configured one as dropped.  A false record on the only
    channel that names lost stereochemistry teaches the caller to ignore it, so the log is asserted empty.
    """
    style = DepictStyle()
    mol = smiles('C/C=C=C/C')
    mol.clean2d()
    with mol.edit():
        mol.set_parity(3, 1)
    unit = next(u for u in mol.stereo_units() if u['kind'] == SU_ALLENE)
    assert mol.parity_of(unit['anchor']), 'premise: the allene is configured'

    plane = mol.coordinates()
    chosen, _ = wedges_for_write(mol, plane=plane)
    assert len(chosen) == 1
    narrow, wide, _ = chosen[0]
    assert narrow != unit['anchor'], \
        'premise: the chooser puts an allene wedge on a chain TERMINAL, never on the anchor'
    assert wide in unit['refs'], 'and its wide end is one of the atoms the unit itself names'

    log = []
    paths, claimed = wedge_paths(mol, plane, labels(mol, plane, style), style, log=log)
    assert len(paths) == 1 and claimed == {(min(narrow, wide), max(narrow, wide))}
    assert log == [], f'the allene IS wedged; nothing to report, got {log}'


def test_a_tetrahedral_centre_drawn_with_its_wedge_is_not_reported_either():
    """the negative half of the log contract: a satisfied centre is silent

    An `any(...)` assertion passes just as well when the record is appended unconditionally, so only this
    one can tell a working check from no check.
    """
    mol, plane, boxes, style = _setup('C[C@H](N)O')
    log = []
    paths, claimed = wedge_paths(mol, plane, boxes, style, log=log)
    assert len(paths) == 1 and claimed, 'premise: this centre does get a wedge'
    assert log == [], f'a centre that was drawn is not a centre that could not be, got {log}'


def test_a_stereogenic_but_unconfigured_centre_is_not_reported():
    """1-aminoethanol: stereogenic, parity 0, and NOT the log's business

    An unset parity states no configuration, so a drawing that states none dropped nothing; the channel
    exists to name information the RENDERER lost.
    """
    mol, plane, boxes, style = _setup('CC(N)O')
    stereogenic = [u for u in mol.stereo_units() if u['stereogenic']]
    assert stereogenic, 'premise: there is a stereogenic centre here'
    assert all(mol.parity_of(u['anchor']) == 0 for u in stereogenic), 'premise: none is configured'
    log = []
    wedge_paths(mol, plane, boxes, style, log=log)
    assert log == [], f'an unstated configuration is not a lost one, got {log}'


def test_a_stored_up_wedge_is_drawn_solid_from_its_narrow_end():
    paths, claimed, log, plane, centre, wide, _, _ = _wedged_all(WEDGE_UP)
    assert len(paths) == 1
    assert paths[0].fill is not None and paths[0].stroke is None, 'a solid wedge is filled, not stroked'
    assert claimed == {(min(centre, wide), max(centre, wide))}
    assert log == [], f'a stored wedge that states this layout is reported nowhere, got {log}'


def test_a_stored_down_wedge_is_drawn_hashed():
    paths, _, log, _, _, _, _, _ = _wedged_all(WEDGE_DOWN)
    assert len(paths) == 1
    assert len(paths[0].subpaths) > 3, 'the rungs are subpaths of ONE path'
    assert log == [], f'a stored wedge that states this layout is reported nowhere, got {log}'


def test_an_either_wedge_is_drawn_as_a_wave():
    """and it is NOT rewedged, however the plane reads back

    A wavy bond says "this drawing does not commit", which is true in every plane, so no plane can
    contradict it: replacing it would assert the configuration the record refused to state.
    """
    paths, _, log, plane, centre, wide, _, mol = _wedged_all(WEDGE_EITHER)
    unit = next(u for u in mol.stereo_units() if u['anchor'] == centre)
    assert _tetra_parity_in(mol, unit, plane, [(centre, wide, WEDGE_EITHER)]) == 0 \
        != mol.parity_of(centre), \
        'premise: an either bond reads back as NO configuration where the arena states one'
    assert any(seg[0] == 'C' for p in paths for sub in p.subpaths for seg in sub)
    assert log == [], f'an either bond is not a contradicted wedge, got {log}'


def test_a_stored_wedge_that_does_not_STATE_this_layout_is_replaced_and_logged():
    """a stored code that misstates the passed plane is a picture of the ENANTIOMER

    A wedge code renders a configuration in ONE plane, not a property of a bond: reflect the layout and
    the same code on the same bond states the opposite centre.  Both ways in are ordinary -- ``plane=`` is
    a public keyword on ``depict()``, and ``clean2d(force=True)`` relays a molecule that already has
    wedges.  The assertion is on the PICTURE, read back through the core's own reader.
    """
    style = DepictStyle()
    mol = smiles('C[C@H](N)O')          # alanine's skeleton, minus the acid
    mol.clean2d()
    unit = next(u for u in mol.stereo_units()
                if u['kind'] == SU_TETRA and mol.parity_of(u['anchor']))
    centre = unit['anchor']
    # The wedges a molfile OF THIS LAYOUT would carry -- i.e. codes that are true where they came from.
    stored, _ = wedges_for_write(mol, plane=mol.coordinates())
    assert stored, 'premise: this centre gets a wedge'
    with mol.edit():
        for n, w, c in stored:
            mol.set_wedge(n, w, c)
    # any OTHER layout: a reflection is the cheapest one that still makes the stored codes state the
    # other centre.
    plane = {sid: (-x, y) for sid, (x, y) in mol.coordinates().items()}
    assert _tetra_parity_in(mol, unit, plane, stored) != mol.parity_of(centre), \
        'premise: in THIS plane the stored codes state the opposite configuration'

    log = []
    paths, claimed = wedge_paths(mol, plane, labels(mol, plane, style), style, log=log)
    assert [(r.rule, r.atoms) for r in log] == [('depict:rewedged', (centre,))], \
        f'the substitution is a fact the caller must be able to see, got {log}'
    drawn = _read_the_picture_back(paths, plane, claimed)
    assert _tetra_parity_in(mol, unit, plane, drawn) == mol.parity_of(centre), \
        'the picture must state the configuration the arena holds, not its mirror image'
    assert drawn != stored, 'and it cannot do that by drawing the stored codes'


def test_a_stored_wedge_THAT_states_this_layout_is_drawn_on_its_own_bond():
    """the stored path stays load-bearing: an AGREEING stored wedge the chooser would place elsewhere

    A centre with three heavy neighbours can be stated by wedging any of several bonds, so "agrees with
    the plane" and "is where the chooser would put it" are different properties.  Both are measured, or
    deleting ``stored`` from the merge would pass the suite.
    """
    style = DepictStyle()
    mol = smiles('C[C@H](N)O')
    mol.clean2d()
    plane = mol.coordinates()
    unit = next(u for u in mol.stereo_units()
                if u['kind'] == SU_TETRA and mol.parity_of(u['anchor']))
    centre = unit['anchor']
    chosen, _ = wedges_for_write(mol, plane=plane)
    chosen_halfedges = {(n, w) for n, w, _ in chosen}
    # measured, not guessed: every (bond, code) at this centre that STATES the arena's parity here and
    # that the chooser did not pick
    elsewhere = [(centre, r, code) for r in sorted(mol.neighbors_of(centre))
                 for code in (WEDGE_UP, WEDGE_DOWN)
                 if (centre, r) not in chosen_halfedges
                 and _tetra_parity_in(mol, unit, plane, [(centre, r, code)]) == mol.parity_of(centre)]
    assert elsewhere, 'premise: this centre can be stated on a bond the chooser did not pick'
    stored = [elsewhere[0]]
    with mol.edit():
        mol.set_wedge(*stored[0])

    log = []
    paths, claimed = wedge_paths(mol, plane, labels(mol, plane, style), style, log=log)
    assert log == [], f'a stored wedge that states this layout is reported nowhere, got {log}'
    assert claimed == {_key(*stored[0][:2])}, \
        'the drawing claims the STORED bond and not the one the chooser would have used'
    assert _read_the_picture_back(paths, plane, claimed) == stored, \
        'and it draws the stored code on it, from the stored end'
    assert _tetra_parity_in(mol, unit, plane, stored) == mol.parity_of(centre), \
        'which is a true picture of the molecule -- honouring it costs no correctness'


def test_a_chosen_wedge_never_lands_on_a_bond_a_STORED_wedge_already_carries():
    """``core.wedge._assign``: "a bond carries at most one wedge, ever.  Never relaxed."

    That guarantee holds WITHIN one ``wedges_for_write`` call, and the chooser is never told what is
    already stored, so the merge's exclusion must be keyed on the BOND and not the narrow atom.  A bond
    wedged from both ends is a picture of nothing, and ``bond_paths(skip=claimed)`` suppresses its plain
    line once either way.  Cholesterol is the construction, because one of its eight centres leaves the
    chooser no bond that avoids another configured centre.  Every atom below is found by measurement.
    """
    style = DepictStyle()
    mol = smiles('CC(C)CCC[C@@H](C)[C@H]1CC[C@H]2[C@@H]3CC=C4C[C@@H](O)CC[C@]4(C)[C@H]3CC[C@]12C')
    mol.clean2d()
    plane = mol.coordinates()
    units = {u['anchor']: u for u in mol.stereo_units()
             if u['kind'] == SU_TETRA and mol.parity_of(u['anchor'])}
    chosen, _ = wedges_for_write(mol, plane=plane)
    pair = next(((n, w) for n, w, _ in chosen if w in units), None)
    assert pair is not None, \
        'premise: the chooser puts a wedge on a bond whose WIDE end is another configured centre'
    chooser_narrow, stored_narrow = pair
    code = next((c for c in (WEDGE_UP, WEDGE_DOWN)
                 if _tetra_parity_in(mol, units[stored_narrow], plane,
                                     [(stored_narrow, chooser_narrow, c)])
                 == mol.parity_of(stored_narrow)), None)
    assert code is not None, 'premise: that bond can state the other centre from the other end'
    with mol.edit():
        mol.set_wedge(stored_narrow, chooser_narrow, code)

    log = []
    paths, claimed = wedge_paths(mol, plane, labels(mol, plane, style), style, log=log)
    assert len(paths) == len(claimed), \
        f'one mark per claimed bond: {len(paths)} paths over {len(claimed)} bonds means a bond ' \
        f'carries two wedges'
    assert _key(chooser_narrow, stored_narrow) in claimed, 'premise: the contested bond IS drawn'
    drawn = _read_the_picture_back(paths, plane, claimed)
    assert [t for t in drawn if _key(*t[:2]) == _key(chooser_narrow, stored_narrow)] == \
        [(stored_narrow, chooser_narrow, code)], 'and the mark on it is the STORED one'
    # the excluded centre is genuinely unwedged now
    assert [(r.rule, r.atoms) for r in log] == [('depict:unwedged', (chooser_narrow,))], \
        f'the excluded centre lost its wedge and nothing else is reported, got {log}'
    # every OTHER configured centre is still stated correctly
    for anchor, unit in sorted(units.items()):
        if anchor == chooser_narrow:
            continue
        assert _tetra_parity_in(mol, unit, plane, drawn) == mol.parity_of(anchor), \
            f'centre {anchor} must still read back as the arena holds it'


@mark.parametrize('code', [WEDGE_UP, WEDGE_DOWN, WEDGE_EITHER])
def test_the_narrow_end_is_the_one_the_arena_says_it_is(code):
    """the direction is the whole content of a wedge; getting it backwards inverts the stereocentre

    ALL THREE CODES, because only the solid wedge makes a reversal obvious: a hashed wedge reversed is
    still a fan of rungs, differing only in which end the narrowest sits at, and that is a drawing of
    the OTHER enantiomer rather than an ugly picture.
    """
    paths, plane, centre, wide, _ = _wedged(code)
    assert len(paths) == 1
    subpaths = paths[0].subpaths

    def from_centre(point):
        return hypot(point[0] - plane[centre][0], point[1] - plane[centre][1])

    def from_wide(point):
        return hypot(point[0] - plane[wide][0], point[1] - plane[wide][1])

    if code == WEDGE_DOWN:
        # no point to look at, so the assertion is on the RUNGS: the narrowest is nearest the anchor and
        # they widen away from it; reversed, both of those flip.
        widths = [hypot(rung[0][1] - rung[1][1], rung[0][2] - rung[1][2]) for rung in subpaths]
        distances = [from_centre(((rung[0][1] + rung[1][1]) / 2., (rung[0][2] + rung[1][2]) / 2.))
                     for rung in subpaths]
        assert len(widths) > 3
        nearest = min(range(len(widths)), key=lambda i: distances[i])
        assert widths[nearest] == min(widths), 'the narrowest rung must be the one nearest the anchor'
        assert distances == sorted(distances), 'the rungs run outward from the narrow end'
        assert widths == sorted(widths), 'and they widen as they go'
    else:
        # the solid wedge's point and the wave's start are both the FIRST segment of the first subpath
        start = (subpaths[0][0][1], subpaths[0][0][2])
        assert from_centre(start) < from_wide(start), 'the drawing starts at the atom it is about'


def test_a_wedge_between_two_bare_vertices_touches_both_of_them():
    """no ink, nothing to clear -- and the point must touch its atom in every case

    The narrow end is trimmed to the label box and NO further -- a floating tip leaves the reader unable
    to say which atom the configuration belongs to, so its clearance is `0.` whether or not there is a
    glyph.  The wide end takes `bond.trim`, but only where there is ink to clear.
    """
    paths, plane, centre, wide, boxes = _wedged(WEDGE_UP)
    assert not has_ink(boxes[centre].box) and not has_ink(boxes[wide].box), \
        'premise: a plain carbon carries no glyph, so neither end of this bond has ink'
    triangle = paths[0].subpaths[0]
    point = (triangle[0][1], triangle[0][2])
    base = [(triangle[i][1], triangle[i][2]) for i in (1, 2)]
    base_middle = ((base[0][0] + base[1][0]) / 2., (base[0][1] + base[1][1]) / 2.)
    assert point == approx(plane[centre], abs=1e-12), \
        'the point is a statement about this atom and must touch it'
    assert base_middle == approx(plane[wide], abs=1e-12), \
        'no ink at the wide end, so no clearance to spend there'


def test_a_claimed_bond_is_not_drawn_twice():
    """the wedge draws the bond, so `bonds.py` must skip it -- or the wedge sits on a plain line"""
    from chython.depict.bonds import bond_paths

    mol, plane, boxes, style = _setup('C[C@H](N)O')
    centre = [a.n for a in mol.atoms() if a.parity][0]
    wide = next(iter(mol.neighbors_of(centre)))
    with mol.edit():
        mol.set_wedge(centre, wide, WEDGE_UP)
    plane = mol.coordinates()
    _, claimed = wedge_paths(mol, plane, boxes, style)
    plain = bond_paths(mol, plane, boxes, style, skip=claimed)
    drawn = sum(len(sub) - 1 for p in plain for sub in p.subpaths)
    assert drawn == sum(1 for _ in mol.bonds()) - 1


def test_a_wedge_is_trimmed_at_a_labelled_end():
    """the wide end of a wedge into an OH must stop at the label like any other bond"""
    mol, plane, boxes, style = _setup('C[C@H](O)N')
    centre = [a.n for a in mol.atoms() if a.parity][0]
    oxygen = [a.n for a in mol.atoms() if a.atomic_symbol == 'O'][0]
    with mol.edit():
        mol.set_wedge(centre, oxygen, WEDGE_UP)
    plane = mol.coordinates()
    boxes = labels(mol, plane, style)
    paths, _ = wedge_paths(mol, plane, boxes, style)
    wide_points = [(s[1], s[2]) for s in paths[0].subpaths[0][1:3]]
    for x, y in wide_points:
        assert not (boxes[oxygen].box.min_x <= x <= boxes[oxygen].box.max_x
                    and boxes[oxygen].box.min_y <= y <= boxes[oxygen].box.max_y)


def test_the_wedge_width_comes_from_the_style():
    """`bond.wedge_width` through `wedge_paths`, not just through the pure-geometry function"""
    mol, plane, boxes, _ = _setup('C[C@H](N)O')
    centre = [a.n for a in mol.atoms() if a.parity][0]
    wide = next(iter(mol.neighbors_of(centre)))
    with mol.edit():
        mol.set_wedge(centre, wide, WEDGE_UP)
    plane = mol.coordinates()

    default_style = DepictStyle()
    wide_style = default_style.tuned(**{'bond.wedge_width': .3})

    def _wide_width(style):
        bx = labels(mol, plane, style)
        paths, _ = wedge_paths(mol, plane, bx, style)
        sub = paths[0].subpaths[0]
        p1, p2 = (sub[1][1], sub[1][2]), (sub[2][1], sub[2][2])
        return hypot(p1[0] - p2[0], p1[1] - p2[1])

    assert _wide_width(default_style) == approx(default_style.bond.wedge_width, abs=1e-9)
    assert _wide_width(wide_style) == approx(.3, abs=1e-9)


def test_the_rung_pitch_comes_from_the_style():
    """`hash_step` reaches the drawing, and reaches it as a PITCH

    The pure-geometry test passes literals to `hashed_wedge`, so it says nothing about which style field
    `wedge_paths` hands it.  The assertion is on the spacing between consecutive rungs rather than on
    their count, because a count is satisfied by any pitch within half a rung of the right one.
    """
    default_style = DepictStyle()
    fine_style = default_style.tuned(**{'bond.hash_step': .045})

    def pitch(style):
        paths, plane, centre, wide, _ = _wedged(WEDGE_DOWN, style)
        middles = [((r[0][1] + r[1][1]) / 2., (r[0][2] + r[1][2]) / 2.) for r in paths[0].subpaths]
        gaps = [hypot(b[0] - a[0], b[1] - a[1]) for a, b in zip(middles, middles[1:])]
        assert len(gaps) > 3
        return gaps

    for gap in pitch(default_style):
        assert gap == approx(default_style.bond.hash_step, abs=1e-9)
    for gap in pitch(fine_style):
        assert gap == approx(.045, abs=1e-9)


def test_the_wave_amplitude_and_period_come_from_the_style():
    """both `either_amplitude` and `either_period`, and neither is the other's

    Two fields wired at one call site, so one test: the amplitude is the largest perpendicular offset any
    control point reaches, the period fixes the number of cubics.  Asserting only one would not catch
    the two arguments being swapped.
    """
    default_style = DepictStyle()
    loud_style = default_style.tuned(**{'bond.either_amplitude': .14, 'bond.either_period': .10})

    def measured(style):
        paths, plane, centre, wide, _ = _wedged(WEDGE_EITHER, style)
        wave = paths[0].subpaths[0]
        start = (wave[0][1], wave[0][2])
        end = (wave[-1][-2], wave[-1][-1])
        length = hypot(end[0] - start[0], end[1] - start[1])
        ux, uy = (end[0] - start[0]) / length, (end[1] - start[1]) / length
        px, py = -uy, ux

        def offset(x, y):
            return (x - start[0]) * px + (y - start[1]) * py

        handles = [offset(seg[i], seg[i + 1]) for seg in wave[1:] for i in (1, 3)]
        return length, len(wave) - 1, max(abs(h) for h in handles)

    for style, amplitude, period in ((default_style, default_style.bond.either_amplitude,
                                      default_style.bond.either_period), (loud_style, .14, .10)):
        length, cubics, reach = measured(style)
        assert reach == approx(amplitude, abs=1e-9), 'the wave is as tall as the style says'
        assert cubics == ceil(length / period), 'and there is one half-wave per period of bond'


def test_the_crossed_bond_takes_its_line_spacing_from_the_style():
    """the two straddle lines of a crossed bond are a DOUBLE bond and share `bond.spacing`

    The value is not re-derived here and must not be: a second constant is how two double bonds end up
    looking unlike each other on one page.
    """
    default_style = DepictStyle()
    open_style = default_style.tuned(**{'bond.spacing': .30})

    def separation(style):
        paths, _, _, plane, anchor, partner, _, _ = _collinear_butene(style)
        (line_a, line_b), _ = _straddle_and_x(paths)
        # p->q is the x-axis here, so the two parallel lines differ only in y
        return abs(line_a[0][1] - line_b[0][1])

    assert separation(default_style) == approx(default_style.bond.spacing, abs=1e-9)
    assert separation(open_style) == approx(.30, abs=1e-9)


def test_the_X_arms_are_sized_from_the_style_wedge_width():
    """the X's arms are `wedge_width/2` out along the axis AND across it, so each arm is that times root 2

    Nothing else on the page sets this mark's size, so a hard-coded arm length looks right at the default
    scale and is invisible at a small one.
    """
    default_style = DepictStyle()
    big_style = default_style.tuned(**{'bond.wedge_width': .32})

    def arm_length(style):
        paths, _, _, _, _, _, _, _ = _collinear_butene(style)
        _, arms = _straddle_and_x(paths)
        lengths = [hypot(b[0] - a[0], b[1] - a[1]) for a, b in arms]
        assert lengths[0] == approx(lengths[1], abs=1e-9), 'an X has two arms of one length'
        return lengths[0]

    assert arm_length(default_style) == approx(default_style.bond.wedge_width * sqrt(2.), abs=1e-9)
    assert arm_length(big_style) == approx(.32 * sqrt(2.), abs=1e-9)


def test_the_X_sits_on_the_middle_of_the_bond_it_crosses():
    """the mark means "this bond", and only its position says which bond

    Counting two arms in two subpaths is satisfied by an X drawn at either atom or at the origin.  Both
    arms are asserted, since one midpoint is also the midpoint of two arms drawn on top of each other.
    """
    paths, _, _, plane, anchor, partner, _, _ = _collinear_butene()
    _, arms = _straddle_and_x(paths)
    middle = ((plane[anchor][0] + plane[partner][0]) / 2.,
              (plane[anchor][1] + plane[partner][1]) / 2.)
    for a, b in arms:
        assert ((a[0] + b[0]) / 2., (a[1] + b[1]) / 2.) == approx(middle, abs=1e-9)


def test_a_crossed_double_bond_is_drawn_for_an_unrepresentable_cis_trans_unit():
    """a cis/trans unit the plane cannot represent is drawn as a crossed double bond

    ``C/C=C/C`` on a collinear plane: every atom at y=0, so ``cis_trans_parity`` returns 0 and
    ``wedge_paths`` emits both straddle lines plus an X on the bond midpoint, claims the bond and logs
    ``depict:crossed`` -- the id for a layout that states NOTHING, which only ``clean2d`` can fix.
    """
    paths, claimed, log, plane, anchor, partner, mol, unit = _collinear_butene()
    assert cis_trans_parity(mol, unit, plane=plane) == 0, \
        'premise: a collinear layout must return parity 0 from cis_trans_parity'
    assert (min(anchor, partner), max(anchor, partner)) in claimed

    # two paths: the straddle lines (2 subpaths) and the X (2 subpaths)
    assert len(paths) == 2
    assert len(paths[0].subpaths) == 2, 'both straddle lines as subpaths of one path'
    assert len(paths[1].subpaths) == 2, 'both arms of the X as subpaths of one path'
    assert [(r.rule, r.atoms) for r in log] == [('depict:crossed', (anchor, partner))], \
        'a picture that dropped a stored configuration must say which bond it dropped it on'


def test_a_plane_that_draws_the_configuration_the_arena_holds_is_left_alone():
    """the negative: an agreeing layout gets a plain double bond and no record

    Without this, "crossed when it disagrees" is indistinguishable from "always crossed".
    ``clean2d`` lays trans-2-butene out as trans, so the read parity equals the stored one and there is
    nothing for this module to draw at all.
    """
    mol, plane, boxes, style = _setup('C/C=C/C')
    unit = next(u for u in mol.stereo_units()
                if u['kind'] == SU_CIS_TRANS and mol.parity_of(u['anchor']))
    assert cis_trans_parity(mol, unit, plane=plane) == mol.parity_of(unit['anchor']), \
        'premise: clean2d drew the configuration the arena holds'
    log = []
    paths, claimed = wedge_paths(mol, plane, boxes, style, log=log)
    assert paths == [] and claimed == set(), 'a plain double bond already asserts the right geometry'
    assert log == []


def test_a_plane_that_draws_the_OPPOSITE_configuration_is_crossed_too():
    """a contradicting layout is crossed too, not drawn plain

    ``C/C=C/C`` is trans in the arena, laid out here by hand as CIS.  A plain double bond asserts
    whatever the plane draws, so it would put cis-2-butene on the page silently -- a contradicting plane
    reads back a perfectly good parity, so a guard on parity 0 alone misses it.  The question is
    AGREEMENT, not readability, and the record is a distinct id: ``depict:crossed-contradiction`` means
    somebody's layout is wrong, which is actionable where a collinear one is not.
    """
    style = DepictStyle()
    mol = smiles('C/C=C/C')
    unit = next(u for u in mol.stereo_units()
                if u['kind'] == SU_CIS_TRANS and mol.parity_of(u['anchor']))
    anchor = unit['anchor']
    partner = next(m for m in mol.neighbors_of(anchor) if unit['refs'][2] in mol.neighbors_of(m))
    # both methyls on the SAME side of the C=C axis: a cis drawing of a trans molecule
    plane = {1: (-.5, .87), 2: (0., 0.), 3: (1., 0.), 4: (1.5, .87)}
    read = cis_trans_parity(mol, unit, plane=plane)
    stored = mol.parity_of(anchor)
    assert read and read != stored, \
        'premise: this layout reads back a real parity, and it is not the stored one'

    log = []
    paths, claimed = wedge_paths(mol, plane, labels(mol, plane, style), style, log=log)
    assert (min(anchor, partner), max(anchor, partner)) in claimed
    assert len(paths) == 2 and len(paths[0].subpaths) == len(paths[1].subpaths) == 2, \
        'straddle lines and an X, exactly as for a collinear plane'
    assert ([(r.rule, r.atoms) for r in log]
            == [('depict:crossed-contradiction', (anchor, partner))]), \
        'a layout that contradicts the arena is a different fact from one that states nothing'


def test_wedge_paths_answers_from_the_plane_parameter_not_stored_coordinates():
    """``plane=`` is why ``cis_trans_parity`` grew the keyword; a renderer must use it

    ``clean2d()`` stores TRANS coordinates that agree with the arena, then a hand-written CIS plane is
    passed: the answer must come from the argument.  Falling back to the stored coordinates would find
    agreement and draw a plain double bond, i.e. the wrong compound, silently.  The only test here with
    ``has_coordinates`` True before a contradicting plane is passed, so the only one that can tell.
    """
    style = DepictStyle()
    mol = smiles('C/C=C/C')
    mol.clean2d()
    assert mol.has_coordinates, 'premise: clean2d() must store coordinates'

    unit = next(u for u in mol.stereo_units()
                if u['kind'] == SU_CIS_TRANS and mol.parity_of(u['anchor']))
    anchor = unit['anchor']
    partner = next(m for m in mol.neighbors_of(anchor) if unit['refs'][2] in mol.neighbors_of(m))

    # premise: clean2d drew it trans, so the stored coordinates agree with the arena
    stored_plane = mol.coordinates()
    assert cis_trans_parity(mol, unit, plane=stored_plane) == mol.parity_of(anchor), \
        'premise: clean2d must draw trans-2-butene in the trans configuration'

    # a hand-written CIS plane: both terminal methyls above the C=C axis.  Atom order is parse order.
    a, b, c, d = sorted(mol)
    cis_plane = {a: (-.5, .87), b: (0., 0.), c: (1., 0.), d: (1.5, .87)}
    assert cis_trans_parity(mol, unit, plane=cis_plane) != mol.parity_of(anchor), \
        'premise: the hand-written plane contradicts the stored trans parity'

    log = []
    paths, claimed = wedge_paths(mol, cis_plane, labels(mol, cis_plane, style), style, log=log)
    assert (min(anchor, partner), max(anchor, partner)) in claimed, \
        'the contradicting plane must claim the double bond -- a crossed bond must be drawn'
    assert len(paths) == 2, 'straddle lines plus X, exactly as for any contradicting plane'
    assert any(r.rule == 'depict:crossed-contradiction' for r in log), \
        'wedge_paths must detect the contradiction against the PASSED plane, not stored XY'


def test_an_allene_whose_stored_wedge_does_not_state_this_layout_is_replaced_and_logged():
    """the allene half of the rewedge check, mirroring the tetrahedral test above

    A stored allene wedge renders the configuration in one plane, so reflecting the layout makes the same
    code state the opposite axial centre; ``_stored_for_plane`` handling only ``SU_TETRA`` lets it through
    and draws the enantiomer.  Record 78 of ``test/stereo.sdf`` is one allene unit with no tetrahedral
    centres competing for the same bonds.  The final assertion is on the PICTURE, read back through
    ``allene_parity`` over the wedge the chooser placed after the drop.
    """
    mol = _read_stereo_record(78)
    unit = next(u for u in mol.stereo_units()
                if u['kind'] == SU_ALLENE and mol.parity_of(u['anchor']))
    anchor = unit['anchor']
    refs = {r for r in unit['refs'] if r is not None}
    stored = list(mol.wedges())
    own_pairs = {(n, w) for n, w, c in stored if w in refs and n not in refs}
    assert own_pairs, 'premise: the fixture carries stored wedges for this allene unit'

    # reflecting x negates the determinant, so the stored codes now read back as the opposite parity
    plane = {sid: (-x, y) for sid, (x, y) in mol.coordinates().items()}

    probe_table = {(n, w): c for n, w, c in stored}

    def probe_stored(a, b):
        return probe_table.get((a, b), WEDGE_NONE)

    read_in_mirror = allene_parity(_Planar(mol, plane), unit, None, probe_stored)
    assert read_in_mirror != mol.parity_of(anchor), \
        'premise: the mirrored plane must read back the opposite parity from the stored codes'

    wedged_units = [u for u in mol.stereo_units()
                    if u['kind'] in (SU_TETRA, SU_ALLENE) and mol.parity_of(u['anchor'])]
    log = []
    result = _stored_for_plane(mol, plane, stored, wedged_units, log)

    remaining_allene_pairs = [(n, w) for n, w, c in result if (n, w) in own_pairs]
    assert remaining_allene_pairs == [], \
        f'the allene unit\'s stored pairs must be dropped; still present: {remaining_allene_pairs}'

    assert any(r.rule == 'depict:rewedged' and anchor in r.atoms for r in log), \
        f'the substitution must be logged with the allene anchor; got {log}'

    # ``wedges_for_write`` on the mirrored plane is exactly what ``wedge_paths`` calls after the drop
    chosen, _ = wedges_for_write(mol, plane=plane)
    allene_chosen = [(n, w, c) for n, w, c in chosen if w in refs and n not in refs]
    assert allene_chosen, 'the chooser must place a wedge for the freed allene centre'

    chosen_table = {(n, w): c for n, w, c in allene_chosen}

    def probe(a, b):
        return chosen_table.get((a, b), WEDGE_NONE)

    read_chosen = allene_parity(_Planar(mol, plane), unit, None, probe)
    assert read_chosen == mol.parity_of(anchor), \
        (f'the chosen wedge must state the arena parity {mol.parity_of(anchor)} in the mirrored '
         f'plane; got {read_chosen}')


def test_an_allene_whose_stored_wedge_does_state_this_layout_is_left_alone():
    """The stored path stays load-bearing: an agreeing allene wedge is not dropped.

    Without this, "drop every allene wedge" passes the preceding test.  Same fixture, own coordinates:
    the stored codes were drawn for THIS plane, so ``allene_parity`` reads back the arena parity.
    """
    mol = _read_stereo_record(78)
    unit = next(u for u in mol.stereo_units()
                if u['kind'] == SU_ALLENE and mol.parity_of(u['anchor']))
    anchor = unit['anchor']
    refs = {r for r in unit['refs'] if r is not None}
    stored = list(mol.wedges())
    own_pairs = {(n, w) for n, w, c in stored if w in refs and n not in refs}
    assert own_pairs, 'premise: the fixture carries stored wedges for this allene unit'

    plane = mol.coordinates()

    probe_table = {(n, w): c for n, w, c in stored}

    def probe_stored(a, b):
        return probe_table.get((a, b), WEDGE_NONE)

    read_in_own = allene_parity(_Planar(mol, plane), unit, None, probe_stored)
    assert read_in_own == mol.parity_of(anchor), \
        'premise: the stored codes must state the arena parity in the molecule\'s own coordinates'

    wedged_units = [u for u in mol.stereo_units()
                    if u['kind'] in (SU_TETRA, SU_ALLENE) and mol.parity_of(u['anchor'])]
    log = []
    result = _stored_for_plane(mol, plane, stored, wedged_units, log)

    remaining_allene_pairs = {(n, w) for n, w, c in result if (n, w) in own_pairs}
    assert remaining_allene_pairs == own_pairs, \
        f'an agreeing stored wedge must survive; dropped: {own_pairs - remaining_allene_pairs}'

    assert not any(r.rule == 'depict:rewedged' and anchor in r.atoms for r in log), \
        f'a correct stored wedge must not be logged as rewedged; got {log}'


def test_a_stored_wedge_on_an_unconfigured_centre_survives_and_nothing_is_logged_about_it():
    """parity 0 is not a claim to verify, so no plane can contradict it

    Comparing a well-drawn layout's read-back (non-zero) against the arena's 0 finds them unequal and
    drops a wedge the file placed, logging a contradiction that does not exist.  Record 15 of
    ``test/stereo.sdf`` is a dichloro diacid with three stereo units: two configured (parities 2 and 1)
    and one unconfigured that still carries a stored wedge.  The unconfigured atom is found by
    ``parity_of(anchor) == 0``, not by a literal id.
    """
    mol = _read_stereo_record(15)
    stored = list(mol.wedges())
    plane = mol.coordinates()  # own coordinates — the bug fires here, not only in a mirror

    # ALL units, not the parity-filtered ones every caller passes: the function promises this itself
    all_tetra = [u for u in mol.stereo_units() if u['kind'] in (SU_TETRA, SU_ALLENE)]
    assert any(not mol.parity_of(u['anchor']) for u in all_tetra), \
        'premise: at least one unit in this record is unconfigured'

    log = []
    result = _stored_for_plane(mol, plane, stored, all_tetra, log)

    assert result == stored, \
        f'all stored wedges must survive in own coordinates; dropped: {set(stored) - set(result)}'
    assert log == [], \
        f'own coordinates with no contradiction means no log lines; got {log}'


def test_the_unconfigured_guard_does_not_exempt_configured_centres_from_being_dropped():
    """the same record, mirrored plane -- configured centres still get dropped

    Without this, returning `stored` unconditionally passes the previous test.  Only the parity guard
    distinguishes the two configured anchors, which must be dropped, from the unconfigured one.
    """
    mol = _read_stereo_record(15)
    stored = list(mol.wedges())
    plane = {sid: (-x, y) for sid, (x, y) in mol.coordinates().items()}  # mirror

    all_tetra = [u for u in mol.stereo_units() if u['kind'] in (SU_TETRA, SU_ALLENE)]
    configured = {u['anchor'] for u in all_tetra if mol.parity_of(u['anchor'])}
    unconfigured = {u['anchor'] for u in all_tetra if not mol.parity_of(u['anchor'])}

    log = []
    result = _stored_for_plane(mol, plane, stored, all_tetra, log)

    for anchor in sorted(configured):
        assert not any(n == anchor for n, w, c in result), \
            f'configured anchor {anchor} must be dropped in a contradicting plane'
    rewedged_anchors = {r.atoms[0] for r in log if r.rule == 'depict:rewedged'}
    assert configured == rewedged_anchors, \
        f'every configured anchor must be logged, and no other; got {rewedged_anchors}'

    for anchor in sorted(unconfigured):
        assert any(n == anchor for n, w, c in result), \
            f'unconfigured anchor {anchor} must survive regardless of plane'
    assert not any(r.rule == 'depict:rewedged' and a in r.atoms
                   for r in log for a in unconfigured), \
        f'unconfigured anchors must never appear in rewedged records'


def test_a_centre_with_two_stored_wedges_has_both_dropped_when_they_misstate_the_layout():
    """EVERY pair of a dropped centre, not just the first: a half-dropped centre states two
    configurations at once.

    The drop accumulator holds ``(narrow, wide)`` pairs so an allene drop cannot remove an unrelated
    tetrahedral wedge sharing a narrow atom, and the tetrahedral branch must contribute all of them.
    Record 21 of ``test/stereo.sdf`` is 2-deuterio-2-butanol: one stereocentre (anchor 6, parity 2) with
    exactly two stored wedges ``(6, 3, 2)`` and ``(6, 4, 1)``, both contradicted by the mirrored plane.
    """
    mol = _read_stereo_record(21)
    stored = list(mol.wedges())
    plane = {sid: (-x, y) for sid, (x, y) in mol.coordinates().items()}  # mirror

    configured_units = [u for u in mol.stereo_units()
                        if u['kind'] in (SU_TETRA, SU_ALLENE) and mol.parity_of(u['anchor'])]
    anchor_units = [u for u in configured_units if u['anchor'] == 6]
    assert len(anchor_units) == 1, 'premise: anchor 6 must be a single configured centre'
    anchor = anchor_units[0]['anchor']

    anchor_pairs = {(n, w) for n, w, c in stored if n == anchor}
    assert len(anchor_pairs) == 2, \
        f'premise: anchor 6 must carry exactly two stored wedges; got {sorted(anchor_pairs)}'

    log = []
    result = _stored_for_plane(mol, plane, stored, configured_units, log)

    remaining = {(n, w) for n, w, c in result if n == anchor}
    assert remaining == set(), \
        f'BOTH stored pairs of anchor {anchor} must be dropped; still present: {remaining}'
    assert any(r.rule == 'depict:rewedged' and anchor in r.atoms for r in log), \
        f'the drop must be logged'
