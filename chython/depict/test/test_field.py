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
"""The scalar field, and the four stages that turn it into smooth closed curves.

Two rules held throughout: an atom with NO value contributes nothing and is not zero (`at()` returns None,
which is drawn as nothing), and a traced vertex is REFINED before it is fitted -- marching squares
interpolates linearly, so its vertices sit off the isoline by O(spacing^2) and the fitted curve wobbles
where two adjacent cells err in opposite directions.  Two Newton steps along the gradient remove it.
"""
from math import exp, hypot, isclose
from pytest import approx, raises
from chython.depict.field import (ScalarField, contour_levels, convex_hull, in_polygon, isolines,
                                  refine, sample, to_cubics, trim_asymptote, _saddle_segments)


def _gaussian_bump():
    """one atom at the origin with value 1 -- the field is then a known radial function"""
    return ScalarField({1: 1.}, {1: (0., 0.)}, sigma=1., cutoff=4.)


def test_the_field_at_an_atom_is_that_atoms_value():
    field = _gaussian_bump()
    assert field.at(0., 0.) == approx(1.)


def test_a_two_atom_field_is_the_sum_of_two_decaying_bumps():
    """f = Σ vᵢ·exp(−‖p−pᵢ‖²/2σ²) — a plain sum, no denominator.

    Normalising would give 0.5 at the midpoint; the un-normalised sum gives exp(−0.5) ≈ 0.607.  Three
    exact points, so a silent re-normalisation cannot sneak back in.
    """
    field = ScalarField({1: 0., 2: 1.}, {1: (0., 0.), 2: (2., 0.)}, sigma=1., cutoff=4.)
    # midpoint: atom1 (v=0) at d=1, atom2 (v=1) at d=1  → 0·exp(−0.5) + 1·exp(−0.5) = exp(−0.5)
    assert field.at(1., 0.) == approx(exp(-0.5))
    # atom1's position: only atom2 contributes (d=2) → 0·exp(0) + 1·exp(−2) = exp(−2)
    assert field.at(0., 0.) == approx(exp(-2.))
    # atom2's position: only atom1 contributes (d=2) → 0·exp(−2) + 1·exp(0) = 1
    assert field.at(2., 0.) == approx(1.)
    # asymmetry: closer to atom2 (v=1) means higher value; closer to atom1 (v=0) means lower
    assert field.at(.5, 0.) < .5
    assert field.at(1.5, 0.) > .5


def test_the_field_is_none_beyond_the_cutoff():
    field = ScalarField({1: 1.}, {1: (0., 0.)}, sigma=1., cutoff=2.)
    assert field.at(5., 0.) is None


def test_an_unvalued_atom_contributes_nothing_and_is_not_zero():
    """the defect this prevents: a partial charge map showing a neutral ring that is really no data"""
    named = ScalarField({1: 1.}, {1: (0., 0.), 2: (1., 0.)}, sigma=1., cutoff=4.)
    alone = ScalarField({1: 1.}, {1: (0., 0.)}, sigma=1., cutoff=4.)
    assert named.at(.9, 0.) == approx(alone.at(.9, 0.))


def test_a_value_naming_an_atom_that_is_not_in_the_plane_is_refused():
    with raises(ValueError, match='7'):
        ScalarField({7: 1.}, {1: (0., 0.)}, sigma=1., cutoff=4.)


def test_the_gradient_is_analytic_and_matches_finite_differences():
    """analytic because Newton refinement calls it per vertex per step, and because a finite difference
    at the .12 grid spacing is itself the error we are trying to remove"""
    field = ScalarField({1: 0., 2: 1.}, {1: (0., 0.), 2: (2., 0.)}, sigma=1., cutoff=4.)
    h = 1e-6
    for x, y in ((.7, .3), (1.4, -.6), (1., 0.)):
        gx, gy = field.gradient(x, y)
        assert gx == approx((field.at(x + h, y) - field.at(x - h, y)) / (2 * h), abs=1e-4)
        assert gy == approx((field.at(x, y + h) - field.at(x, y - h)) / (2 * h), abs=1e-4)


def test_the_gradient_of_a_flat_field_is_zero():
    field = ScalarField({1: .5, 2: .5}, {1: (0., 0.), 2: (2., 0.)}, sigma=1., cutoff=4.)
    assert field.gradient(1., 0.) == approx((0., 0.), abs=1e-9)


def test_bounds_cover_every_named_atom_plus_the_pad():
    field = ScalarField({1: 1., 2: 1.}, {1: (0., 0.), 2: (2., 1.)}, sigma=1., cutoff=4.)
    box = field.bounds(.5)
    assert (box.min_x, box.min_y, box.max_x, box.max_y) == approx((-.5, -.5, 2.5, 1.5))


def test_sampling_covers_the_box_and_records_its_own_origin():
    field = _gaussian_bump()
    grid = sample(field, field.bounds(1.), .25)
    assert grid.spacing == approx(.25)
    assert grid.nx * grid.ny == len(grid.z)
    assert grid.min_x + (grid.nx - 1) * grid.spacing >= field.bounds(1.).max_x


def test_sampling_records_none_outside_the_cutoff():
    field = ScalarField({1: 1.}, {1: (0., 0.)}, sigma=.3, cutoff=.5)
    grid = sample(field, field.bounds(2.), .25)
    assert any(v is None for v in grid.z)
    assert any(v is not None for v in grid.z)


def test_an_isoline_of_a_radial_bump_is_one_closed_ring():
    field = _gaussian_bump()
    grid = sample(field, field.bounds(1.5), .1)
    lines = isolines(grid, .5)
    assert len(lines) == 1
    ring = lines[0]
    assert ring[0] == approx(ring[-1]), 'a ring is closed by repeating its first point'


def test_the_isoline_of_a_radial_bump_is_a_circle_of_the_right_radius():
    """f(r) = exp(-r^2/2) for one atom, so f = .5 at r = sqrt(2 ln 2) ~= 1.1774"""
    field = _gaussian_bump()
    grid = sample(field, field.bounds(1.5), .1)
    ring = refine(field, isolines(grid, .5)[0], .5, 2, spacing=.1)
    radii = [hypot(x, y) for x, y in ring]
    assert min(radii) == approx(1.1774, abs=.01)
    assert max(radii) == approx(1.1774, abs=.01)


def test_refinement_moves_every_vertex_onto_the_isoline():
    field = _gaussian_bump()
    grid = sample(field, field.bounds(1.5), .2)          # deliberately coarse
    raw = isolines(grid, .5)[0]
    fine = refine(field, raw, .5, 2, spacing=.2)
    assert max(abs(field.at(x, y) - .5) for x, y in fine) < \
           max(abs(field.at(x, y) - .5) for x, y in raw) / 10, 'an order of magnitude, at least'


def test_refinement_is_idempotent_on_an_already_exact_vertex():
    field = _gaussian_bump()
    exact = [(1.1774, 0.)]
    assert refine(field, exact, .5, 2, spacing=.1)[0] == approx(exact[0], abs=1e-4)


def test_a_newton_step_longer_than_one_grid_cell_is_abandoned():
    """a step longer than one cell is not refinement, so the traced point is kept EXACTLY as it was

    Newton diverges where the gradient is small and the residual is not: at r = 5.6σ, f = 1.5e-7 and
    |∇f| = 8.7e-7, so one unguarded step is 5.8e6 units long.  The seed must also clear the pre-existing
    |∇f| guard -- here |∇f|² = 7.5e-13, well above its 1e-18 -- or that guard, not this one, holds it.
    """
    field = ScalarField({1: 1.}, {1: (0., 0.)}, sigma=1., cutoff=6.)
    gx, gy = field.gradient(5.6, 0.)
    assert gx * gx + gy * gy > 1e-18, 'the seed must not be caught by the vanishing-gradient guard'
    assert refine(field, [(5.6, 0.)], -5., 2, spacing=.1) == [(5.6, 0.)]


def test_refinement_does_not_move_a_vertex_where_the_gradient_vanishes():
    """dividing by |grad| is the one way this can produce a NaN and put garbage in the SVG"""
    field = ScalarField({1: .5, 2: .5}, {1: (0., 0.), 2: (2., 0.)}, sigma=1., cutoff=4.)
    assert refine(field, [(1., 0.)], .4, 2, spacing=.1)[0] == approx((1., 0.))


def test_two_separated_atoms_give_two_rings_at_a_high_level():
    field = ScalarField({1: 1., 2: 1.}, {1: (0., 0.), 2: (6., 0.)}, sigma=1., cutoff=4.)
    grid = sample(field, field.bounds(1.5), .15)
    assert len(isolines(grid, .5)) == 2


def test_a_saddle_cell_is_resolved_by_its_centre_value_not_arbitrarily():
    """the marching-squares ambiguity: cases 5 and 10 have two valid connections

    This field and grid produce zero case-5 and zero case-10 cells, so the resolution itself is guarded by
    `test_saddle_segments_all_four_cases` instead.
    """
    field = ScalarField({1: 1., 2: 1., 3: -1., 4: -1.},
                        {1: (0., 0.), 2: (2., 2.), 3: (2., 0.), 4: (0., 2.)}, sigma=1., cutoff=6.)
    grid = sample(field, field.bounds(1.), .1)
    lines = isolines(grid, 0.)
    assert lines, 'the zero level of a saddle exists'
    for line in lines:
        for x, y in line:
            assert abs(field.at(x, y)) < .2, 'and every vertex is near it'


def test_saddle_segments_all_four_cases():
    """_saddle_segments is tested directly because the end-to-end path cannot detect it

    Over 35 box × spacing combinations of the saddle field above, `len(isolines(...))` is 2 under the
    correct resolution, the inverted one and a corner-based one alike: one cell in thousands changes a
    segment pairing, not a count.  That is the justification for importing a private name.
    """
    def fs(pairs):
        return {frozenset(p) for p in pairs}

    # case 5: BL and TR above
    # centre above: isolated corners are BR and TL (below) → edges {0,1} and {2,3}
    assert fs(_saddle_segments(5, 1., 0.)) == {frozenset({0, 1}), frozenset({2, 3})}
    # centre below: isolated corners are BL and TR (above) → edges {0,3} and {1,2}
    assert fs(_saddle_segments(5, -1., 0.)) == {frozenset({0, 3}), frozenset({1, 2})}
    # case 10: BR and TL above
    # centre above: isolated corners are BL and TR (below) → edges {0,3} and {1,2}
    assert fs(_saddle_segments(10, 1., 0.)) == {frozenset({0, 3}), frozenset({1, 2})}
    # centre below: isolated corners are BR and TL (above) → edges {0,1} and {2,3}
    assert fs(_saddle_segments(10, -1., 0.)) == {frozenset({0, 1}), frozenset({2, 3})}


def test_an_open_isoline_is_not_closed():
    """a contour that runs off the sampled box is left OPEN; closing one invents a wall.

    The .5 level of this bump is a circle at r = 1.1774 and the box is a square of half-width 1., so the
    circle crosses all four edges and comes back as four corner arcs.  The count is asserted because an
    empty list would otherwise pass with the loop body never running.
    """
    field = ScalarField({1: 1.}, {1: (0., 0.)}, sigma=1., cutoff=4.)
    grid = sample(field, field.bounds(1.), .1)           # a box INSIDE the .5 isoline's radius
    lines = isolines(grid, .5)
    assert len(lines) == 4, 'four corner arcs, not one ring and not nothing'
    for line in lines:
        assert line[0] != approx(line[-1])


def _nodal_field():
    """two atoms of opposite sign, at the shipped sigma -- so level 0 is a real nodal line AND an asymptote

    The numbers are `Oc1ccccc1`'s at a `clean2d()` layout, measured: mean bond .825, so the effective
    sigma is .55 * .825 = .4537, and the two atoms are the phenol oxygen and its ring carbon.
    """
    return ScalarField({1: -.4, 2: .6}, {1: (0., 0.), 2: (.825, 0.)}, sigma=.4537, cutoff=4.)


def test_the_asymptotic_tail_of_a_nodal_line_is_trimmed_and_the_nodal_line_is_kept():
    """level 0 of a mixed-sign field is ONE chain with two natures, and only one of them is a contour

    Measured here: 90 traced vertices, |grad f| running 7.1e-17 .. 1.27 and the distance to the nearer
    atom 0.31 .. 3.94.  Near the atoms it is the nodal line where the two Gaussians cancel; out at 3.9 the
    field is ~4e-18 and floating point decides the sign.  Discarding the chain loses the nodal line.
    """
    field = _nodal_field()
    grid = sample(field, field.bounds(4.), .12)
    chains = isolines(grid, 0.)
    assert len(chains) == 1, 'the premise: one open chain, not a ring and not fragments'
    traced = chains[0]
    assert len(traced) > 60

    runs = trim_asymptote(field, traced)
    kept = [v for run in runs for v in run]
    assert runs, 'the nodal line is a contour and must survive'
    assert len(kept) < len(traced), 'and the asymptote is not one, so something must go'

    # the two sets do not overlap in |grad f| at all
    strongest = max(hypot(*field.gradient(x, y)) for x, y in traced)
    dropped = [v for v in traced if v not in kept]
    assert min(hypot(*field.gradient(x, y)) for x, y in kept) > 1e-6 * strongest
    assert max(hypot(*field.gradient(x, y)) for x, y in dropped) <= 1e-6 * strongest

    # and the trimmed curve no longer runs out to the cutoff boundary
    def near(x, y):
        return min(hypot(x, y), hypot(x - .825, y))

    assert max(near(x, y) for x, y in traced) > 3.9, 'the premise: the traced chain reaches the asymptote'
    assert max(near(x, y) for x, y in kept) < 2.5, 'the drawn one does not'


def test_a_well_conditioned_ring_is_returned_whole_so_its_closure_survives():
    """the trim must be inert on a contour that is one: a ring handed back in pieces cannot be filled, and
    a band that loses its fill because of a guard against noise is the guard doing the damage"""
    field = _gaussian_bump()
    grid = sample(field, field.bounds(1.5), .1)
    ring = isolines(grid, .5)[0]
    assert ring[0] == ring[-1], 'the premise: a closed ring'
    assert trim_asymptote(field, ring) == [list(ring)], 'one run, every vertex, same order'


def test_a_chain_where_the_field_does_not_vary_at_all_is_not_a_contour():
    """all-zero values: level 0 is "satisfied" by every point of the plane, so the trace is an artefact.
    The scale a relative floor divides by is zero here, so the case is answered explicitly, not by 0/0."""
    flat = ScalarField({1: 0., 2: 0.}, {1: (0., 0.), 2: (1., 0.)}, sigma=.5, cutoff=4.)
    assert trim_asymptote(flat, [(0., 0.), (.5, 0.), (1., 0.)]) == []


def test_a_single_vertex_chain_trims_to_nothing_rather_than_to_a_bare_move():
    field = _gaussian_bump()
    assert trim_asymptote(field, [(1., 0.)]) == []
    assert trim_asymptote(field, []) == []


def test_a_pad_that_would_invert_the_box_is_refused_where_the_mistake_is():
    """an inverted Box does not raise: it reaches `sample` as a degenerate grid and comes back as an empty
    contour list, which reads as "this level is not in the field" -- the wrong answer, silently"""
    field = ScalarField({1: 1., 2: 1.}, {1: (0., 0.), 2: (2., 2.)}, sigma=1., cutoff=4.)
    field.bounds(-.9)                                   # the atoms span 2 × 2, so -.9 is still a box
    with raises(ValueError, match='inverts the bounding box'):
        field.bounds(-1.2)


def test_cubics_pass_through_every_polyline_vertex():
    field = _gaussian_bump()
    grid = sample(field, field.bounds(1.5), .15)
    ring = refine(field, isolines(grid, .5)[0], .5, 2, spacing=.15)
    segs = to_cubics(field, ring, .5, closed=True)
    on_curve = [(segs[0][1], segs[0][2])] + [(s[5], s[6]) for s in segs[1:] if s[0] == 'C']
    for x, y in on_curve:
        assert field.at(x, y) == approx(.5, abs=1e-3)


def test_a_closed_contour_emits_a_close_command_and_no_duplicate_point():
    field = _gaussian_bump()
    grid = sample(field, field.bounds(1.5), .15)
    ring = refine(field, isolines(grid, .5)[0], .5, 2, spacing=.15)
    segs = to_cubics(field, ring, .5, closed=True)
    assert segs[-1] == ('Z',)
    assert segs[0][0] == 'M'


def test_the_tangent_at_each_vertex_is_perpendicular_to_the_gradient():
    """this is what "smooth" means here: the curve's direction is the isoline's direction, so two
    adjacent cubics meet with matching tangents and there is no visible crease"""
    field = _gaussian_bump()
    grid = sample(field, field.bounds(1.5), .15)
    ring = refine(field, isolines(grid, .5)[0], .5, 2, spacing=.15)
    segs = to_cubics(field, ring, .5, closed=True)
    for previous, current in zip(segs[1:-2], segs[2:-1]):
        if previous[0] != 'C' or current[0] != 'C':
            continue
        joint = (previous[5], previous[6])
        incoming = (joint[0] - previous[3], joint[1] - previous[4])
        outgoing = (current[1] - joint[0], current[2] - joint[1])
        cross = incoming[0] * outgoing[1] - incoming[1] * outgoing[0]
        dot = incoming[0] * outgoing[0] + incoming[1] * outgoing[1]
        assert dot > 0., 'the handles are collinear and same-signed: G1 continuity'
        assert abs(cross) < 1e-6 * max(1e-9, dot), 'and exactly collinear, not just nearly'
        gx, gy = field.gradient(*joint)
        assert abs(outgoing[0] * gx + outgoing[1] * gy) < 1e-3 * hypot(gx, gy), \
            'perpendicular to the gradient'


def test_a_three_vertex_ring_still_produces_a_valid_closed_path():
    """the degenerate case at a high level in a coarse grid, where a whole ring is three cells"""
    field = _gaussian_bump()
    segs = to_cubics(field, [(.2, 0.), (0., .2), (-.2, 0.), (.2, 0.)], .96, closed=True)
    assert segs[0][0] == 'M' and segs[-1] == ('Z',)


def test_a_two_vertex_line_degrades_to_a_straight_cubic_rather_than_failing():
    field = _gaussian_bump()
    segs = to_cubics(field, [(1.1774, 0.), (1.1, .4)], .5, closed=False)
    assert len(segs) == 2 and segs[1][0] == 'C'


def test_no_zero_length_handle_in_any_contour():
    """guards the chord < 1e-9 skip in to_cubics

    A closed ring whose trace starts beside an almost co-located predecessor gives a wrap-around chord of
    ~1e-16; without the skip the handle length underflows to exactly 0.0 in float64, the G1 dot product is
    0.0, and any caller dividing by handle length divides by zero.  Swept over three fields.
    """
    def _no_zero_handles(segs):
        prev_end = (segs[0][1], segs[0][2])   # M coordinates
        for seg in segs[1:]:
            if seg[0] != 'C':
                continue
            cp1x, cp1y = seg[1], seg[2]
            hx = cp1x - prev_end[0]
            hy = cp1y - prev_end[1]
            assert hx != 0. or hy != 0., (
                f'zero-length outgoing handle at ({prev_end}) → cp1 ({cp1x},{cp1y}); '
                f'chord < 1e-9 skip missing or disabled')
            prev_end = (seg[5], seg[6])

    # the wrap-around degenerate span lives in the radial bump
    field1 = _gaussian_bump()
    grid1 = sample(field1, field1.bounds(1.5), .15)
    ring1 = refine(field1, isolines(grid1, .5)[0], .5, 2, spacing=.15)
    _no_zero_handles(to_cubics(field1, ring1, .5, closed=True))

    # two separated atoms
    field2 = ScalarField({1: 1., 2: 1.}, {1: (0., 0.), 2: (6., 0.)}, sigma=1., cutoff=4.)
    grid2 = sample(field2, field2.bounds(1.5), .15)
    for ring in isolines(grid2, .5):
        _no_zero_handles(to_cubics(field2, refine(field2, ring, .5, 2, spacing=.15), .5,
                                   closed=(ring[0] == ring[-1])))

    # saddle field
    field3 = ScalarField({1: 1., 2: 1., 3: -1., 4: -1.},
                         {1: (0., 0.), 2: (2., 2.), 3: (2., 0.), 4: (0., 2.)}, sigma=1., cutoff=6.)
    grid3 = sample(field3, field3.bounds(1.), .1)
    for ring in isolines(grid3, 0.):
        _no_zero_handles(to_cubics(field3, refine(field3, ring, 0., 2, spacing=.1), 0.,
                                   closed=(ring[0] == ring[-1])))


def test_levels_are_evenly_spaced_inside_the_domain_and_exclude_the_extremes():
    """a contour AT the maximum is a point, and one at the minimum is the whole box outline"""
    field = ScalarField({1: -1., 2: 1.}, {1: (0., 0.), 2: (2., 0.)}, sigma=1., cutoff=4.)
    levels = contour_levels(5, -1., 1.)
    assert len(levels) == 5
    assert all(-1. < level < 1. for level in levels)
    gaps = [b - a for a, b in zip(levels, levels[1:])]
    assert all(gap == approx(gaps[0]) for gap in gaps)


def test_zero_levels_is_allowed_and_yields_nothing():
    field = _gaussian_bump()
    assert contour_levels(0, 0., 1.) == []


def test_the_padded_hull_of_a_square_is_a_bigger_square():
    hull = convex_hull([(0., 0.), (1., 0.), (1., 1.), (0., 1.)], .5)
    assert len(hull) == 4
    xs = sorted(x for x, _ in hull)
    assert xs[0] == approx(-.5) and xs[-1] == approx(1.5)


def test_an_interior_point_is_dropped_from_the_hull():
    hull = convex_hull([(0., 0.), (2., 0.), (1., 2.), (1., .5)], 0.)
    assert (1., .5) not in hull
    assert len(hull) == 3


def test_two_points_still_give_a_usable_hull():
    """a diatomic has no polygon -- the padded hull of a segment is a rectangle, not an empty tuple"""
    hull = convex_hull([(0., 0.), (1., 0.)], .4)
    assert len(hull) >= 4
    assert in_polygon(.5, .3, hull)
    assert not in_polygon(.5, .9, hull)


def test_one_point_gives_a_square_around_it():
    hull = convex_hull([(3., 3.)], .5)
    assert in_polygon(3., 3., hull)
    assert not in_polygon(3., 4., hull)


def test_a_hull_makes_the_field_undefined_outside_it():
    plane = {1: (0., 0.), 2: (4., 0.)}
    hull = convex_hull(list(plane.values()), .5)
    clipped = ScalarField({1: 1., 2: 1.}, plane, sigma=1., cutoff=6., hull=hull)
    unclipped = ScalarField({1: 1., 2: 1.}, plane, sigma=1., cutoff=6.)
    assert clipped.at(2., 0.) == approx(unclipped.at(2., 0.))
    assert unclipped.at(2., 3.) is not None
    assert clipped.at(2., 3.) is None


def test_a_clipped_field_traces_bands_that_stop_at_the_hull():
    """the point of the whole exercise: no band runs off as a rectangle past the structure"""
    plane = {1: (0., 0.), 2: (2., 0.), 3: (1., 1.7)}
    hull = convex_hull(list(plane.values()), .4)
    field = ScalarField({1: 1., 2: .5, 3: 0.}, plane, sigma=1.2, cutoff=8., hull=hull)
    grid = sample(field, field.bounds(1.), .1)
    for line in isolines(grid, .4):
        for x, y in line:
            assert in_polygon(x, y, convex_hull(list(plane.values()), .55)), (x, y)


def test_sigma_and_cutoff_must_be_named():
    """positional (values, plane, sigma, cutoff) freezes an argument order nobody chose"""
    with raises(TypeError):
        ScalarField({1: 1.}, {1: (0., 0.)}, 1., 4.)


def test_the_field_is_defined_out_to_the_stated_cutoff_at_any_sigma():
    """the guard must mean "no atom in range", not "the weights got small" -- at sigma=.35 and
    cutoff=4 the old weight threshold ended the field at r=2.60, 65% of the number it documents"""
    field = ScalarField({1: 1.}, {1: (0., 0.)}, sigma=.35, cutoff=4.)
    # just inside the cutoff: must be defined
    assert field.at(3.99, 0.) is not None
    # just outside the cutoff: must be None
    assert field.at(4.01, 0.) is None


def test_none_corners_are_skipped_and_no_contour_appears_at_the_cutoff_boundary():
    """a cell with a None corner is skipped; substituting 0. invents a hard wall at the cutoff boundary

    THE LEVEL HAS TO BE BELOW f(cutoff) = exp(-2) = .1353: the .5 ring sits at r = 1.1774, far inside the
    boundary, where an invented wall cannot reach it.  At .05 there is honestly no contour at all, and a
    substituted 0. produces one ring at r ≈ 1.97..2.06 -- exactly the cutoff.
    """
    field = ScalarField({1: 1.}, {1: (0., 0.)}, sigma=1., cutoff=2.)
    # box extends well past the cutoff, so the grid has a ring of None cells at the boundary
    box = field.bounds(3.)
    grid = sample(field, box, .1)
    assert any(v is None for v in grid.z), 'sanity: grid must have None cells'

    # .5 > f(cutoff): a real feature, inside the boundary
    ring = isolines(grid, .5)
    assert len(ring) == 1 and ring[0][0] == ring[0][-1]
    assert max(hypot(x, y) for x, y in ring[0]) < field.cutoff

    # .05 < f(cutoff): the field never reaches it inside its own support, so any contour is the wall
    traced = isolines(grid, .05)
    assert traced == [], (
        f'level .05 is below f(cutoff) = {exp(-2.):.4f}, so no vertex of the field equals it; '
        f'{len(traced)} contour(s) were traced at radii '
        f'{[round(hypot(*v), 3) for line in traced for v in line][:6]} -- a None-corner skip is missing '
        f'and the cutoff boundary itself is being contoured')


def test_sigma_squared_is_used_not_sigma():
    """at sigma=1 the square is invisible, so this pins an exact value at sigma=2: one sigma out is always
    exp(-.5), and without the square it would be exp(-1.)"""
    field = ScalarField({1: 1.}, {1: (0., 0.)}, sigma=2., cutoff=10.)
    # at distance r=2 (one sigma), f = exp(-r^2 / 2*sigma^2) = exp(-4/8) = exp(-.5)
    # if sigma² were dropped: exp(-r^2 / 2*sigma) = exp(-4/4) = exp(-1.) ≠ exp(-.5)
    assert field.at(2., 0.) == approx(exp(-.5))


def test_the_gradient_carries_the_one_over_sigma_squared_factor():
    """same blind spot as the field value: at sigma=1 the factor is a division by one, at sigma=2 it is
    a factor of four"""
    field = ScalarField({1: 1.}, {1: (0., 0.)}, sigma=2., cutoff=10.)
    # ∂f/∂x = -x/σ² · exp(-r²/2σ²).  at x=2, r=2: -2/4 · exp(-.5) = -0.303265
    # without 1/σ² factor: -x · exp(-r²/2σ²) = -2 · exp(-.5) = -1.213061
    gx, gy = field.gradient(2., 0.)
    assert gx == approx(-2. / 4. * exp(-.5))
    assert gy == approx(0.)
