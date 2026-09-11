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
"""A colormap is a pure function from a float to an RGB triple, and the domain is the interesting half.

The domain rule: a diverging map fitted to data spanning zero comes out SYMMETRIC.  Fitting `coolwarm` to
charges from -0.3 to +0.8 unsymmetrised puts white at +0.25, so every neutral carbon draws pale blue and
the molecule reads as polarised where it is not.  The map's midpoint means "zero" or it means nothing.
"""
from pytest import approx, raises
from chython.depict.colormap import NAMED_COLORMAPS, Colormap, as_colormap


def test_a_value_at_a_stop_is_that_stops_colour():
    cmap = Colormap(((0., (0., 0., 1.)), (1., (1., 0., 0.))), vmin=0., vmax=1.)
    assert cmap.at(0.) == approx((0., 0., 1.))
    assert cmap.at(1.) == approx((1., 0., 0.))


def test_a_value_between_stops_is_interpolated_channelwise():
    cmap = Colormap(((0., (0., 0., 0.)), (1., (1., 1., 1.))), vmin=0., vmax=1.)
    assert cmap.at(.25) == approx((.25, .25, .25))


def test_interpolation_finds_the_right_span_in_a_multi_stop_map():
    cmap = Colormap(((0., (0., 0., 0.)), (.5, (1., 0., 0.)), (1., (1., 1., 0.))), vmin=0., vmax=1.)
    assert cmap.at(.75) == approx((1., .5, 0.))


def test_a_value_outside_the_domain_is_clamped_not_extrapolated():
    """extrapolating channels leaves the unit cube and produces a colour no format can write"""
    cmap = Colormap(((0., (0., 0., 1.)), (1., (1., 0., 0.))), vmin=0., vmax=1.)
    assert cmap.at(-5.) == approx((0., 0., 1.))
    assert cmap.at(5.) == approx((1., 0., 0.))


def test_normalised_is_the_one_value_to_fraction_map_and_it_clamps():
    """`(v - vmin) / (vmax - vmin)` clamped to [0, 1] is also `AtomHalo`'s radius and `BondScale`'s width:
    a halo whose colour says one thing and whose radius says another cannot be read, so one definition"""
    cmap = Colormap(((0., (0., 0., 1.)), (1., (1., 0., 0.))), vmin=-1., vmax=3.)
    assert cmap.normalised(-1.) == approx(0.)
    assert cmap.normalised(1.) == approx(.5)
    assert cmap.normalised(3.) == approx(1.)
    assert cmap.normalised(-99.) == 0. and cmap.normalised(99.) == 1., 'clamped, not extrapolated'


def test_hex_at_is_the_svg_spelling():
    cmap = Colormap(((0., (0., 0., 0.)), (1., (1., 1., 1.))), vmin=0., vmax=1.)
    assert cmap.hex_at(1.) == '#ffffff'
    assert cmap.hex_at(0.) == '#000000'


def test_a_diverging_map_is_symmetric_by_default():
    fitted = Colormap.named('coolwarm').fitted([-.3, .1, .8])
    assert fitted.vmin == approx(-.8)
    assert fitted.vmax == approx(.8)
    mid = fitted.at(0.)
    assert mid == approx(fitted.stops[len(fitted.stops) // 2][1]), 'zero lands on the middle stop'
    assert mid[0] == approx(mid[2], abs=.1), 'and the middle stop of a diverging map is neutral'


def test_a_diverging_map_fitted_to_one_sided_data_keeps_its_own_range():
    fitted = Colormap.named('coolwarm').fitted([.2, .5, .9])
    assert fitted.vmin == approx(.2)
    assert fitted.vmax == approx(.9)


def test_a_sequential_map_is_never_symmetrised():
    fitted = Colormap.named('viridis').fitted([-.3, .8])
    assert fitted.vmin == approx(-.3)
    assert fitted.vmax == approx(.8)


def test_fitting_constant_data_gives_a_domain_of_nonzero_width():
    """every atom the same value must not divide by zero; it draws one flat colour"""
    base = Colormap.named('viridis')
    fitted = base.fitted([.5, .5, .5])
    assert fitted.vmax > fitted.vmin
    colour = fitted.at(.5)
    assert all(0. <= c <= 1. for c in colour), 'channels must be in unit range'
    # padding is symmetric around the constant value, so 0.5 maps to position 0.5 in the stops
    assert colour == approx(base.at(.5)), 'constant-data padding keeps the value at the map midpoint'


def test_fitting_an_empty_mapping_is_refused():
    with raises(ValueError, match='no values'):
        Colormap.named('viridis').fitted([])


def test_explicit_bounds_survive_fitting():
    cmap = Colormap.named('coolwarm').fitted([-.3, .8], vmin=-1., vmax=1.)
    assert (cmap.vmin, cmap.vmax) == approx((-1., 1.))


def test_the_named_maps_are_all_present_and_well_formed():
    for name, cmap in NAMED_COLORMAPS.items():
        assert cmap.stops[0][0] == approx(0.), name
        assert cmap.stops[-1][0] == approx(1.), name
        assert all(a[0] < b[0] for a, b in zip(cmap.stops, cmap.stops[1:])), name
        assert all(0. <= c <= 1. for _, colour in cmap.stops for c in colour), name


def test_an_unknown_name_is_refused_and_lists_what_exists():
    with raises(KeyError, match='viridis'):
        Colormap.named('inferno')


def test_stops_out_of_order_are_refused():
    with raises(ValueError, match='ascending'):
        Colormap(((1., (0., 0., 0.)), (0., (1., 1., 1.))), vmin=0., vmax=1.)


def test_a_colormap_is_hashable_and_immutable():
    cmap = Colormap.named('viridis')
    assert hash(cmap) == hash(Colormap.named('viridis'))
    with raises(AttributeError):
        cmap.vmin = 3.


def test_coolwarm_is_marked_diverging_and_viridis_is_not():
    assert Colormap.named('coolwarm').diverging
    assert not Colormap.named('viridis').diverging


def test_the_six_named_maps_are_exactly_the_ones_the_spec_lists():
    assert set(NAMED_COLORMAPS) == {'viridis', 'cividis', 'coolwarm', 'RdBu', 'PiYG', 'mono'}
    assert {n for n, c in NAMED_COLORMAPS.items() if c.diverging} == {'coolwarm', 'RdBu', 'PiYG'}


def test_mono_is_grey_end_to_end():
    """the greyscale-journal map: every stop has three equal channels, so it survives a mono print"""
    for _, (r, g, b) in NAMED_COLORMAPS['mono'].stops:
        assert r == approx(g) == approx(b)


def test_a_caller_may_pass_a_list_of_stops():
    cmap = as_colormap([(0., (0., 0., 0.)), (1., (1., 0., 0.))])
    assert cmap.at(cmap.vmax) == approx((1., 0., 0.))


def test_a_bare_colour_list_is_spread_evenly():
    cmap = as_colormap(['#000000', '#808080', '#ffffff'])
    assert [p for p, _ in cmap.stops] == approx([0., .5, 1.])


def test_a_caller_may_pass_a_callable():
    cmap = as_colormap(lambda t: (t, t, t))
    assert len(cmap.stops) == 17
    assert cmap.at(cmap.vmin) == approx((0., 0., 0.))
    assert cmap.at(cmap.vmax) == approx((1., 1., 1.))


def test_a_colormap_passes_through_as_colormap_unchanged():
    cmap = Colormap.named('PiYG')
    assert as_colormap(cmap) is cmap


def test_a_callable_that_does_not_return_a_triple_is_refused_where_it_was_written():
    with raises(ValueError, match='three'):
        as_colormap(lambda t: t)


def test_hex_to_float_preserves_channel_order():
    """a transposed r/b channel in _hex_to_float would silently mirror every user-supplied hex list

    #ff0080 has r=255, g=0, b=128, so a swap returns (128/255, 0, 1) instead of (1, 0, 128/255).
    """
    from chython.depict.colormap import _hex_to_float
    r, g, b = _hex_to_float('#ff0080')
    assert r == approx(1.)
    assert g == approx(0.)
    assert b == approx(128 / 255.)

    # also verify through the public coercion path: red at 0, blue at 1
    cmap = as_colormap(['#ff0000', '#0000ff'])
    assert cmap.hex_at(cmap.vmin) == '#ff0000'


def test_hex_at_rounds_not_truncates():
    """hex_at must use round(), not int() truncation

    At the midpoint of a black-to-white ramp the channel is exactly 0.5, and 0.5 * 255 = 127.5:
    round gives 128 (#808080), int gives 127 (#7f7f7f).
    """
    cmap = Colormap(((0., (0., 0., 0.)), (1., (1., 1., 1.))), vmin=0., vmax=1.)
    assert cmap.hex_at(0.5) == '#808080'
