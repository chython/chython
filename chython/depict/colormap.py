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
"""Colormap: a pure function from a float value to an RGB triple.

The named maps are 9-stop tables from public-domain or attributed sources; each carries its
attribution beside it below.  A diverging map fitted to data that spans zero has its domain
symmetrised to ±max(|vmin|, vmax), so white stays on zero and neutral atoms do not look coloured.
"""
from bisect import bisect_right
from collections.abc import Sequence
from dataclasses import dataclass, field, replace

from .scene import rgb as _scene_rgb, to_hex as _scene_to_hex


__all__ = ['Colormap', 'NAMED_COLORMAPS', 'as_colormap']


def _hex_to_float(colour: str) -> tuple[float, float, float]:
    """Parse a normalised '#rrggbb' hex string to a float (r, g, b) triple in [0, 1]."""
    h = _scene_to_hex(colour)  # normalise / validate
    r = int(h[1:3], 16) / 255.
    g = int(h[3:5], 16) / 255.
    b = int(h[5:7], 16) / 255.
    return r, g, b


@dataclass(frozen=True, slots=True)
class Colormap:
    """A piecewise-linear colormap: a tuple of (position, (r, g, b)) stops in ascending position order.

    Positions are in [0, 1]; channel values are in [0, 1].  The domain [vmin, vmax] maps linearly
    onto [0, 1] before the stops are consulted.  ``diverging`` is set by three of the six named maps and
    defaults to False for a hand-built one.
    """
    stops: tuple[tuple[float, tuple[float, float, float]], ...]
    vmin: float = 0.
    vmax: float = 1.
    diverging: bool = False
    # precomputed stop positions for the bisect in `at`.  Excluded from __init__/__repr__/__hash__/__eq__
    # because it is derived from `stops`; including it would make identical maps hash differently.
    _positions: tuple[float, ...] = field(init=False, repr=False, hash=False, compare=False,
                                          default=())

    def __post_init__(self):
        if self.vmax <= self.vmin:
            raise ValueError(f'vmax ({self.vmax}) must be greater than vmin ({self.vmin})')
        if len(self.stops) < 2:
            raise ValueError(f'a colormap needs at least two stops, got {len(self.stops)}')
        positions = tuple(p for p, _ in self.stops)
        if any(a >= b for a, b in zip(positions, positions[1:])):
            raise ValueError('stop positions must be strictly ascending')
        if abs(positions[0]) > 1e-9:
            raise ValueError(f'first stop must be at position 0, got {positions[0]}')
        if abs(positions[-1] - 1.) > 1e-9:
            raise ValueError(f'last stop must be at position 1, got {positions[-1]}')
        for _, (r, g, b) in self.stops:
            if not (0. <= r <= 1. and 0. <= g <= 1. and 0. <= b <= 1.):
                raise ValueError(f'channel values must be in [0, 1], got ({r}, {g}, {b})')
        object.__setattr__(self, '_positions', positions)

    def normalised(self, value: float) -> float:
        """``value``'s position within [vmin, vmax] as a fraction in [0, 1], clamped at both ends.

        One definition, because three consumers must agree: ``at`` picks a colour with it, ``AtomHalo`` a
        radius, ``BondScale`` a width.  A halo whose colour and radius disagree cannot be read.
        """
        t = (value - self.vmin) / (self.vmax - self.vmin)
        return max(0., min(1., t))

    def at(self, value: float) -> tuple[float, float, float]:
        """Return the interpolated RGB triple for ``value``, clamped to [vmin, vmax]."""
        t = self.normalised(value)

        stops = self.stops
        if t <= stops[0][0]:
            return stops[0][1]
        if t >= stops[-1][0]:
            return stops[-1][1]

        # the first stop past t, so the span containing t is [i-1, i]
        i = bisect_right(self._positions, t)
        p0, (r0, g0, b0) = stops[i - 1]
        p1, (r1, g1, b1) = stops[i]
        f = (t - p0) / (p1 - p0)
        return (r0 + f * (r1 - r0), g0 + f * (g1 - g0), b0 + f * (b1 - b0))

    def hex_at(self, value: float) -> str:
        """Return the hex colour string for ``value``."""
        r, g, b = self.at(value)
        return _scene_rgb(round(r * 255), round(g * 255), round(b * 255))

    def fitted(self, values: Sequence[float], *,
               vmin: float | None = None,
               vmax: float | None = None) -> 'Colormap':
        """Return a copy with a domain derived from ``values``.

        If ``vmin``/``vmax`` are given explicitly they are used unchanged.  Otherwise the domain is
        derived from the data:

        - An empty sequence is refused.
        - A zero-width range (all values equal) is padded by ±0.5 (or ±5 % of the value when
          non-zero), so that ``at`` never divides by zero.
        - When ``self.diverging`` and the data spans zero (raw_min < 0 < raw_max), the domain is
          symmetrised: both ends are set to max(|raw_min|, raw_max).
        """
        values = list(values)
        if not values:
            raise ValueError('no values to fit a colormap to')

        if vmin is not None and vmax is not None:
            return replace(self, vmin=float(vmin), vmax=float(vmax))

        raw_min = min(values)
        raw_max = max(values)

        if raw_min == raw_max:
            v = raw_min
            half = abs(v) * 0.05 if v != 0. else 0.5
            new_min = v - half
            new_max = v + half
        elif self.diverging and raw_min < 0. < raw_max:
            extent = max(abs(raw_min), raw_max)
            new_min = -extent
            new_max = extent
        else:
            new_min = raw_min
            new_max = raw_max

        if vmin is not None:
            new_min = float(vmin)
        if vmax is not None:
            new_max = float(vmax)

        return replace(self, vmin=new_min, vmax=new_max)

    @classmethod
    def named(cls, name: str) -> 'Colormap':
        """Return one of the six named colormaps by name.

        Raises ``KeyError`` listing the available names when the name is unknown.
        """
        try:
            return NAMED_COLORMAPS[name]
        except KeyError:
            available = ', '.join(sorted(NAMED_COLORMAPS))
            raise KeyError(f'{name!r} is not a named colormap; available: {available}') from None


# viridis: matplotlib (CC0).  Perceptually uniform sequential, the default for one-sided data.
_VIRIDIS = Colormap(
    stops=(
        (0.000, (0.267, 0.005, 0.329)),
        (0.125, (0.283, 0.141, 0.458)),
        (0.250, (0.254, 0.265, 0.530)),
        (0.375, (0.207, 0.372, 0.553)),
        (0.500, (0.163, 0.471, 0.558)),
        (0.625, (0.128, 0.567, 0.551)),
        (0.750, (0.197, 0.659, 0.498)),
        (0.875, (0.477, 0.757, 0.345)),
        (1.000, (0.993, 0.906, 0.144)),
    ),
    vmin=0., vmax=1., diverging=False,
)

# cividis: Nuñez, Anderton & Renslow 2018 (CC0 via matplotlib).
# Designed for colour-vision deficiency — reads identically to deuteranopes.
_CIVIDIS = Colormap(
    stops=(
        (0.000, (0.000, 0.135, 0.304)),
        (0.125, (0.093, 0.191, 0.384)),
        (0.250, (0.185, 0.248, 0.427)),
        (0.375, (0.279, 0.309, 0.432)),
        (0.500, (0.377, 0.374, 0.418)),
        (0.625, (0.483, 0.443, 0.390)),
        (0.750, (0.602, 0.519, 0.346)),
        (0.875, (0.737, 0.604, 0.278)),
        (1.000, (0.993, 0.906, 0.144)),
    ),
    vmin=0., vmax=1., diverging=False,
)

# coolwarm: Moreland, K. "Diverging Color Maps for Scientific Visualization" (public domain).
# White is zero, blue is negative, red is positive.
_COOLWARM = Colormap(
    stops=(
        (0.000, (0.230, 0.299, 0.754)),
        (0.125, (0.350, 0.461, 0.858)),
        (0.250, (0.487, 0.607, 0.921)),
        (0.375, (0.628, 0.736, 0.941)),
        (0.500, (0.865, 0.865, 0.865)),
        (0.625, (0.952, 0.704, 0.595)),
        (0.750, (0.918, 0.516, 0.404)),
        (0.875, (0.804, 0.311, 0.288)),
        (1.000, (0.706, 0.016, 0.150)),
    ),
    vmin=0., vmax=1., diverging=True,
)

# RdBu: ColorBrewer 2.0, Cynthia Brewer (http://colorbrewer2.org — free for use with attribution).
# The journal-conventional charge map.
_RDBU = Colormap(
    stops=(
        (0.000, (0.647, 0.000, 0.149)),
        (0.125, (0.839, 0.188, 0.153)),
        (0.250, (0.957, 0.478, 0.357)),
        (0.375, (0.992, 0.733, 0.635)),
        (0.500, (0.969, 0.969, 0.969)),
        (0.625, (0.643, 0.812, 0.894)),
        (0.750, (0.353, 0.627, 0.804)),
        (0.875, (0.137, 0.404, 0.675)),
        (1.000, (0.020, 0.188, 0.380)),
    ),
    vmin=0., vmax=1., diverging=True,
)

# PiYG: ColorBrewer 2.0, Cynthia Brewer (same licence, same attribution as RdBu above).
# A diverging pair that stays distinguishable beside a blue highlight.
_PIYG = Colormap(
    stops=(
        (0.000, (0.557, 0.004, 0.322)),
        (0.125, (0.773, 0.106, 0.490)),
        (0.250, (0.871, 0.467, 0.682)),
        (0.375, (0.945, 0.714, 0.855)),
        (0.500, (0.969, 0.969, 0.969)),
        (0.625, (0.820, 0.902, 0.627)),
        (0.750, (0.576, 0.769, 0.349)),
        (0.875, (0.302, 0.573, 0.129)),
        (1.000, (0.153, 0.392, 0.098)),
    ),
    vmin=0., vmax=1., diverging=True,
)

# mono: a three-stop black-to-white ramp for greyscale journals.
# Every stop has three equal channels, so it survives monochrome print.
_MONO = Colormap(
    stops=(
        (0.0, (0.0, 0.0, 0.0)),
        (0.5, (0.5, 0.5, 0.5)),
        (1.0, (1.0, 1.0, 1.0)),
    ),
    vmin=0., vmax=1., diverging=False,
)

NAMED_COLORMAPS: dict[str, Colormap] = {
    'viridis': _VIRIDIS,
    'cividis': _CIVIDIS,
    'coolwarm': _COOLWARM,
    'RdBu': _RDBU,
    'PiYG': _PIYG,
    'mono': _MONO,
}


def as_colormap(spec) -> Colormap:
    """Coerce whatever the caller passed into a ``Colormap``.

    Four accepted shapes: the map itself, its name, a list of stops, or a function.  A callable is
    sampled at 17 positions rather than stored -- finer than a printer resolves -- so every consumer sees
    one type, ``at`` stays a bisect, and no backend calls user code at draw time.
    """
    if isinstance(spec, Colormap):
        return spec
    if isinstance(spec, str):
        return Colormap.named(spec)
    if callable(spec):
        stops = []
        for i in range(17):
            p = i / 16.
            colour = spec(p)
            try:
                r, g, b = colour
            except (TypeError, ValueError):
                raise ValueError(f'a colormap callable must return three channels, got {colour!r}')
            stops.append((p, (float(r), float(g), float(b))))
        return Colormap(tuple(stops), vmin=0., vmax=1.)

    spec = list(spec)
    if len(spec) < 2:
        raise ValueError(f'a colormap needs at least two stops, got {len(spec)}')
    if isinstance(spec[0], str):  # bare colours: spread them evenly
        n = len(spec) - 1
        stops = tuple((i / n, _hex_to_float(c)) for i, c in enumerate(spec))
    else:
        stops = tuple((float(p), (float(c[0]), float(c[1]), float(c[2]))) for p, c in spec)
    return Colormap(stops, vmin=0., vmax=1.)
