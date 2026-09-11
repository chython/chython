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
"""Five overlay kinds that turn QM/scalar data into scene nodes.

Two rules: an overlay returns (under, over) and never sorts -- fields and halos go UNDER the structure so
bonds and labels stay readable, value labels OVER it.  And overlapping highlights are ONE group with group
opacity, never a boolean union: two 50%-opaque discs drawn apart make a 75%-opaque lens where they meet.
"""
from collections.abc import Sequence
from dataclasses import dataclass, field
from math import fsum, hypot, sqrt
from typing import NamedTuple

from ..core import LogRecord
from .colormap import Colormap, as_colormap
from .field import (ScalarField, contour_levels, convex_hull,
                    isolines, refine, sample, to_cubics, trim_asymptote)
from .label import MINUS
from .metrics import text_box
from .scene import Box, Group, Path, Text, TextRun, circle, close, line, move, polyline, rgb, to_hex
from .style import DepictStyle


__all__ = ['Highlight', 'AtomHalo', 'AtomField', 'BondScale', 'ValueLabels',
           'Overlay', 'Swatch', 'bands_of', 'render_overlays', 'scale_of', 'tiled_swatches']


# Each number below bounds an ALGORITHM rather than describing the page, which is why none is a
# `DepictStyle` field: a style field is for a number a caller can want differently.

#: Pitch between consecutive `Highlight(style='outline')` rings, in outline widths.  Two widths of clear
#: space, so N overlapping outlines read as N concentric rings; at 1.0 the strokes touch and read as one.
_OUTLINE_RING_PITCH = 2.5

#: A `Highlight.label`'s size, and its clearance above the group's box, as fractions of `LabelStyle.size`.
#: Three quarters is the smallest that stays legible at 83 mm beside a full-size atom symbol; the gap
#: clears the group's stroke without floating free of it.
_HIGHLIGHT_LABEL_SCALE = .75
_HIGHLIGHT_LABEL_GAP = .2

#: How much darker a band's boundary stroke is than the band it bounds.  Dark enough to read as a line
#: over its own fill at 83 mm, light enough not to read as a second colour in the sequence.
_BAND_OUTLINE_DARKEN = .65

#: A value label's clearance from the atom's label box, and a bond value label's perpendicular offset
#: from the bond axis, both in label heights.  A little over one height clears the bond's own stroke.
_VALUE_LABEL_GAP = .8
_BOND_LABEL_OFFSET = 1.1

#: `_nudge_clear`'s step, in label heights, and its iteration limit.  20 × .15 is three label heights:
#: past that the label is no longer beside the atom it names, so the loop must end rather than run away.
_NUDGE_STEP = .15
_NUDGE_LIMIT = 20


def _check_atoms(mol, atom_ids):
    """Raise ValueError if any id is not in mol."""
    mol_ids = {a.n for a in mol.atoms()}
    for sid in atom_ids:
        if sid not in mol_ids:
            raise ValueError(f'atom id {sid} is not in this molecule')


def _check_bonds(mol, bond_pairs):
    """Raise ValueError if any pair is not a bond in mol."""
    for n, m in bond_pairs:
        if mol.order_of(n, m) is None:
            raise ValueError(f'atoms {n} and {m} are not bonded in this molecule')


def _mean_bond_length(mol, plane) -> float:
    """Mean Euclidean bond length in `plane`.  Falls back to 1.0 for a single atom.

    `fsum` and not `sum`: this length sets a field's grid spacing, so its last bit is a contour's
    position and a written coordinate.  `sum` accumulates floats in extended precision from 3.12 and
    naively before it, which put a grid line of `docs/images/field-node.svg` on either side of the
    fourth decimal by interpreter version.  A correctly rounded total is the same on all of them.
    """
    lengths = []
    for bond in mol.bonds():
        nx, ny = plane[bond.n]
        mx, my = plane[bond.m]
        lengths.append(hypot(nx - mx, ny - my))
    if not lengths:
        return 1.0
    return fsum(lengths) / len(lengths)


def _free_direction(mol, plane, n) -> tuple[float, float]:
    """Unit vector pointing away from the centroid of atom n's neighbours.

    Places a value label clear of the atom symbol.  A degree-0 atom answers (0, 1), above the atom.
    """
    neighbours = list(mol.neighbors_of(n))
    if not neighbours:
        return 0., 1.
    ax, ay = plane[n]
    cx = fsum(plane[nb][0] for nb in neighbours) / len(neighbours)
    cy = fsum(plane[nb][1] for nb in neighbours) / len(neighbours)
    dx = ax - cx
    dy = ay - cy
    length = hypot(dx, dy)
    if length < 1e-9:
        return 0., 1.
    return dx / length, dy / length


def _darker(hex_colour: str, factor: float) -> str:
    """Return a darker shade of `hex_colour` by multiplying channels by `factor`."""
    r = int(hex_colour[1:3], 16)
    g = int(hex_colour[3:5], 16)
    b = int(hex_colour[5:7], 16)
    return rgb(round(r * factor), round(g * factor), round(b * factor))


def _fitted(overlay) -> Colormap:
    """The fitted colormap `overlay` draws with.  ONE definition, five callers.

    `AtomHalo.render`, `AtomField`'s band tracer, `BondScale.bond_widths`, `BondScale.bond_colours` and
    `scale_of` all need `domain is not None → fitted(vmin, vmax) else fitted()`, and `scale_of` is what
    the colorbar labels the figure by, so the five have to agree exactly.
    """
    cmap = as_colormap(overlay.colormap)
    if overlay.domain is not None:
        return cmap.fitted(overlay.values.values(), vmin=overlay.domain[0], vmax=overlay.domain[1])
    return cmap.fitted(overlay.values.values())


@dataclass(frozen=True, slots=True)
class Highlight:
    """A coloured halo behind one or more atoms and/or bonds.

    `style='fill'`    — filled disc at `style.highlight.opacity` on the group; children are plain.
    `style='outline'` — stroked ring, offset outward by cycle * outline_width * `_OUTLINE_RING_PITCH`
                        so two overlapping outlines read as concentric rings.
    `color=None`      — resolves to `style.highlight.palette[cycle % len(palette)]`.
    `label`           — optional string drawn once at the outside of the group's bounding box.
    """
    atoms: tuple = ()
    bonds: tuple = ()
    color: str | None = None
    style: str = 'fill'
    label: str | None = None

    def __post_init__(self):
        if self.style not in ('fill', 'outline'):
            raise ValueError(f"Highlight.style must be 'fill' or 'outline', got {self.style!r}")
        object.__setattr__(self, 'atoms', tuple(self.atoms))
        object.__setattr__(self, 'bonds', tuple(self.bonds))
        if self.color is not None:
            object.__setattr__(self, 'color', to_hex(self.color))

    def render(self, mol, plane, boxes, style: DepictStyle, *, cycle: int = 0,
               log=None) -> tuple[list, list]:
        _check_atoms(mol, self.atoms)
        _check_bonds(mol, self.bonds)

        colour = self.color if self.color is not None else \
            style.highlight.palette[cycle % len(style.highlight.palette)]

        radius = style.highlight.radius
        bond_width = style.highlight.bond_width

        children = []

        for sid in self.atoms:
            ax, ay = plane[sid]
            if self.style == 'fill':
                r = radius
            else:
                # outline: offset outward for each cycle position
                r = radius + cycle * style.highlight.outline_width * _OUTLINE_RING_PITCH
            path = Path([circle(ax, ay, r)], fill=colour if self.style == 'fill' else None,
                        stroke=colour if self.style == 'outline' else None,
                        width=style.highlight.outline_width if self.style == 'outline' else None)
            children.append(path)

        # Bond capsules: a round-capped stroke at width = bond_width is identical to a full capsule path,
        # since the round cap is exactly the half-disc the capsule ends in.
        for n, m in self.bonds:
            nx, ny = plane[n]
            mx, my = plane[m]
            if self.style == 'fill':
                w = bond_width
            else:
                # twice the radial pitch, because a stroke widens on BOTH sides of the bond axis
                w = bond_width + 2 * cycle * style.highlight.outline_width * _OUTLINE_RING_PITCH
            path = Path([(move(nx, ny), line(mx, my))],
                        fill=None,
                        stroke=colour,
                        width=w,
                        cap='round')
            children.append(path)

        if not children:
            return [], []

        if self.style == 'fill':
            group = Group(children, opacity=style.highlight.opacity)
        else:
            group = Group(children)

        under = [group]
        over = []

        # The label goes in the OVER list: text drawn under the bonds is occluded at print resolution.
        if self.label is not None:
            lbl_size = style.label.size * _HIGHLIGHT_LABEL_SCALE
            b = group.bounds
            lbl_text = Text(
                [TextRun(self.label, family=style.label.family, size=lbl_size)],
                x=(b.min_x + b.max_x) / 2.,
                y=b.max_y + lbl_size * _HIGHLIGHT_LABEL_GAP,
                anchor='middle',
                fill=colour,
            )
            over.append(lbl_text)

        return under, over


@dataclass(frozen=True, slots=True)
class AtomHalo:
    """Per-atom coloured disc, sized and/or coloured by a scalar value.

    `encode='color'`  — fixed radius, fill colour from the colormap.
    `encode='size'`   — radius mapped to [halo_min_radius, halo_max_radius], fill from style default.
    `encode='both'`   — radius AND colour both encoded.
    `radius=None`     — uses `style.highlight.radius` for colour-only, scaled range for size encoding.
    """
    values: dict = field(default_factory=dict)
    colormap: object = 'coolwarm'
    encode: str = 'color'
    radius: float | None = None
    domain: tuple[float, float] | None = None

    def __post_init__(self):
        if self.encode not in ('color', 'size', 'both'):
            raise ValueError(f"AtomHalo.encode must be 'color', 'size' or 'both', got {self.encode!r}")
        object.__setattr__(self, 'values', dict(self.values))

    def render(self, mol, plane, boxes, style: DepictStyle, *, cycle: int = 0,
               log=None) -> tuple[list, list]:
        if not self.values:
            return [], []

        cmap = _fitted(self)

        rmin = style.field.halo_min_radius
        rmax = style.field.halo_max_radius
        default_radius = self.radius if self.radius is not None else style.highlight.radius

        children = []
        for sid, val in self.values.items():
            ax, ay = plane[sid]

            if self.encode in ('size', 'both'):
                r = rmin + cmap.normalised(val) * (rmax - rmin)
            else:
                r = default_radius

            if self.encode in ('color', 'both'):
                fill_colour = cmap.hex_at(val)
            else:
                fill_colour = style.atom.default_colour

            path = Path([circle(ax, ay, r)], fill=fill_colour)
            children.append(path)

        group = Group(children, opacity=style.highlight.opacity)
        return [group], []


@dataclass(frozen=True, slots=True)
class AtomField:
    """Smooth scalar field interpolated from per-atom values, rendered as filled contour bands.

    `sigma=None`   — `style.field.contour.sigma × mean bond length` (σ is a SCALE, not a distance).
    `levels`       — int (number of evenly-spaced bands) or an explicit sequence of level values.
    `clip=None`    — no clip, and the only setting that draws the whole field.  'hull' and 'box', both
                     padded by `contour.pad`, are clip paths drawn through the ATOMS, so they slice every
                     band reaching past them and the contours end in mid-air.
    `fill`/`isolines` — paint filled bands / stroke the boundaries.  Independent; both False draws nothing.
    `opacity=None` — `style.field.contour.fill_opacity`, on the group, since a band and its boundary
                     composite together.  Tracing lives in `_field_bands`, shared with the colorbar.
    """
    values: dict = field(default_factory=dict)
    colormap: object = 'coolwarm'
    levels: object = 9
    sigma: float | None = None
    domain: tuple[float, float] | None = None
    fill: bool = True
    isolines: bool = True
    isoline_labels: bool = False
    clip: str | None = None
    opacity: float | None = None

    def __post_init__(self):
        if self.clip not in ('hull', 'box', None):
            raise ValueError(f"AtomField.clip must be 'hull', 'box' or None, got {self.clip!r}")
        object.__setattr__(self, 'values', dict(self.values))

    def render(self, mol, plane, boxes, style: DepictStyle, *, cycle: int = 0,
               log=None) -> tuple[list, list]:
        if not self.values or (not self.fill and not self.isolines):
            return [], []

        bands = _field_bands(self, mol, plane, style, log=log)
        if not bands:
            return [], []

        # Two accumulation lists so fill+outline pairs always precede open-isoline strokes: children[0] is
        # a fill and children[1] its outline.
        fill_children: list = []
        open_children: list = []
        for swatch, closed_subs, open_subs in bands:
            colour = swatch.colour
            if self.fill:
                if closed_subs:
                    # Two sibling Paths per band, the fill then its boundary stroke: a printer that drops
                    # one of the two still renders the other.
                    fill_children.append(Path(closed_subs, fill=colour))
                    if self.isolines:
                        fill_children.append(
                            Path(closed_subs, stroke=_darker(colour, _BAND_OUTLINE_DARKEN),
                                 width=style.field.contour.line_width))
                # An open contour encloses no region, so it is stroked at the band colour, not filled.
                for segs in open_subs:
                    open_children.append(Path([segs], stroke=colour,
                                              width=style.field.contour.line_width))
            else:
                # Line contour mode: every contour stroked at the band colour, closed or not.
                fill_children.append(Path(closed_subs + open_subs, stroke=colour,
                                          width=style.field.contour.line_width))

        children = fill_children + open_children
        if not children:
            return [], []

        # The field's own boundary wins by default: `ScalarField.at` answers None past `contour.cutoff`, so
        # a band closes where the Gaussians have decayed.  A hull or box is drawn through the ATOMS and cuts
        # the rings reaching past it, hence `clip=None`.  Either stays a renderer clip PATH on the group —
        # handing the hull to `ScalarField` would cut the field itself, and even un-cut bands would then
        # trace along the hull.
        if self.clip == 'hull':
            hull = convex_hull([plane[a] for a in self.values], style.field.contour.pad)
        elif self.clip == 'box':
            b = Box.of(lbl.box for lbl in boxes.values()).inflate(style.field.contour.pad)
            hull = ((b.min_x, b.min_y), (b.max_x, b.min_y),
                    (b.max_x, b.max_y), (b.min_x, b.max_y))
        else:
            hull = None

        clip_path = None
        if hull is not None and len(hull) >= 3:
            hull_segs = [move(*hull[0])] + [line(*p) for p in hull[1:]] + [close()]
            clip_path = Path([tuple(hull_segs)], fill='#000000')

        opacity = self.opacity if self.opacity is not None else style.field.contour.fill_opacity
        return [Group(children, opacity=opacity, clip=clip_path)], []


class Swatch(NamedTuple):
    """One band as the LEGEND sees it, and the whole of what `colorbar` is allowed to know.

    `level` and `colour` are the band's own.  `lo`/`hi` are the VALUE INTERVAL over which that colour is
    what a reader sees: a band is painted over the ones outside it, so its colour stands for the values
    between its own level and the next drawn band's level outward, and the outermost band runs to the
    domain end.  A level the field never reached leaves a GAP, which tells a reader where the bands are.
    `filled=False` says this level drew no closed contour, so no interval of values wears that colour and
    the bar draws a rule at `level` instead of a block; `lo == hi == level` there.
    """
    level: float
    colour: str
    lo: float
    hi: float
    filled: bool


def tiled_swatches(cmap: Colormap, levels: Sequence[float]) -> list[Swatch]:
    """Swatches for a scale that has NO bands to line up with: the domain tiled edge to edge.

    Each level owns the values nearer to it than to either neighbour, so the strip is continuous.  Right
    for an `AtomHalo`, which encodes its value as a RADIUS, and wrong for an `AtomField`, whose gaps must
    stay visible — `_field_bands` builds those.
    """
    ordered = sorted(levels)
    out = []
    for i, level in enumerate(ordered):
        lo = cmap.vmin if i == 0 else (ordered[i - 1] + level) / 2.
        hi = cmap.vmax if i == len(ordered) - 1 else (level + ordered[i + 1]) / 2.
        out.append(Swatch(level, cmap.hex_at(level), lo, hi, True))
    return out


def _field_bands(overlay: AtomField, mol, plane, style: DepictStyle, *,
                 log=None) -> list[tuple[Swatch, list[tuple], list[tuple]]]:
    """Every level of `overlay` that traces any contour at all, in painter's order.

    Returns `(swatch, closed_subpaths, open_subpaths)` per level; a level that traces nothing is absent.
    ONE definition, shared with `bands_of` so the colorbar labels what was drawn.  A closed chain bounds a
    region and can be filled, an open one is stroked; a chain is open where it meets `contour.cutoff`.

    An int `levels` COUNTS BANDS, so the levels come from the range the field was SAMPLED over and not from
    the colormap's domain: a level outside the field's own extremes has no cell with corners either side of
    it and draws nothing.  The sampled range is intersected with the domain, whose ends clamp the colour.
    """
    # σ is a SCALE: multiply by mean bond length so the field tracks the plane's units.  A clean2d() layout
    # has bond ≈ 0.825 mol units; at Å scale (bond ≈ 1.4) the same σ_scale covers the same bond count.
    sigma_scale = overlay.sigma if overlay.sigma is not None else style.field.contour.sigma
    mean_bl = _mean_bond_length(mol, plane)
    effective_sigma = sigma_scale * mean_bl
    cutoff = style.field.contour.cutoff

    # cutoff / 3. is a derivation: a Gaussian carries 99.7% of its mass within 3σ, so a cutoff inside 3σ
    # truncates the kernel while it still carries weight and every atom draws as a hard-edged disc.  The
    # usual cause is a plane in another unit (picometres, say).
    if effective_sigma > cutoff / 3. and log is not None:
        log.append(LogRecord(
            'depict:field-sigma',
            (),
            (f'effective sigma ({effective_sigma:.4g}) exceeds cutoff/3 '
             f'({cutoff / 3.:.4g}); the field is truncated inside its own sigma '
             f'and will draw as hard-edged discs.  '
             f'mean bond length measured: {mean_bl:.4g}, cutoff: {cutoff}.  '
             f'Check whether plane coordinates are in the expected units.'),
        ))

    sf = ScalarField(overlay.values, plane, sigma=effective_sigma, cutoff=cutoff, hull=None)
    cmap = _fitted(overlay)

    # Cutoff as padding, so every closed isoline fits inside the box and none is cut by the grid edge.
    grid = sample(sf, sf.bounds(cutoff), style.field.contour.grid)
    spacing = style.field.contour.grid
    refine_steps = style.field.contour.refine

    # BEFORE the levels, because the levels are of this range (see the docstring).
    sampled = [z for z in grid.z if z is not None]
    if not sampled:
        return []
    if isinstance(overlay.levels, int):
        lo, hi = max(min(sampled), cmap.vmin), min(max(sampled), cmap.vmax)
        if lo >= hi:
            lo, hi = min(sampled), max(sampled)
        level_list = contour_levels(overlay.levels, lo, hi)
    else:
        level_list = list(overlay.levels)
    if not level_list:
        return []

    # Painter's order is largest area first: for a diverging map that is the level closest to the midpoint,
    # for a sequential one the lowest level.
    midpoint = (cmap.vmin + cmap.vmax) / 2.
    if cmap.diverging:
        sorted_levels = sorted(level_list, key=lambda lv: abs(lv - midpoint))
    else:
        sorted_levels = sorted(level_list)

    traced: list[tuple[float, list[tuple], list[tuple]]] = []
    decisions: list[str] = []
    noteworthy = False
    for level in sorted_levels:
        closed_subs: list[tuple] = []
        open_subs: list[tuple] = []
        trimmed = 0
        for poly in isolines(grid, level):
            whole = len(poly)
            was_closed = whole >= 2 and poly[0] == poly[-1]
            runs = trim_asymptote(sf, poly)
            trimmed += whole - sum(len(run) for run in runs)
            for run in runs:
                # A run that IS the whole chain keeps the chain's closure; one that lost vertices to the
                # asymptote is open however it started, since closing it would invent the dropped arc.
                is_closed = was_closed and len(run) == whole
                r = refine(sf, run, level, refine_steps, spacing=spacing)
                if style.field.contour.smooth:
                    segs = to_cubics(sf, r, level, closed=is_closed)
                else:
                    segs = tuple(polyline(r, closed=is_closed))
                if segs:
                    (closed_subs if is_closed else open_subs).append(tuple(segs))
        if closed_subs or open_subs:
            traced.append((level, closed_subs, open_subs))
            decisions.append(f'{level:+.4g}: {len(closed_subs)} filled band(s), '
                             f'{len(open_subs)} open contour(s) stroked' +
                             (f', {trimmed} vertex(es) dropped as asymptotic' if trimmed else ''))
        else:
            decisions.append(f'{level:+.4g}: nothing traced — the field never reaches this level')
        # The record carries EVERY level and fires whenever one came out other than a whole band.
        if trimmed or not closed_subs:
            noteworthy = True

    if log is not None and noteworthy:
        log.append(LogRecord(
            'depict:field-levels',
            tuple(sorted(overlay.values)),
            f'field levels asked for: {len(sorted_levels)}, drawn: {len(traced)}, '
            f'over a sampled range of {min(sampled):+.4g}…{max(sampled):+.4g} '
            f'and a colormap domain of {cmap.vmin:+.4g}…{cmap.vmax:+.4g}.  ' + '; '.join(decisions),
        ))

    # The value interval each band's colour covers: between its own level and that of the band drawn just
    # OUTWARD of it, `filled` only, because an outward level that traced nothing hides nothing and so does
    # not end the interval.  Outward is away from the field's far-field value — the midpoint for a diverging
    # (symmetrized) map, below the domain for a sequential one — the two directions painter's order uses.
    filled_levels = sorted(lv for lv, closed_subs, _ in traced if closed_subs and overlay.fill)
    bands: list[tuple[Swatch, list[tuple], list[tuple]]] = []
    for level, closed_subs, open_subs in traced:
        colour = cmap.hex_at(level)
        if level not in filled_levels:
            bands.append((Swatch(level, colour, level, level, False), closed_subs, open_subs))
        elif not cmap.diverging or level >= midpoint:
            outward = [lv for lv in filled_levels if lv > level]
            # max()/min() against the level itself: an explicit level outside the domain would otherwise
            # give the bar an inverted interval.
            bands.append((Swatch(level, colour, level, min(outward) if outward else max(cmap.vmax, level),
                                 True), closed_subs, open_subs))
        else:
            outward = [lv for lv in filled_levels if lv < level]
            bands.append((Swatch(level, colour, max(outward) if outward else min(cmap.vmin, level), level,
                                 True), closed_subs, open_subs))
    return bands


def bands_of(overlay, mol, plane, style: DepictStyle) -> list[Swatch]:
    """The `Swatch` of every band `overlay` actually draws, in painter's order.

    The colorbar is built from THIS and never from a level count recomputed beside it: the swatches are the
    bands, one for one, or the legend does not label the picture.  An overlay that draws no bands (every
    kind but `AtomField`, or an `AtomField` with no values or levels) answers the empty list, which
    `figure.py` reads as "nothing to line up with".  The traversal is repeated rather than cached, since
    `values` is a mutable dict and a cache keyed on `id()` would outlive the object it described.
    """
    if not isinstance(overlay, AtomField) or not overlay.values:
        return []
    return [swatch for swatch, _, _ in _field_bands(overlay, mol, plane, style)]


@dataclass(frozen=True, slots=True)
class BondScale:
    """Per-bond width and/or colour scaling.

    `render` returns ([], []) because the width and colour of a bond are properties of the bond's own path,
    not a second path drawn on top of it -- two stacked strokes at different widths make a visible outline.
    The caller (figure.py) reads `bond_widths` and `bond_colours` and passes them into `bond_paths`.
    `encode` selects which of the two is populated (`'width'`, `'color'`, `'both'`); `width_range=None`
    takes `(style.field.bond_min_width, style.field.bond_max_width)`.
    """
    values: dict = field(default_factory=dict)
    encode: str = 'width'
    width_range: tuple[float, float] | None = None
    colormap: object = 'viridis'
    domain: tuple[float, float] | None = None

    def __post_init__(self):
        if self.encode not in ('width', 'color', 'both'):
            raise ValueError(f"BondScale.encode must be 'width', 'color' or 'both', "
                             f"got {self.encode!r}")
        if self.width_range is not None and self.encode == 'color':
            raise ValueError("width_range given with encode='color': a parameter that cannot act "
                             "is a typo, not a preference.  Use encode='both' or remove width_range.")
        object.__setattr__(self, 'values', dict(self.values))

    def _normalised_values(self) -> dict[tuple[int, int], float]:
        """Return values keyed by low-first (min(n,m), max(n,m)) pairs."""
        return {(min(n, m), max(n, m)): v for (n, m), v in self.values.items()}

    def bond_widths(self, mol, style: DepictStyle) -> dict[tuple[int, int], float]:
        """Per-bond width overrides keyed by the low-first pair.  Only named bonds are present.

        Always computed; `encode` determines which mapping render() and figure.py apply.
        """
        norm = self._normalised_values()
        if not norm:
            return {}
        cmap = _fitted(self)
        wmin = self.width_range[0] if self.width_range is not None else style.field.bond_min_width
        wmax = self.width_range[1] if self.width_range is not None else style.field.bond_max_width
        return {pair: wmin + cmap.normalised(val) * (wmax - wmin) for pair, val in norm.items()}

    def bond_colours(self, mol, style: DepictStyle) -> dict[tuple[int, int], str]:
        """Per-bond colour overrides keyed by the low-first pair.  Only named bonds are present.

        Always computed; `encode` determines which mapping render() and figure.py apply.
        """
        norm = self._normalised_values()
        if not norm:
            return {}
        cmap = _fitted(self)
        return {pair: cmap.hex_at(val) for pair, val in norm.items()}

    def render(self, mol, plane, boxes, style: DepictStyle, *, cycle: int = 0,
               log=None) -> tuple[list, list]:
        # Width and colour live on the bond's own path (see class docstring).
        return [], []


@dataclass(frozen=True, slots=True)
class ValueLabels:
    """Numeric value labels placed beside each named atom or bond midpoint.

    `on='atom'` — keys are stable atom ids; label is offset away from the mean of neighbours.
    `on='bond'` — keys are (n, m) bond pairs; label sits at the bond midpoint, offset perpendicularly.
    `fmt=None`  — `style.field.value_format`, applied as `fmt.format(value)`.
    `color='auto'` — `style.field.value_colour`; a hex string overrides it.
    `size=None` — `style.label.size * style.field.value_scale`, the same product the colorbar sizes its
                  tick text with, so a figure's numbers are all one size.
    """
    values: dict = field(default_factory=dict)
    on: str = 'atom'
    fmt: str | None = None
    size: float | None = None
    color: str = 'auto'

    def __post_init__(self):
        if self.on not in ('atom', 'bond'):
            raise ValueError(f"ValueLabels.on must be 'atom' or 'bond', got {self.on!r}")
        object.__setattr__(self, 'values', dict(self.values))
        # `on` is inferred as 'bond' when every key is a tuple.
        has_tuple_keys = any(isinstance(k, tuple) for k in self.values)
        has_int_keys = any(isinstance(k, int) for k in self.values)
        if has_tuple_keys and not has_int_keys:
            object.__setattr__(self, 'on', 'bond')
        elif self.on == 'bond' and has_int_keys:
            raise ValueError(
                "ValueLabels.on='bond' but an int key was given; "
                "bond keys must be (n, m) tuples")

    @staticmethod
    def _format(fmt: str, value: float) -> str:
        """Format value and replace a leading ASCII minus with the typographic minus."""
        s = fmt.format(value)
        if s.startswith('-'):
            s = MINUS + s[1:]
        return s

    def render(self, mol, plane, boxes, style: DepictStyle, *, cycle: int = 0,
               log=None) -> tuple[list, list]:
        if not self.values:
            return [], []

        lbl_size = self.size if self.size is not None else style.label.size * style.field.value_scale
        fmt = self.fmt if self.fmt is not None else style.field.value_format
        colour = style.field.value_colour if self.color == 'auto' else to_hex(self.color)

        over = []

        if self.on == 'atom':
            for sid, val in self.values.items():
                ax, ay = plane[sid]
                text_str = self._format(fmt, val)
                dx, dy = _free_direction(mol, plane, sid)
                atom_box = boxes[sid].box if sid in boxes else Box(ax, ay, ax, ay)
                # HALF the larger box side is the box's own radius along the offset direction (the box is
                # centred on the atom), plus `_VALUE_LABEL_GAP` label sizes of air.
                clearance = max(atom_box.width, atom_box.height) * .5 + lbl_size * _VALUE_LABEL_GAP
                lx = ax + dx * clearance
                ly = ay + dy * clearance
                baseline = ly - style.label.baseline_shift * lbl_size
                text = Text([TextRun(text_str, family=style.label.family, size=lbl_size)],
                            x=lx, y=baseline, anchor='middle', fill=colour)
                text = _nudge_clear(text, atom_box, dx, dy, lbl_size)
                over.append(text)

        else:  # on='bond'
            for (n, m), val in self.values.items():
                nx, ny = plane[n]
                mx, my = plane[m]
                mid_x = (nx + mx) / 2.
                mid_y = (ny + my) / 2.
                text_str = self._format(fmt, val)
                bx, by = mx - nx, my - ny
                blen = hypot(bx, by)
                if blen > 1e-9:
                    perp_x, perp_y = -by / blen, bx / blen
                else:
                    perp_x, perp_y = 0., 1.
                offset = lbl_size * _BOND_LABEL_OFFSET
                lx = mid_x + perp_x * offset
                ly = mid_y + perp_y * offset
                baseline = ly - style.label.baseline_shift * lbl_size
                text = Text([TextRun(text_str, family=style.label.family, size=lbl_size)],
                            x=lx, y=baseline, anchor='middle', fill=colour)
                over.append(text)

        return [], over


def _nudge_clear(text: Text, box: Box, dx: float, dy: float, size: float) -> Text:
    """Shift text further along (dx, dy) until its bounds no longer overlap box.

    `_NUDGE_STEP` and `_NUDGE_LIMIT` are algorithm bounds, not style: see their definitions.
    """
    step = size * _NUDGE_STEP
    for _ in range(_NUDGE_LIMIT):
        b = text.bounds
        if (b.max_x <= box.min_x or b.min_x >= box.max_x or
                b.max_y <= box.min_y or b.min_y >= box.max_y):
            break
        text = text.translated(dx * step, dy * step)
    return text


Overlay = Highlight | AtomHalo | AtomField | BondScale | ValueLabels


def render_overlays(overlays: list[Overlay], mol, plane, boxes, style: DepictStyle,
                    *, log=None) -> tuple[list, list]:
    """Render each overlay in the given order and concatenate their (under, over) lists.

    The caller's order is painter's order for the under list: the first overlay is drawn first
    (lowest), the last overlay is drawn last (topmost).  The code never reorders.
    """
    if not overlays:
        return [], []
    all_under: list = []
    all_over: list = []
    for i, overlay in enumerate(overlays):
        u, o = overlay.render(mol, plane, boxes, style, cycle=i, log=log)
        all_under.extend(u)
        all_over.extend(o)
    return all_under, all_over


def scale_of(overlay: Overlay) -> Colormap | None:
    """Return the fitted colormap an overlay draws with, or None for one that carries no scale.

    `colorbar()` labels this scale, through `_fitted` — the SAME call each `render` makes.  An overlay that
    draws no colour answers None: a `Highlight` (its colour is stated, not mapped), a `ValueLabels`, a
    valueless overlay, and a `BondScale` encoding only width, which is the default.
    """
    if not isinstance(overlay, (AtomHalo, AtomField, BondScale)) or not overlay.values:
        return None
    if isinstance(overlay, BondScale) and overlay.encode == 'width':
        return None
    return _fitted(overlay)
