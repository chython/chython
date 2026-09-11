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
"""Every rendering constant the depictor uses, as an immutable tree of frozen dataclasses.

A style is a value: hashable, nested by meaning (bonds, atoms, labels, the page), passed explicitly to
whatever draws, and validated at construction.  `tuned()` is dotted-key sugar over `replace`.
Lengths are in molecule units -- one standard bond is 1.0 -- except `*_mm`, which is page geometry.
"""
from dataclasses import dataclass, fields, replace

from .metrics import FAMILIES
from .scene import to_hex


__all__ = ['AtomStyle', 'BondStyle', 'ContourStyle', 'DepictStyle', 'FieldStyle', 'HighlightStyle',
           'LabelStyle', 'PageStyle', 'PRESETS', 'ReactionStyle', 'get_depict_style',
           'set_depict_style']


def _positive(value, name):
    if value <= 0.:
        raise ValueError(f'{name} must be positive, got {value}')
    return float(value)


def _non_negative(value, name):
    if value < 0.:
        raise ValueError(f'{name} must not be negative, got {value}')
    return float(value)


@dataclass(frozen=True, slots=True)
class BondStyle:
    """Line weights and the geometry of a multiple bond.  Lengths in molecule units."""
    width: float = .04
    spacing: float = .18                  # centre-to-centre of a double bond's two lines
    triple_spacing: float = .16           # offset of each outer line FROM THE AXIS; the spread is 2x
    aromatic: str = 'dashed-inner'        # 'kekule' | 'circle' | 'dashed-inner'
    aromatic_inset: float = .22           # the CIRCLE's inset from the ring bonds
    # Two insets: the dashed inner line sits closer to its bond than the circle does.
    aromatic_dash_inset: float = .14      # the dashed inner line's inset from its ring bond
    aromatic_dashes: tuple[float, ...] = (.15, .05)
    wedges: str = 'stored'                # 'stored': draw the wedge marks the structure carries;
                                          # 'recompute': derive them from the stereo and this layout
    wedge_width: float = .16              # the wide end of a stereo wedge
    hash_step: float = .09                # rung pitch of a hashed wedge
    either_amplitude: float = .07         # the squiggle of an "either" bond
    either_period: float = .20            # and its wavelength
    # A dative bond is a dashed axis with no arrow head: a container does not store which atom donates.
    dative_dashes: tuple[float, ...] = (.2, .1)
    trim: float = .06                     # extra clearance between a bond end and a label's ink box
    join: str = 'miter'
    # Round, because a bond is drawn as one or more separate paths and their ends have to MEET: a chain
    # broken by a label, a ring's inner line, a double bond's second line and a wedge all stop at their
    # own end, and two butt caps arriving at one point from two angles show the notch between them.  A
    # round cap is the same disc whatever the angle, so the seam closes.  It also softens the free end of
    # a terminal bond, which is where a butt cap reads as a cut.  `_dash_pattern` compensates the dashed
    # patterns for it, since a round cap lengthens every dash by half the stroke width at each end.
    cap: str = 'round'
    miter_limit: float = 4.
    colour: str = '#000000'

    def __post_init__(self):
        object.__setattr__(self, 'width', _positive(self.width, 'bond width'))
        object.__setattr__(self, 'spacing', _positive(self.spacing, 'bond spacing'))
        object.__setattr__(self, 'triple_spacing', _positive(self.triple_spacing, 'triple spacing'))
        object.__setattr__(self, 'wedge_width', _positive(self.wedge_width, 'wedge width'))
        object.__setattr__(self, 'hash_step', _positive(self.hash_step, 'hash step'))
        object.__setattr__(self, 'trim', _non_negative(self.trim, 'bond trim'))
        object.__setattr__(self, 'aromatic_inset', _non_negative(self.aromatic_inset,
                                                                 'aromatic inset'))
        object.__setattr__(self, 'aromatic_dash_inset', _non_negative(self.aromatic_dash_inset,
                                                                      'aromatic dash inset'))
        object.__setattr__(self, 'colour', to_hex(self.colour))
        if self.aromatic not in ('kekule', 'circle', 'dashed-inner'):
            raise ValueError(f'aromatic must be kekule, circle or dashed-inner, got {self.aromatic!r}')
        if self.wedges not in ('stored', 'recompute'):
            raise ValueError(f'wedges must be stored or recompute, got {self.wedges!r}')


@dataclass(frozen=True, slots=True)
class AtomStyle:
    """When an atom gets a label at all, and what colour it is."""
    carbon: bool = False                  # label plain carbons?
    hydrogens: bool = True                # write implicit H counts on labelled atoms
    unknown_h_marks: bool = False         # mark an unknown H count with '?' rather than silence
    charges: bool = True
    isotopes: bool = True
    radicals: bool = True
    map_numbers: bool = True              # drawn only where `map_number` is non-zero
    stereo_labels: bool = False           # the stored CIP descriptor beside a centre
    stereo_groups: bool = True            # `&N` (AND), `oN` (OR), `a` (ABS) from `atom.stereo_group`
    query_marks: bool = True              # a query atom's primitives
    colour_by_element: bool = True
    carbon_colour: str = '#000000'
    default_colour: str = '#000000'
    radical_radius: float = .045
    radical_gap: float = .09

    def __post_init__(self):
        object.__setattr__(self, 'carbon_colour', to_hex(self.carbon_colour))
        object.__setattr__(self, 'default_colour', to_hex(self.default_colour))
        object.__setattr__(self, 'radical_radius', _positive(self.radical_radius, 'radical radius'))
        object.__setattr__(self, 'radical_gap', _positive(self.radical_gap, 'radical gap'))


@dataclass(frozen=True, slots=True)
class LabelStyle:
    """Type: family, sizes, and the knock-out that keeps a bond out of a symbol's ink."""
    family: str = 'helvetica'
    size: float = .40
    subscript_scale: float = .70
    superscript_scale: float = .70
    subscript_drop: float = .28           # fraction of `size`, downward
    superscript_rise: float = .40
    stereo_scale: float = .62
    stereo_italic: bool = True
    map_scale: float = .55
    map_colour: str = '#0000cc'
    # Two fixed rows beside the label: the stereo statement above the baseline, the map number below.
    # Fixed, not negotiated -- a map number must not move because the atom also carries a descriptor.
    annotation_rise: float = .40          # fraction of `size`, upward
    annotation_drop: float = .40          # fraction of `size`, downward
    # The knock-out UNDER an annotation: a plate in the background colour between the number and the line
    # it would otherwise sit on.  An annotation STAYS BESIDE ITS ATOM -- a number far enough out to be
    # clear of everything has stopped saying which atom it belongs to -- so the plate, not distance, is
    # what makes the crowded one readable.  Drawn under every annotation and not per number measured
    # against the paths; `_annotation_plates` in `figure.py` gives the reason.
    annotation_plate: str = 'rounded'     # 'rounded' | 'ellipse' | 'none'
    annotation_plate_pad: float = .02     # around the annotation's ink, molecule units
    # None means `page.background`, and white where that is transparent, which is what a knock-out on an
    # unpainted page has to assume: the plate's whole job is to be the colour of what is behind it.
    annotation_plate_colour: str | None = None
    pad: float = .05                      # inflation of the measured ink box, molecule units
    baseline_shift: float = .34           # fraction of `size` that centres a cap on the atom point

    def __post_init__(self):
        if self.family not in FAMILIES:
            raise ValueError(f'unknown label family {self.family!r}: expected one of {sorted(FAMILIES)}')
        object.__setattr__(self, 'size', _positive(self.size, 'label size'))
        object.__setattr__(self, 'pad', _non_negative(self.pad, 'label pad'))
        object.__setattr__(self, 'map_colour', to_hex(self.map_colour))
        object.__setattr__(self, 'annotation_rise',
                           _non_negative(self.annotation_rise, 'annotation rise'))
        object.__setattr__(self, 'annotation_drop',
                           _non_negative(self.annotation_drop, 'annotation drop'))
        if self.annotation_plate not in ('rounded', 'ellipse', 'none'):
            raise ValueError(f'annotation_plate must be rounded, ellipse or none, '
                             f'got {self.annotation_plate!r}')
        object.__setattr__(self, 'annotation_plate_pad',
                           _non_negative(self.annotation_plate_pad, 'annotation plate pad'))
        object.__setattr__(self, 'annotation_plate_colour', to_hex(self.annotation_plate_colour))


@dataclass(frozen=True, slots=True)
class HighlightStyle:
    """Halos around highlighted atoms and ribbons along highlighted bonds."""
    radius: float = .30                   # halo radius around a labelled atom
    bond_width: float = .34               # ribbon width along a highlighted bond
    opacity: float = .45
    outline: str | None = None
    outline_width: float = .02
    # Wong, B. "Points of view: Color blindness."  Nature Methods 8, 441 (2011) -- eight hues that stay
    # distinguishable under deuteranopia, protanopia and tritanopia.  A `Highlight` with no colour of
    # its own reads `palette[i % len(palette)]` for its position `i`; the neutral grey is first.
    palette: tuple[str, ...] = ('#767676', '#e69f00', '#56b4e9', '#009e73',
                                '#f0e442', '#0072b2', '#d55e00', '#cc79a7')

    def __post_init__(self):
        object.__setattr__(self, 'radius', _positive(self.radius, 'highlight radius'))
        object.__setattr__(self, 'bond_width', _positive(self.bond_width, 'highlight bond width'))
        object.__setattr__(self, 'outline_width', _positive(self.outline_width, 'highlight outline width'))
        if not 0. <= self.opacity <= 1.:
            raise ValueError(f'highlight opacity must be within [0, 1], got {self.opacity}')
        if not self.palette:
            raise ValueError('the highlight palette must have at least one colour')
        object.__setattr__(self, 'outline', to_hex(self.outline))
        object.__setattr__(self, 'palette', tuple(to_hex(c) for c in self.palette))


@dataclass(frozen=True, slots=True)
class ContourStyle:
    """How a field's isolines are traced.  The band count is not here: it is `AtomField(levels=...)`."""
    grid: float = .12                     # marching-squares cell size, molecule units
    smooth: bool = True                   # fit tangent-continuous cubics instead of emitting segments
    refine: int = 2                       # Newton steps along the gradient per traced vertex
    line_width: float = .018
    fill_opacity: float = .7              # opacity of the whole field group, read by
                                          # `AtomField.opacity=None`.  There is no `line_opacity`: a
                                          # band and its boundary stroke are sibling Paths in one
                                          # Group, and group opacity is the scene's only compositing
                                          # mechanism.
    sigma: float = .55                    # Gaussian width scale factor: effective sigma = sigma *
                                          # mean bond length of the plane, multiplied in
                                          # `AtomField.render`.  Weight at one bond length is
                                          # exp(-1/(2*.55²)) ≈ 0.19, so neighbours merge but per-atom
                                          # structure survives.
    cutoff: float = 4.                    # beyond this from every named atom the field is UNDEFINED,
                                          # not zero -- `ScalarField.at` returns None and it draws as
                                          # nothing.  An absolute distance (not scaled), so at clean2d
                                          # scale (bond≈0.825) it is ~5 bonds away.
    pad: float = .55                      # how far past the atoms the field is sampled, and the
                                          # padding of the convex hull when `clip='hull'`

    def __post_init__(self):
        if self.refine < 0:
            raise ValueError(f'contour refine steps must not be negative, got {self.refine}')
        object.__setattr__(self, 'grid', _positive(self.grid, 'contour grid'))
        object.__setattr__(self, 'sigma', _positive(self.sigma, 'contour sigma'))
        object.__setattr__(self, 'cutoff', _positive(self.cutoff, 'contour cutoff'))
        object.__setattr__(self, 'line_width', _positive(self.line_width, 'contour line width'))
        if self.pad < 0.:
            raise ValueError(f'contour pad must not be negative, got {self.pad}')
        if self.cutoff < self.sigma:
            raise ValueError(f'contour cutoff {self.cutoff} is inside sigma {self.sigma}: the field '
                             'would be truncated where it is still strong, which draws a visible '
                             'circular edge around every atom.  This is a conservative check on the '
                             'style values themselves; the effective sigma is per-plane (sigma_scale '
                             '× mean bond length) and is checked at render time by AtomField.')
        if not 0. <= self.fill_opacity <= 1.:
            raise ValueError(f'contour fill_opacity must be within [0, 1], got {self.fill_opacity}')


@dataclass(frozen=True, slots=True)
class FieldStyle:
    """Scalar-data channels: the colour map, the value labels, and the contour sub-tree."""
    colormap: str = 'coolwarm'
    contour: ContourStyle = ContourStyle()
    halo_min_radius: float = .10          # read by `AtomHalo.render`
    halo_max_radius: float = .34
    bond_min_width: float = .02           # read by `BondScale.bond_widths` when width_range is None
    bond_max_width: float = .22
    value_scale: float = .52              # label size as a fraction of `LabelStyle.size`; shared by
                                          # `ValueLabels.render` and the colorbar's tick text
    value_format: str = '{:.2f}'          # read by `ValueLabels.render` when fmt is None
    value_colour: str = '#333333'

    def __post_init__(self):
        object.__setattr__(self, 'value_colour', to_hex(self.value_colour))
        if self.halo_min_radius > self.halo_max_radius:
            raise ValueError('halo_min_radius must not exceed halo_max_radius')
        if self.bond_min_width > self.bond_max_width:
            raise ValueError('bond_min_width must not exceed bond_max_width')


@dataclass(frozen=True, slots=True)
class ReactionStyle:
    """Arrow and sign geometry.  The arrangement constants belong to `layout/reaction.py`, not here."""
    arrow_width: float = .05
    head_length: float = .30
    head_width: float = .22
    sign_size: float = .55
    gap: float = .40                      # clearance between a molecule's box and the arrow or sign
    colour: str = '#000000'

    def __post_init__(self):
        object.__setattr__(self, 'arrow_width', _positive(self.arrow_width, 'arrow width'))
        object.__setattr__(self, 'head_length', _positive(self.head_length, 'arrow head length'))
        object.__setattr__(self, 'colour', to_hex(self.colour))


@dataclass(frozen=True, slots=True)
class PageStyle:
    """The one place physical units belong.

    `width_mm` sizes the output; `scale_mm` says how many millimetres one molecule unit becomes.  State
    exactly one -- both is a contradiction and is refused, neither leaves the output with no size.
    """
    width_mm: float | None = None
    scale_mm: float | None = 6.
    margin: float = .35                   # molecule units around the content box
    background: str | None = None         # None means transparent, which is what a figure wants
    legend: str = 'auto'                  # 'auto' | 'right' | 'bottom' | 'none'
    legend_breadth: float = .34           # the short dimension of the swatch strip, molecule units
    legend_fmt: str = '{:+.2f}'           # tick format; a one-sided domain drops the '+' at draw time
    title: bool = False
    title_size: float = .45

    def __post_init__(self):
        if self.width_mm is not None and self.scale_mm is not None:
            raise ValueError('state page width_mm or scale_mm, not both: the two would disagree')
        if self.legend not in ('auto', 'right', 'bottom', 'none'):
            raise ValueError(f'legend must be auto, right, bottom or none, got {self.legend!r}')
        if self.width_mm is None and self.scale_mm is None:
            raise ValueError('state page width_mm or scale_mm: the output needs a physical size')
        if self.width_mm is not None:
            object.__setattr__(self, 'width_mm', _positive(self.width_mm, 'page width_mm'))
        if self.scale_mm is not None:
            object.__setattr__(self, 'scale_mm', _positive(self.scale_mm, 'page scale_mm'))
        object.__setattr__(self, 'margin', _non_negative(self.margin, 'page margin'))
        object.__setattr__(self, 'legend_breadth', _positive(self.legend_breadth, 'legend breadth'))
        try:
            self.legend_fmt.format(-1.5)
        except (ValueError, TypeError):
            raise ValueError(f'legend_fmt {self.legend_fmt!r} cannot format a float')
        object.__setattr__(self, 'background', to_hex(self.background))


@dataclass(frozen=True, slots=True)
class DepictStyle:
    """The root of the tree, and the object every drawing function takes.

    Passed explicitly rather than read from a global, so two pictures in one process can differ.  The
    process default (`get_depict_style`) is only for the convenience entry points -- `mol.depict()`.
    """
    page: PageStyle = PageStyle()
    bond: BondStyle = BondStyle()
    atom: AtomStyle = AtomStyle()
    label: LabelStyle = LabelStyle()
    highlight: HighlightStyle = HighlightStyle()
    field: FieldStyle = FieldStyle()
    reaction: ReactionStyle = ReactionStyle()

    def tuned(self, **settings) -> 'DepictStyle':
        """A copy with dotted keys replaced: `style.tuned(**{'bond.width': .055})`.

        Reaches a sub-branch too, so `field.contour.refine` works.  An unknown key raises `KeyError`
        naming the field it could not find, and every branch is rebuilt through its constructor, so a
        tuned value is validated exactly as a constructed one is.
        """
        if not settings:
            return self
        tree = {}
        for key, value in settings.items():
            if '.' not in key:
                raise KeyError(f'{key!r} is not a dotted style key: write e.g. {"bond." + key!r}')
            head, _, tail = key.partition('.')
            tree.setdefault(head, {})[tail] = value
        return self._apply(tree)

    def _apply(self, tree):
        """The dotted tree from `tuned()`, applied one branch at a time.

        No recursion: the tree is two levels deep, so a sub-branch goes to `_apply_to_leaf`.
        """
        names = {f.name for f in fields(self)}
        changes = {}
        for head, sub in tree.items():
            if head not in names:
                raise KeyError(f'{head} is not a style branch of '
                               f'{type(self).__name__}: expected one of {sorted(names)}')
            branch = getattr(self, head)
            nested = {}
            leaves = {}
            for key, value in sub.items():
                if '.' in key:
                    inner_head, _, inner_tail = key.partition('.')
                    nested.setdefault(inner_head, {})[inner_tail] = value
                else:
                    leaves[key] = value
            if nested:
                branch = _apply_to_leaf(branch, nested, f'{head}.')
            if leaves:
                branch_names = {f.name for f in fields(branch)}
                for key in leaves:
                    if key not in branch_names:
                        raise KeyError(f'{head}.{key} is not a field of '
                                       f'{type(branch).__name__}: expected one of {sorted(branch_names)}')
                branch = replace(branch, **leaves)
            changes[head] = branch
        return replace(self, **changes)

    @classmethod
    def preset(cls, name: str) -> 'DepictStyle':
        """A named starting point, tunable like any other style.

        `DepictStyle.preset('acs').tuned(**{'bond.width': .055})`.  See `PRESETS` for the names.
        """
        try:
            return _PRESETS[name]()
        except KeyError:
            raise ValueError(f'unknown preset {name!r}: expected one of {sorted(_PRESETS)}') from None


def _apply_to_leaf(branch, nested, prefix):
    """Nested `tuned()` one level below a leaf branch -- `field.contour.refine`."""
    changes = {}
    names = {f.name for f in fields(branch)}
    for head, sub in nested.items():
        if head not in names:
            raise KeyError(f'{prefix}{head} is not a field of {type(branch).__name__}: '
                           f'expected one of {sorted(names)}')
        inner = getattr(branch, head)
        inner_names = {f.name for f in fields(inner)}
        for key in sub:
            if key not in inner_names:
                raise KeyError(f'{prefix}{head}.{key} is not a field of {type(inner).__name__}: '
                               f'expected one of {sorted(inner_names)}')
        changes[head] = replace(inner, **sub)
    return replace(branch, **changes)


def _acs() -> DepictStyle:
    """ACS single-column: 83 mm wide (their stated column width), Helvetica, and a bond width above the
    .5 pt line-weight floor their guidelines set."""
    return DepictStyle(page=PageStyle(width_mm=83., scale_mm=None, margin=.25),
                       bond=BondStyle(width=.052, spacing=.20, trim=.07),
                       label=LabelStyle(family='helvetica', size=.42, pad=.055),
                       atom=AtomStyle(stereo_labels=True))


def _print() -> DepictStyle:
    """A Times-set figure for a two-column body, at a fixed scale rather than a fixed width."""
    return DepictStyle(page=PageStyle(width_mm=None, scale_mm=5.5, margin=.3),
                       bond=BondStyle(width=.048),
                       label=LabelStyle(family='times', size=.44))


def _screen() -> DepictStyle:
    """Heavier lines and larger type, for a notebook cell at 96 dpi."""
    return DepictStyle(page=PageStyle(width_mm=None, scale_mm=8., margin=.4),
                       bond=BondStyle(width=.06),
                       label=LabelStyle(size=.45))


def _poster() -> DepictStyle:
    return DepictStyle(page=PageStyle(width_mm=None, scale_mm=14., margin=.5),
                       bond=BondStyle(width=.075, spacing=.22),
                       label=LabelStyle(size=.5, pad=.07))


_PRESETS = {'acs': _acs, 'print': _print, 'screen': _screen, 'poster': _poster}

PRESETS = frozenset(_PRESETS)

# The process default, for the convenience entry points only -- `mol.depict()` with no style argument.
_DEFAULT_STYLE = DepictStyle()


def get_depict_style() -> DepictStyle:
    """The style `mol.depict()` uses when the caller states none."""
    return _DEFAULT_STYLE


def set_depict_style(style: DepictStyle):
    """Set the process default.  Type-checked here, at the assignment, not later at draw time."""
    global _DEFAULT_STYLE

    if not isinstance(style, DepictStyle):
        raise TypeError(f'expected a DepictStyle, got {type(style).__name__}')
    _DEFAULT_STYLE = style
