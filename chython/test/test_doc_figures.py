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
"""A figure in ``docs/`` is output of the sample above it, and this is what holds that true.

``test_doc_samples.py`` proves a sample runs; it says nothing about the picture beside it, and a picture
that no longer matches its sample is the documentation defect a reader cannot detect -- the code says one
thing and the image shows another.  ``docs/figures.py`` renders each page's figures and 3D
scenes from the page itself; here that rendering is compared against the committed files.

``to_svg`` and ``depict3d`` are deterministic and the glyph metrics are shipped TSVs rather than system
fonts, so nothing outside the checkout enters the output -- EXCEPT the plane, when the QuickJS layout
computed it: the engine answers ``Math.sin`` and its neighbours from the platform's libm, so a
coordinate differs in its last bits between hosts and byte equality would be a property of the runner.
The comparison is therefore ``figures.same_drawing``: the text around the numbers must match exactly
and every number to within ``figures.TOLERANCE``.  A failure is fixed by running
``python docs/figures.py``, never by editing an SVG or a scene.
"""

from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
from pytest import approx, mark, skip


def _doc_root():
    """``docs/`` beside the repository's ``chython/``, or ``None`` in an installed package."""
    for parent in Path(__file__).resolve().parents:
        if (parent / 'chython').is_dir() and (parent / 'docs').is_dir():
            return parent / 'docs'
    return None


def _figures():
    """``docs/figures.py`` as a module, loaded by path -- ``docs/`` is not a package and never will be."""
    root = _doc_root()
    if root is None or not (root / 'figures.py').is_file():
        skip('no docs/figures.py beside this package -- an installed copy, not a checkout')
    spec = spec_from_file_location('_doc_figures', root / 'figures.py')
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _pages_with_figures():
    """Cheap text match rather than `items()`: this runs at collection, where a skip is an error.

    `encoding='utf-8'` on every read here and below: the sources are UTF-8 and `read_text` without it
    asks the locale, which is cp1252 on Windows -- and `docs/depiction.rst` holds a `⁻`, so the default
    turns collection of this module into an error rather than a test result.
    """
    root = _doc_root()
    if root is None:
        return []
    return sorted(p for p in root.glob('*.rst')
                  if '.. figure:: images/' in (text := p.read_text(encoding='utf-8'))
                  or ':file: scenes/' in text)


@mark.parametrize('page', [p.name for p in _pages_with_figures()] or ['<no docs/>'])
def test_every_referenced_figure_exists(page):
    """The cheap half: a missing file renders as a broken image and Sphinx only warns.

    Every referenced file, including one whose sample this host skips: a `:skipif:` says the sample
    cannot run here, not that the picture is optional -- the page shows it to every reader.
    """
    root = _doc_root()
    if root is None:
        skip('no docs/ beside this package -- an installed copy, not a checkout')
    figures = _figures()

    missing = [f'{page}:{line} -> {figures.target(kind, name).relative_to(root).as_posix()}'
               for kind, name, line, _ in figures.items(root / page)
               if kind in figures._ASSETS and not figures.target(kind, name).is_file()]
    assert not missing, ('these figures are referenced and not committed; run `python docs/figures.py`: '
                         f'{", ".join(missing)}')


def test_no_committed_image_is_unreferenced():
    """A file no page shows is one no sample can keep true, so it is a stale artefact."""
    root = _doc_root()
    if root is None:
        skip('no docs/ beside this package -- an installed copy, not a checkout')
    figures = _figures()

    orphans = []
    for kind, asset in figures._ASSETS.items():
        directory = root / asset.directory
        if not directory.is_dir():
            continue
        referenced = {name for page in root.glob('*.rst')
                      for k, name, _, _ in figures.items(page) if k == kind}
        orphans.extend(f'{asset.directory}/{p.name}' for p in sorted(directory.glob(f'*{asset.suffix}'))
                       if p.stem not in referenced)
    assert not orphans, f'no page references these; delete them: {", ".join(orphans)}'


def test_the_comparison_still_separates_noise_from_a_different_drawing():
    """The control, and the one thing here whose failure mode is a green module.

    ``same_drawing`` is what every assertion below reads, so a comparison that answered ``None`` for
    anything would pass every page while comparing nothing -- and unlike a byte test it cannot be
    inspected by reading it.  The four cases are the four verdicts it has to give.
    """
    figures = _figures()
    svg = '<path d="M -1.4290 2.0 L 0.5 1.0" fill="#000000"/>'

    # the print boundary, which is what the tolerance exists to absorb: one four-decimal step
    assert figures.same_drawing(svg, svg) is None
    assert figures.same_drawing(svg, svg.replace('-1.4290', '-1.4289')) is None
    # a layout that took the other branch of a tie: most of a bond length, and reported as a number
    assert figures.same_drawing(svg, svg.replace('-1.4290', '-0.6040')) == approx(0.825)
    # not a number at all: a different element, a lost path, a changed style
    assert figures.same_drawing(svg, svg.replace('fill', 'stroke')) == float('inf')
    assert figures.same_drawing(svg, svg.replace(' L 0.5 1.0', '')) == float('inf')


def test_the_content_comparison_drops_the_geometry_and_nothing_else():
    """The second control, for the tier a host that did not render the file is held to.

    ``same_content`` answering ``True`` too easily is the failure mode with no symptom: every page would
    pass while only the geometry differed, which is the one thing it is allowed to ignore.
    """
    figures = _figures()
    svg = ('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 4 2"><path d="M 0 0 L 1 1" '
           'fill="none" stroke="#000000"/><text x="1" y="1" text-anchor="end">O</text></svg>')

    # the geometry, including which side of an atom its label sits on
    assert figures.same_content(svg, svg)
    assert figures.same_content(svg, svg.replace('M 0 0 L 1 1', 'M 0 0 L 0.7 1.3'))
    assert figures.same_content(svg, svg.replace(' text-anchor="end"', ''))
    # everything else: a label, an element, a style, a colour
    assert not figures.same_content(svg, svg.replace('>O<', '>N<'))
    assert not figures.same_content(svg, svg.replace('<text x="1" y="1" text-anchor="end">O</text>', ''))
    assert not figures.same_content(svg, svg.replace('stroke="#000000"', 'stroke="#0000cc"'))
    assert not figures.same_content(svg, svg.replace('fill="none"', 'fill="#000000"'))
    # a scene is an HTML fragment: unparsable is not a match, or two of them would pass as one drawing
    assert not figures.same_content('<div>a scene', '<div>a scene')


@mark.parametrize('page', [p.name for p in _pages_with_figures()] or ['<no docs/>'])
def test_every_figure_matches_its_sample(page):
    """The whole point: re-render the page and compare with the committed file.

    One test per page, because ``render`` runs a page's blocks in one shared namespace -- the same
    all-or-nothing unit ``test_doc_samples.py`` uses, and for the same reason.

    NOT BYTE FOR BYTE, and NOT the geometry on a host that did not render the committed file.  The
    plane comes from the QuickJS layout, which answers ``Math.sin`` and friends from the platform's
    libm: a coordinate written with four decimals lands on the print boundary often enough that byte
    equality is a property of the runner, and a tie inside the engine can take the other branch, which
    moves a fragment about a quarter of a bond and flips the side its label sits on.  So a figure passes
    when it is ``figures.same_drawing`` -- the geometry within ``TOLERANCE``, which is what
    ``docs/figures.py --check`` asserts on the host that renders -- or, failing that, when it is
    ``figures.same_content``: the same elements, labels and colours with the geometry dropped.  A lost
    atom, a changed element, a dropped bond and a recoloured one all still fail, on every host.  The
    message names the geometric magnitude either way, so a drawing that moved is legible in the log.

    A block carrying a ``:skipif:`` that holds draws nothing and is compared to nothing -- the same
    option ``test_doc_samples.py`` honours, so a sample needing an optional extra is skipped on a host
    without it rather than raising the ``ImportError`` the extra exists to name.
    """
    root = _doc_root()
    if root is None:
        skip('no docs/ beside this package -- an installed copy, not a checkout')
    figures = _figures()

    stale = []
    for (kind, name), text in figures.render(root / page).items():
        if text is None:
            continue                # the sample's `:skipif:` holds here, so this host drew nothing
        target = figures.target(kind, name)
        if not target.is_file():
            stale.append('%s (not committed)' % target.relative_to(root).as_posix())
            continue
        committed = target.read_text(encoding='utf-8')
        worst = figures.same_drawing(committed, text)
        if worst is not None and not figures.same_content(committed, text):
            stale.append('%s (%s)' % (target.relative_to(root).as_posix(),
                                      'structure' if worst == float('inf') else 'by %.4g' % worst))
    assert not stale, (f'{page} draws these differently than the committed file; run '
                       f'`python docs/figures.py`: {", ".join(stale)}')
