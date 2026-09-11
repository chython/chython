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
"""Regenerates `docs/images/` and `docs/scenes/` from the samples that claim to have drawn them.

    python docs/figures.py            # write every asset
    python docs/figures.py --check    # exit 1 naming the ones that no longer match

THE PICTURE IS NOT A SECOND SOURCE.  An asset is whatever its variable holds when the `testcode::`
block above the reference finishes, so a sample and its picture cannot drift: editing the sample and
regenerating is the only way to change the image.  The pairing is positional and is exactly what a
reader sees -- one asset per preceding block, no naming rule.

| Asset  | Variable | Committed as    | Referenced by                 |
| ------ | -------- | --------------- | ----------------------------- |
| figure | `svg`    | `images/*.svg`  | `.. figure::` / `.. image::`  |
| scene  | `xml`    | `scenes/*.html` | `.. raw:: html` with `:file:` |

`to_svg` and `depict3d` are deterministic given a plane, so `--check` compares text rather than
images -- but a plane the QuickJS layout computed is NOT bit-identical across hosts, so the comparison
is `same_drawing` below and not equality.  `chython/test/test_doc_figures.py` runs it.
"""
from pathlib import Path
from re import compile as re_compile
from sys import argv, exit as _exit, path as _path
from tempfile import TemporaryDirectory
from typing import NamedTuple


DOC = Path(__file__).resolve().parent

# The script lives beside the page, so `sys.path[0]` is `docs/` and an installed chython would answer
# instead of the checkout the samples are documenting -- and a block runs with the cwd in a scratch
# directory, where `''` would answer with nothing.
if (_root := str(DOC.parent)) not in _path:
    _path.insert(0, _root)

#: The directives whose bodies run, in the order a page's blocks must run -- `test_doc_samples.py`
#: executes the same two, for the same reason.
_EXECUTED = ('testsetup', 'testcode')

#: The directives that reference a figure.  `image` and `figure` differ only in whether a caption
#: follows, and the caption is not this script's business.
_FIGURES = ('image', 'figure')


class _Asset(NamedTuple):
    variable: str       # the local a block leaves it in, so no block has to say which of its locals it is
    directory: str      # under `docs/`
    suffix: str


#: One entry per kind of output a page can carry.  Keyed by the kind `items()` reports.
_ASSETS = {'figure': _Asset('svg', 'images', '.svg'), 'scene': _Asset('xml', 'scenes', '.html')}

#: A scene's `<x3d>` sizes to its parent, so it needs a box with a height -- the same wrapper
#: `JupyterWidget._repr_html_` puts around the same document.  The runtime `<script>` and `<link>` are
#: NOT here: `conf.py` loads them once per page, where a second scene cannot load x3dom twice.
_SCENE_BOX = '<div class="x3d-scene" style="width: 100%%; height: %s">%s</div>\n'
_SCENE_HEIGHT = '360px'


def _reference(kind: str, argument: str, line: int):
    """`(kind, name, line, {})` if `argument` names that kind's committed file, else None."""
    asset = _ASSETS[kind]
    prefix = asset.directory + '/'
    if argument.startswith(prefix) and argument.endswith(asset.suffix):
        return kind, argument[len(prefix):-len(asset.suffix)], line, {}
    return None


def items(path: Path):
    """`[(kind, payload, line, options)]` for one page: `('code', source)` and `(asset kind, 'name')`.

    A directive's body is the indented run after it, dedented to its own first-line indent; its option
    lines become `options`, which is empty for an asset reference.  Only the kinds above are collected;
    everything else is prose.

    `encoding='utf-8'` here and on both ends of the asset comparison below: the sources and the assets
    are UTF-8, and `read_text`/`write_text` without it ask the locale -- cp1252 on Windows, where
    `depiction.rst`'s `⁻` cannot be decoded and a rendered `⁻` cannot be written.
    """
    lines = path.read_text(encoding='utf-8').splitlines()
    out = []
    i = 0
    while i < len(lines):
        stripped = lines[i].lstrip()
        if not stripped.startswith('.. ') or '::' not in stripped:
            i += 1
            continue
        head, _, argument = stripped[3:].partition('::')
        directive = head.strip()
        argument = argument.strip()
        if directive in _FIGURES:
            if (found := _reference('figure', argument, i + 1)) is not None:
                out.append(found)
            i += 1
            continue
        if directive == 'raw' and argument == 'html':
            # The only directive whose reference is an OPTION rather than its argument, and the reason
            # `:file:` is used at all: the generated markup stays out of the page.
            i += 1
            while i < len(lines) and (option := lines[i].strip()).startswith(':'):
                if option.startswith(':file:'):
                    if (found := _reference('scene', option[len(':file:'):].strip(), i + 1)) is not None:
                        out.append(found)
                i += 1
            continue
        if directive not in _EXECUTED:
            i += 1
            continue
        outer = len(lines[i]) - len(stripped)
        start = i + 1
        j = i + 1
        options = {}
        while j < len(lines):                      # option lines and the blanks around them
            option = lines[j].strip()
            if not option:
                j += 1
            elif option.startswith(':') and option.count(':') >= 2:
                key, _, value = option[1:].partition(':')
                options[key.strip()] = value.strip()
                j += 1
            else:
                break
        body, base = [], None
        while j < len(lines):
            line = lines[j]
            if not line.strip():
                body.append('')
                j += 1
                continue
            indent = len(line) - len(line.lstrip())
            if base is None:
                if indent <= outer:
                    break
                base = indent
            elif indent < base:
                break
            body.append(line[base:])
            j += 1
        out.append(('code', '\n'.join(body).rstrip(), start, options))
        i = j
    return out


def target(kind: str, name: str) -> Path:
    """Where that asset is committed.  One place computes this, so a writer and a checker agree."""
    asset = _ASSETS[kind]
    return DOC / asset.directory / f'{name}{asset.suffix}'


_NUMBER = re_compile(r'-?\d+(?:\.\d+)?(?:e-?\d+)?')

#: How far a number may move before the committed file is a different drawing, in the units a
#: coordinate is written in -- a bond is 1.0 there and a bond line is 0.04 wide, so this is a quarter
#: of a line width and nothing a reader can see.  It sits between the two magnitudes it has to
#: separate: a coordinate is printed with four decimals, so a value on a rounding boundary moves by
#: 1e-4 when a summation order changes, while a layout that took the other branch of a tie inside the
#: QuickJS engine moves an atom by a large fraction of a bond length.  Everything that is not a number
#: still has to match exactly, so a lost atom, a changed element or a dropped path is not absorbed, and
#: neither is a tie flip -- `same_drawing` returns the magnitude rather than a verdict, so which of the
#: two a failure is gets read off the message instead of assumed from this constant.
TOLERANCE = 0.01


def same_drawing(committed: str, rendered: str) -> float | None:
    """`None` if the two are the same drawing, else the largest difference that says they are not.

    `inf` reports a difference no tolerance can absorb: the text around the numbers differs, or there
    is a different count of them.
    """
    if _NUMBER.sub('#', committed) != _NUMBER.sub('#', rendered):
        return float('inf')
    a, b = _NUMBER.findall(committed), _NUMBER.findall(rendered)
    worst = max((abs(float(x) - float(y)) for x, y in zip(a, b)), default=0.)
    return None if worst <= TOLERANCE else worst


#: Every attribute whose value the plane decides.  `text-anchor` is one of them: which side of an atom
#: its label sits on is a comparison between two coordinates, so it flips wherever the plane does.
_PLACED = frozenset({'d', 'x', 'y', 'x1', 'y1', 'x2', 'y2', 'cx', 'cy', 'r', 'rx', 'ry', 'points',
                     'transform', 'viewBox', 'width', 'height', 'text-anchor', 'dx', 'dy'})


def same_content(committed: str, rendered: str) -> bool:
    """Whether the two drawings state the same thing, with the geometry dropped.

    What is left is what no host can disagree about: the element sequence, the text of every label, and
    the attributes that are style rather than position.  A lost atom, a changed element, a dropped bond
    or a recoloured one still fails; a fragment that converged elsewhere does not.

    This is the comparison a host that did not render the committed file can make.  Measured on the
    GitHub runners, two of the 54 figures differ across them by more than any usable tolerance --
    `smirks-stereo-racemic.svg` by 0.2848 in 36 of its 770 numbers on Linux and Windows both, and
    `smirks-masked.svg` by 0.6228 on Windows, which also drops a `text-anchor`.  The engine answers
    `Math` from the platform's libm, so one fragment takes the other branch of a tie and lands about a
    quarter of a bond away, and the label side goes with it.  Neither host is the reference: for the
    first figure macOS is the odd one out and for the second Windows is.
    """
    committed_content = _content(committed)
    # `None` on either side is not a match: two scenes nobody could parse are not thereby the same one
    return committed_content is not None and committed_content == _content(rendered)


def _content(svg: str) -> list | None:
    """`svg` as `[(tag, style attributes, text), ...]` in document order, or `None` if it is not XML."""
    from xml.etree.ElementTree import ParseError, fromstring

    try:
        root = fromstring(svg)
    except ParseError:      # a scene is an HTML fragment, which only `same_drawing` speaks about
        return None
    return [(e.tag, tuple(sorted((k, v) for k, v in e.attrib.items() if k not in _PLACED)),
             (e.text or '').strip()) for e in root.iter()]


def body(kind: str, value: str) -> str:
    """What the committed file holds: an SVG verbatim, a scene inside a box that has a height."""
    if kind == 'scene':
        return _SCENE_BOX % (_SCENE_HEIGHT, value)
    return value


def _gated(page: Path, line: int, condition: str, namespace: dict) -> bool:
    """Whether a block's `:skipif:` holds, evaluated in the page's namespace as Sphinx would.

    `test_doc_samples.py` reads the same option the same way, so a page states the condition once and
    both harnesses agree about which hosts run the block.  An unevaluable condition is a defect in the
    page and not a skip: read as one, a mistyped name would quietly stop drawing a figure everywhere.
    """
    try:
        return bool(eval(condition, dict(namespace)))
    except Exception as e:
        raise SystemExit(f'{page.name}:{line} has an unevaluable :skipif: -- {e!r}')


def render(page: Path) -> dict[tuple[str, str], str | None]:
    """`{(kind, name): file body}` for one page, by running its blocks in order in one namespace.

    `None` is an asset this host cannot draw: the block above it carries a `:skipif:` that holds, which
    is how a page states a sample needing an optional extra.  The committed file stays as it is -- a
    host without the extra has nothing to say about it, and neither writing nor calling it stale is that.
    """
    from os import chdir, getcwd

    namespace = {'__name__': f'figures_{page.stem}'}
    drawn: dict[tuple[str, str], str | None] = {}
    pending: dict[str, str] = {}                   # each kind's variable from the last block, if a string
    skipped = False                                # whether that last block is one this host skipped
    was = getcwd()
    with TemporaryDirectory() as scratch:
        try:
            chdir(scratch)                         # a sample that writes a file writes it here
            for kind, payload, line, options in items(page):
                if kind == 'code':
                    pending = {}
                    if skipped := 'skipif' in options and _gated(page, line, options['skipif'], namespace):
                        continue
                    exec(compile(payload, f'{page.name}:{line}', 'exec'), namespace)
                    for other, asset in _ASSETS.items():
                        value = namespace.pop(asset.variable, None)  # popped, so the NEXT block cannot
                        if isinstance(value, str):                   # inherit this picture
                            pending[other] = value
                elif (kind, payload) in drawn:
                    raise SystemExit(f'{page.name}:{line} draws '
                                     f'{target(kind, payload).relative_to(DOC).as_posix()} twice')
                elif kind in pending:
                    drawn[kind, payload] = body(kind, pending.pop(kind))
                elif skipped:
                    drawn[kind, payload] = None
                else:
                    raise SystemExit(f'{page.name}:{line} references '
                                     f'{target(kind, payload).relative_to(DOC).as_posix()} but the block '
                                     f'above it left no `{_ASSETS[kind].variable}` string')
        finally:
            chdir(was)
    return drawn


def main(check: bool) -> int:
    pages = sorted(DOC.glob('*.rst'))
    stale, written, kept, undrawn = [], 0, 0, 0
    for page in pages:
        if not any(kind in _ASSETS for kind, _, _, _ in items(page)):
            continue
        for (kind, name), text in render(page).items():
            path = target(kind, name)
            if text is None:                       # `:skipif:` held: this host draws nothing to compare
                undrawn += 1
                continue
            # a file within `TOLERANCE` is the same drawing, so it is neither stale nor rewritten:
            # rewriting one would put a host's last digits in the diff of every regeneration
            worst = float('inf') if not path.exists() \
                else same_drawing(path.read_text(encoding='utf-8'), text)
            if worst is None:
                kept += 1
            elif check:
                stale.append((path.relative_to(DOC.parent).as_posix(), worst))
            else:
                path.parent.mkdir(exist_ok=True)
                path.write_text(text, encoding='utf-8')
                written += 1
                print(f'{path.relative_to(DOC.parent).as_posix():44s} {len(text) // 1024:4d} KiB')
    # what this host could not draw is reported by count on every run: a figure silently left out of the
    # comparison is how a `:skipif:` widened by mistake would stop being checked anywhere
    missing = f', {undrawn} not drawn here (a `:skipif:` holds)' if undrawn else ''
    if check:
        if stale:
            print('these no longer match their samples; run `python docs/figures.py`:')
            for name, worst in stale:
                print(f'  {name}  ({"structure" if worst == float("inf") else f"by {worst:.4g}"})')
            return 1
        print(f'every figure and scene matches its sample ({kept} within {TOLERANCE}{missing})')
        return 0
    print(f'\n{written} written, {kept} already the same drawing{missing}')
    return 0


if __name__ == '__main__':
    _exit(main('--check' in argv[1:]))
