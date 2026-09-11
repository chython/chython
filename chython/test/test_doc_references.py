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
"""Every name the documentation's **prose** points at exists.

`test_doc_samples.py` executes the samples, which covers the code a page shows and nothing else.  A
`:func:` role, a ``chython.x.y`` in a sentence and a `mol.method()` in a README table are read by no
interpreter: Sphinx renders an unresolvable role as plain text by default, and README is executed by
nothing at all.  So the two references most likely to rot -- a function that moved one package down and
a method that never existed -- rot silently.  Both are checked here by resolving the dotted path and by
asking the container for the attribute.

Scope is `docs/*.rst` plus `README.md`, and the gate is the union: a claim in the release-facing README is
held to the same standard as one in the manual.  Not checked: a bare ``some_function()`` with no owner,
because prose cannot say which namespace it means.

Two documented blind spots, both narrow and both by construction:

* a path whose last segment is private is skipped.  `chython.core._query_boxes` names a ``.pxi``
  fragment, which is a translation-unit member and not a module, so `import` is the wrong question.
* a ``mol.``/``rxn.`` reference whose tail is a file suffix is skipped -- ``molecules.sdf`` in a code
  span and ``mol.state_view()`` in a code span look the same to a regex.
"""

import re
from importlib import import_module
from pathlib import Path
from pytest import mark, skip


#: A dotted path rooted at the package, as prose and roles spell it: ``chython.chemistry.saturate``.
_DOTTED = re.compile(r'\bchython(?:\.[A-Za-z_][A-Za-z_0-9]*)+')

#: A container member by its class name (`:meth:`chython.MoleculeContainer.thiele``, minus the package
#: part, which `_DOTTED` already resolved) or by the variable name every page uses for an instance.
_MEMBER = re.compile(r'\b(MoleculeContainer|ReactionContainer|QueryContainer)\.([A-Za-z_][A-Za-z_0-9]*)'
                     r'|`(mol|rxn|molecule|reaction)\.([A-Za-z_][A-Za-z_0-9]*)')

#: A URL is not an API reference -- ``chython.readthedocs.io`` would resolve to a missing subpackage.
_URL = re.compile(r'https?://\S+')

#: What a member tail may be instead of an attribute: a file name in a code span.  See the docstring.
_SUFFIXES = frozenset(('mol', 'sdf', 'rxn', 'rdf', 'mrv', 'cml', 'smi', 'mol2', 'pdb', 'cif', 'xyz',
                       'json', 'txt', 'gz', 'py', 'rst', 'md', 'svg', 'png'))


def _root():
    """The repository root, or `None` in an installed package."""
    for parent in Path(__file__).resolve().parents:
        if (parent / 'chython').is_dir() and (parent / 'docs').is_dir():
            return parent
    return None


def _sources():
    root = _root()
    if root is None:
        return []
    return sorted(root.glob('docs/*.rst')) + [root / 'README.md']


def _source(root, name):
    """One source by the name the parametrization carries: a page lives in `docs/`, README at the root."""
    return root / name if name == 'README.md' else root / 'docs' / name


def _lines(path):
    """``(line number, text)`` per line, URLs removed.

    `encoding='utf-8'`: the pages are UTF-8 and `docs/depiction.rst` holds a `⁻`, which the locale codec
    Windows hands `read_text` cannot decode.
    """
    return [(i, _URL.sub('', line))
            for i, line in enumerate(path.read_text(encoding='utf-8').splitlines(), 1)]


def _resolves(path):
    """Is `path` an attribute chain from the façade, or a module?

    Both questions, in that order: `chython.chemistry.saturate` is an attribute of a package, while
    `chython.formats.mol2` may not have been imported by the façade yet and is reachable only as a
    module.  A path is what the documentation promises a reader can write, so either answer is a yes.
    """
    import chython

    parts = path.split('.')
    obj = chython
    for i, part in enumerate(parts[1:], 1):
        try:
            obj = getattr(obj, part)
        except AttributeError:
            try:
                obj = import_module('.'.join(parts[:i + 1]))
            except ImportError:
                return False
    return True


@mark.parametrize('source', [p.name for p in _sources()] or ['<no docs/>'])
def test_documented_dotted_paths_exist(source):
    """Every ``chython.x.y`` the prose writes can be imported or reached by attribute.

    [mutant: renaming `chython.interop.conformers.generate_conformers` in `docs/depiction.rst` to the
    path it had before -- this test fails and names the page and the line.]
    """
    root = _root()
    if root is None:
        skip('no docs/ beside this package -- an installed copy, not a checkout')

    bad = []
    for number, text in _lines(_source(root, source)):
        for name in _DOTTED.findall(text):
            if name.rsplit('.', 1)[1].startswith('_'):  # a private tail: see the module docstring
                continue
            if not _resolves(name):
                bad.append(f'{source}:{number} {name}')
    assert not bad, 'documented names that do not exist:\n' + '\n'.join(bad)


@mark.parametrize('source', [p.name for p in _sources()] or ['<no docs/>'])
def test_documented_container_members_exist(source):
    """Every ``mol.method`` and ``MoleculeContainer.method`` the prose names is on the class.

    The instance spellings are the convention the pages already follow -- `mol` is a molecule and `rxn`
    a reaction throughout -- which makes a README table of methods checkable without executing it.

    [mutant: `README.md` claimed a `MoleculeContainer.from_rdkit` classmethod; the conversion is a
    single dispatching `chython.interop.rdkit`.  This test names the line.]
    """
    root = _root()
    if root is None:
        skip('no docs/ beside this package -- an installed copy, not a checkout')

    from chython import MoleculeContainer, QueryContainer, ReactionContainer

    owner = {'MoleculeContainer': MoleculeContainer, 'QueryContainer': QueryContainer,
             'ReactionContainer': ReactionContainer, 'mol': MoleculeContainer,
             'molecule': MoleculeContainer, 'rxn': ReactionContainer, 'reaction': ReactionContainer}

    bad = []
    for number, text in _lines(_source(root, source)):
        for klass, member, variable, attribute in _MEMBER.findall(text):
            klass, member = klass or variable, member or attribute
            if member.startswith('_') or member in _SUFFIXES:
                continue
            if not hasattr(owner[klass], member):
                bad.append(f'{source}:{number} {klass}.{member}')
    assert not bad, 'documented members that do not exist:\n' + '\n'.join(bad)
