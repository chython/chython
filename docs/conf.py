# -*- coding: utf-8 -*-
#
#  Copyright 2021-2026 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from os.path import abspath
from sys import path

parent = abspath('..')
if parent not in path:
    path.insert(0, parent)

project = 'chython'
author = 'Dr. Ramil Nugmanov'
copyright = '2014-2026, Dr. Ramil Nugmanov'
version = '3.x'

needs_sphinx = '7.0'
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.doctest',
    'sphinx.ext.viewcode',
]

# NO `autodoc_mock_imports = ['chython']`.  There is not one autodoc directive in the whole of `docs/`,
# and the first `automodule:: chython` anybody adds must render against the real package: mocked, it
# would render an empty page and still pass.
#
# `sphinx.ext.doctest` is the comparator: every Python sample is a `testcode::` block, so `make doctest`
# executes them.  `chython/test/test_doc_samples.py` runs the same blocks under pytest -- because the
# gate has to fire for a developer who never builds the docs, which is every developer -- and
# additionally forbids `code-block:: python`, the spelling that renders identically and runs never.

# The IO samples write the files they then read back, so run every group in a scratch directory -- the
# same isolation `test_doc_samples.py` gives its side.  Without it `make doctest` drops a dozen
# `output.sdf`-shaped files in whatever directory it was invoked from, which is the repository root.
doctest_global_setup = '''
from os import chdir as _chdir
from tempfile import mkdtemp as _mkdtemp
_chdir(_mkdtemp(prefix='chython-doctest-'))
'''

exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', 'tutorial']
source_suffix = '.rst'
master_doc = 'index'
language = 'en'
pygments_style = 'default'
pygments_dark_style = 'monokai'

# -- 3D scenes ---------------------------------------------------------------
# `depict3d()` returns `<x3d>` markup and the X3DOM runtime is NOT in it: in a notebook
# `JupyterWidget._repr_html_` carries the two tags itself, and here the page carries them, once -- so a
# page with two scenes cannot load x3dom twice.  From x3dom.org rather than vendored into `_static/`: a
# reader offline, or a host with a strict CSP, sees an empty box, which is the failure the notebook
# widget already has.  `docs/figures.py` writes the scenes; `.. raw:: html` with `:file:` inlines one.
html_js_files = [('https://www.x3dom.org/download/x3dom.js', {'defer': 'defer'})]
html_css_files = ['https://www.x3dom.org/download/x3dom.css']

# -- Theme -------------------------------------------------------------------
html_theme = 'furo'
html_title = 'chython'
html_logo = 'logo256.png'
html_favicon = 'logo256.png'
html_show_sourcelink = False
html_show_copyright = True

html_theme_options = {
    'sidebar_hide_name': True,
    'navigation_with_keys': True,
    'source_repository': 'https://github.com/chython/chython',
    'source_branch': 'master',
    'source_directory': 'docs/',
    'light_css_variables': {
        'color-brand-primary': '#2962ff',
        'color-brand-content': '#2962ff',
    },
    'dark_css_variables': {
        'color-brand-primary': '#82b1ff',
        'color-brand-content': '#82b1ff',
    },
}
