# -*- coding: utf-8 -*-
from os.path import abspath
from sys import path

parent = abspath('..')
if parent not in path:
    path.insert(0, parent)

project = 'chython'
author = 'Dr. Ramil Nugmanov'
copyright = '2014-2026, Dr. Ramil Nugmanov'
version = '1.x'

needs_sphinx = '7.0'
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.viewcode',
]

autodoc_mock_imports = ['chython']

exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', 'tutorial']
source_suffix = '.rst'
master_doc = 'index'
language = 'en'
pygments_style = 'default'
pygments_dark_style = 'monokai'

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
    'source_directory': 'doc/',
    'light_css_variables': {
        'color-brand-primary': '#2962ff',
        'color-brand-content': '#2962ff',
    },
    'dark_css_variables': {
        'color-brand-primary': '#82b1ff',
        'color-brand-content': '#82b1ff',
    },
}
