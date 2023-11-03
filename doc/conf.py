# -*- coding: utf-8 -*-
from os.path import abspath
from sys import path
parent = abspath('..')
if parent not in path:
    path.insert(0, parent)
from chython.periodictable import C, QueryC, ListElement, AnyElement, AnyMetal

author = 'Dr. Ramil Nugmanov'
copyright = '2014-2023, Dr. Ramil Nugmanov <nougmanoff@protonmail.com>'
version = '1.x'
project = 'chython'

needs_sphinx = '1.8'
extensions = ['sphinx.ext.autodoc', 'sphinx.ext.autosummary', 'nbsphinx']

exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', '**.ipynb_checkpoints']
templates_path = ['_templates']
source_suffix = '.rst'
master_doc = 'index'

language = 'en'
pygments_style = 'sphinx'
todo_include_todos = False
autoclass_content = 'both'

html_logo = 'logo256.png'
html_favicon = 'logo256.png'
html_theme_options = {'github_user': 'chython', 'github_repo': 'chython', 'show_related': True}
html_show_copyright = True
html_show_sourcelink = False
html_sidebars = {
    '**': [
        'about.html',
        'navigation.html',
        'relations.html',  # needs 'show_related': True theme option to display
        'searchbox.html',
    ]
}

nbsphinx_execute = 'never'
