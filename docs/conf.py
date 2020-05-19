# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# http://www.sphinx-doc.org/en/master/config

import os

from recommonmark.transform import AutoStructify

# check if we are running on readthedocs.org
on_readthedocs = os.environ.get('READTHEDOCS', None) == 'True'

# always build the API docs w/ doxygen on RTD
if on_readthedocs:
    tags.add('use_doxygen')

# -- Project information ------------------------------------------------------

project = 'Acts'
author = 'The Acts authors'
copyright = '2014–2020 CERN for the benefit of the Acts project'
# version = '@PROJECT_VERSION@'
# release = '@PROJECT_VERSION@'

# -- General ------------------------------------------------------------------

extensions = [
    'sphinx.ext.mathjax',
    'recommonmark',
    'sphinx_markdown_tables',
]
source_suffix = {
    '.rst': 'restructuredtext',
    '.md': 'markdown',
}
master_doc = 'index'
# ensure the in-source build directory is ignored
exclude_patterns = [
    '_build',
]
# cpp as default language
primary_domain = 'cpp'
highlight_language = 'cpp'
smartquotes = True

# -- Options for HTML output --------------------------------------------------

# ensure we use the RTD them when building locally
if not on_readthedocs:
    import sphinx_rtd_theme
    html_theme = 'sphinx_rtd_theme'
    html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]

html_theme_options = {
    'collapse_navigation': False,
    'navigation_depth': 2,
    'prev_next_buttons_location': None, # no prev/next links
    'style_external_links': True,
}
html_logo = 'figures/ActsLogo.gif'
html_static_path = []
html_copy_source = False
html_show_sourcelink = False
html_show_sphinx = False

# -- Doxygen integration with Breathe+Exhale ----------------------------------

if tags.has('use_doxygen'):
    extensions += [
        'breathe',
        'exhale',
    ]
    breathe_projects = {
        'Acts': '_build/doxygen-xml',
    }
    breathe_default_project = 'Acts'
    breathe_domain_by_extension = {
        'cpp': 'cpp',
        'hpp': 'cpp',
        'ipp': 'cpp',
    }
    breathe_default_members = ('members', 'undoc-members')
    exhale_args = {
        'containmentFolder': 'api',
        'rootFileName': 'api.rst',
        'rootFileTitle': 'API',
        'createTreeView': True,
        # let exhale handle the doxygen execution
        'exhaleExecutesDoxygen': True,
        'exhaleUseDoxyfile': True,
        # OUTPUT_DIRECTORY in Doxyfile must match breathe default project path
        # this must match STRIP_FROM_PATH in the Doxyfile
        'doxygenStripFromPath': '..',
    }

# -- Markdown bridge setup hook (must come last, not sure why) ----------------

def setup(app):
    app.add_config_value(
        'recommonmark_config',
        {
            'enable_math': True,
            'enable_inline_math': True,
        },
        True)
    app.add_transform(AutoStructify)
