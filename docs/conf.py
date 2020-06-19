# Configuration file for the Sphinx documentation builder.

import os

from recommonmark.transform import AutoStructify

# check if we are running on readthedocs.org
on_readthedocs = os.environ.get('READTHEDOCS', None) == 'True'

# -- Project information ------------------------------------------------------

project = 'Acts'
author = 'The Acts authors'
copyright = '2014–2020 CERN for the benefit of the Acts project'
# version = '@PROJECT_VERSION@'
# release = '@PROJECT_VERSION@'

# -- General ------------------------------------------------------------------

extensions = [
    'breathe',
    'recommonmark',
    'sphinx.ext.mathjax',
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
    html_theme_path = [
        sphinx_rtd_theme.get_html_theme_path(),
    ]

html_theme_options = {
    'collapse_navigation': False,
    'navigation_depth': 2,
    'prev_next_buttons_location': None, # no prev/next links
    'style_external_links': True,
}
html_logo = 'figures/acts_logo_white.svg'
html_static_path = [
    '_static',
]
html_css_files = [
    'custom.css',
]
html_copy_source = False
html_show_sourcelink = False
html_show_sphinx = False

# -- Doxygen integration with Breathe -----------------------------------------

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

# -- Automatic API documentation generation with Exhale -----------------------

exhale_args = {
    'containmentFolder': 'api',
    'rootFileName': 'api.rst',
    'rootFileTitle': 'API',
    'createTreeView': True,
    'exhaleUseDoxyfile': True,
    # note: OUTPUT_DIRECTORY in Doxyfile must match breathe default project path
    # note: this must match STRIP_FROM_PATH in Doxyfile
    'doxygenStripFromPath': '..',
}

if on_readthedocs:
    # if we are running on RTD Doxygen must be run as part of the build
    # let exhale handle the doxygen execution
    extensions.append('exhale')
    exhale_args['exhaleExecutesDoxygen'] = True
elif tags.has('use_exhale'):
    # if exhale is requested manually, we expect the Doxygen has been run for us
    extensions.append('exhale')
    exhale_args['exhaleExecutesDoxygen'] = False

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


    from m2r import MdInclude
    # from m2r to make `mdinclude` work
    app.add_config_value('no_underscore_emphasis', False, 'env')
    app.add_config_value('m2r_parse_relative_links', False, 'env')
    app.add_config_value('m2r_anonymous_references', False, 'env')
    app.add_config_value('m2r_disable_inline_math', False, 'env')
    app.add_directive('mdinclude', MdInclude)
