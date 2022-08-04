# Configuration file for the Sphinx documentation builder.

import os
import sys
import subprocess

# check if we are running on readthedocs.org
on_readthedocs = os.environ.get("READTHEDOCS", None) == "True"

# -- Project information ------------------------------------------------------

project = "Acts"
author = "The Acts authors"
copyright = "2014–2022 CERN for the benefit of the Acts project"
# version = '@PROJECT_VERSION@'
# release = '@PROJECT_VERSION@'

# -- General ------------------------------------------------------------------

extensions = [
    "breathe",
    "myst_parser",
    "sphinx.ext.mathjax",
]
source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}
master_doc = "index"
# ensure the in-source build directory is ignored
exclude_patterns = [
    "_build",
]
# cpp as default language
primary_domain = "cpp"
highlight_language = "cpp"
smartquotes = True
numfig = True

myst_enable_extensions = ["dollarmath", "colon_fence"]
myst_heading_anchors = 3

# -- Options for HTML output --------------------------------------------------

# ensure we use the RTD them when building locally
if not on_readthedocs:
    html_theme = "sphinx_rtd_theme"
    extensions.append("sphinx_rtd_theme")

html_theme_options = {
    "collapse_navigation": False,
    "navigation_depth": 3,
    "prev_next_buttons_location": None,  # no prev/next links
    "style_external_links": True,
}
html_logo = "figures/acts_logo_white.svg"
html_static_path = [
    "_static",
]
html_css_files = [
    "custom.css",
]
html_copy_source = False
html_show_sourcelink = False
html_show_sphinx = False

# -- Doxygen integration with Breathe -----------------------------------------

breathe_projects = {
    "Acts": "_build/doxygen-xml",
}
breathe_default_project = "Acts"
breathe_domain_by_extension = {
    "cpp": "cpp",
    "hpp": "cpp",
    "ipp": "cpp",
}
breathe_default_members = (
    "members",
    "undoc-members",
)

# -- Automatic API documentation ---------------------------------------------

env = os.environ.copy()
env["DOXYGEN_WARN_AS_ERROR"] = "NO"
cwd = os.path.dirname(__file__)

if on_readthedocs or tags.has("run_doxygen"):
    # if we are running on RTD Doxygen must be run as part of the build
    print("Executing doxygen in", cwd)
    print(
        "Doxygen version:",
        subprocess.check_output(["doxygen", "--version"], encoding="utf-8"),
    )
    sys.stdout.flush()
    subprocess.check_call(
        ["doxygen", "Doxyfile"], stdout=subprocess.PIPE, cwd=cwd, env=env
    )

if on_readthedocs or tags.has("run_apidoc"):
    print("Executing breathe apidoc in", cwd)
    subprocess.check_call(
        [sys.executable, "-m", "breathe.apidoc", "_build/doxygen-xml", "-o", "api"],
        stdout=subprocess.PIPE,
        cwd=cwd,
        env=env,
    )

# -- Markdown bridge setup hook (must come last, not sure why) ----------------


def setup(app):
    pass
