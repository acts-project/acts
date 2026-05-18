# Configuration file for the Sphinx documentation builder.

import os
import sys
import subprocess
from pathlib import Path
import shutil
import datetime
import urllib.request
import urllib.error
import json
import html

# check if we are running on readthedocs.org
on_readthedocs = os.environ.get("READTHEDOCS", None) == "True"

# -- Project information ------------------------------------------------------

project = "Acts"
author = "The Acts authors"
copyright = (
    f"2014–{datetime.date.today().year} CERN for the benefit of the Acts project"
)
# version = '@PROJECT_VERSION@'
# release = '@PROJECT_VERSION@'

# -- General ------------------------------------------------------------------

doc_dir = Path(__file__).parent

sys.path.insert(0, str(doc_dir))
sys.path.insert(0, str(doc_dir / "_extensions"))

extensions = [
    "breathe",
    "myst_parser",
    "sphinx.ext.mathjax",
    "sphinx.ext.graphviz",
    "sphinx.ext.todo",
    "warnings_filter",
]

todo_include_todos = True

warnings_filter_config = str(doc_dir / "known-warnings.txt")
warnings_filter_silent = True

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}
master_doc = "index"
# ensure the in-source build directory is ignored
exclude_patterns = ["_build", "api/api_stub.rst", "api/api_index.rst"]
# cpp as default language
primary_domain = "cpp"
highlight_language = "cpp"
smartquotes = True
numfig = True

myst_enable_extensions = ["dollarmath", "colon_fence", "amsmath", "html_image"]
myst_heading_anchors = 3
myst_dmath_allow_labels = True

linkcheck_retries = 5
linkcheck_ignore = []

# Linkcheck ignore patterns are loaded from this URL, so we can
# update without adding pull requests.
linkcheck_ignore_url = (
    "https://raw.githubusercontent.com/acts-project/linkcheck-ignore/main/data.json"
)
try:
    response = urllib.request.urlopen(linkcheck_ignore_url)
    linkcheck_ignore = json.loads(response.read().decode("utf-8"))
except urllib.error.HTTPError:
    print("Error getting linkcheck ignore data, using default")

print("Link check ignore patterns")
print(linkcheck_ignore)


# -- Options for HTML output --------------------------------------------------

html_theme = "sphinx_rtd_theme"
extensions.append("sphinx_rtd_theme")

html_theme_options = {
    "collapse_navigation": False,
    "navigation_depth": 4,
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

nitpicky = True
nitpick_ignore = [
    ("cpp:identifier", "Acts"),
    ("cpp:identifier", "detail"),
    ("cpp:identifier", "eSize"),
    ("cpp:identifier", "eBoundSize"),
    ("cpp:identifier", "eFreeSize"),
    ("cpp:identifier", "open"),
    ("cpp:identifier", "FreeToBoundCorrection"),
]

nitpick_ignore_regex = [
    ("cpp:identifier", r"Eigen.*"),
    ("cpp:identifier", r"boost.*"),
    ("cpp:identifier", r"s_.*"),
    ("cpp:identifier", r"detail::.*"),
    ("cpp:identifier", ".*::Identity"),
    ("cpp:identifier", ".*::Zero"),
    # This blanket ignore only targets the doxygen/breathe auto-generated
    # references. Explicit references should have specific types.
    ("cpp:identifier", r".*"),
]

# -- Automatic API documentation ---------------------------------------------

env = os.environ.copy()

if on_readthedocs or tags.has("run_doxygen"):
    # if we are running on RTD Doxygen must be run as part of the build
    print("Executing doxygen in", doc_dir)
    print(
        "Doxygen version:",
        subprocess.check_output(["doxygen", "--version"], encoding="utf-8"),
    )
    sys.stdout.flush()
    subprocess.check_call(
        ["doxygen", "Doxyfile"], stdout=subprocess.PIPE, cwd=doc_dir, env=env
    )

api_index_target = doc_dir / "api/api.md"

if tags.has("run_apidoc"):
    print("Executing breathe apidoc in", doc_dir)
    subprocess.check_call(
        [sys.executable, "-m", "breathe.apidoc", "_build/doxygen-xml", "-o", "api"],
        stdout=subprocess.DEVNULL,
        cwd=doc_dir,
        env=env,
    )
    if not api_index_target.exists():
        shutil.copyfile(doc_dir / "api/api_index.rst", api_index_target)
    print("breathe apidoc completed")

if tags.has("lazy_autodoc") or on_readthedocs:
    extensions += ["lazy_autodoc"]


if on_readthedocs or tags.has("white_papers"):
    import white_papers

    white_papers.render()

# -- Markdown bridge setup hook (must come last, not sure why) ----------------


def setup(app):
    legacy_redirects = {
        "index": "https://acts-project.github.io/",
        "getting_started": "https://acts-project.github.io/building_acts.html",
        "tracking": "https://acts-project.github.io/tracking.html",
        "versioning": "https://acts-project.github.io/versioning.html",
        "examples/python_bindings": "https://acts-project.github.io/python_bindings.html",
        "contribution/clang_tidy": "https://acts-project.github.io/contribution_clang_tidy.html",
        "contribution/profiling": "https://acts-project.github.io/howto_profiling.html",
        "contribution/release": "https://acts-project.github.io/howto_release.html",
        "contribution/run_formatting": "https://acts-project.github.io/formatting.html",
        "contribution/physmon": "https://acts-project.github.io/physmon.html",
        "contribution/root_hash_checks": "https://acts-project.github.io/python_bindings.html#root_file_hashes",
        "contribution/documentation_build": "https://acts-project.github.io/building_acts.html",
        "misc/spack": "https://acts-project.github.io/howto_spack.html",
    }

    def add_legacy_redirect(app, pagename, templatename, context, doctree):
        target = legacy_redirects.get(pagename)
        if target is None:
            return

        escaped_target = html.escape(target, quote=True)
        redirect_snippet = f"""
<meta http-equiv="refresh" content="0; url={escaped_target}" />
<link rel="canonical" href="{escaped_target}" />
<script>
(function() {{
  var target = "{escaped_target}";
  var destination = target.indexOf("#") === -1 ? target + window.location.hash : target;
  window.location.replace(destination);
}})();
</script>
"""
        context["metatags"] = context.get("metatags", "") + redirect_snippet

    app.connect("html-page-context", add_legacy_redirect)
