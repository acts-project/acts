# Configuration file for the Sphinx documentation builder.

import os
import sys
import re
import subprocess
from pathlib import Path
import shutil
import datetime
from typing import List, Tuple

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

sys.path.insert(0, str(doc_dir / "_extensions"))

extensions = [
    "breathe",
    "myst_parser",
    "sphinx.ext.mathjax",
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

myst_enable_extensions = ["dollarmath", "colon_fence", "amsmath"]
myst_heading_anchors = 3

linkcheck_retries = 5
linkcheck_ignore = [
    r"https://doi.org/.*",
    r"https://cernvm.cern.ch/.*",
    r"http://eigen.tuxfamily.org.*",
    r"https://pythia.org.*",
    r"https://lcginfo.cern.ch/.*",
    r"https://.*\.?intel.com/.*",
]

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
    ("cpp:identifier", "SIZE_MAX"),
    ("cpp:identifier", "M_PI"),
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
    ("cpp:identifier", r".*"),
]

# -- Automatic API documentation ---------------------------------------------

env = os.environ.copy()
env["DOXYGEN_WARN_AS_ERROR"] = "NO"

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
else:
    #  if not api_index_target.exists():
    #  shutil.copyfile(doc_dir / "api/api_stub.rst", api_index_target)
    pass

# Find all roles pointing at code
# @TODO: Member is difficult, unscoped enums don't work
roles = [
    "class",
    "struct",
    "type",
]

role_names = {"class": "Classes", "struct": "Structs", "type": "Types", "enum": "Enums"}

directives = {
    "class": "doxygenclass",
    "struct": "doxygenstruct",
    "type": "doxygentypedef",
    "func": "doxygenfunction",
    "enum": "doxygenenum",
}

role_instances = {k: set() for k in roles}

role_instances["type"] |= {
    "Acts::ActsScalar",
    "Acts::ActsVector",
    "Acts::ActsMatrix",
    "Acts::ActsSquareMatrix",
    "Acts::SquareMatrix2",
    "Acts::SquareMatrix3",
    "Acts::SquareMatrix4",
    "Acts::BoundMatrix",
    "Acts::BoundSquareMatrix",
    "Acts::Vector2",
    "Acts::Vector3",
    "Acts::Vector4",
    "Acts::BoundVector",
    "Acts::BoundTrackParameters",
    "Acts::Transform2",
    "Acts::Transform3",
    "Acts::AngleAxis3",
    "Acts::RotationMatrix2",
    "Acts::RotationMatrix3",
    "Acts::Translation2",
    "Acts::Translation3",
    "Acts::GeometryContext",
    "Acts::FreeVector",
    "Acts::FreeMatrix",
    "Acts::SurfaceVector",
    "Acts::Intersection3D",
    "Acts::OrientedSurface",
    "Acts::OrientedSurfaces",
    "Acts::BoundToFreeMatrix",
    "Acts::FreeToBoundMatrix",
    "Acts::FreeSquareMatrix",
    "Acts::FreeToPathMatrix",
}

role_instances["struct"] |= {
    "Acts::DenseStepperPropagatorOptions",
    "Acts::Experimental::DetectorNavigator::State",
    "Acts::Geant4PhysicalVolumeSelectors::AllSelector",
    "Acts::Geant4PhysicalVolumeSelectors::NameSelector",
}

role_instances["class"] |= {
    "Acts::BinningData",
    "Acts::Direction",
    "Acts::ConstrainedStep",
    "Acts::IAxis",
    "Acts::SeedFilter",
    "Acts::BoundaryCheck",
    "Acts::ConeVolumeBounds",
    "Acts::CuboidVolumeBounds",
    "Acts::CylinderVolumeBounds",
    "Acts::CutoutCylinderVolumeBounds",
    "Acts::GenericCuboidVolumeBounds",
    "Acts::TrapezoidVolumeBounds",
    "Acts::GeometryObject",
    "Acts::TrackContainer",
    "Acts::ConeLayer",
    "Acts::CylinderLayer",
    "Acts::IdentifiedDetectorElement",
    "Acts::DiscLayer",
    "Acts::PlaneLayer",
    "Acts::NullBField",
    "Acts::DiscBounds",
    "Acts::PlanarBounds",
    "Acts::AbstractVolume",
    "Acts::AnnulusBounds",
    "Acts::DiamondBounds",
    "Acts::ConvexPolygonBounds",
    "Acts::ConvexPolygonBoundsBase",
    "Acts::Logging::LevelOutputDecorator",
    "Acts::Logging::NamedOutputDecorator",
    "Acts::Logging::ThreadOutputDecorator",
    "Acts::Logging::TimedOutputDecorator",
    "Acts::Logging::DefaultFilterPolicy",
    "Acts::Logging::DefaultPrintPolicy",
}

role_instances["enum"] = {
    "Acts::BinningValue",
    "Acts::BinningType",
    "Acts::BinningValue",
    "Acts::BoundIndices",
    "Acts::FreeIndices",
}

role_ex = re.compile(r"[{:](" + "|".join(roles) + r")[}:]`(.+?)`")


def process_roles(file: Path) -> List[Tuple[str, str]]:
    text = file.read_text()
    return [m.groups() for m in role_ex.finditer(text)]


for dirpath, _, filenames in os.walk(doc_dir):
    dirpath = Path(dirpath)
    for file in filenames:
        file = dirpath / file
        if file.suffix not in (".rst", ".md"):
            continue
        for role, arg in process_roles(file):
            role_instances[role].add(arg)

# add members to their parents

api_preamble = """

"""

with api_index_target.open("w") as fh:
    fh.write("# API Reference\n\n")
    fh.write(api_preamble)
    for role, instances in sorted(role_instances.items(), key=lambda x: x[0]):
        fh.write(f"## {role_names[role]}\n")
        for instance in sorted(instances):
            fh.write(
                f"""
:::{{{directives[role]}}} {instance}
:::
"""
            )
        fh.write("\n")


# -- Markdown bridge setup hook (must come last, not sure why) ----------------


def setup(app):
    pass
