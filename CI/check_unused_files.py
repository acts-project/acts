#!/usr/bin/env python3

from pathlib import Path
import os
import sys
import subprocess

EXCLUDE_PATHS = (
    ".devcontainer",
    ".git",
    ".github",
    ".idea",
    "CI",
    "cmake",
    # Used by traccc
    "Detray/detectors",
    # CLI tools
    "Detray/tests/tools",
    # TODO: Remove the traccc part.
    "Traccc",
    "git",
    "Python",
    "Scripts",
    # cmake for "DD4hep-tests" looks a bit different
    "Tests/UnitTests/Plugins/DD4hep",
    "thirdparty",
    "white_papers/figures",
)
EXCLUDE_FILES = (
    ".gersemirc",
    ".gitignore",
    ".kodiak.toml",
    ".merge-sentinel.yml",
    ".policy.yml",
    ".pre-commit-config.yaml",
    "acts_logo_colored.svg",
    "CITATION.cff",
    "CMakeLists.txt",
    "CMakePresets.json",
    "CODE_OF_CONDUCT.md",
    "CODEOWNERS",
    "codecov.yml",
    "pytest.ini",
    "README.md",
    "readthedocs.yml",
    "sonar-project.properties",
    # Filename not completed in source
    "vertexing_event_mu20_beamspot.csv",
    "vertexing_event_mu20_tracks.csv",
    "vertexing_event_mu20_vertices_AMVF.csv",
    "event000000001-MuonDriftCircle.csv",
    "event000000001-MuonSimHit.csv",
    # TODO Move the following files to a better place?
    "Magfield.ipynb",
    "SolenoidField.ipynb",
    # TODO Add README next to the following files?
    "generic-input-config.json",
    "generic-alignment-geo.json",
    # TODO Mention these files somewhere?
    "codegen/src/codegen/sympy_common.py",
    "codegen/src/codegen/detray_backend.py",
    "CompressedIO.h",
    "generate_particle_data_table.py",
    "GeometryModule.h",
    "lazy_autodoc.py",
    "runtime_geometry_modules.md",
    # Files for python binding generation
    "acts-version-manager.js",
    "bugs.md",
    "deprecated.md",
    "Python/conftest.py",
    "serve.py",
    "SNIPPETS.md",
    "tex-mml-chtml.js",
    "tgeo_aux.py.in",
    "todo.md",
    # Detray python tests for auto-generated code
    "Detray/codegen/detray-sympy/tests/test_assumptions_D.py",
    "Detray/codegen/detray-sympy/tests/test_matrices.py",
    # Used in traccc
    "Detray/tests/include/detray/test/utils/perigee_stopper.hpp",
    "Detray/tests/include/detray/test/validation/propagation_validation.hpp",
    # Build-time metadata generation
    "Detray/python/detray/detectors/impl/definitions.py",
    "Detray/python/detray/detectors/impl/type_helpers.py",
    # Python uv files
    "Detray/codegen/detray-sympy/uv.lock",
    "Detray/python/detray/uv.lock",
    # TODO: remove after file is gone
    "Core/include/Acts/Utilities/ProtoAxisHelpers.hpp",
)
SUFFIX_CPP = (
    ".hpp",
    ".cuh",
    ".sycl",
    ".hip",
    ".ipp",
    ".cpp",
    ".cu",
)
SUFFIX_IMAGE = (
    ".png",
    ".svg",
    ".jpg",
    ".gif",
)
SUFFIX_PYTHON = (".py",)
SUFFIX_DOC = (
    ".md",
    ".rst",
    ".dox",
    ".html",
    ".bib",
)
SUFFIX_OTHER = (
    "",
    ".C",
    ".csv",
    ".css",
    ".gdml",
    ".hepmc3",
    ".lock",
    ".ico",
    ".in",
    ".ipynb",
    ".json",
    ".j2",
    ".onnx",
    ".root",
    ".toml",
    ".txt",
    ".yml",
    ".xml",
    ".sh",
)


def filter_paths(names, root, exclude_paths=(), exclude_files=()):
    """
    Filter names from os.walk() based on path substrings and file rules.
    Excludes entries if their full path matches exclude_paths or exclude_files.
    """

    def keep(name):
        p = Path(root) / name
        p_str = p.as_posix()
        return not any(ep in p_str for ep in exclude_paths) and not any(
            ef in p_str if "/" in ef else p.name == ef for ef in exclude_files
        )

    return [name for name in names if keep(name)]


def file_can_be_removed(searchstring, scope):
    cmd = "grep -IR '" + searchstring + "' " + " ".join(scope)

    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    output, _ = p.communicate()
    return output == b""


def count_files(path="."):
    count = 0
    for root, dirs, files in os.walk(path):
        dirs[:] = filter_paths(dirs, root, EXCLUDE_PATHS)
        files = filter_paths(files, root, EXCLUDE_PATHS, EXCLUDE_FILES)

        count += len(files)

    return count


def check_wrong_extensions(walk_root):
    """
    Collect files with disallowed suffixes. Returns the number of problematic files.
    """

    suffix_allowed = (
        SUFFIX_CPP + SUFFIX_IMAGE + SUFFIX_PYTHON + SUFFIX_DOC + SUFFIX_OTHER
    )

    wrong_extension = []
    for root, dirs, files in os.walk(walk_root):
        dirs[:] = filter_paths(dirs, root, EXCLUDE_PATHS)
        files = filter_paths(files, root, EXCLUDE_PATHS, EXCLUDE_FILES)

        for f in files:
            p = Path(root) / f
            if p.suffix not in suffix_allowed:
                wrong_extension.append(str(p))

    if len(wrong_extension) != 0:
        print(
            "\n\n\033[31mERROR\033[0m "
            + f"The following {len(wrong_extension)} files have an unsupported extension:\n\n"
            + "\033[31m"
            + "\n".join(wrong_extension)
            + "\033[0m"
            + "\nCheck if you can change the format to one of the following:\n"
            + "\n".join(suffix_allowed)
            + "\nIf you really need that specific extension, add it to the list above.\n"
        )

    return len(wrong_extension)


def find_unused_by_suffix(walk_root, suffixes, search_key, search_scope):
    unused = []

    for root, dirs, files in os.walk(walk_root):
        dirs[:] = filter_paths(dirs, root, EXCLUDE_PATHS)
        files = filter_paths(files, root, EXCLUDE_PATHS, EXCLUDE_FILES)

        for f in files:
            p = Path(root) / f
            if p.suffix in suffixes and file_can_be_removed(
                search_key(p), search_scope
            ):
                unused.append(str(p))

    return unused


def find_unused_python_files(walk_root, dirs_base):
    unused = []

    for root, dirs, files in os.walk(walk_root):
        dirs[:] = filter_paths(dirs, root, EXCLUDE_PATHS)
        files = filter_paths(files, root, EXCLUDE_PATHS, EXCLUDE_FILES)

        for f in files:
            p = Path(root) / f
            if p.suffix not in SUFFIX_PYTHON:
                continue

            if not file_can_be_removed(r"import .*" + p.stem, dirs_base):
                continue

            if not file_can_be_removed(r"from " + p.stem + r" import", dirs_base):
                continue

            if file_can_be_removed(p.name, dirs_base):
                unused.append(str(p))

    return unused


def main():
    print("\033[32mINFO\033[0m Start check_unused_files.py ...")

    exit = 0

    dirs_base = next(os.walk("."))[1]
    dirs_base.append(".")
    dirs_base[:] = filter_paths(dirs_base, Path("."), EXCLUDE_PATHS)
    dirs_base_docs = ("docs",)
    dirs_base_code = filter_paths(dirs_base, Path("."), dirs_base_docs)

    exit += check_wrong_extensions(".")

    # Collector
    unused_files = []

    unused_files += find_unused_by_suffix(
        ".", SUFFIX_CPP, lambda p: p.name, dirs_base_code
    )

    unused_files += find_unused_python_files(".", dirs_base)

    # TODO find more reliable test for this
    unused_files += find_unused_by_suffix(
        ".", SUFFIX_DOC, lambda p: p.stem, dirs_base_docs
    )

    unused_files += find_unused_by_suffix(
        ".", SUFFIX_IMAGE + SUFFIX_OTHER, lambda p: p.name, dirs_base
    )

    if len(unused_files) != 0:
        print(
            "\n\n\033[31mERROR\033[0m "
            + f"The following {len(unused_files)} files seem to be unused:\n"
            + "\033[31m"
            + "\n".join(unused_files)
            + "\033[0m"
            + "\nYou have 3 options:"
            + "\n\t- Remove them"
            + "\n\t- Use them (check proper include)"
            + "\n\t- Modify the ignore list of this check\n"
        )

        exit += 1

    if exit == 0:
        print(
            "\n\n\033[32mINFO\033[0m Finished check_unused_files.py without any errors."
        )

    return exit


if "__main__" == __name__:
    sys.exit(main())
