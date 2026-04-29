#!/usr/bin/env python3

from pathlib import Path
import os
import sys
import subprocess

EXCLUDE_PATHS = (
    "Scripts",
    "thirdparty",
    "CI",
    "Python",
    "git",
    "cmake",
    ".git",
    ".github",
    ".idea",
    ".devcontainer",
    "white_papers/figures",
    # cmake for "DD4hep-tests" looks a bit different
    "Tests/UnitTests/Plugins/DD4hep",
)
EXCLUDE_FILES = (
    "acts_logo_colored.svg",
    ".gitignore",
    "README.md",
    "CMakeLists.txt",
    "pytest.ini",
    ".pre-commit-config.yaml",
    "CITATION.cff",
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
    "odd-digi-smearing-config-notime.json",
    # TODO Mention these files somewhere?
    "generate_particle_data_table.py",
    "lazy_autodoc.py",
    "codegen/src/codegen/sympy_common.py",
    "CompressedIO.h",
    "GeometryModule.h",
    "runtime_geometry_modules.md",
    # Files for python binding generation
    "tgeo_aux.py.in",
    "serve.py",
    "SNIPPETS.md",
    "todo.md",
    "bugs.md",
    "deprecated.md",
    "acts-version-manager.js",
    "tex-mml-chtml.js",
    "Python/conftest.py",
    # Temporarily excluded files. TODO remove in next major release.
    "Core/include/Acts/EventData/detail/ParameterTraits.hpp",
    "Core/include/Acts/Seeding/PathSeeder.hpp",
    "Tests/CommonHelpers/include/ActsTests/CommonHelpers/TestSpacePoint.hpp",
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
    unused_files = ()

    # walk over all files
    for root, dirs, files in os.walk("."):
        dirs[:] = filter_paths(dirs, root, EXCLUDE_PATHS)
        files = filter_paths(files, root, EXCLUDE_PATHS, EXCLUDE_FILES)

        # Skip base-directory
        if str(Path(root)) == ".":
            continue

        # Print progress
        if root[2:] in dirs_base:
            processed_files = 0
            current_base_dir = root
            number_files = count_files(root)
            # print empty to start a new line
            print("")

        root = Path(root)
        for filename in files:
            processed_files += 1
            # get the full path of the file
            filepath = root / filename

            # Check header files and remove
            if filepath.suffix in SUFFIX_CPP:
                if file_can_be_removed(filename, dirs_base_code):
                    unused_files += (str(filepath),)
                    remove_cmd = "rm " + str(filepath)
                    os.system(remove_cmd)

            elif filepath.suffix in SUFFIX_PYTHON:
                if not file_can_be_removed("import .*" + filepath.stem, dirs_base):
                    continue

                if not file_can_be_removed(
                    "from " + filepath.stem + " import", dirs_base
                ):
                    continue

                if file_can_be_removed(filename, dirs_base):
                    unused_files += (str(filepath),)
                    remove_cmd = "rm " + str(filepath)
                    os.system(remove_cmd)

            # Check documentation files (weak tests)
            # TODO find more reliable test for this
            elif filepath.suffix in SUFFIX_DOC:
                if file_can_be_removed(filepath.stem, dirs_base_docs):
                    unused_files += (str(filepath),)
                    remove_cmd = "rm " + str(filepath)
                    os.system(remove_cmd)

            # Check and print other files
            elif filepath.suffix in SUFFIX_IMAGE + SUFFIX_OTHER:
                if file_can_be_removed(filename, dirs_base):
                    unused_files += (str(filepath),)
                    remove_cmd = "rm " + str(filepath)
                    os.system(remove_cmd)

        # Print the progress in place
        progress = int(20 * processed_files / number_files)
        sys.stdout.write("\r")
        sys.stdout.write(
            "Checked [%-20s] %d/%d files in %s"
            % ("=" * progress, processed_files, number_files, current_base_dir)
        )
        sys.stdout.flush()

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
