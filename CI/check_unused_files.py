#!/usr/bin/env python3

from pathlib import Path
import os
import sys
import subprocess


def file_can_be_removed(searchstring, scope):
    cmd = "grep -IR '" + searchstring + "' " + " ".join(scope)

    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    if p.returncode == 1:
        return True

    return False


def main():
    exclude_dirs = (
        "Scripts",
        "thirdparty",
        "CI",
        "git",
        "cmake",
        ".git",
        ".github",
        ".",
        ".idea",
    )
    exclude_files = (
        "acts_logo_colored.svg",
        ".gitignore",
        "README.md",
        "CMakeLists.txt",
        # Filename not completed in source
        "vertexing_event_mu20_beamspot.csv",
        "vertexing_event_mu20_tracks.csv",
        "vertexing_event_mu20_vertices_AMVF.csv",
        # TODO Move the following files to a better place?
        "Magfield.ipynb",
        "SolenoidField.ipynb",
        # TODO Add README next to the following files?
        "default-input-config-generic.json",
        "geoSelection-openDataDetector.json",
        "alignment-geo-contextualDetector.json",
    )

    suffix_header = (
        ".hpp",
        ".cuh",
        ".h",  # TODO discontinue .h over time
    )
    suffix_source = (
        ".ipp",
        ".cpp",
        ".cu",
    )
    suffix_image = (
        ".png",
        ".svg",
        ".jpg",
        ".gif",
    )
    suffix_python = (".py",)
    suffix_doc = (
        ".md",
        ".rst",
    )
    suffix_other = (
        "",
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
    )
    suffix_allowed = (
        suffix_header
        + suffix_source
        + suffix_image
        + suffix_python
        + suffix_doc
        + suffix_other
    )

    exit = 0

    dirs_base = next(os.walk("."))[1]
    dirs_base[:] = [d for d in dirs_base if d not in exclude_dirs]
    dirs_base_docs = ("docs",)
    dirs_base_code = [d for d in dirs_base if d not in dirs_base_docs]
    print(dirs_base)

    # Collectors
    wrong_extension = ()
    unused_files = ()

    # walk over all files
    for root, dirs, files in os.walk("."):  # topdown=True
        dirs[:] = [d for d in dirs if d not in exclude_dirs]
        files[:] = [f for f in files if f not in exclude_files]

        # Skip "white-paper-figures"
        # TODO Find a more elegant way
        if str(root).find("white_papers/figures") != -1:
            continue

        # Skip "DD4hep-tests" since their cmake looks a bit different
        # TODO Find a more elegant way
        if str(root).find("Tests/UnitTests/Plugins/DD4hep") != -1:
            continue

        root = Path(root)

        # Skip base-directory
        if str(root) == ".":
            continue

        for filename in files:
            # get the full path of the file
            filepath = root / filename

            # Check for wrong extensions
            if filepath.suffix not in suffix_allowed:
                wrong_extension += (str(filepath),)

            # Check header files and remove
            if filepath.suffix in suffix_header + suffix_source:
                if file_can_be_removed(filepath.stem, dirs_base_code):
                    unused_files += (str(filepath),)
                    remove_cmd = "rm " + str(filepath)
                    os.system(remove_cmd)

            # TODO Find test to check python files
            if filepath.suffix in suffix_python:
                continue

            # Check documentation files (weak tests)
            # TODO find more reliable test for this
            if filepath.suffix in suffix_doc:
                if file_can_be_removed(filepath.stem, dirs_base_docs):
                    unused_files += (str(filepath),)
                    remove_cmd = "rm " + str(filepath)
                    os.system(remove_cmd)

            # Check and print other files
            if filepath.suffix in suffix_image + suffix_other:
                if file_can_be_removed(filename, dirs_base):
                    unused_files += (str(filepath),)
                    remove_cmd = "rm " + str(filepath)
                    os.system(remove_cmd)

    if len(wrong_extension) != 0:
        print(
            f"ERROR: The following {len(wrong_extension)} files have an unsupported extension:\n\n"
            + "\n".join(wrong_extension)
            + "\nCheck if you can change the format to one of the following:\n"
            + "\n".join(suffix_allowed)
            + "\nIf you really need that specific extension, add it to the list above.\n"
        )

        exit += 1

    if len(unused_files) != 0:
        print(
            f"ERROR: The following {len(unused_files)} files seem to be unused:\n"
            + "\n".join(unused_files)
            + "\nYou have 3 options:"
            + "\n\t- Remove them"
            + "\n\t- Use them (check proper include)"
            + "\n\t- Modify the ignore list of this check\n"
        )

        exit += 1

    if exit == 0:
        print("Finished check_unused_files.py without any errors.")

    return exit


if "__main__" == __name__:
    sys.exit(main())
